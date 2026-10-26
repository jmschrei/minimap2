"""Locate or build minimap2, libminimap2.a, the `make extra` programs, the C
test drivers and the mappy extension module.

Nothing here writes to the source tree. Sources are copied to a cache
directory keyed by a hash of the sources and the compiler, and built there
with the Makefile's own flags (mappy with setup.py's). A file lock
serializes builds, so parallel pytest workers and concurrent test runs do
not race.

Environment variables:
	MINIMAP2_BIN              use this binary instead of building one
	MINIMAP2_TEST_SRC_DIR     source tree to build from (default: repository root)
	MINIMAP2_TEST_CACHE_DIR   build cache (default: $XDG_CACHE_HOME/minimap2-tests
	                          or ~/.cache/minimap2-tests)
	MINIMAP2_TEST_MAKE_JOBS   parallel jobs for make (default: 4)
	MINIMAP2_TEST_MAKE_ARGS   extra make arguments, e.g. "sse2only=1" or "asan=1"
	                          (default: none; part of the cache key)
	CC                        C compiler for make, the C drivers and mappy (default: gcc)
"""

import fcntl
import hashlib
import os
import shlex
import shutil
import subprocess
import sys
import sysconfig
import tempfile

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(TESTS_DIR)
C_DIR = os.path.join(TESTS_DIR, "c")

# Build products that `make`, `make extra` or setup.py would leave in a source
# tree; used to check that the harness never builds in place.
TREE_PRODUCT_PATTERNS = ("*.o", "*.a", "*.so", "*.mmi", "minimap2", "sdust", "minimap2-lite", "a.out",
	"mappy.c", "python/mappy.c", "build", "dist", "mappy.egg-info", "*.tmp")

# The files setup.py reads besides the C sources
MAPPY_FILES = ("setup.py", "python/mappy.pyx", "python/cmappy.pxd", "python/cmappy.h", "python/README.rst",
	"python/minimap2.py")


class BuildError(RuntimeError):
	pass


class MappyUnavailable(RuntimeError):
	"""Cython, setuptools or a C compiler is missing: mappy tests skip."""
	pass


def src_dir():
	return os.path.abspath(os.environ.get("MINIMAP2_TEST_SRC_DIR") or REPO_ROOT)


def compiler():
	return os.environ.get("CC") or "gcc"


def make_jobs():
	try:
		return max(1, int(os.environ.get("MINIMAP2_TEST_MAKE_JOBS", "4")))
	except ValueError:
		return 4


def make_args():
	"""Extra arguments for make from $MINIMAP2_TEST_MAKE_ARGS (shell words), e.g. ["sse2only=1"]."""
	return shlex.split(os.environ.get("MINIMAP2_TEST_MAKE_ARGS", ""))


def sanitizer_flags():
	"""Compiler flags the C drivers need to link against a sanitizer build of libminimap2.a: the
	Makefile's asan=1 and tsan=1 add -fsanitize=address or -fsanitize=thread to CFLAGS and LIBS."""
	out = []
	for a in make_args():
		k, _, v = a.partition("=")
		if k == "asan" and v:
			out.append("-fsanitize=address")
		elif k == "tsan" and v:
			out.append("-fsanitize=thread")
	return out


def cache_root():
	d = os.environ.get("MINIMAP2_TEST_CACHE_DIR")
	if not d:
		base = os.environ.get("XDG_CACHE_HOME") or os.path.join(os.path.expanduser("~"), ".cache")
		d = os.path.join(base, "minimap2-tests")
	os.makedirs(d, exist_ok=True)
	return d


def source_files(root):
	"""Relative paths of the files `make` reads: Makefile, top-level *.c/*.h and sse2neon/."""
	out = []
	for name in sorted(os.listdir(root)):
		p = os.path.join(root, name)
		if os.path.isfile(p) and (name == "Makefile" or name.endswith(".c") or name.endswith(".h")):
			out.append(name)
	for dirpath, dirnames, filenames in os.walk(os.path.join(root, "sse2neon")):
		dirnames.sort()
		for name in sorted(filenames):
			out.append(os.path.relpath(os.path.join(dirpath, name), root))
	if "Makefile" not in out or "minimap.h" not in out:
		raise BuildError("no Makefile or minimap.h in %s; set MINIMAP2_TEST_SRC_DIR to a minimap2 source tree" % root)
	return out


def mappy_source_files(root):
	"""The files setup.py reads: the C sources and headers plus setup.py and python/."""
	out = [f for f in source_files(root) if not f.startswith("sse2neon")] + list(MAPPY_FILES)
	for rel in MAPPY_FILES:
		if not os.path.exists(os.path.join(root, rel)):
			raise BuildError("%s has no %s" % (root, rel))
	return out


def _hash_files(h, root, rels):
	for rel in rels:
		h.update(rel.replace(os.sep, "/").encode() + b"\0")
		with open(os.path.join(root, rel), "rb") as f:
			h.update(hashlib.sha256(f.read()).digest())


def source_hash(root=None):
	root = root or src_dir()
	h = hashlib.sha256()
	_hash_files(h, root, source_files(root))
	return h.hexdigest()


def compiler_id(cc=None):
	"""A string that changes when the compiler does: its resolved path and version banner."""
	cc = cc or compiler()
	exe = shlex.split(cc)[0]
	path = shutil.which(exe) or exe
	try:
		ver = subprocess.run(shlex.split(cc) + ["--version"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
			timeout=60).stdout.decode("utf-8", "replace")
	except (OSError, subprocess.SubprocessError) as e:
		raise BuildError("cannot run compiler %r: %s" % (cc, e))
	return "%s\n%s\n%s" % (cc, os.path.realpath(path), ver)


def build_key(root=None, cc=None):
	h = hashlib.sha256()
	h.update(source_hash(root).encode())
	h.update(compiler_id(cc).encode())
	h.update(("%s-%s" % (sys.platform, os.uname().machine)).encode())
	if make_args():  # keys of builds without extra arguments are unchanged
		h.update(("make:" + "\0".join(make_args())).encode())
	return h.hexdigest()[:20]


class _Lock(object):
	def __init__(self, path):
		self.path = path
		self.fd = None

	def __enter__(self):
		self.fd = os.open(self.path, os.O_CREAT | os.O_RDWR, 0o644)
		fcntl.flock(self.fd, fcntl.LOCK_EX)
		return self

	def __exit__(self, *exc):
		fcntl.flock(self.fd, fcntl.LOCK_UN)
		os.close(self.fd)
		self.fd = None


def _run(cmd, cwd, log, env=None):
	p = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=900, env=env)
	with open(log, "ab") as f:
		f.write(("$ %s\n" % " ".join(shlex.quote(c) for c in cmd)).encode())
		f.write(p.stdout)
	if p.returncode != 0:
		tail = p.stdout.decode("utf-8", "replace")[-4000:]
		raise BuildError("command failed (exit %d): %s\n%s" % (p.returncode, " ".join(cmd), tail))


def _copy_tree(root, rels, dst):
	for rel in rels:
		d = os.path.join(dst, rel)
		os.makedirs(os.path.dirname(d), exist_ok=True)
		shutil.copy2(os.path.join(root, rel), d)


def build_dir(root=None, cc=None):
	"""Build minimap2 and libminimap2.a out of tree; return the build directory.

	The directory holds a copy of the sources plus `minimap2`,
	`libminimap2.a`, the objects and `build.log`. It is reused while the
	sources, the compiler and $MINIMAP2_TEST_MAKE_ARGS stay the same.
	"""
	root = root or src_dir()
	cc = cc or compiler()
	key = build_key(root, cc)
	cache = cache_root()
	final = os.path.join(cache, "build-" + key)
	marker = os.path.join(final, ".complete")
	if os.path.exists(marker):
		return final
	with _Lock(os.path.join(cache, "build-" + key + ".lock")):
		if os.path.exists(marker):
			return final
		if os.path.exists(final):
			shutil.rmtree(final)
		tmp = tempfile.mkdtemp(prefix="tmp-build-" + key + "-", dir=cache)
		try:
			_copy_tree(root, source_files(root), tmp)
			if build_key(tmp, cc) != key:
				raise BuildError("sources changed while being copied from %s; retry" % root)
			log = os.path.join(tmp, "build.log")
			_run(["make", "-j%d" % make_jobs(), "CC=" + cc] + make_args() + ["minimap2", "libminimap2.a"], tmp, log)
			for prod in ("minimap2", "libminimap2.a"):
				if not os.path.exists(os.path.join(tmp, prod)):
					raise BuildError("make did not produce %s in %s" % (prod, tmp))
			with open(os.path.join(tmp, "build.info"), "w") as f:
				f.write("source: %s\nkey: %s\nmake arguments: %s\ncompiler: %s\n" % (root, key, " ".join(make_args()), compiler_id(cc)))
			open(os.path.join(tmp, ".complete"), "w").close()
			os.replace(tmp, final)
		except BaseException:
			shutil.rmtree(tmp, ignore_errors=True)
			raise
	return final


def extra_dir(root=None, cc=None):
	"""Run `make extra` (sdust, minimap2-lite) in the out-of-tree build; return the directory."""
	d = build_dir(root, cc)
	marker = os.path.join(d, ".extra-complete")
	if os.path.exists(marker):
		return d
	with _Lock(os.path.join(d, ".extra.lock")):
		if os.path.exists(marker):
			return d
		cc = cc or compiler()
		_run(["make", "-j%d" % make_jobs(), "CC=" + cc] + make_args() + ["extra"], d, os.path.join(d, "build.log"))
		for prod in ("sdust", "minimap2-lite"):
			if not os.path.exists(os.path.join(d, prod)):
				raise BuildError("make extra did not produce %s in %s" % (prod, d))
		open(marker, "w").close()
	return d


def minimap2_bin():
	"""Path of the binary under test: $MINIMAP2_BIN, else an out-of-tree build."""
	b = os.environ.get("MINIMAP2_BIN")
	if b:
		b = os.path.abspath(b)
		if not (os.path.isfile(b) and os.access(b, os.X_OK)):
			raise BuildError("MINIMAP2_BIN=%s is not an executable file" % b)
		return b
	return os.path.join(build_dir(), "minimap2")


def sdust_bin():
	"""The `sdust` program of `make extra`, built from $MINIMAP2_TEST_SRC_DIR."""
	return os.path.join(extra_dir(), "sdust")


def lite_bin():
	"""The `minimap2-lite` program (example.c) of `make extra`, built from $MINIMAP2_TEST_SRC_DIR."""
	return os.path.join(extra_dir(), "minimap2-lite")


def libminimap2():
	"""(path to libminimap2.a, include directory), both from the out-of-tree build."""
	d = build_dir()
	return os.path.join(d, "libminimap2.a"), d


# ---------------------------------------------------------------------------
# C drivers
# ---------------------------------------------------------------------------


def c_driver_sources(name):
	src = os.path.join(C_DIR, name + ".c")
	if not os.path.exists(src):
		raise BuildError("no C driver %s" % src)
	return [src, os.path.join(C_DIR, "minimap2_test.h")]


def c_driver(name, extra_cflags=()):
	"""Compile tests/c/<name>.c against libminimap2.a; return the executable's path.

	The driver sees minimap2's headers (minimap.h and the internal ones such
	as mmpriv.h, bseq.h, ksw2.h, sdust.h, kalloc.h) and tests/c/minimap2_test.h.
	It is rebuilt when the driver, the test header, the library build or the
	compiler change.
	"""
	lib, inc = libminimap2()
	cc = compiler()
	extra_cflags = list(extra_cflags) + sanitizer_flags()
	h = hashlib.sha256()
	h.update(os.path.basename(inc).encode())
	h.update(" ".join(extra_cflags).encode())
	for p in c_driver_sources(name):
		with open(p, "rb") as f:
			h.update(hashlib.sha256(f.read()).digest())
	out_dir = os.path.join(inc, "drivers")
	os.makedirs(out_dir, exist_ok=True)
	exe = os.path.join(out_dir, "%s-%s" % (name, h.hexdigest()[:16]))
	if os.path.exists(exe):
		return exe
	with _Lock(exe + ".lock"):
		if os.path.exists(exe):
			return exe
		tmp = exe + ".tmp%d" % os.getpid()
		cmd = shlex.split(cc) + ["-std=gnu99", "-g", "-O2", "-Wall"] + list(extra_cflags) + \
			["-I" + inc, "-I" + C_DIR, "-o", tmp, os.path.join(C_DIR, name + ".c"), lib, "-lz", "-lm", "-lpthread"]
		p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=300)
		if p.returncode != 0:
			raise BuildError("failed to compile C driver %s:\n$ %s\n%s" % (name, " ".join(cmd), p.stdout.decode("utf-8", "replace")))
		os.replace(tmp, exe)
	return exe


# ---------------------------------------------------------------------------
# mappy
# ---------------------------------------------------------------------------


def _have_module(name):
	import importlib.util
	try:
		return importlib.util.find_spec(name) is not None
	except (ImportError, ValueError):
		return False


def mappy_missing():
	"""Why mappy cannot be built with this interpreter, or None if it can."""
	missing = [m for m in ("Cython", "setuptools") if not _have_module(m)]
	if missing:
		return ("%s not importable by %s; run the suite with e.g. "
			"`uv run --no-project --with pytest --with cython --with setuptools python3 -m pytest tests`") % (
			" and ".join(missing), sys.executable)
	cc = shlex.split(compiler())[0]
	if shutil.which(cc) is None:
		return "C compiler %r not found" % cc
	return None


def mappy_key(root=None, cc=None):
	root = root or src_dir()
	h = hashlib.sha256()
	_hash_files(h, root, mappy_source_files(root))
	h.update(compiler_id(cc).encode())
	# the base interpreter, not sys.executable: `uv run --with` makes a new virtualenv path each time
	h.update(("%s\0%s\0%s" % (sys.version, sysconfig.get_config_var("EXT_SUFFIX"), os.path.realpath(sys.base_prefix))).encode())
	for m in ("Cython", "setuptools"):
		try:
			mod = __import__(m)
			h.update(("%s=%s" % (m, getattr(mod, "__version__", "?"))).encode())
		except ImportError:
			pass
	return h.hexdigest()[:20]


def mappy_dir(root=None, cc=None):
	"""Build the mappy extension with setup.py out of tree; return the directory to put on sys.path.

	Raises MappyUnavailable when Cython, setuptools or the compiler is
	missing (tests skip), and BuildError when the build itself fails (tests
	fail). The build uses the running interpreter and CC.
	"""
	why = mappy_missing()
	if why:
		raise MappyUnavailable(why)
	root = root or src_dir()
	cc = cc or compiler()
	key = mappy_key(root, cc)
	cache = cache_root()
	final = os.path.join(cache, "mappy-" + key)
	lib = os.path.join(final, "lib")
	marker = os.path.join(final, ".complete")
	if os.path.exists(marker):
		return lib
	with _Lock(os.path.join(cache, "mappy-" + key + ".lock")):
		if os.path.exists(marker):
			return lib
		if os.path.exists(final):
			shutil.rmtree(final)
		tmp = tempfile.mkdtemp(prefix="tmp-mappy-" + key + "-", dir=cache)
		try:
			src = os.path.join(tmp, "src")
			_copy_tree(root, mappy_source_files(root), src)
			env = dict(os.environ)
			env["CC"] = cc
			log = os.path.join(tmp, "build.log")
			_run([sys.executable, "setup.py", "build_ext", "-j", str(make_jobs()), "--build-lib", os.path.join(tmp, "lib"),
				"--build-temp", os.path.join(tmp, "obj")], src, log, env=env)
			ext = sysconfig.get_config_var("EXT_SUFFIX") or ".so"
			if not os.path.exists(os.path.join(tmp, "lib", "mappy" + ext)):
				raise BuildError("setup.py build_ext did not produce mappy%s in %s" % (ext, tmp))
			with open(os.path.join(tmp, "build.info"), "w") as f:
				f.write("source: %s\nkey: %s\npython: %s\ncompiler: %s\n" % (root, key, sys.executable, compiler_id(cc)))
			open(os.path.join(tmp, ".complete"), "w").close()
			os.replace(tmp, final)
		except BaseException:
			shutil.rmtree(tmp, ignore_errors=True)
			raise
	return lib


def import_mappy():
	"""Build (if needed) and import mappy from the cache; returns the module."""
	import importlib.util
	if "mappy" in sys.modules:
		return sys.modules["mappy"]
	lib = mappy_dir()
	ext = sysconfig.get_config_var("EXT_SUFFIX") or ".so"
	path = os.path.join(lib, "mappy" + ext)
	spec = importlib.util.spec_from_file_location("mappy", path)
	mod = importlib.util.module_from_spec(spec)
	spec.loader.exec_module(mod)
	sys.modules["mappy"] = mod
	return mod


# ---------------------------------------------------------------------------
# Source-tree hygiene
# ---------------------------------------------------------------------------


def tree_products(root=None):
	"""Build products present in a source tree, as {relpath: mtime}."""
	import glob
	root = root or src_dir()
	out = {}
	for pat in TREE_PRODUCT_PATTERNS:
		for p in glob.glob(os.path.join(root, pat)):
			out[os.path.relpath(p, root)] = os.stat(p).st_mtime
	return out


if __name__ == "__main__":
	print(minimap2_bin())

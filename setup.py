import configparser
import glob
import mock
import os
import platform
import re
import subprocess
import sys

import setuptools
from distutils import log
from distutils.command.clean import clean as _clean
from setuptools.command.build_ext import build_ext as _build_ext
from setuptools.command.build_clib import build_clib as _build_clib
from setuptools.command.sdist import sdist as _sdist
from setuptools.extension import Extension, Library

try:
    from Cython.Build import cythonize
except ImportError as err:
    cythonize = err

# --- Constants -----------------------------------------------------------------

SYSTEM  = platform.system()
MACHINE = platform.machine()
if re.match("^mips", MACHINE):
    TARGET_CPU = "mips"
elif re.match("^arm", MACHINE):
    TARGET_CPU = "arm"
elif re.match("^aarch64", MACHINE):
    TARGET_CPU = "aarch64"
elif re.match("(x86_64)|(AMD64|amd64)|(^i.86$)", MACHINE):
    TARGET_CPU = "x86"
elif re.match("^(powerpc|ppc)", MACHINE):
    TARGET_CPU = "ppc"

# --- Utils ------------------------------------------------------------------

def _eprint(*args, **kwargs):
    print(*args, **kwargs, file=sys.stderr)

# --- Commands ------------------------------------------------------------------

class sdist(_sdist):
    """A `sdist` that generates a `pyproject.toml` on the fly.
    """

    def run(self):
        # build `pyproject.toml` from `setup.cfg`
        c = configparser.ConfigParser()
        c.add_section("build-system")
        c.set("build-system", "requires", str(self.distribution.setup_requires))
        c.set("build-system", 'build-backend', '"setuptools.build_meta"')
        with open("pyproject.toml", "w") as pyproject:
            c.write(pyproject)
        # run the rest of the packaging
        _sdist.run(self)


class build_ext(_build_ext):
    """A `build_ext` that disables optimizations if compiled in debug mode.
    """

    def finalize_options(self):
        _build_ext.finalize_options(self)
        self._clib_cmd = self.get_finalized_command("build_clib")
        self._clib_cmd.debug = self.debug
        self._clib_cmd.force = self.force

    def build_extension(self, ext):
        # show the compiler being used
        _eprint("building", ext.name, "with", self.compiler.compiler_type, "compiler")

        # add debug flags if we are building in debug mode
        if self.debug:
            if self.compiler.compiler_type in {"unix", "cygwin", "mingw32"}:
                ext.extra_compile_args.append("-O0")
            elif self.compiler.compiler_type == "msvc":
                ext.extra_compile_args.append("/Od")
            if sys.implementation.name == "cpython":
                ext.define_macros.append(("CYTHON_TRACE_NOGIL", 1))
        else:
            ext.define_macros.append(("CYTHON_WITHOUT_ASSERTIONS", 1))

        # update link and include directories
        for name in ext.libraries:
            lib = self._clib_cmd.get_library(name)
            libfile = self.compiler.library_filename(
                lib.name, output_dir=self._clib_cmd.build_clib
            )
            ext.depends.append(libfile)
            ext.extra_objects.append(libfile)

        # build the rest of the extension as normal
        _build_ext.build_extension(self, ext)

    def run(self):
        # check `cythonize` is available
        if isinstance(cythonize, ImportError):
            raise RuntimeError("Cython is required to run `build_ext` command") from cythonize

        # use debug directives with Cython if building in debug mode
        cython_args = {
            "include_path": ["include"],
            "compiler_directives": {},
            "compile_time_env": {
                "SYS_IMPLEMENTATION_NAME": sys.implementation.name,
                "SYS_VERSION_INFO_MAJOR": sys.version_info.major,
                "SYS_VERSION_INFO_MINOR": sys.version_info.minor,
                "SYS_VERSION_INFO_MICRO": sys.version_info.micro,
                "TARGET_CPU": TARGET_CPU,
                "AVX2_BUILD_SUPPORT": False,
                "SSE_BUILD_SUPPORT": False,
            }
        }
        if self.force:
            cython_args["force"] = True
        if self.debug:
            cython_args["annotate"] = True
            cython_args["compiler_directives"]["warn.undeclared"] = True
            cython_args["compiler_directives"]["warn.unreachable"] = True
            cython_args["compiler_directives"]["warn.maybe_uninitialized"] = True
            cython_args["compiler_directives"]["warn.unused"] = True
            cython_args["compiler_directives"]["warn.unused_arg"] = True
            cython_args["compiler_directives"]["warn.unused_result"] = True
            cython_args["compiler_directives"]["warn.multiple_declarators"] = True
        else:
            cython_args["compiler_directives"]["boundscheck"] = False
            cython_args["compiler_directives"]["wraparound"] = False
            cython_args["compiler_directives"]["cdivision"] = True

        # compile the C library
        if not self.distribution.have_run.get("build_clib", False):
            self._clib_cmd.run()

        # check if we can build platform-specific code
        if self._clib_cmd._avx2_supported:
            cython_args["compile_time_env"]["AVX2_BUILD_SUPPORT"] = True
            for ext in self.extensions:
                ext.extra_compile_args.append(self._clib_cmd._avx2_flag())
        if self._clib_cmd._sse2_supported:
            cython_args["compile_time_env"]["SSE2_BUILD_SUPPORT"] = True
            for ext in self.extensions:
                ext.extra_compile_args.append(self._clib_cmd._sse2_flag())

        # cythonize the extensions
        self.extensions = cythonize(self.extensions, **cython_args)

        # patch the extensions (somehow needed for `_build_ext.run` to work)
        for ext in self.extensions:
            ext._needs_stub = False

        # build the extensions as normal
        _build_ext.run(self)


class build_clib(_build_clib):
    """A custom `build_clib` that splits the `training.c` file from Prodigal.
    """

    # --- Autotools-like helpers ---

    def _silent_spawn(self, cmd):
        try:
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as err:
            raise CompileError(err.stderr)

    def _check_simd_generic(self, name, flag, header, vector, set, extract):
        _eprint('checking whether compiler can build', name, 'code', end="... ")

        testfile = os.path.join(self.build_temp, "have_{}.c".format(name))
        binfile = os.path.join(self.build_temp, "have_{}.bin".format(name))
        objects = []

        with open(testfile, "w") as f:
            f.write("""
                #include <{}>
                int main() {{
                    {}      a = {}(1);
                    short   x = {}(a, 1);
                    return (x == 1) ? 0 : 1;
                }}
            """.format(header, vector, set, extract))
        try:
            with mock.patch.object(self.compiler, "spawn", new=self._silent_spawn):
                objects = self.compiler.compile([testfile], debug=self.debug, extra_preargs=[flag])
                self.compiler.link_executable(objects, binfile)
                subprocess.run([binfile], check=True)
        except CompileError:
            _eprint("no")
            return False
        except subprocess.CalledProcessError:
            _eprint("yes, but cannot run code")
            return True  # assume we are cross-compiling, and still build
        else:
            _eprint("yes, with {}".format(flag))
            return True
        finally:
            os.remove(testfile)
            for obj in filter(os.path.isfile, objects):
                os.remove(obj)
            if os.path.isfile(binfile):
                os.remove(binfile)

    def _avx2_flag(self):
        if self.compiler.compiler_type == "msvc":
            return "/arch:AVX2"
        return "-mavx2"

    def _check_avx2(self):
        return self._check_simd_generic(
            "AVX2",
            self._avx2_flag(),
            header="immintrin.h",
            vector="__m256i",
            set="_mm256_set1_epi16",
            extract="_mm256_extract_epi16",
        )

    def _sse2_flag(self):
        if self.compiler.compiler_type == "msvc":
            return "/arch:SSE2"
        return "-msse2"

    def _check_sse2(self):
        return self._check_simd_generic(
            "SSE2",
            self._sse2_flag(),
            header="emmintrin.h",
            vector="__m128i",
            set="_mm_set1_epi16",
            extract="_mm_extract_epi16",
        )

    # --- Compatibility with base `build_clib` command ---

    def check_library_list(self, libraries):
        pass

    def get_library_names(self):
        return [ lib.name for lib in self.libraries ]

    def get_source_files(self):
        return [ source for lib in self.libraries for source in lib.sources ]

    def get_library(self, name):
        return next(lib for lib in self.libraries if lib.name == name)

    # --- Compatibility with `setuptools.Command`

    def initialize_options(self):
        _build_clib.initialize_options(self)
        self._avx2_supported = False
        self._sse2_supported = False

    def finalize_options(self):
        _build_clib.finalize_options(self)
        # extract the training file and the temporary folder where to split
        # the training profiles
        lib = self.libraries[0]
        self.training_file = next(s for s in lib.sources if os.path.basename(s) == "training.c")
        self.training_temp = os.path.join(self.build_temp, "training")

    # --- Build code ---

    def _write_source_split(self, training_temp, index, lines):
        filename = os.path.join(training_temp, "training{:02}.c".format(index))
        with open(filename, "wb") as dst:
            if index != 0:
                dst.write(b'#include "training.h"\n')
            dst.writelines(lines)
        lines.clear()

    def split_training_source(self, training_file, training_temp):
        self.mkpath(training_temp)
        with open(training_file, "rb") as src:
            lines = []
            index = 0
            for line in src:
                if line.startswith(b"void initialize_metagenome"):
                    self._write_source_split(training_temp, index, lines)
                    index += 1
                if line.lstrip().startswith(b"struct _training"):
                    line = line.replace(b"struct _training", b"static const struct _training")
                lines.append(line)
            self._write_source_split(training_temp, index, lines)

    def run(self):
        # split the huge `training.c` file in small chunks with individual
        # functions so that it can compile even on low-memory machines
        self.make_file(
            [self.training_file],
            self.training_temp,
            self.split_training_source,
            (self.training_file, self.training_temp)
        )
        # build the library as normal
        _build_clib.run(self)
        # check which CPU features are supported, now that
        # `_build_clib.run` has properly initialized the C compiler
        self._avx2_supported = self._check_avx2()
        self._sse2_supported = self._check_sse2()

    def build_libraries(self, libraries):
        self.mkpath(self.build_clib)
        for library in libraries:
            self.make_file(
                library.sources,
                self.compiler.library_filename(library.name, output_dir=self.build_clib),
                self.build_library,
                (library,)
            )

    def build_library(self, library):
        # show the compiler being used
        _eprint("building", library.name, "with", self.compiler.compiler_type, "compiler")

        # add debug flags if we are building in debug mode
        if self.debug:
            if self.compiler.compiler_type in {"unix", "cygwin", "mingw32"}:
                library.extra_compile_args.append("-O0")
            elif self.compiler.compiler_type == "msvc":
                library.extra_compile_args.append("/Od")

        # store compile args
        compile_args = (
            library.define_macros,
            library.include_dirs,
            self.debug,
            library.extra_compile_args,
            None,
            library.depends,
        )

        # compile Prodigal
        if library.name == "prodigal":
            # compile sources
            sources_lib = library.sources.copy()
            sources_lib.remove(self.training_file)
            objects_lib = [
                os.path.join(self.build_temp, s.replace(".c", self.compiler.obj_extension))
                for s in sources_lib
            ]
            for source, object in zip(sources_lib, objects_lib):
                self.make_file(
                    [source],
                    object,
                    self.compiler.compile,
                    ([source], self.build_temp, *compile_args)
                )
            # compile `training.c` source files
            sources_training = sorted(glob.iglob(os.path.join(self.training_temp, "*.c")))
            objects_training = [s.replace(".c", self.compiler.obj_extension) for s in sources_training]
            for source, object in zip(sources_training, objects_training):
                self.make_file(
                    [source],
                    object,
                    self.compiler.compile,
                    ([source], None, *compile_args)
                )

        else:
            # compile sources
            sources_lib = library.sources.copy()
            objects_lib = [
                os.path.join(self.build_temp, s.replace(".c", self.compiler.obj_extension))
                for s in sources_lib
            ]
            objects_training = []
            for source, object in zip(sources_lib, objects_lib):
                self.make_file(
                    [source],
                    object,
                    self.compiler.compile,
                    ([source], self.build_temp, *compile_args)
                )

        # link into a static library
        libfile = self.compiler.library_filename(
            library.name,
            output_dir=self.build_clib,
        )
        self.make_file(
            objects_lib + objects_training,
            libfile,
            self.compiler.create_static_lib,
            (objects_lib + objects_training, library.name, self.build_clib, None, self.debug)
        )


class clean(_clean):
    """A `clean` that removes intermediate files created by Cython.
    """

    def run(self):

        source_dir = os.path.join(os.path.dirname(__file__), "pyrodigal")

        patterns = ["*.html"]
        if self.all:
            patterns.extend(["*.so", "*.c"])

        for pattern in patterns:
            for file in glob.glob(os.path.join(source_dir, pattern)):
                log.info("removing {!r}".format(file))
                os.remove(file)

        _clean.run(self)

# --- Setup ---------------------------------------------------------------------

setuptools.setup(
    libraries=[
        Library(
            "prodigal",
            sources=[
                os.path.join("vendor", "Prodigal", "{}.c".format(base))
                for base in [
                    "bitmap",
                    "dprog",
                    "gene",
                    "metagenomic",
                    "node",
                    "sequence",
                    "training"
                ]
            ],
            include_dirs=[os.path.join("vendor", "Prodigal")]
        ),
        Library(
            "cpu_features",
            sources=[
                os.path.join("vendor", "cpu_features", "src", "{}.c".format(base))
                for base in [
                    "cpuinfo_{}".format(TARGET_CPU),
                    "filesystem",
                    "stack_line_reader",
                    "string_view",
                ]
            ],
            include_dirs=[os.path.join("vendor", "cpu_features", "include")],
            define_macros=[("STACK_LINE_READER_BUFFER_SIZE", 1024)]
        )
    ],
    ext_modules=[
        Extension(
            "pyrodigal._pyrodigal",
            sources=[
                "pyrodigal/impl/sse.c",
                "pyrodigal/impl/avx.c",
                "pyrodigal/_pyrodigal.pyx"
            ],
            include_dirs=[
                "pyrodigal",
                os.path.join("vendor", "Prodigal"),
                os.path.join("vendor", "cpu_features", "include"),
            ],
            libraries=[
                "prodigal",
                "cpu_features",
            ],
        ),
    ],
    cmdclass={
        "sdist": sdist,
        "build_ext": build_ext,
        "build_clib": build_clib,
        "clean": clean
    }
)

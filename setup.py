import configparser
import glob
import os
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
        log.info("building {} with {} compiler".format(
            ext.name, self.compiler.compiler_type
        ))

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
            ext.extra_objects.append(self.compiler.library_filename(
                lib.name, output_dir=self._clib_cmd.build_clib
            ))

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

        # cythonize the extensions
        self.extensions = cythonize(self.extensions, **cython_args)

        # patch the extensions (somehow needed for `_build_ext.run` to work)
        for ext in self.extensions:
            ext._needs_stub = False

        # compile the C library
        if not self.distribution.have_run.get("build_clib", False):
            self._clib_cmd.run()

        # build the extensions as normal
        _build_ext.run(self)


class build_clib(_build_clib):
    """A custom `build_clib` that splits the `training.c` file from Prodigal.
    """

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
        # add debug flags if we are building in debug mode
        if self.debug:
            if self.compiler.compiler_type in {"unix", "cygwin", "mingw32"}:
                library.extra_compile_args.append("-O0")
            elif self.compiler.compiler_type == "msvc":
                library.extra_compile_args.append("/Od")

        # compile Prodigal
        sources_lib = library.sources.copy()
        sources_lib.remove(self.training_file)
        objects_lib = self.compiler.compile(
            sources_lib,
            output_dir=self.build_temp,
            include_dirs=library.include_dirs,
            macros=library.define_macros,
            debug=self.debug,
            depends=library.depends,
            extra_preargs=library.extra_compile_args,
        )

        # compile `training.c` source files
        sources_training = sorted(glob.iglob(os.path.join(self.training_temp, "*.c")))
        objects_training = self.compiler.compile(
            sources_training,
            include_dirs=library.include_dirs,
            macros=library.define_macros,
            debug=self.debug,
            depends=library.depends,
            extra_preargs=library.extra_compile_args,
        )

        # link into a static library
        self.compiler.create_static_lib(
            objects_lib + objects_training,
            library.name,
            output_dir=self.build_clib,
            debug=self.debug,
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


setuptools.setup(
    libraries=[
        Library(
            "prodigal",
            sources=[
                os.path.join("Prodigal", base)
                for base in [
                    "bitmap.c",
                    "dprog.c",
                    "gene.c",
                    "metagenomic.c",
                    "node.c",
                    "sequence.c",
                    "training.c"
                ]
            ],
            include_dirs=["Prodigal"],
        )
    ],
    ext_modules=[
        Extension(
            "pyrodigal._pyrodigal",
            sources=["pyrodigal/_pyrodigal.pyx"],
            include_dirs=["Prodigal"],
            libraries=["prodigal"],
        ),
    ],
    cmdclass={
        "sdist": sdist,
        "build_ext": build_ext,
        "build_clib": build_clib,
        "clean": clean
    }
)

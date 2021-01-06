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

    def build_extension(self, ext):
        # show the compiler being used
        log.info("building {} with {} compiler".format(
            ext.name, self.compiler.compiler_type
        ))

        # make sure the C libraries have been built already
        self.run_command("build_clib")
        _clib_cmd = self.get_finalized_command("build_clib")

        # add debug flags if we are building in debug mode
        if self.debug:
            if self.compiler.compiler_type in {"unix", "cygwin", "mingw32"}:
                ext.extra_compile_args.append("-O0")
            elif self.compiler.compiler_type == "msvc":
                ext.extra_compile_args.append("/Od")
            if sys.implementation.name == "cpython":
                ext.define_macros.append(("CYTHON_TRACE_NOGIL", 1))

        # add path to static library as an extra object to make sure linking
        # works on OSX as well
        ext.extra_objects = [
            os.path.join(_clib_cmd.build_clib, self.compiler.library_filename(lib))
            for lib in ext.libraries
        ]

        # build the rest of the extension as normal
        _build_ext.build_extension(self, ext)

    def run(self):
        # check `cythonize` is available
        if isinstance(cythonize, ImportError):
            raise RuntimeError("Cython is required to run `build_ext` command") from cythonize

        # use debug directives with Cython if building in debug mode
        cython_args = {"include_path": ["include"], "compiler_directives": {}}
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

        # patch the extensions (needed for `_build_ext.run` to work)
        for ext in self.extensions:
            ext._needs_stub = False

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
        lib = self.libraries[0]
        training_file = next(s for s in lib.sources if os.path.basename(s) == "training.c")
        training_temp = os.path.join(self.build_temp, "training")
        self.make_file([training_file], training_temp, self.split_training_source, (training_file, training_temp))

        # patch the library sources to use the split training files
        lib.sources.remove(training_file)
        lib.sources.extend(sorted(glob.glob(os.path.join(training_temp, "*.c"))))

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

        # compile and link as usual
        objects = self.compiler.compile(
            library.sources,
            output_dir=self.build_temp,
            include_dirs=library.include_dirs + [self.build_clib],
            macros=library.define_macros,
            debug=self.debug,
            depends=library.depends,
            extra_preargs=library.extra_compile_args,
        )
        self.compiler.create_static_lib(
            objects,
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
            sources=["pyrodigal/__init__.pyx"],
            include_dirs=["Prodigal"],
            libraries=["prodigal"],
        )
    ],
    cmdclass={
        "sdist": sdist,
        "build_ext": build_ext,
        "build_clib": build_clib,
        "clean": clean
    }
)

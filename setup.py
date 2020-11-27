import configparser
import glob
import os
import sys

import setuptools
from Cython.Build import cythonize
from setuptools.command.build_ext import build_ext as _build_ext
from setuptools.command.build_clib import build_clib as _build_clib
from setuptools.command.sdist import sdist as _sdist
from setuptools.extension import Extension, Library


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
        if self.debug:
            if sys.platform == "linux" or sys.platform == "darwin":
                ext.extra_compile_args.append("-O0")
            if sys.implementation.name == "cpython":
                ext.define_macros.append(("CYTHON_TRACE_NOGIL", 1))
        _build_ext.build_extension(self, ext)

    def run(self):
        self.run_command("build_clib")
        _build_ext.run(self)

class build_clib(_build_clib):
    """A custom `build_clib` that compiles out of source.
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

libraries = [
    Library(
        "prodigal",
        sources=glob.glob(os.path.join("Prodigal", "*.c")),
        include_dirs=["Prodigal"],
    )
]

extensions = [
    Extension(
        "pyrodigal._pyrodigal",
        sources=["pyrodigal/__init__.pyx"],
        include_dirs=["Prodigal"],
        libraries=["prodigal"],
    ),
]


setuptools.setup(
    libraries=libraries,
    ext_modules=cythonize(extensions, annotate=True),
    cmdclass=dict(sdist=sdist, build_ext=build_ext, build_clib=build_clib),
)

import configparser
import glob
import os
import sys

import setuptools
from Cython.Build import cythonize
from setuptools.command.build_ext import build_ext as _build_ext
from setuptools.command.sdist import sdist as _sdist


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

    def build_extension(self, ext):
        training_file = ext.sources.pop()
        training_temp = os.path.join(self.build_temp, "training")
        self.make_file([training_file], training_temp, self.split_training_source, (training_file, training_temp))
        ext.sources.extend(sorted(glob.glob(os.path.join(training_temp, "*.c"))))

        if self.debug:
            if sys.platform == "linux" or sys.platform == "darwin":
                ext.extra_compile_args.append("-O0")
            if sys.implementation.name == "cpython":
                ext.define_macros.append(("CYTHON_TRACE_NOGIL", 1))
        _build_ext.build_extension(self, ext)


extensions = [
    setuptools.Extension(
        "pyrodigal._pyrodigal",
        [
            "pyrodigal/__init__.pyx",
            "Prodigal/bitmap.c",
            "Prodigal/dprog.c",
            "Prodigal/gene.c",
            "Prodigal/metagenomic.c",
            "Prodigal/node.c",
            "Prodigal/sequence.c",
            # leave this last
            "Prodigal/training.c"
        ],
        include_dirs=["Prodigal"],
        libraries=[],
        library_dirs=[],
        extra_compile_args=["-Wall"],
    ),
]


setuptools.setup(
    ext_modules=cythonize(extensions, annotate=True),
    cmdclass=dict(sdist=sdist, build_ext=build_ext),
)

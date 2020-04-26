import configparser
import sys

import setuptools
from Cython.Build import cythonize
from setuptools.command.sdist import sdist as _sdist


class sdist(_sdist):
    """A custom `sdist` that generates a `pyproject.toml` on the fly.
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


if "--debug" in sys.argv:
    cflags = ["-pedantic", "-Wall", "-O0"]
else:
    cflags = ["-pedantic", "-Wall"]

extensions = [
    setuptools.Extension(
        "pyrodigal",
        [
            "pyrodigal/__init__.pyx",
            "Prodigal/bitmap.c",
            "Prodigal/dprog.c",
            "Prodigal/gene.c",
            "Prodigal/metagenomic.c",
            "Prodigal/node.c",
            "Prodigal/sequence.c",
            "Prodigal/training.c"
        ],
        include_dirs=["Prodigal"],
        libraries=[],
        library_dirs=[],
        extra_compile_args=cflags,
    ),
]


setuptools.setup(
    ext_modules=cythonize(extensions),
    cmdclass=dict(sdist=sdist),
)

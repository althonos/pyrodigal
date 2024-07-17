Installation
============

.. note::

    Wheels are provided for Linux, MacOS and Windows x86-64, as well
    as Linux and MacOS Aarch64, but other machines will have to build the wheel 
    from the source distribution. Building ``pyrodigal`` involves compiling 
    Prodigal, which requires a C compiler to be available.


PyPi
^^^^

``pyrodigal`` is hosted on GitHub, but the easiest way to install it is to download
the latest release from its `PyPi repository <https://pypi.python.org/pypi/pyrodigal>`_.
It will install all dependencies then install ``pyrodigal`` either from a wheel if
one is available, or from source after compiling the Cython code :

.. code:: console

   $ pip install --user pyrodigal

Conda
^^^^^

Pronto is also available as a `recipe <https://anaconda.org/bioconda/pyrodigal>`_
in the `bioconda <https://bioconda.github.io/>`_ channel. To install, simply
use the ``conda`` installer:

.. code:: console

   $ conda install -c bioconda pyrodigal


Arch User Repository
^^^^^^^^^^^^^^^^^^^^

A package recipe for Arch Linux can be found in the Arch User Repository
under the name `python-pyrodigal <https://aur.archlinux.org/packages/python-pyrodigal>`_.
It will always match the latest release from PyPI.

Steps to install on ArchLinux depend on your `AUR helper <https://wiki.archlinux.org/title/AUR_helpers>`_
(``yaourt``, ``aura``, ``yay``, etc.). For ``aura``, you'll need to run:

.. code:: console

    $ aura -A python-pyrodigal


BioArchLinux
^^^^^^^^^^^^

The `BioArchLinux <https://bioarchlinux.org>`_ project provides pre-compiled packages
based on the AUR recipe. Add the BioArchLinux package repository to ``/etc/pacman.conf``:

.. code:: ini

    [bioarchlinux]
    Server = https://repo.bioarchlinux.org/$arch

Then install the latest version of the package and its dependencies with ``pacman``:

.. code:: console

    $ pacman -S python-pyrodigal


Piwheels
^^^^^^^^

Pyrodigal works on Raspberry Pi computers (with NEON vectorization enabled!), 
and pre-built wheels are compiled for `armv7l` platforms on piwheels.
Run the following command to install these instead of compiling from source:

.. code:: console

   $ pip3 install pyrodigal --extra-index-url https://www.piwheels.org/simple

Check the `piwheels documentation <https://www.piwheels.org/faq.html>`_ for 
more information.


FreeBSD Port
^^^^^^^^^^^^

On BSD systems, Pyrodigal is available as a *port* and can be installed as 
a compiled package with one of the following commands:

.. code:: console

    $ pkg install biology/py-pyrodigal
    $ pkg install py39-pyrodigal

See the `FreshPorts page <https://www.freshports.org/biology/py-pyrodigal/>`_
for more information.


GitHub + ``pip``
^^^^^^^^^^^^^^^^

If, for any reason, you prefer to download the library from GitHub, you can clone
the repository and install the repository by running (with the admin rights):

.. code:: console

   $ git clone --recursive https://github.com/althonos/pyrodigal
   $ pip install --user ./pyrodigal

.. caution::

    Keep in mind this will install always try to install the latest commit,
    which may not even build, so consider using a versioned release instead.


GitHub + ``setuptools``
^^^^^^^^^^^^^^^^^^^^^^^

If you do not want to use ``pip``, you can still clone the repository and
run the ``setup.py`` file manually, although you will need to install the
build dependencies (mainly `Cython <https://pypi.org/project/cython>`_):

.. code:: console

   $ git clone --recursive https://github.com/althonos/pyrodigal
   $ cd pyrodigal
   $ python setup.py build_ext
   # python setup.py install

.. Danger::

    Installing packages without ``pip`` is strongly discouraged, as they can
    only be uninstalled manually, and may damage your system.

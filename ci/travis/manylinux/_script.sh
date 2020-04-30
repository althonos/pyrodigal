#!/bin/sh -e

export PATH="$HOME/.cargo/bin:$PATH"
export PYBIN="$(echo /opt/python/${1}*/bin)"
export PYTHON_SYS_EXECUTABLE="$PYBIN/python"
export PYTHON_LIB=$(${PYBIN}/python -c "import sysconfig; print(sysconfig.get_config_var('LIBDIR'))")
export LIBRARY_PATH="$LIBRARY_PATH:$PYTHON_LIB"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$PYTHON_LIB"

cd /io
$PYTHON_SYS_EXECUTABLE setup.py build_ext --inplace --debug
$PYTHON_SYS_EXECUTABLE -m coverage run -p -m unittest discover -v
$PYTHON_SYS_EXECUTABLE -m coverage combine
$PYTHON_SYS_EXECUTABLE -m coverage xml
$PYTHON_SYS_EXECUTABLE -m coverage report

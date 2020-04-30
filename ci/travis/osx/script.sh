#!/bin/sh -e

. $(dirname $(dirname $0))/functions.sh

# --- Using proper Python executable -----------------------------------------

log Activating pyenv
eval "$(pyenv init -)"
pyenv shell $(pyenv versions --bare)

# --- Test -------------------------------------------------------------------

$PYTHON setup.py build_ext --inplace --debug
$PYTHON setup.py test

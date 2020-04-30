#!/bin/sh -e

. $(dirname $(dirname $0))/functions.sh

# --- Install Python ---------------------------------------------------------

if [ ! -e "$PYENV_ROOT/.git" ]; then
  log Installing pyenv from GitHub
  git clone https://github.com/pyenv/pyenv.git $PYENV_ROOT
else
  log Updating pyenv
  cd "$PYENV_ROOT"
  git pull
  cd "$TRAVIS_BUILD_DIR"
fi

log Activating pyenv
eval "$(pyenv init -)"

log Installing Python v${PYTHON#python}
pyenv install --skip-existing ${PYTHON#python}-dev
pyenv shell ${PYTHON#python}-dev


# --- Check system Python version --------------------------------------------

log Using $(python --version | head -n1 | cut -d' ' -f1,2)


# --- Install Python requirements --------------------------------------------

log Installing Python requirements
$PYTHON -m pip install -U -r "$TRAVIS_BUILD_DIR/ci/requirements.txt"

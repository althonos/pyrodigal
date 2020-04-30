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

if [ "$PYTHON" = "pypy3.7" ]; then
  log Installing PyPy v3.7
  brew unlink python
  brew install pypy3
  rm -f /usr/local/bin/python /usr/local/bin/python3
  ln -s /usr/local/bin/pypy3 /usr/local/bin/python
  ln -s /usr/local/bin/pypy3 /usr/local/bin/python3
elif [ "$PYTHON" = "pypy3.6" ]; then
  log Installing PyPy v3.6
  pyenv install --skip-existing pypy3.6-7.3.0
  pyenv shell pypy3.6-7.3.0
elif [ "$PYTHON" = "pypy3.5" ]; then
  log Installing PyPy v3.5
  pyenv install --skip-existing pypy3.5-7.0.0
  pyenv shell pypy3.5-7.0.0
else
  log Installing Python v${PYTHON#python}
  pyenv install --skip-existing ${PYTHON#python}-dev
  pyenv shell ${PYTHON#python}-dev
fi


# --- Check system Python version --------------------------------------------

log Using $(python --version | head -n1 | cut -d' ' -f1,2)


# --- Install Python requirements --------------------------------------------

log Installing Python requirements
$PYTHON -m pip install -r "$TRAVIS_BUILD_DIR/ci/requirements.txt"

#!/bin/sh -e

. $(dirname $(dirname $0))/functions.sh

# --- Install Python dependencies --------------------------------------------

case $TRAVIS_PYTHON_VERSION in
  pypy3) TAG=pp36-pypy3_*;;
  *)     TAG=cp$(echo $TRAVIS_PYTHON_VERSION | sed 's/\.//');;
esac

log Installing Python dependencies in \`manylinux\` container
docker exec -it manylinux sh /io/ci/travis/manylinux/_install.sh $TAG


# --- Install Python deployment dependencies ---------------------------------

log Installing Python requirements
pip install -r ci/requirements.txt

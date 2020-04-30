#!/bin/sh -e

. $(dirname $(dirname $0))/functions.sh


# --- Build wheel ------------------------------------------------------------

case $TRAVIS_PYTHON_VERSION in
  pypy3) TAG=pp371-pypy3_71;;
  *)     TAG=cp$(echo $TRAVIS_PYTHON_VERSION | sed 's/\.//');;
esac

log Building wheel with $TAG
docker exec -it manylinux sh /io/ci/travis/manylinux/_before_deploy.sh $TAG


# --- Check distribution files -----------------------------------------------

log Checking distribution archives
twine check dist/*

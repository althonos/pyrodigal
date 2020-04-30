#!/bin/sh -e

. $(dirname $(dirname $0))/functions.sh

# --- Coverage ---------------------------------------------------------------

log Uploading coverage data to Codecov
codecov

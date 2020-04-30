#!/bin/sh -e

. $(dirname $(dirname $0))/functions.sh


# --- Cleanup Homebrew cache -------------------------------------------------

brew cleanup

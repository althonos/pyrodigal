#!/bin/sh

. $(dirname $0)/functions.sh

# --- Update GitHub release notes --------------------------------------------

export GEM_PATH="$(ruby -r rubygems -e 'puts Gem.user_dir')"
export PATH="${GEM_PATH}/bin:$PATH"

log Installing chandler gem
gem install --user-install chandler

log Updating GitHub release notes
chandler push --github="$TRAVIS_REPO_SLUG" --changelog="CHANGELOG.md"

# --- Deploy to PyPI ---------------------------------------------------------

log Deploying to PyPI
python3 -m twine upload --skip-existing dist/*.whl dist/*.tar.gz

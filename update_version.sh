#!/usr/bin/bash -l 

#==============================================================================
# BASH Strict mode (i.e. "fail fast" to reduce hard-to-find bugs)
set -e          # EXIT the script if any command returns non-zero exit status.
set -E          # Make ERR trapping work inside functions too.
set -u          # Variables must be pre-defined before using them.
set -o pipefail # If a pipe fails, returns the error code for the failed pipe
                #  even if it isn't the last command in a series of pipes.

message=${1?Commit message required before updating version}

# Date version: Year, month, day, hour, minute
new_version="v$(date +'%Y-%m-%d-%H-%M')"

sed -z "s/# VERSION\n\nv[0-9]\+-[0-9]\+-[0-9]\+-[0-9]\+-[0-9]\+/# VERSION\n\n$new_version/" -i README.md

echo "$new_version" > version.txt

git add README.md
git add version.txt

git commit -m "$message ($new_version)" && git push

git tag -a $new_version -m "$message"

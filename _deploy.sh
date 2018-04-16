#!/bin/sh

set -e

# [ -z "${GITHUB_PAT}" ] && exit 0
# [ "${TRAVIS_BRANCH}" != "master" ] && exit 0
#
# git config --global user.email "martins@gmail.com"
# git config --global user.name "Martin Smith"
#
# git clone -b gh-pages https://${GITHUB_PAT}@github.com/${TRAVIS_REPO_SLUG}.git book-output
# cd book-output
# cp -r ../_book/* ./
# git add --all *
# git commit -m"Update compiled files" || true
# git push -q origin gh-pages

#!/bin/bash

# this is used to stop the script if something fails
set -e

git checkout gh-pages
cp -r docs/build/html/* .
cp -r docs/doxygen/html doxygen/
touch .nojekyll

git add .nojekyll
git add *.html
git add searchindex.js
git add -f _images
git add -f _static
git add -f _sources
git add -f doxygen
git add -f userguide
git add -f devguide
git add -f methods

git commit -am 'update the documentation'
git push origin gh-pages

#git checkout master



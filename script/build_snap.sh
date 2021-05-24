#!/bin/sh
# Downloads, builds and puts in the project tree the Snap library
# version 6.0.
#
# For more info visit: https://snap.stanford.edu/

VERSION="6.0"

cd ../lib
if [ ! -d "snap" ]; then
  wget https://snap.stanford.edu/releases/Snap-"$VERSION".zip
  unzip Snap-"$VERSION".zip && rm Snap-"$VERSION".zip
  mv Snap-"$VERSION" snap
  cd snap
  make --silent
  rm -rf contrib doxygen examples snap-adv snap-exp test tutorials Credits.txt ReadMe.txt Makefile.config Release.txt Makefile Makefile.config .gitignore
fi

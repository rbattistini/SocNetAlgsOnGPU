#!/bin/sh
# Downloads, builds and puts in the project tree the Snap library
#
# For more info visit: http://snap-graph.sourceforge.net/

VERSION="0.4"

cd ../lib
if [ ! -d "snap" ]; then
  if [ ! -f "snap-"$VERSION".tar.gz" ]; then
    wget https://master.dl.sourceforge.net/project/snap-graph/snap-graph/snap-"$VERSION"/snap-"$VERSION".tar.gz
  fi
  tar -xvzf snap-"$VERSION".tar.gz
  mkdir snap
  mkdir build
  cd build
  ../snap-"$VERSION"/configure --prefix="`cd ../snap; pwd`" --enable-openmp
  make
  make install
  cd ../snap
  rm -rf ../build/ bin/ ../snap-"$VERSION"/ ../snap-"$VERSION".tar.gz
fi

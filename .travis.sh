#! /bin/bash

set -e

# install xdrfile library first
wget ftp://ftp.gromacs.org/pub/contrib/xdrfile-1.1.4.tar.gz
tar xzf xdrfile-1.1.4.tar.gz
cd xdrfile-1.1.4
./configure CC=$PLUMED_CC CXX=$PLUMED_CXX
make
sudo make install
cd ../

# these packages take a lot of time and are needed only when
# building documentation
if [ "$MAKEDOC" = yes ] ; then
  sudo apt-get install -y doxygen
  sudo apt-get install -y graphviz
  sudo apt-get install -y texlive-fonts-extra
  sudo apt-get install -y texlive-full
fi


#! /bin/bash

test ! -d "$PLUMED1/common_files" && {
  echo "Cannot find $PLUMED1/common_files"
  echo 'Please set PLUMED1 environment variable'
  exit
}

for file in *.cpp *.h
do
  case "$file" in
  (*Plumed1.cpp) ;;
  (*) rm $file ;;
  esac
done

cd $PLUMED1/common_files/
files="$(echo *.c)"
cd -

for file in $files
do
  ln -s $PLUMED1/common_files/$file ${file}pp
done

ln -s $PLUMED1/common_files/metadyn.h


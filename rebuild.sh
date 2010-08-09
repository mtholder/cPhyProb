#!/bin/sh
set -x
find cPhyProb -name "*.pyc" -exec rm {} \;
rm ./cPhyProb/ccore/dsct_model.so
rm -rf build
python setup.py build --use-beagle --debug --no-inline
cp build/lib.macosx-10.3-i386-2.6/cPhyProb/ccore/dsct_model.so cPhyProb/ccore

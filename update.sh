#!/usr/bin/env bash

set -e
set -x

SRC="tmp"
DEST="src/metis"

git clone https://github.com/KarypisLab/METIS.git $SRC
mkdir -p $DEST
mkdir -p $DEST/libmetis
mkdir -p $DEST/include

cp -f $SRC/LICENSE $DEST/LICENSE
cp -fa $SRC/include/*.h $DEST/include/
cp -fa $SRC/libmetis/*.{h,c} $DEST/libmetis/

# Remove files that are not needed
rm -rf $SRC

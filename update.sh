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
cp -fa $SRC/graphs ./

mkdir -p "app"
mkdir -p "app/cmpfillin"
mkdir -p "app/gpmetis"
mkdir -p "app/graphchk"
mkdir -p "app/m2gmetis"
mkdir -p "app/mpmetis"
mkdir -p "app/ndmetis"
mkdir -p "app/utils"

cp -fa $SRC/programs/*.{h,c} "app/utils"
mv app/utils/cmpfillin.c app/cmpfillin

mv app/utils/gpmetis.c app/gpmetis
mv app/utils/cmdline_gpmetis.c app/gpmetis

mv app/utils/graphchk.c app/graphchk

mv app/utils/m2gmetis.c app/m2gmetis
mv app/utils/cmdline_m2gmetis.c app/m2gmetis

mv app/utils/mpmetis.c app/mpmetis
mv app/utils/cmdline_mpmetis.c app/mpmetis

mv app/utils/ndmetis.c app/ndmetis
mv app/utils/cmdline_ndmetis.c app/ndmetis

# Hard link the remaining files in utils
ln -f app/utils/* app/cmpfillin
ln -f app/utils/* app/gpmetis
ln -f app/utils/* app/graphchk
ln -f app/utils/* app/m2gmetis
ln -f app/utils/* app/mpmetis
ln -f app/utils/* app/ndmetis

# Remove files that are not needed
rm -rf $SRC

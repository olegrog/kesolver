#!/bin/bash

tarfile=$1
#echo $tarfile

filename=$(basename "$tarfile")
#filename="${filename%%.*}"
filename="${filename/.tar.gz/}"
#echo $filename

tempdir=$(mktemp -d /tmp/${filename}.XXXXXXXXXX)
#echo $tempdir

tar -xz -f $tarfile -C $tempdir

mkdir src
cp $tempdir/$filename/src/lib_json/*.{h,inl,cpp} src

mkdir include 
cp -r $tempdir/$filename/include/* include

mkdir lib

rm -rf $tempdir

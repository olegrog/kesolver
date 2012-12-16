#!/bin/bash

msh=$1

kem=${msh%.*}.kem

if [ $# -gt 2 ]; then
    kep=$2
else
    kep=${msh%.*}.kep
fi

if [ $# -gt 3 ]; then
    kei=$3
else
    kei=${kep%.*}.kei
fi

dir=`dirname $0`

$dir/msh2kem.py $msh $kem
$dir/attach_mesh.py $kep $kem $kei


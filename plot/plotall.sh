#!/usr/bin/env bash

echo "[PLOTALL] ...plot principal axes of strain rates..."
./gmtstrainplot.sh -jpg -str strain_info.dat -psta -l -o output_str
echo "[PLOTALL] ...plot rotational rates..."
./gmtstrainplot.sh -jpg -rot strain_info.dat -psta -l -o output_rot
echo "[PLOTALL] ...plot dextral, sinistral maximum shear strain rates..."
./gmtstrainplot.sh -jpg -gtotaxes strain_info.dat -psta -l -o output_gtotaxes
echo "[PLOTALL] ...plot maximum shear strain..."
./gmtstrainplot.sh -jpg -gtot strain_info.dat -psta -l -o output_gtot
echo "[PLOTALL] ...plot dilatation..."
./gmtstrainplot.sh -jpg -dil strain_info.dat -psta -l -o output_dil
echo "[PLOTALL] ...plot second invariant..."
./gmtstrainplot.sh -jpg -secinv strain_info.dat -psta -l -o output_2inv 


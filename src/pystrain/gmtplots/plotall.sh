#!/usr/bin/env bash

echo "[PLOTALL] ...plot principal axes of strain rates..."
./gmtstrainplot.sh -jpg -str strain_info.dat -psta -l -o output_str -max_str_value 300
echo "[PLOTALL] ...plot rotational rates..."
./gmtstrainplot.sh -jpg -rot strain_info.dat -psta -l -o output_rot  -max_str_value 300
echo "[PLOTALL] ...plot dextral, sinistral maximum shear strain rates..."
./gmtstrainplot.sh -jpg -gtotaxes strain_info.dat -psta -l -o output_gtotaxes -max_str_value 300
echo "[PLOTALL] ...plot maximum shear strain..."
./gmtstrainplot.sh -jpg -gtot strain_info.dat -psta -l -o output_gtot -max_str_value 300
echo "[PLOTALL] ...plot dilatation..."
./gmtstrainplot.sh -jpg -dil strain_info.dat -psta -l -o output_dil   -max_str_value 300
echo "[PLOTALL] ...plot second invariant..."
./gmtstrainplot.sh -jpg -secinv strain_info.dat -psta -l -o output_2inv   -max_str_value 300


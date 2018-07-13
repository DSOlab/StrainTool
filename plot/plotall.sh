#!/usr/bin/env bash

echo "[PLOTALL] ...plot principal axes of strain rates..."
./gmtstrainplot.sh -jpg -vhor station_info.dat -psta -l -o output_vel
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

echo "[PLOTALL] ...plot statistics, stations used per cell..."
./gmtstatsplot.sh -jpg -stats strain_stats.dat --stats-stations -leg -o output_stats-stations
echo "[PLOTALL] ...plot statistics, optimal smoothing distance..."
./gmtstatsplot.sh -jpg -stats strain_stats.dat --stats-doptimal -leg -o output_stats-doptimal
echo "[PLOTALL] ...plot statistics, sigma value..."
./gmtstatsplot.sh -jpg -stats strain_stats.dat --stats-sigma -leg -o output_stats-sigma

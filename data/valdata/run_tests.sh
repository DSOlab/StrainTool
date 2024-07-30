#!/bin/bash

## Set bin directory
pth2st=${HOME}/Apps/STRAIN/StrainTool

echo ${pth2bin}
##Run test data


## Cases

## Parameters codes
# 1: [0: reference results 1:product]
# 2: File   [1:midas001 2:midas002, 3:midas003, 4:midas004]
# 3: step grid [0:0.5/0.5, 1:1/1, 2:0.25/0.25, 3:0.1/0.1 4:2/2]
# 4: -c [0:no 1:yes]
# 5: -b barycenter [0:False 1:set]
# 6-7: Wt [<=99]
# 8-9: dmin [< 100 km]
# 10-12: dmax [001 - 999]
# 13-14: dstep [01--99]
# 15: -g generate statistics [0:False, 1:true]
# 16: --multicore [0:False, 1:True]

# 10024015000211: --y-grid-step 0.5 --y-grid-step 0.5 --Wt 24 --dmin 1 --dmax 500 --dstep 2

## Default test for different files
for ivel in midas001.vel midas002.vel midas003.vel midas004.vel
do
    for gstep in 0.5 1
    do
        for Wt in 06 12 24 48
        do
            for dmax in 300 500 700
            do
                testcode=0${ivel:7:1}${gstep:0:1}00${Wt}01${dmax}0210
                ${pth2st}/bin/./StrainTensor.py -i ${pth2st}/data/valdata/${ivel} \
                    --x-grid-step ${gstep} \
                    --y-grid-step ${gstep} \
                    --Wt ${Wt} \
                    --dmax ${dmax} \
                    -g
                mv station_info.dat ${testcode}_station_info.dat
                mv strain_info.dat ${testcode}_strain_info.dat
                mv strain_stats.dat ${testcode}_strain_stats.dat
            done
        done
    done
done

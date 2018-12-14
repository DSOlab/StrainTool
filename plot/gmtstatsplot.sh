#!/usr/bin/env bash

# //////////////////////////////////////////////////////////////////////////////
# ==============================================================================
#
#    |===============================================|
#    |**       DIONYSOS SATELLITE OBSERVATORY      **|
#    |**          HIGHER GEODESY LABORATORY        **|
#    |**  National Technical University of Athens  **|
#    |===============================================|
#
#    filename       : gmtstatsnplot.sh
#                     NAME=gmtstrainplot
#    version        : v-1.0
#                     VERSION=v1.0
#                     RELEASE=rc4.0
#    licence        : MIT
#    created        : JUL-2018
#    usage          :
#    GMT Modules    : gmtset, makecpt, psbasemap, xyz2grd, grdsample, grdimage,
#                     pscoast, psscale, psxy, pstext, psvelo, psconvert, pscontour
#    UNIX progs     : awk 
#    exit code(s)   : 0 -> success
#                   : 1 -> error
#    discription    : 
#    uses           : 
#    notes          :
#    update list    : LAST_UPDATE=DEC-2018
#    contact        : Dimitris Anastasiou (dganastasiou@gmail.com)
#                     Xanthos Papanikolaou (xanthos@mail.ntua.gr)
#    ----------------------------------------------------------------------
# ==============================================================================

##
##  Function to resolve system's python version. The major version (aka 2 or 3)
##+ is stored in a variable called "PYV"
##
resolve_py_version() {
    regex="([1-9])\.[1-9]\.[1-9]+"
    pyv=$(python -c 'import platform; print(platform.python_version())')
    if [[ $pyv =~ $regex ]]
    then
        if test "${BASH_REMATCH[1]}" = 2
        then
            PYV=2
        elif test "${BASH_REMATCH[1]}" = 3
        then
            PYV=3
        else
            >&2 echo "[ERROR] Failed to resolve Python version"
            exit 1
        fi
    else
        >&2 echo "[ERROR] Failed to resolve Python version"
        exit 1
    fi
}
##
##  Alias python call! This is actualy an alias to calling: 'python -c'
##+ depending on the variable PYV; that is if PYV=3, then the call will be
##+ 'python -c ......', else the call will be
##+ 'python -c "from __future__ import print_function; .....'
## ----------------------------------------------------------------------------
##  So, when using pythonc function, just use the Python v3.x print syntax.
## ----------------------------------------------------------------------------
##
##  e.g
##  $>foo=$(pythonc "a=5+5.7; print(a)")
##
pythonc() {
    if test ${PYV} = 3
    then
        python -c "$@"
    else
        python -c "from __future__ import print_function; $@"
    fi
}

# //////////////////////////////////////////////////////////////////////////////
# HELP FUNCTION
function help {
	echo "/*****************************************************************/"
	echo " Program Name : gmtstatplot.sh"
	echo " Version : v-1.0"
	echo " Purpose : Plot statistics generated from straintool"
	echo " Usage   : gmtstatsplot.sh -stats [input] --stats-<option> -o [output]"
	echo " Switches: "
	echo ""
	echo "/*** Basic Plots & Background ***********************************/"
	echo "     -r | --region : region to plot (default Greece)"
	echo "         usage: -r west east south north projscale frame"
	echo ""
	echo "/*** PLOT STATIONS ***********************************************/"
	echo "     -psta [:=stations] plot only stations from input file"
	echo "     -deltr [:= delaunay triangles] plot delaunay triangles"
	echo ""
	echo "/*** PLOT STATISTICS ********************************************/"
	echo "     -stats (input file) set input file"
	echo "     --stats-stations : plot used stations"
	echo "     --stats-doptimal : plot optimal distance (D)"
	echo "     --stats-sigma : plot sigma "
	echo ""
	echo "/*** OTHER OPRTIONS ********************************************/"
	echo "     -o | --output : name of output files"
	echo "     -l | --labels : plot labels"
	echo "     -leg : plot legends"
	echo "     -mt | --map_title <title> : title map default none use quotes"
	echo "     -jpg : convert eps file to jpg"
	echo "     -h | --help : help menu"
	echo " Exit Status:    1 -> help message or error"
	echo " Exit Status:  = 0 -> sucesseful exit"
	echo " run scr: ./gmtstrainplot.sh -jpg -stats strain_stats.dat --stats-stations"
	echo "/*************************************************************/"
	exit 1
}
# //////////////////////////////////////////////////////////////////////////////
# BASH settings
# set -o errexit
set -e
set -o pipefail
# set -o nounset
# set -o xtrace

# //////////////////////////////////////////////////////////////////////////////
# pre define parameters

# program version
VERSION="v.1.0-rc4.0"

# system's Python version
PYV=99
resolve_py_version

##  verbosity level for GMT, see http://gmt.soest.hawaii.edu/doc/latest/gmt.html#v-full
##+ q - Complete silence. n - Normal. c - compatibility warnings. v - progress messages
##+ l - detailed progress messages. d - debugging messages.
export VRBLEVM=n

# //////////////////////////////////////////////////////////////////////////////
# Source function files

# //////////////////////////////////////////////////////////////////////////////
# Pre-defined parameters for bash script
TOPOGRAPHY=0
LABELS=0
OUTJPG=0
LEGEND=0
LOGO=0

PSTA=0
DELTR=0
STATS=0
STATS_STATIONS=0
STATS_DOPTIMAL=0
STATS_SIGMA=0

# //////////////////////////////////////////////////////////////////////////////
# Check default parameters file
if [ ! -f "default-param" ]
then
  echo "default-param file does not exist"
  exit 1
else
  echo "...load default parameters..."
  source default-param
fi

# //////////////////////////////////////////////////////////////////////////////
# GET COMMAND LINE ARGUMENTS
if [ "$#" == "0" ]
then
  help
fi

while [ $# -gt 0 ]
do
  case "$1" in
    -r | --region)
      west=$2
      east=$3
      south=$4
      north=$5
      projscale=$6
      frame=$7
      shift
      shift
      shift
      shift
      shift
      shift
      shift
      ;;
    -mt)
      maptitle=$2
      shift
      shift
      ;;
    -psta)
      pth2sta=${pth2inptf}/station_info.dat
      PSTA=1
      shift
      ;;
    -deltr)
      pth2deltr=${pth2inptf}/delaunay_info.dat
      DELTR=1
      shift
      ;;
    -stats)
      pth2stats=${pth2inptf}/$2
      STATS=1
      shift
      shift
      ;;
    --stats-stations)
      STATS_STATIONS=1
      maptitle="Stations used per grid cell"
      shift
      ;;
    --stats-doptimal)
      STATS_DOPTIMAL=1
      maptitle="Optimal Smoothing Distance (D) per grid cell"
      shift
      ;;
    --stats-sigma)
      STATS_SIGMA=1
      maptitle="sigma (@~\s@~) value estimated per grid cell"
      shift
      ;;
    -o | --output)
      outfile=${2}.ps
      shift
      shift
      ;;
    -l | --labels)
      LABELS=1
      shift
      ;;
    -leg)
      LEGEND=1
      shift
      ;;
    -logo)
      LOGO=1
      shift
      ;;
    -jpg)
      OUTJPG=1
      shift
      ;;
    -h | --help)
      help
      ;;
    -v | --version)
      echo "version: $VERSION"
      exit 1
      shift
      ;;
    *)
      echo "[ERROR] Bad argument structure. argument \"${1}\" is not right"
      echo "[STATUS] Script Finished Unsuccesfully! Exit Status 1"
      exit 1
  esac
done

# //////////////////////////////////////////////////////////////////////////////
# check if files exist

##check statistics input file
if [ "$STATS" -eq 1 ]
then
  if [ ! -f $pth2stats ]
  then
    echo "[ERROR] input file $pth2stats does not exist"
    echo "          please set the correct path and then use this switch"
    echo "[STATUS] Script Finished Unsuccesfully! Exit Status 1"
    exit 1
  fi
fi

##check delaunay triangles file
if [ "$DELTR" -eq 1 ]
then
  if [ ! -f $pth2deltr ]
  then
    echo "[ERROR] input file $pth2deltr does not exist"
    echo "          please set the correct path and then use this switch"
    echo "[STATUS] Script Finished Unsuccesfully! Exit Status 1"
    exit 1
  fi
fi

##check stations info file
if [ "$PSTA" -eq 1 ]
then
  if [ ! -f $pth2sta ]
  then
    echo "[WARNING] input file $pth2sta does not exist"
    echo "          Stations will not printed"
    PSTA=0
  fi
fi

# //////////////////////////////////////////////////////////////////////////////
# READ STATISTICS FILE
if [ "$STATS" -eq 1 ]
then
  # calculate grid size
  west_grd=$(cat $pth2stats | awk "/^Longtitude/,0" | tail -n +3 \
  | awk ' {printf "%+0.3f\n", $1}' | gmt info -El)
  east_grd=$(cat $pth2stats | awk "/^Longtitude/,0" | tail -n +3 \
  | awk ' {printf "%+0.3f\n", $1}' | gmt info -Eh)
  south_grd=$(cat $pth2stats | awk "/^Longtitude/,0" | tail -n +3 \
  | awk ' {printf "%+0.3f\n", $2}' | gmt info -El)
  north_grd=$(cat $pth2stats | awk "/^Longtitude/,0" | tail -n +3 \
  | awk ' {printf "%+0.3f\n", $2}' | gmt info -Eh)
  range_grd="-R${west_grd}/${east_grd}/${south_grd}/${north_grd}"
    
  > .legend
  {
  globFonts="G .5c"
  echo "G 0.2c"
  echo "H 11 Times-Roman StrainTool parameters"
  echo "D 0.3c 1p"
  echo "N 1"

  stat_version=$(grep Version $pth2stats  |awk -F: '{print $2}')
  echo -e "T Version: ${stat_version}\n${globFonts}"
  stat_gps_file=$(tail -n+4 $pth2stats | grep gps_file | awk '{print $3}')
  echo -e "T GPS file: ${stat_gps_file}\n${globFonts}"
  echo "H 11 Times-Roman Interpolation model"
  echo "D 0.3c 1p"
  echo "N 1"
  stat_method=$(tail -n+4 $pth2stats | grep method | awk '{print $3}')
  echo -e "T method: ${stat_method}\n${globFonts}"
  stat_ltype=$(tail -n+4 $pth2stats | grep ltype | awk '{print $3}')
  echo -e "T ltype: ${stat_ltype}\n${globFonts}"
  stat_Wt=$(tail -n+4 $pth2stats | grep Wt | awk '{print $3}')
  echo -e "T Wt: ${stat_Wt}\n${globFonts}"
  stat_dmin=$(tail -n+4 $pth2stats | grep dmin | awk '{print $3}')
  echo -e "T dmin: ${stat_dmin}\n${globFonts}"
  stat_dmax=$(tail -n+4 $pth2stats | grep dmax | awk '{print $3}')
  echo -e "T dmax: ${stat_dmax}\n${globFonts}"
  stat_dstep=$(tail -n+4 $pth2stats | grep dstep | awk '{print $3}')
  echo -e "T dstep: ${stat_dstep}\n${globFonts}"
  echo "H 11 Times-Roman Region parameters"
  echo "D 0.3c 1p"
  echo "N 1"
  echo -e "T @_Grid limits@_:\n${globFonts}"
  echo -e "T west / east: ${west_grd} / ${east_grd}\n${globFonts}"
  echo -e "T south / north: ${south_grd} / ${north_grd}\n${globFonts}"
  stat_x_grid_step=$(tail -n+4 $pth2stats | grep x_grid_step | awk '{print $3}')
  echo -e "T x_grid_step: ${stat_x_grid_step}\n${globFonts}"
  stat_y_grid_step=$(tail -n+4 $pth2stats | grep y_grid_step | awk '{print $3}')
  echo -e "T y_grid_step: ${stat_y_grid_step}\n${globFonts}"
  } >> .legend
  istep_grd=$(pythonc "print(${stat_x_grid_step}*60.)")
fi

# //////////////////////////////////////////////////////////////////////////////
# SET REGION PROPERTIES
tmp_scrate=$(pythonc "print((${projscale}/150000000.)*10.)")
sclat=$(pythonc "print(${south} + ${tmp_scrate})")

tmp_scrate=$(pythonc "print((${projscale}/150000000.)*27.)")
sclon=$(pythonc "print(${east} - ${tmp_scrate})")

tmp_msclat=$(pythonc "print(int((${south} + ${north})/2))")
tmp_msclon=$(pythonc "print(int((${west} + ${east})/2))")
export scale=-Lf${sclon}/${sclat}/${tmp_msclat}:${tmp_msclon}/${sclength}+l+jr
# scale="-Lf20/33.5/36:24/100+l+jr"
range="-R$west/$east/$south/$north"
proj="-Jm24/37/1:$projscale"

# //////////////////////////////////////////////////////////////////////////////
# GMT parameters
gmt gmtset MAP_FRAME_TYPE fancy
gmt gmtset PS_PAGE_ORIENTATION portrait
gmt gmtset FONT_ANNOT_PRIMARY 8 FONT_LABEL 8 MAP_FRAME_WIDTH 0.12c FONT_TITLE 18p,Palatino-BoldItalic
gmt gmtset PS_MEDIA ${PAPER_SIZE}

# ####################### TOPOGRAPHY ###########################
if [ "$TOPOGRAPHY" -eq 0 ]
then
  echo "...plot coastlines..."
  ################## Plot coastlines only ######################	
  gmt	psbasemap $range $proj  -B$frame:."$maptitle": -P -K > $outfile
  gmt	pscoast -R -J -O -K -W0.25 -G225 -Df -Na $scale >> $outfile
fi

# //////////////////////////////////////////////////////////////////////////////
### PLOT ONLY STATIONS ITHOUT ANY OTHER PARAMETER

if [ "$PSTA" -eq 1 ] && [ "$STATS" -eq 0 ] 
then

    awk 'NR > 2 {print $2,$3}' $pth2sta  \
    | gmt psxy -R -J -W.1 -Sc.15c -Gyellow -O -K -V${VRBLEVM} >> $outfile
    
    if [ "$LABELS" -eq 1 ]
    then
      awk 'NR > 2 {print $2,$3, "7,1,black", 0, "RB", $1}' $pth2sta \
      | gmt pstext -R -J -Dj0.1c/0.1c -F+f+a+j -O -K -V${VRBLEVM} >> ${outfile}
    fi
fi

# //////////////////////////////////////////////////////////////////////////////
### PLOT STATIONS USED FOR EACH CELL
if [ "$STATS_STATIONS" -eq 1 ]
then
  echo "...plot stations used for each grid cell..."
  awk 'NR > 24 {print $1,$2,$3}' $pth2stats > tmpstations
  # find min max and create cpt file
  T=`awk '{print $3}' tmpstations | gmt info -Eh `
  # set variables for scale
  if [ $(echo " ${T} <= 10 " | bc -l) == 1 ]
  then
    Tmax_r=0
    Tmax_r_marg=1
    cpt_step=1
    scale_step_r=0
  elif [ $(echo " ${T} > 10 " | bc -l) == 1 ] && [ $(echo " ${T} <= 100 " | bc -l) == 1 ]
  then
    Tmax_r=0
    Tmax_r_marg=5
    cpt_step=1
    scale_step_r=0
  elif [ $(echo " ${T} > 100 " | bc -l) == 1 ] && [ $(echo " ${T} <= 25000 " | bc -l) == 1 ]
  then
    Tmax_r=-1
    Tmax_r_marg=10
    cpt_step=1
    scale_step_r=-1
  else
    echo "ERROR"
    exit 1
  fi
  Tmax=$(pythonc "print(int(round(${T},${Tmax_r})+${Tmax_r_marg}))")
  T=`awk '{print $3}' tmpstations | gmt info -El `
  Tmin=$(pythonc "print(int(round(${T},${Tmax_r})-${Tmax_r_marg}))")
  gmt makecpt -Cjet -T${Tmin}/${Tmax}/${cpt_step} > inx.cpt  
  gmt xyz2grd tmpstations -Gtmpstations.grd ${range_grd} -I${istep_grd}m -V${VRBLEVM}
  gmt grdsample tmpstations.grd -I${istep_grd}m -Gtmpstations_sample.grd -V${VRBLEVM}
  gmt grdimage tmpstations_sample.grd ${proj} ${range} -Cinx.cpt -Q \
	-O -K -V${VRBLEVM}>> $outfile
  scale_step=$(pythonc "print(round(((${Tmax}-${Tmin})/5.),${scale_step_r}))")
  gmt pscoast -R -J -O -K -W0.25 -Df -Na -V${VRBLEVM}>> $outfile
  gmt psscale -Cinx.cpt -D8/-1.1/10/0.3h -B${scale_step}/:"stations": -I -S \
      -O -K -V${VRBLEVM}>> $outfile
  
# plot stations
  if [ "$PSTA" -eq 1 ]
  then

    awk 'NR > 2 {print $2,$3}' $pth2sta  \
    | gmt psxy -R -J -W.1 -Sc.15c -Gyellow -O -K -V${VRBLEVM} >> $outfile
    
    if [ "$LABELS" -eq 1 ]
    then
      awk 'NR > 2 {print $2,$3, "7,1,black", 0, "RB", $1}' $pth2sta \
      | gmt pstext -R -J -Dj0.1c/0.1c -F+f+a+j -O -K -V${VRBLEVM} >> ${outfile}
    fi
  fi
fi

# //////////////////////////////////////////////////////////////////////////////
### PLOT OPTIMAL DISTANCE d
if [ "$STATS_DOPTIMAL" -eq 1 ]
then
  echo "...plot optimal distance D for each grid cell..."
  awk 'NR > 24 {print $1,$2,$4}' $pth2stats > tmpdoptimal
# find min max and create cpt file
  T=`awk '{print $3}' tmpdoptimal | gmt info -Eh `
  if [ $(echo " ${T} <= 10 " | bc -l) == 1 ]
  then
    Tmax_r=0
    Tmax_r_marg=1
    cpt_step=1
    scale_step_r=0
  elif [ $(echo " ${T} > 10 " | bc -l) == 1 ] && [ $(echo " ${T} <= 100 " | bc -l) == 1 ]
  then
    Tmax_r=0
    Tmax_r_marg=5
    cpt_step=1
    scale_step_r=0
  elif [ $(echo " ${T} > 100 " | bc -l) == 1 ] && [ $(echo " ${T} <= 25000 " | bc -l) == 1 ]
  then
    Tmax_r=-1
    Tmax_r_marg=10
    cpt_step=1
    scale_step_r=-1
  else
    echo "ERROR"
    exit 1
  fi
  Tmax=$(pythonc "print(round(${T},${Tmax_r})+${Tmax_r_marg})")
  T=`awk '{print $3}' tmpdoptimal | gmt info -El `
  Tmin=$(pythonc "print(round(${T},${Tmax_r})-${Tmax_r_marg})")
  gmt makecpt -Cjet -T${Tmin}/${Tmax}/${cpt_step} > inx.cpt  
  gmt xyz2grd tmpdoptimal -Gtmpdoptimal.grd ${range_grd} -I${istep_grd}m -V${VRBLEVM}
  gmt grdsample tmpdoptimal.grd -I${istep_grd}m -Gtmpdoptimal_sample.grd -V${VRBLEVM}
  gmt grdimage tmpdoptimal_sample.grd ${proj} ${range} -Cinx.cpt -Q \
	-O -K -V${VRBLEVM}>> $outfile

  scale_step=$(pythonc "print(round(((${Tmax}-${Tmin})/5.),${scale_step_r}))")
  gmt pscoast -R -J -O -K -W0.25 -Df -Na -V${VRBLEVM}>> $outfile
  gmt psscale -Cinx.cpt -D8/-1.1/10/0.3h -B${scale_step}/:"km": -I -S \
      -O -K -V${VRBLEVM}>> $outfile
  
# plot stations
  if [ "$PSTA" -eq 1 ]
  then

    awk 'NR > 2 {print $2,$3}' $pth2sta  \
    | gmt psxy -R -J -W.1 -Sc.15c -Gyellow -O -K -V${VRBLEVM} >> $outfile
    
    if [ "$LABELS" -eq 1 ]
    then
      awk 'NR > 2 {print $2,$3, "7,1,black", 0, "RB", $1}' $pth2sta \
      | gmt pstext -R -J -Dj0.1c/0.1c -F+f+a+j -O -K -V${VRBLEVM} >> ${outfile}
    fi
  fi
fi

# //////////////////////////////////////////////////////////////////////////////
### PLOT SIGMA
if [ "$STATS_SIGMA" -eq 1 ]
then
  echo "...plot sigm estimated for each grid cell..."
  awk 'NR > 24 {print $1,$2,$6}' $pth2stats > tmpsigma
  # find min max and create cpt file
  T=`awk '{print $3}' tmpsigma | gmt info -Eh `
  Tmax=$(pythonc "print(round(${T},3)+.001)")
## comm. I am not sure that this check needed for sigma value
#   if [ $(echo " ${T} <= 10 " | bc -l) == 1 ]
#   then
#     Tmax_r=0
#     Tmax_r_marg=1
#     cpt_step=1
#     scale_step_r=0
#   elif [ $(echo " ${T} > 10 " | bc -l) == 1 ] && [ $(echo " ${T} <= 100 " | bc -l) == 1 ]
#   then
#     Tmax_r=0
#     Tmax_r_marg=5
#     cpt_step=1
#     scale_step_r=0
#   elif [ $(echo " ${T} > 100 " | bc -l) == 1 ] && [ $(echo " ${T} <= 25000 " | bc -l) == 1 ]
#   then
#     Tmax_r=-1
#     Tmax_r_marg=10
#     cpt_step=1
#     scale_step_r=-1
#   else
#     echo "ERROR"
#     exit 1
#   fi
  T=`awk '{print $3}' tmpsigma | gmt info -El `
  Tmin=$(pythonc "print(round(${T},3)-.001)")
  gmt makecpt -Cjet -T${Tmin}/${Tmax}/.001 > inx.cpt  
  
  gmt xyz2grd tmpsigma -Gtmpsigma.grd ${range_grd} -I${istep_grd}m -V${VRBLEVM}
  gmt grdsample tmpsigma.grd -I${istep_grd}m -Gtmpsigma_sample.grd -V${VRBLEVM}
  gmt grdimage tmpsigma_sample.grd ${proj} ${range} -Cinx.cpt -Q \
	-O -K -V${VRBLEVM}>> $outfile

  scale_step=$(pythonc "print(round((${Tmax}/5.),3))")
  gmt pscoast -R -J -O -K -W0.25 -Df -Na -V${VRBLEVM}>> $outfile
  gmt psscale -Cinx.cpt -D8/-1.1/10/0.3h -B${scale_step}/:"sigma": -I -S \
      -O -K -V${VRBLEVM}>> $outfile
  
# plot stations
  if [ "$PSTA" -eq 1 ]
  then

    awk 'NR > 2 {print $2,$3}' $pth2sta  \
    | gmt psxy -R -J -W.1 -Sc.15c -Gyellow -O -K -V${VRBLEVM} >> $outfile
    
    if [ "$LABELS" -eq 1 ]
    then
      awk 'NR > 2 {print $2,$3, "7,1,black", 0, "RB", $1}' $pth2sta \
      | gmt pstext -R -J -Dj0.1c/0.1c -F+f+a+j -O -K -V${VRBLEVM} >> ${outfile}
    fi
  fi
fi

# //////////////////////////////////////////////////////////////////////////////
# plot legend
if [ "$LEGEND" -eq 1 ]
then
  gmt pslegend .legend ${legendc} -C0.1c/0.1c  -O -K -V${VRBLEVM} >> $outfile
fi

# //////////////////////////////////////////////////////////////////////////////
#  Plot custom logo
if [ "$LOGO" -eq 1 ]
then
  gmt psimage $pth2logos -O $logo_pos2 -W1.1c -F0.4  -K -V${VRBLEVM} >>$outfile
fi

# //////////////////////////////////////////////////////////////////////////////
# FINAL SECTION
#################--- Close ps output file ----##################################
echo "$west $south 8,0,black 0 LB This image was produced using" \
  | gmt pstext -Jm ${range} -Dj0.1c/1.1c -F+f+a+j -K  -O -V${VRBLEVM} >> $outfile
echo "$west $south 9,1,white 0 LB STRAINTOOL for EPOS" \
  | gmt pstext -Jm -R -Dj0.2c/.65c -F+f+a+j -G165/0/236 -U$logo_pos -O -V${VRBLEVM} >> $outfile


#################--- Convert to other format ----###############################
if [ "$OUTJPG" -eq 1 ]
then
  echo "...adjust and convert to JPEG format..."
  gmt psconvert ${outfile} -A0.2c -Tj -V${VRBLEVM} 
fi

# clear all teporary files
echo "...remove temporary files..."
rm -rf tmp* gmt.conf gmt.history .legend inx.cpt 2>/dev/null

# Print exit status
echo "[STATUS] Finished. Exit status: $?"

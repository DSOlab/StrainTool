#!/usr/bin/env bash

# program version
VERSION="gmtstrainplot - v1.0"

# //////////////////////////////////////////////////////////////////////////////
# ==============================================================================
#
#    |===============================================|
#    |**       DIONYSOS SATELLITE OBSERVATORY      **|
#    |**          HIGHER GEODESY LABORATORY        **|
#    |**  National Technical University of Athens  **|
#    |===============================================|
#
#    filename       : gmtstrainplot.sh
#                     NAME=gmtstrainplot
#    version        : 1.0
#                     VERSION=1.0
#                     RELEASE=1.0
#    licence        : MIT
#    created        : MAY-2018
#    usage          :
#    GMT Modules    : gmtset, makecpt, psbasemap, xyz2grd, grdsample, grdimage,
#                     pscoast, psscale, psxy, pstext, psvelo, psconvert, pscontour
#    UNIX progs     : awk, tail
#    exit code(s)   : 0 -> success
#                   : 1 -> error
#    discription    : 
#    uses           : 
#    notes          :
#    update list    : LAST_UPDATE=JUN-2019
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

##
##  Function to set scale variables. User must pass in the T variable, and the
##+ function will set the following (global) variables:
##+     * T
##+     * Tmax_r
##+     * Tmax_r_marg
##+     * cpt_step
##+     * scale_step_r
##  Use as: <scalevar_T Tval> where Tvar must be a number (integer or float)
##
scalevar_T() 
{
    ## make sure function i called with a cmd, which is a number
    re="^[+-]?[0-9]+([.][0-9]+)?$"
    if test -z ${1+x} 
    then
        echo "[ERROR] Must supply cmd arg in scalevar_T" && exit 1
    else
        if ! [[ $1 =~ $re ]]
        then
            echo "[ERROR] Must supply numeric cmd arg in scalevar_T" && exit 1
        fi
    fi
    T="${1}"
    ## assign (scale) results to array (via awk)
    local SVARS=($(awk -v T="$T" '
        BEGIN {
            if (T<=1) {
                Tmax_r=1
                Tmax_r_marg=0.1
                cpt_step=0.01
                scale_step_r=1
            } else if (T>1 && T<10) {
                Tmax_r=0
                Tmax_r_marg=1
                cpt_step=0.01
                scale_step_r=0
            } else if (T>=10 && T<100) {
                Tmax_r=-1
                Tmax_r_marg=10
                cpt_step=0.1
                scale_step_r=-1
            } else if (T>=100) {
                Tmax_r=-1
                Tmax_r_marg=10
                cpt_step=1
                scale_step_r=-1
            } else {
                print "[ERROR] Cannot resolve scale variables"
                exit 1
            }
            print Tmax_r, Tmax_r_marg, cpt_step, scale_step_r;
        }
    '))
    ## check and assign to (global) variables
    test "${#SVARS[@]}" -eq 4  \
        || { echo "[ERROR] Failed to resolve scale variables"; exit 1; }
    read Tmax_r Tmax_r_marg cpt_step scale_step_r <<< "${SVARS[@]}"
}

##
##  Check if a (given) variable is a number (either integer or float, in
##+ fixed or scietific format. Returns 0 if variable is number, 1 otherwise.
##
##  Examples:
##      for s in abc 2.a 2 -2 +2.234 2.0e0 -123.45e-9 .25 . " " "" ; do
##          isNumber $s && echo "$s is number"
##      done
##  Result:
##      2 is number
##      -2 is number
##      +2.234 is number
##      2.0e0 is number
##      -123.45e-9 is number
##     .25 is number
##
isNumber() {
    test -z $1 && return 1
    numre="^[+-]?([0-9]+)?([.][0-9]+)?([eE][+-]?[0-9]{1,2})?$"
    [[ "$1" =~ $numre ]] && return 0;
    return 1
}

##
##  This function will read in a variable. If (and only if) the variable is
##+ in the form "a/b/c/d/e/f" (aka 6 numbers, seperated by '/'), then the
##+ function will set the global variables: 
##+ west, east, south, north, projscale, frame to the fields of the string;
##+ that is west=a, east=b, ..., frame=f
##  If the resolution happens successefuly and the variables are assigned,
##+ 0 is returned; in any other case, the variables are not set and the
##+ function returns 1.
##  Note that the string must comprise of numbers, i.e. in the passed in
##+ argument "a/b/c/d/e/f" a, b, c, .., f must be numbers.
##
##  Example:
##      resolve_region "1/2/3//5/6"
##      echo "vars: $west $east $south $north $projscale $frame"
##      resolve_region "1/2/3"
##      echo "vars: $west $east $south $north $projscale $frame"
##      resolve_region "foobar"
##      echo "vars: $west $east $south $north $projscale $frame"
##      resolve_region "1/2/3/4/foo/6"
##      echo "vars: $west $east $south $north $projscale $frame"
##      resolve_region "1/2e0/3/4/5.0/6"
##      echo "vars: $west $east $south $north $projscale $frame"
##  Prints:
##      vars:
##      vars:
##      vars:
##      vars:
##      vars: 1 2e0 3 4 5.0 6
##
resolve_region() {
#    local ar=$(echo $1 | awk -F"/" '{ 
#        if (NF == 6)
#            print $0 
#        else 
#            print "ERROR"
#    }')
#    test "$ar" == "ERROR" && return 1
    local AR=( $(echo $1 | tr '/' ' ') )
    test "${#AR[@]}" -eq 6 || return 1
    for num in "${AR[@]}" ; do isNumber "$num" || return 1 ; done
    read west east south north projscale frame <<< "${AR[@]}"
}

##
##  Check the inputs for the boundary zone that all inuts are ok. That's meean:
##  west < east and (east - west) < 360
##  -90 <= south < north <= 90
##  projscale > 0 and frame > 0 
##
check_region() {
  ## test south north -90 < south < north < 90
  local SVARS=($(awk -v west="$1" -v east="$2" -v south="$3" -v north="$4" -v projscale="$5" -v frame="$6" '
        BEGIN {
            if (-90 <= south && south < north && north <= 90) {
                  south = south
                  north = north
                if (west < east && (east - west) < 360 ) {
                  west = west
                  east = east
				  if (projscale > 0 && frame > 0) {
				    projscale = projscale
				    frame = frame
				  } else {
				    print "[ERROR] check again projscale and frame inputs"
                    exit 1
				  }
                } else {
                  print "[ERROR] check again west east inputs"
                  exit 1
                }
			} else {
                print "[ERROR] check again south north inputs"
                exit 1
            }
            print west, east, south, north, projscale, frame;
        }
    '))
    ## check and assign to (global) variables
    test "${#SVARS[@]}" -eq 6  \
        || { echo "[ERROR] Failed to modify boundary zone"; exit 1; }
    read west east south north projscale frame <<< "${SVARS[@]}"
}

# //////////////////////////////////////////////////////////////////////////////
# HELP FUNCTION
function help {
	echo "/*****************************************************************/"
	echo " Program Name : gmtstrainplot.sh"
	echo " Version : ${VERSION}"
	echo " Purpose : Plot strain maps for StrainTool results"
	echo " Usage   : gmtstrainplot.sh -r  |  | -o [output] | -jpg "
	echo " Switches: "
	echo ""
	echo "/*** Basic Plots & Background ***********************************/"
	echo "     -r | --region : region to plot (default Greece)"
	echo "         usage: -r west east south north projscale frame"
        echo "            or: -r west/east/south/north/projscale/frame"
	echo ""
	echo "/*** PLOT STATIONS ***********************************************/"
	echo "     -psta [:=stations] plot only stations from input file"
	echo "     -deltr [:= delaunay triangles] plot delaunay triangles"
	echo ""
	echo "/*** PLOT VELOCITIES ********************************************/"
	echo "     -vhor (station_file)[:= horizontal velocities]  "
	echo "     -vsc [:=velocity scale] change velocity scale default 0.05"
	echo ""
	echo "/*** PLOT STRAINS **********************************************/"
	echo "     -str (strain file)[:= strains] Plot strain rates "
	echo "     -rot (strain file)[:= rots] Plot rotational rates "
	echo "     -gtot(strain file)[:=shear strain] plot total shear strain rate contours"
	echo "     -gtotaxes (strain file) dextral and sinistral max shear strain rates"
	echo "     -dil (strainfile)[:= dilatation] Plot dilatation and principal axes"
	echo "     -secinv (strain file) [:=2nd invariand] Plot second invariand"
	echo "     -strsc [:=strain scale]"
	echo "     -rotsc [:=rotational scales]"
	echo "  *for -gtot | -dil | -secinv use +grd to plot gridded data"
	echo "        ex:-gtot+grd "
	echo ""
	echo "/*** OTHER OPRTIONS ********************************************/"
	echo "     -o | --output : name of output files"
	echo "     -l | --labels : plot labels"
	echo "     -mt | --map_title <title> : title map default none use quotes"
	echo "     -jpg : convert eps file to jpg"
	echo "     -h | --help : help menu"
	echo " Exit Status:    1 -> help message or error"
	echo " Exit Status:  = 0 -> sucesseful exit"
	echo " run scr: ./gmtstrainplot.sh -jpg -str strain_info.dat -psta -l"
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
FAULTS=0
PCMT=0
TOPOGRAPHY=0
LABELS=0
OUTJPG=0
LEGEND=0
LOGO=0

PSTA=0
DELTR=0
VHORIZONTAL=0
STRAIN=0
STRROT=0
GTOTAL=0
GTOTALAXES=0
DILATATION=0
SECINV=0
GRDDAT=0

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
        ## check the next argument; maybe its w/e/s/n/p/f
        if resolve_region "${2}"
        then
            shift 2
        else
            ## at least 6 (more) arguments should follow, seperated by space(s)
            test "$#" -ge 7 || { echo "[ERROR] Invalid region arg"; exit 1; }
            read west east south north projscale frame <<< $(echo $2 $3 $4 $5 $6 $7)
            shift 7
            ## all arguments should be numeric
            for rarg in "$west" "$east" "$south" "$north" "$projscale" "$frame"
            do
                isNumber "${rarg}" || \
                    { echo "[ERROR] Invalid region arg"; exit 1; }
            done
        fi
        ;;
    -mt)
      maptitle=$2
      shift 2
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
    -vhor)
      pth2stainfo=${pth2inptf}/$2
      maptitle="Velocity field"
      VHORIZONTAL=1
      shift 2
      ;;
    -vsc)
      VSC=$2
      shift 2
      ;;
    -str)
      pth2strinfo=${pth2inptf}/${2}
      maptitle="Principal Axes of Strain Rates"
      STRAIN=1
      shift 2
      ;;
    -strsc)
      STRSC=$2
      shift 2
      ;;
    -rot)
      pth2strinfo=${pth2inptf}/${2}
      maptitle="Rotational Rates"
      STRROT=1
      shift 2
      ;;
    -rotsc)
      ROTSC=${2}
      shift 2
      ;;
    -gtot*)
      if [ "${1:5:8}" == "+grd" ]; then
        pth2strinfo=${pth2inptf}/${2}
        GTOTAL=1
        GRDDAT=1
      elif [ "${1:5:8}" == "axes" ]; then
        pth2strinfo=${pth2inptf}/${2}
        GTOTALAXES=1
        maptitle="Axes of dextral/sinistral shear-strain"
      else
        pth2strinfo=${pth2inptf}/${2}
        GTOTAL=1
        maptitle="Maximum Shear Strain Rates"
      fi
      shift 2
      ;;
    -gtotaxes)
      pth2strinfo=${pth2inptf}/${2}
      GTOTALAXES=1
      maptitle="Axes of dextral/sinistral shear strain"
      shift 2
      ;;
    -dil*)
      pth2strinfo=${pth2inptf}/${2}
      DILATATION=1
      maptitle="Dilatation"
      if [ "${1:4:7}" == "+grd" ]; then
        GRDDAT=1
      fi
      shift 2
      ;;
    -secinv*)
      pth2strinfo=${pth2inptf}/${2}
      SECINV=1
      maptitle="Second invariant of strain rate"
      if [ "${1:7:10}" == "+grd" ]; then
        GRDDAT=1
      fi
      shift 2
      ;;
    -o | --output)
      outfile=${2}.ps
      shift 2
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
      echo "version: ${VERSION}"
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

# Check strain info file exist
if [ "$STRAIN" -eq 1 ] || [ "$STRROT" -eq 1 ] || [ "$GTOTAL" -eq 1 ] || [ "$DILATATION" -eq 1 ] \
   || [ "$GTOTALAXES" -eq 1 ] || [ "$SECINV" -eq 1 ]
then
  if [ ! -f $pth2strinfo ]
  then
    echo "[ERROR] input file $pth2strinfo does not exist"
    echo "          please download it and then use this switch"
    echo "[STATUS] Script Finished Unsuccesfully! Exit Status 1"
    exit 1
  fi
fi

# check delaunay info file
if [ "$DELTR" -eq 1 ]
then
  if [ ! -f $pth2deltr ]
  then
    echo "[ERROR] input file $pth2deltr does not exist"
    echo "          please download it and then use this switch"
    echo "[STATUS] Script Finished Unsuccesfully! Exit Status 1"
    exit 1
  fi
fi

# check station input files
if [ "$PSTA" -eq 1 ]
then
  if [ ! -f $pth2sta ]
  then
    echo "[WARNING] input file $pth2sta does not exist"
    echo "          Stations will not printed"
    PSTA=0
  fi
fi

# Check velocity file
if [ "$VHORIZONTAL" -eq 1 ]
then
  if [ ! -f $pth2stainfo ]
  then
    echo "[WARNING] input file $pth2stainfo does not exist"
    echo "          Horizontal velocities will not plotted"
    VHORIZONTAL=0
  fi
fi

# //////////////////////////////////////////////////////////////////////////////
# SET REGION PROPERTIES
# check region inputs
check_region $west $east $south $north $projscale $frame

#set region parameters
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
if [ "$TOPOGRAPHY" -eq 0 ]  # For v1.0 var $TOPOGRAPHY set only 0
then
  echo "...plot coastlines..."
  ################## Plot coastlines only ######################	
  gmt	psbasemap $range $proj  -B$frame:."$maptitle": -P -K > $outfile
  gmt	pscoast -R -J -O -K -W0.25 -G225 -Df -Na $scale >> $outfile
fi

# //////////////////////////////////////////////////////////////////////////////
### PLOT ONLY STATIONS WITHOUT ANY OTHER PARAMETER

if [ "$PSTA" -eq 1 ] && [ "$STRAIN" -eq 0 ] && [ "$STRROT" -eq 0 ] && [ "$GTOTAL" -eq 0 ] \
    && [ "$DILATATION" -eq 0 ] && [ "$GTOTALAXES" -eq 0 ] && [ "$SECINV" -eq 0 ]
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
### PLOT HORIZONTAL VELOCITIES

if [ "$VHORIZONTAL" -eq 1 ]
then
  echo "...plot horizontal velocities..."
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

  awk 'NR > 2 {print $2,$3,$4,$5,$6,$7,0,$1}' $pth2stainfo \
  | gmt psvelo -R -Jm -Se${VSC}/0.95/0 -W.5p,black -A.05p+e -Gblue \
    -O -K -V${VRBLEVM} >> $outfile

  awk 'NR > 2 {print $2,$3,$4,$5,$6,$7,0,$1}' $pth2stainfo \
  | gmt psvelo -R -Jm -Se${VSC}/0/0 -W2p,blue -A10p+e -Gblue \
    -O -K -V${VRBLEVM} >> $outfile

# plot scale of strain rates
  tmp_scrate=$(pythonc "print((${projscale}/150000000.)*5.)")
  velsclat=$(pythonc "print(${sclat} + ${tmp_scrate})")
  velsclon=$(pythonc "print(${sclon} - ${tmp_scrate})")
  vscmagn_sd=$(pythonc "print(${vscmagn}/20.e0)")

  echo "$velsclon $velsclat ${vscmagn} 0 ${vscmagn_sd} ${vscmagn_sd}  0 " \
  | gmt psvelo -R -Jm -Se${VSC}/0.95/0 -W.5p,black -A.05p+e -Gblue \
    -O -K -V${VRBLEVM} >> $outfile
  
  echo "$velsclon $velsclat ${vscmagn} 0 0 0 0 " \
  | gmt psvelo -R -Jm -Se${VSC}/0/0 -W2p,blue -A10p+e -Gblue \
    -O -K -V${VRBLEVM} >> $outfile
  
  echo "$sclon $velsclat 9,1,black 0 CB ${vscmagn} \261 ${vscmagn_sd} mm/y" \
  | gmt pstext -Jm -R -Dj0c/.5c -F+f+a+j  -O -K -V${VRBLEVM} >> $outfile
  
  sclat=${velsclat}
#   sclon=${velsclon}
fi

# //////////////////////////////////////////////////////////////////////////////
### PLOT SHEAR STRAIN RATES parameters
if [ "$GTOTAL" -eq 1 ]
then
  echo "...plot maximum shear strain rates..."
# plot shear strain rates
  awk 'NR > 2 {print $2,$1,$19*1}' $pth2strinfo > tmpgtot
# find min max and create cpt file
  T=`awk '{print $3}' tmpgtot | gmt info -Eh `
# set variables for scale plot
  scalevar_T ${T}
  Tmax=$(pythonc "print(round(${T},${Tmax_r})+${Tmax_r_marg})")
  
  gmt makecpt -Cjet -T0/${Tmax}/${cpt_step} > inx.cpt
  
  if [ "${GRDDAT}" -eq 0 ]
  then
    gmt pscontour tmpgtot -R -J  -Cinx.cpt -I0.1 -O -K -V${VRBLEVM} >> ${outfile}
  else
    gmt xyz2grd tmpgtot -Gtmpgtot.grd ${range} -I40m= -V
    gmt grdsample tmpgtot.grd -I4s -Gtmpgtot_sample.grd -V${VRBLEVM}
    gmt grdimage tmpgtot_sample.grd ${proj} ${range} -Cinx.cpt -Q \
	-O -K -V${VRBLEVM}>> $outfile
  fi

# PLot scale
  scale_step=$(pythonc "print(round((${Tmax}/5.),${scale_step_r}))")
  gmt pscoast -R -J -O -K -W0.25 -Df -Na -V${VRBLEVM}>> $outfile
  gmt psscale -Cinx.cpt -D8/-1.1/10/0.3h -B${scale_step}/:"nstrain/y": -I -S \
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
### PLOT SHEAR STRAIN RATES parameters
if [ "$DILATATION" -eq 1 ]
then
  echo "...plot dilatation..."
# plot shear strain rates
  awk 'NR > 2 {print $2,$1,$23*1}' $pth2strinfo >tmpdil
  # find min max and create cpt file
  T=`awk '{print $3}' tmpdil | gmt info -Eh `
  # set variables for scale plot
  scalevar_T ${T}
  Tmax=$(pythonc "print(round(${T},${Tmax_r})+${Tmax_r_marg})")
  T=`awk '{print $3}' tmpdil | gmt info -El `
  Tmin=$(pythonc "print(round(${T},${Tmax_r})-${Tmax_r_marg})")
  gmt makecpt -Cjet -T${Tmin}/${Tmax}/${cpt_step} > inx.cpt

  if [ "${GRDDAT}" -eq 0 ]
  then
    gmt pscontour tmpdil -R -J  -Cinx.cpt -I0.1 -O -K -V${VRBLEVM} >> ${outfile}
  else 
    gmt xyz2grd tmpdil -Gtmpdil.grd ${range} -I40m= -V
    gmt grdsample tmpdil.grd -I4s -Gtmpdil_sample.grd -V${VRBLEVM}
    gmt grdimage tmpdil_sample.grd ${proj} ${range} -Cinx.cpt -Q \
	-O -V${VRBLEVM} -K >> $outfile
  fi

  scale_step=$(pythonc "print(round(((${Tmax}-${Tmin})/5.),${scale_step_r}))")
  gmt pscoast -R -J -O -K -W0.25 -Df -Na  >> $outfile
  gmt psscale -Cinx.cpt -D8/-1.1/10/0.3h -B${scale_step}/:"nstrain/y": -I -S \
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
### PLOT SHEAR STRAIN RATES parameters
if [ "$SECINV" -eq 1 ]
then
  echo "...plot 2nd invariant..."
# plot shear strain rates
  awk 'NR > 2 {print $2,$1, $25*1}' $pth2strinfo >tmp2inv
  # find min max and create cpt file
  T=`awk '{print $3}' tmp2inv | gmt info -Eh `
# set variables for scale plot
  scalevar_T ${T}
  
  Tmax=$(pythonc "print(round(${T},${Tmax_r})+${Tmax_r_marg})")
  gmt makecpt -Cjet -T0/${Tmax}/${cpt_step} > inx.cpt
  
  if [ "${GRDDAT}" -eq 0 ]
  then
    gmt pscontour tmp2inv -R -J  -Cinx.cpt -I0.1 -O -K -V${VRBLEVM} >> ${outfile}
  else 
    gmt xyz2grd tmp2inv -Gtmp2inv.grd ${range} -I40m= -V
    gmt grdsample tmp2inv.grd -I4s -Gtmp2inv_sample.grd -V${VRBLEVM}
    gmt grdimage tmp2inv_sample.grd ${proj} ${range} -Cinx.cpt -Q \
	-O -V${VRBLEVM} -K >> $outfile
  fi

  scale_step=$(pythonc "print(round((${Tmax}/5.),${scale_step_r}))")
  gmt pscoast -R -J -O -K -W0.25 -Df -Na >> $outfile
  gmt psscale -Cinx.cpt -D8/-1.1/10/0.3h -B${scale_step}/:"nstrain/y": -I -S \
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
### PLOT STRAIN RATES parameters
if [ "$STRAIN" -eq 1 ]
then
  echo "...plot principal axes of strain rates..."
#   plot delaunay
  if [ "${DELTR}" -eq 1 ]
  then
    gmt psxy ${pth2deltr} -R -J -Wthinner -O -K -V${VRBLEVM} >> $outfile
  fi
  
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

# plot strain rates
  awk 'NR > 2 {print $2,$1,0,$17,$21+90}' $pth2strinfo \
  | gmt psvelo  -Jm $range -Sx${STRSC} -L -A5p+e -Gblue -W1p,blue -O -K -V${VRBLEVM} >> $outfile
  awk 'NR > 2 {print $2,$1,$15,0,$21+90}' $pth2strinfo \
  | gmt psvelo -Jm $range -Sx${STRSC} -L -A5p+e -Gred -W1p,red -O -K -V${VRBLEVM} >> $outfile
  
# plot scale of strain rates
  tmp_scrate=$(pythonc "print((${projscale}/150000000.)*20.)")
  strsclat=$(pythonc "print(${sclat} + ${tmp_scrate})")
  strsclon=$sclon

  echo "$strsclon $strsclat 0 -${strscmagn} 90" \
  | gmt psvelo -Jm $range -Sx${STRSC} -L -A10p+e -Gblue -W1.5p,blue \
        -O -K -V${VRBLEVM} >> $outfile
  echo "$strsclon $strsclat ${strscmagn} 0 90" \
  | gmt psvelo -Jm $range -Sx${STRSC} -L -A10p+e -Gred -W1.5p,red \
        -O -K -V${VRBLEVM} >> $outfile
  echo "$strsclon $strsclat 9 0 1 CB ${strscmagn} nstrain/y" \
  | gmt pstext -Jm -R -Dj0c/1c -Gwhite -O -K -V${VRBLEVM} >> $outfile
  
  sclat=${strasclat}
fi

# //////////////////////////////////////////////////////////////////////////////
### PLOT ROTATIONAL RATES parameters
if [ "$STRROT" -eq 1 ]
then
  echo "...plot rotational rates..."
#   plot delaunay triangles
  if [ "${DELTR}" -eq 1 ]
  then
    gmt psxy ${pth2deltr} -R -J -Wthinner -O -K -V${VRBLEVM} >> $outfile
  fi
  
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

# plot rotational rates
  ROT_wmag_sc=$(pythonc "print(0.206e0/${ROT_wedge_mag})")
  
  awk 'NR > 2 { if ($7 >= 0) print $2,$1,$7,$8}' $pth2strinfo \
  | gmt psvelo -Jm $range -Sw${ROTSC}/${ROT_wmag_sc} -Gred -E255/255/220 -L -A0.02  \
        -O -K -V${VRBLEVM} >> $outfile
  awk 'NR > 2 { if ($7 < 0) print $2,$1,$7,$8}' $pth2strinfo \
  | gmt psvelo -Jm $range -Sw${ROTSC}/${ROT_wmag_sc} -Gblue -E255/255/220 -L -A0.02  \
        -O -K -V${VRBLEVM} >> $outfile
        
# plot scale for rotational rates
  tmp_scrate=$(pythonc "print((${projscale}/150000000.)*20.)")
  rotsclat=$(pythonc "print(${sclat} + ${tmp_scrate})")
  rotsclon=$sclon
  
  ROT_wmagf=$(pythonc "print(${ROT_wedge_mag})")
  ROT_wmagf_sd=$(pythonc "print(${ROT_wedge_mag}/2.)")
  echo "$rotsclon $rotsclat ${ROT_wmagf} ${ROT_wmagf_sd}" \
  | gmt psvelo -Jm $range -Sw${ROTSC}/${ROT_wmag_sc} -Gred -E255/255/220 -L -A0.02  \
        -O -K -V${VRBLEVM} >> $outfile
  echo "$rotsclon $rotsclat -${ROT_wmagf} ${ROT_wmagf_sd}" \
  | gmt psvelo -Jm $range -Sw${ROTSC}/${ROT_wmag_sc} -Gblue -E255/255/220 -L -A0.02 \
        -O -K -V${VRBLEVM} >> $outfile
  echo "$rotsclon $rotsclat 9 0 1 CB ${ROT_wmagf} \261 ${ROT_wmagf_sd} \260/Myr" \
  | gmt pstext -Jm -R -Dj0c/-.6c -Gwhite -O -K -V${VRBLEVM} >> $outfile
  
  sclat=${rotsclat}
fi

# //////////////////////////////////////////////////////////////////////////////
### PLOT DEXTRAL SINISTRAL AXES OF MAXIMUM SHEAR STRAIN RATES parameters
if [ "$GTOTALAXES" -eq 1 ]
then
  echo "...plot dextral and sinistral maximum shear strain axes..."
  #   plot delaunay
  if [ "${DELTR}" -eq 1 ]
  then
    gmt psxy ${pth2deltr} -R -J -Wthinner -O -K -V${VRBLEVM} >> $outfile
  fi
  
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

# plot strain rates
  # dextral 
  awk 'NR > 2 {print $2,$1,0,$19,$21-45+90}' $pth2strinfo \
  | gmt psvelo  -Jm $range -Sx${STRSC} -L -A.1p+e -Gred -W1.5p,red \
        -O -K -V${VRBLEVM} >> $outfile
  # sinistral
  awk 'NR > 2 {print $2,$1,$19,0,$21-45+90}' $pth2strinfo \
  | gmt psvelo -Jm $range -Sx${STRSC} -L -A.1p+e -G0/204/0 -W1.5p,0/204/0 \
        -O -K -V${VRBLEVM} >> $outfile

# plot scale of strain rates
  tmp_scrate=$(pythonc "print((${projscale}/150000000.)*20.)")
  totsclat=$(pythonc "print(${sclat} + ${tmp_scrate})")
  totsclon=$sclon

  echo "$totsclon $totsclat 0 -${strscmagn} -45" \
  | gmt psvelo -Jm $range -Sx${STRSC} -L -A.1p+e -G0/204/0 -W1.5p,0/204/0 \
        -O -K -V${VRBLEVM} >> $outfile
  echo "$totsclon $totsclat ${strscmagn} 0 -45" \
  | gmt psvelo -Jm $range -Sx${STRSC} -L -A.1p+e -Gred -W1.5p,red \
        -O -K -V${VRBLEVM} >> $outfile
  echo "$totsclon $totsclat 9,1,black 45 RB dextral" \
  | gmt pstext -Jm -R -Dj0.2c/0.1c -F+f+a+j  -O -K -V${VRBLEVM} >> $outfile
  echo "$totsclon $totsclat 9,1,black -45 LB sinistral" \
  | gmt pstext -Jm -R -Dj0.2c/0.1c -F+f+a+j  -O -K -V${VRBLEVM} >> $outfile
  echo "$totsclon $totsclat 9,1,black 0 CB ${strscmagn} nstrain/y" \
  | gmt pstext -Jm -R -Dj0c/1c -F+f+a+j  -O -K -V${VRBLEVM} >> $outfile
fi


# //////////////////////////////////////////////////////////////////////////////
#  Plot custom logo
if [ "$LOGO" -eq 1 ]
then
  gmt psimage $pth2logos -O $logo_pos2 -W1.1c -F0.4  -K >>$outfile
fi

# //////////////////////////////////////////////////////////////////////////////
# FINAL SECTION
# add stamp for StrainTool
echo "$west $south 8,0,black 0 LB This image was produced using" \
  | gmt pstext -Jm -R -Dj0.1c/1.1c -F+f+a+j -K  -O -V${VRBLEVM} >> $outfile
echo "$west $south 9,1,white 0 LB STRAINTOOL for EPOS" \
  | gmt pstext -Jm -R -Dj0.2c/.65c -F+f+a+j -G165/0/236 -U$logo_pos -O -V${VRBLEVM} >> $outfile


#################--- Convert to other format ----###############################
if [ "$OUTJPG" -eq 1 ]
then
  echo "...adjust and convert to JPEG format..."
#   gs -sDEVICE=jpeg -dJPEGQ=100 -dNOPAUSE -dBATCH -dSAFER -r300 -sOutputFile=test.jpg ${outfile}
  gmt psconvert ${outfile} -A0.2c -Tj -V${VRBLEVM} 
fi

# clear all teporary files
echo "...remove temporary files..."
rm -rf tmp* gmt.conf gmt.history inx.cpt 2>/dev/null

# Print exit status
echo "[STATUS] Finished. Exit status: $?"

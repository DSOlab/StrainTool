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
#    filename       : gmtsttatsnplot.sh
#                     NAME=gmtstrainplot
#    version        : v-1.0
#                     VERSION=v1.0
#                     RELEASE=beta
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
#    update list    : LAST_UPDATE=JUL-2018
#    contact        : Demitris Anastasiou (dganastasiou@gmail.com)
#                     Xanthos Papanikolaou (xanthos@mail.ntua.gr)
#    ----------------------------------------------------------------------
# ==============================================================================
# //////////////////////////////////////////////////////////////////////////////
# HELP FUNCTION
function help {
	echo "/*****************************************************************/"
	echo " Program Name : gmtstrainplot.sh"
	echo " Version : v-1.0"
	echo " Purpose : Plot maps for gmtstrainplot"
	echo " Usage   : gmtstrainplot.sh -r  |  | -o [output] | -jpg "
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
	echo "/*** PLOT SSTATISTICS ********************************************/"
	echo "     -stats (input file) set input file"
	echo "     --stats-stations : plot used stations"
	echo "     --stats-doptimal : plot optimal distance (D)"
	echo "     --stats-sigma : plot sigma "
# 	echo "     -vhor (station_file)[:= horizontal velocities]  "
# 	echo "     -vsc [:=velocity scale] change velocity scale default 0.05"
# 	echo ""
# 	echo "/*** PLOT STRAINS **********************************************/"
# 	echo "     -str (strain file)[:= strains] Plot strain rates "
# 	echo "     -rot (strain file)[:= rots] Plot rotational rates "
# 	echo "     -gtot(strain file)[:=shear strain] plot total shear strain rate contours"
# 	echo "     -gtotaxes (strain file) dextral and sinistral max shear strain rates"
# 	echo "     -dil (strainfile)[:= dilatation] Plot dilatation and principal axes"
# 	echo "     -secinv (strain file) [:=2nd invariand] Plot second invariand"
# 	echo "     -strsc [:=strain scale]"
# 	echo "     -rotsc [:=rotational scales]"
# 	echo "  *for -gtot | -dil | -secinv use +grd to plot gridded data"
# 	echo "        ex:-gtot+grd "
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
# #BASH settings
# set -o errexit
# set -o pipefail
# set -o nounset
# set -o xtrace

# //////////////////////////////////////////////////////////////////////////////
# pre define parameters

# program version
VERSION="v.1.0-beta4.0"

# verbosity level for GMT, see http://gmt.soest.hawaii.edu/doc/latest/gmt.html#v-full
export VRBLEVM=n

# //////////////////////////////////////////////////////////////////////////////
# Source function files

# //////////////////////////////////////////////////////////////////////////////
# GMT parameters
gmt gmtset MAP_FRAME_TYPE fancy
gmt gmtset PS_PAGE_ORIENTATION portrait
gmt gmtset FONT_ANNOT_PRIMARY 8 FONT_LABEL 8 MAP_FRAME_WIDTH 0.12c FONT_TITLE 18p,Palatino-BoldItalic
gmt gmtset PS_MEDIA 30cx30c

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
# VHORIZONTAL=0
# STRAIN=0
# STRROT=0
# GTOTAL=0
# GTOTALAXES=0
# DILATATION=0
# SECINV=0
# GRDDAT=0
# MAX_STR_VALUE=1000000

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
	shift
	;;
    --stats-doptimal)
	STATS_DOPTIMAL=1
	shift
	;;
    --stats-sigma)
	STATS_SIGMA=1
	shift
	;;
    -topo)
  # switch topo not used in server!
	TOPOGRAPHY=1
	shift
	;;
    -o | --output)
	outfile=${2}.eps
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
	echo "version: "$VERSION
	exit 1
	shift
	;;
    *)
      echo "[ERROR] Bad argument structure. argument \"${1}\" is not right"
      echo "[STATUS] Script Finished Unsuccesful! Exit Status 1"
      exit 1
  esac
done

# //////////////////////////////////////////////////////////////////////////////
# check if files exist

###check dems
if [ "$TOPOGRAPHY" -eq 1 ]
then
  if [ ! -f $inputTopoB ]
  then
    echo "[WARNING] grd file for topography toes not exist, var turn to coastline"
    TOPOGRAPHY=0
  fi
fi

##check inputfiles
if [ "$PSTA" -eq 1 ]
then
  if [ ! -f $pth2sta ]
  then
    echo "[WARNING] input file $pth2sta does not exist"
    echo "          please download it and then use this switch"
    PSTA=0
    exit 1
  fi
fi

##check inputfiles
if [ "$DELTR" -eq 1 ]
then
  if [ ! -f $pth2deltr ]
  then
    echo "[WARNING] input file $pth2deltr does not exist"
    echo "          please download it and then use this switch"
    PSTA=0
    exit 1
  fi
fi

# //////////////////////////////////////////////////////////////////////////////
# READ STATISTICS FILE
if [ "$STATS" -eq 1 ]
then
  stat_x_grid_step=$(grep x_grid_step $pth2stats | awk '{print $3}')
  stat_y_grid_step=$(grep y_grid_step $pth2stats | awk '{print $3}')
  
  
  
#   calculate new variables
  west_grd=$(python -c "print(${west} + (${stat_x_grid_step}/2))")
  south_grd=$(python -c "print(${south} + (${stat_y_grid_step}/2))")
  range_grd="-R$west_grd/$east/$south_grd/$north"

  istep_grd=$(python -c "print(${stat_x_grid_step}*60)")
fi

# //////////////////////////////////////////////////////////////////////////////
# SET REGION PROPERTIES
# tmp_scrate=$(python -c "print((${prjscale}/150000000.)*10.)")
tmp_scrate=$(python -c "print((${projscale}/150000000.)*10.)")
sclat=$(echo print ${south} + ${tmp_scrate} | python)

tmp_scrate=$(python -c "print((${projscale}/150000000.)*27.)")
sclon=$(echo print ${east} - ${tmp_scrate} | python)

tmp_msclat=$(python -c "print int((${south} + ${north})/2)")
tmp_msclon=$(python -c "print int((${west} + ${east})/2)")
export scale=-Lf${sclon}/${sclat}/${tmp_msclat}:${tmp_msclon}/${sclength}+l+jr
# scale="-Lf20/33.5/36:24/100+l+jr"
range="-R$west/$east/$south/$north"
proj="-Jm24/37/1:$projscale"

# ####################### TOPOGRAPHY ###########################
if [ "$TOPOGRAPHY" -eq 0 ]
then
  echo "...plot coatlines..."
  ################## Plot coastlines only ######################	
  gmt	psbasemap $range $proj  -B$frame:."$maptitle": -P -K > $outfile
  gmt	pscoast -R -J -O -K -W0.25 -G225 -Df -Na $scale -U$logo_pos >> $outfile
# 	pscoast -Jm -R -Df -W0.25p,black -G195  -U$logo_pos -K -O -V >> $outfile
# 	psbasemap -R -J -O -K --FONT_ANNOT_PRIMARY=10p $scale --FONT_LABEL=10p >> $outfile
fi

# if [ "$TOPOGRAPHY" -eq 1 ]
# then
#   echo "...plot topography dem..."
#   # ####################### TOPOGRAPHY ###########################
#   # bathymetry
#   gmt makecpt -Cgebco.cpt -T-7000/0/50 -Z > $bathcpt
#   gmt grdimage $inputTopoB $range $proj -C$bathcpt -K > $outfile
#   gmt pscoast $proj -P $range -Df -Gc -K -O >> $outfile
# 	# land
#   gmt makecpt -Cgray.cpt -T-6000/1800/50 -Z > $landcpt
#   gmt grdimage $inputTopoL $range $proj -C$landcpt  -K -O >> $outfile
#   gmt pscoast -R -J -O -K -Q >> $outfile
# 	#------- coastline -------------------------------------------
#   gmt psbasemap -R -J -O -K -B$frame:."$maptitle":  $scale >> $outfile
#   gmt pscoast -J -R -Df -W0.25p,black -K  -O -U$logo_pos >> $outfile
# fi
# 
# # //////////////////////////////////////////////////////////////////////////////
# #  PLOT NOA CATALOGUE FAULTS Ganas et.al, 2013
# if [ "$FAULTS" -eq 1 ]
# then
#   echo "...plot NOA FAULTS CATALOGUE Ganas et.al, 2013 ..."
#   gmt psxy $pth2faults -R -J -O -K  -W.5,204/102/0  >> $outfile
# fi
# 
# # //////////////////////////////////////////////////////////////////////////////
# #  PLOT HARVARD CMT catalogue
# if [ "$PCMT" -eq 1 ]
# then
# # gmt	makecpt -Cseis -T0/150/10 -Z > seis2.cpt
#   gmt psmeca harvardcat.cmt -R -J -Sm0.3c/0 -G247/207/136 -W0.25p -T0 -O -K >> $outfile
#   gmt psmeca papazachos.cmt -R -Jm -Sa0.3/0 -h1 -CP0.25 -G110 -K -O -P >> $outfile
# fi

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

# # //////////////////////////////////////////////////////////////////////////////
# ### PLOT HORIZONTAL VELOCITIES
# 
# if [ "$VHORIZONTAL" -eq 1 ]
# then
#   echo "...plot horizontal velocities..."
# # plot stations
#   if [ "$PSTA" -eq 1 ]
#   then
# 
#     awk 'NR > 2 {print $2,$3}' $pth2sta  \
#     | gmt psxy -R -J -W.1 -Sc.15c -Gyellow -O -K -V${VRBLEVM} >> $outfile
#     
#     if [ "$LABELS" -eq 1 ]
#     then
#       awk 'NR > 2 {print $2,$3, "7,1,black", 0, "RB", $1}' $pth2sta \
#       | gmt pstext -R -J -Dj0.1c/0.1c -F+f+a+j -O -K -V${VRBLEVM} >> ${outfile}
#     fi
#   fi
# 
#   awk 'NR > 2 {print $2,$3,$4,$5,$6,$7,0,$1}' $pth2stainfo \
#   | gmt psvelo -R -Jm -Se${VSC}/0.95/0 -W.5p,black -A.05p+e -Gblue \
#     -O -K -V${VRBLEVM} >> $outfile  # 205/133/63.
# 
#   awk 'NR > 2 {print $2,$3,$4,$5,$6,$7,0,$1}' $pth2stainfo \
#   | gmt psvelo -R -Jm -Se${VSC}/0/0 -W2p,blue -A10p+e -Gblue \
#     -O -K -V${VRBLEVM} >> $outfile  # 205/133/63.
# 
# ###scale
# # plot scale of strain rates
#   tmp_scrate=$(python -c "print((${projscale}/150000000.)*5.)")
#   velsclat=$(echo print ${sclat} + ${tmp_scrate} | python)
#   velsclon=$(echo print ${sclon} - ${tmp_scrate} | python)
# 
#   echo "$velsclon $velsclat ${vscmagn} 0 1 1 0 " \
#   | gmt psvelo -R -Jm -Se${VSC}/0.95/0 -W.5p,black -A.05p+e -Gblue \
#     -O -K -V${VRBLEVM} >> $outfile
#   
#   echo "$velsclon $velsclat ${vscmagn} 0 5 5 0 " \
#   | gmt psvelo -R -Jm -Se${VSC}/0/0 -W2p,blue -A10p+e -Gblue \
#     -O -K -V${VRBLEVM} >> $outfile
#   
#   echo "$sclon $velsclat 9,1,black 0 CB ${vscmagn} \261 1 mm/y" \
#   | gmt pstext -Jm -R -Dj0c/.5c -F+f+a+j  -O -K -V${VRBLEVM} >> $outfile
#   
#   sclat=${velsclat}
# #   sclon=${velsclon}
# fi



# //////////////////////////////////////////////////////////////////////////////
### PLOT STATIONS USED FOR EACH CELL
if [ "$STATS_STATIONS" -eq 1 ]
then
  echo "...plot stations used for each grid cell..."
# plot shear strain rates
  awk 'NR > 24 {print $1-.125,$2-.125,$3}' $pth2stats > tmpstations
  # find min max and create cpt file
  T=`awk '{print $3}' tmpstations | gmt info -Eh `
  Tmax=$(python -c "print(round(${T},-1)+1)")
  T=`awk '{print $3}' tmpstations | gmt info -El `
  Tmin=$(python -c "print(round(${T},-1)-1)")
  gmt makecpt -Cjet -T${Tmin}/${Tmax}/1 > inx.cpt  
  
    gmt xyz2grd tmpstations -Gtmpstations.grd ${range_grd} -I${istep_grd}m -V
    gmt grdsample tmpstations.grd -I30m -Gtmpstations_sample.grd -V${VRBLEVM}
    gmt grdimage tmpstations_sample.grd ${proj} ${range} -Cinx.cpt -Q \
	-O -K -V${VRBLEVM}>> $outfile

  scale_step=$(python -c "print(round(((${Tmax}-${Tmin})/5.),0))")
  gmt pscoast -R -J -O -K -W0.25 -Df -Na -U$logo_pos -V${VRBLEVM}>> $outfile
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
### PLOTOPTIMAL dISTANCE d
if [ "$STATS_DOPTIMAL" -eq 1 ]
then
  echo "...plot optimal distance D for each grid cell..."
# plot shear strain rates
  awk 'NR > 24 {print $1,$2,$4}' $pth2stats > tmpdoptimal
  # find min max and create cpt file
  T=`awk '{print $3}' tmpdoptimal | gmt info -Eh `
  Tmax=$(python -c "print(round(${T},-1)+1)")
  T=`awk '{print $3}' tmpdoptimal | gmt info -El `
  Tmin=$(python -c "print(round(${T},-1)-1)")
  gmt makecpt -Cjet -T${Tmin}/${Tmax}/1 > inx.cpt  
  
    gmt xyz2grd tmpdoptimal -Gtmpdoptimal.grd ${range_grd} -I${istep_grd}m= -V
    gmt grdsample tmpdoptimal.grd -I30m -Gtmpdoptimal_sample.grd -V${VRBLEVM}
    gmt grdimage tmpdoptimal_sample.grd ${proj} ${range} -Cinx.cpt -Q \
	-O -K -V${VRBLEVM}>> $outfile

  scale_step=$(python -c "print(round(((${Tmax}-${Tmin})/5.),-1))")
  gmt pscoast -R -J -O -K -W0.25 -Df -Na -U$logo_pos -V${VRBLEVM}>> $outfile
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
# plot shear strain rates
  awk 'NR > 24 {print $1,$2,$6}' $pth2stats > tmpsigma
  # find min max and create cpt file
  T=`awk '{print $3}' tmpsigma | gmt info -Eh `
  Tmax=$(python -c "print(round(${T},3)+.001)")
  T=`awk '{print $3}' tmpsigma | gmt info -El `
  Tmin=$(python -c "print(round(${T},3)-.001)")
  gmt makecpt -Cjet -T${Tmin}/${Tmax}/.001 > inx.cpt  
  
    gmt xyz2grd tmpsigma -Gtmpsigma.grd ${range_grd} -I${istep_grd}m= -V
    gmt grdsample tmpsigma.grd -I30m -Gtmpsigma_sample.grd -V${VRBLEVM}
    gmt grdimage tmpsigma_sample.grd ${proj} ${range} -Cinx.cpt -Q \
	-O -K -V${VRBLEVM}>> $outfile

  scale_step=$(python -c "print(round((${Tmax}/5.),3))")
  gmt pscoast -R -J -O -K -W0.25 -Df -Na -U$logo_pos -V${VRBLEVM}>> $outfile
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
  gmt pslegend .legend ${legendc} -C0.1c/0.1c -L1.3 -O -K >> $outfile
fi

# //////////////////////////////////////////////////////////////////////////////
#  Plot custom logo
if [ "$LOGO" -eq 1 ]
then
  gmt psimage $pth2logos -O $logo_pos2 -W1.1c -F0.4  -K >>$outfile
fi

# //////////////////////////////////////////////////////////////////////////////
# FINAL SECTION
#################--- Close ps output file ----##################################
#echo "909 909" | gmt psxy -Sc.1 -Jm -R  -W1,red -O -V${VRBLEVM} >> ${outfile}
echo "$west $south 8,0,black 0 LB This image was produced using" \
  | gmt pstext -Jm -R -Dj0.1c/1.1c -F+f+a+j -K  -O -V${VRBLEVM} >> $outfile
echo "$west $south 9,1,white 0 LB STRAINTOOL for EPOS" \
  | gmt pstext -Jm -R -Dj0.2c/.65c -F+f+a+j -G165/0/236 -O -V${VRBLEVM} >> $outfile


#################--- Convert to other format ----###############################
if [ "$OUTJPG" -eq 1 ]
then
  echo "...adjust and convert to JPEG format..."
#   gs -sDEVICE=jpeg -dJPEGQ=100 -dNOPAUSE -dBATCH -dSAFER -r300 -sOutputFile=test.jpg ${outfile}
  gmt psconvert ${outfile} -A0.2c -Tj -V${VRBLEVM} 
fi
# if [ "$OUTPNG" -eq 1 ]
# then
#   echo "...adjust and convert to PNG format..."
#   gmt psconvert ${outfile} -A0.2c -TG -V${VRBLEVM} 	
# fi
# if [ "$OUTEPS" -eq 1 ]
# then
#   echo "...adjust and convert to EPS format..."
#   gmt psconvert ${outfile} -A0.2c -Te -V${VRBLEVM} 
# fi
# if [ "$OUTPDF" -eq 1 ]
# then
#   echo "...adjust and convert to PDF format..."
#   gmt psconvert ${outfile} -A0.2c -Tf -V${VRBLEVM} 
# fi

# clear all teporary files
echo "...remove temporary files..."
rm -rf tmp* gmt.conf gmt.history *.legend inx.cpt

# Print exit status
echo "[STATUS] Finished. Exit status: $?"

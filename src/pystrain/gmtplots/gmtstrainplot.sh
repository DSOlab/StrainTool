#!/usr/bin/env bash

# //////////////////////////////////////////////////////////////////////////////
# ==============================================================================
#
#    |===========================================|
#    |**     DIONYSOS SATELLITE OBSERVATORY    **|
#    |**        HIGHER GEODESY LABORATORY      **|
#    |** National Tecnical University of Athens**|
#    |===========================================|
#
#    filename       : gmtstrainplot.sh
#                     NAME=gmtstrainplot
#    version        : v-1.0
#                     VERSION=v1.0
#                     RELEASE=beta
#    licence        : MIT
#    created        : MAY-2018
#    usage          :
#    GMT Modules    :
#    UNIX progs     :
#    exit code(s)   : 0 -> success
#                   : 1 -> error
#    discription    : 
#    uses           : 
#    notes          :
#    update list    : LAST_UPDATE=OCT-2017
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
	echo "/*** Basic Plots & Background ***********************************/"
	echo "           -r [:= region] region to plot west east south north projscale frame(default Greece)"
	echo "                   use: -r west east south north projscale frame"
# 	echo "           -faults [:= faults] plot NOA fault database"
# 	echo "           -pcmt [:= plot cmt] plot HARVARD cmt and papazachos fm"
# 	echo "           -topo [:= topography] plot topography"
	echo ""
	echo "/*** PLOT STATIONS ***********************************************/"
	echo "           -psta [:=stations] plot only stations from input file"
	echo ""
# 	echo "/*** PLOT VELOCITIES ********************************************/"
# 	echo "           -vhor (gmt_file)[:= horizontal velocities]  "
# 	echo "           -vsc [:=velocity scale] change valocity scale default 0.05"
# 	echo ""
	echo "/*** PLOT STRAINS **********************************************/"
	echo "           -str (gmt_file)[:= strains] Plot strain rates "
	echo "           -rot (gmt_file)[:= rots] Plot rotational rates "
	echo "           -gtot(strain file)[:=shear strain] plot total shear strain rate"
# 	echo "           -dil (gmt_file)[:= dilatation] Plot dilatation and principal axes"
	echo "           -strsc [:=strain scale]"
	echo "           -rotsc [:=rotational scales]"
	echo ""
	echo "/*** OTHER OPRTIONS ********************************************/"
	echo "           -o [:= output] name of output files"
	echo "           -l [:=labels] plot labels"
# 	echo "           -leg [:=legend] insert legends"
	echo "           -logo [:=logo] plot logo"
	echo "           -jpg : convert eps file to jpg"
	echo "           -h [:= help] help menu"
	echo " Exit Status:    1 -> help message or error"
	echo " Exit Status:  = 0 -> sucesseful exit"
	echo " run scr: ./gmtplot.sh -topo -jpg"
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
VERSION="v.1.0-beta1.0"

# verbosity level for GMT, see http://gmt.soest.hawaii.edu/doc/latest/gmt.html#v-full
# 
export VRBLEVM=n

# //////////////////////////////////////////////////////////////////////////////
# Source function files

# //////////////////////////////////////////////////////////////////////////////
# GMT parameters
gmt gmtset MAP_FRAME_TYPE fancy
gmt gmtset PS_PAGE_ORIENTATION portrait
gmt gmtset FONT_ANNOT_PRIMARY 8 FONT_LABEL 8 MAP_FRAME_WIDTH 0.12c FONT_TITLE 18p,Palatino-BoldItalic

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
VHORIZONTAL=0
STRAIN=0
STRROT=0
GTOTAL=0

##//////////////////check default param
if [ ! -f "default-param" ]
then
	echo "default-param file does not exist"
	exit 1
else
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
		-r)
			west=$2
			east=$3
			south=$4
			north=$5
			projscale=$6
			frame=$7
# 			REGION=$2
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
			pth2sta=../input/$2
			PSTA=1
			shift
			shift
			;;
		-vhor)
			pth2vhor=${pth2inptf}/$2
			VHORIZONTAL=1
			shift
			shift
			;;
		-vsc)
			VSC=$2
			shift
			shift
			;;
		-str)
			pth2strinfo=${pth2inptf}/${2}
# 			pth2comp=${pth2inptf}/${2}.comp
# 			pth2ext=${pth2inptf}/${2}.ext
			STRAIN=1
			shift
			shift
			;;
		-strsc)
			STRSC=$2
			shift
			shift
			;;
		-rot)
			pth2strinfo=${pth2inptf}/${2}
# 			pth2strain=${pth2inptf}/${2}par.str
			STRROT=1
			shift
			shift
			;;
		-rotsc)
			ROTSC=${2}
			shift
			shift
			;;
		-gtot)
			pth2strinfo=${pth2inptf}/${2}
# 			pth2strain=${pth2inptf}/${2}par.str
			GTOTAL=1
			shift
			shift
			;;
		-topo)
#                       switch topo not used in server!
			TOPOGRAPHY=1
			shift
			;;
		-faults)
			FAULTS=1
			shift
			;;	
		-pcmt)
			PCMT=1
			shift
			;;
		-o)
			outfile=${2}.eps
			out_jpg=${2}.jpg
			shift
			shift
			;;
		-l)
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
		-h)
			help
			;;
	esac
done


# //////////////////////////////////////////////////////////////////////////////
# check if files exist

###check dems
if [ "$TOPOGRAPHY" -eq 1 ]
then
	if [ ! -f $inputTopoB ]
	then
		echo "grd file for topography toes not exist, var turn to coastline"
		TOPOGRAPHY=0
	fi
fi

##check inputfiles
if [ "$PSTA" -eq 1 ]
then
	if [ ! -f $pth2sta ]
	then
		echo "input file $pth2sta does not exist"
		echo "please download it and then use this switch"
		PSTA=0
		exit 1
	fi
fi

if [ "$VHORIZONTAL" -eq 1 ]
then
	if [ ! -f $pth2vhor ]
	then
		echo "input file $pth2vhor does not exist"
		echo "please download it and then use this switch"
		VHORIZONTAL=0
		exit 1
	fi
fi

if [ "$STRAIN" -eq 1 ]
then
	if [ ! -f $pth2strinfo ]
	then
		echo "input file $pth2strinfo does not exist"
		echo "please download it and then use this switch"
		STRAIN=0
		exit 1
	fi
fi

if [ "$STRROT" -eq 1 ]
then
# echo "no rotation plots yet"
	if [ ! -f $pth2strain ]
	then
		echo "input file $pth2strain does not exist"
		echo "please download it and then use this switch"
		STRROT=0
		exit 1
	fi
fi

###check NOA FAULT catalogue
if [ "$FAULTS" -eq 1 ]
then
	if [ ! -f $pth2faults ]
	then
		echo "NOA Faults database does not exist"
		echo "please download it and then use this switch"
		FAULTS=0
	fi
fi

###check LOGO file
if [ ! -f "$pth2logos" ]
then
	echo "Logo file does not exist"
	LOGO=0
fi




# //////////////////////////////////////////////////////////////////////////////
# SET REGION PROPERTIES
gmt	gmtset PS_MEDIA 22cx22c
# 	scale="-Lf20/33.5/36:24/100+l+jr"
	range="-R$west/$east/$south/$north"
	proj="-Jm24/37/1:$projscale"

# ####################### TOPOGRAPHY ###########################
if [ "$TOPOGRAPHY" -eq 0 ]
then
	################## Plot coastlines only ######################	
gmt	psbasemap $range $proj $scale -B$frame:."$maptitle": -P -K > $outfile
gmt	pscoast -R -J -O -K -W0.25 -G225 -Df -Na -U$logo_pos >> $outfile
# 	pscoast -Jm -R -Df -W0.25p,black -G195  -U$logo_pos -K -O -V >> $outfile
# 	psbasemap -R -J -O -K --FONT_ANNOT_PRIMARY=10p $scale --FONT_LABEL=10p >> $outfile
fi
if [ "$TOPOGRAPHY" -eq 1 ]
then
	# ####################### TOPOGRAPHY ###########################
	# bathymetry
gmt	makecpt -Cgebco.cpt -T-7000/0/50 -Z > $bathcpt
gmt	grdimage $inputTopoB $range $proj -C$bathcpt -K > $outfile
gmt	pscoast $proj -P $range -Df -Gc -K -O >> $outfile
	# land
gmt	makecpt -Cgray.cpt -T-6000/1800/50 -Z > $landcpt
gmt	grdimage $inputTopoL $range $proj -C$landcpt  -K -O >> $outfile
gmt	pscoast -R -J -O -K -Q >> $outfile
	#------- coastline -------------------------------------------
gmt	psbasemap -R -J -O -K -B$frame:."$maptitle":  $scale >> $outfile
gmt	pscoast -J -R -Df -W0.25p,black -K  -O -U$logo_pos >> $outfile
fi

#////////////////////////////////////////////////////////////////
#  PLOT NOA CATALOGUE FAULTS Ganas et.al, 2013
if [ "$FAULTS" -eq 1 ]
then
	echo "plot NOA FAULTS CATALOGUE Ganas et.al, 2013 ..."
gmt	psxy $pth2faults -R -J -O -K  -W.5,204/102/0  >> $outfile
fi

#////////////////////////////////////////////////////////////////
#  PLOT HARVARD CMT catalogue
if [ "$PCMT" -eq 1 ]
then
# gmt	makecpt -Cseis -T0/150/10 -Z > seis2.cpt
gmt	psmeca harvardcat.cmt -R -J -Sm0.3c/0 -G247/207/136 -W0.25p -T0 -O -K >> $outfile

gmt	psmeca papazachos.cmt -R -Jm -Sa0.3/0 -h1 -CP0.25 -G110 -K -O -P >> $outfile
fi
#////////////////////////////////////////////////////////////////
### PLOT STATIONS

if [ "$PSTA" -eq 1 ]
then
# 	awk 'NR != 1 {print $1,$2,$3,$4,0,0,0,$8}' $pth2vhor | gmt psvelo -R -Jm -Se${VSC}/0.95/0 -W2p,blue -A10p+e -Gblue -O -K -L -V >> $outfile  # 205/133/63.
# 	awk -F, '{print $3, $2}' $pth2sta | gmt psxy -Jm -O -R -St0.17c -W0.01c -Gred -K >> $outfile
	awk -F, '{print $1, $2}' $pth2inptf/ionpvel_vhor.sta | gmt psxy -Jm -O -R -St0.27c -W0.01c -Gblue -K >> $outfile
	awk -F, '{print $1, $2}' $pth2inptf/ioncvel_vhor.sta | gmt psxy -Jm -O -R -St0.22c -W0.01c -Gred -K >> $outfile


	if [ "$LABELS" -eq 1 ]
	then
# 		 awk -F, '{print $3,$2,7,0,1,"LM",$1}' $pth2sta | grep -v "PYLO" | gmt pstext -Jm -R -Dj0.2c/0.2c -Gwhite -O -K -V>> $outfile
	gmt pstext $pth2inptf/ionpvel_vhor.sta -Jm -R -Dj0.2c/0.2c -Gwhite -O -K -V>> $outfile
	gmt pstext $pth2inptf/ioncvel_vhor.sta -Jm -R -Dj0.2c/0.2c -Gwhite -O -K -V>> $outfile
	fi

echo "22.3 37.72" | gmt psxy -Jm -O -R -St0.27c -W0.01c -Gblue -K >> $outfile
echo "22.3 37.72 8 0 1 LM Permanent Sites" | gmt pstext -Jm -R -Dj0.4c/0.3c -Gwhite -O -K -V>> $outfile
echo "22.3 37.65" | gmt psxy -Jm -O -R -St0.22c -W0.01c -Gred -K >> $outfile
echo "22.3 37.65 8 0 1 LM Campaign Sites" | gmt pstext -Jm -R -Dj0.4c/0.3c -Gwhite -O -K -V>> $outfile
###scale
# echo "$vsclon $vsclat $vscmagn 0 0 0 0 $vscmagn mm" | gmt psvelo -R -Jm -Se${VSC}/0.95/10 -W2p,blue -A10p+e -Gblue -O -K -L -V >> $outfile

fi



#////////////////////////////////////////////////////////////////
### PLOT HORIZONTAL VELOCITIES

if [ "$VHORIZONTAL" -eq 1 ]
then
	awk -F, '{print $1, $2}' $pth2inptf/ionpvel_vhor.sta | gmt psxy -Jm -O -R -Sc0.17c -W0.001c -G250 -K >> $outfile
	awk -F, '{print $1, $2}' $pth2inptf/ioncvel_vhor.sta | gmt psxy -Jm -O -R -Sc0.17c -W0.001c -G250 -K >> $outfile

# 	awk 'NR != 1 {print $1,$2,$3,$4,0,0,0,$8}' $pth2vhor | gmt psvelo -R -Jm -Se${VSC}/0.95/0 -W2p,blue -A10p+e -Gblue -O -K -L -V >> $outfile  # 205/133/63.
	awk 'NR != 1 {print $1,$2,$3,$4,$5,$6,0,$8}' $pth2inptf/ionpvel_vhor.vel | gmt psvelo -R -Jm -Se${VSC}/0.75/0 -W.5p,50 -A10p+e -Gblue -O -K -L -V >> $outfile  # 205/133/63.
	awk 'NR != 1 {print $1,$2,$3,$4,$5,$6,0,$8}' $pth2inptf/ioncvel_vhor.vel | gmt psvelo -R -Jm -Se${VSC}/0.75/0 -W.5p,50 -A8p+e -Gred -O -K -L -V >> $outfile  # 205/133/63.

	awk 'NR != 1 {print $1,$2,$3,$4,0,0,0,$8}' $pth2inptf/ionpvel_vhor.vel | gmt psvelo -R -Jm -Se${VSC}/0.95/0 -W2p,blue -A10p+e -Gblue -O -K -L -V >> $outfile  # 205/133/63.
	awk 'NR != 1 {print $1,$2,$3,$4,0,0,0,$8}' $pth2inptf/ioncvel_vhor.vel | gmt psvelo -R -Jm -Se${VSC}/0.95/0 -W1.5p,red -A8p+e -Gred -O -K -L -V >> $outfile  # 205/133/63.

	if [ "$LABELS" -eq 1 ]
	then
# 		 awk '{print $1,$2,8,0,1,"LM",$8}' $pth2vhor | gmt pstext -Jm -R -Dj0.2c/0.2c -Gwhite -O -K -V>> $outfile
 	gmt pstext $pth2inptf/ionpvel_vhor.sta -Jm -R -Dj0.2c/0.2c -Gwhite -O -K -V>> $outfile
	gmt pstext $pth2inptf/ioncvel_vhor.sta -Jm -R -Dj0.2c/0.2c -Gwhite -O -K -V>> $outfile

	fi

###scale
# echo "$vsclon $vsclat $vscmagn 0 1 1 0 $vscmagn mm" | gmt psvelo -R -Jm -Se${VSC}/0.95/0 -W.5p,50 -A10p+e -Gblue -O -K -L -V >> $outfile
# echo "$vsclon $vsclat $vscmagn 0 0 0 0 $vscmagn mm" | gmt psvelo -R -Jm -Se${VSC}/0.95/0 -W2p,blue -A10p+e -Gblue -O -K -L -V >> $outfile
# echo "$vsclon $vsclat 9 0 1 LB $vscmagn \261 1 mm/y" | gmt pstext -Jm -R -Dj-.3c/0.5c  -O -K -V>> $outfile
echo "22.48 37.70 $vscmagn 0 1 1 0 $vscmagn mm" | gmt psvelo -R -Jm -Se${VSC}/0.95/0 -W.5p,50 -A10p+e -Gblue -O -K -L -V >> $outfile
echo "22.48 37.70 $vscmagn 0 0 0 0 $vscmagn mm" | gmt psvelo -R -Jm -Se${VSC}/0.95/0 -W2p,blue -A10p+e -Gblue -O -K -L -V >> $outfile

echo "22.48 37.65 $vscmagn 0 1 1 0 $vscmagn mm" | gmt psvelo -R -Jm -Se${VSC}/0.95/0 -W.5p,50 -A10p+e -Gblue -O -K -L -V >> $outfile
echo "22.48 37.65 $vscmagn 0 0 0 0 $vscmagn mm" | gmt psvelo -R -Jm -Se${VSC}/0.95/0 -W2p,red -A10p+e -Gblue -O -K -L -V >> $outfile
echo "22.45 37.69 8 0 1 LB $vscmagn \261 1 mm/y" | gmt pstext -Jm -R -Dj-.3c/0.5c  -O -K -V>> $outfile
echo "22.48 37.70 8 0 1 RM Permanent Sites" | gmt pstext -Jm -R -Dj0.4c/0.3c -Gwhite -O -K -V>> $outfile
echo "22.48 37.65 8 0 1 RM Campaign Sites" | gmt pstext -Jm -R -Dj0.4c/0.3c -Gwhite -O -K -V>> $outfile

fi

#////////////////////////////////////////////////////////////////
### PLOT SHEAR STRAIN RATES parameters
if [ "$GTOTAL" -eq 1 ]
then
awk '{if ($11 < 300) print $2,$1,$11}' $pth2strinfo >tmpgtot
gmt makecpt -Cjet -T0/300/.5 > inx.cpt
# gmt pscontour tmpgtot -R -J -Wthin -Cinx.cpt -G.1i/0 -O -K >> $outfile
gmt pscontour tmpgtot -R -J -I -Cinx.cpt -O -K >> $outfile

# pscoast -J -R -W -Di -O -K -UBL/3.8c/-3.2c/"DSO-HGL/NTUA" >> $ps
gmt pscoast -R -J -O -K -W0.25 -Df -Na -U$logo_pos >> $outfile
gmt psscale -Cinx.cpt -D.9/3/3.2/0.3 -B100/:"nstrain/y": -I -S -O -K >> $outfile
# for i in `ls $pth2work*.sta`;do
#     gmt psxy $i -Jm -O -R -Sc.17c -Gyellow -W0.01c  -K >> $outfile ;   #fill patterns-Gp300/1:F0/0/0B-ls
# #     gmt pstext $i -Jm -R -Dj0.2c/0.2c -Gwhite -O -K -V>> $outfile
# 
# done
# for i in `ls $pth2work*.sta`;do
#     gmt psxy $i -Jm -O -R -Sc.17c -Gyellow -W0.01c  -K >> $outfile ;   #fill patterns-Gp300/1:F0/0/0B-ls
# #     gmt pstext $i -Jm -R -Dj0.2c/0.2c -Gwhite -O -K -V>> $outfile
# 
# done
# 
# for i in `ls $pth2work*.reg`;do
#     gmt psxy $i -Jm -O -R  -W0.02c,90  -K >> $outfile ;   #fill patterns-Gp300/1:F0/0/0B-ls
# done
# 
#     gmt psvelo  $pth2comp -Jm $range -Sx${STRSC} -L -A10p+e -Gblue -W2p,blue -V -K -O>> $outfile
#     gmt psvelo $pth2ext -Jm $range -Sx${STRSC} -L -A10p+e -Gred -W2p,red -V -K -O>> $outfile
#     
#     echo "$strsclon $strsclat 0 -.2 90" | gmt psvelo -Jm $range -Sx${STRSC} -L -A10p+e -Gblue -W2p,blue -V  -K -O>> $outfile
#     echo "$strsclon $strsclat .2 0 90" | gmt psvelo -Jm $range -Sx${STRSC} -L -A10p+e -Gred -W2p,red -V  -K -O>> $outfile
#     echo "$strsclon $strsclat 9 0 1 CB 200 nstrain/y" | gmt pstext -Jm -R -Dj0c/1c -Gwhite -O -K -V>> $outfile

#////////////////////////////////////////////////////////////////
#  PLOT HARVARD CMT catalogue
if [ "$PCMT" -eq 1 ]
then
# gmt	makecpt -Cseis -T0/150/10 -Z > seis2.cpt
gmt	psmeca harvardcat.cmt -R -J -Sm0.3c/0 -G247/207/136 -W0.25p -T0 -O -K >> $outfile

gmt	psmeca papazachos.cmt -R -Jm -Sa0.3/0 -h1 -CP0.25 -G110 -K -O -P >> $outfile
fi
fi

#////////////////////////////////////////////////////////////////
### PLOT STRAIN RATES parameters

if [ "$STRAIN" -eq 1 ]
then
# for i in `ls $pth2work*.sta`;do
#     gmt psxy $i -Jm -O -R -Sc.16c -Gyellow -W0.01c  -K >> $outfile ;   #fill patterns-Gp300/1:F0/0/0B-ls
# #     gmt pstext $i -Jm -R -Dj0.2c/0.2c -Gwhite -O -K -V>> $outfile
# 
# done
# 
# for i in `ls $pth2work*.reg`;do
#     gmt psxy $i -Jm -O -R  -W0.02c,90  -K >> $outfile ;   #fill patterns-Gp300/1:F0/0/0B-ls
# done

    awk 'NR !=1 { if ($9 < 4000 && $10 >-4000) print $2,$1,0,$10,$12+90}' $pth2strinfo |\
    gmt psvelo  -Jm $range -Sx${STRSC} -L -A10p+e -Gblue -W2p,blue -O -K -V${VRBLEVM} >> $outfile
    awk 'NR !=1 { if ($9 < 4000 && $10 >-4000) print $2,$1,$9,0,$12+90}' $pth2strinfo |\
    gmt psvelo -Jm $range -Sx${STRSC} -L -A10p+e -Gred -W2p,red -O -K -V${VRBLEVM} >> $outfile
    
    echo "$strsclon $strsclat 0 -.15 90" | gmt psvelo -Jm $range -Sx${STRSC} -L -A10p+e -Gblue -W2p,blue -V  -K -O>> $outfile
    echo "$strsclon $strsclat .15 0 90" | gmt psvelo -Jm $range -Sx${STRSC} -L -A10p+e -Gred -W2p,red -V  -K -O>> $outfile
    echo "$strsclon $strsclat 9 0 1 CB 150 nstrain/y" | gmt pstext -Jm -R -Dj0c/1c -Gwhite -O -K -V>> $outfile

fi

#////////////////////////////////////////////////////////////////
### PLOT ROTATIONAL RATES parameters
if [ "$STRROT" -eq 1 ]
then
# for i in `ls $pth2work*.sta`;do
#     gmt psxy $i -Jm -O -R -Sc.16c -Gyellow -W0.01c  -K >> $outfile ;   #fill patterns-Gp300/1:F0/0/0B-ls
# #     gmt pstext $i -Jm -R -Dj0.2c/0.2c -Gwhite -O -K -V>> $outfile
# 
# done
# 
# for i in `ls $pth2work*.reg`;do
#     gmt psxy $i -Jm -O -R  -W0.02c  -K >> $outfile ;   #fill patterns-Gp300/1:F0/0/0B-ls
# done

	awk 'NR !=1 {if ($6>=0) print $2,$1,$6/1000000000,0.00000001}' $pth2strinfo |\
	gmt psvelo -Jm $range -Sw${ROTSC}/1.e7 -Gred -E0/0/0/10 -L -A0.02  -O -K -V${VRBLEVM} >> $outfile
	awk 'NR !=1 {if ($6<0) print $2,$1,$6/1000000000,0.00000001}' $pth2strinfo |\
	gmt psvelo -Jm $range -Sw${ROTSC}/1.e7 -Gblue -E0/0/0/10 -L -A0.02  -O -K -V${VRBLEVM} >> $outfile

	echo "$strsclon $strsclat 0.00000005 0.00000001"| gmt psvelo -Jm $range -Sw${ROTSC}/1.e7 -Gred -E0/0/0/10 -L -A0.02  -V -K -O>> $outfile
	echo "$strsclon $strsclat -0.00000005 0.00000001"| gmt psvelo -Jm $range -Sw${ROTSC}/1.e7 -Gblue -E0/0/0/10 -L -A0.02  -V -K -O>> $outfile
	echo "$strsclon $strsclat 9 0 1 CB -0.05 / 0.05 ppm" | gmt pstext -Jm -R -Dj0c/-.6c -Gwhite -O -K -V>> $outfile

# psvelo <<END -Jm $range -Sw2/1.e3 -Gred -E0/0/0/10 -L -A0.05/0/0  -V -K -O>> $outfile
# 20.53726 38.33053 -0.001697 0.0001
# END
fi



# ///////////////// PLOT LEGEND \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/
if [ "$LEGEND" -eq 1 ]
then
gmt	pslegend .legend ${legendc} -C0.1c/0.1c -L1.3 -O -K >> $outfile
fi

#/////////////////PLOT LOGO DSO\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/
if [ "$LOGO" -eq 1 ]
then
gmt	psimage $pth2logos -O $logo_pos2 -W1.1c -F0.4  -K >>$outfile
fi

# //////////////////////////////////////////////////////////////////////////////
# FINAL SECTION
#################--- Close ps output file ----##################################
echo "909 909" | gmt psxy -Sc.1 -Jm -R  -W1,red -O -V${VRBLEVM} >> ${outfile}

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
rm -rf tmp* gmt.conf gmt.history *.legend inx.cpt
# Print exit status
echo "[STATUS] Finished. Exit status: $?"


# NOA FAULTS reference
# Ganas Athanassios, Oikonomou Athanassia I., and Tsimi Christina, 2013. NOAFAULTS: a digital database for active faults in Greece. Bulletin of
#  the Geological Society of Greece, vol. XLVII and Proceedings of the 13th International Congress, Chania, Sept. 2013.
#!/usr/bin/env bash
# program version
VERSION="plotall - v1.0"
# //////////////////////////////////////////////////////////////////////////////
# ==============================================================================
#
#    |===============================================|
#    |**       DIONYSOS SATELLITE OBSERVATORY      **|
#    |**          HIGHER GEODESY LABORATORY        **|
#    |**  National Technical University of Athens  **|
#    |===============================================|
#
#    filename       : plotall.sh
#                     NAME=plotall
#    version        : 1.0
#                     VERSION=1.0
#                     RELEASE=1.0
#    licence        : MIT
#    created        : SEP-2018
#    usage          :
#    GMT Modules    : 
#    UNIX progs     : 
#    exit code(s)   : 0 -> success
#                   : 1 -> error
#    discription    : Plot all maps for StrainTool
#    uses           : 
#    notes          :
#    update list    : LAST_UPDATE=DEC-2018
#    contact        : Dimitris Anastasiou (dganastasiou@gmail.com)
#                     Xanthos Papanikolaou (xanthos@mail.ntua.gr)
#    ----------------------------------------------------------------------
# ==============================================================================

function help {
	echo "/*****************************************************************/"
	echo " Program Name : plotall.sh"
	echo " Version : ${VERSION}"
	echo " Purpose : Plot alla available maps using prefix/suffix work id"
	echo " Usage   : plotall.sh "
	echo " Switches: "
	echo " no switch produce default output names"
	echo "   -p | --prefix <workid>: add prefix work id to output file"
	echo "   -s | --suffix <workid>: add suffix work id to output file"
	echo "   -h | --help: help panel"
    echo "/*****************************************************************/"
exit 0
}

# //////////////////////////////////////////////////////////////////////////////
# BASH settings
# set -o errexit
set -e
set -o pipefail
# set -o nounset
# set -o xtrace

prefix=""
suffix=""

if [ "$#" == "0" ]
then
  prefix=""
  suffix=""
else
while [ $# -gt 0 ]
do
  case "$1" in
    -p | --prefix)
      prefix="${2}-"
      shift
      shift
      ;;
    -s | --suffix)
      suffix="-${2}"
      shift
      shift
      ;;
    -h | --help)
      help
      ;;
    -v | --version)
      echo "Version: "$VERSION
      exit 1
      shift
      ;;
    *)
      echo "[ERROR] Bad argument structure. argument \"${1}\" is not right"
      echo "[STATUS] Script Finished Unsuccesfully! Exit Status 1"
      exit 1
  esac
done
fi

echo "[PLOTALL] ...plot horizontal velocities..."
./gmtstrainplot.sh -jpg -vhor station_info.dat -psta -l -o ${prefix}output_vel${suffix}
echo "[PLOTALL] ...plot principal axes of strain rates..."
./gmtstrainplot.sh -jpg -str strain_info.dat -psta -l -o ${prefix}output_str${suffix}
echo "[PLOTALL] ...plot rotational rates..."
./gmtstrainplot.sh -jpg -rot strain_info.dat -psta -l -o ${prefix}output_rot${suffix}
echo "[PLOTALL] ...plot dextral, sinistral maximum shear strain rates..."
./gmtstrainplot.sh -jpg -gtotaxes strain_info.dat -psta -l -o ${prefix}output_gtotaxes${suffix}
echo "[PLOTALL] ...plot maximum shear strain..."
./gmtstrainplot.sh -jpg -gtot strain_info.dat -psta -l -o ${prefix}output_gtot${suffix}
echo "[PLOTALL] ...plot dilatation..."
./gmtstrainplot.sh -jpg -dil strain_info.dat -psta -l -o ${prefix}output_dil${suffix}
echo "[PLOTALL] ...plot second invariant..."
./gmtstrainplot.sh -jpg -secinv strain_info.dat -psta -l -o ${prefix}output_2inv${suffix}

echo "[PLOTALL] ...plot statistics, stations used per cell..."
./gmtstatsplot.sh -jpg -stats strain_stats.dat --stats-stations -psta -l -leg -o ${prefix}output_stats-stations${suffix}
echo "[PLOTALL] ...plot statistics, optimal smoothing distance..."
./gmtstatsplot.sh -jpg -stats strain_stats.dat --stats-doptimal -psta -l -leg -o ${prefix}output_stats-doptimal${suffix}
echo "[PLOTALL] ...plot statistics, sigma value..."
./gmtstatsplot.sh -jpg -stats strain_stats.dat --stats-sigma -psta -l -leg -o ${prefix}output_stats-sigma${suffix}

echo "--------------------------------------------------------------------------"
echo "[PLOTALL] output file produced from this script:"
echo "--------------------------------------------------------------------------"
echo "Horizontal velocities ..................... ${prefix}output_vel${suffix}.jpg"
echo "Principal axes of strain rates ............ ${prefix}output_str${suffix}.jpg"
echo "Rotational rates .......................... ${prefix}output_rot${suffix}.jpg"
echo "Dextral/sinistral maximum shear strain .... ${prefix}output_gtotaxes${suffix}.jpg"
echo "Maximum shear strain ...................... ${prefix}output_gtot${suffix}.jpg"
echo "Dilatation ................................ ${prefix}output_dil${suffix}.jpg"
echo "Second Invariant .......................... ${prefix}output_2inv${suffix}.jpg"
echo "Statistics, stations used per cell ........ ${prefix}output_stats-stations${suffix}.jpg"
echo "Statistics, optimal smoothing distance .... ${prefix}output_stats-doptimal${suffix}.jpg"
echo "Statistics, sigma value : ................. ${prefix}output_stats-sigma${suffix}.jpg"
echo "--------------------------------------------------------------------------"

# Print exit status
echo "[STATUS] Finished. Exit status: $?"

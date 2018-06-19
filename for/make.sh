gfortran -c voronoi_area_version.f90
gfortran -g -fbacktrace -ffpe-trap=zero,overflow,underflo visr.f -L. voronoi_area_version.o
cat CNRS_midas.vel | awk '{printf "%4s%4s%10.4f%10.4f%7.2f%5.2f%7.2f%5.2f%7.3f%s\n", $1,"_GPS",$2,$3,$4,$6,$5,$7,$8+0.001,"     4   7.8  1994.8"}' | grep -v -i MATG > midas.vel

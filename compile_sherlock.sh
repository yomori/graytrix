#!/bin/bash
#source /home/groups/kipac/yomori/setup_gcc.sh 
DIR_HEALPIX=/home/groups/kipac/yomori/packages/Healpix_3.11_gcc/


#gfortran test.f -o test.x -I${DIR_HEALPIX}/include -DGFORTRAN -fno-second-underscore -fopenmp -fPIC -L${DIR_HEALPIX}/lib -lhealpix -lcfitsio
#gfortran graytrix_v2.2_lite.f -o graytrix_lite_sherlock.x -I${DIR_HEALPIX}/include -DGFORTRAN -fno-second-underscore -fopenmp -fPIC -L${DIR_HEALPIX}/lib -lhealpix -lcfitsio
#gfortran heal2psi_upsample.f90 -o heal2psi_upsample_sherlock.x -I${DIR_HEALPIX}/include -DGFORTRAN -fno-second-underscore -fopenmp -fPIC -L${DIR_HEALPIX}/lib -lhealpix -lcfitsio
gfortran degrade_nsideout.f90 -o degrade_nsideout.x -I${DIR_HEALPIX}/include -DGFORTRAN -fno-second-underscore -fopenmp -fPIC -L${DIR_HEALPIX}/lib -lhealpix -lcfitsio


#ifort heal2psi.f90 -o heal2psi.x -I${DIR_HEALPIX}/include -L${DIR_HEALPIX}/lib -lhealpix -fopenmp -lcfitsio 
#ifort heal2psi_nside16384.f90 -o heal2psi_nside16384.x -I${DIR_HEALPIX}/include -DGFORTRAN -fno-second-underscore -fopenmp -fPIC -L${DIR_HEALPIX}/lib -lhealpix -lcfitsio
#ifort check_psimaps.f90 -o check_psimaps.x -I${DIR_HEALPIX}/include -L${DIR_HEALPIX}/lib -lhealpix -fopenmp -lcfitsio 
#ifort dat2fits.f -o dat2fits.x -I/project2/chihway/packages/Healpix_3.31-intel/include -L/project2/chihway/packages/Healpix_3.31-intel/lib -lhealpix -qopenmp -lcfitsio -lcfitsio -qopt-report=0

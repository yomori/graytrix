#!/bin/bash

DIR_HEALPIX=/project2/chihway/packages/Healpix_3.31-intel/

ifort heal2psi.f90 -o heal2psi.x -I${DIR_HEALPIX}/include -L${DIR_HEALPIX}/lib -lhealpix -qopenmp -lcfitsio -lcfitsio -qopt-report=0

ifort graytrix_v2.2_lite.f -o graytrix_lite.x -I${DIR_HEALPIX}/include -L${DIR_HEALPIX}/lib -lhealpix -qopenmp -lcfitsio -lcfitsio -qopt-report=0

ifort check_psimaps.f90 -o check_psimaps.x -I${DIR_HEALPIX}/include -L${DIR_HEALPIX}/lib -lhealpix -qopenmp -lcfitsio -lcfitsio -qopt-report=0

#ifort dat2fits.f -o dat2fits.x -I/project2/chihway/packages/Healpix_3.31-intel/include -L/project2/chihway/packages/Healpix_3.31-intel/lib -lhealpix -qopenmp -lcfitsio -lcfitsio -qopt-report=0

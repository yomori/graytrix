#!bin/bash

for i in {0..19}
do
rm testmag.${i}
./zzz.x -infile mdpl2_psi2.${i} -outfile testmag.${i}
done

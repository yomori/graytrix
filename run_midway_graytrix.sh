#!/bin/bash
for i in {1..1}
do
echo '#!/bin/bash' > submit_job
echo "#SBATCH -t 8:00:00">> submit_job
echo "#SBATCH --partition=bigmem2" >> submit_job
echo "#SBATCH --nodes=1" >> submit_job
echo '#SBATCH -c 20' >> submit_job
echo "#SBATCH --job-name=graytixmini" >> submit_job
echo "#SBATCH --mem=250GB">> submit_job
echo "cd $PWD" >> submit_job
echo "export OMP_NUM_THREADS=20" >> submit_job
echo "source /project2/chihway/setup/setup_intel.sh" >> submit_job
echo  "./graytrix_lite.x -nres 13 -d /scratch/midway2/yomori/density_maps/ -om 0.307115 -th 25.0 -nt 20 -f raytrace8192test < infile_jul26" >> submit_job
# The lite version only outputs kappa/shear/omega not the source positions
sbatch submit_job
done

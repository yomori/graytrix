#!/bin/bash

echo '#!/bin/bash' > submit_job
echo "#SBATCH -t 48:00:00">> submit_job
echo "#SBATCH --partition=chihway" >> submit_job
echo "#SBATCH --exclusive" >> submit_job
echo "#SBATCH --nodes=1" >> submit_job
echo "#SBATCH --job-name=heal2psi" >> submit_job
echo "#SBATCH --mem=76GB">> submit_job
echo "cd $PWD" >> submit_job
echo "export OMP_NUM_THREADS=40" >> submit_job
echo "source /project2/chihway/setup/setup_intel.sh" >> submit_job
for i in {0..3}
do
echo  "./heal2psi.x -nside 8192 -lmax 24576 -infile /scratch/midway2/yomori/density_maps/mdpl2_graytrix_sdens_${i}.fits -outfile /scratch/midway2/yomori/density_maps/mdpl2_psi2.${i}" >> submit_job
done
sbatch submit_job

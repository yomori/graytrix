#!/bin/bash
dir_in=/scratch/midway2/yomori/density_maps/
dir_out=/scratch/midway2/yomori/density_maps/
cambpar=mdpl2_params.ini
echo '#!/bin/bash' > submit_job
echo "#SBATCH -t 24:00:00">> submit_job
echo "#SBATCH --partition=chihway" >> submit_job
echo "#SBATCH --exclusive" >> submit_job
echo "#SBATCH --nodes=1" >> submit_job
echo "#SBATCH --job-name=heal2psi${i}" >> submit_job
#echo "#SBATCH --array=5-6" >> submit_job
echo "#SBATCH --mem=42GB">> submit_job
#echo "#SBATCH --begin=04:00:00" >> submit_job
echo "cd $PWD" >> submit_job
echo "export OMP_NUM_THREADS=40" >> submit_job
echo "source /project2/chihway/setup/setup_intel.sh" >> submit_job
echo  python make_kappaslices.py ${dir_in} ${dir_out} ${cambpar} >> submit_job
sbatch submit_job

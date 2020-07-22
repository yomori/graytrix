#!/bin/bash
dir_in=/home/users/yomori/oak/mdpl2/dens/
dir_out=/scratch/users/yomori/raytrace2/
cambpar=mdpl2_params.ini
echo '#!/bin/bash' > submit_job
echo "#SBATCH -t 24:00:00">> submit_job
echo "#SBATCH --partition=kipac" >> submit_job
echo "#SBATCH -c 4" >> submit_job
echo "#SBATCH --nodes=1" >> submit_job
echo "#SBATCH --job-name=heal2psi${i}" >> submit_job
echo "#SBATCH --mem=60GB">> submit_job
echo "cd $PWD" >> submit_job
echo "export OMP_NUM_THREADS=4" >> submit_job
echo "source /home/groups/kipac/yomori/setup.sh" >> submit_job
echo  make_densityfiles.py ${dir_in} ${dir_out} ${cambpar} >> submit_job
sbatch submit_job

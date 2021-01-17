#!/bin/bash

echo '#!/bin/bash' > submit_job
echo "#SBATCH -t 12:00:00">> submit_job
echo "#SBATCH --partition=kipac,iric,hns" >> submit_job
echo "#SBATCH --nodes=1" >> submit_job
echo '#SBATCH --cpus-per-task=12' >> submit_job
echo "#SBATCH --job-name=graytix" >> submit_job
echo "#SBATCH --mem=180000">> submit_job
echo "#SBATCH --array=1-1" >> submit_job
echo "cd $PWD" >> submit_job
echo "export OMP_NUM_THREADS=12" >> submit_job
#echo "source activate cobaya" >> submit_job
#for i in {1..105}
#do
#echo  ./degrade_nsideout.x -infile /oak/stanford/orgs/kipac/users/yomori/mdpl2/outputs/raytrace_v4/raytrace_cori/raytrace16384_ip20.zs'$SLURM_ARRAY_TASK_ID' >> submit_job
echo  ./degrade_nsideout.x -infile /oak/stanford/orgs/kipac/users/yomori/mdpl2/outputs/raytrace_v4/raytrace_cori/raytrace16384_ip20_cmbkappa.zs'$SLURM_ARRAY_TASK_ID' >> submit_job
#echo  ./dat2fits.x -d /oak/stanford/orgs/kipac/users/yomori/mdpl2/outputs/raytrace_v4/raytrace_cori/raytrace16384_ip20.zs'$SLURM_ARRAY_TASK_ID'.mag.dat -f /oak/stanford/orgs/kipac/users/yomori/mdpl2/outputs/raytrace_v4/raytrace_cori/raytrace16384_ip20.zs'$SLURM_ARRAY_TASK_ID'.fits>> submit_job
#done
sbatch submit_job

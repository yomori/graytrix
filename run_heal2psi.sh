#!/bin/bash
for i in {145,136,133}
do
echo '#!/bin/bash' > submit_job
echo "#SBATCH -t 2:00:00">> submit_job
echo "#SBATCH --partition=broadwl" >> submit_job
echo "#SBATCH --exclusive" >> submit_job
echo "#SBATCH --nodes=1" >> submit_job
echo "#SBATCH --job-name=heal2psi${i}" >> submit_job
#echo "#SBATCH --array=5-6" >> submit_job
echo "#SBATCH --mem=42GB">> submit_job
#echo "#SBATCH --begin=04:00:00" >> submit_job
echo "cd $PWD" >> submit_job
echo "export OMP_NUM_THREADS=24" >> submit_job
echo "source /project2/chihway/setup/setup_intel.sh" >> submit_job
#echo python project_particles.py ${i} 25 '$SLURM_ARRAY_TASK_ID'>> submit_job
#echo  ./heal2psi.x -nside 8192 -lmax 12000 -infile /scratch/users/yomori/kappa/mdpl2_cmbkappa_'$SLURM_ARRAY_TASK_ID'.fits -outfile /scratch/users/yomori/psifiles/mdpl2_psi.'$SLURM_ARRAY_TASK_ID' >> submit_job
echo  "./heal2psi.x -nside 8192 -lmax 16384 -infile /project2/chihway/yuuki/density/mdpl2_graytrix_sdens_${i}.fits -outfile mdpl2_psi2.${i}" >> submit_job
sbatch submit_job
done

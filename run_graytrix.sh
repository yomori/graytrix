#!/bin/bash
for i in {1..1}
do
echo '#!/bin/bash' > submit_job
echo "#SBATCH -t 4:00:00">> submit_job
echo "#SBATCH --partition=bigmem2" >> submit_job
echo "#SBATCH --nodes=1" >> submit_job
echo '#SBATCH -c 24' >> submit_job
echo "#SBATCH --job-name=graytixmini" >> submit_job
echo "#SBATCH --mem=250GB">> submit_job
echo "cd $PWD" >> submit_job
echo "export OMP_NUM_THREADS=20" >> submit_job
echo "source /project2/chihway/setup/setup_intel.sh" >> submit_job
echo  "./graytrix.x -nres 13 -d /project2/chihway/yuuki/density/ -om 0.307115 -th 25.0 -nt 20 -f raytrace8192 < infile" >> submit_job
sbatch submit_job
done

for i in  0.1
do

cd $i

cat << Eof > sbatch_HEOM.sh
#!/bin/bash
#SBATCH -p action
#SBATCH -o output.log
#SBATCH --mem=60GB
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --job-name=sec
#SBATCH --open-mode=append

time /scratch/smontill/Soderquist_project/test/bin/rhot ./input.json

Eof

chmod +x sbatch_HEOM.sh
sbatch sbatch_HEOM.sh
cd ../

done

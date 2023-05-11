#!/bin/sh -login


#SBATCH --time=03:59:59             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1                   # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks-per-node=1         # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=1           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=25G                   # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name diamond          # you can give your job a name for easier identification (same as -J)
#SBATCH --output=job_reports/%x_%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/anaconda3"

#Change to current directory
cd ${PBS_O_WORKDIR}

#Export paths to conda
export PATH="${conda}/envs/diamond/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/diamond/lib:${LD_LIBRARY_PATH}"

# Arabidopsis
# diamond blastx -d Araport11 \
# -q ../ref/cds_from_genomic.mod.fna \
# -o Vriparia-Athaliana.tsv \
# --iterate \
# --max-target-seqs 1 \
# --unal 0

# Vvinifera
cd Vvinifera

diamond blastx -d Vvinifera_VCost3 \
-q ../cds_from_genomic.mod.fna \
-o Vriparia-Vvinifera.tsv \
--iterate \
--max-target-seqs 1 \
--unal 0

echo "Done"
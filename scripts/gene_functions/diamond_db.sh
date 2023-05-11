#!/bin/sh -login


#SBATCH --time=03:59:59             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1                   # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks-per-node=1         # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=1           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=15G                   # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name diamond_db          # you can give your job a name for easier identification (same as -J)
#SBATCH --output=%x_%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/anaconda3"

#Change to current directory
cd ${PBS_O_WORKDIR}

#Export paths to conda
export PATH="${conda}/envs/diamond/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/diamond/lib:${LD_LIBRARY_PATH}"

# Arabidopsis
# mkdir Arabidopsis
# cd Arabidopis
# diamond makedb --in Araport11_genes.201606.pep.fasta -d Araport11

# cd ../

# Vitis vinifera
mkdir Vvinifera 
cd Vvinifera

wget https://urgi.versailles.inra.fr/files/Vini/Vitis%2012X.2%20annotations/vitviv2.pep.fasta.zip
unzip -a vitviv2.pep.fasta.zip
diamond makedb --in vitviv2.pep.fasta -d Vvinifera_VCost3
...

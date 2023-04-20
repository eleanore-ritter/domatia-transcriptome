#!/bin/sh -login


#SBATCH --time=100:59:59             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=4                   # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks-per-node=1         # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=1           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=60G                   # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name star-index             # you can give your job a name for easier identification (same as -J)
#SBATCH --output=job_reports/%x_%j.SLURMout


#change to working directory
cd $PBS_O_WORKDIR

#Create index files
echo "Making index files"
$HOME/programs/STAR/bin/Linux_x86_64/STAR \
--runThreadN 4 \
--runMode genomeGenerate \
--genomeDir ../ref/star-index \
--genomeFastaFiles ../ref/assembly.fasta \
--sjdbGTFfile ../ref/genomic.gtf \
--genomeSAindexNbases 13

echo "Done"

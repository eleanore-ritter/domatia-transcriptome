#!/bin/sh -login


#SBATCH --time=03:59:59             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1                   # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks-per-node=1         # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=4           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=15G                   # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name star             # you can give your job a name for easier identification (same as -J)
#SBATCH --output=job_reports/%x_%j.SLURMout

#Variables should be set automatically
sample=$(pwd | sed s/^.*\\///)
fastq1=$(ls *.trimmed.1.fastq.gz)
fastq2=$(ls *.trimmed.2.fastq.gz)

#change to working directory
cd $PBS_O_WORKDIR

#Make star directory
mkdir star

#Run star
echo "Running star on ${sample}"
$HOME/programs/STAR/bin/Linux_x86_64/STAR \
--runThreadN 4 \
--genomeDir ../../ref/star-index/ \
--readFilesCommand gunzip -c --readFilesIn ${fastq1} ${fastq2} \
--outFileNamePrefix star/${sample} \
--quantMode GeneCounts

echo "Done"

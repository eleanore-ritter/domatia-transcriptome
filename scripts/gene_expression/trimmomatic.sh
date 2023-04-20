#!/bin/bash --login
#SBATCH --time=03:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=200GB
#SBATCH --job-name trimmomatic
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/anaconda3"

#Set variables
threads="4"
sample=$(pwd | sed 's/.*\///')
r1=$(ls *R1_001.fastq.gz)
r2=$(ls *R2_001.fastq.gz)
t1="${sample}.trimmed.1.fastq.gz"
t2="${sample}.trimmed.2.fastq.gz"
t3="${sample}.trimmed.1.single.fastq.gz"
t4="${sample}.trimmed.2.single.fastq.gz"
adapters="NexteraPE-PE"
path1="fastqc"

#Path to trimmomatic fastas
adapter_path="${conda}/envs/polishing/share/trimmomatic/adapters"

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/polishing/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/polishing/lib:${LD_LIBRARY_PATH}"

#Set temporary directories for large memory operations
export TMPDIR=$(pwd | sed s/data.*/data/)
export TMP=$(pwd | sed s/data.*/data/)
export TEMP=$(pwd | sed s/data.*/data/)

echo "Running trimmomatic on ${r1} and ${r2}"

trimmomatic PE \
-threads ${threads} \
-phred33 \
-trimlog job_reports/trim_log.txt \
-summary job_reports/trim_summary.txt \
${r1} ${r2} ${t1} ${t3} ${t2} ${t4} \
ILLUMINACLIP:${adapter_path}/${adapters}.fa:2:30:10

echo "Running FastQC on trimmed reads"
mkdir ${path1}

fastqc \
-t ${threads} \
-o ${path1} ${t1}

fastqc \
-t ${threads} \
-o ${path1} ${t2}

echo "Done"

#!/bin/bash
#
#SBATCH --job-name=blastx
#SBATCH --mem=20G
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --job-name blastx
#SBATCH --output=job_reports/%x_%j.SLURMout


#change to working directory
cd $PBS_O_WORKDIR

module purge
module load BLAST

blastx \
-query rna.fna \
-db nr \
-evalue 1e-3 \
-max_target_seqs 1 \
-word_size 6 \
-taxids eukaryotes \
-outfmt 5 \
-out Vriparia

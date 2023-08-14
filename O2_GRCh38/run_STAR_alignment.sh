#!/bin/bash
#SBATCH -c 8                    	# Request num cores
#SBATCH -t 0-0:30            		# Runtime in D-HH:MM format
#SBATCH -p short             		# Partition to run in
#SBATCH --mem=32G                	# Memory total in GB (for all cores)
#SBATCH -o Run_STAR_%A_%a.out	        # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e Run_STAR_%A_%a.err   	# File to which STDERR will be written, including job ID (%j)
#SBATCH --mail-user=StuartAi@BroadInstitute.org
#SBATCH --mail-type=ALL
#SBATCH --array=1-48

#####################################################################################################
# User params:

WORKING_DIR="/n/data1/dfci/pedonc/kadoch/Aidan/20230810_Biomek_RNAseq_QC/full_analysis"
IN_DIR="${WORKING_DIR}/fastq"
OUT_DIR="${WORKING_DIR}/aligned"
SAMPLE_SHEET="${WORKING_DIR}/scripts/samples.txt"
GENOME_REF="/n/data1/dfci/pedonc/kadoch/Aidan/REF/STAR/GRCh38_primary_ercc_75"
GTF="/n/data1/dfci/pedonc/kadoch/Aidan/REF/GRCh38_primary_ercc.gtf"
####################################################################################################
# Set working directory
cd ${WORKING_DIR}

module load star/2.7.9a

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_SHEET")

STAR \
	--quantMode GeneCounts \
	--sjdbGTFfile ${GTF} \
	--runThreadN 8 \
	--outFileNamePrefix ${OUT_DIR}/${SAMPLE} \
	--outSAMtype BAM SortedByCoordinate \
	--genomeDir ${GENOME_REF} \
	--readFilesIn ${IN_DIR}/${SAMPLE}_R1.fastq ${IN_DIR}/${SAMPLE}_R2.fastq

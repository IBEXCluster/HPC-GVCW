#!/bin/bash
#######################################################################################################################################################
# SAMtools and GATK HaplotypeCaller workflow 
#  
# Tools used:
# -----------
#    1. Merge Bam files from different Flow cells	- Samtools Merge
#    2. Sort the merged Bam files			- Samtools sort
#    3. Samtools index					- Samtools index
#    4. Haplotype Caller				- GATK HaplotypeCaller
#
#
#     Version 1.0 dated 2 Apr 2020
#     Modification required for any users:
#      1. INPUT (directory location for fastq.gz files)
#      2. PROJECT (directory location for generating the output files)
#
########################################################################################################################################################

## Software Modules 
module load samtools/1.8 tabix/0.2.6
export PATH=/ibex/scratch/projects/c2028/3k_Project/naga/software/gatk-4.1.6.0:$PATH

## Directory Variables 
export REF=/ibex/scratch/projects/c2028/3k_Project/REF/Nipponbare_chr.fasta
export INPUT=/ibex/scratch/projects/c2028/3k_Project/naga/OUTPUT_GATK4.1.6/tmpBAM
export PROJECT=/ibex/scratch/projects/c2028/3k_Project/naga/OUTPUT_GATK4.1.6

export BAM=${PROJECT}/BAM;
export VCF=${PROJECT}/VCF;
export LOGS=${PROJECT}/LOGS;
mkdir -p $BAM;
mkdir -p $VCF;
mkdir -p $LOGS;

export QoS="-A ibex-cs"
## Workflow steps 
for DIR in `ls -d $INPUT/*`;
do 
  PREFIX=${DIR#$INPUT/}

  #### Step 1. Merge the BAM files using Samtools 
    MEM="32gb"
    CORES=16
    JOB1_NAME="samtools-merge"
    JOB1_TYPE="sbatch $QoS --partition=batch --job-name=${JOB1_NAME}.${PREFIX} --time=5:00:00 --output=$LOGS/${JOB1_NAME}.${PREFIX}.%J.out --error=$LOGS/${JOB1_NAME}.${PREFIX}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM}" ;
    JOB1_CMD="time -p samtools merge -f -n --threads $CORES $BAM/${PREFIX}.bam $DIR/${PREFIX}_*.sorted.fxmt.mkdup.addrep.bam"
    JOB1_ID=$(${JOB1_TYPE} --parsable --wrap="${JOB1_CMD}");
    echo "$PREFIX sample with the job id=$JOB1_ID and Job Name=$JOB1_NAME submitted"  

  #### Step 2. Sort the BAM files using Samtools
    MEM="32gb"
    CORES=16
    JOB2_NAME="sam-sort"
    JOB2_TYPE="sbatch $QoS --partition=batch --job-name=${JOB2_NAME}.${PREFIX} --time=10:00:00 --output=$LOGS/${JOB2_NAME}.${PREFIX}.%J.out --error=$LOGS/${JOB2_NAME}.${PREFIX}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM}" ;
    JOB2_CMD="time -p samtools sort --threads $CORES -o $BAM/${PREFIX}.sorted.bam $BAM/${PREFIX}.bam"
    JOB2_ID=$(${JOB2_TYPE} --parsable --dependency=afterok:${JOB1_ID} --wrap="${JOB2_CMD}");
    echo "$PREFIX sample with the job id=$JOB2_ID and Job Name=$JOB2_NAME submitted" 
 
  #### Step 3. Samtools Index
    MEM="16gb"
    CORES=1
    JOB3_NAME="sam-index"
    JOB3_TYPE="sbatch $QoS --partition=batch --job-name=${JOB3_NAME}.${PREFIX} --time=4:00:00 --output=$LOGS/${JOB3_NAME}.${PREFIX}.%J.out --error=$LOGS/${JOB3_NAME}.${PREFIX}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM}" ;
    JOB3_CMD="time -p samtools index $BAM/${PREFIX}.sorted.bam"
    JOB3_ID=$(${JOB3_TYPE} --parsable --dependency=afterok:${JOB2_ID} --wrap="${JOB3_CMD}");
    echo "$PREFIX sample with the job id=$JOB3_ID and Job Name=$JOB3_NAME submitted" 

  #### 4. HaplotypeCaller 
    MEM="64gb"
    CORES=16
    JOB4_NAME="HC"
    JOB4_TYPE="sbatch $QoS --partition=batch --job-name=${JOB4_NAME}.${PREFIX} --time=22:00:00 --output=$LOGS/${JOB4_NAME}.${PREFIX}.%J.out --error=$LOGS/${JOB4_NAME}.${PREFIX}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM}" ;
    JOB4_CMD="time -p gatk HaplotypeCaller --input $BAM/${PREFIX}.sorted.bam --output $VCF/${PREFIX}.snps.indels.g.vcf --reference $REF --emit-ref-confidence GVCF --min-base-quality-score 20 --native-pair-hmm-threads ${CORES} --bam-output $BAM/${PREFIX}.assembledHC.bam"
    JOB4_ID=$(${JOB4_TYPE} --parsable --dependency=afterok:${JOB3_ID} --wrap="${JOB4_CMD}");
    echo "$PREFIX sample with the job id=$JOB4_ID and Job Name=$JOB4_NAME submitted" 

##### 5. Compress the g.VCF file using bgzip
    MEM="32gb"
    CORES=1
    JOB5_NAME="bgzip"
    JOB5_TYPE="sbatch --partition=batch --job-name=${JOB5_NAME}.${PREFIX} --time=7:00:00 --output=$LOGS/${JOB5_NAME}.${PREFIX}.%J.out --error=$LOGS/${JOB5_NAME}.${PREFIX}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM}" ;
    JOB5_CMD="time -p bgzip -o $VCF/$PREFIX.snps.indels.g.vcf" ;
    JOB5_ID=$(${JOB5_TYPE} --parsable --dependency=afterok:${JOB4_ID} --wrap="${JOB5_CMD}");
    echo "$PREFIX sample with the job id=$JOB5_ID and Job Name=$JOB5_NAME submitted"

 ##### 6. Create Tabix-Index for the g.VCF.gz file
    MEM="32gb"
    CORES=1
    JOB6_NAME="tabix"
    JOB6_TYPE="sbatch --partition=batch --job-name=${JOB6_NAME}.${PREFIX} --time=30:00 --output=$LOGS/${JOB6_NAME}.${PREFIX}.%J.out --error=$LOGS/${JOB6_NAME}.${PREFIX}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM}" ;
    JOB6_CMD="time -p tabix -o $VCF/$PREFIX.snps.indels.g.vcf.gz" ;
    JOB6_ID=$(${JOB6_TYPE} --parsable --dependency=afterok:${JOB5_ID} --wrap="${JOB6_CMD}");
    echo "$PREFIX sample with the job id=$JOB6_ID and Job Name=$JOB6_NAME submitted"   
done

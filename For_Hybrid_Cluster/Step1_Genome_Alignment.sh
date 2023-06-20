#!/bin/bash
#######################################################################################################################################################
# BWA, SAMtools and GATK (Picard functionality) for Genome Alignment workflow 
#  
# Tools used:
# -----------
#  1. Sorted Genome Mapping 		- BWA MEM with Samtools
#  2. Fix mate information  		- gatk FixMateInformation
#  3. Mark duplicates			- gatk MarkDuplicates
#  4. Add or replace read groups	- gatk AddOrReplaceReadGroups
#
#     Version 1.0 dated 31 Mar 2020
#     Modification required for any users:
#      1. INPUT (directory location for fastq.gz files)
#      2. PROJECT (directory location for generating the output files)
#
########################################################################################################################################################

## Software Modules 
module load bwa/0.7.17/gnu-6.4.0 samtools/1.8 

export PATH=/ibex/scratch/projects/c2028/3k_Project/naga/software/gatk-4.1.6.0:$PATH

## Directory Variables 
export REF=/ibex/scratch/projects/c2028/3k_Project/REF/Nipponbare_chr.fasta
export INPUT=/ibex/scratch/projects/c2028/3k_Project/naga/INPUT1
export PROJECT=/ibex/scratch/projects/c2028/3k_Project/naga/OUTPUT_GATK4.1.6
export LOGS=${PROJECT}/LOGS;
mkdir -p $LOGS;


#export QoS="-A ibex-cs"
## Workflow steps 
for DATA in `ls $INPUT/*_1.fastq.gz`;
do 
  PREFIX=`basename $DATA _1.fastq.gz` ;
  LOCATION=${DATA%/*};
  SAMPLE=${PREFIX%*_*};
  export BAM=${PROJECT}/tmpBAM/$SAMPLE
  mkdir -p $BAM;
  #echo $PREFIX $LOCATION $SAMPLE

  #### Step 1. Genome Mapping: BWA MEM
    MEM="32gb"
    CORES=8
    JOB1_NAME="bwa-mem"
    JOB1_TYPE="sbatch $QoS --partition=batch --job-name=${JOB1_NAME}.${PREFIX} --time=2:00:00 --output=$LOGS/${JOB1_NAME}.${PREFIX}.%J.out --error=$LOGS/${JOB1_NAME}.${PREFIX}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM}" ;
    JOB1_CMD="time -p bwa mem -M -t $CORES $REF $LOCATION/${PREFIX}_1.fastq.gz $LOCATION/${PREFIX}_2.fastq.gz | samtools view -@ $CORES -b -S -h -q 30 - | samtools sort -T ${PREFIX} - > $BAM/$PREFIX.sorted.bam"
    echo "${JOB1_TYPE} --parsable --wrap="${JOB1_CMD}"" ;
    JOB1_ID=$(${JOB1_TYPE} --parsable --wrap="hostname; module list; which bwa; which samtools; ${JOB1_CMD}");
    echo "$PREFIX sample with the job id=$JOB1_ID and Job Name=$JOB1_NAME submitted"  

  #### Step 2. Fix mate information
    MEM="16gb"
    CORES=1
    JOB2_NAME="Fix-mate"
    JOB2_TYPE="sbatch $QoS --partition=batch --job-name=${JOB2_NAME}.${PREFIX} --time=7:00:00 --output=$LOGS/${JOB2_NAME}.${PREFIX}.%J.out --error=$LOGS/${JOB2_NAME}.${PREFIX}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM}" ;
    JOB2_CMD="time -p gatk FixMateInformation --INPUT $BAM/$PREFIX.sorted.bam --SORT_ORDER coordinate --OUTPUT $BAM/$PREFIX.sorted.fxmt.bam"
    echo "${JOB2_TYPE} --parsable --dependency=afterok:${JOB1_ID} --wrap="${JOB2_CMD}""
    JOB2_ID=$(${JOB2_TYPE} --parsable --dependency=afterok:${JOB1_ID} --wrap="hostname; module list; which gatk; ${JOB2_CMD}");
    echo "$PREFIX sample with the job id=$JOB2_ID and Job Name=$JOB2_NAME submitted" 
 
  #### Step 3. Mark duplicate reads
    MEM="16gb"
    CORES=1
    JOB3_NAME="Mark-dupe"
    JOB3_TYPE="sbatch $QoS --partition=batch --job-name=${JOB3_NAME}.${PREFIX} --time=8:00:00 --output=$LOGS/${JOB3_NAME}.${PREFIX}.%J.out --error=$LOGS/${JOB3_NAME}.${PREFIX}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM}" ;
    JOB3_CMD="time -p gatk MarkDuplicates --INPUT $BAM/$PREFIX.sorted.fxmt.bam --METRICS_FILE $BAM/$PREFIX.metrics --OUTPUT $BAM/$PREFIX.sorted.fxmt.mkdup.bam "
    echo "${JOB3_TYPE} --parsable --dependency=afterok:${JOB2_ID} --wrap="${JOB3_CMD}""
    JOB3_ID=$(${JOB3_TYPE} --parsable --dependency=afterok:${JOB2_ID} --wrap="hostname; module list; which gatk; ${JOB3_CMD}");
    echo "$PREFIX sample with the job id=$JOB3_ID and Job Name=$JOB3_NAME submitted" 

  #### 4. Add or replace read groups
    MEM="16gb"
    CORES=1
    JOB4_NAME="Read-gp"
    JOB4_TYPE="sbatch $QoS --partition=batch --job-name=${JOB4_NAME}.${PREFIX} --time=12:00:00 --output=$LOGS/${JOB4_NAME}.${PREFIX}.%J.out --error=$LOGS/${JOB4_NAME}.${PREFIX}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM}" ;
    JOB4_CMD="time -p gatk AddOrReplaceReadGroups --INPUT $BAM/$PREFIX.sorted.fxmt.mkdup.bam --OUTPUT $BAM/$PREFIX.sorted.fxmt.mkdup.addrep.bam --RGID $SAMPLE --RGPL Illumina --RGSM $SAMPLE --RGLB $SAMPLE --RGPU unit1 --RGCN BGI --SORT_ORDER coordinate" 
     echo "${JOB4_TYPE} --parsable --dependency=afterok:${JOB3_ID} --wrap="${JOB4_CMD}""
    JOB4_ID=$(${JOB4_TYPE} --parsable --dependency=afterok:${JOB3_ID} --wrap="hostname; module list; which gatk; ${JOB4_CMD}");
    echo "$PREFIX sample with the job id=$JOB4_ID and Job Name=$JOB4_NAME submitted" 
done

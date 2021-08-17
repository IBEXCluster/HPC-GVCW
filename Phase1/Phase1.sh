#!/bin/bash
### Environment variables 
export REF=/project/k01/kathirn/3k/ref/Nipponbare_chr.fasta ;
export INPUT=/project/k1514/batch3/cleandata
export PROJECT=/project/k01/kathirn/3k/OUTPUT/batch2 ;
export CORE=1 ;
val=$1;
if [ $val -eq 0 ]
then 
  LINE=1;
else
 LINE=$((val + 1))
fi

## Read one Sample Per Rank 
DATA=`sed -n ${LINE}p  Phase1.txt`;
PREFIX=`basename $DATA _1.fastq.gz` ;
SAMPLE=${PREFIX%*_*};
BAM=${PROJECT}/tmpBAM/$SAMPLE

## Random wait before creating $BAM directory (to safeguard /luster file system)
export MAXWAIT=5
sleep $[ ( $RANDOM % $MAXWAIT )  + 1 ]s
mkdir -p $BAM;


## Step 1
echo "------------------------------------- Step 1 executing: BWA MEM ---------------------------------------------------------------------------------------"
echo "time -p bwa mem -M -t $CORE $REF $INPUT/${PREFIX}_1.fastq.gz $INPUT/${PREFIX}_2.fastq.gz | samtools view -@ $CORE -b -S -h -q 30 - | samtools sort -T /scratch/$USER/tmp/$SLURM_JOB_ID/${PREFIX} - > $BAM/$PREFIX.sorted.bam"
time -p bwa mem -M -t $CORE $REF $INPUT/${PREFIX}_1.fastq.gz $INPUT/${PREFIX}_2.fastq.gz | samtools view -@ $CORE -b -S -h -q 30 - | samtools sort -T /scratch/$USER/tmp/$SLURM_JOB_ID/${PREFIX} - > $BAM/$PREFIX.sorted.bam 
## Step 2
echo "------------------------------------- Step 2 executing: FixMateInformation  -----------------------------------------------------------------------------"
echo "time -p gatk --java-options -Xmx2g -Xms2g FixMateInformation --INPUT $BAM/$PREFIX.sorted.bam --SORT_ORDER coordinate --OUTPUT $BAM/$PREFIX.sorted.fxmt.bam --TMP_DIR /scratch/$USER/tmp/$SLURM_JOB_ID"
time -p gatk --java-options "-Xmx2g -Xms2g" FixMateInformation --INPUT $BAM/$PREFIX.sorted.bam --SORT_ORDER coordinate --OUTPUT $BAM/$PREFIX.sorted.fxmt.bam --TMP_DIR /scratch/$USER/tmp/$SLURM_JOB_ID
## Step 3
echo "------------------------------------- Step 3 executing: MarkDuplicates -----------------------------------------------------------------------------------"
echo "time -p gatk --java-options -Xmx2g -Xms2g MarkDuplicates --INPUT $BAM/$PREFIX.sorted.fxmt.bam --METRICS_FILE $BAM/$PREFIX.metrics --OUTPUT $BAM/$PREFIX.sorted.fxmt.mkdup.bam --TMP_DIR /scratch/$USER/tmp/$SLURM_JOB_ID"
time -p gatk --java-options "-Xmx2g -Xms2g" MarkDuplicates --INPUT $BAM/$PREFIX.sorted.fxmt.bam --METRICS_FILE $BAM/$PREFIX.metrics --OUTPUT $BAM/$PREFIX.sorted.fxmt.mkdup.bam --TMP_DIR /scratch/$USER/tmp/$SLURM_JOB_ID
## Step 4
echo "------------------------------------- Step 4 executing: AddOrReplaceReadGroups ----------------------------------------------------------------------------"
echo "time -p gatk --java-options -Xmx2g -Xms2g AddOrReplaceReadGroups --INPUT $BAM/$PREFIX.sorted.fxmt.mkdup.bam --OUTPUT $BAM/$PREFIX.sorted.fxmt.mkdup.addrep.bam --RGID $SAMPLE --RGPL Illumina --RGSM $SAMPLE --RGLB $SAMPLE --RGPU unit1 --RGCN BGI --SORT_ORDER coordinate --TMP_DIR /scratch/$USER/tmp/$SLURM_JOB_ID"
time -p gatk --java-options "-Xmx2g -Xms2g" AddOrReplaceReadGroups --INPUT $BAM/$PREFIX.sorted.fxmt.mkdup.bam --OUTPUT $BAM/$PREFIX.sorted.fxmt.mkdup.addrep.bam --RGID $SAMPLE --RGPL Illumina --RGSM $SAMPLE --RGLB $SAMPLE --RGPU unit1 --RGCN BGI --SORT_ORDER coordinate --TMP_DIR /scratch/$USER/tmp/$SLURM_JOB_ID

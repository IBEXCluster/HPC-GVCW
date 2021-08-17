#!/bin/bash
### Environment variables 
export REF=/project/k01/kathirn/3k/ref/Nipponbare_chr.fasta ;
export INPUT=/project/k01/kathirn/3k/OUTPUT/tmpBAM ;
export PROJECT=/project/k01/kathirn/3k/OUTPUT ;
export BAM=${PROJECT}/BAM;
export VCF=${PROJECT}/VCF;
export CORE=1 ;
val=$1;
if [ $val -eq 0 ]
then 
  LINE=1;
else
 LINE=$((val + 1))
fi

## Read one Sample Per Rank 
PREFIX=`sed -n ${LINE}p  Phase2.prefix.txt`;
DIRECTORY=`sed -n ${LINE}p  Phase2.directory.txt`;

echo "------------------------------------- Step 1: samtools-merge  ---------------------------------------------------------------------------------------"
time -p samtools merge -f -n $BAM/${PREFIX}.bam $DIRECTORY/${PREFIX}_*.sorted.fxmt.mkdup.addrep.bam
echo "------------------------------------- Step 2: sam-sort-----------------------------------------------------------------------------------------------"
time -p samtools sort -T /scratch/kathirn/tmp/$PREFIX -o $BAM/${PREFIX}.sorted.bam $BAM/${PREFIX}.bam
echo "------------------------------------- Step 3: sam-index ---------------------------------------------------------------------------------------------"
time -p samtools index $BAM/${PREFIX}.sorted.bam
echo "------------------------------------- Step 4: HaplotypeCaller ----------------------------------------------------------------------------------------"
time -p gatk --java-options "-Xmx2g -Xms2g" HaplotypeCaller --input $BAM/${PREFIX}.sorted.bam --output $VCF/${PREFIX}.snps.indels.g.vcf.gz --reference $REF --emit-ref-confidence GVCF --min-base-quality-score 20 --bam-output $BAM/${PREFIX}.assembledHC.bam --tmp-dir /scratch/$USER/tmp/$SLURM_JOB_ID
echo "------------------------------------- Step 5: Tabix --------------------------------------------------------------------------------------------------"  
time -p tabix -f $VCF/$PREFIX.snps.indels.g.vcf.gz

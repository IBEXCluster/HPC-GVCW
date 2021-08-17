#!/bin/bash
## Variables
# export REF=/project/k01/kathirn/3k/ref/Nipponbare_chr.fasta
#export OUTPUT=/project/k01/kathirn/3k/OUTPUT

export REF=/project/k1514/ref/genome6.fa
# export OUTPUT=/project/k1514/scripts/version2/Phase3/testOUTPUT
export OUTPUT=/project/k1514/scripts/version2/Phase3/OUTPUTall

## Rank based value 
val=$1;
if [ $val -eq 0 ]
then 
  LINE=1;
else
 LINE=$((val + 1))
fi

## Read one Split value Per Rank 
CHR=`sed -n ${LINE}p distribution.txt | awk '{print $1}'`
CHUNK=`sed -n ${LINE}p distribution.txt | awk '{print $2}'`

set INPUT_SNP
set INPUT_INDEL
for i in `seq 1 $CHUNK`; 
 do
  INPUT_SNP+="-I $OUTPUT/SNPs/filtered_snps.$CHR.$i.vcf.gz "; 
  INPUT_INDEL+="-I $OUTPUT/INDELs/$CHR.$i.vcf.gz "; 
done

 ## For SNPs
#  time -p gatk GatherVcfs ${INPUT_SNP} -O $OUTPUT/SNPs_per_Chr/$CHR.SNPs.vcf.gz -R $REF
#  time -p gatk VariantsToTable -R $REF -V $OUTPUT/SNPs_per_Chr/$CHR.SNPs.vcf.gz -F CHROM -F POS -F REF -F ALT -GF GT -O $OUTPUT/SNPs_per_Chr/$CHR.table
 ## For INDELs
  time -p gatk GatherVcfs ${INPUT_INDEL} -O $OUTPUT/INDELs_per_Chr/$CHR.INDELs.vcf.gz -R $REF
#  time -p gatk VariantsToTable -R $REF -V $OUTPUT/INDELs_per_Chr/$CHR.INDELs.vcf.gz -F CHROM -F POS -F REF -F ALT -GF GT -O $OUTPUT/INDELs_per_Chr/$CHR.table

#!/bin/bash
export PATH=/ibex/scratch/projects/c2028/3k_Project/naga/software/gatk-4.1.6.0:$PATH
export CHR=$1;
export CHUNK=$2;

export REF=/ibex/scratch/projects/c2028/3k_Project/REF/Nipponbare_chr.fasta
export LOCATION=/ibex/scratch/projects/c2028/3k_Project/naga/3k_combine_gVCF/50sample_chunks
export OUTPUT=/ibex/scratch/projects/c2028/3k_Project/naga/3k_combine_gVCF/50sample_chunks/merge

set INPUT_SNP
for i in `seq 1 $CHUNK`; 
 do
  INPUT_SNP+="-I $LOCATION/SNPs/filtered_snps.$CHR.$i.vcf "; 
done
#echo $INPUT_SNP
time -p gatk GatherVcfs ${INPUT_SNP} -O $OUTPUT/$CHR.SNPs.vcf.gz -R $REF
time -p gatk VariantsToTable -R $REF -V $OUTPUT/$CHR.SNPs.vcf.gz -F CHROM -F POS -F REF -F ALT -GF GT -O $OUTPUT/$CHR.table

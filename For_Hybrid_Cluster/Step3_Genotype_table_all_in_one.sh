#!/bin/bash
#########################################################################################################################################################
# Genotype_Table, version 1.0																#
#  Last update: Nov 1, 2020																#
# 	- Added optimal chunk size calculation    													#
#	- Optimal chunk size is based on any REFERENCE FILE and 											#
#		(a) Optimal number of chunk = {SLURM resource limit / (total no. of chromosomes * no. of job steps)}                                    #
#		(b) Optimal chunk size = {Max. chromosome size / optimal number of chunks}								#
# 	- Bigger Chromosome will be divided into more chunks                                                                                            #
#       - The chunk divident details are stored in the file "submitted_chunks_data.txt"  								#
#																			#
#########################################################################################################################################################
## Software Modules 
export PATH=/ibex/scratch/projects/c2028/3k_Project/naga/software/gatk-4.1.6.0:$PATH

REF=$1;
INPUT_DIR=$2;
PROJECT=$3;
if [ "$#" -eq 3 ];
then
  if [ ! -f $REF ]; then
    printf " Your Genome reference file does not exist \n" ;
    break;
  fi;
  if [ ! -d $INPUT_DIR ];then
    printf " Your INPUT [ $INPUT_DIR ] directory does not exist \n" ;
    break;
  fi
  if [ -d $PROJECT ];then
   echo "###################### WARNING! "
   echo " Please delete $PROJECT and RERUN the script"
  fi
else
 printf "\033c"
 echo " ***************************************************************************************************************************"
 echo ""
 echo " Run this script with 3 arguments:" 
 echo "      (a) Your Reference file"
 echo "      (b) Your gVCF files directory and "
 echo "      (c) Your output file directory "
 echo "  (absulute path required for all these options) "
 echo ""
 echo ""
 echo "  Example: "
 echo "    ./Step3_Genotype_table.sh /ibex/scratch/projects/c2028/3k_Project/REF/Nipponbare_chr.fasta /ibex/scratch/projects/c2028/3k_Project/01_MAGIC16_3KRGP/03_genome1/Stet1_Step2_OUTPUT /ibex/scratch/projects/c2028/3k_Project/naga/3k_combine_gVCF"
 echo ""
 echo ""
 echo " ***************************************************************************************************************************"
 exit;
fi

## Sample/Data Variables 
export gVCF=${PROJECT}/gVCF;
export logs=${PROJECT}/logs;
export SNPs=${PROJECT}/SNPs;
export INDELs=${PROJECT}/INDELs;
mkdir -p $gVCF;
mkdir -p $SNPs;
mkdir -p $logs;


## Calculate the Chromosome split 
MaxSize=0;
MaxJobs=1000;
JobSteps=4;

# Find the Max. size of Chromosome 
while IFS=$'\t' read -r -a myREF
do 
 ChrName=${myREF[0]};
 ChrLen=${myREF[1]};
 if [ $MaxSize -lt $ChrLen ]; then
  MaxSize=$ChrLen
 fi
done < $REF.fai

# Find the total numbers of Chr
TotalChr=`cat $REF.fai | wc -l `

# Find the best possible Chunk size 
Divide=1;
JobSteps=4;
AssinedJobs=$(( $Divide * $TotalChr * $JobSteps )) ;
#echo "$Divide and $AssinedJobs"
while [ ${MaxJobs} -gt $AssinedJobs ];
 do
   Divide=$(( $Divide + 1 ));
   AssinedJobs=$(( $Divide * $TotalChr * $JobSteps )) ;
 done
Divide=$(( Divide - 2 ));
Chunk=$(( $MaxSize / $Divide ));
 echo "		total Chr = $TotalChr, Max Chr length = $MaxSize, Chunk number = $Divide, Chunk size = $Chunk "


# Store the all the Chromosome Chunks for reference in a file "submitted_chunks_data.txt"
rm -Rf submitted_chunks_data.txt; 
while IFS=$'\t' read -r -a myREF
do
 ChrName=${myREF[0]};
 ChrLen=${myREF[1]};
 Part=$(( $ChrLen / $Chunk )) ;
 tmp=$(( $Part * $Chunk )) ;
 if [ $ChrLen -gt $tmp ]; then
    Part=$(( $Part + 1 ));
    echo "$ChrName split into $Part parts" >> submitted_chunks_data.txt
 else 
    echo "$ChrName split into $Part parts" >> submitted_chunks_data.txt
 fi
done < $REF.fai

###
QoS="-A ibex-cs"
####

## Add all the *.g.vcf files for Genotyping 
set INPUT_TMP
for i in `ls -l ${INPUT_DIR}/*.g.vcf.gz | awk '{print $9}'`
do
 INPUT_TMP+="$i -V "; 
done
INPUT=${INPUT_TMP::-4}
 #echo $INPUT
 MEM="115gb"
 CORES=4;
## READ ALL CHROMOSOME and LENGTH from REFERENCE INDEX file 
while IFS=$'\t' read -r -a myREF
do 
 size=1;
 DefSize=$Chunk;
 ChunkSize=$DefSize;
 ChrName=${myREF[0]};
 ChrLen=${myREF[1]};

 ## PREPARE CHUNKS for EACH CHROMOSOME 
  for ((Start=1; Start<$ChrLen; Start+=$DefSize))
   do 
     End=$(( $size*$ChunkSize ))
     if [ $End -lt $ChrLen ]
     then
   echo "		Debug: Chromosome=$ChrName   Chunk=chunk_$size   Start:$Start  End:$End"   
      ## CombineGVCFs (Chunk by Chunk)
       JOB1_NAME="$ChrName.chunk_$size.CombineGVCFs"
       JOB1_TYPE="sbatch --partition=batch $QoS --job-name=${JOB1_NAME} --time=10-00:00:00 --output=$logs/${JOB1_NAME}.%J.out --error=$logs/${JOB1_NAME}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
       JOB1_CMD="time -p gatk CombineGVCFs --variant $INPUT --reference $REF --intervals $ChrName:$Start-$End --output $gVCF/Combine.$ChrName.$size.vcf.gz";
       JOB1_ID=$(${JOB1_TYPE} --parsable --wrap="${JOB1_CMD}");
   echo "$JOB1_NAME with the job id=$JOB1_ID submitted";
      ## CALL GenotypeGVF (Chunk by Chunk)
       JOB2_NAME="$ChrName.chunk_$size.GenotypeGVCFs"
       JOB2_TYPE="sbatch --partition=batch $QoS --job-name=${JOB2_NAME} --time=10-00:00:00 --output=$logs/${JOB2_NAME}.%J.out --error=$logs/${JOB2_NAME}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
       JOB2_CMD="time -p gatk GenotypeGVCFs --variant $gVCF/Combine.$ChrName.$size.vcf.gz --reference $REF --output $gVCF/Genotype.$ChrName.$size.vcf.gz --intervals $ChrName:$Start-$End";
       JOB2_ID=$(${JOB2_TYPE} --parsable --dependency=afterok:${JOB1_ID} --wrap="${JOB2_CMD}");
   echo "$JOB2_NAME with the job id=$JOB2_ID submitted";
      ## Select VARIANT=SNPs (Chunk by Chunk)
       JOB3_NAME="$ChrName.chunk_$size.SNPs.SelectVariants"
       JOB3_TYPE="sbatch --partition=batch $QoS --job-name=${JOB3_NAME} --time=3-00:00:00 --output=$logs/${JOB3_NAME}.%J.out --error=$logs/${JOB3_NAME}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
       JOB3_CMD="time -p gatk SelectVariants --variant $gVCF/Genotype.$ChrName.$size.vcf.gz --reference $REF -select-type SNP --output $SNPs/$ChrName.$size.vcf.gz" ;
       JOB3_ID=$(${JOB3_TYPE} --parsable --dependency=afterok:${JOB2_ID} --wrap="${JOB3_CMD}");
   echo "$JOB3_NAME with the job id=$JOB3_ID submitted";
      ## Call SNPs Filteration (Chunk by Chunk)
       JOB4_NAME="$ChrName.chunk_$size.SNPs.VariantFiltration"
       JOB4_TYPE="sbatch --partition=batch $QoS --job-name=${JOB4_NAME} --time=3-00:00:00 --output=$logs/${JOB4_NAME}.%J.out --error=$logs/${JOB4_NAME}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
       JOB4_CMD="time -p gatk VariantFiltration --variant $SNPs/$ChrName.$size.vcf.gz --reference $REF --filter-expression \"QUAL < 30.0 || QD < 2.5 || MQ < 40.0 || MQRankSum < -12.5 || MQRankSum < -4.0 || DP > 500.0\" --filter-name snp_filter --output $SNPs/filtered_snps.$ChrName.$size.vcf.gz" ;
       JOB4_ID=$(${JOB4_TYPE} --parsable --dependency=afterok:${JOB3_ID} --wrap="${JOB4_CMD}");
   echo "$JOB4_NAME with the job id=$JOB4_ID submitted";
     ## PREPARE for Next CHUNK
     size=$(( $size + 1 ));
    else   ### Last chunk 
      Start=$(( $End - $ChunkSize +1 ));
      End=$ChrLen;
   echo "		Debug: Chromosome=$ChrName   Chunk=chunk_$size   Start:$Start  End:$End"   
      ## CombineGVCFs (Chunk by Chunk)
       JOB1_NAME="$ChrName.chunk_$size.CombineGVCFs"
       JOB1_TYPE="sbatch --partition=batch $QoS --job-name=${JOB1_NAME} --time=10-00:00:00 --output=$logs/${JOB1_NAME}.%J.out --error=$logs/${JOB1_NAME}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
       JOB1_CMD="time -p gatk CombineGVCFs --variant $INPUT --reference $REF --intervals $ChrName:$Start-$End --output $gVCF/Combine.$ChrName.$size.vcf.gz";
       JOB1_ID=$(${JOB1_TYPE} --parsable --wrap="${JOB1_CMD}");
   echo "$JOB1_NAME with the job id=$JOB1_ID submitted";
      ## CALL GenotypeGVF (Chunk by Chunk)
       JOB2_NAME="$ChrName.chunk_$size.GenotypeGVCFs"
       JOB2_TYPE="sbatch --partition=batch $QoS --job-name=${JOB2_NAME} --time=10-00:00:00 --output=$logs/${JOB2_NAME}.%J.out --error=$logs/${JOB2_NAME}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
       JOB2_CMD="time -p gatk GenotypeGVCFs --variant $gVCF/Combine.$ChrName.$size.vcf.gz --reference $REF --output $gVCF/Genotype.$ChrName.$size.vcf.gz --intervals $ChrName:$Start-$End";
       JOB2_ID=$(${JOB2_TYPE} --parsable --dependency=afterok:${JOB1_ID} --wrap="${JOB2_CMD}");
   echo "$JOB2_NAME with the job id=$JOB2_ID submitted";
      ## Select VARIANT=SNPs (Chunk by Chunk)
       JOB3_NAME="$ChrName.chunk_$size.SNPs.SelectVariants"
       JOB3_TYPE="sbatch --partition=batch $QoS --job-name=${JOB3_NAME} --time=3-00:00:00 --output=$logs/${JOB3_NAME}.%J.out --error=$logs/${JOB3_NAME}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
       JOB3_CMD="time -p gatk SelectVariants --variant $gVCF/Genotype.$ChrName.$size.vcf.gz --reference $REF -select-type SNP --output $SNPs/$ChrName.$size.vcf.gz" ;
       JOB3_ID=$(${JOB3_TYPE} --parsable --dependency=afterok:${JOB2_ID} --wrap="${JOB3_CMD}");
   echo "$JOB3_NAME with the job id=$JOB3_ID submitted";
      ## Call SNPs Filteration (Chunk by Chunk)
       JOB4_NAME="$ChrName.chunk_$size.SNPs.VariantFiltration"
       JOB4_TYPE="sbatch --partition=batch $QoS --job-name=${JOB4_NAME} --time=3-00:00:00 --output=$logs/${JOB4_NAME}.%J.out --error=$logs/${JOB4_NAME}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
       JOB4_CMD="time -p gatk VariantFiltration --variant $SNPs/$ChrName.$size.vcf --reference $REF --filter-expression \"QUAL < 30.0 || QD < 2.5 || MQ < 40.0 || MQRankSum < -12.5 || MQRankSum < -4.0 || DP > 500.0\" --filter-name snp_filter --output $SNPs/filtered_snps.$ChrName.$size.vcf.gz" ;
       JOB4_ID=$(${JOB4_TYPE} --parsable --dependency=afterok:${JOB3_ID} --wrap="${JOB4_CMD}");
   echo "$JOB4_NAME with the job id=$JOB4_ID submitted";
    fi ## END OF ALL CHUNKS within the specific CHROMOSOME 
 #  echo "Preparing for the next Chunk:" $ChrName.chunk_$size 
done ## END OF CHROMOSOME (Chr by Chr)
done < $REF.fai
  ## END OF ALL CHROMOSOMEs

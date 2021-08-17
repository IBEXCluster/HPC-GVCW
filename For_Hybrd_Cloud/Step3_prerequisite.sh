#!/bin/bash
REF=$1;
if [ "$#" -eq 1 ];
then
  if [ ! -f $REF ]; then
    printf " Your Genome reference file does not exist \n" ;
    break;
  fi;
else
 printf "\033c"
 echo " ***************************************************************************************************************************"
 echo ""
 echo " Run this script with 1 arguments:" 
 echo "      (a) Your Reference file"
 echo "  Example: "
 echo "    ./Step3_prerequisite.sh /ibex/scratch/projects/c2028/3k_Project/REF/Nipponbare_chr.fasta "
 echo ""
 echo ""
 echo " ***************************************************************************************************************************"
 exit;
fi

rm -Rf submitted_chunks_data.txt split.txt; 
## Calculate the Chromosome split 
MaxSize=0;
MaxJobs=4000;
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
       echo "$ChrName   $size   $Start  $End"  >> split.txt
       size=$(( $size + 1 ));
     else   ### Last chunk 
       Start=$(( $End - $ChunkSize +1 ));
       End=$ChrLen;
       echo "$ChrName  $size  $Start  $End"   >> split.txt
    fi ## END OF ALL CHUNKS within the specific CHROMOSOME 
 #  echo "Preparing for the next Chunk:" $ChrName.chunk_$size 
done ## END OF CHROMOSOME (Chr by Chr)
done < $REF.fai
  ## END OF ALL CHROMOSOMEs

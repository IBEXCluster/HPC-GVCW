#!/bin/bash
### Environment variables 
export INPUT=/project/k01/kathirn/3k/OUTPUT/tmpBAM

if [ $# -ne 1 ]; then
    echo "No arguments provided ..."
    echo " Provide Phase1 output file directory"
    echo " Example: ./Phase2.prerequisite.sh /project/k01/kathirn/3k/OUTPUT/tmpBAM "
    exit 1
else
     DIRECTORY=$1 ;
     rm -Rf Phase2.prefix.txt Phase2.directory.txt;
     for SAMPLE in `ls -d $DIRECTORY/*`;
     do
        echo $SAMPLE >> Phase2.directory.txt;
        PREFIX=${SAMPLE#$INPUT/} ;
        echo $PREFIX >> Phase2.prefix.txt ;
     done
fi
Rank=`cat Phase2.prefix.txt | wc -l`

echo "Set this value in your Phase2.batch file"
Nodes=$((( Rank / 24 ) + 1 ));
echo " -N $Nodes "
echo " --ntasks=$Rank "

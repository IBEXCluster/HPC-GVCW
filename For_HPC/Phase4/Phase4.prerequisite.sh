#!/bin/bash
if [ ! -f ../Phase3/distribution.txt ]; then
    echo "distribution.txt ...File not found!"
    break;
else
   cp ../Phase3/distribution.txt .
   count=`cat distribution.txt | wc -l` ;
   Nodes=$(( ((count+31) / 32) ))
   echo ""
   echo "The following values can be set in your Phase4.batch "
   echo " -N $Nodes and --ntasks=$count "
fi 

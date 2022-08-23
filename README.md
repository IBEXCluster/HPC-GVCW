# Rice-Variant-Calling

## Principal investigators (PI)

<b>Prof. Rod A. Wing</b>, <br> 
Director, Center for Desert Agriculture, <br> 
Professor, Biological and Environmental Science and Engineering, <br> 
4700 King Abdullah University of Science and Technology, <br> 
Thuwal 23955-6900, <br> 
Kingdom of Saudi Arabia <br> 

### Authors:
Nagarajan Kathiresan {nagarajan.kathiresan@kaust.edu.sa} <br> 
Yong Zhou {yong.zhou@kaust.edu.sa} <br>
Luis F. Rivera Serna {luis.riveraserna@kaust.edu.sa}

## Computational systems 

### About Shaheen 
The system has 6,174 dual sockets compute nodes based on 16 core Intel Haswell processors running at 2.3GHz. Each node has 128GB of DDR4 memory running at 2300MHz. Overall the system has a total of 197,568 processor cores and 790TB of aggregate memory. More information is available in https://www.hpc.kaust.edu.sa/content/shaheen-ii 

### About Ibex cluster
Ibex is a heterogeneous group of nodes, a mix of AMD, Intel and Nvidia GPUs with different architectures that gives the users a variety of options to work on. Overall, Ibex is made up of 488+ nodes togeter has a heterogeneous cluster and the workload is managed by the SLURM scheduler. More information is available in https://www.hpc.kaust.edu.sa/ibex

## Workflow for Rice Variant Calling 


![](https://www.hpc.kaust.edu.sa/sites/default/files/files/public/Graphical_abstract.png)

<b> Phase #1 - Data pre-processing</b> <br>
&ensp; The objective of this phase is to get the clean data from the collected rice genome samples. This includes, (a) Genome alignment using BWA MEM algorithm, (b) Update FixMate reads for the same set of genomes, Mark Duplicate and Read grouping using Genome Analysis ToolKit (GATK). <br> <br> 
<b>Phase #2 - Variant discovery </b> <br>
The objective of this phase is to call the variants per sample and generate gVCFs files. Two major steps are required in this variant discovery phase. First, the multiple sorted input files are merged into single BAM file and (re)sorted to the merged BAM using SAMTools. Second step is to call the SNPs and INDELs simultaneously via local denovo-assembly of haplotypes in an active region using GATK called “HaplotypeCaller”. At this end of this phase, we will generate a gVCF output of SNPs and INDELs. <br> <br>
<b>Phase #3 - Callset refinement </b> <br>
<b>Phase #4 - Variant tables </b> <br>

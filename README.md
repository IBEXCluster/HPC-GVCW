# HPC-GVCW

## Principal investigators (PI)

<b>Prof. Rod A. Wing</b>, <br> 
Director, Center for Desert Agriculture, <br> 
Professor, Biological and Environmental Science and Engineering, <br> 
4700 King Abdullah University of Science and Technology, <br> 
Thuwal 23955-6900, <br> 
Kingdom of Saudi Arabia <br> 

## For Pipeline support 
 Contact us: pipeline.cda@gmail.com 

### Authors:
Nagarajan Kathiresan {nagarajan.kathiresan@kaust.edu.sa} <br> 
Yong Zhou {yong.zhou@kaust.edu.sa} <br>
Zhichao Yu {2023317110021@webmail.hzau.edu.cn } <br>
Luis F. Rivera Serna {luis.riveraserna@kaust.edu.sa} <br>
Manjula Thimma {manjula.thimma@kaust.edu.sa} <br> 
Keerthana Manickam {keerthana9811@gmail.com} <br> 
Rod A Wing {rwing@ag.arizona.edu, rod.wing@kaust.edu.sa} 

## Publication: 
 DOI: https://doi.org/10.1186/s12915-024-01820-5 <br>
 PDF available here:  https://link.springer.com/content/pdf/10.1186/s12915-024-01820-5.pdf. <br> 


## Computational systems 

### About Shaheen 
The system has 6,174 dual sockets compute nodes based on 16 core Intel Haswell processors running at 2.3GHz. Each node has 128GB of DDR4 memory running at 2300MHz. Overall the system has a total of 197,568 processor cores and 790TB of aggregate memory. More information is available in https://www.hpc.kaust.edu.sa/content/shaheen-ii 

### About Ibex cluster
Ibex is a heterogeneous group of nodes, a mix of AMD, Intel and Nvidia GPUs with different architectures that gives the users a variety of options to work on. Overall, Ibex is made up of 488+ nodes togeter has a heterogeneous cluster and the workload is managed by the SLURM scheduler. More information is available in https://www.hpc.kaust.edu.sa/ibex

## Workflow for Rice Variant Calling 


![](https://www.hpc.kaust.edu.sa/sites/default/files/files/public/Graphical_abstract.png)

<b> Phase #1 - Data pre-processing</b> <br>
&ensp; The objective of this phase is to get the clean data from the collected rice genome samples. This includes, (a) Genome alignment using BWA MEM algorithm, (b) Update FixMate reads for the same set of genomes, Mark Duplicate and Read grouping using Genome Analysis ToolKit (GATK). 

![](https://www.hpc.kaust.edu.sa/sites/default/files/files/public/Phase1.png)
<br>
<b>Phase #2 - Variant discovery </b> <br>
&ensp; The objective of this phase is to call the variants per sample and generate gVCFs files. Two major steps are required in this variant discovery phase. First, the multiple sorted input files are merged into single BAM file and (re)sorted to the merged BAM using SAMTools. Second step is to call the SNPs and INDELs simultaneously via local denovo-assembly of haplotypes in an active region using GATK called “HaplotypeCaller”. At this end of this phase, we will generate a gVCF output of SNPs and INDELs. 

 ![](https://www.hpc.kaust.edu.sa/sites/default/files/files/public/Phase2.png)
 <br>
<b>Phase #3 - Callset refinement </b> <br>
&ensp; In this phase, we will combine all gVCF files from the HaplotypeCaller and generate joint genotyping across all the samples. This phase is extremely complex because of (i) Multiple samples executed across the cluster of nodes in phase #1 and phase #2 are combined (using GATK CombineGVCFs) into a single file and then, generate multi-sample joint genotyping (using GATK GenotypeGVCFs) and (ii) the CombineGVCFs and GenotypeGVCFs steps are executed in a single core using GATK. 
 ![](https://www.hpc.kaust.edu.sa/sites/default/files/files/public/Phase3a.png)
<br>
&ensp; As we know, the GATK tool is sequential due to programming limitations and the assembling of genotype across multiple samples into a single file takes extremely longer time and required huge memory when the data parallelization is absent. To address these limitations, the latest version of GATK offers variant intervals feature in CombineGVCFs and GenotypeGVCFs calls for data parallelization. 
 ![](https://www.hpc.kaust.edu.sa/sites/default/files/files/public/Phase3b.png)
<br>
<b>Phase #4 - Variant tables </b> <br>
&ensp; In this phase, the quality of genotype is enriched through variant filters and it’s also separated based on SNPs and INDELs from these independent chunks of GenotypeGVCFs files. Once all the chunks of filtered SNPs and INDELs are generated, all these partial chunks can be combined into a single file using GatherVcfs and its recommended to assemble per chromosome. The chromosome-based SNPs and INDELs are converted into variant table.
![](https://www.hpc.kaust.edu.sa/sites/default/files/files/public/Phase4.png)
<br>

## Summary of workflow steps across multiple phases 
&ensp; The below table summarizes various bioinformatics tools used in different stages of the workflow. Additionally, we provided the optimal number of CPUs used, data parallelization methods and input/output file formats are summarized.

![](https://www.hpc.kaust.edu.sa//sites/default/files/files/public/Table.png)

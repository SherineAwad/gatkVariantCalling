

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.0.2-brightgreen.svg)](https://snakemake.github.io)
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)


Snakemake Workflow for Variant Calling
==========================================

This is a GATK variant calling snakemake pipeline written by Sherine Awad. 


We are using GATK4 GVCF mode. To run the pipeline, edit the config file to match your samples file name and reference genome. 
Your samples names should be listed by default in samples.tsv file. You can change this file name in config file if neededi by editing the SAMPLES entry in the config file. 
For example, if your sample name is sample1.s_1.r_1.fq.gz, then your sample name in the samples file should be sample1.s_1. 

You need to update the config file with whether your samples are paired-end or single reads. If your samples are paired-end, then the PAIRD entry in the config file should be set to TRUE, otherwise, set the PAIRED entry in the config file to FALSE. 

You need to update your interval list, by editing the intervals.list file to list only the chromosomes of interest. You can change the name of this file by editing the config file entry INTERVALS. 


You can change the samples.tsv name in the config file. 


    snakemake -jn 

where n is the number of cores for example for 10 cores use:


    snakemake -j10 

### Use conda 

For less froodiness, use conda:


    snakemake -jn --use-conda 


For example, for 10 cores use: 

    snakemake -j10 --use-conda 

This will pull automatically the same versiosn of tools we used. Conda has to be installed in the system, in addition to snakemake. 


### Dry Run


For a dry run use: 
  
  
    snakemake -j1 -n 


and to print command in dry run use: 

  
    snakemake -j1 -n -p 


### Use Corresponding configfile:


Just update your config file to include all your sample names, edit your interval.list file to include your intervals of interest, your path, etc for example: 

  
    snakemake -j1 --configfile config-WES.yaml 
  
or: 


    snakemake -j1 configfile config-WGS.yaml 


### References 

1.  Brouard, Jean-Simon, Flavio Schenkel, Andrew Marete, and Nathalie Bissonnette. "The GATK joint genotyping workflow is appropriate for calling variants in RNA-seq experiments." Journal of animal science and biotechnology 10, no. 1 (2019): 1-6.

2. Van der Auwera, Geraldine A., Mauricio O. Carneiro, Christopher Hartl, Ryan Poplin, Guillermo Del Angel, Ami Levy‐Moonshine, Tadeusz Jordan et al. "From FastQ data to high‐confidence variant calls: the genome analysis toolkit best practices pipeline." Current protocols in bioinformatics 43, no. 1 (2013): 11-10.

3. Poplin, R., Ruano-Rubio, V., DePristo, M. A., Fennell, T. J., Carneiro, M. O., Van der Auwera, G. A., ... & Banks, E. (2018). Scaling accurate genetic variant discovery to tens of thousands of samples. BioRxiv, 201178.

4. Rausch, T., Zichner, T., Schlattl, A., Stütz, A. M., Benes, V., & Korbel, J. O. (2012). DELLY: structural variant discovery by integrated paired-end and split-read analysis. Bioinformatics, 28(18), i333-i339.

5. Li, H. (2011). A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics, 27(21), 2987-2993.

6. Eisfeldt, J., Vezzi, F., Olason, P., Nilsson, D., & Lindstrand, A. (2017). TIDDIT, an efficient and comprehensive structural variant caller for massive parallel sequencing data. F1000Research, 6.

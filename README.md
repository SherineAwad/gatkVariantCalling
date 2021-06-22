This is a GATK variant calling snakemake pipeline written by Sherine Awad. 

To run the pipeline, edit the config file to match your sample names, your reference genome then: 


    snakemake -jn 

where n is the number of cores for example for 10 cores use:


    snakemake -j10 

### Use conda 

For less froodiness, use conda:


    snakemake -jn --use-conda 


For example, for 10 cores use: 

    snakemake -j10 --use-conda 

This will pull automatically the same versiosn of tools we used. Conda has to be installed in the system, in addition to snakemake. 

### Dry run 

For a dry run use: 
  
  
    snakemake -j1 -n 


and to print command in dry run use: 

  
    snakemake -j1 -n -p 


Just update your config file to include all your sample names, edit your interval.list file to include your intervals of interest, your path, etc. 

  


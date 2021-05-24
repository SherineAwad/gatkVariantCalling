This is a GATK variant calling snakemake pipeline written by Sherine Awad. 

To run the pipeline, edit the config file to match your sample names, your reference genome then: 


    snakemake -jn 

where n is the number of cores for example for 10 cores use:


    snakemake -j10 


For a dry run use: 
  
  
    snakemake -j1 -n 


and to print command in dry run use: 

  
    snakemake -j1 -n -p 


Just update your config file to include all your sample names, edit your interval.list file to include your intervals of interest, your path, etc. 

  


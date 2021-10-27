configfile: "config.yaml"

rule all:
      input:
        expand("{genome}.fasta", genome = config['GENOME']),
        config['DBSNP'],
        config['INDELS'],
        config['GOLD_STANDARD'],
        config['Axiom_Exome'],
        config['G_phase1'],
        config['G_omni2'],
        config['hapmap_3'],
        expand("{genome}.fasta.fai", genome = config['GENOME']),
        expand("{genome}.rev.1.bt2", genome = config['GENOME']),
        expand("{genome}.rev.2.bt2", genome = config['GENOME']),
        expand("{genome}.1.bt2", genome = config['GENOME']),
        expand("{genome}.2.bt2", genome = config['GENOME']),
        expand("{genome}.3.bt2", genome = config['GENOME']),
        expand("{genome}.4.bt2", genome = config['GENOME']),
        expand("{sample}.g.vcf", sample=config['SAMPLES']), 
        directory(expand("{my_db}", my_db = config['GENOMEDB'])),
        expand("{cohort}.vcf.gz", cohort=config['ALL_VCF']),
        expand("{COHORT}.{prefix}.vcf", COHORT = config['ALL_VCF'],prefix = config['annovar_prefix'] ),
        output = expand("{COHORT}.{prefix}.filtered.vcf.gz", COHORT=config['ALL_VCF'], prefix=config['annovar_prefix'])


rule download: 
     params:
        config['GENOME']
     output:
        expand("{genome}.fasta", genome = config['GENOME']),
        config['DBSNP'],
        config['INDELS'],
        config['GOLD_STANDARD'],
        config['Axiom_Exome'],
        config['G_phase1'],
        config['G_omni2'],
        config['hapmap_3']

     shell:
         """
          wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
          wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai
          wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
          wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx
          wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz
          wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi
          wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
          wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
          wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz
          wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi
          wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz
          wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi
          wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz
          wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi
          wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz
          wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi
          wget https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.dict
         """ 

rule index:
     input: 
           expand("{genome}.fasta", genome = config['GENOME']) 
     params: 
        config['GENOME'] 
     output:
        expand("{genome}.fasta.fai", genome = config['GENOME']),        
        expand("{genome}.rev.1.bt2", genome = config['GENOME']), 
        expand("{genome}.rev.2.bt2", genome = config['GENOME']),
        expand("{genome}.1.bt2", genome = config['GENOME']),
        expand("{genome}.2.bt2", genome = config['GENOME']),
        expand("{genome}.3.bt2", genome = config['GENOME']),
        expand("{genome}.4.bt2", genome = config['GENOME']), 

     shell: 
         """
          bowtie2-build {input} {params}
          samtools faidx {input} 
         """ 


if config['PAIRED']:
    rule trim:
       input:
           r1 = "{sample}.r_1.fq.gz",
           r2 = "{sample}.r_2.fq.gz"
       output:
           "galore/{sample}.r_1_val_1.fq.gz",
           "galore/{sample}.r_2_val_2.fq.gz"
       conda: 'env/env-trim.yaml'
       shell:
           """
           mkdir -p galore
           mkdir -p fastqc
           trim_galore --gzip --retain_unpaired --trim1 --fastqc --fastqc_args "--outdir fastqc" -o galore --paired {input.r1} {input.r2}
           """
    rule tosam:
       input:
          r1 = "galore/{sample}.r_1_val_1.fq.gz",
          r2 = "galore/{sample}.r_2_val_2.fq.gz"
       params:
          genome = config['GENOME']
       output:
          "{sample}.sam"
       shell:
          "bowtie2 -x {params.genome} -1 {input.r1} -2 {input.r2} -S {output}"
else:
     rule trim:
       input:
           "{sample}.fq.gz",

       output:
           "galore/{sample}_trimmed.fq.gz",
       conda: 'env/env-trim.yaml'
       shell:
           """
           mkdir -p galore
           mkdir -p fastqc
           trim_galore --gzip --retain_unpaired --trim1 --fastqc --fastqc_args "--outdir fastqc" -o galore {input}
           """


     rule tosam:
        input:
           "galore/{sample}_trimmed.fq.gz"
        params:
           genome = config['GENOME']
        output:
           "{sample}.sam"
        shell:
           "bowtie2 -x {params.genome} -U {input} -S {output}"

rule AddRG: 
    input: 
       '{sample}.sam'
    output: 
       '{sample}.RG.sam' 
    params: 
        RG = config['RG']
    shell:
        "picard AddOrReplaceReadGroups I={input} O={output} SO=coordinate RGID=@{params} RGSM={wildcards.sample} RGPL=Illumina RGLB={wildcards.sample} RGPU={params}_{wildcards.sample} VALIDATION_STRINGENCY=SILENT" 


rule dedup: 
     input: 
         '{sample}.RG.sam'
     output:
       '{sample}.dedupped.bam',
       '{sample}.output.metrics'
     shell:
        "picard MarkDuplicates I={input} O={output[0]} CREATE_INDEX=true M={output[1]}"

rule recalibrate: 
    input:  
         sample = '{sample}.dedupped.bam', 
         genome = expand("{genome}.fasta", genome=config['GENOME'])
    output:
         '{sample}.recal_data.table',
         '{sample}.recalibrated.bam'
    params: 
        mem = "-Xmx100g",
        knownsites1 = config['DBSNP'],
        knownsites2 = config['INDELS'],
        knownsites3 = config['GOLD_STANDARD']
    shell:
       """
       gatk --java-options {params.mem} BaseRecalibrator -I {input.sample} -R {input.genome} --known-sites {params.knownsites1} --known-sites {params.knownsites2} --known-sites {params.knownsites3} -O {output[0]}
       gatk --java-options {params.mem} ApplyBQSR  -I {input.sample} -R {input.genome} --bqsr-recal-file {output[0]} -O {output[1]} 
       """ 

rule tovcf:
   input:
      "{sample}.recalibrated.bam",
      expand("{genome}.fasta", genome=config['GENOME'])
   params:
     mem_threads = {"-Xmx100g -XX:ParallelGCThreads=4"},
   output:
       "{sample}.g.vcf" 
   shell:
       """
       gatk --java-options "{params.mem_threads}" HaplotypeCaller -R {input[1]} -I {input[0]} -ERC GVCF -O {output[0]}
       """

rule GenomeDBImport: 
     input:
         sample = expand("{sample}.g.vcf", sample = config['SAMPLES']),
         genome = expand("{genome}.fasta", genome=config['GENOME'])
     params:
         INTERVALS = config['INTERVALS_FILE'],
         DB = config['GENOMEDB'],
         mem = {"-Xmx100g"},
         I =  lambda w: "-V " + " -V ".join(expand("{sample}.g.vcf", sample =config['SAMPLES']))
     output:
        directory(expand("{my_db}", my_db = config['GENOMEDB'])) 
     shell:
         """
         gatk --java-options {params.mem} GenomicsDBImport {params.I} --genomicsdb-workspace-path {params.DB} -L {params.INTERVALS} 
         """

rule jointcall: 
    input:
       DB = directory(expand("{my_db}", my_db = config['GENOMEDB'])), 
       genome = expand("{genome}.fasta", genome=config['GENOME'])
    params:
       mem = {"-Xmx100g"},
    output:
       expand("{cohort}.vcf.gz", cohort=config['ALL_VCF'])
    shell:
        """
          gatk --java-options {params.mem} GenotypeGVCFs -R {input.genome} -V gendb://{input.DB} -O {output}
        """

rule annotate: 
     input: 
         vcf = expand("{COHORT}.vcf.gz", COHORT=config['ALL_VCF'] )
     params: 
          humanDB = config['humanDB'],
          version = config['version'], 
          ANNOVAR = config['ANNOVAR'],
          output = config['ALL_VCF'] 
     output:
         expand("{COHORT}.{prefix}.vcf", COHORT = config['ALL_VCF'], prefix =config['annovar_prefix'] )
     shell:
         """ 
          {params.ANNOVAR}/table_annovar.pl {input.vcf} {params.humanDB} -buildver {params.version} -out {params.output} \
          -remove -protocol refGene,ensGene,cytoBand,exac03,gnomad_exome,avsnp147,dbnsfp33a,clinvar_20170130,revel -operation g,g,f,f,f,f,f,f,f  -nastring . -vcfinput  
         """

rule hard_filter: 
    input: 
         vcf = expand("{COHORT}.{prefix}.vcf", COHORT = config['ALL_VCF'], prefix = config['annovar_prefix'])
    params:
        qd = config['QD'], 
        qual = config['QUAL'],
        sor = config['SOR'],
        fs = config['FS'], 
        mq = config['MQ'], 
        mqranksum = config['MQRankSum'],
        readposranksum = config['ReadPosRankSum']
    output:
         output = expand("{COHORT}.{prefix}.filtered.vcf.gz", COHORT=config['ALL_VCF'], prefix=config['annovar_prefix'])
    shell: 
         """
         gatk VariantFiltration \
         -V {input} \
         -filter "QD < {params.qd}" --filter-name "QD2" \
         -filter "QUAL < {params.qual}" --filter-name "QUAL30" \
         -filter "SOR > {params.sor}" --filter-name "SOR3" \
         -filter "FS > {params.fs}" --filter-name "FS60" \
         -filter "MQ < {params.mq}" --filter-name "MQ40" \
         -filter "MQRankSum < {params.mqranksum}" --filter-name "MQRankSum-12.5" \
         -filter "ReadPosRankSum < {params.readposranksum}" --filter-name "ReadPosRankSum-8" \
         -O {output}
         """




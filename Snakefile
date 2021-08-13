configfile: "config.yaml"

ruleorder: index > trim > tosam > AddRG > dedup > recalibrate > tovcf >  GenomeDBImport > jointcall > hard_filter > annotate  


rule all:
      input:
        expand("{genome}.fasta", genome = config['GENOME']),
        config['DBSNP'],
        expand("{dbsnp}.idx", dbsnp=config['DBSNP']),
        config['INDELS'],
        expand("{indels}.tbi", indels=config['INDELS']),
        config['GOLD_STANDARD'],
        expand("{goldstd}.tbi", goldstd=config['GOLD_STANDARD']),
        config['Axiom_Exome'],
        expand("{axiom}.tbi", axiom=config['Axiom_Exome']),
        config['G_phase1'],
        expand("{gphase}.tbi", gphase=config['G_phase1']),
        config['G_omni2'],
        expand("{gomni}.tbi", gomni=config['G_omni2']),
        config['hapmap_3'],
        expand("{hapmap}.tbi", hapmap=config['hapmap_3']),
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
        expand("{COHORT}.filtered.vcf.gz", COHORT = config['ALL_VCF']),
        output = expand("{COHORT}.{prefix}.vcf", COHORT=config['ALL_VCF'], prefix=config['annovar_prefix'])


rule download: 
     params:
        config['GENOME']
     output:
        expand("{genome}.fasta", genome = config['GENOME']),
        config['DBSNP'],
        expand("{dbsnp}.idx", dbsnp=config['DBSNP']),
        config['INDELS'],
        expand("{indels}.tbi", indels=config['INDELS']),
        config['GOLD_STANDARD'],
        expand("{goldstd}.tbi", goldstd=config['GOLD_STANDARD']), 
        config['Axiom_Exome'],
        expand("{axiom}.tbi", axiom=config['Axiom_Exome']), 
        config['G_phase1'],
        expand("{gphase}.tbi", gphase=config['G_phase1']),
        config['G_omni2'],
        expand("{gomni}.tbi", gomni=config['G_omni2']),
        config['hapmap_3'],
        expand("{hapmap}.tbi", hapmap=config['hapmap_3'])
     benchmark: "logs/download.benchmark"  
     conda: 'env/env-gsutil.yaml' 
     shell:
         """
             gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta .
             gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai .
             gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf .
             gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx .
             gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz .
             gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi .
             gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz .
             gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi .
             gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz .
             gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi .
             gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz .
             gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi .
             gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz .
             gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi .
             gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz .
             gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi .
             gsutil cp gs://genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.dict .
         """ 

rule index:
     input: 
        expand("{genome}.fasta", genome = config['GENOME']),
     params: 
        config['GENOME'] 
     log: "logs/index.log"
     benchmark: "logs/index.benchmark"
     conda: 'env/env-align.yaml'
     threads: 8 
     output:
        expand("{genome}.fasta.fai", genome = config['GENOME']),        
        expand("{genome}.rev.1.bt2", genome = config['GENOME']), 
        expand("{genome}.rev.2.bt2", genome = config['GENOME']),
        expand("{genome}.1.bt2", genome = config['GENOME']),
        expand("{genome}.2.bt2", genome = config['GENOME']),
        expand("{genome}.3.bt2", genome = config['GENOME']),
        expand("{genome}.4.bt2", genome = config['GENOME'])
     shell: 
         """
          bowtie2-build {input} {params} --threads 8
          samtools faidx {input} 
         """ 

rule trim: 
    input: 
       r1 = "{sample}.r_1.fq.gz",
       r2 = "{sample}.r_2.fq.gz"
    log: "logs/{sample}.trim.log"
    benchmark: "logs/{sample}.trim.benchmark" 
    conda :'env/env-trim.yaml' 
    output: 
      val1 = "galore/{sample}.r_1_val_1.fq.gz",
      val2 = "galore/{sample}.r_2_val_2.fq.gz"
    shell: 
        """
         mkdir -p galore
         mkdir -p fastqc
         trim_galore --gzip --retain_unpaired --trim1 --fastqc --fastqc_args "--outdir fastqc" -o galore --paired {input.r1} {input.r2}
        """ 

rule tosam:
    input:
        expand("{genome}.fasta.fai", genome = config['GENOME']),
        expand("{genome}.rev.1.bt2", genome = config['GENOME']),
        expand("{genome}.rev.2.bt2", genome = config['GENOME']),
        expand("{genome}.1.bt2", genome = config['GENOME']),
        expand("{genome}.2.bt2", genome = config['GENOME']),
        expand("{genome}.3.bt2", genome = config['GENOME']),
        expand("{genome}.4.bt2", genome = config['GENOME']),
        r1 = "galore/{sample}.r_1_val_1.fq.gz",
        r2 = "galore/{sample}.r_2_val_2.fq.gz"
    params:
        genome = config['GENOME']
    threads: 4
    log: "logs/{sample}.align.log"
    benchmark: "logs/{sample}.tosam.benchmark"
    conda: 'env/env-align.yaml' 
    output:
        "{sample}.sam"
    shell:
        "bowtie2 -x {params.genome} -1 {input.r1} -2 {input.r2} -S {output} -p 4"

rule AddRG: 
    input: 
       '{sample}.sam'
    output: 
       '{sample}.RG.sam' 
    log: "logs/{sample}.addRG.log"
    benchmark: "logs/{sample}.AddRG.benchmark"
    conda: 'env/env-picard.yaml'  
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
     log: "logs/{sample}.dedup.log"
     benchmark: "logs/{sample}.dedup.benchmark" 
     conda: 'env/env-picard.yaml' 
     shell:
        "picard MarkDuplicates I={input} O={output[0]} CREATE_INDEX=true M={output[1]}"

rule recalibrate: 
    input:  
         sample = '{sample}.dedupped.bam', 
         genome = expand("{genome}.fasta", genome=config['GENOME'])
    log: "logs/{sample}.recalibrate.log" 
    benchmark: "logs/{sample}.recalibrate.benchmark"
    conda: 'env/env-gatk.yaml'
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
     mem_threads = {"-Xmx100g -XX:ParallelGCThreads=4"}
   log: "logs/{sample}.tovcf.log"
   benchmark: "logs/{sample}.tovcf.benchmark"  
   conda: 'env/env-gatk.yaml'
   output:
       "{sample}.g.vcf" 
   shell:
       """
       gatk --java-options "{params.mem_threads}" HaplotypeCaller -R {input[1]} -I {input[0]} -ERC GVCF -O {output[0]} -G StandardAnnotation -G AS_StandardAnnotation
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
     log: "logs/dbimport.log"
     benchmark: "logs/DBImport.benchmark"
     conda: 'env/env-gatk.yaml'
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
       mem = {"-Xmx100g"}
    log: "logs/jointcall.log"
    benchmark: "logs/jointcall.benchmark"
    conda: 'env/env-gatk.yaml' 
    output:
       expand("{cohort}.vcf.gz", cohort=config['ALL_VCF'])
    shell:
        """
          gatk --java-options {params.mem} GenotypeGVCFs -R {input.genome} -V gendb://{input.DB} -O {output} -G StandardAnnotation -G AS_StandardAnnotation
        """

rule annotate: 
     input: 
         vcf = expand("{COHORT}.filtered.vcf.gz", COHORT=config['ALL_VCF'] )
     params: 
          humanDB = config['humanDB'],
          version = config['version'], 
          ANNOVAR = config['ANNOVAR'],
          output = config['ALL_VCF']
     log: "logs/annotate.log" 
     benchmark: "logs/annotate.benchmark" 
     output:
         expand("{COHORT}.{prefix}.vcf", COHORT = config['ALL_VCF'], prefix =config['annovar_prefix'] )
     shell:
         """
          {params.ANNOVAR}/table_annovar.pl {input.vcf} {params.humanDB} -buildver {params.version} -out {params.output} \
          -remove -protocol refGene,ensGene,cytoBand,exac03,gnomad_exome,avsnp147,dbnsfp33a,clinvar_20170130,revel -operation g,g,f,f,f,f,f,f,f  -nastring . -vcfinput  
         """

rule hard_filter: 
    input: 
         vcf = expand("{COHORT}.vcf.gz", COHORT = config['ALL_VCF'])
    params:
        qd = config['QD'], 
        qual = config['QUAL'],
        sor = config['SOR'],
        fs = config['FS'], 
        mq = config['MQ'], 
        mqranksum = config['MQRankSum'],
        readposranksum = config['ReadPosRankSum']
    log: "logs/filter.log"
    benchmark: "logs/filter.benchmark"
    conda: 'env/env-gatk.yaml'
    output:
         output = expand("{COHORT}.filtered.vcf.gz", COHORT=config['ALL_VCF'])
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




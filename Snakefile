ruleorder: index > trim > tosam > AddRG > dedup > recalibrate > tovcf >  GenomeDBImport > jointcall  > select_SNPs > select_INDELs > filter_SNPs > filter_INDELs > annotate
 
with open(config['INTERVALS_FILE']) as fp:
    INTERVALS= fp.read().splitlines()
print(INTERVALS)

with open(config['SAMPLES']) as fp:
    SAMPLES= fp.read().splitlines()
print(SAMPLES)



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
        expand("{sample}.sam", sample = SAMPLES),
        expand("{sample}.RG.sam", sample = SAMPLES), 
        expand("{sample}.dedupped.bam", sample = SAMPLES), 
        expand("{sample}.recalibrated.bam", sample =SAMPLES),
        expand("{sample}.g.vcf", sample=SAMPLES), 
        expand("DB_{interval}", interval =INTERVALS),
        expand("{interval}.vcf", interval= INTERVALS),
        expand("{COHORT}.vcf.gz", COHORT = config['ALL_VCF']),
        expand("{COHORT}_filtered_snps.vcf.gz", COHORT = config['ALL_VCF']), 
        expand("{COHORT}_filtered_indels.vcf.gz", COHORT = config['ALL_VCF']), 
        expand("{COHORT}_SNPs.{prefix}.txt", COHORT = config['ALL_VCF'], prefix =config['annovar_prefix'] ),
        expand("{COHORT}_INDELs.{prefix}.txt", COHORT = config['ALL_VCF'], prefix =config['annovar_prefix'] ),
 

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

if config['PAIRED']:
    rule trim:
       input:
           r1 = "{sample}_r1.fastq.gz",
           r2 = "{sample}_r2.fastq.gz"
       params:
           threads= config['THREADS']
       output:
           "galore/{sample}_r1_val_1.fq.gz",
           "galore/{sample}_r2_val_2.fq.gz"
       conda: 'env/env-trim.yaml'
       shell:
           """
           mkdir -p galore
           mkdir -p fastqc
           trim_galore --gzip --retain_unpaired --trim1 --fastqc --fastqc_args "--outdir fastqc" -o galore --paired {input.r1} {input.r2} --cores {params.threads}
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
          r1 = "galore/{sample}_r1_val_1.fq.gz",
          r2 = "galore/{sample}_r2_val_2.fq.gz" 
       params:
          genome = config['GENOME'],
          threads = config['THREADS']
       benchmark: "logs/{sample}.bowtie2.benchmark"
       conda: 'env/env-align.yaml'
       output:
          "{sample}.sam"
       shell:
          "bowtie2 -x {params.genome} -1 {input.r1} -2 {input.r2} -S {output} -p {params.threads}" 

else:
     rule trim:
       input:
           "{sample}_r1.fastq.gz",
       params:
          threads = config['THREADS']
       output:
           "galore/{sample}_trimmed.fq.gz",
       conda: 'env/env-trim.yaml'
       shell:
           """
           mkdir -p galore
           mkdir -p fastqc
           trim_galore --gzip --retain_unpaired --trim1 --fastqc --fastqc_args "--outdir fastqc" -o galore {input} --cores {params.threads}
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
           "galore/{sample}_trimmed.fq.gz"
        params:
           genome = config['GENOME'],
           threads = config['THREADS']
        benchmark: "logs/{sample}.bowtie2.benchmark"
        conda: 'env/env-align.yaml'
        output:
           "{sample}.sam"
        shell:
           "bowtie2 -x {params.genome} -U {input} -S {output}  -p {params.threads}"


rule AddRG:
    input:
       '{sample}.sam',
    output:
       '{sample}.RG.sam'
    log: "logs/{sample}.addRG.log"
    benchmark: "logs/{sample}.AddRG.benchmark"
    conda: 'env/env-picard.yaml'
    params:
        PL = config['PL'], 
        fastq_file = '{sample}_r1.fastq.gz'
    shell:
        """
         python src/RG.py -s {input} -f {params.fastq_file} -o {output} -p {params.PL}
        """

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
        threads = "-XX:ParallelGCThreads=2",
        knownsites1 = config['DBSNP'],
        knownsites2 = config['INDELS'],
        knownsites3 = config['GOLD_STANDARD']
    shell:
       """
       gatk --java-options {params.mem} {params.threads} BaseRecalibrator -I {input.sample} -R {input.genome} --known-sites {params.knownsites1} --known-sites {params.knownsites2} --known-sites {params.knownsites3} -O {output[0]}
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
   conda: "env/env-gatk.yaml"
   output:
       "{sample}.g.vcf" 
   shell:
       """
       gatk --java-options "{params.mem_threads}" HaplotypeCaller -R {input[1]} -I {input[0]} -ERC GVCF -O {output} -G StandardAnnotation -G AS_StandardAnnotation
       
       """

rule GenomeDBImport:
     input:
        expand("logs/{sample}.tovcf.log", sample = SAMPLES) 
     params: 
         lambda w: "-V " + " -V ".join(expand("{sample}.g.vcf", sample =SAMPLES)),
         interval = "{interval}",
         mem = {"-Xmx100g"}
     log: 
        "logs/dbimport_{interval}.log"
     benchmark: 
        "logs/dbimport_{interval}.benchmark"
     conda: 
        "env/env-gatk.yaml"
     output: 
         directory("DB_{interval}") 
     shell:
         """
           gatk --java-options {params.mem} GenomicsDBImport {params[0]} --genomicsdb-workspace-path {output} -L {params.interval}
         """

rule jointcall: 
    input:
        DB = "DB_{chr}",
        genome = expand("{genome}.fasta", genome=config['GENOME'])
    params:
        lambda w: "-V " + " -V ".join(expand("{sample}.g.vcf", sample =SAMPLES)),
        V = "gendb://DB_{chr}",
        mem = {"-Xmx100g"}
    log: "logs/jointcall_{chr}.log"
    benchmark: "logs/jointcall_{chr}.benchmark"
    conda: 'env/env-gatk.yaml' 
    output:
       "{chr}.vcf"
    shell:
        """
          gatk --java-options {params.mem} GenotypeGVCFs -R {input.genome} -V {params.V} -O {output} -G StandardAnnotation -G AS_StandardAnnotation
        """

rule merge:
     input: 
         expand("{chr}.vcf" , chr=INTERVALS)
     params:  
          I =  lambda w: "I= " + " I= ".join(expand("{chr}.vcf", chr =INTERVALS))
     output: 
          expand("{COHORT}.vcf.gz", COHORT = config['ALL_VCF'])
     shell: 
         """
         picard MergeVcfs {params.I} O= {output}
         """

rule select_SNPs: 
    input: 
        vcf = expand("{COHORT}.vcf.gz", COHORT = config['ALL_VCF']) 
    output: 
         "{COHORT}_SNPs.vcf.gz" 
    shell: 
         """
         gatk SelectVariants -V {input} -select-type SNP -O {output} 
         """ 

rule select_INDELs:
    input:
        vcf = expand("{COHORT}.vcf.gz", COHORT = config['ALL_VCF'])
    output:
         "{COHORT}_INDELs.vcf.gz" 
    shell:
         """
         gatk SelectVariants -V {input} -select-type INDEL -O {output}
         """

rule filter_SNPs: 
    input: 
       "{COHORT}_SNPs.vcf.gz" 
    output:   
       "{COHORT}_filtered_snps.vcf.gz"
    params: 
       qd = config['QD'],
       qual = config['QUAL'],
       sor = config['SNP_SOR'],
       fs = config['SNP_FS'],
       mq = config['SNP_MQ'], 
       mqranksum = config['SNP_MQRankSum'], 
       readposranksum = config['SNP_ReadPosRankSum']
    log: "logs/{COHORT}_snps_filter.log"
    benchmark: "logs/{COHORT}_snps_filter.benchmark"
    conda: 'env/env-gatk.yaml'  
    shell:
        """ 
	gatk VariantFiltration \
        -V {input} \
        -filter "QD < {params.qd}" --filter-name "QD2" \
        -filter "QUAL < {params.qual}" --filter-name "QUAL30" \
        -filter "SOR > {params.sor}" --filter-name "SOR3" \
        -filter "FS > {params.fs}"  --filter-name "FS60" \
        -filter "MQ <  {params.mq}" --filter-name "MQ40" \
        -filter "MQRankSum < {params.mqranksum}" --filter-name "MQRankSum-12.5" \
        -filter "ReadPosRankSum < {params.readposranksum}" --filter-name "ReadPosRankSum-8" \
        -O {output} 
       """



rule filter_INDELs:  
    input:
       "{COHORT}_INDELs.vcf.gz"
    output:
       "{COHORT}_filtered_indels.vcf.gz"
    params: 
       qd = config['QD'],
       qual = config['QUAL'],
       fs = config['INDEL_FS'],
       readposranksum = config['INDEL_ReadPosRankSum']
    log: "logs/{COHORT}_indels_filter.log"
    benchmark: "logs/{COHORT}_indels_filter.benchmark"
    conda: "env/env-gatk.yaml"
    shell: 
       """ 
       gatk VariantFiltration -V {input} \ 
       -filter "QD < {params.qd}" --filter-name "QD2" \
       -filter "QUAL < {params.qual}" --filter-name "QUAL30" \
       -filter "FS > {params.fs}" --filter-name "FS200" \
       -filter "ReadPosRankSum < {params.readposranksum}" --filter-name "ReadPosRankSum-20" \ 
       -O {output} 
       """ 

rule annotate:
     input:
         snps = expand("{COHORT}_filtered_snps.vcf.gz", COHORT=config['ALL_VCF'] ), 
         indels = expand("{COHORT}_filtered_indels.vcf.gz", COHORT = config['ALL_VCF']) 
     params:
          humanDB = config['humanDB'],
          version = config['version'],
          ANNOVAR = config['ANNOVAR'],
          snps_prefix = expand("{COHORT}_SNPs", COHORT= config['ALL_VCF']), 
          indels_prefix =  expand("{COHORT}_INDELs", COHORT= config['ALL_VCF']),
     log: "logs/annotate.log"
     benchmark: "logs/annotate.benchmark"
     output:
         expand("{COHORT}_SNPs.{prefix}.txt", COHORT = config['ALL_VCF'], prefix =config['annovar_prefix'] ),
         expand("{COHORT}_SNPs.{prefix}.vcf", COHORT = config['ALL_VCF'], prefix =config['annovar_prefix'] ),
         expand("{COHORT}_INDELs.{prefix}.txt", COHORT = config['ALL_VCF'], prefix =config['annovar_prefix'] ),
         expand("{COHORT}_INDELs.{prefix}.vcf", COHORT = config['ALL_VCF'], prefix =config['annovar_prefix'] ) 
     shell:
         """
          {params.ANNOVAR}/table_annovar.pl {input.snps} {params.humanDB} -buildver {params.version} -out {params.snps_prefix} \
          -remove -protocol refGene,ensGene,cytoBand,exac03,gnomad_exome,avsnp147,dbnsfp33a,clinvar_20170130,revel -operation g,g,f,f,f,f,f,f,f  -nastring . -vcfinput
          
          {params.ANNOVAR}/table_annovar.pl {input.indels} {params.humanDB} -buildver {params.version} -out {params.indels_prefix} \
          -remove -protocol refGene,ensGene,cytoBand,exac03,gnomad_exome,avsnp147,dbnsfp33a,clinvar_20170130,revel -operation g,g,f,f,f,f,f,f,f  -nastring . -vcfinput
         """







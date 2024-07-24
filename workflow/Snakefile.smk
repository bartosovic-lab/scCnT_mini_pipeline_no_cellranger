rule main:
    input:
        # expand("{out}/{sample}/trimming/{sample}_{R}.fastq.gz", out=config['general']['outdir'], sample=config['samples'],R=['R1','R2','R3']),
        # expand("{out}/{sample}/trimming/{sample}_{R}.fastq.gz", out=config['general']['outdir'], sample=config['samples'],R=['R1','R2']),
        expand("{out}/{sample}/mapping/{sample}_fixed.bam", out=config['general']['outdir'], sample=config['samples']),
# Trim R1 and R3 reads
rule trim_trim_galore:
  input:
    lambda wildcards: config['samples'][wildcards.sample]['R1_PATH'], 
    lambda wildcards: config['samples'][wildcards.sample]['R3_PATH']
  output:
    R1 = temp("{out}/{sample}/trimming/{sample}_R1_temp.fastq.gz"),     # TODO: Mark with temp later on
    R3 = temp("{out}/{sample}/trimming/{sample}_R3_temp.fastq.gz"),     # TODO: Mark with temp later on
  params:
    outdir = "{out}/{sample}/trimming/"
  threads: 8
  resources:
    mem_mb = 8000
  shell:
    "source ~/.zshrc && conda activate trimming && "
    "trim_galore --cores {threads} --fastqc --paired -o {params.outdir} --basename {wildcards.sample} {input} && "
    "mv {wildcards.out}/{wildcards.sample}/trimming/{wildcards.sample}_val_1.fq.gz {output.R1} && "
    "mv {wildcards.out}/{wildcards.sample}/trimming/{wildcards.sample}_val_2.fq.gz {output.R3} "

rule fix_R2_after_trimming:
    input:
        R1 = "{out}/{sample}/trimming/{sample}_R1_temp.fastq.gz",
        R3 = "{out}/{sample}/trimming/{sample}_R3_temp.fastq.gz",
        R2 = lambda wildcards: config['samples'][wildcards.sample]['R2_PATH']
    output:
        R2 = temp("{out}/{sample}/trimming/{sample}_R2_temp.fastq.gz")
    params:
        script = workflow.basedir + "/scripts/filter_R2_according_to_R1.py"
    shell:
        'python3 {params.script} --fastq_R1 {input.R1} --fastq_R2 {input.R2} --out {output.R2}'

 
rule add_R2_barcode_to_R1_R3_read_name:
    input:
        R1 = "{out}/{sample}/trimming/{sample}_R1_temp.fastq.gz",
        R3 = "{out}/{sample}/trimming/{sample}_R3_temp.fastq.gz",
        R2 = "{out}/{sample}/trimming/{sample}_R2_temp.fastq.gz"
    output:
        R1 = "{out}/{sample}/trimming/{sample}_R1.fastq.gz",
        R3 = "{out}/{sample}/trimming/{sample}_R3.fastq.gz"
    params:
        script = workflow.basedir + "/scripts/add_barcode_to_read_name.py",
    shell:
        "python3 {params.script} --fastq_R1 {input.R1} --fastq_R2 {input.R2} --fastq_R3 {input.R3} --out_R1 {output.R1} --out_R3 {output.R3}"

rule map:
    input:
        R1 = "{out}/{sample}/trimming/{sample}_R1.fastq.gz",
        R2 = "{out}/{sample}/trimming/{sample}_R3.fastq.gz"
    output:
        bam = "{out}/{sample}/mapping/{sample}.bam"
    params:
        outdir = "{out}/{sample}/mapping/",
        index  = config['general']['bwa_index']
    threads: 8
    shell:
        "bwa mem -t {threads} {params.index} {input.R1} {input.R2} | samtools view -bS - > {output.bam}"

rule fix_tags_in_bam:
    input:
        bam = "{out}/{sample}/mapping/{sample}.bam"
    output:
        bam = "{out}/{sample}/mapping/{sample}_fixed.bam"
    params:
        script = workflow.basedir + "/scripts/fix_bam_tags.py"
    shell:
        "python {params.script} --bam {input.bam} --out {output.bam}"
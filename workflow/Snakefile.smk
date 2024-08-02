shell.executable("/bin/zsh")

rule main:
    input:
        expand("{out}/{sample}/mapping/{sample}_possorted.bam", out=config['general']['outdir'], sample=config['samples']),
        expand("{out}/{sample}/mapping/{sample}_fragments.tsv.gz", out=config['general']['outdir'], sample=config['samples']),
        expand("{out}/{sample}/mapping/{sample}_mod_fragments.tsv.gz", out=config['general']['outdir'], sample=config['samples']),
        expand("{out}/{sample}/mapping/{sample}.bw", out=config['general']['outdir'], sample=config['samples']),
        expand("{out}/{sample}/peaks/macs2_narrow/{sample}_peaks.narrowPeak", out=config['general']['outdir'], sample=config['samples']),
        expand("{out}/{sample}/peaks/macs2_broad/{sample}_peaks.broadPeak", out=config['general']['outdir'], sample=config['samples']),
        expand("{out}/{sample}/metrics/{sample}_metadata.csv", out=config['general']['outdir'], sample=config['samples']),
        expand("{out}/{sample}/matrix_peaks/", out=config['general']['outdir'], sample=config['samples']),
        expand("{out}/{sample}/matrix_genomic_bins_{binsize}/", out=config['general']['outdir'], sample=config['samples'], binsize=config['general']['binsizes']),
        expand("{out}/{sample}/matrix_GA/", out=config['general']['outdir'], sample=config['samples']),
        "{out}/merged/fragments.tsv.gz".format(out=config['general']['outdir']),
        "{out}/merged/merged_matrix_peaks/".format(out=config['general']['outdir']),
        expand("{out}/merged/merged_matrix_genomic_bins_{binsize}/", out=config['general']['outdir'], binsize=config['general']['binsizes']),
        "{out}/merged/merged_matrix_GA/".format(out=config['general']['outdir']),
        

# Trim R1 and R3 reads
rule trim_trim_galore:
  input:
    lambda wildcards: config['samples'][wildcards.sample]['R1_PATH'], 
    lambda wildcards: config['samples'][wildcards.sample]['R3_PATH']
  output:
    R1 = temp("{out}/{sample}/trimming/{sample}_R1_temp.fastq.gz"), 
    R3 = temp("{out}/{sample}/trimming/{sample}_R3_temp.fastq.gz"),     
  params:
    outdir = "{out}/{sample}/trimming/"
  threads: 4
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
        R1 = temp("{out}/{sample}/trimming/{sample}_R1.fastq.gz"),
        R3 = temp("{out}/{sample}/trimming/{sample}_R3.fastq.gz")
    params:
        script = workflow.basedir + "/scripts/add_barcode_to_read_name.py",
    threads: 4
    shell:
        "python3 {params.script} --fastq_R1 {input.R1} --fastq_R2 {input.R2} --fastq_R3 {input.R3} --out_R1 {output.R1} --out_R3 {output.R3} --threads {threads}"

rule map:
    input:
        R1 = "{out}/{sample}/trimming/{sample}_R1.fastq.gz",
        R2 = "{out}/{sample}/trimming/{sample}_R3.fastq.gz"
    output:
        bam = temp("{out}/{sample}/mapping/{sample}.bam")
    params:
        outdir = "{out}/{sample}/mapping/",
        index  = config['general']['bwa_index']
    threads: 4
    shell:
        "bwa mem -t {threads} {params.index} {input.R1} {input.R2} | samtools view -bS - > {output.bam}"

rule fix_tags_in_bam:
    input:
        bam = "{out}/{sample}/mapping/{sample}.bam"
    output:
        bam = temp("{out}/{sample}/mapping/{sample}_fixed.bam")
    params:
        script = workflow.basedir + "/scripts/fix_bam_tags.py"
    threads: 1
    shell:
        "python {params.script} --bam {input.bam} --out {output.bam}"

rule sort_bam:
    input:
        bam = "{out}/{sample}/mapping/{sample}_fixed.bam"
    output:
        bam   = "{out}/{sample}/mapping/{sample}_possorted.bam",
        index = "{out}/{sample}/mapping/{sample}_possorted.bam.bai",
    threads: 4
    shell:
        "samtools sort -@ {threads} -o {output.bam} {input.bam} && samtools index {output.bam}"

rule bam_to_fragments:
    input:
        bam   = "{out}/{sample}/mapping/{sample}_possorted.bam",
        index = "{out}/{sample}/mapping/{sample}_possorted.bam.bai",
    output:
        fragments = "{out}/{sample}/mapping/{sample}_fragments.tsv.gz"
    threads: 4
    shell:
        'sinto fragments -b {input.bam} -f {output.fragments}_temp -m 30 -p {threads} && '
        'sort -k1,1 -k2,2n {output.fragments}_temp > {output.fragments}_sorted_temp && '
        'bgzip -@ {threads} -c {output.fragments}_sorted_temp > {output.fragments} && '
        'tabix -p bed {output.fragments} && '
        'rm {output.fragments}_temp {output.fragments}_sorted_temp'
    
rule fragments_with_sample_id: 
    input:
        fragments = "{out}/{sample}/mapping/{sample}_fragments.tsv.gz"
    output:
        fragments = "{out}/{sample}/mapping/{sample}_mod_fragments.tsv.gz"
    params:
        script = workflow.basedir + "/scripts/add_sample_name_to_cell_barcode.py"
    threads: 4
    shell:
        'python3 {params.script} --fragments {input.fragments} --out {output.fragments} --sample {wildcards.sample} --threads {threads}'

rule cells_with_sample_id:
    input:
        cells = "{out}/{sample}/metrics/{sample}_cells.txt"
    output:
        cells = "{out}/{sample}/metrics/{sample}_mod_cells.txt"
    threads: 1
    shell:
        """awk -v sample={wildcards.sample} '{{print $0"_"sample}}' {input.cells} > {output.cells}"""

rule bam_to_bw:
    input:
        bam = "{out}/{sample}/mapping/{sample}_possorted.bam",
        index = "{out}/{sample}/mapping/{sample}_possorted.bam.bai",
    output:
        bw = "{out}/{sample}/mapping/{sample}.bw"
    threads: 4
    shell:
        'source ~/.zshrc && conda activate deeptools && '
        'bamCoverage -b {input.bam} -o {output.bw} --binSize 10 --normalizeUsing RPKM --extendReads 200 --centerReads -p {threads}'

rule call_peaks_narrow:
    input:
        bam = "{out}/{sample}/mapping/{sample}_possorted.bam",
        index = "{out}/{sample}/mapping/{sample}_possorted.bam.bai",
    output:
        peaks_narrow = "{out}/{sample}/peaks/macs2_narrow/{sample}_peaks.narrowPeak"
    threads: 1
    shell:
        'source ~/.zshrc && conda activate macs2 && '
        'macs2 callpeak -t {input.bam} -f BAM -g hs -n {wildcards.sample} -q 0.01 --outdir {wildcards.out}/{wildcards.sample}/peaks/macs2_narrow/ --nomodel'

rule call_peaks_broad:
    input:
        bam = "{out}/{sample}/mapping/{sample}_possorted.bam",
        index = "{out}/{sample}/mapping/{sample}_possorted.bam.bai",
    output:
        peaks_broad = "{out}/{sample}/peaks/macs2_broad/{sample}_peaks.broadPeak"
    threads: 1
    shell:
        'source ~/.zshrc && conda activate macs2 && '
        'macs2 callpeak -t {input.bam} -f BAM -g hs -n {wildcards.sample} -q 0.01 --broad --outdir {wildcards.out}/{wildcards.sample}/peaks/macs2_broad/ --llocal 100000 --keep-dup 1 --broad-cutoff 0.1 --max-gap 1000 --nomodel' 


rule barcode_metrics_peaks:
    input:
        bam = "{out}/{sample}/mapping/{sample}_possorted.bam",
        peaks = "{out}/{sample}/peaks/macs2_broad/{sample}_peaks.broadPeak"
    output:
        overlap = "{out}/{sample}/metrics/{sample}_peak_region_barcode_counts.txt"
    params:
        script    = workflow.basedir + "/scripts/summarise_barcode_counts.py",
        whitelist = workflow.basedir + "/../Resources/737K-cratac-v1.txt.gz"
    threads: 1
    shell:
        'bedtools intersect -abam {input.bam} -b {input.peaks} -u > {output.overlap}_temp && '
        'python {params.script} --bam {output.overlap}_temp --out {output.overlap} --whitelist {params.whitelist} &&'
        'rm {output.overlap}_temp'

rule barcode_metrics_overall:
    input:
        bam = "{out}/{sample}/mapping/{sample}_possorted.bam"
    output:
        counts = "{out}/{sample}/metrics/{sample}_barcode_counts.txt"
    params:
        script    = workflow.basedir + "/scripts/summarise_barcode_counts.py",
        whitelist = workflow.basedir + "/../Resources/737K-cratac-v1.txt.gz"
    threads: 1
    shell:
        'python {params.script} --bam {input.bam} --out {output.counts} --whitelist {params.whitelist}'

rule cell_picking:
    input:
        metrics_peak = "{out}/{sample}/metrics/{sample}_peak_region_barcode_counts.txt",
        metrics_all  = "{out}/{sample}/metrics/{sample}_barcode_counts.txt",
    output:
        metadata = "{out}/{sample}/metrics/{sample}_metadata.csv",
        cells    = "{out}/{sample}/metrics/{sample}_cells.txt",
    threads: 1
    shell:
        'Rscript {workflow.basedir}/scripts/pick_cells_light.R --metrics_all {input.metrics_all} --metrics_peak {input.metrics_peak} --out {output.metadata} --cells_list {output.cells}'


rule create_matrix_peaks:
    input:
        cells = "{out}/{sample}/metrics/{sample}_cells.txt",
        fragments = "{out}/{sample}/mapping/{sample}_fragments.tsv.gz",
        peaks = "{out}/{sample}/peaks/macs2_broad/{sample}_peaks.broadPeak"
    output:
        matrix = directory("{out}/{sample}/matrix_peaks/")
    threads: 1
    shell:
        # Requires f2m in $PATH
        # https://github.com/stuart-lab/f2m
        "f2m --fragments {input.fragments} --bed {input.peaks} --cells {input.cells} --outdir {output.matrix}"

rule generate_genomic_bins_annotation:
    input:
        genome = config['general']['chrom_sizes']
    output:
        bed = temp("{out}/genomic_bins_{binsize}.bed")
    threads: 1
    shell:
        'bedtools makewindows -g {input.genome} -w {wildcards.binsize} > {output.bed}'


rule create_matrix_genomic_bins:
    input:
        cells = "{out}/{sample}/metrics/{sample}_cells.txt",
        fragments = "{out}/{sample}/mapping/{sample}_fragments.tsv.gz",
        bins = "{out}/genomic_bins_{binsize}.bed"
    output:
        matrix = directory("{out}/{sample}/matrix_genomic_bins_{binsize}/")
    threads: 1
    shell:
        # Requires f2m in $PATH
        # https://github.com/stuart-lab/f2m
        "f2m --fragments {input.fragments} --bed {input.bins} --cells {input.cells} --outdir {output.matrix}"


rule download_GA_annotation:
    output:
        annotation = "{out}/GA_annotation.bed"
    shell:
        'Rscript {workflow.basedir}/scripts/get_genebody_and_promoter_coordinates.R --out {output.annotation}'

rule create_matrix_GA:
    input:
        cells = "{out}/{sample}/metrics/{sample}_cells.txt",
        fragments = "{out}/{sample}/mapping/{sample}_fragments.tsv.gz",
        annotation = "{out}/GA_annotation.bed"
    output:
        matrix = directory("{out}/{sample}/matrix_GA/")
    threads: 1
    shell:
        # Requires f2m in $PATH
        # https://github.com/stuart-lab/f2m
        "f2m --fragments {input.fragments} --bed {input.annotation} --cells {input.cells} --outdir {output.matrix}"

rule merge_fragment_files:
    input:
        fragments = expand("{out}/{sample}/mapping/{sample}_mod_fragments.tsv.gz", out=config['general']['outdir'], sample=config['samples'])
    output:
        fragments = "{out}/merged/fragments.tsv.gz"
    threads: 4
    shell:
        'gunzip -c {input.fragments} | bedtools sort -i - > {output.fragments}_temp && '
        'bgzip -@ {threads} -c {output.fragments}_temp > {output.fragments} && '
        'tabix -p bed {output.fragments} && '
        'rm {output.fragments}_temp'

rule merge_cells:
    input:
        cells = expand("{out}/{sample}/metrics/{sample}_mod_cells.txt", out=config['general']['outdir'], sample=config['samples'])
    output:
        cells = "{out}/merged/cells.txt"
    threads: 1
    shell:
        'cat {input.cells} > {output.cells}'

rule merge_peaks:
    input:
        peaks = expand("{out}/{sample}/peaks/macs2_broad/{sample}_peaks.broadPeak", out=config['general']['outdir'], sample=config['samples'])
    output:
        peaks = "{out}/merged/peaks.broadPeak"
    threads: 1
    shell:
        'cat {input.peaks} | bedtools sort -i - > {output.peaks}_temp && '
        'bedtools merge -i {output.peaks}_temp > {output.peaks} && '
        'rm {output.peaks}_temp'

rule merged_peak_matrix:
    input:
        fragments = "{out}/merged/fragments.tsv.gz",
        peaks = "{out}/merged/peaks.broadPeak",
        cells = "{out}/merged/cells.txt"
    output:
        matrix = directory("{out}/merged/merged_matrix_peaks/")
    threads: 1
    shell:
        'f2m --fragments {input.fragments} --bed {input.peaks} --cells {input.cells} --outdir {output.matrix}'

rule merged_bin_matrix:
    input:
        fragments = "{out}/merged/fragments.tsv.gz",
        bins = "{out}/genomic_bins_{binsize}.bed",
        cells = "{out}/merged/cells.txt"
    output:
        matrix = directory("{out}/merged/merged_matrix_genomic_bins_{binsize}/")
    threads: 1
    shell:
        'f2m --fragments {input.fragments} --bed {input.bins} --cells {input.cells} --outdir {output.matrix}'

rule merged_GA_matrix:
    input:
        fragments = "{out}/merged/fragments.tsv.gz",
        annotation = "{out}/GA_annotation.bed",
        cells = "{out}/merged/cells.txt"
    output:
        matrix = directory("{out}/merged/merged_matrix_GA/")
    threads: 1
    shell:
        'f2m --fragments {input.fragments} --bed {input.annotation} --cells {input.cells} --outdir {output.matrix}'
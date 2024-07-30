import os

shell.executable("/bin/zsh")


fraction = [i for i in range(5,100,5)]
prefix = '/Users/marek/human/out'

rule all:
    input:
        # expand("{prefix}/{sample}_down_{fraction}/mapping/{sample}_possorted.bam",sample='P32269_1001',fraction = fraction, prefix = prefix),
        # expand('{prefix}/{sample}_down_{fraction}/metrics/{sample}_barcode_counts.tsv',sample='P32269_1001',fraction = fraction, prefix = prefix),
        # expand('{prefix}/{sample}_down_{fraction}/metrics/{sample}_peak_region_barcode_counts.txt',sample='P32269_1001',fraction = fraction, prefix = prefix),
        expand('{prefix}/{sample}_down_{fraction}/mapping/{sample}_fragments.tsv.gz',sample='P32269_1001',fraction = fraction, prefix = prefix),
        expand("{prefix}/{sample}_down_{fraction}/metrics/{sample}_metadata.csv",sample='P32269_1001',fraction = fraction, prefix = prefix),
        expand("{prefix}/{sample}_down_{fraction}/matrix_peaks/",sample='P32269_1001',fraction = fraction, prefix = prefix)


rule downsample_bam:
    input:
        bam = '{prefix}/{sample}/mapping/{sample}_possorted.bam',
    output:
        bam   = temp('{prefix}/{sample}_down_{fraction}/mapping/{sample}_possorted.bam'),
        index = '{prefix}/{sample}_down_{fraction}/mapping/{sample}_possorted.bam.bai',
    params:
        fraction = lambda wildcards: float(wildcards.fraction)/100.0
    threads: 4
    shell:
        """ 
        samtools view  -@ {threads} -s {params.fraction} -b {input} > {output.bam}_temp.bam &&
        samtools sort  -@ {threads} {output.bam}_temp.bam > {output.bam} &&
        samtools index -@ {threads} {output.bam} &&
        rm {output.bam}_temp.bam
        """

rule barcode_metrics_overall:
    input:
        bam = '{prefix}/{sample}_down_{fraction}/mapping/{sample}_possorted.bam',
    output:
        counts = '{prefix}/{sample}_down_{fraction}/metrics/{sample}_barcode_counts.tsv'
    params:
        script    = workflow.basedir + "/scripts/summarise_barcode_counts.py",
        whitelist = workflow.basedir + "/../Resources/737K-cratac-v1.txt.gz"
    threads: 4
    shell:
        'python {params.script} --bam {input.bam} --out {output.counts} --whitelist {params.whitelist}'

rule call_peaks_broad:
    input:
        bam = '{prefix}/{sample}_down_{fraction}/mapping/{sample}_possorted.bam',
        index = '{prefix}/{sample}_down_{fraction}/mapping/{sample}_possorted.bam.bai',
    output:
        peaks_broad = "{prefix}/{sample}_down_{fraction}/peaks/macs2_broad/{sample}_peaks.broadPeak"
    threads: 1
    shell:
        'source ~/.zshrc && conda activate macs2 && '
        'macs2 callpeak -t {input.bam} -f BAM -g hs -n {wildcards.sample} -q 0.01 --broad --outdir {wildcards.prefix}/{wildcards.sample}_down_{wildcards.fraction}/peaks/macs2_broad/ --llocal 100000 --keep-dup 1 --broad-cutoff 0.1 --max-gap 1000 --nomodel --nolambda' 

rule barcode_metrics_peaks:
    input:
        bam = '{prefix}/{sample}_down_{fraction}/mapping/{sample}_possorted.bam',
        peaks = "{prefix}/{sample}_down_{fraction}/peaks/macs2_broad/{sample}_peaks.broadPeak"
    output:
        overlap = "{prefix}/{sample}_down_{fraction}/metrics/{sample}_peak_region_barcode_counts.txt"
    params:
        script    = workflow.basedir + "/scripts/summarise_barcode_counts.py",
        whitelist = workflow.basedir + "/../Resources/737K-cratac-v1.txt.gz"
    threads: 1
    shell:
        'bedtools intersect -abam {input.bam} -b {input.peaks} -u > {output.overlap}_temp && '
        'python {params.script} --bam {output.overlap}_temp --out {output.overlap} --whitelist {params.whitelist} &&'
        'rm {output.overlap}_temp'

rule bam_to_fragments:
    input:
        bam = '{prefix}/{sample}_down_{fraction}/mapping/{sample}_possorted.bam',
        index = '{prefix}/{sample}_down_{fraction}/mapping/{sample}_possorted.bam.bai',
    output:
        fragments = "{prefix}/{sample}_down_{fraction}/mapping/{sample}_fragments.tsv.gz"
    threads: 4
    shell:
        'sinto fragments -b {input.bam} -f {output.fragments}_temp -m 30 -p {threads} && '
        'sort -k1,1 -k2,2n {output.fragments}_temp > {output.fragments}_sorted_temp && '
        'bgzip -@ {threads} -c {output.fragments}_sorted_temp > {output.fragments} && '
        'tabix -p bed {output.fragments} && '
        'rm {output.fragments}_temp {output.fragments}_sorted_temp'

rule cell_picking:
    input:
        metrics_peak = "{prefix}/{sample}_down_{fraction}/metrics/{sample}_peak_region_barcode_counts.txt",
        metrics_all  = '{prefix}/{sample}_down_{fraction}/metrics/{sample}_barcode_counts.tsv',
    output:
        metadata = "{prefix}/{sample}_down_{fraction}/metrics/{sample}_metadata.csv",
        cells    = "{prefix}/{sample}_down_{fraction}/metrics/{sample}_cells.txt",
    threads: 1
    shell:
        'Rscript {workflow.basedir}/scripts/pick_cells_light.R --metrics_all {input.metrics_all} --metrics_peak {input.metrics_peak} --out {output.metadata} --cells_list {output.cells}'

rule create_peak_matrix:
    input:
        peaks = "{prefix}/{sample}_down_{fraction}/peaks/macs2_broad/{sample}_peaks.broadPeak",
        cells = "{prefix}/{sample}_down_{fraction}/metrics/{sample}_cells.txt",
        fragments = "{prefix}/{sample}_down_{fraction}/mapping/{sample}_fragments.tsv.gz"
    output:
        matrix = directory("{prefix}/{sample}_down_{fraction}/matrix_peaks/")
    threads: 2
    shell:
        # Requires f2m in $PATH
        # https://github.com/stuart-lab/f2m
        "f2m --fragments {input.fragments} --bed {input.peaks} --cells {input.cells} --outdir {output.matrix}"
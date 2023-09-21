# CutAndQC performs initial QC on CutAndTag projects
import glob
import os 
from pathlib import Path,PurePath,PurePosixPath
from collections import defaultdict
import pandas as pd
from snakemake.utils import validate, min_version
import plotly as plt
import plotly.graph_objects as go

##### set minimum snakemake version #####
min_version("5.1.2")

include: "src/common.py"

# validate inputs

configfile: "config.yml"
validate(config, schema="schemas/config.schema.yml")

st = pd.read_table('samplesheet.tsv').set_index('sample',drop=False)
validate(st, schema="schemas/samples.schema.yml")

deseq2_md = pd.read_table("src/deseq2_metadata.csv", sep = ",").set_index('sample', drop = False)
validate(deseq2_md, schema="schemas/deseq2.schema.yml")

# parse config files

samps = get_samples()
reads= get_reads()
marks=get_marks()
mark_conditions=get_mark_conditions()

if not os.path.exists("data/fastqc"):
    os.makedirs("data/fastqc")

marks = get_marks()
sample_noigg = [k for k in samps if config["IGG"] not in k]
marks_noigg = [m for m in marks if config["IGG"] not in m]

localrules: frip_plot, fraglength_plot

#singularity: "library://gartician/miniconda-mamba/4.12.0:sha256.7302640e37d37af02dd48c812ddf9c540a7dfdbfc6420468923943651f795591"
singularity: "/home/groups/MaxsonLab/software/singularity-containers/4.12.0_sha256.7302640e37d37af02dd48c812ddf9c540a7dfdbfc6420468923943651f795591.sif"


rule all:
    input:
        expand("data/fastp/{read}.fastq.gz", read=reads) if config["TRIM_ADAPTERS"] else [],
        expand("data/fastqc/{read}_fastqc.{ext}", read=reads, ext = ["html", "zip"]),
        expand("data/fastq_screen/{read}_screen.txt", read=reads),
        expand("data/counts/{mark}_counts.tsv", mark=marks_noigg),
        expand(["data/markd/{sample}.sorted.markd.bam",
                "data/markd/{sample}.sorted.markd.fraglen.tsv",
                "data/tracks/{sample}.bw",
                "data/preseq/lcextrap_{sample}.txt",
                ], sample=samps),
        expand(["data/callpeaks/{sample}_peaks.bed", 
        "data/dtools/fingerprint_{sample}.tsv",
        "data/plotEnrichment/frip_{sample}.tsv",
        ], sample=sample_noigg),
        "data/multiqc/multiqc_report.html",
        expand(["data/deseq2/{mark}/{mark}-rld-pca.pdf",
        "data/deseq2/{mark}/{mark}-vsd-pca.pdf",
        "data/deseq2/{mark}/{mark}-normcounts.csv",
        "data/deseq2/{mark}/{mark}-lognormcounts.csv",
        "data/deseq2/{mark}/{mark}-rld.pdf",
        "data/deseq2/{mark}/{mark}-vsd.pdf",
        "data/deseq2/{mark}/{mark}-vsd-dist.pdf",
        "data/deseq2/{mark}/{mark}-rld-dist.pdf",
        "data/deseq2/{mark}/{mark}-dds.rds",
        "data/homer/{mark}.done"], mark=marks_noigg),
        # quality control plots
        "data/markd/fraglen.html",
        "data/plotEnrichment/frip.html",
        expand("data/mergebw/{mark_condition}.bw", mark_condition=mark_conditions),
        expand("data/highConf/{mark_condition}.highConf.bed", mark_condition=mark_conditions)


if config["TRIM_ADAPTERS"]:
    # trim adapters from reads before alignment
    rule fastp:
        input:
            r1 = "data/raw/{sample}_R1.fastq.gz",
            r2 = "data/raw/{sample}_R2.fastq.gz"
        output:
            r1 = "data/fastp/{sample}_R1.fastq.gz",
            r2 = "data/fastp/{sample}_R2.fastq.gz"
        conda: 
            "envs/fastp.yml"
        log:
            "data/logs/fastp/{sample}.fastp.json"
        threads: 8
        shell:
            "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} --detect_adapter_for_pe --thread {threads} -j {log} -h /dev/null"

# fastqc for each read 
rule fastqc:
    input:
        "data/fastp/{read}.fastq.gz" if config["TRIM_ADAPTERS"] else "data/raw/{read}.fastq.gz"
    output:
        html="data/fastqc/{read}_fastqc.html",
        zip="data/fastqc/{read}_fastqc.zip"
    conda:
        "envs/fastqc.yml"
    log:
        "data/logs/fastqc_{read}.log"
    threads: 4
    shell:
        "fastqc -t {threads} --outdir data/fastqc {input} > {log} 2>&1"

# detect contaminants
rule fastq_screen:
    input:
        "data/raw/{read}.fastq.gz"
    output:
        "data/fastq_screen/{read}_screen.txt"
    conda:
        "envs/fastq_screen.yml"
    log:
        "data/logs/fastq_screen_{read}.log"
    threads: 8
    shell:
        "fastq_screen --aligner bowtie2 --threads {threads} --outdir data/fastq_screen "
        "--conf {config[FASTQ_SCREEN]} --force {input} > {log} 2>&1"

# align samples to genome
rule bowtie2:
    input:
        r1 = "data/fastp/{sample}_R1.fastq.gz" if config["TRIM_ADAPTERS"] else "data/raw/{sample}_R1.fastq.gz",
        r2 = "data/fastp/{sample}_R2.fastq.gz" if config["TRIM_ADAPTERS"] else "data/raw/{sample}_R2.fastq.gz"
    output:
        "data/aligned/{sample}.bam"
    log:
        err="data/logs/bowtie2_{sample}.err"
    conda:
        "envs/align.yml"
    threads: 8
    shell:
        "bowtie2 --local --very-sensitive-local "
        "--no-unal --no-mixed --threads {threads} "
        "--no-discordant --phred33 "
        "-I 10 -X 700 -x {config[GENOME]} "
        "-1 {input.r1} -2 {input.r2} 2>{log.err} | samtools view -@ {threads} -Sbh - > {output}"

rule sort:
    input:
        "data/aligned/{sample}.bam"
    output: 
        temp("data/aligned/{sample}.sort.bam")
    conda:
        "envs/sambamba.yml"
    threads: 4
    log:
        "data/logs/sambamba_sort_{sample}.log"
    shell:
        "sambamba sort --tmpdir=data/aligned -t {threads} -o {output} {input} > {log} 2>&1"

rule markdup:
    input:
        rules.sort.output
    output:
        "data/markd/{sample}.sorted.markd.bam"
    conda:
        "envs/sambamba.yml"
    threads: 4
    log:
        "data/logs/sambamba_markdup_{sample}.log"
    shell:
        "sambamba markdup --tmpdir=data/markd -t {threads} {input} {output} > {log} 2>&1"

rule index:
    input:
        rules.markdup.output
    output:
        "data/markd/{sample}.sorted.markd.bam.bai"
    conda:
        "envs/sambamba.yml"
    threads: 4
    log:
        "data/logs/samtools_index_{sample}.log"
    shell:
        "sambamba index -t {threads} {input} > {log} 2>&1"

rule tracks:
    input:
        bam = rules.markdup.output,
        bai = rules.index.output,
    output:
        "data/tracks/{sample}.bw"
    conda:
        "envs/dtools.yml"
    threads: 8
    shell:
        "bamCoverage -b {input[0]} -o {output} --binSize 10 --smoothLength 50 --normalizeUsing CPM -p {threads} "

rule merge_bw:
    input:
        get_tracks_by_mark_condition
    output:
        "data/mergebw/{mark_condition}.bw"
    conda:
        "envs/mergebw.yml"
    resources:
        mem_mb = 8000
    shell:
        "bash src/mergebw.sh -c {config[CSIZES]} -o {output} {input}"

# merge all peaks to get union peak with at least
# two reps per condition per peak
rule make_high_conf_peaks:
    input:
        get_peaks_by_mark_condition
    output:
        "data/highConf/{mark_condition}.highConf.bed"
    conda:
        "envs/bedtools.yml"
    shell:
        "cat {input} | sort -k1,1 -k2,2n | "
        "bedtools merge | "
        "bedtools intersect -a - -b {input} -c | "
        "awk -v OFS='\t' '$4>=2 {{print}}' > {output}"


rule fraglength:
    input:
        rules.markdup.output
    output:
        "data/markd/{sample}.sorted.markd.fraglen.tsv"
    conda:
        "envs/align.yml"
    shell:
        "src/fraglen-dist.sh {input} {output}"

rule fraglength_plot:
    input:
        expand("data/markd/{sample}.sorted.markd.fraglen.tsv", sample = samps)
    output:
        "data/markd/fraglen.html"
    run:
        pd.options.plotting.backend = "plotly"
        dfs = []
        for i in input:
            cond_marker = [os.path.basename(i).split(".")[0]]
            temp_df = pd.read_csv(i, sep = "\t", index_col = 0, names = cond_marker)
            dfs.append(temp_df)
        df = pd.concat(dfs, axis = 1)
        fraglen = df.plot()
        fraglen.update_layout( 
            title='Fragment Length Distribution', 
            xaxis_title='Fragment Length (bp)', 
            yaxis_title='Counts', 
            legend_title_text='Samples')
        fraglen.write_html(str(output))

rule preseq:
    input:
       rules.markdup.output
    output:
        "data/preseq/estimates_{sample}.txt"
    resources:
        defect_mode = defect_mode
    conda:
        "envs/preseq.yml"
    log:
        "data/logs/preseq_{sample}.log"
    shell:
        "preseq c_curve -B {resources.defect_mode} -l 1000000000 -o {output} {input} > {log} 2>&1"

rule preseq_lcextrap:
    input:
        rules.markdup.output
    output:
        "data/preseq/lcextrap_{sample}.txt"
    resources:
        defect_mode = defect_mode
    conda:
        "envs/preseq.yml"
    log:
        "data/logs/preseq_{sample}.log"
    shell:
        "preseq lc_extrap -B {resources.defect_mode} -l 1000000000 -e 1000000000 -o {output} {input} > {log} 2>&1"

rule plotFinger:
    input:
        "data/markd/{sample}.sorted.markd.bam", "data/markd/{sample}.sorted.markd.bam.bai"
    output:
        "data/dtools/fingerprint_{sample}.tsv"
    conda:
        "envs/dtools.yml"
    log:
        "data/logs/fingerprint_{sample}.log"
    shell:
        "plotFingerprint -b {input[0]} --smartLabels --outRawCounts {output} > {log} 2>&1"

rule callpeaks:
    input:
        get_callpeaks
    output:
        "data/callpeaks/{sample}_peaks.bed"
    conda: 
        "envs/gopeaks.yml"
    log:
        "data/callpeaks/{sample}_gopeaks.json"
    params:
        igg = get_igg,
        params = callpeaks_params
    shell:
        "gopeaks -b {input[0]} {params.igg} -o data/callpeaks/{wildcards.sample} {params.params} > {log} 2>&1"

# get consensus
rule consensus:
    input:
       expand("data/callpeaks/{sample}_peaks.bed", sample=sample_noigg)
    output:
       "data/counts/{mark}_counts.tsv"
    conda:
       "envs/bedtools.yml"
    shell:
       "bash src/consensus_peaks.sh {wildcards.mark} {config[N_INTERSECTS]} {output}"

rule frip:
    input:
        rules.callpeaks.output, "data/markd/{sample}.sorted.markd.bam"
    output:
        "data/plotEnrichment/frip_{sample}.png", "data/plotEnrichment/frip_{sample}.tsv"
    conda:
        "envs/dtools.yml"
    log:
        "data/logs/plotEnrichment_{sample}.log"
    shell:
        "bash src/skip_frip.sh {input[0]} {input[1]} {output[0]} {output[1]} {log}"

rule frip_plot:
    input:
        expand("data/plotEnrichment/frip_{sample}.tsv", sample = sample_noigg)
    output:
        "data/plotEnrichment/frip.html"
    run:
        pd.options.plotting.backend = "plotly"
        dfs = []
        for i in sorted(input):
            cond_marker = "_".join(i.split("_")[1:3])
            temp_df = pd.read_csv(i, sep = "\t", usecols=["percent"]).rename(columns = {'percent': cond_marker})
            dfs.append(temp_df)
        frip_df = pd.concat(dfs, axis = 1)
        frip_df = frip_df.rename(index={0: 'inside'})
        frip_df.loc["outside"] = 100 - frip_df.loc['inside']
        fig = go.Figure(data=[
            go.Bar(name="inside_peaks", x=frip_df.columns, y=frip_df.loc['inside'], marker_color='rgb(255, 201, 57)'),
            go.Bar(name='outside_peaks', x=frip_df.columns, y=frip_df.loc['outside'], marker_color='rgb(0,39, 118)')
        ])
        fig.update_layout(barmode='stack', 
            title='Fraction of Reads in Peaks by Sample', 
            xaxis_tickfont_size=14, yaxis=dict(title='Fraction of reads in peaks', 
            titlefont_size=16, tickfont_size=14), xaxis=dict(title='Samples'))
        fig.write_html(str(output))

rule deseq2:
    input:
        counts="data/counts/{mark}_counts.tsv",
        meta="src/deseq2_metadata.csv",
        genes=config["GENES"]
    output:
        pcaPlot="data/deseq2/{mark}/{mark}-rld-pca.pdf",
        pcaPlotVsd="data/deseq2/{mark}/{mark}-vsd-pca.pdf",
        normCounts="data/deseq2/{mark}/{mark}-normcounts.csv",
        lnormCounts="data/deseq2/{mark}/{mark}-lognormcounts.csv",
        sdMeanRld="data/deseq2/{mark}/{mark}-rld.pdf",
        sdMeanVsd="data/deseq2/{mark}/{mark}-vsd.pdf",
        sampleDistVsd="data/deseq2/{mark}/{mark}-vsd-dist.pdf",
        sampleDistRld="data/deseq2/{mark}/{mark}-rld-dist.pdf",
        rds="data/deseq2/{mark}/{mark}-dds.rds"
    params:
        mark=lambda wildcards: wildcards.mark,
        outdir = "data/deseq2/",
        numclusters = 4
    conda:
        "envs/deseq2.yml"
    script:
        "src/deseq2.R"

rule homer:
    input:
        "data/deseq2/{mark}/{mark}-dds.rds"
    output:
        "data/homer/{mark}.done"
    conda:
        "envs/homer.yml"
    shell:
        "bash src/homer.sh -m {wildcards.mark} -s 1 -p 8 -g {config[FASTA]}"
# this rule submits HOMER runs to SLURM if -s = 1. A run is each unique contrast

rule multiqc:
    input:
        expand("data/fastqc/{read}_fastqc.zip", read=reads),
        expand("data/fastq_screen/{read}_screen.txt", read=reads),
        expand("data/plotEnrichment/frip_{sample}.tsv", sample=sample_noigg),
        # comment out the expand("data/deseq2/...) line if no replicates used in experiment
        expand("data/deseq2/{mark}/{mark}-dds.rds",mark=marks_noigg),
        expand("data/preseq/lcextrap_{sample}.txt", sample=samps)
    output:
        "data/multiqc/multiqc_report.html"
    conda:
        "envs/multiqc.yml"
    log:
        "data/logs/multiqc.log"
    shell:
        # comment out the "export ..." line if not running pipeline through Singularity
        "export LC_ALL=C.UTF-8; export LANG=C.UTF-8; "
        "multiqc data/ -f -c src/multiqc_conf.yml -o data/multiqc --ignore data/homer > {log} 2>&1"

# export different locales for singularity workaround: https://click.palletsprojects.com/en/8.1.x/unicode-support/


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
configfile: "src/config.yml"

st = pd.read_table('samplesheet.tsv').set_index('sample',drop=False)
validate(st, schema="schemas/samples.schema.yml")

samps = get_samples()
reads= get_reads()
marks=get_marks()
mark_conditions=get_mark_conditions()

marks = get_marks()
sample_noigg = [k for k in samps if config["IGG"] not in k]
marks_noigg = [m for m in marks if config["IGG"] not in m]

fastqScreenDict = {
'database': {
   'hg38': {
     'bowtie2': config["BOWTIE2"]["HG38"][0]},
   'mm10': {
     'bowtie2': config["BOWTIE2"]["MM10"][0]}, 
   'ecoli': {
     'bowtie2': config["BOWTIE2"]["ECOLI"][0]}, 
   'myco': {
     'bowtie2': config["BOWTIE2"]["MYCO"][0]}, 
 },
 'aligner_paths': {'bowtie2': 'bowtie2'}
}

localrules: frip_plot, fraglength_plot

rule all:
    input:
        expand("data/fastqc/{read}.html", read=reads),
        expand("data/fastq_screen/{read}.fastq_screen.txt", read=reads),
        expand("data/counts/{mark}_counts.tsv", mark=marks_noigg),
        expand(["data/markd/{sample}.sorted.markd.bam",
                "data/markd/{sample}.sorted.markd.bam",
                "data/markd/{sample}.sorted.markd.fraglen.tsv",
                "data/tracks/{sample}.bw",
                ], sample=samps),
        expand(["data/callpeaks/{sample}_peaks.bed", 
        "data/preseq/lcextrap_{sample}.txt",
        "data/dtools/fingerprint_{sample}.tsv",
        "data/plotEnrichment/frip_{sample}.tsv",
        ], sample=samps),
        "data/multiqc/multiqc_report.html",
        expand(["data/deseq2/{mark}/{mark}-rld-pca.png",
        "data/deseq2/{mark}/{mark}-vsd-pca.png",
        "data/deseq2/{mark}/{mark}-normcounts.csv",
        "data/deseq2/{mark}/{mark}-lognormcounts.csv",
        "data/deseq2/{mark}/{mark}-rld.png",
        "data/deseq2/{mark}/{mark}-vsd.png",
        "data/deseq2/{mark}/{mark}-vsd-dist.png",
        "data/deseq2/{mark}/{mark}-rld-dist.png",
        "data/deseq2/{mark}/{mark}-dds.rds"], mark=marks_noigg),
        # quality control plots
        "data/markd/fraglen.html",
        "data/plotEnrichment/frip.html",
        expand("data/mergebw/{mark_condition}.bw", mark_condition=mark_conditions)

# fastqc for each read 
rule fastqc:
    input:
        "data/raw/{read}.fastq.gz"
    output:
        html="data/fastqc/{read}.html",
        zip="data/fastqc/{read}_fastqc.zip"
    log:
        "data/logs/fastqc_{read}.log"
    threads: 4
    wrapper:
        "0.65.0/bio/fastqc"

# detect contaminants
rule fastq_screen:
    input:
        "data/raw/{read}.fastq.gz"
    output:
        txt="data/fastq_screen/{read}.fastq_screen.txt",
        png="data/fastq_screen/{read}.fastq_screen.png"
    params:
        fastq_screen_config=fastqScreenDict,
        subset=100000,
        aligner='bowtie2'
    log:
        "data/logs/fastq_screen_{read}.log"
    threads: 8
    wrapper:
        "0.65.0/bio/fastq_screen"

# align samples to genome
rule bowtie2:
    input:
        get_bowtie2_input
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
        "-1 {input[0]} -2 {input[1]} 2>{log.err} | samtools view -@ {threads} -Sbh - > {output}"

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
        rules.markdup.output
    output:
        "data/tracks/{sample}.bw"
    conda:
        "envs/dtools.yml"
    threads:
        8
    shell:
        "bamCoverage -p {threads} --binSize 10 --smoothLength 50 --normalizeUsing CPM -b {input} -o {output}"

rule merge_bw:
    input:
        get_tracks_by_mark_condition
    output:
        "data/mergebw/{mark_condition}.bw"
    conda:
        "envs/mergebw.yml"
    shell:
        "bash src/mergebw.sh -c {config[CSIZES]} -o {output} {input}"
    
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
    conda:
        "envs/preseq.yml"
    log:
        "data/logs/preseq_{sample}.log"
    shell:
        "preseq c_curve -B -P -o {output} {input} > {log} 2>&1" 

rule preseq_lcextrap:
    input:
        rules.markdup.output
    output:
        "data/preseq/lcextrap_{sample}.txt"
    conda:
        "envs/preseq.yml"
    log:
        "data/logs/preseq_{sample}.log"
    shell:
        "preseq lc_extrap -B -P -e 1000000000 -o {output} {input} > {log} 2>&1"
    

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
        "data/logs/callpeaks_{sample}.log"
    params:
        igg=get_igg
    shell:
        """
        gopeaks -bam {input[0]} {params.igg} -of {output} > {log} 2>&1
        """

# get consensus
rule consensus:
    input:
       expand("data/callpeaks/{sample}_peaks.bed", sample=sample_noigg)
    output:
       "data/counts/{mark}_counts.tsv"
    conda:
       "envs/bedtools.yml"
    log:
       "data/logs/{mark}_counts.log"
    shell:
       "src/consensus_peaks.sh {wildcards.mark} data/callpeaks data/markd {output} > {log} 2>&1"

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
        "plotEnrichment -b {input[1]} --BED {input[0]} --regionLabels 'frip' --outRawCounts {output[1]} -o {output[0]} > {log} 2>&1"

rule frip_plot:
    input:
        expand("data/plotEnrichment/frip_{sample}.tsv", sample = samps)
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
        pcaPlot="data/deseq2/{mark}/{mark}-rld-pca.png",
        pcaPlotVsd="data/deseq2/{mark}/{mark}-vsd-pca.png",
        normCounts="data/deseq2/{mark}/{mark}-normcounts.csv",
        lnormCounts="data/deseq2/{mark}/{mark}-lognormcounts.csv",
        sdMeanRld="data/deseq2/{mark}/{mark}-rld.png",
        sdMeanVsd="data/deseq2/{mark}/{mark}-vsd.png",
        sampleDistVsd="data/deseq2/{mark}/{mark}-vsd-dist.png",
        sampleDistRld="data/deseq2/{mark}/{mark}-rld-dist.png",
        rds="data/deseq2/{mark}/{mark}-dds.rds"
    params:
        mark=lambda wildcards: wildcards.mark,
        outdir = "data/deseq2/",
        numclusters = 4
    conda:
        "envs/deseq2.yml"
    script:
        "src/deseq2.R"

rule multiqc:
    input:
        expand("data/plotEnrichment/frip_{sample}.tsv", sample=samps),
        expand("data/deseq2/{mark}/{mark}-dds.rds",mark=marks_noigg),
        directory("data/")
    output:
        "data/multiqc/multiqc_report.html"
    conda:
        "envs/multiqc.yml"
    log:
        "data/logs/multiqc.log"
    shell:
        "multiqc -f -c src/multiqc_conf.yml -o data/multiqc {input} > {log} 2>&1"


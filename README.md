# cutTag-pipeline

Snakemake pipeline for Cut&amp;Tag analysis 

# Setup

1. Configure the project directory

```
#clone to your local called 'my-project'
git clone https://github.com/maxsonBraunLab/cutTag-pipeline.git my-project

#create directory for your fastq files
cd my-project
mkdir -p data/raw

# link fastqs to data/raw 
ln -s /path/to/fastq/files/* data/raw

# make scripts executable
chmod +x src/*.py src/*.sh
```


Rename all samples in data/raw to simplify sample information, for example:

This file:
LIB200706TB_M6Q3_RBP1_S93_L001_R1_001.fastq.gz

Is renamed to:
M6Q_3_RBP1_R1.fastq.gz

2. Make the sample sheet and deseq2 metadata.

Activate an environment containing snakemake, and then run the script `make_sample_sheet.py` script from the root of the directory.

```
$ python src/make_sample_sheet.py data/raw
```

This will make a samplesheet for the experiment called samplesheet.tsv in the root of the directory as well as the file `src/deseq2_metadata.csv`, the contents of the samplesheet will be structured like the following example:

```
sample	R1	R2	mark	condition	igg
HoxE_1_IgG	data/raw/HoxE_1_IgG_R1.fastq.gz	data/raw/HoxE_1_IgG_R2.fastq.gz	IgG	HoxE	HoxE_1_IgG
HoxE_1_Rbp1	data/raw/HoxE_1_Rbp1_R1.fastq.gz	data/raw/HoxE_1_Rbp1_R2.fastq.gz	Rbp1	HoxE	HoxE_1_Rbp1
HoxE_2_Rbp1	data/raw/HoxE_2_Rbp1_R1.fastq.gz	data/raw/HoxE_2_Rbp1_R2.fastq.gz	Rbp1	HoxE	HoxE_2_Rbp1
HoxE_3_Rbp1	data/raw/HoxE_3_Rbp1_R1.fastq.gz	data/raw/HoxE_3_Rbp1_R2.fastq.gz	Rbp1	HoxE	HoxE_3_Rbp1
HoxM_1_IgG	data/raw/HoxM_1_IgG_R1.fastq.gz	data/raw/HoxM_1_IgG_R2.fastq.gz	IgG	HoxM	HoxM_1_IgG
HoxM_1_Rbp1	data/raw/HoxM_1_Rbp1_R1.fastq.gz	data/raw/HoxM_1_Rbp1_R2.fastq.gz	Rbp1	HoxM	HoxM_1_Rbp1
HoxM_2_Rbp1	data/raw/HoxM_2_Rbp1_R1.fastq.gz	data/raw/HoxM_2_Rbp1_R2.fastq.gz	Rbp1	HoxM	HoxM_2_Rbp1
HoxM_3_Rbp1	data/raw/HoxM_3_Rbp1_R1.fastq.gz	data/raw/HoxM_3_Rbp1_R2.fastq.gz	Rbp1	HoxM	HoxM_3_Rbp1
HoxW_1_IgG	data/raw/HoxW_1_IgG_R1.fastq.gz	data/raw/HoxW_1_IgG_R2.fastq.gz	IgG	HoxW	HoxW_1_IgG
HoxW_1_Rbp1	data/raw/HoxW_1_Rbp1_R1.fastq.gz	data/raw/HoxW_1_Rbp1_R2.fastq.gz	Rbp1	HoxW	HoxW_1_Rbp1
HoxW_2_Rbp1	data/raw/HoxW_2_Rbp1_R1.fastq.gz	data/raw/HoxW_2_Rbp1_R2.fastq.gz	Rbp1	HoxW	HoxW_2_Rbp1
HoxW_3_Rbp1	data/raw/HoxW_3_Rbp1_R1.fastq.gz	data/raw/HoxW_3_Rbp1_R2.fastq.gz	Rbp1	HoxW	HoxW_3_Rbp1
```

The script splits the file name on the '_' and uses the first split for the condition, and the second split for the mark. The 'igg' column is the same as the 'sample' column and should be manually replaced with the sample name of the IGG or control you would like to use for that sample. If the sample is and IgG it can be the same as it's name, and won't affect peak calling. 

So a fixed version of the table above would look like this:


```
sample	R1	R2	mark	condition	igg
HoxE_1_IgG	data/raw/HoxE_1_IgG_R1.fastq.gz	data/raw/HoxE_1_IgG_R2.fastq.gz	IgG	HoxE	HoxE_1_IgG
HoxE_1_Rbp1	data/raw/HoxE_1_Rbp1_R1.fastq.gz	data/raw/HoxE_1_Rbp1_R2.fastq.gz	Rbp1	HoxE	HoxE_1_IgG
HoxE_2_Rbp1	data/raw/HoxE_2_Rbp1_R1.fastq.gz	data/raw/HoxE_2_Rbp1_R2.fastq.gz	Rbp1	HoxE	HoxE_1_IgG
HoxE_3_Rbp1	data/raw/HoxE_3_Rbp1_R1.fastq.gz	data/raw/HoxE_3_Rbp1_R2.fastq.gz	Rbp1	HoxE	HoxE_1_IgG
HoxM_1_IgG	data/raw/HoxM_1_IgG_R1.fastq.gz	data/raw/HoxM_1_IgG_R2.fastq.gz	IgG	HoxM	HoxM_1_IgG
HoxM_1_Rbp1	data/raw/HoxM_1_Rbp1_R1.fastq.gz	data/raw/HoxM_1_Rbp1_R2.fastq.gz	Rbp1	HoxM	HoxM_1_IgG
HoxM_2_Rbp1	data/raw/HoxM_2_Rbp1_R1.fastq.gz	data/raw/HoxM_2_Rbp1_R2.fastq.gz	Rbp1	HoxM	HoxM_1_IgG
HoxM_3_Rbp1	data/raw/HoxM_3_Rbp1_R1.fastq.gz	data/raw/HoxM_3_Rbp1_R2.fastq.gz	Rbp1	HoxM	HoxM_1_IgG
HoxW_1_IgG	data/raw/HoxW_1_IgG_R1.fastq.gz	data/raw/HoxW_1_IgG_R2.fastq.gz	IgG	HoxW	HoxW_1_IgG
HoxW_1_Rbp1	data/raw/HoxW_1_Rbp1_R1.fastq.gz	data/raw/HoxW_1_Rbp1_R2.fastq.gz	Rbp1	HoxW	HoxW_1_IgG
HoxW_2_Rbp1	data/raw/HoxW_2_Rbp1_R1.fastq.gz	data/raw/HoxW_2_Rbp1_R2.fastq.gz	Rbp1	HoxW	HoxW_1_IgG
HoxW_3_Rbp1	data/raw/HoxW_3_Rbp1_R1.fastq.gz	data/raw/HoxW_3_Rbp1_R2.fastq.gz	Rbp1	HoxW	HoxW_1_IgG
```

For this example there was only one IgG per condition, so the sample name corresponding to that IGG was used for each sample in the condition. In the case that each sample had it's own control file, each entry would correspond to the IgG for that sample. If only one IgG was used in the whole experiment, then it's sample name could be used for each row. If you are not using IgG set config['USEIGG'] to false, and don't modify the samplesheet.


3. Edit configuration files 

Edit runtime configuration in the file src/config.yml:

- Specify the path to the bowtie2 index for the genome you are aligning to.
- Specify the path to the bowtie2 indices for the fastq_screen genomes.

The pipeline detects samples in the subdirectory data/raw with the following assumptions:

 - Paired end reads
 - Read 1 and 2 are designated using "_R1", and "_R2"
 - the epigenetic mark label is the second split of the sample name by _ delimeter. For example M6C3_4Me1_S12_R2.fastq.gz will have the {mark} wildcard set to 4me1. This effects the output files from the calculation of counts tables.

Make sure deseq2_metadata.csv looks right. The file is created when you run `make_sample_sheet.py` and should have the following properties:

 - should have two columns labeled "sample", and "condition"
 - the sample column corresponds to the replicates of the given condition, and should be the same as the first split of the raw file: e.g. M6C3_4Me1_S12_R2.fastq.gz will have "sample" equal to M6C3.
 - the condition should be the name for each sample condition, and doees not have to come from the file name, it could be changed to whatever you would like to have displayed in the deseq2 plots.

 The file src/deseq2_metadata is populated with the following example data:

```
sample,condition
HoxE_1_IgG,HoxE
HoxE_1_Rbp1,HoxE
HoxE_2_Rbp1,HoxE
HoxE_3_Rbp1,HoxE
HoxM_1_IgG,HoxM
HoxM_1_Rbp1,HoxM
HoxM_2_Rbp1,HoxM
HoxM_3_Rbp1,HoxM
HoxW_1_IgG,HoxW
HoxW_1_Rbp1,HoxW
HoxW_2_Rbp1,HoxW
HoxW_3_Rbp1,HoxW
```

# Method

![](rulegraph.svg)

# Execution

Setup snakemake profile to run on compute cluster:

SLURM: follow [these instructions](https://github.com/Snakemake-Profiles/slurm)

Example of how to run pipeline

`
snakemake --use-conda --profile slurm -j 60 --cluster-config src/cluster.yml
`

A "dry-run" can be accomplished to see what files would be generated by using the command:

`
snakemake -nrp
`

If the output of this command is not green, there is something wrong.


# Output


Below is an explanation of each output directory:

```
aligned - sorted and aligned sample bam files
callpeaks - the output of callpeaks.py with peak files and normalized bigwigs for each sample.
counts - the raw counts for each mark over a consensus peak list
deseq2 - the results of deseq2 fore each counts table (by mark)
dtools - fingerprint plot data for multiqc to use
fastqc - fastqc results
fastq_screen - fastq_screen results
logs - runtime logs for each snakemake rule
markd - duplicate marked bam files
multiqc - contains the file multiqc_report.html with a lot of QC info about the experiment.
plotEnrichment - FRIP statistics for each mark
preseq - library complexity data for multiqc report
```

## Deseq2 outputs

Each mark should have the following output files:

```
"data/deseq2/{mark}/{mark}-rld-pca.png" - PCA of counts after reguarlized log transformation. 
"data/deseq2/{mark}/{mark}-vsd-pca.png" - PCA of counts after variance stabilizing transformation.
"data/deseq2/{mark}/{mark}-normcounts.csv" - normalized count for each sample in each consensus peak.
"data/deseq2/{mark}/{mark}-lognormcounts.csv" - log2 normalized counts for each sample in each consensus peak.
"data/deseq2/{mark}/{mark}-rld.png", - the sdVsMean plot using regularized log transformation.
"data/deseq2/{mark}/{mark}-vsd.png" - the sdVsMean plot using variance stabilizing transformation.
"data/deseq2/{mark}/{mark}-vsd-dist.png" - the sample distance matrix after variance stabilizing transformation.
"data/deseq2/{mark}/{mark}-rld-dist.png" - the sample distance matrix using regularized log transformation.
"data/deseq2/{mark}/{mark}-dds.rds" - the R object with the results of running the DEseq() function.
```

For each contrast, the differentially expressed genes are written to a file ending in `-diffexp.tsv` as well as those with an adjusted p-value less than 0.05 with the extension `-sig05-diffexp.tsv`. A summary of the results using an alpha of 0.05 is also written to a file with the extension `-sig05-diffexp-summary.txt`. Additionally two MA plots are written to the file ending in `plotMA.png` that have highlighted differential peaks with an adjusted p-value less than 0.1.

See the following paper for further explanations of the above plots and data transforms:
https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#exporting-results-to-csv-files 

# cutTag-pipeline

Snakemake pipeline for Cut&amp;Tag analysis 

# Setup

Clone this repository into your project directory:

```
#clone to your local called 'my-project'
git clone git@github.com:maxsonBraunLab/cutAndQC.git my-project

#create directory for your fastq files
cd my-project
mkdir -p data/raw

# link your fastqs to here
ln -s /path/to/fastq/files/* data/raw
```


Rename all samples in data/raw to simplify sample information, for example:

This file:
LIB200706TB_M6Q3_RBP1_S93_L001_R1_001.fastq.gz

Is renamed to:
M6Q3_RBP1_S93_R1.fastq.gz


Edit runtime configuration in the file src/config.yml:

- Specify the path to the bowtie2 index for the genome you are aligning to.
- Specify the path to the bowtie2 indices for the fastq_screen genomes.

The pipeline detects samples in the subdirectory data/raw with the following assumptions:

 - Paired end reads
 - Read 1 and 2 are designated using "_R1", and "_R2"


# Execution

Setup snakemake profile to run on compute cluster:

SLURM: follow [these instructions](https://github.com/Snakemake-Profiles/slurm)

Example of how to run pipeline

```
snakemake --use-conda --profile slurm -j 60 --latency-wait 60
```

To change runtime parameters for indvidual rules, you can provide a cluster configuration file via the --cluster-       config flag, for example:

```
snakemake --use-conda --profile slurm --cluster-config src/cluster.yml -j 60 --latency-wait 60
```


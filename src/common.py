# common holds pyhton function to be used in the snakefile
import pandas as pd

# map samples to fastqs
def get_samples():
    """
    return list of samples from samplesheet.tsv
    """
    return list(st.index)

def get_marks():
    """
    return list of marks from samplesheet.tsv
    """
    return list(st['mark'])

def get_bowtie2_input(wildcards):
    """
    returns reads associated with a sample
    """
    r1=st.loc[wildcards.sample]['R1']
    r2=st.loc[wildcards.sample]['R2']
    return r1,r2

def get_reads():
    """
    get list of all reads
    """
    rlist=list(st['R1'])+list(st['R2'])
    rlist=[os.path.basename(f).split('.')[0] for f in rlist]
    return rlist

def get_igg(wildcards):
    """
    Returns the igg file for the sample unless
    the sample is IgG then no control file is used.
    """ 
    if config['USEIGG']:
        igg=st.loc[wildcards.sample]['igg']
        iggbam=f'data/markd/{igg}.sorted.markd.bam'
        isigg=config['IGG'] in wildcards.sample
        if not isigg:
            return f'-cf {iggbam}'
        else:
            return ""
    else:
        return ""

def get_callpeaks(wildcards):
    """
    Returns the callpeaks input files
    """
    bam=f"data/markd/{wildcards.sample}.sorted.markd.bam"
    bai=f"data/markd/{wildcards.sample}.sorted.markd.bam.bai"
    cp="src/gopeaks"
    if config["USEIGG"]:
        igg=st.loc[wildcards.sample]['igg']
        iggbam=f'data/markd/{igg}.sorted.markd.bam'
        iggbam=f'data/markd/{igg}.sorted.markd.bam.bai'
        isigg=config['IGG'] in wildcards.sample
        if not isigg:
            return [bam,bai,cp,iggbam]
        else:
            return [bam,bai,cp]
    else:
        return [bam,bai,cp]

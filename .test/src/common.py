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

def get_mark_conditions():
    """
    return list of samples by condition
    """
    st['mark_condition']=st['mark'].astype(str)+"_"+st['condition']
    return st['mark_condition'].unique().tolist()

def get_tracks_by_mark_condition(wildcards):
    """
    return list of tracks by mark_condition
    """
    st['mark_condition']=st['mark'].astype(str)+"_"+st['condition']
    samps=st.groupby(["mark_condition"])['sample'].apply(list)[wildcards.mark_condition]
    return [f"data/tracks/{s}.bw" for s in samps]
    
def get_peaks_by_mark_condition(wildcards):
    """
    return list of peaks by mark_condition
    """
    st['mark_condition']=st['mark'].astype(str)+"_"+st['condition']
    samps=st.groupby(["mark_condition"])['sample'].apply(list)[wildcards.mark_condition]
    return [f"data/callpeaks/{s}_peaks.bed" for s in samps]

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
            return f'-control {iggbam}'
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
    if config["USEIGG"]:
        igg=st.loc[wildcards.sample]['igg']
        iggbam=f'data/markd/{igg}.sorted.markd.bam'
        iggbam=f'data/markd/{igg}.sorted.markd.bam.bai'
        isigg=config['IGG'] in wildcards.sample
        if not isigg:
            return [bam,bai,iggbam]
        else:
            return [bam,bai]
    else:
        return [bam,bai]

#!/usr/bin/env python3

import pandas as pd
import sys
import os
import re
import yaml

def main():

    # get list of fastq files in raw dir
    raw = sys.argv[1]
    outfile = "samplesheet.tsv"

    samps = [os.path.join(raw,f) for f in os.listdir(raw)]
    regex = re.compile(".+\.fastq\.gz")
    sfilt = [s for s in samps if regex.search(s)]

    # make df 
    df=pd.DataFrame(sfilt, columns=["file"])
    df["sample"] = df["file"].apply(lambda x: "_".join(os.path.basename(x).split("_")[0:3]))
    df["mark"] = df["file"].apply(lambda x: x.split("_")[2])
    df["condition"] = df["sample"].apply(lambda x: x.split("_")[0])

    df["read"] = df["file"].apply(lambda x: "R1" if "R1" in x else "R2")
    df["igg"] = df["sample"]

    df=df.pivot_table(index=['sample','igg','mark','condition'], columns='read', values='file',aggfunc=sum).reset_index()
    df=df[["sample","R1","R2","mark","condition","igg"]]
    df.columns.name = None
    df["gopeaks"] = "-"

    # write sample sheet to file
    df.to_csv(outfile, sep='\t', index=False)

    # define how IgG is spelled
    with open("config.yml", 'r') as f:
        config = yaml.safe_load(f)

    igg = config["IGG"]

    # write deseq2 metadata to file
    dmet=df.copy()
    dmet=dmet[['sample','condition']]
    dmet=dmet[- dmet['sample'].str.contains(igg)]
    dmet.drop_duplicates(inplace=True)
    dmet.to_csv('src/deseq2_metadata.csv', index=False)

    # write diffbind metadata to file
    diffbind_outfile = "src/diffbind_config.csv"

    diffbind = df[['sample','mark','condition']]
    diffbind=diffbind[- diffbind['sample'].str.contains(igg)]
    
    diffbind["replicate"] = diffbind["sample"].apply(lambda x: x.split("_")[1])
    diffbind["bam"] = diffbind["sample"].apply(lambda x: "data/markd/" + x + ".sorted.markd.bam")
    diffbind["peaks"] = diffbind["mark"].apply(lambda x: "data/counts/" + x + "_consensus.bed")

    diffbind["peakcaller"] = ["raw"] * diffbind.shape[0]

    diffbind["peakformat"] = ["bed"] * diffbind.shape[0]
    diffbind.to_csv(diffbind_outfile, sep = ",", header = ["SampleID", "Factor", "Condition", "Replicate", "bamReads", "Peaks", "PeakCaller", "PeakFormat"], index = False)



if __name__ == "__main__": 
    main()

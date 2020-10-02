#!/usr/bin/env python

'''
Call peaks from bam file of aligned CUTTAG reads
'''

import argparse
import os
import sys
import pandas as pd
from shutil import which
from functools import reduce
import subprocess
import numpy as np
from rgt.helper import get_chrom_sizes_as_genomicregionset
from rgt.THOR.get_extension_size import get_extension_size
from rgt.CoverageSet import CoverageSet
from scipy.stats import binom_test
from rgt.GenomicRegionSet import GenomicRegionSet
from rgt.GenomicRegion import GenomicRegion
from rgt.motifanalysis.Statistics import multiple_test_correction


def parseArgs():
    '''Parse command line arguments.'''
    parser = argparse.ArgumentParser(
            prog='callpeaks.py',
            description='Call peaks on CUTTAG bam files')
    parser.add_argument('-b', '--bam',
                        help='Bam file',
                        required=True,
                        type=str)
    parser.add_argument('-o', '--outfile',
                        help='Output prefix (filename without extension)',
                        required=True,
                        type=str)
    parser.add_argument('-cf', '--controlfile',
                        help='control file for experiment or sample',
                        required=False,
                        type=str)
    parser.add_argument('-cs', '--chromsizes',
                        help='Genome chromsizes file',
                        required=True,
                        type=str)
    parser.add_argument('-minreads', '--minreads',
                        help='Minimum number of reads ',
                        required=False,
                        default=15,
                        type=float)
    parser.add_argument('-minsize', '--minsize',
                        help='Only output peaks greater than\
                              or equal to -min-size',
                        required=False,
                        default=300,
                        type=float)
    parser.add_argument('-pv', '--pvalue',
                        help='Pvalue threshold for binomial peak test\
                             (default 0.05)',
                        required=False,
                        type=float,
                        default=0.05)
    parser.add_argument('-cp', '--correct-pval',
                        help='Correct p-values for multiple testing using\
                             Benjamini/Hochberg ("bh") or Benjamini/Yekutieli ("by") method',
                        required=False,
                        type=str,
                        default="bh")

    return parser.parse_args()


def isTool(name):
    '''return true if 'name' is an executible found in PATH'''
    return which(name) is not None


def filter_bins(c, clookup, min_reads, pvalue_theshold=0.05):
    '''
    Check if bin has more reads than expected by chance with binom_test.
    Ignores bins with fewer than min_reads.
    '''
    p = clookup[c] if c > min_reads else 1
    return True if p < pvalue_theshold else False


def norm(cov):
    '''Set 1/0 and 0/0 case to 0'''
    for c in cov:
        c[np.isnan(c)] = 0
        c[np.isinf(c)] = 0
        c[c > 1] = 0
    return cov


def norm_igg(cov, ctrl):
    '''
    Calculates the fraction of igg per coverage bin
    and scales coverage to the complement.
    coverage = coverage * (1-(iggRPM/covRPM))
    '''
    cov_rpm = np.array(cov.coverage, dtype='object') * (1e6 / float(cov.reads))
    ctrl_rpm = np.array(ctrl.coverage, dtype='object') * (1e6 / float(ctrl.reads))
    frac = np.divide(ctrl_rpm, cov_rpm)
    scale = 1-norm(frac)  # complement igg fraction
    for i in range(len(cov.coverage)):
        cov.coverage[i] = np.rint(cov.coverage[i] * scale[i]).astype(int)


def call_peaks(bam, csizes, pval, min_reads, cfile=None):
    '''
    Call peaks on bam file using pvalue and binomial model.
    Returns GenomeRegionSet with peaks, and CoverageSet with signal.
    '''

    # make chromsizes region set
    rs = get_chrom_sizes_as_genomicregionset(csizes)

    print("calculating extension sizes...")
    # calculate ext size
    ext, _ = get_extension_size(bam, start=0, end=300, stepsize=5)

    print("calculating coverage...")
    # calc coverage
    cov = CoverageSet('coverageset', rs)
    cov.coverage_from_bam(bam_file=bam, extension_size=ext)
    if cfile is not None:
        print(f"Using control file: {cfile}")
        control = CoverageSet('contorl', rs)
        control.coverage_from_bam(bam_file=cfile, extension_size=ext)
        with np.errstate(divide='ignore', invalid='ignore'):
            norm_igg(cov, control)

        # recalc overall coverage
        cov.overall_cov = reduce(lambda x, y: np.concatenate((x, y)),
                                 [cov.coverage[i] for i in range(len(cov.genomicRegions))])

    # total coverage
    s = np.sum(cov.overall_cov)
    # probability of event, a read in a bin, (avg reads/bin )/libsize
    p = np.mean(cov.overall_cov[cov.overall_cov > 0]) / s

    # what is the max coverage
    maxcov = np.max(cov.overall_cov)

    # create dict with probability for each count value
    mc = np.arange(0, maxcov+1, dtype="object")
    d = {count: binom_test((count, s-count), p=p) for count in mc}

    # create GenomicRegionSet to hold peaks
    res = GenomicRegionSet('identified_peaks')

    print("calculating peaks...")
    # iterate through bins in genome, store peaks
    for i, c in enumerate(cov.overall_cov):
        if filter_bins(c, d, min_reads):
            chrom, s, e = cov.index2coordinates(i, rs)
            res.add(GenomicRegion(chrom, s, e+1, data=d[c]))

    # merge ol peaks
    res.merge()

    # merge peaks within ext dist
    rc = res.cluster(ext)

    return rc, cov


def write_bed(res, file, minSize):
    '''Write peak result to bed file'''
    with open(file, 'w') as f:
        for i, c in enumerate(res):
            # write minus log10 p-value to conform with BED format
            p = c.data
            if c.initial > c.final:
                st = c.final
                en = c.initial
                if en-st >= minSize:
                    f.write(f"{c.chrom}\t{st}\t{en}\tPeak{i}\t{p}\t.\n")
            elif c.final-c.initial >= minSize:
                f.write(f"{c.chrom}\t{c.initial}\t{c.final}\tPeak{i}\t{p}\t.\n")


def write_wig(cov, filename):
    """Output coverage in wig format.

    *Keyword arguments:*

    - filename -- filepath
    """
    f = open(filename, 'w')
    i = 0
    for region in cov.genomicRegions:
        print('variableStep chrom=' + str(region.chrom)
              + ' span=' + str(cov.stepsize), file=f)
        c = cov.coverage[i]
        i += 1
        for j in range(len(c)):
            if c[j] != 0:
                print(j * cov.stepsize +
                      ((cov.binsize - cov.stepsize) / 2),
                      c[j], file=f)
    f.close()


def write_bigwig(cov, filename, chrom_file, save_wig=False):
    """Output coverage in bigwig format.

    The path to the chromosome size file <chrom_file> is required.
    This file is tab-separated and assigns
    a chromosome to its size.

    *Keyword arguments:*

    - filename -- filepath
    - chrom_file -- chromosome size file
    """

    tmp_path = filename + '.wig'
    cov.write_wig(tmp_path)

    if not isTool("wigToBigWig"):
        print("wigToBigWig not in PATH")
        sys.exit(1)

    c = ['wigToBigWig', "-clip", tmp_path, chrom_file, filename]
    rc = subprocess.run(c, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if rc.returncode != 0:
        print(f"return: {rc}")

    if not save_wig:
        os.remove(tmp_path)


def main():

    # get command line arguments
    args = parseArgs()

    print(f"Using bam: {args.bam}")
    bf = args.bam

    print(f"Will write peaks: {args.outfile}_peaks.bed")
    of = args.outfile

    cf = args.controlfile
    cs = args.chromsizes
    pvalue = args.pvalue
    minreads = args.minreads
    minsize = args.minsize

    corr = args.correct_pval
    if corr not in ["bh", "by", None]:
        print("Invalid correction method (please pass either 'bh' or 'by'")
        sys.exit(1)

    res, cov = call_peaks(bf, cs, pvalue, minreads, cfile=cf)
    # rpm norm the signal before writing to bigwig
    cov.coverage = np.array(cov.coverage, dtype='object') * (1e6 / float(cov.reads))

    bwfile = of+".bw"
    write_bigwig(cov, bwfile, cs)

    outbed = of+"_peaks.bed"
    write_bed(res, outbed, minsize)

    if corr is not None:
        dat = pd.read_csv(outbed,
                          sep="\t",
                          names=["chr", "start",
                                 "end", "name", "score", "strand"])
        if corr == "bh":
            b, corr = multiple_test_correction(dat.score, method="p")
        elif corr == "by":
            b, corr = b, corr = multiple_test_correction(dat.score, method="n")
        dat["score"] = corr
        dat["score"] = dat["score"].apply(lambda x: -np.log(x))
        dat[dat.score > -np.log10(pvalue)].to_csv(outbed, sep="\t",
                                                  header=False, index=False)


if __name__ == "__main__":
    main()

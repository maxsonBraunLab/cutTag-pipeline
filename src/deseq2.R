library("tidyverse")
library("vsn")
library("pheatmap")
library("ggplot2")
library("RColorBrewer")
library("DESeq2")
library("PoiClaClu")
library("gplots")
library("plyranges")
library("tidyr")

message("setting snakemake parameters")

exp_prefix = snakemake@params[["mark"]]
outdir = snakemake@params[["outdir"]]

input = snakemake@input[["counts"]]
meta = snakemake@input[["meta"]]
genes = snakemake@input[["genes"]]
pcaPlot = snakemake@output[["pcaPlot"]]
pcaPlotVsd = snakemake@output[["pcaPlotVsd"]]
normCount = snakemake@output[["normCounts"]]
lnormCount = snakemake@output[["lnormCounts"]]
sdMeanRld = snakemake@output[["sdMeanRld"]]
sdMeanVsd = snakemake@output[["sdMeanVsd"]]
sampleDistPlotVsd = snakemake@output[["sampleDistVsd"]]
sampleDistPlotRld = snakemake@output[["sampleDistRld"]]
rds = snakemake@output[["rds"]]


# read in genes file
genestab = read_tsv(genes, col_names=c("seqnames", "start", "end", "name", "score", "strand")) %>% GRanges()

# counts join
peaks = read_tsv(input) %>% select(chrom, start, end, peak) %>% rename(peak='name')

# read in counts table
counts = read.delim(input, header=T, stringsAsFactors = F)
rownames(counts) = counts$peak
counts = counts[,7:ncol(counts)]

# read in metadata
meta <- read.csv(meta, header = T, stringsAsFactors = F)

# make the meta reflect the available sample names
meta = meta[meta$sample %in% colnames(counts),]

message('calculating deseq')
# make deseqdataset
ddsMat <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~condition)
dds <- DESeq(ddsMat)

# write normalized counts tables
normCounts <- counts(dds, normalized=TRUE)
normCounts <- 0.00001+(normCounts) # ensure nonzero 
lnormCounts <- log2(normCounts)
lnormCounts <- as.data.frame(lnormCounts)
lnormCounts$ID <- counts$peak
write.csv(normCounts, normCount)
write.csv(lnormCounts, lnormCount) 

# calculate rlog transform
rld <- rlog(dds, blind=FALSE)
vsd <- varianceStabilizingTransformation(dds, blind = TRUE,
  fitType = "parametric")
ntd <- normTransform(dds)

message('plotting count transforms')
# plot the transform for the first two samples
png(sdMeanRld)
rldCounts <- meanSdPlot(assay(rld))
rldCounts$gg + labs(title="meanSdPlot - rlog transformed")
dev.off()

png(sdMeanVsd)
vsdCounts <- meanSdPlot(assay(vsd))
vsdCounts$gg + labs(title="meanSdPlot - vsd transformed")
dev.off()

message('plotting poisson distance sample cross correlation')
# plot sample distances 
png(sampleDistPlotVsd)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(meta$sample, meta$condition, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,
         main=paste(exp_prefix, "vsd sample distance matrix"))
dev.off()

png(sampleDistPlotRld)
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(meta$sample, meta$condition, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,
         main=paste(exp_prefix, "rld sample distance matrix"))
dev.off()

# plot PCA rld
png(pcaPlot)
plotPCA(rld, intgroup = c("condition")) + labs(title=paste0(exp_prefix,"-rld"))
dev.off()

# plot PCA vsd
png(pcaPlotVsd)
plotPCA(rld, intgroup = c("condition")) + labs(title=paste0(exp_prefix,"-vsd"))
dev.off()

# make contrasts
c = combn(unique(meta$condition), 2)
lsc=split(c, rep(1:ncol(c), each = nrow(c)))

message('writing contrast results')
# write differential results for each contrast
for (k in 1:length(lsc)) {
    cl=lsc[[k]]
    rname=paste0(cl[1],"vs",cl[2])
    res=results(dds, 
            independentFiltering = TRUE,
            cooksCutoff = FALSE,
            contrast = c("condition", cl[1], cl[2]),
            alpha=0.05)

    # estimate shrinkage 
    resLFC <- lfcShrink(dds, coef=2, type="apeglm")

    #plotMA
    maplot=paste0(outdir,"/",exp_prefix,"/",exp_prefix,"-",rname,"plotMA.png")
    png(maplot)
    par(mfrow=c(1,2), mar=c(4,4,2,1))
    xlim <- c(1,1e5); ylim <- c(-2,2)
    plotMA(resLFC, xlim=xlim, ylim=ylim, alpha=0.05, main=paste(rname, "apeglm"))
    plotMA(res, xlim=xlim, ylim=ylim, alpha=0.05, main=paste(rname, "normal"))
    dev.off()

    tableName=paste0(outdir,"/",exp_prefix,"/",cl[1],'-',cl[2],'-',exp_prefix,'-','peaks-annot.tsv')
   
    # annottate peaks
    diffexp=res %>% 
        as_tibble(rownames='name') %>%
        left_join(peaks) %>%
        GRanges() %>%
        join_nearest(genestab, suffix=c(".peak", ".gene")) %>%
        as_tibble() 
    write_tsv(diffexp, path=tableName)

    # differential peak files lfc > 0 | lfc > 0  padj < 0.05 | padj < 0.01
    deup05=paste0(outdir,"/",exp_prefix,"/",cl[1],'-',cl[2],'-',exp_prefix,'-','differential-up-05.bed')
    dedown05=paste0(outdir,"/",exp_prefix,"/",cl[1],'-',cl[2],'-',exp_prefix,'-','differential-down-05.bed')
    deup01=paste0(outdir,"/",exp_prefix,"/",cl[1],'-',cl[2],'-',exp_prefix,'-','differential-up-01.bed')
    dedown01=paste0(outdir,"/",exp_prefix,"/",cl[1],'-',cl[2],'-',exp_prefix,'-','differential-down-01.bed')
    # sig up peaks
    diffexp %>%
        filter(log2FoldChange > 0, padj < 0.05) %>%
        select(seqnames, start, end, name.peak,score,strand) %>%
        write_tsv(deup05, col_names=FALSE)

    # sig down peaks
    diffexp %>%
        filter(log2FoldChange < 0, padj < 0.05) %>%
        select(seqnames, start, end, name.peak,score,strand) %>%
        write_tsv(dedown05, col_names=FALSE)

    # sig up peaks
    diffexp %>%
        filter(log2FoldChange > 0, padj < 0.01) %>%
        select(seqnames, start, end, name.peak,score,strand) %>%
        write_tsv(deup01, col_names=FALSE)

    # sig down peaks
    diffexp %>%
        filter(log2FoldChange < 0, padj < 0.01) %>%
        select(seqnames, start, end, name.peak,score,strand) %>%
        write_tsv(dedown01, col_names=FALSE)

    summaryName=paste0(outdir,"/",exp_prefix,"/",cl[1],"-",cl[2],"-",exp_prefix,"-diffexp-summary.txt")
    diffexp %>%
        summarise(Condition=paste0(exp_prefix,"-",rname),
                  DP05=nrow(filter(.,padj<0.05)),
                  DP01=nrow(filter(.,padj<0.01))) %>%
        write_tsv(summaryName)

}

# save RDS of deseq object
saveRDS(dds, rds)

library(pheatmap)
library(dply)
library(tidyr)
library(dendextend)
library(randomcoloR)

# the number of clusters for the heatmap
numclusters=4

# heatmap title
heattitle="H3K27Ac-C-vs-D"

heatmap_filename=paste0(heattitle, "-heatmap.svg")

# peaks annot
table <- "/home/groups/MaxsonLab/vancampe/cutTag-pipeline-molm6-2020-10-08/data/deseq2/27Ac/C-D-27Ac-peaks-annot.tsv"
tab <- read_delim(table, delim = "\t")
sig <- tab %>% filter(padj < 0.05)

annot_filename=gsub(".tsv", "-05-clust.tsv", basename(table))

direction <- sig %>% 
  select(name.peak,log2FoldChange) %>% 
  mutate(direction=ifelse(log2FoldChange>0,"up","down")) %>%
  data.frame() 

rownames(direction) = direction[,1]
direction[,1]=NULL

# read in counts
counts <- "/home/groups/MaxsonLab/vancampe/cutTag-pipeline-molm6-2020-10-08/data/deseq2/27Ac/27Ac-lognormcounts.csv"
counts <- read_csv(counts) %>% rename(X1 = "peak")

# read in metadata file with condition information
meta <- "/home/groups/MaxsonLab/vancampe/cutTag-pipeline-molm6-2020-10-08/src/deseq2_metadata.csv"
metatab <- read_csv(meta)

# subset counts by sig peaks
sigcounts <- counts %>% filter(peak %in% sig$name.peak)

# organize samples by condition
meta_conditions <- metatab %>% group_by(condition) %>%
  summarise(samples=paste(sample, collapse = ","))

# rowMean counts for each sample condition
for (c in meta_conditions$condition) {
  samples = filter(meta_conditions, condition==c)$samples
  s=unlist(strsplit(samples, ","))
  sigcounts <- sigcounts %>% mutate(!!c := rowMeans(sigcounts[s])) 
}

# retain mean counts per condition 
sigcounts <- sigcounts[append("peak", meta_conditions$condition)]

# scale the mean counts
scaled <- t(apply(as.matrix(sigcounts[,-1]), 1, scale))
colnames(scaled) = colnames(sigcounts[-1])

# re-assign rownames
rownames(scaled) = sigcounts$peak

# generate clusters
clust=hclust(dist(scaled), method="ward.D2")
clustdf=data.frame(cluster=cutree(as.dendrogram(clust), k = numclusters))
annots=merge.data.frame(clustdf, direction, by=0, all = T) %>% select(-log2FoldChange)
rownames(annots)=annots[,1]
annots[,1]=NULL

# get a numebr of unique colors
cols=RColorBrewer::brewer.pal(numclusters, name="Set1")

# name the colors by cluster
names(cols) = unique(clustdf$cluster)

# set cluster colors
colors = list(
  cluster=cols
)

# create heatmap
heat <- pheatmap(scaled,
         main=heattitle,
         show_rownames = F,
         show_colnames = T,
         cluster_rows = T,
         cluster_method = "ward.D2",
         annotation_row = annots,
         annotation_colors = colors)

save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
  svg(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png(heat, heatmap_filename)

# annotate cluster and direction
annots[["name.peak"]] = rownames(annots)
sig %>% left_join(annots) %>% write_delim(annot_filename)



#!/usr/bin/env Rscript

# This script generates plots from Diffbind data in an DBA object from an RDS file

# load libraries
library(DiffBind)
# library(rtracklayer)
library(GenomicRanges)
# library(ggplot2)
library(dplyr)
library(stringr)
library(optparse)

# function to parse command line arguments
get_args <- function(){
	option_list <- list(
		make_option(c("-i", "--dba_rds_file"), action = "store", default = NA, type = "character", help = "Character string specifying path to an RDS file containing a DBA object outputted by Diffbind."),
		make_option(c("-p", "--pval_cutoff"), action = "store", default = 0.05, type = "numeric", help = "Number specifying adjusted p-value threshold to be used by Diffbind [default %default]."),
		make_option(c("-o", "--outdir"), action = "store", default = NA, type = "character", help = "Character string specifying path to folder in which to store Diffbind output.")
	)
	opt <- parse_args(OptionParser(option_list = option_list))
	return(opt)
}


# main driver function
diffbind_main <- function(){
	#-------- general setup --------
	# get command line arguments
	message(date())
	message("Getting command line arguments...")

	args <- get_args()
	dba_rds_file <- args$dba_rds_file
	pval_cutoff <- args$pval_cutoff
	diffbind_outdir <- args$outdir
	
	# if output directory doesn't exist, create it
	if (!dir.exists(diffbind_outdir)) {
		message("Creating output directory: ", diffbind_outdir)
		dir.create(diffbind_outdir, recursive = TRUE, mode = "774")	
	}
	
	message(date())
	message("Reading in DBA object from RDS file...")
	dba_obj <- readRDS(dba_rds_file)

	# get mark from input file name
	mark <- stringr::str_split(
		string = basename(dba_rds_file),
		pattern = "_",
		simplify = TRUE)[,1]


	#------- generate plots -------
	message(date())
	message("Generating plots for all samples in mark/protein: ", mark)

	# sample correlation plot w/all mark samples
	corr_plot_outfile <- file.path(diffbind_outdir, paste0(mark, "_sample_correlation.pdf"))
	pdf(
		file = corr_plot_outfile,
		width = 10,
		height = 12
	)
	dba.plotHeatmap(dba_obj, correlations = TRUE)
	dev.off()

	#  PCA plot w/all mark samples, color by DBA_CONDITION
	pca_outfile <- file.path(diffbind_outdir, paste0(mark, "_pca.pdf"))
	pdf(
		file = pca_outfile,
		width = 10,
		height = 10
	)
	dba.plotPCA(
		dba_obj, 
		attributes = DBA_CONDITION, 
		label = DBA_ID
	)
	dev.off()

	# heatmap w/all mark samples
	heatmap_outfile <- file.path(diffbind_outdir, paste0(mark, "_heatmap.pdf"))
	pdf(file = heatmap_outfile,
		width = 10,
		height = 32
	)
	dba.plotHeatmap(dba_obj, correlations = FALSE)
	dev.off()


	# loop through all contrasts and export plots
	dba_meta <- dba.show(dba_obj, bContrasts = TRUE)

	for (i in 1:nrow(dba_meta)) {
		de_counts <- dba_meta[i, "DB.DESeq2"]
		group1 <- dba_meta[i, "Group"]
		group2 <- dba_meta[i, "Group2"]
		contrast <- paste0(group1, "-vs-", group2)

		# create output directory for contrast
		diffbind_contrast_outdir <- file.path(diffbind_outdir, contrast)

		if (!dir.exists(diffbind_contrast_outdir)) {
			message(date())
			message("Creating output directory: ", diffbind_contrast_outdir)
			dir.create(diffbind_contrast_outdir, recursive = TRUE, mode = "774")	
		}

		#------- generate plots -------
		message(date())
		message("Generating plots for contrast: ", contrast)
		
		# MA plot per contrast
		ma_plot_outfile <- file.path(diffbind_contrast_outdir, paste0(mark, "_", contrast, "_MAplot.pdf"))
		pdf(
			file = ma_plot_outfile,
			width = 10,
			height = 10
		)
		# use FDR for significance cutoff (bUsePval = FALSE)
		dba.plotMA(
			dba_obj, 
			contrast = i, 
			th = pval_cutoff, 
			bUsePval = FALSE, 
			factor = paste0(mark, " ")
		)
		dev.off()

		# volcano plot per contrast
		volcano_plot_outfile <- file.path(diffbind_contrast_outdir, paste0(mark, "_", contrast, "_volcano.pdf"))
		pdf(
			file = volcano_plot_outfile,
			width = 10,
			height = 10
		)
		dba.plotVolcano(
			dba_obj, 
			contrast = i, 
			th = pval_cutoff, 
			bUsePval = FALSE, 
			factor = paste0(mark, " ")
		)
		dev.off()

	}
}


# run main function
diffbind_main()

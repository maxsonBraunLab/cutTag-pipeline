#!/usr/bin/env Rscript

# This script runs Diffbind for a given input mark/protein of interest to normalize counts in peaks and perform differential analysis
# Adapted from diffbind.R script in MaxsonBraunLab dba_obj_seq pipeline

# load libraries
library(DiffBind)
# library(rtracklayer)
library(GenomicRanges)
library(plyranges)
# library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(optparse)

# function to parse command line arguments
get_args <- function(){
	option_list <- list(
		make_option(c("-m", "--metadata_file"), action = "store", default = NA, type = "character", help = "Character string specifying path to a comma-delimited file (.csv) containing metadata required for Diffbind."),
		make_option(c("-c", "--consensus_peak_file"), action = "store", default = NA, type = "character", help = "Character string specifying path to a tab-delimited BED file containing consensus peaks to be used by Diffbind."),
		make_option(c("-g", "--refseq_file"), action = "store", default = NA, type = "character", help = "Character string specifying path to a BED file containing RefSeq gene coordinates. Used for annotating differential peaks."),
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
	metadata_df <- read.csv(file = args$metadata_file)
	consensus_peak_file <- args$consensus_peak_file
	refseq_file <- args$refseq_file
	pval_cutoff <- args$pval_cutoff
	diffbind_outdir <- args$outdir
	
	# if output directory doesn't exist, create it
	if (!dir.exists(diffbind_outdir)) {
		message("Creating output directory: ", diffbind_outdir)
		dir.create(diffbind_outdir, recursive = TRUE, mode = "774")	
	}

	# read in gene coordinates from refseq file
	refseq_granges <- 
		readr::read_tsv(
			file = refseq_file,
			col_names = c("seqnames", "start", "end", "name", "score", "strand")) |>
		GRanges()
	

	#-------- diffbind setup --------
	message(date())
	message("Subsetting metadata for mark/protein of interest...")

	# get mark from input consensus peak file name
	mark <- stringr::str_split(
		string = basename(consensus_peak_file),
		pattern = "_",
		simplify = TRUE)[,1]

	metadata_df_subset <- metadata_df |>
		dplyr::filter(Factor == mark)

	message(date())
	message("Creating GRanges object from consensus peaks...")
	consensus_peaks_granges <- 
		read.table(file = consensus_peak_file, col.names=c("seqnames", "start", "end")) |>
		dplyr::mutate(name = paste0("peak", dplyr::row_number())) |>
		GRanges()


	#-------- run diffbind --------
	message(date())
	message("Running Diffbind...")

	# create DBA object from metadata
	message(date())
	message("Creating DBA object and counting reads in regions...")
	dba_obj <- dba(sampleSheet = metadata_df_subset)
	# count reads in regions
	dba_obj <- dba.count(dba_obj, peaks = consensus_peaks_granges)

	# normalize by full library size
	message(date())
	message("Normalizing...")
	dba_obj <- dba.normalize(
		dba_obj, 
		normalize = DBA_NORM_LIB, 
		library = DBA_LIBSIZE_FULL
	)
	# set up contrasts based on conditions for differential analysis
	message(date())
	message("Setting up contrasts...")
	dba_obj <- dba.contrast(
		dba_obj,
		minMembers = 2,
		categories = DBA_CONDITION
	)
	# run differential analysis
	message(date())
	message("Running differential analysis...")
	dba_obj <- dba.analyze(dba_obj, bParallel = TRUE)


	#-------- export diffbind data --------
	# to do: 
	# - turn this into separate function(export_diffbind_data()), passing as arguments: dba_obj, pval_cutoff, diffbind_outdir, mark
	# - add nearest gene to each consensus region before outputting final files

	# message(date())
	# message("Exporting Diffbind data...")
	# export_diffbind_data(
	# 	dba_object = dba_obj,
	# 	pvalue_cutoff = pval_cutoff,
	# 	mark_of_interest = mark,
	# 	output_dir = diffbind_outdir
	# )

	message(date())
	message("Exporting Diffbind DBA object to RDS file...")
	rds_outfile <- file.path(diffbind_outdir, paste0(mark, "_DBAobj.rds"))
	saveRDS(object = dba_obj, file = rds_outfile)
	message("DBA object saved to: ", rds_outfile)

	message(date())
	message("Exporting normalized counts table...")
	normcounts_df <- dba.peakset(dba_obj, bRetrieve = TRUE, DataType = DBA_DATA_FRAME)

	normcounts_outfile <- file.path(diffbind_outdir, paste0(mark, "_normcounts.csv"))
	write.csv(x = normcounts_df, file = normcounts_outfile, row.names = FALSE)
	message("Normalized counts saved to: ", normcounts_outfile)

	# #------- generate plots -------
	# message(date())
	# message("Generating plots for all samples in mark/protein: ", mark)

	# # sample correlation plot w/all mark samples
	# corr_plot_outfile <- file.path(diffbind_outdir, paste0(mark, "_sample_correlation.pdf"))
	# pdf(
	# 	file = corr_plot_outfile,
	# 	width = 12,
	# 	height = 12
	# )
	# dba.plotHeatmap(dba_obj, correlations = TRUE)
	# dev.off()

	# #  PCA plot w/all mark samples, color by DBA_CONDITION
	# pca_outfile <- file.path(diffbind_outdir, paste0(mark, "_pca.pdf"))
	# pdf(
	# 	file = pca_outfile,
	# 	width = 12,
	# 	height = 12
	# )
	# dba.plotPCA(
	# 	dba_obj, 
	# 	attributes = DBA_CONDITION, 
	# 	label = DBA_ID
	# )
	# dev.off()

	# # heatmap w/all mark samples
	# heatmap_outfile <- file.path(diffbind_outdir, paste0(mark, "_heatmap.pdf"))
	# pdf(file = heatmap_outfile,
	# 	width = 12,
	# 	height = 34
	# )
	# dba.plotHeatmap(dba_obj, correlations = FALSE)
	# dev.off()


	# loop through all contrasts and export differential analysis results and plots
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

		# #------- generate plots -------
		# message(date())
		# message("Generating plots for contrast: ", contrast)
		
		# # MA plot per contrast
		# ma_plot_outfile <- file.path(diffbind_contrast_outdir, paste0(mark, "_", contrast, "_MAplot.pdf"))
		# pdf(
		# 	file = ma_plot_outfile,
		# 	width = 12,
		# 	height = 12
		# )
		# # use FDR for significance cutoff (bUsePval = FALSE)
		# dba.plotMA(
		# 	dba_obj, 
		# 	contrast = i, 
		# 	th = pval_cutoff, 
		# 	bUsePval = FALSE, 
		# 	factor = paste0(mark, " ")
		# )
		# dev.off()

		# # volcano plot per contrast
		# volcano_plot_outfile <- file.path(diffbind_contrast_outdir, paste0(mark, "_", contrast, "_volcano.pdf"))
		# pdf(
		# 	file = volcano_plot_outfile,
		# 	width = 12,
		# 	height = 12
		# )
		# dba.plotVolcano(
		# 	dba_obj, 
		# 	contrast = i, 
		# 	th = pval_cutoff, 
		# 	bUsePval = FALSE, 
		# 	factor = paste0(mark, " ")
		# )
		# dev.off()

		#------- generate output tables -------
		message(date())
		message("Exporting differential analysis results for contrast: ", contrast)

		# if de_counts is null, need to put in separate if() statement to check for null
		# combining the first if/else-if statements will fail if de_counts is null
		if (is.null(de_counts) | length(de_counts) == 0) {
			message("WARNING: ", contrast, " did not have any significant differential peaks.")
			next
		} else if (is.na(de_counts) | de_counts == 0) {
			message("WARNING: ", contrast, " did not have any significant differential peaks.")
			next
		} else {
			# use FDR for significance cutoff (bUsePval = FALSE)
			dba_intervals <- dba.report(dba_obj, contrast = i, th = pval_cutoff, bUsePval = FALSE) |>
				as.data.frame() |>
				GRanges() |>
				join_nearest(y = refseq_granges, suffix = c(".peak", ".gene"), distance = TRUE) |>
        		as.data.frame()

			# rename Conc_ prefix in columns to "log2MeanReads"
			colnames(dba_intervals) <- stringr::str_replace(
				string = colnames(dba_intervals), 
				pattern = "Conc", 
				replacement = "log2MeanReads")

			upregulated_intervals <- dba_intervals |> 
				dplyr::filter(Fold > 0)
			downregulated_intervals <- dba_intervals |> 
				dplyr::filter(Fold < 0)

			# define output files
			all_sig <- file.path(diffbind_contrast_outdir, paste0(contrast, "-all_sig.txt"))
			up_output <- file.path(diffbind_contrast_outdir, paste0(contrast, "-up-", pval_cutoff, ".txt"))
			dn_output <- file.path(diffbind_contrast_outdir, paste0(contrast, "-dn-", pval_cutoff, ".txt"))
			up_bed <- file.path(diffbind_contrast_outdir, paste0(contrast, "-up-", pval_cutoff, ".bed"))
			dn_bed <- file.path(diffbind_contrast_outdir, paste0(contrast, "-dn-", pval_cutoff, ".bed"))

			write.table(
				x = dba_intervals, 
				file = all_sig, 
				sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
			
			write.table(
				x = upregulated_intervals, 
				file = up_output, 
				sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
			
			write.table(
				x = downregulated_intervals, 
				file = dn_output, 
				sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
			
			upregulated_intervals |> 
				dplyr::select(seqnames, start, end) |> 
				write.table(file = up_bed, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
			
			downregulated_intervals |> 
				dplyr::select(seqnames, start, end) |> 
				write.table(file = dn_bed, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
		}
	}
}


# run main function
diffbind_main()

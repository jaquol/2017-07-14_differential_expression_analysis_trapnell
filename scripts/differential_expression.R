#==================================================================================================
# Created on: 2017-07-21
# Usage: ./differential_expression.R
# Author: Javier Quilez (GitHub: jaquol)
# Goal: Performs differential expression analysis between 2 conditions using sleuth
#==================================================================================================



#==================================================================================================
# Configuration variables and paths
#==================================================================================================

# installation --only needs to be done once
#source("http://bioconductor.org/biocLite.R")
#biocLite("rhdf5")
#install.packages("devtools")
#devtools::install_github("pachterlab/sleuth")
#install.packages("RSQLite")
#biocLite("biomaRt")
#install.packages("DESeq2")
#install.packages("pheatmap")
#install.packages("cowplot")

# load R libraries
library("sleuth")
library("RSQLite")
library("dplyr")
library("tidyr")
library('cowplot')

# variables
data_type <- "rnaseq"
project <- "misc"
analysis <- "2017-07-14_differential_expression_analysis_trapnell"
species = "homo_sapiens"
assembly_version = "hg38_mmtv"
sequencing_type <- 'paired_end'

# paths
HOME <- "/users/GR/mb/jquilez"
SAMPLES <- paste0(HOME, "/data/", data_type, "/samples")
ANALYSIS <- paste0(HOME, "/projects/", project, "/analysis/", analysis)
ifile_metadata <- paste0(ANALYSIS, "/tables/sample_id_to_srr.txt")



#==================================================================================================
# Import metadata
#==================================================================================================

# import metadata
metadata <- read.delim(ifile_metadata, header = FALSE)
metadata <- metadata[, c("V1", "V3")]
names(metadata) <- c('sample', 'condition')
samples <- as.vector(metadata$sample)



#==================================================================================================
# Import GENCODE V24 annotation
#==================================================================================================

# transcript quantifications are generated with kallisto using GENCODE V24 annotation
# gene and transcript ids can be obtained from the kallisto's quantification file for any of the samples

# read in file
itsv <- paste0(SAMPLES, '/', samples[1], '/quantifications/kallisto/', assembly_version, '/', sequencing_type, '/abundance.tsv')
target_mapping <- read.delim(itsv)

# extract information from target_id
target_mapping$target_id_raw <- target_mapping$target_id
target_id_names <- c('transcript_id', 'gene_id', 'ott1', 'ott2', 'transcript_name', 'gene_name', 'transcript_length', 'gene_type', 'empty')
target_mapping <- target_mapping %>% separate(target_id_raw, target_id_names, "\\|")

# select specific columns
target_mapping <- target_mapping[, c("target_id", "gene_id", 'gene_name')]



#==================================================================================================
# Principal component analysis (PCA) using the transcript-level counts
#==================================================================================================

# transcript-level counts estimated with kallisto
# dPCA performed with sleuth

diff_exp_sleuth <- function(condition1, condition2, level) {

	# re-format metadata to meet sleuth's required format
	# later on, the 'sleuth_prep' function of sleuth will assign 0=untreated and 1=treated condition
	# importantly, by default such function will set '0' to first alphabetical match in the values
	# found in the 'condition' field so here we force that 0=condition1
	tab_metadata <- metadata
	rownames(tab_metadata) <- tab_metadata$sample

	conditions <- c()
	for (i in 1:nrow(tab_metadata)) {
		tag <- tab_metadata[i, 'condition']
		if (tag == condition1) {
			conditions <- c(conditions, 0)
		}
		else if (tag == condition2) {
			conditions <- c(conditions, 1)
		}
		else {
			conditions <- c(conditions, -1)
		}
	}
	tab_metadata$condition_name <- tab_metadata$condition
 	tab_metadata$condition <- conditions

	# subset samples to include only those from the 2 conditions compared
	cond1 <- tab_metadata$condition == 0
  	cond2 <- tab_metadata$condition == 1
 	tab_metadata <- tab_metadata[cond1 | cond2, ]

 	# get paths to data
 	samples <- tab_metadata$sample
	paths <- c()
	for (s in samples) {
		p <- paste0(SAMPLES, "/", s, "/quantifications/kallisto/", assembly_version, '/', sequencing_type)
		paths <- c(paths, p)
	}
	tab_metadata$path <- paths

	# (1) load the kallisto processed data and make a regression model using 'condition' as the dependent variable
	# also, `target_mapping` adds annotation information (e.g. gene_id)
	if (level == 'gene') {
	  	so <- sleuth_prep(tab_metadata, ~ condition, target_mapping = target_mapping, aggregation_column = 'gene_id')	
	} else if (level == 'transcript') {
	  	so <- sleuth_prep(tab_metadata, ~ condition, target_mapping = target_mapping)
	}
  
	# (2) estimate parameters for the sleuth response error measurement (full) model as responding to the 'condition' factor
  	so <- sleuth_fit(so)
  
	# (3) Create another model where the gene expression is not dependent on any factor.
  	so <- sleuth_fit(so, ~1, 'reduced')
  
	# (4.1) Run a likelihood ratio test (LRT) between the two models to see what transcripts appear 
	# to really be affected by the time factor value
	so <- sleuth_lrt(so, 'reduced', 'full')

	# (4.2) Run the Wald test (WT), a statistical tests which:
	# - is somewhat related to the LRT and is also used to test for differential expression
	# - LRT is considerd a better test than the WT but
	# - WT is used becase it generates the beta statistic, which approximates to the fold change in expression between
	# the 2 condition tested, which is typically reported in differential expression analysis
	so <- sleuth_wt(so, paste0('condition'))

	# export normalised transcript abundance values (no needed for level = gene)
	if (level == 'transcript') {
		condition1_name <- tolower(condition1)
		condition2_name <- tolower(condition2)
		otab = paste0(ANALYSIS, "/tables/normalized_abundance_", level, "_level_sleuth_", condition1_name, "_", condition2_name, ".tsv")
 		write.table(kallisto_table(so), otab, sep = "\t", quote = FALSE, row.names = FALSE)
 	}

	# add beta (b), beta's standard error (se_b) and the mean expression in the samples (mean_obs)
 	res_lrt <- sleuth_results(so, 'reduced:full', test_type = 'lrt')
 	res_wt <- sleuth_results(so, 'condition')
 	res <- merge(res_lrt, res_wt[, c('target_id', 'b', 'se_b', 'mean_obs')], on = 'target_id', sort = FALSE)
	
	# export
	condition1_name <- tolower(condition1)
	condition2_name <- tolower(condition2)
	otab = paste0(ANALYSIS, "/tables/differential_expression_analysis_", level, "_level_sleuth_", condition1_name, "_", condition2_name, ".tsv")
 	write.table(res, otab, sep = "\t", quote = FALSE, row.names = FALSE)

 	# perform principal component analysis (PCA)
	png(paste0(ANALYSIS, "/figures/pca_tpm_", level, ".png"), res = 150, w = 1500, h = 500)
	par(mfrow = c(1, 2))
	p1 <- plot_pca(so, pc_x = 1L, pc_y = 2L, color_by = 'condition_name')
	p2 <- plot_pc_variance(so, PC_relative = TRUE)
	plot_grid(p1, p2, align = "h")
	dev.off()

 	return(so)

}

sleuth_object_transcript <- diff_exp_sleuth("scramble", "HOXA1KD", 'transcript')
sleuth_object = diff_exp_sleuth("scramble", "HOXA1KD", 'gene')




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
#library("RSQLite")
#library("dplyr")
#library("tidyr")
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

pca_sleuth <- function() {

	# re-format metadata to meet sleuth's required format
	# later on, the 'sleuth_prep' function of sleuth will assign 0=untreated and 1=treated condition
	# importantly, by default such function will set '0' to first alphabetical match in the values
	# found in the 'condition' field so here we force that 0=condition1
	tab_metadata <- metadata
	rownames(tab_metadata) <- tab_metadata$sample

 	# get paths to data
	paths <- c()
	for (s in samples) {
		p <- paste0(SAMPLES, "/", s, "/quantifications/kallisto/", assembly_version, '/', sequencing_type)
		paths <- c(paths, p)
	}
	tab_metadata$path <- paths

	# (1) load the kallisto processed data and make a regression model using 'condition' as the dependent variable
  	so <- sleuth_prep(tab_metadata, ~ condition, target_mapping = target_mapping)

	png(paste0(ANALYSIS, "/figures/pca_tpm.png"), res = 150, w = 1500, h = 500)
	par(mfrow = c(1, 2))
	p1 <- plot_pca(so, pc_x = 1L, pc_y = 2L, units = 'tpm', color_by = 'condition')
	p2 <- plot_pc_variance(so, units = 'tpm', PC_relative = TRUE)
	plot_grid(p1, p2, align = "h")
	dev.off()

}

pca_sleuth()



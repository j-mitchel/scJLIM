
#' Example simulated dataset containing a true colocalization in one of two cell types
#'
#' @format A list of various data objects:
#' \describe{
#'   \item{1}{scRNA seq normalized and log transformed counts (genes x cells)}
#'   \item{2}{scRNA seq metadata with donors, library size, and/or other covariates (cells x covariates)}
#'   \item{3}{scRNA seq expression PCs (cells x PCs)}
#'   \item{4}{scRNA donor genotypes, maf > 0.1, but can be different (SNPs x donors)}
#'   \item{5}{GWAS summary statistics (SNPs as rows) (CHR, BP, P as columns)}
#'   \item{6}{reference LD matrix}
#'   \item{7}{the preprogrammed cell type labels (unknown in real data)}
#'   \item{8}{the preprogrammed causal GWAS and eQTL SNP (unknown in real data)}
#'   \item{9}{the preprogrammed cell type with the eQTL (unknown in real data)}
#'   \item{10}{the preprogrammed eGene (unknown in real data)}
#' }
#' @source Generated using 'data-raw/prep_sim.R'
"sim_dat_all"

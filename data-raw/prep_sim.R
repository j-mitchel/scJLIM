
## saving one example simulation for use with the tutorial

library(VariantAnnotation)
library(scater)
library(cowplot)
library(umap)

# # a positive simulation with two cell types of equal sizes
# sims_type <- 1
# sims_subtype <- 1
# 
# # which sim to pick 
# # sim_ndx <- 3
# # gw_ndx <- 3
# sim_ndx <- 2
# gw_ndx <- 1
# 
# base_dir <- '/home/jmitchel/data/GWAS_sim_v3'
# 
# loci_coords <- read.csv(paste0(base_dir,'/loci_coords.csv'),row.names = 1,skip = 1, header = FALSE)
# 
# sc_sims_all <- readRDS(file=paste0(base_dir,'/splat_sims/sc_sims_all6.rds'))
# # sc_sims_all <- readRDS(file=paste0(base_dir,'/splat_sims/sc_sims_all7.rds'))
# 
# x <- sc_sims_all[[sim_ndx]]
# sims_iter <- x[[sims_type]]
# sim_select <- sims_iter[[gw_ndx]][[sims_subtype]]
# 
# # load the full eqtl genotypes matrix
# vcf_file <- paste0(base_dir,'/hapgen_gwas/sim_',sim_ndx,'/eqtl_full.vcf')
# donor_geno_eqtl <- readVcf(vcf_file, "hg19")
# donor_geno_eqtl <- as.data.frame(t(as(genotypeToSnpMatrix(donor_geno_eqtl)$genotype, "numeric"))) 
# snps_keep <- which(!is.na(donor_geno_eqtl[,1]))
# donor_geno_eqtl <- donor_geno_eqtl[snps_keep,]
# 
# # load the GWAS sumstats
# main_tr_lst <- readRDS(file=paste0(base_dir,'/hapgen_gwas/sim_',sim_ndx,'/gw_sumstats.rds'))
# 
# # load the reference panel for LD
# ref_matrix <- readRDS(file=paste0(base_dir,'/hapgen_gwas/sim_',sim_ndx,'/ref_panel_jlim.rds'))
# 
# # get the causal snp
# causal_snp <- read.table(file=paste0('/home/jmitchel/data/GWAS_sim_v3/hapgen_gwas/sim_',sim_ndx,'/causal_rh_snps.txt'),header=FALSE)
# causal_snp <- strsplit(as.character(causal_snp),split=',')[[1]][[1]]
# 
# chr <- loci_coords[sim_ndx,1]
# 
# main_tr <- main_tr_lst[[gw_ndx]]
# 
# rownames(main_tr) <- main_tr$all_snp_rsid
# colnames(main_tr) <- c('SNP','pval','position','beta','varbeta','snp_ndx')
# 
# eQTL_ndx <- which(sim_select@rowRanges@elementMetadata$eSNP.ID==causal_snp)
# eQTL_ndx <- eQTL_ndx[1]
# eQTL_group <- sim_select@rowRanges@elementMetadata$eQTL.group[eQTL_ndx]
# 
# eQTL_gene <- sim_select@rowRanges@partitioning@NAMES[eQTL_ndx]
# eQTL_gene_split <- strsplit(eQTL_gene,split='_')[[1]]
# eQTL_gene <- paste0(eQTL_gene_split[1],'-',eQTL_gene_split[2])
# 
# # make cell names unique
# new_cellnames <- sapply(1:nrow(sim_select@colData),function(x){
#   paste0('cell',x)
# })
# rownames(sim_select@colData) <- new_cellnames
# sim_select@colData$Cell <- new_cellnames
# colnames(counts(sim_select)) <- new_cellnames
# 
# # normalize expression counts
# sim_select <- logNormCounts(sim_select)
# 
# # extract the normalized counts
# tmp_norm <- sim_select@assays@data@listData[["logcounts"]]
# 
# # renaming genes
# rownames(tmp_norm) <- sapply(rownames(tmp_norm),function(x){
#   mysplit=strsplit(x,split='_')[[1]]
#   return(paste0(mysplit[1],'-',mysplit[2]))
# })
# rownames(sim_select) <- rownames(tmp_norm)
# 
# ### adding some more variation between cells so the PCs don't capture the eQTL effect
# g_alter <- sample(rownames(tmp_norm),50)
# cell_var_vals <- .2*rnorm(ncol(tmp_norm))
# for (g in g_alter) {
#     tmp_norm[g,] <- tmp_norm[g,] + cell_var_vals
# }
# sim_select@assays@data@listData[["logcounts"]] <- tmp_norm
# 
# # get principal components of expression data
# sim_select <- runPCA(sim_select, ncomponents = 5, scale=TRUE)
# pcs <- reducedDims(sim_select)@listData[["PCA"]]
# # ures <- umap(pcs)
# 
# # extract cell to donor information
# covar_df <- sim_select@colData[,c('Sample'),drop=FALSE] # group is true cell type label
# colnames(covar_df) <- c('donors')
# 
# # compute cell library sizes to use as a covariate in the mixed models
# lib_sizes <- colSums(counts(sim_select))
# 
# # log transform and scale lib_sizes
# lib_sizes_scaled <- scale(log(lib_sizes[rownames(covar_df)]))[,1]
# covar_df$lib_sizes <- lib_sizes_scaled
# 
# true_ct <- sim_select@colData[,c('Group'),drop=FALSE]
# colnames(true_ct) <- 'cell_type'
# 
# saveRDS(list(tmp_norm,covar_df,pcs,donor_geno_eqtl,main_tr,ref_matrix,chr,true_ct,causal_snp,eQTL_group,eQTL_gene),
#         file='/home/jmitchel/scJLIM/data/sim_data.rds',compress = 'xz')
# 
# 
# 
# 
# 
# 
# 
# 
# tmp <- cbind.data.frame(covar_df,true_ct,pcs)
# p1 <- ggplot(tmp,aes(x=PC1,y=PC2,color=cell_type)) +
#   geom_point(alpha=.3) +
#   ggtitle('True cell type labels') +
#   theme_bw() +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# p2 <- ggplot(tmp,aes(x=PC1,y=PC2,color=donors)) +
#   geom_point(alpha=.3) +
#   ggtitle('Individuals') +
#   theme_bw() +
#   theme(legend.position = "none",plot.title = element_text(hjust = 0.5))
# plot_grid(p1,p2,nrow=1,rel_widths=c(.57,.45))
# 
# 
# # tmp <- cbind.data.frame(covar_df,true_ct,ures[["layout"]])
# # colnames(tmp)[4:5] <- c('UMAP1','UMAP2')
# # p1 <- ggplot(tmp,aes(x=UMAP1,y=UMAP2,color=cell_type)) +
# #   geom_point(alpha=.3) +
# #   ggtitle('True cell type labels') +
# #   theme_bw() +
# #   theme(plot.title = element_text(hjust = 0.5))
# # 
# # p2 <- ggplot(tmp,aes(x=UMAP1,y=UMAP2,color=donors)) +
# #   geom_point(alpha=.3) +
# #   ggtitle('Individuals') +
# #   theme_bw() +
# #   theme(legend.position = "none",plot.title = element_text(hjust = 0.5))
# # plot_grid(p1,p2,nrow=1,rel_widths=c(.57,.45))














base_dir <- '/home/jmitchel/data/GWAS_sim_v3'

loci_coords <- read.csv(paste0(base_dir,'/loci_coords.csv'),row.names = 1,skip = 1, header = FALSE)

sc_sims_all <- readRDS(file=paste0(base_dir,'/splat_sims/sc_sims_all6.rds'))
# sc_sims_all <- readRDS(file=paste0(base_dir,'/splat_sims/sc_sims_all7.rds'))

sims_type <- 1
sims_subtype <- 1

sc_sims_sub <- lapply(1:length(sc_sims_all),function(sim_ndx) {
  x <- sc_sims_all[[sim_ndx]]
  if (class(x)=='character') {
    return(NULL)
  }
  sims_iter <- x[[sims_type]]
  if (sims_type==3) {
    return(sims_iter)
  } else {
    sims_iter_sub <- list()
    for (gw_ndx in 1:length(sims_iter)) {
      sim_select <- sims_iter[[gw_ndx]][[sims_subtype]]
      # for null cases, need to save the name of the eQTL gene and group of the corresponding positive case
      if (sims_type==2) {
        if (sims_subtype==1) {
          corresponding_pos <- x[[1]][[gw_ndx]][[1]]
        } else if (sims_subtype==2) {
          corresponding_pos <- x[[1]][[gw_ndx]][[4]]
        } else {
          stop('sims_subtype can only be c(1,2)')
        }
        # get the causal snp
        causal_snp <- read.table(file=paste0('/home/jmitchel/data/GWAS_sim_v3/hapgen_gwas/sim_',sim_ndx,'/causal_rh_snps.txt'),header=FALSE)
        causal_snp <- strsplit(as.character(causal_snp),split=',')[[1]][[1]]
        eQTL_ndx <- which(corresponding_pos@rowRanges@elementMetadata$eSNP.ID==causal_snp)
        eQTL_gene <- corresponding_pos@rowRanges@partitioning@NAMES[eQTL_ndx][1]
        eQTL_gene_split <- strsplit(eQTL_gene,split='_')[[1]]
        eQTL_gene <- paste0(eQTL_gene_split[1],'-',eQTL_gene_split[2])
        eQTL_group <- corresponding_pos@rowRanges@elementMetadata$eQTL.group[eQTL_ndx][1]
        sim_select$eQTL_gene <- eQTL_gene
        sim_select$eQTL_group <- eQTL_group
      }
      
      sims_iter_sub[[gw_ndx]] <- sim_select
      
    } 
    return(sims_iter_sub)
  }
})



n_g_sims <- 100

### splitting sims up into batches for more memory efficient parallelization
n_batches <- 4
by_len <- n_g_sims/n_batches
batch_start_ndx <- seq(1, (n_g_sims+1), by = by_len)
sc_sims_sub_batches <- list()
for (i in 1:(length(batch_start_ndx)-1)) {
  batch_ndx <- batch_start_ndx[i]:(batch_start_ndx[i+1]-1)
  sims_subset <- sc_sims_sub[batch_ndx]
  names(sims_subset) <- paste0('sim_',batch_ndx)
  sc_sims_sub_batches[[i]] <- sims_subset
}

sc_sims_batch <- sc_sims_sub_batches[[1]]

batch_sim_ndx <- 2

sim_ndx <- strsplit(names(sc_sims_batch)[batch_sim_ndx],split='_')[[1]][[2]]

chr <- loci_coords[sim_ndx,1]

# load the full eqtl genotypes matrix
vcf_file <- paste0(base_dir,'/hapgen_gwas/sim_',sim_ndx,'/eqtl_full.vcf')
donor_geno_eqtl <- readVcf(vcf_file, "hg19")
donor_geno_eqtl <- as.data.frame(t(as(genotypeToSnpMatrix(donor_geno_eqtl)$genotype, "numeric"))) 
snps_keep <- which(!is.na(donor_geno_eqtl[,1]))
donor_geno_eqtl <- donor_geno_eqtl[snps_keep,]

# keep snps with allele frequency over .1
af <- rowSums(donor_geno_eqtl)/(ncol(donor_geno_eqtl)*2)
maf <- pmin(af,1-af)
snps_keep <- names(maf)[maf>=.1]
donor_geno_eqtl <- donor_geno_eqtl[snps_keep,]

# load the GWAS sumstats
main_tr_lst <- readRDS(file=paste0(base_dir,'/hapgen_gwas/sim_',sim_ndx,'/gw_sumstats.rds'))

# load the reference panel for LD
ref_matrix <- readRDS(file=paste0(base_dir,'/hapgen_gwas/sim_',sim_ndx,'/ref_panel_jlim.rds'))

# get the causal snp
causal_snp <- read.table(file=paste0('/home/jmitchel/data/GWAS_sim_v3/hapgen_gwas/sim_',sim_ndx,'/causal_rh_snps.txt'),header=FALSE)
causal_snp <- strsplit(as.character(causal_snp),split=',')[[1]][[1]]

gw_ndx <- 1

main_tr <- main_tr_lst[[gw_ndx]]

rownames(main_tr) <- main_tr$all_snp_rsid
colnames(main_tr) <- c('SNP','pval','position','beta','varbeta','snp_ndx')

# reduce both datasets to the intersection of variants
snps_both <- intersect(rownames(donor_geno_eqtl),rownames(main_tr))
donor_geno_eqtl <- donor_geno_eqtl[snps_both,]
main_tr <- main_tr[snps_both,]
main_tr$CHR <- chr
main_tr <- main_tr[,c('CHR','position','pval')]
colnames(main_tr) <- c('CHR','BP','P')
sim.sc.gr2 <- sc_sims_batch[[batch_sim_ndx]][[gw_ndx]]

ne_ndx=NULL

# get the gene that the gwas variant eQTL corresponds to
if ('eQTL_gene' %in% colnames(colData(sim.sc.gr2))) { # for negative cases
  eQTL_gene <- colData(sim.sc.gr2)$eQTL_gene[1]
  eQTL_group <- colData(sim.sc.gr2)$eQTL_group[1]
} else {
  eQTL_ndx <- which(sim.sc.gr2@rowRanges@elementMetadata$eSNP.ID==causal_snp)
  if (length(eQTL_ndx)==0 | !is.null(ne_ndx)) {
    if (!is.null(ne_ndx)) { # for null extra cases
      eQTL_ndx <- which(!is.na(sim.sc.gr2@rowRanges@elementMetadata$eSNP.ID))[ne_ndx]
      eQTL_group <- sim.sc.gr2@rowRanges@elementMetadata$eQTL.group[eQTL_ndx]
    } else if ('eSNP.ID2' %in% colnames(sim.sc.gr2@rowRanges@elementMetadata)) { # otherwise, it's red herring w causal snp in second column
      eQTL_ndx <- which(sim.sc.gr2@rowRanges@elementMetadata$eSNP.ID2==causal_snp)
      eQTL_group <- sim.sc.gr2@rowRanges@elementMetadata$eQTL.group2[eQTL_ndx]
    } 
  } else {
    eQTL_ndx <- eQTL_ndx[1]
    eQTL_group <- sim.sc.gr2@rowRanges@elementMetadata$eQTL.group[eQTL_ndx]
  }
  eQTL_gene <- sim.sc.gr2@rowRanges@partitioning@NAMES[eQTL_ndx] # for new sims
  eQTL_gene_split <- strsplit(eQTL_gene,split='_')[[1]]
  eQTL_gene <- paste0(eQTL_gene_split[1],'-',eQTL_gene_split[2])
}

# make cell names unique
new_cellnames <- sapply(1:nrow(sim.sc.gr2@colData),function(x){
  paste0('cell',x)
})
rownames(sim.sc.gr2@colData) <- new_cellnames
sim.sc.gr2@colData$Cell <- new_cellnames
colnames(counts(sim.sc.gr2)) <- new_cellnames

# normalize counts
sim.sc.gr2 <- logNormCounts(sim.sc.gr2)



### temporarily testing adding some more variation between cells
tmp_norm <- sim.sc.gr2@assays@data@listData[["logcounts"]]
g_alter <- sample(rownames(tmp_norm),50)
cell_var_vals <- .2*rnorm(ncol(tmp_norm))
for (g in g_alter) {
  tmp_norm[g,] <- tmp_norm[g,] + cell_var_vals
}
sim.sc.gr2@assays@data@listData[["logcounts"]] <- tmp_norm
###

pcell_counts_norm <- sim.sc.gr2@assays@data@listData[["logcounts"]]
rownames(pcell_counts_norm) <- sapply(rownames(pcell_counts_norm),function(x){
  mysplit=strsplit(x,split='_')[[1]]
  return(paste0(mysplit[1],'-',mysplit[2]))
})
rownames(sim.sc.gr2) <- rownames(pcell_counts_norm)

lib_sizes <- colSums(counts(sim.sc.gr2))

n_PCs <- 2
sim.sc.gr2 <- runPCA(sim.sc.gr2, ncomponents = n_PCs, scale=TRUE)

pseudo_pcs <- reducedDims(sim.sc.gr2)@listData[["PCA"]]

covar_df <- sim.sc.gr2@colData[,'Sample',drop=FALSE]
colnames(covar_df)[1] <- c('donors')

# log transform and scale lib_sizes
lib_sizes_scaled <- scale(log(lib_sizes[rownames(covar_df)]))[,1]
covar_df$lib_sizes <- lib_sizes_scaled

# extract true cell type labels
true_ct <- sim.sc.gr2@colData[,'Group',drop=FALSE]

saveRDS(list(pcell_counts_norm,covar_df,pseudo_pcs,donor_geno_eqtl,main_tr,ref_matrix,true_ct,causal_snp,eQTL_group,eQTL_gene),
        file='/home/jmitchel/scJLIM/data/sim_data.rds',compress = 'xz')







# ### testing just using one i know worked before
# sim.sc.gr2 <- readRDS(file='/home/jmitchel/tmp_testdat.rds')
# 
# pcell_counts_norm <- sim.sc.gr2@assays@data@listData[["logcounts"]]
# rownames(pcell_counts_norm) <- sapply(rownames(pcell_counts_norm),function(x){
#   mysplit=strsplit(x,split='_')[[1]]
#   return(paste0(mysplit[1],'-',mysplit[2]))
# })
# rownames(sim.sc.gr2) <- rownames(pcell_counts_norm)
# lib_sizes <- colSums(counts(sim.sc.gr2))
# pseudo_pcs <- reducedDims(sim.sc.gr2)@listData[["PCA"]]
# 
# covar_df <- sim.sc.gr2@colData[,'Sample',drop=FALSE]
# colnames(covar_df)[1] <- c('donors')
# 
# # log transform and scale lib_sizes
# lib_sizes_scaled <- scale(log(lib_sizes[rownames(covar_df)]))[,1]
# covar_df$lib_sizes <- lib_sizes_scaled
# 
# # extract true cell type labels
# true_ct <- sim.sc.gr2@colData[,'Group',drop=FALSE]
# 
# saveRDS(list(pcell_counts_norm,covar_df,pseudo_pcs,donor_geno_eqtl,main_tr,ref_matrix,true_ct,causal_snp,eQTL_group,eQTL_gene),
#         file='/home/jmitchel/scJLIM/data/sim_data.rds',compress = 'xz')








# # create formatted objects for scJLIM
# jlim_vars <- prep_jlim(main_tr,refLD.mat=ref_matrix,min.MAF=0.05)
# main_tr <- jlim_vars[[1]]
#
# # create a null distribution from permuted genotypes
# nperm <- 100000
# null_dist <- get_null_dist(jlim_vars,
#                            sectr.sample.size=ncol(donor_geno_eqtl),
#                            nperm=nperm,r2res=.8,
#                            n.cores=1)
#
#
# snp_res_mat <- get_eQTL_res(gene_test = eQTL_gene, norm_counts=pcell_counts_norm,
#                             cell_meta=covar_df, cell_pcs=pseudo_pcs[,1:n_PCs], n_PCs,
#                             geno_mat=donor_geno_eqtl, main_tr=main_tr,n.cores=1,
#                             use_ivw=FALSE,progress=FALSE,
#                             julia_dir="/home/jmitchel/.julia/juliaup/julia-1.10.4+0.x64.linux.gnu/bin")
#
# jlim_res1 <- jlim_main(snp_res_mat, jlim_vars, null_dist, sectr.sample.size=ncol(donor_geno_eqtl),
#                        min.SNPs.count=15, n.cores=1,global_adjust_method='cauchy',
#                        adjust_method='none',cluster_pvals=FALSE,n_eff=NULL,
#                        progress=TRUE,censor_cells=TRUE,censor_type='min_pval')
#
#
# ## testing reshuffling high pvals
# per_cell_pvals <- jlim_res1[[2]]
# per_cell_pvals[per_cell_pvals>.5] <- runif(sum(per_cell_pvals>.5),min=.5,max=1)
# per_cell_pvals[per_cell_pvals==0] <- 1/length(null_dist[[1]])
# cauchy_all <- ACAT(per_cell_pvals)



















# #### trying to debug why my sims look weird
# suppressPackageStartupMessages({
#   library(devtools)
#   load_all('/home/jmitchel/scJLIM/')
#   library(SingleCellExperiment)
#   library(scater)
#   library(cowplot)
#   library(parallel)
#   library(ACAT)
#   library(viridis)
#   library(ggrastr)
# })
# 
# # first need to setup juliacall package for faster model fits
# library(JuliaCall)
# julia <- julia_setup(verbose = TRUE)
# julia_install_package_if_needed("MixedModels")
# julia_install_package_if_needed("RCall")
# julia$library("MixedModels")
# julia$library("RCall")
# 
# julia_dir="/home/jmitchel/.julia/juliaup/julia-1.10.4+0.x64.linux.gnu/bin"
# 
# # load small simulated dataset with a positive colocalization in one of two cell types
# sim_dat_all <- readRDS(file='/home/jmitchel/scJLIM/data/sim_data_old.rds')
# 
# counts_norm <- sim_dat_all[[1]] # scRNA seq normalized and log transformed counts
# covar_df <- sim_dat_all[[2]] # scRNA seq metadata with donor and library size covariates (can contain more)
# pcs <- sim_dat_all[[3]] # scRNA seq expression PCs
# donor_geno_eqtl <- sim_dat_all[[4]] # scRNA donor genotypes (maf > 0.1, but can be different)
# main_tr <- sim_dat_all[[5]] # GWAS summary statistics
# ref_matrix <- sim_dat_all[[6]] # reference LD matrix
# true_ct <- sim_dat_all[[7]] # the preprogrammed cell type labels (unknown in real data)
# causal_snp <- sim_dat_all[[8]] # the preprogrammed causal GWAS and eQTL SNP (unknown in real data)
# eQTL_group <- sim_dat_all[[9]] # the preprogrammed cell type with the eQTL (unknown in real data)
# eQTL_gene <- sim_dat_all[[10]] # the preprogrammed eGene (unknown in real data)
# 
# # create formatted objects for scJLIM
# jlim_vars <- prep_jlim(main_tr,refLD.mat=ref_matrix,min.MAF=0.05)
# main_tr <- jlim_vars[[1]]
# 
# # create a null distribution from permuted genotypes
# nperm <- 100000
# null_dist <- get_null_dist(jlim_vars,
#                            sectr.sample.size=ncol(donor_geno_eqtl),
#                            nperm=nperm,r2res=.8,
#                            n.cores=1)
# 
# # fit mixed models for each SNP in the locus
# n_PCs <- 2 # number of PC interaction terms to use (must be 2 or greater)
# 
# # snp_res_mat <- get_eQTL_res(gene_test = eQTL_gene, norm_counts=counts_norm,
# #                             cell_meta=covar_df, cell_pcs=pcs, n_PCs,
# #                             geno_mat=donor_geno_eqtl, main_tr=main_tr,n.cores=1,
# #                             progress=FALSE,
# #                             julia_dir=julia_dir)
# # compute colocalization p-values per cell and a cauchy (global) p-value over all cells
# # jlim_res1 <- jlim_main(snp_res_mat, jlim_vars, null_dist, sectr.sample.size=ncol(donor_geno_eqtl),
# #                        min.SNPs.count=15, n.cores=1)
# ## extract results
# # cauchy_all <- jlim_res1[[1]] # global p-value over all cells
# # per_cell_pvals <- jlim_res1[[2]] # p-values for all individual cells
# 
# 
# snp_res_mat <- get_eQTL_res(gene_test = eQTL_gene, norm_counts=counts_norm,
#                             cell_meta=covar_df, cell_pcs=pcs[,1:n_PCs], n_PCs,
#                             geno_mat=donor_geno_eqtl, main_tr=main_tr,n.cores=1,
#                             use_ivw=FALSE,progress=FALSE,
#                             julia_dir="/home/jmitchel/.julia/juliaup/julia-1.10.4+0.x64.linux.gnu/bin")
# 
# jlim_res1 <- jlim_main(snp_res_mat, jlim_vars, null_dist, sectr.sample.size=ncol(donor_geno_eqtl),
#                        min.SNPs.count=15, n.cores=1,global_adjust_method='cauchy',
#                        adjust_method='none',cluster_pvals=FALSE,n_eff=NULL,
#                        progress=TRUE,censor_cells=TRUE,censor_type='min_pval')
# 
# 
# ## testing reshuffling high pvals
# per_cell_pvals <- jlim_res1[[2]]
# per_cell_pvals[per_cell_pvals>.5] <- runif(sum(per_cell_pvals>.5),min=.5,max=1)
# per_cell_pvals[per_cell_pvals==0] <- 1/length(null_dist[[1]])
# cauchy_all <- ACAT(per_cell_pvals)
# 
# 
# 
# 
# 
# # looking at the scRNA-seq cell population and donor structure
# tmp <- cbind.data.frame(covar_df,true_ct,pcs)
# tmp2 <- cbind.data.frame(tmp,per_cell_pvals)
# tmp2$true_coloc_cells <- tmp2$Group==eQTL_group
# 
# p1 <- ggplot(tmp2,aes(x=PC1,y=PC2,color=true_coloc_cells)) + 
#   geom_point_rast(size = .5) +
#   ggtitle('Ground truth') +
#   theme_bw() +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   guides(color=guide_legend(title="Has coloc."))
# 
# p2 <- ggplot(tmp,aes(x=PC1,y=PC2,color=-log10(per_cell_pvals))) +
#   geom_point_rast(size = .5) +
#   scale_colour_gradientn(
#     colors = viridis(256,option = "plasma"),
#     values = scales::rescale(c(0,.15,.35,.5,1)),  # normalize your breakpoints
#     limits = c(0, 5)) +
#   theme_bw() +
#   labs(title = 'scJLIM per-cell pvalues',
#        color = "-log10(p)") +
#   theme(plot.title = element_text(hjust = 0.5),
#         legend.spacing.y = unit(0.1, "cm"),   # reduces vertical spacing between legend items
#         legend.key.height = unit(0.4, "cm"),    # decreases the height of the legend key
#         legend.key.width = unit(0.6, "cm"))
# 
# plot_grid(p1,p2,nrow=1,rel_widths=c(.5,.5))
# 
# 
# 
# 
# 
# 
# 
# ### picking out some cells with low pc1 and high pc2
# print(causal_snp)
# pc_sub <- pcs[pcs[,'PC1']<0&pcs[,'PC2']>0,]
# head(pc_sub)
# cell_select <- rownames(pc_sub)[4]
# 
# # looking at the simulated GWAS locus
# sec_tr <- main_tr
# sec_tr$P <- snp_res_mat[[cell_select]]
# sec_tr$causal_snp <- rownames(sec_tr)==causal_snp
# 
# p <- ggplot(sec_tr,aes(x=BP,y=-log10(P),color=causal_snp,size=causal_snp)) +
#   geom_point(data=sec_tr[which(!sec_tr$causal_snp),],alpha=.7) +
#   geom_point(data=sec_tr[which(sec_tr$causal_snp),],alpha=1) +
#   scale_size_manual(values = c("TRUE" = 3, "FALSE"=1, 1)) +
#   scale_color_manual(values = c("blue2", "red1")) +
#   ggtitle('eQTLs for a cell') +
#   theme_bw() +
#   theme(plot.title = element_text(hjust = 0.5)) + 
#   guides(color=guide_legend(title="Causal SNP"),size=guide_legend(title="Causal SNP"))
# p
















## saving one example simulation for use with the tutorial

library(VariantAnnotation)
library(scater)
library(cowplot)
library(umap)

base_dir <- '/home/jmitchel/data/GWAS_sim_v3'

loci_coords <- read.csv(paste0(base_dir,'/loci_coords.csv'),row.names = 1,skip = 1, header = FALSE)

sc_sims_all <- readRDS(file=paste0(base_dir,'/splat_sims/sc_sims_all6.rds'))

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


### adding some more variation between cells so eQTL isn't a PC itself
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

# saveRDS(list(pcell_counts_norm,covar_df,pseudo_pcs,donor_geno_eqtl,main_tr,ref_matrix,true_ct,causal_snp,eQTL_group,eQTL_gene),
#         file='/home/jmitchel/scJLIM/data/sim_data.rds',compress = 'xz')



















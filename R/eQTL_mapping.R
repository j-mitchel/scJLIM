
lme_helper <- function(data,geno_mat,snp_name,n_PCs) {
  
  # need to add genotype of lead variant to a column
  gt <- unlist(geno_mat[snp_name,])
  data <- cbind.data.frame(data,gt[as.character(data[,'donors'])]) # needs to be character because indexing doesnt work if it's a factor
  colnames(data)[ncol(data)] <- 'geno'
  
  if ('batch' %in% colnames(data)) {
    base_eqn <- paste0('expr ~ geno + (1 | batch) + (1 | donors)')
  } else {
    base_eqn <- paste0('expr ~ geno + (1 | donors)')
  }
  
  covars_vec <- colnames(data)
  covars_vec <- covars_vec[!(covars_vec %in% c('expr','donors','batch','geno'))]
  covar_eqn <- paste0(covars_vec, collapse = " + ")
  pc_g_intr_eqn <- paste0("geno:PC", 1:n_PCs, collapse = " + ")
  
  # create formula with all variables
  f_mod_vars <- paste(base_eqn,covar_eqn,pc_g_intr_eqn,sep = " + ")
  f_mod_form <- as.formula(f_mod_vars)
  
  full_model <- lmerTest::lmer(formula=f_mod_form, data=data)
  coef_mat <- summary(full_model)$coefficients[c('geno',paste0("geno:PC", 1:n_PCs)),c('Estimate','Std. Error','Pr(>|t|)')]
  
  return(coef_mat)
}

get_per_cell_pv <- function(data,coef_mat,n_PCs,use_ivw=TRUE) {
  intr_pc_names <- paste0("geno:PC", 1:n_PCs)
  pc_names <- paste0("PC", 1:n_PCs)
  
  coefs <- coef_mat[intr_pc_names,'Estimate']
  std_err <- coef_mat[intr_pc_names,'Std. Error']
  
  pcs_select <- as.matrix(data[,pc_names,drop=FALSE])
  
  if (use_ivw) {
    ### computing inverse variance weighted average genotypic effects per cell
    
    # expanding coefficient values by pc value per cell
    expanded_coefs <- pcs_select %*% diag(coefs) # dimensions are cells x pc_interaction terms
    
    # append the persistent genotypic effect term as a column
    expanded_coefs <- cbind(expanded_coefs,coef_mat['geno','Estimate'])
    
    # now weight the coefs by their corresponding variance
    variances <- c(std_err^2,coef_mat['geno','Std. Error']^2)
    coefs_weighted <- sweep(expanded_coefs,2,variances,'/')
    
    ### The "full" way to calculate Z-scores
    # # now get the average genotypic effect across all terms involving genotype
    # total_weighted_effects <- rowSums(coefs_weighted)
    # ivw_effects <- total_weighted_effects / sum(1/variances)
    # 
    # ### now calculating standard errors per cell
    # # append persistent effect error
    # numerator_squared <- rowSums(sweep(pcs_select^2,2,std_err^2,'/'))
    # numerator_squared <- numerator_squared + (1/(coef_mat['geno','Std. Error']^2))
    # numerator <- sqrt(numerator_squared)
    # denominator <- sum(1/variances)
    # ivw_error <- numerator/denominator
    # zsc <- ivw_effects/ivw_error
    
    ### the slightly simplified way of calculating Z-scores which gives identical results
    total_weighted_effects <- rowSums(coefs_weighted)
    se_numerator_squared <- rowSums(sweep(pcs_select^2,2,std_err^2,'/'))
    se_numerator_squared <- se_numerator_squared + (1/(coef_mat['geno','Std. Error']^2))
    se_numerator <- sqrt(se_numerator_squared)
    zsc <- total_weighted_effects / se_numerator
    pval <- pnorm(abs(zsc), mean = 0, sd = 1, lower.tail = FALSE) * 2
  } else {
    expanded_coefs <- pcs_select %*% diag(coefs)
    expanded_errors <- pcs_select %*% diag(std_err) # sd gets multiplied by the constant
    total_effect <- rowSums(expanded_coefs) + coef_mat['geno','Estimate']
    total_effect_error <- sqrt(rowSums(expanded_errors^2) + (coef_mat['geno','Std. Error']^2)) # variances added across terms
    zsc <- total_effect / total_effect_error
    pval <- pnorm(abs(zsc), mean = 0, sd = 1, lower.tail = FALSE) * 2
  }
  
  return(pval)
}


#' Main function to run eQTL mapping
#' 
#' @param gene_test character The gene to test for cis-eQTLs with
#' @param norm_counts matrix Normalized UMI counts matrix (genes x cells)
#' @param cell_meta data.frame Cell-level metadata. Columns must include 'donors', as 
#' the model includes donor random effects. Rownames must be cell names. Optionally include a 
#' 'batch' column to also incorporate batch random effects. Include all other covariates to be 
#' used in the model. We suggest pre-scaling and log transforming variables with long tails such as library sizes.
#' @param cell_pcs matrix Cell by PC matrix, used for genotype-interaction testing. Columns
#' should be labeled PC1, PC2, etc. Rownames should be cell names.
#' @param n_PCs numeric The number of cell_pcs to use for interaction testing
#' @param geno_mat matrix Imputed genotype matrix (SNPs by donors), with values representing
#' the number of alternate alleles per donor (0, 1, or 2)
#' @param main_tr data.frame The GWAS locus data with columns as 'CHR','BP','P', corresponding
#' to chromosome, position, and p-value
#' @param geno_pcs matrix Optional donor by PC matrix to include as covariates in the model. 
#' Columns should be labeled as geno_PC1, geno_PC2, etc.
#'
#' @return
#' @export
#'
#' @examples
get_eQTL_res <- function(gene_test,norm_counts,cell_meta,cell_pcs,geno_mat,main_tr,geno_pcs=NULL,n.cores=4) {
  # check everything is the right dimensions
  cell_names <- colnames(norm_counts)
  if (!all(cell_names %in% rownames(cell_pcs))) {
    stop('Not all cells are in the PC matrix')
  } else if (!all(cell_names %in% rownames(cell_meta))) {
    stop('Not all cells are in the cell_meta data.frame')
  } else if (!(gene_test %in% rownames(norm_counts))) {
    stop('gene_test is not in rownames of norm_counts')
  }
  
  # check that a 'donors' column exists in cell_meta
  if (!('donors' %in% colnames(cell_meta))) {
    stop('cell_meta ust a column called donors')
  }
  
  # check that pc names are specified correctly
  expected_pc_names <- paste0('PC',1:ncol(cell_pcs))
  if (!all.equal(expected_pc_names,colnames(cell_pcs))) {
    stop('Columns of cell_pcs should be named PC1, PC2, ...')
  }
  
  # check that all donors are represented in geno_mat
  donors_all <- unique(cell_meta$donors)
  if (!all(donors_all %in% colnames(geno_mat))) {
    stop('Not all donors are in geno_mat')
  }
  
  # first make a matrix with all covariates to be reused for mapping each snp
  data <- cbind.data.frame(norm_counts[gene_test,],
                                cell_pcs[cell_names,],
                                cell_meta[cell_names,])
  colnames(data)[1] <- 'expr'
  
  if (!is.null(geno_pcs)) {
    expected_geno_pc_names <- paste0('geno_PC',1:ncol(geno_pcs))
    if (!all.equal(expected_geno_pc_names,colnames(geno_pcs))) {
      stop('Columns of geno_pcs should be named geno_PC1, geno_PC2, ...')
    }
    
    ## append genotype PC covariates
    data <- cbind.data.frame(data,geno_pcs[as.character(data[,'donors']),])
  }
  
  ### quantile normalize expression to gaussian
  norm_expr <- normalize.quantiles.use.target(as.matrix(data$expr),rnorm(10000),copy=TRUE,subset=NULL)[,1]
  data$expr <- norm_expr
  
  # run eQTL mapping for each SNP and get total effects and errors
  snp_res <- plapply(1:nrow(main_tr),function(snp_j) tryCatch({
    snp_name <- rownames(main_tr)[snp_j]
    coef_mat <- lme_helper(data,geno_mat,snp_name,n_PCs)
    per_cell_pvals <- get_per_cell_pv(data,coef_mat,n_PCs)
    return(per_cell_pvals)
  },error=function(e) paste0('error_index_',snp_j)),progress = TRUE,n.cores = n.cores,mc.preschedule = TRUE)
  
  snp_res_mat <- do.call(rbind.data.frame, snp_res)
  colnames(snp_res_mat) <- names(snp_res[[1]])
  rownames(snp_res_mat) <- rownames(main_tr)
  return(snp_res_mat)
}





lme_helper <- function(data_in,geno_mat,snp_name,n_PCs) {

  # need to add genotype of lead variant to a column
  gt <- unlist(geno_mat[snp_name,])
  data_in <- cbind.data.frame(data_in,gt[as.character(data_in[,'donors'])]) # needs to be character because indexing doesnt work if it's a factor
  colnames(data_in)[ncol(data_in)] <- 'geno'
  
  if ('batch' %in% colnames(data_in)) {
    base_eqn <- paste0('expr ~ geno + (1 | batch) + (1 | donors)')
  } else {
    base_eqn <- paste0('expr ~ geno + (1 | donors)')
  }
  
  covars_vec <- colnames(data_in)
  covars_vec <- covars_vec[!(covars_vec %in% c('expr','donors','batch','geno'))]
  covar_eqn <- paste0(covars_vec, collapse = " + ")
  pc_g_intr_eqn <- paste0("geno:PC", 1:n_PCs, collapse = " + ")
  
  # create formula with all variables
  f_mod_vars <- paste(base_eqn,covar_eqn,pc_g_intr_eqn,sep = " + ")
  f_mod_form <- as.formula(f_mod_vars)
  
  
  ## trying mixedmodels in julia to see if it's any faster
  julia$assign("lmm_data_in", data_in)
  julia$assign("form", f_mod_form)
  full_model <- julia$eval("full_model = fit(LinearMixedModel, form, lmm_data_in)",need_return = c("Julia"))
  
  beta.fixed <- julia_eval("coef(full_model)")
  std_err <- julia_eval("full_model.stderror")
  pvals <- julia_eval("full_model.pvalues")
  vcov_mat <- julia_eval("vcov(full_model)")
  
  coef_mat <- cbind.data.frame(beta.fixed, std_err, pvals)
  colnames(coef_mat) <- c('Estimate','Std. Error','Pr(>|t|)')
  fixed_coef_nms <- julia_eval("fixefnames(full_model)")
  fixed_coef_nms[(length(fixed_coef_nms)-n_PCs+1):length(fixed_coef_nms)] <- paste0("geno:PC", 1:n_PCs)
  rownames(coef_mat) <- fixed_coef_nms
  coef_mat <- coef_mat[c('geno',paste0("geno:PC", 1:n_PCs)),]

  colnames(vcov_mat) <- fixed_coef_nms
  rownames(vcov_mat) <- fixed_coef_nms
  
  # full_model <- lmerTest::lmer(formula=f_mod_form, data=data_in, control = lmerControl(calc.derivs = FALSE))
  # coef_mat <- summary(full_model)$coefficients[c('geno',paste0("geno:PC", 1:n_PCs)),c('Estimate','Std. Error','Pr(>|t|)')]
  # vcov_mat <- vcov(full_model)
  return(list(coef_mat,vcov_mat))
}

get_per_cell_pv <- function(data_in,coef_mat,n_PCs,use_ivw=TRUE,ivw_type='orig_var') {
  intr_pc_names <- paste0("geno:PC", 1:n_PCs)
  pc_names <- paste0("PC", 1:n_PCs)
  
  coefs <- coef_mat[intr_pc_names,'Estimate']
  std_err <- coef_mat[intr_pc_names,'Std. Error']
  
  pcs_select <- as.matrix(data_in[,pc_names,drop=FALSE])
  
  if (use_ivw) {
    if (ivw_type=='orig_var') {
      ##### trying to inverse variance weight only the interaction terms and
      # then take the average of this with the persistent term
      expanded_coefs <- pcs_select %*% diag(coefs) # dimensions are cells x pc_interaction terms
      variances <- c(std_err^2)
      coefs_weighted <- sweep(expanded_coefs,2,variances,'/')
      total_weighted_effects <- rowSums(coefs_weighted)
      ivw_effects <- total_weighted_effects / sum(1/variances)

      # now combine the average PC modification to the persistent effect
      estimated_effect <- coef_mat['geno','Estimate'] + ivw_effects

      # now calculate variance of these estimated effects
      numerator <- rowSums(sweep(pcs_select^2,2,std_err^2,'/'))
      denominator <- sum(1/variances)^2
      pc_var <- numerator / denominator
      total_var <- (coef_mat['geno','Std. Error']^2) + pc_var
      total_se <- sqrt(total_var)
      zsc <- estimated_effect / total_se
      pval <- pnorm(abs(zsc), mean = 0, sd = 1, lower.tail = FALSE) * 2

    } else if (ivw_type=='pc_expanded_var') {

      ##### trying the interaction term weighting approach but with pc expansions
      expanded_coefs <- pcs_select %*% diag(coefs) # dimensions are cells x pc_interaction terms
      expanded_errors <- pcs_select %*% diag(std_err)
      coefs_weighted <- expanded_coefs / (expanded_errors^2)
      total_weighted_effects <- rowSums(coefs_weighted)
      ivw_effects <- total_weighted_effects / rowSums(1/(expanded_errors^2))

      # now combine the average PC modification to the persistent effect
      estimated_effect <- coef_mat['geno','Estimate'] + ivw_effects

      weighted_variances <- (expanded_errors^2) / ((expanded_errors^2)^2)
      av_pc_weighted_variances <- rowSums(weighted_variances) / (rowSums(1/(expanded_errors^2))^2)
      total_var <- (coef_mat['geno','Std. Error']^2) + av_pc_weighted_variances
      total_se <- sqrt(total_var)
      zsc <- estimated_effect / total_se
      pval <- pnorm(abs(zsc), mean = 0, sd = 1, lower.tail = FALSE) * 2
      
    } else if (ivw_type=='z_var') {
      
      ### trying with the square root of se instead of variance
      expanded_coefs <- pcs_select %*% diag(coefs) # dimensions are cells x pc_interaction terms
      
      # append the persistent genotypic effect term as a column
      expanded_coefs <- cbind(coef_mat['geno','Estimate'],expanded_coefs)
      
      # now weight the coefs by their corresponding se
      total_weighted_effects <- expanded_coefs %*% (1/sqrt(coef_mat[,'Std. Error']))
      
      ### the slightly simplified way of calculating Z-scores which gives identical results
      se_numerator_squared <- (pcs_select^2) %*% std_err
      se_numerator_squared <- se_numerator_squared + (1/sqrt(coef_mat['geno','Std. Error']))
      se_numerator <- sqrt(se_numerator_squared)
      zsc <- total_weighted_effects / se_numerator
      pval <- pnorm(abs(zsc), mean = 0, sd = 1, lower.tail = FALSE) * 2
      pval <- c(pval)
    } else {
      stop('The requested ivw_type not implemented')
    }
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


get_per_cell_pv_covar <- function(data_in,coef_mat,n_PCs,vcov_mat) {
  intr_pc_names <- paste0("geno:PC", 1:n_PCs)
  pc_names <- paste0("PC", 1:n_PCs)
  
  vcov_mat <- vcov_mat[c('geno',intr_pc_names),c('geno',intr_pc_names)]
  vcov_mat[lower.tri(vcov_mat)] <- 0
  diag(vcov_mat) <- 0
  vcov_mat <- vcov_mat[,intr_pc_names] # remove geno column
  
  coefs <- coef_mat[intr_pc_names,'Estimate']
  std_err <- coef_mat[intr_pc_names,'Std. Error']
  
  pcs_select <- as.matrix(data_in[,pc_names,drop=FALSE])
  
  expanded_coefs <- pcs_select %*% diag(coefs)
  expanded_errors <- pcs_select %*% diag(std_err) # sd gets multiplied by the constant
  
  # now computing the expanded variance with covariance terms using formula from case 3 of: https://mattgolder.com/wp-content/uploads/2015/05/standarderrors1.png
  covar_list <- list()
  for (myterm in c('geno',intr_pc_names)) {
    term_covars <- vcov_mat[myterm,]
    expanded_term_covars <- pcs_select %*% diag(term_covars)
    if (myterm=='geno') {
      expanded_term_covars <- 2 * expanded_term_covars
    } else {
      expanded_term_covars <- 2 * expanded_term_covars * pcs_select[,strsplit(myterm,split=':')[[1]][[2]]]
    }
    covar_list[[myterm]] <- expanded_term_covars
  }
  
  covar_list2 <- do.call(cbind,covar_list)
  
  total_effect <- rowSums(expanded_coefs) + coef_mat['geno','Estimate']
  total_effect_error <- sqrt(rowSums(expanded_errors^2) + (coef_mat['geno','Std. Error']^2) + rowSums(covar_list2)) # variances added across terms

  zsc <- total_effect / total_effect_error
  
  pval <- pnorm(abs(zsc), mean = 0, sd = 1, lower.tail = FALSE) * 2
  
  return(pval)
}


get_per_cell_pv_sig_only <- function(data_in,coef_mat,n_PCs,pc_ndx_keep) {
  intr_pc_names <- paste0("geno:PC", 1:n_PCs)
  pc_names <- paste0("PC", 1:n_PCs)
  
  intr_pc_names <- intr_pc_names[pc_ndx_keep]
  pc_names <- pc_names[pc_ndx_keep]
  
  if (length(pc_ndx_keep)==0) {
    zsc <- coef_mat['geno','Estimate'] / coef_mat['geno','Std. Error']
    pval <- pnorm(abs(zsc), mean = 0, sd = 1, lower.tail = FALSE) * 2
    pval <- rep(pval,nrow(data_in))
  } else {
    if (length(pc_ndx_keep)==1) {
      coefs <- coef_mat[intr_pc_names,'Estimate',drop=FALSE]
      std_err <- coef_mat[intr_pc_names,'Std. Error',drop=FALSE]
    } else {
      coefs <- coef_mat[intr_pc_names,'Estimate']
      std_err <- coef_mat[intr_pc_names,'Std. Error']
    }
    
    pcs_select <- as.matrix(data_in[,pc_names,drop=FALSE])
    expanded_coefs <- pcs_select %*% diag(coefs)
    expanded_errors <- pcs_select %*% diag(std_err) # sd gets multiplied by the constant
    total_effect <- rowSums(expanded_coefs) + coef_mat['geno','Estimate']
    total_effect_error <- sqrt(rowSums(expanded_errors^2) + (coef_mat['geno','Std. Error']^2)) # variances added across terms
    zsc <- total_effect / total_effect_error
    pval <- pnorm(abs(zsc), mean = 0, sd = 1, lower.tail = FALSE) * 2
  }

  return(pval)
}


get_per_cell_pv_covar_sig_only <- function(data_in,coef_mat,n_PCs,vcov_mat,pc_ndx_keep) {
  intr_pc_names <- paste0("geno:PC", 1:n_PCs)
  pc_names <- paste0("PC", 1:n_PCs)
  
  intr_pc_names <- intr_pc_names[pc_ndx_keep]
  pc_names <- pc_names[pc_ndx_keep]
  
  if (length(pc_ndx_keep)==0) {
    zsc <- coef_mat['geno','Estimate'] / coef_mat['geno','Std. Error']
    pval <- pnorm(abs(zsc), mean = 0, sd = 1, lower.tail = FALSE) * 2
    pval <- rep(pval,nrow(data_in))
    return(pval)
  }
    
  vcov_mat <- vcov_mat[c('geno',intr_pc_names),c('geno',intr_pc_names)]
  vcov_mat[lower.tri(vcov_mat)] <- 0
  diag(vcov_mat) <- 0
  vcov_mat <- vcov_mat[,intr_pc_names,drop=FALSE] # remove geno column
  
  if (length(pc_ndx_keep)==1) {
    coefs <- coef_mat[intr_pc_names,'Estimate',drop=FALSE]
    std_err <- coef_mat[intr_pc_names,'Std. Error',drop=FALSE]
  } else {
    coefs <- coef_mat[intr_pc_names,'Estimate']
    std_err <- coef_mat[intr_pc_names,'Std. Error']
  }
  
  pcs_select <- as.matrix(data_in[,pc_names,drop=FALSE])
  
  expanded_coefs <- pcs_select %*% diag(coefs)
  expanded_errors <- pcs_select %*% diag(std_err) # sd gets multiplied by the constant
  
  # now computing the expanded variance with covariance terms using formula from case 3 of: https://mattgolder.com/wp-content/uploads/2015/05/standarderrors1.png
  covar_list <- list()
  for (myterm in c('geno',intr_pc_names)) {
    term_covars <- vcov_mat[myterm,]
    if (length(intr_pc_names)==1) {
      expanded_term_covars <- pcs_select * term_covars
    } else {
      expanded_term_covars <- pcs_select %*% diag(term_covars)
    }
    if (myterm=='geno') {
      expanded_term_covars <- 2 * expanded_term_covars
    } else {
      expanded_term_covars <- 2 * expanded_term_covars * pcs_select[,strsplit(myterm,split=':')[[1]][[2]]]
    }
    covar_list[[myterm]] <- expanded_term_covars
  }
  
  covar_list2 <- do.call(cbind,covar_list)
  
  total_effect <- rowSums(expanded_coefs) + coef_mat['geno','Estimate']
  total_effect_error <- sqrt(rowSums(expanded_errors^2) + (coef_mat['geno','Std. Error']^2) + rowSums(covar_list2)) # variances added across terms
  
  zsc <- total_effect / total_effect_error
  
  pval <- pnorm(abs(zsc), mean = 0, sd = 1, lower.tail = FALSE) * 2
  
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
#' @param n.cores numeric Number of cores to use
#' @param use_ivw logical Set to TRUE to use inverse variance weighting
#' @param ivw_type character Set to 'orig_var' to use one over the interaction term variances as 
#' the weights. Set to 'pc_expanded_var' to use one over the PC expanded interaction term variances
#' as the weights.
#' 
#' @return
#' @export
#'
#' @examples
get_eQTL_res <- function(gene_test,norm_counts,cell_meta,cell_pcs,n_PCs,geno_mat,main_tr,
                         geno_pcs=NULL,n.cores=4,use_ivw=TRUE,ivw_type='orig_var',progress=TRUE) {
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
  
  ### scale expression
  scaled_expr <- scale(data[,'expr',drop=FALSE])
  data$expr <- c(scaled_expr)
  
  # check for any column names with '.' in it
  col_to_fix <- sapply(colnames(data),function(x){
    return(grepl( '.', x, fixed = TRUE))
  })
  col_to_fix <-  which(col_to_fix)
  new_col_nm <- sapply(colnames(data)[col_to_fix],function(x){
    new_splt <- strsplit(x,split='.',fixed=TRUE)[[1]]
    new_nm <- paste0(new_splt[[1]],'_',new_splt[[2]])
    return(new_nm)
  })
  colnames(data)[col_to_fix] <- new_col_nm
  
  # # run eQTL mapping for each SNP and get total effects and errors
  # snp_res <- plapply(1:nrow(main_tr),function(snp_j) tryCatch({
  # # snp_res <- lapply(1:nrow(main_tr),function(snp_j) {
  #   snp_name <- rownames(main_tr)[snp_j]
  #   lme_out <- lme_helper(data,geno_mat,snp_name,n_PCs)
  #   coef_mat <- lme_out[[1]]
  #   vcov_mat <- lme_out[[2]]
  #   # per_cell_pvals <- get_per_cell_pv(data,coef_mat,n_PCs,use_ivw=use_ivw,ivw_type=ivw_type)
  #   per_cell_pvals <- get_per_cell_pv_covar(data,coef_mat,n_PCs,vcov_mat)
  #   return(per_cell_pvals)
  # # })
  # },error=function(e) paste0('error_index_',snp_j)),progress = progress,n.cores = n.cores,mc.preschedule = TRUE)
  
  # ####### testing only 
  # plan(multicore, workers = 30)
  # plan(multisession, workers = 10)
  # snp_res <- future_lapply(1:50,function(snp_j) {
  #   julia$library("MixedModels")
  #   snp_name <- rownames(main_tr)[snp_j]
  #   lme_out <- lme_helper(data,geno_mat,snp_name,n_PCs)
  #   return(lme_out)
  # }, future.seed = 1)
  # 
  # snp_res <- plapply(1:50,function(snp_j) tryCatch({
  #   # julia$library("MixedModels")
  #   snp_name <- rownames(main_tr)[snp_j]
  #   lme_out <- lme_helper(data,geno_mat,snp_name,n_PCs)
  #   return(lme_out)
  # },error=function(e) paste0('error_index_',snp_j)),progress = progress,n.cores = 30,mc.preschedule = TRUE)
  # ########
  
  snp_res <- plapply(1:nrow(main_tr),function(snp_j) tryCatch({
    snp_name <- rownames(main_tr)[snp_j]
    lme_out <- lme_helper(data,geno_mat,snp_name,n_PCs)
    return(lme_out)
  },error=function(e) paste0('error_index_',snp_j)),progress = progress,n.cores = n.cores,mc.preschedule = TRUE)
  
  ## testing using acat here
  sig_vals <- get_PC_sig(snp_res)
  pc_pv <- c()
  for (i in 2:ncol(sig_vals)) {
    per_sig <- sig_vals[,i]
    per_sig_locus <- ACAT(per_sig)
    pc_pv <- c(pc_pv,per_sig_locus)
  }
  
  snp_res <- plapply(1:length(snp_res),function(i) tryCatch({
    coef_mat <- snp_res[[i]][[1]]
    vcov_mat <- snp_res[[i]][[2]]
    per_cell_pvals <- get_per_cell_pv_covar(data,coef_mat,n_PCs,vcov_mat)
    return(per_cell_pvals)
  },error=function(e) paste0('error_index_',snp_j)),progress = progress,n.cores = n.cores,mc.preschedule = TRUE)
  
  snp_res_mat <- do.call(rbind.data.frame, snp_res)
  colnames(snp_res_mat) <- names(snp_res[[1]])
  rownames(snp_res_mat) <- rownames(main_tr)
  return(list(snp_res_mat,pc_pv))
}



get_eQTL_res_sig_any <- function(gene_test,norm_counts,cell_meta,cell_pcs,n_PCs,geno_mat,main_tr,
                         geno_pcs=NULL,n.cores=4,progress=TRUE) {
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
  
  # ### scale expression
  scaled_expr <- scale(data[,'expr',drop=FALSE])
  data$expr <- c(scaled_expr)
  
  snp_res <- plapply(1:nrow(main_tr),function(snp_j) tryCatch({
    snp_name <- rownames(main_tr)[snp_j]
    lme_out <- lme_helper(data,geno_mat,snp_name,n_PCs)
    return(lme_out)
  },error=function(e) paste0('error_index_',snp_j)),progress = progress,n.cores = n.cores,mc.preschedule = TRUE)

  ## testing using acat here
  sig_vals <- get_PC_sig(snp_res)
  pc_pv <- c()
  for (i in 2:ncol(sig_vals)) {
    per_sig <- sig_vals[,i]
    per_sig_locus <- ACAT(per_sig)
    pc_pv <- c(pc_pv,per_sig_locus)
  }
  # pc_pv <- p.adjust(pc_pv,method='fdr')
  # pc_ndx_keep <- which(pc_pv<.05) # keep interaction terms that are significant
  pc_ndx_keep <- which(pc_pv<.05) # keep interaction terms that are significant
  
  
  ## testing not running the test if persistent term isn't significant across the locus
  # sig_vals <- get_PC_sig(snp_res)
  # per_sig <- sig_vals[,1]
  # per_sig_locus <- ACAT(per_sig)
  # if (per_sig_locus >.05) {
  #   return(NA)
  # }

  snp_res <- plapply(1:length(snp_res),function(i) tryCatch({
    coef_mat <- snp_res[[i]][[1]]
    vcov_mat <- snp_res[[i]][[2]]
    # per_cell_pvals <- get_per_cell_pv_covar(data,coef_mat,n_PCs,vcov_mat)
    per_cell_pvals <- get_per_cell_pv_covar_sig_only(data,coef_mat,n_PCs,vcov_mat,pc_ndx_keep)
    return(per_cell_pvals)
  },error=function(e) paste0('error_index_',snp_j)),progress = progress,n.cores = n.cores,mc.preschedule = TRUE)
  
  snp_res_mat <- do.call(rbind.data.frame, snp_res)
  colnames(snp_res_mat) <- names(snp_res[[1]])
  rownames(snp_res_mat) <- rownames(main_tr)
  return(snp_res_mat)
}

# get PC interaction term max significance
get_PC_max_sig <- function(snp_res) {
  sig_vals_all <- list()
  for (i in 1:length(snp_res)) {
    sig_vals <- snp_res[[i]][,'Pr(>|t|)']
    sig_vals_all[[i]] <- sig_vals
  }
  df <- do.call("rbind.data.frame", sig_vals_all)
  min_pvals <- colMins(as.matrix(df),useNames = FALSE)
  return(min_pvals)
}

# get PC interaction term max significance
get_PC_sig <- function(snp_res) {
  snp_res <- lapply(snp_res,function(x){
    return(x[[1]])
  })
  sig_vals_all <- list()
  for (i in 1:length(snp_res)) {
    sig_vals <- snp_res[[i]][,'Pr(>|t|)']
    sig_vals_all[[i]] <- sig_vals
  }
  df <- do.call("rbind.data.frame", sig_vals_all)
  return(df)
}

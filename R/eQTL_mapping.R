
#' Fit a mixed model for one SNP and extract coefficients
#'
#' @param data_in data.frame A data object with cells as rows. Expression values
#' are in the first column with PCs in the next columns, followed by donors and
#' other covariate data.
#' @param geno_mat data.frame Donor genotype data (SNPs x donors)
#' @param snp_name character The name of the SNP to test. Should match the name
#' of a row from geno_mat.
#' @param n_PCs numeric The number of PCs and PC-interaction terms to use
#' @param julia environment The output of julia_setup()
#'
#' @return a list with the coefficients matrix in the first element and the
#' variance-covariance matrix in the second element
#' @export
lme_helper <- function(data_in,geno_mat,snp_name,n_PCs,julia) {

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

  return(list(coef_mat,vcov_mat))
}

#' Compute total effects, errors, and p-values per cell for one SNP
#'
#' @param data_in data.frame A data object with cells as rows. Expression values
#' are in the first column with PCs in the next columns, followed by donors and
#' other covariate data.
#' @param coef_mat data.frame The first list element output from lme_helper(), containing
#' the mixed model term coefficients, standard errors, and p-values.
#' @param n_PCs numeric The number of PCs and PC interaction terms used.
#' @param vcov_mat matrix The second list element output from lme_helper(), containing
#' the variance and covariances between each model term.
#'
#' @return a numeric vector containing the p-values for all cells
#' @export
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
#' @param julia_dir character The directory of your Julia installation
#' @param geno_pcs matrix Optional donor by PC matrix to include as covariates in the model.
#' Columns should be labeled as geno_PC1, geno_PC2, etc. (default=NULL)
#' @param n.cores numeric Number of cores to use (default=4)
#' @param progress logical Whether to show a progress bar (default=TRUE)
#'
#' @return A list with one element per cell. Each element contains a vector with the p-values
#' for each tested SNP for the given cell.
#' @export
get_eQTL_res <- function(gene_test,norm_counts,cell_meta,cell_pcs,n_PCs,geno_mat,main_tr,julia_dir,
                         geno_pcs=NULL,n.cores=4,progress=TRUE) {

  if (n_PCs>ncol(cell_pcs)) {
    stop('n_PCs must be >=ncol(cell_pcs)')
  }

  cell_pcs <- cell_pcs[,1:n_PCs]
  
  # ensure geno_mat and main_tr have same snps
  snps_int <- intersect(rownames(geno_mat),rownames(main_tr))
  if (length(snps_int)==0) {
    stop('Need to have the same SNPs in rownames(geno_mat) and rownames(main_tr)')
  } else {
    geno_mat <- geno_mat[snps_int,]
    main_tr <- main_tr[snps_int,]
  }

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

  if (n.cores > 1) {
    # using parLapply for parallel processes
    cl <- makeCluster(n.cores)
    clusterExport(cl, c("julia","julia_dir","lme_helper","get_per_cell_pv_covar","main_tr","data","geno_mat","n_PCs"),
                  envir=environment())
    # clusterEvalQ(cl, c(library(JuliaCall),julia$library("MixedModels")))
    clusterEvalQ(cl,{
      library(JuliaCall)
      julia <- julia_setup(JULIA_HOME = julia_dir,
                           rebuild = TRUE, force = TRUE)
      julia$library("MixedModels")
      julia$library("RCall")
      NULL
    })

    snp_res <- parLapply(cl,X=1:nrow(main_tr),fun = function(snp_j) {
      snp_name <- rownames(main_tr)[snp_j]
      lme_out <- lme_helper(data,geno_mat,snp_name,n_PCs,julia)
      coef_mat <- lme_out[[1]]
      vcov_mat <- lme_out[[2]]
      per_cell_pvals <- get_per_cell_pv_covar(data,coef_mat,n_PCs,vcov_mat)
      return(per_cell_pvals)
    })
    stopCluster(cl)
  } else {
    julia <- julia_setup(JULIA_HOME = julia_dir,
                         rebuild = TRUE, force = TRUE)
    julia$library("MixedModels")
    julia$library("RCall")

    # slightly faster to use lapply if only using one core per gene test
    snp_res <- lapply(1:nrow(main_tr),FUN = function(snp_j) {
      snp_name <- rownames(main_tr)[snp_j]
      lme_out <- lme_helper(data,geno_mat,snp_name,n_PCs,julia)
      coef_mat <- lme_out[[1]]
      vcov_mat <- lme_out[[2]]
      per_cell_pvals <- get_per_cell_pv_covar(data,coef_mat,n_PCs,vcov_mat)
      return(per_cell_pvals)
    })
  }

  snp_res_t <- data.table::transpose(snp_res)
  names(snp_res_t) <- names(snp_res[[1]])
  # snp names match rownames(main_tr)

  return(snp_res_t)
}


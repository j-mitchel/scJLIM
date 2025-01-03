



# For each PC, check if the cells have stronger coloc toward PC center compared to 
# the extremes (test both sides). The side where such pattern exists also must be
# more significant compared to the other side to flag it
detect_pattern <- function(jlim_res,cell_pcs,sig_thresh_pv=.001) {
  ## first looks to see if your per cell jlim pvalues are concentrated in the middle of any individual PC
  # loop through PCs
  sig_thresh <- -log10(sig_thresh_pv)
  sig_vals <- jlim_res[[2]]
  pcs_flag <- c()
  pc_side_flagged <- c()
  for (pc in 1:ncol(cell_pcs)) {
    mypcs <- cell_pcs[,pc]
    names(mypcs) <- rownames(cell_pcs)
    test_bounds <- quantile(abs(mypcs),c(.1,.9))
    all_sides <- c('left','right')
    for (side in all_sides) {
      if (side=='left') {
        cells_inner <- c(names(mypcs)[mypcs<0 & mypcs>(-1*test_bounds[1])])
        cells_middle <- c(names(mypcs)[mypcs<(-1*test_bounds[1]) & mypcs>(-1*test_bounds[2])])
        cells_outer <- c(names(mypcs)[mypcs<(-1*test_bounds[2])])
        cells_side <- c(names(mypcs)[mypcs<0])
        cells_other_side <- c(names(mypcs)[mypcs>0])
      } else if (side=='right') {
        cells_inner <- c(names(mypcs)[mypcs>0 & mypcs<test_bounds[1]])
        cells_middle <- c(names(mypcs)[mypcs>test_bounds[1] & mypcs<test_bounds[2]])
        cells_outer <- c(names(mypcs)[mypcs>test_bounds[2]])
        cells_side <- c(names(mypcs)[mypcs>0])
        cells_other_side <- c(names(mypcs)[mypcs<0])
      }
      pv_inner <- -log10(sig_vals[cells_inner])
      pv_middle <- -log10(sig_vals[cells_middle])
      pv_outer <- -log10(sig_vals[cells_outer])
      pv_side <- -log10(sig_vals[cells_side])
      pv_other_side <- -log10(sig_vals[cells_other_side])
      
      max_inner <- max(pv_inner)
      max_middle <- max(pv_middle)
      max_outer <- max(pv_outer)
      mean_side <- mean(pv_side)
      mean_otherside <- mean(pv_other_side)
      
      if (max_middle>max_inner & max_middle>max_outer & max_middle>sig_thresh & mean_side>mean_otherside) {
        pcs_flag <- c(pcs_flag,pc)
        pc_side_flagged <- c(pc_side_flagged,side)
      }
    }
  }
  return(list(pcs_flag,pc_side_flagged))
}

# helper for below fn
snp_id_helper <- function(jlim_vars, null_dist, sec_tr, sectr.sample.size,
                          min.SNPs.count) {
  # unpack variables
  main_tr <- jlim_vars[[1]]
  assoc1 <- jlim_vars[[2]]
  maf_vec <- jlim_vars[[3]]
  refgt.org <- jlim_vars[[4]]
  refgt_num <- jlim_vars[[5]]
  ld_cormat <- jlim_vars[[6]]
  indexSNP <- jlim_vars[[7]]
  
  NULLDIST <- null_dist[[1]]
  permmat <- null_dist[[2]]
  r2res <- null_dist[[3]]
  
  results.allgene <- matrix(ncol=13, nrow=0)
  colnames(results.allgene) <- c("userIdxBP"," actualIdxBP","STAT", "pvalue",
                                 "usedSNPsNo", "startBP","endBP", "sectrSampleSize",
                                 "sectrGeneName","sectrIdxSNPAssocPvalue", "sectrMinAssocPvalue",
                                 "executedPerm" ,"desc")
  
  # copy the matrix, replace rownames, append Z-scores from the p-values
  assoc2 <- sec_tr
  rownames(assoc2) <- 1:nrow(assoc2)
  if (!('Z' %in% colnames(assoc2))) {
    assoc2 <- cbind(assoc2, Z=PtoZ(assoc2$P))
  }
  
  if (!('P' %in% colnames(assoc2))) {
    assoc2 <- cbind(assoc2, P=pnorm(abs(assoc2$Z), mean = 0, sd = 1, lower.tail = TRUE))
  }
  
  if (!all(assoc1$BP==assoc2$BP)) {
    stop('SNPs in secondary trait are not the same or in the same order as the first trait')
  }
  
  refgt <-  refgt.org
  ld1 <-  ld_cormat
  ld2 <-  ld_cormat
  ld0.maf <-  maf_vec
  
  if(nrow(assoc1) < min.SNPs.count){
    if (nrow(assoc1)< min.SNPs.count)
      stop("too few common SNPs to run JLIM")
  }
  
  refgt0 <- refgt_num
  
  assoc1.sel <- assoc1$BP[assoc1$P <= 0.1]
  assoc2.sel <- assoc2$BP[assoc2$P <= 0.1]
  
  markers.t <- union(assoc1.sel, assoc2.sel)
  
  ASSERT(length(markers.t) > 1)
  ASSERT(assoc1$BP[which.min(assoc1$P)] %in% markers.t)
  
  ld1.t <- ld1[colnames(ld1) %in% markers.t, colnames(ld1)  %in% markers.t]
  ld2.t <- ld2[colnames(ld2) %in% markers.t, colnames(ld2)  %in% markers.t]
  
  assoc1.t <- assoc1[assoc1$BP %in% markers.t, ]
  assoc2.t <- assoc2[assoc2$BP %in% markers.t, ]
  
  if (!is.null(ld0.maf)) {
    ld0.maf.t <- ld0.maf[assoc1$BP %in% markers.t]
  }
  
  best1 <- which.max(abs(assoc1.t$Z))
  sectrIndSNPpvalue <- assoc2$P[assoc2$BP==indexSNP]
  sectrMinpvalue <- min(assoc2$P)
  sectrSNPWithMinpvalue <- assoc2$BP[ assoc2$P==min(assoc2$P)][1]
  jlim.res <- new("jlim",
                  userIdxBP=assoc1.t$BP[best1],
                  actualIdxBP=assoc1.t$BP[best1],
                  STAT=NA_real_, pvalue=NA_real_,
                  usedSNPsNo=nrow(assoc1.t),
                  startBP= min(assoc1.t$BP),
                  endBP= max(assoc1.t$BP),
                  sectrSampleSize=sectr.sample.size,
                  sectrGeneName="",
                  sectrIndSNPpvalue=sectrIndSNPpvalue,
                  sectrMinpvalue=sectrMinpvalue,
                  sectrSNPWithMinpvalue=sectrSNPWithMinpvalue,
                  desc="", executedPerm=0,
                  permmat=permmat)
  
  # check the number of remaining snps in the assoc1
  if(nrow(assoc1.t) < min.SNPs.count ){
    stop("too few SNPs to run JLIM")
  }
  
  assoc1 <- assoc1.t
  assoc2 <- assoc2.t
  ld0 <- ld1.t
  ld2 <- ld2.t
  R2thr <- r2res
  
  logP1 <- (abs(assoc1$Z) ** 2) / 2
  logP2 <- (abs(assoc2$Z) ** 2) / 2
  
  ASSERT(sum(is.na(assoc1$Z)) == 0)
  ASSERT(sum(is.na(assoc2$Z)) == 0)
  
  ### PEAK SELECTION (local)
  maxI1 <- which.max(logP1)
  relP1 <- exp(logP1 - max(logP1))
  postP1 <- exp(logP1) / sum(exp(logP1))
  local <- which(ld0[maxI1, ]**2 >= R2thr)
  
  gap <- 0
  
  gap_per_snp <- c()
  for (I in local) {
    mystat <- relP1[I] * (logP2[I] - max(logP2[ld2[I, ]**2 < R2thr]))
    gap_per_snp <- c(gap_per_snp,mystat)
  }
  snp_ndx <- which(gap_per_snp==max(gap_per_snp))[1]
  bp_best <- names(local)[snp_ndx]
  return(bp_best)
}

# identify the top snp that is most contributing to the colocalization
identify_snps_contrib <- function(jlim_res,pc_flags,cell_pcs,snp_res_mat,jlim_vars,null_dist,
                                  sectr.sample.size,min.SNPs.count,sig_thresh_pv=.001) {
  main_tr <- jlim_vars[[1]]
  sig_vals <- jlim_res[[2]]
  pcs_flag <- pc_flags[[1]]
  pc_side_flagged <- pc_flags[[2]]
  pc_snps <- c()
  for (pc_ndx in 1:length(pcs_flag)) {
    pc <- pcs_flag[pc_ndx]
    side <- pc_side_flagged[pc_ndx]
    mypcs <- cell_pcs[,pc]
    names(mypcs) <- rownames(cell_pcs)
    if (side=='left') {
      cells_sub <- names(mypcs)[mypcs<0]
    } else if (side=='right') {
      cells_sub <- names(mypcs)[mypcs>0]
    }
    top_cells_sig <- names(sig_vals)[sig_vals<sig_thresh_pv]
    top_cells_select <- intersect(cells_sub,top_cells_sig)
    all_best_bp <- c()
    for (mycell in top_cells_select) {
      ## now, for these cells, need to regenerate the jlim stats to see which snps contribute to significance
      snp_res_un <- snp_res_mat[[mycell]]
      # select pvalues for all snps for a cell
      names(snp_res_un) <- rownames(main_tr)
      ## make sec_tr df
      sec_tr <- cbind.data.frame(main_tr$CHR,main_tr$BP,snp_res_un)
      colnames(sec_tr) <- c('CHR','BP','P')
      
      best_bp <- snp_id_helper(jlim_vars, null_dist, sec_tr, sectr.sample.size, min.SNPs.count)
      all_best_bp <- c(all_best_bp,best_bp)
    }
    best_bp_counts <- table(all_best_bp)
    best_bp <- names(best_bp_counts)[order(best_bp_counts,decreasing=TRUE)][1]
    pc_snps <- c(pc_snps,best_bp)
  }
  return(pc_snps)
}

# check if the magnitude of the effect (or significance) of the tested snp gets smaller as the pc gets bigger
check_snp_direc <- function(snp_res_mat,jlim_vars,cell_pcs,pc_flags,pc_snps) {
  main_tr <- jlim_vars[[1]]
  pcs_flag <- pc_flags[[1]]
  pc_side_flagged <- pc_flags[[2]]
  final_flag <- c()
  for (pc_ndx in 1:length(pc_snps)) {
    bp <- pc_snps[pc_ndx]
    snp_ndx <- which(main_tr$BP==bp)
    snp_pvals <- sapply(snp_res_mat,function(x){
      return(x[snp_ndx])
    })
    snp_pvals <- -log10(snp_pvals)
    tmp <- cbind.data.frame(cell_pcs,snp_pvals)
    pc_nm <- paste0('PC',pcs_flag[pc_ndx])
    myform <- as.formula(paste0('snp_pvals ~ ',pc_nm))
    lmres <- summary(lm(myform,data=tmp))
    lmeffect <- lmres$coefficients[pc_nm,'Estimate']
    side <- pc_side_flagged[pc_ndx]
    
    # also generally checking significance levels of outer bounds on either side
    outer_bounds <- quantile(tmp[,pc_nm],probs = c(.1,.9))
    upper_bound_pv <- max(tmp[tmp[,pc_nm]>outer_bounds[2],'snp_pvals'])
    lower_bound_pv <- max(tmp[tmp[,pc_nm]<outer_bounds[1],'snp_pvals'])
    
    inner_bounds <- quantile(tmp[,pc_nm],probs = c(.4,.6))
    middle_pv <- max(tmp[(tmp[,pc_nm]>inner_bounds[1] & tmp[,pc_nm]<inner_bounds[2]),'snp_pvals'])
    
    # checking snp is more significant opposite from side of jlim significnace
    # and that the SNPs significance trend is monotonic
    if (side=='left') {
      if (lmeffect>0 & upper_bound_pv>middle_pv & middle_pv>lower_bound_pv) {
        flag_log <- TRUE
      } else {
        flag_log <- FALSE
      }
    } else if (side=='right') {
      if (lmeffect<0 & lower_bound_pv>middle_pv & middle_pv>upper_bound_pv) {
        flag_log <- TRUE
      } else {
        flag_log <- FALSE
      }
    }
    final_flag <- c(final_flag,flag_log)
  }
  return(final_flag)
}
















































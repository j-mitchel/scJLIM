

#' jlim class to keep the result of JLIM single test
#'
#'#' @useDynLib scJLIM, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' 
#' @slot userIdxBP user specified index SNP
#' @slot actualIdxBP actual found index SNP
#' @slot STAT jlim statistic lambda
#' @slot pvalue permutation pvalue
#' @slot startBP start position of the tested locus
#' @slot endBP end position of the tested locus
#' @slot sectrSampleSize end position of the tested locus
#' @slot sectrGeneName name of the Gene in case of multiple gene in the second trait
#' @slot sectrIndSNPpvalue pvalue of the indexSNP in the second trait.
#' @slot sectrMinpvalue minimum pvalue of in the second trait.
#' @slot sectrSNPWithMinpvalue  SNP with the minimum pvalue in the second trait
#' @slot desc  status of the JLIM test
#' @slot executedPerm number of the executed permutations
#' @slot permmat the permutation matrix used for the JLIM test
#' @export
#' @import methods

setClass("jlim",
         slots = list(
           userIdxBP="numeric",
           actualIdxBP="numeric",
           STAT="numeric",
           pvalue="numeric",
           usedSNPsNo="numeric",
           startBP="numeric",
           endBP="numeric",
           sectrSampleSize="numeric",
           sectrGeneName="character",
           sectrIndSNPpvalue="numeric",
           sectrMinpvalue="numeric",
           sectrSNPWithMinpvalue="numeric",
           desc="character",
           executedPerm="numeric",
           permmat="matrix"
         ))

#getVec.jlim <- function(object){}
setGeneric("getVec.jlim", function(object) standardGeneric("getVec.jlim"))
setMethod("getVec.jlim",
          "jlim",
          function(object) {
            c(object@userIdxBP, object@actualIdxBP, object@STAT, object@pvalue,
              object@usedSNPsNo, object@startBP, object@endBP, object@sectrSampleSize,
              object@sectrGeneName, object@sectrIndSNPpvalue, object@sectrMinpvalue,
              object@executedPerm, object@desc)
          })



# loads the reference panel and filters out low MAF variants
loadRefLD <- function(refgt, start.bp, end.bp, assoc2.genes.org,
                       assoc2.genes, min.MAF ){
  ld0.maf <- NULL
  
  #  reference gt from 1000 genomes
  colnames(refgt)[1:5] <- c("CHROM", "POS", "REF", "ALT", "GT")
  
  # remove individuals with missing data
  refgt <- refgt[complete.cases(refgt[ , 6:ncol(refgt)]),]
  
  # subset reference data to the snps we have in the data
  refgt <- refgt[refgt$POS >= start.bp & refgt$POS <= end.bp, ]
  
  # filter rare variant from refLD
  ld0.maf <- calcMAF.refPanels(refgt, rep(TRUE, nrow(refgt)))
  
  refgt <- refgt[ld0.maf >= min.MAF, ]
  ld0.maf <- ld0.maf[ld0.maf >= min.MAF]
  
  # exclude overlapping variants
  if (sum(duplicated(refgt$POS)) > 0) {
    dup.POS <- refgt$POS[duplicated(refgt$POS)]
    
    MSG("\nRemoving multiple common alleles at same BP in refgt: BP =",
        paste(dup.POS, collapse=", "))
    
    ld0.maf <- ld0.maf[!(refgt$POS %in% dup.POS)]
    refgt <- refgt[!(refgt$POS %in% dup.POS), ]
  }
  
  reslist <- list(refgt, ld0.maf)
  return (reslist)
}

# convert matrix of genotypes into matrix of 1s and 2s for ref/alt
# also getting the LD matrix of pearson correlations
process.refPanels <- function(refgt.org) {
  
  refgt.org <- refgt.org[order(refgt.org$POS), ]
  
  refgt.mat <- refgt.org[, 6:ncol(refgt.org)]
  refgt.mat <- as.matrix(refgt.mat)
  
  ASSERT(sum(refgt.mat == ".") == 0)
  # ref allele := 1
  refgt.mat <-
    t(sapply(1:nrow(refgt.mat),
             function(I) as.numeric(refgt.mat[I, ] == refgt.org$REF[I]),
             simplify=TRUE))
  # alt allele := 2
  refgt.mat[refgt.mat != 1] <- 2
  refgt.mat0 <- toGT(refgt.mat)
  rownames(refgt.mat0) <- refgt.org$POS
  rownames(refgt.mat) <- refgt.org$POS

  mode(refgt.mat0) = "numeric"
  gt0 <- stdGT2(t(refgt.mat0))
  # product of standardized genotype vectors to get the ld matrix, essentially a dot product for each pair of snps
  # if a pair is highly correlated across individuals, it will have a large dot product
  # it's just the covariance matrix
  ld0 <- (t(gt0) %*% gt0) # should be dimensions of snps x snps
  # next, just normalizing the ld matrix so the values are pearson correlations
  var0 <- apply(gt0, 2, function(v) sqrt(sum(v*v))) # getting a vector of norms for each variant
  for (I in 1:nrow(ld0)) {
    for (J in 1:nrow(ld0)) {
      ld0[I, J] <- ld0[I, J]/var0[I]/var0[J]
    }
  }
  
  return(list(ld0,refgt.mat))
}

calc.stat <- function (assoc1, assoc2, ld0, ld2, R2thr) {
  
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
  sumRelP1 <- sum(relP1[local])
  
  for (I in local) {
    gap <- gap +
      relP1[I] * (logP2[I] - max(logP2[ld2[I, ]**2 < R2thr]))
  }
  
  gap <- gap / sumRelP1
  gap
  
  return(gap)
}

#' @export
perm.test <- function (assoc1, permmat, ld0, ld2,
                       r2res, lambda.t, n.cores, progress=progress) {
  
  thresholdingZ <- PtoZ(0.1)
  
  logP1 <- (abs(assoc1$Z) ** 2) / 2 # this is the non-centrality parameter
  
  # SNP selection part 1
  markers.p1 <- abs(assoc1$Z) >= thresholdingZ
  
  # peak selection
  maxlogP1 <- max(logP1) # max non-centrality param
  maxI1 <- which.max(logP1) # index of the max
  relP1 <- exp(logP1 - maxlogP1) # e^difference for each snp. 
  # ~0 if much less than max, higher if closer
  
  simNo <- 1
  NULLGAP <- c()
  ASSERT(ncol(permmat) == nrow(ld0))
  ASSERT(ncol(permmat) == nrow(ld2))
  
  cl <- makeCluster(n.cores)
  clusterExport(cl, c("ASSERT","abs","max","sum","is.na","which","permmat","markers.p1","thresholdingZ","maxI1","ld0","r2res","relP1"), 
                envir=environment())
  NULLGAP <- parLapply(cl, X = 1:nrow(permmat), fun = function(simNo) {
    assoc2n.Z <- permmat[simNo, ]
    
    ASSERT(sum(is.na(assoc2n.Z)) == 0)
    
    # SNP selection part 2
    markers.p <-
      markers.p1 | (abs(assoc2n.Z) >= thresholdingZ)
    
    ASSERT(markers.p[maxI1]) # always include maxI1
    
    # need multiple markers
    local <- intersect(which(ld0[maxI1, ]**2 >= r2res),
                       which(markers.p))
    
    logP2n <- (abs(assoc2n.Z) ** 2) / 2
    
    # gap is the joint likelihood statistic
    gap <-
      sum(relP1[local] *
            sapply(local, function (I)
              (logP2n[I] -
                 max(logP2n[(ld2[I, ]**2 < r2res) & markers.p])),
              simplify=TRUE))
    
    gap_norm <- gap / sum(relP1[local]) # statistic gets normalized
    return(gap_norm)
  })
  stopCluster(cl)
  
  NULLGAP <- unlist(NULLGAP)
  
  return(NULLGAP)
}

# refgt_num is the reference panel with ref/alt as 1s and 2s
get_permmat <- function(refgt_num, sectr.sample.size, nperm, n.cores, progress) {
  snpIds <- rownames(refgt_num)
  refgt_num <- toGT(refgt_num)
  refgt_num <- stdGT2(t(refgt_num))
  refgt_num <- t(refgt_num)
  rownames(refgt_num) <- snpIds
  
  y <- matrix(rnorm(sectr.sample.size * nperm), ncol=sectr.sample.size, nrow=nperm) # represent gene expression values sampled from the null (random normal distribution)
  # permmat <- matrix(0, nrow=nperm, ncol=nrow(refgt_num))
  refgt_num <- t(refgt_num)
  
  cl <- makeCluster(n.cores)
  clusterExport(cl, c("refgt_num","sample","nrow","sectr.sample.size","y"), 
                envir=environment())
  permmat_list <- parLapply(cl, X=1:nperm, fun = function(IP) {
    sampledgt <- refgt_num[sample(1:nrow(refgt_num), sectr.sample.size, replace=TRUE),]
    
    z_vec <- c(y[IP, , drop=FALSE] %*% sampledgt)
    
    return(z_vec)
  })
  stopCluster(cl)
  
  permmat <- do.call(rbind,permmat_list)
  permmat <- permmat / sqrt(sectr.sample.size)
  
  # permmat has dimensions of nperm by snps, values are dot product of expression x genotype_snp_i
  # thus, if expression is correlated with genotype, we get a large value and vice vera
  # however, since the expression values sampled from the null, we should get a null distribution of covariance values
  # denominator is standard deviation of sum of product of two rnorms
  # for example:
  # The variance of the product, var_prod, of the two standardized variables is 1
  # the sum of this product will have variance var_sum = n_samples*var_prod
  # thus, the sd of the sum of the product is sd_sum_prod = sqrt(n_samples*var_prod)
  # so, sd_sum_prod = sqrt(n_samples)
  # therefore, dividing the values in permmat by this, converts the values into z-scores
  return(permmat)
}

# converts p-values to Z-scores
PtoZ <- function(pv) {
  -qnorm(pv/2)
}

# creates genotypes from haplotypes by summing number of alternate alleles from two random haplotypes
toGT <- function(gt) {
  M <- ncol(gt)/2
  gt.alt <- (gt != 1) + 0 # nonref
  matrix(apply(gt.alt, 1, function(v) (v[seq(1, 2*M, by=2)] + v[seq(2, 2*M, by=2)])),
         byrow=TRUE, ncol=M)
}

# standardize 0s and 1s matrix of genotypes
stdGT2 <- function(gt) {
  f <- apply(gt, 2, mean)/2
  for (I in 1:ncol(gt)) {
    gt[, I] <- (gt[, I] - 2*f[I]) # for each donor, subtract of their genotype mean
  }
  
  sd0 <- apply(gt, 2, sd)
  for (I in 1:ncol(gt)) {
    gt[, I] <- gt[, I] / sd0[I] # gets a z-score per snp per donor
  }
  
  return(gt)
}

ASSERT <- function(test) {
  if (!test) {
    stop(paste("ASSERT fail:", deparse(substitute(test))),"\n")
  }
}

# find the proper ref-panel from the list of the files in the specified path
find.panel <-function( refLD.dir, start.bp, end.bp, CHR ) {
  
  ref.LD.list <- vector()
  ref.LD.list0 <-
    list.files(refLD.dir,
               paste("[^.]+.*.txt.gz", sep=''))

  ref.LD.list <- gsub("_", ".", ref.LD.list0)
  
  panels.List<-
    sapply(ref.LD.list,
           function(str) strsplit(str, '.', fixed=TRUE),
           simplify=TRUE)
  
  panels <- as.data.frame(matrix(unlist(panels.List), ncol  = 5, byrow = TRUE))
  panels[,2] <- as.numeric(as.character(panels[,2]))
  panels[,3] <- as.numeric(as.character(panels[,3]))
  
  panel.se <- which(panels[,1]==paste("chr",CHR,sep="") & panels[,2] <= start.bp & panels[,3]>=end.bp)
  
  if(length(panel.se)==0){
    cat("\nThere is not any matching panel for the given interval: ", start.bp," - ",end.bp,"\n" )
    q(status=1)
  }else if(length(panel.se)==1){
    panel.toUse <- panel.se
  }else{
    panel.toUse <- panel.se[1]
  }
  
  return(paste(refLD.dir,"/",ref.LD.list0[panel.toUse],sep=""))
}

calcMAF.refPanels <- function(refgt, refgt.all.sel) {
  
  refgt0 <- refgt[refgt.all.sel, ]
  refgt0 <- refgt0[order(refgt0$POS), ]
  
  ASSERT(sum(refgt.all.sel) == nrow(refgt0))
  
  refgt.mat <- refgt0[, 6:ncol(refgt0)]
  refgt.mat <- as.matrix(refgt.mat)
  
  # ref allele := 1
  refgt.mat0 <-
    t(sapply(1:nrow(refgt.mat),
             function(I) as.numeric(refgt.mat[I, ] == refgt0$REF[I]),
             simplify=TRUE))
  # alt allele := 2
  refgt.mat0[refgt.mat0 != 1] <- 2
  rownames(refgt.mat0) <- refgt0$POS
  
  refgt.mat <- toGT(refgt.mat0)
  rownames(refgt.mat) <- refgt0$POS
  
  meanF0 <- apply(refgt.mat, 1, mean)/2
  meanF0 <- pmin(meanF0, 1 - meanF0)
  
  return(meanF0)
}

# function to generate a single null distribution and other needed inputs for main
# jlim.test function
# this only gets run once per full gene-cell type test (not per cell)
# should be run before colocalization
# should not depend on sec_tr because this changes for each cell in coloc test
prep_jlim <- function(main_tr,refLD.dir=NULL,refLD.mat=NULL,min.MAF=.05) {
  if (is.null(refLD.dir) & is.null(refLD.mat)) {
    stop('need to supply either refLD.dir or refLD.mat')
  }
  verbose <- FALSE
  start.bp <- min(main_tr$BP)
  end.bp <- max(main_tr$BP)
  CHR <- main_tr$CHR[1]
  
  if (!is.null(refLD.dir)) {
    refld.file <- find.panel(refLD.dir, start.bp=start.bp,
                             end.bp=end.bp, CHR=CHR)
    refgt <- read.delim(file=refld.file, header=FALSE, sep="\t", stringsAsFactors=FALSE) # matrix of variants x individuals from reference
  } else {
    refgt <- refLD.mat
    refgt[,1] <- as.integer(refgt[,1])
    refgt[,2] <- as.integer(refgt[,2])
  }

  # process the refgt matrix, remove low MAF variants, and get out maf vector
  refLD.res <- loadRefLD(refgt, start.bp, end.bp, NULL,
                         NULL, min.MAF )
  refgt.org <- refLD.res[[1]] # maf cleaned refgt
  maf_vec <- refLD.res[[2]] # vector of mafs - names are BP. previously named ld0.maf.org
  
  # get SNPs in both ref matrix and main_tr
  bp_both <- intersect(as.numeric(refgt.org[,'POS']),main_tr$BP)
  
  # subset ref matrix, ld vector, and main_tr to these snps
  refgt.org <- refgt.org[refgt.org$POS %in% bp_both,]
  maf_vec <- maf_vec[as.character(bp_both)]
  main_tr <- main_tr[main_tr$BP %in% bp_both,]
  
  # sort variants by pos
  main_tr <- main_tr[match(sort(main_tr$BP),main_tr$BP),]

  # get the index snp (lead GWAS SNP)
  indexSNP <- main_tr[main_tr$P==min(main_tr$P),'BP'][1]
  
  # copy the matrix, replace rownames, append Z-scores from the p-values
  assoc1 <- main_tr
  rownames(assoc1) <- 1:nrow(assoc1)
  assoc1 <- cbind(assoc1, Z=PtoZ(assoc1$P))
  assoc1 <- assoc1[!is.na(assoc1$P), ]
  
  ld_cormat_refgt_num <- process.refPanels(refgt.org)
  ld_cormat <- ld_cormat_refgt_num[[1]] # previously called ld0.org
  refgt_num <- ld_cormat_refgt_num[[2]] # previously called refgt0.org

  return(list(main_tr,assoc1,maf_vec,refgt.org,refgt_num,ld_cormat,indexSNP))
}

get_null_dist <- function(jlim_vars,sectr.sample.size,nperm,n.cores,r2res=.8,progress=TRUE) {
  # unpack variables
  main_tr <- jlim_vars[[1]]
  assoc1 <- jlim_vars[[2]]
  maf_vec <- jlim_vars[[3]]
  refgt.org <- jlim_vars[[4]]
  refgt_num <- jlim_vars[[5]]
  ld_cormat <- jlim_vars[[6]]
  indexSNP <- jlim_vars[[7]]
  
  ## get permutation matrix
  permmat <- get_permmat(refgt_num, sectr.sample.size, nperm, n.cores=n.cores, progress=progress)
  
  # subset permutation matrix to same variants as in assoc1
  NULLDIST <- perm.test(assoc1, permmat=permmat,ld0=ld_cormat, ld2=ld_cormat,
                        r2res=r2res, lambda.t=Inf, n.cores=n.cores, progress=progress)
  
  return(list(NULLDIST,r2res))
}

jlim.test <- function(jlim_vars, null_dist, sec_tr, sectr.sample.size, min.SNPs.count){
  
  # unpack variables
  main_tr <- jlim_vars[[1]]
  assoc1 <- jlim_vars[[2]]
  maf_vec <- jlim_vars[[3]]
  refgt.org <- jlim_vars[[4]]
  refgt_num <- jlim_vars[[5]]
  ld_cormat <- jlim_vars[[6]]
  indexSNP <- jlim_vars[[7]]
  
  NULLDIST <- null_dist[[1]]
  r2res <- null_dist[[2]]
  
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
  
  jlim.res <- SNPselction(assoc1, assoc2, ld1=ld1, ld2=ld2, ld0.maf, r2res = r2res,
                             refgt_num, sectr.sample.size, min.SNPs.count,
                          indexSNP=indexSNP, NULLDIST)
  jlim.res@userIdxBP <- indexSNP

  # cat("\nJLIM results:",colnames(results.allgene),"\n",sep = "   ")
  # cat("\n",getVec.jlim(jlim.res))
  results.allgene <- rbind (results.allgene, getVec.jlim(jlim.res))
  
  return(results.allgene)
}


## adapted from JLIM package
SNPselction <- function(assoc1, assoc2, ld1, ld2, ld0.maf, r2res,
                        refgt0, sectr.sample.size,
                        min.SNPs.count, indexSNP, NULLDIST){
  
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
                  desc="", executedPerm=0)
  
  # check the number of remaining snps in the assoc1
  if(nrow(assoc1.t) < min.SNPs.count ){
    stop("too few SNPs to run JLIM")
  }
  
  lambda.t <- calc.stat(assoc1.t, assoc2.t, ld1.t, ld2.t, r2res)
  
  #########################################################################
  executedPerm <- length(NULLDIST)
  
  permP <- sum(NULLDIST >= lambda.t, na.rm=TRUE)/sum(!is.na(NULLDIST))
  
  jlim.res@desc <- "executed"
  jlim.res@STAT=lambda.t
  jlim.res@pvalue=permP
  jlim.res@executedPerm=executedPerm
  return(jlim.res)
}

jlim_main <- function(snp_res_mat, jlim_vars, null_dist, sectr.sample.size, 
                      min.SNPs.count=15, n.cores=20, global_adjust_method='cauchy',
                      adjust_method='n_eff_cauchy',cluster_pvals=FALSE,n_eff=NULL,
                      progress=TRUE,censor_cells=FALSE,censor_type='min_pval') {
  main_tr <- jlim_vars[[1]]
  
  cl <- makeCluster(n.cores)
  clusterExport(cl, c("calc.stat","ASSERT","jlim.test","PtoZ","SNPselction",
                      "as.numeric","cbind.data.frame","ACAT","main_tr",
                      "jlim_vars","null_dist","sectr.sample.size",
                      "min.SNPs.count","censor_cells","censor_type"), 
                envir=environment())
  clusterEvalQ(cl, c(library(ACAT)))
  
  per_cell_jlim <- parLapply(cl,snp_res_mat,function(snp_res_un) {
    setClass("jlim",
             slots = list(
               userIdxBP="numeric",
               actualIdxBP="numeric",
               STAT="numeric",
               pvalue="numeric",
               usedSNPsNo="numeric",
               startBP="numeric",
               endBP="numeric",
               sectrSampleSize="numeric",
               sectrGeneName="character",
               sectrIndSNPpvalue="numeric",
               sectrMinpvalue="numeric",
               sectrSNPWithMinpvalue="numeric",
               desc="character",
               executedPerm="numeric",
               permmat="matrix"
             ))
    
    #getVec.jlim <- function(object){}
    setGeneric("getVec.jlim", function(object) standardGeneric("getVec.jlim"))
    setMethod("getVec.jlim",
              "jlim",
              function(object) {
                c(object@userIdxBP, object@actualIdxBP, object@STAT, object@pvalue,
                  object@usedSNPsNo, object@startBP, object@endBP, object@sectrSampleSize,
                  object@sectrGeneName, object@sectrIndSNPpvalue, object@sectrMinpvalue,
                  object@executedPerm, object@desc)
              })
    
    if (is.na(snp_res_un[1])) {
      return(list(FALSE,1))
    }
    
    # select pvalues for all snps for a cell
    names(snp_res_un) <- rownames(main_tr)
    ## make sec_tr df
    
    # if there are any 0 eQTL pvalues, set to the largest non-zero pvalue
    if (any(snp_res_un==0)) {
      snp_res_un[snp_res_un==0] <- min(snp_res_un[snp_res_un!=0])
    }
    
    sec_tr <- cbind.data.frame(main_tr$CHR,main_tr$BP,snp_res_un)
    colnames(sec_tr) <- c('CHR','BP','P')
    
    if (censor_type=='min_pval') {
      eGene_sig <- min(sec_tr$P)
      eGene_thresh <- .01
    } else if (censor_type=='cauchy_pval') {
      eGene_sig <- ACAT(sec_tr$P)
      eGene_thresh <- .10
    }
    
    if (censor_cells) {
      if (eGene_sig>eGene_thresh) {
        c_test <- FALSE
      } else {
        c_test <- TRUE
      }
    } else {
      c_test <- TRUE
    }
    
    jlim_res <- jlim.test(jlim_vars, null_dist, sec_tr, sectr.sample.size=sectr.sample.size,
                          min.SNPs.count=min.SNPs.count)
    pval <- as.numeric(jlim_res[1,'pvalue'])
    return(list(c_test,pval))
    # return(pval)
  })
  stopCluster(cl)
  
  cells_test_ind <- sapply(per_cell_jlim,function(i){
    return(i[[1]])
  })
  per_cell_jlim_un <- sapply(per_cell_jlim,function(i){
    return(i[[2]])
  })
  
  names(per_cell_jlim_un) <- names(snp_res_mat)
  # per_cell_jlim_un <- per_cell_jlim_un[!is.na(per_cell_jlim_un)]

  # if (length(per_cell_jlim_un)==0) {
  #   return(NA)
  # }
  
  #### to run cauchy combination test, can't have pvals exactly 0 or 1
  # # set pvals==0 to be equal to 1/number of permuatations
  # per_cell_jlim_un[per_cell_jlim_un==0] <- 1/length(null_dist[[1]])
  # # set pvals==1 to be 1-(1/number of permuatations)
  # per_cell_jlim_un[per_cell_jlim_un==1] <- 1-(1/length(null_dist[[1]]))
  
  if (sum(cells_test_ind)==0) {
    return(list(NA,per_cell_jlim_un,c()))
  }
  cells_test_ind <- which(cells_test_ind) # logical to indices
  
  ## removing pvals exactly equal to 0 or 1
  ndx_rm1 <- which(per_cell_jlim_un==0)
  ndx_rm2 <- which(per_cell_jlim_un==1)
  cells_test_ind <- cells_test_ind[!(cells_test_ind %in% c(ndx_rm1,ndx_rm2))]

  if (global_adjust_method %in% c('n_eff_bon','cauchy_omni','all_methods')) {
    if (is.null(n_eff)) {
      stop('need to set n_eff with this adjust method')
    }
  }
  
  if (global_adjust_method=='cauchy') {
    ## compute p-value globally with cauchy combination test
    global_p <- ACAT(per_cell_jlim_un[cells_test_ind])
  } else if (global_adjust_method=='n_eff_bon') {
    global_p <- min(c(min(per_cell_jlim_un) * n_eff,1))
  } else {
    stop('use one of the available global adjust methods')
  }

  if (adjust_method=='n_eff_cauchy') {
    # ## the following logic used to adjust individual cell pvals
    # pval_adj <- pval_orig * effective_n  #for bonferroni correction
    # pval_global <- min(pval_orig * effective_n)
    # # this can be rewritten as the following since it works monotonically
    # pval_global <- min(pval_orig) * effective_n
    # # thus, we can compute the effective n with the global pval
    # take min orig pval
    min_orig_pval <- min(per_cell_jlim_un)
    # now divide global pval by this to get effective_n
    effective_n <- global_p / min_orig_pval
    # now adjust all individual pvalues by this
    pval_adj <- pmin(per_cell_jlim_un * effective_n,1)
  } else if (adjust_method=='n_eff_bon') {
    if (!is.null(n_eff)) {
      pval_adj <- pmin(per_cell_jlim_un * n_eff,1)
    } else {
      stop('need to provide n_eff with n_eff_pre correction method')
    }
  } else if (adjust_method=='fdr') {
    pval_adj <- p.adjust(per_cell_jlim_un,method='fdr')
  } else if (adjust_method=='none') {
    pval_adj <- per_cell_jlim_un
  } else {
    stop('pick an available adjust_method')
  }

  if (cluster_pvals) {
    clust_out <- kmeans(-log10(per_cell_jlim_un), 2, iter.max = 10, nstart = 1)
    return(list(global_p,pval_adj,clust_out,cells_test_ind))
  } else {
    return(list(global_p,pval_adj,cells_test_ind))
  }
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector indicesGreaterThan(NumericVector x, double threshold) {
  int n = x.size();
  IntegerVector indices;
  
  for (int i = 0; i < n; ++i) {
    if (x[i] >= threshold) {
      indices.push_back(i);
    }
  }
  
  return indices;
}

// [[Rcpp::export]]
IntegerVector indicesIsTrue(LogicalVector x) {
  int n = x.size();
  IntegerVector indices;
  bool true_val = 1;
  
  for (int i = 0; i < n; ++i) {
    if (x[i] == true_val) {
      indices.push_back(i);
    }
  }
  
  return indices;
}

//' Perform permutation test
//'
//' This function performs a permutation test using the given parameters.
//'
//' @param assoc1 A numeric vector representing assoc1.
//' @param permmat A numeric matrix representing permmat.
//' @param ld0 A numeric matrix representing ld0.
//' @param ld2 A numeric matrix representing ld2.
//' @param R2thr A numeric value representing R2thr.
//' @param lambda_t A numeric value representing lambda_t.
//' @return A numeric vector containing the results of the permutation test.
//' @export
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericVector perm_test(DataFrame assoc1, NumericMatrix permmat, NumericMatrix ld0,
                        NumericMatrix ld2, double R2thr, double lambda_t) {
  
  // Convert assoc1$Z to NumericVector
  NumericVector assoc1_Z = assoc1["Z"];
  
  // Z threshold - from PtoZ(0.1)
  double thresholdingZ = 1.644854;

  // Calculate logP1
  NumericVector logP1 = pow(abs(assoc1_Z), 2) / 2.0;

  // SNP selection part 1
  LogicalVector markers_p1 = abs(assoc1_Z) >= thresholdingZ;
  
  double maxlogP1 = max(logP1);
  
  int maxI1 = which_max(logP1);
  
  NumericVector relP1_diff = logP1 - maxlogP1;

  NumericVector relP1 = exp(relP1_diff);
  
  NumericVector NULLGAP;

  // Iterate over permutations
  for (int simNo = 0; simNo < permmat.nrow(); simNo++) {
    // Extract assoc2n.Z from permmat
    NumericVector assoc2n_Z = permmat(simNo, _);
    
    // SNP selection part 2
    LogicalVector markers_p = markers_p1 | (abs(assoc2n_Z) >= thresholdingZ);
    
    // NumericVector local = intersect(indicesGreaterThan(pow(ld0(maxI1, _), 2), R2thr), which(markers_p)) - 1;
    
    NumericVector ld_r_to_lead = ld0(maxI1, _);
    NumericVector ld_rsq_to_lead = pow(ld_r_to_lead, 2);
    IntegerVector markers_in_ld_to_lead_ndx = indicesGreaterThan(ld_rsq_to_lead,R2thr);
    IntegerVector markers_p_ndx = indicesIsTrue(markers_p);
    IntegerVector local = intersect(markers_in_ld_to_lead_ndx,markers_p_ndx);
    
    NumericVector logP2n = pow(abs(assoc2n_Z), 2) / 2.0;
    
    // to calculate gap stat first loop through local
    int n = local.size();
    NumericVector t2_diffs;
    for (int i = 0; i < n; ++i) {
      int I = local[i];
      double t2_stat = logP2n[I];
      NumericVector ld_r_to_snp_i = ld0(I, _);
      NumericVector ld_rsq_to_snp_i = pow(ld_r_to_snp_i, 2);
      LogicalVector markers_in_ld_to_snp_i = (ld_rsq_to_snp_i < R2thr) & markers_p;
      NumericVector t2_stats_outside = logP2n[markers_in_ld_to_snp_i];
      double max_outside = max(t2_stats_outside);
      double mydiff = t2_stat - max_outside;
      t2_diffs.push_back(mydiff);
    }
    
    NumericVector relP1_local = relP1[local];
    NumericVector gap_product = relP1_local * t2_diffs;
    double gap = sum(gap_product);
    double gap_norm = gap / sum(relP1_local);
    NULLGAP.push_back(gap_norm);
  }
  
  return NULLGAP;
}



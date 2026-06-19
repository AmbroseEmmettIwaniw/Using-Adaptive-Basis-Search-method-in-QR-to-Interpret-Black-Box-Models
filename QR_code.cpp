#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector d_register(IntegerVector r, int k) {
  int d = r.size();
  
  for (int i = (d-1); i >= 0; i--)
  {
    if (r[i] < k-1)
    {
      r[i]++;
      break;
    }
    else if (r[i] == k-1)
    {
      r[i] = 0;
    }
  }
  return r;
}

// [[Rcpp::export]]
std::vector<IntegerVector> interpretable_sets(int d, int B_0, int B_1, int B_inf)
{
  int k = B_inf + 1;
  std::vector<IntegerVector> U;
  IntegerVector r(d);
  
  for (int i = 0; i < pow(k,d); i++)
  {
    if (sum(r) <= B_1 && sum(r > 0) <= B_0)
    {
      U.push_back(clone(r));
    }
    r = clone(d_register(r,k));
  }
  
  return U;
}

// [[Rcpp::export]]
double orthonormal_basis(double z, int degree)
{
  double u = z - 1.0/2.0;
  double phi_d = 0;
  
  switch(degree)
  {
  case 0:
    phi_d = 1.0;
    break;
  case 1:
    phi_d = sqrt(12.0)*u;
    break;
  case 2:
    phi_d = sqrt(180.0)*(pow(u,2) - 1.0/12.0);
    break;
  case 3:
    phi_d = sqrt(2800.0)*(pow(u,3) - 3.0/20*u);
    break;
  case 4:
    phi_d = 210*(pow(u,4) - 3.0/14.0*pow(u,2) + 3.0/560.0);
    break;
  case 5:
    phi_d = 252*sqrt(11.0)*(pow(u,5) -5.0/18.0*pow(u,3) + 5.0/336.0*u);
    break;
  case 6:
    phi_d = 924*sqrt(13.0)*(pow(u,6) - 15.0/44.0*pow(u,4) + 5.0/176.0*pow(u,2) - 5.0/14784.0);
    break;
  default:
    Rcout << "The degree of 7 or greater is not defined." << std::endl;
  break;
  }
  return phi_d;
}

// [[Rcpp::export]]
double mulidimensional_basis(const IntegerVector& r, const NumericVector& x)
{
  double phi_tensor_d = 1;
  int d = r.size();
  
  for (int i = 0; i < d; i++)
  {
    phi_tensor_d *= orthonormal_basis(x[i],r[i]);
  }
  
  return phi_tensor_d;
}

// [[Rcpp::export]]
std::vector<IntegerVector> interpretable_poly_helper(int d, int B_1, int B_inf)
{
  int max_r = B_inf + 1;
  
  std::vector<IntegerVector> U;
  IntegerVector r(d);
  
  for (int i = 0; i < pow(max_r,d); i++)
  {
    if (sum(r) <= B_1 && sum(r==0) == 0)
    {
      U.push_back(clone(r));
    }
    r = clone(d_register(r,max_r));
  }
  
  return U;
}

// [[Rcpp::export]]
std::vector<IntegerVector> interpretable_poly_sets(int d, int B_0, int B_1, int B_inf)
{
  int max_r = B_inf;
  
  std::vector<IntegerVector> U;
  
  IntegerVector r(d);
  U.push_back(clone(r)); // The empty element.
  
  if (B_0 >= 1) // Single factors.
  {
    for (int i = 0; i < d; i++)
    {
      IntegerVector r1(d);
      for (int j = 0; j < max_r; j++)
      {
        r1[i]++;
        U.push_back(clone(r1));
      }
    }
  }
  
  if (B_0 >= 2)
  {
    std::vector<IntegerVector> U2 = interpretable_poly_helper(2,B_1,B_inf);
    int U2_size = U2.size();
    
    for (int k = 0; k < U2_size; k++)
    {
      IntegerVector v = U2[k];
      
      for (int i = 0; i < d; i++)
      {
        for (int j = i+1; j < d; j++)
        {
          IntegerVector r1(d);
          r1[i] = v[0];
          r1[j] = v[1];
          U.push_back(clone(r1));
        }
      }
    }
  }
  
  if (B_0 >= 3)
  {
    std::vector<IntegerVector> U3 = interpretable_poly_helper(3,B_1,B_inf);
    int U3_size = U3.size();
    
    for (int k = 0; k < U3_size; k++)
    {
      IntegerVector v = U3[k];
      
      for (int i = 0; i < d; i++)
      {
        for (int j = i+1; j < d; j++)
        {
          for (int m = j+1; m < d; m++)
          {
            IntegerVector r1(d);
            r1[i] = v[0];
            r1[j] = v[1];
            r1[m] = v[2];
            U.push_back(clone(r1));
          }
        }
      }
    }
  }
  if (B_0 >= 4)
  {
    std::vector<IntegerVector> U4 = interpretable_poly_helper(4,B_1,B_inf);
    int U4_size = U4.size();
    
    for (int k = 0; k < U4_size; k++)
    {
      IntegerVector v = U4[k];
      
      for (int i = 0; i < d; i++)
      {
        for (int j = i+1; j < d; j++)
        {
          for (int m = j+1; m < d; m++)
          {
            for (int p = m+1; p < d; p++)
            {
              IntegerVector r1(d);
              r1[i] = v[0];
              r1[j] = v[1];
              r1[m] = v[2];
              r1[p] = v[3];
              U.push_back(clone(r1));
            }
          }
        }
      }
    }
  }
  
  return U;
}

// [[Rcpp::export]]
NumericMatrix gen_unif_mat(int nrow, int ncol)
{
  NumericMatrix unif_mat(nrow, ncol, Rcpp::runif(nrow*ncol).begin());
  return unif_mat;
}

// [[Rcpp::export]]
double artifical_dense_function(const IntegerMatrix& U, NumericVector x)
{
  double adf_out = 0.0;
  
  for (int j = 0; j < U.nrow(); j++)
  {
    adf_out += pow(0.8,sum(U(j,_)))*mulidimensional_basis(U(j,_), x);
  }
  return adf_out;
}

// [[Rcpp::export]]
double artifical_sparse_function(const IntegerMatrix& U, NumericVector x)
{
  double adf_out = 0.0;
  
  for (int j = 0; j < U.nrow(); j++)
  {
    if (sum(U(j,_) > 0) <= 2 && sum(U(j,_)) <= 2 && max(U(j,_)) <= 2)
    {
      adf_out += pow(0.8,sum(U(j,_)))*mulidimensional_basis(U(j,_), x);
    }
  }
  return adf_out;
}


// USED

//[[Rcpp::export]]
NumericVector apply_mean_col(NumericMatrix X)
{
  int numRows = X.nrow();
  NumericVector f_out(numRows);
  
  for(int i = 0; i < numRows; i++)
  {
    f_out[i]=mean(X(i,_));
  }    
  return f_out;
}

// USED

//[[Rcpp::export]]
NumericVector apply_var_col(NumericMatrix X)
{
  int numRows = X.nrow();
  NumericVector f_out(numRows);
  
  for(int i = 0; i < numRows; i++)
  {
    f_out[i]=var(X(i,_));
  }    
  return f_out;
}

// USED

//[[Rcpp::export]]
double sigma_hat_QM(int Q, const std::vector<NumericVector>& f_X_Q)
{
  NumericVector f_sqr_X_Q_mean(Q);
  NumericVector mu_hat_rqmc(Q);
  
  for (int i = 0; i < Q; i++)
  {
    mu_hat_rqmc[i] = mean(f_X_Q[i]);
    f_sqr_X_Q_mean[i] = mean(pow(f_X_Q[i],2));
  }
  
  double cov_mu_hat_QM_prod = 0;
  
  for (int i = 0; i < Q-1; i++)
  {
    for (int j = i+1; j < Q; j++)
    {
      cov_mu_hat_QM_prod += mu_hat_rqmc[i]*mu_hat_rqmc[j];
    }
  }
  
  return mean(f_sqr_X_Q_mean) - (2.0 / (Q*(Q-1))) * cov_mu_hat_QM_prod;
}

///// For Interpreting Black-Box Functions. //////

// [[Rcpp::export]]
std::vector<NumericVector> QR_adap(int n, const IntegerMatrix& U, 
                                   const NumericMatrix& X, const NumericVector& f_X)
{
  std::vector<NumericVector> outList;
  
  int numRows = U.nrow();
  
  NumericVector betas(numRows);
  NumericVector Var_betas(numRows);
  
  NumericVector beta_errors(n);
  
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < numRows; j++)
    {
      double diff_res = f_X[i]*mulidimensional_basis(U(j,_),X(i,_)) - betas[j];
      betas[j] += (1.0/(i+1))*(diff_res);
      Var_betas[j] = (1.0-2.0/(i+1))*Var_betas[j] + 1.0 / ((i+1.0)*(i+1.0)) * diff_res * diff_res;
    }
    
    double k = log2(i+1);
    
    if (k == floor(k) && k > 9)
    {
      outList.push_back(Rcpp::clone(betas));
      outList.push_back(Rcpp::clone(Var_betas));
    }
  }
  
  return outList;
}

// [[Rcpp::export]]
std::vector<NumericVector> QR_email(int n, const IntegerMatrix& U, 
                                   const NumericMatrix& X, const NumericVector& f_X)
{
  std::vector<NumericVector> outList;
  
  int numRows = U.nrow();
  
  NumericVector betas(numRows);
  NumericVector Var_betas(numRows);
  
  NumericVector beta_errors(n);
  
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < numRows; j++)
    {
      double diff_res = f_X[i]*mulidimensional_basis(U(j,_),X(i,_)) - betas[j];
      betas[j] += (1.0/(i+1))*(diff_res);
      Var_betas[j] = (1.0-2.0/(i+1))*Var_betas[j] + 1.0 / ((i+1.0)*(i+1.0)) * diff_res * diff_res;
    }
    
    double k = log2(i+1);
    
    if (k == floor(k) && k > 9)
    {
      outList.push_back(Rcpp::clone(betas));
      outList.push_back(Rcpp::clone(Var_betas));
    }
  }

  //outList.push_back(Rcpp::clone(betas));
  //outList.push_back(Rcpp::clone(Var_betas));
  
  return outList;
}

// [[Rcpp::export]]
std::vector<NumericVector> QRS_adap(int n, const IntegerMatrix& U,
                                    const NumericMatrix& X, const NumericVector& f_X, 
                                    int burn_in = 500, bool oracle_shrinkage=true)
{
  std::vector<NumericVector> outList;
  
  int numRows = U.nrow();
  
  NumericVector betas(numRows);
  NumericVector Var_betas(numRows);
  NumericVector gammas(numRows);
  
  // Burn in value.
  int B = burn_in;
  
  for (int i = 0; i < n; i++)
  {
    double shrinkage_total = 0;
    
    NumericVector psi_basis(numRows);
    
    for (int j=0; j < numRows; j++)
    {
      psi_basis[j] = mulidimensional_basis(U(j,_),X(i,_));
    }
    
    // Burn-in of 500.
    if (i >= B)
    {
      for (int j=0; j < numRows; j++)
      {
        shrinkage_total += gammas[j]*betas[j]*psi_basis[j];
      }
    }
    
    for (int j = 0; j < numRows; j++)
    {
      double shrinkage_minus_j = 0;
      if (i >= B)
      {
        shrinkage_minus_j = shrinkage_total -
          gammas[j]*betas[j]*psi_basis[j];
      }
      
      double T_ij = psi_basis[j]*(f_X[i] - shrinkage_minus_j);
      Var_betas[j] = (1-2.0/(i+1))*Var_betas[j] + 1.0/pow(i+1,2)*pow(T_ij - betas[j],2);
      betas[j] = (1.0/(i+1))*(i*betas[j] + T_ij);
    } 
    
    if (i >= B)
    {
      double med_var_betas = median(Var_betas);
      
      if (oracle_shrinkage)
      {
        for (int j = 0; j < numRows; j++)
        {
          gammas[j] = pow(betas[j],2) / (pow(betas[j],2) + med_var_betas);
        }
      }
      else
      {
        double ti = sqrt(2.0*log(numRows)*med_var_betas);
        for (int j = 0; j < numRows; j++)
        {
          gammas[j] = static_cast<double>(fabs(betas[j]) > ti);
        }
      }
    }
    
    double k = log2(i+1);
    
    if (k == floor(k) && k > 9)
    {
      outList.push_back(Rcpp::clone(betas));
      outList.push_back(Rcpp::clone(Var_betas));
      outList.push_back(Rcpp::clone(gammas));
    }
  }
  
  return outList;
}

// [[Rcpp::export]]
std::vector<NumericVector> QRS_email(int n, const IntegerMatrix& U,
                                    const NumericMatrix& X, const NumericVector& f_X, 
                                    int burn_in = 500, bool oracle_shrinkage=true)
{
  std::vector<NumericVector> outList;
  
  int numRows = U.nrow();
  
  NumericVector betas(numRows);
  NumericVector Var_betas(numRows);
  NumericVector gammas(numRows);
  
  // Burn in value.
  int B = burn_in;
  
  for (int i = 0; i < n; i++)
  {
    double shrinkage_total = 0;
    
    NumericVector psi_basis(numRows);
    
    for (int j=0; j < numRows; j++)
    {
      psi_basis[j] = mulidimensional_basis(U(j,_),X(i,_));
    }
    
    // Burn-in of 500.
    if (i >= B)
    {
      for (int j=0; j < numRows; j++)
      {
        shrinkage_total += gammas[j]*betas[j]*psi_basis[j];
      }
    }
    
    for (int j = 0; j < numRows; j++)
    {
      double shrinkage_minus_j = 0;
      if (i >= B)
      {
        shrinkage_minus_j = shrinkage_total -
          gammas[j]*betas[j]*psi_basis[j];
      }
      
      double T_ij = psi_basis[j]*(f_X[i] - shrinkage_minus_j);
      Var_betas[j] = (1-2.0/(i+1))*Var_betas[j] + 1.0/pow(i+1,2)*pow(T_ij - betas[j],2);
      betas[j] = (1.0/(i+1))*(i*betas[j] + T_ij);
    } 
    
    if (i >= B)
    {
      double med_var_betas = median(Var_betas);
      
      if (oracle_shrinkage)
      {
        for (int j = 0; j < numRows; j++)
        {
          gammas[j] = pow(betas[j],2) / (pow(betas[j],2) + med_var_betas);
        }
      }
      else
      {
        double ti = sqrt(2.0*log(numRows)*med_var_betas);
        for (int j = 0; j < numRows; j++)
        {
          gammas[j] = static_cast<double>(fabs(betas[j]) > ti);
        }
      }
    }
    
    double k = log2(i+1);
    
    if (k == floor(k) && k > 9)
    {
      outList.push_back(Rcpp::clone(betas));
      outList.push_back(Rcpp::clone(Var_betas));
      outList.push_back(Rcpp::clone(gammas));
    }
  }
  
  //outList.push_back(Rcpp::clone(betas));
  //outList.push_back(Rcpp::clone(Var_betas));
  //outList.push_back(Rcpp::clone(gammas));
  
  return outList;
}

// [[Rcpp::export]]
std::vector<NumericVector> QRS_email_P(int n, const IntegerMatrix& U,
                                     const NumericMatrix& X, const NumericVector& f_X, 
                                     int burn_in = 500, bool oracle_shrinkage=true)
{
  std::vector<NumericVector> outList;
  
  int numRows = U.nrow();
  
  NumericVector betas(numRows);
  NumericVector Var_betas(numRows);
  NumericVector gammas(numRows);
  
  // Burn in value.
  int B = burn_in;
  
  for (int i = 0; i < n; i++)
  {
    double shrinkage_total = 0;
    
    NumericVector psi_basis(numRows);
    
    for (int j=0; j < numRows; j++)
    {
      psi_basis[j] = mulidimensional_basis(U(j,_),X(i,_));
    }
    
    // Burn-in of 500.
    if (i >= B)
    {
      for (int j=0; j < numRows; j++)
      {
        shrinkage_total += gammas[j]*betas[j]*psi_basis[j];
      }
    }
    
    for (int j = 0; j < numRows; j++)
    {
      double shrinkage_minus_j = 0;
      if (i >= B)
      {
        shrinkage_minus_j = shrinkage_total -
          gammas[j]*betas[j]*psi_basis[j];
      }
      
      double T_ij = psi_basis[j]*(f_X[i] - shrinkage_minus_j);
      Var_betas[j] = (1-2.0/(i+1))*Var_betas[j] + 1.0/pow(i+1,2)*pow(T_ij - betas[j],2);
      betas[j] = (1.0/(i+1))*(i*betas[j] + T_ij);
    } 
    
    if (i >= B)
    {
      double med_var_betas = median(Var_betas);
      
      if (oracle_shrinkage)
      {
        for (int j = 0; j < numRows; j++)
        {
          gammas[j] = pow(betas[j],2) / (pow(betas[j],2) + med_var_betas);
        }
      }
      else
      {
        double ti = sqrt(2.0*log(numRows)*med_var_betas);
        for (int j = 0; j < numRows; j++)
        {
          gammas[j] = static_cast<double>(fabs(betas[j]) > ti);
        }
      }
    }
  }
  
  outList.push_back(Rcpp::clone(betas));
  outList.push_back(Rcpp::clone(Var_betas));
  outList.push_back(Rcpp::clone(gammas));
  
  return outList;
}

// [[Rcpp::export]]
std::vector<NumericVector> QR_email_P(int n, const IntegerMatrix& U, 
                                    const NumericMatrix& X, const NumericVector& f_X)
{
  std::vector<NumericVector> outList;
  
  int numRows = U.nrow();
  
  NumericVector betas(numRows);
  NumericVector Var_betas(numRows);
  
  NumericVector beta_errors(n);
  
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < numRows; j++)
    {
      double diff_res = f_X[i]*mulidimensional_basis(U(j,_),X(i,_)) - betas[j];
      betas[j] += (1.0/(i+1))*(diff_res);
      Var_betas[j] = (1.0-2.0/(i+1))*Var_betas[j] + 1.0 / ((i+1.0)*(i+1.0)) * diff_res * diff_res;
    }
  }
  
  outList.push_back(Rcpp::clone(betas));
  outList.push_back(Rcpp::clone(Var_betas));
  
  return outList;
}

// [[Rcpp::export]]
NumericVector LOF_EST_cpp(int d, const NumericMatrix& X, NumericVector f_X, double var_f, IntegerMatrix U, NumericVector betas, NumericVector gamma, int n) {
  NumericVector Int_sqr_err(n);
  
  int numRow = U.nrow();
  
  for (int i = 0; i < n; i++) {
    double fi_X_i = 0.0;
    for (int j = 0; j < numRow; j++) {
      fi_X_i += gamma[j] * betas[j] * mulidimensional_basis(U(j,_),X(i,_));
    }
    Int_sqr_err[i] = pow(f_X[i] - fi_X_i, 2);
  }
  
  return Int_sqr_err;
}

///// Comparison to true beta. //////

// [[Rcpp::export]]
std::vector<NumericVector> QR_true(int n, const IntegerMatrix& U, 
                                   const NumericMatrix& X, const NumericVector& f_X,
                                   const NumericVector& true_beta)
{
  std::vector<NumericVector> outList;
  
  int numRows = U.nrow();
  
  NumericVector betas(numRows);
  NumericVector Var_betas(numRows);
  
  NumericVector beta_errors(n);
  
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < numRows; j++)
    {
      double diff_res = f_X[i]*mulidimensional_basis(U(j,_),X(i,_)) - betas[j];
      betas[j] += (1.0/(i+1))*(diff_res);
      Var_betas[j] = (1.0-2.0/(i+1))*Var_betas[j] + 1.0 / ((i+1.0)*(i+1.0)) * diff_res * diff_res;
      double beta_diff = (true_beta[j] != 0) ? (betas[j] - true_beta[j])*(betas[j] - true_beta[j]) / (true_beta[j]*true_beta[j]) : betas[j]*betas[j];
      beta_errors[i] = (1.0/(j+1))*(j*beta_errors[i] + beta_diff);
    }
    
    double k = log2(i+1);
    
    if (k == floor(k) && k > 9)
    {
      outList.push_back(Rcpp::clone(betas));
      outList.push_back(Rcpp::clone(Var_betas));
    }
  }
  
  outList.push_back(beta_errors);
  
  return outList;
}

// Code up deep shrinkage.

// USED

// [[Rcpp::export]]
std::vector<NumericVector> QRS_true(int n, const IntegerMatrix& U,
                                    const NumericMatrix& X, const NumericVector& f_X,
                                    const NumericVector& true_beta, int burn_in = 500, bool oracle_shrinkage=true)
{
  std::vector<NumericVector> outList;
  
  int numRows = U.nrow();
  
  NumericVector betas(numRows);
  NumericVector Var_betas(numRows);
  NumericVector gammas(numRows);
  
  NumericVector beta_errors(n);
  
  // Burn in value.
  int B = burn_in;
  
  for (int i = 0; i < n; i++)
  {
    double shrinkage_total = 0;
    
    NumericVector psi_basis(numRows);
    
    for (int j=0; j < numRows; j++)
    {
      psi_basis[j] = mulidimensional_basis(U(j,_),X(i,_));
    }
    
    // Burn-in of 500.
    if (i >= B)
    {
      for (int j=0; j < numRows; j++)
      {
        shrinkage_total += gammas[j]*betas[j]*psi_basis[j];
      }
    }
    
    for (int j = 0; j < numRows; j++)
    {
      double shrinkage_minus_j = 0;
      if (i >= B)
      {
        shrinkage_minus_j = shrinkage_total -
          gammas[j]*betas[j]*psi_basis[j];
      }
      
      double T_ij = psi_basis[j]*(f_X[i] - shrinkage_minus_j);
      Var_betas[j] = (1-2.0/(i+1))*Var_betas[j] + 1.0/pow(i+1,2)*pow(T_ij - betas[j],2);
      betas[j] = (1.0/(i+1))*(i*betas[j] + T_ij);
      double beta_diff = (true_beta[j] != 0) ? (betas[j] - true_beta[j])*(betas[j] - true_beta[j]) / (true_beta[j]*true_beta[j]) : betas[j]*betas[j];
      beta_errors[i] = (1.0/(j+1))*(j*beta_errors[i] + beta_diff);
    } 
    
    if (i >= B)
    {
      double med_var_betas = median(Var_betas);
      
      if (oracle_shrinkage)
      {
        for (int j = 0; j < numRows; j++)
        {
          gammas[j] = pow(betas[j],2) / (pow(betas[j],2) + med_var_betas);
        }
      }
      else
      {
        double ti = sqrt(2.0*log(numRows)*med_var_betas);
        for (int j = 0; j < numRows; j++)
        {
          gammas[j] = static_cast<double>(fabs(betas[j]) > ti);
        }
      }
    }
    
    double k = log2(i+1);
    
    if (k == floor(k) && k > 9)
    {
      outList.push_back(Rcpp::clone(betas));
      outList.push_back(Rcpp::clone(Var_betas));
      outList.push_back(Rcpp::clone(gammas));
    }
  }
  
  outList.push_back(beta_errors);
  
  return outList;
}

// [[Rcpp::export]]
std::vector<NumericVector> OQRS_true(int n, const IntegerMatrix& U, int Q, 
                                     const std::vector<NumericMatrix>& X_Q, const std::vector<NumericVector>& f_X_Q,
                                     const NumericVector& true_beta, int burn_in=500, bool oracle_shrinkage=true)
{
  std::vector<NumericVector> outList;
  
  int numRows = U.nrow();
  
  NumericMatrix betas_Q(numRows, Q);
  NumericVector gammas(numRows);
  
  NumericMatrix beta_errors(n,Q);
  
  // Burn in value.
  int B = burn_in;
  
  for (int i = 0; i < n; i++)
  {
    for (int k = 0; k < Q; k++)
    {
      double shrinkage_total = 0;
      
      NumericVector psi_basis(numRows);
      
      for (int j = 0; j < numRows; j++)
      {
        psi_basis[j] = mulidimensional_basis(U(j,_),X_Q[k](i,_));
      }
      
      // Burn-in of 500.
      if (i >= B)
      {
        for (int j=0; j < numRows; j++)
        {
          shrinkage_total += gammas[j]*betas_Q(j,k)*psi_basis[j];
        }
      }
      
      for (int j = 0; j < numRows; j++)
      {
        double shrinkage_minus_j = 0;
        
        if (i >= B)
        {
          shrinkage_minus_j = shrinkage_total -
            gammas[j]*betas_Q(j,k)*psi_basis[j];
        }
        
        double T_ij = psi_basis[j]*(f_X_Q[k][i] - shrinkage_minus_j);
        betas_Q(j,k) = (1.0/(i+1))*(i*betas_Q(j,k) + T_ij);
        double beta_diff = (true_beta[j] != 0) ? (betas_Q(j,k) - true_beta[j])*(betas_Q(j,k) - true_beta[j]) / (true_beta[j]*true_beta[j]) : betas_Q(j,k)*betas_Q(j,k);
        beta_errors(i, k) = (1.0/(j+1))*(j*beta_errors(i, k) + beta_diff);
      }
    }
    
    if (i >= B)
    {
      NumericVector beta_Q_avg = apply_mean_col(betas_Q);
      NumericVector beta_Q_var = apply_var_col(betas_Q);
      double med_beta_Q_var = median(beta_Q_var);
      
      if (oracle_shrinkage)
      {
        for (int j = 0; j < numRows; j++)
        {
          gammas[j] = pow(beta_Q_avg[j],2) / (pow(beta_Q_avg[j],2) + med_beta_Q_var);
        }
      }
      else
      {
        double ti = sqrt(2.0*log(numRows)*med_beta_Q_var);
        for (int j = 0; j < numRows; j++)
        {
          gammas[j] = static_cast<double>(fabs(beta_Q_avg[j]) > ti);
        }
      }
      
      double k = log2(i+1);
      
      if (k == floor(k) && k > 9)
      {
        outList.push_back(Rcpp::clone(beta_Q_avg));
        outList.push_back(Rcpp::clone(beta_Q_var));
        outList.push_back(Rcpp::clone(gammas));
      }
    }
  }

  outList.push_back(apply_mean_col(beta_errors));
  
  return outList;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/

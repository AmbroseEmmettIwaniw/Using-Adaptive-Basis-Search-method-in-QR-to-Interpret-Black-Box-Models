#include <Rcpp.h>
using namespace Rcpp;


// Used a register as in hardware to move through the possible basis efficiently.
// e.g k = 3, (0,1,1,2) -> (0,1,2,0).
// [[Rcpp::export]]
IntegerVector d_register(IntegerVector r, int k) {
  int d = r.size();
  
  for (int i = (d-1); i >= 0; i--)
  {
    if (r[i] < k-1)
    {
      r[i]++;
      
      for (int j = i+1; j < d; j++)
      {
        r[j] = r[i];
      }
      break;
    }
    else if (r[i] == k-1)
    {
      r[i] = 0;
    }
  }
  return r;
}

// The first 6 degree orthonormal polynomial basis.
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

// In ABS we know the basis already so we can efficiently compute the coefficient (don't need to take products of 0 degree basis)
// [[Rcpp::export]]
double mulidim_basis(const IntegerVector& r, const NumericVector& x, 
                             const IntegerVector& r_indices)
{
  double phi_tensor_d = 1;
  int num_non_zero = r_indices.size();
  
  for (int i = 0; i < num_non_zero; i++)
  {
    phi_tensor_d *= orthonormal_basis(x[r_indices[i]],r[r_indices[i]]);
  }
  
  return phi_tensor_d;
}

// Multidimensional basis of orthonomial polynomials. r is basis index and x is input [0,1]^d.
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

// Faster way in Rcpp to generate a uniform matrix of points.
// [[Rcpp::export]]
NumericMatrix gen_unif_mat(int nrow, int ncol)
{
  NumericMatrix unif_mat(nrow, ncol, Rcpp::runif(nrow*ncol).begin());
  return unif_mat;
}

// This computes ...
// [[Rcpp::export]]
std::vector<IntegerVector> poly_interpretable_helper(int d, int B_1, int B_inf)
{
  int max_r = B_inf + 1;
  
  std::vector<IntegerVector> U;
  IntegerVector r(d);
  
  for (int i = 0; i < pow(max_r,d); i++)
  {
    int r_1_norm = sum(r);
    
    if (r_1_norm <= B_1 && sum(r==0) == 0)
    {
      U.push_back(clone(r));
    }
    
    if (r_1_norm == d*B_inf)
      break;
    
    r = clone(d_register(r,max_r));
  }
  
  return U;
}

// [[Rcpp::export]]
std::vector<IntegerVector> non_zero_indices_pairs(const IntegerVector& r)
{
  std::vector<IntegerVector> indices_pairs;
  
  for (int i = 0; i < r.size(); i++)
  {
    if (r[i] > 0)
    {
      IntegerVector pair(2);
      pair[0] = i;
      pair[1] = r[i];
      indices_pairs.push_back(pair);
    }
  }
  
  return indices_pairs;
}

//[[Rcpp::export]] 
int IsSubsetSet(IntegerVector sub, IntegerVector sup)
{
  for (int i = 0; i < sub.size(); i++)
  {
    for (int j = 0; j < sup.size(); j++)
    {
      if (sub[i] == sup[j])
      {
        sup.erase(j);
        break;
      }
    }
  }
  
  if (sup.size() == 1)
    return sup[0];
  else
    return -1;
}

// [[Rcpp::export]]
std::vector<IntegerVector> Pot_inter_poly_fcts(IntegerVector r,
                                                int B_1, int B_inf)
{
  int max_r = B_inf;
  int d = r.size();
  
  std::vector<IntegerVector> U;
  
  if (sum(r>0) == 0) // Single factors.
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
  else if (sum(r>0) == 1) // two-factors
  {
    IntegerVector nz_pair = non_zero_indices_pairs(r)[0];
    std::vector<IntegerVector> U2 = poly_interpretable_helper(2,B_1,B_inf);
    int U2_size = U2.size();
    
    for (int k = 0; k < U2_size; k++)
    {
      IntegerVector v(1,nz_pair[1]);
      
      int subset_val = IsSubsetSet(v, U2[k]);
      
      if (subset_val != -1)
      {
        for (int i = 0; i < d; i++)
        {
          if (i != nz_pair[0])
          {
            IntegerVector r1(d);
            r1[i] = subset_val;
            r1[nz_pair[0]] = nz_pair[1];
            U.push_back(clone(r1));
          }
        }
      }
    }
  }
  else if (sum(r>0) == 2) // three-factors
  {
    // TODO: For 3-factors.
    std::vector<IntegerVector> nz_pairs = non_zero_indices_pairs(r);
    std::vector<IntegerVector> U3 = poly_interpretable_helper(3,B_1,B_inf);
    int U3_size = U3.size();
    
    for (int k = 0; k < U3_size; k++)
    {
      IntegerVector v = IntegerVector::create(nz_pairs[0][1], nz_pairs[1][1]);
      
      int subset_val = IsSubsetSet(v, U3[k]);
      
      if (subset_val != -1)
      {
        for (int i = 0; i < d; i++)
        {
          if (i != nz_pairs[0][0] && i != nz_pairs[1][0])
          {
            IntegerVector r1(d);
            r1[i] = subset_val;
            r1[nz_pairs[0][0]] = nz_pairs[0][1];
            r1[nz_pairs[1][0]] = nz_pairs[1][1];
            U.push_back(clone(r1));
          }
        }
      }
    }
  }
  else if (sum(r>0) == 3) // three-factors
  {
    std::vector<IntegerVector> nz_pairs = non_zero_indices_pairs(r);
    std::vector<IntegerVector> U4 = poly_interpretable_helper(4,B_1,B_inf);
    int U4_size = U4.size();
    
    for (int k = 0; k < U4_size; k++)
    {
      IntegerVector v = IntegerVector::create(nz_pairs[0][1], nz_pairs[1][1], nz_pairs[2][1]);
      
      int subset_val = IsSubsetSet(v, U4[k]);
      
      if (subset_val != -1)
      {
        for (int i = 0; i < d; i++)
        {
          if (i != nz_pairs[0][0] && i != nz_pairs[1][0] && i != nz_pairs[2][0])
          {
            IntegerVector r1(d);
            r1[i] = subset_val;
            r1[nz_pairs[0][0]] = nz_pairs[0][1];
            r1[nz_pairs[1][0]] = nz_pairs[1][1];
            r1[nz_pairs[2][0]] = nz_pairs[2][1];
            U.push_back(clone(r1));
          }
        }
      }
    }
  }
  
  return U;
}

//////// Haar Basis. 

// [[Rcpp::export]]
IntegerVector d_register_haar(IntegerVector r, int k) {
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

//[[Rcpp::export]]
double haar_wavelet_01(double x)
{
  if (x < 0)
    return 0;
  else if (x <= 0.5)
    return -1;
  else if (x <= 1)
    return 1;
  else
    return 0;
}

//[[Rcpp::export]]
double haar_wavelet_basis(double x, int k)
{
  if (k == 0)
    return 1;
  else
  {
    double tk = floor(log(k) / log(2)) + 1;
    double lk = k - pow(2,tk-1);
    
    return pow(2, (tk-1)/2.0)*haar_wavelet_01(pow(2,tk-1)*x-lk);
  }
}

//[[Rcpp::export]]
NumericVector haar_deg(IntegerVector r, int max_deg=0)
{
  NumericVector result(r.size());
  
  for (int i = 0; i < r.size(); i++)
  {
    result[i] = (r[i] != 0) ? (floor(log(r[i]) / log(2)) + 1) : max_deg;
  }
  return result;
}

// [[Rcpp::export]]
double mulidim_haar_basis(const IntegerVector& r, const NumericVector& x, 
                          const IntegerVector& r_indices)
{
  double phi_tensor_d = 1;
  int num_non_zero = r_indices.size();
  
  for (int i = 0; i < num_non_zero; i++)
  {
    phi_tensor_d *= haar_wavelet_basis(x[r_indices[i]],r[r_indices[i]]);
  }
  
  return phi_tensor_d;
}

// [[Rcpp::export]]
double mulidimensional_haar_basis(const IntegerVector& r, const NumericVector& x)
{
  double phi_tensor_d = 1;
  int d = r.size();
  
  for (int i = 0; i < d; i++)
  {
    phi_tensor_d *= haar_wavelet_basis(x[i],r[i]);
  }
  
  return phi_tensor_d;
}

// [[Rcpp::export]]
std::vector<IntegerVector> haar_interpretable_helper(int d, int B_1, int B_inf)
{
  int max_r = pow(2,B_inf);
  int num_haar_basis = max_r;
  
  if (d == 2)
    num_haar_basis = (max_r + 1) * max_r / 2;
  else if (d==3)
    num_haar_basis =  (max_r + 1) * max_r / 4 + (max_r + 1)*max_r*(2*max_r+1) / 12;
  
  std::vector<IntegerVector> U;
  IntegerVector r(d);
  
  for (int i = 0; i < num_haar_basis; i++)
  {
    if (sum(haar_deg(r, B_1)) <= B_1)
      U.push_back(clone(r));
    
    r = clone(d_register_haar(r,max_r));
  }
  
  return U;
}

// [[Rcpp::export]]
std::vector<IntegerVector> interpretable_haar_helper(int d, int B_1, int B_inf)
{
  int max_r = pow(2,B_inf);
  
  std::vector<IntegerVector> U;
  IntegerVector r(d);
  
  for (int i = 0; i < pow(max_r,d); i++)
  {
    if (sum(haar_deg(r, B_1)) <= B_1)
    {
      U.push_back(clone(r));
    }
    r = clone(d_register_haar(r,max_r));
  }
  
  return U;
}

// [[Rcpp::export]]
std::vector<IntegerVector> interpretable_haar_sets(int d, int B_0, int B_1, int B_inf)
{
  int max_r = pow(2,B_inf)-1;
  
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
    std::vector<IntegerVector> U2 = interpretable_haar_helper(2,B_1,B_inf);
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
    std::vector<IntegerVector> U3 = interpretable_haar_helper(3,B_1,B_inf);
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
  
  return U;
}

// [[Rcpp::export]]
std::vector<IntegerVector> Pot_inter_haar_fcts(IntegerVector r,
                                               int B_1, int B_inf)
{
  int max_r = pow(2,B_inf);
  int d = r.size();
  
  std::vector<IntegerVector> U;
  
  if (sum(r>0) == 0) // Single factors.
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
  else if (sum(r>0) == 1)
  {
    IntegerVector nz_pair = non_zero_indices_pairs(r)[0];
    std::vector<IntegerVector> U2 = haar_interpretable_helper(2,B_1,B_inf);
    int U2_size = U2.size();
    
    for (int k = 0; k < U2_size; k++)
    {
      IntegerVector v(1,nz_pair[1]);
      
      int subset_val = IsSubsetSet(v, U2[k]);
      
      if (subset_val != -1)
      {
        for (int i = 0; i < d; i++)
        {
          if (i != nz_pair[0])
          {
            IntegerVector r1(d);
            r1[i] = subset_val;
            r1[nz_pair[0]] = nz_pair[1];
            U.push_back(clone(r1));
          }
        }
      }
    }
  }
  else if (sum(r>0) == 2) //  For 3-factors.
  {
    std::vector<IntegerVector> nz_pairs = non_zero_indices_pairs(r);
    std::vector<IntegerVector> U3 = haar_interpretable_helper(3,B_1,B_inf);
    int U3_size = U3.size();
    
    for (int k = 0; k < U3_size; k++)
    {
      IntegerVector v = IntegerVector::create(nz_pairs[0][1], nz_pairs[1][1]);
      
      int subset_val = IsSubsetSet(v, U3[k]);
      
      if (subset_val != -1)
      {
        for (int i = 0; i < d; i++)
        {
          if (i != nz_pairs[0][0] && i != nz_pairs[1][0])
          {
            IntegerVector r1(d);
            r1[i] = subset_val;
            r1[nz_pairs[0][0]] = nz_pairs[0][1];
            r1[nz_pairs[1][0]] = nz_pairs[1][1];
            U.push_back(clone(r1));
          }
        }
      }
    }
  }
  else if (sum(r>0) == 3) //  For 4-factors.
  {
    std::vector<IntegerVector> nz_pairs = non_zero_indices_pairs(r);
    std::vector<IntegerVector> U4 = haar_interpretable_helper(4,B_1,B_inf);
    int U4_size = U4.size();
    
    for (int k = 0; k < U4_size; k++)
    {
      IntegerVector v = IntegerVector::create(nz_pairs[0][1], nz_pairs[1][1], nz_pairs[2][1]);
      
      int subset_val = IsSubsetSet(v, U4[k]);
      
      if (subset_val != -1)
      {
        for (int i = 0; i < d; i++)
        {
          if (i != nz_pairs[0][0] && i != nz_pairs[1][0] && i != nz_pairs[2][0])
          {
            IntegerVector r1(d);
            r1[i] = subset_val;
            r1[nz_pairs[0][0]] = nz_pairs[0][1];
            r1[nz_pairs[1][0]] = nz_pairs[1][1];
            r1[nz_pairs[2][0]] = nz_pairs[2][1];
            U.push_back(clone(r1));
          }
        }
      }
    }
  }
  
  return U;
}

// Compute the beta coefficients.

// [[Rcpp::export]]
IntegerVector non_zero_indices(const IntegerVector& r)
{
  IntegerVector indices;
  
  for (int i = 0; i < r.size(); i++)
  {
    if (r[i] > 0)
      indices.push_back(i);
  }
  
  return indices;
}

// [[Rcpp::export]]
double beta_r(int n, const IntegerVector& r, 
              const NumericMatrix& X, const NumericVector& f_X)
{
  double beta_est=0;
  IntegerVector r_NZ_ind = non_zero_indices(r);
  
  for (int i = 0; i < n; i++)
  {
    double diff_res = f_X[i]*mulidim_basis(r,X(i,_), r_NZ_ind) - beta_est;
    beta_est += (1.0/(i+1))*(diff_res);
  }
  
  return beta_est;
}

// [[Rcpp::export]]
std::vector<double> fit_beta_r_inactive(int n, const std::vector<IntegerVector>& new_inactive_r, 
                         const NumericMatrix& X, const NumericVector& f_X)
{
  std::vector<double> betas;
  int inactive_size = new_inactive_r.size();
  
  for (int i = 0; i < inactive_size; i++)
  {
    betas.push_back(beta_r(n, new_inactive_r[i], X, f_X));
  }
  
  return betas;
}

// [[Rcpp::export]]
double haar_beta_r(int n, const IntegerVector& r, 
              const NumericMatrix& X, const NumericVector& f_X)
{
  double beta_est=0;
  IntegerVector r_NZ_ind = non_zero_indices(r);
  
  for (int i = 0; i < n; i++)
  {
    double diff_res = f_X[i]*mulidim_haar_basis(r,X(i,_), r_NZ_ind) - beta_est;
    beta_est += (1.0/(i+1))*(diff_res);
  }
  
  return beta_est;
}

// [[Rcpp::export]]
std::vector<double> haar_beta_r_inactive(int n, const std::vector<IntegerVector>& new_inactive_r, 
                                        const NumericMatrix& X, const NumericVector& f_X)
{
  std::vector<double> betas;
  int inactive_size = new_inactive_r.size();
  
  for (int i = 0; i < inactive_size; i++)
  {
    betas.push_back(haar_beta_r(n, new_inactive_r[i], X, f_X));
  }
  
  return betas;
}

// [[Rcpp::export]]
double ResidualErrorEstimate(int n, const std::vector<IntegerVector>& activeSet,
                             const std::vector<double>& betas_active, 
                             const IntegerVector& r_pot, double beta_pot,
                             const NumericMatrix& X, const NumericVector& f_X)
{
  double resid_error = 0;
  
  for (int i = 0; i < n; i++)
  {
    double active_res = 0;
    int active_size = activeSet.size();
    
    for (int j = 0; j < active_size; j++)
    {
      active_res += betas_active[j] * mulidimensional_basis(activeSet[j], X(i,_));
    }
    
    double res_diff = f_X[i] - active_res - 
      beta_pot*mulidimensional_basis(r_pot, X(i,_)) - resid_error;
    
    resid_error += (1.0/(i+1))*res_diff;
  }
  
  return resid_error;
}

int argmin_REE(int n, const std::vector<IntegerVector>& activeSet,
               const std::vector<double>& betas_active, 
               const std::vector<IntegerVector>& inactiveSet,
               const std::vector<double>& betas_inactive,
               const NumericMatrix& X, const NumericVector& f_X)
{
  int num_pot_r = inactiveSet.size();
  NumericVector REE_ests(num_pot_r);
  
  for (int i = 0; i < num_pot_r; i++)
  {
    REE_ests[i] = ResidualErrorEstimate(n, activeSet, betas_active,
                                        inactiveSet[i], betas_inactive[i],
                                        X, f_X);
  }
  
  return which_min(REE_ests);
}


// [[Rcpp::export]]
std::vector<IntegerVector> MAQR(int n, int B_0, int B_1, int B_inf, 
                              const NumericMatrix& X, const NumericVector& f_X,
                              size_t M=100)
{
  int d = X.ncol();
  
  std::vector<IntegerVector> activeSet;
  std::vector<IntegerVector> inactiveSet;
  std::vector<double> betas_active;
  std::vector<double> betas_inactive;
  
  // Stage 0: Fit the intercept.
  IntegerVector intercept(d,0);
  activeSet.push_back(intercept);
  betas_active.push_back(beta_r(n, intercept, X, f_X));
  
  inactiveSet = Pot_inter_poly_fcts(intercept, B_1, B_inf);
  betas_inactive = fit_beta_r_inactive(n, inactiveSet, X, f_X);
  
  
  while(activeSet.size() < M)
  {
    std::vector<IntegerVector> new_inactiveSet;
    std::vector<double> new_betas_inactive;
    
    int new_r_ind = argmin_REE(n, activeSet, betas_active, 
                               inactiveSet, betas_inactive,
                               X, f_X);
    
    IntegerVector new_r = inactiveSet[new_r_ind];
    double new_beta_r = betas_inactive[new_r_ind];
    
    activeSet.push_back(new_r);
    betas_active.push_back(new_beta_r);
    
    inactiveSet.erase(inactiveSet.begin() + new_r_ind);
    betas_inactive.erase(betas_inactive.begin()+new_r_ind);
    
    if (sum(new_r > 0) < B_0)
    {
      new_inactiveSet = Pot_inter_poly_fcts(new_r, B_1, B_inf);
      new_betas_inactive = fit_beta_r_inactive(n, new_inactiveSet, X, f_X);
      
      inactiveSet.insert(inactiveSet.end(), new_inactiveSet.begin(), new_inactiveSet.end());
      betas_inactive.insert(betas_inactive.end(), new_betas_inactive.begin(), new_betas_inactive.end());
    }
  }
  
  return activeSet;
}

// [[Rcpp::export]]
double SEE(int n, const std::vector<IntegerVector>& activeSet,
           const std::vector<double>& betas_active,
           const NumericMatrix& X, const NumericVector& f_X)
{
  int active_size = activeSet.size();
  double resid_error = 0;
  
  for (int i = 0; i < n; i++)
  {
    double active_res = 0;
    
    for (int j = 0; j < active_size; j++)
    {
      active_res += betas_active[j] * mulidimensional_basis(activeSet[j], X(i,_));
    }
    
    double res_diff = (f_X[i] - active_res)*(f_X[i] - active_res);
    
    resid_error = (1.0/(i+1))*(i*resid_error + res_diff);
  }
  
  return sqrt(resid_error * ((double)n / (n - active_size)));
}

//[[Rcpp::export]]
std::vector<IntegerVector> SetDiff(int d, const std::vector<IntegerVector>& A, 
                                   const std::vector<IntegerVector>& B)
{
  // Return C = A\B.
  std::vector<IntegerVector> C;
  
  int A_size = A.size();
  int B_size = B.size();
  
  for (int i = 0; i < A_size; i++)
  {
    bool inBFlag = false;
    
    for (int j = 0; j < B_size; j++)
    {
      if (sum(B[j] == A[i])  == d) 
      {
        inBFlag = true;
        break;
      }
    }
    
    if (!inBFlag)
    {
      C.push_back(A[i]);
    }
  }
  
  return C;
}

// [[Rcpp::export]]
std::vector<IntegerVector> MAQRNS(int n, int B_0, int B_1, int B_inf, 
                                const NumericMatrix& X, const NumericVector& f_X,
                                double alpha = 0.05)
{
  int d = X.ncol();
  double norm_quant = R::qnorm(1-alpha/2,0,1,true,false);
  
  std::vector<IntegerVector> activeSet;
  std::vector<double> betas_active;
  
  // Stage 0: Fit the intercept.
  IntegerVector intercept(d,0);
  activeSet.push_back(intercept);
  betas_active.push_back(beta_r(n, intercept, X, f_X));
  
  for(int i = 0; i < B_0; i++)
  {
    std::vector<IntegerVector> new_activeSet;
    std::vector<double> new_betas_active;
    
    int active_size = activeSet.size();
    
    for (int j = 0; j < active_size; j++)
    {
      if (sum(activeSet[j] > 0) == i)
      {
        new_activeSet = Pot_inter_poly_fcts(activeSet[j], B_1, B_inf);
        new_activeSet = SetDiff(d, new_activeSet, activeSet);
        new_betas_active = fit_beta_r_inactive(n, new_activeSet, X, f_X);
        
        activeSet.insert(activeSet.end(), new_activeSet.begin(), new_activeSet.end());
        betas_active.insert(betas_active.end(), new_betas_active.begin(), new_betas_active.end());
      }
    }
    
    double sigma_est = SEE(n, activeSet, betas_active, X, f_X);
    double noise_bound = norm_quant*sigma_est/sqrt(n);
    
    active_size = activeSet.size();
    
    for (int k = active_size-1; k >= 0; k--)
    {
      if (fabs(betas_active[k]) < noise_bound)
      {
        activeSet.erase(activeSet.begin() + k);
        betas_active.erase(betas_active.begin()+k);
      }
    }
    
  }
  
  return activeSet;
}

// [[Rcpp::export]]
std::vector<IntegerVector> MAQRLS(int n, int B_0, int B_1, int B_inf, 
                                  const NumericMatrix& X, const NumericVector& f_X)
{
  int d = X.ncol();
  
  std::vector<IntegerVector> activeSet;
  std::vector<double> betas_active;
  
  // Stage 0: Fit the intercept.
  IntegerVector intercept(d,0);
  activeSet.push_back(intercept);
  betas_active.push_back(beta_r(n, intercept, X, f_X));
  
  for(int i = 0; i < B_0; i++)
  {
    std::vector<IntegerVector> new_activeSet;
    std::vector<double> new_betas_active;
    
    int active_size = activeSet.size();
    
    for (int j = 0; j < active_size; j++)
    {
      if (sum(activeSet[j] > 0) == i)
      {
        new_activeSet = Pot_inter_poly_fcts(activeSet[j], B_1, B_inf);
        new_activeSet = SetDiff(d, new_activeSet, activeSet);
        new_betas_active = fit_beta_r_inactive(n, new_activeSet, X, f_X);
        
        activeSet.insert(activeSet.end(), new_activeSet.begin(), new_activeSet.end());
        betas_active.insert(betas_active.end(), new_betas_active.begin(), new_betas_active.end());
      }
    }
    
    double sigma_est = SEE(n, activeSet, betas_active, X, f_X);
    // log(n) or log(number of basis in active set)?
    double shrinkage_bound = sigma_est/sqrt(n) * sqrt(2*log(activeSet.size()));
    Rcpp::Rcerr << shrinkage_bound << "\n";
    
    active_size = activeSet.size();
    
    for (int k = active_size-1; k >= 0; k--)
    {
      if (fabs(betas_active[k]) < shrinkage_bound)
      {
        activeSet.erase(activeSet.begin() + k);
        betas_active.erase(betas_active.begin()+k);
      }
    }
    
  }
  
  return activeSet;
}

// [[Rcpp::export]]
double quantile_rcpp(const std::vector<double> x, double prob)
{
  std::vector<double> y = x;
  std::sort(y.begin(), y.end());
  return x[ceil(y.size()*prob)-1];
}

// [[Rcpp::export]]
std::vector<IntegerVector> MAQRQS(int n, int B_0, int B_1, int B_inf, 
                                  const NumericMatrix& X, const NumericVector& f_X,
                                  double alpha = 0.05)
{
  int d = X.ncol();
  
  std::vector<IntegerVector> activeSet;
  std::vector<double> betas_active;
  
  // Stage 0: Fit the intercept.
  IntegerVector intercept(d,0);
  activeSet.push_back(intercept);
  betas_active.push_back(beta_r(n, intercept, X, f_X));
  
  for(int i = 0; i < B_0; i++)
  {
    std::vector<IntegerVector> new_activeSet;
    std::vector<double> new_betas_active;
    
    int active_size = activeSet.size();
    
    for (int j = 0; j < active_size; j++)
    {
      if (sum(activeSet[j] > 0) == i)
      {
        new_activeSet = Pot_inter_poly_fcts(activeSet[j], B_1, B_inf);
        new_activeSet = SetDiff(d, new_activeSet, activeSet);
        new_betas_active = fit_beta_r_inactive(n, new_activeSet, X, f_X);
        
        activeSet.insert(activeSet.end(), new_activeSet.begin(), new_activeSet.end());
        betas_active.insert(betas_active.end(), new_betas_active.begin(), new_betas_active.end());
      }
    }
    
    double noise_bound = quantile_rcpp(betas_active, alpha);
    
    active_size = activeSet.size();
    
    for (int k = active_size-1; k >= 0; i--)
    {
      if (fabs(betas_active[k]) < noise_bound)
      {
        activeSet.erase(activeSet.begin() + k);
        betas_active.erase(betas_active.begin()+k);
      }
    }
    
  }
  
  return activeSet;
}

// [[Rcpp::export]]
std::vector<IntegerVector> MAQRSA(int n, int B_0, int B_1, int B_inf, 
                                  const NumericMatrix& X, const NumericVector& f_X,
                                  double alpha = 0.05)
{
  int d = X.ncol();
  
  std::vector<IntegerVector> activeSet;
  std::vector<double> betas_active;
  
  // Stage 0: Fit the intercept.
  IntegerVector intercept(d,0);
  activeSet.push_back(intercept);
  betas_active.push_back(beta_r(n, intercept, X, f_X));
  
  for(int i = 0; i < B_0; i++)
  {
    std::vector<IntegerVector> new_activeSet;
    std::vector<double> new_betas_active;
    
    int active_size = activeSet.size();
    
    for (int j = 0; j < active_size; j++)
    {
      if (sum(activeSet[j] > 0) == i)
      {
        new_activeSet = Pot_inter_poly_fcts(activeSet[j], B_1, B_inf);
        new_activeSet = SetDiff(d, new_activeSet, activeSet);
        new_betas_active = fit_beta_r_inactive(n, new_activeSet, X, f_X);
        
        activeSet.insert(activeSet.end(), new_activeSet.begin(), new_activeSet.end());
        betas_active.insert(betas_active.end(), new_betas_active.begin(), new_betas_active.end());
      }
    }
    
    std::vector<IntegerVector> imp_activeSet;
    imp_activeSet.push_back(intercept);
    std::vector<double> imp_betas_active;
    imp_betas_active.push_back(betas_active[0]);
    std::vector<double> squ_betas;
    double total_variability = 0;
    double curr_variability = 0;
    
    active_size = activeSet.size();
    
    for (int k = 0; k < active_size; k++)
    {
      if (k == 0)
      {
        squ_betas.push_back(-betas_active[k]*betas_active[k]);
      }
      else
      {
        total_variability += betas_active[k]*betas_active[k];
        squ_betas.push_back(betas_active[k]*betas_active[k]);
      }
    }
    
    for (int k = 0; k < active_size-1; k++)
    {
      int max_index = std::distance(squ_betas.begin(),std::max_element(squ_betas.begin(), squ_betas.end()));
      curr_variability += squ_betas[max_index];
      squ_betas[max_index] = -squ_betas[max_index];
      imp_betas_active.push_back(betas_active[max_index]);
      imp_activeSet.push_back(activeSet[max_index]);
      
      if (curr_variability / total_variability > 1- alpha)
      {
        break;
      }
    }
    activeSet = imp_activeSet;
    betas_active = imp_betas_active;
    
  }
  
  return activeSet;
}

// [[Rcpp::export]]
std::vector<double> beta_r_ests(int n, const IntegerVector& r, 
              const NumericMatrix& X, const NumericVector& f_X)
{
  double beta_est=0;
  double var_beta_est = 0;
  
  IntegerVector r_NZ_ind = non_zero_indices(r);
  
  for (int i = 0; i < n; i++)
  {
    double diff_res = f_X[i]*mulidim_basis(r,X(i,_), r_NZ_ind) - beta_est;
    beta_est += (1.0/(i+1))*(diff_res);
    if (i >= 1)
      var_beta_est = (i-1.0) / (i+1.0) * var_beta_est + 1.0 / ((i+1.0)*(i+1.0))*(diff_res)*(diff_res);
  }
  
  std::vector<double> out;
  out.push_back(beta_est);
  out.push_back(var_beta_est);
  
  return out;
}

// [[Rcpp::export]]
std::vector<std::vector<double> > fit_beta_r_ests_inactive(int n, const std::vector<IntegerVector>& new_inactive_r, 
                                        const NumericMatrix& X, const NumericVector& f_X)
{
  std::vector<std::vector<double> > betas;
  int inactive_size = new_inactive_r.size();
  
  for (int i = 0; i < inactive_size; i++)
  {
    betas.push_back(beta_r_ests(n, new_inactive_r[i], X, f_X));
  }
  
  return betas;
}

// [[Rcpp::export]]
std::vector<IntegerVector> MAQRNNS(int n, int B_0, int B_1, int B_inf, 
                                  const NumericMatrix& X, const NumericVector& f_X, double alpha)
{
  int d = X.ncol();
  double norm_quant = R::qnorm(1-alpha/2,0,1,true,false);
  
  std::vector<IntegerVector> activeSet;
  std::vector<std::vector<double> > betas_active;
  
  // Stage 0: Fit the intercept.
  IntegerVector intercept(d,0);
  activeSet.push_back(intercept);
  betas_active.push_back(beta_r_ests(n, intercept, X, f_X));
  
  for(int i = 0; i < B_0; i++)
  {
    std::vector<IntegerVector> new_activeSet;
    std::vector<std::vector<double> > new_betas_active;
    
    int active_size = activeSet.size();
    
    for (int j = 0; j < active_size; j++)
    {
      if (sum(activeSet[j] > 0) == i)
      {
        new_activeSet = Pot_inter_poly_fcts(activeSet[j], B_1, B_inf);
        new_activeSet = SetDiff(d, new_activeSet, activeSet);
        new_betas_active = fit_beta_r_ests_inactive(n, new_activeSet, X, f_X);
        
        activeSet.insert(activeSet.end(), new_activeSet.begin(), new_activeSet.end());
        betas_active.insert(betas_active.end(), new_betas_active.begin(), new_betas_active.end());
      }
    }
    
    active_size = activeSet.size();
    
    for (int k = active_size-1; k >= 0; k--)
    {
      if (fabs(betas_active[k][0]) < norm_quant*sqrt(betas_active[k][1]))
      {
        activeSet.erase(activeSet.begin() + k);
        betas_active.erase(betas_active.begin()+k);
      }
    }
    
  }
  
  return activeSet;
}

// [[Rcpp::export]]
std::vector<IntegerVector> MAQROS(int n, int B_0, int B_1, int B_inf, 
                                  const NumericMatrix& X, const NumericVector& f_X)
{
  int d = X.ncol();
  
  // TODO: Remove this. Maybe add RQMC instead. n=2^14. m=2^3.
  
  std::vector<IntegerVector> activeSet;
  std::vector<std::vector<double> > betas_active;
  
  // Stage 0: Fit the intercept.
  IntegerVector intercept(d,0);
  activeSet.push_back(intercept);
  betas_active.push_back(beta_r_ests(n, intercept, X, f_X));
  
  for(int i = 0; i < B_0; i++)
  {
    std::vector<IntegerVector> new_activeSet;
    std::vector<std::vector<double> > new_betas_active;
    
    int active_size = activeSet.size();
    
    for (int j = 0; j < active_size; j++)
    {
      if (sum(activeSet[j] > 0) == i)
      {
        new_activeSet = Pot_inter_poly_fcts(activeSet[j], B_1, B_inf);
        new_activeSet = SetDiff(d, new_activeSet, activeSet);
        new_betas_active = fit_beta_r_ests_inactive(n, new_activeSet, X, f_X);
        
        activeSet.insert(activeSet.end(), new_activeSet.begin(), new_activeSet.end());
        betas_active.insert(betas_active.end(), new_betas_active.begin(), new_betas_active.end());
      }
    }
    
    active_size = activeSet.size();
    
    for (int k = active_size-1; k >= 0; k--)
    {
      if (fabs(betas_active[k][0]) < sqrt(2*log(activeSet.size())*betas_active[k][1]))
      {
        activeSet.erase(activeSet.begin() + k);
        betas_active.erase(betas_active.begin()+k);
      }
    }
    
  }
  
  return activeSet;
}

// [[Rcpp::export]]
std::vector<double> beta_r_ests_shrink(int n, const IntegerVector& r, 
                                const NumericMatrix& X, const NumericVector& f_X)
{
  double beta_est=0;
  double var_beta_est = 0;
  
  IntegerVector r_NZ_ind = non_zero_indices(r);
  
  for (int i = 0; i < n; i++)
  {
    double diff_res = f_X[i]*mulidim_basis(r,X(i,_), r_NZ_ind) - beta_est;
    beta_est += (1.0/(i+1))*(diff_res);
    var_beta_est = (i-1.0) / (i+1.0) * var_beta_est + 1.0 / ((i+1.0)*(i+1.0))*(diff_res)*(diff_res);
  }
  
  std::vector<double> out;
  out.push_back(beta_est);
  out.push_back(var_beta_est);
  
  return out;
}

// [[Rcpp::export]]
std::vector<std::vector<double> > fit_beta_r_ests_shrink_inactive(int n, const std::vector<IntegerVector>& new_inactive_r, 
                                                           const NumericMatrix& X, const NumericVector& f_X)
{
  std::vector<std::vector<double> > betas;
  int inactive_size = new_inactive_r.size();
  
  for (int i = 0; i < inactive_size; i++)
  {
    betas.push_back(beta_r_ests_shrink(n, new_inactive_r[i], X, f_X));
  }
  
  return betas;
}

// [[Rcpp::export]]
std::vector<IntegerVector> MAQRSJ(int n, int B_0, int B_1, int B_inf, 
                                 const NumericMatrix& X, const NumericVector& f_X)
{
  int d = X.ncol();
  
  std::vector<IntegerVector> activeSet;
  std::vector<std::vector<double> > betas_active;
  
  // Stage 0: Fit the intercept.
  IntegerVector intercept(d,0);
  activeSet.push_back(intercept);
  betas_active.push_back(beta_r_ests(n, intercept, X, f_X));
  
  for(int i = 0; i < B_0; i++)
  {
    std::vector<IntegerVector> new_activeSet;
    std::vector<std::vector<double> > new_betas_active;
    
    int active_size = activeSet.size();
    
    for (int j = 0; j < active_size; j++)
    {
      if (sum(activeSet[j] > 0) == i)
      {
        new_activeSet = Pot_inter_poly_fcts(activeSet[j], B_1, B_inf);
        new_activeSet = SetDiff(d, new_activeSet, activeSet);
        new_betas_active = fit_beta_r_ests_inactive(n, new_activeSet, X, f_X);
        
        activeSet.insert(activeSet.end(), new_activeSet.begin(), new_activeSet.end());
        betas_active.insert(betas_active.end(), new_betas_active.begin(), new_betas_active.end());
      }
    }
    
    active_size = activeSet.size();
    
    for (int k = active_size-1; k >= 0; k--)
    {
      if (fabs(betas_active[k][0]) < sqrt(2*log(activeSet.size())*betas_active[k][1]))
      {
        activeSet.erase(activeSet.begin() + k);
        betas_active.erase(betas_active.begin()+k);
      }
    }
    
  }
  
  return activeSet;
}


///// Practical Methods //////

// [[Rcpp::export]]
List MARQRNS(int n, int B_0, int B_1, int B_inf, 
                                  const NumericMatrix& X, const NumericVector& f_X,
                                  double alpha = 0.05)
{
  int d = X.ncol();
  double norm_quant = R::qnorm(1-alpha/2,0,1,true,false);
  
  std::vector<IntegerVector> activeSet;
  std::vector<double> betas_active;
  
  // Stage 0: Fit the intercept.
  IntegerVector intercept(d,0);
  activeSet.push_back(intercept);
  betas_active.push_back(beta_r(n, intercept, X, f_X));
  
  for(int i = 0; i < B_0; i++)
  {
    std::vector<IntegerVector> new_activeSet;
    std::vector<double> new_betas_active;
    
    int active_size = activeSet.size();
    
    for (int j = 0; j < active_size; j++)
    {
      if (sum(activeSet[j] > 0) == i)
      {
        new_activeSet = Pot_inter_poly_fcts(activeSet[j], B_1, B_inf);
        new_activeSet = SetDiff(d, new_activeSet, activeSet);
        new_betas_active = fit_beta_r_inactive(n, new_activeSet, X, f_X);
        
        activeSet.insert(activeSet.end(), new_activeSet.begin(), new_activeSet.end());
        betas_active.insert(betas_active.end(), new_betas_active.begin(), new_betas_active.end());
      }
    }
    
    double sigma_est = SEE(n, activeSet, betas_active, X, f_X);
    double noise_bound = norm_quant*sigma_est/sqrt(n);
    
    active_size = activeSet.size();
    
    for (int k = active_size-1; k >= 0; k--)
    {
      if (fabs(betas_active[k]) < noise_bound)
      {
        activeSet.erase(activeSet.begin() + k);
        betas_active.erase(betas_active.begin()+k);
      }
    }
    
  }
  
  List output_vals = List::create(activeSet, betas_active);
  
  return output_vals;
}

// [[Rcpp::export]]
List MARQRSA(int n, int B_0, int B_1, int B_inf, 
                                  const NumericMatrix& X, const NumericVector& f_X,
                                  double alpha = 0.05)
{
  int d = X.ncol();
  
  std::vector<IntegerVector> activeSet;
  std::vector<double> betas_active;
  
  // Stage 0: Fit the intercept.
  IntegerVector intercept(d,0);
  activeSet.push_back(intercept);
  betas_active.push_back(beta_r(n, intercept, X, f_X));
  
  for(int i = 0; i < B_0; i++)
  {
    std::vector<IntegerVector> new_activeSet;
    std::vector<double> new_betas_active;
    
    int active_size = activeSet.size();
    
    for (int j = 0; j < active_size; j++)
    {
      if (sum(activeSet[j] > 0) == i)
      {
        new_activeSet = Pot_inter_poly_fcts(activeSet[j], B_1, B_inf);
        new_activeSet = SetDiff(d, new_activeSet, activeSet);
        new_betas_active = fit_beta_r_inactive(n, new_activeSet, X, f_X);
        
        activeSet.insert(activeSet.end(), new_activeSet.begin(), new_activeSet.end());
        betas_active.insert(betas_active.end(), new_betas_active.begin(), new_betas_active.end());
      }
    }
    
    std::vector<IntegerVector> imp_activeSet;
    imp_activeSet.push_back(intercept);
    std::vector<double> imp_betas_active;
    imp_betas_active.push_back(betas_active[0]);
    std::vector<double> squ_betas;
    double total_variability = 0;
    double curr_variability = 0;
    
    active_size = activeSet.size();
    
    for (int k = 0; k < active_size; k++)
    {
      if (k == 0)
      {
        squ_betas.push_back(-betas_active[k]*betas_active[k]);
      }
      else
      {
        total_variability += betas_active[k]*betas_active[k];
        squ_betas.push_back(betas_active[k]*betas_active[k]);
      }
    }
    
    for (int k = 0; k < active_size-1; k++)
    {
      int max_index = std::distance(squ_betas.begin(),std::max_element(squ_betas.begin(), squ_betas.end()));
      curr_variability += squ_betas[max_index];
      squ_betas[max_index] = -squ_betas[max_index];
      imp_betas_active.push_back(betas_active[max_index]);
      imp_activeSet.push_back(activeSet[max_index]);
      
      if (curr_variability / total_variability > 1- alpha)
      {
        break;
      }
    }
    activeSet = imp_activeSet;
    betas_active = imp_betas_active;
    
  }
  
  List output_vals = List::create(activeSet, betas_active);
  
  return output_vals;
}


// [[Rcpp::export]]
List MARQRLS(int n, int B_0, int B_1, int B_inf, 
                                  const NumericMatrix& X, const NumericVector& f_X)
{
  int d = X.ncol();
  
  std::vector<IntegerVector> activeSet;
  std::vector<double> betas_active;
  
  // Stage 0: Fit the intercept.
  IntegerVector intercept(d,0);
  activeSet.push_back(intercept);
  betas_active.push_back(beta_r(n, intercept, X, f_X));
  
  for(int i = 0; i < B_0; i++)
  {
    std::vector<IntegerVector> new_activeSet;
    std::vector<double> new_betas_active;
    
    int active_size = activeSet.size();
    
    for (int j = 0; j < active_size; j++)
    {
      if (sum(activeSet[j] > 0) == i)
      {
        new_activeSet = Pot_inter_poly_fcts(activeSet[j], B_1, B_inf);
        new_activeSet = SetDiff(d, new_activeSet, activeSet);
        new_betas_active = fit_beta_r_inactive(n, new_activeSet, X, f_X);
        
        activeSet.insert(activeSet.end(), new_activeSet.begin(), new_activeSet.end());
        betas_active.insert(betas_active.end(), new_betas_active.begin(), new_betas_active.end());
      }
    }
    
    double sigma_est = SEE(n, activeSet, betas_active, X, f_X);
    double shrinkage_bound = sigma_est/sqrt(n) * sqrt(2*log(activeSet.size()));
    
    Rcpp::Rcerr << shrinkage_bound << "\n";
    
    active_size = activeSet.size();
    
    for (int k = active_size-1; k >= 0; k--)
    {
      if (fabs(betas_active[k]) < shrinkage_bound)
      {
        activeSet.erase(activeSet.begin() + k);
        betas_active.erase(betas_active.begin()+k);
      }
    }
    
  }
  
  List output_vals = List::create(activeSet, betas_active);
  
  return output_vals;
}

// [[Rcpp::export]]
List MARQRNNS(int n, int B_0, int B_1, int B_inf, 
             const NumericMatrix& X, const NumericVector& f_X, double alpha)
{
  int d = X.ncol();
  
  double norm_quant = R::qnorm(1-alpha/2,0,1,true,false);
  
  std::vector<IntegerVector> activeSet;
  std::vector<std::vector<double>> betas_active;
  
  // Stage 0: Fit the intercept.
  IntegerVector intercept(d,0);
  activeSet.push_back(intercept);
  betas_active.push_back(beta_r_ests(n, intercept, X, f_X));
  
  for(int i = 0; i < B_0; i++)
  {
    std::vector<IntegerVector> new_activeSet;
    std::vector<std::vector<double>> new_betas_active;
    
    int active_size = activeSet.size();
    
    for (int j = 0; j < active_size; j++)
    {
      if (sum(activeSet[j] > 0) == i)
      {
        new_activeSet = Pot_inter_poly_fcts(activeSet[j], B_1, B_inf);
        new_activeSet = SetDiff(d, new_activeSet, activeSet);
        new_betas_active = fit_beta_r_ests_inactive(n, new_activeSet, X, f_X);
        
        activeSet.insert(activeSet.end(), new_activeSet.begin(), new_activeSet.end());
        betas_active.insert(betas_active.end(), new_betas_active.begin(), new_betas_active.end());
      }
    }
    
    active_size = activeSet.size();
    
    for (int k = active_size-1; k >= 0; k--)
    {
      if (fabs(betas_active[k][0]) < norm_quant*sqrt(betas_active[k][1]))
      {
        activeSet.erase(activeSet.begin() + k);
        betas_active.erase(betas_active.begin()+k);
      }
    }
    
  }
  
  List output_vals = List::create(activeSet, betas_active);
  
  return output_vals;
}

// [[Rcpp::export]]
List MARQROS(int n, int B_0, int B_1, int B_inf, 
             const NumericMatrix& X, const NumericVector& f_X)
{
  int d = X.ncol();
  
  double var_factor = 1;
  
  std::vector<IntegerVector> activeSet;
  std::vector<std::vector<double>> betas_active;
  
  // Stage 0: Fit the intercept.
  IntegerVector intercept(d,0);
  activeSet.push_back(intercept);
  betas_active.push_back(beta_r_ests(n, intercept, X, f_X));
  
  for(int i = 0; i < B_0; i++)
  {
    std::vector<IntegerVector> new_activeSet;
    std::vector<std::vector<double>> new_betas_active;
    
    int active_size = activeSet.size();
    
    for (int j = 0; j < active_size; j++)
    {
      if (sum(activeSet[j] > 0) == i)
      {
        new_activeSet = Pot_inter_poly_fcts(activeSet[j], B_1, B_inf);
        new_activeSet = SetDiff(d, new_activeSet, activeSet);
        new_betas_active = fit_beta_r_ests_inactive(n, new_activeSet, X, f_X);
        
        activeSet.insert(activeSet.end(), new_activeSet.begin(), new_activeSet.end());
        betas_active.insert(betas_active.end(), new_betas_active.begin(), new_betas_active.end());
      }
    }
    
    active_size = activeSet.size();
    
    for (int k = active_size-1; k >= 0; k--)
    {
      if (fabs(betas_active[k][0]) < sqrt(2*log(active_size)*betas_active[k][1] / var_factor))
      {
        activeSet.erase(activeSet.begin() + k);
        betas_active.erase(betas_active.begin()+k);
      }
    }
    
  }
  
  List output_vals = List::create(activeSet, betas_active);
  
  return output_vals;
}

///// HAAR ////

// [[Rcpp::export]]
double SEE_haar(int n, const std::vector<IntegerVector>& activeSet,
           const std::vector<double>& betas_active,
           const NumericMatrix& X, const NumericVector& f_X)
{
  double resid_error = 0;
  
  for (int i = 0; i < n; i++)
  {
    double active_res = 0;
    int active_size = activeSet.size();
    
    for (int j = 0; j < active_size; j++)
    {
      active_res += betas_active[j] * mulidimensional_haar_basis(activeSet[j], X(i,_));
    }
    
    double res_diff = (f_X[i] - active_res)*(f_X[i] - active_res) - resid_error;
    
    resid_error += (1.0/(i+1))*res_diff;
  }
  
  return sqrt(resid_error * ((double)n / (n - activeSet.size())));
}

// [[Rcpp::export]]
std::vector<IntegerVector> MAQRNS_haar(int n, int B_0, int B_1, int B_inf, 
                                  const NumericMatrix& X, const NumericVector& f_X,
                                  double alpha = 0.05)
{
  int d = X.ncol();
  double norm_quant = R::qnorm(1-alpha/2,0,1,true,false);
  
  std::vector<IntegerVector> activeSet;
  std::vector<double> betas_active;
  
  // Stage 0: Fit the intercept.
  IntegerVector intercept(d,0);
  activeSet.push_back(intercept);
  betas_active.push_back(haar_beta_r(n, intercept, X, f_X));
  
  for(int i = 0; i < B_0; i++)
  {
    std::vector<IntegerVector> new_activeSet;
    std::vector<double> new_betas_active;
    
    int active_size = activeSet.size();
    
    for (int j = 0; j < active_size; j++)
    {
      if (sum(activeSet[j] > 0) == i)
      {
        new_activeSet = Pot_inter_haar_fcts(activeSet[j], B_1, B_inf);
        new_activeSet = SetDiff(d, new_activeSet, activeSet);
        new_betas_active = haar_beta_r_inactive(n, new_activeSet, X, f_X);
        
        activeSet.insert(activeSet.end(), new_activeSet.begin(), new_activeSet.end());
        betas_active.insert(betas_active.end(), new_betas_active.begin(), new_betas_active.end());
      }
    }
    
    double sigma_est = SEE_haar(n, activeSet, betas_active, X, f_X);
    double noise_bound = norm_quant*sigma_est/sqrt(n);
    
    Rcpp::Rcerr << noise_bound << "\n";
    
    active_size = activeSet.size();
    
    for (int i = active_size-1; i >= 0; i--)
    {
      if (fabs(betas_active[i]) < noise_bound)
      {
        activeSet.erase(activeSet.begin() + i);
        betas_active.erase(betas_active.begin()+i);
      }
    }
    
  }
  
  return activeSet;
}

// [[Rcpp::export]]
std::vector<IntegerVector> MAQRLS_haar(int n, int B_0, int B_1, int B_inf, 
                                  const NumericMatrix& X, const NumericVector& f_X)
{
  int d = X.ncol();
  
  std::vector<IntegerVector> activeSet;
  std::vector<double> betas_active;
  
  // Stage 0: Fit the intercept.
  IntegerVector intercept(d,0);
  activeSet.push_back(intercept);
  betas_active.push_back(haar_beta_r(n, intercept, X, f_X));
  
  for(int i = 0; i < B_0; i++)
  {
    std::vector<IntegerVector> new_activeSet;
    std::vector<double> new_betas_active;
    
    int active_size = activeSet.size();
    
    for (int j = 0; j < active_size; j++)
    {
      if (sum(activeSet[j] > 0) == i)
      {
        new_activeSet = Pot_inter_haar_fcts(activeSet[j], B_1, B_inf);
        new_activeSet = SetDiff(d, new_activeSet, activeSet);
        new_betas_active = haar_beta_r_inactive(n, new_activeSet, X, f_X);
        
        activeSet.insert(activeSet.end(), new_activeSet.begin(), new_activeSet.end());
        betas_active.insert(betas_active.end(), new_betas_active.begin(), new_betas_active.end());
      }
    }
    
    double sigma_est = SEE_haar(n, activeSet, betas_active, X, f_X);
    double shrinkage_bound = sigma_est/sqrt(n) * sqrt(2*log(activeSet.size()));
    //double shrinkage_bound = sigma_est/sqrt(n) * sqrt(2*log(n));
    
    Rcpp::Rcerr << shrinkage_bound << "\n";
    Rcpp::Rcerr << activeSet.size() << "\n";
    
    active_size = activeSet.size();
    
    for (int i = active_size-1; i >= 0; i--)
    {
      if (fabs(betas_active[i]) < shrinkage_bound)
      {
        activeSet.erase(activeSet.begin() + i);
        betas_active.erase(betas_active.begin()+i);
      }
    }
    Rcpp::Rcerr << activeSet.size() << "\n";
  }
  
  return activeSet;
}

// [[Rcpp::export]]
std::vector<NumericVector> QR_email_haar_fast(int n, const IntegerMatrix& U, 
                                         const NumericMatrix& X, const NumericVector& f_X)
{
  std::vector<NumericVector> outList;
  
  int numRows = U.nrow();
  int d = U.ncol();
  
  NumericVector betas(numRows);
  NumericVector Var_betas(numRows);
  
  IntegerMatrix non_zero_terms(numRows, d);
  
  for (int j = 1; j < numRows; j++)
  {
    non_zero_terms(j,_) = non_zero_indices(U(j,_));
  }
  
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < numRows; j++)
    {
      double psi_out = 1;
      
      if (j != 0)
      {
        psi_out = mulidim_haar_basis(U(j,_),X(i,_),non_zero_terms(j,_));
      }
      
      double diff_res = f_X[i]*psi_out - betas[j];
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
  
  outList.push_back(Rcpp::clone(betas));
  outList.push_back(Rcpp::clone(Var_betas));
  
  return outList;
}

// [[Rcpp::export]]
std::vector<NumericVector> QR_email_haar(int n, const IntegerMatrix& U, 
                                    const NumericMatrix& X, const NumericVector& f_X)
{
  std::vector<NumericVector> outList;
  
  int numRows = U.nrow();
  
  NumericVector betas(numRows);
  NumericVector Var_betas(numRows);
  
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < numRows; j++)
    {
      double diff_res = f_X[i]*mulidimensional_haar_basis(U(j,_),X(i,_)) - betas[j];
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
  
  outList.push_back(Rcpp::clone(betas));
  outList.push_back(Rcpp::clone(Var_betas));
  
  return outList;
}

// [[Rcpp::export]]
std::vector<NumericVector> QRS_email_haar(int n, const IntegerMatrix& U,
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
      psi_basis[j] = mulidimensional_haar_basis(U(j,_),X(i,_));
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
  
  outList.push_back(Rcpp::clone(betas));
  outList.push_back(Rcpp::clone(Var_betas));
  outList.push_back(Rcpp::clone(gammas));
  
  return outList;
}


// [[Rcpp::export]]
NumericVector LOF_EST_haar_cpp(int d, const NumericMatrix& X, NumericVector f_X, double var_f, IntegerMatrix U, NumericVector betas, NumericVector gamma, int n) {
  NumericVector Int_sqr_err(n);
  
  int numRow = U.nrow();
  
  for (int i = 0; i < n; i++) {
    double fi_X_i = 0.0;
    for (int j = 0; j < numRow; j++) {
      fi_X_i += gamma[j] * betas[j] * mulidimensional_haar_basis(U(j,_),X(i,_));
    }
    Int_sqr_err[i] = pow(f_X[i] - fi_X_i, 2);
  }
  
  return Int_sqr_err;
}

// [[Rcpp::export]]
std::vector<NumericVector> QR_haar(int n, const IntegerMatrix& U, 
                              const NumericMatrix& X, const NumericVector& f_X,
                              const double& var_f)
{
  int numRows = U.nrow();
  
  NumericVector betas(numRows);
  NumericVector LOF_fit(n);
  
  NumericVector squared_errors(n);
  
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < numRows; j++)
    {
      double diff_res = f_X[i]*mulidimensional_haar_basis(U(j,_),X(i,_)) - betas[j];
      betas[j] += (1.0/(i+1))*(diff_res);
    }
    
    double fi_X_ipi = 0;
    
    for (int j = 0; j < numRows; j++)
    {
      fi_X_ipi += betas[j]*mulidimensional_haar_basis(U(j,_),X(i+1,_));
    }
    
    squared_errors[i] = pow(f_X[i+1] - fi_X_ipi, 2);
    int m = trunc(pow(i, 2.0/3));
    
    if (i > 0)
    {
      LOF_fit[i] = mean(squared_errors[Rcpp::Range(i-m+1, i)]) / var_f;
    }
  }
  
  std::vector<NumericVector> outList;
  outList.push_back(betas);
  outList.push_back(LOF_fit);
  
  return outList;
}

/*** R
*/

/*  
    This file's goal is to figure out why bcgd is returning NAs with certain parameter inputs
    All functions have been copied over from drmdelLasso.c
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include "basisFuncs.h"
#include "utilities.h"


void lp_val(unsigned long m, unsigned long d, double * restrict h, /*inputs*/
            double *restrict* restrict par_mat, /*inputs*/
            double *restrict lp /*outputs*/)
  /* lp (linear predictor) is a vector of length m for the m non-baseline
   *   samples. */
{
  
  /* loop indices */
  unsigned long i, j;
  
  /* calculating lp[i] */
  for (i = 0; i < m; ++i) {
    lp[i] = par_mat[i][0];
    for (j = 1; j < d + 1; ++j) {
      lp[i] = lp[i] + par_mat[i][j] * h[j-1];
    }
  }
  
}


double logDualLGLCache(unsigned long n_total, /*inputs*/
                       unsigned long * restrict n_samples, /*inputs*/
                       unsigned long m, unsigned long d, /*inputs*/
                       double *restrict* restrict par_mat, /*inputs*/
                       double *restrict* restrict* restrict x_mat_H /*inputs*/,
                       double lambda, double * restrict pen_g /*inputs*/)
                       /* Calculating log dual empirical likelihood (+ n \log n) at a given parameter value.
                        * This function assumes the values of the basis function have been cached.
                        * Inputs:
                        *   n_total -- total sample size;
                        *   n_samples -- a vector of length m+1 specifying the size of each sample;
                        *   m -- number of samples - 1;
                        *   d -- dimension of h(x); 
                        *   par_mat -- values of parameters (pointer matrix of dimension m by (d+1));
                        *   x_matH -- 3-D pointer array of data, the first two dimensions are organized as x_0, x_1, ..., x_m,
                        *    with the third dimension being the cached value of the basis function at x_ij
                        *   lambda -- penalty threshold for the group lasso penalty
                        *    pen_g -- a pointer vector of length d for the group specific penalty
                        * Outputs:
                        *   ldlGL_val -- value of group lasso objective function at a given "par" value.
                        */
{
  /* loop indices */
  unsigned long i, j, k;
  
  double * restrict lp;
  lp = (double * restrict) malloc((size_t) (m*sizeof(double)));
  if (lp == NULL) errMsg("malloc() allocation failure for lp!");
  for (i = 0; i < m; ++i) {
    lp[i] = 0;
  }
  
  /* Create vector for l2 norms in penalty*/
  double * restrict l2;
  l2 = (double * restrict) malloc((size_t) (d*sizeof(double)));
  if (l2 == NULL) errMsg("malloc() allocation failure for l2!");
  for (i = 0; i < d; ++i) {
    l2[i] = 0;
  }
  
  double * restrict rho;
  rho = (double * restrict) malloc((size_t) ((m+1)*sizeof(double)));
  if (rho == NULL) errMsg("malloc() allocation failure for rho!");
  for (i = 0; i < m+1; ++i) {
    rho[i] = (double)n_samples[i]/(double)n_total;
  }
  
  
  /*other variables*/
  double S;
  
  /*define output*/
  double ldlGL_val;
  ldlGL_val = 0.0;
  
  for (i = 0; i < m+1; ++i) {
    for (j = 0; j < n_samples[i]; ++j) {
      
      lp_val(m, d, x_mat_H[i][j], par_mat, lp); /*update lp*/
  
  /* calculating q_i */
  S = rho[0];
  for (k = 0; k < m; ++k) {
    S += rho[k+1] * exp(lp[k]);
  }
  
  if (i == 0) {
    ldlGL_val = ldlGL_val - log(S);  
  } else {
    ldlGL_val = ldlGL_val + lp[i-1] - log(S);
  }
  
    }
    if (i==0){
      continue;
    }
    pen_val(d, par_mat[i-1], l2);/*update l2*/
  }
  
  /*Add pen to ldlGL_val*/
  for (i=0; i<d; ++i){
    ldlGL_val -= lambda*pen_g[i]*sqrt(l2[i]);/*subtracting because we will return the negative*/
  }
  
  /* free arrays */
  free((void *) lp);
  free((void *) rho);
  free((void *) l2);
  
  return ldlGL_val;
  
}

void R_val(unsigned long m, unsigned long d, double * restrict h, /*inputs*/
  double *restrict* restrict par_mat, double * restrict n_samples, /*inputs*/
  double *restrict R /*outputs*/)
  /* R is a vector of length m for the m non-baseline samples.
   * R[i] = n_samples[i] * exp(alpha_i + beta_i^T h(x)), i = 1, 2, ..., m*/
{
  
  /* loop indices */
  unsigned long i, j;
  
  /* calculating r_i */
  for (i = 0; i < m; ++i) {
    R[i] = par_mat[i][0];
    for (j = 1; j < d + 1; ++j) {
      R[i] = R[i] + par_mat[i][j] * h[j-1];
    }
    R[i] = n_samples[i+1] * exp(R[i]);
  }
  
}

void grad_Gt(unsigned long n_total, /*inputs*/
  unsigned long * restrict n_samples, /*inputs*/
  unsigned long m, unsigned long d, /*inputs*/
  double *restrict* restrict par_mat, /*inputs*/
  double *restrict* restrict* restrict x_mat_H /*inputs*/,
  unsigned long g, /*inputs*/
  double * restrict grad_G /*output*/)
  /* Calculating the gradient for a particular group of entries
   * parameter value.
   * Inputs:
   *   n_total -- total sample size;
   *   n_samples -- a vector of length m+1 specifying the size of each sample;
   *   m -- number of samples - 1;
   *   d -- dimension of h(x); 
   *   par_mat -- values of parameters (pointer matrix of dimension m by (d+1));
   *   x_matH -- 3-D pointer array of data, the first two dimensions are organized as x_0, x_1, ..., x_m,
   *    with the third dimension being the cached value of the basis function at x_ij
   *   g -- the lasso group
   * Outputs:
   *   grad_G -- a pointer array of dimension m; value of the gradient of the ldl (log dual empirical likelihood) at a given "par" value and at the group g indices.
   */
{
  /* loop indices */
  unsigned long i, j, k;
  
  double * restrict R;
  R = (double * restrict) malloc((size_t) (m*sizeof(double)));
  if (R == NULL) errMsg("malloc() allocation failure for R!");
  for (i = 0; i < m; ++i) {
    R[i] = 0.0;
  }
  
  double * restrict rho;
  rho = (double * restrict) malloc((size_t) ((m+1)*sizeof(double)));
  if (rho == NULL) errMsg("malloc() allocation failure for rho!");
  for (i = 0; i < m+1; ++i) {
    rho[i] = (double)n_samples[i]/(double)n_total;
  }
  
  /*other variables*/
  double S;
  double tmp_double;
  
  /*initializing grad_G as safeguard*/
  for (i = 0; i < m; ++i) {
    grad_G[i] = 0.0;
  }
  
  for (i = 0; i < m+1; ++i) {
    
    for (j = 0; j < n_samples[i]; ++j) {
      
      /*update H = (1, h^T)^T*/
      //(*h_func)(x_mat[i][j], H+1); /*update H*/
      
      R_val(m, d, x_mat_H[i][j], par_mat, rho, R); /*update R*/
      
      /*calculating S*/
      S = rho[0];
      for (k = 0; k < m; ++k) {
        S += R[k];
      }
      
      /* calculating the gradient of ldl */
      for (k = 0; k < m; ++k) {
        
        tmp_double = -R[k]/S;
        
        if(g==0){ // group is the constant group, so factor is 1
          grad_G[k] += tmp_double;
        } else{
          grad_G[k] += tmp_double * x_mat_H[i][j][g-1]; //use g-1 as index since index starts at first term of basis function
        }
        
      }
      
      if (i > 0) {
        if(g==0){
          grad_G[i-1] = 1 + grad_G[i-1];
        } else{
          grad_G[i-1] = x_mat_H[i][j][g-1] + grad_G[i-1];
        }
      }
      
    }
    
  }
  
  /* free arrays */
  free((void *) R);
  free((void *) rho);
  
}

double hG(unsigned long n_total, /*inputs*/
  double * restrict n_samples, /*inputs*/
  unsigned long m, unsigned long d, /*inputs*/
  double *restrict* restrict par_mat, /*inputs*/
  double *restrict* restrict* restrict x_mat_H /*inputs*/,
  unsigned long g /*inputs*/)
  /* Calculate the maximum diagonal of the negative Hessian for a given group.
   * Inputs:
   *   n_total -- length of data;
   *   n_samples -- a (double, not unsigned long!) vector of length m+1
   *     specifying the size of each sample;
   *   m -- number of samples - 1;
   *   d -- dimension of h(x); 
   *   par_mat -- values of parameters (pointer matrix of dimension m by (d+1));
   *   x_matH -- 3-D pointer array of data, the first two dimensions are organized as x_0, x_1, ..., x_m,
   *    with the third dimension being the cached value of the basis function at x_ij
   *   g -- the lasso group
   * Outputs:
   *   h_g -- the maximum of the diagonal along the group's entries.
   */
{
  /* loop indices */
  unsigned long i, j, k;
  
  /* Create pointer array of length m for the diagonal of the group entries*/
  /* initialize as 0 */
  double * restrict HDiagG;
  HDiagG = (double * restrict) malloc((size_t) (m*sizeof(double)));
  if (HDiagG == NULL) errMsg("malloc() allocation failure for HDiagG!");
  for (i = 0; i < m; ++i) {
    HDiagG[i] = 0.0;
  }
  
  double * restrict R;
  R = (double * restrict) malloc((size_t) (m*sizeof(double)));
  if (R == NULL) errMsg("malloc() allocation failure for R!");
  for (i = 0; i < m; ++i) {
    R[i] = 0.0;
  }
  
  /*qaa matrix*/
  double *restrict* restrict qaa;
  qaa = (double *restrict* restrict) malloc((size_t) (m*sizeof(double*)));
  if (qaa == NULL) errMsg("malloc() allocation failure for qaa!");
  
  qaa[0] = (double * restrict) malloc((size_t) ((m*m)*sizeof(double)));
  if (qaa[0] == NULL) errMsg("malloc() allocation failure for qaa[0]!");
  for(i = 1; i < m; ++i) {
    qaa[i] = qaa[i-1] + m;
  }
  
  for (i = 0; i < m; ++i) {
    for (k = 0; k < m; ++k) {
      qaa[i][k] = 0.0;
    }
  }
  
  /*other variables*/
  double S;
  
  /*output variable*/
  double h_g;
  
  for (i = 0; i < m+1; ++i) {
    for (j = 0; j < n_samples[i]; ++j) {
      
      R_val(m, d, x_mat_H[i][j], par_mat, n_samples, R); /*update R*/
  
  /*calculating S*/
  S = n_samples[0];
  for (k = 0; k < m; ++k) {
    S += R[k];
  }
  
  for (k = 0; k < m; ++k) {
    qaa[k][k] = R[k] * R[k] / (S * S) - R[k] / S;
  }
  
  // Add relevant values to diagonal vector
  for (k = 0; k < m; ++k) {
    if(g==0){ // g represents with the constants, and we should simply multiply qaa[k][k] by 1
      HDiagG[k] += qaa[k][k]*1;
    } else{ // g represents a basis function term, so we should use the square of x_mat_H[i][j][g-1]
      // g-1 corrects for the fact that the constant is not represented in this matrix, so the 0th group would be the first basis function term
      HDiagG[k] += qaa[k][k]*x_mat_H[i][j][g-1]*x_mat_H[i][j][g-1];
    }
  }
    }
  }
  
  // Initialize output as the negative of the group's first diagonal entry
  h_g = (-1)*HDiagG[0];
  
  /*calculating the group's diagonal's of the negative hessian, and saving the maximum.*/
  for (k = 1; k < m; ++k) {
    if(h_g < (-1)*HDiagG[k]){
      h_g = (-1)*HDiagG[k];
    }
  }
  
  /* bound h_g from below to ensure convergence */
  if (h_g < 0.01){
    h_g = 0.01;
  } else if(h_g > 1000000000.0){ // bound from above as well
    h_g = 1000000000.0;
  }
  
  /* free arrays */
  free((void *) R);
  free((void *) HDiagG);
  
  free((void *) qaa[0]);
  free((void *) qaa);
  
  /*return h_g*/
  return h_g;
}


void bcgd(
    double *restrict n_total, /* total pooled sample size */
    double *restrict n_samples, /* vector samples sizes*/
    double *restrict m, double *restrict d, /* number of samples and length of basis function*/
    double *restrict model, double *restrict x, /*vector of sample values*/
    double *restrict lambda, double *restrict pen_g, /* penalty values*/
    double *restrict omega_0, double *restrict psi, double *restrict sigma, /* optimization hyperparameters*/
    double *restrict threshold, double *restrict max_iters, /*convergence criteria*/
    double *restrict theta_0, /*input/output: initial -> optimized value of parameter vector*/
    double *restrict opt_val, /*output: the minimized value of the negLDLGL function */
    double *restrict total_iters /*output: number of iterations to convergence*/){
  //Create loop indices
  unsigned long i, j;
  
  //Create vector of sample sizes
  unsigned long * restrict n_samples_use;
  n_samples_use = (unsigned long * restrict) malloc((size_t) (((unsigned long)*m + 1)*sizeof(unsigned long)));                                                             
  if (n_samples_use == NULL) {
    errMsg("malloc() allocation failure for m_samples_use!");
  }
  for (i = 0; i < ((unsigned long)*m + 1); ++i) {
    n_samples_use[i] = (unsigned long)n_samples[i];
  }
  
  // Setup X as matrix for computations
  double *restrict* restrict x_mat;
  x_mat = (double *restrict* restrict) malloc((size_t) (((unsigned
                                                            long)*m+1)*sizeof(double*)));
  if (x_mat == NULL) errMsg("malloc() allocation failure for x_mat!");
  x_mat[0] = x;
  for (i = 1; i < ((unsigned long)*m + 1); ++i){
    x_mat[i] = x_mat[i-1] + n_samples_use[i - 1];
  }
  
  // Assign initial values to parameter matrix
  double *restrict* restrict par_mat;
  par_mat = (double *restrict* restrict) malloc((size_t)
                                                  (((unsigned long)*m)*sizeof(double*)));
  if (par_mat == NULL) errMsg("malloc() allocation failure for par_mat!");
  par_mat[0] = theta_0;
  for (i = 1; i < (unsigned long)*m; ++i){
    par_mat[i] = par_mat[i-1] + ((unsigned long)*d + 1);
  }
  
  // Assign the proper basis function
  void (*h_func)(double, double *restrict) = 0; // define to avoid warning
  
  switch ((unsigned long)*model)
  {
  case 1 :
    if ((unsigned long)*d != 1) {
      errMsg("For model 1, h(x) = x, d must be 1!");
    }
    h_func = &h1x;
    break;
    
  case 2 :
    if ((unsigned long)*d != 1) {
      errMsg("For model 2, h(x) = log(x), d must be 1!");
    }
    h_func = &h1logx;
    break;
    
  case 3 :
    if ((unsigned long)*d != 1) {
      errMsg("For model 3, h(x) = sqrt(x), d must be 1!");
    }
    h_func = &h1sqrtx;
    break;
    
  case 4 :
    if ((unsigned long)*d != 1) {
      errMsg("For model 4, h(x) = x^2, d must be 1!");
    }
    h_func = &h1xSquare;
    break;
    
  case 5 :
    if ((unsigned long)*d != 2) {
      errMsg("For model 5 (Normal model), h(x) = (x, x^2), d must be 2!");
    }
    h_func = &h2Normal;
    break;
    
  case 6 :
    if ((unsigned long)*d != 2) {
      errMsg("For model 6 (Gamma model), h(x) = (x, log(x)), d must be 2!");
    }
    h_func = &h2Gamma;
    break;
    
  case 7 :
    if ((unsigned long)*d != 3) {
      errMsg("For model 7, h(x) = (log(x), sqrt(x), x), d must be 3!");
    }
    h_func = &h3a;
    break;
    
  case 8 :
    if ((unsigned long)*d != 3) {
      errMsg("For model 8, h(x) = (log(x), sqrt(x), x^2), d must be 3!");
    }
    h_func = &h3b;
    break;
    
  case 9 :
    if ((unsigned long)*d != 3) {
      errMsg("For model 9, h(x) = (log(x), x, x^2), d must be 3!");
    }
    h_func = &h3c;
    break;
    
  case 10 :
    if ((unsigned long)*d != 3) {
      errMsg("For model 10, h(x) = (sqrt(x), x, x^2), d must be 3!");
    }
    h_func = &h3d;
    break;
    
  case 11 :
    if ((unsigned long)*d != 4) {
      errMsg("For model 11, h(x) = (log(x), sqrt(x), x, x^2), d must be 4!");
    }
    h_func = &h4a;
    break;
    
  case 12 :
    if ((unsigned long)*d != 5) {
      errMsg("For model 11, h(x) = (log(x), log(x)^2, sqrt(x), x, x^2), d must be 5!");
    }
    h_func = &h5a;
    break;
    
  default :
    errMsg("'Model' must be an integer between 1 and 12 or a function of a single data point");
  break;
  }
  
  // Setup x_mat_H as matrix for computations, which will cache the results of the basis function for each x_ij
  double *restrict* restrict* restrict x_mat_H;
  x_mat_H = (double *restrict* restrict* restrict) malloc((size_t) (((unsigned
                                                                        long)*m+1)*sizeof(double**)));
  if (x_mat_H == NULL) errMsg("malloc() allocation failure for x_mat_H!");
  // define idx to track index of x in flat pointer
  unsigned long idx = 0;
  
  for (i = 0; i < ((unsigned long)*m + 1); ++i){
    //update idx
    if(i > 0){
      idx += n_samples_use[i-1];
    }
    
    x_mat_H[i] = (double *restrict* restrict) malloc((size_t) ((n_samples_use[i])*sizeof(double*)));
    if (x_mat_H[i] == NULL) errMsg("malloc() allocation failure for x_mat_H[i]!");
    for(j = 0; j < n_samples_use[i]; j++){
      x_mat_H[i][j] = (double * restrict) malloc((size_t) (((unsigned long)*d)*sizeof(double)));
      if (x_mat_H[i][j] == NULL) errMsg("malloc() allocation failure for x_mat_H[i][j]!");
      (*h_func)(x[idx+j], x_mat_H[i][j]); // update x_mat_H[i][j]
    }
  }
  
  // Create inner loop index for groups
  unsigned long g;
  
  // Create initial value of function
  double initial_ldlGL;
  
  // Create final value of function
  double final_ldlGL = 0;
  
  // Create h_g
  double h_g;
  
  // Set threshold for omega_t during loop
  double omega_threshold = pow(10.0, -30.0);
  
  // Begin outer loop
  for(i = 1;i<=(unsigned long)*max_iters;i++){
    // set initial value of the function
    if(i == 1){
      initial_ldlGL = (-1)*logDualLGLCache((unsigned long)*n_total, n_samples_use, (unsigned long)*m,
                       (unsigned long)*d, par_mat, x_mat_H, *lambda, pen_g);
    } else {
      initial_ldlGL = final_ldlGL; // can do this here because it will be equal to the most recent evaluation
    }
    
    // Begin inner loop over groups
    for(g=0;g<(unsigned long)*d+1;g++){
      
      // Compute grad_gt
      // Since we cannot return arrays in C, must setup as pointer
      double *restrict grad_G = (double *restrict) malloc((size_t) (((unsigned long)*m)*sizeof(double)));
      if (grad_G == NULL) {
        errMsg("malloc() allocation failure for grad_G!");
      }
      
      grad_Gt((unsigned long)*n_total, n_samples_use, (unsigned long)*m,
              (unsigned long)*d, par_mat, x_mat_H, g, /*inputs*/
              grad_G /*output*/);
      
      // Compute h_g
      h_g = hG((unsigned long)*n_total, n_samples, (unsigned long)*m,
               (unsigned long)*d, par_mat, x_mat_H, g)/*inputs*/;
      
      
      // Initialize d_g, the direction vector
      double *restrict d_g = (double *restrict) malloc((size_t) (((unsigned long)*m)*sizeof(double)));
      if (d_g == NULL) {
        errMsg("malloc() allocation failure for d_g!");
      }
      
      // Set value for delta_t, which is computed conditionally on the g != 0
      double delta_t = 0.0;
      // initialize the initial value of the objective function
      double init_ldlGl_iter;
      
      // Begin cases for d_g computation
      
      if(g==0){ // 0 is the group of intercepts
        // Loop through d_g and set values
        for(j = 0; j<(unsigned long)*m; j++){
          d_g[j] = (1.0/h_g) * grad_G[j];
          // Compute delta that will be used in line search
          delta_t += (-1.0)*d_g[j]*grad_G[j];
        }
        
        // set initial value for objective
        init_ldlGl_iter = initial_ldlGL;
      } else{
        // First we want to compute the norm of grad_g + h_g*theta_g
        double norm = 0.0;
        for(j = 0; j<(unsigned long)*m; j++){
          // compute the vector entry
          double vecEntry = grad_G[j] + h_g*par_mat[j][g];
          norm += vecEntry * vecEntry;
        }
        norm = sqrt(norm);
        if(norm <= *lambda * pen_g[g-1]){ // we set d equal to the negative of theta. g offset since it is only of length d, and g goes from 0 to d 
          for(j = 0; j<(unsigned long)*m; j++){
            d_g[j] = (-1.0)*par_mat[j][g];
          }
        } else{
          for(j = 0; j<(unsigned long)*m; j++){
            d_g[j] = (1.0/h_g) * (grad_G[j] - *lambda * pen_g[g-1] * (grad_G[j] + h_g * par_mat[j][g]) / norm);
          }
        }
        
        // Compute delta for coefficient group
        double delta_norm_1 = 0.0;
        double delta_norm_2 = 0.0;
        
        for(j = 0; j<(unsigned long)*m; j++){
          // first add -d_g*theta_g to delta
          delta_t += (-1.0)*d_g[j]*grad_G[j];
          // update the two norms
          // compute sum of entries first
          double delta_norm_1_entry = par_mat[j][g] + d_g[j];
          delta_norm_1 += delta_norm_1_entry * delta_norm_1_entry;
          delta_norm_2 += par_mat[j][g] * par_mat[j][g];
        }
        
        // Take root of norms
        delta_norm_1 = sqrt(delta_norm_1);
        delta_norm_2 = sqrt(delta_norm_2);
        delta_t = delta_t + *lambda * pen_g[g-1] * (delta_norm_1 - delta_norm_2);
        
        // set initial value for objective
        init_ldlGl_iter = final_ldlGL; // can do this here because it will be equal to the most recent evaluation
      }
      
      // initialize omega_t
      double omega_t = *omega_0;
      
      // before starting line search, copy original values of gth column for efficient updates
      double *restrict original_col_g = (double *restrict) malloc((size_t) (((unsigned long)*m)*sizeof(double)));
      if (original_col_g == NULL) {
        errMsg("malloc() allocation failure for original_col_g!");
      }
      for (j = 0; j < (unsigned long)*m; ++j) {
        original_col_g[j] = par_mat[j][g];
      }
      
      
      // do line search for omega_t
      while (omega_t > omega_threshold){
        // update par_mat
        for(j = 0; j<(unsigned long)*m; j++){
          par_mat[j][g] = original_col_g[j] + omega_t * d_g[j];
        }
        // update final_ldlGL
        final_ldlGL = (-1)*logDualLGLCache((unsigned long)*n_total, n_samples_use, (unsigned long)*m,
                       (unsigned long)*d, par_mat, x_mat_H, *lambda, pen_g);
        // Compute difference
        double sub = final_ldlGL - init_ldlGl_iter;
        // Check if step size is optimal
        if(sub <= delta_t * omega_t * (*sigma)){
          break;
        } else{
          omega_t= omega_t* (*psi);
        }
      }
      
      // Free group loop pointers
      free(grad_G);
      free(d_g);
      free(original_col_g);
    }
    // Check if alg has converged
    if(fabs(final_ldlGL - initial_ldlGL) < *threshold){
      *total_iters = i;
      break;
    }else if(i == *max_iters){
      *total_iters = i;
    }
  }
  // Set the minimum value to return
  // Free allocated memory for x-mat_H
  for (unsigned long i = 0; i < ((unsigned long)*m + 1); ++i) {
    for (unsigned long j = 0; j < n_samples_use[i]; ++j) {
      free(x_mat_H[i][j]);  // Free each double* (the actual data storage)
    }
    free((void *) x_mat_H[i]);  // Free each double** (array of double*)
  }
  free((void *) x_mat_H);  // Free the top-level double*** pointer
  
  *opt_val = final_ldlGL;
  
  // free other variables
  free(n_samples_use);
  free((void *) par_mat);
}

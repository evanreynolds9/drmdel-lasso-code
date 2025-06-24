/* This software is written by Song Cai and published under GPLv3.
 *
 * Version 1.3, April 08, 2014.
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

void r_val(unsigned long m, unsigned long d, double * restrict h, /*inputs*/
    double *restrict* restrict par_mat, /*inputs*/
    double *restrict r /*outputs*/)
/* r is a vector of length m for the m non-baseline samples. */
{

  /* loop indices */
  unsigned long i, j;

  /* calculating r_i */
  for (i = 0; i < m; ++i) {
    r[i] = par_mat[i][0];
    for (j = 1; j < d + 1; ++j) {
      r[i] = r[i] + par_mat[i][j] * h[j-1];
    }
    r[i] = exp(r[i]);
  }

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

void pen_val(unsigned long d, /*inputs*/
  double *restrict par_mat_i, /*inputs*/
  double *restrict l2 /*outputs*/)
/* l2 is a d length vector that stores the squared sum
 * of the coefficients of each basis function term */
{
  unsigned long i;
  for (i=1; i<d+1; ++i){
    l2[i-1] += par_mat_i[i]*par_mat_i[i];
  }
}

double logDualL(unsigned long n_total, /*inputs*/
    unsigned long * restrict n_samples, /*inputs*/
    unsigned long m, unsigned long d, /*inputs*/
    double *restrict* restrict par_mat, /*inputs*/
    void (*h_func)(double, double * restrict), /*input*/
    double *restrict* restrict x_mat /*inputs*/)
/* Calculating log dual empirical likelihood (+ n \log n) at a given parameter value.
 * Inputs:
 *   n_total -- total sample size;
 *   n_samples -- a vector of length m+1 specifying the size of each sample;
 *   m -- number of samples - 1;
 *   d -- dimension of h(x); 
 *   par_mat -- values of parameters (pointer matrix of dimension m by (d+1));
 *   h_func -- the basis function of the DRM;
 *   x_mat -- 2-D pointer array of data, organized as x_0, x_1, ..., x_m.
 * Outputs:
 *   ldl_val -- value of ldl (log dual empirical likelihood) at a given "par" value.
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

  double * restrict h;
  h = (double * restrict) malloc((size_t) (d*sizeof(double)));
  if (h == NULL) errMsg("malloc() allocation failure for h!");
  for (i = 0; i < d; ++i) {
    h[i] = 0;
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
  double ldl_val;
  ldl_val = 0.0;

  for (i = 0; i < m+1; ++i) {

    for (j = 0; j < n_samples[i]; ++j) {

      (*h_func)(x_mat[i][j], h); /*update h*/

      lp_val(m, d, h, par_mat, lp); /*update lp*/

      /* calculating q_i */
      S = rho[0];
      for (k = 0; k < m; ++k) {
        S += rho[k+1] * exp(lp[k]);
      }

      if (i == 0) {
        ldl_val = ldl_val - log(S);  
      } else {
        ldl_val = ldl_val + lp[i-1] - log(S);
      }

    }

  }

  /* free arrays */
  free((void *) lp);
  free((void *) h);
  free((void *) rho);

  return ldl_val;

}

void logDualLWrapper(double * restrict n_total, /*inputs*/
    double * restrict n_samples, /*inputs*/
    double * restrict m, double * restrict d,
    double * restrict par, /*inputs*/
    double * restrict model, double * restrict x, /*inputs*/
    double * restrict ldl_val /*output*/)
/* Calculating log dual empirical likelihood (+ n \log n) at a given parameter value.
 * Inputs:
 *   n_total -- total sample size;
 *   n_samples -- a vector of length m+1 specifying the size of each sample;
 *   m -- number of samples - 1;
 *   d -- dimension of h(x); 
 *   par -- values of parameters (length of m(d+1)), organized as \theta_1, \cdots, \theta_m;
 *   h_func -- the basis function of the DRM;
 *   x -- data, organized as x_0, x_1, ..., x_m.
 * Outputs:
 *   ldl_val -- value of ldl (log dual empirical likelihood) at a given "par" value.
 */
{
  /* loop indices */
  unsigned long i;

  unsigned long * restrict n_samples_use;
  n_samples_use = (unsigned long * restrict) malloc((size_t) (((unsigned
            long)*m + 1)*sizeof(unsigned long)));
  if (n_samples_use == NULL) {
    errMsg("malloc() allocation failure for m_samples_use!");
  }
  for (i = 0; i < ((unsigned long)*m + 1); ++i) {
    n_samples_use[i] = (unsigned long)n_samples[i];
  }

  double *restrict* restrict par_mat;
  double *restrict* restrict x_mat;

  /* converting x and par to matrices */
  par_mat = (double *restrict* restrict) malloc((size_t)
      (((unsigned long)*m)*sizeof(double*)));
  if (par_mat == NULL) errMsg("malloc() allocation failure for par_mat!");
  par_mat[0] = par;
  for (i = 1; i < (unsigned long)*m; ++i){
    par_mat[i] = par_mat[i-1] + ((unsigned long)*d + 1);
  }

  x_mat = (double *restrict* restrict) malloc((size_t) (((unsigned
            long)*m+1)*sizeof(double*)));
  if (x_mat == NULL) errMsg("malloc() allocation failure for x_mat!");
  x_mat[0] = x;
  for (i = 1; i < ((unsigned long)*m + 1); ++i){
    x_mat[i] = x_mat[i-1] + n_samples_use[i - 1];
  }

  /* calculating log dual empirical likelihood value at 'par' */
  switch ((unsigned long)*model)
  {
    case 1 :
      /*printf("h(x) = (x)");*/
      if ((unsigned long)*d != 1) {
        errMsg("For model 1, h(x) = x, d must be 1!");
      }
      *ldl_val = logDualL((unsigned long)*n_total,
          /*n_samples_use, (unsigned long)*m, (unsigned long)*d,*/
          n_samples_use, (unsigned long)*m, 1,
          par_mat, &h1x, x_mat);
      break;

    case 2 :
      /*printf("h(x) = (log(x))");*/
      if ((unsigned long)*d != 1) {
        errMsg("For model 2, h(x) = log(x), d must be 1!");
      }
      *ldl_val = logDualL((unsigned long)*n_total,
          /*n_samples_use, (unsigned long)*m, (unsigned long)*d,*/
          n_samples_use, (unsigned long)*m, 1,
          par_mat, &h1logx, x_mat);
      break;

    case 3 :
      /*printf("h(x) = (sqrt(x))");*/
      if ((unsigned long)*d != 1) {
        errMsg("For model 3, h(x) = sqrt(x), d must be 1!");
      }
      *ldl_val = logDualL((unsigned long)*n_total,
          /*n_samples_use, (unsigned long)*m, (unsigned long)*d,*/
          n_samples_use, (unsigned long)*m, 1,
          par_mat, &h1sqrtx, x_mat);
      break;

    case 4 :
      /*printf("h(x) = (x^2)");*/
      if ((unsigned long)*d != 1) {
        errMsg("For model 4, h(x) = x^2, d must be 1!");
      }
      *ldl_val = logDualL((unsigned long)*n_total,
          /*n_samples_use, (unsigned long)*m, (unsigned long)*d,*/
          n_samples_use, (unsigned long)*m, 1,
          par_mat, &h1xSquare, x_mat);
      break;

    case 5 :
      /*printf("h(x) = (x, x^2) -- Normal model");*/
      if ((unsigned long)*d != 2) {
        errMsg("For model 5 (Normal model), h(x) = (x, x^2), d must be 2!");
      }
      *ldl_val = logDualL((unsigned long)*n_total,
          /*n_samples_use, (unsigned long)*m, (unsigned long)*d,*/
          n_samples_use, (unsigned long)*m, 2,
          par_mat, &h2Normal, x_mat);
      break;

    case 6 :
      /*printf("h(x) = (x, log(x)) -- Gamma model");*/
      if ((unsigned long)*d != 2) {
        errMsg("For model 6 (Gamma model), h(x) = (x, log(x)), d must be 2!");
      }
      *ldl_val = logDualL((unsigned long)*n_total,
          /*n_samples_use, (unsigned long)*m, (unsigned long)*d,*/
          n_samples_use, (unsigned long)*m, 2,
          par_mat, &h2Gamma, x_mat);
      break;

    case 7 :
      /*printf("h(x) = (log(x), sqrt(x), x)");*/
      if ((unsigned long)*d != 3) {
        errMsg("For model 7, h(x) = (log(x), sqrt(x), x), d must be 3!");
      }
      *ldl_val = logDualL((unsigned long)*n_total,
          /*n_samples_use, (unsigned long)*m, (unsigned long)*d,*/
          n_samples_use, (unsigned long)*m, 3,
          par_mat, &h3a, x_mat);
      break;

    case 8 :
      /*printf("h(x) = (log(x), sqrt(x), x^2)");*/
      if ((unsigned long)*d != 3) {
        errMsg("For model 8, h(x) = (log(x), sqrt(x), x^2), d must be 3!");
      }
      *ldl_val = logDualL((unsigned long)*n_total,
          /*n_samples_use, (unsigned long)*m, (unsigned long)*d,*/
          n_samples_use, (unsigned long)*m, 3,
          par_mat, &h3b, x_mat);
      break;

    case 9 :
      /*printf("h(x) = (log(x), x, x^2)");*/
      if ((unsigned long)*d != 3) {
        errMsg("For model 9, h(x) = (log(x), x, x^2), d must be 3!");
      }
      *ldl_val = logDualL((unsigned long)*n_total,
          /*n_samples_use, (unsigned long)*m, (unsigned long)*d,*/
          n_samples_use, (unsigned long)*m, 3,
          par_mat, &h3c, x_mat);
      break;

    case 10 :
      /*printf("h(x) = (sqrt(x), x, x^2)");*/
      if ((unsigned long)*d != 3) {
        errMsg("For model 10, h(x) = (sqrt(x), x, x^2), d must be 3!");
      }
      *ldl_val = logDualL((unsigned long)*n_total,
          /*n_samples_use, (unsigned long)*m, (unsigned long)*d,*/
          n_samples_use, (unsigned long)*m, 3,
          par_mat, &h3d, x_mat);
      break;

    case 11 :
      /*printf("h(x) = (log(x), sqrt(x), x, x^2)");*/
      if ((unsigned long)*d != 4) {
        errMsg("For model 11, h(x) = (log(x), sqrt(x), x, x^2), d must be 4!");
      }
      *ldl_val = logDualL((unsigned long)*n_total,
          /*n_samples_use, (unsigned long)*m, (unsigned long)*d,*/
          n_samples_use, (unsigned long)*m, 4,
          par_mat, &h4a, x_mat);
      break;
    
    case 12 :
      /*printf("h(x) = (log(x), log(x)^2, sqrt(x), x, x^2)");*/
      if ((unsigned long)*d != 5) {
        errMsg("For model 11, h(x) = (log(x), log(x)^2, sqrt(x), x, x^2), d must be 5!");
      }
      *ldl_val = logDualL((unsigned long)*n_total,
                        /*n_samples_use, (unsigned long)*m, (unsigned long)*d,*/
                        n_samples_use, (unsigned long)*m, 5,
                        par_mat, &h5a, x_mat);
    break;

    default :
      errMsg("'Model' must be an integer between 1 and 12 or a function of a single data point");
      break;

  }

  /* free arrays */
  free((void *) n_samples_use);
  free((void *) x_mat);
  free((void *) par_mat);

}

double logDualLGL(unsigned long n_total, /*inputs*/
  unsigned long * restrict n_samples, /*inputs*/
  unsigned long m, unsigned long d, /*inputs*/
  double *restrict* restrict par_mat, /*inputs*/
  void (*h_func)(double, double * restrict), /*input*/
  double *restrict* restrict x_mat /*inputs*/,
  double lambda, double * restrict pen_g /*inputs*/)
  /* Calculating log dual empirical likelihood (+ n \log n) at a given parameter value.
   * Inputs:
   *   n_total -- total sample size;
   *   n_samples -- a vector of length m+1 specifying the size of each sample;
   *   m -- number of samples - 1;
   *   d -- dimension of h(x); 
   *   par_mat -- values of parameters (pointer matrix of dimension m by (d+1));
   *   h_func -- the basis function of the DRM;
   *   x_mat -- 2-D pointer array of data, organized as x_0, x_1, ..., x_m.
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
  
  double * restrict h;
  h = (double * restrict) malloc((size_t) (d*sizeof(double)));
  if (h == NULL) errMsg("malloc() allocation failure for h!");
  /* Create vector for l2 norms in penalty*/
  double * restrict l2;
  l2 = (double * restrict) malloc((size_t) (d*sizeof(double)));
  if (l2 == NULL) errMsg("malloc() allocation failure for l2!");
  for (i = 0; i < d; ++i) {
    h[i] = 0;
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
      
      (*h_func)(x_mat[i][j], h); /*update h*/
  
      lp_val(m, d, h, par_mat, lp); /*update lp*/
  
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
  free((void *) h);
  free((void *) rho);
  free((void *) l2);
  
  return ldlGL_val;
  
}

void logDualLGLWrapper( double * restrict n_total, /*inputs*/
     double * restrict n_samples, /*inputs*/
     double * restrict m, double * restrict d,
     double * restrict par, /*inputs*/
     double * restrict model, double * restrict x, /*inputs*/
     double * restrict lambda, double * restrict pen_g, /*inputs*/
     double * restrict ldlGL_val /*output*/)
  /* Calculating objective group lasso DEL function at a given parameter value.
   * Inputs:
   *   n_total -- total sample size;
   *   n_samples -- a vector of length m+1 specifying the size of each sample;
   *   m -- number of samples - 1;
   *   d -- dimension of h(x); 
   *   par_mat -- values of parameters (pointer matrix of dimension m by (d+1));
   *   h_func -- the basis function of the DRM;
   *   x_mat -- 2-D pointer array of data, organized as x_0, x_1, ..., x_m.
   *   lambda -- a positive threshold for the group lasso penalty.
   *   pen_g -- a pointer vector of length d for the group specific penalty
   * Outputs:
   *   ldlGL_val -- value of the objective function at a given "par" value.
   */
{
  /* loop indices */
  unsigned long i;
  
  unsigned long * restrict n_samples_use;
  n_samples_use = (unsigned long * restrict) malloc((size_t) (((unsigned
                                                                  long)*m + 1)*sizeof(unsigned long)));
  if (n_samples_use == NULL) {
    errMsg("malloc() allocation failure for m_samples_use!");
  }
  for (i = 0; i < ((unsigned long)*m + 1); ++i) {
    n_samples_use[i] = (unsigned long)n_samples[i];
  }
  
  double *restrict* restrict par_mat;
  double *restrict* restrict x_mat;
  
  /* converting x and par to matrices */
  par_mat = (double *restrict* restrict) malloc((size_t)
                                                  (((unsigned long)*m)*sizeof(double*)));
  if (par_mat == NULL) errMsg("malloc() allocation failure for par_mat!");
  par_mat[0] = par;
  for (i = 1; i < (unsigned long)*m; ++i){
    par_mat[i] = par_mat[i-1] + ((unsigned long)*d + 1);
  }
  
  x_mat = (double *restrict* restrict) malloc((size_t) (((unsigned
                                                            long)*m+1)*sizeof(double*)));
  if (x_mat == NULL) errMsg("malloc() allocation failure for x_mat!");
  x_mat[0] = x;
  for (i = 1; i < ((unsigned long)*m + 1); ++i){
    x_mat[i] = x_mat[i-1] + n_samples_use[i - 1];
  }
  
  /* calculating objective function at 'par' */
  switch ((unsigned long)*model)
  {
  case 1 :
    /*printf("h(x) = (x)");*/
    if ((unsigned long)*d != 1) {
      errMsg("For model 1, h(x) = x, d must be 1!");
    }
    *ldlGL_val = logDualLGL((unsigned long)*n_total,
                        /*n_samples_use, (unsigned long)*m, (unsigned long)*d,*/
                        n_samples_use, (unsigned long)*m, 1,
                        par_mat, &h1x, x_mat,
                        *lambda, pen_g);
    break;
    
  case 2 :
    /*printf("h(x) = (log(x))");*/
    if ((unsigned long)*d != 1) {
      errMsg("For model 2, h(x) = log(x), d must be 1!");
    }
    *ldlGL_val = logDualLGL((unsigned long)*n_total,
                        /*n_samples_use, (unsigned long)*m, (unsigned long)*d,*/
                        n_samples_use, (unsigned long)*m, 1,
                        par_mat, &h1logx, x_mat,
                        *lambda, pen_g);
    break;
    
  case 3 :
    /*printf("h(x) = (sqrt(x))");*/
    if ((unsigned long)*d != 1) {
      errMsg("For model 3, h(x) = sqrt(x), d must be 1!");
    }
    *ldlGL_val = logDualLGL((unsigned long)*n_total,
                        /*n_samples_use, (unsigned long)*m, (unsigned long)*d,*/
                        n_samples_use, (unsigned long)*m, 1,
                        par_mat, &h1sqrtx, x_mat,
                        *lambda, pen_g);
    break;
    
  case 4 :
    /*printf("h(x) = (x^2)");*/
    if ((unsigned long)*d != 1) {
      errMsg("For model 4, h(x) = x^2, d must be 1!");
    }
    *ldlGL_val = logDualLGL((unsigned long)*n_total,
                        /*n_samples_use, (unsigned long)*m, (unsigned long)*d,*/
                        n_samples_use, (unsigned long)*m, 1,
                        par_mat, &h1xSquare, x_mat,
                        *lambda, pen_g);
    break;
    
  case 5 :
    /*printf("h(x) = (x, x^2) -- Normal model");*/
    if ((unsigned long)*d != 2) {
      errMsg("For model 5 (Normal model), h(x) = (x, x^2), d must be 2!");
    }
    *ldlGL_val = logDualLGL((unsigned long)*n_total,
                        /*n_samples_use, (unsigned long)*m, (unsigned long)*d,*/
                        n_samples_use, (unsigned long)*m, 2,
                        par_mat, &h2Normal, x_mat,
                        *lambda, pen_g);
    break;
    
  case 6 :
    /*printf("h(x) = (x, log(x)) -- Gamma model");*/
    if ((unsigned long)*d != 2) {
      errMsg("For model 6 (Gamma model), h(x) = (x, log(x)), d must be 2!");
    }
    *ldlGL_val = logDualLGL((unsigned long)*n_total,
                        /*n_samples_use, (unsigned long)*m, (unsigned long)*d,*/
                        n_samples_use, (unsigned long)*m, 2,
                        par_mat, &h2Gamma, x_mat,
                        *lambda, pen_g);
    break;
    
  case 7 :
    /*printf("h(x) = (log(x), sqrt(x), x)");*/
    if ((unsigned long)*d != 3) {
      errMsg("For model 7, h(x) = (log(x), sqrt(x), x), d must be 3!");
    }
    *ldlGL_val = logDualLGL((unsigned long)*n_total,
                        /*n_samples_use, (unsigned long)*m, (unsigned long)*d,*/
                        n_samples_use, (unsigned long)*m, 3,
                        par_mat, &h3a, x_mat,
                        *lambda, pen_g);
    break;
    
  case 8 :
    /*printf("h(x) = (log(x), sqrt(x), x^2)");*/
    if ((unsigned long)*d != 3) {
      errMsg("For model 8, h(x) = (log(x), sqrt(x), x^2), d must be 3!");
    }
    *ldlGL_val = logDualLGL((unsigned long)*n_total,
                        /*n_samples_use, (unsigned long)*m, (unsigned long)*d,*/
                        n_samples_use, (unsigned long)*m, 3,
                        par_mat, &h3b, x_mat,
                        *lambda, pen_g);
    break;
    
  case 9 :
    /*printf("h(x) = (log(x), x, x^2)");*/
    if ((unsigned long)*d != 3) {
      errMsg("For model 9, h(x) = (log(x), x, x^2), d must be 3!");
    }
    *ldlGL_val = logDualLGL((unsigned long)*n_total,
                        /*n_samples_use, (unsigned long)*m, (unsigned long)*d,*/
                        n_samples_use, (unsigned long)*m, 3,
                        par_mat, &h3c, x_mat,
                        *lambda, pen_g);
    break;
    
  case 10 :
    /*printf("h(x) = (sqrt(x), x, x^2)");*/
    if ((unsigned long)*d != 3) {
      errMsg("For model 10, h(x) = (sqrt(x), x, x^2), d must be 3!");
    }
    *ldlGL_val = logDualLGL((unsigned long)*n_total,
                        /*n_samples_use, (unsigned long)*m, (unsigned long)*d,*/
                        n_samples_use, (unsigned long)*m, 3,
                        par_mat, &h3d, x_mat,
                        *lambda, pen_g);
    break;
    
  case 11 :
    /*printf("h(x) = (log(x), sqrt(x), x, x^2)");*/
    if ((unsigned long)*d != 4) {
      errMsg("For model 11, h(x) = (log(x), sqrt(x), x, x^2), d must be 4!");
    }
    *ldlGL_val = logDualLGL((unsigned long)*n_total,
                        /*n_samples_use, (unsigned long)*m, (unsigned long)*d,*/
                        n_samples_use, (unsigned long)*m, 4,
                        par_mat, &h4a, x_mat,
                        *lambda, pen_g);
    break;
    
  case 12 :
    /*printf("h(x) = (log(x), log(x)^2, sqrt(x), x, x^2)");*/
    if ((unsigned long)*d != 5) {
      errMsg("For model 11, h(x) = (log(x), log(x)^2, sqrt(x), x, x^2), d must be 5!");
    }
    *ldlGL_val = logDualLGL((unsigned long)*n_total,
                        /*n_samples_use, (unsigned long)*m, (unsigned long)*d,*/
                        n_samples_use, (unsigned long)*m, 5,
                        par_mat, &h5a, x_mat,
                        *lambda, pen_g);
    break;
    
  default :
    errMsg("'Model' must be an integer between 1 and 12 or a function of a single data point");
  break;
  
  }
  
  /* free arrays */
  free((void *) n_samples_use);
  free((void *) x_mat);
  free((void *) par_mat);
  
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

void logDualLGr(unsigned long n_total, /*inputs*/
    unsigned long * restrict n_samples, /*inputs*/
    unsigned long m, unsigned long d, /*inputs*/
    double *restrict* restrict par_mat, /*inputs*/
    void (*h_func)(double, double * restrict), /*input*/
    double *restrict* restrict x_mat, /*inputs*/
    double *restrict* restrict ldl_gr_mat /*outputs*/)
/* Calculating the gradient of log dual empirical likelihood (+ n \log n) at a given
 * parameter value.
 * Inputs:
 *   n_total -- total sample size;
 *   n_samples -- a vector of length m+1 specifying the size of each sample;
 *   m -- number of samples - 1;
 *   d -- dimension of h(x); 
 *   par_mat -- values of parameters (pointer matrix of dimension m by (d+1));
 *   h_func -- the basis function of the DRM;
 *   x_mat -- 2-D pointer array of data, organized as x_0, x_1, ..., x_m.
 * Outputs:
 *   ldl_gr_mat -- a 2-D pointer array of dimension m by (d+1); value of the gradient of the ldl (log dual empirical likelihood) at a given "par" value.
 */
{
  /* loop indices */
  unsigned long i, j, k, l;

  double * restrict R;
  R = (double * restrict) malloc((size_t) (m*sizeof(double)));
  if (R == NULL) errMsg("malloc() allocation failure for R!");
  for (i = 0; i < m; ++i) {
    R[i] = 0.0;
  }

  double * restrict H;
  H = (double * restrict) malloc((size_t) ((d+1)*sizeof(double)));
  if (H == NULL) errMsg("malloc() allocation failure for H!");
  H[0] = 1.0;
  for (i = 1; i < (d+1); ++i) {
    H[i] = 0.0;
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

  /*initializing ldl_gr_mat as safeguard*/
  for (i = 0; i < m; ++i) {
    for (j = 0; j < (d+1); ++j) {
      ldl_gr_mat[i][j] = 0.0;
    }
  }

  for (i = 0; i < m+1; ++i) {

    for (j = 0; j < n_samples[i]; ++j) {

      /*update H = (1, h^T)^T*/
      (*h_func)(x_mat[i][j], H+1); /*update H*/

      R_val(m, d, H+1, par_mat, rho, R); /*update R*/

      /*calculating S*/
      S = rho[0];
      for (k = 0; k < m; ++k) {
        S += R[k];
      }

      /* calculating the gradient of ldl */
      for (k = 0; k < m; ++k) {

        tmp_double = -R[k]/S;

        for (l = 0; l < (d+1); ++l) {
          ldl_gr_mat[k][l] += tmp_double * H[l];
        }

      }

      if (i > 0) {

        for (l = 0; l < (d+1); ++l) {
          ldl_gr_mat[i-1][l] = H[l] + ldl_gr_mat[i-1][l];
        }

      }

    }

  }

  /* free arrays */
  free((void *) R);
  free((void *) H);
  free((void *) rho);

}

void logDualLGrWrapper(double * restrict n_total, /*inputs*/
    double * restrict n_samples, /*inputs*/
    double * restrict m, double * restrict d,
    double * restrict par, /*inputs*/
    double * restrict model, double * restrict x, /*inputs*/
    double * restrict ldl_gr /*output*/)
/* Calculating log dual empirical likelihood (+ n \log n) at a given parameter value.
 * Inputs:
 *   n_total -- total sample size;
 *   n_samples -- a vector of length m+1 specifying the size of each sample;
 *   m -- number of samples - 1;
 *   d -- dimension of h(x); 
 *   par -- values of parameters (length of m(d+1)), organized as \theta_1, \cdots, \theta_m;
 *   h_func -- the basis function of the DRM;
 *   x -- data, organized as x_0, x_1, ..., x_m.
 * Outputs:
 *   ldl_gr -- a vector of length m*(d+1); value of the gradient of the ldl (log dual empirical likelihood) at a given "par" value.
 */
{
  /* loop indices */
  unsigned long i;

  unsigned long * restrict n_samples_use;
  n_samples_use = (unsigned long * restrict) malloc((size_t) (((unsigned
            long)*m + 1)*sizeof(unsigned long)));
  if (n_samples_use == NULL) {
    errMsg("malloc() allocation failure for m_samples_use!");
  }
  for (i = 0; i < ((unsigned long)*m + 1); ++i) {
    n_samples_use[i] = (unsigned long)n_samples[i];
  }

  double *restrict* restrict par_mat;
  double *restrict* restrict x_mat;
  double *restrict* restrict ldl_gr_mat;

  /* converting x, par and ldl_gr to matrices */
  par_mat = (double *restrict* restrict) malloc((size_t)
      (((unsigned long)*m)*sizeof(double*)));
  if (par_mat == NULL) errMsg("malloc() allocation failure for par_mat!");
  par_mat[0] = par;
  for (i = 1; i < (unsigned long)*m; ++i){
    par_mat[i] = par_mat[i-1] + ((unsigned long)*d + 1);
  }

  x_mat = (double *restrict* restrict) malloc((size_t) (((unsigned
            long)*m+1)*sizeof(double*)));
  if (x_mat == NULL) errMsg("malloc() allocation failure for x_mat!");
  x_mat[0] = x;
  for (i = 1; i < ((unsigned long)*m + 1); ++i){
    x_mat[i] = x_mat[i-1] + n_samples_use[i - 1];
  }

  ldl_gr_mat = (double *restrict* restrict) malloc((size_t)
      (((unsigned long)*m)*sizeof(double*)));
  if (ldl_gr_mat == NULL){
    errMsg("malloc() allocation failure for ldl_gr_mat!");
  }
  ldl_gr_mat[0] = ldl_gr;
  for (i = 1; i < (unsigned long)*m; ++i){
    ldl_gr_mat[i] = ldl_gr_mat[i-1] + ((unsigned long)*d + 1);
  }


  /* calculating log dual empirical likelihood value at 'par' */
  switch ((unsigned long)*model)
  {
    case 1 :
      /*printf("h(x) = (x)");*/
      if ((unsigned long)*d != 1) {
        errMsg("For model 1, h(x) = x, d must be 1!");
      }
      logDualLGr((unsigned long)*n_total, /*inputs*/
          /*n_samples_use, (unsigned long)*m, (unsigned long)*d, [>inputs<]*/
          n_samples_use, (unsigned long)*m, 1, /*inputs*/
          par_mat, &h1x, x_mat, /*inputs*/
          ldl_gr_mat /*outputs*/);
      break;

    case 2 :
      /*printf("h(x) = (log(x))");*/
      if ((unsigned long)*d != 1) {
        errMsg("For model 2, h(x) = log(x), d must be 1!");
      }
      logDualLGr((unsigned long)*n_total, /*inputs*/
          /*n_samples_use, (unsigned long)*m, (unsigned long)*d, [>inputs<]*/
          n_samples_use, (unsigned long)*m, 1, /*inputs*/
          par_mat, &h1logx, x_mat, /*inputs*/
          ldl_gr_mat /*outputs*/);
      break;

    case 3 :
      /*printf("h(x) = (sqrt(x))");*/
      if ((unsigned long)*d != 1) {
        errMsg("For model 3, h(x) = sqrt(x), d must be 1!");
      }
      logDualLGr((unsigned long)*n_total, /*inputs*/
          /*n_samples_use, (unsigned long)*m, (unsigned long)*d, [>inputs<]*/
          n_samples_use, (unsigned long)*m, 1, /*inputs*/
          par_mat, &h1sqrtx, x_mat, /*inputs*/
          ldl_gr_mat /*outputs*/);
      break;

    case 4 :
      /*printf("h(x) = (x^2)");*/
      if ((unsigned long)*d != 1) {
        errMsg("For model 4, h(x) = x^2, d must be 1!");
      }
      logDualLGr((unsigned long)*n_total, /*inputs*/
          /*n_samples_use, (unsigned long)*m, (unsigned long)*d, [>inputs<]*/
          n_samples_use, (unsigned long)*m, 1, /*inputs*/
          par_mat, &h1xSquare, x_mat, /*inputs*/
          ldl_gr_mat /*outputs*/);
      break;

    case 5 :
      /*printf("h(x) = (x, x^2) -- Normal model");*/
      if ((unsigned long)*d != 2) {
        errMsg("For model 5 (Normal model), h(x) = (x, x^2), d must be 2!");
      }
      logDualLGr((unsigned long)*n_total, /*inputs*/
          /*n_samples_use, (unsigned long)*m, (unsigned long)*d, [>inputs<]*/
          n_samples_use, (unsigned long)*m, 2, /*inputs*/
          par_mat, &h2Normal, x_mat, /*inputs*/
          ldl_gr_mat /*outputs*/);
      break;

    case 6 :
      /*printf("h(x) = (x, log(x)) -- Gamma model");*/
      if ((unsigned long)*d != 2) {
        errMsg("For model 6 (Gamma model), h(x) = (x, log(x)), d must be 2!");
      }
      logDualLGr((unsigned long)*n_total, /*inputs*/
          /*n_samples_use, (unsigned long)*m, (unsigned long)*d, [>inputs<]*/
          n_samples_use, (unsigned long)*m, 2, /*inputs*/
          par_mat, &h2Gamma, x_mat, /*inputs*/
          ldl_gr_mat /*outputs*/);
      break;

    case 7 :
      /*printf("h(x) = (log(x), sqrt(x), x)");*/
      if ((unsigned long)*d != 3) {
        errMsg("For model 7, h(x) = (log(x), sqrt(x), x), d must be 3!");
      }
      logDualLGr((unsigned long)*n_total, /*inputs*/
          /*n_samples_use, (unsigned long)*m, (unsigned long)*d, [>inputs<]*/
          n_samples_use, (unsigned long)*m, 3, /*inputs*/
          par_mat, &h3a, x_mat, /*inputs*/
          ldl_gr_mat /*outputs*/);
      break;

    case 8 :
      /*printf("h(x) = (log(x), sqrt(x), x^2)");*/
      if ((unsigned long)*d != 3) {
        errMsg("For model 8, h(x) = (log(x), sqrt(x), x^2), d must be 3!");
      }
      logDualLGr((unsigned long)*n_total, /*inputs*/
          /*n_samples_use, (unsigned long)*m, (unsigned long)*d, [>inputs<]*/
          n_samples_use, (unsigned long)*m, 3, /*inputs*/
          par_mat, &h3b, x_mat, /*inputs*/
          ldl_gr_mat /*outputs*/);
      break;

    case 9 :
      /*printf("h(x) = (log(x), x, x^2)");*/
      if ((unsigned long)*d != 3) {
        errMsg("For model 9, h(x) = (log(x), x, x^2), d must be 3!");
      }
      logDualLGr((unsigned long)*n_total, /*inputs*/
          /*n_samples_use, (unsigned long)*m, (unsigned long)*d, [>inputs<]*/
          n_samples_use, (unsigned long)*m, 3, /*inputs*/
          par_mat, &h3c, x_mat, /*inputs*/
          ldl_gr_mat /*outputs*/);
      break;

    case 10 :
      /*printf("h(x) = (sqrt(x), x, x^2)");*/
      if ((unsigned long)*d != 3) {
        errMsg("For model 10, h(x) = (sqrt(x), x, x^2), d must be 3!");
      }
      logDualLGr((unsigned long)*n_total, /*inputs*/
          /*n_samples_use, (unsigned long)*m, (unsigned long)*d, [>inputs<]*/
          n_samples_use, (unsigned long)*m, 3, /*inputs*/
          par_mat, &h3d, x_mat, /*inputs*/
          ldl_gr_mat /*outputs*/);
      break;

    case 11 :
      /*printf("h(x) = (log(x), sqrt(x), x, x^2)");*/
      if ((unsigned long)*d != 4) {
        errMsg("For model 11, h(x) = (log(x), sqrt(x), x, x^2), d must be 4!");
      }
      logDualLGr((unsigned long)*n_total, /*inputs*/
          /*n_samples_use, (unsigned long)*m, (unsigned long)*d, [>inputs<]*/
          n_samples_use, (unsigned long)*m, 4, /*inputs*/
          par_mat, &h4a, x_mat, /*inputs*/
          ldl_gr_mat /*outputs*/);
      break;
      
    case 12 :
      /*printf("h(x) = (log(x), log(x)^2, sqrt(x), x, x^2)");*/
      if ((unsigned long)*d != 5) {
        errMsg("For model 11, h(x) = (log(x), log(x)^2, sqrt(x), x, x^2), d must be 5!");
      }
      logDualLGr((unsigned long)*n_total, /*inputs*/
        /*n_samples_use, (unsigned long)*m, (unsigned long)*d, [>inputs<]*/
        n_samples_use, (unsigned long)*m, 5, /*inputs*/
        par_mat, &h5a, x_mat, /*inputs*/
        ldl_gr_mat /*outputs*/);
    break;

    default :
      errMsg("'Model' must be an integer between 1 and 11 or a function of a single data point");
      break;

  }

  /* free arrays */
  free((void *) n_samples_use);
  free((void *) x_mat);
  free((void *) par_mat);
  free((void *) ldl_gr_mat);

}

void logDualLHessian(unsigned long n_total, /*inputs*/
    double * restrict n_samples, /*inputs*/
    unsigned long m, unsigned long d, /*inputs*/
    double *restrict* restrict par_mat, /*inputs*/
    void (*h_func)(double, double * restrict), /*input*/
    double * restrict x, /*inputs*/
    double *restrict* restrict ldl_hessian_mat /*outputs*/)
/* Calculating the Hessian of log dual empirical likelihood (+ n \log n) at a given
 * parameter value.
 * Inputs:
 *   n_total -- length of data;
 *   n_samples -- a (double, not unsigned long!) vector of length m+1
 *     specifying the size of each sample;
 *   m -- number of samples - 1;
 *   d -- dimension of h(x); 
 *   par_mat -- values of parameters (pointer matrix of dimension m by (d+1));
 *   h_func -- the basis function of the DRM;
 *   x -- data, a long vector in order of as x_0, x_1, ..., x_m.
 * Outputs:
 *   ldl_hessian_mat -- the m(d+1) by m(d+1) Hessian matrix evaluated at the
 *     parameter value par_mat.
 */
{
  /* loop indices */
  unsigned long i, k, l, kh, lh;

  double * restrict R;
  R = (double * restrict) malloc((size_t) (m*sizeof(double)));
  if (R == NULL) errMsg("malloc() allocation failure for R!");
  for (i = 0; i < m; ++i) {
    R[i] = 0.0;
  }

  double * restrict H;
  H = (double * restrict) malloc((size_t) ((d+1)*sizeof(double)));
  if (H == NULL) errMsg("malloc() allocation failure for H!");
  H[0] = 1.0;
  for (i = 1; i < (d+1); ++i) {
    H[i] = 0.0;
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

  /*initializing ldl_hessian_mat as safeguard*/
  for (i = 0; i < m*(d+1); ++i) {
    for (k = 0; k < m*(d+1); ++k) {
      ldl_hessian_mat[i][k] = 0.0;
    }
  }

  for (i = 0; i < n_total; ++i) {

    /*update H = (1, h^T)^T*/
    (*h_func)(x[i], H+1); /*update H*/

    R_val(m, d, H+1, par_mat, n_samples, R); /*update R*/

    /*calculating S*/
    S = n_samples[0];
    for (k = 0; k < m; ++k) {
      S += R[k];
    }

    for (k = 0; k < m; ++k) {

      for (l = 0; l < m; ++l) {
        qaa[k][l] = R[k]*R[l]/(S*S);
      }

    }
    for (k = 0; k < m; ++k) {
        qaa[k][k] -= R[k]/S;
    }

    /*calculating the Hessian matrix m(d+1) by m(d+1)*/
    /*set very entry of ldl_hessian_mat to be 0;*/
    /*k and l are for vertical (row) indices; kh and lh are for horizontal
     * (column) indices.*/
    for (k = 0; k < m; ++k) {
      for (l = 0; l < (d+1); ++l) {

        for (kh = 0; kh < m; ++kh) {
          for (lh = 0; lh < (d+1); ++lh) {

            ldl_hessian_mat[(k*(d+1) + l)][(kh*(d+1) + lh)] +=
              qaa[k][kh]*H[l]*H[lh];

          }
        }

      }
    }

  }

  /* free arrays */
  free((void *) R);
  free((void *) H);

  free((void *) qaa[0]);
  free((void *) qaa);

}

void logDualLHessianWrapper(double * restrict n_total, /*inputs*/
    double * restrict n_samples, /*inputs*/
    double * restrict m, double * restrict d, /*inputs*/
    double * restrict par, double * restrict model, /*inputs*/
    double * restrict x, /*inputs*/
    double * restrict ldl_hessian /*output*/)
/* Calculating the Hessian of log dual empirical likelihood (+ n \log n) at a given
 * parameter value.
 * Inputs:
 *   n_total -- length of data;
 *   n_samples -- a vector of length m+1 specifying the size of each sample;
 *   m -- number of samples - 1;
 *   d -- dimension of h(x); 
 *   par -- values of parameters (length of m(d+1)), organized as \theta_1,
 *     \cdots, \theta_m;
 *   h_func -- the basis function of the DRM;
 *   x -- data, a long vector in order of as x_0, x_1, ..., x_m;
 * Outputs:
 *   ldl_hessian -- the vector form (row by row) of the m(d+1) by m(d+1)
 *     Hessian matrix evaluated at the parameter value par_mat.
 */
{
  /* loop indices */
  unsigned long i;

  double *restrict* restrict par_mat;
  double *restrict* restrict ldl_hessian_mat;

  /* converting par and ldl_hessian to matrices */
  par_mat = (double *restrict* restrict) malloc((size_t)
      (((unsigned long)*m)*sizeof(double*)));
  if (par_mat == NULL) errMsg("malloc() allocation failure for par_mat!");
  par_mat[0] = par;
  for (i = 1; i < (unsigned long)*m; ++i){
    par_mat[i] = par_mat[i-1] + ((unsigned long)*d + 1);
  }

  ldl_hessian_mat = (double *restrict* restrict) malloc((size_t)
      ((((unsigned long)*m) * ((unsigned long)*d + 1))*sizeof(double*)));
  if (ldl_hessian_mat == NULL){
    errMsg("malloc() allocation failure for ldl_hessian_mat!");
  }
  ldl_hessian_mat[0] = ldl_hessian;
  for (i = 1; i < (((unsigned long)*m) * ((unsigned long)*d + 1)); ++i){
    ldl_hessian_mat[i] = ldl_hessian_mat[i-1] + 
      (((unsigned long)*m) * ((unsigned long)*d + 1));
  }

  /* calculating log dual empirical likelihood value at 'par' */
  switch ((unsigned long)*model)
  {
    case 1 :
      /*printf("h(x) = (x)");*/
      if ((unsigned long)*d != 1) {
        errMsg("For model 1, h(x) = x, d must be 1!");
      }
      logDualLHessian((unsigned long)*n_total, n_samples, /*inputs*/
          (unsigned long)*m, 1, /*inputs*/
          par_mat, &h1x, x, /*inputs*/
          ldl_hessian_mat /*outputs*/);
      break;

    case 2 :
      /*printf("h(x) = (log(x))");*/
      if ((unsigned long)*d != 1) {
        errMsg("For model 2, h(x) = log(x), d must be 1!");
      }
      logDualLHessian((unsigned long)*n_total, n_samples, /*inputs*/
          (unsigned long)*m, 1, /*inputs*/
          par_mat, &h1logx, x, /*inputs*/
          ldl_hessian_mat /*outputs*/);
      break;

    case 3 :
      /*printf("h(x) = (sqrt(x))");*/
      if ((unsigned long)*d != 1) {
        errMsg("For model 3, h(x) = sqrt(x), d must be 1!");
      }
      logDualLHessian((unsigned long)*n_total, n_samples, /*inputs*/
          (unsigned long)*m, 1, /*inputs*/
          par_mat, &h1sqrtx, x, /*inputs*/
          ldl_hessian_mat /*outputs*/);
      break;

    case 4 :
      /*printf("h(x) = (x^2)");*/
      if ((unsigned long)*d != 1) {
        errMsg("For model 4, h(x) = x^2, d must be 1!");
      }
      logDualLHessian((unsigned long)*n_total, n_samples, /*inputs*/
          (unsigned long)*m, 1, /*inputs*/
          par_mat, &h1xSquare, x, /*inputs*/
          ldl_hessian_mat /*outputs*/);
      break;

    case 5 :
      /*printf("h(x) = (x, x^2) -- Normal model");*/
      if ((unsigned long)*d != 2) {
        errMsg("For model 5 (Normal model), h(x) = (x, x^2), d must be 2!");
      }
      logDualLHessian((unsigned long)*n_total, n_samples, /*inputs*/
          (unsigned long)*m, 2, /*inputs*/
          par_mat, &h2Normal, x, /*inputs*/
          ldl_hessian_mat /*outputs*/);
      break;

    case 6 :
      /*printf("h(x) = (x, log(x)) -- Gamma model");*/
      if ((unsigned long)*d != 2) {
        errMsg("For model 6 (Gamma model), h(x) = (x, log(x)), d must be 2!");
      }
      logDualLHessian((unsigned long)*n_total, n_samples, /*inputs*/
          (unsigned long)*m, 2, /*inputs*/
          par_mat, &h2Gamma, x, /*inputs*/
          ldl_hessian_mat /*outputs*/);
      break;

    case 7 :
      /*printf("h(x) = (log(x), sqrt(x), x)");*/
      if ((unsigned long)*d != 3) {
        errMsg("For model 7, h(x) = (log(x), sqrt(x), x), d must be 3!");
      }
      logDualLHessian((unsigned long)*n_total, n_samples, /*inputs*/
          (unsigned long)*m, 3, /*inputs*/
          par_mat, &h3a, x, /*inputs*/
          ldl_hessian_mat /*outputs*/);
      break;

    case 8 :
      /*printf("h(x) = (log(x), sqrt(x), x^2)");*/
      if ((unsigned long)*d != 3) {
        errMsg("For model 8, h(x) = (log(x), sqrt(x), x^2), d must be 3!");
      }
      logDualLHessian((unsigned long)*n_total, n_samples, /*inputs*/
          (unsigned long)*m, 3, /*inputs*/
          par_mat, &h3b, x, /*inputs*/
          ldl_hessian_mat /*outputs*/);
      break;

    case 9 :
      /*printf("h(x) = (log(x), x, x^2)");*/
      if ((unsigned long)*d != 3) {
        errMsg("For model 9, h(x) = (log(x), x, x^2), d must be 3!");
      }
      logDualLHessian((unsigned long)*n_total, n_samples, /*inputs*/
          (unsigned long)*m, 3, /*inputs*/
          par_mat, &h3c, x, /*inputs*/
          ldl_hessian_mat /*outputs*/);
      break;

    case 10 :
      /*printf("h(x) = (sqrt(x), x, x^2)");*/
      if ((unsigned long)*d != 3) {
        errMsg("For model 10, h(x) = (sqrt(x), x, x^2), d must be 3!");
      }
      logDualLHessian((unsigned long)*n_total, n_samples, /*inputs*/
          (unsigned long)*m, 3, /*inputs*/
          par_mat, &h3d, x, /*inputs*/
          ldl_hessian_mat /*outputs*/);
      break;

    case 11 :
      /*printf("h(x) = (log(x), sqrt(x), x, x^2)");*/
      if ((unsigned long)*d != 4) {
        errMsg("For model 11, h(x) = (log(x), sqrt(x), x, x^2), d must be 4!");
      }
      logDualLHessian((unsigned long)*n_total, n_samples, /*inputs*/
          (unsigned long)*m, 4, /*inputs*/
          par_mat, &h4a, x, /*inputs*/
          ldl_hessian_mat /*outputs*/);
      break;
      
    case 12 :
      /*printf("h(x) = (log(x), log(x)^2, sqrt(x), x, x^2)");*/
      if ((unsigned long)*d != 5) {
        errMsg("For model 11, h(x) = (log(x), log(x)^2, sqrt(x), x, x^2), d must be 5!");
      }
      logDualLHessian((unsigned long)*n_total, n_samples, /*inputs*/
        (unsigned long)*m, 5, /*inputs*/
        par_mat, &h5a, x, /*inputs*/
        ldl_hessian_mat /*outputs*/);
      break;

    default :
      errMsg("'Model' must be an integer between 1 and 11 or a function of a single data point");
      break;

  }

  /* free arrays */
  free((void *) par_mat);
  free((void *) ldl_hessian_mat);

}

void probEst(unsigned long n_total, /*inputs*/
    double * restrict n_samples, /*inputs*/
    unsigned long m, unsigned long d, /*inputs*/
    double *restrict* restrict par_mat, /*inputs*/
    void (*h_func)(double, double * restrict), /*input*/
    double * restrict x, /*inputs*/
    unsigned short normalize,
    double *restrict* restrict pEst_mat /*output*/)
/* Estimating \hat{p} for all samples in density ratio models. \hat{p}'s are
 * useful for estimating quantiles, distribution function and density of each
 * population.
 * Inputs:
 *   n_total -- total sample size;
 *   n_samples -- a (double, not unsigned long!) vector of length m+1
 *     specifying the size of each sample;
 *   m -- number of samples - 1;
 *   d -- dimension of h(x); 
 *   par_mat -- values of parameters (pointer matrix of dimension m by (d+1));
 *   h_func -- the basis function of the DRM;
 *   x -- data, a long vector in order of as x_0, x_1, ..., x_m.
 *   normalize -- indicator; wether normalize the probabilies so that the sum
 *     is strictly one: 1 for yes, 0 for no.
 * Outputs:
 *   pEst_mat -- a (m+1) by n_total matrix; estimated probability for x's at
 *     parameter value 'par_mat' for each sample.
 */
{
  /* loop indices */
  unsigned long i, j;

  double * restrict r;
  r = (double * restrict) malloc((size_t) (m*sizeof(double)));
  if (r == NULL) errMsg("malloc() allocation failure for r!");
  for (i = 0; i < m; ++i) {
    r[i] = 0;
  }

  double * restrict h;
  h = (double * restrict) malloc((size_t) (d*sizeof(double)));
  if (h == NULL) errMsg("malloc() allocation failure for h!");
  for (i = 0; i < d; ++i) {
    h[i] = 0;
  }

  double p_bl_tmp;

  /* auxiliary variable, used only when normalize == 1 */
  double * restrict pEst_sum;
  pEst_sum = (double * restrict) malloc((size_t) ((m+1)*sizeof(double)));
  if (pEst_sum == NULL) errMsg("malloc() allocation failure for pEst_sum!");
  for (i = 0; i < m+1; ++i) {
    pEst_sum[i] = 0.0;
  }

  /* initializing output variable pEst_mat to 0 as safeguard */
  for (i = 0; i < (m + 1); ++i) {
    for (j = 0; j < n_total; ++j) {
        pEst_mat[i][j] = 0.0;
    }
  }

  for (i = 0; i < n_total; ++i) {

    (*h_func)(x[i], h); /*update h*/

    r_val(m, d, h, par_mat, r); /*update r*/

    /* calculate baseline probabilities */
    p_bl_tmp = n_samples[0];
    for (j = 0; j < m; ++j) {
      p_bl_tmp += n_samples[j+1] * r[j];
    }
    p_bl_tmp = 1.0/p_bl_tmp;

    /* calculate probabilities for all m+1 sampels*/
    pEst_mat[0][i] = p_bl_tmp; /*baseline probablities*/
    if (normalize==1) pEst_sum[0] += pEst_mat[0][i];
    for (j = 1; j < m+1; ++j) {
      pEst_mat[j][i] = r[j-1]*p_bl_tmp;
      if (normalize==1) pEst_sum[j] += pEst_mat[j][i];
    }

  }

  if (normalize==1) {

    for (i = 0; i < m+1; ++i) {

      for (j = 0; j < n_total; ++j) {
          pEst_mat[i][j] /= pEst_sum[i];
      }

    }

  }

  /* free arrays */
  free((void *) r);
  free((void *) h);
  free((void *) pEst_sum);

}

void probEstWrapper(double * restrict n_total, /*inputs*/
    double * restrict n_samples, /*inputs*/
    double * restrict m, double * restrict d,
    double * restrict par, /*inputs*/
    double * restrict model, double * restrict x, /*inputs*/
    double * restrict normalize,
    double * restrict pEst /*output*/)
/* Estimating \hat{p} for all samples in density ratio models. \hat{p}'s are
 * useful for estimating quantiles, distribution function and density of each
 * population.
 * Inputs:
 *   n_total -- total sample size;
 *   n_samples -- a vector of length m+1 specifying the size of each sample;
 *   m -- number of samples - 1;
 *   d -- dimension of h(x); 
 *   par_mat -- values of parameters (pointer matrix of dimension m by (d+1));
 *   h_func -- the basis function of the DRM;
 *   x_mat -- 2-D pointer array of data, organized as x_0, x_1, ..., x_m.
 *   normalize -- indicator; wether normalize the probabilies so that the sum
 *     is strictly one: 1 for yes, 0 for no.
 * Outputs:
 *   pEst -- a vector of length (m+1)*n_total; estimated probability for x's at
 *     parameter value 'par_mat' for each sample.
 */
{
  /* loop indices */
  unsigned long i;

  double *restrict* restrict par_mat;
  /* converting x and par to matrices */
  par_mat = (double *restrict* restrict) malloc((size_t)
      (((unsigned long)*m)*sizeof(double*)));
  if (par_mat == NULL) errMsg("malloc() allocation failure for par_mat!");
  par_mat[0] = par;
  for (i = 1; i < (unsigned long)*m; ++i){
    par_mat[i] = par_mat[i-1] + ((unsigned long)*d + 1);
  }

  double *restrict* restrict pEst_mat;
  pEst_mat = (double *restrict* restrict) malloc((size_t)
      (((unsigned long)*m+1)*sizeof(double*)));
  if (pEst_mat == NULL) errMsg("malloc() allocation failure for pEst_mat!");
  pEst_mat[0] = pEst;
  for (i = 1; i < ((unsigned long)*m + 1); ++i){
    pEst_mat[i] = pEst_mat[i-1] + (unsigned long)*n_total;
  }

  /* calculating estimated probability for baseline at a given 'par' */
  switch ((unsigned long)*model)
  {
    case 1 :
      /*printf("h(x) = (x)");*/
      if ((unsigned long)*d != 1) {
        errMsg("For model 1, h(x) = x, d must be 1!");
      }
      probEst((unsigned long)*n_total, n_samples, /*inputs*/
          (unsigned long) *m, 1, /*inputs*/
          par_mat, &h1x, x, (unsigned short) *normalize, /*inputs*/
          pEst_mat /*output*/);
      break;

    case 2 :
      /*printf("h(x) = (log(x))");*/
      if ((unsigned long)*d != 1) {
        errMsg("For model 2, h(x) = log(x), d must be 1!");
      }
      probEst((unsigned long)*n_total, n_samples, /*inputs*/
          (unsigned long) *m, 1, /*inputs*/
          par_mat, &h1logx, x, (unsigned short) *normalize, /*inputs*/
          pEst_mat /*output*/);
      break;

    case 3 :
      /*printf("h(x) = (sqrt(x))");*/
      if ((unsigned long)*d != 1) {
        errMsg("For model 3, h(x) = sqrt(x), d must be 1!");
      }
      probEst((unsigned long)*n_total, n_samples, /*inputs*/
          (unsigned long) *m, 1, /*inputs*/
          par_mat, &h1sqrtx, x, (unsigned short) *normalize, /*inputs*/
          pEst_mat /*output*/);
      break;

    case 4 :
      /*printf("h(x) = (x^2)");*/
      if ((unsigned long)*d != 1) {
        errMsg("For model 4, h(x) = x^2, d must be 1!");
      }
      probEst((unsigned long)*n_total, n_samples, /*inputs*/
          (unsigned long) *m, 1, /*inputs*/
          par_mat, &h1xSquare, x, (unsigned short) *normalize, /*inputs*/
          pEst_mat /*output*/);
      break;

    case 5 :
      /*printf("h(x) = (x, x^2) -- Normal model");*/
      if ((unsigned long)*d != 2) {
        errMsg("For model 5 (Normal model), h(x) = (x, x^2), d must be 2!");
      }
      probEst((unsigned long)*n_total, n_samples, /*inputs*/
          (unsigned long) *m, 2, /*inputs*/
          par_mat, &h2Normal, x, (unsigned short) *normalize, /*inputs*/
          pEst_mat /*output*/);
      break;

    case 6 :
      /*printf("h(x) = (x, log(x)) -- Gamma model");*/
      if ((unsigned long)*d != 2) {
        errMsg("For model 6 (Gamma model), h(x) = (x, log(x)), d must be 2!");
      }
      probEst((unsigned long)*n_total, n_samples, /*inputs*/
          (unsigned long) *m, 2, /*inputs*/
          par_mat, &h2Gamma, x, (unsigned short) *normalize, /*inputs*/
          pEst_mat /*output*/);
      break;

    case 7 :
      /*printf("h(x) = (log(x), sqrt(x), x)");*/
      if ((unsigned long)*d != 3) {
        errMsg("For model 7, h(x) = (log(x), sqrt(x), x), d must be 3!");
      }
      probEst((unsigned long)*n_total, n_samples, /*inputs*/
          (unsigned long) *m, 3, /*inputs*/
          par_mat, &h3a, x, (unsigned short) *normalize, /*inputs*/
          pEst_mat /*output*/);
      break;

    case 8 :
      /*printf("h(x) = (log(x), sqrt(x), x^2)");*/
      if ((unsigned long)*d != 3) {
        errMsg("For model 8, h(x) = (log(x), sqrt(x), x^2), d must be 3!");
      }
      probEst((unsigned long)*n_total, n_samples, /*inputs*/
          (unsigned long) *m, 3, /*inputs*/
          par_mat, &h3b, x, (unsigned short) *normalize, /*inputs*/
          pEst_mat /*output*/);
      break;

    case 9 :
      /*printf("h(x) = (log(x), x, x^2)");*/
      if ((unsigned long)*d != 3) {
        errMsg("For model 9, h(x) = (log(x), x, x^2), d must be 3!");
      }
      probEst((unsigned long)*n_total, n_samples, /*inputs*/
          (unsigned long) *m, 3, /*inputs*/
          par_mat, &h3c, x, (unsigned short) *normalize, /*inputs*/
          pEst_mat /*output*/);
      break;

    case 10 :
      /*printf("h(x) = (sqrt(x), x, x^2)");*/
      if ((unsigned long)*d != 3) {
        errMsg("For model 10, h(x) = (sqrt(x), x, x^2), d must be 3!");
      }
      probEst((unsigned long)*n_total, n_samples, /*inputs*/
          (unsigned long) *m, 3, /*inputs*/
          par_mat, &h3d, x, (unsigned short) *normalize, /*inputs*/
          pEst_mat /*output*/);
      break;

    case 11 :
      /*printf("h(x) = (log(x), sqrt(x), x, x^2)");*/
      if ((unsigned long)*d != 4) {
        errMsg("For model 11, h(x) = (log(x), sqrt(x), x, x^2), d must be 4!");
      }
      probEst((unsigned long)*n_total, n_samples, /*inputs*/
          (unsigned long) *m, 4, /*inputs*/
          par_mat, &h4a, x, (unsigned short) *normalize, /*inputs*/
          pEst_mat /*output*/);
      break;
    
    case 12 :
      /*printf("h(x) = (log(x), log(x)^2, sqrt(x), x, x^2)");*/
      if ((unsigned long)*d != 5) {
        errMsg("For model 11, h(x) = (log(x), log(x)^2, sqrt(x), x, x^2), d must be 5!");
      }
      probEst((unsigned long)*n_total, n_samples, /*inputs*/
        (unsigned long) *m, 5, /*inputs*/
        par_mat, &h5a, x, (unsigned short) *normalize, /*inputs*/
        pEst_mat /*output*/);
      break;

    default :
      errMsg("'Model' must be an integer between 1 and 11 or a function of a single data point");
      break;

  }

  /* free arrays */
  free((void *) par_mat);
  free((void *) pEst_mat);

}

/* ***** User specified basis function (Uf) version ***** */

double logDualLUf(unsigned long n_total, /*inputs*/
    unsigned long * restrict n_samples, /*inputs*/
    unsigned long m, unsigned long d, /*inputs*/
    double *restrict* restrict par_mat, /*inputs*/
    SEXP h_func, SEXP env, /*input*/
    double *restrict* restrict x_mat /*inputs*/)
/* Calculating log dual empirical likelihood (+ n \log n) at a given parameter value.
 * Inputs:
 *   n_total -- total sample size;
 *   n_samples -- a vector of length m+1 specifying the size of each sample;
 *   m -- number of samples - 1;
 *   d -- dimension of h(x); 
 *   par_mat -- values of parameters (pointer matrix of dimension m by (d+1));
 *   h_func -- the basis function of the DRM;
 *   x_mat -- 2-D pointer array of data, organized as x_0, x_1, ..., x_m.
 * Outputs:
 *   ldl_val -- value of ldl (log dual empirical likelihood) at a given "par" value.
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

  double * restrict h;
  h = (double * restrict) malloc((size_t) (d*sizeof(double)));
  if (h == NULL) errMsg("malloc() allocation failure for h!");
  for (i = 0; i < d; ++i) {
    h[i] = 0;
  }
  SEXP h_args, h_r_call, h_r_call_result;
  PROTECT(h_args = allocVector(REALSXP, 1));

  double * restrict rho;
  rho = (double * restrict) malloc((size_t) ((m+1)*sizeof(double)));
  if (rho == NULL) errMsg("malloc() allocation failure for rho!");
  for (i = 0; i < m+1; ++i) {
    rho[i] = (double)n_samples[i]/(double)n_total;
  }

  /*other variables*/
  double S;

  /*define output*/
  double ldl_val;
  ldl_val = 0.0;

  for (i = 0; i < m+1; ++i) {

    for (j = 0; j < n_samples[i]; ++j) {

      /*update h*/
      REAL(h_args)[0] = x_mat[i][j];
      PROTECT(h_r_call = lang2(h_func, h_args));
      PROTECT(h_r_call_result = eval(h_r_call, env));
      for (k = 0; k < d; ++k) {
        h[k] = REAL(h_r_call_result)[k];
      }
      UNPROTECT(2);

      lp_val(m, d, h, par_mat, lp); /*update lp*/

      /* calculating q_i */
      S = rho[0];
      for (k = 0; k < m; ++k) {
        S += rho[k+1] * exp(lp[k]);
      }

      if (i == 0) {
        ldl_val = ldl_val - log(S);  
      } else {
        ldl_val = ldl_val + lp[i-1] - log(S);
      }

    }

  }

  UNPROTECT(1);

  /* free arrays */
  free((void *) lp);
  free((void *) h);
  free((void *) rho);

  return ldl_val;

}

SEXP logDualLUfWrapper(SEXP n_total, SEXP n_samples, SEXP m, SEXP d, SEXP par,
    SEXP h_func, SEXP env, SEXP x)
/* Calculating log dual empirical likelihood (+ n \log n) at a given parameter value.
 * Inputs:
 *   n_total -- total sample size;
 *   n_samples -- a vector of length m+1 specifying the size of each sample;
 *   m -- number of samples - 1;
 *   d -- dimension of h(x); 
 *   par -- values of parameters (length of m(d+1)), organized as \theta_1, \cdots, \theta_m;
 *   h_func -- the basis function of the DRM;
 *   x -- data, organized as x_0, x_1, ..., x_m.
 * Outputs:
 *   ldl_val -- value of ldl (log dual empirical likelihood) at a given "par" value.
 */
{
  double * restrict n_total_c;
  n_total_c = REAL(n_total);
  double * restrict n_samples_c;
  n_samples_c = REAL(n_samples);
  double * restrict m_c;
  m_c = REAL(m);
  double * restrict d_c;
  d_c = REAL(d);
  double * restrict par_c;
  par_c = REAL(par);
  double * restrict x_c;
  x_c = REAL(x);

  /* loop indices */
  unsigned long i;

  double ldl_val_c;
  SEXP ldl_val;

  unsigned long * restrict n_samples_use;
  n_samples_use = (unsigned long * restrict) malloc((size_t) (((unsigned
            long)*m_c + 1)*sizeof(unsigned long)));
  if (n_samples_use == NULL) {
    errMsg("malloc() allocation failure for m_samples_use!");
  }
  for (i = 0; i < ((unsigned long)*m_c + 1); ++i) {
    n_samples_use[i] = (unsigned long)n_samples_c[i];
  }

  double *restrict* restrict par_mat;
  double *restrict* restrict x_mat;

  /* converting x and par to matrices */
  par_mat = (double *restrict* restrict) malloc((size_t)
      (((unsigned long)*m_c)*sizeof(double*)));
  if (par_mat == NULL) errMsg("malloc() allocation failure for par_mat!");
  par_mat[0] = par_c;
  for (i = 1; i < (unsigned long)*m_c; ++i){
    par_mat[i] = par_mat[i-1] + ((unsigned long)*d_c + 1);
  }

  x_mat = (double *restrict* restrict) malloc((size_t) (((unsigned
            long)*m_c+1)*sizeof(double*)));
  if (x_mat == NULL) errMsg("malloc() allocation failure for x_mat!");
  x_mat[0] = x_c;
  for (i = 1; i < ((unsigned long)*m_c + 1); ++i){
    x_mat[i] = x_mat[i-1] + n_samples_use[i - 1];
  }

  ldl_val_c = logDualLUf((unsigned long)*n_total_c, n_samples_use,
      (unsigned long)*m_c, (unsigned long)*d_c, par_mat,
      h_func, env, x_mat);

  PROTECT(ldl_val = allocVector(REALSXP, 1));
  REAL(ldl_val)[0] = ldl_val_c;

  /* free arrays */
  free((void *) n_samples_use);
  free((void *) x_mat);
  free((void *) par_mat);

  UNPROTECT(1);

  return(ldl_val);

}

void logDualLGrUf(unsigned long n_total, /*inputs*/
    unsigned long * restrict n_samples, /*inputs*/
    unsigned long m, unsigned long d, /*inputs*/
    double *restrict* restrict par_mat, /*inputs*/
    SEXP h_func, SEXP env, /*input*/
    double *restrict* restrict x_mat, /*inputs*/
    double *restrict* restrict ldl_gr_mat /*outputs*/)
/* Calculating the gradient of log dual empirical likelihood (+ n \log n) at a given
 * parameter value.
 * Inputs:
 *   n_total -- total sample size;
 *   n_samples -- a vector of length m+1 specifying the size of each sample;
 *   m -- number of samples - 1;
 *   d -- dimension of h(x); 
 *   par_mat -- values of parameters (pointer matrix of dimension m by (d+1));
 *   h_func -- the basis function of the DRM;
 *   x_mat -- 2-D pointer array of data, organized as x_0, x_1, ..., x_m.
 * Outputs:
 *   ldl_gr_mat -- a 2-D pointer array of dimension m by (d+1); value of the gradient of the ldl (log dual empirical likelihood) at a given "par" value.
 */
{
  /* loop indices */
  unsigned long i, j, k, l;

  double * restrict R;
  R = (double * restrict) malloc((size_t) (m*sizeof(double)));
  if (R == NULL) errMsg("malloc() allocation failure for R!");
  for (i = 0; i < m; ++i) {
    R[i] = 0.0;
  }

  double * restrict H;
  H = (double * restrict) malloc((size_t) ((d+1)*sizeof(double)));
  if (H == NULL) errMsg("malloc() allocation failure for H!");
  H[0] = 1.0;
  for (i = 1; i < (d+1); ++i) {
    H[i] = 0.0;
  }
  SEXP h_args, h_r_call, h_r_call_result;
  PROTECT(h_args = allocVector(REALSXP, 1));

  double * restrict rho;
  rho = (double * restrict) malloc((size_t) ((m+1)*sizeof(double)));
  if (rho == NULL) errMsg("malloc() allocation failure for rho!");
  for (i = 0; i < m+1; ++i) {
    rho[i] = (double)n_samples[i]/(double)n_total;
  }

  /*other variables*/
  double S;
  double tmp_double;

  /*initializing ldl_gr_mat as safeguard*/
  for (i = 0; i < m; ++i) {
    for (j = 0; j < (d+1); ++j) {
      ldl_gr_mat[i][j] = 0.0;
    }
  }

  for (i = 0; i < m+1; ++i) {

    for (j = 0; j < n_samples[i]; ++j) {

      /*update H = (1, h^T)^T*/
      REAL(h_args)[0] = x_mat[i][j];
      PROTECT(h_r_call = lang2(h_func, h_args));
      PROTECT(h_r_call_result = eval(h_r_call, env));
      for (k = 0; k < d; ++k) {
        H[k+1] = REAL(h_r_call_result)[k];
      }
      UNPROTECT(2);

      R_val(m, d, H+1, par_mat, rho, R); /*update R*/

      /*calculating S*/
      S = rho[0];
      for (k = 0; k < m; ++k) {
        S += R[k];
      }

      /* calculating the gradient of ldl */
      for (k = 0; k < m; ++k) {

        tmp_double = -R[k]/S;

        for (l = 0; l < (d+1); ++l) {
          ldl_gr_mat[k][l] += tmp_double * H[l];
        }

      }

      if (i > 0) {

        for (l = 0; l < (d+1); ++l) {
          ldl_gr_mat[i-1][l] = H[l] + ldl_gr_mat[i-1][l];
        }

      }

    }

  }

  UNPROTECT(1);

  /* free arrays */
  free((void *) R);
  free((void *) H);
  free((void *) rho);

}

SEXP logDualLGrUfWrapper(SEXP n_total, SEXP n_samples, SEXP m, SEXP d,
    SEXP par, SEXP h_func, SEXP env, SEXP x /*input*/)
/* Calculating log dual empirical likelihood (+ n \log n) at a given parameter value.
 * Inputs:
 *   n_total -- total sample size;
 *   n_samples -- a vector of length m+1 specifying the size of each sample;
 *   m -- number of samples - 1;
 *   d -- dimension of h(x); 
 *   par -- values of parameters (length of m(d+1)), organized as \theta_1, \cdots, \theta_m;
 *   h_func -- the basis function of the DRM;
 *   x -- data, organized as x_0, x_1, ..., x_m.
 */
{
  double * restrict n_total_c;
  n_total_c = REAL(n_total);
  double * restrict n_samples_c;
  n_samples_c = REAL(n_samples);
  double * restrict m_c;
  m_c = REAL(m);
  double * restrict d_c;
  d_c = REAL(d);
  double * restrict par_c;
  par_c = REAL(par);
  double * restrict x_c;
  x_c = REAL(x);

  /* Outputs:
   *   ldl_gr -- a vector of length m*(d+1); value of the gradient of the ldl (log dual empirical likelihood) at a given "par" value.
   */
  SEXP ldl_gr;
  PROTECT(ldl_gr = allocVector(REALSXP, ((unsigned long)*m_c) * ((unsigned
            long)*d_c + 1)));
  double * restrict ldl_gr_c;
  ldl_gr_c = REAL(ldl_gr);

  /* loop indices */
  unsigned long i;

  unsigned long * restrict n_samples_use;
  n_samples_use = (unsigned long * restrict) malloc((size_t) (((unsigned
            long)*m_c + 1)*sizeof(unsigned long)));
  if (n_samples_use == NULL) {
    errMsg("malloc() allocation failure for m_samples_use!");
  }
  for (i = 0; i < ((unsigned long)*m_c + 1); ++i) {
    n_samples_use[i] = (unsigned long)n_samples_c[i];
  }

  double *restrict* restrict par_mat;
  double *restrict* restrict x_mat;
  double *restrict* restrict ldl_gr_mat;

  /* converting x, par and ldl_gr to matrices */
  par_mat = (double *restrict* restrict) malloc((size_t)
      (((unsigned long)*m_c)*sizeof(double*)));
  if (par_mat == NULL) errMsg("malloc() allocation failure for par_mat!");
  par_mat[0] = par_c;
  for (i = 1; i < (unsigned long)*m_c; ++i){
    par_mat[i] = par_mat[i-1] + ((unsigned long)*d_c + 1);
  }

  x_mat = (double *restrict* restrict) malloc((size_t) (((unsigned
            long)*m_c+1)*sizeof(double*)));
  if (x_mat == NULL) errMsg("malloc() allocation failure for x_mat!");
  x_mat[0] = x_c;
  for (i = 1; i < ((unsigned long)*m_c + 1); ++i){
    x_mat[i] = x_mat[i-1] + n_samples_use[i - 1];
  }

  ldl_gr_mat = (double *restrict* restrict) malloc((size_t)
      (((unsigned long)*m_c)*sizeof(double*)));
  if (ldl_gr_mat == NULL){
    errMsg("malloc() allocation failure for ldl_gr_mat!");
  }
  ldl_gr_mat[0] = ldl_gr_c;
  for (i = 1; i < (unsigned long)*m_c; ++i){
    ldl_gr_mat[i] = ldl_gr_mat[i-1] + ((unsigned long)*d_c + 1);
  }

  logDualLGrUf((unsigned long)*n_total_c, /*inputs*/
      n_samples_use, (unsigned long)*m_c, (unsigned long)*d_c, /*inputs*/
      par_mat, h_func, env, x_mat, /*inputs*/
      ldl_gr_mat /*outputs*/);

  /* free arrays */
  free((void *) n_samples_use);
  free((void *) x_mat);
  free((void *) par_mat);
  free((void *) ldl_gr_mat);

  UNPROTECT(1);

  return(ldl_gr);

}

void logDualLHessianUf(unsigned long n_total, /*inputs*/
    double * restrict n_samples, /*inputs*/
    unsigned long m, unsigned long d, /*inputs*/
    double *restrict* restrict par_mat, /*inputs*/
    SEXP h_func, SEXP env, /*input*/
    double * restrict x, /*inputs*/
    double *restrict* restrict ldl_hessian_mat /*outputs*/)
/* Calculating the Hessian of log dual emprical likelihood (+ n \log n) at a given
 * parameter value.
 * Inputs:
 *   n_total -- length of data;
 *   n_samples -- a (double, not unsigned long!) vector of length m+1
 *     specifying the size of each sample;
 *   m -- number of samples - 1;
 *   d -- dimension of h(x); 
 *   par_mat -- values of parameters (pointer matrix of dimension m by (d+1));
 *   h_func -- the basis function of the DRM;
 *   x -- data, a long vector in order of as x_0, x_1, ..., x_m.
 * Outputs:
 *   ldl_hessian_mat -- the m(d+1) by m(d+1) Hessian matrix evaluated at the
 *     parameter value par_mat.
 */
{
  /* loop indices */
  unsigned long i, k, l, kh, lh;

  double * restrict R;
  R = (double * restrict) malloc((size_t) (m*sizeof(double)));
  if (R == NULL) errMsg("malloc() allocation failure for R!");
  for (i = 0; i < m; ++i) {
    R[i] = 0.0;
  }

  double * restrict H;
  H = (double * restrict) malloc((size_t) ((d+1)*sizeof(double)));
  if (H == NULL) errMsg("malloc() allocation failure for H!");
  H[0] = 1.0;
  for (i = 1; i < (d+1); ++i) {
    H[i] = 0.0;
  }
  SEXP h_args, h_r_call, h_r_call_result;
  PROTECT(h_args = allocVector(REALSXP, 1));

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

  /*initializing ldl_hessian_mat as safeguard*/
  for (i = 0; i < m*(d+1); ++i) {
    for (k = 0; k < m*(d+1); ++k) {
      ldl_hessian_mat[i][k] = 0.0;
    }
  }

  for (i = 0; i < n_total; ++i) {

    /*update H = (1, h^T)^T*/
    REAL(h_args)[0] = x[i];
    PROTECT(h_r_call = lang2(h_func, h_args));
    PROTECT(h_r_call_result = eval(h_r_call, env));
    for (k = 0; k < d; ++k) {
      H[k+1] = REAL(h_r_call_result)[k];
    }
    UNPROTECT(2);

    R_val(m, d, H+1, par_mat, n_samples, R); /*update R*/

    /*calculating S*/
    S = n_samples[0];
    for (k = 0; k < m; ++k) {
      S += R[k];
    }

    for (k = 0; k < m; ++k) {

      for (l = 0; l < m; ++l) {
        qaa[k][l] = R[k]*R[l]/(S*S);
      }

    }
    for (k = 0; k < m; ++k) {
        qaa[k][k] -= R[k]/S;
    }

    /*calculating the Hessian matrix m(d+1) by m(d+1)*/
    /*set very entry of ldl_hessian_mat to be 0;*/
    /*k and l are for vertical (row) indices; kh and lh are for horizontal
     * (column) indices.*/
    for (k = 0; k < m; ++k) {
      for (l = 0; l < (d+1); ++l) {

        for (kh = 0; kh < m; ++kh) {
          for (lh = 0; lh < (d+1); ++lh) {

            ldl_hessian_mat[(k*(d+1) + l)][(kh*(d+1) + lh)] += 
              qaa[k][kh]*H[l]*H[lh];

          }
        }

      }
    }

  }

  UNPROTECT(1);

  /* free arrays */
  free((void *) R);
  free((void *) H);

  free((void *) qaa[0]);
  free((void *) qaa);

}


SEXP logDualLHessianUfWrapper(SEXP n_total, SEXP n_samples, SEXP m, SEXP d,
    SEXP par, SEXP h_func, SEXP env, SEXP x /*input*/)
/* Calculating log dual empirical likelihood (+ n \log n) at a given parameter value.
 * Inputs:
 *   n_total -- total sample size;
 *   n_samples -- a vector of length m+1 specifying the size of each sample;
 *   m -- number of samples - 1;
 *   d -- dimension of h(x); 
 *   par -- values of parameters (length of m(d+1)), organized as \theta_1, \cdots, \theta_m;
 *   h_func -- the basis function of the DRM;
 *   x -- data, organized as x_0, x_1, ..., x_m.
 */
{
  double * restrict n_total_c;
  n_total_c = REAL(n_total);
  double * restrict n_samples_c;
  n_samples_c = REAL(n_samples);
  double * restrict m_c;
  m_c = REAL(m);
  double * restrict d_c;
  d_c = REAL(d);
  double * restrict par_c;
  par_c = REAL(par);
  double * restrict x_c;
  x_c = REAL(x);

  /* Outputs:
   *   ldl_hessian -- a m(d+1) x m(d+1) vector for hessian matrix of ldl at the parameter value par.
   */
  SEXP ldl_hessian;
  PROTECT(ldl_hessian = allocMatrix(REALSXP,
        ((unsigned long)*m_c) * ((unsigned long)*d_c + 1),
        ((unsigned long)*m_c) * ((unsigned long)*d_c + 1)
        ));
  double * restrict ldl_hessian_c;
  ldl_hessian_c = REAL(ldl_hessian);

  /* loop indices */
  unsigned long i;

  double *restrict* restrict par_mat;
  double *restrict* restrict ldl_hessian_mat;

  /* converting x, par and ldl_hessian to matrices */
  par_mat = (double *restrict* restrict) malloc((size_t)
      (((unsigned long)*m_c)*sizeof(double*)));
  if (par_mat == NULL) errMsg("malloc() allocation failure for par_mat!");
  par_mat[0] = par_c;
  for (i = 1; i < (unsigned long)*m_c; ++i){
    par_mat[i] = par_mat[i-1] + ((unsigned long)*d_c + 1);
  }

  ldl_hessian_mat = (double *restrict* restrict) malloc((size_t)
      ((((unsigned long)*m_c) * ((unsigned long)*d_c + 1))*sizeof(double*)));
  if (ldl_hessian_mat == NULL){
    errMsg("malloc() allocation failure for ldl_hessian_mat!");
  }
  ldl_hessian_mat[0] = ldl_hessian_c;
  for (i = 1; i < (((unsigned long)*m_c) * ((unsigned long)*d_c + 1)); ++i){
    ldl_hessian_mat[i] = ldl_hessian_mat[i-1] + 
      (((unsigned long)*m_c) * ((unsigned long)*d_c + 1));
  }

  logDualLHessianUf((unsigned long)*n_total_c, /*inputs*/
      n_samples_c, (unsigned long)*m_c, (unsigned long)*d_c, /*inputs*/
      par_mat, h_func, env, x_c, /*inputs*/
      ldl_hessian_mat /*outputs*/);

  /* free arrays */
  free((void *) par_mat);
  free((void *) ldl_hessian_mat);

  UNPROTECT(1);

  return(ldl_hessian);

}

void probEstUf(unsigned long n_total, /*inputs*/
    double * restrict n_samples, /*inputs*/
    unsigned long m, unsigned long d, /*inputs*/
    double *restrict* restrict par_mat, /*inputs*/
    SEXP h_func, SEXP env, /*inputs*/
    double * restrict x, /*inputs*/
    unsigned short normalize, /*inputs*/
    double *restrict* restrict pEst_mat /*output*/)
/* Estimating \hat{p} for all samples in density ratio models. \hat{p}'s are
 * useful for estimating quantiles, distribution function and density of each
 * population.
 * Inputs:
 *   n_total -- total sample size;
 *   n_samples -- a (double, not unsigned long!) vector of length m+1
 *     specifying the size of each sample;
 *   m -- number of samples - 1;
 *   d -- dimension of h(x); 
 *   par_mat -- values of parameters (pointer matrix of dimension m by (d+1));
 *   h_func -- the basis function of the DRM;
 *   x -- data, a long vector in order of as x_0, x_1, ..., x_m.
 *   normalize -- indicator; wether normalize the probabilies so that the sum
 *     is strictly one: 1 for yes, 0 for no.
 * Outputs:
 *   pEst_mat -- a (m+1) by n_total matrix; estimated probability for x's at
 *     parameter value 'par_mat' for each sample.
 */
{
  /* loop indices */
  unsigned long i, j;

  double * restrict r;
  r = (double * restrict) malloc((size_t) (m*sizeof(double)));
  if (r == NULL) errMsg("malloc() allocation failure for r!");
  for (i = 0; i < m; ++i) {
    r[i] = 0;
  }

  double * restrict h;
  h = (double * restrict) malloc((size_t) (d*sizeof(double)));
  if (h == NULL) errMsg("malloc() allocation failure for h!");
  for (i = 0; i < d; ++i) {
    h[i] = 0;
  }
  SEXP h_args, h_r_call, h_r_call_result;
  PROTECT(h_args = allocVector(REALSXP, 1));

  double p_bl_tmp;

  /* auxiliary variable, used only when normalize == 1 */
  double * restrict pEst_sum;
  pEst_sum = (double * restrict) malloc((size_t) ((m+1)*sizeof(double)));
  if (pEst_sum == NULL) errMsg("malloc() allocation failure for pEst_sum!");
  for (i = 0; i < m+1; ++i) {
    pEst_sum[i] = 0.0;
  }

  /* initializing output variable pEst_mat to 0 as safeguard */
  for (i = 0; i < (m + 1); ++i) {
    for (j = 0; j < n_total; ++j) {
        pEst_mat[i][j] = 0.0;
    }
  }

  for (i = 0; i < n_total; ++i) {

    /*update h*/
    REAL(h_args)[0] = x[i];
    PROTECT(h_r_call = lang2(h_func, h_args));
    PROTECT(h_r_call_result = eval(h_r_call, env));
    for (j = 0; j < d; ++j) {
      h[j] = REAL(h_r_call_result)[j];
    }
    UNPROTECT(2);

    r_val(m, d, h, par_mat, r); /*update r*/

    /* calculate baseline probabilities */
    p_bl_tmp = n_samples[0];
    for (j = 0; j < m; ++j) {
      p_bl_tmp += n_samples[j+1] * r[j];
    }
    p_bl_tmp = 1.0/p_bl_tmp;

    /* calculate probabilities for all m+1 sampels*/
    pEst_mat[0][i] = p_bl_tmp; /*baseline probablities*/
    if (normalize==1) pEst_sum[0] += pEst_mat[0][i];
    for (j = 1; j < m+1; ++j) {
      pEst_mat[j][i] = r[j-1]*p_bl_tmp;
      if (normalize==1) pEst_sum[j] += pEst_mat[j][i];
    }

  }

  if (normalize==1) {

    for (i = 0; i < m+1; ++i) {

      for (j = 0; j < n_total; ++j) {
          pEst_mat[i][j] /= pEst_sum[i];
      }

    }

  }

  UNPROTECT(1);

  /* free arrays */
  free((void *) r);
  free((void *) h);
  free((void *) pEst_sum);

}

SEXP probEstUfWrapper(SEXP n_total, SEXP n_samples, SEXP m, SEXP d,
    SEXP par, SEXP h_func, SEXP env, SEXP x, SEXP normalize /*input*/)
/* Estimating \hat{p} for all samples in density ratio models. \hat{p}'s are
 * useful for estimating quantiles, distribution function and density of each
 * population.
 * Inputs:
 *   n_total -- total sample size;
 *   n_samples -- a vector of length m+1 specifying the size of each sample;
 *   m -- number of samples - 1;
 *   d -- dimension of h(x); 
 *   par_mat -- values of parameters (pointer matrix of dimension m by (d+1));
 *   h_func -- the basis function of the DRM;
 *   x_mat -- 2-D pointer array of data, organized as x_0, x_1, ..., x_m.
 *   normalize -- indicator; wether normalize the probabilies so that the sum
 *     is strictly one: 1 for yes, 0 for no.
 */
{
  double * restrict n_total_c;
  n_total_c = REAL(n_total);
  double * restrict n_samples_c;
  n_samples_c = REAL(n_samples);
  double * restrict m_c;
  m_c = REAL(m);
  double * restrict d_c;
  d_c = REAL(d);
  double * restrict par_c;
  par_c = REAL(par);
  double * restrict x_c;
  x_c = REAL(x);
  double * restrict normalize_c;
  normalize_c = REAL(normalize);

  /* Outputs:
   *   pEst -- a vector of length (m+1)*n_total; estimated probability for x's at parameter value 'par_mat' for each sample.
   */
  SEXP pEst;
  PROTECT(pEst = allocVector(REALSXP, ((unsigned long)*m_c + 1) * (unsigned
          long)*n_total_c));
  double * restrict pEst_c;
  pEst_c = REAL(pEst);

  /* loop indices */
  unsigned long i;

  /* converting x and par to matrices */
  double *restrict* restrict par_mat;
  par_mat = (double *restrict* restrict) malloc((size_t)
      (((unsigned long)*m_c)*sizeof(double*)));
  if (par_mat == NULL) errMsg("malloc() allocation failure for par_mat!");
  par_mat[0] = par_c;
  for (i = 1; i < (unsigned long)*m_c; ++i){
    par_mat[i] = par_mat[i-1] + ((unsigned long)*d_c + 1);
  }

  double *restrict* restrict pEst_mat;
  pEst_mat = (double *restrict* restrict) malloc((size_t)
      (((unsigned long)*m_c+1)*sizeof(double*)));
  if (pEst_mat == NULL) errMsg("malloc() allocation failure for pEst_mat!");
  pEst_mat[0] = pEst_c;
  for (i = 1; i < ((unsigned long)*m_c + 1); ++i){
    pEst_mat[i] = pEst_mat[i-1] + (unsigned long)*n_total_c;
  }

  probEstUf((unsigned long)*n_total_c, /*inputs*/
      n_samples_c, (unsigned long)*m_c, (unsigned long)*d_c, /*inputs*/
      par_mat, h_func, env, x_c, (unsigned short)*normalize_c, /*inputs*/
      pEst_mat /*outputs*/);

  /* free arrays */
  free((void *) par_mat);
  free((void *) pEst_mat);

  UNPROTECT(1);

  return(pEst);

}

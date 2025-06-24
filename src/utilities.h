/* This software is written by Song Cai and published under GPLv3.
 *
 * Version 1.3, April 08, 2014.
 */

/* **********Error messages handling********** */
/*void errMsg(char err_text[]);*/

void errMsg(char err_text[]);
/* ******************** */

/* **********Matrix operations********** */
void kroneckerProd(double *restrict* restrict A, /*inputs*/
    unsigned long m, unsigned long n, /*inputs*/
    double *restrict* restrict B, /*inputs*/
    unsigned long p, unsigned long q, /*inputs*/
    double *restrict* restrict C /*output*/);
/* ******************** */

/* **********C wrappers to LAPACK functions********** */
void dgetrfCWrapper(double * restrict m, /*input*/
    double * restrict n, /*input*/
    double * restrict a, /*input and output*/
    double * restrict lda, /*input*/
    double * restrict ipiv, /*output*/
    double * restrict info /*output*/);
/* ******************** */

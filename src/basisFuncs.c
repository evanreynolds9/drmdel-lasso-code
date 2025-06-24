/* This software is written by Song Cai and published under GPLv3.
 *
 * Version 1.3, April 08, 2014.
 */

#include <math.h>
#include "basisFuncs.h"

/* Basis functions, h(x), for density ratio models */

/* model 1 */
void h1x(double x, /*input*/
    double * restrict h /*output*/)
{
  h[0] = x;
}

/* model 2 */
void h1logx(double x, /*input*/
    double * restrict h /*output*/)
{
  h[0] = log(fabs(x));
}

/* model 3 */
void h1sqrtx(double x, /*input*/
    double * restrict h /*output*/)
{
  h[0] = sqrt(fabs(x));
}

/* model 4 */
void h1xSquare(double x, /*input*/
    double * restrict h /*output*/)
{
  h[0] = x*x;
}

/* model 5 */
void h2Normal(double x, /*input*/
    double * restrict h /*output*/)
/* h(x) function for Normal family */
{
  h[0] = x; h[1] = x*x;
}

/* model 6 */
void h2Gamma(double x, /*input*/
    double * restrict h /*output*/)
/* h(x) function for Gamma family */
{
  h[0] = x; h[1] = log(fabs(x));
}

/* model 7 */
void h3a(double x, /*input*/
    double * restrict h /*output*/)
{
  h[0] = log(fabs(x)); h[1] = sqrt(fabs(x)); h[2] = x;
}

/* model 8 */
void h3b(double x, /*input*/
    double * restrict h /*output*/)
{
  h[0] = log(fabs(x)); h[1] = sqrt(fabs(x)); h[2] = x*x;
}

/* model 9 */
void h3c(double x, /*input*/
    double * restrict h /*output*/)
{
  h[0] = log(fabs(x)); h[1] = x; h[2] = x*x;
}

/* model 10 */
void h3d(double x, /*input*/
    double * restrict h /*output*/)
{
  h[0] = sqrt(fabs(x)); h[1] = x; h[2] = x*x;
}

/* model 11 */
void h4a(double x, /*input*/
    double * restrict h /*output*/)
{
  h[0] = log(fabs(x)); h[1] = sqrt(fabs(x)); h[2] = x; h[3] = x*x;
}

/* model 12 */
void h5a(double x, /*input*/
double * restrict h /*output*/)
{
  h[0] = log(fabs(x)); h[1] = pow(log(fabs(x)), 2); h[2] = sqrt(fabs(x)); h[3] = x; h[4] = x*x;
}

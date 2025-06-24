/* This software is written by Song Cai and published under GPLv3.
 *
 * Version 1.3, April 08, 2014.
 */

/* Basis functions, h(x), for dentity ratio models */

/* model 1 */
void h1x(double x, /*input*/
    double * restrict h /*output*/);

/* model 2 */
void h1logx(double x, /*input*/
    double * restrict h /*output*/);

/* model 3 */
void h1sqrtx(double x, /*input*/
    double * restrict h /*output*/);

/* model 4 */
void h1xSquare(double x, /*input*/
    double * restrict h /*output*/);

/* model 5 */
/* h(x) function for Normal family */
void h2Normal(double x, /*input*/
    double * restrict h /*output*/);

/* model 6 */
/* h(x) function for Gamma family */
void h2Gamma(double x, /*input*/
    double * restrict h /*output*/);

/* model 7 */
void h3a(double x, /*input*/
    double * restrict h /*output*/);

/* model 8 */
void h3b(double x, /*input*/
    double * restrict h /*output*/);

/* model 9 */
void h3c(double x, /*input*/
    double * restrict h /*output*/);

/* model 10 */
void h3d(double x, /*input*/
    double * restrict h /*output*/);

/* model 11 */
void h4a(double x, /*input*/
    double * restrict h /*output*/);

/* model 12 */
void h5a(double x, /*input*/
double * restrict h /*output*/);

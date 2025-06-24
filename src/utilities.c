/* This software is written by Song Cai and published under GPLv3.
 *
 * Version 1.3.2, October 23, 2021.
 */

#include <stdlib.h>
#include <R_ext/Error.h>
#include <R_ext/Lapack.h>
#include "utilities.h"

/* **********Error messages handling********** */
/*void errMsg(char err_text[])*/
/*[> Standard error handler <]*/
/*{*/
    /*fprintf(stderr, "Program run-time error:\n");*/
    /*fprintf(stderr, "    %s\n", err_text);*/
    /*fprintf(stderr, "...now exiting to system...\n");*/
    /*exit(1);*/
/*}*/

void errMsg(char err_text[])
/* Standard error handler; R version */
{
    error(err_text);
}
/* ******************** */


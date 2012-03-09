/*! \file    print.h
    \ingroup optking
    \brief header for print functions
*/

#ifndef _opt_print_h_
#define _opt_print_h_

#include <cstdlib>
#include <cstdio>

namespace opt {

void print_matrix(const FILE *fp, double **A, const int x, const int y);

void print_array(const FILE *fp, double *A, const int x);

void print_geom_array(const FILE *fp, double *A, const int natom);

}

#endif


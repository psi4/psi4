/*###########################################################################*/
/*! \file
    \ingroup EXTREMA
  \brief Small utility functions. */ 
/*						Joseph P. Kenny 11/29/01
  ###########################################################################*/

#include<cstdlib>
#include<cstdio>
#include<cmath>

#include <libciomr/libciomr.h>
#include <psifiles.h> 

namespace psi { namespace extrema {


void stop_io();
void punt(const char *mess);
double **symm_matrix_invert(double **_A, int dim, int print_det, 
			    int redundant); 

/*--------------------------------------------------------------------------*/
/*! \fn gprgid()
  \brief The obligatory gprgid function. */
/*---------------------------------------------------------------------------*/
}



/*--------------------------------------------------------------------------*/
/*! \fn punt(const char *mess)
  \brief Handles unexpected termination. */
/*--------------------------------------------------------------------------*/

void punt(const char *mess)
{
  fprintf(outfile, "  error: %s\n", mess);
  fprintf(stderr, "  EXTREMA error: %s\n", mess);
  stop_io();
  throw PsiException("extrema error", __FILE__, __LINE__);
}



/*---------------------------------------------------------------------------*/
/*! \fn symm_matrix_invert(double **_A, int dim, int print_det, int redundant)
  \brief Inverts a symmetric matrix. 
  \param **_A the matrix to invert
  \param dim dimension of the matrix
  \param print_det 1 to print determinant
  \param redundant 1 if redundant coordinates */
/*                                                      written by Rollin King
----------------------------------------------------------------------------*/

double **symm_matrix_invert(double **_A, int dim, int print_det, int redundant) {
  int i;
  double **_A_inv, **_A_vects, *_A_vals, **_A_temp, det=1.0;

  _A_inv   = init_matrix(dim,dim);
  _A_temp  = init_matrix(dim,dim);
  _A_vects = init_matrix(dim,dim);
  _A_vals  = init_array(dim);

  sq_rsp(dim,dim,_A,_A_vals,1,_A_vects,1E-10);

  if (redundant == 0) {
     for (i=0;i<dim;++i) {
        det *= _A_vals[i];
        _A_inv[i][i] = 1.0/_A_vals[i];
     }
     if (print_det)
        fprintf(outfile,"Determinant: %10.6e\n",det);
     if (fabs(det) < 1E-10) {
        fprintf(outfile,"Determinant: %10.6e\n",det);
        fprintf(outfile,"Determinant is too small...aborting.\n");
        fclose(outfile);
        abort();
     }
  }
  else {
     for (i=0;i<dim;++i) {
        det *= _A_vals[i];
        if (fabs(_A_vals[i]) > 1E-10)
           _A_inv[i][i] = 1.0/_A_vals[i];
     }
     if (print_det)
        fprintf(outfile,"Determinant: %10.6e\n",det);
  }

  mmult(_A_inv,0,_A_vects,1,_A_temp,0,dim,dim,dim,0);
  mmult(_A_vects,0,_A_temp,0,_A_inv,0,dim,dim,dim,0);

  free(_A_vals);
  free_matrix(_A_vects,dim);
  free_matrix(_A_temp,dim);
  return _A_inv;
}

}} // namespace psi::extrema

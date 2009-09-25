#ifndef _psi_bin_extrema_mathtools_h_
#define _psi_bin_extrema_mathtools_h_

namespace psi { namespace extrema {

/*###########################################################################*/
/*! 
  \file
  \ingroup EXTREMA
  \brief Math tools concrete class header file. 
*/

/*! \class math_tools
   \brief A general mathematics routine class.

   General implementations of basic mathematical algorithms.  The only 
   dependency other than standard libraries is the <b>PSI 3.0</b>
   library <b>libciomr</b>. */
/*                                                  Joseph P. Kenny 11/29/01
  ##########################################################################*/

class math_tools{

  protected:
    math_tools() { return; }
    ~math_tools() { return; }

    /* optimization tools */
    double* newton_step(int dim, double **_Hi, double *g);
    double** update_bfgs(int dim, double *_var_dif, 
			 double *grad_dif, double **_Hi_old);
    double** update_ms(int dim, double *_var_dif, 
		       double *_grad_dif, double **_Hi_old);

    /* symmetry tools */
    double** rep_reduce(char* label, double **rep_matrix, int reps);
    double** rep_project(char *label, int dim_vec, double **result_vecs, 
			 int *irrep_proj);
    double** orthogonalize( int nvecs, int dimvecs, double **vecs, 
			    int normalize, double norm_tol, int *nindep);

};

/*! \fn math_tools::math_tools()
  \brief Default constructor, does nothing. */
/*! \fn math_tools::~math_tools()
  \brief Default destructor, does nothing. */



/*##########################################################################*/
/*! \class char_table
  \brief Full character table.

  The character class contains a set of functions for access to point
  group information. */
/*##########################################################################*/

class char_table {

  private:
    char *ptgrp;
    void get_char_table();
    int get_num_irreps();
    const char **get_sym_ops();
    const char **get_irrep_labels();
    int *get_ops_coeffs();
    int get_num_ops();
    int get_num_classes();

  public:
    char_table(char*);
    ~char_table();
    int **ctable;
    int num_irreps;
    const char **sym_ops;
    const char **irrep_labels;
    int *ops_coeffs;
    int num_ops;
    int num_classes;
    
};

}} // namespace psi::extrema

#endif // header guard

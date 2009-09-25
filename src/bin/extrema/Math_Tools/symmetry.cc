/*###########################################################################*/
/*! 
  \file
  \ingroup EXTREMA
  \brief Symmetry related functions. 

  Joseph P. Kenny 12/10/01
  ########################################################################### */

#define EXTERN
#include "extrema.h"

using namespace psi::extrema;

/*---------------------------------------------------------------------------*/
/*! \fn math_tools::rep_reduce(char *label, double **rep_matrix, int num_reps)
  \brief Reduces a set of representation vectors. 
  \param label the point group label
  \param rep_matrix the matrix of representations, each representation 
  should be a row of this matrix whose length = number of irreps
  Parses input */
/*---------------------------------------------------------------------------*/

double **math_tools::rep_reduce(char *label,double **rep_matrix,int num_reps) {
    
    int i,j, vec, rep, pg_order=0;
    char_table ctab(label);
    double **coef_matrix, red_coef;

    coef_matrix = init_matrix(num_reps,ctab.num_irreps);

    for(i=0;i<ctab.num_irreps;++i) 
	pg_order += ctab.ops_coeffs[i];

    int prod;
    for(rep=0;rep<num_reps;++rep) { 
	for(i=0;i<ctab.num_irreps;++i) {
	    red_coef = 0.0;
	    for(j=0;j<ctab.num_irreps;++j) {
		prod = ctab.ctable[i][j] * ctab.ops_coeffs[j];
		red_coef += rep_matrix[rep][j] * (double) prod;
	    }
	    coef_matrix[rep][i] = red_coef/(double)pg_order;
	}
    }

    return coef_matrix;
}

    

/*--------------------------------------------------------------------------*/
/*! \fn math_tools::rep_project(char* label, int dim_vec,
  double **result_vecs, int *irrep_proj)
  \brief Projects out irrep components from a reducible representation.
  \param label point group label
  \param num_vars dimension of the vector
  \param vectors resulting from each symmetry operation
  \param irrep_proj 1 if irrep should be projected out and returned */
/*--------------------------------------------------------------------------*/

double **math_tools::rep_project(char *label, int dim_vec, 
				 double **result_vecs, int *irrep_proj) {
    
    int i,j,ir, proj_num, num_projections, pg_order;
    double **projected;
    char_table c(label);
    
    num_projections=0;
    for(ir=0;ir<c.num_irreps;++ir) 
	if( irrep_proj[ir] )
	    ++num_projections;
    projected = init_matrix(num_projections,dim_vec);

    pg_order=0;
    for(i=0;i<c.num_irreps;++i) 
	pg_order += c.ops_coeffs[i];

    /* project out desired irreps */
    proj_num = 0;
    for(ir=0;ir<c.num_irreps;++ir) 
	if( irrep_proj[ir] ) {
	    for(i=0;i<c.num_irreps;++i) 
		for(j=0;j<dim_vec;++j) {
		    projected[proj_num][j] += c.ctable[ir][i]*
			result_vecs[i][j];
		}
	    for(i=0;i<dim_vec;++i)
		projected[proj_num][i] /= (double) pg_order;
	    ++proj_num;
	}

    return projected;
}



/*------------------------------------------------------------------------*/
/*! \fn math_tools::orthogonalize(
  int nvecs, int dimvecs, double **vecs, int normalize, double norm_tol,
  int *nindep)
  \brief Othogonalizes a set of vectors.
  \param nvecs number of vectors to orthogonalize
  \param dimvecs length of each vector
  \param vec nvec x dimvecs matrix containing vectors
  \param nindep where number of independent vectors is saved
  \param norm_tol norm of vectors big enough to keep
  \param normalize 1 if vectors should be normalized as well
/*------------------------------------------------------------------------*/

double** math_tools::orthogonalize( int nvecs, int dimvecs, double **vecs, 
				    int normalize, double norm_tol,
				    int *nindep) {

    int i,j,k,p, num_big=0;
    double **omat_temp, **omat, norm, dot1, dot2, *temp;

    temp = init_array(dimvecs);

    /* max number of independent vectors is nvecs */
    omat_temp = init_matrix(nvecs,dimvecs);
    
    /* start with first vector */
    for(i=0;i<dimvecs;++i)
	omat_temp[0][i] = vecs[0][i]; 

    /* make each vector orthogonal to all previous vectors */
    for(i=1;i<nvecs;++i) {
	for(k=0;k<dimvecs;++k)
	    omat_temp[i][k] = vecs[i][k]; 
	for(j=0;j<i;++j) {
	    dot1=dot2=0.0;
	    for(k=0;k<dimvecs;++k) {
		dot1 += vecs[i][k]*omat_temp[j][k];
		dot2 += omat_temp[j][k]*omat_temp[j][k];
	    }
	    if( (fabs(dot1)>ALMOST_ZERO) && (dot2>ALMOST_ZERO) ) {
		for(k=0;k<dimvecs;++k) 
		    omat_temp[i][k] -= dot1/dot2 * omat_temp[j][k];
	    }
	}
    }

    /* if norm greater than tolerance, copy over to omat */
    for(i=0;i<nvecs;++i) {
	norm = 0.0;
	for(j=0;j<dimvecs;++j)
	    norm += omat_temp[i][j]*omat_temp[i][j];
	if(norm>norm_tol) 
	    ++num_big;
    }

    omat = init_matrix(num_big,dimvecs);

    p=0;
    for(i=0;i<nvecs;++i) {
	norm = 0.0;
	for(j=0;j<dimvecs;++j)
	    norm += omat_temp[i][j]*omat_temp[i][j];
	if(norm>norm_tol) {
	    for(j=0;j<dimvecs;++j)
		omat[p][j] = omat_temp[i][j];
	    ++p;
	}
    }   

    (*nindep) = num_big;

    if(normalize) 
	for(i=0;i<num_big;++i) {
	    norm=0;
	    for(j=0;j<dimvecs;++j)
		norm += omat[i][j] * omat[i][j];
	    for(j=0;j<dimvecs;++j)
		omat[i][j] /= sqrt(norm);
	}

    return omat;
}

/*###########################################################################*/
/*! 
** \file
** \ingroup EXTREMA
** \brief Member functions for math tools concrete class. */
/*						Joseph P. Kenny 11/29/01
  ###########################################################################*/

#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include"math_tools.h"

using namespace psi::extrema;

/*---------------------------------------------------------------------------*/
/*! \fn math_tools::newton_step(int dim, double **_Hi, double *g)
  \brief Computes the Newton-Raphson optimization step.
  \param dim number of variables
  \param **_Hi inverse hessian matrix
  \param *g gradient vector 
  \return  *_s displacement vector */
/*---------------------------------------------------------------------------*/

double* math_tools :: newton_step(int dim, double **_Hi, double *g) {
    
    int i, j;
    double *_s;
	
    _s = init_array(dim);
    
    for(i=0;i<dim;++i) 
	for(j=0;j<dim;++j) 
	    _s[i] += -1.0 * _Hi[i][j] * g[j];

    return _s;
}



/*---------------------------------------------------------------------------*/
/*! \fn math_tools::update_bfgs(int dim, double *_var_dif, double *_grad_dif,
  double **_Hi_old)
  \brief Performs bfgs update on inverse of hessian matrix.
  \param dim dimension of hessian inverse
  \param *_var_dif difference of current and previous variable
  \param *grad_dif difference of current and previous gradient
  \param **_Hi_old inverse hessian from previous iteration
  \return **_Hi_new updated hessian inverse */
/*---------------------------------------------------------------------------*/

double** math_tools :: update_bfgs(int dim, double *_var_dif, 
				   double *_grad_dif, double **_Hi_old) {

    int i,j;
    double **temp_mat, **_Hi_new, **mat1, **mat2, **mat3, num1=0.0, num2=0.0;
    _Hi_new = init_matrix(dim,dim);
    temp_mat = init_matrix(1,dim);
    mat1 = init_matrix(dim,dim);
    mat2 = init_matrix(dim,dim);
    mat3 = init_matrix(dim,dim);


    mmult(&_grad_dif,0,_Hi_old,0,temp_mat,0,1,dim,dim,0);
    for(i=0;i<dim;++i) 
	num1 += temp_mat[0][i]*_grad_dif[i];
 
    for(i=0;i<dim;++i)
	num2 += _var_dif[i]*_grad_dif[i];

    for(i=0;i<dim;++i)
	for(j=0;j<dim;++j)
	    mat1[i][j] = _var_dif[i]*_var_dif[j];

    mmult(&_grad_dif,0,_Hi_old,0,temp_mat,0,1,dim,dim,0);
    for(i=0;i<dim;++i)
	for(j=0;j<dim;++j)
	    mat2[i][j] = _var_dif[i]*temp_mat[0][j];

    free_matrix(temp_mat,1);
    temp_mat = init_matrix(dim,dim);
    for(i=0;i<dim;++i) 
	for(j=0;j<dim;++j)
	    temp_mat[i][j] = _grad_dif[i]*_var_dif[j];
    mmult(_Hi_old,0,temp_mat,0,mat3,0,dim,dim,dim,0);

    for(i=0;i<dim;++i)
	for(j=0;j<dim;++j)
	    _Hi_new[i][j] = _Hi_old[i][j] + (1 + num1/num2)*(mat1[i][j]/num2)
		- mat2[i][j]/num2 - mat3[i][j]/num2;
    
    free_matrix(temp_mat,dim);
    free_matrix(mat1,dim);
    free_matrix(mat2,dim);
    free_matrix(mat3,dim);

    return _Hi_new;
}



/*---------------------------------------------------------------------------*/
/*! \fn math_tools::update_ms(int dim, double *_var_dif, double *_grad_dif,
  double **_Hi_old)
  \brief Performs update on inverse of hessian matrix attributed to Murtah and
  Sargent among many others.
  \param dim dimension of hessian inverse
  \param *_var_dif difference of current and previous variable
  \param *_grad_dif difference of current and previous gradient
  \param **_Hi_old inverse hessian from previous iteration
  \return **_Hi_new = updated hessian inverse */
/*---------------------------------------------------------------------------*/

double** math_tools :: update_ms(int dim, double *_var_dif, 
				   double *_grad_dif, double **_Hi_old) {

    int i, j;
    double div, *temp_arr, **temp_mat, **_Hi_new;

    /*allocate memory*/
    temp_arr = init_array(dim);
    temp_mat = init_matrix(dim,dim);
    _Hi_new = init_matrix(dim,dim);
    
    for(i=0;i<dim;++i) {
	for(j=0;j<dim;++j) {
	    temp_arr[i] += _Hi_old[i][j] * _grad_dif[j];
	}
    }
    
    for(i=0;i<dim;++i) {
	temp_arr[i] = 0.0;
	for(j=0;j<dim;++j)
	    temp_arr[i] += _Hi_old[i][j] * _grad_dif[j]; 
    }
    for(i=0;i<dim;++i)
	temp_arr[i] = _var_dif[i] - temp_arr[i];
    
    for(i=0;i<dim;++i)
	for(j=0;j<dim;++j)
	    temp_mat[i][j] = temp_arr[i]*temp_arr[j];
    
    div = 0.0;
    for(i=0;i<dim;++i)
	div += temp_arr[i] * _grad_dif[i];
    
    for(i=0;i<dim;++i)
	for(j=0;j<dim;++j)
	    _Hi_new[i][j] = _Hi_old[i][j] + temp_mat[i][j]/div;
    
    /*free up memory*/
    free(temp_arr);
    free_matrix(temp_mat,dim);
    
    return _Hi_new;
}



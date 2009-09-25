/*###########################################################################*/
/*! 
  \file
  \ingroup EXTREMA
  \brief Constructor and member functions for coordinate base class. */
/*						Joseph P. Kenny 11/29/01
  ###########################################################################*/

#define EXTERN
#include "extrema.h"

/*
namespace psi { namespace extrema {
void stop_io();
}}
*/

using namespace psi::extrema;

/*---------------------------------------------------------------------------*/
/*! \fn coord_base::coord_base()
  \brief Coord_base constructor. 

  Parses input */
/*---------------------------------------------------------------------------*/

coord_base :: coord_base() : coord_base_carts(), math_tools() {

    coord_base::parse_input();

    atomic_nums = chkpt_rd_zvals();
    
    return; 
}



/*---------------------------------------------------------------------------*/
/*! \fn coord_base::mem_alloc() 
  \brief Allocates memory.

  Must be called once number of optimized coordinates
  is determined. */
/*---------------------------------------------------------------------------*/

void coord_base :: mem_alloc() {
    
    coords = init_array(num_coords);
    grads = init_array(num_coords);
    Hi = init_matrix(num_coords,num_coords);
    coords_old = init_array(num_coords);
    grads_old = init_array(num_coords);  
    Hi_old = init_matrix(num_coords,num_coords);
    coord_write = init_array(num_coords);
   
    return;
}



/*---------------------------------------------------------------------------*/
/*! \fn coord_base::update_Hi()
  \brief Updates inverse hessian.

  Interface for <b>math_tools</b> update functions. */
/*---------------------------------------------------------------------------*/

void coord_base :: update_Hi() {

    int i;
    double *coord_dif, *grad_dif;
    coord_dif = init_array(num_coords);
    grad_dif = init_array(num_coords);
    for(i=0;i<num_coords;++i) {
	coord_dif[i] = coords[i] - coords_old[i];
	grad_dif[i] = grads[i] - grads_old[i];
    }

    if(print_lvl >= RIDICULOUS_PRINT) {
	fprintf(outfile,"\n  Hi matrix from previous itertation:\n");
	print_mat(Hi_old,num_coords,num_coords,outfile);
	fprintf(outfile,"\n  Coordinate differences:\n");
	for(i=0;i<num_coords;++i)
	    fprintf(outfile,"  %lf - %lf = %lf\n",
		    coords[i],coords_old[i],coord_dif[i]);
	fprintf(outfile,"\n  Gradient differences:\n");
	for(i=0;i<num_coords;++i) 
	    fprintf(outfile,"  %lf - %lf = %lf\n",
		    grads[i],grads_old[i],grad_dif[i]);
    }

    if(!strcmp(update,"MS")) {
	fprintf(outfile,"\n  Performing ms update of inverse hessian\n");
	Hi = math_tools :: update_ms(num_coords, coord_dif, grad_dif, Hi_old);
    }
    else { 
	fprintf(outfile,"\n  Performing bfgs update of inverse hessian\n");
	Hi = math_tools ::update_bfgs(num_coords, coord_dif, grad_dif, Hi_old);
    }

    free(coord_dif);
    free(grad_dif);

    return;
}



/*---------------------------------------------------------------------------*/
/*! \fn coord_base::grad_test()
  \brief Tests for convergence of cartesian gradients. */
/*---------------------------------------------------------------------------*/

void coord_base :: grad_test() {
    
    int i, conv=1;
    double sum=0.0;

    for(i=0;i<3*num_entries;++i) 
	sum += fabs(c_grads[i]);
    
    sum /= (3*num_atoms);

    fprintf(outfile,"\n  RMS gradient = %.10lf\n", sum);

    for(i=0;i<(3*num_entries);++i) 
	if(fabs(c_grads[i]) > (1.0/pow(10.0,(double)grad_max)))
	    conv = 0;
    if(conv) {
	fprintf(outfile,"\n  All gradients below convergence criteria of");
	fprintf(outfile," 10^-%d",grad_max);
	fprintf(outfile,"\n  Optimization completed\n");
	fprintf(stdout,"\n Optimization completed\n");  
        stop_io();
	/* stop_io() calls psi_stop
           psi_stop(infile,outfile,psi_file_prefix) closes infile and outfile */
	exit(PSI_RETURN_ENDLOOP);
    }
    
    return;
}



/*--------------------------------------------------------------------------*/
/*! \fn coord_base::H_test()
  \brief Computes hessian and its eigenvalues. */
/*-------------------------------------------------------------------------*/

void coord_base :: H_test() {

    int i;
    double *evals, **evecs, **H;

    evals = init_array(num_coords);
    evecs = init_matrix(num_coords,num_coords);

    H = symm_matrix_invert(Hi,num_coords,0,0);

    sq_rsp( num_coords, num_coords, H, evals, 1, evecs, 1.0e-14); 
    
    if(print_lvl > NORMAL_PRINT) {
	fprintf(outfile,"\n  H matrix (a.u.):\n");
	print_mat(H, num_coords, num_coords, outfile);

	fprintf(outfile,"\n  H eigenvalues:");
	for(i=0;i<num_coords;++i)
	    fprintf(outfile,"\n  %lf",evals[i]);
	fprintf(outfile,"\n");
    }

    free(evals);
    free_matrix(evecs, num_coords);
    free_matrix(H, num_coords);

    return;
}



/*---------------------------------------------------------------------------*/
/*! \fn coord_base::print_Hi()
  \brief prints inverse of hessian matrix. */
/*---------------------------------------------------------------------------*/

void coord_base :: print_Hi() {

    fprintf(outfile,"\n  Inverse hessian matrix(a.u.):\n");
    print_mat(Hi,num_coords,num_coords,outfile);

    return;
}



/*---------------------------------------------------------------------------*/
/*! \fn coord_base::print_H()
  \brief prints hessian matrix. */
/*---------------------------------------------------------------------------*/

void coord_base :: print_H() {

    double **H_mat;

    H_mat = symm_matrix_invert(Hi, num_coords, 0, 0);

    fprintf(outfile,"\n  Hessian matrix(a.u.):\n");
    print_mat(H_mat,num_coords,num_coords,outfile);

    free_matrix(H_mat,num_coords);

    return;
}







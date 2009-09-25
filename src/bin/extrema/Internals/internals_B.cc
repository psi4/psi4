/*###########################################################################*/
/*! 
   \file
   \ingroup EXTREMA
   \brief Member functions of internals derived class  
   associated with the B matrix. */
/*						Joseph P. Kenny 11/29/01
  ###########################################################################*/

#define EXTERN
#include "extrema.h"
#include "inline.h"

using namespace psi::extrema;

/*---------------------------------------------------------------------------*/
/*! \fn internals::B_row_bond(double *carr, int atom1, int atom2)
   \brief Computes a row of B corrresponding to a simple bonding coordinate.
   \param *carr full cartesian coordinate array
   \param atom1 reference atom 1
   \param atom2 atom bonded to 1 */
/*---------------------------------------------------------------------------*/

double* internals :: B_row_bond(double* c_arr, int atom1, int atom2 ) {

  int i;
  double *row, *unit12;

  row = init_array(num_entries*3);

  unit12 = unit_vec(c_arr, atom1, atom2);

  for(i=0;i<3;++i) {
      row[3*atom1+i] = -unit12[i];
      row[3*atom2+i] =  unit12[i];
    }
  return row;
}



/*--------------------------------------------------------------------------*/
/*! \fn internals::B_row_angle(double *c_arr, int atom1, int atom3, int atom2)
  \brief Computes a row of B corresponding to a simple bending coordinate.
  \param *c_arr full cartesian coordinate array
  \param atom1 reference atom 1
  \param atom3 atom bonded to 1
  \param atom2 atom defining angle 1-3-2 */
/*---------------------------------------------------------------------------*/

double* internals :: B_row_angle(double* c_arr, int atom1, int atom3, 
				 int atom2 ) {

  int i;
  double *row, *unit31, *unit32, r31, r32, cosine, sine;
 
  row = init_array(num_entries*3);

  unit31 = unit_vec(c_arr,atom3, atom1);
  unit32 = unit_vec(c_arr,atom3, atom2);
  cosine = dot_pdt(unit31,unit32);
  sine = sin(acos(cosine));
  r31 = vec_norm(c_arr,atom3,atom1);
  r32 = vec_norm(c_arr,atom3,atom2);

  for(i=0;i<3;++i) {
      row[3*atom1+i] = ( cosine * unit31[i] - unit32[i] ) / ( r31 * sine );
      row[3*atom2+i] = ( cosine * unit32[i] - unit31[i] ) / ( r32 * sine );
      row[3*atom3+i] = ((r31 - r32 * cosine) * unit31[i] + (r32 - r31 * cosine)
			* unit32[i]) / (r31 * r32 * sine); 
    }

  return row;
} 



/*---------------------------------------------------------------------------*/
/*! \fn internals::B_row_tors(double* c_arr,int atom1, int atom2,
                              int atom3, int atom4 )
  \brief Computes row of B matrix corresponding to a simple torsion coordinate.
  \param *c_arr full cartesian array
  \param atom1 reference atom 1
  \param atom2 atom bonded to 1
  \param atom3 atom defining angle 1-2-3
  \param atom4 atom defining torsion 1-2-3-4 */
/*---------------------------------------------------------------------------*/

double* internals :: B_row_tors(double* c_arr,int atom1, int atom2, 
				int atom3, int atom4 ) {

  int i;
  double *row, *unit12, *unit23, *unit43, *unit32, *cross_12_23, *cross_43_32,
      cosine_an2, sine_an2, sine2_an2, 
      cosine_an3, sine_an3, sine2_an3, r12, r23;

  row = init_array(num_entries*3);

  unit12 = unit_vec(c_arr,atom1, atom2);
  unit23 = unit_vec(c_arr,atom2, atom3);
  unit43 = unit_vec(c_arr,atom4, atom3);
  unit32 = unit_vec(c_arr,atom3, atom2);
  cross_12_23 = cross_pdt(unit12, unit23);
  cross_43_32 = cross_pdt(unit43, unit32);
  cosine_an2 = dot_pdt( unit_vec(c_arr, atom2,atom1), 
			unit_vec(c_arr, atom2,atom3));
  sine_an2 = sin(acos( cosine_an2 ));
  sine2_an2 = sine_an2 * sine_an2;
  cosine_an3 = dot_pdt( unit_vec(c_arr, atom3,atom2), 
			unit_vec(c_arr, atom3,atom4));
  sine_an3 = sin(acos( cosine_an3 ));
  sine2_an3 = sine_an3 * sine_an3;
  r12 = vec_norm(c_arr,atom1,atom2);
  r23 = vec_norm(c_arr,atom2,atom3);

  for(i=0;i<3;++i) {
      row[3*atom1+i] = - cross_12_23[i] / ( r12 * sine2_an2 );
      row[3*atom2+i] = (r23 - r12 * cosine_an2) * cross_12_23[i] / 
	  ( r23 * r12 * sine2_an2 ) +
	  cosine_an3 * cross_43_32[i] / ( r23 * sine2_an3 );
    }

  
  /* entries for atom3 same as those for atom2 with permutation of 
     1 with 4 and 2 with 3, entries for atom4 same as those for 
     atom1 with same permutations*/

  /* I permute everything but the angles here */
  unit12 = unit_vec(c_arr,atom4, atom3);
  unit23 = unit_vec(c_arr,atom3, atom2);
  unit43 = unit_vec(c_arr,atom1, atom2);
  unit32 = unit_vec(c_arr,atom2, atom3);
  cross_12_23 = cross_pdt(unit12, unit23);
  cross_43_32 = cross_pdt(unit43, unit32);
  r12 = vec_norm(c_arr,atom4,atom3);
  r23 = vec_norm(c_arr,atom3,atom2);

  /* I permute the angles here */
  for(i=0;i<3;++i) {
      row[3*atom4+i] = - cross_12_23[i] / ( r12 * sine2_an3 );
      row[3*atom3+i] = (r23 - r12 * cosine_an3) * cross_12_23[i] / 
	  ( r23 * r12 * sine2_an3 ) +
	  cosine_an2 * cross_43_32[i] / ( r23 * sine2_an2 );
  }
  
  return row;
} 



/*---------------------------------------------------------------------------*/
/*! \fn internals::compute_G()
  \brief Computes the G=B.B^t matrix. */
/*---------------------------------------------------------------------------*/
  
void internals :: compute_G() {

    double **temp1;
    temp1 = init_matrix(fnum_coords,3*num_entries);

    mmult(B, 0, u, 0, temp1, 0, fnum_coords, 3*num_entries, 3*num_entries, 0);
    mmult(temp1, 0, B, 1, G, 0, fnum_coords, 3*num_entries, fnum_coords, 0);

    free_matrix(temp1,fnum_coords);

    if(print_lvl > RIDICULOUS_PRINT) {
	fprintf(outfile,"\n  G matrix:\n");
	print_mat(G,fnum_coords,fnum_coords,outfile);
    }

    return;
}



/*---------------------------------------------------------------------------*/
/*! \fn internals::compute_A()
  \brief Computes the A=u.B^t.G matrix.
  ---------------------------------------------------------------------------*/
  
void internals :: compute_A() {

    double **temp1, **temp2;
    temp1 = init_matrix(fnum_coords,fnum_coords);
    temp2 = init_matrix(3*num_entries,fnum_coords);

    temp1 = symm_matrix_invert(G, fnum_coords, 0, 1);
    mmult(B,1,temp1,0,temp2,0,3*num_entries,fnum_coords,fnum_coords,0);
    mmult(u,0,temp2,0,A,0,3*num_entries,3*num_entries,fnum_coords,0);

    free_matrix(temp1,fnum_coords);
    free_matrix(temp2,3*num_entries);

    if(print_lvl > RIDICULOUS_PRINT) {
	fprintf(outfile,"\n  A matrix:\n");
	print_mat(A,3*num_entries,fnum_coords,outfile);
    }

    return;
}



/*---------------------------------------------------------------------------*/
/*! \fn internals::print_B()
  \brief Prints the B matrix. */
/*---------------------------------------------------------------------------*/

void internals :: print_B() {

    if(print_lvl > NORMAL_PRINT) {
	fprintf(outfile,"\n  B matrix:\n");
	print_mat(B,fnum_coords,3*num_entries,outfile);
    }
    return;
}



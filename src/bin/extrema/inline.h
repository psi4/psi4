/*###########################################################################*/
/*! \file
    \ingroup EXTREMA
  \brief Inline functions, used to compute simple internal values. */
/*						Joseph P. Kenny 11/29/01
  ###########################################################################*/

#ifndef _psi_bin_extrema_inline_h_
#define _psi_bin_extrema_inline_h_

#include <cmath>

namespace psi { namespace extrema {

inline double *unit_vec(double* cart_arr, int atom1, int atom2 ) {

  int i;
  double *u_vec;

  u_vec = init_array(3);

  for(i=0;i<3;++i) {
      u_vec[i] = (cart_arr[3*atom2+i]-cart_arr[3*atom1+i]) 
	  / vec_norm(cart_arr,atom1,atom2);
    }

  return u_vec;
}



inline double vec_norm(double* cart_arr, int atom1, int atom2 ) {

  int i;
  double norm, temp1, temp2;

  norm = 0;
  for(i=0;i<3;++i) {
      temp1 = cart_arr[3*atom1+i];
      temp2 = cart_arr[3*atom2+i];
      norm += (temp2 - temp1) * (temp2 - temp1);
    }

  norm = sqrt(norm);

  return norm;
}



inline double dot_pdt( double* vec1, double* vec2 ) {

  int i;
  double val;

  val=0;
  for(i=0;i<3;++i) {
      val += vec1[i] * vec2[i];
    }

  return val;
}



/*  double dot_pdt(double *vec1, double *vec2, int num) { */

/*    int i; */
/*    double result=0; */

/*    for(i=0;i<num;++i) { */
/*        result += vec1[i] * vec2[i]; */
/*      } */

/*    return result; */
/*  }      */
  

   
inline double *cross_pdt( double* vec1, double* vec2 ) {

  double* result_vec;

  result_vec = init_array(3);

  result_vec[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
  result_vec[1] = -(vec1[0]*vec2[2] - vec1[2]*vec2[0]);
  result_vec[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];

  return result_vec;
}



inline double norm(double* cart_arr, int atom1, int atom2 ) {

  int i;
  double norm;

  norm = 0;
  for(i=0;i<3;++i) {
      norm += ((cart_arr[3*atom2+i] - cart_arr[3*atom1+i]) * 
	       (cart_arr[3*atom2+i] - cart_arr[3*atom1+i]));
    }

  norm = sqrt(norm);

   return norm;
}



inline double compute_bond(double *car, int atm, int bnd) {
    double l_val;
    l_val = norm(car,atm,bnd);
    return l_val;
}



inline double compute_angle(double *car, int atm, int bnd, int ang) {
    double a_val;
    a_val = acos ( dot_pdt( unit_vec(car, bnd, atm ), 
			    unit_vec(car, bnd, ang ) ) );
    return a_val;
}



inline double compute_torsion(double *car, int atm, int bnd, int ang, int tor){
    
    double t_val, *temp1, *temp2, temp_num;

    temp1 = cross_pdt( unit_vec( car, bnd, ang ),
		       unit_vec( car, ang, tor) ); 
		
    temp2 = cross_pdt( unit_vec( car, atm , bnd), 
		       unit_vec( car, bnd, ang ) );

    temp_num = dot_pdt(temp1,temp2);

    temp_num /= sin( acos( dot_pdt( unit_vec( car, ang, tor ), 
				    unit_vec( car, ang, bnd ) )))
	* sin( acos( dot_pdt( unit_vec( car, bnd, ang ),
			      unit_vec( car, bnd, atm ) )));

    if(temp_num>(1-ALMOST_ONE))
	t_val = 0.0;
    else if(temp_num<(-1+ALMOST_ONE))
	t_val = _pi;
    else 
	t_val = acos(temp_num);
    
    free(temp1);
    free(temp2);

    return t_val;
}

}} // namespace psi::extrema

#endif

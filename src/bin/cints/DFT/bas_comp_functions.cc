/*! \file bas_comp_functions.cc
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>


namespace psi { namespace cints {

  double calc_exp_basis(int shell_num, double rr){
    
    int i;
    int shell_type;
    int shell_start;
    int shell_end;
    double expon,coeff;
    double bastmp;
    
    shell_type = BasisSet.shells[shell_num].am;
    shell_start = BasisSet.shells[shell_num].fprim-1;
    shell_end = shell_start+BasisSet.shells[shell_num].n_prims;
    
    bastmp = 0.0;
    for(i=shell_start;i<shell_end;i++){
      expon = -BasisSet.cgtos[i].exp;
      coeff = BasisSet.cgtos[i].ccoeff[shell_type-1];
      bastmp += coeff*exp(expon*rr);
    }
    
    return bastmp;
  }

  double calc_radial_bas(int shell_num, double rr, double r){
    
    int i;
    int shell_type;
    double bastmp;
    
    shell_type = BasisSet.shells[BasisSet.am2shell[shell_num]].am-1;
    bastmp = calc_exp_basis(shell_num,rr);
    for(i=0;i<shell_type;i++)
      bastmp *= r;
    
    return bastmp;
  }
  
}

}

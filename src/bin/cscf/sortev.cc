/*! \file 
    \ingroup (CSCF)
    \brief Enter brief description of file here 
*/
/*-------------------------------------------------------------

  sortev.c-a routine for sorting the energies and symmetries of
  of orbitals into two arrays containing the symmetries of the
  orbitals and the other containing the energies

  ----------------------------------------------------------*/
#define EXTERN
#include "includes.h"
#include "common.h"

namespace psi { namespace cscf {

void sortev()
{
  int i,j,k;
  static int *ct;
  static int *pt;
  int rep;
  double lowvalue;
  int elect;
  struct symm *s;
  int num_arr,flag;
  int num_mo = 0;

  ct = (int *) init_int_array(num_ir);/*counter map array*/
  pt = (int *) init_int_array(num_ir);/*pointer map array*/
  elect = 1;
  num_arr = num_ir;
 
  for(i=0; i < num_arr;i++) {
    ct[i]=0;
    pt[i]=i;
    num_mo += scf_info[i].num_mo;
  }

  /* if this function is used anytime besides the first iteration
     the eigenvalues must be converted to the arrays scf_info.hevals*/
  
  if(iter != 0 && iopen !=1 || converged)
    for(i=0; i < num_ir; i++)
      for(j=0; j < scf_info[i].num_mo; j++)
	scf_info[i].hevals[j] = scf_info[i].fock_evals[j];

/* excluding irreps that have zero so's out of consideration */
  for(i=0; i < num_arr; i++)
    if(ct[pt[i]]==scf_info[pt[i]].num_mo) {
      if(i != num_arr-1) {
	for(j=i;j<num_arr-1;j++)
	  pt[j]=pt[j+1];
	num_arr--;
	i--;
      }
      else
	num_arr--;
    }

/* Sort eigenvalues and put them in order in the ener_tot array*/
  
  do{
      lowvalue = scf_info[pt[0]].hevals[ct[pt[0]]];
      rep = pt[0];
      flag = 0;
      
      for(j = 1;j < num_arr; j++)
	  if(lowvalue > scf_info[pt[j]].hevals[ct[pt[j]]])
	      {
		lowvalue = scf_info[pt[j]].hevals[ct[pt[j]]];
		rep = pt[j];
		flag = j;
	      }
      
      /* Checking if have any SOs left in irrep rep */
      ct[rep]++;
      if(ct[rep]==scf_info[rep].num_mo) {
	  if(flag != num_arr-1) {
	    for(j=flag;j<num_arr-1;j++)
	      pt[j]=pt[j+1];
	    num_arr--;
	  }
	  else
	    num_arr--;
      }
      ener_tot[elect-1]=lowvalue;
      symm_tot[elect-1]=rep;
      elect = elect + 1;
      
  } while (elect != num_mo + 1);

}

}} // namespace psi::cscf

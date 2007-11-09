/*! \file 
    \ingroup (CSCF)
    \brief Enter brief description of file here 
*/
/* ------------------------------------

   cmatsplit - takes cmat and splits
   the scf_info struct into two spin
   struct cmats.  Also intend to add
   a mixing here to force convergence
   to UHF solutions.

   -------------------------------------*/

#define EXTERN
#include "includes.h"
#include "common.h"

namespace psi { namespace cscf {

void cmatsplit(void){
    int i,j,k,l;
    int nn, nmo;
    double temp;
    
    for(i=0;i<num_ir;i++){
	nn=scf_info[i].num_so;
        nmo=scf_info[i].num_mo;
	for(j=0;j<nn;j++){
	    for(k=0;k<nmo;k++){
		temp = scf_info[i].cmat[j][k];
		spin_info[0].scf_spin[i].cmat[j][k] = temp;
		spin_info[1].scf_spin[i].cmat[j][k] = temp;
	    }
	}
	/*fprintf(outfile,"\n\nCmat for spin = %d and irrep = %s\n"
		,0,scf_info[i].irrep_label);
	print_mat(spin_info[0].scf_spin[i].cmat,nn,nmo,outfile);
	fprintf(outfile,"\n\nCmat for spin = %d and irrep = %s\n"
		,1,scf_info[i].irrep_label);
		print_mat(spin_info[1].scf_spin[i].cmat,nn,nmo,outfile);*/
    }
}
		
}} // namespace psi::cscf

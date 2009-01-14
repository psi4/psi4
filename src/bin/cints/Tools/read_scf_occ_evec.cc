/*! \file read_scf_occ_evec.cc
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libint/libint.h>

#include "defines.h"
#define EXTERN
#include "global.h"

namespace psi { namespace CINTS {

void read_scf_occ_evec(void)
{ 
    int i, j, k, l, jj, ij;
    int nstri;
    int nstria,nstrib;
    int num_so,num_ao,num_mo,ndocc;
    int aocc,bocc;
    int bas_off;
    int shell_start,shell_end,shell_type;
    double **SO_cmat, **SO_cmato;
    double **Cocc_un;
    
    /* Read All information necessary to do DFT procedure with Eigenvector */
    
    psio_open(IOUnits.itapDSCF, PSIO_OPEN_OLD);
    
    psio_read_entry(IOUnits.itapDSCF, "Number of MOs", 
		    (char *) &(MOInfo.num_mo), sizeof(int));
    if(UserOptions.reftype == uhf){
	psio_read_entry(IOUnits.itapDSCF, "Number of Alpha DOCC", 
		    (char *) &(MOInfo.alpha_occ),sizeof(int));
	psio_read_entry(IOUnits.itapDSCF, "Number of Beta DOCC", 
		    (char *) &(MOInfo.beta_occ),sizeof(int));
    }
    else{
	psio_read_entry(IOUnits.itapDSCF, "Number of DOCC", 
			(char *) &(MOInfo.ndocc),sizeof(int));
    }
    
/*-------------------
  Variables needed
  -------------------*/
    
    num_ao = BasisSet.num_ao;
    num_so = Symmetry.num_so;
    num_mo = MOInfo.num_mo;
    ndocc = MOInfo.ndocc;
    nstri = num_mo*ndocc;
    aocc = MOInfo.alpha_occ;
    bocc = MOInfo.beta_occ;
    nstria = num_mo*aocc;
    nstrib = num_mo*bocc;
    
    
    if (UserOptions.reftype == uhf){
	  SO_cmat = block_matrix(num_so,aocc);
	  SO_cmato = block_matrix(num_so,bocc);
   
	  psio_read_entry(IOUnits.itapDSCF, "Alpha Occupied SCF Eigenvector", 
		      (char *) &(SO_cmat[0][0]), sizeof(double)*nstria);
	  psio_read_entry(IOUnits.itapDSCF, "Beta Occupied SCF Eigenvector", 
			  (char *) &(SO_cmato[0][0]), sizeof(double)*nstrib);
      }
      else {
	  SO_cmat = block_matrix(num_so,ndocc);
	  psio_read_entry(IOUnits.itapDSCF, "Occupied SCF Eigenvector", 
			  (char *) &(SO_cmat[0][0]), sizeof(double)*nstri);
      }
      psio_close(IOUnits.itapDSCF, 1);
      
  
   /*----------------------
    transform to AO basis
    ----------------------*/
      /*fprintf(outfile,"\nUSOTAO matrix");
	print_mat(Symmetry.usotao,num_so,num_ao,outfile);
	fprintf(outfile,"\nSO Cmat");
	print_mat(SO_cmat,num_so,num_mo,outfile);
      */
      if(UserOptions.reftype == uhf){
	  Cocca = block_matrix(num_ao,aocc);
	  mmult(Symmetry.usotao,1,SO_cmat,0,Cocca,0,num_ao,num_so,aocc,0);
	  free_block(SO_cmat);
	  Coccb = block_matrix(num_ao,bocc);
	  mmult(Symmetry.usotao,1,SO_cmato,0,Coccb,0,num_ao,num_so,bocc,0);
	  free_block(SO_cmato);
      }
      else{
	  Cocc = block_matrix(num_ao,ndocc);
	  mmult(Symmetry.usotao,1,SO_cmat,0,Cocc,0,num_ao,num_so,ndocc,0);
	  free_block(SO_cmat);
      }
      
      /*--------------------------
	Remove after done testing
	--------------------------*/   
      
      /*  fprintf(outfile,"\nCocca");
      print_mat(Cocca,num_ao,aocc,outfile);
      fprintf(outfile,"\nCoccb");
      print_mat(Coccb,num_ao,bocc,outfile);*/
      
      return;
}


}}

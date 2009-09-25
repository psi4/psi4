#define MAX(i,j) ((i>j) ? i : j)
#define MIN(i,j) ((i>j) ? j : i)
#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))

#include <iostream>

#include <psifiles.h>
#include <libmoinfo/libmoinfo.h>
#include <liboptions/liboptions.h>
#include <libiwl/iwl.h>
#include <libutil/libutil.h>

#include "scf.h"

extern FILE* outfile;

namespace psi{ namespace mcscf{

void SCF::read_so_oei()
{
  // Read all the SO one electron integrals in Pitzer order
  double* buffer = new double[nso*(nso+1)/2];

  // Read the kinetic energy integrals
  for(int k=0; k<nso*(nso+1)/2;++k) buffer[k] = 0.0;
  iwl_rdone(PSIF_OEI,const_cast<char*>(PSIF_SO_T),buffer,nso*(nso+1)/2,0,0,outfile);

  for(int h=0;h<nirreps;h++){
    for(int i=0; i < H->get_rows(h); i++){
      for(int j=0; j < H->get_cols(h); j++){
        int ij = INDEX(H->get_abs_row(h,i),H->get_abs_col(h,j));
        H->set(h,i,j,buffer[ij]);
      }
    }
  }
  // Read the potential energy integrals
  for(int k=0; k<nso*(nso+1)/2;++k) buffer[k] = 0.0;
  iwl_rdone(PSIF_OEI,const_cast<char*>(PSIF_SO_V),buffer,nso*(nso+1)/2,0,0,outfile);

  for(int h=0;h<nirreps;h++){
    for(int i=0; i < H->get_rows(h); i++){
      for(int j=0; j < H->get_cols(h); j++){
        int ij = INDEX(H->get_abs_row(h,i),H->get_abs_col(h,j));
        H->add(h,i,j,buffer[ij]);
      }
    }
  }

  // Read the overlap integrals
  for(int k=0; k<nso*(nso+1)/2;++k) buffer[k] = 0.0;
  iwl_rdone(PSIF_OEI,const_cast<char*>(PSIF_SO_S),buffer,nso*(nso+1)/2,0,0,outfile);

  for(int h=0;h<nirreps;h++){
    for(int i=0; i < S->get_rows(h); i++){
      for(int j=0; j < S->get_rows(h); j++){
        int ij = INDEX(H->get_abs_row(h,i),H->get_abs_col(h,j));
        S->set(h,i,j, buffer[ij] );
      }
    }
  }

  delete[] buffer;

  if(options_get_int("DEBUG") > 4){
    S->print();
    H->print();
  }
}

}} /* End Namespaces */


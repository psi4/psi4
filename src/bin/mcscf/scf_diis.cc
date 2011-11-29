#include <iostream>
#include <cstdio>

#include <liboptions/liboptions.h>
#include <libutil/libutil.h>

#include "scf.h"
#include "algebra_interface.h"

extern FILE* outfile;

namespace psi{ namespace mcscf{

extern MemoryManager* memory_manager;

void SCF::diis(int cycle)
{
  // Transform Feff from MO to orthogonal AO
  SBlockMatrix CFeffC("CFeffC",nirreps,sopi,sopi);

  transform(Feff_t,CFeffC,C_T);
  transform(CFeffC,Feff_oAO,S_sqrt);

  diis_F[current_diis] = Feff_oAO;

  // Build and transform the error vector from MO to orthogonal AO
  SBlockMatrix CeC("CeC",nirreps,sopi,sopi);
  SBlockMatrix e_oAO("e_oAO",nirreps,sopi,sopi);

  e = Feff_t;
  e->zero_diagonal();
  transform(e,CeC,C_T);
  transform(CeC,e_oAO,S_sqrt);

  diis_e[current_diis] = e_oAO;

  if(reference == tcscf){
    for(int I = 0 ; I < nci; ++I){
      diis_ci[I][current_diis] = ci[I];
    }
  }

  fprintf(outfile," S");

  if(cycle >= ndiis){
    fprintf(outfile,"/E");
    int  matrix_size = ndiis + 1;

    double** diis_B;
    double*  diis_A;
    allocate1(double,diis_A,matrix_size);
    allocate2(double,diis_B,matrix_size,matrix_size);

    // Zero A and B
    for(int i=0;i<ndiis;i++){
      diis_A[i]=0.0;
      diis_B[i][ndiis]=diis_B[ndiis][i]=-1.0;
      for(int j=0;j<ndiis;j++)
        diis_B[i][j]=0.0;
    }
    diis_B[ndiis][ndiis]=0.0;
    diis_A[ndiis]=-1.0;

    // Build the diis_B matrix
    for(int i=0;i<ndiis;i++){
      for(int j=i;j<ndiis;j++){
        diis_B[i][j] += dot(diis_e[i],diis_e[j]);
        diis_B[j][i] = diis_B[i][j];
      }
    }

    // Solve B x = A
    int* IPIV = new int[matrix_size];
    int nrhs = 1;
    int info = 0;
    F_DGESV(&matrix_size, &nrhs, &(diis_B[0][0]),&matrix_size, &(IPIV[0]), &(diis_A[0]),&matrix_size, &info);
    delete[] IPIV;

    // Update F = sum diis_F(i) * A(i);
    if(!info){
      Feff_oAO->zero();
      for(int i=0;i<ndiis;i++){
        e = diis_F[i];
        e->scale(diis_A[i]);
        Feff_oAO += e;
      }

      if(reference == tcscf && options_.get_bool("CI_DIIS") ){
        for(int I = 0 ; I < nci; ++I){
          ci[I] = 0.0;
          for(int i=0; i < ndiis;i++){
            ci[I] += diis_ci[I][i] * diis_A[i];
          }
        }
      }

    }else{
      fprintf(outfile," (singularities found)");
    }

    release1(diis_A);
    release2(diis_B);
  }
  current_diis++;
  if(current_diis == ndiis)
    current_diis  = 0;
}

}} /* End Namespaces */

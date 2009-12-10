/*! \file
    \ingroup LMP2
    \brief localized the SCF MO's
*/
#include "mpi.h"
//#include <iostream>
//#include <fstream>              // file I/O support
//#include <cstdio>
//#include <cstdlib>
//#include <cstring>
//#include <cmath>
//#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
//#include <libchkpt/chkpt.h>
//#include <libpsio/psio.h>
//#include <libqt/qt.h>
//#include <libint/libint.h>
#include <libmints/basisset.h>
//#include <libmints/onebody.h>
#include <libmints/twobody.h>
#include <libmints/integral.h>
#include <libmints/factory.h>
#include <libmints/symmetry.h>
#include <libmints/wavefunction.h>
//#include <psifiles.h>
#define EXTERN
#include "globals.h"

namespace psi{

extern int myid;
extern int nprocs;

namespace lmp2{

extern int myid_lmp2;
extern int nprocs_lmp2;


void LMP2::direct_transformation() {

  using namespace std;

  int i, j, a, b, ab, t, ij;
  int k, v;
  int count;
  int *ij_owner, *ij_local, *mn_owner;
  int **MN_shell;
  int M, R, N, S;
  double ****eri_1, ***eri_2, ***eri_2_mn, ***eri_3, ***eri_4;
  MPI_Status stat;

  //  ****  These are required to utilize libmints  ****

  // Required for libmints, allocates and computes the following:
  // ioff, fac, df, bc
  Wavefunction::initialize_singletons();

  // Create a basis set object and initialize it using the checkpoint file.
  shared_ptr<BasisSet> basis(new BasisSet(chkpt));

  // Initialize an integral factory
  shared_ptr<IntegralFactory> integral(new IntegralFactory(basis, basis, basis, basis));

  // Create an integral object for ERIs
  shared_ptr<TwoBodyInt> eri(integral->eri());

  // Get the storage buffer from the eri object
  const double *buffer = eri->buffer();

  // **** Integral direct transformation ****

  //  **** Allocate memory for the first quarter integral transformation ****
//  double *B1 = (double *) calloc(nso*nocc*nso*nocc , sizeof(double));

  ij_owner = get_ij_owner();
  ij_local = get_ij_local();

  if(ij_pairs%nprocs == 0) {
    eri_2 = (double ***) malloc((ij_pairs/nprocs) * sizeof(double **));
    for(ij=0; ij < ij_pairs; ij++) {
      if(myid == ij_owner[ij])
        eri_2[ij_local[ij]] = block_matrix(nso, nso);
    }
  }
  else {
    if(myid < ij_pairs%nprocs) {
      eri_2 = (double ***) malloc(((ij_pairs/nprocs) + 1) * sizeof(double **));
      for(ij=0; ij < ij_pairs; ij++) {
        if(myid == ij_owner[ij])
          eri_2[ij_local[ij]] = block_matrix(nso, nso);
      }
    }
    else {
      eri_2 = (double ***) malloc(ij_pairs/nprocs * sizeof(double **));
      for(ij=0; ij < ij_pairs; ij++) {
        if(myid == ij_owner[ij])
          eri_2[ij_local[ij]] = block_matrix(nso, nso);
      }
    }
  }


  num_unique_shells=0;
  for (M=0; M < nshell; M++) {
    for (N=0; N < nshell; N++) {
      num_unique_shells++;
    }
  }

  MN_shell = init_int_matrix(4, num_unique_shells);

  count = 0;
  for (M=0; M < nshell; M++) {
    int numm = basis->shell(M)->nfunction();
    for (N=0; N < nshell; N++, count++) {
      int numn = basis->shell(N)->nfunction();
      MN_shell[0][count] = M;
      MN_shell[1][count] = numm;
      MN_shell[2][count] = N;
      MN_shell[3][count] = numn;
    }
  }

  sort_shell(MN_shell, num_unique_shells);

  eri_1 = (double ****)  calloc(MN_shell[1][0] , sizeof(double ***));
  for(i=0; i < MN_shell[1][0]; i++) {
    eri_1[i] = (double ***)  calloc(nso , sizeof(double **));
    for(j=0; j < nso; j++) { 
      eri_1[i][j] = block_matrix(MN_shell[1][0],nocc);
    } 
  } 

  eri_2_mn = (double ***)  calloc(ij_pairs , sizeof(double **));
  for(ij=0; ij < ij_pairs; ij++) {
    eri_2_mn[ij] = block_matrix(MN_shell[1][0],MN_shell[1][0]);
  }   

  double **rec;
  rec = block_matrix(MN_shell[1][0],MN_shell[1][0]);

  mn_owner = get_mn_owner(num_unique_shells);

  //  **** Integral direct half AO-MO integral transformation
  v=0;
  for(count=0; count < num_unique_shells; count++) { 
    M = MN_shell[0][count];
    int numm = MN_shell[1][count];
    N = MN_shell[2][count];
    int numn = MN_shell[3][count];

      if(myid == mn_owner[count]) {

        for(i=0; i < MN_shell[1][0]; i++) {
          for(j=0; j < nso; j++) {
            zero_mat(eri_1[i][j], MN_shell[1][0], nocc);
          }
        }

        for (R=0; R < nshell; R++) {
          int numr = basis->shell(R)->nfunction(); 
          for (S=0; S < nshell; S++) {
            int nums = basis->shell(S)->nfunction();

            eri->compute_shell(M, R, N, S);

            for(j=0; j < nocc; j++) {
              //  **** First quarter integral transformation ****
              int index = 0;
              for(int m=0; m < numm; m++) {
                for(int r=0; r < numr; r++) {
                  int oor = basis->shell(R)->function_index()+r;
                  for(int n=0; n < numn; n++) {
                    for(int s=0; s < nums; s++, index++) {
                      int os = basis->shell(S)->function_index()+s;

                      if (fabs(buffer[index]) > 1.0e-14) 
                        eri_1[m][oor][n][j] += C[os][j] * buffer[index];

                    } // End of s loop
                  } // End of n loop
                } // End of r loop
              } // End of m loop
            } // End of j loop
          } // End of S loop
        } // Enf of N loop

        for(ij=0; ij < ij_pairs; ij++) {
          zero_mat(eri_2_mn[ij], MN_shell[1][0], MN_shell[1][0]);
        }

        //  **** Second quarter integral transformation ****
        for(i=0, ij=0; i < nocc; i++) {
          for(j=0; j <= i; j++, ij++) {
            for(int m=0; m < numm; m++) {
              for(int n=0; n < numn; n++) {
                for (R=0; R < nshell; R++) {
                  int numr = basis->shell(R)->nfunction();
                  for(int r=0; r < numr; r++) {
                    int oor = basis->shell(R)->function_index()+r;

                    if (fabs(eri_1[m][oor][n][j]) > 1.0e-14) 
                      eri_2_mn[ij][m][n] += C[oor][i] * eri_1[m][oor][n][j];

                  } // End of n loop
                } // End of N loop
              } // End of r loop
            } // End of m loop
          } // End of j loop
        } // End of i loop
      } // End if myid loop

      for(ij=0; ij < ij_pairs; ij++) {
        zero_mat(rec, MN_shell[1][0], MN_shell[1][0]);
        if ((myid == ij_owner[ij]) && (myid == mn_owner[count])) {
          for(int m=0; m < numm; m++) {
            int om = basis->shell(M)->function_index()+m;
            for(int n=0; n < numn; n++) {
              int oon = basis->shell(N)->function_index()+n;

              if (fabs(eri_2_mn[ij][m][n]) > 1.0e-14) 
                eri_2[ij_local[ij]][om][oon] = eri_2_mn[ij][m][n];

            }
          }
        }
        else if ((myid == ij_owner[ij]) && (myid != mn_owner[count])) {
          MPI::COMM_WORLD.Recv(&(rec[0][0]), (MN_shell[1][0])*(MN_shell[1][0]), MPI::DOUBLE, mn_owner[count], ij);
          for(int m=0; m < numm; m++) {
            int om = basis->shell(M)->function_index()+m;
            for(int n=0; n < numn; n++) {
              int oon = basis->shell(N)->function_index()+n;

              if (fabs(rec[m][n]) > 1.0e-14)
                eri_2[ij_local[ij]][om][oon] = rec[m][n];

            }
          }
        }
        else if ((myid != ij_owner[ij]) && (myid == mn_owner[count]))
          MPI::COMM_WORLD.Send(&(eri_2_mn[ij][0][0]), (MN_shell[1][0])*(MN_shell[1][0]), MPI::DOUBLE, ij_owner[ij], ij);
      }

      v++;
  } // End of M loop

  //  **** Free the memory used for the first quarter integral transformation ****
  for(i=0; i < MN_shell[1][0]; i++) {
    for(j=0; j < nso; j++) {
      free_block(eri_1[i][j]);
    }
    free(eri_1[i]);
  }
  free(eri_1);


  if(ij_pairs%nprocs == 0) {
    eri_3 = (double ***) malloc((ij_pairs/nprocs) * sizeof(double **));
    for(ij=0; ij < ij_pairs; ij++) {
      if(myid == ij_owner[ij])
        eri_3[ij_local[ij]] = block_matrix(nso, pairdom_len[ij]);
    }
  }
  else {
    if(myid < ij_pairs%nprocs) {
      eri_3 = (double ***) malloc(((ij_pairs/nprocs) + 1) * sizeof(double **));
      for(ij=0; ij < ij_pairs; ij++) {
        if(myid == ij_owner[ij])
          eri_3[ij_local[ij]] = block_matrix(nso, pairdom_len[ij]);
      }
    }
    else {
      eri_3 = (double ***) malloc(ij_pairs/nprocs * sizeof(double **));
      for(ij=0; ij < ij_pairs; ij++) {
        if(myid == ij_owner[ij])
          eri_3[ij_local[ij]] = block_matrix(nso, pairdom_len[ij]);
      }
    }
  }



  //  **** Third quarter integral transformation ****
  v=0;
  for(i=0, ij=0; i < nocc; i++) {
    for(j=0; j <= i; j++, ij++) {
      if(v%nprocs == myid) {
        for(k=0, b=0; k < natom; k++) {
          if(pairdomain[ij][k]) {
            for(t=aostart[k]; t <= aostop[k]; t++, b++) {
              for (M=0; M < nshell; M++) {
                int numm = basis->shell(M)->nfunction();
                for(int m=0; m < numm; m++) {
                  int om = basis->shell(M)->function_index()+m;
                  for (N=0; N < nshell; N++) {
                    int numn = basis->shell(N)->nfunction();
                    for(int n=0; n < numn; n++) {
                      int oon = basis->shell(N)->function_index()+n;

                      if (fabs(eri_2[ij_local[ij]][om][oon]) > 1.0e-14) {

                        eri_3[ij_local[ij]][om][b] += Rt_full[oon][t] * eri_2[ij_local[ij]][om][oon];

                      } 
                    } // End of n loop
                  } // End of N loop
                } // End of m loop
              } // End of M loop
            } // End of b loop
          } // End of if pairdomain
        } // End of k loop 
      } // End of if myid loop
      v++;
    } // End of j loop
  } // End of i loop


  if(ij_pairs%nprocs == 0) {
    Ktilde = (double ***) malloc((ij_pairs/nprocs) * sizeof(double **));
    for(ij=0; ij < ij_pairs; ij++) {
      if(myid == ij_owner[ij])
        Ktilde[ij_local[ij]] = block_matrix(pairdom_len[ij], pairdom_len[ij]);
    }
  }
  else {
    if(myid < ij_pairs%nprocs) {
      Ktilde = (double ***) malloc(((ij_pairs/nprocs) + 1) * sizeof(double **));
      for(ij=0; ij < ij_pairs; ij++) {
        if(myid == ij_owner[ij])
          Ktilde[ij_local[ij]] = block_matrix(pairdom_len[ij], pairdom_len[ij]);
      }
    }
    else {
      Ktilde = (double ***) malloc(ij_pairs/nprocs * sizeof(double **));
      for(ij=0; ij < ij_pairs; ij++) {
        if(myid == ij_owner[ij])
          Ktilde[ij_local[ij]] = block_matrix(pairdom_len[ij], pairdom_len[ij]);
      }
    }
  }

  //  **** Fourth quarter integral transformation ****
  v=0;
  for(i=0, ij=0; i < nocc; i++) {
    for(j=0; j <= i; j++, ij++) {
      if(v%nprocs == myid) {
        for(k=0, a=0; k < natom; k++) {
          if(pairdomain[ij][k]) {
            for(t=aostart[k]; t <= aostop[k]; t++, a++) {
              for(b=0; b < pairdom_len[ij]; b++) {
//                Ktilde[ij_local[ij]][a][b] = 0;
                for (M=0; M < nshell; M++) {
                  int numm = basis->shell(M)->nfunction();
                  for(int m=0; m < numm; m++) {
                    int om = basis->shell(M)->function_index()+m;
        
                    if (fabs(eri_3[ij_local[ij]][om][b]) > 1.0e-14) {
//                      Ktilde[ij][a][b] += Rt_full[om][t] * eri_3[ij][om][b];
                      Ktilde[ij_local[ij]][a][b] += Rt_full[om][t] * eri_3[ij_local[ij]][om][b];
                    }

                  } // End of m loop
                } // End of M loop
              } // End of b loop
            } // End of a loop
          } // End of if pairdomain
        } // End of k loop
      } // End of if myid loop
      v++;
    } // End of j loop
  } // End of i loop

  if(print > 2) {
    v=0;
    for(i=0, ij=0; i < nocc; i++) {
      for(j=0; j <= i; j++, ij++) {
        if(v%nprocs == myid) {
          fprintf(outfile, "Ktilde[%d] matrix", ij);
          print_mat(Ktilde[ij_local[ij]], pairdom_len[ij], pairdom_len[ij], outfile);
        }
        v++;
      }
    }
  }


  free(MN_shell);

}

}} // namespace psi::lmp2

/*! \file
    \ingroup LMP2
    \brief Construct MO integrals from density-fitted 3-center quantities
 */

#include <psi4-dec.h>
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <libipv1/ip_lib.h>
#include <libiwl/iwl.h>
#include <libqt/qt.h>

#include <iostream>
#include <fstream>              // file I/O support
#include <libciomr/libciomr.h>
#include <libmints/basisset.h>
#include <libmints/onebody.h>
#include <libmints/twobody.h>
#include <libmints/integral.h>
#include <libmints/factory.h>
//#include <libmints/symmetry.h>
#include <libmints/wavefunction.h>
#include <libparallel/parallel.h>
//#include <psifiles.h>
#define EXTERN
#include "globals.h"

#define TIME_DF_LMP2 1

namespace psi {

namespace lmp2 {

void LMP2::direct_df_transformation2() {

  using namespace std;

  // Required for libmints, allocates and computes:
  // ioff, fac, df, bc
  Wavefunction::initialize_singletons();

  // Create a basis set object and initialize it using the checkpoint file.
  shared_ptr<BasisSet>basis = shared_ptr<BasisSet > (new BasisSet(chkpt));
  shared_ptr<BasisSet>ribasis = shared_ptr<BasisSet > (new BasisSet(chkpt, "DF_BASIS"));
  shared_ptr<BasisSet>zero = BasisSet::zero_basis_set();

  // Create integral factory
  IntegralFactory rifactory(ribasis, zero, basis, basis);
  IntegralFactory rifactory_J(ribasis, zero, ribasis, zero);

  // Create an integral object for ERIs
  
  TwoBodyInt* eri = rifactory.eri();
  TwoBodyInt* Jint = rifactory_J.eri();

  int i, j, ij, k, l, m, n, a, b, s, t, u, v;

  int *ij_owner, *ij_local;
  int **ij_map, *abs_ij_map;
  int **pairdomain, *pairdom_len;

  ij_owner = get_ij_owner();
  ij_local = get_ij_local();
 
  ij_map = get_ij_map();
  abs_ij_map = original_ij_map();

  pairdomain = compute_pairdomain(ij_map);
  pairdom_len = compute_pairdomlen(ij_map);

  // Set up the indexing arrays
  // auxstart[atom] = the AO number for the first aux bf on this atom
  // auxstop[atom]  = the AO number for the last aux bf on this atom
  int *auxstart = init_int_array(natom);
  int *auxstop = init_int_array(natom);
  int *auxstart_shell = init_int_array(natom);
  int *auxstop_shell = init_int_array(natom);
  int rinshell = ribasis->nshell();
  int *aux_stype = get_aux_stype("DF_BASIS",rinshell);
  int *aux_snuc = get_aux_snuc("DF_BASIS",rinshell);
  int atom, offset, am, shell_length;
  for(i=0,atom=-1,offset=0; i < ribasis->nshell(); i++) {
    am = aux_stype[i] - 1;    // am is the angular momentum of the orbital
    // shell_length is the number of obritals in each shell
    shell_length = l_length[am];        
    // aux_snuc is the nucleus that the shell belongs to
    if(atom != aux_snuc[i]-1) {             
      if(atom != -1) {
        auxstop[atom] = offset-1;
        auxstop_shell[atom] = i-1;
      }
      atom = aux_snuc[i]-1;
      auxstart[atom] = offset;
      auxstart_shell[atom] = i;
    }
    offset += shell_length;
  }
  auxstop[atom] = offset-1;
  auxstop_shell[atom] = i-1;

  int *auxsize = init_int_array(natom);
  for (i=0; i<natom; i++) {
    auxsize[i] = auxstop[i] - auxstart[i] + 1;
  }

  int *aosize = init_int_array(natom);
  for (i=0; i<natom; i++) {
    aosize[i] = aostop[i] - aostart[i] + 1;
  }

  // aux2atom[i] = the atom number which aux bf i is on
  int* aux2atom = init_int_array(ribasis->nbf());
  for(i=0; i < natom; i++) {
    for(j=auxstart[i]; j <= auxstop[i]; j++) {
      aux2atom[j] = i;
    }
  }

  // if we use "method 1" for the fitting, we'll need to get an 
  // aux_pairdom_len here...

  /* Build the "united pair domains" (Werner JCP 118, 8149 (2003) */
  /* This will only work for now, where we have assumed that all
     pairs are strong pairs.  This assumption will change, then
     we may need to change this code --CDS 12/09
   */
  int** uniteddomain = init_int_matrix(nocc, natom);
  int** uniteddomain_len2 = init_int_matrix(nocc, natom);
  int* uniteddomain_len = init_int_array(nocc);
  int* fit_len = init_int_array(nocc);
  int* fit_atoms = init_int_array(nocc);
  int** fit_atom_id = init_int_matrix(nocc, natom);
  for(i = 0; i < nocc; i++) {
    for(k = 0; k < natom; k++) {
      for(j = 0; j < nocc; j++) {
        ij = INDEX(i, j);
        ij = abs_ij_map[ij]; // map org ij to new ij
        if(pairdomain[ij][k] && uniteddomain[i][k] == 0) {
          uniteddomain[i][k] = 1;
          uniteddomain_len[i] += aostop[k] - aostart[k] + 1;
          // method 2 fit domains, not extended (Rd=0)
          fit_len[i] += auxsize[k];
          fit_atom_id[i][fit_atoms[i]] = k;
          fit_atoms[i]++;
        }
      }
    }
  }

  int** unitedfitdomain = init_int_matrix(nocc, natom);
  for (i=0,ij=0; i<nocc; i++) {
    for (j=0; j<=i; j++, ij++) {
      if (pairdom_exist[ij]) {
        for (k=0; k<natom; k++) {
          if (uniteddomain[i][k] || uniteddomain[j][k]) {
            unitedfitdomain[i][k] = 1;
            unitedfitdomain[j][k] = 1;
        }
      }
    }
  }
 
  int** unitedfit_start = init_int_matrix(nocc, natom);
  int counter;
  for (i=0; i<nocc; i++) {
    counter = 0;
    for (k=0; k<natom; k++) {
      unitedfit_start[i][k] = -1;
      if (unitedfitdomain[i][k]) {
        unitedfit_start[i][k] = counter;
        unitedfit_len[i] += auxsize[k];
        counter += auxsize[k];
      }
    }
  }

  int **uniteddomain_abs2rel = init_int_matrix(nocc, nso);
  for(i = 0; i < nocc; i++) {
    for(k = 0, a = 0; k < natom; k++) {
      if(uniteddomain[i][k]) {
        for(t = aostart[k]; t <= aostop[k]; t++, a++) {
          uniteddomain_abs2rel[i][t] = a;
        }
      }
    }
  }


  // Schwartz Screening 
 
  IntegralFactory ao_eri_factory(basis, basis, basis, basis);
  TwoBodyInt* ao_eri = ao_eri_factory.eri();
  const double *ao_buffer = ao_eri->buffer();

  double *Schwartz = init_array(basis->nshell() * (basis->nshell()+1) / 2);
  double *DFSchwartz = init_array(ribasis->nshell());

  int numw, numx, P,Q,PQ,w,x,index;
  double tei, max;
  for(P=0,PQ=0;P<basis->nshell();P++) {
    numw = basis->shell(P)->nfunction();
    for(int Q=0;Q<=P;Q++,PQ++) {
      numx = basis->shell(Q)->nfunction();
      max=0.0;

      ao_eri->compute_shell(P, Q, P, Q);

      for(w=0;w<numw;w++) {
        for(x=0;x<numx;x++) {
          index = ( ( (w*numx + x) * numw + w) * numx + x);
          tei = ao_buffer[index];
          if(fabs(tei) > max) max = fabs(tei);
        }
      }
      Schwartz[PQ] = max;
    }
  }

  for(P=0;P<ribasis->nshell();P++) {
    numw = ribasis->shell(P)->nfunction();
    max=0.0;

    Jint->compute_shell(P, 0, P, 0);

    for(w=0;w<numw;w++) {
      tei = Jbuffer[w];
      if(fabs(tei) > max) max = fabs(tei);
    }
    DFSchwartz[P] = max;
  }

  double **SchwartzBlock = block_matrix(natom,natom);
  for(m = 0; m < natom; m++) {
    for(n = 0; n < natom; n++) { // NU block
      max = 0.0;
      for(MU = aostart_shell[m]; MU <= aostop_shell[m]; ++MU) {
        for(NU = aostart_shell[n]; NU <= aostop_shell[n]; ++NU) {

        if (NU>MU) continue;
        mn = INDEX(MU,NU); 
 
        MUNUmax = Schwartz[mn];
        if(fabs(MUNUmax) > max) max = fabs(MUNUmax);
        } 
      }
      SchwartzBlock[m][n] = max;
    }
  }

  double *DFSchwartzBlock = init_array(natom);
  for(k = 0; k < natom; k++) {
    max = 0.0;
    for(Pshell = auxstart_shell[k]; Pshell <= auxstop_shell[k]; ++Pshell){ 
      MUNUmax = DFSchwartz[Pshell];
      if(fabs(MUNUmax) > max) max = fabs(MUNUmax);
    }
    DFSchwartzBlock[k] = max;
  }

  double **Cmax = block_matrix(natom, nocc);
  double Cval, max1;
  for(i=0; i < nocc; i++) {
    for(k=0; k < natom; k++) {
      max1 = 0.0;
      for(t=aostart[k]; l <= aostop[k]; t++) {
        Cval = C[t][i];
        if(fabs(Cval) > max1) max1 = fabs(Cval);
      }
      Cmax[k][i] = max1;
    }
  }

  // this is really unefficient 
  double **C_t = block_matrix(nocc, nso);
  for(i=0; i < nocc; i++)
    for(t=0; r < nso; t++)
       C_t[i][t] = C[t][i];


  int numPshell, Pshell, MU, NU, P, oP, Q, oQ, mu, nu, nummu, numnu, omu, onu, mn;

  // find out the max number of P's in a P shell
  int screened=0;
  int max_aoblock_len  = 0;
  int max_auxblock_len = 0;
  for(k = 0; k < natom; k++){
    if (auxsize[k] > max_auxblock_len) max_auxblock_len = auxsize[k];
    if (aosize[k] > max_aoblock_len) max_aoblock_len = aosize[k];
  }

  double** temp = block_matrix(max_auxblock_len*max_aoblock_len,max_aoblock_len);
  double*  I1 = init_array(max_auxblock_len*max_aoblock_len);
  double** I2 = block_matrix(max_auxblock_len,max_aoblock_len);

  double ***I3 = (double ***) malloc(sizeof (double **) * nocc);
  for(i = 0; i < nocc; i++)
    I3[i] = (double **) malloc(sizeof (double *) * max_auxblock_len);
  for(i = 0; i < nocc; i++) {
    for(j = 0; j < max_auxblock_len; j++) {
      I3[i][j] = (double *) malloc(sizeof (double) * max_aoblock_len);
      memset(I3[i][j], '\0', sizeof (double) * max_aoblock_len);
    }
  }
           
  double* eigval = init_array(max_auxblock_len);
  int lwork = max_auxblock_len * 3;
  double* work = init_array(lwork);
  int stat, p;

  double Imax = 0.0;
  const double *buffer = eri->buffer();
  for(k = 0; k < natom; k++) { // Pshell block
    for(m = 0; m < natom; m++) { // MU block
      for(n = 0; n < natom; n++) { // NU block

        for(Pshell = auxstart_shell[k]; Pshell <= auxstop_shell[k]; ++Pshell) {
          numPshell = ribasis->shell(Pshell)->nfunction();

          for(MU = aostart_shell[m]; MU <= aostop_shell[m]; ++MU) {
            //mu_block_len = aostop[l]-aostart[l]+1; 
            nummu = basis->shell(MU)->nfunction();

            for(NU = aostart_shell[n]; NU <= aostop_shell[n]; ++NU) {
              //nu_block_len = aostop[m]-aostart[m]+1; 
              numnu = basis->shell(NU)->nfunction();

              if (NU>MU) continue;
              mn = INDEX(MU,NU);

              Imax = sqrt(Schwartz[mn]*DFSchwartz[Pshell]); 
              if( Imax > tol) {
                eri->compute_shell(Pshell, 0, MU, NU);

                for(P = 0, index = 0; P < numPshell; ++P) { 
                  oP = ribasis->shell(Pshell)->function_index() + P - auxstart[k];

                  for(mu = 0; mu < nummu; ++mu) {
                    omu = basis->shell(MU)->function_index() + mu - aostart[m];

                    for(nu = 0; nu < numnu; ++nu, ++index) {
                      onu = basis->shell(NU)->function_index() + nu - aostart[n];

                      temp[oP*aosize[n]+onu][omu] = buffer[index]; // (oP | omu onu) integral
                      temp[oP*aosize[n]+omu][onu] = buffer[index]; // (oP | omu onu) integral
                    }
                  }
                } // end loop over P in Pshell
              } // end Schwartz inequality
              else screened++;
            } // end loop over NU
          } // end lopp over MU
        } // end loop over Pshell
   
        Imax = sqrt(SchwartzBlock[l][m]*DFSchwartzBlock[k]);

        for(i = 0, ij=0; i < nocc; i++) {

          if( Cmax[l][i] * Imax > 1.0E-7 ) {
            if ( unitedfit_start[i][k] >= 0 ) {

            // first transformation (1/3 rd done)
            C_DGEMV('n',auxsize[k]*aosize[n],aosize[m],1.0,&temp[0][0],max_aoblock_len, 
                    C_t[i]+aostart[m],1,0.0,&I1[0],1);
            
            for ( p=0; p<auxsize[k]; p++) { 
              for (nu=0; nu<aosize[n]; nu++){
                 I2[p][nu] = I1[p * aosize[n] + nu];  
              }
            }

            // Second Transformation (2/3 rd done)
            for(l = 0; l < natom; l++) {
              if(uniteddomain[i][n]) {
                for(p=0;p<auxsize[k];p++) {              
                  for(t = aostart[l], a = 0; t <= aostop[l]; t++, a++) {
                    for(u = aostart[n], nu=0; u <= aostop[n]; u++, nu++) { 
                      I3[i][p][a] += Rt_full[u][t] * I2[p][nu];
                    }
                  }
                }
              }
            }
            
            } // end if ( unitedfit_start[i][k] >=0 )
          } // if Cmax[l][i] * Imax
        } // end loop over i

      }  // end loop over NUblock
    }  // end loop over MUblock
  } // end loop over Pblock


  double **J = block_matrix(max_auxblock_len,max_auxblock_len);
  double **J_inv = block_matrix(max_auxblock_len,max_auxblock_len);

//  double ***J_inv = (double ***) malloc(sizeof (double **) * nocc);
//  for(i = 0; i < nocc; i++)
//    J_inv[i] = (double **) malloc(sizeof (double *) * max_auxblock_len);
//  for(i = 0; i < nocc; i++) {
//    for(j = 0; j < max_auxblock_len; j++) {
//      J_inv[i][j] = (double *) malloc(sizeof (double) * max_auxblock_len);
//      memset(J_inv[i][j], '\0', sizeof (double) * max_auxblock_len);
//    }
//  }


  const double *Jbuffer = Jint->buffer();
  int m2, n2;
  for(i = 0; i < nocc; i++) {  

    for(m = 0; m < fit_atoms[i]; m++) { // MU block
      m2 = fit_atom_id[i][m];
      for(n = 0; n < fit_atoms[i]; n++) { // NU block
        n2 = fit_atom_id[i][n];      

        for(MU = auxstart_shell[m2]; MU <= auxstop_shell[m2]; ++MU) {
          nummu = ribasis->shell(MU)->nfunction();

          for(NU = auxstart_shell[n2]; NU <= auxstop_shell[n2]; ++NU) {
            numnu = ribasis->shell(NU)->nfunction();

            Jint->compute_shell(MU, 0, NU, 0);

            for(mu = 0, index=0; mu < nummu; ++mu) {
              omu = ribasis->shell(MU)->function_index() + mu - auxstart[m2];

              for(nu = 0; nu < numnu; ++nu, ++index) {
                onu = ribasis->shell(NU)->function_index() + nu - auxstart[n2];

                J[omu][onu] = Jbuffer[index];

              }
            }
          }  
        }
      }
    } 
    // Form J^-1/2
    // First, diagonalize J
    // the C_DSYEV call replaces the original matrix J with its eigenvectors
    //
    //
    //fprintf(outfile, "\n\tJ_mat\n");
    //print_mat(J,omu+1,onu+1, outfile);
    //fflush(outfile);              

    m2 = fit_atom_id[i][fit_atoms[i]];
 
    stat = C_DSYEV('v', 'u', m2, &J[0],max_auxblock_len, eigval, work, lwork);
    if(stat != 0) {
      throw PsiException("Error in DGESV in RI-LMP2 J-matrix construction", __FILE__, __LINE__);
    }

    // Now J contains the eigenvectors of the original J
    // Copy J to J_copy
    //J_copy = block_matrix(omu, omu);
    C_DCOPY(max_auxblock_len * max_auxblock_len, J[0], 1, J_copy[0], 1);

    // Now form J^-1 = U(T)*j^-1*U,
    // where j^-1 is the diagonal matrix of the inverse 
    // of the eigenvalues, and U is the matrix of eigenvectors of J
    for(int j = 0; j < m2; j++) {
      if(eigval[j] < 1.0E-10)
        eigval[j] = 0.0;
      else {
        eigval[j] = 1.0 / eigval[j];
      }
      // scale one set of eigenvectors by the diagonal elements j^{-1}
      C_DSCAL(m2, eigval[j], J[j], 1);
    }
        //free(eigval);

    // J_inv = J_copy(T) * J
    C_DGEMM('t', 'n', auxsize[k],auxsize[k],auxsize[k],
             1.0, J_copy[0],max_auxblock_len,,J[0],max_auxblock_len, 
             0.0, J_inv,,max_auxblock_len,);

    // D[i][a][p] = I[i][a][q] * J_inv[q][p]  
    C_DGEMM('n', 'n', auxsize[k],aosize[k],auxsize[k],
             1.0, I3[i],max_aoblock_len,,J_inv[0],max_auxblock_len, 
             0.0, d[i],max_aoblock_len,);
 
  } // end loop over i

  //Allocating the memory needed by Ktilde
  if(ij_pairs % nprocs == 0) {
    Ktilde = (double ***) malloc((ij_pairs / nprocs) * sizeof (double **));
    for(ij = 0; ij < ij_pairs; ij++) {
      if(myid == ij_owner[ij])
        Ktilde[ij_local[ij]] = block_matrix(pairdom_len[ij], pairdom_len[ij]);
    }
  } 
  else {
    if(myid < ij_pairs % nprocs) {
      Ktilde = (double ***) malloc(((ij_pairs / nprocs) + 1) * sizeof (double **));
      for(ij = 0; ij < ij_pairs; ij++) {
        if(myid == ij_owner[ij])
          Ktilde[ij_local[ij]] = block_matrix(pairdom_len[ij], pairdom_len[ij]);
      }
    } 
    else {
      Ktilde = (double ***) malloc(ij_pairs / nprocs * sizeof (double **));
      for(ij = 0; ij < ij_pairs; ij++) {
        if(myid == ij_owner[ij])
          Ktilde[ij_local[ij]] = block_matrix(pairdom_len[ij], pairdom_len[ij]);
      }
    }
  }

  //
  // construct Ktilde[ij][a][b] = (ia|jb)
  //
  v = 0;
  int a2, b2;
  /* This now loops over all ij_pairs even if we are neglecting
   * distant pairs ... um, no, ij_pairs would subtract distant ones, yes? */
  for(ij = 0; ij < ij_pairs; ij++, v++) {
    i = ij_map[ij][0];
    j = ij_map[ij][1];
    if(v % nprocs == myid) {
      for(k = 0, a = 0; k < natom; k++) {
        if(pairdomain[ij][k]) {
          for(t = aostart[k]; t <= aostop[k]; t++, a++) {
            a2 = uniteddomain_abs2rel[i][t];
            for(l = 0, b = 0; l < natom; l++) {
              if(pairdomain[ij][l]) {
                for(u = aostart[l]; u <= aostop[l]; u++, b++) {
                  b2 = uniteddomain_abs2rel[j][u];
                  // fit_atoms[i] returns the number of atoms in the
                  // fitting basis i_{fit}
                  for (int m=0, q=0; m<fit_atoms[i]; m++) {
                    // fit_atom_id[i][m] returns the absolute atom number
                    // for the mth atom belonging to the [i]fit list for
                    // atom i
                    m2 = fit_atom_id[i][m];
                    // auxsize[m2] returns the number of aux fns on atom m2
                    // ufitstart[j][m2] returns the starting location in
                    // the united orbital fit domain [j]_{fit}^u for
                    // auxiliary functions belonging to atom m2
                    Ktilde[ij_local[ij]][a][b] = 
                      C_DDOT(auxsize[m2], &(d[i][a2][q]),
                      &(I[j][b2][ufitstart[j][m2]]), 1);
                }
              }
            }
          }
        }
      }
    }
  } // End of ij_pair loop


#ifdef TIME_DF_LMP2
if(myid == 0) timer_off("Compute DF-LMP2");
#endif

}

}} // end namespaces


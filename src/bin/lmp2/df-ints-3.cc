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

void LMP2::direct_df_transformation3() {

  using namespace std;

  // Required for libmints, allocates and computes:
  // ioff, fac, df, bc
  Wavefunction::initialize_singletons();

  // Create a basis set object and initialize it using the checkpoint file.
  shared_ptr<BasisSet>basis = shared_ptr<BasisSet > (new BasisSet(chkpt));
  shared_ptr<BasisSet>ribasis = shared_ptr<BasisSet > (new BasisSet(chkpt, "DF_BASIS_MP2"));
  shared_ptr<BasisSet>zero = BasisSet::zero_basis_set();

  // Create integral factory
  IntegralFactory rifactory(ribasis, zero, basis, basis);
  IntegralFactory rifactory_J(ribasis, zero, ribasis, zero);

  // Create an integral object for ERIs
  
  TwoBodyInt* eri = rifactory.eri();
  TwoBodyInt* Jint = rifactory_J.eri();
  const double *Jbuffer = Jint->buffer();

  int i, j, ij, k, l, m, n, p, a, b, s, t, u, v;
  int MU, NU, P, oP, Q, oQ, mu, nu, nummu, numnu, omu, onu, mn;
  int numPshell, Pshell, index,stat,counter;
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
  int *aux_stype = get_aux_stype("DF_BASIS_MP2",rinshell);
  int *aux_snuc = get_aux_snuc("DF_BASIS_MP2",rinshell);
  int atom, offset, am, shell_length;

  // Compute the length of each AM block
  l_length = init_int_array(LIBINT_MAX_AM);
  l_length[0] = 1;
  for(l=1; l < (LIBINT_MAX_AM); l++) {
    if(puream) l_length[l] = 2 * l + 1;
    else l_length[l] = l_length[l-1] + l + 1;
  }

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
  
  fprintf(outfile, "\n\n\t\t NATOM aostart aostop aosize\n");
  fprintf(outfile, "\t\t===============\n");
  for(i=0;i<natom;i++)
    fprintf(outfile,"\t\t  %5d  %5d  %5d  %5d\n",i,aostart[i], aostop[i], aosize[i]);

  fprintf(outfile, "\n\n\t\t NATOM auxstart auxstop auxsize\n");
  fprintf(outfile, "\t\t================\n");
  for(i=0;i<natom;i++)
    fprintf(outfile,"\t\t  %5d  %5d  %5d  %5d\n",i,auxstart[i], auxstop[i], auxsize[i]);
  fflush(outfile);

  // aux2atom[i] = the atom number which aux bf i is on
  int nribasis = ribasis->nbf();
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
     we may need to change this code 
  */
  int** uniteddomain = init_int_matrix(nocc, natom);
  int* uniteddomain_len = init_int_array(nocc);
  int* fit_len = init_int_array(nocc);
  int* fit_atoms = init_int_array(nocc);
  int** fit_atom_id = init_int_matrix(nocc, natom);

  for(i = 0; i < nocc; i++) {
    for(k = 0; k < natom; k++) {
      counter = 1; 
      for(j = 0; j < nocc; j++) {
        ij = INDEX(i,j);
        ij = abs_ij_map[ij]; // map org ij to new ij
        if (ij == -1)
		continue;
	if( pairdomain[ij][k] == 1 && uniteddomain[i][k] == 0 ) {
          uniteddomain[i][k] = 1;
          uniteddomain_len[i] += aosize[k];
          // method 2 fit domains, not extended (Rd=0)
          fit_len[i] += auxsize[k];
          fit_atoms[i] = counter;
          fit_atom_id[i][fit_atoms[i]] = k;
          counter++;
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
  }
    
  int** unitedfit_start = init_int_matrix(nocc, natom);
  int* unitedfit_len = init_int_array(nocc);
  for (i=0; i<nocc; i++) {
    counter = 0;
    for (k=0; k<natom; k++) {
      unitedfit_start[i][k] = -1;
      if (unitedfitdomain[i][k]) {
	fprintf(outfile,"  i = %d, k = %d, counter = %d, auxsize = %d\n",i,k,counter,auxsize[k]);
        unitedfit_start[i][k] = counter;
        unitedfit_len[i] += auxsize[k];
      }
      counter += auxsize[k];
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
  

fprintf(outfile, "\n\n NOCC NATOM uniteddomain_len[i] fit_len[i] fit_atom_id fit_atoms\n");
fprintf(outfile, "\t\t================\n");

for(i = 0; i < nocc; i++) {
  for(k = 0; k < natom; k++) {
    fprintf(outfile,"\t\t  %5d  %5d  %5d  %5d  %5d  %5d\n",
    i,k,uniteddomain_len[i],fit_len[i],
    fit_atom_id[i][fit_atoms[i]], fit_atoms[i]);
  }
}
fflush(outfile);

fprintf(outfile, "\n\n NOCC NATOM unitedfitdomain[i][k] unitedfit_start[i][k] unitedfit_len[i]\n");
fprintf(outfile, "\t\t================\n");
for(i = 0; i < nocc; i++) {
  for(k = 0; k < natom; k++) {
    fprintf(outfile,"\t\t  %5d  %5d  %5d  %5d  %5d \n",
           i,k,unitedfitdomain[i][k],unitedfit_start[i][k],unitedfit_len[i]);
  }
}


   fprintf(outfile, "\nBook keeping done...\n");
   fflush(outfile);

  // Schwartz Screening 
  
  IntegralFactory ao_eri_factory(basis, basis, basis, basis);
  TwoBodyInt* ao_eri = ao_eri_factory.eri();
  const double *ao_buffer = ao_eri->buffer();
  
  double *Schwartz = init_array(basis->nshell() * (basis->nshell()+1) / 2);
  double *DFSchwartz = init_array(ribasis->nshell());
  
  int numw, numx, PQ,w,x;
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
  
   fprintf(outfile, "\nAO Schwartz Screening   done...\n");
   fflush(outfile);

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
  
   fprintf(outfile, "\nAUX Schwartz Screening   done...\n");
   fflush(outfile);
  double **SchwartzBlock = block_matrix(natom,natom);
  double MUNUmax;
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
      for(t=aostart[k]; t <= aostop[k]; t++) {
        Cval = C[t][i];
        if(fabs(Cval) > max1) max1 = fabs(Cval);
      }
      Cmax[k][i] = max1;
    }
  }
  
  
  // this is really unefficient 
  double **C_t = block_matrix(nocc, nso);
  for(i=0; i < nocc; i++)
    for(t=0; t < nso; t++)
      C_t[i][t] = C[t][i];
  
   fprintf(outfile, "\nBook keeping done...\n");
   fflush(outfile);
  
  
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
    I3[i] = block_matrix(uniteddomain_len[i],fit_len[i]);
  //for(i = 0; i < nocc; i++)
    //I3[i] = (double **) malloc(sizeof (double *) * uniteddomain_len[i]);
  //for(i = 0; i < nocc; i++) {
  //  for(j = 0; j < uniteddomain_len[i]; j++) {
  //    I3[i][j] = (double *) malloc(sizeof (double) * fit_len[i]);
  //    memset(I3[i][j], '\0', sizeof (double) * fit_len[i]);
  //  }
  //}
  

  double* eigval = init_array(max_auxblock_len);
  int lwork = max_auxblock_len * 3;
  double* work = init_array(lwork);
  int* AOsize = init_int_array(nocc);
  
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
	
        Imax = sqrt(SchwartzBlock[m][n]*DFSchwartzBlock[k]);
	
        for(i = 0, ij=0; i < nocc; i++) {
	  
          if( Cmax[m][i] * Imax > 1.0E-7 ) {
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
	      
              for(p=0;p<auxsize[k];p++) {              
	        for(l = 0,a=0; l < natom; l++) {
		  if(uniteddomain[i][l]) {
		    for(t = aostart[l]; t <= aostop[l]; t++, a++) {
		      for(u = aostart[n], nu=0; u <= aostop[n]; u++, nu++) { 
  		        I3[i][a][p+unitedfit_start[i][k]] += Rt_full[u][t] * I2[p][nu];
		      }
		    } 
		  }
		}
	      } 
              AOsize[i] = a;
	      fprintf(outfile,"  i = %d AO Size %d\n",i,a); 
            } // end if ( unitedfit_start[i][k] >=0 )
          } // if Cmax[l][i] * Imax
        } // end loop over i
	
      }  // end loop over NUblock
    }  // end loop over MUblock
  } // end loop over Pblock
  

  int max_fit_len = 0;
  for(i = 0; i < nocc; i++){
    if (fit_len[i] >  max_fit_len) max_fit_len = fit_len[i];
  }

  double **J = block_matrix(max_fit_len,max_fit_len);
  double **J_inv = block_matrix(max_fit_len,max_fit_len);
  
  double ***D = (double ***) malloc(sizeof (double **) * nocc);
  for(i = 0; i < nocc; i++)
    D[i] = block_matrix(uniteddomain_len[i],fit_len[i]);
  //for(i = 0; i < nocc; i++)
    //D[i] = (double **) malloc(sizeof (double *) * uniteddomain_len[i]);
  //for(i = 0; i < nocc; i++) {
    //for(j = 0; j < uniteddomain_len[i]; j++) {
      //D[i][j] = (double *) malloc(sizeof (double) * fit_len[i]);
      //memset(D[i][j], '\0', sizeof (double) * fit_len[i]);
    //}
  //}
  
  int m2, n2, Jsize;
  for(i = 0; i < nocc; i++) {  
   
    for(m = 0; m < fit_atoms[i]; m++) { // MU block
      m2 = fit_atom_id[i][m];

      for(MU = auxstart_shell[m2]; MU <= auxstop_shell[m2]; ++MU) {
        nummu = ribasis->shell(MU)->nfunction();

        for(n = 0; n < fit_atoms[i]; n++) { // NU block
          n2 = fit_atom_id[i][n];      
	  
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
            
          //onu += numnu;
          } // end loop over NU shell
          
        } 
      //omu += nummu;
      } // end loop over MU shell
    } 

    // Form J^-1
    // The plan is to perform an LU decomposition via dgetrf
    // and them calculate the inverse via dgetri
    //m2 = fit_atom_id[i][fit_atoms[i]];
    //Jsize = auxsize[m2];
    Jsize = omu+1;
    
    //fprintf(outfile,"  J matrix i = %d ,Jsize = %d: \n\n",i, Jsize);
    //print_mat(J,Jsize,Jsize,outfile);
    
    int M = max_aoblock_len;
    int N = max_fit_len;
    int* ipiv = init_int_array(N);
    //memset(ipiv, '\0', sizeof (int) * N);
    stat = C_DGETRF(Jsize,Jsize,&J[0][0],N,ipiv);
    if(stat != 0) {
      fflush(outfile);
      throw PsiException("Error in DGESV in RI-LMP2 J-matrix at C_DGETRF() ",
			 __FILE__, __LINE__);
    }

    int lwork =4 * N;
    double* work = init_array(lwork);
    stat = C_DGETRI(Jsize,&J[0][0],N,ipiv,work,lwork);
    if(stat != 0) {
      fflush(outfile);
      throw PsiException("Error in DGESV in RI-LMP2 J-matrix at C_DGETRI() ",
			 __FILE__, __LINE__);
    }
    
    //fprintf(outfile,"  J matrix Inverse i = %d: \n\n",i);
    //print_mat(J,Jsize,Jsize,outfile);

    fprintf(outfile, "\nBook keeping done... %d\n",i);
    fflush(outfile);
    //fprintf(outfile,"  I3 Matrix, AOsize = %d fit_len[%d] = %d ,max_fit_len = %d \n\n",AOsize[i],i, fit_len[i],max_fit_len);
    //print_mat(I3[i], AOsize[i], Jsize, outfile);

    //Pointer based writeout
    //for(int kk = 0; kk<AOsize[i]*Jsize; kk++)
    //    fprintf(outfile,"  %d = %14.10f\n",kk, *(I3[i][0]+kk));

    // D[i][a][p] = I[i][a][q] * J_inv[q][p]  

   //for (int AA = 0; AA<AOsize[i]; AA++)
   //	for (int Q = 0; Q<Jsize; Q++)
   //	    for (int P = 0; P<Jsize; P++)
   //		D[i][AA][Q] += I3[i][AA][P]*J[P][Q];


    C_DGEMM('n', 'n', AOsize[i], Jsize, Jsize,
            1.0, &I3[i][0][0], fit_len[i], &J[0][0],max_fit_len, 
             0.0, &D[i][0][0], fit_len[i]);
    //C_DGEMM('n', 'n', 5, 9, 9,
    //        1.0, I3[i][0], 9, J[0],9, 
    //         0.0, D[i][0], 9);

    
    //fprintf(outfile,"  D Matrix, AOsize = %d\n\n",AOsize[i]);
    //print_mat(D[i], AOsize[i], Jsize, outfile);
    //fflush(outfile);

  } // end loop over i




  
  //Allocating the memory needed by Ktilde
  Ktilde = (double ***) malloc(pairs_per_proc * sizeof (double **));
  for(ij = 0; ij < ij_pairs; ij++) {
    if (Communicator::world->me() == ij_owner[ij])
      Ktilde[ij_local[ij]] = block_matrix(pairdom_len[ij], pairdom_len[ij]);
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
    if (Communicator::world->me() == ij_owner[ij]) {
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
                             C_DDOT(auxsize[m2], &(D[i][a2][q]),1,
			     &(I3[j][b2][unitedfit_start[j][m2]]), 1);
                     q += auxsize[m2];
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  } // End of ij_pair loop
  
  
//#ifdef TIME_DF_LMP2
//  if(Communicator::world->me() == 0) 
//    timer_off("Compute DF-LMP2");
//#endif

}

}} // end namespaces


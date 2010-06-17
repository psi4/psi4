/*
 *
 */

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>

#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.hpp>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>
#include <psifiles.h>

#include "dfmp2.h"

#include <libmints/basisset.h>
#include <libmints/onebody.h>
#include <libmints/twobody.h>
#include <libmints/integral.h>
#include <libmints/molecule.h>

using namespace boost;
using namespace std;
using namespace psi;

namespace psi { namespace dfmp2 {

DFMP2::DFMP2(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt_)
    : Wavefunction(options, psio, chkpt_)
{
}

DFMP2::~DFMP2()
{
}
double DFMP2::compute_E()
{
  fprintf(outfile, "\t\t\t*************************\n");
  fprintf(outfile, "\t\t\t*                       *\n");
  fprintf(outfile, "\t\t\t*         DF-MP2        *\n");
  fprintf(outfile, "\t\t\t*                       *\n");
  fprintf(outfile, "\t\t\t*************************\n");
  fflush(outfile);

  // Required for libmints, allocates and computes:

  // Create a new matrix factory
  MatrixFactory factory;

  // Initialize the factory with data from checkpoint
  factory.init_with_chkpt(chkpt_);

  // Form ribasis object:
  shared_ptr<BasisSet> ribasis =shared_ptr<BasisSet>(new BasisSet(chkpt_, "DF_BASIS_MP2"));
  int naux = ribasis->nbf();
  shared_ptr<BasisSet> zero = BasisSet::zero_basis_set();
  
  // Taken from mp2 module, get_params.cc 
  // get parameters related to SCS-MP2 or SCS-N-MP2 
  // see papers by S. Grimme or J. Platz 
  int scs = options_.get_bool("SCS");
  int scsn= options_.get_bool("SCS_N");
  double scs_scale_os = 6.0/5.0;
  double scs_scale_ss = 1.0/3.0;
  if (scs == 1) {
    scs_scale_os = options_.get_double("SCALE_OS");
    scs_scale_ss = options_.get_double("SCALE_SS");
    fprintf(outfile,"\n  Spin-Component Scaled RI-MP2 requested\n");
    fprintf(outfile,"    Opposite-spin scaled by %10.4lf\n",scs_scale_os);
    fprintf(outfile,"    Same-spin scaled by     %10.4lf\n",scs_scale_ss);
    }
  if (scsn == 1) {
    fprintf(outfile,"\n  SCSN-RI-MP2 energies will be printed\n");
  }
    
  int nirreps = chkpt_->rd_nirreps();
  int *clsdpi = chkpt_->rd_clsdpi();
  int *orbspi = chkpt_->rd_orbspi();
  int *frzcpi = chkpt_->rd_frzcpi();
  int *frzvpi = chkpt_->rd_frzvpi();
  int ndocc = 0;
  int nvirt = 0;
  int nfocc = 0;
  int nfvir = 0;
  int norbs = 0;
  int nact_docc = 0;
  int nact_virt = 0;
  for(int h=0; h < nirreps; ++h){
      nfocc     += frzcpi[h];
      nfvir     += frzvpi[h];
      ndocc     += clsdpi[h];
      nact_docc += clsdpi[h] - frzcpi[h];
      nvirt     += orbspi[h] - clsdpi[h];
      nact_virt += orbspi[h] - frzvpi[h] - clsdpi[h];
      norbs     += orbspi[h];
  }

  fprintf(outfile, "\n\t\t==============================================\n");
  fprintf(outfile, "\t\t #ORBITALS #RI  FOCC DOCC AOCC AVIR VIRT FVIR \n");
  fprintf(outfile, "\t\t----------------------------------------------\n");
  fprintf(outfile, "\t\t  %5d  %5d  %4d %4d %4d %4d %4d %4d\n",
          norbs,ribasis->nbf(),nfocc,ndocc,nact_docc,nact_virt,nvirt,nfvir);
  fprintf(outfile, "\t\t==============================================\n");

  double escf = chkpt_->rd_escf();

  // Read in C coefficients

  SimpleMatrix *C_so = 
    factory.create_simple_matrix("MO coefficients (SO basisset_)");
  double **vectors = chkpt_->rd_scf();
  if (vectors == NULL) {
    fprintf(stderr, "Could not find MO coefficients. Run cscf first.\n");
    return EXIT_FAILURE;
  }
  C_so->set(vectors);
  free_block(vectors);
  double **Umat = chkpt_->rd_usotbf();
  SimpleMatrix *U = factory.create_simple_matrix("SO->BF");
  U->set(Umat);
  free_block(Umat);
  // Transfrom the eigenvectors from the SO to the AO basisset_
  SimpleMatrix *C = 
    factory.create_simple_matrix("MO coefficients (AO basisset_)");
  C->gemm(true, false, 1.0, U, C_so, 0.0);

  // Load in orbital energies
  double *orbital_energies = chkpt_->rd_evals();

  // Create integral factory
  IntegralFactory rifactory(ribasis, zero, basisset_, basisset_);
  IntegralFactory rifactory_J(ribasis, zero, ribasis, zero);
    
  TwoBodyInt* eri = rifactory.eri();
  TwoBodyInt* Jint = rifactory_J.eri();
  double **J = block_matrix(ribasis->nbf(), ribasis->nbf());
  double **J_mhalf = block_matrix(ribasis->nbf(), ribasis->nbf());
  const double *Jbuffer = Jint->buffer();

  timer_on("Form J");

  int index = 0;
    
  for (int MU=0; MU < ribasis->nshell(); ++MU) {
    int nummu = ribasis->shell(MU)->nfunction();
        
    for (int NU=0; NU < ribasis->nshell(); ++NU) {
      int numnu = ribasis->shell(NU)->nfunction();
    
      Jint->compute_shell(MU, 0, NU, 0);
            
      index = 0;
      for (int mu=0; mu < nummu; ++mu) {
        int omu = ribasis->shell(MU)->function_index() + mu;
                
        for (int nu=0; nu < numnu; ++nu, ++index) {
          int onu = ribasis->shell(NU)->function_index() + nu;

          J[omu][onu] = Jbuffer[index];
        }
      }
    }
  }
    
  // fprintf(outfile, "J:\n");
  // print_mat(J, ribasis->nbf(), ribasis->nbf(), outfile);    
  // Invert J matrix
  // invert_matrix(J, J_inverse, ribasis->nbf(), outfile);
  // fprintf(outfile, "J^-1:\n");
  // print_mat(J_inverse, ribasis->nbf(), ribasis->nbf(), outfile);
  // Form J^-1/2
  // First, diagonalize J
  // the C_DSYEV call replaces the original matrix J with its eigenvectors
  double* eigval = init_array(ribasis->nbf());
  int lwork = ribasis->nbf() * 3;
  double* work = init_array(lwork);
  int stat = C_DSYEV('v','u',ribasis->nbf(),J[0],ribasis->nbf(),eigval,
    work,lwork);
  if (stat != 0) {
    fprintf(outfile, "C_DSYEV failed\n");
    exit(PSI_RETURN_FAILURE);
  }
  free(work);

  // Now J contains the eigenvectors of the original J
  // Copy J to J_copy
  double **J_copy = block_matrix(ribasis->nbf(), ribasis->nbf());
  C_DCOPY(ribasis->nbf()*ribasis->nbf(),J[0],1,J_copy[0],1); 
  
  // Now form J^{-1/2} = U(T)*j^{-1/2}*U,
  // where j^{-1/2} is the diagonal matrix of the inverse square roots
  // of the eigenvalues, and U is the matrix of eigenvectors of J
  for (int i=0; i<ribasis->nbf(); i++) {
    if (eigval[i] < 1.0E-10)
      eigval[i] = 0.0;
    else {
      eigval[i] = 1.0 / sqrt(eigval[i]);
    }
    // scale one set of eigenvectors by the diagonal elements j^{-1/2}
    C_DSCAL(ribasis->nbf(), eigval[i], J[i], 1);
  }
  free(eigval);

  // J_mhalf = J_copy(T) * J
  C_DGEMM('t','n',ribasis->nbf(),ribasis->nbf(),ribasis->nbf(),1.0,
    J_copy[0],ribasis->nbf(),J[0],ribasis->nbf(),0.0,J_mhalf[0],ribasis->nbf());

  free_block(J);
  free_block(J_copy);

  timer_off("Form J");

  const double *buffer = eri->buffer();
  
  double** Co   = block_matrix(norbs, nact_docc);
  double** Cv   = block_matrix(norbs, nact_virt);
  double** half = block_matrix(nact_docc, norbs);
  double* epsilon_act_docc = new double[nact_docc];
  double* epsilon_act_virt = new double[nact_virt];
  int*    docc_sym = new int[nact_docc];
  int*    virt_sym = new int[nact_virt];
  int offset = 0;
  int act_docc_count  = 0;
  int act_virt_count  = 0;
  for(int h=0; h<nirreps; ++h){
    // Skip over the frozen core orbitals in this irrep
    offset += frzcpi[h];
    // Copy over the info for active occupied orbitals
    for(int i=0; i<clsdpi[h]-frzcpi[h]; ++i){
      for (int mu=0; mu<norbs; ++mu)
        Co[mu][act_docc_count] = C->get(mu, offset);
      epsilon_act_docc[act_docc_count] = orbital_energies[offset];
      docc_sym[act_docc_count] = h;
      ++act_docc_count;
      ++offset;
    }
    // Copy over the info for active virtual orbitals
    for(int a=0; a<orbspi[h]-clsdpi[h]-frzvpi[h]; ++a){
      for (int mu=0; mu<norbs; ++mu)
        Cv[mu][act_virt_count] = C->get(mu, offset);
      epsilon_act_virt[act_virt_count] = orbital_energies[offset];
      virt_sym[act_virt_count] = h;
      ++offset;
      ++act_virt_count;
    }
    // Skip over the frozen virtual orbitals in this irrep
    offset += frzvpi[h];
  }
  //Chkpt::free(orbital_energies);

  //fprintf(outfile, "Co:\n");
  //print_mat(Co, norbs, nact_docc, outfile);
  //fprintf(outfile, "Cv:\n");
  //print_mat(Cv, norbs, nact_virt, outfile);

  //How many doubles do we have?
  double doubles = memory_/(1.0*sizeof(double))*0.9;
  double doubles_available = doubles-(1.0)*naux*naux;

  //What is the size of the virtual block?
  int vir_blk_size = (int) floor(doubles_available/(3.0*naux*nact_docc));
  if (vir_blk_size > nact_virt)
      vir_blk_size = nact_virt;
  if (vir_blk_size < 1)
      vir_blk_size = 1;
  int full_blks = nact_virt/vir_blk_size;
  int gimp_blk_size = nact_virt- full_blks*vir_blk_size;

  // find out the max number of P's in a P shell
  int maxPshell = 0;
  for (int Pshell=0; Pshell < ribasis->nshell(); ++Pshell) {
    int numPshell = ribasis->shell(Pshell)->nfunction();
    if (numPshell > maxPshell) maxPshell = numPshell;
  }

  //ao, three-index tensor (A|mn)
  double** ao_buffer = block_matrix(maxPshell,norbs*norbs);
  //half-tansformed three-index tensor (A|mi) buffer  
  double** half_buffer = block_matrix(maxPshell,norbs*nact_docc);
 
  //Transformed DF Integrals by blocks of a 
  double** Pia = block_matrix(naux,nact_docc*nact_virt);
  double** Pjb = block_matrix(naux,nact_docc*nact_virt);

  timer_on("Form mo_p_ia");


int numPshell,Pshell,MU,NU,P,mu,nu,nummu,numnu,omu,onu,omuonu;
for (Pshell=0; Pshell < ribasis->nshell(); ++Pshell) {
    numPshell = ribasis->shell(Pshell)->nfunction();
    //Zero that guy out!
    memset((void*)&ao_buffer[0][0],'\0',maxPshell*norbs*norbs); 
 
    for (MU=0; MU < basisset_->nshell(); ++MU) {
      nummu = basisset_->shell(MU)->nfunction();
      for (NU=0; NU <= MU; ++NU) {
        numnu = basisset_->shell(NU)->nfunction();
        eri->compute_shell(Pshell, 0, MU, NU);
        for (P=0, index=0; P < numPshell; ++P) {
          for (mu=0; mu < nummu; ++mu) {
            omu = basisset_->shell(MU)->function_index() + mu;
            for (nu=0; nu < numnu; ++nu, ++index) {
              onu = basisset_->shell(NU)->function_index() + nu;
              ao_buffer[P][omu*norbs+onu] = buffer[index]; // (oP | omu onu) integral
              ao_buffer[P][onu*norbs+omu] = buffer[index]; // (oP | onu omu) integral
            }
          }
        } // end loop over P in Pshell
      } // end loop over NU shell
    } // end loop over MU shell
    // now we've gone through all P, mu, nu for a given Pshell
    // transform the integrals for all P in the given P shell

    // Do transform
    C_DGEMM('N', 'N', numPshell*norbs, nact_docc, norbs, 1.0, &(ao_buffer[0][0]),
        norbs, &(Co[0][0]), nact_docc, 0.0, &(half_buffer[0][0]), nact_docc);

    for (int P=0, index=0; P < numPshell; ++P) {
      int oP = ribasis->shell(Pshell)->function_index() + P;

      C_DGEMM('T', 'N', nact_virt, nact_docc, norbs, 1.0, &(Cv[0][0]),
        nact_virt, &(half_buffer[P][0]), nact_docc, 0.0, &(Pia[oP][0]), nact_docc);
      C_DGEMM('T', 'N', nact_virt, nact_docc, norbs, 1.0, &(Cv[0][0]),
        nact_virt, &(half_buffer[P][0]), nact_docc, 0.0, &(Pjb[oP][0]), nact_docc);
                    
      // Do transform
//      C_DGEMM('T', 'N', nact_docc, nact_virt, norbs, 1.0, &(half_buffer[P][0]),
//        nact_docc, &(Cv[0][0]), nact_virt, 0.0, &(Pia[oP][0]), nact_virt);
//      C_DGEMM('T', 'N', nact_docc, nact_virt, norbs, 1.0, &(half_buffer[P][0]),
//        nact_docc, &(Cv[0][0]), nact_virt, 0.0, &(Pjb[oP][0]), nact_virt);
         
    }
  } // end loop over P shells; done with forming MO basisset_ (P|ia)'s
   
  free_block(ao_buffer);
  free_block(half_buffer);
 
  timer_off("Form mo_p_ia");

  // fprintf(outfile, "mo_p_ia:\n");
  // print_mat(mo_p_ia, ribasis->nbf(), nact_docc*nact_virt, outfile);
    
  timer_on("Form B_ia^P");

  // mo_p_ia has integrals
  // B_ia^P = Sum_Q (i a | Q) (J^-1/2)_QP
  double **Qia = block_matrix(nact_docc * nact_virt, ribasis->nbf());

  C_DGEMM('T','N',nact_docc*nact_virt,ribasis->nbf(),ribasis->nbf(),
    1.0, Pia[0], nact_docc*nact_virt, J_mhalf[0], ribasis->nbf(),
    0.0, Qia[0], ribasis->nbf());

  double **Qjb = block_matrix(nact_docc * nact_virt, ribasis->nbf());

  C_DGEMM('T','N',nact_docc*nact_virt,ribasis->nbf(),ribasis->nbf(),
    1.0, Pjb[0], nact_docc*nact_virt, J_mhalf[0], ribasis->nbf(),
    0.0, Qjb[0], ribasis->nbf());

  free_block(J_mhalf);

  timer_off("Form B_ia^P");

  // fprintf(outfile, "B_ia_p:\n");
  // print_mat(B_ia_p, nact_docc*nact_virt, ribasis->nbf(), outfile);


  timer_on("Compute EMP2");

  double emp2 = 0.0; 

  double **I = block_matrix(nact_docc,nact_docc);

  for (int a=0; a < nact_virt; ++a) {
    for (int b=0; b < nact_virt; ++b) {

      C_DGEMM('N','T',nact_docc,nact_docc,ribasis->nbf(),1.0,&(Qia[a*nact_docc][0]),
        ribasis->nbf(),&(Qjb[b*nact_docc][0]),ribasis->nbf(),0.0,&(I[0][0]),nact_docc);

      double iajb, ibja, tval, denom;
      for (int i=0; i < nact_docc; ++i) {
        for (int j=0; j < nact_docc; ++j) {
          iajb = I[i][j];
          ibja = I[j][i];
          denom = 1.0 /
            (epsilon_act_docc[i] + epsilon_act_docc[j] -
             epsilon_act_virt[a] - epsilon_act_virt[b]);

          tval = (2.0*iajb*iajb - iajb*ibja)*denom;
          emp2 += tval;
        }
      }
    } // end loop over j<=i
  } // end loop over i

/*
  double **I = block_matrix(nact_virt,nact_virt);

      for (int i=0; i < nact_docc; ++i) {
        for (int j=0; j < nact_docc; ++j) {

      C_DGEMM('N','T',nact_virt,nact_virt,ribasis->nbf(),1.0,&(Qia[i*nact_virt][0]),
        ribasis->nbf(),&(Qjb[j*nact_virt][0]),ribasis->nbf(),0.0,&(I[0][0]),nact_virt);

      double iajb, ibja, tval, denom;
  for (int a=0; a < nact_virt; ++a) {
    for (int b=0; b < nact_virt; ++b) {
          iajb = I[a][b];
          ibja = I[b][a];
          denom = 1.0 /
            (epsilon_act_docc[i] + epsilon_act_docc[j] -
             epsilon_act_virt[a] - epsilon_act_virt[b]);

          tval = (2.0*iajb*iajb - iajb*ibja)*denom;
          emp2 += tval;
        }
      }
    } // end loop over j<=i
  } // end loop over i
*/
  timer_off("Compute EMP2");
  
  free_block(I);

  fprintf(outfile,"\tRI-MP2 correlation energy         = %20.15f\n",emp2);
  fprintf(outfile,"      * RI-MP2 total energy               = %20.15f\n\n",
    escf + emp2);
//fprintf(outfile,"\tOpposite-Spin correlation energy  = %20.15f\n",os_mp2);
//fprintf(outfile,"\tSame-Spin correlation energy      = %20.15f\n\n",ss_mp2);
//fprintf(outfile,"      * SCS-RI-MP2 total energy           = %20.15f\n\n",
//  escf + scs_scale_os*os_mp2 + scs_scale_ss*ss_mp2);
//if (scsn) {
//  fprintf(outfile,"      * SCSN-RI-MP2 total energy          = %20.15f\n\n",
//    escf + 1.76*ss_mp2);
//  }
  
  
  // Shut down psi
  fflush(outfile);
  return emp2;

}

}}

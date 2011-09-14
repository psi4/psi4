/****************************************************************************
*
*       dfmp2.cc in psi4/src/bin/dfmp2
*       By Rob Parrish, Ed Hohenstein, David Sherrill, CCMST Georgia Tech
*       Justin Turney and Andy Simmonett CCQC University of Georgia
*       contact: robparrish@gmail.com
*       28 June 2010
*
*       Canonical (delocalized) density fitting routines for MP2 and SCS-MP2
*
*****************************************************************************/

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
#include <float.h>

#include "dfmp2.h"

#include <libmints/mints.h>
#include <libparallel/parallel.h>

//MKL Header
#ifdef _MKL
#include <mkl.h>
#endif

//OpenMP Header
//_OPENMP is defined by the compiler if it exists
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace boost;
using namespace std;
using namespace psi;

namespace psi { namespace dfmp2 {

DFMP2::DFMP2(Options& options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt_)
    : Wavefunction(options, psio, chkpt_)
{
    setup();
}
void  DFMP2::setup()
{
  timer_on("Setup");

  boost::shared_ptr<Wavefunction> ref = Process::environment.reference_wavefunction();
  boost::shared_ptr<Matrix> C;
  boost::shared_ptr<Vector> epsilon;

  if (ref.get() != NULL) {
      E_scf_ = Process::environment.globals["SCF ENERGY"];
      nirrep_ = ref->nirrep();
      nso_ = ref->nso();
      clsdpi_ = ref->doccpi();
      orbspi_ = ref->nmopi();
      frzcpi_ = ref->frzcpi();
      frzvpi_ = ref->frzvpi();

  } else {
      // Form primary basis indexing
//      if (Communicator::world->me() == 0) {
          E_scf_ = chkpt_->rd_escf();
          nirrep_ = chkpt_->rd_nirreps();
          clsdpi_ = chkpt_->rd_clsdpi();
          orbspi_ = chkpt_->rd_orbspi();
          frzcpi_ = chkpt_->rd_frzcpi();
          frzvpi_ = chkpt_->rd_frzvpi();
          nso_ = chkpt_->rd_nso();

//      }	else {
//          // Leak a little memory in parallel
//          clsdpi_ = new int[8];
//          orbspi_ = new int[8];
//          frzcpi_ = new int[8];
//          frzvpi_ = new int[8];
//      }
//      Communicator::world->bcast(E_scf_, 0);
//      Communicator::world->bcast(nirrep_, 0);
//      Communicator::world->bcast(clsdpi_, nirrep_, 0);
//      Communicator::world->bcast(orbspi_, nirrep_, 0);
//      Communicator::world->bcast(frzcpi_, nirrep_, 0);
//      Communicator::world->bcast(frzvpi_, nirrep_, 0);
//      Communicator::world->bcast(nso_, 0);
  }

  ndocc_ = 0;
  nvirt_ = 0;
  nf_docc_ = 0;
  nf_virt_ = 0;
  nact_docc_ = 0;
  nact_virt_ = 0;
  for(int h=0; h < nirrep_; ++h){
      nf_docc_     += frzcpi_[h];
      nf_virt_     += frzvpi_[h];
      ndocc_     += clsdpi_[h];
      nact_docc_ += clsdpi_[h] - frzcpi_[h];
      nvirt_     += orbspi_[h] - clsdpi_[h];
      nact_virt_ += orbspi_[h] - frzvpi_[h] - clsdpi_[h];
  }
 // printf("nact_docc_ is %d before anything happens\n", nact_docc_);
  nmo_ = ndocc_+nvirt_;

  // If the user doesn't spec a basis name, pick it yourself
  // TODO: Verify that the basis assign does not messs this up
  if (options_.get_str("RI_BASIS_MP2") == "") {
      molecule_->set_basis_all_atoms(options_.get_str("BASIS") + "-RI", "RI_BASIS_MP2");
      fprintf(outfile, "  No auxiliary basis selected, defaulting to %s-RI\n\n", options_.get_str("BASIS").c_str());
      fflush(outfile);
  }

  // Form ribasis object and auxiliary basis indexing:
  boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
  ribasis_ = BasisSet::construct(parser, molecule_, "RI_BASIS_MP2");
  naux_raw_ = ribasis_->nbf();
  naux_fin_ = ribasis_->nbf(); //For now, may be pared later
  zerobasis_ = BasisSet::zero_ao_basis_set();

  // Form initial energy values
  E_ss_ = 0.0;
  E_os_ = 0.0;
  E_ = 0.0;
  E_tot_ = 0.0;
  E_sss_ = 0.0;
  E_sos_ = 0.0;
  E_scs_ = 0.0;
  E_tot_scs_ = 0.0;

  //Form scaling parameters
  os_scale_ = options_.get_double("SCALE_OS");
  ss_scale_ = options_.get_double("SCALE_SS");

  // Form print level
  print_ = options_.get_int("PRINT");

  algorithm_type_ = options_.get_str("DFMP2_TYPE");

  fitting_symmetry_ = options_.get_str("FITTING_SYMMETRY");
  fitting_conditioning_ = options_.get_str("FITTING_CONDITIONING");
  fitting_inversion_ = options_.get_str("FITTING_INVERSION");

 // Form the initial memory usage estimate
  df_memory_ = (unsigned long int)(options_.get_double("DFMP2_MEM_FACTOR")*memory_/sizeof(double));

  //Required core memory
  unsigned long int core_memory = naux_raw_*nact_docc_*(ULI)nact_virt_ + naux_raw_*(ULI)naux_raw_;

  if (algorithm_type_ == "DEFAULT")
    if (core_memory < df_memory_)
        algorithm_type_ = "CORE";
    else
        algorithm_type_ = "DISK";

  if (ref.get() != NULL) {
      C = ref->Ca();
      epsilon = ref->epsilon_a();
  } else {
      C = boost::shared_ptr<Matrix>(new Matrix("C Matrix", nso_, nmo_));
      epsilon = boost::shared_ptr<Vector>(new Vector(nmo_));

//      if (Communicator::world->me() == 0) {
        double **vectors;
        double *orbital_energies;

        // Read in C coefficients
        vectors = chkpt_->rd_scf();

        // Load in orbital energies
        orbital_energies = chkpt_->rd_evals();

        C_DCOPY(nso_*nmo_, vectors[0], 1, C->pointer()[0], 1);
        C_DCOPY(nmo_, orbital_energies, 1, epsilon->pointer(), 1);

        // Ironically, I leak less memory
        free(orbital_energies);
        free_block(vectors);

//      }
//      Communicator::world->bcast(C->pointer()[0], nmo_*nso_, 0);
//      Communicator::world->bcast(epsilon->pointer(), nmo_, 0);

  }

  //Put the orbitals and C matrix into some nice arrays
  C_docc_   = block_matrix(nso_, nact_docc_);
  C_virt_   = block_matrix(nso_, nact_virt_);
  eps_docc_ = init_array(nact_docc_);
  eps_virt_ = init_array(nact_virt_);
  sym_docc_ = init_int_array(nact_docc_);
  sym_virt_ = init_int_array(nact_virt_);
  int act_docc_count  = 0;
  int act_virt_count  = 0;
  for(int h=0; h<nirrep_; ++h){
    // Skip over the frozen core orbitals in this irrep
    // Copy over the info for active occupied orbitals
    for(int i=0; i<clsdpi_[h]-frzcpi_[h]; ++i){
      for (int mu=0; mu<nso_; ++mu)
        C_docc_[mu][act_docc_count] = C->get(h,mu,i+frzcpi_[h]);;
      eps_docc_[act_docc_count] = epsilon->get(h,i+frzcpi_[h]);
      sym_docc_[act_docc_count] = h;
      ++act_docc_count;
    }
    // Copy over the info for active virtual orbitals
    for(int a=0; a<orbspi_[h]-clsdpi_[h]-frzvpi_[h]; ++a){
      for (int mu=0; mu<nso_; ++mu)
        C_virt_[mu][act_virt_count] = C->get(h,mu,a+clsdpi_[h]);
      eps_virt_[act_virt_count] = epsilon->get(h,a+clsdpi_[h]);
      sym_virt_[act_virt_count] = h;
      ++act_virt_count;
    }
    // Skip over the frozen virtual orbitals in this irrep
  }

  schwarz_ = options_.get_double("SCHWARZ_CUTOFF");
  schwarz_shell_pairs_ = NULL;
  sig_shell_pairs_ = 0;

  print_header();

  //
  //diagnostics
  if (print_>3) {
    fprintf(outfile, "  C_occ:\n");
    print_mat(C_docc_,nso_,nact_docc_,outfile);
    fprintf(outfile, "  C_virt:\n");
    print_mat(C_virt_,nso_,nact_virt_,outfile);
    fprintf(outfile, "  Eps_occ:\n");
    for (int i = 0; i<nact_docc_; i++)
        fprintf(outfile,"  %5d: %14.10f\n", i+1, eps_docc_[i]);
    fprintf(outfile, "  Eps_virt:\n");
    for (int i = 0; i<nact_virt_; i++)
        fprintf(outfile,"  %5d: %14.10f\n", i+1+nact_docc_, eps_virt_[i]);

  }
  timer_off("Setup");

}
void DFMP2::print_header()
{
    int nthreads = 1;
    #ifdef _OPENMP
        nthreads = omp_get_max_threads();
    #endif

    fprintf(outfile, "\t********************************************************\n");
    fprintf(outfile, "\t*                                                      *\n");
    fprintf(outfile, "\t*                        DF-MP2                        *\n");
    fprintf(outfile, "\t*    2nd-Order Density-Fitted Moller-Plesset Theory    *\n");
    fprintf(outfile, "\t*             %4s Algorithm, %3d Threads              *\n",algorithm_type_.c_str(),nthreads);
    fprintf(outfile, "\t*                                                      *\n");
    fprintf(outfile, "\t*    Rob Parrish, Justin Turney, C. David Sherrill,    *\n");
    fprintf(outfile, "\t*        Andy Simmonett, and Edward Hohenstein         *\n");
    fprintf(outfile, "\t*                                                      *\n");
    fprintf(outfile, "\t********************************************************\n");
    fprintf(outfile, "\n\t  ==================================================\n");
    fprintf(outfile, "\t    NSO  NAUX  FOCC  DOCC  AOCC  AVIR  VIRT  FVIR \n");
    fprintf(outfile, "\t   ------------------------------------------------\n");
    fprintf(outfile, "\t   %5d %5d %5d %5d %5d %5d %5d %5d\n",
          nso_,naux_raw_,nf_docc_,ndocc_,nact_docc_,nact_virt_,nvirt_,nf_virt_);
    fprintf(outfile, "\t  ==================================================\n");
    fflush(outfile);
    fflush(outfile);
}
DFMP2::~DFMP2()
{
   free(sym_docc_);
   free(sym_virt_);
   free(eps_docc_);
   free(eps_virt_);
   free(schwarz_shell_pairs_);
   free_block(C_docc_);
   free_block(C_virt_);
}
double DFMP2::compute_energy()
{
    if (algorithm_type_ == "DISK")
        E_ = compute_E_disk();
    else if (algorithm_type_ == "CORE")
        E_ = compute_E_core();
    else if (algorithm_type_ == "OLD")
        E_ = compute_E_old();


    //Finish summing energies
    E_ = E_ss_ + E_os_;
    E_sss_ = ss_scale_*E_ss_;
    E_sos_ = os_scale_*E_os_;
    E_scs_ = E_sss_ + E_sos_;
    E_tot_ = E_ + E_scf_;
    E_tot_scs_ = E_scs_ + E_scf_;

    //Display energies
    fprintf(outfile,"\n   Canonical DF-MP2:\n");
    fprintf(outfile,"    Same-Spin Energy:                %20.16Lf [H]\n",E_ss_);
    fprintf(outfile,"    Opposite-Spin Energy:            %20.16Lf [H]\n",E_os_);
    fprintf(outfile,"    DF-MP2 Correlation Energy:       %20.16f [H]\n",E_);
    fprintf(outfile,"    DF-MP2 Total Energy:             %20.16f [H]\n",E_tot_);
    fprintf(outfile,"   SCS-DF-MP2 (Pss = %14.10f, Pos = %14.10f):\n",ss_scale_,os_scale_);
    fprintf(outfile,"    SCS-Same-Spin Energy:            %20.16f [H]\n",E_sss_);
    fprintf(outfile,"    SCS-Opposite-Spin Energy:        %20.16f [H]\n",E_sos_);
    fprintf(outfile,"    SCS-DF-MP2 Correlation Energy:   %20.16f [H]\n",E_scs_);
    fprintf(outfile,"    SCS-DF-MP2 Total Energy:         %20.16f [H]\n",E_tot_scs_);
    fflush(outfile);

    Process::environment.globals["CURRENT ENERGY"] = E_tot_;
    Process::environment.globals["DF-MP2 ENERGY"] = E_tot_;
    Process::environment.globals["SCS-DF-MP2 ENERGY"] = E_tot_scs_;

    return E_tot_;
}
double DFMP2::compute_E_core()
{
    if (options_.get_bool("PARALLEL_DFMP2")) {
        form_Qia_core_parallel();
        evaluate_contributions_core_parallel();
        free_Qia_core();
    } else {
        form_Qia_core();
        evaluate_contributions_core_sym();
        free_Qia_core();
    }

    return E_;
}
double DFMP2::compute_E_disk()
{
    form_Aia_disk();
    form_Qia_disk();
    evaluate_contributions_disk();

    return E_;
}
void DFMP2::form_Schwarz()
{
    timer_on("Schwarz");
    if (schwarz_ > 0.0) {

        schwarz_shell_pairs_ = init_int_array(basisset_->nshell()*(basisset_->nshell()+1)/2);
        double* max_shell_val = init_array(basisset_->nshell()*(basisset_->nshell()+1)/2);;
        double max_global_val = 0.0;

        IntegralFactory schwarzfactory(basisset_,basisset_,basisset_,basisset_);
        boost::shared_ptr<TwoBodyAOInt> eri = boost::shared_ptr<TwoBodyAOInt>(schwarzfactory.eri());
        const double *buffer = eri->buffer();

        int MU, NU, mu, nu,omu,onu, nummu, numnu, index;
        int MUNU = 0;
        int munu = 0;
        for (MU=0; MU < basisset_->nshell(); ++MU) {
            nummu = basisset_->shell(MU)->nfunction();
            for (NU=0; NU <= MU; ++NU, ++MUNU) {
                numnu = basisset_->shell(NU)->nfunction();
                eri->compute_shell(MU,NU,MU,NU);
                for (mu=0; mu < nummu; ++mu) {
                    omu = basisset_->shell(MU)->function_index() + mu;
                    for (nu=0; nu < numnu; ++nu) {
                        onu = basisset_->shell(NU)->function_index() + nu;

                        if (omu>=onu) {
                            index = mu*(numnu*nummu*numnu+numnu)+nu*(nummu*numnu+1);
                            if (max_global_val<abs(buffer[index]))
                                max_global_val = abs(buffer[index]);
                            if (max_shell_val[MUNU]<abs(buffer[index]))
                                max_shell_val[MUNU] = abs(buffer[index]);
                        }
                    }
                }
            }
        }

        for (int ij = 0; ij < basisset_->nshell()*(basisset_->nshell()+1)/2; ij ++)
            if (max_shell_val[ij]*max_global_val>=schwarz_*schwarz_){
                schwarz_shell_pairs_[ij] = 1;
                sig_shell_pairs_++;
            }

        free(max_shell_val);
    } else {
        schwarz_shell_pairs_ = init_int_array(basisset_->nshell()*(basisset_->nshell()+1)/2);
        for (int ij = 0; ij < basisset_->nshell()*(basisset_->nshell()+1)/2; ij++) {
            schwarz_shell_pairs_[ij] = 1;
            sig_shell_pairs_++;
        }
    }
    timer_off("Schwarz");

}
double** DFMP2::form_W(boost::shared_ptr<BasisSet> bas)
{
  boost::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();

  IntegralFactory rifactory_J(bas, zero, bas, zero);

  int naux = bas->nbf();

  boost::shared_ptr<TwoBodyAOInt> Jint = (boost::shared_ptr<TwoBodyAOInt>)rifactory_J.eri();
  const double *Jbuffer = Jint->buffer();

  double **J = block_matrix(naux, naux);

  int index = 0;

  for (int MU=0; MU < bas->nshell(); ++MU) {
    int nummu = bas->shell(MU)->nfunction();

    for (int NU=0; NU < bas->nshell(); ++NU) {
      int numnu = bas->shell(NU)->nfunction();

      Jint->compute_shell(MU, 0, NU, 0);

      index = 0;
      for (int mu=0; mu < nummu; ++mu) {
        int omu = bas->shell(MU)->function_index() + mu;

        for (int nu=0; nu < numnu; ++nu, ++index) {
          int onu = bas->shell(NU)->function_index() + nu;

          J[omu][onu] = Jbuffer[index];
        }
      }
    }
  }

  return J;
}
double** DFMP2::form_W_overlap(boost::shared_ptr<BasisSet> bas)
{

  boost::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();
  int naux = bas->nbf();

  IntegralFactory rifactory_JS(bas, bas, zero, zero);
  boost::shared_ptr<OneBodyAOInt> S = (boost::shared_ptr<OneBodyAOInt>)rifactory_JS.ao_overlap();

  //Put the integrals in a good old double**
  double** SJ;
  boost::shared_ptr<SimpleMatrix> S_J(new SimpleMatrix("J Overlap", naux, naux));
  //Compute those integrals
  S->compute(S_J);
  //S_J->print(outfile);
  SJ = S_J->to_block_matrix();

  return SJ;
}
int DFMP2::matrix_power(double** J, int n, double power, double max_cond)
{

  double** J1 = block_matrix(n,n);
  C_DCOPY(n*(ULI)n,J[0],1,J1[0],1);

  // the C_DSYEV call replaces the original matrix J with its eigenvectors
  double* eigval = init_array(n);
  int lwork = n * 3;
  double* work = init_array(lwork);
  int stat = C_DSYEV('v','u',n,J1[0],n,eigval,
    work,lwork);
  if (stat != 0) {
    fprintf(outfile, "C_DSYEV failed\n");
    exit(PSI_RETURN_FAILURE);
  }
  free(work);

  double min_J = eigval[0];
  double max_J = eigval[n-1];

  double min_allowed = (max_cond == DBL_MAX ? max_J/max_cond : -DBL_MAX);

  // Now J contains the eigenvectors of the original J
  // Copy J to J_copy
  double **J2 = block_matrix(n, n);
  C_DCOPY(n*(ULI)n,J1[0],1,J2[0],1);

  // Now form J^{-1/2} = U(T)*j^{-1/2}*U,
  // where j^{-1/2} is the diagonal matrix of the inverse square roots
  // of the eigenvalues, and U is the matrix of eigenvectors of J
  int clipped = 0;
  #pragma omp parallel for schedule (static, 1) reduction(+: clipped)
  for (int i=0; i<n; i++) {
    if (eigval[i] < min_allowed) {
      eigval[i] = 0.0;
      clipped++;
    }
    else {
      eigval[i] = pow(eigval[i],power);
    }
    // scale one set of eigenvectors by the diagonal elements j^{-1/2}
    C_DSCAL(n, eigval[i], J1[i], 1);
  }
  free(eigval);

  // J_mhalf = J_copy(T) * J
  C_DGEMM('t','n',n,n,n,1.0,
    J2[0],n,J1[0],n,0.0,J[0],n);

  free_block(J1);
  free_block(J2);

  if (print_>1) {
    fprintf(outfile,"  Smallest eigenvalue in the matrix is %14.10E.\n",min_J);
    fprintf(outfile,"  Largest eigenvalue in the matrix is %14.10E.\n",max_J);
    fprintf(outfile,"  Condition number of the matrix is %14.10E.\n",max_J/min_J);
    fprintf(outfile,"  %d columns removed to maintain a condition number of %14.10E.\n",clipped,max_cond);
    fflush(outfile);
  }
  return clipped;
}
void DFMP2::cholesky_decomposition(double** J, int n)
{
  //Choleskify J (results returned in place, should not be touched)
  int stat = C_DPOTRF('U',n,J[0],n);
  if (stat < 0) {
    fprintf(outfile, "Cholesky argument %d is bad.\n",stat);
    exit(PSI_RETURN_FAILURE);
  }
  else if (stat > 0) {
    fprintf(outfile, "W^+1/2 is not positive definite.\n");
    exit(PSI_RETURN_FAILURE);
  }
}
void DFMP2::cholesky_inverse(double** J, int n)
{
  //Choleskify J (results returned in place, should not be touched)
  int stat = C_DPOTRF('U',n,J[0],n);
  if (stat < 0) {
    fprintf(outfile, "Cholesky argument %d is bad.\n",stat);
    exit(PSI_RETURN_FAILURE);
  }
  else if (stat > 0) {
    fprintf(outfile, "W^+1/2 is not positive definite.\n");
    exit(PSI_RETURN_FAILURE);
  }
  //Inversion (in place)
  int IError = C_DPOTRI('U',n,J[0],n);
  if (IError !=0 )
    throw std::domain_error("J Matrix Inversion Failed!");

  //LAPACK is smart and all, only uses half of the thing
  for (int m = 0; m<n; m++)
    for (int n = 0; n<m; n++)
      J[n][m] = J[m][n];
}
double** DFMP2::form_X(double** W, int n, double max_cond, int* clipped)
{
    int naux_raw = n;
    int naux_fin = n;

    double** SJ = block_matrix(naux_raw,naux_raw);
    C_DCOPY(naux_raw*(ULI)naux_raw,W[0],1,SJ[0],1);

    // the C_DSYEV call replaces the original matrix SJ with its eigenvectors
    double* s_eigval = init_array(naux_raw);
    int lswork = naux_raw * 3;
    double* swork = init_array(lswork);
    int stat = C_DSYEV('v','u',naux_raw,SJ[0],naux_raw,s_eigval, swork,lswork);
    if (stat != 0) {
        fprintf(outfile, "C_DSYEV failed\n");
        exit(PSI_RETURN_FAILURE);
    }
    free(swork);

    //Find the condition number and largest/smallest eigenvalues of SJ
    //The largest one is what is important
    double min_SJ = s_eigval[0];
    double max_SJ = s_eigval[naux_raw-1];
    //Print some stuff
    double cond_SJ = max_SJ/min_SJ;
    if (print_>1) {
        fprintf(outfile,"  Min overlap eigenvalue is %14.10E.\n",min_SJ);
        fprintf(outfile,"  Max overlap eigenvalue is %14.10E.\n",max_SJ);
        fprintf(outfile,"  Overlap condition number is %14.10E.\n",cond_SJ);
        fflush(outfile);
    }

    //Find the user's allowed condition number
    double cond_SJ_tol = max_cond;
    //and compute the smallest allowed eigenvlue in the overlap matrix
    double cutoff_SJ = max_SJ/cond_SJ_tol;

    //Go through and scale columns of eigenvectors by allowed s^-1/2
    int clip = 0;
    double new_min_SJ = s_eigval[naux_raw-1]; //should be the biggest
    for (int i = 0; i<naux_raw; i++) {
        if (s_eigval[i] <= cutoff_SJ) {
            //OK this eigenvalue does not contribute to the subspace completeness
            //and moreover it kills my condition number
            s_eigval[i] = 0.0; //This isn't officially a big deal, but it's proper etiquette
            naux_fin--;
            clip++;
        } else {
            if (new_min_SJ > s_eigval[i])
                new_min_SJ = s_eigval[i];
            //OK, this eigenvalue is approved for launch
            s_eigval[i] = 1.0 / sqrt(s_eigval[i]);
        }
    }

    int ind;
    #pragma omp parallel for private (ind) schedule (dynamic)
    for (int ind = 0; ind<naux_raw; ind++) {

        //Scale the column out to produce the porous X matrix in place
        //(Columns are rows in this messed up LAPACK-ness
        if (s_eigval[ind] > 0.0)
            C_DSCAL(naux_raw, s_eigval[ind], SJ[ind], 1);
    }

    //So how'd that go?
    if (print_>1) {
        fprintf(outfile, "  %d of %d finished auxiliary functions eliminated by canonical orthogonalization.\n",clip, naux_raw);
        fprintf(outfile, "  Resultant condition number in overlap matrix is %14.10E.\n",max_SJ/new_min_SJ);
        fflush(outfile);
    }
    //The famed transformation matrix
    double **X = block_matrix(naux_raw,naux_fin);

    //We'll put it in canonical (descending eigenvalue) order
    #pragma omp parallel for schedule (static, 1)
    for (int m = 0; m<naux_raw; m++)
        for (int n = 0; n<naux_fin; n++)
            X[m][n] = SJ[naux_raw-n-1][m];

    free_block(SJ);
    free(s_eigval);

    *clipped = clip;
    return X;
}
void DFMP2::form_Wm12_fin()
{
    if (print_>1)
        fprintf(outfile,"  Fitting metric conditioned by canonical orthogonalization.\n");

    timer_on("Form W_S");
    double** SJ = form_W_overlap(ribasis_);
    timer_off("Form W_S");

    int clipped = 0;
    double max_cond = options_.get_double("RI_MAX_COND");
    timer_on("Form X");
    double** X = form_X(SJ, naux_raw_,max_cond, &clipped);
    timer_off("Form X");

    naux_fin_ -= clipped;
    free_block(SJ);

    timer_on("Form W");
    double **J = form_W(ribasis_);
    timer_off("Form W");

    timer_on("Form W'");
    //Tansform the J matrix to J'
    double **Jtemp1 = block_matrix(naux_fin_,naux_raw_);
    double **Jp = block_matrix(naux_fin_,naux_fin_);

    //Temp = X'*J
    C_DGEMM('T','N',naux_fin_,naux_raw_,naux_raw_,1.0,X[0],naux_fin_,J[0],naux_raw_,0.0,Jtemp1[0],naux_raw_);
    free_block(J);

    //J' = Temp*X
    C_DGEMM('N','N',naux_fin_,naux_fin_,naux_raw_,1.0,Jtemp1[0],naux_raw_,X[0],naux_fin_,0.0,Jp[0],naux_fin_);
    free_block(Jtemp1);
    timer_off("Form W'");


    //Find J'^-1/2
    timer_on("W'^-1/2");
    matrix_power(Jp, naux_fin_, -1.0/2.0);
    timer_off("W'^-1/2");

    timer_on("XW'^-1/2");
    //W_AB' = X_AC'*J_CB'^-1/2
    W_ = block_matrix(naux_raw_,naux_fin_);

    C_DGEMM('N','N',naux_raw_,naux_fin_,naux_fin_,1.0,X[0],naux_fin_,Jp[0],naux_fin_,0.0,W_[0],naux_fin_);

    free_block(Jp);
    free_block(X);
    timer_off("XW'^-1/2");

    if (debug_) {
        fprintf(outfile,"  W:\n");
        print_mat(W_,naux_raw_,naux_fin_,outfile);
    }
}
void DFMP2::form_Wm12_raw()
{
  if (print_>1)
    fprintf(outfile,"  WARNING: No conditioning algorithm selected for fitting metric.\n");

  timer_on("Form W");
  W_ = form_W(ribasis_);
  timer_off("Form W");

  double max_cond = options_.get_double("RI_MAX_COND");

  timer_on("W^+1/2");
  matrix_power(W_, naux_raw_, -1.0/2.0,max_cond);
  timer_off("W^+1/2");


  if (debug_) {
    fprintf(outfile,"  W:\n");
    print_mat(W_,naux_raw_,naux_fin_,outfile);
  }
}
void DFMP2::form_Wp12_chol()
{
  if (print_>1)
    fprintf(outfile,"  Fitting metric conditioned by Cholesky decomposition.\n");

  timer_on("Form W");
  W_ = form_W(ribasis_);
  timer_off("Form W");

  timer_on("W^+1/2");
  matrix_power(W_, naux_raw_, 1.0/2.0);
  timer_off("W^+1/2");

  timer_on("Chol(W^+1/2)");
  cholesky_decomposition(W_,naux_raw_);
  timer_off("Chol(W^+1/2)");

  if (debug_) {
    fprintf(outfile,"  W:\n");
    print_mat(W_,naux_raw_,naux_fin_,outfile);
  }
}
void DFMP2::form_Aia_disk()
{
  form_Schwarz();

  //Convention
  int norbs =  nso_;
  int nact_docc =  nact_docc_;
  int nact_virt =  nact_virt_;

  //Storage required for each A
  unsigned long int storage_per_row = norbs*norbs+nact_docc*norbs+nact_docc*nact_virt;

  int maxPshell = 0;
  for (int P = 0; P<ribasis_->nshell(); P++)
    if (maxPshell < ribasis_->shell(P)->nfunction())
      maxPshell = ribasis_->shell(P)->nfunction();

  //Max rows of A to pursue
  int max_rows = df_memory_/storage_per_row;
  if (max_rows > naux_raw_)
      max_rows = naux_raw_;
  if (max_rows < maxPshell )
      max_rows = maxPshell;

  //Block sizes
  int nblocks = 2*naux_raw_/max_rows;
  int* p_starts = init_int_array(nblocks);
  int* p_sizes = init_int_array(nblocks);
  int* block_starts = init_int_array(nblocks);
  int* block_stops = init_int_array(nblocks);
  int* block_sizes = init_int_array(nblocks);

  //Determine block sizes
  nblocks = 0;
  int fun_counter = 0;
  for (int A=0; A<ribasis_->nshell(); A++) {
    if (A == ribasis_->nshell() - 1) {
      block_sizes[nblocks]++;
      p_sizes[nblocks] += ribasis_->shell(A)->nfunction();
      block_stops[nblocks] = ribasis_->nshell();
      nblocks++;
      break;
    }
    if (fun_counter+ribasis_->shell(A)->nfunction() > max_rows) {
      block_stops[nblocks] = A;
      block_starts[nblocks+1] = A;
      block_sizes[nblocks+1] = 1;
      p_sizes[nblocks+1] = ribasis_->shell(A)->nfunction();
      p_starts[nblocks+1] = ribasis_->shell(A)->function_index();
      nblocks++;
      fun_counter = ribasis_->shell(A)->nfunction();
      continue;
    }
    p_sizes[nblocks] += ribasis_->shell(A)->nfunction();
    block_sizes[nblocks]++;
    fun_counter += ribasis_->shell(A)->nfunction();
  }

  if (print_>1) {

    fprintf(outfile,"\n  Disk Formation of Aia, using %d blocks.\n",nblocks);
    //for (int b = 0; b<nblocks; b++)
    //  fprintf(outfile,"  block_starts = %d, block_stops  %d, block_sizes = %d, p_start = %d, p_sizes = %d.\n",
    //    block_starts[b], block_stops[b], block_sizes[b], p_starts[b], p_sizes[b]);
  }

  //Buffer tensors
  double** Amn = block_matrix(max_rows,norbs*norbs);
  double** Ami = block_matrix(max_rows,norbs*nact_docc);
  double** Aia = block_matrix(max_rows,nact_docc*nact_virt);

  int nthread = 1;
  #ifdef _OPENMP
  if (options_.get_int("RI_INTS_NUM_THREADS") == 0)
    nthread = omp_get_max_threads();
  else
    nthread = options_.get_int("RI_INTS_NUM_THREADS");
  #endif
  int rank = 0;

  //Get an ERI object for the AO three-index integrals
  IntegralFactory rifactory(ribasis_,zerobasis_, basisset_, basisset_);
  //Get a TEI for each thread
  const double **buffer = new const double*[nthread];
  boost::shared_ptr<TwoBodyAOInt> *eri = new boost::shared_ptr<TwoBodyAOInt>[nthread];
  for (int Q = 0; Q<nthread; Q++) {
    eri[Q] = boost::shared_ptr<TwoBodyAOInt>(rifactory.eri());
    buffer[Q] = eri[Q]->buffer();
  }

  //File handlers
  //Fast index is ia, not A, as less storage is required in the next step
  psio_->open(PSIF_DFMP2_AIA,PSIO_OPEN_NEW);
  psio_address next_PSIF_DFMP2_AIA = PSIO_ZERO;

  //indices for three-index integrals
  int P, MU, NU, nump, nummu, numnu, p, mu, nu, op, omu, onu, index;

  //Loop over blocks of A
  for (int block = 0; block<nblocks; block++) {
    //Zero that guy out!
    memset((void*)&Amn[0][0],'\0',p_sizes[block]*norbs*(ULI)norbs);

    //Form Amn ints
    timer_on("(A|mn)");
    #pragma omp parallel for private (P, MU, NU, p, mu, nu, nump, nummu, numnu, op, omu, onu, index, rank) schedule (dynamic) num_threads(nthread)
    for (P=block_starts[block]; P < block_stops[block]; ++P) {
      #ifdef _OPENMP
         rank = omp_get_thread_num();
      #endif
      nump = ribasis_->shell(P)->nfunction();
      for (MU=0; MU < basisset_->nshell(); ++MU) {
        nummu = basisset_->shell(MU)->nfunction();
        for (NU=0; NU <= MU; ++NU) {
          numnu = basisset_->shell(NU)->nfunction();
          if (schwarz_shell_pairs_[MU*(MU+1)/2+NU] == 1) {
            eri[rank]->compute_shell(P, 0, MU, NU);
            for (p=0, index=0; p < nump; ++p) {
              op = ribasis_->shell(P)->function_index() + p;
              for (mu=0; mu < nummu; ++mu) {
                omu = basisset_->shell(MU)->function_index() + mu;
                for (nu=0; nu < numnu; ++nu, ++index) {
                  onu = basisset_->shell(NU)->function_index() + nu;
                  Amn[op-p_starts[block]][omu*norbs+onu] = buffer[rank][index]; // (op | omu onu) integral
                  Amn[op-p_starts[block]][onu*norbs+omu] = buffer[rank][index]; // (op | onu omu) integral
                }
              }
            }
          }
        }
      }
    }
    timer_off("(A|mn)");

    //fprintf(outfile, "  Amn\n");
    //print_mat(Amn,max_rows,norbs*norbs, outfile);

    //Transform to Ami
    // (A|mi) = (Amn)C_ni
    timer_on("(A|mi)");
    C_DGEMM('N', 'N', p_sizes[block]*norbs, nact_docc, norbs, 1.0, &(Amn[0][0]),
        norbs, &(C_docc_[0][0]), nact_docc, 0.0, &(Ami[0][0]), nact_docc);
    timer_off("(A|mi)");


    //fprintf(outfile, "  Ami\n");
    //print_mat(Ami,max_rows,nact_docc*norbs, outfile);

    #ifdef _MKL
       int mkl_nthreads = mkl_get_max_threads();
       mkl_set_num_threads(1);
    #endif

    //Transform to Aia
    // (A|ia) = C_ma(A|mi)
    timer_on("(A|ia)");
    #pragma omp parallel for
    for (int A = 0; A<p_sizes[block]; A++) {
      C_DGEMM('T', 'N', nact_docc, nact_virt, norbs, 1.0, &(Ami[A][0]),
        nact_docc, &(C_virt_[0][0]), nact_virt, 0.0, &(Aia[A][0]), nact_virt);
    }
    timer_off("(A|ia)");

    #ifdef _MKL
       mkl_set_num_threads(mkl_nthreads);
    #endif

    //fprintf(outfile, "  Aia\n");
    //print_mat(Aia,max_rows,nact_docc*nact_virt, outfile);

    //Stripe to disk
    timer_on("(A|ia) Write");
    psio_->write(PSIF_DFMP2_AIA,"Aia Three-Index Integrals",(char *)(&Aia[0][0]),sizeof(double)*p_sizes[block]*nact_docc*(ULI)nact_virt,next_PSIF_DFMP2_AIA,&next_PSIF_DFMP2_AIA);
    timer_off("(A|ia) Write");

  }

  psio_->close(PSIF_DFMP2_AIA,1);
  delete[] buffer;
  delete[]eri;
  free(p_starts);
  free(p_sizes);
  free(block_starts);
  free(block_stops);
  free(block_sizes);
  free_block(Amn);
  free_block(Ami);
  free_block(Aia);
}
void DFMP2::form_Qia_disk()
{
  //Form fitting metric
  if (options_.get_str("RI_FITTING_TYPE") == "RAW")
      form_Wm12_raw();
  else if (options_.get_str("RI_FITTING_TYPE") == "FINISHED")
      form_Wm12_fin();
  else if (options_.get_str("RI_FITTING_TYPE") == "CHOLESKY")
      form_Wp12_chol();

  //Available memory is lower due to W
  unsigned long int available_memory = df_memory_-naux_raw_*naux_fin_;

  //max_cols to attempt
  int max_cols = available_memory/((long) naux_raw_+naux_fin_);
  if (max_cols > nact_virt_*nact_docc_)
      max_cols = nact_virt_*nact_docc_;
  if (max_cols < 1)
      max_cols = 1;

  //How many blocks are there?
  int nblocks = nact_docc_*nact_virt_/max_cols;
  if (max_cols*nblocks < nact_docc_*nact_virt_)
    nblocks++;

  //File handlers
  psio_->open(PSIF_DFMP2_AIA,PSIO_OPEN_OLD);
  psio_address next_PSIF_DFMP2_AIA = PSIO_ZERO;
  psio_->open(PSIF_DFMP2_QIA,PSIO_OPEN_NEW);
  psio_address next_PSIF_DFMP2_QIA = PSIO_ZERO;

  //Prestripe
  timer_on("(Q|ia) Prestripe");
  double *Prestripe = init_array(nact_docc_*nact_virt_);
  for (int Q = 0; Q < naux_fin_; Q++) {
       psio_->write(PSIF_DFMP2_QIA,"Qia Three-Index Integrals",(char *) &(Prestripe[0]),sizeof(double)*nact_docc_*nact_virt_,next_PSIF_DFMP2_QIA,&next_PSIF_DFMP2_QIA);
  }
  timer_off("(Q|ia) Prestripe");
  free(Prestripe);
  next_PSIF_DFMP2_QIA = PSIO_ZERO;

  //Allocate buffers
  double** Aia = block_matrix(naux_raw_, max_cols);
  double** Qia = block_matrix(max_cols, naux_fin_);

  int current_column = 0;
  int current_columns;
  for (int block = 0; block<nblocks; block++) {

    //How many columns in this pass?
    current_columns = max_cols;
    if (current_column+current_columns>=nact_docc_*nact_virt_)
        current_columns = nact_docc_*nact_virt_-current_column;

    //Read Aia in
    timer_on("(A|ia) Read");
    for (int A = 0; A<naux_raw_; A++) {
        next_PSIF_DFMP2_AIA = psio_get_address(PSIO_ZERO,(ULI)(A*(ULI)nact_docc_*nact_virt_*sizeof(double)+current_column*sizeof(double)));
        psio_->read(PSIF_DFMP2_AIA,"Aia Three-Index Integrals",(char *)(&Aia[A][0]),sizeof(double)*(ULI)current_columns,next_PSIF_DFMP2_AIA,&next_PSIF_DFMP2_AIA);
    }
    timer_off("(A|ia) Read");

    //Embed fitting
    timer_on("(Q|ia)");
    if (options_.get_str("RI_FITTING_TYPE") == "CHOLESKY") {

        #pragma omp parallel for schedule(static)
        for (int A = 0; A< naux_raw_; A++) {
            C_DCOPY(current_columns, &Aia[A][0], 1, &Qia[0][A], naux_raw_);
        }

        int info = C_DPOTRS('U',naux_raw_,current_columns,&W_[0][0],naux_raw_,&Qia[0][0],naux_raw_);

    } else {
        C_DGEMM('T','N',current_columns,naux_fin_,naux_raw_,1.0, Aia[0], max_cols, W_[0], naux_fin_,0.0, Qia[0],naux_fin_);
    }
    timer_off("(Q|ia)");

    //fprintf(outfile,"  Nblocks = %d, Max cols = %d, current_columns = %d, current_column = %d\n",nblocks, max_cols, current_columns, current_column);

    //Write Qia out
    timer_on("(Q|ia) Write");
    psio_->write(PSIF_DFMP2_QIA,"Qia Three-Index Integrals",(char *)(&Qia[0][0]),sizeof(double)*naux_fin_*(ULI)current_columns,next_PSIF_DFMP2_QIA,&next_PSIF_DFMP2_QIA);
    timer_off("(Q|ia) Write");

    // Increment column offset
    current_column += current_columns;

  }
  psio_->close(PSIF_DFMP2_AIA,1);
  psio_->close(PSIF_DFMP2_QIA,1);

  //Free fitting metric
  free_block(W_);
  free_block(Aia);
  free_block(Qia);
}
void DFMP2::evaluate_contributions_disk()
{
    unsigned long int available_memory = df_memory_;

    //How many occupieds per block
    unsigned long int N_I = available_memory/(2*nact_virt_*(ULI)naux_fin_);
    if (N_I > nact_docc_)
        N_I = nact_docc_;
    if (N_I < 1)
        N_I = 1;

    //How many blocks
    unsigned long int nblocks = nact_docc_/N_I;
    if (nact_docc_ > N_I*nblocks)
        nblocks++;

    //Block indexing
    int* block_starts = init_int_array(nblocks);
    int* block_sizes = init_int_array(nblocks);
    int* block_stops = init_int_array(nblocks);

    for (int block = 0; block < nblocks; block++) {
        block_starts[block] = N_I*block;
        if (block == nblocks - 1)
            block_sizes[block] = nact_docc_-block_starts[block];
        else
            block_sizes[block] = N_I;
    }

    for (int block = 0; block < nblocks; block++) {
        block_stops[block] = block_starts[block] + block_sizes[block];
    }

    if (debug_) {
        fprintf(outfile, "  Nblocks = %ld, N_I = %ld\n", nblocks, N_I);
        for (int block = 0; block < nblocks; block++)
            fprintf(outfile,"  Block %d, start = %d, size = %d, stop = %d\n",block+1, block_starts[block], block_sizes[block], block_stops[block]);
        fflush(outfile);
    }

    //Integral buffers
    double **Qia = block_matrix(N_I*(ULI)nact_virt_, naux_fin_);
    double **Qjb = block_matrix(N_I*(ULI)nact_virt_, naux_fin_);

    //Threads
    int nthreads = 1;
    #ifdef _OPENMP
        nthreads = omp_get_max_threads();
    #endif

    // MKL trick
    #ifdef _MKL
        int mkl_nthreads = mkl_get_max_threads();
        mkl_set_num_threads(1);
    #endif

    //zipped integrals buffer
    double*** I = new double**[nthreads];
    for (int r = 0 ; r<nthreads; r++) {
        I[r] = block_matrix(nact_virt_,nact_virt_);
    }

    //File handler to fitted integrals
    psio_->open(PSIF_DFMP2_QIA,PSIO_OPEN_OLD);
    psio_address next_PSIF_DFMP2_QIA = PSIO_ZERO;

    //Loop over blocks
    for (int block1 = 0; block1 < nblocks; block1++) {

        //Startup the first block
        if (block1 == 0) {
            timer_on("(Q|ia) Read");
            next_PSIF_DFMP2_QIA = psio_get_address(PSIO_ZERO,(ULI)(block_starts[block1]*nact_virt_*naux_fin_*(ULI)sizeof(double)));
            psio_->read(PSIF_DFMP2_QIA,"Qia Three-Index Integrals",(char *)(&Qjb[0][0]),sizeof(double)*naux_fin_*nact_virt_*(ULI)block_sizes[block1],next_PSIF_DFMP2_QIA,&next_PSIF_DFMP2_QIA);
            timer_off("(Q|ia) Read");
        }

        //Copy Qjb up to Qia
        C_DCOPY(naux_fin_*nact_virt_*(ULI)block_sizes[block1],Qjb[0],1,Qia[0],1);

        //Same block contribution (ai|Q)(Q|ai)
        find_disk_contributions(Qia,Qjb,I,block_starts,block_sizes,block_stops,block1,block1);

        //Different block contribution (ai|Q)(Q|bj)
        for (int block2 = nblocks-1; block2 > block1; block2--) {
            //Read in the second block
            timer_on("(Q|ia) Read");
            next_PSIF_DFMP2_QIA = psio_get_address(PSIO_ZERO,(ULI)(block_starts[block2]*nact_virt_*naux_fin_*(ULI)sizeof(double)));
            psio_->read(PSIF_DFMP2_QIA,"Qia Three-Index Integrals",(char *)(&Qjb[0][0]),sizeof(double)*naux_fin_*nact_virt_*(ULI)block_sizes[block2],next_PSIF_DFMP2_QIA,&next_PSIF_DFMP2_QIA);
            timer_off("(Q|ia) Read");

            //Find the contributions
            find_disk_contributions(Qia,Qjb,I,block_starts,block_sizes,block_stops,block1,block2);

        }
    }

    //Close
    psio_->close(PSIF_DFMP2_QIA,1);

    //Frees
    for (int r = 0; r<nthreads; r++)
        free_block(I[r]);
    delete[] I;

    free_block(Qia);
    free_block(Qjb);
    free(block_starts);
    free(block_stops);
    free(block_sizes);

    //Reset MKL
    #ifdef _MKL
        mkl_set_num_threads(mkl_nthreads);
    #endif

}
void DFMP2::find_disk_contributions(double** Qia, double** Qjb, double***I, int* starts, int* sizes, int* stops, int block1, int block2)
{
  int rank = 0;
  double iajb, ibja, denom, perm_scale;
  long double e_ss, e_os;
  int a, b, i, j;

  e_ss = 0.0;
  e_os = 0.0;

  #pragma omp parallel for private (rank, a, b, i, j, perm_scale, iajb, ibja, denom) reduction (+: e_ss, e_os) schedule (dynamic)
  for (j=starts[block2]; j < stops[block2]; ++j) {

    #ifdef _OPENMP
      rank = omp_get_thread_num();
    #endif

    for (i=starts[block1]; i < stops[block1] && i<= j ; ++i) {

      if (rank == 0) timer_on("(ia|jb)");
      C_DGEMM('N','T',nact_virt_,nact_virt_,naux_fin_,1.0,&(Qia[(i-starts[block1])*nact_virt_][0]),
        naux_fin_,&(Qjb[(j-starts[block2])*nact_virt_][0]),naux_fin_,0.0,&(I[rank][0][0]),nact_virt_);
      if (rank == 0) timer_off("(ia|jb)");

      //fprintf(outfile, "  ij = (%d,%d), I:\n",i,j);
      //print_mat(I[0],nact_virt_,nact_virt_,outfile);

      if (rank == 0) timer_on("T_MP2");
      for (a=0; a < nact_virt_; ++a) {
        for (b=0; b < nact_virt_; ++b) {
          iajb = I[rank][a][b];
          ibja = I[rank][b][a];
          denom = 1.0 /
            (eps_docc_[i] + eps_docc_[j] -
             eps_virt_[a] - eps_virt_[b]);

          perm_scale = ((i == j)? 1.0 : 2.0 );

          e_ss += perm_scale*(iajb*iajb-iajb*ibja)*denom;
          e_os += perm_scale*(iajb*iajb)*denom;
        }
      }
      if (rank == 0) timer_off("T_MP2");
    }
  }

  if (debug_)
    fprintf(outfile,"  Block1 = %d, Block2 = %d, e_ss = %20Lf, e_os = %20Lf\n",block1,block2, e_ss, e_os);

  E_ss_ += e_ss;
  E_os_ += e_os;

}
double** DFMP2::form_Aia_core()
{
  //Convention
  int norbs =  nso_;
  int nact_docc =  nact_docc_;
  int nact_virt =  nact_virt_;

  //Setup Schwarz sieve
  form_Schwarz();

  //How much memory does Aia take?
  unsigned long int three_mem = naux_raw_*nact_docc_*(ULI)nact_virt_;
  //How much memory is left for buffers
  unsigned long int available_mem = df_memory_-three_mem;

  //Setup blocks
  //Storage required for each A
  unsigned long int storage_per_row = norbs*(ULI)norbs+nact_docc*(ULI)norbs;

  int maxPshell = 0;
  for (int P = 0; P<ribasis_->nshell(); P++)
    if (maxPshell < ribasis_->shell(P)->nfunction())
      maxPshell = ribasis_->shell(P)->nfunction();

  //Max rows of A to pursue
  int max_rows = available_mem/storage_per_row;
  if (max_rows > naux_raw_ + maxPshell)
      max_rows = naux_raw_ + maxPshell;
  if (max_rows < maxPshell )
      max_rows = maxPshell;

  //Block sizes
  unsigned long int nblocks = 2*naux_raw_/max_rows;
  int* p_starts = init_int_array(nblocks);
  int* p_sizes = init_int_array(nblocks);
  int* block_starts = init_int_array(nblocks);
  int* block_stops = init_int_array(nblocks);
  int* block_sizes = init_int_array(nblocks);

  // How many functions per shell?
  //for (int P = 0; P < ribasis_->nshell(); P++) {
  //  printf("RI shell %d has %d functions\n", P, ribasis_->shell(P)->nfunction());
  //}

  //Determine block sizes
  nblocks = 0;
  int fun_counter = 0;
  for (int A=0; A<ribasis_->nshell(); A++) {
    if (A == ribasis_->nshell() - 1) {
      block_sizes[nblocks]++;
      p_sizes[nblocks] += ribasis_->shell(A)->nfunction();
      block_stops[nblocks] = ribasis_->nshell();
      nblocks++;
      break;
    }
    if (fun_counter+ribasis_->shell(A)->nfunction() > max_rows) {
      block_stops[nblocks] = A;
      block_starts[nblocks+1] = A;
      block_sizes[nblocks+1] = 1;
      p_sizes[nblocks+1] = ribasis_->shell(A)->nfunction();
      p_starts[nblocks+1] = ribasis_->shell(A)->function_index();
      nblocks++;
      fun_counter = ribasis_->shell(A)->nfunction();
      continue;
    }
    p_sizes[nblocks] += ribasis_->shell(A)->nfunction();
    block_sizes[nblocks]++;
    fun_counter += ribasis_->shell(A)->nfunction();
  }

  if (print_>1) {

    fprintf(outfile,"\n  Core Formation of Aia, using %ld blocks.\n",nblocks);
    if (debug_) {
      for (int b = 0; b<nblocks; b++)
        fprintf(outfile,"  block_starts = %d, block_stops  %d, block_sizes = %d, p_start = %d, p_sizes = %d.\n",
          block_starts[b], block_stops[b], block_sizes[b], p_starts[b], p_sizes[b]);
    }
  }

  double** Aia = block_matrix(naux_raw_,nact_docc_*(ULI)nact_virt_);

  //Buffer tensors
  double** Amn = block_matrix(max_rows,norbs*(ULI)norbs);
  double** Ami = block_matrix(max_rows,norbs*(ULI)nact_docc);

  int nthread = 1;
  #ifdef _OPENMP
  if (options_.get_int("RI_INTS_NUM_THREADS") == 0)
    nthread = omp_get_max_threads();
  else
    nthread = options_.get_int("RI_INTS_NUM_THREADS");
  #endif
  int rank = 0;

  //Get an ERI object for the AO three-index integrals
  IntegralFactory rifactory(ribasis_,zerobasis_, basisset_, basisset_);
  //Get a TEI for each thread
  const double **buffer = new const double*[nthread];
  boost::shared_ptr<TwoBodyAOInt> *eri = new boost::shared_ptr<TwoBodyAOInt>[nthread];
  for (int Q = 0; Q<nthread; Q++) {
    eri[Q] = boost::shared_ptr<TwoBodyAOInt>(rifactory.eri());
    buffer[Q] = eri[Q]->buffer();
  }

  //indices for three-index integrals
  int P, MU, NU, nump, nummu, numnu, p, mu, nu, op, omu, onu, index;

  //Loop over blocks of A
  for (int block = 0; block<nblocks; block++) {
    //Zero that guy out!
    timer_on("(A|mn)");
    memset((void*)&Amn[0][0],'\0',block_sizes[block]*norbs*norbs);

    //Form Amn ints
    #pragma omp parallel for private (P, MU, NU, p, mu, nu, nump, nummu, numnu, op, omu, onu, index, rank) schedule (dynamic) num_threads(nthread)
    for (P=block_starts[block]; P < block_stops[block]; ++P) {
      #ifdef _OPENMP
         rank = omp_get_thread_num();
      #endif

      nump = ribasis_->shell(P)->nfunction();
      for (MU=0; MU < basisset_->nshell(); ++MU) {
        nummu = basisset_->shell(MU)->nfunction();
        for (NU=0; NU <= MU; ++NU) {
          numnu = basisset_->shell(NU)->nfunction();
          if (schwarz_shell_pairs_[MU*(MU+1)/2+NU] == 1) {
            eri[rank]->compute_shell(P, 0, MU, NU);
            for (p=0, index=0; p < nump; ++p) {
              op = ribasis_->shell(P)->function_index() + p;
              for (mu=0; mu < nummu; ++mu) {
                omu = basisset_->shell(MU)->function_index() + mu;
                for (nu=0; nu < numnu; ++nu, ++index) {
                  onu = basisset_->shell(NU)->function_index() + nu;
                  Amn[op-p_starts[block]][omu*norbs+onu] = buffer[rank][index]; // (op | omu onu) integral
                  Amn[op-p_starts[block]][onu*norbs+omu] = buffer[rank][index]; // (op | onu omu) integral
                }
              }
            }
          }
        }
      }
    }

    if (debug_) {
      fprintf(outfile, "  Amn\n");
      print_mat(Amn,max_rows,norbs*norbs, outfile);
    }
    timer_off("(A|mn)");

    //Transform to Ami
    // (A|mi) = (Amn)C_ni
    timer_on("(A|mi)");
    C_DGEMM('N', 'N', p_sizes[block]*norbs, nact_docc, norbs, 1.0, &(Amn[0][0]),
        norbs, &(C_docc_[0][0]), nact_docc, 0.0, &(Ami[0][0]), nact_docc);
    timer_off("(A|mi)");

    if (debug_) {
      fprintf(outfile, "  Ami\n");
      print_mat(Ami,max_rows,nact_docc*norbs, outfile);
    }

    timer_on("(A|ia)");
    #ifdef _MKL
       int mkl_nthreads = mkl_get_max_threads();
       mkl_set_num_threads(1);
    #endif

    //Transform to Aia
    // (A|ia) = C_ma(A|mi)
    #pragma omp parallel for
    for (int A = 0; A<p_sizes[block]; A++) {
      C_DGEMM('T', 'N', nact_docc, nact_virt, norbs, 1.0, &(Ami[A][0]),
        nact_docc, &(C_virt_[0][0]), nact_virt, 0.0, &(Aia[p_starts[block]+A][0]), nact_virt_);
    }

    #ifdef _MKL
       mkl_set_num_threads(mkl_nthreads);
    #endif
    timer_off("(A|ia)");

  }
  if (debug_) {
    fprintf(outfile, "  Aia\n");
    print_mat(Aia,naux_raw_,nact_docc*nact_virt, outfile);
  }

  delete[] buffer;
  delete[] eri;
  free(p_starts);
  free(p_sizes);
  free(block_starts);
  free(block_stops);
  free(block_sizes);
  free_block(Amn);
  free_block(Ami);

  return Aia;
}
void DFMP2::form_Qia_core()
{

  //Convention
  int norbs =  nso_;
  int nact_docc =  nact_docc_;
  int nact_virt =  nact_virt_;
  //printf("nact_docc_ is %d in form_Qia_core()\n", nact_docc_);
  Qia_ = form_Aia_core();

  //Setup fitting tensor
  if (options_.get_str("RI_FITTING_TYPE") == "RAW")
      form_Wm12_raw();
  else if (options_.get_str("RI_FITTING_TYPE") == "FINISHED")
      form_Wm12_fin();
  else if (options_.get_str("RI_FITTING_TYPE") == "CHOLESKY")
      form_Wp12_chol();

  unsigned long int three_mem = naux_raw_*nact_docc_*(ULI)nact_virt_;
  unsigned long int available_mem = df_memory_-three_mem-naux_raw_*(ULI)naux_fin_;
  unsigned long int max_cols = available_mem/(ULI)naux_raw_;
  if (max_cols > nact_docc_*(ULI)nact_virt_)
    max_cols = nact_docc_*(ULI)nact_virt;
  if (max_cols < 1)
    max_cols = 1;

  unsigned long int nblocks = nact_docc_*(ULI)nact_virt_/max_cols;
  if (nblocks*max_cols < nact_docc_*nact_virt_)
    nblocks++;

  int* column_starts = init_int_array(nblocks);
  int* column_sizes = init_int_array(nblocks);

  for (int block = 0; block < nblocks; block++) {
    if (block == nblocks - 1) {
      column_sizes[block] = nact_docc_*(ULI)nact_virt_ - column_starts[block];
      break;
    }
    column_sizes[block] = max_cols;
    column_starts[block+1] = column_starts[block] + column_sizes[block];
  }

  //Apply fitting by blocks
  double** Abuffer;
  if (options_.get_str("RI_FITTING_TYPE") == "CHOLESKY") {
      Abuffer = block_matrix(max_cols, naux_raw_);
  } else {
      Abuffer = block_matrix(naux_raw_, max_cols);
  }
  timer_on("(Q|ia)");
  for (int block = 0; block < nblocks; block++) {

     if (options_.get_str("RI_FITTING_TYPE") == "CHOLESKY") {

      #pragma omp parallel for schedule(static)
      for (int A = 0; A< naux_raw_; A++) {
        C_DCOPY(column_sizes[block], &Qia_[A][column_starts[block]], 1, &Abuffer[0][A], naux_raw_);
      }

      int info = C_DPOTRS('U',naux_raw_,column_sizes[block],&W_[0][0],naux_raw_,&Abuffer[0][0],naux_raw_);

      #pragma omp parallel for schedule(static)
      for (int A = 0; A< naux_raw_; A++) {
        C_DCOPY(column_sizes[block], &Abuffer[0][A], naux_raw_, &Qia_[A][column_starts[block]], 1);
      }
    } else {
      #pragma omp parallel for schedule(static)
      for (int A = 0; A< naux_raw_; A++) {
        C_DCOPY(column_sizes[block], &Qia_[A][column_starts[block]], 1, &Abuffer[A][0], 1);
      }

      //Multiply back into Qia_ with tight striping
      C_DGEMM('T','N',naux_fin_, column_sizes[block], naux_raw_, 1.0, W_[0], naux_fin_, Abuffer[0], max_cols, 0.0, &Qia_[0][column_starts[block]], nact_docc_*(ULI) nact_virt_);
    }
  }
  timer_off("(Q|ia)");

  if (debug_) {
    fprintf(outfile,"  Qia");
    print_mat(Qia_,naux_fin_,nact_docc_*nact_virt_,outfile);
  }

  free_block(W_);
  free_block(Abuffer);
  free(column_sizes);
  free(column_starts);
}
void DFMP2::evaluate_contributions_core_sym()
{
  int nthreads = 1;
  int rank = 0;

  #ifdef _OPENMP
    nthreads = omp_get_max_threads();
  #endif

  #ifdef _MKL
    int mkl_nthreads = mkl_get_max_threads();
    mkl_set_num_threads(1);
  #endif

  double*** I = new double**[nthreads];
  for (int r = 0 ; r<nthreads; r++) {
    I[r] = block_matrix(nact_virt_,nact_virt_);
  }

  double iajb, ibja, denom, perm_scale;
  long double e_ss, e_os;
  int a, b, i, j;

  e_ss = 0.0;
  e_os = 0.0;

  #pragma omp parallel for private (rank, a, b, i, j, perm_scale, iajb, ibja, denom) reduction (+: e_ss, e_os) schedule (dynamic)
  for (i=0; i < nact_docc_; ++i) {

    #ifdef _OPENMP
      rank = omp_get_thread_num();
    #endif

    for (j=0; j <= i; ++j) {

      if (rank == 0) timer_on("(ia|jb)");
      C_DGEMM('T','N',nact_virt_,nact_virt_,naux_fin_,1.0,&(Qia_[0][i*nact_virt_]),
        nact_docc_*(ULI) nact_virt_,&(Qia_[0][j*nact_virt_]),nact_docc_*(ULI) nact_virt_,0.0,&(I[rank][0][0]),nact_virt_);
      if (rank == 0) timer_off("(ia|jb)");

      //fprintf(outfile, "  ij = (%d,%d), I:\n",i,j);
      //print_mat(I[0],nact_virt_,nact_virt_,outfile);

      if (rank == 0) timer_on("T_MP2");
      for (a=0; a < nact_virt_; ++a) {
        for (b=0; b < nact_virt_; ++b) {
          iajb = I[rank][a][b];
          ibja = I[rank][b][a];
          denom = 1.0 /
            (eps_docc_[i] + eps_docc_[j] -
             eps_virt_[a] - eps_virt_[b]);

          perm_scale = ((i == j)? 1.0 : 2.0 );

          e_ss += perm_scale*(iajb*iajb-iajb*ibja)*denom;
          e_os += perm_scale*(iajb*iajb)*denom;
        }
      }
      if (rank == 0) timer_off("T_MP2");
    }
  }

  for (int r = 0; r<nthreads; r++)
    free_block(I[r]);
  delete[] I;

  #ifdef _MKL
    mkl_set_num_threads(mkl_nthreads);
  #endif

  E_ss_ += e_ss;
  E_os_ += e_os;
}
double** DFMP2::form_Aia_core_parallel()
{
  //Convention
  int norbs =  nso_;
  int nact_docc =  nact_docc_;
  int nact_virt =  nact_virt_;

  //Setup Schwarz sieve
  form_Schwarz();

  //How much memory does Aia take?
  unsigned long int three_mem = naux_raw_*nact_docc_*(ULI)nact_virt_;
  //How much memory is left for buffers
  unsigned long int available_mem = df_memory_-three_mem;

  //Setup blocks
  //Storage required for each A
  unsigned long int storage_per_row = norbs*(ULI)norbs+nact_docc*(ULI)norbs;

  int maxPshell = 0;
  for (int P = 0; P<ribasis_->nshell(); P++)
    if (maxPshell < ribasis_->shell(P)->nfunction())
      maxPshell = ribasis_->shell(P)->nfunction();

  //Max rows of A to pursue
  int max_rows = available_mem/storage_per_row;
  if (max_rows > naux_raw_ + maxPshell)
      max_rows = naux_raw_ + maxPshell;
  if (max_rows < maxPshell )
      max_rows = maxPshell;

  //Block sizes
  unsigned long int nblocks = 2*naux_raw_/max_rows;
  int* p_starts = init_int_array(nblocks);
  int* p_sizes = init_int_array(nblocks);
  int* block_starts = init_int_array(nblocks);
  int* block_stops = init_int_array(nblocks);
  int* block_sizes = init_int_array(nblocks);

  // How many functions per shell?
  for (int P = 0; P < ribasis_->nshell(); P++) {
    printf("RI shell %d has %d functions\n", P, ribasis_->shell(P)->nfunction());
  }

  //Determine block sizes
  nblocks = 0;
  int fun_counter = 0;
  for (int A=0; A<ribasis_->nshell(); A++) {
    if (A == ribasis_->nshell() - 1) {
      block_sizes[nblocks]++;
      p_sizes[nblocks] += ribasis_->shell(A)->nfunction();
      block_stops[nblocks] = ribasis_->nshell();
      nblocks++;
      break;
    }
    if (fun_counter+ribasis_->shell(A)->nfunction() > max_rows) {
      block_stops[nblocks] = A;
      block_starts[nblocks+1] = A;
      block_sizes[nblocks+1] = 1;
      p_sizes[nblocks+1] = ribasis_->shell(A)->nfunction();
      p_starts[nblocks+1] = ribasis_->shell(A)->function_index();
      nblocks++;
      fun_counter = ribasis_->shell(A)->nfunction();
      continue;
    }
    p_sizes[nblocks] += ribasis_->shell(A)->nfunction();
    block_sizes[nblocks]++;
    fun_counter += ribasis_->shell(A)->nfunction();
  }

  if (print_>1) {

    fprintf(outfile,"\n  Core Formation of Aia, using %ld blocks.\n",nblocks);
    if (debug_) {
      for (int b = 0; b<nblocks; b++)
        fprintf(outfile,"  block_starts = %d, block_stops  %d, block_sizes = %d, p_start = %d, p_sizes = %d.\n",
          block_starts[b], block_stops[b], block_sizes[b], p_starts[b], p_sizes[b]);
    }
  }

  //Final tensor (will be fitted later)
  double** Aia = block_matrix(naux_raw_,nact_docc_*(ULI)nact_virt_);

  //Buffer tensors
  double** Amn = block_matrix(max_rows,norbs*(ULI)norbs);
  double** Ami = block_matrix(max_rows,norbs*(ULI)nact_docc);

  int nthread = 1;
  #ifdef _OPENMP
  if (options_.get_int("RI_INTS_NUM_THREADS") == 0)
    nthread = omp_get_max_threads();
  else
    nthread = options_.get_int("RI_INTS_NUM_THREADS");
  #endif
  int rank = 0;

  //Get an ERI object for the AO three-index integrals
  IntegralFactory rifactory(ribasis_,zerobasis_, basisset_, basisset_);
  //Get a TEI for each thread
  const double **buffer = new const double*[nthread];
  boost::shared_ptr<TwoBodyAOInt> *eri = new boost::shared_ptr<TwoBodyAOInt>[nthread];
  for (int Q = 0; Q<nthread; Q++) {
    eri[Q] = boost::shared_ptr<TwoBodyAOInt>(rifactory.eri());
    buffer[Q] = eri[Q]->buffer();
  }

  //indices for three-index integrals
  int P, MU, NU, nump, nummu, numnu, p, mu, nu, op, omu, onu, index;

  //Loop over blocks of A
  for (int block = 0; block<nblocks; block++) {
    //Zero that guy out!
    timer_on("(A|mn)");
    memset((void*)&Amn[0][0],'\0',block_sizes[block]*norbs*norbs);

    //Form Amn ints
    #pragma omp parallel for private (P, MU, NU, p, mu, nu, nump, nummu, numnu, op, omu, onu, index, rank) schedule (dynamic) num_threads(nthread)
    for (P=block_starts[block]; P < block_stops[block]; ++P) {
      #ifdef _OPENMP
         rank = omp_get_thread_num();
      #endif

      nump = ribasis_->shell(P)->nfunction();
      for (MU=0; MU < basisset_->nshell(); ++MU) {
        nummu = basisset_->shell(MU)->nfunction();
        for (NU=0; NU <= MU; ++NU) {
          numnu = basisset_->shell(NU)->nfunction();
          if (schwarz_shell_pairs_[MU*(MU+1)/2+NU] == 1) {
            eri[rank]->compute_shell(P, 0, MU, NU);
            for (p=0, index=0; p < nump; ++p) {
              op = ribasis_->shell(P)->function_index() + p;
              for (mu=0; mu < nummu; ++mu) {
                omu = basisset_->shell(MU)->function_index() + mu;
                for (nu=0; nu < numnu; ++nu, ++index) {
                  onu = basisset_->shell(NU)->function_index() + nu;
                  Amn[op-p_starts[block]][omu*norbs+onu] = buffer[rank][index]; // (op | omu onu) integral
                  Amn[op-p_starts[block]][onu*norbs+omu] = buffer[rank][index]; // (op | onu omu) integral
                }
              }
            }
          }
        }
      }
    }

    if (debug_) {
      fprintf(outfile, "  Amn\n");
      print_mat(Amn,max_rows,norbs*norbs, outfile);
    }
    timer_off("(A|mn)");

    //Transform to Ami
    // (A|mi) = (Amn)C_ni
    timer_on("(A|mi)");
    C_DGEMM('N', 'N', p_sizes[block]*norbs, nact_docc, norbs, 1.0, &(Amn[0][0]),
        norbs, &(C_docc_[0][0]), nact_docc, 0.0, &(Ami[0][0]), nact_docc);
    timer_off("(A|mi)");

    if (debug_) {
      fprintf(outfile, "  Ami\n");
      print_mat(Ami,max_rows,nact_docc*norbs, outfile);
    }

    timer_on("(A|ia)");
    #ifdef _MKL
       int mkl_nthreads = mkl_get_max_threads();
       mkl_set_num_threads(1);
    #endif

    //Transform to Aia
    // (A|ia) = C_ma(A|mi)
    #pragma omp parallel for
    for (int A = 0; A<p_sizes[block]; A++) {
      C_DGEMM('T', 'N', nact_docc, nact_virt, norbs, 1.0, &(Ami[A][0]),
        nact_docc, &(C_virt_[0][0]), nact_virt, 0.0, &(Aia[p_starts[block]+A][0]), nact_virt_);
    }

    #ifdef _MKL
       mkl_set_num_threads(mkl_nthreads);
    #endif
    timer_off("(A|ia)");

  }
  if (debug_) {
    fprintf(outfile, "  Aia\n");
    print_mat(Aia,naux_raw_,nact_docc*nact_virt, outfile);
  }

  delete[] buffer;
  delete[] eri;
  free(p_starts);
  free(p_sizes);
  free(block_starts);
  free(block_stops);
  free(block_sizes);
  free_block(Amn);
  free_block(Ami);

  return Aia;
}
void DFMP2::form_Qia_core_parallel()
{

  //Convention
  int norbs =  nso_;
  int nact_docc =  nact_docc_;
  int nact_virt =  nact_virt_;
  //printf("nact_docc_ is %d in form_Qia_core()\n", nact_docc_);
  Qia_ = form_Aia_core_parallel();

  //Setup fitting tensor
  if (options_.get_str("RI_FITTING_TYPE") == "RAW")
      form_Wm12_raw();
  else if (options_.get_str("RI_FITTING_TYPE") == "FINISHED")
      form_Wm12_fin();
  else if (options_.get_str("RI_FITTING_TYPE") == "CHOLESKY")
      form_Wp12_chol();

  unsigned long int three_mem = naux_raw_*nact_docc_*(ULI)nact_virt_;
  unsigned long int available_mem = df_memory_-three_mem-naux_raw_*(ULI)naux_fin_;
  unsigned long int max_cols = available_mem/(ULI)naux_raw_;
  if (max_cols > nact_docc_*(ULI)nact_virt_)
    max_cols = nact_docc_*(ULI)nact_virt;
  if (max_cols < 1)
    max_cols = 1;

  unsigned long int nblocks = nact_docc_*(ULI)nact_virt_/max_cols;
  if (nblocks*max_cols < nact_docc_*nact_virt_)
    nblocks++;

  int* column_starts = init_int_array(nblocks);
  int* column_sizes = init_int_array(nblocks);

  for (int block = 0; block < nblocks; block++) {
    if (block == nblocks - 1) {
      column_sizes[block] = nact_docc_*(ULI)nact_virt_ - column_starts[block];
      break;
    }
    column_sizes[block] = max_cols;
    column_starts[block+1] = column_starts[block] + column_sizes[block];
  }

  //Apply fitting by blocks
  double** Abuffer;
  if (options_.get_str("RI_FITTING_TYPE") == "CHOLESKY") {
      Abuffer = block_matrix(max_cols, naux_raw_);
  } else {
      Abuffer = block_matrix(naux_raw_, max_cols);
  }
  timer_on("(Q|ia)");
  for (int block = 0; block < nblocks; block++) {

     if (options_.get_str("RI_FITTING_TYPE") == "CHOLESKY") {

      #pragma omp parallel for schedule(static)
      for (int A = 0; A< naux_raw_; A++) {
        C_DCOPY(column_sizes[block], &Qia_[A][column_starts[block]], 1, &Abuffer[0][A], naux_raw_);
      }

      int info = C_DPOTRS('U',naux_raw_,column_sizes[block],&W_[0][0],naux_raw_,&Abuffer[0][0],naux_raw_);

      #pragma omp parallel for schedule(static)
      for (int A = 0; A< naux_raw_; A++) {
        C_DCOPY(column_sizes[block], &Abuffer[0][A], naux_raw_, &Qia_[A][column_starts[block]], 1);
      }
    } else {
      #pragma omp parallel for schedule(static)
      for (int A = 0; A< naux_raw_; A++) {
        C_DCOPY(column_sizes[block], &Qia_[A][column_starts[block]], 1, &Abuffer[A][0], 1);
      }

      //Multiply back into Qia_ with tight striping
      C_DGEMM('T','N',naux_fin_, column_sizes[block], naux_raw_, 1.0, W_[0], naux_fin_, Abuffer[0], max_cols, 0.0, &Qia_[0][column_starts[block]], nact_docc_*(ULI) nact_virt_);
    }
  }
  timer_off("(Q|ia)");

  if (debug_) {
    fprintf(outfile,"  Qia");
    print_mat(Qia_,naux_fin_,nact_docc_*nact_virt_,outfile);
  }

  free_block(W_);
  free_block(Abuffer);
  free(column_sizes);
  free(column_starts);
}
void DFMP2::evaluate_contributions_core_parallel()
{
   int nthreads = 1;
  int rank = 0;

    // Ben Seibert, your stuff goes here
#ifdef HAVE_MPI
  int nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  #ifdef _OPENMP
    nthreads = 1;
  #endif

  #ifdef _MKL
    int mkl_nthreads = mkl_get_max_threads();
    mkl_set_num_threads(mkl_nthreads);
  #endif

  double*** I = new double**[nprocs];
  for (int r = 0 ; r<nprocs; r++) {
    I[r] = block_matrix(nact_virt_,nact_virt_);
  }

  double iajb, ibja, denom, perm_scale;
  long double e_ss, e_os;
  int a, b, i, j;

  e_ss = 0.0;
  e_os = 0.0;

  for (i=rank; i < nact_docc_; i+=nprocs) {
    for(j=0;j<=i;j++){
      if (rank == 0) timer_on("(ia|jb)");
      C_DGEMM('T','N',nact_virt_,nact_virt_,naux_fin_,1.0,&(Qia_[0][i*nact_virt_]),
        nact_docc_*(ULI) nact_virt_,&(Qia_[0][j*nact_virt_]),nact_docc_*(ULI) nact_virt_,0.0,&(I[rank][0][0]),nact_virt_);
      if (rank == 0) timer_off("(ia|jb)");
      if (rank == 0) timer_on("T_MP2");
      for (a=0; a < nact_virt_; ++a) {
        for (b=0; b < nact_virt_; ++b) {
          iajb = I[rank][a][b];
          ibja = I[rank][b][a];
          denom = 1.0 /
            (eps_docc_[i] + eps_docc_[j] -
             eps_virt_[a] - eps_virt_[b]);

          perm_scale = ((i == j)? 1.0 : 2.0 );

          e_ss += perm_scale*(iajb*iajb-iajb*ibja)*denom;
          e_os += perm_scale*(iajb*iajb)*denom;
        }
      }
      if (rank == 0) timer_off("T_MP2");
    }
  }

   double *in_buf_ess, *out_buf_ess, *in_buf_eos, *out_buf_eos;
   in_buf_ess = (double *) malloc(1*sizeof(double));
   in_buf_eos = (double *) malloc(1*sizeof(double));
   out_buf_ess = (double *) malloc(1*sizeof(double));
   out_buf_eos = (double *) malloc(1*sizeof(double));

   *in_buf_ess = e_ss;
   *in_buf_eos = e_os;


   MPI_Reduce(in_buf_ess, out_buf_ess, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(in_buf_eos, out_buf_eos, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

   e_ss = *out_buf_ess;
   e_os = *out_buf_eos;

  for (int r = 0; r<nthreads; r++)
  {
    printf("Freeing block %d of I\n", r);
    free_block(I[r]);
  }
  delete[] I;

  #ifdef _MKL
    mkl_set_num_threads(mkl_nthreads);
  #endif

  E_ss_ += e_ss;
  E_os_ += e_os;
#endif
}
void DFMP2::free_Qia_core()
{
    free_block(Qia_);
}
double DFMP2::compute_E_old()
{
        /**g
  fprintf(outfile, "\t\t*************************\n");
  fprintf(outfile, "\t\t*                       *\n");
  fprintf(outfile, "\t\t*         DF-MP2        *\n");
  fprintf(outfile, "\t\t*                       *\n");
  fprintf(outfile, "\t\t*                       *\n");
  fprintf(outfile, "\t\t*************************\n");
  fflush(outfile);

  //How many processes are we running on?
  int nproc = Communicator::world->nproc();
  //What process number is this?
  int rank = Communicator::world->me();

  fprintf(outfile,"\n\n  Running on %d processes.\n",nproc);
  fflush(outfile);
  // Required for libmints, allocates and computes:

  int nirreps;
  int *clsdpi;
  int *orbspi;
  int *frzcpi;
  int *frzvpi;
  nirreps = chkpt_->rd_nirreps();
  clsdpi = chkpt_->rd_clsdpi();
  orbspi = chkpt_->rd_orbspi();
  frzcpi = chkpt_->rd_frzcpi();
  frzvpi = chkpt_->rd_frzvpi();
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


  double escf;
  escf = chkpt_->rd_escf();

  // Create a new matrix factory
  boost::shared_ptr<MatrixFactory> factory = boost::shared_ptr<MatrixFactory>(new MatrixFactory());

  // Initialize the factory with data from checkpoint
  factory->init_with_chkpt(chkpt_);
  //factory->init_with(nirreps,orbspi,orbspi);

  // Read in C coefficients
  double **vectors;
  vectors = chkpt_->rd_scf();

  //Do some transform
  SimpleMatrix *C_so =
    factory->create_simple_matrix("MO coefficients (SO basisset_)");
  if (vectors == NULL) {
    fprintf(stderr, "Could not find MO coefficients. Run cscf first.\n");
    return EXIT_FAILURE;
  }
  C_so->set(vectors);
  free_block(vectors);
  double **Umat;
  Umat = chkpt_->rd_usotbf();
  SimpleMatrix *U = factory->create_simple_matrix("SO->BF");
  U->set(Umat);
  free_block(Umat);
  // Transfrom the eigenvectors from the SO to the AO basisset_
  SimpleMatrix *C =
    factory->create_simple_matrix("MO coefficients (AO basisset_)");
  C->gemm(true, false, 1.0, U, C_so, 0.0);

  // Load in orbital energies
  double *orbital_energies;
  orbital_energies = chkpt_->rd_evals();

  // Form ribasis object:
  boost::shared_ptr<BasisSet> ribasis =boost::shared_ptr<BasisSet>(new BasisSet(chkpt_, "DF_BASIS_MP2"));
  int naux = ribasis->nbf();
  boost::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();

  //Put the orbitals and C matrix into some nice arrays
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
    //fprintf(outfile, "  C_occ:\n");
    //print_mat(Co,norbs,nact_docc,outfile);
    //fprintf(outfile, "  C_virt:\n");
    //print_mat(Cv,norbs,nact_virt,outfile);
    //fprintf(outfile, "  Eps_occ:\n");
    //for (int i = 0; i<nact_docc; i++)
    //    fprintf(outfile,"  %5d: %14.10f\n", i+1, epsilon_act_docc[i]);
    //fprintf(outfile, "  Eps_virt:\n");
    //for (int i = 0; i<nact_virt; i++)
    //    fprintf(outfile,"  %5d: %14.10f\n", i+1+nact_docc, epsilon_act_virt[i]);

  fprintf(outfile, "\n\t\t==============================================\n");
  fprintf(outfile, "\t\t #ORBITALS #RI  FOCC DOCC AOCC AVIR VIRT FVIR \n");
  fprintf(outfile, "\t\t----------------------------------------------\n");
  fprintf(outfile, "\t\t  %5d  %5d  %4d %4d %4d %4d %4d %4d\n",
          norbs,naux,nfocc,ndocc,nact_docc,nact_virt,nvirt,nfvir);
  fprintf(outfile, "\t\t==============================================\n");
  fflush(outfile);

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
  fflush(outfile);


  // Create integral factory
  IntegralFactory rifactory(ribasis, zero, basisset_, basisset_);
  IntegralFactory rifactory_J(ribasis, zero, ribasis, zero);

  TwoBodyAOInt* eri = rifactory.eri();
  TwoBodyAOInt* Jint = rifactory_J.eri();
  double **J = block_matrix(naux, naux);
  double **J_mhalf = block_matrix(naux, naux);
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
  // print_mat(J, naux, naux, outfile);
  // Invert J matrix
  // invert_matrix(J, J_inverse, naux, outfile);
  // fprintf(outfile, "J^-1:\n");
  // print_mat(J_inverse, naux, naux, outfile);
  // Form J^-1/2
  // First, diagonalize J
  // the C_DSYEV call replaces the original matrix J with its eigenvectors
  double* eigval = init_array(naux);
  int lwork = naux * 3;
  double* work = init_array(lwork);
  int stat = C_DSYEV('v','u',naux,J[0],naux,eigval,
    work,lwork);
  if (stat != 0) {
    fprintf(outfile, "C_DSYEV failed\n");
    exit(PSI_RETURN_FAILURE);
  }
  free(work);

  // Now J contains the eigenvectors of the original J
  // Copy J to J_copy
  double **J_copy = block_matrix(naux, naux);
  C_DCOPY(naux*naux,J[0],1,J_copy[0],1);

  // Now form J^{-1/2} = U(T)*j^{-1/2}*U,
  // where j^{-1/2} is the diagonal matrix of the inverse square roots
  // of the eigenvalues, and U is the matrix of eigenvectors of J
  for (int i=0; i<naux; i++) {
    if (eigval[i] < 1.0E-10)
      eigval[i] = 0.0;
    else {
      eigval[i] = 1.0 / sqrt(eigval[i]);
    }
    // scale one set of eigenvectors by the diagonal elements j^{-1/2}
    C_DSCAL(naux, eigval[i], J[i], 1);
  }
  free(eigval);

  // J_mhalf = J_copy(T) * J
  C_DGEMM('t','n',naux,naux,naux,1.0,
    J_copy[0],naux,J[0],naux,0.0,J_mhalf[0],naux);

  free_block(J);
  free_block(J_copy);

  timer_off("Form J");

  fprintf(outfile," J is formed\n");
  fflush(outfile);

  //fprintf(outfile,"  W:\n");
  //print_mat(J_mhalf,naux,naux,outfile);

  const double *buffer = eri->buffer();


  //How many doubles do we have?
  double doubles = memory_/(1.0*sizeof(double))*0.9;
  double doubles_available = doubles-(1.0)*naux*naux;

  //What is the size of the virtual block?
  int vir_blk_size = (int) floor(doubles_available/(4.0*naux*nact_docc));
  if (vir_blk_size > nact_virt)
      vir_blk_size = nact_virt;
  if (vir_blk_size < 1)
      vir_blk_size = 1;
  int n_full_blks = nact_virt/vir_blk_size;
  int gimp_blk_size = nact_virt - n_full_blks*vir_blk_size;

  //How many actual blocks are there in one dimension?
  int n_virt_blks = n_full_blks;
  //Is there a gimp block?
  if (gimp_blk_size != 0)
    n_virt_blks++;

  if (nproc > n_virt_blks*(n_virt_blks+1)/2) {
    n_full_blks = (int) floor((-1.0 + sqrt((1.0+8.0*((double)nproc)))) / 2.0);
    vir_blk_size = nact_virt/n_full_blks;
    gimp_blk_size = nact_virt - n_full_blks*vir_blk_size;
    n_virt_blks = n_full_blks;
  if (gimp_blk_size != 0) {
    n_virt_blks++;
    }
  }
  fprintf(outfile,"  %d virtual/virtual blocks of dimension %d will be used.\n",n_full_blks,vir_blk_size);
  fprintf(outfile,"  The gimp index dimension is %d\n",gimp_blk_size);
  fprintf(outfile,"\n");
  fflush(outfile);


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
  double** Pia = block_matrix(naux,nact_docc*vir_blk_size);
  double** Pjb = block_matrix(naux,nact_docc*vir_blk_size);

  double **Qia = block_matrix(nact_docc * vir_blk_size, naux);
  double **Qjb = block_matrix(nact_docc * vir_blk_size, naux);

  double **I = block_matrix(nact_docc,nact_docc);

  //The energy!
  register double emp2 = 0.0;

//Outer loop
int block_A, block_B, hash;
//Actual block start and sizes
int a_start, b_start, a_size, b_size, gimped_a=0;
//Scale due to permutational symmetry
double perm_scale;
//AO integral indexing
int numPshell,Pshell,MU,NU,P,mu,nu,nummu,numnu,omu,onu;
for (block_A = 0; block_A < n_virt_blks; block_A++) {
  for (block_B = 0; block_B <= block_A; block_B++) {
    double emp2_ab = 0.0;

  //hash value
  hash = block_A*(block_A+1)/2+block_B;
  //Should we do this block?
  if (hash%nproc != rank) {
    //GTFO
    continue;
  }

  //so where is the job?
  a_start = vir_blk_size*block_A;
  b_start = vir_blk_size*block_B;
  a_size = vir_blk_size;
  b_size = vir_blk_size;
  if (block_A == n_full_blks) {
    a_size = gimp_blk_size;
    if (gimped_a == 0) {
      free_block(Pia);
      Pia = block_matrix(naux,nact_docc*a_size);
      gimped_a = 1;
    }
//    free_block(Qia);
//    Qia = block_matrix(nact_docc * a_size, naux);
  }
  if (block_B == n_full_blks) {
    b_size = gimp_blk_size;
    free_block(Pjb);
    Pjb = block_matrix(naux,nact_docc*b_size);
//    free_block(Qjb);
//    Qjb = block_matrix(nact_docc * b_size, naux);
  }

  //How much is this block worth due to permutational symmetry?
  perm_scale = (block_A == block_B)? 1.0 : 2.0;

  //fprintf(outfile,"  Thread %d: a_start = %d, a_size = %d, b_start = %d, b_size = %d\n", rank, a_start, a_size, b_start, b_size);

  //Inner loop
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
    if (debug_) {
      fprintf(outfile, "  Amn for shell %d\n",Pshell);
      print_mat(ao_buffer,numPshell,norbs*norbs,outfile);
    }
    // Do transform
    C_DGEMM('N', 'N', numPshell*norbs, nact_docc, norbs, 1.0, &(ao_buffer[0][0]),
        norbs, &(Co[0][0]), nact_docc, 0.0, &(half_buffer[0][0]), nact_docc);

    if (debug_) {
      fprintf(outfile, "  Ami for shell %d\n",Pshell);
      print_mat(half_buffer,numPshell,norbs*ndocc,outfile);
    }

    for (int P=0, index=0; P < numPshell; ++P) {
      int oP = ribasis->shell(Pshell)->function_index() + P;

      C_DGEMM('T', 'N', a_size, nact_docc, norbs, 1.0, &(Cv[0][a_start]),
        nact_virt, &(half_buffer[P][0]), nact_docc, 0.0, &(Pia[oP][0]), nact_docc);
      C_DGEMM('T', 'N', b_size, nact_docc, norbs, 1.0, &(Cv[0][b_start]),
        nact_virt, &(half_buffer[P][0]), nact_docc, 0.0, &(Pjb[oP][0]), nact_docc);

    }
  } // end loop over P shells; done with forming MO basisset_ (P|ia)'s

  if (debug_) {
    fprintf(outfile,"  Aai:\n");
    print_mat(Pia,naux,nact_docc*nact_virt,outfile);
  }

  //Embed the fitting

  C_DGEMM('T','N',nact_docc*a_size,naux,naux,
    1.0, Pia[0], nact_docc*a_size, J_mhalf[0], naux,
    0.0, Qia[0], naux);

  C_DGEMM('T','N',nact_docc*b_size,naux,naux,
    1.0, Pjb[0], nact_docc*b_size, J_mhalf[0], naux,
    0.0, Qjb[0], naux);

  if (debug_) {
    fprintf(outfile,"  Qai:\n");
    print_mat(Qia,nact_docc*nact_virt,naux,outfile);
  }

  for (int a=0; a < a_size; ++a) {
    for (int b=0; b < b_size; ++b) {

      C_DGEMM('N','T',nact_docc,nact_docc,naux,1.0,&(Qia[a*nact_docc][0]),
        naux,&(Qjb[b*nact_docc][0]),naux,0.0,&(I[0][0]),nact_docc);

      double iajb, ibja, tval, denom;
      for (int i=0; i < nact_docc; ++i) {
        for (int j=0; j < nact_docc; ++j) {
          iajb = I[i][j];
          ibja = I[j][i];
          denom = 1.0 /
            (epsilon_act_docc[i] + epsilon_act_docc[j] -
             epsilon_act_virt[a+a_start] - epsilon_act_virt[b+b_start]);

          tval = (2.0*iajb*iajb - iajb*ibja)*denom;
          emp2 += perm_scale*tval;
          emp2_ab += perm_scale*tval;
        }
      }
    } // end loop over j<=i
  } // end loop over i
    //fprintf(outfile,"  Thread %d: Block A = %d, Block B = %d, Hash = %d, energy contribution = %14.10f\n",rank, block_A, block_B, hash, emp2_ab);
    //fflush(outfile);
    //fprintf(outfile,"  Thread %d: Block A = %d, Block B = %d, Hash = %d, energy contribution = %14.10f\n",rank, block_A, block_B, hash, emp2_ab);
  //end inner loop
}} // end outer loop (parallelization)

  free_block(Qjb);
  free_block(Qia);
  free_block(Pia);
  free_block(Pjb);
  free_block(I);
  free_block(ao_buffer);
  free_block(half_buffer);
  free_block(J_mhalf);
  //Reduce the energy
  Communicator::world->sum(emp2);

  //Print the results
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

  // Get out
  fflush(outfile);
  return emp2;
**/
return 0.0;

}

}}

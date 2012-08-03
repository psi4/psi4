#include "sapt2.h"
#include <physconst.h>

namespace psi { namespace sapt {

SAPT2::SAPT2(Options& options, boost::shared_ptr<PSIO> psio,
  boost::shared_ptr<Chkpt> chkpt) : SAPT(options, psio, chkpt),
  e_elst10_(0.0),
  e_elst12_(0.0),
  e_exch10_(0.0),
  e_exch10_s2_(0.0),
  e_exch11_(0.0),
  e_exch12_(0.0),
  e_ind20_(0.0),
  e_ind22_(0.0),
  e_exch_ind20_(0.0),
  e_exch_ind22_(0.0),
  e_disp20_(0.0),
  e_no_disp20_(0.0),
  e_exch_disp20_(0.0),
  e_sapt0_(0.0),
  e_sapt2_(0.0)
{
  psio_->open(PSIF_SAPT_AA_DF_INTS,PSIO_OPEN_NEW);
  psio_->open(PSIF_SAPT_BB_DF_INTS,PSIO_OPEN_NEW);
  psio_->open(PSIF_SAPT_AB_DF_INTS,PSIO_OPEN_NEW);
  psio_->open(PSIF_SAPT_AMPS,PSIO_OPEN_NEW);

  maxiter_ = options_.get_int("MAXITER");
  e_conv_ = options_.get_double("E_CONVERGENCE");
  d_conv_ = options_.get_double("D_CONVERGENCE");

  nat_orbs_ = options.get_bool("NAT_ORBS");
  nat_orbs_t2_ = options.get_bool("NAT_ORBS_T2");
  occ_cutoff_ = options.get_double("OCC_TOLERANCE");

  ioff_ = (int *) malloc (sizeof(int) * (nso_*(nso_+1)/2));
  index2i_ = (int *) malloc (sizeof(int) * (nso_*(nso_+1)/2));
  index2j_ = (int *) malloc (sizeof(int) * (nso_*(nso_+1)/2));

  ioff_[0] = 0;

  for (int i=1; i < (nso_*(nso_+1)/2); i++) {
    ioff_[i] = ioff_[i-1] + i;
  }

  for (int i=0, ij=0; i<nso_; i++) {
    for (int j=0; j<=i; j++, ij++) {
      index2i_[ij] = i;
      index2j_[ij] = j;
  }}

  wBAR_ = NULL;
  wABS_ = NULL;
  wBAA_ = NULL;
  wBRR_ = NULL;
  wABB_ = NULL;
  wASS_ = NULL;
  no_evalsA_ = NULL;
  no_evalsB_ = NULL;
  no_CA_ = NULL;
  no_CB_ = NULL;
}

SAPT2::~SAPT2()
{
  if (wBAR_ != NULL) free_block(wBAR_);
  if (wABS_ != NULL) free_block(wABS_);

  if (wBAA_ != NULL) free_block(wBAA_);
  if (wBRR_ != NULL) free_block(wBRR_);
  if (wABB_ != NULL) free_block(wABB_);
  if (wASS_ != NULL) free_block(wASS_);

  if (nat_orbs_) {
    if (no_evalsA_ != NULL) free(no_evalsA_);
    if (no_evalsB_ != NULL) free(no_evalsB_);
    if (no_CA_ != NULL) free_block(no_CA_);
    if (no_CB_ != NULL) free_block(no_CB_);
  }

  free(ioff_);
  free(index2i_);
  free(index2j_);

  psio_->close(PSIF_SAPT_AA_DF_INTS,1);
  psio_->close(PSIF_SAPT_BB_DF_INTS,1);
  psio_->close(PSIF_SAPT_AB_DF_INTS,1);
  psio_->close(PSIF_SAPT_AMPS,1);
}

double SAPT2::compute_energy()
{
  print_header();

  timer_on("DF Integrals       ");
    df_integrals();
  timer_off("DF Integrals       ");
  timer_on("Omega Integrals    ");
    w_integrals();
  timer_off("Omega Integrals    ");
  timer_on("Amplitudes         ");
    amplitudes();
  timer_off("Amplitudes         ");
  timer_on("Elst10             ");
    elst10();
  timer_off("Elst10             ");
  timer_on("Exch10 S^2         ");
    exch10_s2();
  timer_off("Exch10 S^2         ");
  timer_on("Exch10             ");
    exch10();
  timer_off("Exch10             ");
  timer_on("Ind20,r            ");
    ind20r();
  timer_off("Ind20,r            ");
  timer_on("Exch-Ind20,r       ");
    exch_ind20r();
  timer_off("Exch-Ind20,r       ");
  timer_on("Disp20             ");
    disp20();
  timer_off("Disp20             ");
  timer_on("Exch-Disp20        ");
    exch_disp20();
  timer_off("Exch-Disp20        ");
  timer_on("Elst12             ");
    elst12();
  timer_off("Elst12             ");
  timer_on("Exch11             ");
    exch11();
  timer_off("Exch11             ");
  timer_on("Exch12             ");
    exch12();
  timer_off("Exch12             ");
  timer_on("Ind22              ");
    ind22();
  timer_off("Ind22              ");

  print_results();

  return (e_sapt0_);
}

void SAPT2::print_header()
{
  fprintf(outfile,"        SAPT2  \n");
  fprintf(outfile,"    Ed Hohenstein\n") ;
  fprintf(outfile,"     6 June 2009\n") ;
  fprintf(outfile,"\n");
  fprintf(outfile,"      Orbital Information\n");
  fprintf(outfile,"  --------------------------\n");
  fprintf(outfile,"    NSO        = %9d\n",nso_);
  if (nsoA_ != nsoB_) {
    fprintf(outfile,"    NSO A      = %9d\n",nsoA_);
    fprintf(outfile,"    NSO B      = %9d\n",nsoB_);
  }
  fprintf(outfile,"    NMO        = %9d\n",nmo_);
  if (nmoA_ != nmoB_) {
    fprintf(outfile,"    NMO A      = %9d\n",nmoA_);
    fprintf(outfile,"    NMO B      = %9d\n",nmoB_);
  }
  fprintf(outfile,"    NRI        = %9d\n",ndf_);
  fprintf(outfile,"    NOCC A     = %9d\n",noccA_);
  fprintf(outfile,"    NOCC B     = %9d\n",noccB_);
  fprintf(outfile,"    FOCC A     = %9d\n",foccA_);
  fprintf(outfile,"    FOCC B     = %9d\n",foccB_);
  fprintf(outfile,"    NVIR A     = %9d\n",nvirA_);
  fprintf(outfile,"    NVIR B     = %9d\n",nvirB_);
  fprintf(outfile,"\n");

  long int mem = (long int) memory_;
  mem /= 8L;
  long int occ = noccA_;
  if (noccB_ > noccA_)
    occ = noccB_;
  long int vir = nvirA_;
  if (nvirB_ > nvirA_)
    vir = nvirB_;
  long int ovov = occ*occ*vir*vir;
  long int vvnri = vir*vir*ndf_;
  double memory = 8.0*(vvnri + ovov*3L)/1000000.0;
  if (print_) {
    fprintf(outfile,"    Estimated memory usage: %.1lf MB\n\n",memory);
    fflush(outfile);
  }
  if (options_.get_bool("SAPT_MEM_CHECK"))
    if (mem < vvnri + ovov*3L)
      throw PsiException("Not enough memory", __FILE__,__LINE__);

  fflush(outfile);
}

void SAPT2::print_results()
{
  e_sapt0_ = eHF_ + e_disp20_ + e_exch_disp20_;
  e_sapt2_ = e_sapt0_ + e_elst12_ + e_exch11_ + e_exch12_  + e_ind22_
    + e_exch_ind22_;
  double dHF = eHF_ - (e_elst10_ + e_exch10_ + e_ind20_ + e_exch_ind20_);

  double tot_elst = e_elst10_ + e_elst12_;
  double tot_exch = e_exch10_ + e_exch11_ + e_exch12_;
  double tot_ind = e_ind20_ + e_exch_ind20_ + dHF + e_ind22_ + e_exch_ind22_;
  double tot_ct = e_ind20_ + e_exch_ind20_ + e_ind22_ + e_exch_ind22_;
  double tot_disp = e_disp20_ + e_exch_disp20_;

  fprintf(outfile,"\n    SAPT Results  \n");
  fprintf(outfile,"  -----------------------------------------------------------------------\n");
  fprintf(outfile,"    Electrostatics     %16.8lf mH %16.8lf kcal mol^-1\n",
    tot_elst*1000.0,tot_elst*_hartree2kcalmol);
  fprintf(outfile,"      Elst10,r         %16.8lf mH %16.8lf kcal mol^-1\n",
    e_elst10_*1000.0,e_elst10_*_hartree2kcalmol);
  fprintf(outfile,"      Elst12,r         %16.8lf mH %16.8lf kcal mol^-1\n\n",
    e_elst12_*1000.0,e_elst12_*_hartree2kcalmol);
  fprintf(outfile,"    Exchange           %16.8lf mH %16.8lf kcal mol^-1\n",
    tot_exch*1000.0,tot_exch*_hartree2kcalmol);
  fprintf(outfile,"      Exch10           %16.8lf mH %16.8lf kcal mol^-1\n",
    e_exch10_*1000.0,e_exch10_*_hartree2kcalmol);
  fprintf(outfile,"      Exch10(S^2)      %16.8lf mH %16.8lf kcal mol^-1\n",
    e_exch10_s2_*1000.0,e_exch10_s2_*_hartree2kcalmol);
  fprintf(outfile,"      Exch11(S^2)      %16.8lf mH %16.8lf kcal mol^-1\n",
    e_exch11_*1000.0,e_exch11_*_hartree2kcalmol);
  fprintf(outfile,"      Exch12(S^2)      %16.8lf mH %16.8lf kcal mol^-1\n\n",
    e_exch12_*1000.0,e_exch12_*_hartree2kcalmol);
  fprintf(outfile,"    Induction          %16.8lf mH %16.8lf kcal mol^-1\n",
    tot_ind*1000.0,tot_ind*_hartree2kcalmol);
  fprintf(outfile,"      Ind20,r          %16.8lf mH %16.8lf kcal mol^-1\n",
    e_ind20_*1000.0,e_ind20_*_hartree2kcalmol);
  fprintf(outfile,"      Ind22            %16.8lf mH %16.8lf kcal mol^-1\n",
    e_ind22_*1000.0,e_ind22_*_hartree2kcalmol);
  fprintf(outfile,"      Exch-Ind20,r     %16.8lf mH %16.8lf kcal mol^-1\n",
    e_exch_ind20_*1000.0,e_exch_ind20_*_hartree2kcalmol);
  fprintf(outfile,"      Exch-Ind22       %16.8lf mH %16.8lf kcal mol^-1\n",
    e_exch_ind22_*1000.0,e_exch_ind22_*_hartree2kcalmol);
  fprintf(outfile,"      delta HF,r       %16.8lf mH %16.8lf kcal mol^-1\n\n",
    dHF*1000.0,dHF*_hartree2kcalmol);
  fprintf(outfile,"    Dispersion         %16.8lf mH %16.8lf kcal mol^-1\n",
    tot_disp*1000.0,tot_disp*_hartree2kcalmol);
  fprintf(outfile,"      Disp20           %16.8lf mH %16.8lf kcal mol^-1\n",
    e_disp20_*1000.0,e_disp20_*_hartree2kcalmol);
  fprintf(outfile,"      Exch-Disp20      %16.8lf mH %16.8lf kcal mol^-1\n\n",
    e_exch_disp20_*1000.0,e_exch_disp20_*_hartree2kcalmol);

  fprintf(outfile,"    Total HF           %16.8lf mH %16.8lf kcal mol^-1\n",
    eHF_*1000.0,eHF_*_hartree2kcalmol);
  fprintf(outfile,"    Total SAPT0        %16.8lf mH %16.8lf kcal mol^-1\n",
    e_sapt0_*1000.0,e_sapt0_*_hartree2kcalmol);
  fprintf(outfile,"    Total SAPT2        %16.8lf mH %16.8lf kcal mol^-1\n",
    e_sapt2_*1000.0,e_sapt2_*_hartree2kcalmol);

  Process::environment.globals["SAPT ELST ENERGY"] = tot_elst;
  Process::environment.globals["SAPT EXCH ENERGY"] = tot_exch;
  Process::environment.globals["SAPT IND ENERGY"] = tot_ind;
  Process::environment.globals["SAPT CT ENERGY"] = tot_ct;
  Process::environment.globals["SAPT DISP ENERGY"] = tot_disp;
  Process::environment.globals["SAPT SAPT0 ENERGY"] = e_sapt0_;
  Process::environment.globals["SAPT SAPT2 ENERGY"] = e_sapt2_;
  Process::environment.globals["SAPT ENERGY"] = e_sapt2_;
  Process::environment.globals["CURRENT ENERGY"] = Process::environment.globals["SAPT ENERGY"];
}

void SAPT2::df_integrals()
{
  // Get fitting metric
  boost::shared_ptr<FittingMetric> metric = boost::shared_ptr<FittingMetric>(
    new FittingMetric(ribasis_));
  metric->form_eig_inverse();
  double **J_temp = metric->get_metric()->pointer();
  double **J_mhalf = block_matrix(ndf_,ndf_);
  C_DCOPY(ndf_*ndf_,J_temp[0],1,J_mhalf[0],1);
  metric.reset();

  boost::shared_ptr<IntegralFactory> rifactory_J(new IntegralFactory(ribasis_,
    zero_, ribasis_, zero_));
  boost::shared_ptr<TwoBodyAOInt> Jint = boost::shared_ptr<TwoBodyAOInt>(
    rifactory_J->eri());
  const double *Jbuffer = Jint->buffer();

  boost::shared_ptr<IntegralFactory> ao_eri_factory(new IntegralFactory(
    basisset_, basisset_, basisset_, basisset_));
  boost::shared_ptr<TwoBodyAOInt> ao_eri = boost::shared_ptr<TwoBodyAOInt>(
    ao_eri_factory->eri());
  const double *ao_buffer = ao_eri->buffer();

  double *Schwartz = init_array(basisset_->nshell()*(basisset_->nshell()+1)/2);
  double *DFSchwartz = init_array(ribasis_->nshell());

  for(int P=0,PQ=0;P<basisset_->nshell();P++) {
    int numw = basisset_->shell(P).nfunction();
    for(int Q=0;Q<=P;Q++,PQ++) {
      int numx = basisset_->shell(Q).nfunction();
      double tei, max=0.0;

      ao_eri->compute_shell(P, Q, P, Q);

      for(int w=0;w<numw;w++) {
        for(int x=0;x<numx;x++) {
          int index = ( (w*numx + x) * numw + w) * numx + x;
          tei = ao_buffer[index];
          if(fabs(tei) > max) max = fabs(tei);
        }
      }
      Schwartz[PQ] = max;
    }
  }

  for(int P=0;P<ribasis_->nshell();P++) {
    int numw = ribasis_->shell(P).nfunction();
    double tei, max=0.0;

    Jint->compute_shell(P, 0, P, 0);

    for(int w=0;w<numw;w++) {
      tei = Jbuffer[w];
      if(fabs(tei) > max) max = fabs(tei);
    }
    DFSchwartz[P] = max;
  }

  boost::shared_ptr<IntegralFactory> rifactory(new IntegralFactory(ribasis_,
    zero_, basisset_, basisset_));

  int nthreads = 1;
  #ifdef _OPENMP
    nthreads = omp_get_max_threads();
  #endif
  int rank = 0;

  boost::shared_ptr<TwoBodyAOInt> *eri =
    new boost::shared_ptr<TwoBodyAOInt>[nthreads];
  const double **buffer = new const double*[nthreads];
  for(int i = 0;i < nthreads;++i){
    eri[i] = boost::shared_ptr<TwoBodyAOInt>(rifactory->eri());
    buffer[i] = eri[i]->buffer();
  }

  int mn;
  int maxPshell = 0;
  for (int Pshell=0; Pshell < ribasis_->nshell(); ++Pshell) {
    int numPshell = ribasis_->shell(Pshell).nfunction();
    if (numPshell > maxPshell) maxPshell = numPshell;
  }

  psio_->open(PSIF_SAPT_TEMP,0);

  double** AO_RI = block_matrix(maxPshell,nso_*nso_);
  double* halftrans = init_array(nmoA_*nso_);
  double** MO_RI = block_matrix(maxPshell,nmoA_*nmoA_);

  psio_address next_DF_MO = PSIO_ZERO;
  psio_address next_bare_AR = PSIO_ZERO;

  int nshelltri = basisset_->nshell()*(basisset_->nshell()+1)/2;

  for (int Pshell=0; Pshell < ribasis_->nshell(); ++Pshell) {
    int numPshell = ribasis_->shell(Pshell).nfunction();

    #pragma omp for private(mn) schedule(guided)
    for (mn=0; mn < nshelltri; ++mn) {
      #ifdef _OPENMP
        rank = omp_get_thread_num();
      #endif
      int MU = index2i_[mn];
      int NU = index2j_[mn];
      int nummu = basisset_->shell(MU).nfunction();
      int numnu = basisset_->shell(NU).nfunction();

      if (sqrt(Schwartz[mn]*DFSchwartz[Pshell])>schwarz_) {
        eri[rank]->compute_shell(Pshell, 0, MU, NU);

        for (int P=0, index=0; P < numPshell; ++P) {

          for (int mu=0; mu < nummu; ++mu) {
            int omu = basisset_->shell(MU).function_index() + mu;

            for (int nu=0; nu < numnu; ++nu, ++index) {
              int onu = basisset_->shell(NU).function_index() + nu;

              AO_RI[P][omu*nso_+onu] = buffer[rank][index];
              AO_RI[P][onu*nso_+omu] = buffer[rank][index];
            }
          }
        }

      } // end Schwartz inequality
    } // end loop over MU,NU shells

    for (int P=0; P < numPshell; ++P) {
      C_DGEMM('T', 'N', nmoA_, nso_, nso_, 1.0, CA_[0], nmoA_, AO_RI[P], nso_,
        0.0, halftrans, nso_);
      C_DGEMM('N', 'N', nmoA_, nmoA_, nso_, 1.0, halftrans, nso_, CA_[0],
        nmoA_, 0.0, MO_RI[P], nmoA_);
    }

    psio_->write(PSIF_SAPT_TEMP,"MO AA RI Integrals",(char *) &(MO_RI[0][0]),
      sizeof(double)*numPshell*nmoA_*nmoA_,next_DF_MO,&next_DF_MO);

  }

  free_block(AO_RI);
  free(halftrans);
  free_block(MO_RI);

  zero_disk(PSIF_SAPT_AA_DF_INTS,"AA RI Integrals",noccA_*noccA_,ndf_+3);
  zero_disk(PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",noccA_*nvirA_,ndf_+3);
  zero_disk(PSIF_SAPT_AA_DF_INTS,"RR RI Integrals",nvirA_*nvirA_,ndf_+3);

  long int numP;
  long int temp_size = mem_ / (2*ndf_);

  if (temp_size > nmoA_*nmoA_)
    temp_size = nmoA_*nmoA_;

  double** temp;
  double** temp_J;

  psio_address next_DF_AA = PSIO_ZERO;
  psio_address next_DF_AR = PSIO_ZERO;
  psio_address next_DF_RR = PSIO_ZERO;

  for (int i=0,oP=0; oP < nmoA_*nmoA_; ++i, oP += numP) {
    if ((i+1)*temp_size < nmoA_*nmoA_) numP = temp_size;
    else numP = nmoA_*nmoA_ - oP;

    temp = block_matrix(ndf_,numP);

    next_DF_MO = psio_get_address(PSIO_ZERO,sizeof(double)*oP);
    for (int P=0; P < ndf_; ++P) {
      psio_->read(PSIF_SAPT_TEMP,"MO AA RI Integrals",(char *) &(temp[P][0]),
        sizeof(double)*numP,next_DF_MO,&next_DF_MO);
      next_DF_MO = psio_get_address(next_DF_MO,sizeof(double)*
        (nmoA_*nmoA_-numP));
    }

    temp_J = block_matrix(numP,ndf_+3);

    C_DGEMM('T','N',numP,ndf_,ndf_,1.0,temp[0],numP,
      J_mhalf[0],ndf_,0.0,temp_J[0],ndf_+3);

    free_block(temp);

    for (int ij=0; ij<numP; ij++) {
      int i = (ij+oP)/nmoA_;
      int j = (ij+oP)%nmoA_;
      if (i < noccA_ && j < noccA_) {
        psio_->write(PSIF_SAPT_AA_DF_INTS,"AA RI Integrals",(char *)
          &(temp_J[ij][0]),sizeof(double)*(ndf_+3),next_DF_AA,&next_DF_AA);
      }
    }

    for (int ij=0; ij<numP; ij++) {
      int i = (ij+oP)/nmoA_;
      int j = (ij+oP)%nmoA_;
      if (i < noccA_ && j >= noccA_) {
        psio_->write(PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",(char *)
          &(temp_J[ij][0]),sizeof(double)*(ndf_+3),next_DF_AR,&next_DF_AR);
      }
    }

    for (int ij=0; ij<numP; ij++) {
      int i = (ij+oP)/nmoA_;
      int j = (ij+oP)%nmoA_;
      if (i >= noccA_ && j >= noccA_) {
        psio_->write(PSIF_SAPT_AA_DF_INTS,"RR RI Integrals",(char *)
          &(temp_J[ij][0]),sizeof(double)*(ndf_+3),next_DF_RR,&next_DF_RR);
      }
    }

    free_block(temp_J);
  }

  AO_RI = block_matrix(maxPshell,nso_*nso_);
  halftrans = init_array(nmoB_*nso_);
  MO_RI = block_matrix(maxPshell,nmoB_*nmoB_);

  next_DF_MO = PSIO_ZERO;
  psio_address next_bare_BS = PSIO_ZERO;

  for (int Pshell=0; Pshell < ribasis_->nshell(); ++Pshell) {
    int numPshell = ribasis_->shell(Pshell).nfunction();

    #pragma omp for private(mn) schedule(guided)
    for (mn=0; mn < nshelltri; ++mn) {
      #ifdef _OPENMP
        rank = omp_get_thread_num();
      #endif
      int MU = index2i_[mn];
      int NU = index2j_[mn];
      int nummu = basisset_->shell(MU).nfunction();
      int numnu = basisset_->shell(NU).nfunction();

      if (sqrt(Schwartz[mn]*DFSchwartz[Pshell])>schwarz_) {
        eri[rank]->compute_shell(Pshell, 0, MU, NU);

        for (int P=0, index=0; P < numPshell; ++P) {

          for (int mu=0; mu < nummu; ++mu) {
            int omu = basisset_->shell(MU).function_index() + mu;

            for (int nu=0; nu < numnu; ++nu, ++index) {
              int onu = basisset_->shell(NU).function_index() + nu;

              AO_RI[P][omu*nso_+onu] = buffer[rank][index];
              AO_RI[P][onu*nso_+omu] = buffer[rank][index];
            }
          }
        }

      } // end Schwartz inequality
    } // end loop over MU,NU shells

    for (int P=0; P < numPshell; ++P) {
      C_DGEMM('T', 'N', nmoB_, nso_, nso_, 1.0, CB_[0], nmoB_, AO_RI[P], nso_,
        0.0, halftrans, nso_);
      C_DGEMM('N', 'N', nmoB_, nmoB_, nso_, 1.0, halftrans, nso_, CB_[0],
        nmoB_, 0.0, MO_RI[P], nmoB_);
    }

    psio_->write(PSIF_SAPT_TEMP,"MO BB RI Integrals",(char *) &(MO_RI[0][0]),
      sizeof(double)*numPshell*nmoB_*nmoB_,next_DF_MO,&next_DF_MO);

  }

  free_block(AO_RI);
  free(halftrans);
  free_block(MO_RI);

  zero_disk(PSIF_SAPT_BB_DF_INTS,"BB RI Integrals",noccB_*noccB_,ndf_+3);
  zero_disk(PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",noccB_*nvirB_,ndf_+3);
  zero_disk(PSIF_SAPT_BB_DF_INTS,"SS RI Integrals",nvirB_*nvirB_,ndf_+3);

  if (temp_size > nmoB_*nmoB_)
    temp_size = nmoB_*nmoB_;

  psio_address next_DF_BB = PSIO_ZERO;
  psio_address next_DF_BS = PSIO_ZERO;
  psio_address next_DF_SS = PSIO_ZERO;

  for (int i=0,oP=0; oP < nmoB_*nmoB_; ++i, oP += numP) {
    if ((i+1)*temp_size < nmoB_*nmoB_) numP = temp_size;
    else numP = nmoB_*nmoB_ - oP;

    temp = block_matrix(ndf_,numP);

    next_DF_MO = psio_get_address(PSIO_ZERO,sizeof(double)*oP);
    for (int P=0; P < ndf_; ++P) {
      psio_->read(PSIF_SAPT_TEMP,"MO BB RI Integrals",(char *) &(temp[P][0]),
        sizeof(double)*numP,next_DF_MO,&next_DF_MO);
      next_DF_MO = psio_get_address(next_DF_MO,sizeof(double)*(nmoB_*nmoB_
        -numP));
    }

    temp_J = block_matrix(numP,ndf_+3);

    C_DGEMM('T','N',numP,ndf_,ndf_,1.0,temp[0],numP,J_mhalf[0],ndf_,0.0,
      temp_J[0],ndf_+3);

    free_block(temp);

    for (int ij=0; ij<numP; ij++) {
      int i = (ij+oP)/nmoB_;
      int j = (ij+oP)%nmoB_;
      if (i < noccB_ && j < noccB_) {
        psio_->write(PSIF_SAPT_BB_DF_INTS,"BB RI Integrals",(char *)
          &(temp_J[ij][0]),sizeof(double)*(ndf_+3),next_DF_BB,&next_DF_BB);
      }
    }

    for (int ij=0; ij<numP; ij++) {
      int i = (ij+oP)/nmoB_;
      int j = (ij+oP)%nmoB_;
      if (i < noccB_ && j >= noccB_) {
        psio_->write(PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",(char *)
          &(temp_J[ij][0]),sizeof(double)*(ndf_+3),next_DF_BS,&next_DF_BS);
      }
    }

    for (int ij=0; ij<numP; ij++) {
      int i = (ij+oP)/nmoB_;
      int j = (ij+oP)%nmoB_;
      if (i >= noccB_ && j >= noccB_) {
        psio_->write(PSIF_SAPT_BB_DF_INTS,"SS RI Integrals",(char *)
          &(temp_J[ij][0]),sizeof(double)*(ndf_+3),next_DF_SS,&next_DF_SS);
      }
    }

    free_block(temp_J);
  }

  AO_RI = block_matrix(maxPshell,nso_*nso_);
  halftrans = init_array(nmoA_*nso_);
  MO_RI = block_matrix(maxPshell,nmoA_*nmoB_);

  next_DF_MO = PSIO_ZERO;

  for (int Pshell=0; Pshell < ribasis_->nshell(); ++Pshell) {
    int numPshell = ribasis_->shell(Pshell).nfunction();

    #pragma omp for private(mn) schedule(guided)
    for (mn=0; mn < nshelltri; ++mn) {
      #ifdef _OPENMP
        rank = omp_get_thread_num();
      #endif
      int MU = index2i_[mn];
      int NU = index2j_[mn];
      int nummu = basisset_->shell(MU).nfunction();
      int numnu = basisset_->shell(NU).nfunction();

      if (sqrt(Schwartz[mn]*DFSchwartz[Pshell])>schwarz_) {
        eri[rank]->compute_shell(Pshell, 0, MU, NU);

        for (int P=0, index=0; P < numPshell; ++P) {

          for (int mu=0; mu < nummu; ++mu) {
            int omu = basisset_->shell(MU).function_index() + mu;

            for (int nu=0; nu < numnu; ++nu, ++index) {
              int onu = basisset_->shell(NU).function_index() + nu;

              AO_RI[P][omu*nso_+onu] = buffer[rank][index];
              AO_RI[P][onu*nso_+omu] = buffer[rank][index];
            }
          }
        }

      } // end Schwartz inequality
    } // end loop over MU,NU shells

    for (int P=0; P < numPshell; ++P) {
      C_DGEMM('T', 'N', nmoA_, nso_, nso_, 1.0, CA_[0], nmoA_, AO_RI[P], nso_,
        0.0, halftrans, nso_);
      C_DGEMM('N', 'N', nmoA_, nmoB_, nso_, 1.0, halftrans, nso_, CB_[0],
        nmoB_, 0.0, MO_RI[P], nmoB_);
    }

    psio_->write(PSIF_SAPT_TEMP,"MO AB RI Integrals",(char *) &(MO_RI[0][0]),
      sizeof(double)*numPshell*nmoA_*nmoB_,next_DF_MO,&next_DF_MO);

  }

  free_block(AO_RI);
  free_block(MO_RI);

  zero_disk(PSIF_SAPT_AB_DF_INTS,"AB RI Integrals",noccA_*noccB_,ndf_+3);
  zero_disk(PSIF_SAPT_AB_DF_INTS,"AS RI Integrals",noccA_*nvirB_,ndf_+3);
  zero_disk(PSIF_SAPT_AB_DF_INTS,"RB RI Integrals",nvirA_*noccB_,ndf_+3);

  if (temp_size > nmoA_*nmoB_)
    temp_size = nmoA_*nmoB_;

  psio_address next_DF_AB = PSIO_ZERO;
  psio_address next_DF_AS = PSIO_ZERO;
  psio_address next_DF_RB = PSIO_ZERO;

  for (int i=0,oP=0; oP < nmoA_*nmoB_; ++i, oP += numP) {
    if ((i+1)*temp_size < nmoA_*nmoB_) numP = temp_size;
    else numP = nmoA_*nmoB_ - oP;

    temp = block_matrix(ndf_,numP);

    next_DF_MO = psio_get_address(PSIO_ZERO,sizeof(double)*oP);
    for (int P=0; P < ndf_; ++P) {
      psio_->read(PSIF_SAPT_TEMP,"MO AB RI Integrals",(char *) &(temp[P][0]),
        sizeof(double)*numP,next_DF_MO,&next_DF_MO);
      next_DF_MO = psio_get_address(next_DF_MO,sizeof(double)*(nmoA_*nmoB_
        -numP));
    }

    temp_J = block_matrix(numP,ndf_+3);

    C_DGEMM('T','N',numP,ndf_,ndf_,1.0,temp[0],numP,J_mhalf[0],ndf_,
      0.0,temp_J[0],ndf_+3);

    free_block(temp);

    for (int ij=0; ij<numP; ij++) {
      int i = (ij+oP)/nmoB_;
      int j = (ij+oP)%nmoB_;
      if (i < noccA_ && j < noccB_) {
        psio_->write(PSIF_SAPT_AB_DF_INTS,"AB RI Integrals",(char *)
          &(temp_J[ij][0]),sizeof(double)*(ndf_+3),next_DF_AB,&next_DF_AB);
      }
    }

    for (int ij=0; ij<numP; ij++) {
      int i = (ij+oP)/nmoB_;
      int j = (ij+oP)%nmoB_;
      if (i < noccA_ && j >= noccB_) {
        psio_->write(PSIF_SAPT_AB_DF_INTS,"AS RI Integrals",(char *)
          &(temp_J[ij][0]),sizeof(double)*(ndf_+3),next_DF_AS,&next_DF_AS);
      }
    }

    for (int ij=0; ij<numP; ij++) {
      int i = (ij+oP)/nmoB_;
      int j = (ij+oP)%nmoB_;
      if (i >= noccA_ && j < noccB_) {
        psio_->write(PSIF_SAPT_AB_DF_INTS,"RB RI Integrals",(char *)
          &(temp_J[ij][0]),sizeof(double)*(ndf_+3),next_DF_RB,&next_DF_RB);
      }
    }

    free_block(temp_J);
  }

  free(halftrans);

  free_block(J_mhalf);

  psio_->close(PSIF_SAPT_TEMP,0);

  free(Schwartz);
  free(DFSchwartz);

  fflush(outfile);
}

void SAPT2::w_integrals()
{
  double **B_p_A = get_diag_AA_ints(1);
  double **B_p_B = get_diag_BB_ints(1);

  diagAA_ = init_array(ndf_+3);
  diagBB_ = init_array(ndf_+3);

  for(int a=0; a<noccA_; a++){
    C_DAXPY(ndf_+3,1.0,&(B_p_A[a][0]),1,diagAA_,1);
  }

  for(int b=0; b<noccB_; b++){
    C_DAXPY(ndf_+3,1.0,&(B_p_B[b][0]),1,diagBB_,1);
  }

  free_block(B_p_A);
  free_block(B_p_B);

  wBAR_ = block_matrix(noccA_,nvirA_);

  for(int a=0; a<noccA_; a++){
    C_DAXPY(nvirA_,1.0,&(vBAA_[a][noccA_]),1,&(wBAR_[a][0]),1);
  }

  double **B_p_AR = get_AR_ints(0);

  C_DGEMV('n',noccA_*nvirA_,ndf_,2.0,&(B_p_AR[0][0]),ndf_+3,diagBB_,1,1.0,
    &(wBAR_[0][0]),1);

  free_block(B_p_AR);

  wABS_ = block_matrix(noccB_,nvirB_);

  for(int b=0; b<noccB_; b++){
    C_DAXPY(nvirB_,1.0,&(vABB_[b][noccB_]),1,&(wABS_[b][0]),1);
  }

  double **B_p_BS = get_BS_ints(0);

  C_DGEMV('n',noccB_*nvirB_,ndf_,2.0,&(B_p_BS[0][0]),ndf_+3,diagAA_,1,1.0,
    &(wABS_[0][0]),1);

  free_block(B_p_BS);

  wBAA_ = block_matrix(noccA_,noccA_);

  for(int a=0; a<noccA_; a++){
    C_DAXPY(noccA_,1.0,&(vBAA_[a][0]),1,&(wBAA_[a][0]),1);
  }

  double **B_p_AA = get_AA_ints(0);

  C_DGEMV('n',noccA_*noccA_,ndf_,2.0,&(B_p_AA[0][0]),ndf_+3,diagBB_,1,1.0,
    &(wBAA_[0][0]),1);

  free_block(B_p_AA);

  wABB_ = block_matrix(noccB_,noccB_);

  for(int b=0; b<noccB_; b++){
    C_DAXPY(noccB_,1.0,&(vABB_[b][0]),1,&(wABB_[b][0]),1);
  }

  double **B_p_BB = get_BB_ints(0);

  C_DGEMV('n',noccB_*noccB_,ndf_,2.0,&(B_p_BB[0][0]),ndf_+3,diagAA_,1,1.0,
    &(wABB_[0][0]),1);

  free_block(B_p_BB);

  wBRR_ = block_matrix(nvirA_,nvirA_);

  for(int r=0; r<nvirA_; r++){
    C_DAXPY(nvirA_,1.0,&(vBAA_[r+noccA_][noccA_]),1,&(wBRR_[r][0]),1);
  }

  double **B_p_RR = get_RR_ints(0);

  C_DGEMV('n',nvirA_*nvirA_,ndf_,2.0,&(B_p_RR[0][0]),ndf_+3,diagBB_,1,1.0,
    &(wBRR_[0][0]),1);

  free_block(B_p_RR);

  wASS_ = block_matrix(nvirB_,nvirB_);

  for(int s=0; s<nvirB_; s++){
    C_DAXPY(nvirB_,1.0,&(vABB_[s+noccB_][noccB_]),1,&(wASS_[s][0]),1);
  }

  double **B_p_SS = get_SS_ints(0);

  C_DGEMV('n',nvirB_*nvirB_,ndf_,2.0,&(B_p_SS[0][0]),ndf_+3,diagAA_,1,1.0,
    &(wASS_[0][0]),1);

  free_block(B_p_SS);
}

void SAPT2::natural_orbitalify(int ampfile, const char *VV_opdm,
  double *evals, int foccA, int noccA, int nvirA, const char monomer)
{
  double **P = block_matrix(nvirA,nvirA);

  psio_->read_entry(ampfile,VV_opdm,(char *) P[0],
    sizeof(double)*nvirA*nvirA);

  C_DSCAL(nvirA*nvirA,2.0,P[0],1);

  double *occnum = init_array(nvirA);
  double **nat_orbs_MO = block_matrix(nvirA,nvirA);

  sq_rsp(nvirA,nvirA,P,occnum,3,nat_orbs_MO,1.0e-14);

  int num_no_vir = 0;

  for (int i=0; i<nvirA; i++) {
    if (occnum[i] > occ_cutoff_) {
      num_no_vir++;
    }
    else break;
  }

  if (print_) {
    fprintf(outfile,"    Monomer %c: %d virtual orbitals dropped\n",monomer,
          nvirA-num_no_vir);
    fflush(outfile);
  }

  double **Fock_MO = block_matrix(nvirA,nvirA);

  for (int i=0; i<nvirA; i++) {
    Fock_MO[i][i] = evals[i+noccA];
  }

  double **tempmat = block_matrix(num_no_vir,nvirA);
  double **Fock_NO = block_matrix(num_no_vir,num_no_vir);

  C_DGEMM('T','N',num_no_vir,nvirA,nvirA,1.0,&(nat_orbs_MO[0][0]),nvirA,
    &(Fock_MO[0][0]),nvirA,0.0,&(tempmat[0][0]),nvirA);
  C_DGEMM('N','N',num_no_vir,num_no_vir,nvirA,1.0,&(tempmat[0][0]),nvirA,
    &(nat_orbs_MO[0][0]),nvirA,0.0,&(Fock_NO[0][0]),num_no_vir);

  double *epsilon = init_array(num_no_vir);
  double **X = block_matrix(num_no_vir,num_no_vir);
  sq_rsp(num_no_vir,num_no_vir,Fock_NO,epsilon,1,X,1.0e-14);

  double **MO_MVO = block_matrix(nvirA,num_no_vir);

  C_DGEMM('N','N',nvirA,num_no_vir,num_no_vir,1.0,&(nat_orbs_MO[0][0]),
    nvirA,&(X[0][0]),num_no_vir,0.0,&(MO_MVO[0][0]),num_no_vir);

  if (monomer == 'A') {
    no_CA_ = block_matrix(nvirA,num_no_vir);
    no_evalsA_ = init_array(noccA+num_no_vir);
    no_nvirA_ = num_no_vir;

    C_DCOPY(nvirA*num_no_vir,&(MO_MVO[0][0]),1,
      &(no_CA_[0][0]),1);
    C_DCOPY(noccA,evals,1,no_evalsA_,1);
    C_DCOPY(num_no_vir,epsilon,1,&(no_evalsA_[noccA]),1);
  }

  if (monomer == 'B') {
    no_CB_ = block_matrix(nvirA,num_no_vir);
    no_evalsB_ = init_array(noccA+num_no_vir);
    no_nvirB_ = num_no_vir;

    C_DCOPY(nvirA*num_no_vir,&(MO_MVO[0][0]),1,
      &(no_CB_[0][0]),1);
    C_DCOPY(noccA,evals,1,no_evalsB_,1);
    C_DCOPY(num_no_vir,epsilon,1,&(no_evalsB_[noccA]),1);
  }

  free(epsilon);
  free(occnum);
  free_block(P);
  free_block(nat_orbs_MO);
  free_block(tempmat);
  free_block(Fock_MO);
  free_block(Fock_NO);
  free_block(X);
  free_block(MO_MVO);
}

void SAPT2::natural_orbitalify_df_ints()
{
  double **B_p_AR = get_DF_ints(PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",
    0,noccA_,0,nvirA_);
  double **C_p_AR = block_matrix(noccA_*no_nvirA_,ndf_+3);

  for (int a=0; a<noccA_; a++) {
    C_DGEMM('T','N',no_nvirA_,ndf_+3,nvirA_,1.0,no_CA_[0],no_nvirA_,
      B_p_AR[a*nvirA_],ndf_+3,0.0,C_p_AR[a*no_nvirA_],ndf_+3);
  }

  psio_->write_entry(PSIF_SAPT_AA_DF_INTS,"AR NO RI Integrals",
    (char *) C_p_AR[0],sizeof(double)*noccA_*no_nvirA_*(ndf_+3));

  free_block(B_p_AR);
  free_block(C_p_AR);

  double **B_p_BS = get_DF_ints(PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",
    0,noccB_,0,nvirB_);
  double **C_p_BS = block_matrix(noccB_*no_nvirB_,ndf_+3);

  for (int b=0; b<noccB_; b++) {
    C_DGEMM('T','N',no_nvirB_,ndf_+3,nvirB_,1.0,no_CB_[0],no_nvirB_,
      B_p_BS[b*nvirB_],ndf_+3,0.0,C_p_BS[b*no_nvirB_],ndf_+3);
  }

  psio_->write_entry(PSIF_SAPT_BB_DF_INTS,"BS NO RI Integrals",
    (char *) C_p_BS[0],sizeof(double)*noccB_*no_nvirB_*(ndf_+3));

  free_block(B_p_BS);
  free_block(C_p_BS);

  double **B_p_RR = get_DF_ints(PSIF_SAPT_AA_DF_INTS,"RR RI Integrals",
    0,nvirA_,0,nvirA_);
  double **C_p_RR = block_matrix(no_nvirA_*nvirA_,ndf_+3);

  C_DGEMM('T','N',no_nvirA_,nvirA_*(ndf_+3),nvirA_,1.0,no_CA_[0],no_nvirA_,
    B_p_RR[0],nvirA_*(ndf_+3),0.0,C_p_RR[0],nvirA_*(ndf_+3));

  free_block(B_p_RR);

  double **D_p_RR = block_matrix(no_nvirA_*no_nvirA_,ndf_+3);

  for (int r=0; r<no_nvirA_; r++) {
    C_DGEMM('T','N',no_nvirA_,ndf_+3,nvirA_,1.0,no_CA_[0],no_nvirA_,
      C_p_RR[r*nvirA_],ndf_+3,0.0,D_p_RR[r*no_nvirA_],ndf_+3);
  }

  psio_->write_entry(PSIF_SAPT_AA_DF_INTS,"RR NO RI Integrals",
    (char *) D_p_RR[0],sizeof(double)*no_nvirA_*no_nvirA_*(ndf_+3));

  free_block(C_p_RR);
  free_block(D_p_RR);

  double **B_p_SS = get_DF_ints(PSIF_SAPT_BB_DF_INTS,"SS RI Integrals",
    0,nvirB_,0,nvirB_);
  double **C_p_SS = block_matrix(no_nvirB_*nvirB_,ndf_+3);

  C_DGEMM('T','N',no_nvirB_,nvirB_*(ndf_+3),nvirB_,1.0,no_CB_[0],no_nvirB_,
    B_p_SS[0],nvirB_*(ndf_+3),0.0,C_p_SS[0],nvirB_*(ndf_+3));

  free_block(B_p_SS);

  double **D_p_SS = block_matrix(no_nvirB_*no_nvirB_,ndf_+3);

  for (int s=0; s<no_nvirB_; s++) {
    C_DGEMM('T','N',no_nvirB_,ndf_+3,nvirB_,1.0,no_CB_[0],no_nvirB_,
      C_p_SS[s*nvirB_],ndf_+3,0.0,D_p_SS[s*no_nvirB_],ndf_+3);
  }

  psio_->write_entry(PSIF_SAPT_BB_DF_INTS,"SS NO RI Integrals",
    (char *) D_p_SS[0],sizeof(double)*no_nvirB_*no_nvirB_*(ndf_+3));

  free_block(C_p_SS);
  free_block(D_p_SS);
}

}}

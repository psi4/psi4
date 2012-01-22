#include "sapt0.h"

namespace psi { namespace sapt {

SAPT0::SAPT0(Options& options, boost::shared_ptr<PSIO> psio, 
  boost::shared_ptr<Chkpt> chkpt) : SAPT(options, psio, chkpt),
  e_elst10_(0.0),
  e_exch10_(0.0),
  e_exch10_s2_(0.0),
  e_ind20_(0.0),
  e_exch_ind20_(0.0),
  e_disp20_(0.0),
  e_exch_disp20_(0.0),
  e_disp20_ss_(0.0),
  e_disp20_os_(0.0),
  e_exch_disp20_ss_(0.0),
  e_exch_disp20_os_(0.0),
  e_sapt0_(0.0),
  e_sapt0_scs_(0.0)
{
  print_header();

  maxiter_ = options_.get_int("MAXITER");
  e_conv_ = options_.get_double("E_CONVERGENCE");
  d_conv_ = options_.get_double("D_CONVERGENCE");
  no_response_ = options_.get_bool("NO_RESPONSE");
  aio_cphf_ = options_.get_bool("AIO_CPHF");
  aio_dfints_ = options_.get_bool("AIO_DF_INTS");

  wBAR_ = NULL;
  wABS_ = NULL;
}

SAPT0::~SAPT0()
{
  if (wBAR_ != NULL) free_block(wBAR_);
  if (wABS_ != NULL) free_block(wABS_);
  psio_->close(PSIF_SAPT_AA_DF_INTS,1);
  psio_->close(PSIF_SAPT_BB_DF_INTS,1);
  psio_->close(PSIF_SAPT_AB_DF_INTS,1);
}

double SAPT0::compute_energy()
{
  check_memory();

  if (elst_basis_) 
    first_order_terms();

  psio_->open(PSIF_SAPT_AA_DF_INTS,PSIO_OPEN_NEW);
  psio_->open(PSIF_SAPT_BB_DF_INTS,PSIO_OPEN_NEW);
  psio_->open(PSIF_SAPT_AB_DF_INTS,PSIO_OPEN_NEW);

  timer_on("DF Integrals       ");
    if (aio_dfints_)
      df_integrals_aio();
    else
      df_integrals();
  timer_off("DF Integrals       ");
  timer_on("W Integrals        ");
    w_integrals();
  timer_off("W Integrals        ");
  if (!elst_basis_) {
    timer_on("Elst10             ");
      elst10();
    timer_off("Elst10             ");
    timer_on("Exch10             ");
      exch10();
    timer_off("Exch10             ");
    timer_on("Exch10 S^2         ");
      exch10_s2();
    timer_off("Exch10 S^2         ");
  }
  timer_on("Ind20              ");
    if (debug_ || no_response_) ind20();
    if (!no_response_) ind20r();
  timer_off("Ind20              ");
  timer_on("Exch-Ind20         ");
    exch_ind20A_B();
    exch_ind20B_A();
  timer_off("Exch-Ind20         ");
  if (debug_) disp20();
  timer_on("Exch-Disp20 N^5    ");
    psio_->open(PSIF_SAPT_TEMP,PSIO_OPEN_NEW);
    exch_disp20_n5();
  timer_off("Exch-Disp20 N^5    ");
  timer_on("Exch-Disp20 N^4    ");
    exch_disp20_n4();
    psio_->close(PSIF_SAPT_TEMP,0);
  timer_off("Exch-Disp20 N^4    ");

  print_results();

  return (e_sapt0_);
}

void SAPT0::print_header()
{
  fprintf(outfile,"        SAPT0  \n");
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
  if (elst_basis_) {
    fprintf(outfile,"    NRI        = %9d\n",ribasis_->nbf());
    fprintf(outfile,"    NRI (Elst) = %9d\n",elstbasis_->nbf());
  }
  else {
    fprintf(outfile,"    NRI        = %9d\n",ndf_);
  }
  fprintf(outfile,"    NOCC A     = %9d\n",noccA_);
  fprintf(outfile,"    NOCC B     = %9d\n",noccB_);
  fprintf(outfile,"    FOCC A     = %9d\n",foccA_);
  fprintf(outfile,"    FOCC B     = %9d\n",foccB_);
  fprintf(outfile,"    NVIR A     = %9d\n",nvirA_);
  fprintf(outfile,"    NVIR B     = %9d\n",nvirB_);
  fprintf(outfile,"\n");
  fflush(outfile);
}

void SAPT0::print_results()
{
  e_sapt0_ = eHF_ + e_disp20_ + e_exch_disp20_;
  double SOS = options_.get_double("SAPT_OS_SCALE");
  double SSS = options_.get_double("SAPT_SS_SCALE");
  e_sapt0_scs_ = eHF_ + SOS*(e_disp20_os_ + e_exch_disp20_os_) 
    + SSS*(e_disp20_ss_ + e_exch_disp20_ss_);
  double dHF = eHF_ - (e_elst10_ + e_exch10_ + e_ind20_ + e_exch_ind20_);

  double tot_elst = e_elst10_;
  double tot_exch = e_exch10_;
  double tot_ind = e_ind20_ + e_exch_ind20_ + dHF;
  double tot_disp = e_disp20_ + e_exch_disp20_;
  double tot_scs_disp = e_sapt0_scs_ - eHF_;

  fprintf(outfile,"\n    SAPT Results  \n");
  fprintf(outfile,"  -----------------------------------------------------------------------\n");
  fprintf(outfile,"    Electrostatics     %16.8lf mH %16.8lf kcal mol^-1\n",
    tot_elst*1000.0,tot_elst*627.5095);
  fprintf(outfile,"      Elst10,r         %16.8lf mH %16.8lf kcal mol^-1\n\n",
    e_elst10_*1000.0,e_elst10_*627.5095);
  fprintf(outfile,"    Exchange           %16.8lf mH %16.8lf kcal mol^-1\n",
    tot_exch*1000.0,tot_exch*627.5095);
  fprintf(outfile,"      Exch10           %16.8lf mH %16.8lf kcal mol^-1\n",
    e_exch10_*1000.0,e_exch10_*627.5095);
  fprintf(outfile,"      Exch10(S^2)      %16.8lf mH %16.8lf kcal mol^-1\n\n",
    e_exch10_s2_*1000.0,e_exch10_s2_*627.5095);
  fprintf(outfile,"    Induction          %16.8lf mH %16.8lf kcal mol^-1\n",
    tot_ind*1000.0,tot_ind*627.5095);
  if (no_response_) {
    fprintf(outfile,"      Ind20            %16.8lf mH %16.8lf kcal mol^-1\n",
      e_ind20_*1000.0,e_ind20_*627.5095);
    fprintf(outfile,"      Exch-Ind20       %16.8lf mH %16.8lf kcal mol^-1\n",
      e_exch_ind20_*1000.0,e_exch_ind20_*627.5095);
    fprintf(outfile,"      delta HF         %16.8lf mH %16.8lf kcal mol^-1\n\n",
      dHF*1000.0,dHF*627.5095);
  } else {
    fprintf(outfile,"      Ind20,r          %16.8lf mH %16.8lf kcal mol^-1\n",
      e_ind20_*1000.0,e_ind20_*627.5095);
    fprintf(outfile,"      Exch-Ind20,r     %16.8lf mH %16.8lf kcal mol^-1\n",
      e_exch_ind20_*1000.0,e_exch_ind20_*627.5095);
    fprintf(outfile,"      delta HF,r       %16.8lf mH %16.8lf kcal mol^-1\n\n",
      dHF*1000.0,dHF*627.5095);
  }
  fprintf(outfile,"    Dispersion         %16.8lf mH %16.8lf kcal mol^-1\n",
    tot_disp*1000.0,tot_disp*627.5095);
  fprintf(outfile,"      Disp20           %16.8lf mH %16.8lf kcal mol^-1\n",
    e_disp20_*1000.0,e_disp20_*627.5095);
  fprintf(outfile,"      Exch-Disp20      %16.8lf mH %16.8lf kcal mol^-1\n\n",
    e_exch_disp20_*1000.0,e_exch_disp20_*627.5095);
  fprintf(outfile,"    SCS Dispersion     %16.8lf mH %16.8lf kcal mol^-1\n",
    tot_scs_disp*1000.0,tot_scs_disp*627.5095);
  fprintf(outfile,"      Disp20 (SS)      %16.8lf mH %16.8lf kcal mol^-1\n",
    e_disp20_ss_*1000.0,e_disp20_ss_*627.5095);
  fprintf(outfile,"      Disp20 (OS)      %16.8lf mH %16.8lf kcal mol^-1\n",
    e_disp20_os_*1000.0,e_disp20_os_*627.5095);
  fprintf(outfile,"      Exch-Disp20 (SS) %16.8lf mH %16.8lf kcal mol^-1\n",
    e_exch_disp20_ss_*1000.0,e_exch_disp20_ss_*627.5095);
  fprintf(outfile,"      Exch-Disp20 (OS) %16.8lf mH %16.8lf kcal mol^-1\n\n",
    e_exch_disp20_os_*1000.0,e_exch_disp20_os_*627.5095);

  fprintf(outfile,"    Same-Spin Scale        %11.3E\n", SSS);
  fprintf(outfile,"    Opposite-Spin Scale    %11.3E\n\n", SOS);

  fprintf(outfile,"    Total HF           %16.8lf mH %16.8lf kcal mol^-1\n",
    eHF_*1000.0,eHF_*627.5095);
  fprintf(outfile,"    Total SAPT0        %16.8lf mH %16.8lf kcal mol^-1\n",
    e_sapt0_*1000.0,e_sapt0_*627.5095);
  fprintf(outfile,"    Total SCS-SAPT0    %16.8lf mH %16.8lf kcal mol^-1\n",
    e_sapt0_scs_*1000.0,e_sapt0_scs_*627.5095);

  Process::environment.globals["SAPT ELST ENERGY"] = tot_elst;
  Process::environment.globals["SAPT EXCH ENERGY"] = tot_exch;
  Process::environment.globals["SAPT IND ENERGY"] = tot_ind;
  Process::environment.globals["SAPT CT ENERGY"] = e_ind20_ + e_exch_ind20_;
  Process::environment.globals["SAPT DISP ENERGY"] = tot_disp;
  Process::environment.globals["SAPT SCS-DISP ENERGY"] = tot_scs_disp;
  Process::environment.globals["SAPT SAPT0 ENERGY"] = e_sapt0_;
  Process::environment.globals["SAPT SCS-SAPT0 ENERGY"] = e_sapt0_scs_;
  Process::environment.globals["SAPT ENERGY"] = e_sapt0_;
  Process::environment.globals["CURRENT ENERGY"] = Process::environment.globals["SAPT ENERGY"];
}

void SAPT0::check_memory()
{
  double memory = 8.0*mem_/1000000.0;

  if (debug_) {
    fprintf(outfile,"    Using %8.1lf MB Memory\n\n",memory);
    fflush(outfile);
  }

  bool fail = false;

  int max_func_per_shell = basisset_->max_function_per_shell(); 
  int nsotri = nso_*(nso_+1)/2;

  long int dfint = ndf_*ndf_ + 2*max_func_per_shell*ndf_;
  long int indices = nsotri + noccA_*noccA_ + noccA_*nvirA_ + nvirA_*nvirA_
    + noccB_*noccB_ + noccB_*nvirB_ + nvirB_*nvirB_ + noccB_*noccB_
    + noccA_*nvirB_ + noccB_*nvirA_;
  long int exchdisp = 2L*nvirB_*(ndf_+3) + 1L*nvirA_*(ndf_) 
    + 2L*nvirA_*(ndf_+3) + 1L*nvirB_*(ndf_);

  if (dfint > mem_) fail = true;
  if (indices > mem_) fail = true;
  if (exchdisp > mem_) fail = true;

  if (fail) throw PsiException("Not enough memory", __FILE__,__LINE__);
}

void SAPT0::first_order_terms()
{
  ndf_ = elstbasis_->nbf();

  psio_->open(PSIF_SAPT_AA_DF_INTS,PSIO_OPEN_NEW);
  psio_->open(PSIF_SAPT_BB_DF_INTS,PSIO_OPEN_NEW);
  psio_->open(PSIF_SAPT_AB_DF_INTS,PSIO_OPEN_NEW);

  timer_on("OO DF Integrals    ");
    oo_df_integrals();
  timer_off("OO DF Integrals    ");
  timer_on("Elst10             ");
    elst10();
  timer_off("Elst10             ");
  timer_on("Exch10             ");
    exch10();
  timer_off("Exch10             ");
  timer_on("Exch10 S^2         ");
    exch10_s2();
  timer_off("Exch10 S^2         ");

  psio_->close(PSIF_SAPT_AA_DF_INTS,1);
  psio_->close(PSIF_SAPT_BB_DF_INTS,1);
  psio_->close(PSIF_SAPT_AB_DF_INTS,1);

  free(diagAA_);
  free(diagBB_);

  ndf_ = ribasis_->nbf();
}

void SAPT0::df_integrals()
{
  psio_->open(PSIF_SAPT_TEMP,PSIO_OPEN_NEW);

  // Get fitting metric
  boost::shared_ptr<FittingMetric> metric = boost::shared_ptr<FittingMetric>(
    new FittingMetric(ribasis_));
  metric->form_eig_inverse();
  double **J_temp = metric->get_metric()->pointer();
  double **J_mhalf = block_matrix(ndf_,ndf_);
  C_DCOPY(ndf_*ndf_,J_temp[0],1,J_mhalf[0],1);
  metric.reset();

  // Get Schwartz screening arrays
  double maxSchwartz = 0.0;
  double *Schwartz = init_array(basisset_->nshell()*(basisset_->nshell()+1)/2);

  boost::shared_ptr<IntegralFactory> ao_eri_factory = 
    boost::shared_ptr<IntegralFactory>(new IntegralFactory(basisset_, 
    basisset_, basisset_, basisset_));
  boost::shared_ptr<TwoBodyAOInt> ao_eri = boost::shared_ptr<TwoBodyAOInt>(
    ao_eri_factory->eri());
  const double *ao_buffer = ao_eri->buffer();

  for(int P=0,PQ=0;P<basisset_->nshell();P++) {
    int numw = basisset_->shell(P)->nfunction();
    for(int Q=0;Q<=P;Q++,PQ++) {
      int numx = basisset_->shell(Q)->nfunction();
      double tei, max=0.0;

      ao_eri->compute_shell(P, Q, P, Q);

      for(int w=0;w<numw;w++) {
        for(int x=0;x<numx;x++) {
          int index = ( ( (w*numx + x) * numw + w) * numx + x);
          tei = ao_buffer[index];
          if(fabs(tei) > max) max = fabs(tei);
        }
      }
      Schwartz[PQ] = max;
      if (max > maxSchwartz) maxSchwartz = max;
    }
  }

  ao_eri.reset();
  ao_eri_factory.reset();

  double *DFSchwartz = init_array(ribasis_->nshell());

  boost::shared_ptr<IntegralFactory> df_eri_factory = 
    boost::shared_ptr<IntegralFactory>(new IntegralFactory(ribasis_, zero_, 
    ribasis_, zero_));
  boost::shared_ptr<TwoBodyAOInt> df_eri = boost::shared_ptr<TwoBodyAOInt>(
    df_eri_factory->eri());
  const double *df_buffer = df_eri->buffer();

  for(int P=0;P<ribasis_->nshell();P++) {
    int numw = ribasis_->shell(P)->nfunction();
    double tei, max=0.0;

    df_eri->compute_shell(P, 0, P, 0);

    for(int w=0;w<numw;w++) {
      tei = df_buffer[w];
      if(fabs(tei) > max) max = fabs(tei);
    }
    DFSchwartz[P] = max;
  }

  df_eri.reset();
  df_eri_factory.reset();

  long int nsotri = nso_*(nso_+1)/2;
  long int nsotri_screened = 0;

  for(int P=0,PQ=0;P<basisset_->nshell();P++) {
    int numw = basisset_->shell(P)->nfunction();
    for(int Q=0;Q<=P;Q++,PQ++) {
      int numx = basisset_->shell(Q)->nfunction();
      int numPQ = numw*numx;
      if (P == Q) numPQ = numw*(numw+1)/2;
      if (sqrt(Schwartz[PQ]*maxSchwartz)>schwarz_)
        nsotri_screened += numPQ;
  }}

  long int avail_mem = mem_ - (long int) ndf_*ndf_; 
  long int mem_tot = (long int) 2*ndf_*nsotri_screened + (long int) ndf_*ndf_;
  if (avail_mem < (long int) 2*ndf_)
    throw PsiException("Not enough memory", __FILE__,__LINE__);
  long int max_size = avail_mem / ((long int) 2*ndf_);
  if (max_size > nsotri_screened)
    max_size = nsotri_screened;

  if (debug_) {
    fprintf(outfile,"Requires storage of %ld doubles\n",mem_tot);
    fprintf(outfile,"Max nso x nso block is %ld\n\n",max_size);
    fflush(outfile);
  }

  int size = 0;
  int num_blocks = 1;

  for(int P=0,PQ=0;P<basisset_->nshell();P++) {
    int numw = basisset_->shell(P)->nfunction();
    for(int Q=0;Q<=P;Q++,PQ++) {
      int numx = basisset_->shell(Q)->nfunction();
      int numPQ = numw*numx;
      if (P == Q) numPQ = numw*(numw+1)/2;

      if (sqrt(Schwartz[PQ]*maxSchwartz)>schwarz_) {
        size += numPQ;
        if (max_size < size) {
          if (debug_)
            fprintf(outfile,"Block %d : %d\n",num_blocks,size-numPQ);
          num_blocks++;
          size = numPQ;
        }
      }
  }} 

  if (debug_)
    fprintf(outfile,"Block %d : %d\n\n",num_blocks,size);

  int *PQ_start = init_int_array(num_blocks);
  int *PQ_stop = init_int_array(num_blocks);
  int *block_length =  init_int_array(num_blocks);

  int block_num = 0;
  size = 0;

  PQ_start[0] = 0;

  for(int P=0,PQ=0;P<basisset_->nshell();P++) {
    int numw = basisset_->shell(P)->nfunction();
    for(int Q=0;Q<=P;Q++,PQ++) {
      int numx = basisset_->shell(Q)->nfunction();
      int numPQ = numw*numx;
      if (P == Q) numPQ = numw*(numw+1)/2;
      if (sqrt(Schwartz[PQ]*maxSchwartz)>schwarz_) {
        size += numPQ;
        if (max_size < size) {
          PQ_stop[block_num] = PQ - 1;
          block_length[block_num] = size-numPQ;
          block_num++;
          PQ_start[block_num] = PQ;
          size = numPQ;
        }
      }
  }} 

  PQ_stop[num_blocks-1] = basisset_->nshell()*(basisset_->nshell()+1)/2;
  block_length[num_blocks-1] = size;

  int max_func_per_shell = basisset_->max_function_per_shell();

  if (max_func_per_shell*max_func_per_shell > max_size) 
    throw PsiException("Not enough memory", __FILE__,__LINE__);

  if (debug_) {
    for (int i=0; i<num_blocks; i++) 
      fprintf(outfile,"Block %2d : PQ %4d - %4d : %d\n",i,PQ_start[i],
        PQ_stop[i],block_length[i]);
    fprintf(outfile,"\n");
    fflush(outfile);
  }

  boost::shared_ptr<IntegralFactory> rifactory = 
    boost::shared_ptr<IntegralFactory>(new IntegralFactory(ribasis_, zero_, 
    basisset_, basisset_));

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

  zero_disk(PSIF_SAPT_TEMP,"AO RI Integrals",ndf_,nsotri_screened);

  psio_address next_DF_AO = PSIO_ZERO;

  double** AO_RI = block_matrix(max_size,ndf_);
  double** J_AO_RI = block_matrix(ndf_,max_size);

  int munu_offset = 0;
  int curr_block = 0;
  int offset = 0;
  int Pshell;

  for(int MU=0,MUNU=0;MU<basisset_->nshell();MU++) {
    int nummu = basisset_->shell(MU)->nfunction();
    for(int NU=0;NU<=MU;NU++,MUNU++) {
      int numnu = basisset_->shell(NU)->nfunction();

      if (sqrt(Schwartz[MUNU]*maxSchwartz)>schwarz_ ) {

#pragma omp parallel
{
        #pragma omp for private(Pshell,rank) schedule(dynamic)
        for (Pshell=0; Pshell < ribasis_->nshell(); ++Pshell) {
          int numPshell = ribasis_->shell(Pshell)->nfunction();
  
          #ifdef _OPENMP
            rank = omp_get_thread_num();
          #endif
  
          if (sqrt(Schwartz[MUNU]*DFSchwartz[Pshell])>schwarz_) {
            eri[rank]->compute_shell(Pshell, 0, MU, NU);
 
            if (MU != NU) {
              for (int P=0, index=0; P < numPshell; ++P) {
                int oP = ribasis_->shell(Pshell)->function_index() + P;
    
                for (int mu=0,munu=0; mu < nummu; ++mu) {
                  int omu = basisset_->shell(MU)->function_index() + mu;
    
                  for (int nu=0; nu < numnu; ++nu, ++index, ++munu) {
                    int onu = basisset_->shell(NU)->function_index() + nu;
    
                    AO_RI[munu+munu_offset][oP] = buffer[rank][index];
                  }
                }
              }
            }
            else {
              for (int P=0; P < numPshell; ++P) {
                int oP = ribasis_->shell(Pshell)->function_index() + P;
  
                for (int mu=0,munu=0; mu < nummu; ++mu) {
                  int omu = basisset_->shell(MU)->function_index() + mu;
    
                  for (int nu=0; nu <= mu; ++nu, ++munu) {
                    int onu = basisset_->shell(NU)->function_index() + nu;
                    int index = P*nummu*nummu + mu*nummu + nu;
    
                    AO_RI[munu+munu_offset][oP] = buffer[rank][index];
                  }
                }
              }
            }
          }
        }
}
        if (MU != NU) {
          munu_offset += nummu*numnu;
        }
        else {
          munu_offset += nummu*(nummu+1)/2;
        }
      }
      if (PQ_stop[curr_block] == MUNU) {

        C_DGEMM('N','T',ndf_,block_length[curr_block],ndf_,1.0,J_mhalf[0],
          ndf_,&(AO_RI[0][0]),ndf_,0.0,&(J_AO_RI[0][0]),max_size);

        for (int P=0; P < ndf_; ++P) {
          next_DF_AO = psio_get_address(PSIO_ZERO,
            sizeof(double)*P*nsotri_screened+sizeof(double)*offset);
          psio_->write(PSIF_SAPT_TEMP,"AO RI Integrals",(char *) 
            &(J_AO_RI[P][0]),sizeof(double)*block_length[curr_block],
            next_DF_AO,&next_DF_AO);
        }

        memset(&(AO_RI[0][0]),'\0',sizeof(double)*max_size*ndf_);
        memset(&(J_AO_RI[0][0]),'\0',sizeof(double)*max_size*ndf_);

        offset += block_length[curr_block];
        munu_offset = 0;
        curr_block++;
      }

  }}

  C_DGEMM('N','T',ndf_,block_length[curr_block],ndf_,1.0,J_mhalf[0],
    ndf_,&(AO_RI[0][0]),ndf_,0.0,&(J_AO_RI[0][0]),max_size);

  for (int P=0; P < ndf_; ++P) {
    next_DF_AO = psio_get_address(PSIO_ZERO,
      sizeof(double)*P*nsotri_screened+sizeof(double)*offset);
    psio_->write(PSIF_SAPT_TEMP,"AO RI Integrals",(char *)
      &(J_AO_RI[P][0]),sizeof(double)*block_length[curr_block],
      next_DF_AO,&next_DF_AO);
  }

  free_block(J_mhalf);
  free_block(AO_RI);
  free_block(J_AO_RI);

  avail_mem = mem_;
  long int indices = nsotri_screened + noccA_*noccA_ + noccA_*nvirA_ 
    + nvirA_*(nvirA_+1)/2 + noccB_*noccB_ + noccB_*nvirB_ 
    + nvirB_*(nvirB_+1)/2 + noccB_*noccB_ + noccA_*nvirB_ + noccB_*nvirA_;
  mem_tot = (long int) ndf_*indices;
  if (indices > avail_mem)
    throw PsiException("Not enough memory", __FILE__,__LINE__);
  max_size = avail_mem / (indices);
  if (max_size > ndf_)
    max_size = ndf_;

  int Pblocks = ndf_/max_size;
  int gimp = ndf_%max_size;

  if (gimp) Pblocks++;

  int Plength = max_size;

  double **B_p_munu = block_matrix(Plength,nsotri_screened);
  double **B_p_AA = block_matrix(Plength,noccA_*noccA_);
  double **B_p_AR = block_matrix(Plength,noccA_*nvirA_);
  double **B_p_RR = block_matrix(Plength,nvirA_*(nvirA_+1)/2);
  double **B_p_BB = block_matrix(Plength,noccB_*noccB_);
  double **B_p_BS = block_matrix(Plength,noccB_*nvirB_);
  double **B_p_SS = block_matrix(Plength,nvirB_*(nvirB_+1)/2);
  double **B_p_AB = block_matrix(Plength,noccA_*noccB_);
  double **B_p_AS = block_matrix(Plength,noccA_*nvirB_);
  double **B_p_RB = block_matrix(Plength,noccB_*nvirA_);

  double **munu_temp = block_matrix(nthreads,nso_*nso_);
  double **Inu_temp = block_matrix(nthreads,nmo_*nso_);
  double **IJ_temp = block_matrix(nthreads,nmo_*nmo_);
 
  next_DF_AO = PSIO_ZERO;
  psio_address next_DF_AA = PSIO_ZERO;
  psio_address next_DF_AR = PSIO_ZERO;
  psio_address next_DF_RR = PSIO_ZERO;
  psio_address next_DF_BB = PSIO_ZERO;
  psio_address next_DF_BS = PSIO_ZERO;
  psio_address next_DF_SS = PSIO_ZERO;
  psio_address next_DF_AB = PSIO_ZERO;
  psio_address next_DF_AS = PSIO_ZERO;
  psio_address next_DF_RB = PSIO_ZERO;

  zero_disk(PSIF_SAPT_AA_DF_INTS,"AA RI Integrals",ndf_,noccA_*noccA_);
  zero_disk(PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",ndf_,noccA_*nvirA_);
  zero_disk(PSIF_SAPT_AA_DF_INTS,"RR RI Integrals",ndf_,nvirA_*(nvirA_+1)/2);
  zero_disk(PSIF_SAPT_BB_DF_INTS,"BB RI Integrals",ndf_,noccB_*noccB_);
  zero_disk(PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",ndf_,noccB_*nvirB_);
  zero_disk(PSIF_SAPT_BB_DF_INTS,"SS RI Integrals",ndf_,nvirB_*(nvirB_+1)/2);
  zero_disk(PSIF_SAPT_AB_DF_INTS,"AB RI Integrals",ndf_,noccA_*noccB_);
  zero_disk(PSIF_SAPT_AB_DF_INTS,"AS RI Integrals",ndf_,noccA_*nvirB_);
  zero_disk(PSIF_SAPT_AB_DF_INTS,"RB RI Integrals",ndf_,nvirA_*noccB_);

  int Prel;

  for (int Pbl=0; Pbl<Pblocks; Pbl++) {

    int length = max_size;
    if (gimp && Pbl == Pblocks-1) length = gimp;

    psio_->read(PSIF_SAPT_TEMP,"AO RI Integrals",(char *) &(B_p_munu[0][0]),
      sizeof(double)*length*nsotri_screened,next_DF_AO,&next_DF_AO);

#pragma omp parallel
{
    #pragma omp for private(Prel,rank) schedule(dynamic)
    for (Prel=0; Prel<length; Prel++) {

      #ifdef _OPENMP
        rank = omp_get_thread_num();
      #endif

      memset(&(munu_temp[rank][0]),'\0',sizeof(double)*nso_*nso_);

      int PQoff = 0;
      for(int MU=0,MUNU=0;MU<basisset_->nshell();MU++) {
        int nummu = basisset_->shell(MU)->nfunction();
        for(int NU=0;NU<=MU;NU++,MUNU++) {
          int numnu = basisset_->shell(NU)->nfunction();
          if (sqrt(Schwartz[MUNU]*maxSchwartz)>schwarz_) { 

            if (MU != NU) {
              for (int mu=0,munu=0; mu < nummu; ++mu) {
                int omu = basisset_->shell(MU)->function_index() + mu;
  
                for (int nu=0; nu < numnu; ++nu, ++munu) {
                  int onu = basisset_->shell(NU)->function_index() + nu;
  
                  munu_temp[rank][omu*nso_+onu] = 
                    B_p_munu[Prel][munu+PQoff];
                  munu_temp[rank][onu*nso_+omu] = 
                    B_p_munu[Prel][munu+PQoff];
                }
              }
              PQoff += nummu*numnu;
            }
            else {
              for (int mu=0,munu=0; mu < nummu; ++mu) {
                int omu = basisset_->shell(MU)->function_index() + mu;
  
                for (int nu=0; nu <= mu; ++nu, ++munu) {
                  int onu = basisset_->shell(NU)->function_index() + nu;
  
                  munu_temp[rank][omu*nso_+onu] = 
                    B_p_munu[Prel][munu+PQoff];
                  munu_temp[rank][onu*nso_+omu] = 
                    B_p_munu[Prel][munu+PQoff];
                }
              }
              PQoff += nummu*(nummu+1)/2;
            }
          }
        }
      }
  
      C_DGEMM('T', 'N', nmoA_, nso_, nso_, 1.0, &(CA_[0][0]), nmoA_,
        munu_temp[rank], nso_, 0.0, Inu_temp[rank], nso_);
      C_DGEMM('N', 'N', nmoA_, nmoA_, nso_, 1.0, Inu_temp[rank], nso_,
        &(CA_[0][0]), nmoA_, 0.0, IJ_temp[rank], nmoA_);
  
      for (int a=0; a<noccA_; a++) {
        C_DCOPY(noccA_,&(IJ_temp[rank][a*nmoA_]),1,&(B_p_AA[Prel][a*noccA_]),1);
        C_DCOPY(nvirA_,&(IJ_temp[rank][a*nmoA_+noccA_]),1,
          &(B_p_AR[Prel][a*nvirA_]),1);
      }
      for (int r=0; r<nvirA_; r++) {
        C_DCOPY(r+1,&(IJ_temp[rank][(r+noccA_)*nmoA_+noccA_]),1,
          &(B_p_RR[Prel][r*(r+1)/2]),1);
      }
  
      C_DGEMM('N', 'N', nmoA_, nmoB_, nso_, 1.0, Inu_temp[rank], nso_,
        &(CB_[0][0]), nmoB_, 0.0, IJ_temp[rank], nmoB_);
  
      for (int a=0; a<noccA_; a++) {
        C_DCOPY(noccB_,&(IJ_temp[rank][a*nmoB_]),1,&(B_p_AB[Prel][a*noccB_]),1);
        C_DCOPY(nvirB_,&(IJ_temp[rank][a*nmoB_+noccB_]),1,
          &(B_p_AS[Prel][a*nvirB_]),1);
      }
      for (int r=0; r<nvirA_; r++) {
        C_DCOPY(noccB_,&(IJ_temp[rank][(r+noccA_)*nmoB_]),1,
          &(B_p_RB[Prel][r*noccB_]),1);
      } 
  
      C_DGEMM('T', 'N', nmoB_, nso_, nso_, 1.0, &(CB_[0][0]), nmoB_,
        munu_temp[rank], nso_, 0.0, Inu_temp[rank], nso_);
      C_DGEMM('N', 'N', nmoB_, nmoB_, nso_, 1.0, Inu_temp[rank], nso_,
        &(CB_[0][0]), nmoB_, 0.0, IJ_temp[rank], nmoB_);
  
      for (int b=0; b<noccB_; b++) {
        C_DCOPY(noccB_,&(IJ_temp[rank][b*nmoB_]),1,&(B_p_BB[Prel][b*noccB_]),1);
        C_DCOPY(nvirB_,&(IJ_temp[rank][b*nmoB_+noccB_]),1,
          &(B_p_BS[Prel][b*nvirB_]),1);
      }
      for (int s=0; s<nvirB_; s++) {
        C_DCOPY(s+1,&(IJ_temp[rank][(s+noccB_)*nmoB_+noccB_]),1,
          &(B_p_SS[Prel][s*(s+1)/2]),1);
      } 
 
    }
}
    psio_->write(PSIF_SAPT_AA_DF_INTS,"AA RI Integrals",(char *)
      &(B_p_AA[0][0]),sizeof(double)*length*noccA_*noccA_,
      next_DF_AA,&next_DF_AA);
    psio_->write(PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",(char *)
      &(B_p_AR[0][0]),sizeof(double)*length*noccA_*nvirA_,
      next_DF_AR,&next_DF_AR);
    psio_->write(PSIF_SAPT_AA_DF_INTS,"RR RI Integrals",(char *)
      &(B_p_RR[0][0]),sizeof(double)*length*(nvirA_*(nvirA_+1)/2),
      next_DF_RR,&next_DF_RR);

    psio_->write(PSIF_SAPT_BB_DF_INTS,"BB RI Integrals",(char *)
      &(B_p_BB[0][0]),sizeof(double)*length*noccB_*noccB_,
      next_DF_BB,&next_DF_BB);
    psio_->write(PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",(char *)
      &(B_p_BS[0][0]),sizeof(double)*length*noccB_*nvirB_,
      next_DF_BS,&next_DF_BS);
    psio_->write(PSIF_SAPT_BB_DF_INTS,"SS RI Integrals",(char *)
      &(B_p_SS[0][0]),sizeof(double)*length*(nvirB_*(nvirB_+1)/2),
      next_DF_SS,&next_DF_SS);

    psio_->write(PSIF_SAPT_AB_DF_INTS,"AB RI Integrals",(char *)
      &(B_p_AB[0][0]),sizeof(double)*length*noccA_*noccB_,
      next_DF_AB,&next_DF_AB);
    psio_->write(PSIF_SAPT_AB_DF_INTS,"AS RI Integrals",(char *)
      &(B_p_AS[0][0]),sizeof(double)*length*noccA_*nvirB_,
      next_DF_AS,&next_DF_AS);
    psio_->write(PSIF_SAPT_AB_DF_INTS,"RB RI Integrals",(char *)
      &(B_p_RB[0][0]),sizeof(double)*length*nvirA_*noccB_,
      next_DF_RB,&next_DF_RB);

  }

  free_block(B_p_munu); 
  free_block(B_p_AA);
  free_block(B_p_AR);
  free_block(B_p_RR);
  free_block(B_p_BB);
  free_block(B_p_BS);
  free_block(B_p_SS);
  free_block(B_p_AB);
  free_block(B_p_AS);
  free_block(B_p_RB);
  free_block(munu_temp);
  free_block(Inu_temp);
  free_block(IJ_temp);
  free(Schwartz);
  free(DFSchwartz);
  free(PQ_start);
  free(PQ_stop);
  free(block_length);

  for(int i = 0; i < nthreads; ++i)
    eri[i].reset();

  psio_->close(PSIF_SAPT_TEMP,0);
}

void SAPT0::df_integrals_aio()
{
  boost::shared_ptr<AIOHandler> aio(new AIOHandler(psio_));

  psio_->open(PSIF_SAPT_TEMP,PSIO_OPEN_NEW);

  // Get fitting metric
  boost::shared_ptr<FittingMetric> metric = boost::shared_ptr<FittingMetric>(
    new FittingMetric(ribasis_));
  metric->form_eig_inverse();
  double **J_temp = metric->get_metric()->pointer();
  double **J_mhalf = block_matrix(ndf_,ndf_);
  C_DCOPY(ndf_*ndf_,J_temp[0],1,J_mhalf[0],1);
  metric.reset();

  // Get Schwartz screening arrays
  double maxSchwartz = 0.0;
  double *Schwartz = init_array(basisset_->nshell()*(basisset_->nshell()+1)/2);

  boost::shared_ptr<IntegralFactory> ao_eri_factory = 
    boost::shared_ptr<IntegralFactory>(new IntegralFactory(basisset_, 
    basisset_, basisset_, basisset_));
  boost::shared_ptr<TwoBodyAOInt> ao_eri = boost::shared_ptr<TwoBodyAOInt>(
    ao_eri_factory->eri());
  const double *ao_buffer = ao_eri->buffer();

  for(int P=0,PQ=0;P<basisset_->nshell();P++) {
    int numw = basisset_->shell(P)->nfunction();
    for(int Q=0;Q<=P;Q++,PQ++) {
      int numx = basisset_->shell(Q)->nfunction();
      double tei, max=0.0;

      ao_eri->compute_shell(P, Q, P, Q);

      for(int w=0;w<numw;w++) {
        for(int x=0;x<numx;x++) {
          int index = ( ( (w*numx + x) * numw + w) * numx + x);
          tei = ao_buffer[index];
          if(fabs(tei) > max) max = fabs(tei);
        }
      }
      Schwartz[PQ] = max;
      if (max > maxSchwartz) maxSchwartz = max;
    }
  }

  ao_eri.reset();
  ao_eri_factory.reset();

  double *DFSchwartz = init_array(ribasis_->nshell());

  boost::shared_ptr<IntegralFactory> df_eri_factory = 
    boost::shared_ptr<IntegralFactory>(new IntegralFactory(ribasis_, zero_, 
    ribasis_, zero_));
  boost::shared_ptr<TwoBodyAOInt> df_eri = boost::shared_ptr<TwoBodyAOInt>(
    df_eri_factory->eri());
  const double *df_buffer = df_eri->buffer();

  for(int P=0;P<ribasis_->nshell();P++) {
    int numw = ribasis_->shell(P)->nfunction();
    double tei, max=0.0;

    df_eri->compute_shell(P, 0, P, 0);

    for(int w=0;w<numw;w++) {
      tei = df_buffer[w];
      if(fabs(tei) > max) max = fabs(tei);
    }
    DFSchwartz[P] = max;
  }

  df_eri.reset();
  df_eri_factory.reset();

  long int nsotri = nso_*(nso_+1)/2;
  long int nsotri_screened = 0;

  for(int P=0,PQ=0;P<basisset_->nshell();P++) {
    int numw = basisset_->shell(P)->nfunction();
    for(int Q=0;Q<=P;Q++,PQ++) {
      int numx = basisset_->shell(Q)->nfunction();
      int numPQ = numw*numx;
      if (P == Q) numPQ = numw*(numw+1)/2;
      if (sqrt(Schwartz[PQ]*maxSchwartz)>schwarz_)
        nsotri_screened += numPQ;
  }}

  long int avail_mem = mem_ - (long int) ndf_*ndf_; 
  long int mem_tot = (long int) 4*ndf_*nsotri_screened + (long int) ndf_*ndf_;
  if (avail_mem < (long int) 4*ndf_)
    throw PsiException("Not enough memory", __FILE__,__LINE__);
  long int max_size = avail_mem / ((long int) 4*ndf_);
  if (max_size > nsotri_screened)
    max_size = nsotri_screened;

  if (debug_) {
    fprintf(outfile,"Requires storage of %ld doubles\n",mem_tot);
    fprintf(outfile,"Max nso x nso block is %ld\n\n",max_size);
    fflush(outfile);
  }

  int size = 0;
  int num_blocks = 1;

  for(int P=0,PQ=0;P<basisset_->nshell();P++) {
    int numw = basisset_->shell(P)->nfunction();
    for(int Q=0;Q<=P;Q++,PQ++) {
      int numx = basisset_->shell(Q)->nfunction();
      int numPQ = numw*numx;
      if (P == Q) numPQ = numw*(numw+1)/2;

      if (sqrt(Schwartz[PQ]*maxSchwartz)>schwarz_) {
        size += numPQ;
        if (max_size < size) {
          if (debug_)
            fprintf(outfile,"Block %d : %d\n",num_blocks,size-numPQ);
          num_blocks++;
          size = numPQ;
        }
      }
  }} 

  if (debug_)
    fprintf(outfile,"Block %d : %d\n\n",num_blocks,size);

  int *PQ_start = init_int_array(num_blocks);
  int *PQ_stop = init_int_array(num_blocks);
  int *block_length =  init_int_array(num_blocks);

  int block_num = 0;
  size = 0;

  PQ_start[0] = 0;

  for(int P=0,PQ=0;P<basisset_->nshell();P++) {
    int numw = basisset_->shell(P)->nfunction();
    for(int Q=0;Q<=P;Q++,PQ++) {
      int numx = basisset_->shell(Q)->nfunction();
      int numPQ = numw*numx;
      if (P == Q) numPQ = numw*(numw+1)/2;
      if (sqrt(Schwartz[PQ]*maxSchwartz)>schwarz_) {
        size += numPQ;
        if (max_size < size) {
          PQ_stop[block_num] = PQ - 1;
          block_length[block_num] = size-numPQ;
          block_num++;
          PQ_start[block_num] = PQ;
          size = numPQ;
        }
      }
  }} 

  PQ_stop[num_blocks-1] = basisset_->nshell()*(basisset_->nshell()+1)/2;
  block_length[num_blocks-1] = size;

  int max_func_per_shell = basisset_->max_function_per_shell();

  if (max_func_per_shell*max_func_per_shell > max_size) 
    throw PsiException("Not enough memory", __FILE__,__LINE__);

  if (debug_) {
    for (int i=0; i<num_blocks; i++) 
      fprintf(outfile,"Block %2d : PQ %4d - %4d : %d\n",i,PQ_start[i],
        PQ_stop[i],block_length[i]);
    fprintf(outfile,"\n");
    fflush(outfile);
  }

  boost::shared_ptr<IntegralFactory> rifactory = 
    boost::shared_ptr<IntegralFactory>(new IntegralFactory(ribasis_, zero_, 
    basisset_, basisset_));

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

  zero_disk(PSIF_SAPT_TEMP,"AO RI Integrals",ndf_,nsotri_screened);

  psio_address next_DF_AO = PSIO_ZERO;

  double** AO_RI[2]; 
  double** J_AO_RI[2]; 

  AO_RI[0] = block_matrix(max_size,ndf_);
  J_AO_RI[0] = block_matrix(ndf_,max_size);
  AO_RI[1] = block_matrix(max_size,ndf_);
  J_AO_RI[1] = block_matrix(ndf_,max_size);

  int munu_offset = 0;
  int curr_block = 0;
  int offset = 0;
  int Pshell;

  for(int MU=0,MUNU=0;MU<basisset_->nshell();MU++) {
    int nummu = basisset_->shell(MU)->nfunction();
    for(int NU=0;NU<=MU;NU++,MUNU++) {
      int numnu = basisset_->shell(NU)->nfunction();

      if (sqrt(Schwartz[MUNU]*maxSchwartz)>schwarz_ ) {

#pragma omp parallel
{
        #pragma omp for private(Pshell,rank) schedule(dynamic)
        for (Pshell=0; Pshell < ribasis_->nshell(); ++Pshell) {
          int numPshell = ribasis_->shell(Pshell)->nfunction();
  
          #ifdef _OPENMP
            rank = omp_get_thread_num();
          #endif
  
          if (sqrt(Schwartz[MUNU]*DFSchwartz[Pshell])>schwarz_) { 
            eri[rank]->compute_shell(Pshell, 0, MU, NU);
  
            if (MU != NU) {
              for (int P=0, index=0; P < numPshell; ++P) {
                int oP = ribasis_->shell(Pshell)->function_index() + P;
    
                for (int mu=0,munu=0; mu < nummu; ++mu) {
                  int omu = basisset_->shell(MU)->function_index() + mu;
    
                  for (int nu=0; nu < numnu; ++nu, ++index, ++munu) {
                    int onu = basisset_->shell(NU)->function_index() + nu;
    
                    AO_RI[curr_block%2][munu+munu_offset][oP] 
                      = buffer[rank][index];
                  }
                }
              }
            }
            else {
              for (int P=0; P < numPshell; ++P) {
                int oP = ribasis_->shell(Pshell)->function_index() + P;
  
                for (int mu=0,munu=0; mu < nummu; ++mu) {
                  int omu = basisset_->shell(MU)->function_index() + mu;
    
                  for (int nu=0; nu <= mu; ++nu, ++munu) {
                    int onu = basisset_->shell(NU)->function_index() + nu;
                    int index = P*nummu*nummu + mu*nummu + nu;
    
                    AO_RI[curr_block%2][munu+munu_offset][oP] 
                      = buffer[rank][index];
                  }
                }
              }
            }
          }
        }
}
        if (MU != NU) {
          munu_offset += nummu*numnu;
        }
        else {
          munu_offset += nummu*(nummu+1)/2;
        }
      }

      if (PQ_stop[curr_block] == MUNU) {

        C_DGEMM('N','T',ndf_,block_length[curr_block],ndf_,1.0,J_mhalf[0],
          ndf_,&(AO_RI[curr_block%2][0][0]),ndf_,0.0,
          &(J_AO_RI[curr_block%2][0][0]),max_size);

        if (curr_block > 0)
          aio->synchronize();

        next_DF_AO = psio_get_address(PSIO_ZERO,sizeof(double)*offset);
        aio->write_discont(PSIF_SAPT_TEMP,"AO RI Integrals",
          J_AO_RI[curr_block%2],ndf_,block_length[curr_block],
          nsotri_screened-block_length[curr_block],next_DF_AO);

        offset += block_length[curr_block];
        munu_offset = 0;
        curr_block++;

        memset(&(AO_RI[curr_block%2][0][0]),'\0',sizeof(double)*max_size
          *ndf_);
        memset(&(J_AO_RI[curr_block%2][0][0]),'\0',sizeof(double)*max_size
          *ndf_);
      }

  }}

  C_DGEMM('N','T',ndf_,block_length[curr_block],ndf_,1.0,J_mhalf[0],
    ndf_,&(AO_RI[curr_block%2][0][0]),ndf_,0.0,
    &(J_AO_RI[curr_block%2][0][0]),max_size);

  if (curr_block > 0)
    aio->synchronize();

  next_DF_AO = psio_get_address(PSIO_ZERO,sizeof(double)*offset);
  for (int P=0; P < ndf_; ++P) {
    psio_->write(PSIF_SAPT_TEMP,"AO RI Integrals",(char *)
      &(J_AO_RI[curr_block%2][P][0]),sizeof(double)*block_length[curr_block],
      next_DF_AO,&next_DF_AO);
    next_DF_AO = psio_get_address(next_DF_AO,
      sizeof(double)*(nsotri_screened-block_length[curr_block]));
  }

  free_block(J_mhalf);
  free_block(AO_RI[0]);
  free_block(J_AO_RI[0]);
  free_block(AO_RI[1]);
  free_block(J_AO_RI[1]);

  avail_mem = mem_;
  long int indices = nsotri_screened + noccA_*noccA_ + noccA_*nvirA_ 
    + nvirA_*(nvirA_+1)/2 + noccB_*noccB_ + noccB_*nvirB_ 
    + nvirB_*(nvirB_+1)/2 + noccB_*noccB_ + noccA_*nvirB_ + noccB_*nvirA_;
  mem_tot = (long int) ndf_*indices*2;
  if (indices*2 > avail_mem)
    throw PsiException("Not enough memory", __FILE__,__LINE__);
  max_size = avail_mem / (indices*2);
  if (max_size > ndf_)
    max_size = ndf_;

  int Pblocks = ndf_/max_size;
  int gimp = ndf_%max_size;

  if (gimp) Pblocks++;

  int Plength = max_size;

  double **B_p_munu[2];
  double **B_p_AA[2];
  double **B_p_AR[2];
  double **B_p_RR[2];
  double **B_p_BB[2];
  double **B_p_BS[2];
  double **B_p_SS[2];
  double **B_p_AB[2];
  double **B_p_AS[2];
  double **B_p_RB[2];

  B_p_munu[0] = block_matrix(Plength,nsotri_screened);
  B_p_AA[0] = block_matrix(Plength,noccA_*noccA_);
  B_p_AR[0] = block_matrix(Plength,noccA_*nvirA_);
  B_p_RR[0] = block_matrix(Plength,nvirA_*(nvirA_+1)/2);
  B_p_BB[0] = block_matrix(Plength,noccB_*noccB_);
  B_p_BS[0] = block_matrix(Plength,noccB_*nvirB_);
  B_p_SS[0] = block_matrix(Plength,nvirB_*(nvirB_+1)/2);
  B_p_AB[0] = block_matrix(Plength,noccA_*noccB_);
  B_p_AS[0] = block_matrix(Plength,noccA_*nvirB_);
  B_p_RB[0] = block_matrix(Plength,noccB_*nvirA_);

  B_p_munu[1] = block_matrix(Plength,nsotri_screened);
  B_p_AA[1] = block_matrix(Plength,noccA_*noccA_);
  B_p_AR[1] = block_matrix(Plength,noccA_*nvirA_);
  B_p_RR[1] = block_matrix(Plength,nvirA_*(nvirA_+1)/2);
  B_p_BB[1] = block_matrix(Plength,noccB_*noccB_);
  B_p_BS[1] = block_matrix(Plength,noccB_*nvirB_);
  B_p_SS[1] = block_matrix(Plength,nvirB_*(nvirB_+1)/2);
  B_p_AB[1] = block_matrix(Plength,noccA_*noccB_);
  B_p_AS[1] = block_matrix(Plength,noccA_*nvirB_);
  B_p_RB[1] = block_matrix(Plength,noccB_*nvirA_);

  double **munu_temp = block_matrix(nthreads,nso_*nso_);
  double **Inu_temp = block_matrix(nthreads,nmo_*nso_);
  double **IJ_temp = block_matrix(nthreads,nmo_*nmo_);
 
  next_DF_AO = PSIO_ZERO;
  psio_address next_DF_AA = PSIO_ZERO;
  psio_address next_DF_AR = PSIO_ZERO;
  psio_address next_DF_RR = PSIO_ZERO;
  psio_address next_DF_BB = PSIO_ZERO;
  psio_address next_DF_BS = PSIO_ZERO;
  psio_address next_DF_SS = PSIO_ZERO;
  psio_address next_DF_AB = PSIO_ZERO;
  psio_address next_DF_AS = PSIO_ZERO;
  psio_address next_DF_RB = PSIO_ZERO;

  zero_disk(PSIF_SAPT_AA_DF_INTS,"AA RI Integrals",ndf_,noccA_*noccA_);
  zero_disk(PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",ndf_,noccA_*nvirA_);
  zero_disk(PSIF_SAPT_AA_DF_INTS,"RR RI Integrals",ndf_,nvirA_*(nvirA_+1)/2);
  zero_disk(PSIF_SAPT_BB_DF_INTS,"BB RI Integrals",ndf_,noccB_*noccB_);
  zero_disk(PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",ndf_,noccB_*nvirB_);
  zero_disk(PSIF_SAPT_BB_DF_INTS,"SS RI Integrals",ndf_,nvirB_*(nvirB_+1)/2);
  zero_disk(PSIF_SAPT_AB_DF_INTS,"AB RI Integrals",ndf_,noccA_*noccB_);
  zero_disk(PSIF_SAPT_AB_DF_INTS,"AS RI Integrals",ndf_,noccA_*nvirB_);
  zero_disk(PSIF_SAPT_AB_DF_INTS,"RB RI Integrals",ndf_,nvirA_*noccB_);

  psio_->read(PSIF_SAPT_TEMP,"AO RI Integrals",(char *) &(B_p_munu[0][0][0]),
    sizeof(double)*max_size*nsotri_screened,next_DF_AO,&next_DF_AO);

  int Prel;

  for (int Pbl=0; Pbl<Pblocks; Pbl++) {

    int length = max_size;
    if (gimp && Pbl == Pblocks-1) length = gimp;

      if (Pbl < Pblocks-1) {
        int read_length = max_size;
        if (gimp && Pbl == Pblocks-2) read_length = gimp;
        aio->read(PSIF_SAPT_TEMP,"AO RI Integrals",(char *) 
          &(B_p_munu[(Pbl+1)%2][0][0]),sizeof(double)*read_length
          *nsotri_screened,next_DF_AO,&next_DF_AO);
      }

#pragma omp parallel
{
    #pragma omp for private(Prel,rank) schedule(dynamic)
    for (Prel=0; Prel<length; Prel++) {

      #ifdef _OPENMP
        rank = omp_get_thread_num();
      #endif

      memset(&(munu_temp[rank][0]),'\0',sizeof(double)*nso_*nso_);

      int PQoff = 0;
      for(int MU=0,MUNU=0;MU<basisset_->nshell();MU++) {
        int nummu = basisset_->shell(MU)->nfunction();
        for(int NU=0;NU<=MU;NU++,MUNU++) {
          int numnu = basisset_->shell(NU)->nfunction();
          if (sqrt(Schwartz[MUNU]*maxSchwartz)>schwarz_) {
 
            if (MU != NU) {
              for (int mu=0,munu=0; mu < nummu; ++mu) {
                int omu = basisset_->shell(MU)->function_index() + mu;
  
                for (int nu=0; nu < numnu; ++nu, ++munu) {
                  int onu = basisset_->shell(NU)->function_index() + nu;
  
                  munu_temp[rank][omu*nso_+onu] = 
                    B_p_munu[Pbl%2][Prel][munu+PQoff];
                  munu_temp[rank][onu*nso_+omu] = 
                    B_p_munu[Pbl%2][Prel][munu+PQoff];
                }
              }
              PQoff += nummu*numnu;
            }
            else {
              for (int mu=0,munu=0; mu < nummu; ++mu) {
                int omu = basisset_->shell(MU)->function_index() + mu;
  
                for (int nu=0; nu <= mu; ++nu, ++munu) {
                  int onu = basisset_->shell(NU)->function_index() + nu;
  
                  munu_temp[rank][omu*nso_+onu] = 
                    B_p_munu[Pbl%2][Prel][munu+PQoff];
                  munu_temp[rank][onu*nso_+omu] = 
                    B_p_munu[Pbl%2][Prel][munu+PQoff];
                }
              }
              PQoff += nummu*(nummu+1)/2;
            }
          }
        }
      }
  
      C_DGEMM('T', 'N', nmoA_, nso_, nso_, 1.0, &(CA_[0][0]), nmoA_,
        munu_temp[rank], nso_, 0.0, Inu_temp[rank], nso_);
      C_DGEMM('N', 'N', nmoA_, nmoA_, nso_, 1.0, Inu_temp[rank], nso_,
        &(CA_[0][0]), nmoA_, 0.0, IJ_temp[rank], nmoA_);
  
      for (int a=0; a<noccA_; a++) {
        C_DCOPY(noccA_,&(IJ_temp[rank][a*nmoA_]),1,
          &(B_p_AA[Pbl%2][Prel][a*noccA_]),1);
        C_DCOPY(nvirA_,&(IJ_temp[rank][a*nmoA_+noccA_]),1,
          &(B_p_AR[Pbl%2][Prel][a*nvirA_]),1);
      }
      for (int r=0; r<nvirA_; r++) {
        C_DCOPY(r+1,&(IJ_temp[rank][(r+noccA_)*nmoA_+noccA_]),1,
          &(B_p_RR[Pbl%2][Prel][r*(r+1)/2]),1);
      }
  
      C_DGEMM('N', 'N', nmoA_, nmoB_, nso_, 1.0, Inu_temp[rank], nso_,
        &(CB_[0][0]), nmoB_, 0.0, IJ_temp[rank], nmoB_);
  
      for (int a=0; a<noccA_; a++) {
        C_DCOPY(noccB_,&(IJ_temp[rank][a*nmoB_]),1,
          &(B_p_AB[Pbl%2][Prel][a*noccB_]),1);
        C_DCOPY(nvirB_,&(IJ_temp[rank][a*nmoB_+noccB_]),1,
          &(B_p_AS[Pbl%2][Prel][a*nvirB_]),1);
      }
      for (int r=0; r<nvirA_; r++) {
        C_DCOPY(noccB_,&(IJ_temp[rank][(r+noccA_)*nmoB_]),1,
          &(B_p_RB[Pbl%2][Prel][r*noccB_]),1);
      } 
  
      C_DGEMM('T', 'N', nmoB_, nso_, nso_, 1.0, &(CB_[0][0]), nmoB_,
        munu_temp[rank], nso_, 0.0, Inu_temp[rank], nso_);
      C_DGEMM('N', 'N', nmoB_, nmoB_, nso_, 1.0, Inu_temp[rank], nso_,
        &(CB_[0][0]), nmoB_, 0.0, IJ_temp[rank], nmoB_);
  
      for (int b=0; b<noccB_; b++) {
        C_DCOPY(noccB_,&(IJ_temp[rank][b*nmoB_]),1,
          &(B_p_BB[Pbl%2][Prel][b*noccB_]),1);
        C_DCOPY(nvirB_,&(IJ_temp[rank][b*nmoB_+noccB_]),1,
          &(B_p_BS[Pbl%2][Prel][b*nvirB_]),1);
      }
      for (int s=0; s<nvirB_; s++) {
        C_DCOPY(s+1,&(IJ_temp[rank][(s+noccB_)*nmoB_+noccB_]),1,
          &(B_p_SS[Pbl%2][Prel][s*(s+1)/2]),1);
      } 
 
    }
}
    if (Pblocks > 1)
      aio->synchronize();

    aio->write(PSIF_SAPT_AA_DF_INTS,"AA RI Integrals",(char *)
      &(B_p_AA[Pbl%2][0][0]),sizeof(double)*length*noccA_*noccA_,
      next_DF_AA,&next_DF_AA);
    aio->write(PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",(char *)
      &(B_p_AR[Pbl%2][0][0]),sizeof(double)*length*noccA_*nvirA_,
      next_DF_AR,&next_DF_AR);
    aio->write(PSIF_SAPT_AA_DF_INTS,"RR RI Integrals",(char *)
      &(B_p_RR[Pbl%2][0][0]),sizeof(double)*length*(nvirA_*(nvirA_+1)/2),
      next_DF_RR,&next_DF_RR);

    aio->write(PSIF_SAPT_BB_DF_INTS,"BB RI Integrals",(char *)
      &(B_p_BB[Pbl%2][0][0]),sizeof(double)*length*noccB_*noccB_,
      next_DF_BB,&next_DF_BB);
    aio->write(PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",(char *)
      &(B_p_BS[Pbl%2][0][0]),sizeof(double)*length*noccB_*nvirB_,
      next_DF_BS,&next_DF_BS);
    aio->write(PSIF_SAPT_BB_DF_INTS,"SS RI Integrals",(char *)
      &(B_p_SS[Pbl%2][0][0]),sizeof(double)*length*(nvirB_*(nvirB_+1)/2),
      next_DF_SS,&next_DF_SS);

    aio->write(PSIF_SAPT_AB_DF_INTS,"AB RI Integrals",(char *)
      &(B_p_AB[Pbl%2][0][0]),sizeof(double)*length*noccA_*noccB_,
      next_DF_AB,&next_DF_AB);
    aio->write(PSIF_SAPT_AB_DF_INTS,"AS RI Integrals",(char *)
      &(B_p_AS[Pbl%2][0][0]),sizeof(double)*length*noccA_*nvirB_,
      next_DF_AS,&next_DF_AS);
    aio->write(PSIF_SAPT_AB_DF_INTS,"RB RI Integrals",(char *)
      &(B_p_RB[Pbl%2][0][0]),sizeof(double)*length*nvirA_*noccB_,
      next_DF_RB,&next_DF_RB);

  }

  aio->synchronize();

  free_block(B_p_munu[0]); 
  free_block(B_p_AA[0]);
  free_block(B_p_AR[0]);
  free_block(B_p_RR[0]);
  free_block(B_p_BB[0]);
  free_block(B_p_BS[0]);
  free_block(B_p_SS[0]);
  free_block(B_p_AB[0]);
  free_block(B_p_AS[0]);
  free_block(B_p_RB[0]);
  free_block(B_p_munu[1]); 
  free_block(B_p_AA[1]);
  free_block(B_p_AR[1]);
  free_block(B_p_RR[1]);
  free_block(B_p_BB[1]);
  free_block(B_p_BS[1]);
  free_block(B_p_SS[1]);
  free_block(B_p_AB[1]);
  free_block(B_p_AS[1]);
  free_block(B_p_RB[1]);
  free_block(munu_temp);
  free_block(Inu_temp);
  free_block(IJ_temp);
  free(Schwartz);
  free(DFSchwartz);
  free(PQ_start);
  free(PQ_stop);
  free(block_length);

  for(int i = 0; i < nthreads; ++i)
    eri[i].reset();

  psio_->close(PSIF_SAPT_TEMP,0);
}

void SAPT0::w_integrals()
{
  diagAA_ = init_array(ndf_+3);
  SAPTDFInts B_p_AA = set_A_AA();
  Iterator AA_iter = get_iterator(mem_,&B_p_AA);

  for (int i=0,off=0; i<AA_iter.num_blocks; i++) {
    read_block(&AA_iter,&B_p_AA);

    for (int a=0; a<noccA_; a++){
      C_DAXPY(AA_iter.curr_size,1.0,&(B_p_AA.B_p_[0][a*noccA_+a]),
        noccA_*noccA_,&(diagAA_[off]),1);
    }

    off += AA_iter.curr_size;
  }

  B_p_AA.done();

  diagBB_ = init_array(ndf_+3);
  SAPTDFInts B_p_BB = set_B_BB();
  Iterator BB_iter = get_iterator(mem_,&B_p_BB);

  for (int i=0,off=0; i<BB_iter.num_blocks; i++) {
    read_block(&BB_iter,&B_p_BB);

    for (int b=0; b<noccB_; b++){
      C_DAXPY(BB_iter.curr_size,1.0,&(B_p_BB.B_p_[0][b*noccB_+b]),
        noccB_*noccB_,&(diagBB_[off]),1);
    }

    off += BB_iter.curr_size;
  }

  B_p_BB.done();

  wBAR_ = block_matrix(noccA_,nvirA_);
  SAPTDFInts B_p_AR = set_C_AR();
  Iterator AR_iter = get_iterator(mem_,&B_p_AR);

  for(int a=0; a<noccA_; a++){
    C_DAXPY(nvirA_,1.0,&(vBAA_[a][noccA_]),1,&(wBAR_[a][0]),1);
  }

  for (int i=0,off=0; i<AR_iter.num_blocks; i++) {
    read_block(&AR_iter,&B_p_AR);

    C_DGEMV('t',AR_iter.curr_size,noccA_*nvirA_,2.0,&(B_p_AR.B_p_[0][0]),
      noccA_*nvirA_,&(diagBB_[off]),1,1.0,&(wBAR_[0][0]),1);

    off += AR_iter.curr_size;
  }

  B_p_AR.done();

  wABS_ = block_matrix(noccB_,nvirB_);
  SAPTDFInts B_p_BS = set_C_BS();
  Iterator BS_iter = get_iterator(mem_,&B_p_BS);

  for(int b=0; b<noccB_; b++){
    C_DAXPY(nvirB_,1.0,&(vABB_[b][noccB_]),1,&(wABS_[b][0]),1);
  }
  
  for (int i=0,off=0; i<BS_iter.num_blocks; i++) {
    read_block(&BS_iter,&B_p_BS);
    
    C_DGEMV('t',BS_iter.curr_size,noccB_*nvirB_,2.0,&(B_p_BS.B_p_[0][0]),
      noccB_*nvirB_,&(diagAA_[off]),1,1.0,&(wABS_[0][0]),1);

    off += BS_iter.curr_size;
  }

  B_p_BS.done();
}

void SAPT0::oo_df_integrals()
{
  psio_->open(PSIF_SAPT_TEMP,PSIO_OPEN_NEW);
  
  // Get Schwartz screening arrays
  double maxSchwartz = 0.0;
  int nshelltri = basisset_->nshell()*(basisset_->nshell()+1)/2;
  double *Schwartz = init_array(basisset_->nshell()*(basisset_->nshell()+1)/2);
  
  boost::shared_ptr<IntegralFactory> ao_eri_factory = 
    boost::shared_ptr<IntegralFactory>(new IntegralFactory(basisset_, 
    basisset_, basisset_, basisset_));
  boost::shared_ptr<TwoBodyAOInt> ao_eri = boost::shared_ptr<TwoBodyAOInt>(
    ao_eri_factory->eri()); 
  const double *ao_buffer = ao_eri->buffer();
  
  for(int P=0,PQ=0;P<basisset_->nshell();P++) {
    int numw = basisset_->shell(P)->nfunction();
    for(int Q=0;Q<=P;Q++,PQ++) {
      int numx = basisset_->shell(Q)->nfunction();
      double tei, max=0.0;

      ao_eri->compute_shell(P, Q, P, Q);

      for(int w=0;w<numw;w++) { 
        for(int x=0;x<numx;x++) { 
          int index = ( ( (w*numx + x) * numw + w) * numx + x);
          tei = ao_buffer[index]; 
          if(fabs(tei) > max) max = fabs(tei);
        }
      }
      Schwartz[PQ] = max;
      if (max > maxSchwartz) maxSchwartz = max;
    }
  }
  
  ao_eri.reset();
  ao_eri_factory.reset();

  double *DFSchwartz = init_array(elstbasis_->nshell());
    
  boost::shared_ptr<IntegralFactory> df_eri_factory = 
    boost::shared_ptr<IntegralFactory>(new IntegralFactory(elstbasis_, zero_, 
    elstbasis_, zero_));
  boost::shared_ptr<TwoBodyAOInt> df_eri = boost::shared_ptr<TwoBodyAOInt>(
    df_eri_factory->eri());
  const double *df_buffer = df_eri->buffer();
  
  for(int P=0;P<elstbasis_->nshell();P++) {
    int numw = elstbasis_->shell(P)->nfunction();
    double tei, max=0.0;

    df_eri->compute_shell(P, 0, P, 0);

    for(int w=0;w<numw;w++) {
      tei = df_buffer[w];
      if(fabs(tei) > max) max = fabs(tei);
    }
    DFSchwartz[P] = max;
  }

  df_eri.reset();
  df_eri_factory.reset();

  int maxPshell = elstbasis_->max_function_per_shell();

  boost::shared_ptr<IntegralFactory> rifactory =
    boost::shared_ptr<IntegralFactory>(new IntegralFactory(elstbasis_, zero_,
    basisset_, basisset_));

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

  int *MUNUtoMU = init_int_array(nshelltri);
  int *MUNUtoNU = init_int_array(nshelltri);

  for(int MU=0, MUNU=0; MU < basisset_->nshell(); MU++) { 
    for(int NU=0; NU <= MU; NU++, MUNU++) {
      MUNUtoMU[MUNU] = MU;
      MUNUtoNU[MUNU] = NU;
  }}

  long int avail_mem = mem_ - (long int) maxPshell*(noccA_*noccA_ 
    + noccA_*noccB_ + noccB_*noccB_ + nso_*nso_);

  if (0 > avail_mem)
    throw PsiException("Not enough memory", __FILE__,__LINE__);

  double **B_p_AA = block_matrix(maxPshell,noccA_*noccA_);
  double **B_p_AB = block_matrix(maxPshell,noccA_*noccB_);
  double **B_p_BB = block_matrix(maxPshell,noccB_*noccB_);

  double **temp = block_matrix(maxPshell,nso_*nso_);
  double *tempA = init_array(noccA_*nso_);
  double *tempB = init_array(noccB_*nso_);

  psio_address next_DF_AA = PSIO_ZERO;
  psio_address next_DF_AB = PSIO_ZERO;
  psio_address next_DF_BB = PSIO_ZERO;

  zero_disk(PSIF_SAPT_TEMP,"AA RI Integrals",ndf_,noccA_*noccA_);
  zero_disk(PSIF_SAPT_TEMP,"AB RI Integrals",ndf_,noccA_*noccB_);
  zero_disk(PSIF_SAPT_TEMP,"BB RI Integrals",ndf_,noccB_*noccB_);

  for (int Pshell=0; Pshell<elstbasis_->nshell(); Pshell++) {
    int numPshell = elstbasis_->shell(Pshell)->nfunction();

#pragma omp parallel
{
    #pragma omp for private(rank) schedule(dynamic)
    for(int MUNU=0; MUNU < nshelltri; MUNU++) {
      #ifdef _OPENMP
        rank = omp_get_thread_num();
      #endif

      int MU = MUNUtoMU[MUNU];
      int NU = MUNUtoNU[MUNU];
      int nummu = basisset_->shell(MU)->nfunction();
      int numnu = basisset_->shell(NU)->nfunction();
      if (sqrt(Schwartz[MUNU]*maxSchwartz)>schwarz_ && 
        sqrt(Schwartz[MUNU]*DFSchwartz[Pshell])>schwarz_) {

        eri[rank]->compute_shell(Pshell, 0, MU, NU);
 
        for (int P=0, index=0; P < numPshell; ++P) {
          for (int mu=0; mu < nummu; ++mu) {
            int omu = basisset_->shell(MU)->function_index() + mu;
            for (int nu=0; nu < numnu; ++nu, ++index) {
              int onu = basisset_->shell(NU)->function_index() + nu;

              temp[P][omu*nso_+onu] = buffer[rank][index];
              temp[P][onu*nso_+omu] = buffer[rank][index];
            }
          }
        }
      }
    }
}
    for (int P=0, index=0; P < numPshell; ++P) {
      C_DGEMM('T', 'N', noccA_, nso_, nso_, 1.0, &(CA_[0][0]), nmoA_,
        temp[P], nso_, 0.0, tempA, nso_);
      C_DGEMM('N', 'N', noccA_, noccA_, nso_, 1.0, tempA, nso_,
        &(CA_[0][0]), nmoA_, 0.0, B_p_AA[P], noccA_);
      C_DGEMM('N', 'N', noccA_, noccB_, nso_, 1.0, tempA, nso_,
        &(CB_[0][0]), nmoB_, 0.0, B_p_AB[P], noccB_);
      C_DGEMM('T', 'N', noccB_, nso_, nso_, 1.0, &(CB_[0][0]), nmoB_,
        temp[P], nso_, 0.0, tempB, nso_);
      C_DGEMM('N', 'N', noccB_, noccB_, nso_, 1.0, tempB, nso_,
        &(CB_[0][0]), nmoB_, 0.0, B_p_BB[P], noccB_);
    }

    psio_->write(PSIF_SAPT_TEMP,"AA RI Integrals",(char *)
      &(B_p_AA[0][0]),sizeof(double)*numPshell*noccA_*noccA_,
      next_DF_AA,&next_DF_AA);
    psio_->write(PSIF_SAPT_TEMP,"AB RI Integrals",(char *)
      &(B_p_AB[0][0]),sizeof(double)*numPshell*noccA_*noccB_,
      next_DF_AB,&next_DF_AB);
    psio_->write(PSIF_SAPT_TEMP,"BB RI Integrals",(char *)
      &(B_p_BB[0][0]),sizeof(double)*numPshell*noccB_*noccB_,
      next_DF_BB,&next_DF_BB);
  }

  free(Schwartz);
  free(DFSchwartz);
  free_block(temp);
  free(tempA);
  free(tempB);
  free_block(B_p_AA);
  free_block(B_p_AB);
  free_block(B_p_BB);

  // Get fitting metric
  boost::shared_ptr<FittingMetric> metric = boost::shared_ptr<FittingMetric>(
    new FittingMetric(elstbasis_));
  metric->form_eig_inverse();
  double **J_temp = metric->get_metric()->pointer();
  double **J_mhalf = block_matrix(ndf_,ndf_);
  C_DCOPY(ndf_*ndf_,J_temp[0],1,J_mhalf[0],1);
  metric.reset();

  avail_mem = mem_ - (long int) ndf_*ndf_;

  if (2L*ndf_ > avail_mem)
    throw PsiException("Not enough memory", __FILE__,__LINE__);

  long int max_size = avail_mem/(2L*ndf_);
  if (max_size > noccA_*noccA_)
    max_size = noccA_*noccA_;

  B_p_AA = block_matrix(ndf_,max_size);
  double **B_q_AA = block_matrix(ndf_,max_size);

  next_DF_AA = PSIO_ZERO;
  psio_address next_DFJ_AA = PSIO_ZERO;

  int blocks = (noccA_*noccA_)/max_size;
  if ((noccA_*noccA_)%max_size)
    blocks++;

  zero_disk(PSIF_SAPT_AA_DF_INTS,"AA RI Integrals",ndf_,noccA_*noccA_);

  for (int n=0; n<blocks; n++) { 

    int start = n*max_size;
    int size = max_size;
    if (n == blocks-1)
      size = noccA_*noccA_ - start;

    for (int P=0; P<ndf_; P++) {
      next_DF_AA = psio_get_address(PSIO_ZERO,sizeof(double)*P*noccA_*noccA_
        + sizeof(double)*start);
      psio_->read(PSIF_SAPT_TEMP,"AA RI Integrals",(char *) &(B_p_AA[P][0]),
        sizeof(double)*size,next_DF_AA,&next_DF_AA);
    }
  
    C_DGEMM('N','N',ndf_,size,ndf_,1.0,J_mhalf[0],ndf_,B_p_AA[0],max_size,
      0.0,B_q_AA[0],max_size);
  
    for (int P=0; P<ndf_; P++) {
      next_DFJ_AA = psio_get_address(PSIO_ZERO,sizeof(double)*P*noccA_*noccA_
        + sizeof(double)*start);
      psio_->write(PSIF_SAPT_AA_DF_INTS,"AA RI Integrals",(char *) 
        &(B_q_AA[P][0]),sizeof(double)*size,next_DFJ_AA,
        &next_DFJ_AA);
    }
  }

  free_block(B_p_AA);
  free_block(B_q_AA);

  max_size = avail_mem/(2L*ndf_);
  if (max_size > noccA_*noccB_)
    max_size = noccA_*noccB_;

  B_p_AB = block_matrix(ndf_,max_size);
  double **B_q_AB = block_matrix(ndf_,max_size);

  next_DF_AB = PSIO_ZERO;
  psio_address next_DFJ_AB = PSIO_ZERO;

  blocks = (noccA_*noccB_)/max_size;
  if ((noccA_*noccB_)%max_size)
    blocks++;

  zero_disk(PSIF_SAPT_AB_DF_INTS,"AB RI Integrals",ndf_,noccA_*noccB_);

  for (int n=0; n<blocks; n++) {

    int start = n*max_size;
    int size = max_size;
    if (n == blocks-1)
      size = noccA_*noccB_ - start;

    for (int P=0; P<ndf_; P++) {
      next_DF_AB = psio_get_address(PSIO_ZERO,sizeof(double)*P*noccA_*noccB_
        + sizeof(double)*start);
      psio_->read(PSIF_SAPT_TEMP,"AB RI Integrals",(char *) &(B_p_AB[P][0]),
        sizeof(double)*size,next_DF_AB,&next_DF_AB);
    }

    C_DGEMM('N','N',ndf_,size,ndf_,1.0,J_mhalf[0],
      ndf_,B_p_AB[0],max_size,0.0,B_q_AB[0],max_size);

    for (int P=0; P<ndf_; P++) {
      next_DFJ_AB = psio_get_address(PSIO_ZERO,sizeof(double)*P*noccA_*noccB_
        + sizeof(double)*start);
      psio_->write(PSIF_SAPT_AB_DF_INTS,"AB RI Integrals",(char *)
        &(B_q_AB[P][0]),sizeof(double)*size,next_DFJ_AB,
        &next_DFJ_AB);
    }
  }

  free_block(B_p_AB);
  free_block(B_q_AB);

  max_size = avail_mem/(2L*ndf_);
  if (max_size > noccB_*noccB_)
    max_size = noccB_*noccB_;

  B_p_BB = block_matrix(ndf_,max_size);
  double **B_q_BB = block_matrix(ndf_,max_size);

  next_DF_BB = PSIO_ZERO;
  psio_address next_DFJ_BB = PSIO_ZERO;

  blocks = (noccB_*noccB_)/max_size;
  if ((noccB_*noccB_)%max_size)
    blocks++;

  zero_disk(PSIF_SAPT_BB_DF_INTS,"BB RI Integrals",ndf_,noccB_*noccB_);

  for (int n=0; n<blocks; n++) {

    int start = n*max_size;
    int size = max_size;
    if (n == blocks-1)
      size = noccB_*noccB_ - start;

    for (int P=0; P<ndf_; P++) {
      next_DF_BB = psio_get_address(PSIO_ZERO,sizeof(double)*P*noccB_*noccB_
        + sizeof(double)*start);
      psio_->read(PSIF_SAPT_TEMP,"BB RI Integrals",(char *) &(B_p_BB[P][0]),
        sizeof(double)*size,next_DF_BB,&next_DF_BB);
    }

    C_DGEMM('N','N',ndf_,size,ndf_,1.0,J_mhalf[0],
      ndf_,B_p_BB[0],max_size,0.0,B_q_BB[0],max_size);

    for (int P=0; P<ndf_; P++) {
      next_DFJ_BB = psio_get_address(PSIO_ZERO,sizeof(double)*P*noccB_*noccB_
        + sizeof(double)*start);
      psio_->write(PSIF_SAPT_BB_DF_INTS,"BB RI Integrals",(char *)
        &(B_q_BB[P][0]),sizeof(double)*size,next_DFJ_BB,
        &next_DFJ_BB);
    }
  }

  free_block(B_p_BB);
  free_block(B_q_BB);

  free_block(J_mhalf);
  free(MUNUtoMU);
  free(MUNUtoNU);

  psio_->close(PSIF_SAPT_TEMP,0);

  diagAA_ = init_array(ndf_+3);
  SAPTDFInts C_p_AA = set_A_AA();
  Iterator AA_iter = get_iterator(mem_,&C_p_AA);

  for (int i=0,off=0; i<AA_iter.num_blocks; i++) {
    read_block(&AA_iter,&C_p_AA);

    for (int a=0; a<noccA_; a++){
      C_DAXPY(AA_iter.curr_size,1.0,&(C_p_AA.B_p_[0][a*noccA_+a]),
        noccA_*noccA_,&(diagAA_[off]),1);
    }

    off += AA_iter.curr_size;
  }

  C_p_AA.done();

  diagBB_ = init_array(ndf_+3);
  SAPTDFInts C_p_BB = set_B_BB();
  Iterator BB_iter = get_iterator(mem_,&C_p_BB);

  for (int i=0,off=0; i<BB_iter.num_blocks; i++) {
    read_block(&BB_iter,&C_p_BB);

    for (int b=0; b<noccB_; b++){
      C_DAXPY(BB_iter.curr_size,1.0,&(C_p_BB.B_p_[0][b*noccB_+b]),
        noccB_*noccB_,&(diagBB_[off]),1);
    }

    off += BB_iter.curr_size;
  }

  C_p_BB.done();
}

}}

#include "sapt0.h"

using namespace boost;
using namespace std;
using namespace psi;

namespace psi { namespace sapt {

SAPT0::SAPT0(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt)
    : SAPT(options, psio, chkpt)
{
  maxiter_ = options_.get_int("MAXITER");
  e_conv_ = pow(10.0,-options_.get_int("E_CONVERGE"));
  d_conv_ = pow(10.0,-options_.get_int("D_CONVERGE"));

  psio_->open(PSIF_SAPT_AA_DF_INTS,PSIO_OPEN_NEW);
  psio_->open(PSIF_SAPT_BB_DF_INTS,PSIO_OPEN_NEW);
  psio_->open(PSIF_SAPT_AB_DF_INTS,PSIO_OPEN_NEW);

  print_header();
  df_integrals();
  w_integrals();
}

SAPT0::~SAPT0()
{
  free_block(wBAR_);
  free_block(wABS_);
  psio_->close(PSIF_SAPT_AA_DF_INTS,1);
  psio_->close(PSIF_SAPT_BB_DF_INTS,1);
  psio_->close(PSIF_SAPT_AB_DF_INTS,1);
}

double SAPT0::compute_energy()
{
  elst10();
  exch10();
  exch10_s2();
  if (debug_) ind20();
  ind20r();
  exch_ind20A_B();
  exch_ind20B_A();
  disp20();

  print_results();

  return (e_sapt0_);
}

void SAPT0::print_header()
{
  fprintf(outfile,"        SAPT0  \n");
  fprintf(outfile,"    Ed Hohenstein\n") ;
  fprintf(outfile,"     6 June 2009\n") ;
  fprintf(outfile,"\n");
  fprintf(outfile,"    Orbital Information\n");
  fprintf(outfile,"  -----------------------\n");
  fprintf(outfile,"    NSO     = %9d\n",nso_);
  fprintf(outfile,"    NMO     = %9d\n",nmo_);
  fprintf(outfile,"    NRI     = %9d\n",ndf_);
  fprintf(outfile,"    NOCC A  = %9d\n",noccA_);
  fprintf(outfile,"    NOCC B  = %9d\n",noccB_);
  fprintf(outfile,"    FOCC A  = %9d\n",foccA_);
  fprintf(outfile,"    FOCC B  = %9d\n",foccB_);
  fprintf(outfile,"    NVIR A  = %9d\n",nvirA_);
  fprintf(outfile,"    NVIR B  = %9d\n",nvirB_);
  fprintf(outfile,"\n");
  fflush(outfile);
}

void SAPT0::print_results()
{
  e_sapt0_ = eHF_ + e_disp20_ + e_exch_disp20_;
  double dHF = eHF_ - (e_elst10_ + e_exch10_ + e_ind20_ + e_exch_ind20_);

  fprintf(outfile,"\n    SAPT Results  \n");
  fprintf(outfile,"  ------------------------------------------------------------------\n");
  fprintf(outfile,"    E_HF          %16.8lf mH %16.8lf kcal mol^-1\n",
    eHF_*1000.0,eHF_*627.5095);
  fprintf(outfile,"    Elst10        %16.8lf mH %16.8lf kcal mol^-1\n",
    e_elst10_*1000.0,e_elst10_*627.5095);
  fprintf(outfile,"    Exch10        %16.8lf mH %16.8lf kcal mol^-1\n",
    e_exch10_*1000.0,e_exch10_*627.5095);
  fprintf(outfile,"    Exch10(S^2)   %16.8lf mH %16.8lf kcal mol^-1\n",
    e_exch10_s2_*1000.0,e_exch10_s2_*627.5095);
  fprintf(outfile,"    Ind20,r       %16.8lf mH %16.8lf kcal mol^-1\n",
    e_ind20_*1000.0,e_ind20_*627.5095);
  fprintf(outfile,"    Exch-Ind20,r  %16.8lf mH %16.8lf kcal mol^-1\n",
    e_exch_ind20_*1000.0,e_exch_ind20_*627.5095);
  fprintf(outfile,"    delta HF,r    %16.8lf mH %16.8lf kcal mol^-1\n",
    dHF*1000.0,dHF*627.5095);
  fprintf(outfile,"    Disp20        %16.8lf mH %16.8lf kcal mol^-1\n",
    e_disp20_*1000.0,e_disp20_*627.5095);
  fprintf(outfile,"    Exch-Disp20   %16.8lf mH %16.8lf kcal mol^-1\n\n",
    e_exch_disp20_*1000.0,e_exch_disp20_*627.5095);
  fprintf(outfile,"    Total SAPT0   %16.8lf mH %16.8lf kcal mol^-1\n",
    e_sapt0_*1000.0,e_sapt0_*627.5095);

  double tot_elst = e_elst10_;
  double tot_exch = e_exch10_;
  double tot_ind = e_ind20_ + e_exch_ind20_ + dHF;
  double tot_disp = e_disp20_ + e_exch_disp20_;

  Process::environment.globals["SAPT ELST ENERGY"] = tot_elst;
  Process::environment.globals["SAPT EXCH ENERGY"] = tot_exch;
  Process::environment.globals["SAPT IND ENERGY"] = tot_ind;
  Process::environment.globals["SAPT DISP ENERGY"] = tot_disp;
  Process::environment.globals["SAPT SAPT0 ENERGY"] = e_sapt0_;
  Process::environment.globals["SAPT ENERGY"] = e_sapt0_;
}

void SAPT0::df_integrals()
{
  psio_->open(PSIF_SAPT_TEMP,PSIO_OPEN_NEW);

  // Get fitting metric
  shared_ptr<FittingMetric> metric(new FittingMetric(ribasis_));
  metric->form_eig_inverse();
  double **J_temp = metric->get_metric()->pointer();
  double **J_mhalf = block_matrix(ndf_,ndf_);
  C_DCOPY(ndf_*ndf_,J_temp[0],1,J_mhalf[0],1);
  metric.reset();

  // Get Schwartz screening arrays
  double *Schwartz = init_array(basisset_->nshell()*(basisset_->nshell()+1)/2);

  shared_ptr<IntegralFactory> ao_eri_factory(new IntegralFactory(basisset_,
    basisset_, basisset_, basisset_));
  shared_ptr<TwoBodyAOInt> ao_eri = shared_ptr<TwoBodyAOInt>(
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
    }
  }

  ao_eri.reset();
  ao_eri_factory.reset();

  double *DFSchwartz = init_array(ribasis_->nshell());

  shared_ptr<IntegralFactory> df_eri_factory(new IntegralFactory(ribasis_,
    zero_, ribasis_, zero_));
  shared_ptr<TwoBodyAOInt> df_eri = shared_ptr<TwoBodyAOInt>(
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
  long int avail_mem = mem_ - (long int) ndf_*ndf_; 
  long int mem_tot = (long int) 2*ndf_*nsotri + (long int) ndf_*ndf_;
  if (avail_mem < (long int) 2*ndf_)
    throw PsiException("Not enough memory", __FILE__,__LINE__);
  long int max_size = avail_mem / ((long int) 2*ndf_);
  if (max_size > nsotri)
    max_size = nsotri;

//fprintf(outfile,"Requires storage of %ld doubles\n",mem_tot);
//fprintf(outfile,"Max nso x nso block is %ld\n\n",max_size);

  int size = 0;
  int num_blocks = 1;

  for(int P=0,PQ=0;P<basisset_->nshell();P++) {
    int numw = basisset_->shell(P)->nfunction();
    for(int Q=0;Q<=P;Q++,PQ++) {
      int numx = basisset_->shell(Q)->nfunction();
      int numPQ = numw*numx;
      if (P == Q) numPQ = numw*(numw+1)/2;

      size += numPQ;
      if (max_size < size) {
//      fprintf(outfile,"Block %d : %d\n",num_blocks,size-numPQ);
        num_blocks++;
        size = numPQ;
      }
  }} 

//fprintf(outfile,"Block %d : %d\n\n",num_blocks,size);

  int *PQ_start = init_int_array(num_blocks);
  int *PQ_stop = init_int_array(num_blocks);
  int *block_length =  init_int_array(num_blocks);
  int *PQ_offset = init_int_array(basisset_->nshell()*(basisset_->nshell()+1)/2);

  int block_num = 0;
  int totalPQ = 0;
  size = 0;

  PQ_start[0] = 0;

  for(int P=0,PQ=0;P<basisset_->nshell();P++) {
    int numw = basisset_->shell(P)->nfunction();
    for(int Q=0;Q<=P;Q++,PQ++) {
      int numx = basisset_->shell(Q)->nfunction();
      int numPQ = numw*numx;
      if (P == Q) numPQ = numw*(numw+1)/2;
      PQ_offset[PQ] = totalPQ;
      totalPQ += numPQ;

      size += numPQ;
      if (max_size < size) {
        PQ_stop[block_num] = PQ - 1;
        block_length[block_num] = size-numPQ;
        block_num++;
        PQ_start[block_num] = PQ;
        size = numPQ;
      }
  }} 

  PQ_stop[num_blocks-1] = basisset_->nshell()*(basisset_->nshell()+1)/2;
  block_length[num_blocks-1] = size;

//for (int i=0; i<num_blocks; i++) 
//  fprintf(outfile,"Block %2d : PQ %4d - %4d : %d\n",i,PQ_start[i],
//    PQ_stop[i],block_length[i]);
//fprintf(outfile,"\n");

  shared_ptr<IntegralFactory> rifactory(new IntegralFactory(ribasis_, zero_,
    basisset_, basisset_));

  int nthreads = 1;
  #ifdef _OPENMP
    nthreads = omp_get_max_threads();
  #endif
  int rank = 0;
  
  shared_ptr<TwoBodyAOInt> *eri = new shared_ptr<TwoBodyAOInt>[nthreads];
  const double **buffer = new const double*[nthreads];
  for(int i = 0;i < nthreads;++i){
    eri[i] = shared_ptr<TwoBodyAOInt>(rifactory->eri());
    buffer[i] = eri[i]->buffer();
  }

  zero_disk(PSIF_SAPT_TEMP,"AO RI Integrals",ndf_,nsotri);

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

      if (MU != NU) {
        munu_offset += nummu*numnu;
      }
      else {
        munu_offset += nummu*(nummu+1)/2;
      }

      if (PQ_stop[curr_block] == MUNU) {

        C_DGEMM('N','T',ndf_,block_length[curr_block],ndf_,1.0,J_mhalf[0],
          ndf_,&(AO_RI[0][0]),ndf_,0.0,&(J_AO_RI[0][0]),max_size);

        for (int P=0; P < ndf_; ++P) {
          next_DF_AO = psio_get_address(PSIO_ZERO,sizeof(double)*P*nsotri+
            sizeof(double)*offset);
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
    next_DF_AO = psio_get_address(PSIO_ZERO,sizeof(double)*P*nsotri+
      sizeof(double)*offset);
    psio_->write(PSIF_SAPT_TEMP,"AO RI Integrals",(char *)
      &(J_AO_RI[P][0]),sizeof(double)*block_length[curr_block],
      next_DF_AO,&next_DF_AO);
  }

  free_block(J_mhalf);
  free_block(AO_RI);
  free_block(J_AO_RI);
  free(Schwartz);
  free(DFSchwartz);

  avail_mem = mem_;
  long int indices = nsotri + noccA_*noccA_ + noccA_*nvirA_ + nvirA_*nvirA_
    + noccB_*noccB_ + noccB_*nvirB_ + nvirB_*nvirB_ + noccB_*noccB_ 
    + noccA_*nvirB_ + noccB_*nvirA_;
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

  double **B_p_munu = block_matrix(Plength,nsotri);
  double **B_p_AA = block_matrix(Plength,noccA_*noccA_);
  double **B_p_AR = block_matrix(Plength,noccA_*nvirA_);
  double **B_p_RR = block_matrix(Plength,nvirA_*nvirA_);
  double **B_p_BB = block_matrix(Plength,noccB_*noccB_);
  double **B_p_BS = block_matrix(Plength,noccB_*nvirB_);
  double **B_p_SS = block_matrix(Plength,nvirB_*nvirB_);
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
  zero_disk(PSIF_SAPT_AA_DF_INTS,"RR RI Integrals",ndf_,nvirA_*nvirA_);
  zero_disk(PSIF_SAPT_BB_DF_INTS,"BB RI Integrals",ndf_,noccB_*noccB_);
  zero_disk(PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",ndf_,noccB_*nvirB_);
  zero_disk(PSIF_SAPT_BB_DF_INTS,"SS RI Integrals",ndf_,nvirB_*nvirB_);
  zero_disk(PSIF_SAPT_AB_DF_INTS,"AB RI Integrals",ndf_,noccA_*noccB_);
  zero_disk(PSIF_SAPT_AB_DF_INTS,"AS RI Integrals",ndf_,noccA_*nvirB_);
  zero_disk(PSIF_SAPT_AB_DF_INTS,"RB RI Integrals",ndf_,nvirA_*noccB_);

  #ifdef HAVE_MKL
    int mkl_nthread = mkl_get_max_threads();
    mkl_set_num_threads(1);
  #endif

  int Prel;

  for (int Pbl=0; Pbl<Pblocks; Pbl++) {

    int length = max_size;
    if (gimp && Pbl == Pblocks-1) length = gimp;

    psio_->read(PSIF_SAPT_TEMP,"AO RI Integrals",(char *) &(B_p_munu[0][0]),
      sizeof(double)*length*nsotri,next_DF_AO,&next_DF_AO);

    #pragma omp for private(Prel,rank) schedule(dynamic)
    for (Prel=0; Prel<length; Prel++) {

      #ifdef _OPENMP
        rank = omp_get_thread_num();
      #endif

      munu_offset = 0;
      for(int MU=0,MUNU=0;MU<basisset_->nshell();MU++) {
        int nummu = basisset_->shell(MU)->nfunction();
        for(int NU=0;NU<=MU;NU++,MUNU++) {
          int numnu = basisset_->shell(NU)->nfunction();
  
          if (MU != NU) {
            for (int mu=0,munu=0; mu < nummu; ++mu) {
              int omu = basisset_->shell(MU)->function_index() + mu;
  
              for (int nu=0; nu < numnu; ++nu, ++munu) {
                int onu = basisset_->shell(NU)->function_index() + nu;
  
                munu_temp[rank][omu*nso_+onu] = 
                  B_p_munu[Prel][munu+PQ_offset[MUNU]];
                munu_temp[rank][onu*nso_+omu] = 
                  B_p_munu[Prel][munu+PQ_offset[MUNU]];
              }
            }
          }
          else {
            for (int mu=0,munu=0; mu < nummu; ++mu) {
              int omu = basisset_->shell(MU)->function_index() + mu;
  
              for (int nu=0; nu <= mu; ++nu, ++munu) {
                int onu = basisset_->shell(NU)->function_index() + nu;
  
                munu_temp[rank][omu*nso_+onu] = 
                  B_p_munu[Prel][munu+PQ_offset[MUNU]];
                munu_temp[rank][onu*nso_+omu] = 
                  B_p_munu[Prel][munu+PQ_offset[MUNU]];
              }
            }
          }
        }
      }
  
      C_DGEMM('T', 'N', nmo_, nso_, nso_, 1.0, &(CA_[0][0]), nmo_,
        munu_temp[rank], nso_, 0.0, Inu_temp[rank], nso_);
      C_DGEMM('N', 'N', nmo_, nmo_, nso_, 1.0, Inu_temp[rank], nso_,
        &(CA_[0][0]), nmo_, 0.0, IJ_temp[rank], nmo_);
  
      for (int a=0; a<noccA_; a++) {
        C_DCOPY(noccA_,&(IJ_temp[rank][a*nmo_]),1,&(B_p_AA[Prel][a*noccA_]),1);
        C_DCOPY(nvirA_,&(IJ_temp[rank][a*nmo_+noccA_]),1,
          &(B_p_AR[Prel][a*nvirA_]),1);
      }
      for (int r=0; r<nvirA_; r++) {
        C_DCOPY(nvirA_,&(IJ_temp[rank][(r+noccA_)*nmo_+noccA_]),1,
          &(B_p_RR[Prel][r*nvirA_]),1);
      }
  
      C_DGEMM('N', 'N', nmo_, nmo_, nso_, 1.0, Inu_temp[rank], nso_,
        &(CB_[0][0]), nmo_, 0.0, IJ_temp[rank], nmo_);
  
      for (int a=0; a<noccA_; a++) {
        C_DCOPY(noccB_,&(IJ_temp[rank][a*nmo_]),1,&(B_p_AB[Prel][a*noccB_]),1);
        C_DCOPY(nvirB_,&(IJ_temp[rank][a*nmo_+noccB_]),1,
          &(B_p_AS[Prel][a*nvirB_]),1);
      }
      for (int r=0; r<nvirA_; r++) {
        C_DCOPY(noccB_,&(IJ_temp[rank][(r+noccA_)*nmo_]),1,
          &(B_p_RB[Prel][r*noccB_]),1);
      } 
  
      C_DGEMM('T', 'N', nmo_, nso_, nso_, 1.0, &(CB_[0][0]), nmo_,
        munu_temp[rank], nso_, 0.0, Inu_temp[rank], nso_);
      C_DGEMM('N', 'N', nmo_, nmo_, nso_, 1.0, Inu_temp[rank], nso_,
        &(CB_[0][0]), nmo_, 0.0, IJ_temp[rank], nmo_);
  
      for (int b=0; b<noccB_; b++) {
        C_DCOPY(noccB_,&(IJ_temp[rank][b*nmo_]),1,&(B_p_BB[Prel][b*noccB_]),1);
        C_DCOPY(nvirB_,&(IJ_temp[rank][b*nmo_+noccB_]),1,
          &(B_p_BS[Prel][b*nvirB_]),1);
      }
      for (int s=0; s<nvirB_; s++) {
        C_DCOPY(nvirB_,&(IJ_temp[rank][(s+noccB_)*nmo_+noccB_]),1,
          &(B_p_SS[Prel][s*nvirB_]),1);
      } 
 
    }

    psio_->write(PSIF_SAPT_AA_DF_INTS,"AA RI Integrals",(char *)
      &(B_p_AA[0][0]),sizeof(double)*length*noccA_*noccA_,
      next_DF_AA,&next_DF_AA);
    psio_->write(PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",(char *)
      &(B_p_AR[0][0]),sizeof(double)*length*noccA_*nvirA_,
      next_DF_AR,&next_DF_AR);
    psio_->write(PSIF_SAPT_AA_DF_INTS,"RR RI Integrals",(char *)
      &(B_p_RR[0][0]),sizeof(double)*length*nvirA_*nvirA_,
      next_DF_RR,&next_DF_RR);

    psio_->write(PSIF_SAPT_BB_DF_INTS,"BB RI Integrals",(char *)
      &(B_p_BB[0][0]),sizeof(double)*length*noccB_*noccB_,
      next_DF_BB,&next_DF_BB);
    psio_->write(PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",(char *)
      &(B_p_BS[0][0]),sizeof(double)*length*noccB_*nvirB_,
      next_DF_BS,&next_DF_BS);
    psio_->write(PSIF_SAPT_BB_DF_INTS,"SS RI Integrals",(char *)
      &(B_p_SS[0][0]),sizeof(double)*length*nvirB_*nvirB_,
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

  #ifdef HAVE_MKL
    mkl_set_num_threads(mkl_nthread);
  #endif

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
  free(PQ_start);
  free(PQ_stop);
  free(block_length);
  free(PQ_offset);

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

}}

#include "drpa.h"

#include <time.h>
#include <libqt/qt.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libmints/mints.h>

using namespace boost;
using namespace psi;

namespace psi { namespace dfcc {

dRPA::dRPA(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt)
  : CC(options, psio, chkpt)
{
  print_header();
  psio_->open(DFCC_INT_FILE,PSIO_OPEN_NEW);
  df_integrals();
}

dRPA::~dRPA()
{
  psio_->close(DFCC_INT_FILE,1);
}

double dRPA::compute_energy()
{
  double energy;

  if (options_.get_str("RPA_ALGORITHM") == "DF")
    energy = df_compute_energy();
  if (options_.get_str("RPA_ALGORITHM") == "CD")
    energy = cd_compute_energy();

  return (energy);
}

double dRPA::df_compute_energy()
{
  time_t start = time(NULL);
  time_t stop;

  if (options_.get_bool("DIIS"))
    diis_ = shared_ptr<DFCCDIIS>(new DFCCDIIS(DFCC_DIIS_FILE,naocc_*naocc_*
      navir_*navir_,options_.get_int("MAX_DIIS_VECS"),psio_));

  B_p_IA_ = block_matrix(naocc_*navir_,ndf_);
  Th_p_IA_ = block_matrix(naocc_*navir_,ndf_);
  tIAJB_ = init_array((long int) naocc_*navir_*naocc_*navir_);
  t2IAJB_ = init_array((long int) naocc_*navir_*naocc_*navir_);
  xIAJB_ = init_array((long int) naocc_*navir_*naocc_*navir_);

  psio_->read_entry(DFCC_INT_FILE,"OV DF Integrals",(char *) &(B_p_IA_[0][0]),
    sizeof(double)*naocc_*navir_*ndf_);

  C_DGEMM('N','T',naocc_*navir_,naocc_*navir_,ndf_,1.0,B_p_IA_[0],ndf_,
    B_p_IA_[0],ndf_,0.0,tIAJB_,naocc_*navir_);

  apply_denom();
  C_DCOPY((long int) naocc_*navir_*naocc_*navir_,tIAJB_,1,t2IAJB_,1);

  C_DGEMM('N','N',naocc_*navir_,ndf_,naocc_*navir_,1.0,tIAJB_,naocc_*navir_,
    B_p_IA_[0],ndf_,0.0,Th_p_IA_[0],ndf_); 

  double emp2 = df_energy();

  fprintf(outfile,"  Reference Energy            %18.10lf\n",Eref_);
  fprintf(outfile,"  Correlation Energy          %18.10lf\n",emp2);
  fprintf(outfile,"  Total DF-MP2 (J) Energy     %18.10lf\n\n",Eref_+emp2);
  fprintf(outfile,"  Iter       Energy (H)          dE (H)             RMS (H)     Time (s)\n");
  fflush(outfile);

  int iter = 1;
  int done = 0;
  double e_new;
  double e_old = emp2;
  double rms = 1.0;

  do {

    C_DGEMM('N','T',naocc_*navir_,naocc_*navir_,ndf_,1.0,B_p_IA_[0],ndf_,
      B_p_IA_[0],ndf_,0.0,tIAJB_,naocc_*navir_);
    C_DGEMM('N','T',naocc_*navir_,naocc_*navir_,ndf_,1.0,B_p_IA_[0],ndf_,
      Th_p_IA_[0],ndf_,1.0,tIAJB_,naocc_*navir_);
    C_DGEMM('N','T',naocc_*navir_,naocc_*navir_,ndf_,1.0,Th_p_IA_[0],ndf_,
      B_p_IA_[0],ndf_,1.0,tIAJB_,naocc_*navir_);
    C_DGEMM('N','T',naocc_*navir_,naocc_*navir_,ndf_,1.0,Th_p_IA_[0],ndf_,
      Th_p_IA_[0],ndf_,1.0,tIAJB_,naocc_*navir_);
  
    apply_denom();
    
    C_DGEMM('N','N',naocc_*navir_,ndf_,naocc_*navir_,1.0,tIAJB_,naocc_*navir_,
      B_p_IA_[0],ndf_,0.0,Th_p_IA_[0],ndf_);
    
    e_new = df_energy();

    rms = df_store_error_vecs();

    stop = time(NULL);
    fprintf(outfile,"  %4d %16.8lf %17.9lf %17.9lf %12ld",iter,e_new,
      e_old-e_new,rms,stop-start);
    fflush(outfile);

    if (options_.get_int("MIN_DIIS_VECS") <= iter &&
        options_.get_bool("DIIS")) {
      diis_->get_new_vector(tIAJB_,xIAJB_);
      fprintf(outfile,"  DIIS\n");
      C_DGEMM('N','N',naocc_*navir_,ndf_,naocc_*navir_,1.0,tIAJB_,
        naocc_*navir_,B_p_IA_[0],ndf_,0.0,Th_p_IA_[0],ndf_);
    }
    else {
      fprintf(outfile,"\n");
    }
    fflush(outfile);

    iter++;

    if (iter > options_.get_int("MAXITER")) done = 1;
    if (fabs(e_old-e_new) < pow(10.0,-(double) options_.get_int("E_CONVERGE"))
      && rms < pow(10.0,-(double) options_.get_int("T_CONVERGE"))) done = 1;

    e_old = e_new;
  }
  while(!done);

  fprintf(outfile,"\n");
  fprintf(outfile,"  Reference Energy            %18.10lf\n",Eref_);
  fprintf(outfile,"  Correlation Energy          %18.10lf\n",e_old);
  fprintf(outfile,"  Total DF-dRPA Energy        %18.10lf\n\n",Eref_+e_old);

  free_block(B_p_IA_);
  free_block(Th_p_IA_);
  free(tIAJB_);
  free(xIAJB_);
  free(t2IAJB_);

  return(0.0);
}

double dRPA::cd_compute_energy()
{
  time_t start = time(NULL);
  time_t stop;

  double **B_p_IA = block_matrix(naocc_*navir_,ndf_);

  psio_->read_entry(DFCC_INT_FILE,"OV DF Integrals",(char *)
      &(B_p_IA[0][0]),naocc_*navir_*ndf_*sizeof(double));

  double emp2 = 0.0;

  fprintf(outfile,"  Reference Energy            %18.10lf\n",Eref_);
  fprintf(outfile,"  Correlation Energy          %18.10lf\n",emp2);
  fprintf(outfile,"  Total DF-MP2 (J) Energy     %18.10lf\n\n",Eref_+emp2);
  fprintf(outfile,"  Iter       Energy (H)          dE (H)             RMS (H)     Time (s)\n");
  fflush(outfile);

  int iter = 1;
  int done = 0;
  double e_new;
  double e_old = emp2;
  double rms;

  do {

    stop = time(NULL);
    fprintf(outfile,"  %4d %16.8lf %17.9lf %17.9lf %12ld\n",iter,e_new,
      e_old-e_new,rms,stop-start);
    fflush(outfile);
    iter++;

    if (iter > options_.get_int("MAXITER")) done = 1;
    if (fabs(e_old-e_new) < pow(10.0,-(double) options_.get_int("E_CONVERGE")))
      done = 1;
    if (rms < pow(10.0,-(double) options_.get_int("T_CONVERGE"))) done = 1;

    e_old = e_new;
  }
  while(!done);

  fprintf(outfile,"\n");
  fprintf(outfile,"  Reference Energy            %18.10lf\n",Eref_);
  fprintf(outfile,"  Correlation Energy          %18.10lf\n",e_old);
  fprintf(outfile,"  Total DF-dRPA Energy         %18.10lf\n\n",Eref_+e_old);

  return(0.0);
}

void dRPA::print_header()
{
    fprintf(outfile, "\t********************************************************\n");
    fprintf(outfile, "\t*                                                      *\n");
    fprintf(outfile, "\t*                        dRPA                          *\n");
    fprintf(outfile, "\t*           Direct Random Phase Approximation          *\n");
    fprintf(outfile, "\t*                with all sorts of shit                *\n");
    fprintf(outfile, "\t*                                                      *\n");
    fprintf(outfile, "\t*            Rob Parrish and Ed Hohenstein             *\n");
    fprintf(outfile, "\t*                                                      *\n");
    fprintf(outfile, "\t********************************************************\n");
    fprintf(outfile, "\n");
    CC::print_header();

}

void dRPA::apply_denom()
{
  for (int i=0,ia=0; i<naocc_; i++) {
  for (int a=0; a<navir_; a++,ia++) {
    for (int j=0,jb=0; j<naocc_; j++) {
    for (int b=0; b<navir_; b++,jb++) {
      long int iajb = (long int) ia*naocc_*navir_ + jb;
      tIAJB_[iajb] /= evals_aoccp_[i] + evals_aoccp_[j] - evals_avirp_[a]
        - evals_avirp_[b];
    }}
  }}
}

double dRPA::df_energy()
{
  return(2.0*C_DDOT(naocc_*navir_*ndf_,Th_p_IA_[0],1,B_p_IA_[0],1));
}

double dRPA::df_store_error_vecs()
{
  double rms;

  C_DAXPY((long int) naocc_*navir_*naocc_*navir_,-1.0,tIAJB_,1,t2IAJB_,1);

  if (options_.get_bool("DIIS"))
    diis_->store_vectors(tIAJB_,t2IAJB_);

  rms = C_DDOT((long int) naocc_*navir_*naocc_*navir_,t2IAJB_,1,t2IAJB_,1);
  rms /= (double) naocc_*navir_*naocc_*navir_;

  C_DCOPY((long int) naocc_*navir_*naocc_*navir_,tIAJB_,1,t2IAJB_,1);

  return(sqrt(rms));
}

}}


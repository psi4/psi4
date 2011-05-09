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

}}


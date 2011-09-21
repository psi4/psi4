/*! \file
    \ingroup ccresponse
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include <libpsio/psio.h>
#include <psifiles.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccresponse {

void init_X(const char *pert, int irrep, double omega);
void sort_X(const char *pert, int irrep, double omega);
void cc2_sort_X(const char *pert, int irrep, double omega);
void X1_build(const char *pert, int irrep, double omega);
void X2_build(const char *pert, int irrep, double omega);
void cc2_X1_build(const char *pert, int irrep, double omega);
void cc2_X2_build(const char *pert, int irrep, double omega);
double converged(const char *pert, int irrep, double omega);
void save_X(const char *pert, int irrep, double omega);
void print_X(const char *pert, int irrep, double omega);
void update_X(const char *pert, int irrep, double omega);
void diis(int iter, const char *pert, int irrep, double omega);
double pseudopolar(const char *pert, int irrep, double omega);
void cleanup(void);
void exit_io(void);
void amp_write(const char *pert, int irrep, double omega);

void analyze(const char *pert, int irrep, double omega);

void compute_X(const char *pert, int irrep, double omega)
{
  int i, iter=0, done=0;
  double rms, polar, X2_norm;
  char lbl[32];
  dpdbuf4 X2;

  timer_on("compute_X");

  fprintf(outfile, "\n\tComputing %s-Perturbed Wave Function (%5.3f E_h).\n", pert, omega);
  init_X(pert, irrep, omega);
  fprintf(outfile, "\tIter   Pseudopolarizability       RMS \n");
  fprintf(outfile, "\t----   --------------------   -----------\n");
  fflush(outfile);

  if (params.wfn == "CC2")
    cc2_sort_X(pert, irrep, omega);
  else
    sort_X(pert, irrep, omega);
  polar = -2.0*pseudopolar(pert, irrep, omega);
  fprintf(outfile, "\t%4d   %20.12f\n", iter, polar);
  fflush(outfile);

  for(iter=1; iter <= params.maxiter; iter++) {

    if (params.wfn == "CC2") {
      cc2_sort_X(pert, irrep, omega);
      cc2_X1_build(pert, irrep, omega);
      cc2_X2_build(pert, irrep, omega);
    }
    else {
      sort_X(pert, irrep, omega);
      X1_build(pert, irrep, omega);
      X2_build(pert, irrep, omega);
    }
    update_X(pert, irrep, omega);
    rms = converged(pert, irrep, omega);
    if(rms <= params.convergence) {
      done = 1;
      save_X(pert, irrep, omega);
      if (params.wfn == "CC2")
	cc2_sort_X(pert, irrep, omega);
      else
	sort_X(pert, irrep, omega);
      fprintf(outfile, "\t-----------------------------------------\n");
      fprintf(outfile, "\tConverged %s-Perturbed Wfn to %4.3e\n", pert, rms);
      if(params.print & 2) {
	sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
	dpd_buf4_init(&X2, CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
	X2_norm = dpd_buf4_dot_self(&X2);
	dpd_buf4_close(&X2);
	X2_norm = sqrt(X2_norm);
	fprintf(outfile, "\tNorm of the converged X2 amplitudes %20.15f\n", X2_norm);
	amp_write(pert, irrep, omega);
      }
      fflush(outfile);
      break;
    }
    if(params.diis) diis(iter, pert, irrep, omega);
    save_X(pert, irrep, omega);
    if (params.wfn == "CC2")
      cc2_sort_X(pert, irrep, omega);
    else
      sort_X(pert, irrep, omega);

    polar = -2.0*pseudopolar(pert, irrep, omega);
    fprintf(outfile, "\t%4d   %20.12f    %4.3e\n", iter, polar, rms);
    fflush(outfile);

  }
  if(!done) {
    fflush(outfile);
    dpd_close(0);
    cleanup();
    exit_io();
    throw PsiException("Failed to converge perturbed wavefunction",__FILE__,__LINE__);
  }

  /* Clean up disk space */
  psio_close(CC_DIIS_AMP, 0);
  psio_close(CC_DIIS_ERR, 0);

  psio_open(CC_DIIS_AMP, 0);
  psio_open(CC_DIIS_ERR, 0);

  for(i=CC_TMP; i <= CC_TMP11; i++) {
    psio_close(i,0);
    psio_open(i,0);
  }

  if(params.analyze) analyze(pert, irrep, omega);

  /*  print_X(pert, irrep, omega); */

  timer_off("compute_X");
}

}} // namespace psi::ccresponse

/*! \file
    \ingroup CPHF
    \brief Enter brief description of file here 
*/
/*
** CPHF: Program to solve the Coupled Perturbed Hartree-Fock equations
** for nuclear and electric field perturbations and to compute
** electric polarizabilities, harmonic vibrational frequencies, and IR
** intensities.
**
** Limitations of and future plans for this code:
**
** (1) Spin-restricted closed-shell Hartree-Fock (RHF) wave functions
** only.  Extension to ROHF and UHF cases is needed.
**
** (2) All two-electron integrals are held in core and used in the MO
** basis.  Out-of-core and AO-based algorithms are needed in mohess.c,
** cphf_X.c, and build_hessian.c in order to handle larger basis sets
** and to avoid the two-electron integral transformation.
**
** (3) Symmetry-blocking is used in most of the loops, but is not
** actually used to improve storage or computational order.  I've put
** this off because the nuclear perturbations (x, y, and z on each
** nucleus) are not yet symmetry-adapted.  Some effort in this area is
** needed.
**
** (4) Thermodynamic functions, including enthalpies, entropies, heat
** capacities, free energies, and partition functions can be computed,
** given the vibrational data computed in vibration.c.
**
** TDC, December 2001 and October 2002
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libiwl/iwl.h>
#include <libchkpt/chkpt.h>
#include <libqt/qt.h>
#include <psifiles.h>
#include "globals.h"

namespace psi { namespace cphf {

void init_io(int argc, char *argv[]);
void exit_io(void);
void title(void);
void init_ioff(void);

void setup(void);
void cleanup(void);

void out_of_core(double ***, double ***, double ***, double **);
void sort_B(double ***, double **);
void sort_A(double **, double **);
void mohess(double **);
void cphf_X(double ***, double **, double **, double ***);
void cphf_F(double **, double ***);
void polarize(double ***);
void build_hessian(double ***, double ***, double **, double ***, double **);
void build_dipder(double ***);
void vibration(double **, double **);
void cphf_B(double ***, double **);

}} // namespace psi::cphf

int main(int argc, char *argv[])
{
  using namespace psi::cphf;
  int errcod = 0;
  int coord = 0;
  double ***F;
  double ***S;
  double ***B;
  double **A;
  double **Baijk;
  double **Aaibj;
  double ***UX;
  double ***UF;
  double **hessian;
  double **L;
  
  init_io(argc,argv);
  title();
  init_ioff();

  timer_init();
  timer_on("CPHF Main");

  print_lvl = 0;
  errcod = ip_data("PRINT", "%d", &(print_lvl), 0);

  setup();

  F = (double ***) malloc(natom*3 * sizeof(double **));
  S = (double ***) malloc(natom*3 * sizeof(double **));
  B = (double ***) malloc(natom*3 * sizeof(double **));
  for(coord=0; coord < natom*3; coord++) {
    F[coord] = block_matrix(nmo, nmo);
    S[coord] = block_matrix(nmo, nmo);
    B[coord] = block_matrix(nmo, nmo);
  }

  A = block_matrix(num_pq,num_pq);

  out_of_core(F, S, B, A);

  Baijk = block_matrix(natom*3, num_ai);
  
  sort_B(B, Baijk);

  for(coord=0; coord < natom*3; coord++) {
    free_block(B[coord]);
  }
  free(B);

  Aaibj = block_matrix(num_ai,num_ai);
  
  sort_A(A, Aaibj);
  
  mohess(Aaibj);
  
  UX = (double ***) malloc(natom*3 * sizeof(double **));
  for(coord=0; coord < natom*3; coord++) {
    UX[coord] = block_matrix(nmo,nmo);
  }
  
  cphf_X(S, Baijk, Aaibj, UX);

  if (X_only) { /* only compute cphf coefficients, then quit */
    free_block(A);
    free_block(Baijk);
    free_block(Aaibj);

    for(coord=0; coord < natom*3; coord++) {
      free_block(UX[coord]);
      free_block(F[coord]);
      free_block(S[coord]);
    }
    free(UX); free(F); free(S);
    timer_off("CPHF Main");
    timer_done();
    exit_io();
    exit(PSI_RETURN_SUCCESS);
  }
  
  UF = (double ***) malloc(3 * sizeof(double **));
  for(coord=0; coord < 3; coord++)
    UF[coord] = block_matrix(nmo,nmo);

  cphf_F(Aaibj, UF);
  
  polarize(UF);
  
  hessian = block_matrix(natom*3, natom*3);
  build_hessian(F, S, A, UX, hessian);
  
  dipder = block_matrix(3, natom*3);
  dipder_q = block_matrix(3, natom*3);
  build_dipder(UX);
  
  L = block_matrix(natom*3, natom*3);
  vibration(hessian,L);

  cphf_B(UX,L);
  
  cleanup(); 

  free_block(A);
  free_block(Baijk);
  free_block(Aaibj);

  for(coord=0; coord < natom*3; coord++) { 
    free_block(UX[coord]);
    free_block(F[coord]);
    free_block(S[coord]);
  }
  for(coord=0; coord < 3; coord++) { 
    free_block(UF[coord]);
  }
  free(UX); free(UF); free(F); free(S);
  free_block(hessian);
  free_block(dipder);
  free_block(L);
  
  timer_off("CPHF Main");
  timer_done();

  exit_io();
  exit(PSI_RETURN_SUCCESS);
}

namespace psi { namespace cphf {

void init_io(int argc, char *argv[])
{
  int i, num_unparsed;
  char *progid, *argv_unparsed[100];

  X_only = 0;
  for (i=1, num_unparsed=0; i<argc; ++i) {
    if (!strcmp(argv[i],"--X_only")) /* only compute cphf coefficients, then quit */
      X_only = 1;
    else
      argv_unparsed[num_unparsed++] = argv[i];
  }

  psi_start(&infile,&outfile,&psi_file_prefix,num_unparsed,argv_unparsed,0);
  ip_cwk_add(progid);
  free(progid);
  tstart(outfile);

}

void title(void)
{
  fprintf(outfile, "\t\t\t**************************\n");
  fprintf(outfile, "\t\t\t*                        *\n");
  fprintf(outfile, "\t\t\t*          CPHF          *\n");
  fprintf(outfile, "\t\t\t*                        *\n");
  fprintf(outfile, "\t\t\t**************************\n");
}

void exit_io(void)
{
  int i;
  psio_done();
  tstop(outfile);
  psi_stop(infile,outfile,psi_file_prefix);
}

void init_ioff(void)
{
  int i;
  ioff = (int *) malloc(IOFF_MAX * sizeof(int));
  ioff[0] = 0;
  for(i=1; i < IOFF_MAX; i++) {
      ioff[i] = ioff[i-1] + i;
  }
}

}} // namespace psi::cphf

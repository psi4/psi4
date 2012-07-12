/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <stdio.h>
#include <math.h>
#include "psi4-dec.h"
#include <libmints/wavefunction.h>
#include <libmints/basisset.h>
#include <libmints/matrix.h>
#include <libmints/vector.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libiwl/iwl.h>
#include <libchkpt/chkpt.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include <psifiles.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

/*
** densgrid_RHF(): Compute the values of the ground-state one-particle 
** density at a set of grid points.
**
** TDC, 7/2012
*/
void compute_delta(double **delta, double x, double y, double z);
void setup_delta(void);
int nmo, nao; // global
double **scf, **u; // global

void densgrid_RHF(void)
{
  int nmo = moinfo.nmo;
  int nfzv = moinfo.nfzv;
  int nirreps = moinfo.nirreps;
  int nactive = nmo - nfzv;
  double dens;
  double **D, **delta;
  double xmin, xmax, ymin, ymax, zmin, zmax;
  double xstep, ystep, zstep;

  boost::shared_ptr<Molecule> molecule = Process::environment.molecule();

  D = moinfo.opdm; // A block matrix
  delta = block_matrix(nactive, nactive); // Dirac delta function 

  xmin = -3.0; xmax = +3.0; xstep = 0.1;
  ymin = -3.0; ymax = +3.0; ystep = 0.1;
  zmin = -3.0; zmax = +3.0; zstep = 0.1;

  // Loop over points
  for(double x=xmin; x <= xmax; x += xstep) {
    for(double y=ymin; y <= ymax; y += ystep) {
      for(double z=zmin; z <= zmax; z += zstep) {

        // Compute delta function in Gaussian basis
        compute_delta(delta, x, y, z);

        dens = 0.0;
        for(int i=0; i < nactive; i++)
          for(int j=0; j < nactive; j++)
            dens += delta[i][j] * D[i][j];

      } // z
    }  // y
  } // x

  free_block(delta);

}

void compute_delta(double **delta, double x, double y, double z)
{
  int i, j;
  double *phi_ao, *phi_so, *phi_mo;

  setup_delta();

  phi_ao = init_array(nao);  /* AO function values */
  phi_so = init_array(nmo);  /* SO function values */
  phi_mo = init_array(nmo);  /* MO function values */

  boost::shared_ptr<BasisSet> basis = Process::environment.reference_wavefunction()->basisset();
  basis->compute_phi(phi_ao, x, y, z);

  /*  for(i=0; i < nao; i++) printf("%d %20.10f\n", i, phi_ao[i]); */

  /* Transform the basis function values to the MO basis */
  C_DGEMV('n', nmo, nao, 1.0, &(u[0][0]), nao, &(phi_ao[0]), 1,
          0.0, &(phi_so[0]), 1);

  C_DGEMV('t', nmo, nmo, 1.0, &(scf[0][0]), nmo, &(phi_so[0]), 1,
          0.0, &(phi_mo[0]), 1);

  /* for(i=0; i < nmo; i++) printf("%d %20.10f\n", i, phi_mo[i]); */


  /* Build the MO-basis delta function */
  for(i=0; i < nmo; i++)
    for(j=0; j < nmo; j++)
      delta[i][j] = phi_mo[i] * phi_mo[j];

  free(phi_ao);
  free(phi_so);
  free(phi_mo);
}

void setup_delta(void)
{
  static int done=0;
  int i, I, j;
  int nirreps, nfzc, nfzv;
  int *order, *clsdpi, *openpi, *orbspi, *fruocc, *frdocc;
  double **scf_pitzer;

  if(done) return;

  chkpt_init(PSIO_OPEN_OLD);
  nmo = chkpt_rd_nmo();
  nao = chkpt_rd_nao();
  nirreps = chkpt_rd_nirreps();
  clsdpi = chkpt_rd_clsdpi();
  openpi = chkpt_rd_openpi();
  orbspi = chkpt_rd_orbspi();
  scf_pitzer = chkpt_rd_scf();
  u = chkpt_rd_usotao();
  chkpt_close();

  moinfo.frdocc = Process::environment.reference_wavefunction()->frzcpi();
  moinfo.fruocc = Process::environment.reference_wavefunction()->frzvpi();

  nfzc = nfzv = 0;
  for(i=0; i < nirreps; i++) {
    nfzc += frdocc[i];
    nfzv += fruocc[i];
  }

  /*** Get the Pitzer -> QT reordering array ***/
  order = init_int_array(nmo);
  reorder_qt(clsdpi, openpi, frdocc, fruocc, order, orbspi, nirreps);

  /*** Arrange the SCF eigenvectors into QT ordering ***/
  scf = block_matrix(nmo, nmo);
  for(i=0; i < nmo; i++) {
      I = order[i];  /* Pitzer --> QT */
      for(j=0; j < nmo; j++) scf[j][I] = scf_pitzer[j][i];
    }

  free(order);
  free(clsdpi);
  free(openpi);
  free(orbspi);
  free(fruocc);
  free(frdocc);
  free_block(scf_pitzer);

  done = 1;

  return;
}

}} // namespace psi::ccdensity

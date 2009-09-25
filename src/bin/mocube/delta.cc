/*! \file
    \ingroup MOCUBE
    \brief Enter brief description of file here 
*/

/* * returns values of mos at a grid point, RAK, Nov. 2002
** hacked from cusp: delta() ** TDC, June 2001
*/

#include <cstdio>
#include <cstdlib>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libqt/qt.h>
#define EXTERN
#include "mocube.h"

namespace psi { namespace mocube {

void compute_phi(double *phi, double x, double y, double z);

void compute_mos(double *movals, double x, double y, double z,
    double **scf, double **u)
{
  int i, j, nmo, nao;
  double *phi_ao, *phi_so, *phi_mo, tval;

  nmo = params.nmo;
  nao = params.nao;

  /* setup_delta(); */

  phi_ao = init_array(nao);  /* AO function values */
  phi_so = init_array(nmo);  /* SO function values */
  phi_mo = init_array(nmo);  /* MO function values */

  compute_phi(phi_ao, x, y, z);

  /* Transform the basis function values to the MO basis */
  C_DGEMV('n', nmo, nao, 1.0, &(u[0][0]), nao, &(phi_ao[0]), 1,
	  0.0, &(phi_so[0]), 1);

  C_DGEMV('t', nmo, nmo, 1.0, &(scf[0][0]), nmo, &(phi_so[0]), 1,
	  0.0, &(phi_mo[0]), 1);

  free(phi_ao);
  free(phi_so);

  for (i=0; i<cube.nmo_to_plot; ++i)
    movals[i] = phi_mo[ cube.mos_to_plot[i] ];

  free(phi_mo);
  return;
}

}} // namespace psi::mocube

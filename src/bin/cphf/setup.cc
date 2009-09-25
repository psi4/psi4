/*! \file
    \ingroup CPHF
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#define EXTERN
#include "globals.h"

namespace psi { namespace cphf {

void setup(void)
{
  int i, coord, h, foffset, loffset;

  chkpt_init(PSIO_OPEN_OLD);
  natom = chkpt_rd_natom();
  nmo = chkpt_rd_nmo();
  nso = chkpt_rd_nso();
  nao = chkpt_rd_nao();
  evals = chkpt_rd_evals();
  nirreps = chkpt_rd_nirreps();
  orbspi = chkpt_rd_orbspi();
  clsdpi = chkpt_rd_clsdpi();
  zvals = chkpt_rd_zvals();
  scf = chkpt_rd_scf();
  usotao = chkpt_rd_usotao();
  geom = chkpt_rd_geom();
  rottype = chkpt_rd_rottype();
  chkpt_close();

  // Number of normal coordinates
  if(rottype==3) nnc = (3*natom-5);
  else nnc = (3*natom-6);

  asymbol = (char**)malloc(natom*3*sizeof(char*));
  for(i=0; i<natom*3; i++)
    asymbol[i] = (char*)malloc(3*sizeof(char));

  for(coord=0; coord<natom; coord++) {
    for(i=0; i<3; i++) {
      zval_to_symbol(zvals[coord],asymbol[coord*3+i]);
    }
  }

  ntri = nmo * (nmo + 1)/2;
  ntei = ntri * (ntri + 1)/2;
  noei = nso * (nso + 1)/2;
  noei_ao = nao * (nao + 1)/2;

  uoccpi = init_int_array(nirreps);
  for(h=0; h < nirreps; h++) uoccpi[h] = orbspi[h] - clsdpi[h];

  ndocc = 0; nuocc = 0;
  for(h=0; h < nirreps; h++) { ndocc += clsdpi[h]; nuocc += uoccpi[h]; }
  
  /* For Pitzer to QTS reordering array */
  openpi = init_int_array(nirreps);
  frdoccpi = init_int_array(nirreps);
  fruoccpi = init_int_array(nirreps);
  for(h=0; h < nirreps; h++) {
    openpi[h] = 0; 
    frdoccpi[h] = 0;
    fruoccpi[h] = 0;
  }
  
  num_ai = nuocc * ndocc;
  num_pq = nmo * (nmo +1)/2;
  num_pi = nmo * ndocc;
  
  /* Build the first and last lookup arrays */
  first = init_int_array(nirreps);
  last = init_int_array(nirreps);
  foffset = 0;
  loffset = orbspi[0] - 1;
  first[0] = foffset;
  last[0] = loffset;
  for(h=1; h < nirreps; h++) {
    foffset += orbspi[h-1];
    loffset += orbspi[h];
    first[h] = foffset;
    last[h] = loffset;
  }

  ofirst = init_int_array(nirreps);
  olast = init_int_array(nirreps);
  foffset = 0;
  loffset = clsdpi[0] - 1;
  ofirst[0] = foffset;
  olast[0] = loffset;
  for(h=1; h < nirreps; h++) {
    foffset += orbspi[h-1];
    loffset += uoccpi[h-1] + clsdpi[h];
    ofirst[h] = foffset;
    olast[h] = loffset;
  }

  vfirst = init_int_array(nirreps);
  vlast = init_int_array(nirreps);
  foffset = clsdpi[0];
  loffset = orbspi[0] - 1;
  vfirst[0] = foffset;
  vlast[0] = loffset;
  for(h=1; h < nirreps; h++) {
    foffset += uoccpi[h-1] + clsdpi[h];
    loffset += orbspi[h];
    vfirst[h] = foffset;
    vlast[h] = loffset;
  }
}

void cleanup(void)
{
  free(orbspi);
  free(clsdpi);
  free(uoccpi);
  free(evals);
  free(zvals);
  free(first);
  free(last);
  free(ofirst);
  free(olast);
  free(vfirst);
  free(vlast);
  free(ioff);
  free_block(scf);
  free_block(geom);
}

}} // namespace psi::cphf

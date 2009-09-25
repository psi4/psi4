/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libiwl/iwl.h>
#include <libchkpt/chkpt.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include <psifiles.h>
#include <physconst.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

#define _au2cgs 471.44353920

void transdip(void);
void transp(double sign);
void transL(double sign);

void rotational_strength(struct TD_Params *S)
{
  int i, j, k;
  int no, nv, nt;
  double lt_x, lt_y, lt_z;
  double rt_x, rt_y, rt_z;
  double rs_lx, rs_ly, rs_lz;
  double rs_rx, rs_ry, rs_rz;
  double rs_x, rs_y, rs_z;
  double rs;
  double conv;
  int nmo = moinfo.nmo;

  transdip();

  fprintf(outfile,"\n\tLength-Gauge Rotational Strength for %d%3s\n",S->root+1,
          moinfo.labels[S->irrep]);
  fprintf(outfile,"\t                              X    \t       Y    \t       Z\n");

  lt_x = lt_y = lt_z = 0.0;
  rt_x = rt_y = rt_z = 0.0;
  rs_lx = rs_ly = rs_lz = 0.0;
  rs_rx = rs_ry = rs_rz = 0.0;
  rs_x = rs_y = rs_z = 0.0;

  for(i=0; i < nmo; i++)
    for(j=0; j < nmo; j++) {
      lt_x += moinfo.ltd[i][j] * moinfo.dip[0][i][j];
      lt_y += moinfo.ltd[i][j] * moinfo.dip[1][i][j];
      lt_z += moinfo.ltd[i][j] * moinfo.dip[2][i][j];
    }

  transL(+1.0);

  for(i=0; i < nmo; i++)
    for(j=0; j < nmo; j++) {
      rt_x += moinfo.rtd[i][j] * moinfo.L[0][i][j];
      rt_y += moinfo.rtd[i][j] * moinfo.L[1][i][j];
      rt_z += moinfo.rtd[i][j] * moinfo.L[2][i][j];
    }

  rs_lx = lt_x * rt_x;
  rs_ly = lt_y * rt_y;
  rs_lz = lt_z * rt_z;

  fprintf(outfile,"\t<0|mu_e|n>              %11.8lf \t %11.8lf \t %11.8lf\n",
          lt_x,lt_y,lt_z);
  fprintf(outfile,"\t<n|mu_m|0>              %11.8lf \t %11.8lf \t %11.8lf\n",
          rt_x,rt_y,rt_z);

  // Complex Conjugate

  lt_x = lt_y = lt_z = 0.0;
  rt_x = rt_y = rt_z = 0.0;

  for(i=0; i < 3 ; i++)
    for(j=0; j < nmo; j++)
      for(k=0; k < nmo; k++)
        moinfo.L[i][j][k] = 0.0;

  transL(-1.0);

  for(i=0; i < nmo; i++)
    for(j=0; j < nmo; j++) {
      lt_x += moinfo.ltd[i][j] * moinfo.L[0][i][j];
      lt_y += moinfo.ltd[i][j] * moinfo.L[1][i][j];
      lt_z += moinfo.ltd[i][j] * moinfo.L[2][i][j];
    }

  for(i=0; i < nmo; i++)
    for(j=0; j < nmo; j++) {
      rt_x += moinfo.rtd[i][j] * moinfo.dip[0][i][j];
      rt_y += moinfo.rtd[i][j] * moinfo.dip[1][i][j];
      rt_z += moinfo.rtd[i][j] * moinfo.dip[2][i][j];
    }

  rs_rx = lt_x * rt_x;
  rs_ry = lt_y * rt_y;
  rs_rz = lt_z * rt_z;

  fprintf(outfile,"\t<0|mu_m|n>*             %11.8lf \t %11.8lf \t %11.8lf\n",
          lt_x,lt_y,lt_z);
  fprintf(outfile,"\t<n|mu_e|0>*             %11.8lf \t %11.8lf \t %11.8lf\n",
          rt_x,rt_y,rt_z);

  rs_x = 0.5 * ( rs_lx + rs_rx);
  rs_y = 0.5 * ( rs_ly + rs_ry);
  rs_z = 0.5 * ( rs_lz + rs_rz);

  rs = rs_x + rs_y + rs_z;
  S->RS_length = rs;

  fprintf(outfile,"\n");
  fprintf(outfile,"\tRotational Strength (au)                 %11.8lf\n",rs);
  fprintf(outfile,"\tRotational Strength (10^-40 esu^2 cm^2)  %11.8lf\n",rs*_au2cgs);
  fflush(outfile);

  fprintf(outfile,"\n\tVelocity-Gauge Rotational Strength for %d%3s\n",S->root+1,
          moinfo.labels[S->irrep]);
  fprintf(outfile,"\t                              X    \t       Y    \t       Z\n");

  lt_x = lt_y = lt_z = 0.0;
  rt_x = rt_y = rt_z = 0.0;
  rs_lx = rs_ly = rs_lz = 0.0;
  rs_rx = rs_ry = rs_rz = 0.0;
  rs_x = rs_y = rs_z = 0.0;

  transp(+1.0);

  for(i=0; i < nmo; i++)
    for(j=0; j < nmo; j++) {
      lt_x += moinfo.ltd[i][j] * moinfo.nabla[0][i][j];
      lt_y += moinfo.ltd[i][j] * moinfo.nabla[1][i][j];
      lt_z += moinfo.ltd[i][j] * moinfo.nabla[2][i][j];
    }

  transL(+1.0);

  for(i=0; i < nmo; i++)
    for(j=0; j < nmo; j++) {
      rt_x += moinfo.rtd[i][j] * moinfo.L[0][i][j];
      rt_y += moinfo.rtd[i][j] * moinfo.L[1][i][j];
      rt_z += moinfo.rtd[i][j] * moinfo.L[2][i][j];
    }

  rs_lx = lt_x * rt_x;
  rs_ly = lt_y * rt_y;
  rs_lz = lt_z * rt_z;

  fprintf(outfile,"\t<0|mu_e|n>              %11.8lf \t %11.8lf \t %11.8lf\n",
          lt_x,lt_y,lt_z);
  fprintf(outfile,"\t<n|mu_m|0>              %11.8lf \t %11.8lf \t %11.8lf\n",
          rt_x,rt_y,rt_z);

  // Complex Conjugate

  lt_x = lt_y = lt_z = 0.0;
  rt_x = rt_y = rt_z = 0.0;

  for(i=0; i < 3 ; i++)
    for(j=0; j < nmo; j++)
      for(k=0; k < nmo; k++) {
        moinfo.nabla[i][j][k] = 0.0;
        moinfo.L[i][j][k] = 0.0;
      }

  transL(-1.0);

  for(i=0; i < nmo; i++)
    for(j=0; j < nmo; j++) {
      lt_x += moinfo.ltd[i][j] * moinfo.L[0][i][j];
      lt_y += moinfo.ltd[i][j] * moinfo.L[1][i][j];
      lt_z += moinfo.ltd[i][j] * moinfo.L[2][i][j];
    }

  transp(-1.0);

  for(i=0; i < nmo; i++)
    for(j=0; j < nmo; j++) {
      rt_x += moinfo.rtd[i][j] * moinfo.nabla[0][i][j];
      rt_y += moinfo.rtd[i][j] * moinfo.nabla[1][i][j];
      rt_z += moinfo.rtd[i][j] * moinfo.nabla[2][i][j];
    }

  rs_rx = lt_x * rt_x;
  rs_ry = lt_y * rt_y;
  rs_rz = lt_z * rt_z;

  rs_x = 0.5 * ( rs_lx + rs_rx);
  rs_y = 0.5 * ( rs_ly + rs_ry);
  rs_z = 0.5 * ( rs_lz + rs_rz);

  fprintf(outfile,"\t<0|mu_m|n>*             %11.8lf \t %11.8lf \t %11.8lf\n",
          lt_x,lt_y,lt_z);
  fprintf(outfile,"\t<n|mu_e|0>*             %11.8lf \t %11.8lf \t %11.8lf\n",
          rt_x,rt_y,rt_z);

  rs_x = rs_x / S->cceom_energy;
  rs_y = rs_y / S->cceom_energy;
  rs_z = rs_z / S->cceom_energy;

  rs = rs_x + rs_y + rs_z;
  S->RS_velocity = rs;

  fprintf(outfile,"\n");
  fprintf(outfile,"\tRotational Strength (au)                 %11.8lf\n",rs);
  fprintf(outfile,"\tRotational Strength (10^-40 esu^2 cm^2)  %11.8lf\n",rs*_au2cgs);
  fflush(outfile);

  return;
}

}} // namespace psi::ccdensity

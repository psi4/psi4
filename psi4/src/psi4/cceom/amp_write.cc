/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \file
    \ingroup CCEOM
    \brief print out converged right-hand eigenvectors
        amp_write_RHF()
        amp_write_UHF()
        amp_write_ROHF()
      int namps = number of amplitudes to be printed
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "psi4/libdpd/dpd.h"
#include <vector>

#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

// minimum magnitude of amplitude to include in output
#define MIN_TO_SHOW (1e-5)

namespace psi { namespace cceom {

using std::vector;

class R1_amp {
  double value;
  int i, a, Gi, Ga;
  public:
    R1_amp() {
      i = a = Gi = Ga = 0;
      value = 0.0;
    }
    void zero(void) {
      i = a = Gi = Ga = 0;
      value = 0.0;
    }
  friend void amp_write_RHF(dpdfile2 *RIA, dpdbuf4 *RIjAb, int namps);
  friend void amp_write_UHF(dpdfile2 *RIA, dpdfile2 *Ria, dpdbuf4 *RIJAB,
    dpdbuf4 *Rijab, dpdbuf4 *RIjAb, int namps);
  friend void amp_write_ROHF(dpdfile2 *RIA, dpdfile2 *Ria, dpdbuf4 *RIJAB,
    dpdbuf4 *Rijab, dpdbuf4 *RIjAb, int namps);
  friend void get_largest_R1_amps(dpdfile2 *R1, int namps, vector<R1_amp> & R1_stack);
};

class R2_amp {
  double value;
  int i, j, a, b, Gi, Gj, Ga, Gb;
  public:
    R2_amp() {
      i = j = a = b = Gi = Gj = Ga = Gb = 0;
      value = 0.0;
    }
    void zero(void) {
      i = j = a = b = Gi = Gj = Ga = Gb = 0;
      value = 0.0;
    }

  friend void amp_write_RHF(dpdfile2 *RIA, dpdbuf4 *RIjAb, int namps);
  friend void amp_write_UHF(dpdfile2 *RIA, dpdfile2 *Ria, dpdbuf4 *RIJAB,
    dpdbuf4 *Rijab, dpdbuf4 *RIjAb, int namps);
  friend void amp_write_ROHF(dpdfile2 *RIA, dpdfile2 *Ria, dpdbuf4 *RIJAB,
    dpdbuf4 *Rijab, dpdbuf4 *RIjAb, int namps);
  friend void get_largest_R2_amps(dpdbuf4 *R2, int namps, vector<R2_amp> & R2_stack);
};

void get_largest_R1_amps(dpdfile2 *R1, int namps, vector<R1_amp> & R1s);
void get_largest_R2_amps(dpdbuf4 *R2, int namps, vector<R2_amp> & R2_stack);

void amp_write_RHF(dpdfile2 *RIA, dpdbuf4 *RIjAb, int namps) {
  int m, i, j, a, b, Gi, Gj, Ga, Gb, *frdocc, *clsdpi;
  char lbli[10], lblj[10], lbla[10], lblb[10];
  frdocc = moinfo.frdocc;
  clsdpi = moinfo.clsdpi;

  // Do RIA
  vector<R1_amp> R1_stack;
  get_largest_R1_amps(RIA, namps, R1_stack);

  outfile->Printf(" RIA (libdpd indices) : (cscf notation)\n");
  for(m=0; m < R1_stack.size(); m++) {
    if(fabs(R1_stack[m].value) > MIN_TO_SHOW) {
      Gi = R1_stack[m].Gi;
      Ga = R1_stack[m].Ga;
      i = frdocc[Gi] + R1_stack[m].i + 1;
      a = frdocc[Ga] + clsdpi[Ga] + R1_stack[m].a + 1;
      sprintf(lbli,"%d%s", i, moinfo.irr_labs_lowercase[Gi]);
      sprintf(lbla,"%d%s", a, moinfo.irr_labs_lowercase[Ga]);
      outfile->Printf( "       %3d > %3d      :    %6s > %6s : %15.10f\n",
        R1_stack[m].i, R1_stack[m].a, lbli, lbla, R1_stack[m].value);
    }
  }
  R1_stack.clear();

  // Do RIjAb
  vector<R2_amp> R2_stack;
  get_largest_R2_amps(RIjAb, namps, R2_stack);

  outfile->Printf(" RIjAb (libdpd indices)     : (cscf notation)\n");
  for(m=0; m < R2_stack.size(); m++) {
    if(fabs(R2_stack[m].value) > MIN_TO_SHOW) {
      Gi = R2_stack[m].Gi;
      Gj = R2_stack[m].Gj;
      Ga = R2_stack[m].Ga;
      Gb = R2_stack[m].Gb;

      i = frdocc[Gi] + R2_stack[m].i + 1;
      j = frdocc[Gj] + R2_stack[m].j + 1;
      a = frdocc[Ga] + clsdpi[Ga] + R2_stack[m].a + 1;
      b = frdocc[Gb] + clsdpi[Gb] + R2_stack[m].b + 1;

      sprintf(lbli,"%d%s", i, moinfo.irr_labs_lowercase[Gi]);
      sprintf(lblj,"%d%s", j, moinfo.irr_labs_lowercase[Gj]);
      sprintf(lbla,"%d%s", a, moinfo.irr_labs_lowercase[Ga]);
      sprintf(lblb,"%d%s", b, moinfo.irr_labs_lowercase[Gb]);
      outfile->Printf( "      %3d %3d > %3d %3d     : %6s %6s > %6s %6s : %15.10f\n",
        R2_stack[m].i, R2_stack[m].j, R2_stack[m].a, R2_stack[m].b,
          lbli, lblj, lbla, lblb, R2_stack[m].value);
    }
  }
  R2_stack.clear();
}

// ** amp_write_UHF()
void amp_write_UHF(dpdfile2 *RIA, dpdfile2 *Ria, dpdbuf4 *RIJAB,
    dpdbuf4 *Rijab, dpdbuf4 *RIjAb, int namps) {
  int m, i, j, a, b, Gi, Gj, Ga, Gb, *frdocc, *clsdpi, *openpi;
  char lbli[10], lblj[10], lbla[10], lblb[10];

  frdocc = moinfo.frdocc;
  clsdpi = moinfo.clsdpi;
  openpi = moinfo.openpi;

  // Do RIA
  vector<R1_amp> R1_stack;
  get_largest_R1_amps(RIA, namps, R1_stack);

  outfile->Printf(" RIA (libdpd indices) : (cscf notation)\n");
  for(m=0; m < R1_stack.size(); m++) {
    if(fabs(R1_stack[m].value) > MIN_TO_SHOW) {
      Gi = R1_stack[m].Gi;
      Ga = R1_stack[m].Ga;
      i = frdocc[Gi] + R1_stack[m].i + 1;
      a = frdocc[Ga] + clsdpi[Ga] + openpi[Ga] + R1_stack[m].a + 1;
      sprintf(lbli,"%d%s", i, moinfo.irr_labs_lowercase[Gi]);
      sprintf(lbla,"%d%s", a, moinfo.irr_labs_lowercase[Ga]);
      outfile->Printf( "       %3d > %3d      :    %6s > %6s : %15.10f\n",
        R1_stack[m].i, R1_stack[m].a, lbli, lbla, R1_stack[m].value);
    }
  }
  R1_stack.clear();

  // Do Ria
  get_largest_R1_amps(Ria, namps, R1_stack);

  outfile->Printf(" Ria (libdpd indices) : (cscf notation)\n");
  for(m=0; m < R1_stack.size(); m++) {
    if(fabs(R1_stack[m].value) > MIN_TO_SHOW) {
      Gi = R1_stack[m].Gi;
      Ga = R1_stack[m].Ga;
      i = frdocc[Gi] + R1_stack[m].i + 1;
      a = frdocc[Ga] + clsdpi[Ga] + R1_stack[m].a + 1;
      sprintf(lbli,"%d%s", i, moinfo.irr_labs_lowercase[Gi]);
      sprintf(lbla,"%d%s", a, moinfo.irr_labs_lowercase[Ga]);
      outfile->Printf( "       %3d > %3d      :    %6s > %6s : %15.10f\n",
        R1_stack[m].i, R1_stack[m].a, lbli, lbla, R1_stack[m].value);
    }
  }
  R1_stack.clear();

  // Do RIjAb
  vector<R2_amp> R2_stack;
  get_largest_R2_amps(RIjAb, namps, R2_stack);

  outfile->Printf(" RIjAb (libdpd indices)     : (cscf notation)\n");
  for(m=0; m < R2_stack.size(); m++) {
    if(fabs(R2_stack[m].value) > MIN_TO_SHOW) {
      Gi = R2_stack[m].Gi;
      Gj = R2_stack[m].Gj;
      Ga = R2_stack[m].Ga;
      Gb = R2_stack[m].Gb;

      i = frdocc[Gi] + R2_stack[m].i + 1;
      j = frdocc[Gj] + R2_stack[m].j + 1;
      a = frdocc[Ga] + clsdpi[Ga] + openpi[Ga] + R2_stack[m].a + 1;
      b = frdocc[Gb] + clsdpi[Gb] + R2_stack[m].b + 1;

      sprintf(lbli,"%d%s", i, moinfo.irr_labs_lowercase[Gi]);
      sprintf(lblj,"%d%s", j, moinfo.irr_labs_lowercase[Gj]);
      sprintf(lbla,"%d%s", a, moinfo.irr_labs_lowercase[Ga]);
      sprintf(lblb,"%d%s", b, moinfo.irr_labs_lowercase[Gb]);
      outfile->Printf( "      %3d %3d > %3d %3d     : %6s %6s > %6s %6s : %15.10f\n",
        R2_stack[m].i, R2_stack[m].j, R2_stack[m].a, R2_stack[m].b,
          lbli, lblj, lbla, lblb, R2_stack[m].value);
    }
  }
  R2_stack.clear();

  // Do RIJAB
  get_largest_R2_amps(RIJAB, namps, R2_stack);

  outfile->Printf(" RIJAB (libdpd indices)     : (cscf notation)\n");
  for(m=0; m < R2_stack.size(); m++) {
    if(fabs(R2_stack[m].value) > MIN_TO_SHOW) {
      Gi = R2_stack[m].Gi;
      Gj = R2_stack[m].Gj;
      Ga = R2_stack[m].Ga;
      Gb = R2_stack[m].Gb;

      i = frdocc[Gi] + R2_stack[m].i + 1;
      j = frdocc[Gj] + R2_stack[m].j + 1;
      a = frdocc[Ga] + clsdpi[Ga] + openpi[Ga] + R2_stack[m].a + 1;
      b = frdocc[Gb] + clsdpi[Gb] + openpi[Gb] + R2_stack[m].b + 1;

      sprintf(lbli,"%d%s", i, moinfo.irr_labs_lowercase[Gi]);
      sprintf(lblj,"%d%s", j, moinfo.irr_labs_lowercase[Gj]);
      sprintf(lbla,"%d%s", a, moinfo.irr_labs_lowercase[Ga]);
      sprintf(lblb,"%d%s", b, moinfo.irr_labs_lowercase[Gb]);
      outfile->Printf( "      %3d %3d > %3d %3d     : %6s %6s > %6s %6s : %15.10f\n",
        R2_stack[m].i, R2_stack[m].j, R2_stack[m].a, R2_stack[m].b,
          lbli, lblj, lbla, lblb, R2_stack[m].value);
    }
  }
  R2_stack.clear();

  // Do Rijab
  get_largest_R2_amps(Rijab, namps, R2_stack);

  outfile->Printf(" Rijab (libdpd indices)     : (cscf notation)\n");
  for(m=0; m < R2_stack.size(); m++) {
    if(fabs(R2_stack[m].value) > MIN_TO_SHOW) {
      Gi = R2_stack[m].Gi;
      Gj = R2_stack[m].Gj;
      Ga = R2_stack[m].Ga;
      Gb = R2_stack[m].Gb;

      i = frdocc[Gi] + R2_stack[m].i + 1;
      j = frdocc[Gj] + R2_stack[m].j + 1;
      a = frdocc[Ga] + clsdpi[Ga] + R2_stack[m].a + 1;
      b = frdocc[Gb] + clsdpi[Gb] + R2_stack[m].b + 1;

      sprintf(lbli,"%d%s", i, moinfo.irr_labs_lowercase[Gi]);
      sprintf(lblj,"%d%s", j, moinfo.irr_labs_lowercase[Gj]);
      sprintf(lbla,"%d%s", a, moinfo.irr_labs_lowercase[Ga]);
      sprintf(lblb,"%d%s", b, moinfo.irr_labs_lowercase[Gb]);
      outfile->Printf( "      %3d %3d > %3d %3d     : %6s %6s > %6s %6s : %15.10f\n",
        R2_stack[m].i, R2_stack[m].j, R2_stack[m].a, R2_stack[m].b,
          lbli, lblj, lbla, lblb, R2_stack[m].value);
    }
  }
  R2_stack.clear();
}

// ** amp_write_ROHF()
void amp_write_ROHF(dpdfile2 *RIA, dpdfile2 *Ria, dpdbuf4 *RIJAB,
    dpdbuf4 *Rijab, dpdbuf4 *RIjAb, int namps) {
  int m, i, j, a, b, Gi, Gj, Ga, Gb, *virtpi, *clsdpi, *openpi, *frdocc;
  char lbli[10], lblj[10], lbla[10], lblb[10];

  virtpi = moinfo.virtpi;
  clsdpi = moinfo.clsdpi;
  openpi = moinfo.openpi;
  frdocc = moinfo.frdocc;

  // Do RIA
  vector<R1_amp> R1_stack;
  get_largest_R1_amps(RIA, namps, R1_stack);

  outfile->Printf(" RIA (libdpd indices) : (cscf notation)\n");
  for(m=0; m < R1_stack.size(); m++) {
    if(fabs(R1_stack[m].value) > MIN_TO_SHOW) {
      Gi = R1_stack[m].Gi;
      Ga = R1_stack[m].Ga;
      i = frdocc[Gi] + R1_stack[m].i + 1;
      a = frdocc[Ga] + clsdpi[Ga] + openpi[Ga] + R1_stack[m].a + 1;
      sprintf(lbli,"%d%s", i, moinfo.irr_labs_lowercase[Gi]);
      sprintf(lbla,"%d%s", a, moinfo.irr_labs_lowercase[Ga]);
      outfile->Printf( "       %3d > %3d      :    %6s > %6s : %15.10f\n",
        R1_stack[m].i, R1_stack[m].a, lbli, lbla, R1_stack[m].value);
    }
  }
  R1_stack.clear();

  // Do Ria
  get_largest_R1_amps(Ria, namps, R1_stack);

  outfile->Printf(" Ria (libdpd indices) : (cscf notation)\n");
  for(m=0; m < R1_stack.size(); m++) {
    if(fabs(R1_stack[m].value) > MIN_TO_SHOW) {
      Gi = R1_stack[m].Gi;
      Ga = R1_stack[m].Ga;

      i = frdocc[Gi] + R1_stack[m].i + 1;
      if (R1_stack[m].a < (virtpi[Ga]-openpi[Ga]) )
        a = frdocc[Ga] + clsdpi[Ga] + openpi[Ga] + R1_stack[m].a + 1;
      else
        a = frdocc[Ga] + clsdpi[Ga] + (R1_stack[m].a - (virtpi[Ga]-openpi[Ga])) + 1;

      sprintf(lbli,"%d%s", i, moinfo.irr_labs_lowercase[Gi]);
      sprintf(lbla,"%d%s", a, moinfo.irr_labs_lowercase[Ga]);
      outfile->Printf( "       %3d > %3d      :    %6s > %6s : %15.10f\n",
        R1_stack[m].i, R1_stack[m].a, lbli, lbla, R1_stack[m].value);
    }
  }
  R1_stack.clear();

  // Do RIjAb
  vector<R2_amp> R2_stack;
  get_largest_R2_amps(RIjAb, namps, R2_stack);

  outfile->Printf(" RIjAb (libdpd indices)     : (cscf notation)\n");
  for(m=0; m < R2_stack.size(); m++) {
    if(fabs(R2_stack[m].value) > MIN_TO_SHOW) {
      Gi = R2_stack[m].Gi;
      Gj = R2_stack[m].Gj;
      Ga = R2_stack[m].Ga;
      Gb = R2_stack[m].Gb;

      i = frdocc[Gi] + R2_stack[m].i + 1;
      j = frdocc[Gj] + R2_stack[m].j + 1;
      a = frdocc[Ga] + clsdpi[Ga] + openpi[Ga] + R2_stack[m].a + 1;
      if (R2_stack[m].b < (virtpi[Gb]-openpi[Gb]) )
        b = frdocc[Gb] + clsdpi[Gb] + openpi[Gb] + R2_stack[m].b + 1;
      else
        b = frdocc[Gb] + clsdpi[Gb] + (R2_stack[m].b - (virtpi[Gb]-openpi[Gb])) + 1;

      sprintf(lbli,"%d%s", i, moinfo.irr_labs_lowercase[Gi]);
      sprintf(lblj,"%d%s", j, moinfo.irr_labs_lowercase[Gj]);
      sprintf(lbla,"%d%s", a, moinfo.irr_labs_lowercase[Ga]);
      sprintf(lblb,"%d%s", b, moinfo.irr_labs_lowercase[Gb]);
      outfile->Printf( "      %3d %3d > %3d %3d     : %6s %6s > %6s %6s : %15.10f\n",
        R2_stack[m].i, R2_stack[m].j, R2_stack[m].a, R2_stack[m].b,
          lbli, lblj, lbla, lblb, R2_stack[m].value);
    }
  }
  R2_stack.clear();

  // Do RIJAB
  get_largest_R2_amps(RIJAB, namps, R2_stack);

  outfile->Printf(" RIJAB (libdpd indices)     : (cscf notation)\n");
  for(m=0; m < R2_stack.size(); m++) {
    if(fabs(R2_stack[m].value) > MIN_TO_SHOW) {
      Gi = R2_stack[m].Gi;
      Gj = R2_stack[m].Gj;
      Ga = R2_stack[m].Ga;
      Gb = R2_stack[m].Gb;

      i = frdocc[Gi] + R2_stack[m].i + 1;
      j = frdocc[Gj] + R2_stack[m].j + 1;
      a = frdocc[Ga] + clsdpi[Ga] + openpi[Ga] + R2_stack[m].a + 1;
      b = frdocc[Gb] + clsdpi[Gb] + openpi[Gb] + R2_stack[m].b + 1;

      sprintf(lbli,"%d%s", i, moinfo.irr_labs_lowercase[Gi]);
      sprintf(lblj,"%d%s", j, moinfo.irr_labs_lowercase[Gj]);
      sprintf(lbla,"%d%s", a, moinfo.irr_labs_lowercase[Ga]);
      sprintf(lblb,"%d%s", b, moinfo.irr_labs_lowercase[Gb]);
      outfile->Printf( "      %3d %3d > %3d %3d     : %6s %6s > %6s %6s : %15.10f\n",
        R2_stack[m].i, R2_stack[m].j, R2_stack[m].a, R2_stack[m].b,
          lbli, lblj, lbla, lblb, R2_stack[m].value);
    }
  }
  R2_stack.clear();

  // Do Rijab
  get_largest_R2_amps(Rijab, namps, R2_stack);

  outfile->Printf(" Rijab (libdpd indices)     : (cscf notation)\n");
  for(m=0; m < R2_stack.size(); m++) {
    if(fabs(R2_stack[m].value) > MIN_TO_SHOW) {
      Gi = R2_stack[m].Gi;
      Gj = R2_stack[m].Gj;
      Ga = R2_stack[m].Ga;
      Gb = R2_stack[m].Gb;

      i = frdocc[Gi] + R2_stack[m].i + 1;
      j = frdocc[Gj] + R2_stack[m].j + 1;

      if (R2_stack[m].a < (virtpi[Ga]-openpi[Ga]) )
        a = frdocc[Ga] + clsdpi[Ga] + openpi[Ga] + R2_stack[m].a + 1;
      else
        a = frdocc[Ga] + clsdpi[Ga] + (R2_stack[m].a - (virtpi[Ga]-openpi[Ga])) + 1;

      if (R2_stack[m].b < (virtpi[Gb]-openpi[Gb]) )
        b = frdocc[Gb] + clsdpi[Gb] + openpi[Gb] + R2_stack[m].b + 1;
      else
        b = frdocc[Gb] + clsdpi[Gb] + (R2_stack[m].b - (virtpi[Gb]-openpi[Gb])) + 1;

      sprintf(lbli,"%d%s", i, moinfo.irr_labs_lowercase[Gi]);
      sprintf(lblj,"%d%s", j, moinfo.irr_labs_lowercase[Gj]);
      sprintf(lbla,"%d%s", a, moinfo.irr_labs_lowercase[Ga]);
      sprintf(lblb,"%d%s", b, moinfo.irr_labs_lowercase[Gb]);
      outfile->Printf( "      %3d %3d > %3d %3d     : %6s %6s > %6s %6s : %15.10f\n",
        R2_stack[m].i, R2_stack[m].j, R2_stack[m].a, R2_stack[m].b,
          lbli, lblj, lbla, lblb, R2_stack[m].value);
    }
  }
  R2_stack.clear();
}


/* builds std::vector of R1_amps of largest magnitude for dpd_file2 structures */
void get_largest_R1_amps(dpdfile2 *R1, int namps, vector<R1_amp> & R1_stack) {
  int h, m, i, a, Gia, nirreps;
  R1_amp one_R1;

  nirreps = R1->params->nirreps;
  Gia = R1->my_irrep;
  global_dpd_->file2_mat_init(R1);
  global_dpd_->file2_mat_rd(R1);
  R1_stack.push_back(one_R1); // create stack with 1 zero entry

  for(h=0; h < nirreps; h++) {
    one_R1.Gi = h;
    one_R1.Ga = h^Gia;
    for(i=0; i < R1->params->rowtot[h]; i++) {
      one_R1.i = i;
      for(a=0; a < R1->params->coltot[h^Gia]; a++) {
        one_R1.a = a;
        one_R1.value = R1->matrix[h][i][a];
        for(m=0; m < R1_stack.size() ; m++) {
          if((fabs(one_R1.value) - fabs(R1_stack[m].value)) > 1e-12) {
            R1_stack.insert(R1_stack.begin()+m, one_R1);
            if (R1_stack.size() > namps)
              R1_stack.erase(R1_stack.end()-1);
            break;
          }
        }
      }
    }
  }
  global_dpd_->file2_mat_close(R1);
  return;
}

/* builds std::vector of R2_amps of largest magnitude for dpd_file2 structures */
void get_largest_R2_amps(dpdbuf4 *R2, int namps, vector<R2_amp> & R2_stack) {
  int a, b, i, j, h, ij, ab, m, nirreps, Gijab;
  R2_amp one_R2;

  nirreps = R2->params->nirreps;
  Gijab = R2->file.my_irrep;
  R2_stack.push_back(one_R2); // create stack with 1 zero entry

  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(R2, h);
    global_dpd_->buf4_mat_irrep_rd(R2, h);

    for(ij=0; ij < R2->params->rowtot[h]; ij++) {
      i = R2->params->roworb[h][ij][0];
      j = R2->params->roworb[h][ij][1];
      one_R2.Gi = R2->params->psym[i];
      one_R2.Gj = R2->params->qsym[j];
      one_R2.i = i - R2->params->poff[one_R2.Gi];
      one_R2.j = j - R2->params->qoff[one_R2.Gj];

      for(ab=0; ab < R2->params->coltot[h^Gijab]; ab++) {
        a = R2->params->colorb[h^Gijab][ab][0];
        b = R2->params->colorb[h^Gijab][ab][1];
        one_R2.Ga = R2->params->rsym[a];
        one_R2.Gb = R2->params->ssym[b];
        one_R2.a = a - R2->params->roff[one_R2.Ga];
        one_R2.b = b - R2->params->soff[one_R2.Gb];

        one_R2.value = R2->matrix[h][ij][ab];
        for(m=0; m < R2_stack.size() ; m++) {
          if((fabs(one_R2.value) - fabs(R2_stack[m].value)) > 1e-12) {
            R2_stack.insert(R2_stack.begin()+m, one_R2);
            if (R2_stack.size() > namps)
              R2_stack.erase(R2_stack.end()-1);
            break;
          }
        }
      }
    }
    global_dpd_->buf4_mat_irrep_close(R2, h);
  }
}

}} // namespace psi::cceom

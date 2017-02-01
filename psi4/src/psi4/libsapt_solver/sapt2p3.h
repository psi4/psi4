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

#ifndef SAPT2p3_H
#define SAPT2p3_H

#include "sapt2p.h"

namespace psi { namespace sapt {

class SAPT2p3 : public SAPT2p {
private:
  virtual void print_header();
  virtual void print_results();

  bool third_order_;

protected:
  double e_elst13_;
  double e_ind30_;
  double e_exch_ind30_;
  double e_ind30r_;
  double e_exch_ind30r_;
  double e_ind_disp30_;
  double e_exch_ind_disp30_;
  double e_disp30_;
  double e_exch_disp30_;
  double e_sapt2pp3_;
  double e_sapt2p3_;
  double e_sapt2pp3_ccd_;
  double e_sapt2p3_ccd_;

  void Y3(int, const char *, const char *, const char *, int, const char *,
    const char *, const char *, const char *, const char *, const char *,
    int, int, int, double *, int, const char *);
  void Y3_1(double **, int, const char *, const char *, int, const char *,
    int, int, int);
  void Y3_2(double **, int, const char *, const char *, int, const char *,
    int, int, int);
  void Y3_3(double **, int, const char *, const char *, const char *, int,
    const char *, int, int, int);
  void Y3_4(double **, int, const char *, const char *, const char *, int,
    const char *, int, int, int);

  double elst130(double **, double **, double **, int, const char *,
    const char *, const char *, int, int, int);

  void ind30_amps(int, const char *, int, const char *, double **, double **,
    double **, double **, int, int, double *, int, int, double *, int,
    const char *);

  void inddisp30_amps();
  void inddisp30_ov(int, const char *, const char *, int, const char *,
    int, int, int, double *, int, const char *);
  void inddisp30_ovov();

  double disp30_1(int, const char *, int, const char *, int, const char *,
    int, int, int, int, int, int);
  double disp30_2(int, const char *, int, const char *, const char *, int,
    const char *, const char *, int, int, int, int, int, int);

  void disp30_amps(int, const char *, int, const char *, const char *,
    int, const char *, const char *, int, int, int, double *,
    int, int, int, double *, int, const char *);

  double exch_ind30_1(double **, double **);
  double exch_ind30_2(double **);
  double exch_ind30_3(double **);

  double exch_ind_disp30_21(double **);
  double exch_ind_disp30_12(double **);

  double exch_disp30_20();
  double exch_disp30_02();
  double exch_disp30_22();

  double ind30r_1(double **, double **, double **, double **, int,
    const char *, const char *, const char *, int, const char *, int,
    int, int, int);

public:
  SAPT2p3(SharedWavefunction Dimer, SharedWavefunction MonomerA,
          SharedWavefunction MonomerB, Options& options,
          std::shared_ptr<PSIO>psio);
  virtual ~SAPT2p3();

  virtual double compute_energy();

  virtual void amplitudes();

  void elst13();
  void ind30();
  void ind30r();
  void exch_ind30();
  void ind_disp30();
  void exch_ind_disp30();
  void disp30();
  void exch_disp30();

};

}}

#endif

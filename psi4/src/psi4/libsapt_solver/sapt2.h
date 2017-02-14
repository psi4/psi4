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

#ifndef SAPT2_H
#define SAPT2_H

#include "sapt.h"

namespace psi { namespace sapt {

class SAPT2 : public SAPT {
private:
  virtual void print_header();
  virtual void print_results();

protected:
  size_t no_nvirA_;
  size_t no_nvirB_;

  double *no_evalsA_;
  double *no_evalsB_;

  double **no_CA_;
  double **no_CB_;

// These ints should not overflow below 46000 basis functions
// They contain indices up to nso_ * (nso_ + 1) / 2
  int *ioff_;
  int *index2i_;
  int *index2j_;

  int maxiter_;
  double e_conv_;
  double d_conv_;

  bool nat_orbs_t3_;
  bool nat_orbs_t2_;
  bool nat_orbs_v4_;
  double occ_cutoff_;

  double e_elst10_;
  double e_elst12_;
  double e_exch10_;
  double e_exch10_s2_;
  double e_exch11_;
  double e_exch12_;
  double e_ind20_;
  double e_ind22_;
  double e_exch_ind20_;
  double e_exch_ind22_;
  double e_disp20_;
  double e_no_disp20_;
  double e_exch_disp20_;
  double e_sapt0_;
  double e_sapt2_;

  double **wBAA_;
  double **wBAR_;
  double **wBRR_;

  double **wABB_;
  double **wABS_;
  double **wASS_;

  double** get_AA_ints(const int, int=0, int=0);
  double** get_diag_AA_ints(const int);
  double** get_AR_ints(const int, int=0);
  double** get_RR_ints(const int);
  double** get_BB_ints(const int, int=0, int=0);
  double** get_diag_BB_ints(const int);
  double** get_BS_ints(const int, int=0);
  double** get_SS_ints(const int);
  double** get_AB_ints(const int, int=0, int=0);
  double** get_AS_ints(const int, int=0);
  double** get_RB_ints(const int, int=0);

  void df_integrals();
  void w_integrals();

  double **get_DF_ints(int, const char *, int, int, int, int);
  double **get_DF_ints_nongimp(int, const char *, int, int, int, int);
  void antisym(double *, int, int);
  void antisym(double **, int, int);

  void cphf_solver(double**, double **, double *, int, const char *,
    const char *, const char *, int, int);

  void exch_ind20rA_B();
  void exch_ind20rB_A();

  void tOVOV(int, const char *, int, int, int, double *, int, const char *,
    int, int, int, double *, int, const char *);
  void pOOpVV(int, const char *, const char *, int, int, int, const char *,
    const char *);
  void theta(int, const char *, const char, bool, int, int, int, int,
    const char *, int, const char *);

  void Y2(int, const char *, const char *, const char *, int, const char *,
    const char *, const char *, int, int, int, double *, int, const char *,
    const char *);
  void Y2_1(double **, int, const char *, const char *, int, const char *,
    int, int, int);
  void Y2_2(double **, int, const char *, const char *, int, const char *,
    int, int, int);
  void Y2_3(double **, int, const char *, const char *, int, const char *,
    int, int, int);

  void t2OVOV(int, const char *, const char *, int, const char *,
    const char *, const char *, int, int, int, double *, int, const char *);
  void t2OVOV(int, const char *, const char *, const char *, int,
    const char *, const char *, const char *, const char *, int, int, int,
    int, double *, double **, int, const char *);

  void OVOpVp_to_OVpOpV(double *, int, int);
  void ijkl_to_ikjl(double *, int, int, int, int);
  void symmetrize(double *, int, int);

  void natural_orbitalify(int, const char *, double *evals, int, int, int,
    const char);
  void natural_orbitalify_df_ints();

  double elst120(double **, double **, double **, int, const char *,
    const char *, const char *, int, int, int);

  double exch110(int, const char *);
  double exch101(int, const char *);
  double exch111();
  double exch120_k2f();
  double exch102_k2f();
  double exch120_k11u_1();
  double exch102_k11u_1();
  double exch120_k11u_2();
  double exch102_k11u_2();
  double exch120_k11u_3();
  double exch102_k11u_3();
  double exch120_k11u_4();
  double exch102_k11u_4();
  double exch120_k11u_5();
  double exch102_k11u_5();
  double exch120_k11u_6();
  double exch102_k11u_6();

  double ind220();
  double ind202();
  double ind220_1(int, const char *, const char *, const char *, int,
    const char *, double **, double **, double **, int, int, int, double *);
  double ind220_2(int, const char *, double **, double **, double **, int,
    int, int);
  double ind220_3(int, const char *, const char *, double **, double **,
    int, int, int);
  double ind220_4(int, const char *, int, const char *, double **, int,
    int, int);
  double ind220_5(int, const char *, double **, int, int, int, double *);
  double ind220_6(int, const char *, const char *, const char *, int,
    const char *, double **, int, int, int);
  double ind220_7(int, const char *, const char *, const char *, int,
    const char *, int, const char *, const char *, const char *, double **,
    int, int, int, int, int, int);

public:
  SAPT2(SharedWavefunction Dimer, SharedWavefunction MonomerA,
        SharedWavefunction MonomerB, Options& options,
        std::shared_ptr<PSIO>psio);
  virtual ~SAPT2();

  virtual double compute_energy();

  virtual void amplitudes();

  void elst10();
  void exch10_s2();
  void exch10();
  void ind20r();
  void exch_ind20r();
  void disp20();
  void exch_disp20();
  void elst12();
  void exch11();
  void exch12();
  void ind22();

};

}}

#endif

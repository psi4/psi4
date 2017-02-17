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

#ifndef SAPT_H
#define SAPT_H


#include "psi4/psifiles.h"

#ifdef _OPENMP
  #include <omp.h>
#endif

#ifdef USING_LAPACK_MKL
  #include <mkl.h>
#endif

#define INDEX(i,j) ((i>=j) ? (ioff_[i] + j) : (ioff_[j] + i))


#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsio/aiohandler.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"
#include "psi4/lib3index/3index.h"
#include "psi4/libmints/wavefunction.h"


namespace psi { namespace sapt {

class SAPT : public Wavefunction {

private:
  void initialize(SharedWavefunction MonomerA, SharedWavefunction MonomerB);
  void get_denom();

protected:
  std::shared_ptr<BasisSet> ribasis_;
  std::shared_ptr<BasisSet> elstbasis_;
  std::shared_ptr<BasisSet> zero_;

  size_t nsoA_;
  size_t nmoA_;
  size_t nsoB_;
  size_t nmoB_;
  size_t ndf_;
  size_t noccA_;
  size_t foccA_;
  size_t aoccA_;
  size_t noccB_;
  size_t foccB_;
  size_t aoccB_;
  size_t nvirA_;
  size_t nvirB_;
  int NA_;
  int NB_;
  int natomsA_;
  int natomsB_;

  bool elst_basis_;

  long int mem_;

  // Alpha exponent for exchange scaling
  double exch_scale_alpha_;

  double enuc_;
  double eHF_;
  double schwarz_;

  double *evalsA_;
  double *evalsB_;
  double *diagAA_;
  double *diagBB_;

  double **CA_;
  double **CB_;
  double **CHFA_;
  double **CHFB_;
  double **sAB_;
  double **vABB_;
  double **vBAA_;
  double **vAAB_;
  double **vBAB_;

  std::shared_ptr<SAPTDenominator> denom_;

  size_t nvec_;

  double **dAR_;
  double **dBS_;

  void zero_disk(int, const char *, int, int);

public:
  SAPT(SharedWavefunction Dimer, SharedWavefunction MonomerA,
       SharedWavefunction MonomerB, Options& options,
       std::shared_ptr<PSIO> psio);
  virtual ~SAPT();

  virtual double compute_energy()=0;
};

class CPHFDIIS {

private:
  int max_diis_vecs_;
  size_t vec_length_;

  int curr_vec_;
  int num_vecs_;

  double **t_vecs_;
  double **err_vecs_;

protected:

public:
  CPHFDIIS(int, int);
  ~CPHFDIIS();

  void store_vectors(double *, double *);
  void get_new_vector(double *);
};

}}

#endif

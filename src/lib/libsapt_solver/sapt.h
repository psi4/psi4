/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

#ifndef SAPT_H
#define SAPT_H

#include <psiconfig.h>
#include <psifiles.h>

#ifdef _OPENMP
  #include <omp.h>
#endif

#ifdef HAVE_MKL
  #include <mkl.h>
#endif

#define INDEX(i,j) ((i>=j) ? (ioff_[i] + j) : (ioff_[j] + i))

#include <libmints/mints.h>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>
#include <libpsio/aiohandler.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <lib3index/3index.h>

namespace boost {

template<class T> class shared_ptr;

}

namespace psi { namespace sapt {

class SAPT : public Wavefunction {

private:
  void initialize();
  void get_denom();

protected:
  boost::shared_ptr<BasisSet> ribasis_;
  boost::shared_ptr<BasisSet> elstbasis_;
  boost::shared_ptr<BasisSet> zero_;

  int nso_;
  int nmo_;
  int nsoA_;
  int nmoA_;
  int nsoB_;
  int nmoB_;
  int ndf_;
  int noccA_;
  int foccA_;
  int aoccA_;
  int noccB_;
  int foccB_;
  int aoccB_;
  int nvirA_;
  int nvirB_;
  int NA_;
  int NB_;
  int natomsA_;
  int natomsB_;

  bool elst_basis_;

  long int mem_;

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

  boost::shared_ptr<SAPTDenominator> denom_;

  int nvec_;

  double **dAR_;
  double **dBS_;

  void zero_disk(int, const char *, int, int);

public:
  virtual bool same_a_b_orbs() const { return reference_wavefunction_->same_a_b_orbs(); }
  virtual bool same_a_b_dens() const { return reference_wavefunction_->same_a_b_dens(); }

  SAPT(Options& options, boost::shared_ptr<PSIO> psio,
    boost::shared_ptr<Chkpt> chkpt);
  virtual ~SAPT();

  virtual double compute_energy()=0;
};

class CPHFDIIS {

private:
  int max_diis_vecs_;
  int vec_length_;

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

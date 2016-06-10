/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
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

#ifndef _psi_libfock_link_h
#define _psi_libfock_link_h

#include <cfloat>

#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <sstream>

#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libiwl/iwl.h>
#include <libqt/qt.h>

#include <libfock/apps.h>
#include <libfock/jk.h>

#include <libmints/typedefs.h>
#include <libmints/sieve.h>


namespace psi {

  class LinK {
    
  protected:
    
    boost::shared_ptr<BasisSet> basis_;
    
    bool do_J_;
    bool do_K_;
    bool do_wK_;
    
    //std::vector<SharedMatrix >& D_;
    std::vector<SharedMatrix> D_;
    
    // the coulomb and exchange matrices, to be filled in and passed back out
    std::vector<SharedMatrix> K_;
    //std::vector<SharedMatrix> wK_;
    
    // ERIs
    std::vector<boost::shared_ptr<TwoBodyAOInt> > eri_;
    /// Integral factory (must be retained for Spherical Transforms)
    boost::shared_ptr<IntegralFactory> factory_;
    /// ERI Sieve
    boost::shared_ptr<ERISieve> sieve_;
    
    
    // for each mu, we have a vector of pairs (d_munu * numax, nu_ind)
    std::vector<std::vector<std::pair<double, int> > > nu_for_mu_indices_;
    
    // shell_to_shell[mu_ind] is pairs Q_munu, nu_ind for each shell mu_ind,
    // sorted by decreasing values of Q
    std::vector<std::vector<std::pair<double, int> > > shell_to_shell_;
    
    // (mu_max | mu_max)^{1/2} for each mu
    std::vector<double> shell_max_q_sqr_;
    
    double shell_pair_threshold_sqr_;
    double integral_threshold_sqr_;
    
    long long int num_integrals_;
    long long int total_num_integrals_;
    
    
    ///////////////////// functions ////////////////////////////
    
    
    void FindSignificantNuforMu_(int mu_ind);
    // sort these lists (for each mu) when finished
    
    // for each bra shell pair, form list of ket shell pairs
    void FormSignificantShellPairList_();
    
    void ContractIntegrals_(int mu_ind, int lambda_ind,
                            std::vector<std::pair<int,int> >& ml_integrals);

    
  public:
    
    LinK(boost::shared_ptr<BasisSet> basis_in,
         std::vector<SharedMatrix>& density_in);
    
    ~LinK();
    
    void Compute();
    
    // Assuming this always gets called before Compute()
    void Update(const std::vector<SharedMatrix>& D_new);
    
    std::vector<SharedMatrix>& J();
    
    std::vector<SharedMatrix>& K();
    
    void set_do_J(bool do_it);
    
    void set_do_K(bool do_it);
    
    void print_header() const;

    
    
  }; // class

} // namespace

#endif
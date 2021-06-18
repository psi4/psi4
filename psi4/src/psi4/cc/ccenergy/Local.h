/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2021 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \file
    \ingroup CCENERGY
    \brief Local_cc class definition for PNO simulation
*/

#ifndef _psi_src_bin_ccenergy_local_h
#define _psi_src_bin_ccenergy_local_h

#include <string>
#include "psi4/libmints/wavefunction.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/local.h"
#include "psi4/libmints/basisset.h"

namespace psi {

class Local_cc {
    public:
        Local_cc();
        ~Local_cc();
        void local_init();
        void init_pno();
        void init_pnopp(const double omega, bool combined=false);
        void init_cpnopp(const double omega);
        void init_filter_T2();
        void local_filter_T1(dpdfile2 *T1);
        void local_filter_T2(dpdbuf4 *T2);
        void local_done();

        std::string method;
        std::string weakp;
        std::string pairdef;
        std::string pert;
        int filter_singles;
        int nocc;
        int nvir;
        double cutoff;
        double unpert_cutoff;
        double weak_pair_cutoff;
        double weak_pair_energy;

    private:
        void get_matvec(dpdbuf4 *buf_obj, std::vector<SharedMatrix> *matvec);
        void get_semicanonical_transforms(std::vector<SharedMatrix> Q);
        std::vector<SharedMatrix> build_PNO_lists(double cutoff, std::vector<SharedMatrix> D);
        std::vector<SharedMatrix> build_cPNO_lists(double cutoff, std::vector<SharedMatrix> D);
        void mp2_pair_energy();
        void pair_perturbation();
        
        int npairs;

};

} //namespace psi

#endif  // _psi_src_bin_ccenergy_local_h

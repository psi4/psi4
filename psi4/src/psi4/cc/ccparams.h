/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
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

#pragma once

#include <string>

namespace psi {

class Options;

namespace cc {

enum class Reference { RHF, ROHF, UHF };
enum class DerivativeType { NONE = 0, FIRST = 1, RESPONSE = 3 };
enum class CacheType { LRU, LOW };

/*! Holds coupled cluster calculation parameters
 *  FIXME
 *  memory, just_energy, just_residuals have not been ported from Params.h
 */
struct CCParams final {
    CCParams() {}
    CCParams(Options options);

    /*! Which coupled wavefunction? */
    std::string wfn;
    /*! Use new triples? */
    bool newtrips;
    /*! Use Brueckner? */
    bool brueckner;
    /*! Use density-fitting? */
    bool df;
    /*! Semicanonical calculation? */
    bool semicanonical;
    /*! Reference determinant */
    Reference ref;
    /*! Analyze T2 amplitudes? */
    bool analyze;
    /*! Derivative type */
    DerivativeType dertype;
    /*! Print level */
    int print;
    /*! Maximum number of CC iterations */
    int maxiter;
    /*! Convergence threshold on the residual */
    double convergence;
    /*! Convergence threshold on the energy */
    double e_convergence;
    /*! Is this a restart? */
    bool restart;
    std::string aobasis;
    int cachelevel;
    CacheType cachetype;
    /*! Number of threads.
     * TODO This is only used in cc3.cc and can probably be inferred in some other way.
     */
    int nthreads;
    bool diis;
    bool t2_coupled;
    std::string prop;
    std::string abcd;
    /*! Number of amplitudes to print */
    int num_amps;
    double bconv;

    bool print_mp2_amps;
    bool print_pair_energies;
    bool spinadapt_energies;
    bool t3_Ws_incore;

    bool scsn;
    bool scs;
    bool scscc;
    double scsmp2_scale_os;
    double scsmp2_scale_ss;
    double scscc_scale_os;
    double scscc_scale_ss;

    bool local;
    /*! @{ Local coupled cluster calculation set up options. Previously in Local */
    double local_cutoff;
    std::string local_method;
    std::string local_weakp;
    double local_cphf_cutoff;
    bool local_freeze_core;
    std::string local_pairdef;
};

void print_parameters(const CCParams &params, size_t memory);

}  // namespace cc
}  // namespace psi

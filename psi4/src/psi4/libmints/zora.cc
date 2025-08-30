/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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

#include "psi4/physconst.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libmints/zora.h"
#include "psi4/libqt/qt.h"

#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/basisset.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libfock/points.h"
#include "psi4/libfock/cubature.h"

#ifdef _OPENMP
#include <omp.h>
#include "psi4/libpsi4util/process.h"
#endif

#include <string>
#include <cmath>

namespace psi {

ZORA::ZORA(std::shared_ptr<Molecule> molecule, std::shared_ptr<BasisSet> basis, Options& options)
    : options_(options), molecule_(molecule), primary_(basis) {}

ZORA::~ZORA() {}

void ZORA::setup() {
    outfile->Printf("\n        ____  ____  ___  ___\n");
    outfile->Printf("       /_  / / __ \\/ _ \\/ _ |\n");
    outfile->Printf("        / /_/ /_/ / , _/ __ |\n");
    outfile->Printf("       /___/\\____/_/|_/_/ |_|\n");
    outfile->Printf("        by Nathan Gillispie\n");
    outfile->Printf("      ========================\n");

    timer_on("Make Grid");
    // Initialize grid with options
    std::map<std::string, std::string> grid_str_options = {
        {"DFT_PRUNING_SCHEME", options_.get_str("ZORA_PRUNING_SCHEME")},
        {"DFT_RADIAL_SCHEME", "TREUTLER"},
        {"DFT_NUCLEAR_SCHEME", "TREUTLER"},
        {"DFT_GRID_NAME", ""},
        {"DFT_BLOCK_SCHEME", "OCTREE"},
    };

    std::map<std::string, int> grid_int_options = {
        {"DFT_SPHERICAL_POINTS", options_.get_int("ZORA_SPHERICAL_POINTS")},
        {"DFT_RADIAL_POINTS", options_.get_int("ZORA_RADIAL_POINTS")},
        {"DFT_BLOCK_MIN_POINTS", 100},
        {"DFT_BLOCK_MAX_POINTS", 256},
    };

    std::map<std::string, double> grid_double_options = {
        {"DFT_BASIS_TOLERANCE", options_.get_double("ZORA_BASIS_TOLERANCE")},
        {"DFT_BS_RADIUS_ALPHA", 1.0},
        {"DFT_PRUNING_ALPHA", 1.0},
        {"DFT_BLOCK_MAX_RADIUS", 3.0},
        {"DFT_WEIGHTS_TOLERANCE", 1e-15},
    };

    grid_ = std::make_shared<DFTGrid>(primary_->molecule(), primary_, grid_int_options, grid_str_options,
                                      grid_double_options, options_);
    timer_off("Make Grid");

    outfile->Printf("\n  ==> ZORA Grid Details <==\n");
    outfile->Printf("    Basis: %s\n", primary_->name().c_str());
    outfile->Printf("    Total number of grid points: %d \n", grid_->npoints());
    outfile->Printf("    Total number of blocks: %d \n", grid_->blocks().size());
}

void ZORA::compute(SharedMatrix T_SR) {
    timer_on("ZORA");

    setup();

    int nblocks = grid_->blocks().size();
    int max_points = grid_->max_points();
    int max_funcs = grid_->max_functions();

    int nthreads = 1;
#ifdef _OPENMP
    nthreads = Process::environment.get_n_threads();
#endif
    outfile->Printf("    Using %d thread(s)\n\n", nthreads);

    // Basis function computer on each thread
    std::vector<std::shared_ptr<BasisFunctions>> pworkers;
    for (int i = 0; i < nthreads; i++) {
        auto p_tmp = std::make_shared<BasisFunctions>(primary_, max_points, max_funcs);
        p_tmp->set_deriv(1);
        pworkers.push_back(p_tmp);
    }

    timer_on("Effective Potential");
    veff_ = std::make_shared<std::map<int, SharedVector>>();
    if (options_.get_bool("ZORA_NR_DEBUG")) {
        compute_debug_veff();
    } else {
        compute_veff();
    }
    timer_off("Effective Potential");

    timer_on("Scalar Relativistic Kinetic");
    compute_TSR(pworkers, T_SR);
    timer_off("Scalar Relativistic Kinetic");

    timer_off("ZORA");
}

void ZORA::compute_debug_veff() {
    for (const auto &block : grid_->blocks()) {
        int index = block->index();
        int npoints = block->npoints();

        auto veff_block = std::make_shared<Vector>(npoints);
        veff_block->zero();
        veff_->insert({index, veff_block});
    }
}

void ZORA::compute_veff() {
    int natoms = molecule_->natom();

    int nthreads = 1;
#ifdef _OPENMP
    nthreads = Process::environment.get_n_threads();
#endif

#pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (const auto &block : grid_->blocks()) {
        int index = block->index();
        int npoints = block->npoints();

        double* x = block->x();
        double* y = block->y();
        double* z = block->z();

        auto veff_block = std::make_shared<Vector>(npoints);
        veff_block->zero(); // just in case

        for (int atom_idx = 0; atom_idx < natoms; atom_idx++) {
            int Z = molecule_->Z(atom_idx);
            if (Z > 104) throw PSIEXCEPTION("Molecule contains atom with atomic number above 104.");
            auto pos_a = molecule_->xyz(atom_idx);

            const double* coef_a  = &coeffs[c_aIndex[Z-1]];
            const double* alpha_a = &alphas[c_aIndex[Z-1]];
            int nc_a = c_aIndex[Z] - c_aIndex[Z-1];

            //einsums("i,ip->p", 𝕔[i], erf(α[i]⊗ r[p]))/r[p]
            for (int p = 0; p < npoints; p++) {
                double dist = std::hypot(pos_a[0]-x[p], pos_a[1]-y[p], pos_a[2]-z[p]);
                double outer = 0;
                for (int i = 0; i < nc_a; i++) {
                    outer += std::erf(dist * alpha_a[i]) * coef_a[i];
                }
                outer /= dist;
                outer -= Z/dist;
                veff_block->add(p, outer);
            }
        }
#pragma omp critical
        {
            veff_->insert({index, veff_block});
        }
    }
}


//Scalar Relativistic Kinetic Energy Matrix
void ZORA::compute_TSR(std::vector<std::shared_ptr<BasisFunctions>> pworkers, SharedMatrix &T_SR) {
    // Speed of light in atomic units
    double C = pc_c_au;

    int nthreads = 1;
#ifdef _OPENMP
    nthreads = Process::environment.get_n_threads();
#endif

#pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (const auto &block : grid_->blocks()) {
        const auto &bf_map = block->functions_local_to_global();
        auto local_nbf = bf_map.size();
        int npoints = block->npoints();

        auto veff_block = veff_->at(block->index());

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif
        pworkers[thread]->compute_functions(block);
        auto phi_x = pworkers[thread]->basis_value("PHI_X");
        auto phi_y = pworkers[thread]->basis_value("PHI_Y");
        auto phi_z = pworkers[thread]->basis_value("PHI_Z");


        // Preprocess kernel c²/(2c²-veff) * weight
        double* w = block->w();
        double kernel[npoints];
        for (int p = 0; p < npoints; p++) {
            kernel[p] = C *C /(2.*C *C - veff_block->get(p)) * w[p];
        }

        auto tmp = std::make_shared<Matrix>(T_SR->ncol(), T_SR->nrow());

        // Compute kinetic integral using kernel above.
        // T_SR --> non-relativistic T when veff --> 0.
        // bf_map is needed because the basis funcions differ from that of a given block
        for (int l_mu = 0; l_mu < local_nbf; l_mu++) {
            int mu = bf_map[l_mu];
            for (int l_nu = l_mu; l_nu < local_nbf; l_nu++) {
                int nu = bf_map[l_nu];
                for (int p = 0; p < npoints; p++) {
                    // ∇²φ(r)*kernel
                    tmp->add(mu, nu,
                         (phi_x->get(p, l_mu) * phi_x->get(p, l_nu) +
                          phi_y->get(p, l_mu) * phi_y->get(p, l_nu) +
                          phi_z->get(p, l_mu) * phi_z->get(p, l_nu)) * kernel[p]);
                }
            }
        }

        // Lock the T_SR matrix before adding
#pragma omp critical
        {
            T_SR->add(tmp);
        }
    }

    T_SR->copy_upper_to_lower();
}

}  // namespace psi

/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
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

#include "jk.h"
#include "SplitJK.h"
#include "psi4/libqt/qt.h"
#include "psi4/libfock/cubature.h"
#include "psi4/libfock/points.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/electrostatic.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/integral.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/lib3index/dftensor.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#include <unordered_set>
#include <vector>
#include <map>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace psi;

namespace psi {

// convers a Psi4::Molecule object to a GauXC::Molecule object
GauXC::Molecule snLinK::psi4_to_gauxc_molecule(std::shared_ptr<Molecule> psi4_molecule) {
    GauXC::Molecule gauxc_molecule;
 
    // TODO: Check if Bohr/Angstrom conversion is needed
    for (size_t iatom = 0; iatom != psi4_molecule->natom(); ++iatom) {
        gauxc_molecule.emplace_back(
            GauXC::AtomicNumber(psi4_molecule->true_atomic_number(iatom)), 
            psi4_molecule->x(iatom),
            psi4_molecule->y(iatom),
            psi4_molecule->z(iatom)
        );
    }

    return gauxc_molecule;
}

// converts a Psi4::Matrix to an Eigen::MatrixXd
Eigen::MatrixXd snLinK::psi4_to_eigen_matrix(SharedMatrix psi4_matrix) {
    if (psi4_matrix->nirrep() != 1) {
        throw PSIEXCEPTION("Psi4:: Matrix must be in C1 symmetry to be transformed into Eigen::MatrixXd!");
    }

    Eigen::MatrixXd eigen_matrix(psi4_matrix->nrow(), psi4_matrix->ncol());

    auto psi4_matrix_p = psi4_matrix->pointer(); 
    for (size_t irow = 0; irow != psi4_matrix->nrow(); ++irow) {
        for (size_t icol = 0; icol != psi4_matrix->ncol(); ++icol) {
            eigen_matrix(irow, icol) = psi4_matrix_p[irow][icol]; //TODO: check that Psi4::Matrix is column-major
        }
    }

    return eigen_matrix;
}

snLinK::snLinK(std::shared_ptr<BasisSet> primary, Options& options) : SplitJK(primary, options) {

    // always executing on host (i.e., CPU) for now
    ex_ = std::make_unique<GauXC::ExecutionSpace>(GauXC::ExecutionSpace::Host); 
    rt_ = std::make_unique<GauXC::RuntimeEnvironment>(GAUXC_MPI_CODE(MPI_COMM_WORLD));

    // convert Psi4 Molecule to GauXC molecule
    gauxc_mol_ = psi4_to_gauxc_molecule(primary_->molecule());
    // convert Psi4 basis set to GauXC basis set
    
    // construct load balancer
    const std::string load_balancer_kernel = "Default";
    gauxc_load_balancer_factory_ = std::make_unique<GauXC::LoadBalancerFactory>(*ex_, load_balancer_kernel);
    // gauxc_load_balancer_ =  
    //GauXC::LoadBalancer gauxc_load_balancer_;
    
    // construct weights module
    const std::string mol_weights_kernel = "Default";
    gauxc_mol_weights_factory_ = std::make_unique<GauXC::MolecularWeightsFactory>(*ex_, mol_weights_kernel, GauXC::MolecularWeightsSettings{});
    //GauXC::MolecularWeights gauxc_mol_weights_;
    
    // actually construct integrator
    const std::string integrator_input_type = "Replicated";
    const std::string integrator_kernel = "Default";
    const std::string reduction_kernel = "Default";
    const std::string lwd_kernel = "Default";
    const std::string dummy_func = "B3LYP";

    integrator_factory_ = std::make_unique<GauXC::XCIntegratorFactory<matrix_type> >(*ex_, integrator_input_type, integrator_kernel, lwd_kernel, reduction_kernel);
    //integrator_ = integrator_factory.get_instance(ExchCXX::functional_map.value(dummy_func_), gauxc_load_balancer_); 
    return;
}

snLinK::~snLinK() {}

size_t snLinK::num_computed_shells() {
    return num_computed_shells_;
}

void snLinK::print_header() const {
    if (print_) {
        outfile->Printf("\n");
        outfile->Printf("  ==> snLinK: Semi-Numerical Linear Exchange K <==\n\n");

        //outfile->Printf("    K Screening Cutoff: %11.0E\n", kscreen_);
        //outfile->Printf("    K Density Cutoff:   %11.0E\n", dscreen_);
        //outfile->Printf("    K Basis Cutoff:     %11.0E\n", basis_tol_);
        //outfile->Printf("    K Overlap Fitting:  %11s\n", (overlap_fitted_ ? "Yes" : "No"));
    }
}

// build the K matrix using Neeses's Chain-of-Spheres Exchange algorithm
// algorithm is originally proposed in https://doi.org/10.1016/j.chemphys.2008.10.036
// overlap fitting is discussed in https://doi.org/10.1063/1.3646921
void snLinK::build_G_component(std::vector<std::shared_ptr<Matrix>>& D, std::vector<std::shared_ptr<Matrix>>& K,
    std::vector<std::shared_ptr<TwoBodyAOInt> >& eri_computers) {
    return;
}

}  // namespace psi

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
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/basisset.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libfock/points.h"
#include "psi4/libfock/cubature.h"

#include <map>
#include <string>
#include <cmath>

namespace psi {

ZORA::ZORA(std::shared_ptr<Molecule> molecule, std::shared_ptr<BasisSet> basis, Options& options)
	: options_(options), molecule_(molecule), primary_(basis) {}

ZORA::~ZORA() {}

void ZORA::setup() {
    outfile->Printf("\n");
    outfile->Printf("          ███████╗ ██████╗ ██████╗  █████╗\n");
    outfile->Printf("          ╚══███╔╝██╔═══██╗██╔══██╗██╔══██╗\n");
    outfile->Printf("            ███╔╝ ██║   ██║██████╔╝███████║\n");
    outfile->Printf("           ███╔╝  ██║   ██║██╔══██╗██╔══██║\n");
    outfile->Printf("          ███████╗╚██████╔╝██║  ██║██║  ██║\n");
    outfile->Printf("          ╚══════╝ ╚═════╝ ╚═╝  ╚═╝╚═╝  ╚═╝\n");
    outfile->Printf("                 by Nathan Gillispie\n");
    outfile->Printf("         ===================================\n");

	
    // Print options
    outfile->Printf("\n  ==> ZORA Options <==\n");
    outfile->Printf("\n    Basis: %s", primary_->name().c_str());

	timer_on("Make Grid");

	// Initialize grid with options
	// TODO: add options_.get_str("ZORA_PRUNING_SCHEME")
	std::map<std::string, std::string> grid_str_options = {
		{"DFT_RADIAL_SCHEME",  "TREUTLER"},
		{"DFT_PRUNING_SCHEME", "TREUTLER"},
		{"DFT_NUCLEAR_SCHEME", "TREUTLER"},
		{"DFT_GRID_NAME",      ""},
		{"DFT_BLOCK_SCHEME",   "OCTREE"},
	};
	
	// TODO: add options_.get_int("ZORA_SPHERICAL_POINTS" + gridname_uppercase)
	// add options_.get_int("ZORA_RADIAL_POINTS" + gridname_uppercase
	std::map<std::string, int> grid_int_options = {
		{"DFT_BLOCK_MAX_POINTS", 256},
		{"DFT_BLOCK_MIN_POINTS", 100},
		{"DFT_SPHERICAL_POINTS", 1202},
		{"DFT_RADIAL_POINTS",    100},
	};

	std::map<std::string, double> grid_float_options = {
		{"DFT_BS_RADIUS_ALPHA",   1.0},
		{"DFT_PRUNING_ALPHA",     1.0},
		{"DFT_WEIGHTS_TOLERANCE", 1e-15},
		{"DFT_BLOCK_MAX_RADIUS",  3.0},
		{"DFT_BASIS_TOLERANCE",   1e-12},
	};

	grid_ = std::make_shared<DFTGrid>(primary_->molecule(), primary_, grid_int_options, grid_str_options, grid_float_options, options_);

	timer_off("Make Grid");

#if ZORADEBUG
	outfile->Printf("\n  ==> ZORA: Grid Details\n");
	auto npoints = grid_->npoints();
	auto nblocks = grid_->blocks().size();
	outfile->Printf("    Total number of grid points: %d \n", npoints);
	outfile->Printf("    Total number of batches: %d \n", nblocks);
	grid_->print_details();
#endif
}

void ZORA::compute(SharedMatrix T_SR) {
	timer_on("ZORA");

	setup();

	int nblocks = grid_->blocks().size();
	int max_points = grid_->max_points();
	int max_funcs = grid_->max_functions();

	timer_on("Effective Potential");
	auto veff = std::make_shared<Matrix>("Effective potential", nblocks, max_points);
	compute_veff(veff);
	timer_off("Effective Potential");

	timer_on("Scalar Relativistic Kinetic");
	BasisFunctions bf_computer(primary_, max_points, max_funcs);
	bf_computer.set_deriv(1);
	compute_TSR(bf_computer, veff, T_SR);
	timer_off("Scalar Relativistic Kinetic");

	timer_off("ZORA");
#if ZORADEBUG
	T_SR->print();
#endif
}

void ZORA::compute_veff(SharedMatrix veff)
{
	// Speed of light in atomic units
	double C = pc_c_au;

	int natoms = molecule_->natom();
	for (int a = 0; a < natoms; a++) {
		int Z = molecule_->Z(a);
		if (Z > 104) throw PSIEXCEPTION("Z too big. Max value 104");
		auto pos_a = molecule_->xyz(a);

		const double* coef_a  = &coeffs[c_aIndex[Z-1]];
		const double* alpha_a = &alphas[c_aIndex[Z-1]];
		int nc_a = c_aIndex[Z] - c_aIndex[Z-1];

		int index = 0;
		for (const auto &block : grid_->blocks()) {
			int npoints = block->npoints();

			double* x = block->x();
			double* y = block->y();
			double* z = block->z();

			//einsums("i,ip->p", 𝕔[i], erf(α[i]⊗ r[p]))/r[p]
			for (int p = 0; p < npoints; p++) {
				double dist = std::hypot(pos_a[0]-x[p], pos_a[1]-y[p], pos_a[2]-z[p]);
				double outer = 0;
				for (int i = 0; i < nc_a; i++) {
					outer += std::erf(dist * alpha_a[i]) * coef_a[i];
				}
				outer /= dist;
				outer -= Z/dist;
				veff->add(index, p, outer);
			}
			index++;
		}
	}
}

//Scalar Relativistic Kinetic Energy Matrix
void ZORA::compute_TSR(BasisFunctions &props, SharedMatrix veff, SharedMatrix &T_SR)
{
	// Speed of light in atomic units
	double C = pc_c_au;
	double** T_SRp = T_SR->pointer();

	double* kernel = new double[props.max_points()];

	int index = 0;
	for (const auto &block : grid_->blocks()) {
		const auto &bf_map = block->functions_local_to_global();
		auto local_nbf = bf_map.size();
		int npoints = block->npoints();

		props.compute_functions(block);
		auto phi_x = props.basis_value("PHI_X");
		auto phi_y = props.basis_value("PHI_Y");
		auto phi_z = props.basis_value("PHI_Z");

		double* w = block->w();

		//preprocess kernel c²/(2c²-veff) * weight
		for (int p = 0; p < npoints; p++) {
			kernel[p] = C *C /(2.*C *C - veff->get(index,p)) * w[p];
		}

		for (int l_mu = 0; l_mu < local_nbf; l_mu++) {
			int mu = bf_map[l_mu];
			for (int l_nu = l_mu; l_nu < local_nbf; l_nu++) {
				int nu = bf_map[l_nu];
				for (int p = 0; p < npoints; p++) {
					T_SRp[mu][nu] += kernel[p] * (
						phi_x->get(p,l_mu)*phi_x->get(p,l_nu) +
						phi_y->get(p,l_mu)*phi_y->get(p,l_nu) +
						phi_z->get(p,l_mu)*phi_z->get(p,l_nu) );

				}
			}
		}
		index++;
	}

	T_SR->copy_upper_to_lower();
	delete kernel;
}

}  // namespace psi

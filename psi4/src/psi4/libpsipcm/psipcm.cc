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
#ifdef USING_PCMSolver

#include <algorithm>
#include <array>
#include <memory>
#include <tuple>
#include <vector>

#include "psipcm.h"

#include "psi4/psi4-dec.h"

#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/petitelist.h"
#include "psi4/libmints/potentialint.h"
#include "psi4/libmints/vector.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#include <PCMSolver/PCMInput.h>

namespace psi {

namespace detail {
std::tuple<std::vector<double>, std::vector<double>> collect_atoms(std::shared_ptr<Molecule> molecule) {
    int nat = molecule->natom();
    std::vector<double> charges(nat);
    for (int i = 0; i < nat; ++i) {
        charges[i] = molecule->fZ(i);
    }

    Matrix geom = molecule->geometry();
    std::vector<double> centers(3 * nat);
    for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < nat; ++i) {
            // Eigen stores matrices in Column-Major order
            centers[j + i * 3] = geom.get(i, j);
        }
    }

    return std::make_tuple(charges, centers);
}

PCMInput pcmsolver_input() {
    PCMInput host_input;

    // These parameters would be set by the host input reading
    // Length and area parameters are all assumed to be in Angstrom,
    // the module will convert to Bohr internally
    strcpy(host_input.cavity_type, "gepol");
    host_input.patch_level = 2;
    host_input.coarsity = 0.5;
    host_input.area = 0.2;
    host_input.min_distance = 0.1;
    host_input.der_order = 4;
    host_input.scaling = true;
    strcpy(host_input.radii_set, "bondi");
    strcpy(host_input.restart_name, "cavity.npz");
    host_input.min_radius = 100.0;

    strcpy(host_input.solver_type, "iefpcm");
    strcpy(host_input.solvent, "water");
    strcpy(host_input.equation_type, "secondkind");
    host_input.correction = 0.0;
    host_input.probe_radius = 1.0;

    strcpy(host_input.inside_type, "vacuum");
    host_input.outside_epsilon = 1.0;
    strcpy(host_input.outside_type, "uniformdielectric");

    return host_input;
}

void host_writer(const char *message) {
    outfile->Printf(message);
    outfile->Printf("\n");
}
}

PCM::PCM(int print_level, std::shared_ptr<BasisSet> basisset) : pcm_print_(print_level), basisset_(basisset) {
    if (!pcmsolver_is_compatible_library()) throw PSIEXCEPTION("Incompatible PCMSolver library version.");

    std::shared_ptr<Molecule> molecule = basisset_->molecule();
    int natom = molecule->natom();

    auto integrals = std::make_shared<IntegralFactory>(basisset, basisset, basisset, basisset);

    PetiteList petite(basisset, integrals, true);
    my_aotoso_ = petite.aotoso();

    potential_int_ = static_cast<PCMPotentialInt *>(integrals->pcm_potentialint());

    /* PCMSolver needs to know who has to parse the input.
     * We should have something like this here:
     * int PSI4_provides_input = false;
     * if (PSI4_has_pcmsolver_input) {
     *    PSI4_provides_input = true;
     * }
     */
    std::vector<double> charges(natom), coordinates(3 * natom);
    std::tie(charges, coordinates) = detail::collect_atoms(molecule);
    std::array<int, 4> symmetry_info{{0, 0, 0, 0}};
    int PSI4_provides_input = false;
    PCMInput host_input;
    if (PSI4_provides_input) {
        host_input = detail::pcmsolver_input();
        context_ = std::shared_ptr<pcmsolver_context_t>(
            pcmsolver_new(PCMSOLVER_READER_HOST, natom, charges.data(), coordinates.data(), symmetry_info.data(),
                          &host_input, detail::host_writer),
            pcmsolver_delete);
    } else {
        context_ = std::shared_ptr<pcmsolver_context_t>(
            pcmsolver_new(PCMSOLVER_READER_OWN, natom, charges.data(), coordinates.data(), symmetry_info.data(),
                          &host_input, detail::host_writer),
            pcmsolver_delete);
    }
    outfile->Printf("  **PSI4:PCMSOLVER Interface Active**\n");
    pcmsolver_print(context_.get());
    ntess_ = pcmsolver_get_cavity_size(context_.get());
    ntessirr_ = pcmsolver_get_irreducible_cavity_size(context_.get());
    outfile->Printf("  There are %d tesserae, %d of which irreducible.\n\n", ntess_, ntessirr_);
    tess_pot_ = new double[ntess_];
    tess_pot_e_ = new double[ntess_];
    tess_pot_n_ = new double[ntess_];
    tess_charges_e_ = new double[ntess_];
    tess_charges_n_ = new double[ntess_];
    tess_charges_ = new double[ntess_];

    // The charge and {x,y,z} coordinates (in bohr) for each tessera
    tess_Zxyz_ = std::make_shared<Matrix>("Tess Zxyz", ntess_, 4);
    double **ptess_Zxyz = tess_Zxyz_->pointer();
    // Set the tesserae's coordinates (note the loop bounds; this function is 1-based)
    for (int tess = 1; tess <= ntess_; ++tess) pcmsolver_get_center(context_.get(), tess, &(ptess_Zxyz[tess - 1][1]));

    // Compute the nuclear potentials at the tesserae
    std::fill(tess_pot_n_, tess_pot_n_ + ntess_, 0.0);
    for (int atom = 0; atom < natom; ++atom) {
        double Z = charges[atom];
        for (int tess = 0; tess < ntess_; ++tess) {
            double dx = ptess_Zxyz[tess][1] - coordinates[atom * 3];
            double dy = ptess_Zxyz[tess][2] - coordinates[atom * 3 + 1];
            double dz = ptess_Zxyz[tess][3] - coordinates[atom * 3 + 2];
            double r = std::sqrt(dx * dx + dy * dy + dz * dz);
            tess_pot_n_[tess] += Z / r;
            if (r < 1.0E-3) outfile->Printf("Warning! Tessera %d is only %.3f bohr from atom %d!\n", tess, r, atom + 1);
        }
    }

    // A little debug info
    if (pcm_print_ > 2) {
        outfile->Printf("Nuclear MEP at each tessera:\n");
        outfile->Printf("----------------------------\n");
        for (int tess = 0; tess < ntess_; ++tess) outfile->Printf("tess[%4d] -> %16.10f\n", tess, tess_pot_n_[tess]);
    }

    // Compute the nuclear charges, since they don't change
    std::fill(tess_charges_n_, tess_charges_n_ + ntess_, 0.0);
    const char *potential_name = "NucMEP";
    const char *charge_name = "NucASC";
    pcmsolver_set_surface_function(context_.get(), ntess_, tess_pot_n_, potential_name);
    int irrep = 0;
    pcmsolver_compute_asc(context_.get(), potential_name, charge_name, irrep);
    pcmsolver_get_surface_function(context_.get(), ntess_, tess_charges_n_, charge_name);

    // A little debug info
    if (pcm_print_ > 2) {
        outfile->Printf("Nuclear ASC at each tessera:\n");
        outfile->Printf("----------------------------\n");
        for (int tess = 0; tess < ntess_; ++tess)
            outfile->Printf("tess[%4d] -> %16.10f\n", tess, tess_charges_n_[tess]);
    }

}  // PCM()

PCM::~PCM() {
    delete[] tess_pot_;
    delete[] tess_pot_e_;
    delete[] tess_pot_n_;
    delete[] tess_charges_;
    delete[] tess_charges_e_;
    delete[] tess_charges_n_;
}

double PCM::compute_E(SharedMatrix &D, CalcType type) {
    switch (type) {
        case CalcType::Total:
            return compute_E_total(D);
        case CalcType::NucAndEle:
            return compute_E_separate(D);
        case CalcType::EleOnly:
            return compute_E_electronic(D);
        default:
            throw PSIEXCEPTION("Unknown PCM calculation type.");
    }
}

double PCM::compute_E_total(SharedMatrix &D) {
    double **ptess_Zxyz = tess_Zxyz_->pointer();
    std::fill(tess_pot_e_, tess_pot_e_ + ntess_, 0.0);
    std::fill(tess_charges_e_, tess_charges_e_ + ntess_, 0.0);
    for (int tess = 0; tess < ntess_; ++tess) ptess_Zxyz[tess][0] = 1.0;
    potential_int_->set_charge_field(tess_Zxyz_);

    SharedMatrix D_carts;
    if (basisset_->has_puream()) {
        D_carts = std::make_shared<Matrix>("D carts", basisset_->nao(), basisset_->nao());
        D_carts->back_transform(D, my_aotoso_);
    } else
        D_carts = D;

    ContractOverDensityFunctor contract_density_functor(ntess_, tess_pot_e_, D_carts);
    // Add in the electronic contribution to the potential at each tessera
    potential_int_->compute(contract_density_functor);

    // A little debug info
    if (pcm_print_ > 2) {
        outfile->Printf("Electronic MEP at each tessera:\n");
        outfile->Printf("-------------------------------\n");
        for (int tess = 0; tess < ntess_; ++tess) outfile->Printf("tess[%4d] -> %16.10f\n", tess, tess_pot_e_[tess]);
    }

    // Combine the nuclear and electronic potentials at each tessera
    for (int tess = 0; tess < ntess_; ++tess) tess_pot_[tess] = tess_pot_n_[tess] + tess_pot_e_[tess];

    const char *tot_potential_name = "TotMEP";
    const char *tot_charge_name = "TotASC";
    pcmsolver_set_surface_function(context_.get(), ntess_, tess_pot_, tot_potential_name);
    int irrep = 0;
    pcmsolver_compute_asc(context_.get(), tot_potential_name, tot_charge_name, irrep);
    pcmsolver_get_surface_function(context_.get(), ntess_, tess_charges_, tot_charge_name);

    if (pcm_print_ > 2) {
        outfile->Printf("Total MEP & ASC at each tessera:\n");
        outfile->Printf("-------------------------------------------------\n");
        outfile->Printf("Tessera# Total MEP       Total ASC\n");
        outfile->Printf("----------------------------------------------------------------------------\n");
        for (int tess = 0; tess < ntess_; ++tess)
            outfile->Printf("%4d  %16.10f  %16.10f\n", tess, tess_pot_[tess], tess_charges_[tess]);
    }

    // Grab the polarization energy from PCMSolver
    double Epol = pcmsolver_compute_polarization_energy(context_.get(), tot_potential_name, tot_charge_name);
    outfile->Printf("   PCM polarization energy = %16.14f\n", Epol);

    return Epol;
}

double PCM::compute_E_separate(SharedMatrix &D) {
    double **ptess_Zxyz = tess_Zxyz_->pointer();
    std::fill(tess_pot_e_, tess_pot_e_ + ntess_, 0.0);
    std::fill(tess_charges_, tess_charges_ + ntess_, 0.0);
    std::fill(tess_charges_e_, tess_charges_e_ + ntess_, 0.0);
    for (int tess = 0; tess < ntess_; ++tess) ptess_Zxyz[tess][0] = 1.0;
    potential_int_->set_charge_field(tess_Zxyz_);

    SharedMatrix D_carts;
    if (basisset_->has_puream()) {
        D_carts = std::make_shared<Matrix>("D carts", basisset_->nao(), basisset_->nao());
        D_carts->back_transform(D, my_aotoso_);
    } else
        D_carts = D;

    ContractOverDensityFunctor contract_density_functor(ntess_, tess_pot_e_, D_carts);
    // Add in the electronic contribution to the potential at each tessera
    potential_int_->compute(contract_density_functor);

    // A little debug info
    if (pcm_print_ > 2) {
        outfile->Printf("Electronic MEP at each tessera:\n");
        outfile->Printf("-------------------------------\n");
        for (int tess = 0; tess < ntess_; ++tess) outfile->Printf("tess[%4d] -> %16.10f\n", tess, tess_pot_e_[tess]);
    }

    const char *e_potential_name = "EleMEP";
    const char *e_charge_name = "EleASC";
    pcmsolver_set_surface_function(context_.get(), ntess_, tess_pot_e_, e_potential_name);
    int irrep = 0;
    pcmsolver_compute_asc(context_.get(), e_potential_name, e_charge_name, irrep);
    pcmsolver_get_surface_function(context_.get(), ntess_, tess_charges_e_, e_charge_name);

    // Combine the nuclear and electronic potentials at each tessera
    for (int tess = 0; tess < ntess_; ++tess) tess_pot_[tess] = tess_pot_n_[tess] + tess_pot_e_[tess];

    if (pcm_print_ > 2) {
        outfile->Printf("Nuclear and Electronic MEP & ASC at each tessera:\n");
        outfile->Printf("-------------------------------------------------\n");
        outfile->Printf("Tessera# Nuclear MEP       Nuclear ASC       Elec. MEP         Elec. ASC:\n");
        outfile->Printf("----------------------------------------------------------------------------\n");
        for (int tess = 0; tess < ntess_; ++tess)
            outfile->Printf("%4d  %16.10f  %16.10f  %16.10f  %16.10f\n", tess, tess_pot_n_[tess], tess_charges_n_[tess],
                            tess_pot_e_[tess], tess_charges_e_[tess]);
    }

    const char *tot_potential_name = "TotMEP";
    const char *tot_charge_name = "TotASC";
    pcmsolver_set_surface_function(context_.get(), ntess_, tess_pot_, tot_potential_name);
    pcmsolver_compute_asc(context_.get(), tot_potential_name, tot_charge_name, irrep);
    pcmsolver_get_surface_function(context_.get(), ntess_, tess_charges_, tot_charge_name);

    // Grab the polarization energy from PCMSolver
    const char *n_potential_name = "NucMEP";
    const char *n_charge_name = "NucASC";
    double U_NN = pcmsolver_compute_polarization_energy(context_.get(), n_potential_name, n_charge_name);
    double U_eN = pcmsolver_compute_polarization_energy(context_.get(), n_potential_name, e_charge_name);
    double U_Ne = pcmsolver_compute_polarization_energy(context_.get(), e_potential_name, n_charge_name);
    double U_ee = pcmsolver_compute_polarization_energy(context_.get(), e_potential_name, e_charge_name);
    outfile->Printf("  Polarization energy components\n");
    outfile->Printf("  U_ee = %16.14f\n", 2.0 * U_ee);
    outfile->Printf("  U_eN = %16.14f\n", 2.0 * U_eN);
    outfile->Printf("  U_Ne = %16.14f\n", 2.0 * U_Ne);
    outfile->Printf("  U_NN = %16.14f\n", 2.0 * U_NN);
    outfile->Printf("  U_eN - U_Ne = %16.14f\n", U_eN - U_Ne);
    double Epol = U_NN + U_eN + U_Ne + U_ee;
    outfile->Printf("   PCM polarization energy = %16.14f\n", Epol);
    return Epol;
}

double PCM::compute_E_electronic(SharedMatrix &D) {
    double **ptess_Zxyz = tess_Zxyz_->pointer();
    std::fill(tess_pot_e_, tess_pot_e_ + ntess_, 0.0);
    std::fill(tess_charges_e_, tess_charges_e_ + ntess_, 0.0);
    for (int tess = 0; tess < ntess_; ++tess) ptess_Zxyz[tess][0] = 1.0;
    potential_int_->set_charge_field(tess_Zxyz_);

    SharedMatrix D_carts;
    if (basisset_->has_puream()) {
        D_carts = std::make_shared<Matrix>("D carts", basisset_->nao(), basisset_->nao());
        D_carts->back_transform(D, my_aotoso_);
    } else
        D_carts = D;

    ContractOverDensityFunctor contract_density_functor(ntess_, tess_pot_e_, D_carts);
    // Add in the electronic contribution to the potential at each tessera
    potential_int_->compute(contract_density_functor);

    // A little debug info
    if (pcm_print_ > 2) {
        outfile->Printf("Electronic MEP at each tessera:\n");
        outfile->Printf("-------------------------------\n");
        for (int tess = 0; tess < ntess_; ++tess) outfile->Printf("tess[%4d] -> %16.10f\n", tess, tess_pot_e_[tess]);
    }

    const char *e_potential_name = "EleMEP";
    const char *e_charge_name = "EleASC";
    pcmsolver_set_surface_function(context_.get(), ntess_, tess_pot_e_, e_potential_name);
    int irrep = 0;
    pcmsolver_compute_asc(context_.get(), e_potential_name, e_charge_name, irrep);
    pcmsolver_get_surface_function(context_.get(), ntess_, tess_charges_e_, e_charge_name);

    if (pcm_print_ > 2) {
        outfile->Printf("Electronic MEP & ASC at each tessera:\n");
        outfile->Printf("-------------------------------------------------\n");
        outfile->Printf("Tessera# Elec. MEP         Elec. ASC:\n");
        outfile->Printf("----------------------------------------------------------------------------\n");
        for (int tess = 0; tess < ntess_; ++tess)
            outfile->Printf("%4d  %16.10f  %16.10f\n", tess, tess_pot_e_[tess], tess_charges_e_[tess]);
    }

    // Grab the polarization energy from PCMSolver
    double Epol = pcmsolver_compute_polarization_energy(context_.get(), e_potential_name, e_charge_name);
    outfile->Printf("   PCM polarization energy (electronic only) = %16.14f\n", Epol);
    return Epol;
}

SharedMatrix PCM::compute_V() {
    auto V_pcm_cart = std::make_shared<Matrix>("PCM potential cart", basisset_->nao(), basisset_->nao());
    ContractOverChargesFunctor contract_charges_functor(tess_charges_, V_pcm_cart);
    potential_int_->compute(contract_charges_functor);
    // The potential might need to be transformed to the spherical harmonic basis
    SharedMatrix V_pcm_pure;
    if (basisset_->has_puream()) {
        V_pcm_pure = std::make_shared<Matrix>("PCM potential pure", basisset_->nbf(), basisset_->nbf());
        V_pcm_pure->transform(V_pcm_cart, my_aotoso_);
    }
    if (basisset_->has_puream())
        return V_pcm_pure;
    else
        return V_pcm_cart;
}

SharedMatrix PCM::compute_V_electronic() {
    auto V_pcm_cart = std::make_shared<Matrix>("PCM potential cart", basisset_->nao(), basisset_->nao());
    ContractOverChargesFunctor contract_charges_functor(tess_charges_e_, V_pcm_cart);
    potential_int_->compute(contract_charges_functor);
    // The potential might need to be transformed to the spherical harmonic basis
    SharedMatrix V_pcm_pure;
    if (basisset_->has_puream()) {
        V_pcm_pure = std::make_shared<Matrix>("PCM potential pure", basisset_->nbf(), basisset_->nbf());
        V_pcm_pure->transform(V_pcm_cart, my_aotoso_);
    }
    if (basisset_->has_puream())
        return V_pcm_pure;
    else
        return V_pcm_cart;
}

}  // psi namespace
#endif

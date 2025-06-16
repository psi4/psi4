/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
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

#include "psipcm.h"

#include "psi4/psi4-dec.h"

#include "psi4/libmints/basisset.h"
#include "psi4/libmints/dimension.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/petitelist.h"
#include "psi4/libmints/potentialint.h"
#include "psi4/libmints/vector.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#include <PCMSolver/PCMInput.h>

#include <algorithm>
#include <array>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace psi {

namespace detail {
std::pair<std::vector<double>, std::vector<double>> collect_atoms(std::shared_ptr<Molecule> molecule) {
    int nat = molecule->natom();
    std::vector<double> charges(nat);
    for (int i = 0; i < nat; ++i) {
        charges[i] = static_cast<double>(molecule->true_atomic_number(i));
    }

    Matrix geom = molecule->geometry();
    std::vector<double> centers(3 * nat);
    for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < nat; ++i) {
            // Eigen stores matrices in Column-Major order
            centers[j + i * 3] = geom.get(i, j);
        }
    }

    return std::make_pair(charges, centers);
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

std::shared_ptr<pcmsolver_context_t> init_PCMSolver(const std::string &pcmsolver_parsed_fname,
                                                    const std::shared_ptr<Molecule> &molecule) {
    /* PCMSolver needs to know who has to parse the input.
     * We should have something like this here:
     * int PSI4_provides_input = false;
     * if (PSI4_has_pcmsolver_input) {
     *    PSI4_provides_input = true;
     * }
     */
    int natom = molecule->natom();
    std::vector<double> charges(natom), coordinates(3 * natom);
    std::tie(charges, coordinates) = detail::collect_atoms(molecule);
    // TODO grab point group generators into this array
    std::array<int, 4> symmetry_info{{0, 0, 0, 0}};
    int PSI4_provides_input = false;
    PCMInput host_input;
    if (PSI4_provides_input) {
        host_input = detail::pcmsolver_input();
        return std::shared_ptr<pcmsolver_context_t>(
            pcmsolver_new(PCMSOLVER_READER_HOST, natom, charges.data(), coordinates.data(), symmetry_info.data(),
                          &host_input, detail::host_writer),
            pcmsolver_delete);
    } else {
        return std::shared_ptr<pcmsolver_context_t>(
            pcmsolver_new_v1112(PCMSOLVER_READER_OWN, natom, charges.data(), coordinates.data(), symmetry_info.data(),
                                pcmsolver_parsed_fname.c_str(), &host_input, detail::host_writer),
            pcmsolver_delete);
    }
}
}  // namespace detail

PCM::PCM(const std::string &pcmsolver_parsed_fname, int print_level, std::shared_ptr<BasisSet> basisset)
    : pcmsolver_parsed_fname_(pcmsolver_parsed_fname), pcm_print_(print_level), basisset_(basisset) {
    if (!pcmsolver_is_compatible_library()) throw PSIEXCEPTION("Incompatible PCMSolver library version.");

    std::shared_ptr<Molecule> molecule = basisset_->molecule();

    auto integrals = std::make_shared<IntegralFactory>(basisset, basisset, basisset, basisset);

    PetiteList petite(basisset, integrals, true);
    my_aotoso_ = petite.aotoso();

    potential_int_ = static_cast<PCMPotentialInt *>(integrals->pcm_potentialint().release());

    context_ = detail::init_PCMSolver(pcmsolver_parsed_fname_, molecule);
    outfile->Printf("  **PSI4:PCMSOLVER Interface Active**\n");
    pcmsolver_citation(detail::host_writer);
    pcmsolver_print(context_.get());
    ntess_ = pcmsolver_get_cavity_size(context_.get());
    ntessirr_ = pcmsolver_get_irreducible_cavity_size(context_.get());
    std::vector<int> tmp(int(ntess_ / ntessirr_));
    std::fill(tmp.begin(), tmp.end(), ntessirr_);
    tesspi_ = Dimension(tmp);

    // The charge and {x,y,z} coordinates (in bohr) for each tessera
    tess_Zxyz_ = std::make_shared<Matrix>("Tess Zxyz", ntess_, 4);
    double **ptess_Zxyz = tess_Zxyz_->pointer();
    // Set the tesserae's coordinates (note the loop bounds; this function is 1-based)
    for (int tess = 1; tess <= ntess_; ++tess) pcmsolver_get_center(context_.get(), tess, &(ptess_Zxyz[tess - 1][1]));

    int natom = molecule->natom();
    std::vector<double> charges(natom), coordinates(3 * natom);
    std::tie(charges, coordinates) = detail::collect_atoms(molecule);
    MEP_n_ = std::make_shared<Vector>(tesspi_);
    // Compute the nuclear potentials at the tesserae
    for (int atom = 0; atom < natom; ++atom) {
        double Z = static_cast<double>(molecule->Z(atom));
        for (int tess = 0; tess < ntess_; ++tess) {
            double dx = ptess_Zxyz[tess][1] - coordinates[atom * 3];
            double dy = ptess_Zxyz[tess][2] - coordinates[atom * 3 + 1];
            double dz = ptess_Zxyz[tess][3] - coordinates[atom * 3 + 2];
            double r = std::sqrt(dx * dx + dy * dy + dz * dz);
            MEP_n_->add(0, tess, Z / r);
        }
    }

    // Compute the nuclear charges, since they don't change
    std::string MEP_n_label("NucMEP");
    std::string ASC_n_label("NucASC");
    pcmsolver_set_surface_function(context_.get(), ntess_, MEP_n_->pointer(0), MEP_n_label.c_str());
    int irrep = 0;
    pcmsolver_compute_asc(context_.get(), MEP_n_label.c_str(), ASC_n_label.c_str(), irrep);

    if (pcm_print_ > 2) {
        auto tess_charges_n = std::make_shared<Vector>(tesspi_);
        pcmsolver_get_surface_function(context_.get(), ntess_, tess_charges_n->pointer(0), ASC_n_label.c_str());
        outfile->Printf("Nuclear MEP & ASC at each tessera:\n");
        outfile->Printf("-------------------------------------------------\n");
        outfile->Printf("Tessera# Nuclear MEP       Nuclear ASC\n");
        outfile->Printf("----------------------------------------------------------------------------\n");
        for (int tess = 0; tess < ntess_; ++tess)
            outfile->Printf("%4d  %16.10f  %16.10f\n", tess, MEP_n_->get(0, tess), tess_charges_n->get(0, tess));
    }
}

PCM::PCM(const PCM *other) {
    ntess_ = other->ntess_;
    ntessirr_ = other->ntessirr_;
    tesspi_ = other->tesspi_;
    tess_Zxyz_ = other->tess_Zxyz_->clone();
    MEP_n_ = std::make_shared<Vector>(std::move(other->MEP_n_->clone()));
    basisset_ = other->basisset_;
    my_aotoso_ = other->my_aotoso_->clone();
    potential_int_ = other->potential_int_;
    context_ = detail::init_PCMSolver(other->pcmsolver_parsed_fname_, basisset_->molecule());
    pcm_print_ = other->pcm_print_;
}

SharedVector PCM::compute_electronic_MEP(const SharedMatrix &D) const {
    double **pZxyz = tess_Zxyz_->pointer();
    std::vector<std::pair<double, std::array<double, 3>>> field;
    for (int tess = 0; tess < ntess_; ++tess) {
        field.push_back({1.0, {pZxyz[tess][1], pZxyz[tess][2], pZxyz[tess][3]}});
    }
    potential_int_->set_charge_field(field);

    auto MEP = std::make_shared<Vector>(tesspi_);
    ContractOverDensityFunctor contract_density_functor(ntess_, MEP->pointer(0), D);
    // Add in the electronic contribution to the potential at each tessera
    potential_int_->compute(contract_density_functor);

    // A little debug info
    if (pcm_print_ > 2) {
        outfile->Printf("Electronic MEP at each tessera:\n");
        outfile->Printf("-------------------------------\n");
        for (int tess = 0; tess < ntess_; ++tess) outfile->Printf("tess[%4d] -> %16.10f\n", tess, MEP->get(0, tess));
    }

    return MEP;
}

std::pair<double, SharedMatrix> PCM::compute_PCM_terms(const SharedMatrix &D, CalcType type) const {
    double upcm = 0.0;
    auto MEP_e = compute_electronic_MEP(D);
    auto ASC = std::make_shared<Vector>(tesspi_);
    switch (type) {
        case CalcType::Total:
            upcm = compute_E_total(MEP_e);
            pcmsolver_get_surface_function(context_.get(), ntess_, ASC->pointer(0), "TotASC");
            break;
        case CalcType::NucAndEle:
            upcm = compute_E_separate(MEP_e);
            pcmsolver_get_surface_function(context_.get(), ntess_, ASC->pointer(0), "TotASC");
            break;
        case CalcType::EleOnly:
            upcm = compute_E_electronic(MEP_e);
            pcmsolver_get_surface_function(context_.get(), ntess_, ASC->pointer(0), "EleASC");
            break;
        default:
            throw PSIEXCEPTION("Unknown PCM calculation type.");
    }
    return std::make_pair(upcm, compute_Vpcm(ASC));
}

SharedMatrix PCM::compute_V(const SharedMatrix &D) {
    // returns V_elec for a given D with CalcType=total and response_asc if needed
    auto MEP_e = compute_electronic_MEP(D);
    auto ASC = std::make_shared<Vector>(tesspi_);
    std::string MEP_label("OITMEP");  // one-index-transformed MEP
    std::string ASC_label("OITASC");  // one-index-transformed ASC
    pcmsolver_set_surface_function(context_.get(), ntess_, MEP_e->pointer(0), MEP_label.c_str());
    int irrep = 0;
    pcmsolver_compute_response_asc(context_.get(), MEP_label.c_str(), ASC_label.c_str(), irrep);
    pcmsolver_get_surface_function(context_.get(), ntess_, ASC->pointer(0), ASC_label.c_str());
    return compute_Vpcm(ASC);
}

double PCM::compute_E_total(const SharedVector &MEP_e) const {
    // Combine the nuclear and electronic potentials at each tessera
    MEP_e->add(*MEP_n_);
    std::string MEP_label("TotMEP");
    std::string ASC_label("TotASC");
    pcmsolver_set_surface_function(context_.get(), ntess_, MEP_e->pointer(0), MEP_label.c_str());
    int irrep = 0;
    pcmsolver_compute_asc(context_.get(), MEP_label.c_str(), ASC_label.c_str(), irrep);

    if (pcm_print_ > 2) {
        auto ASC = std::make_shared<Vector>(tesspi_);
        pcmsolver_get_surface_function(context_.get(), ntess_, ASC->pointer(0), ASC_label.c_str());
        outfile->Printf("Total MEP & ASC at each tessera:\n");
        outfile->Printf("-------------------------------------------------\n");
        outfile->Printf("Tessera# Total MEP       Total ASC\n");
        outfile->Printf("----------------------------------------------------------------------------\n");
        for (int tess = 0; tess < ntess_; ++tess)
            outfile->Printf("%4d  %16.10f  %16.10f\n", tess, MEP_e->get(0, tess), ASC->get(0, tess));
    }

    // Grab the polarization energy from PCMSolver
    double Epol = pcmsolver_compute_polarization_energy(context_.get(), MEP_label.c_str(), ASC_label.c_str());
    return Epol;
}

double PCM::compute_E_separate(const SharedVector &MEP_e) const {
    // Surface functions labels
    std::string MEP_n_label("NucMEP");
    std::string ASC_n_label("NucASC");
    std::string MEP_e_label("EleMEP");
    std::string ASC_e_label("EleASC");
    std::string MEP_label("TotMEP");
    std::string ASC_label("TotASC");

    pcmsolver_set_surface_function(context_.get(), ntess_, MEP_e->pointer(0), MEP_e_label.c_str());
    int irrep = 0;
    pcmsolver_compute_asc(context_.get(), MEP_e_label.c_str(), ASC_e_label.c_str(), irrep);

    if (pcm_print_ > 2) {
        auto ASC_n = std::make_shared<Vector>(tesspi_);
        pcmsolver_get_surface_function(context_.get(), ntess_, ASC_n->pointer(0), ASC_n_label.c_str());
        auto ASC_e = std::make_shared<Vector>(tesspi_);
        pcmsolver_get_surface_function(context_.get(), ntess_, ASC_e->pointer(0), ASC_e_label.c_str());
        outfile->Printf("Nuclear and Electronic MEP & ASC at each tessera:\n");
        outfile->Printf("-------------------------------------------------\n");
        outfile->Printf("Tessera# Nuclear MEP       Nuclear ASC       Elec. MEP         Elec. ASC:\n");
        outfile->Printf("----------------------------------------------------------------------------\n");
        for (int tess = 0; tess < ntess_; ++tess)
            outfile->Printf("%4d  %16.10f  %16.10f  %16.10f  %16.10f\n", tess, MEP_n_->get(0, tess),
                            ASC_n->get(0, tess), MEP_e->get(0, tess), ASC_e->get(0, tess));
    }

    MEP_e->add(*MEP_n_);
    pcmsolver_set_surface_function(context_.get(), ntess_, MEP_e->pointer(0), MEP_label.c_str());
    pcmsolver_compute_asc(context_.get(), MEP_label.c_str(), ASC_label.c_str(), irrep);

    // Grab the polarization energy from PCMSolver
    double U_NN = pcmsolver_compute_polarization_energy(context_.get(), MEP_n_label.c_str(), ASC_n_label.c_str());
    double U_eN = pcmsolver_compute_polarization_energy(context_.get(), MEP_n_label.c_str(), ASC_e_label.c_str());
    double U_Ne = pcmsolver_compute_polarization_energy(context_.get(), MEP_e_label.c_str(), ASC_n_label.c_str());
    double U_ee = pcmsolver_compute_polarization_energy(context_.get(), MEP_e_label.c_str(), ASC_e_label.c_str());
    outfile->Printf("  Polarization energy components\n");
    outfile->Printf("  U_ee = %16.14f\n", 2.0 * U_ee);
    outfile->Printf("  U_eN = %16.14f\n", 2.0 * U_eN);
    outfile->Printf("  U_Ne = %16.14f\n", 2.0 * U_Ne);
    outfile->Printf("  U_NN = %16.14f\n", 2.0 * U_NN);
    outfile->Printf("  U_eN - U_Ne = %16.14f\n", U_eN - U_Ne);
    double Epol = U_NN + U_eN + U_Ne + U_ee;
    return Epol;
}

double PCM::compute_E_electronic(const SharedVector &MEP_e) const {
    std::string MEP_e_label("EleMEP");
    std::string ASC_e_label("EleASC");
    pcmsolver_set_surface_function(context_.get(), ntess_, MEP_e->pointer(0), MEP_e_label.c_str());
    int irrep = 0;
    pcmsolver_compute_asc(context_.get(), MEP_e_label.c_str(), ASC_e_label.c_str(), irrep);

    if (pcm_print_ > 2) {
        auto ASC_e = std::make_shared<Vector>(tesspi_);
        pcmsolver_get_surface_function(context_.get(), ntess_, ASC_e->pointer(0), ASC_e_label.c_str());
        outfile->Printf("Electronic MEP & ASC at each tessera:\n");
        outfile->Printf("-------------------------------------------------\n");
        outfile->Printf("Tessera# Elec. MEP         Elec. ASC:\n");
        outfile->Printf("----------------------------------------------------------------------------\n");
        for (int tess = 0; tess < ntess_; ++tess)
            outfile->Printf("%4d  %16.10f  %16.10f\n", tess, MEP_e->get(0, tess), ASC_e->get(0, tess));
    }

    // Grab the polarization energy from PCMSolver
    double Epol = pcmsolver_compute_polarization_energy(context_.get(), MEP_e_label.c_str(), ASC_e_label.c_str());
    outfile->Printf("   PCM polarization energy (electronic only) = %16.14f\n", Epol);
    return Epol;
}

SharedMatrix PCM::compute_Vpcm(const SharedVector &ASC) const {
    auto V_pcm = std::make_shared<Matrix>("PCM potential cart", basisset_->nbf(), basisset_->nbf());
    ContractOverChargesFunctor contract_charges_functor(ASC->pointer(0), V_pcm);
    potential_int_->compute(contract_charges_functor);
    return V_pcm;
}
}  // namespace psi

#endif

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
#include "psi4/libmints/wavefunction.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libparallel/parallel.h"
#include "psi4/psifiles.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libiwl/iwl.h"
#include "psi4/libqt/qt.h"
#include "psi4/libparallel/parallel.h"
#include "psi4/libmints/view.h"
#include "psi4/libmints/orbitalspace.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/petitelist.h"
#include "psi4/libmints/sobasis.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/factory.h"
#include "psi4/libmints/sointegral_onebody.h"
#include "psi4/libmints/corrtab.h"
#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/exception.h"

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <tuple>

using namespace psi;

// Globals (Seriously? This is where we instantiate these.)
size_t ioff[MAX_IOFF];
double df[MAX_DF];
double bc[MAX_BC][MAX_BC];
double fac[MAX_FAC];

Wavefunction::Wavefunction(std::shared_ptr<Molecule> molecule,
                           std::shared_ptr<BasisSet> basis,
                           Options &options) :
        options_(options), basisset_(basis), molecule_(molecule)
{
    common_init();
}

Wavefunction::Wavefunction(std::shared_ptr<Molecule> molecule,
                           std::shared_ptr<BasisSet> basis,
                           std::shared_ptr<BasisSet> ecpbasis) :
        options_(Process::environment.options), basisset_(basis), ecpbasisset_(ecpbasis), molecule_(molecule)
{
    common_init();
}

Wavefunction::Wavefunction(std::shared_ptr<Molecule> molecule, std::shared_ptr<BasisSet> basis) :
        options_(Process::environment.options), basisset_(basis), molecule_(molecule)
{

    common_init();
}

Wavefunction::Wavefunction(Options &options) :
        options_(options)
{
}

Wavefunction::~Wavefunction()
{
}

void Wavefunction::shallow_copy(SharedWavefunction other)
{
    shallow_copy(other.get());
}

void Wavefunction::shallow_copy(const Wavefunction *other)
{

    name_ = other->name_;
    basisset_ = other->basisset_;
    ecpbasisset_ = other->ecpbasisset_;
    basissets_ = other->basissets_;
    sobasisset_ = other->sobasisset_;
    AO2SO_ = other->AO2SO_;
    S_ = other->S_;
    molecule_ = other->molecule_;

    psio_ = other->psio_;
    integral_ = other->integral_;
    factory_ = other->factory_;
    memory_ = other->memory_;
    nalpha_ = other->nalpha_;
    nbeta_ = other->nbeta_;
    nfrzc_ = other->nfrzc_;

    print_ = other->print_;
    debug_ = other->debug_;
    density_fitted_ = other->density_fitted_;

    energy_ = other->energy_;
    efzc_ = other->efzc_;

    doccpi_ = other->doccpi_;
    soccpi_ = other->soccpi_;
    frzcpi_ = other->frzcpi_;
    frzvpi_ = other->frzvpi_;
    nalphapi_ = other->nalphapi_;
    nbetapi_ = other->nbetapi_;
    nsopi_ = other->nsopi_;
    nmopi_ = other->nmopi_;

    nso_ = other->nso_;
    nmo_ = other->nmo_;
    nirrep_ = other->nirrep_;

    oeprop_ = other->oeprop_;

    same_a_b_dens_ = other->same_a_b_dens_;
    same_a_b_orbs_ = other->same_a_b_orbs_;

    // Set by other classes
    H_ = other->H_;
    Ca_ = other->Ca_;
    Cb_ = other->Cb_;
    Da_ = other->Da_;
    Db_ = other->Db_;
    Fa_ = other->Fa_;
    Fb_ = other->Fb_;
    epsilon_a_ = other->epsilon_a_;
    epsilon_b_ = other->epsilon_b_;

    gradient_ = other->gradient_;
    hessian_ = other->hessian_;
    tpdm_gradient_contribution_ = other->tpdm_gradient_contribution_;
}

void Wavefunction::deep_copy(SharedWavefunction other)
{
    deep_copy(other.get());
}

void Wavefunction::deep_copy(const Wavefunction *other)
{
    if (!S_) {
        throw PSIEXCEPTION("Wavefunction::deep_copy must copy an initialized wavefunction.");
    }

    /// From typical constructor
    /// Some member data is not clone-able so we will copy
    name_ = other->name_;
    molecule_ = std::shared_ptr<Molecule>(new Molecule(other->molecule_->clone()));
    basisset_ = basisset_;
    ecpbasisset_ = other->ecpbasisset_;
    basissets_ = other->basissets_; // Still cannot copy basissets
    if(ecpbasisset_)
        integral_ = std::shared_ptr<IntegralFactory>(new IntegralFactory(basisset_, basisset_, basisset_, basisset_, ecpbasisset_));
    else
        integral_ = std::shared_ptr<IntegralFactory>(new IntegralFactory(basisset_, basisset_, basisset_, basisset_));
    sobasisset_ = std::shared_ptr<SOBasisSet>(new SOBasisSet(basisset_, integral_));
    factory_ = std::shared_ptr<MatrixFactory>(new MatrixFactory);
    factory_->init_with(other->nsopi_, other->nsopi_);
    AO2SO_ = other->AO2SO_->clone();
    S_ = other->S_->clone();

    psio_ = other->psio_; // We dont actually copy psio
    memory_ = other->memory_;
    nalpha_ = other->nalpha_;
    nbeta_ = other->nbeta_;
    nfrzc_ = other->nfrzc_;

    print_ = other->print_;
    debug_ = other->debug_;
    density_fitted_ = other->density_fitted_;

    energy_ = other->energy_;
    efzc_ = other->efzc_;
    variables_ = other->variables_;

    // How should we handle OEProp? oeprop_ = std::shared_ptr<OEProp>(this);

    doccpi_ = other->doccpi_;
    soccpi_ = other->soccpi_;
    frzcpi_ = other->frzcpi_;
    frzvpi_ = other->frzvpi_;
    nalphapi_ = other->nalphapi_;
    nbetapi_ = other->nbetapi_;
    nsopi_ = other->nsopi_;
    nmopi_ = other->nmopi_;

    nso_ = other->nso_;
    nmo_ = other->nmo_;
    nirrep_ = other->nirrep_;

    same_a_b_dens_ = other->same_a_b_dens_;
    same_a_b_orbs_ = other->same_a_b_orbs_;

    /// Below is not set in the typical constructor
    if (other->H_) H_ = other->H_->clone();
    if (other->Ca_) Ca_ = other->Ca_->clone();
    if (other->Cb_) Cb_ = other->Cb_->clone();
    if (other->Da_) Da_ = other->Da_->clone();
    if (other->Db_) Db_ = other->Db_->clone();
    if (other->Fa_) Fa_ = other->Fa_->clone();
    if (other->Fb_) Fb_ = other->Fb_->clone();
    if (other->epsilon_a_) epsilon_a_ = SharedVector(other->epsilon_a_->clone());
    if (other->epsilon_b_) epsilon_b_ = SharedVector(other->epsilon_b_->clone());

    if (other->gradient_) gradient_ = other->gradient_->clone();
    if (other->hessian_) hessian_ = other->hessian_->clone();
    if (other->tpdm_gradient_contribution_)
        tpdm_gradient_contribution_ = other->tpdm_gradient_contribution_->clone();

}

void Wavefunction::common_init()
{
    Wavefunction::initialize_singletons();
    if (!basisset_)
        throw PSIEXCEPTION("You can't initialize a Wavefunction that doesn't "
                                   "have a basis set");

    // Check the point group of the molecule. If it is not set, set it.
    if (!molecule_->point_group()) {
        molecule_->set_point_group(molecule_->find_point_group());
    }

    // Create an SO basis...we need the point group for this part.
    if(ecpbasisset_)
        integral_ = std::shared_ptr<IntegralFactory>(new IntegralFactory(basisset_, basisset_, basisset_, basisset_, ecpbasisset_));
    else
        integral_ = std::shared_ptr<IntegralFactory>(new IntegralFactory(basisset_, basisset_, basisset_, basisset_));
    sobasisset_ = std::shared_ptr<SOBasisSet>(new SOBasisSet(basisset_, integral_));

    std::shared_ptr<PetiteList> pet(new PetiteList(basisset_, integral_));
    AO2SO_ = pet->aotoso();

    // Obtain the dimension object to initialize the factory.
    nsopi_ = sobasisset_->dimension();
    nsopi_.set_name("SOs per irrep");

    factory_ = std::shared_ptr<MatrixFactory>(new MatrixFactory);
    factory_->init_with(nsopi_, nsopi_);

    nirrep_ = nsopi_.n();

    S_ = factory_->create_shared_matrix("S");
    std::shared_ptr<OneBodySOInt> Sint(integral_->so_overlap());
    Sint->compute(S_);

    // Initialize array that hold dimensionality information
    nmopi_ = Dimension(nirrep_, "MOs per irrep");
    nalphapi_ = Dimension(nirrep_, "Alpha electrons per irrep");
    nbetapi_ = Dimension(nirrep_, "Beta electrons per irrep");
    doccpi_ = Dimension(nirrep_, "Doubly occupied orbitals per irrep");
    soccpi_ = Dimension(nirrep_, "Singly occupied orbitals per irrep");
    frzcpi_ = Dimension(nirrep_, "Frozen core orbitals per irrep");
    frzvpi_ = Dimension(nirrep_, "Frozen virtual orbitals per irrep");

    // Obtain memory amount from the environment
    memory_ = Process::environment.get_memory();

    nso_ = basisset_->nbf();
    nmo_ = basisset_->nbf();
    for (int k = 0; k < nirrep_; k++) {
        nmopi_[k] = 0;
        doccpi_[k] = 0;
        soccpi_[k] = 0;
        nalphapi_[k] = 0;
        nbetapi_[k] = 0;
    }

    density_fitted_ = false;
    energy_ = 0.0;
    efzc_ = 0.0;
    same_a_b_dens_ = true;
    same_a_b_orbs_ = false;

    // Read in the debug flag
    debug_ = options_.get_int("DEBUG");
    print_ = options_.get_int("PRINT");

    // Determine the number of electrons in the system
    int nelectron = 0;
    for (int i = 0; i < molecule_->natom(); ++i)
        nelectron += (int) molecule_->Z(i);
    nelectron -= molecule_->molecular_charge();

    // If the user told us the multiplicity, read it from the input
    int multiplicity;
    if (molecule_->multiplicity_specified()) {
        multiplicity = molecule_->multiplicity();
    } else {
        if (nelectron % 2) {
            multiplicity = 2;
            molecule_->set_multiplicity(2);
            // There are an odd number of electrons
            outfile->Printf("    There are an odd number of electrons - assuming doublet.\n"
                                    "    Specify the multiplicity in the molecule input block.\n\n");
        } else {
            multiplicity = 1;
            molecule_->set_multiplicity(1);
            // There are an even number of electrons
            outfile->Printf("    There are an even number of electrons - assuming singlet.\n"
                                    "    Specify the multiplicity in the molecule input block.\n\n");
        }
    }

    // Make sure that the multiplicity is reasonable
    if (multiplicity - 1 > nelectron) {
        char *str = new char[100];
        sprintf(str, "There are not enough electrons for multiplicity = %d.\n"
                "Please check your input", multiplicity);
        throw SanityCheckError(str, __FILE__, __LINE__);
        delete[] str;
    }
    if (multiplicity % 2 == nelectron % 2) {
        char *str = new char[100];
        sprintf(str, "A multiplicity of %d with %d electrons is impossible.\n"
                        "Please check your input",
                multiplicity, nelectron);
        throw SanityCheckError(str, __FILE__, __LINE__);
        delete[] str;
    }

    nbeta_ = (nelectron - multiplicity + 1) / 2;
    nalpha_ = nbeta_ + multiplicity - 1;
}

Dimension Wavefunction::map_irreps(const Dimension& dimpi)
{
    std::shared_ptr<PointGroup> full = Process::environment.parent_symmetry();
    // If the parent symmetry hasn't been set, no displacements have been made
    if (!full) return dimpi;
    std::shared_ptr<PointGroup> sub = molecule_->point_group();

    // If the point group between the full and sub are the same return
    if (full->symbol() == sub->symbol())
        return dimpi;

    // Build the correlation table between full, and subgroup
    CorrelationTable corrtab(full, sub);
    Dimension mapped_dimpi(sub->char_table().nirrep());
    for (int h = 0; h < full->char_table().nirrep(); ++h) {
        int target = corrtab.gamma(h, 0);
        mapped_dimpi[target] += dimpi[h];
    }

    return mapped_dimpi;
}

void Wavefunction::initialize_singletons()
{
    static bool done = false;

    if (done)
        return;

    ioff[0] = 0;
    for (size_t i = 1; i < MAX_IOFF; ++i)
        ioff[i] = ioff[i - 1] + i;

    df[0] = 1.0;
    df[1] = 1.0;
    df[2] = 1.0;
    for (int i = 3; i < MAX_DF; ++i)
        df[i] = (i - 1) * df[i - 2];

    for (int i = 0; i < MAX_BC; ++i)
        for (int j = 0; j <= i; ++j)
            bc[i][j] = combinations(i, j);

    fac[0] = 1.0;
    for (int i = 1; i < MAX_FAC; ++i)
        fac[i] = i * fac[i - 1];

    done = true;
}

std::shared_ptr<Molecule> Wavefunction::molecule() const
{
    return molecule_;
}

std::shared_ptr<PSIO> Wavefunction::psio() const
{
    return psio_;
}

Options &Wavefunction::options() const
{
    return options_;
}

std::shared_ptr<IntegralFactory> Wavefunction::integral() const
{
    return integral_;
}

std::shared_ptr<BasisSet> Wavefunction::basisset() const
{
    return basisset_;
}

std::shared_ptr<BasisSet> Wavefunction::ecpbasisset() const
{
    return ecpbasisset_;
}

std::shared_ptr<BasisSet> Wavefunction::get_basisset(std::string label)
{
    // This may be slightly confusing, but better than changing this in 800 other places
    if (label == "ORBITAL"){
        return basisset_;
    } else if (basissets_.count(label) == 0){
        outfile->Printf("Could not find requested basisset (%s).", label.c_str());
        throw PSIEXCEPTION("Wavefunction::get_basisset: Requested basis set (" + label + ") was not set!\n");
    } else {
        return basissets_[label];
    }
}
void Wavefunction::set_basisset(std::string label, std::shared_ptr<BasisSet> basis)
{
    if (label == "ORBITAL"){
        throw PSIEXCEPTION("Cannot set the ORBITAL basis after the Wavefunction is built!");
    } else {
        basissets_[label] = basis;
    }
}

bool Wavefunction::basisset_exists(std::string label)
{
    if (basissets_.count(label) == 0)
        return false;
    else
        return true;
}

std::shared_ptr<SOBasisSet> Wavefunction::sobasisset() const
{
    return sobasisset_;
}

std::shared_ptr<MatrixFactory> Wavefunction::matrix_factory() const
{
    return factory_;
}

std::shared_ptr<Wavefunction> Wavefunction::reference_wavefunction() const
{
    return reference_wavefunction_;
}

void Wavefunction::set_reference_wavefunction(const std::shared_ptr<Wavefunction> wfn)
{
    reference_wavefunction_ = wfn;
}

SharedMatrix Wavefunction::Ca() const
{
    if (!Ca_) {
        if (!reference_wavefunction_)
            throw PSIEXCEPTION("Wavefunction::Ca: Unable to obtain MO coefficients.");
        else
            return reference_wavefunction_->Ca();
    }

    return Ca_;
}

SharedMatrix Wavefunction::Cb() const
{
    if (!Cb_) {
        if (!reference_wavefunction_)
            throw PSIEXCEPTION("Wavefunction::Cb: Unable to obtain MO coefficients.");
        else
            return reference_wavefunction_->Cb();
    }

    return Cb_;
}

std::vector <std::vector<int>> Wavefunction::subset_occupation(const Dimension &noccpi, const std::string &subset)
{
    if (!(subset == "FROZEN_OCC" ||
          subset == "FROZEN_VIR" ||
          subset == "ACTIVE_OCC" ||
          subset == "ACTIVE_VIR" ||
          subset == "FROZEN" ||
          subset == "ACTIVE" ||
          subset == "OCC" ||
          subset == "VIR" ||
          subset == "ALL"))
        throw PSIEXCEPTION("Orbital subset is not defined, should be ALL, OCC, VIR, FROZEN, ACTIVE, or look at doxygen");

    // Vector of relevent positions by irrep
    std::vector <std::vector<int>> positions;
    positions.resize(nirrep_);

    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < frzcpi_[h]; i++) {
            if (subset == "ALL" || subset == "FROZEN_OCC" || subset == "FROZEN" || subset == "OCC")
                positions[h].push_back(i);
        }
        for (int i = frzcpi_[h]; i < noccpi[h]; i++) {
            if (subset == "ALL" || subset == "ACTIVE_OCC" || subset == "ACTIVE" || subset == "OCC")
                positions[h].push_back(i);
        }
        for (int i = noccpi[h]; i < nmopi_[h] - frzvpi_[h]; i++) {
            if (subset == "ALL" || subset == "ACTIVE_VIR" || subset == "ACTIVE" || subset == "VIR")
                positions[h].push_back(i);
        }
        for (int i = nmopi_[h] - frzvpi_[h]; i < nmopi_[h]; i++) {
            if (subset == "ALL" || subset == "FROZEN_VIR" || subset == "FROZEN" || subset == "VIR")
                positions[h].push_back(i);
        }
    }
    return positions;
}

SharedMatrix Wavefunction::C_subset_helper(SharedMatrix C, const Dimension &noccpi, SharedVector epsilon, const std::string &basis, const std::string &subset)
{
    std::vector <std::vector<int>> positions = subset_occupation(noccpi, subset);

    Dimension nmopi(nirrep_);
    for (int h = 0; h < (int) positions.size(); h++) {
        nmopi[h] = positions[h].size();
    }
    SharedMatrix C2(new Matrix("C " + basis + " " + subset, nsopi_, nmopi));
    for (int h = 0; h < (int) positions.size(); h++) {
        for (int i = 0; i < (int) positions[h].size(); i++) {
            C_DCOPY(nsopi_[h], &C->pointer(h)[0][positions[h][i]], nmopi_[h], &C2->pointer(h)[0][i], nmopi[h]);
        }
    }

    if (basis == "AO") {

        SharedMatrix C3(new Matrix("C " + basis + " " + subset, nso_, nmopi.sum()));
        std::swap(C2, C3);

        std::vector <std::tuple<double, int, int>> order;
        for (int h = 0; h < nirrep_; h++) {
            for (int i = 0; i < (int) positions[h].size(); i++) {
                order.push_back(std::tuple<double, int, int>(epsilon->get(h, positions[h][i]), i, h));
            }
        }

        std::sort(order.begin(), order.end(), std::less < std::tuple < double, int, int > > ());

        for (int index = 0; index < (int) order.size(); index++) {
            int i = std::get<1>(order[index]);
            int h = std::get<2>(order[index]);

            int nao = nso_;
            int nso = nsopi_[h];

            if (!nso) continue;

            C_DGEMV('N', nao, nso, 1.0, AO2SO_->pointer(h)[0], nso, &C3->pointer(h)[0][i], nmopi[h], 0.0, &C2->pointer()[0][index], nmopi.sum());
        }

    } else if (basis == "SO" || basis == "MO") {
        // Already done
    } else {
        throw PSIEXCEPTION("Invalid basis requested, use AO, SO, or MO");
    }

    return C2;
}

SharedVector Wavefunction::epsilon_subset_helper(SharedVector epsilon, const Dimension &noccpi, const std::string &basis, const std::string &subset)
{
    std::vector <std::vector<int>> positions = subset_occupation(noccpi, subset);

    Dimension nmopi(nirrep_);
    for (int h = 0; h < (int) positions.size(); h++) {
        nmopi[h] = positions[h].size();
    }

    SharedVector C2;

    if (basis == "AO") {

        C2 = SharedVector(new Vector("Epsilon " + basis + " " + subset, nmopi.sum()));

        std::vector <std::tuple<double, int, int>> order;
        for (int h = 0; h < nirrep_; h++) {
            for (int i = 0; i < (int) positions[h].size(); i++) {
                order.push_back(std::tuple<double, int, int>(epsilon->get(h, positions[h][i]), i, h));
            }
        }

        std::sort(order.begin(), order.end(), std::less < std::tuple < double, int, int > > ());

        for (int index = 0; index < (int) order.size(); index++) {
            C2->set(0, index, std::get<0>(order[index]));
        }

    } else if (basis == "SO" || basis == "MO") {

        C2 = SharedVector(new Vector("Epsilon " + basis + " " + subset, nmopi));
        for (int h = 0; h < (int) positions.size(); h++) {
            for (int i = 0; i < (int) positions[h].size(); i++) {
                C2->set(h, i, epsilon->get(h, positions[h][i]));
            }
        }

    } else {
        throw PSIEXCEPTION("Invalid basis requested, use AO, SO, or MO");
    }

    return C2;
}

SharedMatrix Wavefunction::F_subset_helper(SharedMatrix F, SharedMatrix C, const std::string &basis)
{
    if (basis == "AO") {
        double *temp = new double[AO2SO_->max_ncol() * AO2SO_->max_nrow()];
        SharedMatrix F2 = SharedMatrix(new Matrix("Fock (AO basis)", basisset_->nbf(), basisset_->nbf()));
        int symm = F->symmetry();
        for (int h = 0; h < AO2SO_->nirrep(); ++h) {
            int nao = AO2SO_->rowspi()[0];
            int nsol = AO2SO_->colspi()[h];
            int nsor = AO2SO_->colspi()[h ^ symm];
            if (!nsol || !nsor) continue;
            double **Ulp = AO2SO_->pointer(h);
            double **Urp = AO2SO_->pointer(h ^ symm);
            double **FSOp = F->pointer(h ^ symm);
            double **FAOp = F2->pointer();
            C_DGEMM('N', 'T', nsol, nao, nsor, 1.0, FSOp[0], nsor, Urp[0], nsor, 0.0, temp, nao);
            C_DGEMM('N', 'N', nao, nao, nsol, 1.0, Ulp[0], nsol, temp, nao, 1.0, FAOp[0], nao);
        }
        delete[] temp;
        return F2;
    } else if (basis == "SO") {
        return SharedMatrix(F->clone());
    } else if (basis == "MO") {
        SharedMatrix F2(new Matrix("Fock (MO Basis)", C->colspi(), C->colspi()));

        int symm = F->symmetry();
        int nirrep = F->nirrep();

        double *SC = new double[C->max_ncol() * C->max_nrow()];
        double *temp = new double[C->max_ncol() * C->max_nrow()];
        for (int h = 0; h < nirrep; h++) {
            int nmol = C->colspi()[h];
            int nmor = C->colspi()[h ^ symm];
            int nsol = C->rowspi()[h];
            int nsor = C->rowspi()[h ^ symm];
            if (!nmol || !nmor || !nsol || !nsor) continue;
            double **Slp = S_->pointer(h);
            double **Srp = S_->pointer(h ^ symm);
            double **Clp = C->pointer(h);
            double **Crp = C->pointer(h ^ symm);
            double **Fmop = F2->pointer(h);
            double **Fsop = F->pointer(h);

            C_DGEMM('N', 'N', nsor, nmor, nsor, 1.0, Srp[0], nsor, Crp[0], nmor, 0.0, SC, nmor);
            C_DGEMM('N', 'N', nsol, nmor, nsor, 1.0, Fsop[0], nsor, SC, nmor, 0.0, temp, nmor);
            C_DGEMM('N', 'N', nsol, nmol, nsol, 1.0, Slp[0], nsol, Clp[0], nmol, 0.0, SC, nmol);
            C_DGEMM('T', 'N', nmol, nmor, nsol, 1.0, SC, nmol, temp, nmor, 0.0, Fmop[0], nmor);
        }
        delete[] temp;
        delete[] SC;
        return F2;
    } else {
        throw PSIEXCEPTION("Invalid basis requested, use AO, SO, or MO");
    }
}


SharedMatrix Wavefunction::D_subset_helper(SharedMatrix D, SharedMatrix C, const std::string &basis)
{
    if (basis == "AO") {
        double *temp = new double[AO2SO_->max_ncol() * AO2SO_->max_nrow()];
        SharedMatrix D2 = SharedMatrix(new Matrix("D (AO basis)", basisset_->nbf(), basisset_->nbf()));
        int symm = D->symmetry();
        for (int h = 0; h < AO2SO_->nirrep(); ++h) {
            int nao = AO2SO_->rowspi()[0];
            int nsol = AO2SO_->colspi()[h];
            int nsor = AO2SO_->colspi()[h ^ symm];
            if (!nsol || !nsor) continue;
            double **Ulp = AO2SO_->pointer(h);
            double **Urp = AO2SO_->pointer(h ^ symm);
            double **DSOp = D->pointer(h ^ symm);
            double **DAOp = D2->pointer();
            C_DGEMM('N', 'T', nsol, nao, nsor, 1.0, DSOp[0], nsor, Urp[0], nsor, 0.0, temp, nao);
            C_DGEMM('N', 'N', nao, nao, nsol, 1.0, Ulp[0], nsol, temp, nao, 1.0, DAOp[0], nao);
        }
        delete[] temp;
        return D2;
    } else if (basis == "CartAO") {
        /*
         * Added by ACS. Rob's definition of AO is simply a desymmetrized SO (i.e. using spherical basis
         * functions).  In cases like EFP and PCM where many OE integral evaluations are needed, we want
         * to avoid the spherical transformation, so we need to back transform the density matrix all the
         * way back to Cartesian AOs.
         */

        PetiteList petite(basisset_, integral_, true);
        SharedMatrix my_aotoso = petite.aotoso();
        double *temp = new double[my_aotoso->max_ncol() * my_aotoso->max_nrow()];
        SharedMatrix D2 = SharedMatrix(new Matrix("D (ao basis)", basisset_->nao(), basisset_->nao()));
        int symm = D->symmetry();
        for (int h = 0; h < my_aotoso->nirrep(); ++h) {
            int nao = my_aotoso->rowspi()[0];
            int nsol = my_aotoso->colspi()[h];
            int nsor = my_aotoso->colspi()[h ^ symm];
            if (!nsol || !nsor) continue;
            double **Ulp = my_aotoso->pointer(h);
            double **Urp = my_aotoso->pointer(h ^ symm);
            double **DSOp = D->pointer(h ^ symm);
            double **DAOp = D2->pointer();
            C_DGEMM('N', 'T', nsol, nao, nsor, 1.0, DSOp[0], nsor, Urp[0], nsor, 0.0, temp, nao);
            C_DGEMM('N', 'N', nao, nao, nsol, 1.0, Ulp[0], nsol, temp, nao, 1.0, DAOp[0], nao);
        }
        delete[] temp;
        return D2;
    } else if (basis == "SO") {
        return SharedMatrix(D->clone());
    } else if (basis == "MO") {
        SharedMatrix D2(new Matrix("D (MO Basis)", C->colspi(), C->colspi()));

        int symm = D->symmetry();
        int nirrep = D->nirrep();

        double *SC = new double[C->max_ncol() * C->max_nrow()];
        double *temp = new double[C->max_ncol() * C->max_nrow()];
        for (int h = 0; h < nirrep; h++) {
            int nmol = C->colspi()[h];
            int nmor = C->colspi()[h ^ symm];
            int nsol = C->rowspi()[h];
            int nsor = C->rowspi()[h ^ symm];
            if (!nmol || !nmor || !nsol || !nsor) continue;
            double **Slp = S_->pointer(h);
            double **Srp = S_->pointer(h ^ symm);
            double **Clp = C->pointer(h);
            double **Crp = C->pointer(h ^ symm);
            double **Dmop = D2->pointer(h);
            double **Dsop = D->pointer(h);

            C_DGEMM('N', 'N', nsor, nmor, nsor, 1.0, Srp[0], nsor, Crp[0], nmor, 0.0, SC, nmor);
            C_DGEMM('N', 'N', nsol, nmor, nsor, 1.0, Dsop[0], nsor, SC, nmor, 0.0, temp, nmor);
            C_DGEMM('N', 'N', nsol, nmol, nsol, 1.0, Slp[0], nsol, Clp[0], nmol, 0.0, SC, nmol);
            C_DGEMM('T', 'N', nmol, nmor, nsol, 1.0, SC, nmol, temp, nmor, 0.0, Dmop[0], nmor);
        }
        delete[] temp;
        delete[] SC;
        return D2;
    } else {
        throw PSIEXCEPTION("Invalid basis requested, use AO, CartAO, SO, or MO");
    }
}
SharedMatrix Wavefunction::basis_projection(SharedMatrix C_A, Dimension noccpi,
                                            std::shared_ptr<BasisSet> old_basis,
                                            std::shared_ptr<BasisSet> new_basis) {

    //Based on Werner's method from Mol. Phys. 102, 21-22, 2311
    std::shared_ptr<IntegralFactory> newfactory(new IntegralFactory(new_basis,new_basis,new_basis,new_basis));
    std::shared_ptr<IntegralFactory> hybfactory(new IntegralFactory(old_basis,new_basis,old_basis,new_basis));
    std::shared_ptr<OneBodySOInt> intBB(newfactory->so_overlap());
    std::shared_ptr<OneBodySOInt> intAB(hybfactory->so_overlap());

    std::shared_ptr<PetiteList> pet(new PetiteList(new_basis, newfactory));
    SharedMatrix AO2USO(pet->aotoso());

    SharedMatrix SAB(new Matrix("S_AB", C_A->nirrep(), C_A->rowspi(), AO2USO->colspi()));
    SharedMatrix SBB(new Matrix("S_BB", C_A->nirrep(), AO2USO->colspi(), AO2USO->colspi()));

    intAB->compute(SAB);
    intBB->compute(SBB);

    // C_A->print();
    // noccpi.print();
    // SAB->print();
    // SBB->print();

    newfactory.reset();
    hybfactory.reset();
    intAB.reset();
    intBB.reset();
    pet.reset();

    // Constrained to the same symmetry at the moment, we can relax this soon
    SharedMatrix C_B(new Matrix("C_B", C_A->nirrep(), AO2USO->colspi(), noccpi));

    // Block over irreps (soon united irreps)
    for (int h = 0; h < C_A->nirrep(); h++) {

        int nocc = noccpi[h];
        int na = C_A->rowspi()[h];
        int nb = AO2USO->colspi()[h];
        int nc = C_A->colspi()[h];

        if (nocc == 0 || na == 0 || nb == 0) continue;

        double** Ca = C_A->pointer(h);
        double** Cb = C_B->pointer(h);
        double** Sab = SAB->pointer(h);
        double** Sbb = SBB->pointer(h);

        int CholError = C_DPOTRF('L',nb,Sbb[0],nb);
        if (CholError !=0 )
            throw std::domain_error("S_BB Matrix Cholesky failed!");

        //Inversion (in place)
        int IError = C_DPOTRI('L',nb,Sbb[0],nb);
        if (IError !=0 )
            throw std::domain_error("S_BB Inversion Failed!");

        //LAPACK is smart and all, only uses half of the thing
        for (int m = 0; m<nb; m++)
            for (int n = 0; n<m; n++)
                Sbb[m][n] = Sbb[n][m];

        //Form T
        double** Temp1 = block_matrix(nb,nocc);
        C_DGEMM('T','N',nb,nocc,na,1.0,Sab[0],nb,Ca[0],nc,0.0,Temp1[0],nocc);

        // outfile->Printf(" Temp1:\n");
        // print_mat(Temp1,nb,nocc,"outfile");

        double** Temp2 = block_matrix(nb,nocc);
        C_DGEMM('N','N',nb,nocc,nb,1.0,Sbb[0],nb,Temp1[0],nocc,0.0,Temp2[0],nocc);

        //outfile->Printf(" Temp2:\n");
        //print_mat(Temp2,nb,nocc,outfile);

        double** Temp3 = block_matrix(na,nocc);
        C_DGEMM('N','N',na,nocc,nb,1.0,Sab[0],nb,Temp2[0],nocc,0.0,Temp3[0],nocc);

        //outfile->Printf(" Temp3:\n");
        //print_mat(Temp3,na,nocc,outfile);

        double** T = block_matrix(nocc,nocc);
        C_DGEMM('T','N',nocc,nocc,na,1.0,Ca[0],nc,Temp3[0],nocc,0.0,T[0],nocc);

        // outfile->Printf(" T:\n");
        // print_mat(T,nocc,nocc,"outfile");

        //Find T^-1/2
        // First, diagonalize T
        // the C_DSYEV call replaces the original matrix T with its eigenvectors
        double* eigval = init_array(nocc);
        int lwork = nocc * 3;
        double* work = init_array(lwork);
        int stat = C_DSYEV('v','u',nocc,T[0],nocc,eigval, work,lwork);
        if (stat != 0) {
            outfile->Printf( "C_DSYEV failed\n");
            exit(PSI_RETURN_FAILURE);
        }
        free(work);

        // Now T contains the eigenvectors of the original T
        // Copy T to T_copy
        double **T_mhalf = block_matrix(nocc, nocc);
        double **T_copy = block_matrix(nocc, nocc);
        C_DCOPY(nocc*nocc,T[0],1,T_copy[0],1);

        // Now form T^{-1/2} = U(T)*t^{-1/2}*U,
        // where t^{-1/2} is the diagonal matrix of the inverse square roots
        // of the eigenvalues, and U is the matrix of eigenvectors of T
        for (int i=0; i<nocc; i++) {
            if (eigval[i] < 1.0E-10)
                eigval[i] = 0.0;
            else
                eigval[i] = 1.0 / sqrt(eigval[i]);

            // scale one set of eigenvectors by the diagonal elements t^{-1/2}
            C_DSCAL(nocc, eigval[i], T[i], 1);
        }
        free(eigval);

        // T_mhalf = T_copy(T) * T
        C_DGEMM('t','n',nocc,nocc,nocc,1.0,
                T_copy[0],nocc,T[0],nocc,0.0,T_mhalf[0],nocc);

        //Form CB
        C_DGEMM('N','N',nb,nocc,nocc,1.0,Temp2[0],nocc,T_mhalf[0],nocc,0.0,Cb[0],nocc);

        free_block(Temp1);
        free_block(Temp2);
        free_block(Temp3);
        free_block(T);
        free_block(T_copy);
        free_block(T_mhalf);

    }
    return C_B;
}

OrbitalSpace Wavefunction::alpha_orbital_space(const std::string &id, const std::string &basis, const std::string &subset)
{
    return OrbitalSpace(id, subset, Ca_subset(basis, subset), epsilon_a_subset(basis, subset), basisset_, integral_);
}

OrbitalSpace Wavefunction::beta_orbital_space(const std::string &id, const std::string &basis, const std::string &subset)
{
    return OrbitalSpace(id, subset, Cb_subset(basis, subset), epsilon_b_subset(basis, subset), basisset_, integral_);
}

SharedMatrix Wavefunction::Ca_subset(const std::string &basis, const std::string &subset)
{
    return C_subset_helper(Ca_, nalphapi_, epsilon_a_, basis, subset);
}

SharedMatrix Wavefunction::Cb_subset(const std::string &basis, const std::string &subset)
{
    return C_subset_helper(Cb_, nbetapi_, epsilon_b_, basis, subset);
}

SharedMatrix Wavefunction::Da_subset(const std::string &basis)
{
    return D_subset_helper(Da_, Ca_, basis);
}

SharedMatrix Wavefunction::Db_subset(const std::string &basis)
{
    return D_subset_helper(Db_, Cb_, basis);
}

SharedVector Wavefunction::epsilon_a_subset(const std::string &basis, const std::string &subset)
{
    return epsilon_subset_helper(epsilon_a_, nalphapi_, basis, subset);
}

SharedVector Wavefunction::epsilon_b_subset(const std::string &basis, const std::string &subset)
{
    return epsilon_subset_helper(epsilon_b_, nbetapi_, basis, subset);
}

SharedMatrix Wavefunction::Fa() const
{
    return Fa_;
}

SharedMatrix Wavefunction::Fb() const
{
    return Fb_;
}

std::shared_ptr<Matrix> Wavefunction::Lagrangian() const
{
    return Lagrangian_;
}

std::shared_ptr<Matrix> Wavefunction::tpdm_gradient_contribution() const
{
    throw PSIEXCEPTION("This type of wavefunction has not defined a TPDM gradient contribution!");
}

std::shared_ptr<Vector> Wavefunction::epsilon_a() const
{
    return epsilon_a_;
}

std::shared_ptr<Vector> Wavefunction::epsilon_b() const
{
    return epsilon_b_;
}

const SharedMatrix Wavefunction::Da() const
{
    return Da_;
}

SharedMatrix Wavefunction::Db() const
{
    return Db_;
}

SharedMatrix Wavefunction::X() const
{
    return Lagrangian_;
}

SharedMatrix Wavefunction::gradient() const
{
    return gradient_;
}

void Wavefunction::set_gradient(SharedMatrix &grad)
{
    gradient_ = grad;
}

SharedMatrix Wavefunction::hessian() const
{
    return hessian_;
}

void Wavefunction::set_hessian(SharedMatrix &hess)
{
    hessian_ = hess;
}

std::shared_ptr<Vector> Wavefunction::frequencies() const
{
    return frequencies_;
}

std::shared_ptr<Vector> Wavefunction::normalmodes() const
{
    return normalmodes_;
}

void Wavefunction::set_frequencies(std::shared_ptr<Vector> &freqs)
{
    frequencies_ = freqs;
}

void Wavefunction::set_normalmodes(std::shared_ptr<Vector> &norms)
{
    normalmodes_ = norms;
}

void Wavefunction::save() const
{
}

std::shared_ptr<Vector> Wavefunction::get_atomic_point_charges() const
{
    std::shared_ptr<std::vector<double>> q = atomic_point_charges();

    int n = molecule_->natom();
    std::shared_ptr<Vector> q_vector(new Vector(n));
    for (int i = 0; i < n; ++i)
        q_vector->set(i, (*q)[i]);
    return q_vector;
}
double Wavefunction::get_variable(std::string label)
{
    std::string uc_label = label;
    std::transform(uc_label.begin(), uc_label.end(), uc_label.begin(), ::toupper);

    if (variables_.count(uc_label) == 0){
        throw PSIEXCEPTION("Wavefunction::get_variable: Requested variable was not set!\n");
    } else {
        return variables_[uc_label];
    }
}
SharedMatrix Wavefunction::get_array(std::string label)
{
    if (arrays_.count(label) == 0){
        throw PSIEXCEPTION("Wavefunction::get_array: Requested array was not set!\n");
    } else {
        return arrays_[label];
    }
}

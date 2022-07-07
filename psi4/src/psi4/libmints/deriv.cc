/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
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

/*
 * deriv.cc
 *
 *  Created on: Feb 24, 2009
 *      Author: jturney
 */

#include "psi4/libmints/sointegral_twobody.h"
#include "psi4/libmints/deriv.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/factory.h"
#include "psi4/libmints/sointegral_onebody.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"

#include <cstdio>
#include <cstdlib>
#include <string>
#include <memory>
#include <array>

namespace psi {

size_t counter;

class PSI_API CorrelatedFunctor {
    /// The buffer to hold the TPDM
    double *tpdm_buffer_;
    /// Pointer to the current TPDM element
    double *tpdm_ptr_;
    /// How large the buffer is, for each shell pair
    size_t *buffer_sizes_;
    /// The PSIO object to use for disk I/O
    std::shared_ptr<PSIO> psio_;

   public:
    int nthread;
    std::vector<SharedVector> result;

    CorrelatedFunctor() {
        throw PSIEXCEPTION("CorrelatedRestrictedFunctor(): Default constructor called. This shouldn't happen.");
    }
    CorrelatedFunctor(SharedVector results) : psio_(_default_psio_lib_) {
        nthread = Process::environment.get_n_threads();
        result.push_back(results);
        for (int i = 1; i < nthread; ++i) result.push_back(std::make_shared<Vector>(std::move(result[0]->clone())));
        size_t num_pairs = 0;
        psio_->read_entry(PSIF_AO_TPDM, "Num. Pairs", (char *)&num_pairs, sizeof(size_t));
        buffer_sizes_ = new size_t[num_pairs];
        psio_->read_entry(PSIF_AO_TPDM, "TPDM Buffer Sizes", (char *)buffer_sizes_, num_pairs * sizeof(size_t));
        size_t max_size = 0;
        for (size_t i = 0; i < num_pairs; ++i) max_size = max_size > buffer_sizes_[i] ? max_size : buffer_sizes_[i];
        tpdm_buffer_ = new double[max_size];
        tpdm_ptr_ = tpdm_buffer_;
    }

    void finalize() {
        // Do summation over threads
        for (int i = 1; i < nthread; ++i) {
            result[0]->add(*result[i]);
        }
        delete[] tpdm_buffer_;
        delete[] buffer_sizes_;
    }

    void load_tpdm(size_t id) {
        // TODO, make this work with threads (each thread needs its own buffer)
        auto *toc = new char[40];
        sprintf(toc, "SO_TPDM_FOR_PAIR_%zd", id);
        size_t buffer_size = buffer_sizes_[id];
        psio_->read_entry(PSIF_AO_TPDM, toc, (char *)tpdm_buffer_, buffer_size * sizeof(double));
        delete[] toc;
        tpdm_ptr_ = tpdm_buffer_;
    }

    void next_tpdm_element() { ++tpdm_ptr_; }

    void operator()(int salc, int pabs, int qabs, int rabs, int sabs, int /*pirrep*/, int /*pso*/, int /*qirrep*/,
                    int /*qso*/, int /*rirrep*/, int /*rso*/, int /*sirrep*/, int /*sso*/, double value) {
        int thread = 0;
        // This was the old line, but it always returned 0 anyways...
        // WorldComm->thread_id(pthread_self());

        double prefactor = 8.0;
        if (pabs == qabs) prefactor *= 0.5;
        if (rabs == sabs) prefactor *= 0.5;
        if (pabs == rabs && qabs == sabs) prefactor *= 0.5;
        result[thread]->add(salc, prefactor * (*tpdm_ptr_) * value);
    }
};

class PSI_API ScfRestrictedFunctor {
    SharedMatrix D_;

   public:
    int nthread;
    std::vector<SharedVector> result;

    ScfRestrictedFunctor() {
        throw PSIEXCEPTION("ScfRestrictedFunctor(): Default constructor called. This shouldn't happen.");
    }

    // Added for debugging MADNESS.
    //    ScfRestrictedFunctor(const ScfRestrictedFunctor&) {
    //        throw PSIEXCEPTION("ScfRestrictedFunctor(): Copy constructor called.\n");
    //    }

    //    ScfRestrictedFunctor& operator=(const ScfRestrictedFunctor&) {
    //        throw PSIEXCEPTION("ScfRestrictedFunctor(): Assignment operator called. This shouldn't happen.");
    //        return *this;
    //    }

    ScfRestrictedFunctor(SharedVector results, std::shared_ptr<Matrix> D) : D_(D) {
        counter = 0;
        nthread = Process::environment.get_n_threads();
        result.push_back(results);

        for (int i = 1; i < nthread; ++i) result.push_back(std::make_shared<Vector>(std::move(results->clone())));
    }
    ~ScfRestrictedFunctor() {}

    void finalize() {
        // Do summation over threads
        for (int i = 1; i < nthread; ++i) {
            result[0]->add(*result[i]);
        }
    }

    void load_tpdm(size_t /*id*/) {}
    void next_tpdm_element() {}

    void operator()(int salc, int pabs, int qabs, int rabs, int sabs, int pirrep, int pso, int qirrep, int qso,
                    int rirrep, int rso, int sirrep, int sso, double value) {
        int thread = 0;
        // Old call::WorldComm->thread_id(pthread_self());

        // Previously, we applied a factor of 4 after the fact...apply it from the beginning now.
        double prefactor = 4.0;

        if (pabs == qabs) prefactor *= 0.5;
        if (rabs == sabs) prefactor *= 0.5;
        if (pabs == rabs && qabs == sabs) prefactor *= 0.5;

        double four_index_D = 0.0;

        if (pirrep == qirrep && rirrep == sirrep)
            four_index_D = 4.0 * D_->get(pirrep, pso, qso) * D_->get(rirrep, rso, sso);
        if (pirrep == rirrep && qirrep == sirrep) four_index_D -= D_->get(pirrep, pso, rso) * D_->get(qirrep, qso, sso);
        if (pirrep == sirrep && qirrep == rirrep) four_index_D -= D_->get(pirrep, pso, sso) * D_->get(qirrep, qso, rso);

        four_index_D *= prefactor;

        result[thread]->add(salc, four_index_D * value);
        counter++;
    }
};

class PSI_API ScfAndDfCorrelationRestrictedFunctor {
    SharedMatrix D_ref_;
    SharedMatrix D_;
    ScfRestrictedFunctor scf_functor_;
    std::vector<SharedVector> result_vec_;
    SharedVector results_;

   public:
    int nthread;

    //    ScfAndDfCorrelationRestrictedFunctor() {
    //        throw PSIEXCEPTION("SCFAndDFCorrelationRestrictedFunctor(): Default constructor called. This shouldn't
    //        happen.");
    //    }

    ScfAndDfCorrelationRestrictedFunctor(SharedVector results, ScfRestrictedFunctor &scf_functor,
                                         std::shared_ptr<Matrix> D, std::shared_ptr<Matrix> D_ref)
        : D_ref_(D_ref), D_(D), scf_functor_(scf_functor), results_(results) {
        counter = 0;
        nthread = Process::environment.get_n_threads();
        result_vec_.push_back(results);

        for (int i = 1; i < nthread; ++i) result_vec_.push_back(std::make_shared<Vector>(std::move(results->clone())));
    }

    ScfAndDfCorrelationRestrictedFunctor() {
        throw PSIEXCEPTION(
            "ScfAndDfCorrelationRestrictedFunctor(): Default constructor called. This shouldn't happen.");
    }

    ~ScfAndDfCorrelationRestrictedFunctor() {}

    void load_tpdm(size_t /*id*/) {}
    void next_tpdm_element() {}

    void finalize() {
        // Make sure the SCF code is done
        scf_functor_.finalize();
        // Do summation over threads
        for (int i = 1; i < nthread; ++i) {
            result_vec_[0]->add(*result_vec_[i]);
        }
    }

    void operator()(int salc, int pabs, int qabs, int rabs, int sabs, int pirrep, int pso, int qirrep, int qso,
                    int rirrep, int rso, int sirrep, int sso, double value) {
        int thread = 0;
        // Old call WorldComm->thread_id(pthread_self());

        bool braket = pabs != rabs || qabs != sabs;
        bool bra = pabs != qabs;
        bool ket = rabs != sabs;

        double four_index_D = 0.0;

        double Coulomb1 = 0.0;
        double Coulomb2 = 0.0;
        double Exchange1 = 0.0;
        if (pirrep == qirrep && rirrep == sirrep) {
            Coulomb1 = 2.0 * D_->get(pirrep, pso, qso) * D_ref_->get(rirrep, rso, sso);
            Coulomb2 = 2.0 * D_->get(rirrep, rso, sso) * D_ref_->get(pirrep, pso, qso);
        }
        if (pirrep == rirrep && qirrep == sirrep) Exchange1 = D_->get(pirrep, pso, rso) * D_ref_->get(qirrep, qso, sso);
        // (pq|rs)
        four_index_D = Coulomb1 - Exchange1;

        if (bra && ket && braket) {
            four_index_D += 3 * Coulomb1;
            four_index_D += 4 * Coulomb2;
            // (qp|rs) and (rs|qp)
            if (qirrep == rirrep && pirrep == sirrep)
                four_index_D -= 2.0 * D_->get(qirrep, qso, rso) * D_ref_->get(pirrep, pso, sso);
            // (pq|sr) and (sr|pq)
            if (pirrep == sirrep && qirrep == rirrep)
                four_index_D -= 2.0 * D_->get(pirrep, pso, sso) * D_ref_->get(qirrep, qso, rso);
            // (qp|sr) and (sr|qp)
            if (qirrep == sirrep && pirrep == rirrep)
                four_index_D -= 2.0 * D_->get(qirrep, qso, sso) * D_ref_->get(pirrep, pso, rso);
            // (rs|pq)
            four_index_D -= Exchange1;
        } else if (bra && ket) {
            four_index_D += 3 * Coulomb1;
            // (qp|rs)
            if (qirrep == rirrep && pirrep == sirrep)
                four_index_D -= D_->get(qirrep, qso, rso) * D_ref_->get(pirrep, pso, sso);
            // (pq|sr)
            if (pirrep == sirrep && qirrep == rirrep)
                four_index_D -= D_->get(pirrep, pso, sso) * D_ref_->get(qirrep, qso, rso);
            // (qp|sr)
            if (qirrep == sirrep && pirrep == rirrep)
                four_index_D -= D_->get(qirrep, qso, sso) * D_ref_->get(pirrep, pso, rso);
        } else if (bra) {
            four_index_D += Coulomb1;
            four_index_D += 2 * Coulomb2;
            // (qp|rs)
            if (qirrep == rirrep && pirrep == sirrep)
                four_index_D -= D_->get(qirrep, qso, rso) * D_ref_->get(pirrep, pso, sso);
            // (rs|pq)
            four_index_D -= Exchange1;
            // (rs|qp)
            if (rirrep == qirrep && sirrep == pirrep)
                four_index_D -= D_->get(rirrep, rso, qso) * D_ref_->get(sirrep, sso, pso);
        } else if (ket) {
            four_index_D += Coulomb1;
            four_index_D += 2 * Coulomb2;
            // (pq|sr)
            if (pirrep == sirrep && qirrep == rirrep)
                four_index_D -= D_->get(pirrep, pso, sso) * D_ref_->get(qirrep, qso, rso);
            // (rs|pq)
            four_index_D -= Exchange1;
            // (sr|qp)
            if (sirrep == qirrep && rirrep == pirrep)
                four_index_D -= D_->get(sirrep, sso, qso) * D_ref_->get(rirrep, rso, pso);
        } else if (braket) {
            four_index_D += Coulomb2;
            // (rs|pq)
            four_index_D -= Exchange1;
        }

        result_vec_[thread]->add(salc, four_index_D * value);

        // Make sure the SCF contribution is computed.
        scf_functor_(salc, pabs, qabs, rabs, sabs, pirrep, pso, qirrep, qso, rirrep, rso, sirrep, sso, value);
        counter++;
    }
};

class PSI_API ScfUnrestrictedFunctor {
    SharedMatrix Da_;
    SharedMatrix Db_;

   public:
    int nthread;
    std::vector<SharedVector> result;

    ScfUnrestrictedFunctor() { throw PSIEXCEPTION("ScfUnrestrictedFunctor(): Oh come on!!!"); }

    ScfUnrestrictedFunctor(SharedVector results, std::shared_ptr<Matrix> Da, std::shared_ptr<Matrix> Db)
        : Da_(Da), Db_(Db) {
        nthread = Process::environment.get_n_threads();
        result.push_back(results);
        for (int i = 1; i < nthread; ++i) result.push_back(std::make_shared<Vector>(std::move(result[0]->clone())));
    }
    ~ScfUnrestrictedFunctor() {}

    void load_tpdm(size_t /*id*/) {}
    void next_tpdm_element() {}

    void finalize() {
        // Do summation over threads
        for (int i = 1; i < nthread; ++i) result[0]->add(*result[i]);
    }

    void operator()(int salc, int pabs, int qabs, int rabs, int sabs, int pirrep, int pso, int qirrep, int qso,
                    int rirrep, int rso, int sirrep, int sso, double value) {
        int thread = 0;
        // Old call: WorldComm->thread_id(pthread_self());
        double prefactor = 1.0;

        if (pabs == qabs) prefactor *= 0.5;
        if (rabs == sabs) prefactor *= 0.5;
        if (pabs == rabs && qabs == sabs) prefactor *= 0.5;

        double four_index_D = 0.0;

        if (pirrep == qirrep && rirrep == sirrep) {
            four_index_D = 4.0 * (Da_->get(pirrep, pso, qso) + Db_->get(pirrep, pso, qso)) *
                           (Da_->get(rirrep, rso, sso) + Db_->get(rirrep, rso, sso));
        }
        if (pirrep == rirrep && qirrep == sirrep) {
            four_index_D -= 2.0 * ((Da_->get(pirrep, pso, rso) * Da_->get(qirrep, qso, sso)) +
                                   (Db_->get(pirrep, pso, rso) * Db_->get(qirrep, qso, sso)));
        }
        if (pirrep == sirrep && qirrep == rirrep) {
            four_index_D -= 2.0 * ((Da_->get(pirrep, pso, sso) * Da_->get(rirrep, rso, qso)) +
                                   (Db_->get(pirrep, pso, sso) * Db_->get(rirrep, rso, qso)));
        }
        value *= prefactor;

        result[thread]->add(salc, four_index_D * value);
    }
};

Deriv::Deriv(const std::shared_ptr<Wavefunction> &wave, char needed_irreps, bool project_out_translations,
             bool project_out_rotations)
    : wfn_(wave), cdsalcs_(wave->molecule(), needed_irreps, project_out_translations, project_out_rotations) {
    integral_ = wave->integral();
    basis_ = wave->basisset();
    sobasis_ = wave->sobasisset();
    factory_ = wave->matrix_factory();
    molecule_ = wave->molecule();
    natom_ = molecule_->natom();
    tpdm_presorted_ = false;
    deriv_density_backtransformed_ = false;
    ignore_reference_ = false;

    // Results go here.
    tpdm_contr_ = factory_->create_shared_matrix("Two-electron contribution to gradient", natom_, 3);
    gradient_ = factory_->create_shared_matrix("Total gradient", natom_, 3);

    if (wfn_->options().get_int("PRINT") > 2) cdsalcs_.print();
}

/* Technical Documentation
 * This requires five intermediates.
 * 1. One-Particle Density Matrix
 *  Analagous to <Ψ|a^p_q|Ψ>. The (pq) element contracts against the integral derivative h^x_(pq).
 *  This quantity separates into alpha and beta blocks. Each block must be back-transformed, and then
 *  stored on the wavefunction as Da_ and Db_.
 * 2. Lagrangian, AKA, Energy-Weighted Density Matrix, AKA, Generalized Fock Matrix
 *  The (pq) element contracts against the overlap integral derivative -S^x_(pq).
 *  This quantity separates into alpha and beta blocks. Each block must be back-transformed, the blocks
 *  summed together, and the result stored on the wavefunction as Lagangian_.
 *  Note that this is NOT the behavior of the Lagrangian for RHF! As defined in Psi4, for Lagrangians
 *  for RHF wavefunctions, the blocks are not added together.
 * 3. Metric Density
 *  The (pq) element contracts against the metric integral derivative -J^x_(pq).
 *  Again, this quantity separates into alpha and beta blocks. Build the blocks, and add them
 *  together. Each auxiliary basis set has its own density. These must be stored in PSIF_AO_TPDM in
 *  LowerTriangular format under the names "Metric Reference Density" and "Metric Correlation Density".
 * 4. 3-Center Density
 *  The (Q|pq) element contracts against the 3-center integral derivative (Q|pq)^x.
 *  These must be stored in PSIF_AO_TPDM under the names "3-Center Reference Density" and
 *  "3-Center Correlation Density". First looping over the auxiliary index and then store
 *  the lower triangular block of the relevant slice of metric elements.
 * Note (Back-Transformation):
 *  An MO basis quantity is back-transformed to the AO index when each index is contracted against the C
 *  matrix to give a quantity in the AO basis. A back-transformed density matrix contracted against the AO
 *  basis derivative integrals gives the same result as the density matrix contracted against the MO basis
 *  derivative integrals. This single back-transformation allows us to avoid back-transforming the derivative
 *  integrals for each perturbation.
 * Note (Densities):
 *  Analytic gradient theory can always be formulated as involving a derivative of the two-electron integrals
 *  against contracted some other quantity, known as the (relaxed) two-particle density matrix. In the density-fitting
 *  approximation, inserting the expression for the two-electron intermediates and applying the chain rules shows that
 *  with suitably defined intermediates, all derivative integral contractions are against two or three index quantities.
 *  This is the one and only difference between conventional integral and density-fitted integral analytic gradient
 *  theory. But it's important, as the explicit, nao^4 two-particle density matrix reduces to a naux^2 and nao^2*naux
 *  intermediate.
 * Note (Basis Sets):
 *  It is standard for electronic structure methods to use one basis set for density-fitting the SCF reference and a
 *  different basis set for a correlated correction. In this case, the derivatives of both basis sets appear in the
 *  gradient formula, and each have their associated density matrix. Contracting a term against derivatives of integrals
 *  from the wrong auxiliary basis can lead to gradient errors on the order of 10^-5. When implementing density-fitted
 *  gradients, you must correctly assign which TEI terms in the conventional integral theory come from which basis set.
 *  For SCF reference theories, these terms always come from the orbital Z-vector or the Fock matrix. For non-SCF
 *  reference theories, proceed with extra caution.
 */
SharedMatrix Deriv::compute_df(const std::string& ref_aux_name, const std::string& cor_aux_name) {
    molecule_->print_in_bohr();
    
    if (!wfn_) throw("In Deriv: The wavefunction passed in is empty!");

    if (natom_ == 1) {
        // This is an atom...there is no gradient.
        outfile->Printf("    A single atom has no gradient.\n");
        // Save the gradient to the wavefunction so that optking can optimize with it
        wfn_->set_gradient(gradient_);
        return gradient_;
    }

    // The indirection of letting this be a pointer allows us to define tei_terms later.
    auto gradient_terms = std::make_shared<std::vector<SharedMatrix>>();

    // Obtain nuclear repulsion contribution from the wavefunction
    auto enuc = std::make_shared<Matrix>(molecule_->nuclear_repulsion_energy_deriv1(wfn_->get_dipole_field_strength()));
    gradient_terms->push_back(enuc);

    const auto& mints = wfn_->mintshelper();

    // One-electron derivatives.
    auto Da = wfn_->Da();
    auto Db = wfn_->Db();
    auto Dtot_AO = std::make_shared<Matrix>("AO basis total D", wfn_->nso(), wfn_->nso());
    auto Dtot = Da->clone();
    Dtot->add(Db);
    Dtot_AO->remove_symmetry(Dtot, wfn_->aotoso()->transpose());
    gradient_terms->push_back(mints->core_hamiltonian_grad(Dtot_AO));

    // Overlap derivatives
    auto s_deriv = cdsalcs_.create_matrices("S'", *factory_);
    auto s_int = integral_->so_overlap(1);
    s_int->compute_deriv1(s_deriv, cdsalcs_);
    auto X = wfn_->lagrangian();
    std::vector<double> Xcont;
    for (size_t cd = 0; cd < cdsalcs_.ncd(); ++cd) {
        Xcont.push_back(-X->vector_dot(s_deriv[cd]));
    }
    // B^t g_q^t = g_x^t -> g_q B = g_x
    // In einsum notation i,ijk -> jk
    auto st = cdsalcs_.matrix();
    auto B = st->pointer(0);
    auto *cart = new double[3 * natom_];
    C_DGEMM('n', 'n', 1, 3 * natom_, cdsalcs_.ncd(), 1.0, Xcont.data(), cdsalcs_.ncd(), B[0], 3 * natom_, 0.0, cart,
            3 * natom_);
    auto x_contr = factory_->create_shared_matrix("Lagrangian contribution to gradient", natom_, 3);
    for (int a = 0; a < natom_; ++a)
        for (int xyz = 0; xyz < 3; ++xyz) x_contr->set(a, xyz, cart[3 * a + xyz]);
    gradient_terms->push_back(x_contr);

    // DF TEI derivatives
    std::vector<std::pair<std::string, std::string>> aux_data{{ref_aux_name, "Reference"}, {cor_aux_name, "Correlation"}};
    _default_psio_lib_->open(PSIF_AO_TPDM, PSIO_OPEN_OLD);
    bool separate_tei = wfn_->options().get_int("PRINT") > 2;
    auto tei_terms = separate_tei ? gradient_terms : std::make_shared<std::vector<SharedMatrix>>();
    for (const auto& aux_datum: aux_data) {
        auto naux = wfn_->get_basisset(aux_datum.first)->nbf();
        auto metric_density = std::make_shared<Matrix>("Metric " + aux_datum.second + " Density", naux, naux);
        metric_density->load(_default_psio_lib_, PSIF_AO_TPDM, Matrix::SaveType::LowerTriangle);
        metric_density->scale(2);
        std::map<std::string, SharedMatrix> densities;
        densities["Metric " + aux_datum.second] = metric_density;
        auto results = mints->metric_grad(densities, aux_datum.first);
        for (const auto& kv: results) {
            tei_terms->push_back(kv.second);
        }
        auto result = mints->three_idx_grad(aux_datum.first, "3-Center " + aux_datum.second + " Density" , "3-Center " + aux_datum.second);
        tei_terms->push_back(result);
    }
    _default_psio_lib_->close(PSIF_AO_TPDM, 1); // 1 = keep contents of PSIF_AO_TPDM

    if (!separate_tei) {
        for (const auto& tei_term : *tei_terms) {
            tpdm_contr_->add(tei_term);
        }
        gradient_terms->push_back(tpdm_contr_);
    }

    for (auto gradient: *gradient_terms) {
        gradient->symmetrize_gradient(molecule_);
        if (wfn_->options().get_int("PRINT") > 1) gradient->print_atom_vector();
        gradient_->add(gradient);
    }

    // Print the atom vector
    gradient_->print_atom_vector();

    // Save the gradient to the wavefunction so that optking can optimize with it
    wfn_->set_gradient(gradient_);

    return gradient_;
}

SharedMatrix Deriv::compute(DerivCalcType deriv_calc_type) {
    molecule_->print_in_bohr();

    if (!wfn_) throw("In Deriv: The wavefunction passed in is empty!");

    if (natom_ == 1) {
        // This is an atom...there is no gradient.
        outfile->Printf("    A single atom has no gradient.\n");
        // Save the gradient to the wavefunction so that optking can optimize with it
        wfn_->set_gradient(gradient_);
        return gradient_;
    }

    // Initialize an ERI object requesting derivatives.
    int nthread = Process::environment.get_n_threads();
    std::vector<std::shared_ptr<TwoBodyAOInt> > ao_eri(nthread);
    ao_eri[0] = std::shared_ptr<TwoBodyAOInt>(integral_->eri(1));
    for (int i = 1; i < nthread; ++i)
        ao_eri[i] = std::shared_ptr<TwoBodyAOInt>(ao_eri.front()->clone());
    TwoBodySOInt so_eri(ao_eri, integral_, cdsalcs_);

    // A certain optimization can be used if we know we only need totally symmetric
    // derivatives.
    so_eri.set_only_totally_symmetric(true);

    // Compute one-electron derivatives.
    auto s_deriv = cdsalcs_.create_matrices("S'", *factory_);
    auto s_int = integral_->so_overlap(1);

    s_int->compute_deriv1(s_deriv, cdsalcs_);

    auto mints = std::make_shared<MintsHelper>(wfn_->basisset(), wfn_->options());
    int ncd = cdsalcs_.ncd();
    auto TPDMcont_vector = std::make_shared<Vector>(ncd);
    auto Xcont_vector = std::make_shared<Vector>(ncd);
    auto Dcont_vector = std::make_shared<Vector>(ncd);
    SharedVector TPDM_ref_cont_vector;
    SharedVector X_ref_cont_vector;
    auto Xcont = Xcont_vector->pointer();
    auto TPDMcont = TPDMcont_vector->pointer();
    double *TPDM_ref_cont = nullptr;
    double *X_ref_cont = nullptr;

    // Try and grab the OPDM and lagrangian from the wavefunction
    auto Da = wfn_->Da();
    auto Db = wfn_->Db();
    auto X = wfn_->lagrangian();

    // The current wavefunction's reference wavefunction, nullptr for SCF/DFT
    std::shared_ptr<Wavefunction> ref_wfn = wfn_->reference_wavefunction();
    // Whether the SCF contribution is separate from the correlated terms
    // This is currently a hack and should be improved.
    bool reference_separate = (X) && ref_wfn && wfn_->module() != "dct";

    // If the function is called without specifying the type of derivative computation
    // then use information from the wave function to determine it
    if (deriv_calc_type == DerivCalcType::Default) {
        if (!ref_wfn) {
            // If wavefunction doesn't have a reference wavefunction
            // itself, we assume that we're dealing with SCF.
            deriv_calc_type = DerivCalcType::SCF;
        } else {
            if (wfn_->density_fitted()) {
                deriv_calc_type = DerivCalcType::SCFandDF;
            } else {
                deriv_calc_type = DerivCalcType::Correlated;
            }
        }
    }

    if (deriv_calc_type == DerivCalcType::SCF) {
        // Derivatives code for SCF computations
        if (!Da || !Db) throw PSIEXCEPTION("Deriv::compute: Unable to access OPDM.");
        if (!X) throw PSIEXCEPTION("Deriv::compute: Unable to access Lagrangian.");

        if (wfn_->same_a_b_dens()) {  // RHF
            // We need to account for spin integration
            X->scale(2.0);
            ScfRestrictedFunctor functor(TPDMcont_vector, Da);
            so_eri.compute_integrals_deriv1(functor);
            functor.finalize();
        } else {  // ROHF and UHF
            ScfUnrestrictedFunctor functor(TPDMcont_vector, Da, Db);
            so_eri.compute_integrals_deriv1(functor);
            functor.finalize();
        }
        for (size_t cd = 0; cd < cdsalcs_.ncd(); ++cd) TPDMcont[cd] = TPDMcont_vector->get(cd);
    }

    if (deriv_calc_type == DerivCalcType::SCFandDF) {
        // Derivatives code for correlated calculations using density fitting
        // The Lagrangian and density matrices are held as member variables, but these contain only
        // the correlated part. The reference contributions are harvested from the
        // reference_wavefunction member of wfn_.
        // If density fitting was used, we don't want to compute two electron contributions here

        X_ref_cont_vector = std::make_shared<Vector>(ncd);
        TPDM_ref_cont_vector = std::make_shared<Vector>(ncd);
        X_ref_cont = X_ref_cont_vector->pointer();
        TPDM_ref_cont = TPDM_ref_cont_vector->pointer();
        x_ref_contr_ = factory_->create_shared_matrix("Reference Lagrangian contribution to gradient", natom_, 3);
        tpdm_ref_contr_ = factory_->create_shared_matrix("Reference two-electron contribution to gradient", natom_, 3);

        // Here we need to extract the reference contributions
        SharedMatrix X_ref = ref_wfn->lagrangian();
        SharedMatrix Da_ref = ref_wfn->Da();

        for (size_t cd = 0; cd < cdsalcs_.ncd(); ++cd) {
            double temp = -X_ref->vector_dot(s_deriv[cd]);
            X_ref_cont[cd] = temp;
        }

        if (wfn_->same_a_b_orbs()) {
            // In the restricted case, the alpha D is really the total D.  Undefine the beta one, so
            // that the one-electron contribution, computed below, is correct.
            Db = factory_->create_shared_matrix("nullptr");
            ScfRestrictedFunctor scf_functor(TPDM_ref_cont_vector, Da_ref);
            ScfAndDfCorrelationRestrictedFunctor functor(Dcont_vector, scf_functor, Da, Da_ref);
            so_eri.compute_integrals_deriv1(functor);
            functor.finalize();
        } else
            throw PSIEXCEPTION("Unrestricted DF gradient not implemented yet.");
    }

    if (deriv_calc_type == DerivCalcType::Correlated) {
        // Derivatives code for correlated calculations for older codes (CI/CC)
        // This is the part of the code reached from CI/CC.  In this case, the total (alpha+beta) density
        // matrices are backtransformed to the SO basis and dumped to disk.  The one particle terms are
        // just combined into the alpha density (with the beta OPDM set to zero, so that the one-particle
        // terms below are computed correctly.  The two-particle terms are computed the same in both cases
        // as all spin cases have been collapsed into the a single SO TPDM.

        if (!deriv_density_backtransformed_) {
            // Dial up an integral transformation object to backtransform the OPDM, TPDM and Lagrangian
            std::vector<std::shared_ptr<MOSpace> > spaces;
            spaces.push_back(MOSpace::all);
            std::shared_ptr<IntegralTransform> ints_transform =
                std::shared_ptr<IntegralTransform>(new IntegralTransform(
                    wfn_, spaces,
                    wfn_->same_a_b_orbs() ? IntegralTransform::TransformationType::Restricted
                                          : IntegralTransform::TransformationType::Unrestricted,  // Transformation type
                    IntegralTransform::OutputType::DPDOnly,                                       // Output buffer
                    IntegralTransform::MOOrdering::QTOrder,                                       // MO ordering
                    IntegralTransform::FrozenOrbitals::None));                                    // Frozen orbitals?
            dpd_set_default(ints_transform->get_dpd_id());

            // Some codes already presort the tpdm, do not follow this as an example
            if (tpdm_presorted_) ints_transform->set_tpdm_already_presorted(true);

            ints_transform->backtransform_density();

            Da = factory_->create_shared_matrix("SO-basis OPDM");
            Db = factory_->create_shared_matrix("nullptr");
            Da->load(_default_psio_lib_, PSIF_AO_OPDM);
            X = factory_->create_shared_matrix("SO-basis Lagrangian");
            X->load(_default_psio_lib_, PSIF_AO_OPDM);
            // The CC lagrangian is defined with a different prefactor to SCF / MP2, so we account for it here
            X->scale(0.5);
        }

        _default_psio_lib_->open(PSIF_AO_TPDM, PSIO_OPEN_OLD);
        CorrelatedFunctor functor(TPDMcont_vector);
        so_eri.compute_integrals_deriv1(functor);
        functor.finalize();
        _default_psio_lib_->close(PSIF_AO_TPDM, 1);

        for (size_t cd = 0; cd < cdsalcs_.ncd(); ++cd) TPDMcont[cd] = TPDMcont_vector->get(cd);
    }

    outfile->Printf("\n");

    // Now, compute the one electron terms
    auto Dtot_AO = std::make_shared<Matrix>("AO basis total D", wfn_->nso(), wfn_->nso());
    if (wfn_->density_fitted()) {
        Dtot_AO->add(wfn_->Da_subset("AO"));
        Dtot_AO->add(wfn_->Db_subset("AO"));
    }
    auto Dtot = Da->clone();
    Dtot->add(Db);
    Dtot_AO->remove_symmetry(Dtot, wfn_->aotoso()->transpose());
    auto opdm_contr = mints->core_hamiltonian_grad(Dtot_AO);
    for (size_t cd = 0; cd < cdsalcs_.ncd(); ++cd) {
        Xcont[cd] = -X->vector_dot(s_deriv[cd]);
    }

    // Transform the SALCs back to cartesian space
    auto st = cdsalcs_.matrix();
    auto B = st->pointer(0);
    auto *cart = new double[3 * natom_];

    if (TPDM_ref_cont) {
        // B^t g_q^t = g_x^t -> g_q B = g_x
        C_DGEMM('n', 'n', 1, 3 * natom_, cdsalcs_.ncd(), 1.0, TPDM_ref_cont, cdsalcs_.ncd(), B[0], 3 * natom_, 0.0,
                cart, 3 * natom_);

        for (int a = 0; a < natom_; ++a)
            for (int xyz = 0; xyz < 3; ++xyz) tpdm_ref_contr_->set(a, xyz, cart[3 * a + xyz]);
    } else {
        // B^t g_q^t = g_x^t -> g_q B = g_x
        C_DGEMM('n', 'n', 1, 3 * natom_, cdsalcs_.ncd(), 1.0, TPDMcont, cdsalcs_.ncd(), B[0], 3 * natom_, 0.0, cart,
                3 * natom_);

        for (int a = 0; a < natom_; ++a)
            for (int xyz = 0; xyz < 3; ++xyz) tpdm_contr_->set(a, xyz, cart[3 * a + xyz]);
    }

    // B^t g_q^t = g_x^t -> g_q B = g_x
    C_DGEMM('n', 'n', 1, 3 * natom_, cdsalcs_.ncd(), 1.0, Xcont, cdsalcs_.ncd(), B[0], 3 * natom_, 0.0, cart,
            3 * natom_);
    
    auto x_contr = factory_->create_shared_matrix("Lagrangian contribution to gradient", natom_, 3);
    for (int a = 0; a < natom_; ++a)
        for (int xyz = 0; xyz < 3; ++xyz) x_contr->set(a, xyz, cart[3 * a + xyz]);

    if (X_ref_cont) {
        // B^t g_q^t = g_x^t -> g_q B = g_x
        C_DGEMM('n', 'n', 1, 3 * natom_, cdsalcs_.ncd(), 1.0, X_ref_cont, cdsalcs_.ncd(), B[0], 3 * natom_, 0.0, cart,
                3 * natom_);

        for (int a = 0; a < natom_; ++a)
            for (int xyz = 0; xyz < 3; ++xyz) x_ref_contr_->set(a, xyz, cart[3 * a + xyz]);
    }

    // Obtain nuclear repulsion contribution from the wavefunction
    auto enuc = std::make_shared<Matrix>(molecule_->nuclear_repulsion_energy_deriv1(wfn_->get_dipole_field_strength()));

    // Print things out, after making sure that each component is properly symmetrized
    enuc->symmetrize_gradient(molecule_);
    opdm_contr->symmetrize_gradient(molecule_);
    x_contr->symmetrize_gradient(molecule_);
    tpdm_contr_->symmetrize_gradient(molecule_);

    if (wfn_->options().get_int("PRINT") > 1) {
        enuc->print_atom_vector();
        opdm_contr->print_atom_vector();
        x_contr->print_atom_vector();
        tpdm_contr_->print_atom_vector();
    }

    if (x_ref_contr_) {
        x_ref_contr_->symmetrize_gradient(molecule_);
        if (wfn_->options().get_int("PRINT") > 1) x_ref_contr_->print_atom_vector();
    }
    if (tpdm_ref_contr_) {
        tpdm_ref_contr_->symmetrize_gradient(molecule_);
        if (wfn_->options().get_int("PRINT") > 1) tpdm_ref_contr_->print_atom_vector();
    }

    // Add everything up into a temp.
    auto corr = std::make_shared<Matrix>("Correlation contribution to gradient", molecule_->natom(), 3);
    gradient_->add(enuc);
    corr->add(opdm_contr);
    corr->add(x_contr);
    corr->add(tpdm_contr_);
    if (reference_separate && !ignore_reference_) {
        gradient_->add(x_ref_contr_);
        gradient_->add(tpdm_ref_contr_);
        SharedMatrix scf_gradient(gradient_->clone());
        scf_gradient->set_name("Reference Gradient");
        scf_gradient->print_atom_vector();
        wfn_->reference_wavefunction()->set_gradient(scf_gradient);
        corr->print_atom_vector();
    }
    gradient_->add(corr);

    // Print the atom vector
    gradient_->print_atom_vector();

    // Save the gradient to the wavefunction so that optking can optimize with it
    wfn_->set_gradient(gradient_);

    return gradient_;
}

}  // namespace psi

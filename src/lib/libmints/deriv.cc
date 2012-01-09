/*
 * deriv.cc
 *
 *  Created on: Feb 24, 2009
 *      Author: jturney
 */

#include <stdio.h>
#include <stdlib.h>
#include <string>

#include <boost/foreach.hpp>

#include <libtrans/integraltransform.h>
#include <libdpd/dpd.h>
#include <boost/shared_ptr.hpp>

#include "mints.h"
#include "deriv.h"

using namespace std;

namespace psi {

size_t counter;

class CorrelatedFunctor
{
    /// The buffer to hold the TPDM
    double *tpdm_buffer_;
    /// Pointer to the current TPDM element
    double *tpdm_ptr_;
    /// How large the buffer is, for each shell pair
    size_t *buffer_sizes_;
    /// The PSIO object to use for disk I/O
    boost::shared_ptr<PSIO> psio_;
public:
    int nthread;
    std::vector<SharedVector> result;

    CorrelatedFunctor() {
        throw PSIEXCEPTION("CorrelatedRestrictedFunctor(): Default constructor called. This shouldn't happen.");
    }
    CorrelatedFunctor(SharedVector results) : psio_(_default_psio_lib_)
    {
        nthread = Communicator::world->nthread();
        result.push_back(results);
        for (int i=1; i<nthread; ++i)
            result.push_back(SharedVector(result[0]->clone()));
        size_t num_pairs = 0;
        psio_->read_entry(PSIF_AO_TPDM, "Num. Pairs", (char*)&num_pairs, sizeof(size_t));
        buffer_sizes_ = new size_t[num_pairs];
        psio_->read_entry(PSIF_AO_TPDM, "TPDM Buffer Sizes", (char*)buffer_sizes_, num_pairs*sizeof(size_t));
        size_t max_size = 0;
        for(int i = 0; i < num_pairs; ++i)
            max_size = max_size > buffer_sizes_[i] ? max_size : buffer_sizes_[i];
        tpdm_buffer_ = new double[max_size];
        tpdm_ptr_ = tpdm_buffer_;
    }

    void finalize() {
        // Do summation over threads
        for (int i=1; i<nthread; ++i) {
            result[0]->add(result[i]);
        }
        // Do MPI global summation
        result[0]->sum();
        delete [] tpdm_buffer_;
        delete [] buffer_sizes_;
    }

    void load_tpdm(size_t id){
        // TODO, make this work with threads (each thread needs its own buffer)
        char *toc = new char[40];
        sprintf(toc, "SO_TPDM_FOR_PAIR_%zd", id);
        size_t buffer_size = buffer_sizes_[id];
        psio_->read_entry(PSIF_AO_TPDM, toc, (char*)tpdm_buffer_, buffer_size*sizeof(double));
        delete [] toc;
        tpdm_ptr_ = tpdm_buffer_;
    }

    void next_tpdm_element(){
        ++tpdm_ptr_;
    }

    void operator()(int salc, int pabs, int qabs, int rabs, int sabs,
                    int pirrep, int pso,
                    int qirrep, int qso,
                    int rirrep, int rso,
                    int sirrep, int sso,
                    double value)
    {
        int thread = Communicator::world->thread_id(pthread_self());

        double prefactor = 8.0;
        if (pabs == qabs)
            prefactor *= 0.5;
        if (rabs == sabs)
            prefactor *= 0.5;
        if (pabs == rabs && qabs == sabs)
            prefactor *= 0.5;
        result[thread]->add(salc, prefactor * (*tpdm_ptr_) * value);
    }
};

class ScfRestrictedFunctor
{
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

    ScfRestrictedFunctor(SharedVector results, boost::shared_ptr<Matrix> D)
        : D_(D)
    {
        counter=0;
        nthread = Communicator::world->nthread();
        result.push_back(results);

        for (int i=1; i<nthread; ++i)
            result.push_back(SharedVector(result[0]->clone()));
    }
    ~ScfRestrictedFunctor() {
    }

    void finalize() {
        // Do summation over threads
        for (int i=1; i<nthread; ++i) {
            result[0]->add(result[i]);
        }
        // Do MPI global summation
        result[0]->sum();
    }

    void load_tpdm(size_t id) {}
    void next_tpdm_element(){}

    void operator()(int salc, int pabs, int qabs, int rabs, int sabs,
                    int pirrep, int pso,
                    int qirrep, int qso,
                    int rirrep, int rso,
                    int sirrep, int sso,
                    double value)
    {
        int thread = Communicator::world->thread_id(pthread_self());

        // Previously, we applied a factor of 4 after the fact...apply it from the beginning now.
        double prefactor = 4.0;

        if (pabs == qabs)
            prefactor *= 0.5;
        if (rabs == sabs)
            prefactor *= 0.5;
        if (pabs == rabs && qabs == sabs)
            prefactor *= 0.5;

        double four_index_D = 0.0;

        if (pirrep == qirrep && rirrep == sirrep)
            four_index_D = 4.0 * D_->get(pirrep, pso, qso) * D_->get(rirrep, rso, sso);
        if (pirrep == rirrep && qirrep == sirrep)
            four_index_D -= D_->get(pirrep, pso, rso) * D_->get(qirrep, qso, sso);
        if (pirrep == sirrep && qirrep == rirrep)
            four_index_D -= D_->get(pirrep, pso, sso) * D_->get(qirrep, qso, rso);

        four_index_D *= prefactor;
//        fprintf(outfile, "Salc %d (%d %d | %d %d) = %16.10f\n", salc, pabs, qabs, rabs, sabs, value);

        result[thread]->add(salc, four_index_D * value);
        counter++;
    }
};

class ScfAndDfCorrelationRestrictedFunctor
{
    SharedMatrix D_ref_;
    SharedMatrix D_;
    ScfRestrictedFunctor scf_functor_;
    std::vector<SharedVector> result_vec_;
    SharedVector results_;

public:
    int nthread;

//    ScfAndDfCorrelationRestrictedFunctor() {
//        throw PSIEXCEPTION("SCFAndDFCorrelationRestrictedFunctor(): Default constructor called. This shouldn't happen.");
//    }

    ScfAndDfCorrelationRestrictedFunctor(SharedVector results,
                                         ScfRestrictedFunctor& scf_functor,
                                         boost::shared_ptr<Matrix> D,
                                         boost::shared_ptr<Matrix> D_ref)
        : D_ref_(D_ref), D_(D), scf_functor_(scf_functor), results_(results)
    {
        counter=0;
        nthread = Communicator::world->nthread();
        result_vec_.push_back(results);

        for (int i=1; i<nthread; ++i)
            result_vec_.push_back(SharedVector(results->clone()));
    }

    ScfAndDfCorrelationRestrictedFunctor() {
        throw PSIEXCEPTION("ScfAndDfCorrelationRestrictedFunctor(): Default constructor called. This shouldn't happen.");
    }

    ~ScfAndDfCorrelationRestrictedFunctor() {
    }

    void load_tpdm(size_t id) {}
    void next_tpdm_element(){}

    void finalize() {
        // Make sure the SCF code is done
        scf_functor_.finalize();
        // Do summation over threads
        for (int i=1; i<nthread; ++i) {
            result_vec_[0]->add(result_vec_[i]);
        }
        // Do MPI global summation
        result_vec_[0]->sum();
    }

    void operator()(int salc, int pabs, int qabs, int rabs, int sabs,
                    int pirrep, int pso,
                    int qirrep, int qso,
                    int rirrep, int rso,
                    int sirrep, int sso,
                    double value)
    {
        int thread = Communicator::world->thread_id(pthread_self());

        bool braket = pabs!=rabs || qabs!=sabs;
        bool bra    = pabs!=qabs;
        bool ket    = rabs!=sabs;

        double four_index_D = 0.0;

        double Coulomb1  = 0.0;
        double Coulomb2  = 0.0;
        double Exchange1 = 0.0;
        if (pirrep == qirrep && rirrep == sirrep){
            Coulomb1 = 2.0 * D_->get(pirrep, pso, qso) * D_ref_->get(rirrep, rso, sso);
            Coulomb2 = 2.0 * D_->get(rirrep, rso, sso) * D_ref_->get(pirrep, pso, qso);
        }
        if (pirrep == rirrep && qirrep == sirrep)
            Exchange1 = D_->get(pirrep, pso, rso) * D_ref_->get(qirrep, qso, sso);
        // (pq|rs)
        four_index_D = Coulomb1 - Exchange1;

        if(bra && ket && braket){
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
        }else if(bra && ket){
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
        }else if(bra){
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
        }else if(ket){
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
        }else if(braket){
            four_index_D += Coulomb2;
            // (rs|pq)
            four_index_D -= Exchange1;
        }

        result_vec_[thread]->add(salc, four_index_D * value);

        // Make sure the SCF contribution is computed.
        scf_functor_(salc, pabs, qabs, rabs, sabs,  pirrep, pso, qirrep, qso,
                     rirrep, rso, sirrep, sso, value);
        counter++;
    }
};


class ScfUnrestrictedFunctor
{
    SharedMatrix Da_;
    SharedMatrix Db_;

public:
    int nthread;
    std::vector<SharedVector> result;

    ScfUnrestrictedFunctor() { throw PSIEXCEPTION("ScfUnrestrictedFunctor(): Oh come on!!!"); }

    ScfUnrestrictedFunctor(SharedVector results, boost::shared_ptr<Matrix> Da, boost::shared_ptr<Matrix> Db)
        : Da_(Da),
          Db_(Db)
    {
        nthread = Communicator::world->nthread();
        result.push_back(results);
        for (int i=1; i<nthread; ++i) {
            result.push_back(SharedVector(result[0]->clone()));
        }
    }
    ~ScfUnrestrictedFunctor() {
    }

    void load_tpdm(size_t id) {}
    void next_tpdm_element(){}

    void finalize() {
        // Do summation over threads
        for (int i=1; i<nthread; ++i) {
            result[0]->add(result[i]);
        }
        // Do MPI global summation
        result[0]->sum();
    }

    void operator()(int salc, int pabs, int qabs, int rabs, int sabs,
                    int pirrep, int pso,
                    int qirrep, int qso,
                    int rirrep, int rso,
                    int sirrep, int sso,
                    double value)
    {
        int thread = Communicator::world->thread_id(pthread_self());
        double prefactor = 1.0;

        if (pabs == qabs)
            prefactor *= 0.5;
        if (rabs == sabs)
            prefactor *= 0.5;
        if (pabs == rabs && qabs == sabs)
            prefactor *= 0.5;

        double four_index_D = 0.0;

        if (pirrep == qirrep && rirrep == sirrep) {
            four_index_D = 4.0 * (Da_->get(pirrep, pso, qso) + Db_->get(pirrep, pso, qso)) *
                                 (Da_->get(rirrep, rso, sso) + Db_->get(rirrep, rso, sso));
        }
        if (pirrep == rirrep && qirrep == sirrep) {
            four_index_D -= 2.0 * ((Da_->get(pirrep, pso, rso) * Da_->get(qirrep, qso, sso))
                                 + (Db_->get(pirrep, pso, rso) * Db_->get(qirrep, qso, sso)));
        }
        if (pirrep == sirrep && qirrep == rirrep) {
            four_index_D -= 2.0 * ((Da_->get(pirrep, pso, sso) * Da_->get(rirrep, rso, qso))
                                 + (Db_->get(pirrep, pso, sso) * Db_->get(rirrep, rso, qso)));
        }
//        four_index_D *= prefactor;
        value *= prefactor;

        result[thread]->add(salc, four_index_D * value);
    }
};

Deriv::Deriv(const boost::shared_ptr<Wavefunction>& wave,
             char needed_irreps,
             bool project_out_translations,
             bool project_out_rotations)
    : wfn_(wave),
      cdsalcs_(wave->molecule(),
          wave->matrix_factory(),
          needed_irreps,
          project_out_translations,
          project_out_rotations)
{
    integral_ = wave->integral();
    basis_    = wave->basisset();
    sobasis_  = wave->sobasisset();
    factory_  = wave->matrix_factory();
    molecule_ = wave->molecule();
    natom_    = molecule_->natom();

    // Results go here.
    opdm_contr_ = factory_->create_shared_matrix("One-electron contribution to gradient", natom_, 3);
    x_contr_    = factory_->create_shared_matrix("Lagrangian contribution to gradient", natom_, 3);
    tpdm_contr_ = factory_->create_shared_matrix("Two-electron contribution to gradient", natom_, 3);
    gradient_   = factory_->create_shared_matrix("Total gradient", natom_, 3);

    cdsalcs_.print();
}

SharedMatrix Deriv::compute()
{
    molecule_->print_in_bohr();

    if (natom_ == 1) {
        // This is an atom...there is no gradient.
        fprintf(outfile, "    A single atom has no gradient.\n");
        // Save the gradient to the wavefunction so that optking can optimize with it
        wfn_->set_gradient(gradient_);
        return gradient_;
    }

    // Initialize an ERI object requesting derivatives.
    std::vector<boost::shared_ptr<TwoBodyAOInt> > ao_eri;
    for (int i=0; i<Communicator::world->nthread(); ++i)
        ao_eri.push_back(boost::shared_ptr<TwoBodyAOInt>(integral_->eri(1)));
    TwoBodySOInt so_eri(ao_eri, integral_, cdsalcs_);

    // A certain optimization can be used if we know we only need totally symmetric
    // derivatives.
    so_eri.set_only_totally_symmetric(true);

    // Compute one-electron derivatives.
    vector<SharedMatrix> s_deriv = cdsalcs_.create_matrices("S'");
    vector<SharedMatrix> h_deriv = cdsalcs_.create_matrices("H'");

    boost::shared_ptr<OneBodySOInt> s_int(integral_->so_overlap(1));
    boost::shared_ptr<OneBodySOInt> t_int(integral_->so_kinetic(1));
    boost::shared_ptr<OneBodySOInt> v_int(integral_->so_potential(1));

    s_int->compute_deriv1(s_deriv, cdsalcs_);
    t_int->compute_deriv1(h_deriv, cdsalcs_);
    v_int->compute_deriv1(h_deriv, cdsalcs_);

    int ncd = cdsalcs_.ncd();
    SharedVector TPDMcont_vector(new Vector(1, &ncd));
    SharedVector Xcont_vector(new Vector(1, &ncd));
    SharedVector Dcont_vector(new Vector(1, &ncd));
    SharedVector TPDM_ref_cont_vector;
    SharedVector X_ref_cont_vector;
    SharedVector D_ref_cont_vector;
    double *Dcont           = Dcont_vector->pointer();
    double *Xcont           = Xcont_vector->pointer();
    double *TPDMcont        = TPDMcont_vector->pointer();
    double *TPDM_ref_cont   = 0;
    double *X_ref_cont      = 0;
    double *D_ref_cont      = 0;

    if(!wfn_)
        throw("In Deriv: The wavefunction passed in is empty!");

    // Try and grab the OPDM and lagrangian from the wavefunction
    SharedMatrix Da = wfn_->Da();
    SharedMatrix Db = wfn_->Db();
    SharedMatrix X = wfn_->X();

    // The current wavefunction's reference wavefunction, NULL for SCF/DFT
    boost::shared_ptr<Wavefunction> ref_wfn = wfn_->reference_wavefunction();
    // Whether the SCF contribution is separate from the correlated terms
    bool reference_separate = (Da || Db || X) && ref_wfn;

    if(!ref_wfn){
        // If wavefunction doesn't have a reference wavefunction
        // itself, we assume that we're dealing with SCF.
        if (!Da || !Db)
            throw PSIEXCEPTION("Deriv::compute: Unable to access OPDM.");
        if (!X)
            throw PSIEXCEPTION("Deriv::compute: Unable to access Lagrangian.");

        if(wfn_->same_a_b_orbs()){
            // We need to account for spin integration
            X->scale(2.0);
            ScfRestrictedFunctor functor(TPDMcont_vector, Da);
            so_eri.compute_integrals_deriv1(functor);
            functor.finalize();
        }else{
            ScfUnrestrictedFunctor functor(TPDMcont_vector, Da, Db);
            so_eri.compute_integrals_deriv1(functor);
            functor.finalize();
        }
        for (int cd=0; cd < cdsalcs_.ncd(); ++cd) {
            TPDMcont[cd] = TPDMcont_vector->get(cd);
//            fprintf(outfile, "    SALC #%d TPDM contribution:         %+lf\n", cd, TPDMcont[cd]);
        }
        fflush(outfile);
    } else {
        /* For correlated calculations, we have two different types.  The older CI/CC codes dump the
           Lagrangian to disk and density matrices to disk, and these both include the reference
           contributions.  The newer codes hold these quantities as member variables, but these contain only
           the correlated part.  The reference contributions must be harvested from the reference_wavefunction
           member.  If density fitting was used, we don't want to compute two electron contributions here*/
        if(wfn_->density_fitted()){
            X_ref_cont_vector    = SharedVector(new Vector(1, &ncd));
            D_ref_cont_vector    = SharedVector(new Vector(1, &ncd));
            TPDM_ref_cont_vector = SharedVector(new Vector(1, &ncd));
            X_ref_cont           = X_ref_cont_vector->pointer();
            D_ref_cont           = D_ref_cont_vector->pointer();
            TPDM_ref_cont        = TPDM_ref_cont_vector->pointer();
            x_ref_contr_         = factory_->create_shared_matrix("Reference Lagrangian contribution to gradient", natom_, 3);
            opdm_ref_contr_      = factory_->create_shared_matrix("Reference one-electron contribution to gradient", natom_, 3);
            tpdm_ref_contr_      = factory_->create_shared_matrix("Reference two-electron contribution to gradient", natom_, 3);

            // Here we need to extract the reference contributions
            SharedMatrix X_ref  = ref_wfn->Lagrangian();
            SharedMatrix Da_ref = ref_wfn->Da();
            SharedMatrix Db_ref = ref_wfn->Db();

            for (int cd=0; cd < cdsalcs_.ncd(); ++cd) {
                double temp = Da_ref->vector_dot(h_deriv[cd]);
                temp += Db_ref->vector_dot(h_deriv[cd]);
                D_ref_cont[cd] = temp;
//                fprintf(outfile, "    SALC #%d Reference One-electron contribution: %+lf\n", cd, temp);
            }
            fprintf(outfile, "\n");

            for (int cd=0; cd < cdsalcs_.ncd(); ++cd) {
                double temp = -X_ref->vector_dot(s_deriv[cd]);
                X_ref_cont[cd] = temp;
//                fprintf(outfile, "    SALC #%d Reference Lagrangian contribution:   %+lf\n", cd, temp);
            }
            fprintf(outfile, "\n");

            if(wfn_->same_a_b_orbs()){
                // In the restricted case, the alpha D is really the total D.  Undefine the beta one, so
                // that the one-electron contribution, computed below, is correct.
                Db = factory_->create_shared_matrix("NULL");
                ScfRestrictedFunctor scf_functor(TPDM_ref_cont_vector, Da_ref);
                ScfAndDfCorrelationRestrictedFunctor functor(Dcont_vector, scf_functor, Da, Da_ref);
                so_eri.compute_integrals_deriv1(functor);
                functor.finalize();
                tpdm_contr_ = wfn_->tpdm_gradient_contribution();
            }else{
                throw PSIEXCEPTION("Unrestricted DF gradient not implemented yet.");
            }
        }else{
            /* This is the part of the code reached from CI/CC.  In this case, the total (alpha+beta) density
               matrices are backtransformed to the SO basis and dumped to disk.  The one particle terms are
               just combined into the alpha density (with the beta OPDM set to zero, so that the one-particle
               terms below are computed correctly.  The two-particle terms are computed the same in both cases
               as all spin cases have been collapsed into the a single SO TPDM. */

            // Dial up an integral tranformation object to backtransform the OPDM, TPDM and Lagrangian
            vector<boost::shared_ptr<MOSpace> > spaces;
            spaces.push_back(MOSpace::all);
            boost::shared_ptr<IntegralTransform> ints_transform = boost::shared_ptr<IntegralTransform>(
                        new IntegralTransform(wfn_,
                                              spaces,
                                              wfn_->same_a_b_orbs() ? IntegralTransform::Restricted : IntegralTransform::Unrestricted, // Transformation type
                                              IntegralTransform::DPDOnly,    // Output buffer
                                              IntegralTransform::QTOrder,    // MO ordering
                                              IntegralTransform::None));     // Frozen orbitals?
            dpd_set_default(ints_transform->get_dpd_id());
            ints_transform->backtransform_density();

            Da = factory_->create_shared_matrix("SO-basis OPDM");
            Db = factory_->create_shared_matrix("NULL");
            Da->load(_default_psio_lib_, PSIF_AO_OPDM);
            X = factory_->create_shared_matrix("SO-basis Lagrangian");
            X->load(_default_psio_lib_, PSIF_AO_OPDM);
            // The CC lagrangian is defined with a different prefactor to SCF / MP2, so we account for it here
            X->scale(0.5);

            _default_psio_lib_->open(PSIF_AO_TPDM, PSIO_OPEN_OLD);
            CorrelatedFunctor functor(TPDMcont_vector);
            so_eri.compute_integrals_deriv1(functor);
            functor.finalize();
            _default_psio_lib_->close(PSIF_AO_TPDM, 1);

            for (int cd=0; cd < cdsalcs_.ncd(); ++cd) {
                TPDMcont[cd] = TPDMcont_vector->get(cd);
//                fprintf(outfile, "    SALC #%d TPDM contribution:         %+lf\n", cd, TPDMcont[cd]);
            }
            fflush(outfile);
        }

        fprintf(outfile, "\n");
    }

    // Now, compute the one electron terms
    for (int cd=0; cd < cdsalcs_.ncd(); ++cd) {
        double temp = Dcont[cd]; // In the df case, the HxP2 terms are already in here
        temp += Da->vector_dot(h_deriv[cd]);
        temp += Db->vector_dot(h_deriv[cd]);
        Dcont[cd] = temp;
//        fprintf(outfile, "    SALC #%d One-electron contribution: %+lf\n", cd, temp);
    }
    fprintf(outfile, "\n");

    for (int cd=0; cd < cdsalcs_.ncd(); ++cd) {
        double temp = X->vector_dot(s_deriv[cd]);
        Xcont[cd] = -temp;
//        fprintf(outfile, "    SALC #%d Lagrangian contribution:   %+lf\n", cd, temp);
    }
    fprintf(outfile, "\n");

    // Transform the SALCs back to cartesian space
    SharedMatrix st = cdsalcs_.matrix();
    double **B = st->pointer(0);
    double *cart = new double[3*natom_];

    // B^t g_q^t = g_x^t -> g_q B = g_x
    C_DGEMM('n', 'n', 1, 3*natom_, cdsalcs_.ncd(),
            1.0, Dcont, cdsalcs_.ncd(), B[0],
            3*natom_, 0.0, cart, 3*natom_);

    for (int a=0; a<natom_; ++a)
        for (int xyz=0; xyz<3; ++xyz)
            opdm_contr_->set(a, xyz, cart[3*a+xyz]);

    if(D_ref_cont){
        // B^t g_q^t = g_x^t -> g_q B = g_x
        C_DGEMM('n', 'n', 1, 3*natom_, cdsalcs_.ncd(),
                1.0, D_ref_cont, cdsalcs_.ncd(), B[0],
                3*natom_, 0.0, cart, 3*natom_);

        for (int a=0; a<natom_; ++a)
            for (int xyz=0; xyz<3; ++xyz)
                opdm_ref_contr_->set(a, xyz, cart[3*a+xyz]);
    }

    if(TPDM_ref_cont){
        // B^t g_q^t = g_x^t -> g_q B = g_x
        C_DGEMM('n', 'n', 1, 3*natom_, cdsalcs_.ncd(),
                1.0, TPDM_ref_cont, cdsalcs_.ncd(), B[0],
                3*natom_, 0.0, cart, 3*natom_);

        for (int a=0; a<natom_; ++a)
            for (int xyz=0; xyz<3; ++xyz)
                tpdm_ref_contr_->set(a, xyz, cart[3*a+xyz]);
    }else{
        // B^t g_q^t = g_x^t -> g_q B = g_x
        C_DGEMM('n', 'n', 1, 3*natom_, cdsalcs_.ncd(),
                1.0, TPDMcont, cdsalcs_.ncd(), B[0],
                3*natom_, 0.0, cart, 3*natom_);

        for (int a=0; a<natom_; ++a)
            for (int xyz=0; xyz<3; ++xyz)
                tpdm_contr_->set(a, xyz, cart[3*a+xyz]);
    }

    // B^t g_q^t = g_x^t -> g_q B = g_x
    C_DGEMM('n', 'n', 1, 3*natom_, cdsalcs_.ncd(),
            1.0, Xcont, cdsalcs_.ncd(), B[0],
            3*natom_, 0.0, cart, 3*natom_);

    for (int a=0; a<natom_; ++a)
        for (int xyz=0; xyz<3; ++xyz)
            x_contr_->set(a, xyz, cart[3*a+xyz]);

    if(X_ref_cont){
        // B^t g_q^t = g_x^t -> g_q B = g_x
        C_DGEMM('n', 'n', 1, 3*natom_, cdsalcs_.ncd(),
                1.0, X_ref_cont, cdsalcs_.ncd(), B[0],
                3*natom_, 0.0, cart, 3*natom_);

        for (int a=0; a<natom_; ++a)
            for (int xyz=0; xyz<3; ++xyz)
                x_ref_contr_->set(a, xyz, cart[3*a+xyz]);
    }

    // Obtain nuclear repulsion contribution from the wavefunction
    SharedMatrix enuc(new Matrix(molecule_->nuclear_repulsion_energy_deriv1()));

    // Print things out, after making sure that each component is properly symmetrized
    symmetrize_gradient(enuc)->print_atom_vector();
    symmetrize_gradient(opdm_contr_)->print_atom_vector();
    symmetrize_gradient(x_contr_)->print_atom_vector();
    symmetrize_gradient(tpdm_contr_)->print_atom_vector();
    if(x_ref_contr_)
        symmetrize_gradient(x_ref_contr_)->print_atom_vector();
    if(opdm_ref_contr_)
        symmetrize_gradient(opdm_ref_contr_)->print_atom_vector();
    if(tpdm_ref_contr_)
        symmetrize_gradient(tpdm_ref_contr_)->print_atom_vector();

    // Add everything up into a temp.
    SharedMatrix corr(new Matrix("Correlation contribution to gradient", molecule_->natom(), 3));
    gradient_->add(enuc);
    corr->add(opdm_contr_);
    corr->add(x_contr_);
    corr->add(tpdm_contr_);
    if(reference_separate){
        gradient_->add(x_ref_contr_);
        gradient_->add(opdm_ref_contr_);
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

SharedMatrix
Deriv::symmetrize_gradient(SharedMatrix grad)
{
    // Make a temporary storage object
    SharedMatrix temp(grad->clone());
    temp->zero();

    CharacterTable ct = molecule_->point_group()->char_table();

    // Obtain atom mapping of atom * symm op to atom
    int **atom_map = compute_atom_map(molecule_);

    // Symmetrize the gradients to remove any noise
    for (int atom=0; atom<molecule_->natom(); ++atom) {
        for (int g=0; g<ct.order(); ++g) {

            int Gatom = atom_map[atom][g];

            SymmetryOperation so = ct.symm_operation(g);

            temp->add(atom, 0, so(0, 0) * grad->get(Gatom, 0) / ct.order());
            temp->add(atom, 1, so(1, 1) * grad->get(Gatom, 1) / ct.order());
            temp->add(atom, 2, so(2, 2) * grad->get(Gatom, 2) / ct.order());
        }
    }
    // Delete the atom map.
    delete_atom_map(atom_map, molecule_);

    grad->copy(temp);
    return grad;
}


}

#ifdef HAVE_MADNESS
namespace madness {
namespace archive {

template <class Archive>
struct ArchiveStoreImpl< Archive, psi::ScfRestrictedFunctor> {
    static void store(const Archive &ar, const psi::ScfRestrictedFunctor &t) {
    }
};

template <class Archive>
struct ArchiveStoreImpl< Archive, psi::ScfUnrestrictedFunctor> {
    static void store(const Archive &ar, const psi::ScfUnrestrictedFunctor &t) {
    }
};

template <class Archive>
struct ArchiveStoreImpl< Archive, psi::CorrelatedFunctor> {
    static void store(const Archive &ar, const psi::CorrelatedFunctor &t) {
    }
};

template <class Archive>
struct ArchiveStoreImpl< Archive, psi::ScfAndDfCorrelationRestrictedFunctor> {
    static void store(const Archive &ar, const psi::ScfAndDfCorrelationRestrictedFunctor &t) {
    }
};

template <class Archive>
struct ArchiveLoadImpl< Archive, psi::ScfRestrictedFunctor> {
    static void load(const Archive &ar, const psi::ScfRestrictedFunctor &t) {
    }
};

template <class Archive>
struct ArchiveLoadImpl< Archive, psi::ScfUnrestrictedFunctor> {
    static void load(const Archive &ar, const psi::ScfUnrestrictedFunctor &t) {
    }
};

template <class Archive>
struct ArchiveLoadImpl< Archive, psi::CorrelatedFunctor> {
    static void load(const Archive &ar, const psi::CorrelatedFunctor &t) {
    }
};

template <class Archive>
struct ArchiveLoadImpl< Archive, psi::ScfAndDfCorrelationRestrictedFunctor> {
    static void load(const Archive &ar, const psi::ScfAndDfCorrelationRestrictedFunctor &t) {
    }
};

}
}
#endif // HAVE_MADNESS

#ifndef INTEGRALFUNCTORS_H
#define INTEGRALFUNCTORS_H

#include "hf.h"
#include <libiwl/iwl.hpp>
#include <libmints/mints.h>
#include "pseudospectral.h"
#include "pkintegrals.h"
#include "df.h"

namespace boost {
template<class T> class shared_ptr;
}

// We're assuming that for (PQ|RS) P>=Q, R>=S and PQ>=RS

/*
 * All functors MUST implement the following, or the code won't compile.
 * bool k_required() const - just a boolean stating whether this functor requires K
 * void initialize() this will always be called at the beginning and should do things
 *                   like clear memory, if needed.  It can just do nothing.
 * void finalize()   this will always be called at the end and should do things like
 *                   free memory and copy data, if needed.  It can just do nothing.
 * void operator()(int pabs, int qabs, int rabs, int sabs,
 *                 int psym, int prel, int qsym, int qrel,
 *                 int rsym, int rrel, int ssym, int srel, double value)
 *                   this is called by the out-of-core and direct codes, and will take
 *                   the relative (within irrep) and absolute values of p, q, r and s, along
 *                   with their symmetries, to generate the J and K matrices.
 * void operator()(boost::shared_ptr<DFHF> dfhf, boost::shared_ptr<PseudospectralHF> psHF) which will set up a PSHF object
 * void operator()(boost::shared_ptr<DSHF> dfHF) which will set up a DFHF object
 * void operator()(boost::shared_ptr<PKIntegrals> pk_integrals) which is used in the PK algorithms.
 */


namespace psi{ namespace scf{

/**
  @brief This can be passed into the templated HF::process_tei() function
         and will compute Coulomb and exchange matrices for restricted orbitals
*/
class J_K_Functor
{
    /// The density matrix
    const boost::shared_ptr<Matrix> D_;
    /// The occupation matrix
    const boost::shared_ptr<Matrix> C_;
    /// The number of alpha (and beta) electrons
    const int* N_;
    std::vector<SharedMatrix> J_;
    std::vector<SharedMatrix> K_;
    /// The communicator
    std::string comm;
    /// The scf type
    std::string scf_type;
    /// Number of threads
    int nthread_;

public:
    bool k_required() const {return true;}
    void initialize(std::string scf_type_){
        comm = Communicator::world->communicator();
        nthread_ = Communicator::world->nthread();
        scf_type = scf_type_;
        J_[0]->zero();
        K_[0]->zero();

        for (int i=1; i<nthread_; ++i) {
            J_.push_back(SharedMatrix(J_[0]->clone()));
            K_.push_back(SharedMatrix(K_[0]->clone()));
        }
    }
    void finalize(){

        // Sum up all the threads
        for (int i=1; i<nthread_; ++i) {
            J_[0]->add(J_[i]);
            K_[0]->add(K_[i]);
        }

        // Perform MPI global sum
        if (comm != "LOCAL" && scf_type == "DIRECT") {
            for (int i=0; i < J_[0]->nirrep(); i++)
                Communicator::world->sum(J_[0]->get_pointer(i), J_[0]->rowdim(i)*J_[0]->coldim(i));
            for (int i=0; i < K_[0]->nirrep(); i++)
                Communicator::world->sum(K_[0]->get_pointer(i), K_[0]->rowdim(i)*K_[0]->coldim(i));
        }

        J_[0]->copy_lower_to_upper();
        K_[0]->copy_lower_to_upper();
    }

    J_K_Functor() {
        throw PSIEXCEPTION("J_K_Functor(): Default constructor called. This shouldn't happen.");
    }

    // This was added for testing.
//    J_K_Functor(const J_K_Functor&) {
//        throw PSIEXCEPTION("J_K_Functor(): Copy constructor called.");
//    }
//
//    J_K_Functor& operator=(const J_K_Functor&) {
//        throw PSIEXCEPTION("J_K_Functor(): Assignment operator called. This shouldn't happen.");
//        return *this;
//    }

    J_K_Functor(boost::shared_ptr<Matrix> J, boost::shared_ptr<Matrix> K, const boost::shared_ptr<Matrix> D,
        const boost::shared_ptr<Matrix> C, const int* N)
        :  D_(D), C_(C), N_(N)
    {
        J_.push_back(J);
        K_.push_back(K);
    }

    void operator()(boost::shared_ptr<DFHF> dfhf, boost::shared_ptr<PseudospectralHF> pshf) {
        dfhf->set_restricted(true);
        dfhf->set_jk(false);
        dfhf->set_J(J_[0]);
        dfhf->set_Da(D_);
        pshf->set_restricted(true);
        pshf->set_Da(D_);
        pshf->set_Ka(K_[0]);
    }

    void operator()(boost::shared_ptr<PKIntegrals> pk_integrals) {
        pk_integrals->setup(J_[0], K_[0], D_, D_);
    }

    void operator()(boost::shared_ptr<DFHF> dfhf) {
        dfhf->set_restricted(true);
        dfhf->set_jk(true);
        dfhf->set_J(J_[0]);
        dfhf->set_Ka(K_[0]);
        dfhf->set_Da(D_);
        dfhf->set_Ca(C_);
        dfhf->set_Na(N_);
    }

    void operator()(int pabs, int qabs, int rabs, int sabs,
                int psym, int prel, int qsym, int qrel,
                int rsym, int rrel, int ssym, int srel, double value) {
        int thread = Communicator::world->thread_id(pthread_self());

        /* (pq|rs) */
        if(rsym == ssym){
            J_[thread]->add(rsym, rrel, srel, 2.0 * D_->get(psym, prel, qrel) * value);
        }
        if(qabs >= rabs){
            if(qsym == rsym){
                K_[thread]->add(qsym, qrel, rrel, D_->get(psym, prel, srel) * value);
            }
        }

        if(pabs!=qabs && rabs!=sabs && (pabs!=rabs || qabs!=sabs)){
            /* (pq|sr) */
            if(qabs >= sabs){
                if(qsym == ssym){
                    K_[thread]->add(qsym, qrel, srel, D_->get(psym, prel, rrel) * value);
                }
            }

            /* (qp|rs) */
            if(rsym == ssym){
                J_[thread]->add(rsym, rrel, srel, 2.0 * D_->get(qsym, qrel, prel) * value);
            }

            if(pabs >= rabs){
                if(psym == rsym){
                    K_[thread]->add(psym, prel, rrel, D_->get(qsym, qrel, srel) * value);
                }
            }

            /* (qp|sr) */
            if(pabs >= sabs){
                if(psym == ssym){
                    K_[thread]->add(psym, prel, srel, D_->get(qsym, qrel, rrel) * value);
                }
            }

            /* (rs|pq) */
            if(psym == qsym){
                J_[thread]->add(psym, prel, qrel, 2.0 * D_->get(rsym, rrel, srel) * value);
            }

            if(sabs >= pabs){
                if(ssym == psym){
                    K_[thread]->add(ssym, srel, prel, D_->get(rsym, rrel, qrel) * value);
                }
            }

            /* (sr|pq) */
            if(psym == qsym){
                J_[thread]->add(psym, prel, qrel, 2.0 * D_->get(ssym, srel, rrel) * value);
            }

            if(rabs >= pabs){
                if(rsym == psym){
                    K_[thread]->add(rsym, rrel, prel, D_->get(ssym, srel, qrel) * value);
                }
            }

            /* (rs|qp) */
            if(sabs >= qabs){
                if(ssym == qsym){
                    K_[thread]->add(ssym, srel, qrel, D_->get(rsym, rrel, prel) * value);
                }
            }

            /* (sr|qp) */
            if(rabs >= qabs){
                if(rsym == qsym){
                    K_[thread]->add(rsym, rrel, qrel, D_->get(ssym, srel, prel) * value);
                }
            }
        }else if(pabs!=qabs && rabs!=sabs && pabs==rabs && qabs==sabs){
            /* (pq|sr) */
            if(qabs >= sabs){
                if(qsym == ssym){
                    K_[thread]->add(qsym, qrel, srel, D_->get(psym, prel, rrel) * value);
                }
            }
            /* (qp|rs) */
            if(rsym == ssym){
                J_[thread]->add(rsym, rrel, srel, 2.0 * D_->get(qsym, qrel, prel) * value);
            }
            if(pabs >= rabs){
                if(psym == rsym){
                    K_[thread]->add(psym, prel, rrel, D_->get(qsym, qrel, srel) * value);
                }
            }

            /* (qp|sr) */
            if(pabs >= sabs){
                if(psym == ssym){
                    K_[thread]->add(psym, prel, srel, D_->get(qsym, qrel, rrel) * value);
                }
            }
        }else if(pabs!=qabs && rabs==sabs){
            /* (qp|rs) */
            if(rsym == ssym){
                J_[thread]->add(rsym, rrel, srel, 2.0 * D_->get(qsym, qrel, prel) * value);
            }

            if(pabs >= rabs){
                if(qsym == rsym){
                    K_[thread]->add(psym, prel, rrel, D_->get(qsym, qrel, srel) * value);
                }
            }

            /* (rs|pq) */
            if(psym == qsym){
                J_[thread]->add(psym, prel, qrel, 2.0 * D_->get(rsym, rrel, srel) * value);
            }

            if(sabs >= pabs){
                if(ssym == psym){
                    K_[thread]->add(ssym, srel, prel, D_->get(rsym, rrel, qrel) * value);
                }
            }

            /* (rs|qp) */
            if(sabs >= qabs){
                if(ssym == qsym){
                    K_[thread]->add(ssym, srel, qrel, D_->get(rsym, rrel, prel) * value);
                }
            }
        }else if(pabs==qabs && rabs!=sabs){
            /* (pq|sr) */
            if(qabs >= sabs){
                if(qsym == ssym){
                    K_[thread]->add(qsym, qrel, srel, D_->get(psym, prel, rrel) * value);
                }
            }

            /* (rs|pq) */
            if(psym == qsym){
                J_[thread]->add(psym, prel, qrel, 2.0 * D_->get(rsym, rrel, srel) * value);
            }

            if(sabs >= pabs){
                if(ssym == psym){
                    K_[thread]->add(ssym, srel, prel, D_->get(rsym, rrel, qrel) * value);
                }
            }

            /* (sr|pq) */
            if(psym == qsym){
                J_[thread]->add(psym, prel, qrel, 2.0 * D_->get(ssym, srel, rrel) * value);
            }

            if(rabs >= pabs){
                if(rsym == psym){
                    K_[thread]->add(rsym, rrel, prel, D_->get(ssym, srel, qrel) * value);
                }
            }
        }else if(pabs==qabs && rabs==sabs && (pabs!=rabs || qabs!=sabs)){
            /* (rs|pq) */
            if(psym == qsym){
                J_[thread]->add(psym, prel, qrel, 2.0 * D_->get(rsym, rrel, srel) * value);
            }

            if(sabs >= pabs){
                if(ssym == psym){
                    K_[thread]->add(ssym, srel, prel, D_->get(rsym, rrel, qrel) * value);
                }
            }
        }
    }
};


/**
  @brief This can be passed into the templated HF::process_tei() function
         and will compute only the Coulomb matrix for restricted orbitals
*/
class J_Functor
{
    /// The alpha density matrix
    const boost::shared_ptr<Matrix> Da_;
    /// The beta density matrix
    const boost::shared_ptr<Matrix> Db_;
    /// The Coulomb matrix
    std::vector<SharedMatrix> J_;
    /// Whether this is restricted or not
    bool restricted_;
    /// The communicator
    std::string comm;
    int nthread_;
    /// The scf type
    std::string scf_type;

public:
    bool k_required() const {return false;}

    void initialize(std::string scf_type_){
        comm = Communicator::world->communicator();
        nthread_ = Communicator::world->nthread();
        scf_type = scf_type_;
        J_[0]->zero();

        for (int i=1; i<nthread_; ++i)
            J_.push_back(SharedMatrix(J_[0]->clone()));
    }
    void finalize(){
        for (int i=1; i<nthread_; ++i)
            J_[0]->add(J_[i]);

        if (comm != "LOCAL" && scf_type == "DIRECT") {
            for (int i=0; i < J_[0]->nirrep(); i++)
                Communicator::world->sum(J_[0]->get_pointer(i), J_[0]->rowdim(i)*J_[0]->coldim(i));
        }

        J_[0]->copy_lower_to_upper();
    }

    // NEVER CALL THIS FUNCTION
    J_Functor() { throw PSIEXCEPTION("J_Functor(): This should never have been called, idiot."); }

    J_Functor(boost::shared_ptr<Matrix> J, const boost::shared_ptr<Matrix> Da)
        : Da_(Da), Db_(Da)
    {
        J_.push_back(J);
        restricted_ = true;
    }

    J_Functor(boost::shared_ptr<Matrix> J, const boost::shared_ptr<Matrix> Da, const boost::shared_ptr<Matrix> Db)
        : Da_(Da), Db_(Db)
    {
        J_.push_back(J);
        restricted_ = false;
    }

    void operator()(boost::shared_ptr<DFHF> dfhf, boost::shared_ptr<PseudospectralHF> pshf) {
        // Pseudospectral is exchange only, this does nothing
    }

    void operator()(boost::shared_ptr<DFHF> dfhf) {
        dfhf->set_restricted(restricted_);
        dfhf->set_jk(false);
        dfhf->set_J(J_[0]);
        dfhf->set_Da(Da_);
        dfhf->set_Db(Db_);
    }

    void operator()(boost::shared_ptr<PKIntegrals> pk_integrals) {
        pk_integrals->setup(J_[0], Da_, Db_);
    }

    void operator()(int pabs, int qabs, int rabs, int sabs,
                    int psym, int prel, int qsym, int qrel,
                    int rsym, int rrel, int ssym, int srel, double value) {
        int thread = Communicator::world->thread_id(pthread_self());
        /* (pq|rs) */
        if(rsym == ssym){
            J_[thread]->add(rsym, rrel, srel, (Da_->get(psym, prel, qrel) + Db_->get(psym, prel, qrel)) * value);
        }

        if(pabs!=qabs && rabs!=sabs && (pabs!=rabs || qabs!=sabs)){
            /* (pq|sr) */

            /* (qp|rs) */
            if(rsym == ssym){
                J_[thread]->add(rsym, rrel, srel, (Da_->get(qsym, qrel, prel) + Db_->get(qsym, qrel, prel)) * value);
            }

            /* (qp|sr) */

            /* (rs|pq) */
            if(psym == qsym){
                J_[thread]->add(psym, prel, qrel, (Da_->get(rsym, rrel, srel) + Db_->get(rsym, rrel, srel)) * value);
            }

            /* (sr|pq) */
            if(psym == qsym){
                J_[thread]->add(psym, prel, qrel, (Da_->get(ssym, srel, rrel) + Db_->get(ssym, srel, rrel)) * value);
            }

            /* (rs|qp) */

            /* (sr|qp) */
        }else if(pabs!=qabs && rabs!=sabs && pabs==rabs && qabs==sabs){
            /* (pq|sr) */

            /* (qp|rs) */
            if(rsym == ssym){
                J_[thread]->add(rsym, rrel, srel, (Da_->get(qsym, qrel, prel) + Db_->get(qsym, qrel, prel)) * value);
            }

            /* (qp|sr) */

        }else if(pabs!=qabs && rabs==sabs){
            /* (qp|rs) */
            if(rsym == ssym){
                J_[thread]->add(rsym, rrel, srel, (Da_->get(qsym, qrel, prel) + Db_->get(qsym, qrel, prel)) * value);
            }

            /* (rs|pq) */
            if(psym == qsym){
                J_[thread]->add(psym, prel, qrel, (Da_->get(rsym, rrel, srel) + Db_->get(rsym, rrel, srel)) * value);
            }
            /* (rs|qp) */

        }else if(pabs==qabs && rabs!=sabs){
            /* (pq|sr) */

            /* (rs|pq) */
            if(psym == qsym){
                J_[thread]->add(psym, prel, qrel, (Da_->get(rsym, rrel, srel) + Db_->get(rsym, rrel, srel)) * value);
            }

            /* (sr|pq) */
            if(psym == qsym){
                J_[thread]->add(psym, prel, qrel, (Da_->get(ssym, srel, rrel) + Db_->get(ssym, srel, rrel)) * value);
            }
        }else if(pabs==qabs && rabs==sabs && (pabs!=rabs || qabs!=sabs)){
            /* (rs|pq) */
            if(psym == qsym){
                J_[thread]->add(psym, prel, qrel, (Da_->get(rsym, rrel, srel) + Db_->get(rsym, rrel, srel)) * value);
            }
        }
    }
};


/**
  @brief This can be passed into the templated HF::process_tei() function
         and will compute alpha and beta Coulomb and exchange matrices
*/
class J_Ka_Kb_Functor
{
    /// The alpha density matrix
    const boost::shared_ptr<Matrix> Da_;
    /// The beta density matrix
    const boost::shared_ptr<Matrix> Db_;
    /// The alpha occupation matrix
    const boost::shared_ptr<Matrix> Ca_;
    /// The beta occupation matrix
    const boost::shared_ptr<Matrix> Cb_;
    /// The number of alpha electrons
    const int* Na_;
    /// The number of beta electrons
    const int* Nb_;
    /// The alpha Coulomb matrix
    std::vector<SharedMatrix> J_;
    /// The alpha exchange matrix
    std::vector<SharedMatrix> Ka_;
    /// The beta exchange matrix
    std::vector<SharedMatrix> Kb_;
    /// The communicator
    std::string comm;
    int nthread_;
    /// The scf type
    std::string scf_type;

public:
    bool k_required() const {return true;}

    void initialize(std::string scf_type_){
        comm = Communicator::world->communicator();
        nthread_ = Communicator::world->nthread();
        scf_type = scf_type_;
        J_[0]->zero();
        Ka_[0]->zero();
        Kb_[0]->zero();

        for (int i=1; i<nthread_; ++i) {
            J_.push_back(SharedMatrix(J_[0]->clone()));
            Ka_.push_back(SharedMatrix(Ka_[0]->clone()));
            Kb_.push_back(SharedMatrix(Kb_[0]->clone()));
        }
    }
    void finalize(){

        for (int i=1; i<nthread_; ++i) {
            J_[0]->add(J_[i]);
            Ka_[0]->add(Ka_[i]);
            Kb_[0]->add(Kb_[i]);
        }
        if (comm != "LOCAL" && scf_type == "DIRECT") {
            for (int i=0; i < J_[0]->nirrep(); i++)
                Communicator::world->sum(J_[0]->get_pointer(i), J_[0]->rowdim(i)*J_[0]->coldim(i));
            for (int i=0; i < Ka_[0]->nirrep(); i++)
                Communicator::world->sum(Ka_[0]->get_pointer(i), Ka_[0]->rowdim(i)*Ka_[0]->coldim(i));
            for (int i=0; i < Kb_[0]->nirrep(); i++)
                Communicator::world->sum(Kb_[0]->get_pointer(i), Kb_[0]->rowdim(i)*Kb_[0]->coldim(i));
        }

        J_[0]->copy_lower_to_upper();
        Ka_[0]->copy_lower_to_upper();
        Kb_[0]->copy_lower_to_upper();
    }

    // NEVER CALL THIS FUNCTION
    J_Ka_Kb_Functor() { throw PSIEXCEPTION("J_Ka_Kb_Functor(): This should never have been called, idiot."); }

    J_Ka_Kb_Functor(boost::shared_ptr<Matrix> J, boost::shared_ptr<Matrix> Ka,
               boost::shared_ptr<Matrix> Kb, const boost::shared_ptr<Matrix> Da, const boost::shared_ptr<Matrix> Db, const boost::shared_ptr<Matrix> Ca, const boost::shared_ptr<Matrix> Cb, const int* Na, const int* Nb)
        : Da_(Da), Db_(Db), Ca_(Ca), Cb_(Cb), Na_(Na), Nb_(Nb)
    {
        J_.push_back(J);
        Ka_.push_back(Ka);
        Kb_.push_back(Kb);
    }

    void operator()(boost::shared_ptr<DFHF> dfhf, boost::shared_ptr<PseudospectralHF> pshf) {
        pshf->set_restricted(false);
        pshf->set_Da(Da_);
        pshf->set_Db(Db_);
        pshf->set_Ka(Ka_[0]);
        pshf->set_Kb(Kb_[0]);
        dfhf->set_restricted(false);
        dfhf->set_jk(false);
        dfhf->set_J(J_[0]);
        dfhf->set_Da(Da_);
        dfhf->set_Db(Db_);
    }

    void operator()(boost::shared_ptr<DFHF> dfhf) {
        dfhf->set_restricted(false);
        dfhf->set_jk(true);
        dfhf->set_J(J_[0]);
        dfhf->set_Da(Da_);
        dfhf->set_Db(Db_);
        dfhf->set_Ca(Ca_);
        dfhf->set_Cb(Cb_);
        dfhf->set_Na(Na_);
        dfhf->set_Nb(Nb_);
        dfhf->set_Ka(Ka_[0]);
        dfhf->set_Kb(Kb_[0]);
    }

    void operator()(boost::shared_ptr<PKIntegrals> pk_integrals) {
        pk_integrals->setup(J_[0], Ka_[0], Kb_[0], Da_, Db_);
    }


    void operator()(int pabs, int qabs, int rabs, int sabs,
                    int psym, int prel, int qsym, int qrel,
                    int rsym, int rrel, int ssym, int srel, double value) {
        int thread = Communicator::world->thread_id(pthread_self());
        /* (pq|rs) */
        if(rsym == ssym){
            J_[thread]->add(rsym, rrel, srel, (Da_->get(psym, prel, qrel) + Db_->get(psym, prel, qrel)) * value);
        }
        if(qabs >= rabs){
            if(qsym == rsym){
                Ka_[thread]->add(qsym, qrel, rrel, Da_->get(psym, prel, srel) * value);
                Kb_[thread]->add(qsym, qrel, rrel, Db_->get(psym, prel, srel) * value);
            }
        }

        if(pabs!=qabs && rabs!=sabs && (pabs!=rabs || qabs!=sabs)){
            /* (pq|sr) */
            if(qabs >= sabs){
                if(qsym == ssym){
                    Ka_[thread]->add(qsym, qrel, srel, Da_->get(psym, prel, rrel) * value);
                    Kb_[thread]->add(qsym, qrel, srel, Db_->get(psym, prel, rrel) * value);
                }
            }

            /* (qp|rs) */
            if(rsym == ssym){
                J_[thread]->add(rsym, rrel, srel, (Da_->get(qsym, qrel, prel) + Db_->get(qsym, qrel, prel)) * value);
            }

            if(pabs >= rabs){
                if(psym == rsym){
                    Ka_[thread]->add(psym, prel, rrel, Da_->get(qsym, qrel, srel) * value);
                    Kb_[thread]->add(psym, prel, rrel, Db_->get(qsym, qrel, srel) * value);
                }
            }

            /* (qp|sr) */
            if(pabs >= sabs){
                if(psym == ssym){
                    Ka_[thread]->add(psym, prel, srel, Da_->get(qsym, qrel, rrel) * value);
                    Kb_[thread]->add(psym, prel, srel, Db_->get(qsym, qrel, rrel) * value);
                }
            }

            /* (rs|pq) */
            if(psym == qsym){
                J_[thread]->add(psym, prel, qrel, (Da_->get(rsym, rrel, srel) + Db_->get(rsym, rrel, srel)) * value);
            }
            if(sabs >= pabs){
                if(ssym == psym){
                    Ka_[thread]->add(ssym, srel, prel, Da_->get(rsym, rrel, qrel) * value);
                    Kb_[thread]->add(ssym, srel, prel, Db_->get(rsym, rrel, qrel) * value);
                }
            }

            /* (sr|pq) */
            if(psym == qsym){
                J_[thread]->add(psym, prel, qrel, (Da_->get(ssym, srel, rrel) + Db_->get(ssym, srel, rrel)) * value);
            }
            if(rabs >= pabs){
                if(rsym == psym){
                    Ka_[thread]->add(rsym, rrel, prel, Da_->get(ssym, srel, qrel) * value);
                    Kb_[thread]->add(rsym, rrel, prel, Db_->get(ssym, srel, qrel) * value);
                }
            }

            /* (rs|qp) */
            if(sabs >= qabs){
                if(ssym == qsym){
                    Ka_[thread]->add(ssym, srel, qrel, Da_->get(rsym, rrel, prel) * value);
                    Kb_[thread]->add(ssym, srel, qrel, Db_->get(rsym, rrel, prel) * value);
                }
            }

            /* (sr|qp) */
            if(rabs >= qabs){
                if(rsym == qsym){
                    Ka_[thread]->add(rsym, rrel, qrel, Da_->get(ssym, srel, prel) * value);
                    Kb_[thread]->add(rsym, rrel, qrel, Db_->get(ssym, srel, prel) * value);
                }
            }
        }else if(pabs!=qabs && rabs!=sabs && pabs==rabs && qabs==sabs){
            /* (pq|sr) */
            if(qabs >= sabs){
                if(qsym == ssym){
                    Ka_[thread]->add(qsym, qrel, srel, Da_->get(psym, prel, rrel) * value);
                    Kb_[thread]->add(qsym, qrel, srel, Db_->get(psym, prel, rrel) * value);
                }
            }
            /* (qp|rs) */
            if(rsym == ssym){
                J_[thread]->add(rsym, rrel, srel, (Da_->get(qsym, qrel, prel) + Db_->get(qsym, qrel, prel)) * value);
            }
            if(pabs >= rabs){
                if(psym == rsym){
                    Ka_[thread]->add(psym, prel, rrel, Da_->get(qsym, qrel, srel) * value);
                    Kb_[thread]->add(psym, prel, rrel, Db_->get(qsym, qrel, srel) * value);
                }
            }

            /* (qp|sr) */
            if(pabs >= sabs){
                if(psym == ssym){
                    Ka_[thread]->add(psym, prel, srel, Da_->get(qsym, qrel, rrel) * value);
                    Kb_[thread]->add(psym, prel, srel, Db_->get(qsym, qrel, rrel) * value);
                }
            }
        }else if(pabs!=qabs && rabs==sabs){
            /* (qp|rs) */
            if(rsym == ssym){
                J_[thread]->add(rsym, rrel, srel, (Da_->get(qsym, qrel, prel) + Db_->get(qsym, qrel, prel)) * value);
            }
            if(pabs >= rabs){
                if(qsym == rsym){
                    Ka_[thread]->add(psym, prel, rrel, Da_->get(qsym, qrel, srel) * value);
                    Kb_[thread]->add(psym, prel, rrel, Db_->get(qsym, qrel, srel) * value);
                }
            }

            /* (rs|pq) */
            if(psym == qsym){
                J_[thread]->add(psym, prel, qrel, (Da_->get(rsym, rrel, srel) + Db_->get(rsym, rrel, srel)) * value);
            }
            if(sabs >= pabs){
                if(ssym == psym){
                    Ka_[thread]->add(ssym, srel, prel, Da_->get(rsym, rrel, qrel) * value);
                    Kb_[thread]->add(ssym, srel, prel, Db_->get(rsym, rrel, qrel) * value);
                }
            }

            /* (rs|qp) */
            if(sabs >= qabs){
                if(ssym == qsym){
                    Ka_[thread]->add(ssym, srel, qrel, Da_->get(rsym, rrel, prel) * value);
                    Kb_[thread]->add(ssym, srel, qrel, Db_->get(rsym, rrel, prel) * value);
                }
            }
        }else if(pabs==qabs && rabs!=sabs){
            /* (pq|sr) */
            if(qabs >= sabs){
                if(qsym == ssym){
                    Ka_[thread]->add(qsym, qrel, srel, Da_->get(psym, prel, rrel) * value);
                    Kb_[thread]->add(qsym, qrel, srel, Db_->get(psym, prel, rrel) * value);
                }
            }

            /* (rs|pq) */
            if(psym == qsym){
                J_[thread]->add(psym, prel, qrel, (Da_->get(rsym, rrel, srel) + Db_->get(rsym, rrel, srel)) * value);
            }
            if(sabs >= pabs){
                if(ssym == psym){
                    Ka_[thread]->add(ssym, srel, prel, Da_->get(rsym, rrel, qrel) * value);
                    Kb_[thread]->add(ssym, srel, prel, Db_->get(rsym, rrel, qrel) * value);
                }
            }

            /* (sr|pq) */
            if(psym == qsym){
                J_[thread]->add(psym, prel, qrel, (Da_->get(ssym, srel, rrel) + Db_->get(ssym, srel, rrel)) * value);
            }
            if(rabs >= pabs){
                if(rsym == psym){
                    Ka_[thread]->add(rsym, rrel, prel, Da_->get(ssym, srel, qrel) * value);
                    Kb_[thread]->add(rsym, rrel, prel, Db_->get(ssym, srel, qrel) * value);
                }
            }
        }else if(pabs==qabs && rabs==sabs && (pabs!=rabs || qabs!=sabs)){
            /* (rs|pq) */
            if(psym == qsym){
                J_[thread]->add(psym, prel, qrel, (Da_->get(rsym, rrel, srel) + Db_->get(rsym, rrel, srel)) * value);
            }
            if(sabs >= pabs){
                if(ssym == psym){
                    Ka_[thread]->add(ssym, srel, prel, Da_->get(rsym, rrel, qrel) * value);
                    Kb_[thread]->add(ssym, srel, prel, Db_->get(rsym, rrel, qrel) * value);
                }
            }
        }
    }
};

template <class JKFunctor>
void HF::process_tei(JKFunctor & functor)
{
    if(scf_type_ == "OUT_OF_CORE"){
        functor.initialize(scf_type_);
        IWL *iwl = new IWL(psio_.get(), PSIF_SO_TEI, integral_threshold_, 1, 1);
        Label *lblptr = iwl->labels();
        Value *valptr = iwl->values();
        int labelIndex, pabs, qabs, rabs, sabs, prel, qrel, rrel, srel, psym, qsym, rsym, ssym;
        double value;
        bool lastBuffer;
        do{
            lastBuffer = iwl->last_buffer();
            for(int index = 0; index < iwl->buffer_count(); ++index){
                labelIndex = 4*index;
                pabs  = abs((int) lblptr[labelIndex++]);
                qabs  = (int) lblptr[labelIndex++];
                rabs  = (int) lblptr[labelIndex++];
                sabs  = (int) lblptr[labelIndex++];
                prel  = so2index_[pabs];
                qrel  = so2index_[qabs];
                rrel  = so2index_[rabs];
                srel  = so2index_[sabs];
                psym  = so2symblk_[pabs];
                qsym  = so2symblk_[qabs];
                rsym  = so2symblk_[rabs];
                ssym  = so2symblk_[sabs];
                value = (double) valptr[index];
                functor(pabs,qabs,rabs,sabs,psym,prel,qsym,qrel,rsym,rrel,ssym,srel,value);
            } /* end loop through current buffer */
            if(!lastBuffer) iwl->fetch();
        }while(!lastBuffer);
        iwl->set_keep_flag(1);
        functor.finalize();
        Communicator::world->sync();
        delete iwl;
    }else if (scf_type_ == "DIRECT"){
        functor.initialize(scf_type_);
        eri_->compute_integrals(functor);  // parallelized
        functor.finalize();
    }else if (scf_type_ == "PSEUDOSPECTRAL"){
        functor.initialize(scf_type_);
        functor(df_, pseudospectral_);
        df_->form_J_DF();
        if(functor.k_required())
            pseudospectral_->form_K_PS();
    }else if (scf_type_ == "DF"){
        functor.initialize(scf_type_);
        functor(df_);
        if(functor.k_required())
            df_->form_JK_DF();
        else
            df_->form_J_DF();
    }else if(scf_type_ == "PK"){
        functor.initialize(scf_type_);
        functor(pk_integrals_);
        if(functor.k_required())
            pk_integrals_->compute_J_and_K();
        else
            pk_integrals_->compute_J();
        functor.finalize();
    }else{
        throw PSIEXCEPTION("SCF_TYPE " + scf_type_ + " is not supported in HF::process_tei");
    }
}

}} // Namespaces

#if HAVE_MADNESS
namespace madness {
namespace archive {

template <class Archive>
struct ArchiveStoreImpl< Archive, psi::scf::J_Functor> {
    static void store(const Archive &ar, const psi::scf::J_Functor &t) {
    }
};

template <class Archive>
struct ArchiveStoreImpl< Archive, psi::scf::J_K_Functor> {
    static void store(const Archive &ar, const psi::scf::J_K_Functor &t) {
    }
};

template <class Archive>
struct ArchiveStoreImpl< Archive, psi::scf::J_Ka_Kb_Functor> {
    static void store(const Archive &ar, const psi::scf::J_Ka_Kb_Functor &t) {
    }
};


template <class Archive>
struct ArchiveLoadImpl< Archive, psi::scf::J_Functor> {
    static void load(const Archive &ar, const psi::scf::J_Functor &t) {
    }
};

template <class Archive>
struct ArchiveLoadImpl< Archive, psi::scf::J_K_Functor> {
    static void load(const Archive &ar, const psi::scf::J_K_Functor &t) {
    }
};

template <class Archive>
struct ArchiveLoadImpl< Archive, psi::scf::J_Ka_Kb_Functor> {
    static void load(const Archive &ar, const psi::scf::J_Ka_Kb_Functor &t) {
    }
};

}
}
#endif // HAVE_MADNESS

#endif // INTEGRALFUNCTORS_H

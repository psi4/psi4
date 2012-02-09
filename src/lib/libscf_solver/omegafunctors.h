#ifndef OMEGAFUNCTORS_H
#define OMEGAFUNCTORS_H

#include "ks.h"
#include <libiwl/iwl.hpp>
#include <libmints/mints.h>

namespace boost {
template<class T> class shared_ptr;
}

// We're assuming that for (PQ|RS) P>=Q, R>=S and PQ>=RS

/*
 * All functors MUST implement the following, or the code won't compile.
 * void initialize() this will always be called at the beginning and should do things
 *                   like clear memory, if needed.  It can just do nothing.
 * void finalize()   this will always be called at the end and should do things like
 *                   free memory and copy data, if needed.  It can just do nothing.
 * void operator()(int pabs, int qabs, int rabs, int sabs,
 *                 int psym, int prel, int qsym, int qrel,
 *                 int rsym, int rrel, int ssym, int srel, double value)
 *                   this is called by the out-of-core and direct codes, and will take
 *                   the relative (within irrep) and absolute values of p, q, r and s, along
 *                   with their symmetries, to generate the long-range wK matrices.
 */


namespace psi{ namespace scf{

/**
  @brief This can be passed into the templated wKS::process_omega_tei() function
         and will compute the long-range exchange matrix for restricted orbitals
         in range-corrected functionals (a.k.a. omega functionals)
*/
class Omega_K_Functor
{
    /// The density matrix
    const SharedMatrix D_;
    /// The occupation matrix
    const SharedMatrix C_;
    /// The number of alpha (and beta) electrons
    const int* N_;
    /// The omega for the iteration
    double omega_;
    /// Dummy matrix for DFHF
    std::vector<SharedMatrix> J_;
    /// The long-range exchange matrix
    std::vector<SharedMatrix> wK_;
    /// The communicator
    std::string comm;
    int nthread_;
    /// The scf type
    std::string scf_type;

public:
    void initialize(std::string scf_type_){
        comm = Communicator::world->communicator();
        nthread_ = Communicator::world->nthread();
        scf_type = scf_type_;
        wK_[0]->zero();
        for (int i=1; i<nthread_; ++i)
            wK_.push_back(SharedMatrix(wK_[0]->clone()));
    }
    void finalize(){
        for (int i=1; i<nthread_; ++i)
            wK_[0]->add(wK_[i]);

        if (comm != "LOCAL" && scf_type == "DIRECT") {
            for (int i=0; i < wK_[0]->nirrep(); i++)
                Communicator::world->sum(wK_[0]->get_pointer(i), wK_[0]->rowdim(i)*wK_[0]->coldim(i));
        }

        wK_[0]->copy_lower_to_upper();
    }

    Omega_K_Functor() { throw PSIEXCEPTION("Omega_K_Functor(): Don't call me, idiot."); }

    Omega_K_Functor(double omega, SharedMatrix wK, const SharedMatrix D,
        const SharedMatrix C, const int* N)
        : omega_(omega), D_(D), C_(C), N_(N)
    {
        J_.clear();
        J_.push_back(SharedMatrix(D->clone()));
        wK_.push_back(wK);
    }

    void operator()(int pabs, int qabs, int rabs, int sabs,
                int psym, int prel, int qsym, int qrel,
                int rsym, int rrel, int ssym, int srel, double value) {
        int thread = Communicator::world->thread_id(pthread_self());
        /* (pq|rs) */
        if(qabs >= rabs){
            if(qsym == rsym){
                wK_[thread]->add(qsym, qrel, rrel, D_->get(psym, prel, srel) * value);
            }
        }

        if(pabs!=qabs && rabs!=sabs && (pabs!=rabs || qabs!=sabs)){
            /* (pq|sr) */
            if(qabs >= sabs){
                if(qsym == ssym){
                    wK_[thread]->add(qsym, qrel, srel, D_->get(psym, prel, rrel) * value);
                }
            }

            /* (qp|rs) */
            if(pabs >= rabs){
                if(psym == rsym){
                    wK_[thread]->add(psym, prel, rrel, D_->get(qsym, qrel, srel) * value);
                }
            }

            /* (qp|sr) */
            if(pabs >= sabs){
                if(psym == ssym){
                    wK_[thread]->add(psym, prel, srel, D_->get(qsym, qrel, rrel) * value);
                }
            }

            /* (rs|pq) */
            if(sabs >= pabs){
                if(ssym == psym){
                    wK_[thread]->add(ssym, srel, prel, D_->get(rsym, rrel, qrel) * value);
                }
            }

            /* (sr|pq) */
            if(rabs >= pabs){
                if(rsym == psym){
                    wK_[thread]->add(rsym, rrel, prel, D_->get(ssym, srel, qrel) * value);
                }
            }

            /* (rs|qp) */
            if(sabs >= qabs){
                if(ssym == qsym){
                    wK_[thread]->add(ssym, srel, qrel, D_->get(rsym, rrel, prel) * value);
                }
            }

            /* (sr|qp) */
            if(rabs >= qabs){
                if(rsym == qsym){
                    wK_[thread]->add(rsym, rrel, qrel, D_->get(ssym, srel, prel) * value);
                }
            }
        }else if(pabs!=qabs && rabs!=sabs && pabs==rabs && qabs==sabs){
            /* (pq|sr) */
            if(qabs >= sabs){
                if(qsym == ssym){
                    wK_[thread]->add(qsym, qrel, srel, D_->get(psym, prel, rrel) * value);
                }
            }
            /* (qp|rs) */
            if(pabs >= rabs){
                if(psym == rsym){
                    wK_[thread]->add(psym, prel, rrel, D_->get(qsym, qrel, srel) * value);
                }
            }

            /* (qp|sr) */
            if(pabs >= sabs){
                if(psym == ssym){
                    wK_[thread]->add(psym, prel, srel, D_->get(qsym, qrel, rrel) * value);
                }
            }
        }else if(pabs!=qabs && rabs==sabs){
            /* (qp|rs) */
            if(pabs >= rabs){
                if(qsym == rsym){
                    wK_[thread]->add(psym, prel, rrel, D_->get(qsym, qrel, srel) * value);
                }
            }

            /* (rs|pq) */
            if(sabs >= pabs){
                if(ssym == psym){
                    wK_[thread]->add(ssym, srel, prel, D_->get(rsym, rrel, qrel) * value);
                }
            }

            /* (rs|qp) */
            if(sabs >= qabs){
                if(ssym == qsym){
                    wK_[thread]->add(ssym, srel, qrel, D_->get(rsym, rrel, prel) * value);
                }
            }
        }else if(pabs==qabs && rabs!=sabs){
            /* (pq|sr) */
            if(qabs >= sabs){
                if(qsym == ssym){
                    wK_[thread]->add(qsym, qrel, srel, D_->get(psym, prel, rrel) * value);
                }
            }

            /* (rs|pq) */
            if(sabs >= pabs){
                if(ssym == psym){
                    wK_[thread]->add(ssym, srel, prel, D_->get(rsym, rrel, qrel) * value);
                }
            }

            /* (sr|pq) */
            if(rabs >= pabs){
                if(rsym == psym){
                    wK_[thread]->add(rsym, rrel, prel, D_->get(ssym, srel, qrel) * value);
                }
            }
        }else if(pabs==qabs && rabs==sabs && (pabs!=rabs || qabs!=sabs)){
            /* (rs|pq) */
            if(sabs >= pabs){
                if(ssym == psym){
                    wK_[thread]->add(ssym, srel, prel, D_->get(rsym, rrel, qrel) * value);
                }
            }
        }
    }
};


/**
  @brief This can be passed into the templated KS::process_omega_tei() function
         and will compute the long-range exchange matrices for restricted orbitals
         in range-corrected functionals (a.k.a. omega functionals)
*/
class Omega_Ka_Kb_Functor
{
    /// The alpha density matrix
    const SharedMatrix Da_;
    /// The beta density matrix
    const SharedMatrix Db_;
    /// The alpha occupation matrix
    const SharedMatrix Ca_;
    /// The beta occupation matrix
    const SharedMatrix Cb_;
    /// The number of alpha electrons
    const int* Na_;
    /// The number of beta electrons
    const int* Nb_;
    /// The omega for the iteration
    double omega_;
    /// Dummy for DFHF
    std::vector<SharedMatrix> J_;
    /// The alpha exchange matrix
    std::vector<SharedMatrix> wKa_;
    /// The beta exchange matrix
    std::vector<SharedMatrix> wKb_;
    /// The communicator
    std::string comm;
    int nthread_;
    /// The scf type
    std::string scf_type;

public:

    void initialize(std::string scf_type_){
        comm = Communicator::world->communicator();
        nthread_ = Communicator::world->nthread();
        scf_type = scf_type_;
        wKa_[0]->zero();
        wKb_[0]->zero();

        for (int i=1; i<nthread_; ++i) {
            wKa_.push_back(SharedMatrix(wKa_[0]->clone()));
            wKb_.push_back(SharedMatrix(wKb_[0]->clone()));
        }
    }

    void finalize(){
        for (int i=1; i<nthread_; ++i) {
            wKa_[0]->add(wKa_[i]);
            wKb_[0]->add(wKb_[i]);
        }
        if (comm != "LOCAL" && scf_type == "DIRECT") {
            for (int i=0; i < wKa_[0]->nirrep(); i++)
                Communicator::world->sum(wKa_[0]->get_pointer(i), wKa_[0]->rowdim(i)*wKa_[0]->coldim(i));
            for (int i=0; i < wKb_[0]->nirrep(); i++)
                Communicator::world->sum(wKb_[0]->get_pointer(i), wKb_[0]->rowdim(i)*wKb_[0]->coldim(i));
        }

        wKa_[0]->copy_lower_to_upper();
        wKb_[0]->copy_lower_to_upper();
    }

    Omega_Ka_Kb_Functor() { throw PSIEXCEPTION("Omega_Ka_Kb_Functor(): Really? You want to do this?"); }

    Omega_Ka_Kb_Functor(double omega, SharedMatrix wKa,
               SharedMatrix wKb, const SharedMatrix Da, const SharedMatrix Db, const SharedMatrix Ca, const SharedMatrix Cb, const int* Na, const int* Nb)
        : omega_(omega), Da_(Da), Db_(Db), Ca_(Ca), Cb_(Cb), Na_(Na), Nb_(Nb)
    {
        J_.clear();
        J_.push_back(SharedMatrix(Da->clone()));
        wKa_.push_back(wKa);
        wKb_.push_back(wKb);
    }

    void operator()(int pabs, int qabs, int rabs, int sabs,
                    int psym, int prel, int qsym, int qrel,
                    int rsym, int rrel, int ssym, int srel, double value) {
        int thread = Communicator::world->thread_id(pthread_self());
        /* (pq|rs) */
        if(qabs >= rabs){
            if(qsym == rsym){
                wKa_[thread]->add(qsym, qrel, rrel, Da_->get(psym, prel, srel) * value);
                wKb_[thread]->add(qsym, qrel, rrel, Db_->get(psym, prel, srel) * value);
            }
        }

        if(pabs!=qabs && rabs!=sabs && (pabs!=rabs || qabs!=sabs)){
            /* (pq|sr) */
            if(qabs >= sabs){
                if(qsym == ssym){
                    wKa_[thread]->add(qsym, qrel, srel, Da_->get(psym, prel, rrel) * value);
                    wKb_[thread]->add(qsym, qrel, srel, Db_->get(psym, prel, rrel) * value);
                }
            }

            /* (qp|rs) */
            if(pabs >= rabs){
                if(psym == rsym){
                    wKa_[thread]->add(psym, prel, rrel, Da_->get(qsym, qrel, srel) * value);
                    wKb_[thread]->add(psym, prel, rrel, Db_->get(qsym, qrel, srel) * value);
                }
            }

            /* (qp|sr) */
            if(pabs >= sabs){
                if(psym == ssym){
                    wKa_[thread]->add(psym, prel, srel, Da_->get(qsym, qrel, rrel) * value);
                    wKb_[thread]->add(psym, prel, srel, Db_->get(qsym, qrel, rrel) * value);
                }
            }

            /* (rs|pq) */
            if(sabs >= pabs){
                if(ssym == psym){
                    wKa_[thread]->add(ssym, srel, prel, Da_->get(rsym, rrel, qrel) * value);
                    wKb_[thread]->add(ssym, srel, prel, Db_->get(rsym, rrel, qrel) * value);
                }
            }

            /* (sr|pq) */
            if(rabs >= pabs){
                if(rsym == psym){
                    wKa_[thread]->add(rsym, rrel, prel, Da_->get(ssym, srel, qrel) * value);
                    wKb_[thread]->add(rsym, rrel, prel, Db_->get(ssym, srel, qrel) * value);
                }
            }

            /* (rs|qp) */
            if(sabs >= qabs){
                if(ssym == qsym){
                    wKa_[thread]->add(ssym, srel, qrel, Da_->get(rsym, rrel, prel) * value);
                    wKb_[thread]->add(ssym, srel, qrel, Db_->get(rsym, rrel, prel) * value);
                }
            }

            /* (sr|qp) */
            if(rabs >= qabs){
                if(rsym == qsym){
                    wKa_[thread]->add(rsym, rrel, qrel, Da_->get(ssym, srel, prel) * value);
                    wKb_[thread]->add(rsym, rrel, qrel, Db_->get(ssym, srel, prel) * value);
                }
            }
        }else if(pabs!=qabs && rabs!=sabs && pabs==rabs && qabs==sabs){
            /* (pq|sr) */
            if(qabs >= sabs){
                if(qsym == ssym){
                    wKa_[thread]->add(qsym, qrel, srel, Da_->get(psym, prel, rrel) * value);
                    wKb_[thread]->add(qsym, qrel, srel, Db_->get(psym, prel, rrel) * value);
                }
            }
            /* (qp|rs) */
            if(pabs >= rabs){
                if(psym == rsym){
                    wKa_[thread]->add(psym, prel, rrel, Da_->get(qsym, qrel, srel) * value);
                    wKb_[thread]->add(psym, prel, rrel, Db_->get(qsym, qrel, srel) * value);
                }
            }

            /* (qp|sr) */
            if(pabs >= sabs){
                if(psym == ssym){
                    wKa_[thread]->add(psym, prel, srel, Da_->get(qsym, qrel, rrel) * value);
                    wKb_[thread]->add(psym, prel, srel, Db_->get(qsym, qrel, rrel) * value);
                }
            }
        }else if(pabs!=qabs && rabs==sabs){
            /* (qp|rs) */
            if(pabs >= rabs){
                if(qsym == rsym){
                    wKa_[thread]->add(psym, prel, rrel, Da_->get(qsym, qrel, srel) * value);
                    wKb_[thread]->add(psym, prel, rrel, Db_->get(qsym, qrel, srel) * value);
                }
            }

            /* (rs|pq) */
            if(sabs >= pabs){
                if(ssym == psym){
                    wKa_[thread]->add(ssym, srel, prel, Da_->get(rsym, rrel, qrel) * value);
                    wKb_[thread]->add(ssym, srel, prel, Db_->get(rsym, rrel, qrel) * value);
                }
            }

            /* (rs|qp) */
            if(sabs >= qabs){
                if(ssym == qsym){
                    wKa_[thread]->add(ssym, srel, qrel, Da_->get(rsym, rrel, prel) * value);
                    wKb_[thread]->add(ssym, srel, qrel, Db_->get(rsym, rrel, prel) * value);
                }
            }
        }else if(pabs==qabs && rabs!=sabs){
            /* (pq|sr) */
            if(qabs >= sabs){
                if(qsym == ssym){
                    wKa_[thread]->add(qsym, qrel, srel, Da_->get(psym, prel, rrel) * value);
                    wKb_[thread]->add(qsym, qrel, srel, Db_->get(psym, prel, rrel) * value);
                }
            }

            /* (rs|pq) */
            if(sabs >= pabs){
                if(ssym == psym){
                    wKa_[thread]->add(ssym, srel, prel, Da_->get(rsym, rrel, qrel) * value);
                    wKb_[thread]->add(ssym, srel, prel, Db_->get(rsym, rrel, qrel) * value);
                }
            }

            /* (sr|pq) */
            if(rabs >= pabs){
                if(rsym == psym){
                    wKa_[thread]->add(rsym, rrel, prel, Da_->get(ssym, srel, qrel) * value);
                    wKb_[thread]->add(rsym, rrel, prel, Db_->get(ssym, srel, qrel) * value);
                }
            }
        }else if(pabs==qabs && rabs==sabs && (pabs!=rabs || qabs!=sabs)){
            /* (rs|pq) */
            if(sabs >= pabs){
                if(ssym == psym){
                    wKa_[thread]->add(ssym, srel, prel, Da_->get(rsym, rrel, qrel) * value);
                    wKb_[thread]->add(ssym, srel, prel, Db_->get(rsym, rrel, qrel) * value);
                }
            }
        }
    }
};

template <class OmegaKFunctor>
void KS::process_omega_tei(OmegaKFunctor & functor)
{
    if (options_.get_str("SCF_TYPE") == "DIRECT"){
        functor.initialize(options_.get_str("SCF_TYPE"));
        omega_eri_->compute_integrals(functor);  // parallelized
        functor.finalize();
    }else{
        throw PSIEXCEPTION("SCF_TYPE " + options_.get_str("SCF_TYPE") + " is not supported in KS::process_omega_tei");
    }
}

}} // Namespaces

#if HAVE_MADNESS
namespace madness {
namespace archive {

template <class Archive>
struct ArchiveStoreImpl< Archive, psi::scf::Omega_K_Functor> {
    static void store(const Archive &ar, const psi::scf::Omega_K_Functor &t) {
    }
};

template <class Archive>
struct ArchiveStoreImpl< Archive, psi::scf::Omega_Ka_Kb_Functor> {
    static void store(const Archive &ar, const psi::scf::Omega_Ka_Kb_Functor &t) {
    }
};


template <class Archive>
struct ArchiveLoadImpl< Archive, psi::scf::Omega_K_Functor> {
    static void load(const Archive &ar, const psi::scf::Omega_K_Functor &t) {
    }
};

template <class Archive>
struct ArchiveLoadImpl< Archive, psi::scf::Omega_Ka_Kb_Functor> {
    static void load(const Archive &ar, const psi::scf::Omega_Ka_Kb_Functor &t) {
    }
};

}
}
#endif // HAVE_MADNESS

#endif // OMEGAFUNCTORS_H

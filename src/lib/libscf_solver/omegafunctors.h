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
    const boost::shared_ptr<Matrix> &D_;
    /// The occupation matrix
    const boost::shared_ptr<Matrix> &C_;
    /// The number of alpha (and beta) electrons
    const int* N_;
    /// The omega for the iteration
    double omega_; 
    /// The long-range exchange matrix
    boost::shared_ptr<Matrix> &wK_;
    /// The communicator
    std::string comm;
    /// The scf type
    std::string scf_type;

public:
    void initialize(std::string comm_, std::string scf_type_){
        comm = comm_;
        scf_type = scf_type_;
        wK_->zero();
    }
    void finalize(){

        if (comm != "LOCAL" && scf_type == "DIRECT") {
            for (int i=0; i < wK_->nirrep(); i++)
                Communicator::world->sum(wK_->get_pointer(i), wK_->rowdim(i)*wK_->coldim(i));
        }

        wK_->copy_lower_to_upper();
    }

    Omega_K_Functor(double omega, boost::shared_ptr<Matrix> wK, const boost::shared_ptr<Matrix> D,
        const boost::shared_ptr<Matrix> C, const int* N)
        : omega_(omega), wK_(wK), D_(D), C_(C), N_(N)
    { }

    void operator()(int pabs, int qabs, int rabs, int sabs,
                int psym, int prel, int qsym, int qrel,
                int rsym, int rrel, int ssym, int srel, double value) {
        /* (pq|rs) */
        if(qabs >= rabs){
            if(qsym == rsym){
                wK_->add(qsym, qrel, rrel, D_->get(psym, prel, srel) * value);
            }
        }

        if(pabs!=qabs && rabs!=sabs && (pabs!=rabs || qabs!=sabs)){
            /* (pq|sr) */
            if(qabs >= sabs){
                if(qsym == ssym){
                    wK_->add(qsym, qrel, srel, D_->get(psym, prel, rrel) * value);
                }
            }

            /* (qp|rs) */
            if(pabs >= rabs){
                if(psym == rsym){
                    wK_->add(psym, prel, rrel, D_->get(qsym, qrel, srel) * value);
                }
            }

            /* (qp|sr) */
            if(pabs >= sabs){
                if(psym == ssym){
                    wK_->add(psym, prel, srel, D_->get(qsym, qrel, rrel) * value);
                }
            }

            /* (rs|pq) */
            if(sabs >= pabs){
                if(ssym == psym){
                    wK_->add(ssym, srel, prel, D_->get(rsym, rrel, qrel) * value);
                }
            }

            /* (sr|pq) */
            if(rabs >= pabs){
                if(rsym == psym){
                    wK_->add(rsym, rrel, prel, D_->get(ssym, srel, qrel) * value);
                }
            }

            /* (rs|qp) */
            if(sabs >= qabs){
                if(ssym == qsym){
                    wK_->add(ssym, srel, qrel, D_->get(rsym, rrel, prel) * value);
                }
            }

            /* (sr|qp) */
            if(rabs >= qabs){
                if(rsym == qsym){
                    wK_->add(rsym, rrel, qrel, D_->get(ssym, srel, prel) * value);
                }
            }
        }else if(pabs!=qabs && rabs!=sabs && pabs==rabs && qabs==sabs){
            /* (pq|sr) */
            if(qabs >= sabs){
                if(qsym == ssym){
                    wK_->add(qsym, qrel, srel, D_->get(psym, prel, rrel) * value);
                }
            }
            /* (qp|rs) */
            if(pabs >= rabs){
                if(psym == rsym){
                    wK_->add(psym, prel, rrel, D_->get(qsym, qrel, srel) * value);
                }
            }

            /* (qp|sr) */
            if(pabs >= sabs){
                if(psym == ssym){
                    wK_->add(psym, prel, srel, D_->get(qsym, qrel, rrel) * value);
                }
            }
        }else if(pabs!=qabs && rabs==sabs){
            /* (qp|rs) */
            if(pabs >= rabs){
                if(qsym == rsym){
                    wK_->add(psym, prel, rrel, D_->get(qsym, qrel, srel) * value);
                }
            }

            /* (rs|pq) */
            if(sabs >= pabs){
                if(ssym == psym){
                    wK_->add(ssym, srel, prel, D_->get(rsym, rrel, qrel) * value);
                }
            }

            /* (rs|qp) */
            if(sabs >= qabs){
                if(ssym == qsym){
                    wK_->add(ssym, srel, qrel, D_->get(rsym, rrel, prel) * value);
                }
            }
        }else if(pabs==qabs && rabs!=sabs){
            /* (pq|sr) */
            if(qabs >= sabs){
                if(qsym == ssym){
                    wK_->add(qsym, qrel, srel, D_->get(psym, prel, rrel) * value);
                }
            }

            /* (rs|pq) */
            if(sabs >= pabs){
                if(ssym == psym){
                    wK_->add(ssym, srel, prel, D_->get(rsym, rrel, qrel) * value);
                }
            }

            /* (sr|pq) */
            if(rabs >= pabs){
                if(rsym == psym){
                    wK_->add(rsym, rrel, prel, D_->get(ssym, srel, qrel) * value);
                }
            }
        }else if(pabs==qabs && rabs==sabs && (pabs!=rabs || qabs!=sabs)){
            /* (rs|pq) */
            if(sabs >= pabs){
                if(ssym == psym){
                    wK_->add(ssym, srel, prel, D_->get(rsym, rrel, qrel) * value);
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
    /// The omega for the iteration
    double omega_; 
    /// The alpha exchange matrix
    boost::shared_ptr<Matrix> &wKa_;
    /// The beta exchange matrix
    boost::shared_ptr<Matrix> &wKb_;
    /// The communicator
    std::string comm;
    /// The scf type
    std::string scf_type;

public:

    void initialize(std::string comm_, std::string scf_type_){
        comm = comm_;
        scf_type = scf_type_;
        wKa_->zero();
        wKb_->zero();
    }
    void finalize(){

        if (comm != "LOCAL" && scf_type == "DIRECT") {
            for (int i=0; i < wKa_->nirrep(); i++)
                Communicator::world->sum(wKa_->get_pointer(i), wKa_->rowdim(i)*wKa_->coldim(i));
            for (int i=0; i < wKb_->nirrep(); i++)
                Communicator::world->sum(wKb_->get_pointer(i), wKb_->rowdim(i)*wKb_->coldim(i));
        }

        wKa_->copy_lower_to_upper();
        wKb_->copy_lower_to_upper();
    }


    Omega_Ka_Kb_Functor(double omega, boost::shared_ptr<Matrix> wKa,
               boost::shared_ptr<Matrix> wKb, const boost::shared_ptr<Matrix> Da, const boost::shared_ptr<Matrix> Db, const boost::shared_ptr<Matrix> Ca, const boost::shared_ptr<Matrix> Cb, const int* Na, const int* Nb)
        : omega_(omega), wKa_(wKa), wKb_(wKb), Da_(Da), Db_(Db), Ca_(Ca), Cb_(Cb), Na_(Na), Nb_(Nb)
    { }

    void operator()(int pabs, int qabs, int rabs, int sabs,
                    int psym, int prel, int qsym, int qrel,
                    int rsym, int rrel, int ssym, int srel, double value) {
        /* (pq|rs) */
        if(qabs >= rabs){
            if(qsym == rsym){
                wKa_->add(qsym, qrel, rrel, Da_->get(psym, prel, srel) * value);
                wKb_->add(qsym, qrel, rrel, Db_->get(psym, prel, srel) * value);
            }
        }

        if(pabs!=qabs && rabs!=sabs && (pabs!=rabs || qabs!=sabs)){
            /* (pq|sr) */
            if(qabs >= sabs){
                if(qsym == ssym){
                    wKa_->add(qsym, qrel, srel, Da_->get(psym, prel, rrel) * value);
                    wKb_->add(qsym, qrel, srel, Db_->get(psym, prel, rrel) * value);
                }
            }

            /* (qp|rs) */
            if(pabs >= rabs){
                if(psym == rsym){
                    wKa_->add(psym, prel, rrel, Da_->get(qsym, qrel, srel) * value);
                    wKb_->add(psym, prel, rrel, Db_->get(qsym, qrel, srel) * value);
                }
            }

            /* (qp|sr) */
            if(pabs >= sabs){
                if(psym == ssym){
                    wKa_->add(psym, prel, srel, Da_->get(qsym, qrel, rrel) * value);
                    wKb_->add(psym, prel, srel, Db_->get(qsym, qrel, rrel) * value);
                }
            }

            /* (rs|pq) */
            if(sabs >= pabs){
                if(ssym == psym){
                    wKa_->add(ssym, srel, prel, Da_->get(rsym, rrel, qrel) * value);
                    wKb_->add(ssym, srel, prel, Db_->get(rsym, rrel, qrel) * value);
                }
            }

            /* (sr|pq) */
            if(rabs >= pabs){
                if(rsym == psym){
                    wKa_->add(rsym, rrel, prel, Da_->get(ssym, srel, qrel) * value);
                    wKb_->add(rsym, rrel, prel, Db_->get(ssym, srel, qrel) * value);
                }
            }

            /* (rs|qp) */
            if(sabs >= qabs){
                if(ssym == qsym){
                    wKa_->add(ssym, srel, qrel, Da_->get(rsym, rrel, prel) * value);
                    wKb_->add(ssym, srel, qrel, Db_->get(rsym, rrel, prel) * value);
                }
            }

            /* (sr|qp) */
            if(rabs >= qabs){
                if(rsym == qsym){
                    wKa_->add(rsym, rrel, qrel, Da_->get(ssym, srel, prel) * value);
                    wKb_->add(rsym, rrel, qrel, Db_->get(ssym, srel, prel) * value);
                }
            }
        }else if(pabs!=qabs && rabs!=sabs && pabs==rabs && qabs==sabs){
            /* (pq|sr) */
            if(qabs >= sabs){
                if(qsym == ssym){
                    wKa_->add(qsym, qrel, srel, Da_->get(psym, prel, rrel) * value);
                    wKb_->add(qsym, qrel, srel, Db_->get(psym, prel, rrel) * value);
                }
            }
            /* (qp|rs) */
            if(pabs >= rabs){
                if(psym == rsym){
                    wKa_->add(psym, prel, rrel, Da_->get(qsym, qrel, srel) * value);
                    wKb_->add(psym, prel, rrel, Db_->get(qsym, qrel, srel) * value);
                }
            }

            /* (qp|sr) */
            if(pabs >= sabs){
                if(psym == ssym){
                    wKa_->add(psym, prel, srel, Da_->get(qsym, qrel, rrel) * value);
                    wKb_->add(psym, prel, srel, Db_->get(qsym, qrel, rrel) * value);
                }
            }
        }else if(pabs!=qabs && rabs==sabs){
            /* (qp|rs) */
            if(pabs >= rabs){
                if(qsym == rsym){
                    wKa_->add(psym, prel, rrel, Da_->get(qsym, qrel, srel) * value);
                    wKb_->add(psym, prel, rrel, Db_->get(qsym, qrel, srel) * value);
                }
            }

            /* (rs|pq) */
            if(sabs >= pabs){
                if(ssym == psym){
                    wKa_->add(ssym, srel, prel, Da_->get(rsym, rrel, qrel) * value);
                    wKb_->add(ssym, srel, prel, Db_->get(rsym, rrel, qrel) * value);
                }
            }

            /* (rs|qp) */
            if(sabs >= qabs){
                if(ssym == qsym){
                    wKa_->add(ssym, srel, qrel, Da_->get(rsym, rrel, prel) * value);
                    wKb_->add(ssym, srel, qrel, Db_->get(rsym, rrel, prel) * value);
                }
            }
        }else if(pabs==qabs && rabs!=sabs){
            /* (pq|sr) */
            if(qabs >= sabs){
                if(qsym == ssym){
                    wKa_->add(qsym, qrel, srel, Da_->get(psym, prel, rrel) * value);
                    wKb_->add(qsym, qrel, srel, Db_->get(psym, prel, rrel) * value);
                }
            }

            /* (rs|pq) */
            if(sabs >= pabs){
                if(ssym == psym){
                    wKa_->add(ssym, srel, prel, Da_->get(rsym, rrel, qrel) * value);
                    wKb_->add(ssym, srel, prel, Db_->get(rsym, rrel, qrel) * value);
                }
            }

            /* (sr|pq) */
            if(rabs >= pabs){
                if(rsym == psym){
                    wKa_->add(rsym, rrel, prel, Da_->get(ssym, srel, qrel) * value);
                    wKb_->add(rsym, rrel, prel, Db_->get(ssym, srel, qrel) * value);
                }
            }
        }else if(pabs==qabs && rabs==sabs && (pabs!=rabs || qabs!=sabs)){
            /* (rs|pq) */
            if(sabs >= pabs){
                if(ssym == psym){
                    wKa_->add(ssym, srel, prel, Da_->get(rsym, rrel, qrel) * value);
                    wKb_->add(ssym, srel, prel, Db_->get(rsym, rrel, qrel) * value);
                }
            }
        }
    }
};

template <class OmegaKFunctor>
void KS::process_omega_tei(OmegaKFunctor & functor)
{
    if (options_.get_str("SCF_TYPE") == "DIRECT"){
        SOShellCombinationsIterator shellIter(sobasisset_, sobasisset_, sobasisset_, sobasisset_);
        std::string comm_ = Process::environment("COMMUNICATOR");
        functor.initialize(comm_, options_.get_str("SCF_TYPE"));
        omega_eri_ao_->setOmega(functional_->getOmega());

        if (comm_ == "LOCAL" || comm_ == "MPI" || comm_ == "GA") {
            int v=0;
            for (shellIter.first(); shellIter.is_done() == false; shellIter.next()) {
                if (Communicator::world->me() == v%Communicator::world->nproc())
                    omega_eri_->compute_shell(shellIter, functor);
                v++;
            }
            timer_on("KS::OmegaFunctor Barrier");
            Communicator::world->sync();
            timer_off("KS::OmegaFunctor Barrier");
        }

        functor.finalize();
    }else{
        throw PSIEXCEPTION("SCF_TYPE " + options_.get_str("SCF_TYPE") + " is not supported in KS::process_omega_tei");
    }
}




}} // Namespaces

#endif // OMEGAFUNCTORS_H

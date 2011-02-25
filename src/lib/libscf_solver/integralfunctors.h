#include  "hf.h"
#include "libiwl/iwl.hpp"
#include "libmints/mints.h"

#ifndef INTEGRALFUNCTORS_H
#define INTEGRALFUNCTORS_H

// If you don't want the restriction that for (PQ|RS) P>=Q, R>=S and
// PQ>=RS, simply define this value to be nonzero
#define NONSTANDARD_ORDERING 0

namespace psi{ namespace scf{

/**
  @brief This can be passed into the templated HF::process_tei() function
         and will compute Coulomb and exchange matrices for restricted orbitals
*/
class J_K_Functor
{
    /// The density matrix
    const shared_ptr<Matrix> D_;
    /// The Coulomb matrix
    shared_ptr<Matrix> &J_;
    /// The exchange matrix
    shared_ptr<Matrix> &K_;

public:
    void initialize(){
        J_->zero();
        K_->zero();
    }
    void finalize(){
        J_->copy_lower_to_upper();
        K_->copy_lower_to_upper();
    }


    J_K_Functor(shared_ptr<Matrix> J, shared_ptr<Matrix> K, const shared_ptr<Matrix> D)
        : J_(J), K_(K), D_(D)
    { }

    void operator()(int pabs, int qabs, int rabs, int sabs,
                    int psym, int prel, int qsym, int qrel,
                    int rsym, int rrel, int ssym, int srel, double value) {
        /* (pq|rs) */
        if(rsym == ssym){
            J_->add(rsym, rrel, srel, D_->get(psym, prel, qrel) * value);
        }
        if(qabs >= rabs){
            if(qsym == rsym){
                K_->add(qsym, qrel, rrel, D_->get(psym, prel, srel) * value);
            }
        }

        if(pabs!=qabs && rabs!=sabs && (pabs!=rabs || qabs!=sabs)){
            /* (pq|sr) */
#if NONSTANDARD_ORDERING
            if(sabs >= rabs){
                if(ssym == rsym){
                    J_->add(ssym, srel, rrel, D_->get(psym, prel, qrel) * value);
                }
            }
#endif
            if(qabs >= sabs){
                if(qsym == ssym){
                    K_->add(qsym, qrel, srel, D_->get(psym, prel, rrel) * value);
                }
            }

            /* (qp|rs) */
#if NONSTANDARD_ORDERING
            if(rabs >= sabs){
#endif
                if(rsym == ssym){
                    J_->add(rsym, rrel, srel, D_->get(qsym, qrel, prel) * value);
                }
#if NONSTANDARD_ORDERING
            }
#endif
            if(pabs >= rabs){
                if(psym == rsym){
                    K_->add(psym, prel, rrel, D_->get(qsym, qrel, srel) * value);
                }
            }

            /* (qp|sr) */
#if NONSTANDARD_ORDERING
            if(sabs >= rabs){
                if(ssym == rsym){
                    J_->add(ssym, srel, rrel,  D_->get(qsym, qrel, prel) * value);
                }
            }
#endif
            if(pabs >= sabs){
                if(psym == ssym){
                    K_->add(psym, prel, srel, D_->get(qsym, qrel, rrel) * value);
                }
            }

            /* (rs|pq) */
#if NONSTANDARD_ORDERING
            if(pabs >= qabs){
#endif
                if(psym == qsym){
                    J_->add(psym, prel, qrel, D_->get(rsym, rrel, srel) * value);
                }
#if NONSTANDARD_ORDERING
            }
#endif
            if(sabs >= pabs){
                if(ssym == psym){
                    K_->add(ssym, srel, prel, D_->get(rsym, rrel, qrel) * value);
                }
            }

            /* (sr|pq) */
#if NONSTANDARD_ORDERING
            if(pabs >= qabs){
#endif
                if(psym == qsym){
                    J_->add(psym, prel, qrel, D_->get(ssym, srel, rrel) * value);
                }
#if NONSTANDARD_ORDERING
            }
#endif
            if(rabs >= pabs){
                if(rsym == psym){
                    K_->add(rsym, rrel, prel, D_->get(ssym, srel, qrel) * value);
                }
            }

            /* (rs|qp) */
#if NONSTANDARD_ORDERING
            if(qabs >= pabs){
                if(qsym == psym){
                    J_->add(qsym, qrel, prel, D_->get(rsym, rrel, srel) * value);
                }
            }
#endif
            if(sabs >= qabs){
                if(ssym == qsym){
                    K_->add(ssym, srel, qrel, D_->get(rsym, rrel, prel) * value);
                }
            }

            /* (sr|qp) */
#if NONSTANDARD_ORDERING
            if(qabs >= pabs){
                if(qsym == psym){
                    J_->add(qsym, qrel, prel, D_->get(ssym, srel, rrel) * value);
                }
            }
#endif
            if(rabs >= qabs){
                if(rsym == qsym){
                    K_->add(rsym, rrel, qrel, D_->get(ssym, srel, prel) * value);
                }
            }
        }else if(pabs!=qabs && rabs!=sabs && pabs==rabs && qabs==sabs){
            /* (pq|sr) */
#if NONSTANDARD_ORDERING
            if(sabs >= rabs){
                if(ssym == rsym){
                    J_->add(ssym, srel, rrel, D_->get(psym, prel, qrel) * value);
                }
            }
#endif
            if(qabs >= sabs){
                if(qsym == ssym){
                    K_->add(qsym, qrel, srel, D_->get(psym, prel, rrel) * value);
                }
            }
            /* (qp|rs) */
#if NONSTANDARD_ORDERING
            if(rabs >= sabs){
#endif
                if(rsym == ssym){
                    J_->add(rsym, rrel, srel, D_->get(qsym, qrel, prel) * value);
                }
#if NONSTANDARD_ORDERING
            }
#endif
            if(pabs >= rabs){
                if(psym == rsym){
                    K_->add(psym, prel, rrel, D_->get(qsym, qrel, srel) * value);
                }
            }

            /* (qp|sr) */
#if NONSTANDARD_ORDERING
            if(sabs >= rabs){
                if(ssym == rsym){
                    J_->add(ssym, srel, rrel, D_->get(qsym, qrel, prel) * value);
                }
            }
#endif
            if(pabs >= sabs){
                if(psym == ssym){
                    K_->add(psym, prel, srel, D_->get(qsym, qrel, rrel) * value);
                }
            }
        }else if(pabs!=qabs && rabs==sabs){
            /* (qp|rs) */
#if NONSTANDARD_ORDERING
            if(rabs >= sabs){
#endif
                if(rsym == ssym){
                    J_->add(rsym, rrel, srel, D_->get(qsym, qrel, prel) * value);
                }
#if NONSTANDARD_ORDERING
            }
#endif
            if(pabs >= rabs){
                if(qsym == rsym){
                    K_->add(psym, prel, rrel, D_->get(qsym, qrel, srel) * value);
                }
            }

            /* (rs|pq) */
#if NONSTANDARD_ORDERING
            if(pabs >= qabs){
#endif
                if(psym == qsym){
                    J_->add(psym, prel, qrel, D_->get(rsym, rrel, srel) * value);
                }
#if NONSTANDARD_ORDERING
            }
#endif
            if(sabs >= pabs){
                if(ssym == psym){
                    K_->add(ssym, srel, prel, D_->get(rsym, rrel, qrel) * value);
                }
            }

            /* (rs|qp) */
#if NONSTANDARD_ORDERING
            if(qabs >= pabs){
                if(qsym == psym){
                    J_->add(qsym, qrel, prel, D_->get(rsym, rrel, srel) * value);
                }
            }
#endif
            if(sabs >= qabs){
                if(ssym == qsym){
                    K_->add(ssym, srel, qrel, D_->get(rsym, rrel, prel) * value);
                }
            }
        }else if(pabs==qabs && rabs!=sabs){
            /* (pq|sr) */
#if NONSTANDARD_ORDERING
            if(sabs >= rabs){
                if(ssym == rsym){
                    J_->add(ssym, srel, rrel, D_->get(psym, prel, qrel) * value);
                }
            }
#endif
            if(qabs >= sabs){
                if(qsym == ssym){
                    K_->add(qsym, qrel, srel, D_->get(psym, prel, rrel) * value);
                }
            }

            /* (rs|pq) */
#if NONSTANDARD_ORDERING
            if(pabs >= qabs){
#endif
                if(psym == qsym){
                    J_->add(psym, prel, qrel, D_->get(rsym, rrel, srel) * value);
                }
#if NONSTANDARD_ORDERING
            }
#endif
            if(sabs >= pabs){
                if(ssym == psym){
                    K_->add(ssym, srel, prel, D_->get(rsym, rrel, qrel) * value);
                }
            }

            /* (sr|pq) */
#if NONSTANDARD_ORDERING
            if(pabs >= qabs){
#endif
                if(psym == qsym){
                    J_->add(psym, prel, qrel, D_->get(ssym, srel, rrel) * value);
                }
#if NONSTANDARD_ORDERING
            }
#endif
            if(rabs >= pabs){
                if(rsym == psym){
                    K_->add(rsym, rrel, prel, D_->get(ssym, srel, qrel) * value);
                }
            }
        }else if(pabs==qabs && rabs==sabs && (pabs!=rabs || qabs!=sabs)){
            /* (rs|pq) */
#if NONSTANDARD_ORDERING
            if(pabs >= qabs){
#endif
                if(psym == qsym){
                    J_->add(psym, prel, qrel, D_->get(rsym, rrel, srel) * value);
                }
#if NONSTANDARD_ORDERING
            }
#endif
            if(sabs >= pabs){
                if(ssym == psym){
                    K_->add(ssym, srel, prel, D_->get(rsym, rrel, qrel) * value);
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
    /// The density matrix
    const shared_ptr<Matrix> D_;
    /// The Coulomb matrix
    shared_ptr<Matrix> &J_;

public:
    void initialize(){
        J_->zero();
    }
    void finalize(){
        J_->copy_lower_to_upper();
    }


    J_Functor(shared_ptr<Matrix> J, const shared_ptr<Matrix> D)
        : J_(J), D_(D)
    { }

    void operator()(int pabs, int qabs, int rabs, int sabs,
                    int psym, int prel, int qsym, int qrel,
                    int rsym, int rrel, int ssym, int srel, double value) {
        /* (pq|rs) */
        if(rsym == ssym){
            J_->add(rsym, rrel, srel, D_->get(psym, prel, qrel) * value);
        }

        if(pabs!=qabs && rabs!=sabs && (pabs!=rabs || qabs!=sabs)){
            /* (pq|sr) */
#if NONSTANDARD_ORDERING
            if(sabs >= rabs){
                if(ssym == rsym){
                    J_->add(ssym, srel, rrel, D_->get(psym, prel, qrel) * value);
                }
            }
#endif

            /* (qp|rs) */
#if NONSTANDARD_ORDERING
            if(rabs >= sabs){
#endif
                if(rsym == ssym){
                    J_->add(rsym, rrel, srel, D_->get(qsym, qrel, prel) * value);
                }
#if NONSTANDARD_ORDERING
            }
#endif

            /* (qp|sr) */
#if NONSTANDARD_ORDERING
            if(sabs >= rabs){
                if(ssym == rsym){
                    J_->add(ssym, srel, rrel,  D_->get(qsym, qrel, prel) * value);
                }
            }
#endif

            /* (rs|pq) */
#if NONSTANDARD_ORDERING
            if(pabs >= qabs){
#endif
                if(psym == qsym){
                    J_->add(psym, prel, qrel, D_->get(rsym, rrel, srel) * value);
                }
#if NONSTANDARD_ORDERING
            }
#endif

            /* (sr|pq) */
#if NONSTANDARD_ORDERING
            if(pabs >= qabs){
#endif
                if(psym == qsym){
                    J_->add(psym, prel, qrel, D_->get(ssym, srel, rrel) * value);
                }
#if NONSTANDARD_ORDERING
            }
#endif

            /* (rs|qp) */
#if NONSTANDARD_ORDERING
            if(qabs >= pabs){
                if(qsym == psym){
                    J_->add(qsym, qrel, prel, D_->get(rsym, rrel, srel) * value);
                }
            }
#endif

            /* (sr|qp) */
#if NONSTANDARD_ORDERING
            if(qabs >= pabs){
                if(qsym == psym){
                    J_->add(qsym, qrel, prel, D_->get(ssym, srel, rrel) * value);
                }
            }
#endif
        }else if(pabs!=qabs && rabs!=sabs && pabs==rabs && qabs==sabs){
            /* (pq|sr) */
#if NONSTANDARD_ORDERING
            if(sabs >= rabs){
                if(ssym == rsym){
                    J_->add(ssym, srel, rrel, D_->get(psym, prel, qrel) * value);
                }
            }
#endif

            /* (qp|rs) */
#if NONSTANDARD_ORDERING
            if(rabs >= sabs){
#endif
                if(rsym == ssym){
                    J_->add(rsym, rrel, srel, D_->get(qsym, qrel, prel) * value);
                }
#if NONSTANDARD_ORDERING
            }
#endif

            /* (qp|sr) */
#if NONSTANDARD_ORDERING
            if(sabs >= rabs){
                if(ssym == rsym){
                    J_->add(ssym, srel, rrel, D_->get(qsym, qrel, prel) * value);
                }
            }
#endif

        }else if(pabs!=qabs && rabs==sabs){
            /* (qp|rs) */
#if NONSTANDARD_ORDERING
            if(rabs >= sabs){
#endif
                if(rsym == ssym){
                    J_->add(rsym, rrel, srel, D_->get(qsym, qrel, prel) * value);
                }
#if NONSTANDARD_ORDERING
            }
#endif

            /* (rs|pq) */
#if NONSTANDARD_ORDERING
            if(pabs >= qabs){
#endif
                if(psym == qsym){
                    J_->add(psym, prel, qrel, D_->get(rsym, rrel, srel) * value);
                }
#if NONSTANDARD_ORDERING
            }
#endif

            /* (rs|qp) */
#if NONSTANDARD_ORDERING
            if(qabs >= pabs){
                if(qsym == psym){
                    J_->add(qsym, qrel, prel, D_->get(rsym, rrel, srel) * value);
                }
            }
#endif

        }else if(pabs==qabs && rabs!=sabs){
            /* (pq|sr) */
#if NONSTANDARD_ORDERING
            if(sabs >= rabs){
                if(ssym == rsym){
                    J_->add(ssym, srel, rrel, D_->get(psym, prel, qrel) * value);
                }
            }
#endif

            /* (rs|pq) */
#if NONSTANDARD_ORDERING
            if(pabs >= qabs){
#endif
                if(psym == qsym){
                    J_->add(psym, prel, qrel, D_->get(rsym, rrel, srel) * value);
                }
#if NONSTANDARD_ORDERING
            }
#endif

            /* (sr|pq) */
#if NONSTANDARD_ORDERING
            if(pabs >= qabs){
#endif
                if(psym == qsym){
                    J_->add(psym, prel, qrel, D_->get(ssym, srel, rrel) * value);
                }
#if NONSTANDARD_ORDERING
            }
#endif

        }else if(pabs==qabs && rabs==sabs && (pabs!=rabs || qabs!=sabs)){
            /* (rs|pq) */
#if NONSTANDARD_ORDERING
            if(pabs >= qabs){
#endif
                if(psym == qsym){
                    J_->add(psym, prel, qrel, D_->get(rsym, rrel, srel) * value);
                }
#if NONSTANDARD_ORDERING
            }
#endif
        }
    }
};


/**
  @brief This can be passed into the templated HF::process_tei() function
         and will compute alpha and beta Coulomb and exchange matrices
*/
class Ja_Jb_Ka_Kb_Functor
{
    /// The alpha density matrix
    const shared_ptr<Matrix> Da_;
    /// The beta density matrix
    const shared_ptr<Matrix> Db_;
    /// The alpha Coulomb matrix
    shared_ptr<Matrix> &Ja_;
    /// The beta Coulomb matrix
    shared_ptr<Matrix> &Jb_;
    /// The alpha exchange matrix
    shared_ptr<Matrix> &Ka_;
    /// The beta exchange matrix
    shared_ptr<Matrix> &Kb_;

public:
    void initialize(){
        Ja_->zero();
        Jb_->zero();
        Ka_->zero();
        Kb_->zero();
    }
    void finalize(){
        Ja_->copy_lower_to_upper();
        Jb_->copy_lower_to_upper();
        Ka_->copy_lower_to_upper();
        Kb_->copy_lower_to_upper();
    }


    Ja_Jb_Ka_Kb_Functor(shared_ptr<Matrix> Ja, shared_ptr<Matrix> Jb, shared_ptr<Matrix> Ka,
               shared_ptr<Matrix> Kb, const shared_ptr<Matrix> Da, const shared_ptr<Matrix> Db)
        : Ja_(Ja), Jb_(Jb), Ka_(Ka), Kb_(Kb), Da_(Da), Db_(Db)
    { }

    void operator()(int pabs, int qabs, int rabs, int sabs,
                    int psym, int prel, int qsym, int qrel,
                    int rsym, int rrel, int ssym, int srel, double value) {
        double temp;

        /* (pq|rs) */
        if(rsym == ssym){
            temp = (Da_->get(psym, prel, qrel) + Db_->get(psym, prel, qrel)) * value;
            Ja_->add(rsym, rrel, srel, temp);
            Jb_->add(rsym, rrel, srel, temp);
        }
        if(qabs >= rabs){
            if(qsym == rsym){
                Ka_->add(qsym, qrel, rrel, Da_->get(psym, prel, srel) * value);
                Kb_->add(qsym, qrel, rrel, Db_->get(psym, prel, srel) * value);
            }
        }

        if(pabs!=qabs && rabs!=sabs && (pabs!=rabs || qabs!=sabs)){
            /* (pq|sr) */
#if NONSTANDARD_ORDERING
            if(sabs >= rabs){
                if(ssym == rsym){
                    temp = (Da_->get(psym, prel, qrel) + Db_->get(psym, prel, qrel)) * value;
                    Ja_->add(ssym, srel, rrel, temp);
                    Jb_->add(ssym, srel, rrel, temp);
                }
            }
#endif
            if(qabs >= sabs){
                if(qsym == ssym){
                    Ka_->add(qsym, qrel, srel, Da_->get(psym, prel, rrel) * value);
                    Kb_->add(qsym, qrel, srel, Db_->get(psym, prel, rrel) * value);
                }
            }

            /* (qp|rs) */
#if NONSTANDARD_ORDERING
            if(rabs >= sabs){
#endif
                if(rsym == ssym){
                    temp = (Da_->get(qsym, qrel, prel) + Db_->get(qsym, qrel, prel)) * value;
                    Ja_->add(rsym, rrel, srel, temp);
                    Jb_->add(rsym, rrel, srel, temp);
                }
#if NONSTANDARD_ORDERING
            }
#endif
            if(pabs >= rabs){
                if(psym == rsym){
                    Ka_->add(psym, prel, rrel, Da_->get(qsym, qrel, srel) * value);
                    Kb_->add(psym, prel, rrel, Db_->get(qsym, qrel, srel) * value);
                }
            }

            /* (qp|sr) */
#if NONSTANDARD_ORDERING
            if(sabs >= rabs){
                if(ssym == rsym){
                    temp = (Da_->get(qsym, qrel, prel) + Db_->get(qsym, qrel, prel)) * value;
                    Ja_->add(ssym, srel, rrel, temp);
                    Jb_->add(ssym, srel, rrel, temp);
                }
            }
#endif
            if(pabs >= sabs){
                if(psym == ssym){
                    Ka_->add(psym, prel, srel, Da_->get(qsym, qrel, rrel) * value);
                    Kb_->add(psym, prel, srel, Db_->get(qsym, qrel, rrel) * value);
                }
            }

            /* (rs|pq) */
#if NONSTANDARD_ORDERING
            if(pabs >= qabs){
#endif
                if(psym == qsym){
                    temp = (Da_->get(rsym, rrel, srel) + Db_->get(rsym, rrel, srel)) * value;
                    Ja_->add(psym, prel, qrel, temp);
                    Jb_->add(psym, prel, qrel, temp);
                }
#if NONSTANDARD_ORDERING
            }
#endif
            if(sabs >= pabs){
                if(ssym == psym){
                    Ka_->add(ssym, srel, prel, Da_->get(rsym, rrel, qrel) * value);
                    Kb_->add(ssym, srel, prel, Db_->get(rsym, rrel, qrel) * value);
                }
            }

            /* (sr|pq) */
#if NONSTANDARD_ORDERING
            if(pabs >= qabs){
#endif
                if(psym == qsym){
                    temp = (Da_->get(ssym, srel, rrel) + Db_->get(ssym, srel, rrel)) * value;
                    Ja_->add(psym, prel, qrel, temp);
                    Jb_->add(psym, prel, qrel, temp);
                }
#if NONSTANDARD_ORDERING
            }
#endif
            if(rabs >= pabs){
                if(rsym == psym){
                    Ka_->add(rsym, rrel, prel, Da_->get(ssym, srel, qrel) * value);
                    Kb_->add(rsym, rrel, prel, Db_->get(ssym, srel, qrel) * value);
                }
            }

            /* (rs|qp) */
#if NONSTANDARD_ORDERING
            if(qabs >= pabs){
                if(qsym == psym){
                    temp = (Da_->get(rsym, rrel, srel) + Db_->get(rsym, rrel, srel)) * value;
                    Ja_->add(qsym, qrel, prel, temp);
                    Jb_->add(qsym, qrel, prel, temp);
                }
            }
#endif
            if(sabs >= qabs){
                if(ssym == qsym){
                    Ka_->add(ssym, srel, qrel, Da_->get(rsym, rrel, prel) * value);
                    Kb_->add(ssym, srel, qrel, Db_->get(rsym, rrel, prel) * value);
                }
            }

            /* (sr|qp) */
#if NONSTANDARD_ORDERING
            if(qabs >= pabs){
                if(qsym == psym){
                    temp = (Da_->get(ssym, srel, rrel) + Db_->get(ssym, srel, rrel)) * value;
                    Ja_->add(qsym, qrel, prel, temp);
                    Jb_->add(qsym, qrel, prel, temp);
                }
            }
#endif
            if(rabs >= qabs){
                if(rsym == qsym){
                    Ka_->add(rsym, rrel, qrel, Da_->get(ssym, srel, prel) * value);
                    Kb_->add(rsym, rrel, qrel, Db_->get(ssym, srel, prel) * value);
                }
            }
        }else if(pabs!=qabs && rabs!=sabs && pabs==rabs && qabs==sabs){
            /* (pq|sr) */
#if NONSTANDARD_ORDERING
            if(sabs >= rabs){
                if(ssym == rsym){
                    temp = (Da_->get(psym, prel, qrel) + Db_->get(psym, prel, qrel)) * value;
                    Ja_->add(ssym, srel, rrel, temp);
                    Jb_->add(ssym, srel, rrel, temp);
                }
            }
#endif
            if(qabs >= sabs){
                if(qsym == ssym){
                    Ka_->add(qsym, qrel, srel, Da_->get(psym, prel, rrel) * value);
                    Kb_->add(qsym, qrel, srel, Db_->get(psym, prel, rrel) * value);
                }
            }
            /* (qp|rs) */
#if NONSTANDARD_ORDERING
            if(rabs >= sabs){
#endif
                if(rsym == ssym){
                    temp = (Da_->get(qsym, qrel, prel) + Db_->get(qsym, qrel, prel)) * value;
                    Ja_->add(rsym, rrel, srel, temp);
                    Jb_->add(rsym, rrel, srel, temp);
                }
#if NONSTANDARD_ORDERING
            }
#endif
            if(pabs >= rabs){
                if(psym == rsym){
                    Ka_->add(psym, prel, rrel, Da_->get(qsym, qrel, srel) * value);
                    Kb_->add(psym, prel, rrel, Db_->get(qsym, qrel, srel) * value);
                }
            }

            /* (qp|sr) */
#if NONSTANDARD_ORDERING
            if(sabs >= rabs){
                if(ssym == rsym){
                    temp = (Da_->get(qsym, qrel, prel) + Db_->get(qsym, qrel, prel)) * value;
                    Ja_->add(ssym, srel, rrel, temp);
                    Jb_->add(ssym, srel, rrel, temp);
                }
            }
#endif
            if(pabs >= sabs){
                if(psym == ssym){
                    Ka_->add(psym, prel, srel, Da_->get(qsym, qrel, rrel) * value);
                    Kb_->add(psym, prel, srel, Db_->get(qsym, qrel, rrel) * value);
                }
            }
        }else if(pabs!=qabs && rabs==sabs){
            /* (qp|rs) */
#if NONSTANDARD_ORDERING
            if(rabs >= sabs){
#endif
                if(rsym == ssym){
                    temp = (Da_->get(qsym, qrel, prel) + Db_->get(qsym, qrel, prel)) * value;
                    Ja_->add(rsym, rrel, srel, temp);
                    Jb_->add(rsym, rrel, srel, temp);
                }
#if NONSTANDARD_ORDERING
            }
#endif
            if(pabs >= rabs){
                if(qsym == rsym){
                    Ka_->add(psym, prel, rrel, Da_->get(qsym, qrel, srel) * value);
                    Kb_->add(psym, prel, rrel, Db_->get(qsym, qrel, srel) * value);
                }
            }

            /* (rs|pq) */
#if NONSTANDARD_ORDERING
            if(pabs >= qabs){
#endif
                if(psym == qsym){
                    temp = (Da_->get(rsym, rrel, srel) + Db_->get(rsym, rrel, srel)) * value;
                    Ja_->add(psym, prel, qrel, temp);
                    Jb_->add(psym, prel, qrel, temp);
                }
#if NONSTANDARD_ORDERING
            }
#endif
            if(sabs >= pabs){
                if(ssym == psym){
                    Ka_->add(ssym, srel, prel, Da_->get(rsym, rrel, qrel) * value);
                    Kb_->add(ssym, srel, prel, Db_->get(rsym, rrel, qrel) * value);
                }
            }

            /* (rs|qp) */
#if NONSTANDARD_ORDERING
            if(qabs >= pabs){
                if(qsym == psym){
                    temp = (Da_->get(rsym, rrel, srel) + Db_->get(rsym, rrel, srel)) * value;
                    Ja_->add(qsym, qrel, prel, temp);
                    Jb_->add(qsym, qrel, prel, temp);
                }
            }
#endif
            if(sabs >= qabs){
                if(ssym == qsym){
                    Ka_->add(ssym, srel, qrel, Da_->get(rsym, rrel, prel) * value);
                    Kb_->add(ssym, srel, qrel, Db_->get(rsym, rrel, prel) * value);
                }
            }
        }else if(pabs==qabs && rabs!=sabs){
            /* (pq|sr) */
#if NONSTANDARD_ORDERING
            if(sabs >= rabs){
                if(ssym == rsym){
                    temp = (Da_->get(psym, prel, qrel) + Db_->get(psym, prel, qrel)) * value;
                    Ja_->add(ssym, srel, rrel, temp);
                    Jb_->add(ssym, srel, rrel, temp);
                }
            }
#endif
            if(qabs >= sabs){
                if(qsym == ssym){
                    Ka_->add(qsym, qrel, srel, Da_->get(psym, prel, rrel) * value);
                    Kb_->add(qsym, qrel, srel, Db_->get(psym, prel, rrel) * value);
                }
            }

            /* (rs|pq) */
#if NONSTANDARD_ORDERING
            if(pabs >= qabs){
#endif
                if(psym == qsym){
                    temp = (Da_->get(rsym, rrel, srel) + Db_->get(rsym, rrel, srel)) * value;
                    Ja_->add(psym, prel, qrel, temp);
                    Jb_->add(psym, prel, qrel, temp);
                }
#if NONSTANDARD_ORDERING
            }
#endif
            if(sabs >= pabs){
                if(ssym == psym){
                    Ka_->add(ssym, srel, prel, Da_->get(rsym, rrel, qrel) * value);
                    Kb_->add(ssym, srel, prel, Db_->get(rsym, rrel, qrel) * value);
                }
            }

            /* (sr|pq) */
#if NONSTANDARD_ORDERING
            if(pabs >= qabs){
#endif
                if(psym == qsym){
                    temp = (Da_->get(ssym, srel, rrel) + Db_->get(ssym, srel, rrel)) * value;
                    Ja_->add(psym, prel, qrel, temp);
                    Jb_->add(psym, prel, qrel, temp);
                }
#if NONSTANDARD_ORDERING
            }
#endif
            if(rabs >= pabs){
                if(rsym == psym){
                    Ka_->add(rsym, rrel, prel, Da_->get(ssym, srel, qrel) * value);
                    Kb_->add(rsym, rrel, prel, Db_->get(ssym, srel, qrel) * value);
                }
            }
        }else if(pabs==qabs && rabs==sabs && (pabs!=rabs || qabs!=sabs)){
            /* (rs|pq) */
#if NONSTANDARD_ORDERING
            if(pabs >= qabs){
#endif
                if(psym == qsym){
                    temp = (Da_->get(rsym, rrel, srel) + Db_->get(rsym, rrel, srel)) * value;
                    Ja_->add(psym, prel, qrel, temp);
                    Jb_->add(psym, prel, qrel, temp);
                }
#if NONSTANDARD_ORDERING
            }
#endif
            if(sabs >= pabs){
                if(ssym == psym){
                    Ka_->add(ssym, srel, prel, Da_->get(rsym, rrel, qrel) * value);
                    Kb_->add(ssym, srel, prel, Db_->get(rsym, rrel, qrel) * value);
                }
            }
        }
    }
};

/**
  @brief This can be passed into the templated HF::process_tei() function
         and will compute just the alpha and beta Coulomb matrices
*/
class Ja_Jb_Functor
{
    /// The alpha density matrix
    const shared_ptr<Matrix> Da_;
    /// The beta density matrix
    const shared_ptr<Matrix> Db_;
    /// The alpha Coulomb matrix
    shared_ptr<Matrix> &Ja_;
    /// The beta Coulomb matrix
    shared_ptr<Matrix> &Jb_;

public:
    void initialize(){
        Ja_->zero();
        Jb_->zero();
    }
    void finalize(){
        Ja_->copy_lower_to_upper();
        Jb_->copy_lower_to_upper();
    }
Ja_Jb_Functor(shared_ptr<Matrix> Ja, shared_ptr<Matrix> Jb,
              const shared_ptr<Matrix> Da, const shared_ptr<Matrix> Db)
    : Ja_(Ja), Jb_(Jb), Da_(Da), Db_(Db)
{ }

void operator()(int pabs, int qabs, int rabs, int sabs,
                int psym, int prel, int qsym, int qrel,
                int rsym, int rrel, int ssym, int srel, double value) {
    double temp;

    /* (pq|rs) */
    if(rsym == ssym){
        temp = (Da_->get(psym, prel, qrel) + Db_->get(psym, prel, qrel)) * value;
        Ja_->add(rsym, rrel, srel, temp);
        Jb_->add(rsym, rrel, srel, temp);
    }

    if(pabs!=qabs && rabs!=sabs && (pabs!=rabs || qabs!=sabs)){
        /* (pq|sr) */
#if NONSTANDARD_ORDERING
        if(sabs >= rabs){
            if(ssym == rsym){
                temp = (Da_->get(psym, prel, qrel) + Db_->get(psym, prel, qrel)) * value;
                Ja_->add(ssym, srel, rrel, temp);
                Jb_->add(ssym, srel, rrel, temp);
            }
        }
#endif

        /* (qp|rs) */
#if NONSTANDARD_ORDERING
        if(rabs >= sabs){
#endif
            if(rsym == ssym){
                temp = (Da_->get(qsym, qrel, prel) + Db_->get(qsym, qrel, prel)) * value;
                Ja_->add(rsym, rrel, srel, temp);
                Jb_->add(rsym, rrel, srel, temp);
            }
#if NONSTANDARD_ORDERING
        }
#endif

        /* (qp|sr) */
#if NONSTANDARD_ORDERING
        if(sabs >= rabs){
            if(ssym == rsym){
                temp = (Da_->get(qsym, qrel, prel) + Db_->get(qsym, qrel, prel)) * value;
                Ja_->add(ssym, srel, rrel, temp);
                Jb_->add(ssym, srel, rrel, temp);
            }
        }
#endif

        /* (rs|pq) */
#if NONSTANDARD_ORDERING
        if(pabs >= qabs){
#endif
            if(psym == qsym){
                temp = (Da_->get(rsym, rrel, srel) + Db_->get(rsym, rrel, srel)) * value;
                Ja_->add(psym, prel, qrel, temp);
                Jb_->add(psym, prel, qrel, temp);
            }
#if NONSTANDARD_ORDERING
        }
#endif

        /* (sr|pq) */
#if NONSTANDARD_ORDERING
        if(pabs >= qabs){
#endif
            if(psym == qsym){
                temp = (Da_->get(ssym, srel, rrel) + Db_->get(ssym, srel, rrel)) * value;
                Ja_->add(psym, prel, qrel, temp);
                Jb_->add(psym, prel, qrel, temp);
            }
#if NONSTANDARD_ORDERING
        }
#endif

        /* (rs|qp) */
#if NONSTANDARD_ORDERING
        if(qabs >= pabs){
            if(qsym == psym){
                temp = (Da_->get(rsym, rrel, srel) + Db_->get(rsym, rrel, srel)) * value;
                Ja_->add(qsym, qrel, prel, temp);
                Jb_->add(qsym, qrel, prel, temp);
            }
        }
#endif

        /* (sr|qp) */
#if NONSTANDARD_ORDERING
        if(qabs >= pabs){
            if(qsym == psym){
                temp = (Da_->get(ssym, srel, rrel) + Db_->get(ssym, srel, rrel)) * value;
                Ja_->add(qsym, qrel, prel, temp);
                Jb_->add(qsym, qrel, prel, temp);
            }
        }
#endif

    }else if(pabs!=qabs && rabs!=sabs && pabs==rabs && qabs==sabs){
        /* (pq|sr) */
#if NONSTANDARD_ORDERING
        if(sabs >= rabs){
            if(ssym == rsym){
                temp = (Da_->get(psym, prel, qrel) + Db_->get(psym, prel, qrel)) * value;
                Ja_->add(ssym, srel, rrel, temp);
                Jb_->add(ssym, srel, rrel, temp);
            }
        }
#endif

        /* (qp|rs) */
#if NONSTANDARD_ORDERING
        if(rabs >= sabs){
#endif
            if(rsym == ssym){
                temp = (Da_->get(qsym, qrel, prel) + Db_->get(qsym, qrel, prel)) * value;
                Ja_->add(rsym, rrel, srel, temp);
                Jb_->add(rsym, rrel, srel, temp);
            }
#if NONSTANDARD_ORDERING
        }
#endif

        /* (qp|sr) */
#if NONSTANDARD_ORDERING
        if(sabs >= rabs){
            if(ssym == rsym){
                temp = (Da_->get(qsym, qrel, prel) + Db_->get(qsym, qrel, prel)) * value;
                Ja_->add(ssym, srel, rrel, temp);
                Jb_->add(ssym, srel, rrel, temp);
            }
        }
#endif

    }else if(pabs!=qabs && rabs==sabs){
        /* (qp|rs) */
#if NONSTANDARD_ORDERING
        if(rabs >= sabs){
#endif
            if(rsym == ssym){
                temp = (Da_->get(qsym, qrel, prel) + Db_->get(qsym, qrel, prel)) * value;
                Ja_->add(rsym, rrel, srel, temp);
                Jb_->add(rsym, rrel, srel, temp);
            }
#if NONSTANDARD_ORDERING
        }
#endif

        /* (rs|pq) */
#if NONSTANDARD_ORDERING
        if(pabs >= qabs){
#endif
            if(psym == qsym){
                temp = (Da_->get(rsym, rrel, srel) + Db_->get(rsym, rrel, srel)) * value;
                Ja_->add(psym, prel, qrel, temp);
                Jb_->add(psym, prel, qrel, temp);
            }
#if NONSTANDARD_ORDERING
        }
#endif

        /* (rs|qp) */
#if NONSTANDARD_ORDERING
        if(qabs >= pabs){
            if(qsym == psym){
                temp = (Da_->get(rsym, rrel, srel) + Db_->get(rsym, rrel, srel)) * value;
                Ja_->add(qsym, qrel, prel, temp);
                Jb_->add(qsym, qrel, prel, temp);
            }
        }
#endif

    }else if(pabs==qabs && rabs!=sabs){
        /* (pq|sr) */
#if NONSTANDARD_ORDERING
        if(sabs >= rabs){
            if(ssym == rsym){
                temp = (Da_->get(psym, prel, qrel) + Db_->get(psym, prel, qrel)) * value;
                Ja_->add(ssym, srel, rrel, temp);
                Jb_->add(ssym, srel, rrel, temp);
            }
        }
#endif

        /* (rs|pq) */
#if NONSTANDARD_ORDERING
        if(pabs >= qabs){
#endif
            if(psym == qsym){
                temp = (Da_->get(rsym, rrel, srel) + Db_->get(rsym, rrel, srel)) * value;
                Ja_->add(psym, prel, qrel, temp);
                Jb_->add(psym, prel, qrel, temp);
            }
#if NONSTANDARD_ORDERING
        }
#endif

        /* (sr|pq) */
#if NONSTANDARD_ORDERING
        if(pabs >= qabs){
#endif
            if(psym == qsym){
                temp = (Da_->get(ssym, srel, rrel) + Db_->get(ssym, srel, rrel)) * value;
                Ja_->add(psym, prel, qrel, temp);
                Jb_->add(psym, prel, qrel, temp);
            }
#if NONSTANDARD_ORDERING
        }
#endif

    }else if(pabs==qabs && rabs==sabs && (pabs!=rabs || qabs!=sabs)){
        /* (rs|pq) */
#if NONSTANDARD_ORDERING
        if(pabs >= qabs){
#endif
            if(psym == qsym){
                temp = (Da_->get(rsym, rrel, srel) + Db_->get(rsym, rrel, srel)) * value;
                Ja_->add(psym, prel, qrel, temp);
                Jb_->add(psym, prel, qrel, temp);
            }
#if NONSTANDARD_ORDERING
        }
#endif
    }
}
};


template <class JKFunctor>
void HF::process_tei(JKFunctor & functor)
{
    if(scf_type_ == "OUT_OF_CORE"){
        functor.initialize();
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
        delete iwl;
    }else if (scf_type_ == "DIRECT"){
        if (!eri_) {
            shared_ptr<IntegralFactory> integral = shared_ptr<IntegralFactory>(new IntegralFactory(basisset_, basisset_, basisset_, basisset_));
            shared_ptr<TwoBodyAOInt> aoeri = shared_ptr<TwoBodyAOInt>(integral->eri());
            eri_ = shared_ptr<TwoBodySOInt>(new TwoBodySOInt(aoeri, integral));
        }

        SOShellCombinationsIterator shellIter(sobasisset_, sobasisset_, sobasisset_, sobasisset_);
        functor.initialize();
        for (shellIter.first(); shellIter.is_done() == false; shellIter.next()) {
            eri_->compute_shell(shellIter, functor);
        }
        functor.finalize();
    }else{
        throw PSIEXCEPTION("SCF_TYPE " + scf_type_ + " is not supported in HF::process_tei");
    }
}




}} // Namespaces


#endif // INTEGRALFUNCTORS_H

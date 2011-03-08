#ifndef INTEGRALFUNCTORS_H
#define INTEGRALFUNCTORS_H

#include  "hf.h"
#include "libiwl/iwl.hpp"
#include "libmints/mints.h"
#include "pseudospectral.h"
#include "pkintegrals.h"
#include "df.h"

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
 * void operator()(shared_ptr<DFHF> dfhf, shared_ptr<PseudospectralHF> psHF) which will set up a PSHF object
 * void operator()(shared_ptr<DSHF> dfHF) which will set up a DFHF object
 * void operator()(shared_ptr<PKIntegrals> pk_integrals) which is used in the PK algorithms.
 */


namespace psi{ namespace scf{

/**
  @brief This can be passed into the templated HF::process_tei() function
         and will compute Coulomb and exchange matrices for restricted orbitals
*/
class J_K_Functor
{
    /// The density matrix
    const shared_ptr<Matrix> D_;
    /// The occupation matrix
    const shared_ptr<Matrix> &C_;
    /// The number of alpha (and beta) electrons 
    const int* N_;
    /// The Coulomb matrix
    shared_ptr<Matrix> &J_;
    /// The exchange matrix
    shared_ptr<Matrix> &K_;

public:
    bool k_required() const {return true;}
    void initialize(){
        J_->zero();
        K_->zero();
    }
    void finalize(){
        J_->copy_lower_to_upper();
        K_->copy_lower_to_upper();
    }


    J_K_Functor(shared_ptr<Matrix> J, shared_ptr<Matrix> K, const shared_ptr<Matrix> D,
        const shared_ptr<Matrix> C, const int* N)
        : J_(J), K_(K), D_(D), C_(C), N_(N)
    { }

    void operator()(shared_ptr<DFHF> dfhf, shared_ptr<PseudospectralHF> pshf) {
        dfhf->set_restricted(true);
        dfhf->set_jk(false);
        dfhf->set_J(J_);
        dfhf->set_Da(D_);
        pshf->set_restricted(true);
        pshf->set_Da(D_);
        pshf->set_Ka(K_);
    }

    void operator()(shared_ptr<PKIntegrals> pk_integrals) {
        pk_integrals->setup(J_, K_, D_, D_);
    }

    void operator()(shared_ptr<DFHF> dfhf) {
        dfhf->set_restricted(true);
        dfhf->set_jk(true);
        dfhf->set_J(J_);
        dfhf->set_Ka(K_);
        dfhf->set_Da(D_);
        dfhf->set_Ca(C_);
        dfhf->set_Na(N_);
    }

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
            if(qabs >= sabs){
                if(qsym == ssym){
                    K_->add(qsym, qrel, srel, D_->get(psym, prel, rrel) * value);
                }
            }

            /* (qp|rs) */
            if(rsym == ssym){
                J_->add(rsym, rrel, srel, D_->get(qsym, qrel, prel) * value);
            }

            if(pabs >= rabs){
                if(psym == rsym){
                    K_->add(psym, prel, rrel, D_->get(qsym, qrel, srel) * value);
                }
            }

            /* (qp|sr) */
            if(pabs >= sabs){
                if(psym == ssym){
                    K_->add(psym, prel, srel, D_->get(qsym, qrel, rrel) * value);
                }
            }

            /* (rs|pq) */
            if(psym == qsym){
                J_->add(psym, prel, qrel, D_->get(rsym, rrel, srel) * value);
            }

            if(sabs >= pabs){
                if(ssym == psym){
                    K_->add(ssym, srel, prel, D_->get(rsym, rrel, qrel) * value);
                }
            }

            /* (sr|pq) */
            if(psym == qsym){
                J_->add(psym, prel, qrel, D_->get(ssym, srel, rrel) * value);
            }

            if(rabs >= pabs){
                if(rsym == psym){
                    K_->add(rsym, rrel, prel, D_->get(ssym, srel, qrel) * value);
                }
            }

            /* (rs|qp) */
            if(sabs >= qabs){
                if(ssym == qsym){
                    K_->add(ssym, srel, qrel, D_->get(rsym, rrel, prel) * value);
                }
            }

            /* (sr|qp) */
            if(rabs >= qabs){
                if(rsym == qsym){
                    K_->add(rsym, rrel, qrel, D_->get(ssym, srel, prel) * value);
                }
            }
        }else if(pabs!=qabs && rabs!=sabs && pabs==rabs && qabs==sabs){
            /* (pq|sr) */
            if(qabs >= sabs){
                if(qsym == ssym){
                    K_->add(qsym, qrel, srel, D_->get(psym, prel, rrel) * value);
                }
            }
            /* (qp|rs) */
            if(rsym == ssym){
                J_->add(rsym, rrel, srel, D_->get(qsym, qrel, prel) * value);
            }
            if(pabs >= rabs){
                if(psym == rsym){
                    K_->add(psym, prel, rrel, D_->get(qsym, qrel, srel) * value);
                }
            }

            /* (qp|sr) */
            if(pabs >= sabs){
                if(psym == ssym){
                    K_->add(psym, prel, srel, D_->get(qsym, qrel, rrel) * value);
                }
            }
        }else if(pabs!=qabs && rabs==sabs){
            /* (qp|rs) */
            if(rsym == ssym){
                J_->add(rsym, rrel, srel, D_->get(qsym, qrel, prel) * value);
            }

            if(pabs >= rabs){
                if(qsym == rsym){
                    K_->add(psym, prel, rrel, D_->get(qsym, qrel, srel) * value);
                }
            }

            /* (rs|pq) */
            if(psym == qsym){
                J_->add(psym, prel, qrel, D_->get(rsym, rrel, srel) * value);
            }

            if(sabs >= pabs){
                if(ssym == psym){
                    K_->add(ssym, srel, prel, D_->get(rsym, rrel, qrel) * value);
                }
            }

            /* (rs|qp) */
            if(sabs >= qabs){
                if(ssym == qsym){
                    K_->add(ssym, srel, qrel, D_->get(rsym, rrel, prel) * value);
                }
            }
        }else if(pabs==qabs && rabs!=sabs){
            /* (pq|sr) */
            if(qabs >= sabs){
                if(qsym == ssym){
                    K_->add(qsym, qrel, srel, D_->get(psym, prel, rrel) * value);
                }
            }

            /* (rs|pq) */
            if(psym == qsym){
                J_->add(psym, prel, qrel, D_->get(rsym, rrel, srel) * value);
            }

            if(sabs >= pabs){
                if(ssym == psym){
                    K_->add(ssym, srel, prel, D_->get(rsym, rrel, qrel) * value);
                }
            }

            /* (sr|pq) */
            if(psym == qsym){
                J_->add(psym, prel, qrel, D_->get(ssym, srel, rrel) * value);
            }

            if(rabs >= pabs){
                if(rsym == psym){
                    K_->add(rsym, rrel, prel, D_->get(ssym, srel, qrel) * value);
                }
            }
        }else if(pabs==qabs && rabs==sabs && (pabs!=rabs || qabs!=sabs)){
            /* (rs|pq) */
            if(psym == qsym){
                J_->add(psym, prel, qrel, D_->get(rsym, rrel, srel) * value);
            }

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
    /// The alpha density matrix
    const shared_ptr<Matrix> Da_;
    /// The beta density matrix
    const shared_ptr<Matrix> Db_;
    /// The Coulomb matrix
    shared_ptr<Matrix> &J_;
    /// Whether this is restricted or not
    bool restricted_;

public:
    bool k_required() const {return false;}

    void initialize(){
        J_->zero();
    }
    void finalize(){
        J_->copy_lower_to_upper();
    }


    J_Functor(shared_ptr<Matrix> J, const shared_ptr<Matrix> Da)
        : J_(J), Da_(Da), Db_(Da)
    {
        restricted_ = true;
    }

    J_Functor(shared_ptr<Matrix> J, const shared_ptr<Matrix> Da, const shared_ptr<Matrix> Db)
        : J_(J), Da_(Da), Db_(Db)
    {
        restricted_ = false;
    }

    void operator()(shared_ptr<DFHF> dfhf, shared_ptr<PseudospectralHF> pshf) {
        // Pseudospectral is exchange only, this does nothing
    }

    void operator()(shared_ptr<DFHF> dfhf) {
        dfhf->set_restricted(restricted_);
        dfhf->set_jk(false);
        dfhf->set_J(J_);
        dfhf->set_Da(Da_);
        dfhf->set_Db(Db_);
    }

    void operator()(shared_ptr<PKIntegrals> pk_integrals) {
        pk_integrals->setup(J_, Da_, Db_);
    }

    void operator()(int pabs, int qabs, int rabs, int sabs,
                    int psym, int prel, int qsym, int qrel,
                    int rsym, int rrel, int ssym, int srel, double value) {
        /* (pq|rs) */
        if(rsym == ssym){
            J_->add(rsym, rrel, srel, Da_->get(psym, prel, qrel) * value);
        }

        if(pabs!=qabs && rabs!=sabs && (pabs!=rabs || qabs!=sabs)){
            /* (pq|sr) */

            /* (qp|rs) */
            if(rsym == ssym){
                J_->add(rsym, rrel, srel, (Da_->get(qsym, qrel, prel) + Db_->get(qsym, qrel, prel)) * value);
            }

            /* (qp|sr) */

            /* (rs|pq) */
            if(psym == qsym){
                J_->add(psym, prel, qrel, (Da_->get(rsym, rrel, srel) + Db_->get(rsym, rrel, srel)) * value);
            }

            /* (sr|pq) */
            if(psym == qsym){
                J_->add(psym, prel, qrel, (Da_->get(ssym, srel, rrel) + Db_->get(ssym, srel, rrel)) * value);
            }

            /* (rs|qp) */

            /* (sr|qp) */
        }else if(pabs!=qabs && rabs!=sabs && pabs==rabs && qabs==sabs){
            /* (pq|sr) */

            /* (qp|rs) */
            if(rsym == ssym){
                J_->add(rsym, rrel, srel, (Da_->get(qsym, qrel, prel) + Db_->get(qsym, qrel, prel)) * value);
            }

            /* (qp|sr) */

        }else if(pabs!=qabs && rabs==sabs){
            /* (qp|rs) */
            if(rsym == ssym){
                J_->add(rsym, rrel, srel, (Da_->get(qsym, qrel, prel) + Db_->get(qsym, qrel, prel)) * value);
            }

            /* (rs|pq) */
            if(psym == qsym){
                J_->add(psym, prel, qrel, (Da_->get(rsym, rrel, srel) + Db_->get(rsym, rrel, srel)) * value);
            } 
            /* (rs|qp) */

        }else if(pabs==qabs && rabs!=sabs){
            /* (pq|sr) */

            /* (rs|pq) */
            if(psym == qsym){
                J_->add(psym, prel, qrel, (Da_->get(rsym, rrel, srel) + Db_->get(rsym, rrel, srel)) * value);
            }

            /* (sr|pq) */
            if(psym == qsym){
                J_->add(psym, prel, qrel, (Da_->get(ssym, srel, rrel) + Db_->get(ssym, srel, rrel)) * value);
            }
        }else if(pabs==qabs && rabs==sabs && (pabs!=rabs || qabs!=sabs)){
            /* (rs|pq) */
            if(psym == qsym){
                J_->add(psym, prel, qrel, (Da_->get(rsym, rrel, srel) + Db_->get(rsym, rrel, srel)) * value);
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
    const shared_ptr<Matrix> Da_;
    /// The beta density matrix
    const shared_ptr<Matrix> Db_;
    /// The alpha occupation matrix
    const shared_ptr<Matrix> Ca_;
    /// The beta occupation matrix
    const shared_ptr<Matrix> Cb_;
    /// The number of alpha electrons 
    const int* Na_;
    /// The number of beta electrons 
    const int* Nb_;
    /// The alpha Coulomb matrix
    shared_ptr<Matrix> &J_;
    /// The alpha exchange matrix
    shared_ptr<Matrix> &Ka_;
    /// The beta exchange matrix
    shared_ptr<Matrix> &Kb_;

public:
    bool k_required() const {return true;}

    void initialize(){
        J_->zero();
        Ka_->zero();
        Kb_->zero();
    }
    void finalize(){
        J_->copy_lower_to_upper();
        Ka_->copy_lower_to_upper();
        Kb_->copy_lower_to_upper();
    }


    J_Ka_Kb_Functor(shared_ptr<Matrix> J, shared_ptr<Matrix> Ka,
               shared_ptr<Matrix> Kb, const shared_ptr<Matrix> Da, const shared_ptr<Matrix> Db, const shared_ptr<Matrix> Ca, const shared_ptr<Matrix> Cb, const int* Na, const int* Nb)
        : J_(J), Ka_(Ka), Kb_(Kb), Da_(Da), Db_(Db), Ca_(Ca), Cb_(Cb), Na_(Na), Nb_(Nb)
    { }

    void operator()(shared_ptr<DFHF> dfhf, shared_ptr<PseudospectralHF> pshf) {
        pshf->set_restricted(false);
        pshf->set_Da(Da_);
        pshf->set_Db(Db_);
        pshf->set_Ka(Ka_);
        pshf->set_Kb(Kb_);
        dfhf->set_restricted(false);
        dfhf->set_jk(false);
        dfhf->set_J(J_);
        dfhf->set_Da(Da_);
        dfhf->set_Db(Db_);
    }

    void operator()(shared_ptr<DFHF> dfhf) {
        dfhf->set_restricted(false);
        dfhf->set_jk(true);
        dfhf->set_J(J_);
        dfhf->set_Da(Da_);
        dfhf->set_Db(Db_);
        dfhf->set_Ca(Ca_);
        dfhf->set_Cb(Cb_);
        dfhf->set_Na(Na_);
        dfhf->set_Nb(Nb_);
        dfhf->set_Ka(Ka_);
        dfhf->set_Kb(Kb_);
    }

    void operator()(shared_ptr<PKIntegrals> pk_integrals) {
        pk_integrals->setup(J_, Ka_, Kb_, Da_, Db_);
    }


    void operator()(int pabs, int qabs, int rabs, int sabs,
                    int psym, int prel, int qsym, int qrel,
                    int rsym, int rrel, int ssym, int srel, double value) {
        /* (pq|rs) */
        if(rsym == ssym){
            J_->add(rsym, rrel, srel, (Da_->get(psym, prel, qrel) + Db_->get(psym, prel, qrel)) * value);
        }
        if(qabs >= rabs){
            if(qsym == rsym){
                Ka_->add(qsym, qrel, rrel, Da_->get(psym, prel, srel) * value);
                Kb_->add(qsym, qrel, rrel, Db_->get(psym, prel, srel) * value);
            }
        }

        if(pabs!=qabs && rabs!=sabs && (pabs!=rabs || qabs!=sabs)){
            /* (pq|sr) */
            if(qabs >= sabs){
                if(qsym == ssym){
                    Ka_->add(qsym, qrel, srel, Da_->get(psym, prel, rrel) * value);
                    Kb_->add(qsym, qrel, srel, Db_->get(psym, prel, rrel) * value);
                }
            }

            /* (qp|rs) */
            if(rsym == ssym){
                J_->add(rsym, rrel, srel, (Da_->get(qsym, qrel, prel) + Db_->get(qsym, qrel, prel)) * value);
            }

            if(pabs >= rabs){
                if(psym == rsym){
                    Ka_->add(psym, prel, rrel, Da_->get(qsym, qrel, srel) * value);
                    Kb_->add(psym, prel, rrel, Db_->get(qsym, qrel, srel) * value);
                }
            }

            /* (qp|sr) */
            if(pabs >= sabs){
                if(psym == ssym){
                    Ka_->add(psym, prel, srel, Da_->get(qsym, qrel, rrel) * value);
                    Kb_->add(psym, prel, srel, Db_->get(qsym, qrel, rrel) * value);
                }
            }

            /* (rs|pq) */
            if(psym == qsym){
                J_->add(psym, prel, qrel, (Da_->get(rsym, rrel, srel) + Db_->get(rsym, rrel, srel)) * value);
            }
            if(sabs >= pabs){
                if(ssym == psym){
                    Ka_->add(ssym, srel, prel, Da_->get(rsym, rrel, qrel) * value);
                    Kb_->add(ssym, srel, prel, Db_->get(rsym, rrel, qrel) * value);
                }
            }

            /* (sr|pq) */
            if(psym == qsym){
                J_->add(psym, prel, qrel, (Da_->get(ssym, srel, rrel) + Db_->get(ssym, srel, rrel)) * value);
            }
            if(rabs >= pabs){
                if(rsym == psym){
                    Ka_->add(rsym, rrel, prel, Da_->get(ssym, srel, qrel) * value);
                    Kb_->add(rsym, rrel, prel, Db_->get(ssym, srel, qrel) * value);
                }
            }

            /* (rs|qp) */
            if(sabs >= qabs){
                if(ssym == qsym){
                    Ka_->add(ssym, srel, qrel, Da_->get(rsym, rrel, prel) * value);
                    Kb_->add(ssym, srel, qrel, Db_->get(rsym, rrel, prel) * value);
                }
            }

            /* (sr|qp) */
            if(rabs >= qabs){
                if(rsym == qsym){
                    Ka_->add(rsym, rrel, qrel, Da_->get(ssym, srel, prel) * value);
                    Kb_->add(rsym, rrel, qrel, Db_->get(ssym, srel, prel) * value);
                }
            }
        }else if(pabs!=qabs && rabs!=sabs && pabs==rabs && qabs==sabs){
            /* (pq|sr) */
            if(qabs >= sabs){
                if(qsym == ssym){
                    Ka_->add(qsym, qrel, srel, Da_->get(psym, prel, rrel) * value);
                    Kb_->add(qsym, qrel, srel, Db_->get(psym, prel, rrel) * value);
                }
            }
            /* (qp|rs) */
            if(rsym == ssym){
                J_->add(rsym, rrel, srel, (Da_->get(qsym, qrel, prel) + Db_->get(qsym, qrel, prel)) * value);
            }
            if(pabs >= rabs){
                if(psym == rsym){
                    Ka_->add(psym, prel, rrel, Da_->get(qsym, qrel, srel) * value);
                    Kb_->add(psym, prel, rrel, Db_->get(qsym, qrel, srel) * value);
                }
            }

            /* (qp|sr) */
            if(pabs >= sabs){
                if(psym == ssym){
                    Ka_->add(psym, prel, srel, Da_->get(qsym, qrel, rrel) * value);
                    Kb_->add(psym, prel, srel, Db_->get(qsym, qrel, rrel) * value);
                }
            }
        }else if(pabs!=qabs && rabs==sabs){
            /* (qp|rs) */
            if(rsym == ssym){
                J_->add(rsym, rrel, srel, (Da_->get(qsym, qrel, prel) + Db_->get(qsym, qrel, prel)) * value);
            }
            if(pabs >= rabs){
                if(qsym == rsym){
                    Ka_->add(psym, prel, rrel, Da_->get(qsym, qrel, srel) * value);
                    Kb_->add(psym, prel, rrel, Db_->get(qsym, qrel, srel) * value);
                }
            }

            /* (rs|pq) */
            if(psym == qsym){
                J_->add(psym, prel, qrel, (Da_->get(rsym, rrel, srel) + Db_->get(rsym, rrel, srel)) * value);
            }
            if(sabs >= pabs){
                if(ssym == psym){
                    Ka_->add(ssym, srel, prel, Da_->get(rsym, rrel, qrel) * value);
                    Kb_->add(ssym, srel, prel, Db_->get(rsym, rrel, qrel) * value);
                }
            }

            /* (rs|qp) */
            if(sabs >= qabs){
                if(ssym == qsym){
                    Ka_->add(ssym, srel, qrel, Da_->get(rsym, rrel, prel) * value);
                    Kb_->add(ssym, srel, qrel, Db_->get(rsym, rrel, prel) * value);
                }
            }
        }else if(pabs==qabs && rabs!=sabs){
            /* (pq|sr) */
            if(qabs >= sabs){
                if(qsym == ssym){
                    Ka_->add(qsym, qrel, srel, Da_->get(psym, prel, rrel) * value);
                    Kb_->add(qsym, qrel, srel, Db_->get(psym, prel, rrel) * value);
                }
            }

            /* (rs|pq) */
            if(psym == qsym){
                J_->add(psym, prel, qrel, (Da_->get(rsym, rrel, srel) + Db_->get(rsym, rrel, srel)) * value);
            }
            if(sabs >= pabs){
                if(ssym == psym){
                    Ka_->add(ssym, srel, prel, Da_->get(rsym, rrel, qrel) * value);
                    Kb_->add(ssym, srel, prel, Db_->get(rsym, rrel, qrel) * value);
                }
            }

            /* (sr|pq) */
            if(psym == qsym){
                J_->add(psym, prel, qrel, (Da_->get(ssym, srel, rrel) + Db_->get(ssym, srel, rrel)) * value);
            }
            if(rabs >= pabs){
                if(rsym == psym){
                    Ka_->add(rsym, rrel, prel, Da_->get(ssym, srel, qrel) * value);
                    Kb_->add(rsym, rrel, prel, Db_->get(ssym, srel, qrel) * value);
                }
            }
        }else if(pabs==qabs && rabs==sabs && (pabs!=rabs || qabs!=sabs)){
            /* (rs|pq) */
            if(psym == qsym){
                J_->add(psym, prel, qrel, (Da_->get(rsym, rrel, srel) + Db_->get(rsym, rrel, srel)) * value);
            }
            if(sabs >= pabs){
                if(ssym == psym){
                    Ka_->add(ssym, srel, prel, Da_->get(rsym, rrel, qrel) * value);
                    Kb_->add(ssym, srel, prel, Db_->get(rsym, rrel, qrel) * value);
                }
            }
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
        SOShellCombinationsIterator shellIter(sobasisset_, sobasisset_, sobasisset_, sobasisset_);
        functor.initialize();
        for (shellIter.first(); shellIter.is_done() == false; shellIter.next()) {
            eri_->compute_shell(shellIter, functor);
        }
        functor.finalize();
    }else if (scf_type_ == "PSEUDOSPECTRAL"){
        functor.initialize();
        functor(df_, pseudospectral_);
        df_->form_J_DF();
        if(functor.k_required())
            pseudospectral_->form_K_PS();
    }else if (scf_type_ == "DF"){
        functor.initialize();
        functor(df_);
        if(functor.k_required())
            df_->form_JK_DF();
        else 
            df_->form_J_DF();
    }else if(scf_type_ == "PK"){
        functor.initialize();
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


#endif // INTEGRALFUNCTORS_H

#include "pkintegrals.h"
#include <libiwl/iwl.hpp>
#include <liboptions/liboptions.h>
#include <exception.h>
#include <libmints/matrix.h>
#include <libpsio/psio.hpp>
#include <boost/shared_ptr.hpp>
#include <psifiles.h>

#define INDEX2(i,j) (i < j ? (j)*((j)+1)/2 + i : (i)*((i)+1)/2 + j)

using namespace boost;

namespace psi{ namespace scf{

PKIntegrals::PKIntegrals(size_t memory, boost::shared_ptr<PSIO> psio, Options& options,
                         int nirreps, const int* sopi, const int *so2index, const int *so2symblk):
        memory_(memory),
        psio_(psio),
        options_(options),
        nirreps_(nirreps),
        sopi_(sopi),
        so2index_(so2index),
        so2symblk_(so2symblk),
        pk_initialized_(false),
        p_k_(0),
        p_j_(0)
{
    print_ = options_.get_int("PRINT");
}

void
PKIntegrals::finalize()
{
    if(p_j_) delete [] p_j_;
    if(p_k_) delete [] p_k_;
    p_j_ = p_k_ = 0;
}

PKIntegrals::~PKIntegrals()
{
    finalize();
}

/**
 * Sets this object up for use with a J functor.
 * @param J a shared pointer to a Matrix, which will hold the resulting J
 * @param Da a shared pointer to an alpha density matrix
 * @param Db a shared pointer to a beta density matrix
 */
void
PKIntegrals::setup(boost::shared_ptr<Matrix> J, const boost::shared_ptr<Matrix> Da, const boost::shared_ptr<Matrix> Db)
{
    J_  = J;
    Da_ = Da;
    Db_ = Db;
    restricted_ = false;
    setup_arrays(false);
}

/**
 * Sets this object up for use with a J_K functor.
 * @param J a shared pointer to a Matrix, which will hold the resulting J
 * @param K a shared pointer to a Matrix, which will hold the resulting K
 * @param D a shared pointer to the density matrix
 */
void
PKIntegrals::setup(boost::shared_ptr<Matrix> J, boost::shared_ptr<Matrix> K,  const boost::shared_ptr<Matrix> Da, const boost::shared_ptr<Matrix> Db)
{
    J_  = J;
    Ka_ = K;
    Da_ = Da;
    Db_ = Db;
    restricted_ = true;
    setup_arrays(true);
}

/**
 * Sets this object up for use with a J_Ka_Kb functor.
 * @param J a shared pointer to a Matrix, which will hold the resulting J
 * @param Ka a shared pointer to a Matrix, which will hold the resulting alpha K
 * @param Kb a shared pointer to a Matrix, which will hold the resulting beta K
 * @param Da a shared pointer to an alpha density matrix
 * @param Db a shared pointer to a beta density matrix
 */
void
PKIntegrals::setup(boost::shared_ptr<Matrix> J, boost::shared_ptr<Matrix> Ka, boost::shared_ptr<Matrix> Kb,
                   const boost::shared_ptr<Matrix> Da, const boost::shared_ptr<Matrix> Db)
{
    J_  = J;
    Ka_ = Ka;
    Kb_ = Kb;
    Da_ = Da;
    Db_ = Db;
    restricted_ = false;
    setup_arrays(true);
}

void
PKIntegrals::process_J_Ka_Kb_block(double *J_block, double *K_block, size_t start, size_t end,
                                   double *&Da_vector,  double *&Db_vector,
                                   double *&J_vector, double *&Ka_vector, double *&Kb_vector)
{
    double* Da_rs;
    double* Db_rs;
    double* J_rs;
    double* Ka_rs;
    double* Kb_rs;
    double Da_pq;
    double Db_pq;
    double J_pq;
    double Ka_pq;
    double Kb_pq;
    for (size_t pq = start; pq < end; ++pq) {
        Da_pq = Da_vector[pq];
        Da_rs = Da_vector;
        Db_pq = Db_vector[pq];
        Db_rs = Db_vector;
        J_pq = 0.0;
        J_rs = J_vector;
        Ka_pq = 0.0;
        Kb_pq = 0.0;
        Ka_rs = Ka_vector;
        Kb_rs = Kb_vector;
        for (size_t rs = 0; rs <= pq; ++rs) {
            J_pq   += *J_block * (*Da_rs + *Db_rs);
            *J_rs  += *J_block * (Da_pq + Db_pq);
            Ka_pq  += *K_block * (*Da_rs);
            *Ka_rs += *K_block * Da_pq;
            Kb_pq  += *K_block * (*Db_rs);
            *Kb_rs += *K_block * Db_pq;
            ++Da_rs;
            ++Db_rs;
            ++J_rs;
            ++Ka_rs;
            ++Kb_rs;
            ++J_block;
            ++K_block;
        }
        J_vector[pq]  += J_pq;
        Ka_vector[pq] += Ka_pq;
        Kb_vector[pq] += Kb_pq;
    }
}

void
PKIntegrals::process_J_K_block(double *J_block, double *K_block, size_t start, size_t end,
                               double *&Da_vector, double *&J_vector,  double *&Ka_vector)
{
    double* Da_rs;
    double* J_rs;
    double* Ka_rs;
    double Da_pq;
    double J_pq;
    double Ka_pq;
    for (size_t pq = start; pq < end; ++pq) {
        Da_pq = Da_vector[pq];
        Da_rs = Da_vector;
        J_pq = 0.0;
        J_rs = J_vector;
        Ka_pq = 0.0;
        Ka_rs = Ka_vector;
        for (size_t rs = 0; rs <= pq; ++rs) {
            J_pq   += *J_block * (*Da_rs);
            *J_rs  += *J_block * Da_pq;
            Ka_pq  += *K_block * (*Da_rs);
            *Ka_rs += *K_block * Da_pq;
            ++Da_rs;
            ++J_rs;
            ++Ka_rs;
            ++J_block;
            ++K_block;
        }
        J_vector[pq]  += J_pq;
        Ka_vector[pq] += Ka_pq;
    }
}

void
PKIntegrals::process_J_block(double *J_block, size_t start, size_t end, double *&Da_vector,
                             double *&Db_vector, double *&J_vector)
{
    double* Da_rs;
    double* Db_rs;
    double* J_rs;
    double Da_pq;
    double Db_pq;
    double J_pq;
    for (size_t pq = start; pq < end; ++pq) {
        Da_pq = Da_vector[pq];
        Da_rs = Da_vector;
        Db_pq = Db_vector[pq];
        Db_rs = Db_vector;
        J_pq = 0.0;
        J_rs = J_vector;
        for (size_t rs = 0; rs <= pq; ++rs) {
            J_pq   += *J_block * (*Da_rs + *Db_rs);
            *J_rs  += *J_block * (Da_pq + Db_pq);
            ++Da_rs;
            ++Db_rs;
            ++J_rs;
            ++J_block;
        }
        J_vector[pq]  += J_pq;
    }
}


/**
 * Computes the J matrix, using the P_J matrix.  Assumes that the appropriate setup() method
 * has already been called, to set up the matrices needed.
 */
void
PKIntegrals::compute_J()
{
    double *J_vector;
    double *Da_vector;
    double *Db_vector;
    Da_vector = new double[pk_pairs_];
    Db_vector = new double[pk_pairs_];
    J_vector  = new double[pk_pairs_];
    ::memset(J_vector,  0, pk_pairs_ * sizeof(double));

    // The off-diagonal terms need to be doubled here
    size_t pqval = 0;
    for (int h = 0; h < nirreps_; ++h) {
        for (int p = 0; p < sopi_[h]; ++p) {
            for (int q = 0; q <= p; ++q) {
                if (p != q) {
                    Da_vector[pqval] = 2.0 * Da_->get(h, p, q);
                    Db_vector[pqval] = 2.0 * Db_->get(h, p, q);
                }else{
                    Da_vector[pqval] =  Da_->get(h, p, q);
                    Db_vector[pqval] =  Db_->get(h, p, q);
                }
                ++pqval;
            }
        }
    }

    process_J_block(p_j_, 0, pk_pairs_, Da_vector, Db_vector, J_vector);

    pqval = 0;
    for (int h = 0; h < nirreps_; ++h) {
        for (int p = 0; p < sopi_[h]; ++p) {
            for (int q = 0; q <= p; ++q) {
                J_->set(h, p, q, J_vector[pqval]);
                ++pqval;
            }
        }
    }

    delete [] J_vector;
}


/**
 * Computes the J and K matrices, using the P_J and P_K matrices.  Assumes that the appropriate
 * setup() method has already been called, to set up the matrices needed.
 */
void
PKIntegrals::compute_J_and_K()
{
    double *J_vector;
    double *Ka_vector;
    double *Kb_vector;
    double *Da_vector;
    double *Db_vector;
    Da_vector = new double[pk_pairs_];
    Db_vector = new double[pk_pairs_];
    J_vector = new double[pk_pairs_];
    ::memset(J_vector,  0, pk_pairs_ * sizeof(double));
    Ka_vector = new double[pk_pairs_];
    ::memset(Ka_vector, 0, pk_pairs_ * sizeof(double));
    if(!restricted_){
        Kb_vector = new double[pk_pairs_];
        ::memset(Kb_vector, 0, pk_pairs_ * sizeof(double));
    }

    // The off-diagonal terms need to be doubled here
    size_t pqval = 0;
    for (int h = 0; h < nirreps_; ++h) {
        for (int p = 0; p < sopi_[h]; ++p) {
            for (int q = 0; q <= p; ++q) {
                if (p != q) {
                    Da_vector[pqval] = 2.0 * Da_->get(h, p, q);
                    Db_vector[pqval] = 2.0 * Db_->get(h, p, q);
                }else{
                    Da_vector[pqval] =  Da_->get(h, p, q);
                    Db_vector[pqval] =  Db_->get(h, p, q);
                }
                ++pqval;
            }
        }
    }

    if(restricted_)
        process_J_K_block(p_j_, p_k_, 0, pk_pairs_, Da_vector, J_vector, Ka_vector);
    else
        process_J_Ka_Kb_block(p_j_, p_k_, 0, pk_pairs_, Da_vector, Db_vector, J_vector, Ka_vector, Kb_vector);

    double prefactor = restricted_ ? 2.0 : 1.0;
    pqval = 0;
    for (int h = 0; h < nirreps_; ++h) {
        for (int p = 0; p < sopi_[h]; ++p) {
            for (int q = 0; q <= p; ++q) {
                J_->set(h, p, q, prefactor * J_vector[pqval]);
                Ka_->set(h, p, q, Ka_vector[pqval]);
                if(!restricted_) Kb_->set(h, p, q, Kb_vector[pqval]);
                ++pqval;
            }
        }
    }

    delete [] J_vector;
    delete [] Ka_vector;
    if(!restricted_)
        delete [] Kb_vector;
}

/**
 * Makes the PJ (Coulomb-like) and, optionally, the PK (exchange-like)
 * matrices needed in SCF procedures.
 * @param Whether the K matrix is to be built along with J
 */
void PKIntegrals::setup_arrays(bool build_k)
{
    if(pk_initialized_)
        return;
    pk_initialized_ = true;

    // Compute PK symmetry mapping
    int *pk_symoffset = new int [nirreps_];

    pk_size_ = 0; pk_pairs_ = 0;
    for (int h=0; h<nirreps_; ++h) {
        pk_symoffset[h] = pk_pairs_;
        // Add up possible pair combinations that yield A1 symmetry
        pk_pairs_ += sopi_[h]*(sopi_[h] + 1)/2;
    }
    // Compute the number of pairs in PK
    pk_size_ = INDEX2(pk_pairs_-1, pk_pairs_-1) + 1;

    size_t memory_needed = pk_size_;
    if(build_k) memory_needed *= 2;

    // Allocate memory for the PJ and PK matrices (using a vector)
    if (memory_needed < 0.8 *memory_) {
        p_j_ = new double[pk_size_];
        ::memset(p_j_, 0 , pk_size_ * sizeof(double));
        if(build_k){
            p_k_ = new double[pk_size_];
            ::memset(p_k_, 0 , pk_size_ * sizeof(double));
        }

        if (p_j_ == NULL || p_k_ == NULL) {
            if(p_j_) delete [] p_j_;
            if(p_k_) delete [] p_k_;
            throw PSIEXCEPTION("Not enough memory for PK matrix, try scf_type = out_of_core instead");
        } else {
            // PK and PJ are zeroed out just before filling, not here
            if(print_ > 2){
                fprintf(outfile,
                "  Allocated %lu elements (%lu pairs) for PJ. (%5f MiB)\n",
                (unsigned long)pk_size_, (unsigned long)pk_pairs_, pk_size_ * 8.0 / 1048576.0);
                fprintf(outfile,
                "  Allocated %lu elements (%lu pairs) for PK. (%5f MiB)\n\n",
                (unsigned long)pk_size_, (unsigned long)pk_pairs_, pk_size_ * 8.0 / 1048576.0);
            }
        }
    } else {
        throw PSIEXCEPTION("Not enough memory for PK matrix, try scf_type = out_of_core instead");
        // TODO just write an out of core code instead. Shouldn't be too hard.
    }

    // struct iwlbuf ERIIN;
    int last_buffer, nbuf;
    int i, j, k, l;
    int ii, jj, kk, ll;
    int is, js, ks, ls;
    int fi;
    size_t bra, ket, braket;
    int idx;
    double value;

    // PK zeroed out during allocation
    if(print_ > 1){
        fprintf(outfile, "  Forming PK and K matrices.\n");
        fflush(outfile);
    }

    IWL iwl(psio_.get(), PSIF_SO_TEI, 0.0, 1, 1);

    do {
        last_buffer = iwl.last_buffer();
        nbuf = iwl.buffer_count();

        fi = 0;
        for (idx = 0; idx < nbuf; ++idx) {
            i = iwl.labels()[fi] > 0 ? iwl.labels()[fi] : -iwl.labels()[fi];
            j = iwl.labels()[fi+1];
            k = iwl.labels()[fi+2];
            l = iwl.labels()[fi+3];
            value = iwl.values()[idx];
            fi += 4;

            // Get the symmetries
            is = so2symblk_[i];
            js = so2symblk_[j];
            ks = so2symblk_[k];
            ls = so2symblk_[l];

            // Get the offset of the SO index in its symblock
            ii = so2index_[i];
            jj = so2index_[j];
            kk = so2index_[k];
            ll = so2index_[l];

            // J
            if ((is == js) && (ks == ls)) {
                bra = INDEX2(ii, jj);
                ket = INDEX2(kk, ll);
                // pk_symoffset_ corrects for the symmetry offset in the pk_ vector
                braket = INDEX2(bra + pk_symoffset[is], ket + pk_symoffset[ks]);
                p_j_[braket] += value;

                // K (2nd sort)
                if (build_k && (ii != jj) && (kk != ll)) {
                    if ((is == ls) && (js == ks)) {
                        bra = INDEX2(ii, ll);
                        ket = INDEX2(jj, kk);
                        braket = INDEX2(bra + pk_symoffset[is], ket + pk_symoffset[js]);
                        if ((ii == ll) || (jj == kk)) {
                            p_k_[braket] += value;
                        } else {
                            p_k_[braket] += 0.5 * value;
                        }
                    }
                }
            }

            // K (1st sort)
            if (build_k && (is == ks) && (js == ls)) {
                bra = INDEX2(ii, kk);
                ket = INDEX2(jj, ll);
                braket = INDEX2(bra + pk_symoffset[is], ket + pk_symoffset[js]);
                if ((ii == kk) || (jj == ll)) {
                    p_k_[braket] += value;
                } else {
                    p_k_[braket] += 0.5 * value;
                }
            }
        }

        if (!last_buffer)
            iwl.fetch();
    } while (!last_buffer);


    // Going out of scope will close the buffer and kill the integrals
    iwl.set_keep_flag(1);

    delete [] pk_symoffset;
    // After stage two is complete, the elements of P must be halved for the case IJ=KL.
    for (size_t ij=0; ij < pk_pairs_; ++ij) {
        size_t address = INDEX2(ij,ij);
        p_j_[address] *= 0.5;
        if(build_k) p_k_[address] *= 0.5;
    }
}

}}// Namespaces

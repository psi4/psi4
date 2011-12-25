#ifndef INTEGRALTRANSFORM_FUNCTORS_H
#define INTEGRALTRANSFORM_FUNCTORS_H

#include "libdpd/dpd.h"
#include "libiwl/iwl.hpp"
#include "libpsio/psio.hpp"

namespace psi{

class FrozenCoreAndFockRestrictedFunctor {
private:
    /// The density matrix (stored as a lower triangular array)
    const double *D_;
    /// The frozen core density matrix (stored as a lower triangular array)
    const double *FzD_;
    /// The Fock matrix (stored as a lower triangular array)
    double *F_;
    /// The frozen core operator (stored as a lower triangular array)
    double *Fz_;
    /// Some temporary arrays for handling permutations
    int al[8], bl[8], cl[8], dl[8], dum, found;
public:
    /*
     * Generates the fock and frozen core operators, given the density matrices as input
     *
     * @param D: The density matrix (stored as a lower triangular array)
     * @param FzD: The frozen core density matrix (stored as a lower triangular array)
     * @param F: The Fock matrix (stored as a lower triangular array)
     * @param Fz: The frozen core operator (stored as a lower triangular array)
     */
    FrozenCoreAndFockRestrictedFunctor(const double *D,
                                       const double *FzD,
                                       double *F,
                                       double *Fz)
        : D_(D), FzD_(FzD), F_(F), Fz_(Fz)
    {}

    void operator()(int pabs, int qabs, int rabs, int sabs,
                    int psym, int prel, int qsym, int qrel,
                    int rsym, int rrel, int ssym, int srel, double value)
    {
        int a = al[0] = pabs;
        int b = bl[0] = qabs;
        int c = cl[0] = rabs;
        int d = dl[0] = sabs;
        int ab = INDEX(a,b);
        int cd = INDEX(c,d);
        int bc = INDEX(b,c);
        int ad = INDEX(a,d);

        Fz_[cd] += 2.0 * FzD_[ab] * value;
        F_[cd]  += 2.0 * D_[ab] * value;
        if(b >= c){
            Fz_[bc] -= FzD_[ad] * value;
            F_[bc]  -= D_[ad] * value;
        }

        a = al[1] = qabs;
        b = bl[1] = pabs;
        c = cl[1] = rabs;
        d = dl[1] = sabs;
        if(!(a==al[0] && b==bl[0] && c==cl[0] && d==dl[0])) {
            ab = INDEX(a,b);
            cd = INDEX(c,d);
            bc = INDEX(b,c);
            ad = INDEX(a,d);
            if(c >= d){
                Fz_[cd] += 2.0 * FzD_[ab] * value;
                F_[cd]  += 2.0 * D_[ab] * value;
            }
            if(b >= c){
                Fz_[bc] -= FzD_[ad] * value;
                F_[bc]  -= D_[ad] * value;
            }
        }

        a = al[2] = pabs;
        b = bl[2] = qabs;
        c = cl[2] = sabs;
        d = dl[2] = rabs;
        for(dum=0,found=0; dum < 2 && !found; dum++)
            if(a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
        if(!found){
            ab = INDEX(a,b);
            cd = INDEX(c,d);
            bc = INDEX(b,c);
            ad = INDEX(a,d);
            if(c >= d){
                Fz_[cd] += 2.0 * FzD_[ab] * value;
                F_[cd]  += 2.0 * D_[ab] * value;
            }
            if(b >= c){
                Fz_[bc] -= FzD_[ad] * value;
                F_[bc]  -= D_[ad] * value;
            }
        }

        a = al[3] = qabs;
        b = bl[3] = pabs;
        c = cl[3] = sabs;
        d = dl[3] = rabs;
        for(dum=0,found=0; dum < 3 && !found; dum++)
            if(a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
        if(!found){
            ab = INDEX(a,b);
            cd = INDEX(c,d);
            bc = INDEX(b,c);
            ad = INDEX(a,d);
            if(c >= d){
                Fz_[cd] += 2.0 * FzD_[ab] * value;
                F_[cd]  += 2.0 * D_[ab] * value;
            }
            if(b >= c){
                Fz_[bc] -= FzD_[ad] * value;
                F_[bc]  -= D_[ad] * value;
            }
        }

        a = al[4] = rabs;
        b = bl[4] = sabs;
        c = cl[4] = pabs;
        d = dl[4] = qabs;
        for(dum=0,found=0; dum < 4 && !found; dum++)
            if(a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
        if(!found){
            ab = INDEX(a,b);
            cd = INDEX(c,d);
            bc = INDEX(b,c);
            ad = INDEX(a,d);
            if(c >= d){
                Fz_[cd] += 2.0 * FzD_[ab] * value;
                F_[cd]  += 2.0 * D_[ab] * value;
            }
            if(b >= c){
                Fz_[bc] -= FzD_[ad] * value;
                F_[bc]  -= D_[ad] * value;
            }
        }

        a = al[5] = rabs;
        b = bl[5] = sabs;
        c = cl[5] = qabs;
        d = dl[5] = pabs;
        for(dum=0,found=0; dum < 5 && !found; dum++)
            if(a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
        if(!found) {
            ab = INDEX(a,b);
            cd = INDEX(c,d);
            bc = INDEX(b,c);
            ad = INDEX(a,d);
            if(c >= d){
                Fz_[cd] += 2.0 * FzD_[ab] * value;
                F_[cd]  += 2.0 * D_[ab] * value;
            }
            if(b >= c){
                Fz_[bc] -= FzD_[ad] * value;
                F_[bc]  -= D_[ad] * value;
            }
        }

        a = al[6] = sabs;
        b = bl[6] = rabs;
        c = cl[6] = pabs;
        d = dl[6] = qabs;
        for(dum=0, found=0; dum < 6 && !found; ++dum)
            if(a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
        if(!found){
            ab = INDEX(a,b);
            cd = INDEX(c,d);
            bc = INDEX(b,c);
            ad = INDEX(a,d);
            if(c >= d){
                Fz_[cd] += 2.0 * FzD_[ab] * value;
                F_[cd]  += 2.0 * D_[ab] * value;
            }
            if(b >= c){
                Fz_[bc] -= FzD_[ad] * value;
                F_[bc]  -= D_[ad] * value;
            }
        }

        a = al[7] = sabs;
        b = bl[7] = rabs;
        c = cl[7] = qabs;
        d = dl[7] = pabs;
        for(dum=0,found=0; dum < 7 && !found; dum++)
            if(a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
        if(!found){
            ab = INDEX(a,b);
            cd = INDEX(c,d);
            bc = INDEX(b,c);
            ad = INDEX(a,d);
            if(c >= d){
                Fz_[cd] += 2.0 * FzD_[ab] * value;
                F_[cd]  += 2.0 * D_[ab] * value;
            }
            if(b >= c){
                Fz_[bc] -= FzD_[ad] * value;
                F_[bc]  -= D_[ad] * value;
            }
        }
    }
};

class FrozenCoreAndFockUnrestrictedFunctor {
private:
    /// The alpha density matrix (stored as a lower triangular array)
    const double *Da_;
    /// The beta density matrix (stored as a lower triangular array)
    const double *Db_;
    /// The alpha frozen core density matrix (stored as a lower triangular array)
    const double *FzDa_;
    /// The beta frozen core density matrix (stored as a lower triangular array)
    const double *FzDb_;
    /// The alpha Fock matrix (stored as a lower triangular array)
    double *Fa_;
    /// The beta Fock matrix (stored as a lower triangular array)
    double *Fb_;
    /// The alpha frozen core operator (stored as a lower triangular array)
    double *Fza_;
    /// The beta frozen core operator (stored as a lower triangular array)
    double *Fzb_;
    /// Some temporary arrays for handling permutations
    int al[8], bl[8], cl[8], dl[8], dum, found;
public:
    /*
     * Generates the fock and frozen core operators, given the density matrices as input
     *
     * @param Da: The alpha density matrix (stored as a lower triangular array)
     * @param Db: The beta density matrix (stored as a lower triangular array)
     * @param FzDa: The alpha frozen core density matrix (stored as a lower triangular array)
     * @param FzDb: The beta frozen core density matrix (stored as a lower triangular array)
     * @param Fa: The alpha Fock matrix (stored as a lower triangular array)
     * @param Fa: The beta Fock matrix (stored as a lower triangular array)
     * @param Fza: The alpha frozen core operator (stored as a lower triangular array)
     * @param Fzb: The beta frozen core operator (stored as a lower triangular array)
     */
    FrozenCoreAndFockUnrestrictedFunctor(const double *Da, const double *Db,
                                              const double *FzDa, const double *FzDb,
                                              double *Fa, double *Fb,
                                              double *Fza, double *Fzb)
        : Da_(Da), Db_(Db), FzDa_(FzDa), FzDb_(FzDb), Fa_(Fa), Fb_(Fb), Fza_(Fza), Fzb_(Fzb)
    {}

    void operator()(int pabs, int qabs, int rabs, int sabs,
                    int psym, int prel, int qsym, int qrel,
                    int rsym, int rrel, int ssym, int srel, double value)
    {
        int a = al[0] = pabs;
        int b = bl[0] = qabs;
        int c = cl[0] = rabs;
        int d = dl[0] = sabs;
        int ab = INDEX(a,b);
        int cd = INDEX(c,d);
        int bc = INDEX(b,c);
        int ad = INDEX(a,d);
        Fza_[cd] += (FzDa_[ab] + FzDb_[ab]) * value;
        Fzb_[cd] += (FzDa_[ab] + FzDb_[ab]) * value;
        Fa_[cd]  += (Da_[ab] + Db_[ab]) * value;
        Fb_[cd]  += (Da_[ab] + Db_[ab]) * value;
        if(b >= c) {
            Fza_[bc] -= FzDa_[ad] * value;
            Fa_[bc]  -= Da_[ad] * value;
            Fzb_[bc] -= FzDb_[ad] * value;
            Fb_[bc]  -= Db_[ad] * value;
        }

        a = al[1] = qabs;
        b = bl[1] = pabs;
        c = cl[1] = rabs;
        d = dl[1] = sabs;
        if(!(a==al[0] && b==bl[0] && c==cl[0] && d==dl[0])){
            ab = INDEX(a,b);
            cd = INDEX(c,d);
            bc = INDEX(b,c);
            ad = INDEX(a,d);
            if(c >= d){
                Fza_[cd] += (FzDa_[ab] + FzDb_[ab]) * value;
                Fa_[cd] += (Da_[ab] + Db_[ab]) * value;
                Fzb_[cd] += (FzDa_[ab] + FzDb_[ab]) * value;
                Fb_[cd] += (Da_[ab] + Db_[ab]) * value;
            }
            if(b >= c){
                Fza_[bc] -= FzDa_[ad] * value;
                Fa_[bc]  -= Da_[ad] * value;
                Fzb_[bc] -= FzDb_[ad] * value;
                Fb_[bc]  -= Db_[ad] * value;
            }
        }

        a = al[2] = pabs;
        b = bl[2] = qabs;
        c = cl[2] = sabs;
        d = dl[2] = rabs;
        for(dum=0,found=0; dum < 2 && !found; dum++)
            if(a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
        if(!found){
            ab = INDEX(a,b);
            cd = INDEX(c,d);
            bc = INDEX(b,c);
            ad = INDEX(a,d);
            if(c >= d){
                Fza_[cd] += (FzDa_[ab] + FzDb_[ab]) * value;
                Fa_[cd] += (Da_[ab] + Db_[ab]) * value;
                Fzb_[cd] += (FzDa_[ab] + FzDb_[ab]) * value;
                Fb_[cd] += (Da_[ab] + Db_[ab]) * value;
            }
            if(b >= c){
                Fza_[bc] -= FzDa_[ad] * value;
                Fa_[bc]  -= Da_[ad] * value;
                Fzb_[bc] -= FzDb_[ad] * value;
                Fb_[bc]  -= Db_[ad] * value;
            }
        }

        a = al[3] = qabs;
        b = bl[3] = pabs;
        c = cl[3] = sabs;
        d = dl[3] = rabs;
        for(dum=0,found=0; dum < 3 && !found; dum++)
            if(a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
        if(!found){
            ab = INDEX(a,b);
            cd = INDEX(c,d);
            bc = INDEX(b,c);
            ad = INDEX(a,d);
            if(c >= d){
                Fza_[cd] += (FzDa_[ab] + FzDb_[ab]) * value;
                Fa_[cd]  += (Da_[ab] + Db_[ab]) * value;
                Fzb_[cd] += (FzDa_[ab] + FzDb_[ab]) * value;
                Fb_[cd]  += (Da_[ab] + Db_[ab]) * value;
            }
            if(b >= c){
                Fza_[bc] -= FzDa_[ad] * value;
                Fa_[bc]  -= Da_[ad] * value;
                Fzb_[bc] -= FzDb_[ad] * value;
                Fb_[bc]  -= Db_[ad] * value;
            }
        }

        a = al[4] = rabs;
        b = bl[4] = sabs;
        c = cl[4] = pabs;
        d = dl[4] = qabs;
        for(dum=0,found=0; dum < 4 && !found; dum++)
            if(a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
        if(!found){
            ab = INDEX(a,b);
            cd = INDEX(c,d);
            bc = INDEX(b,c);
            ad = INDEX(a,d);
            if(c >= d){
                Fza_[cd] += (FzDa_[ab] + FzDb_[ab]) * value;
                Fa_[cd]  += (Da_[ab] + Db_[ab]) * value;
                Fzb_[cd] += (FzDa_[ab] + FzDb_[ab]) * value;
                Fb_[cd]  += (Da_[ab] + Db_[ab]) * value;
            }
            if(b >= c){
                Fza_[bc] -= FzDa_[ad] * value;
                Fa_[bc]  -= Da_[ad] * value;
                Fzb_[bc] -= FzDb_[ad] * value;
                Fb_[bc]  -= Db_[ad] * value;
            }
        }

        a = al[5] = rabs;
        b = bl[5] = sabs;
        c = cl[5] = qabs;
        d = dl[5] = pabs;
        for(dum=0,found=0; dum < 5 && !found; dum++)
            if(a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
        if(!found){
            ab = INDEX(a,b);
            cd = INDEX(c,d);
            bc = INDEX(b,c);
            ad = INDEX(a,d);
            if(c >= d){
                Fza_[cd] += (FzDa_[ab] + FzDb_[ab]) * value;
                Fa_[cd]  += (Da_[ab] + Db_[ab]) * value;
                Fzb_[cd] += (FzDa_[ab] + FzDb_[ab]) * value;
                Fb_[cd]  += (Da_[ab] + Db_[ab]) * value;
            }
            if(b >= c){
                Fza_[bc] -= FzDa_[ad] * value;
                Fa_[bc]  -= Da_[ad] * value;
                Fzb_[bc] -= FzDb_[ad] * value;
                Fb_[bc]  -= Db_[ad] * value;
            }
        }

        a = al[6] = sabs;
        b = bl[6] = rabs;
        c = cl[6] = pabs;
        d = dl[6] = qabs;
        for(dum=0,found=0; dum < 6 && !found; dum++)
            if(a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
        if(!found) {
            ab = INDEX(a,b);
            cd = INDEX(c,d);
            bc = INDEX(b,c);
            ad = INDEX(a,d);
            if(c >= d){
                Fza_[cd] += (FzDa_[ab] + FzDb_[ab]) * value;
                Fa_[cd]  += (Da_[ab] + Db_[ab]) * value;
                Fzb_[cd] += (FzDa_[ab] + FzDb_[ab]) * value;
                Fb_[cd]  += (Da_[ab] + Db_[ab]) * value;
            }
            if(b >= c){
                Fza_[bc] -= FzDa_[ad] * value;
                Fa_[bc]  -= Da_[ad] * value;
                Fzb_[bc] -= FzDb_[ad] * value;
                Fb_[bc]  -= Db_[ad] * value;
            }
        }

        a = al[7] = sabs;
        b = bl[7] = rabs;
        c = cl[7] = qabs;
        d = dl[7] = pabs;
        for(dum=0,found=0; dum < 7 && !found; dum++)
            if(a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
        if(!found){
            ab = INDEX(a,b);
            cd = INDEX(c,d);
            bc = INDEX(b,c);
            ad = INDEX(a,d);
            if(c >= d){
                Fza_[cd] += (FzDa_[ab] + FzDb_[ab]) * value;
                Fa_[cd]  += (Da_[ab] + Db_[ab]) * value;
                Fzb_[cd] += (FzDa_[ab] + FzDb_[ab]) * value;
                Fb_[cd]  += (Da_[ab] + Db_[ab]) * value;
            }
            if(b >= c){
                Fza_[bc] -= FzDa_[ad] * value;
                Fa_[bc]  -= Da_[ad] * value;
                Fzb_[bc] -= FzDb_[ad] * value;
                Fb_[bc]  -= Db_[ad] * value;
            }
        }
    }
};

class DPDFillerFunctor {
private:
    dpdfile4 *file_;
    dpdparams4 *params_;
    int this_bucket_;
    int **bucket_map_;
    int **bucket_offset_;
    bool symmetrize_;
    bool have_bra_ket_sym_;
public:
    DPDFillerFunctor(dpdfile4 *file, int this_bucket, int **bucket_map,
                     int **bucket_offset, bool symmetrize, bool have_bra_ket_sym):
        file_(file), this_bucket_(this_bucket), bucket_map_(bucket_map),
        bucket_offset_(bucket_offset), symmetrize_(symmetrize),
        have_bra_ket_sym_(have_bra_ket_sym)
    {
        params_ = file_->params;
    }
    void operator()(int p, int q, int r, int s, double value)
    {
        if(symmetrize_){
            // Symmetrize the quantity (used in density matrix processing)
            if(p!=q) value *= 0.5;
            if(r!=s) value *= 0.5;
        }

        bool bra_ket_different = !(p==r && q==s);

        /* Get the orbital symmetries */
        int p_sym = params_->psym[p];
        int q_sym = params_->qsym[q];
        int r_sym = params_->rsym[r];
        int s_sym = params_->ssym[s];
        int pq_sym = p_sym^q_sym;
        int rs_sym = r_sym^s_sym;

        /* The allowed (Mulliken) permutations are very simple in this case */
        if(bucket_map_[p][q] == this_bucket_) {
            /* Get the row and column indices and assign the value */
            int pq = params_->rowidx[p][q];
            int rs = params_->colidx[r][s];
            int offset = bucket_offset_[this_bucket_][pq_sym];
            if((pq-offset >= params_->rowtot[pq_sym]) || (rs >= params_->coltot[rs_sym]))
                error("MP Params_make: pq, rs", p,q,r,s,pq,rs,pq_sym,rs_sym);
            file_->matrix[pq_sym][pq-offset][rs] += value;
        }

        /*
         * We also add in the bra-ket transposed value, as a result of the matrix
         * storage, but we need to make sure we don't duplicate "diagonal" values.
         * We don't do this if the quantity does not have bra-ket symmetry, like
         * in the Alpha-Beta TPDM.
         */
        if(bucket_map_[r][s] == this_bucket_ && bra_ket_different && have_bra_ket_sym_) {
            int rs = params_->rowidx[r][s];
            int pq = params_->colidx[p][q];
            int offset = bucket_offset_[this_bucket_][rs_sym];
            if((rs-offset >= params_->rowtot[rs_sym])||(pq >= params_->coltot[pq_sym]))
                error("MP Params_make: rs, pq", p,q,r,s,rs,pq,rs_sym,pq_sym);
            file_->matrix[rs_sym][rs-offset][pq] += value;
        }
    }

private:
    void error(const char *message, int p, int q, int r, int s,
               int pq, int rs, int pq_sym, int rs_sym)
    {

        fprintf(outfile, "\n\tDPD Parameter Error in %s\n", message);
        fprintf(outfile,"\t-------------------------------------------------\n");
        fprintf(outfile,"\t    p      q      r      s  [   pq]  [   rs] pq_symm rs_symm\n");
        fprintf(outfile,"\t%5d  %5d  %5d  %5d  [%5d]  [%5d]   %1d   %1d\n", p,q,r,s,
                pq,rs,pq_sym,rs_sym);
        throw PsiException("DPD idx failure.", __FILE__, __LINE__);
    }
};

class NullFunctor {
public:
    /*
     * Just an empty functor that will be compiled away (hopefully)
     */
    void operator()(int pabs, int qabs, int rabs, int sabs,
                    int psym, int prel, int qsym, int qrel,
                    int rsym, int rrel, int ssym, int srel, double value)
    {
    }
    void operator()(int p, int q, int r, int s, double value)
    {
    }
};

template <class DPDFunctor, class FockFunctor>
void iwl_integrals(IWL* iwl, DPDFunctor &dpd, FockFunctor &fock)
{
    Label *lblptr = iwl->labels();
    Value *valptr = iwl->values();
    int labelIndex, p, q, r, s;
    double value;
    bool lastBuffer;
    do{
        lastBuffer = iwl->last_buffer();
        for(int index = 0; index < iwl->buffer_count(); ++index){
            labelIndex = 4*index;
            p = abs((int) lblptr[labelIndex++]);
            q = (int) lblptr[labelIndex++];
            r = (int) lblptr[labelIndex++];
            s = (int) lblptr[labelIndex++];
            value = (double) valptr[index];
            dpd(p,q,r,s,value);
            fock(p,q,r,s,0,0,0,0,0,0,0,0,value);
        } /* end loop through current buffer */
        if(!lastBuffer) iwl->fetch();
    }while(!lastBuffer);
    iwl->set_keep_flag(1);
}

} // Namespaces
#endif // INTEGRALTRANSFORM_FUNCTORS_H

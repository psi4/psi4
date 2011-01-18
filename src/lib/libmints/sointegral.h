#ifndef _psi_src_lib_libmints_sointegral_h_
#define _psi_src_lib_libmints_sointegral_h_

#include "onebody.h"
#include "twobody.h"
#include "basisset.h"
#include "integral.h"
#include "sobasis.h"
#include "gshell.h"
#include "petitelist.h"

namespace boost {
template <class T>
class shared_ptr;
}

namespace psi {

class Matrix;

class OneBodySOInt
{
protected:
    boost::shared_ptr<OneBodyInt> ob_;
    boost::shared_ptr<IntegralFactory> integral_;

    boost::shared_ptr<SOBasis> b1_;
    boost::shared_ptr<SOBasis> b2_;

    double *buffer_;

    int only_totally_symmetric_;
public:
    OneBodySOInt(const boost::shared_ptr<OneBodyInt>& , const boost::shared_ptr<IntegralFactory> &);
    virtual ~OneBodySOInt();

    boost::shared_ptr<SOBasis> basis() const;
    boost::shared_ptr<SOBasis> basis1() const;
    boost::shared_ptr<SOBasis> basis2() const;

    /**
     * Resulting integral buffer.
     */
    const double* buffer() const { return buffer_; }

    /**
     * Computes a one-electron integral matrix. Only works for symmetric operators
     * (multipole operators will not work).
     */
    void compute(boost::shared_ptr<Matrix> result);

    /**
     * Computes a one-electron integral matrix. Only works for symmetric operators
     * (multipole operators will not work).
     */
    void compute_petitelist(boost::shared_ptr<Matrix> result);

    /**
     * Computes a one-electron integral matrix. Only works for symmetric operators
     * (multipole operators will not work).
     */
    void compute_dcr(boost::shared_ptr<Matrix> result);

    /**
     * Compute a given SO shell doublet placing the result into buffer_
     */
    virtual void compute_shell(int, int);

    /**
     * Currently these settings are ignored, placed in for future work.
     */
    int only_totally_symmetric() const { return only_totally_symmetric_; }
    void set_only_totally_symmetric(int i) { only_totally_symmetric_ = i; }
};

// Only include the following function if Doxygen is running to generate appropriate
// documentation.
#ifdef DOXYGEN
class TwoBodySOIntFunctor
{
public:
    void operator()(int pirrep, int pso, int qirrep, int qso, int rirrep, int rso, int sirrep, int sso, double value);
};
#endif

class TwoBodySOInt
{
protected:
    boost::shared_ptr<TwoBodyInt> tb_;
    boost::shared_ptr<IntegralFactory> integral_;

    boost::shared_ptr<SOBasis> b1_;
    boost::shared_ptr<SOBasis> b2_;
    boost::shared_ptr<SOBasis> b3_;
    boost::shared_ptr<SOBasis> b4_;

    double *buffer_;

public:
    TwoBodySOInt(const boost::shared_ptr<TwoBodyInt>&,
                 const boost::shared_ptr<IntegralFactory>&);
    virtual ~TwoBodySOInt();

    boost::shared_ptr<SOBasis> basis() const;
    boost::shared_ptr<SOBasis> basis1() const;
    boost::shared_ptr<SOBasis> basis2() const;
    boost::shared_ptr<SOBasis> basis3() const;
    boost::shared_ptr<SOBasis> basis4() const;

    const double *buffer() const { return buffer_; }

    template<typename TwoBodySOIntFunctor>
    void compute_shell(int, int, int, int, TwoBodySOIntFunctor& body);

    template<typename TwoBodySOIntFunctor>
    void compute_shell_smart(int, int, int, int, TwoBodySOIntFunctor& body);
};

template<typename TwoBodySOFunctor>
void TwoBodySOInt::compute_shell(int ish, int jsh, int ksh, int lsh, TwoBodySOFunctor& body)
{
//    fprintf(outfile, "computing shell (%d %d %d %d)\n", ish, jsh, ksh, lsh);

    const double *aobuff = tb_->buffer();

    const SOTransform &t1 = b1_->trans(ish);
    const SOTransform &t2 = b2_->trans(jsh);
    const SOTransform &t3 = b3_->trans(ksh);
    const SOTransform &t4 = b4_->trans(lsh);

    int nso1 = b1_->nfunction(ish);
    int nso2 = b2_->nfunction(jsh);
    int nso3 = b3_->nfunction(ksh);
    int nso4 = b4_->nfunction(lsh);

    int nao1 = b1_->naofunction(ish);
    int nao2 = b2_->naofunction(jsh);
    int nao3 = b3_->naofunction(ksh);
    int nao4 = b4_->naofunction(lsh);

//    fprintf(outfile, "nao1 = %d nao2 = %d nao3 = %d nao4 = %d\n", nao1, nao2, nao3, nao4);

    memset(buffer_, 0, 16*nao1*nao2*nao3*nao4*sizeof(double));

    // loop through the ao shells that make up this so shell
    for (int i=0; i<t1.naoshell; i++) {
        const SOTransformShell &s1 = t1.aoshell[i];
        for (int j=0; j<t2.naoshell; j++) {
            const SOTransformShell &s2 = t2.aoshell[j];
            for (int k=0; k<t3.naoshell; k++) {
                const SOTransformShell &s3 = t3.aoshell[k];
                for (int l=0; l<t4.naoshell; l++) {
                    const SOTransformShell &s4 = t4.aoshell[l];
                    tb_->compute_shell(s1.aoshell, s2.aoshell, s3.aoshell, s4.aoshell);

//                    fprintf(outfile, "ao: (%d %d %d %d)\n", s1.aoshell, s2.aoshell,
//                            s3.aoshell, s4.aoshell);

//                    for (int z=0; z < INT_NPURE(tb_->basis1()->shell(s1.aoshell)->am()) *
//                                      INT_NPURE(tb_->basis2()->shell(s2.aoshell)->am()) *
//                                      INT_NPURE(tb_->basis3()->shell(s3.aoshell)->am()) *
//                                      INT_NPURE(tb_->basis4()->shell(s4.aoshell)->am()); ++z) {
//                        fprintf(outfile, "raw: %d -> %8.5f\n", z, aobuff[z]);
//                    }

                    for (int itr=0; itr<s1.nfunc; itr++) {
                        const SOTransformFunction &ifunc = s1.func[itr];
                        double icoef = ifunc.coef;
                        int iaofunc = ifunc.aofunc;
                        int isofunc = b1_->function_offset_within_shell(ish,
                                                                        ifunc.irrep)
                                + ifunc.sofunc;
                        int iaooff = iaofunc;
                        int isooff = isofunc;
                        for (int jtr=0; jtr<s2.nfunc; jtr++) {
                            const SOTransformFunction &jfunc = s2.func[jtr];
                            double jcoef = jfunc.coef * icoef;
                            int jaofunc = jfunc.aofunc;
                            int jsofunc = b2_->function_offset_within_shell(jsh,
                                                                            jfunc.irrep)
                                    + jfunc.sofunc;
                            int jaooff = iaooff*nao2 + jaofunc;
                            int jsooff = isooff*nso2 + jsofunc;
                            for (int ktr=0; ktr<s3.nfunc; ktr++) {
                                const SOTransformFunction &kfunc = s3.func[ktr];
                                double kcoef = kfunc.coef * jcoef;
                                int kaofunc = kfunc.aofunc;
                                int ksofunc = b3_->function_offset_within_shell(ksh,
                                                                                kfunc.irrep)
                                        + kfunc.sofunc;
                                int kaooff = jaooff*nao3 + kaofunc;
                                int ksooff = jsooff*nso3 + ksofunc;
                                for (int ltr=0; ltr<s4.nfunc; ltr++) {
                                    const SOTransformFunction &lfunc = s4.func[ltr];
                                    double lcoef = lfunc.coef * kcoef;
                                    int laofunc = lfunc.aofunc;
                                    int lsofunc = b4_->function_offset_within_shell(lsh,
                                                                                    lfunc.irrep)
                                            + lfunc.sofunc;
                                    int laooff = kaooff*nao4 + laofunc;
                                    int lsooff = ksooff*nso4 + lsofunc;
                                    buffer_[lsooff] += lcoef * aobuff[laooff];
//                                    if (fabs(aobuff[laooff]*lcoef) > 1.0e-10) {
//                                        fprintf(outfile, "(%d %d|%d %d) += %8.5lf * (%d %d|%d %d): %8.5lf (laoff = %d) -> %8.5lf (lsoff = %d)\n",
//                                                isofunc, jsofunc, ksofunc, lsofunc, lcoef,
//                                                iaofunc, jaofunc, kaofunc, laofunc,
//                                                aobuff[laooff], laooff,
//                                                buffer_[lsooff], lsooff);
//                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    const SOTransformShell &s1 = t1.aoshell[0];
    const SOTransformShell &s2 = t2.aoshell[0];
    const SOTransformShell &s3 = t3.aoshell[0];
    const SOTransformShell &s4 = t4.aoshell[0];

    for (int itr=0; itr<s1.nfunc; itr++) {
        const SOTransformFunction &ifunc = s1.func[itr];
        double icoef = ifunc.coef;
        int iaofunc = ifunc.aofunc;
        int isofunc = b1_->function_offset_within_shell(ish,
                                                        ifunc.irrep)
                + ifunc.sofunc;
        int iaooff = iaofunc;
        int isooff = isofunc;
        for (int jtr=0; jtr<s2.nfunc; jtr++) {
            const SOTransformFunction &jfunc = s2.func[jtr];
            double jcoef = jfunc.coef * icoef;
            int jaofunc = jfunc.aofunc;
            int jsofunc = b2_->function_offset_within_shell(jsh,
                                                            jfunc.irrep)
                    + jfunc.sofunc;
            int jaooff = iaooff*nao2 + jaofunc;
            int jsooff = isooff*nso2 + jsofunc;
            for (int ktr=0; ktr<s3.nfunc; ktr++) {
                const SOTransformFunction &kfunc = s3.func[ktr];
                double kcoef = kfunc.coef * jcoef;
                int kaofunc = kfunc.aofunc;
                int ksofunc = b3_->function_offset_within_shell(ksh,
                                                                kfunc.irrep)
                        + kfunc.sofunc;
                int kaooff = jaooff*nao3 + kaofunc;
                int ksooff = jsooff*nso3 + ksofunc;
                for (int ltr=0; ltr<s4.nfunc; ltr++) {
                    const SOTransformFunction &lfunc = s4.func[ltr];
                    double lcoef = lfunc.coef * kcoef;
                    int laofunc = lfunc.aofunc;
                    int lsofunc = b4_->function_offset_within_shell(lsh,
                                                                    lfunc.irrep)
                            + lfunc.sofunc;
                    int laooff = kaooff*nao4 + laofunc;
                    int lsooff = ksooff*nso4 + lsofunc;

//                    buffer_[lsooff] += lcoef * aobuff[laooff];
                    if (fabs(buffer_[lsooff]) > 1.0e-16)
                        body(ifunc.irrep, isofunc,
                             jfunc.irrep, jsofunc,
                             kfunc.irrep, ksofunc,
                             lfunc.irrep, lsofunc,
                             buffer_[lsooff]);
                }
            }
        }
    }
}

template<typename TwoBodySOFunctor>
void TwoBodySOInt::compute_shell_smart(int ish, int jsh, int ksh, int lsh, TwoBodySOFunctor& body)
{
    // ish ... are unique shells.
    fprintf(outfile, "\tcomputing SO shell quartet (%d %d %d %d)\n", ish, jsh, ksh, lsh);

    const double *aobuff = tb_->buffer();
    boost::shared_ptr<PetiteList> petite = b1_->petitelist();

    const SOTransform &t1 = b1_->trans(ish);
    const SOTransform &t2 = b2_->trans(jsh);
    const SOTransform &t3 = b3_->trans(ksh);
    const SOTransform &t4 = b4_->trans(lsh);

    int nso1 = b1_->nfunction(ish);
    int nso2 = b2_->nfunction(jsh);
    int nso3 = b3_->nfunction(ksh);
    int nso4 = b4_->nfunction(lsh);

    int nao1 = b1_->naofunction(ish);
    int nao2 = b2_->naofunction(jsh);
    int nao3 = b3_->naofunction(ksh);
    int nao4 = b4_->naofunction(lsh);

//    fprintf(outfile, "nao1 = %d nao2 = %d nao3 = %d nao4 = %d\n", nao1, nao2, nao3, nao4);

    memset(buffer_, 0, 16*nao1*nao2*nao3*nao4*sizeof(double));

    // Map unique shell to shell
    int si = t1.aoshell[t1.naoshell-1].aoshell;
    int sj = t2.aoshell[t2.naoshell-1].aoshell;
    int sk = t3.aoshell[t3.naoshell-1].aoshell;
    int sl = t4.aoshell[t4.naoshell-1].aoshell;

    // Order of the group.
    double order = (double)petite->order();

    fprintf(outfile, "si = %d sj = %d sk = %d sl = %d\n", si, sj, sk, sl);

    for (int i=0; i<t1.naoshell; ++i) {
        const SOTransformShell& s1 = t1.aoshell[i];

        if (!petite->in_p1(s1.aoshell)) {
            fprintf(outfile, "skipping i = %d\n", i);
            continue;
        }

        int atom_i = integral_->basis1()->shell(s1.aoshell)->ncenter();
        int stab_i = petite->stabilizer(atom_i);

        for (int j=0; j<t2.naoshell; ++j) {
            const SOTransformShell& s2 = t2.aoshell[j];

            if (!petite->in_p2(s1.aoshell, s2.aoshell))
                continue;

            int atom_j = integral_->basis2()->shell(s2.aoshell)->ncenter();
            int stab_j = petite->stabilizer(atom_j);

            int stab_ij = petite->GnG(stab_i, stab_j);

            for (int k=0; k<t3.naoshell; ++k) {
                const SOTransformShell& s3 = t3.aoshell[k];

                int atom_k = integral_->basis3()->shell(s3.aoshell)->ncenter();
                int stab_k = petite->stabilizer(atom_k);

                for (int l=0; l<t4.naoshell; ++l) {
                    const SOTransformShell& s4 = t4.aoshell[l];

                    if (!petite->in_p4(s1.aoshell, s2.aoshell, s3.aoshell, s4.aoshell))
                        continue;

                    int atom_l = integral_->basis4()->shell(s4.aoshell)->ncenter();
                    int stab_l = petite->stabilizer(atom_l);

                    int stab_kl = petite->GnG(stab_k, stab_l);

                    double lambda_T = order / (double)petite->dcr_degeneracy(stab_ij, stab_kl);

                    fprintf(outfile, "i = %d j = %d k = %d l = %d lambda_T = %8.5f\n",
                            i, j, k, l, lambda_T);

                    tb_->compute_shell(s1.aoshell, s2.aoshell, s3.aoshell, s4.aoshell);

                    for (int itr=0; itr<s1.nfunc; itr++) {
                        const SOTransformFunction &ifunc = s1.func[itr];
                        double icoef = ifunc.coef;
                        int iaofunc = ifunc.aofunc;
                        int isofunc = b1_->function_offset_within_shell(ish,
                                                                        ifunc.irrep)
                                + ifunc.sofunc;
                        int iaooff = iaofunc;
                        int isooff = isofunc;
                        for (int jtr=0; jtr<s2.nfunc; jtr++) {
                            const SOTransformFunction &jfunc = s2.func[jtr];
                            double jcoef = jfunc.coef * icoef;
                            int jaofunc = jfunc.aofunc;
                            int jsofunc = b2_->function_offset_within_shell(jsh,
                                                                            jfunc.irrep)
                                    + jfunc.sofunc;
                            int jaooff = iaooff*nao2 + jaofunc;
                            int jsooff = isooff*nso2 + jsofunc;
                            for (int ktr=0; ktr<s3.nfunc; ktr++) {
                                const SOTransformFunction &kfunc = s3.func[ktr];
                                double kcoef = kfunc.coef * jcoef;
                                int kaofunc = kfunc.aofunc;
                                int ksofunc = b3_->function_offset_within_shell(ksh,
                                                                                kfunc.irrep)
                                        + kfunc.sofunc;
                                int kaooff = jaooff*nao3 + kaofunc;
                                int ksooff = jsooff*nso3 + ksofunc;
                                for (int ltr=0; ltr<s4.nfunc; ltr++) {
                                    const SOTransformFunction &lfunc = s4.func[ltr];
                                    double lcoef = lfunc.coef * kcoef;
                                    int laofunc = lfunc.aofunc;
                                    int lsofunc = b4_->function_offset_within_shell(lsh,
                                                                                    lfunc.irrep)
                                            + lfunc.sofunc;
                                    int laooff = kaooff*nao4 + laofunc;
                                    int lsooff = ksooff*nso4 + lsofunc;

                                    buffer_[lsooff] += lambda_T * lcoef * aobuff[laooff];
//                                    if (fabs(aobuff[laooff]*lcoef) > 1.0e-10) {
//                                        fprintf(outfile, "(%d %d|%d %d) += %8.5lf * (%d %d|%d %d): %8.5lf (laoff = %d) -> %8.5lf (lsoff = %d)\n",
//                                                isofunc, jsofunc, ksofunc, lsofunc, lcoef,
//                                                iaofunc, jaofunc, kaofunc, laofunc,
//                                                aobuff[laooff], laooff,
//                                                buffer_[lsooff], lsooff);
//                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    const SOTransformShell &s1 = t1.aoshell[0];
    const SOTransformShell &s2 = t2.aoshell[0];
    const SOTransformShell &s3 = t3.aoshell[0];
    const SOTransformShell &s4 = t4.aoshell[0];

    for (int itr=0; itr<s1.nfunc; itr++) {
        const SOTransformFunction &ifunc = s1.func[itr];
        double icoef = ifunc.coef;
        int iaofunc = ifunc.aofunc;
        int isofunc = b1_->function_offset_within_shell(ish,
                                                        ifunc.irrep)
                + ifunc.sofunc;
        int iaooff = iaofunc;
        int isooff = isofunc;
        for (int jtr=0; jtr<s2.nfunc; jtr++) {
            const SOTransformFunction &jfunc = s2.func[jtr];
            double jcoef = jfunc.coef * icoef;
            int jaofunc = jfunc.aofunc;
            int jsofunc = b2_->function_offset_within_shell(jsh,
                                                            jfunc.irrep)
                    + jfunc.sofunc;
            int jaooff = iaooff*nao2 + jaofunc;
            int jsooff = isooff*nso2 + jsofunc;
            for (int ktr=0; ktr<s3.nfunc; ktr++) {
                const SOTransformFunction &kfunc = s3.func[ktr];
                double kcoef = kfunc.coef * jcoef;
                int kaofunc = kfunc.aofunc;
                int ksofunc = b3_->function_offset_within_shell(ksh,
                                                                kfunc.irrep)
                        + kfunc.sofunc;
                int kaooff = jaooff*nao3 + kaofunc;
                int ksooff = jsooff*nso3 + ksofunc;
                for (int ltr=0; ltr<s4.nfunc; ltr++) {
                    const SOTransformFunction &lfunc = s4.func[ltr];
                    double lcoef = lfunc.coef * kcoef;
                    int laofunc = lfunc.aofunc;
                    int lsofunc = b4_->function_offset_within_shell(lsh,
                                                                    lfunc.irrep)
                            + lfunc.sofunc;
                    int laooff = kaooff*nao4 + laofunc;
                    int lsooff = ksooff*nso4 + lsofunc;

//                    buffer_[lsooff] += lcoef * aobuff[laooff];
                    if (fabs(buffer_[lsooff]) > 1.0e-16)
                        body(ifunc.irrep, isofunc,
                             jfunc.irrep, jsofunc,
                             kfunc.irrep, ksofunc,
                             lfunc.irrep, lsofunc,
                             buffer_[lsooff]);
                }
            }
        }
    }
}

}

#endif // _psi_src_lib_libmints_sointegral_h_

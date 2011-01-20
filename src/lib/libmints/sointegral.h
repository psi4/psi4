#ifndef _psi_src_lib_libmints_sointegral_h_
#define _psi_src_lib_libmints_sointegral_h_

#include "onebody.h"
#include "twobody.h"
#include "basisset.h"
#include "integral.h"
#include "sobasis.h"
#include "gshell.h"
#include "petitelist.h"
#include "wavefunction.h"

namespace boost {
template <class T>
class shared_ptr;
}

namespace psi {

class Matrix;

template<typename T>
void swap_index(T& a, T& b) {
    T temp;
    temp = b;
    b = a;
    a = temp;
}

#define SWAP_INDEX(a, b) swap_index(a ## abs, b ## abs); swap_index(a ## rel, b ## rel); swap_index(a ## irrep, b ## irrep);

class OneBodySOInt
{
protected:
    boost::shared_ptr<OneBodyAOInt> ob_;
    const IntegralFactory* integral_;

    boost::shared_ptr<SOBasisSet> b1_;
    boost::shared_ptr<SOBasisSet> b2_;

    double *buffer_;

    int only_totally_symmetric_;
public:
    OneBodySOInt(const boost::shared_ptr<OneBodyAOInt>& , const boost::shared_ptr<IntegralFactory> &);
    OneBodySOInt(const boost::shared_ptr<OneBodyAOInt>& , const IntegralFactory*);
    virtual ~OneBodySOInt();

    boost::shared_ptr<SOBasisSet> basis() const;
    boost::shared_ptr<SOBasisSet> basis1() const;
    boost::shared_ptr<SOBasisSet> basis2() const;

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
    boost::shared_ptr<TwoBodyAOInt> tb_;
    boost::shared_ptr<IntegralFactory> integral_;

    boost::shared_ptr<SOBasisSet> b1_;
    boost::shared_ptr<SOBasisSet> b2_;
    boost::shared_ptr<SOBasisSet> b3_;
    boost::shared_ptr<SOBasisSet> b4_;

    double *buffer_;

    template<typename TwoBodySOIntFunctor>
    void provide_IIII(int, int, int, int, TwoBodySOIntFunctor& body);

    template<typename TwoBodySOIntFunctor>
    void provide_IIJJ(int, int, int, int, TwoBodySOIntFunctor& body);

    template<typename TwoBodySOIntFunctor>
    void provide_IJKL(int, int, int, int, TwoBodySOIntFunctor& body);

public:
    TwoBodySOInt(const boost::shared_ptr<TwoBodyAOInt>&,
                 const boost::shared_ptr<IntegralFactory>&);
    virtual ~TwoBodySOInt();

    boost::shared_ptr<SOBasisSet> basis() const;
    boost::shared_ptr<SOBasisSet> basis1() const;
    boost::shared_ptr<SOBasisSet> basis2() const;
    boost::shared_ptr<SOBasisSet> basis3() const;
    boost::shared_ptr<SOBasisSet> basis4() const;

    const double *buffer() const { return buffer_; }

    template<typename TwoBodySOIntFunctor>
    void compute_shell(int, int, int, int, TwoBodySOIntFunctor& body);

    template<typename TwoBodySOIntFunctor>
    void compute_shell_smart(int, int, int, int, TwoBodySOIntFunctor& body);
};

template<typename TwoBodySOIntFunctor>
void TwoBodySOInt::compute_shell(int ish, int jsh, int ksh, int lsh, TwoBodySOIntFunctor& body)
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
    int irrepoff[8];

    memset(irrepoff, 0, sizeof(int) * 8);

    for (int h=1; h<b1_->nirrep(); ++h) {
        irrepoff[h] = irrepoff[h-1] + b1_->nfunction_in_irrep(h-1);
    }

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
//                        int irel = b1_->function_within_irrep(ish, isofunc);
//                        int iabs = irrepoff[ifunc.irrep] + irel;

                        for (int jtr=0; jtr<s2.nfunc; jtr++) {
                            const SOTransformFunction &jfunc = s2.func[jtr];
                            double jcoef = jfunc.coef * icoef;
                            int jaofunc = jfunc.aofunc;
                            int jsofunc = b2_->function_offset_within_shell(jsh,
                                                                            jfunc.irrep)
                                    + jfunc.sofunc;
                            int jaooff = iaooff*nao2 + jaofunc;
                            int jsooff = isooff*nso2 + jsofunc;
//                            int jrel = b2_->function_within_irrep(jsh, jsofunc);
//                            int jabs = irrepoff[jfunc.irrep] + jrel;

                            for (int ktr=0; ktr<s3.nfunc; ktr++) {
                                const SOTransformFunction &kfunc = s3.func[ktr];
                                double kcoef = kfunc.coef * jcoef;
                                int kaofunc = kfunc.aofunc;
                                int ksofunc = b3_->function_offset_within_shell(ksh,
                                                                                kfunc.irrep)
                                        + kfunc.sofunc;
                                int kaooff = jaooff*nao3 + kaofunc;
                                int ksooff = jsooff*nso3 + ksofunc;
//                                int krel = b3_->function_within_irrep(ksh, ksofunc);
//                                int kabs = irrepoff[kfunc.irrep] + krel;

                                for (int ltr=0; ltr<s4.nfunc; ltr++) {
                                    const SOTransformFunction &lfunc = s4.func[ltr];
                                    double lcoef = lfunc.coef * kcoef;
                                    int laofunc = lfunc.aofunc;
                                    int lsofunc = b4_->function_offset_within_shell(lsh,
                                                                                    lfunc.irrep)
                                            + lfunc.sofunc;
                                    int laooff = kaooff*nao4 + laofunc;
                                    int lsooff = ksooff*nso4 + lsofunc;
//                                    int lrel = b4_->function_within_irrep(lsh, lsofunc);
//                                    int labs = irrepoff[lfunc.irrep] + lrel;

                                    // If you're doing the two-stage SO integral uncomment the next line
                                    buffer_[lsooff] += lcoef * aobuff[laooff];

                                    // ---- Begin new code ----
//                                    double partial_value = lcoef * aobuff[laooff];
//                                    if (fabs(partial_value) > 1.0e-14) {
//                                        int iiabs = iabs;
//                                        int jjabs = jabs;
//                                        int kkabs = kabs;
//                                        int llabs = labs;

//                                        int iiirrep = ifunc.irrep;
//                                        int jjirrep = jfunc.irrep;
//                                        int kkirrep = kfunc.irrep;
//                                        int llirrep = lfunc.irrep;

//                                        int iirel = irel;
//                                        int jjrel = jrel;
//                                        int kkrel = krel;
//                                        int llrel = lrel;

////                                        fprintf(outfile, "original tr=%d %d %d %d (%d %d %d %d)\n", itr, jtr, ktr, ltr, iiabs, jjabs, kkabs, llabs);
//                                        if (ish == jsh) {
//                                            if (iabs < jabs)
//                                            continue;

//                                            if (ksh == lsh) {
//                                                if (kabs < labs)
//                                                continue;
//                                                if (INDEX2(iabs, jabs) < INDEX2(kabs, labs)) {
//                                                    if (ish == ksh)   // IIII case
//                                                    continue;
//                                                    else {            // IIJJ case
//                                                        SWAP_INDEX(ii, kk);
//                                                        SWAP_INDEX(jj, ll);
//                                                    }
//                                                }
//                                            }
//                                            else{                     // IIJK case
//                                                if (labs > kabs) {
//                                                    SWAP_INDEX(kk, ll);
//                                                }
//                                                if (INDEX2(iabs, jabs) < INDEX2(kabs, labs)) {
//                                                    SWAP_INDEX(ii, kk);
//                                                    SWAP_INDEX(jj, ll);
//                                                }
//                                            }
//                                        }
//                                        else {
//                                            if (ksh == lsh) {         // IJKK case
//                                                if (kabs < labs)
//                                                continue;
//                                                if (iabs < jabs) {
//                                                    SWAP_INDEX(ii, jj);
//                                                }
//                                                if (INDEX2(iabs, jabs) < INDEX2(kabs, labs)) {
//                                                    SWAP_INDEX(ii, kk);
//                                                    SWAP_INDEX(jj, ll);
//                                                }
//                                            }
//                                            else {                   // IJIJ case
//                                                if (ish == ksh && jsh == lsh && INDEX2(iabs, jabs) < INDEX2(kabs, labs))
//                                                continue;
//                                                // IJKL case
//                                                if (iabs < jabs) {
//                                                    SWAP_INDEX(ii, jj);
//                                                }
//                                                if (kabs < labs) {
//                                                    SWAP_INDEX(kk, ll);
//                                                }
//                                                if (INDEX2(iabs, jabs) < INDEX2(kabs, labs)) {
//                                                    SWAP_INDEX(ii, kk);
//                                                    SWAP_INDEX(jj, ll);
//                                                }
//                                            }
//                                        }

//                                        // func off/on
//                                        body(iiabs, jjabs, kkabs, llabs,
//                                             iiirrep, iirel,
//                                             jjirrep, jjrel,
//                                             kkirrep, kkrel,
//                                             llirrep, llrel,
//                                             partial_value);
//                                    }

//                                    // ---- End new code   ----

//                                    if (fabs(aobuff[laooff]*lcoef) > 1.0e-10) {
//                                        fprintf(outfile, "!(%d %d|%d %d) += %8.5lf * (%d %d|%d %d): %8.5lf (laoff = %d) -> %8.5lf (lsoff = %d)\n",
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

    // OK, I'm playing around here. The following ONLY WORKS if the basis sets are the same
//    if (ish == jsh && ish == ksh && ish == lsh)
//        provide_IIII(ish, jsh, ksh, lsh, body);
//    else if (ish == jsh && ksh == lsh)
//        provide_IIJJ(ish, jsh, ksh, lsh, body);
    provide_IJKL(ish, jsh, ksh, lsh, body);
}

#if 0
template<typename TwoBodySOIntFunctor>
void TwoBodySOInt::provide_IIII(int ish, int jsh, int ksh, int lsh, TwoBodySOIntFunctor& body)
{
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

    const SOTransformShell &s1 = t1.aoshell[0];

    int iirrepoff[8];

    memset(iirrepoff, 0, sizeof(int) * 8);

    for (int h=1; h<b1_->nirrep(); ++h) {
        iirrepoff[h] = iirrepoff[h-1] + b1_->nfunction_in_irrep(h-1);
    }

    for (int itr=0; itr<s1.nfunc; itr++) {
        const SOTransformFunction &ifunc = s1.func[itr];
        int isofunc = b1_->function_offset_within_shell(ish,
                                                        ifunc.irrep);
        int irel = b1_->function_within_irrep(ish, isofunc);
        isofunc += ifunc.sofunc;
        int isooff = isofunc;
        int iabs = iirrepoff[ifunc.irrep] + irel;

        for (int jtr=0; jtr<=itr; jtr++) {
            const SOTransformFunction &jfunc = s1.func[jtr];
            int jsofunc = b2_->function_offset_within_shell(jsh,
                                                            jfunc.irrep);
            int jrel = b2_->function_within_irrep(jsh, jsofunc);
            jsofunc += jfunc.sofunc;
            int jsooff = isooff*nso2 + jsofunc;
            int jabs = iirrepoff[jfunc.irrep] + jrel;

            for (int ktr=0; ktr<=itr; ktr++) {
                const SOTransformFunction &kfunc = s1.func[ktr];
                int ksofunc = b3_->function_offset_within_shell(ksh,
                                                                kfunc.irrep);
                int krel = b3_->function_within_irrep(ksh, ksofunc);
                ksofunc += kfunc.sofunc;
                int ksooff = jsooff*nso3 + ksofunc;
                int kabs = iirrepoff[kfunc.irrep] + krel;

                for (int ltr=0; ltr<=ktr; ltr++) {
                    const SOTransformFunction &lfunc = s1.func[ltr];
                    int lsofunc = b4_->function_offset_within_shell(lsh,
                                                                    lfunc.irrep);
                    int lrel = b4_->function_within_irrep(lsh, lsofunc);
                    lsofunc += lfunc.sofunc;
                    int lsooff = ksooff*nso4 + lsofunc;
                    int labs = iirrepoff[lfunc.irrep] + lrel;

                    int iiabs = iabs;
                    int jjabs = jabs;
                    int kkabs = kabs;
                    int llabs = labs;

                    int iiirrep = ifunc.irrep;
                    int jjirrep = jfunc.irrep;
                    int kkirrep = kfunc.irrep;
                    int llirrep = lfunc.irrep;

                    int iirel = irel;
                    int jjrel = jrel;
                    int kkrel = krel;
                    int llrel = lrel;

                    if (fabs(buffer_[lsooff]) > 1.0e-16) {
//                        fprintf(outfile, "!(%d %d|%d %d) : %8.5lf (lsoff = %d)\n",
//                                isofunc, jsofunc, ksofunc, lsofunc,
//                                buffer_[lsooff], lsooff);
                        // func off/on
                        body(iiabs, jjabs, kkabs, llabs,
                             iiirrep, iirel,
                             jjirrep, jjrel,
                             kkirrep, kkrel,
                             llirrep, llrel,
                             buffer_[lsooff]);
                    }
                }
            }
        }
    }
}
#endif

#if 0
template<typename TwoBodySOIntFunctor>
void TwoBodySOInt::provide_IIJJ(int ish, int jsh, int ksh, int lsh, TwoBodySOIntFunctor& body)
{
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

    const SOTransformShell &s1 = t1.aoshell[0];
    const SOTransformShell &s3 = t3.aoshell[0];

    int iirrepoff[8];
    int jirrepoff[8];

    memset(iirrepoff, 0, sizeof(int) * 8);
    memset(jirrepoff, 0, sizeof(int) * 8);

    for (int h=1; h<b1_->nirrep(); ++h) {
        iirrepoff[h] = iirrepoff[h-1] + b1_->nfunction_in_irrep(h-1);
        jirrepoff[h] = jirrepoff[h-1] + b3_->nfunction_in_irrep(h-1);
    }

    for (int itr=0; itr<s1.nfunc; itr++) {
        const SOTransformFunction &ifunc = s1.func[itr];
        int isofunc = b1_->function_offset_within_shell(ish,
                                                        ifunc.irrep) + ifunc.sofunc;
        int irel = b1_->function_within_irrep(ish, isofunc);
//        isofunc += ifunc.sofunc;
        int isooff = isofunc;
        int iabs = iirrepoff[ifunc.irrep] + irel;

        for (int jtr=0; jtr<=itr; jtr++) {
            const SOTransformFunction &jfunc = s1.func[jtr];
            int jsofunc = b2_->function_offset_within_shell(jsh,
                                                            jfunc.irrep) + jfunc.sofunc;
            int jrel = b2_->function_within_irrep(jsh, jsofunc);
//            jsofunc += jfunc.sofunc;
            int jsooff = isooff*nso2 + jsofunc;
            int jabs = iirrepoff[jfunc.irrep] + jrel;

            for (int ktr=0; ktr<s3.nfunc; ktr++) {
                const SOTransformFunction &kfunc = s3.func[ktr];
                int ksofunc = b3_->function_offset_within_shell(ksh,
                                                                kfunc.irrep) + kfunc.sofunc;
                int krel = b3_->function_within_irrep(ksh, ksofunc);
//                ksofunc += kfunc.sofunc;
                int ksooff = jsooff*nso3 + ksofunc;
                int kabs = jirrepoff[kfunc.irrep] + krel;

                for (int ltr=0; ltr<=ktr; ltr++) {
                    const SOTransformFunction &lfunc = s3.func[ltr];
                    int lsofunc = b4_->function_offset_within_shell(lsh,
                                                                    lfunc.irrep) + lfunc.sofunc;
                    int lrel = b4_->function_within_irrep(lsh, lsofunc);
//                    lsofunc += lfunc.sofunc;
                    int lsooff = ksooff*nso4 + lsofunc;
                    int labs = jirrepoff[lfunc.irrep] + lrel;

                    int iiabs = iabs;
                    int jjabs = jabs;
                    int kkabs = kabs;
                    int llabs = labs;

                    int iiirrep = ifunc.irrep;
                    int jjirrep = jfunc.irrep;
                    int kkirrep = kfunc.irrep;
                    int llirrep = lfunc.irrep;

                    int iirel = irel;
                    int jjrel = jrel;
                    int kkrel = krel;
                    int llrel = lrel;

//                    fprintf(outfile, "ifunc.sofunc = %d\n", ifunc.sofunc);
//                    fprintf(outfile, "jfunc.sofunc = %d\n", jfunc.sofunc);
//                    fprintf(outfile, "kfunc.sofunc = %d\n", kfunc.sofunc);
//                    fprintf(outfile, "lfunc.sofunc = %d\n", lfunc.sofunc);

                    if (fabs(buffer_[lsooff]) > 1.0e-16) {
//                        fprintf(outfile, "!(%d %d|%d %d) : %8.5lf (lsoff = %d)\n",
//                                isofunc, jsofunc, ksofunc, lsofunc,
//                                buffer_[lsooff], lsooff);
                        if (iiabs >= kkabs) {
                            // func off/on
                            body(iiabs, jjabs, kkabs, llabs,
                                 iiirrep, iirel,
                                 jjirrep, jjrel,
                                 kkirrep, kkrel,
                                 llirrep, llrel,
                                 buffer_[lsooff]);
                        }
                        else {
                            body(kkabs, llabs, iiabs, jjabs,
                                 kkirrep, kkrel,
                                 llirrep, llrel,
                                 iiirrep, iirel,
                                 jjirrep, jjrel,
                                 buffer_[lsooff]);
                        }
                    }
                }
            }
        }
    }
}
#endif

template<typename TwoBodySOIntFunctor>
void TwoBodySOInt::provide_IJKL(int ish, int jsh, int ksh, int lsh, TwoBodySOIntFunctor& body)
{
    const double *aobuff = tb_->buffer();

    const SOTransform &t1 = b1_->trans(ish);
    const SOTransform &t2 = b2_->trans(jsh);
    const SOTransform &t3 = b3_->trans(ksh);
    const SOTransform &t4 = b4_->trans(lsh);

    int nso1 = b1_->nfunction(ish);
    int nso2 = b2_->nfunction(jsh);
    int nso3 = b3_->nfunction(ksh);
    int nso4 = b4_->nfunction(lsh);

    int n1 = b1_->nfunction(ish);
    int n2 = b2_->nfunction(jsh);
    int n3 = b3_->nfunction(ksh);
    int n4 = b4_->nfunction(lsh);

//    fprintf(outfile, "n1 = %d, n2 = %d, n3 = %d, n4 = %d\n", n1, n2, n3, n4); fflush(outfile);

    const SOTransformShell &s1 = t1.aoshell[0];
    const SOTransformShell &s2 = t2.aoshell[0];
    const SOTransformShell &s3 = t3.aoshell[0];
    const SOTransformShell &s4 = t4.aoshell[0];

    int iirrepoff[8];
    int jirrepoff[8];
    int kirrepoff[8];
    int lirrepoff[8];

    memset(iirrepoff, 0, sizeof(int) * 8);
    memset(jirrepoff, 0, sizeof(int) * 8);
    memset(kirrepoff, 0, sizeof(int) * 8);
    memset(lirrepoff, 0, sizeof(int) * 8);

    for (int h=1; h<b1_->nirrep(); ++h) {
        iirrepoff[h] = iirrepoff[h-1] + b1_->nfunction_in_irrep(h-1);
        jirrepoff[h] = jirrepoff[h-1] + b2_->nfunction_in_irrep(h-1);
        kirrepoff[h] = kirrepoff[h-1] + b3_->nfunction_in_irrep(h-1);
        lirrepoff[h] = lirrepoff[h-1] + b4_->nfunction_in_irrep(h-1);
    }

    int itr, itrfunc;
    int jtr, jtrfunc;
    int ktr, ktrfunc;
    int ltr, ltrfunc;

    for (itr=0, itrfunc=0; itr<n1; itr++, itrfunc++) {

        int ifunc = b1_->function(ish) + itr;
        int isym = b1_->irrep(ifunc);
        int irel = b1_->function_within_irrep(ifunc);
        int iabs = iirrepoff[isym] + irel;
        int isooff = itr;

        for (jtr=0; jtr<n2; jtr++) {

            int jfunc = b2_->function(jsh) + jtr;
            int jsym = b2_->irrep(jfunc);
            int jrel = b2_->function_within_irrep(jfunc);
            int jabs = jirrepoff[jsym] + jrel;
            int jsooff = isooff*nso2 + jtr;

            for (ktr=0; ktr<n3; ktr++) {

                int kfunc = b3_->function(ksh) + ktr;
                int ksym = b3_->irrep(kfunc);
                int krel = b3_->function_within_irrep(kfunc);
                int kabs = kirrepoff[ksym] + krel;
                int ksooff = jsooff*nso3 + ktr;

                for (ltr=0; ltr<n4; ltr++) {

                    int lfunc = b4_->function(lsh) + ltr;
                    int lsym = b4_->irrep(lfunc);
                    int lrel = b4_->function_within_irrep(lfunc);
                    int labs = lirrepoff[lsym] + lrel;
                    int lsooff = ksooff*nso4 + ltr;

                    int iiabs = iabs;
                    int jjabs = jabs;
                    int kkabs = kabs;
                    int llabs = labs;

                    int iiirrep = isym;
                    int jjirrep = jsym;
                    int kkirrep = ksym;
                    int llirrep = lsym;

                    int iirel = irel;
                    int jjrel = jrel;
                    int kkrel = krel;
                    int llrel = lrel;

//                    fprintf(outfile, "ifunc.sofunc = %d\n", ifunc.sofunc);
//                    fprintf(outfile, "jfunc.sofunc = %d\n", jfunc.sofunc);
//                    fprintf(outfile, "kfunc.sofunc = %d\n", kfunc.sofunc);
//                    fprintf(outfile, "lfunc.sofunc = %d\n", lfunc.sofunc);

                    if (fabs(buffer_[lsooff]) > 1.0e-14) {
//                        fprintf(outfile, "original tr=%d %d %d %d (%d %d %d %d)\n", itr, jtr, ktr, ltr, iiabs, jjabs, kkabs, llabs);
                        if (ish == jsh) {
                            if (iabs < jabs)
                            continue;

                            if (ksh == lsh) {
                                if (kabs < labs)
                                continue;
                                if (INDEX2(iabs, jabs) < INDEX2(kabs, labs)) {
                                    if (ish == ksh)   // IIII case
                                    continue;
                                    else {            // IIJJ case
                                        SWAP_INDEX(ii, kk);
                                        SWAP_INDEX(jj, ll);
                                    }
                                }
                            }
                            else{                     // IIJK case
                                if (labs > kabs) {
                                    SWAP_INDEX(kk, ll);
                                }
                                if (INDEX2(iabs, jabs) < INDEX2(kabs, labs)) {
                                    SWAP_INDEX(ii, kk);
                                    SWAP_INDEX(jj, ll);
                                }
                            }
                        }
                        else {
                            if (ksh == lsh) {         // IJKK case
                                if (kabs < labs)
                                continue;
                                if (iabs < jabs) {
                                    SWAP_INDEX(ii, jj);
                                }
                                if (INDEX2(iabs, jabs) < INDEX2(kabs, labs)) {
                                    SWAP_INDEX(ii, kk);
                                    SWAP_INDEX(jj, ll);
                                }
                            }
                            else {                   // IJIJ case
                                if (ish == ksh && jsh == lsh && INDEX2(iabs, jabs) < INDEX2(kabs, labs))
                                continue;
                                // IJKL case
                                if (iabs < jabs) {
                                    SWAP_INDEX(ii, jj);
                                }
                                if (kabs < labs) {
                                    SWAP_INDEX(kk, ll);
                                }
                                if (INDEX2(iabs, jabs) < INDEX2(kabs, labs)) {
                                    SWAP_INDEX(ii, kk);
                                    SWAP_INDEX(jj, ll);
                                }
                            }
                        }

                        // func off/on
                        body(iiabs, jjabs, kkabs, llabs,
                             iiirrep, iirel,
                             jjirrep, jjrel,
                             kkirrep, kkrel,
                             llirrep, llrel,
                             buffer_[lsooff]);
                    }
                }
            }
        }
    }
}

//template<typename TwoBodySOIntFunctor>
//void TwoBodySOInt::compute_shell_smart(int ish, int jsh, int ksh, int lsh, TwoBodySOIntFunctor& body)
//{
//    // ish ... are unique shells.
//    fprintf(outfile, "\tcomputing SO shell quartet (%d %d %d %d)\n", ish, jsh, ksh, lsh);

//    const double *aobuff = tb_->buffer();
//    boost::shared_ptr<PetiteList> petite = b1_->petitelist();

//    const SOTransform &t1 = b1_->trans(ish);
//    const SOTransform &t2 = b2_->trans(jsh);
//    const SOTransform &t3 = b3_->trans(ksh);
//    const SOTransform &t4 = b4_->trans(lsh);

//    int nso1 = b1_->nfunction(ish);
//    int nso2 = b2_->nfunction(jsh);
//    int nso3 = b3_->nfunction(ksh);
//    int nso4 = b4_->nfunction(lsh);

//    int nao1 = b1_->naofunction(ish);
//    int nao2 = b2_->naofunction(jsh);
//    int nao3 = b3_->naofunction(ksh);
//    int nao4 = b4_->naofunction(lsh);

////    fprintf(outfile, "nao1 = %d nao2 = %d nao3 = %d nao4 = %d\n", nao1, nao2, nao3, nao4);

//    memset(buffer_, 0, 16*nao1*nao2*nao3*nao4*sizeof(double));

//    // Map unique shell to shell
//    int si = t1.aoshell[t1.naoshell-1].aoshell;
//    int sj = t2.aoshell[t2.naoshell-1].aoshell;
//    int sk = t3.aoshell[t3.naoshell-1].aoshell;
//    int sl = t4.aoshell[t4.naoshell-1].aoshell;

//    // Order of the group.
//    double order = (double)petite->order();

//    fprintf(outfile, "si = %d sj = %d sk = %d sl = %d\n", si, sj, sk, sl);

//    for (int i=0; i<t1.naoshell; ++i) {
//        const SOTransformShell& s1 = t1.aoshell[i];

//        if (!petite->in_p1(s1.aoshell)) {
//            fprintf(outfile, "skipping i = %d\n", i);
//            continue;
//        }

//        int atom_i = integral_->basis1()->shell(s1.aoshell)->ncenter();
//        int stab_i = petite->stabilizer(atom_i);

//        for (int j=0; j<t2.naoshell; ++j) {
//            const SOTransformShell& s2 = t2.aoshell[j];

//            if (!petite->in_p2(s1.aoshell, s2.aoshell))
//                continue;

//            int atom_j = integral_->basis2()->shell(s2.aoshell)->ncenter();
//            int stab_j = petite->stabilizer(atom_j);

//            int stab_ij = petite->GnG(stab_i, stab_j);

//            for (int k=0; k<t3.naoshell; ++k) {
//                const SOTransformShell& s3 = t3.aoshell[k];

//                int atom_k = integral_->basis3()->shell(s3.aoshell)->ncenter();
//                int stab_k = petite->stabilizer(atom_k);

//                for (int l=0; l<t4.naoshell; ++l) {
//                    const SOTransformShell& s4 = t4.aoshell[l];

//                    if (!petite->in_p4(s1.aoshell, s2.aoshell, s3.aoshell, s4.aoshell))
//                        continue;

//                    int atom_l = integral_->basis4()->shell(s4.aoshell)->ncenter();
//                    int stab_l = petite->stabilizer(atom_l);

//                    int stab_kl = petite->GnG(stab_k, stab_l);

//                    double lambda_T = order / (double)petite->dcr_degeneracy(stab_ij, stab_kl);

//                    fprintf(outfile, "i = %d j = %d k = %d l = %d lambda_T = %8.5f\n",
//                            i, j, k, l, lambda_T);

//                    tb_->compute_shell(s1.aoshell, s2.aoshell, s3.aoshell, s4.aoshell);

//                    for (int itr=0; itr<s1.nfunc; itr++) {
//                        const SOTransformFunction &ifunc = s1.func[itr];
//                        double icoef = ifunc.coef;
//                        int iaofunc = ifunc.aofunc;
//                        int isofunc = b1_->function_offset_within_shell(ish,
//                                                                        ifunc.irrep)
//                                + ifunc.sofunc;
//                        int iaooff = iaofunc;
//                        int isooff = isofunc;
//                        for (int jtr=0; jtr<s2.nfunc; jtr++) {
//                            const SOTransformFunction &jfunc = s2.func[jtr];
//                            double jcoef = jfunc.coef * icoef;
//                            int jaofunc = jfunc.aofunc;
//                            int jsofunc = b2_->function_offset_within_shell(jsh,
//                                                                            jfunc.irrep)
//                                    + jfunc.sofunc;
//                            int jaooff = iaooff*nao2 + jaofunc;
//                            int jsooff = isooff*nso2 + jsofunc;
//                            for (int ktr=0; ktr<s3.nfunc; ktr++) {
//                                const SOTransformFunction &kfunc = s3.func[ktr];
//                                double kcoef = kfunc.coef * jcoef;
//                                int kaofunc = kfunc.aofunc;
//                                int ksofunc = b3_->function_offset_within_shell(ksh,
//                                                                                kfunc.irrep)
//                                        + kfunc.sofunc;
//                                int kaooff = jaooff*nao3 + kaofunc;
//                                int ksooff = jsooff*nso3 + ksofunc;
//                                for (int ltr=0; ltr<s4.nfunc; ltr++) {
//                                    const SOTransformFunction &lfunc = s4.func[ltr];
//                                    double lcoef = lfunc.coef * kcoef;
//                                    int laofunc = lfunc.aofunc;
//                                    int lsofunc = b4_->function_offset_within_shell(lsh,
//                                                                                    lfunc.irrep)
//                                            + lfunc.sofunc;
//                                    int laooff = kaooff*nao4 + laofunc;
//                                    int lsooff = ksooff*nso4 + lsofunc;

//                                    buffer_[lsooff] += lambda_T * lcoef * aobuff[laooff];
////                                    if (fabs(aobuff[laooff]*lcoef) > 1.0e-10) {
////                                        fprintf(outfile, "(%d %d|%d %d) += %8.5lf * (%d %d|%d %d): %8.5lf (laoff = %d) -> %8.5lf (lsoff = %d)\n",
////                                                isofunc, jsofunc, ksofunc, lsofunc, lcoef,
////                                                iaofunc, jaofunc, kaofunc, laofunc,
////                                                aobuff[laooff], laooff,
////                                                buffer_[lsooff], lsooff);
////                                    }
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//    }

//    const SOTransformShell &s1 = t1.aoshell[0];
//    const SOTransformShell &s2 = t2.aoshell[0];
//    const SOTransformShell &s3 = t3.aoshell[0];
//    const SOTransformShell &s4 = t4.aoshell[0];

//    for (int itr=0; itr<s1.nfunc; itr++) {
//        const SOTransformFunction &ifunc = s1.func[itr];
//        double icoef = ifunc.coef;
//        int iaofunc = ifunc.aofunc;
//        int isofunc = b1_->function_offset_within_shell(ish,
//                                                        ifunc.irrep)
//                + ifunc.sofunc;
//        int iaooff = iaofunc;
//        int isooff = isofunc;
//        for (int jtr=0; jtr<s2.nfunc; jtr++) {
//            const SOTransformFunction &jfunc = s2.func[jtr];
//            double jcoef = jfunc.coef * icoef;
//            int jaofunc = jfunc.aofunc;
//            int jsofunc = b2_->function_offset_within_shell(jsh,
//                                                            jfunc.irrep)
//                    + jfunc.sofunc;
//            int jaooff = iaooff*nao2 + jaofunc;
//            int jsooff = isooff*nso2 + jsofunc;
//            for (int ktr=0; ktr<s3.nfunc; ktr++) {
//                const SOTransformFunction &kfunc = s3.func[ktr];
//                double kcoef = kfunc.coef * jcoef;
//                int kaofunc = kfunc.aofunc;
//                int ksofunc = b3_->function_offset_within_shell(ksh,
//                                                                kfunc.irrep)
//                        + kfunc.sofunc;
//                int kaooff = jaooff*nao3 + kaofunc;
//                int ksooff = jsooff*nso3 + ksofunc;
//                for (int ltr=0; ltr<s4.nfunc; ltr++) {
//                    const SOTransformFunction &lfunc = s4.func[ltr];
//                    double lcoef = lfunc.coef * kcoef;
//                    int laofunc = lfunc.aofunc;
//                    int lsofunc = b4_->function_offset_within_shell(lsh,
//                                                                    lfunc.irrep)
//                            + lfunc.sofunc;
//                    int laooff = kaooff*nao4 + laofunc;
//                    int lsooff = ksooff*nso4 + lsofunc;

////                    buffer_[lsooff] += lcoef * aobuff[laooff];
//                    if (fabs(buffer_[lsooff]) > 1.0e-16)
//                        body(ifunc.irrep, isofunc,
//                             jfunc.irrep, jsofunc,
//                             kfunc.irrep, ksofunc,
//                             lfunc.irrep, lsofunc,
//                             buffer_[lsooff]);
//                }
//            }
//        }
//    }
//}

}

#endif // _psi_src_lib_libmints_sointegral_h_

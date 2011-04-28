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
#include "cdsalclist.h"

#ifdef MINTS_TIMER
#include <libqt/qt.h>
#endif

#include <vector>

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
    int deriv_;

    boost::shared_ptr<SOBasisSet> b1_;
    boost::shared_ptr<SOBasisSet> b2_;

    size_t size_;
    double *buffer_;

    void common_init();

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
     *
     * \param result Where the integrals are going.
     */
    void compute(boost::shared_ptr<Matrix> result);

    /**
     * Computes one-electron integral matrices. Should be able to handle multipole operators
     *
     * \param results Where the integrals are going.
     */
    void compute(std::vector<boost::shared_ptr<Matrix> > results);

    /**
     * Computes one-electron integral derivative matrices.
     *
     * \param result Where the integral derivatives are going.
     * \param cdsalcs The Cartesian displacement SALCs that you are interested in.
     */
    void compute_deriv1(std::vector<boost::shared_ptr<Matrix> > result,
                        const CdSalcList& cdsalcs);
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

    size_t size_;
    double *buffer_;

    int iirrepoff_[8], jirrepoff_[8], kirrepoff_[8], lirrepoff_[8];

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
    void compute_shell(const SOShellCombinationsIterator& shellIter, TwoBodySOIntFunctor& body) {
        compute_shell(shellIter.p(), shellIter.q(), shellIter.r(), shellIter.s(),
                      body);
    }

    template<typename TwoBodySOIntFunctor>
    void compute_shell(int, int, int, int, TwoBodySOIntFunctor& body);

    template<typename TwoBodySOIntFunctor>
    void compute_shell_smart(int, int, int, int, TwoBodySOIntFunctor& body);
};

template<typename TwoBodySOIntFunctor>
void TwoBodySOInt::compute_shell(int ish, int jsh, int ksh, int lsh, TwoBodySOIntFunctor& body)
{
#ifdef MINTS_TIMER
    timer_on("TwoBodySOInt::compute_shell overall");
#endif

#ifdef MINTS_TIMER
    timer_on("TwoBodySOInt::compute_shell setup");
#endif // MINTS_TIMER

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

#ifdef MINTS_TIMER
    timer_on("TwoBodySOInt::compute_shell zero buffer");
#endif // MINTS_TIMER

    memset(buffer_, 0, nso1*nso2*nso3*nso4*sizeof(double));

#ifdef MINTS_TIMER
    timer_off("TwoBodySOInt::compute_shell zero buffer");
#endif // MINTS_TIMER

#ifdef MINTS_TIMER
    timer_off("TwoBodySOInt::compute_shell setup");
#endif // MINTS_TIMER

#ifdef MINTS_TIMER
    timer_on("TwoBodySOInt::compute_shell full shell transform");
#endif // MINTS_TIMER

    // loop through the ao shells that make up this so shell
    for (int i=0; i<t1.naoshell; i++) {
        const SOTransformShell &s1 = t1.aoshell[i];
        for (int j=0; j<t2.naoshell; j++) {
            const SOTransformShell &s2 = t2.aoshell[j];
            for (int k=0; k<t3.naoshell; k++) {
                const SOTransformShell &s3 = t3.aoshell[k];
                for (int l=0; l<t4.naoshell; l++) {
                    const SOTransformShell &s4 = t4.aoshell[l];

#ifdef MINTS_TIMER
    timer_on("TwoBodySOInt::compute_shell AO eri overhead");
#endif // MINTS_TIMER
                    tb_->compute_shell(s1.aoshell, s2.aoshell, s3.aoshell, s4.aoshell);
#ifdef MINTS_TIMER
    timer_off("TwoBodySOInt::compute_shell AO eri overhead");
    timer_on("TwoBodySOInt::compute_shell AO->SO transform");
#endif // MINTS_TIMER

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

                                    // If you're doing the two-stage SO integral uncomment the next line
                                    buffer_[lsooff] += lcoef * aobuff[laooff];
                                }
                            }
                        }
                    }

#ifdef MINTS_TIMER
    timer_off("TwoBodySOInt::compute_shell AO->SO transform");
#endif // MINTS_TIMER

                }
            }
        }
    }

#ifdef MINTS_TIMER
    timer_off("TwoBodySOInt::compute_shell full shell transform");
#endif // MINTS_TIMER

    provide_IJKL(ish, jsh, ksh, lsh, body);

#ifdef MINTS_TIMER
    timer_off("TwoBodySOInt::compute_shell overall");
#endif
}

template<typename TwoBodySOIntFunctor>
void TwoBodySOInt::provide_IJKL(int ish, int jsh, int ksh, int lsh, TwoBodySOIntFunctor& body)
{
#ifdef MINTS_TIMER
    timer_on("TwoBodySOInt::provide_IJKL overall");
#endif

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

    const SOTransformShell &s1 = t1.aoshell[0];
    const SOTransformShell &s2 = t2.aoshell[0];
    const SOTransformShell &s3 = t3.aoshell[0];
    const SOTransformShell &s4 = t4.aoshell[0];

    int itr, itrfunc;
    int jtr, jtrfunc;
    int ktr, ktrfunc;
    int ltr, ltrfunc;

    for (itr=0, itrfunc=0; itr<n1; itr++, itrfunc++) {

        int ifunc = b1_->function(ish) + itr;
        int isym = b1_->irrep(ifunc);
        int irel = b1_->function_within_irrep(ifunc);
        int iabs = iirrepoff_[isym] + irel;
        int isooff = itr;

        for (jtr=0; jtr<n2; jtr++) {

            int jfunc = b2_->function(jsh) + jtr;
            int jsym = b2_->irrep(jfunc);
            int jrel = b2_->function_within_irrep(jfunc);
            int jabs = jirrepoff_[jsym] + jrel;
            int jsooff = isooff*nso2 + jtr;

            for (ktr=0; ktr<n3; ktr++) {

                int kfunc = b3_->function(ksh) + ktr;
                int ksym = b3_->irrep(kfunc);
                int krel = b3_->function_within_irrep(kfunc);
                int kabs = kirrepoff_[ksym] + krel;
                int ksooff = jsooff*nso3 + ktr;

                for (ltr=0; ltr<n4; ltr++) {

                    int lfunc = b4_->function(lsh) + ltr;
                    int lsym = b4_->irrep(lfunc);
                    int lrel = b4_->function_within_irrep(lfunc);
                    int labs = lirrepoff_[lsym] + lrel;
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

                    if (fabs(buffer_[lsooff]) > 1.0e-14) {
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

#ifdef MINTS_TIMER
    timer_on("TwoBodySOInt::provide_IJKL functor");
#endif
                        // func off/on
                        body(iiabs, jjabs, kkabs, llabs,
                             iiirrep, iirel,
                             jjirrep, jjrel,
                             kkirrep, kkrel,
                             llirrep, llrel,
                             buffer_[lsooff]);
#ifdef MINTS_TIMER
    timer_off("TwoBodySOInt::provide_IJKL functor");
#endif
                    }
                }
            }
        }
    }
#ifdef MINTS_TIMER
    timer_off("TwoBodySOInt::provide_IJKL overall");
#endif
}
}

#endif // _psi_src_lib_libmints_sointegral_h_

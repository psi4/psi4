/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#ifndef _psi_src_lib_libmints_sointegral_h_
#define _psi_src_lib_libmints_sointegral_h_

#include "typedefs.h"

#include "onebody.h"
#include "twobody.h"
#include "basisset.h"
#include "integral.h"
#include "sobasis.h"
#include "gshell.h"
#include "petitelist.h"
#include "wavefunction.h"
#include "cdsalclist.h"
#include "dcd.h"


#include <libqt/qt.h>
#include <vector>

//#define DebugPrint 1

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef DebugPrint
#   define DebugPrint 0
#endif

#if DebugPrint
#define dprintf(...) outfile->Printf( __VA_ARGS__)
#else
#define dprintf(...)
#endif

namespace psi {

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
    std::vector<boost::shared_ptr<TwoBodyAOInt> > tb_;
    boost::shared_ptr<IntegralFactory> integral_;

    boost::shared_ptr<SOBasisSet> b1_;
    boost::shared_ptr<SOBasisSet> b2_;
    boost::shared_ptr<SOBasisSet> b3_;
    boost::shared_ptr<SOBasisSet> b4_;

    size_t size_;
    std::vector<double *> buffer_;
    std::vector<double *> temp_;
    std::vector<double *> temp2_;
    std::vector<double **> deriv_;

    int iirrepoff_[8], jirrepoff_[8], kirrepoff_[8], lirrepoff_[8];

    boost::shared_ptr<PetiteList> petite1_;
    boost::shared_ptr<PetiteList> petite2_;
    boost::shared_ptr<PetiteList> petite3_;
    boost::shared_ptr<PetiteList> petite4_;

    boost::shared_ptr<PointGroup> pg_;

    boost::shared_ptr<DCD> dcd_;

    bool only_totally_symmetric_;
    int nthread_;
    std::string comm_;
    int nproc_;
    int me_;

    double cutoff_;

    const CdSalcList* cdsalcs_;

    template<typename TwoBodySOIntFunctor>
    void provide_IJKL(int, int, int, int, TwoBodySOIntFunctor& body);

    template<typename TwoBodySOIntFunctor>
    void provide_IJKL_deriv1(int ish, int jsh, int ksh, int lsh, TwoBodySOIntFunctor& body);

    void common_init();
public:
    // Constructor, assuming 1 thread
    TwoBodySOInt(const boost::shared_ptr<TwoBodyAOInt>&,
                 const boost::shared_ptr<IntegralFactory>&);
    // Constructor, using vector of AO objects for threading
    TwoBodySOInt(const std::vector<boost::shared_ptr<TwoBodyAOInt> > &tb,
                 const boost::shared_ptr<IntegralFactory>& integral);
    TwoBodySOInt(const boost::shared_ptr<TwoBodyAOInt>& aoint,
                 const boost::shared_ptr<IntegralFactory>& intfac,
                 const CdSalcList& cdsalcs);
    TwoBodySOInt(const std::vector<boost::shared_ptr<TwoBodyAOInt> >& tb,
                 const boost::shared_ptr<IntegralFactory>& integral,
                 const CdSalcList& cdsalcs);

    virtual ~TwoBodySOInt();

    bool only_totally_symmetric() const { return only_totally_symmetric_; }
    void set_only_totally_symmetric(bool ots) { only_totally_symmetric_ = ots; }

    boost::shared_ptr<SOBasisSet> basis() const;
    boost::shared_ptr<SOBasisSet> basis1() const;
    boost::shared_ptr<SOBasisSet> basis2() const;
    boost::shared_ptr<SOBasisSet> basis3() const;
    boost::shared_ptr<SOBasisSet> basis4() const;

    const double *buffer(int thread=0) const { return buffer_[thread]; }

    void set_cutoff(double ints_tolerance) { cutoff_ = ints_tolerance; }

    // Normal integrals
    template<typename TwoBodySOIntFunctor>
    void compute_shell(const SOShellCombinationsIterator& shellIter, TwoBodySOIntFunctor& body) {
        compute_shell(shellIter.p(), shellIter.q(), shellIter.r(), shellIter.s(),
                      body);
    }

    template<typename TwoBodySOIntFunctor>
    void compute_shell(int, int, int, int, TwoBodySOIntFunctor& body);

    // User provides an iterator object and this function will walk through it.
    // Assumes serial run (nthread = 1)
    template<typename ShellIter, typename TwoBodySOIntFunctor>
    void compute_quartets(ShellIter &shellIter, TwoBodySOIntFunctor &body) {
        for (shellIter->first(); shellIter->is_done() == false; shellIter->next()) {
            this->compute_shell(shellIter->p(), shellIter->q(), shellIter->r(), shellIter->s(), body);
        }
    }

    // Compute integrals in parallel
    template<typename TwoBodySOIntFunctor>
    void compute_integrals(TwoBodySOIntFunctor &functor);

    // Derivative integrals
    // User provides an iterator object and this function will walk through it.
    template<typename ShellIter, typename TwoBodySOIntFunctor>
    void compute_quartets_deriv1(ShellIter &shellIter, TwoBodySOIntFunctor &body) {
        for (shellIter->first(); shellIter->is_done() == false; shellIter->next()) {
            compute_shell_deriv1(shellIter->p(), shellIter->q(), shellIter->r(), shellIter->s(), body);
        }
    }

    template<typename TwoBodySOIntFunctor>
    void compute_shell_deriv1(const SOShellCombinationsIterator& shellIter, TwoBodySOIntFunctor& body) {
        compute_shell_deriv1(shellIter.p(), shellIter.q(), shellIter.r(), shellIter.s(),
                             body);
    }

    template<typename TwoBodySOIntFunctor>
    void compute_shell_deriv1(int, int, int, int, TwoBodySOIntFunctor& body);

    // Compute integrals in parallel
    template<typename TwoBodySOIntFunctor>
    void compute_integrals_deriv1(TwoBodySOIntFunctor &functor);

    template<typename TwoBodySOIntFunctor>
    int compute_pq_pair_deriv1(const int &p, const int &q, const size_t &pair_number, const TwoBodySOIntFunctor &body) {

        const_cast<TwoBodySOIntFunctor &>(body).load_tpdm(pair_number);
        boost::shared_ptr<SO_RS_Iterator> shellIter(
                    new SO_RS_Iterator(p, q,
                                       b1_, b2_, b3_, b4_));

        compute_quartets_deriv1(shellIter, const_cast<TwoBodySOIntFunctor &>(body));

        return 0;
    }
};

template<typename TwoBodySOIntFunctor>
void TwoBodySOInt::compute_shell(int uish, int ujsh, int uksh, int ulsh, TwoBodySOIntFunctor& body)
{
    dprintf("uish %d, ujsh %d, uksh %d, ulsh %d\n", uish, ujsh, uksh, ulsh);

    int thread = 0;
    #ifdef _OPENMP
      thread = omp_get_thread_num();
    #endif
    //Old call WorldComm->thread_id(pthread_self());

    //timer_on("TwoBodySOInt::compute_shell overall");
    mints_timer_on("TwoBodySOInt::compute_shell overall");
    mints_timer_on("TwoBodySOInt::compute_shell setup");

    const double *aobuff = tb_[thread]->buffer();

    const SOTransform &t1 = b1_->sotrans(uish);
    const SOTransform &t2 = b2_->sotrans(ujsh);
    const SOTransform &t3 = b3_->sotrans(uksh);
    const SOTransform &t4 = b4_->sotrans(ulsh);

    const int nso1 = b1_->nfunction(uish);
    const int nso2 = b2_->nfunction(ujsh);
    const int nso3 = b3_->nfunction(uksh);
    const int nso4 = b4_->nfunction(ulsh);
    const size_t nso = nso1*nso2*nso3*nso4;

    const int nao2 = b2_->naofunction(ujsh);
    const int nao3 = b3_->naofunction(uksh);
    const int nao4 = b4_->naofunction(ulsh);

    const int iatom = tb_[thread]->basis1()->shell(t1.aoshell[0].aoshell).ncenter();
    const int jatom = tb_[thread]->basis2()->shell(t2.aoshell[0].aoshell).ncenter();
    const int katom = tb_[thread]->basis3()->shell(t3.aoshell[0].aoshell).ncenter();
    const int latom = tb_[thread]->basis4()->shell(t4.aoshell[0].aoshell).ncenter();

    int nirrep = b1_->nirrep();

    mints_timer_on("TwoBodySOInt::compute_shell zero buffer");

    ::memset(buffer_[thread], 0, nso*sizeof(double));

    mints_timer_off("TwoBodySOInt::compute_shell zero buffer");
    mints_timer_off("TwoBodySOInt::compute_shell setup");

    mints_timer_on("TwoBodySOInt::compute_shell full shell transform");

    // Get the atomic stablizer (the first symmetry operation that maps the atom
    // onto itself.

    // These 3 sections are not shell specific so we can just use petite1_
    const unsigned short istablizer = petite1_->stablizer(iatom);
    const unsigned short jstablizer = petite1_->stablizer(jatom);
    const unsigned short kstablizer = petite1_->stablizer(katom);
    const unsigned short lstablizer = petite1_->stablizer(latom);

    const int istabdense = dcd_->bits_to_dense_numbering(istablizer);
    const int jstabdense = dcd_->bits_to_dense_numbering(jstablizer);
    const int kstabdense = dcd_->bits_to_dense_numbering(kstablizer);
    const int lstabdense = dcd_->bits_to_dense_numbering(lstablizer);

    const int ijstablizer = dcd_->intersection(istabdense, jstabdense);
    const int klstablizer = dcd_->intersection(kstabdense, lstabdense);
    const int ijklstablizer = dcd_->intersection(ijstablizer, klstablizer);

    const int* R_list = dcd_->dcr(istabdense, jstabdense);
    const int* S_list = dcd_->dcr(kstabdense, lstabdense);
    const int* T_list = dcd_->dcr(ijstablizer, klstablizer);

    const int R_size = R_list[0];
    const int S_size = S_list[0];
    const int T_size = T_list[0];

    // Check with Andy on this:
    int lambda_T = petite1_->nirrep() / dcd_->subgroup_dimensions(ijklstablizer);

    std::vector<int> sj_arr, sk_arr, sl_arr;

    int si = petite1_->unique_shell_map(uish, 0);
    const int siatom = tb_[thread]->basis1()->shell(si).ncenter();

    dprintf("dcd %d", petite1_->group());
    dprintf("istab %d, jstab %d, kstab %d, lstab %d, ijstab %d, klstab %d\n", istabdense, jstabdense, kstabdense, lstabdense, ijstablizer, klstablizer);
    dprintf("R_size %d, S_size %d, T_size %d\n", R_size, S_size, T_size);

    for (int ij=1; ij <= R_size; ++ij) {
        int sj = petite2_->unique_shell_map(ujsh, R_list[ij]);
        const int sjatom = tb_[thread]->basis2()->shell(sj).ncenter();

        for (int ijkl=1; ijkl <= T_size; ++ijkl) {
            int sk = petite3_->unique_shell_map(uksh, T_list[ijkl]);
            int llsh = petite4_->unique_shell_map(ulsh, T_list[ijkl]);
            const int skatom = tb_[thread]->basis3()->shell(sk).ncenter();

            for (int kl=1; kl <= S_size; ++kl) {
                int sl = petite4_->shell_map(llsh, S_list[kl]);
                const int slatom = tb_[thread]->basis4()->shell(sl).ncenter();

                // Check AM
                int total_am = tb_[thread]->basis1()->shell(si).am() +
                               tb_[thread]->basis2()->shell(sj).am() +
                               tb_[thread]->basis3()->shell(sk).am() +
                               tb_[thread]->basis4()->shell(sl).am();

                if (!(total_am % 2) ||
                    (siatom != sjatom) ||
                    (sjatom != skatom) ||
                    (skatom != slatom)) {
                    sj_arr.push_back(sj);
                    sk_arr.push_back(sk);
                    sl_arr.push_back(sl);
                }
            }
        }
    }

    // Compute integral using si, sj_arr, sk_arr, sl_arr
    // Loop over unique quartets
    const AOTransform& s1 = b1_->aotrans(si);

    const unsigned short *ifuncpi = s1.nfuncpi;

    for (size_t n=0; n<sj_arr.size(); ++n) {
        int sj = sj_arr[n];
        int sk = sk_arr[n];
        int sl = sl_arr[n];

        const AOTransform& s2 = b2_->aotrans(sj);
        const AOTransform& s3 = b3_->aotrans(sk);
        const AOTransform& s4 = b4_->aotrans(sl);

        const unsigned short *jfuncpi = s2.nfuncpi;
        const unsigned short *kfuncpi = s3.nfuncpi;
        const unsigned short *lfuncpi = s4.nfuncpi;

        // Compute this unique AO shell
        //timer_on("Computing the AO shell");
        tb_[thread]->compute_shell(si, sj, sk, sl);
        //timer_off("Computing the AO shell");

        mints_timer_on("TwoBodySOInt::compute_shell actual transform");

        for (int isym=0; isym<nirrep; ++isym) {
            unsigned short nifunc = ifuncpi[isym];
            for (int itr=0; itr<nifunc; itr++) {
                mints_timer_on("itr");

                const AOTransformFunction &ifunc = s1.soshellpi[isym][itr];
                double icoef = ifunc.coef;
                int iaofunc = ifunc.aofunc;
                int isofunc = ifunc.sofunc;
                int iaooff = iaofunc;
                int isooff = isofunc;

                mints_timer_off("itr");

                for (int jsym=0; jsym<nirrep; ++jsym) {
                    unsigned short njfunc = jfuncpi[jsym];
                    for (int jtr=0; jtr<njfunc; jtr++) {

                        mints_timer_on("jtr");

                        const AOTransformFunction &jfunc = s2.soshellpi[jsym][jtr];
                        double jcoef = jfunc.coef * icoef;
                        int jaofunc = jfunc.aofunc;
                        int jsofunc = jfunc.sofunc;
                        int jaooff = iaooff*nao2 + jaofunc;
                        int jsooff = isooff*nso2 + jsofunc;

                        mints_timer_off("jtr");

                        for (int ksym=0; ksym<nirrep; ++ksym) {
                            unsigned short nkfunc = kfuncpi[ksym];
                            for (int ktr=0; ktr<nkfunc; ktr++) {

                                mints_timer_on("ktr");

                                const AOTransformFunction &kfunc = s3.soshellpi[ksym][ktr];
                                double kcoef = kfunc.coef * jcoef;
                                int kaofunc = kfunc.aofunc;
                                int ksofunc = kfunc.sofunc;
                                int kaooff = jaooff*nao3 + kaofunc;
                                int ksooff = jsooff*nso3 + ksofunc;

                                mints_timer_off("ktr");

                                int lsym = isym ^ jsym ^ ksym;
                                unsigned short nlfunc = lfuncpi[lsym];
                                for (int ltr=0; ltr<nlfunc; ltr++) {

                                    mints_timer_on("ltr");

                                    const AOTransformFunction &lfunc = s4.soshellpi[lsym][ltr];
                                    double lcoef = lfunc.coef * kcoef;
                                    int laofunc = lfunc.aofunc;
                                    int lsofunc = lfunc.sofunc;
                                    int laooff = kaooff*nao4 + laofunc;
                                    int lsooff = ksooff*nso4 + lsofunc;

                                    mints_timer_off("ltr");
                                    mints_timer_on("transform");

                                    dprintf("lambda_T %d, lcoef %lf, aobuff %lf\n", lambda_T, lcoef, aobuff[laooff]);
                                    buffer_[thread][lsooff] += lambda_T * lcoef * aobuff[laooff];

                                    mints_timer_off("transform");
                                }
                            }
                        }
                    }
                }
            }
        }

        mints_timer_off("TwoBodySOInt::compute_shell actual transform");
    }

    mints_timer_off("TwoBodySOInt::compute_shell full shell transform");

    provide_IJKL(uish, ujsh, uksh, ulsh, body);

    mints_timer_off("TwoBodySOInt::compute_shell overall");
    //timer_off("TwoBodySOInt::compute_shell overall");
}

template<typename TwoBodySOIntFunctor>
void TwoBodySOInt::provide_IJKL(int ish, int jsh, int ksh, int lsh, TwoBodySOIntFunctor& body)
{
    int thread = 0;
    //Old call WorldComm->thread_id(pthread_self());
    #ifdef _OPENMP
      thread = omp_get_thread_num();
    #endif

    mints_timer_on("TwoBodySOInt::provide_IJKL overall");
    //timer_on("TwoBodySOInt::provide_IJKL overall");

    int nso2 = b2_->nfunction(jsh);
    int nso3 = b3_->nfunction(ksh);
    int nso4 = b4_->nfunction(lsh);

    int n1 = b1_->nfunction(ish);
    int n2 = b2_->nfunction(jsh);
    int n3 = b3_->nfunction(ksh);
    int n4 = b4_->nfunction(lsh);

    int itr;
    int jtr;
    int ktr;
    int ltr;

    for (itr=0; itr<n1; itr++) {

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

                    if (fabs(buffer_[thread][lsooff]) > cutoff_) {
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

                        mints_timer_on("TwoBodySOInt::provide_IJKL functor");

                        // func off/on
                        body(iiabs, jjabs, kkabs, llabs,
                             iiirrep, iirel,
                             jjirrep, jjrel,
                             kkirrep, kkrel,
                             llirrep, llrel,
                             buffer_[thread][lsooff]);

                        mints_timer_off("TwoBodySOInt::provide_IJKL functor");
                    }
                }
            }
        }
    }
    //timer_off("TwoBodySOInt::provide_IJKL overall");
    mints_timer_off("TwoBodySOInt::provide_IJKL overall");
}

template<typename TwoBodySOIntFunctor>
void TwoBodySOInt::compute_shell_deriv1(int uish, int ujsh, int uksh, int ulsh, TwoBodySOIntFunctor& body)
{
    int thread = 0;
    //Old call: WorldComm->thread_id(pthread_self());

    const double *aobuffer = tb_[thread]->buffer();

    const SOTransform &t1 = b1_->sotrans(uish);
    const SOTransform &t2 = b2_->sotrans(ujsh);
    const SOTransform &t3 = b3_->sotrans(uksh);
    const SOTransform &t4 = b4_->sotrans(ulsh);

    const int nso1 = b1_->nfunction(uish);
    const int nso2 = b2_->nfunction(ujsh);
    const int nso3 = b3_->nfunction(uksh);
    const int nso4 = b4_->nfunction(ulsh);
    const size_t nso = nso1*nso2*nso3*nso4;

    // These are needed, because the buffers for each perturbation are spaced by
    // quartets expressed in terms of Cartesian basis functions, regardless of puream.
    const int nao1 = b1_->naofunction(uish);
    const int nao2 = b2_->naofunction(ujsh);
    const int nao3 = b3_->naofunction(uksh);
    const int nao4 = b4_->naofunction(ulsh);
    const size_t nao = nao1*nao2*nao3*nao4;

    const int iatom = tb_[thread]->basis1()->shell(t1.aoshell[0].aoshell).ncenter();
    const int jatom = tb_[thread]->basis2()->shell(t2.aoshell[0].aoshell).ncenter();
    const int katom = tb_[thread]->basis3()->shell(t3.aoshell[0].aoshell).ncenter();
    const int latom = tb_[thread]->basis4()->shell(t4.aoshell[0].aoshell).ncenter();

    // These 3 sections are not shell specific so we can just use petite1_
    const unsigned short istablizer = petite1_->stablizer(iatom);
    const unsigned short jstablizer = petite1_->stablizer(jatom);
    const unsigned short kstablizer = petite1_->stablizer(katom);
    const unsigned short lstablizer = petite1_->stablizer(latom);

    const int istabdense = dcd_->bits_to_dense_numbering(istablizer);
    const int jstabdense = dcd_->bits_to_dense_numbering(jstablizer);
    const int kstabdense = dcd_->bits_to_dense_numbering(kstablizer);
    const int lstabdense = dcd_->bits_to_dense_numbering(lstablizer);

    const int ijstablizer = dcd_->intersection(istabdense, jstabdense);
    const int klstablizer = dcd_->intersection(kstabdense, lstabdense);
    const int ijklstablizer = dcd_->intersection(ijstablizer, klstablizer);

    const int* R_list = dcd_->dcr(istabdense, jstabdense);
    const int* S_list = dcd_->dcr(kstabdense, lstabdense);
    const int* T_list = dcd_->dcr(ijstablizer, klstablizer);

    const int R_size = R_list[0];
    const int S_size = S_list[0];
    const int T_size = T_list[0];

    int lambda_T = petite1_->nirrep() / dcd_->subgroup_dimensions(ijklstablizer);

    std::vector<int> sj_arr, sk_arr, sl_arr;

    int si = petite1_->unique_shell_map(uish, 0);
    const int siatom = tb_[thread]->basis1()->shell(si).ncenter();

    for (int ij=1; ij <= R_size; ++ij) {
        int sj = petite2_->unique_shell_map(ujsh, R_list[ij]);
        const int sjatom = tb_[thread]->basis2()->shell(sj).ncenter();

        for (int ijkl=1; ijkl <= T_size; ++ijkl) {
            int sk = petite3_->unique_shell_map(uksh, T_list[ijkl]);
            int llsh = petite4_->unique_shell_map(ulsh, T_list[ijkl]);
            const int skatom = tb_[thread]->basis3()->shell(sk).ncenter();

            for (int kl=1; kl <= S_size; ++kl) {
                int sl = petite4_->shell_map(llsh, S_list[kl]);
                const int slatom = tb_[thread]->basis4()->shell(sl).ncenter();

                // Check AM
                int total_am = tb_[thread]->basis1()->shell(si).am() +
                               tb_[thread]->basis2()->shell(sj).am() +
                               tb_[thread]->basis3()->shell(sk).am() +
                               tb_[thread]->basis4()->shell(sl).am();

                //                if(!(total_am%2)||
                //                             (BasisSet.shells[si].center!=BasisSet.shells[sj].center)||
                //                             (BasisSet.shells[sj].center!=BasisSet.shells[sk].center)||
                //                             (BasisSet.shells[sk].center!=BasisSet.shells[sl].center)) {

                //                outfile->Printf( "total_am %d siatom %d sjatom %d skatom %d slatom %d\n",
                //                        total_am, siatom, sjatom, skatom, slatom);
                if (!(total_am % 2) ||
                    (siatom != sjatom) ||
                    (sjatom != skatom) ||
                    (skatom != slatom)) {
                    //                    outfile->Printf( "\tadding\n");
                    sj_arr.push_back(sj);
                    sk_arr.push_back(sk);
                    sl_arr.push_back(sl);
                }
            }
        }
    }

    // Obtain SALC transformation objects
    // This probably won't work. I'll probably need
    // the SALC for the petite list of shells.
    // If so, this will move into the for loop below.

    // Compute integral using si, sj_arr, sk_arr, sl_arr
    // Loop over unique quartets
    const AOTransform& s1 = b1_->aotrans(si);
    const unsigned short *ifuncpi = s1.nfuncpi;
    const CdSalcWRTAtom& c1 = cdsalcs_->atom_salc(siatom);

    // Zero out SALC memory
    for (size_t i=0; i<cdsalcs_->ncd(); ++i)
        ::memset(deriv_[thread][i], 0, sizeof(double)*nso);


    double pfac = 1.0;
    //    if (uish == ujsh)
    //      pfac *= 0.5;
    //    if (uksh == ulsh)
    //      pfac *= 0.5;
    //    if (uish == uksh && ujsh == ulsh || uish == ulsh && ujsh == uksh)
    //      pfac *= 0.5;

    for (size_t n=0; n<sj_arr.size(); ++n) {
        int sj = sj_arr[n];
        int sk = sk_arr[n];
        int sl = sl_arr[n];


        const AOTransform& s2 = b2_->aotrans(sj);
        const AOTransform& s3 = b3_->aotrans(sk);
        const AOTransform& s4 = b4_->aotrans(sl);

        const unsigned short *jfuncpi = s2.nfuncpi;
        const unsigned short *kfuncpi = s3.nfuncpi;
        const unsigned short *lfuncpi = s4.nfuncpi;

        const CdSalcWRTAtom& c2 = cdsalcs_->atom_salc(tb_[thread]->basis2()->shell(sj).ncenter());
        const CdSalcWRTAtom& c3 = cdsalcs_->atom_salc(tb_[thread]->basis3()->shell(sk).ncenter());
        const CdSalcWRTAtom& c4 = cdsalcs_->atom_salc(tb_[thread]->basis4()->shell(sl).ncenter());

        // Compute this unique AO shell
        tb_[thread]->compute_shell_deriv1(si, sj, sk, sl);

        int nirrep = b1_->nirrep();

        for (int isym=0; isym<nirrep; ++isym) {
            unsigned short nifunc = ifuncpi[isym];
            for (int itr=0; itr<nifunc; itr++) {
                mints_timer_on("itr");

                const AOTransformFunction &ifunc = s1.soshellpi[isym][itr];
                double icoef = ifunc.coef;
                int iaofunc = ifunc.aofunc;
                int isofunc = ifunc.sofunc;
                int iaooff = iaofunc;
                int isooff = isofunc;

                mints_timer_off("itr");

                for (int jsym=0; jsym<nirrep; ++jsym) {
                    unsigned short njfunc = jfuncpi[jsym];
                    for (int jtr=0; jtr<njfunc; jtr++) {

                        mints_timer_on("jtr");

                        const AOTransformFunction &jfunc = s2.soshellpi[jsym][jtr];
                        double jcoef = jfunc.coef * icoef;
                        int jaofunc = jfunc.aofunc;
                        int jsofunc = jfunc.sofunc;
                        int jaooff = iaooff*nao2 + jaofunc;
                        int jsooff = isooff*nso2 + jsofunc;

                        mints_timer_off("jtr");

                        for (int ksym=0; ksym<nirrep; ++ksym) {
                            unsigned short nkfunc = kfuncpi[ksym];
                            for (int ktr=0; ktr<nkfunc; ktr++) {

                                mints_timer_on("ktr");

                                const AOTransformFunction &kfunc = s3.soshellpi[ksym][ktr];
                                double kcoef = kfunc.coef * jcoef;
                                int kaofunc = kfunc.aofunc;
                                int ksofunc = kfunc.sofunc;
                                int kaooff = jaooff*nao3 + kaofunc;
                                int ksooff = jsooff*nso3 + ksofunc;

                                mints_timer_off("ktr");

                                int lsym = isym ^ jsym ^ ksym;
                                unsigned short nlfunc = lfuncpi[lsym];
                                for (int ltr=0; ltr<nlfunc; ltr++) {

                                    mints_timer_on("ltr");

                                    const AOTransformFunction &lfunc = s4.soshellpi[lsym][ltr];
                                    double lcoef = lfunc.coef * kcoef;
                                    int laofunc = lfunc.aofunc;
                                    int lsofunc = lfunc.sofunc;
                                    int laooff = kaooff*nao4 + laofunc;
                                    int lsooff = ksooff*nso4 + lsofunc;

                                    int total_symmetry = ifunc.irrep ^ jfunc.irrep ^ kfunc.irrep ^ lfunc.irrep;

                                    // If we're only interested in totally symmetric derivatives, skip all others.
                                    if (only_totally_symmetric_ == true && total_symmetry != 0)
                                        continue;

                                    // OK, the integral at 12 aobuff[laooff] needs
                                    // to have lambda_T * lcoef * salcCoef applied to
                                    // it and saved to salc index
                                    // Special case is center B which must be determined
                                    // from centers A, C, and D.
                                    double A[3], B[3], C[3], D[3];

                                    A[0] = aobuffer[0*nao+laooff]; A[1] = aobuffer[1*nao+laooff]; A[2] = aobuffer[2*nao+laooff];
                                    C[0] = aobuffer[3*nao+laooff]; C[1] = aobuffer[4*nao+laooff]; C[2] = aobuffer[5*nao+laooff];
                                    D[0] = aobuffer[6*nao+laooff]; D[1] = aobuffer[7*nao+laooff]; D[2] = aobuffer[8*nao+laooff];

                                    // Use translational invariance to determine B
                                    B[0] = -(A[0] + C[0] + D[0]);  B[1] = -(A[1] + C[1] + D[1]);  B[2] = -(A[2] + C[2] + D[2]);

                                    A[0] *= lambda_T * pfac * lcoef;
                                    A[1] *= lambda_T * pfac * lcoef;
                                    A[2] *= lambda_T * pfac * lcoef;
                                    B[0] *= lambda_T * pfac * lcoef;
                                    B[1] *= lambda_T * pfac * lcoef;
                                    B[2] *= lambda_T * pfac * lcoef;
                                    C[0] *= lambda_T * pfac * lcoef;
                                    C[1] *= lambda_T * pfac * lcoef;
                                    C[2] *= lambda_T * pfac * lcoef;
                                    D[0] *= lambda_T * pfac * lcoef;
                                    D[1] *= lambda_T * pfac * lcoef;
                                    D[2] *= lambda_T * pfac * lcoef;

                                    dprintf("so' derivatives: A[0] %+lf A[1] %+lf A[2] %+lf\n"
                                            "                 B[0] %+lf B[1] %+lf B[2] %+lf\n"
                                            "                 C[0] %+lf C[1] %+lf C[2] %+lf\n"
                                            "                 D[0] %+lf D[1] %+lf D[2] %+lf\n"
                                            "lsooff: %d\n"
                                            "iirrep: %d jirrep: %d kirrep: %d lirrep: %d combined: %d\n",
                                            A[0], A[1], A[2], B[0], B[1], B[2], C[0], C[1], C[2], D[0], D[1], D[2], lsooff,
                                            ifunc.irrep, jfunc.irrep, kfunc.irrep, lfunc.irrep,
                                            ifunc.irrep ^ jfunc.irrep ^ kfunc.irrep ^ lfunc.irrep);

                                    // For each center apply the so transform and salc coef.
                                    // Ax
                                    for (int nx=0; nx<c1.nx(); ++nx) {
                                        const CdSalcWRTAtom::Component element = c1.x(nx);
                                        double temp = element.coef * A[0];
                                        dprintf("Ax SALC#%d pfac %lf, A[0] %lf, contr %lf\n", element.salc, element.coef, A[0], temp);
                                        if (total_symmetry == element.irrep && fabs(temp) > cutoff_)
                                            deriv_[thread][element.salc][lsooff] += temp;
                                        dprintf(" val: %lf\n", deriv_[thread][element.salc][lsooff]);
                                    }

                                    // Ay
                                    for (int ny=0; ny<c1.ny(); ++ny) {
                                        const CdSalcWRTAtom::Component element = c1.y(ny);
                                        double temp = element.coef * A[1];
                                        dprintf("Ay SALC#%d pfac %lf, A[1] %lf, contr %lf\n", element.salc, element.coef, A[1], temp);
                                        if (total_symmetry == element.irrep && fabs(temp) > cutoff_)
                                            deriv_[thread][element.salc][lsooff] += temp;
                                        dprintf(" val: %lf\n", deriv_[thread][element.salc][lsooff]);
                                    }

                                    // Az
                                    for (int nz=0; nz<c1.nz(); ++nz) {
                                        const CdSalcWRTAtom::Component element = c1.z(nz);
                                        double temp = element.coef * A[2];
                                        dprintf("Az SALC#%d pfac %lf, A[2] %lf, contr %lf\n", element.salc, element.coef, A[2], temp);
                                        if (total_symmetry == element.irrep && fabs(temp) > cutoff_)
                                            deriv_[thread][element.salc][lsooff] += temp;
                                        dprintf(" val: %lf\n", deriv_[thread][element.salc][lsooff]);
                                    }

                                    // Bx
                                    for (int nx=0; nx<c2.nx(); ++nx) {
                                        const CdSalcWRTAtom::Component element = c2.x(nx);
                                        double temp = element.coef * B[0];
                                        dprintf("Bx SALC#%d pfac %lf, B[0] %lf, contr %lf\n", element.salc, element.coef, B[0], temp);
                                        if (total_symmetry == element.irrep && fabs(temp) > cutoff_)
                                            deriv_[thread][element.salc][lsooff] += temp;
                                        dprintf(" val: %lf\n", deriv_[thread][element.salc][lsooff]);
                                    }

                                    // By
                                    for (int ny=0; ny<c2.ny(); ++ny) {
                                        const CdSalcWRTAtom::Component element = c2.y(ny);
                                        double temp = element.coef * B[1];
                                        dprintf("By SALC#%d pfac %lf, B[1] %lf, contr %lf\n", element.salc, element.coef, B[1], temp);
                                        if (total_symmetry == element.irrep && fabs(temp) > cutoff_)
                                            deriv_[thread][element.salc][lsooff] += temp;
                                        dprintf(" val: %lf\n", deriv_[thread][element.salc][lsooff]);
                                    }

                                    // Bz
                                    for (int nz=0; nz<c2.nz(); ++nz) {
                                        const CdSalcWRTAtom::Component element = c2.z(nz);
                                        double temp = element.coef * B[2];
                                        dprintf("Bz SALC#%d pfac %lf, B[2] %lf, contr %lf\n", element.salc, element.coef, B[2], temp);
                                        if (total_symmetry == element.irrep && fabs(temp) > cutoff_)
                                            deriv_[thread][element.salc][lsooff] += temp;
                                        dprintf(" val: %lf\n", deriv_[thread][element.salc][lsooff]);
                                    }

                                    // Cx
                                    for (int nx=0; nx<c3.nx(); ++nx) {
                                        const CdSalcWRTAtom::Component element = c3.x(nx);
                                        double temp = element.coef * C[0];
                                        dprintf("Cx SALC#%d pfac %lf, C[0] %lf, contr %lf\n", element.salc, element.coef, C[0], temp);
                                        if (total_symmetry == element.irrep && fabs(temp) > cutoff_)
                                            deriv_[thread][element.salc][lsooff] += temp;
                                        dprintf(" val: %lf\n", deriv_[thread][element.salc][lsooff]);
                                    }

                                    // Cy
                                    for (int ny=0; ny<c3.ny(); ++ny) {
                                        const CdSalcWRTAtom::Component element = c3.y(ny);
                                        double temp = element.coef * C[1];
                                        dprintf("Cy SALC#%d pfac %lf, C[1] %lf, contr %lf\n", element.salc, element.coef, C[1], temp);
                                        if (total_symmetry == element.irrep && fabs(temp) > cutoff_)
                                            deriv_[thread][element.salc][lsooff] += temp;
                                        dprintf(" val: %lf\n", deriv_[thread][element.salc][lsooff]);
                                    }

                                    // Cz
                                    for (int nz=0; nz<c3.nz(); ++nz) {
                                        const CdSalcWRTAtom::Component element = c3.z(nz);
                                        double temp = element.coef * C[2];
                                        dprintf("Cz SALC#%d pfac %lf, C[2] %lf, contr %lf\n", element.salc, element.coef, C[2], temp);
                                        if (total_symmetry == element.irrep && fabs(temp) > cutoff_)
                                            deriv_[thread][element.salc][lsooff] += temp;
                                        dprintf(" val: %lf\n", deriv_[thread][element.salc][lsooff]);
                                    }

                                    // Dx
                                    for (int nx=0; nx<c4.nx(); ++nx) {
                                        const CdSalcWRTAtom::Component element = c4.x(nx);
                                        double temp = element.coef * D[0];
                                        dprintf("Dx SALC#%d pfac %lf, D[0] %lf, contr %lf\n", element.salc, element.coef, D[0], temp);
                                        if (total_symmetry == element.irrep && fabs(temp) > cutoff_)
                                            deriv_[thread][element.salc][lsooff] += temp;
                                        dprintf(" val: %lf\n", deriv_[thread][element.salc][lsooff]);
                                    }

                                    // Dy
                                    for (int ny=0; ny<c4.ny(); ++ny) {
                                        const CdSalcWRTAtom::Component element = c4.y(ny);
                                        double temp = element.coef * D[1];
                                        dprintf("Dy SALC#%d pfac %lf, D[1] %lf, contr %lf\n", element.salc, element.coef, D[1], temp);
                                        if (total_symmetry == element.irrep && fabs(temp) > cutoff_)
                                            deriv_[thread][element.salc][lsooff] += temp;
                                        dprintf(" val: %lf\n", deriv_[thread][element.salc][lsooff]);
                                    }

                                    // Dz
                                    for (int nz=0; nz<c4.nz(); ++nz) {
                                        const CdSalcWRTAtom::Component element = c4.z(nz);
                                        double temp = element.coef * D[2];
                                        dprintf("Dz SALC#%d pfac %lf, D[2] %lf, contr %lf\n", element.salc, element.coef, D[2], temp);
                                        if (total_symmetry == element.irrep && fabs(temp) > cutoff_)
                                            deriv_[thread][element.salc][lsooff] += temp;
                                        dprintf(" val: %lf\n", deriv_[thread][element.salc][lsooff]);
                                    }

                                    dprintf("\n");
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    provide_IJKL_deriv1(uish, ujsh, uksh, ulsh, body);

} // function

template<typename TwoBodySOIntFunctor>
void TwoBodySOInt::provide_IJKL_deriv1(int ish, int jsh, int ksh, int lsh, TwoBodySOIntFunctor& body)
{
    int thread = 0;
    //Old call: WorldComm->thread_id(pthread_self());

    mints_timer_on("TwoBodySOInt::provide_IJKL overall");

    int nso2 = b2_->nfunction(jsh);
    int nso3 = b3_->nfunction(ksh);
    int nso4 = b4_->nfunction(lsh);

    int n1 = b1_->nfunction(ish);
    int n2 = b2_->nfunction(jsh);
    int n3 = b3_->nfunction(ksh);
    int n4 = b4_->nfunction(lsh);

    int itr;
    int jtr;
    int ktr;
    int ltr;

    for (itr=0; itr<n1; itr++) {

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

                    // Only totally symmetric pertubations are considered here!
                    if(isym^jsym^ksym^lsym) continue;

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

                    mints_timer_on("TwoBodySOInt::provide_IJKL functor");
                    for (size_t i=0; i<cdsalcs_->ncd(); ++i) {
                        if (fabs(deriv_[thread][i][lsooff]) > cutoff_)
                            body(i, iiabs, jjabs, kkabs, llabs,
                                 iiirrep, iirel,
                                 jjirrep, jjrel,
                                 kkirrep, kkrel,
                                 llirrep, llrel,
                                 deriv_[thread][i][lsooff]);
                    }
                    body.next_tpdm_element();

                    mints_timer_off("TwoBodySOInt::provide_IJKL functor");
                }
            }
        }
    }
    mints_timer_off("TwoBodySOInt::provide_IJKL overall");
}

template<typename TwoBodySOIntFunctor>
void TwoBodySOInt::compute_integrals(TwoBodySOIntFunctor &functor)
{
    if (comm_ == "MADNESS") {
        throw PSIEXCEPTION("PSI4 was not built with MADNESS. "
                           "Please rebuild PSI4 with MADNESS, or "
                           "change your COMMUNICATOR "
                           "environment variable to MPI or LOCAL.\n");
    }
    else {
        boost::shared_ptr<SOShellCombinationsIterator>
        shellIter(new SOShellCombinationsIterator(b1_, b2_, b3_, b4_));
        this->compute_quartets(shellIter, functor);
    }
}

template<typename TwoBodySOIntFunctor>
void TwoBodySOInt::compute_integrals_deriv1(TwoBodySOIntFunctor &functor)
{
    if(!only_totally_symmetric_)
        throw PSIEXCEPTION("The way the TPDM is stored and iterated enables only totally symmetric"
                           " perturbations to be considered right now!");

    if (comm_ == "MADNESS") {  }
    else {
        boost::shared_ptr<SO_PQ_Iterator> PQIter(new SO_PQ_Iterator(b1_));
        size_t pair_number = 0;
        for (PQIter->first(); PQIter->is_done() == false; PQIter->next()) {
            compute_pq_pair_deriv1<TwoBodySOIntFunctor>(
                        PQIter->p(), PQIter->q(), pair_number, functor);
            pair_number++;
        }
    }
}

typedef boost::shared_ptr<OneBodySOInt> SharedOneBodySOInt;
typedef boost::shared_ptr<TwoBodySOInt> SharedTwoBodySOInt;

}

#endif // _psi_src_lib_libmints_sointegral_h_

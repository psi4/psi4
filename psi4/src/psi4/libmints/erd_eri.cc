/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
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


#ifdef USING_erd

#include "psi4/psi4-dec.h"
#include "erd_eri.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/wavefunction.h"

#define OLDCODE 0
#define DEBUG 0

// Name mangling
#include <ERD/ERD_MANGLE.h>
#define C_ERD__GENER_ERI_BATCH  ERD_MANGLE_GLOBAL_(erd__gener_eri_batch,  ERD__GENER_ERI_BATCH)
#define C_ERD__MEMORY_CSGTO     ERD_MANGLE_GLOBAL_(erd__memory_csgto,     ERD__MEMORY_CSGTO)
#define C_ERD__MEMORY_ERI_BATCH ERD_MANGLE_GLOBAL_(erd__memory_eri_batch, ERD__MEMORY_ERI_BATCH)


extern "C" {
void C_ERD__GENER_ERI_BATCH(const F_INT &imax, const F_INT &zmax, const F_INT &nalpha, const F_INT &ncoeff,
                            const F_INT &ncsum, const F_INT &ncgto1, const F_INT &ncgto2,
                            const F_INT &ncgto3, const F_INT &ncgto4, const F_INT &npgto1,
                            const F_INT &npgto2, const F_INT &npgto3, const F_INT &npgto4,
                            const F_INT &shell1, const F_INT &shell2, const F_INT &shell3,
                            const F_INT &shell4, const double &x1, const double &y1, const double &z1,
                            const double &x2,const double &y2,const double &z2, const double &x3,
                            const double &y3,const double &z3,const double &x4, const double &y4, const double &z4,
                            const double *alpha, const double *cc, const F_INT *ccbeg, const F_INT *ccend,
                            const F_BOOL &spheric,  const F_BOOL &screen, F_INT *icore,
                            F_INT &nbatch, F_INT & nfirst, double *zcore );
void C_ERD__MEMORY_ERI_BATCH(const F_INT &nalpha, const F_INT &ncoeff,
                             const F_INT &ncgto1, const F_INT &ncgto2, const F_INT &ncgto3, const F_INT &ncgto4,
                             const F_INT &npgto1, const F_INT &npgto2, const F_INT &npgto3, const F_INT &npgto4,
                             const F_INT &shell1, const F_INT &shell2, const F_INT &shell3, const F_INT &shell4,
                             const double &x1, const double &y1, const double &z1, const double &x2, const double &y2,
                             const double &z2, const double &x3, const double &y3, const double &z3, const double &x4,
                             const double &y4, const double &z4, const double *alpha, const double *cc, const F_BOOL &spheric,
                             F_INT &imin, F_INT &iopt, F_INT &zmin, F_INT &zopt);
}

namespace psi{

ERDTwoElectronInt::ERDTwoElectronInt(const IntegralFactory* integral, int deriv, bool use_shell_pairs):
    TwoBodyAOInt(integral, deriv),
    d_buffer_size_(0L),
    i_buffer_size_(0L),
    buffer_offset_(0L)
{
    bs1_ = original_bs1_;
    bs2_ = original_bs2_;
    bs3_ = original_bs3_;
    bs4_ = original_bs4_;
    same_bs_ = (bs1_ == bs2_ && bs1_ == bs3_ && bs1_ == bs4_);

    has_puream_ = bs1_->has_puream() || bs2_->has_puream() ||
                  bs3_->has_puream() || bs4_->has_puream();

    spheric_ = 0; //bs1_->has_puream() && bs2_->has_puream() &&
                  //bs3_->has_puream() && bs4_->has_puream();

    size_t max_cart = INT_NCART(basis1()->max_am()) * INT_NCART(basis2()->max_am()) *
                      INT_NCART(basis3()->max_am()) * INT_NCART(basis4()->max_am());

    try {
        tformbuf_ = new double[max_cart];
    }
    catch (std::bad_alloc& e) {
        outfile->Printf( "Error allocating tformbuf_.\n%s\n", e.what());
        exit(EXIT_FAILURE);
    }
    memset(tformbuf_, 0, sizeof(double)*max_cart);


    try {
        target_ = new double[max_cart];
    }
    catch (std::bad_alloc& e) {
        outfile->Printf( "Error allocating target_.\n%s\n", e.what());
        exit(EXIT_FAILURE);
    }
    memset(target_, 0, sizeof(double)*max_cart);

    size_t max_nprim = basis1()->max_nprimitive() +
                       basis2()->max_nprimitive() +
                       basis3()->max_nprimitive() +
                       basis4()->max_nprimitive();
    cc_ = new double[max_nprim];
    alpha_ = new double[max_nprim];

#if OLDCODE
    // Basis set dependent info
    new_cc_1_ = new double[original_bs1_->nprimitive()];
    alpha_1_ = new double[original_bs1_->nprimitive()];
    npgto_1_ = new int[original_bs1_->nshell()];
    xyz_1_ = new double[3*original_bs1_->nshell()];
    am_1_ = new int[original_bs1_->nshell()];
    pgto_offsets_1_ = new int[original_bs1_->nshell()];
    if(same_bs_){
        new_cc_2_ = new_cc_1_;
        new_cc_3_ = new_cc_1_;
        new_cc_4_ = new_cc_1_;
        alpha_2_ = alpha_1_;
        alpha_3_ = alpha_1_;
        alpha_4_ = alpha_1_;
        xyz_2_ = xyz_1_;
        xyz_3_ = xyz_1_;
        xyz_4_ = xyz_1_;
        npgto_2_ = npgto_1_;
        npgto_3_ = npgto_1_;
        npgto_4_ = npgto_1_;
        am_2_ = am_1_;
        am_3_ = am_1_;
        am_4_ = am_1_;
        pgto_offsets_2_ = pgto_offsets_1_;
        pgto_offsets_3_ = pgto_offsets_1_;
        pgto_offsets_4_ = pgto_offsets_1_;
    }else{
        new_cc_2_ = new double[original_bs2_->nprimitive()];
        new_cc_3_ = new double[original_bs3_->nprimitive()];
        new_cc_4_ = new double[original_bs4_->nprimitive()];
        alpha_2_ = new double[original_bs2_->nprimitive()];
        alpha_3_ = new double[original_bs3_->nprimitive()];
        alpha_4_ = new double[original_bs4_->nprimitive()];
        xyz_2_ = new double[3*original_bs2_->nshell()];
        xyz_3_ = new double[3*original_bs3_->nshell()];
        xyz_4_ = new double[3*original_bs4_->nshell()];
        npgto_2_ = new int[original_bs2_->nshell()];
        npgto_3_ = new int[original_bs3_->nshell()];
        npgto_4_ = new int[original_bs4_->nshell()];
        am_2_ = new int[original_bs2_->nshell()];
        am_3_ = new int[original_bs3_->nshell()];
        am_4_ = new int[original_bs4_->nshell()];
        pgto_offsets_2_ = new int[original_bs2_->nshell()];
        pgto_offsets_3_ = new int[original_bs3_->nshell()];
        pgto_offsets_4_ = new int[original_bs4_->nshell()];
    }
#endif

    screen_ = 0;

    // Ask ERD for the maximum amount of scratch it'll need
    compute_scratch_size();
    try {
        dscratch_ = new double[d_buffer_size_];
    }
    catch (std::bad_alloc& e) {
        outfile->Printf( "Error allocating dscratch_.\n%s\n", e.what());
        exit(EXIT_FAILURE);
    }
    memset(dscratch_, 0, sizeof(double)*d_buffer_size_);

    try {
        iscratch_ = new F_INT[i_buffer_size_];
    }
    catch (std::bad_alloc& e) {
        outfile->Printf( "Error allocating iscratch.\n%s\n", e.what());
        exit(EXIT_FAILURE);
    }

    normalize_basis();
}


ERDTwoElectronInt::~ERDTwoElectronInt()
{
    delete[] alpha_;
    delete[] cc_;
#if OLDCODE
    delete[] alpha_1_;
    delete[] am_1_;
    delete[] npgto_1_;
    delete[] xyz_1_;
    delete[] pgto_offsets_1_;
    delete[] new_cc_1_;
    if(!same_bs_){
        delete[] alpha_2_;
        delete[] alpha_3_;
        delete[] alpha_4_;
        delete[] am_2_;
        delete[] am_3_;
        delete[] am_4_;
        delete[] npgto_2_;
        delete[] npgto_3_;
        delete[] npgto_4_;
        delete[] xyz_2_;
        delete[] xyz_3_;
        delete[] xyz_4_;
        delete[] pgto_offsets_2_;
        delete[] pgto_offsets_3_;
        delete[] pgto_offsets_4_;
        delete[] new_cc_2_;
        delete[] new_cc_3_;
        delete[] new_cc_4_;
    }
#endif
    delete[] tformbuf_;
    delete[] target_;
    delete[] dscratch_;
    delete[] iscratch_;
}


size_t ERDTwoElectronInt::compute_shell(const psi::AOShellCombinationsIterator& shellIter)
{
    return compute_shell(shellIter.p(), shellIter.q(), shellIter.r(), shellIter.s());
}


void ERDTwoElectronInt::compute_scratch_size()
{
    int shell1;
    int npgto1 = 0;
    for(int shell = 0; shell < original_bs1_->nshell(); ++shell)
        if(original_bs1_->shell(shell).nprimitive() > npgto1){
            npgto1 = original_bs1_->shell(shell).nprimitive();
            shell1 = shell;
        }
    int shell2;
    int npgto2 = 0;
    for(int shell = 0; shell < original_bs2_->nshell(); ++shell)
        if(original_bs2_->shell(shell).nprimitive() > npgto2){
            npgto2 = original_bs2_->shell(shell).nprimitive();
            shell2 = shell;
        }
    int shell3;
    int npgto3 = 0;
    for(int shell = 0; shell < original_bs3_->nshell(); ++shell)
        if(original_bs3_->shell(shell).nprimitive() > npgto3){
            npgto3 = original_bs3_->shell(shell).nprimitive();
            shell3 = shell;
        }
    int shell4;
    int npgto4 = 0;
    for(int shell = 0; shell < original_bs4_->nshell(); ++shell)
        if(original_bs4_->shell(shell).nprimitive() > npgto4){
            npgto4 = original_bs4_->shell(shell).nprimitive();
            shell4 = shell;
        }
    const GaussianShell &gs1 = original_bs1_->shell(shell1);
    const GaussianShell &gs2 = original_bs2_->shell(shell2);
    const GaussianShell &gs3 = original_bs3_->shell(shell3);
    const GaussianShell &gs4 = original_bs4_->shell(shell4);
    double x1 = 1.0;
    double y1 = 1.0;
    double z1 = 1.0;
    double x2 = 2.0;
    double y2 = 2.0;
    double z2 = 2.0;
    double x3 = 3.0;
    double y3 = 3.0;
    double z3 = 3.0;
    double x4 = 4.0;
    double y4 = 4.0;
    double z4 = 4.0;
    F_INT ncgto1 = 1;
    F_INT ncgto2 = 1;
    F_INT ncgto3 = 1;
    F_INT ncgto4 = 1;
    F_INT ncgto = 4;
    F_INT npgto = npgto1 + npgto2 + npgto3 + npgto4;
    F_INT am1 = original_bs1_->max_am();
    F_INT am2 = original_bs2_->max_am();
    F_INT am3 = original_bs3_->max_am();
    F_INT am4 = original_bs4_->max_am();
    F_INT imin = 0;
    F_INT iopt = 0;
    F_INT zmin = 0;
    F_INT zopt = 0;
    F_INT last_pgto = 0;
    for(int pgto1 = 0; pgto1 < npgto1; ++pgto1){
        cc_[last_pgto] = gs1.coef(pgto1);
        alpha_[last_pgto] = gs1.exp(pgto1);
        ++last_pgto;
    }
    for(int pgto2 = 0; pgto2 < npgto2; ++pgto2){
        cc_[last_pgto] = gs2.coef(pgto2);
        alpha_[last_pgto] = gs2.exp(pgto2);
        ++last_pgto;
    }
    for(int pgto3 = 0; pgto3 < npgto3; ++pgto3){
        cc_[last_pgto] = gs3.coef(pgto3);
        alpha_[last_pgto] = gs3.exp(pgto3);
        ++last_pgto;
    }
    for(int pgto4 = 0; pgto4 < npgto4; ++pgto4){
        cc_[last_pgto] = gs4.coef(pgto4);
        alpha_[last_pgto] = gs4.exp(pgto4);
        ++last_pgto;
    }
    long int nbatch = 0;
    // Compute the amount of memory needed for the largest quartet
    C_ERD__MEMORY_ERI_BATCH(npgto, npgto, ncgto1, ncgto2, ncgto3, ncgto4,
                            npgto1, npgto2, npgto3, npgto4, am1, am2, am3, am4,
                            x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4,
                            alpha_, cc_, spheric_, imin, iopt, zmin, zopt);
#if DEBUG
    outfile->Printf( "\timin %d iopt %d zmin %d zopt %d\n", imin, iopt, zmin, zopt);
#endif
    i_buffer_size_ = iopt;
    d_buffer_size_ = zopt;

    // We always treat basis sets as segmented, right now.
    ccbeg_[0] = 1;
    ccbeg_[1] = 1;
    ccbeg_[2] = 1;
    ccbeg_[3] = 1;
}

/**
 * Applies the normalization needed for ERD, with psi4.  Note that all calls to
 * erd__normalize_cartesian must be eliminated in the ERD code for this to work.
 * While we're at it, we populate the flattened data arrays containing basis info.
 */
void ERDTwoElectronInt::normalize_basis()
{
#if OLDCODE
    // Basis set 1
    int count = 0;
    for(int shell = 0; shell < original_bs1_->nshell(); ++shell){
        pgto_offsets_1_[shell] = count;
        const GaussianShell &gs = original_bs1_->shell(shell);
        int L = gs.am();
        double sum = 0.0;
        for(int j = 0; j < gs.nprimitive(); j++){
            for(int k = 0; k <= j; k++){
                double a1 = gs.exp(j);
                double a2 = gs.exp(k);
                double temp = (gs.original_coef(j) * gs.original_coef(k));
                double temp2 = ((double) L + 1.5);
                double temp3 = (2.0 * sqrt(a1 * a2) / (a1 + a2));
                temp3 = pow(temp3, temp2);
                temp = temp * temp3;
                sum = sum + temp;
                if(j != k)
                    sum = sum + temp;
            }
        }
        double prefac = 1.0;
        if(L > 1)
            prefac =pow(2.0, 2*L) / df[2*L];
        double norm = sqrt(prefac / sum);
        for(int j = 0; j < gs.nprimitive(); j++){
            alpha_1_[count] = gs.exp(j);
            new_cc_1_[count++] = gs.original_coef(j) * norm;
        }
        Vector3 xyz = gs.center();
        xyz_1_[3*shell+0] = xyz[0];
        xyz_1_[3*shell+1] = xyz[1];
        xyz_1_[3*shell+2] = xyz[2];
        npgto_1_[shell] = gs.nprimitive();
        am_1_[shell] = gs.am();
    }
    if(!same_bs_) return;

    // Basis set 2
    count = 0;
    for(int shell = 0; shell < original_bs2_->nshell(); ++shell){
        pgto_offsets_2_[shell] = count;
        const GaussianShell &gs = original_bs2_->shell(shell);
        int L = gs.am();
        double sum = 0.0;
        for(int j = 0; j < gs.nprimitive(); j++){
            for(int k = 0; k <= j; k++){
                double a1 = gs.exp(j);
                double a2 = gs.exp(k);
                double temp = (gs.original_coef(j) * gs.original_coef(k));
                double temp2 = ((double) L + 1.5);
                double temp3 = (2.0 * sqrt(a1 * a2) / (a1 + a2));
                temp3 = pow(temp3, temp2);
                temp = temp * temp3;
                sum = sum + temp;
                if(j != k)
                    sum = sum + temp;
            }
        }
        double prefac = 1.0;
        if(L > 1)
            prefac =pow(2.0, 2*L) / df[2*L];
        double norm = sqrt(prefac / sum);
        for(int j = 0; j < gs.nprimitive(); j++){
            alpha_2_[count] = gs.exp(j);
            new_cc_2_[count++] = gs.original_coef(j) * norm;
        }
        Vector3 xyz = gs.center();
        xyz_2_[3*shell+0] = xyz[0];
        xyz_2_[3*shell+1] = xyz[1];
        xyz_2_[3*shell+2] = xyz[2];
        npgto_2_[shell] = gs.nprimitive();
        am_2_[shell] = gs.am();
    }
    // Basis set 3
    count = 0;
    for(int shell = 0; shell < original_bs3_->nshell(); ++shell){
        pgto_offsets_3_[shell] = count;
        const GaussianShell &gs = original_bs3_->shell(shell);
        int L = gs.am();
        double sum = 0.0;
        for(int j = 0; j < gs.nprimitive(); j++){
            for(int k = 0; k <= j; k++){
                double a1 = gs.exp(j);
                double a2 = gs.exp(k);
                double temp = (gs.original_coef(j) * gs.original_coef(k));
                double temp2 = ((double) L + 1.5);
                double temp3 = (2.0 * sqrt(a1 * a2) / (a1 + a2));
                temp3 = pow(temp3, temp2);
                temp = temp * temp3;
                sum = sum + temp;
                if(j != k)
                    sum = sum + temp;
            }
        }
        double prefac = 1.0;
        if(L > 1)
            prefac =pow(2.0, 2*L) / df[2*L];
        double norm = sqrt(prefac / sum);
        for(int j = 0; j < gs.nprimitive(); j++){
            alpha_3_[count] = gs.exp(j);
            new_cc_3_[count++] = gs.original_coef(j) * norm;
        }
        Vector3 xyz = gs.center();
        xyz_3_[3*shell+0] = xyz[0];
        xyz_3_[3*shell+1] = xyz[1];
        xyz_3_[3*shell+2] = xyz[2];
        npgto_3_[shell] = gs.nprimitive();
        am_3_[shell] = gs.am();
    }
    // Basis set 4
    count = 0;
    for(int shell = 0; shell < original_bs4_->nshell(); ++shell){
        pgto_offsets_4_[shell] = count;
        const GaussianShell &gs = original_bs4_->shell(shell);
        int L = gs.am();
        double sum = 0.0;
        for(int j = 0; j < gs.nprimitive(); j++){
            for(int k = 0; k <= j; k++){
                double a1 = gs.exp(j);
                double a2 = gs.exp(k);
                double temp = (gs.original_coef(j) * gs.original_coef(k));
                double temp2 = ((double) L + 1.5);
                double temp3 = (2.0 * sqrt(a1 * a2) / (a1 + a2));
                temp3 = pow(temp3, temp2);
                temp = temp * temp3;
                sum = sum + temp;
                if(j != k)
                    sum = sum + temp;
            }
        }
        double prefac = 1.0;
        if(L > 1)
            prefac =pow(2.0, 2*L) / df[2*L];
        double norm = sqrt(prefac / sum);
        for(int j = 0; j < gs.nprimitive(); j++){
            alpha_4_[count] = gs.exp(j);
            new_cc_4_[count++] = gs.original_coef(j) * norm;
        }
        Vector3 xyz = gs.center();
        xyz_4_[3*shell+0] = xyz[0];
        xyz_4_[3*shell+1] = xyz[1];
        xyz_4_[3*shell+2] = xyz[2];
        npgto_4_[shell] = gs.nprimitive();
        am_4_[shell] = gs.am();
    }

#endif
}

size_t ERDTwoElectronInt::compute_shell(int shell_i, int shell_j, int shell_k, int shell_l)
{
#if OLDCODE
    double *xyzptr1 = &xyz_1_[3*shell_i];
    double x1 = *(xyzptr1++);
    double y1 = *(xyzptr1++);
    double z1 = *(xyzptr1);
    double *xyzptr2 = &xyz_2_[3*shell_j];
    double x2 = *(xyzptr2++);
    double y2 = *(xyzptr2++);
    double z2 = *(xyzptr2);
    double *xyzptr3 = &xyz_3_[3*shell_k];
    double x3 = *(xyzptr3++);
    double y3 = *(xyzptr3++);
    double z3 = *(xyzptr3);
    double *xyzptr4 = &xyz_4_[3*shell_l];
    double x4 = *(xyzptr4++);
    double y4 = *(xyzptr4++);
    double z4 = *(xyzptr4);
    F_INT ncgto1 = 1;
    F_INT ncgto2 = 1;
    F_INT ncgto3 = 1;
    F_INT ncgto4 = 1;
    F_INT ncgto = 4;
    F_INT npgto1 = npgto_1_[shell_i];
    F_INT npgto2 = npgto_2_[shell_j];
    F_INT npgto3 = npgto_3_[shell_k];
    F_INT npgto4 = npgto_4_[shell_l];
    F_INT npgto = npgto1 + npgto2 + npgto3 + npgto4;
    F_INT am1 = am_1_[shell_i];
    F_INT am2 = am_2_[shell_j];
    F_INT am3 = am_3_[shell_k];
    F_INT am4 = am_4_[shell_l];
    ccend_[0] = npgto4;
    ccend_[1] = npgto3;
    ccend_[2] = npgto2;
    ccend_[3] = npgto1;
    int offset_i = 0;
    int offset_j = offset_i + npgto4;
    int offset_k = offset_j + npgto3;
    int offset_l = offset_k + npgto2;

    // Copy exponents and coefficients over
    ::memcpy(&(alpha_[offset_i]), &(alpha_4_[pgto_offsets_4_[shell_l]]), sizeof(double)*npgto4);
    ::memcpy(&(alpha_[offset_j]), &(alpha_3_[pgto_offsets_3_[shell_k]]), sizeof(double)*npgto3);
    ::memcpy(&(alpha_[offset_k]), &(alpha_2_[pgto_offsets_2_[shell_j]]), sizeof(double)*npgto2);
    ::memcpy(&(alpha_[offset_l]), &(alpha_1_[pgto_offsets_1_[shell_i]]), sizeof(double)*npgto1);
    ::memcpy(&(cc_[offset_i]), &(new_cc_4_[pgto_offsets_4_[shell_l]]), sizeof(double)*npgto4);
    ::memcpy(&(cc_[offset_j]), &(new_cc_3_[pgto_offsets_3_[shell_k]]), sizeof(double)*npgto3);
    ::memcpy(&(cc_[offset_k]), &(new_cc_2_[pgto_offsets_2_[shell_j]]), sizeof(double)*npgto2);
    ::memcpy(&(cc_[offset_l]), &(new_cc_1_[pgto_offsets_1_[shell_i]]), sizeof(double)*npgto1);
#else
    const GaussianShell &gs1 = original_bs1_->shell(shell_i);
    const GaussianShell &gs2 = original_bs2_->shell(shell_j);
    const GaussianShell &gs3 = original_bs3_->shell(shell_k);
    const GaussianShell &gs4 = original_bs4_->shell(shell_l);
    const double *xyzptr1 = gs1.center();
    double x1 = *(xyzptr1++);
    double y1 = *(xyzptr1++);
    double z1 = *(xyzptr1);
    const double *xyzptr2 = gs2.center();
    double x2 = *(xyzptr2++);
    double y2 = *(xyzptr2++);
    double z2 = *(xyzptr2);
    const double *xyzptr3 = gs3.center();
    double x3 = *(xyzptr3++);
    double y3 = *(xyzptr3++);
    double z3 = *(xyzptr3);
    const double *xyzptr4 = gs4.center();
    double x4 = *(xyzptr4++);
    double y4 = *(xyzptr4++);
    double z4 = *(xyzptr4);
    F_INT ncgto1 = 1;
    F_INT ncgto2 = 1;
    F_INT ncgto3 = 1;
    F_INT ncgto4 = 1;
    F_INT ncgto = 4;
    F_INT npgto1 = gs1.nprimitive();
    F_INT npgto2 = gs2.nprimitive();
    F_INT npgto3 = gs3.nprimitive();
    F_INT npgto4 = gs4.nprimitive();
    F_INT npgto = npgto1 + npgto2 + npgto3 + npgto4;
    F_INT am1 = gs1.am();
    F_INT am2 = gs2.am();
    F_INT am3 = gs3.am();
    F_INT am4 = gs4.am();
    ccend_[0] = npgto4;
    ccend_[1] = npgto3;
    ccend_[2] = npgto2;
    ccend_[3] = npgto1;
    int offset_i = 0;
    int offset_j = offset_i + npgto4;
    int offset_k = offset_j + npgto3;
    int offset_l = offset_k + npgto2;

    // Copy exponents and coefficients over
    ::memcpy(&(alpha_[offset_i]), gs4.exps(), sizeof(double)*npgto4);
    ::memcpy(&(alpha_[offset_j]), gs3.exps(), sizeof(double)*npgto3);
    ::memcpy(&(alpha_[offset_k]), gs2.exps(), sizeof(double)*npgto2);
    ::memcpy(&(alpha_[offset_l]), gs1.exps(), sizeof(double)*npgto1);
    ::memcpy(&(cc_[offset_i]), gs4.erd_coefs(), sizeof(double)*npgto4);
    ::memcpy(&(cc_[offset_j]), gs3.erd_coefs(), sizeof(double)*npgto3);
    ::memcpy(&(cc_[offset_k]), gs2.erd_coefs(), sizeof(double)*npgto2);
    ::memcpy(&(cc_[offset_l]), gs1.erd_coefs(), sizeof(double)*npgto1);
#endif

#if DEBUG
    outfile->Printf( "\n\nShell (%2d %2d | %2d %2d) - center (%2d %2d | %2d %2d) - angular momentum (%d %d | %d %d)\n",
                    shell_i, shell_j, shell_k, shell_l,
                    bs1_->shell(shell_i).ncenter(), bs2_->shell(shell_j).ncenter(), bs3_->shell(shell_k).ncenter(), bs4_->shell(shell_l).ncenter(),
                    am1, am2, am3, am4);
    outfile->Printf( "XYZ1: %16.10f %16.10f %16.10f\n", x1, y1, z1);
    outfile->Printf( "XYZ2: %16.10f %16.10f %16.10f\n", x2, y2, z2);
    outfile->Printf( "XYZ3: %16.10f %16.10f %16.10f\n", x3, y3, z3);
    outfile->Printf( "XYZ4: %16.10f %16.10f %16.10f\n", x4, y4, z4);
    outfile->Printf( "Indices  -> %d, %d\n", ccbeg_[0], ccend_[0]);

    outfile->Printf( "Number of primitives: %d %d %d %d\n", npgto1, npgto2, npgto3, npgto4);
    outfile->Printf( "Coefficients: ");
    for(int n = 0; n < npgto; ++n)
        outfile->Printf( "%14.10f ", cc_[n]);
    outfile->Printf( "\n");
    outfile->Printf( "Exponents:    ");
    for(int n = 0; n < npgto; ++n)
        outfile->Printf( "%14.10f ", alpha_[n]);
    outfile->Printf( "\n");
#endif

    F_INT nbatch;
    // Call ERD.  N.B. We reverse the shell ordering, because the first index is
    // the fastest running index in the buffer, which should be l for us.
    C_ERD__GENER_ERI_BATCH(i_buffer_size_, d_buffer_size_, npgto, npgto, ncgto,
                           ncgto4, ncgto3, ncgto2, ncgto1,
                           npgto4, npgto3, npgto2, npgto1,
                           am4, am3, am2, am1,
                           x4, y4, z4, x3, y3, z3, x2, y2, z2, x1, y1, z1,
                           alpha_, cc_, ccbeg_, ccend_, spheric_, screen_,
                           iscratch_, nbatch, buffer_offset_, dscratch_);

#if DEBUG
    outfile->Printf( "Buffer offset is %d\n", buffer_offset_-1);
    outfile->Printf( "%d integrals were computed\n", nbatch);
#endif

    int n1, n2, n3, n4;
    if (force_cartesian_) {
        n1 = gs1.ncartesian();
        n2 = gs2.ncartesian();
        n3 = gs3.ncartesian();
        n4 = gs4.ncartesian();
    } else {
        n1 = gs1.nfunction();
        n2 = gs2.nfunction();
        n3 = gs3.nfunction();
        n4 = gs4.nfunction();
    }
    const int n1234 = n1 * n2 * n3 * n4;

    if(nbatch == 0) {
        // The code should check the return value, and ignore the integrals in the buffer if we get here
        //TODO: LOTS OF CODE DOESN'T DO THIS
        ::memset(target_, 0, sizeof(double)*n1234);
    }
    else if(has_puream_ && !spheric_ && !force_cartesian_) {
        source_ = &(dscratch_[buffer_offset_-1]);
        pure_transform(shell_i, shell_j, shell_k, shell_l, 1);
    }
    else {
        ::memcpy(target_, &(dscratch_[buffer_offset_-1]), sizeof(double)*nbatch);
    }
    return n1234;
}


size_t ERDTwoElectronInt::compute_shell_deriv1(int, int, int, int)
{
    throw PSIEXCEPTION("Derivatives for ERD are NYI!");
}


size_t ERDTwoElectronInt::compute_shell_deriv2(int, int, int, int)
{
    throw PSIEXCEPTION("Derivatives for ERD are NYI!");
}

ERDERI::ERDERI(const IntegralFactory *integral, int deriv, bool use_shell_pairs)
    : ERDTwoElectronInt(integral, deriv, use_shell_pairs)
{
}

ERDERI::~ERDERI()
{
}

} // Namespace

#endif // USING_erd

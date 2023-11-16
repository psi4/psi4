/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#ifndef MP2F12_H
#define MP2F12_H

#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integralparameters.h"
#include "psi4/libmints/onebody.h"
#include "psi4/libmints/orbitalspace.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/liboptions/liboptions.h"

#include "einsums.hpp"

namespace psi { namespace mp2f12 {

class MP2F12 : public Wavefunction {
   public: 
    MP2F12(SharedWavefunction reference_wavefunction, Options& options);
    ~MP2F12() override;

    /* Compute the total MP2-F12/3C(FIX) Energy */
    virtual double compute_energy();

   protected:
    /* Print level */
    int print_;

    /* Number of OMP_THREADS */
    int nthreads_;

    /* Choose to compute CABS Singles correction */
    bool singles_;

    /* Choose CONV or DF F12 computation */
    std::string f12_type_;

    /* Density-fitting Basis Set (DFBS) */
    std::shared_ptr<BasisSet> DFBS_;

    /* Bool to turn on DF */
    bool use_df_ = false;

    /* Bool to read in precomputed F12 integrals */
    bool f12_restart_ = false;

    /* List of orbital spaces: Orbital Basis Set (OBS) 
       and Complimentary Auxiliary Basis Set (CABS) */
    std::vector<OrbitalSpace> bs_;

    /* Number of basis functions in OBS */
    int nobs_;

    /* Number of basis functions in CABS */
    int ncabs_;

    /* Number of basis functions in total */
    int nri_;

    /* Number of occupied orbitals */
    int nocc_;

    /* Number of virtual orbitals */
    int nvir_;

    /* Number of basis functions in DFBS */
    int naux_;

    /* Number of frozen core orbitals */
    int nfrzn_;

    /* F12 Correlation Factor, Contracted Gaussian-Type Geminal */
    std::vector<std::pair<double, double>> cgtg_;

    /* $\beta$ in F12 CGTG */
    double beta_;

    /* F12 energy */
    double E_f12_ = 0.0;

    /* CABS Singles Correction */
    double E_singles_ = 0.0;

    /* Total MP2-F12/3C(FIX) Energy */
    double E_f12_total_ = 0.0;

    void common_init();

    void print_header();

    /* Form the basis sets OBS and CABS */
    void form_basissets();

    /* Form the energy denominator */
    virtual void form_D(einsums::Tensor<double, 4> *D, einsums::Tensor<double, 2> *f);

    /* Form the CABS Singles correction $\frac{|f^{a'}_{i}}|^2}{e_{a'} - e_{i}}$ */
    virtual void form_cabs_singles(einsums::Tensor<double,2> *f);

    /* Form the F12/3C(FIX) correlation energy */
    virtual void form_f12_energy(einsums::Tensor<double,4> *V, einsums::Tensor<double,4> *X,
                                 einsums::Tensor<double,4> *C, einsums::Tensor<double,4> *B,
                                 einsums::Tensor<double,2> *f, einsums::Tensor<double,4> G_,
                                 einsums::Tensor<double,4> *D);
   
    /* Form the one-electron integrals H = T + V */
    virtual void form_oeints(einsums::Tensor<double, 2> *h);

    /* Form the convetional two-electron integrals */
    virtual void form_teints(const std::string& int_type, einsums::Tensor<double, 4> *ERI);
    
    /* Form the density-fitted two-electron integrals */
    virtual void form_df_teints(const std::string& int_type, einsums::Tensor<double, 4> *ERI,
                                einsums::Tensor<double, 3> *Metric);
    
    /* Form the Fock matrix */
    virtual void form_fock(einsums::Tensor<double, 2> *f, einsums::Tensor<double, 2> *k,
                   einsums::Tensor<double, 2> *fk, einsums::Tensor<double, 2> *h);

    /* Form the DF Fock matrix */
    virtual void form_df_fock(einsums::Tensor<double, 2> *f, einsums::Tensor<double, 2> *k, 
                      einsums::Tensor<double, 2> *fk, einsums::Tensor<double, 2> *h);

    /* Form the $V^{ij}_{kl}$ or $X^{ij}_{kl}$ tensor */
    virtual void form_V_or_X(einsums::Tensor<double, 4> *VX, einsums::Tensor<double, 4> *F,
                           einsums::Tensor<double, 4> *G_F, einsums::Tensor<double, 4> *FG_F2);
    
    /* Form the $C^{kl}_{ab}$ tensor */
    virtual void form_C(einsums::Tensor<double, 4> *C, einsums::Tensor<double, 4> *F,
                einsums::Tensor<double, 2> *f);
    
    /* Form the $B^{kl}_{mn}$ tensor */
    virtual void form_B(einsums::Tensor<double, 4> *B, einsums::Tensor<double, 4> *Uf,
                        einsums::Tensor<double, 4> *F2, einsums::Tensor<double, 4> *F,
                        einsums::Tensor<double, 2> *f, einsums::Tensor<double, 2> *fk,
                        einsums::Tensor<double, 2> *kk);

    void print_results();

    /* Returns the fixed amplitudes value */
    double t_(const int& p, const int& q, const int& r, const int& s);

    /* Form the $T^{ij}_{ij}\Tilde{V}^{ij}_{ij}$ contirbution to the energy */
    virtual std::pair<double, double> V_Tilde(einsums::Tensor<double, 2>& V_, einsums::Tensor<double, 4> *C,
                                      einsums::TensorView<double, 2>& K_ij, einsums::TensorView<double, 2>& D_ij,
                                      const int& i, const int& j);

    /* Form the $T^{ij}_{ij}\Tilde{B}^{ij}_{ij}T^{ij}_{ij}$ contirbution to the energy */
    virtual std::pair<double, double> B_Tilde(einsums::Tensor<double, 4>& B, einsums::Tensor<double, 4> *C,
                                      einsums::TensorView<double, 2>& D_ij,
                                      const int& i, const int& j);

    /* Converts the AO to MO matrices to einsum::Tensors */
    void convert_C(einsums::Tensor<double,2> *C, OrbitalSpace bs, const int& dim1, const int& dim2);

    /* Places the computed integral in the einsum::Tensor */
    virtual void set_ERI(einsums::TensorView<double, 4>& ERI_Slice, einsums::Tensor<double, 4> *Slice);
    void set_ERI(einsums::TensorView<double, 3>& ERI_Slice, einsums::Tensor<double, 3> *Slice);

    /* Computes the conventional two-body integrals */
    void two_body_ao_computer(const std::string& int_type, einsums::Tensor<double, 4> *GAO,
                              std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2,
                              std::shared_ptr<BasisSet> bs3, std::shared_ptr<BasisSet> bs4);

    /* Computes the DF three-index integrals */
    void three_index_ao_computer(const std::string& int_type, einsums::Tensor<double, 3> *Bpq,
                                 std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2);

    /* Form the integrals containing the DF metric [J_AB]^{-1}(B|PQ) */
    void form_metric_ints(einsums::Tensor<double, 3> *DF_ERI, bool is_fock);
    
    /* Form the integrals containing the explicit correlation (B|\hat{A}_{12}|PQ) */
    void form_oper_ints(const std::string& int_type, einsums::Tensor<double, 3> *DF_ERI);
    void form_oper_ints(const std::string& int_type, einsums::Tensor<double, 2> *DF_ERI);
};

class DiskMP2F12 : public MP2F12 {
   public:
    DiskMP2F12(SharedWavefunction reference_wavefunction, Options& options);
    ~DiskMP2F12() override;

    /* Compute the total MP2-F12/3C(FIX) Energy */
    double compute_energy() override;

   protected:
    /* Form the energy denominator */
    void form_D(einsums::DiskTensor<double, 4> *D, einsums::DiskTensor<double, 2> *f);

   //  /* Form the CABS Singles correction $\frac{|f^{a'}_{i}}|^2}{e_{a'} - e_{i}}$ */
    void form_cabs_singles(einsums::DiskTensor<double,2> *f);

    /* Form the F12/3C(FIX) correlation energy */
    void form_f12_energy(einsums::DiskTensor<double,4> *V, einsums::DiskTensor<double,4> *X,
                         einsums::DiskTensor<double,4> *C, einsums::DiskTensor<double,4> *B,
                         einsums::DiskTensor<double,2> *f, einsums::DiskTensor<double,4> *G,
                         einsums::DiskTensor<double,4> *D);

    /* Form the one-electron integrals H = T + V */
    void form_oeints(einsums::DiskTensor<double, 2> *h);

    /* Form the convetional two-electron integrals */
    void form_teints(const std::string& int_type, einsums::DiskTensor<double, 4> *ERI);

    /* Form the density-fitted two-electron integrals */
    void form_df_teints(const std::string& int_type, einsums::DiskTensor<double, 4> *ERI,
                        einsums::Tensor<double, 3> *Metric);

    /* Form the Fock matrix */
    void form_fock(einsums::DiskTensor<double, 2> *f, einsums::DiskTensor<double, 2> *k,
                   einsums::DiskTensor<double, 2> *fk, einsums::DiskTensor<double, 2> *h);

    /* Form the DF Fock matrix */
    void form_df_fock(einsums::DiskTensor<double, 2> *f, einsums::DiskTensor<double, 2> *k,
                      einsums::DiskTensor<double, 2> *fk, einsums::DiskTensor<double, 2> *h);

    /* Form the $V^{ij}_{kl}$ or $X^{ij}_{kl}$ tensor */
    void form_V_or_X(einsums::DiskTensor<double, 4> *VX, einsums::DiskTensor<double, 4> *F,
                     einsums::DiskTensor<double, 4> *G_F, einsums::DiskTensor<double, 4> *FG_F2);

    /* Form the $C^{kl}_{ab}$ tensor */
    void form_C(einsums::DiskTensor<double, 4> *C, einsums::DiskTensor<double, 4> *F,
                einsums::DiskTensor<double, 2> *f);

    /* Form the $B^{kl}_{mn}$ tensor */
    void form_B(einsums::DiskTensor<double, 4> *B, einsums::DiskTensor<double, 4> *Uf,
                einsums::DiskTensor<double, 4> *F2, einsums::DiskTensor<double, 4> *F,
                einsums::DiskTensor<double, 2> *f, einsums::DiskTensor<double, 2> *fk,
                einsums::DiskTensor<double, 2> *kk);


    /* Form the $T^{ij}_{ij}\Tilde{V}^{ij}_{ij}$ contirbution to the energy */
    std::pair<double, double> V_Tilde(einsums::Tensor<double, 2>& V_ij, einsums::DiskTensor<double, 4> *C,
                              einsums::DiskView<double, 2, 4>& K_ij, einsums::DiskView<double, 2, 4>& D_ij,
                              const int& i, const int& j);

    /* Form the $T^{ij}_{ij}\Tilde{B}^{ij}_{ij}T^{ij}_{ij}$ contirbution to the energy */
    std::pair<double, double> B_Tilde(einsums::Tensor<double, 4>& B_ij, einsums::DiskTensor<double, 4> *C,
                                      einsums::DiskView<double, 2, 4>& D_ij, const int& i, const int& j);

    /* Places the computed integral in the einsum::DiskTensor */
    void set_ERI(einsums::DiskView<double, 2, 4>& ERI_Slice, einsums::TensorView<double, 2>& Slice);
};

}} // end namespaces
#endif

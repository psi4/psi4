/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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

#ifndef SAPT2p_H
#define SAPT2p_H

#include "sapt2.h"

namespace psi {
namespace sapt {

class SAPT2p : public SAPT2 {
   private:
    void print_header() override;
    void print_results() override;

   protected:
    double e_disp21_;
    double e_disp22sdq_;
    double e_disp22t_;
    double e_est_disp22t_;
    double e_sapt2p_;

    void gARARxtARBS(int, const char *, const char, int, const char *, const char *, const char *, size_t, size_t, size_t, size_t,
                     size_t, size_t, size_t, const char *);

    double disp21_1(int, const char *, const char *, size_t, size_t, size_t, size_t);
    double disp21_2(int, const char *, const char *, size_t, size_t);

    double disp211();
    double disp220s(int, const char *, const char *, int, const char *, const char *, size_t, size_t, size_t);
    double disp220d_1(int, const char *, const char *, int, const char *, size_t, size_t, size_t);
    double disp220d_2(int, const char *, const char *, int, const char *, size_t, size_t, size_t, size_t, size_t, size_t, double *,
                      double *, const char);
    double disp220q_1(int, const char *, const char *, const char *, size_t, size_t);
    double disp220q_2(int, const char *, const char *, const char *, int, const char *, size_t, size_t, size_t);
    double disp220q_3(int, const char *, const char *, const char, int, const char *, size_t, size_t, size_t, size_t, size_t, size_t);
    double disp220q_4(int, const char *, const char *, const char, int, const char *, size_t, size_t, size_t, size_t, size_t, size_t);

    double disp220t(int, const char *, const char *, const char *, int, const char *, int, const char *, size_t, size_t, size_t,
                    size_t, size_t, size_t, double *, double *);

    // CCD Dispersion Values
    double e_disp2d_ccd_;
    double e_disp22s_ccd_;
    double e_disp22t_ccd_;
    double e_est_disp22t_ccd_;
    double e_sapt2p_ccd_;

    // CCD Dispersion Parameters
    bool ccd_disp_;
    int ccd_maxiter_;
    int min_ccd_vecs_;
    int max_ccd_vecs_;
    double ccd_e_conv_;
    double ccd_t_conv_;

    // Do MBPT and CCD dispersion?
    bool mbpt_disp_;

    // CCD Dispersion Methods
    void r_ccd_prep(const char *, const char *, const char *, const char *, const char *, const char *, const char *,
                    const char *, const char *, const char *, const char *, const char *, const char *, const char *,
                    const char *, const char *, const char *, int, const char *, int, const char *, double *, double *,
                    size_t, size_t, size_t, size_t, size_t, size_t);
    double r_ccd_energy(const char *, const char *, size_t, size_t, size_t, size_t);
    double r_ccd_iterate(const char *, const char *, const char *, const char *, const char *, const char *,
                         const char *, const char *, const char *, const char *, const char *, const char *, double *,
                         double *, int, int, int, int, int, int);
    double r_ccd_amplitudes(const char *, const char *, const char *, const char *, const char *, const char *,
                            const char *, const char *, const char *, const char *, const char *, double *, double *,
                            size_t, size_t, size_t, size_t, size_t, size_t);

    void s_ccd_prep(const char *, const char *, const char *, const char *, const char *, const char *, const char *,
                    double *, size_t, size_t, size_t, size_t, size_t, size_t);
    double s_ccd_iterate(const char *, const char *, const char *, const char *, const char *, const char *,
                         const char *, const char *, const char *, const char *, const char *, const char *,
                         const char *, const char *, const char *, double *, int, int, int, std::shared_ptr<Matrix>);
    double s_ccd_amplitudes(const char *, const char *, const char *, const char *, const char *, const char *,
                            const char *, const char *, const char *, const char *, const char *, const char *,
                            const char *, const char *, const char *, double *, size_t, size_t, size_t, std::shared_ptr<Matrix>);

    void ccd_prep(const char *, const char *, const char *, const char *, const char *, const char *, const char *,
                  const char *, const char *, int, const char *, const char *, const char *, double *, size_t, size_t, size_t,
                  std::shared_ptr<Matrix>, const char *);
    double ccd_energy(const char *, const char *, size_t, size_t);
    void ccd_iterate(const char *, const char *, const char *, const char *, const char *, const char *, const char *,
                     const char *, const char *, const char *, double *, int, int, int, std::shared_ptr<Matrix>);
    double ccd_amplitudes(const char *, const char *, const char *, const char *, const char *, const char *,
                          const char *, const char *, const char *, const char *, double *, size_t, size_t, size_t,
                          std::shared_ptr<Matrix>);

    void vvvv_prep(const char *, const char *, double **, int, int, std::shared_ptr<Matrix>);
    double **vvvv_ccd(const char *, const char *, const char *, size_t, size_t, std::shared_ptr<Matrix>);
    std::shared_ptr<Matrix> mo2no(int ampfile, const char *VV_opdm, int nvir, double cutoff);

    double **read_IJKL(int, char *, int, int);
    void write_IJKL(double **, int, const char *, int, int);

    // CCD (S)
    void disp_s_prep(const char *, const char *, const char *, const char *, int, const char *, const char *,
                     const char *, int, const char *, double *, size_t, size_t, size_t, size_t, size_t, size_t);

    // CCD (T)
    void natural_orbitalify_ccd();
    double disp220tccd(int, const char *, int, const char *, const char *, int, const char *, int, const char *,
                       const char *, double *, double *, size_t, size_t, size_t, size_t, size_t, size_t);

   public:
    SAPT2p(SharedWavefunction Dimer, SharedWavefunction MonomerA, SharedWavefunction MonomerB, Options &options,
           std::shared_ptr<PSIO> psio);
    ~SAPT2p() override;

    double compute_energy() override;

    void amplitudes() override;

    // PT Dispersion

    void disp21();
    void disp22sdq();
    void disp22t();

    // CCD Dispersion

    void disp2ccd();
    void disp22tccd();
};

/**
 * SAPTDIIS is a legacy helper for CCD
 **/
class SAPTDIIS {
   private:
    int filenum_;
    const char *vec_label_;
    const char *err_label_;
    int max_diis_vecs_;

    int diis_file_;
    size_t vec_length_;

    int curr_vec_;
    int num_vecs_;

    std::string get_err_label(const int64_t);
    std::string get_vec_label(const int64_t);

   protected:
    std::shared_ptr<PSIO> psio_;

   public:
    SAPTDIIS(int, const char *, const char *, size_t, int, std::shared_ptr<PSIO>);
    ~SAPTDIIS();

    void store_vectors();
    void get_new_vector();
};
}  // namespace sapt
}  // namespace psi

#endif

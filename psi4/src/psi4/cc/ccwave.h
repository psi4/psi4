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

#ifndef CCWAVE_H
#define CCWAVE_H

#include <array>

#include "psi4/libmints/wavefunction.h"
#include "psi4/libdpd/dpd.h"

#include "ccenergy/MOInfo.h"
#include "ccenergy/Params.h"
#include "ccenergy/Local.h"

namespace psi {
class Options;
struct dpdfile2;
struct dpdbuf4;
struct iwlbuf;
}  // namespace psi

namespace psi {
namespace ccenergy {

class CCEnergyWavefunction : public Wavefunction {
   public:
    CCEnergyWavefunction(std::shared_ptr<Wavefunction> reference_wavefunction, Options &options);
    ~CCEnergyWavefunction() override;
    double compute_energy() override;

    std::map<std::string, SharedMatrix> get_amplitudes();
    // (index within irrep, irrep) -> global index
    std::map<std::tuple<int, int>, int> total_indices;

   private:
    /* setup, info and teardown */
    void init();
    void init_io();
    void exit_io();
    void cleanup();
    void status(const char *, std::string);
    void title();

    /* calculation info */
    void get_moinfo();
    void get_params(Options &);

    /* amplitude handling */
    void init_amps();
    void sort_amps();
    void tsave();
    void ccdump();
    void spinad_amps();
    void amp_write();
    void checkpoint();

    /* intermediates */
    void update();
    void Fae_build();
    void Fmi_build();
    void Fme_build();
    void Wmnij_build();
    void Wmbej_build();
    void purge_Wabei();
    void purge_Wmnij();
    void purge_Wmnie();
    void purge_Wmbij();
    void purge_Wamef();
    void FaetT2();
    void FmitT2();
    void WmnijT2();
    void WmbejT2();
    void BT2();
    void CT2();
    void DT2();
    void ET2();
    void FT2();
    void ZT2();
    void dijabT2();
    void t1_build();
    void tau_build();
    void taut_build();
    void Z_build();
    void Y_build();
    void X_build();
    void t2_build();
    int converged(double);

    /* DPD cache */
    void init_priority_list();
    int **cacheprep_uhf(int level, int *cachefiles);
    int **cacheprep_rhf(int level, int *cachefiles);
    void cachedone_rhf(int **cachelist);
    void cachedone_uhf(int **cachelist);

    /* Brueckner */
    bool rotate();
    double **fock_build(double **D);

    /* CC2 / CC3 */
    void cc2_fmiT2();
    void cc2_faeT2();
    void cc2_WmbijT2();
    void cc2_WabeiT2();
    void cc2_WabijT2();
    void cc2_Wmnij_build();
    void cc2_Wmbij_build();
    void cc2_Wabei_build();
    void cc2_t2_build();
    void FT2_CC2();
    void purge_cc2_Wmnij();
    void purge_cc2_Wmbij();
    void purge_cc2_Wabei();
    void t1_ijab();
    void cc3_Wmnie();
    void cc3_Wamef();
    void cc3_Wmnij();
    void cc3_Wmbij();
    void cc3_Wabei();
    void cc3();

    /* energies */
    double energy();
    double mp2_energy();
    double uhf_mp2_energy();
    double rhf_mp2_energy();
    void one_step();
    void denom();
    // Grab pair energies from the OOVV block of C2 (not T2), storing in the input vectors.
    // For HF orbitals, the correlation energy is the sum of pair energies.
    void pair_energies(std::vector<double>& epair_aa, std::vector<double>& epair_ab) const;
    // Use input AA and AB pair energies to save pair energies to the wavefunction and print them.
    void print_pair_energies(const std::vector<double>& emp2_aa, const std::vector<double>& emp2_ab, const std::vector<double>& ecc_aa, const std::vector<double>& ecc_ab);
    // Form density-fitted integrals and connect them to DPD. Only implemented for RHF.
    void form_df_ints(Options &options, int **cachelist, int *cachefiles);

    /* diagnostics */
    void analyze();
    double diagnostic();
    double d1diag();
    double new_d1diag();
    double new_d1diag_t1_rohf();
    double d2diag();
    double d1diag_t1_rhf();
    double d1diag_t1_rohf();
    double d2diag_rhf();

    /* local correlation */
    void lmp2();
    void local_filter_T1(dpdfile2 *T1);
    void local_filter_T2(dpdbuf4 *T2);
    void local_init();
    void local_done();

    /* AO basis */
    void BT2_AO();
    void halftrans(dpdbuf4 *Buf1, int dpdnum1, dpdbuf4 *Buf2, int dpdnum2, double ***C1, double ***C2, int nirreps,
                   int **mo_row, int **so_row, int *mospi_left, int *mospi_right, int *sospi, int type, double alpha,
                   double beta);
    int AO_contribute(struct iwlbuf *InBuf, dpdbuf4 *tau1_AO, dpdbuf4 *tau2_AO);

    double rhf_energy();
    double uhf_energy();
    double rohf_energy();
    void rhf_fock_build(double **fock, double **D);
    void uhf_fock_build(double **fock_a, double **fock_b, double **D_a, double **D_b);

    /* DIIS */
    void diis(int iter);
    void diis_RHF(int);
    void diis_ROHF(int);
    void diis_UHF(int);
    void diis_invert_B(double **B, double *C, int dimension, double tolerance);

    /* member variables */
    Dimension act_occpi_;
    Dimension act_virpi_;
    MOInfo moinfo_;
    Params params_;
    Local local_;
    std::array<dpd_file4_cache_entry, 113> cache_priority_list_;
};

}  // namespace ccenergy
}  // namespace psi

#endif  // CCWAVE_H

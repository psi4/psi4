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

#ifndef CCWAVE_H
#define CCWAVE_H

#include "psi4/libmints/wavefunction.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"

namespace psi {
class Wavefunction;
class Options;
struct dpd_file4_cache_entry;
struct dpdfile2;
struct dpdbuf4;
}

namespace psi { namespace ccenergy {

class CCEnergyWavefunction : public Wavefunction
{
public:
    CCEnergyWavefunction(std::shared_ptr<Wavefunction> reference_wavefunction, Options &options);
    virtual ~CCEnergyWavefunction();

    double compute_energy();

private:

    /* setup, info and teardown */
    void init();
    void init_io();
    void init_ioff(void);
    void exit_io(void);
    void cleanup(void);
    void status(const char *, std::string );
    void title(void);

    /* calculation info */
    void get_moinfo(void);
    void get_params(Options &);

    /* amplitude handling */
    void init_amps(void);
    void sort_amps(void);
    void tsave(void);
    void ccdump(void);
    void spinad_amps(void);
    void amp_write(void);
    void checkpoint(void);

    /* intermediates */
    void update(void);
    void Fae_build(void);
    void Fmi_build(void);
    void Fme_build(void);
    void Wmnij_build(void);
    void Wmbej_build(void);
    void purge_Wabei(void);
    void purge_Wmnij(void);
    void purge_Wmnie(void);
    void purge_Wmbij(void);
    void purge_Wamef(void);
    void FaetT2(void);
    void FmitT2(void);
    void WmnijT2(void);
    void WmbejT2(void);
    void BT2(void);
    void CT2(void);
    void DT2(void);
    void ET2(void);
    void FT2(void);
    void ZT2(void);
    void dijabT2(void);
    void t1_build(void);
    void tau_build(void);
    void taut_build(void);
    void Z_build(void);
    void Y_build(void);
    void X_build(void);
    void t2_build(void);
    int converged(double);

    /* DPD cache */
    void init_priority_list(void);
    int **cacheprep_uhf(int level, int *cachefiles);
    int **cacheprep_rhf(int level, int *cachefiles);
    void cachedone_rhf(int **cachelist);
    void cachedone_uhf(int **cachelist);

    /* Brueckner */
    int rotate(void);
    double **fock_build(double **D);

    /* CC2 / CC3 */
    void cc2_fmiT2(void);
    void cc2_faeT2(void);
    void cc2_WmbijT2(void);
    void cc2_WabeiT2(void);
    void cc2_WabijT2(void);
    void cc2_Wmnij_build(void);
    void cc2_Wmbij_build(void);
    void cc2_Wabei_build(void);
    void cc2_t2_build(void);
    void FT2_CC2(void);
    void purge_cc2_Wmnij(void);
    void purge_cc2_Wmbij(void);
    void purge_cc2_Wabei(void);
    void t1_ijab(void);
    void cc3_Wmnie(void);
    void cc3_Wamef(void);
    void cc3_Wmnij(void);
    void cc3_Wmbij(void);
    void cc3_Wabei(void);
    void cc3(void);

    /* energies */
    double energy(void);
    double mp2_energy(void);
    double uhf_mp2_energy(void);
    double rhf_mp2_energy(void);
    void one_step(void);
    void denom(void);
    void pair_energies(double** epair_aa, double** epair_ab);
    void print_pair_energies(double* emp2_aa, double* emp2_ab, double* ecc_aa,
                             double* ecc_ab);

    void form_df_ints(Options &options, int **cachelist, int *cachefiles, dpd_file4_cache_entry *priority);

    /* diagnostics */
    void analyze(void);
    double diagnostic(void);
    double d1diag(void);
    double new_d1diag(void);
    double new_d1diag_t1_rohf(void);
    double d2diag(void);
    double d1diag_t1_rhf(void);
    double d1diag_t1_rohf(void);
    double d2diag_rhf(void);

    /* local correlation */
    void lmp2(void);
    void local_filter_T1(dpdfile2 *T1);
    void local_filter_T2(dpdbuf4 *T2);
    void local_init(void);
    void local_done(void);


    /* AO basis */
    void BT2_AO(void);
    void halftrans(dpdbuf4 *Buf1, int dpdnum1, dpdbuf4 *Buf2, int dpdnum2, double ***C1, double ***C2,
                   int nirreps, int **mo_row, int **so_row, int *mospi_left, int *mospi_right,
                   int *sospi, int type, double alpha, double beta);
    int AO_contribute(struct iwlbuf *InBuf, dpdbuf4 *tau1_AO, dpdbuf4 *tau2_AO);


    double rhf_energy(void);
    double uhf_energy(void);
    double rohf_energy(void);
    void rhf_fock_build(double **fock, double  **D);
    void uhf_fock_build(double **fock_a, double **fock_b, double **D_a, double **D_b);

    /* DIIS */
    void diis(int iter);
    void diis_RHF(int);
    void diis_ROHF(int);
    void diis_UHF(int);
    void diis_invert_B(double** B, double* C, int dimension, double tolerance);

    /* member variables */
    int *ioff_;
    MOInfo moinfo_;
    Params params_;
    Local local_;
    dpd_file4_cache_entry *cache_priority_list_;
};

}}

#endif // CCWAVE_H

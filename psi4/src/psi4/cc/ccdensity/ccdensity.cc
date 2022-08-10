/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
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

/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here
*/
/*
**  CCDENSITY: Program to calculate the coupled-cluster one- and
**             two-particle densities.
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "ccdensity.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#include "globals.h"
#include "psi4/cc/ccwave.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libiwl/iwl.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libqt/qt.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/matrix.h"
#include "psi4/psi4-dec.h"
#include "psi4/psifiles.h"
namespace psi {

class Molecule;
namespace ccdensity {

void init_io();
void title();
void get_moinfo(std::shared_ptr<Wavefunction> wfn);
void get_frozen();
void get_params(Options &options);
void exit_io();
void onepdm(const struct RHO_Params&);
void sortone(const struct RHO_Params&);
void twopdm();
void energy(const struct RHO_Params&);
// void resort_tei();
// void resort_gamma();
void lag(const struct RHO_Params& rho_params);
void build_X();
void build_A();
void build_Z();
void relax_I();
void relax_D(const struct RHO_Params& rho_params);
void sortI();
void fold(const struct RHO_Params& rho_params);
void deanti(const struct RHO_Params& rho_params);
void add_ref_RHF(struct iwlbuf *);
void add_ref_ROHF(struct iwlbuf *);
void add_ref_UHF(struct iwlbuf *, struct iwlbuf *, struct iwlbuf *);
void add_core_ROHF(struct iwlbuf *);
// void add_core_UHF(struct iwlbuf *, struct iwlbuf *, struct iwlbuf *);
void dump_RHF(struct iwlbuf *, const struct RHO_Params& rho_params);
void dump_ROHF(struct iwlbuf *, const struct RHO_Params& rho_params);
void dump_UHF(struct iwlbuf *, struct iwlbuf *, struct iwlbuf *, const struct RHO_Params& rho_params);
void kinetic(std::shared_ptr<Wavefunction> wfn);
void probable();
int **cacheprep_rhf(int level, int *cachefiles);
int **cacheprep_uhf(int level, int *cachefiles);
void cachedone_rhf(int **cachelist);
void cachedone_uhf(int **cachelist);
void setup_LR(const struct RHO_Params&);
void G_build();
void x_oe_intermediates(const struct RHO_Params&);
void x_onepdm(const struct RHO_Params&);
void x_te_intermediates();
void x_Gijkl();
void x_Gabcd();
void x_Gibja();
void x_Gijka();
void x_Gijab();
void x_Gciab();
void V_build_x();
void x_xi1();
void x_xi_zero();
void x_xi2();
void x_xi_oe_intermediates();
// void G_norm();
void zero_onepdm(const struct RHO_Params& rho_params);
void zero_twopdm();
void get_rho_params(Options &options);
void get_td_params(Options &options);
void td_setup(const struct TD_Params& S);
void tdensity(const struct TD_Params& S);
void td_print();
void oscillator_strength(ccenergy::CCEnergyWavefunction& wfn, struct TD_Params *S);
void rotational_strength(ccenergy::CCEnergyWavefunction& wfn, struct TD_Params *S);
void cleanup();
void td_cleanup();
void x_oe_intermediates_rhf(const struct RHO_Params& rho_params);
void x_te_intermediates_rhf();
void x_xi_intermediates();
void V_build();
void V_cc2();
void ex_tdensity(char hand, const struct TD_Params& S, const struct TD_Params& U);
void ex_td_setup(const struct TD_Params& S, const struct TD_Params& U);
void ex_td_cleanup();
void ex_oscillator_strength(ccenergy::CCEnergyWavefunction& wfn, struct TD_Params *S, struct TD_Params *U,
                            struct XTD_Params *xtd_data);
void ex_rotational_strength(ccenergy::CCEnergyWavefunction& wfn, struct TD_Params *S, struct TD_Params *U, struct XTD_Params *xtd_data);
void ex_td_print(std::vector<struct XTD_Params>);
SharedMatrix block_to_matrix(double **);

PsiReturnType ccdensity(std::shared_ptr<ccenergy::CCEnergyWavefunction> ref_wfn, Options &options) {
    int i;
    int **cachelist, *cachefiles;
    struct iwlbuf OutBuf;
    struct iwlbuf OutBuf_AA, OutBuf_BB, OutBuf_AB;
    dpdfile2 D;
    double tval;

    init_io();
    title();
    /*  get_frozen(); */
    get_params(options);
    get_moinfo(ref_wfn);
    get_rho_params(options);

    if ((moinfo.nfzc || moinfo.nfzv) && params.relax_opdm) {
        outfile->Printf("\n\tGradients/orbital relaxation involving frozen orbitals not yet available.\n");
        throw PsiException("ccdensity: error", __FILE__, __LINE__);
    }

    cachefiles = init_int_array(PSIO_MAXUNIT);

    if (params.ref == 0 || params.ref == 1) { /** RHF or ROHF **/
        cachelist = cacheprep_rhf(params.cachelev, cachefiles);

        std::vector<int *> spaces;
        spaces.push_back(moinfo.occpi);
        spaces.push_back(moinfo.occ_sym.data());
        spaces.push_back(moinfo.virtpi);
        spaces.push_back(moinfo.vir_sym.data());
        delete dpd_list[0];
        dpd_list[0] = new DPD(0, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, nullptr, 2, spaces);
        dpd_set_default(0);

    } else if (params.ref == 2) { /** UHF **/
        cachelist = cacheprep_uhf(params.cachelev, cachefiles);

        std::vector<int *> spaces;
        spaces.push_back(moinfo.aoccpi);
        spaces.push_back(moinfo.aocc_sym.data());
        spaces.push_back(moinfo.avirtpi);
        spaces.push_back(moinfo.avir_sym.data());
        spaces.push_back(moinfo.boccpi);
        spaces.push_back(moinfo.bocc_sym.data());
        spaces.push_back(moinfo.bvirtpi);
        spaces.push_back(moinfo.bvir_sym.data());
        dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, nullptr, 4, spaces);
    }

    for (i = 0; i < params.nstates; ++i) {
        /* CC_GLG will contain L, or R0*L + Zeta, if relaxed and zeta is available */
        /* CC_GL will contain L */
        setup_LR(rho_params[i]);

        /* Calculate Xi, put Xi in EOM_XI, and quit */
        if (params.calc_xi) {
            /* these intermediates go into EOM_TMP and are used to compute Xi;
               they may be reused to compute the excited-state density matrix */
            if (params.ref == 0) {
                x_oe_intermediates_rhf(rho_params[i]);
                x_te_intermediates_rhf();
            } else {
                x_oe_intermediates(rho_params[i]);
                x_te_intermediates();
            }
            x_xi_intermediates(); /*Xi intermediates put in EOM_TMP_XI */
            x_xi_zero();          /* make blank Xi */
            x_xi1();
            x_xi2();
            dpd_close(0);
            if (params.ref == 2)
                cachedone_uhf(cachelist);
            else
                cachedone_rhf(cachelist);
            free(cachefiles);
            cleanup();
            psio_close(PSIF_EOM_TMP_XI, 0); /* delete EOM_TMP_XI */
            psio_open(PSIF_EOM_TMP_XI, PSIO_OPEN_NEW);
            exit_io();
            return Success;
        }

        /* compute ground state parts of onepdm or put zeroes there */
        if (((rho_params[i].L_irr == rho_params[i].G_irr) || (params.use_zeta))) {
            zero_onepdm(rho_params[i]);
            onepdm(rho_params[i]);
        } else
            zero_onepdm(rho_params[i]);

        /* if the one-electron excited-state intermediates are not already on disk (from a Xi
           calculation, compute them.  They are nearly all necessary to compute the excited-state
           onepdm. Then complete excited-state onepdm.*/
        if (!rho_params[i].R_ground) {
            x_oe_intermediates(rho_params[i]); /* change to x_oe_intermediates_rhf() when rho gets spin-adapted */
            x_onepdm(rho_params[i]);
        }

        /* begin construction of twopdm */
        if (!params.onepdm) {
            /* Compute intermediates for construction of ground-state twopdm */
            if ((params.L_irr == params.G_irr) || (params.use_zeta)) {
                if (params.wfn == "CC2" && params.dertype == 1)
                    V_cc2();
                else {
                    V_build(); /* uses CC_GLG, writes tau2*L2 to CC_MISC */
                    G_build(); /* uses CC_GLG, writes t2*L2 to CC_GLG */
                }
            }
            /* Compute ground-state twopdm or ground-state-like contributions to the excited twodpm */
            if ((params.L_irr == params.G_irr) || (params.use_zeta))
                twopdm();
            else
                zero_twopdm();

            /* Compute intermediates for construction of excited-state twopdm */
            if (!params.ground) {
                x_te_intermediates(); /* change to x_te_intermediates_rhf() when rho gets spin-adapted */
                V_build_x();          /* uses CC_GL, writes t2*L2 to EOM_TMP */

                /* add in non-R0 parts of onepdm and twopdm */
                x_Gijkl();
                x_Gabcd();
                x_Gibja();
                x_Gijka();
                x_Gciab();
                x_Gijab();
            }
        }

        sortone(rho_params[i]); /* puts full 1-pdm into moinfo.opdm */
        if (!params.onepdm) {
            if (!params.aobasis && params.debug_) energy(rho_params[i]);

            kinetic(ref_wfn); /* puts kinetic energy integrals into MO basis */

            lag(rho_params[i]); /* builds the orbital lagrangian pieces, I */

            build_X(); /* builds orbital rotation gradient X */
            build_A(); /* construct MO Hessian A */
            build_Z(); /* solves the orbital Z-vector equations */

            relax_I(); /* adds orbital response contributions to Lagrangian */

            if (params.relax_opdm) {
                relax_D(rho_params[i]); /* adds orbital response contributions to onepdm */
            }
            sortone(rho_params[i]); /* builds large moinfo.opdm matrix */
            sortI();                /* builds large lagrangian matrix I */
            fold(rho_params[i]);
            deanti(rho_params[i]);
        }

        if (params.ref == 0) { /** RHF **/

            iwl_buf_init(&OutBuf, PSIF_MO_TPDM, params.tolerance, 0, 0);

            add_core_ROHF(&OutBuf);
            add_ref_RHF(&OutBuf);

            if (params.onepdm_grid_dump) dx_write(ref_wfn, options, moinfo.opdm);

            dump_RHF(&OutBuf, rho_params[i]);

            iwl_buf_flush(&OutBuf, 1);
            iwl_buf_close(&OutBuf, 1);
        } else if (params.ref == 1) { /** ROHF **/

            iwl_buf_init(&OutBuf, PSIF_MO_TPDM, params.tolerance, 0, 0);

            add_core_ROHF(&OutBuf);
            add_ref_ROHF(&OutBuf);

            dump_ROHF(&OutBuf, rho_params[i]);

            iwl_buf_flush(&OutBuf, 1);
            iwl_buf_close(&OutBuf, 1);
        } else if (params.ref == 2) { /** UHF **/

            iwl_buf_init(&OutBuf_AA, PSIF_MO_AA_TPDM, params.tolerance, 0, 0);
            iwl_buf_init(&OutBuf_BB, PSIF_MO_BB_TPDM, params.tolerance, 0, 0);
            iwl_buf_init(&OutBuf_AB, PSIF_MO_AB_TPDM, params.tolerance, 0, 0);

            /*    add_core_UHF(&OutBuf_AA, &OutBuf_BB, &OutBuf_AB); */
            add_ref_UHF(&OutBuf_AA, &OutBuf_BB, &OutBuf_AB);

            dump_UHF(&OutBuf_AA, &OutBuf_BB, &OutBuf_AB, rho_params[i]);

            iwl_buf_flush(&OutBuf_AA, 1);
            iwl_buf_flush(&OutBuf_BB, 1);
            iwl_buf_flush(&OutBuf_AB, 1);
            iwl_buf_close(&OutBuf_AA, 1);
            iwl_buf_close(&OutBuf_BB, 1);
            iwl_buf_close(&OutBuf_AB, 1);
        }
        std::shared_ptr<Matrix> Ca = ref_wfn->Ca();
        std::shared_ptr<Matrix> Cb = ref_wfn->Cb();

        Dimension nmopi = ref_wfn->nmopi();
        Dimension frzvpi = ref_wfn->frzvpi();

        // Grab the GS OPDM and set it in the ref_wfn object
        auto Pa = std::make_shared<Matrix>("P alpha", Ca->colspi(), Ca->colspi());
        auto Pb = std::make_shared<Matrix>("P beta", Cb->colspi(), Cb->colspi());
        int mo_offset = 0;

        /*
        for (int h = 0; h < Ca->nirrep(); h++) {
            int nmo = nmopi[h];
            int nfv = frzvpi[h];
            int nmor = nmo - nfv;
            if (!nmo || !nmor) continue;

            // Loop over QT, convert to Pitzer
            double **Pap = Pa->pointer(h);
            double **Pbp = Pb->pointer(h);
            for (int i = 0; i < nmor; i++) {
                for (int j = 0; j < nmor; j++) {
                    int I = moinfo.pitzer2qt[i + mo_offset];
                    int J = moinfo.pitzer2qt[j + mo_offset];
                    if (ref_wfn->same_a_b_dens())
                        Pap[i][j] = moinfo.opdm[I][J];
                    else {
                        Pap[i][j] = moinfo.opdm_a[I][J];
                        Pbp[i][j] = moinfo.opdm_b[I][J];
                    }
                }
            }
            mo_offset += nmo;
        }
        */

        SharedMatrix Pa_x, Pb_x;
        if (ref_wfn->same_a_b_dens()) {
            Pa_x = block_to_matrix(moinfo.opdm);
        } else {
            Pa_x = block_to_matrix(moinfo.opdm_a);
            Pb_x = block_to_matrix(moinfo.opdm_b);
        }

        std::string short_name;
        if (params.wfn == "EOM_CC2") {
            short_name = "CC2";
        } else if (params.wfn == "EOM_CC3") {
            short_name = "CC3";
        } else if (params.wfn == "EOM_CCSD") {
            short_name = "CCSD";
        }

        /* Transform Da/b to so basis and set in wfn */
        // If this becomes a wavefunction subclass someday, just redefine the densities directly.
        if (ref_wfn->same_a_b_dens()) {
            Pa_x->scale(0.5);
            auto Pa_so = linalg::triplet(ref_wfn->Ca(), Pa_x, ref_wfn->Ca(), false, false, true);
            if (i == 0) {
                auto ref_Da_so = ref_wfn->Da();
                ref_Da_so->copy(Pa_so);
            } else {
                density_saver(*ref_wfn, &(rho_params[i]), "ALPHA", Pa_so);
            }
        } else {
            auto Pa_so = linalg::triplet(ref_wfn->Ca(), Pa_x, ref_wfn->Ca(), false, false, true);
            auto Pb_so = linalg::triplet(ref_wfn->Cb(), Pb_x, ref_wfn->Cb(), false, false, true);
            if (i == 0) {
                auto ref_Da_so = ref_wfn->Da();
                auto ref_Db_so = ref_wfn->Db();
                ref_Da_so->copy(Pa_so);
                ref_Db_so->copy(Pb_so);
            } else {
                density_saver(*ref_wfn, &(rho_params[i]), "ALPHA", Pa_so);
                density_saver(*ref_wfn, &(rho_params[i]), "BETA", Pb_so);
            }
        }

        // For psivar scraper

        // Process::environment.globals["CCname ROOT n DIPOLE"]
        // Process::environment.globals["CCname ROOT n DIPOLE - h TRANSITION"]
        // Process::environment.globals["CCname ROOT n (h) DIPOLE"]
        // Process::environment.globals["CCname ROOT n QUADRUPOLE"]
        // Process::environment.globals["CCname ROOT n QUADRUPOLE - h TRANSITION"]
        // Process::environment.globals["CCname ROOT n (h) QUADRUPOLE"]

        free_block(moinfo.opdm);
        free_block(moinfo.opdm_a);
        free_block(moinfo.opdm_b);

        psio_close(PSIF_CC_TMP, 0);
        psio_open(PSIF_CC_TMP, PSIO_OPEN_NEW);
        psio_close(PSIF_EOM_TMP0, 0);
        psio_open(PSIF_EOM_TMP0, PSIO_OPEN_NEW);
        psio_close(PSIF_EOM_TMP1, 0);
        psio_open(PSIF_EOM_TMP1, PSIO_OPEN_NEW);
        psio_close(PSIF_CC_GLG, 0);
        psio_open(PSIF_CC_GLG, PSIO_OPEN_NEW);
        psio_close(PSIF_CC_GL, 0);
        psio_open(PSIF_CC_GL, PSIO_OPEN_NEW);
        psio_close(PSIF_CC_GR, 0);
        psio_open(PSIF_CC_GR, PSIO_OPEN_NEW);
        if (!params.calc_xi) {
            psio_close(PSIF_EOM_TMP, 0);
            psio_open(PSIF_EOM_TMP, PSIO_OPEN_NEW);
        }
    }

    if (params.transition) {
        MintsHelper mints(ref_wfn->basisset(), options, 0);
        get_td_params(options);
        for (i = 0; i < params.nstates; i++) {
            td_setup(td_params[i]);
            tdensity(td_params[i]);
            outfile->Printf("Doing transition\n");
            oscillator_strength(*ref_wfn, &(td_params[i]));
            outfile->Printf("Doing transition\n");
            if (params.ref == 0) {
                rotational_strength(*ref_wfn, &(td_params[i]));
            }
            outfile->Printf("Doing transition\n");
            td_cleanup();
        }
        outfile->Printf("Doing transition\n");
        td_print();
        outfile->Printf("Doing transition\n");

        /* Excited State Transition Data */
        //  The convention is that the transition is one of absorption.
        //  That is to say - we always go from a lower excited state
        //  to a higher one - which maintains a defintion for
        //  labeling the LTD and the RTD.
        int j;

        if (params.nstates > 1) {  // Can't do this with one excited state.
            outfile->Printf("\n\t*********************************************************************\n");
            outfile->Printf("\t*********************************************************************\n");
            outfile->Printf("\t************                                             ************\n");
            outfile->Printf("\t************ Excited State-Excited State Transition Data ************\n");
            outfile->Printf("\t************                                             ************\n");
            outfile->Printf("\t*********************************************************************\n");
            outfile->Printf("\t*********************************************************************\n\n");

            std::vector<struct XTD_Params> xtd_params;
            struct XTD_Params xtd_data;
            int state1;
            int state2;

            for (i = 0; i < (params.nstates - 1); i++) {
                for (j = 0; j <= i; j++) {
                    //- Set States
                    if (td_params[j].cceom_energy <= td_params[i + 1].cceom_energy) {
                        state1 = j;
                        state2 = i + 1;
                    } else {
                        state1 = i + 1;
                        state2 = j;
                    }
                    /*
                              outfile->Printf( "State %d%s Energy = %20.12lf\n",
                                    td_params[state1].root+1,moinfo.labels[td_params[state1].irrep],td_params[state1].cceom_energy);
                              outfile->Printf( "State %d%s Energy = %20.12lf\n",
                                    td_params[state1].root+1,moinfo.labels[td_params[state2].irrep],td_params[state2].cceom_energy);
                    */

                    //- <Lx|O|Ry> (y>x)
                    outfile->Printf("\n\t*** Computing <%d%2s|X{pq}}|%d%2s> (LEFT) Transition Density ***\n\n",
                                    td_params[state1].root + 1, moinfo.labels[td_params[state1].irrep].c_str(),
                                    td_params[state2].root + 1, moinfo.labels[td_params[state2].irrep].c_str());
                    ex_td_setup(td_params[state1], td_params[state2]);
                    outfile->Printf("\t\t*** LTD Setup complete.\n");

                    ex_tdensity('l', td_params[state1], td_params[state2]);

                    //- Clean out Amp Files (Might not be necessary, but seems to be.)
                    ex_td_cleanup();

                    //- <Ly|O|Rx> (y>x)
                    outfile->Printf("\n\t*** Computing <%d%2s|X{pq}}|%d%2s> (RIGHT) Transition Density ***\n\n",
                                    td_params[state2].root + 1, moinfo.labels[td_params[state2].irrep].c_str(),
                                    td_params[state1].root + 1, moinfo.labels[td_params[state1].irrep].c_str());
                    ex_td_setup(td_params[state2], td_params[state1]);
                    outfile->Printf("\t\t*** RTD Setup complete.\n");

                    ex_tdensity('r', td_params[state2], td_params[state1]);

                    outfile->Printf("\n\t*** %d%s -> %d%s transition densities complete.\n", td_params[state1].root + 1,
                                    moinfo.labels[td_params[state1].irrep].c_str(), td_params[state2].root + 1,
                                    moinfo.labels[td_params[state2].irrep].c_str());

                    // ex_oscillator_strength(&(td_params[j]),&(td_params[i+1]), &xtd_data);
                    ex_oscillator_strength(*ref_wfn, &(td_params[state1]), &(td_params[state2]), &xtd_data);
                    if (params.ref == 0) {
                        // ex_rotational_strength(&(td_params[j]),&(td_params[i+1]), &xtd_data);
                        ex_rotational_strength(*ref_wfn, &(td_params[state1]), &(td_params[state2]), &xtd_data);
                    }

                    xtd_params.push_back(xtd_data);

                    td_cleanup();
                }
            }
            td_print();
            ex_td_print(xtd_params);
        }

    }  // End params.transition IF loop

    dpd_close(0);

    if (params.ref == 2)
        cachedone_uhf(cachelist);
    else
        cachedone_rhf(cachelist);
    free(cachefiles);

    cleanup();
    exit_io();
    return Success;
}

// must be fixed with options for excited state densities
void init_io() {
    int i, num_unparsed;
    char *argv_unparsed[100];

    params.onepdm = 0;
    params.prop_all = 0;
    params.calc_xi = 0;
    params.restart = 0;
    params.use_zeta = 0;
    params.transition = 0;

    /*
      for (i=1, num_unparsed=0; i<argc; ++i) {
        if(!strcmp(argv[i], "--onepdm")) {
          params.onepdm = 1; // generate ONLY the onepdm (for one-electron properties)
        }
        else if (!strcmp(argv[i],"--use_zeta")) {
          params.use_zeta = 1;
          params.ground = 0;
          params.restart = 1;
        }
        else if (!strcmp(argv[i],"--calc_xi")) {
          params.calc_xi = 1;
          params.ground = 0;
          params.restart = 0;
        }
        else if (!strcmp(argv[i],"--prop_all")) {
          params.prop_all = 1;
        }
        else if (!strcmp(argv[i],"--transition")) {
          params.transition = 1;
          params.relax_opdm = 0;
          params.ground = 0;
        }
        else {
          argv_unparsed[num_unparsed++] = argv[i];
        }
      }
    */

    timer_on("ccdensity");

    for (i = PSIF_CC_MIN; i <= PSIF_CC_MAX; i++) psio_open(i, PSIO_OPEN_OLD);
    // erase old files
    psio_close(PSIF_CC_GR, 0);
    psio_close(PSIF_CC_GL, 0);
    psio_close(PSIF_EOM_TMP0, 0);
    psio_open(PSIF_CC_GR, PSIO_OPEN_NEW);
    psio_open(PSIF_CC_GL, PSIO_OPEN_NEW);
    psio_open(PSIF_EOM_TMP0, PSIO_OPEN_NEW);
}

void title() {
    outfile->Printf("\n");
    outfile->Printf("\t\t\t**************************\n");
    outfile->Printf("\t\t\t*                        *\n");
    outfile->Printf("\t\t\t*        CCDENSITY       *\n");
    outfile->Printf("\t\t\t*                        *\n");
    outfile->Printf("\t\t\t**************************\n");
    outfile->Printf("\n");
}

void exit_io() {
    int i;

    /* delete temporary EOM files */
    psio_close(PSIF_EOM_TMP0, 0);
    psio_close(PSIF_EOM_TMP1, 0);
    psio_close(PSIF_CC_GLG, 0);
    psio_open(PSIF_EOM_TMP0, PSIO_OPEN_NEW);
    psio_open(PSIF_EOM_TMP1, PSIO_OPEN_NEW);
    psio_open(PSIF_CC_GLG, PSIO_OPEN_NEW);
    if (!params.calc_xi) {
        psio_close(PSIF_EOM_TMP, 0);
        psio_open(PSIF_EOM_TMP, PSIO_OPEN_NEW);
    }
    if (params.use_zeta) { /* we're done with Xi amplitudes */
        psio_close(PSIF_EOM_XI, 0);
        psio_open(PSIF_EOM_XI, PSIO_OPEN_NEW);
    }

    /* Close all dpd data files here */
    for (i = PSIF_CC_MIN; i <= PSIF_CC_MAX; i++) psio_close(i, 1);

    timer_off("ccdensity");
}

}  // namespace ccdensity
}  // namespace psi

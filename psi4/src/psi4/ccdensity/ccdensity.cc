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
#include "psi4/libpsio/psio.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libiwl/iwl.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/psi4-dec.h"
#include <cmath>
#include "psi4/psifiles.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#include "globals.h"
#include "psi4/libqt/qt.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/oeprop.h"
#include "psi4/libmints/writer.h"
#include "psi4/libmints/molecule.h"

#ifdef USING_PCMSolver
#include "psi4/libpsipcm/psipcm.h"
#endif

namespace psi {
namespace ccdensity {

void init_io(void);
void title(void);
void get_moinfo(std::shared_ptr<Wavefunction> wfn);
void get_frozen(void);
void get_params(Options &options);
void exit_io(void);
void onepdm(struct RHO_Params);
void sortone(struct RHO_Params);
void twopdm(void);
void energy(struct RHO_Params);
// void resort_tei(void);
// void resort_gamma(void);
void lag(struct RHO_Params rho_params);
void build_X(void);
void build_A(void);
void build_Z(void);
void relax_I(void);
void relax_D(struct RHO_Params rho_params);
void sortI(void);
void fold(struct RHO_Params rho_params);
void deanti(struct RHO_Params rho_params);
void add_ref_RHF(struct iwlbuf *);
void add_ref_ROHF(struct iwlbuf *);
void add_ref_UHF(struct iwlbuf *, struct iwlbuf *, struct iwlbuf *);
void add_core_ROHF(struct iwlbuf *);
// void add_core_UHF(struct iwlbuf *, struct iwlbuf *, struct iwlbuf *);
void dump_RHF(struct iwlbuf *, struct RHO_Params rho_params);
void dump_ROHF(struct iwlbuf *, struct RHO_Params rho_params);
void dump_UHF(struct iwlbuf *, struct iwlbuf *, struct iwlbuf *, struct RHO_Params rho_params);
void kinetic(std::shared_ptr<Wavefunction> wfn);
void probable(void);
int **cacheprep_rhf(int level, int *cachefiles);
int **cacheprep_uhf(int level, int *cachefiles);
void cachedone_rhf(int **cachelist);
void cachedone_uhf(int **cachelist);
void setup_LR(struct RHO_Params);
void G_build(void);
void x_oe_intermediates(struct RHO_Params);
void x_onepdm(struct RHO_Params);
void x_te_intermediates(void);
void x_Gijkl(void);
void x_Gabcd(void);
void x_Gibja(void);
void x_Gijka(void);
void x_Gijab(void);
void x_Gciab(void);
void V_build_x(void);
void x_xi1(void);
void x_xi_zero(void);
void x_xi2(void);
void x_xi_oe_intermediates(void);
// void G_norm(void);
void zero_onepdm(struct RHO_Params rho_params);
void zero_twopdm(void);
void get_rho_params(Options &options);
void get_td_params(Options &options);
void td_setup(struct TD_Params S);
void tdensity(struct TD_Params S);
void td_print(void);
void oscillator_strength(std::shared_ptr<Wavefunction> wfn, struct TD_Params *S);
void rotational_strength(MintsHelper &mints, struct TD_Params *S);
void ael(struct RHO_Params *rho_params);
void cleanup(void);
void td_cleanup(void);
void x_oe_intermediates_rhf(struct RHO_Params rho_params);
void x_te_intermediates_rhf(void);
void x_xi_intermediates(void);
void V_build(void);
void V_cc2(void);
void update_F_pcm_rhf(SharedMatrix &MO_PCM_potential);
void update_F_pcm_uhf(SharedMatrix &MO_PCM_potential_A, SharedMatrix &MO_PCM_potential_B);
void ex_tdensity(char hand, struct TD_Params S, struct TD_Params U);
void ex_td_setup(struct TD_Params S, struct TD_Params U);
void ex_td_cleanup();
void ex_oscillator_strength(std::shared_ptr<Wavefunction> wfn, struct TD_Params *S, struct TD_Params *U,
                            struct XTD_Params *xtd_data);
void ex_rotational_strength(MintsHelper &mints, struct TD_Params *S, struct TD_Params *U, struct XTD_Params *xtd_data);
void ex_td_print(std::vector<struct XTD_Params>);

PsiReturnType ccdensity(std::shared_ptr<Wavefunction> ref_wfn, Options &options) {
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
        spaces.push_back(moinfo.occ_sym);
        spaces.push_back(moinfo.virtpi);
        spaces.push_back(moinfo.vir_sym);
        delete dpd_list[0];
        dpd_list[0] = new DPD(0, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, NULL, 2, spaces);
        dpd_set_default(0);

    } else if (params.ref == 2) { /** UHF **/
        cachelist = cacheprep_uhf(params.cachelev, cachefiles);

        std::vector<int *> spaces;
        spaces.push_back(moinfo.aoccpi);
        spaces.push_back(moinfo.aocc_sym);
        spaces.push_back(moinfo.avirtpi);
        spaces.push_back(moinfo.avir_sym);
        spaces.push_back(moinfo.boccpi);
        spaces.push_back(moinfo.bocc_sym);
        spaces.push_back(moinfo.bvirtpi);
        spaces.push_back(moinfo.bvir_sym);
        dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, NULL, 4, spaces);
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

        /* PCM-CC stuff */
        if (options.get_bool("PCM") &&
            (options.get_str("PCM_CC_TYPE") == "PTED")) {  // PTED scheme: R. Cammi, JCP, 131, 164104 (2009)
            //   Build the PCM polarization charges from the CC density matrix at the current iteration.
            //   The polarization charges are then contracted with the electrostatic potential integrals
            //   to build the PCM potential to be added in the CC T and Lambda equations.
            SharedPCM cc_pcm = std::make_shared<PCM>(options, ref_wfn->psio(), moinfo.nirreps, ref_wfn->basisset());
            if (params.ref == 0 || params.ref == 1) /** RHF or ROHF **/
            {
                SharedMatrix MO_OPDM(new Matrix("MO_OPDM", moinfo.nirreps, moinfo.orbspi, moinfo.orbspi));
                SharedMatrix MO_PCM_potential(
                    new Matrix("MO_PCM_potential", moinfo.nirreps, moinfo.orbspi, moinfo.orbspi));

                int offset = 0;
                for (int h = 0; h < moinfo.nirreps; ++h) {
                    for (int i = 0; i < moinfo.orbspi[h]; ++i) {
                        int iabs = moinfo.pitzer2qt[i + offset];
                        for (int j = 0; j < moinfo.orbspi[h]; ++j) {
                            int jabs = moinfo.pitzer2qt[j + offset];
                            MO_OPDM->set(h, i, j, moinfo.opdm[iabs][jabs]);
                        }
                    }
                    offset += moinfo.orbspi[h];
                }
                // MO_OPDM->print();
                SharedMatrix C = ref_wfn->Ca();
                SharedMatrix SO_OPDM = SharedMatrix(new Matrix("SO_OPDM", ref_wfn->nsopi(), ref_wfn->nsopi()));
                SO_OPDM->back_transform(MO_OPDM, C);
                // SO_OPDM->print();
                // ref_wfn->Da()->print();
                // We have to get the 0.5 * (QV) energy contribution
                double Epcm_correlated = cc_pcm->compute_E(SO_OPDM, psi::PCM::EleOnly);
                Process::environment.globals["PCM-CC-PTED CORRELATED POLARIZATION ENERGY"] = Epcm_correlated;
                double E_correlation = Process::environment.globals["CURRENT CORRELATION ENERGY"];
                E_correlation -= Epcm_correlated;
                Process::environment.globals["CURRENT CORRELATION ENERGY"] = E_correlation;
                double E = Process::environment.globals["CURRENT ENERGY"];
                E -= Epcm_correlated;
                Process::environment.globals["CURRENT ENERGY"] = E;
                // outfile->Printf("Epol_correlated = %20.12f\n", Epol_correlated);
                SharedMatrix SO_PCM_potential = cc_pcm->compute_V_electronic();  // This is in SO basis
                // We now transform it to MO basis...
                MO_PCM_potential->transform(SO_PCM_potential, C);
                // MO_PCM_potential->print();
                update_F_pcm_rhf(MO_PCM_potential);
            } else if (params.ref == 2) { /** UHF case **/
                SharedMatrix MO_OPDM_A(new Matrix("MO_OPDM_A", moinfo.nirreps, moinfo.orbspi, moinfo.orbspi));
                SharedMatrix MO_OPDM_B(new Matrix("MO_OPDM_B", moinfo.nirreps, moinfo.orbspi, moinfo.orbspi));

                SharedMatrix MO_PCM_potential_A(
                    new Matrix("MO_PCM_potential_A", moinfo.nirreps, moinfo.orbspi, moinfo.orbspi));
                SharedMatrix MO_PCM_potential_B(
                    new Matrix("MO_PCM_potential_B", moinfo.nirreps, moinfo.orbspi, moinfo.orbspi));

                int offset = 0;
                for (int h = 0; h < moinfo.nirreps; ++h) {
                    for (int i = 0; i < moinfo.orbspi[h]; ++i) {
                        int iabs_A = moinfo.pitzer2qt[i + offset];
                        int iabs_B = moinfo.pitzer2qt[i + offset];
                        for (int j = 0; j < moinfo.orbspi[h]; ++j) {
                            int jabs_A = moinfo.pitzer2qt[j + offset];
                            int jabs_B = moinfo.pitzer2qt[j + offset];
                            MO_OPDM_A->set(h, i, j, moinfo.opdm_a[iabs_A][jabs_A]);
                            MO_OPDM_B->set(h, i, j, moinfo.opdm_b[iabs_B][jabs_B]);
                        }
                    }
                    offset += moinfo.orbspi[h];
                }
                // MO_OPDM->print();
                SharedMatrix Ca = ref_wfn->Ca();
                SharedMatrix Cb = ref_wfn->Cb();
                SharedMatrix SO_OPDM_A = SharedMatrix(new Matrix("SO_OPDM_A", ref_wfn->nsopi(), ref_wfn->nsopi()));
                SharedMatrix SO_OPDM_B = SharedMatrix(new Matrix("SO_OPDM_B", ref_wfn->nsopi(), ref_wfn->nsopi()));
                SO_OPDM_A->back_transform(MO_OPDM_A, Ca);
                SO_OPDM_B->back_transform(MO_OPDM_B, Cb);
                SO_OPDM_A->add(SO_OPDM_B);
                // SO_OPDM->print();
                // ref_wfn->Da()->print();
                // We have to get the 0.5 * (QV) energy contribution
                double Epcm_correlated = cc_pcm->compute_E(SO_OPDM_A, PCM::EleOnly);
                Process::environment.globals["PCM-CC-PTED CORRELATED POLARIZATION ENERGY"] = Epcm_correlated;
                double E_correlation = Process::environment.globals["CURRENT CORRELATION ENERGY"];
                E_correlation -= Epcm_correlated;
                Process::environment.globals["CURRENT CORRELATION ENERGY"] = E_correlation;
                double E = Process::environment.globals["CURRENT ENERGY"];
                E -= Epcm_correlated;
                Process::environment.globals["CURRENT ENERGY"] = E;
                // outfile->Printf("Epol_correlated = %20.12f\n", Epol_correlated);
                SharedMatrix SO_PCM_potential = cc_pcm->compute_V_electronic();  // This is in SO basis
                // We now transform it to MO basis...
                MO_PCM_potential_A->transform(SO_PCM_potential, Ca);
                MO_PCM_potential_B->transform(SO_PCM_potential, Cb);
                // MO_PCM_potential->print();
                update_F_pcm_uhf(MO_PCM_potential_A, MO_PCM_potential_B);
            }
            double Epte = Process::environment.globals["PCM-CC-PTE CORRELATION ENERGY"];
            double Epte_s = Process::environment.globals["PCM-CC-PTE(S) CORRELATED POLARIZATION ENERGY"];
            outfile->Printf("\tSCF energy       (chkpt)              = %20.15f\n", moinfo.escf);
            outfile->Printf("\tReference energy (file100)            = %20.15f\n", moinfo.eref);
            outfile->Printf("\tPTE correlation energy                = %20.15f\n", Epte);
            outfile->Printf("\tPTE(S) correlated polarization energy = %20.15f\n", Epte_s);
            outfile->Printf("\tPTED correlated polarization energy   = %20.15f\n",
                            Process::environment.globals["PCM-CC-PTED CORRELATED POLARIZATION ENERGY"]);
            outfile->Printf("\tCCSD correlation energy               = %20.15f\n", moinfo.ecc);
            outfile->Printf("\tPCM-PTE-CCSD correlation energy       = %20.15f\n", Epte);
            outfile->Printf("\tPCM-PTE(S)-CCSD correlation energy    = %20.15f\n", Epte + Epte_s);
            outfile->Printf("\tPCM-PTED-CCSD correlation energy      = %20.15f\n",
                            Process::environment.globals["CURRENT CORRELATION ENERGY"]);
            outfile->Printf("      * PCM-PTE-CCSD total energy             = %20.15f\n", moinfo.eref + Epte);
            outfile->Printf("      * PCM-PTE(S)-CCSD total energy          = %20.15f\n", moinfo.eref + Epte + Epte_s);
            outfile->Printf("      * PCM-PTED-CCSD total energy            = %20.15f\n",
                            Process::environment.globals["CURRENT ENERGY"]);
        } /* PCM-CC stuff */

        if (!params.onepdm) {
            if (!params.aobasis && params.debug_) energy(rho_params[i]);

            kinetic(ref_wfn); /* puts kinetic energy integrals into MO basis */

            lag(rho_params[i]); /* builds the orbital lagrangian pieces, I */

            /* dpd_init(1, moinfo.nirreps, params.memory, 2, frozen.occpi, frozen.occ_sym,
               frozen.virtpi, frozen.vir_sym); */

            /*  if(moinfo.nfzc || moinfo.nfzv) {
                resort_gamma();
                resort_tei();
                } */

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

        /*  dpd_close(0); dpd_close(1); */

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
        SharedMatrix Pa(new Matrix("P alpha", Ca->colspi(), Ca->colspi()));
        SharedMatrix Pb(new Matrix("P beta", Cb->colspi(), Cb->colspi()));
        int mo_offset = 0;

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
        /*Call OEProp for each root opdm */
        std::shared_ptr<OEProp> oe(new OEProp(ref_wfn));
        if (ref_wfn->same_a_b_dens()) {
            Pa->scale(0.5);
            oe->set_Da_mo(Pa);
            Pb = Pa;
        } else {
            oe->set_Da_mo(Pa);
            oe->set_Db_mo(Pb);
        }
        oe->add("DIPOLE");
        oe->add("QUADRUPOLE");
        oe->add("MULLIKEN_CHARGES");
        oe->add("NO_OCCUPATIONS");
        std::string cc_prop_label("CC");
        if (i == 0) {  // ground state
            oe->set_title(cc_prop_label);
            oe->compute();
            // set OPDM in ref_wfn for GS only
            SharedMatrix cc_Da = oe->Da_so();
            SharedMatrix ref_Da = ref_wfn->Da();
            ref_Da->copy(cc_Da);
            if (params.nstates > 1) {
                Process::environment.globals["CC ROOT 0 DIPOLE X"] = Process::environment.globals["CC DIPOLE X"];
                Process::environment.globals["CC ROOT 0 DIPOLE Y"] = Process::environment.globals["CC DIPOLE Y"];
                Process::environment.globals["CC ROOT 0 DIPOLE Z"] = Process::environment.globals["CC DIPOLE Z"];
                Process::environment.globals["CC ROOT 0 QUADRUPOLE XX"] =
                    Process::environment.globals["CC QUADRUPOLE XX"];
                Process::environment.globals["CC ROOT 0 QUADRUPOLE XY"] =
                    Process::environment.globals["CC QUADRUPOLE XY"];
                Process::environment.globals["CC ROOT 0 QUADRUPOLE XZ"] =
                    Process::environment.globals["CC QUADRUPOLE XZ"];
                Process::environment.globals["CC ROOT 0 QUADRUPOLE YY"] =
                    Process::environment.globals["CC QUADRUPOLE YY"];
                Process::environment.globals["CC ROOT 0 QUADRUPOLE YZ"] =
                    Process::environment.globals["CC QUADRUPOLE YZ"];
                Process::environment.globals["CC ROOT 0 QUADRUPOLE ZZ"] =
                    Process::environment.globals["CC QUADRUPOLE ZZ"];
            }

            // Get the NOs/occupation numbers
            std::pair<SharedMatrix, SharedVector> NOa_pair = oe->Na_mo();
            std::pair<SharedMatrix, SharedVector> NOb_pair = NOa_pair;
            if (!ref_wfn->same_a_b_dens()) {
                SharedMatrix cc_Db = oe->Db_so();
                SharedMatrix ref_Db = ref_wfn->Db();
                ref_Db->copy(cc_Db);
                // if alpha/beta are distinct get the beta NO info
                NOb_pair = oe->Nb_mo();
            }
            // write molden files for NOs?
            if (params.write_nos) {
                MoldenWriter nowriter(ref_wfn);
                std::string mol_name = ref_wfn->molecule()->name();
                SharedMatrix NO_Ca = Matrix::doublet(Ca, NOa_pair.first, false, false);
                SharedMatrix NO_Cb = Matrix::doublet(Cb, NOb_pair.first, false, false);
                nowriter.write(mol_name + "NO.molden", NO_Ca, NO_Cb, NOa_pair.second, NOb_pair.second, NOa_pair.second,
                               NOb_pair.second, options.get_bool("MOLDEN_WITH_VIRTUAL"));
            }
        } else {
            // this should set psivars correctly for root Properties
            oe->set_title(cc_prop_label + " ROOT " + std::to_string(i));
            oe->compute();
            /*- Process::environment.globals["CC ROOT n DIPOLE X"] -*/
            /*- Process::environment.globals["CC ROOT n DIPOLE Y"] -*/
            /*- Process::environment.globals["CC ROOT n DIPOLE Z"] -*/
            /*- Process::environment.globals["CC ROOT n QUADRUPOLE XX"] -*/
            /*- Process::environment.globals["CC ROOT n QUADRUPOLE XY"] -*/
            /*- Process::environment.globals["CC ROOT n QUADRUPOLE XZ"] -*/
            /*- Process::environment.globals["CC ROOT n QUADRUPOLE YY"] -*/
            /*- Process::environment.globals["CC ROOT n QUADRUPOLE YZ"] -*/
            /*- Process::environment.globals["CC ROOT n QUADRUPOLE ZZ"] -*/
        }

        free_block(moinfo.opdm);

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

    /*
      if ( params.ael && (params.nstates > 1) )
        ael(rho_params);
    */
    // outfile->Printf("I am here\n");

    if (params.transition) {
        MintsHelper mints(ref_wfn->basisset(), options, 0);
        get_td_params(options);
        for (i = 0; i < params.nstates; i++) {
            td_setup(td_params[i]);
            tdensity(td_params[i]);
            outfile->Printf("Doing transition\n");
            oscillator_strength(ref_wfn, &(td_params[i]));
            outfile->Printf("Doing transition\n");
            if (params.ref == 0) {
                rotational_strength(mints, &(td_params[i]));
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
                    ex_oscillator_strength(ref_wfn, &(td_params[state1]), &(td_params[state2]), &xtd_data);
                    if (params.ref == 0) {
                        // ex_rotational_strength(&(td_params[j]),&(td_params[i+1]), &xtd_data);
                        ex_rotational_strength(mints, &(td_params[state1]), &(td_params[state2]), &xtd_data);
                    }

                    xtd_params.push_back(xtd_data);

                    td_cleanup();
                }
            }
            td_print();
            ex_td_print(xtd_params);
        }

    }  // End params.transition IF loop

    // outfile->Printf("I am here\n");
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
void init_io(void) {
    int i, num_unparsed;
    char *argv_unparsed[100];
    ;

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

    tstart();

    for (i = PSIF_CC_MIN; i <= PSIF_CC_MAX; i++) psio_open(i, PSIO_OPEN_OLD);
    // erase old files
    psio_close(PSIF_CC_GR, 0);
    psio_close(PSIF_CC_GL, 0);
    psio_close(PSIF_EOM_TMP0, 0);
    psio_open(PSIF_CC_GR, PSIO_OPEN_NEW);
    psio_open(PSIF_CC_GL, PSIO_OPEN_NEW);
    psio_open(PSIF_EOM_TMP0, PSIO_OPEN_NEW);
}

void title(void) {
    outfile->Printf("\n");
    outfile->Printf("\t\t\t**************************\n");
    outfile->Printf("\t\t\t*                        *\n");
    outfile->Printf("\t\t\t*        CCDENSITY       *\n");
    outfile->Printf("\t\t\t*                        *\n");
    outfile->Printf("\t\t\t**************************\n");
    outfile->Printf("\n");
}

void exit_io(void) {
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

    tstop();
}

void update_F_pcm_rhf(SharedMatrix &MO_PCM_potential) {
    dpdfile2 fIJ, fij, fAB, fab, fIA, fia;
    if (!(_default_psio_lib_->tocentry_exists(PSIF_CC_OEI,
                                              "fIJ original")))  // We check just one of the six blocks for existence
    {
        // If it doesn't exist load the various blocks...
        global_dpd_->file2_init(&fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ");
        global_dpd_->file2_init(&fij, PSIF_CC_OEI, 0, 0, 0, "fij");
        global_dpd_->file2_init(&fAB, PSIF_CC_OEI, 0, 1, 1, "fAB");
        global_dpd_->file2_init(&fab, PSIF_CC_OEI, 0, 1, 1, "fab");
        global_dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
        global_dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 0, 1, "fia");
        // ...copy them in the "original" file...
        global_dpd_->file2_copy(&fIJ, PSIF_CC_OEI, "fIJ original");
        global_dpd_->file2_copy(&fij, PSIF_CC_OEI, "fij original");
        global_dpd_->file2_copy(&fAB, PSIF_CC_OEI, "fAB original");
        global_dpd_->file2_copy(&fab, PSIF_CC_OEI, "fab original");
        global_dpd_->file2_copy(&fIA, PSIF_CC_OEI, "fIA original");
        global_dpd_->file2_copy(&fia, PSIF_CC_OEI, "fia original");
        // ...and close.
        global_dpd_->file2_close(&fIJ);
        global_dpd_->file2_close(&fij);
        global_dpd_->file2_close(&fAB);
        global_dpd_->file2_close(&fab);
        global_dpd_->file2_close(&fIA);
        global_dpd_->file2_close(&fia);
    }
    // Read from the original Fock matrix (the one coming from the PCM-SCF step)
    global_dpd_->file2_init(&fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ original");
    global_dpd_->file2_mat_init(&fIJ);
    global_dpd_->file2_mat_rd(&fIJ);
    global_dpd_->file2_init(&fij, PSIF_CC_OEI, 0, 0, 0, "fij original");
    global_dpd_->file2_mat_init(&fij);
    global_dpd_->file2_mat_rd(&fij);
    global_dpd_->file2_init(&fAB, PSIF_CC_OEI, 0, 1, 1, "fAB original");
    global_dpd_->file2_mat_init(&fAB);
    global_dpd_->file2_mat_rd(&fAB);
    global_dpd_->file2_init(&fab, PSIF_CC_OEI, 0, 1, 1, "fab original");
    global_dpd_->file2_mat_init(&fab);
    global_dpd_->file2_mat_rd(&fab);
    global_dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA original");
    global_dpd_->file2_mat_init(&fIA);
    global_dpd_->file2_mat_rd(&fIA);
    global_dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 0, 1, "fia original");
    global_dpd_->file2_mat_init(&fia);
    global_dpd_->file2_mat_rd(&fia);
    // Open a series of buffers to put the various blocks in
    dpdfile2 fIJ2, fij2, fAB2, fab2, fIA2, fia2;
    global_dpd_->file2_init(&fIJ2, PSIF_CC_OEI, 0, 0, 0, "fIJ");
    global_dpd_->file2_mat_init(&fIJ2);
    global_dpd_->file2_mat_rd(&fIJ2);
    global_dpd_->file2_init(&fij2, PSIF_CC_OEI, 0, 0, 0, "fij");
    global_dpd_->file2_mat_init(&fij2);
    global_dpd_->file2_mat_rd(&fij2);
    global_dpd_->file2_init(&fAB2, PSIF_CC_OEI, 0, 1, 1, "fAB");
    global_dpd_->file2_mat_init(&fAB2);
    global_dpd_->file2_mat_rd(&fAB2);
    global_dpd_->file2_init(&fab2, PSIF_CC_OEI, 0, 1, 1, "fab");
    global_dpd_->file2_mat_init(&fab2);
    global_dpd_->file2_mat_rd(&fab2);
    global_dpd_->file2_init(&fIA2, PSIF_CC_OEI, 0, 0, 1, "fIA");
    global_dpd_->file2_mat_init(&fIA2);
    global_dpd_->file2_mat_rd(&fIA2);
    global_dpd_->file2_init(&fia2, PSIF_CC_OEI, 0, 0, 1, "fia");
    global_dpd_->file2_mat_init(&fia2);
    global_dpd_->file2_mat_rd(&fia2);

    // Add MO_PCM_potential on top of the original Fock matrix
    for (int h = 0; h < moinfo.nirreps; ++h) {
        int upper_bound_ij = moinfo.occpi[h] - moinfo.openpi[h];
        int upper_bound_IA = moinfo.virtpi[h] - moinfo.openpi[h];
        for (int i = 0; i < moinfo.occpi[h]; ++i) {
            int I = i + moinfo.frdocc[h];
            for (int j = 0; j < moinfo.occpi[h]; ++j) {
                int J = j + moinfo.frdocc[h];
                fIJ2.matrix[h][i][j] = fIJ.matrix[h][i][j] + MO_PCM_potential->get(h, I, J);
            }
        }
        for (int i = 0; i < upper_bound_ij; ++i) {
            int I = i + moinfo.frdocc[h];
            for (int j = 0; j < upper_bound_ij; ++j) {
                int J = j + moinfo.frdocc[h];
                fij2.matrix[h][i][j] = fij.matrix[h][i][j] + MO_PCM_potential->get(h, I, J);
            }
        }
        for (int i = 0; i < moinfo.occpi[h]; ++i) {
            int I = i + moinfo.frdocc[h];
            for (int a = 0; a < upper_bound_IA; ++a) {
                int A = a + moinfo.frdocc[h] + moinfo.occpi[h];
                fIA2.matrix[h][i][a] = fIA.matrix[h][i][a] + MO_PCM_potential->get(h, I, A);
            }
        }
        for (int i = 0; i < upper_bound_ij; ++i) {
            int I = i + moinfo.frdocc[h];
            for (int a = 0; a < upper_bound_IA; ++a) {
                int A = a + moinfo.frdocc[h] + moinfo.occpi[h];
                fia2.matrix[h][i][a] = fia.matrix[h][i][a] + MO_PCM_potential->get(h, I, A);
            }
        }
        for (int i = 0; i < upper_bound_ij; ++i) {
            int I = i + moinfo.frdocc[h];
            for (int a = 0; a < moinfo.openpi[h]; ++a) {
                int A = a + upper_bound_IA;
                int aA = a + moinfo.frdocc[h] + upper_bound_ij;
                fia2.matrix[h][i][A] = fia.matrix[h][i][A] + MO_PCM_potential->get(h, I, aA);
            }
        }
        for (int a = 0; a < upper_bound_IA; ++a) {
            int A = a + moinfo.frdocc[h] + moinfo.occpi[h];
            for (int b = 0; b < upper_bound_IA; ++b) {
                int B = b + moinfo.frdocc[h] + moinfo.occpi[h];
                fAB2.matrix[h][a][b] = fAB.matrix[h][a][b] + MO_PCM_potential->get(h, A, B);
                fab2.matrix[h][a][b] = fab.matrix[h][a][b] + MO_PCM_potential->get(h, A, B);
            }
        }
        for (int a = 0; a < moinfo.openpi[h]; ++a) {
            int A = a + upper_bound_IA;
            int aA = a + moinfo.frdocc[h] + upper_bound_ij;
            for (int b = 0; b < moinfo.openpi[h]; ++b) {
                int B = b + upper_bound_IA;
                int bB = b + moinfo.frdocc[h] + upper_bound_ij;
                fab2.matrix[h][A][B] = fab.matrix[h][A][B] + MO_PCM_potential->get(h, aA, bB);
            }
        }
        for (int a = 0; a < upper_bound_IA; ++a) {
            int A = a + moinfo.frdocc[h] + moinfo.occpi[h];
            for (int b = 0; b < moinfo.openpi[h]; ++b) {
                int B = b + upper_bound_IA;
                int bB = b + moinfo.frdocc[h] + upper_bound_ij;
                fab2.matrix[h][a][B] = fab.matrix[h][a][B] + MO_PCM_potential->get(h, A, bB);
            }
        }
        for (int a = 0; a < moinfo.openpi[h]; ++a) {
            int A = a + upper_bound_IA;
            int aA = a + moinfo.frdocc[h] + upper_bound_ij;
            for (int b = 0; b < upper_bound_IA; ++b) {
                int B = b + moinfo.frdocc[h] + moinfo.occpi[h];
                fab2.matrix[h][A][b] = fab.matrix[h][A][b] + MO_PCM_potential->get(h, aA, B);
            }
        }
    }
    // Write the updated Fock matrix
    global_dpd_->file2_mat_wrt(&fIJ2);
    global_dpd_->file2_mat_wrt(&fij2);
    global_dpd_->file2_mat_wrt(&fAB2);
    global_dpd_->file2_mat_wrt(&fab2);
    global_dpd_->file2_mat_wrt(&fIA2);
    global_dpd_->file2_mat_wrt(&fia2);

    global_dpd_->file2_mat_print(&fIA2, "outfile");
    global_dpd_->file2_mat_print(&fAB2, "outfile");
    global_dpd_->file2_mat_print(&fIJ2, "outfile");
    // Close the original Fock matrix
    global_dpd_->file2_close(&fIJ);
    global_dpd_->file2_close(&fij);
    global_dpd_->file2_close(&fAB);
    global_dpd_->file2_close(&fab);
    global_dpd_->file2_close(&fIA);
    global_dpd_->file2_close(&fia);
    // Close the updated Fock matrix
    global_dpd_->file2_close(&fIJ2);
    global_dpd_->file2_close(&fij2);
    global_dpd_->file2_close(&fAB2);
    global_dpd_->file2_close(&fab2);
    global_dpd_->file2_close(&fIA2);
    global_dpd_->file2_close(&fia2);
}

void update_F_pcm_uhf(SharedMatrix &MO_PCM_potential_A, SharedMatrix &MO_PCM_potential_B) {
    dpdfile2 fIJ, fij, fAB, fab, fIA, fia;
    if (!(_default_psio_lib_->tocentry_exists(PSIF_CC_OEI,
                                              "fIJ original")))  // We check just one of the six blocks for existence
    {
        // If it doesn't exist load the various blocks...
        global_dpd_->file2_init(&fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ");
        global_dpd_->file2_init(&fij, PSIF_CC_OEI, 0, 2, 2, "fij");
        global_dpd_->file2_init(&fAB, PSIF_CC_OEI, 0, 1, 1, "fAB");
        global_dpd_->file2_init(&fab, PSIF_CC_OEI, 0, 3, 3, "fab");
        global_dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
        global_dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 2, 3, "fia");
        // ...copy them in the "original" file...
        global_dpd_->file2_copy(&fIJ, PSIF_CC_OEI, "fIJ original");
        global_dpd_->file2_copy(&fij, PSIF_CC_OEI, "fij original");
        global_dpd_->file2_copy(&fAB, PSIF_CC_OEI, "fAB original");
        global_dpd_->file2_copy(&fab, PSIF_CC_OEI, "fab original");
        global_dpd_->file2_copy(&fIA, PSIF_CC_OEI, "fIA original");
        global_dpd_->file2_copy(&fia, PSIF_CC_OEI, "fia original");
        // ...and close.
        global_dpd_->file2_close(&fIJ);
        global_dpd_->file2_close(&fij);
        global_dpd_->file2_close(&fAB);
        global_dpd_->file2_close(&fab);
        global_dpd_->file2_close(&fIA);
        global_dpd_->file2_close(&fia);
    }
    // Read from the original Fock matrix (the one coming from the PCM-SCF step)
    global_dpd_->file2_init(&fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ original");
    global_dpd_->file2_mat_init(&fIJ);
    global_dpd_->file2_mat_rd(&fIJ);
    global_dpd_->file2_init(&fij, PSIF_CC_OEI, 0, 2, 2, "fij original");
    global_dpd_->file2_mat_init(&fij);
    global_dpd_->file2_mat_rd(&fij);
    global_dpd_->file2_init(&fAB, PSIF_CC_OEI, 0, 1, 1, "fAB original");
    global_dpd_->file2_mat_init(&fAB);
    global_dpd_->file2_mat_rd(&fAB);
    global_dpd_->file2_init(&fab, PSIF_CC_OEI, 0, 3, 3, "fab original");
    global_dpd_->file2_mat_init(&fab);
    global_dpd_->file2_mat_rd(&fab);
    global_dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA original");
    global_dpd_->file2_mat_init(&fIA);
    global_dpd_->file2_mat_rd(&fIA);
    global_dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 2, 3, "fia original");
    global_dpd_->file2_mat_init(&fia);
    global_dpd_->file2_mat_rd(&fia);
    // Open a series of buffers to put the various blocks in
    dpdfile2 fIJ2, fij2, fAB2, fab2, fIA2, fia2;
    global_dpd_->file2_init(&fIJ2, PSIF_CC_OEI, 0, 0, 0, "fIJ");
    global_dpd_->file2_mat_init(&fIJ2);
    global_dpd_->file2_mat_rd(&fIJ2);
    global_dpd_->file2_init(&fij2, PSIF_CC_OEI, 0, 2, 2, "fij");
    global_dpd_->file2_mat_init(&fij2);
    global_dpd_->file2_mat_rd(&fij2);
    global_dpd_->file2_init(&fAB2, PSIF_CC_OEI, 0, 1, 1, "fAB");
    global_dpd_->file2_mat_init(&fAB2);
    global_dpd_->file2_mat_rd(&fAB2);
    global_dpd_->file2_init(&fab2, PSIF_CC_OEI, 0, 3, 3, "fab");
    global_dpd_->file2_mat_init(&fab2);
    global_dpd_->file2_mat_rd(&fab2);
    global_dpd_->file2_init(&fIA2, PSIF_CC_OEI, 0, 0, 1, "fIA");
    global_dpd_->file2_mat_init(&fIA2);
    global_dpd_->file2_mat_rd(&fIA2);
    global_dpd_->file2_init(&fia2, PSIF_CC_OEI, 0, 2, 3, "fia");
    global_dpd_->file2_mat_init(&fia2);
    global_dpd_->file2_mat_rd(&fia2);
    // Add MO_PCM_potential on top of the original Fock matrix
    for (int h = 0; h < moinfo.nirreps; ++h) {
        for (int i = 0; i < moinfo.aoccpi[h]; ++i) {
            int I = i + moinfo.frdocc[h];
            for (int j = 0; j < moinfo.aoccpi[h]; ++j) {
                int J = j + moinfo.frdocc[h];
                fIJ2.matrix[h][i][j] = fIJ.matrix[h][i][j] + MO_PCM_potential_A->get(h, I, J);
            }
        }
        for (int i = 0; i < moinfo.boccpi[h]; ++i) {
            int I = i + moinfo.frdocc[h];
            for (int j = 0; j < moinfo.boccpi[h]; ++j) {
                int J = j + moinfo.frdocc[h];
                fij2.matrix[h][i][j] = fij.matrix[h][i][j] + MO_PCM_potential_B->get(h, I, J);
            }
        }
        for (int i = 0; i < moinfo.aoccpi[h]; ++i) {
            int I = i + moinfo.frdocc[h];
            for (int a = 0; a < moinfo.avirtpi[h]; ++a) {
                int A = a + moinfo.frdocc[h] + moinfo.aoccpi[h];
                fIA2.matrix[h][i][a] = fIA.matrix[h][i][a] + MO_PCM_potential_A->get(h, I, A);
            }
        }
        for (int i = 0; i < moinfo.boccpi[h]; ++i) {
            int I = i + moinfo.frdocc[h];
            for (int a = 0; a < moinfo.bvirtpi[h]; ++a) {
                int A = a + moinfo.frdocc[h] + moinfo.boccpi[h];
                fia2.matrix[h][i][a] = fia.matrix[h][i][a] + MO_PCM_potential_B->get(h, I, A);
            }
        }
        for (int a = 0; a < moinfo.avirtpi[h]; ++a) {
            int A = a + moinfo.frdocc[h] + moinfo.aoccpi[h];
            for (int b = 0; b < moinfo.avirtpi[h]; ++b) {
                int B = b + moinfo.frdocc[h] + moinfo.aoccpi[h];
                fAB2.matrix[h][a][b] = fAB.matrix[h][a][b] + MO_PCM_potential_A->get(h, A, B);
            }
        }
        for (int a = 0; a < moinfo.bvirtpi[h]; ++a) {
            int A = a + moinfo.frdocc[h] + moinfo.boccpi[h];
            for (int b = 0; b < moinfo.bvirtpi[h]; ++b) {
                int B = b + moinfo.frdocc[h] + moinfo.boccpi[h];
                fab2.matrix[h][a][b] = fab.matrix[h][a][b] + MO_PCM_potential_B->get(h, A, B);
            }
        }
    }
    global_dpd_->file2_mat_wrt(&fIJ2);
    global_dpd_->file2_mat_wrt(&fij2);
    global_dpd_->file2_mat_wrt(&fAB2);
    global_dpd_->file2_mat_wrt(&fab2);
    global_dpd_->file2_mat_wrt(&fIA2);
    global_dpd_->file2_mat_wrt(&fia2);

    global_dpd_->file2_close(&fIJ);
    global_dpd_->file2_close(&fij);
    global_dpd_->file2_close(&fAB);
    global_dpd_->file2_close(&fab);
    global_dpd_->file2_close(&fIA);
    global_dpd_->file2_close(&fia);
    global_dpd_->file2_close(&fIJ2);
    global_dpd_->file2_close(&fij2);
    global_dpd_->file2_close(&fAB2);
    global_dpd_->file2_close(&fab2);
    global_dpd_->file2_close(&fIA2);
    global_dpd_->file2_close(&fia2);
}
}
}  // namespace psi::ccdensity

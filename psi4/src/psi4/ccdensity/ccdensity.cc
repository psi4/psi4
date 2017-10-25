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
void update_Fock_matrix_rhf(SharedMatrix &);
void update_Fock_matrix_uhf(SharedMatrix &, SharedMatrix &);
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

        // Prepare density matrices for PCM-CC and OEprop
        SharedMatrix Ca = ref_wfn->Ca();
        SharedMatrix Cb = ref_wfn->Cb();

        Dimension nmopi = ref_wfn->nmopi();
        Dimension frzvpi = ref_wfn->frzvpi();

        // MO basis alpha and beta correlated OPDMs
        SharedMatrix MO_OPDM_a =
            std::make_shared<Matrix>("Alpha density matrix in MO basis", Ca->rowspi(), Ca->colspi());
        SharedMatrix MO_OPDM_b =
            std::make_shared<Matrix>("Beta density matrix in MO basis", Cb->rowspi(), Cb->colspi());

        int mo_offset = 0;
        for (int h = 0; h < Ca->nirrep(); h++) {
            int nmo = nmopi[h];
            int nfv = frzvpi[h];
            int nmor = nmo - nfv;
            if (!nmo || !nmor) continue;

            // Loop over QT, convert to Pitzer
            double **MO_OPDM_ap = MO_OPDM_a->pointer(h);
            double **MO_OPDM_bp = MO_OPDM_b->pointer(h);
            for (int i = 0; i < nmor; i++) {
                for (int j = 0; j < nmor; j++) {
                    int I = moinfo.pitzer2qt[i + mo_offset];
                    int J = moinfo.pitzer2qt[j + mo_offset];
                    if (ref_wfn->same_a_b_dens())
                        MO_OPDM_ap[i][j] = moinfo.opdm[I][J];
                    else {
                        MO_OPDM_ap[i][j] = moinfo.opdm_a[I][J];
                        MO_OPDM_bp[i][j] = moinfo.opdm_b[I][J];
                    }
                }
            }
            mo_offset += nmo;
        }

        // SO basis alpha and beta correlated OPDMs
        SharedMatrix SO_OPDM_a =
            std::make_shared<Matrix>("Alpha density matrix in SO basis", ref_wfn->nsopi(), ref_wfn->nsopi());
        SO_OPDM_a->back_transform(MO_OPDM_a, Ca);
        SharedMatrix SO_OPDM_b =
            std::make_shared<Matrix>("Beta density matrix in SO basis", ref_wfn->nsopi(), ref_wfn->nsopi());
        SO_OPDM_b->back_transform(MO_OPDM_b, Cb);

        /* PCM-CC stuff */
        if (options.get_bool("PCM") &&
            (options.get_str("PCM_CC_TYPE") == "PTED")) {  // PTED scheme: R. Cammi, JCP, 131, 164104 (2009)
            //   Build the PCM polarization charges from the CC density matrix at the current iteration.
            //   The polarization charges are then contracted with the electrostatic potential integrals
            //   to build the PCM potential to be added in the CC T and Lambda equations.
            SharedPCM cc_pcm = std::make_shared<PCM>(options, ref_wfn->psio(), moinfo.nirreps, ref_wfn->basisset());
            if (params.ref == 0 || params.ref == 1) { /* RHF or ROHF */
                SharedMatrix MO_PCM_potential =
                    std::make_shared<Matrix>("PCM potential in MO basis", Ca->rowspi(), Ca->colspi());
                // We have to get the 0.5 * (QV) energy contribution
                double Epcm_correlated = cc_pcm->compute_E(SO_OPDM_a, psi::PCM::EleOnly);
                Process::environment.globals["PCM-CC-PTED CORRELATED POLARIZATION ENERGY"] = Epcm_correlated;
                double E_correlation = Process::environment.globals["CURRENT CORRELATION ENERGY"];
                E_correlation -= Epcm_correlated;
                Process::environment.globals["CURRENT CORRELATION ENERGY"] = E_correlation;
                double E = Process::environment.globals["CURRENT ENERGY"];
                E -= Epcm_correlated;
                Process::environment.globals["CURRENT ENERGY"] = E;
                SharedMatrix SO_PCM_potential = cc_pcm->compute_V_electronic();  // This is in SO basis
                // We now transform it to MO basis...
                MO_PCM_potential->transform(SO_PCM_potential, Ca);
                update_Fock_matrix_rhf(MO_PCM_potential);
            } else if (params.ref == 2) { /* UHF case */
                SharedMatrix MO_PCM_potential_a =
                    std::make_shared<Matrix>("Alpha PCM potential in MO basis", Ca->rowspi(), Ca->colspi());
                SharedMatrix MO_PCM_potential_b =
                    std::make_shared<Matrix>("Beta PCM potential in MO basis", Cb->rowspi(), Cb->colspi());

                SharedMatrix SO_OPDM(SO_OPDM_a->clone());
                SO_OPDM->add(SO_OPDM_b);
                // We have to get the 0.5 * (QV) energy contribution
                double Epcm_correlated = cc_pcm->compute_E(SO_OPDM, PCM::EleOnly);
                Process::environment.globals["PCM-CC-PTED CORRELATED POLARIZATION ENERGY"] = Epcm_correlated;
                double E_correlation = Process::environment.globals["CURRENT CORRELATION ENERGY"];
                E_correlation -= Epcm_correlated;
                Process::environment.globals["CURRENT CORRELATION ENERGY"] = E_correlation;
                double E = Process::environment.globals["CURRENT ENERGY"];
                E -= Epcm_correlated;
                Process::environment.globals["CURRENT ENERGY"] = E;
                SharedMatrix SO_PCM_potential = cc_pcm->compute_V_electronic();  // This is in SO basis
                // We now transform it to MO basis...
                MO_PCM_potential_a->transform(SO_PCM_potential, Ca);
                MO_PCM_potential_b->transform(SO_PCM_potential, Cb);
                update_Fock_matrix_uhf(MO_PCM_potential_a, MO_PCM_potential_b);
            }
            double Epte = Process::environment.globals["PCM-CC-PTE CORRELATION ENERGY"];
            outfile->Printf("\tSCF energy       (chkpt)              = %20.15f\n", moinfo.escf);
            outfile->Printf("\tReference energy (file100)            = %20.15f\n", moinfo.eref);
            outfile->Printf("\tPTE correlation energy                = %20.15f\n", Epte);
            outfile->Printf("\tPTED correlated polarization energy   = %20.15f\n",
                            Process::environment.globals["PCM-CC-PTED CORRELATED POLARIZATION ENERGY"]);
            outfile->Printf("\tCCSD correlation energy               = %20.15f\n", moinfo.ecc);
            outfile->Printf("\tPCM-PTE-CCSD correlation energy       = %20.15f\n", Epte);
            outfile->Printf("\tPCM-PTED-CCSD correlation energy      = %20.15f\n",
                            Process::environment.globals["CURRENT CORRELATION ENERGY"]);
            outfile->Printf("      * PCM-PTE-CCSD total energy             = %20.15f\n", moinfo.eref + Epte);
            outfile->Printf("      * PCM-PTED-CCSD total energy            = %20.15f\n",
                            Process::environment.globals["CURRENT ENERGY"]);
        } /* PCM-CC stuff */

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

        // Update contents of MO_OPDM_a and MO_OPDM_b
        // with the total (reference + correlated) OPDM
        mo_offset = 0;
        for (int h = 0; h < Ca->nirrep(); h++) {
            int nmo = nmopi[h];
            int nfv = frzvpi[h];
            int nmor = nmo - nfv;
            if (!nmo || !nmor) continue;

            // Loop over QT, convert to Pitzer
            double **MO_OPDM_ap = MO_OPDM_a->pointer(h);
            double **MO_OPDM_bp = MO_OPDM_b->pointer(h);
            for (int i = 0; i < nmor; i++) {
                for (int j = 0; j < nmor; j++) {
                    int I = moinfo.pitzer2qt[i + mo_offset];
                    int J = moinfo.pitzer2qt[j + mo_offset];
                    if (ref_wfn->same_a_b_dens())
                        MO_OPDM_ap[i][j] = moinfo.opdm[I][J];
                    else {
                        MO_OPDM_ap[i][j] = moinfo.opdm_a[I][J];
                        MO_OPDM_bp[i][j] = moinfo.opdm_b[I][J];
                    }
                }
            }
            mo_offset += nmo;
        }
        /*Call OEProp for each root opdm */
        std::shared_ptr<OEProp> oe = std::make_shared<OEProp>(ref_wfn);
        if (ref_wfn->same_a_b_dens()) {
            MO_OPDM_a->scale(0.5);
            oe->set_Da_mo(MO_OPDM_a);
            MO_OPDM_b = MO_OPDM_a;
        } else {
            oe->set_Da_mo(MO_OPDM_a);
            oe->set_Db_mo(MO_OPDM_b);
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

        if (ref_wfn->same_a_b_dens())
            free_block(moinfo.opdm);
        else {
            free_block(moinfo.opdm_a);
            free_block(moinfo.opdm_b);
        }

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

void update_Fock_matrix_rhf(SharedMatrix &MO_OE_operator) {
    dpdfile2 dpd_fIJ, dpd_fij, dpd_fAB, dpd_fab, dpd_fIA, dpd_fia;
    // We save the Fock matrix from the SCF step under a "original" label.
    // The MO_OE_operator will always be added on top of that.
    if (!(_default_psio_lib_->tocentry_exists(PSIF_CC_OEI, "fIJ original"))) {
        global_dpd_->file2_init(&dpd_fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ");
        global_dpd_->file2_init(&dpd_fij, PSIF_CC_OEI, 0, 0, 0, "fij");
        global_dpd_->file2_init(&dpd_fAB, PSIF_CC_OEI, 0, 1, 1, "fAB");
        global_dpd_->file2_init(&dpd_fab, PSIF_CC_OEI, 0, 1, 1, "fab");
        global_dpd_->file2_init(&dpd_fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
        global_dpd_->file2_init(&dpd_fia, PSIF_CC_OEI, 0, 0, 1, "fia");
        global_dpd_->file2_copy(&dpd_fIJ, PSIF_CC_OEI, "fIJ original");
        global_dpd_->file2_copy(&dpd_fij, PSIF_CC_OEI, "fij original");
        global_dpd_->file2_copy(&dpd_fAB, PSIF_CC_OEI, "fAB original");
        global_dpd_->file2_copy(&dpd_fab, PSIF_CC_OEI, "fab original");
        global_dpd_->file2_copy(&dpd_fIA, PSIF_CC_OEI, "fIA original");
        global_dpd_->file2_copy(&dpd_fia, PSIF_CC_OEI, "fia original");
        global_dpd_->file2_close(&dpd_fIJ);
        global_dpd_->file2_close(&dpd_fij);
        global_dpd_->file2_close(&dpd_fAB);
        global_dpd_->file2_close(&dpd_fab);
        global_dpd_->file2_close(&dpd_fIA);
        global_dpd_->file2_close(&dpd_fia);
    }
    // Generate Dimensions for slices
    Dimension zero(moinfo.nirreps);

    std::vector<int> tmp(moinfo.nirreps);
    tmp.assign(moinfo.frdocc, moinfo.frdocc + moinfo.nirreps);
    Dimension frdoccpi(tmp);

    tmp.assign(moinfo.occpi, moinfo.occpi + moinfo.nirreps);
    Dimension occpi(tmp);

    tmp.assign(moinfo.openpi, moinfo.openpi + moinfo.nirreps);
    Dimension openpi(tmp);

    tmp.assign(moinfo.virtpi, moinfo.virtpi + moinfo.nirreps);
    Dimension virtpi(tmp);

    /* Only defined for UHF
    tmp.assign(moinfo.avirtpi, moinfo.avirtpi + moinfo.nirreps);
    Dimension avirtpi(tmp);

    tmp.assign(moinfo.bvirtpi, moinfo.bvirtpi + moinfo.nirreps);
    Dimension bvirtpi(tmp);
    */

    global_dpd_->file2_init(&dpd_fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ original");
    SharedMatrix fIJ = std::make_shared<Matrix>(&dpd_fIJ);
    global_dpd_->file2_close(&dpd_fIJ);
    fIJ->set_name("fIJ");
    fIJ->add(MO_OE_operator->get_block({frdoccpi, frdoccpi + occpi}, {frdoccpi, frdoccpi + occpi}));
    global_dpd_->file2_init(&dpd_fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ");
    fIJ->write_to_dpdfile2(&dpd_fIJ);
    global_dpd_->file2_close(&dpd_fIJ);

    global_dpd_->file2_init(&dpd_fij, PSIF_CC_OEI, 0, 0, 0, "fij original");
    SharedMatrix fij = std::make_shared<Matrix>(&dpd_fij);
    global_dpd_->file2_close(&dpd_fij);
    fij->set_name("fij");
    fij->add(
        MO_OE_operator->get_block({frdoccpi, frdoccpi + (occpi - openpi)}, {frdoccpi, frdoccpi + (occpi - openpi)}));
    global_dpd_->file2_init(&dpd_fij, PSIF_CC_OEI, 0, 0, 0, "fij");
    fij->write_to_dpdfile2(&dpd_fij);
    global_dpd_->file2_close(&dpd_fij);

    global_dpd_->file2_init(&dpd_fAB, PSIF_CC_OEI, 0, 1, 1, "fAB original");
    SharedMatrix fAB = std::make_shared<Matrix>(&dpd_fAB);
    global_dpd_->file2_close(&dpd_fAB);
    fAB->set_name("fAB");
    fAB->add(MO_OE_operator->get_block({(frdoccpi + occpi), (frdoccpi + occpi) + (virtpi - openpi)},
                                       {(frdoccpi + occpi), (frdoccpi + occpi) + (virtpi - openpi)}));
    global_dpd_->file2_init(&dpd_fAB, PSIF_CC_OEI, 0, 1, 1, "fAB");
    fAB->write_to_dpdfile2(&dpd_fAB);
    global_dpd_->file2_close(&dpd_fAB);

    global_dpd_->file2_init(&dpd_fab, PSIF_CC_OEI, 0, 1, 1, "fab original");
    SharedMatrix fab = std::make_shared<Matrix>(&dpd_fab);
    global_dpd_->file2_close(&dpd_fab);
    fab->set_name("fab");
    fab->add(MO_OE_operator->get_block({(frdoccpi + occpi), (frdoccpi + occpi) + (virtpi - openpi)},
                                       {(frdoccpi + occpi), (frdoccpi + occpi) + (virtpi - openpi)}));
    (fab->get_block({virtpi - openpi, virtpi}, {virtpi - openpi, virtpi}))
        ->add(MO_OE_operator->get_block({frdoccpi + occpi - openpi, frdoccpi + occpi},
                                        {frdoccpi + occpi - openpi, frdoccpi + occpi}));
    (fab->get_block({zero, virtpi - openpi}, {virtpi - openpi, virtpi}))
        ->add(MO_OE_operator->get_block({frdoccpi + occpi, (frdoccpi + occpi + (virtpi - openpi))},
                                        {frdoccpi + occpi - openpi, frdoccpi + occpi}));
    (fab->get_block({virtpi - openpi, virtpi}, {zero, virtpi - openpi}))
        ->add(MO_OE_operator->get_block({frdoccpi + occpi - openpi, frdoccpi + occpi},
                                        {frdoccpi + occpi, frdoccpi + occpi + (virtpi - openpi)}));
    global_dpd_->file2_init(&dpd_fab, PSIF_CC_OEI, 0, 1, 1, "fab");
    fab->write_to_dpdfile2(&dpd_fab);
    global_dpd_->file2_close(&dpd_fab);

    global_dpd_->file2_init(&dpd_fIA, PSIF_CC_OEI, 0, 0, 1, "fIA original");
    SharedMatrix fIA = std::make_shared<Matrix>(&dpd_fIA);
    global_dpd_->file2_close(&dpd_fIA);
    fIA->set_name("fIA");
    fIA->add(MO_OE_operator->get_block({frdoccpi, frdoccpi + occpi},
                                       {frdoccpi + occpi, frdoccpi + occpi + (virtpi - openpi)}));
    global_dpd_->file2_init(&dpd_fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
    fIA->write_to_dpdfile2(&dpd_fIA);
    global_dpd_->file2_close(&dpd_fIA);

    global_dpd_->file2_init(&dpd_fia, PSIF_CC_OEI, 0, 0, 1, "fia original");
    SharedMatrix fia = std::make_shared<Matrix>(&dpd_fia);
    global_dpd_->file2_close(&dpd_fia);
    fia->set_name("fia");
    fia->add(MO_OE_operator->get_block({frdoccpi, frdoccpi + (occpi - openpi)},
                                       {frdoccpi + occpi, frdoccpi + occpi + (virtpi - openpi)}));
    (fia->get_block({zero, occpi - openpi}, {virtpi - openpi, virtpi}))
        ->add(MO_OE_operator->get_block({frdoccpi, frdoccpi + (occpi - openpi)},
                                        {frdoccpi + (occpi - openpi), frdoccpi + occpi}));
    global_dpd_->file2_init(&dpd_fia, PSIF_CC_OEI, 0, 0, 1, "fia");
    fia->write_to_dpdfile2(&dpd_fia);
    global_dpd_->file2_close(&dpd_fia);
}

void update_Fock_matrix_uhf(SharedMatrix &MO_OE_operator_a, SharedMatrix &MO_OE_operator_b) {
    dpdfile2 dpd_fIJ, dpd_fij, dpd_fAB, dpd_fab, dpd_fIA, dpd_fia;
    // We save the Fock matrix from the SCF step under a "original" label.
    // The MO_OE_operator will always be added on top of that.
    if (!(_default_psio_lib_->tocentry_exists(PSIF_CC_OEI, "fIJ original"))) {
        global_dpd_->file2_init(&dpd_fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ");
        global_dpd_->file2_init(&dpd_fij, PSIF_CC_OEI, 0, 2, 2, "fij");
        global_dpd_->file2_init(&dpd_fAB, PSIF_CC_OEI, 0, 1, 1, "fAB");
        global_dpd_->file2_init(&dpd_fab, PSIF_CC_OEI, 0, 3, 3, "fab");
        global_dpd_->file2_init(&dpd_fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
        global_dpd_->file2_init(&dpd_fia, PSIF_CC_OEI, 0, 2, 3, "fia");
        global_dpd_->file2_copy(&dpd_fIJ, PSIF_CC_OEI, "fIJ original");
        global_dpd_->file2_copy(&dpd_fij, PSIF_CC_OEI, "fij original");
        global_dpd_->file2_copy(&dpd_fAB, PSIF_CC_OEI, "fAB original");
        global_dpd_->file2_copy(&dpd_fab, PSIF_CC_OEI, "fab original");
        global_dpd_->file2_copy(&dpd_fIA, PSIF_CC_OEI, "fIA original");
        global_dpd_->file2_copy(&dpd_fia, PSIF_CC_OEI, "fia original");
        global_dpd_->file2_close(&dpd_fIJ);
        global_dpd_->file2_close(&dpd_fij);
        global_dpd_->file2_close(&dpd_fAB);
        global_dpd_->file2_close(&dpd_fab);
        global_dpd_->file2_close(&dpd_fIA);
        global_dpd_->file2_close(&dpd_fia);
    }
    // Generate Dimensions for slices
    Dimension zero(moinfo.nirreps);

    std::vector<int> tmp(moinfo.nirreps);
    tmp.assign(moinfo.frdocc, moinfo.frdocc + moinfo.nirreps);
    Dimension frdoccpi(tmp);

    tmp.assign(moinfo.aoccpi, moinfo.aoccpi + moinfo.nirreps);
    Dimension aoccpi(tmp);

    tmp.assign(moinfo.boccpi, moinfo.boccpi + moinfo.nirreps);
    Dimension boccpi(tmp);

    tmp.assign(moinfo.avirtpi, moinfo.avirtpi + moinfo.nirreps);
    Dimension avirtpi(tmp);

    tmp.assign(moinfo.bvirtpi, moinfo.bvirtpi + moinfo.nirreps);
    Dimension bvirtpi(tmp);

    global_dpd_->file2_init(&dpd_fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ original");
    SharedMatrix fIJ = std::make_shared<Matrix>(&dpd_fIJ);
    global_dpd_->file2_close(&dpd_fIJ);
    fIJ->set_name("fIJ");
    fIJ->add(MO_OE_operator_a->get_block({frdoccpi, frdoccpi + aoccpi}, {frdoccpi, frdoccpi + aoccpi}));
    global_dpd_->file2_init(&dpd_fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ");
    fIJ->write_to_dpdfile2(&dpd_fIJ);
    global_dpd_->file2_close(&dpd_fIJ);

    global_dpd_->file2_init(&dpd_fij, PSIF_CC_OEI, 0, 2, 2, "fij original");
    SharedMatrix fij = std::make_shared<Matrix>(&dpd_fij);
    global_dpd_->file2_close(&dpd_fij);
    fij->set_name("fij");
    fij->add(MO_OE_operator_b->get_block({frdoccpi, frdoccpi + boccpi}, {frdoccpi, frdoccpi + boccpi}));
    global_dpd_->file2_init(&dpd_fij, PSIF_CC_OEI, 0, 2, 2, "fij");
    fij->write_to_dpdfile2(&dpd_fij);
    global_dpd_->file2_close(&dpd_fij);

    global_dpd_->file2_init(&dpd_fAB, PSIF_CC_OEI, 0, 1, 1, "fAB original");
    SharedMatrix fAB = std::make_shared<Matrix>(&dpd_fAB);
    global_dpd_->file2_close(&dpd_fAB);
    fAB->set_name("fAB");
    fAB->add(MO_OE_operator_a->get_block({frdoccpi + aoccpi, frdoccpi + aoccpi + avirtpi},
                                         {frdoccpi + aoccpi, frdoccpi + aoccpi + avirtpi}));
    global_dpd_->file2_init(&dpd_fAB, PSIF_CC_OEI, 0, 1, 1, "fAB");
    fAB->write_to_dpdfile2(&dpd_fAB);
    global_dpd_->file2_close(&dpd_fAB);

    global_dpd_->file2_init(&dpd_fab, PSIF_CC_OEI, 0, 3, 3, "fab original");
    SharedMatrix fab = std::make_shared<Matrix>(&dpd_fab);
    global_dpd_->file2_close(&dpd_fab);
    fab->set_name("fab");
    fab->add(MO_OE_operator_b->get_block({frdoccpi + boccpi, frdoccpi + boccpi + bvirtpi},
                                         {frdoccpi + boccpi, frdoccpi + boccpi + bvirtpi}));
    global_dpd_->file2_init(&dpd_fab, PSIF_CC_OEI, 0, 3, 3, "fab");
    fab->write_to_dpdfile2(&dpd_fab);
    global_dpd_->file2_close(&dpd_fab);

    global_dpd_->file2_init(&dpd_fIA, PSIF_CC_OEI, 0, 0, 1, "fIA original");
    SharedMatrix fIA = std::make_shared<Matrix>(&dpd_fIA);
    global_dpd_->file2_close(&dpd_fIA);
    fIA->set_name("fIA");
    fIA->add(
        MO_OE_operator_a->get_block({frdoccpi, frdoccpi + aoccpi}, {frdoccpi + aoccpi, frdoccpi + aoccpi + avirtpi}));
    global_dpd_->file2_init(&dpd_fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
    fIA->write_to_dpdfile2(&dpd_fIA);
    global_dpd_->file2_close(&dpd_fIA);

    global_dpd_->file2_init(&dpd_fia, PSIF_CC_OEI, 0, 2, 3, "fia original");
    SharedMatrix fia = std::make_shared<Matrix>(&dpd_fia);
    global_dpd_->file2_close(&dpd_fia);
    fia->set_name("fia");
    fia->add(
        MO_OE_operator_b->get_block({frdoccpi, frdoccpi + boccpi}, {frdoccpi + boccpi, frdoccpi + boccpi + bvirtpi}));
    global_dpd_->file2_init(&dpd_fia, PSIF_CC_OEI, 0, 2, 3, "fia");
    fia->write_to_dpdfile2(&dpd_fia);
    global_dpd_->file2_close(&dpd_fia);
}
}
}  // namespace psi::ccdensity

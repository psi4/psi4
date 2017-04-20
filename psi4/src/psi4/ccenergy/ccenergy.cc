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
    \ingroup CCENERGY
    \brief Enter brief description of file here
*/
/*
**  CCENERGY: Program to calculate coupled cluster energies.
*/

#include "Params.h"
#include "MOInfo.h"
#include "Local.h"
#include "ccwave.h"

#include "psi4/libciomr/libciomr.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/psifiles.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/liboptions/liboptions.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <sys/types.h>
#include <unistd.h>
#include <sys/types.h>

#ifdef USING_PCMSolver
#include "psi4/libpsipcm/psipcm.h"
#endif

namespace psi { namespace ccenergy {

#define IOFF_MAX 32641

}} //namespace psi::ccenergy

// Forward declaration to call cctriples
namespace psi { namespace cctriples {
PsiReturnType cctriples(std::shared_ptr<Wavefunction> reference_wavefunction, Options &options);
}}

namespace psi { namespace ccenergy {

CCEnergyWavefunction::CCEnergyWavefunction(std::shared_ptr<Wavefunction> reference_wavefunction, Options &options)
    : Wavefunction(options)
{
    set_reference_wavefunction(reference_wavefunction);
    init();
    #define NUM_ENTRIES 113
    cache_priority_list_ = new dpd_file4_cache_entry[NUM_ENTRIES];
}

CCEnergyWavefunction::~CCEnergyWavefunction()
{
    if(cache_priority_list_)
        delete [] cache_priority_list_;
}

void CCEnergyWavefunction::init()
{
    shallow_copy(reference_wavefunction_);
}

double CCEnergyWavefunction::compute_energy()
{
    int done=0, brueckner_done=0;
    int h, i, j, a, b, row, col, natom;
    double **geom, *zvals, value;
    FILE *efile;
    int **cachelist, *cachefiles;
    dpdfile2 t1;
    dpdbuf4 t2;
    double *emp2_aa, *emp2_ab, *ecc_aa, *ecc_ab, tval;

    moinfo_.iter=0;

    init_io();
    init_ioff();
    title();

#ifdef TIME_CCENERGY
    timer_on("CCEnergy");
#endif

    get_moinfo();
    get_params(options_);

    cachefiles = init_int_array(PSIO_MAXUNIT);

    if(params_.ref == 2) { /** UHF **/
        cachelist = cacheprep_uhf(params_.cachelev, cachefiles);

        std::vector<int*> spaces;
        spaces.push_back(moinfo_.aoccpi);
        spaces.push_back(moinfo_.aocc_sym);
        spaces.push_back(moinfo_.avirtpi);
        spaces.push_back(moinfo_.avir_sym);
        spaces.push_back(moinfo_.boccpi);
        spaces.push_back(moinfo_.bocc_sym);
        spaces.push_back(moinfo_.bvirtpi);
        spaces.push_back(moinfo_.bvir_sym);
        delete[] dpd_list[0];
        dpd_list[0] = new DPD(0, moinfo_.nirreps, params_.memory, 0, cachefiles,
                              cachelist, NULL, 4, spaces);
        dpd_set_default(0);

        if( params_.df ){
            form_df_ints(options_, cachelist, cachefiles, cache_priority_list_);
        }else if( params_.aobasis != "NONE" ) { /* Set up new DPD's for AO-basis algorithm */
            std::vector<int*> aospaces;
            aospaces.push_back(moinfo_.aoccpi);
            aospaces.push_back(moinfo_.aocc_sym);
            aospaces.push_back(moinfo_.sopi);
            aospaces.push_back(moinfo_.sosym);
            aospaces.push_back(moinfo_.boccpi);
            aospaces.push_back(moinfo_.bocc_sym);
            aospaces.push_back(moinfo_.sopi);
            aospaces.push_back(moinfo_.sosym);
            dpd_init(1, moinfo_.nirreps, params_.memory, 0, cachefiles, cachelist, NULL, 4, aospaces);
            dpd_set_default(0);
        }

    }
    else { /** RHF or ROHF **/
        cachelist = cacheprep_rhf(params_.cachelev, cachefiles);

        init_priority_list();
        std::vector<int*> spaces;
        spaces.push_back(moinfo_.occpi);
        spaces.push_back(moinfo_.occ_sym);
        spaces.push_back(moinfo_.virtpi);
        spaces.push_back(moinfo_.vir_sym);

        dpd_init(0, moinfo_.nirreps, params_.memory, params_.cachetype, cachefiles, cachelist, cache_priority_list_, 2, spaces);

        if( params_.df ){
            form_df_ints(options_, cachelist, cachefiles, cache_priority_list_);
        }else if( params_.aobasis != "NONE") { /* Set up new DPD for AO-basis algorithm */
            std::vector<int*> aospaces;
            aospaces.push_back(moinfo_.occpi);
            aospaces.push_back(moinfo_.occ_sym);
            aospaces.push_back(moinfo_.sopi);
            aospaces.push_back(moinfo_.sosym);
            dpd_init(1, moinfo_.nirreps, params_.memory, 0, cachefiles, cachelist, NULL, 2, aospaces);
            dpd_set_default(0);
        }

    }

    if ( (params_.just_energy) || (params_.just_residuals) ) {
        one_step();
        if(params_.ref == 2) cachedone_uhf(cachelist); else cachedone_rhf(cachelist);
        free(cachefiles);
        cleanup();
        free(ioff_);
        exit_io();
        return Success;
    }

    if(params_.local) {
        local_init();
        if(local_.weakp=="MP2") lmp2();
    }

    init_amps();

    /* Compute the MP2 energy while we're here */
    if(params_.ref == 0 || params_.ref == 2) {
        moinfo_.emp2 = mp2_energy();
        outfile->Printf("MP2 correlation energy %4.16f\n", moinfo_.emp2);
        psio_write_entry(PSIF_CC_INFO, "MP2 Energy", (char *) &(moinfo_.emp2),sizeof(double));
        Process::environment.globals["MP2 CORRELATION ENERGY"] = moinfo_.emp2;
        Process::environment.globals["MP2 TOTAL ENERGY"] = moinfo_.emp2 + moinfo_.eref;
    }

    if(params_.print_mp2_amps) amp_write();

    tau_build();
    taut_build();
    outfile->Printf( "                Solving CC Amplitude Equations\n");
    outfile->Printf( "                ------------------------------\n");
    outfile->Printf( "  Iter             Energy              RMS        T1Diag      D1Diag    New D1Diag    D2Diag\n");
    outfile->Printf( "  ----     ---------------------    ---------   ----------  ----------  ----------   --------\n");
    moinfo_.ecc = energy();
    pair_energies(&emp2_aa, &emp2_ab);
    double last_energy;

    moinfo_.t1diag = diagnostic();
    moinfo_.d1diag = d1diag();
    moinfo_.new_d1diag = new_d1diag();

    moinfo_.d2diag = d2diag();
    update();
    checkpoint();
    for(moinfo_.iter=1; moinfo_.iter <= params_.maxiter; moinfo_.iter++) {

        sort_amps();

        if (options_.get_bool("PCM") && (options_.get_str("PCM_CC_TYPE") == "PTED")) {
          Process::environment.globals["T MICROITERATIONS"] = (moinfo_.iter + 1);
        }

#ifdef TIME_CCENERGY
        timer_on("F build");
#endif
        /* If the user asked for a PCM-CC-PTES calculation, we have to update
         * the Fock matrix by adding the PCM potential at this point in the
         * setup of the T equations.
         * The Fock matrix update is done as in ccdensity.cc, by first saving the
         * original Fock matrix.
         */
        if (options_.get_bool("PCM") && options_.get_str("PCM_CC_TYPE") == "PTES") {
          // PTES scheme: M. Caricato, JCP, 135, 074113 (2011)
          //   Build the PCM polarization charges from the t_i^a CC singles amplitudes at the current iteration.
          //   The polarization charges are then contracted with the electrostatic potential integrals
          //   to build the PCM potential to be added in the CC T equations.
          SharedPCM cc_pcm = std::make_shared<PCM>(options_,
              reference_wavefunction_->psio(), moinfo_.nirreps,
              reference_wavefunction_->basisset());
          if (params_.ref == 0) { /** RHF **/
            SharedMatrix MO_t1_A(new Matrix("MO_t1_A", moinfo_.nirreps, moinfo_.orbspi, moinfo_.orbspi));
            SharedMatrix SO_t1_A(new Matrix("SO_t1_A", moinfo_.nirreps, reference_wavefunction_->nsopi(), reference_wavefunction_->nsopi()));

            SharedMatrix MO_PCM_potential(new Matrix("MO_PCM_potential", moinfo_.nirreps, moinfo_.orbspi, moinfo_.orbspi));

            get_t1_rhf(MO_t1_A);

            SharedMatrix C = reference_wavefunction_->Ca();
            SO_t1_A->back_transform(MO_t1_A, C);

            // We have to get the 0.5 * (QV) energy contribution
            double Epcm_correlated = cc_pcm->compute_E(SO_t1_A, PCM::EleOnly);
            Process::environment.globals["PCM-CC-PTES CORRELATED POLARIZATION ENERGY"] = Epcm_correlated;
            double E_correlation = moinfo_.ecc; // Just the CCSD part
            E_correlation += Epcm_correlated; // We add the PCM contribution on top
            Process::environment.globals["CURRENT CORRELATION ENERGY"] = E_correlation; // Save into the globals array
            double E = Process::environment.globals["CURRENT ENERGY"];
            E += Epcm_correlated;
            Process::environment.globals["CURRENT ENERGY"] = E;
            //outfile->Printf("Epol_correlated = %20.12f\n", Epol_correlated);
            SharedMatrix SO_PCM_potential = cc_pcm->compute_V_electronic(); // This is in SO basis
            // We now transform it to MO basis...
            MO_PCM_potential->transform(SO_PCM_potential, C);
            //MO_PCM_potential->print();
            update_F_pcm_rhf(MO_PCM_potential);
          } else if (params_.ref == 1) {/** ROHF case **/
            SharedMatrix MO_t1_A(new Matrix("MO_t1_A", moinfo_.nirreps, moinfo_.orbspi, moinfo_.orbspi));
            SharedMatrix MO_t1_B(new Matrix("MO_t1_B", moinfo_.nirreps, moinfo_.orbspi, moinfo_.orbspi));
            SharedMatrix SO_t1_A(new Matrix("SO_t1_A", moinfo_.nirreps, reference_wavefunction_->nsopi(), reference_wavefunction_->nsopi()));
            SharedMatrix SO_t1_B(new Matrix("SO_t1_B", moinfo_.nirreps, reference_wavefunction_->nsopi(), reference_wavefunction_->nsopi()));

            SharedMatrix MO_PCM_potential(new Matrix("MO_PCM_potential", moinfo_.nirreps, moinfo_.orbspi, moinfo_.orbspi));

            get_t1_rohf(MO_t1_A, MO_t1_B);

            SharedMatrix C = reference_wavefunction_->Ca();
            SO_t1_A->back_transform(MO_t1_A, C);
            SO_t1_B->back_transform(MO_t1_B, C);

            SO_t1_A->add(SO_t1_B);

            // We have to get the 0.5 * (QV) energy contribution
            double Epcm_correlated = cc_pcm->compute_E(SO_t1_A, PCM::EleOnly);
            Process::environment.globals["PCM-CC-PTES CORRELATED POLARIZATION ENERGY"] = Epcm_correlated;
            double E_correlation = moinfo_.ecc; // Just the CCSD part
            E_correlation += Epcm_correlated; // We add the PCM contribution on top
            Process::environment.globals["CURRENT CORRELATION ENERGY"] = E_correlation; // Save into the globals array
            double E = Process::environment.globals["CURRENT ENERGY"];
            E += Epcm_correlated;
            Process::environment.globals["CURRENT ENERGY"] = E;
            //outfile->Printf("Epol_correlated = %20.12f\n", Epol_correlated);
            SharedMatrix SO_PCM_potential = cc_pcm->compute_V_electronic(); // This is in SO basis
            // We now transform it to MO basis...
            MO_PCM_potential->transform(SO_PCM_potential, C);
            //MO_PCM_potential->print();
            update_F_pcm_rhf(MO_PCM_potential);
          } else if (params_.ref == 2) {/** UHF case **/
            SharedMatrix MO_t1_A(new Matrix("MO_t1_A", moinfo_.nirreps, moinfo_.orbspi, moinfo_.orbspi));
            SharedMatrix MO_t1_B(new Matrix("MO_t1_B", moinfo_.nirreps, moinfo_.orbspi, moinfo_.orbspi));
            SharedMatrix SO_t1_A(new Matrix("SO_t1_A", moinfo_.nirreps, reference_wavefunction_->nsopi(), reference_wavefunction_->nsopi()));
            SharedMatrix SO_t1_B(new Matrix("SO_t1_B", moinfo_.nirreps, reference_wavefunction_->nsopi(), reference_wavefunction_->nsopi()));

            get_t1_uhf(MO_t1_A, MO_t1_B);

            SharedMatrix MO_PCM_potential_A(new Matrix("MO_PCM_potential_A", moinfo_.nirreps, moinfo_.orbspi, moinfo_.orbspi));
            SharedMatrix MO_PCM_potential_B(new Matrix("MO_PCM_potential_B", moinfo_.nirreps, moinfo_.orbspi, moinfo_.orbspi));

            SharedMatrix Ca = reference_wavefunction_->Ca();
            SharedMatrix Cb = reference_wavefunction_->Cb();
            SO_t1_A->back_transform(MO_t1_A, Ca);
            SO_t1_B->back_transform(MO_t1_B, Cb);
            SO_t1_A->add(SO_t1_B);

            // We have to get the 0.5 * (QV) energy contribution
            double Epcm_correlated = cc_pcm->compute_E(SO_t1_A, PCM::EleOnly);
            Process::environment.globals["PCM-CC-PTES CORRELATED POLARIZATION ENERGY"] = Epcm_correlated;
            double E_correlation = moinfo_.ecc; // Just the CCSD part
            E_correlation += Epcm_correlated; // We add the PCM contribution on top
            Process::environment.globals["CURRENT CORRELATION ENERGY"] = E_correlation; // Save into the globals array
            double E = Process::environment.globals["CURRENT ENERGY"];
            E += Epcm_correlated;
            Process::environment.globals["CURRENT ENERGY"] = E;
            SharedMatrix SO_PCM_potential = cc_pcm->compute_V_electronic(); // This is in SO basis
            // We now transform it to MO basis...
            MO_PCM_potential_A->transform(SO_PCM_potential, Ca);
            MO_PCM_potential_B->transform(SO_PCM_potential, Cb);
            //MO_PCM_potential->print();
            update_F_pcm_uhf(MO_PCM_potential_A, MO_PCM_potential_B);
          }
        }

        Fme_build(); Fae_build(); Fmi_build();
        if(params_.print & 2) status("F intermediates", "outfile");
#ifdef TIME_CCENERGY
        timer_off("F build");
#endif

        t1_build();
        if(params_.print & 2) status("T1 amplitudes", "outfile");

        if( params_.wfn == "CC2"  || params_.wfn == "EOM_CC2" ) {

            cc2_Wmnij_build();
            if(params_.print & 2) status("Wmnij", "outfile");

#ifdef TIME_CCENERGY
            timer_on("Wmbij build");
#endif
            cc2_Wmbij_build();
            if(params_.print & 2) status("Wmbij", "outfile");
#ifdef TIME_CCENERGY
            timer_off("Wmbij build");
#endif

#ifdef TIME_CCENERGY
            timer_on("Wabei build");
#endif
            cc2_Wabei_build();
            if(params_.print & 2) status("Wabei", "outfile");
#ifdef TIME_CCENERGY
            timer_off("Wabei build");
#endif

#ifdef TIME_CCENERGY
            timer_on("T2 Build");
#endif
            cc2_t2_build();
            if(params_.print & 2) status("T2 amplitudes", "outfile");
#ifdef TIME_CCENERGY
            timer_off("T2 Build");
#endif

        }

        else {

#ifdef TIME_CCENERGY
            timer_on("Wmbej build");
#endif
            Wmbej_build();
            if(params_.print & 2) status("Wmbej", "outfile");
#ifdef TIME_CCENERGY
            timer_off("Wmbej build");
#endif

            Z_build();
            if(params_.print & 2) status("Z", "outfile");
            Wmnij_build();
            if(params_.print & 2) status("Wmnij", "outfile");

#ifdef TIME_CCENERGY
            timer_on("T2 Build");
#endif
            t2_build();
            if(params_.print & 2) status("T2 amplitudes", "outfile");
#ifdef TIME_CCENERGY
            timer_off("T2 Build");
#endif

            if( params_.wfn == "CC3" || params_.wfn == "EOM_CC3" ) {

                /* step1: build cc3 intermediates, Wabei, Wmnie, Wmbij, Wamef */
                cc3_Wmnij();
                cc3_Wmbij();
                cc3_Wmnie();
                cc3_Wamef();
                cc3_Wabei();

                /* step2: loop over T3's and add contributions to T1 and T2 as you go */
                cc3();
            }
        }

        if (!params_.just_residuals)
            denom(); /* apply denominators to T1 and T2 */

        if(converged(last_energy - moinfo_.ecc)) {
          done = 1;

          tsave();
          tau_build(); taut_build();
          last_energy = moinfo_.ecc;
          moinfo_.ecc = energy();
          moinfo_.t1diag = diagnostic();
          moinfo_.d1diag = d1diag();
          moinfo_.new_d1diag = new_d1diag();
          moinfo_.d2diag = d2diag();
          sort_amps();
          update();
          outfile->Printf( "\n    Iterations converged.\n");

          outfile->Printf( "\n");
          amp_write();
          /* What to save in the globals array.
           * PCM-CC-PTED: the PCM-CC-PTE and PCM-CC-PTE(S) correlation energies,
           * the PCM-CC-PTE(S) correlated polarization energy. Only at the first MACROiteration.
           * PCM-CC-PTES: the PCM-CC-PTES correlation energy and correlated polarization energy.
           * Within this scheme we don't have access to the PTE and PTE(S) quantities.
           * PCM-CC-PTE(S): the PCM-CC-PTE and PCM-CC-PTE(S) correlation energies,
           * the PCM-CC-PTE(S) correlated polarization energy.
           * PCM-CC-PTE: the PCM-CC-PTE correlation energy.
           */
          if (options_.get_bool("PCM") && (options_.get_str("PCM_CC_TYPE") == "PTE(S)" || options_.get_str("PCM_CC_TYPE") == "PTED")) {
            bool macroiter_0 = (int(Process::environment.globals["MACROITERATION"]) == 0);
            bool pte_s_calc = (options_.get_str("PCM_CC_TYPE") == "PTE(S)");
            bool pted_calc = (options_.get_str("PCM_CC_TYPE") == "PTED");
            SharedMatrix SO_converged_t1_A(new Matrix("SO_converged_t1_A", reference_wavefunction_->nsopi(), reference_wavefunction_->nsopi()));
            if (pte_s_calc || (pted_calc && macroiter_0)) {
              SharedMatrix MO_converged_t1_A(new Matrix("MO_converged_t1_A", moinfo_.nirreps, moinfo_.orbspi, moinfo_.orbspi));
              if (params_.ref == 0) { /** RHF **/
                get_t1_rhf(MO_converged_t1_A);

                SharedMatrix C = reference_wavefunction_->Ca();
                SO_converged_t1_A->back_transform(MO_converged_t1_A, C);
              } else if (params_.ref == 1) { /** ROHF **/
                SharedMatrix MO_converged_t1_B(new Matrix("MO_converged_t1_B", moinfo_.nirreps, moinfo_.orbspi, moinfo_.orbspi));
                get_t1_rohf(MO_converged_t1_A, MO_converged_t1_B);

                SharedMatrix C = reference_wavefunction_->Ca();
                SharedMatrix SO_converged_t1_A(new Matrix("SO_converged_t1_A", reference_wavefunction_->nsopi(), reference_wavefunction_->nsopi()));
                SharedMatrix SO_converged_t1_B(new Matrix("SO_converged_t1_B", reference_wavefunction_->nsopi(), reference_wavefunction_->nsopi()));
                SO_converged_t1_A->back_transform(MO_converged_t1_A, C);
                SO_converged_t1_B->back_transform(MO_converged_t1_B, C);

                SO_converged_t1_A->add(SO_converged_t1_B);
              } else if (params_.ref == 2) { /** UHF **/
                SharedMatrix MO_converged_t1_B(new Matrix("MO_converged_t1_B", moinfo_.nirreps, moinfo_.orbspi, moinfo_.orbspi));
                get_t1_uhf(MO_converged_t1_A, MO_converged_t1_B);

                SharedMatrix Ca = reference_wavefunction_->Ca();
                SharedMatrix Cb = reference_wavefunction_->Cb();
                SharedMatrix SO_converged_t1_A(new Matrix("SO_converged_t1_A", reference_wavefunction_->nsopi(), reference_wavefunction_->nsopi()));
                SharedMatrix SO_converged_t1_B(new Matrix("SO_converged_t1_B", reference_wavefunction_->nsopi(), reference_wavefunction_->nsopi()));
                SO_converged_t1_A->back_transform(MO_converged_t1_A, Ca);
                SO_converged_t1_B->back_transform(MO_converged_t1_B, Cb);

                SO_converged_t1_A->add(SO_converged_t1_B);
              }
            }
            // We have to get the 0.5 * (QV) energy contribution
            std::shared_ptr<PCM> cc_pcm = std::make_shared<PCM>(options_,
                reference_wavefunction_->psio(), moinfo_.nirreps,
                reference_wavefunction_->basisset());
            double Epcm_correlated = cc_pcm->compute_E(SO_converged_t1_A, PCM::EleOnly);
            if (pte_s_calc) {
              Process::environment.globals["PCM-CC-PTE(S) CORRELATED POLARIZATION ENERGY"] = Epcm_correlated;
              double E_correlation = moinfo_.ecc; //Process::environment.globals["CURRENT CORRELATION ENERGY"];
              E_correlation += Epcm_correlated;
              Process::environment.globals["CURRENT CORRELATION ENERGY"] = E_correlation;
              double E = moinfo_.eref + moinfo_.ecc; //Process::environment.globals["CURRENT ENERGY"];
              E += Epcm_correlated;
              Process::environment.globals["CURRENT ENERGY"] = E;
              outfile->Printf("PCM-CC-PTE(S) correlated polarization energy %16.14f\n", Epcm_correlated);
            } else if (pted_calc && macroiter_0) {
              Process::environment.globals["PCM-CC-PTE(S) CORRELATED POLARIZATION ENERGY"] = Epcm_correlated;
              double E_correlation = moinfo_.ecc; //Process::environment.globals["CURRENT CORRELATION ENERGY"];
              Process::environment.globals["PCM-CC-PTE CORRELATION ENERGY"] = E_correlation;
            }
          }
          if (params_.analyze != 0) analyze();
          break;
        }
        if(params_.diis) diis(moinfo_.iter);
        tsave();
        tau_build(); taut_build();
        last_energy = moinfo_.ecc;
        moinfo_.ecc = energy();
        moinfo_.t1diag = diagnostic();
        moinfo_.d1diag = d1diag();
        moinfo_.new_d1diag = new_d1diag();
        moinfo_.d2diag = d2diag();
        update();
        checkpoint();
    }  // end loop over iterations

    // DGAS Edit
    Process::environment.globals["CC T1 DIAGNOSTIC"] = moinfo_.t1diag;
    Process::environment.globals["CC D1 DIAGNOSTIC"] = moinfo_.d1diag;
    Process::environment.globals["CC NEW D1 DIAGNOSTIC"] = moinfo_.new_d1diag;
    Process::environment.globals["CC D2 DIAGNOSTIC"] = moinfo_.d2diag;

    outfile->Printf( "\n");
    if(!done) {
        outfile->Printf( "     ** Wave function not converged to %2.1e ** \n",
                params_.convergence);

        if( params_.aobasis != "NONE" ) dpd_close(1);
        dpd_close(0);
        cleanup();
#ifdef TIME_CCENERGY
        timer_off("CCEnergy");
#endif
        free(ioff_);
        exit_io();
        return Failure;
    }

    outfile->Printf( "    SCF energy       (wfn)                    = %20.15f\n", moinfo_.escf);
    outfile->Printf( "    Reference energy (file100)                = %20.15f\n", moinfo_.eref);

    //Process::environment.globals["SCF TOTAL ENERGY (CHKPT)"] = moinfo_.escf;
    //Process::environment.globals["SCF TOTAL ENERGY"] = moinfo_.eref;

    if(params_.ref == 0 || params_.ref == 2) {
        if (params_.scs) {
            outfile->Printf( "\n    OS SCS-MP2 correlation energy             = %20.15f\n", moinfo_.emp2_os*params_.scsmp2_scale_os);
            outfile->Printf( "    SS SCS-MP2 correlation energy             = %20.15f\n", moinfo_.emp2_ss*params_.scsmp2_scale_ss);
            outfile->Printf( "    SCS-MP2 correlation energy                = %20.15f\n", moinfo_.emp2_os*params_.scsmp2_scale_os
                    + moinfo_.emp2_ss*params_.scsmp2_scale_ss);
            outfile->Printf( "      * SCS-MP2 total energy                  = %20.15f\n", moinfo_.eref
                    + moinfo_.emp2_os*params_.scsmp2_scale_os + moinfo_.emp2_ss*params_.scsmp2_scale_ss);

            Process::environment.globals["SCS-MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = moinfo_.emp2_os*params_.scsmp2_scale_os;
            Process::environment.globals["SCS-MP2 SAME-SPIN CORRELATION ENERGY"] = moinfo_.emp2_ss*params_.scsmp2_scale_ss;
            Process::environment.globals["SCS-MP2 CORRELATION ENERGY"] = moinfo_.emp2_os*params_.scsmp2_scale_os +
                    moinfo_.emp2_ss*params_.scsmp2_scale_ss;
            Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = moinfo_.eref +
                    moinfo_.emp2_os*params_.scsmp2_scale_os + moinfo_.emp2_ss*params_.scsmp2_scale_ss;
        }
        if (params_.scsn) {
            outfile->Printf( "\nOS SCSN-MP2 correlation energy            = %20.15f\n", 0.0);
            outfile->Printf( "    SS SCSN-MP2 correlation energy        = %20.15f\n", moinfo_.emp2_ss*1.76);
            outfile->Printf( "    SCSN-MP2 correlation energy           = %20.15f\n", moinfo_.emp2_ss*1.76);
            outfile->Printf( "      * SCSN-MP2 total energy             = %20.15f\n", moinfo_.eref + moinfo_.emp2_ss*1.76);

            Process::environment.globals["SCSN-MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = 0.0;
            Process::environment.globals["SCSN-MP2 SAME-SPIN CORRELATION ENERGY"] = moinfo_.emp2_ss*1.76;
            Process::environment.globals["SCSN-MP2 CORRELATION ENERGY"] = moinfo_.emp2_ss*1.76;
            Process::environment.globals["SCSN-MP2 TOTAL ENERGY"] = moinfo_.eref + moinfo_.emp2_ss*1.76;
        }

        outfile->Printf( "\n    Opposite-spin MP2 correlation energy      = %20.15f\n", moinfo_.emp2_os);
        outfile->Printf( "    Same-spin MP2 correlation energy          = %20.15f\n", moinfo_.emp2_ss);
        outfile->Printf( "    MP2 correlation energy                    = %20.15f\n", moinfo_.emp2);
        outfile->Printf( "      * MP2 total energy                      = %20.15f\n", moinfo_.eref + moinfo_.emp2);

        Process::environment.globals["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = moinfo_.emp2_os;
        Process::environment.globals["MP2 SAME-SPIN CORRELATION ENERGY"] = moinfo_.emp2_ss;
        Process::environment.globals["MP2 CORRELATION ENERGY"] = moinfo_.emp2;
        Process::environment.globals["MP2 TOTAL ENERGY"] = moinfo_.eref + moinfo_.emp2;
    }
    if( params_.wfn == "CC3"  || params_.wfn == "EOM_CC3" ) {
        outfile->Printf( "    CC3 correlation energy                    = %20.15f\n", moinfo_.ecc);
        outfile->Printf( "      * CC3 total energy                      = %20.15f\n", moinfo_.eref + moinfo_.ecc);

        Process::environment.globals["CC3 CORRELATION ENERGY"] = moinfo_.ecc;
        Process::environment.globals["CC3 TOTAL ENERGY"] = moinfo_.eref + moinfo_.ecc;
    }
    else if( params_.wfn == "CC2" || params_.wfn == "EOM_CC2" )  {
        outfile->Printf( "    CC2 correlation energy                    = %20.15f\n", moinfo_.ecc);
        outfile->Printf( "      * CC2 total energy                      = %20.15f\n", moinfo_.eref + moinfo_.ecc);

        Process::environment.globals["CC2 CORRELATION ENERGY"] = moinfo_.ecc;
        Process::environment.globals["CC2 TOTAL ENERGY"] = moinfo_.eref + moinfo_.ecc;

        if(params_.local && local_.weakp == "MP2" )
            outfile->Printf( "      * LCC2 (+LMP2) total energy             = %20.15f\n",
                    moinfo_.eref + moinfo_.ecc + local_.weak_pair_energy);
        Process::environment.globals["LCC2 (+LMP2) TOTAL ENERGY"] =
                moinfo_.eref + moinfo_.ecc + local_.weak_pair_energy;
    }
    else {
        if (params_.scscc) {
            outfile->Printf( "\n    OS SCS-CCSD correlation energy            = %20.15f\n", moinfo_.ecc_os*params_.scscc_scale_os);
            outfile->Printf( "    SS SCS-CCSD correlation energy            = %20.15f\n", moinfo_.ecc_ss*params_.scscc_scale_ss);
            outfile->Printf( "    SCS-CCSD correlation energy               = %20.15f\n", moinfo_.ecc_os*params_.scscc_scale_os
                    + moinfo_.ecc_ss*params_.scscc_scale_ss);
            outfile->Printf( "      * SCS-CCSD total energy                 = %20.15f\n", moinfo_.eref
                    + moinfo_.ecc_os*params_.scscc_scale_os + moinfo_.ecc_ss*params_.scscc_scale_ss);

            // LAB TODO  reconsider variable names for ss/os cc
            Process::environment.globals["SCS-CCSD OPPOSITE-SPIN CORRELATION ENERGY"] = moinfo_.ecc_os*params_.scscc_scale_os;
            Process::environment.globals["SCS-CCSD SAME-SPIN CORRELATION ENERGY"] = moinfo_.ecc_ss*params_.scscc_scale_ss;
            Process::environment.globals["SCS-CCSD CORRELATION ENERGY"] = moinfo_.ecc_os*params_.scscc_scale_os +
                    moinfo_.ecc_ss*params_.scscc_scale_ss;
            Process::environment.globals["SCS-CCSD TOTAL ENERGY"] = moinfo_.eref +
                    moinfo_.ecc_os*params_.scscc_scale_os + moinfo_.ecc_ss*params_.scscc_scale_ss;
        }

        outfile->Printf( "\n    Opposite-spin CCSD correlation energy     = %20.15f\n", moinfo_.ecc_os);
        outfile->Printf( "    Same-spin CCSD correlation energy         = %20.15f\n", moinfo_.ecc_ss);

        if (options_.get_bool("PCM") && (options_.get_str("PCM_CC_TYPE") == "PTES")) {
          double ptes_energy = Process::environment.globals["PCM-CC-PTES CORRELATED POLARIZATION ENERGY"];
          outfile->Printf("\tPTES correlated polarization energy   = %20.15f\n", ptes_energy);
          outfile->Printf("\tCCSD correlation energy               = %20.15f\n", moinfo_.ecc);
          outfile->Printf("\tPCM-PTES-CCSD correlation energy      = %20.15f\n", moinfo_.ecc + ptes_energy);
          outfile->Printf("      * PCM-PTES-CCSD total energy            = %20.15f\n", moinfo_.eref + moinfo_.ecc + ptes_energy);
        } else if (options_.get_bool("PCM") && (options_.get_str("PCM_CC_TYPE") == "PTE(S)")) {
          double pte_s_energy = Process::environment.globals["PCM-CC-PTE(S) CORRELATED POLARIZATION ENERGY"];
          outfile->Printf("\tPTE(S) correlated polarization energy   = %20.15f\n", pte_s_energy);
          outfile->Printf("\tCCSD correlation energy                 = %20.15f\n", moinfo_.ecc);
          outfile->Printf("\tPCM-PTE(S)-CCSD correlation energy      = %20.15f\n", moinfo_.ecc + pte_s_energy);
          outfile->Printf("      * PCM-PTE(S)-CCSD total energy      = %20.15f\n", moinfo_.eref + moinfo_.ecc + pte_s_energy);
        } else {
          outfile->Printf( "    CCSD correlation energy                   = %20.15f\n", moinfo_.ecc);
          outfile->Printf( "      * CCSD total energy                     = %20.15f\n", moinfo_.eref + moinfo_.ecc);
        }

        Process::environment.globals["CCSD OPPOSITE-SPIN CORRELATION ENERGY"] = moinfo_.ecc_os;
        Process::environment.globals["CCSD SAME-SPIN CORRELATION ENERGY"] = moinfo_.ecc_ss;

        if (options_.get_bool("PCM") && (options_.get_str("PCM_CC_TYPE") == "PTES")) {
          double ptes_energy = Process::environment.globals["PCM-CC-PTES CORRELATED POLARIZATION ENERGY"];
          Process::environment.globals["CCSD CORRELATION ENERGY"] = moinfo_.ecc + ptes_energy;
          Process::environment.globals["CCSD TOTAL ENERGY"] = moinfo_.ecc + moinfo_.eref + ptes_energy;
        } else if (options_.get_bool("PCM") && (options_.get_str("PCM_CC_TYPE") == "PTE(S)")) {
          double pte_s_energy = Process::environment.globals["PCM-CC-PTE(S) CORRELATED POLARIZATION ENERGY"];
          Process::environment.globals["CCSD CORRELATION ENERGY"] = moinfo_.ecc + pte_s_energy;
          Process::environment.globals["CCSD TOTAL ENERGY"] = moinfo_.ecc + moinfo_.eref + pte_s_energy;
        } else {
          Process::environment.globals["CCSD CORRELATION ENERGY"] = moinfo_.ecc;
          Process::environment.globals["CCSD TOTAL ENERGY"] = moinfo_.ecc + moinfo_.eref;
        }

        if(params_.local && local_.weakp == "MP2" )
            outfile->Printf( "      * LCCSD (+LMP2) total energy            = %20.15f\n",
                    moinfo_.eref + moinfo_.ecc + local_.weak_pair_energy);
        Process::environment.globals["LCCSD (+LMP2) TOTAL ENERGY"] =
                moinfo_.eref + moinfo_.ecc + local_.weak_pair_energy;
    }
    outfile->Printf( "\n");

    /* Generate the spin-adapted RHF amplitudes for later codes */
    if(params_.ref == 0) spinad_amps();

    /* Compute pair energies */
    if(params_.print_pair_energies) {
        pair_energies(&ecc_aa, &ecc_ab);
        print_pair_energies(emp2_aa, emp2_ab, ecc_aa, ecc_ab);
    }

    if(params_.wfn == "CC2" && params_.dertype ==1) t1_ijab();

    if( (params_.wfn == "CC3" || params_.wfn == "EOM_CC3" )
            && (params_.dertype == 1 || params_.dertype == 3) && params_.ref == 0) {
        params_.ref = 1;
        /* generate the ROHF versions of the He^T1 intermediates */
        cc3_Wmnij();
        cc3_Wmbij();
        cc3_Wmnie();
        cc3_Wamef();
        cc3_Wabei();
        //    params_.ref == 0;
    }

    if(params_.local) {
        /*    local_print_T1_norm(); */
        local_done();
    }

    if(params_.brueckner)
        Process::environment.globals["BRUECKNER CONVERGED"] = rotate();

    if( params_.aobasis != "NONE" ) dpd_close(1);
    dpd_close(0);

    if(params_.ref == 2) cachedone_uhf(cachelist);
    else cachedone_rhf(cachelist);
    free(cachefiles);

    cleanup();

#ifdef TIME_CCENERGY
    timer_off("CCEnergy");
#endif

    energy_ = moinfo_.ecc + moinfo_.eref;
    name_  = "CCSD";
    Process::environment.globals["CURRENT ENERGY"] = moinfo_.ecc+moinfo_.eref;
    Process::environment.globals["CURRENT CORRELATION ENERGY"] = moinfo_.ecc;
    //Process::environment.globals["CC TOTAL ENERGY"] = moinfo_.ecc+moinfo_.eref;
    //Process::environment.globals["CC CORRELATION ENERGY"] = moinfo_.ecc;

    free(ioff_);
    exit_io();
    //  if(params_.brueckner && brueckner_done)
    //     throw FeatureNotImplemented("CCENERGY", "Brueckner end loop", __FILE__, __LINE__);
    //else

    if ((options_.get_str("WFN") == "CCSD_T")) {
        // Run cctriples
        if (psi::cctriples::cctriples(reference_wavefunction_, options_) == Success)
            energy_ = Process::environment.globals["CURRENT ENERGY"];
        else
            energy_ = 0.0;
    }

    return energy_;
}


void CCEnergyWavefunction::init_io()
{
    params_.just_energy = 0;
    params_.just_residuals = 0;
    tstart();
    for(int i =PSIF_CC_MIN; i <= PSIF_CC_MAX; i++) psio_open(i,1);
}

void CCEnergyWavefunction::title(void)
{
    outfile->Printf( "            **************************\n");
    outfile->Printf( "            *                        *\n");
    outfile->Printf( "            *        CCENERGY        *\n");
    outfile->Printf( "            *                        *\n");
    outfile->Printf( "            **************************\n");
}

void CCEnergyWavefunction::exit_io(void)
{
    int i;
    for(i=PSIF_CC_MIN; i < PSIF_CC_TMP; i++) psio_close(i,1);
    for(i=PSIF_CC_TMP; i <= PSIF_CC_TMP11; i++) psio_close(i,0); /* delete CC_TMP files */
    for(i=PSIF_CC_TMP11+1; i <= PSIF_CC_MAX; i++) psio_close(i,1);
    tstop();

}

void CCEnergyWavefunction::init_ioff(void)
{
    int i;
    ioff_ = init_int_array(IOFF_MAX);
    ioff_[0] = 0;
    for(i=1; i < IOFF_MAX; i++) ioff_[i] = ioff_[i-1] + i;
}


void CCEnergyWavefunction::checkpoint(void)
{
    int i;

    for(i=PSIF_CC_MIN; i <= PSIF_CC_MAX; i++) psio_close(i,1);
    for(i=PSIF_CC_MIN; i <= PSIF_CC_MAX; i++) psio_open(i,1);
}

/* just use T's on disk and don't iterate */
void CCEnergyWavefunction::one_step(void) {
    dpdfile2 t1;
    dpdbuf4 t2;
    double tval;

    moinfo_.ecc = energy();
    outfile->Printf("\n    Values computed from T amplitudes on disk.\n");
    outfile->Printf("Reference expectation value computed: %20.15lf\n", moinfo_.ecc);
    psio_write_entry(PSIF_CC_HBAR, "Reference expectation value", (char *) &(moinfo_.ecc), sizeof(double));

    if (params_.just_residuals) {
        Fme_build(); Fae_build(); Fmi_build();
        t1_build();
        Wmbej_build();
        Z_build();
        Wmnij_build();
        t2_build();
        if ( (params_.ref == 0) || (params_.ref == 1) ) {
            global_dpd_->file2_init(&t1, PSIF_CC_OEI, 0, 0, 1, "New tIA");
            global_dpd_->file2_copy(&t1, PSIF_CC_OEI, "FAI residual");
            global_dpd_->file2_close(&t1);
            global_dpd_->file2_init(&t1, PSIF_CC_OEI, 0, 0, 1, "FAI residual");
            tval = global_dpd_->file2_dot_self(&t1);
            global_dpd_->file2_close(&t1);
            outfile->Printf("    Norm squared of <Phi_I^A|Hbar|0> = %20.15lf\n",tval);
        }
        if (params_.ref == 1) {
            global_dpd_->file2_init(&t1, PSIF_CC_OEI, 0, 0, 1, "New tia");
            global_dpd_->file2_copy(&t1, PSIF_CC_OEI, "Fai residual");
            global_dpd_->file2_close(&t1);
        }
        else if (params_.ref == 2) {
            global_dpd_->file2_init(&t1, PSIF_CC_OEI, 0, 2, 3, "New tia");
            global_dpd_->file2_copy(&t1, PSIF_CC_OEI, "Fai residual");
            global_dpd_->file2_close(&t1);
        }
        if (params_.ref == 0) {
            global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
            global_dpd_->buf4_copy(&t2, PSIF_CC_HBAR, "WAbIj residual");
            global_dpd_->buf4_close(&t2);
            global_dpd_->buf4_init(&t2, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "WAbIj residual");
            tval = global_dpd_->buf4_dot_self(&t2);
            outfile->Printf("    Norm squared of <Phi^Ij_Ab|Hbar|0>: %20.15lf\n",tval);
            global_dpd_->buf4_close(&t2);
        }
        else if (params_.ref == 1) {
            global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tIJAB");
            global_dpd_->buf4_copy(&t2, PSIF_CC_HBAR, "WABIJ residual");
            global_dpd_->buf4_close(&t2);
            global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tijab");
            global_dpd_->buf4_copy(&t2, PSIF_CC_HBAR, "Wabij residual");
            global_dpd_->buf4_close(&t2);
            global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
            global_dpd_->buf4_copy(&t2, PSIF_CC_HBAR, "WAbIj residual");
            global_dpd_->buf4_close(&t2);
        }
        else if(params_.ref ==2) {
            global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tIJAB");
            global_dpd_->buf4_copy(&t2, PSIF_CC_HBAR, "WABIJ residual");
            global_dpd_->buf4_close(&t2);
            global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "New tijab");
            global_dpd_->buf4_copy(&t2, PSIF_CC_HBAR, "Wabij residual");
            global_dpd_->buf4_close(&t2);
            global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");
            global_dpd_->buf4_copy(&t2, PSIF_CC_HBAR, "WAbIj residual");
            global_dpd_->buf4_close(&t2);
        }
    }
    return;
}

void CCEnergyWavefunction::get_t1_rhf(SharedMatrix &_t1)
{
  dpdfile2 tIA;
  global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->file2_mat_init(&tIA);
  global_dpd_->file2_mat_rd(&tIA);

  for(int h = 0; h < moinfo_.nirreps; ++h)
  {
    int upper_bound_ij = moinfo_.occpi[h] - moinfo_.openpi[h];
    int upper_bound_IA = moinfo_.virtpi[h] - moinfo_.openpi[h];

    for(int i = 0; i < moinfo_.occpi[h]; ++i)
    {
      int I = i + moinfo_.frdocc[h];
      for(int a = 0; a < upper_bound_IA; ++a)
      {
        int A = a + moinfo_.frdocc[h] + moinfo_.occpi[h];
        _t1->set(h, I, A, tIA.matrix[h][i][a]);
      }
    }
  }

  global_dpd_->file2_close(&tIA);
}

void CCEnergyWavefunction::get_t1_rohf(SharedMatrix &_t1_A, SharedMatrix &_t1_B)
{
  dpdfile2 tIA, tia;
  global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->file2_mat_init(&tIA);
  global_dpd_->file2_mat_rd(&tIA);

  global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");
  global_dpd_->file2_mat_init(&tia);
  global_dpd_->file2_mat_rd(&tia);

  for(int h = 0; h < moinfo_.nirreps; ++h)
  {
    int upper_bound_ij = moinfo_.occpi[h] - moinfo_.openpi[h];
    int upper_bound_IA = moinfo_.virtpi[h] - moinfo_.openpi[h];

    for(int i = 0; i < moinfo_.occpi[h]; ++i)
    {
      int I = i + moinfo_.frdocc[h];
      for(int a = 0; a < upper_bound_IA; ++a)
      {
        int A = a + moinfo_.frdocc[h] + moinfo_.occpi[h];
        _t1_A->set(h, I, A, tIA.matrix[h][i][a]);
      }
    }
    for(int i = 0; i < upper_bound_ij; ++i)
    {
      int I = i + moinfo_.frdocc[h];
      for(int a = 0; a < upper_bound_IA; ++a)
      {
        int A = a + moinfo_.frdocc[h] + moinfo_.occpi[h];
        _t1_B->set(h, I, A, tia.matrix[h][i][a]);
      }
    }
    for(int i = 0; i < upper_bound_ij; ++i)
    {
      int I = i + moinfo_.frdocc[h];
      for(int a = 0; a < moinfo_.openpi[h]; ++a)
      {
        int A = a + upper_bound_IA;
        int aA = a + moinfo_.frdocc[h] + upper_bound_ij;
        _t1_B->set(h, I, aA, tia.matrix[h][i][A]);
      }
    }
  }

  global_dpd_->file2_close(&tIA);
  global_dpd_->file2_close(&tia);
}

void CCEnergyWavefunction::get_t1_uhf(SharedMatrix &_t1_A, SharedMatrix &_t1_B)
{
  dpdfile2 tIA, tia;
  global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->file2_mat_init(&tIA);
  global_dpd_->file2_mat_rd(&tIA);

  global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->file2_mat_init(&tia);
  global_dpd_->file2_mat_rd(&tia);

  for(int h = 0; h < moinfo_.nirreps; ++h)
  {
    for(int i = 0; i < moinfo_.aoccpi[h]; ++i)
    {
      int I = i + moinfo_.frdocc[h];
      for(int a = 0; a < moinfo_.avirtpi[h]; ++a)
      {
        int A = a + moinfo_.frdocc[h] + moinfo_.aoccpi[h];
        _t1_A->set(h, I, A, tIA.matrix[h][i][a]);
      }
    }
    for(int i = 0; i < moinfo_.boccpi[h]; ++i)
    {
      int I = i + moinfo_.frdocc[h];
      for(int a = 0; a < moinfo_.bvirtpi[h]; ++a)
      {
        int A = a + moinfo_.frdocc[h] + moinfo_.boccpi[h];
        _t1_B->set(h, I, A, tia.matrix[h][i][a]);
      }
    }
  }

  global_dpd_->file2_close(&tIA);
  global_dpd_->file2_close(&tia);
}

void CCEnergyWavefunction::update_F_pcm_rhf(SharedMatrix &MO_PCM_potential)
{
  dpdfile2 fIJ, fij, fAB, fab, fIA, fia;
  if (!(_default_psio_lib_->tocentry_exists(PSIF_CC_OEI, "fIJ original"))) // We check just one of the six blocks for existence
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
  for(int h = 0; h < moinfo_.nirreps; ++h)
  {
    int upper_bound_ij = moinfo_.occpi[h] - moinfo_.openpi[h];
    int upper_bound_IA = moinfo_.virtpi[h] - moinfo_.openpi[h];
    for(int i = 0; i < moinfo_.occpi[h]; ++i)
    {
      int I = i + moinfo_.frdocc[h];
      for(int j = 0; j < moinfo_.occpi[h]; ++j)
      {
        int J = j + moinfo_.frdocc[h];
        fIJ2.matrix[h][i][j] = fIJ.matrix[h][i][j] + MO_PCM_potential->get(h, I, J);
      }
    }
    for(int i = 0; i < upper_bound_ij; ++i)
    {
      int I = i + moinfo_.frdocc[h];
      for(int j = 0; j < upper_bound_ij; ++j)
      {
        int J = j + moinfo_.frdocc[h];
        fij2.matrix[h][i][j] = fij.matrix[h][i][j] + MO_PCM_potential->get(h, I, J);
      }
    }
    for(int i = 0; i < moinfo_.occpi[h]; ++i)
    {
      int I = i + moinfo_.frdocc[h];
      for(int a = 0; a < upper_bound_IA; ++a)
      {
        int A = a + moinfo_.frdocc[h] + moinfo_.occpi[h];
        fIA2.matrix[h][i][a] = fIA.matrix[h][i][a] + MO_PCM_potential->get(h, I, A);
      }
    }
    for(int i = 0; i < upper_bound_ij; ++i)
    {
      int I = i + moinfo_.frdocc[h];
      for(int a = 0; a < upper_bound_IA; ++a)
      {
        int A = a + moinfo_.frdocc[h] + moinfo_.occpi[h];
        fia2.matrix[h][i][a] = fia.matrix[h][i][a] + MO_PCM_potential->get(h, I, A);
      }
    }
    for(int i = 0; i < upper_bound_ij; ++i)
    {
      int I = i + moinfo_.frdocc[h];
      for(int a = 0; a < moinfo_.openpi[h]; ++a)
      {
        int A = a + upper_bound_IA;
        int aA = a + moinfo_.frdocc[h] + upper_bound_ij;
        fia2.matrix[h][i][A] = fia.matrix[h][i][A] + MO_PCM_potential->get(h, I, aA);
      }
    }
    for(int a = 0; a < upper_bound_IA; ++a)
    {
      int A = a + moinfo_.frdocc[h] + moinfo_.occpi[h];
      for(int b = 0; b < upper_bound_IA; ++b)
      {
        int B = b + moinfo_.frdocc[h] + moinfo_.occpi[h];
        fAB2.matrix[h][a][b] = fAB.matrix[h][a][b] + MO_PCM_potential->get(h, A, B);
        fab2.matrix[h][a][b] = fab.matrix[h][a][b] + MO_PCM_potential->get(h, A, B);
      }
    }
    for(int a = 0; a < moinfo_.openpi[h]; ++a)
    {
      int A = a + upper_bound_IA;
      int aA = a + moinfo_.frdocc[h] + upper_bound_ij;
      for(int b = 0; b < moinfo_.openpi[h]; ++b)
      {
        int B = b + upper_bound_IA;
        int bB = b + moinfo_.frdocc[h] + upper_bound_ij;
        fab2.matrix[h][A][B] = fab.matrix[h][A][B] + MO_PCM_potential->get(h, aA, bB);
      }
    }
    for(int a = 0; a < upper_bound_IA; ++a)
    {
      int A = a + moinfo_.frdocc[h] + moinfo_.occpi[h];
      for(int b = 0; b < moinfo_.openpi[h]; ++b)
      {
        int B = b + upper_bound_IA;
        int bB = b + moinfo_.frdocc[h] + upper_bound_ij;
        fab2.matrix[h][a][B] = fab.matrix[h][a][B] + MO_PCM_potential->get(h, A, bB);
      }
    }
    for(int a = 0; a < moinfo_.openpi[h]; ++a)
    {
      int A = a + upper_bound_IA;
      int aA = a + moinfo_.frdocc[h] + upper_bound_ij;
      for(int b = 0; b < upper_bound_IA; ++b)
      {
        int B = b + moinfo_.frdocc[h] + moinfo_.occpi[h];
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

void CCEnergyWavefunction::update_F_pcm_uhf(SharedMatrix &MO_PCM_potential_A, SharedMatrix &MO_PCM_potential_B)
{
  dpdfile2 fIJ, fij, fAB, fab, fIA, fia;
  if (!(_default_psio_lib_->tocentry_exists(PSIF_CC_OEI, "fIJ original"))) // We check just one of the six blocks for existence
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
  for(int h = 0; h < moinfo_.nirreps; ++h)
  {
    for(int i = 0; i < moinfo_.aoccpi[h]; ++i)
    {
      int I = i + moinfo_.frdocc[h];
      for(int j = 0; j < moinfo_.aoccpi[h]; ++j)
      {
        int J = j + moinfo_.frdocc[h];
        fIJ2.matrix[h][i][j] = fIJ.matrix[h][i][j] + MO_PCM_potential_A->get(h, I, J);
      }
    }
    for(int i = 0; i < moinfo_.boccpi[h]; ++i)
    {
      int I = i + moinfo_.frdocc[h];
      for(int j = 0; j < moinfo_.boccpi[h]; ++j)
      {
        int J = j + moinfo_.frdocc[h];
        fij2.matrix[h][i][j] = fij.matrix[h][i][j] + MO_PCM_potential_B->get(h, I, J);
      }
    }
    for(int i = 0; i < moinfo_.aoccpi[h]; ++i)
    {
      int I = i + moinfo_.frdocc[h];
      for(int a = 0; a < moinfo_.avirtpi[h]; ++a)
      {
        int A = a + moinfo_.frdocc[h] + moinfo_.aoccpi[h];
        fIA2.matrix[h][i][a] = fIA.matrix[h][i][a] + MO_PCM_potential_A->get(h, I, A);
      }
    }
    for(int i = 0; i < moinfo_.boccpi[h]; ++i)
    {
      int I = i + moinfo_.frdocc[h];
      for(int a = 0; a < moinfo_.bvirtpi[h]; ++a)
      {
        int A = a + moinfo_.frdocc[h] + moinfo_.boccpi[h];
        fia2.matrix[h][i][a] = fia.matrix[h][i][a] + MO_PCM_potential_B->get(h, I, A);
      }
    }
    for(int a = 0; a < moinfo_.avirtpi[h]; ++a)
    {
      int A = a + moinfo_.frdocc[h] + moinfo_.aoccpi[h];
      for(int b = 0; b < moinfo_.avirtpi[h]; ++b)
      {
        int B = b + moinfo_.frdocc[h] + moinfo_.aoccpi[h];
        fAB2.matrix[h][a][b] = fAB.matrix[h][a][b] + MO_PCM_potential_A->get(h, A, B);
      }
    }
    for(int a = 0; a < moinfo_.bvirtpi[h]; ++a)
    {
      int A = a + moinfo_.frdocc[h] + moinfo_.boccpi[h];
      for(int b = 0; b < moinfo_.bvirtpi[h]; ++b)
      {
        int B = b + moinfo_.frdocc[h] + moinfo_.boccpi[h];
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
}} // namespace psi::ccenergy

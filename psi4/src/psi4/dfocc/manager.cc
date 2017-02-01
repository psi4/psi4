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

#include "psi4/libqt/qt.h"
#include "dfocc.h"
#include "psi4/libciomr/libciomr.h"

using namespace psi;
using namespace std;


namespace psi{ namespace dfoccwave{

//======================================================================
//             OMP2 Manager
//======================================================================
void DFOCC::omp2_manager()
{
        time4grad = 0;// means I will not compute the gradient
	mo_optimized = 0;// means MOs are not optimized
	orbs_already_opt = 0;// means orbitals are not optimized yet.
	orbs_already_sc = 0;// menas orbitals are not semicanonical yet.
        timer_on("DF CC Integrals");
        df_corr();
        trans_corr();
        df_ref();
        trans_ref();
        outfile->Printf("\tNumber of basis functions in the DF-HF basis: %3d\n", nQ_ref);
        outfile->Printf("\tNumber of basis functions in the DF-CC basis: %3d\n", nQ);
        timer_off("DF CC Integrals");

        // memalloc for density intermediates
        Jc = SharedTensor1d(new Tensor1d("DF_BASIS_SCF J_Q", nQ_ref));
        g1Qc = SharedTensor1d(new Tensor1d("DF_BASIS_SCF G1_Q", nQ_ref));
        g1Qt = SharedTensor1d(new Tensor1d("DF_BASIS_SCF G1t_Q", nQ_ref));
        g1Q = SharedTensor1d(new Tensor1d("DF_BASIS_CC G1_Q", nQ));
        g1Qt2 = SharedTensor1d(new Tensor1d("DF_BASIS_CC G1t_Q", nQ));

     if (reference_ == "RESTRICTED") {
        // mem for amplitudes
        cost_ampAA = 0.0;
        cost_ampAA = nocc2AA * nvir2AA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        cost_amp = 3.0 * cost_ampAA;
        memory = Process::environment.get_memory();
        memory_mb = (double)memory/(1024.0 * 1024.0);
        outfile->Printf("\n\tAvailable memory                      : %9.2lf MB \n", memory_mb);
        outfile->Printf("\tMinimum required memory for amplitudes: %9.2lf MB \n", cost_amp);

        // memory requirements
        // DF-HF B(Q,mn)
        cost_ampAA = 0;
        cost_ampAA = nQ_ref * nso2_;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        outfile->Printf("\tMemory requirement for B-HF (Q|mu nu) : %9.2lf MB \n", cost_ampAA);

        // DF-HF B(Q,ab)
        cost_ampAA = 0;
        cost_ampAA = nQ_ref * nvir2AA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        outfile->Printf("\tMemory requirement for B-HF (Q|ab)    : %9.2lf MB \n", cost_ampAA);

        // Cost of Integral transform for DF-HF B(Q,ab)
        cost_ampAA = 0.0;
        cost_ampAA = nQ_ref * nso2_;
        cost_ampAA += nQ_ref * navirA * navirA;
        cost_ampAA += nQ_ref * nso_ * navirA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        outfile->Printf("\tMemory requirement for DF-HF int trans: %9.2lf MB \n", cost_ampAA);

        // DF-CC B(Q,mn)
        cost_ampAA = 0.0;
        cost_ampAA = nQ * nso2_;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        outfile->Printf("\tMemory requirement for B-CC (Q|mu nu) : %9.2lf MB \n", cost_ampAA);

        // DF-CC B(Q,ab)
        cost_ampAA = 0.0;
        cost_ampAA = nQ * navirA * navirA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        outfile->Printf("\tMemory requirement for B-CC (Q|ab)    : %9.2lf MB \n", cost_ampAA);

        // Cost of Integral transform for DF-CC B(Q,ab)
        cost_ampAA = 0.0;
        cost_ampAA = nQ * nso2_;
        cost_ampAA += nQ * navirA * navirA;
        cost_ampAA += nQ * nso_ * navirA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        outfile->Printf("\tMemory requirement for DF-CC int trans: %9.2lf MB \n", cost_ampAA);
     }  // end if (reference_ == "RESTRICTED")

     else if (reference_ == "UNRESTRICTED") {
        // memory requirements
        cost_ampAA = 0.0;
        cost_ampAA = nocc2AA * nvir2AA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        cost_ampBB = nocc2BB * nvir2BB;
        cost_ampBB /= 1024.0 * 1024.0;
        cost_ampBB *= sizeof(double);
        cost_ampAB = nocc2AB * nvir2AB;
        cost_ampAB /= 1024.0 * 1024.0;
        cost_ampAB *= sizeof(double);
        cost_amp = MAX0(cost_ampAA, cost_ampBB);
        cost_amp = MAX0(cost_amp, cost_ampAB);
        cost_amp = 3.0 * cost_amp;
        memory = Process::environment.get_memory();
        memory_mb = (double)memory/(1024.0 * 1024.0);
        outfile->Printf("\n\tAvailable memory                      : %9.2lf MB \n", memory_mb);
        outfile->Printf("\tMinimum required memory for amplitudes: %9.2lf MB \n", cost_amp);
     }  // end else if (reference_ == "UNRESTRICTED")

        // Fock
        fock();

        // QCHF
        if (qchf_ == "TRUE") qchf();

        // ROHF REF
        if (reference == "ROHF") t1_1st_sc();
	t2_1st_sc();
	Emp2L=Emp2;
        EcorrL=Emp2L-Escf;
	Emp2L_old=Emp2;

	outfile->Printf("\n");
	if (reference == "ROHF") outfile->Printf("\tComputing DF-MP2 energy using SCF MOs (DF-ROHF-MP2)... \n");
	else outfile->Printf("\tComputing DF-MP2 energy using SCF MOs (Canonical DF-MP2)... \n");
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	outfile->Printf("\tDF-HF Energy (a.u.)                : %20.14f\n", Escf);
	outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA+Escsmp2BB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SCS-MP2 Total Energy (a.u.)     : %20.14f\n", Escsmp2);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SOS-MP2 Total Energy (a.u.)     : %20.14f\n", Esosmp2);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SCSN-MP2 Total Energy (a.u.)    : %20.14f\n", Escsnmp2);
	if (reference == "ROHF") outfile->Printf("\tDF-MP2 Singles Energy (a.u.)       : %20.14f\n", Emp2_t1);
	if (reference == "ROHF") outfile->Printf("\tDF-MP2 Doubles Energy (a.u.)       : %20.14f\n", Ecorr - Emp2_t1);
	outfile->Printf("\tDF-MP2 Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
	outfile->Printf("\tDF-MP2 Total Energy (a.u.)         : %20.14f\n", Emp2);
	outfile->Printf("\t======================================================================= \n");

	Process::environment.globals["CURRENT ENERGY"] = Emp2;
	Process::environment.globals["MP2 TOTAL ENERGY"] = Emp2;
	Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = Escsmp2;
	Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
	Process::environment.globals["SCSN-MP2 TOTAL ENERGY"] = Escsnmp2;

        Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
        Process::environment.globals["CURRENT CORRELATION ENERGY"] = Emp2 - Escf;
        Process::environment.globals["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
        Process::environment.globals["SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
        Process::environment.globals["SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
        Process::environment.globals["SCSN-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;

        Process::environment.globals["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
        Process::environment.globals["MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA+Emp2BB;

        // S2
        if (comput_s2_ == "TRUE" && reference_ == "UNRESTRICTED") s2_response();

	omp2_opdm();
	omp2_tpdm();
        mp2l_energy();
        separable_tpdm();
	gfock_vo();
	gfock_ov();
        gfock_oo();
        gfock_vv();
	idp();
	mograd();
        occ_iterations();

        // main if
        if (rms_wog <= tol_grad && fabs(DE) >= tol_Eod && regularization == "FALSE") {
           orbs_already_opt = 1;
	   if (conver == 1) outfile->Printf("\n\tOrbitals are optimized now.\n");
	   else if (conver == 0) {
                    outfile->Printf("\n\tMAX MOGRAD did NOT converged, but RMS MOGRAD converged!!!\n");
	            outfile->Printf("\tI will consider the present orbitals as optimized.\n");
           }
	   outfile->Printf("\tTransforming MOs to the semicanonical basis... \n");

	   semi_canonic();
	   outfile->Printf("\tSwitching to the standard DF-MP2 computation... \n");

           trans_corr();
           trans_ref();
           fock();
	   t2_1st_sc();
           conver = 1;
           if (dertype == "FIRST") {
	       omp2_opdm();
	       omp2_tpdm();
               separable_tpdm();
               gfock_vo();
               gfock_ov();
               gfock_oo();
               gfock_vv();
           }
        }// end main if

        else if (rms_wog <= tol_grad && fabs(DE) >= tol_Eod && regularization == "TRUE") {
	   outfile->Printf("\tOrbital gradient converged, but energy did not... \n");
	   outfile->Printf("\tA tighter rms_mograd_convergence tolerance is recommended... \n");
           throw PSIEXCEPTION("A tighter rms_mograd_convergence tolerance is recommended.");
        }

  if (conver == 1) {
        ref_energy();
	mp2_energy();
        if (orbs_already_opt == 1) Emp2L = Emp2;

	outfile->Printf("\n");
	outfile->Printf("\tComputing MP2 energy using optimized MOs... \n");
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	outfile->Printf("\tDF-HF Energy (a.u.)                : %20.14f\n", Escf);
	outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA+Escsmp2BB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SCS-MP2 Total Energy (a.u.)     : %20.14f\n", Escsmp2);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SOS-MP2 Total Energy (a.u.)     : %20.14f\n", Esosmp2);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SCSN-MP2 Total Energy (a.u.)    : %20.14f\n", Escsnmp2);
	outfile->Printf("\tDF-MP2 Correlation Energy (a.u.)   : %20.14f\n", Emp2 - Escf);
	outfile->Printf("\tDF-MP2 Total Energy (a.u.)         : %20.14f\n", Emp2);
	outfile->Printf("\t======================================================================= \n");


	outfile->Printf("\n");
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\t================ DF-OMP2 FINAL RESULTS ================================ \n");
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	outfile->Printf("\tDF-HF Energy (a.u.)                : %20.14f\n", Escf);
	outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SCS-OMP2 Total Energy (a.u.)    : %20.14f\n", Escsmp2);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SOS-OMP2 Total Energy (a.u.)    : %20.14f\n", Esosmp2);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SCSN-OMP2 Total Energy (a.u.)   : %20.14f\n", Escsnmp2);
	outfile->Printf("\tDF-OMP2 Correlation Energy (a.u.)  : %20.14f\n", Emp2L-Escf);
	outfile->Printf("\tEdfomp2 - Eref (a.u.)              : %20.14f\n", Emp2L-Eref);
	outfile->Printf("\tDF-OMP2 Total Energy (a.u.)        : %20.14f\n", Emp2L);
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\n");


	// Set the global variables with the energies
	Process::environment.globals["OMP2 TOTAL ENERGY"] = Emp2L;
	Process::environment.globals["SCS-OMP2 TOTAL ENERGY"] =  Escsmp2;
	Process::environment.globals["SOS-OMP2 TOTAL ENERGY"] =  Esosmp2;
	Process::environment.globals["SCSN-OMP2 TOTAL ENERGY"] = Escsnmp2;
	Process::environment.globals["CURRENT ENERGY"] = Emp2L;
	Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
	Process::environment.globals["CURRENT CORRELATION ENERGY"] = Emp2L-Escf;

        Process::environment.globals["OMP2 CORRELATION ENERGY"] = Emp2L - Escf;
        Process::environment.globals["SCS-OMP2 CORRELATION ENERGY"] =  Escsmp2 - Escf;
        Process::environment.globals["SOS-OMP2 CORRELATION ENERGY"] =  Esosmp2 - Escf;
        Process::environment.globals["SCSN-OMP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;

        // if scs on
	if (do_scs == "TRUE") {
	    if (scs_type_ == "SCS") {
	       Process::environment.globals["CURRENT ENERGY"] = Escsmp2;
	       Process::environment.globals["CURRENT CORRELATION ENERGY"] = Escsmp2 - Escf;
            }

	    else if (scs_type_ == "SCSN") {
	       Process::environment.globals["CURRENT ENERGY"] = Escsnmp2;
	       Process::environment.globals["CURRENT CORRELATION ENERGY"] = Escsnmp2 - Escf;
            }
	}

        // else if sos on
	else if (do_sos == "TRUE") {
	     if (sos_type_ == "SOS") {
	         Process::environment.globals["CURRENT ENERGY"] = Esosmp2;
  	         Process::environment.globals["CURRENT CORRELATION ENERGY"] = Esosmp2 - Escf;
             }
	}

	//if (natorb == "TRUE") nbo();
	//if (occ_orb_energy == "TRUE") semi_canonic();

        // OEPROP
        if (oeprop_ == "TRUE") oeprop();

        // S2
        if (comput_s2_ == "TRUE" && reference_ == "UNRESTRICTED") s2_response();
        //if (comput_s2_ == "TRUE" && reference_ == "UNRESTRICTED") s2_lagrangian();

        // Compute Analytic Gradients
        if (dertype == "FIRST") dfgrad();

	// Save MOs to wfn
	save_mo_to_wfn();

  }// end if (conver == 1)
}// end omp2_manager

//======================================================================
//             MP2 Manager
//======================================================================
void DFOCC::mp2_manager()
{
        time4grad = 0;// means i will not compute the gradient
	mo_optimized = 0;// means MOs are not optimized
        timer_on("DF CC Integrals");
        df_corr();
        if (dertype == "NONE" && oeprop_ == "FALSE" && ekt_ip_ == "FALSE" && comput_s2_ == "FALSE" && qchf_ == "FALSE") {
            trans_mp2();
        }
        else {
            trans_corr();
            df_ref();
            trans_ref();
            outfile->Printf("\tNumber of basis functions in the DF-HF basis: %3d\n", nQ_ref);

            // memalloc for density intermediates
            Jc = SharedTensor1d(new Tensor1d("DF_BASIS_SCF J_Q", nQ_ref));
            g1Qc = SharedTensor1d(new Tensor1d("DF_BASIS_SCF G1_Q", nQ_ref));
            g1Qt = SharedTensor1d(new Tensor1d("DF_BASIS_SCF G1t_Q", nQ_ref));
            g1Q = SharedTensor1d(new Tensor1d("DF_BASIS_CC G1_Q", nQ));
            g1Qt2 = SharedTensor1d(new Tensor1d("DF_BASIS_CC G1t_Q", nQ));
            if (reference == "ROHF") {
                g1Qp = SharedTensor1d(new Tensor1d("DF_BASIS_SCF G1p_Q", nQ_ref));
            }
        }
        outfile->Printf("\tNumber of basis functions in the DF-CC basis: %3d\n", nQ);

        timer_off("DF CC Integrals");

     if (reference_ == "RESTRICTED") {
        // mem for amplitudes
        cost_ampAA = 0.0;
        cost_ampAA = nocc2AA * nvir2AA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        cost_amp = 3.0 * cost_ampAA;
        memory = Process::environment.get_memory();
        memory_mb = (double)memory/(1024.0 * 1024.0);
        outfile->Printf("\n\tAvailable memory                      : %9.2lf MB \n", memory_mb);
        outfile->Printf("\tMinimum required memory for amplitudes: %9.2lf MB \n", cost_amp);

        // memory requirements
        /*
        // DF-HF B(Q,mn)
        cost_ampAA = 0;
        cost_ampAA = nQ_ref * nso2_;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        outfile->Printf("\n\tMemory requirement for B-HF (Q|mu nu) : %9.2lf MB \n", cost_ampAA);

        // DF-HF B(Q,ab)
        cost_ampAA = 0;
        cost_ampAA = nQ_ref * nvir2AA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        outfile->Printf("\tMemory requirement for B-HF (Q|ab)    : %9.2lf MB \n", cost_ampAA);

        // Cost of Integral transform for DF-HF B(Q,ab)
        cost_ampAA = 0.0;
        cost_ampAA = nQ_ref * nso2_;
        cost_ampAA += nQ_ref * navirA * navirA;
        cost_ampAA += nQ_ref * nso_ * navirA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        outfile->Printf("\tMemory requirement for DF-HF int trans: %9.2lf MB \n", cost_ampAA);
        */

        // DF-CC B(Q,mn)
        cost_ampAA = 0.0;
        cost_ampAA = nQ * nso2_;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        outfile->Printf("\tMemory requirement for B-CC (Q|mu nu) : %9.2lf MB \n", cost_ampAA);

        // DF-CC B(Q,ab)
        cost_ampAA = 0.0;
        cost_ampAA = nQ * navirA * navirA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        outfile->Printf("\tMemory requirement for B-CC (Q|ab)    : %9.2lf MB \n", cost_ampAA);

        // Cost of Integral transform for DF-CC B(Q,ab)
        cost_ampAA = 0.0;
        cost_ampAA = nQ * nso2_;
        cost_ampAA += nQ * navirA * navirA;
        cost_ampAA += nQ * nso_ * navirA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        outfile->Printf("\tMemory requirement for DF-CC int trans: %9.2lf MB \n", cost_ampAA);
     }  // end if (reference_ == "RESTRICTED")

     else if (reference_ == "UNRESTRICTED") {
        // memory requirements
        cost_ampAA = 0.0;
        cost_ampAA = nocc2AA * nvir2AA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        cost_ampBB = nocc2BB * nvir2BB;
        cost_ampBB /= 1024.0 * 1024.0;
        cost_ampBB *= sizeof(double);
        cost_ampAB = nocc2AB * nvir2AB;
        cost_ampAB /= 1024.0 * 1024.0;
        cost_ampAB *= sizeof(double);
        cost_amp = MAX0(cost_ampAA, cost_ampBB);
        cost_amp = MAX0(cost_amp, cost_ampAB);
        cost_amp = 3.0 * cost_amp;
        memory = Process::environment.get_memory();
        memory_mb = (double)memory/(1024.0 * 1024.0);
        outfile->Printf("\n\tAvailable memory                      : %9.2lf MB \n", memory_mb);
        outfile->Printf("\tMinimum required memory for amplitudes: %9.2lf MB \n", cost_amp);
     }  // end else if (reference_ == "UNRESTRICTED")

        // QCHF
        if (qchf_ == "TRUE") qchf();

        // ROHF REF
        //outfile->Printf("\tI am here.\n");
        if (reference == "ROHF") t1_1st_sc();
        if (dertype == "NONE" && oeprop_ == "FALSE" && ekt_ip_ == "FALSE" && comput_s2_ == "FALSE") mp2_direct();
        else {
             fock();
             if (mp2_amp_type_ == "DIRECT") mp2_direct();
             else {
	         t2_1st_sc();
                 mp2_energy();
             }
        }
	Emp2L=Emp2;
        EcorrL=Emp2L-Escf;

	outfile->Printf("\n");
	if (reference == "ROHF") outfile->Printf("\tComputing DF-MP2 energy using SCF MOs (DF-ROHF-MP2)... \n");
	else outfile->Printf("\tComputing DF-MP2 energy using SCF MOs (Canonical DF-MP2)... \n");
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	outfile->Printf("\tDF-HF Energy (a.u.)                : %20.14f\n", Escf);
	outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA+Escsmp2BB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SCS-MP2 Total Energy (a.u.)     : %20.14f\n", Escsmp2);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SOS-MP2 Total Energy (a.u.)     : %20.14f\n", Esosmp2);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SCSN-MP2 Total Energy (a.u.)    : %20.14f\n", Escsnmp2);
	if (reference_ == "ROHF") outfile->Printf("\tDF-MP2 Singles Energy (a.u.)       : %20.14f\n", Emp2_t1);
	if (reference_ == "ROHF") outfile->Printf("\tDF-MP2 Doubles Energy (a.u.)       : %20.14f\n", Ecorr - Emp2_t1);
	outfile->Printf("\tDF-MP2 Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
	outfile->Printf("\tDF-MP2 Total Energy (a.u.)         : %20.14f\n", Emp2);
	outfile->Printf("\t======================================================================= \n");

	Process::environment.globals["CURRENT ENERGY"] = Emp2;
	Process::environment.globals["MP2 TOTAL ENERGY"] = Emp2;
	Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = Escsmp2;
	Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
	Process::environment.globals["SCSN-MP2 TOTAL ENERGY"] = Escsnmp2;

        Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
        Process::environment.globals["CURRENT CORRELATION ENERGY"] = Emp2 - Escf;
        Process::environment.globals["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
        Process::environment.globals["SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
        Process::environment.globals["SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
        Process::environment.globals["SCSN-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;

        Process::environment.globals["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
        Process::environment.globals["MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA+Emp2BB;

        // S2
        //if (comput_s2_ == "TRUE" && reference_ == "UNRESTRICTED" && dertype == "NONE") s2_response();
        if (comput_s2_ == "TRUE" && reference_ == "UNRESTRICTED") {
            if (reference == "UHF" || reference == "UKS") s2_response();
        }

        // Compute Analytic Gradients
        if (dertype == "FIRST" || oeprop_ == "TRUE" || ekt_ip_ == "TRUE") {
            outfile->Printf("\n\tComputing unrelaxed response density matrices...\n");

 	    omp2_opdm();
	    omp2_tpdm();
	    //mp2l_energy();
            prepare4grad();
            if (oeprop_ == "TRUE") oeprop();
            if (dertype == "FIRST") dfgrad();
            //if (ekt_ip_ == "TRUE") ekt_ip();

        }// if (dertype == "FIRST" || ekt_ip_ == "TRUE")

}// end mp2_manager

//======================================================================
//             CCSD Manager
//======================================================================
void DFOCC::ccsd_manager()
{

        time4grad = 0;// means i will not compute the gradient
	mo_optimized = 0;// means MOs are not optimized

        timer_on("DF CC Integrals");
        df_corr();
        trans_corr();
        timer_off("DF CC Integrals");
        outfile->Printf("\n\tNumber of basis functions in the DF-CC basis: %3d\n", nQ);

        if (dertype == "FIRST" || oeprop_ == "TRUE" || ekt_ip_ == "TRUE" || qchf_ == "TRUE") {
            timer_on("DF REF Integrals");
            df_ref();
            trans_ref();
            timer_off("DF REF Integrals");
            outfile->Printf("\tNumber of basis functions in the DF-HF basis: %3d\n", nQ_ref);
            Jc = SharedTensor1d(new Tensor1d("DF_BASIS_SCF J_Q", nQ_ref));
        }

        // Memory allocation
        T1c = SharedTensor1d(new Tensor1d("DF_BASIS_CC T1_Q", nQ));

     if (reference_ == "RESTRICTED") {
        t1A = SharedTensor2d(new Tensor2d("T1 <I|A>", naoccA, navirA));
        t1newA = SharedTensor2d(new Tensor2d("New T1 <I|A>", naoccA, navirA));
        FiaA = SharedTensor2d(new Tensor2d("Fint <I|A>", naoccA, navirA));
        FtijA = SharedTensor2d(new Tensor2d("Ftilde <I|J>", naoccA, naoccA));
        FtabA = SharedTensor2d(new Tensor2d("Ftilde <A|B>", navirA, navirA));

        // avaliable mem
        memory = Process::environment.get_memory();
        memory_mb = (double)memory/(1024.0 * 1024.0);
        outfile->Printf("\n\tAvailable memory                      : %9.2lf MB \n", memory_mb);

        // memory requirements

        // DF-CC B(Q,ab) + B(Q,ia) + B(Q,ij)
        cost_df = 0.0;
        cost_df = (navirA * navirA) + (navirA * naoccA) + (naoccA * naoccA);
        cost_df *= nQ;
        cost_df /= 1024.0 * 1024.0;
        cost_df *= sizeof(double);
        outfile->Printf("\tMemory requirement for 3-index ints   : %9.2lf MB \n", cost_df);

        // Cost of Integral transform for B(Q,ab)
        cost_ampAA = 0.0;
        cost_ampAA = nQ * nso2_;
        cost_ampAA += nQ * navirA * navirA;
        cost_ampAA += nQ * nso_ * navirA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        outfile->Printf("\tMemory requirement for DF-CC int trans: %9.2lf MB \n", cost_ampAA);

        // Mem for amplitudes
        cost_ampAA = 0.0;
        cost_ampAA = nocc2AA * nvir2AA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        cost_3amp = 3.0 * cost_ampAA;
        cost_4amp = 4.0 * cost_ampAA;
        cost_5amp = 5.0 * cost_ampAA;
        /*
        if ((cost_5amp+cost_df) <= memory_mb) {
             outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_5amp);
             outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_5amp+cost_df);
             nincore_amp = 5;
             t2_incore = true;
             df_ints_incore = true;
        }
        */
        if ((cost_4amp+cost_df) <= memory_mb) {
             outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_4amp);
             outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_4amp+cost_df);
             nincore_amp = 4;
             t2_incore = true;
             df_ints_incore = true;
        }
        else if ((cost_3amp+cost_df) <= memory_mb) {
             outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_3amp);
             //outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_3amp+cost_df);
             outfile->Printf("\tWarning: T2 amplitudes will be stored on the disk!\n");
             nincore_amp = 3;
             t2_incore = false;
             df_ints_incore = false;
        }
        else if (cost_3amp < memory_mb && cost_df < memory_mb ) {
             outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_3amp);
             outfile->Printf("\tWarning: T2 amplitudes will be stored on the disk!\n");
             nincore_amp = 3;
             t2_incore = false;
             df_ints_incore = false;
        }
        else {
             outfile->Printf("\tWarning: There is NOT enough memory for CC contractions!\n");
             outfile->Printf("\tIncrease memory by                    : %9.2lf MB \n", cost_3amp+cost_df-memory_mb);
             throw PSIEXCEPTION("There is NOT enough memory for CC contractions!");
        }

        // W_abef term
        double cost_amp1 = 0.0;
        cost_amp1 = 2.5 * naoccA * naoccA * navirA * navirA;
        cost_amp1 += nQ * navirA * navirA;
        cost_amp1 /= 1024.0 * 1024.0;
        cost_amp1 *= sizeof(double);
        double cost_amp2 = 0.0;
        cost_amp2 = 1.5 * naoccA * naoccA * navirA * navirA;
        cost_amp2 += 4.0 * nQ * navirA * navirA;
        cost_amp2 /= 1024.0 * 1024.0;
        cost_amp2 *= sizeof(double);
        double cost_amp3 = 0.0;
        cost_amp3 = 2.0 * naoccA * naoccA * navirA * navirA;
        cost_amp3 += 3.0 * nQ * navirA * navirA;
        cost_amp3 += 2.0 * navirA * navirA * navirA;
        cost_amp3 /= 1024.0 * 1024.0;
        cost_amp3 *= sizeof(double);
        cost_amp = MAX0(cost_amp1, cost_amp2);
        cost_amp = MAX0(cost_amp, cost_amp3);
        outfile->Printf("\tMemory requirement for Wabef term (T2): %9.2lf MB \n", cost_amp);

        if (cc_lambda_ == "TRUE") {
            cost_amp1 = navirA * navirA * navirA;
            cost_amp1 /= 1024.0 * 1024.0;
            cost_amp += cost_amp1;
            outfile->Printf("\tMemory requirement for Wefab term (L2): %9.2lf MB \n", cost_amp);
        }

 	// cost_ppl_hm
	//cost_ppl_hm = ntri_abAA * ntri_abAA;
        //cost_ppl_hm /= 1024.0 * 1024.0;
	cost_ppl_hm = ntri_abAA / 1024.0;
        cost_ppl_hm *= cost_ppl_hm;
        cost_ppl_hm *= sizeof(double);
        cost_ppl_hm += cost_amp;
        outfile->Printf("\tMemory for high mem Wabef algorithm   : %9.2lf MB \n", cost_ppl_hm);
	if (cost_ppl_hm > memory_mb && Wabef_type_ == "AUTO") {
	    do_ppl_hm = false;
            outfile->Printf("\tI will use the LOW_MEM Wabef algorithm! \n");
	}
	else if (cost_ppl_hm <= memory_mb && Wabef_type_ == "AUTO") {
	    do_ppl_hm = true;
            outfile->Printf("\tI will use the HIGH_MEM Wabef algorithm! \n");
	}

        // Mem alloc for DF ints
        if (df_ints_incore) {
            bQijA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA));
            bQiaA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA));
            bQabA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA));
            bQijA->read(psio_, PSIF_DFOCC_INTS);
            bQiaA->read(psio_, PSIF_DFOCC_INTS);
            bQabA->read(psio_, PSIF_DFOCC_INTS, true, true);
        }

        //  Malloc
        if (t2_incore) {
            t2 = SharedTensor2d(new Tensor2d("T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
        }

     }  // end if (reference_ == "RESTRICTED")

     else if (reference_ == "UNRESTRICTED") {
        t1A = SharedTensor2d(new Tensor2d("T1 <I|A>", naoccA, navirA));
        t1B = SharedTensor2d(new Tensor2d("T1 <i|a>", naoccB, navirB));
        FiaA = SharedTensor2d(new Tensor2d("Fint <I|A>", naoccA, navirA));
        FiaB = SharedTensor2d(new Tensor2d("Fint <i|a>", naoccB, navirB));
        FtijA = SharedTensor2d(new Tensor2d("Ftilde <I|J>", naoccA, naoccA));
        FtabA = SharedTensor2d(new Tensor2d("Ftilde <A|B>", navirA, navirA));
        FtijB = SharedTensor2d(new Tensor2d("Ftilde <i|j>", naoccB, naoccB));
        FtabB = SharedTensor2d(new Tensor2d("Ftilde <a|b>", navirB, navirB));

        // memory requirements
        cost_ampAA = 0.0;
        cost_ampAA = nocc2AA * nvir2AA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        cost_ampBB = nocc2BB * nvir2BB;
        cost_ampBB /= 1024.0 * 1024.0;
        cost_ampBB *= sizeof(double);
        cost_ampAB = nocc2AB * nvir2AB;
        cost_ampAB /= 1024.0 * 1024.0;
        cost_ampAB *= sizeof(double);
        cost_amp = MAX0(cost_ampAA, cost_ampBB);
        cost_amp = MAX0(cost_amp, cost_ampAB);
        cost_amp = 3.0 * cost_amp;
        memory = Process::environment.get_memory();
        memory_mb = (double)memory/(1024.0 * 1024.0);
        outfile->Printf("\n\tAvailable memory                      : %9.2lf MB \n", memory_mb);
        outfile->Printf("\tMinimum required memory for amplitudes: %9.2lf MB \n", cost_amp);
     }// else if (reference_ == "UNRESTRICTED")

        // memalloc for density intermediates
        if (qchf_ == "TRUE" || dertype == "FIRST") {
            g1Qc = SharedTensor1d(new Tensor1d("DF_BASIS_SCF G1_Q", nQ_ref));
            g1Qt = SharedTensor1d(new Tensor1d("DF_BASIS_SCF G1t_Q", nQ_ref));
            g1Qp = SharedTensor1d(new Tensor1d("DF_BASIS_SCF G1p_Q", nQ_ref));
            g1Q = SharedTensor1d(new Tensor1d("DF_BASIS_CC G1_Q", nQ));
            g1Qt2 = SharedTensor1d(new Tensor1d("DF_BASIS_CC G1t_Q", nQ));
        }

        // QCHF
        if (qchf_ == "TRUE") qchf();

        // Fock
        if (dertype == "FIRST" || oeprop_ == "TRUE" || ekt_ip_ == "TRUE") fock();

        // Compute MP2 energy
        if (reference == "ROHF") t1_1st_sc();
        if (t2_incore) ccsd_mp2();
        else ccsd_mp2_low();

	outfile->Printf("\n");
	if (reference == "ROHF") outfile->Printf("\tComputing DF-MP2 energy (DF-ROHF-MP2)... \n");
	else outfile->Printf("\tComputing DF-MP2 energy ... \n");
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	outfile->Printf("\tDF-HF Energy (a.u.)                : %20.14f\n", Escf);
	outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA+Escsmp2BB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SCS-MP2 Total Energy (a.u.)     : %20.14f\n", Escsmp2);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SOS-MP2 Total Energy (a.u.)     : %20.14f\n", Esosmp2);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SCSN-MP2 Total Energy (a.u.)    : %20.14f\n", Escsnmp2);
	if (reference_ == "ROHF") outfile->Printf("\tDF-MP2 Singles Energy (a.u.)       : %20.14f\n", Emp2_t1);
	if (reference_ == "ROHF") outfile->Printf("\tDF-MP2 Doubles Energy (a.u.)       : %20.14f\n", Ecorr - Emp2_t1);
	outfile->Printf("\tDF-MP2 Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
	outfile->Printf("\tDF-MP2 Total Energy (a.u.)         : %20.14f\n", Emp2);
	outfile->Printf("\t======================================================================= \n");

	Process::environment.globals["MP2 TOTAL ENERGY"] = Emp2;
	Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = Escsmp2;
	Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
	Process::environment.globals["SCSN-MP2 TOTAL ENERGY"] = Escsnmp2;
        Process::environment.globals["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
        Process::environment.globals["SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
        Process::environment.globals["SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
        Process::environment.globals["SCSN-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;
        Process::environment.globals["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
        Process::environment.globals["MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA+Emp2BB;

        // Perform CCSD iterations
        timer_on("CCSD");
        if (t2_incore) ccsd_iterations();
        else ccsd_iterations_low();
        timer_off("CCSD");

	outfile->Printf("\n");
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\t================ CCSD FINAL RESULTS =================================== \n");
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	outfile->Printf("\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	outfile->Printf("\tDF-CCSD Correlation Energy (a.u.)  : %20.14f\n", Ecorr);
	outfile->Printf("\tDF-CCSD Total Energy (a.u.)        : %20.14f\n", Eccsd);
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\n");

	Process::environment.globals["CURRENT ENERGY"] = Eccsd;
        Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
        Process::environment.globals["CURRENT CORRELATION ENERGY"] = Eccsd - Escf;
	Process::environment.globals["CCSD TOTAL ENERGY"] = Eccsd;
        Process::environment.globals["CCSD CORRELATION ENERGY"] = Eccsd - Escf;

        // CCSDL
        if (dertype == "FIRST" || cc_lambda_ == "TRUE") {
	    // memalloc
            if (dertype == "FIRST") {
                GtijA = SharedTensor2d(new Tensor2d("Gtilde Intermediate <I|J>", naoccA, naoccA));
                GtabA = SharedTensor2d(new Tensor2d("Gtilde Intermediate <A|B>", navirA, navirA));
                L1c = SharedTensor1d(new Tensor1d("DF_BASIS_CC L1_Q", nQ));
	        gQt = SharedTensor1d(new Tensor1d("CCSD PDM G_Qt", nQ));
            }

            timer_on("CCSDL");
            if (t2_incore) {
                tstop();
                tstart();
		lambda_title();
                outfile->Printf("\tSolving Lambda amplitude equations...\n");
		ccsdl_iterations();
	    }
            else throw PSIEXCEPTION("There is NOT enough memory for Lambda equations!");
            timer_off("CCSDL");
        }

        // Compute Analytic Gradients
        if (dertype == "FIRST" || ekt_ip_ == "TRUE") {
            tstop();
            tstart();
            pdm_title();

	    // memalloc
	    G1c_ov = SharedTensor2d(new Tensor2d("Correlation OPDM <O|V>", noccA, nvirA));
	    G1c_vo = SharedTensor2d(new Tensor2d("Correlation OPDM <V|O>", nvirA, noccA));

            outfile->Printf("\tComputing unrelaxed response density matrices...\n");
 	    ccsd_opdm();
	    ccsd_tpdm();
	    //ccl_energy();
            prepare4grad();
            if (oeprop_ == "TRUE") oeprop();
            if (dertype == "FIRST") dfgrad();
            //if (ekt_ip_ == "TRUE") ekt_ip();
        }// if (dertype == "FIRST" || ekt_ip_ == "TRUE")

}// end ccsd_manager

//======================================================================
//             CCSD(T) Manager
//======================================================================
void DFOCC::ccsd_t_manager()
{

        time4grad = 0;// means i will not compute the gradient
	mo_optimized = 0;// means MOs are not optimized

        timer_on("DF CC Integrals");
        df_corr();
        trans_corr();
        timer_off("DF CC Integrals");
        outfile->Printf("\n\tNumber of basis functions in the DF-CC basis: %3d\n", nQ);

        if (dertype == "FIRST" || oeprop_ == "TRUE" || ekt_ip_ == "TRUE" || qchf_ == "TRUE") {
            timer_on("DF REF Integrals");
            df_ref();
            trans_ref();
            timer_off("DF REF Integrals");
            outfile->Printf("\tNumber of basis functions in the DF-HF basis: %3d\n", nQ_ref);
            Jc = SharedTensor1d(new Tensor1d("DF_BASIS_SCF J_Q", nQ_ref));
        }

        // Memory allocation
        T1c = SharedTensor1d(new Tensor1d("DF_BASIS_CC T1_Q", nQ));

     if (reference_ == "RESTRICTED") {
        t1A = SharedTensor2d(new Tensor2d("T1 <I|A>", naoccA, navirA));
        t1newA = SharedTensor2d(new Tensor2d("New T1 <I|A>", naoccA, navirA));
        FiaA = SharedTensor2d(new Tensor2d("Fint <I|A>", naoccA, navirA));
        FtijA = SharedTensor2d(new Tensor2d("Ftilde <I|J>", naoccA, naoccA));
        FtabA = SharedTensor2d(new Tensor2d("Ftilde <A|B>", navirA, navirA));

        // avaliable mem
        memory = Process::environment.get_memory();
        memory_mb = (double)memory/(1024.0 * 1024.0);
        outfile->Printf("\n\tAvailable memory                      : %9.2lf MB \n", memory_mb);

        // memory requirements

        // DF-CC B(Q,ab) + B(Q,ia) + B(Q,ij)
        cost_df = 0.0;
        cost_df = (navirA * navirA) + (navirA * naoccA) + (naoccA * naoccA);
        cost_df *= nQ;
        cost_df /= 1024.0 * 1024.0;
        cost_df *= sizeof(double);
        outfile->Printf("\tMemory requirement for 3-index ints   : %9.2lf MB \n", cost_df);

        // Cost of Integral transform for B(Q,ab)
        cost_ampAA = 0.0;
        cost_ampAA = nQ * nso2_;
        cost_ampAA += nQ * navirA * navirA;
        cost_ampAA += nQ * nso_ * navirA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        outfile->Printf("\tMemory requirement for DF-CC int trans: %9.2lf MB \n", cost_ampAA);

        // Mem for amplitudes
        cost_ampAA = 0.0;
        cost_ampAA = nocc2AA * nvir2AA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        cost_3amp = 3.0 * cost_ampAA;
        cost_4amp = 4.0 * cost_ampAA;
        cost_5amp = 5.0 * cost_ampAA;

        if ((cost_4amp+cost_df) <= memory_mb) {
             outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_4amp);
             outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_4amp+cost_df);
             nincore_amp = 4;
             t2_incore = true;
             df_ints_incore = true;
        }
        else if ((cost_3amp+cost_df) <= memory_mb) {
             outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_3amp);
             //outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_3amp+cost_df);
             outfile->Printf("\tWarning: T2 amplitudes will be stored on the disk!\n");
             nincore_amp = 3;
             t2_incore = false;
             df_ints_incore = false;
        }
        else if (cost_3amp < memory_mb && cost_df < memory_mb ) {
             outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_3amp);
             outfile->Printf("\tWarning: T2 amplitudes will be stored on the disk!\n");
             nincore_amp = 3;
             t2_incore = false;
             df_ints_incore = false;
        }
        else {
             outfile->Printf("\tWarning: There is NOT enough memory for CC contractions!\n");
             outfile->Printf("\tIncrease memory by                    : %9.2lf MB \n", cost_3amp+cost_df-memory_mb);
             throw PSIEXCEPTION("There is NOT enough memory for CC contractions!");
        }

        // W_abef term
        double cost_amp1 = 0.0;
        cost_amp1 = 2.5 * naoccA * naoccA * navirA * navirA;
        cost_amp1 += nQ * navirA * navirA;
        cost_amp1 /= 1024.0 * 1024.0;
        cost_amp1 *= sizeof(double);
        double cost_amp2 = 0.0;
        cost_amp2 = 1.5 * naoccA * naoccA * navirA * navirA;
        cost_amp2 += 4.0 * nQ * navirA * navirA;
        cost_amp2 /= 1024.0 * 1024.0;
        cost_amp2 *= sizeof(double);
        double cost_amp3 = 0.0;
        cost_amp3 = 2.0 * naoccA * naoccA * navirA * navirA;
        cost_amp3 += 3.0 * nQ * navirA * navirA;
        cost_amp3 += 2.0 * navirA * navirA * navirA;
        cost_amp3 /= 1024.0 * 1024.0;
        cost_amp3 *= sizeof(double);
        cost_amp = MAX0(cost_amp1, cost_amp2);
        cost_amp = MAX0(cost_amp, cost_amp3);
        outfile->Printf("\tMemory requirement for Wabef term (T2): %9.2lf MB \n", cost_amp);

        if (cc_lambda_ == "TRUE") {
            cost_amp1 = navirA * navirA * navirA;
            cost_amp1 /= 1024.0 * 1024.0;
            cost_amp += cost_amp1;
            outfile->Printf("\tMemory requirement for Wefab term (L2): %9.2lf MB \n", cost_amp);
        }

 	// cost_ppl_hm
	cost_ppl_hm = ntri_abAA / 1024.0;
        cost_ppl_hm *= cost_ppl_hm;
        cost_ppl_hm *= sizeof(double);
        cost_ppl_hm += cost_amp;
        outfile->Printf("\tMemory for high mem Wabef algorithm   : %9.2lf MB \n", cost_ppl_hm);
	if (cost_ppl_hm > memory_mb && Wabef_type_ == "AUTO") {
	    do_ppl_hm = false;
            outfile->Printf("\tI will use the LOW_MEM Wabef algorithm! \n");
	}
	else if (cost_ppl_hm <= memory_mb && Wabef_type_ == "AUTO") {
	    do_ppl_hm = true;
            outfile->Printf("\tI will use the HIGH_MEM Wabef algorithm! \n");
	}

	// Memory for triples: 2*O^2V^2 + 5*V^3 + O^3V + V^2N + V^3/2
        cost_amp1 = 0.0;
        cost_amp1 = 2.0 * naoccA * naoccA * navirA * navirA;
        cost_amp1 += 5.0 * naoccA * navirA * navirA;
        cost_amp1 += naoccA * naoccA * naoccA * navirA;
        cost_amp1 += nQ * navirA * navirA;
        cost_amp1 += navirA * ntri_abAA;
        cost_amp1 /= 1024.0 * 1024.0;
        cost_amp1 *= sizeof(double);
	// Memory: OV^3 + 2*O^2V^2 + 2*V^3 + O^3V + V^2N
        cost_triples_iabc = 0.0;
        cost_triples_iabc = 2.0 * naoccA * naoccA * navirA * navirA;
        cost_triples_iabc += 5.0 * naoccA * navirA * navirA;
        cost_triples_iabc += naoccA * naoccA * naoccA * navirA;
        cost_triples_iabc += nQ * navirA * navirA;
        cost_triples_iabc /= 1024.0 * 1024.0;
        cost_amp2 = 0.0;
        cost_amp2 = (navirA * navirA) / 1024.0;
        cost_amp2 *= (naoccA * navirA) / 1024.0;
        cost_triples_iabc += cost_amp2;
        cost_triples_iabc *= sizeof(double);

	if (triples_iabc_type_ == "DISK") {
	    do_triples_hm = false;
            //outfile->Printf("\n\tI will use a DISK algorithm for (ia|bc) in (T)! \n");
            outfile->Printf("\tMemory requirement for (T) correction : %9.2lf MB \n", cost_amp1);
	}
	else if (triples_iabc_type_ == "DIRECT") {
	    do_triples_hm = false;
            outfile->Printf("\n\tI will use a DIRECT algorithm for (ia|bc) in (T)! \n");
            outfile->Printf("\tMemory requirement for (T) correction : %9.2lf MB \n", cost_amp1);
	}
	else if (cost_triples_iabc > memory_mb && triples_iabc_type_ == "AUTO") {
	    do_triples_hm = false;
            outfile->Printf("\n\tI will use a DIRECT algorithm for (ia|bc) in (T)! \n");
            outfile->Printf("\tMemory requirement for (T) correction : %9.2lf MB \n", cost_amp1);
	}
	else if (cost_triples_iabc <= memory_mb && triples_iabc_type_ == "AUTO") {
	    do_triples_hm = true;
            outfile->Printf("\n\tI will use an INCORE algorithm for (ia|bc) in (T)! \n");
            outfile->Printf("\tMemory requirement for (T) correction : %9.2lf MB \n", cost_triples_iabc);
	}

        // Mem alloc for DF ints
        if (df_ints_incore) {
            bQijA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA));
            bQiaA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA));
            bQabA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA));
            bQijA->read(psio_, PSIF_DFOCC_INTS);
            bQiaA->read(psio_, PSIF_DFOCC_INTS);
            bQabA->read(psio_, PSIF_DFOCC_INTS, true, true);
        }

        //  Malloc
        if (t2_incore) {
            t2 = SharedTensor2d(new Tensor2d("T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
        }

     }  // end if (reference_ == "RESTRICTED")

     else if (reference_ == "UNRESTRICTED") {
        t1A = SharedTensor2d(new Tensor2d("T1 <I|A>", naoccA, navirA));
        t1B = SharedTensor2d(new Tensor2d("T1 <i|a>", naoccB, navirB));
        FiaA = SharedTensor2d(new Tensor2d("Fint <I|A>", naoccA, navirA));
        FiaB = SharedTensor2d(new Tensor2d("Fint <i|a>", naoccB, navirB));
        FtijA = SharedTensor2d(new Tensor2d("Ftilde <I|J>", naoccA, naoccA));
        FtabA = SharedTensor2d(new Tensor2d("Ftilde <A|B>", navirA, navirA));
        FtijB = SharedTensor2d(new Tensor2d("Ftilde <i|j>", naoccB, naoccB));
        FtabB = SharedTensor2d(new Tensor2d("Ftilde <a|b>", navirB, navirB));

        // memory requirements
        cost_ampAA = 0.0;
        cost_ampAA = nocc2AA * nvir2AA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        cost_ampBB = nocc2BB * nvir2BB;
        cost_ampBB /= 1024.0 * 1024.0;
        cost_ampBB *= sizeof(double);
        cost_ampAB = nocc2AB * nvir2AB;
        cost_ampAB /= 1024.0 * 1024.0;
        cost_ampAB *= sizeof(double);
        cost_amp = MAX0(cost_ampAA, cost_ampBB);
        cost_amp = MAX0(cost_amp, cost_ampAB);
        cost_amp = 3.0 * cost_amp;
        memory = Process::environment.get_memory();
        memory_mb = (double)memory/(1024.0 * 1024.0);
        outfile->Printf("\n\tAvailable memory                      : %9.2lf MB \n", memory_mb);
        outfile->Printf("\tMinimum required memory for amplitudes: %9.2lf MB \n", cost_amp);
     }// else if (reference_ == "UNRESTRICTED")

        // memalloc for density intermediates
        if (qchf_ == "TRUE" || dertype == "FIRST") {
            g1Qc = SharedTensor1d(new Tensor1d("DF_BASIS_SCF G1_Q", nQ_ref));
            g1Qt = SharedTensor1d(new Tensor1d("DF_BASIS_SCF G1t_Q", nQ_ref));
            g1Qp = SharedTensor1d(new Tensor1d("DF_BASIS_SCF G1p_Q", nQ_ref));
            g1Q = SharedTensor1d(new Tensor1d("DF_BASIS_CC G1_Q", nQ));
            g1Qt2 = SharedTensor1d(new Tensor1d("DF_BASIS_CC G1t_Q", nQ));
        }

        // QCHF
        if (qchf_ == "TRUE") qchf();

        // Fock
        if (dertype == "FIRST" || oeprop_ == "TRUE" || ekt_ip_ == "TRUE") fock();

        // Compute MP2 energy
        if (reference == "ROHF") t1_1st_sc();
        if (t2_incore) ccsd_mp2();
        else ccsd_mp2_low();

	outfile->Printf("\n");
	if (reference == "ROHF") outfile->Printf("\tComputing DF-MP2 energy (DF-ROHF-MP2)... \n");
	else outfile->Printf("\tComputing DF-MP2 energy ... \n");
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	outfile->Printf("\tDF-HF Energy (a.u.)                : %20.14f\n", Escf);
	outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA+Escsmp2BB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SCS-MP2 Total Energy (a.u.)     : %20.14f\n", Escsmp2);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SOS-MP2 Total Energy (a.u.)     : %20.14f\n", Esosmp2);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SCSN-MP2 Total Energy (a.u.)    : %20.14f\n", Escsnmp2);
	if (reference_ == "ROHF") outfile->Printf("\tDF-MP2 Singles Energy (a.u.)       : %20.14f\n", Emp2_t1);
	if (reference_ == "ROHF") outfile->Printf("\tDF-MP2 Doubles Energy (a.u.)       : %20.14f\n", Ecorr - Emp2_t1);
	outfile->Printf("\tDF-MP2 Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
	outfile->Printf("\tDF-MP2 Total Energy (a.u.)         : %20.14f\n", Emp2);
	outfile->Printf("\t======================================================================= \n");

	Process::environment.globals["MP2 TOTAL ENERGY"] = Emp2;
	Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = Escsmp2;
	Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
	Process::environment.globals["SCSN-MP2 TOTAL ENERGY"] = Escsnmp2;
        Process::environment.globals["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
        Process::environment.globals["SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
        Process::environment.globals["SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
        Process::environment.globals["SCSN-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;
        Process::environment.globals["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
        Process::environment.globals["MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA+Emp2BB;

        // Perform CCSD iterations
        timer_on("CCSD");
        if (t2_incore) ccsd_iterations();
        else ccsd_iterations_low();
        timer_off("CCSD");

	outfile->Printf("\n");
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\t================ CCSD FINAL RESULTS =================================== \n");
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	outfile->Printf("\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	outfile->Printf("\tDF-CCSD Correlation Energy (a.u.)  : %20.14f\n", Ecorr);
	outfile->Printf("\tDF-CCSD Total Energy (a.u.)        : %20.14f\n", Eccsd);
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\n");
	Process::environment.globals["CCSD TOTAL ENERGY"] = Eccsd;
        Process::environment.globals["CCSD CORRELATION ENERGY"] = Eccsd - Escf;

	// CCSD(T)
        tstop();
        tstart();
	pt_title();
        outfile->Printf("\tComputing (T) correction...\n");
        timer_on("(T)");
        if (triples_iabc_type_ == "DISK") ccsd_canonic_triples_disk();
	else if (triples_iabc_type_ == "AUTO") {
	    if (do_triples_hm) ccsd_canonic_triples_hm();
	    else ccsd_canonic_triples();
        }
        else if (triples_iabc_type_ == "INCORE") ccsd_canonic_triples_hm();
        else if (triples_iabc_type_ == "DIRECT") ccsd_canonic_triples();
        else if (triples_iabc_type_ == "DISK") ccsd_canonic_triples_disk();
        timer_off("(T)");
	outfile->Printf("\t(T) Correction (a.u.)              : %20.14f\n", E_t);
	outfile->Printf("\tDF-CCSD(T) Total Energy (a.u.)     : %20.14f\n", Eccsd_t);

	Process::environment.globals["CURRENT ENERGY"] = Eccsd_t;
        Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
        Process::environment.globals["CURRENT CORRELATION ENERGY"] = Eccsd_t - Escf;
	Process::environment.globals["CCSD(T) TOTAL ENERGY"] = Eccsd_t;
	Process::environment.globals["(T) CORRECTION ENERGY"] = E_t;

	/*
        // CCSDL
        if (dertype == "FIRST" || cc_lambda_ == "TRUE") {
	    // memalloc
            if (dertype == "FIRST") {
                GtijA = SharedTensor2d(new Tensor2d("Gtilde Intermediate <I|J>", naoccA, naoccA));
                GtabA = SharedTensor2d(new Tensor2d("Gtilde Intermediate <A|B>", navirA, navirA));
                L1c = SharedTensor1d(new Tensor1d("DF_BASIS_CC L1_Q", nQ));
	        gQt = SharedTensor1d(new Tensor1d("CCSD PDM G_Qt", nQ));
            }

            timer_on("CCSDL");
            if (t2_incore) ccsdl_iterations();
            else throw PSIEXCEPTION("There is NOT enough memory for Lambda equations!");
            timer_off("CCSDL");
        }

        // Compute Analytic Gradients
        if (dertype == "FIRST" || ekt_ip_ == "TRUE") {
	    // memalloc
	    G1c_ov = SharedTensor2d(new Tensor2d("Correlation OPDM <O|V>", noccA, nvirA));
	    G1c_vo = SharedTensor2d(new Tensor2d("Correlation OPDM <V|O>", nvirA, noccA));

            outfile->Printf("\tComputing unrelaxed response density matrices...\n");
 	    ccsd_opdm();
	    ccsd_tpdm();
	    //ccl_energy();
            prepare4grad();
            if (oeprop_ == "TRUE") oeprop();
            if (dertype == "FIRST") dfgrad();
            //if (ekt_ip_ == "TRUE") ekt_ip();
        }// if (dertype == "FIRST" || ekt_ip_ == "TRUE")
	*/

}// end ccsd_t_manager

//======================================================================
//             Lambda-CCSD(T) Manager
//======================================================================
void DFOCC::ccsdl_t_manager()
{

        time4grad = 0;// means i will not compute the gradient
	mo_optimized = 0;// means MOs are not optimized

        timer_on("DF CC Integrals");
        df_corr();
        trans_corr();
        timer_off("DF CC Integrals");
        outfile->Printf("\n\tNumber of basis functions in the DF-CC basis: %3d\n", nQ);

        if (dertype == "FIRST" || oeprop_ == "TRUE" || ekt_ip_ == "TRUE" || qchf_ == "TRUE") {
            timer_on("DF REF Integrals");
            df_ref();
            trans_ref();
            timer_off("DF REF Integrals");
            outfile->Printf("\tNumber of basis functions in the DF-HF basis: %3d\n", nQ_ref);
            Jc = SharedTensor1d(new Tensor1d("DF_BASIS_SCF J_Q", nQ_ref));
        }

        // Memory allocation
        T1c = SharedTensor1d(new Tensor1d("DF_BASIS_CC T1_Q", nQ));

     if (reference_ == "RESTRICTED") {
        t1A = SharedTensor2d(new Tensor2d("T1 <I|A>", naoccA, navirA));
        t1newA = SharedTensor2d(new Tensor2d("New T1 <I|A>", naoccA, navirA));
        FiaA = SharedTensor2d(new Tensor2d("Fint <I|A>", naoccA, navirA));
        FtijA = SharedTensor2d(new Tensor2d("Ftilde <I|J>", naoccA, naoccA));
        FtabA = SharedTensor2d(new Tensor2d("Ftilde <A|B>", navirA, navirA));

        // avaliable mem
        memory = Process::environment.get_memory();
        memory_mb = (double)memory/(1024.0 * 1024.0);
        outfile->Printf("\n\tAvailable memory                      : %9.2lf MB \n", memory_mb);

        // memory requirements

        // DF-CC B(Q,ab) + B(Q,ia) + B(Q,ij)
        cost_df = 0.0;
        cost_df = (navirA * navirA) + (navirA * naoccA) + (naoccA * naoccA);
        cost_df *= nQ;
        cost_df /= 1024.0 * 1024.0;
        cost_df *= sizeof(double);
        outfile->Printf("\tMemory requirement for 3-index ints   : %9.2lf MB \n", cost_df);

        // Cost of Integral transform for B(Q,ab)
        cost_ampAA = 0.0;
        cost_ampAA = nQ * nso2_;
        cost_ampAA += nQ * navirA * navirA;
        cost_ampAA += nQ * nso_ * navirA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        outfile->Printf("\tMemory requirement for DF-CC int trans: %9.2lf MB \n", cost_ampAA);

        // Mem for amplitudes
        cost_ampAA = 0.0;
        cost_ampAA = nocc2AA * nvir2AA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        cost_3amp = 3.0 * cost_ampAA;
        cost_4amp = 4.0 * cost_ampAA;
        cost_5amp = 5.0 * cost_ampAA;
        /*
        if ((cost_5amp+cost_df) <= memory_mb) {
             outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_5amp);
             outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_5amp+cost_df);
             nincore_amp = 5;
             t2_incore = true;
             df_ints_incore = true;
        }
        */
        if ((cost_4amp+cost_df) <= memory_mb) {
             outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_4amp);
             outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_4amp+cost_df);
             nincore_amp = 4;
             t2_incore = true;
             df_ints_incore = true;
        }
        else if ((cost_3amp+cost_df) <= memory_mb) {
             outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_3amp);
             //outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_3amp+cost_df);
             outfile->Printf("\tWarning: T2 amplitudes will be stored on the disk!\n");
             nincore_amp = 3;
             t2_incore = false;
             df_ints_incore = false;
        }
        else if (cost_3amp < memory_mb && cost_df < memory_mb ) {
             outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_3amp);
             outfile->Printf("\tWarning: T2 amplitudes will be stored on the disk!\n");
             nincore_amp = 3;
             t2_incore = false;
             df_ints_incore = false;
        }
        else {
             outfile->Printf("\tWarning: There is NOT enough memory for CC contractions!\n");
             outfile->Printf("\tIncrease memory by                    : %9.2lf MB \n", cost_3amp+cost_df-memory_mb);
             throw PSIEXCEPTION("There is NOT enough memory for CC contractions!");
        }

        // W_abef term
        double cost_amp1 = 0.0;
        cost_amp1 = 2.5 * naoccA * naoccA * navirA * navirA;
        cost_amp1 += nQ * navirA * navirA;
        cost_amp1 /= 1024.0 * 1024.0;
        cost_amp1 *= sizeof(double);
        double cost_amp2 = 0.0;
        cost_amp2 = 1.5 * naoccA * naoccA * navirA * navirA;
        cost_amp2 += 4.0 * nQ * navirA * navirA;
        cost_amp2 /= 1024.0 * 1024.0;
        cost_amp2 *= sizeof(double);
        double cost_amp3 = 0.0;
        cost_amp3 = 2.0 * naoccA * naoccA * navirA * navirA;
        cost_amp3 += 3.0 * nQ * navirA * navirA;
        cost_amp3 += 2.0 * navirA * navirA * navirA;
        cost_amp3 /= 1024.0 * 1024.0;
        cost_amp3 *= sizeof(double);
        cost_amp = MAX0(cost_amp1, cost_amp2);
        cost_amp = MAX0(cost_amp, cost_amp3);
        outfile->Printf("\tMemory requirement for Wabef term (T2): %9.2lf MB \n", cost_amp);

        if (cc_lambda_ == "TRUE") {
            cost_amp1 = navirA * navirA * navirA;
            cost_amp1 /= 1024.0 * 1024.0;
            cost_amp += cost_amp1;
            outfile->Printf("\tMemory requirement for Wefab term (L2): %9.2lf MB \n", cost_amp);
        }

 	// cost_ppl_hm
	//cost_ppl_hm = ntri_abAA * ntri_abAA;
        //cost_ppl_hm /= 1024.0 * 1024.0;
	cost_ppl_hm = ntri_abAA / 1024.0;
        cost_ppl_hm *= cost_ppl_hm;
        cost_ppl_hm *= sizeof(double);
        cost_ppl_hm += cost_amp;
        outfile->Printf("\tMemory for high mem Wabef algorithm   : %9.2lf MB \n", cost_ppl_hm);
	if (cost_ppl_hm > memory_mb && Wabef_type_ == "AUTO") {
	    do_ppl_hm = false;
            outfile->Printf("\tI will use the LOW_MEM Wabef algorithm! \n");
	}
	else if (cost_ppl_hm <= memory_mb && Wabef_type_ == "AUTO") {
	    do_ppl_hm = true;
            outfile->Printf("\tI will use the HIGH_MEM Wabef algorithm! \n");
	}

	// Memory for triples: 3*O^2V^2 + 5*V^3 + O^3V + V^2N + V^3/2
        cost_amp1 = 0.0;
        cost_amp1 = 3.0 * naoccA * naoccA * navirA * navirA;
        cost_amp1 += 5.0 * naoccA * navirA * navirA;
        cost_amp1 += naoccA * naoccA * naoccA * navirA;
        cost_amp1 += nQ * navirA * navirA;
        cost_amp1 += navirA * ntri_abAA;
        cost_amp1 /= 1024.0 * 1024.0;
        cost_amp1 *= sizeof(double);
        outfile->Printf("\tMemory requirement for (AT) correction : %9.2lf MB \n", cost_amp1);

	/*
	// Memory: OV^3 + 2*O^2V^2 + 2*V^3 + O^3V + V^2N
        cost_triples_iabc = 0.0;
        cost_triples_iabc = 2.0 * naoccA * naoccA * navirA * navirA;
        cost_triples_iabc += 5.0 * naoccA * navirA * navirA;
        cost_triples_iabc += naoccA * naoccA * naoccA * navirA;
        cost_triples_iabc += nQ * navirA * navirA;
        cost_triples_iabc /= 1024.0 * 1024.0;
        cost_amp2 = 0.0;
        cost_amp2 = (navirA * navirA) / 1024.0;
        cost_amp2 *= (naoccA * navirA) / 1024.0;
        cost_triples_iabc += cost_amp2;
        cost_triples_iabc *= sizeof(double);

	if (triples_iabc_type_ == "DISK") {
	    do_triples_hm = false;
            //outfile->Printf("\n\tI will use a DISK algorithm for (ia|bc) in (T)! \n");
            outfile->Printf("\tMemory requirement for (T) correction : %9.2lf MB \n", cost_amp1);
	}
	else if (triples_iabc_type_ == "DIRECT") {
	    do_triples_hm = false;
            outfile->Printf("\n\tI will use a DIRECT algorithm for (ia|bc) in (T)! \n");
            outfile->Printf("\tMemory requirement for (T) correction : %9.2lf MB \n", cost_amp1);
	}
	else if (cost_triples_iabc > memory_mb && triples_iabc_type_ == "AUTO") {
	    do_triples_hm = false;
            outfile->Printf("\n\tI will use a DIRECT algorithm for (ia|bc) in (T)! \n");
            outfile->Printf("\tMemory requirement for (T) correction : %9.2lf MB \n", cost_amp1);
	}
	else if (cost_triples_iabc <= memory_mb && triples_iabc_type_ == "AUTO") {
	    do_triples_hm = true;
            outfile->Printf("\n\tI will use an INCORE algorithm for (ia|bc) in (T)! \n");
            outfile->Printf("\tMemory requirement for (T) correction : %9.2lf MB \n", cost_triples_iabc);
	}
	*/


        // Mem alloc for DF ints
        if (df_ints_incore) {
            bQijA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA));
            bQiaA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA));
            bQabA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA));
            bQijA->read(psio_, PSIF_DFOCC_INTS);
            bQiaA->read(psio_, PSIF_DFOCC_INTS);
            bQabA->read(psio_, PSIF_DFOCC_INTS, true, true);
        }

        //  Malloc
        if (t2_incore) {
            t2 = SharedTensor2d(new Tensor2d("T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
        }

     }  // end if (reference_ == "RESTRICTED")

     else if (reference_ == "UNRESTRICTED") {
        t1A = SharedTensor2d(new Tensor2d("T1 <I|A>", naoccA, navirA));
        t1B = SharedTensor2d(new Tensor2d("T1 <i|a>", naoccB, navirB));
        FiaA = SharedTensor2d(new Tensor2d("Fint <I|A>", naoccA, navirA));
        FiaB = SharedTensor2d(new Tensor2d("Fint <i|a>", naoccB, navirB));
        FtijA = SharedTensor2d(new Tensor2d("Ftilde <I|J>", naoccA, naoccA));
        FtabA = SharedTensor2d(new Tensor2d("Ftilde <A|B>", navirA, navirA));
        FtijB = SharedTensor2d(new Tensor2d("Ftilde <i|j>", naoccB, naoccB));
        FtabB = SharedTensor2d(new Tensor2d("Ftilde <a|b>", navirB, navirB));

        // memory requirements
        cost_ampAA = 0.0;
        cost_ampAA = nocc2AA * nvir2AA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        cost_ampBB = nocc2BB * nvir2BB;
        cost_ampBB /= 1024.0 * 1024.0;
        cost_ampBB *= sizeof(double);
        cost_ampAB = nocc2AB * nvir2AB;
        cost_ampAB /= 1024.0 * 1024.0;
        cost_ampAB *= sizeof(double);
        cost_amp = MAX0(cost_ampAA, cost_ampBB);
        cost_amp = MAX0(cost_amp, cost_ampAB);
        cost_amp = 3.0 * cost_amp;
        memory = Process::environment.get_memory();
        memory_mb = (double)memory/(1024.0 * 1024.0);
        outfile->Printf("\n\tAvailable memory                      : %9.2lf MB \n", memory_mb);
        outfile->Printf("\tMinimum required memory for amplitudes: %9.2lf MB \n", cost_amp);
     }// else if (reference_ == "UNRESTRICTED")

        // memalloc for density intermediates
        if (qchf_ == "TRUE" || dertype == "FIRST") {
            g1Qc = SharedTensor1d(new Tensor1d("DF_BASIS_SCF G1_Q", nQ_ref));
            g1Qt = SharedTensor1d(new Tensor1d("DF_BASIS_SCF G1t_Q", nQ_ref));
            g1Qp = SharedTensor1d(new Tensor1d("DF_BASIS_SCF G1p_Q", nQ_ref));
            g1Q = SharedTensor1d(new Tensor1d("DF_BASIS_CC G1_Q", nQ));
            g1Qt2 = SharedTensor1d(new Tensor1d("DF_BASIS_CC G1t_Q", nQ));
        }

        // QCHF
        if (qchf_ == "TRUE") qchf();

        // Fock
        if (dertype == "FIRST" || oeprop_ == "TRUE" || ekt_ip_ == "TRUE") fock();

        // Compute MP2 energy
        if (reference == "ROHF") t1_1st_sc();
        if (t2_incore) ccsd_mp2();
        else ccsd_mp2_low();

	outfile->Printf("\n");
	if (reference == "ROHF") outfile->Printf("\tComputing DF-MP2 energy (DF-ROHF-MP2)... \n");
	else outfile->Printf("\tComputing DF-MP2 energy ... \n");
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	outfile->Printf("\tDF-HF Energy (a.u.)                : %20.14f\n", Escf);
	outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA+Escsmp2BB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SCS-MP2 Total Energy (a.u.)     : %20.14f\n", Escsmp2);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SOS-MP2 Total Energy (a.u.)     : %20.14f\n", Esosmp2);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SCSN-MP2 Total Energy (a.u.)    : %20.14f\n", Escsnmp2);
	if (reference_ == "ROHF") outfile->Printf("\tDF-MP2 Singles Energy (a.u.)       : %20.14f\n", Emp2_t1);
	if (reference_ == "ROHF") outfile->Printf("\tDF-MP2 Doubles Energy (a.u.)       : %20.14f\n", Ecorr - Emp2_t1);
	outfile->Printf("\tDF-MP2 Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
	outfile->Printf("\tDF-MP2 Total Energy (a.u.)         : %20.14f\n", Emp2);
	outfile->Printf("\t======================================================================= \n");

	Process::environment.globals["MP2 TOTAL ENERGY"] = Emp2;
	Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = Escsmp2;
	Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
	Process::environment.globals["SCSN-MP2 TOTAL ENERGY"] = Escsnmp2;
        Process::environment.globals["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
        Process::environment.globals["SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
        Process::environment.globals["SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
        Process::environment.globals["SCSN-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;
        Process::environment.globals["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
        Process::environment.globals["MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA+Emp2BB;

        // Perform CCSD iterations
        timer_on("CCSD");
        if (t2_incore) ccsd_iterations();
        else ccsd_iterations_low();
        timer_off("CCSD");

	outfile->Printf("\n");
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\t================ CCSD FINAL RESULTS =================================== \n");
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	outfile->Printf("\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	outfile->Printf("\tDF-CCSD Correlation Energy (a.u.)  : %20.14f\n", Ecorr);
	outfile->Printf("\tDF-CCSD Total Energy (a.u.)        : %20.14f\n", Eccsd);
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\n");
	Process::environment.globals["CCSD TOTAL ENERGY"] = Eccsd;
        Process::environment.globals["CCSD CORRELATION ENERGY"] = Eccsd - Escf;

        // CCSDL
        timer_on("CCSDL");
        if (t2_incore) {
            tstop();
            tstart();
            lambda_title();
            outfile->Printf("\tSolving Lambda amplitude equations...\n");
            ccsdl_iterations();
	}
        else throw PSIEXCEPTION("There is NOT enough memory for Lambda equations!");
        timer_off("CCSDL");

	// CCSD(AT)
        tstop();
        tstart();
	pat_title();
        outfile->Printf("\tComputing asymmetric triples (AT) correction...\n");
        timer_on("(AT)");
        ccsdl_canonic_triples_disk();
        timer_off("(AT)");
	outfile->Printf("\t(AT) Correction (a.u.)             : %20.14f\n", E_at);
	outfile->Printf("\tDF-CCSD(AT) Total Energy (a.u.)    : %20.14f\n", Eccsd_at);

	Process::environment.globals["CURRENT ENERGY"] = Eccsd_at;
        Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
        Process::environment.globals["CURRENT CORRELATION ENERGY"] = Eccsd_at - Escf;
	Process::environment.globals["CCSD(AT) TOTAL ENERGY"] = Eccsd_at;
	Process::environment.globals["(AT) CORRECTION ENERGY"] = E_at;

	/*
        // Compute Analytic Gradients
        if (dertype == "FIRST" || ekt_ip_ == "TRUE") {
            tstop();
            tstart();
            pdm_title();

	    // memalloc
	    G1c_ov = SharedTensor2d(new Tensor2d("Correlation OPDM <O|V>", noccA, nvirA));
	    G1c_vo = SharedTensor2d(new Tensor2d("Correlation OPDM <V|O>", nvirA, noccA));

            outfile->Printf("\tComputing unrelaxed response density matrices...\n");
 	    ccsd_opdm();
	    ccsd_tpdm();
	    //ccl_energy();
            prepare4grad();
            if (oeprop_ == "TRUE") oeprop();
            if (dertype == "FIRST") dfgrad();
            //if (ekt_ip_ == "TRUE") ekt_ip();
        }// if (dertype == "FIRST" || ekt_ip_ == "TRUE")
	*/

}// end ccsdl_t_manager

//======================================================================
//             CCD Manager
//======================================================================
void DFOCC::ccd_manager()
{

        time4grad = 0;// means i will not compute the gradient
	mo_optimized = 0;// means MOs are not optimized

        timer_on("DF CC Integrals");
        df_corr();
        trans_corr();
        timer_off("DF CC Integrals");
        outfile->Printf("\n\tNumber of basis functions in the DF-CC basis: %3d\n", nQ);

        if (dertype == "FIRST" || oeprop_ == "TRUE" || ekt_ip_ == "TRUE" || qchf_ == "TRUE") {
            timer_on("DF REF Integrals");
            df_ref();
            trans_ref();
            timer_off("DF REF Integrals");
            outfile->Printf("\tNumber of basis functions in the DF-HF basis: %3d\n", nQ_ref);
            Jc = SharedTensor1d(new Tensor1d("DF_BASIS_SCF J_Q", nQ_ref));
        }

        // Memory allocation
        //T1c = SharedTensor1d(new Tensor1d("DF_BASIS_CC T1_Q", nQ));

     if (reference_ == "RESTRICTED") {
        // avaliable mem
        memory = Process::environment.get_memory();
        memory_mb = (double)memory/(1024.0 * 1024.0);
        outfile->Printf("\n\tAvailable memory                      : %9.2lf MB \n", memory_mb);

        // memory requirements
        // DF-CC B(Q,ab) + B(Q,ia) + B(Q,ij)
        cost_df = 0.0;
        cost_df = (navirA * navirA) + (navirA * naoccA) + (naoccA * naoccA);
        cost_df *= nQ;
        cost_df /= 1024.0 * 1024.0;
        cost_df *= sizeof(double);
        outfile->Printf("\tMemory requirement for 3-index ints   : %9.2lf MB \n", cost_df);

        // Cost of Integral transform for B(Q,ab)
        cost_ampAA = 0.0;
        cost_ampAA = nQ * nso2_;
        cost_ampAA += nQ * navirA * navirA;
        cost_ampAA += nQ * nso_ * navirA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        outfile->Printf("\tMemory requirement for DF-CC int trans: %9.2lf MB \n", cost_ampAA);

        // Mem for amplitudes
        cost_ampAA = 0.0;
        cost_ampAA = nocc2AA * nvir2AA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        cost_3amp = 3.0 * cost_ampAA;
        cost_4amp = 4.0 * cost_ampAA;
        cost_5amp = 5.0 * cost_ampAA;
        /*
        if ((cost_5amp+cost_df) <= memory_mb) {
             outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_5amp);
             outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_5amp+cost_df);
             nincore_amp = 5;
             t2_incore = true;
             df_ints_incore = true;
        }
        */
        if ((cost_4amp+cost_df) <= memory_mb) {
             outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_4amp);
             outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_4amp+cost_df);
             nincore_amp = 4;
             t2_incore = true;
             df_ints_incore = true;
        }
        else if ((cost_3amp+cost_df) <= memory_mb) {
             outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_3amp);
             //outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_3amp+cost_df);
             outfile->Printf("\tWarning: T2 amplitudes will be stored on the disk!\n");
             nincore_amp = 3;
             t2_incore = false;
             df_ints_incore = false;
        }
        else if (cost_3amp < memory_mb && cost_df < memory_mb ) {
             outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_3amp);
             outfile->Printf("\tWarning: T2 amplitudes will be stored on the disk!\n");
             nincore_amp = 3;
             t2_incore = false;
             df_ints_incore = false;
        }
        else {
             outfile->Printf("\tWarning: There is NOT enough memory for CC contractions!\n");
             outfile->Printf("\tIncrease memory by                    : %9.2lf MB \n", cost_3amp+cost_df-memory_mb);
             throw PSIEXCEPTION("There is NOT enough memory for CC contractions!");
        }

        // W_abef term
        cost_ampAA = 0.0;
        cost_ampAA = naoccA * naoccA * navirA * navirA;
        cost_ampAA += 2.0 * nQ * navirA * navirA;
        cost_ampAA += navirA * navirA * navirA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        double cost_ampAA2 = 0.0;
        cost_ampAA2 = naoccA * naoccA * navirA * navirA;
        cost_ampAA2 += nQ * navirA * navirA;
        cost_ampAA2 += 3.0 * navirA * navirA * navirA;
        cost_ampAA2 /= 1024.0 * 1024.0;
        cost_ampAA2 *= sizeof(double);
        cost_amp = MAX0(cost_ampAA, cost_ampAA2);
        outfile->Printf("\tMemory requirement for Wabef term     : %9.2lf MB \n", cost_amp);

	// cost_ppl_hm
	cost_ppl_hm = ntri_abAA / 1024.0;
        cost_ppl_hm *= cost_ppl_hm;
        cost_ppl_hm *= sizeof(double);
        cost_ppl_hm += cost_amp;
        outfile->Printf("\tMemory for high mem Wabef algorithm   : %9.2lf MB \n", cost_ppl_hm);
	if (cost_ppl_hm > memory_mb && Wabef_type_ == "AUTO") {
	    do_ppl_hm = false;
            outfile->Printf("\tI will use the LOW_MEM Wabef algorithm! \n");
	}
	else if (cost_ppl_hm <= memory_mb && Wabef_type_ == "AUTO") {
	    do_ppl_hm = true;
            outfile->Printf("\tI will use the HIGH_MEM Wabef algorithm! \n");
	}

        // Mem alloc for DF ints
        if (df_ints_incore) {
            bQijA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA));
            bQiaA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA));
            bQabA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA));
            bQijA->read(psio_, PSIF_DFOCC_INTS);
            bQiaA->read(psio_, PSIF_DFOCC_INTS);
            bQabA->read(psio_, PSIF_DFOCC_INTS, true, true);
        }

        //  Malloc
        if (t2_incore) {
            t2 = SharedTensor2d(new Tensor2d("T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
        }

     }  // end if (reference_ == "RESTRICTED")

     else if (reference_ == "UNRESTRICTED") {
        // memory requirements
        cost_ampAA = 0.0;
        cost_ampAA = nocc2AA * nvir2AA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        cost_ampBB = nocc2BB * nvir2BB;
        cost_ampBB /= 1024.0 * 1024.0;
        cost_ampBB *= sizeof(double);
        cost_ampAB = nocc2AB * nvir2AB;
        cost_ampAB /= 1024.0 * 1024.0;
        cost_ampAB *= sizeof(double);
        cost_amp = MAX0(cost_ampAA, cost_ampBB);
        cost_amp = MAX0(cost_amp, cost_ampAB);
        cost_amp = 3.0 * cost_amp;
        memory = Process::environment.get_memory();
        memory_mb = (double)memory/(1024.0 * 1024.0);
        outfile->Printf("\n\tAvailable memory                      : %9.2lf MB \n", memory_mb);
        outfile->Printf("\tMinimum required memory for amplitudes: %9.2lf MB \n", cost_amp);
     }// else if (reference_ == "UNRESTRICTED")

        // memalloc for density intermediates
        if (qchf_ == "TRUE" || dertype == "FIRST") {
            g1Qc = SharedTensor1d(new Tensor1d("DF_BASIS_SCF G1_Q", nQ_ref));
            g1Qt = SharedTensor1d(new Tensor1d("DF_BASIS_SCF G1t_Q", nQ_ref));
            g1Qp = SharedTensor1d(new Tensor1d("DF_BASIS_SCF G1p_Q", nQ_ref));
            g1Q = SharedTensor1d(new Tensor1d("DF_BASIS_CC G1_Q", nQ));
            g1Qt2 = SharedTensor1d(new Tensor1d("DF_BASIS_CC G1t_Q", nQ));
        }

        // QCHF
        if (qchf_ == "TRUE") qchf();

        // Fock
        if (dertype == "FIRST" || oeprop_ == "TRUE" || ekt_ip_ == "TRUE") fock();

        // Compute MP2 energy
        if (reference == "ROHF") t1_1st_sc();
        if (t2_incore) ccd_mp2();
        else ccd_mp2_low();

	outfile->Printf("\n");
	if (reference == "ROHF") outfile->Printf("\tComputing DF-MP2 energy (DF-ROHF-MP2)... \n");
	else outfile->Printf("\tComputing DF-MP2 energy ... \n");
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	outfile->Printf("\tDF-HF Energy (a.u.)                : %20.14f\n", Escf);
	outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA+Escsmp2BB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SCS-MP2 Total Energy (a.u.)     : %20.14f\n", Escsmp2);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SOS-MP2 Total Energy (a.u.)     : %20.14f\n", Esosmp2);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SCSN-MP2 Total Energy (a.u.)    : %20.14f\n", Escsnmp2);
	if (reference_ == "ROHF") outfile->Printf("\tDF-MP2 Singles Energy (a.u.)       : %20.14f\n", Emp2_t1);
	if (reference_ == "ROHF") outfile->Printf("\tDF-MP2 Doubles Energy (a.u.)       : %20.14f\n", Ecorr - Emp2_t1);
	outfile->Printf("\tDF-MP2 Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
	outfile->Printf("\tDF-MP2 Total Energy (a.u.)         : %20.14f\n", Emp2);
	outfile->Printf("\t======================================================================= \n");

	Process::environment.globals["MP2 TOTAL ENERGY"] = Emp2;
	Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = Escsmp2;
	Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
	Process::environment.globals["SCSN-MP2 TOTAL ENERGY"] = Escsnmp2;
        Process::environment.globals["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
        Process::environment.globals["SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
        Process::environment.globals["SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
        Process::environment.globals["SCSN-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;
        Process::environment.globals["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
        Process::environment.globals["MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA+Emp2BB;

        // Perform CCD iterations
        timer_on("CCD");
        if (t2_incore) ccd_iterations();
        else ccd_iterations_low();
        timer_off("CCD");

	outfile->Printf("\n");
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\t================ CCD FINAL RESULTS ==================================== \n");
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	outfile->Printf("\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	outfile->Printf("\tDF-CCD Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
	outfile->Printf("\tDF-CCD Total Energy (a.u.)         : %20.14f\n", Eccd);
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\n");

	Process::environment.globals["CURRENT ENERGY"] = Eccd;
        Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
        Process::environment.globals["CURRENT CORRELATION ENERGY"] = Eccd - Escf;
	Process::environment.globals["CCD TOTAL ENERGY"] = Eccd;
        Process::environment.globals["CCD CORRELATION ENERGY"] = Eccd - Escf;

        // CCDL
        if (dertype == "FIRST" || cc_lambda_ == "TRUE") {
            // memalloc
            if (dertype == "FIRST") {
	        gQt = SharedTensor1d(new Tensor1d("CCD PDM G_Qt", nQ));
            }

            timer_on("CCDL");
            if (t2_incore) {
                tstop();
                tstart();
		lambda_title();
                outfile->Printf("\tSolving Lambda amplitude equations...\n");
		ccdl_iterations();
	    }
            else throw PSIEXCEPTION("There is NOT enough memory for Lambda equations!");
            timer_off("CCDL");
        }

        // Compute Analytic Gradients
        if (dertype == "FIRST" || ekt_ip_ == "TRUE") {
            tstop();
            tstart();
            pdm_title();

	    // memalloc
	    G1c_ov = SharedTensor2d(new Tensor2d("Correlation OPDM <O|V>", noccA, nvirA));
	    G1c_vo = SharedTensor2d(new Tensor2d("Correlation OPDM <V|O>", nvirA, noccA));

            outfile->Printf("\tComputing unrelaxed response density matrices...\n");
 	    ccd_opdm();
	    ccd_tpdm();
	    //ccl_energy();
            prepare4grad();
            if (oeprop_ == "TRUE") oeprop();
            if (dertype == "FIRST") dfgrad();
            //if (ekt_ip_ == "TRUE") ekt_ip();
        }// if (dertype == "FIRST" || ekt_ip_ == "TRUE")

}// end ccd_manager

//======================================================================
//             OMP3 Manager
//======================================================================
void DFOCC::omp3_manager()
{

        time4grad = 0;// means I will not compute the gradient
	mo_optimized = 0;// means MOs are not optimized
	orbs_already_opt = 0;// means orbitals are not optimized yet.
	orbs_already_sc = 0;// means orbitals are not semicanonical yet.

        timer_on("DF CC Integrals");
        df_corr();
        trans_corr();
        df_ref();
        trans_ref();
        outfile->Printf("\tNumber of basis functions in the DF-HF basis: %3d\n", nQ_ref);
        outfile->Printf("\tNumber of basis functions in the DF-CC basis: %3d\n", nQ);
        timer_off("DF CC Integrals");

        // memalloc for density intermediates
        Jc = SharedTensor1d(new Tensor1d("DF_BASIS_SCF J_Q", nQ_ref));
        g1Qc = SharedTensor1d(new Tensor1d("DF_BASIS_SCF G1_Q", nQ_ref));
        g1Qt = SharedTensor1d(new Tensor1d("DF_BASIS_SCF G1t_Q", nQ_ref));
        g1Qp = SharedTensor1d(new Tensor1d("DF_BASIS_SCF G1p_Q", nQ_ref));
        g1Q = SharedTensor1d(new Tensor1d("DF_BASIS_CC G1_Q", nQ));
        g1Qt2 = SharedTensor1d(new Tensor1d("DF_BASIS_CC G1t_Q", nQ));

        // avaliable mem
        memory = Process::environment.get_memory();
        memory_mb = (double)memory/(1024.0 * 1024.0);
        outfile->Printf("\n\tAvailable memory                      : %9.2lf MB \n", memory_mb);

        // memory requirements
        // DF-CC B(Q,ab) + B(Q,ia) + B(Q,ij)
        cost_df = 0.0;
        cost_df = (navirA * navirA) + (navirA * naoccA) + (naoccA * naoccA);
        cost_df *= nQ;
        cost_df /= 1024.0 * 1024.0;
        cost_df *= sizeof(double);
        if (reference_ == "RESTRICTED") outfile->Printf("\tMemory requirement for 3-index ints   : %9.2lf MB \n", cost_df);
	else if (reference_ == "UNRESTRICTED") outfile->Printf("\tMemory requirement for 3-index ints   : %9.2lf MB \n", 2.0*cost_df);

        // Cost of Integral transform for B(Q,ab)
        cost_ampAA = 0.0;
        cost_ampAA = nQ * nso2_;
        cost_ampAA += nQ * navirA * navirA;
        cost_ampAA += nQ * nso_ * navirA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        outfile->Printf("\tMemory requirement for DF-CC int trans: %9.2lf MB \n", cost_ampAA);

        // Mem for amplitudes
        cost_ampAA = 0.0;
        cost_ampAA = nocc2AA * nvir2AA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        cost_3amp = 3.0 * cost_ampAA;
        cost_4amp = 4.0 * cost_ampAA;
        cost_5amp = 5.0 * cost_ampAA;

        if ((cost_4amp+cost_df) <= memory_mb) {
             outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_4amp);
             outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_4amp+cost_df);
             nincore_amp = 4;
             t2_incore = true;
             df_ints_incore = true;
        }
        else if ((cost_3amp+cost_df) <= memory_mb) {
             outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_3amp);
             //outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_3amp+cost_df);
             outfile->Printf("\tWarning: T2 amplitudes will be stored on the disk!\n");
             nincore_amp = 3;
             t2_incore = false;
             df_ints_incore = false;
        }
        else if (cost_3amp < memory_mb && cost_df < memory_mb ) {
             outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_3amp);
             outfile->Printf("\tWarning: T2 amplitudes will be stored on the disk!\n");
             nincore_amp = 3;
             t2_incore = false;
             df_ints_incore = false;
        }
        else {
             outfile->Printf("\tWarning: There is NOT enough memory for CC contractions!\n");
             outfile->Printf("\tIncrease memory by                    : %9.2lf MB \n", cost_3amp+cost_df-memory_mb);
             throw PSIEXCEPTION("There is NOT enough memory for CC contractions!");
        }

        // W_abef term
        cost_ampAA = 0.0;
        cost_ampAA = naoccA * naoccA * navirA * navirA;
        cost_ampAA += 2.0 * nQ * navirA * navirA;
        cost_ampAA += navirA * navirA * navirA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        double cost_ampAA2 = 0.0;
        cost_ampAA2 = naoccA * naoccA * navirA * navirA;
        cost_ampAA2 += nQ * navirA * navirA;
        cost_ampAA2 += 3.0 * navirA * navirA * navirA;
        cost_ampAA2 /= 1024.0 * 1024.0;
        cost_ampAA2 *= sizeof(double);
        cost_amp = MAX0(cost_ampAA, cost_ampAA2);
        outfile->Printf("\tMemory requirement for Wabef term     : %9.2lf MB \n", cost_amp);

        // QCHF
        if (qchf_ == "TRUE") qchf();

        // Fock
        fock();

        // ROHF REF
        if (reference == "ROHF") {
            t1A = SharedTensor2d(new Tensor2d("T1_1 <I|A>", naoccA, navirA));
            t1B = SharedTensor2d(new Tensor2d("T1_1 <i|a>", naoccB, navirB));
	    t1_1st_sc();
	}
	mp3_t2_1st_sc();
	Emp2L=Emp2;
        EcorrL=Emp2L-Escf;
	Emp2L_old=Emp2;

	outfile->Printf("\n");
	if (reference == "ROHF") outfile->Printf("\tComputing DF-MP2 energy using SCF MOs (DF-ROHF-MP2)... \n");
	else outfile->Printf("\tComputing DF-MP2 energy using SCF MOs (Canonical DF-MP2)... \n");
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	outfile->Printf("\tDF-HF Energy (a.u.)                : %20.14f\n", Escf);
	outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA+Escsmp2BB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SCS-MP2 Total Energy (a.u.)     : %20.14f\n", Escsmp2);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SOS-MP2 Total Energy (a.u.)     : %20.14f\n", Esosmp2);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SCSN-MP2 Total Energy (a.u.)    : %20.14f\n", Escsnmp2);
	if (reference == "ROHF") outfile->Printf("\tDF-MP2 Singles Energy (a.u.)       : %20.14f\n", Emp2_t1);
	if (reference == "ROHF") outfile->Printf("\tDF-MP2 Doubles Energy (a.u.)       : %20.14f\n", Ecorr - Emp2_t1);
	outfile->Printf("\tDF-MP2 Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
	outfile->Printf("\tDF-MP2 Total Energy (a.u.)         : %20.14f\n", Emp2);
	outfile->Printf("\t======================================================================= \n");

	Process::environment.globals["MP2 TOTAL ENERGY"] = Emp2;
	Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = Escsmp2;
	Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
	Process::environment.globals["SCSN-MP2 TOTAL ENERGY"] = Escsnmp2;
        Process::environment.globals["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
        Process::environment.globals["SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
        Process::environment.globals["SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
        Process::environment.globals["SCSN-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;
        Process::environment.globals["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
        Process::environment.globals["MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA+Emp2BB;

        // Perform MP3 iterations
        timer_on("MP3");
	t2_2nd_sc();
        timer_off("MP3");

	outfile->Printf("\n");
	outfile->Printf("\tComputing DF-MP3 energy using SCF MOs (Canonical DF-MP3)... \n");
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	outfile->Printf("\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp3AA);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp3AB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp3BB);
	outfile->Printf("\t3rd Order Energy (a.u.)            : %20.14f\n", Emp3-Emp2);
	outfile->Printf("\tDF-MP2.5 Correlation Energy (a.u.) : %20.14f\n", (Emp2 - Escf) + 0.5 * (Emp3-Emp2));
	outfile->Printf("\tDF-MP2.5 Total Energy (a.u.)       : %20.14f\n", 0.5 * (Emp3+Emp2));
	outfile->Printf("\tDF-MP3 Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
	outfile->Printf("\tDF-MP3 Total Energy (a.u.)         : %20.14f\n", Emp3);
	outfile->Printf("\t======================================================================= \n");

	Process::environment.globals["MP3 TOTAL ENERGY"] = Emp3;
        Process::environment.globals["MP3 CORRELATION ENERGY"] = Emp3 - Escf;
	Emp3L=Emp3;
        EcorrL=Emp3L-Escf;
	Emp3L_old=Emp3;

        // Malloc for PDMs
	gQt = SharedTensor1d(new Tensor1d("CCD PDM G_Qt", nQ));
	if (reference_ == "RESTRICTED") {
	    G1c_ov = SharedTensor2d(new Tensor2d("Correlation OPDM <O|V>", noccA, nvirA));
	    G1c_vo = SharedTensor2d(new Tensor2d("Correlation OPDM <V|O>", nvirA, noccA));
	}
	else if (reference_ == "UNRESTRICTED") {
	    G1c_ovA = SharedTensor2d(new Tensor2d("Correlation OPDM <O|V>", noccA, nvirA));
	    G1c_ovB = SharedTensor2d(new Tensor2d("Correlation OPDM <o|v>", noccB, nvirB));
	    G1c_voA = SharedTensor2d(new Tensor2d("Correlation OPDM <V|O>", nvirA, noccA));
	    G1c_voB = SharedTensor2d(new Tensor2d("Correlation OPDM <v|o>", nvirB, noccB));
	}

	mp3_pdm_3index_intr();
	omp3_opdm();
	omp3_tpdm();
	//ccl_energy();
	sep_tpdm_cc();
	gfock_cc_vo();
	gfock_cc_ov();
        gfock_cc_oo();
        gfock_cc_vv();
	idp();
	mograd();
        occ_iterations();

        // main if
        if (rms_wog <= tol_grad && fabs(DE) >= tol_Eod) {
           orbs_already_opt = 1;
	   if (conver == 1) outfile->Printf("\n\tOrbitals are optimized now.\n");
	   else if (conver == 0) {
                    outfile->Printf("\n\tMAX MOGRAD did NOT converged, but RMS MOGRAD converged!!!\n");
	            outfile->Printf("\tI will consider the present orbitals as optimized.\n");
           }
	   outfile->Printf("\tTransforming MOs to the semicanonical basis... \n");
	   semi_canonic();
	   outfile->Printf("\tSwitching to the standard DF-MP3 computation... \n");
           trans_corr();
           trans_ref();
           fock();
           ref_energy();
	   mp3_t2_1st_sc();
	   t2_2nd_sc();
           conver = 1;
           if (dertype == "FIRST") {
	       mp3_pdm_3index_intr();
	       omp3_opdm();
	       omp3_tpdm();
	       sep_tpdm_cc();
	       gfock_cc_vo();
	       gfock_cc_ov();
               gfock_cc_oo();
               gfock_cc_vv();
           }
        }// end main if

        else if (rms_wog <= tol_grad && fabs(DE) >= tol_Eod && regularization == "TRUE") {
	   outfile->Printf("\tOrbital gradient converged, but energy did not... \n");
	   outfile->Printf("\tA tighter rms_mograd_convergence tolerance is recommended... \n");
           throw PSIEXCEPTION("A tighter rms_mograd_convergence tolerance is recommended.");
        }


  if (conver == 1) {
        if (orbs_already_opt == 1) Emp3L = Emp3;

	outfile->Printf("\n");
	outfile->Printf("\tComputing DF-MP3 energy using optimized MOs... \n");
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	outfile->Printf("\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp3AA);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp3AB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp3BB);
	outfile->Printf("\t3rd Order Energy (a.u.)            : %20.14f\n", Emp3-Emp2);
	outfile->Printf("\tDF-MP2.5 Correlation Energy (a.u.) : %20.14f\n", (Emp2 - Escf) + 0.5 * (Emp3-Emp2));
	outfile->Printf("\tDF-MP2.5 Total Energy (a.u.)       : %20.14f\n", 0.5 * (Emp3+Emp2));
	outfile->Printf("\tDF-MP3 Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
	outfile->Printf("\tDF-MP3 Total Energy (a.u.)         : %20.14f\n", Emp3);
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\n");

	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\t================ DF-OMP3 FINAL RESULTS ================================ \n");
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	outfile->Printf("\tDF-HF Energy (a.u.)                : %20.14f\n", Escf);
	outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	outfile->Printf("\tDF-OMP3 Correlation Energy (a.u.)  : %20.14f\n", Emp3L-Escf);
	outfile->Printf("\tEdfomp3 - Eref (a.u.)              : %20.14f\n", Emp3L-Eref);
	outfile->Printf("\tDF-OMP3 Total Energy (a.u.)        : %20.14f\n", Emp3L);
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\n");

	// Set the global variables with the energies
	//Emp3L=Emp3;
	Process::environment.globals["CURRENT ENERGY"] = Emp3L;
        Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
        Process::environment.globals["CURRENT CORRELATION ENERGY"] = Emp3L - Escf;
	Process::environment.globals["OMP3 TOTAL ENERGY"] = Emp3L;
        Process::environment.globals["OMP3 CORRELATION ENERGY"] = Emp3L - Escf;

        // OEPROP
        if (oeprop_ == "TRUE") oeprop();

        // Compute Analytic Gradients
        if (dertype == "FIRST") dfgrad();

	// Save MOs to wfn
	save_mo_to_wfn();

  }// end if (conver == 1)


}// end omp3_manager

//======================================================================
//             MP3 Manager
//======================================================================
void DFOCC::mp3_manager()
{

        time4grad = 0;// means i will not compute the gradient
	mo_optimized = 0;// means MOs are not optimized

        timer_on("DF CC Integrals");
        df_corr();
        trans_corr();
        timer_off("DF CC Integrals");
        outfile->Printf("\n\tNumber of basis functions in the DF-CC basis: %3d\n", nQ);

        if (dertype == "FIRST" || oeprop_ == "TRUE" || ekt_ip_ == "TRUE" || qchf_ == "TRUE") {
            timer_on("DF REF Integrals");
            df_ref();
            trans_ref();
            timer_off("DF REF Integrals");
            outfile->Printf("\tNumber of basis functions in the DF-HF basis: %3d\n", nQ_ref);
            Jc = SharedTensor1d(new Tensor1d("DF_BASIS_SCF J_Q", nQ_ref));
        }

        // avaliable mem
        memory = Process::environment.get_memory();
        memory_mb = (double)memory/(1024.0 * 1024.0);
        outfile->Printf("\n\tAvailable memory                      : %9.2lf MB \n", memory_mb);

        // memory requirements
        // DF-CC B(Q,ab) + B(Q,ia) + B(Q,ij)
        cost_df = 0.0;
        cost_df = (navirA * navirA) + (navirA * naoccA) + (naoccA * naoccA);
        cost_df *= nQ;
        cost_df /= 1024.0 * 1024.0;
        cost_df *= sizeof(double);
        if (reference_ == "RESTRICTED") outfile->Printf("\tMemory requirement for 3-index ints   : %9.2lf MB \n", cost_df);
	else if (reference_ == "UNRESTRICTED") outfile->Printf("\tMemory requirement for 3-index ints   : %9.2lf MB \n", 2.0*cost_df);

        // Cost of Integral transform for B(Q,ab)
        cost_ampAA = 0.0;
        cost_ampAA = nQ * nso2_;
        cost_ampAA += nQ * navirA * navirA;
        cost_ampAA += nQ * nso_ * navirA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        outfile->Printf("\tMemory requirement for DF-CC int trans: %9.2lf MB \n", cost_ampAA);

        // Mem for amplitudes
        cost_ampAA = 0.0;
        cost_ampAA = nocc2AA * nvir2AA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        cost_3amp = 3.0 * cost_ampAA;
        cost_4amp = 4.0 * cost_ampAA;
        cost_5amp = 5.0 * cost_ampAA;

        if ((cost_4amp+cost_df) <= memory_mb) {
             outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_4amp);
             outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_4amp+cost_df);
             nincore_amp = 4;
             t2_incore = true;
             df_ints_incore = true;
        }
        else if ((cost_3amp+cost_df) <= memory_mb) {
             outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_3amp);
             //outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_3amp+cost_df);
             outfile->Printf("\tWarning: T2 amplitudes will be stored on the disk!\n");
             nincore_amp = 3;
             t2_incore = false;
             df_ints_incore = false;
        }
        else if (cost_3amp < memory_mb && cost_df < memory_mb ) {
             outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_3amp);
             outfile->Printf("\tWarning: T2 amplitudes will be stored on the disk!\n");
             nincore_amp = 3;
             t2_incore = false;
             df_ints_incore = false;
        }
        else {
             outfile->Printf("\tWarning: There is NOT enough memory for CC contractions!\n");
             outfile->Printf("\tIncrease memory by                    : %9.2lf MB \n", cost_3amp+cost_df-memory_mb);
             throw PSIEXCEPTION("There is NOT enough memory for CC contractions!");
        }

        // W_abef term
        cost_ampAA = 0.0;
        cost_ampAA = naoccA * naoccA * navirA * navirA;
        cost_ampAA += 2.0 * nQ * navirA * navirA;
        cost_ampAA += navirA * navirA * navirA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        double cost_ampAA2 = 0.0;
        cost_ampAA2 = naoccA * naoccA * navirA * navirA;
        cost_ampAA2 += nQ * navirA * navirA;
        cost_ampAA2 += 3.0 * navirA * navirA * navirA;
        cost_ampAA2 /= 1024.0 * 1024.0;
        cost_ampAA2 *= sizeof(double);
        cost_amp = MAX0(cost_ampAA, cost_ampAA2);
        outfile->Printf("\tMemory requirement for Wabef term     : %9.2lf MB \n", cost_amp);

        // Mem alloc for DF ints
	/*
        if (df_ints_incore) {
            bQijA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA));
            bQiaA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA));
            bQabA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA));
            bQijA->read(psio_, PSIF_DFOCC_INTS);
            bQiaA->read(psio_, PSIF_DFOCC_INTS);
            bQabA->read(psio_, PSIF_DFOCC_INTS, true, true);
        }
	*/

	/*
        if (t2_incore) {
            t2 = SharedTensor2d(new Tensor2d("T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
        }
	*/

        // memalloc for density intermediates
        if (qchf_ == "TRUE" || dertype == "FIRST") {
            g1Qc = SharedTensor1d(new Tensor1d("DF_BASIS_SCF G1_Q", nQ_ref));
            g1Qt = SharedTensor1d(new Tensor1d("DF_BASIS_SCF G1t_Q", nQ_ref));
            g1Qp = SharedTensor1d(new Tensor1d("DF_BASIS_SCF G1p_Q", nQ_ref));
            g1Q = SharedTensor1d(new Tensor1d("DF_BASIS_CC G1_Q", nQ));
            g1Qt2 = SharedTensor1d(new Tensor1d("DF_BASIS_CC G1t_Q", nQ));
        }

        // QCHF
        if (qchf_ == "TRUE") qchf();

        // Compute MP2 energy
        if (reference == "ROHF") {
            t1A = SharedTensor2d(new Tensor2d("T1_1 <I|A>", naoccA, navirA));
            t1B = SharedTensor2d(new Tensor2d("T1_1 <i|a>", naoccB, navirB));
	    t1_1st_sc();
	}
	mp3_t2_1st_sc();

	outfile->Printf("\n");
	if (reference == "ROHF") outfile->Printf("\tComputing DF-MP2 energy (DF-ROHF-MP2)... \n");
	else outfile->Printf("\tComputing DF-MP2 energy ... \n");
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	outfile->Printf("\tDF-HF Energy (a.u.)                : %20.14f\n", Escf);
	outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA+Escsmp2BB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SCS-MP2 Total Energy (a.u.)     : %20.14f\n", Escsmp2);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SOS-MP2 Total Energy (a.u.)     : %20.14f\n", Esosmp2);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SCSN-MP2 Total Energy (a.u.)    : %20.14f\n", Escsnmp2);
	if (reference == "ROHF") outfile->Printf("\tDF-MP2 Singles Energy (a.u.)       : %20.14f\n", Emp2_t1);
	if (reference == "ROHF") outfile->Printf("\tDF-MP2 Doubles Energy (a.u.)       : %20.14f\n", Ecorr - Emp2_t1);
	outfile->Printf("\tDF-MP2 Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
	outfile->Printf("\tDF-MP2 Total Energy (a.u.)         : %20.14f\n", Emp2);
	outfile->Printf("\t======================================================================= \n");

	Process::environment.globals["MP2 TOTAL ENERGY"] = Emp2;
	Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = Escsmp2;
	Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
	Process::environment.globals["SCSN-MP2 TOTAL ENERGY"] = Escsnmp2;
        Process::environment.globals["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
        Process::environment.globals["SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
        Process::environment.globals["SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
        Process::environment.globals["SCSN-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;
        Process::environment.globals["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
        Process::environment.globals["MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA+Emp2BB;

        // Perform MP3 iterations
        timer_on("MP3");
	t2_2nd_sc();
        timer_off("MP3");

	outfile->Printf("\n");
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\t================ MP3 FINAL RESULTS ==================================== \n");
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	outfile->Printf("\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp3AA);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp3AB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp3BB);
	outfile->Printf("\t3rd Order Energy (a.u.)            : %20.14f\n", Emp3-Emp2);
	outfile->Printf("\tDF-MP2.5 Correlation Energy (a.u.) : %20.14f\n", (Emp2 - Escf) + 0.5 * (Emp3-Emp2));
	outfile->Printf("\tDF-MP2.5 Total Energy (a.u.)       : %20.14f\n", 0.5 * (Emp3+Emp2));
	outfile->Printf("\tDF-MP3 Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
	outfile->Printf("\tDF-MP3 Total Energy (a.u.)         : %20.14f\n", Emp3);
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\n");

	Process::environment.globals["CURRENT ENERGY"] = Emp3;
        Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
        Process::environment.globals["CURRENT CORRELATION ENERGY"] = Emp3 - Escf;
	Process::environment.globals["MP3 TOTAL ENERGY"] = Emp3;
        Process::environment.globals["MP3 CORRELATION ENERGY"] = Emp3 - Escf;
	Emp3L=Emp3;

        // Compute Analytic Gradients
        if (dertype == "FIRST" || ekt_ip_ == "TRUE") {
	    // memalloc
	    gQt = SharedTensor1d(new Tensor1d("CCD PDM G_Qt", nQ));
	    if (reference_ == "RESTRICTED") {
	        G1c_ov = SharedTensor2d(new Tensor2d("Correlation OPDM <O|V>", noccA, nvirA));
	        G1c_vo = SharedTensor2d(new Tensor2d("Correlation OPDM <V|O>", nvirA, noccA));
	    }
	    else if (reference_ == "UNRESTRICTED") {
	        G1c_ovA = SharedTensor2d(new Tensor2d("Correlation OPDM <O|V>", noccA, nvirA));
	        G1c_ovB = SharedTensor2d(new Tensor2d("Correlation OPDM <o|v>", noccB, nvirB));
	        G1c_voA = SharedTensor2d(new Tensor2d("Correlation OPDM <V|O>", nvirA, noccA));
	        G1c_voB = SharedTensor2d(new Tensor2d("Correlation OPDM <v|o>", nvirB, noccB));
	    }

            outfile->Printf("\tComputing unrelaxed response density matrices...\n");
	    mp3_pdm_3index_intr();
 	    omp3_opdm();
	    omp3_tpdm();
	    //ccl_energy();
            prepare4grad();
            if (oeprop_ == "TRUE") oeprop();
            if (dertype == "FIRST") dfgrad();
            //if (ekt_ip_ == "TRUE") ekt_ip();
        }// if (dertype == "FIRST" || ekt_ip_ == "TRUE")

}// end mp3_manager

//======================================================================
//             OMP2.5 Manager
//======================================================================
void DFOCC::omp2_5_manager()
{
        time4grad = 0;// means I will not compute the gradient
	mo_optimized = 0;// means MOs are not optimized
	orbs_already_opt = 0;// means orbitals are not optimized yet.
	orbs_already_sc = 0;// means orbitals are not semicanonical yet.

        timer_on("DF CC Integrals");
        df_corr();
        trans_corr();
        df_ref();
        trans_ref();
        outfile->Printf("\tNumber of basis functions in the DF-HF basis: %3d\n", nQ_ref);
        outfile->Printf("\tNumber of basis functions in the DF-CC basis: %3d\n", nQ);
        timer_off("DF CC Integrals");

        // memalloc for density intermediates
        Jc = SharedTensor1d(new Tensor1d("DF_BASIS_SCF J_Q", nQ_ref));
        g1Qc = SharedTensor1d(new Tensor1d("DF_BASIS_SCF G1_Q", nQ_ref));
        g1Qt = SharedTensor1d(new Tensor1d("DF_BASIS_SCF G1t_Q", nQ_ref));
        g1Qp = SharedTensor1d(new Tensor1d("DF_BASIS_SCF G1p_Q", nQ_ref));
        g1Q = SharedTensor1d(new Tensor1d("DF_BASIS_CC G1_Q", nQ));
        g1Qt2 = SharedTensor1d(new Tensor1d("DF_BASIS_CC G1t_Q", nQ));

        // avaliable mem
        memory = Process::environment.get_memory();
        memory_mb = (double)memory/(1024.0 * 1024.0);
        outfile->Printf("\n\tAvailable memory                      : %9.2lf MB \n", memory_mb);

        // memory requirements
        // DF-CC B(Q,ab) + B(Q,ia) + B(Q,ij)
        cost_df = 0.0;
        cost_df = (navirA * navirA) + (navirA * naoccA) + (naoccA * naoccA);
        cost_df *= nQ;
        cost_df /= 1024.0 * 1024.0;
        cost_df *= sizeof(double);
        if (reference_ == "RESTRICTED") outfile->Printf("\tMemory requirement for 3-index ints   : %9.2lf MB \n", cost_df);
	else if (reference_ == "UNRESTRICTED") outfile->Printf("\tMemory requirement for 3-index ints   : %9.2lf MB \n", 2.0*cost_df);

        // Cost of Integral transform for B(Q,ab)
        cost_ampAA = 0.0;
        cost_ampAA = nQ * nso2_;
        cost_ampAA += nQ * navirA * navirA;
        cost_ampAA += nQ * nso_ * navirA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        outfile->Printf("\tMemory requirement for DF-CC int trans: %9.2lf MB \n", cost_ampAA);

        // Mem for amplitudes
        cost_ampAA = 0.0;
        cost_ampAA = nocc2AA * nvir2AA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        cost_3amp = 3.0 * cost_ampAA;
        cost_4amp = 4.0 * cost_ampAA;
        cost_5amp = 5.0 * cost_ampAA;

        if ((cost_4amp+cost_df) <= memory_mb) {
             outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_4amp);
             outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_4amp+cost_df);
             nincore_amp = 4;
             t2_incore = true;
             df_ints_incore = true;
        }
        else if ((cost_3amp+cost_df) <= memory_mb) {
             outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_3amp);
             //outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_3amp+cost_df);
             outfile->Printf("\tWarning: T2 amplitudes will be stored on the disk!\n");
             nincore_amp = 3;
             t2_incore = false;
             df_ints_incore = false;
        }
        else if (cost_3amp < memory_mb && cost_df < memory_mb ) {
             outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_3amp);
             outfile->Printf("\tWarning: T2 amplitudes will be stored on the disk!\n");
             nincore_amp = 3;
             t2_incore = false;
             df_ints_incore = false;
        }
        else {
             outfile->Printf("\tWarning: There is NOT enough memory for CC contractions!\n");
             outfile->Printf("\tIncrease memory by                    : %9.2lf MB \n", cost_3amp+cost_df-memory_mb);
             throw PSIEXCEPTION("There is NOT enough memory for CC contractions!");
        }

        // W_abef term
        cost_ampAA = 0.0;
        cost_ampAA = naoccA * naoccA * navirA * navirA;
        cost_ampAA += 2.0 * nQ * navirA * navirA;
        cost_ampAA += navirA * navirA * navirA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        double cost_ampAA2 = 0.0;
        cost_ampAA2 = naoccA * naoccA * navirA * navirA;
        cost_ampAA2 += nQ * navirA * navirA;
        cost_ampAA2 += 3.0 * navirA * navirA * navirA;
        cost_ampAA2 /= 1024.0 * 1024.0;
        cost_ampAA2 *= sizeof(double);
        cost_amp = MAX0(cost_ampAA, cost_ampAA2);
        outfile->Printf("\tMemory requirement for Wabef term     : %9.2lf MB \n", cost_amp);

        // Fock
        fock();

        // QCHF
        if (qchf_ == "TRUE") qchf();

        // ROHF REF
        if (reference == "ROHF") {
            t1A = SharedTensor2d(new Tensor2d("T1_1 <I|A>", naoccA, navirA));
            t1B = SharedTensor2d(new Tensor2d("T1_1 <i|a>", naoccB, navirB));
	    t1_1st_sc();
	}
	mp3_t2_1st_sc();
	Emp2L=Emp2;
        EcorrL=Emp2L-Escf;
	Emp2L_old=Emp2;

	outfile->Printf("\n");
	if (reference == "ROHF") outfile->Printf("\tComputing DF-MP2 energy using SCF MOs (DF-ROHF-MP2)... \n");
	else outfile->Printf("\tComputing DF-MP2 energy using SCF MOs (Canonical DF-MP2)... \n");
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	outfile->Printf("\tDF-HF Energy (a.u.)                : %20.14f\n", Escf);
	outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA+Escsmp2BB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SCS-MP2 Total Energy (a.u.)     : %20.14f\n", Escsmp2);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SOS-MP2 Total Energy (a.u.)     : %20.14f\n", Esosmp2);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SCSN-MP2 Total Energy (a.u.)    : %20.14f\n", Escsnmp2);
	if (reference == "ROHF") outfile->Printf("\tDF-MP2 Singles Energy (a.u.)       : %20.14f\n", Emp2_t1);
	if (reference == "ROHF") outfile->Printf("\tDF-MP2 Doubles Energy (a.u.)       : %20.14f\n", Ecorr - Emp2_t1);
	outfile->Printf("\tDF-MP2 Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
	outfile->Printf("\tDF-MP2 Total Energy (a.u.)         : %20.14f\n", Emp2);
	outfile->Printf("\t======================================================================= \n");

	Process::environment.globals["MP2 TOTAL ENERGY"] = Emp2;
	Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = Escsmp2;
	Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
	Process::environment.globals["SCSN-MP2 TOTAL ENERGY"] = Escsnmp2;
        Process::environment.globals["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
        Process::environment.globals["SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
        Process::environment.globals["SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
        Process::environment.globals["SCSN-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;
        Process::environment.globals["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
        Process::environment.globals["MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA+Emp2BB;

        // Perform MP3 iterations
        timer_on("MP3");
	t2_2nd_sc();
        timer_off("MP3");

	outfile->Printf("\n");
	outfile->Printf("\tComputing DF-MP2.5 energy using SCF MOs (Canonical DF-MP2.5)... \n");
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	outfile->Printf("\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp3AA);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp3AB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp3BB);
	outfile->Printf("\tDF-MP3 Correlation Energy (a.u.)   : %20.14f\n", (Emp2-Escf) + 2.0*(Emp3-Emp2));
	outfile->Printf("\tDF-MP3 Total Energy (a.u.)         : %20.14f\n", Emp2 + 2.0*(Emp3-Emp2));
	outfile->Printf("\tDF-MP2.5 Correlation Energy (a.u.) : %20.14f\n", Ecorr);
	outfile->Printf("\tDF-MP2.5 Total Energy (a.u.)       : %20.14f\n", Emp3);
	outfile->Printf("\t======================================================================= \n");

	Process::environment.globals["MP2.5 TOTAL ENERGY"] = Emp3;
        Process::environment.globals["MP2.5 CORRELATION ENERGY"] = Emp3 - Escf;
	Emp3L=Emp3;
        EcorrL=Emp3L-Escf;
	Emp3L_old=Emp3;

        // Malloc for PDMs
	gQt = SharedTensor1d(new Tensor1d("CCD PDM G_Qt", nQ));
	if (reference_ == "RESTRICTED") {
	    G1c_ov = SharedTensor2d(new Tensor2d("Correlation OPDM <O|V>", noccA, nvirA));
	    G1c_vo = SharedTensor2d(new Tensor2d("Correlation OPDM <V|O>", nvirA, noccA));
	}
	else if (reference_ == "UNRESTRICTED") {
	    G1c_ovA = SharedTensor2d(new Tensor2d("Correlation OPDM <O|V>", noccA, nvirA));
	    G1c_ovB = SharedTensor2d(new Tensor2d("Correlation OPDM <o|v>", noccB, nvirB));
	    G1c_voA = SharedTensor2d(new Tensor2d("Correlation OPDM <V|O>", nvirA, noccA));
	    G1c_voB = SharedTensor2d(new Tensor2d("Correlation OPDM <v|o>", nvirB, noccB));
	}

	mp3_pdm_3index_intr();
	omp3_opdm();
	omp3_tpdm();
	//ccl_energy();
	sep_tpdm_cc();
	gfock_cc_vo();
	gfock_cc_ov();
        gfock_cc_oo();
        gfock_cc_vv();
	idp();
	mograd();
        occ_iterations();

        // main if
        if (rms_wog <= tol_grad && fabs(DE) >= tol_Eod) {
           orbs_already_opt = 1;
	   if (conver == 1) outfile->Printf("\n\tOrbitals are optimized now.\n");
	   else if (conver == 0) {
                    outfile->Printf("\n\tMAX MOGRAD did NOT converged, but RMS MOGRAD converged!!!\n");
	            outfile->Printf("\tI will consider the present orbitals as optimized.\n");
           }
	   outfile->Printf("\tTransforming MOs to the semicanonical basis... \n");
	   semi_canonic();
	   outfile->Printf("\tSwitching to the standard DF-MP3 computation... \n");
           trans_corr();
           trans_ref();
           fock();
           ref_energy();
	   mp3_t2_1st_sc();
	   t2_2nd_sc();
           conver = 1;
           if (dertype == "FIRST") {
	       mp3_pdm_3index_intr();
	       omp3_opdm();
	       omp3_tpdm();
	       sep_tpdm_cc();
	       gfock_cc_vo();
	       gfock_cc_ov();
               gfock_cc_oo();
               gfock_cc_vv();
           }
        }// end main if

        else if (rms_wog <= tol_grad && fabs(DE) >= tol_Eod && regularization == "TRUE") {
	   outfile->Printf("\tOrbital gradient converged, but energy did not... \n");
	   outfile->Printf("\tA tighter rms_mograd_convergence tolerance is recommended... \n");
           throw PSIEXCEPTION("A tighter rms_mograd_convergence tolerance is recommended.");
        }


  if (conver == 1) {
        if (orbs_already_opt == 1) Emp3L = Emp3;

	outfile->Printf("\n");
	outfile->Printf("\tComputing DF-MP2.5 energy using optimized MOs... \n");
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	outfile->Printf("\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp3AA);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp3AB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp3BB);
	outfile->Printf("\tDF-MP3 Correlation Energy (a.u.)   : %20.14f\n", (Emp2-Escf) + 2.0*(Emp3-Emp2));
	outfile->Printf("\tDF-MP3 Total Energy (a.u.)         : %20.14f\n", Emp2 + 2.0*(Emp3-Emp2));
	outfile->Printf("\tDF-MP2.5 Correlation Energy (a.u.) : %20.14f\n", Ecorr);
	outfile->Printf("\tDF-MP2.5 Total Energy (a.u.)       : %20.14f\n", Emp3);
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\n");

	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\t================ DF-OMP2.5 FINAL RESULTS ============================== \n");
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	outfile->Printf("\tDF-HF Energy (a.u.)                : %20.14f\n", Escf);
	outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	outfile->Printf("\tDF-OMP2.5 Correlation Energy (a.u.): %20.14f\n", Emp3L-Escf);
	outfile->Printf("\tEdfomp2.5 - Eref (a.u.)            : %20.14f\n", Emp3L-Eref);
	outfile->Printf("\tDF-OMP2.5 Total Energy (a.u.)      : %20.14f\n", Emp3L);
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\n");

	// Set the global variables with the energies
	//Emp3L=Emp3;
	Process::environment.globals["CURRENT ENERGY"] = Emp3L;
        Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
        Process::environment.globals["CURRENT CORRELATION ENERGY"] = Emp3L - Escf;
	Process::environment.globals["OMP2.5 TOTAL ENERGY"] = Emp3L;
        Process::environment.globals["OMP2.5 CORRELATION ENERGY"] = Emp3L - Escf;

        // OEPROP
        if (oeprop_ == "TRUE") oeprop();

        // Compute Analytic Gradients
        if (dertype == "FIRST") dfgrad();

	// Save MOs to wfn
	save_mo_to_wfn();

  }// end if (conver == 1)


}// end omp2_5_manager

//======================================================================
//             MP2.5 Manager
//======================================================================
void DFOCC::mp2_5_manager()
{

        time4grad = 0;// means i will not compute the gradient
	mo_optimized = 0;// means MOs are not optimized

        timer_on("DF CC Integrals");
        df_corr();
        trans_corr();
        timer_off("DF CC Integrals");
        outfile->Printf("\n\tNumber of basis functions in the DF-CC basis: %3d\n", nQ);

        if (dertype == "FIRST" || oeprop_ == "TRUE" || ekt_ip_ == "TRUE" || qchf_ == "TRUE") {
            timer_on("DF REF Integrals");
            df_ref();
            trans_ref();
            timer_off("DF REF Integrals");
            outfile->Printf("\tNumber of basis functions in the DF-HF basis: %3d\n", nQ_ref);
            Jc = SharedTensor1d(new Tensor1d("DF_BASIS_SCF J_Q", nQ_ref));
        }

        // avaliable mem
        memory = Process::environment.get_memory();
        memory_mb = (double)memory/(1024.0 * 1024.0);
        outfile->Printf("\n\tAvailable memory                      : %9.2lf MB \n", memory_mb);

        // memory requirements
        // DF-CC B(Q,ab) + B(Q,ia) + B(Q,ij)
        cost_df = 0.0;
        cost_df = (navirA * navirA) + (navirA * naoccA) + (naoccA * naoccA);
        cost_df *= nQ;
        cost_df /= 1024.0 * 1024.0;
        cost_df *= sizeof(double);
        if (reference_ == "RESTRICTED") outfile->Printf("\tMemory requirement for 3-index ints   : %9.2lf MB \n", cost_df);
	else if (reference_ == "UNRESTRICTED") outfile->Printf("\tMemory requirement for 3-index ints   : %9.2lf MB \n", 2.0*cost_df);

        // Cost of Integral transform for B(Q,ab)
        cost_ampAA = 0.0;
        cost_ampAA = nQ * nso2_;
        cost_ampAA += nQ * navirA * navirA;
        cost_ampAA += nQ * nso_ * navirA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        outfile->Printf("\tMemory requirement for DF-CC int trans: %9.2lf MB \n", cost_ampAA);

        // Mem for amplitudes
        cost_ampAA = 0.0;
        cost_ampAA = nocc2AA * nvir2AA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        cost_3amp = 3.0 * cost_ampAA;
        cost_4amp = 4.0 * cost_ampAA;
        cost_5amp = 5.0 * cost_ampAA;

        if ((cost_4amp+cost_df) <= memory_mb) {
             outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_4amp);
             outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_4amp+cost_df);
             nincore_amp = 4;
             t2_incore = true;
             df_ints_incore = true;
        }
        else if ((cost_3amp+cost_df) <= memory_mb) {
             outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_3amp);
             //outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_3amp+cost_df);
             outfile->Printf("\tWarning: T2 amplitudes will be stored on the disk!\n");
             nincore_amp = 3;
             t2_incore = false;
             df_ints_incore = false;
        }
        else if (cost_3amp < memory_mb && cost_df < memory_mb ) {
             outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_3amp);
             outfile->Printf("\tWarning: T2 amplitudes will be stored on the disk!\n");
             nincore_amp = 3;
             t2_incore = false;
             df_ints_incore = false;
        }
        else {
             outfile->Printf("\tWarning: There is NOT enough memory for CC contractions!\n");
             outfile->Printf("\tIncrease memory by                    : %9.2lf MB \n", cost_3amp+cost_df-memory_mb);
             throw PSIEXCEPTION("There is NOT enough memory for CC contractions!");
        }

        // W_abef term
        cost_ampAA = 0.0;
        cost_ampAA = naoccA * naoccA * navirA * navirA;
        cost_ampAA += 2.0 * nQ * navirA * navirA;
        cost_ampAA += navirA * navirA * navirA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        double cost_ampAA2 = 0.0;
        cost_ampAA2 = naoccA * naoccA * navirA * navirA;
        cost_ampAA2 += nQ * navirA * navirA;
        cost_ampAA2 += 3.0 * navirA * navirA * navirA;
        cost_ampAA2 /= 1024.0 * 1024.0;
        cost_ampAA2 *= sizeof(double);
        cost_amp = MAX0(cost_ampAA, cost_ampAA2);
        outfile->Printf("\tMemory requirement for Wabef term     : %9.2lf MB \n", cost_amp);

        // memalloc for density intermediates
        if (qchf_ == "TRUE" || dertype == "FIRST") {
            g1Qc = SharedTensor1d(new Tensor1d("DF_BASIS_SCF G1_Q", nQ_ref));
            g1Qt = SharedTensor1d(new Tensor1d("DF_BASIS_SCF G1t_Q", nQ_ref));
            g1Qp = SharedTensor1d(new Tensor1d("DF_BASIS_SCF G1p_Q", nQ_ref));
            g1Q = SharedTensor1d(new Tensor1d("DF_BASIS_CC G1_Q", nQ));
            g1Qt2 = SharedTensor1d(new Tensor1d("DF_BASIS_CC G1t_Q", nQ));
        }

        // QCHF
        if (qchf_ == "TRUE") qchf();

        // Compute MP2 energy
        if (reference == "ROHF") {
            t1A = SharedTensor2d(new Tensor2d("T1_1 <I|A>", naoccA, navirA));
            t1B = SharedTensor2d(new Tensor2d("T1_1 <i|a>", naoccB, navirB));
	    t1_1st_sc();
	}
	mp3_t2_1st_sc();

	outfile->Printf("\n");
	if (reference == "ROHF") outfile->Printf("\tComputing DF-MP2 energy (DF-ROHF-MP2)... \n");
	else outfile->Printf("\tComputing DF-MP2 energy ... \n");
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	outfile->Printf("\tDF-HF Energy (a.u.)                : %20.14f\n", Escf);
	outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA+Escsmp2BB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SCS-MP2 Total Energy (a.u.)     : %20.14f\n", Escsmp2);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SOS-MP2 Total Energy (a.u.)     : %20.14f\n", Esosmp2);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SCSN-MP2 Total Energy (a.u.)    : %20.14f\n", Escsnmp2);
	if (reference == "ROHF") outfile->Printf("\tDF-MP2 Singles Energy (a.u.)       : %20.14f\n", Emp2_t1);
	if (reference == "ROHF") outfile->Printf("\tDF-MP2 Doubles Energy (a.u.)       : %20.14f\n", Ecorr - Emp2_t1);
	outfile->Printf("\tDF-MP2 Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
	outfile->Printf("\tDF-MP2 Total Energy (a.u.)         : %20.14f\n", Emp2);
	outfile->Printf("\t======================================================================= \n");

	Process::environment.globals["MP2 TOTAL ENERGY"] = Emp2;
	Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = Escsmp2;
	Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
	Process::environment.globals["SCSN-MP2 TOTAL ENERGY"] = Escsnmp2;
        Process::environment.globals["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
        Process::environment.globals["SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
        Process::environment.globals["SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
        Process::environment.globals["SCSN-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;
        Process::environment.globals["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
        Process::environment.globals["MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA+Emp2BB;

        // Perform MP3 iterations
        timer_on("MP3");
	t2_2nd_sc();
        timer_off("MP3");

	outfile->Printf("\n");
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\t================ MP2.5 FINAL RESULTS ================================== \n");
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	outfile->Printf("\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp3AA);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp3AB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp3BB);
	outfile->Printf("\tDF-MP3 Correlation Energy (a.u.)   : %20.14f\n", (Emp2-Escf) + 2.0*(Emp3-Emp2));
	outfile->Printf("\tDF-MP3 Total Energy (a.u.)         : %20.14f\n", Emp2 + 2.0*(Emp3-Emp2));
	outfile->Printf("\tDF-MP2.5 Correlation Energy (a.u.) : %20.14f\n", Ecorr);
	outfile->Printf("\tDF-MP2.5 Total Energy (a.u.)       : %20.14f\n", Emp3);
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\n");

	Process::environment.globals["CURRENT ENERGY"] = Emp3;
        Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
        Process::environment.globals["CURRENT CORRELATION ENERGY"] = Emp3 - Escf;
	Process::environment.globals["MP2.5 TOTAL ENERGY"] = Emp3;
        Process::environment.globals["MP2.5 CORRELATION ENERGY"] = Emp3 - Escf;
	Emp3L=Emp3;

        // Compute Analytic Gradients
        if (dertype == "FIRST" || ekt_ip_ == "TRUE") {
	    // memalloc
	    gQt = SharedTensor1d(new Tensor1d("CCD PDM G_Qt", nQ));
	    if (reference_ == "RESTRICTED") {
	        G1c_ov = SharedTensor2d(new Tensor2d("Correlation OPDM <O|V>", noccA, nvirA));
	        G1c_vo = SharedTensor2d(new Tensor2d("Correlation OPDM <V|O>", nvirA, noccA));
	    }
	    else if (reference_ == "UNRESTRICTED") {
	        G1c_ovA = SharedTensor2d(new Tensor2d("Correlation OPDM <O|V>", noccA, nvirA));
	        G1c_ovB = SharedTensor2d(new Tensor2d("Correlation OPDM <o|v>", noccB, nvirB));
	        G1c_voA = SharedTensor2d(new Tensor2d("Correlation OPDM <V|O>", nvirA, noccA));
	        G1c_voB = SharedTensor2d(new Tensor2d("Correlation OPDM <v|o>", nvirB, noccB));
	    }

            outfile->Printf("\tComputing unrelaxed response density matrices...\n");
	    mp3_pdm_3index_intr();
 	    omp3_opdm();
	    omp3_tpdm();
	    //ccl_energy();
            prepare4grad();
            if (oeprop_ == "TRUE") oeprop();
            if (dertype == "FIRST") dfgrad();
            //if (ekt_ip_ == "TRUE") ekt_ip();
        }// if (dertype == "FIRST" || ekt_ip_ == "TRUE")

}// end mp2_5_manager

//======================================================================
//             OLCCD Manager
//======================================================================
void DFOCC::olccd_manager()
{

        time4grad = 0;// means I will not compute the gradient
	mo_optimized = 0;// means MOs are not optimized
	orbs_already_opt = 0;// means orbitals are not optimized yet.
	orbs_already_sc = 0;// means orbitals are not semicanonical yet.

        timer_on("DF CC Integrals");
        df_corr();
        trans_corr();
        df_ref();
        trans_ref();
        outfile->Printf("\tNumber of basis functions in the DF-HF basis: %3d\n", nQ_ref);
        outfile->Printf("\tNumber of basis functions in the DF-CC basis: %3d\n", nQ);
        timer_off("DF CC Integrals");

        // memalloc for density intermediates
        Jc = SharedTensor1d(new Tensor1d("DF_BASIS_SCF J_Q", nQ_ref));
        g1Qc = SharedTensor1d(new Tensor1d("DF_BASIS_SCF G1_Q", nQ_ref));
        g1Qt = SharedTensor1d(new Tensor1d("DF_BASIS_SCF G1t_Q", nQ_ref));
        g1Qp = SharedTensor1d(new Tensor1d("DF_BASIS_SCF G1p_Q", nQ_ref));
        g1Q = SharedTensor1d(new Tensor1d("DF_BASIS_CC G1_Q", nQ));
        g1Qt2 = SharedTensor1d(new Tensor1d("DF_BASIS_CC G1t_Q", nQ));

        // avaliable mem
        memory = Process::environment.get_memory();
        memory_mb = (double)memory/(1024.0 * 1024.0);
        outfile->Printf("\n\tAvailable memory                      : %9.2lf MB \n", memory_mb);

        // memory requirements
        // DF-CC B(Q,ab) + B(Q,ia) + B(Q,ij)
        cost_df = 0.0;
        cost_df = (navirA * navirA) + (navirA * naoccA) + (naoccA * naoccA);
        cost_df *= nQ;
        cost_df /= 1024.0 * 1024.0;
        cost_df *= sizeof(double);
        if (reference_ == "RESTRICTED") outfile->Printf("\tMemory requirement for 3-index ints   : %9.2lf MB \n", cost_df);
	else if (reference_ == "UNRESTRICTED") outfile->Printf("\tMemory requirement for 3-index ints   : %9.2lf MB \n", 2.0*cost_df);

        // Cost of Integral transform for B(Q,ab)
        cost_ampAA = 0.0;
        cost_ampAA = nQ * nso2_;
        cost_ampAA += nQ * navirA * navirA;
        cost_ampAA += nQ * nso_ * navirA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        outfile->Printf("\tMemory requirement for DF-CC int trans: %9.2lf MB \n", cost_ampAA);

        // Mem for amplitudes
        cost_ampAA = 0.0;
        cost_ampAA = nocc2AA * nvir2AA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        cost_3amp = 3.0 * cost_ampAA;
        cost_4amp = 4.0 * cost_ampAA;
        cost_5amp = 5.0 * cost_ampAA;

        if ((cost_4amp+cost_df) <= memory_mb) {
             outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_4amp);
             outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_4amp+cost_df);
             nincore_amp = 4;
             t2_incore = true;
             df_ints_incore = true;
        }
        else if ((cost_3amp+cost_df) <= memory_mb) {
             outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_3amp);
             //outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_3amp+cost_df);
             outfile->Printf("\tWarning: T2 amplitudes will be stored on the disk!\n");
             nincore_amp = 3;
             t2_incore = false;
             df_ints_incore = false;
        }
        else if (cost_3amp < memory_mb && cost_df < memory_mb ) {
             outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_3amp);
             outfile->Printf("\tWarning: T2 amplitudes will be stored on the disk!\n");
             nincore_amp = 3;
             t2_incore = false;
             df_ints_incore = false;
        }
        else {
             outfile->Printf("\tWarning: There is NOT enough memory for CC contractions!\n");
             outfile->Printf("\tIncrease memory by                    : %9.2lf MB \n", cost_3amp+cost_df-memory_mb);
             throw PSIEXCEPTION("There is NOT enough memory for CC contractions!");
        }

        // W_abef term
        cost_ampAA = 0.0;
        cost_ampAA = naoccA * naoccA * navirA * navirA;
        cost_ampAA += 2.0 * nQ * navirA * navirA;
        cost_ampAA += navirA * navirA * navirA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        double cost_ampAA2 = 0.0;
        cost_ampAA2 = naoccA * naoccA * navirA * navirA;
        cost_ampAA2 += nQ * navirA * navirA;
        cost_ampAA2 += 3.0 * navirA * navirA * navirA;
        cost_ampAA2 /= 1024.0 * 1024.0;
        cost_ampAA2 *= sizeof(double);
        cost_amp = MAX0(cost_ampAA, cost_ampAA2);
        outfile->Printf("\tMemory requirement for Wabef term     : %9.2lf MB \n", cost_amp);

        // Fock
        fock();

        // QCHF
        if (qchf_ == "TRUE") qchf();

        // ROHF REF
        if (reference == "ROHF") {
            t1A = SharedTensor2d(new Tensor2d("T1_1 <I|A>", naoccA, navirA));
            t1B = SharedTensor2d(new Tensor2d("T1_1 <i|a>", naoccB, navirB));
	    t1_1st_sc();
	}
	lccd_t2_1st_sc();
	Elccd=Emp2;
	ElccdL=Emp2;
        EcorrL=Emp2-Escf;
	ElccdL_old=Emp2;

	outfile->Printf("\n");
	if (reference == "ROHF") outfile->Printf("\tComputing DF-MP2 energy using SCF MOs (DF-ROHF-MP2)... \n");
	else outfile->Printf("\tComputing DF-MP2 energy using SCF MOs (Canonical DF-MP2)... \n");
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	outfile->Printf("\tDF-HF Energy (a.u.)                : %20.14f\n", Escf);
	outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA+Escsmp2BB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SCS-MP2 Total Energy (a.u.)     : %20.14f\n", Escsmp2);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SOS-MP2 Total Energy (a.u.)     : %20.14f\n", Esosmp2);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SCSN-MP2 Total Energy (a.u.)    : %20.14f\n", Escsnmp2);
	if (reference == "ROHF") outfile->Printf("\tDF-MP2 Singles Energy (a.u.)       : %20.14f\n", Emp2_t1);
	if (reference == "ROHF") outfile->Printf("\tDF-MP2 Doubles Energy (a.u.)       : %20.14f\n", Ecorr - Emp2_t1);
	outfile->Printf("\tDF-MP2 Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
	outfile->Printf("\tDF-MP2 Total Energy (a.u.)         : %20.14f\n", Emp2);
	outfile->Printf("\t======================================================================= \n");

	Process::environment.globals["MP2 TOTAL ENERGY"] = Emp2;
	Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = Escsmp2;
	Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
	Process::environment.globals["SCSN-MP2 TOTAL ENERGY"] = Escsnmp2;
        Process::environment.globals["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
        Process::environment.globals["SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
        Process::environment.globals["SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
        Process::environment.globals["SCSN-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;
        Process::environment.globals["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
        Process::environment.globals["MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA+Emp2BB;

        // Malloc for PDMs
	gQt = SharedTensor1d(new Tensor1d("CCD PDM G_Qt", nQ));
	if (reference_ == "RESTRICTED") {
	    G1c_ov = SharedTensor2d(new Tensor2d("Correlation OPDM <O|V>", noccA, nvirA));
	    G1c_vo = SharedTensor2d(new Tensor2d("Correlation OPDM <V|O>", nvirA, noccA));
	}
	else if (reference_ == "UNRESTRICTED") {
	    G1c_ovA = SharedTensor2d(new Tensor2d("Correlation OPDM <O|V>", noccA, nvirA));
	    G1c_ovB = SharedTensor2d(new Tensor2d("Correlation OPDM <o|v>", noccB, nvirB));
	    G1c_voA = SharedTensor2d(new Tensor2d("Correlation OPDM <V|O>", nvirA, noccA));
	    G1c_voB = SharedTensor2d(new Tensor2d("Correlation OPDM <v|o>", nvirB, noccB));
	}

	lccd_pdm_3index_intr();
	omp3_opdm();
	olccd_tpdm();
	sep_tpdm_cc();
	gfock_cc_vo();
	gfock_cc_ov();
        gfock_cc_oo();
        gfock_cc_vv();
	idp();
	mograd();
        occ_iterations();
	Elccd=ElccdL;

        // main if
        if (rms_wog <= tol_grad && fabs(DE) >= tol_Eod) {
           orbs_already_opt = 1;
	   mo_optimized = 1;
	   if (conver == 1) outfile->Printf("\n\tOrbitals are optimized now.\n");
	   else if (conver == 0) {
                    outfile->Printf("\n\tMAX MOGRAD did NOT converged, but RMS MOGRAD converged!!!\n");
	            outfile->Printf("\tI will consider the present orbitals as optimized.\n");
           }
	   outfile->Printf("\tTransforming MOs to the semicanonical basis... \n");
	   semi_canonic();
	   outfile->Printf("\tSwitching to the standard DF-LCCD computation... \n");
           trans_corr();
           trans_ref();
           fock();
           ref_energy();
           lccd_iterations();
           conver = 1;
           if (dertype == "FIRST") {
	       lccd_pdm_3index_intr();
	       omp3_opdm();
	       olccd_tpdm();
	       sep_tpdm_cc();
	       gfock_cc_vo();
	       gfock_cc_ov();
               gfock_cc_oo();
               gfock_cc_vv();
           }
        }// end main if

        else if (rms_wog <= tol_grad && fabs(DE) >= tol_Eod && regularization == "TRUE") {
	   outfile->Printf("\tOrbital gradient converged, but energy did not... \n");
	   outfile->Printf("\tA tighter rms_mograd_convergence tolerance is recommended... \n");
           throw PSIEXCEPTION("A tighter rms_mograd_convergence tolerance is recommended.");
        }


  if (conver == 1) {
        if (orbs_already_opt == 1) ElccdL = Elccd;

	outfile->Printf("\n");
	outfile->Printf("\tComputing DF-LCCD energy using optimized MOs... \n");
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	outfile->Printf("\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", ElccdAA);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", ElccdAB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", ElccdBB);
	outfile->Printf("\tDF-LCCD Correlation Energy (a.u.)  : %20.14f\n", Ecorr);
	outfile->Printf("\tDF-LCCD Total Energy (a.u.)        : %20.14f\n", Elccd);
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\n");

	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\t================ DF-OLCCD FINAL RESULTS =============================== \n");
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	outfile->Printf("\tDF-HF Energy (a.u.)                : %20.14f\n", Escf);
	outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	outfile->Printf("\tDF-OLCCD Correlation Energy (a.u.) : %20.14f\n", ElccdL-Escf);
	outfile->Printf("\tEdfolccd - Eref (a.u.)             : %20.14f\n", ElccdL-Eref);
	outfile->Printf("\tDF-OLCCD Total Energy (a.u.)       : %20.14f\n", ElccdL);
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\n");

	// Set the global variables with the energies
	Process::environment.globals["CURRENT ENERGY"] = ElccdL;
        Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
        Process::environment.globals["CURRENT CORRELATION ENERGY"] = ElccdL - Escf;
	Process::environment.globals["OLCCD TOTAL ENERGY"] = ElccdL;
        Process::environment.globals["OLCCD CORRELATION ENERGY"] = ElccdL - Escf;

        // OEPROP
        if (oeprop_ == "TRUE") oeprop();

        // Compute Analytic Gradients
        if (dertype == "FIRST") dfgrad();

	// Save MOs to wfn
	save_mo_to_wfn();

  }// end if (conver == 1)

}// end olccd_manager

//======================================================================
//             LCCD Manager
//======================================================================
void DFOCC::lccd_manager()
{

        time4grad = 0;// means i will not compute the gradient
	mo_optimized = 0;// means MOs are not optimized

        timer_on("DF CC Integrals");
        df_corr();
        trans_corr();
        timer_off("DF CC Integrals");
        outfile->Printf("\n\tNumber of basis functions in the DF-CC basis: %3d\n", nQ);

        if (dertype == "FIRST" || oeprop_ == "TRUE" || ekt_ip_ == "TRUE" || qchf_ == "TRUE") {
            timer_on("DF REF Integrals");
            df_ref();
            trans_ref();
            timer_off("DF REF Integrals");
            outfile->Printf("\tNumber of basis functions in the DF-HF basis: %3d\n", nQ_ref);
            Jc = SharedTensor1d(new Tensor1d("DF_BASIS_SCF J_Q", nQ_ref));
        }

        // avaliable mem
        memory = Process::environment.get_memory();
        memory_mb = (double)memory/(1024.0 * 1024.0);
        outfile->Printf("\n\tAvailable memory                      : %9.2lf MB \n", memory_mb);

        // memory requirements
        // DF-CC B(Q,ab) + B(Q,ia) + B(Q,ij)
        cost_df = 0.0;
        cost_df = (navirA * navirA) + (navirA * naoccA) + (naoccA * naoccA);
        cost_df *= nQ;
        cost_df /= 1024.0 * 1024.0;
        cost_df *= sizeof(double);
        if (reference_ == "RESTRICTED") outfile->Printf("\tMemory requirement for 3-index ints   : %9.2lf MB \n", cost_df);
	else if (reference_ == "UNRESTRICTED") outfile->Printf("\tMemory requirement for 3-index ints   : %9.2lf MB \n", 2.0*cost_df);

        // Cost of Integral transform for B(Q,ab)
        cost_ampAA = 0.0;
        cost_ampAA = nQ * nso2_;
        cost_ampAA += nQ * navirA * navirA;
        cost_ampAA += nQ * nso_ * navirA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        outfile->Printf("\tMemory requirement for DF-CC int trans: %9.2lf MB \n", cost_ampAA);

        // Mem for amplitudes
        cost_ampAA = 0.0;
        cost_ampAA = nocc2AA * nvir2AA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        cost_3amp = 3.0 * cost_ampAA;
        cost_4amp = 4.0 * cost_ampAA;
        cost_5amp = 5.0 * cost_ampAA;

        if ((cost_4amp+cost_df) <= memory_mb) {
             outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_4amp);
             outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_4amp+cost_df);
             nincore_amp = 4;
             t2_incore = true;
             df_ints_incore = true;
        }
        else if ((cost_3amp+cost_df) <= memory_mb) {
             outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_3amp);
             //outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_3amp+cost_df);
             outfile->Printf("\tWarning: T2 amplitudes will be stored on the disk!\n");
             nincore_amp = 3;
             t2_incore = false;
             df_ints_incore = false;
        }
        else if (cost_3amp < memory_mb && cost_df < memory_mb ) {
             outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_3amp);
             outfile->Printf("\tWarning: T2 amplitudes will be stored on the disk!\n");
             nincore_amp = 3;
             t2_incore = false;
             df_ints_incore = false;
        }
        else {
             outfile->Printf("\tWarning: There is NOT enough memory for CC contractions!\n");
             outfile->Printf("\tIncrease memory by                    : %9.2lf MB \n", cost_3amp+cost_df-memory_mb);
             throw PSIEXCEPTION("There is NOT enough memory for CC contractions!");
        }

        // W_abef term
        cost_ampAA = 0.0;
        cost_ampAA = naoccA * naoccA * navirA * navirA;
        cost_ampAA += 2.0 * nQ * navirA * navirA;
        cost_ampAA += navirA * navirA * navirA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        double cost_ampAA2 = 0.0;
        cost_ampAA2 = naoccA * naoccA * navirA * navirA;
        cost_ampAA2 += nQ * navirA * navirA;
        cost_ampAA2 += 3.0 * navirA * navirA * navirA;
        cost_ampAA2 /= 1024.0 * 1024.0;
        cost_ampAA2 *= sizeof(double);
        cost_amp = MAX0(cost_ampAA, cost_ampAA2);
        outfile->Printf("\tMemory requirement for Wabef term     : %9.2lf MB \n", cost_amp);

        // memalloc for density intermediates
        if (qchf_ == "TRUE" || dertype == "FIRST") {
            g1Qc = SharedTensor1d(new Tensor1d("DF_BASIS_SCF G1_Q", nQ_ref));
            g1Qt = SharedTensor1d(new Tensor1d("DF_BASIS_SCF G1t_Q", nQ_ref));
            g1Qp = SharedTensor1d(new Tensor1d("DF_BASIS_SCF G1p_Q", nQ_ref));
            g1Q = SharedTensor1d(new Tensor1d("DF_BASIS_CC G1_Q", nQ));
            g1Qt2 = SharedTensor1d(new Tensor1d("DF_BASIS_CC G1t_Q", nQ));
        }

        // QCHF
        if (qchf_ == "TRUE") qchf();

        // Compute MP2 energy
        if (reference == "ROHF") {
            t1A = SharedTensor2d(new Tensor2d("T1_1 <I|A>", naoccA, navirA));
            t1B = SharedTensor2d(new Tensor2d("T1_1 <i|a>", naoccB, navirB));
	    t1_1st_sc();
	}
	lccd_t2_1st_sc();

	outfile->Printf("\n");
	if (reference == "ROHF") outfile->Printf("\tComputing DF-MP2 energy (DF-ROHF-MP2)... \n");
	else outfile->Printf("\tComputing DF-MP2 energy ... \n");
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	outfile->Printf("\tDF-HF Energy (a.u.)                : %20.14f\n", Escf);
	outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA+Escsmp2BB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SCS-MP2 Total Energy (a.u.)     : %20.14f\n", Escsmp2);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SOS-MP2 Total Energy (a.u.)     : %20.14f\n", Esosmp2);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SCSN-MP2 Total Energy (a.u.)    : %20.14f\n", Escsnmp2);
	if (reference == "ROHF") outfile->Printf("\tDF-MP2 Singles Energy (a.u.)       : %20.14f\n", Emp2_t1);
	if (reference == "ROHF") outfile->Printf("\tDF-MP2 Doubles Energy (a.u.)       : %20.14f\n", Ecorr - Emp2_t1);
	outfile->Printf("\tDF-MP2 Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
	outfile->Printf("\tDF-MP2 Total Energy (a.u.)         : %20.14f\n", Emp2);
	outfile->Printf("\t======================================================================= \n");

	Process::environment.globals["MP2 TOTAL ENERGY"] = Emp2;
	Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = Escsmp2;
	Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
	Process::environment.globals["SCSN-MP2 TOTAL ENERGY"] = Escsnmp2;
        Process::environment.globals["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
        Process::environment.globals["SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
        Process::environment.globals["SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
        Process::environment.globals["SCSN-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;
        Process::environment.globals["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
        Process::environment.globals["MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA+Emp2BB;

        // Perform LCCD iterations
        timer_on("LCCD");
        lccd_iterations();
        timer_off("LCCD");

	outfile->Printf("\n");
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\t================ LCCD FINAL RESULTS =================================== \n");
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	outfile->Printf("\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", ElccdAA);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", ElccdBB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", ElccdAB);
	outfile->Printf("\tDF-LCCD Correlation Energy (a.u.)  : %20.14f\n", Ecorr);
	outfile->Printf("\tDF-LCCD Total Energy (a.u.)        : %20.14f\n", Elccd);
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\n");

	Process::environment.globals["CURRENT ENERGY"] = Elccd;
        Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
        Process::environment.globals["CURRENT CORRELATION ENERGY"] = Elccd - Escf;
	Process::environment.globals["LCCD TOTAL ENERGY"] = Elccd;
        Process::environment.globals["LCCD CORRELATION ENERGY"] = Elccd - Escf;
	ElccdL = Elccd;

        // Compute Analytic Gradients
        if (dertype == "FIRST" || ekt_ip_ == "TRUE") {
	    // memalloc
	    gQt = SharedTensor1d(new Tensor1d("CCD PDM G_Qt", nQ));
	    if (reference_ == "RESTRICTED") {
	        G1c_ov = SharedTensor2d(new Tensor2d("Correlation OPDM <O|V>", noccA, nvirA));
	        G1c_vo = SharedTensor2d(new Tensor2d("Correlation OPDM <V|O>", nvirA, noccA));
	    }
	    else if (reference_ == "UNRESTRICTED") {
	        G1c_ovA = SharedTensor2d(new Tensor2d("Correlation OPDM <O|V>", noccA, nvirA));
	        G1c_ovB = SharedTensor2d(new Tensor2d("Correlation OPDM <o|v>", noccB, nvirB));
	        G1c_voA = SharedTensor2d(new Tensor2d("Correlation OPDM <V|O>", nvirA, noccA));
	        G1c_voB = SharedTensor2d(new Tensor2d("Correlation OPDM <v|o>", nvirB, noccB));
	    }

            outfile->Printf("\tComputing unrelaxed response density matrices...\n");
	    lccd_pdm_3index_intr();
 	    omp3_opdm();
	    olccd_tpdm();
	    //ccl_energy();
            prepare4grad();
            if (oeprop_ == "TRUE") oeprop();
            if (dertype == "FIRST") dfgrad();
            //if (ekt_ip_ == "TRUE") ekt_ip();
        }// if (dertype == "FIRST" || ekt_ip_ == "TRUE")

}// end lccd_manager

//======================================================================
//             HF Manager
//======================================================================
void DFOCC::qchf_manager()
{
        time4grad = 0;// means i will not compute the gradient
	mo_optimized = 0;// means MOs are not optimized

	/*
        timer_on("DF CC Integrals");
        df_corr();
        trans_corr();
        timer_off("DF CC Integrals");
	*/

        df_ref();
        //trans_ref();
        outfile->Printf("\tNumber of basis functions in the DF-HF basis: %3d\n", nQ_ref);

        // memalloc for density intermediates
        Jc = SharedTensor1d(new Tensor1d("DF_BASIS_SCF J_Q", nQ_ref));
	/*
        g1Qc = SharedTensor1d(new Tensor1d("DF_BASIS_SCF G1_Q", nQ_ref));
        g1Qt = SharedTensor1d(new Tensor1d("DF_BASIS_SCF G1t_Q", nQ_ref));
        g1Q = SharedTensor1d(new Tensor1d("DF_BASIS_CC G1_Q", nQ));
        g1Qt2 = SharedTensor1d(new Tensor1d("DF_BASIS_CC G1t_Q", nQ));
	*/
        //outfile->Printf("\tNumber of basis functions in the DF-CC basis: %3d\n", nQ);


     if (reference_ == "RESTRICTED") {
        // memory requirements
        // DF-CC B(Q,mn)
        cost_ampAA = 0.0;
        cost_ampAA = nQ_ref * nso2_;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        outfile->Printf("\tMemory requirement for B-CC (Q|mu nu) : %9.2lf MB \n", cost_ampAA);

        // DF-CC B(Q,ab)
        cost_ampAA = 0.0;
        cost_ampAA = nQ_ref * navirA * navirA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        outfile->Printf("\tMemory requirement for B-CC (Q|ab)    : %9.2lf MB \n", cost_ampAA);

        // Cost of Integral transform for DF-CC B(Q,ab)
        cost_ampAA = 0.0;
        cost_ampAA = nQ_ref * nso2_;
        cost_ampAA += nQ_ref * navirA * navirA;
        cost_ampAA += nQ_ref * nso_ * navirA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        outfile->Printf("\tMemory requirement for DF-CC int trans: %9.2lf MB \n", cost_ampAA);
     }  // end if (reference_ == "RESTRICTED")

     else if (reference_ == "UNRESTRICTED") {
        // memory requirements
        memory = Process::environment.get_memory();
        memory_mb = (double)memory/(1024.0 * 1024.0);
        outfile->Printf("\n\tAvailable memory                      : %9.2lf MB \n", memory_mb);
     }  // end else if (reference_ == "UNRESTRICTED")

        // QCHF
        qchf();

}// end qchf_manager


}} // End Namespaces



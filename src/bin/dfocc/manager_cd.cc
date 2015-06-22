/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

#include <libqt/qt.h>
#include "dfocc.h"

using namespace psi;
using namespace std;


namespace psi{ namespace dfoccwave{
  
//======================================================================
//             CD-OMP2 Manager
//======================================================================             
void DFOCC::cd_omp2_manager()
{
        do_cd = "TRUE";
        time4grad = 0;// means i will not compute the gradient
	mo_optimized = 0;// means MOs are not optimized
	orbs_already_opt = 0;// means orbitals are not optimized yet.
	orbs_already_sc = 0;// menas orbitals are not semicanonical yet.
        timer_on("CD Integrals");
        cd_ints();
        timer_off("CD Integrals");
        timer_on("CD Trans");
        trans_cd();
        timer_off("CD Trans");

        // memalloc for density intermediates
        Jc = SharedTensor1d(new Tensor1d("DF_BASIS_SCF J_Q", nQ_ref));
        g1Qc = SharedTensor1d(new Tensor1d("DF_BASIS_SCF G1_Q", nQ_ref));
        g1Qt = SharedTensor1d(new Tensor1d("DF_BASIS_SCF G1t_Q", nQ_ref));
        g1Q = SharedTensor1d(new Tensor1d("DF_BASIS_CC G1_Q", nQ));
        g1Qt2 = SharedTensor1d(new Tensor1d("DF_BASIS_CC G1t_Q", nQ));

        // Fock
        fock();

        // ROHF REF
        if (reference == "ROHF") t1_1st_sc();
	t2_1st_sc();
	Emp2L=Emp2;
        EcorrL=Emp2L-Escf;
	Emp2L_old=Emp2;
	
	outfile->Printf("\n");
	if (reference == "ROHF") outfile->Printf("\tComputing CD-MP2 energy using SCF MOs (CD-ROHF-MP2)... \n"); 
	else outfile->Printf("\tComputing CD-MP2 energy using SCF MOs (Canonical CD-MP2)... \n"); 
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	outfile->Printf("\tCD-HF Energy (a.u.)                : %20.14f\n", Escf);
	outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA+Escsmp2BB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SCS-MP2 Total Energy (a.u.)     : %20.14f\n", Escsmp2);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SOS-MP2 Total Energy (a.u.)     : %20.14f\n", Esosmp2);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SCSN-MP2 Total Energy (a.u.)    : %20.14f\n", Escsnmp2);
	if (reference == "ROHF") outfile->Printf("\tCD-MP2 Singles Energy (a.u.)       : %20.14f\n", Emp2_t1);
	if (reference == "ROHF") outfile->Printf("\tCD-MP2 Doubles Energy (a.u.)       : %20.14f\n", Ecorr - Emp2_t1);
	outfile->Printf("\tCD-MP2 Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
	outfile->Printf("\tCD-MP2 Total Energy (a.u.)         : %20.14f\n", Emp2);
	outfile->Printf("\t======================================================================= \n");
	
	Process::environment.globals["CURRENT ENERGY"] = Emp2;
	Process::environment.globals["CD-MP2 TOTAL ENERGY"] = Emp2;
	Process::environment.globals["CD-SCS-MP2 TOTAL ENERGY"] = Escsmp2;
	Process::environment.globals["CD-SOS-MP2 TOTAL ENERGY"] = Esosmp2;
	Process::environment.globals["CD-SCSN-MP2 TOTAL ENERGY"] = Escsnmp2;

        Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
        Process::environment.globals["CURRENT CORRELATION ENERGY"] = Emp2 - Escf;
        Process::environment.globals["CD-MP2 CORRELATION ENERGY"] = Emp2 - Escf;
        Process::environment.globals["CD-SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
        Process::environment.globals["CD-SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
        Process::environment.globals["CD-SCSN-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;

        Process::environment.globals["CD-MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
        Process::environment.globals["CD-MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA+Emp2BB;

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
	
        if (rms_wog <= tol_grad && fabs(DE) >= tol_Eod) {
           orbs_already_opt = 1;
	   if (conver == 1) outfile->Printf("\n\tOrbitals are optimized now.\n");
	   else if (conver == 0) { 
                    outfile->Printf("\n\tMAX MOGRAD did NOT converged, but RMS MOGRAD converged!!!\n");
	            outfile->Printf("\tI will consider the present orbitals as optimized.\n");
           }
	   outfile->Printf("\tSwitching to the standard CD-MP2 computation after semicanonicalization of the MOs... \n");
	   
	   semi_canonic();
           trans_cd();
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
        } 

  if (conver == 1) {
        ref_energy();
	mp2_energy();
        if (orbs_already_opt == 1) Emp2L = Emp2;
	
	outfile->Printf("\n");
	outfile->Printf("\tComputing MP2 energy using optimized MOs... \n");
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	outfile->Printf("\tCD-HF Energy (a.u.)                : %20.14f\n", Escf);
	outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA+Escsmp2BB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SCS-MP2 Total Energy (a.u.)     : %20.14f\n", Escsmp2);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SOS-MP2 Total Energy (a.u.)     : %20.14f\n", Esosmp2);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SCSN-MP2 Total Energy (a.u.)    : %20.14f\n", Escsnmp2);
	outfile->Printf("\tCD-MP2 Correlation Energy (a.u.)   : %20.14f\n", Emp2 - Escf);
	outfile->Printf("\tCD-MP2 Total Energy (a.u.)         : %20.14f\n", Emp2);
	outfile->Printf("\t======================================================================= \n");
	

	outfile->Printf("\n");
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\t================ CD-OMP2 FINAL RESULTS ================================ \n");
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	outfile->Printf("\tCD-HF Energy (a.u.)                : %20.14f\n", Escf);
	outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SCS-OMP2 Total Energy (a.u.)    : %20.14f\n", Escsmp2);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SOS-OMP2 Total Energy (a.u.)    : %20.14f\n", Esosmp2);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SCSN-OMP2 Total Energy (a.u.)   : %20.14f\n", Escsnmp2);
	outfile->Printf("\tCD-OMP2 Correlation Energy (a.u.)  : %20.14f\n", Emp2L-Escf);
	outfile->Printf("\tEcdomp2 - Eref (a.u.)              : %20.14f\n", Emp2L-Eref);
	outfile->Printf("\tCD-OMP2 Total Energy (a.u.)        : %20.14f\n", Emp2L);
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\n");
	

	// Set the global variables with the energies
	Process::environment.globals["CD-OMP2 TOTAL ENERGY"] = Emp2L;
	Process::environment.globals["CD-SCS-OMP2 TOTAL ENERGY"] =  Escsmp2;
	Process::environment.globals["CD-SOS-OMP2 TOTAL ENERGY"] =  Esosmp2;
	Process::environment.globals["CD-SCSN-OMP2 TOTAL ENERGY"] = Escsnmp2;
	Process::environment.globals["CURRENT ENERGY"] = Emp2L;
	Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
	Process::environment.globals["CURRENT CORRELATION ENERGY"] = Emp2L-Escf;

        Process::environment.globals["CD-OMP2 CORRELATION ENERGY"] = Emp2L - Escf;
        Process::environment.globals["CD-SCS-OMP2 CORRELATION ENERGY"] =  Escsmp2 - Escf;
        Process::environment.globals["CD-SOS-OMP2 CORRELATION ENERGY"] =  Esosmp2 - Escf;
        Process::environment.globals["CD-SCSN-OMP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;

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

        // OEPROP
        if (oeprop_ == "TRUE") oeprop();
 
	//if (natorb == "TRUE") nbo();
	//if (occ_orb_energy == "TRUE") semi_canonic(); 

        // Compute Analytic Gradients
        /*
        if (dertype == "FIRST") {
	    outfile->Printf("\tAnalytic gradient computation is starting...\n");
	    
            coord_grad();
	    outfile->Printf("\tNecessary information has been sent to DERIV, which will take care of the rest.\n");
	    
        }
        */

  }// end if (conver == 1)
}// end omp2_manager 

//======================================================================
//             CD-MP2 Manager
//======================================================================             
void DFOCC::cd_mp2_manager()
{
        do_cd = "TRUE";
        time4grad = 0;// means i will not compute the gradient
	mo_optimized = 0;// means MOs are not optimized
        timer_on("CD Integrals");
        cd_ints();
        if (dertype == "NONE" && ekt_ip_ == "FALSE") {
            trans_cd_mp2();
        }
        else trans_cd();
        timer_off("CD Integrals");

        // ROHF REF
        //outfile->Printf("\tI am here.\n"); 
        if (reference == "ROHF") t1_1st_sc();

	//t2_1st_sc();
        //mp2_energy();
        if (dertype == "NONE" && ekt_ip_ == "FALSE") mp2_direct();
        else {
	     t2_1st_sc();
             mp2_energy();
        }
	Emp2L=Emp2;
        EcorrL=Emp2L-Escf;
	
	outfile->Printf("\n");
	if (reference == "ROHF") outfile->Printf("\tComputing CD-MP2 energy using SCF MOs (CD-ROHF-MP2)... \n"); 
	else outfile->Printf("\tComputing CD-MP2 energy using SCF MOs (Canonical CD-MP2)... \n"); 
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	outfile->Printf("\tCD-HF Energy (a.u.)                : %20.14f\n", Escf);
	outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA+Escsmp2BB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SCS-MP2 Total Energy (a.u.)     : %20.14f\n", Escsmp2);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SOS-MP2 Total Energy (a.u.)     : %20.14f\n", Esosmp2);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SCSN-MP2 Total Energy (a.u.)    : %20.14f\n", Escsnmp2);
	if (reference_ == "ROHF") outfile->Printf("\tCD-MP2 Singles Energy (a.u.)       : %20.14f\n", Emp2_t1);
	if (reference_ == "ROHF") outfile->Printf("\tCD-MP2 Doubles Energy (a.u.)       : %20.14f\n", Ecorr - Emp2_t1);
	outfile->Printf("\tCD-MP2 Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
	outfile->Printf("\tCD-MP2 Total Energy (a.u.)         : %20.14f\n", Emp2);
	outfile->Printf("\t======================================================================= \n");
	
	Process::environment.globals["CURRENT ENERGY"] = Emp2;
	Process::environment.globals["CD-MP2 TOTAL ENERGY"] = Emp2;
	Process::environment.globals["CD-SCS-MP2 TOTAL ENERGY"] = Escsmp2;
	Process::environment.globals["CD-SOS-MP2 TOTAL ENERGY"] = Esosmp2;
	Process::environment.globals["CD-SCSN-MP2 TOTAL ENERGY"] = Escsnmp2;

        Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
        Process::environment.globals["CURRENT CORRELATION ENERGY"] = Emp2 - Escf;
        Process::environment.globals["CD-MP2 CORRELATION ENERGY"] = Emp2 - Escf;
        Process::environment.globals["CD-SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
        Process::environment.globals["CD-SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
        Process::environment.globals["CD-SCSN-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;

        Process::environment.globals["CD-MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
        Process::environment.globals["CD-MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA+Emp2BB;


}// end mp2_manager 

//======================================================================
//             CCSD Manager
//======================================================================             
void DFOCC::ccsd_manager_cd()
{

        do_cd = "TRUE";
        time4grad = 0;// means i will not compute the gradient
	mo_optimized = 0;// means MOs are not optimized

        timer_on("CD Integrals");
        cd_ints();
        trans_cd();
        timer_off("CD Integrals");

        // Memory allocation
        T1c = SharedTensor1d(new Tensor1d("DF_BASIS_CC T1_Q", nQ));
	Jc = SharedTensor1d(new Tensor1d("DF_BASIS_SCF J_Q", nQ_ref));

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
	if (reference == "ROHF") outfile->Printf("\tComputing CD-MP2 energy (CD-ROHF-MP2)... \n"); 
	else outfile->Printf("\tComputing DF-MP2 energy ... \n"); 
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	outfile->Printf("\tCD-HF Energy (a.u.)                : %20.14f\n", Escf);
	outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA+Escsmp2BB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SCS-MP2 Total Energy (a.u.)     : %20.14f\n", Escsmp2);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SOS-MP2 Total Energy (a.u.)     : %20.14f\n", Esosmp2);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SCSN-MP2 Total Energy (a.u.)    : %20.14f\n", Escsnmp2);
	if (reference_ == "ROHF") outfile->Printf("\tCD-MP2 Singles Energy (a.u.)       : %20.14f\n", Emp2_t1);
	if (reference_ == "ROHF") outfile->Printf("\tCD-MP2 Doubles Energy (a.u.)       : %20.14f\n", Ecorr - Emp2_t1);
	outfile->Printf("\tCD-MP2 Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
	outfile->Printf("\tCD-MP2 Total Energy (a.u.)         : %20.14f\n", Emp2);
	outfile->Printf("\t======================================================================= \n");
	
	Process::environment.globals["CD-MP2 TOTAL ENERGY"] = Emp2;
	Process::environment.globals["CD-SCS-MP2 TOTAL ENERGY"] = Escsmp2;
	Process::environment.globals["CD-SOS-MP2 TOTAL ENERGY"] = Esosmp2;
	Process::environment.globals["CD-SCSN-MP2 TOTAL ENERGY"] = Escsnmp2;
        Process::environment.globals["CD-MP2 CORRELATION ENERGY"] = Emp2 - Escf;
        Process::environment.globals["CD-SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
        Process::environment.globals["CD-SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
        Process::environment.globals["CD-SCSN-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;
        Process::environment.globals["CD-MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
        Process::environment.globals["CD-MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA+Emp2BB;

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
	outfile->Printf("\tCD-CCSD Correlation Energy (a.u.)  : %20.14f\n", Ecorr);
	outfile->Printf("\tCD-CCSD Total Energy (a.u.)        : %20.14f\n", Eccsd);
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\n");
	
	Process::environment.globals["CURRENT ENERGY"] = Eccsd;
        Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
        Process::environment.globals["CURRENT CORRELATION ENERGY"] = Eccsd - Escf;
	Process::environment.globals["CD-CCSD TOTAL ENERGY"] = Eccsd;
        Process::environment.globals["CD-CCSD CORRELATION ENERGY"] = Eccsd - Escf;

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
        if (dertype == "FIRST" || oeprop_ == "TRUE" || ekt_ip_ == "TRUE") {
	    // memalloc
	    G1c_ov = SharedTensor2d(new Tensor2d("Correlation OPDM <O|V>", noccA, nvirA));
	    G1c_vo = SharedTensor2d(new Tensor2d("Correlation OPDM <V|O>", nvirA, noccA));

            outfile->Printf("\tComputing unrelaxed response density matrices...\n");
 	    ccsd_opdm();
	    ccsd_tpdm();
	    ccl_energy();
            prepare4grad();
            //if (oeprop_ == "TRUE") oeprop();
            //if (dertype == "FIRST") dfgrad();
            //if (ekt_ip_ == "TRUE") ekt_ip(); 
        }// if (dertype == "FIRST" || ekt_ip_ == "TRUE") 

}// end ccsd_manager_cd 

//======================================================================
//             CCD Manager
//======================================================================             
void DFOCC::ccd_manager_cd()
{

	do_cd = "TRUE";
        time4grad = 0;// means i will not compute the gradient
	mo_optimized = 0;// means MOs are not optimized

        timer_on("CD Integrals");
        cd_ints();
        trans_cd();
        timer_off("CD Integrals");

        // Memory allocation
        T1c = SharedTensor1d(new Tensor1d("DF_BASIS_CC T1_Q", nQ));
	Jc = SharedTensor1d(new Tensor1d("DF_BASIS_SCF J_Q", nQ_ref));


     if (reference_ == "RESTRICTED") {
        // avaliable mem
        memory = Process::environment.get_memory();
        memory_mb = (double)memory/(1024.0 * 1024.0);
        outfile->Printf("\n\tAvailable memory                      : %9.2lf MB \n", memory_mb);

        // memory requirements
        /*
        // DF-HF B(Q,mn)
        cost_ampAA = 0;
        cost_ampAA = nQ_ref * nso2_;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        cost_amp = 3.0 * cost_ampAA;
        outfile->Printf("\n\tMemory requirement for B-HF (Q|mu nu) : %9.2lf MB \n", cost_ampAA);
        outfile->Printf("\tMemory requirement for 3*B-HF (Q|mu nu) %9.2lf MB \n", cost_amp);
 
        // DF-HF B(Q,ab)
        cost_ampAA = 0;
        cost_ampAA = nQ_ref * nvir2AA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        cost_amp = 3.0 * cost_ampAA;
        outfile->Printf("\tMemory requirement for B-HF (Q|ab)    : %9.2lf MB \n", cost_ampAA);
        outfile->Printf("\tMemory requirement for 3*B-HF (Q|ab)  : %9.2lf MB \n", cost_amp);

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
        */

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
	if (reference == "ROHF") outfile->Printf("\tComputing CD-MP2 energy (CD-ROHF-MP2)... \n"); 
	else outfile->Printf("\tComputing CD-MP2 energy ... \n"); 
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	outfile->Printf("\tCD-HF Energy (a.u.)                : %20.14f\n", Escf);
	outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA+Escsmp2BB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SCS-MP2 Total Energy (a.u.)     : %20.14f\n", Escsmp2);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SOS-MP2 Total Energy (a.u.)     : %20.14f\n", Esosmp2);
	if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SCSN-MP2 Total Energy (a.u.)    : %20.14f\n", Escsnmp2);
	if (reference_ == "ROHF") outfile->Printf("\tCD-MP2 Singles Energy (a.u.)       : %20.14f\n", Emp2_t1);
	if (reference_ == "ROHF") outfile->Printf("\tCD-MP2 Doubles Energy (a.u.)       : %20.14f\n", Ecorr - Emp2_t1);
	outfile->Printf("\tCD-MP2 Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
	outfile->Printf("\tCD-MP2 Total Energy (a.u.)         : %20.14f\n", Emp2);
	outfile->Printf("\t======================================================================= \n");
	
	Process::environment.globals["CD-MP2 TOTAL ENERGY"] = Emp2;
	Process::environment.globals["CD-SCS-MP2 TOTAL ENERGY"] = Escsmp2;
	Process::environment.globals["CD-SOS-MP2 TOTAL ENERGY"] = Esosmp2;
	Process::environment.globals["CD-SCSN-MP2 TOTAL ENERGY"] = Escsnmp2;
        Process::environment.globals["CD-MP2 CORRELATION ENERGY"] = Emp2 - Escf;
        Process::environment.globals["CD-SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
        Process::environment.globals["CD-SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
        Process::environment.globals["CD-SCSN-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;
        Process::environment.globals["CD-MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
        Process::environment.globals["CD-MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA+Emp2BB;

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
	outfile->Printf("\tCD-CCD Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
	outfile->Printf("\tCD-CCD Total Energy (a.u.)         : %20.14f\n", Eccd);
	outfile->Printf("\t======================================================================= \n");
	outfile->Printf("\n");
	
	Process::environment.globals["CURRENT ENERGY"] = Eccd;
        Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
        Process::environment.globals["CURRENT CORRELATION ENERGY"] = Eccd - Escf;
	Process::environment.globals["CD-CCD TOTAL ENERGY"] = Eccd;
        Process::environment.globals["CD-CCD CORRELATION ENERGY"] = Eccd - Escf;

        // CCDL 
        if (dertype == "FIRST" || cc_lambda_ == "TRUE") {
	    // memalloc
            if (dertype == "FIRST") {
	        gQt = SharedTensor1d(new Tensor1d("CCSD PDM G_Qt", nQ));
            }

            timer_on("CCDL");
            if (t2_incore) ccdl_iterations();
            else throw PSIEXCEPTION("There is NOT enough memory for Lambda equations!");
            timer_off("CCDL");
        }

        // Compute Analytic Gradients
        if (dertype == "FIRST" || oeprop_ == "TRUE" || ekt_ip_ == "TRUE") {
	    // memalloc
	    G1c_ov = SharedTensor2d(new Tensor2d("Correlation OPDM <O|V>", noccA, nvirA));
	    G1c_vo = SharedTensor2d(new Tensor2d("Correlation OPDM <V|O>", nvirA, noccA));

            outfile->Printf("\tComputing unrelaxed response density matrices...\n");
 	    ccd_opdm();
	    ccd_tpdm();
	    ccl_energy();
            prepare4grad();
            //if (oeprop_ == "TRUE") oeprop();
            //if (dertype == "FIRST") dfgrad();
            //if (ekt_ip_ == "TRUE") ekt_ip(); 
        }// if (dertype == "FIRST" || ekt_ip_ == "TRUE") 

}// end ccd_manager_cd 


}} // End Namespaces




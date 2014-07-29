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
        psi::fprintf(outfile,"\tNumber of basis functions in the DF-HF basis: %3d\n", nQ_ref);
        psi::fprintf(outfile,"\tNumber of basis functions in the DF-CC basis: %3d\n", nQ);
        fflush(outfile);
        timer_off("DF CC Integrals");

        // memalloc for density intermediates
        Jc = SharedTensor1d(new Tensor1d("DF_BASIS_SCF J_Q", nQ_ref));
        g1Qc = SharedTensor1d(new Tensor1d("DF_BASIS_SCF G1_Q", nQ_ref));
        g1Qt = SharedTensor1d(new Tensor1d("DF_BASIS_SCF G1t_Q", nQ_ref));
        g1Q = SharedTensor1d(new Tensor1d("DF_BASIS_CC G1_Q", nQ));
        g1Qt2 = SharedTensor1d(new Tensor1d("DF_BASIS_CC G1t_Q", nQ));

        if (conv_tei_type == "DISK") { 
           tei_oooo_chem_ref();
           tei_ooov_chem_ref();
           tei_oovv_chem_ref();
           tei_ovov_chem_ref(); 
           if (reference_ == "UNRESTRICTED") {
            tei_oooo_phys_ref();
            tei_ooov_phys_ref();
            tei_oovv_phys_ref();
            tei_ovov_phys_ref();
            tei_oooo_anti_symm_ref();
            tei_ooov_anti_symm_ref();
            tei_oovv_anti_symm_ref();
            tei_ovov_anti_symm_ref();
           }

           tei_iajb_chem();
           //tei_ijab_chem();// for Hessian
           if (reference_ == "UNRESTRICTED") {
               tei_ijab_phys();
               tei_ijab_anti_symm();
           }
        }// if (conv_tei_type == "DISK")  
           fock();

        // ROHF REF
        if (reference == "ROHF") t1_1st_sc();
	t2_1st_sc();
	Emp2L=Emp2;
        EcorrL=Emp2L-Escf;
	Emp2L_old=Emp2;
	
	psi::fprintf(outfile,"\n");
	if (reference == "ROHF") psi::fprintf(outfile,"\tComputing DF-MP2 energy using SCF MOs (DF-ROHF-MP2)... \n"); 
	else psi::fprintf(outfile,"\tComputing DF-MP2 energy using SCF MOs (Canonical DF-MP2)... \n"); 
	psi::fprintf(outfile,"\t======================================================================= \n");
	psi::fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	psi::fprintf(outfile,"\tDF-HF Energy (a.u.)                : %20.14f\n", Escf);
	psi::fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	if (reference_ == "UNRESTRICTED") psi::fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
	if (reference_ == "UNRESTRICTED") psi::fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
	if (reference_ == "UNRESTRICTED") psi::fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
	if (reference_ == "UNRESTRICTED") psi::fprintf(outfile,"\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA+Escsmp2BB);
	if (reference_ == "UNRESTRICTED") psi::fprintf(outfile,"\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
	if (reference_ == "UNRESTRICTED") psi::fprintf(outfile,"\tDF-SCS-MP2 Total Energy (a.u.)     : %20.14f\n", Escsmp2);
	if (reference_ == "UNRESTRICTED") psi::fprintf(outfile,"\tDF-SOS-MP2 Total Energy (a.u.)     : %20.14f\n", Esosmp2);
	if (reference_ == "UNRESTRICTED") psi::fprintf(outfile,"\tDF-SCSN-MP2 Total Energy (a.u.)    : %20.14f\n", Escsnmp2);
	if (reference == "ROHF") psi::fprintf(outfile,"\tDF-MP2 Singles Energy (a.u.)       : %20.14f\n", Emp2_t1);
	if (reference == "ROHF") psi::fprintf(outfile,"\tDF-MP2 Doubles Energy (a.u.)       : %20.14f\n", Ecorr - Emp2_t1);
	psi::fprintf(outfile,"\tDF-MP2 Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
	psi::fprintf(outfile,"\tDF-MP2 Total Energy (a.u.)         : %20.14f\n", Emp2);
	psi::fprintf(outfile,"\t======================================================================= \n");
	fflush(outfile);
	Process::environment.globals["CURRENT ENERGY"] = Emp2;
	Process::environment.globals["DF-MP2 TOTAL ENERGY"] = Emp2;
	Process::environment.globals["DF-SCS-MP2 TOTAL ENERGY"] = Escsmp2;
	Process::environment.globals["DF-SOS-MP2 TOTAL ENERGY"] = Esosmp2;
	Process::environment.globals["DF-SCSN-MP2 TOTAL ENERGY"] = Escsnmp2;

        Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
        Process::environment.globals["CURRENT CORRELATION ENERGY"] = Emp2 - Escf;
        Process::environment.globals["DF-MP2 CORRELATION ENERGY"] = Emp2 - Escf;
        Process::environment.globals["DF-SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
        Process::environment.globals["DF-SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
        Process::environment.globals["DF-SCSN-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;

        Process::environment.globals["DF-MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
        Process::environment.globals["DF-MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA+Emp2BB;

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
	
        if (rms_wog <= tol_grad && fabs(DE) >= tol_Eod) {
           orbs_already_opt = 1;
	   if (conver == 1) psi::fprintf(outfile,"\n\tOrbitals are optimized now.\n");
	   else if (conver == 0) { 
                    psi::fprintf(outfile,"\n\tMAX MOGRAD did NOT converged, but RMS MOGRAD converged!!!\n");
	            psi::fprintf(outfile,"\tI will consider the present orbitals as optimized.\n");
           }
	   psi::fprintf(outfile,"\tTransforming MOs to the semicanonical basis... \n");
	   fflush(outfile);
	   semi_canonic();
	   psi::fprintf(outfile,"\tSwitching to the standard DF-MP2 computation... \n");
	   fflush(outfile);
           trans_corr();
        if (conv_tei_type == "DISK") { 
           tei_iajb_chem();
           if (reference_ == "UNRESTRICTED") {
               tei_ijab_phys();
               tei_ijab_anti_symm();
           }
        }// if (conv_tei_type == "DISK") 
           trans_ref();
        if (conv_tei_type == "DISK") { 
           tei_oooo_chem_ref();
           tei_ooov_chem_ref();
           tei_oovv_chem_ref();
           tei_ovov_chem_ref();
           if (reference_ == "UNRESTRICTED") {
            tei_oooo_phys_ref();
            tei_ooov_phys_ref();
            tei_oovv_phys_ref();
            tei_ovov_phys_ref();
            tei_oooo_anti_symm_ref();
            tei_ooov_anti_symm_ref();
            tei_oovv_anti_symm_ref();
            tei_ovov_anti_symm_ref();
           }
        }// if (conv_tei_type == "DISK") 
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
	
	psi::fprintf(outfile,"\n");
	psi::fprintf(outfile,"\tComputing MP2 energy using optimized MOs... \n");
	psi::fprintf(outfile,"\t======================================================================= \n");
	psi::fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	psi::fprintf(outfile,"\tDF-HF Energy (a.u.)                : %20.14f\n", Escf);
	psi::fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	if (reference_ == "UNRESTRICTED") psi::fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
	if (reference_ == "UNRESTRICTED") psi::fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
	if (reference_ == "UNRESTRICTED") psi::fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
	if (reference_ == "UNRESTRICTED") psi::fprintf(outfile,"\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA+Escsmp2BB);
	if (reference_ == "UNRESTRICTED") psi::fprintf(outfile,"\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
	if (reference_ == "UNRESTRICTED") psi::fprintf(outfile,"\tDF-SCS-MP2 Total Energy (a.u.)     : %20.14f\n", Escsmp2);
	if (reference_ == "UNRESTRICTED") psi::fprintf(outfile,"\tDF-SOS-MP2 Total Energy (a.u.)     : %20.14f\n", Esosmp2);
	if (reference_ == "UNRESTRICTED") psi::fprintf(outfile,"\tDF-SCSN-MP2 Total Energy (a.u.)    : %20.14f\n", Escsnmp2);
	psi::fprintf(outfile,"\tDF-MP2 Correlation Energy (a.u.)   : %20.14f\n", Emp2 - Escf);
	psi::fprintf(outfile,"\tDF-MP2 Total Energy (a.u.)         : %20.14f\n", Emp2);
	psi::fprintf(outfile,"\t======================================================================= \n");
	fflush(outfile);

	psi::fprintf(outfile,"\n");
	psi::fprintf(outfile,"\t======================================================================= \n");
	psi::fprintf(outfile,"\t================ DF-OMP2 FINAL RESULTS ================================ \n");
	psi::fprintf(outfile,"\t======================================================================= \n");
	psi::fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	psi::fprintf(outfile,"\tDF-HF Energy (a.u.)                : %20.14f\n", Escf);
	psi::fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	if (reference_ == "UNRESTRICTED") psi::fprintf(outfile,"\tDF-SCS-OMP2 Total Energy (a.u.)    : %20.14f\n", Escsmp2);
	if (reference_ == "UNRESTRICTED") psi::fprintf(outfile,"\tDF-SOS-OMP2 Total Energy (a.u.)    : %20.14f\n", Esosmp2);
	if (reference_ == "UNRESTRICTED") psi::fprintf(outfile,"\tDF-SCSN-OMP2 Total Energy (a.u.)   : %20.14f\n", Escsnmp2);
	psi::fprintf(outfile,"\tDF-OMP2 Correlation Energy (a.u.)  : %20.14f\n", Emp2L-Escf);
	psi::fprintf(outfile,"\tEdfomp2 - Eref (a.u.)              : %20.14f\n", Emp2L-Eref);
	psi::fprintf(outfile,"\tDF-OMP2 Total Energy (a.u.)        : %20.14f\n", Emp2L);
	psi::fprintf(outfile,"\t======================================================================= \n");
	psi::fprintf(outfile,"\n");
	fflush(outfile);

	// Set the global variables with the energies
	Process::environment.globals["DF-OMP2 TOTAL ENERGY"] = Emp2L;
	Process::environment.globals["DF-SCS-OMP2 TOTAL ENERGY"] =  Escsmp2;
	Process::environment.globals["DF-SOS-OMP2 TOTAL ENERGY"] =  Esosmp2;
	Process::environment.globals["DF-SCSN-OMP2 TOTAL ENERGY"] = Escsnmp2;
	Process::environment.globals["CURRENT ENERGY"] = Emp2L;
	Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
	Process::environment.globals["CURRENT CORRELATION ENERGY"] = Emp2L-Escf;

        Process::environment.globals["DF-OMP2 CORRELATION ENERGY"] = Emp2L - Escf;
        Process::environment.globals["DF-SCS-OMP2 CORRELATION ENERGY"] =  Escsmp2 - Escf;
        Process::environment.globals["DF-SOS-OMP2 CORRELATION ENERGY"] =  Esosmp2 - Escf;
        Process::environment.globals["DF-SCSN-OMP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;

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
        if (dertype == "NONE" && oeprop_ == "FALSE" && ekt_ip_ == "FALSE" && ekt_ea_ == "FALSE" && comput_s2_ == "FALSE") {
            trans_mp2();
        }
        else {
            trans_corr();
            df_ref();
            trans_ref();
            psi::fprintf(outfile,"\tNumber of basis functions in the DF-HF basis: %3d\n", nQ_ref);

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
        psi::fprintf(outfile,"\tNumber of basis functions in the DF-CC basis: %3d\n", nQ);
        fflush(outfile);
        timer_off("DF CC Integrals");

        // ROHF REF
        //psi::fprintf(outfile,"\tI am here.\n"); fflush(outfile);
        if (reference == "ROHF") t1_1st_sc();

        if (conv_tei_type == "DISK") { 
            tei_iajb_chem();
            if (reference_ == "UNRESTRICTED") {
                tei_ijab_phys();
                tei_ijab_anti_symm();
            }
        }// if (conv_tei_type == "DISK")  
        if (dertype == "NONE" && oeprop_ == "FALSE" && ekt_ip_ == "FALSE" && ekt_ea_ == "FALSE" && comput_s2_ == "FALSE") mp2_direct();
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
	
	psi::fprintf(outfile,"\n");
	if (reference == "ROHF") psi::fprintf(outfile,"\tComputing DF-MP2 energy using SCF MOs (DF-ROHF-MP2)... \n"); 
	else psi::fprintf(outfile,"\tComputing DF-MP2 energy using SCF MOs (Canonical DF-MP2)... \n"); 
	psi::fprintf(outfile,"\t======================================================================= \n");
	psi::fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	psi::fprintf(outfile,"\tDF-HF Energy (a.u.)                : %20.14f\n", Escf);
	psi::fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	if (reference_ == "UNRESTRICTED") psi::fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
	if (reference_ == "UNRESTRICTED") psi::fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
	if (reference_ == "UNRESTRICTED") psi::fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
	if (reference_ == "UNRESTRICTED") psi::fprintf(outfile,"\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA+Escsmp2BB);
	if (reference_ == "UNRESTRICTED") psi::fprintf(outfile,"\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
	if (reference_ == "UNRESTRICTED") psi::fprintf(outfile,"\tDF-SCS-MP2 Total Energy (a.u.)     : %20.14f\n", Escsmp2);
	if (reference_ == "UNRESTRICTED") psi::fprintf(outfile,"\tDF-SOS-MP2 Total Energy (a.u.)     : %20.14f\n", Esosmp2);
	if (reference_ == "UNRESTRICTED") psi::fprintf(outfile,"\tDF-SCSN-MP2 Total Energy (a.u.)    : %20.14f\n", Escsnmp2);
	if (reference_ == "ROHF") psi::fprintf(outfile,"\tDF-MP2 Singles Energy (a.u.)       : %20.14f\n", Emp2_t1);
	if (reference_ == "ROHF") psi::fprintf(outfile,"\tDF-MP2 Doubles Energy (a.u.)       : %20.14f\n", Ecorr - Emp2_t1);
	psi::fprintf(outfile,"\tDF-MP2 Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
	psi::fprintf(outfile,"\tDF-MP2 Total Energy (a.u.)         : %20.14f\n", Emp2);
	psi::fprintf(outfile,"\t======================================================================= \n");
	fflush(outfile);
	Process::environment.globals["CURRENT ENERGY"] = Emp2;
	Process::environment.globals["DF-MP2 TOTAL ENERGY"] = Emp2;
	Process::environment.globals["DF-SCS-MP2 TOTAL ENERGY"] = Escsmp2;
	Process::environment.globals["DF-SOS-MP2 TOTAL ENERGY"] = Esosmp2;
	Process::environment.globals["DF-SCSN-MP2 TOTAL ENERGY"] = Escsnmp2;

        Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
        Process::environment.globals["CURRENT CORRELATION ENERGY"] = Emp2 - Escf;
        Process::environment.globals["DF-MP2 CORRELATION ENERGY"] = Emp2 - Escf;
        Process::environment.globals["DF-SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
        Process::environment.globals["DF-SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
        Process::environment.globals["DF-SCSN-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;

        Process::environment.globals["DF-MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
        Process::environment.globals["DF-MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA+Emp2BB;

        // S2
        //if (comput_s2_ == "TRUE" && reference_ == "UNRESTRICTED" && dertype == "NONE") s2_response();
        if (comput_s2_ == "TRUE" && reference_ == "UNRESTRICTED") {
            if (reference == "UHF" || reference == "UKS") s2_response();
        }

        // Compute Analytic Gradients
        if (dertype == "FIRST" || oeprop_ == "TRUE" || ekt_ip_ == "TRUE" || ekt_ea_ == "TRUE") {
            psi::fprintf(outfile,"\n\tComputing unrelaxed response density matrices...\n");
            fflush(outfile);
 	    omp2_opdm();
	    omp2_tpdm();
            prepare4grad();
            if (oeprop_ == "TRUE") oeprop();
            if (dertype == "FIRST") dfgrad();

            /*
            if (ekt_ip_ == "TRUE" && ekt_ea_ == "TRUE") {
                ekt_ip();
                ekt_ea();
            }

            else if (ekt_ip_ == "TRUE" && ekt_ea_ == "FALSE") {
                ekt_ip();
            }

            else if (ekt_ip_ == "FALSE" && ekt_ea_ == "TRUE") {
                ekt_ea();
            }
            */
        }// if (dertype == "FIRST" || ekt_ip_ == "TRUE" || ekt_ea_ == "TRUE") 

}// end mp2_manager 


//======================================================================
//             OMP3 Manager
//======================================================================             
void DFOCC::omp3_manager()
{
        /*
	mo_optimized = 0;
	orbs_already_opt = 0;
	orbs_already_sc = 0;
        time4grad = 0;// means i will not compute the gradient
        timer_on("trans_ints");
	if (reference_ == "RESTRICTED") trans_ints_rhf();  
	else if (reference_ == "UNRESTRICTED") trans_ints_uhf();  
        timer_off("trans_ints");
        timer_on("REF Energy");
	ref_energy();
        timer_off("REF Energy");
        timer_on("T2(1)");
	omp3_t2_1st_sc();
        timer_off("T2(1)");
        timer_on("MP2 Energy");
	omp3_mp2_energy();
        timer_off("MP2 Energy");

	psi::fprintf(outfile,"\n"); 
	psi::fprintf(outfile,"\tComputing MP2 energy using SCF MOs (Canonical MP2)... \n"); 
	psi::fprintf(outfile,"\t============================================================================== \n");
	psi::fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	psi::fprintf(outfile,"\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	psi::fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	psi::fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
	psi::fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
	psi::fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
	psi::fprintf(outfile,"\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA+Escsmp2BB);
	psi::fprintf(outfile,"\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
	psi::fprintf(outfile,"\tSCS-MP2 Total Energy (a.u.)        : %20.14f\n", Escsmp2);
	psi::fprintf(outfile,"\tSOS-MP2 Total Energy (a.u.)        : %20.14f\n", Esosmp2);
	psi::fprintf(outfile,"\tSCSN-MP2 Total Energy (a.u.)       : %20.14f\n", Escsnmp2);
	psi::fprintf(outfile,"\tSCS-MP2-VDW Total Energy (a.u.)    : %20.14f\n", Escsmp2vdw);
	psi::fprintf(outfile,"\tSOS-PI-MP2 Total Energy (a.u.)     : %20.14f\n", Esospimp2);
	psi::fprintf(outfile,"\tMP2 Correlation Energy (a.u.)      : %20.14f\n", Ecorr);
	psi::fprintf(outfile,"\tMP2 Total Energy (a.u.)            : %20.14f\n", Emp2);
	psi::fprintf(outfile,"\t============================================================================== \n");
	psi::fprintf(outfile,"\n"); 
	fflush(outfile);
	Process::environment.globals["MP2 TOTAL ENERGY"] = Emp2;
	Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = Escsmp2;
	Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
	Process::environment.globals["SCSN-MP2 TOTAL ENERGY"] = Escsnmp2;
	Process::environment.globals["SCS-MP2-VDW TOTAL ENERGY"] = Escsmp2vdw;
	Process::environment.globals["SOS-PI-MP2 TOTAL ENERGY"] = Esospimp2;

        Process::environment.globals["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
        Process::environment.globals["SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
        Process::environment.globals["SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
        Process::environment.globals["SCSN-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;
        Process::environment.globals["SCS-MP2-VDW CORRELATION ENERGY"] = Escsmp2vdw - Escf;
        Process::environment.globals["SOS-PI-MP2 CORRELATION ENERGY"] = Esospimp2 - Escf;

        Process::environment.globals["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
        Process::environment.globals["MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA+Emp2BB;

        timer_on("T2(2)");
	t2_2nd_sc();
        timer_off("T2(2)");
        timer_on("MP3 Energy");
	mp3_energy();
        timer_off("MP3 Energy");
	Emp3L=Emp3;
        EcorrL=Emp3L-Escf;
	Emp3L_old=Emp3;
        if (ip_poles == "TRUE") omp3_ip_poles();
	
	psi::fprintf(outfile,"\n"); 
	psi::fprintf(outfile,"\tComputing MP3 energy using SCF MOs (Canonical MP3)... \n"); 
	psi::fprintf(outfile,"\t============================================================================== \n");
	psi::fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	psi::fprintf(outfile,"\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	psi::fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	psi::fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp3AA);
	psi::fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp3AB);
	psi::fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp3BB);
	psi::fprintf(outfile,"\tMP2.5 Correlation Energy (a.u.)    : %20.14f\n", (Emp2 - Escf) + 0.5 * (Emp3-Emp2));
	psi::fprintf(outfile,"\tMP2.5 Total Energy (a.u.)          : %20.14f\n", 0.5 * (Emp3+Emp2));
	psi::fprintf(outfile,"\tSCS-MP3 Total Energy (a.u.)        : %20.14f\n", Escsmp3);
	psi::fprintf(outfile,"\tSOS-MP3 Total Energy (a.u.)        : %20.14f\n", Esosmp3);
	psi::fprintf(outfile,"\tSCSN-MP3 Total Energy (a.u.)       : %20.14f\n", Escsnmp3);
	psi::fprintf(outfile,"\tSCS-MP3-VDW Total Energy (a.u.)    : %20.14f\n", Escsmp3vdw);
	psi::fprintf(outfile,"\tSOS-PI-MP3 Total Energy (a.u.)     : %20.14f\n", Esospimp3);
	psi::fprintf(outfile,"\t3rd Order Energy (a.u.)            : %20.14f\n", Emp3-Emp2);
	psi::fprintf(outfile,"\tMP3 Correlation Energy (a.u.)      : %20.14f\n", Ecorr);
	psi::fprintf(outfile,"\tMP3 Total Energy (a.u.)            : %20.14f\n", Emp3);
	psi::fprintf(outfile,"\t============================================================================== \n");
	psi::fprintf(outfile,"\n"); 
	fflush(outfile);
	Process::environment.globals["MP3 TOTAL ENERGY"] = Emp3;
	Process::environment.globals["SCS-MP3 TOTAL ENERGY"] = Escsmp3;
	Process::environment.globals["SOS-MP3 TOTAL ENERGY"] = Esosmp3;
	Process::environment.globals["SCSN-MP3 TOTAL ENERGY"] = Escsnmp3;
	Process::environment.globals["SCS-MP3-VDW TOTAL ENERGY"] = Escsmp3vdw;
	Process::environment.globals["SOS-PI-MP3 TOTAL ENERGY"] = Esospimp3;

	Process::environment.globals["MP2.5 CORRELATION ENERGY"] = (Emp2 - Escf) + 0.5 * (Emp3-Emp2);
	Process::environment.globals["MP2.5 TOTAL ENERGY"] = 0.5 * (Emp3+Emp2);
	Process::environment.globals["MP3 CORRELATION ENERGY"] = Emp3 - Escf;
	Process::environment.globals["SCS-MP3 CORRELATION ENERGY"] = Escsmp3 - Escf;
	Process::environment.globals["SOS-MP3 CORRELATION ENERGY"] = Esosmp3 - Escf;
	Process::environment.globals["SCSN-MP3 CORRELATION ENERGY"] = Escsnmp3 - Escf;
	Process::environment.globals["SCS-MP3-VDW CORRELATION ENERGY"] = Escsmp3vdw - Escf;
	Process::environment.globals["SOS-PI-MP3 CORRELATION ENERGY"] = Esospimp3 - Escf;

	omp3_response_pdms();
	gfock();
	idp();
	mograd();
        occ_iterations();
	
        if (rms_wog <= tol_grad && fabs(DE) >= tol_Eod) {
           orbs_already_opt = 1;
	   if (conver == 1) psi::fprintf(outfile,"\n\tOrbitals are optimized now.\n");
	   else if (conver == 0) { 
                    psi::fprintf(outfile,"\n\tMAX MOGRAD did NOT converged, but RMS MOGRAD converged!!!\n");
	            psi::fprintf(outfile,"\tI will consider the present orbitals as optimized.\n");
           }
	   psi::fprintf(outfile,"\tSwitching to the standard MP3 computation after semicanonicalization of the MOs... \n");
	   fflush(outfile);
	   semi_canonic();
	   if (reference_ == "RESTRICTED") trans_ints_rhf();  
	   else if (reference_ == "UNRESTRICTED") trans_ints_uhf();  
	   omp3_t2_1st_sc();
	   t2_2nd_sc();
           conver = 1;
           if (dertype == "FIRST") {
               omp3_response_pdms();
	       gfock();
           }
        }     

  if (conver == 1) {
        ref_energy();
	omp3_mp2_energy();
	mp3_energy();
        if (orbs_already_opt == 1) Emp3L = Emp3;

        if (ip_poles == "TRUE") {
	   if (orbs_already_sc == 0) {
               semi_canonic();
	       if (reference_ == "RESTRICTED") trans_ints_rhf();  
	       else if (reference_ == "UNRESTRICTED") trans_ints_uhf();  
	       omp3_t2_1st_sc();
	       t2_2nd_sc();
           }
           omp3_ip_poles();
        }

        // EKT
        if (ekt_ip_ == "TRUE" || ekt_ea_ == "TRUE") { 
            if (orbs_already_sc == 1) {
	        omp3_response_pdms();
	        gfock();
            }
            gfock_diag();
            if (ekt_ip_ == "TRUE") ekt_ip();
            if (ekt_ea_ == "TRUE") ekt_ea();
        }

        psi::fprintf(outfile,"\n"); 
	psi::fprintf(outfile,"\tComputing MP2 energy using optimized MOs... \n"); 
	psi::fprintf(outfile,"\t============================================================================== \n");
	psi::fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	psi::fprintf(outfile,"\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	psi::fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	psi::fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
	psi::fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
	psi::fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
	psi::fprintf(outfile,"\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA+Escsmp2BB);
	psi::fprintf(outfile,"\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
	psi::fprintf(outfile,"\tSCS-MP2 Total Energy (a.u.)        : %20.14f\n", Escsmp2);
	psi::fprintf(outfile,"\tSOS-MP2 Total Energy (a.u.)        : %20.14f\n", Esosmp2);
	psi::fprintf(outfile,"\tSCSN-MP2 Total Energy (a.u.)       : %20.14f\n", Escsnmp2);
	psi::fprintf(outfile,"\tSCS-MP2-VDW Total Energy (a.u.)    : %20.14f\n", Escsmp2vdw);
	psi::fprintf(outfile,"\tSOS-PI-MP2 Total Energy (a.u.)     : %20.14f\n", Esospimp2);
	psi::fprintf(outfile,"\tMP2 Correlation Energy (a.u.)      : %20.14f\n", Ecorr);
	psi::fprintf(outfile,"\tMP2 Total Energy (a.u.)            : %20.14f\n", Emp2);
	psi::fprintf(outfile,"\t============================================================================== \n");
	psi::fprintf(outfile,"\n"); 
	fflush(outfile);

	psi::fprintf(outfile,"\n"); 
	psi::fprintf(outfile,"\tComputing MP3 energy using optimized MOs... \n"); 
	psi::fprintf(outfile,"\t============================================================================== \n");
	psi::fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	psi::fprintf(outfile,"\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	psi::fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	psi::fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp3AA);
	psi::fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp3AB);
	psi::fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp3BB);
	psi::fprintf(outfile,"\tMP2.5 Correlation Energy (a.u.)    : %20.14f\n", (Emp2 - Escf) + 0.5 * (Emp3-Emp2));
	psi::fprintf(outfile,"\tMP2.5 Total Energy (a.u.)          : %20.14f\n", 0.5 * (Emp3+Emp2));
	psi::fprintf(outfile,"\tSCS-MP3 Total Energy (a.u.)        : %20.14f\n", Escsmp3);
	psi::fprintf(outfile,"\tSOS-MP3 Total Energy (a.u.)        : %20.14f\n", Esosmp3);
	psi::fprintf(outfile,"\tSCSN-MP3 Total Energy (a.u.)       : %20.14f\n", Escsnmp3);
	psi::fprintf(outfile,"\tSCS-MP3-VDW Total Energy (a.u.)    : %20.14f\n", Escsmp3vdw);
	psi::fprintf(outfile,"\tSOS-PI-MP3 Total Energy (a.u.)     : %20.14f\n", Esospimp3);
	psi::fprintf(outfile,"\t3rd Order Energy (a.u.)            : %20.14f\n", Emp3-Emp2);
	psi::fprintf(outfile,"\tMP3 Correlation Energy (a.u.)      : %20.14f\n", Ecorr);
	psi::fprintf(outfile,"\tMP3 Total Energy (a.u.)            : %20.14f\n", Emp3);
	psi::fprintf(outfile,"\t============================================================================== \n");
	psi::fprintf(outfile,"\n"); 
	fflush(outfile);


	psi::fprintf(outfile,"\n");
	psi::fprintf(outfile,"\t============================================================================== \n");
	psi::fprintf(outfile,"\t================ OMP3 FINAL RESULTS ========================================== \n");
	psi::fprintf(outfile,"\t============================================================================== \n");
	psi::fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	psi::fprintf(outfile,"\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	psi::fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	psi::fprintf(outfile,"\tSCS-OMP3 Total Energy (a.u.)       : %20.14f\n", Escsmp3);
	psi::fprintf(outfile,"\tSOS-OMP3 Total Energy (a.u.)       : %20.14f\n", Esosmp3);
	psi::fprintf(outfile,"\tSCSN-OMP3 Total Energy (a.u.)      : %20.14f\n", Escsnmp3);
	psi::fprintf(outfile,"\tSCS-OMP3-VDW Total Energy (a.u.    : %20.14f\n", Escsmp3vdw);
	psi::fprintf(outfile,"\tSOS-PI-OMP3 Total Energy (a.u.)    : %20.14f\n", Esospimp3);
	psi::fprintf(outfile,"\tOMP3 Correlation Energy (a.u.)     : %20.14f\n", Emp3L-Escf);
	psi::fprintf(outfile,"\tEomp3 - Eref (a.u.)                : %20.14f\n", Emp3L-Eref);
	psi::fprintf(outfile,"\tOMP3 Total Energy (a.u.)           : %20.14f\n", Emp3L);
	psi::fprintf(outfile,"\t============================================================================== \n");
	psi::fprintf(outfile,"\n");
	fflush(outfile);
	
	// Set the global variables with the energies
	Process::environment.globals["OMP3 TOTAL ENERGY"] = Emp3L;
	Process::environment.globals["SCS-OMP3 TOTAL ENERGY"] =  Escsmp3;
	Process::environment.globals["SOS-OMP3 TOTAL ENERGY"] =  Esosmp3;
	Process::environment.globals["SCSN-OMP3 TOTAL ENERGY"] = Escsnmp3;
	Process::environment.globals["SCS-OMP3-VDW TOTAL ENERGY"] = Escsmp3vdw;
	Process::environment.globals["SOS-PI-OMP3 TOTAL ENERGY"] = Esospimp3;
	Process::environment.globals["CURRENT ENERGY"] = Emp3L;
	Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
	Process::environment.globals["CURRENT CORRELATION ENERGY"] = Emp3L-Escf;

	Process::environment.globals["OMP3 CORRELATION ENERGY"] = Emp3L - Escf;
	Process::environment.globals["SCS-OMP3 CORRELATION ENERGY"] =  Escsmp3 - Escf;
	Process::environment.globals["SOS-OMP3 CORRELATION ENERGY"] =  Esosmp3 - Escf;
	Process::environment.globals["SCSN-OMP3 CORRELATION ENERGY"] = Escsnmp3 - Escf;
	Process::environment.globals["SCS-OMP3-VDW CORRELATION ENERGY"] = Escsmp3vdw - Escf;
	Process::environment.globals["SOS-PI-OMP3 CORRELATION ENERGY"] = Esospimp3 - Escf;

        // if scs on	
	if (do_scs == "TRUE") {
	    if (scs_type_ == "SCS") {
	       Process::environment.globals["CURRENT ENERGY"] = Escsmp3;
	       Process::environment.globals["CURRENT CORRELATION ENERGY"] = Escsmp3 - Escf;
            }

	    else if (scs_type_ == "SCSN") {
	       Process::environment.globals["CURRENT ENERGY"] = Escsnmp3;
	       Process::environment.globals["CURRENT CORRELATION ENERGY"] = Escsnmp3 - Escf;
            }

	    else if (scs_type_ == "SCSVDW") {
	       Process::environment.globals["CURRENT ENERGY"] = Escsmp3vdw;
	       Process::environment.globals["CURRENT CORRELATION ENERGY"] = Escsmp3vdw - Escf;
            }
	}
    
        // else if sos on	
	else if (do_sos == "TRUE") {
	     if (sos_type_ == "SOS") {
	         Process::environment.globals["CURRENT ENERGY"] = Esosmp3;
	         Process::environment.globals["CURRENT CORRELATION ENERGY"] = Esosmp3 - Escf;
             }

	     else if (sos_type_ == "SOSPI") {
	             Process::environment.globals["CURRENT ENERGY"] = Esospimp3;
	             Process::environment.globals["CURRENT CORRELATION ENERGY"] = Esospimp3 - Escf;
             }
	}

	if (natorb == "TRUE") nbo();
	if (occ_orb_energy == "TRUE") semi_canonic(); 

        // Compute Analytic Gradients
        if (dertype == "FIRST") {
            time4grad = 1;
	    psi::fprintf(outfile,"\tAnalytic gradient computation is starting...\n");
	    fflush(outfile);
            coord_grad();
	    psi::fprintf(outfile,"\tNecessary information has been sent to DERIV, which will take care of the rest.\n");
	    fflush(outfile);
        }

  }// end if (conver == 1)
  */
}// end omp3_manager 

//======================================================================
//             MP3 Manager
//======================================================================             
void DFOCC::mp3_manager()
{
        /*
        time4grad = 0;// means i will not compute the gradient
        timer_on("trans_ints");
	if (reference_ == "RESTRICTED") trans_ints_rhf();  
	else if (reference_ == "UNRESTRICTED") trans_ints_uhf();  
        timer_off("trans_ints");
        Eref = Escf;
        timer_on("T2(1)");
	omp3_t2_1st_sc();
        timer_off("T2(1)");
        timer_on("MP2 Energy");
	omp3_mp2_energy();
        timer_off("MP2 Energy");

	psi::fprintf(outfile,"\n"); 
	psi::fprintf(outfile,"\tComputing MP2 energy using SCF MOs (Canonical MP2)... \n"); 
	psi::fprintf(outfile,"\t============================================================================== \n");
	psi::fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	psi::fprintf(outfile,"\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	psi::fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	psi::fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
	psi::fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
	psi::fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
	psi::fprintf(outfile,"\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA+Escsmp2BB);
	psi::fprintf(outfile,"\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
	psi::fprintf(outfile,"\tSCS-MP2 Total Energy (a.u.)        : %20.14f\n", Escsmp2);
	psi::fprintf(outfile,"\tSOS-MP2 Total Energy (a.u.)        : %20.14f\n", Esosmp2);
	psi::fprintf(outfile,"\tSCSN-MP2 Total Energy (a.u.)       : %20.14f\n", Escsnmp2);
	psi::fprintf(outfile,"\tSCS-MP2-VDW Total Energy (a.u.)    : %20.14f\n", Escsmp2vdw);
	psi::fprintf(outfile,"\tSOS-PI-MP2 Total Energy (a.u.)     : %20.14f\n", Esospimp2);
	psi::fprintf(outfile,"\tMP2 Correlation Energy (a.u.)      : %20.14f\n", Ecorr);
	psi::fprintf(outfile,"\tMP2 Total Energy (a.u.)            : %20.14f\n", Emp2);
	psi::fprintf(outfile,"\t============================================================================== \n");
	psi::fprintf(outfile,"\n"); 
	fflush(outfile);
	Process::environment.globals["MP2 TOTAL ENERGY"] = Emp2;
	Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = Escsmp2;
	Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
	Process::environment.globals["SCSN-MP2 TOTAL ENERGY"] = Escsnmp2;
	Process::environment.globals["SCS-MP2-VDW TOTAL ENERGY"] = Escsmp2vdw;
	Process::environment.globals["SOS-PI-MP2 TOTAL ENERGY"] = Esospimp2;

        Process::environment.globals["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
        Process::environment.globals["SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
        Process::environment.globals["SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
        Process::environment.globals["SCSN-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;
        Process::environment.globals["SCS-MP2-VDW CORRELATION ENERGY"] = Escsmp2vdw - Escf;
        Process::environment.globals["SOS-PI-MP2 CORRELATION ENERGY"] = Esospimp2 - Escf;

        Process::environment.globals["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
        Process::environment.globals["MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA+Emp2BB;

        timer_on("T2(2)");
	t2_2nd_sc();
        timer_off("T2(2)");
        timer_on("MP3 Energy");
	mp3_energy();
        timer_off("MP3 Energy");
	Emp3L=Emp3;
        EcorrL=Emp3L-Escf;
	Emp3L_old=Emp3;
        if (ip_poles == "TRUE") omp3_ip_poles();
	
	psi::fprintf(outfile,"\n"); 
	psi::fprintf(outfile,"\tComputing MP3 energy using SCF MOs (Canonical MP3)... \n"); 
	psi::fprintf(outfile,"\t============================================================================== \n");
	psi::fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	psi::fprintf(outfile,"\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	psi::fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	psi::fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp3AA);
	psi::fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp3AB);
	psi::fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp3BB);
	psi::fprintf(outfile,"\tMP2.5 Correlation Energy (a.u.)    : %20.14f\n", (Emp2 - Escf) + 0.5 * (Emp3-Emp2));
	psi::fprintf(outfile,"\tMP2.5 Total Energy (a.u.)          : %20.14f\n", 0.5 * (Emp3+Emp2));
	psi::fprintf(outfile,"\tSCS-MP3 Total Energy (a.u.)        : %20.14f\n", Escsmp3);
	psi::fprintf(outfile,"\tSOS-MP3 Total Energy (a.u.)        : %20.14f\n", Esosmp3);
	psi::fprintf(outfile,"\tSCSN-MP3 Total Energy (a.u.)       : %20.14f\n", Escsnmp3);
	psi::fprintf(outfile,"\tSCS-MP3-VDW Total Energy (a.u.)    : %20.14f\n", Escsmp3vdw);
	psi::fprintf(outfile,"\tSOS-PI-MP3 Total Energy (a.u.)     : %20.14f\n", Esospimp3);
	psi::fprintf(outfile,"\t3rd Order Energy (a.u.)            : %20.14f\n", Emp3-Emp2);
	psi::fprintf(outfile,"\tMP3 Correlation Energy (a.u.)      : %20.14f\n", Ecorr);
	psi::fprintf(outfile,"\tMP3 Total Energy (a.u.)            : %20.14f\n", Emp3);
	psi::fprintf(outfile,"\t============================================================================== \n");
	psi::fprintf(outfile,"\n"); 
	fflush(outfile);
	Process::environment.globals["CURRENT ENERGY"] = Emp3;
	Process::environment.globals["CURRENT CORRELATION ENERGY"] = Emp3 - Escf;
	Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;

	Process::environment.globals["MP3 TOTAL ENERGY"] = Emp3;
	Process::environment.globals["SCS-MP3 TOTAL ENERGY"] = Escsmp3;
	Process::environment.globals["SOS-MP3 TOTAL ENERGY"] = Esosmp3;
	Process::environment.globals["SCSN-MP3 TOTAL ENERGY"] = Escsnmp3;
	Process::environment.globals["SCS-MP3-VDW TOTAL ENERGY"] = Escsmp3vdw;
	Process::environment.globals["SOS-PI-MP3 TOTAL ENERGY"] = Esospimp3;

	Process::environment.globals["MP2.5 CORRELATION ENERGY"] = (Emp2 - Escf) + 0.5 * (Emp3-Emp2);
	Process::environment.globals["MP2.5 TOTAL ENERGY"] = 0.5 * (Emp3+Emp2);
	Process::environment.globals["MP3 CORRELATION ENERGY"] = Emp3 - Escf;
	Process::environment.globals["SCS-MP3 CORRELATION ENERGY"] = Escsmp3 - Escf;
	Process::environment.globals["SOS-MP3 CORRELATION ENERGY"] = Esosmp3 - Escf;
	Process::environment.globals["SCSN-MP3 CORRELATION ENERGY"] = Escsnmp3 - Escf;
	Process::environment.globals["SCS-MP3-VDW CORRELATION ENERGY"] = Escsmp3vdw - Escf;
	Process::environment.globals["SOS-PI-MP3 CORRELATION ENERGY"] = Esospimp3 - Escf;

        // if scs on	
	if (do_scs == "TRUE") {
	    if (scs_type_ == "SCS") {
	       Process::environment.globals["CURRENT ENERGY"] = Escsmp3;
	       Process::environment.globals["CURRENT CORRELATION ENERGY"] = Escsmp3 - Escf;
            }

	    else if (scs_type_ == "SCSN") {
	       Process::environment.globals["CURRENT ENERGY"] = Escsnmp3;
	       Process::environment.globals["CURRENT CORRELATION ENERGY"] = Escsnmp3 - Escf;
            }

	    else if (scs_type_ == "SCSVDW") {
	       Process::environment.globals["CURRENT ENERGY"] = Escsmp3vdw;
	       Process::environment.globals["CURRENT CORRELATION ENERGY"] = Escsmp3vdw - Escf;
            }
	}
    
        // else if sos on	
	else if (do_sos == "TRUE") {
	     if (sos_type_ == "SOS") {
	         Process::environment.globals["CURRENT ENERGY"] = Esosmp3;
	         Process::environment.globals["CURRENT CORRELATION ENERGY"] = Esosmp3 - Escf;
             }

	     else if (sos_type_ == "SOSPI") {
	             Process::environment.globals["CURRENT ENERGY"] = Esospimp3;
	             Process::environment.globals["CURRENT CORRELATION ENERGY"] = Esospimp3 - Escf;
             }
	}

        // Compute Analytic Gradients
        if (dertype == "FIRST" || ekt_ip_ == "TRUE" || ekt_ea_ == "TRUE") {
            time4grad = 1;
	    psi::fprintf(outfile,"\tAnalytic gradient computation is starting...\n");
            psi::fprintf(outfile,"\tComputing response density matrices...\n");
            fflush(outfile);
	    omp3_response_pdms();
            psi::fprintf(outfile,"\tComputing off-diagonal blocks of GFM...\n");
            fflush(outfile);
	    gfock();
            psi::fprintf(outfile,"\tForming independent-pairs...\n");
            fflush(outfile);
	    idp2();
            psi::fprintf(outfile,"\tComputing orbital gradient...\n");
            fflush(outfile);
	    mograd();
            coord_grad();

            if (ekt_ip_ == "TRUE" && ekt_ea_ == "TRUE") {
                ekt_ip();
                ekt_ea();
            }

            else if (ekt_ip_ == "TRUE" && ekt_ea_ == "FALSE") {
                ekt_ip();
            }

            else if (ekt_ip_ == "FALSE" && ekt_ea_ == "TRUE") {
                ekt_ea();
            }

            else if (ekt_ip_ == "FALSE" && ekt_ea_ == "FALSE") {
	        psi::fprintf(outfile,"\tNecessary information has been sent to DERIV, which will take care of the rest.\n");
	        fflush(outfile);
            }

        }
        */

}// end mp3_manager 


//======================================================================
//             OCEPA Manager
//======================================================================             
void DFOCC::ocepa_manager()
{
        /*
	mo_optimized = 0;// means MOs are not optimized yet.
	orbs_already_opt = 0;
	orbs_already_sc = 0;
        time4grad = 0;// means i will not compute the gradient
        timer_on("trans_ints");
	if (reference_ == "RESTRICTED") trans_ints_rhf();  
	else if (reference_ == "UNRESTRICTED") trans_ints_uhf();  
        timer_off("trans_ints");
        timer_on("REF Energy");
	ref_energy();
        timer_off("REF Energy");
        timer_on("T2(1)");
	ocepa_t2_1st_sc();
        timer_off("T2(1)");
        timer_on("MP2 Energy");
	ocepa_mp2_energy();
        timer_off("MP2 Energy");
        Ecepa = Emp2;
	EcepaL = Ecepa;
        EcorrL = Ecorr;
	EcepaL_old = Ecepa;

	psi::fprintf(outfile,"\n"); 
	psi::fprintf(outfile,"\tComputing MP2 energy using SCF MOs (Canonical MP2)... \n"); 
	psi::fprintf(outfile,"\t============================================================================== \n");
	psi::fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	psi::fprintf(outfile,"\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	psi::fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	psi::fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
	psi::fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
	psi::fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
	psi::fprintf(outfile,"\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA+Escsmp2BB);
	psi::fprintf(outfile,"\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
	psi::fprintf(outfile,"\tSCS-MP2 Total Energy (a.u.)        : %20.14f\n", Escsmp2);
	psi::fprintf(outfile,"\tSOS-MP2 Total Energy (a.u.)        : %20.14f\n", Esosmp2);
	psi::fprintf(outfile,"\tSCSN-MP2 Total Energy (a.u.)       : %20.14f\n", Escsnmp2);
	psi::fprintf(outfile,"\tSCS-MP2-VDW Total Energy (a.u.)    : %20.14f\n", Escsmp2vdw);
	psi::fprintf(outfile,"\tSOS-PI-MP2 Total Energy (a.u.)     : %20.14f\n", Esospimp2);
	psi::fprintf(outfile,"\tMP2 Correlation Energy (a.u.)      : %20.14f\n", Ecorr);
	psi::fprintf(outfile,"\tMP2 Total Energy (a.u.)            : %20.14f\n", Emp2);
	psi::fprintf(outfile,"\t============================================================================== \n");
	fflush(outfile);
	Process::environment.globals["MP2 TOTAL ENERGY"] = Emp2;
	Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = Escsmp2;
	Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
	Process::environment.globals["SCSN-MP2 TOTAL ENERGY"] = Escsnmp2;
	Process::environment.globals["SCS-MP2-VDW TOTAL ENERGY"] = Escsmp2vdw;
	Process::environment.globals["SOS-PI-MP2 TOTAL ENERGY"] = Esospimp2;

	Process::environment.globals["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
	Process::environment.globals["SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
	Process::environment.globals["SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
	Process::environment.globals["SCSN-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;
	Process::environment.globals["SCS-MP2-VDW CORRELATION ENERGY"] = Escsmp2vdw - Escf;
	Process::environment.globals["SOS-PI-MP2 CORRELATION ENERGY"] = Esospimp2 - Escf;

        Process::environment.globals["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
        Process::environment.globals["MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA+Emp2BB;

	ocepa_response_pdms();
	gfock();
	idp();
	mograd();
        if (rms_wog > tol_grad) occ_iterations();
        else {
           orbs_already_opt = 1;
	   psi::fprintf(outfile,"\n\tOrbitals are already optimized, switching to the canonical CEPA computation... \n");
	   fflush(outfile);
           cepa_iterations();
        }
	
        if (rms_wog <= tol_grad && fabs(DE) >= tol_Eod) {
           orbs_already_opt = 1;
	   if (conver == 1) psi::fprintf(outfile,"\n\tOrbitals are optimized now.\n");
	   else if (conver == 0) { 
                    psi::fprintf(outfile,"\n\tMAX MOGRAD did NOT converged, but RMS MOGRAD converged!!!\n");
	            psi::fprintf(outfile,"\tI will consider the present orbitals as optimized.\n");
           }
	   psi::fprintf(outfile,"\tSwitching to the standard CEPA computation... \n");
	   fflush(outfile);
           ref_energy();
	   cepa_energy();
           Ecepa_old = EcepaL;
           cepa_iterations();
           if (dertype == "FIRST") {
               ocepa_response_pdms();
	       gfock();
           }
        } 

  if (conver == 1) {
        ref_energy();
	cepa_energy();
        if (orbs_already_opt == 1) EcepaL = Ecepa;

        // EKT
        if (ekt_ip_ == "TRUE" || ekt_ea_ == "TRUE") { 
            if (orbs_already_opt == 1) {
	        ocepa_response_pdms();
	        gfock();
            }
            gfock_diag();
            if (ekt_ip_ == "TRUE") ekt_ip();
            if (ekt_ea_ == "TRUE") ekt_ea();
        }

	psi::fprintf(outfile,"\n");
	psi::fprintf(outfile,"\t============================================================================== \n");
	psi::fprintf(outfile,"\t================ OCEPA FINAL RESULTS ========================================= \n");
	psi::fprintf(outfile,"\t============================================================================== \n");
	psi::fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	psi::fprintf(outfile,"\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	psi::fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	psi::fprintf(outfile,"\tSCS-OCEPA(0) Total Energy (a.u.)   : %20.14f\n", Escscepa);
	psi::fprintf(outfile,"\tSOS-OCEPA(0) Total Energy (a.u.)   : %20.14f\n", Esoscepa);
	psi::fprintf(outfile,"\tOCEPA(0) Correlation Energy (a.u.) : %20.14f\n", EcepaL-Escf);
	psi::fprintf(outfile,"\tEocepa - Eref (a.u.)               : %20.14f\n", EcepaL-Eref);
	psi::fprintf(outfile,"\tOCEPA(0) Total Energy (a.u.)       : %20.14f\n", EcepaL);
	psi::fprintf(outfile,"\t============================================================================== \n");
	psi::fprintf(outfile,"\n");
	fflush(outfile);
	
	// Set the global variables with the energies
	Process::environment.globals["OCEPA(0) TOTAL ENERGY"] = EcepaL;
	Process::environment.globals["SCS-OCEPA(0) TOTAL ENERGY"] =  Escscepa;
	Process::environment.globals["SOS-OCEPA(0) TOTAL ENERGY"] =  Esoscepa;
	Process::environment.globals["CURRENT ENERGY"] = EcepaL;
	Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
	Process::environment.globals["CURRENT CORRELATION ENERGY"] = EcepaL-Escf;

	Process::environment.globals["OCEPA(0) CORRELATION ENERGY"] = EcepaL - Escf;
	Process::environment.globals["SCS-OCEPA(0) CORRELATION ENERGY"] =  Escscepa - Escf;
	Process::environment.globals["SOS-OCEPA(0) CORRELATION ENERGY"] =  Esoscepa - Escf;

        // if scs on	
	if (do_scs == "TRUE") {
	    Process::environment.globals["CURRENT ENERGY"] = Escscepa;
	    Process::environment.globals["CURRENT CORRELATION ENERGY"] = Escscepa - Escf;

	}
    
        // else if sos on	
	else if (do_sos == "TRUE") {
	         Process::environment.globals["CURRENT ENERGY"] = Esoscepa;
	         Process::environment.globals["CURRENT CORRELATION ENERGY"] = Esoscepa - Escf;
	}

	if (natorb == "TRUE") nbo();
	if (occ_orb_energy == "TRUE") semi_canonic(); 
	
        // Compute Analytic Gradients
        if (dertype == "FIRST") {
            time4grad = 1;
	    psi::fprintf(outfile,"\tAnalytic gradient computation is starting...\n");
	    fflush(outfile);
            coord_grad();
	    psi::fprintf(outfile,"\tNecessary information has been sent to DERIV, which will take care of the rest.\n");
	    fflush(outfile);
        }

  }// end if (conver == 1)
  */
}// end ocepa_manager 


//======================================================================
//             CEPA Manager
//======================================================================             
void DFOCC::cepa_manager()
{
        /*
        timer_on("trans_ints");
	if (reference_ == "RESTRICTED") trans_ints_rhf();  
	else if (reference_ == "UNRESTRICTED") trans_ints_uhf();  
        timer_off("trans_ints");
        Eref = Escf;
        timer_on("T2(1)");
	ocepa_t2_1st_sc();
        timer_off("T2(1)");
        timer_on("MP2 Energy");
	ocepa_mp2_energy();
        timer_off("MP2 Energy");
        Ecepa = Emp2;
	Ecepa_old = Emp2;

	psi::fprintf(outfile,"\n"); 
	psi::fprintf(outfile,"\tComputing MP2 energy using SCF MOs (Canonical MP2)... \n"); 
	psi::fprintf(outfile,"\t============================================================================== \n");
	psi::fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	psi::fprintf(outfile,"\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	psi::fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	psi::fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
	psi::fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
	psi::fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
	psi::fprintf(outfile,"\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA+Escsmp2BB);
	psi::fprintf(outfile,"\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
	psi::fprintf(outfile,"\tSCS-MP2 Total Energy (a.u.)        : %20.14f\n", Escsmp2);
	psi::fprintf(outfile,"\tSOS-MP2 Total Energy (a.u.)        : %20.14f\n", Esosmp2);
	psi::fprintf(outfile,"\tSCSN-MP2 Total Energy (a.u.)       : %20.14f\n", Escsnmp2);
	psi::fprintf(outfile,"\tSCS-MP2-VDW Total Energy (a.u.)    : %20.14f\n", Escsmp2vdw);
	psi::fprintf(outfile,"\tSOS-PI-MP2 Total Energy (a.u.)     : %20.14f\n", Esospimp2);
	psi::fprintf(outfile,"\tMP2 Correlation Energy (a.u.)      : %20.14f\n", Ecorr);
	psi::fprintf(outfile,"\tMP2 Total Energy (a.u.)            : %20.14f\n", Emp2);
	psi::fprintf(outfile,"\t============================================================================== \n");
	fflush(outfile);
	Process::environment.globals["MP2 TOTAL ENERGY"] = Emp2;
	Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = Escsmp2;
	Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
	Process::environment.globals["SCSN-MP2 TOTAL ENERGY"] = Escsnmp2;
	Process::environment.globals["SCS-MP2-VDW TOTAL ENERGY"] = Escsmp2vdw;
	Process::environment.globals["SOS-PI-MP2 TOTAL ENERGY"] = Esospimp2;

	Process::environment.globals["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
	Process::environment.globals["SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
	Process::environment.globals["SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
	Process::environment.globals["SCSN-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;
	Process::environment.globals["SCS-MP2-VDW CORRELATION ENERGY"] = Escsmp2vdw - Escf;
	Process::environment.globals["SOS-PI-MP2 CORRELATION ENERGY"] = Esospimp2 - Escf;

        Process::environment.globals["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
        Process::environment.globals["MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA+Emp2BB;

        // Perform CEPA iterations
        cepa_iterations();

	psi::fprintf(outfile,"\n");
	psi::fprintf(outfile,"\t============================================================================== \n");
	psi::fprintf(outfile,"\t================ CEPA FINAL RESULTS ========================================== \n");
	psi::fprintf(outfile,"\t============================================================================== \n");
	psi::fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	psi::fprintf(outfile,"\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	psi::fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	psi::fprintf(outfile,"\tCEPA(0) Correlation Energy (a.u.)  : %20.14f\n", Ecorr);
	psi::fprintf(outfile,"\tCEPA(0) Total Energy (a.u.)        : %20.14f\n", Ecepa);
	psi::fprintf(outfile,"\t============================================================================== \n");
	psi::fprintf(outfile,"\n");
	fflush(outfile);
	
	// Set the global variables with the energies
	Process::environment.globals["CEPA(0) TOTAL ENERGY"] = Ecepa;
	Process::environment.globals["CURRENT ENERGY"] = Ecepa;
	Process::environment.globals["CURRENT REFERENCE ENERGY"] = Eref;
	Process::environment.globals["CURRENT CORRELATION ENERGY"] = Ecorr;
        //EcepaL = Ecepa;
        
        // Compute Analytic Gradients
        if (dertype == "FIRST" || ekt_ip_ == "TRUE" || ekt_ea_ == "TRUE") {
            time4grad = 1;
	    psi::fprintf(outfile,"\tAnalytic gradient computation is starting...\n");
            psi::fprintf(outfile,"\tComputing response density matrices...\n");
            fflush(outfile);
	    ocepa_response_pdms();
            psi::fprintf(outfile,"\tComputing off-diagonal blocks of GFM...\n");
            fflush(outfile);
	    gfock();
            psi::fprintf(outfile,"\tForming independent-pairs...\n");
            fflush(outfile);
	    idp2();
            psi::fprintf(outfile,"\tComputing orbital gradient...\n");
            fflush(outfile);
	    mograd();
            coord_grad();

            if (ekt_ip_ == "TRUE" && ekt_ea_ == "TRUE") {
                ekt_ip();
                ekt_ea();
            }

            else if (ekt_ip_ == "TRUE" && ekt_ea_ == "FALSE") {
                ekt_ip();
            }

            else if (ekt_ip_ == "FALSE" && ekt_ea_ == "TRUE") {
                ekt_ea();
            }

            else if (ekt_ip_ == "FALSE" && ekt_ea_ == "FALSE") {
	        psi::fprintf(outfile,"\tNecessary information has been sent to DERIV, which will take care of the rest.\n");
	        fflush(outfile);
            }

        }
        */
}// end cepa_manager 


//======================================================================
//             OMP2.5 Manager
//======================================================================             
void DFOCC::omp2_5_manager()
{
        /*
	mo_optimized = 0;
	orbs_already_opt = 0;
	orbs_already_sc = 0;
        time4grad = 0;// means i will not compute the gradient
        timer_on("trans_ints");
	if (reference_ == "RESTRICTED") trans_ints_rhf();  
	else if (reference_ == "UNRESTRICTED") trans_ints_uhf();  
        timer_off("trans_ints");
        timer_on("REF Energy");
	ref_energy();
        timer_off("REF Energy");
        timer_on("T2(1)");
	omp3_t2_1st_sc();
        timer_off("T2(1)");
        timer_on("MP2 Energy");
	omp3_mp2_energy();
        timer_off("MP2 Energy");

	psi::fprintf(outfile,"\n"); 
	psi::fprintf(outfile,"\tComputing MP2 energy using SCF MOs (Canonical MP2)... \n"); 
	psi::fprintf(outfile,"\t============================================================================== \n");
	psi::fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	psi::fprintf(outfile,"\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	psi::fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	psi::fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
	psi::fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
	psi::fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
	psi::fprintf(outfile,"\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA+Escsmp2BB);
	psi::fprintf(outfile,"\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
	psi::fprintf(outfile,"\tSCS-MP2 Total Energy (a.u.)        : %20.14f\n", Escsmp2);
	psi::fprintf(outfile,"\tSOS-MP2 Total Energy (a.u.)        : %20.14f\n", Esosmp2);
	psi::fprintf(outfile,"\tSCSN-MP2 Total Energy (a.u.)       : %20.14f\n", Escsnmp2);
	psi::fprintf(outfile,"\tSCS-MP2-VDW Total Energy (a.u.)    : %20.14f\n", Escsmp2vdw);
	psi::fprintf(outfile,"\tSOS-PI-MP2 Total Energy (a.u.)     : %20.14f\n", Esospimp2);
	psi::fprintf(outfile,"\tMP2 Correlation Energy (a.u.)      : %20.14f\n", Ecorr);
	psi::fprintf(outfile,"\tMP2 Total Energy (a.u.)            : %20.14f\n", Emp2);
	psi::fprintf(outfile,"\t============================================================================== \n");
	psi::fprintf(outfile,"\n"); 
	fflush(outfile);
	Process::environment.globals["MP2 TOTAL ENERGY"] = Emp2;
	Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = Escsmp2;
	Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
	Process::environment.globals["SCSN-MP2 TOTAL ENERGY"] = Escsnmp2;
	Process::environment.globals["SCS-MP2-VDW TOTAL ENERGY"] = Escsmp2vdw;
	Process::environment.globals["SOS-PI-MP2 TOTAL ENERGY"] = Esospimp2;

	Process::environment.globals["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
	Process::environment.globals["SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
	Process::environment.globals["SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
	Process::environment.globals["SCSN-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;
	Process::environment.globals["SCS-MP2-VDW CORRELATION ENERGY"] = Escsmp2vdw - Escf;
	Process::environment.globals["SOS-PI-MP2 CORRELATION ENERGY"] = Esospimp2 - Escf;

        Process::environment.globals["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
        Process::environment.globals["MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA+Emp2BB;

        timer_on("T2(2)");
	t2_2nd_sc();
        timer_off("T2(2)");
        timer_on("MP3 Energy");
	mp3_energy();
        timer_off("MP3 Energy");
	Emp3L=Emp3;
        EcorrL=Emp3L-Escf;
	Emp3L_old=Emp3;
        if (ip_poles == "TRUE") omp3_ip_poles();
	
	psi::fprintf(outfile,"\n"); 
	psi::fprintf(outfile,"\tComputing MP2.5 energy using SCF MOs (Canonical MP2.5)... \n"); 
	psi::fprintf(outfile,"\t============================================================================== \n");
	psi::fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	psi::fprintf(outfile,"\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	psi::fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	psi::fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp3AA);
	psi::fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp3AB);
	psi::fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp3BB);
	psi::fprintf(outfile,"\t0.5 Energy Correction (a.u.)       : %20.14f\n", Emp3-Emp2);
	psi::fprintf(outfile,"\tMP2.5 Correlation Energy (a.u.)    : %20.14f\n", Ecorr);
	psi::fprintf(outfile,"\tMP2.5 Total Energy (a.u.)          : %20.14f\n", Emp3);
	psi::fprintf(outfile,"\t============================================================================== \n");
	psi::fprintf(outfile,"\n"); 
	fflush(outfile);
	Process::environment.globals["MP2.5 TOTAL ENERGY"] = Emp3;

	omp3_response_pdms();
	gfock();
        //ccl_energy();
	idp();
	mograd();
        occ_iterations();
	
        if (rms_wog <= tol_grad && fabs(DE) >= tol_Eod) {
           orbs_already_opt = 1;
	   if (conver == 1) psi::fprintf(outfile,"\n\tOrbitals are optimized now.\n");
	   else if (conver == 0) { 
                    psi::fprintf(outfile,"\n\tMAX MOGRAD did NOT converged, but RMS MOGRAD converged!!!\n");
	            psi::fprintf(outfile,"\tI will consider the present orbitals as optimized.\n");
           }
	   psi::fprintf(outfile,"\tSwitching to the standard MP2.5 computation after semicanonicalization of the MOs... \n");
	   fflush(outfile);
	   semi_canonic();
	   if (reference_ == "RESTRICTED") trans_ints_rhf();  
	   else if (reference_ == "UNRESTRICTED") trans_ints_uhf();  
	   omp3_t2_1st_sc();
	   t2_2nd_sc();
           conver = 1;
           if (dertype == "FIRST") {
               omp3_response_pdms();
	       gfock();
           }
        }     

  if (conver == 1) {
        ref_energy();
	omp3_mp2_energy();
	mp3_energy();
        if (orbs_already_opt == 1) Emp3L = Emp3;

        if (ip_poles == "TRUE") {
	   if (orbs_already_sc == 0) {
               semi_canonic();
	       if (reference_ == "RESTRICTED") trans_ints_rhf();  
	       else if (reference_ == "UNRESTRICTED") trans_ints_uhf();  
	       omp3_t2_1st_sc();
	       t2_2nd_sc();
           }
           omp3_ip_poles();
        }

        // EKT
        if (ekt_ip_ == "TRUE" || ekt_ea_ == "TRUE") { 
            if (orbs_already_sc == 1) {
	        omp3_response_pdms();
	        gfock();
            }
            gfock_diag();
            if (ekt_ip_ == "TRUE") ekt_ip();
            if (ekt_ea_ == "TRUE") ekt_ea();
        }

        psi::fprintf(outfile,"\n"); 
	psi::fprintf(outfile,"\tComputing MP2 energy using optimized MOs... \n"); 
	psi::fprintf(outfile,"\t============================================================================== \n");
	psi::fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	psi::fprintf(outfile,"\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	psi::fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	psi::fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
	psi::fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
	psi::fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
	psi::fprintf(outfile,"\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA+Escsmp2BB);
	psi::fprintf(outfile,"\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
	psi::fprintf(outfile,"\tSCS-MP2 Total Energy (a.u.)        : %20.14f\n", Escsmp2);
	psi::fprintf(outfile,"\tSOS-MP2 Total Energy (a.u.)        : %20.14f\n", Esosmp2);
	psi::fprintf(outfile,"\tSCSN-MP2 Total Energy (a.u.)       : %20.14f\n", Escsnmp2);
	psi::fprintf(outfile,"\tSCS-MP2-VDW Total Energy (a.u.)    : %20.14f\n", Escsmp2vdw);
	psi::fprintf(outfile,"\tSOS-PI-MP2 Total Energy (a.u.)     : %20.14f\n", Esospimp2);
	psi::fprintf(outfile,"\tMP2 Correlation Energy (a.u.)      : %20.14f\n", Ecorr);
	psi::fprintf(outfile,"\tMP2 Total Energy (a.u.)            : %20.14f\n", Emp2);
	psi::fprintf(outfile,"\t============================================================================== \n");
	psi::fprintf(outfile,"\n"); 
	fflush(outfile);

	psi::fprintf(outfile,"\n"); 
	psi::fprintf(outfile,"\tComputing MP2.5 energy using optimized MOs... \n"); 
	psi::fprintf(outfile,"\t============================================================================== \n");
	psi::fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	psi::fprintf(outfile,"\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	psi::fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	psi::fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp3AA);
	psi::fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp3AB);
	psi::fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp3BB);
	psi::fprintf(outfile,"\t0.5 Energy Correction (a.u.)       : %20.14f\n", Emp3-Emp2);
	psi::fprintf(outfile,"\tMP2.5 Correlation Energy (a.u.)    : %20.14f\n", Ecorr);
	psi::fprintf(outfile,"\tMP2.5 Total Energy (a.u.)          : %20.14f\n", Emp3);
	psi::fprintf(outfile,"\t============================================================================== \n");
	psi::fprintf(outfile,"\n"); 
	fflush(outfile);


	psi::fprintf(outfile,"\n");
	psi::fprintf(outfile,"\t============================================================================== \n");
	psi::fprintf(outfile,"\t================ OMP2.5 FINAL RESULTS ======================================== \n");
	psi::fprintf(outfile,"\t============================================================================== \n");
	psi::fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	psi::fprintf(outfile,"\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	psi::fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	psi::fprintf(outfile,"\tOMP2.5 Correlation Energy (a.u.)   : %20.14f\n", Emp3L-Escf);
	psi::fprintf(outfile,"\tEomp2.5 - Eref (a.u.)              : %20.14f\n", Emp3L-Eref);
	psi::fprintf(outfile,"\tOMP2.5 Total Energy (a.u.)         : %20.14f\n", Emp3L);
	psi::fprintf(outfile,"\t============================================================================== \n");
	psi::fprintf(outfile,"\n");
	fflush(outfile);
	
	// Set the global variables with the energies
	Process::environment.globals["OMP2.5 TOTAL ENERGY"] = Emp3L;
	Process::environment.globals["OMP2.5 CORRELATION ENERGY"] = Emp3L - Escf;
	Process::environment.globals["CURRENT ENERGY"] = Emp3L;
	Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
	Process::environment.globals["CURRENT CORRELATION ENERGY"] = Emp3L-Escf;

	if (natorb == "TRUE") nbo();
	if (occ_orb_energy == "TRUE") semi_canonic(); 

        // Compute Analytic Gradients
        if (dertype == "FIRST") {
            time4grad = 1;
	    psi::fprintf(outfile,"\tAnalytic gradient computation is starting...\n");
	    fflush(outfile);
            coord_grad();
	    psi::fprintf(outfile,"\tNecessary information has been sent to DERIV, which will take care of the rest.\n");
	    fflush(outfile);
        }

  }// end if (conver == 1)
  */
}// end omp2.5_manager 


//======================================================================
//             MP2.5 Manager
//======================================================================             
void DFOCC::mp2_5_manager()
{
        /*
        time4grad = 0;// means i will not compute the gradient
        timer_on("trans_ints");
	if (reference_ == "RESTRICTED") trans_ints_rhf();  
	else if (reference_ == "UNRESTRICTED") trans_ints_uhf();  
        timer_off("trans_ints");
        Eref = Escf;
        timer_on("T2(1)");
	omp3_t2_1st_sc();
        timer_off("T2(1)");
        timer_on("MP2 Energy");
	omp3_mp2_energy();
        timer_off("MP2 Energy");

	psi::fprintf(outfile,"\n"); 
	psi::fprintf(outfile,"\tComputing MP2 energy using SCF MOs (Canonical MP2)... \n"); 
	psi::fprintf(outfile,"\t============================================================================== \n");
	psi::fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	psi::fprintf(outfile,"\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	psi::fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	psi::fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
	psi::fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
	psi::fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
	psi::fprintf(outfile,"\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA+Escsmp2BB);
	psi::fprintf(outfile,"\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
	psi::fprintf(outfile,"\tSCS-MP2 Total Energy (a.u.)        : %20.14f\n", Escsmp2);
	psi::fprintf(outfile,"\tSOS-MP2 Total Energy (a.u.)        : %20.14f\n", Esosmp2);
	psi::fprintf(outfile,"\tSCSN-MP2 Total Energy (a.u.)       : %20.14f\n", Escsnmp2);
	psi::fprintf(outfile,"\tSCS-MP2-VDW Total Energy (a.u.)    : %20.14f\n", Escsmp2vdw);
	psi::fprintf(outfile,"\tSOS-PI-MP2 Total Energy (a.u.)     : %20.14f\n", Esospimp2);
	psi::fprintf(outfile,"\tMP2 Correlation Energy (a.u.)      : %20.14f\n", Ecorr);
	psi::fprintf(outfile,"\tMP2 Total Energy (a.u.)            : %20.14f\n", Emp2);
	psi::fprintf(outfile,"\t============================================================================== \n");
	psi::fprintf(outfile,"\n"); 
	fflush(outfile);
	Process::environment.globals["MP2 TOTAL ENERGY"] = Emp2;
	Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = Escsmp2;
	Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
	Process::environment.globals["SCSN-MP2 TOTAL ENERGY"] = Escsnmp2;
	Process::environment.globals["SCS-MP2-VDW TOTAL ENERGY"] = Escsmp2vdw;
	Process::environment.globals["SOS-PI-MP2 TOTAL ENERGY"] = Esospimp2;

	Process::environment.globals["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
	Process::environment.globals["SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
	Process::environment.globals["SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
	Process::environment.globals["SCSN-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;
	Process::environment.globals["SCS-MP2-VDW CORRELATION ENERGY"] = Escsmp2vdw - Escf;
	Process::environment.globals["SOS-PI-MP2 CORRELATION ENERGY"] = Esospimp2 - Escf;

        Process::environment.globals["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
        Process::environment.globals["MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA+Emp2BB;

        timer_on("T2(2)");
	t2_2nd_sc();
        timer_off("T2(2)");
        timer_on("MP3 Energy");
	mp3_energy();
        timer_off("MP3 Energy");
	Emp3L=Emp3;
        EcorrL=Emp3L-Escf;
	Emp3L_old=Emp3;
        if (ip_poles == "TRUE") omp3_ip_poles();
	
	psi::fprintf(outfile,"\n"); 
	psi::fprintf(outfile,"\tComputing MP2.5 energy using SCF MOs (Canonical MP2.5)... \n"); 
	psi::fprintf(outfile,"\t============================================================================== \n");
	psi::fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	psi::fprintf(outfile,"\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	psi::fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	psi::fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp3AA);
	psi::fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp3AB);
	psi::fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp3BB);
	psi::fprintf(outfile,"\t0.5 Energy Correction (a.u.)       : %20.14f\n", Emp3-Emp2);
	psi::fprintf(outfile,"\tMP2.5 Correlation Energy (a.u.)    : %20.14f\n", Ecorr);
	psi::fprintf(outfile,"\tMP2.5 Total Energy (a.u.)          : %20.14f\n", Emp3);
	psi::fprintf(outfile,"\t============================================================================== \n");
	psi::fprintf(outfile,"\n"); 
	fflush(outfile);
	Process::environment.globals["MP2.5 TOTAL ENERGY"] = Emp3;
	Process::environment.globals["MP2.5 CORRELATION ENERGY"] = Emp3 - Escf;
	Process::environment.globals["CURRENT ENERGY"] = Emp3L;
	Process::environment.globals["CURRENT REFERENCE ENERGY"] = Eref;
	Process::environment.globals["CURRENT CORRELATION ENERGY"] = Emp3L-Escf;

        // Compute Analytic Gradients
        if (dertype == "FIRST" || ekt_ip_ == "TRUE" || ekt_ea_ == "TRUE") {
            time4grad = 1;
	    psi::fprintf(outfile,"\tAnalytic gradient computation is starting...\n");
            psi::fprintf(outfile,"\tComputing response density matrices...\n");
            fflush(outfile);
	    omp3_response_pdms();
            psi::fprintf(outfile,"\tComputing off-diagonal blocks of GFM...\n");
            fflush(outfile);
	    gfock();
            psi::fprintf(outfile,"\tForming independent-pairs...\n");
            fflush(outfile);
	    idp2();
            psi::fprintf(outfile,"\tComputing orbital gradient...\n");
            fflush(outfile);
	    mograd();
            coord_grad();

            if (ekt_ip_ == "TRUE" && ekt_ea_ == "TRUE") {
                ekt_ip();
                ekt_ea();
            }

            else if (ekt_ip_ == "TRUE" && ekt_ea_ == "FALSE") {
                ekt_ip();
            }

            else if (ekt_ip_ == "FALSE" && ekt_ea_ == "TRUE") {
                ekt_ea();
            }

            else if (ekt_ip_ == "FALSE" && ekt_ea_ == "FALSE") {
	        psi::fprintf(outfile,"\tNecessary information has been sent to DERIV, which will take care of the rest.\n");
	        fflush(outfile);
            }

        }
        */

}// end omp2.5_manager 


}} // End Namespaces



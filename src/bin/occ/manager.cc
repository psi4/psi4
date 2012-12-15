/** Standard library includes */
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <sstream>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>

/** Required PSI3 includes */ 
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <libiwl/iwl.h>
#include <libqt/qt.h>
#include <libtrans/mospace.h>
#include <libtrans/integraltransform.h>

/** Required libmints includes */
#include <libmints/mints.h>
#include <libmints/factory.h>
#include <libmints/wavefunction.h>

#include "defines.h"
#include "occwave.h"

using namespace boost;
using namespace psi;
using namespace std;


namespace psi{ namespace occwave{
  
void OCCWave::omp2_manager()
{
	mo_optimized = 0;
	orbs_already_opt = 0;
        timer_on("trans_ints");
	if (reference_ == "RESTRICTED") trans_ints_rhf();  
	else if (reference_ == "UNRESTRICTED") trans_ints_uhf();  
        timer_off("trans_ints");
        timer_on("T2(1)");
	omp2_t2_1st_sc();
        timer_off("T2(1)");
        timer_on("REF Energy");
	ref_energy();
        timer_off("REF Energy");
        timer_on("MP2 Energy");
	omp2_mp2_energy();
        timer_off("MP2 Energy");
	Emp2L=Emp2;
        EcorrL=Emp2L-Escf;
	Emp2L_old=Emp2;
	
	fprintf(outfile,"\n"); 
	fprintf(outfile,"\tComputing MP2 energy using SCF MOs (Canonical MP2)... \n"); 
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %12.14f\n", Enuc);
	fprintf(outfile,"\tSCF Energy (a.u.)                  : %12.14f\n", Escf);
	fprintf(outfile,"\tREF Energy (a.u.)                  : %12.14f\n", Eref);
	fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %12.14f\n", Emp2AA);
	fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %12.14f\n", Emp2AB);
	fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %12.14f\n", Emp2BB);
	fprintf(outfile,"\tScaled_SS Correlation Energy (a.u.): %12.14f\n", Escsmp2AA+Escsmp2BB);
	fprintf(outfile,"\tScaled_OS Correlation Energy (a.u.): %12.14f\n", Escsmp2AB);
	fprintf(outfile,"\tSCS-MP2 Total Energy (a.u.)        : %12.14f\n", Escsmp2);
	fprintf(outfile,"\tSOS-MP2 Total Energy (a.u.)        : %12.14f\n", Esosmp2);
	fprintf(outfile,"\tSCSN-MP2 Total Energy (a.u.)       : %12.14f\n", Escsnmp2);
	fprintf(outfile,"\tSCS-MI-MP2 Total Energy (a.u.)     : %12.14f\n", Escsmimp2);
	fprintf(outfile,"\tSCS-MP2-VDW Total Energy (a.u.)    : %12.14f\n", Escsmp2vdw);
	fprintf(outfile,"\tSOS-PI-MP2 Total Energy (a.u.)     : %12.14f\n", Esospimp2);
	fprintf(outfile,"\tMP2 Correlation Energy (a.u.)      : %12.14f\n", Ecorr);
	fprintf(outfile,"\tMP2 Total Energy (a.u.)            : %12.14f\n", Emp2);
	fprintf(outfile,"\t============================================================================== \n");
	fflush(outfile);
	Process::environment.globals["MP2 TOTAL ENERGY"] = Emp2;
	Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = Escsmp2;
	Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
	Process::environment.globals["SCSN-MP2 TOTAL ENERGY"] = Escsnmp2;
	Process::environment.globals["SCS-MI-MP2 TOTAL ENERGY"] = Escsmimp2;
	Process::environment.globals["SCS-MP2-VDW TOTAL ENERGY"] = Escsmp2vdw;
	Process::environment.globals["SOS-PI-MP2 TOTAL ENERGY"] = Esospimp2;

	omp2_response_pdms();
	gfock();
	idp();
	mograd();
        occ_iterations();
	
        if (rms_wog <= tol_grad && fabs(DE) >= tol_Eod) {
           orbs_already_opt = 1;
	   if (conver == 1) fprintf(outfile,"\n\tOrbitals are optimized now.\n");
	   else if (conver == 0) { 
                    fprintf(outfile,"\n\tMAX MOGRAD did NOT converged, but RMS MOGRAD converged!!!\n");
	            fprintf(outfile,"\tI will consider the present orbitals as optimized.\n");
           }
	   fprintf(outfile,"\tSwitching to the standard MP2 computation after semicanonicalization of the MOs... \n");
	   fflush(outfile);
	   semi_canonic();
	   if (reference_ == "RESTRICTED") trans_ints_rhf();  
	   else if (reference_ == "UNRESTRICTED") trans_ints_uhf();  
	   omp2_t2_1st_sc();
           conver = 1;
           if (dertype == "FIRST") {
               omp2_response_pdms();
	       gfock();
           }
        }     

  if (conver == 1) {
        ref_energy();
	omp2_mp2_energy();
        if (orbs_already_opt == 1) Emp2L = Emp2;
	
	fprintf(outfile,"\n"); 
	fprintf(outfile,"\tComputing MP2 energy using optimized MOs... \n"); 
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %12.14f\n", Enuc);
	fprintf(outfile,"\tSCF Energy (a.u.)                  : %12.14f\n", Escf);
	fprintf(outfile,"\tREF Energy (a.u.)                  : %12.14f\n", Eref);
	fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %12.14f\n", Emp2AA);
	fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %12.14f\n", Emp2AB);
	fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %12.14f\n", Emp2BB);
	fprintf(outfile,"\tScaled_SS Correlation Energy (a.u.): %12.14f\n", Escsmp2AA+Escsmp2BB);
	fprintf(outfile,"\tScaled_OS Correlation Energy (a.u.): %12.14f\n", Escsmp2AB);
	fprintf(outfile,"\tSCS-MP2 Total Energy (a.u.)        : %12.14f\n", Escsmp2);
	fprintf(outfile,"\tSOS-MP2 Total Energy (a.u.)        : %12.14f\n", Esosmp2);
	fprintf(outfile,"\tSCSN-MP2 Total Energy (a.u.)       : %12.14f\n", Escsnmp2);
	fprintf(outfile,"\tSCS-MI-MP2 Total Energy (a.u.)     : %12.14f\n", Escsmimp2);
	fprintf(outfile,"\tSCS-MP2-VDW Total Energy (a.u.)    : %12.14f\n", Escsmp2vdw);
	fprintf(outfile,"\tSOS-PI-MP2 Total Energy (a.u.)     : %12.14f\n", Esospimp2);
	fprintf(outfile,"\tMP2 Correlation Energy (a.u.)      : %12.14f\n", Ecorr);
	fprintf(outfile,"\tMP2 Total Energy (a.u.)            : %12.14f\n", Emp2);
	fprintf(outfile,"\t============================================================================== \n");
	fflush(outfile);


	fprintf(outfile,"\n");
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\t================ OMP2 FINAL RESULTS ========================================== \n");
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %12.14f\n", Enuc);
	fprintf(outfile,"\tSCF Energy (a.u.)                  : %12.14f\n", Escf);
	fprintf(outfile,"\tREF Energy (a.u.)                  : %12.14f\n", Eref);
	fprintf(outfile,"\tSCS-OMP2 Total Energy (a.u.)       : %12.14f\n", Escsmp2);
	fprintf(outfile,"\tSOS-OMP2 Total Energy (a.u.)       : %12.14f\n", Esosmp2);
	fprintf(outfile,"\tSCSN-OMP2 Total Energy (a.u.)      : %12.14f\n", Escsnmp2);
	fprintf(outfile,"\tSCS-MI-OMP2 Total Energy (a.u.)    : %12.14f\n", Escsmimp2);
	fprintf(outfile,"\tSCS-OMP2-VDW Total Energy (a.u.)   : %12.14f\n", Escsmp2vdw);
	fprintf(outfile,"\tSOS-PI-OMP2 Total Energy (a.u.)    : %12.14f\n", Esospimp2);
	fprintf(outfile,"\tOMP2 Correlation Energy (a.u.)     : %12.14f\n", Emp2L-Escf);
	fprintf(outfile,"\tEomp2 - Eref (a.u.)                : %12.14f\n", Emp2L-Eref);
	fprintf(outfile,"\tOMP2 Total Energy (a.u.)           : %12.14f\n", Emp2L);
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\n");
	fflush(outfile);

	// Set the global variables with the energies
	Process::environment.globals["OMP2 TOTAL ENERGY"] = Emp2L;
	Process::environment.globals["SCS-OMP2 TOTAL ENERGY"] =  Escsmp2;
	Process::environment.globals["SOS-OMP2 TOTAL ENERGY"] =  Esosmp2;
	Process::environment.globals["SCSN-OMP2 TOTAL ENERGY"] = Escsnmp2;
	Process::environment.globals["SCS-MI-OMP2 TOTAL ENERGY"] = Escsmimp2;
	Process::environment.globals["SCS-OMP2-VDW TOTAL ENERGY"] = Escsmp2vdw;
	Process::environment.globals["SOS-PI-OMP2 TOTAL ENERGY"] = Esospimp2;
	Process::environment.globals["CURRENT ENERGY"] = Emp2L;
	Process::environment.globals["CURRENT REFERENCE ENERGY"] = Eref;
	Process::environment.globals["CURRENT CORRELATION ENERGY"] = Emp2L-Escf;

        // if scs on	
	if (do_scs == "TRUE") {
	    if (scs_type_ == "SCS") {
	       Process::environment.globals["CURRENT ENERGY"] = Escsmp2;
            }

	    else if (scs_type_ == "SCSN") {
	       Process::environment.globals["CURRENT ENERGY"] = Escsnmp2;
            }

	    else if (scs_type_ == "SCSMI") {
	       Process::environment.globals["CURRENT ENERGY"] = Escsmimp2;
            }

	    else if (scs_type_ == "SCSVDW") {
	       Process::environment.globals["CURRENT ENERGY"] = Escsmp2vdw;
            }
	}
    
        // else if sos on	
	else if (do_sos == "TRUE") {
	     if (sos_type_ == "SOS") {
	         Process::environment.globals["CURRENT ENERGY"] = Esosmp2;
             }

	     else if (sos_type_ == "SOSPI") {
	             Process::environment.globals["CURRENT ENERGY"] = Esospimp2;
             }
	}
 
	if (natorb == "TRUE") nbo();
	if (occ_orb_energy == "TRUE") semi_canonic(); 

        // Compute Analytic Gradients
        if (dertype == "FIRST") coord_grad();

  }// end if (conver == 1)
}// end omp2_manager 


void OCCWave::omp3_manager()
{
	mo_optimized = 0;
	orbs_already_opt = 0;
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

	fprintf(outfile,"\n"); 
	fprintf(outfile,"\tComputing MP2 energy using SCF MOs (Canonical MP2)... \n"); 
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %12.14f\n", Enuc);
	fprintf(outfile,"\tSCF Energy (a.u.)                  : %12.14f\n", Escf);
	fprintf(outfile,"\tREF Energy (a.u.)                  : %12.14f\n", Eref);
	fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %12.14f\n", Emp2AA);
	fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %12.14f\n", Emp2AB);
	fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %12.14f\n", Emp2BB);
	fprintf(outfile,"\tScaled_SS Correlation Energy (a.u.): %12.14f\n", Escsmp2AA+Escsmp2BB);
	fprintf(outfile,"\tScaled_OS Correlation Energy (a.u.): %12.14f\n", Escsmp2AB);
	fprintf(outfile,"\tSCS-MP2 Total Energy (a.u.)        : %12.14f\n", Escsmp2);
	fprintf(outfile,"\tSOS-MP2 Total Energy (a.u.)        : %12.14f\n", Esosmp2);
	fprintf(outfile,"\tSCSN-MP2 Total Energy (a.u.)       : %12.14f\n", Escsnmp2);
	fprintf(outfile,"\tSCS-MI-MP2 Total Energy (a.u.)     : %12.14f\n", Escsmimp2);
	fprintf(outfile,"\tSCS-MP2-VDW Total Energy (a.u.)    : %12.14f\n", Escsmp2vdw);
	fprintf(outfile,"\tSOS-PI-MP2 Total Energy (a.u.)     : %12.14f\n", Esospimp2);
	fprintf(outfile,"\tMP2 Correlation Energy (a.u.)      : %12.14f\n", Ecorr);
	fprintf(outfile,"\tMP2 Total Energy (a.u.)            : %12.14f\n", Emp2);
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\n"); 
	fflush(outfile);
	Process::environment.globals["MP2 TOTAL ENERGY"] = Emp2;
	Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = Escsmp2;
	Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
	Process::environment.globals["SCSN-MP2 TOTAL ENERGY"] = Escsnmp2;
	Process::environment.globals["SCS-MI-MP2 TOTAL ENERGY"] = Escsmimp2;
	Process::environment.globals["SCS-MP2-VDW TOTAL ENERGY"] = Escsmp2vdw;
	Process::environment.globals["SOS-PI-MP2 TOTAL ENERGY"] = Esospimp2;

        timer_on("T2(2)");
	t2_2nd_sc();
        timer_off("T2(2)");
        timer_on("MP3 Energy");
	mp3_energy();
        timer_off("MP3 Energy");
	Emp3L=Emp3;
        EcorrL=Emp3L-Escf;
	Emp3L_old=Emp3;
	
	fprintf(outfile,"\n"); 
	fprintf(outfile,"\tComputing MP3 energy using SCF MOs (Canonical MP3)... \n"); 
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %12.14f\n", Enuc);
	fprintf(outfile,"\tSCF Energy (a.u.)                  : %12.14f\n", Escf);
	fprintf(outfile,"\tREF Energy (a.u.)                  : %12.14f\n", Eref);
	fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %12.14f\n", Emp3AA);
	fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %12.14f\n", Emp3AB);
	fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %12.14f\n", Emp3BB);
	fprintf(outfile,"\tSCS-MP3 Total Energy (a.u.)        : %12.14f\n", Escsmp3);
	fprintf(outfile,"\tSOS-MP3 Total Energy (a.u.)        : %12.14f\n", Esosmp3);
	fprintf(outfile,"\tSCSN-MP3 Total Energy (a.u.)       : %12.14f\n", Escsnmp3);
	fprintf(outfile,"\tSCS-MI-MP3 Total Energy (a.u.)     : %12.14f\n", Escsmimp3);
	fprintf(outfile,"\tSCS-MP3-VDW Total Energy (a.u.)    : %12.14f\n", Escsmp3vdw);
	fprintf(outfile,"\tSOS-PI-MP3 Total Energy (a.u.)     : %12.14f\n", Esospimp3);
	fprintf(outfile,"\t3rd Order Energy (a.u.)            : %12.14f\n", Emp3-Emp2);
	fprintf(outfile,"\tMP3 Correlation Energy (a.u.)      : %12.14f\n", Ecorr);
	fprintf(outfile,"\tMP3 Total Energy (a.u.)            : %12.14f\n", Emp3);
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\n"); 
	fflush(outfile);
	Process::environment.globals["MP3 TOTAL ENERGY"] = Emp3;
	Process::environment.globals["SCS-MP3 TOTAL ENERGY"] = Escsmp3;
	Process::environment.globals["SOS-MP3 TOTAL ENERGY"] = Esosmp3;
	Process::environment.globals["SCSN-MP3 TOTAL ENERGY"] = Escsnmp3;
	Process::environment.globals["SCS-MI-MP3 TOTAL ENERGY"] = Escsmimp3;
	Process::environment.globals["SCS-MP3-VDW TOTAL ENERGY"] = Escsmp3vdw;
	Process::environment.globals["SOS-PI-MP3 TOTAL ENERGY"] = Esospimp3;

	omp3_response_pdms();
	gfock();
	
	idp();
	mograd();
        occ_iterations();
	
        if (rms_wog <= tol_grad && fabs(DE) >= tol_Eod) {
           orbs_already_opt = 1;
	   if (conver == 1) fprintf(outfile,"\n\tOrbitals are optimized now.\n");
	   else if (conver == 0) { 
                    fprintf(outfile,"\n\tMAX MOGRAD did NOT converged, but RMS MOGRAD converged!!!\n");
	            fprintf(outfile,"\tI will consider the present orbitals as optimized.\n");
           }
	   fprintf(outfile,"\tSwitching to the standard MP3 computation after semicanonicalization of the MOs... \n");
	   fflush(outfile);
	   semi_canonic();
	   if (reference_ == "RESTRICTED") trans_ints_rhf();  
	   else if (reference_ == "UNRESTRICTED") trans_ints_uhf();  
	   omp3_t2_1st_sc();
	   t2_2nd_sc();
           conver = 1;
        }     

  if (conver == 1) {
        ref_energy();
	omp3_mp2_energy();
	mp3_energy();
        if (orbs_already_opt == 1) Emp3L = Emp3;

        fprintf(outfile,"\n"); 
	fprintf(outfile,"\tComputing MP2 energy using optimized MOs... \n"); 
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %12.14f\n", Enuc);
	fprintf(outfile,"\tSCF Energy (a.u.)                  : %12.14f\n", Escf);
	fprintf(outfile,"\tREF Energy (a.u.)                  : %12.14f\n", Eref);
	fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %12.14f\n", Emp2AA);
	fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %12.14f\n", Emp2AB);
	fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %12.14f\n", Emp2BB);
	fprintf(outfile,"\tScaled_SS Correlation Energy (a.u.): %12.14f\n", Escsmp2AA+Escsmp2BB);
	fprintf(outfile,"\tScaled_OS Correlation Energy (a.u.): %12.14f\n", Escsmp2AB);
	fprintf(outfile,"\tSCS-MP2 Total Energy (a.u.)        : %12.14f\n", Escsmp2);
	fprintf(outfile,"\tSOS-MP2 Total Energy (a.u.)        : %12.14f\n", Esosmp2);
	fprintf(outfile,"\tSCSN-MP2 Total Energy (a.u.)       : %12.14f\n", Escsnmp2);
	fprintf(outfile,"\tSCS-MI-MP2 Total Energy (a.u.)     : %12.14f\n", Escsmimp2);
	fprintf(outfile,"\tSCS-MP2-VDW Total Energy (a.u.)    : %12.14f\n", Escsmp2vdw);
	fprintf(outfile,"\tSOS-PI-MP2 Total Energy (a.u.)     : %12.14f\n", Esospimp2);
	fprintf(outfile,"\tMP2 Correlation Energy (a.u.)      : %12.14f\n", Ecorr);
	fprintf(outfile,"\tMP2 Total Energy (a.u.)            : %12.14f\n", Emp2);
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\n"); 
	fflush(outfile);

	fprintf(outfile,"\n"); 
	fprintf(outfile,"\tComputing MP3 energy using optimized MOs... \n"); 
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %12.14f\n", Enuc);
	fprintf(outfile,"\tSCF Energy (a.u.)                  : %12.14f\n", Escf);
	fprintf(outfile,"\tREF Energy (a.u.)                  : %12.14f\n", Eref);
	fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %12.14f\n", Emp3AA);
	fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %12.14f\n", Emp3AB);
	fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %12.14f\n", Emp3BB);
	fprintf(outfile,"\tSCS-MP3 Total Energy (a.u.)        : %12.14f\n", Escsmp3);
	fprintf(outfile,"\tSOS-MP3 Total Energy (a.u.)        : %12.14f\n", Esosmp3);
	fprintf(outfile,"\tSCSN-MP3 Total Energy (a.u.)       : %12.14f\n", Escsnmp3);
	fprintf(outfile,"\tSCS-MI-MP3 Total Energy (a.u.)     : %12.14f\n", Escsmimp3);
	fprintf(outfile,"\tSCS-MP3-VDW Total Energy (a.u.)    : %12.14f\n", Escsmp3vdw);
	fprintf(outfile,"\tSOS-PI-MP3 Total Energy (a.u.)     : %12.14f\n", Esospimp3);
	fprintf(outfile,"\t3rd Order Energy (a.u.)            : %12.14f\n", Emp3-Emp2);
	fprintf(outfile,"\tMP3 Correlation Energy (a.u.)      : %12.14f\n", Ecorr);
	fprintf(outfile,"\tMP3 Total Energy (a.u.)            : %12.14f\n", Emp3);
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\n"); 
	fflush(outfile);


	fprintf(outfile,"\n");
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\t================ OMP3 FINAL RESULTS ========================================== \n");
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %12.14f\n", Enuc);
	fprintf(outfile,"\tSCF Energy (a.u.)                  : %12.14f\n", Escf);
	fprintf(outfile,"\tREF Energy (a.u.)                  : %12.14f\n", Eref);
	fprintf(outfile,"\tSCS-OMP3 Total Energy (a.u.)       : %12.14f\n", Escsmp3);
	fprintf(outfile,"\tSOS-OMP3 Total Energy (a.u.)       : %12.14f\n", Esosmp3);
	fprintf(outfile,"\tSCSN-OMP3 Total Energy (a.u.)      : %12.14f\n", Escsnmp3);
	fprintf(outfile,"\tSCS-MI-OMP3 Total Energy (a.u.)    : %12.14f\n", Escsmimp3);
	fprintf(outfile,"\tSCS-OMP3-VDW Total Energy (a.u.    : %12.14f\n", Escsmp3vdw);
	fprintf(outfile,"\tSOS-PI-OMP3 Total Energy (a.u.)    : %12.14f\n", Esospimp3);
	fprintf(outfile,"\tOMP3 Correlation Energy (a.u.)     : %12.14f\n", Emp3L-Escf);
	fprintf(outfile,"\tEomp3 - Eref (a.u.)                : %12.14f\n", Emp3L-Eref);
	fprintf(outfile,"\tOMP3 Total Energy (a.u.)           : %12.14f\n", Emp3L);
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\n");
	fflush(outfile);
	
	// Set the global variables with the energies
	Process::environment.globals["OMP3 TOTAL ENERGY"] = Emp3L;
	Process::environment.globals["SCS-OMP3 TOTAL ENERGY"] =  Escsmp3;
	Process::environment.globals["SOS-OMP3 TOTAL ENERGY"] =  Esosmp3;
	Process::environment.globals["SCSN-OMP3 TOTAL ENERGY"] = Escsnmp3;
	Process::environment.globals["SCS-MI-OMP3 TOTAL ENERGY"] = Escsmimp3;
	Process::environment.globals["SCS-OMP3-VDW TOTAL ENERGY"] = Escsmp3vdw;
	Process::environment.globals["SOS-PI-OMP3 TOTAL ENERGY"] = Esospimp3;
	Process::environment.globals["CURRENT ENERGY"] = Emp3L;
	Process::environment.globals["CURRENT REFERENCE ENERGY"] = Eref;
	Process::environment.globals["CURRENT CORRELATION ENERGY"] = Emp3L-Escf;

        // if scs on	
	if (do_scs == "TRUE") {
	    if (scs_type_ == "SCS") {
	       Process::environment.globals["CURRENT ENERGY"] = Escsmp3;
            }

	    else if (scs_type_ == "SCSN") {
	       Process::environment.globals["CURRENT ENERGY"] = Escsnmp3;
            }

	    else if (scs_type_ == "SCSMI") {
	       Process::environment.globals["CURRENT ENERGY"] = Escsmimp3;
            }

	    else if (scs_type_ == "SCSVDW") {
	       Process::environment.globals["CURRENT ENERGY"] = Escsmp3vdw;
            }
	}
    
        // else if sos on	
	else if (do_sos == "TRUE") {
	     if (sos_type_ == "SOS") {
	         Process::environment.globals["CURRENT ENERGY"] = Esosmp3;
             }

	     else if (sos_type_ == "SOSPI") {
	             Process::environment.globals["CURRENT ENERGY"] = Esospimp3;
             }
	}

	if (natorb == "TRUE") nbo();
	if (occ_orb_energy == "TRUE") semi_canonic(); 
  }// end if (conver == 1)
}// end omp3_manager 


void OCCWave::ocepa_manager()
{
	mo_optimized = 0;// means MOs are not optimized yet.
	orbs_already_opt = 0;
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

	fprintf(outfile,"\n"); 
	fprintf(outfile,"\tComputing MP2 energy using SCF MOs (Canonical MP2)... \n"); 
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %12.14f\n", Enuc);
	fprintf(outfile,"\tSCF Energy (a.u.)                  : %12.14f\n", Escf);
	fprintf(outfile,"\tREF Energy (a.u.)                  : %12.14f\n", Eref);
	fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %12.14f\n", Emp2AA);
	fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %12.14f\n", Emp2AB);
	fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %12.14f\n", Emp2BB);
	fprintf(outfile,"\tScaled_SS Correlation Energy (a.u.): %12.14f\n", Escsmp2AA+Escsmp2BB);
	fprintf(outfile,"\tScaled_OS Correlation Energy (a.u.): %12.14f\n", Escsmp2AB);
	fprintf(outfile,"\tSCS-MP2 Total Energy (a.u.)        : %12.14f\n", Escsmp2);
	fprintf(outfile,"\tSOS-MP2 Total Energy (a.u.)        : %12.14f\n", Esosmp2);
	fprintf(outfile,"\tSCSN-MP2 Total Energy (a.u.)       : %12.14f\n", Escsnmp2);
	fprintf(outfile,"\tSCS-MI-MP2 Total Energy (a.u.)     : %12.14f\n", Escsmimp2);
	fprintf(outfile,"\tSCS-MP2-VDW Total Energy (a.u.)    : %12.14f\n", Escsmp2vdw);
	fprintf(outfile,"\tSOS-PI-MP2 Total Energy (a.u.)     : %12.14f\n", Esospimp2);
	fprintf(outfile,"\tMP2 Correlation Energy (a.u.)      : %12.14f\n", Ecorr);
	fprintf(outfile,"\tMP2 Total Energy (a.u.)            : %12.14f\n", Emp2);
	fprintf(outfile,"\t============================================================================== \n");
	fflush(outfile);
	Process::environment.globals["MP2 TOTAL ENERGY"] = Emp2;
	Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = Escsmp2;
	Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
	Process::environment.globals["SCSN-MP2 TOTAL ENERGY"] = Escsnmp2;
	Process::environment.globals["SCS-MI-MP2 TOTAL ENERGY"] = Escsmimp2;
	Process::environment.globals["SCS-MP2-VDW TOTAL ENERGY"] = Escsmp2vdw;
	Process::environment.globals["SOS-PI-MP2 TOTAL ENERGY"] = Esospimp2;

	ocepa_response_pdms();
	gfock();
	idp();
	mograd();
        if (rms_wog > tol_grad) occ_iterations();
        else {
	   fprintf(outfile,"\n\tOrbitals are already optimized, switching to the canonical CEPA computation... \n");
	   fflush(outfile);
           cepa_iterations();
        }
	
        if (rms_wog <= tol_grad && fabs(DE) >= tol_Eod) {
           orbs_already_opt = 1;
	   if (conver == 1) fprintf(outfile,"\n\tOrbitals are optimized now.\n");
	   else if (conver == 0) { 
                    fprintf(outfile,"\n\tMAX MOGRAD did NOT converged, but RMS MOGRAD converged!!!\n");
	            fprintf(outfile,"\tI will consider the present orbitals as optimized.\n");
           }
	   fprintf(outfile,"\tSwitching to the standard CEPA computation... \n");
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

	fprintf(outfile,"\n");
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\t================ OCEPA FINAL RESULTS ========================================= \n");
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %12.14f\n", Enuc);
	fprintf(outfile,"\tSCF Energy (a.u.)                  : %12.14f\n", Escf);
	fprintf(outfile,"\tREF Energy (a.u.)                  : %12.14f\n", Eref);
	fprintf(outfile,"\tSCS-OCEPA(0) Total Energy (a.u.)   : %12.14f\n", Escscepa);
	fprintf(outfile,"\tSOS-OCEPA(0) Total Energy (a.u.)   : %12.14f\n", Esoscepa);
	fprintf(outfile,"\tOCEPA(0) Correlation Energy (a.u.) : %12.14f\n", EcepaL-Escf);
	fprintf(outfile,"\tEocepa - Eref (a.u.)               : %12.14f\n", EcepaL-Eref);
	fprintf(outfile,"\tOCEPA(0) Total Energy (a.u.)       : %12.14f\n", EcepaL);
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\n");
	fflush(outfile);
	
	// Set the global variables with the energies
	Process::environment.globals["OCEPA(0) TOTAL ENERGY"] = EcepaL;
	Process::environment.globals["SCS-OCEPA(0) TOTAL ENERGY"] =  Escscepa;
	Process::environment.globals["SOS-OCEPA(0) TOTAL ENERGY"] =  Esoscepa;
	Process::environment.globals["CURRENT ENERGY"] = EcepaL;
	Process::environment.globals["CURRENT REFERENCE ENERGY"] = Eref;
	Process::environment.globals["CURRENT CORRELATION ENERGY"] = EcepaL-Escf;

        // if scs on	
	if (do_scs == "TRUE") {
	    Process::environment.globals["CURRENT ENERGY"] = Escscepa;

	}
    
        // else if sos on	
	else if (do_sos == "TRUE") {
	         Process::environment.globals["CURRENT ENERGY"] = Esoscepa;
	}

	if (natorb == "TRUE") nbo();
	if (occ_orb_energy == "TRUE") semi_canonic(); 
	
        // Compute Analytic Gradients
        if (dertype == "FIRST") {
            time4grad = 1;
            coord_grad();
        }

  }// end if (conver == 1)
}// end ocepa_manager 


void OCCWave::cepa_manager()
{
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
	Ecepa_old = Emp2;

	fprintf(outfile,"\n"); 
	fprintf(outfile,"\tComputing MP2 energy using SCF MOs (Canonical MP2)... \n"); 
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %12.14f\n", Enuc);
	fprintf(outfile,"\tSCF Energy (a.u.)                  : %12.14f\n", Escf);
	fprintf(outfile,"\tREF Energy (a.u.)                  : %12.14f\n", Eref);
	fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %12.14f\n", Emp2AA);
	fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %12.14f\n", Emp2AB);
	fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %12.14f\n", Emp2BB);
	fprintf(outfile,"\tScaled_SS Correlation Energy (a.u.): %12.14f\n", Escsmp2AA+Escsmp2BB);
	fprintf(outfile,"\tScaled_OS Correlation Energy (a.u.): %12.14f\n", Escsmp2AB);
	fprintf(outfile,"\tSCS-MP2 Total Energy (a.u.)        : %12.14f\n", Escsmp2);
	fprintf(outfile,"\tSOS-MP2 Total Energy (a.u.)        : %12.14f\n", Esosmp2);
	fprintf(outfile,"\tSCSN-MP2 Total Energy (a.u.)       : %12.14f\n", Escsnmp2);
	fprintf(outfile,"\tSCS-MI-MP2 Total Energy (a.u.)     : %12.14f\n", Escsmimp2);
	fprintf(outfile,"\tSCS-MP2-VDW Total Energy (a.u.)    : %12.14f\n", Escsmp2vdw);
	fprintf(outfile,"\tSOS-PI-MP2 Total Energy (a.u.)     : %12.14f\n", Esospimp2);
	fprintf(outfile,"\tMP2 Correlation Energy (a.u.)      : %12.14f\n", Ecorr);
	fprintf(outfile,"\tMP2 Total Energy (a.u.)            : %12.14f\n", Emp2);
	fprintf(outfile,"\t============================================================================== \n");
	fflush(outfile);
	Process::environment.globals["MP2 TOTAL ENERGY"] = Emp2;
	Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = Escsmp2;
	Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
	Process::environment.globals["SCSN-MP2 TOTAL ENERGY"] = Escsnmp2;
	Process::environment.globals["SCS-MI-MP2 TOTAL ENERGY"] = Escsmimp2;
	Process::environment.globals["SCS-MP2-VDW TOTAL ENERGY"] = Escsmp2vdw;
	Process::environment.globals["SOS-PI-MP2 TOTAL ENERGY"] = Esospimp2;

        // Perform CEPA iterations
        cepa_iterations();

	fprintf(outfile,"\n");
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\t================ CEPA FINAL RESULTS ========================================== \n");
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %12.14f\n", Enuc);
	fprintf(outfile,"\tSCF Energy (a.u.)                  : %12.14f\n", Escf);
	fprintf(outfile,"\tREF Energy (a.u.)                  : %12.14f\n", Eref);
	fprintf(outfile,"\tCEPA(0) Correlation Energy (a.u.)  : %12.14f\n", Ecorr);
	fprintf(outfile,"\tCEPA(0) Total Energy (a.u.)        : %12.14f\n", Ecepa);
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\n");
	fflush(outfile);
	
	// Set the global variables with the energies
	Process::environment.globals["CEPA(0) TOTAL ENERGY"] = Ecepa;
	Process::environment.globals["CURRENT ENERGY"] = Ecepa;
	Process::environment.globals["CURRENT REFERENCE ENERGY"] = Eref;
	Process::environment.globals["CURRENT CORRELATION ENERGY"] = Ecorr;
        //EcepaL = Ecepa;
}// end cepa_manager 

}} // End Namespaces



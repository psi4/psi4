#include <libqt/qt.h>
#include <libtrans/integraltransform.h>

#include "occwave.h"
#include "defines.h"

using namespace boost;
using namespace psi;
using namespace std;

namespace psi{ namespace plugin_occ{

void OCCWave::trans_ints_rhf()
{    
    //fprintf(outfile,"\n trans_ints is starting... \n"); fflush(outfile);
/********************************************************************************************/
/************************** Transform 2-electron int. to MO space ***************************/
/********************************************************************************************/  
    ints->update_orbitals();  
    ints->set_print(print_ - 2 >= 0 ? print_ - 2 : 0);
    ints->set_keep_dpd_so_ints(1);
    
    // Trans (OO|OO)
    timer_on("Trans (OO|OO)");
    ints->transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::occ, MOSpace::occ, IntegralTransform::MakeAndKeep);
    timer_off("Trans (OO|OO)");
    
    // Trans (OO|OV)
    timer_on("Trans (OO|OV)");
    ints->transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::occ, MOSpace::vir, IntegralTransform::ReadAndKeep);
    timer_off("Trans (OO|OV)");
    
    // Trans (OO|VV)
    timer_on("Trans (OO|VV)");
    ints->transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::vir, MOSpace::vir, IntegralTransform::ReadAndNuke);
    timer_off("Trans (OO|VV)");

    // Trans (OV|OV)
    timer_on("Trans (OV|OV)");
    ints->transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::occ, MOSpace::vir, IntegralTransform::MakeAndKeep);
    timer_off("Trans (OV|OV)");
    
    // Trans (OV|VV)
    timer_on("Trans (OV|VV)");
    ints->transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::vir, MOSpace::vir, IntegralTransform::ReadAndNuke);
    timer_off("Trans (OV|VV)");
    
//if (wfn_type_ == "OMP3" || wfn_type_ == "OCEPA" || wfn_type_ == "CEPA") { 
if (wfn_type_ != "OMP2") { 
    // Trans (VV|VV)
    timer_on("Trans (VV|VV)");
    ints->transform_tei(MOSpace::vir, MOSpace::vir, MOSpace::vir, MOSpace::vir);
    timer_off("Trans (VV|VV)");
} 
    
/********************************************************************************************/
/************************** sort chem -> phys ***********************************************/
/********************************************************************************************/  
     dpdbuf4 K, G;
     
     psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
     
     // Build MO ints    
     timer_on("Sort chem -> phys");
     timer_on("Sort (OO|OO) -> <OO|OO>");
     // (OO|OO) -> <OO|OO>
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O>=O]+"), ID("[O>=O]+"), 0, "MO Ints (OO|OO)");
     dpd_buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[O,O]"), ID("[O,O]"), "MO Ints <OO|OO>");
     dpd_buf4_close(&K);
     timer_off("Sort (OO|OO) -> <OO|OO>");
 
  
     timer_on("Sort (OO|OV) -> <OO|OV>");
     // (OO|OV) -> <OO|OV>
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,V]"),
                  ID("[O>=O]+"), ID("[O,V]"), 0, "MO Ints (OO|OV)");
     dpd_buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[O,O]"), ID("[O,V]"), "MO Ints <OO|OV>");
     dpd_buf4_close(&K);
     timer_off("Sort (OO|OV) -> <OO|OV>");
     
     
     timer_on("Sort (OV|OV) -> <OO|VV>");
     // (OV|OV) -> <OO|VV>
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");
     dpd_buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[O,O]"), ID("[V,V]"), "MO Ints <OO|VV>");
     dpd_buf4_close(&K);
     timer_off("Sort (OV|OV) -> <OO|VV>");
     
    
     timer_on("Sort (OO|VV) -> <OV|OV>");
     // (OO|VV) -> <OV|OV>
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>=O]+"), ID("[V>=V]+"), 0, "MO Ints (OO|VV)");
     dpd_buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[O,V]"), ID("[O,V]"), "MO Ints <OV|OV>");
     dpd_buf4_close(&K);
     timer_off("Sort (OO|VV) -> <OV|OV>");
     
     
     timer_on("Sort (OV|VV) -> <OV|VV>");
     // (OV|VV) -> <OV|VV>
if (wfn_type_ == "OMP2" && incore_iabc_ == 0) { 
     tei_sort_iabc();
}

else {
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V>=V]+"), 0, "MO Ints (OV|VV)");
     dpd_buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[O,V]"), ID("[V,V]"), "MO Ints <OV|VV>");
     dpd_buf4_close(&K);
}
     timer_off("Sort (OV|VV) -> <OV|VV>");
     
       
//if (wfn_type_ == "OMP3" || wfn_type_ == "OCEPA" || wfn_type_ == "CEPA") { 
if (wfn_type_ != "OMP2") { 
     timer_on("Sort (VV|VV) -> <VV|VV>");
     // (VV|VV) -> <VV|VV>
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                 ID("[V>=V]+"), ID("[V>=V]+"), 0, "MO Ints (VV|VV)");
     dpd_buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[V,V]"), ID("[V,V]"), "MO Ints <VV|VV>");
     dpd_buf4_close(&K);
     timer_off("Sort (VV|VV) -> <VV|VV>");
}// end if (wfn_type_ == "OMP3" || wfn_type_ == "OCEPA") { 
     timer_off("Sort chem -> phys");
     
/********************************************************************************************/
/************************** Transform 1-electron int. to MO space ***************************/
/********************************************************************************************/        
      // Trans H matrix
      timer_on("Trans OEI");
      HmoA->copy(Hso);
      HmoA->transform(Ca_);
      timer_off("Trans OEI");
      
      if (print_ > 2) {
	HmoA->print();
      }
      
      // Trans Fock matrix    
      timer_on("Build Fock");
      fock_alpha();      
      timer_off("Build Fock");

      timer_on("Build Denominators");
      denominators_rhf();
      timer_off("Build Denominators");
      psio_->close(PSIF_LIBTRANS_DPD, 1);
      //fprintf(outfile,"\n trans_ints done. \n"); fflush(outfile);
 
}//


void OCCWave::denominators_rhf()
{
    //fprintf(outfile,"\n denominators is starting... \n"); fflush(outfile);
    dpdbuf4 D;
    dpdfile2 Fo,Fv;
    
    double *aOccEvals = new double [nacooA];
    double *aVirEvals = new double [nacvoA];
    
    // Pick out the diagonal elements of the Fock matrix, making sure that they are in the order
    // used by the DPD library, i.e. starting from zero for each space and ordering by irrep
    
    int aOccCount = 0, bOccCount = 0, aVirCount = 0, bVirCount = 0;
    
    //Diagonal elements of the Fock matrix
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0; i < aoccpiA[h]; ++i) aOccEvals[aOccCount++] = FockA->get(h, i + frzcpi_[h], i + frzcpi_[h]);
        for(int a = 0; a < avirtpiA[h]; ++a) aVirEvals[aVirCount++] = FockA->get(h, occpiA[h] + a, occpiA[h] + a); 
    }
    
    // Build denominators
    dpd_buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "D <OO|VV>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&D, h);
        for(int row = 0; row < D.params->rowtot[h]; ++row){
            int i = D.params->roworb[h][row][0];
            int j = D.params->roworb[h][row][1];
            for(int col = 0; col < D.params->coltot[h]; ++col){
                int a = D.params->colorb[h][col][0];
                int b = D.params->colorb[h][col][1];
                D.matrix[h][row][col] = 1.0/(aOccEvals[i] + aOccEvals[j] - aVirEvals[a] - aVirEvals[b]);
            }
        }
        dpd_buf4_mat_irrep_wrt(&D, h);
        dpd_buf4_mat_irrep_close(&D, h);
    }
    if (print_ > 2) dpd_buf4_print(&D, outfile, 1);
    dpd_buf4_close(&D);
    
    
    //Print
    if(print_ > 1){
      fprintf(outfile,"\n \n"); fflush(outfile);
      for(int i = 0; i<nacooA; i++) {
	fprintf(outfile,"\taOccEvals[%1d]: %20.14f\n", i, aOccEvals[i]); 
	fflush(outfile);
      }
      
      fprintf(outfile,"\n \n"); fflush(outfile);
      for(int i = 0; i<nacvoA; i++) {
	fprintf(outfile,"\taVirEvals[%1d]: %20.14f\n", i, aVirEvals[i]); 
	fflush(outfile);
      }
    }
    
    delete [] aOccEvals;
    delete [] aVirEvals;

 
    // Off-diagonal elements of the Fock matrix    
    // Build Occupied-Occupied block
    // The alpha-alpha spin case 
    dpd_file2_init(&Fo, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "F <O|O>");
    dpd_file2_mat_init(&Fo);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < aoccpiA[h]; ++i){
            for(int j = 0 ; j < aoccpiA[h]; ++j){
		if (i != j) Fo.matrix[h][i][j] = FockA->get(h, i + frzcpi_[h], j + frzcpi_[h]);
		else Fo.matrix[h][i][j] = 0.0;
            }
        }
    }
    dpd_file2_mat_wrt(&Fo);
    dpd_file2_close(&Fo);
    
    if (print_ > 2) {
      dpd_file2_init(&Fo, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "F <O|O>");
      dpd_file2_mat_init(&Fo);
      dpd_file2_mat_print(&Fo, outfile);
      dpd_file2_close(&Fo);
    }
    
    // Build Virtual-Virtual block
    // The alpha-alpha spin case 
    dpd_file2_init(&Fv, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "F <V|V>");
    dpd_file2_mat_init(&Fv);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < avirtpiA[h]; ++i){
            for(int j = 0 ; j < avirtpiA[h]; ++j){
                if (i != j) Fv.matrix[h][i][j] = FockA->get(h, i + occpiA[h], j + occpiA[h]);
		else Fv.matrix[h][i][j] = 0.0;
            }
        }
    }
    dpd_file2_mat_wrt(&Fv);
    dpd_file2_close(&Fv);
    
    if (print_ > 2) {
      dpd_file2_init(&Fv, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "F <V|V>");
      dpd_file2_mat_init(&Fv);
      dpd_file2_mat_print(&Fv, outfile);
      dpd_file2_close(&Fv);
    }

//fprintf(outfile,"\n denominators done. \n"); fflush(outfile);
}// end denominators
}} // End Namespaces


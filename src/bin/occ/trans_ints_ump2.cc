#include <libqt/qt.h>
#include <libtrans/integraltransform.h>

#include "occwave.h"
#include "defines.h"

using namespace boost;
using namespace psi;
using namespace std;

namespace psi{ namespace occwave{

void OCCWave::trans_ints_ump2()
{    
    //fprintf(outfile,"\n trans_ints is starting... \n"); fflush(outfile);
/********************************************************************************************/
/************************** Transform 2-electron int. to MO space ***************************/
/********************************************************************************************/  
    ints->update_orbitals();  
    ints->set_print(print_ - 2 >= 0 ? print_ - 2 : 0);
    ints->set_keep_dpd_so_ints(1);

    // Trans (OV|OV)
    timer_on("Trans (OV|OV)");
    ints->transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::occ, MOSpace::vir);
    timer_off("Trans (OV|OV)");
    
/********************************************************************************************/
/************************** sort chem -> phys ***********************************************/
/********************************************************************************************/  
     dpdbuf4 K, G;
     
     psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
     
     // Build MO ints    
     timer_on("Sort chem -> phys");
     // (OV|OV) -> <OO|VV>
     timer_on("Sort (OV|OV) -> <OO|VV>");
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");
     dpd_buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[O,O]"), ID("[V,V]"), "MO Ints <OO|VV>");
     dpd_buf4_close(&K);
     
     // (ov|ov) -> <oo|vv>
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints (ov|ov)");
     dpd_buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[o,o]"), ID("[v,v]"), "MO Ints <oo|vv>");
     dpd_buf4_close(&K);
     
     // (OV|ov) -> <Oo|Vv>
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "MO Ints (OV|ov)");
     dpd_buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[O,o]"), ID("[V,v]"), "MO Ints <Oo|Vv>");
     dpd_buf4_close(&K);
     
     // (OV|ov) -> <Ov|Vo>: <Ia||Bj> = <Ia|Bj> = (IB|ja)
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "MO Ints (OV|ov)");
     dpd_buf4_sort(&K, PSIF_LIBTRANS_DPD , psqr, ID("[O,v]"), ID("[V,o]"), "MO Ints <Ov|Vo>");
     dpd_buf4_close(&K);
     timer_off("Sort (OV|OV) -> <OO|VV>");
     timer_off("Sort chem -> phys");
     
/********************************************************************************************/
/************************** Antisymmetrized Ints ********************************************/
/********************************************************************************************/
      timer_on("Antisymmetrize integrals");
     // <OO||VV>:  <IJ||AB> = <IJ|AB> - <IJ|BA> 
     timer_on("Make <OO||VV>");
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
     dpd_buf4_copy(&K, PSIF_LIBTRANS_DPD, "MO Ints <OO||VV>");
     dpd_buf4_close(&K);
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO||VV>");
     dpd_buf4_init(&G, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
     for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&K, h);
	dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&K, h);
        dpd_buf4_mat_irrep_rd(&G, h);
        for(int ij = 0; ij < K.params->rowtot[h]; ++ij){
            for(int ab = 0; ab < K.params->coltot[h]; ++ab){
                int a = K.params->colorb[h][ab][0];
                int b = K.params->colorb[h][ab][1];
		int ba = G.params->colidx[b][a];
		K.matrix[h][ij][ab] -= G.matrix[h][ij][ba];                
            }
        }
        dpd_buf4_mat_irrep_wrt(&K, h);
        dpd_buf4_mat_irrep_close(&K, h);
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&K);
    dpd_buf4_close(&G);

     
     // <oo||vv>:  <ij||ab> = <ij|ab> - <ij|ba> 
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo|vv>");
     dpd_buf4_copy(&K, PSIF_LIBTRANS_DPD, "MO Ints <oo||vv>");
     dpd_buf4_close(&K);
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo||vv>");
     dpd_buf4_init(&G, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo|vv>");
     for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&K, h);
	dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&K, h);
        dpd_buf4_mat_irrep_rd(&G, h);
        for(int ij = 0; ij < K.params->rowtot[h]; ++ij){
            for(int ab = 0; ab < K.params->coltot[h]; ++ab){
                int a = K.params->colorb[h][ab][0];
                int b = K.params->colorb[h][ab][1];
		int ba = G.params->colidx[b][a];
		K.matrix[h][ij][ab] -= G.matrix[h][ij][ba];                
            }
        }
        dpd_buf4_mat_irrep_wrt(&K, h);
        dpd_buf4_mat_irrep_close(&K, h);
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&K);
    dpd_buf4_close(&G);
    timer_off("Make <OO||VV>");
    timer_off("Antisymmetrize integrals");

/********************************************************************************************/
/************************** Transform 1-electron int. to MO space ***************************/
/********************************************************************************************/        
      // Trans H matrix
      timer_on("Trans OEI");
      HmoA->copy(Hso);
      HmoB->copy(Hso);
      HmoA->transform(Ca_);
      HmoB->transform(Cb_);
      timer_off("Trans OEI");
      
      if (print_ > 1) {
	HmoA->print();
	HmoB->print();
      }
      
         for(int h = 0; h < nirrep_; ++h){
             for(int i = 0; i < occpiA[h]; ++i) FockA->set(h, i, i, epsilon_a_->get(h,i));
             for(int i = 0; i < occpiB[h]; ++i) FockB->set(h, i, i, epsilon_b_->get(h,i));
             for(int a = 0; a < virtpiA[h]; ++a) FockA->set(h, a + occpiA[h], a + occpiA[h], epsilon_a_->get(h, a + occpiA[h]));
             for(int a = 0; a < virtpiB[h]; ++a) FockB->set(h, a + occpiB[h], a + occpiB[h], epsilon_b_->get(h, a + occpiB[h]));
         }

      timer_on("Build Denominators");
      denominators_ump2();
      timer_off("Build Denominators");
      psio_->close(PSIF_LIBTRANS_DPD, 1);
      //fprintf(outfile,"\n trans_ints done. \n"); fflush(outfile);
 
}//


void OCCWave::denominators_ump2()
{
    //fprintf(outfile,"\n denominators is starting... \n"); fflush(outfile);
    dpdbuf4 D;
    dpdfile2 Fo,Fv;
    
    double *aOccEvals = new double [nacooA];
    double *bOccEvals = new double [nacooB];
    double *aVirEvals = new double [nacvoA];
    double *bVirEvals = new double [nacvoB];
    memset(aOccEvals, 0, sizeof(double)*nacooA);    
    memset(bOccEvals, 0, sizeof(double)*nacooB);    
    memset(aVirEvals, 0, sizeof(double)*nacvoA);    
    memset(bVirEvals, 0, sizeof(double)*nacvoB);    
    
    // Pick out the diagonal elements of the Fock matrix, making sure that they are in the order
    // used by the DPD library, i.e. starting from zero for each space and ordering by irrep
    
    int aOccCount = 0, bOccCount = 0, aVirCount = 0, bVirCount = 0;
    
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0; i < aoccpiA[h]; ++i) aOccEvals[aOccCount++] = epsilon_a_->get(h, i + frzcpi_[h]);
	for(int i = 0; i < aoccpiB[h]; ++i) bOccEvals[bOccCount++] = epsilon_b_->get(h, i + frzcpi_[h]);
        for(int a = 0; a < avirtpiA[h]; ++a) aVirEvals[aVirCount++] = epsilon_a_->get(h, a + occpiA[h]); 
	for(int a = 0; a < avirtpiB[h]; ++a) bVirEvals[bVirCount++] = epsilon_b_->get(h, a + occpiB[h]); 
    }
    
    // Build denominators
    // The alpha-alpha spin case
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
    

    // The beta-beta spin case 
    dpd_buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "D <oo|vv>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&D, h);
        for(int row = 0; row < D.params->rowtot[h]; ++row){
            int i = D.params->roworb[h][row][0];
            int j = D.params->roworb[h][row][1];
            for(int col = 0; col < D.params->coltot[h]; ++col){
                int a = D.params->colorb[h][col][0];
                int b = D.params->colorb[h][col][1];
                D.matrix[h][row][col] = 1.0/(bOccEvals[i] + bOccEvals[j] - bVirEvals[a] - bVirEvals[b]);
            }
        }
        dpd_buf4_mat_irrep_wrt(&D, h);
        dpd_buf4_mat_irrep_close(&D, h);
    }
    if (print_ > 2) dpd_buf4_print(&D, outfile, 1);
    dpd_buf4_close(&D);
    
    
    // The alpha-beta spin case 
    dpd_buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "D <Oo|Vv>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&D, h);
        for(int row = 0; row < D.params->rowtot[h]; ++row){
            int i = D.params->roworb[h][row][0];
            int j = D.params->roworb[h][row][1];
            for(int col = 0; col < D.params->coltot[h]; ++col){
                int a = D.params->colorb[h][col][0];
                int b = D.params->colorb[h][col][1];
                D.matrix[h][row][col] = 1.0/(aOccEvals[i] + bOccEvals[j] - aVirEvals[a] - bVirEvals[b]);
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
      for(int i = 0; i<nacooB; i++) {
	fprintf(outfile,"\tbOccEvals[%1d]: %20.14f\n", i, bOccEvals[i]); 
	fflush(outfile);
      }
      
      fprintf(outfile,"\n \n"); fflush(outfile);
      for(int i = 0; i<nacvoA; i++) {
	fprintf(outfile,"\taVirEvals[%1d]: %20.14f\n", i, aVirEvals[i]); 
	fflush(outfile);
      }
      
      fprintf(outfile,"\n \n"); fflush(outfile);
      for(int i = 0; i<nacvoB; i++) {
	fprintf(outfile,"\tbVirEvals[%1d]: %20.14f\n", i, bVirEvals[i]);
	fflush(outfile);
      }      
    }
    
    delete [] aOccEvals;
    delete [] bOccEvals;
    delete [] aVirEvals;
    delete [] bVirEvals;
 
//fprintf(outfile,"\n denominators done. \n"); fflush(outfile);
}// end denominators
}} // End Namespaces


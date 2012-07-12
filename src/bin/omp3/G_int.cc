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

/** Required PSI4 includes */ 
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

#include "omp3wave.h"
#include "defines.h"

using namespace boost;
using namespace psi;
using namespace std;

namespace psi{ namespace omp3wave{

void OMP3Wave::G_int()
{  
        //fprintf(outfile,"\n G_int is starting... \n"); fflush(outfile);
	GooA->zero();
	GooB->zero();
	GvvA->zero();
	GvvB->zero();

	dpdbuf4 T2_1AA, T2_1AB, T2_1BB, L2_1AA, L2_1AB, L2_1BB;
	dpdbuf4 T2_2AA, T2_2AB, T2_2BB, L2_2AA, L2_2AB, L2_2BB;
	dpdfile2 G;
	
	psio_->open(PSIF_OMP3_DPD, PSIO_OPEN_OLD);  
        psio_->open(PSIF_OMP3_DENSITY, PSIO_OPEN_OLD);
	
/********************************************************************************************/
/************************** Build G^(2) ints ************************************************/
/********************************************************************************************/
	// Open 1st order amplitude files
	dpd_buf4_init(&T2_1AA, PSIF_OMP3_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
	dpd_buf4_init(&T2_1BB, PSIF_OMP3_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2_1 <oo|vv>");
	dpd_buf4_init(&T2_1AB, PSIF_OMP3_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2_1 <Oo|Vv>");
	dpd_buf4_init(&L2_1AA, PSIF_OMP3_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
	dpd_buf4_init(&L2_1BB, PSIF_OMP3_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2_1 <oo|vv>");
	dpd_buf4_init(&L2_1AB, PSIF_OMP3_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2_1 <Oo|Vv>");
	
	// Open 2nd order amplitude files
	dpd_buf4_init(&T2_2AA, PSIF_OMP3_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_2 <OO|VV>");
	dpd_buf4_init(&T2_2BB, PSIF_OMP3_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2_2 <oo|vv>");
	dpd_buf4_init(&T2_2AB, PSIF_OMP3_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2_2 <Oo|Vv>");
	dpd_buf4_init(&L2_2AA, PSIF_OMP3_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_2 <OO|VV>");
	dpd_buf4_init(&L2_2BB, PSIF_OMP3_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2_2 <oo|vv>");
	dpd_buf4_init(&L2_2AB, PSIF_OMP3_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2_2 <Oo|Vv>");
	
	
	
	
	// Occupied-Occupied block
	// Alpha-Alpha spin case
	// G_IM = 1/2 \sum{N,E,F} t_IN^EF(1) * l_EF^MN(1) = 1/2 \sum{N,E,F} t_IN^EF(1) * t_MN^EF(1)
	dpd_file2_init(&G, PSIF_OMP3_DENSITY, 0, ID('O'), ID('O'), "G <O|O>");  
	dpd_contract442(&T2_1AA, &L2_1AA, &G, 0, 0, 0.5, 0.0);
	
	// G_IM += 1/2 \sum{N,E,F} t_IN^EF(2) * l_EF^MN(1) = 1/2 \sum{N,E,F} t_IN^EF(2) * t_MN^EF(1)
	dpd_contract442(&T2_2AA, &L2_1AA, &G, 0, 0, 0.5, 1.0);
	
	// G_IM += 1/2 \sum{N,E,F} t_IN^EF(1) * l_EF^MN(2) = 1/2 \sum{N,E,F} t_IN^EF(1) * t_MN^EF(2)
	dpd_contract442(&T2_1AA, &L2_2AA, &G, 0, 0, 0.5, 1.0);
	
	
	// G_IM += \sum{n,E,f} t_In^Ef(1) * l_Ef^Mn(1) = \sum{N,E,F} t_In^Ef(1) * t_Mn^Ef(1)
	dpd_contract442(&T2_1AB, &L2_1AB, &G, 0, 0, 1.0, 1.0);
	
	// G_IM += \sum{n,E,f} t_In^Ef(2) * l_Ef^Mn(1) = \sum{N,E,F} t_In^Ef(2) * t_Mn^Ef(1)
	dpd_contract442(&T2_2AB, &L2_1AB, &G, 0, 0, 1.0, 1.0);
	
	// G_IM += \sum{n,E,f} t_In^Ef(1) * l_Ef^Mn(2) = \sum{N,E,F} t_In^Ef(1) * t_Mn^Ef(2)
	dpd_contract442(&T2_1AB, &L2_2AB, &G, 0, 0, 1.0, 1.0);
	dpd_file2_close(&G);
	
	
	// Beta-Beta spin case
	// G_im = 1/2 \sum{n,e,f} t_in^ef(1) * l_ef^mn(1) = 1/2 \sum{n,e,f} t_in^ef(1) * t_mn^ef(1)
	dpd_file2_init(&G, PSIF_OMP3_DENSITY, 0, ID('o'), ID('o'), "G <o|o>");  
	dpd_contract442(&T2_1BB, &L2_1BB, &G, 0, 0, 0.5, 0.0);
	
	// G_im += 1/2 \sum{n,e,f} t_in^ef(2) * l_ef^mn(1) = 1/2 \sum{n,e,f} t_in^ef(2) * t_mn^ef(1)
	dpd_contract442(&T2_2BB, &L2_1BB, &G, 0, 0, 0.5, 1.0);
	
	// G_im += 1/2 \sum{n,e,f} t_in^ef(1) * l_ef^mn(2) = 1/2 \sum{n,e,f} t_in^ef(1) * t_mn^ef(2)
	dpd_contract442(&T2_1BB, &L2_2BB, &G, 0, 0, 0.5, 1.0);
	
	// G_im  += \sum{N,e,F} t_Ni^Fe(1) * l_Fe^Nm(1) = \sum{N,e,F} t_Ni^Fe(1) * t_Nm^Fe(1)
	dpd_contract442(&T2_1AB, &L2_1AB, &G, 1, 1, 1.0, 1.0);
	
	// G_im  += \sum{N,e,F} t_Ni^Fe(2) * l_Fe^Nm(1) = \sum{N,e,F} t_Ni^Fe(2) * t_Nm^Fe(1)
	dpd_contract442(&T2_2AB, &L2_1AB, &G, 1, 1, 1.0, 1.0);
	
	// G_im  += \sum{N,e,F} t_Ni^Fe(1) * l_Fe^Nm(2) = \sum{N,e,F} t_Ni^Fe(1) * t_Nm^Fe(2)
	dpd_contract442(&T2_1AB, &L2_2AB, &G, 1, 1, 1.0, 1.0);
	dpd_file2_close(&G);
	
	
	
	
	
	// Virtual-Virtual block
	// Alpha-Alpha spin case
	// G_EA = -1/2 \sum{M,N,F} t_MN^AF(1) * l_EF^MN(1) = -1/2 \sum{M,N,F} t_MN^AF(1) * t_MN^EF(1)
	dpd_file2_init(&G, PSIF_OMP3_DENSITY, 0, ID('V'), ID('V'), "G <V|V>");  
	dpd_contract442(&T2_1AA, &L2_1AA, &G, 2, 2, -0.5, 0.0); 
	
	// G_EA += -1/2 \sum{M,N,F} t_MN^AF(2) * l_EF^MN(1) = -1/2 \sum{M,N,F} t_MN^AF(2) * t_MN^EF(1)  
	dpd_contract442(&T2_2AA, &L2_1AA, &G, 2, 2, -0.5, 1.0); 
	
	// G_EA += -1/2 \sum{M,N,F} t_MN^AF(1) * l_EF^MN(2) = -1/2 \sum{M,N,F} t_MN^AF(1) * t_MN^EF(2)
	dpd_contract442(&T2_1AA, &L2_2AA, &G, 2, 2, -0.5, 1.0); 
	
	// G_EA += - \sum{M,n,f} t_Mn^Af(1) * l_Ef^Mn(1) = - \sum{M,n,f} t_Mn^Af(1) * t_Mn^Ef(1)
	dpd_contract442(&T2_1AB, &L2_1AB, &G, 2, 2, -1.0, 1.0); 
	
	// G_EA += - \sum{M,n,f} t_Mn^Af(2) * l_Ef^Mn(1) = - \sum{M,n,f} t_Mn^Af(2) * t_Mn^Ef(1)
	dpd_contract442(&T2_2AB, &L2_1AB, &G, 2, 2, -1.0, 1.0); 
	
	// G_EA += - \sum{M,n,f} t_Mn^Af(1) * l_Ef^Mn(2) = - \sum{M,n,f} t_Mn^Af(1) * t_Mn^Ef(2)
	dpd_contract442(&T2_1AB, &L2_2AB, &G, 2, 2, -1.0, 1.0); 
	dpd_file2_close(&G);
	

	// Beta-Beta spin case
	// G_ea = -1/2 \sum{m,n,f} t_mn^af(1) * l_ef^mn(1) = -1/2 \sum{m,n,f} t_mn^af(1) * t_mn^ef(1)
	dpd_file2_init(&G, PSIF_OMP3_DENSITY, 0, ID('v'), ID('v'), "G <v|v>");  
	dpd_contract442(&T2_1BB, &L2_1BB, &G, 2, 2, -0.5, 0.0); 
	
	// G_ea += -1/2 \sum{m,n,f} t_mn^af(2) * l_ef^mn(1) = -1/2 \sum{m,n,f} t_mn^af(2) * t_mn^ef(1) 
	dpd_contract442(&T2_2BB, &L2_1BB, &G, 2, 2, -0.5, 1.0); 
	
	// G_ea += -1/2 \sum{m,n,f} t_mn^af(1) * l_ef^mn(2) = -1/2 \sum{m,n,f} t_mn^af(1) * t_mn^ef(2)
	dpd_contract442(&T2_1BB, &L2_2BB, &G, 2, 2, -0.5, 1.0); 
	
	// G_ea += - \sum{M,n,F} t_Mn^Fa(1) * l_Fe^Mn(1) = - \sum{M,n,F} t_Mn^Fa(1) * t_Mn^Fe(1)
	dpd_contract442(&T2_1AB, &L2_1AB, &G, 3, 3, -1.0, 1.0); 
	
	// G_ea += - \sum{M,n,F} t_Mn^Fa(2) * l_Fe^Mn(1) = - \sum{M,n,F} t_Mn^Fa(2) * t_Mn^Fe(1)
	dpd_contract442(&T2_2AB, &L2_1AB, &G, 3, 3, -1.0, 1.0); 
	
	// G_ea += - \sum{M,n,F} t_Mn^Fa(1) * l_Fe^Mn(2) = - \sum{M,n,F} t_Mn^Fa(1) * t_Mn^Fe(2)
	dpd_contract442(&T2_1AB, &L2_2AB, &G, 3, 3, -1.0, 1.0); 
	dpd_file2_close(&G);
	
	
	// Close amplitude files
	dpd_buf4_close(&T2_1AA);
	dpd_buf4_close(&T2_1BB);
	dpd_buf4_close(&T2_1AB);
	dpd_buf4_close(&L2_1AA);
	dpd_buf4_close(&L2_1BB);
	dpd_buf4_close(&L2_1AB);
	
	dpd_buf4_close(&T2_2AA);
	dpd_buf4_close(&T2_2BB);
	dpd_buf4_close(&T2_2AB);
	dpd_buf4_close(&L2_2AA);
	dpd_buf4_close(&L2_2BB);
	dpd_buf4_close(&L2_2AB);
	
	
/********************************************************************************************/
/************************** Load dpd_file2 to Matrix ****************************************/
/********************************************************************************************/
	// Load dpd_file2 to Matrix (Goo)
	// Alpha-Alpha spin case
	dpd_file2_init(&G, PSIF_OMP3_DENSITY, 0, ID('O'), ID('O'), "G <O|O>");  
	dpd_file2_mat_init(&G);
	dpd_file2_mat_rd(&G);
	for(int h = 0; h < nirreps; ++h){
	  for(int i = 0 ; i < aoccpiA[h]; ++i){
            for(int j = 0 ; j < aoccpiA[h]; ++j){
                GooA->set(h, i, j, G.matrix[h][i][j]);
            }
	  }
	}
	dpd_file2_close(&G);
	
	// Beta-Beta spin case
	dpd_file2_init(&G, PSIF_OMP3_DENSITY, 0, ID('o'), ID('o'), "G <o|o>");  
	dpd_file2_mat_init(&G);
	dpd_file2_mat_rd(&G);
	for(int h = 0; h < nirreps; ++h){
	  for(int i = 0 ; i < aoccpiB[h]; ++i){
            for(int j = 0 ; j < aoccpiB[h]; ++j){
                GooB->set(h, i, j, G.matrix[h][i][j]);
            }
	  }
	}
	dpd_file2_close(&G);
	
	
	
	// Load dpd_file2 to Matrix (Gvv)
	// Alpha-Alpha spin case
	dpd_file2_init(&G, PSIF_OMP3_DENSITY, 0, ID('V'), ID('V'), "G <V|V>"); 
	dpd_file2_mat_init(&G);
	dpd_file2_mat_rd(&G);
	for(int h = 0; h < nirreps; ++h){
	  for(int i = 0 ; i < avirtpiA[h]; ++i){
            for(int j = 0 ; j < avirtpiA[h]; ++j){
                GvvA->set(h, i, j, G.matrix[h][i][j]);
            }
	  }
	}
	dpd_file2_close(&G);
	
	// Beta-Beta spin case
	dpd_file2_init(&G, PSIF_OMP3_DENSITY, 0, ID('v'), ID('v'), "G <v|v>");  
	dpd_file2_mat_init(&G);
	dpd_file2_mat_rd(&G);
	for(int h = 0; h < nirreps; ++h){
	  for(int i = 0 ; i < avirtpiB[h]; ++i){
            for(int j = 0 ; j < avirtpiB[h]; ++j){
                GvvB->set(h, i, j, G.matrix[h][i][j]);
            }
	  }
	}
	dpd_file2_close(&G);
	
	psio_->close(PSIF_OMP3_DPD, 1);  
        psio_->close(PSIF_OMP3_DENSITY, 1);

	if (print_ > 1) {
	  GooA->print();
	  GooB->print();
	  GvvA->print();
	  GvvB->print();
	}
	
	//printf(outfile,"\n G_int done... \n"); fflush(outfile);

} // end of G_int
}} // End Namespaces


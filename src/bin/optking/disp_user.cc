/*! \file
    \ingroup OPTKING
    \brief DISP_USER only performs input-specified displacements
*/

#define EXTERN
#include "globals.h"
#undef EXTERN
#include "cartesians.h"
#include "simples.h"
#include "salc.h"
#include "opt.h"

#include <libipv1/ip_lib.h>

namespace psi { namespace optking {

int disp_user(const cartesians &carts, simples_class & simples, const salc_set &all_salcs) {
  int i,j,a,b,success,h;
  int  num_disps = 0, disp_length = 0, restart_geom_file, line_length_count;
  double *geom, *djunk, *dq, **displacements, disp = 0;
  char *disp_label, *ch_temp, *salc_lbl;
  int *irrep_per_disp;

  disp_label = new char[MAX_LINELENGTH];
  //djunk = new double[3*carts.get_natom()];

  dq = init_array(all_salcs.get_num());

  int *irrep_salcs;
  irrep_salcs = init_int_array(all_salcs.get_num());
  for (i=0; i<all_salcs.get_num(); ++i) {
    salc_lbl = all_salcs.get_label(i);
    for (h=0; h<syminfo.nirreps; ++h) {
      if (!strcmp(salc_lbl, syminfo.irrep_lbls[h]))
        irrep_salcs[i] = h;
     }
  }

  ip_count("DISPLACEMENTS",&num_disps,0);
  irrep_per_disp = init_int_array(num_disps);

  if (num_disps < 1) {
    punt("No DISPLACEMENTS vector found in input.");
  }

  displacements = init_matrix(num_disps,all_salcs.get_num());
  for (i=0;i<num_disps;++i) {
    disp_length = 0;
    ip_count("DISPLACEMENTS",&disp_length,1,i);
    if (div_int(disp_length,2) && disp_length != 0) {
      for (j=0;j<disp_length;j+=2) {
        ip_data("DISPLACEMENTS","%d",&a,2,i,j);
        ip_data("DISPLACEMENTS","%lf",&disp,2,i,j+1);
        if ((a>0) && (a<=all_salcs.get_num())) {
          displacements[i][a-1] = disp; 
          irrep_per_disp[i] = irrep_salcs[a-1];
        }
        else {
          punt("internal to be displaced does not exist");
        }
      }
    }
    else {
      punt("displacement vector has wrong number of elements");
      exit(1);
    }
  }

  free_int_array(irrep_salcs);

  fprintf(outfile,"Displacement Matrix\n");
  print_mat5(displacements, num_disps, all_salcs.get_num(), outfile);


  /*** generate and store Micro_iteration cartesian geometries ***/
  double **micro_geoms;
  micro_geoms = init_matrix(num_disps, 3*carts.get_natom());
  for (i=0;i<num_disps;++i)  {
    //sprintf(disp_label,"Displaced geometry %d in a.u.\n",i+1);

    // print to text file
    sprintf(disp_label, "Disp:");
    ch_temp = disp_label + 5 ;
    line_length_count = 5;
    for (j=0; j<all_salcs.get_num(); ++j) {
      if (fabs(displacements[i][j]) > 1.0E-8) {
        if (line_length_count < (MAX_LINELENGTH - 10)) {
          sprintf(ch_temp, " (%d %5.3lf)", j+1, displacements[i][j]);
          ch_temp += 10;
          line_length_count += 10;
        }
      }
    }

    success = new_geom(carts,simples,all_salcs,displacements[i],PRINT_TO_GEOM,
        0, disp_label, i, 0, micro_geoms[i]);

    //
    if (!success) {
      fprintf(outfile,"Unable to generate displaced geometry.\n");
      fclose(outfile);
      exit(PSI_RETURN_FAILURE);
    }
  }
  free_matrix(displacements);

  open_PSIF();
  psio_write_entry(PSIF_OPTKING, "OPT: Displaced geometries",
      (char *) &(micro_geoms[0][0]), num_disps*3*carts.get_natom()*sizeof(double));

  optinfo.disp_num = 0;
  psio_write_entry(PSIF_OPTKING, "OPT: Current disp_num",
      (char *) &(optinfo.disp_num),sizeof(int));

  psio_write_entry(PSIF_OPTKING, "OPT: Num. of disp.",
      (char *) &(num_disps), sizeof(int));

  psio_write_entry(PSIF_OPTKING, "OPT: Irrep per disp",
    (char *) &(irrep_per_disp[0]), num_disps*sizeof(int));
  free_int_array(irrep_per_disp);

  double *disp_e;
  disp_e = init_array(num_disps);
  psio_write_entry(PSIF_OPTKING, "OPT: Displaced energies",
      (char *) &(disp_e[0]), num_disps*sizeof(double));
  free_array(disp_e);

  close_PSIF();

  free_matrix(micro_geoms);
  free(disp_label);
  //free(djunk);

  //return num_disps;
  return 100;
}

}} /* namespace psi::optking */


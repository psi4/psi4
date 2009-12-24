/*! \file
    \ingroup OPTKING
    \brief GRAD_ENERGY computes a file11 entry from energies in chkpt Rollin King, 2002
*/

#define EXTERN
#include "globals.h"
#undef EXTERN
#include "cartesians.h"
#include "simples.h"
#include "salc.h"
#include "opt.h"

#include <libchkpt/chkpt.h>
#include <libipv1/ip_lib.h>
#include <libpsio/psio.h>

namespace psi { namespace optking {

void grad_energy(cartesians &carts, simples_class &simples, const salc_set &symm) {

  int i,j,a,b, dim, dim_carts, num_disps, cnt, natom;
  double **B, *geom, *forces;
  double energy, *energies, **micro_geoms, **displacements;
  double *f, *f_q, *dq, *q, tval, **geom2D;
  char *disp_label, *line1;
  FILE *fp_energy_dat;

  disp_label = new char[MAX_LINELENGTH];
  natom = carts.get_natom();
  dim_carts = 3*carts.get_natom();
  natom = carts.get_natom();

  if (symm.get_num() == 0) {
    punt("No symmetric internal coordinates to optimize.\n");
  }

  if (optinfo.points == 3)
    num_disps = 2 * symm.get_num();
  else if (optinfo.points == 5)
    num_disps = 4 * symm.get_num();

  /* read in energies */
  energies = new double[num_disps];
  for (i=0;i<num_disps;++i) energies[i] = 0;
  if (optinfo.energy_dat) { /* read from energy.dat text file */
    fp_energy_dat = fopen("energy.dat", "r");
    rewind (fp_energy_dat);
    line1 = new char[MAX_LINELENGTH+1];
    // ACS (11/06) Allow external program to be used to compute energies
    if(optinfo.external_energies){
      /* Read the first energy and dump it as the reference energy */
      double temp;
      fscanf(fp_energy_dat,"%lf",&temp);
      open_PSIF();
      psio_write_entry(PSIF_OPTKING, "OPT: Reference energy",(char *) &(temp), sizeof(double));
      close_PSIF();
    }
    for (i=0; i<num_disps; i++) {
      fscanf(fp_energy_dat,"%lf",&(energies[i]));
    }
    fclose(fp_energy_dat);
    delete [] line1;
  }
  else { /* read from checkpoint file (default) */
    open_PSIF();
    psio_read_entry(PSIF_OPTKING, "OPT: Displaced energies",
        (char *) &(energies[0]), num_disps*sizeof(double));
    close_PSIF();
  }

  fprintf(outfile,"Energies of displaced geometries. Check for precision!\n");
  cnt = -1;
  for (i=0;i<symm.get_num();++i) {
    fprintf(outfile,"Coordinate %d: ",i);
    for (j=0;j<optinfo.points-1;++j)
      fprintf(outfile,"%15.10lf ",energies[++cnt]);
    fprintf(outfile,"\n");
  }
  fflush(outfile);

  // Calculate forces in internal coordinates
  f_q = new double[symm.get_num()];
  if (optinfo.points == 3) {
    for (i=0;i<symm.get_num();++i) {
      f_q[i] = (energies[2*i+1]-energies[2*i]) / (2.0 * optinfo.disp_size);
      f_q[i] = -1.0 * f_q[i] * _hartree2J * 1.0E18 ;
    }
  }
  else if (optinfo.points == 5) {
    for (i=0;i<symm.get_num();++i) {
      f_q[i] = ( energies[4*i]-8.0*energies[4*i+1]+8.0*energies[4*i+2]-energies[4*i+3])
                  / (12.0 * optinfo.disp_size);
      f_q[i] = -1.0 * f_q[i] * _hartree2J * 1.0E18 ;
    }
  }
  free_array(energies);

  // Print internal coordinate forces
fprintf(outfile,"\nInternal coordinate forces\n");
for (i=0;i<symm.get_num();++i)
   fprintf(outfile,"%13.10lf\n",f_q[i]);

  // write out approximate file11.dat
  geom = new double[dim_carts];
  open_PSIF();
  psio_read_entry(PSIF_OPTKING, "OPT: Reference geometry",
      (char *) &(geom[0]), dim_carts*sizeof(double));
  psio_read_entry(PSIF_OPTKING, "OPT: Reference energy",
      (char *) &(energy), sizeof(double));
  close_PSIF();

  // Transform forces to cartesian coordinates
  simples.compute(geom);
  simples.compute_s(geom);
  B = compute_B(simples, symm);
  f = new double[dim_carts];
  opt_mmult(B,1,&f_q,1,&f,1,dim_carts,symm.get_num(),1,0);
  free_matrix(B);

  // change forces to gradient for writing a file11 entry
  for(i=0;i<dim_carts;++i)
    f[i] = -1.0 * f[i] / _hartree2J / 1.0E18 * _bohr2angstroms;

  opt_ffile(&fp_11, "file11.dat", 1);
  char *wfn,*dertype;
  int errcod = ip_string("WFN", &wfn, 0);
  if (errcod != IPE_OK)
    punt("Keyword WFN not found in input file");
  errcod = ip_string("DERTYPE", &dertype, 0);
  if (errcod != IPE_OK)
    dertype = strdup("NONE");
  chkpt_init(PSIO_OPEN_OLD);
  char* label = chkpt_rd_label();
  chkpt_close();
  sprintf(disp_label,"%-59.59s %-10.10s%-8.8s",label,wfn,dertype);
  free(label); free(wfn); free(dertype);
  carts.set_energy(energy);
  carts.set_coord(geom);
  carts.set_grad(f);
  carts.print(11,fp_11,0,disp_label, 0);
  fclose(fp_11);

  // write out geometry, gradient and energy to chkpt file
  cnt = -1;
  geom2D = init_matrix(carts.get_natom(),3);
  for (i=0; i<carts.get_natom(); ++i)
    for (j=0; j<3; ++j)
      geom2D[i][j] = geom[++cnt];

  chkpt_init(PSIO_OPEN_OLD);
  chkpt_wt_geom(geom2D);
  chkpt_wt_grad(f);
  chkpt_wt_etot(energy);
  chkpt_close();
  free_matrix(geom2D);
  free_array(f);
  free_array(geom);

  // recompute values of internals and s vectors -- too late!
//  simples.compute(carts.get_coord());
//  simples.compute_s(carts.get_coord() );

  // use optking --opt_step to take a step
  // opt_step(carts, simples, symm);

  // reset microiteration value in disp_all
  /*
  optinfo.micro_iteration = 0;
  open_PSIF();
  psio_write_entry(PSIF_OPTKING, "Micro_iteration",
      (char *) &(optinfo.micro_iteration),sizeof(int));
  close_PSIF();
  */
}

}} /* namespace psi::optking */


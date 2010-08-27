/*! \file optking.cc
    \ingroup opt
    \brief optimization */

#include <cstdio>
#include <fstream>
#include "globals.h"
#include "molecule.h"
#include "print.h"

namespace opt {
  void print_title(FILE *fp);
  void print_end(FILE *fp);
  void set_params(void); // set optimization parameters

  //Package = PSI; // designate package

#ifdef PSI4
using psi::Options;
OptReturnType optking(Options & options) {
#elif "QCHEM4"
OptReturnType optking(void) {
outfile = fopen("output.dat","a");
#endif

  using namespace std;

  // open output file and write title header
  print_title(outfile);

  // set optimization parameters
  set_params();

  // read in geometry, gradient, energy from text input file
  ifstream fin;
  fin.open(FILENAME_GEOM_GRAD_IN, ios_base::in);
  MOLECULE mol1(fin);

  // automatically generate coordinates
  mol1.update_connectivity_by_distances();
  mol1.add_simples_by_connectivity();
  mol1.compute_intco_values();

  // print to output file
  mol1.print_connectivity(outfile);
  mol1.print_geom_grad(outfile);
  mol1.print_intcos(outfile);

  // read binary file for previous step data
  p_Opt_data = new OPT_DATA(mol1.g_nintco(), 3*mol1.g_natom());

  // save geometry, gradient, and energy
  double * x = mol1.g_geom_array();
  p_Opt_data->save_geom_energy(x, mol1.g_energy());
  free_array(x);

  // compute forces in internal coordinates from cartesian gradient
  mol1.forces(); // puts forces in p_Opt_data->step[last one]

  if (p_Opt_data->g_iteration() == 1)
    mol1.H_guess(); // puts Hessian guess in p_Opt_data->step[last one]
  else
    p_Opt_data->H_update(mol1);

  mol1.project_f_and_H();

  if (Opt_params.step_type == opt::OPT_PARAMS::NR)
    mol1.nr_step(); // puts dq in p_Opt_data->step
  else if (Opt_params.step_type == opt::OPT_PARAMS::RFO)
    mol1.rfo_step(); // puts dq in p_Opt_data->step

  if ( p_Opt_data->conv_check() ) {
    printf("\t *** Geometry is converged!\n");
    fprintf(outfile,"\t *** Geometry is converged!\n");
    p_Opt_data->summary();
    p_Opt_data->write();
#ifdef PSI4
  mol1.write_geom_chkpt(); // write geometry for next step to chkpt file
#endif
    delete p_Opt_data;
    print_end(outfile);
#ifdef QCHEM
    fclose(outfile);
#endif
    return OptReturnEndloop;
  }

  p_Opt_data->write();
  delete p_Opt_data;

  mol1.write_geom(); // write geometry for next step to output file

#ifdef PSI4
  mol1.write_geom_chkpt(); // write geometry for next step to chkpt file
#endif

  print_end(outfile);
#ifdef QCHEM
  fclose(outfile);
#endif
  return OptReturnSuccess;
}

void print_title(FILE *fp) {
  fprintf(fp, "\n\t\t\t----------------------------------\n");
  fprintf(fp, "\t\t\t OPTKING: for geometry optimizations  \n");
  fprintf(fp, "\t\t\t  - R.A. King,  Bethel University   \n");
  fprintf(fp, "\t\t\t------------------------------------\n");
}
void print_end(FILE *fp) {
  fprintf(fp, "\t\t\t--------------------------\n");
  fprintf(fp, "\t\t\t OPTKING Finished Execution \n");
  fprintf(fp, "\t\t\t--------------------------\n");
}
}


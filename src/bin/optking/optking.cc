/*! \file optking.cc
    \ingroup opt
    \brief optimization */

#include <cstdio>
#include <fstream>
#include "globals.h"
#include "molecule.h"
#include "print.h"

#include "io.h"

// Define the return types for optking.
#if defined(OPTKING_PACKAGE_PSI)
  typedef psi::PsiReturnType OptReturnType;
  #define OptReturnEndloop (psi::Endloop)
  #define OptReturnSuccess (psi::Success)
  #define OptReturnFailure (psi::Failure)
#elif QCHEM4
  typedef int OptReturnType;
  #define OptReturnEndloop 1
  #define OptReturnSuccess 0
#endif

namespace opt {
  void open_output_dat(void); // open/link outfile to text output
  void close_output_dat(void);// close above
  void print_title(FILE *fp); // print header
  void print_end(FILE *fp);   // print footer
  void set_params(void);      // set optimization parameters
  void init_ioff(void);

#if defined(OPTKING_PACKAGE_PSI)
OptReturnType optking(psi::Options & options) {
#elif defined(QCHEM4)
OptReturnType optking(void) {
#endif

  open_output_dat();   // assign output.dat file pointer

  using namespace std;
  MOLECULE *mol1;
  bool newly_generated_coordinates; // Are internal coordinates produced or read-in?

  print_title(outfile); // print header
  //init_ioff(); // not used at present

  set_params(); // set optimization parameters

  // try to open old internal coordinates
  std::ifstream if_intco(FILENAME_INTCO_DAT, ios_base::in);

  if (if_intco.is_open()) { // old internal coordinates are present

    fprintf(outfile,"\n\tPrevous internal coordinate definitions found.\n");
    newly_generated_coordinates = false;

    mol1 = new MOLECULE(0);    // make an empty molecule

    // read internal coordinate and fragment definitions
    // create, allocate, and add fragment objects
    mol1->read_intcos(if_intco);
    if_intco.close();

    mol1->update_connectivity_by_bonds();

    // read geometry and gradient into the existing fragments
    mol1->read_geom_grad();

  }
  else { // automatically generate coordinates

    fprintf(outfile,"\n\tInternal coordinates to be generated automatically.\n");
    newly_generated_coordinates = true;

    // read number of atoms ; make one fragment of that size ;
    mol1 = new MOLECULE( read_natoms() );

    // read geometry and gradient into fragment
    mol1->read_geom_grad();

    // use covalent radii to define bonds
    mol1->update_connectivity_by_distances();

    // if fragment_mode == SINGLE, connects all separated groups of atoms by modifying frag.connectivity
    // if fragment_mode == MULTI, splits into fragments and makes interfragment coordinates
    mol1->fragmentize();

    if (Opt_params.fragment_mode == OPT_PARAMS::SINGLE) {
      mol1->add_intrafragment_simples_by_connectivity();
    }
      
    // newly constructed fragments need connectivity generated
    if (Opt_params.fragment_mode == OPT_PARAMS::MULTI) {
      mol1->update_connectivity_by_distances();
      mol1->add_intrafragment_simples_by_connectivity();
      mol1->add_interfragment();
    }

    // print out internal coordinates for future steps
    FILE *fp_intco = fopen(FILENAME_INTCO_DAT, "w");
    mol1->print_intco_dat(fp_intco);
    fclose(fp_intco);

    // only generate coordinates and print them out
    if (Opt_params.generate_intcos_only)
      return OptReturnSuccess;

  }

  // print geometry and gradient
  mol1->print_geom_grad(outfile);

mol1->print_connectivity(outfile); fflush(outfile);

  // print internal coordinate definitions and values
  mol1->print_intcos(outfile);
  fflush(outfile);

  if (Opt_params.test_B)
    mol1->test_B();

  if (Opt_params.test_derivative_B)
    mol1->test_derivative_B();

  // read binary file for previous steps
  p_Opt_data = new OPT_DATA(mol1->g_nintco(), 3*mol1->g_natom());

  // save geometry and energy
  double * x = mol1->g_geom_array();
  p_Opt_data->save_geom_energy(x, mol1->g_energy());
  free_array(x);

  p_Opt_data->previous_step_report();

  //double **H_int = mol1->cartesian_H_to_internals();
  //free_matrix(H_int);

  // compute forces in internal coordinates from cartesian gradient
  mol1->forces(); // puts forces in p_Opt_data->step[last one]

  if (p_Opt_data->g_iteration() == 1)
    mol1->H_guess(); // puts Hessian guess in p_Opt_data->step[last one]
  else {
    try {
      p_Opt_data->H_update(*mol1);
    } catch (const char * str) {
      fprintf(stderr, "%s\n", str);
      fprintf(outfile, "H_update failed\n");
      fprintf(outfile, "%s\n", str);
      return OptReturnFailure;
    }
  }

  mol1->project_f_and_H();

  if (Opt_params.step_type == opt::OPT_PARAMS::NR)
    mol1->nr_step(); // puts dq in p_Opt_data->step
  else if (Opt_params.step_type == opt::OPT_PARAMS::RFO)
    mol1->rfo_step(); // puts dq in p_Opt_data->step

  if ( p_Opt_data->conv_check() ) {
    printf("\t *** Geometry is converged!\n");
    fprintf(outfile,"\t *** Geometry is converged!\n");
    p_Opt_data->summary();
    p_Opt_data->write();
    mol1->write_geom(); // geometry to chkpt in PSI
    delete p_Opt_data;
    fprintf(outfile,"\tFinal (next step) structure:\n");
    mol1->print_geom(); // geometry to output file
    print_end(outfile);
    close_output_dat();
    return OptReturnEndloop;
  }

  p_Opt_data->write();
  delete p_Opt_data;

  mol1->write_geom(); // write geometry to package

  fprintf(outfile,"\tStructure for next step:\n"); fflush(outfile);
  mol1->print_geom(); // write geometry for next step to output file

  print_end(outfile);
  close_output_dat();
  //free_int_array(ioff); // not used at present
  return OptReturnSuccess;
}

// Functions to set and release (possibly open and close) the file pointer for
// the standard text output file
void open_output_dat(void) {
#if defined (OPTKING_PACKAGE_PSI)
  opt::outfile = psi::outfile;
#elif defined (OPTKING_PACKAGE_QCHEM)
  opt::outfile = fopen(FILENAME_OUTPUT_DAT,"a");
#endif
}

void close_output_dat(void) {
#if defined(OPTKING_PACKAGE_QCHEM)
  fclose(opt::outfile);
#endif
}

void print_title(FILE *fp) {
  fprintf(fp, "\n\t\t\t----------------------------------\n");
  fprintf(fp, "\t\t\t OPTKING: for geometry optimizations  \n");
  fprintf(fp, "\t\t\t  - R.A. King,  Bethel University   \n");
  fprintf(fp, "\t\t\t------------------------------------\n");
  fflush(outfile);
}

void print_end(FILE *fp) {
  fprintf(fp, "\t\t\t--------------------------\n");
  fprintf(fp, "\t\t\t OPTKING Finished Execution \n");
  fprintf(fp, "\t\t\t--------------------------\n");
  fflush(outfile);
}

void init_ioff(void)
{
  int i;
  ioff = init_int_array(IOFF_MAX);
  ioff[0] = 0;
  for(i=1; i < IOFF_MAX; i++) ioff[i] = ioff[i-1] + i;
}

}


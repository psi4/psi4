/*! 
    \defgroup optking
    \file optking.cc : main optimizer
    \ingroup optking
*/
#include <cstdio>
#include <fstream>
#include "globals.h"
#include "molecule.h"
#include "print.h"
#include "io.h"
#include <libparallel/parallel.h>


// Define the return types for optking.
#if defined(OPTKING_PACKAGE_PSI)
  typedef psi::PsiReturnType OptReturnType;
  #define OptReturnEndloop (psi::EndLoop)
  #define OptReturnSuccess (psi::Success)
  #define OptReturnFailure (psi::Failure)
#elif defined(OPTKING_PACKAGE_QCHEM)
  typedef int OptReturnType;
  #define OptReturnFailure 2
  #define OptReturnEndloop 1
  #define OptReturnSuccess 0
#endif

#if defined(OPTKING_PACKAGE_PSI)
namespace psi { void psiclean(void); }
#endif

#if defined(OPTKING_PACKAGE_QCHEM)
 #include "qchem.h"
 namespace opt {
   OptReturnType optking(void); // declare optking
 }
 OptReturnType opt2man_main(void) { opt::optking(); } // QCHEM wrapper/alias
#endif

namespace opt {
  void open_output_dat(void); // open/link outfile to text output
  void close_output_dat(void);// close above
  void print_title(void); // print header
  void print_end(void);   // print footer
  void init_ioff(void);
  bool INTCO_EXCEPT::already_tried_other_intcos = false;
  bool INTCO_EXCEPT::override_fragment_mode = false; // to override MULTI setting by exception algorithm

#if defined(OPTKING_PACKAGE_PSI)
  void set_params(psi::Options & options);      // set optimization parameters
#else
  void set_params(void);
#endif

#if defined(OPTKING_PACKAGE_PSI)
OptReturnType optking(psi::Options & options) {
#elif defined(OPTKING_PACKAGE_QCHEM)
OptReturnType optking(void) {
#endif

  using namespace std;
  using namespace opt;

  open_output_dat();   // assign output.dat file pointer

  MOLECULE *mol1;
  bool newly_generated_coordinates; // Are internal coordinates produced or read-in?

  print_title(); // print header

  try {

#if defined(OPTKING_PACKAGE_PSI)
  set_params(options);
#else
  set_params(); // set optimization parameters
#endif
  if (INTCO_EXCEPT::override_fragment_mode && (Opt_params.fragment_mode == OPT_PARAMS::MULTI))
    Opt_params.fragment_mode = OPT_PARAMS::SINGLE;

  // try to open old internal coordinates
  std::ifstream if_intco(FILENAME_INTCO_DAT, ios_base::in);

  if (if_intco.is_open()) { // old internal coordinates are present

    fprintf(outfile,"\n\tPrevious internal coordinate definitions found.\n"); fflush(outfile);
    newly_generated_coordinates = false;

    mol1 = new MOLECULE(0);    // make an empty molecule

    // read internal coordinate and fragment definitions
    // create, allocate, and add fragment objects
    mol1->read_intcos(if_intco);
#if defined(OPTKING_PACKAGE_PSI)
    psi::Communicator::world->sync();
#endif
    if_intco.close();

    mol1->update_connectivity_by_bonds();

    // read geometry and gradient into the existing fragments
    mol1->read_geom_grad();
    mol1->set_masses();

  }
  else { // automatically generate coordinates

    fprintf(outfile,"\n\tInternal coordinates to be generated automatically.\n"); fflush(outfile);
    newly_generated_coordinates = true;

    // read number of atoms ; make one fragment of that size ;
    mol1 = new MOLECULE( read_natoms() );

    // read geometry and gradient into fragment
    mol1->read_geom_grad();
    mol1->set_masses();

    // use covalent radii to define bonds
    mol1->update_connectivity_by_distances();

    // if fragment_mode == SINGLE, connects all separated groups of atoms by modifying frag.connectivity
    // if fragment_mode == MULTI, splits into fragments and makes interfragment coordinates
    mol1->fragmentize();
mol1->print_connectivity(outfile);

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
    psi::Communicator::world->sync();
    fclose(fp_intco);

    // only generate coordinates and print them out
    if (Opt_params.generate_intcos_only) {
      fprintf(outfile,"\tGenerating intcos and halting.");
      close_output_dat();
#if defined(OPTKING_PACKAGE_PSI)
      psi::psiclean();
#endif
      return OptReturnEndloop;
    }

  }

#if defined(OPTKING_PACKAGE_QCHEM)
  if (Opt_params.efp_fragments_only)
    mol1 = new MOLECULE(0);       // construct an empty molecule
  if (Opt_params.efp_fragments)
    mol1->add_efp_fragments();    // add EFP fragments
#endif

  // print geometry and gradient
  mol1->print_geom_grad(outfile);
  //mol1->print_connectivity(outfile); fflush(outfile);

  // read binary file for previous steps ; history needed to compute EFP values
  p_Opt_data = new OPT_DATA(mol1->g_nintco(), 3*mol1->g_natom());
  if (p_Opt_data->g_iteration() == 1 && Opt_params.opt_type == OPT_PARAMS::IRC) {
    p_irc_data = new IRC_DATA();
    fprintf(stdout,"IRC data object created\n");
  }

#if defined(OPTKING_PACKAGE_QCHEM)
  mol1->update_efp_values(); // EFP values calculated from old opt_data
#endif

  // print internal coordinate definitions and values
  mol1->print_intcos(outfile);

  if (Opt_params.test_B)
    mol1->test_B();

  if (Opt_params.test_derivative_B)
    mol1->test_derivative_B();

  // save geometry and energy
  double * x = mol1->g_geom_array();
  p_Opt_data->save_geom_energy(x, mol1->g_energy());
  if (x!=NULL) free_array(x);

  // print out report on progress
  p_Opt_data->previous_step_report();
  p_Opt_data->reset_consecutive_backsteps();

  // compute forces in internal coordinates from cartesian gradient
  mol1->forces(); // puts forces in p_Opt_data->step[last one]

  if (p_Opt_data->g_iteration() == 1) { // 1st iteration -> put initial Hessian in p_Opt_data
    bool read_H_worked = false;
    if (Opt_params.read_cartesian_H) {
      read_H_worked = mol1->cartesian_H_to_internals(); // read and transform cartesian Hessian
      if (read_H_worked)
        fprintf(outfile,"\tRead in cartesian Hessian and transformed it.\n");
      else 
        fprintf(outfile,"\tUnable to read and transform Hessian.\n");
    }

    if ( (!Opt_params.read_cartesian_H) || (!read_H_worked) )
      mol1->H_guess(); // empirical model guess Hessian
  }
  else { // do Hessian update
    try {
//      if(Opt_params.opt_type != OPT_PARAMS::IRC)
        p_Opt_data->H_update(*mol1);
    } catch (const char * str) {
      fprintf(stderr, "%s\n", str);
      fprintf(outfile, "H_update failed\n");
      fprintf(outfile, "%s\n", str);
      return OptReturnFailure;
    }
  }

//mol1->project_f_and_H();

  // step functions put dq in p_Opt_data->step
  if (Opt_params.opt_type == OPT_PARAMS::IRC)
    mol1->irc_step();
  else {
    if (Opt_params.step_type == OPT_PARAMS::NR)
      mol1->nr_step();
    else if (Opt_params.step_type == OPT_PARAMS::RFO)
      mol1->rfo_step(); 
    else if (Opt_params.step_type == OPT_PARAMS::P_RFO)
      mol1->prfo_step();
    else if (Opt_params.step_type == OPT_PARAMS::SD)
      mol1->sd_step();
  }

  bool converged = p_Opt_data->conv_check(*mol1);

#if defined(OPTKING_PACKAGE_QCHEM)
  rem_write((int) converged, REM_GEOM_OPT_CONVERGED); // tell QChem if converged (return value ignored for now)
  rem_write(p_Opt_data->g_iteration(), REM_GEOM_OPT_CYCLE); // tell QChem current iteration number
#endif
  if ( converged ) {
    if (Opt_params.opt_type == OPT_PARAMS::IRC)
    {
cout << "Converged point!\nSize of opt_data is: " << p_Opt_data->nsteps() << "\n";
//   TODO : could delete old opt_data entries
      //delete all entries but those on reaction path
      //assuming coord has already been incremented; is >=1
//      while(p_Opt_data->nsteps() > 2)
//        p_Opt_data->erase_step(1);
//      p_Opt_data->H_update(*mol1);
//      p_Opt_data->erase_step(0);

      p_irc_data->point_converged(*mol1);
    }
    else
    {
      fprintf(outfile,"\n\t **** Optimization is complete! ****\n");
      p_Opt_data->summary();
      p_Opt_data->write(); // save data to optimization binary file

      fprintf(outfile,"\tFinal energy is %20.13lf\n", p_Opt_data->g_energy());

      if (Opt_params.write_final_step_geometry) {
        fprintf(outfile,"\tFinal (next step) structure:\n");
        mol1->print_geom();  // write geometry -> output file
        fprintf(outfile,"\tSaving final (next step) structure.\n");
      }
      else { // default - get last geometry and write that one
        double *x = p_Opt_data->g_geom_const_pointer(p_Opt_data->nsteps()-1);
        mol1->set_geom_array(x);
        fprintf(outfile,"\tFinal (previous) structure:\n");
        mol1->print_geom();  // write geometry -> output file
        fprintf(outfile,"\tSaving final (previous) structure.\n");
      }

      delete p_Opt_data;
      mol1->write_geom();  // write geometry -> chkpt file (also output for QChem)
      print_end();

      close_output_dat();
      return OptReturnEndloop;
    }
  }

  p_Opt_data->write();
  delete p_Opt_data;

//#if defined(OPTKING_PACKAGE_PSI)
  fprintf(outfile,"\tStructure for next step:\n");
  mol1->print_geom(); // write geometry for next step to output file
//#endif

  mol1->write_geom(); // write geometry -> chkpt (also output for QChem)
  print_end();

  } // end big try
  catch (INTCO_EXCEPT exc) {

    if (exc.try_again() && !exc.already_tried_other_intcos) {

      fprintf(outfile,"\tThe optimizer encountered the following error:\n\t%s\n", exc.g_message());
      fprintf(outfile,"\tWill attempt to restart optimization with redefined internal coordinates.\n");

      opt_intco_dat_remove(); // rm intco definitions
      opt_io_remove(); // rm optimization data

#if defined(OPTKING_PACKAGE_QCHEM)
      rem_write(0, REM_GEOM_OPT_CYCLE); // reset iteration counter
#endif
      // if multi mode has failed, for now try single mode.
      if (Opt_params.fragment_mode == OPT_PARAMS::MULTI)
        exc.override_fragment_mode = true;

      exc.already_tried_other_intcos = true;
      close_output_dat();
      return OptReturnSuccess;
    }
    else {
      fprintf(stderr, "%s\n", exc.g_message());
      fprintf(outfile, "%s\n", exc.g_message());
#if defined (OPTKING_PACKAGE_QCHEM)
      QCrash(exc.g_message());
#elif defined (OPTKING_PACKAGE_PSI)
      abort();
#endif
    }
  }
  catch (BAD_STEP_EXCEPT exc) {

    fprintf(outfile,"\tA bad-step exception has been caught.\n");
    fprintf(outfile,"\t%s", exc.g_message());

    p_Opt_data->decrease_trust_radius();

    mol1->backstep();

    p_Opt_data->write();
    delete p_Opt_data;
    fprintf(outfile,"\tStructure for next step:\n");
    mol1->print_geom(); // write geometry for next step to output file

    mol1->write_geom(); // write geometry -> chkpt (also output for QChem)
    print_end();

    close_output_dat();
    return OptReturnSuccess;
  }
  catch (...) {
#if defined (OPTKING_PACKAGE_QCHEM)
      QCrash(exc.g_message());
#elif defined (OPTKING_PACKAGE_PSI)
      abort();
#endif
  }

  close_output_dat();
  return OptReturnSuccess;
}

// Functions to set and release (possibly open and close) the file pointer for
// the standard text output file
void open_output_dat(void) {
#if defined (OPTKING_PACKAGE_PSI)
  outfile = psi::outfile;
#elif defined (OPTKING_PACKAGE_QCHEM)
  outfile = stdout;
#endif
}

void close_output_dat(void) {
#if defined(OPTKING_PACKAGE_QCHEM)
  fflush(outfile);
#endif
}

void print_title(void) {
  fprintf(outfile, "\n\t\t\t----------------------------------\n");
  fprintf(outfile, "\t\t\t OPTKING: for geometry optimizations  \n");
  fprintf(outfile, "\t\t\t  - R.A. King,  Bethel University   \n");
  fprintf(outfile, "\t\t\t------------------------------------\n");
  fflush(outfile);
}

void print_end(void) {
#if defined (OPTKING_PACKAGE_PSI)
  fprintf(outfile, "\t\t\t--------------------------\n");
  fprintf(outfile, "\t\t\t OPTKING Finished Execution \n");
  fprintf(outfile, "\t\t\t--------------------------\n");
#endif
}

void init_ioff(void)
{
  int i;
  ioff = init_int_array(IOFF_MAX);
  ioff[0] = 0;
  for(i=1; i < IOFF_MAX; i++) ioff[i] = ioff[i-1] + i;
}

}


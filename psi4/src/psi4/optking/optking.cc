/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */

/*!
    \defgroup optking
    \file optking.cc : main optimizer
    \ingroup optking
*/
#include <cstdio>
#include <fstream>
#include <iostream>
#include "globals.h"
#include "molecule.h"
#include "print.h"
#include "io.h"

#if defined(OPTKING_PACKAGE_PSI)
  #include "psi4/libpsi4util/exception.h"
  #include "psi4/libparallel/parallel.h"
#endif

// Define the return types for optking.
#if defined(OPTKING_PACKAGE_PSI)
  typedef psi::PsiReturnType OptReturnType;
  #define OptReturnSuccess (psi::Success)
  #define OptReturnFailure (psi::Failure)
  #define OptReturnEndloop (psi::EndLoop)
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
 OptReturnType opt2man_main(void) { return opt::optking(); } // QCHEM wrapper/alias
#endif

namespace opt {
  void open_output_dat(void); // open/link outfile to text output
  void close_output_dat(void);// close above
  void print_title_out(void); // print header
  void print_end_out(void);   // print footer
  void init_ioff(void);
  int INTCO_EXCEPT::dynamic_level = 0;          // initialized only once
  std::vector<int> INTCO_EXCEPT::linear_angles; // static class members must be defined
  //bool INTCO_EXCEPT::already_tried_other_intcos = false;
  //bool INTCO_EXCEPT::override_fragment_mode = false; // to override MULTI setting by exception algorithm

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

  print_title_out(); // print header

  try { // STEP2 ; do other exceptions
  try { // STEP1 ; check bad_step exception first

#if defined(OPTKING_PACKAGE_PSI)
  set_params(options);
#else
  set_params(); // set optimization parameters
#endif
  //if (INTCO_EXCEPT::override_fragment_mode && (Opt_params.fragment_mode == OPT_PARAMS::MULTI))
  //  Opt_params.fragment_mode = OPT_PARAMS::SINGLE;

  // try to open old internal coordinates
  std::ifstream if_intco(FILENAME_INTCO_DAT, ios_base::in);

  if (if_intco.is_open()) { // old internal coordinates are present

    oprintf_out("\n\tPrevious internal coordinate definitions found.\n");
    newly_generated_coordinates = false;

    mol1 = new MOLECULE(0);    // make an empty molecule

    // read internal coordinate and fragment definitions
    // create, allocate, and add fragment objects
    mol1->read_coords(if_intco);

    if_intco.close();
    //mol1->form_intrafragment_coord_combinations();  // skip if already present

    // If we only read simples, and there are no combinations indicated, then
    // assume 'trivial' (no combos);
    if (mol1->Ncoord() == 0)
      mol1->form_trivial_coord_combinations();

    mol1->update_connectivity_by_bonds();

    // read geometry and gradient into the existing fragments
    mol1->read_geom_grad();
    mol1->set_masses();

  }
  else { // automatically generate coordinates

    oprintf_out("\n\tInternal coordinates to be generated automatically.\n");
    newly_generated_coordinates = true;

    // read number of atoms ; make one fragment of that size ;
    mol1 = new MOLECULE( read_natoms() );

    // read geometry and gradient into fragment
    mol1->read_geom_grad();

    // Quit nicely if there is only one atom present
    if (mol1->g_natom() == 1) {
      oprintf_out("\tThere is only one atom present, so your optimization is complete!\n");
      close_output_dat();
      return OptReturnEndloop;
    }

    mol1->set_masses();

    // use covalent radii to define connectivity
    mol1->update_connectivity_by_distances();
    if ( Opt_params.coordinates != OPT_PARAMS::CARTESIAN )
      mol1->add_intrafragment_hbonds();

    // Critical function here.
    // If fragment_mode == SINGLE, it connects all separated groups of atoms with additional connections.
    // If fragment_mode == MULTI, splits into fragments based on connectivity.
    mol1->fragmentize();
    mol1->print_connectivity(psi_outfile, qc_outfile);

    // General simple internal coordinates
    if (Opt_params.coordinates == OPT_PARAMS::REDUNDANT ||
        Opt_params.coordinates == OPT_PARAMS::DELOCALIZED ||
        Opt_params.coordinates == OPT_PARAMS::NATURAL ||
        Opt_params.coordinates == OPT_PARAMS::BOTH) {

          mol1->add_intrafragment_simples_by_connectivity();

          if (Opt_params.add_auxiliary_bonds) { // adds auxiliary bonds (only) within existing fragments
            if (Opt_params.coordinates == OPT_PARAMS::NATURAL)
              oprintf_out("\tAuxiliary bonds not added, as considered contrary to natural internals.\n");
            else
              mol1->add_intrafragment_auxiliary_bonds();
          }
    }

    if (mol1->g_nfragment() > 1) { // >1 EFP fragment need connected
      mol1->add_interfragment();
      mol1->freeze_interfragment_asymm(); // Try to remove some assymetric problematic ones.
    }

    if ( Opt_params.coordinates == OPT_PARAMS::CARTESIAN )
      mol1->add_cartesians();

    // In the future, which all types of constraints will be support?
    mol1->apply_input_constraints();

    if (Opt_params.coordinates == OPT_PARAMS::DELOCALIZED)
      mol1->form_delocalized_coord_combinations();
    else if (Opt_params.coordinates == OPT_PARAMS::NATURAL)
      mol1->form_natural_coord_combinations();
    else
      mol1->form_trivial_coord_combinations(); // Add trivial combinations.

    if (Opt_params.freeze_intrafragment) {
      mol1->freeze_intrafragments();
      mol1->freeze_intrafragment_coords();
    }

    if ( Opt_params.coordinates == OPT_PARAMS::BOTH )
      mol1->add_cartesians(); // also adds trivial combos

    // print out internal coordinates for future steps
    FILE *qc_intco = NULL;
    std::string psi_intco = FILENAME_INTCO_DAT;
#if defined(OPTKING_PACKAGE_QCHEM)
    qc_intco = fopen(FILENAME_INTCO_DAT, "w");
#endif
    mol1->print_intco_dat(psi_intco, qc_intco);
#if defined(OPTKING_PACKAGE_QCHEM)
    fclose(qc_intco);
#endif

    // only generate coordinates and print them out
    if (Opt_params.intcos_generate_exit) {
      oprintf_out("\tUpon request, generating intcos and halting.\n");
      close_output_dat();
      return OptReturnEndloop;
    }

  }

  if (Opt_params.fb_fragments_only)
    mol1 = new MOLECULE(0);       // construct an empty molecule
  if (Opt_params.fb_fragments)
    mol1->add_fb_fragments();    // add EFP fragments

  // print geometry and gradient
  mol1->print_geom_grad(psi_outfile, qc_outfile);
  //mol1->print_connectivity(psi_outfile, qc_outfile);

  // read binary file for previous steps ; history needed to compute EFP values
  p_Opt_data = new OPT_DATA(mol1->Ncoord(), 3*mol1->g_natom());
  if (p_Opt_data->g_iteration() == 1 && Opt_params.opt_type == OPT_PARAMS::IRC) {
    p_irc_data = new IRC_DATA();
    oprintf_out("IRC data object created\n");
  }

  // If first iteration, start optional xyz trajectory file.
  if (p_Opt_data->g_iteration() == 1)
    if (Opt_params.print_trajectory_xyz_file && Opt_params.opt_type != OPT_PARAMS::IRC)
      mol1->print_xyz(-1);

  mol1->update_fb_values(); // EFP values calculated from old opt_data

  // print internal coordinate definitions and values
  mol1->print_coords(psi_outfile, qc_outfile);

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

  bool read_H_worked = false;

  if (Opt_params.step_type != OPT_PARAMS::SD) {  // ignore all hessian stuff if SD

    if (Opt_params.H_guess_every) { // ignore Hessian already present
        mol1->H_guess(); // empirical model guess Hessian
    }  /* if one wants to read_cartesian_H every time, then user can set this to true and H_update to none */
    else if (Opt_params.read_cartesian_H) {
      double **H_cart = p_Opt_data->read_cartesian_H();       // read Cartesian Hessian
      read_H_worked = mol1->cartesian_H_to_internals(H_cart); // transform to internal coordinates
      free_matrix(H_cart);                                    // free Cartesian Hessian
      if (read_H_worked) {
        oprintf_out("\tRead in cartesian Hessian and transformed it.\n");
        p_Opt_data->reset_steps_since_last_H();
      }
      else {
        oprintf_out("\tUnable to read and transform cartesian Hessian.\n");
        mol1->H_guess(); // empirical model guess Hessian
      }
    }
    else if (p_Opt_data->g_iteration() == 1) {
        mol1->H_guess(); // empirical model guess Hessian
    }
    else { // do Hessian update
      try {
          p_Opt_data->H_update(*mol1);
      } catch (const char * str) {
        fprintf(stderr, "%s\n", str);
        oprintf_out( "H_update failed\n");
        oprintf_out( "%s\n", str);
        return OptReturnFailure;
      }
    }
  } // end !steepest descent

  // Increase number of steps since last Hessian by 1 unless we read in a new one
  if (read_H_worked)
    p_Opt_data->reset_steps_since_last_H();
  else
    p_Opt_data->increment_steps_since_last_H();

  // apply user-defined constraints
  mol1->apply_constraint_forces();
  // project out constraints for fixed intcos and unphysical displacements
  mol1->project_f_and_H();

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
    else if (Opt_params.step_type == OPT_PARAMS::LINESEARCH_STATIC) {
      // compute geometries and then quit
      mol1->linesearch_step();
      delete p_Opt_data;
      print_end_out();
      close_output_dat();
      return OptReturnEndloop;
    }
  }

  bool converged = p_Opt_data->conv_check(*mol1);

#if defined(OPTKING_PACKAGE_QCHEM)
  rem_write((int) converged, REM_GEOM_OPT_CONVERGED); // tell QChem if converged (return value ignored for now)
  rem_write(p_Opt_data->g_iteration(), REM_GEOM_OPT_CYCLE); // tell QChem current iteration number
#endif

  if(Opt_params.opt_type == OPT_PARAMS::IRC && p_Opt_data->g_iteration() == Opt_params.geom_maxiter)
    p_irc_data->progress_report(*mol1);

  if ( converged ) {
    if (Opt_params.opt_type == OPT_PARAMS::IRC) {

      if(!p_irc_data->go) {
        oprintf_out("\n\t **** Optimization is complete! ****\n");
        p_Opt_data->write(); // save data to optimization binary file
        oprintf_out("\tFinal energy is %20.13lf\n", p_Opt_data->g_energy());

        if (Opt_params.write_final_step_geometry) {
          oprintf_out("\tFinal (next step) structure:\n");
          mol1->print_geom_out();  // write geometry -> output file
          oprintf_out("\tSaving final (next step) structure.\n");
        }
        else { // default - get last geometry and write that one
          double *x = p_irc_data->g_x();
          mol1->set_geom_array(x);
          oprintf_out("\tFinal (previous) structure:\n");
          mol1->print_geom_out();  // write geometry -> output file
          oprintf_out("\tSaving final (previous) structure.\n");
        }
        p_irc_data->progress_report(*mol1);

        p_Opt_data->reset_trust_radius();
        delete p_Opt_data;
        INTCO_EXCEPT::dynamic_level = options.get_int("DYNAMIC_LEVEL"); // reset for future optimizations
        mol1->write_geom();  // write geometry -> chkpt file (also output for QChem)
        print_end_out();

        close_output_dat();
        return OptReturnEndloop;
      }

      INTCO_EXCEPT::dynamic_level = options.get_int("DYNAMIC_LEVEL"); // reset for future optimizations
      p_irc_data->point_converged(*mol1);
    }
    else
    {
      oprintf_out("\n  **** Optimization is complete! (in %d steps) ****\n", p_Opt_data->nsteps());
      p_Opt_data->summary();
      p_Opt_data->write(); // save data to optimization binary file

      oprintf_out("\tFinal energy is %20.13lf\n", p_Opt_data->g_energy());

      if (Opt_params.write_final_step_geometry) {
        oprintf_out("\tFinal (next step) structure:\n");
        mol1->print_geom_out();  // write geometry -> output file
        if (Opt_params.print_trajectory_xyz_file)
          mol1->print_xyz();
        oprintf_out("\tSaving final (next step) structure.\n");
      }
      else { // default - get last geometry and write that one
        double *x = p_Opt_data->g_geom_const_pointer(p_Opt_data->nsteps()-1);
        mol1->set_geom_array(x);
        oprintf_out("\tFinal (previous) structure:\n");
        mol1->print_geom_out();  // write geometry -> output file
        if (Opt_params.print_trajectory_xyz_file) mol1->print_xyz();
        oprintf_out("\tSaving final (previous) structure.\n");
      }

      p_Opt_data->reset_trust_radius();
      delete p_Opt_data;
      INTCO_EXCEPT::dynamic_level = options.get_int("DYNAMIC_LEVEL"); // reset for future optimizations
      opt_intco_dat_remove(); // rm intco definitions
      opt_io_remove();        // rm optimization data
      mol1->write_geom();  // write geometry -> chkpt file (also output for QChem)
      print_end_out();

      close_output_dat();
      return OptReturnEndloop;
    }
  }

  if(p_Opt_data->g_iteration() == Opt_params.geom_maxiter) {
    oprintf_out("\n  **** Optimization has failed! (in %d steps) ****\n", p_Opt_data->nsteps());
      p_Opt_data->summary();
      p_Opt_data->write(); // save data to optimization binary file

/*  We could print out this info if it didn't confuse users.
      oprintf_out("\tFinal energy is %20.13lf\n", p_Opt_data->g_energy());
      oprintf_out("\tFinal (next step) structure:\n");
      mol1->print_geom_out();  // write geometry -> output file
      if (Opt_params.print_trajectory_xyz_file) mol1->print_xyz();
      oprintf_out("\tSaving final (next step) structure.\n");
*/

      p_Opt_data->reset_trust_radius();
      delete p_Opt_data;
      INTCO_EXCEPT::dynamic_level = options.get_int("DYNAMIC_LEVEL"); // reset for future optimizations
      mol1->write_geom();  // write geometry -> chkpt file (also output for QChem)
      print_end_out();
      close_output_dat();
      return OptReturnSuccess;
  }

  p_Opt_data->write();
  delete p_Opt_data;

  oprintf_out("\tStructure for next step:\n");
  mol1->print_geom_out(); // write geometry for next step to output file
  if (Opt_params.print_trajectory_xyz_file && Opt_params.opt_type != OPT_PARAMS::IRC)
    mol1->print_xyz();

  mol1->write_geom(); // write geometry -> chkpt (also output for QChem)
  print_end_out();

  } // end try STEP1
  catch (BAD_STEP_EXCEPT exc) {
    oprintf_out("\tThe BAD_STEP_EXCEPTion handler:\n\t%s\n", exc.g_message());
    oprintf_out("\tDynamic level is %d.\n", Opt_params.dynamic);
    oprintf_out("\tConsecutive backsteps is %d.\n", p_Opt_data->g_consecutive_backsteps());

    if ( (!Opt_params.dynamic || (p_Opt_data->g_consecutive_backsteps() < Opt_params.consecutive_backsteps_allowed))
        && (p_Opt_data->g_iteration() > 1) ) {
      // Do backward step.  backstep() function increments counter in opt_data.
      p_Opt_data->decrease_trust_radius();
      mol1->backstep();
      p_Opt_data->write();
      delete p_Opt_data;

      oprintf_out("\tStructure for next step:\n");
      mol1->print_geom_out(); // write geometry for next step to output file
      if (Opt_params.print_trajectory_xyz_file && Opt_params.opt_type != OPT_PARAMS::IRC)
        mol1->print_xyz();

      mol1->write_geom(); // write geometry -> chkpt (also output for QChem)
      print_end_out();
      close_output_dat();
    }
    else if (Opt_params.dynamic) {
      throw(INTCO_EXCEPT("Too many bad steps."));
    }
    else {
      oprintf_out("This exception should not be called when neither dynamic nor backsteps are allowed.\n");
      throw("This exception should not be called when neither dynamic nor backsteps are allowed.");
    }
    return OptReturnSuccess;
  } // end BAD_STEP_EXCEPT
  } // end try STEP2
  catch (BROKEN_SYMMETRY_EXCEPT exc) {
    p_Opt_data->reset_trust_radius();
    oprintf_out("\n  **** Optimization has failed! (in %d steps) ****\n", p_Opt_data->nsteps());
    delete p_Opt_data;
    print_end_out();
    close_output_dat();
    return OptReturnFailure;
  }
  catch (INTCO_EXCEPT exc) {
    oprintf_out("\tThe INTCO_EXCEPTion handler:\n\t%s\n", exc.g_message());
    fprintf(stderr, "%s\n", exc.g_message());
    oprintf_out("\tDynamic level is %d.\n", Opt_params.dynamic);
    oprintf_out("\texc.g_really_quit() is %d.\n", exc.g_really_quit());

    if (exc.g_really_quit()) { // quit no matter what is indicated by exception
      oprintf_out("\n  **** Optimization has failed! (in %d steps) ****\n", p_Opt_data->nsteps());
      delete p_Opt_data;
      print_end_out();
      close_output_dat();
      return OptReturnFailure;
    }

    // There is no indication that new linear coordinates need to be added AND we are
    // either not doing a dynamic optimization, or we have tried all available levels.
    if ( (Opt_params.dynamic == 0 || Opt_params.dynamic == 7) && exc.linear_angles.empty()) {
      oprintf_out("\n  **** Optimization has failed! (in %d steps) ****\n", p_Opt_data->nsteps());
      delete p_Opt_data;
      print_end_out();
      close_output_dat();
      return OptReturnFailure;
/*
#if defined (OPTKING_PACKAGE_QCHEM)
      QCrash(exc.g_message());
#elif defined (OPTKING_PACKAGE_PSI)
      abort();
#endif
*/
    }
    else {
      if (exc.linear_angles.empty()) {
        oprintf_out("\tRaising dynamic level to %d. ~\n", Opt_params.dynamic+1);
        exc.increment_dynamic_level();
      }
      else
        oprintf_out("\tRestarting with new internal coordinates at same level.");

      int use_geom = 0;
      if (p_Opt_data->nsteps()>2)
        use_geom = p_Opt_data->nsteps()-2;
      else if (p_Opt_data->nsteps()>1)
        use_geom = p_Opt_data->nsteps()-1;
      else
        use_geom = 0; // probably need to check this

      double *x = p_Opt_data->g_geom_const_pointer(use_geom);
      mol1->set_geom_array(x);
      mol1->write_geom();

      delete p_Opt_data;
      opt_intco_dat_remove(); // rm intco definitions
      opt_io_remove();        // rm optimization data

      //if (Opt_params.fragment_mode == OPT_PARAMS::MULTI)
      //  exc.override_fragment_mode = true;

#if defined(OPTKING_PACKAGE_QCHEM)
      rem_write(0, REM_GEOM_OPT_CYCLE); // reset iteration counter
#endif
      close_output_dat();
      return OptReturnSuccess;
    }
  }
#if defined (OPTKING_PACKAGE_PSI)
  catch (psi::PsiException e){
      oprintf_out("\t%s", e.what());
  }
#endif
  catch (...) {
#if defined (OPTKING_PACKAGE_QCHEM)
      QCrash("Exception thrown in optking() leading to abort.");
#elif defined (OPTKING_PACKAGE_PSI)
      abort();
#endif
  }

  close_output_dat();
  return OptReturnSuccess;
}

// Standard text output file string (psi) or file pointer (qchem)
// Interpreted by functions in print.cc
void open_output_dat(void) {
#if defined (OPTKING_PACKAGE_PSI)
  psi_outfile = "outfile";
#elif defined (OPTKING_PACKAGE_QCHEM)
  qc_outfile = stdout;
#endif
}

void close_output_dat(void) {
#if defined(OPTKING_PACKAGE_QCHEM)

#endif
}

void print_title_out(void) {
  oprintf_out( "\n\t\t\t-----------------------------------------\n");
  oprintf_out(   "\t\t\t OPTKING 2.0: for geometry optimizations \n");
  oprintf_out(   "\t\t\t  - R.A. King,  Bethel University        \n");
  oprintf_out(   "\t\t\t-----------------------------------------\n");

}

void print_end_out(void) {
#if defined (OPTKING_PACKAGE_PSI)
  oprintf_out( "\t\t\t--------------------------\n");
  oprintf_out( "\t\t\t OPTKING Finished Execution \n");
  oprintf_out( "\t\t\t--------------------------\n");
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

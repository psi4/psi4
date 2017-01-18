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

/*! \file set_params.cc
    \ingroup optking
    \brief set optimization parameters
*/

#include "print.h"
#define EXTERN
#include "globals.h"

#if defined(OPTKING_PACKAGE_PSI)
 #include "psi4/psi4-dec.h"
#elif defined(OPTKING_PACKAGE_QCHEM)
 #include <qchem.h>
#endif

#include <sstream>

namespace opt {

void print_params_out(void);

#if defined(OPTKING_PACKAGE_PSI)
void set_params(psi::Options & options)
#else
void set_params(void)
#endif
{

  std::string s;

#if defined(OPTKING_PACKAGE_PSI)

// optimization type
    s = options.get_str("OPT_TYPE");
    if(s == "MIN")  Opt_params.opt_type = OPT_PARAMS::MIN;
    else if (s == "TS")  Opt_params.opt_type = OPT_PARAMS::TS;
    else if (s == "IRC")  Opt_params.opt_type = OPT_PARAMS::IRC;

    if (options["STEP_TYPE"].has_changed()) {
      s = options.get_str("STEP_TYPE");
      if (s == "RFO") {
        if (Opt_params.opt_type == OPT_PARAMS::MIN)
          Opt_params.step_type = OPT_PARAMS::RFO;
        else if (Opt_params.opt_type == OPT_PARAMS::TS)
          Opt_params.step_type = OPT_PARAMS::P_RFO;
      }
      else if (s == "NR") Opt_params.step_type = OPT_PARAMS::NR;
      else if (s == "SD") Opt_params.step_type = OPT_PARAMS::SD;
      else if (s == "LINESEARCH_STATIC") Opt_params.step_type = OPT_PARAMS::LINESEARCH_STATIC;
   }
   else { // Set defaults for step type.
     if (Opt_params.opt_type == OPT_PARAMS::MIN)
       Opt_params.step_type = OPT_PARAMS::RFO;
     else if (Opt_params.opt_type == OPT_PARAMS::TS)
       Opt_params.step_type = OPT_PARAMS::P_RFO;
     // else if (Opt_params.opt_type == OPT_PARAMS::IRC) options?
   }

   s = options.get_str("OPT_COORDINATES");
   if (s == "INTERNAL" || s == "REDUNDANT")
     Opt_params.coordinates = OPT_PARAMS::REDUNDANT;
   else if (s == "DELOCALIZED")
     Opt_params.coordinates = OPT_PARAMS::DELOCALIZED;
   else if (s == "NATURAL")
     Opt_params.coordinates = OPT_PARAMS::NATURAL;
   else if (s == "CARTESIAN")
     Opt_params.coordinates = OPT_PARAMS::CARTESIAN;
   else if (s == "BOTH")
     Opt_params.coordinates = OPT_PARAMS::BOTH;

// Maximum step size in bohr or radian along an internal coordinate {double}
   // For fixed coordinate optimizations, also see below.
//  Opt_params.intrafragment_step_limit = 0.4;
    Opt_params.intrafragment_step_limit = options.get_double("INTRAFRAG_STEP_LIMIT");
    Opt_params.intrafragment_step_limit_min = options.get_double("INTRAFRAG_STEP_LIMIT_MIN");
    Opt_params.intrafragment_step_limit_max = options.get_double("INTRAFRAG_STEP_LIMIT_MAX");

    Opt_params.interfragment_step_limit = options.get_double("INTERFRAG_STEP_LIMIT");

// For now, the initial Hessian guess for cartesians with BOTH is stupid, so don't scale
// step size down too much.
    if (!options["INTRAFRAG_STEP_LIMIT_MIN"].has_changed() && Opt_params.coordinates == OPT_PARAMS::BOTH)
      Opt_params.intrafragment_step_limit_min = Opt_params.intrafragment_step_limit / 2.0;

// Steepest descent has no good hessian, so don't restrict step much
    if (!options["INTRAFRAG_STEP_LIMIT_MIN"].has_changed() && Opt_params.step_type == OPT_PARAMS::SD)
      Opt_params.intrafragment_step_limit_min = Opt_params.intrafragment_step_limit;

// Reduce step size to ensure convergence of back-transformation of internal coordinate
// step to cartesians.
    Opt_params.ensure_bt_convergence = options.get_bool("ENSURE_BT_CONVERGENCE");

// do stupid, linear scaling of internal coordinates to step limit (not RS-RFO);
    Opt_params.simple_step_scaling = options.get_bool("SIMPLE_STEP_SCALING");

// Whether to 'follow' the initial RFO vector after the first step {true, false}
    Opt_params.rfo_follow_root = options.get_bool("RFO_FOLLOW_ROOT");
// Which RFO root to follow; internally 0 (externally 1) indicates minimum; {integer}
//  Opt_params.rfo_root = 0;
    Opt_params.rfo_root = options.get_int("RFO_ROOT");

// When determining connectivity, a bond is assigned if interatomic distance
// is less than (this number) * sum of covalent radii {double}
//  Opt_params.scale_connectivity = 1.3;
    Opt_params.scale_connectivity = options.get_double("COVALENT_CONNECT");

// When determining connectivity BETWEEN FRAGMENTS when fragment mode is set to
// SIMPLE, distance coordinates are created if atoms on different fragments
// are at a distance less than (this number) * sum of covalent radii {double}
// The criterion is gradually increased until all fragments are connected.
//  Opt_params.interfragment_scale_connectivity = 1.8;
    Opt_params.interfragment_scale_connectivity = options.get_double("INTERFRAGMENT_CONNECT");

// Whether to treat multiple molecule fragments as a single bonded molecule;
// or via interfragment coordinates ; a primary difference is that in MULTI mode,
// the interfragment coordinates are not redundant. {SINGLE, MULTI}
    s = options.get_str("FRAG_MODE");
    if (s == "SINGLE")     Opt_params.fragment_mode = OPT_PARAMS::SINGLE;
    else if (s == "MULTI") Opt_params.fragment_mode = OPT_PARAMS::MULTI;

// whether to use fixed linear combinations of atoms as reference points for
//interfragment coordinates or whether to use principal axes {FIXED, PRINCIPAL_AXES}
    s = options.get_str("INTERFRAG_MODE");
    if (s == "FIXED")               Opt_params.interfragment_mode = OPT_PARAMS::FIXED;
    else if (s == "PRINCIPAL_AXES") Opt_params.interfragment_mode = OPT_PARAMS::PRINCIPAL_AXES;

//  Atoms to define reference points on fragments
    // If already specified, (probably because iter == 0); don't add them.
    if (options["FRAG_REF_ATOMS"].has_changed() && !Opt_params.frag_ref_atoms.size()) {
        int nfrag = options["FRAG_REF_ATOMS"].size();
        for (int F=0; F<nfrag; ++F) {
            std::vector<std::vector<int> > frag;
            Opt_params.frag_ref_atoms.push_back(frag);
            int nref = options["FRAG_REF_ATOMS"][F].size();
            if (nref > 3)
                throw("Fragment can have only 3 reference atoms.\n");
            for (int R=0; R<nref; ++R) {
                std::vector<int> ref;
                Opt_params.frag_ref_atoms[F].push_back(ref);
                for (int A=0; A<(int) options["FRAG_REF_ATOMS"][F][R].size(); ++A) {
                    int atom = options["FRAG_REF_ATOMS"][F][R][A].to_integer();
                    Opt_params.frag_ref_atoms[F][R].push_back(atom); // will need decremented on use!
                }
            }
        }
        Opt_params.interfragment_mode = OPT_PARAMS::FIXED;
    }

// Whether to only generate the internal coordinates and then stop {true, false}
    Opt_params.intcos_generate_exit = options.get_bool("INTCOS_GENERATE_EXIT");

// What model Hessian to use to guess intrafragment force constants {SCHLEGEL, FISCHER, SIMPLE, LINDH}
    s = options.get_str("INTRAFRAG_HESS");
    if (s == "FISCHER")       Opt_params.intrafragment_H = OPT_PARAMS::FISCHER;
    else if (s == "SCHLEGEL") Opt_params.intrafragment_H = OPT_PARAMS::SCHLEGEL;
    else if (s == "SIMPLE") Opt_params.intrafragment_H = OPT_PARAMS::SIMPLE;
    else if (s == "LINDH") Opt_params.intrafragment_H = OPT_PARAMS::LINDH;
    else if (s == "LINDH_SIMPLE") Opt_params.intrafragment_H = OPT_PARAMS::LINDH_SIMPLE;

// Re-estimate the hessian every step.  Usually default is false.
    Opt_params.H_guess_every = options.get_bool("H_GUESS_EVERY");

// Original Lindh specification was to redo at every step.
    if (Opt_params.intrafragment_H == OPT_PARAMS::LINDH) {
      if (!options["H_GUESS_EVERY"].has_changed())
        Opt_params.H_guess_every = true;
    }

//  The default for cartesian coordinates will be to use Lindh force field
//  for initial guess, and then switch to BFGS.  This should be a flexible approach for
//  difficult problems.
    if (Opt_params.coordinates == OPT_PARAMS::CARTESIAN) {
      if (!options["INTRAFRAG_HESS"].has_changed())
        Opt_params.intrafragment_H = OPT_PARAMS::LINDH;
      if (!options["H_GUESS_EVERY"].has_changed())
        Opt_params.H_guess_every = false;
    }

// Whether to use the default of FISCHER_LIKE force constants for the initial guess {DEFAULT, FISCHER_LIKE}
    s = options.get_str("INTERFRAG_HESS");
    if (s == "DEFAULT")           Opt_params.interfragment_H = OPT_PARAMS::DEFAULT;
    else if (s == "FISCHER_LIKE") Opt_params.interfragment_H = OPT_PARAMS::FISCHER_LIKE;

// Whether to freeze all fragments rigid
    Opt_params.freeze_intrafragment = options.get_bool("FREEZE_INTRAFRAG");

// Whether to freeze all interfragment modes
    Opt_params.freeze_interfragment = options.get_bool("FREEZE_INTERFRAG");

// Add auxiliary bonds for non-bonded (but nearby) atoms.
    Opt_params.add_auxiliary_bonds = options.get_bool("ADD_AUXILIARY_BONDS");

// Covalent distance times this factor is used to choose extra stretch coordinates
    Opt_params.auxiliary_bond_factor = options.get_double("AUXILIARY_BOND_FACTOR");

// By default, optking prints and saves the last (previous) geometry at the end of an
// optimization, i.e., the one at which a gradient was computed.  If this keyword is
// set to true, then the structure obtained from the last projected step is printed out and saved instead.
//  Opt_params.write_final_step_geometry = false;
    Opt_params.write_final_step_geometry = options.get_bool("FINAL_GEOM_WRITE");

// Choose from supported Hessian updates {NONE, BFGS, MS, POWELL, BOFILL}
    s = options.get_str("HESS_UPDATE");
    if (s == "NONE")        Opt_params.H_update = OPT_PARAMS::NONE;
    else if (s == "BFGS")   Opt_params.H_update = OPT_PARAMS::BFGS;
    else if (s == "MS")     Opt_params.H_update = OPT_PARAMS::MS;
    else if (s == "POWELL") Opt_params.H_update = OPT_PARAMS::POWELL;
    else if (s == "BOFILL") Opt_params.H_update = OPT_PARAMS::BOFILL;

// Set Bofill as default for TS optimizations
    if ( Opt_params.opt_type == OPT_PARAMS::TS || Opt_params.opt_type == OPT_PARAMS::IRC)
      if (!options["HESS_UPDATE"].has_changed())
        Opt_params.H_update = OPT_PARAMS::BOFILL;

//  How many previous steps' data to use in Hessian update; 0=use them all ; {integer}
//  Opt_params.H_update_use_last = 6;
    Opt_params.H_update_use_last = options.get_int("HESS_UPDATE_USE_LAST");

// Whether to limit the magnitutde of changes caused by the Hessian update {true, false}
//  Opt_params.H_update_limit = true;
    Opt_params.H_update_limit = options.get_bool("HESS_UPDATE_LIMIT");

// If the above is true, changes to the Hessian from the update are limited to the larger of
// (H_update_limit_scale)*(the previous value) and H_update_limit_max (in au).
//  Opt_params.H_update_limit_scale = 0.50;
//  Opt_params.H_update_limit_max  = 1.0;
    Opt_params.H_update_limit_scale = options.get_double("HESS_UPDATE_LIMIT_SCALE");
    Opt_params.H_update_limit_max = options.get_double("HESS_UPDATE_LIMIT_MAX");

// Whether to use 1/R(AB) for stretching coordinate between fragments (or just R(AB))
//  Opt_params.interfragment_distance_inverse = false;
    Opt_params.interfragment_distance_inverse = options.get_bool("INTERFRAG_DIST_INV");

// For now, this is a general maximum distance for the definition of H-bonds
//  Opt_params.maximum_H_bond_distance = 4.3;
    Opt_params.maximum_H_bond_distance = options.get_double("H_BOND_CONNECT");

// User-specified direction for IRC mapping
    s = options.get_str("IRC_DIRECTION");
    if (s == "FORWARD")        Opt_params.IRC_direction = OPT_PARAMS::FORWARD;
    else if (s == "BACKWARD")  Opt_params.IRC_direction = OPT_PARAMS::BACKWARD;

// Decide how to quit IRC mapping
// if 'ASK' ==> user will be prompted to decide whether to continue when
// the gradient dotted with the previous gradient goes negative
   s = options.get_str("IRC_STOP");
   if (s == "STOP")        Opt_params.IRC_stop = OPT_PARAMS::STOP;
   else if (s == "ASK")  Opt_params.IRC_stop = OPT_PARAMS::ASK;
   else if (s == "GO")    Opt_params.IRC_stop = OPT_PARAMS::GO;

// General optimization criteria
    Opt_params.i_max_force = false;
    Opt_params.i_rms_force = false;
    Opt_params.i_max_DE = false;
    Opt_params.i_max_disp = false;
    Opt_params.i_rms_disp = false;
    Opt_params.i_untampered = false;

    Opt_params.general_conv = options.get_str("G_CONVERGENCE");
    if (Opt_params.general_conv == "QCHEM") {
      Opt_params.i_untampered = true;
      Opt_params.conv_max_force = 3.0e-4;  Opt_params.i_max_force = true;
      Opt_params.conv_max_DE    = 1.0e-6;  Opt_params.i_max_DE = true;
      Opt_params.conv_max_disp  = 1.2e-3;  Opt_params.i_max_disp = true;
    }
    else if (Opt_params.general_conv == "MOLPRO") {
      Opt_params.i_untampered = true;
      Opt_params.conv_max_force = 3.0e-4;  Opt_params.i_max_force = true;
      Opt_params.conv_max_DE    = 1.0e-6;  Opt_params.i_max_DE = true;
      Opt_params.conv_max_disp  = 3.0e-4;  Opt_params.i_max_disp = true;
    }
    else if (Opt_params.general_conv == "GAU") {
      Opt_params.i_untampered = true;
      Opt_params.conv_max_force = 4.5e-4;  Opt_params.i_max_force = true;
      Opt_params.conv_rms_force = 3.0e-4;  Opt_params.i_rms_force = true;
      Opt_params.conv_max_disp  = 1.8e-3;  Opt_params.i_max_disp = true;
      Opt_params.conv_rms_disp  = 1.2e-3;  Opt_params.i_rms_disp = true;
    }
    else if (Opt_params.general_conv == "GAU_TIGHT") {
      Opt_params.i_untampered = true;
      Opt_params.conv_max_force = 1.5e-5;  Opt_params.i_max_force = true;
      Opt_params.conv_rms_force = 1.0e-5;  Opt_params.i_rms_force = true;
      Opt_params.conv_max_disp  = 6.0e-5;  Opt_params.i_max_disp = true;
      Opt_params.conv_rms_disp  = 4.0e-5;  Opt_params.i_rms_disp = true;
    }
    else if (Opt_params.general_conv == "GAU_VERYTIGHT") {
      Opt_params.i_untampered = true;
      Opt_params.conv_max_force = 2.0e-6;  Opt_params.i_max_force = true;
      Opt_params.conv_rms_force = 1.0e-6;  Opt_params.i_rms_force = true;
      Opt_params.conv_max_disp  = 6.0e-6;  Opt_params.i_max_disp = true;
      Opt_params.conv_rms_disp  = 4.0e-6;  Opt_params.i_rms_disp = true;
    }
    else if (Opt_params.general_conv == "GAU_LOOSE") {
      Opt_params.i_untampered = true;
      Opt_params.conv_max_force = 2.5e-3;  Opt_params.i_max_force = true;
      Opt_params.conv_rms_force = 1.7e-3;  Opt_params.i_rms_force = true;
      Opt_params.conv_max_disp  = 1.0e-2;  Opt_params.i_max_disp = true;
      Opt_params.conv_rms_disp  = 6.7e-3;  Opt_params.i_rms_disp = true;
    }
    else if (Opt_params.general_conv == "TURBOMOLE") {
      Opt_params.i_untampered = true;
      Opt_params.conv_max_force = 1.0e-3;  Opt_params.i_max_force = true;
      Opt_params.conv_rms_force = 5.0e-4;  Opt_params.i_rms_force = true;
      Opt_params.conv_max_DE    = 1.0e-6;  Opt_params.i_max_DE = true;
      Opt_params.conv_max_disp  = 1.0e-3;  Opt_params.i_max_disp = true;
      Opt_params.conv_rms_disp  = 5.0e-4;  Opt_params.i_rms_disp = true;
    }
    else if (Opt_params.general_conv == "CFOUR") {
      Opt_params.i_untampered = true;
      Opt_params.conv_rms_force = 1.0e-4;  Opt_params.i_rms_force = true;
    }
    else if (Opt_params.general_conv == "NWCHEM_LOOSE") {
      Opt_params.i_untampered = true;
      Opt_params.conv_max_force = 4.5e-3;  Opt_params.i_max_force = true;
      Opt_params.conv_rms_force = 3.0e-3;  Opt_params.i_rms_force = true;
      Opt_params.conv_max_disp  = 5.4e-3;  Opt_params.i_max_disp = true;
      Opt_params.conv_rms_disp  = 3.6e-3;  Opt_params.i_rms_disp = true;
    }

// Specific optimization criteria
    if (options["MAX_FORCE_G_CONVERGENCE"].has_changed()) {
      Opt_params.i_untampered = false;
      Opt_params.i_max_force = true;
      Opt_params.conv_max_force = fabs(options.get_double("MAX_FORCE_G_CONVERGENCE"));
    }
    if (options["RMS_FORCE_G_CONVERGENCE"].has_changed()) {
      Opt_params.i_untampered = false;
      Opt_params.i_rms_force = true;
      Opt_params.conv_rms_force = fabs(options.get_double("RMS_FORCE_G_CONVERGENCE"));
    }
    if (options["MAX_ENERGY_G_CONVERGENCE"].has_changed()) {
      Opt_params.i_untampered = false;
      Opt_params.i_max_DE = true;
      Opt_params.conv_max_DE = fabs(options.get_double("MAX_ENERGY_G_CONVERGENCE"));
    }
    if (options["MAX_DISP_G_CONVERGENCE"].has_changed()) {
      Opt_params.i_untampered = false;
      Opt_params.i_max_disp = true;
      Opt_params.conv_max_disp = fabs(options.get_double("MAX_DISP_G_CONVERGENCE"));
    }
    if (options["RMS_DISP_G_CONVERGENCE"].has_changed()) {
      Opt_params.i_untampered = false;
      Opt_params.i_rms_disp = true;
      Opt_params.conv_rms_disp = fabs(options.get_double("RMS_DISP_G_CONVERGENCE"));
    }

    // even if a specific threshold were given, allow for Molpro/Qchem/G03 flex criteria
    if (options.get_bool("FLEXIBLE_G_CONVERGENCE"))
      Opt_params.i_untampered = true;

// Whether to test B matrix and derivative B matrix numerically
//  Opt_params.test_B = false;
//  Opt_params.test_derivative_B = false;
    Opt_params.test_B = options.get_bool("TEST_B");
    Opt_params.test_derivative_B = options.get_bool("TEST_DERIVATIVE_B");

//  Opt_params.print_lvl = 1;
    Opt_params.print_lvl = options.get_int("PRINT");

    Opt_params.print_params = options.get_bool("PRINT_OPT_PARAMS");

// make trajectory file printing the default for IRC.
    if ((Opt_params.opt_type == OPT_PARAMS::IRC) &&
        (options["PRINT_TRAJECTORY_XYZ_FILE"].has_changed() == false))
      Opt_params.print_trajectory_xyz_file = true;
    else
      Opt_params.print_trajectory_xyz_file = options.get_bool("PRINT_TRAJECTORY_XYZ_FILE");

// default


// Read cartesian Hessian.  Make reading the default for IRC.
    if ((Opt_params.opt_type == OPT_PARAMS::IRC) &&
        (options["CART_HESS_READ"].has_changed() == 0))
      Opt_params.read_cartesian_H = true;
    else
      Opt_params.read_cartesian_H = options.get_bool("CART_HESS_READ");

// only treating "dummy fragments"
    // These are not found in psi4/read_options.cc
    // Not sure if we need these.
  Opt_params.fb_fragments = false;
  Opt_params.fb_fragments_only = false;

//IRC stepsize
  Opt_params.IRC_step_size = options.get_double("IRC_STEP_SIZE");

  // keep internal coordinate definitions file after optimization
  Opt_params.keep_intcos = options.get_bool("KEEP_INTCOS");

  // if we are running only to generate them, then we assume we'll keep them
  if (Opt_params.intcos_generate_exit && !options["KEEP_INTCOS"].has_changed())
    Opt_params.keep_intcos = true;

  // for coordinates with user-specified equilibrium values - this is the force constant
  Opt_params.fixed_coord_force_constant = options.get_double("FIXED_COORD_FORCE_CONSTANT");

  // Currently, a static line search merely displaces along the gradient in internal
  // coordinates generating LINESEARCH_STATIC_N geometries.  The other two keywords
  // control the min and the max of the largest internal coordinate displacement.
  Opt_params.linesearch_static_N   = options.get_int("LINESEARCH_STATIC_N");
  Opt_params.linesearch_static_min = options.get_double("LINESEARCH_STATIC_MIN");
  Opt_params.linesearch_static_max = options.get_double("LINESEARCH_STATIC_MAX");

// consecutive number of backsteps allowed before giving up
  Opt_params.consecutive_backsteps_allowed = options.get_int("CONSECUTIVE_BACKSTEPS");

  // if steepest-descent, then make much larger default
  if (Opt_params.step_type == OPT_PARAMS::SD && !(options["CONSECUTIVE_BACKSTEPS"].has_changed()))
    Opt_params.consecutive_backsteps_allowed = 10;

  Opt_params.geom_maxiter = options.get_int("GEOM_MAXITER");

  // For RFO step, eigenvectors of augmented Hessian are divided by the last
  // element unless it is smaller than this value {double}.  Can be used to eliminate
  // asymmetric steps not otherwise detected (e.g. in degenerate point groups).
  // For multi-fragment modes, we presume that smaller Delta-E's are possible, and
  // this threshold should be made larger.
  Opt_params.rfo_normalization_max = options.get_double("RFO_NORMALIZATION_MAX");
  if (Opt_params.fragment_mode == OPT_PARAMS::MULTI &&
      !(options["RFO_NORMALIZATION_MAX"].has_changed()))
    Opt_params.rfo_normalization_max = 1.0e5;

// Hessian update is avoided if the denominators (Dq*Dq) or (Dq*Dg) are smaller than this
  Opt_params.H_update_den_tol = options.get_double("H_UPDATE_DEN_TOL");

// Symmetry tolerance for whether to follow non-symmetric modes.
  Opt_params.symm_tol = options.get_double("SYMM_TOL");

  // Absolute maximum for value of alpha in RS-RFO
  Opt_params.rsrfo_alpha_max = options.get_double("RSRFO_ALPHA_MAX");


#elif defined(OPTKING_PACKAGE_QCHEM)

  int i;

// MIN = 0 ; TS = 1 ; IRC = 2      (default 0)
  i = rem_read(REM_GEOM_OPT2_OPT_TYPE);
  if (i == 0)      Opt_params.opt_type = OPT_PARAMS::MIN;
  else if (i == 1) Opt_params.opt_type = OPT_PARAMS::TS;
  else if (i == 2) Opt_params.opt_type = OPT_PARAMS::IRC;

// RFO = 0 ; NR = 1 ; P_RFO = 2      (default 0)
// defaults should be RFO for MIN; P_RFO for TS
  if (Opt_params.opt_type == OPT_PARAMS::MIN)
    Opt_params.step_type = OPT_PARAMS::RFO;
  else if (Opt_params.opt_type == OPT_PARAMS::TS)
    Opt_params.step_type = OPT_PARAMS::P_RFO;

  // How to check if user specified this?
  i = rem_read(REM_GEOM_OPT2_STEP_TYPE);
  if (i == 0)      Opt_params.step_type = OPT_PARAMS::RFO;
  else if (i == 1) Opt_params.step_type = OPT_PARAMS::NR;
  else if (i == 2) Opt_params.step_type = OPT_PARAMS::P_RFO;

// coordinates for optimization; default is 0
  i = rem_read(REM_GEOM_OPT2_OPT_COORDINATES);
  if (i == 0)
    Opt_params.coordinates = OPT_PARAMS::INTERNAL;
  else if (i == 1)
    Opt_params.coordinates = OPT_PARAMS::CARTESIAN;
  else if (i == 2)
    Opt_params.coordinates = OPT_PARAMS::BOTH;

  // Maximum change in an internal coordinate is au; limits on steps are rem / 1000
  i = rem_read(REM_GEOM_OPT2_INTRAFRAG_STEP_LIMIT);     // default is  400 -> 0.4
  Opt_params.intrafragment_step_limit =     i / 1000.0;
  i = rem_read(REM_GEOM_OPT2_INTRAFRAG_STEP_LIMIT_MIN); // default is    1 -> 0.001
  Opt_params.intrafragment_step_limit_min = i / 1000.0;
  i = rem_read(REM_GEOM_OPT2_INTRAFRAG_STEP_LIMIT_MAX); // default is 1000 -> 1.0
  Opt_params.intrafragment_step_limit_max = i / 1000.0;
  i = rem_read(REM_GEOM_OPT2_INTERFRAG_STEP_LIMIT);
  Opt_params.interfragment_step_limit     = i / 1000.0; // default is  400 -> 0.4

  // Reduce step size to ensure convergence of back-transformation of internal coordinate
  // step to cartesians.
  Opt_params.ensure_bt_convergence = rem_read("REM_GEOM_OPT2_ENSURE_BT_CONVERGENCE");

// follow root   (default 0)
  Opt_params.rfo_follow_root = rem_read(REM_GEOM_OPT2_RFO_FOLLOW_ROOT);

// which root    (default is 0 for minimum) - should be corrected in qchem default
  Opt_params.rfo_root = rem_read(REM_GEOM_OPT2_RFO_ROOT);

// scale = i / 10   (default 13)
  i = rem_read(REM_GEOM_OPT2_SCALE_CONNECTIVITY);
  Opt_params.scale_connectivity = i / 10.0;

// scale = i / 10   (default 18) // not yet implemented in QChem
  i = rem_read(REM_GEOM_OPT2_INTERFRAGMENT_SCALE_CONNECTIVITY);
  Opt_params.interfragment_scale_connectivity = i/ 10;

// multi-mode (0=single ; 1= multi) (default 0)
  i = rem_read(REM_GEOM_OPT2_FRAGMENT_MODE);
  if (i == 0)      Opt_params.fragment_mode = OPT_PARAMS::SINGLE;
  else if (i == 1) Opt_params.fragment_mode = OPT_PARAMS::MULTI;

  i = rem_read(REM_GEOM_OPT2_INTERFRAGMENT_MODE);
  if (i == 0)      Opt_params.interfragment_mode = OPT_PARAMS::FIXED;
  else if (i == 1) Opt_params.interfragment_mode = OPT_PARAMS::PRINCIPAL_AXES;

// only generate intcos
  Opt_params.intcos_generate_exit = rem_read(REM_GEOM_OPT2_GENERATE_INTCOS_ONLY);

// model 0=FISCHER ; 1 = SCHLEGEL (default 0) ; 2 = simple
  i = rem_read(REM_GEOM_OPT2_INTRAFRAGMENT_H);
  if (i == 0)      Opt_params.intrafragment_H = OPT_PARAMS::FISCHER;
  else if (i == 1) Opt_params.intrafragment_H = OPT_PARAMS::SCHLEGEL;
  else if (i == 2) Opt_params.intrafragment_H = OPT_PARAMS::SIMPLE;
  else if (i == 3) Opt_params.intrafragment_H = OPT_PARAMS::LINDH;

// interfragment 0=DEFAULT ; 1=FISCHER_LIKE
  i = rem_read(REM_GEOM_OPT2_INTERFRAGMENT_H);
  if (i == 0)      Opt_params.interfragment_H = OPT_PARAMS::DEFAULT;
  else if (i == 1) Opt_params.interfragment_H = OPT_PARAMS::FISCHER_LIKE;

// Whether to freeze all fragments rigid (default 0);
  Opt_params.freeze_intrafragment = rem_read(REM_GEOM_OPT2_FREEZE_INTRAFRAGMENT);

// Needs added to QChem
  Opt_params.H_guess_every = rem_read("REM_GEOM_OPT2_H_GUESS_EVERY");

// Whether to freeze all interfragment modes:
  Opt_params.freeze_interfragment = rem_read(REM_GEOM_OPT2_FREEZE_INTERFRAGMENT);

// write final step ; default (false);
  Opt_params.write_final_step_geometry = rem_read(REM_GEOM_OPT2_WRITE_FINAL_STEP_GEOMETRY);

// {NONE, BFGS, MS, POWELL, BOFILL} (default 1)
  i = rem_read(REM_GEOM_OPT2_H_UPDATE);
  if (i == 0) Opt_params.H_update = OPT_PARAMS::NONE;
  else if (i == 1) Opt_params.H_update = OPT_PARAMS::BFGS;
  else if (i == 2) Opt_params.H_update = OPT_PARAMS::MS;
  else if (i == 3) Opt_params.H_update = OPT_PARAMS::POWELL;
  else if (i == 4) Opt_params.H_update = OPT_PARAMS::BOFILL;

// previous steps to use ; (0=all) ; default (6)
  Opt_params.H_update_use_last = rem_read(REM_GEOM_OPT2_H_UPDATE_USE_LAST);

// limit hessian changes (default true)
  Opt_params.H_update_limit = rem_read(REM_GEOM_OPT2_H_UPDATE_LIMIT);

// scale is i / 10 (default 5 -> 0.5)  ; max is i / 10 (default 10 -> 1.0)
  i = rem_read(REM_GEOM_OPT2_H_UPDATE_LIMIT_SCALE);
  Opt_params.H_update_limit_scale = i / 10;

  i = rem_read(REM_GEOM_OPT2_H_UPDATE_LIMIT_MAX);
  Opt_params.H_update_limit_max  = i / 10;

// use 1/R(AB) ; (default 0)
  Opt_params.interfragment_distance_inverse = rem_read(REM_GEOM_OPT2_INTERFRAGMENT_DISTANCE_INVERSE);

// H-bond distance = i / 10 ; default (43 -> 4.3)
  i = rem_read(REM_GEOM_OPT2_MAXIMUM_H_BOND_DISTANCE);
  Opt_params.maximum_H_bond_distance = i / 10;

// QCHEM optimization criteria ; names adapted from old QCHEM
  i = rem_read(REM_GEOM_OPT2_TOL_GRADIENT);
  Opt_params.conv_max_force = i / 1.0e6; // default (300 -> 3e-4)

  i = rem_read(REM_GEOM_OPT2_TOL_DISPLACEMENT);
  Opt_params.conv_max_disp  = i / 1.0e6; // default (1200 -> 1.2e-3)

  i = rem_read(REM_GEOM_OPT2_TOL_ENERGY);
  Opt_params.conv_max_DE    = i / 1.0e8; // default (100 -> 1.0e-6)

// Turn "on" these convergence criteria; At least for now, QChem, doesn't
// support all the special string names for convergence criteria.
  Opt_params.i_untampered = true; // allow flex between force and displacement
  Opt_params.i_max_force = true;
  Opt_params.i_max_disp = true;
  Opt_params.i_max_DE = true;

// test B (default 0)
  Opt_params.test_B = rem_read(REM_GEOM_OPT2_TEST_B);
  Opt_params.test_derivative_B = rem_read(REM_GEOM_OPT2_TEST_DERIVATIVE_B);

// (default 1)
  Opt_params.print_lvl = rem_read(REM_GEOM_OPT2_PRINT_LVL);

  Opt_params.print_trajectory_xyz_file = rem_read(REM_GEOM_OPT2_PRINT_TRAJECTORY_XYZ_FILE);

// read Hessian (default 0)
  Opt_params.read_cartesian_H = rem_read(REM_GEOM_OPT2_READ_CARTESIAN_H);

// This optimizer will not work unless only EFP fragments are present
// Last I tried, I can't even get geometry data when running EFP_opt.in
  Opt_params.fb_fragments = rem_read(REM_EFP);

// are ONLY EFP fragments present
  if(Opt_params.fb_fragments)
    Opt_params.fb_fragments_only = rem_read(REM_EFP_FRAGMENTS_ONLY);
  else {
    Opt_params.fb_fragments_only = false;
  }

  // for coordinates with user-specified equilibrium values - this is the force constant
  //i = rem_read(REM_FIXED_COORD_FORCE_CONSTANT);
  //Opt_params.fixed_coord_force_constant = i / 10; // default (20 -> 2 au)

  Opt_params.consecutive_backsteps_allowed = rem_read(REM_GEOM_OPT2_CONSECUTIVE_BACKSTEPS);

  // if steepest-descent, then make much larger default
  if (Opt_params.step_type == OPT_PARAMS::SD && RemUninitialized(REM_GEOM_OPT2_CONSECUTIVE_BACKSTEPS))
    Opt_params.consecutive_backsteps_allowed = 1;

//TO DO: initialize IRC_step_size for Q-Chem

#endif

// Strings that carry user-specified constraints
// "frozen" means unchanging, while "fixed" means eq. value is specified
#if defined(OPTKING_PACKAGE_PSI)
  Opt_params.frozen_distance_str = options.get_str("FROZEN_DISTANCE");
  Opt_params.frozen_bend_str     = options.get_str("FROZEN_BEND");
  Opt_params.frozen_dihedral_str = options.get_str("FROZEN_DIHEDRAL");
  Opt_params.frozen_cartesian_str = options.get_str("FROZEN_CARTESIAN");

  Opt_params.fixed_distance_str = options.get_str("FIXED_DISTANCE");
  Opt_params.fixed_bend_str     = options.get_str("FIXED_BEND");
  Opt_params.fixed_dihedral_str = options.get_str("FIXED_DIHEDRAL");

  if (!Opt_params.fixed_distance_str.empty() ||
      !Opt_params.fixed_bend_str.empty()     ||
      !Opt_params.fixed_dihedral_str.empty()) {
    if (!options["INTRAFRAG_STEP_LIMIT"].has_changed())     Opt_params.intrafragment_step_limit = 0.1;
    if (!options["INTRAFRAG_STEP_LIMIT_MIN"].has_changed()) Opt_params.intrafragment_step_limit_min = 0.1;
    if (!options["INTRAFRAG_STEP_LIMIT_MAX"].has_changed()) Opt_params.intrafragment_step_limit_max = 0.1;
  }

#elif defined(OPTKING_PACKAGE_QCHEM)
  // Read QChem input and write all the frozen distances into a string
  if (rem_read(REM_GEOM_OPT2_FROZEN_DISTANCES) > 0) {
    INTEGER n_frozen = rem_read(REM_GEOM_OPT2_FROZEN_DISTANCES);
    double* fd = init_array(2*n_frozen);
    FileMan(FM_READ,FILE_FROZEN_DISTANCES,FM_DP,2*n_frozen,0,FM_BEG,fd);

    std::stringstream atoms;
    for (int i=0; i<2*n_frozen; ++i)
      atoms << (int) fd[i] << ' ';
    Opt_params.frozen_distance_str = atoms.str();

    free_array(fd);
  }
  // Read QChem input and write all the frozen bends into a string
  if (rem_read(REM_GEOM_OPT2_FROZEN_BENDS) > 0) {
    INTEGER n_frozen = rem_read(REM_GEOM_OPT2_FROZEN_BENDS);
    double* fd = init_array(3*n_frozen);
    FileMan(FM_READ,FILE_FROZEN_BENDS,FM_DP,3*n_frozen,0,FM_BEG,fd);

    std::stringstream atoms;
    for (int i=0; i<3*n_frozen; ++i)
      atoms << (int) fd[i] << ' ';
    Opt_params.frozen_bend_str = atoms.str();

    free_array(fd);
  }
  // Read QChem input and write all the frozen dihedrals into a string
  if (rem_read(REM_GEOM_OPT2_FROZEN_DIHEDRALS) > 0) {
    INTEGER n_frozen = rem_read(REM_GEOM_OPT2_FROZEN_DIHEDRALS);
    double* fd = init_array(4*n_frozen);
    FileMan(FM_READ,FILE_FROZEN_DIHEDRALS,FM_DP,4*n_frozen,0,FM_BEG,fd);

    std::stringstream atoms;
    for (int i=0; i<4*n_frozen; ++i)
      atoms << (int) fd[i] << ' ';
    Opt_params.frozen_dihedral_str = atoms.str();

    free_array(fd);
  }
  // Read QChem input and write all the frozen dihedrals into a string
  if (rem_read(REM_GEOM_OPT2_FROZEN_CARTESIANS) > 0) {
    INTEGER n_frozen = rem_read(REM_GEOM_OPT2_FROZEN_CARTESIANS);
    double* fd = init_array(2*n_frozen);
    FileMan(FM_READ,FILE_FROZEN_CARTESIANS,FM_DP,2*n_frozen,0,FM_BEG,fd);

    // TODO: will have to add code for "xyz" format to integer format for
    // subsequent reads of frozen cartesians from file.
    std::stringstream atoms;
    for (int i=0; i<2*n_frozen; ++i)
      atoms << (int) fd[i] << ' ';
    Opt_params.frozen_cartesian_str = atoms.str();

    free_array(fd);
  }

#endif

// ** Items are below unlikely to need modified

// how close to pi should a torsion be to assume it may have passed through 180
  Opt_params.fix_tors_near_pi = _pi / 2;

// torsional angles will not be computed if the contained bond angles are within
// this many radians of zero or 180. (< ~1 and > ~179 degrees)
  Opt_params.tors_angle_lim = 0.017;

// only used for determining which atoms in a fragment are acceptable for use
// as reference atoms.  We avoid collinear sets.
// angle is 0/pi if the bond angle is within this fraction of pi from 0/pi
  Opt_params.interfrag_collinear_tol = 0.01;

// cos(torsional angle) must be this close to -1/+1 for angle to count as 0/pi
  Opt_params.tors_cos_tol = 1e-10;

// if bend exceeds this value, then also create linear bend complement
  Opt_params.linear_bend_threshold = 3.05; // about 175 degrees

// If bend is smaller than this value, then never fix its associated vectors
// this allows iterative steps through and near zero degrees.
  Opt_params.small_bend_fix_threshold = 0.35;

// threshold for which entries in diagonalized redundant matrix are kept and inverted
// while computing a generalized inverse of a matrix
  Opt_params.redundant_eval_tol = 1.0e-10;

// Parameters that control the backtransformation to cartesian coordinates
  Opt_params.bt_max_iter = 25;
  Opt_params.bt_dx_conv = 1.0e-6;
  Opt_params.bt_dx_conv_rms_change = 1.0e-12;
  //Opt_params.bt_dx_conv = 1.0e-10;
  //Opt_params.bt_dx_conv_rms_change = 1.0e-14;


// Hessian update is avoided if any internal coordinate has changed by more than this in radians/au
  Opt_params.H_update_dq_tol = 0.5;

// Some parameter error-checking / modification
  if (Opt_params.fb_fragments_only) {
    Opt_params.test_B = false;
    Opt_params.test_derivative_B = false;
  }

/*  if dynamic mode is on, then other settings are overridden.
* step_type = step
* intrafragment_step_limit = step_limit
* consecutive_backsteps = backsteps
* RI = redundant internals; D = default anyway

*dynamic  step   coord   step_limit      backsteps              criteria
* level                                               for downmove    for upmove
*  0      RFO    RI      dynamic         no           none            none
*
*  1      RFO    RI      dynamic(D)      no           1 bad step
*
*  2      RFO    RI      small initial   yes (1)      1 bad step
*                        dynamic(D)
*
*  3      SD     RI      large(D)        yes (1)      1 bad step
*
*  4      SD     RI+XYZ  large(D)        yes (1)      1 bad step
*
*  5      SD     XYZ     large(D)        yes (1)      1 bad step
*
*  6      SD     XYZ     small           yes (1)      1 bad step
*
*  7  abort
*
*  BackStep:
*   DE > 0 in minimization
*
*  BadStep:
*   DE > 0 and backsteps exceeded and iterations > 5  ** OR **
*   badly defined internal coordinate or derivative
*
* */
  if (options.get_int("DYNAMIC_LEVEL") == 0) // if 0, then not dynamic
    Opt_params.dynamic = 0;
  else if (INTCO_EXCEPT::dynamic_level != 0) // already set
    Opt_params.dynamic = INTCO_EXCEPT::dynamic_level;
  else
    Opt_params.dynamic = options.get_int("DYNAMIC_LEVEL");

  Opt_params.sd_hessian = 1.0; // small step

  switch(Opt_params.dynamic) {
    case 0: // not dynamic
      break;
    case 1:
      Opt_params.coordinates = OPT_PARAMS::REDUNDANT;
      Opt_params.consecutive_backsteps_allowed = 0;
      Opt_params.step_type = OPT_PARAMS::RFO;
             printf("At level 1: Red. Int., RFO, no backsteps, dynamic trust\n");
      oprintf_out("\tAt level 1: Red. Int., RFO, no backsteps, dynamic trust\n");
      break;
    case 2:
      Opt_params.coordinates = OPT_PARAMS::REDUNDANT;
      Opt_params.consecutive_backsteps_allowed = 1;
      Opt_params.step_type = OPT_PARAMS::RFO;
      Opt_params.intrafragment_step_limit = 0.2;
      Opt_params.intrafragment_step_limit_min = 0.2; //this code overwrites changes anyway
      Opt_params.intrafragment_step_limit_max = 0.2;
             printf("At level 2: Red. Int., RFO, backsteps, smaller trust.\n");
      oprintf_out("\tAt level 2: Red. Int., RFO, backsteps, smaller trust.\n");
      break;
    case 3:
      Opt_params.coordinates = OPT_PARAMS::BOTH;
      Opt_params.consecutive_backsteps_allowed = 1;
      Opt_params.step_type = OPT_PARAMS::RFO;
      Opt_params.intrafragment_step_limit = 0.1;
      Opt_params.intrafragment_step_limit_min = 0.1; //this code overwrites changes anyway
      Opt_params.intrafragment_step_limit_max = 0.1;
             printf("At level 3: Red. Int. + XYZ, RFO, backsteps, smaller trust.\n");
      oprintf_out("\tAt level 3: Red. Int. + XYZ, RFO, backsteps, smaller trust.\n");
      break;
    case 4:
      Opt_params.coordinates = OPT_PARAMS::CARTESIAN;
      Opt_params.consecutive_backsteps_allowed = 1;
      Opt_params.step_type = OPT_PARAMS::RFO;
      Opt_params.intrafragment_H = OPT_PARAMS::LINDH;
      Opt_params.intrafragment_step_limit = 0.3;
      Opt_params.intrafragment_step_limit_min = 0.3; //this code overwrites changes anyway
      Opt_params.intrafragment_step_limit_max = 0.3;
             printf("At level 4: XYZ, RFO, backsteps, larger trust.\n");
      oprintf_out("\tAt level 4: XYZ, RFO, backsteps, larger trust.\n");
      break;
    case 5:   // Try repeating level 4 which is working well
      Opt_params.coordinates = OPT_PARAMS::CARTESIAN;
      Opt_params.consecutive_backsteps_allowed = 1;
      Opt_params.step_type = OPT_PARAMS::RFO;
      Opt_params.intrafragment_H = OPT_PARAMS::LINDH;
      Opt_params.intrafragment_step_limit = 0.2;
      Opt_params.intrafragment_step_limit_min = 0.2; //this code overwrites changes anyway
      Opt_params.intrafragment_step_limit_max = 0.2;
             printf("At level 5: XYZ, RFO, backsteps, medium trust.\n");
      oprintf_out("\tAt level 5: XYZ, RFO, backsteps, medium trust.\n");
/* Opt_params.coordinates = OPT_PARAMS::CARTESIAN;
      Opt_params.consecutive_backsteps_allowed = 1;
      Opt_params.step_type = OPT_PARAMS::SD;
      Opt_params.sd_hessian = 0.1;
      Opt_params.intrafragment_step_limit = 0.3;
      Opt_params.intrafragment_step_limit_min = 0.3; //this code overwrites changes anyway
      Opt_params.intrafragment_step_limit_max = 0.3;
             printf("At level 5: XYZ, SD, backsteps, larger trust.\n");
      oprintf_out("\tAt level 5: XYZ, SD, backsteps, larger trust.\n"); */
      break;
    case 6:
      Opt_params.coordinates = OPT_PARAMS::CARTESIAN;
      Opt_params.consecutive_backsteps_allowed = 1;
      Opt_params.step_type = OPT_PARAMS::SD;
      Opt_params.sd_hessian = 0.3;
      Opt_params.intrafragment_step_limit = 0.3;
      Opt_params.intrafragment_step_limit_min = 0.3; //this code overwrites changes anyway
      Opt_params.intrafragment_step_limit_max = 0.3;
             printf("At level 5: XYZ, SD, backsteps, larger trust.\n");
      oprintf_out("\tAt level 5: XYZ, SD, backsteps, larger trust.\n");
      break;
    case 7:
      Opt_params.coordinates = OPT_PARAMS::CARTESIAN;
      Opt_params.consecutive_backsteps_allowed = 1;
      Opt_params.step_type = OPT_PARAMS::SD;
      Opt_params.sd_hessian = 0.6;
      Opt_params.intrafragment_step_limit = 0.1;
      Opt_params.intrafragment_step_limit_min = 0.1; //this code overwrites changes anyway
      Opt_params.intrafragment_step_limit_max = 0.1;
             printf("Moving to level 6: XYZ, SD, backsteps, small trust, smaller steps.\n");
      oprintf_out("\tMoving to level 6: XYZ, SD, backsteps, small trust, smaller steps.\n");
      break;
  default:
      oprintf_out("Unknown value of Opt_params.dynamic variable.\n");
  }

  if (Opt_params.print_lvl > 1 || Opt_params.print_params) print_params_out();

}

void print_params_out(void) {

  oprintf_out( "dynamic level          = %18d\n", Opt_params.dynamic);
  oprintf_out( "conv_max_force         = %18.2e\n", Opt_params.conv_max_force);
  oprintf_out( "conv_rms_force         = %18.2e\n", Opt_params.conv_rms_force);
  oprintf_out( "conv_max_DE            = %18.2e\n", Opt_params.conv_max_DE);
  oprintf_out( "conv_max_disp          = %18.2e\n", Opt_params.conv_max_disp);
  oprintf_out( "conv_rms_disp          = %18.2e\n", Opt_params.conv_rms_disp);
  oprintf_out( "SD Hessian             = %18.2e\n", Opt_params.sd_hessian);

  oprintf_out( "scale_connectivity     = %18.2e\n", Opt_params.scale_connectivity);
  oprintf_out( "interfragment_scale_connectivity = %18.2e\n",
    Opt_params.interfragment_scale_connectivity);
  //oprintf_out( "fixed_coord_force_constant = %14.2e\n", Opt_params.fixed_coord_force_constant);

  if (Opt_params.fragment_mode == OPT_PARAMS::SINGLE)
  oprintf_out( "fragment_mode          = %18s\n", "single");
  else if (Opt_params.fragment_mode == OPT_PARAMS::MULTI)
  oprintf_out( "fragment_mode          = %18s\n", "multi");

  if (Opt_params.interfragment_mode == OPT_PARAMS::FIXED)
  oprintf_out( "interfragment_mode        = %18s\n", "fixed");
  else if (Opt_params.interfragment_mode == OPT_PARAMS::PRINCIPAL_AXES)
  oprintf_out( "interfragment_mode        = %18s\n", "principal axes");

  for (int i=0; i<(int) Opt_params.frag_ref_atoms.size(); ++i) {
      if (i == 0) oprintf_out( "Reference points specified for fragments:\n");
      oprintf_out( "Fragment %d\n", i);
      for (int j=0; j<(int) Opt_params.frag_ref_atoms[i].size(); ++j) {
          oprintf_out( "Reference atom %d: ", j);
          for (int k=0; k<(int) Opt_params.frag_ref_atoms[i][j].size(); ++k)
              oprintf_out( "%d ", Opt_params.frag_ref_atoms[i][j][k]);
          oprintf_out("\n");
      }
  }

  if (Opt_params.intcos_generate_exit)
    oprintf_out( "intcos_generate_exit   = %18s\n", "true");
  else
    oprintf_out( "intcos_generate_exit   = %18s\n", "false");

  //if (Opt_params.print_params)
  //oprintf_out( "print_params           = %18s\n", "true");
  //else
  //oprintf_out( "print_params           = %18s\n", "false");

  oprintf_out( "print_params           = %18s\n", Opt_params.print_params ? "true" : "false");

  oprintf_out( "print_lvl              = %d\n", Opt_params.print_lvl);

  if (Opt_params.ensure_bt_convergence)
    oprintf_out("ensure_bt_convergence = %17s\n", "true");
  else
    oprintf_out("ensure_bt_convergence = %17s\n", "false");

  if (Opt_params.rfo_follow_root)
  oprintf_out( "rfo_follow_root        = %18s\n", "true");
  else
  oprintf_out( "rfo_follow_root        = %18s\n", "false");

  oprintf_out( "rfo_root               = %18d\n", Opt_params.rfo_root);

  oprintf_out( "rfo_normalization_max  = %18.2e\n", Opt_params.rfo_normalization_max);
  oprintf_out( "rsrfo_alpha_max        = %18.3e\n", Opt_params.rsrfo_alpha_max);

  if (Opt_params.step_type == OPT_PARAMS::NR)
  oprintf_out( "step_type              = %18s\n", "N-R");
  else if (Opt_params.step_type == OPT_PARAMS::RFO)
  oprintf_out( "step_type              = %18s\n", "RFO");
  else if (Opt_params.step_type == OPT_PARAMS::P_RFO)
  oprintf_out( "step_type              = %18s\n", "P_RFO");
  else if (Opt_params.step_type == OPT_PARAMS::LINESEARCH_STATIC)
  oprintf_out( "step_type              = %18s\n", "Static linesearch");

  if (Opt_params.coordinates == OPT_PARAMS::REDUNDANT)
  oprintf_out( "opt. coordinates       = %18s\n", "Redundant Internals");
  else if (Opt_params.coordinates == OPT_PARAMS::DELOCALIZED)
  oprintf_out( "opt. coordinates       = %18s\n", "Delocalized");
  else if (Opt_params.coordinates == OPT_PARAMS::NATURAL)
  oprintf_out( "opt. coordinates       = %18s\n", "Natural");
  else if (Opt_params.coordinates == OPT_PARAMS::CARTESIAN)
  oprintf_out( "opt. coordinates       = %18s\n", "Cartesian");
  else if (Opt_params.coordinates == OPT_PARAMS::BOTH)
  oprintf_out( "opt. coordinates       = %18s\n", "Add Cartesians");

  oprintf_out( "linesearch_static_N    = %18d\n", Opt_params.linesearch_static_N);
  oprintf_out( "linesearch_static_min  = %18.3e\n", Opt_params.linesearch_static_min);
  oprintf_out( "linesearch_static_max  = %18.3e\n", Opt_params.linesearch_static_max);

  oprintf_out( "consecutive_backsteps  = %18d\n",  Opt_params.consecutive_backsteps_allowed);

  if (Opt_params.intrafragment_H == OPT_PARAMS::FISCHER)
  oprintf_out( "intrafragment_H        = %18s\n", "Fischer");
  else if (Opt_params.intrafragment_H == OPT_PARAMS::SCHLEGEL)
  oprintf_out( "intrafragment_H        = %18s\n", "Schlegel");
  else if (Opt_params.intrafragment_H == OPT_PARAMS::SIMPLE)
  oprintf_out( "intrafragment_H        = %18s\n", "Simple");
  else if (Opt_params.intrafragment_H == OPT_PARAMS::LINDH)
  oprintf_out( "intrafragment_H        = %18s\n", "Lindh");
  else if (Opt_params.intrafragment_H == OPT_PARAMS::LINDH_SIMPLE)
  oprintf_out( "intrafragment_H        = %18s\n", "Lindh - Simple");

  if (Opt_params.interfragment_H == OPT_PARAMS::DEFAULT)
  oprintf_out( "interfragment_H        = %18s\n", "Default");
  else if (Opt_params.interfragment_H == OPT_PARAMS::FISCHER_LIKE)
  oprintf_out( "interfragment_H        = %18s\n", "Fischer_like");

  if (Opt_params.H_update == OPT_PARAMS::NONE)
  oprintf_out( "H_update               = %18s\n", "None");
  else if (Opt_params.H_update == OPT_PARAMS::BFGS)
  oprintf_out( "H_update               = %18s\n", "BFGS");
  else if (Opt_params.H_update == OPT_PARAMS::MS)
  oprintf_out( "H_update               = %18s\n", "MS");
  else if (Opt_params.H_update == OPT_PARAMS::POWELL)
  oprintf_out( "H_update               = %18s\n", "Powell");
  else if (Opt_params.H_update == OPT_PARAMS::BOFILL)
  oprintf_out( "H_update               = %18s\n", "Bofill");

  oprintf_out( "H_update_use_last      = %18d\n", Opt_params.H_update_use_last);

  if (Opt_params.freeze_intrafragment)
  oprintf_out( "freeze_intrafragment   = %18s\n", "true");
  else
  oprintf_out( "freeze_intrafragment   = %18s\n", "false");

  oprintf_out( "intrafragment_step_limit=%18.2e\n", Opt_params.intrafragment_step_limit);

  oprintf_out( "interfragment_step_limit=%18.2e\n", Opt_params.interfragment_step_limit);

  if (Opt_params.add_auxiliary_bonds)
    oprintf_out( "add_auxiliary_bonds   = %18s\n", "true");
  else
    oprintf_out( "add_auxiliary_bonds   = %18s\n", "false");

  if (Opt_params.H_guess_every)
    oprintf_out( "H_guess_every         = %18s\n", "true");
  else
    oprintf_out( "H_guess_every         = %18s\n", "false");

  oprintf_out( "auxiliary_bond_factor =%18.2e\n", Opt_params.auxiliary_bond_factor);

  if (Opt_params.H_update_limit)
    oprintf_out( "H_update_limit         = %18s\n", "true");
  else
    oprintf_out( "H_update_limit         = %18s\n", "false");

  oprintf_out( "H_update_limit_scale   = %18.2e\n", Opt_params.H_update_limit_scale);
  oprintf_out( "H_update_limit_max     = %18.2e\n", Opt_params.H_update_limit_max);
  oprintf_out( "H_update_den_tol       = %18.2e\n", Opt_params.H_update_den_tol);

  if (Opt_params.interfragment_distance_inverse)
  oprintf_out( "interfragment_distance_inverse=%12s\n", "true");
  else
  oprintf_out( "interfragment_distance_inverse=%12s\n", "false");

  if (Opt_params.write_final_step_geometry)
  oprintf_out( "write_final_step_geometry= %16s\n", "true");
  else
  oprintf_out( "write_final_step_geometry= %16s\n", "false");

  oprintf_out( "maximum_H_bond_distance= %18.2e\n", Opt_params.maximum_H_bond_distance);

  if (Opt_params.read_cartesian_H)
  oprintf_out( "read_cartesian_H       = %18s\n", "true");
  else
  oprintf_out( "read_cartesian_H       = %18s\n", "false");

  if (Opt_params.fb_fragments)
  oprintf_out( "fb_fragments          = %18s\n", "true");
  else
  oprintf_out( "fb_fragments          = %18s\n", "false");

  if (Opt_params.fb_fragments_only)
  oprintf_out( "fb_fragments_only     = %18s\n", "true");
  else
  oprintf_out( "fb_fragments_only     = %18s\n", "false");

  oprintf_out( "frozen_distance: \n");
  if (!Opt_params.frozen_distance_str.empty())
    oprintf_out( "%s\n", Opt_params.frozen_distance_str.c_str());

  oprintf_out( "frozen_bend: \n");
  if (!Opt_params.frozen_bend_str.empty())
    oprintf_out( "%s\n", Opt_params.frozen_bend_str.c_str());

  oprintf_out( "frozen_dihedral: \n");
  if (!Opt_params.frozen_dihedral_str.empty())
    oprintf_out( "%s\n", Opt_params.frozen_dihedral_str.c_str());

  oprintf_out( "frozen_cartesian: \n");
  if (!Opt_params.frozen_cartesian_str.empty())
    oprintf_out( "%s\n", Opt_params.frozen_cartesian_str.c_str());

  oprintf_out( "fixed_distance: \n");
  if (!Opt_params.fixed_distance_str.empty())
    oprintf_out( "%s\n", Opt_params.fixed_distance_str.c_str());

  oprintf_out( "fixed_bend: \n");
  if (!Opt_params.fixed_bend_str.empty())
    oprintf_out( "%s\n", Opt_params.fixed_bend_str.c_str());

  oprintf_out( "fixed_dihedral: \n");
  if (!Opt_params.fixed_dihedral_str.empty())
    oprintf_out( "%s\n", Opt_params.fixed_dihedral_str.c_str());

  if (Opt_params.print_trajectory_xyz_file)
    oprintf_out("print_trajectory_xyz_file = %18s\n", "true");
  else
    oprintf_out("print_trajectory_xyz_file = %18s\n", "false");
}

}

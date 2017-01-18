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

/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/

#include <string>

namespace psi { namespace ccdensity {

/* Input parameters for cclambda */
struct Params {
  double tolerance;
  long int memory;
  int cachelev;
  int aobasis;
  int ref;
  int onepdm; /* produce ONLY the onepdm for properties */
  int onepdm_grid_dump; // dump the onepdm on a grid to a dx file
  int relax_opdm;
  int use_zeta;
  int calc_xi;
  int connect_xi;
  int restart;
  int ground;
  int transition; 
  int dertype;
  std::string wfn;
  int nstates;
  int prop_sym;
  int prop_root;
  int prop_all;
  std::string gauge;
  bool write_nos;

  /* these are used by Xi and twopdm code */
  int G_irr; 
  int R_irr; 
  int L_irr;
  double R0;   
  double L0;
  int ael;
  double cceom_energy;
  double overlap1; /* <L1|R1> */
  double overlap2; /* <L2|R2> */
  double RD_overlap; /* Rmnef <mn||ef> */
  double RZ_overlap; /* <R|zeta> */
};

struct RHO_Params {
  int L_irr;  
  int R_irr; 
  int G_irr;
  int L_root; 
  int R_root;
  int L_ground; 
  int R_ground;
  double R0;   
  double L0;
  char L1A_lbl[32];
  char L1B_lbl[32];
  char L2AA_lbl[32];
  char L2BB_lbl[32];
  char L2AB_lbl[32];
  char R1A_lbl[32];
  char R1B_lbl[32];
  char R2AA_lbl[32];
  char R2BB_lbl[32];
  char R2AB_lbl[32];
  double cceom_energy;
  double overlap1; /* <L1|R1> */
  double overlap2; /* <L2|R2> */
  double RD_overlap; /* Rmnef <mn||ef> */
  char DIJ_lbl[10];
  char Dij_lbl[10];
  char DAB_lbl[10];
  char Dab_lbl[10];
  char DIA_lbl[10];
  char Dia_lbl[10];
  char DAI_lbl[10];
  char Dai_lbl[10];
  char opdm_lbl[32];
  char opdm_a_lbl[32];
  char opdm_b_lbl[32];
};

struct TD_Params {
  int irrep;
  int root;
  double R0;
  double cceom_energy;
  char L1A_lbl[32];
  char L1B_lbl[32];
  char L2AA_lbl[32];
  char L2BB_lbl[32];
  char L2AB_lbl[32];
  char R1A_lbl[32];
  char R1B_lbl[32];
  char R2AA_lbl[32];
  char R2BB_lbl[32];
  char R2AB_lbl[32];
  double OS;
  double RS_length;
  double RS_velocity;
  double einstein_a;
  double einstein_b;
};

struct XTD_Params {
  int irrep1;
  int irrep2;
  int root1;
  int root2;
  double cceom_energy;
  double OS;
  double RS_length;
  double RS_velocity;
  double einstein_a;
  double einstein_b;
};


}} // namespace psi::ccdensity
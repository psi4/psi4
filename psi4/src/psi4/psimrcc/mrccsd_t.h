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

#ifndef _psi_src_bin_psimrcc_mrccsd_t_h_
#define _psi_src_bin_psimrcc_mrccsd_t_h_

#include <vector>

namespace psi{ namespace psimrcc{

class CCIndex;
class BlockMatrix;
class IndexMatrix;
class Hamiltonian;

class MRCCSD_T
{
enum TriplesAlgorithm {UnrestrictedTriples,RestrictedTriples,SpinAdaptedTriples};
public:
  // Constructor and destructor
  MRCCSD_T(Options& options,Hamiltonian* h_eff_);
  ~MRCCSD_T();
private:
  void startup();
  void cleanup();
  void compute();
  void compute_ooo_triples();
  void compute_OOO_triples();
  void compute_ooO_triples();
  void compute_oOO_triples();

  void compute_restricted();
  void compute_ooo_triples_restricted();
  void compute_ooO_triples_restricted();
  void compute_oOO_triples_restricted();
  void compute_OOO_triples_restricted();

  void compute_spin_adapted();
  void compute_ooo_triples_spin_adapted();
  void compute_ooO_triples_spin_adapted();
  void compute_oOO_triples_spin_adapted();
  void compute_OOO_triples_spin_adapted();

  void compute_ooo_contribution_to_Heff(int i,int j,int k,int mu,BlockMatrix* T3);
  void compute_ooO_contribution_to_Heff(int i,int j,int k,int mu,BlockMatrix* T3);
  void compute_oOO_contribution_to_Heff(int i,int j,int k,int mu,BlockMatrix* T3);
  void compute_OOO_contribution_to_Heff(int i,int j,int k,int mu,BlockMatrix* T3);

  double compute_A_ooo_contribution_to_Heff(int u_abs,int x_abs,int i_abs,int j_abs,int k_abs,int mu,BlockMatrix* T3);
  double compute_A_ooO_contribution_to_Heff(int u_abs,int x_abs,int i_abs,int j_abs,int k_abs,int mu,BlockMatrix* T3);
  double compute_A_oOO_contribution_to_Heff(int u_abs,int x_abs,int i_abs,int j_abs,int k_abs,int mu,BlockMatrix* T3);

  double compute_B_ooO_contribution_to_Heff(int U_abs,int X_abs,int i_abs,int j_abs,int k_abs,int mu,BlockMatrix* T3);
  double compute_B_oOO_contribution_to_Heff(int U_abs,int X_abs,int i_abs,int j_abs,int k_abs,int mu,BlockMatrix* T3);
  double compute_B_OOO_contribution_to_Heff(int U_abs,int X_abs,int i_abs,int j_abs,int k_abs,int mu,BlockMatrix* T3);

  double compute_AB_ooO_contribution_to_Heff(int u_abs,int V_abs,int x_abs,int Y_abs,int i_abs,int j_abs,int k_abs,int mu,BlockMatrix* T3);
  double compute_AB_oOO_contribution_to_Heff(int u_abs,int V_abs,int x_abs,int Y_abs,int i_abs,int j_abs,int k_abs,int mu,BlockMatrix* T3);

  void compute_ooo_contribution_to_Heff_restricted(int i,int j,int k,int mu,BlockMatrix* T3);
  void compute_ooO_contribution_to_Heff_restricted(int i,int j,int k,int mu,BlockMatrix* T3);
  void compute_oOO_contribution_to_Heff_restricted(int i,int j,int k,int mu,BlockMatrix* T3);
  void compute_OOO_contribution_to_Heff_restricted(int i,int j,int k,int mu,BlockMatrix* T3);

  double compute_A_ooo_contribution_to_Heff_restricted(int u_abs,int x_abs,int i_abs,int j_abs,int k_abs,int mu,BlockMatrix* T3);
  double compute_A_ooO_contribution_to_Heff_restricted(int u_abs,int x_abs,int i_abs,int j_abs,int k_abs,int mu,BlockMatrix* T3);
  double compute_A_oOO_contribution_to_Heff_restricted(int u_abs,int x_abs,int i_abs,int j_abs,int k_abs,int mu,BlockMatrix* T3);

  double compute_B_ooO_contribution_to_Heff_restricted(int U_abs,int X_abs,int i_abs,int j_abs,int k_abs,int mu,BlockMatrix* T3);
  double compute_B_oOO_contribution_to_Heff_restricted(int U_abs,int X_abs,int i_abs,int j_abs,int k_abs,int mu,BlockMatrix* T3);
  double compute_B_OOO_contribution_to_Heff_restricted(int U_abs,int X_abs,int i_abs,int j_abs,int k_abs,int mu,BlockMatrix* T3);

  double compute_AB_ooO_contribution_to_Heff_restricted(int u_abs,int V_abs,int x_abs,int Y_abs,int i_abs,int j_abs,int k_abs,int mu,BlockMatrix* T3);
  double compute_AB_oOO_contribution_to_Heff_restricted(int u_abs,int V_abs,int x_abs,int Y_abs,int i_abs,int j_abs,int k_abs,int mu,BlockMatrix* T3);



  void form_T2_ij_a_b(IndexMatrix* T2_ij_a_b,bool spin1,bool spin2,bool transpose);
  void form_T2_i_ab_j(IndexMatrix* T2_i_ab_j,bool spin1,bool spin2,bool transpose);
  void form_V_k_bc_e(IndexMatrix* V_k_bc_e,double direct,double exchange);
  void form_V_jk_c_m(IndexMatrix* V_jk_c_m,double direct,double exchange);

  void build_W_intermediates();
  void check_intruders();

  Options &options_;

  int nirreps;
  int nrefs;

  double threshold;

  TriplesAlgorithm triples_algorithm;

  Hamiltonian* h_eff;

  std::vector<std::vector<bool> > is_aocc;
  std::vector<std::vector<bool> > is_bocc;
  std::vector<std::vector<bool> > is_avir;
  std::vector<std::vector<bool> > is_bvir;

  // Denominators
  std::vector<std::vector<double> > e_oo;
  std::vector<std::vector<double> > e_OO;
  std::vector<std::vector<double> > e_vv;
  std::vector<std::vector<double> > e_VV;

  std::vector<std::vector<double> > Mk_factor;
  std::vector<double>               Mk_shift;

  std::vector<double***> F_ov;
  std::vector<double***> F_OV;

  std::vector<double***> F2_ov;
  std::vector<double***> F2_OV;

  std::vector<double***> T1_ov;
  std::vector<double***> T1_OV;

  std::vector<double***> W_ooov;
  std::vector<double***> W_oOoV;
  std::vector<double***> W_OoOv;
  std::vector<double***> W_OOOV;

  std::vector<double***> W_vovv;
  std::vector<double***> W_vOvV;
  std::vector<double***> W_VoVv;
  std::vector<double***> W_VOVV;

  double*** V_oovv;
  double*** V_oOvV;
//  double*** V_ooov;
//  double*** V_oOoV;
//  double*** V_vovv;
//  double*** V_vOvV;

  std::vector<double***> T2_oovv;
  std::vector<double***> T2_oOvV;
  std::vector<double***> T2_OOVV;

  CCIndex* o;
  CCIndex* oo;
  CCIndex* v;
  CCIndex* vv;
  CCIndex* vvv;
  CCIndex* vo;
  CCIndex* ov;
  CCIndex* ovv;
  CCIndex* ooo;

  BlockMatrix*** Z;
  BlockMatrix*** W;
  BlockMatrix*** W_ijk;
  BlockMatrix*** W_ikj;
  BlockMatrix*** W_jki;
  BlockMatrix*** T;

  IndexMatrix* T2_ij_a_b;
  IndexMatrix* T2_iJ_a_B;
  IndexMatrix* T2_iJ_B_a;
  IndexMatrix* T2_IJ_A_B;

  IndexMatrix* T2_i_ab_j;
  IndexMatrix* T2_i_aB_J;
  IndexMatrix* T2_J_aB_i;
  IndexMatrix* T2_I_AB_J;

  IndexMatrix* V_k_bc_e;
  IndexMatrix* V_K_bC_e;
  IndexMatrix* V_k_bC_E;

  IndexMatrix* V_jk_c_m;
  IndexMatrix* V_jK_c_M;
  IndexMatrix* V_jK_C_m;

  std::vector<double> e4T;
  std::vector<double> e4ST;
  std::vector<double> e4DT;

  double E4,E4T,E4ST,E4DT;

  std::vector<double> E4T_ooo;
  std::vector<double> E4T_ooO;
  std::vector<double> E4T_oOO;
  std::vector<double> E4T_OOO;
  std::vector<double> E4ST_ooo;
  std::vector<double> E4ST_ooO;
  std::vector<double> E4ST_oOO;
  std::vector<double> E4ST_OOO;
  std::vector<double> E4DT_ooo;
  std::vector<double> E4DT_ooO;
  std::vector<double> E4DT_oOO;
  std::vector<double> E4DT_OOO;
  std::vector<double> E4_ooo;
  std::vector<double> E4_ooO;
  std::vector<double> E4_oOO;
  std::vector<double> E4_OOO;

  std::vector<std::vector<double> > d_h_eff;
};

}} /* End Namespaces */

#endif // _psi_src_bin_psimrcc_mrccsd_t_h_
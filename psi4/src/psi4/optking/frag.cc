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
   \file frag.cc
   \ingroup optking
   \brief fragment (molecule) class
*/

#include "frag.h"

#include "mem.h"
#include "v3d.h"
#include "atom_data.h"
#include "cov_radii.h"
#include "opt_data.h"
#include "psi4/optking/physconst.h"
#include "linear_algebra.h"
#include "psi4/psi4-dec.h"
#include "print.h"
#define EXTERN
#include "globals.h"

#if defined(OPTKING_PACKAGE_PSI)
 #include <cmath>
#elif defined (OPTKING_PACKAGE_QCHEM)
 #include "qcmath.h"
#endif
#include "psi4/libparallel/ParallelPrinter.h"
namespace opt {

using namespace v3d;

// Memory for Z and geom is provided by calling function
FRAG::FRAG(int natom_in, double *Z_in, double **geom_in) {
  natom = natom_in;
  Z = Z_in;
  geom = geom_in;

  frozen = false;

  grad = init_matrix(natom,3);
  connectivity = init_bool_matrix(natom,natom);
  mass = init_array(natom);
}

FRAG::FRAG(int natom_in) {
  natom = natom_in;

  frozen = false;

  Z = init_array(natom);
  geom = init_matrix(natom,3);
  grad = init_matrix(natom,3);
  connectivity = init_bool_matrix(natom,natom);
  mass = init_array(natom);
}

FRAG::~FRAG() {
  free_array(Z);
  free_matrix(geom);
  free_matrix(grad);
  free_array(mass);
  free_bool_matrix(connectivity);
  coords.clear_combos();
  for (ULI i=0; i<coords.simples.size(); ++i)
    delete coords.simples[i];
  coords.simples.clear();
}

// for now just set default masses - fix later
void FRAG::set_masses(void) {
  int i;
  for (i=0; i<natom; ++i)
    mass[i] = Z_to_mass[(int) Z[i]];
  return;
}

// automatically determine bond connectivity by comparison of interatomic distance
//with scale_connectivity * sum of covalent radii
void FRAG::update_connectivity_by_distances(void) {
  int i, j, *Zint;
  double Rij;
  double scale = Opt_params.scale_connectivity;

  Zint = new int [natom];
  for (i=0; i<natom; ++i) {
    Zint[i] = (int) Z[i];
    if ( Zint[i] > LAST_COV_RADII_INDEX )
      throw(INTCO_EXCEPT("Warning: cannot automatically bond atom with strange atomic number"));
  }

  for (i=0; i<natom; ++i)
    for (j=0; j<natom; ++j)
      connectivity[i][j] = false;

  for (i=0; i<natom; ++i) {
    for (j=0; j<i; ++j) {
      Rij = v3d_dist(geom[i], geom[j]);
      if (Rij < scale * (cov_radii[Zint[i]] + cov_radii[Zint[j]])/_bohr2angstroms)
        connectivity[i][j] = connectivity[j][i] = true;
    }
  }
  delete [] Zint;
}

//build connectivity matrix from the current set of bonds
void FRAG::update_connectivity_by_bonds(void) {
  for (int i=0; i<natom; ++i)
    for (int j=0; j<natom; ++j)
      connectivity[i][j] = false;

  for (std::size_t i=0; i<coords.simples.size(); ++i) {
    if (coords.simples.at(i)->g_type() == stre_type) {
      int a = coords.simples.at(i)->g_atom(0);
      int b = coords.simples.at(i)->g_atom(1);
      connectivity[a][b] = connectivity[b][a] = true;
    }
  }
}

// automatically add bond stretch coodinates based on connectivity matrix
// return number of added coordinates
int FRAG::add_stre_by_connectivity(void) {
  int nadded = 0;
  int i,j;
  for (i=0; i<natom; ++i) {
    for (j=i+1; j<natom; ++j) {
      if (connectivity[i][j]) {
        STRE *one_stre = new STRE(i,j);
        if (!present(one_stre)) {
          coords.simples.push_back(one_stre);
          ++nadded;
        }
        else
          delete one_stre;
      }
    }
  }
  return nadded;
}

// Add missing hydrogen bond stretches - return number added
// defined as [O,N,F,Cl]-H ... [O,N,F,Cl] with distance < 2.3 Angstroms
// and angle greater than 90 degrees
int FRAG::add_hbonds(void) {
  int nadded = 0;
  double dist, ang;
  const double pi = acos(-1);

  bool *is_X = init_bool_array(natom);
  bool *is_H = init_bool_array(natom);
  for (int i=0; i<natom; ++i) {
    if (Z[i] == 1)
      is_H[i] = true;
    else if (Z[i] == 7 || Z[i] == 8 || Z[i] == 9 || Z[i] == 17)
      is_X[i] = true;
  }

  double cov_scale = Opt_params.interfragment_scale_connectivity;
  double cov_H     = cov_radii[1]/_bohr2angstroms;

  for (int x=0; x<natom; ++x) {
    if (is_X[x]) { // electronegative atom
      for (int h=0; h<natom; ++h) {
        if (connectivity[x][h] && is_H[h]) { // x is bonded to h
          for (int y=0; y<natom; ++y) {
            if (y != x && is_X[y]) {    // eligible y exists
              dist = v3d_dist(geom[h], geom[y]);  // check distance
              if (dist < Opt_params.maximum_H_bond_distance) { // y is close enough to h
                if (v3d_angle(geom[x], geom[h], geom[y], ang)) { // now check bond angle
                  if (ang > pi/2) {
                    STRE *one_stre = new STRE(h,y);
                    int index = find(one_stre);
                    if (index == (int) coords.simples.size()) { // H-bond is absent
                      one_stre->set_hbond(true);
                      coords.simples.push_back(one_stre);
                      ++nadded;
                    }
                    else { // X-H ... Y stretch already exists,
                      // if it is not a covalent bond, make it a H bond
                      double cov_Y = cov_radii[ (int) Z[y] ]/_bohr2angstroms;
                      if (dist > cov_scale * (cov_H + cov_Y)) {
                        coords.simples[index]->set_hbond(true);
                      }
                      delete one_stre;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return nadded;
}

// Add auxiliary bonds; distance is < 2.5 times sum of covalent radii
int FRAG::add_auxiliary_bonds(void) {
  int nadded = 0;

  int *Zint = new int [natom];
  for (int i=0; i<natom; ++i)
    Zint[i] = (int) Z[i];

  for (int a=0; a<natom; ++a) {
    for (int b=a+1; b<natom; ++b) {
      if (connectivity[a][b]) continue; // already joined by regular bond

      // Omit auxiliary bonds involving H atoms
      if (Zint[a] == 1 || Zint[b] == 1) continue;

      double R = v3d_dist(geom[a], geom[b]);
      double Rcov = (cov_radii[Zint[a]] + cov_radii[Zint[b]])/_bohr2angstroms;

      if (R < Rcov * Opt_params.auxiliary_bond_factor) {

        bool omit = false;
        // Omit auxiliary bonds between a and b, if a-c-b
        for (int c=0; c<natom; ++c)
          if (c != a && c != b)
            if (connectivity[a][c] && connectivity[b][c])
              omit = true;

        // Omit auxiliary bonds between a and b, if a-c-d-b
        for (int c=0; c<natom; ++c)
          if (c != a && c != b)
            if (connectivity[c][a])
              for (int d=0; d<natom; ++d)
                if (d != a && d != b && d !=c)
                  if (connectivity[d][c] && connectivity[d][b])
                    omit = true;

        if (!omit) {
          STRE *one_stre = new STRE(a,b);
          if (!present(one_stre)) {
            coords.simples.push_back(one_stre);
            ++nadded;
          }
          else delete one_stre;
        }
        //else {
        //  oprintf_out("\tOmitting auxiliary bond %d %d bc of connectivity.\n", a+1, b+1);
        //}
      }
    }
  }
  delete [] Zint;
  return nadded;
}

// angles for all bonds present; return number added
int FRAG::add_bend_by_connectivity(void) {
  int nadded = 0;
  double phi;

  for (int i=0; i<natom; ++i)
    for (int j=0; j<natom; ++j)
      if (connectivity[i][j])
        for (int k=i+1; k<natom; ++k)
          if (connectivity[j][k]) {
            if (v3d_angle(geom[i], geom[j], geom[k], phi)) { // can be computed

              BEND *one_bend = new BEND(i,j,k);
              if (phi < Opt_params.linear_bend_threshold) { // < ~175 degrees
                if (!present(one_bend)) {
                  coords.simples.push_back(one_bend);
                  ++nadded;
                }
                else
                  delete one_bend;
              }
              else { // linear angle
                one_bend->make_lb_normal();
                if (!present(one_bend)) {
                  coords.simples.push_back(one_bend);
                  ++nadded;
                }
                else
                  delete one_bend;

                one_bend = new BEND(i,j,k);
                one_bend->make_lb_complement();
                if (!present(one_bend)) {
                  coords.simples.push_back(one_bend);
                  ++nadded;
                }
                else
                  delete one_bend;
              }
            }
          } // ijk

  for (ULI i=0; i<opt::INTCO_EXCEPT::linear_angles.size(); i+=3) {
    int A = opt::INTCO_EXCEPT::linear_angles[i];
    int B = opt::INTCO_EXCEPT::linear_angles[i+1];
    int C = opt::INTCO_EXCEPT::linear_angles[i+2];

    // Erase regular bend; add 2 linear bends
    BEND *old_bend = new BEND(A,B,C);
    if (present(old_bend)) {
      int i = find(old_bend);
      delete coords.simples[i];
      coords.simples.erase(coords.simples.begin() + i);
    }

    BEND *one_bend = new BEND(A,B,C);
    one_bend->make_lb_normal();
    if (!present(one_bend)) {
      oprintf_out("\tException forcing addition of linear bend (%d,%d,%d)\n", A+1, B+1, C+1);
      coords.simples.push_back(one_bend);
      ++nadded;
    }
    else
      delete one_bend;

    one_bend = new BEND(A,B,C);
    one_bend->make_lb_complement();
    if (!present(one_bend)) {
      oprintf_out("\tException forcing addition of linear bend complement (%d,%d,%d)\n", A+1, B+1, C+1);
      coords.simples.push_back(one_bend);
      ++nadded;
    }
    else
      delete one_bend;
  }
  opt::INTCO_EXCEPT::linear_angles.clear();

  return nadded;
}

// torsions for all bonds present; return number added
int FRAG::add_tors_by_connectivity(void) {
  int nadded = 0;
  int i,j,k,l;
  double phi;
  //double const phi_lim = Opt_params.linear_bend_threshold;
  // logic changed to use existence of linear bends to determine linearity
  // in this function.

  // bonding i-j-k-l but i-j-k && j-k-l are not collinear
  // use presence of linear bend coordinate to judge collinearity
  for (i=0; i<natom; ++i)
    for (j=0; j<natom; ++j)
      if (connectivity[i][j])
        for (k=0; k<natom; ++k)
          if ( connectivity[k][j] && k!=i ) {

            // ensure i-j-k is not collinear
            BEND *one_bend = new BEND(i,j,k);
            one_bend->make_lb_normal();
            if (present(one_bend)) {
              delete one_bend;
              continue;
            }
            delete one_bend;
            //if (!v3d_angle(geom[i], geom[j], geom[k], phi)) continue;
            //if (phi > phi_lim) continue;

            for (l=i+1; l<natom; ++l) {
              if (connectivity[l][k] && l!=j) {

                // ensure j-k-l is not collinear
                BEND *one_bend = new BEND(j,k,l);
                one_bend->make_lb_normal();
                if (present(one_bend)) {
                  delete one_bend;
                  continue;
                }
                delete one_bend;
                //if (!v3d_angle(geom[j], geom[k], geom[l], phi)) continue; // can't compute
                //if (phi > phi_lim) continue;

                TORS *one_tors = new TORS(i,j,k,l);
                if (!present(one_tors)) {
                  coords.simples.push_back(one_tors);
                  ++nadded;
                }
                else
                  delete one_tors;
              }
           }
         }

  // search for additional torsions around collinear segments
  int I,J,K,L,m;
  int nbonds;

  // find collinear fragment j-m-k
  for (j=0; j<natom; ++j)
    for (m=0; m<natom; ++m)
      if (connectivity[j][m])       // j!=m
        for (k=j+1; k<natom; ++k)   // k!=j
          if (connectivity[k][m]) { // m!=k

            BEND *one_bend = new BEND(j,m,k);
            one_bend->make_lb_normal();
            if (!present(one_bend)) {  // j-m-k are not collinear
              delete one_bend;
              continue;
            }
            delete one_bend;
            // Found unique, collinear j-m-k

            nbonds = 0; // count atoms bonded to m
            for (int n=0; n<natom; ++n)
              if (connectivity[n][m]) ++nbonds;

            if (nbonds == 2) { // nothing else is bonded to m

              // look for an 'I' for I-J-[m]-k-L such that I-J-K is not collinear
              J = j;
              for (i=0; i<natom; ++i) {
                if (connectivity[i][J] && i!=m) { // i!=J i!=m
                  BEND *one_bend = new BEND(i,j,k);
                  one_bend->make_lb_normal();
                  if (present(one_bend)) { // i,J,k is collinear
                    delete one_bend;
                    J = i;
                    i = 0;
                    continue;
                  }
                  else { // have I-J-[m]-k. Look for L
                    delete one_bend;
                    I = i;
                    K = k;
                    for (l=0; l<natom; ++l) {
                      if (connectivity[l][K] && l!=m && l!=j && l!=i) { // l!=K l!=m
                        BEND *one_bend = new BEND(l,K,J);
                        one_bend->make_lb_normal();
                        if (present(one_bend)) { // J-k-l is collinear
                          delete one_bend;
                          K = l;
                          continue;
                        }
                        else { // have found I-J-K-L
                          delete one_bend;
                          L = l;

oprintf_out("trying %d %d %d %d\n", I+1, J+1, K+1, L+1);
                          // ensure torsion is defined, i.e., angles not linear
                          if (!v3d_tors(geom[I], geom[J], geom[K], geom[L], phi))
                            continue;
oprintf_out("passed phi is %10.5lf\n", phi);

                          TORS *one_tors = new TORS(I,J,K,L);
                          if (!present(one_tors)) {
                            coords.simples.push_back(one_tors);
                            ++nadded;
                          }
                          else
                            delete one_tors;
                        } // end have IJKL
                      } // end l atom found
                    } // loop over l
                  } // end have IJK
                }
              } // look for i
            } // nbonds = 2 on central atom
          }
  return nadded;
}

// is simple already present in list ?
bool FRAG::present(const SIMPLE_COORDINATE *one) const {
  for (ULI k=0; k<coords.simples.size(); ++k) {
    if (*one == *(coords.simples[k]))
      return true;
  }
  return false;
}

void FRAG::add_combination_coord(vector<int> ids, vector<double> coeffs) {
  coords.index.push_back(ids);
  coords.coeff.push_back(coeffs);;
}

void FRAG::add_trivial_coord_combination(int simple_id) {
  std::vector<int> i1;
  i1.push_back(simple_id);
  coords.index.push_back(i1);

  std::vector<double> c1;
  c1.push_back(1.0);
  coords.coeff.push_back(c1);
}

int FRAG::form_trivial_coord_combinations(void) {
  coords.clear_combos();
  for (ULI s=0; s<coords.simples.size(); ++s)
    add_trivial_coord_combination(s);
  return coords.simples.size();
}

// Determine initial delocalized coordinate coefficients.
int FRAG::form_delocalized_coord_combinations(void) {
  // Get B matrix for simples
  int Nsimples = form_trivial_coord_combinations();
  double **B = compute_B();
  coords.clear_combos();

  oprintf_out("\n\tDiagonalizing (B B^t) to form delocalized coordinates for fragment.\n");
  oprintf_out("\tStarting with %d simple coordinates.\n", Nsimples);

  // Diagonalize B*B^t to get coordinate rows
  double **BBt = init_matrix(Nsimples, Nsimples);
  opt_matrix_mult(B, 0, B, 1, BBt, 0, Nsimples, 3*g_natom(), Nsimples, 0);
  free_matrix(B);

  double *evals = init_array(Nsimples);
  opt_symm_matrix_eig(BBt, Nsimples, evals);

  if (Opt_params.print_lvl > 2) {
    oprintf_out("Eigenvectors of BBt\n");
    oprint_matrix_out(BBt, Nsimples, Nsimples);
    oprintf_out("Eigenvalues of BBt\n");
    oprint_array_out(evals, Nsimples);
  }

  double eval_threshold  = 1e-8; // ??
  double evect_threshold = 1e-5; // ??

  for (int i=0; i<Nsimples; ++i) {
    if (fabs(evals[i]) < eval_threshold) {
      if (Opt_params.print_lvl > 2)
        oprintf_out("Eigenvector %d removed for low eigenvalue.\n",i+1);
    }
    else {
      // Delete tiny components to get cleaner functions.
      for (int j=0; j<Nsimples; ++j)
        if (fabs(BBt[i][j]) < evect_threshold)
          BBt[i][j] = 0.0;

      // Make largest component positive.
      double sign = array_max(BBt[i], Nsimples) / array_abs_max(BBt[i], Nsimples);
      if (sign < 0.99)
        array_scm(BBt[i], -1.0, Nsimples);

      // Normalize
      array_normalize(BBt[i], Nsimples);

      // Add to fragment coordinate set.
      vector<int> one_index;
      vector<double> one_coeff;

      for (int j=0; j<Nsimples; ++j) {
        if (fabs(BBt[i][j]) > 1.0e-14) {
          one_index.push_back(j);
          one_coeff.push_back(BBt[i][j]);
        }
      }
      coords.index.push_back(one_index);
      coords.coeff.push_back(one_coeff);
    }
  }
  free_matrix(BBt);
  free_array(evals);

  oprintf_out("\tInitially, formed %d delocalized coordinates for fragment.\n", coords.index.size());
  return coords.index.size();
}

int FRAG::add_cartesians(void) {
  int nadded = 0;

  for (int i=0; i<natom; ++i)
    for (int xyz=0; xyz<3;++xyz) {
      CART *one_cart = new CART(i,xyz);
      if (!present(one_cart)) {
        coords.simples.push_back(one_cart);
        ++nadded;
        vector<int> one_id;
        vector<double> one_coeff;
        one_id.push_back(coords.index.size());
        one_coeff.push_back(1.0);
        coords.index.push_back(one_id);
        coords.coeff.push_back(one_coeff);
      }
    }
  return nadded;
}

bool FRAG::is_noncart_present(void) const {
  for (ULI k=0; k<coords.simples.size(); ++k) {
    if (coords.simples[k]->g_type() != INTCO_TYPE::cart_type)
      return true;
  }

  return false;
}

// given any simple coordinate - this function looks to see if that coordinate
// is already present in the set.  If so, it returns the index.
// If not, it returns the index of the end + 1.
int FRAG::find(const SIMPLE_COORDINATE *one) const {
  for (ULI k=0; k<coords.simples.size(); ++k) {
    if (*one == *(coords.simples[k]))
      return k;
  }
  return coords.simples.size();
}

// Return values of all coordinates.
double * FRAG::coord_values(void) const {
  return coords.values(geom);
}
double * FRAG::coord_values(GeomType new_geom) const {
  return coords.values(new_geom);
}

// return the value of a single internal
double FRAG::coord_value(int coord_index) const {
  return coords.value(geom, coord_index);
}
double FRAG::coord_value(GeomType new_geom, int coord_index) const {
  return coords.value(new_geom, coord_index);
}

// Fills in B matrix of coordinates for this fragment.
// Uses passed-in memory.  Optional offsets allow molecule to build B matrix directly.
void FRAG::compute_B(double **B, int coord_offset, int atom_offset) const {

  for (int cc=0; cc<Ncoord(); ++cc)
    for (int x = 0; x <3*natom; ++x)
      B[coord_offset+cc][3*atom_offset + x] = 0.0;

  for (int cc=0; cc<Ncoord(); ++cc)
    coords.DqDx(geom, cc, B[coord_offset+cc], atom_offset);
}


// Returns B matrix of only the simple coordinates for this fragment.
/*
void FRAG::compute_B_simples(double **B, int coord_offset, int atom_offset) const {

  for (int s=0; s<coords.simples.size(); ++s) {
    double **dqdx_simple = coords.simples[s]->DqDx(geom);

    for (int j=0; j < coords.simples[s]->g_natom(); ++j) { // loop over atoms in s vector
      int atom = atom_offset + coords.simples[s]->g_atom(j);

      for (int xyz=0; xyz<3; ++xyz)
        B[s+coord_offset][3*atom + xyz] = dqdx_simple[j][xyz];
    }
  }
  return;
}
*/

// Returns B matrix of coordinates only for this fragment.
double ** FRAG::compute_B(void) const {
  double **newB = init_matrix(Ncoord(), 3*natom);
  compute_B(newB, 0, 0);
  return newB;
}

// Returns B' matrix for one internal coordinate computed with given geometry.
void FRAG::compute_derivative_B(GeomType new_geom, int coord_index, double **Bprime,
    int frag_atom_offset) const {
  coords.Dq2Dx2(new_geom, coord_index, Bprime, frag_atom_offset);
  return;
}

// Use present geometry
void FRAG::compute_derivative_B(int coord_index, double **Bprime, int frag_atom_offset) const {
  compute_derivative_B(geom, coord_index, Bprime, frag_atom_offset);
}

// Allocates and return B' matrix only for this fragment.
double ** FRAG::compute_derivative_B(int coord_index) const {
  double **Bprime =init_matrix(3*natom, 3*natom);
  compute_derivative_B(geom, coord_index, Bprime, 0);
  return Bprime;
}

// returns matrix of constraints (for now, 1's on diagonals to be frozen)
double ** FRAG::compute_constraints(void) const {
  double **C = init_matrix(coords.simples.size(), coords.simples.size());

  for (ULI i=0; i<coords.simples.size(); ++i)
    if (coords.simples[i]->is_frozen())
      C[i][i] = 1.0;

  return C;
}

// freeze coords within fragments
void FRAG::freeze_coords(void) {
  for (ULI i=0; i<coords.simples.size(); ++i)
    coords.simples[i]->freeze();
}

//// returns G matrix, mass-weighted or not
void FRAG::compute_G(double **G, bool use_masses) const {
  double **B = compute_B();

  if (use_masses) {
    for (int i=0; i<Ncoord(); ++i)
      for (int a=0; a<natom; ++a)
        for(int xyz=0; xyz<3; ++xyz)
          B[i][3*a+xyz] /= sqrt(mass[a]);
  }

  opt_matrix_mult(B, 0, B, 1, G, 0, Ncoord(), 3*natom, Ncoord(), 0);
  free_matrix(B);
  return;
}

void FRAG::fix_tors_near_180(void) {
  for (ULI i=0; i<coords.simples.size(); ++i)
    if (coords.simples[i]->g_type() == tors_type)
      coords.simples[i]->fix_tors_near_180(geom);
}

// Compute axes for bends, then mark as fixed.
void FRAG::fix_bend_axes(void) {
  BEND * a_bend;
  for (ULI i=0; i<coords.simples.size(); ++i)
    if (coords.simples[i]->g_type() == bend_type) {
      a_bend = static_cast<BEND*>(coords.simples[i]);
      // If value is small, then don't fix so that value and Bmatrix
      // element are correct for displacements through zero.
      if (a_bend->value(geom) > Opt_params.small_bend_fix_threshold) {
        a_bend->compute_axes(geom);
        a_bend->fix_axes();
      }
    }
}

void FRAG::unfix_bend_axes(void) {
  BEND * a_bend;
  for (ULI i=0; i<coords.simples.size(); ++i)
    if (coords.simples[i]->g_type() == bend_type) {
      a_bend = static_cast<BEND*>(coords.simples[i]);
      a_bend->unfix_axes();
    }
}

void FRAG::fix_oofp_near_180(void) {
  for (ULI i=0; i<coords.simples.size(); ++i)
    if (coords.simples[i]->g_type() == oofp_type)
      coords.simples[i]->fix_oofp_near_180(geom);
}

void FRAG::set_geom_array(double * geom_array_in) {
  int xyz, i, cnt = 0;
  for (i=0; i<natom; ++i)
    for (xyz=0; xyz<3; ++xyz)
      geom[i][xyz] = geom_array_in[cnt++];
}

void FRAG::set_geom(double ** geom_in) {
  for (int i=0; i<natom; ++i)
    for (int xyz=0; xyz<3; ++xyz)
      geom[i][xyz] = geom_in[i][xyz];
}

void FRAG::set_grad(double ** grad_in) {
  for (int i=0; i<natom; ++i) {
    grad[i][0] = grad_in[i][0];
    grad[i][1] = grad_in[i][1];
    grad[i][2] = grad_in[i][2];
  }
}

double ** FRAG::g_geom(void) const {
  double **g = matrix_return_copy(geom,natom,3);
  return g;
}

double * FRAG::g_geom_array(void) {
  int cnt=0;
  double * geom_array = init_array(3*natom);
  for (int i=0; i<natom; ++i)
    for (int xyz=0; xyz<3; ++xyz)
      geom_array[cnt++] = geom[i][xyz];
  return geom_array;
}

double * FRAG::g_grad_array(void) {
  int i, xyz, cnt=0;
  double * grad_array = init_array(3*natom);
  for (i=0; i<natom; ++i)
    for (xyz=0; xyz<3; ++xyz)
      grad_array[cnt++] = grad[i][xyz];
  return grad_array;
}

double ** FRAG::g_grad(void) {
  double **g = matrix_return_copy(grad,natom,3);
  return g;
}

double * FRAG::g_Z(void) const {
  double *z = init_array(natom);
  for (int i=0; i<natom; ++i) z[i] = Z[i];
  return z;
}

// Identify simple bond angles passing through 0 or going to 180.
//  returns a std::vector<int> containing quadruplets of frag #, A, B, C
// for atoms that should have linear coordinates
std::vector<int> FRAG::validate_angles(double const * const dq, int atom_offset) {

  // Compute change in simple coordinates.
  double *dq_simple = init_array(coords.simples.size());
  for (ULI cc=0; cc<coords.index.size(); ++cc)
    for (ULI s=0; s<coords.index[cc].size(); ++s)
      dq_simple[ coords.index[cc][s] ] += dq[cc] * coords.coeff[cc][s];

  std::vector<int> lin_angle;

  for (std::size_t s=0; s<coords.simples.size(); ++s)
    if (coords.simples[s]->g_type() == bend_type) {

      int A = coords.simples[s]->g_atom(0)+atom_offset;
      int B = coords.simples[s]->g_atom(1)+atom_offset;
      int C = coords.simples[s]->g_atom(2)+atom_offset;
      double new_bend_val = coords.simples[s]->value(geom) + dq_simple[s];

//oprintf_out("bend (%d,%d,%d) old: %10.5f new: %10.5f\n",A+1,B+1,C+1,coords.simples[s]->value(geom), new_bend_val);

      // angle is passing through 0.
      if (new_bend_val < 0.0) { // < ABC<0.  A-C-B should be linear bends.
        if (A < B) { lin_angle.push_back(A); lin_angle.push_back(C); lin_angle.push_back(B); }
        else       { lin_angle.push_back(B); lin_angle.push_back(C); lin_angle.push_back(A); }
      }

      // angle is passing through 180.0.  // < ABC~pi. Add A-B-C linear bends.
      if (new_bend_val > Opt_params.linear_bend_threshold) {
        BEND *one_bend = new BEND(A, B, C);
        one_bend->make_lb_normal();
        int loc = find(one_bend);
        if (((std::size_t) loc) != coords.simples.size()) // linear bend is already there.  No problem.
          continue;
        else {
          lin_angle.push_back(A); lin_angle.push_back(B); lin_angle.push_back(C);
        }
        delete one_bend;
      }

    }
  return lin_angle;
}

bool ** FRAG::g_connectivity(void) const {
  bool **c = init_bool_matrix(natom, natom);
  for (int i=0; i<natom; ++i)
    for (int j=0; j<=i; ++j)
      c[i][j] = c[j][i] = connectivity[i][j];
  return c;
}

const bool * const * FRAG::g_connectivity_pointer(void) const {
  return connectivity;
}

// compute COM
double * FRAG::com(GeomType in_geom) {
  double *center = init_array(3);
  double sum = 0.0;
  for (int i=0; i<g_natom(); ++i) {
    sum += mass[i];
    for (int xyz=0; xyz<3; ++xyz)
      center[xyz] += mass[i] * in_geom[i][xyz];
  }
  for (int xyz=0; xyz<3; ++xyz)
    center[xyz] /= sum;
  return center;
}

// compute intertia tensor
double ** FRAG::inertia_tensor (GeomType in_geom) {
  int xyz, xyz2;
  double tval;
  double *center = com(in_geom);
  double **I = init_matrix(3,3);

  for (int atom=0; atom<g_natom(); ++atom) {
    tval = 0.0;
    for (xyz=0; xyz<3; ++xyz) {
      for (xyz2=0; xyz2<3; ++xyz2) {
        if (xyz == xyz2)
          tval += (in_geom[atom][xyz] - center[xyz]) * (in_geom[atom][xyz] - center[xyz]);
        tval -= (in_geom[atom][xyz] - center[xyz]) * (in_geom[atom][xyz2] - center[xyz2]);
      }
      I[xyz][xyz2] = tval;
    }
  }
  free_array(center);
  return I;
}

// Compute principal axes.
int FRAG::principal_axes(GeomType in_geom, double **axes, double *evals) {

  double **I = inertia_tensor(in_geom);
  double *I_evals = init_array(3);

  opt_symm_matrix_eig(I, 3, I_evals);

  axes = init_matrix(3,3);
  evals = init_array(3);

  int cnt = 0;
  for (int i=0; i<3; ++i) {
    if (fabs(I_evals[i]) > 1.0e-14) {
      evals[cnt] = I_evals[i];
      axes[cnt][0] = I[i][0];
      axes[cnt][1] = I[i][1];
      axes[cnt][2] = I[i][2];
      ++cnt;
    }
  }
  free_array(I_evals);
  free_matrix(I);
  return cnt;
}

INTCO_TYPE FRAG::get_simple_type(int simple_index) {
  return coords.simples.at(simple_index)->g_type();
}

}

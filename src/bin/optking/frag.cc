/*!
   \file frag.cc
   \ingroup optking
   \brief fragment (molecule) class
*/

#include "frag.h"

#include <cmath>

#include "mem.h"
#include "v3d.h"
#include "atom_data.h"
#include "cov_radii.h"
#include "print.h"
#include "opt_data.h"
#include "physconst.h"
#include "linear_algebra.h"

#define EXTERN
#include "globals.h"

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
  for (int i=0; i<intcos.size(); ++i)
  delete intcos[i];
  intcos.clear();
}

// for now just set default masses - fix later
void FRAG::set_masses(void) {
  int i;
  for (i=0; i<natom; ++i)
    mass[i] = Z_to_mass[(int) Z[i]];
  return;
}

void FRAG::print_geom(FILE *fp, const int id, bool print_masses) {
  int i;
  fprintf(fp,"\t---Fragment %d Geometry---\n", id+1);
  if (print_masses) {
    for (i=0; i<natom; ++i)
      fprintf(fp,"\t %-4s%20.10lf%20.10lf%20.10lf%20.10lf\n",
        Z_to_symbol[(int) Z[i]], mass[i], geom[i][0], geom[i][1], geom[i][2]);
  }
  else {
    for (i=0; i<natom; ++i)
      fprintf(fp,"\t %-4s%20.10lf%20.10lf%20.10lf\n",
        Z_to_symbol[(int) Z[i]], geom[i][0], geom[i][1], geom[i][2]);
  }
  fprintf(fp, "\n");
  fflush(fp);
}

void FRAG::print_geom_grad(FILE *fp, const int id, bool print_masses) {
  int i;
  fprintf(fp,"\t---Fragment %d Geometry and Gradient---\n", id+1);
  if (print_masses) {
    for (i=0; i<natom; ++i)
      fprintf(fp,"\t %-4s%20.10lf%20.10lf%20.10lf%20.10lf\n",
        Z_to_symbol[(int) Z[i]], mass[i], geom[i][0], geom[i][1], geom[i][2]);
  }
  else {
    for (i=0; i<natom; ++i)
      fprintf(fp,"\t %-4s%20.10lf%20.10lf%20.10lf\n",
        Z_to_symbol[(int) Z[i]], geom[i][0], geom[i][1], geom[i][2]);
  }
  for (i=0; i<natom; ++i)
    fprintf(fp,"\t %24.10lf%20.10lf%20.10lf\n", grad[i][0], grad[i][1], grad[i][2]);
  fprintf(fp, "\n");
  fflush(fp);
}

void FRAG::print_geom(FILE *fp_geom) {
  for (int i=0; i<natom; ++i)
    fprintf(fp_geom, "\t  %3s  %15.10lf%15.10lf%15.10lf\n",
      Z_to_symbol[(int) Z[i]], geom[i][0], geom[i][1], geom[i][2]);
  fflush(fp_geom);
}


void FRAG::print_intcos(FILE *fp, int atom_offset) {
  fprintf(fp,"\t * Coordinate *           * BOHR/RAD *       * ANG/DEG *\n");
  for (int i=0; i<intcos.size(); ++i)
    intcos.at(i)->print(fp,geom,atom_offset);
  fprintf(fp, "\n");
  fflush(fp);
}

void FRAG::print_intco_dat(FILE *fp, int atom_offset) {
  //fprintf(fp,"\t---Fragment %d Intrafragment Coordinates---\n", id+1);
  for (int i=0; i<intcos.size(); ++i)
    intcos.at(i)->print_intco_dat(fp,atom_offset);
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

  for (int i=0; i<intcos.size(); ++i) {
    if (intcos.at(i)->g_type() == stre_type) {
      int a = intcos.at(i)->g_atom(0);
      int b = intcos.at(i)->g_atom(1);
      connectivity[a][b] = connectivity[b][a] = true;
    }
  }
}

void FRAG::print_connectivity(FILE *fp, const int id, const int offset) const {
  fprintf(fp,"\t---Fragment %d Bond Connectivity---\n", id+1);
  int i,j;
  for (i=0; i<natom; ++i) {
    fprintf(fp,"\t %d :", i+1+offset);
      for (j=0; j<natom; ++j)
        if (connectivity[i][j]) fprintf(fp," %d", j+1+offset);
    fprintf(fp,"\n");
  }
  fprintf(fp,"\n");
  fflush(fp);
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
          intcos.push_back(one_stre);
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
                    one_stre->make_hbond();
                    if (!present(one_stre)) {
                      intcos.push_back(one_stre);
                      ++nadded;
                    } 
                    else delete one_stre;
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

// angles for all bonds present; return number added
int FRAG::add_bend_by_connectivity(void) {
  int nadded = 0;
  double phi;

  for (int i=0; i<natom; ++i)
    for (int j=0; j<natom; ++j)
      if (connectivity[i][j])
        for (int k=i+1; k<natom; ++k)
          if (connectivity[j][k]) {
            BEND *one_bend = new BEND(i,j,k);
            if (!present(one_bend)) {
              intcos.push_back(one_bend);
              ++nadded;
            }
            else
              delete one_bend;

            // add linear bend complement if necessary
            if (v3d_angle(geom[i], geom[j], geom[k], phi)) { // can be computed
              if (phi > Opt_params.linear_bend_threshold) { // ~175 degrees
                one_bend = new BEND(i,j,k);
                one_bend->make_linear_bend();
                if (!present(one_bend)) {
                  intcos.push_back(one_bend);
                  ++nadded;
                }
                else
                  delete one_bend;
              }
            }
          } // ijk

  return nadded;
}

// torsions for all bonds present; return number added
int FRAG::add_tors_by_connectivity(void) {
  int nadded = 0;
  int i,j,k,l;
  double phi;
  double const pi = acos(-1);

  // bonding i-j-k-l but i-j-k && j-k-l are not collinear
  for (i=0; i<natom; ++i)
    for (j=0; j<natom; ++j)
      if (connectivity[i][j])
        for (k=0; k<natom; ++k)
          if ( connectivity[k][j] && k!=i ) {
            // ensure i-j-k is not collinear
            if ( !v3d_angle(geom[i], geom[j], geom[k], phi) ) continue;
            if (phi == pi) continue;
            for (l=i+1; l<natom; ++l) {
              if ( connectivity[l][k] && l!=j) {
                // ensure j-k-l is not collinear
                if ( !v3d_angle(geom[j], geom[k], geom[l], phi) ) continue; // can't compute
                if (phi == pi) continue;
                TORS *one_tors = new TORS(i,j,k,l);
                if (!present(one_tors)) {
                  intcos.push_back(one_tors);
                  ++nadded;
                }
                else
                  delete one_tors;
              }
           }
         }

  // search for additional torsions around collinear segments
  bool I_found, L_found, more_found;
  int I,J,K,L,m;
  int nbonds;

  // find collinear fragment j-m-k
  for (j=0; j<natom; ++j)
    for (m=0; m<natom; ++m)
      if (connectivity[j][m])
        for (k=j+1; k<natom; ++k)
          if (connectivity[k][m]) {
            if ( !v3d_angle(geom[j], geom[m], geom[k], phi) ) continue;
            if (phi == pi) { // found j-m-k collinear

              nbonds = 0;
              for (int n=0; n<natom; ++n)
                if (connectivity[n][m]) ++nbonds;
              if (nbonds == 2) { // nothing else is bonded to j

                // look for an 'I' for I-j-[m]-k-L such that I-J-K is not collinear
                J = j;
                for (i=0; i<natom; ++i) {
                  if (connectivity[i][J] && i!=m) {
                    if ( !v3d_angle(geom[i], geom[J], geom[k], phi) ) continue;
                    if (phi == pi) {
                      J = i;
                      i = 0;
                      continue;
                    }
                    else { // have I-J-K. Look for L
                      I = i;
                      K = k;
                      for (l=0; l<natom; ++l) {
                        if (connectivity[l][K] && l!=m) {
                          if ( !v3d_angle(geom[l], geom[K], geom[J], phi) ) continue;
                          if (phi == pi) {
                            K = l;
                            continue;
                          }
                          else { // have IJKL
                            L = l;

                            TORS *one_tors = new TORS(I,J,K,L);
                            if ( !v3d_tors(geom[I], geom[J], geom[K], geom[L], phi) ) continue;
                            if (!present(one_tors)) {
                              intcos.push_back(one_tors);
                              ++nadded;
                            }
                            else
                              delete one_tors;
                          } // end have IJKL
                        } // end l atom found
                      } // loop over l
                    } // end have IJK
                  }
                }
              }
            }
          }
  return nadded;
}

// is simple already present in list ?
bool FRAG::present(const SIMPLE *one) const {
  int k;
  for (k=0; k<intcos.size(); ++k) {
    if (*one == *(intcos[k]))
      return true;
  }
  return false;
}

// given any simple coordinate - this function looks to see if that coordinate
// is already present in the set.  If so, it returns the index.
// If not, it returns the index of the end + 1.
int FRAG::find(const SIMPLE *one) const {
  int k;
  for (k=0; k<intcos.size(); ++k) {
    if (*one == *(intcos[k]))
      return k;
  }
  return intcos.size();
}

// returns values of internal coordinates - using member geometry
double * FRAG::intco_values(void) const {
  double * q = init_array(intcos.size());
  for (int i=0; i<intcos.size(); ++i)
    q[i] = intcos.at(i)->value(geom);
  return q;
}

// returns values of internal coordinates - using given geometry
double * FRAG::intco_values(GeomType new_geom) const {
  double * q = init_array(intcos.size());
  for (int i=0; i<intcos.size(); ++i)
    q[i] = intcos.at(i)->value(new_geom);
  return q;
}

// returns B' matrix for one internal coordinate
double ** FRAG::compute_derivative_B(int intco_index) const {
  return compute_derivative_B(intco_index,geom);
}

// returns B' matrix for one internal coordinate computed with given geometry
double ** FRAG::compute_derivative_B(int intco_index, GeomType new_geom) const {
  double **frag_dq2dx2 = intcos.at(intco_index)->Dq2Dx2(new_geom);
  return frag_dq2dx2;
}

// returns B matrix of internal coordinates
double ** FRAG::compute_B(void) const {
  double **Bintco;
  double **B = init_matrix(intcos.size(), 3*natom);

  for (int i=0; i<intcos.size(); ++i) {
    Bintco = intcos.at(i)->DqDx(geom);

    for (int j=0; j < intcos.at(i)->g_natom(); ++j)
      for (int xyz=0; xyz<3; ++xyz)
        B[i][3*intcos.at(i)->g_atom(j) + xyz] = Bintco[j][xyz];

    free_matrix(Bintco);
  }
  return B;
}

// returns matrix of constraints (for now, 1's on diagonals to be frozen)
double ** FRAG::compute_constraints(void) const {
  double **C = init_matrix(intcos.size(), intcos.size());

  for (int i=0; i<intcos.size(); ++i)
    if (intcos[i]->is_frozen())
      C[i][i] = 1.0;

  return C;
}


// returns B matrix of internal coordinates; use previously allocated memory
void FRAG::compute_B(double **B) const {
  double **Bintco;

  zero_matrix(B, intcos.size(), 3*natom);  
  for (int i=0; i<intcos.size(); ++i) {
    Bintco = intcos.at(i)->DqDx(geom);

    for (int j=0; j < intcos.at(i)->g_natom(); ++j)
      for (int xyz=0; xyz<3; ++xyz)
        B[i][3*intcos.at(i)->g_atom(j) + xyz] = Bintco[j][xyz];

    free_matrix(Bintco);
  }
  return;
}

//// returns G matrix, mass-weighted or not
void FRAG::compute_G(double **G, bool use_masses) const {
  double **B = compute_B();

  if (use_masses) {
    for (int i=0; i<intcos.size(); ++i)
      for (int a=0; a<natom; ++a)
        for(int xyz=0; xyz<3; ++xyz)
          B[i][3*a+xyz] /= sqrt(mass[a]);
  }

  opt_matrix_mult(B, 0, B, 1, G, 0, intcos.size(), 3*natom, intcos.size(), 0);
  free_matrix(B);
  return;
}

// computes and print B matrix
void FRAG::print_B(FILE *fp) const {
  double **B = compute_B();
  fprintf(fp,"\t---B matrix---\n");
  print_matrix(fp, B, intcos.size(), 3*natom);
  fprintf(fp,"\n");
  fflush(fp);
  free_matrix(B);
}

void FRAG::fix_tors_near_180(void) {
  for (int i=0; i<intcos.size(); ++i) {
    if (intcos[i]->g_type() == tors_type)
      intcos[i]->fix_tors_near_180(geom);
  }
}

/*bool FRAG::check_tors_for_bad_angles(void) {
  for (int i=0; i<intcos.size(); ++i) {
    if (intcos[i]->g_type() == tors_type)
      intcos[i]->check_tors_for_bad_angles(geom);
  }
}*/

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

// don't let bends pass through zero - assumes dq is in right order for fragment
void FRAG::check_zero_angles(double const * const dq) {
  int i, cnt = 0;
  for (i=0; i<intcos.size(); ++i) {
    if (intcos[i]->g_type() == bend_type) {
      if (intcos[i]->value(geom) + dq[cnt++] < 0.0) 
        throw(INTCO_EXCEPT("Bond angle passing through zero", true));
    }
  }
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
  int atom, xyz, xyz2;
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

}


/*! \file frag.cc
    \ingroup OPT10
    \brief fragment (molecule) class
*/

#include "frag.h"

#include "mem.h"
#include "v3d.h"
#include "atom_data.h"
#include "cov_radii.h"
#include "print.h"
#include "opt_data.h"
#include "physconst.h"

#define EXTERN
#include "globals.h"

namespace opt {

using namespace v3d;

// Memory for Z and geom is provided by calling function
FRAG::FRAG(int natom_in, double *Z_in, double **geom_in) {
  natom = natom_in;
  Z = Z_in;
  geom = geom_in;

  grad = init_matrix(natom,3);
  connectivity = init_bool_matrix(natom,natom);
  mass = init_array(natom);
}

FRAG::~FRAG() {
  //printf("Destructing fragment\n");
  free_array(Z);
  free_matrix(geom);
  free_matrix(grad);
  free_array(mass);
  free_bool_matrix(connectivity);
  for (int i=0; i<intcos.size(); ++i)
  delete intcos[i];
  intcos.clear();
}

void FRAG::set_default_masses(void) {
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
      fprintf(fp,"\t %24.10lf%20.10lf%20.10lf\n",
        grad[i][0], grad[i][1], grad[i][2]);
  fprintf(fp, "\n");
}

void FRAG::write_geom(FILE *fp_geom) {
  for (int i=0; i<natom; ++i)
    fprintf(fp_geom, "( %15.10lf %15.10lf %15.10lf %15.10lf) \n",
      Z[i], geom[i][0], geom[i][1], geom[i][2]);
}


void FRAG::print_intcos(FILE *fp, const int id) {
  fprintf(fp,"\t---Fragment %d Intrafragment Coordinates---\n", id+1);
  fprintf(fp,"\t * Coordinate *           * BOHR/RAD *       * ANG/DEG *\n");
  //fprintf(fp,"\t ----------                 --------         ----------\n");
  //fprintf(fp,"\t Coordinate                 BOHR/RAD           ANG/DEG\n");
  //fprintf(fp,"\t ----------                 --------         ----------\n");
  int i;
  for (i=0; i<intcos.size(); ++i)
    intcos.at(i)->print(fp);
  fprintf(fp, "\n");
}

void FRAG::compute_intco_values(void) {
  int i;
  for (i=0; i<intcos.size(); ++i)
    intcos.at(i)->compute_val(geom);
}

// automatically determine bond connectivity by comparison of interatomic distance
//with scale_connectivity * sum of covalent radii
void FRAG::update_connectivity_by_distances(double scale) {
  int i, j, *Zint;
  double Rij;

  if (scale == -1) // not specified by argument, so use default
    scale = Opt_params.scale_connectivity;

  Zint = new int [natom];
  for (i=0; i<natom; ++i) {
    Zint[i] = (int) Z[i];
    if ( Zint[i] > LAST_COV_RADII_INDEX )
      throw("Warning: cannot automatically bond atom with strange atomic number");
  }

  for (i=0; i<natom; ++i) {
    for (j=0; j<i; ++j) {
      Rij = v3d_dist(geom[i], geom[j]);
      if (Rij < scale * (cov_radii[Zint[i]] + cov_radii[Zint[j]])/_bohr2angstroms)
        connectivity[i][j] = connectivity[j][i] = true;
    }
  }
  delete [] Zint;
}

void FRAG::print_connectivity(FILE *fp, const int id) const {
  fprintf(fp,"\t---Fragment %d Bond Connectivity---\n", id+1);
  int i,j;
  for (i=0; i<natom; ++i) {
    fprintf(fp,"\t %d :", i+1);
      for (j=0; j<natom; ++j)
        if (connectivity[i][j]) fprintf(fp," %d", j+1);
    fprintf(fp,"\n");
  }
  fprintf(fp,"\n");
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


// angles for all bonds present; return number added
int FRAG::add_bend_by_connectivity(void) {
  int nadded = 0;
  int i,j,k;

  for (i=0; i<natom; ++i)
    for (j=0; j<natom; ++j)
      if (connectivity[i][j])
        for (k=i+1; k<natom; ++k)
          if (connectivity[j][k]) {
            BEND *one_bend = new BEND(i,j,k);
            if (!present(one_bend)) {
              intcos.push_back(one_bend);
              ++nadded;
            }
            else
              delete one_bend;
          }

  return nadded;
}

// torsions for all bonds present; return number added
int FRAG::add_tors_by_connectivity(void) {
  int nadded = 0;
  int i,j,k,l;

  for (i=0; i<natom; ++i)
    for (j=0; j<natom; ++j)
      if (connectivity[i][j])
        for (k=0; k<natom; ++k)
          if ( connectivity[k][j] && k!=i )
            for (l=i+1; l<natom; ++l)
              if ( connectivity[l][k] && l!=j) {
                TORS *one_tors = new TORS(i,j,k,l);
                if (!present(one_tors)) {
                  intcos.push_back(one_tors);
                  ++nadded;
                }
                else
                  delete one_tors;
              }

  return nadded;
}

// is simple already present in list ?
bool FRAG::present(const SIMPLE *one) const {
  int k;
  for (k=0; k<intcos.size(); ++k) {
    if (intcos.at(k) == one)
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
    if (intcos.at(k) == one)
      return k;
  }
  return intcos.size();
}

// returns values of internal coordinates
double * FRAG::g_intco_values(void) const {
  double * q = init_array(intcos.size());
  for (int i=0; i<intcos.size(); ++i)
      q[i] = intcos.at(i)->g_val();
  return q;
}

// returns B matrix of internal coordinates
double ** FRAG::compute_B(void) const {
  int i, j, xyz;
  double *s1;
  double **B = init_matrix(intcos.size(), 3*natom);

  for (i=0; i<intcos.size(); ++i) {
    for (j=0; j < intcos.at(i)->g_natom(); ++j) {
      s1 = intcos.at(i)->g_s(geom, j);
      for (xyz=0; xyz<3; ++xyz)
        B[i][3*intcos.at(i)->g_atom(j) + xyz] = s1[xyz];
      free_array(s1);
    }
  }
  return B;
}

// returns B matrix of internal coordinates; use previously allocated memory
void FRAG::compute_B(double **B) const {
  int i, j, xyz;
  double *s1;

  for (i=0; i<intcos.size(); ++i) {
    for (j=0; j < intcos.at(i)->g_natom(); ++j) {
      s1 = intcos.at(i)->g_s(geom, j);
      for (xyz=0; xyz<3; ++xyz)
        B[i][3*intcos.at(i)->g_atom(j) + xyz] = s1[xyz];
      free_array(s1);
    }
  }
  return;
}


// computes and print B matrix
void FRAG::print_B(FILE *fp) const {
  double **B = compute_B();
  fprintf(fp,"\t---B matrix---\n");
  print_matrix(fp, B, intcos.size(), 3*natom);
  fprintf(fp,"\n");
}

void FRAG::fix_tors_near_180(void) {
  for (int i=0; i<intcos.size(); ++i) {
    if (intcos[i]->g_type() == tors_type)
      intcos[i]->fix_near_180();
  }
}

void FRAG::set_geom_array(double * geom_array_in) {
  int xyz, i, cnt = 0;
  for (i=0; i<natom; ++i)
    for (xyz=0; xyz<3; ++xyz)
      geom[i][xyz] = geom_array_in[cnt++];
}

void FRAG::set_grad(double ** grad_in) {
  for (int i=0; i<natom; ++i) {
    grad[i][0] = grad_in[i][0];
    grad[i][1] = grad_in[i][1];
    grad[i][2] = grad_in[i][2];
  }
} 

double * FRAG::g_geom_array(void) {
  int i, xyz, cnt=0;
  double * geom_array = init_array(3*natom);
  for (i=0; i<natom; ++i)
    for (xyz=0; xyz<3; ++xyz)
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

// don't let bends pass through zero - assumes dq is in right order for fragment
void FRAG::check_zero_angles(double const * const dq) {
  int i, cnt = 0;
  for (i=0; i<intcos.size(); ++i) {
    if (intcos[i]->g_type() == bend_type) {
      if (intcos[i]->g_val() + dq[cnt++] < 0.0) 
          throw("Bond angle passing through zero. Try new bonds, angles, or considering higher symmetry");
    }
  }
}

}


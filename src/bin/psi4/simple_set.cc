/*! \file simple.h
    \ingroup OPT09
    \brief simple internal coordinate base class
*/

#include "simple_set.h"
#include "cov_radii.h"

#include <iostream>

namespace psi { namespace opt09 {

// takes pointers to 2 simple coordinates - tells if they are identical
bool equiv(const SIMPLE * s1, const SIMPLE * s2) {
  if (s1->get_itype() != s2->get_itype())
    return false;
  else if (s1->get_itype() == stre)
    return (*dynamic_cast<const STRETCH *>(s1) == *dynamic_cast<const STRETCH *>(s2));
  else if (s1->get_itype() == bend)
    return (*dynamic_cast<const BEND *>(s1) == *dynamic_cast<const BEND *>(s2));
  else
    return false;
}

// given any simple coordinate - this function tells if that coordinate is already present
bool SIMPLE_SET::present(const SIMPLE * s2) const {
  for(int i=0; i<simples.size(); ++i) {
    if ( equiv(simples.at(i), s2) )
      return true;
  }
  return false;
}

// given any simple coordinate - this function looks to see if that coordinate
// is already present in the object.  If so, it returns the index.  If not, it returns -1.
int SIMPLE_SET::find_index(const SIMPLE * s2) const {
  for(int i=0; i<simples.size(); ++i) {
    if ( equiv(simples.at(i), s2) )
      return i;
  }
  return -1;
}

// tells whether an id number is already in use
bool SIMPLE_SET::id_present(int a) const {
  for(int i=0; i<simples.size(); ++i) {
    if (simples.at(i)->get_id() == a)
      return true;
  }
  return false;
}

// returns number of added coordinates
int SIMPLE_SET::add_simples_from_input(void) {
  int num=0,tot=0, i,j,a,b,c,d;
  char error[80];

  int size_initial = simples.size();

  try {
  ip_cwk_add(":INTCO");

  if (ip_exist("STRE",0)) {
    ip_count("STRE",&num,0);
    for(i=0;i<num;++i) {
      ip_count("STRE",&j,1,i);
      if (j != 3) {
        sprintf(error,"Stretch %d is of wrong dimension.",i+1);
        throw(error);
      }
      ip_data("STRE","%d",&(a),2,i,0);
      ip_data("STRE","%d",&(b),2,i,1);
      ip_data("STRE","%d",&(c),2,i,2);

      if (id_present(a)) {
        sprintf(error,"Simple internal id number %d is being used more than once.",a);
        throw(error);
      }

      STRETCH *one_stre = new STRETCH(a,b-1,c-1);
      if (!present(one_stre))
        simples.push_back(one_stre);
      else {
        fprintf(outfile,"Ignoring redundant stretch %d\n",a);
        delete one_stre;
      }
    }
  }

  if (ip_exist("BEND",0)) {
    ip_count("BEND",&num,0);
    for(i=0;i<num;++i) {
      ip_count("BEND",&j,1,i);
      if (j != 4) {
        sprintf(error,"Bend %d is of wrong dimension.",i+1);
        throw(error);
      }
      ip_data("BEND","%d",&(a),2,i,0);
      ip_data("BEND","%d",&(b),2,i,1);
      ip_data("BEND","%d",&(c),2,i,2);
      ip_data("BEND","%d",&(d),2,i,3);

      if (id_present(a)) {
        sprintf(error,"Simple internal id number %d is being used more than once.",a);
        throw(error);
      }

      BEND *one_bend = new BEND(a,b-1,c-1,d-1);
      if (!present(one_bend))
        simples.push_back(one_bend);
      else {
        fprintf(outfile,"Ignoring redundant bend %d\n",a);
        delete one_bend;
      }
    }
  }

  } // end try
  catch(char const *str) {
    printf("%s\n",str);
    printf("Trouble reading simple internal coordinates\n");
    fprintf(outfile,"%s\n",str);
    fprintf(outfile,"Trouble reading simple internal coordinates\n");
    abort();
  }

  return (simples.size() - size_initial);
}

// returns values of internal coordinates
double * SIMPLE_SET::get_q(void) const {
  double * q = new double [simples.size()];
  for (int i=0; i<simples.size(); ++i)
    q[i] = simples.at(i)->get_val();
  return q;
}

// returns B matrix of internal coordinates
double ** SIMPLE_SET::get_B(int natom) const {
  int i, j, xyz;
  double **B = block_matrix(simples.size(),3*natom);

  for (i=0; i<simples.size(); ++i)
    for (j=0; j<simples.at(i)->get_na(); ++j)
      for (xyz=0; xyz<3; ++xyz)
        B[i][3*simples.at(i)->get_atom(j)+xyz] = simples.at(i)->get_s(j,xyz) ;

  return B;
}

/* automatically generate simples from bond distances
natom = # of atoms;
Z = atomic numbers;
geom = geometry in au;
scale_connectivity = factor to multiply times cov_radii
*/
int SIMPLE_SET::add_simples_by_distance(int natom, double *Z,
  double *geom, double scale_connectivity) {

/*
  int i, j, Z1, Z2, already;
  double **R = block_matrix(natom,natom);
  for (i=0; i<natom; ++i)
    for (j=0; j<natom; ++j)
      R[i][j] = v3d_dist(geom[i],geom[j]);

  int **bonds = init_int_matrix(natom, natom);
  for (i=0; i<natom; ++i) {
    Z1 = (int) Z[i];
    for (j=0; j<i; ++j) {
      Z2 = (int) Z[j];
      if ( Z1>LAST_COV_RADII_INDEX || Z2>LAST_COV_RADII_INDEX)
        throw("Warning: cannot automatically bond atom with strange atomic number");
    if (R[i][j] < scale_connectivity * (cov_radii[Z1] + cov_radii[Z2]))
      bonds[i][j] = bonds[j][i] = 1;
    }
  }

  num = simples.size();
  for (i=0; i<natom; ++i) {
    for (j=i+1; j<natom; ++j) {
      if (bonds[i][j]) {
        STRETCH *one_stre = new STRETCH(num+1,i,j);

        if (find(stre.begin(),stre.end(),lstretch) == stre.end()) {
          stre.push_back(lstretch);
          ++num;
          ++num_new;
        }
        simples.push_back(one_stre);
      }
    }
  }
*/

}

}}

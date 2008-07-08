/*! \file
    \ingroup OPTKING
    \brief Class for stretches
*/

#include "intcos.h"

namespace psi { namespace optking {

extern void read_masses(int natom, double *zvals);

inline char * c_string(const string s);

Intcos::~Intcos(void) {
 printf("Destructing Intcos\n");
 stre.clear(); // calls Stretch destructors
 bend.clear();
 if (Bsimp_present) free_block(Bsimp);
}

Intcos::Intcos(const Intcos & s) {
 printf("using empty intcos copy constructor.\n");
}

int Intcos::size(void) const {
  return stre.size() + bend.size();
}

// default - read cartesian coordinates from chkpt file
Intcos::Intcos() {
 int i,xyz;

 chkpt_init(PSIO_OPEN_OLD);
 natom = chkpt_rd_natom();
 nallatom = chkpt_rd_nallatom();
 geom = chkpt_rd_geom();
 fgeom = chkpt_rd_fgeom();
 zvals = chkpt_rd_zvals();

 read_masses(natom,zvals); // move to input
 masses = chkpt_rd_masses();

 chkpt_close(); 

fprintf(outfile,"Masses:");
for (i=0;i<natom;++i)
  fprintf(outfile,"%15.10lf\n",masses[i]);

 for (i=0; i<natom; ++i)
   for (xyz=0; xyz<natom; ++xyz)
     geom[i][xyz] *= _bohr2angstroms;

 for (i=0; i<nallatom; ++i)
   for (xyz=0; xyz<nallatom; ++xyz)
     fgeom[i][xyz] *= _bohr2angstroms;

 Bsimp_present = 0;
}

// appends stretches from input, returns # found
int Intcos::add_intcos_from_input(void) {
  int nnew=0;
  // char *keyword = "BONDS";
nnew += add_stretches_from_input();
//  nnew += add_stretches_from_input("BONDS");
  nnew += add_bends_from_input();
  return nnew;
}

// return number of new stretches found
// default keyword is "STRE"
int Intcos::add_stretches_from_input(string key_in) {
  int i, nrow=0, j, a, b, num_new=0, max_atom=100; //fix later
  //char* keyword = &key_in[0];
  char *keyword = c_string(key_in);
  if (!ip_exist(keyword,0)) return 0;

  int num = stre.size();
  Stretch lstretch; // "l" is for local

  ip_count(keyword,&nrow,0);
  for(i=0; i<nrow;++i) {
    j=0;
    ip_count(keyword,&j,1,i);
    if (j != 2) {
      fprintf(outfile,"Warning: Stretch in row %d has wrong dimension.\n",i+1);
      continue;
    }
    a = b = -1;
    ip_data(keyword,"%d",&(a),2,i,0);
    ip_data(keyword,"%d",&(b),2,i,1);
    a-=1;
    b-=1;
    if ( a<0 || a>max_atom || b<0 || b>max_atom) {
      fprintf(outfile,"Warning: Stretch in row %d has bad atom number.\n",i+1);
      continue;
    }
printf("num_1,a,b: %d, %d, %d\n",num+1,a,b);
    lstretch.set(num+1,a,b);
    if (find(stre.begin(),stre.end(),lstretch) == stre.end()) {
      stre.push_back(lstretch);
      ++num;
      ++num_new;
    }
  }
  free(keyword);
  return num_new;
}

// return number of new bends found
// default keyword is "BEND"
int Intcos::add_bends_from_input(string key_in) {
  int i, nrow=0, j, a, b, c, num_new=0, max_atom=100; //fix later
  char* keyword = &key_in[0];
  if (!ip_exist(keyword,0)) return 0;

  int num = bend.size();
  fprintf(outfile,"Number of already bends: %d\n",num);
  Bend lbend; // local bend

  ip_count(keyword,&nrow,0);
  for(i=0; i<nrow;++i) {
    j=0;
    ip_count(keyword,&j,1,i);
    if (j != 3) {
      fprintf(outfile,"Warning: Bend in row %d has wrong dimension.\n",i+1);
      continue;
    }
    a = b = c = -1;
    ip_data(keyword,"%d",&(a),2,i,0);
    ip_data(keyword,"%d",&(b),2,i,1);
    ip_data(keyword,"%d",&(c),2,i,2);
    a-=1;
    b-=1;
    c-=1;
    if ( a<0 || a>max_atom || b<0 || b>max_atom || c<0 || c>max_atom) {
      fprintf(outfile,"Warning: Bend in row %d has bad atom number.\n",i+1);
      continue;
    }
    lbend.set(num+1,a,b,c);
    if (find(bend.begin(),bend.end(),lbend) == bend.end()) {
      bend.push_back(lbend);
      ++num;
      ++num_new;
    }
  }
  return num_new;
}


// return number of new stretches found
int Intcos::add_stretches_by_distance(double scale_connectivity) {
  int i, j, Z1, Z2, already;

  double **R;
  R = block_matrix(natom,natom);
  for (i=0; i<natom; ++i)
    for (j=0; j<natom; ++j)
      R[i][j] = v3d_dist(geom[i],geom[j]);

  // determine bond connectivity matrix
  int **bonds;
  bonds = init_int_matrix(natom, natom);
  for (i=0; i<natom; ++i) {
    Z1 = (int) zvals[i];
    for (j=0; j<i; ++j) {
      Z2 = (int) zvals[j];
      if (Z1==0 || Z2==0)  continue; // ghost
      if ( Z1>LAST_COV_RADII_INDEX || Z2>LAST_COV_RADII_INDEX) {
        fprintf(outfile,"Warning: cannot automatically bond atom %d.\n",i+1);
        continue;
      }
    if (R[i][j] < scale_connectivity * (cov_radii[Z1] + cov_radii[Z2]))
      bonds[i][j] = bonds[j][i] = 1;
    }
  }

  int num_new = 0;
  int num = stre.size();
  fprintf(outfile,"Number of already stretches: %d\n",num);
  Stretch lstretch; //local stretch
  for (i=0; i<natom; ++i) {
    for (j=i+1; j<natom; ++j) {
      if (bonds[i][j]) {
        lstretch.set(num+1,i,j);
        if (find(stre.begin(),stre.end(),lstretch) == stre.end()) {
          stre.push_back(lstretch);
          ++num;
          ++num_new;
        }
      }
    }
  }
  return num_new;
}

// adds angles between sequentially bonded atoms
int Intcos::add_angles_by_bonds(void) {
  int i,j,ia,ib,ja,jb,new_bend = 0;
  int num_bend = bend.size();
  Stretch lstretch;
  Bend lbend;
  for (i=0; i<stre.size(); ++i) {
    ia = stre.at(i).get_A();
    ib = stre.at(i).get_B();
    for (j=0; j<stre.size(); ++j) {
      ja = stre.at(j).get_A();
      jb = stre.at(j).get_B();
      if ((ia==ja) && (ib!=jb))
        lbend.set(num_bend+1,ib,ia,jb);
      else if ((ib==jb) && (ia!=ja))
        lbend.set(num_bend+1,ia,ib,ja);
      else if ((ib==ja) && (ia!=jb))
        lbend.set(num_bend+1,ia,ib,jb);
      else
        continue;
      if (find(bend.begin(),bend.end(),lbend) == bend.end()) {
        ++num_bend; ++new_bend;
        bend.push_back(lbend);
      }
    }
  }
  return new_bend;
}

void Intcos::print_Bsimp(void) const {
  if (size() && Bsimp_present)
    print_mat(Bsimp, size(), 3*natom, outfile);
}

// compute values
void Intcos::compute(void) {
  int i;
  for (i=0; i<stre.size(); ++i) stre.at(i).compute(geom);
  for (i=0; i<bend.size(); ++i) bend.at(i).compute(geom);
}

// compute values and s vectors
void Intcos::compute_s(void) {
  int i;
  compute();
  for (i=0; i<stre.size(); ++i) stre.at(i).compute_s(geom);
  for (i=0; i<bend.size(); ++i) bend.at(i).compute_s(geom);
}

// print intcos, print_flag if values desired
void Intcos::print(int print_flag) const {
  int i;
  for (i=0; i<stre.size(); ++i) stre.at(i).print(print_flag);
  for (i=0; i<bend.size(); ++i) bend.at(i).print(print_flag);
}

// print s vectors
void Intcos::print_s(void) const {
  int i;
  for (i=0; i<stre.size(); ++i) stre.at(i).print_s();
  for (i=0; i<bend.size(); ++i) bend.at(i).print_s();
}

//The following functions are generated by templates in tmpl.h
int Intcos::build_Bsimp(void) {
  if (Bsimp_present) free_block(Bsimp);
  Bsimp = block_matrix(size(),3*natom);
  Bsimp_present = 1;

  compute_s();
  int row=0;
  row = build_B_vtype(stre,row);
  row = build_B_vtype(bend,row);
printf("created Bsimp: size %d 3*natom %d\n",size(),3*natom);
}

inline char * c_string(const string s) {
  char * buf = new char [s.size()+1];
  for (int i=0; i<s.size(); ++i) buf[i] = s[i];
  buf[s.size()] = '\0';
  return buf;
}

}}


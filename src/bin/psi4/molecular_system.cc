#include <molecular_system.h>

namespace psi {

int Molecular_system::get_natoms(void) const {
  int i, num=0;
  for (i=0; i<fragment.size(); ++i)
    num += fragment.at(i).get_natoms();
  return num;
}

// copy constructor
Molecular_system::Molecular_system(const Molecular_system & sys) {
  fragment = sys.fragment;
  charge = sys.charge;
}

// default constructor: get molecule from input file
Molecular_system::Molecular_system(double conv_factor) {
  int nfragments = 1, i, errcod = 0;
  string geom_label;

  errcod = ip_data("NUM_FRAGMENTS","%d",&i,0);
  if (errcod == IPE_OK) nfragments = i;

  bool cart = false;
  if (ip_exist("GEOMETRY",0))
    cart = true;
  else if (ip_exist("ZMAT",0))
    throw("zmatrix not yet implemented");
  else
    throw("could not find GEOMETRY or ZMAT in input!");

  if (cart) {
    Fragment lfrag;
    ostringstream outstr;
    for (i=0; i<nfragments; ++i) {
      if (i==0)
        geom_label = "GEOMETRY";
      else {
        outstr << "GEOMETRY" << i;
        geom_label = outstr.str();
        outstr.clear();
      }
      lfrag.read_cartesian_from_input(geom_label,conv_factor);
      if (lfrag.read_masses_from_input(i)) {;}
      else (lfrag.read_default_masses());
      fragment.push_back(lfrag);
    }
  }
}

double **Molecular_system::get_geom(void) const {
  int j,xyz,cnt=0;
  double **lgeom = block_matrix(get_natoms(),3);
  vector<Fragment>::const_iterator it;

  for(it=fragment.begin(); it!=fragment.end(); ++it) {
    for (j=0; j<it->natoms; ++j) {
      for (xyz=0; xyz<3; ++xyz)
        lgeom[cnt][xyz] = it->geom[j][xyz];
      ++cnt;
    }
  }

  return lgeom;
}

double *Molecular_system::get_Z(void) const {
  int f, j, cnt=0;
  double *lZ = new double [get_natoms()];
  vector<Fragment>::const_iterator it;

  for(it=fragment.begin(); it!=fragment.end(); ++it)
    for (j=0; j<it->natoms; ++j)
      lZ[cnt++] = it->Z[j];

  return lZ;
}

string *Molecular_system::get_atom_label(void) const {
  int j,cnt=0;
  string *latom_label = new string [get_natoms()];
  vector<Fragment>::const_iterator it;

  for(it=fragment.begin(); it!=fragment.end(); ++it)
    for (j=0; j<it->natoms; ++j)
      latom_label[cnt++] = it->atom_label[j];

  return latom_label;
}

void Molecular_system::print(void) const {
  int j,xyz,i=0;

  vector<Fragment>::const_iterator it;
  for(it=fragment.begin(); it!=fragment.end(); ++it) {
    fprintf(outfile,"\n\tGeometry for fragment %d:\n",++i);
    it->print();
/*
    if (it->geom != NULL) {
      for (j=0; j<it->natoms; ++j) {
        fprintf(outfile,"\t\t%s",it->atom_label[j].c_str());
        for (xyz=0; xyz<3; ++xyz)
          fprintf(outfile,"%15.10lf",it->geom[j][xyz]);
      fprintf(outfile,"\n");
      }
    }
    else {
      fprintf(outfile,"Geometry is NULL.\n");
    }
*/
  }
  fprintf(outfile,"\n");
}

}

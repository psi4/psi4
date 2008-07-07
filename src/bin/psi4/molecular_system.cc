#include <Molecular_system.h>
#include <psi4-dec.h> // for outfile

namespace psi {

using std::vector;
using std::string;
using std::ostringstream;

// deep copy constructor
Molecular_system::Molecular_system(const Molecular_system & sys)
{
  fragment = sys.fragment;
  charge = sys.charge;
}

// default constructor: get molecule from input file
Molecular_system::Molecular_system(double conv_factor)
{
  int num_fragments=1, i, errcod = 0;
  string geom_label;

  // read how many fragments to try to read
  errcod = ip_data("NUM_FRAGMENTS","%d",&i,0);
  if (errcod == IPE_OK) num_fragments = i;

  // Cartesian or z-matrix input?
  bool cart = false;
  if (ip_exist("GEOMETRY",0))
    cart = true;
  else if (ip_exist("ZMAT",0))
    throw("zmatrix not yet implemented");
  else
    throw("could not find GEOMETRY or ZMAT in input!");

  // add fragment from input.dat
  if (cart) {
    Fragment lfrag;
    ostringstream outstr;
    for (i=0; i<num_fragments; ++i) {
      if (i==0)
        geom_label = "GEOMETRY";
      else {
        outstr << "GEOMETRY" << i;
        geom_label = outstr.str();
        outstr.clear();
      }
      lfrag.read_cartesian_from_input(geom_label,conv_factor);
      fragment.push_back(lfrag);
    }
  }
}

double **Molecular_system::get_geom(void)
{
  int j,xyz,cnt=0;
  double **lgeom = block_matrix(get_num_atoms(),3);
  vector<Fragment>::iterator it;

  for(it=fragment.begin(); it!=fragment.end(); ++it) {
    for (j=0; j<(*it).num_atoms; ++j) {
      for (xyz=0; xyz<3; ++xyz)
        lgeom[cnt][xyz] = (*it).geom[j][xyz];
      ++cnt;
    }
  }
  return lgeom;
}

double *Molecular_system::get_Z(void)
{
  int j,cnt=0;
  double *lZ = init_array(get_num_atoms());
  vector<Fragment>::iterator it;

  for(it=fragment.begin(); it!=fragment.end(); ++it)
    for (j=0; j<(*it).num_atoms; ++j)
      lZ[cnt++] = (*it).Z[j];

  return lZ;
}

string *Molecular_system::get_atom_label(void)
{
  int j,cnt=0;
  string *latom_label = new string [get_num_atoms()];
  vector<Fragment>::iterator it;

  for(it=fragment.begin(); it!=fragment.end(); ++it)
    for (j=0; j<(*it).num_atoms; ++j)
      latom_label[cnt++] = (*it).atom_label[j];

  return latom_label;
}

void Molecular_system::print(void) const
{
  int j,xyz,i=0;
  vector<Fragment>::const_iterator it;

  for(it=fragment.begin(); it!=fragment.end(); ++it) {
    fprintf(outfile,"\n\tGeometry for fragment %d:\n",++i);
    for (j=0; j<(*it).num_atoms; ++j) {
      for (xyz=0; xyz<3; ++xyz)
        fprintf(outfile,"\t\t%15.10lf",(*it).geom[j][xyz]);
      fprintf(outfile,"\n");
    }
  }
  fprintf(outfile,"\n");
}


} // psi


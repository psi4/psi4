#include <Molecular_system.h>
#include <psi4-dec.h>

namespace psi {

using std::vector;
using std::string;
using std::ostringstream;

// deep copy constructor
Molecular_system::Molecular_system(const Molecular_system & sys)
{
  num_atoms = sys.num_atoms;
  num_fragments = sys.num_fragments;
  fragment = sys.fragment;
  charge = sys.charge;
}

// default constructor: get molecule from input file
Molecular_system::Molecular_system()
{
  int i, errcod = 0;
  string geom_label;

  // read how many fragments there are
  num_fragments = 1;
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

  num_atoms = 0;

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
      lfrag.read_cartesian_from_input(geom_label);
      fragment.push_back(lfrag);
      // keep _system::num_atoms updated
      num_atoms += fragment[i].num_atoms;
    }
  }
}

double **Molecular_system::get_geom(void)
{
  int i,j,xyz,cnt=0;
  double **lgeom = block_matrix(num_atoms,3);

  for(i=0; i<num_fragments; ++i)
    for (j=0; j<fragment[i].get_num_atoms(); ++j)
      for (xyz=0; xyz<3; ++xyz)
        lgeom[cnt++][xyz] = fragment[i].geom[j][xyz];

  return lgeom;
}

double *Molecular_system::get_Z(void)
{
  int i,j,cnt=0;
  double *lZ = init_array(3*num_atoms);

  for(i=0; i<num_fragments; ++i)
    for (j=0; j<fragment[i].get_num_atoms(); ++j)
      lZ[cnt++] = fragment[i].Z[j];

  return lZ;
}

string *Molecular_system::get_atom_label(void)
{
  int i,j,cnt=0;

  string *latom_label = new string [num_atoms];
  for(i=0; i<num_fragments; ++i)
    for (j=0; j<fragment[i].num_atoms; ++j)
      latom_label[cnt++] = fragment[i].atom_label[j];

  return latom_label;
}

void Molecular_system::print(void) const
{
  int i,j,xyz;
  for(i=0; i<num_fragments; ++i) {
    fprintf(outfile,"Geometry for fragment %d:\n",i);
    for (j=0; j<fragment[i].num_atoms; ++j) {
      for (xyz=0; xyz<3; ++xyz)
        fprintf(outfile,"\t%20.10lf",fragment[i].geom[j][xyz]);
      fprintf(outfile,"\n");
    }
  }
}


} // psi


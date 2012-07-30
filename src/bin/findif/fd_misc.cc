/*! \file util.cc
    \ingroup OPTKING
    \brief miscellaneous
*/

#include "findif.h"

#include <physconst.h>

namespace psi { namespace findif {

bool ascending(const VIBRATION *vib1, const VIBRATION *vib2) {
  if (vib1->km < vib2->km)
    return true;
  else
    return false;
}

// function to print out (frequencies and normal modes) vector of vibrations
void print_vibrations(std::vector<VIBRATION *> modes) {

  const boost::shared_ptr<Molecule> mol = psi::Process::environment.molecule();
  char **irrep_lbls = mol->irrep_labels();
  int Natom = mol->natom();

  // compute harmonic frequencies, +/- in wavenumbers
  /* Convert evals from H/(kg bohr^2) to J/(kg m^2) = 1/s^2 */
  /* v = 1/(2 pi c) sqrt( eval ) */
  const double k_convert = _hartree2J/(_bohr2m * _bohr2m * _amu2kg);
  const double cm_convert = 1.0/(2.0 * _pi * _c * 100.0);

  for (int i=0; i<modes.size(); ++i) {
    if(modes[i]->km < 0.0)
      modes[i]->cm = -1*cm_convert * sqrt(-k_convert * modes[i]->km);
    else
      modes[i]->cm =    cm_convert * sqrt( k_convert * modes[i]->km);
  }

  // Sort modes by increasing eigenvalues.
  sort(modes.begin(), modes.end(), ascending);

  // Print out frequencies and irreps to output file.
  fprintf(outfile, "\n\t  Irrep      Harmonic Frequency   \n");
  fprintf(outfile,   "\t                  (cm-1)          \n");
  fprintf(outfile,   "\t-----------------------------------------------\n");

  for(int i=0; i<modes.size(); ++i) {
    if(modes[i]->cm < 0.0)
      fprintf(outfile, "\t  %5s   %15.4fi \n", irrep_lbls[modes[i]->irrep], -modes[i]->cm);
    else
      fprintf(outfile, "\t  %5s   %15.4f  \n", irrep_lbls[modes[i]->irrep], modes[i]->cm);
  }

  fprintf(outfile,   "\t-----------------------------------------------\n");
  fflush(outfile);

  // Return list of frequencies to wavefunction object.
  boost::shared_ptr<Vector> freq_vector(new Vector(modes.size()));
  for (int i=0; i<modes.size(); ++i)
    freq_vector->set(i, modes[i]->cm);

  //freq_vector->print_out();
  Process::environment.wavefunction()->set_frequencies(freq_vector);

  double sum = 0.0;
  for (int a=0; a<Natom; ++a)
     sum += mol->mass(a);

  // print out normal modes in format that WebMO likes
  fprintf(outfile, "\n\tNormal Modes (mass-weighted).\n");
  fprintf(outfile, "\tMolecular mass is %10.5f amu.\n", sum);
  fprintf(outfile, "\tFrequencies in cm^-1; force constants in au.\n");

  for(int i=0; i<modes.size(); ++i) { // print descending order
    if (fabs(cm_convert * sqrt(k_convert * fabs(modes[i]->km))) < 5.0) continue;
    fprintf(outfile,"\n");
    if (modes[i]->km < 0.0)
      fprintf(outfile, "   Frequency:      %8.2fi\n", cm_convert * sqrt(-k_convert * modes[i]->km));
    else
      fprintf(outfile, "   Frequency:      %8.2f\n", cm_convert * sqrt(k_convert * modes[i]->km));

    fprintf(outfile,   "   Force constant: %8.4f\n", modes[i]->km);

    //fprintf(outfile,   "   IR Intensity: %8.2f\n", irint[i]*ir_prefactor);

    fprintf(outfile, "\t     X       Y       Z           mass\t\n");
    for (int a=0; a<Natom; a++) {
      fprintf(outfile, "  %s \t", mol->symbol(a).c_str() );

      for (int xyz=0; xyz<3; ++xyz)
        fprintf(outfile, "%8.3f", modes[i]->lx[3*a+xyz]);

      fprintf(outfile,"%15.6f", mol->mass(a));

      fprintf(outfile, "\n");
    }
  }

  // awkward, but need nirrep to free labels
  int Nirrep = mol->point_group()->char_table().nirrep();

  for (int i=0; i<Nirrep; ++i)
    free(irrep_lbls[i]);
  free(irrep_lbls);
}

// displaces from a reference geometry: geom += salclist[salc_i] * disp_i * disp_size
// disp_size is in mass-weighted coordinates; cartesian displacment is DX/sqrt(mass)
void displace_cart(SharedMatrix geom, const CdSalcList & salclist,
  int salc_i, int disp_factor, double disp_size) {

  geom->set_name("Coord: " + to_string(salc_i) + ", Disp: " + to_string(disp_factor));

  const boost::shared_ptr<Molecule> mol = psi::Process::environment.molecule();

  int nc = salclist[salc_i].ncomponent();

  for (int c=0; c<nc; ++c) {
    int a          = salclist[salc_i].component(c).atom;
    int xyz        = salclist[salc_i].component(c).xyz;
    double coef    = salclist[salc_i].component(c).coef;

    geom->add(0, a, xyz, disp_factor * disp_size * coef / sqrt(mol->mass(a)));
  }

  return;
}

// displaces from a reference geometry.
// geom += salclist[salc_i] * disp_i * disp_size + salclist[salc_j] * disp_j * disp_size
// disp_size is in mass-weighted coordinates; cartesian displacment is DX/sqrt(mass)
void displace_cart(SharedMatrix geom, const CdSalcList & salclist,
  int salc_i, int salc_j, int disp_factor_i, int disp_factor_j, double disp_size) {

  geom->set_name("Coord: " + to_string(salc_i) + ", Disp: " + to_string(disp_factor_i)
    + "Coord: " + to_string(salc_j) + ", Disp: " + to_string(disp_factor_j));

  const boost::shared_ptr<Molecule> mol = psi::Process::environment.molecule();

  int a, xyz;
  double coef;

  for (int c=0; c<salclist[salc_i].ncomponent(); ++c) {
    a    = salclist[salc_i].component(c).atom;
    xyz  = salclist[salc_i].component(c).xyz;
    coef = salclist[salc_i].component(c).coef;

    geom->add(0, a, xyz, disp_factor_i * disp_size * coef / sqrt(mol->mass(a)));
  }

  for (int c=0; c<salclist[salc_j].ncomponent(); ++c) {
    a    = salclist[salc_j].component(c).atom;
    xyz  = salclist[salc_j].component(c).xyz;
    coef = salclist[salc_j].component(c).coef;

    geom->add(0, a, xyz, disp_factor_j * disp_size * coef / sqrt(mol->mass(a)));
  }

  return;
}

// it's assumed columns are cartesian dimensions
void mass_weight_columns_plus_one_half(SharedMatrix B) {
  const boost::shared_ptr<Molecule> mol = psi::Process::environment.molecule();
  double u;

  for (int col=0; col<B->ncol(); ++col) {
    u = sqrt(mol->mass(col/3));
    for (int row=0; row<B->nrow(); ++row)
      B->set(row, col, B->get(row,col) * u);
  }
}

}}


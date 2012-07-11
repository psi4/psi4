/*! \file thermo.cc
    \ingroup thermo
    \brief Compute thermodynamic quantities.
*/

#include <psi4-dec.h>
#include <libmints/mints.h>
#include <libmints/molecule.h>

#include <physconst.h>

/* thermo: Computes thermodynamic quantities.
 *  by Rollin King, 2012
*/

namespace psi { namespace thermo {
  void title(void);

PsiReturnType thermo(Options &options) {
  title();
  fflush(outfile);

  double T = options.get_double("T"); // T in K
  double P = options.get_double("P"); // P in Pascals

  // Read in essential data
  const boost::shared_ptr<Wavefunction> wf = psi::Process::environment.reference_wavefunction();

  const boost::shared_ptr<Molecule> mol = psi::Process::environment.molecule();

  // use this one?
  double E_elec = psi::Process::environment.globals["CURRENT ENERGY"];

  int Natom = mol->natom();
  int multiplicity = mol->multiplicity();
  const boost::shared_ptr<Vector> vib_freqs = wf->frequencies();

  Vector rot_const = mol->rotational_constants();
  RotorType rot_type = mol->rotor_type();

  mol->set_full_point_group();
  std::string pg = mol->full_point_group_with_n();
  std::string pg_n_replaced = mol->full_point_group();
  fprintf(outfile,"\tFull point group: %s (%s)\n", pg.c_str(), pg_n_replaced.c_str());
  int full_pg_n = mol->full_pg_n();

  int rot_symm_num;

  if (pg == "ATOM" || pg == "C1" || pg == "Ci" || pg == "Cs" || pg == "C_inf_v")
    rot_symm_num = 1;
  else if (pg == "D_inf_h")
    rot_symm_num = 2;
  else if (pg == "Td")
    rot_symm_num = 12;
  else if (pg == "Oh")
    rot_symm_num = 24;
  else if (pg == "Ih")
    rot_symm_num = 60;
  else if (pg == "Cn" || pg == "Cnv" || pg == "Cnh")
    rot_symm_num = full_pg_n;
  else if (pg == "Dn" || pg == "Dnd" || pg == "Dnh")
    rot_symm_num = 2*full_pg_n;
  else if (pg == "Sn")
    rot_symm_num = full_pg_n / 2;
  else 
    throw PsiException("thermo(): Could not interpret molecular point group.", __FILE__, __LINE__);

  // Set number of vibrational frequencies.
  int nvib_freqs;
  if (rot_type == 6) nvib_freqs = 0; //atom
  else if (rot_type == 3) nvib_freqs = 3*Natom-5; //linear
  else nvib_freqs = 3*Natom-6;

  if (vib_freqs->dim() != nvib_freqs) {
    fprintf(outfile, "\tERROR: Number of vibrational frequencies provided inconsistent with rotor type \
and number of atoms.\n");
    throw PsiException("thermo(): Wrong number of vibrational frequencies provided.", __FILE__, __LINE__);
  }

  fprintf(outfile,"\n\tData used to determine thermochemical information:\n");
  fprintf(outfile,  "\tTemperature (K): %15.2lf\n",T);
  fprintf(outfile,  "\tPressure (Pa)  : %15.2lf\n",P);
  fprintf(outfile,  "\tMultiplicity   : %15d\n",multiplicity);

  fprintf(outfile,  "\tRotor type     : %15s\n", RotorTypeList[rot_type].c_str());
  fprintf(outfile,  "\tRotational symmetry number : %3d\n",rot_symm_num);

  fprintf(outfile, "\n\tRotational constants:\n");
  fprintf(outfile, "\t\t   wavenumbers      GHz\n");
  if (rot_type < 4) {
    fprintf(outfile,"\t\tA:  %10.6lf   %10.5lf\n",rot_const[0],_c*rot_const[0]/1e7);
    fprintf(outfile,"\t\tB:  %10.6lf   %10.5lf\n",rot_const[1],_c*rot_const[1]/1e7);
    fprintf(outfile,"\t\tC:  %10.6lf   %10.5lf\n",rot_const[2],_c*rot_const[2]/1e7);
  }

  // Can eliminate when debugged - should be printed out in freq. code
  fprintf(outfile, "\n\tNuclear masses:\n");
  int cnt=0;
  for (int i=0; i<Natom; ++i) {
    if (cnt == 0) fprintf(outfile,"\t\t");
    fprintf(outfile,"%10.6f", mol->mass(i));
    ++cnt;
    if (cnt == 6 || i == (Natom-1)) {
      cnt = 0;
      fprintf(outfile,"\n");
    }
  }

  for (int i=0; i<nvib_freqs; ++i)
    if (vib_freqs->get(i) < 0) {
      fprintf(outfile, "\tWARNING: At least one vibrational frequency is imaginary!\n");
    }

  Vector vib_temp(nvib_freqs);
  Vector vib_energy(nvib_freqs);
  Vector q_vib(nvib_freqs);
  Vector s_vib(nvib_freqs);

  double ZPVE = 0, molecular_mass = 0;
  double  Etrans,  Eelec,  Evib, Erot,  Etotal;
  double Cvtrans, Cvelec, Cvvib, Cvrot, Cvtotal;
  double  qtrans,  qelec,  qvib, qrot,  qtotal;
  double  Strans,  Selec,  Svib, Srot,  Stotal;

  for(int i=0; i < Natom; i++)
    molecular_mass += mol->mass(i);

  const double kT = _kb * T;
  double phi_A, phi_B, phi_C;

  // variables Etrans, Cvtrans, and Strans are divided by R
  Etrans = 1.5 * T;
  Cvtrans = 1.5;
  qtrans = pow(_twopi * molecular_mass * _amu2kg * kT / (_h * _h), 1.5) * _na * kT / P ;
  Strans = 5.0/2.0 + log(qtrans/_na);

  // electronic part
  Eelec = 0;
  Cvelec = 0;
  qelec = multiplicity;
  Selec = log(qelec);

  // rotational part
  if(rot_type == 6) { // atom 
    Erot = Cvrot = Srot = 0;
  }
  else if(rot_type == 3) { // linear molecule
    Erot = T;
    Cvrot = 1.0;
    qrot = kT / (rot_symm_num * 100 * _c * _h * rot_const[1]); // B goes from cm^-1 to 1/s
    Srot = 1.0 + log(qrot);
  }
  else {
    Erot = 1.5 * T; 
    Cvrot = 1.5;
    phi_A = rot_const[0] * 100 * _h * _c / _kb;
    phi_B = rot_const[1] * 100 * _h * _c / _kb;
    phi_C = rot_const[2] * 100 * _h * _c / _kb;
    qrot = sqrt(_pi) * pow(T,1.5) / (rot_symm_num * sqrt(phi_A*phi_B*phi_C));
    Srot = 1.5 + log(qrot);
  }

  // vibrational part
  for(int i=0; i < nvib_freqs; i++)
    vib_temp[i] = 100 * _h * _c * vib_freqs->get(i) / _kb;

  if (nvib_freqs)
    fprintf(outfile,"\n\tVibrational frequencies (cm^-1)      Temperatures (K):\n");
  for (int i=0; i<nvib_freqs; ++i)
    fprintf(outfile, "\t%20.3f           %20.3f\n", vib_freqs->get(i), vib_temp[i]);

  for(int i=0; i < nvib_freqs; i++) {
    double rT = vib_temp[i] / T; // reduced T
    if (vib_temp[i] < 900)
      fprintf(outfile,"\tWarning: used thermodynamic relations are not appropriate for low frequency modes.");
    Evib += vib_temp[i] * (0.5 + 1.0 / (exp(rT) - 1));
    Svib += rT/(exp(rT) - 1) - log(1 - exp(-rT));
    Cvvib += exp(rT) * pow(rT/(exp(rT)-1), 2);
    // q_vib[i] = (exp(-vib_temp[i] / (2*T))) / (1 - exp(-vib_temp[i] / T));
    ZPVE += vib_freqs->get(i) / 2.0; //in cm^-2
  }

  // convert quantities in units of R into units of cal/mol
  double R_to_cal = _R / _cal2J;

  Eelec *= R_to_cal/1000; // go to kcal/mol
  Etrans *= R_to_cal/1000;
  Erot *= R_to_cal/1000;
  Evib *= R_to_cal/1000;
  Etotal = Eelec + Etrans + Erot + Evib;

  Selec *= R_to_cal;
  Strans *= R_to_cal;
  Srot *= R_to_cal;
  Svib *= R_to_cal;
  Stotal = Selec + Strans + Srot + Svib;

  Cvelec *= R_to_cal;
  Cvtrans *= R_to_cal;
  Cvrot *= R_to_cal;
  Cvvib *= R_to_cal;
  Cvtotal = Cvelec + Cvtrans + Cvrot + Cvvib;

  fprintf(outfile,"\n");
  fprintf(outfile,"\t                * Thermal Energy *       * Cv *          * S *\n");
  fprintf(outfile,"\t                        kcal/mol      cal/(mol K)    cal/(mol K) \n");
  fprintf(outfile,"\tElectronic    : %15.3lf%15.3lf%15.3lf\n", Eelec,  Cvelec,  Selec);
  fprintf(outfile,"\tTranslational : %15.3lf%15.3lf%15.3lf\n", Etrans, Cvtrans, Strans);
  fprintf(outfile,"\tRotational    : %15.3lf%15.3lf%15.3lf\n", Erot,   Cvrot,   Srot);
  fprintf(outfile,"\tVibrational   : %15.3lf%15.3lf%15.3lf\n", Evib,   Cvvib,   Svib);
  fprintf(outfile,"\tTotal         : %15.3lf%15.3lf%15.3lf\n", Etotal, Cvtotal, Stotal);

  double ZPVE_au = ZPVE * 100 * _h * _c / _hartree2J ; // cm^-1 -> au/particle

  fprintf(outfile,"\n\t                          cm^(-1)             au         kJ/mol\n");
  fprintf(outfile,"\tZero-point energy:    %15.10lf%15.10lf%15.10lf\n", ZPVE, ZPVE_au,
    ZPVE_au * _hartree2J / 1000 * _na);

  double DU = Etotal * 1000.0 * _cal2J / _na / _hartree2J ;
  double DH = DU + _kb * T / _hartree2J ;
  double DG = DH - T * Stotal * _cal2J / _na / _hartree2J ;

  fprintf(outfile,"\tEnergies in Hartree/particle:       Correction         Total\n");
  fprintf(outfile,"\t\tEnergy (0 K)         = %15.8lf  %15.8lf\n", ZPVE_au, E_elec+ZPVE_au);
  fprintf(outfile,"\t\tInternal energy      = %15.8lf  %15.8lf\n",  DU, E_elec + DU);
  fprintf(outfile,"\t\tEnthalpy             = %15.8lf  %15.8lf\n",  DH, E_elec + DH);
  fprintf(outfile,"\t\tGibbs Free Energy    = %15.8lf  %15.8lf\n",  DG, E_elec + DG);

  Process::environment.globals["GIBBS FREE ENERGY"] = E_elec + DG;

  return Success;
}

void title(void)
{
  fprintf(outfile, "\t\t\t*********************************************\n");
  fprintf(outfile, "\t\t\t* Thermodynamic Analysis by R.A. King, 2012 *\n");
  fprintf(outfile, "\t\t\t*********************************************\n");
}

}}


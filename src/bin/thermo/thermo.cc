/*! \file thermo.cc
    \ingroup thermo
    \brief Compute thermodynamic quantities.
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libqt/qt.h>
#include <libchkpt/chkpt.h>
#include <psifiles.h>
#include <physconst.h>
#include <masses.h>
#include <psi4-dec.h>

/* thermo: Computes thermodynamic quantities given:
 *
 *  rot_type : rotational type
 *  rotconst : rotational constants
 *  rot_symm_num : rotational symmetry number
 *  vib_freqs : vibrational frequencies
 *  by Taylor Mach and Rollin King, 2012
*/

namespace psi { namespace thermo {
  void title(void);

PsiReturnType thermo(Options &options) {

  title();

  // Read in thermochemical data
  double  E_elec 
  int     rot_symm_num 
  int     rot_type 
  double *rotconst = chkpt_rd_rotconst();
  double *vib_freqs = chkpt_rd_vib_freqs();
  double *zvals

  chkpt_init(PSIO_OPEN_OLD);
  int Natom = chkpt_rd_natom();
  double E_elec = chkpt_rd_etot(); // use psi::Process::environment.globals["CURRENT ENERGY"] ???
  int rot_symm_num = chkpt_rd_rot_symm_num();
  int rot_type = chkpt_rd_rottype();
  double *rotconst = chkpt_rd_rotconst();
  double *vib_freqs = chkpt_rd_vib_freqs();
  double *zvals = chkpt_rd_zvals(); // for default masses
  chkpt_close();

  double *masses;

  // Set default T, P
  double T = 298.15; // T in K
  double P = 101325; // P in Pascals
  int multiplicity = 1;

  // TODO: Add options for specifying T and P.

  // We could also have user options for rot_type, vib_freqs, such as
  /* options.add("ROTATIONAL_SYMMETRY_TYPE", "ASYMMETRIC_TOP",
       "ASYMMETRIC_TOP SYMMETRIC_TOP SPHERICAL_TOP LINEAR ATOM");
    "ASYMMETRIC_TOP" rot_type = 0;
    "SYMMETRIC_TOP"  rot_type = 1;
    "SPHERICAL_TOP"  rot_type = 2;
    "LINEAR"         rot_type = 3;
    "ATOM"           rot_type = 6;

      // asymmetric, symmetric and spherical tops
      if ( rot_type < 3 )
        expect 3 rotational constants
      else if (rot_type == 3) // linear molecule
        expect 1 rotational constant
      else if (rot_type == 6) // atom
        expect 0 rotational constants

    // Read vibrational frequencies
      Could read in vector of vibrational frequencies if provided.
*/

  // Set rotor string for output.
  std::string srotor; 
  if      (rot_type == 0) srotor = "ASYMMETRIC_TOP";
  else if (rot_type == 1) srotor = "SYMMETRIC_TOP";
  else if (rot_type == 2) srotor = "SPHERICAL_TOP";
  else if (rot_type == 3) srotor = "LINEAR";
  else if (rot_type == 6) srotor = "ATOM";

  // Set number of vibrational frequencies.
  int nvib_freqs;
  if (rot_type == 6) nvib_freqs = 0; //atom
  else if (rot_type == 3) nvib_freqs = 3*Natom-5; //linear
  else nvib_freqs = 3*Natom-6;

  // later, read this from checkpoint file
  masses = new double[Natom];
  for (i=0; i<Natom; ++i)
    masses[i] = an2masses[(int) zvals[i]]; 

  fprintf(outfile, "\tData used to determine thermochemical information:\n");
  fprintf(outfile, "\t\tRotor type: %s\n", srotor);
  fprintf(outfile, "\t\tRotational symmetry number: %d\n",rot_symm_num);
  fprintf(outfile, "\t\tRotational constants:\n");
  fprintf(outfile, "\t\t\t   wavenumbers,  GHz\n");
  if (rot_type < 4) {
    fprintf(outfile,"\t\t\tA: %10.6lf , %10.5lf\n",rotconst[0],_c*rotconst[0]/1e7);
    fprintf(outfile,"\t\t\tB: %10.6lf , %10.5lf\n",rotconst[1],_c*rotconst[1]/1e7);
    fprintf(outfile,"\t\t\tC: %10.6lf , %10.5lf\n",rotconst[2],_c*rotconst[2]/1e7);
  }
  if (nvib_freqs) fprintf(outfile,"\t\tVibrational frequencies:\n");
  for (i=0; i<nvib_freqs; ++i)
    fprintf(outfile, "\t\t\t%10.3f\n", vib_freqs[i]);
  fprintf(outfile, "\t\tTemperature (K): %10.2lf\n",T);
  fprintf(outfile, "\t\tPressure (Pa)  : %10.2lf\n",P);
  fprintf(outfile, "\t\tMultiplicity   : %10d\n",multiplicity);
  fprintf(outfile, "\t\tNuclear masses :\n");
  for (i=0; i<Natom; ++i)
    fprintf(outfile,"\t\t\t%10.6f\n", masses[i]);

  if (vib_freqs[i] < 0.0) {
    fprintf(outfile, "\tWARNING: A vibrational frequency is imaginary. The thermodynamic results are invalid!\n");
  }

  double *vib_temp   = init_array(nvib_freqs);
  double *vib_energy = init_array(nvib_freqs);
  double *q_vib      = init_array(nvib_freqs);
  double *s_vib      = init_array(nvib_freqs);

  double ZPVE, molecular_mass = 0;
  double  Etrans,  Eelec,  Evib, Erot,  Etotal;
  double Cvtrans, Cvelec, Cvvib, Cvrot, Cvtotal;
  double  qtrans,  qelec,  qvib, qrot,  qtotal;
  double  Strans,  Selec,  Svib, Srot,  Stotal;

  for(i=0; i < Natom; i++)
    molecular_mass += masses[i];

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
    qrot = kT / (rot_symm_num * 100 * _c * _h * rotconst[1]); // B goes from cm^-1 to 1/s
    Srot = 1.0 + log(qrot);
  }
  else {
    Erot = 1.5 * T; 
    Cvrot = 1.5;
    phi_A = rotconst[0] * 100 * _h * _c / _kb;
    phi_B = rotconst[1] * 100 * _h * _c / _kb;
    phi_C = rotconst[2] * 100 * _h * _c / _kb;
    qrot = sqrt(_pi) * pow(T,1.5) / (rot_symm_num * sqrt(phi_A*phi_B*phi_C));
    Srot = 1.5 + log(qrot);
  }

  // vibrational part
  for(i=0; i < nvib_freqs; i++)
    vib_temp[i] = 100 * _h * _c * vib_freqs[i] / _kb;

  fprintf(outfile,"\tVibrational temperatures (K):\n\t");
  for(i=0, cnt=0; i < nvib_freqs; i++,cnt++) {
    fprintf(outfile,"%12.3lf",vib_temp[i]);
    if (cnt == 6) { cnt=-1; fprintf(outfile,"\n\t"); }
  }
  fprintf(outfile,"\n");

  for(i=0; i < nvib_freqs; i++) {
    double rT = vib_temp[i] / T; // reduced T
    if (vib_temp[i] < 900)
      fprintf(outfile,"\tWarning: used thermodynamic relations are not appropriate for low frequency modes.");
    Evib += vib_temp[i] * (0.5 + 1.0 / (exp(rT) - 1));
    Svib += rT/(exp(rT) - 1) - log(1 - exp(-rT));
    Cvvib += exp(rT) * pow(rT/(exp(rT)-1), 2);
    // q_vib[i] = (exp(-vib_temp[i] / (2*T))) / (1 - exp(-vib_temp[i] / T));
    ZPVE += vib_freqs[i] / 2.0; //in cm^-2
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
  fprintf(outfile,"\tThermodynamic data\n");
  fprintf(outfile,"\t\t                Thermal Energy            Cv              S \n");
  fprintf(outfile,"\t\t                   kcal/mol           cal/(mol K)    cal/(mol K) \n");
  fprintf(outfile,"\t\tElectronic    : %15.3lf%15.3lf%15.3lf\n", Eelec,  Cvelec,  Selec);
  fprintf(outfile,"\t\tTranslational : %15.3lf%15.3lf%15.3lf\n", Etrans, Cvtrans, Strans);
  fprintf(outfile,"\t\tRotational    : %15.3lf%15.3lf%15.3lf\n", Erot,   Cvrot,   Srot);
  fprintf(outfile,"\t\tVibrational   : %15.3lf%15.3lf%15.3lf\n", Evib,   Cvvib,   Svib);
  fprintf(outfile,"\t\tTotal         : %15.3lf%15.3lf%15.3lf\n", Etotal, Cvtotal, Stotal);

  ZPVE *= 100 * _h * _c / _hartree2J ; // cm^-1 -> au/particle

  double U, H, G;
  U = E_elec + Etotal * 1000.0 * _cal2J / _na / _hartree2J ;
  H = U + _kb * T / _hartree2J ;
  G = H - T * Stotal * _cal2J / _na / _hartree2J ;

  fprintf(outfile,"\n\tTotal energies in Hartree/particle\n");
  fprintf(outfile,"\t\tTotal energy (0 K) = %15.7lf\n", E_elec + ZPVE);
  fprintf(outfile,"\t\tTotal energy       = %15.7lf\n", U);
  fprintf(outfile,"\t\tEnthalpy           = %15.7lf\n", H);
  fprintf(outfile,"\t\tFree Energy        = %15.7lf\n", G);

  free(vib_freqs);
  free(vib_temp);
  free(vib_energy);
  free(q_vib);
  free(s_vib);
  free(rotconst);
  free(zvals);

  tstop();

}

void title(void)
{
  fprintf(outfile, "\t\t\t*********************************\n");
  fprintf(outfile, "\t\t\t* Thermodynamic Analysis        *\n");
  fprintf(outfile, "\t\t\t* Taylor Mach & Rollin King '12 *\n");
  fprintf(outfile, "\t\t\t*********************************\n");
}

}}


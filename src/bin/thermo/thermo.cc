/*! \file thermo.cc
    \ingroup thermo
    \brief Compute thermodynamic quantities
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
 *  rottype : rotational type
 *  rotconst : rotational constants
 *  rot_symm_num : rotational symmetry number
 *  vib_freqs : vibrational frequencies
*/

namespace psi { namespace thermo {
  void init_io(int argc, char *argv[]);
  void exit_io(void);
  void title(void);

int thermo(int argc, char *argv[])
{
  using namespace psi::thermo;
  double tval;
  int errcod, a, i, cnt;
  char *junk;

  init_io(argc,argv);
  title();

  // Read in thermochemical data from checkpoint file
  chkpt_init(PSIO_OPEN_OLD);
  int rot_symm_num = chkpt_rd_rot_symm_num();
  int rottype = chkpt_rd_rottype();
  double *rotconst = chkpt_rd_rotconst();
  double *vib_freqs = chkpt_rd_vib_freqs();
  int natom = chkpt_rd_natom();
  double *zvals = chkpt_rd_zvals(); // for default masses
  double chkpt_energy = chkpt_rd_etot();
  chkpt_close();
  double *masses;

  // Set default T, P
  double T = 298.15; // T in K
  double P = 101325; // P in Pascals
  int multiplicity = 1;

  // Set number of vibrational frequencies consistent with chkpt file
  int nvib_freqs;
  if (rottype == 6) nvib_freqs = 0; //atom
  else if (rottype == 3) nvib_freqs = 3*natom-5; //linear
  else nvib_freqs = 3*natom-6;

  try { // Read overrides from input file
    // Read rotational symmetry type
    // integers are defined in chkpt_rd_rottype()
    if (ip_exist("ROT_TYPE",0)) {
      errcod = ip_string("ROT_TYPE", &(junk),0);
      if(errcod != IPE_OK) throw("Unable to read ROT_TYPE from input");

      if      (!strcmp(junk,"ASYMMETRIC_TOP")) rottype = 0;
      else if (!strcmp(junk,"SYMMETRIC_TOP")) rottype = 1;
      else if (!strcmp(junk,"SPHERICAL_TOP")) rottype = 2;
      else if (!strcmp(junk,"LINEAR")) rottype = 3;
      else if (!strcmp(junk,"ATOM")) rottype = 6;
      else throw("ROT_TYPE string was read but not recognized");

      free(junk);
    }

    // Update number of vibrations
    if (rottype == 6) nvib_freqs = 0; //atom
    else if (rottype == 3) nvib_freqs = 3*natom-5; //linear
    else nvib_freqs = 3*natom-6;

    // Read rotational symmetry number
    fprintf(outfile,"\tChecking input file for parameters.\n");
    if (ip_exist("ROT_SYMM_NUM",0)) {
      errcod = ip_data("ROT_SYMM_NUM", "%d", &a, 0);
      if (errcod != IPE_OK) throw("Unable to read ROT_SYMM_NUM");
      rot_symm_num = a;
    }

    // Read rotational constants
    if (ip_exist("ROT_CONST",0)) {
      errcod = ip_count("ROT_CONST",&a,0);
      if (errcod != IPE_OK) throw("Unable to count ROT_CONST entries");

      // asymmetric, symmetric and spherical tops
      if ( rottype < 3 ) {
        if (a != 3) throw("Expected 3 rotational constants with ROT_CONST in input");
        for (i=0; i<3; ++i) {
          errcod = ip_data("ROT_CONST", "%f", &(tval), 1, i);
          if (errcod != IPE_OK) throw("Unable to read rotational constant with ROT_CONST in input");
          rotconst[i] = tval;
        }
      }
      else if (rottype == 3) { // linear molecule
        if (a != 1) throw("Expected only 1 rotational constant with ROT_CONST in input");
        rotconst[0] = 0; rotconst[1] = rotconst[2] = tval;
      }
      else if (rottype == 6) { // atom
        throw("Don't specify any constants with ROT_CONST for an atom");
      }
    }

    // Read vibrational frequencies
    if (ip_exist("VIB_FREQS",0)) {
      errcod = ip_count("VIB_FREQS",&a,0);
      if (errcod != IPE_OK) throw("Unable to count VIB_FREQS entries");
      if (a != nvib_freqs) throw("Wrong number of entries in VIB_FREQS");
      for (i=0; i<nvib_freqs; ++i) {
        errcod = ip_data("VIB_FREQS", "%f", &(tval), 1, i);
        if (errcod != IPE_OK) throw("Unable to read vibrational frequency with VIB_FREQS in input");
        vib_freqs[i] = tval;
      }
    }

    // default multiplicity is 1 above
    a=0;
    if (ip_exist("MULTIPLICITY",0))
      errcod = ip_data("MULTIPLICITY", "%d", &a, 0);
    else if (ip_exist("MULTI",0))
      errcod = ip_data("MULTI", "%d", &a, 0);
    else if (ip_exist("MULTP",0))
      errcod = ip_data("MULTP", "%d", &a, 0);
    else if (ip_exist("SOCC",0)) {
      int num_ir, open, mguess;
      ip_count("SOCC",&(num_ir),0);
      for(i=0;i<num_ir;i++){
        errcod = ip_data("SOCC","%d",&open,1,i);
        mguess += open;
      }

      if(mguess == 1) a = 2;
      else if(mguess == 0) a = 1;
      else if(mguess == 2) throw("You must specify the MULTP keyword");
      else a = mguess + 1;
    }
    if (errcod != IPE_OK) throw("Unable to read multiplicity");
    if (a != 0) multiplicity = a;

    // Read in temperature and pressure
    if (ip_exist("T",0)) {
      errcod = ip_data("T", "%d",&(tval),0);
      if (errcod != IPE_OK) throw("Unable to read temperature in K");
      T = tval;
    }
    if (ip_exist("P",0)) {
      errcod = ip_data("P", "%d",&(tval),0);
      if (errcod != IPE_OK) throw("Unable to read pressure in Pa");
      P = tval;
    }
  }
  catch (const char *message) {
    printf("\n\t** Error: trouble reading parameters from input\n");
    fprintf(outfile,"\n\t** Error: trouble reading parameters from input\n");
    printf("\t\t%s",message);
    fprintf(outfile,"\t\t%s",message);
    exit_io();
    exit(PSI_RETURN_FAILURE);
  }

  // later, read this from checkpoint file
  masses = new double[natom];
  for (i=0; i<natom; ++i)
    masses[i] = an2masses[(int) zvals[i]]; 

  // set rotor string for output
  char srotor[80];   //string to store rotor type
  if     (rottype == 0) strcpy(srotor,"ASYMMETRIC_TOP");
  else if(rottype == 1) strcpy(srotor,"SYMMETRIC_TOP");
  else if(rottype == 2) strcpy(srotor,"SPHERICAL_TOP");
  else if(rottype == 3) strcpy(srotor,"LINEAR");
  else if(rottype == 6) strcpy(srotor,"ATOM");

  fprintf(outfile,"\tData used to determine thermochemical information:\n");
  fprintf(outfile,"\t\tRotor type: %s\n", srotor);
  fprintf(outfile,"\t\tRotational symmetry number: %d\n",rot_symm_num);
fprintf(outfile,"\t ** Check rotational symmetry number - not yet debugged!\n");
fprintf(outfile,"\t\t set rot_symm_num = 1 heteronuclear diatomics\n");
fprintf(outfile,"\t\t set rot_symm_num = 2 homonuclear diatomics\n");
  fprintf(outfile,"\t\tRotational constants:\n");
  fprintf(outfile,"\t\t\t   wavenumbers,  GHz\n");
  if (rottype < 4) {
    fprintf(outfile,"\t\t\tA: %10.6lf , %10.5lf\n",rotconst[0],_c*rotconst[0]/1e7);
    fprintf(outfile,"\t\t\tB: %10.6lf , %10.5lf\n",rotconst[1],_c*rotconst[1]/1e7);
    fprintf(outfile,"\t\t\tC: %10.6lf , %10.5lf\n",rotconst[2],_c*rotconst[2]/1e7);
  }
  if (nvib_freqs) fprintf(outfile,"\t\tVibrational frequencies:\n");
  for (i=0; i<nvib_freqs; ++i)
    fprintf(outfile,"\t\t\t%10.3f\n", vib_freqs[i]);
  fprintf(outfile,"\t\tTemperature (K): %10.2lf\n",T);
  fprintf(outfile,"\t\tPressure (Pa)  : %10.2lf\n",P);
  fprintf(outfile,"\t\tMultiplicity   : %10d\n",multiplicity);
  fprintf(outfile,"\t\tNuclear masses:\n");
  for (i=0; i<natom; ++i)
    fprintf(outfile,"\t\t\t%10.6f\n", masses[i]);

  for (i=0; i<nvib_freqs; ++i) {
    if (vib_freqs[i] < 0.0) {
      fprintf(outfile,"Imaginary vibrational frequencies are not allowed.\n");
      exit_io();
      exit(PSI_RETURN_FAILURE);
    }
  }

  double *vib_temp = init_array(nvib_freqs);
  double *vib_energy = init_array(nvib_freqs);
  double *q_vib = init_array(nvib_freqs);
  double *s_vib = init_array(nvib_freqs);

  double ZPVE = 0.0;
  double Etrans=0,  Eelec=0, Evib=0, Erot=0, Etotal;
  double Cvtrans=0, Cvelec=0, Cvvib=0, Cvrot=0, Cvtotal;
  double qtrans=0,  qelec=0, qvib=0, qrot=0, qtotal;
  double Strans=0,  Selec=0, Snuc=0, Svib=0, Srot=0, Stotal;
  double molecmass = 0.0;

  for(i=0; i < natom; i++)
    molecmass += masses[i];

  double kT = _kb * T;
  double rT;
  double phi_A, phi_B, phi_C;

  // variables Etrans, Cvtrans, and Strans are divided by R
  Etrans = 1.5 * T;
  Cvtrans = 1.5;
  qtrans = pow(_twopi * molecmass * _amu2kg * kT / (_h * _h), 1.5) * _na * kT / P ;
  Strans = 5.0/2.0 + log(qtrans/_na);

  // electronic part
  Eelec = 0;
  Cvelec = 0;
  qelec = multiplicity;
  Selec = log(qelec);

  // rotational part
  if(rottype == 6) { // atom 
    Erot = Cvrot = Srot = 0;
  }
  else if(rottype == 3) { // linear molecule
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
    rT = vib_temp[i] / T; // reduced T
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
  U = chkpt_energy + Etotal * 1000.0 * _cal2J / _na / _hartree2J ;
  H = U + _kb * T / _hartree2J ;
  G = H - T * Stotal * _cal2J / _na / _hartree2J ;

  fprintf(outfile,"\n\tTotal energies in Hartree/particle\n");
  fprintf(outfile,"\t\tTotal energy (0 K) = %15.7lf\n", chkpt_energy + ZPVE);
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


void init_io(int argc, char *argv[]) {
  int i, num_unparsed;
  char *progid, *argv_unparsed[100];

  progid = (char *) malloc(strlen(gprgid())+2);
  sprintf(progid, ":%s",gprgid());

  /* this code shows how to parse and extrace command-line arguments as needed */
/* X_only = 0;
  for (i=1, num_unparsed=0; i<argc; ++i) {
    if (!strcmp(argv[i],"--X_only"))
      X_only = 1;
    else
      argv_unparsed[num_unparsed++] = argv[i];
  }
*/
  for (i=1, num_unparsed=0; i<argc; ++i)
    argv_unparsed[num_unparsed++] = argv[i];

  ip_cwk_add(progid);
  free(progid);
  tstart();
}

void title(void)
{
  fprintf(outfile, "\t\t\t*********************************\n");
  fprintf(outfile, "\t\t\t*         THERMO                *\n");
  fprintf(outfile, "\t\t\t* Taylor Mach & Rollin King '09 *\n");
  fprintf(outfile, "\t\t\t*********************************\n");
}

void exit_io(void) {
  tstop();
  psi_stop(infile,outfile,psi_file_prefix);
}

}} // namespace psi::cphf


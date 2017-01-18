/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \file thermo.cc
    \ingroup thermo
    \brief Compute thermodynamic quantities.
*/

#include "psi4/psi4-dec.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/molecule.h"

#include "psi4/physconst.h"

/* thermo: Computes thermodynamic quantities.
 *  by Rollin King, 2012
*/

namespace psi { namespace thermo {
  void title(void);

PsiReturnType thermo(SharedWavefunction ref_wfn, SharedVector vib_freqs, Options &options) {
  title();


  double T = options.get_double("T"); // T in K
  double P = options.get_double("P"); // P in Pascals

  const std::shared_ptr<Molecule> mol = ref_wfn->molecule();

  double E_elec = psi::Process::environment.globals["CURRENT ENERGY"];

  int Natom = mol->natom();
  int multiplicity = mol->multiplicity();

  Vector rot_const = mol->rotational_constants();
  RotorType rot_type = mol->rotor_type();

  mol->set_full_point_group();
  std::string pg = mol->full_point_group_with_n();
  std::string pg_n_replaced = mol->full_point_group();
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

  if (options["ROTATIONAL_SYMMETRY_NUMBER"].has_changed())
    rot_symm_num = options.get_int("ROTATIONAL_SYMMETRY_NUMBER");

  // Set number of vibrational frequencies.
  int nvib_freqs;
  if (rot_type == RT_ATOM)
    nvib_freqs = 0;
  else if (rot_type == RT_LINEAR)
    nvib_freqs = 3 * Natom - 5;
  else
    nvib_freqs = 3 * Natom - 6;
  if (vib_freqs->dim() != nvib_freqs) {
    outfile->Printf("\n");
    outfile->Printf("    Not all frequencies have been computed, skipping thermodynamic analysis.\n");
    return Failure;
    //throw PsiException("thermo(): Wrong number of vibrational frequencies provided.", __FILE__, __LINE__);
    //outfile->Printf( "    ERROR: Number of vibrational frequencies provided inconsistent with rotor type and number of atoms.\n");
  }
  outfile->Printf("  ==> Inputs to Thermochemistry <==\n\n");

  outfile->Printf("  Temperature:              %15.2lf [K]\n",T);
  outfile->Printf("  Pressure:                 %15.2lf [Pa]\n",P);
  outfile->Printf("  Multiplicity:              %14d\n",multiplicity);
  outfile->Printf("  Full point group:                 %s (%s)\n", pg.c_str(), pg_n_replaced.c_str());
  outfile->Printf("  Rotor type:               %15s\n", RotorTypeList[rot_type].c_str());
  outfile->Printf("  Rotational symmetry no.:              %3d\n",rot_symm_num);
  if ((rot_type == RT_ASYMMETRIC_TOP) ||
      (rot_type == RT_SYMMETRIC_TOP) ||
      (rot_type == RT_SPHERICAL_TOP)) {
    mol->print_rotational_constants();
  }

  // Can eliminate when debugged - should be printed out in freq. code
  //outfile->Printf( "\n  Nuclear masses:\n");
  //int cnt=0;
  //for (int i=0; i<Natom; ++i) {
  //  if (cnt == 0) outfile->Printf("      ");
  //  outfile->Printf("%10.6f", mol->mass(i));
  //  ++cnt;
  //  if (cnt == 6 || i == (Natom-1)) {
  //    cnt = 0;
  //    outfile->Printf("\n");
  //  }
  //}

  //Flag to only print low/imaginary frequency warnings once
  bool WarnLowImag; //True if the warning was already printed

  WarnLowImag = 0; //Reset to false before checking imaginary frequencies
  for (int i=0; i<nvib_freqs; ++i)
    if ((vib_freqs->get(i) < 0) and (!WarnLowImag)) {
      outfile->Printf( "    WARNING: At least one vibrational frequency is imaginary!\n");
      WarnLowImag = 1; //Do not print a second message
    }

  Vector vib_temp(nvib_freqs);
  Vector vib_energy(nvib_freqs);
  Vector q_vib(nvib_freqs);
  Vector s_vib(nvib_freqs);

  double ZPVE = 0.0, molecular_mass = 0;
  double  Etrans = 0.0,  Eelec = 0.0,  Evib = 0.0, Erot = 0.0,  Etotal = 0.0;
  double Cvtrans = 0.0, Cvelec = 0.0, Cvvib = 0.0, Cvrot = 0.0, Cvtotal = 0.0;
  double Cptrans = 0.0, Cpelec = 0.0, Cpvib = 0.0, Cprot = 0.0, Cptotal = 0.0;
  double  qtrans = 0.0,  qelec = 0.0,  qvib = 0.0, qrot = 0.0,  qtotal = 0.0;
  double  Strans = 0.0,  Selec = 0.0,  Svib = 0.0, Srot = 0.0,  Stotal = 0.0;
  double  Htrans = 0.0,  Helec = 0.0,  Hvib = 0.0, Hrot = 0.0,  Htotal = 0.0;
  double  Gtrans = 0.0,  Gelec = 0.0,  Gvib = 0.0, Grot = 0.0,  Gtotal = 0.0;

  for(int i=0; i < Natom; i++)
    molecular_mass += mol->mass(i);

  const double kT = pc_kb * T;
  double phi_A = 0.0, phi_B = 0.0, phi_C = 0.0;

  // variables E*, Cv*, Cp*, and S* are divided by R

  // translation part
  qtrans = pow(pc_twopi * molecular_mass * pc_amu2kg * kT /
               (pc_h * pc_h), 1.5) * pc_na * kT / P;
  Strans = 2.5 + log(qtrans/pc_na);
  Cvtrans = 1.5;
  Cptrans = 2.5;
  Etrans = 1.5 * T;
  Htrans = 2.5 * T;

  Strans *= pc_R;  // [R] --> [J/mol.K]
  Cvtrans *= pc_R;  // [R] --> [J/mol.K]
  Cptrans *= pc_R;  // [R] --> [J/mol.K]
  Etrans *= 0.001 * pc_R; // [K.R] --> [kJ/mol]
  Htrans *= 0.001 * pc_R; // [K.R] --> [kJ/mol]
  Gtrans = Htrans - 0.001 * T * Strans;  // [kJ/mol]

  // electronic part
  qelec = multiplicity;
  Selec = log(qelec);
  Cvelec = 0;
  Cpelec = 0;
  Eelec = 0;
  Helec = 0;

  Selec *= pc_R;  // [R] --> [J/mol.K]
  Gelec = Helec - 0.001 * T * Selec;  // [kJ/mol]

  // rotational part
  if (rot_type == RT_ATOM) {
    Erot = Cvrot = Cprot = Srot = 0;
  }
  else if(rot_type == RT_LINEAR) {
    qrot = kT / (rot_symm_num * 100 * pc_c * pc_h * rot_const[1]); // B goes from cm^-1 to 1/s
    Srot = 1.0 + log(qrot);
    Cvrot = 1.0;
    Cprot = 1.0;
    Erot = T;
  }
  else {
    phi_A = rot_const[0] * 100 * pc_h * pc_c / pc_kb;
    phi_B = rot_const[1] * 100 * pc_h * pc_c / pc_kb;
    phi_C = rot_const[2] * 100 * pc_h * pc_c / pc_kb;
    qrot = sqrt(pc_pi) * pow(T,1.5) / (rot_symm_num * sqrt(phi_A*phi_B*phi_C));
    Srot = 1.5 + log(qrot);
    Cvrot = 1.5;
    Cprot = 1.5;
    Erot = 1.5 * T;
  }

  Srot *= pc_R;  // [R] --> [J/mol.K]
  Cvrot *= pc_R;  // [R] --> [J/mol.K]
  Cprot *= pc_R;  // [R] --> [J/mol.K]
  Erot *= 0.001 * pc_R;  // [K.R] --> [kJ/mol]
  Hrot = Erot;  // [kJ/mol]
  Grot = Hrot - 0.001 * T * Srot;  // [kJ/mol]

  // vibrational part
  for(int i=0; i < nvib_freqs; i++)
    vib_temp[i] = 100 * pc_h * pc_c * vib_freqs->get(i) / pc_kb;

  if (nvib_freqs)
    outfile->Printf("\n  No.    Vib. Freq. [cm^-1]      Vib. Temp. [K]\n");
  for (int i=0; i<nvib_freqs; ++i)
    outfile->Printf( "  %3i  %20.3f          %10.3f\n", i+1,vib_freqs->get(i), vib_temp[i]);

  WarnLowImag = 0; //Reset to false before checking low frequencies
  for(int i=0; i < nvib_freqs; i++) {
    double rT = vib_temp[i] / T; // reduced T
    if ((vib_temp[i] < 900) and (!WarnLowImag)) {
      outfile->Printf("    Warning: used thermodynamic relations are not appropriate for low frequency modes.\n");
      WarnLowImag = 1; //Do not print a second message
    }
    if (vib_temp[i] < 0) {
      outfile->Printf("    Warning: vibration with imaginary frequency neglected in vibrational contributions.\n");
      continue;
    }
    Evib += vib_temp[i] * (0.5 + 1.0 / (exp(rT) - 1));  // first term zpve
    Svib += rT/(exp(rT) - 1) - log(1 - exp(-rT));
    Cvvib += exp(rT) * pow(rT/(exp(rT)-1), 2);
    // q_vib[i] = (exp(-vib_temp[i] / (2*T))) / (1 - exp(-vib_temp[i] / T));
    ZPVE += vib_freqs->get(i) / 2.0;  // [cm^-1]
  }

  Svib *= pc_R;  // [R] --> [J/mol.K]
  Cvvib *= pc_R;  // [R] --> [J/mol.K]
  Cpvib = Cvvib;  // [J/mol.K]
  Evib *= 0.001 * pc_R;  // [K.R] --> [kJ/mol]
  Hvib = Evib;  // [kJ/mol]
  Gvib = Hvib - 0.001 * T * Svib;  // [kJ/mol]

  Stotal = Selec + Strans + Srot + Svib;  // [J/mol.K]
  Cvtotal = Cvelec + Cvtrans + Cvrot + Cvvib;  // [J/mol.K]
  Cptotal = Cpelec + Cptrans + Cprot + Cpvib;  // [J/mol.K]
  Etotal = Eelec + Etrans + Erot + Evib;  // [kJ/mol]
  Htotal = Helec + Htrans + Hrot + Hvib;  // [kJ/mol]
  Gtotal = Gelec + Gtrans + Grot + Gvib;  // [kJ/mol]

  outfile->Printf("\n  ==> Components <==\n\n");

  outfile->Printf("  Entropy, S\n");
  outfile->Printf("    Electronic S    %11.3lf [cal/(mol K)] %11.3lf [J/(mol K)] %15.8lf [mEh/K] (multiplicity = %d)\n",
    Selec / pc_cal2J, Selec, Selec / pc_hartree2kJmol, multiplicity);
  outfile->Printf("    Translational S %11.3lf [cal/(mol K)] %11.3lf [J/(mol K)] %15.8lf [mEh/K] (mol. weight = %.4f [u], P = %.2f [Pa])\n",
    Strans / pc_cal2J, Strans, Strans / pc_hartree2kJmol, molecular_mass, P);
  outfile->Printf("    Rotational S    %11.3lf [cal/(mol K)] %11.3lf [J/(mol K)] %15.8lf [mEh/K] (symmetry no. = %d)\n",
    Srot / pc_cal2J, Srot, Srot / pc_hartree2kJmol, rot_symm_num);
  outfile->Printf("    Vibrational S   %11.3lf [cal/(mol K)] %11.3lf [J/(mol K)] %15.8lf [mEh/K]\n",
    Svib / pc_cal2J, Svib, Svib / pc_hartree2kJmol);
  outfile->Printf("  Total S           %11.3lf [cal/(mol K)] %11.3lf [J/(mol K)] %15.8lf [mEh/K]\n\n",
    Stotal / pc_cal2J, Stotal, Stotal / pc_hartree2kJmol);

  outfile->Printf("  Constant volume heat capacity, Cv\n");
  outfile->Printf("    Electronic Cv   %11.3lf [cal/(mol K)] %11.3lf [J/(mol K)] %15.8lf [mEh/K]\n",
    Cvelec / pc_cal2J, Cvelec, Cvelec / pc_hartree2kJmol);
  outfile->Printf("    Translational Cv%11.3lf [cal/(mol K)] %11.3lf [J/(mol K)] %15.8lf [mEh/K]\n",
    Cvtrans / pc_cal2J, Cvtrans, Cvtrans / pc_hartree2kJmol);
  outfile->Printf("    Rotational Cv   %11.3lf [cal/(mol K)] %11.3lf [J/(mol K)] %15.8lf [mEh/K]\n",
    Cvrot / pc_cal2J, Cvrot, Cvrot / pc_hartree2kJmol);
  outfile->Printf("    Vibrational Cv  %11.3lf [cal/(mol K)] %11.3lf [J/(mol K)] %15.8lf [mEh/K]\n",
    Cvvib / pc_cal2J, Cvvib, Cvvib / pc_hartree2kJmol);
  outfile->Printf("  Total Cv          %11.3lf [cal/(mol K)] %11.3lf [J/(mol K)] %15.8lf [mEh/K]\n\n",
    Cvtotal / pc_cal2J, Cvtotal, Cvtotal / pc_hartree2kJmol);

  outfile->Printf("  Constant pressure heat capacity, Cp\n");
  outfile->Printf("    Electronic Cp   %11.3lf [cal/(mol K)] %11.3lf [J/(mol K)] %15.8lf [mEh/K]\n",
    Cpelec / pc_cal2J, Cpelec, Cpelec / pc_hartree2kJmol);
  outfile->Printf("    Translational Cp%11.3lf [cal/(mol K)] %11.3lf [J/(mol K)] %15.8lf [mEh/K]\n",
    Cptrans / pc_cal2J, Cptrans, Cptrans / pc_hartree2kJmol);
  outfile->Printf("    Rotational Cp   %11.3lf [cal/(mol K)] %11.3lf [J/(mol K)] %15.8lf [mEh/K]\n",
    Cprot / pc_cal2J, Cprot, Cprot / pc_hartree2kJmol);
  outfile->Printf("    Vibrational Cp  %11.3lf [cal/(mol K)] %11.3lf [J/(mol K)] %15.8lf [mEh/K]\n",
    Cpvib / pc_cal2J, Cpvib, Cpvib / pc_hartree2kJmol);
  outfile->Printf("  Total Cp          %11.3lf [cal/(mol K)] %11.3lf [J/(mol K)] %15.8lf [mEh/K]\n\n",
    Cptotal / pc_cal2J, Cptotal, Cptotal / pc_hartree2kJmol);

  outfile->Printf("  ==> Energy Analysis <== \n\n");

  double ZPVE_au = ZPVE / pc_hartree2wavenumbers;
  double DU = Etotal / pc_hartree2kJmol;
  double DH = Htotal / pc_hartree2kJmol;
  double DG = Gtotal / pc_hartree2kJmol;

  Process::environment.globals["ZPVE"] = ZPVE_au;
  Process::environment.globals["THERMAL ENERGY CORRECTION"] = DU;
  Process::environment.globals["ENTHALPY CORRECTION"] = DH;
  Process::environment.globals["GIBBS FREE ENERGY CORRECTION"] = DG;

  Process::environment.globals["ZERO K ENTHALPHY"] = E_elec + ZPVE_au;
  Process::environment.globals["THERMAL ENERGY"] = E_elec + DU;
  Process::environment.globals["ENTHALPY"] = E_elec + DH;
  Process::environment.globals["GIBBS FREE ENERGY"] = E_elec + DG;

  outfile->Printf("  Raw electronic energy, E0\n");
  outfile->Printf("  Total E0, Electronic energy at well bottom at 0 [K]           %15.8lf [Eh]\n\n",
    E_elec);

  outfile->Printf("  Zero-point energy, ZPE_vib = Sum_i nu_i / 2\n");
  outfile->Printf("    Electronic ZPE  %11.3lf [kcal/mol] %11.3lf [kJ/mol] %15.8lf [Eh]\n",
    0.0, 0.0, 0.0);
  outfile->Printf("    Translational ZPE %9.3lf [kcal/mol] %11.3lf [kJ/mol] %15.8lf [Eh]\n",
    0.0, 0.0, 0.0);
  outfile->Printf("    Rotational ZPE  %11.3lf [kcal/mol] %11.3lf [kJ/mol] %15.8lf [Eh]\n",
    0.0, 0.0, 0.0);
  outfile->Printf("    Vibrational ZPE %11.3lf [kcal/mol] %11.3lf [kJ/mol] %15.8lf [Eh] %15.3lf [cm^-1]\n",
    ZPVE_au * pc_hartree2kcalmol, ZPVE_au * pc_hartree2kJmol, ZPVE_au, ZPVE);
  outfile->Printf("    Correction ZPE  %11.3lf [kcal/mol] %11.3lf [kJ/mol] %15.8lf [Eh] %15.3lf [cm^-1]\n",
    ZPVE_au * pc_hartree2kcalmol, ZPVE_au * pc_hartree2kJmol, ZPVE_au, ZPVE);
  outfile->Printf("  Total ZPE, Electronic energy at 0 [K]                         %15.8lf [Eh]\n\n",
    E_elec + ZPVE_au);

  outfile->Printf("  Thermal Energy, E (includes ZPE)\n");
  outfile->Printf("    Electronic E    %11.3lf [kcal/mol] %11.3lf [kJ/mol] %15.8lf [Eh]\n",
    Eelec / pc_cal2J, Eelec, Eelec / pc_hartree2kJmol);
  outfile->Printf("    Translational E %11.3lf [kcal/mol] %11.3lf [kJ/mol] %15.8lf [Eh]\n",
    Etrans / pc_cal2J, Etrans, Etrans / pc_hartree2kJmol);
  outfile->Printf("    Rotational E    %11.3lf [kcal/mol] %11.3lf [kJ/mol] %15.8lf [Eh]\n",
    Erot / pc_cal2J, Erot, Erot / pc_hartree2kJmol);
  outfile->Printf("    Vibrational E   %11.3lf [kcal/mol] %11.3lf [kJ/mol] %15.8lf [Eh]\n",
    Evib / pc_cal2J, Evib, Evib / pc_hartree2kJmol);
  outfile->Printf("  Correction E      %11.3lf [kcal/mol] %11.3lf [kJ/mol] %15.8lf [Eh]\n",
    Etotal / pc_cal2J, Etotal, Etotal / pc_hartree2kJmol);
  outfile->Printf("  Total E, Electronic energy at %7.2f [K]                     %15.8lf [Eh]\n\n",
    T, E_elec + DU);

  outfile->Printf("  Enthalpy, H_trans = E_trans + k_B * T\n");
  outfile->Printf("    Electronic H    %11.3lf [kcal/mol] %11.3lf [kJ/mol] %15.8lf [Eh]\n",
    Helec / pc_cal2J, Helec, Helec / pc_hartree2kJmol);
  outfile->Printf("    Translational H %11.3lf [kcal/mol] %11.3lf [kJ/mol] %15.8lf [Eh]\n",
    Htrans / pc_cal2J, Htrans, Htrans / pc_hartree2kJmol);
  outfile->Printf("    Rotational H    %11.3lf [kcal/mol] %11.3lf [kJ/mol] %15.8lf [Eh]\n",
    Hrot / pc_cal2J, Hrot, Hrot / pc_hartree2kJmol);
  outfile->Printf("    Vibrational H   %11.3lf [kcal/mol] %11.3lf [kJ/mol] %15.8lf [Eh]\n",
    Hvib / pc_cal2J, Hvib, Hvib / pc_hartree2kJmol);
  outfile->Printf("  Correction H      %11.3lf [kcal/mol] %11.3lf [kJ/mol] %15.8lf [Eh]\n",
    Htotal / pc_cal2J, Htotal, Htotal / pc_hartree2kJmol);
  outfile->Printf("  Total H, Enthalpy at %7.2f [K]                              %15.8lf [Eh]\n\n",
    T, E_elec + DH);

  outfile->Printf("  Gibbs free energy, G = H - T * S\n");
  outfile->Printf("    Electronic G    %11.3lf [kcal/mol] %11.3lf [kJ/mol] %15.8lf [Eh]\n",
    Gelec / pc_cal2J, Gelec, Gelec / pc_hartree2kJmol);
  outfile->Printf("    Translational G %11.3lf [kcal/mol] %11.3lf [kJ/mol] %15.8lf [Eh]\n",
    Gtrans / pc_cal2J, Gtrans, Gtrans / pc_hartree2kJmol);
  outfile->Printf("    Rotational G    %11.3lf [kcal/mol] %11.3lf [kJ/mol] %15.8lf [Eh]\n",
    Grot / pc_cal2J, Grot, Grot / pc_hartree2kJmol);
  outfile->Printf("    Vibrational G   %11.3lf [kcal/mol] %11.3lf [kJ/mol] %15.8lf [Eh]\n",
    Gvib / pc_cal2J, Gvib, Gvib / pc_hartree2kJmol);
  outfile->Printf("  Correction G      %11.3lf [kcal/mol] %11.3lf [kJ/mol] %15.8lf [Eh]\n",
    Gtotal / pc_cal2J, Gtotal, Gtotal / pc_hartree2kJmol);
  outfile->Printf("  Total G, Free enthalpy at %7.2f [K]                         %15.8lf [Eh]\n\n",
    T, E_elec + DG);

  return Success;
}

void title(void)
{
  outfile->Printf("\n");
  outfile->Printf( "            *********************************************\n");
  outfile->Printf( "            * Thermodynamic Analysis by R.A. King, 2012 *\n");
  outfile->Printf( "            *********************************************\n");
  outfile->Printf("\n");
}

}}

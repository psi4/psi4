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

/*! \file
    \ingroup ccresponse
    \brief Compute the three tensors needed for Raman Optical Activity.

    ROA requires the following polarizability tensors:
      (1) electric-dipole/electric-dipole;
      (2) electric-dipole/electric-quadrupole; and
      (3) electric-dipole/magnetic-dipole.

  -TDC, August 2009
*/
//#include <stdio.h>
//#include <stdlib.h>
//#include <math.h>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libqt/qt.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"
#include <vector>
#include "psi4/psi4-dec.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/molecule.h"
#include "psi4/physconst.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libparallel/ParallelPrinter.h"

//#include "oldphysconst.h"
//#include "mass.h"


using namespace psi;

int levi(int a, int b, int c);
double tensor_mean(SharedMatrix alpha);
double beta_alpha2(SharedMatrix alpha);
double beta_G2(SharedMatrix alpha, SharedMatrix G);
double beta_A2(SharedMatrix alpha, double ***A, double omega);

double raman_linear(double alpha, double beta2);
double depolar_linear(double alpha, double beta2);
double raman_circular(double alpha, double beta2);
double depolar_circular(double alpha, double beta2);
void rs(int nm, int n, double **array, double *e_vals, int matz,
        double **e_vecs, double toler);


namespace psi { namespace ccresponse {

void print_tensor_der(std::shared_ptr<OutFile> myfile, std::vector<SharedMatrix> my_tensor_list);

//void scatter(double step, std::vector <SharedMatrix> pol, std::vector <SharedMatrix> rot, std::vector <SharedMatrix> quad)
void scatter(std::shared_ptr<Molecule> molecule, Options &options, double step, std::vector <SharedMatrix> pol, std::vector <SharedMatrix> rot, std::vector <SharedMatrix> quad)
{
    double mstep = options.get_double("DISP_SIZE");
    //-> This is troublesome, need to decide if option should be set as global or local-"FINDIF"
    outfile->Printf("STEPSIZE from Options Object = %lf bohr\n\n",mstep);
    outfile->Printf("STEPSIZE passed from roa.py and used = %lf bohr\n\n",step);

    int print = options.get_int("PRINT");
    outfile->Printf("Print Level = %d\n", print);
    int i,j,k;
    int a,b,c,d;

    // grab the field freqs from input -- a few units are converted to E_h
    //-> Do I want to handle an array of omega values,
    //-> OR do a separate scatter call for each one??
    //double omega = 0.085645;
    //double omega = params.omega[0];
    int count, nomega;
    double omega;
    count = options["OMEGA"].size();
    if(count == 0) { // Assume 0.0 E_h for field energy
      nomega = 1;
      params.omega = init_array(1);
      params.omega[0] = 0.0;
      omega = params.omega[0];
    }
    else if(count == 1) { // Assume E_h for field energy and read value
      params.nomega = 1;
      params.omega = init_array(1);
      params.omega[0] = options["OMEGA"][0].to_double();
      omega = params.omega[0];
    }
    else if(count == 2) {
      params.nomega = count-1;
      params.omega = init_array(params.nomega);
      std::string units = options["OMEGA"][count-1].to_string();
      for(i=0; i < (count-1); i++) {
        params.omega[i] = options["OMEGA"][i].to_double();
        if(units == "HZ" || units == "Hz" || units == "hz")
          params.omega[i] *= pc_h / pc_hartree2J;
        else if(units == "AU" || units == "Au" || units == "au") continue; // do nothing
        else if(units == "NM" || units == "nm")
          params.omega[i] = (pc_c*pc_h*1e9)/(params.omega[i]*pc_hartree2J);
        else if(units == "EV" || units == "ev" || units == "eV")
          params.omega[i] /= pc_hartree2ev;
        else
          throw PsiException("Error in unit for input field frequencies, should be au, Hz, nm, or eV", __FILE__,__LINE__);
      }
      omega = params.omega[0];
    }
    else if (count > 2) {
      throw PsiException("ROA Scattering only working for one wavelength at a time", __FILE__,__LINE__);
    }
    // Print the Wavelength
    outfile->Printf("\t Wavelength (in au): = %20.12f\n\n", omega);

    // COMPUTE TENSOR DERIVATIVES
    // Replicate the part of the roa.pl code that does this here in C++
    // Dipole Polarizability Tensor Gradients

    SharedMatrix denom_pol(new Matrix(3,3));
    denom_pol->set(2.0 * step);
    std::vector <SharedMatrix> pol_grad;
    for(std::vector<SharedMatrix>::iterator it_pol=pol.begin(); it_pol != pol.end(); ++it_pol) {
        SharedMatrix grad_mat(new Matrix(3,3));
        grad_mat->add(*it_pol);
        ++it_pol;
        grad_mat->subtract(*it_pol);
        grad_mat->apply_denominator(denom_pol);
        pol_grad.push_back(grad_mat);
    }

    SharedMatrix denom_rot(new Matrix(3,3));
    denom_rot->set(2.0 * step);
	/*
	 *  PSI4's OR Tensor is opposite in sign compared to PSI3.
	 *  If we tried the trick below to match PSI3 OR tensor signs,
     *  AND read in the normal coord. transform, we could match TDC's
	 *  numbers spot on, basically.
	 */
    //denom_rot->set(-2.0 * step);
    std::vector <SharedMatrix> rot_grad;
    for(std::vector<SharedMatrix>::iterator it_rot=rot.begin(); it_rot != rot.end(); ++it_rot) {
        SharedMatrix grad_mat(new Matrix(3,3));
        grad_mat->add(*it_rot);
        ++it_rot;
        grad_mat->subtract(*it_rot);
        grad_mat->apply_denominator(denom_rot);
        rot_grad.push_back(grad_mat);
    }

    SharedMatrix denom_quad(new Matrix(9,3));
    denom_quad->set(2.0 * step);
    std::vector <SharedMatrix> quad_grad;
    for(std::vector<SharedMatrix>::iterator it_quad=quad.begin(); it_quad != quad.end(); ++it_quad) {
        SharedMatrix grad_mat(new Matrix(9,3));
        grad_mat->add(*it_quad);
        ++it_quad;
        grad_mat->subtract(*it_quad);
        grad_mat->apply_denominator(denom_quad);
        quad_grad.push_back(grad_mat);
    }

    // Write Out the Tensor Derivatives to File tender.dat //
    // Outfile derivs("tender.dat", "w");
    std::shared_ptr<OutFile> derivs(new OutFile("tender.dat", TRUNCATE));
    derivs->Printf( "******************************************************\n");
    derivs->Printf( "**********                                  **********\n");
    derivs->Printf( "**********        TENSOR DERIVATIVES        **********\n");
    derivs->Printf( "**********   FOR COMPUTING ROA SCATTERING   **********\n");
    derivs->Printf( "**********                                  **********\n");
    derivs->Printf( "******************************************************\n\n\n");
    derivs->Printf( "\t*** Dipole Polarizability Derivative Tensors ***\n\n");
    print_tensor_der(derivs, pol_grad);
    derivs->Printf( "*********************************************************\n\n");
    derivs->Printf( "\t*** Optical Rotation Derivative Tensors ***\n\n");
    print_tensor_der(derivs, rot_grad);
    derivs->Printf( "*********************************************************\n\n");
    derivs->Printf( "\t*** Dipole/Quadrupole Derivative Tensors ***\n\n");
    print_tensor_der(derivs, quad_grad);

    int natom = molecule->natom();
    SharedMatrix geom(new Matrix(natom,3));

    // Reading in the Hessian //
    FILE* hessian;
    FILE* dipole_moment;
    hessian=fopen("file15.dat","r");
    SharedMatrix F(new Matrix(natom*3,natom*3));
    double Fval;
    SharedMatrix M(new Matrix(natom*3,natom*3));
    for(i=0; i < (3*natom); i++)
    {
      for(j=0; j < (3*natom); j++)
      {
        int statusvalue=fscanf(hessian,"%lf",&F->pointer()[i][j]);
      }
    }
    fclose(hessian);

	// Read in the Dipole-Moment Derivatives //
    dipole_moment=fopen("file17.dat","r");
    SharedMatrix dipder(new Matrix(3, natom*3));
    for(i=0; i < 3; i++)
    {
      for(j=0; j < natom; j++)
      {
        int statusvalue=fscanf(dipole_moment, "%lf %lf %lf", &dipder->pointer()[i][j*3], &dipder->pointer()[i][j*3+1], &dipder->pointer()[i][j*3+2]);
      }
    }
    fclose(dipole_moment);

	// Convert Vectors of SharedMatrices to Single SuperMatrix //
    SharedMatrix polder(new Matrix(natom*3, 9));
    for (i=0;i<3*natom;i++)
    {
      for (j=0;j<9;j++)
      {
        polder->set(i,j,pol_grad[i]->get(j/3,j%3));
      }
    }

    SharedMatrix optder(new Matrix(natom*3, 9));
    for (i=0;i<3*natom;i++)
    {
      for (j=0;j<9;j++)
      {
        optder->set(i,j,rot_grad[i]->get(j/3,j%3));
      }
    }

    SharedMatrix quadder(new Matrix(natom*3, 27));
    for (i=0;i<3*natom;i++)
    {
      for (j=0;j<27;j++)
      {
        quadder->set(i,j,quad_grad[i]->get(j/3,j%3));
      }
    }

	// Probably delete all these?
    //quadder->print(stdout);
    //for (i=0;i<pol_grad.size();i++)
    //{
      //pol_grad[i]->print(stdout);
    //}
    //polder->print(stdout);
    //F->print(stdout);

    geom->copy(molecule->geometry());
    SharedMatrix geom_orig(new Matrix(natom,3));
    geom_orig->copy(geom);

	// Translate Molecule to Center of Mass //
	if(print >= 1)  {
      outfile->Printf("\tInput coordinates (bohr):\n");
	  molecule->geometry().print();
    }
	molecule->move_to_com();
	geom_orig->copy(molecule->geometry());

    // Reading the Z-vals //
    //-> Just use the Molecule's mass(x) function to get appropriate mass
    //   of atom.
/*
    int zvals[natom];
    for(i=0; i<natom; i++)
    {
        zvals[i]=molecule->Z(i);
    }
*/

    // Mass-weighting the co-ordinates //
    //double massi[natom];
    outfile->Printf("\tAtomic Masses for Raman Computation:\n");
    for(i=0; i < natom; i++) outfile->Printf("\t%d %12.8f\n", i, molecule->mass(i));
    outfile->Printf("\n");

    for(i=0; i<natom; i++)
    {
        //massi[i] = molecule->mass(i);
        for(j=0; j<3; j++)
        {
            //geom->set(i,j, geom_orig->get(i,j) * sqrt(massi[i]));
            geom->set(i,j, geom_orig->get(i,j) * sqrt(molecule->mass(i)));
        }
    }
	if(print >= 1)  {
	  outfile->Printf("\tMass-Weighted coordinates relative to the center of mass:\n");
	  geom->print();
    }
    //delete[] massi;

    //Generating the inertia tensor

    SharedMatrix I(new Matrix(3,3));
    I->copy(molecule->inertia_tensor());
    //I = molecule->inertia_tensor();
	if(print >= 2)  {
	  outfile->Printf("\tMoment of Inertia Tensor:\n");
      I->print();
    }

/*
    for(i=0; i < 3; i++)
    {
      for(j=0; j <= i; j++)
      {
        if(i==j)
          for(k=0; k < natom; k++)
            I->add(i,j,(geom->get(k,(i+1)%3)*geom->get(k,(i+1)%3) + geom->get(k,(i+2)%3)*geom->get(k,(i+2)%3)));
        else
        {
          for(k=0; k < natom; k++)
            I->add(i,j,-1.0 * (geom->get(k,i)*geom->get(k,j)));
        I->set(j,i,I->get(i,j));
        }
      }
    }
    //I->print();
*/

    // Diagonalizing the inertia tensor //

    SharedMatrix Ievecs(new Matrix("Inertia Eigenvectors",3,3));
    SharedVector Ievals(new Vector("Inertia Eigenvalues",3));

    I->diagonalize(Ievecs,Ievals);
	//rs(3,3,I->pointer(),Ievals->pointer(),1,Ievecs->pointer(),1e-12);

    // Constructing I-inverse matrix //

    SharedMatrix Iinv(new Matrix("I-Inverse",3,3));
    SharedMatrix Itmp(new Matrix(3,3));

    Iinv->zero();
    for(i=0;i<3;i++)
    {
      Iinv->set(i,i,(1.0/Ievals->get(i)));
    }

    Itmp->gemm(0,1,1.0,Iinv,Ievecs,0.0);
    Iinv->gemm(0,0,1.0,Ievecs,Itmp,0.0);

	if(print >= 2)  {
	  outfile->Printf("\tInertia Tensor EigenData\n\n");
	  Ievecs->print();
	  Ievals->print();
      Iinv->print();
    }

    // Generating the 6 pure rotation and translation vectors //
    SharedMatrix P(new Matrix(natom*3,natom*3));
    int icart,jcart,iatom,jatom;
    double imass,jmass,total_mass;

    total_mass=0.0;
    for(i=0;i<natom;i++)
    {
      //total_mass+=an2mass[(int)zvals[i]];
      total_mass+=molecule->mass(i);
    }

    for(i=0; i < natom*3; i++)
    {
      icart = i % 3;
      iatom = i/3;
      //imass = an2mass[(int)zvals[iatom]];
      imass = molecule->mass(iatom);

      P->set(i,i,1.0);

      for(j=0; j < natom*3; j++)
      {
        jcart = j % 3;
        jatom = j/3;
        //jmass = an2mass[(int)zvals[jatom]];
        jmass = molecule->mass(jatom);

        P->add(i,j,-1.0*sqrt(imass*jmass)/total_mass*(icart==jcart));

        for(a=0; a < 3; a++)
          for(b=0; b < 3; b++)
            for(c=0; c < 3; c++)
              for(d=0; d < 3; d++)
              {
                P->add(i,j,-1.0*levi(a,b,icart)*geom->get(iatom,b)*Iinv->get(a,c)*levi(c,d,jcart)*geom->get(jatom,d));
              }

      }
    }

    // Generate mass-weighted Hessian matrix [Eh/(bohr^2 amu)] //
    M->zero();
    SharedMatrix T(new Matrix(natom*3,natom*3));
    for(i=0; i < natom; i++)
    {
      for(j=0; j < 3; j++)
      {
        //M->set((i*3+j),(i*3+j),1/sqrt(an2mass[(int)zvals[i]]/_au2amu));
        M->set((i*3+j),(i*3+j),1/sqrt((molecule->mass(i))/pc_au2amu));
      }
    }
    //printf("Mass-Weighting Matrix (for Hessian):\n");
	//M->print(stdout);

    T->gemm(0,0,1.0,M,F,0.0);
    F->gemm(0,0,1.0,T,M,0.0);

    // Project out rotational and translational degrees of freedom from mass-weighted Hessian //
    T->zero();
    T->gemm(0,0,1.0,F,P,0.0);
    F->gemm(0,0,1.0,P,T,0.0);
	if(print >= 2)  {
	  outfile->Printf("\tProjected, Mass-Weighted Hessian:\n");
	  F->print();
    }

    // Diagonalize projected mass-weighted Hessian //
    SharedMatrix Fevecs(new Matrix(3*natom,3*natom));
    SharedVector Fevals(new Vector("Feigenval",3*natom));
    SharedMatrix Lx(new Matrix("Normal Transform Matrix",3*natom,3*natom));
    SharedVector redmass(new Vector("ReducedMass",3*natom));
    double norm=0.0;

    F->diagonalize(Fevecs,Fevals);
	//rs(3*natom,3*natom,F->pointer(),Fevals->pointer(),1,Fevecs->pointer(),1e-12);
	if(print >= 3)  {
	  Fevals->print();
	  Fevecs->print();
    }

    Lx->gemm(0,0,1.0,M,Fevecs,0.0);
/*
	FILE* lxf = fopen("lx.dat","r");
    for(i=0; i < (3*natom); i++)
    {
      for(j=0; j < (3*natom); j++)
      {
        fscanf(lxf,"%lf",&Lx->pointer()[i][j]);
      }
    }
    fclose(lxf);
*/
	if(print >= 2)  {
	  Lx->print();
    }

    for(i=0; i < 3*natom; i++)
    {
       norm = 0.0;
       for(j=0; j < 3*natom; j++)
       {
         norm += Lx->get(j,i)*Lx->get(j,i)/pc_au2amu;
       }
       if(norm > 1e-3)
       {
       redmass->set(i,1.0/norm);
       }
    }

   //redmass->print();

   // Transform dipole-moment derivatives to normal coordinates //
   SharedMatrix dipder_q(new Matrix(3,3*natom));
   dipder_q->gemm(0,0,1.0,dipder,Lx,0.0);

   // Compute IR intensities in projected normal coordinates //
   double dipder_conv;
   dipder_conv = pc_dipmom_debye2si*pc_dipmom_debye2si/(1e-20 * pc_amu2kg * pc_au2amu);
   dipder_conv *= pc_na * pc_pi/(3.0 * pc_c * pc_c * 4.0 * pc_pi * pc_e0 * 1000.0);
   SharedVector IRint(new Vector("IRint",3*natom));
   for(i=0; i < natom*3; i++)
   {
     for(j=0; j < 3; j++)
     {
       IRint->add(i,dipder_conv * dipder_q->get(j,i) * dipder_q->get(j,i));
     }
   }


   // Transform polarizability derivatives to normal coordinates //
   SharedMatrix polder_q(new Matrix(9,3*natom));
   polder_q->gemm(1,0,1.0,polder,Lx,0.0);
   if(print >=2) {
     outfile->Printf("\n\tPolarizability Derivatives in Cart. Coord.\n");
     polder->print();
     outfile->Printf("\n\tPolarizability Derivatives in Normal Coord.\n");
     polder_q->print();
   }
   int jk=0;
   std::vector<SharedMatrix> alpha_der(3*natom);
   for(i=0; i < alpha_der.size(); ++i)
   {
     alpha_der[i] = SharedMatrix(new Matrix(3,3));
   }
   SharedMatrix alpha_der_mat(new Matrix(3,3));

   for(i=0; i < natom*3; i++)
   {
     for(j=0,jk=0; j < 3; j++)
     {
       for(k=0; k < 3; k++,jk++)
       {
         alpha_der[i]->set(j,k,polder_q->get(jk,i));

       }
     }
   }

   // Transform optical rotation tensor derivatives to normal coordinates //
   SharedMatrix optder_q(new Matrix(9,3*natom));
   optder_q->gemm(1,0,1.0,optder,Lx,0.0);
   if(print >=2) {
     outfile->Printf("\n\tOptical Rotation Tensor Derivatives in Cart. Coord..\n");
     optder->print();
     outfile->Printf("\n\tOptical Rotation Tensor Derivatives in Normal Coord.\n");
     optder_q->print();
   }

   std::vector<SharedMatrix> G_der(3*natom);
   for(i=0; i < G_der.size(); ++i)
   {
     G_der[i] = SharedMatrix(new Matrix(3,3));
   }
   SharedMatrix G_der_mat(new Matrix(3,3));

   for(i=0; i < natom*3; i++)
   {
     for(j=0,jk=0; j < 3; j++)
     {
       for(k=0; k < 3; k++,jk++)
       {
         G_der[i]->set(j,k,optder_q->get(jk,i));
       }
     }
   }


   // Transform dipole/quarupole tensor derivatives to normal coordinates //
   SharedMatrix quadder_q(new Matrix(27,3*natom));
   quadder_q->gemm(1,0,1.0,quadder,Lx,0.0);
   if(print >=2) {
     outfile->Printf("\n\tDipole/Quadrupole Tensor Derivatives in Cart. Coord..\n");
     quadder->print();
     outfile->Printf("\n\tDipole/Quadrupole Tensor Derivatives in Normal Coord.\n");
     quadder_q->print();
   }

   int jkl,l;
   double**** Q_der = (double ****) malloc(natom*3*sizeof(double ***));
   for(i=0; i < natom*3; i++)
   {
     Q_der[i] = (double ***) malloc(3*sizeof(double **));
     for(j=0,jkl=0; j < 3; j++)
     {
       Q_der[i][j] = (double **) malloc(3*sizeof(double *));
       for(k=0; k < 3; k++)
       {
       Q_der[i][j][k] = (double *) malloc(3*sizeof(double));
         for(l=0; l < 3; l++,jkl++)
         {
           Q_der[i][j][k][l] = quadder_q->get(jkl,i);
         }
       }
     }
   }

   // Compute the Raman scattering activity //
   double raman_conv = pc_bohr2angstroms * pc_bohr2angstroms * pc_bohr2angstroms * pc_bohr2angstroms / pc_au2amu;
   double km_convert = pc_hartree2J/(pc_bohr2m * pc_bohr2m * pc_amu2kg * pc_au2amu);
   double cm_convert = 1.0/(2.0 * pc_pi * pc_c * 100.0);

   SharedVector alpha(new Vector("Alpha",3*natom));
   SharedVector betaalpha2(new Vector("BetaAlpha2",3*natom));
   SharedVector ramint_linear(new Vector("RamIntLinear",3*natom));
   SharedVector depol_linear(new Vector("DepolLinear",3*natom));
   SharedVector ramint_circular(new Vector("RamIntCircular",3*natom));
   SharedVector depol_circular(new Vector("DepolCircular",3*natom));

   for(i=0; i < natom*3; i++)
   {
     alpha->set(i,tensor_mean(alpha_der[i]));
     betaalpha2->set(i,beta_alpha2(alpha_der[i]));
     ramint_linear->set(i,raman_linear(alpha->get(i), betaalpha2->get(i)));
     depol_linear->set(i,depolar_linear(alpha->get(i), betaalpha2->get(i)));

     ramint_circular->set(i,raman_circular(alpha->get(i), betaalpha2->get(i)));
     depol_circular->set(i,depolar_circular(alpha->get(i), betaalpha2->get(i)));
   }

  /* compute the frequencies and spit them out in a nice table in the output file*/
  outfile->Printf("\n\t     Harmonic Freq.  IR Intensity   Red. Mass       Alpha^2        Beta^2    Raman Int.  Depol. Ratio\n");
  outfile->Printf("\t        (cm-1)         (km/mol)       (amu) \n");
  outfile->Printf("\t---------------------------------------------------------------------------------------------------\n");
  for(i=natom*3-1; i >= 0; i--)
  {
    if(Fevals->get(i) < 0.0)
      outfile->Printf("\t  %3d  %9.3fi    %9.4f       %7.4f      %9.4f     %9.4f    %9.4f    %9.4f\n",
       (natom*3-i), cm_convert*sqrt(-km_convert*Fevals->get(i)), IRint->get(i),
       redmass->get(i), alpha->get(i)*alpha->get(i)*raman_conv, betaalpha2->get(i)*raman_conv,
       ramint_linear->get(i)*raman_conv, depol_linear->get(i));
    else
      outfile->Printf("\t  %3d  %9.3f     %9.4f       %7.4f      %9.4f     %9.4f    %9.4f    %9.4f\n",
       (natom*3-i), cm_convert*sqrt(km_convert*Fevals->get(i)), IRint->get(i),
       redmass->get(i), alpha->get(i)*alpha->get(i)*raman_conv, betaalpha2->get(i)*raman_conv,
       ramint_linear->get(i)*raman_conv, depol_linear->get(i));
  }
  outfile->Printf("\t---------------------------------------------------------------------------------------------------\n");

  /* compute the frequencies and spit them out in a nice table */
  outfile->Printf("----------------------------------------------------------------------------------------------\n");
  outfile->Printf("                                Raman Scattering Parameters\n");
  outfile->Printf("----------------------------------------------------------------------------------------------\n");
  outfile->Printf("     Harmonic Freq.  Alpha^2   Beta^2    Raman Act.   Dep. Ratio  Raman Act.  Dep. Ratio\n");
  outfile->Printf("        (cm-1)                           (linear)      (linear)   (natural)   (natural)\n");
  outfile->Printf("----------------------------------------------------------------------------------------------\n");
  for(i=natom*3-1; i >= 0; i--)
  {
    if(Fevals->get(i) < 0.0)
      outfile->Printf("  %3d  %9.3fi %9.4f   %9.4f  %9.4f  %9.4f    %9.4f  %9.4f\n",
       (natom*3-i), cm_convert*sqrt(-km_convert*Fevals->get(i)),
       alpha->get(i)*alpha->get(i)*raman_conv, betaalpha2->get(i)*raman_conv,
       ramint_linear->get(i)*raman_conv, depol_linear->get(i),
       ramint_circular->get(i)*raman_conv, depol_circular->get(i));
    else
      outfile->Printf("  %3d  %9.3f  %9.4f   %9.4f  %9.4f  %9.4f    %9.4f  %9.4f\n",
       (natom*3-i), cm_convert*sqrt(km_convert*Fevals->get(i)),
       alpha->get(i)*alpha->get(i)*raman_conv, betaalpha2->get(i)*raman_conv,
       ramint_linear->get(i)*raman_conv, depol_linear->get(i),
       ramint_circular->get(i)*raman_conv, depol_circular->get(i));
  }
  outfile->Printf("----------------------------------------------------------------------------------------------\n");

  SharedVector G(new Vector("G",3*natom));
  SharedVector betaG2(new Vector("betaG2",3*natom));
  SharedVector betaA2(new Vector("betaA2",3*natom));
  G->zero();
  betaG2->zero();
  betaA2->zero();

  for(i=0; i < natom*3; i++) {
    G->set(i,tensor_mean(G_der[i]));
    betaG2->set(i,beta_G2(alpha_der[i], G_der[i]));
    betaA2->set(i,beta_A2(alpha_der[i], Q_der[i], omega));
  }


  double cvel = pc_c * pc_me * pc_bohr2angstroms * 1e-10 / (pc_h/(2.0*pc_pi));
  double roa_conv = raman_conv * 1e6 / cvel;
//  roa_conv /= _c * _me * _bohr2angstroms * 1e-10;
//  roa_conv *= _h /(2.0 * _pi);
  outfile->Printf("-----------------------------------------------------\n");
  outfile->Printf("               ROA Scattering Invariants\n");
  outfile->Printf("-----------------------------------------------------\n");
  outfile->Printf("     Harmonic Freq.  AlphaG   Beta(G)^2  Beta(A)^2\n");
  outfile->Printf("        (cm-1)\n");
  outfile->Printf("-----------------------------------------------------\n");
  for(i=natom*3-1; i >= 0; i--) {
    if(Fevals->get(i) < 0.0)
      outfile->Printf("  %3d  %9.3fi %9.4f   %10.4f  %10.4f\n",
       (natom*3-i), cm_convert*sqrt(-km_convert*Fevals->get(i)),
       alpha->get(i)*G->get(i)*roa_conv, betaG2->get(i)*roa_conv, betaA2->get(i)*roa_conv);
    else
      outfile->Printf("  %3d  %9.3f  %9.4f   %10.4f  %10.4f\n",
       (natom*3-i), cm_convert*sqrt(km_convert*Fevals->get(i)),
       alpha->get(i)*G->get(i)*roa_conv, betaG2->get(i)*roa_conv, betaA2->get(i)*roa_conv);
  }
  outfile->Printf("-----------------------------------------------------\n");

  outfile->Printf("\n----------------------------------------------------------------------\n");
  outfile->Printf("         ROA Difference Parameter R-L (Angstrom^4/amu * 1000) \n");
  outfile->Printf("----------------------------------------------------------------------\n");
  outfile->Printf("     Harmonic Freq. Delta_z(90) Delta_x(90)   Delta(0)  Delta(180)\n");
  outfile->Printf("        (cm-1)\n");
  outfile->Printf("----------------------------------------------------------------------\n");

  double delta_0,delta_180,delta_x,delta_z;
  for(i=natom*3-1; i >= 0; i--)
  {

    delta_0 = 4.0 * (180.0*alpha->get(i)*G->get(i) + 4.0*betaG2->get(i) - 4.0*betaA2->get(i));
    delta_0 *= 1.0 / cvel;
    delta_180 = 4.0 * (24.0 * betaG2->get(i) + 8.0 * betaA2->get(i));
    delta_180 *= 1.0 / cvel;
    delta_x = 4.0 * (45.0 * alpha->get(i) * G->get(i) + 7.0 * betaG2->get(i) + betaA2->get(i));
    delta_x *= 1.0 / cvel;
    delta_z = 4.0 * (6.0 * betaG2->get(i) - 2.0 * betaA2->get(i));
    delta_z *= 1.0 / cvel;
    if(Fevals->get(i) < 0.0)
      outfile->Printf("  %3d  %9.3fi %10.4f   %10.4f   %10.4f   %10.4f\n",
       (natom*3-i), cm_convert*sqrt(-km_convert*Fevals->get(i)),
       delta_z*raman_conv*1e3, delta_x*raman_conv*1e3, delta_0*raman_conv*1e3, delta_180*raman_conv*1e3);
    else
      outfile->Printf("  %3d  %9.3f  %10.4f   %10.4f   %10.4f   %10.4f\n",
       (natom*3-i), cm_convert*sqrt(km_convert*Fevals->get(i)),
       delta_z*raman_conv*1e3, delta_x*raman_conv*1e3, delta_0*raman_conv*1e3, delta_180*raman_conv*1e3);
  }

}

// Handy Tensor Derivative Array Printer
void print_tensor_der(std::shared_ptr<OutFile> myfile, std::vector<SharedMatrix> my_tensor_list)
{
  for(int i=0; i < my_tensor_list.size(); ++i)  {
    int atom_num  = i/3;
    int xyz       = i%3;
    if(xyz==0) myfile->Printf( "\tAtom #%d, X-coord.:\n", atom_num);
    if(xyz==1) myfile->Printf( "\tAtom #%d, Y-coord.:\n", atom_num);
    if(xyz==2) myfile->Printf( "\tAtom #%d, Z-coord.:\n", atom_num);
    my_tensor_list[i]->print("myfile");
  }
}

}} // namespace psi::ccresponse


/* The Levi-Civitas evaluator */

int levi(int a, int b, int c)
{
  int val=0;
  int x=0, y=1, z=2;

  if(a==x && b==y && c==z) val=1;
  else if(a==y && b==z && c==x) val=1;
  else if(a==z && b==x && c==y) val=1;
  else if(a==x && b==z && c==y) val=-1;
  else if(a==y && b==x && c==z) val=-1;
  else if(a==z && b==y && c==x) val=-1;
  else val=0;

  return val;
}

/* compute mean of a property tensor: alpha = 1/3 alpha_ii */
double tensor_mean(SharedMatrix alpha)
{
  double mean=0.0;
  int i;
  for(i=0; i < 3; i++)
    mean += alpha->get(i,i);
  mean /= 3.0;
  return mean;
}

/* compute beta(alpha)^2 = 1/2 [ 3 * alpha_ij*alpha_ij - alpha_ii*alpha_jj */
double beta_alpha2(SharedMatrix alpha)
{
  double value = 0.0;
  int i,j;
  for(i=0; i < 3; i++)
    for(j=0; j < 3; j++)
      value += 0.5*(3.0*alpha->get(i,j)*alpha->get(i,j) - alpha->get(i,i)*alpha->get(j,j));

  return value;
}

/* compute beta(G')^2 = 1/2[3 * alpha_ij*G'_ij - alpha_ii*G'_jj */
double beta_G2(SharedMatrix alpha, SharedMatrix G)
{
  double value = 0.0;
  int i,j;
  for(i=0; i < 3; i++)
    for(j=0; j < 3; j++)
      value += 0.5*(3.0*alpha->get(i,j)*G->get(i,j) - alpha->get(i,i)*G->get(j,j));
  return value;
}

/* compute beta(A)^2 = 1/2 omega * alpha_ij epsilon_ikl * A_klj */
double beta_A2(SharedMatrix alpha, double ***A, double omega)
{
  double value=0.0;
  int i,j,k,l;
  for(i=0; i < 3; i++)
    for(j=0; j < 3; j++)
      for(k=0; k < 3; k++)
        for(l=0; l < 3; l++)
          value += 0.5 * omega * alpha->get(i,j) * levi(i,k,l) * A[k][l][j];

  return value;
}

/* compute Raman intensity for linearly polarized radiation:
    A = 45 (alpha^2) + 4 (beta^2)
*/
double raman_linear(double alpha, double beta2)
{
  double value = 0.0;
  value = 45.0 * alpha * alpha + 4.0 * beta2;
  return value;
}

/* compute Raman depolarization ratio for 90-degree scattering of linearly
   polarized radiation:

  ratio = [ 3 * beta(alpha)^2)/(45 * alpha^2 + 4 * beta(alpha)^2) ]
*/
double depolar_linear(double alpha, double beta2)
{
  double numer, denom;

  numer = 3.0 * beta2;
  denom = (45.0 * alpha * alpha) + (4.0 * beta2);

  if(denom > 1e-6) return numer/denom;
  else return 0.0;
}

/* compute Raman intensity for circularly polarized radiation:
    A = 45 (alpha^2) + 7 (beta^2);
*/
double raman_circular(double alpha, double beta2)
{
  double value = 0.0;
  value = 45.0 * alpha * alpha + 7.0 * beta2;
  return value;
}

// compute Raman depolarization ratio for 90-degree scattering of circularly polarized radiation://
// ratio = [ 6 * beta(alpha)^2)/(45 * alpha^2 + 7 * beta(alpha)^2) ] //

double depolar_circular(double alpha, double beta2)
{
  double numer, denom;

  numer = 6.0 * beta2;
  denom = (45.0 * alpha * alpha) + (7.0 * beta2);

  if(denom > 1e-6) return numer/denom;
  else return 0.0;
}

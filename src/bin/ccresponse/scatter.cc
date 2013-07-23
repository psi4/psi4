/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
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
#include <libciomr/libciomr.h>
#include <libpsio/psio.h> 
#include <libqt/qt.h>
#include "physconst.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"
#include <vector>
#include <psi4-dec.h>
#include <libmints/mints.h>
#include "mass.h"

using namespace boost;
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

void print_tensor_der(FILE *myfile, std::vector<SharedMatrix> my_tensor_list);

//void scatter(Options &options, std::vector <SharedMatrix> pol, std::vector <SharedMatrix> rot, std::vector <SharedMatrix> quad)
void scatter(double step, std::vector <SharedMatrix> pol, std::vector <SharedMatrix> rot, std::vector <SharedMatrix> quad)
{
    //double step = options.get_double("DISP_SIZE");
    printf("STEPSIZE = %lf bohr\n\n",step);
	int quiet = 1;
    int i,j,k;
    int a,b,c,d;
    double omega = 0.085645;

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
	FILE *derivs = fopen("tender.dat", "w");
	fprintf(derivs, "******************************************************\n");
	fprintf(derivs, "**********                                  **********\n");
	fprintf(derivs, "**********        TENSOR DERIVATIVES        **********\n");
	fprintf(derivs, "**********   FOR COMPUTING ROA SCATTERING   **********\n");
	fprintf(derivs, "**********                                  **********\n");
	fprintf(derivs, "******************************************************\n\n\n");
	fprintf(derivs, "\t*** Dipole Polarizability Derivative Tensors ***\n\n");
	print_tensor_der(derivs, pol_grad);
	fprintf(derivs, "*********************************************************\n\n");
	fprintf(derivs, "\t*** Optical Rotation Derivative Tensors ***\n\n");
	print_tensor_der(derivs, rot_grad);
	fprintf(derivs, "*********************************************************\n\n");
	fprintf(derivs, "\t*** Dipole/Quadrupole Derivative Tensors ***\n\n");
	print_tensor_der(derivs, quad_grad);
	fclose(derivs);

    boost::shared_ptr<Molecule> molecule = Process::environment.molecule();
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
        fscanf(hessian,"%lf",&F->pointer()[i][j]);
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
        fscanf(dipole_moment, "%lf %lf %lf", &dipder->pointer()[i][j*3], &dipder->pointer()[i][j*3+1], &dipder->pointer()[i][j*3+2]);
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
	if(!quiet)  {
      printf("\tInput coordinates (bohr):\n");
	  molecule->geometry().print(stdout);
    }
	molecule->move_to_com();
	geom_orig->copy(molecule->geometry());
    
    // Reading the Z-vals //
    int zvals[natom];
    for(i=0; i<natom; i++)
    {
        zvals[i]=molecule->Z(i);
    }
    
    // Mass-weighting the co-ordinates //
    double massi[natom];

    for(i=0; i<natom; i++)
    {
        massi[i] = molecule->mass(i);
        for(j=0; j<3; j++)
        {
            geom->set(i,j, geom_orig->get(i,j) * sqrt(massi[i]));
        }
    }
	if(!quiet)  {
	  printf("\tMass-Weighted coordinates relative to the center of mass:\n");
	  geom->print(stdout);
    }
 
    //Generating the inertia tensor
   
    SharedMatrix I(new Matrix(3,3));
    I->copy(molecule->inertia_tensor());
    //I = molecule->inertia_tensor();
	if(!quiet)  {
	  printf("\tMoment of Inertia Tensor:\n");
      I->print(stdout);
	  fflush(stdout);
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
   
    SharedMatrix Ievecs(new Matrix(3,3));
    SharedVector Ievals(new Vector("Ieigenval",3));

    I->diagonalize(Ievecs,Ievals);
	//rs(3,3,I->pointer(),Ievals->pointer(),1,Ievecs->pointer(),1e-12);

    // Constructing I-inverse matrix //

    SharedMatrix Iinv(new Matrix(3,3));
    SharedMatrix Itmp(new Matrix(3,3));
     
    Iinv->zero();
    for(i=0;i<3;i++)
    {
      Iinv->set(i,i,(1.0/Ievals->get(i)));
    }

    Itmp->gemm(0,1,1.0,Iinv,Ievecs,0.0);
    Iinv->gemm(0,0,1.0,Ievecs,Itmp,0.0);

	if(!quiet)  {
	  printf("\tIegv, Ievec, I^-1 matrix\n");
	  Ievecs->print(stdout);
	  Ievals->print(stdout);
      Iinv->print(stdout);
	  fflush(stdout);
    }
 
    // Generating the 6 pure rotation and translation vectors //
    SharedMatrix P(new Matrix(natom*3,natom*3));
    int icart,jcart,iatom,jatom;
    double imass,jmass,total_mass;
  
    total_mass=0.0;
    for(i=0;i<natom;i++)
    {
      //total_mass+=an2mass[(int)zvals[iatom]];
      total_mass+=an2mass[(int)zvals[i]];
    }

    for(i=0; i < natom*3; i++)
    {
      icart = i % 3;
      iatom = i/3;
      imass = an2mass[(int)zvals[iatom]];

      P->set(i,i,1.0);

      for(j=0; j < natom*3; j++)
      {
        jcart = j % 3;
        jatom = j/3;
        jmass = an2mass[(int)zvals[jatom]];

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

    // Generate mass-weighted Hessian matrix (Hartree/(bohr^2 amu)) //
    M->zero();    
    SharedMatrix T(new Matrix(natom*3,natom*3));
    for(i=0; i < natom; i++)
    {
      for(j=0; j < 3; j++)
      {
        M->set((i*3+j),(i*3+j),1/sqrt(an2mass[(int)zvals[i]]/_au2amu));
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
	if(!quiet)  {
	  printf("\tProjected, Mass-Weighted Hessian:\n");
	  F->print(stdout);
	  fflush(stdout);
    }

    // Diagonalize projected mass-weighted Hessian //
    SharedMatrix Fevecs(new Matrix(3*natom,3*natom));
    SharedVector Fevals(new Vector("Feigenval",3*natom));
    SharedMatrix Lx(new Matrix("Normal Transform Matrix",3*natom,3*natom));
    SharedVector redmass(new Vector("ReducedMass",3*natom));
    double norm=0.0;

    F->diagonalize(Fevecs,Fevals);
	//rs(3*natom,3*natom,F->pointer(),Fevals->pointer(),1,Fevecs->pointer(),1e-12);
	if(!quiet)  {
	  Fevals->print(stdout);
	  Fevecs->print(stdout);
      fflush(stdout);
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
	if(!quiet)  {
	  Lx->print(stdout);
      fflush(stdout);
    }

    for(i=0; i < 3*natom; i++)
    {
       norm = 0.0;
       for(j=0; j < 3*natom; j++)
       {
         norm += Lx->get(j,i)*Lx->get(j,i)/_au2amu;
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
   dipder_conv = _dipmom_debye2si*_dipmom_debye2si/(1e-20 * _amu2kg * _au2amu);
   dipder_conv *= _na * _pi/(3.0 * _c * _c * 4.0 * _pi * _e0 * 1000.0);
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
   fprintf(outfile,"\n\tPolarizability Derivatives, normal\n");
   //polder->print(stdout);
   polder_q->print();
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
   fprintf(outfile,"\tOptical Rotation Tensor Derivatives, normal\n");
   //optder->print(stdout);
   optder_q->print();
  
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
   fprintf(outfile,"\tDipole/Quadrupole Tensor Derivatives, normal\n");
   //quadder->print(stdout);
   quadder_q->print();
 
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
   double raman_conv = _bohr2angstroms * _bohr2angstroms * _bohr2angstroms * _bohr2angstroms / _au2amu;
   double km_convert = _hartree2J/(_bohr2m * _bohr2m * _amu2kg * _au2amu);
   double cm_convert = 1.0/(2.0 * _pi * _c * 100.0);

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

  /* compute the frequencies and spit them out in a nice table */
  printf("\n\t     Harmonic Freq.  IR Intensity   Red. Mass       Alpha^2        Beta^2    Raman Int.  Depol. Ratio\n");
  printf(  "\t        (cm-1)         (km/mol)       (amu) \n");
  printf("\t---------------------------------------------------------------------------------------------------\n");
  for(i=natom*3-1; i >= 0; i--)
  {
    if(Fevals->get(i) < 0.0)
      printf("\t  %3d  %9.3fi    %9.4f       %7.4f      %9.4f     %9.4f    %9.4f    %9.4f\n",
       (natom*3-i), cm_convert*sqrt(-km_convert*Fevals->get(i)), IRint->get(i),
       redmass->get(i), alpha->get(i)*alpha->get(i)*raman_conv, betaalpha2->get(i)*raman_conv,
       ramint_linear->get(i)*raman_conv, depol_linear->get(i));
    else
      printf("\t  %3d  %9.3f     %9.4f       %7.4f      %9.4f     %9.4f    %9.4f    %9.4f\n",
       (natom*3-i), cm_convert*sqrt(km_convert*Fevals->get(i)), IRint->get(i),
       redmass->get(i), alpha->get(i)*alpha->get(i)*raman_conv, betaalpha2->get(i)*raman_conv,
       ramint_linear->get(i)*raman_conv, depol_linear->get(i));
  }
  printf("\t---------------------------------------------------------------------------------------------------\n");

  /* compute the frequencies and spit them out in a nice table */
  printf("----------------------------------------------------------------------------------------------\n");
  printf("                                Raman Scattering Parameters\n");
  printf("----------------------------------------------------------------------------------------------\n");
  printf("     Harmonic Freq.  Alpha^2   Beta^2    Raman Act.   Dep. Ratio  Raman Act.  Dep. Ratio\n");
  printf("        (cm-1)                           (linear)      (linear)   (natural)   (natural)\n");
  printf("----------------------------------------------------------------------------------------------\n");
  for(i=natom*3-1; i >= 0; i--)
  {
    if(Fevals->get(i) < 0.0)
      printf("  %3d  %9.3fi %9.4f   %9.4f  %9.4f  %9.4f    %9.4f  %9.4f\n",
       (natom*3-i), cm_convert*sqrt(-km_convert*Fevals->get(i)),
       alpha->get(i)*alpha->get(i)*raman_conv, betaalpha2->get(i)*raman_conv,
       ramint_linear->get(i)*raman_conv, depol_linear->get(i),
       ramint_circular->get(i)*raman_conv, depol_circular->get(i));
    else
      printf("  %3d  %9.3f  %9.4f   %9.4f  %9.4f  %9.4f    %9.4f  %9.4f\n",
       (natom*3-i), cm_convert*sqrt(km_convert*Fevals->get(i)),
       alpha->get(i)*alpha->get(i)*raman_conv, betaalpha2->get(i)*raman_conv,
       ramint_linear->get(i)*raman_conv, depol_linear->get(i),
       ramint_circular->get(i)*raman_conv, depol_circular->get(i));
  }
  printf("----------------------------------------------------------------------------------------------\n");

  SharedVector G(new Vector("G",3*natom));
  SharedVector betaG2(new Vector("betaG2",3*natom));
  SharedVector betaA2(new Vector("betaA2",3*natom));
  
  for(i=0; i < natom*3; i++) {
    G->set(i,tensor_mean(G_der[i]));
    betaG2->set(i,beta_G2(alpha_der[i], G_der[i]));
    betaA2->set(i,beta_A2(alpha_der[i], Q_der[i], omega));
  }

  double cvel;
  double roa_conv;

  cvel = _c * _me * _bohr2angstroms * 1e-10 / (_h/(2.0*_pi));
  printf("cvel = %20.14f\n", cvel);
  roa_conv = raman_conv * 1e6 / cvel;
//  roa_conv /= _c * _me * _bohr2angstroms * 1e-10;
//  roa_conv *= _h /(2.0 * _pi);
  printf("-----------------------------------------------------\n");
  printf("               ROA Scattering Invariants\n");
  printf("-----------------------------------------------------\n");
  printf("     Harmonic Freq.  AlphaG   Beta(G)^2  Beta(A)^2\n");
  printf("        (cm-1)\n");
  printf("-----------------------------------------------------\n");
  for(i=natom*3-1; i >= 0; i--) {
    if(Fevals->get(i) < 0.0)
      printf("  %3d  %9.3fi %9.4f   %10.4f  %10.4f\n",
       (natom*3-i), cm_convert*sqrt(-km_convert*Fevals->get(i)),
       alpha->get(i)*G->get(i)*roa_conv, betaG2->get(i)*roa_conv, betaA2->get(i)*roa_conv);
    else
      printf("  %3d  %9.3f  %9.4f   %10.4f  %10.4f\n",
       (natom*3-i), cm_convert*sqrt(km_convert*Fevals->get(i)),
       alpha->get(i)*G->get(i)*roa_conv, betaG2->get(i)*roa_conv, betaA2->get(i)*roa_conv);
  }
  printf("-----------------------------------------------------\n");

  printf("\n----------------------------------------------------------------------\n");
  printf("         ROA Difference Parameter R-L (Angstrom^4/amu * 1000) \n");
  printf("----------------------------------------------------------------------\n");
  printf("     Harmonic Freq. Delta_z(90) Delta_x(90)   Delta(0)  Delta(180)\n");
  printf("        (cm-1)\n");
  printf("----------------------------------------------------------------------\n");
  
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
      printf("  %3d  %9.3fi %10.4f   %10.4f   %10.4f   %10.4f\n",
       (natom*3-i), cm_convert*sqrt(-km_convert*Fevals->get(i)),
       delta_z*raman_conv*1e3, delta_x*raman_conv*1e3, delta_0*raman_conv*1e3, delta_180*raman_conv*1e3);
    else
      printf("  %3d  %9.3f  %10.4f   %10.4f   %10.4f   %10.4f\n",
       (natom*3-i), cm_convert*sqrt(km_convert*Fevals->get(i)),
       delta_z*raman_conv*1e3, delta_x*raman_conv*1e3, delta_0*raman_conv*1e3, delta_180*raman_conv*1e3);
  }

}

	void print_tensor_der(FILE *myfile, std::vector<SharedMatrix> my_tensor_list)
	{
        for(int i=0; i < my_tensor_list.size(); ++i)  {
            int atom_num  = i/3;
            int xyz       = i%3;
            if(xyz==0) fprintf(myfile, "\tAtom #%d, X-coord.:\n", atom_num);
            if(xyz==1) fprintf(myfile, "\tAtom #%d, Y-coord.:\n", atom_num);
            if(xyz==2) fprintf(myfile, "\tAtom #%d, Z-coord.:\n", atom_num);
            my_tensor_list[i]->print(myfile);
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



#define DSIGN(a,b) ((b) >= 0.0) ? (fabs(a)) : (-fabs(a))
void tred2(int, double **, double *, double *, int);
//int tqli(int, double *, double **, double *, int, double);
void tqli(int, double *, double **, double *, int, double);
void new_eigsort(double *, double **, int);

double *init_array(int n)
{
  int i;
  double *A;
  A = (double*) malloc (n * sizeof(double));
  for(i=0; i < n; i++) A[i] = 0.0;
  return(A);
}

double **init_matrix(int n, int m)
{
  int i, j;
  double **A;
  A = (double**) malloc (n * sizeof(double*));
  for (i=0; i < n; i++) {
    A[i] = (double*) malloc (m * sizeof(double));
    for(j=0; j < m; j++) A[i][j] = 0.0;
  }
  return(A);
}

void free_matrix(double **A, int n)
{
  int i;
  for (i=0; i < n; i++) free(A[i]);
  free(A);
}


/*
** rs(): Diagonalize a symmetric matrix and return its eigenvalues and 
** eigenvectors.
**
** int nm: The row/column dimension of the matrix.
** int n:  Ditto.
** double **array: The matrix
** double *e_vals: A vector, which will contain the eigenvalues.
** int matz: Boolean for returning the eigenvectors.
** double **e_vecs: A matrix whose columns are the eigenvectors.
** double toler: A tolerance limit for the iterative procedure.  Rec: 1e-13
**
*/


/* translation into c of a translation into FORTRAN77 of the EISPACK */
/* matrix diagonalization routines */

void rs(int nm, int n, double **array, double *e_vals, int matz,
        double **e_vecs, double toler)
   {
      int i, j, ii, ij, ierr;
      int ascend_order;
      double *fv1, **temp;
      double zero = 0.0;
      double one = 1.0;

/* Modified by Ed - matz can have the values 0 through 3 */
      
      if ((matz > 3) || (matz < 0)) {
        matz = 0;
        ascend_order = 1;
        }
      else
        if (matz < 2)
          ascend_order = 1;	/* Eigenvalues in ascending order */
        else {
          matz -= 2;
          ascend_order = 0;	/* Eigenvalues in descending order */
          }

      fv1 = init_array(n);
      temp = init_matrix(n,n);

      if (n > nm) {
         ierr = 10*n;
         fprintf(stderr,"n = %d is greater than nm = %d in rsp\n",n,nm);
         exit(ierr);
         }

      for (i=0; i < n; i++) {
         for (j=0; j < n; j++) {
            e_vecs[i][j] = array[i][j];
            }
          }

      tred2(n,e_vecs,e_vals,fv1,matz);

      for (i=0; i < n; i++)
         for (j=0; j < n; j++)
            temp[i][j]=e_vecs[j][i];

      tqli(n,e_vals,temp,fv1,matz,toler);

      for (i=0; i < n; i++)
         for (j=0; j < n; j++)
            e_vecs[i][j]=temp[j][i];

      if (ascend_order)
        new_eigsort(e_vals,e_vecs,n);
      else
        new_eigsort(e_vals,e_vecs,(-1)*n);

      free(fv1);
      free_matrix(temp,n);

      }


/* converts symmetric matrix to a tridagonal form for use in tqli */
/* if matz = 0, only find eigenvalues, else find both eigenvalues and */
/* eigenvectors */

void tred2(int n, double **a, double *d, double *e, int matz)
   {
      int i,j,k,l,il,ik,jk,kj;
      double f,g,h,hh,scale,scale_inv,h_inv;
      double temp;

      if (n == 1) return;

      for (i=n-1; i > 0; i--) {
         l = i-1;
         h = 0.0;
         scale = 0.0;
         if (l) {
            for (k=0; k <= l; k++) {
               scale += fabs(a[i][k]);
               }
            if (scale == 0.0) {
               e[i] = a[i][l];
               }
            else {
               scale_inv=1.0/scale;
               for (k=0; k <= l; k++) {
                  a[i][k] *= scale_inv;
                  h += a[i][k]*a[i][k];
                  }
               f=a[i][l];
               g= -(DSIGN(sqrt(h),f));
               e[i] = scale*g;
               h -= f*g;
               a[i][l] = f-g;
               f = 0.0;
               h_inv=1.0/h;
               for (j=0; j <= l; j++) {
                  if (matz) a[j][i] = a[i][j]*h_inv;
                  g = 0.0;
                  for (k=0; k <= j; k++) {
                     g += a[j][k]*a[i][k];
                     }
                  if (l > j) {
                     for (k=j+1; k <= l; k++) {
                        g += a[k][j]*a[i][k];
                        }
                     }
                  e[j] = g*h_inv;
                  f += e[j]*a[i][j];
                  }
               hh = f/(h+h);
               for (j=0; j <= l; j++) {
                  f = a[i][j];
                  g = e[j] - hh*f;
                  e[j] = g;
                  for (k=0; k <= j; k++) {
                     a[j][k] -= (f*e[k] + g*a[i][k]);
                     }
                  }
               }
            }
         else {
            e[i] = a[i][l];
            }
         d[i] = h;
         }
      if (matz) d[0] = 0.0;
      e[0] = 0.0;

      for (i=0; i < n; i++) {
         l = i-1;
         if (matz) {
            if (d[i]) {
               for (j=0; j <= l; j++) {
                  g = 0.0;
                  for (k=0; k <= l; k++) {
                     g += a[i][k]*a[k][j];
                     }
                  for (k=0; k <= l; k++) {
                     a[k][j] -= g*a[k][i];
                     }
                  }
               }
            }
         d[i] = a[i][i];
         if (matz) {
            a[i][i] = 1.0;
            if (l >= 0) {
               for (j=0; j<= l; j++) {
                  a[i][j] = 0.0;
                  a[j][i] = 0.0;
                  }
               }
            }
         }
      }

/* diagonalizes tridiagonal matrix output by tred2 */
/* gives only eigenvalues if matz = 0, both eigenvalues and eigenvectors */
/* if matz = 1 */

//int tqli(int n, double *d, double **z, double *e, int matz, double toler)
void tqli(int n, double *d, double **z, double *e, int matz, double toler)
   {
      register int k;
      int i,j,l,m,iter;
      double dd,g,r,s,c,p,f,b,h;
      double azi;

      f=0.0;
      if (n == 1) {
         d[0]=z[0][0];
         z[0][0] = 1.0;
         return;
         }

      for (i=1; i < n ; i++) {
         e[i-1] = e[i];
         }
      e[n-1] = 0.0;
      for (l=0; l < n; l++) {
         iter = 0;
L1:
         for (m=l; m < n-1;m++) {
            dd = fabs(d[m]) + fabs(d[m+1]);
#if 0
            if (fabs(e[m])+dd == dd) goto L2;
#else
            if (fabs(e[m]) < toler) goto L2;
#endif
            }
         m=n-1;
L2:
         if (m != l) {
            if (iter++ == 30) {
               fprintf (stderr,"tqli not converging\n");
                continue;
#if 0
               exit(30);
#endif
               }

            g = (d[l+1]-d[l])/(2.0*e[l]);
            r = sqrt(g*g + 1.0);
            g = d[m] - d[l] + e[l]/((g + DSIGN(r,g)));
            s=1.0;
            c=1.0;
            p=0.0;
            for (i=m-1; i >= l; i--) {
               f = s*e[i];
               b = c*e[i];
               if (fabs(f) >= fabs(g)) {
                  c = g/f;
                  r = sqrt(c*c + 1.0);
                  e[i+1] = f*r;
                  s=1.0/r;
                  c *= s;
                  }
               else {
                  s = f/g;
                  r = sqrt(s*s + 1.0);
                  e[i+1] = g*r;
                  c = 1.0/r;
                  s *= c;
                  }
               g = d[i+1] - p;
               r = (d[i]-g)*s + 2.0*c*b;
               p = s*r;
               d[i+1] = g+p;
               g = c*r-b;

               if (matz) {
                  double *zi = z[i];
                  double *zi1 = z[i+1];
                  for (k=n; k ; k--,zi++,zi1++) {
                     azi = *zi;
                     f = *zi1;
                     *zi1 = azi*s + c*f;
                     *zi = azi*c - s*f;
                     }
                  }
               }

            d[l] -= p;
            e[l] = g;
            e[m] = 0.0;
            goto L1;
            }
         }
   }

void new_eigsort(double *d, double **v, int n)
{
      int i,j,k;
      double p;

/// Modified by Ed - if n is negative - sort eigenvalues in descending order

      if (n >= 0) {
        for (i=0; i < n-1 ; i++) {
           k=i;
           p=d[i];
           for (j=i+1; j < n; j++) {
              if (d[j] < p) {
                 k=j;
                 p=d[j];
                 }
              }
           if (k != i) {
              d[k]=d[i];
              d[i]=p;
              for (j=0; j < n; j++) {
                 p=v[j][i];
                 v[j][i]=v[j][k];
                 v[j][k]=p;
                 }
               }
            }
        }
      else {
        n = abs(n);
        for (i=0; i < n-1 ; i++) {
           k=i;
           p=d[i];
           for (j=i+1; j < n; j++) {
              if (d[j] > p) {
                 k=j;
                 p=d[j];
                 }
              }
           if (k != i) {
              d[k]=d[i];
              d[i]=p;
              for (j=0; j < n; j++) {
                  p=v[j][i];
                  v[j][i]=v[j][k];
                  v[j][k]=p;
                  }
              }
           }
        }
   }

void print_mat(double **a, int m, int n,FILE *out)
{
      int ii,jj,kk,nn,ll;
      int i,j,k;

      ii=0;jj=0;
L200:
      ii++;
      jj++;
      kk=10*jj;
      nn=n;
      if (nn > kk) nn=kk;
      ll = 2*(nn-ii+1)+1;
      fprintf (out,"\n");
      for (i=ii; i <= nn; i++) fprintf(out,"       %5d",i);
      fprintf (out,"\n");
      for (i=0; i < m; i++) {
         fprintf (out,"\n%5d",i+1);
         for (j=ii-1; j < nn; j++) {
            fprintf (out,"%12.7f",a[i][j]);
            }
         }
      fprintf (out,"\n");
      if (n <= kk) {
         fflush(out);
         return;
         }
      ii=kk; goto L200;
}

/*!
** eigout(): Print out eigenvectors and eigenvalues.  Prints an n x m
** matrix of eigenvectors.  Under each eigenvector, the corresponding
** elements of an array, b, will also be printed.  This is useful for
** printing, for example, the SCF eigenvectors with their associated
** eigenvalues (orbital energies).
*/
void eigout(double **a, double *b, int m, int n, FILE *out)
   {
      int ii,jj,kk,nn;
      int i,j;

      ii=0;jj=0;
L200:
      ii++;
      jj++;
      kk=10*jj;
      nn=n;
      if (nn > kk) nn=kk;
      fprintf (out,"\n");
      for (i=ii; i <= nn; i++) fprintf(out,"       %5d",i);
      fprintf (out,"\n");
      for (i=0; i < m; i++) {
         fprintf (out,"\n%5d",i+1);
         for (j=ii-1; j < nn; j++) {
            fprintf (out,"%12.7f",a[i][j]);
            }
         }
      fprintf (out,"\n");
      fprintf (out,"\n     ");
      for (j=ii-1; j < nn; j++) {
         fprintf(out,"%12.7f",b[j]);
         }
      fprintf (out,"\n");
      if (n <= kk) {
         fflush(out);
         return;
         }
      ii=kk; goto L200;
      }

double dot(double *A, double *B, int len) 
{
  int i;
  double value = 0.0;

  for(i=0; i < len; i++) value += A[i] * B[i];
  return value;
}


/*!
**                                                             
** mmult():
** a reasonably fast matrix multiply (at least on the DEC3100) 
** written by ETS                                              
**                                                             
** AF,BF,and CF are fortran arrays                             
**                                                             
** ta,tb and tc indicate whether the corresponding arrays are  
**              to be converted to their transpose             
**                                                             
** nr,nl,nc are the number of rows,links,and columns in the    
**          final matrices to be multiplied together           
**          if ta=0 AF should have the dimensions nr x nl      
**          if ta=1 AF should have the dimensions nl x nr      
**          if tb=0 BF should have the dimensions nl x nc      
**          if tb=1 BF should have the dimensions nc x nl      
**          if tc=0 CF should have the dimensions nr x nc      
**          if tc=1 CF should have the dimensions nc x nr      
**                                                             
** add is 1 if this matrix is to be added to the one passed    
**        in as CF, 0 otherwise                                
**
** \ingroup (CIOMR)
*/
void mmult(double **AF, int ta, double **BF, int tb, double **CF, int tc,
	   int nr, int nl, int nc, int add)
{
  int odd_nr,odd_nc,odd_nl;
  int i,j,k,ij;
  double t00,t01,t10,t11;
  double **a,**b;
  double *att,*bt;
  double *at1,*bt1;
  static int keep_nr=0;
  static int keep_nl=0;
  static int keep_nc=0;
  static double **aa,**bb;

  if(!aa) {
    aa = (double **) init_matrix(nr,nl);
    bb = (double **) init_matrix(nc,nl);
    keep_nr = nr;
    keep_nl = nl;
    keep_nc = nc;
  }

  if(nl > keep_nl) {
    free_matrix(aa,keep_nr);
    free_matrix(bb,keep_nc);
    keep_nl = nl;
    keep_nr = (nr > keep_nr) ? nr : keep_nr;
    keep_nc = (nc > keep_nc) ? nc : keep_nc;
    aa = (double **) init_matrix(keep_nr,keep_nl);
    bb = (double **) init_matrix(keep_nc,keep_nl);
  }
  if(nr > keep_nr) {
    free_matrix(aa,keep_nr);
    keep_nr = nr;
    aa = (double **) init_matrix(keep_nr,keep_nl);
  }
  if(nc > keep_nc) {
    free_matrix(bb,keep_nc);
    keep_nc = nc;
    bb = (double **) init_matrix(keep_nc,keep_nl);
  }

  odd_nr = (nr)%2;
  odd_nc = (nc)%2;
  odd_nl = (nl)%2;

  a=aa;
  if(ta)
    for(i=0; i < nr ; i++)
      for(j=0; j < nl ; j++)
	a[i][j] = AF[j][i];
  else
    a=AF;

  b=bb;
  if(tb)
    b=BF;
  else
    for(i=0; i < nc ; i++)
      for(j=0; j < nl ; j++)
	b[i][j] = BF[j][i];
      
  for(j=0; j < nc-1 ; j+=2) {
    for(i=0; i < nr-1 ; i+=2) {
      att=a[i]; bt=b[j];
      at1=a[i+1]; bt1=b[j+1];
      if(add) {
	if(tc) {
	  t00 = CF[j][i];
	  t01 = CF[j+1][i];
	  t10 = CF[j][i+1];
	  t11 = CF[j+1][i+1];
	}
	else {
	  t00 = CF[i][j];
	  t01 = CF[i][j+1];
	  t10 = CF[i+1][j];
	  t11 = CF[i+1][j+1];
	}
      }
      else t00=t01=t10=t11=0.0;
      for(k=nl; k ; k--,att++,bt++,at1++,bt1++) {
	t00 += *att * *bt;
	t01 += *att * *bt1;
	t10 += *at1 * *bt;
	t11 += *at1 * *bt1;
      }
      if(tc) {
	CF[j][i]=t00;
	CF[j+1][i]=t01;
	CF[j][i+1]=t10;
	CF[j+1][i+1]=t11;
      }
      else {
	CF[i][j]=t00;
	CF[i][j+1]=t01;
	CF[i+1][j]=t10;
	CF[i+1][j+1]=t11;
      }
    }
    if(odd_nr) {
      att=a[i]; bt=b[j];
      bt1=b[j+1];
      if(add) {
	if(tc) {
	  t00 = CF[j][i];
	  t01 = CF[j+1][i];
	}
	else {
	  t00 = CF[i][j];
	  t01 = CF[i][j+1];
	}
      }
      else t00=t01=0.0;
      for(k= nl; k ; k--,att++,bt++,bt1++) {
	t00 += *att * *bt;
	t01 += *att * *bt1;
      }
      if(tc) {
	CF[j][i]=t00;
	CF[j+1][i]=t01;
      }
      else {
	CF[i][j]=t00;
	CF[i][j+1]=t01;
      }
    }
  }
  if(odd_nc) {
    for(i=0; i < nr-1 ; i+=2) {
      att=a[i]; bt=b[j];
      at1=a[i+1];
      if(add) {
	if(tc) {
	  t00 = CF[j][i];
	  t10 = CF[j][i+1];
	}
	else {
	  t00 = CF[i][j];
	  t10 = CF[i+1][j];
	}
      }
      else t00=t10=0.0;
      for(k= nl; k ; k--,att++,bt++,at1++) {
	t00 += *att * *bt;
	t10 += *at1 * *bt;
      }
      if(tc) {
	CF[j][i]=t00;
	CF[j][i+1]=t10;
      }
      else {
	CF[i][j]=t00;
	CF[i+1][j]=t10;
      }
    }
    if(odd_nr) {
      att=a[i]; bt=b[j];
      if(add)
	t00 = (tc) ? CF[j][i] : CF[i][j];
      else t00=0.0;
      for(k=nl; k ; k--,att++,bt++)
	t00 += *att * *bt;
      if(tc) CF[j][i]=t00;
      else CF[i][j]=t00;
    }
  }
}

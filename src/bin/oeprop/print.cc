/*! \file
    \ingroup OEPROP
    \brief Enter brief description of file here 
*/
#define EXTERN
#include "includes.h"
#include "prototypes.h"
#include "globals.h"

namespace psi { namespace oeprop {

double F_Zalpha(double Z_alpha);

void print_intro()
{ 
  fprintf(outfile,"    **********************************************\n");
  fprintf(outfile,"    *                    OEPROP                  *\n");
  fprintf(outfile,"    *          A simple property program         *\n");
  fprintf(outfile,"    *              by a big TOOL fan             *\n");
  fprintf(outfile,"    **********************************************\n\n");
}



void print_tasks()
{ 
   fprintf(outfile,"\n  TASKS to be performed :\n");
   
   if (read_opdm) {
     fprintf(outfile,"    $One-particle density in %s basis in %s form will be read from file%d",
             opdm_basis.c_str(),opdm_format.c_str(),opdm_file);
     if (asymm_opdm)
       fprintf(outfile," and symmetrized.\n");
     else
       fprintf(outfile,".\n");

     if (wrtnos)
       fprintf(outfile,"    $Natural orbitals will be written to checkpoint.\n");
   }
   else
     fprintf(outfile,"    $One-particle density will be computed from the eigenvector in checkpoint.\n");

   if (spin_prop)
     fprintf(outfile,"    $Spin properties will be evaluated.\n");
   
   switch (mpmax) {
     
     case 1:
       fprintf(outfile,"    $Only electric dipole moment will be computed.\n");
       break;
     
     case 2:
       fprintf(outfile,"    $Electric dipole and quadrupole moments will be computed.\n");
       break;
     
     case 3:
       fprintf(outfile,"    $Electric dipole, quadrupole, and octopole moments will be computed.\n");
       break;
   }
   
   fprintf(outfile,"    $Reference point for the electric multipole moments calculation is ");
   switch (mp_ref) {
     
     case 1: fprintf(outfile,"\n      the center of mass.\n");
             break;
             
     case 2: fprintf(outfile,"\n      the origin of the coordinate system.\n");
             break;

     case 3: fprintf(outfile,"\n      the center of electronic charge computed from a Mulliken analysis.\n");
             break;
     
     case 4: fprintf(outfile,"\n      the center of nuclear charge.\n");
             break;
     
     case 5: fprintf(outfile,"\n      the center of net charge computed from a Mulliken analysis.\n");
             break;
     
     default: fprintf(outfile,"\n      at (%lf %lf %lf)\n",
                      mp_ref_xyz[0],mp_ref_xyz[1],mp_ref_xyz[2]);
   }
  if (corr)
    fprintf(outfile,"    $Correlation corrections to electric multipole moments will be computed.\n");
  fprintf(outfile,"    $Reference point for the electric angular momentum calculation is ");
  fprintf(outfile,"\n      at (%lf %lf %lf)\n",
	  Lm_ref_xyz[0],Lm_ref_xyz[1],Lm_ref_xyz[2]);
  if (nuc_esp)
    fprintf(outfile,"    $Electrostatic properties at the nuclei will be evaluated.\n");
  if (grid) {
    switch(grid) {
      case 1: fprintf(outfile,"    $Electrostatic potential ");
  	      break;
      case 2: if (!spin_prop)
	        fprintf(outfile,"    $Electron density ");
              else
	        fprintf(outfile,"    $Spin density ");
	      break;
      case 3: if (!spin_prop)
	        fprintf(outfile,"    $Electron density gradient ");
              else
	        fprintf(outfile,"    $Spin density gradient ");
	      break;
      case 4: if (!spin_prop)
	        fprintf(outfile,"    $Laplacian of electron density ");
              else
	        fprintf(outfile,"    $Laplacian of spin density ");
	      break;
    }
    fprintf(outfile,"will be evaluated over a rectangular %dx%d grid.\n",nix+1,niy+1);
  }

}



void print_params()
{
   int i;

   fprintf(outfile,"\n  Title : '%s'\n",title.c_str());
   fprintf(outfile,"\n  List of PARAMETERS :\n");
   fprintf(outfile,"    # of atoms                 = %6d\n",natom);
   fprintf(outfile,"    # of molecular orbitals    = %6d\n",nmo);
   fprintf(outfile,"    # of basis functions       = %6d\n",nbfso);
   fprintf(outfile,"    # of atomic orbitals       = %6d\n",nbfao);
   fprintf(outfile,"    # of irreps                = %6d\n",nirreps);
   fprintf(outfile,"    Total charge               = %6d\n",charge);
   fprintf(outfile,"    # of unique shells         = %6d\n",nshell);
   fprintf(outfile,"    # of primitives            = %6d\n",nprim);
   fprintf(outfile,"    Print level                = %6d\n",print_lvl);
   if (fine_structure_alpha != 1.0)
     fprintf(outfile,"    Fine-structure a/a_0       =\t  %8.5lf\n",fine_structure_alpha);
   if (grid3d == 0) {
     fprintf(outfile,"\n  List of GRID PARAMETERS :\n");
     fprintf(outfile,"    GRID_ORIGIN                =\t( %8.5lf %8.5lf %8.5lf )\n",grid_origin[0],grid_origin[1],grid_origin[2]);
     fprintf(outfile,"    GRID_STEP_X                =\t( %8.5lf %8.5lf %8.5lf )\n",grid_step_x[0],grid_step_x[1],grid_step_x[2]);
     fprintf(outfile,"    GRID_STEP_Y                =\t( %8.5lf %8.5lf %8.5lf )\n",grid_step_y[0],grid_step_y[1],grid_step_y[2]);
     fprintf(outfile,"    NIX                        =\t  %d\n",nix+1);
     fprintf(outfile,"    NIY                        =\t  %d\n",niy+1);
     fprintf(outfile,"    GRID_ZMIN                  =\t  %8.5lf\n",grid_zmin);
     fprintf(outfile,"    GRID_ZMAX                  =\t  %8.5lf\n",grid_zmax);
   }
   else if (grid3d) {
     fprintf(outfile,"\n  List of GRID PARAMETERS :\n");
     fprintf(outfile,"    GRID_ORIGIN                =\t( %8.5lf %8.5lf %8.5lf )\n",grid_origin[0],grid_origin[1],grid_origin[2]);
     fprintf(outfile,"    GRID_STEP_X                =\t( %8.5lf %8.5lf %8.5lf )\n",grid_step_x[0],grid_step_x[1],grid_step_x[2]);
     fprintf(outfile,"    GRID_STEP_Y                =\t( %8.5lf %8.5lf %8.5lf )\n",grid_step_y[0],grid_step_y[1],grid_step_y[2]);
     fprintf(outfile,"    GRID_STEP_Z                =\t( %8.5lf %8.5lf %8.5lf )\n",grid_step_z[0],grid_step_z[1],grid_step_z[2]);
     fprintf(outfile,"    NIX                        =\t  %d\n",nix+1);
     fprintf(outfile,"    NIY                        =\t  %d\n",niy+1);
     fprintf(outfile,"    NIZ                        =\t  %d\n",niz+1);
   }
   fprintf(outfile,"\n");
}


void print_pop_header()
{
  fprintf(outfile," --------------------------------------------------------------\n");
  fprintf(outfile,"   ** Mulliken population analysis of one-particle density **\n");
  fprintf(outfile," --------------------------------------------------------------\n\n");
}


void print_mp()
{
  FILE *file_dipmom;

  fprintf(outfile," --------------------------------------------------------------\n");
  fprintf(outfile,"                *** Electric multipole moments ***\n");
  fprintf(outfile," --------------------------------------------------------------\n\n");
  if (charge != 0) {
    fprintf(outfile,"  CAUTION : The system has non-vanishing charge, therefore dipole\n");
    fprintf(outfile,"    and higher moments depend on the reference point. \n\n");
  }
  else
  if ((dtot > 1.0E-15) && (mpmax > 1)) {
    fprintf(outfile,"  CAUTION : The system has non-vanishing dipole moment, therefore\n");
    fprintf(outfile,"    quadrupole and higher moments depend on the reference point.\n\n");
  }
  else
  if (((qvals[0]*qvals[0] + qvals[1]*qvals[1] + qvals[2]*qvals[2]) > 1.0E-15) 
      && (mpmax > 2)) {
    fprintf(outfile,"  CAUTION : The system has non-vanishing quadrupole moment, therefore\n");
    fprintf(outfile,"    octopole and higher moments depend on the reference point.\n\n");
  }

  fprintf(outfile," -Coordinates of the reference point (a.u.) :\n");
  fprintf(outfile,"           x                     y                     z\n");
  fprintf(outfile,"  --------------------  --------------------  --------------------\n");
  fprintf(outfile,"  %20.10lf  %20.10lf  %20.10lf\n\n",
          mp_ref_xyz[0],mp_ref_xyz[1],mp_ref_xyz[2]);
  if (print_lvl >= PRINTDIPCOMPLEVEL) {
    fprintf(outfile," -Contributions to electric dipole moment (a.u.) :\n\n");
    fprintf(outfile,"   -Electronic part :\n\n");
    fprintf(outfile,"    mu(X) =  %11.8lf  mu(Y) =  %11.8lf  mu(Z) = %11.8lf\n\n",
            dx_e,dy_e,dz_e);
    fprintf(outfile,"   -Nuclear part :\n\n");
    fprintf(outfile,"    mu(X) =  %11.8lf  mu(Y) =  %11.8lf  mu(Z) = %11.8lf\n\n",
            dx_n,dy_n,dz_n);
  }
  fprintf(outfile," -Electric dipole moment (expectation values) :\n\n");
  fprintf(outfile,"    mu(X)  =  %8.5lf D  =  %15.8e C*m  =  %11.8lf a.u.\n",
          dx*_dipmom_au2debye,dx*_dipmom_au2si,dx);
  fprintf(outfile,"    mu(Y)  =  %8.5lf D  =  %15.8e C*m  =  %11.8lf a.u.\n",
          dy*_dipmom_au2debye,dy*_dipmom_au2si,dy);
  fprintf(outfile,"    mu(Z)  =  %8.5lf D  =  %15.8e C*m  =  %11.8lf a.u.\n",
          dz*_dipmom_au2debye,dz*_dipmom_au2si,dz);
  fprintf(outfile,"    |mu|   =  %8.5lf D  =  %15.8e C*m  =  %11.8lf a.u.\n",
          dtot*_dipmom_au2debye,dtot*_dipmom_au2si,dtot);
  if (corr) {
    fprintf(outfile,"\n -Correlation correction to electric dipole moment :\n\n");
    fprintf(outfile,"    cc(X)  =  %8.5lf D  =  %15.8e C*m  =  %11.8lf a.u.\n",
            dxcc*_dipmom_au2debye,dxcc*_dipmom_au2si,dxcc);
    fprintf(outfile,"    cc(Y)  =  %8.5lf D  =  %15.8e C*m  =  %11.8lf a.u.\n",
            dycc*_dipmom_au2debye,dycc*_dipmom_au2si,dycc);
    fprintf(outfile,"    cc(Z)  =  %8.5lf D  =  %15.8e C*m  =  %11.8lf a.u.\n",
            dzcc*_dipmom_au2debye,dzcc*_dipmom_au2si,dzcc);
    fprintf(outfile,"\n -Corrected electric dipole moment :\n\n");
    fprintf(outfile,"   mu(X)   =  %8.5lf D  =  %15.8e C*m  =  %11.8lf a.u.\n",
            (dx+dxcc)*_dipmom_au2debye,(dx+dxcc)*_dipmom_au2si,(dx+dxcc));
    fprintf(outfile,"   mu(Y)   =  %8.5lf D  =  %15.8e C*m  =  %11.8lf a.u.\n",
            (dy+dycc)*_dipmom_au2debye,(dy+dycc)*_dipmom_au2si,(dy+dycc));
    fprintf(outfile,"   mu(Z)   =  %8.5lf D  =  %15.8e C*m  =  %11.8lf a.u.\n",
            (dz+dzcc)*_dipmom_au2debye,(dz+dzcc)*_dipmom_au2si,(dz+dzcc));
    dtot = sqrt((dx+dxcc)*(dx+dxcc) + (dy+dycc)*(dy+dycc) + (dz+dzcc)*(dz+dzcc));
    fprintf(outfile,"    |mu|   =  %8.5lf D  =  %15.8e C*m  =  %11.8lf a.u.\n",
            dtot*_dipmom_au2debye,dtot*_dipmom_au2si,dtot);
  }

  /* write the total dipole moment to an ASCII file */
  if (wrt_dipmom) { 
    ffile(&file_dipmom,"dipmom.dat",1);
    fprintf(file_dipmom,"%20.10lf%20.10lf%20.10lf\n",
      (dx+dzcc)*_dipmom_au2debye,(dy+dycc)*_dipmom_au2debye,
      (dz+dzcc)*_dipmom_au2debye);
    fclose(file_dipmom);
  }


  if (mpmax > 1) {
    fprintf(outfile,"\n -Components of electric quadrupole moment (expectation values) (a.u.) :\n\n");
    fprintf(outfile,"     Q(XX) =  %12.8lf   Q(YY) =  %12.8lf   Q(ZZ) =  %12.8lf\n",
            qxx,qyy,qzz);
    fprintf(outfile,"     Q(XY) =  %12.8lf   Q(XZ) =  %12.8lf   Q(YZ) =  %12.8lf\n",
            qxy,qxz,qyz);
    if (corr) {
      fprintf(outfile,"\n -Correlation correction to electric quadrupole moment (a.u.) :\n\n");
      fprintf(outfile,"    cc(XX) =  %12.8lf  cc(YY) =  %12.8lf  cc(ZZ) =  %12.8lf\n",
              qxxcc,qyycc,qzzcc);
      fprintf(outfile,"    cc(XY) =  %12.8lf  cc(XZ) =  %12.8lf  cc(YZ) =  %12.8lf\n",
              qxycc,qxzcc,qyzcc);
      fprintf(outfile,"\n -Principal values (a.u.) and axis of corrected electric quadrupole moment :\n\n");
    }
    else
      fprintf(outfile,"\n -Principal values (a.u.) and axis of electric quadrupole moment :\n\n");
    fprintf(outfile,"    Q1     =  %12.8lf      V1 = (%11.8lf %11.8lf %11.8lf)\n",
            qvals[0],qvecs[0][0],qvecs[1][0],qvecs[2][0]);
    fprintf(outfile,"    Q2     =  %12.8lf      V2 = (%11.8lf %11.8lf %11.8lf)\n",
            qvals[1],qvecs[0][1],qvecs[1][1],qvecs[2][1]);
    fprintf(outfile,"    Q3     =  %12.8lf      V3 = (%11.8lf %11.8lf %11.8lf)\n",
            qvals[2],qvecs[0][2],qvecs[1][2],qvecs[2][2]);
  }
  if (mpmax > 2) {
    fprintf(outfile,"\n -Components of electric octopole moment (expectation values) (a.u.) :\n\n");
    fprintf(outfile,"    O(XXX) =  %12.8lf  O(XXY) =  %12.8lf  O(XXZ) =  %12.8lf\n",
            oxxx,oxxy,oxxz);
    fprintf(outfile,"    O(YYY) =  %12.8lf  O(XYY) =  %12.8lf  O(YYZ) =  %12.8lf\n",
            oyyy,oxyy,oyyz);
    fprintf(outfile,"    O(ZZZ) =  %12.8lf  O(XZZ) =  %12.8lf  O(YZZ) =  %12.8lf\n",
            ozzz,oxzz,oyzz);
    fprintf(outfile,"                            O(XYZ) =  %12.8lf \n",oxyz);
    if (corr) {
      fprintf(outfile,"\n -Correlation correction to electric octopole moment (a.u.) :\n\n");
      fprintf(outfile,"   cc(XXX) =  %12.8lf cc(XXY) =  %12.8lf cc(XXZ) =  %12.8lf\n",
              oxxxcc,oxxycc,oxxzcc);
      fprintf(outfile,"   cc(YYY) =  %12.8lf cc(XYY) =  %12.8lf cc(YYZ) =  %12.8lf\n",
              oyyycc,oxyycc,oyyzcc);
      fprintf(outfile,"   cc(ZZZ) =  %12.8lf cc(XZZ) =  %12.8lf cc(YZZ) =  %12.8lf\n",
              ozzzcc,oxzzcc,oyzzcc);
      fprintf(outfile,"                           cc(XYZ) =  %12.8lf \n",oxyzcc);
      fprintf(outfile,"\n -Corrected electric octopole moment (a.u.) :\n\n");
      fprintf(outfile,"    O(XXX) =  %12.8lf  O(XXY) =  %12.8lf  O(XXZ) =  %12.8lf\n",
              oxxx+oxxxcc,oxxy+oxxycc,oxxz+oxxzcc);
      fprintf(outfile,"    O(YYY) =  %12.8lf  O(XYY) =  %12.8lf  O(YYZ) =  %12.8lf\n",
              oyyy+oyyycc,oxyy+oxyycc,oyyz+oyyzcc);
      fprintf(outfile,"    O(ZZZ) =  %12.8lf  O(XZZ) =  %12.8lf  O(YZZ) =  %12.8lf\n",
              ozzz+ozzzcc,oxzz+oxzzcc,oyzz+oyzzcc);
      fprintf(outfile,"                            O(XYZ) =  %12.8lf \n",oxyz+oxyzcc);
    }
  }
  fprintf(outfile,"\n\n");
}


void print_lm()
{
  fprintf(outfile," --------------------------------------------------------------\n");
  fprintf(outfile,"              *** Electronic angular momentum ***\n");
  fprintf(outfile," --------------------------------------------------------------\n\n");

  fprintf(outfile," -Coordinates of the reference point (a.u.) :\n");
  fprintf(outfile,"           x                     y                     z\n");
  fprintf(outfile,"  --------------------  --------------------  --------------------\n");
  fprintf(outfile,"  %20.10lf  %20.10lf  %20.10lf\n\n",
          Lm_ref_xyz[0],Lm_ref_xyz[1],Lm_ref_xyz[2]);
  fprintf(outfile," -Electronic angular momentum (expectation values) :\n\n");
  fprintf(outfile,"    Lx  =  %11.8lf a.u.\n", Lx);
  fprintf(outfile,"    Ly  =  %11.8lf a.u.\n", Ly);
  fprintf(outfile,"    Lz  =  %11.8lf a.u.\n", Lz);
  if (mpmax > 1) {
    fprintf(outfile,"    Lx^2 = %11.8lf a.u.\n", Lx2);
    fprintf(outfile,"    Ly^2 = %11.8lf a.u.\n", Ly2);
    fprintf(outfile,"    Lz^2 = %11.8lf a.u.\n", Lz2);
    fprintf(outfile,"    L^2 =  %11.8lf a.u.\n", Lx2+Ly2+Lz2);
  }
  fprintf(outfile,"\n\n");
}


void print_esp()
{
  int i;
  
  fprintf(outfile," --------------------------------------------------------------\n");
  fprintf(outfile,"      *** Electrostatic  properties at atomic centers ***\n");
  fprintf(outfile," --------------------------------------------------------------\n\n");
  fprintf(outfile," -Coordinates of atomic centers (a.u.):\n");
  fprintf(outfile,"    #   Charge           x");
  fprintf(outfile,"                     y                     z\n");
  fprintf(outfile,"   ---  ------  --------------------");
  fprintf(outfile,"  --------------------  --------------------\n");
  for(i=0;i<natom;i++)
    fprintf(outfile,"%5d%7d    %20.10lf  %20.10lf  %20.10lf\n",
                    i+1,(int)zvals[i],geom[i][0]+mp_ref_xyz[0],
                        geom[i][1]+mp_ref_xyz[1],
                        geom[i][2]+mp_ref_xyz[2]);
  fprintf(outfile,"\n\n");
  fprintf(outfile," -Electrostatic potential and electric field (a.u.) :\n\n");
  fprintf(outfile,"    Center         phi            Ex             Ey           Ez\n");
  fprintf(outfile,"    ------    ------------   ------------  ------------  ------------\n");
  for(i=0;i<natom;i++)
    fprintf(outfile,"%8d      %12.8lf   %12.8lf  %12.8lf  %12.8lf\n",
            i+1,phi[i],ex[i],ey[i],ez[i]);
  fprintf(outfile,"\n\n");
  fprintf(outfile," -Electric field gradient (regular form) (a.u.):\n\n");
  fprintf(outfile,"    Center           XX                    YY                    ZZ\n");
  fprintf(outfile,"    ------  --------------------  --------------------  --------------------\n");
  for(i=0;i<natom;i++)
    fprintf(outfile,"%8d    %20.8lf  %20.8lf  %20.8lf\n",
            i+1,dexx[i],deyy[i],dezz[i]);
  fprintf(outfile,"\n    Center           XY                    XZ                    YZ\n");
  fprintf(outfile,"    ------  --------------------  --------------------  --------------------\n");
  for(i=0;i<natom;i++)
    fprintf(outfile,"%8d    %20.8lf  %20.8lf  %20.8lf\n",
            i+1,dexy[i],dexz[i],deyz[i]);
  fprintf(outfile,"\n\n");
  fprintf(outfile," -Electric field gradient (traceless tensor form) (a.u.):\n\n");
  fprintf(outfile,"    Center        XX - RR/3             YY - RR/3             ZZ - RR/3\n");
  fprintf(outfile,"    ------  --------------------  --------------------  --------------------\n");
  for(i=0;i<natom;i++)
    fprintf(outfile,"%8d    %20.8lf  %20.8lf  %20.8lf\n",i+1,
            (2*dexx[i]-deyy[i]-dezz[i])/3,
            (2*deyy[i]-dexx[i]-dezz[i])/3,
            (2*dezz[i]-dexx[i]-deyy[i])/3);
  fprintf(outfile,"\n    Center           XY                    XZ                    YZ\n");
  fprintf(outfile,"    ------  --------------------  --------------------  --------------------\n");
  for(i=0;i<natom;i++)
    fprintf(outfile,"%8d    %20.8lf  %20.8lf  %20.8lf\n",
            i+1,dexy[i],dexz[i],deyz[i]);
  fprintf(outfile,"\n\n");

  if (spin_prop) {
    fprintf(outfile," -Dipole-dipole term in hyperfine coupling constant (tensor form) (a.u.):\n\n");
    fprintf(outfile,"    Center        XX - RR/3             YY - RR/3             ZZ - RR/3\n");
    fprintf(outfile,"    ------  --------------------  --------------------  --------------------\n");
    for(i=0;i<natom;i++)
      fprintf(outfile,"%8d    %20.8lf  %20.8lf  %20.8lf\n",
              i+1,(2*ahfsxx[i]-ahfsyy[i]-ahfszz[i])/3,(2*ahfsyy[i]-ahfsxx[i]-ahfszz[i])/3,(2*ahfszz[i]-ahfsxx[i]-ahfsyy[i])/3);
    fprintf(outfile,"\n    Center           XY                    XZ                    YZ\n");
    fprintf(outfile,"    ------  --------------------  --------------------  --------------------\n");
    for(i=0;i<natom;i++)
      fprintf(outfile,"%8d    %20.8lf  %20.8lf  %20.8lf\n",
	      i+1,ahfsxy[i],ahfsxz[i],ahfsyz[i]);
    fprintf(outfile,"\n\n");
  }

  if (spin_prop) {
    fprintf(outfile," -Electron and spin densities (a.u.):\n\n");
    fprintf(outfile,"    Center    Electron density         Spin density\n");
    fprintf(outfile,"    ------  --------------------  --------------------\n");
    for(i=0;i<natom;i++)
      fprintf(outfile,"%8d    %20.8lf  %20.8lf\n",i+1,edens[i],sdens[i]);
  }
  else {
    fprintf(outfile," -Electron density (a.u.):\n\n");
    fprintf(outfile,"    Center           rho\n");
    fprintf(outfile,"    ------  --------------------\n");
    for(i=0;i<natom;i++)
      fprintf(outfile,"%8d    %20.8lf\n",i+1,edens[i]);
  }
  fprintf(outfile,"\n\n");
}


void print_misc()
{
  int i,j,k;
  double tval, energy;

  fprintf(outfile," --------------------------------------------------------------\n");
  fprintf(outfile,"                *** Miscellaneous properties ***\n");
  fprintf(outfile," --------------------------------------------------------------\n\n");

  fprintf(outfile," -Relativistic MVD one-electron corrections to the energy (a.u.):\n\n");
  fprintf(outfile,"    Mass-velocity (p^4) term     :   %12.15lf\n",massveloc);
  fprintf(outfile,"    One-electron Darwin term     :   %12.15lf\n",darw);
  fprintf(outfile,"    Total one-electron MVD terms :   %12.15lf\n",massveloc+darw);
  fprintf(outfile,"\n");

  if ((print_lvl >= PRINTDARWINCOMPLEVEL) || (QED_darwin)) {
    fprintf(outfile," -One-electron Darwin term per atom : \n\n");
    for (i=0; i<natom; ++i)
      fprintf(outfile,"    Atom %d: %12.15lf\n", i+1, darw_per_atom[i]);
    fprintf(outfile,"\n");
  }

  if (QED_darwin) {
    fprintf(outfile," -QED correction relative to one-electron Darwin term per atom : \n\n");
    for(i=0;i<natom;i++) {
      tval = 2 * (fine_structure_alpha/_c_au) / _pi * \
        (F_Zalpha(zvals[i]*fine_structure_alpha) - 4.0/15.0);
      fprintf(outfile,"    Atom %d: %12.15lf\n", i+1, tval);
      darw_per_atom[i] += darw_per_atom[i] * tval;
    }
    fprintf(outfile,"\n");
    darw = 0.0;
    for(i=0;i<natom;i++) /* recompute darwin term to add in QED */
      darw += darw_per_atom[i];
    fprintf(outfile," -Updating darwin energy (%15.10lf) to include QED correction\n",darw);
  }
  free(darw_per_atom);

  if (update_energy_with_MVD) {
    chkpt_init(PSIO_OPEN_OLD);
    energy = chkpt_rd_etot();
    energy = energy + massveloc + darw;
    chkpt_wt_etot(energy);
    chkpt_close();
    fprintf(outfile," -Updating total energy (%15.10lf) in chkpt file with MVD correction\n", energy);
    if (QED_darwin) fprintf(outfile," -(yes, including QED term)\n");
  }

  if (mpmax > 1) {
    fprintf(outfile,"  NOTE : Spatial extents are computed with respect to the same reference point\n");
    fprintf(outfile,"         as multipole moments.\n\n");
    fprintf(outfile," -Electronic spatial extents (a.u.) :\n\n");
    fprintf(outfile,"     <X^2> =  %11.4lf    <Y^2> =  %11.4lf    <Z^2> =  %11.4lf\n",
            exp_x2,exp_y2,exp_z2);
    fprintf(outfile,"                             <R^2> =  %11.4lf\n",
            exp_x2+exp_y2+exp_z2);
    fprintf(outfile,"\n -Orbital spatial extents ");
    if (read_opdm && wrtnos) 
      fprintf(outfile,"of NOs contructed from onepdm in file%d ",opdm_file);
    else
      fprintf(outfile,"of MOs in checkpoint ");
    fprintf(outfile,"(a.u.) :\n\n");
    fprintf(outfile,"    MO #   Symm     <X^2>        <Y^2>        <Z^2>        <R^2>\n");
    fprintf(outfile,"   ------  ----  -----------  -----------  -----------  -----------\n");
    k = 0;
    for(i=0;i<nirreps;i++)
      for(j=0;j<orbspi[i];j++)
        fprintf(outfile,"   %4d     %3s   %9.4lf    %9.4lf    %9.4lf    %9.4lf\n",
                k+1,irr_labs[i],MOXX[k],MOYY[k],MOZZ[k],MOXX[k]+MOYY[k]+MOZZ[k++]);
    fprintf(outfile,"\n");
  }
}

}} // namespace psi::oeprop

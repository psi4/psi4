/*! \file
    \ingroup OEPROP
    \brief Enter brief description of file here 
*/
#define EXTERN
#include "includes.h"
#include "prototypes.h"
#include "globals.h"
#include <rgb.h>

namespace {
  void print_2d_summary();
  void print_3d_summary();
  void print_grid_plotmtv();
  void print_grid_megapovray();
  void create_megapovray_file();
  void print_grid_gausscube();
  FILE* grid_file;
}

namespace psi { namespace oeprop {

void print_grid()
{
  if (grid_format == "MEGAPOVPLUS")
    print_grid_megapovray();
  else if (grid_format == "PLOTMTV")
    print_grid_plotmtv();
  else if (grid_format == "GAUSSCUBE")
    print_grid_gausscube();

  return;
}

}} // namespace psi::oeprop

namespace {
using namespace psi::oeprop;
using namespace psi;

void print_2d_summary()
{
  fprintf(outfile," --------------------------------------------------------------\n");
  fprintf(outfile,"    *** Evaluating properties over a rectangular 2D grid ***\n");
  fprintf(outfile," --------------------------------------------------------------\n\n");
  fprintf(outfile," -Coordinates of the lower left, lower right, and upper left corners of\n");
  fprintf(outfile,"  the grid rectangle (a.u.):\n");
  fprintf(outfile,"    **            x");
  fprintf(outfile,"                     y                     z\n");
  fprintf(outfile,"   ----  --------------------");
  fprintf(outfile,"  --------------------  --------------------\n");
  fprintf(outfile,"    LL   %20.10lf  %20.10lf  %20.10lf\n",
	  grid_origin[0],grid_origin[1],grid_origin[2]);
  fprintf(outfile,"    LR   %20.10lf  %20.10lf  %20.10lf\n",
	  grid_origin[0]+grid_step_x[0]*nix,
	  grid_origin[1]+grid_step_x[1]*nix,
	  grid_origin[2]+grid_step_x[2]*nix);
  fprintf(outfile,"    UL   %20.10lf  %20.10lf  %20.10lf\n\n\n",
	  grid_origin[0]+grid_step_y[0]*niy,
	  grid_origin[1]+grid_step_y[1]*niy,
	  grid_origin[2]+grid_step_y[2]*niy);
}

void print_3d_summary()
{
  fprintf(outfile," --------------------------------------------------------------\n");
  fprintf(outfile,"    *** Evaluating properties over a rectangular 3D grid ***\n");
  fprintf(outfile," --------------------------------------------------------------\n\n");
  fprintf(outfile," -Coordinates of the lower left and upper right corners of\n");
  fprintf(outfile,"  the grid box (a.u.):\n");
  fprintf(outfile,"    **            x");
  fprintf(outfile,"                     y                     z\n");
  fprintf(outfile,"   ----  --------------------");
  fprintf(outfile,"  --------------------  --------------------\n");
  fprintf(outfile,"    LL   %20.10lf  %20.10lf  %20.10lf\n",
	  grid_origin[0],grid_origin[1],grid_origin[2]);
  fprintf(outfile,"    UR   %20.10lf  %20.10lf  %20.10lf\n",
	  grid_origin[0]+grid_step_x[0]*nix+grid_step_y[0]*niy+grid_step_z[0]*niz,
	  grid_origin[1]+grid_step_x[1]*nix+grid_step_y[1]*niy+grid_step_z[1]*niz,
	  grid_origin[2]+grid_step_x[2]*nix+grid_step_y[2]*niy+grid_step_z[2]*niz);
}

void print_grid_plotmtv()
{
  int i,j,k;
  double step_x, step_y, x, y;

  print_2d_summary();
  
  switch (grid) {
    case 1:
      grid_file = fopen("esp.dat","w");
      break;

    case 2:
      if (!spin_prop)
        grid_file = fopen("edens.dat","w");
      else
	grid_file = fopen("sdens.dat","w");
      break;
    
    case 3:
      if (!spin_prop)
        grid_file = fopen("edgrad.dat","w");
      else
	grid_file = fopen("sdgrad.dat","w");
      break;

    case 4:
      if (!spin_prop)
        grid_file = fopen("edlapl.dat","w");
      else
	grid_file = fopen("sdlapl.dat","w");
      break;

  }

/*    fprintf(grid_file,"%8.4lf %d %8.4lf  %8.4lf %d %8.4lf  %8.4lf %d %8.4lf\n",
            grid_xyz0[0],nix,grid_xyz1[0],
            grid_xyz0[1],niy,grid_xyz1[1],
            grid_xyz0[2],niz,grid_xyz1[2]); */
  switch (grid) {
    case 1:
    case 2:
    case 4:
      fprintf(grid_file,"$DATA = CONTOUR\n");
      fprintf(grid_file,"%% xmin = %lf xmax = %lf nx = %d\n",grid_xy0[0],grid_xy1[0],nix+1);
      fprintf(grid_file,"%% ymin = %lf ymax = %lf ny = %d\n",grid_xy0[1],grid_xy1[1],niy+1);
      fprintf(grid_file,"%% zmin = %lf zmax = %lf\n",grid_zmin,grid_zmax);
      fprintf(grid_file,"%% contfill = T meshplot = T\n");
      for(i=0;i<=niy;i++) {
        for(j=0;j<=nix;j++)
	  if (grid_pts[j][i] < grid_zmin)
            fprintf(grid_file," %lf ",grid_zmin);
	  else
	    if (grid_pts[j][i] > grid_zmax)
	      fprintf(grid_file," %lf ",grid_zmax);
	    else
	      fprintf(grid_file," %lf ",grid_pts[j][i]);
        fprintf(grid_file,"\n");
      }
      break; 
      
    case 3:
      fprintf(grid_file,"$DATA = VECTOR\n");
      fprintf(grid_file,"%% xmin = %lf xmax = %lf\n",grid_xy0[0],grid_xy1[0]);
      fprintf(grid_file,"%% ymin = %lf ymax = %lf\n",grid_xy0[1],grid_xy1[1]);
      fprintf(grid_file,"%% zmin = %lf zmax = %lf\n",0.0,0.0);
      fprintf(grid_file,"%% linecolor = 3\n");
      fprintf(grid_file,"%% xlog = off vscale = 0.2\n\n");
      step_x = (grid_xy1[0]-grid_xy0[0])/nix;
      step_y = (grid_xy1[1]-grid_xy0[1])/niy;
      for(i=0;i<=nix;i++) {
	x = grid_xy0[0] + step_x*i;
	for(j=0;j<=niy;j++) {
        y = grid_xy0[1] + step_y*j;
        if (fabs(grid_pts[i][j]) <= MAXDENSGRAD)
          fprintf(grid_file,"%9.5lf  %9.5lf  %9.5lf  %12.8lf  %12.8lf  %12.8lf\n",x,y,0.0,
                  grid_vecX[i][j],grid_vecY[i][j],0.0);
        }
      }
      break;
  }
      
  fprintf(grid_file,"$END\n");
  fclose(grid_file);
}


void print_grid_megapovray()
{
  int i,j,k;
  double step_x, step_y, step_z, x, y, z;

  print_3d_summary();

  /*--- Write out a data file ---*/
  switch (grid) {
  case 5:
  case 6:
      grid_file = fopen("mo.dat","w");
      break;
  }

  switch (grid) {
  case 5:
  case 6:
    for(k=0;k<=niz;k++)
      for(j=0;j<=niy;j++)
	for(i=0;i<=nix;i++) {
	  fprintf(grid_file,"%15.8lf\n",grid3d_pts[0][i][j][k]);
	}
      break;
  }

  fclose(grid_file);

  /*--- Write out a command file ---*/
  create_megapovray_file();
  
  return;
}

void create_megapovray_file()
{
  int i, j, z;
  double radius, midx, midy, midz;
  double camera_distance, frame_width, frame_height, frame_center[3];
  double dimx, dimy, dimz, maxdim;
  double **grid_geom, **Rgrid;
  FILE *mpvfile;

  dimx = grid_xyz1[0] - grid_xyz0[0];
  dimy = grid_xyz1[1] - grid_xyz0[1];
  dimz = grid_xyz1[2] - grid_xyz0[2];
#define MAX(x,y) (((x) > (y)) ? (x) : (y))
  maxdim = MAX(dimz,MAX(dimx,dimy));
  frame_width = maxdim + 2.0;
  frame_height = maxdim + 2.0;
  frame_center[0] = grid_origin[0] + 0.5*dimx;
  frame_center[1] = grid_origin[1] + 0.5*dimy;
  frame_center[2] = grid_origin[2] + 0.5*dimz;
  camera_distance = 3.0*frame_height;

  /* compute atomic coordinates in the grid coordinate system,
     since MegaPov simply assumes a simple XYZ grid */
  grid_geom = block_matrix(natom,3);
  Rgrid = block_matrix(3,3);
  for(i=0; i<3; i++) {
    Rgrid[i][0] = grid_unit_x[i];
    Rgrid[i][1] = grid_unit_y[i];
    Rgrid[i][2] = grid_unit_z[i];
  }
  mmult(geom,0,Rgrid,0,grid_geom,0,natom,3,3,0);
  free_block(Rgrid);
    
  compute_connectivity();

  mpvfile = fopen("mo.pov","w");
  fprintf(mpvfile,"// File: mo.pov\n");
  fprintf(mpvfile,"// Creator: oeprop (Psi 3.2)\n");
  fprintf(mpvfile,"// Version: MegaPov 0.5\n\n");

  fprintf(mpvfile,"#version unofficial MegaPov 0.5;\n\n");

  fprintf(mpvfile,"// Camera\n");
  fprintf(mpvfile,"camera {\n");
  fprintf(mpvfile,"    orthographic\n");
  fprintf(mpvfile,"    location <%lf, %lf, %lf>\n",
	  camera_distance,frame_center[1],frame_center[2]);  
  fprintf(mpvfile,"    look_at <%lf, %lf, %lf>\n",
	  frame_center[0],frame_center[1],frame_center[2]);
  fprintf(mpvfile,"    up <%lf, %lf, %lf>\n",
	  frame_center[0],frame_center[1],frame_center[2]+frame_height);
  fprintf(mpvfile,"    right <%lf, %lf, %lf>\n",
	  frame_center[0],frame_center[1]+frame_width,frame_center[2]);
  fprintf(mpvfile,"}\n\n");

  fprintf(mpvfile,"// Light\n");
  fprintf(mpvfile,"light_source {<%lf, %lf, %lf> color rgb <1.0, 1.0, 1.0>}\n\n",
	  camera_distance*10.0,camera_distance*(-10.0),camera_distance*(-10.0));

  fprintf(mpvfile,"// Background\n");
  fprintf(mpvfile,"background { color rgb <1.0, 1.0, 1.0>}\n\n");

  fprintf(mpvfile,"union {\n");
  fprintf(mpvfile,"// Objects\n");
  fprintf(mpvfile,"  // union finish\n");
  fprintf(mpvfile,"  #declare F = finish {specular 0.4 roughness 0.005 diffuse 0.8 ambient 0.2}\n");
  fprintf(mpvfile,"  // transparency of atomic surfaces \n");
  fprintf(mpvfile,"  #declare T = 0;\n");
  fprintf(mpvfile,"  // no_shadow\n");
  fprintf(mpvfile,"  // hollow\n\n");
  
  fprintf(mpvfile,"  // Atoms\n");
  for(i=0;i<natom;i++) {
    z = (int)zvals[i];
    fprintf(mpvfile,"  // Atom # %d, charge %d\n",i,z);
    fprintf(mpvfile,"  object {\n");
    fprintf(mpvfile,"    sphere { < %lf, %lf, %lf > ",grid_geom[i][0],grid_geom[i][1],grid_geom[i][2]);
    if (z <= 2)
      radius = 0.4;
    else if (z <= 10)
      radius = 0.5;
    else
      radius = 0.6;
    fprintf(mpvfile,"%lf }\n",radius);
    fprintf(mpvfile,"     texture {\n");
    if (z <= LAST_RGB_INDEX)
      fprintf(mpvfile,"       pigment {color rgbt <%lf, %lf, %lf, T>}\n",atomic_rgb[z][0],atomic_rgb[z][1],atomic_rgb[z][2]);
    else
      fprintf(mpvfile,"       pigment {color rgbt <%lf, %lf, %lf, T>}\n",atomic_rgb[0][0],atomic_rgb[0][1],atomic_rgb[0][2]);
    fprintf(mpvfile,"       finish {F}\n");
    fprintf(mpvfile,"     }\n");
    fprintf(mpvfile,"  }\n");
  }
  fprintf(mpvfile,"\n");

  fprintf(mpvfile,"  // Bonds\n");
  for(i=0;i<natom;i++)
    for(j=0;j<i;j++)
      if (connectivity[i][j]) {
	fprintf(mpvfile,"  // Bond between atoms %d and %d\n",i,j);
	midx = 0.5*(grid_geom[i][0] + grid_geom[j][0]);
	midy = 0.5*(grid_geom[i][1] + grid_geom[j][1]);
	midz = 0.5*(grid_geom[i][2] + grid_geom[j][2]);
	fprintf(mpvfile,"  object {\n");
	fprintf(mpvfile,"    cylinder { < %lf, %lf, %lf > < %lf, %lf, %lf > 0.15 }\n",
		grid_geom[i][0],grid_geom[i][1],grid_geom[i][2],
		midx,midy,midz);
	fprintf(mpvfile,"     texture {\n");
	z = (int) zvals[i];
	if (z <= LAST_RGB_INDEX)
	    fprintf(mpvfile,"       pigment {color rgbt <%lf, %lf, %lf, T>}\n",atomic_rgb[z][0],atomic_rgb[z][1],atomic_rgb[z][2]);
	else
	    fprintf(mpvfile,"       pigment {color rgbt <%lf, %lf, %lf, T>}\n",atomic_rgb[0][0],atomic_rgb[0][1],atomic_rgb[0][2]);
	fprintf(mpvfile,"       finish {F}\n");
	fprintf(mpvfile,"     }\n");
	fprintf(mpvfile,"  }\n");
	fprintf(mpvfile,"  object {\n");
	fprintf(mpvfile,"    cylinder { < %lf, %lf, %lf > < %lf, %lf, %lf > 0.15 }\n",
		grid_geom[j][0],grid_geom[j][1],grid_geom[j][2],
		midx,midy,midz);
	fprintf(mpvfile,"     texture {\n");
	z = (int) zvals[j];
	if (z <= LAST_RGB_INDEX)
	    fprintf(mpvfile,"       pigment {color rgbt <%lf, %lf, %lf, T>}\n",atomic_rgb[z][0],atomic_rgb[z][1],atomic_rgb[z][2]);
	else
	    fprintf(mpvfile,"       pigment {color rgbt <%lf, %lf, %lf, T>}\n",atomic_rgb[0][0],atomic_rgb[0][1],atomic_rgb[0][2]);
	fprintf(mpvfile,"       finish {F}\n");
	fprintf(mpvfile,"     }\n");
	fprintf(mpvfile,"  }\n");
      }
  
  fprintf(mpvfile,"  #declare ISO = 1;\n");
  fprintf(mpvfile,"  #declare NX = %d;\n",nix);
  fprintf(mpvfile,"  #declare NY = %d;\n",niy);
  fprintf(mpvfile,"  #declare NZ = %d;\n",niz);
  fprintf(mpvfile,"  #declare RX = %lf;\n",grid_origin[0]);
  fprintf(mpvfile,"  #declare RY = %lf;\n",grid_origin[1]);
  fprintf(mpvfile,"  #declare RZ = %lf;\n",grid_origin[2]);
  fprintf(mpvfile,"  #declare LX = %lf;\n",grid_xyz1[0] - grid_xyz0[0]);
  fprintf(mpvfile,"  #declare LY = %lf;\n",grid_xyz1[1] - grid_xyz0[1]);
  fprintf(mpvfile,"  #declare LZ = %lf;\n",grid_xyz1[2] - grid_xyz0[2]);
  fprintf(mpvfile,"  #declare LEVSCALE = 11;\n\n");

  fprintf(mpvfile,"  #ifdef (ISO)\n");
  fprintf(mpvfile,"  #declare wfun1 = function {\"data_3D_3\",\n");
  fprintf(mpvfile,"     <-LEVSCALE>, library \"i_dat3d\",  \"mo.dat\", <NX+1,NY+1,NZ+1,0>}\n");
  fprintf(mpvfile,"  #declare wfun2 = function {\"data_3D_3\",\n");
  fprintf(mpvfile,"      <LEVSCALE>, library \"i_dat3d\",  \"mo.dat\", <NX+1,NY+1,NZ+1,0>}\n");
  fprintf(mpvfile,"  union {\n");
  fprintf(mpvfile,"    isosurface {\n");
  fprintf(mpvfile,"      function {wfun1}\n");
  fprintf(mpvfile,"      contained_by { box {<0,0,0>,<NX,NY,NZ>} }\n");
  fprintf(mpvfile,"      threshold -1\n");
  fprintf(mpvfile,"      max_gradient 1.000\n");
  fprintf(mpvfile,"      eval\n");
  fprintf(mpvfile,"      texture { pigment { color rgbt <1.0,0.0,0.0,0.8> } }\n");
  fprintf(mpvfile,"    }\n");
  fprintf(mpvfile,"    isosurface {\n");
  fprintf(mpvfile,"      function {wfun2}\n");
  fprintf(mpvfile,"      contained_by { box {<0,0,0>,<NX,NY,NZ>} }\n");
  fprintf(mpvfile,"      threshold -1\n");
  fprintf(mpvfile,"      max_gradient 1.000\n");
  fprintf(mpvfile,"      eval\n");
  fprintf(mpvfile,"      texture { pigment { color rgbt <0.0,1.0,0.0,0.8> } }\n");
  fprintf(mpvfile,"    }\n");
  fprintf(mpvfile,"    translate <RX*NX/LX,RY*NY/LY,RZ*NZ/LZ>\n");
  fprintf(mpvfile,"    scale <LX/NX,LY/NY,LZ/NZ>\n");
  fprintf(mpvfile,"  }\n");
  fprintf(mpvfile,"  #end\n\n");
  
  fprintf(mpvfile,"  // Rotate\n");
  fprintf(mpvfile,"  rotate 0*z\n");
  fprintf(mpvfile,"}\n");

  fclose(mpvfile);

  free_block(grid_geom);

  return;
}


void print_grid_gausscube()
{
  int i,j,k;
  int atom, mo;
  int n_per_line;
  double step_x, step_y, step_z, x, y, z;

  print_3d_summary();

  /* file name */
  switch (grid) {
  case 5:
      grid_file = fopen("mo.cube","w");
      break;
  case 6:
      grid_file = fopen("dens.cube","w");
      break;
  }

  /*------------------------------------------------------------------------
    Gaussian Cube standard overhead as described in G94 programmer's manual
   ------------------------------------------------------------------------*/
  /* Comment and subcomment */
  fprintf(grid_file,"Gaussian Cube file created by OEPROP (Psi 3.2)\n");
  fprintf(grid_file,"Calculation title: %s\n",title.c_str());
  if (grid == 5)
    fprintf(grid_file,"%5d",-1*natom);
  else if (grid == 6)
    fprintf(grid_file,"%5d",natom);
  fprintf(grid_file,"%12.6lf%12.6lf%12.6lf\n", grid_origin[0], grid_origin[1], grid_origin[2]);
  fprintf(grid_file,"%5d%12.6lf%12.6lf%12.6lf\n", nix+1, grid_step_x[0], grid_step_x[1], grid_step_x[2]);
  fprintf(grid_file,"%5d%12.6lf%12.6lf%12.6lf\n", niy+1, grid_step_y[0], grid_step_y[1], grid_step_y[2]);
  fprintf(grid_file,"%5d%12.6lf%12.6lf%12.6lf\n", niz+1, grid_step_z[0], grid_step_z[1], grid_step_z[2]);
  for(atom=0;atom<natom;atom++) {
    fprintf(grid_file,"%5d%12.6lf%12.6lf%12.6lf%12.6lf\n", (int)zvals[atom], zvals[atom],
	    geom[atom][0], geom[atom][1], geom[atom][2]);
  }
  if (grid == 5) {
    fprintf(grid_file,"%5d",num_mos_to_plot);
    for(mo=0;mo<num_mos_to_plot;mo++)
      fprintf(grid_file,"%5d",mos_to_plot[mo]);
    fprintf(grid_file,"\n");
  }

  switch (grid) {
  case 5:
    for(i=0;i<=nix;i++)
      for(j=0;j<=niy;j++) {
	n_per_line = 0;
	for(k=0;k<=niz;k++)
	  for(mo=0;mo<num_mos_to_plot;mo++) {
	    fprintf(grid_file,"%13.5E",grid3d_pts[mo][i][j][k]);
	    n_per_line++;
	    if (n_per_line == 6) {
	      fprintf(grid_file,"\n");
	      n_per_line = 0;
	    }
	  }
	if (n_per_line != 0)
	  fprintf(grid_file,"\n");
      }
    break;

  case 6:
    for(i=0;i<=nix;i++)
      for(j=0;j<=niy;j++) {
	n_per_line = 0;
	for(k=0;k<=niz;k++) {
	  fprintf(grid_file,"%13.5E",grid3d_pts[0][i][j][k]);
	  n_per_line++;
	  if (n_per_line == 6) {
	    fprintf(grid_file,"\n");
	    n_per_line = 0;
	  }
	}
	if (n_per_line != 0)
	  fprintf(grid_file,"\n");
      }
    break;
  }

  fclose(grid_file);

  return;
}

} // namespace

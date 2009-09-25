/*! \file
    \ingroup OEPROP
    \brief Enter brief description of file here
*/
#define EXTERN
#include "includes.h"
#include "prototypes.h"
#include "globals.h"

namespace {
  int* fock_to_pitzer(int, int*);
}

namespace psi { namespace oeprop {

void grid_parse()
{
  int i, errcod, atom;
  int mo, *fock_mos;
  char *tmpstring;
  double xmin, ymin, zmin, xmax, ymax, zmax;    /* molecular dimensions */

  errcod = ip_data("GRID","%d",&grid,0);
  if (grid < 0)
    punt("GRID type must be positive");

  if (grid == 5 || grid == 6)
    grid3d = 1;

  /* which MOs to plot? */
  if (grid == 5) {

    if (read_opdm)
      punt("Correlated WFN but asked to plot orbitals. Use WFN=SCF.");

    if (ip_exist("MO_TO_PLOT",0)) {
      errcod = ip_count("MO_TO_PLOT",&num_mos_to_plot,0);
      if (errcod == IPE_NOT_AN_ARRAY)
	punt("MO_TO_PLOT must be an array. See the manual.");

      errcod = ip_string("MO_TO_PLOT",&tmpstring,1,0);
      /* If signed integers are used - Fock ordering is used, convert indices back to Pitzer ... */
      if (tmpstring[0] == '+' || tmpstring[0] == '-') {
	fock_mos = init_int_array(num_mos_to_plot);
	for(mo=0;mo<num_mos_to_plot;mo++) {
	  errcod = ip_data("MO_TO_PLOT","%d",&fock_mos[mo],1,mo);
	  if (errcod != IPE_OK)
	    punt("Problem reading elements of MO_TO_PLOT. See the manual.");
	}
	mos_to_plot = fock_to_pitzer(num_mos_to_plot,fock_mos);
	free(fock_mos);
      }
      /* ... else indices are already Pitzer indices */
      else {
	mos_to_plot = init_int_array(num_mos_to_plot);
	for(mo=0;mo<num_mos_to_plot;mo++) {
	  errcod = ip_data("MO_TO_PLOT","%d",&mos_to_plot[mo],1,mo);
	  if (mos_to_plot[mo] <= 0 || mos_to_plot[mo] > nmo)
	    punt("One of the elements of MO_TO_PLOT out of range");
	  mos_to_plot[mo]--;
	  if (errcod != IPE_OK)
	    punt("Problem reading elements of MO_TO_PLOT. See the manual.");
	}
      }
      free(tmpstring);
    }
    /* If MO_TO_PLOT is not specified just use HOMO and LUMO */
    else {
      num_mos_to_plot = 2;
      fock_mos = init_int_array(num_mos_to_plot);
      fock_mos[0] = -1;
      fock_mos[1] = 1;
      mos_to_plot = fock_to_pitzer(num_mos_to_plot,fock_mos);
      free(fock_mos);
    }
  }

  /*--------------------------------------------------------------------
    If GRID_ORIGIN is present read in user's specification for the grid
    else create a default grid using molecular dimensions
   --------------------------------------------------------------------*/
  if (ip_exist("GRID_ORIGIN",0)) {

    ip_count("GRID_ORIGIN",&i,0);
    if (i != 3)
      punt("GRID_ORIGIN must have 3 components");
    for (i=0;i<3;i++) {
      errcod = ip_data("GRID_ORIGIN","%lf",&grid_origin[i],1,i);
      if (errcod != IPE_OK)
	punt("Error in the definition of GRID_ORIGIN");
    }

    if (ip_exist("GRID_UNIT_X",0)) {
      ip_count("GRID_UNIT_X",&i,0);
      if (i != 3)
	punt("GRID_UNIT_X must have 3 components");
      for (i=0;i<3;i++) {
	errcod = ip_data("GRID_UNIT_X","%lf",&grid_unit_x[i],1,i);
	if (errcod != IPE_OK)
	  punt("Problem in parsing GRID_UNIT_X");
      }
    }
    else
      punt("GRID_UNIT_X must be defined when GRID_ORIGIN is given");

    if (ip_exist("GRID_UNIT_Y",0)) {
      ip_count("GRID_UNIT_Y",&i,0);
      if (i != 3)
	punt("GRID_UNIT_Y must have 3 components");
      for (i=0;i<3;i++) {
	errcod = ip_data("GRID_UNIT_Y","%lf",&grid_unit_y[i],1,i);
	if (errcod != IPE_OK)
	  punt("Problem in parsing GRID_UNIT_Y");
      }
    }
    else
      punt("GRID_UNIT_Y must be defined when GRID_ORIGIN is given");

    if (grid3d == 0) {
      if (ip_exist("GRID_XY0",0)) {
	ip_count("GRID_XY0",&i,0);
	if (i != 2)
	  punt("GRID_XY0 must have 2 components");
	for (i=0;i<2;i++) {
	  errcod = ip_data("GRID_XY0","%lf",&grid_xy0[i],1,i);
	  if (errcod != IPE_OK)
	    punt("Error in the definition of GRID_XY0");
	}
      }
      else
	punt("GRID_XY0 is not defined");

      if (ip_exist("GRID_XY1",0)) {
	ip_count("GRID_XY1",&i,0);
	if (i != 2)
	  punt("GRID_XY1 must have 2 components");
	for (i=0;i<2;i++) {
	  errcod = ip_data("GRID_XY1","%lf",&grid_xy1[i],1,i);
	  if (errcod != IPE_OK)
	    punt("Error in the definition of GRID_XY1");
	  else if (grid_xy1[i] <= grid_xy0[i])
	    punt("GRID_XY1 must point to the upper right corner of the grid");
	}
      }
      else
	punt("GRID_XY1 is not defined");
    }
    else {
      if (ip_exist("GRID_XYZ0",0)) {
	ip_count("GRID_XYZ0",&i,0);
	if (i != 3)
	  punt("GRID_XYZ0 must have 3 components");
	for (i=0;i<3;i++) {
	  errcod = ip_data("GRID_XYZ0","%lf",&grid_xyz0[i],1,i);
	  if (errcod != IPE_OK)
	    punt("Error in the definition of GRID_XYZ0");
	}
      }
      else
	punt("GRID_XYZ0 is not defined");

      if (ip_exist("GRID_XYZ1",0)) {
	ip_count("GRID_XYZ1",&i,0);
	if (i != 3)
	  punt("GRID_XYZ1 must have 3 components");
	for (i=0;i<3;i++) {
	  errcod = ip_data("GRID_XYZ1","%lf",&grid_xyz1[i],1,i);
	  if (errcod != IPE_OK)
	    punt("Error in the definition of GRID_XYZ1");
	  else if (grid_xyz1[i] <= grid_xyz0[i])
	    punt("GRID_XYZ1 must point to the upper right corner of the grid parallelepiped");
	}
      }
      else
	punt("GRID_XYZ1 is not defined");

    }
  }
  /* GRID_ORIGIN not specified -- use molecular dimensions in 3D case */
  else if (grid3d) {

    xmax = ymax = zmax = 0.0;
    xmin = ymin = zmin = 0.0;

    for(atom=0; atom<natom; atom++) {
      if (geom[atom][0] > xmax)
	xmax = geom[atom][0];
      if (geom[atom][0] < xmin)
	xmin = geom[atom][0];
      if (geom[atom][1] > ymax)
	ymax = geom[atom][1];
      if (geom[atom][1] < ymin)
	ymin = geom[atom][1];
      if (geom[atom][2] > zmax)
	zmax = geom[atom][2];
      if (geom[atom][2] < zmin)
	zmin = geom[atom][2];
    }

    grid_origin[0] = 0.0;
    grid_origin[1] = 0.0;
    grid_origin[2] = 0.0;

    grid_unit_x[0] = 1.0;
    grid_unit_x[1] = 0.0;
    grid_unit_x[2] = 0.0;
    grid_unit_y[0] = 0.0;
    grid_unit_y[1] = 1.0;
    grid_unit_y[2] = 0.0;

    grid_xyz0[0] = xmin - 2.0;
    grid_xyz0[1] = ymin - 2.0;
    grid_xyz0[2] = zmin - 2.0;
    grid_xyz1[0] = xmax + 2.0;
    grid_xyz1[1] = ymax + 2.0;
    grid_xyz1[2] = zmax + 2.0;
  }
  /* It's a 2D grid and no GRID_ORIGIN is given -- fail */
  else
    punt("GRID_ORIGIN must be specified for two-dimentional grids");

  i = 0; errcod = ip_data("NIX","%d",&i,0);
  nix = i - 1;
  if (i <= 1)
    punt("NIX must be greater than 1");
  i = 0; errcod = ip_data("NIY","%d",&i,0);
  niy = i - 1;
  if (i <= 1)
    punt("NIY must be greater than 1");
  if (grid3d) {
    i = 0; errcod = ip_data("NIZ","%d",&i,0);
    niz = i - 1;
    if (i <= 1)
      punt("NIZ must be greater than 1");
  }

  /* orthonormalize grid unit vectors */
  grid_unitvec();

  /* Which format to use for the 3D grid */
  grid_format = NULL;
  errcod = ip_string("GRID_FORMAT", &grid_format, 0);
  if (grid_format == NULL) {
    if (grid3d)
      grid_format = strdup("GAUSSCUBE");
    else
      grid_format = strdup("PLOTMTV");
  }
  else if (strcmp(grid_format,"GAUSSCUBE") &&
	   strcmp(grid_format,"PLOTMTV") &&
	   strcmp(grid_format,"MEGAPOVPLUS"))
    punt("Invalid value for GRID_FORMAT");

  /* Check if the grid output format is compatibe with the type of grid */
  if (!strcmp(grid_format,"PLOTMTV") && grid3d)
    punt("GRID_FORMAT=PLOTMTV can only be used for 2-d grids");
  if (strcmp(grid_format,"PLOTMTV") && !grid3d)
    punt("Only GRID_FORMAT=PLOTMTV can be used for 2-d grids");

  if (grid3d == 0) {
    errcod = ip_data("GRID_ZMIN","%lf",&grid_zmin,0);
    errcod = ip_data("GRID_ZMAX","%lf",&grid_zmax,0);
    if (grid_zmin >= grid_zmax)
      punt("GRID_ZMIN must be less than GRID_ZMAX");
    errcod = ip_data("EDGRAD_LOGSCALE","%d",&edgrad_logscale,0);
  }
}

}} //namespace psi::oeprop
namespace {
  using namespace psi::oeprop;

int* fock_to_pitzer(int num_mo, int* fock_mos)
{
  int mo, irrep, o, v;
  int nvirt, nocc, ndocc, nsocc;
  int nv, no, na;
  int current_highest, current_lowest;
  int *f2p_virt, *f2p_occ;
  int *mo_offset, *occ_offset, *virt_offset;
  int *pitzer_mos;
  double eval;

  ndocc = nvirt = nsocc = 0;
  for(irrep=0;irrep<nirreps;irrep++) {
    ndocc += clsdpi[irrep];
    nsocc += openpi[irrep];
    nvirt += orbspi[irrep] - clsdpi[irrep] - openpi[irrep];
  }
  nocc = ndocc + nsocc;
  f2p_occ = init_int_array(nocc);
  f2p_virt = init_int_array(nvirt);
  pitzer_mos = init_int_array(num_mo);
  virt_offset = init_int_array(nirreps);
  occ_offset = init_int_array(nirreps);
  mo_offset = init_int_array(nirreps);

  virt_offset[0] = clsdpi[0] + openpi[0];
  occ_offset[0] = clsdpi[0] + openpi[0] - 1;
  mo_offset[0] = 0;
  for(irrep=1;irrep<nirreps;irrep++) {
    na = orbspi[irrep];
    mo_offset[irrep] = mo_offset[irrep-1] + orbspi[irrep-1];
    no = clsdpi[irrep] + openpi[irrep];
    nv = na - no;
    if (no)
      occ_offset[irrep] = mo_offset[irrep] + no - 1;
    else
      occ_offset[irrep] = -1;
    if(nv)
      virt_offset[irrep] = mo_offset[irrep] + no;
    else
      virt_offset[irrep] = -1;
  }


  /*---------------------------------------------------------------
    construct Fock->Pitzer mapping for virtual orbitals (f2p_virt)
   ---------------------------------------------------------------*/
  for(v=0;v<nvirt;v++) {
    eval = 1E300;    /* to be greater than any eigenvalue to be encountered */
    for(irrep=0;irrep<nirreps;irrep++) {
      if (virt_offset[irrep] >= 0)
	if (eval > scf_evals[virt_offset[irrep]]) {
	  current_lowest = irrep;
	  eval = scf_evals[virt_offset[irrep]];
	}
    }
    f2p_virt[v] = virt_offset[current_lowest];
    virt_offset[current_lowest]++;
    /* If all virtual MOs in this irrep have been exhausted - set virt_offset for
       that irrep to -1 so that it's skipped on subsequent passes */
    if (current_lowest == nirreps-1) {
      if (virt_offset[current_lowest] == nmo)
	virt_offset[current_lowest] = -1;
    }
    else {
      if (virt_offset[current_lowest] == mo_offset[current_lowest+1])
	virt_offset[current_lowest] = -1;
    }
  }

  /*----------------------------------------------------------------
    construct Fock->Pitzer mapping for occupied orbitals (f2o_occ)
   ----------------------------------------------------------------*/
  for(o=0;o<nocc;o++) {
    eval = -1E300;    /* to be smaller than any eigenvalue to be encountered */
    for(irrep=0;irrep<nirreps;irrep++) {
      if (occ_offset[irrep] >= 0)
	if (eval < scf_evals[occ_offset[irrep]]) {
	  current_highest = irrep;
	  eval = scf_evals[occ_offset[irrep]];
	}
    }
    f2p_occ[o] = occ_offset[current_highest];
    occ_offset[current_highest]--;
    /* If all occupied MOs in this irrep have been exhausted - set occ_offset for
       that irrep to -1 so that it's skipped on subsequent passes */
    if (occ_offset[current_highest] == mo_offset[current_highest]-1)
      occ_offset[current_highest] = -1;
  }

  for(mo=0;mo<num_mo;mo++)
    if (fock_mos[mo] < 0) {
      o = -fock_mos[mo] - 1;
      if (o >= nocc)
	punt("MOS_TO_PLOT contains occupied indices outside of allowed range");
      else
	pitzer_mos[mo] = f2p_occ[o];
    }
    else {
      v = fock_mos[mo] - 1;
      if (v >= nvirt)
	punt("MOS_TO_PLOT contains unoccupied indices outside of allowed range");
      else
	pitzer_mos[mo] = f2p_virt[v];
    }

  free(virt_offset);
  free(occ_offset);
  free(mo_offset);
  free(f2p_occ);
  free(f2p_virt);

  return pitzer_mos;
}

} // namespace

namespace psi { namespace oeprop {

void grid_unitvec()
{
  int i;
  double sum1, sum2, dot;
  double step_x, step_y, step_z;

  /*-------------------------------------------------
    First orthonormalize grid_unit_x and grid_unit_y
    and compute grid_unit_z
   -------------------------------------------------*/

  /* Normalizing the grid_unit's */
  dot_arr(grid_unit_x,grid_unit_x,3,&sum1);
  dot_arr(grid_unit_y,grid_unit_y,3,&sum2);
  sum1 = sqrt(sum1); sum2 = sqrt(sum2);
  for(i=0;i<3;i++) {
    grid_unit_x[i] /= sum1;
    grid_unit_y[i] /= sum2;
  }

  /* Checking if vectors are parallel */
  dot_arr(grid_unit_x,grid_unit_y,3,&dot);
  if (1.0 - fabs(dot) < ADOTB_ORTHOGONAL) /* Vectors are parallel - abort */
    punt("Vectors GRID_UNIT_X and GRID_UNIT_Y are parallel");
  else if (fabs(dot) > ADOTB_ORTHOGONAL) {  /* Vectors are not orthogonal - orthonormalizing */
    for(i=0;i<3;i++)
      grid_unit_y[i] -= dot*grid_unit_x[i];
    dot_arr(grid_unit_y,grid_unit_y,3,&sum1);
    sum1 = sqrt(sum1);
    for(i=0;i<3;i++)
      grid_unit_y[i] /= sum1;
  }

  /* Get grid_unit_z as a vector product */
  if (grid3d) {
    grid_unit_z[0] = grid_unit_x[1]*grid_unit_y[2] - grid_unit_x[2]*grid_unit_y[1];
    grid_unit_z[1] = grid_unit_x[2]*grid_unit_y[0] - grid_unit_x[0]*grid_unit_y[2];
    grid_unit_z[2] = grid_unit_x[0]*grid_unit_y[1] - grid_unit_x[1]*grid_unit_y[0];
  }

  /* grid_origin will now contain the origin of the grid rectangle/box rather
     than the grid coordinate systems */
  if (grid3d)
    for(i=0;i<3;i++)
      grid_origin[i] += grid_xyz0[0]*grid_unit_x[i] + grid_xyz0[1]*grid_unit_y[i] + grid_xyz0[2]*grid_unit_z[i];
  else
    for(i=0;i<3;i++)
      grid_origin[i] += grid_xy0[0]*grid_unit_x[i] + grid_xy0[1]*grid_unit_y[i];


  /* Compute grid size and unit cell vectors */
  if (grid3d) {
    step_x = (grid_xyz1[0] - grid_xyz0[0])/nix;
    step_y = (grid_xyz1[1] - grid_xyz0[1])/niy;
    step_z = (grid_xyz1[2] - grid_xyz0[2])/niz;
    for(i=0;i<3;i++) {
      grid_step_x[i] = grid_unit_x[i]*step_x;
      grid_step_y[i] = grid_unit_y[i]*step_y;
      grid_step_z[i] = grid_unit_z[i]*step_z;
    }
  }
  else {
    step_x = (grid_xy1[0] - grid_xy0[0])/nix;
    step_y = (grid_xy1[1] - grid_xy0[1])/niy;
    for(i=0;i<3;i++) {
      grid_step_x[i] = grid_unit_x[i]*step_x;
      grid_step_y[i] = grid_unit_y[i]*step_y;
    }
  }

  return;
}

}} // namespace psi::oeprop

/*! \file
    \ingroup OEPROP
    \brief Enter brief description of file here
*/
#define EXTERN
#include "includes.h"
#include "prototypes.h"
#include "globals.h"

namespace psi { namespace oeprop {

int* fock_to_pitzer(int, int*);
void grid_unitvec(void);

void grid_parse(Options & options)
{
  int i, errcod, atom;
  int mo, *fock_mos;
  std::string tmpstring;
  double xmin, ymin, zmin, xmax, ymax, zmax;    /* molecular dimensions */
  double *v3;

  grid = options.get_int("GRID");
  if (grid < 0)
    throw PsiException("GRID type must be positive", __FILE__, __LINE__);

  if (grid == 5 || grid == 6)
    grid3d = 1;

  /* which MOs to plot? */
  if (grid == 5) {

    if (read_opdm)
      throw PsiException("Correlated WFN but asked to plot orbitals. Use WFN=SCF.", __FILE__, __LINE__);

    if(options["MO_TO_PLOT"].has_changed()) {
      num_mos_to_plot = options["MO_TO_PLOT"].size();

      tmpstring = options.get_str("MO_TO_PLOT");

      /* If signed integers are used - Fock ordering is used, convert indices back to Pitzer ... */
      if (tmpstring[0] == '+' || tmpstring[0] == '-') {
        fock_mos = options.get_int_array("MO_TO_PLOT");
        mos_to_plot = fock_to_pitzer(num_mos_to_plot,fock_mos);
        delete [] fock_mos;
      }
      /* ... else indices are already Pitzer indices */
      else {
        mos_to_plot = options.get_int_array("MO_TO_PLOT");
        for(mo=0;mo<num_mos_to_plot;mo++) {
          if (mos_to_plot[mo] <= 0 || mos_to_plot[mo] > nmo)
            throw PsiException("One of the elements of MO_TO_PLOT out of range", __FILE__, __LINE__);
          mos_to_plot[mo]--;
        }
      }
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
  if(options["GRID_ORIGIN"].has_changed()) {

    if (options["GRID_ORIGIN"].size() != 3)
      throw PsiException("GRID_ORIGIN must have 3 components", __FILE__, __LINE__);
    v3 = options.get_double_array("GRID_ORIGIN");
    memcpy(grid_origin, v3, 3*sizeof(double));
    delete [] v3;

    if(options["GRID_UNIT_X"].has_changed()) {
      if (options["GRID_UNIT_X"].size() != 3)
        throw PsiException("GRID_UNIT_X must have 3 components", __FILE__, __LINE__);
      v3 = options.get_double_array("GRID_UNIT_X");
      memcpy(grid_unit_x, v3, 3*sizeof(double));
      delete [] v3;
    }
    else
      throw PsiException("GRID_UNIT_X must be defined when GRID_ORIGIN is given", __FILE__, __LINE__);

    if(options["GRID_UNIT_Y"].has_changed()) {
      if (options["GRID_UNIT_Y"].size() != 3)
        throw PsiException("GRID_UNIT_Y must have 3 components", __FILE__, __LINE__);
      v3 = options.get_double_array("GRID_UNIT_Y");
      memcpy(grid_unit_y, v3, 3*sizeof(double));
      delete [] v3;
    }
    else
      throw PsiException("GRID_UNIT_Y must be defined when GRID_ORIGIN is given", __FILE__, __LINE__);

    if (grid3d == 0) {
      if(options["GRID_XY0"].has_changed()) {
        if (options["GRID_XY0"].size() != 2)
          throw PsiException("GRID_XY0 must have 2 components", __FILE__, __LINE__);
        v3 = options.get_double_array("GRID_XY0");
        memcpy(grid_xy0, v3, 2*sizeof(double));
        delete [] v3;
      }
      else throw PsiException("GRID_XY0 is not defined", __FILE__, __LINE__);
  
      if(options["GRID_XY1"].has_changed()) {
        if (options["GRID_XY1"].size() != 2)
          throw PsiException("GRID_XY1 must have 2 components", __FILE__, __LINE__);
        v3 = options.get_double_array("GRID_XY1");
        memcpy(grid_xy1, v3, 2*sizeof(double));
        delete [] v3;
      }
      else throw PsiException("GRID_XY1 is not defined", __FILE__, __LINE__);
  
      for (i=0; i<2; ++i)
	    if (grid_xy1[i] <= grid_xy0[i])
          throw("GRID_XY1 must point to the upper right corner of the grid");
    }
    else { // grid3d
      if(options["GRID_XYZ0"].has_changed()) {
        if (options["GRID_XYZ0"].size() != 3)
          throw PsiException("GRID_XYZ0 must have 3 components", __FILE__, __LINE__);
        v3 = options.get_double_array("GRID_XYZ0");
        memcpy(grid_xyz0, v3, 3*sizeof(double));
        delete [] v3;
      }
      else throw PsiException("GRID_XYZ0 is not defined", __FILE__, __LINE__);
  
      if(options["GRID_XYZ1"].has_changed()) {
        if (options["GRID_XYZ1"].size()  != 3)
          throw PsiException("GRID_XYZ1 must have 3 components", __FILE__, __LINE__);
        v3 = options.get_double_array("GRID_XYZ1");
        memcpy(grid_xyz1, v3, 3*sizeof(double));
        delete [] v3;
      }
      else throw PsiException("GRID_XYZ1 is not defined", __FILE__, __LINE__);

      for (i=0; i<3; ++i)
	    if (grid_xyz1[i] <= grid_xyz0[i])
	      throw PsiException("GRID_XYZ1 must point to the upper right corner of the grid parallelepiped",
            __FILE__, __LINE__);
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
    throw PsiException("GRID_ORIGIN must be specified for two-dimentional grids", __FILE__, __LINE__);

  i = options.get_int("NIX"); // default 0
  nix = i - 1;
  if (i <= 1)
    throw PsiException("NIX must be greater than 1", __FILE__, __LINE__);

  i = options.get_int("NIY");
  niy = i - 1;
  if (i <= 1)
    throw PsiException("NIY must be greater than 1", __FILE__, __LINE__);
  if (grid3d) {
    i = options.get_int("NIZ");
    niz = i - 1;
    if (i <= 1)
      throw PsiException("NIZ must be greater than 1", __FILE__, __LINE__);
  }

  /* orthonormalize grid unit vectors */
  grid_unitvec();

  /* Which format to use for the 3D grid */
  grid_format = "";
  if(options["GRID_FORMAT"].has_changed())
    grid_format = options.get_str("GRID_FORMAT");
  else {
    if (grid3d)
      grid_format = "GAUSSCUBE";
    else
      grid_format = "PLOTMTV";
  }

  if (grid_format != "GAUSSCUBE" &&
      grid_format != "PLOTMTV" &&
      grid_format != "MEGAPOVPLUS")
    throw PsiException("Invalid value for GRID_FORMAT", __FILE__, __LINE__);

  /* Check if the grid output format is compatibe with the type of grid */
  if(grid_format == "PLOTMTV" && grid3d)
    throw PsiException("GRID_FORMAT=PLOTMTV can only be used for 2-d grids", __FILE__, __LINE__);
  if(grid_format != "PLOTMTV" && !grid3d)
    throw PsiException("Only GRID_FORMAT=PLOTMTV can be used for 2-d grids", __FILE__, __LINE__);

  if (grid3d == 0) {
    grid_zmin = options.get_double("GRID_ZMIN");
    grid_zmax = options.get_double("GRID_ZMAX");
    if (grid_zmin >= grid_zmax)
      throw PsiException("GRID_ZMIN must be less than GRID_ZMAX", __FILE__, __LINE__);
    edgrad_logscale = options.get_int("EDGRAD_LOGSCALE");
  }
}

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
	throw("MOS_TO_PLOT contains occupied indices outside of allowed range");
      else
	pitzer_mos[mo] = f2p_occ[o];
    }
    else {
      v = fock_mos[mo] - 1;
      if (v >= nvirt)
	throw("MOS_TO_PLOT contains unoccupied indices outside of allowed range");
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

void grid_unitvec(void)
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
    throw PsiException("Vectors GRID_UNIT_X and GRID_UNIT_Y are parallel", __FILE__, __LINE__);

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

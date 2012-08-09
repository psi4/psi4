/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <stdio.h>
#include <math.h>
#include "psi4-dec.h"
#include <libmints/mints.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libiwl/iwl.h>
#include <libchkpt/chkpt.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include <psifiles.h>
#include <physconst.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

/*
** densgrid_RHF(): Compute the values of the ground-state one-particle 
** density at a set of grid points.
**
** TDC, 7/2012
*/
void compute_delta(double **delta, double x, double y, double z);
int nmo, nso, nao; // global
double **scf, **u; // global
boost::shared_ptr<Molecule> molecule;
boost::shared_ptr<BasisSet> basis;
boost::shared_ptr<Wavefunction> wfn;

void densgrid_RHF(Options& options)
{
  double dens;
  double **D, **delta;
  double x, y, z;
  double xmin, xmax, ymin, ymax, zmin, zmax;
  double xstep, ystep, zstep;
  int *order;
  double **scf_pitzer;

  wfn = Process::environment.wavefunction();
  molecule = wfn->molecule();
  basis = wfn->basisset();

  nao = basis->nao();
  nso = moinfo.nso;
  nmo = moinfo.nmo;
  chkpt_init(PSIO_OPEN_OLD);
  scf_pitzer = chkpt_rd_scf();
  chkpt_close();

  D = moinfo.opdm; // A block matrix
  delta = block_matrix(nmo, nmo); // Dirac delta function 

  // Set up AO->SO transformation matrix (u)
  MintsHelper helper(options, 0);
  SharedMatrix aotoso = helper.petite_list(true)->aotoso();
  int *col_offset = new int[wfn->nirrep()];
  col_offset[0] = 0;
  for(int h=1; h < wfn->nirrep(); h++)
    col_offset[h] = col_offset[h-1] + aotoso->coldim(h-1);

  u = block_matrix(nao, nso);
  for(int h=0; h < wfn->nirrep(); h++)
    for(int j=0; j < aotoso->coldim(h); j++)
      for(int i=0; i < nao; i++)
        u[i][j+col_offset[h]] = aotoso->get(h, i, j);
  delete[] col_offset;

  /*** Arrange the SCF eigenvectors into QT ordering ***/
  order = moinfo.pitzer2qt;
  scf = block_matrix(nso, nmo);
  for(int i=0; i < nmo; i++) {
      int I = order[i];  /* Pitzer --> QT */
      for(int j=0; j < nso; j++) scf[j][I] = scf_pitzer[j][i];
    }

  // Scan along Cartesian axes to determine dimensions of box 
  molecule->print();
  fprintf(outfile, "  Grid domain:\n");
  xmin = xmax = molecule->xyz(0, 0);
  ymin = ymax = molecule->xyz(0, 1);
  zmin = zmax = molecule->xyz(0, 2);
  for(int atom=1; atom < molecule->natom(); atom++) {
    if(molecule->xyz(atom, 0) < xmin) xmin = molecule->xyz(atom, 0);
    if(molecule->xyz(atom, 1) < ymin) ymin = molecule->xyz(atom, 1);
    if(molecule->xyz(atom, 2) < zmin) zmin = molecule->xyz(atom, 2);
    if(molecule->xyz(atom, 0) > xmax) xmax = molecule->xyz(atom, 0);
    if(molecule->xyz(atom, 1) > ymax) ymax = molecule->xyz(atom, 1);
    if(molecule->xyz(atom, 2) > zmax) zmax = molecule->xyz(atom, 2);
  }

  xmin *= _bohr2angstroms;
  xmax *= _bohr2angstroms;
  ymin *= _bohr2angstroms;
  ymax *= _bohr2angstroms;
  zmin *= _bohr2angstroms;
  zmax *= _bohr2angstroms;

  double b2a3 = _bohr2angstroms * _bohr2angstroms * _bohr2angstroms;

  do {
    xmin -= 0.1;
    compute_delta(delta, xmin/_bohr2angstroms, 0, 0);
    dens = 0.0;
    for(int i=0; i < nmo; i++)
      for(int j=0; j < nmo; j++)
        dens += delta[i][j] * D[i][j];
  } while((dens/b2a3) > options.get_double("ONEPDM_GRID_CUTOFF"));
  fprintf(outfile, "  xmin = %8.6f (Angstrom);  density(xmin,0,0) (e/Ang^3) = %8.6e\n", xmin, dens);

  do {
    xmax += 0.1;
    compute_delta(delta, xmax/_bohr2angstroms, 0, 0);
    dens = 0.0;
    for(int i=0; i < nmo; i++)
      for(int j=0; j < nmo; j++)
        dens += delta[i][j] * D[i][j];
  } while((dens/b2a3) > options.get_double("ONEPDM_GRID_CUTOFF"));
  fprintf(outfile, "  xmax = %8.6f (Angstrom);   density(xmax,0,0) (e/Ang^3) = %8.6e\n", xmax, dens);

  do {
    ymin -= 0.1;
    compute_delta(delta, 0, ymin/_bohr2angstroms, 0);
    dens = 0.0;
    for(int i=0; i < nmo; i++)
      for(int j=0; j < nmo; j++)
        dens += delta[i][j] * D[i][j];
  } while((dens/b2a3) > options.get_double("ONEPDM_GRID_CUTOFF"));
  fprintf(outfile, "  ymin = %8.6f (Angstrom);  density(0,ymin,0) (e/Ang^3) = %8.6e\n", ymin, dens);

  do {
    ymax += 0.1;
    compute_delta(delta, 0, ymax/_bohr2angstroms, 0);
    dens = 0.0;
    for(int i=0; i < nmo; i++)
      for(int j=0; j < nmo; j++)
        dens += delta[i][j] * D[i][j];
  } while((dens/b2a3) > options.get_double("ONEPDM_GRID_CUTOFF"));
  fprintf(outfile, "  ymax = %8.6f (Angstrom);   density(0,ymax,0) (e/Ang^3) = %8.6e\n", ymax, dens);

  do {
    zmin -= 0.1;
    compute_delta(delta, 0, 0, zmin/_bohr2angstroms);
    dens = 0.0;
    for(int i=0; i < nmo; i++)
      for(int j=0; j < nmo; j++)
        dens += delta[i][j] * D[i][j];
  } while((dens/b2a3) > options.get_double("ONEPDM_GRID_CUTOFF"));
  fprintf(outfile, "  zmin = %8.6f (Angstrom);  density(0,0,zmin) (e/Ang^3) = %8.6e\n", zmin, dens);

  do {
    zmax += 0.1;
    compute_delta(delta, 0, 0, zmax/_bohr2angstroms);
    dens = 0.0;
    for(int i=0; i < nmo; i++)
      for(int j=0; j < nmo; j++)
        dens += delta[i][j] * D[i][j];
  } while((dens/b2a3) > options.get_double("ONEPDM_GRID_CUTOFF"));
  fprintf(outfile, "  zmax = %8.6f (Angstrom);   density(0,0,zmax) (e/Ang^3) = %8.6e\n", zmax, dens);

  // Compute density at the nuclei
  fprintf(outfile, "  Density at nuclei:\n");
  for(int atom=0; atom < molecule->natom(); atom++) {
    x = molecule->xyz(atom, 0);
    y = molecule->xyz(atom, 1);
    z = molecule->xyz(atom, 2);
    compute_delta(delta, x, y, z);
    dens = 0.0;
    for(int i=0; i < nmo; i++)
      for(int j=0; j < nmo; j++)
        dens += delta[i][j] * D[i][j];

    fprintf(outfile, "  Atom %d (%8.6f, %8.6f, %8.5f), dens = %20.12f (e/Ang^3)\n", atom, x*_bohr2angstroms,
y*_bohr2angstroms, z*_bohr2angstroms, dens/b2a3);
   }

  double step_size = options.get_double("ONEPDM_GRID_STEPSIZE");
  int xsteps = (int) ((xmax - xmin)/step_size + 1);
  int ysteps = (int) ((ymax - ymin)/step_size + 1);
  int zsteps = (int) ((zmax - zmin)/step_size + 1);

  // Prep .dx file
  FILE *dxfile;
  ffile(&dxfile, "density.dx", 0);
  fprintf(dxfile, "#  Output from PSI4 calculation\n");
  fprintf(dxfile, "#  Electronic density (in e/ang^3) for: \n");
  fprintf(dxfile, "object 1 class gridpositions counts %d %d %d\n", xsteps, ysteps, zsteps);
  fprintf(dxfile, "origin %8.6E  %8.6E  %8.6E\n", 0.0, 0.0, 0.0);
  fprintf(dxfile, "delta %8.6E  %8.6E  %8.6E\n", step_size, 0.0, 0.0);
  fprintf(dxfile, "delta %8.6E  %8.6E  %8.6E\n", 0.0, step_size, 0.0);
  fprintf(dxfile, "delta %8.6E  %8.6E  %8.6E\n", 0.0, 0.0, step_size);
  fprintf(dxfile, "object 1 class gridpositions counts %d %d %d\n", xsteps, ysteps, zsteps);
  fprintf(dxfile, "object 3 class array double rank 0 items %d data follows\n", xsteps*ysteps*zsteps);

  // Loop over points and integrate along the way
  double charge = 0;
  int count=0;
  for(x=xmin; x <= xmax; x += step_size) {
    for(y=ymin; y <= ymax; y += step_size) {
      for(z=zmin; z <= zmax; z += step_size) {

        // Compute delta function in Gaussian basis
        compute_delta(delta, x/_bohr2angstroms, y/_bohr2angstroms, z/_bohr2angstroms);

        dens = 0.0; // e/bohr^3
        for(int i=0; i < nmo; i++)
          for(int j=0; j < nmo; j++)
            dens += delta[i][j] * D[i][j];

        dens /= b2a3; // convert to e/Ang^3

        fprintf(dxfile, "  %8.6E", dens);
        count++;
        if(count % 3 == 0) fprintf(dxfile, "\n");

        charge += dens * step_size * step_size * step_size;

      } // z
    }  // y
  } // x

  fprintf(outfile, "    Number of electrons = %20.12f?\n", charge);

  if(count % 3 != 0) fprintf(dxfile, "\n");
  fprintf(dxfile, "attribute \"dep\" string \"positions\"\n");
  fprintf(dxfile, "object \"regular positions regular connections\" class field\n");
  fprintf(dxfile, "component \"positions\" value 1\n");
  fprintf(dxfile, "component \"connections\" value 2\n");
  fprintf(dxfile, "component \"data\" value 3\n");
  fprintf(dxfile, "\n");
  fprintf(dxfile, "end");
  fclose(dxfile);

  ffile(&dxfile, "molecule.dx", 0);
  fprintf(dxfile, "%d\n", molecule->natom());
  fprintf(dxfile, "Initial atomic coordinates\n");
  for(int i=0; i < molecule->natom(); i++) {
    fprintf(dxfile, "%2s  ", molecule->symbol(i).c_str());
    fprintf(dxfile, "  %9.6f  %9.6f  %9.6f\n", molecule->x(i)*_bohr2angstroms, molecule->y(i)*_bohr2angstroms, molecule->z(i)*_bohr2angstroms);
  }
  fflush(dxfile);
  fclose(dxfile);

  free_block(delta);
  free_block(scf);
}

void compute_delta(double **delta, double x, double y, double z)
{
  int i, j;
  double *phi_ao, *phi_so, *phi_mo;

  phi_ao = init_array(nao);  /* AO function values */
  phi_so = init_array(nso);  /* SO function values */
  phi_mo = init_array(nmo);  /* MO function values */

  basis->compute_phi(phi_ao, x, y, z);

  /*  for(i=0; i < nao; i++) printf("%d %20.10f\n", i, phi_ao[i]); */

  /* Transform the basis function values to the MO basis */
  C_DGEMV('t', nao, nso, 1.0, u[0], nso, phi_ao, 1, 0.0, phi_so, 1);

  C_DGEMV('t', nmo, nso, 1.0, scf[0], nmo, phi_so, 1, 0.0, phi_mo, 1);

  /* for(i=0; i < nmo; i++) printf("%d %20.10f\n", i, phi_mo[i]); */


  /* Build the MO-basis delta function */
  for(i=0; i < nmo; i++)
    for(j=0; j < nmo; j++)
      delta[i][j] = phi_mo[i] * phi_mo[j];

  free(phi_ao);
  free(phi_so);
  free(phi_mo);
}

}} // namespace psi::ccdensity

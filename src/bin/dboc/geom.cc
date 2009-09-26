/*! \file
    \ingroup DBOC
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include "molecule.h"
#include "params.h"

namespace {
  void append_geom(FILE *geomdat, double **geom, int disp);
}

namespace psi { namespace DBOC {

extern Molecule_t Molecule;
extern Params_t Params;

void setup_geoms()
{
  const int disp_per_coord = Params.disp_per_coord;
  const double delta = Params.delta;
  int d;
  int coord, disp, atom, xyz;
  FILE *geometry;

  /*--- Open geom.dat for writing ---*/
  ffile(&geometry, "geom.dat", 0);
  fprintf(geometry,"input:(\n");

  /*--- make a local copy of reference geometry ---*/
  double **geom_copy = block_matrix(Molecule.natom,3);
  for(atom=0; atom<Molecule.natom; atom++)
    for(xyz=0; xyz<3; xyz++)
      geom_copy[atom][xyz] = Molecule.geom[atom][xyz];

  //
  // Setup Params.coords if user didn't specify displacements
  //
  if (Params.coords == NULL) {
    Params.ncoord = Molecule.natom * 3;
    Params.coords = new Params_t::Coord_t[Params.ncoord];
    for(int coord=0; coord<Params.ncoord; coord++) {
      Params.coords[coord].index = coord;
      Params.coords[coord].coeff = 1.0;
    }
  }

  for(coord=0,disp=1; coord<Params.ncoord; coord++) {
    Params_t::Coord_t* c = &(Params.coords[coord]);
    int atom = c->index/3;
    int xyz = c->index%3;
    for(d=-1; d<=1; d+=2,disp++) {
      geom_copy[atom][xyz] += d*delta;
      append_geom(geometry,geom_copy,disp);
      geom_copy[atom][xyz] -= d*delta;
    }
    if (Params.disp_per_coord == 4)
      for(d=-2; d<=2; d+=4,disp++) {
        geom_copy[atom][xyz] += d*delta;
        append_geom(geometry,geom_copy,disp);
        geom_copy[atom][xyz] -= d*delta;
      }
  }
  
  free_block(geom_copy);
  fprintf(geometry,")\n");
  fclose(geometry);

}

}} // namespace psi::DBOC

namespace {
  using namespace psi::DBOC;
  void append_geom(FILE *geomdat, double **geom, int disp)
  {
    fprintf(geomdat,"%% DBOC cartesian displacement %d\n",disp);
    fprintf(geomdat,"geometry%d = (\n",disp);
    for(int atom=0; atom<Molecule.natom; atom++)
      fprintf(geomdat,"  (%4.2lf %20.15lf %20.15lf %20.15lf)\n",
              Molecule.zvals[atom],geom[atom][0],geom[atom][1],geom[atom][2]);
    fprintf(geomdat,")\n",disp);
  }
}


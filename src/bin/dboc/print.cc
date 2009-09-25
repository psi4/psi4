/*! \file
    \ingroup DBOC
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstring>
#include "params.h"
#include "molecule.h"
#include "float.h"

#define PRINT_INTRO 1
#define PRINT_PARAMS 1
#define PRINT_DISP 2

namespace psi { namespace dboc {

extern Params_t Params;
extern Molecule_t Molecule;
extern "C" FILE *outfile;

void print_intro()
{
  if (Params.print_lvl >= PrintLevels::print_intro) {
     fprintf(outfile,"                  ----------------------------------------------\n");
     fprintf(outfile,"                    DBOC: diagonal Born-Oppenheimer correction\n");
     fprintf(outfile,"                      evaluation by numerical differentiation\n");
     fprintf(outfile,"                  ----------------------------------------------\n\n");
  }
  fflush(outfile);
}

void print_params()
{
  if (Params.print_lvl >= PrintLevels::print_params) {
    fprintf(outfile,"\n  -OPTIONS:\n");
    fprintf(outfile,"    Label                       = %s\n",Params.label);
    fprintf(outfile,"    # of disp. per coord        = %d\n",Params.disp_per_coord);
    fprintf(outfile,"    Displacement size           = %lf a.u.\n",Params.delta);
    if (strcmp(Params.wfn,"SCF"))
      fprintf(outfile,"    Wave function               = %s\n",Params.wfn);
    else if (Params.reftype == Params_t::rhf)
      fprintf(outfile,"    Wave function               = RHF SCF\n");
    else if (Params.reftype == Params_t::rohf)
      fprintf(outfile,"    Wave function               = ROHF SCF\n");
    else if (Params.reftype == Params_t::uhf)
      fprintf(outfile,"    Wave function               = UHF SCF\n");
    fprintf(outfile,"    Print level                 = %d\n",Params.print_lvl);
    fprintf(outfile,"    Memory                      = %ld MB\n",Params.max_memory/(1L << 20));
    fprintf(outfile,"    Number of threads           = %d\n",Params.num_threads);
    fprintf(outfile,"    sizeof(double)              = %d\n",sizeof(double));
    fprintf(outfile,"    sizeof(long double)         = %d\n",sizeof(long double));
    fprintf(outfile,"    sizeof(FLOAT)               = %d\n",sizeof(FLOAT));
    fprintf(outfile,"\n");

    fprintf(outfile,"    Cartesian displacements:\n");
    for(int coord=0; coord<Params.ncoord; coord++) {
      int atom = Params.coords[coord].index/3;
      int xyz = Params.coords[coord].index%3;
      char dir;
      switch (xyz) {
      case 0:  dir = 'x'; break;
      case 1:  dir = 'y'; break;
      case 2:  dir = 'z'; break;
      }
      fprintf(outfile,"      atom %d  %c  %7s   degen %5.1lf\n",
	      atom+1, dir, (Params.coords[coord].symm ? "symm" : "nonsymm"),
	      Params.coords[coord].coeff);
      fprintf(outfile,"\n");
    }

    if (Params.isotopes) {
      fprintf(outfile,"    User-specified isotope labels:\n");
      for(int atom=0; atom<Params.nisotope; atom++)
	fprintf(outfile,"      %s\n",Params.isotopes[atom]);
      fprintf(outfile,"\n");
    }
    else
      fprintf(outfile,"    By default, using the most abundant isotope for each element\n\n");
  }

  fflush(outfile);
}
  
void print_geom()
{
  if (Params.print_lvl >= PrintLevels::print_intro) {
    fprintf(outfile, "  -Reference Geometry:\n");
    for(int i=0; i < Molecule.natom; i++) {
      fprintf(outfile, "\n   %1.0f ", Molecule.zvals[i]);
      for(int j=0; j < 3; j++)
	fprintf(outfile, "%20.10f  ", Molecule.geom[i][j]);
    }
    fprintf(outfile, "\n\n");
  }
}

}} // namespace psi::dboc

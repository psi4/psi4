/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680


  This was written by Judith Hill.
*/


#include <mra/mra.h>

void doplotptk(World& world, int step, const functionT& psi, double Lplot, long numpt, const char* fname) {
    double start = wall_time();
    Tensor<double> cell(3,2);
    std::vector<long> npt(3, numpt);
    cell(_,0) = -Lplot;
    cell(_,1) =  Lplot;
    plotvtk(psi, fname, cell, npt, false);
    if (world.rank() == 0) print("plotting used", wall_time()-start);
}

template <typename T, int NDIM>
void plotvtk(const Function<T,NDIM>& function,
            const char* filename,
            const Tensor<double>& cell,
            const std::vector<long>& npt,
            bool binary) {
    PROFILE_FUNC;
    MADNESS_ASSERT(NDIM<=6);

    function.verify();
    World& world = const_cast< Function<T,NDIM>& >(function).world();
    FILE *f=0;
    if (world.rank() == 0) {
        f = fopen(filename, "w");
        if (!f) MADNESS_EXCEPTION("plotvtk: failed to open the plot file", 0);

        double spacex = (cell(0,1) - cell(0,0))/(npt[0] - 1) ; 
        double spacey = (cell(1,1) - cell(1,0))/(npt[1] - 1) ; 
        double spacez = (cell(2,1) - cell(2,0))/(npt[2] - 1) ; 
  
        fprintf(f, "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n") ;
        fprintf(f, "  <StructuredGrid WholeExtent=\"0 %ld 0 %ld 0 %ld\">\n",
                       npt[0]-1, npt[1]-1, npt[2]-1) ;
        fprintf(f, "  <Piece Extent=\"0 %ld 0 %ld 0 %ld\">\n",
                       npt[0]-1, npt[1]-1, npt[2]-1) ;
        fprintf(f, "      <PointData>\n") ;
        fprintf(f, "        <DataArray Name=\"data\" format=\"ascii\" type=\"Float32\" NumberOfComponents=\"1\">\n") ;

        world.gop.fence();
        Tensor<T> tmpr = function.eval_cube(cell, npt);
            for (IndexIterator it(npt); it; ++it) {
                fprintf (f, "%.6e\n", tmpr(*it)) ;
            }
        fprintf(f, "        </DataArray>\n") ;
        fprintf(f, "      </PointData>\n") ;

        fprintf(f, "      <Points>\n") ;
        fprintf(f, "        <DataArray NumberOfComponents=\"3\" type=\"Float32\" format=\"ascii\">\n") ;

        double coordx = cell(0,0) ;
        double coordy = cell(1,0) ;
        double coordz = cell(2,0) ;

        for (int i=0; i<npt[0]; i++)
        {
          coordy = cell(1,0) ;
          for (int j=0; j<npt[1]; j++)
          {
            coordz = cell(2,0) ;
            for (int k=0; k<npt[2]; k++)
            {
              fprintf(f, "%f %f %f\n", coordx, coordy, coordz) ;
              coordz = coordz + spacez ;
            }
            coordy = coordy + spacey ;
          }
          coordx = coordx + spacex ;
        }
        fprintf(f, "        </DataArray>\n") ;
        fprintf(f, "      </Points>\n") ;

        fprintf(f, "      <CellData>\n") ;
        fprintf(f, "      </CellData>\n") ;
        fprintf(f, "    </Piece>\n") ;
        fprintf(f, "  </StructuredGrid>\n") ;
        fprintf(f, "</VTKFile>\n") ;
        fclose(f);
    }
    world.gop.fence();
}

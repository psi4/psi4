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

/*
 *
 Dx_read.cc Code modifed by Kevin Hannon
11/13/2012
 */

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <string>
#include <sstream>
#include <vector>
#include <utility>

#include "psi4/libmints/basisset.h"
#include "psi4/libqt/qt.h"
#include "psi4/physconst.h"
#include "psi4/psi4-dec.h"



using namespace std;
using namespace psi;

namespace psi {

/*():Reads dx files
**
** KPH, November 2012
**
** \param V_eff      = matrix to diagonalize
** \param phi_ao      = dimension of A
** \param phi_so      = number of roots desired
** \param nmo    = # of molecule orbitals
** \param nso      = # of symmetry orbitals
** \param u = transformation matrix
*/
         void dx_read(double **V_eff, double *phi_ao,double *phi_so,int nao,int nso,double **u){

          int delta_count = 0;
          std::shared_ptr<BasisSet> basis;
          bool data_ready = false;
          int data_read = 0;
          double xstep, ystep, zstep;
          int xsteps, ysteps, zsteps;
          double xmin, ymin, zmin;
          int xinc = 0;
          int yinc = 0;
          int zinc = 0;
          int num_steps = 0;
          int total;
          double b2a3 = pc_bohr2angstroms * pc_bohr2angstroms * pc_bohr2angstroms;
          ifstream input;
          input.open("potential.dx");
          if (!input.good()) throw PSIEXCEPTION("Error opening potential.dx.");
          while (!input.eof())
          {
            char buf[512];
            input.getline(buf, 512);
            stringstream cppbuf(buf);
            string token;
            vector <string> tokens;
            while(cppbuf >> token) tokens.push_back(token);
            if(tokens.size()) { // skip blank lines

              if(data_ready && data_read <= total) { // data line
                for(int i=0; i < tokens.size(); i++) {
                  double x = xmin + ((double) xinc) * xstep;
                  double y = ymin + ((double) yinc) * ystep;
                  double z = zmin + ((double) zinc) * zstep;
                  double pot_val = atof(tokens[i].c_str());

                  basis->compute_phi(phi_ao, x, y, z);
                  // Transform phi_ao to SO basis
                  C_DGEMV('t', nao, nso, 1.0, &(u[0][0]), nso, &(phi_ao[0]), 1, 0.0, &(phi_so[0]), 1);
                  for(int i=0; i < nso; i++)
                    for(int j=0; j < nso; j++)
                      V_eff[i][j] += pot_val * xstep * ystep * zstep * phi_so[i] * phi_so[j] / b2a3;

//                  outfile->Printf( "x = %f; y = %f; z = %f; v = %f\n", x, y, z, pot_val);
                  num_steps++;
                  zinc++;
                  if(zinc == zsteps) { yinc++; zinc = 0; }
                  if(yinc == ysteps) { xinc++; yinc = 0; }
                  if(xinc == xsteps) {
                    outfile->Printf( "Total points read: %d\n", num_steps);
                    data_ready = false;
                  }
                }
              }

              if(tokens[0] == "#") // comment lines
                outfile->Printf( "%s\n", buf);
              if(tokens[0] == "origin") {
                xmin = atof(tokens[1].c_str());
                ymin = atof(tokens[2].c_str());
                zmin = atof(tokens[3].c_str());
                outfile->Printf( "%f %f %f\n", xmin, ymin, zmin);
              }

             if(tokens[0] == "delta") {
               if(delta_count == 0)
                 xstep = atof(tokens[1].c_str());
               else if(delta_count == 1)
                 ystep = atof(tokens[2].c_str());
               else if(delta_count == 2) {
                 zstep = atof(tokens[3].c_str());
                 outfile->Printf( "Step sizes: %f %f %f\n", xstep, ystep, zstep);
               }
               delta_count++;
             }

             if(tokens[0] == "object") {
               if(tokens[1] == "1") {
                 xsteps = atoi(tokens[5].c_str());
                 ysteps = atoi(tokens[6].c_str());
                 zsteps = atoi(tokens[7].c_str());
                 outfile->Printf( "%d %d %d\n", xsteps, ysteps, zsteps);
               }
               else if(tokens[1] == "3") {
                 total = atoi(tokens[9].c_str());
                 outfile->Printf( "%d\n", total);
                 data_ready = true;
               }
             }

            }
          }
          input.close();
        } // dx file

}

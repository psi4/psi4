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

#include <libmints/mints.h>
#include <libmints/basisset.h>
#include <libfunctional/superfunctional.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.hpp>
#include <libparallel/parallel.h>
#include <libiwl/iwl.hpp>
#include "qt.h"
#include <liboptions/python.h>
#include <psifiles.h>
#include <libfock/jk.h>
#include <physconst.h>
#include <psi4-dec.h>


using namespace boost;
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
          boost::shared_ptr<BasisSet> basis;
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
          double b2a3 = _bohr2angstroms * _bohr2angstroms * _bohr2angstroms;
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

//                  fprintf(outfile, "x = %f; y = %f; z = %f; v = %f\n", x, y, z, pot_val);
                  num_steps++;
                  zinc++;
                  if(zinc == zsteps) { yinc++; zinc = 0; }
                  if(yinc == ysteps) { xinc++; yinc = 0; }
                  if(xinc == xsteps) {
                    fprintf(outfile, "Total points read: %d\n", num_steps);
                    data_ready = false;
                  }
                }
              }

              if(tokens[0] == "#") // comment lines
                fprintf(outfile, "%s\n", buf);
              if(tokens[0] == "origin") {
                xmin = atof(tokens[1].c_str());
                ymin = atof(tokens[2].c_str());
                zmin = atof(tokens[3].c_str());
                fprintf(outfile, "%f %f %f\n", xmin, ymin, zmin);
              }

             if(tokens[0] == "delta") {
               if(delta_count == 0) 
                 xstep = atof(tokens[1].c_str());
               else if(delta_count == 1) 
                 ystep = atof(tokens[2].c_str());
               else if(delta_count == 2) {
                 zstep = atof(tokens[3].c_str());
                 fprintf(outfile, "Step sizes: %f %f %f\n", xstep, ystep, zstep);
               }
               delta_count++;
             }

             if(tokens[0] == "object") {
               if(tokens[1] == "1") {
                 xsteps = atoi(tokens[5].c_str());
                 ysteps = atoi(tokens[6].c_str());
                 zsteps = atoi(tokens[7].c_str());
                 fprintf(outfile, "%d %d %d\n", xsteps, ysteps, zsteps);
               }              
               else if(tokens[1] == "3") {
                 total = atoi(tokens[9].c_str());
                 fprintf(outfile, "%d\n", total);
                 data_ready = true;
               }
             }

            }
          }
          input.close();
        } // dx file
            
}


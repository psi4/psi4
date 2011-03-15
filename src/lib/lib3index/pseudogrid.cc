#include "3index.h"
#include <libmints/mints.h>
#include <libmints/integrator_defines.h>
#include <libqt/qt.h>
#include <boost/regex.hpp>
#include <boost/xpressive/xpressive.hpp>
#include <boost/xpressive/regex_actions.hpp>
#include <boost/algorithm/string.hpp>

#include <string>
#include <sstream>
#include <iostream>
#include <cstdio>
#include <fstream>
#include <algorithm>
#include <utility>
#include <ctype.h>

using namespace boost;
using namespace std;
using namespace psi;

namespace psi {
// the third parameter of from_string() should be
// one of std::hex, std::dec or std::oct
template <class T>
bool from_string(T& t,
                 const std::string& s,
                 std::ios_base& (*f)(std::ios_base&))
{
    std::istringstream iss(s);
    return !(iss >> f >> t).fail();
}

PseudoGrid::PseudoGrid(shared_ptr<Molecule> mol, const std::string& name) :
    molecule_(mol), name_(name)
{
}
PseudoGrid::~PseudoGrid()
{
    if (grid_.get() != NULL && grid_->getX() != NULL) {
        delete[] grid_->getX();
        delete[] grid_->getY();
        delete[] grid_->getZ();
        delete[] grid_->getWeights();
    }
}
void PseudoGrid::parse(const std::string& filename)
{
    std::string PSIDATADIR = Process::environment("PSIDATADIR");
    std::string gridname = PSIDATADIR + "grids/" + filename + ".grid";

    cout << gridname;

    // Phase I: read the grid in
    std::vector<std::vector<std::pair<int, std::pair<double, double> > > > array;
    array.resize(molecule_->natom());

    regex comment("^\\s*\\!.*");                                       // line starts with !
    regex separator("^\\*\\*\\*\\*");                                  // line starts with ****
    regex atom_array("^\\s*([A-Za-z]+)\\s+0.*");                       // array of atomic symbols terminated by 0

#define NUMBER "((?:[-+]?\\d*\\.\\d+(?:[DdEe][-+]?\\d+)?)|(?:[-+]?\\d+\\.\\d*(?:[DdEe][-+]?\\d+)?))"

    regex grid_shell("^\\s*(\\d+)\\s+" NUMBER ".*"); // match
    // Hold the result of a regex_match
    smatch what;

    std::vector<std::string> lines;
    std::string text;
    ifstream infile(gridname.c_str());
    if (!infile)
        throw PSIEXCEPTION("PseudoGridParser::parse: Unable to open pseudospectral grid file: " + filename);
    while (infile.good()) {
        getline(infile, text);
        lines.push_back(text);
    }

    for (int atom=0; atom<molecule_->natom(); ++atom) {

        int lineno = 0;
        bool found = false;
        while (lineno < lines.size()) {
            string line = lines[lineno++];

            // Spit out the line for debugging.
            //cout << line << endl;

            // Ignore blank lines
            if (line.empty())
                continue;

            // Do some matches
            if (regex_match(line, what, comment)) {
                //cout << " matched a comment line.\n";
                continue;
            }
            if (regex_match(line, what, separator)) {
                //cout << " matched a separator line.\n";
                continue;
            }
            // Match: H    0
            // or:    H    O...     0
            if (regex_match(line, what, atom_array)) {
                //cout << " matched a atom array line.\n";
                //cout << " what.captures(1)" << what[1].str() << "\n";

                // Check the captures and see if this basis set is for the atom we need.
                if (iequals(molecule_->label(atom), what[1].str())) {
                    found = true;
                    line = lines[lineno++];

                    // Need to do the following until we match a "****" which is the end of the basis set
                    while (!regex_match(line, what, separator)) {
                        // cout << " Atom line " << line;
                        if (regex_match(line, what, grid_shell)) {
                            int L;
                            double r;
                            if (!from_string<int>(L, what[1], std::dec))
                                throw PSIEXCEPTION("PseudoGridParser::parse: Unable to convert number of points (order):\n" + line);
                            if (!from_string<double>(r, what[2], std::dec))
                                throw PSIEXCEPTION("PseudoGridParser::parse: Unable to convert grid shell radius:\n" + line);
                            array[atom].push_back(make_pair(L, make_pair(r,1.0)));
                        }
                        line = lines[lineno++];
                    }
                    break;
                } else {
                    continue;
                }
            } else {
                continue;
            }

        }
        if (found == false)
            throw PSIEXCEPTION("PseudoGridParser::parser: Unable to find the grid for " + molecule_->label(atom));
    }

    int npoints = 0;
    for (int atom = 0; atom < molecule_->natom(); atom++) {

        for (std::vector<std::pair<int, std::pair<double, double> > >::iterator it = array[atom].begin(); it != array[atom].end(); it++) {
            std::pair<int, std::pair<double,double> > row = (*it);
            int L = row.first;
            // Defined in integrator_defines.h
            int max_grid = n_lebedev_;
            int ngrid;
            for (int index = 0; index < max_grid; index++) {
                if (lebedev_orders_[index] == L)
                    ngrid = lebedev_points_[index];
            }
            if (ngrid == 0)
                throw PSIEXCEPTION("Grid order does not match available Lebedev grid orders");
            npoints += ngrid;
        }
    }

    if (npoints == 0)
        throw PSIEXCEPTION("PseudoGridParser: No Grid points in this molecule");

    // Phase II: build the grid
    shared_ptr<Integrator> integrator(new Integrator(molecule_, shared_ptr<PSIO>(new PSIO())));
    grid_ = shared_ptr<GridBlock>(new GridBlock());
    grid_->setMaxPoints(npoints);
    grid_->setTruePoints(npoints);

    double* xp = new double[npoints];
    double* yp = new double[npoints];
    double* zp = new double[npoints];
    double* wp = new double[npoints];
    grid_->setGrid(xp,yp,zp,wp);

    int counter = 0;
    for (int atom = 0; atom < molecule_->natom(); atom++) {

        Vector3 center = molecule_->xyz(atom);
        double xc = center[0];
        double yc = center[1];
        double zc = center[2];

        for (std::vector<std::pair<int, std::pair<double, double> > >::iterator it = array[atom].begin(); it != array[atom].end(); it++) {
            std::pair<int, std::pair<double,double> > row = (*it);
            int L = row.first;
            double r = row.second.first;
            double w = row.second.second;

            int ngrid = 0;

            // Defined in integrator_defines.h
            int max_grid = n_lebedev_;
            for (int index = 0; index < max_grid; index++) {
                if (lebedev_orders_[index] == L)
                    ngrid = lebedev_points_[index];
            }
            if (ngrid == 0)
                throw PSIEXCEPTION("Grid order does not match available Lebedev grid orders");

            SphericalQuadrature quad = integrator->getLebedevSpherical(ngrid);

            for (int p = 0; p < ngrid; p++) {
                xp[counter + p] = quad.x[p];
                yp[counter + p] = quad.y[p];
                zp[counter + p] = quad.z[p];
                wp[counter + p] = quad.w[p];
            }
            free(quad.x);
            free(quad.y);
            free(quad.z);
            free(quad.w);

            // Scale to radius
            for (int l = 0; l < ngrid; l++) {
                xp[counter + l] *= r;
                yp[counter + l] *= r;
                zp[counter + l] *= r;
//                wp[counter + l] *= r*r;
            }

            // Center
            for (int l = 0; l < ngrid; l++) {
                xp[counter + l] += xc;
                yp[counter + l] += yc;
                zp[counter + l] += zc;
            }

            // TODO rotate to standard orientation
            // TODO scale for radial weight
            // TODO account for atomic cells

            counter += ngrid;
        }
    }

    fprintf(outfile, "  Pseudospectral Grid Points:\n");
    for (int P = 0; P < npoints; P++) {
        fprintf(outfile, "   P = %5d: (%24.16E, %24.16E, %24.16E) x %24.16E\n", P, xp[P], yp[P], zp[P], wp[P]);
    }

}


}

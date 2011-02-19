#include "3index.h"
#include <libmints/mints.h>
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
    std::vector<std::vector<std::pair<int, std::pair<double, double> > > > array;
    array.resize(molecule_->natom());
    
    regex comment("^\\s*\\!.*");                                       // line starts with !
    regex separator("^\\*\\*\\*\\*");                                  // line starts with ****
    regex atom_array("^\\s*([A-Za-z]+)\\s+0.*");                       // array of atomic symbols terminated by 0

#define NUMBER "((?:[-+]?\\d*\\.\\d+(?:[DdEe][-+]?\\d+)?)|(?:[-+]?\\d+\\.\\d*(?:[DdEe][-+]?\\d+)?))"

    regex grid_shell("^\\s*(\\d+)\\s+" NUMBER "\\s+" NUMBER ".*"); // match  
    // Hold the result of a regex_match
    smatch what;

    std::vector<std::string> lines;
    std::string text;
    ifstream infile(filename.c_str());
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
                            double w;
                            if (!from_string<int>(L, what[1], std::dec))
                                throw PSIEXCEPTION("PseudoGridParser::parse: Unable to convert number of points (order):\n" + line);
                            if (!from_string<double>(r, what[2], std::dec))
                                throw PSIEXCEPTION("PseudoGridParser::parse: Unable to convert grid shell radius:\n" + line);
                            if (!from_string<double>(w, what[3], std::dec))
                                throw PSIEXCEPTION("PseudoGridParser::parse: Unable to convert grid shell weight:\n" + line);
                            array[atom].push_back(make_pair(L, make_pair(r,w)));
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
            npoints += L;
        }
    }

    if (npoints == 0)
        throw PSIEXCEPTION("PseudoGridParser: No Grid points in this molecule");

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

            if (L == 6) {
                xp[counter +   0] =  1.0000000000000000E+000; yp[counter +   0] =  0.0000000000000000E+000; zp[counter +   0] =  0.0000000000000000E+000; wp[counter +   0] = w;
                xp[counter +   1] = -1.0000000000000000E+000; yp[counter +   1] =  0.0000000000000000E+000; zp[counter +   1] =  0.0000000000000000E+000; wp[counter +   1] = w;
                xp[counter +   2] =  0.0000000000000000E+000; yp[counter +   2] =  1.0000000000000000E+000; zp[counter +   2] =  0.0000000000000000E+000; wp[counter +   2] = w;
                xp[counter +   3] =  0.0000000000000000E+000; yp[counter +   3] = -1.0000000000000000E+000; zp[counter +   3] =  0.0000000000000000E+000; wp[counter +   3] = w;
                xp[counter +   4] =  0.0000000000000000E+000; yp[counter +   4] =  0.0000000000000000E+000; zp[counter +   4] =  1.0000000000000000E+000; wp[counter +   4] = w;
                xp[counter +   5] =  0.0000000000000000E+000; yp[counter +   5] =  0.0000000000000000E+000; zp[counter +   5] = -1.0000000000000000E+000; wp[counter +   5] = w;
            } else if (L == 8) {
                xp[counter +   0] =  5.7735026918962573E-001; yp[counter +   0] =  5.7735026918962573E-001; zp[counter +   0] =  5.7735026918962573E-001; wp[counter +   0] = w;
                xp[counter +   1] = -5.7735026918962573E-001; yp[counter +   1] =  5.7735026918962573E-001; zp[counter +   1] =  5.7735026918962573E-001; wp[counter +   1] = w;
                xp[counter +   2] =  5.7735026918962573E-001; yp[counter +   2] = -5.7735026918962573E-001; zp[counter +   2] =  5.7735026918962573E-001; wp[counter +   2] = w;
                xp[counter +   3] =  5.7735026918962573E-001; yp[counter +   3] =  5.7735026918962573E-001; zp[counter +   3] = -5.7735026918962573E-001; wp[counter +   3] = w;
                xp[counter +   4] = -5.7735026918962573E-001; yp[counter +   4] = -5.7735026918962573E-001; zp[counter +   4] =  5.7735026918962573E-001; wp[counter +   4] = w;
                xp[counter +   5] =  5.7735026918962573E-001; yp[counter +   5] = -5.7735026918962573E-001; zp[counter +   5] = -5.7735026918962573E-001; wp[counter +   5] = w;
                xp[counter +   6] = -5.7735026918962573E-001; yp[counter +   6] =  5.7735026918962573E-001; zp[counter +   6] = -5.7735026918962573E-001; wp[counter +   6] = w;
                xp[counter +   7] = -5.7735026918962573E-001; yp[counter +   7] = -5.7735026918962573E-001; zp[counter +   7] = -5.7735026918962573E-001; wp[counter +   7] = w;
            } else if (L == 14) {                                                                                                                                 
                xp[counter +   0] =  1.0000000000000000E+000; yp[counter +   0] =  0.0000000000000000E+000; zp[counter +   0] =  0.0000000000000000E+000; wp[counter +   0] = w;
                xp[counter +   1] = -1.0000000000000000E+000; yp[counter +   1] =  0.0000000000000000E+000; zp[counter +   1] =  0.0000000000000000E+000; wp[counter +   1] = w;
                xp[counter +   2] =  0.0000000000000000E+000; yp[counter +   2] =  1.0000000000000000E+000; zp[counter +   2] =  0.0000000000000000E+000; wp[counter +   2] = w;
                xp[counter +   3] =  0.0000000000000000E+000; yp[counter +   3] = -1.0000000000000000E+000; zp[counter +   3] =  0.0000000000000000E+000; wp[counter +   3] = w;
                xp[counter +   4] =  0.0000000000000000E+000; yp[counter +   4] =  0.0000000000000000E+000; zp[counter +   4] =  1.0000000000000000E+000; wp[counter +   4] = w;
                xp[counter +   5] =  0.0000000000000000E+000; yp[counter +   5] =  0.0000000000000000E+000; zp[counter +   5] = -1.0000000000000000E+000; wp[counter +   5] = w;
                xp[counter +   6] =  5.7735026918962573E-001; yp[counter +   6] =  5.7735026918962573E-001; zp[counter +   6] =  5.7735026918962573E-001; wp[counter +   6] = w;
                xp[counter +   7] = -5.7735026918962573E-001; yp[counter +   7] =  5.7735026918962573E-001; zp[counter +   7] =  5.7735026918962573E-001; wp[counter +   7] = w;
                xp[counter +   8] =  5.7735026918962573E-001; yp[counter +   8] = -5.7735026918962573E-001; zp[counter +   8] =  5.7735026918962573E-001; wp[counter +   8] = w;
                xp[counter +   9] =  5.7735026918962573E-001; yp[counter +   9] =  5.7735026918962573E-001; zp[counter +   9] = -5.7735026918962573E-001; wp[counter +   9] = w;
                xp[counter +  10] = -5.7735026918962573E-001; yp[counter +  10] = -5.7735026918962573E-001; zp[counter +  10] =  5.7735026918962573E-001; wp[counter +  10] = w;
                xp[counter +  11] =  5.7735026918962573E-001; yp[counter +  11] = -5.7735026918962573E-001; zp[counter +  11] = -5.7735026918962573E-001; wp[counter +  11] = w;
                xp[counter +  12] = -5.7735026918962573E-001; yp[counter +  12] =  5.7735026918962573E-001; zp[counter +  12] = -5.7735026918962573E-001; wp[counter +  12] = w;
                xp[counter +  13] = -5.7735026918962573E-001; yp[counter +  13] = -5.7735026918962573E-001; zp[counter +  13] = -5.7735026918962573E-001; wp[counter +  13] = w;
            } else if (L == 26) {
                xp[counter +   0] =  1.0000000000000000E+000; yp[counter +   0] =  0.0000000000000000E+000; zp[counter +   0] =  0.0000000000000000E+000; wp[counter +   0] = w;
                xp[counter +   1] = -1.0000000000000000E+000; yp[counter +   1] =  0.0000000000000000E+000; zp[counter +   1] =  0.0000000000000000E+000; wp[counter +   1] = w;
                xp[counter +   2] =  0.0000000000000000E+000; yp[counter +   2] =  1.0000000000000000E+000; zp[counter +   2] =  0.0000000000000000E+000; wp[counter +   2] = w;
                xp[counter +   3] =  0.0000000000000000E+000; yp[counter +   3] = -1.0000000000000000E+000; zp[counter +   3] =  0.0000000000000000E+000; wp[counter +   3] = w;
                xp[counter +   4] =  0.0000000000000000E+000; yp[counter +   4] =  0.0000000000000000E+000; zp[counter +   4] =  1.0000000000000000E+000; wp[counter +   4] = w;
                xp[counter +   5] =  0.0000000000000000E+000; yp[counter +   5] =  0.0000000000000000E+000; zp[counter +   5] = -1.0000000000000000E+000; wp[counter +   5] = w;
                xp[counter +   6] =  0.0000000000000000E+000; yp[counter +   6] =  7.0710678118654757E-001; zp[counter +   6] =  7.0710678118654757E-001; wp[counter +   6] = w;
                xp[counter +   7] =  0.0000000000000000E+000; yp[counter +   7] = -7.0710678118654757E-001; zp[counter +   7] =  7.0710678118654757E-001; wp[counter +   7] = w;
                xp[counter +   8] =  0.0000000000000000E+000; yp[counter +   8] =  7.0710678118654757E-001; zp[counter +   8] = -7.0710678118654757E-001; wp[counter +   8] = w;
                xp[counter +   9] =  0.0000000000000000E+000; yp[counter +   9] = -7.0710678118654757E-001; zp[counter +   9] = -7.0710678118654757E-001; wp[counter +   9] = w;
                xp[counter +  10] =  7.0710678118654757E-001; yp[counter +  10] =  0.0000000000000000E+000; zp[counter +  10] =  7.0710678118654757E-001; wp[counter +  10] = w;
                xp[counter +  11] =  7.0710678118654757E-001; yp[counter +  11] =  0.0000000000000000E+000; zp[counter +  11] = -7.0710678118654757E-001; wp[counter +  11] = w;
                xp[counter +  12] = -7.0710678118654757E-001; yp[counter +  12] =  0.0000000000000000E+000; zp[counter +  12] =  7.0710678118654757E-001; wp[counter +  12] = w;
                xp[counter +  13] = -7.0710678118654757E-001; yp[counter +  13] =  0.0000000000000000E+000; zp[counter +  13] = -7.0710678118654757E-001; wp[counter +  13] = w;
                xp[counter +  14] =  7.0710678118654757E-001; yp[counter +  14] =  7.0710678118654757E-001; zp[counter +  14] =  0.0000000000000000E+000; wp[counter +  14] = w;
                xp[counter +  15] = -7.0710678118654757E-001; yp[counter +  15] =  7.0710678118654757E-001; zp[counter +  15] =  0.0000000000000000E+000; wp[counter +  15] = w;
                xp[counter +  16] =  7.0710678118654757E-001; yp[counter +  16] = -7.0710678118654757E-001; zp[counter +  16] =  0.0000000000000000E+000; wp[counter +  16] = w;
                xp[counter +  17] = -7.0710678118654757E-001; yp[counter +  17] = -7.0710678118654757E-001; zp[counter +  17] =  0.0000000000000000E+000; wp[counter +  17] = w;
                xp[counter +  18] =  5.7735026918962573E-001; yp[counter +  18] =  5.7735026918962573E-001; zp[counter +  18] =  5.7735026918962573E-001; wp[counter +  18] = w;
                xp[counter +  19] = -5.7735026918962573E-001; yp[counter +  19] =  5.7735026918962573E-001; zp[counter +  19] =  5.7735026918962573E-001; wp[counter +  19] = w;
                xp[counter +  20] =  5.7735026918962573E-001; yp[counter +  20] = -5.7735026918962573E-001; zp[counter +  20] =  5.7735026918962573E-001; wp[counter +  20] = w;
                xp[counter +  21] =  5.7735026918962573E-001; yp[counter +  21] =  5.7735026918962573E-001; zp[counter +  21] = -5.7735026918962573E-001; wp[counter +  21] = w;
                xp[counter +  22] = -5.7735026918962573E-001; yp[counter +  22] = -5.7735026918962573E-001; zp[counter +  22] =  5.7735026918962573E-001; wp[counter +  22] = w;
                xp[counter +  23] =  5.7735026918962573E-001; yp[counter +  23] = -5.7735026918962573E-001; zp[counter +  23] = -5.7735026918962573E-001; wp[counter +  23] = w;
                xp[counter +  24] = -5.7735026918962573E-001; yp[counter +  24] =  5.7735026918962573E-001; zp[counter +  24] = -5.7735026918962573E-001; wp[counter +  24] = w;
                xp[counter +  25] = -5.7735026918962573E-001; yp[counter +  25] = -5.7735026918962573E-001; zp[counter +  25] = -5.7735026918962573E-001; wp[counter +  25] = w;
            } else if (L == 38) {
                xp[counter +   0] =  1.0000000000000000E+000; yp[counter +   0] =  0.0000000000000000E+000; zp[counter +   0] =  0.0000000000000000E+000; wp[counter +   0] = w;
                xp[counter +   1] = -1.0000000000000000E+000; yp[counter +   1] =  0.0000000000000000E+000; zp[counter +   1] =  0.0000000000000000E+000; wp[counter +   1] = w;
                xp[counter +   2] =  0.0000000000000000E+000; yp[counter +   2] =  1.0000000000000000E+000; zp[counter +   2] =  0.0000000000000000E+000; wp[counter +   2] = w;
                xp[counter +   3] =  0.0000000000000000E+000; yp[counter +   3] = -1.0000000000000000E+000; zp[counter +   3] =  0.0000000000000000E+000; wp[counter +   3] = w;
                xp[counter +   4] =  0.0000000000000000E+000; yp[counter +   4] =  0.0000000000000000E+000; zp[counter +   4] =  1.0000000000000000E+000; wp[counter +   4] = w;
                xp[counter +   5] =  0.0000000000000000E+000; yp[counter +   5] =  0.0000000000000000E+000; zp[counter +   5] = -1.0000000000000000E+000; wp[counter +   5] = w;
                xp[counter +   6] =  5.7735026918962573E-001; yp[counter +   6] =  5.7735026918962573E-001; zp[counter +   6] =  5.7735026918962573E-001; wp[counter +   6] = w;
                xp[counter +   7] = -5.7735026918962573E-001; yp[counter +   7] =  5.7735026918962573E-001; zp[counter +   7] =  5.7735026918962573E-001; wp[counter +   7] = w;
                xp[counter +   8] =  5.7735026918962573E-001; yp[counter +   8] = -5.7735026918962573E-001; zp[counter +   8] =  5.7735026918962573E-001; wp[counter +   8] = w;
                xp[counter +   9] =  5.7735026918962573E-001; yp[counter +   9] =  5.7735026918962573E-001; zp[counter +   9] = -5.7735026918962573E-001; wp[counter +   9] = w;
                xp[counter +  10] = -5.7735026918962573E-001; yp[counter +  10] = -5.7735026918962573E-001; zp[counter +  10] =  5.7735026918962573E-001; wp[counter +  10] = w;
                xp[counter +  11] =  5.7735026918962573E-001; yp[counter +  11] = -5.7735026918962573E-001; zp[counter +  11] = -5.7735026918962573E-001; wp[counter +  11] = w;
                xp[counter +  12] = -5.7735026918962573E-001; yp[counter +  12] =  5.7735026918962573E-001; zp[counter +  12] = -5.7735026918962573E-001; wp[counter +  12] = w;
                xp[counter +  13] = -5.7735026918962573E-001; yp[counter +  13] = -5.7735026918962573E-001; zp[counter +  13] = -5.7735026918962573E-001; wp[counter +  13] = w;
                xp[counter +  14] =  4.5970084338098310E-001; yp[counter +  14] =  8.8807383397711526E-001; zp[counter +  14] =  0.0000000000000000E+000; wp[counter +  14] = w;
                xp[counter +  15] = -4.5970084338098310E-001; yp[counter +  15] =  8.8807383397711526E-001; zp[counter +  15] =  0.0000000000000000E+000; wp[counter +  15] = w;
                xp[counter +  16] =  4.5970084338098310E-001; yp[counter +  16] = -8.8807383397711526E-001; zp[counter +  16] =  0.0000000000000000E+000; wp[counter +  16] = w;
                xp[counter +  17] = -4.5970084338098310E-001; yp[counter +  17] = -8.8807383397711526E-001; zp[counter +  17] =  0.0000000000000000E+000; wp[counter +  17] = w;
                xp[counter +  18] =  8.8807383397711526E-001; yp[counter +  18] =  4.5970084338098310E-001; zp[counter +  18] =  0.0000000000000000E+000; wp[counter +  18] = w;
                xp[counter +  19] = -8.8807383397711526E-001; yp[counter +  19] =  4.5970084338098310E-001; zp[counter +  19] =  0.0000000000000000E+000; wp[counter +  19] = w;
                xp[counter +  20] =  8.8807383397711526E-001; yp[counter +  20] = -4.5970084338098310E-001; zp[counter +  20] =  0.0000000000000000E+000; wp[counter +  20] = w;
                xp[counter +  21] = -8.8807383397711526E-001; yp[counter +  21] = -4.5970084338098310E-001; zp[counter +  21] =  0.0000000000000000E+000; wp[counter +  21] = w;
                xp[counter +  22] =  4.5970084338098310E-001; yp[counter +  22] =  0.0000000000000000E+000; zp[counter +  22] =  8.8807383397711526E-001; wp[counter +  22] = w;
                xp[counter +  23] = -4.5970084338098310E-001; yp[counter +  23] =  0.0000000000000000E+000; zp[counter +  23] =  8.8807383397711526E-001; wp[counter +  23] = w;
                xp[counter +  24] =  4.5970084338098310E-001; yp[counter +  24] =  0.0000000000000000E+000; zp[counter +  24] = -8.8807383397711526E-001; wp[counter +  24] = w;
                xp[counter +  25] = -4.5970084338098310E-001; yp[counter +  25] =  0.0000000000000000E+000; zp[counter +  25] = -8.8807383397711526E-001; wp[counter +  25] = w;
                xp[counter +  26] =  8.8807383397711526E-001; yp[counter +  26] =  0.0000000000000000E+000; zp[counter +  26] =  4.5970084338098310E-001; wp[counter +  26] = w;
                xp[counter +  27] = -8.8807383397711526E-001; yp[counter +  27] =  0.0000000000000000E+000; zp[counter +  27] =  4.5970084338098310E-001; wp[counter +  27] = w;
                xp[counter +  28] =  8.8807383397711526E-001; yp[counter +  28] =  0.0000000000000000E+000; zp[counter +  28] = -4.5970084338098310E-001; wp[counter +  28] = w;
                xp[counter +  29] = -8.8807383397711526E-001; yp[counter +  29] =  0.0000000000000000E+000; zp[counter +  29] = -4.5970084338098310E-001; wp[counter +  29] = w;
                xp[counter +  30] =  0.0000000000000000E+000; yp[counter +  30] =  4.5970084338098310E-001; zp[counter +  30] =  8.8807383397711526E-001; wp[counter +  30] = w;
                xp[counter +  31] =  0.0000000000000000E+000; yp[counter +  31] = -4.5970084338098310E-001; zp[counter +  31] =  8.8807383397711526E-001; wp[counter +  31] = w;
                xp[counter +  32] =  0.0000000000000000E+000; yp[counter +  32] =  4.5970084338098310E-001; zp[counter +  32] = -8.8807383397711526E-001; wp[counter +  32] = w;
                xp[counter +  33] =  0.0000000000000000E+000; yp[counter +  33] = -4.5970084338098310E-001; zp[counter +  33] = -8.8807383397711526E-001; wp[counter +  33] = w;
                xp[counter +  34] =  0.0000000000000000E+000; yp[counter +  34] =  8.8807383397711526E-001; zp[counter +  34] =  4.5970084338098310E-001; wp[counter +  34] = w;
                xp[counter +  35] =  0.0000000000000000E+000; yp[counter +  35] = -8.8807383397711526E-001; zp[counter +  35] =  4.5970084338098310E-001; wp[counter +  35] = w;
                xp[counter +  36] =  0.0000000000000000E+000; yp[counter +  36] =  8.8807383397711526E-001; zp[counter +  36] = -4.5970084338098310E-001; wp[counter +  36] = w;
                xp[counter +  37] =  0.0000000000000000E+000; yp[counter +  37] = -8.8807383397711526E-001; zp[counter +  37] = -4.5970084338098310E-001; wp[counter +  37] = w;
            } else if (L == 50) {
                xp[counter +   0] =  1.0000000000000000E+000; yp[counter +   0] =  0.0000000000000000E+000; zp[counter +   0] =  0.0000000000000000E+000; wp[counter +   0] = w;
                xp[counter +   1] = -1.0000000000000000E+000; yp[counter +   1] =  0.0000000000000000E+000; zp[counter +   1] =  0.0000000000000000E+000; wp[counter +   1] = w;
                xp[counter +   2] =  0.0000000000000000E+000; yp[counter +   2] =  1.0000000000000000E+000; zp[counter +   2] =  0.0000000000000000E+000; wp[counter +   2] = w;
                xp[counter +   3] =  0.0000000000000000E+000; yp[counter +   3] = -1.0000000000000000E+000; zp[counter +   3] =  0.0000000000000000E+000; wp[counter +   3] = w;
                xp[counter +   4] =  0.0000000000000000E+000; yp[counter +   4] =  0.0000000000000000E+000; zp[counter +   4] =  1.0000000000000000E+000; wp[counter +   4] = w;
                xp[counter +   5] =  0.0000000000000000E+000; yp[counter +   5] =  0.0000000000000000E+000; zp[counter +   5] = -1.0000000000000000E+000; wp[counter +   5] = w;
                xp[counter +   6] =  0.0000000000000000E+000; yp[counter +   6] =  7.0710678118654757E-001; zp[counter +   6] =  7.0710678118654757E-001; wp[counter +   6] = w;
                xp[counter +   7] =  0.0000000000000000E+000; yp[counter +   7] = -7.0710678118654757E-001; zp[counter +   7] =  7.0710678118654757E-001; wp[counter +   7] = w;
                xp[counter +   8] =  0.0000000000000000E+000; yp[counter +   8] =  7.0710678118654757E-001; zp[counter +   8] = -7.0710678118654757E-001; wp[counter +   8] = w;
                xp[counter +   9] =  0.0000000000000000E+000; yp[counter +   9] = -7.0710678118654757E-001; zp[counter +   9] = -7.0710678118654757E-001; wp[counter +   9] = w;
                xp[counter +  10] =  7.0710678118654757E-001; yp[counter +  10] =  0.0000000000000000E+000; zp[counter +  10] =  7.0710678118654757E-001; wp[counter +  10] = w;
                xp[counter +  11] =  7.0710678118654757E-001; yp[counter +  11] =  0.0000000000000000E+000; zp[counter +  11] = -7.0710678118654757E-001; wp[counter +  11] = w;
                xp[counter +  12] = -7.0710678118654757E-001; yp[counter +  12] =  0.0000000000000000E+000; zp[counter +  12] =  7.0710678118654757E-001; wp[counter +  12] = w;
                xp[counter +  13] = -7.0710678118654757E-001; yp[counter +  13] =  0.0000000000000000E+000; zp[counter +  13] = -7.0710678118654757E-001; wp[counter +  13] = w;
                xp[counter +  14] =  7.0710678118654757E-001; yp[counter +  14] =  7.0710678118654757E-001; zp[counter +  14] =  0.0000000000000000E+000; wp[counter +  14] = w;
                xp[counter +  15] = -7.0710678118654757E-001; yp[counter +  15] =  7.0710678118654757E-001; zp[counter +  15] =  0.0000000000000000E+000; wp[counter +  15] = w;
                xp[counter +  16] =  7.0710678118654757E-001; yp[counter +  16] = -7.0710678118654757E-001; zp[counter +  16] =  0.0000000000000000E+000; wp[counter +  16] = w;
                xp[counter +  17] = -7.0710678118654757E-001; yp[counter +  17] = -7.0710678118654757E-001; zp[counter +  17] =  0.0000000000000000E+000; wp[counter +  17] = w;
                xp[counter +  18] =  5.7735026918962573E-001; yp[counter +  18] =  5.7735026918962573E-001; zp[counter +  18] =  5.7735026918962573E-001; wp[counter +  18] = w;
                xp[counter +  19] = -5.7735026918962573E-001; yp[counter +  19] =  5.7735026918962573E-001; zp[counter +  19] =  5.7735026918962573E-001; wp[counter +  19] = w;
                xp[counter +  20] =  5.7735026918962573E-001; yp[counter +  20] = -5.7735026918962573E-001; zp[counter +  20] =  5.7735026918962573E-001; wp[counter +  20] = w;
                xp[counter +  21] =  5.7735026918962573E-001; yp[counter +  21] =  5.7735026918962573E-001; zp[counter +  21] = -5.7735026918962573E-001; wp[counter +  21] = w;
                xp[counter +  22] = -5.7735026918962573E-001; yp[counter +  22] = -5.7735026918962573E-001; zp[counter +  22] =  5.7735026918962573E-001; wp[counter +  22] = w;
                xp[counter +  23] =  5.7735026918962573E-001; yp[counter +  23] = -5.7735026918962573E-001; zp[counter +  23] = -5.7735026918962573E-001; wp[counter +  23] = w;
                xp[counter +  24] = -5.7735026918962573E-001; yp[counter +  24] =  5.7735026918962573E-001; zp[counter +  24] = -5.7735026918962573E-001; wp[counter +  24] = w;
                xp[counter +  25] = -5.7735026918962573E-001; yp[counter +  25] = -5.7735026918962573E-001; zp[counter +  25] = -5.7735026918962573E-001; wp[counter +  25] = w;
                xp[counter +  26] =  3.0151134457776357E-001; yp[counter +  26] =  3.0151134457776357E-001; zp[counter +  26] =  9.0453403373329089E-001; wp[counter +  26] = w;
                xp[counter +  27] = -3.0151134457776357E-001; yp[counter +  27] =  3.0151134457776357E-001; zp[counter +  27] =  9.0453403373329089E-001; wp[counter +  27] = w;
                xp[counter +  28] =  3.0151134457776357E-001; yp[counter +  28] = -3.0151134457776357E-001; zp[counter +  28] =  9.0453403373329089E-001; wp[counter +  28] = w;
                xp[counter +  29] =  3.0151134457776357E-001; yp[counter +  29] =  3.0151134457776357E-001; zp[counter +  29] = -9.0453403373329089E-001; wp[counter +  29] = w;
                xp[counter +  30] = -3.0151134457776357E-001; yp[counter +  30] = -3.0151134457776357E-001; zp[counter +  30] =  9.0453403373329089E-001; wp[counter +  30] = w;
                xp[counter +  31] = -3.0151134457776357E-001; yp[counter +  31] =  3.0151134457776357E-001; zp[counter +  31] = -9.0453403373329089E-001; wp[counter +  31] = w;
                xp[counter +  32] =  3.0151134457776357E-001; yp[counter +  32] = -3.0151134457776357E-001; zp[counter +  32] = -9.0453403373329089E-001; wp[counter +  32] = w;
                xp[counter +  33] = -3.0151134457776357E-001; yp[counter +  33] = -3.0151134457776357E-001; zp[counter +  33] = -9.0453403373329089E-001; wp[counter +  33] = w;
                xp[counter +  34] = -3.0151134457776357E-001; yp[counter +  34] =  9.0453403373329089E-001; zp[counter +  34] =  3.0151134457776357E-001; wp[counter +  34] = w;
                xp[counter +  35] =  3.0151134457776357E-001; yp[counter +  35] = -9.0453403373329089E-001; zp[counter +  35] =  3.0151134457776357E-001; wp[counter +  35] = w;
                xp[counter +  36] =  3.0151134457776357E-001; yp[counter +  36] =  9.0453403373329089E-001; zp[counter +  36] = -3.0151134457776357E-001; wp[counter +  36] = w;
                xp[counter +  37] = -3.0151134457776357E-001; yp[counter +  37] = -9.0453403373329089E-001; zp[counter +  37] =  3.0151134457776357E-001; wp[counter +  37] = w;
                xp[counter +  38] = -3.0151134457776357E-001; yp[counter +  38] =  9.0453403373329089E-001; zp[counter +  38] = -3.0151134457776357E-001; wp[counter +  38] = w;
                xp[counter +  39] =  3.0151134457776357E-001; yp[counter +  39] = -9.0453403373329089E-001; zp[counter +  39] = -3.0151134457776357E-001; wp[counter +  39] = w;
                xp[counter +  40] = -3.0151134457776357E-001; yp[counter +  40] = -9.0453403373329089E-001; zp[counter +  40] = -3.0151134457776357E-001; wp[counter +  40] = w;
                xp[counter +  41] =  3.0151134457776357E-001; yp[counter +  41] =  9.0453403373329089E-001; zp[counter +  41] =  3.0151134457776357E-001; wp[counter +  41] = w;
                xp[counter +  42] =  9.0453403373329089E-001; yp[counter +  42] =  3.0151134457776357E-001; zp[counter +  42] =  3.0151134457776357E-001; wp[counter +  42] = w;
                xp[counter +  43] = -9.0453403373329089E-001; yp[counter +  43] =  3.0151134457776357E-001; zp[counter +  43] =  3.0151134457776357E-001; wp[counter +  43] = w;
                xp[counter +  44] =  9.0453403373329089E-001; yp[counter +  44] = -3.0151134457776357E-001; zp[counter +  44] =  3.0151134457776357E-001; wp[counter +  44] = w;
                xp[counter +  45] =  9.0453403373329089E-001; yp[counter +  45] =  3.0151134457776357E-001; zp[counter +  45] = -3.0151134457776357E-001; wp[counter +  45] = w;
                xp[counter +  46] = -9.0453403373329089E-001; yp[counter +  46] = -3.0151134457776357E-001; zp[counter +  46] =  3.0151134457776357E-001; wp[counter +  46] = w;
                xp[counter +  47] = -9.0453403373329089E-001; yp[counter +  47] =  3.0151134457776357E-001; zp[counter +  47] = -3.0151134457776357E-001; wp[counter +  47] = w;
                xp[counter +  48] =  9.0453403373329089E-001; yp[counter +  48] = -3.0151134457776357E-001; zp[counter +  48] = -3.0151134457776357E-001; wp[counter +  48] = w;
                xp[counter +  49] = -9.0453403373329089E-001; yp[counter +  49] = -3.0151134457776357E-001; zp[counter +  49] = -3.0151134457776357E-001; wp[counter +  49] = w;
            } else if (L == 74) {
                xp[counter +   0] =  1.0000000000000000E+000; yp[counter +   0] =  0.0000000000000000E+000; zp[counter +   0] =  0.0000000000000000E+000; wp[counter +   0] = w;
                xp[counter +   1] = -1.0000000000000000E+000; yp[counter +   1] =  0.0000000000000000E+000; zp[counter +   1] =  0.0000000000000000E+000; wp[counter +   1] = w;
                xp[counter +   2] =  0.0000000000000000E+000; yp[counter +   2] =  1.0000000000000000E+000; zp[counter +   2] =  0.0000000000000000E+000; wp[counter +   2] = w;
                xp[counter +   3] =  0.0000000000000000E+000; yp[counter +   3] = -1.0000000000000000E+000; zp[counter +   3] =  0.0000000000000000E+000; wp[counter +   3] = w;
                xp[counter +   4] =  0.0000000000000000E+000; yp[counter +   4] =  0.0000000000000000E+000; zp[counter +   4] =  1.0000000000000000E+000; wp[counter +   4] = w;
                xp[counter +   5] =  0.0000000000000000E+000; yp[counter +   5] =  0.0000000000000000E+000; zp[counter +   5] = -1.0000000000000000E+000; wp[counter +   5] = w;
                xp[counter +   6] =  0.0000000000000000E+000; yp[counter +   6] =  7.0710678118654757E-001; zp[counter +   6] =  7.0710678118654757E-001; wp[counter +   6] = w;
                xp[counter +   7] =  0.0000000000000000E+000; yp[counter +   7] = -7.0710678118654757E-001; zp[counter +   7] =  7.0710678118654757E-001; wp[counter +   7] = w;
                xp[counter +   8] =  0.0000000000000000E+000; yp[counter +   8] =  7.0710678118654757E-001; zp[counter +   8] = -7.0710678118654757E-001; wp[counter +   8] = w;
                xp[counter +   9] =  0.0000000000000000E+000; yp[counter +   9] = -7.0710678118654757E-001; zp[counter +   9] = -7.0710678118654757E-001; wp[counter +   9] = w;
                xp[counter +  10] =  7.0710678118654757E-001; yp[counter +  10] =  0.0000000000000000E+000; zp[counter +  10] =  7.0710678118654757E-001; wp[counter +  10] = w;
                xp[counter +  11] =  7.0710678118654757E-001; yp[counter +  11] =  0.0000000000000000E+000; zp[counter +  11] = -7.0710678118654757E-001; wp[counter +  11] = w;
                xp[counter +  12] = -7.0710678118654757E-001; yp[counter +  12] =  0.0000000000000000E+000; zp[counter +  12] =  7.0710678118654757E-001; wp[counter +  12] = w;
                xp[counter +  13] = -7.0710678118654757E-001; yp[counter +  13] =  0.0000000000000000E+000; zp[counter +  13] = -7.0710678118654757E-001; wp[counter +  13] = w;
                xp[counter +  14] =  7.0710678118654757E-001; yp[counter +  14] =  7.0710678118654757E-001; zp[counter +  14] =  0.0000000000000000E+000; wp[counter +  14] = w;
                xp[counter +  15] = -7.0710678118654757E-001; yp[counter +  15] =  7.0710678118654757E-001; zp[counter +  15] =  0.0000000000000000E+000; wp[counter +  15] = w;
                xp[counter +  16] =  7.0710678118654757E-001; yp[counter +  16] = -7.0710678118654757E-001; zp[counter +  16] =  0.0000000000000000E+000; wp[counter +  16] = w;
                xp[counter +  17] = -7.0710678118654757E-001; yp[counter +  17] = -7.0710678118654757E-001; zp[counter +  17] =  0.0000000000000000E+000; wp[counter +  17] = w;
                xp[counter +  18] =  5.7735026918962573E-001; yp[counter +  18] =  5.7735026918962573E-001; zp[counter +  18] =  5.7735026918962573E-001; wp[counter +  18] = w;
                xp[counter +  19] = -5.7735026918962573E-001; yp[counter +  19] =  5.7735026918962573E-001; zp[counter +  19] =  5.7735026918962573E-001; wp[counter +  19] = w;
                xp[counter +  20] =  5.7735026918962573E-001; yp[counter +  20] = -5.7735026918962573E-001; zp[counter +  20] =  5.7735026918962573E-001; wp[counter +  20] = w;
                xp[counter +  21] =  5.7735026918962573E-001; yp[counter +  21] =  5.7735026918962573E-001; zp[counter +  21] = -5.7735026918962573E-001; wp[counter +  21] = w;
                xp[counter +  22] = -5.7735026918962573E-001; yp[counter +  22] = -5.7735026918962573E-001; zp[counter +  22] =  5.7735026918962573E-001; wp[counter +  22] = w;
                xp[counter +  23] =  5.7735026918962573E-001; yp[counter +  23] = -5.7735026918962573E-001; zp[counter +  23] = -5.7735026918962573E-001; wp[counter +  23] = w;
                xp[counter +  24] = -5.7735026918962573E-001; yp[counter +  24] =  5.7735026918962573E-001; zp[counter +  24] = -5.7735026918962573E-001; wp[counter +  24] = w;
                xp[counter +  25] = -5.7735026918962573E-001; yp[counter +  25] = -5.7735026918962573E-001; zp[counter +  25] = -5.7735026918962573E-001; wp[counter +  25] = w;
                xp[counter +  26] =  4.8038446141526142E-001; yp[counter +  26] =  4.8038446141526142E-001; zp[counter +  26] =  7.3379938570534275E-001; wp[counter +  26] = w;
                xp[counter +  27] = -4.8038446141526142E-001; yp[counter +  27] =  4.8038446141526142E-001; zp[counter +  27] =  7.3379938570534275E-001; wp[counter +  27] = w;
                xp[counter +  28] =  4.8038446141526142E-001; yp[counter +  28] = -4.8038446141526142E-001; zp[counter +  28] =  7.3379938570534275E-001; wp[counter +  28] = w;
                xp[counter +  29] =  4.8038446141526142E-001; yp[counter +  29] =  4.8038446141526142E-001; zp[counter +  29] = -7.3379938570534275E-001; wp[counter +  29] = w;
                xp[counter +  30] = -4.8038446141526142E-001; yp[counter +  30] = -4.8038446141526142E-001; zp[counter +  30] =  7.3379938570534275E-001; wp[counter +  30] = w;
                xp[counter +  31] = -4.8038446141526142E-001; yp[counter +  31] =  4.8038446141526142E-001; zp[counter +  31] = -7.3379938570534275E-001; wp[counter +  31] = w;
                xp[counter +  32] =  4.8038446141526142E-001; yp[counter +  32] = -4.8038446141526142E-001; zp[counter +  32] = -7.3379938570534275E-001; wp[counter +  32] = w;
                xp[counter +  33] = -4.8038446141526142E-001; yp[counter +  33] = -4.8038446141526142E-001; zp[counter +  33] = -7.3379938570534275E-001; wp[counter +  33] = w;
                xp[counter +  34] = -4.8038446141526142E-001; yp[counter +  34] =  7.3379938570534275E-001; zp[counter +  34] =  4.8038446141526142E-001; wp[counter +  34] = w;
                xp[counter +  35] =  4.8038446141526142E-001; yp[counter +  35] = -7.3379938570534275E-001; zp[counter +  35] =  4.8038446141526142E-001; wp[counter +  35] = w;
                xp[counter +  36] =  4.8038446141526142E-001; yp[counter +  36] =  7.3379938570534275E-001; zp[counter +  36] = -4.8038446141526142E-001; wp[counter +  36] = w;
                xp[counter +  37] = -4.8038446141526142E-001; yp[counter +  37] = -7.3379938570534275E-001; zp[counter +  37] =  4.8038446141526142E-001; wp[counter +  37] = w;
                xp[counter +  38] = -4.8038446141526142E-001; yp[counter +  38] =  7.3379938570534275E-001; zp[counter +  38] = -4.8038446141526142E-001; wp[counter +  38] = w;
                xp[counter +  39] =  4.8038446141526142E-001; yp[counter +  39] = -7.3379938570534275E-001; zp[counter +  39] = -4.8038446141526142E-001; wp[counter +  39] = w;
                xp[counter +  40] = -4.8038446141526142E-001; yp[counter +  40] = -7.3379938570534275E-001; zp[counter +  40] = -4.8038446141526142E-001; wp[counter +  40] = w;
                xp[counter +  41] =  4.8038446141526142E-001; yp[counter +  41] =  7.3379938570534275E-001; zp[counter +  41] =  4.8038446141526142E-001; wp[counter +  41] = w;
                xp[counter +  42] =  7.3379938570534275E-001; yp[counter +  42] =  4.8038446141526142E-001; zp[counter +  42] =  4.8038446141526142E-001; wp[counter +  42] = w;
                xp[counter +  43] = -7.3379938570534275E-001; yp[counter +  43] =  4.8038446141526142E-001; zp[counter +  43] =  4.8038446141526142E-001; wp[counter +  43] = w;
                xp[counter +  44] =  7.3379938570534275E-001; yp[counter +  44] = -4.8038446141526142E-001; zp[counter +  44] =  4.8038446141526142E-001; wp[counter +  44] = w;
                xp[counter +  45] =  7.3379938570534275E-001; yp[counter +  45] =  4.8038446141526142E-001; zp[counter +  45] = -4.8038446141526142E-001; wp[counter +  45] = w;
                xp[counter +  46] = -7.3379938570534275E-001; yp[counter +  46] = -4.8038446141526142E-001; zp[counter +  46] =  4.8038446141526142E-001; wp[counter +  46] = w;
                xp[counter +  47] = -7.3379938570534275E-001; yp[counter +  47] =  4.8038446141526142E-001; zp[counter +  47] = -4.8038446141526142E-001; wp[counter +  47] = w;
                xp[counter +  48] =  7.3379938570534275E-001; yp[counter +  48] = -4.8038446141526142E-001; zp[counter +  48] = -4.8038446141526142E-001; wp[counter +  48] = w;
                xp[counter +  49] = -7.3379938570534275E-001; yp[counter +  49] = -4.8038446141526142E-001; zp[counter +  49] = -4.8038446141526142E-001; wp[counter +  49] = w;
                xp[counter +  50] =  3.2077264898077640E-001; yp[counter +  50] =  9.4715622136258792E-001; zp[counter +  50] =  0.0000000000000000E+000; wp[counter +  50] = w;
                xp[counter +  51] = -3.2077264898077640E-001; yp[counter +  51] =  9.4715622136258792E-001; zp[counter +  51] =  0.0000000000000000E+000; wp[counter +  51] = w;
                xp[counter +  52] =  3.2077264898077640E-001; yp[counter +  52] = -9.4715622136258792E-001; zp[counter +  52] =  0.0000000000000000E+000; wp[counter +  52] = w;
                xp[counter +  53] = -3.2077264898077640E-001; yp[counter +  53] = -9.4715622136258792E-001; zp[counter +  53] =  0.0000000000000000E+000; wp[counter +  53] = w;
                xp[counter +  54] =  9.4715622136258792E-001; yp[counter +  54] =  3.2077264898077640E-001; zp[counter +  54] =  0.0000000000000000E+000; wp[counter +  54] = w;
                xp[counter +  55] = -9.4715622136258792E-001; yp[counter +  55] =  3.2077264898077640E-001; zp[counter +  55] =  0.0000000000000000E+000; wp[counter +  55] = w;
                xp[counter +  56] =  9.4715622136258792E-001; yp[counter +  56] = -3.2077264898077640E-001; zp[counter +  56] =  0.0000000000000000E+000; wp[counter +  56] = w;
                xp[counter +  57] = -9.4715622136258792E-001; yp[counter +  57] = -3.2077264898077640E-001; zp[counter +  57] =  0.0000000000000000E+000; wp[counter +  57] = w;
                xp[counter +  58] =  3.2077264898077640E-001; yp[counter +  58] =  0.0000000000000000E+000; zp[counter +  58] =  9.4715622136258792E-001; wp[counter +  58] = w;
                xp[counter +  59] = -3.2077264898077640E-001; yp[counter +  59] =  0.0000000000000000E+000; zp[counter +  59] =  9.4715622136258792E-001; wp[counter +  59] = w;
                xp[counter +  60] =  3.2077264898077640E-001; yp[counter +  60] =  0.0000000000000000E+000; zp[counter +  60] = -9.4715622136258792E-001; wp[counter +  60] = w;
                xp[counter +  61] = -3.2077264898077640E-001; yp[counter +  61] =  0.0000000000000000E+000; zp[counter +  61] = -9.4715622136258792E-001; wp[counter +  61] = w;
                xp[counter +  62] =  9.4715622136258792E-001; yp[counter +  62] =  0.0000000000000000E+000; zp[counter +  62] =  3.2077264898077640E-001; wp[counter +  62] = w;
                xp[counter +  63] = -9.4715622136258792E-001; yp[counter +  63] =  0.0000000000000000E+000; zp[counter +  63] =  3.2077264898077640E-001; wp[counter +  63] = w;
                xp[counter +  64] =  9.4715622136258792E-001; yp[counter +  64] =  0.0000000000000000E+000; zp[counter +  64] = -3.2077264898077640E-001; wp[counter +  64] = w;
                xp[counter +  65] = -9.4715622136258792E-001; yp[counter +  65] =  0.0000000000000000E+000; zp[counter +  65] = -3.2077264898077640E-001; wp[counter +  65] = w;
                xp[counter +  66] =  0.0000000000000000E+000; yp[counter +  66] =  3.2077264898077640E-001; zp[counter +  66] =  9.4715622136258792E-001; wp[counter +  66] = w;
                xp[counter +  67] =  0.0000000000000000E+000; yp[counter +  67] = -3.2077264898077640E-001; zp[counter +  67] =  9.4715622136258792E-001; wp[counter +  67] = w;
                xp[counter +  68] =  0.0000000000000000E+000; yp[counter +  68] =  3.2077264898077640E-001; zp[counter +  68] = -9.4715622136258792E-001; wp[counter +  68] = w;
                xp[counter +  69] =  0.0000000000000000E+000; yp[counter +  69] = -3.2077264898077640E-001; zp[counter +  69] = -9.4715622136258792E-001; wp[counter +  69] = w;
                xp[counter +  70] =  0.0000000000000000E+000; yp[counter +  70] =  9.4715622136258792E-001; zp[counter +  70] =  3.2077264898077640E-001; wp[counter +  70] = w;
                xp[counter +  71] =  0.0000000000000000E+000; yp[counter +  71] = -9.4715622136258792E-001; zp[counter +  71] =  3.2077264898077640E-001; wp[counter +  71] = w;
                xp[counter +  72] =  0.0000000000000000E+000; yp[counter +  72] =  9.4715622136258792E-001; zp[counter +  72] = -3.2077264898077640E-001; wp[counter +  72] = w;
                xp[counter +  73] =  0.0000000000000000E+000; yp[counter +  73] = -9.4715622136258792E-001; zp[counter +  73] = -3.2077264898077640E-001; wp[counter +  73] = w;
            } else if (L == 86) {
                xp[counter +   0] =  1.0000000000000000E+000; yp[counter +   0] =  0.0000000000000000E+000; zp[counter +   0] =  0.0000000000000000E+000; wp[counter +   0] = w;
                xp[counter +   1] = -1.0000000000000000E+000; yp[counter +   1] =  0.0000000000000000E+000; zp[counter +   1] =  0.0000000000000000E+000; wp[counter +   1] = w;
                xp[counter +   2] =  0.0000000000000000E+000; yp[counter +   2] =  1.0000000000000000E+000; zp[counter +   2] =  0.0000000000000000E+000; wp[counter +   2] = w;
                xp[counter +   3] =  0.0000000000000000E+000; yp[counter +   3] = -1.0000000000000000E+000; zp[counter +   3] =  0.0000000000000000E+000; wp[counter +   3] = w;
                xp[counter +   4] =  0.0000000000000000E+000; yp[counter +   4] =  0.0000000000000000E+000; zp[counter +   4] =  1.0000000000000000E+000; wp[counter +   4] = w;
                xp[counter +   5] =  0.0000000000000000E+000; yp[counter +   5] =  0.0000000000000000E+000; zp[counter +   5] = -1.0000000000000000E+000; wp[counter +   5] = w;
                xp[counter +   6] =  5.7735026918962573E-001; yp[counter +   6] =  5.7735026918962573E-001; zp[counter +   6] =  5.7735026918962573E-001; wp[counter +   6] = w;
                xp[counter +   7] = -5.7735026918962573E-001; yp[counter +   7] =  5.7735026918962573E-001; zp[counter +   7] =  5.7735026918962573E-001; wp[counter +   7] = w;
                xp[counter +   8] =  5.7735026918962573E-001; yp[counter +   8] = -5.7735026918962573E-001; zp[counter +   8] =  5.7735026918962573E-001; wp[counter +   8] = w;
                xp[counter +   9] =  5.7735026918962573E-001; yp[counter +   9] =  5.7735026918962573E-001; zp[counter +   9] = -5.7735026918962573E-001; wp[counter +   9] = w;
                xp[counter +  10] = -5.7735026918962573E-001; yp[counter +  10] = -5.7735026918962573E-001; zp[counter +  10] =  5.7735026918962573E-001; wp[counter +  10] = w;
                xp[counter +  11] =  5.7735026918962573E-001; yp[counter +  11] = -5.7735026918962573E-001; zp[counter +  11] = -5.7735026918962573E-001; wp[counter +  11] = w;
                xp[counter +  12] = -5.7735026918962573E-001; yp[counter +  12] =  5.7735026918962573E-001; zp[counter +  12] = -5.7735026918962573E-001; wp[counter +  12] = w;
                xp[counter +  13] = -5.7735026918962573E-001; yp[counter +  13] = -5.7735026918962573E-001; zp[counter +  13] = -5.7735026918962573E-001; wp[counter +  13] = w;
                xp[counter +  14] =  3.6960284645415020E-001; yp[counter +  14] =  3.6960284645415020E-001; zp[counter +  14] =  8.5251831170126757E-001; wp[counter +  14] = w;
                xp[counter +  15] = -3.6960284645415020E-001; yp[counter +  15] =  3.6960284645415020E-001; zp[counter +  15] =  8.5251831170126757E-001; wp[counter +  15] = w;
                xp[counter +  16] =  3.6960284645415020E-001; yp[counter +  16] = -3.6960284645415020E-001; zp[counter +  16] =  8.5251831170126757E-001; wp[counter +  16] = w;
                xp[counter +  17] =  3.6960284645415020E-001; yp[counter +  17] =  3.6960284645415020E-001; zp[counter +  17] = -8.5251831170126757E-001; wp[counter +  17] = w;
                xp[counter +  18] = -3.6960284645415020E-001; yp[counter +  18] = -3.6960284645415020E-001; zp[counter +  18] =  8.5251831170126757E-001; wp[counter +  18] = w;
                xp[counter +  19] = -3.6960284645415020E-001; yp[counter +  19] =  3.6960284645415020E-001; zp[counter +  19] = -8.5251831170126757E-001; wp[counter +  19] = w;
                xp[counter +  20] =  3.6960284645415020E-001; yp[counter +  20] = -3.6960284645415020E-001; zp[counter +  20] = -8.5251831170126757E-001; wp[counter +  20] = w;
                xp[counter +  21] = -3.6960284645415020E-001; yp[counter +  21] = -3.6960284645415020E-001; zp[counter +  21] = -8.5251831170126757E-001; wp[counter +  21] = w;
                xp[counter +  22] = -3.6960284645415020E-001; yp[counter +  22] =  8.5251831170126757E-001; zp[counter +  22] =  3.6960284645415020E-001; wp[counter +  22] = w;
                xp[counter +  23] =  3.6960284645415020E-001; yp[counter +  23] = -8.5251831170126757E-001; zp[counter +  23] =  3.6960284645415020E-001; wp[counter +  23] = w;
                xp[counter +  24] =  3.6960284645415020E-001; yp[counter +  24] =  8.5251831170126757E-001; zp[counter +  24] = -3.6960284645415020E-001; wp[counter +  24] = w;
                xp[counter +  25] = -3.6960284645415020E-001; yp[counter +  25] = -8.5251831170126757E-001; zp[counter +  25] =  3.6960284645415020E-001; wp[counter +  25] = w;
                xp[counter +  26] = -3.6960284645415020E-001; yp[counter +  26] =  8.5251831170126757E-001; zp[counter +  26] = -3.6960284645415020E-001; wp[counter +  26] = w;
                xp[counter +  27] =  3.6960284645415020E-001; yp[counter +  27] = -8.5251831170126757E-001; zp[counter +  27] = -3.6960284645415020E-001; wp[counter +  27] = w;
                xp[counter +  28] = -3.6960284645415020E-001; yp[counter +  28] = -8.5251831170126757E-001; zp[counter +  28] = -3.6960284645415020E-001; wp[counter +  28] = w;
                xp[counter +  29] =  3.6960284645415020E-001; yp[counter +  29] =  8.5251831170126757E-001; zp[counter +  29] =  3.6960284645415020E-001; wp[counter +  29] = w;
                xp[counter +  30] =  8.5251831170126757E-001; yp[counter +  30] =  3.6960284645415020E-001; zp[counter +  30] =  3.6960284645415020E-001; wp[counter +  30] = w;
                xp[counter +  31] = -8.5251831170126757E-001; yp[counter +  31] =  3.6960284645415020E-001; zp[counter +  31] =  3.6960284645415020E-001; wp[counter +  31] = w;
                xp[counter +  32] =  8.5251831170126757E-001; yp[counter +  32] = -3.6960284645415020E-001; zp[counter +  32] =  3.6960284645415020E-001; wp[counter +  32] = w;
                xp[counter +  33] =  8.5251831170126757E-001; yp[counter +  33] =  3.6960284645415020E-001; zp[counter +  33] = -3.6960284645415020E-001; wp[counter +  33] = w;
                xp[counter +  34] = -8.5251831170126757E-001; yp[counter +  34] = -3.6960284645415020E-001; zp[counter +  34] =  3.6960284645415020E-001; wp[counter +  34] = w;
                xp[counter +  35] = -8.5251831170126757E-001; yp[counter +  35] =  3.6960284645415020E-001; zp[counter +  35] = -3.6960284645415020E-001; wp[counter +  35] = w;
                xp[counter +  36] =  8.5251831170126757E-001; yp[counter +  36] = -3.6960284645415020E-001; zp[counter +  36] = -3.6960284645415020E-001; wp[counter +  36] = w;
                xp[counter +  37] = -8.5251831170126757E-001; yp[counter +  37] = -3.6960284645415020E-001; zp[counter +  37] = -3.6960284645415020E-001; wp[counter +  37] = w;
                xp[counter +  38] =  6.9435400660266644E-001; yp[counter +  38] =  6.9435400660266644E-001; zp[counter +  38] =  1.8906355288539498E-001; wp[counter +  38] = w;
                xp[counter +  39] = -6.9435400660266644E-001; yp[counter +  39] =  6.9435400660266644E-001; zp[counter +  39] =  1.8906355288539498E-001; wp[counter +  39] = w;
                xp[counter +  40] =  6.9435400660266644E-001; yp[counter +  40] = -6.9435400660266644E-001; zp[counter +  40] =  1.8906355288539498E-001; wp[counter +  40] = w;
                xp[counter +  41] =  6.9435400660266644E-001; yp[counter +  41] =  6.9435400660266644E-001; zp[counter +  41] = -1.8906355288539498E-001; wp[counter +  41] = w;
                xp[counter +  42] = -6.9435400660266644E-001; yp[counter +  42] = -6.9435400660266644E-001; zp[counter +  42] =  1.8906355288539498E-001; wp[counter +  42] = w;
                xp[counter +  43] = -6.9435400660266644E-001; yp[counter +  43] =  6.9435400660266644E-001; zp[counter +  43] = -1.8906355288539498E-001; wp[counter +  43] = w;
                xp[counter +  44] =  6.9435400660266644E-001; yp[counter +  44] = -6.9435400660266644E-001; zp[counter +  44] = -1.8906355288539498E-001; wp[counter +  44] = w;
                xp[counter +  45] = -6.9435400660266644E-001; yp[counter +  45] = -6.9435400660266644E-001; zp[counter +  45] = -1.8906355288539498E-001; wp[counter +  45] = w;
                xp[counter +  46] = -6.9435400660266644E-001; yp[counter +  46] =  1.8906355288539498E-001; zp[counter +  46] =  6.9435400660266644E-001; wp[counter +  46] = w;
                xp[counter +  47] =  6.9435400660266644E-001; yp[counter +  47] = -1.8906355288539498E-001; zp[counter +  47] =  6.9435400660266644E-001; wp[counter +  47] = w;
                xp[counter +  48] =  6.9435400660266644E-001; yp[counter +  48] =  1.8906355288539498E-001; zp[counter +  48] = -6.9435400660266644E-001; wp[counter +  48] = w;
                xp[counter +  49] = -6.9435400660266644E-001; yp[counter +  49] = -1.8906355288539498E-001; zp[counter +  49] =  6.9435400660266644E-001; wp[counter +  49] = w;
                xp[counter +  50] = -6.9435400660266644E-001; yp[counter +  50] =  1.8906355288539498E-001; zp[counter +  50] = -6.9435400660266644E-001; wp[counter +  50] = w;
                xp[counter +  51] =  6.9435400660266644E-001; yp[counter +  51] = -1.8906355288539498E-001; zp[counter +  51] = -6.9435400660266644E-001; wp[counter +  51] = w;
                xp[counter +  52] = -6.9435400660266644E-001; yp[counter +  52] = -1.8906355288539498E-001; zp[counter +  52] = -6.9435400660266644E-001; wp[counter +  52] = w;
                xp[counter +  53] =  6.9435400660266644E-001; yp[counter +  53] =  1.8906355288539498E-001; zp[counter +  53] =  6.9435400660266644E-001; wp[counter +  53] = w;
                xp[counter +  54] =  1.8906355288539498E-001; yp[counter +  54] =  6.9435400660266644E-001; zp[counter +  54] =  6.9435400660266644E-001; wp[counter +  54] = w;
                xp[counter +  55] = -1.8906355288539498E-001; yp[counter +  55] =  6.9435400660266644E-001; zp[counter +  55] =  6.9435400660266644E-001; wp[counter +  55] = w;
                xp[counter +  56] =  1.8906355288539498E-001; yp[counter +  56] = -6.9435400660266644E-001; zp[counter +  56] =  6.9435400660266644E-001; wp[counter +  56] = w;
                xp[counter +  57] =  1.8906355288539498E-001; yp[counter +  57] =  6.9435400660266644E-001; zp[counter +  57] = -6.9435400660266644E-001; wp[counter +  57] = w;
                xp[counter +  58] = -1.8906355288539498E-001; yp[counter +  58] = -6.9435400660266644E-001; zp[counter +  58] =  6.9435400660266644E-001; wp[counter +  58] = w;
                xp[counter +  59] = -1.8906355288539498E-001; yp[counter +  59] =  6.9435400660266644E-001; zp[counter +  59] = -6.9435400660266644E-001; wp[counter +  59] = w;
                xp[counter +  60] =  1.8906355288539498E-001; yp[counter +  60] = -6.9435400660266644E-001; zp[counter +  60] = -6.9435400660266644E-001; wp[counter +  60] = w;
                xp[counter +  61] = -1.8906355288539498E-001; yp[counter +  61] = -6.9435400660266644E-001; zp[counter +  61] = -6.9435400660266644E-001; wp[counter +  61] = w;
                xp[counter +  62] =  3.7424303909034118E-001; yp[counter +  62] =  9.2733065715117247E-001; zp[counter +  62] =  0.0000000000000000E+000; wp[counter +  62] = w;
                xp[counter +  63] = -3.7424303909034118E-001; yp[counter +  63] =  9.2733065715117247E-001; zp[counter +  63] =  0.0000000000000000E+000; wp[counter +  63] = w;
                xp[counter +  64] =  3.7424303909034118E-001; yp[counter +  64] = -9.2733065715117247E-001; zp[counter +  64] =  0.0000000000000000E+000; wp[counter +  64] = w;
                xp[counter +  65] = -3.7424303909034118E-001; yp[counter +  65] = -9.2733065715117247E-001; zp[counter +  65] =  0.0000000000000000E+000; wp[counter +  65] = w;
                xp[counter +  66] =  9.2733065715117247E-001; yp[counter +  66] =  3.7424303909034118E-001; zp[counter +  66] =  0.0000000000000000E+000; wp[counter +  66] = w;
                xp[counter +  67] = -9.2733065715117247E-001; yp[counter +  67] =  3.7424303909034118E-001; zp[counter +  67] =  0.0000000000000000E+000; wp[counter +  67] = w;
                xp[counter +  68] =  9.2733065715117247E-001; yp[counter +  68] = -3.7424303909034118E-001; zp[counter +  68] =  0.0000000000000000E+000; wp[counter +  68] = w;
                xp[counter +  69] = -9.2733065715117247E-001; yp[counter +  69] = -3.7424303909034118E-001; zp[counter +  69] =  0.0000000000000000E+000; wp[counter +  69] = w;
                xp[counter +  70] =  3.7424303909034118E-001; yp[counter +  70] =  0.0000000000000000E+000; zp[counter +  70] =  9.2733065715117247E-001; wp[counter +  70] = w;
                xp[counter +  71] = -3.7424303909034118E-001; yp[counter +  71] =  0.0000000000000000E+000; zp[counter +  71] =  9.2733065715117247E-001; wp[counter +  71] = w;
                xp[counter +  72] =  3.7424303909034118E-001; yp[counter +  72] =  0.0000000000000000E+000; zp[counter +  72] = -9.2733065715117247E-001; wp[counter +  72] = w;
                xp[counter +  73] = -3.7424303909034118E-001; yp[counter +  73] =  0.0000000000000000E+000; zp[counter +  73] = -9.2733065715117247E-001; wp[counter +  73] = w;
                xp[counter +  74] =  9.2733065715117247E-001; yp[counter +  74] =  0.0000000000000000E+000; zp[counter +  74] =  3.7424303909034118E-001; wp[counter +  74] = w;
                xp[counter +  75] = -9.2733065715117247E-001; yp[counter +  75] =  0.0000000000000000E+000; zp[counter +  75] =  3.7424303909034118E-001; wp[counter +  75] = w;
                xp[counter +  76] =  9.2733065715117247E-001; yp[counter +  76] =  0.0000000000000000E+000; zp[counter +  76] = -3.7424303909034118E-001; wp[counter +  76] = w;
                xp[counter +  77] = -9.2733065715117247E-001; yp[counter +  77] =  0.0000000000000000E+000; zp[counter +  77] = -3.7424303909034118E-001; wp[counter +  77] = w;
                xp[counter +  78] =  0.0000000000000000E+000; yp[counter +  78] =  3.7424303909034118E-001; zp[counter +  78] =  9.2733065715117247E-001; wp[counter +  78] = w;
                xp[counter +  79] =  0.0000000000000000E+000; yp[counter +  79] = -3.7424303909034118E-001; zp[counter +  79] =  9.2733065715117247E-001; wp[counter +  79] = w;
                xp[counter +  80] =  0.0000000000000000E+000; yp[counter +  80] =  3.7424303909034118E-001; zp[counter +  80] = -9.2733065715117247E-001; wp[counter +  80] = w;
                xp[counter +  81] =  0.0000000000000000E+000; yp[counter +  81] = -3.7424303909034118E-001; zp[counter +  81] = -9.2733065715117247E-001; wp[counter +  81] = w;
                xp[counter +  82] =  0.0000000000000000E+000; yp[counter +  82] =  9.2733065715117247E-001; zp[counter +  82] =  3.7424303909034118E-001; wp[counter +  82] = w;
                xp[counter +  83] =  0.0000000000000000E+000; yp[counter +  83] = -9.2733065715117247E-001; zp[counter +  83] =  3.7424303909034118E-001; wp[counter +  83] = w;
                xp[counter +  84] =  0.0000000000000000E+000; yp[counter +  84] =  9.2733065715117247E-001; zp[counter +  84] = -3.7424303909034118E-001; wp[counter +  84] = w;
                xp[counter +  85] =  0.0000000000000000E+000; yp[counter +  85] = -9.2733065715117247E-001; zp[counter +  85] = -3.7424303909034118E-001; wp[counter +  85] = w;
            } else if (L == 110) {
                xp[counter +   0] =  1.0000000000000000E+000; yp[counter +   0] =  0.0000000000000000E+000; zp[counter +   0] =  0.0000000000000000E+000; wp[counter +   0] = w;
                xp[counter +   1] = -1.0000000000000000E+000; yp[counter +   1] =  0.0000000000000000E+000; zp[counter +   1] =  0.0000000000000000E+000; wp[counter +   1] = w;
                xp[counter +   2] =  0.0000000000000000E+000; yp[counter +   2] =  1.0000000000000000E+000; zp[counter +   2] =  0.0000000000000000E+000; wp[counter +   2] = w;
                xp[counter +   3] =  0.0000000000000000E+000; yp[counter +   3] = -1.0000000000000000E+000; zp[counter +   3] =  0.0000000000000000E+000; wp[counter +   3] = w;
                xp[counter +   4] =  0.0000000000000000E+000; yp[counter +   4] =  0.0000000000000000E+000; zp[counter +   4] =  1.0000000000000000E+000; wp[counter +   4] = w;
                xp[counter +   5] =  0.0000000000000000E+000; yp[counter +   5] =  0.0000000000000000E+000; zp[counter +   5] = -1.0000000000000000E+000; wp[counter +   5] = w;
                xp[counter +   6] =  5.7735026918962573E-001; yp[counter +   6] =  5.7735026918962573E-001; zp[counter +   6] =  5.7735026918962573E-001; wp[counter +   6] = w;
                xp[counter +   7] = -5.7735026918962573E-001; yp[counter +   7] =  5.7735026918962573E-001; zp[counter +   7] =  5.7735026918962573E-001; wp[counter +   7] = w;
                xp[counter +   8] =  5.7735026918962573E-001; yp[counter +   8] = -5.7735026918962573E-001; zp[counter +   8] =  5.7735026918962573E-001; wp[counter +   8] = w;
                xp[counter +   9] =  5.7735026918962573E-001; yp[counter +   9] =  5.7735026918962573E-001; zp[counter +   9] = -5.7735026918962573E-001; wp[counter +   9] = w;
                xp[counter +  10] = -5.7735026918962573E-001; yp[counter +  10] = -5.7735026918962573E-001; zp[counter +  10] =  5.7735026918962573E-001; wp[counter +  10] = w;
                xp[counter +  11] =  5.7735026918962573E-001; yp[counter +  11] = -5.7735026918962573E-001; zp[counter +  11] = -5.7735026918962573E-001; wp[counter +  11] = w;
                xp[counter +  12] = -5.7735026918962573E-001; yp[counter +  12] =  5.7735026918962573E-001; zp[counter +  12] = -5.7735026918962573E-001; wp[counter +  12] = w;
                xp[counter +  13] = -5.7735026918962573E-001; yp[counter +  13] = -5.7735026918962573E-001; zp[counter +  13] = -5.7735026918962573E-001; wp[counter +  13] = w;
                xp[counter +  14] =  1.8511563534473621E-001; yp[counter +  14] =  1.8511563534473621E-001; zp[counter +  14] =  9.6512403508659406E-001; wp[counter +  14] = w;
                xp[counter +  15] = -1.8511563534473621E-001; yp[counter +  15] =  1.8511563534473621E-001; zp[counter +  15] =  9.6512403508659406E-001; wp[counter +  15] = w;
                xp[counter +  16] =  1.8511563534473621E-001; yp[counter +  16] = -1.8511563534473621E-001; zp[counter +  16] =  9.6512403508659406E-001; wp[counter +  16] = w;
                xp[counter +  17] =  1.8511563534473621E-001; yp[counter +  17] =  1.8511563534473621E-001; zp[counter +  17] = -9.6512403508659406E-001; wp[counter +  17] = w;
                xp[counter +  18] = -1.8511563534473621E-001; yp[counter +  18] = -1.8511563534473621E-001; zp[counter +  18] =  9.6512403508659406E-001; wp[counter +  18] = w;
                xp[counter +  19] = -1.8511563534473621E-001; yp[counter +  19] =  1.8511563534473621E-001; zp[counter +  19] = -9.6512403508659406E-001; wp[counter +  19] = w;
                xp[counter +  20] =  1.8511563534473621E-001; yp[counter +  20] = -1.8511563534473621E-001; zp[counter +  20] = -9.6512403508659406E-001; wp[counter +  20] = w;
                xp[counter +  21] = -1.8511563534473621E-001; yp[counter +  21] = -1.8511563534473621E-001; zp[counter +  21] = -9.6512403508659406E-001; wp[counter +  21] = w;
                xp[counter +  22] = -1.8511563534473621E-001; yp[counter +  22] =  9.6512403508659406E-001; zp[counter +  22] =  1.8511563534473621E-001; wp[counter +  22] = w;
                xp[counter +  23] =  1.8511563534473621E-001; yp[counter +  23] = -9.6512403508659406E-001; zp[counter +  23] =  1.8511563534473621E-001; wp[counter +  23] = w;
                xp[counter +  24] =  1.8511563534473621E-001; yp[counter +  24] =  9.6512403508659406E-001; zp[counter +  24] = -1.8511563534473621E-001; wp[counter +  24] = w;
                xp[counter +  25] = -1.8511563534473621E-001; yp[counter +  25] = -9.6512403508659406E-001; zp[counter +  25] =  1.8511563534473621E-001; wp[counter +  25] = w;
                xp[counter +  26] = -1.8511563534473621E-001; yp[counter +  26] =  9.6512403508659406E-001; zp[counter +  26] = -1.8511563534473621E-001; wp[counter +  26] = w;
                xp[counter +  27] =  1.8511563534473621E-001; yp[counter +  27] = -9.6512403508659406E-001; zp[counter +  27] = -1.8511563534473621E-001; wp[counter +  27] = w;
                xp[counter +  28] = -1.8511563534473621E-001; yp[counter +  28] = -9.6512403508659406E-001; zp[counter +  28] = -1.8511563534473621E-001; wp[counter +  28] = w;
                xp[counter +  29] =  1.8511563534473621E-001; yp[counter +  29] =  9.6512403508659406E-001; zp[counter +  29] =  1.8511563534473621E-001; wp[counter +  29] = w;
                xp[counter +  30] =  9.6512403508659406E-001; yp[counter +  30] =  1.8511563534473621E-001; zp[counter +  30] =  1.8511563534473621E-001; wp[counter +  30] = w;
                xp[counter +  31] = -9.6512403508659406E-001; yp[counter +  31] =  1.8511563534473621E-001; zp[counter +  31] =  1.8511563534473621E-001; wp[counter +  31] = w;
                xp[counter +  32] =  9.6512403508659406E-001; yp[counter +  32] = -1.8511563534473621E-001; zp[counter +  32] =  1.8511563534473621E-001; wp[counter +  32] = w;
                xp[counter +  33] =  9.6512403508659406E-001; yp[counter +  33] =  1.8511563534473621E-001; zp[counter +  33] = -1.8511563534473621E-001; wp[counter +  33] = w;
                xp[counter +  34] = -9.6512403508659406E-001; yp[counter +  34] = -1.8511563534473621E-001; zp[counter +  34] =  1.8511563534473621E-001; wp[counter +  34] = w;
                xp[counter +  35] = -9.6512403508659406E-001; yp[counter +  35] =  1.8511563534473621E-001; zp[counter +  35] = -1.8511563534473621E-001; wp[counter +  35] = w;
                xp[counter +  36] =  9.6512403508659406E-001; yp[counter +  36] = -1.8511563534473621E-001; zp[counter +  36] = -1.8511563534473621E-001; wp[counter +  36] = w;
                xp[counter +  37] = -9.6512403508659406E-001; yp[counter +  37] = -1.8511563534473621E-001; zp[counter +  37] = -1.8511563534473621E-001; wp[counter +  37] = w;
                xp[counter +  38] =  6.9042104838229224E-001; yp[counter +  38] =  6.9042104838229224E-001; zp[counter +  38] =  2.1595729184584844E-001; wp[counter +  38] = w;
                xp[counter +  39] = -6.9042104838229224E-001; yp[counter +  39] =  6.9042104838229224E-001; zp[counter +  39] =  2.1595729184584844E-001; wp[counter +  39] = w;
                xp[counter +  40] =  6.9042104838229224E-001; yp[counter +  40] = -6.9042104838229224E-001; zp[counter +  40] =  2.1595729184584844E-001; wp[counter +  40] = w;
                xp[counter +  41] =  6.9042104838229224E-001; yp[counter +  41] =  6.9042104838229224E-001; zp[counter +  41] = -2.1595729184584844E-001; wp[counter +  41] = w;
                xp[counter +  42] = -6.9042104838229224E-001; yp[counter +  42] = -6.9042104838229224E-001; zp[counter +  42] =  2.1595729184584844E-001; wp[counter +  42] = w;
                xp[counter +  43] = -6.9042104838229224E-001; yp[counter +  43] =  6.9042104838229224E-001; zp[counter +  43] = -2.1595729184584844E-001; wp[counter +  43] = w;
                xp[counter +  44] =  6.9042104838229224E-001; yp[counter +  44] = -6.9042104838229224E-001; zp[counter +  44] = -2.1595729184584844E-001; wp[counter +  44] = w;
                xp[counter +  45] = -6.9042104838229224E-001; yp[counter +  45] = -6.9042104838229224E-001; zp[counter +  45] = -2.1595729184584844E-001; wp[counter +  45] = w;
                xp[counter +  46] = -6.9042104838229224E-001; yp[counter +  46] =  2.1595729184584844E-001; zp[counter +  46] =  6.9042104838229224E-001; wp[counter +  46] = w;
                xp[counter +  47] =  6.9042104838229224E-001; yp[counter +  47] = -2.1595729184584844E-001; zp[counter +  47] =  6.9042104838229224E-001; wp[counter +  47] = w;
                xp[counter +  48] =  6.9042104838229224E-001; yp[counter +  48] =  2.1595729184584844E-001; zp[counter +  48] = -6.9042104838229224E-001; wp[counter +  48] = w;
                xp[counter +  49] = -6.9042104838229224E-001; yp[counter +  49] = -2.1595729184584844E-001; zp[counter +  49] =  6.9042104838229224E-001; wp[counter +  49] = w;
                xp[counter +  50] = -6.9042104838229224E-001; yp[counter +  50] =  2.1595729184584844E-001; zp[counter +  50] = -6.9042104838229224E-001; wp[counter +  50] = w;
                xp[counter +  51] =  6.9042104838229224E-001; yp[counter +  51] = -2.1595729184584844E-001; zp[counter +  51] = -6.9042104838229224E-001; wp[counter +  51] = w;
                xp[counter +  52] = -6.9042104838229224E-001; yp[counter +  52] = -2.1595729184584844E-001; zp[counter +  52] = -6.9042104838229224E-001; wp[counter +  52] = w;
                xp[counter +  53] =  6.9042104838229224E-001; yp[counter +  53] =  2.1595729184584844E-001; zp[counter +  53] =  6.9042104838229224E-001; wp[counter +  53] = w;
                xp[counter +  54] =  2.1595729184584844E-001; yp[counter +  54] =  6.9042104838229224E-001; zp[counter +  54] =  6.9042104838229224E-001; wp[counter +  54] = w;
                xp[counter +  55] = -2.1595729184584844E-001; yp[counter +  55] =  6.9042104838229224E-001; zp[counter +  55] =  6.9042104838229224E-001; wp[counter +  55] = w;
                xp[counter +  56] =  2.1595729184584844E-001; yp[counter +  56] = -6.9042104838229224E-001; zp[counter +  56] =  6.9042104838229224E-001; wp[counter +  56] = w;
                xp[counter +  57] =  2.1595729184584844E-001; yp[counter +  57] =  6.9042104838229224E-001; zp[counter +  57] = -6.9042104838229224E-001; wp[counter +  57] = w;
                xp[counter +  58] = -2.1595729184584844E-001; yp[counter +  58] = -6.9042104838229224E-001; zp[counter +  58] =  6.9042104838229224E-001; wp[counter +  58] = w;
                xp[counter +  59] = -2.1595729184584844E-001; yp[counter +  59] =  6.9042104838229224E-001; zp[counter +  59] = -6.9042104838229224E-001; wp[counter +  59] = w;
                xp[counter +  60] =  2.1595729184584844E-001; yp[counter +  60] = -6.9042104838229224E-001; zp[counter +  60] = -6.9042104838229224E-001; wp[counter +  60] = w;
                xp[counter +  61] = -2.1595729184584844E-001; yp[counter +  61] = -6.9042104838229224E-001; zp[counter +  61] = -6.9042104838229224E-001; wp[counter +  61] = w;
                xp[counter +  62] =  3.9568947305594188E-001; yp[counter +  62] =  3.9568947305594188E-001; zp[counter +  62] =  8.2876998125259227E-001; wp[counter +  62] = w;
                xp[counter +  63] = -3.9568947305594188E-001; yp[counter +  63] =  3.9568947305594188E-001; zp[counter +  63] =  8.2876998125259227E-001; wp[counter +  63] = w;
                xp[counter +  64] =  3.9568947305594188E-001; yp[counter +  64] = -3.9568947305594188E-001; zp[counter +  64] =  8.2876998125259227E-001; wp[counter +  64] = w;
                xp[counter +  65] =  3.9568947305594188E-001; yp[counter +  65] =  3.9568947305594188E-001; zp[counter +  65] = -8.2876998125259227E-001; wp[counter +  65] = w;
                xp[counter +  66] = -3.9568947305594188E-001; yp[counter +  66] = -3.9568947305594188E-001; zp[counter +  66] =  8.2876998125259227E-001; wp[counter +  66] = w;
                xp[counter +  67] = -3.9568947305594188E-001; yp[counter +  67] =  3.9568947305594188E-001; zp[counter +  67] = -8.2876998125259227E-001; wp[counter +  67] = w;
                xp[counter +  68] =  3.9568947305594188E-001; yp[counter +  68] = -3.9568947305594188E-001; zp[counter +  68] = -8.2876998125259227E-001; wp[counter +  68] = w;
                xp[counter +  69] = -3.9568947305594188E-001; yp[counter +  69] = -3.9568947305594188E-001; zp[counter +  69] = -8.2876998125259227E-001; wp[counter +  69] = w;
                xp[counter +  70] = -3.9568947305594188E-001; yp[counter +  70] =  8.2876998125259227E-001; zp[counter +  70] =  3.9568947305594188E-001; wp[counter +  70] = w;
                xp[counter +  71] =  3.9568947305594188E-001; yp[counter +  71] = -8.2876998125259227E-001; zp[counter +  71] =  3.9568947305594188E-001; wp[counter +  71] = w;
                xp[counter +  72] =  3.9568947305594188E-001; yp[counter +  72] =  8.2876998125259227E-001; zp[counter +  72] = -3.9568947305594188E-001; wp[counter +  72] = w;
                xp[counter +  73] = -3.9568947305594188E-001; yp[counter +  73] = -8.2876998125259227E-001; zp[counter +  73] =  3.9568947305594188E-001; wp[counter +  73] = w;
                xp[counter +  74] = -3.9568947305594188E-001; yp[counter +  74] =  8.2876998125259227E-001; zp[counter +  74] = -3.9568947305594188E-001; wp[counter +  74] = w;
                xp[counter +  75] =  3.9568947305594188E-001; yp[counter +  75] = -8.2876998125259227E-001; zp[counter +  75] = -3.9568947305594188E-001; wp[counter +  75] = w;
                xp[counter +  76] = -3.9568947305594188E-001; yp[counter +  76] = -8.2876998125259227E-001; zp[counter +  76] = -3.9568947305594188E-001; wp[counter +  76] = w;
                xp[counter +  77] =  3.9568947305594188E-001; yp[counter +  77] =  8.2876998125259227E-001; zp[counter +  77] =  3.9568947305594188E-001; wp[counter +  77] = w;
                xp[counter +  78] =  8.2876998125259227E-001; yp[counter +  78] =  3.9568947305594188E-001; zp[counter +  78] =  3.9568947305594188E-001; wp[counter +  78] = w;
                xp[counter +  79] = -8.2876998125259227E-001; yp[counter +  79] =  3.9568947305594188E-001; zp[counter +  79] =  3.9568947305594188E-001; wp[counter +  79] = w;
                xp[counter +  80] =  8.2876998125259227E-001; yp[counter +  80] = -3.9568947305594188E-001; zp[counter +  80] =  3.9568947305594188E-001; wp[counter +  80] = w;
                xp[counter +  81] =  8.2876998125259227E-001; yp[counter +  81] =  3.9568947305594188E-001; zp[counter +  81] = -3.9568947305594188E-001; wp[counter +  81] = w;
                xp[counter +  82] = -8.2876998125259227E-001; yp[counter +  82] = -3.9568947305594188E-001; zp[counter +  82] =  3.9568947305594188E-001; wp[counter +  82] = w;
                xp[counter +  83] = -8.2876998125259227E-001; yp[counter +  83] =  3.9568947305594188E-001; zp[counter +  83] = -3.9568947305594188E-001; wp[counter +  83] = w;
                xp[counter +  84] =  8.2876998125259227E-001; yp[counter +  84] = -3.9568947305594188E-001; zp[counter +  84] = -3.9568947305594188E-001; wp[counter +  84] = w;
                xp[counter +  85] = -8.2876998125259227E-001; yp[counter +  85] = -3.9568947305594188E-001; zp[counter +  85] = -3.9568947305594188E-001; wp[counter +  85] = w;
                xp[counter +  86] =  4.7836902881215021E-001; yp[counter +  86] =  8.7815891060406615E-001; zp[counter +  86] =  0.0000000000000000E+000; wp[counter +  86] = w;
                xp[counter +  87] = -4.7836902881215021E-001; yp[counter +  87] =  8.7815891060406615E-001; zp[counter +  87] =  0.0000000000000000E+000; wp[counter +  87] = w;
                xp[counter +  88] =  4.7836902881215021E-001; yp[counter +  88] = -8.7815891060406615E-001; zp[counter +  88] =  0.0000000000000000E+000; wp[counter +  88] = w;
                xp[counter +  89] = -4.7836902881215021E-001; yp[counter +  89] = -8.7815891060406615E-001; zp[counter +  89] =  0.0000000000000000E+000; wp[counter +  89] = w;
                xp[counter +  90] =  8.7815891060406615E-001; yp[counter +  90] =  4.7836902881215021E-001; zp[counter +  90] =  0.0000000000000000E+000; wp[counter +  90] = w;
                xp[counter +  91] = -8.7815891060406615E-001; yp[counter +  91] =  4.7836902881215021E-001; zp[counter +  91] =  0.0000000000000000E+000; wp[counter +  91] = w;
                xp[counter +  92] =  8.7815891060406615E-001; yp[counter +  92] = -4.7836902881215021E-001; zp[counter +  92] =  0.0000000000000000E+000; wp[counter +  92] = w;
                xp[counter +  93] = -8.7815891060406615E-001; yp[counter +  93] = -4.7836902881215021E-001; zp[counter +  93] =  0.0000000000000000E+000; wp[counter +  93] = w;
                xp[counter +  94] =  4.7836902881215021E-001; yp[counter +  94] =  0.0000000000000000E+000; zp[counter +  94] =  8.7815891060406615E-001; wp[counter +  94] = w;
                xp[counter +  95] = -4.7836902881215021E-001; yp[counter +  95] =  0.0000000000000000E+000; zp[counter +  95] =  8.7815891060406615E-001; wp[counter +  95] = w;
                xp[counter +  96] =  4.7836902881215021E-001; yp[counter +  96] =  0.0000000000000000E+000; zp[counter +  96] = -8.7815891060406615E-001; wp[counter +  96] = w;
                xp[counter +  97] = -4.7836902881215021E-001; yp[counter +  97] =  0.0000000000000000E+000; zp[counter +  97] = -8.7815891060406615E-001; wp[counter +  97] = w;
                xp[counter +  98] =  8.7815891060406615E-001; yp[counter +  98] =  0.0000000000000000E+000; zp[counter +  98] =  4.7836902881215021E-001; wp[counter +  98] = w;
                xp[counter +  99] = -8.7815891060406615E-001; yp[counter +  99] =  0.0000000000000000E+000; zp[counter +  99] =  4.7836902881215021E-001; wp[counter +  99] = w;
                xp[counter + 100] =  8.7815891060406615E-001; yp[counter + 100] =  0.0000000000000000E+000; zp[counter + 100] = -4.7836902881215021E-001; wp[counter + 100] = w;
                xp[counter + 101] = -8.7815891060406615E-001; yp[counter + 101] =  0.0000000000000000E+000; zp[counter + 101] = -4.7836902881215021E-001; wp[counter + 101] = w;
                xp[counter + 102] =  0.0000000000000000E+000; yp[counter + 102] =  4.7836902881215021E-001; zp[counter + 102] =  8.7815891060406615E-001; wp[counter + 102] = w;
                xp[counter + 103] =  0.0000000000000000E+000; yp[counter + 103] = -4.7836902881215021E-001; zp[counter + 103] =  8.7815891060406615E-001; wp[counter + 103] = w;
                xp[counter + 104] =  0.0000000000000000E+000; yp[counter + 104] =  4.7836902881215021E-001; zp[counter + 104] = -8.7815891060406615E-001; wp[counter + 104] = w;
                xp[counter + 105] =  0.0000000000000000E+000; yp[counter + 105] = -4.7836902881215021E-001; zp[counter + 105] = -8.7815891060406615E-001; wp[counter + 105] = w;
                xp[counter + 106] =  0.0000000000000000E+000; yp[counter + 106] =  8.7815891060406615E-001; zp[counter + 106] =  4.7836902881215021E-001; wp[counter + 106] = w;
                xp[counter + 107] =  0.0000000000000000E+000; yp[counter + 107] = -8.7815891060406615E-001; zp[counter + 107] =  4.7836902881215021E-001; wp[counter + 107] = w;
                xp[counter + 108] =  0.0000000000000000E+000; yp[counter + 108] =  8.7815891060406615E-001; zp[counter + 108] = -4.7836902881215021E-001; wp[counter + 108] = w;
                xp[counter + 109] =  0.0000000000000000E+000; yp[counter + 109] = -8.7815891060406615E-001; zp[counter + 109] = -4.7836902881215021E-001; wp[counter + 109] = w;
            } else {
                throw PSIEXCEPTION("PseudoGrid: Allowed grid L values are 6, 8, 14, 26, 38, 50, 74, 86, 100");
            }            

            // Scale to radius
            for (int l = 0; l < L; l++) {
                xp[counter + l] *= r;
                yp[counter + l] *= r;
                zp[counter + l] *= r;
            }

            // Center
            for (int l = 0; l < L; l++) {
                xp[counter + l] += xc;
                yp[counter + l] += yc;
                zp[counter + l] += zc;
            }
           
            // TODO rotate to standard orientation  
 
            counter += L;
 
        }
    }

    fprintf(outfile, "  Pseudospectral Grid Points:\n");
    for (int P = 0; P < npoints; P++) {
        fprintf(outfile, "   P = %5d: (%24.16E, %24.16E, %24.16E) x %24.16E\n", P, xp[P], yp[P], zp[P], wp[P]);
    }

}


}

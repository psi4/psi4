#include "mints.h"

#include <psi4-dec.h>

#include <boost/regex.hpp>
#include <boost/xpressive/xpressive.hpp>
#include <boost/xpressive/regex_actions.hpp>
#include <boost/algorithm/string.hpp>

#include <cstdio>
#include <fstream>
#include <algorithm>
#include <ctype.h>

using namespace psi;
using namespace boost;
using namespace std;

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

BasisSetParser::BasisSetParser(const std::string& _searchpath)
   : searchpath_(_searchpath)
{
    // If the search path is empty use either PSIDATADIR or INSTALLEDPSIDATADIR
    if (searchpath_.empty()) {
        searchpath_ = Process::environment("PSIDATADIR") + "/basis";
    }
}

BasisSetParser::~BasisSetParser()
{
}

void Gaussian94BasisSetParser::parse(shared_ptr<BasisSet>& basisSet, const vector<string> &basisnames)
{
    // Obtain hold of molecule object from the basis set, we'll need to walk through this
    shared_ptr<Molecule> molecule = basisSet->molecule();

    // Ensure that the number of atoms match the number of basis names
    if (molecule->natom() != basisnames.size()) {
        fprintf(outfile, "Gaussian94BasisSetParser::parse: Number of atoms does not match number of basis set names.\n");
        fprintf(outfile, "molecule->natom() = %d     basisnames.size() = %lu\n", molecule->natom(), basisnames.size());
        fflush(outfile);
        throw PSIEXCEPTION("Gaussian94BasisSetParser::parse: Number of atoms does not match number of basis set names.");
    }

    // Regular expressions that we'll be checking for.
    regex cartesian("^cartesian", regbase::icase);
    regex spherical("^spherical", regbase::icase);
    regex comment("^\\s*\\!.*");                                       // line starts with !
    regex separator("^\\*\\*\\*\\*");                                  // line starts with ****
    regex atom_array("^\\s*([A-Za-z]+)\\s+0.*");                       // array of atomic symbols terminated by 0
    regex shell("^\\s*(\\w+)\\s*(\\d+)\\s*(-?\\d+\\.\\d+)");           // Match beginning of contraction

#define NUMBER "((?:-?\\d*\\.\\d+(?:[Ee][-+]\\d+)?)|(?:-?\\d+\\.\\d*(?:[Ee][-+]\\d+)?))"

    regex primitives1("^\\s*" NUMBER "\\s+" NUMBER ".*");    // Match s, p, d, f, g, ... functions
    regex primitives2("^\\s*" NUMBER "\\s+" NUMBER "\\s+" NUMBER ".*"); // match sp functions
    regex primitives3("^\\s*" NUMBER "\\s+" NUMBER "\\s+" NUMBER "\\s+" NUMBER ".*"); // match spd functions

    // s, p and s, p, d can be grouped together in Pople-style basis sets
    const string sp("SP"), spd("SPD");

    //                     a  b  c  d  e  f  g  h  i  j  k  l  m  n  o  p  q  r  s  t  u  v  w  x  y  z
    char shell_to_am[] = {-1,-1,-1, 2,-1, 3, 4, 5, 6,-1, 7, 8, 9,10,11, 1,12,13, 0,14,15,16,17,18,19,20};

    // Hold the result of a regex_match
    smatch what;

    for (int atom=0; atom<molecule->natom(); ++atom) {
        // Basis type.
        GaussianType gaussian_type = Pure;

        // Pull atomic center location
        Vector3 center = molecule->xyz(atom);

        // Modify the name of the basis set to generate a filename: STO-3G -> sto-3g
        string basisname = basisnames[atom];
        // First make it lower case
        transform(basisname.begin(), basisname.end(), basisname.begin(), ::tolower);

        string format_underscore("_"); // empty string
        // Replace all '(' with '_'
        xpressive::sregex match_format = xpressive::as_xpr("(");
        basisname = regex_replace(basisname, match_format, format_underscore);

        // Replace all ')' with '_'
        match_format = xpressive::as_xpr(")");
        basisname = regex_replace(basisname, match_format, format_underscore);

        // Replace all ',' with '_'
        match_format = xpressive::as_xpr(",");
        basisname = regex_replace(basisname, match_format, format_underscore);

        // Replace all '*' with 's'
        match_format = xpressive::as_xpr("*");
        string format_star("s");
        basisname = regex_replace(basisname, match_format, format_star);

        // Replace all '+' with 'p'
        match_format = xpressive::as_xpr("+");
        string format_plus("p");
        basisname = regex_replace(basisname, match_format, format_plus);

//        cout << " basisname with '-' removed: " << basisname << endl;
        string basis_filename = searchpath() + "/" + basisname + ".gbs";
//        cout << " will attempt to read from: " << basis_filename << endl;

        // Load in entire file.
        vector<string> lines;
        string text;
        ifstream infile(basis_filename.c_str());
        if (!infile)
            throw PSIEXCEPTION("Gaussian94BasisSetParser::parse: Unable to open basis set file: " + basis_filename);
        while (infile.good()) {
            getline(infile, text);
            lines.push_back(text);
        }

        int lineno = 0;
        bool found = false;
        while (lineno < lines.size()) {
            string line = lines[lineno++];

            // Spit out the line for debugging.
//            cout << line << endl;

            // Ignore blank lines
            if (line.empty())
                continue;

            // Look for Cartesian or Spherical
            if (regex_match(line, what, cartesian)) {
                gaussian_type = Cartesian;
                if (Process::environment.options.get_global("PUREAM").has_changed()) {
                    gaussian_type = ((Process::environment.options.get_global("PUREAM").to_integer()) ? Pure : Cartesian); 
                }
                continue;
            } else if (regex_match(line, what, spherical)) {
                gaussian_type = Pure;
                if (Process::environment.options.get_global("PUREAM").has_changed()) {
                    gaussian_type = ((Process::environment.options.get_global("PUREAM").to_integer()) ? Pure : Cartesian); 
                }
                continue;
            }

            // Do some matches
            if (regex_match(line, what, comment)) {
//                cout << " matched a comment line.\n";
                continue;
            }
            if (regex_match(line, what, separator)) {
//                cout << " matched a separator line.\n";
                continue;
            }
            // Match: H    0
            // or:    H    O...     0
            if (regex_match(line, what, atom_array)) {
//                cout << " matched a atom array line.\n";
//                cout << " what.captures(1)" << what[1].str() << "\n";

                // Check the captures and see if this basis set is for the atom we need.
                found = false;
                if (iequals(molecule->label(atom), what[1].str())) {
                    found = true;

//                    cout << "found\n";

                    // Read in the next line
                    line = lines[lineno++];

                    // Need to do the following until we match a "****" which is the end of the basis set
                    while (!regex_match(line, what, separator)) {
                        // Match shell information
                        if (regex_match(line, what, shell)) {
                            string shell_type(what[1].first, what[1].second);
                            std::transform(shell_type.begin(), shell_type.end(), shell_type.begin(), ::toupper);
                            int nprimitive;
                            double scale;
                            double exponent, contraction;

                            if (!from_string<int>(nprimitive, what[2], std::dec))
                                throw PSIEXCEPTION("Gaussian94BasisSetParser::parse: Unable to convert number of primitives:\n" + line);
                            if (!from_string<double>(scale, what[3], std::dec))
                                throw PSIEXCEPTION("Gaussian94BasisSetParser::parse: Unable to convert scale factor:\n" + line);

                            if (shell_type.size() == 1) {
                                int am = (int)shell_to_am[shell_type[0] - 'A'];

                                double *exponents = new double[nprimitive];
                                double *contractions = new double[nprimitive];

                                for (int p=0; p<nprimitive; ++p) {
                                    line = lines[lineno++];

                                    // Must match primitives1; will work on the others later
                                    if (!regex_match(line, what, primitives1))
                                        throw PSIEXCEPTION("Gaussian94BasisSetParser::parse: Unable to match an exponent with one contraction:\n" + line);

                                    if (!from_string<double>(exponent, what[1], std::dec))
                                        throw PSIEXCEPTION("Gaussian94BasisSetParser::parse: Unable to convert exponent:\n" + line);
                                    if (!from_string<double>(contraction, what[2], std::dec))
                                        throw PSIEXCEPTION("Gaussian94BasisSetParser::parse: Unable to convert contraction:\n" + line);

                                    // Scale the contraction
                                    contraction *= scale;

                                    // Save the information.
                                    exponents[p] = exponent;
                                    contractions[p] = contraction;
                                }

//                                printf("Adding new shell. nprimitive = %d\n", nprimitive);
                                // We have a full shell, push it to the basis set
                                shared_ptr<GaussianShell> new_shell(new GaussianShell);
                                new_shell->init(nprimitive,
                                                exponents,
                                                am,
                                                gaussian_type,
                                                contractions,
                                                atom,
                                                center,
                                                0);
                                new_shell->normalize_shell();
                                basisSet->shells_.push_back(new_shell);

                                delete[] exponents;
                                delete[] contractions;
                            }
                            else if (shell_type.size() == 2) {
                                // This is to handle instances of SP, PD, DF, FG, ...
                                int am1 = (int)shell_to_am[shell_type[0] - 'A'];
                                int am2 = (int)shell_to_am[shell_type[1] - 'A'];

                                double *exponents = new double[nprimitive];
                                double *contractions1 = new double[nprimitive];
                                double *contractions2 = new double[nprimitive];

                                for (int p=0; p<nprimitive; ++p) {
                                    line = lines[lineno++];

                                    // Must match primitivies2;
                                    if (!regex_match(line, what, primitives2))
                                        throw PSIEXCEPTION("Gaussian94BasisSetParser::parse: Unable to match an exponent with two contractions:\n" + line);

                                    if (!from_string<double>(exponent, what[1], std::dec))
                                        throw PSIEXCEPTION("Gaussian94BasisSetParser::parse: Unable to convert exponent:\n" + line);
                                    if (!from_string<double>(contraction, what[2], std::dec))
                                        throw PSIEXCEPTION("Gaussian94BasisSetParser::parse: Unable to convert first contraction:\n" + line);

                                    // Scale the contraction
                                    contraction *= scale;

                                    // Save the information
                                    exponents[p] = exponent;
                                    contractions1[p] = contraction;

                                    // Do the other contraction
                                    if (!from_string<double>(contraction, what[3], std::dec))
                                        throw PSIEXCEPTION("Gaussian94BasisSetParser::parse: Unable to convert second contraction:\n" + line);

                                    // Scale the contraction
                                    contraction *= scale;

                                    // Save the information
                                    contractions2[p] = contraction;
                                }

//                                printf("Adding 2 new shells. nprimitive = %d\n", nprimitive);
                                shared_ptr<GaussianShell> new_shell(new GaussianShell);
                                new_shell->init(nprimitive,
                                                exponents,
                                                am1,
                                                gaussian_type,
                                                contractions1,
                                                atom,
                                                center,
                                                0);
                                new_shell->normalize_shell();
                                basisSet->shells_.push_back(new_shell);

                                new_shell = shared_ptr<GaussianShell>(new GaussianShell);
                                new_shell->init(nprimitive,
                                                exponents,
                                                am2,
                                                gaussian_type,
                                                contractions2,
                                                atom,
                                                center,
                                                0);
                                new_shell->normalize_shell();
                                basisSet->shells_.push_back(new_shell);

                                delete[] exponents;
                                delete[] contractions1;
                                delete[] contractions2;
                            }
                            else {
                                throw PSIEXCEPTION("Gaussian94BasisSetParser::parse: Unable to parse basis sets with spd, or higher grouping\n");
                            }
                        } else {
                            throw PSIEXCEPTION("Gaussian94BasisSetParser::parse: Expected shell information, but got:\n" + line);
                        }
                        line = lines[lineno++];
                    }
                    break;
                }
            }
        }
        if (found == false)
            throw PSIEXCEPTION("Gaussian94BasisSetParser::parser: Unable to find the basis set for " + molecule->label(atom));
    }

    // Have the basis set object refresh its internal data.
    basisSet->refresh();
}


#include "mints.h"

#include <psi4-dec.h>

#include <boost/shared_ptr.hpp>
#include <boost/regex.hpp>
#include <boost/xpressive/xpressive.hpp>
#include <boost/xpressive/regex_actions.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>

#include <cstdio>
#include <fstream>
#include <algorithm>
#include <ctype.h>

using namespace psi;
using namespace boost;
using namespace std;

boost::regex basis_separator("^\\s*\\[\\s*(.*?)\\s*\\]\\s*$");

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

BasisSetFileNotFound::BasisSetFileNotFound(string message,
                                   const char* _file,
                                   int _line) throw()
    : PsiException(message, _file, _line)
{
    stringstream sstr;
    sstr << "sanity check failed! " << message;
    rewrite_msg(sstr.str());
}

BasisSetFileNotFound::~BasisSetFileNotFound() throw()
{
}

BasisSetNotFound::BasisSetNotFound(string message,
                                   const char* _file,
                                   int _line) throw()
    : PsiException(message, _file, _line)
{
    stringstream sstr;
    sstr << "sanity check failed! " << message;
    rewrite_msg(sstr.str());
}

BasisSetNotFound::~BasisSetNotFound() throw()
{
}

BasisSetParser::BasisSetParser()
{
}

BasisSetParser::~BasisSetParser()
{
}

vector<string> BasisSetParser::load_file(const std::string& filename,
                                         const std::string& basisname)
{
    int me = Communicator::world->me();

    // Loads an entire file.
    vector<string> lines;

    if (Communicator::world->me() == 0) {
        smatch what;

        // temp variable
        string text;

        // Stream to use
        ifstream infile(filename.c_str());

        if (!infile)
            throw BasisSetFileNotFound("BasisSetParser::parse: Unable to open basis set file: " + filename, __FILE__, __LINE__);

        bool given_basisname = basisname.empty() ? false : true;
        bool found_basisname = false;

        while (infile.good()) {
            getline(infile, text);

            // If no basisname was given always save the line.
            if (given_basisname == false)
                lines.push_back(text);

            if (found_basisname) {

                // If we find another [*] we're done.
                if (regex_match(text, what, basis_separator))
                    break;

                lines.push_back(text);
                continue;
            }

            // If the user gave a basisname AND text matches the basisname we want to trigger to retain
            if (given_basisname && regex_match(text, what, basis_separator)) {
                if (boost::iequals(what[1].str(), basisname))
                    found_basisname = true;
            }
        }
    }

    size_t nlines = lines.size();
    if (Communicator::world->nproc() > 1) {
        Communicator::world->bcast(&nlines, 1);

        if (me > 0)
            lines.resize(nlines);

        for (size_t i=0; i<nlines; ++i) {
            Communicator::world->bcast(lines[i]);
        }
    }

    return lines;
}

vector<string> BasisSetParser::string_to_vector(const std::string &data)
{
    vector<string> lines;
    boost::split(lines, data, boost::is_any_of("\n"));
    return lines;
}

std::vector<boost::shared_ptr<GaussianShell> >
Gaussian94BasisSetParser::parse(const string& symbol, const std::vector<std::string> &lines)
{
    // Regular expressions that we'll be checking for.
    regex cartesian("^\\s*cartesian\\s*", regbase::icase);
    regex spherical("^\\s*spherical\\s*", regbase::icase);
    regex comment("^\\s*\\!.*");                                       // line starts with !
    regex separator("^\\s*\\*\\*\\*\\*");                                  // line starts with ****
    regex atom_array("^\\s*([A-Za-z]+)\\s+0.*");                       // array of atomic symbols terminated by 0
    regex shell("^\\s*(\\w+)\\s*(\\d+)\\s*(-?\\d+\\.\\d+)");           // Match beginning of contraction

    // NUMBER is in psi4-dec.h
    regex primitives1("^\\s*" NUMBER "\\s+" NUMBER ".*");    // Match s, p, d, f, g, ... functions
    regex primitives2("^\\s*" NUMBER "\\s+" NUMBER "\\s+" NUMBER ".*"); // match sp functions
    regex primitives3("^\\s*" NUMBER "\\s+" NUMBER "\\s+" NUMBER "\\s+" NUMBER ".*"); // match spd functions

    // s, p and s, p, d can be grouped together in Pople-style basis sets
    const string sp("SP"), spd("SPD");

    //                     a  b  c  d  e  f  g  h  i  j  k  l  m  n  o  p  q  r  s  t  u  v  w  x  y  z
    char shell_to_am[] = {-1,-1,-1, 2,-1, 3, 4, 5, 6,-1, 7, 8, 9,10,11, 1,12,13, 0,14,15,16,17,18,19,20};

    // Hold the result of a regex_match
    smatch what;

    // Basis type.
    GaussianType gaussian_type = Pure;

    // Need a dummy center for the shell.
    Vector3 center;

    vector<boost::shared_ptr<GaussianShell> > shell_list;

    int lineno = 0;
    bool found = false;

    while (lineno < lines.size()) {
        string line = lines[lineno++];

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
            continue;
        }
        if (regex_match(line, what, separator)) {
            continue;
        }

        // Match: H    0
        // or:    H    O...     0
        if (regex_match(line, what, atom_array)) {
            // Check the captures and see if this basis set is for the atom we need.
            found = false;
            if (iequals(symbol, what[1].str())) {
                found = true;

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

                                int idx;
                                while((idx=line.find_first_of('D')) >= 0 ) {
                                    line.replace( idx, 1, "e" );
                                }
                                while((idx=line.find_first_of('d')) >= 0 ) {
                                    line.replace( idx, 1, "e" );
                                }

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
                            boost::shared_ptr<GaussianShell> new_shell(new GaussianShell);
                            new_shell->init(nprimitive,
                                            exponents,
                                            am,
                                            gaussian_type,
                                            contractions,
                                            0,
                                            center,
                                            0);
                            new_shell->normalize_shell();
                            shell_list.push_back(new_shell);

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

                                int idx;
                                while((idx=line.find_first_of('D')) >= 0 ) {
                                    line.replace( idx, 1, "e" );
                                }
                                while((idx=line.find_first_of('d')) >= 0 ) {
                                    line.replace( idx, 1, "e" );
                                }

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
                            boost::shared_ptr<GaussianShell> new_shell(new GaussianShell);
                            new_shell->init(nprimitive,
                                            exponents,
                                            am1,
                                            gaussian_type,
                                            contractions1,
                                            0,
                                            center,
                                            0);
                            new_shell->normalize_shell();
                            shell_list.push_back(new_shell);

                            new_shell = boost::shared_ptr<GaussianShell>(new GaussianShell);
                            new_shell->init(nprimitive,
                                            exponents,
                                            am2,
                                            gaussian_type,
                                            contractions2,
                                            0,
                                            center,
                                            0);
                            new_shell->normalize_shell();
                            shell_list.push_back(new_shell);

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
        throw BasisSetNotFound("Gaussian94BasisSetParser::parser: Unable to find the basis set for " + symbol, __FILE__, __LINE__);

    // The constructor, or the caller, should refresh the basis set.
    return shell_list;
}


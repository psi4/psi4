/*!
 \file
 \ingroup MINTS
 */

#include <boost/regex.hpp>
#include <cstdlib>
#include <unistd.h>
#include <cstring>
#include <boost/shared_ptr.hpp>
#include <psi4-dec.h>
#include <libmints/molecule.h>

using namespace boost;

namespace psi {


/*!
** get_writer_file_prefix()
**
** For archival files such as the Hessian, Molden files, etc., we want
** to write a file in the current working directory with a name like
** benzene.molden, benzene.hess, etc.  This function returns the basename
** of such files.  The basename is determined as follows: If the option
** "WRITER_FILE_LABEL" exists, we will use that.  It could be user-specified, 
** or created by a python driver doing a complex task (e.g., 
** S22-1-dimer-mp2-ccpvdz).  If the label option is blank (default), 
** then we will just make up a prefix using the basename of the output
** file along with the active molecule name (e.g., test.h2o).
**
** C. David Sherrill and Lori Burns
**
** Returns: C++ standard string of the filename prefix
**
** \ingroup MINTS
*/

std::string get_writer_file_prefix(void)
{

  std::string label = Process::environment.options.get_str("WRITER_FILE_LABEL");
  if (label != "") {
    return(label);
  }
 
  // If no available options WRITER_FILE_LABEL, then we build a defult:
  // Get the basename of the output filename, append any active molecule name
  // to it, and return the resulting string

  boost::regex outfileBase("(\\w+)(\\.out|\\.dat)", boost::regbase::normal | boost::regbase::icase);
  boost::smatch reMatches;

  // outfile_name is in psi4-dec.h and is a global std::string with the
  // name of the output file
  std::string prefix;
  if (regex_match(outfile_name, reMatches, outfileBase)) {
    prefix = reMatches[1].str();
  }
  else {
    prefix = outfile_name;
  }

  std::string molecule_name = Process::environment.molecule()->name();
  if (molecule_name != "") {
    prefix += "." + molecule_name;
  }

  return(prefix);
}


}


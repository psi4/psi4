#include <cstdio>

namespace psi {

/*! Print PSI version information that was set in configure.ac */
void print_version(FILE *myout)
{
  const char *PSI_VERSION = "0.1";
  fprintf(myout, "    -----------------------------------------------------------------------\n");
  fprintf(myout, "          PSI4: An Open-Source Ab Initio Electronic Structure Package\n");
  fprintf(myout, "                              PSI %s Driver\n", PSI_VERSION);
  fprintf(myout, "    T. D. Crawford, C. D. Sherrill, E. F. Valeev, J. T. Fermann, R. A. King,\n");
  fprintf(myout, "    M. L. Leininger, S. T. Brown, C. L. Janssen, E. T. Seidl, J. P. Kenny,\n");
  fprintf(myout, "    and W. D. Allen, J. Comput. Chem. 28, 1610-1616 (2007)\n");
  fprintf(myout, "\n");
  fprintf(myout, "                         Additional Contributions by\n");
  fprintf(myout, "    Francesco Evangelista, Andrew Simmonett, Justin Turney, Jeremy Wilke\n");
  fprintf(myout, "    -----------------------------------------------------------------------\n\n");
}

}

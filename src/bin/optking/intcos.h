/*! \file
    \ingroup OPTKING
    \brief Class for stretches
*/

#ifndef _psi3_src_bin_optking_intcos_h_
#define _psi3_src_bin_optking_intcos_h_

#include <iostream>
#include <vector>
#include <string>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libipv1/ip_lib.h>
#include <cov_radii.h>
#include <physconst.h>
#include <functional>
#include <algorithm>

extern "C" {
  extern FILE *infile, *outfile;
  extern char *psi_file_prefix;
}

#include "v3d.h"
using namespace psi::v3d;
#include "stretch.h"
#include "bend.h"

namespace psi { namespace optking {

using std::vector;
using std::string;
using std::bind1st;

class Intcos {

  private:
    int natom;      // # of atoms in geom
    double **geom;  //cartesian coordinates in Angstroms
    double *zvals;  // atomic numbers
    double *masses;

    int nallatom;   // # of atoms in geom
    double **fgeom; //cartesian coordinates in Angstroms

    vector<Stretch> stre;
    vector<Bend> bend;

    int add_stretches_from_input(string keyword = "STRE");
    int add_bends_from_input(string keyword = "BEND");

#include "tmpl.h" // intcos templates

  public:

    double **Bsimp;     //B matrix

    Intcos();
    ~Intcos();
    Intcos(const Intcos & s);

    int size(void) const;
    void print(int print_flag) const;
    void print_s(void) const;
    void print_Bsimp(void) const;
    int Bsimp_present;

    void compute(void);
    void compute_s(void);
    int add_intcos_from_input(void);
    int add_stretches_by_distance(double scale_connectivity = 1.3);
    int add_angles_by_bonds(void);
    int build_Bsimp(void);
    double **get_Bsimp(void) {return Bsimp;}
//    friend class Geom_opt;
};

}}

#endif

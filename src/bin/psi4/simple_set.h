/*! \file simple.h
    \ingroup OPT09
    \brief simple internal coordinate base class
*/

#ifndef _psi4_src_bin_opt09_simple_set_h_
#define _psi4_src_bin_opt09_simple_set_h_

#include <psi4-dec.h>

#include "simple.h"
#include "stretch.h"
#include "bend.h"
#include "stretch.h"

#include <vector>
//#include <typeinfo>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>

using std::vector;

namespace psi {

class SIMPLE_SET {

  vector<SIMPLE *> simples;

  public: 

    // default constructor

    ~SIMPLE_SET() {
      for(int i=0; i<simples.size(); ++i)
        delete simples.at(i);
    }

    int size(void) const { return simples.size(); };

    //compute values of each simple coordinate - input geometry is in au
    void compute(double **geom) {
      for(int i=0; i<simples.size(); ++i)
        simples.at(i)->compute(geom);
    }

    //compute s vectors of each simple coordinate - input geometry is in au
    void compute_s(double **geom) {
      for(int i=0; i<simples.size(); ++i)
        simples.at(i)->compute_s(geom);
    }

    //print values of each simple coordinate
    void print(int print_flag) const {
      for(int i=0; i<simples.size(); ++i)
        simples.at(i)->print(print_flag);
    }

    //print s vectors to output file
    void print_s(void) const {
      for(int i=0; i<simples.size(); ++i)
        simples.at(i)->print_s();
    }

    // read simples from intco: section
    int add_simples_from_input(void);

    // automatically generate simples from bond distances
    int add_simples_by_distance(int natom, double *Z, double *geom, double scale_connectivity);

    // tells whether a given simple is already present in object
    bool SIMPLE_SET::present(const SIMPLE * s2) const;

    // tells whether a given id number is alreaady in use
    bool SIMPLE_SET::id_present(int a) const;

    //returns index of simple that matches the given one, if not present returns -1
    int find_index(const SIMPLE * s2) const;

    // returns if two simple coordinates are identical
    friend bool equiv(const SIMPLE * s1, const SIMPLE * s2) ;

    // returns array of internal coordinate values
    double * get_q(void) const ;

    // returns B matrix
    double ** get_B(int natom) const;

};

}

#endif

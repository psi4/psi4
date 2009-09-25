#ifndef _psi_bin_extrema_simple_h_
#define _psi_bin_extrema_simple_h_

namespace psi { namespace extrema {

/*##########################################################################*/
/*! 
  \file
  \ingroup EXTREMA
  \brief Simple coordinate concrete class declaration and definition.
  
  definition for the internal type:
  <ul>   
  <li> bond length = 0
  <li> valence angle = 1
  <li> torsional angle = 2
  </ul> */

/*! \class simple
  \brief Class to hold simple definition and interface with other classes. */

/*						Joseph P. Kenny 11/29/01
  ###########################################################################*/

extern void punt(const char*);
extern void stop_io();

class simple {
    int type, equiv_grp;
    double val;
    int atom, bond, angle, tors;
    int opt;
    char *label;
  public:
     void set_simple(int ty, double value, int at, 
		     int bd, int an, int tr, int op) {
	 type=ty; val=value; atom=at; bond=bd; angle=an; tors=tr; opt=op;
	 return;
     }    
    void set_equiv_grp(int grp_num) { equiv_grp = grp_num; return; }
    int get_equiv_grp() {return equiv_grp;}
    void set_val(double value) { val = value; return; }
    double get_val() { return val; }
    void set_label(char *lab) { label = lab; return; }
    char *get_label() { return label; }
    int get_type() { return type; }
    int get_atom() { return atom; }
    int get_bond() { return bond; }
    int get_opt() { return opt; }
    int get_angle() { 
	if(type==0)
	    punt("class simple error: non-existent angle");
	return angle;
    }
    int get_tors() {
	if(type<2)
	    punt("class simple error: non-existent torsion");
	return tors;
    }  
};

/*! \fn simple::set_simple(int ty, double value, int at, 
		   int bd, int an, int tr, int op) 
  \brief Sets definition for a simple internal coordinate.
  \param ty type of coordinate
  \param value value of coordinate
  \param at reference atom 1
  \param bd atom (2) bonded to 1 
  \param an atom (3) defining angle 1-2-3
  \param tr atom (4) defining torsion 1-2-3-4
  \param op 1 if optimized, 0 otherwise */

/*! \fn simple::set_equiv_grp(int grp_num)
  \brief Sets the equivalent group of a simple.
  
  Simple internals are given a positive integer, coordinates with 
  matching integers were input as symmetrically equivalent by the user.
  \param grp_num group number */

/*! \fn simple::get_equiv_grp()
  \brief Returns number of equivalent group. */

/*! \fn simple::set_val(double value)
  \brief Sets value of simple. */

/*! \fn simple::get_val()
  \brief Returns value of simple. */

/*! \fn simple::set_label(char *lab)
  \brief Sets label of simple. */

/*! \fn simple::get_label()
  \brief Returns label of simple. */

/*! \fn simple::get_type()
  \brief Returns simple type. */

/*! \fn simple::get_atom()
  \brief Returns reference atom 1. */

/*! \fn simple::get_bond()
  \brief Returns atom bonded to reference atom 1. */

/*! \fn simple::get_opt()
  \brief Returns 1 if variable flagged to be optimized by user, 0 otherwise. */

/*! \fn simple::get_angle()
  \brief Returns atom (3) which defines angle 1-2-3. */

/*! \fn simple::get_tors()
  \brief Returns atom (4) which defines torsion 1-2-3-4 */

}} // namespace psi::extrema

#endif // header guard

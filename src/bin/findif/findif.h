#ifndef _findif_h_
#define _findif_h_

#include <sstream>
#include <vector>
#include <libmints/mints.h>
#include <psi4-dec.h>
#include <libciomr/libciomr.h>

namespace psi { namespace findif {

// functions to generate displacements
std::vector< SharedMatrix > fd_geoms_1_0(Options &options);
std::vector< SharedMatrix > fd_geoms_2_0(Options &options);
std::vector< SharedMatrix > fd_geoms_freq_0(Options &options, int irrep=-1);
std::vector< SharedMatrix > fd_geoms_freq_1(Options &options, int irrep=-1);

// functions to carry out finite-differences
PsiReturnType fd_1_0(Options &options, const boost::python::list& energies);
PsiReturnType fd_2_0(Options &options, const boost::python::list& energies);
PsiReturnType fd_freq_0(Options &options, const boost::python::list& energies, int irrep=-1);
PsiReturnType fd_freq_1(Options &options, const boost::python::list& E_list, int irrep=-1);

// class to accumulate and print vibrations
class VIBRATION {
  int irrep;       // irrep
  double km;    // force constant
  double *lx;   // normal mode in mass-weighted cartesians
  double cm;    // harmonic frequency in wavenumbers
  
  public:
    friend PsiReturnType fd_freq_0(Options &options, const boost::python::list& energies, int irrep);
    friend PsiReturnType fd_freq_1(Options &options, const boost::python::list& gradients, int irrep);
    friend bool ascending(const VIBRATION *, const VIBRATION *);
    friend void print_vibrations(std::vector<VIBRATION *> modes);
  
    VIBRATION(int irrep_in, double km_in, double *lx_in) { irrep = irrep_in; km = km_in; lx = lx_in; }
    ~VIBRATION() { free(lx); }
};

// function to print vibrations
void print_vibrations(std::vector<VIBRATION *> modes);
  
// to order vibrations
bool ascending(const VIBRATION *vib1, const VIBRATION *vib2);

// for displacing along a salc
void displace_cart(SharedMatrix geom, const CdSalcList & salclist,
  int salc_i, int disp_factor, double disp_size);

void displace_cart(SharedMatrix geom, const CdSalcList & salclist,
  int salc_i, int salc_j, int disp_factor_i, int disp_factor_j, double disp_size);

// to massweight columns of a shared matrix
void mass_weight_columns_plus_one_half(SharedMatrix B);


template <class T>
inline std::string to_string (const T& t)
{
std::stringstream ss;
ss << t;
return ss.str();
}

}}

#endif


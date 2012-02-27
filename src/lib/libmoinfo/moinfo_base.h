#ifndef _psi_src_lib_libmoinfo_moinfo_base_h_
#define _psi_src_lib_libmoinfo_moinfo_base_h_

/*! \file    moinfo_base.h
    \ingroup LIBMOINFO
    \brief   This class stores all the basic info regarding MOs
*/

#define PSI_NULL(args) args = NULL;
#define PSI_FREE(args) if(args != NULL) free(args);
#define PSI_DELETE(args) if(args != NULL) delete args;
#define PSI_DELETE_ARRAY(args) if(args != NULL) delete[] args;
#define IOFF 5000000

#include <string>
#include <libutil/libutil.h>

typedef std::vector<int>                    intvec;
typedef std::vector<bool>                   boolvec;

namespace psi {

class Options;

class MOInfoBase{
public:
  MOInfoBase(Options& options_, bool silent_ = false);
  ~MOInfoBase();

  double      get_nuclear_energy()               const {return(nuclear_energy);}

  char**      get_irr_labs()                     const {return(irr_labs);}
  char*       get_irr_labs(int i)                const {return(irr_labs[i]);}

  int         get_nirreps()                      const {return(nirreps);}
  int         get_nso()                          const {return(nso);}

  size_t*     get_ioff()                         const {return(ioff);}
  intvec      get_sopi()                         const {return(sopi);}
  intvec      get_docc()                         const {return(docc);}
  intvec      get_actv()                         const {return(actv);}
  bool        get_guess_occupation()             const {return(guess_occupation);}
  int         get_ndocc()                        const {return(ndocc);}
  int         get_nactv()                        const {return(nactv);}

  int         get_nael()                         const {return(nael);}                // # of alpha electrons including frozen
  int         get_nbel()                         const {return(nbel);}                // # of  beta electrons including frozen

  double**    get_scf_mos()                      const {return(scf);}
  double**    get_scf_mos(int i)                 const {return(scf_irrep[i]);}
  double      get_scf_mos(int i,int j)           const {if((i<nmo)&&(j<nso)) return(scf[i][j]); else return(0.0);}
  void        write_chkpt_mos();
protected:
  void        read_data();
  void        compute_number_of_electrons();
  void        correlate(char *ptgrp, int irrep, int& nirreps_old, int& nirreps_new,int*& correlation);
  void        read_mo_space(int nirreps_ref, int& n, intvec& mo, std::string labels);
  void        print_mo_space(int& nmo, intvec& mo, std::string labels);
  intvec      convert_int_array_to_vector(int n, const int* array);

  void        startup();
  void        cleanup();
  void        compute_ioff();

  Options&    options;
  int         nirreps;
  int         wfn_sym;
  int         charge;
  int         multiplicity;

  int         nso;              // PSI nso (number of symmetry-adapted atomic orbitals)
  int         nmo;              // Psi nmo (number of molecular orbitals, including frozen core and frozen virtual)
  int         ndocc;
  int         nactv;
  int         nael;
  int         nbel;
  int         nactive_ael;
  int         nactive_bel;

  size_t*     ioff;
  intvec      sopi;
  intvec      docc;
  intvec      actv;
  bool        guess_occupation;
  bool        silent;

  double      nuclear_energy;

  double**    scf;                                   // MO coefficients
  double***   scf_irrep;                             // MO coefficients

  char**      irr_labs;
};

}

#endif // _psi_src_lib_libmoinfo_moinfo_base_h_

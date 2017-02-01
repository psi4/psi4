/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#ifndef _psi_src_lib_libmoinfo_moinfo_h_
#define _psi_src_lib_libmoinfo_moinfo_h_

#include <bitset>
#include <string>
#include <vector>
#include <utility>

#include "moinfo_base.h"

#define size_det 2048

namespace psi {

enum    ReferenceType {AllRefs,UniqueRefs,ClosedShellRefs,UniqueOpenShellRefs};

class MOInfo : public MOInfoBase
{

  typedef std::vector<std::string>            strvec;
  typedef std::vector<std::pair<int,int> >    intpairvec;
public:
  /*********************************************************
    SlaterDeterminant Class
    1) Purpose
      This class is used to store all the information that
      belongs to a Slater Determinant
    2) Use
    3) Details
      The MOs that describe the reference are stored in the
      arrays aocc,bocc,avir and bvir. These refer to the MOs
      in Pitzer order with the frozen occupied and virtual
      already excluded. Therefore this code assumes that the
      integral transformation code has already eliminated
      the frozen integrals and relabeled them.

      type stores the type of reference
      type = 0 -> closed-shell determinant
      type = 2 -> open-shell   determinant

      number stores the ID of this reference

    4) Uses
      STL vector
  *********************************************************/
  class SlaterDeterminant{
      const MOInfo* moinfo;
  public:
    typedef std::bitset<size_det> bitdet;
    SlaterDeterminant(const MOInfo *);
    ~SlaterDeterminant();
    void        set(int n)                               {bits.set(n);}
    bool        test(int n)                        const {return(bits.test(n));}
    bool        is_closed_shell();
    bool        is_spin_flipped(SlaterDeterminant& det);
    bitdet&     get_bits()                               {return(bits);}
    void        get_internal_excitations(SlaterDeterminant& det, double& sign,
                                                  std::vector<std::pair<int,int> >& alpha_operators,
                                                  std::vector<std::pair<int,int> >& beta_operators);
    char        get_occupation_symbol(int i);
    std::string get_label();
    intvec      get_aocc();
    intvec      get_bocc();
    intvec      get_avir();
    intvec      get_bvir();
    boolvec     get_is_aocc();
    boolvec     get_is_bocc();
    boolvec     get_is_avir();
    boolvec     get_is_bvir();

  private:
    double      annihilate(bitdet& bits_det,int so);
    double      create(bitdet& bits_det,int so);
    bitdet      bits;
    std::string type;
  };
public:
  friend class SlaterDeterminant;
  MOInfo(Wavefunction& ref_wf_, Options& options_, bool silent_ = false);
  ~MOInfo();

  // DGEMM timing
  void        set_dgemm_timing(double value)           {dgemm_timing=value;}
  void        add_dgemm_timing(double value)           {dgemm_timing+=value;}
  double      get_dgemm_timing()                 const {return(dgemm_timing);}

  // Convergence Options
  double      get_no_damp_convergence()          const {return(no_damp_convergence);}

  intvec      get_mo_sym()                       const {return(all_sym);}
  int         get_mo_sym(int i)                  const {return(all_sym[i]);}
  int         get_wfn_sym()                      const {return(wfn_sym);}
  int         get_root()                         const {return(root);}

  int         get_nmo()                          const {return(nmo);}
  int         get_nactive_ael()                  const {return(nactive_ael);}
  int         get_nactive_bel()                  const {return(nactive_bel);}


  int         get_nall()                         const {return(nall);}
  int         get_nfocc()                        const {return(nfocc);}
  int         get_nextr()                        const {return(nextr);}
  int         get_nfvir()                        const {return(nfvir);}
  int         get_nocc()                         const {return(nocc);}
  int         get_nvir()                         const {return(nvir);}

  intvec      get_sopi()                         const {return(sopi);}
  intvec      get_mopi()                         const {return(mopi);}
  intvec      get_docc()                         const {return(docc);}
  intvec      get_actv()                         const {return(actv);}
  intvec      get_focc()                         const {return(focc);}
  intvec      get_extr()                         const {return(extr);}
  intvec      get_fvir()                         const {return(fvir);}
  intvec      get_occ()                          const {return(occ);}
  intvec      get_vir()                          const {return(vir);}
  intvec      get_all()                          const {return(all);}

  int         get_sopi(int i)                    const {return(sopi[i]);}
  int         get_mopi(int i)                    const {return(mopi[i]);}
  int         get_focc(int i)                    const {return(focc[i]);}
  int         get_docc(int i)                    const {return(docc[i]);}
  int         get_actv(int i)                    const {return(actv[i]);}
  int         get_extr(int h)                    const {return(extr[h]);}
  int         get_fvir(int i)                    const {return(fvir[i]);}

  // Mapping functions
  intvec      get_focc_to_mo()                   const {return(focc_to_mo);}
  intvec      get_docc_to_mo()                   const {return(docc_to_mo);}
  intvec      get_actv_to_mo()                   const {return(actv_to_mo);}
  intvec      get_extr_to_mo()                   const {return(extr_to_mo);}
  intvec      get_fvir_to_mo()                   const {return(fvir_to_mo);}
  intvec      get_occ_to_mo()                    const {return(occ_to_mo);}
  intvec      get_vir_to_mo()                    const {return(vir_to_mo);}
  intvec      get_all_to_mo()                    const {return(all_to_mo);}
  intvec      get_mo_to_all()                    const {return(mo_to_all);}
  intvec      get_actv_to_occ()                  const {return(actv_to_occ);}
  intvec      get_actv_to_vir()                  const {return(actv_to_vir);}
  intvec      get_occ_to_actv()                  const {return(occ_to_actv);}
  intvec      get_vir_to_actv()                  const {return(vir_to_actv);}
  boolvec     get_is_actv_in_occ()               const {return(is_actv_in_occ);}
  boolvec     get_is_actv_in_vir()               const {return(is_actv_in_vir);}

  int         get_all_to_occ(int i)              const {return(all_to_occ[i]);}
  int         get_all_to_vir(int i)              const {return(all_to_vir[i]);}
  int         get_all_to_mo(int i)               const {return(all_to_mo[i]);}



  double      get_scf_energy()                   const {return(scf_energy);}
  double      get_fzcore_energy()                const {return(fzcore_energy);}
  void        set_fzcore_energy(double efzc)           {fzcore_energy=efzc;}

  // Model space functions
  void        setup_model_space();
  int         get_nrefs()                              {return(all_refs.size());};
  int         get_nunique()                            {return(unique_refs.size());};
  int         get_ref_number(int n, ReferenceType ref_type = AllRefs);
  int         get_ref_size(ReferenceType ref_type);
  std::string get_determinant_label(int i);

  strvec      get_matrix_names(std::string str);
  intvec      get_aocc(int i, ReferenceType ref_type);
  intvec      get_bocc(int i, ReferenceType ref_type);
  intvec      get_avir(int i, ReferenceType ref_type);
  intvec      get_bvir(int i, ReferenceType ref_type);

  boolvec     get_is_aocc(int i, ReferenceType ref_type);
  boolvec     get_is_bocc(int i, ReferenceType ref_type);
  boolvec     get_is_avir(int i, ReferenceType ref_type);
  boolvec     get_is_bvir(int i, ReferenceType ref_type);

  intvec      get_determinant(int i);  // Array with occupation of reference i (in all ordering)

  intpairvec  get_alpha_internal_excitation(int i,int j);
  intpairvec  get_beta_internal_excitation(int i,int j);
  double      get_sign_internal_excitation(int i,int j);

private:
  void        tuning();
  void        read_info();
  void        read_mo_spaces();
  void        read_mo_spaces2();
  void        compute_mo_mappings();
  void        print_info();
  void        print_mo();
  void        free_memory();

  // Model space functions
  void        print_model_space();
  void        build_model_space();
  void        make_internal_excitations();

  /////////////////////////////////////////////////////////////////////////////////////////////////
  // MOInfo variables
  int         root;


  double      scf_energy;
  double      fzcore_energy;

  /////////////////////////////////////////////////////////////////////////////////////////////////
  double      dgemm_timing;
  // In-core/Out-of-core
  double      no_damp_convergence;

  int         nel;
  int         reference;

  // Total number of orbitals in each space
  int         nfocc;                                // Frozen doubly-occupied MOs
  int         nfvir;                                // Frozen external MOs
  int         nactv_docc;                           // Number of active ???
  int         nocc;                                 // Generalized occupied (docc + actv)
  int         nvir;                                 // Generalized virtual (actv + extr)
  int         nall;                                 // Non-frozen MOs (docc + actv + extr)
  int         nextr;                                // Non-frozen external orbitals (extr)

  // Orbitals arrays
  intvec      focc;
  intvec      fvir;
  intvec      occ;
  intvec      vir;
  intvec      all;
  intvec      extr;
  intvec      actv_docc;
  intvec      mopi;

  // Mapping arrays
  intvec      all_to_mo;
  intvec      mo_to_all;
  intvec      orbs_to_mo;
  intvec      focc_to_mo;
  intvec      docc_to_mo;
  intvec      actv_to_mo;
  intvec      extr_to_mo;
  intvec      fvir_to_mo;
  intvec      occ_to_mo;
  intvec      vir_to_mo;
  intvec      mo_to_occ_act;
  intvec      mo_to_act_vir;
  intvec      occ_to_vir;
  intvec      all_to_occ;
  intvec      all_to_vir;
  intvec      actv_to_occ;
  intvec      actv_to_vir;
  intvec      occ_to_actv;
  intvec      vir_to_actv;
  intvec      occ_to_all;
  intvec      extr_to_all;
  boolvec     is_actv_in_occ;
  boolvec     is_actv_in_vir;

  // Symmetry
  intvec      all_sym;

  // Model space
  std::vector<SlaterDeterminant> references;
  std::vector<std::vector<std::vector<std::pair<int,int> > > > alpha_internal_excitations;
  std::vector<std::vector<std::vector<std::pair<int,int> > > > beta_internal_excitations;
  std::vector<std::vector<double> >                            sign_internal_excitations;
  std::vector<int> all_refs;
  std::vector<int> unique_refs;
  std::vector<int> closed_shell_refs;
  std::vector<int> unique_open_shell_refs;
};

}

#endif // _psi_src_lib_libmoinfo_moinfo_h_
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

#include <cstdlib>

#include "psi4/libmoinfo/libmoinfo.h"
#include "psi4/libpsi4util/libpsi4util.h"
#include "psi4/liboptions/liboptions.h"

#include "blas.h"
#include "mp2_ccsd.h"
#include "debugging.h"
#include "matrix.h"
#include "sort.h"

extern FILE* outfile;

namespace psi{ namespace psimrcc{

using namespace std;

MP2_CCSD::MP2_CCSD(SharedWavefunction ref_wfn, Options &options):
        CCManyBody(ref_wfn, options)
{
  triples_type = pt2;
  add_matrices();
}

MP2_CCSD::~MP2_CCSD()
{
}

void MP2_CCSD::compute_mp2_ccsd_energy()
{


  generate_integrals();
  generate_denominators();
  compute_reference_energy();

  build_offdiagonal_F();

  blas->diis_add("t2[oO][vV]{u}","t2_delta[oO][vV]{u}");

  print_method("  MP2");

  outfile->Printf("\n  ------------------------------------------------------------------------------");
  outfile->Printf("\n     MP2      Cycle        Energy            Delta E    DIIS");
  outfile->Printf("\n     MP2                    [Eh]              [Eh]  ");
  outfile->Printf("\n  ------------------------------------------------------------------------------");

  // Start the MP2 cycle
  bool converged = false;
  int  cycle     = 0;
  delta_energy = 0.0;
  current_energy = compute_energy();
  while(!converged){
    outfile->Printf("\n    @MP2      %5d   %20.15f  %11.4e",cycle,current_energy,delta_energy);
    build_mp2_t2_iJaB_amplitudes();
    blas->diis_save_t_amps(cycle);
    blas->diis(cycle,delta_energy,DiisEachCycle);

    blas->solve("t2[oo][vv]{u}  = t2[oO][vV]{u}");
    blas->solve("t2[oo][vv]{u} += #2134# - t2[oO][vV]{u}");
    blas->solve("t2[OO][VV]{u}  = t2[oo][vv]{u}");

    synchronize_amps(); // TODO: make this more efficient
    build_tau();

    current_energy = compute_energy();
    delta_energy   = current_energy - old_energy;
    old_energy = current_energy;

    converged = (fabs(delta_energy) < options_.get_double("E_CONVERGENCE"));

    cycle++;

  }

  outfile->Printf("\n  ------------------------------------------------------------------------------");

  outfile->Printf("\n\n   * MP2@       =%25.15f\n",current_energy);

  // Compute the singlet and triplet MP2 contribution to the energy
  compute_mp2_components();


  print_method("  MP2-CCSD");

  outfile->Printf("\n  ------------------------------------------------------------------------------");
  outfile->Printf("\n     MP2-CCSD Cycle        Energy            Delta E    DIIS");
  outfile->Printf("\n     MP2-CCSD               [Eh]              [Eh]  ");
  outfile->Printf("\n  ------------------------------------------------------------------------------");

  blas->diis_add("t1[o][v]{u}","t1_delta[o][v]{u}");

  // Start the MP2-CCSD cycle
  converged = false;
  cycle     = 0;
  delta_energy = 0.0;
  while(!converged){
    outfile->Printf("\n    @MP2-CCSD %5d   %20.15f  %11.4e",cycle,current_energy,delta_energy);

    // These two go together before updating any other intermediate
    build_F_intermediates();
    build_W_intermediates();
    build_Z_intermediates();

    build_amplitudes();
    blas->diis_save_t_amps(cycle);
    blas->diis(cycle,delta_energy,DiisEachCycle);

    blas->solve("t2[oo][vv]{u}  = t2[oO][vV]{u}");
    blas->solve("t2[oo][vv]{u} += #2134# - t2[oO][vV]{u}");
    blas->solve("t2[OO][VV]{u}  = t2[oo][vv]{u}");
    blas->solve("t1[O][V]{u} = t1[o][v]{u}");

    synchronize_amps();
    build_tau();

    current_energy = compute_energy();

    delta_energy = current_energy-old_energy;
    converged = (fabs(delta_energy) < options_.get_double("E_CONVERGENCE"));
    old_energy=current_energy;

    if(cycle>options_.get_int("MAXITER")){
      outfile->Printf("\n\n\tThe calculation did not converge in %d cycles\n\tQuitting PSIMRCC\n",options_.get_int("MAXITER"));

      exit(1);
    }
    cycle++;

  }
  outfile->Printf("\n  ------------------------------------------------------------------------------");

  outfile->Printf("\n\n   * MP2-CCSD total energy = %25.15f\n",current_energy);

  compute_mp2_ccsd_components();


}

double MP2_CCSD::compute_energy()
{
  // Compute the energy using a simple UCCSD energy expression
  blas->solve("Eaa{u}   = t1[o][v]{u} . fock[o][v]{u}");
  blas->solve("Ebb{u}   = t1[O][V]{u} . fock[O][V]{u}");

  blas->solve("Eaaaa{u} = 1/4 tau[oo][vv]{u} . <[oo]:[vv]>");
  blas->solve("Eabab{u} =     tau[oO][vV]{u} . <[oo]|[vv]>");
  blas->solve("Ebbbb{u} = 1/4 tau[OO][VV]{u} . <[oo]:[vv]>");

  blas->solve("EPT2{u}  = Eaa{u} + Ebb{u} + Eaaaa{u} + Eabab{u} + Ebbbb{u} + ERef{u}");

  return(blas->get_scalar("EPT2",0));
}

void MP2_CCSD::read_mp2_ccsd_integrals()
{
  START_TIMER(1,"Reading the integrals required by MP2-CCSD");

  // CCSort reads the one and two electron integrals
  // and creates the Fock matrices
  sorter = new CCSort(ref_wfn_, out_of_core_sort);

  END_TIMER(1);
}

void MP2_CCSD::build_amplitudes()
{
  // These are required by the t1 amplitude equations
  build_t1_ia_amplitudes();
  build_t1_IA_amplitudes();
  build_t2_iJaB_amplitudes();
  build_t2_ijab_amplitudes();
  build_t2_IJAB_amplitudes();
}

void MP2_CCSD::compute_mp2_components()
{
  blas->solve("Eaaaa{u} = 1/4 tau[oo][vv]{u} . <[oo]:[vv]>");
  blas->solve("Eabab{u} =     tau[oO][vV]{u} . <[oo]|[vv]>");
  blas->solve("Ebbbb{u} = 1/4 tau[OO][VV]{u} . <[oo]:[vv]>");

  double mp2_triplet = blas->get_scalar("Eaaaa",0) + blas->get_scalar("Ebbbb",0);
  double mp2_singlet = blas->get_scalar("Eabab",0);

  outfile->Printf("\n   * MP2 Singlet correlation energy = %20.15f",mp2_singlet);
  outfile->Printf("\n   * MP2 Triplet correlation energy = %20.15f",mp2_triplet);
}

void MP2_CCSD::compute_mp2_ccsd_components()
{
  blas->solve("Eaa{u}   = t1[o][v]{u} . fock[o][v]{u}");
  blas->solve("Ebb{u}   = t1[O][V]{u} . fock[O][V]{u}");

  blas->solve("Eaaaa{u} = 1/4 tau[oo][vv]{u} . <[oo]:[vv]>");
  blas->solve("Eabab{u} =     tau[oO][vV]{u} . <[oo]|[vv]>");
  blas->solve("Ebbbb{u} = 1/4 tau[OO][VV]{u} . <[oo]:[vv]>");

  double mp2_ccsd_singles = blas->get_scalar("Eaa",0) + blas->get_scalar("Ebb",0);
  double mp2_ccsd_triplet = blas->get_scalar("Eaaaa",0) + blas->get_scalar("Ebbbb",0);
  double mp2_ccsd_singlet = blas->get_scalar("Eabab",0);

  outfile->Printf("\n   * MP2-CCSD  Singles                    = %20.15f",mp2_ccsd_singles);
  outfile->Printf("\n   * MP2-CCSD  Singlet correlation energy = %20.15f",mp2_ccsd_singlet);
  outfile->Printf("\n   * MP2-CCSD  Triplet correlation energy = %20.15f\n",mp2_ccsd_triplet);


  /////////////////////////////////
  // Compute the CCSD contribution
  /////////////////////////////////

  // Save the MP2-CCSD Hbar in t2_delta
  blas->solve("t2_delta[oO][vV]{u} = t2_eqns[oO][vV]{u}");

  blas->zero("t2_eqns[oO][vV]{u}");

  // Eliminate the (oa,aa) and (aa,va) blocks from the amplitudes
  if(options_.get_str("MP2_CCSD_METHOD")=="II"){
    blas->expand_spaces("HiJaB[oA][aA]{u}","t2_eqns[oO][vV]{u}");
    blas->expand_spaces("HiJaB[aO][aA]{u}","t2_eqns[oO][vV]{u}");
    blas->expand_spaces("HiJaB[aA][vA]{u}","t2_eqns[oO][vV]{u}");
    blas->expand_spaces("HiJaB[aA][aV]{u}","t2_eqns[oO][vV]{u}");
  }
  // Add the (aa,aa) block from the amplitudes
  blas->expand_spaces("HiJaB[aA][aA]{u}","t2_eqns[oO][vV]{u}");

  // Compute CCSD amplitudes
  blas->solve("t2[oO][vV]{u}  = t2_eqns[oO][vV]{u} / d2[oO][vV]{u}");

  blas->solve("t2_eqns[oo][vv]{u}  = t2_eqns[oO][vV]{u}");
  blas->solve("t2_eqns[oo][vv]{u} += #2134# - t2_eqns[oO][vV]{u}");
  blas->solve("t2[oo][vv]{u}  = t2_eqns[oo][vv]{u} / d2[oo][vv]{u}");

  blas->solve("t2[OO][VV]{u}  = t2[oo][vv]{u}");

  build_tau();

  blas->solve("Eaaaa{u} = 1/4 tau[oo][vv]{u} . <[oo]:[vv]>");
  blas->solve("Eabab{u} =     tau[oO][vV]{u} . <[oo]|[vv]>");
  blas->solve("Ebbbb{u} = 1/4 tau[OO][VV]{u} . <[oo]:[vv]>");

  double ccsd_term_singlet = blas->get_scalar("Eabab",0);
  double ccsd_term_triplet = blas->get_scalar("Eaaaa",0) + blas->get_scalar("Ebbbb",0);


  ////////////////////////////////
  // Compute the MP2 contribution
  ////////////////////////////////

  // Load the MP2-CCSD Hbar from t2_delta
  blas->solve("t2_eqns[oO][vV]{u} = t2_delta[oO][vV]{u}");

  // Eliminate the (oa,aa) and (aa,va) blocks from the amplitudes
  if(options_.get_str("MP2_CCSD_METHOD")=="II"){
    blas->zero("HiJaB[oA][aA]{u}");
    blas->zero("HiJaB[aO][aA]{u}");
    blas->zero("HiJaB[aA][vA]{u}");
    blas->zero("HiJaB[aA][aV]{u}");

    blas->expand_spaces("HiJaB[oA][aA]{u}","t2_eqns[oO][vV]{u}");
    blas->expand_spaces("HiJaB[aO][aA]{u}","t2_eqns[oO][vV]{u}");

    blas->expand_spaces("HiJaB[aA][vA]{u}","t2_eqns[oO][vV]{u}");
    blas->expand_spaces("HiJaB[aA][aV]{u}","t2_eqns[oO][vV]{u}");
  }

  // Eliminate the (aa,aa) block from the amplitudes
  blas->zero("HiJaB[aA][aA]{u}");

  blas->expand_spaces("HiJaB[aA][aA]{u}","t2_eqns[oO][vV]{u}");

  // Compute MP2 amplitudes
  blas->solve("t2[oO][vV]{u}  = t2_eqns[oO][vV]{u} / d2[oO][vV]{u}");

  blas->solve("t2_eqns[oo][vv]{u}  = t2_eqns[oO][vV]{u}");
  blas->solve("t2_eqns[oo][vv]{u} += #2134# - t2_eqns[oO][vV]{u}");
  blas->solve("t2[oo][vv]{u}  = t2_eqns[oo][vv]{u} / d2[oo][vv]{u}");

  blas->solve("t2[OO][VV]{u}  = t2[oo][vv]{u}");

  blas->solve("Eaaaa{u} = 1/4 t2[oo][vv]{u} . <[oo]:[vv]>");
  blas->solve("Eabab{u} =     t2[oO][vV]{u} . <[oo]|[vv]>");
  blas->solve("Ebbbb{u} = 1/4 t2[OO][VV]{u} . <[oo]:[vv]>");


  double mp2_term_singlet = blas->get_scalar("Eabab",0);
  double mp2_term_triplet = blas->get_scalar("Eaaaa",0) + blas->get_scalar("Ebbbb",0);

  outfile->Printf("\n   * MP2  Term Singlet correlation energy = %20.15f",mp2_term_singlet);
  outfile->Printf("\n   * MP2  Term Triplet correlation energy = %20.15f\n",mp2_term_triplet);
  outfile->Printf("\n   * CCSD Term Singlet correlation energy = %20.15f",ccsd_term_singlet);
  outfile->Printf("\n   * CCSD Term Triplet correlation energy = %20.15f",ccsd_term_triplet);
}

}} /* End Namespaces */

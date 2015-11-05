/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */


#include <libqt/qt.h>
#include <libciomr/libciomr.h>
#include <libmints/mints.h>

#include <libfock/soscf.h>
#include <libdiis/diismanager.h>
#include <libdiis/diisentry.h>

#include "ciwave.h"
#include "structs.h"

namespace psi { namespace detci {

/*
** compute_mcscf(); Optimizes MO and CI coefficients
**
** Parameters:
**    options = options object
**    alplist = list of alpha strings
**    betlist = list of beta strings
**
** Returns: none
*/
void CIWavefunction::compute_mcscf()
{

  Parameters_->print_lvl = 0;

  boost::shared_ptr<SOMCSCF> somcscf;
  if (MCSCF_Parameters_->mcscf_type == "DF"){
    transform_dfmcscf_ints(!MCSCF_Parameters_->orbital_so);
    somcscf = boost::shared_ptr<SOMCSCF>(new DFSOMCSCF(jk_, dferi_, AO2SO_, H_));
  }
  else {
    transform_mcscf_ints(!MCSCF_Parameters_->orbital_so);
    somcscf = boost::shared_ptr<SOMCSCF>(new DiskSOMCSCF(jk_, ints_, AO2SO_, H_));
  }

  // We assume some kind of ras here.
  if (Parameters_->wfn != "CASSCF"){
    std::vector<Dimension> ras_spaces;

    // We only have four spaces currently
    for (int nras = 0; nras < 4; nras++){
      Dimension rasdim = Dimension(nirrep_, "RAS" + psi::to_string(nras));
      for (int h = 0; h < nirrep_; h++){
          rasdim[h] = CalcInfo_->ras_opi[nras][h];
      }
      ras_spaces.push_back(rasdim);
    }
    somcscf->set_ras(ras_spaces);

  }

  // Set fzc energy
  SharedMatrix Cfzc = get_orbitals("FZC");
  somcscf->set_frozen_orbitals(Cfzc);


  /// => Start traditional two-step MCSCF <= //
  // Parameters
  int conv = 0;
  double ediff, grad_rms, current_energy;
  double old_energy = CalcInfo_->escf;
  std::string itertype = "Initial CI";
  std::string mcscf_type;
  if (MCSCF_Parameters_->mcscf_type == "DF") mcscf_type = "   @DF-MCSCF";
  else mcscf_type = "      @MCSCF";


  // Setup the DIIS manager
  Dimension dim_oa = get_dimension("DOCC") + get_dimension("ACT");
  Dimension dim_av = get_dimension("VIR") + get_dimension("ACT");

  SharedMatrix x(new Matrix("Rotation Matrix", nirrep_, dim_oa, dim_av));
  SharedMatrix xstep;
  SharedMatrix original_orbs = get_orbitals("ROT");

  boost::shared_ptr<DIISManager> diis_manager(new DIISManager(MCSCF_Parameters_->diis_max_vecs,
                      "MCSCF DIIS", DIISManager::OldestAdded, DIISManager::InCore));
  diis_manager->set_error_vector_size(1, DIISEntry::Matrix, x.get());
  diis_manager->set_vector_size(1, DIISEntry::Matrix, x.get());
  int diis_count = 0;

  // Energy header
  outfile->Printf("\n                             "
                    "Total Energy         Delta E       RMS Grad\n\n");

  // Iterate
  for (int iter=1; iter<(MCSCF_Parameters_->max_iter + 1); iter++){

    // Run CI and set quantities
    diag_h();
    //outfile->Printf("Diag h\n");
    form_opdm();
    //outfile->Printf("opdm\n");
    set_opdm();
    //outfile->Printf("tpdm\n");
    form_tpdm();
    //outfile->Printf("Formed matrices\n");

    current_energy = Process::environment.globals["MCSCF TOTAL ENERGY"];
    ediff = current_energy - old_energy;

    //// Get orbitals for new update
    SharedMatrix Cdocc = get_orbitals("DOCC");
    SharedMatrix Cact  = get_orbitals("ACT");
    SharedMatrix Cvir  = get_orbitals("VIR");
    SharedMatrix actOPDM = get_active_opdm();
    SharedMatrix actTPDM = get_active_tpdm();

    somcscf->update(Cdocc, Cact, Cvir, actOPDM, actTPDM);
    grad_rms = somcscf->gradient_rms();

    outfile->Printf("%s Iter %3d:  % 5.16lf   % 1.5e  % 1.5e %s\n", mcscf_type.c_str(), iter,
                    current_energy, ediff, grad_rms, itertype.c_str());

    if (grad_rms < MCSCF_Parameters_->rms_grad_convergence &&
        (fabs(ediff) < fabs(MCSCF_Parameters_->energy_convergence)) &&
        (iter > 3)){
      outfile->Printf("\n       MCSCF has converged!\n");
      break;
    }
    old_energy = current_energy;

    bool do_so_orbital = (((grad_rms < MCSCF_Parameters_->so_start_grad) &&
                           (ediff < MCSCF_Parameters_->so_start_e) &&
                           (iter >= 2))
                           || (itertype == "SOMCSCF"));

    if (do_so_orbital && MCSCF_Parameters_->orbital_so){
      itertype = "SOMCSCF";
      xstep = somcscf->solve(3, 1.e-10, false);
    }
    else {
      itertype = "APPROX";
      xstep = somcscf->approx_solve();
    }

    // Scale x if needed
    double maxx = 0.0;
    for (int h=0; h<xstep->nirrep(); ++h) {
        for (int i=0; i<xstep->rowspi()[h]; i++) {
            for (int j=0; j<xstep->colspi()[h]; j++) {
                if (fabs(xstep->get(h,i,j)) > maxx) maxx = fabs(xstep->get(h,i,j));
            }
        }
    }

    if (maxx > MCSCF_Parameters_->max_rot){
      xstep->scale((MCSCF_Parameters_->max_rot)/maxx);
    }

    // Add step to overall rotation
    x->add(xstep);

    // Do we add diis?
    if (iter > MCSCF_Parameters_->diis_start){
      diis_manager->add_entry(2, xstep.get(), x.get());
      diis_count++;
    }

    // Do we do diis?
    if ((itertype == "APPROX") && !(diis_count % MCSCF_Parameters_->diis_freq) && (iter > MCSCF_Parameters_->diis_start)){
      diis_manager->extrapolate(1, x.get());
      itertype = "APPROX, DIIS";
    }

    SharedMatrix new_orbs = somcscf->Ck(original_orbs, x);
    set_orbitals("ROT", new_orbs);

    // Transform integrals
    if (MCSCF_Parameters_->mcscf_type == "DF"){
      transform_dfmcscf_ints(!MCSCF_Parameters_->orbital_so);
    }
    else {
      transform_mcscf_ints(!MCSCF_Parameters_->orbital_so);
    }

  }// End MCSCF
  diis_manager->delete_diis_file();
  diis_manager.reset();


}

}} // End namespace

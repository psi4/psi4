// Just scratch


#include <libqt/qt.h>         
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>
#include <libmints/wavefunction.h>
#include <libmints/matrix.h>
#include <libtrans/integraltransform.h>
#include <libtrans/mospace.h>
#include "MCSCF.h"
#include "globaldefs.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace detci {

}
void MCSCF::mo_transform()
{
    std::vector<shared_ptr<MOSpace> > spaces;
    spaces.push_back(MOSpace::all);

    // Build MO spaces
    double *transorbs = init_int_array(CalcInfo.nirreps);
    for (int h=0; h<CalcInfo.nirreps; h++){
        transorbs[h] = CalcInfo.orbs_per_irr[h] - MCSCF_CalcInfo.frozen_docc[h] - MCSCF_CalcInfo.frozen_uocc[h]
    }

    shared_ptr<Wavefunction> wfn(new Wavefunction());
    SharedMatrix c(new Matrix("C", CalcInfo.nirreps, CalcInfo.nso, MCSCF_CalcInfo.frozen_docc));
    SharedMatrix i(new Matrix("I", CalcInfo.nirreps, CalcInfo.nso, transorbs));
    SharedMatrix a(new Matrix("A", CalcInfo.nirreps, CalcInfo.nso, 0));
    SharedMatrix v(new Matrix("V", CalcInfo.nirreps, CalcInfo.nso, MCSCF_CalcInfo.frozen_uocc));
    


    boost::shared_ptr<IntegralTransform> ints(new IntegralTransform(c, i, a, v, spaces,
                                              IntegralTransform::Restricted,
                                              IntegralTransform::IWLOnly,
                                              IntegralTransform::PitzerOrder,
                                              IntegralTransform::OccOnly
                                              false));
    ints->set_keep_iwl_so_ints(true);
    ints->set_print(3);
    ints->initialize();
    ints->transform_tei(MOSpace::all, MOSpace::all, MOSpace::all, MOSpace::all);

    CalcInfo.efzc = ints->get_frozen_core_energy();
}
}} // end namespace psi::detci

  ciwfn->psio()->tocprint(PSIF_OEI);
  // Libtrans
  std::vector<boost::shared_ptr<MOSpace> > spaces;
  spaces.push_back(MOSpace::all);
  outfile->Printf("Here\n");
  boost::shared_ptr<IntegralTransform> ints(new IntegralTransform((static_cast<boost::shared_ptr<Wavefunction> > (ciwfn)),
                                            spaces,
                                            IntegralTransform::Restricted,
                                            IntegralTransform::IWLOnly,
                                            IntegralTransform::PitzerOrder,
                                            IntegralTransform::OccOnly,
                                            false));
  ints->set_keep_iwl_so_ints(true);
  ints->set_keep_dpd_so_ints(true);
  outfile->Printf("Finished libtrans init\n");
  ints->initialize();
  ciwfn->psio()->tocprint(PSIF_OEI);
  outfile->Printf("V exists: %d\n", ciwfn->psio()->tocentry_exists(PSIF_OEI, PSIF_SO_V));
  outfile->Printf("Built OEI\n");
  ints->transform_oei(MOSpace::all, MOSpace::all, "");
  outfile->Printf("Trans OEI done\n");
  ints->transform_tei(MOSpace::all, MOSpace::all, MOSpace::all, MOSpace::all, IntegralTransform::MakeAndKeep);
  outfile->Printf("Finished Transform\n");
  CalcInfo.efzc = ints->get_frozen_core_energy();
  outfile->Printf("Frozen core energy done\n");


    ints->update_orbitals();
    outfile->Printf("Here\n");
    ints->transform_oei(MOSpace::all, MOSpace::all, "");
    outfile->Printf("Trans OEI done\n");
    ints->transform_tei(MOSpace::all, MOSpace::all, MOSpace::all, MOSpace::all, IntegralTransform::ReadAndKeep);
    outfile->Printf("Finished Transform\n");
    CalcInfo.efzc = ints->get_frozen_core_energy();
    outfile->Printf("Frozen core energy done\n");


   // DGAS This will eventually get nuked
   // CalcInfo.nirreps = reference_wavefunction()->nirrep();
   // CalcInfo.nso = reference_wavefunction()->nso();
   // CalcInfo.nmo = reference_wavefunction()->nmo();
   // CalcInfo.nmo = reference_wavefunction()->nmo();
   // CalcInfo.docc = reference_wavefunction()->doccpi(); 
   // CalcInfo.socc = reference_wavefunction()->soccpi(); 
   // CalcInfo.frozen_docc = reference_wavefunction()->frzcpi();
   // CalcInfo.frozen_uocc = reference_wavefunction()->frzvpi();

   // CalcInfo.enuc = reference_wavefunction()->molecule()->nuclear_repulsion_energy();
   // CalcInfo.escf = reference_wavefunction()->reference_energy();
   // CalcInfo.iopen = CalcInfo.nso * (CalcInfo.nso * 2) / 2;
   // CalcInfo.labels = reference_wavefunction()->molecule()->irrep_labels();
   // CalcInfo.orbs_per_irr = reference_wavefunction()->nmopi();
   // CalcInfo.so_per_irr = reference_wavefunction()->nsopi();
   
   //eig_unsrt = reference_wavefunction()->epsilon_a()->to_block_vector(); 

   // DGAS Edit libtrans will catch this next
   // Most test should fail
   //CalcInfo.efzc = 0;

#include<psifiles.h>
#include<libpsio/psio.hpp>
#include"blas.h"
#include"cim.h"
#include<libmints/wavefunction.h>
#include<libmints/matrix.h>
#include<libtrans/integraltransform.h>


using namespace psi;
using namespace boost;

namespace boost {
template<class T> class shared_ptr;
}

namespace psi{

/* 
 * transform integrals
 */
void CIM::TransformIntegrals(int cluster,SharedMatrix Clmo){
  SharedMatrix Cc (new Matrix("fzc",nso,nfzc_));
  SharedMatrix Ci (new Matrix("occ",nso,ndoccact_));
  SharedMatrix Ca (new Matrix("vir",nso,nvirt_));
  SharedMatrix Cv (new Matrix("fzv",nso,nfzv_));

  double**Clmo_pointer = Clmo->pointer();
  double**pointer;
  if (nfzc>0){
     pointer = Cc->pointer();
     for (int i=0; i<nso; i++){
         for (int j=0; j<nfzc_; j++){
             pointer[i][j] = Clmo_pointer[i][j];
         }
         //F_DCOPY(nfzc,&Clmo_pointer[i][0],1,&pointer[i][0],1);
     }
  }
  pointer = Ci->pointer();
  for (int i=0; i<nso; i++){
      for (int j=0; j<ndoccact_; j++){
          pointer[i][j] = Clmo_pointer[i][j+nfzc];
      }
      //F_DCOPY(ndoccact,&Clmo_pointer[i][nfzc],1,&pointer[i][0],1);
  }
  pointer = Ca->pointer();
  for (int i=0; i<nso; i++){
      for (int j=0; j<nvirt_; j++){
          pointer[i][j] = Clmo_pointer[i][j+ndocc];
      }
      //F_DCOPY(nvirt,&Clmo_pointer[i][ndocc],1,&pointer[i][0],1);
  }
  if (nfzv_>0){
     pointer = Cv->pointer();
     for (int i=0; i<nso; i++){
         for (int j=0; j<nfzv_; j++){
             pointer[i][j] = Clmo_pointer[i][j+ndocc+nvirt_];
         }
          //F_DCOPY(nfzv,&Clmo_pointer[i][nmo],1,&pointer[i][0],1);
     }
  }

  std::vector<shared_ptr<MOSpace> > spaces;
  spaces.push_back(MOSpace::all);
  spaces.push_back(MOSpace::occ);
  spaces.push_back(MOSpace::vir);

  /*boost::shared_ptr<Wavefunction> ref =
       Process::environment.reference_wavefunction();

  boost::shared_ptr<Matrix> Cc = ref->Ca_subset("SO","FROZEN_OCC");
  boost::shared_ptr<Matrix> Ci = ref->Ca_subset("SO","ACTIVE_OCC");
  boost::shared_ptr<Matrix> Ca = ref->Ca_subset("SO","ACTIVE_VIR");
  boost::shared_ptr<Matrix> Cv = ref->Ca_subset("SO","FROZEN_VIR");*/

  // so->no integral transformation
  IntegralTransform ints(Cc, Ci, Ca, Cv, spaces, IntegralTransform::Restricted,
           IntegralTransform::IWLOnly, IntegralTransform::QTOrder, IntegralTransform::OccAndVir, true);
  ints.set_keep_dpd_so_ints(1);
  ints.set_keep_iwl_so_ints(1);
  ints.initialize();
  ints.transform_tei(MOSpace::all, MOSpace::all, MOSpace::all, MOSpace::all);
}

}

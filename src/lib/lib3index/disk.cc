#include "3index.h"

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>

#include <psifiles.h>
#include <libpsio/psio.h>
#include <libqt/qt.h>
#include <libciomr/libciomr.h>
#include <libmints/mints.h>

//MKL Header
#ifdef HAVE_MKL
#include <mkl.h>
#endif

//OpenMP Header
//_OPENMP is defined by the compiler if it exists
#ifdef _OPENMP
#include <omp.h>
#endif


using namespace std;
using namespace psi;

namespace psi { 

void DFTensor::disk_tensor(shared_ptr<Matrix> C1, shared_ptr<Matrix> C2, bool useC1, bool useC2, 
    bool apply_fitting, const std::string & entry)
{
}
void DFTensor::form_Qmn_disk()
{
    shared_ptr<Matrix> null(new Matrix("Null", 1, 1));
    disk_tensor(null, null, false, false, true, "Qmn DF Integrals");
}
void DFTensor::form_Qmi_disk(shared_ptr<Matrix> C_act_occ) 
{
    shared_ptr<Matrix> null(new Matrix("Null", 1, 1));
    disk_tensor(C_act_occ, null, true, false, true, "Qmi DF Integrals");
}
void DFTensor::form_Qma_disk(shared_ptr<Matrix> C_act_virt) 
{
    shared_ptr<Matrix> null(new Matrix("Null", 1, 1));
    disk_tensor(C_act_virt, null, true, false, true, "Qma DF Integrals");
}
void DFTensor::form_Qii_disk(shared_ptr<Matrix> C1_act_occ, shared_ptr<Matrix> C2_act_occ) 
{
    disk_tensor(C1_act_occ, C2_act_occ, true, true, true, "Qii DF Integrals");
}
void DFTensor::form_Qia_disk(shared_ptr<Matrix> C_act_occ, shared_ptr<Matrix> C_act_virt) 
{
    disk_tensor(C_act_occ, C_act_virt, true, true, true, "Qia DF Integrals");
}
void DFTensor::form_Qaa_disk(shared_ptr<Matrix> C1_act_virt, shared_ptr<Matrix> C2_act_virt) 
{
    disk_tensor(C1_act_virt, C2_act_virt, true, true, true, "Qaa DF Integrals");
}
void DFTensor::form_Amn_disk()
{
    shared_ptr<Matrix> null(new Matrix("Null", 1, 1));
    disk_tensor(null, null, false, false, false, "Amn DF Integrals");
}
void DFTensor::form_Ami_disk(shared_ptr<Matrix> C_act_occ) 
{
    shared_ptr<Matrix> null(new Matrix("Null", 1, 1));
    disk_tensor(C_act_occ, null, true, false, false, "Ami DF Integrals");
}
void DFTensor::form_Ama_disk(shared_ptr<Matrix> C_act_virt) 
{
    shared_ptr<Matrix> null(new Matrix("Null", 1, 1));
    disk_tensor(C_act_virt, null, true, false, false, "Ama DF Integrals");
}
void DFTensor::form_Aii_disk(shared_ptr<Matrix> C1_act_occ, shared_ptr<Matrix> C2_act_occ) 
{
    disk_tensor(C1_act_occ, C2_act_occ, true, true, false, "Aii DF Integrals");
}
void DFTensor::form_Aia_disk(shared_ptr<Matrix> C_act_occ, shared_ptr<Matrix> C_act_virt) 
{
    disk_tensor(C_act_occ, C_act_virt, true, true, false, "Aia DF Integrals");
}
void DFTensor::form_Aaa_disk(shared_ptr<Matrix> C1_act_virt, shared_ptr<Matrix> C2_act_virt) 
{
    disk_tensor(C1_act_virt, C2_act_virt, true, true, false, "Aaa DF Integrals");
}

}

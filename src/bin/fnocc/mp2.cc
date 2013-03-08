#include"psi4-dec.h"
#include<libmints/vector.h>
#include<libmints/matrix.h>
#include<libmints/wavefunction.h>
#include<libqt/qt.h>
#include<sys/times.h>
#include<libtrans/mospace.h>
#include<libtrans/integraltransform.h>
#include<libiwl/iwl.h>
#include<psifiles.h>
#ifdef _OPENMP
    #include<omp.h>
#else
    #define omp_get_wtime() 0.0
#endif

#include"blas.h"
#include"ccsd.h"

using namespace psi;
namespace psi{ namespace fnocc{
void SortOVOV(struct iwlbuf *Buf,int nfzc,int nfzv,int norbs,int ndoccact,int nvirt);
}}

namespace psi{ namespace fnocc{

/**
  * canonical MP2
  */
void CoupledCluster::MP2(){
    int o = ndoccact;
    int v = nvirt;

    // transform integrals
    fprintf(outfile,"        ==> Transform (OV|OV) integrals <==\n");
    fprintf(outfile,"\n");
    boost::shared_ptr<psi::Wavefunction> wfn = Process::environment.wavefunction();
    std::vector<boost::shared_ptr<MOSpace> > spaces;
    spaces.push_back(MOSpace::occ);
    spaces.push_back(MOSpace::vir);
    boost::shared_ptr<IntegralTransform>
        ints(new IntegralTransform(wfn,
                                   spaces,
                                   IntegralTransform::Restricted,
                                   IntegralTransform::IWLOnly,
                                   IntegralTransform::QTOrder,
                                   IntegralTransform::None,
                                   false));
    ints->set_keep_dpd_so_ints(1);
    ints->set_keep_iwl_so_ints(1);
    ints->initialize();
    ints->transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::occ, MOSpace::vir);
    ints.reset();

    // sort integrals
    fprintf(outfile,"\n");
    fprintf(outfile,"        ==> Sort (OV|OV) integrals <==\n");
    fprintf(outfile,"\n");
    struct iwlbuf Buf;
    iwl_buf_init(&Buf,PSIF_MO_TEI,0.0,1,1);
    SortOVOV(&Buf,nfzc,nfzv,nfzc+nfzv+ndoccact+nvirt,ndoccact,nvirt);
    iwl_buf_close(&Buf,1);


    // energy
    double * v2 = (double*)malloc(o*o*v*v*sizeof(double));

    boost::shared_ptr<PSIO> psio(new PSIO());
    psio->open(PSIF_DCC_IAJB,PSIO_OPEN_OLD);
    psio->read_entry(PSIF_DCC_IAJB,"E2iajb",(char*)&v2[0],o*o*v*v*sizeof(double));
    psio->close(PSIF_DCC_IAJB,0);

    emp2_os = emp2_ss = 0.0;
    for (int i = 0; i < o; i++) {
        double di = eps[i];
        for (int j = 0; j < o; j++) {
            double dij = di + eps[j];
            for (int a = 0; a < v; a++) {
                double dija = dij - eps[a+o];
                for (int b = 0; b < v; b++) {
                    double dijab = dija - eps[b+o];
                    double val1 = v2[i*o*v*v+a*o*v+j*v+b];
                    emp2_os += val1 * val1 / dijab;
                    emp2_ss += (val1 - v2[i*o*v*v+b*o*v+j*v+a]) * val1 / dijab;
                }
            }
        }
    }
    emp2 = emp2_os + emp2_ss;

    fprintf(outfile,"        OS MP2 correlation energy:       %20.12lf\n",emp2_os);
    fprintf(outfile,"        SS MP2 correlation energy:       %20.12lf\n",emp2_ss);
    fprintf(outfile,"        MP2 correlation energy:          %20.12lf\n",emp2);
    fprintf(outfile,"      * MP2 total energy:                %20.12lf\n",emp2+escf);
    fprintf(outfile,"\n");

    free(v2);
}

}} // end of namespaces

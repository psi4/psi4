/**
  * Frozen natural orbitals
  * Eugene DePrince
  * April 2013
  *
  */

#include"psi4-dec.h"
#include<psifiles.h>
#include<libmints/mints.h>
#include<libmints/mintshelper.h>
#include<libmints/wavefunction.h>
#include<libmints/matrix.h>
#include<libtrans/mospace.h>
#include<libtrans/integraltransform.h>
#include<libiwl/iwl.h>
#include"ccsd.h"
#include"blas.h"
#include"frozen_natural_orbitals.h"
#include<libciomr/libciomr.h>
#include<lib3index/dftensor.h>
#include<lib3index/cholesky.h>
#include <libmints/sieve.h>

#include<libdpd/dpd.h>
#define ID(x) ints->DPD_ID(x)


using namespace psi;
using namespace boost;

namespace psi{namespace fnocc{

FrozenNO::FrozenNO(boost::shared_ptr<Wavefunction>wfn,Options&options):
  Wavefunction(options, _default_psio_lib_)
{
    // copy wave function.
    copy(wfn);
    common_init();
}
FrozenNO::~FrozenNO()
{
}

void FrozenNO::common_init() {
    nso = nmo = ndocc = nvirt = nfzc = nfzv = 0;
    for (int h=0; h<nirrep_; h++){
        nfzc   += frzcpi_[h];
        nfzv   += frzvpi_[h];
        nso    += nsopi_[h];
        nmo    += nmopi_[h];
        ndocc  += doccpi_[h];
    }
    ndoccact = ndocc - nfzc;
    nvirt    = nmo - ndocc;

    // quit if number of virtuals is less than number of doubly occupied
    if (nvirt<ndoccact){
       throw PsiException("ndocc must be less than nvirt",__FILE__,__LINE__);
    }

}
// use this function to return the mp2 energy in the full basis.
double FrozenNO::compute_energy(){
  return emp2;
}

/*
 * build natural orbitals
 */
void FrozenNO::ComputeNaturalOrbitals(){
    tstart();

    fflush(outfile);
    fprintf(outfile,"\n\n");
    fprintf(outfile, "        *******************************************************\n");
    fprintf(outfile, "        *                                                     *\n");
    fprintf(outfile, "        *               Frozen Natural Orbitals               *\n");
    fprintf(outfile, "        *                                                     *\n");
    fprintf(outfile, "        *******************************************************\n");
    fprintf(outfile,"\n\n");
    fflush(outfile);

    fprintf(outfile,"        ==> Transform (OV|OV) integrals <==\n");
    fprintf(outfile,"\n");

    std::vector<shared_ptr<MOSpace> > spaces;
    spaces.push_back(MOSpace::occ);
    spaces.push_back(MOSpace::vir);
    boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
    boost::shared_ptr<IntegralTransform> ints(new IntegralTransform(wfn, spaces, IntegralTransform::Restricted,
               IntegralTransform::DPDOnly, IntegralTransform::QTOrder, IntegralTransform::OccAndVir, false));
    ints->set_dpd_id(0);
    ints->set_keep_iwl_so_ints(true);
    ints->set_keep_dpd_so_ints(true);
    ints->initialize();
    ints->transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::occ, MOSpace::vir);

    fprintf(outfile,"\n");
    fprintf(outfile,"        ==> Build MP2 amplitudes, OPDM, and NOs <==\n");
    fprintf(outfile,"\n");

    dpdbuf4 amps1,amps2;
    boost::shared_ptr<PSIO> psio = _default_psio_lib_;
    psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    // Use the IntegralTransform object's DPD instance, for convenience
    dpd_set_default(ints->get_dpd_id());

    // orbital energies:
    int numAOcc = 0, numBOcc = 0, numAVir = 0, numBVir = 0;
    int aOccCount = 0, bOccCount = 0, aVirCount = 0, bVirCount = 0;
    int * aOccOrbsPI = new int[nirrep_];
    int * bOccOrbsPI = new int[nirrep_];
    int * aVirOrbsPI = new int[nirrep_];
    int * bVirOrbsPI = new int[nirrep_];
    for(int h = 0; h < nirrep_; ++h){
        aOccOrbsPI[h] = doccpi_[h] + soccpi_[h] - frzcpi_[h];
        bOccOrbsPI[h] = doccpi_[h] - frzcpi_[h];
        aVirOrbsPI[h] = nmopi_[h] - doccpi_[h] - soccpi_[h] - frzvpi_[h];
        bVirOrbsPI[h] = nmopi_[h] - doccpi_[h] - frzvpi_[h];
        numAOcc += aOccOrbsPI[h];
        numBOcc += bOccOrbsPI[h];
        numAVir += aVirOrbsPI[h];
        numBVir += bVirOrbsPI[h];
    }
    double * aOccEvals = new double[numAOcc];
    double * bOccEvals = new double[numBOcc];
    double * aVirEvals = new double[numAVir];
    double * bVirEvals = new double[numBVir];

    boost::shared_ptr<Vector> epsA = reference_wavefunction_->epsilon_a();
    boost::shared_ptr<Vector> epsB = reference_wavefunction_->epsilon_b();
    for(int h = 0; h < nirrep_; ++h){
        for(int a = frzcpi_[h]; a < doccpi_[h] + soccpi_[h]; ++a) aOccEvals[aOccCount++] = epsA->get(h, a);
        for(int b = frzcpi_[h]; b < doccpi_[h]; ++b)              bOccEvals[bOccCount++] = epsB->get(h, b);
        for(int a = doccpi_[h] + soccpi_[h]; a < nmopi_[h]; ++a)  aVirEvals[aVirCount++] = epsA->get(h, a);
        for(int b = doccpi_[h]; b < nmopi_[h]; ++b)               bVirEvals[bVirCount++] = epsB->get(h, b);
    }

    dpd_buf4_init(&amps1, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");
    dpd_buf4_sort(&amps1, PSIF_LIBTRANS_DPD, prqs, ID("[O,O]"), ID("[V,V]"), "MO Ints <OO|VV>");

    // T(ijab) -> T(jiab)
    dpd_buf4_sort(&amps1, PSIF_LIBTRANS_DPD, prqs, ID("[O,O]"), ID("[V,V]"), "Tijab <OO|VV>");
    dpd_buf4_sort(&amps1, PSIF_LIBTRANS_DPD, rpqs, ID("[O,O]"), ID("[V,V]"), "Tjiab <OO|VV>");
    dpd_buf4_close(&amps1);

    // only worry about alpha-beta
    dpd_buf4_init(&amps1, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
    dpd_buf4_init(&amps2, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),ID("[O,O]"), ID("[V,V]"), 0, "Tjiab <OO|VV>");

    double emp2_os = 0.0;
    double emp2_ss = 0.0;
    for(int h = 0; h < nirrep_; ++h){

        dpd_buf4_mat_irrep_init(&amps1, h);
        dpd_buf4_mat_irrep_rd(&amps1, h);

        dpd_buf4_mat_irrep_init(&amps2, h);
        dpd_buf4_mat_irrep_rd(&amps2, h);

        for(int ij = 0; ij < amps1.params->rowtot[h]; ++ij){
            int i = amps1.params->roworb[h][ij][0];
            int j = amps1.params->roworb[h][ij][1];

            for(int ab = 0; ab < amps1.params->coltot[h]; ++ab){
                int a = amps1.params->colorb[h][ab][0];
                int b = amps1.params->colorb[h][ab][1];
                double val1 = amps1.matrix[h][ij][ab];
                double val2 = amps2.matrix[h][ij][ab];
                double denom = aOccEvals[i] + bOccEvals[j] - aVirEvals[a] - bVirEvals[b];

                emp2_os += val1 * val1 / denom;
                emp2_ss += val1 * ( val1 - val2 ) / denom;
            }
        }
        dpd_buf4_mat_irrep_close(&amps1, h);
        dpd_buf4_mat_irrep_close(&amps2, h);
    }
    dpd_buf4_close(&amps1);

    double escf = Process::environment.globals["SCF TOTAL ENERGY"];
    Process::environment.globals["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = emp2_os;
    Process::environment.globals["MP2 SAME-SPIN CORRELATION ENERGY"] = emp2_ss;
    Process::environment.globals["MP2 CORRELATION ENERGY"] = emp2_os + emp2_ss;
    Process::environment.globals["MP2 TOTAL ENERGY"] = emp2_os + emp2_ss + escf;

    // build amps1(ij,ab) = 2*T(ij,ab) - T(ji,ab)
    dpd_buf4_init(&amps1, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
    dpd_buf4_scm(&amps1, 2.0);
    dpd_buf4_axpy(&amps2, &amps1, -1.0);

    dpd_buf4_close(&amps1);
    dpd_buf4_close(&amps2);

    fprintf(outfile,"        OS MP2 correlation energy:       %20.12lf\n",emp2_os);
    fprintf(outfile,"        SS MP2 correlation energy:       %20.12lf\n",emp2_ss);
    fprintf(outfile,"        MP2 correlation energy:          %20.12lf\n",emp2_os+emp2_ss);
    fprintf(outfile,"      * MP2 total energy:                %20.12lf\n",emp2_os+emp2_ss+escf);
    fprintf(outfile,"\n");

    // scale amps by denominator
    dpd_buf4_init(&amps2, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
    dpd_buf4_init(&amps1, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),ID("[O,O]"), ID("[V,V]"), 0, "Tijab <OO|VV>");
    for(int h = 0; h < nirrep_; ++h){

        dpd_buf4_mat_irrep_init(&amps1, h);
        dpd_buf4_mat_irrep_rd(&amps1, h);

        dpd_buf4_mat_irrep_init(&amps2, h);
        dpd_buf4_mat_irrep_rd(&amps2, h);

        for(int ij = 0; ij < amps1.params->rowtot[h]; ++ij){
            int i = amps1.params->roworb[h][ij][0];
            int j = amps1.params->roworb[h][ij][1];
            for(int ab = 0; ab < amps1.params->coltot[h]; ++ab){
                int a = amps1.params->colorb[h][ab][0];
                int b = amps1.params->colorb[h][ab][1];
                double denom = aOccEvals[i] + bOccEvals[j] - aVirEvals[a] - bVirEvals[b];
                amps1.matrix[h][ij][ab] /= ( aOccEvals[i] + bOccEvals[j] - aVirEvals[a] - bVirEvals[b]);
                amps2.matrix[h][ij][ab] /= ( aOccEvals[i] + bOccEvals[j] - aVirEvals[a] - bVirEvals[b]);
            }
        }

        dpd_buf4_mat_irrep_wrt(&amps1, h);
        dpd_buf4_mat_irrep_close(&amps1, h);
        dpd_buf4_mat_irrep_wrt(&amps2, h);
        dpd_buf4_mat_irrep_close(&amps2, h);
    }
    dpd_buf4_close(&amps1);
    dpd_buf4_close(&amps2);

    // build virtual-virtual block of opdm: sum(ijc) 2.0 * [ 2 t(ij,ac) - t(ji,ac) ] * t(ij,bc)
    dpdfile2 Dab;
    dpd_file2_init(&Dab, PSIF_LIBTRANS_DPD,    0, 1, 1, "Dab");
    dpd_buf4_init(&amps2, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
    dpd_buf4_init(&amps1, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),ID("[O,O]"), ID("[V,V]"), 0, "Tijab <OO|VV>");
    dpd_contract442(&amps1, &amps2, &Dab, 3, 3, 2.0, 0.0);
    dpd_buf4_close(&amps1);
    dpd_buf4_close(&amps2);
    dpd_file2_close(&Dab);

    // diagonalize virtual-virtual block of opdm
    int symmetry = Ca_->symmetry();
    boost::shared_ptr<Matrix> D (new Matrix("Dab",nirrep_,aVirOrbsPI,aVirOrbsPI,symmetry));

    dpd_file2_init(&Dab, PSIF_LIBTRANS_DPD,    0, 1, 1, "Dab");
    dpd_file2_mat_init(&Dab);
    dpd_file2_mat_rd(&Dab);
    for (int h = 0; h < nirrep_; h++) {
        int v = Dab.params->rowtot[h];

        double ** mat = D->pointer(h);
        for (int a = 0; a < v; a++) {
            for (int b = 0; b < v; b++) {
                mat[a][b] = Dab.matrix[h][a][b];
            }
        }
    }
    dpd_file2_close(&Dab);

    // done with dpd and ints ... reset
    psio->close(PSIF_LIBTRANS_DPD, 1);
    ints.reset();

    boost::shared_ptr<Matrix> eigvec (new Matrix("Dab eigenvectors",nirrep_,aVirOrbsPI,aVirOrbsPI,symmetry));
    boost::shared_ptr<Vector> eigval (new Vector("Dab eigenvalues",nirrep_,aVirOrbsPI));
    D->diagonalize(eigvec,eigval,descending);

    // overwrite ao/mo C matrix with ao/no transformation
    boost::shared_ptr<Matrix> temp (new Matrix("temp",nirrep_,nsopi_,aVirOrbsPI,symmetry));
    for (int h = 0; h < nirrep_; h++) {

        int v = aVirOrbsPI[h];

        double ** c_newv = eigvec->pointer(h);
        double ** c_oldv = Ca_->pointer(h);
        double ** tp     = temp->pointer(h);

        for (int mu = 0; mu < nsopi_[h]; mu++) {
            for (int a = 0; a < v; a++) {
                double dum = 0.0;
                for (int b = 0; b < v; b++) {
                    dum += c_oldv[mu][doccpi_[h] + b] * c_newv[b][a];
                }
                tp[mu][a] = dum;
            }
        }
        for (int mu = 0; mu < nsopi_[h]; mu++) {
            for (int a = 0; a < v; a++) {
                c_oldv[mu][doccpi_[h]+a] = tp[mu][a];
            }
        }
    }

    // determine how many orbitals will be retained
    double cutoff = options_.get_double("OCC_TOLERANCE");
    int * newVirOrbsPI = new int[nirrep_];
    for (int h = 0; h < nirrep_; h++) {
        newVirOrbsPI[h] = 0;
        double * vec = eigval->pointer(h);
        for (int a = 0; a < aVirOrbsPI[h]; a++) {
            if ( vec[a] > cutoff ) newVirOrbsPI[h]++;
        }
    }

    fprintf(outfile,"        Cutoff for significant NO occupancy: %5.3le\n",cutoff);
    fprintf(outfile,"\n");
    fprintf(outfile,"        No. virtuals per irrep (original):  [");
    for (int h = 0; h < nirrep_; h++) fprintf(outfile,"%4i",aVirOrbsPI[h]);
    fprintf(outfile," ]\n");
    fprintf(outfile,"        No. virtuals per irrep (truncated): [");
    for (int h = 0; h < nirrep_; h++) fprintf(outfile,"%4i",newVirOrbsPI[h]);
    fprintf(outfile," ]\n");
    fprintf(outfile,"\n");

    int nvirt_no = 0;
    int nvirt    = 0;
    for (int h = 0; h < nirrep_; h++) {
        nvirt    += aVirOrbsPI[h];
        nvirt_no += newVirOrbsPI[h];
    }
    fprintf(outfile,"        Retaining %i of %i virtual orbitals.\n",nvirt_no,nvirt);
    fprintf(outfile,"\n");

    // transform Fock matrix to truncated NO basis

    boost::shared_ptr<Matrix> Fab (new Matrix("Fab(NO)",nirrep_,newVirOrbsPI,newVirOrbsPI,symmetry));
    for (int h = 0; h < nirrep_; h++) {
        int o    = doccpi_[h];
        int vnew = newVirOrbsPI[h];
        int vold = aVirOrbsPI[h];
        double ** Fnew = Fab->pointer(h);
        double ** cp   = eigvec->pointer(h);
        double * Fold  = epsA->pointer(h);
        for (int a = 0; a < vnew; a++) {
            for (int b = 0; b < vnew; b++) {
                double dum = 0.0;
                for (int c = 0; c < vold; c++) {
                    dum += Fold[o+c] * cp[c][a] * cp[c][b];
                }
                Fnew[a][b] = dum;
            }
        }
    }

    // semicanonicalize orbitals:
    boost::shared_ptr<Matrix> eigvecF (new Matrix("Fab eigenvectors",nirrep_,newVirOrbsPI,newVirOrbsPI,symmetry));
    boost::shared_ptr<Vector> eigvalF (new Vector("Fab eigenvalues",nirrep_,newVirOrbsPI));
    Fab->diagonalize(eigvecF,eigvalF);

    // overwrite ao/no C matrix with ao/semicanonical no transformation:
    for (int h = 0; h < nirrep_; h++) {

        int v = newVirOrbsPI[h];

        double ** c_newv = eigvecF->pointer(h);
        double ** c_oldv = Ca_->pointer(h);
        double ** tp     = temp->pointer(h);

        for (int mu = 0; mu < nsopi_[h]; mu++) {
            for (int a = 0; a < v; a++) {
                double dum = 0.0;
                for (int b = 0; b < v; b++) {
                    dum += c_oldv[mu][doccpi_[h] + b] * c_newv[b][a];
                }
                tp[mu][a] = dum;
            }
        }
        for (int mu = 0; mu < nsopi_[h]; mu++) {
            for (int a = 0; a < v; a++) {
                c_oldv[mu][doccpi_[h]+a] = tp[mu][a];
            }
        }
    }

    // adjust number of frozen virtual orbitals:
    for (int h = 0; h < nirrep_; h++) {
        frzvpi_[h] += aVirOrbsPI[h] - newVirOrbsPI[h];
    }

    // put modified orbital energies back into epsilon_a
    boost::shared_ptr<Vector> eps = epsilon_a_;
    for (int h = 0; h < nirrep_; h++) {
        double * epsp = eps->pointer(h);
        double * eigp = eigvalF->pointer(h);
        for (int a = 0; a < newVirOrbsPI[h]; a++) {
            epsp[doccpi_[h]+a] = eigp[a];
        }
    }

    // free memory
    delete newVirOrbsPI;

    delete aOccOrbsPI;
    delete bOccOrbsPI;
    delete aVirOrbsPI;
    delete bVirOrbsPI;
    delete aOccEvals;
    delete bOccEvals;
    delete aVirEvals;
    delete bVirEvals;

    tstop();
}

// DF FNO class members

DFFrozenNO::DFFrozenNO(boost::shared_ptr<Wavefunction>wfn,Options&options):
  FrozenNO(wfn,options)
{
}
DFFrozenNO::~DFFrozenNO()
{
}

// use this function to return the mp2 energy in the full basis.  
// this isn't used anymore ... i'm reluctant to remove it, though
double DFFrozenNO::compute_energy(){
  return emp2;
}

void DFFrozenNO::ThreeIndexIntegrals() {

  fprintf(outfile,"  ==> 3-index integrals <==\n");
  fprintf(outfile,"\n");

  long int o = ndoccact;
  long int v = nvirt;
  long int nQ;

  // 1.  read scf 3-index integrals from disk

  // get ntri from sieve
  boost::shared_ptr<ERISieve> sieve (new ERISieve(basisset_, options_.get_double("INTS_TOLERANCE")));
  const std::vector<std::pair<int, int> >& function_pairs = sieve->function_pairs();
  long int ntri = function_pairs.size();

  // read integrals that were written to disk in the scf
  long int nQ_scf = Process::environment.globals["NAUX (SCF)"];
  if ( options_.get_str("SCF_TYPE") == "DF" ) {
      boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
      boost::shared_ptr<BasisSet> auxiliary = BasisSet::construct(parser, molecule(), "DF_BASIS_SCF");
      nQ_scf = auxiliary->nbf();
      Process::environment.globals["NAUX (SCF)"] = nQ_scf;
  }

  boost::shared_ptr<Matrix> Qmn = SharedMatrix(new Matrix("Qmn Integrals",nQ_scf,ntri));
  double** Qmnp = Qmn->pointer();
  boost::shared_ptr<PSIO> psio(new PSIO());
  psio->open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DFSCF_BJ, "(Q|mn) Integrals", (char*) Qmnp[0], sizeof(double) * ntri * nQ_scf);
  psio->close(PSIF_DFSCF_BJ,1);

  // unpack and write again in my format
  boost::shared_ptr<Matrix>L (new Matrix("3-index ERIs (SCF)", nQ_scf , nso*nso));
  double ** Lp = L->pointer();
  for (long int mn = 0; mn < ntri; mn++) {
      long int m = function_pairs[mn].first;
      long int n = function_pairs[mn].second;
      for (long int P = 0; P < nQ_scf; P++) {
          Lp[P][m*nso+n] = Qmnp[P][mn];
          Lp[P][n*nso+m] = Qmnp[P][mn];
      }
  }
  psio->open(PSIF_DCC_QSO,PSIO_OPEN_NEW);
  psio->write_entry(PSIF_DCC_QSO,"Qso SCF",(char*)&Lp[0][0],nQ_scf*nso*nso*sizeof(double));
  psio->close(PSIF_DCC_QSO,1);
  L.reset();
  Qmn.reset();

  // for DFCC, assume that the DF basis differs between the SCF and CC (TODO generalize)
  if ( ( options_.get_str("DF_BASIS_CC") != "CHOLESKY" ) ){

      boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
      boost::shared_ptr<BasisSet> auxiliary = BasisSet::construct(parser, molecule(), "DF_BASIS_CC");

      boost::shared_ptr<DFTensor> DF (new DFTensor(basisset(),auxiliary,Ca(),ndocc,nvirt+nfzv,ndoccact,nvirt,options_));
      nQ = auxiliary->nbf();
      boost::shared_ptr<Matrix> tmp = DF->Qso();
      double ** Qso = tmp->pointer();

      // write Qso to disk
      psio->open(PSIF_DCC_QSO,PSIO_OPEN_OLD);
      psio->write_entry(PSIF_DCC_QSO,"Qso CC",(char*)&Qso[0][0],nQ*nso*nso*sizeof(double));
      psio->close(PSIF_DCC_QSO,1);
      fprintf(outfile,"        Number of auxiliary functions:       %5li\n",nQ);

      // stick nQ in process environment so ccsd can know it
      Process::environment.globals["NAUX (CC)"] = (double)nQ;
  }else{

      // Cholesky

      // read integrals from disk if they were generated in the SCF
      if ( options_.get_str("SCF_TYPE") == "CD" ) {
          fprintf(outfile,"        Reading Cholesky vectors from disk ...\n");
          nQ = Process::environment.globals["NAUX (SCF)"];
          fprintf(outfile,"        Cholesky decomposition threshold: %8.2le\n", options_.get_double("CHOLESKY_TOLERANCE"));
          fprintf(outfile,"        Number of Cholesky vectors:          %5li\n",nQ);

          // ntri comes from sieve above
          boost::shared_ptr<Matrix> Qmn = SharedMatrix(new Matrix("Qmn Integrals",nQ,ntri));
          double** Qmnp = Qmn->pointer();
          // TODO: use my 3-index integral file in SCF for DFCC jobs
          psio->open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
          psio->read_entry(PSIF_DFSCF_BJ, "(Q|mn) Integrals", (char*) Qmnp[0], sizeof(double) * ntri * nQ);
          psio->close(PSIF_DFSCF_BJ,1);

          boost::shared_ptr<Matrix>L (new Matrix("CD Integrals", nQ , nso*nso));
          double ** Lp = L->pointer();
          for (long int mn = 0; mn < ntri; mn++) {
              long int m = function_pairs[mn].first;
              long int n = function_pairs[mn].second;
              for (long int P = 0; P < nQ; P++) {
                  Lp[P][m*nso+n] = Qmnp[P][mn];
                  Lp[P][n*nso+m] = Qmnp[P][mn];
              }
          }
          psio->open(PSIF_DCC_QSO,PSIO_OPEN_OLD);
          psio->write_entry(PSIF_DCC_QSO,"Qso CC",(char*)&Lp[0][0],nQ*nso*nso*sizeof(double));
          psio->close(PSIF_DCC_QSO,1);
      }else {

          // generate Cholesky 3-index integrals
          fprintf(outfile,"        Generating Cholesky vectors ...\n");
          boost::shared_ptr<BasisSet> primary = basisset();
          boost::shared_ptr<IntegralFactory> integral (new IntegralFactory(primary,primary,primary,primary));
          double tol = options_.get_double("CHOLESKY_TOLERANCE");
          boost::shared_ptr<CholeskyERI> Ch (new CholeskyERI(boost::shared_ptr<TwoBodyAOInt>(integral->eri()),0.0,tol,Process::environment.get_memory()));
          Ch->choleskify();
          nQ  = Ch->Q();
          boost::shared_ptr<Matrix> L = Ch->L();
          double ** Lp = L->pointer();

          // write Qso to disk
          psio->open(PSIF_DCC_QSO,PSIO_OPEN_OLD);
          psio->write_entry(PSIF_DCC_QSO,"Qso CC",(char*)&Lp[0][0],nQ*nso*nso*sizeof(double));
          psio->close(PSIF_DCC_QSO,1);
          fprintf(outfile,"        Cholesky decomposition threshold: %8.2le\n", tol);
          fprintf(outfile,"        Number of Cholesky vectors:          %5li\n",nQ);

      }

      // stick nQ in process environment so ccsd can know it
      Process::environment.globals["NAUX (CC)"] = (double)nQ;
  }
  fprintf(outfile,"\n");

}

/*
 * build natural orbitals and transform TEIs
 */
void DFFrozenNO::ComputeNaturalOrbitals(){

  fflush(outfile);
  fprintf(outfile, "  ==> Frozen Natural Orbitals <==\n");
  fprintf(outfile,"\n");
  fflush(outfile);

  long int o      = ndoccact;
  long int v      = nvirt;
  long int nQ     = Process::environment.globals["NAUX (CC)"];
  long int nQ_scf = Process::environment.globals["NAUX (SCF)"];
  long int memory = Process::environment.get_memory();

  if ( memory < 8L*(3L*nso*nso+nso*nso*nQ+o*v*nQ) ) {
      throw PsiException("not enough memory (fno)",__FILE__,__LINE__);
  }

  boost::shared_ptr<PSIO> psio(new PSIO());

  // read in 3-index integrals specific to the CC method:
  double * tmp2 = (double*)malloc(nso*nso*nQ * sizeof(double));
  psio->open(PSIF_DCC_QSO,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_QSO,"Qso CC",(char*)&tmp2[0],nQ*nso*nso*sizeof(double));
  psio->close(PSIF_DCC_QSO,1);

  // transform Qso -> Qov:
  TransformQ(nQ,tmp2);
  double * Qov = (double*)malloc(o*v*nQ*sizeof(double));
  F_DCOPY(o*v*nQ,tmp2,1,Qov,1);
  free(tmp2);

  if ( memory < 8L*(o*o*v*v+o*v*nQ) ) {
      throw PsiException("not enough memory (fno)",__FILE__,__LINE__);
  }

  // allocate memory for a couple of buffers
  double * amps2 = (double*)malloc(o*o*v*v*sizeof(double));

  // build (ia|jb) integrals
  F_DGEMM('n','t',o*v,o*v,nQ,1.0,Qov,o*v,Qov,o*v,0.0,amps2,o*v);
  free(Qov);

  if ( memory < 16L*o*o*v*v ) {
      throw PsiException("not enough memory (fno)",__FILE__,__LINE__);
  }

  double * amps1 = (double*)malloc(o*o*v*v*sizeof(double));

  nvirt_no = nvirt;

  boost::shared_ptr<Vector> eps_test = epsilon_a();
  double * tempeps = eps_test->pointer();
  double * F       = tempeps + nfzc;
  double * Dab     = (double*)malloc(v*v*sizeof(double));
  double * temp    = (double*)malloc(nso*v*sizeof(double));
  double * newFock = (double*)malloc(v*v*sizeof(double));
  double * neweps  = (double*)malloc(nvirt_no*sizeof(double));

  // build mp2 amplitudes for mp2 density
  long int ijab = 0;
  emp2 = 0.0;
  double emp2_os = 0.0;
  double emp2_ss = 0.0;
  for (long int a=o; a<o+v; a++){
      double da = F[a];
      for (long int b=o; b<o+v; b++){
          double dab = da + F[b];
          for (long int i=0; i<o; i++){
              double dabi = dab - F[i];
              for (long int j=0; j<o; j++){
                  long int iajb = i*v*v*o+(a-o)*v*o+j*v+(b-o);
                  double dijab = dabi-F[j];
                  amps1[ijab++] = - amps2[iajb]/dijab;
                  emp2_os -= amps2[iajb] * amps2[iajb]/dijab;
                  emp2_ss -= amps2[iajb] * (amps2[iajb] - amps2[j*o*v*v+(a-o)*o*v+i*v+(b-o)])/dijab;
              }
          }
      }
  }
  emp2 = emp2_os + emp2_ss;

  fprintf(outfile,"        Doubles contribution to MP2 energy in full space: %20.12lf\n",emp2);
  fprintf(outfile,"\n");

  Process::environment.globals["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = emp2_os;
  Process::environment.globals["MP2 SAME-SPIN CORRELATION ENERGY"] = emp2_ss;
  Process::environment.globals["MP2 CORRELATION ENERGY"] = emp2;
  Process::environment.globals["MP2 TOTAL ENERGY"] = emp2 + Process::environment.globals["SCF TOTAL ENERGY"];

  ijab = 0;
  for (long int a=o; a<o+v; a++){
      for (long int b=o; b<o+v; b++){
          for (long int i=0; i<o; i++){
              for (long int j=0; j<o; j++){
                  long int ijba = (b-o)*o*o*v+(a-o)*o*o+i*o+j;
                  amps2[ijab] = 2.0*amps1[ijab] - amps1[ijba];
                  ijab++;
              }
          }
      }
  }

  // build ab block of the density:
  F_DGEMM('t','n',v,v,v*o*o,2.0,amps1,v*o*o,amps2,v*o*o,0.0,Dab,v);

  // diagonalize Dab
  double*eigvalDab=(double*)malloc(v*sizeof(double));
  Diagonalize(v,Dab,eigvalDab);

  // reorder transformation matrix:
  for (long int i=0; i<v; i++){
      F_DCOPY(v,Dab+(v-1-i)*v,1,temp+i*v,1);
  }

  // establish cutoff for frozen virtuals
  double cutoff = options_.get_double("OCC_TOLERANCE");
  nvirt_no = 0;
  for (long int i=0; i<v; i++) if (eigvalDab[i]>cutoff) nvirt_no++;

  fprintf(outfile,"        Cutoff for significant NO occupancy: %5.3le\n",cutoff);
  fprintf(outfile,"\n");
  fprintf(outfile,"        Number of virtual orbitals in original space:  %5li\n",v);
  fprintf(outfile,"        Number of virtual orbitals in truncated space: %5li\n",nvirt_no);
  fprintf(outfile,"\n");

  // transform Fock matrix to MP2 NO basis
  memset((void*)newFock,'\0',v*v*sizeof(double));
  F_DCOPY(v,F+ndoccact,1,newFock,v+1);
  F_DGEMM('n','n',v,nvirt_no,v,1.0,newFock,v,temp,v,0.0,Dab,v);
  F_DGEMM('t','n',nvirt_no,nvirt_no,v,1.0,temp,v,Dab,v,0.0,newFock,nvirt_no);

  // diagonalize new Fock matrix for semi-canonical orbitals
  Diagonalize(nvirt_no,newFock,neweps);

  // construct full mo -> no transformation matrix
  F_DGEMM('n','n',v,nvirt_no,nvirt_no,1.0,temp,v,newFock,nvirt_no,0.0,Dab,v);

  // put orbital energies back in F - doesn't matter in this implementation
  F_DCOPY(nvirt_no,neweps,1,F+ndoccact,1);

  // free memory before using libtrans
  free(temp);
  free(neweps);
  free(amps1);
  free(amps2);
  free(eigvalDab);
  free(newFock);

  // change number of frozen virtual orbitals
  frzvpi_[0] = nfzv+(nvirt-nvirt_no);

  // modify c matrix
  ModifyCa(Dab);

  free(Dab);
}

void DFFrozenNO::ModifyCa(double*Dab){

  long int v = nvirt;

  boost::shared_ptr<psi::Wavefunction> ref = Process::environment.wavefunction();

  boost::shared_ptr<Matrix> Caomo = ref->Ca();

  double**Capointer = Caomo->pointer();

  // so->no
  double*temp = (double*)malloc(nso*nvirt_no*sizeof(double));
  for (long int i=0; i<nso; i++){
      for (long int j=0; j<nvirt_no; j++){
          double dum = 0.0;
          for (long int k=0; k<v; k++){
              dum += Capointer[i][ndocc+k] * Dab[j*v+k];
          }
          temp[i*nvirt_no+j] = dum;
      }
  }
  for (long int i=0; i<nso; i++){
      for (long int j=0; j<nvirt_no; j++){
          Capointer[i][ndocc+j] = temp[i*nvirt_no+j];
      }
  }
  free(temp);
}
void DFFrozenNO::ModifyCa_occ(double*Dij){

  long int o = ndoccact;

  boost::shared_ptr<psi::Wavefunction> ref = Process::environment.wavefunction();

  boost::shared_ptr<Matrix> Caomo = ref->Ca();

  double**Capointer = Caomo->pointer();

  // so->no
  double*temp = (double*)malloc(nso*o*sizeof(double));
  for (long int i=0; i<nso; i++){
      for (long int j=0; j<o; j++){
          double dum = 0.0;
          for (long int k=0; k<o; k++){
              dum += Capointer[i][nfzc+k] * Dij[j*o+k];
          }
          temp[i*o+j] = dum;
      }
  }
  for (long int i=0; i<nso; i++){
      for (long int j=0; j<o; j++){
          Capointer[i][nfzc+j] = temp[i*o+j];
      }
  }
  free(temp);
}

void DFFrozenNO::TransformQ(long int nQ,double*Qso) {

    long int o = ndoccact;
    long int v = nvirt;

    double ** Cap = Ca()->pointer();
    double * tmp = (double*)malloc(nso*o*nQ*sizeof(double));

    for (int q = 0; q < nQ; q++) {
        for (int mu = 0; mu < nso; mu++) {
            for (int i = 0; i < o; i++) {
                double dum = 0.0;
                for (int nu = 0; nu < nso; nu++) {
                    dum += Qso[q*nso*nso+mu*nso+nu] * Cap[nu][nfzc+i];
                }
                tmp[q*nso*o+i*nso+mu] = dum;
            }
        }
    }
    for (int q = 0; q < nQ; q++) {
        for (int i = 0; i < o; i++) {
            for (int a = 0; a < v; a++) {
                double dum = 0.0;
                for (int nu = 0; nu < nso; nu++) {
                    dum += tmp[q*nso*o+i*nso+nu] * Cap[nu][nfzc+o+a];
                }
                Qso[q*v*o+i*v+a] = dum;
            }
        }
    }

    free(tmp);
}

void DFFrozenNO::BuildFock(long int nQ,double*Qso,double*F) {

    double ** Cap = Ca()->pointer();

    // Transform Qso to MO basis:
    double * tmp = (double*)malloc(nso*nso*nQ*sizeof(double));
    F_DCOPY(nso*nso*nQ,Qso,1,tmp,1);
    F_DGEMM('n','n',nmo,nso*nQ,nso,1.0,&Cap[0][0],nmo,tmp,nso,0.0,Qso,nmo);
    #pragma omp parallel for schedule (static)
    for (long int q = 0; q < nQ; q++) {
        for (long int mu = 0; mu < nso; mu++) {
            F_DCOPY(nmo,Qso+q*nso*nmo+mu*nmo,1,tmp+q*nso*nmo+mu,nmo);
        }
    }
    F_DGEMM('n','n',nmo,nmo*nQ,nso,1.0,&Cap[0][0],nmo,tmp,nso,0.0,Qso,nmo);

    // build Fock matrix:

    // transform H
    // one-electron integrals
    boost::shared_ptr<MintsHelper> mints(new MintsHelper());
    SharedMatrix H = mints->so_kinetic();
    H->add(mints->so_potential());

    long int max = nQ > nso*nso ? nQ : nso*nso;
    double * temp2 = (double*)malloc(max*sizeof(double));
    double * temp3 = (double*)malloc(nso*nso*sizeof(double));
    double * h     = (double*)malloc(nmo*nmo*sizeof(double));
    double ** hp   = H->pointer();

    F_DGEMM('n','t',nso,nmo,nso,1.0,&hp[0][0],nso,&Cap[0][0],nmo,0.0,temp2,nso);
    F_DGEMM('n','n',nmo,nmo,nso,1.0,&Cap[0][0],nmo,temp2,nso,0.0,h,nmo);

    // build Fock matrix:  sum_k (Q|kk)
    #pragma omp parallel for schedule (static)
    for (long int q = 0; q < nQ; q++) {
        double dum = 0.0;
        for (long int k = 0; k < ndocc; k++) {
            dum += Qso[q*nmo*nmo + k*nmo + k];
        }
        temp2[q] = 2.0 * dum;
    }

    // use temp and Qso as storage for Qmo( q, r, k)
    #pragma omp parallel for schedule (static)
    for (long int q = 0; q < nQ; q++) {
        for (long int r = 0; r < nmo; r++) {
            for (long int k = 0; k < ndocc; k++) {
                tmp[q*nmo*ndocc+k*nmo+r]  = Qso[q*nmo*nmo+k*nmo+r];
            }
        }
    }
    // I(r,s) = sum q k (q|ks)(q|rk)
    F_DGEMM('n','t',nmo,nmo,nQ*ndocc,1.0,tmp,nmo,tmp,nmo,0.0,temp3,nmo);

    // Fock matrix
    #pragma omp parallel for schedule (static)
    for (long int i = 0; i < nmo; i++) {
        for (long int j = 0; j < nmo; j++) {
            double dum = h[i*nmo+j] - temp3[i*nmo+j];
            dum += F_DDOT(nQ,temp2,1,Qso + i*nmo + j , nmo*nmo);
            F[i*nmo+j] = dum;
        }
    }

    free(h);
    free(tmp);
    free(temp2);
    free(temp3);
}

}} // end of namespaces

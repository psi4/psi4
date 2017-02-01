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

/**
  * Frozen natural orbitals
  * Eugene DePrince
  * April 2013
  *
  */

#include "psi4/psi4-dec.h"
#include "psi4/psifiles.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libtrans/mospace.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libmints/integral.h"
#include "psi4/libiwl/iwl.hpp"
#include "ccsd.h"
#include "blas.h"
#include "frozen_natural_orbitals.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/lib3index/dftensor.h"
#include "psi4/lib3index/cholesky.h"
#include "psi4/libmints/sieve.h"
#include "psi4/libqt/qt.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libdpd/dpd.h"
#define ID(x) ints->DPD_ID(x)



namespace psi{namespace fnocc{

FrozenNO::FrozenNO(SharedWavefunction wfn, Options& options):
  Wavefunction(options)
{
    // copy wave function.
    shallow_copy(wfn);
    reference_wavefunction_ = wfn;
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

    // quit if not RHF
    if ( options_.get_str("REFERENCE")!="RHF") {
        throw PsiException("FNOs only implemented for reference=rhf",__FILE__,__LINE__);
    }

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


    outfile->Printf("\n\n");
    outfile->Printf( "        *******************************************************\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *               Frozen Natural Orbitals               *\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *******************************************************\n");
    outfile->Printf("\n\n");


    outfile->Printf("        ==> Transform (OV|OV) integrals <==\n");
    outfile->Printf("\n");

    std::vector<std::shared_ptr<MOSpace> > spaces;
    spaces.push_back(MOSpace::occ);
    spaces.push_back(MOSpace::vir);
    std::shared_ptr<Wavefunction> wfn = reference_wavefunction_;
    std::shared_ptr<IntegralTransform> ints(new IntegralTransform(wfn, spaces, IntegralTransform::Restricted,
               IntegralTransform::DPDOnly, IntegralTransform::QTOrder, IntegralTransform::OccAndVir, false));
    ints->set_dpd_id(0);
    ints->set_keep_iwl_so_ints(true);
    ints->set_keep_dpd_so_ints(true);
    ints->initialize();
    ints->transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::occ, MOSpace::vir);

    outfile->Printf("\n");
    outfile->Printf("        ==> Build MP2 amplitudes, OPDM, and NOs <==\n");
    outfile->Printf("\n");

    dpdbuf4 amps1,amps2;
    std::shared_ptr<PSIO> psio = _default_psio_lib_;
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

    std::shared_ptr<Vector> epsA = reference_wavefunction_->epsilon_a();
    std::shared_ptr<Vector> epsB = reference_wavefunction_->epsilon_b();
    for(int h = 0; h < nirrep_; ++h){
        for(int a = frzcpi_[h]; a < doccpi_[h] + soccpi_[h]; ++a) aOccEvals[aOccCount++] = epsA->get(h, a);
        for(int b = frzcpi_[h]; b < doccpi_[h]; ++b)              bOccEvals[bOccCount++] = epsB->get(h, b);
        for(int a = doccpi_[h] + soccpi_[h]; a < nmopi_[h]; ++a)  aVirEvals[aVirCount++] = epsA->get(h, a);
        for(int b = doccpi_[h]; b < nmopi_[h]; ++b)               bVirEvals[bVirCount++] = epsB->get(h, b);
    }

    global_dpd_->buf4_init(&amps1, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");
    global_dpd_->buf4_sort(&amps1, PSIF_LIBTRANS_DPD, prqs, ID("[O,O]"), ID("[V,V]"), "MO Ints <OO|VV>");

    // T(ijab) -> T(jiab)
    global_dpd_->buf4_sort(&amps1, PSIF_LIBTRANS_DPD, prqs, ID("[O,O]"), ID("[V,V]"), "Tijab <OO|VV>");
    global_dpd_->buf4_sort(&amps1, PSIF_LIBTRANS_DPD, rpqs, ID("[O,O]"), ID("[V,V]"), "Tjiab <OO|VV>");
    global_dpd_->buf4_close(&amps1);

    // only worry about alpha-beta
    global_dpd_->buf4_init(&amps1, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
    global_dpd_->buf4_init(&amps2, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),ID("[O,O]"), ID("[V,V]"), 0, "Tjiab <OO|VV>");

    double emp2_os = 0.0;
    double emp2_ss = 0.0;
    for(int h = 0; h < nirrep_; ++h){

        global_dpd_->buf4_mat_irrep_init(&amps1, h);
        global_dpd_->buf4_mat_irrep_rd(&amps1, h);

        global_dpd_->buf4_mat_irrep_init(&amps2, h);
        global_dpd_->buf4_mat_irrep_rd(&amps2, h);

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
        global_dpd_->buf4_mat_irrep_close(&amps1, h);
        global_dpd_->buf4_mat_irrep_close(&amps2, h);
    }
    global_dpd_->buf4_close(&amps1);

    double escf = Process::environment.globals["SCF TOTAL ENERGY"];
    Process::environment.globals["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = emp2_os;
    Process::environment.globals["MP2 SAME-SPIN CORRELATION ENERGY"] = emp2_ss;
    Process::environment.globals["MP2 CORRELATION ENERGY"] = emp2_os + emp2_ss;
    Process::environment.globals["MP2 TOTAL ENERGY"] = emp2_os + emp2_ss + escf;

    // build amps1(ij,ab) = 2*T(ij,ab) - T(ji,ab)
    global_dpd_->buf4_init(&amps1, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
    global_dpd_->buf4_scm(&amps1, 2.0);
    global_dpd_->buf4_axpy(&amps2, &amps1, -1.0);

    global_dpd_->buf4_close(&amps1);
    global_dpd_->buf4_close(&amps2);

    outfile->Printf("        OS MP2 correlation energy:       %20.12lf\n",emp2_os);
    outfile->Printf("        SS MP2 correlation energy:       %20.12lf\n",emp2_ss);
    outfile->Printf("        MP2 correlation energy:          %20.12lf\n",emp2_os+emp2_ss);
    outfile->Printf("      * MP2 total energy:                %20.12lf\n",emp2_os+emp2_ss+escf);
    outfile->Printf("\n");

    // scale amps by denominator
    global_dpd_->buf4_init(&amps2, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
    global_dpd_->buf4_init(&amps1, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),ID("[O,O]"), ID("[V,V]"), 0, "Tijab <OO|VV>");
    for(int h = 0; h < nirrep_; ++h){

        global_dpd_->buf4_mat_irrep_init(&amps1, h);
        global_dpd_->buf4_mat_irrep_rd(&amps1, h);

        global_dpd_->buf4_mat_irrep_init(&amps2, h);
        global_dpd_->buf4_mat_irrep_rd(&amps2, h);

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

        global_dpd_->buf4_mat_irrep_wrt(&amps1, h);
        global_dpd_->buf4_mat_irrep_close(&amps1, h);
        global_dpd_->buf4_mat_irrep_wrt(&amps2, h);
        global_dpd_->buf4_mat_irrep_close(&amps2, h);
    }
    global_dpd_->buf4_close(&amps1);
    global_dpd_->buf4_close(&amps2);

    // build virtual-virtual block of opdm: sum(ijc) 2.0 * [ 2 t(ij,ac) - t(ji,ac) ] * t(ij,bc)
    dpdfile2 Dab;
    global_dpd_->file2_init(&Dab, PSIF_LIBTRANS_DPD,    0, 1, 1, "Dab");
    global_dpd_->buf4_init(&amps2, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
    global_dpd_->buf4_init(&amps1, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),ID("[O,O]"), ID("[V,V]"), 0, "Tijab <OO|VV>");
    global_dpd_->contract442(&amps1, &amps2, &Dab, 3, 3, 2.0, 0.0);
    global_dpd_->buf4_close(&amps1);
    global_dpd_->buf4_close(&amps2);
    global_dpd_->file2_close(&Dab);

    // diagonalize virtual-virtual block of opdm
    int symmetry = Ca_->symmetry();
    std::shared_ptr<Matrix> D (new Matrix("Dab",nirrep_,aVirOrbsPI,aVirOrbsPI,symmetry));

    global_dpd_->file2_init(&Dab, PSIF_LIBTRANS_DPD,    0, 1, 1, "Dab");
    global_dpd_->file2_mat_init(&Dab);
    global_dpd_->file2_mat_rd(&Dab);
    for (int h = 0; h < nirrep_; h++) {
        int v = Dab.params->rowtot[h];

        double ** mat = D->pointer(h);
        for (int a = 0; a < v; a++) {
            for (int b = 0; b < v; b++) {
                mat[a][b] = Dab.matrix[h][a][b];
            }
        }
    }
    global_dpd_->file2_close(&Dab);

    // done with dpd and ints ... reset
    psio->close(PSIF_LIBTRANS_DPD, 1);
    ints.reset();

    std::shared_ptr<Matrix> eigvec (new Matrix("Dab eigenvectors",nirrep_,aVirOrbsPI,aVirOrbsPI,symmetry));
    std::shared_ptr<Vector> eigval (new Vector("Dab eigenvalues",nirrep_,aVirOrbsPI));
    D->diagonalize(eigvec,eigval,descending);

    // overwrite ao/mo C matrix with ao/no transformation
    std::shared_ptr<Matrix> temp (new Matrix("temp",nirrep_,nsopi_,aVirOrbsPI,symmetry));
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

    if (!options_["ACTIVE_NAT_ORBS"].has_changed()) {

        if ( !options_["OCC_PERCENTAGE"].has_changed() ) {
            // use occupancy tolerance:
            for (int h = 0; h < nirrep_; h++) {
                newVirOrbsPI[h] = 0;
                double * vec = eigval->pointer(h);
                for (int a = 0; a < aVirOrbsPI[h]; a++) {
                    if ( vec[a] > cutoff ) newVirOrbsPI[h]++;
                }
            }
            outfile->Printf("        Cutoff for significant NO occupancy: %5.3le\n",cutoff);
            outfile->Printf("\n");
        }else {
            // retain orbitals with some percentage of total virtual occupation

            // calculate total trace of vv OPDM
            double occ_total = 0.0;
            for (int h = 0; h < nirrep_; h++) {
                newVirOrbsPI[h] = 0;
                double * vec = eigval->pointer(h);
                for (int a = 0; a < aVirOrbsPI[h]; a++) {
                    occ_total += vec[a];
                }
            }

            // choose virtuals from greatest occupation to least until
            // ratio of occupations matches occ_fraction
            double frac = options_.get_double("OCC_PERCENTAGE")/100.0;

            int * skip = (int*)malloc(nvirt*sizeof(int));
            memset((void*)skip,'\0',nvirt*sizeof(int));

            double my_occ_total = 0.0;
            for (int h = 0; h < nirrep_; h++) {
                for (int a = 0; a < aVirOrbsPI[h]; a++) {

                    int maxcount = -999;
                    int maxh     = -999;
                    int maxb     = -999;
                    double max   = -9e99;

                    int count = 0;
                    for (int h2 = 0; h2 < nirrep_; h2++) {
                        double * vec = eigval->pointer(h2);
                        for (int b = 0; b < aVirOrbsPI[h2]; b++) {
                            if ( !skip[count] && vec[b] > max ) {
                                max      = vec[b];
                                maxh     = h2;
                                maxb     = b;
                                maxcount = count;
                            }
                            count++;
                        }
                    }
                    if ( maxcount < 0 ) {
                        throw PsiException("having trouble sorting virtual NOs",__FILE__,__LINE__);
                    }
                    skip[maxcount] = 1;

                    double my_occ = eigval->pointer(maxh)[maxb];
                    if ( ( my_occ_total / occ_total ) < frac ) {
                        my_occ_total += my_occ;
                        newVirOrbsPI[maxh]++;
                    }
                }
            }
            free(skip);
            outfile->Printf("        Cutoff for retaining NOs is %5.2lf%% occupancy\n",frac*100.0);
            outfile->Printf("\n");
        }
    }else{
        // use user-specified number of virtuals
        for (int h = 0; h < nirrep_; h++) {
            newVirOrbsPI[h] = (long int)options_["ACTIVE_NAT_ORBS"][h].to_double();
        }
    }

    outfile->Printf("        No. virtuals per irrep (original):  [");
    for (int h = 0; h < nirrep_; h++) outfile->Printf("%4i",aVirOrbsPI[h]);
    outfile->Printf(" ]\n");
    outfile->Printf("        No. virtuals per irrep (truncated): [");
    for (int h = 0; h < nirrep_; h++) outfile->Printf("%4i",newVirOrbsPI[h]);
    outfile->Printf(" ]\n");
    outfile->Printf("\n");

    int nvirt_no = 0;
    int nvirt    = 0;
    for (int h = 0; h < nirrep_; h++) {
        nvirt    += aVirOrbsPI[h];
        nvirt_no += newVirOrbsPI[h];
    }
    outfile->Printf("        Retaining %i of %i virtual orbitals.\n",nvirt_no,nvirt);
    outfile->Printf("\n");

    // transform Fock matrix to truncated NO basis

    std::shared_ptr<Matrix> Fab (new Matrix("Fab(NO)",nirrep_,newVirOrbsPI,newVirOrbsPI,symmetry));
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
    std::shared_ptr<Matrix> eigvecF (new Matrix("Fab eigenvectors",nirrep_,newVirOrbsPI,newVirOrbsPI,symmetry));
    std::shared_ptr<Vector> eigvalF (new Vector("Fab eigenvalues",nirrep_,newVirOrbsPI));
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
    std::shared_ptr<Vector> eps = epsilon_a_;
    for (int h = 0; h < nirrep_; h++) {
        double * epsp = eps->pointer(h);
        double * eigp = eigvalF->pointer(h);
        for (int a = 0; a < newVirOrbsPI[h]; a++) {
            epsp[doccpi_[h]+a] = eigp[a];
        }
    }

    // free memory
    delete[] newVirOrbsPI;

    delete[] aOccOrbsPI;
    delete[] bOccOrbsPI;
    delete[] aVirOrbsPI;
    delete[] bVirOrbsPI;
    delete[] aOccEvals;
    delete[] bOccEvals;
    delete[] aVirEvals;
    delete[] bVirEvals;

    tstop();
}

// DF FNO class members

DFFrozenNO::DFFrozenNO(std::shared_ptr<Wavefunction>wfn,Options&options):
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

  outfile->Printf("  ==> 3-index integrals <==\n");
  outfile->Printf("\n");

  long int o = ndoccact;
  long int v = nvirt;
  long int nQ;

  // 1.  read scf 3-index integrals from disk

  // get ntri from sieve
  std::shared_ptr<ERISieve> sieve (new ERISieve(basisset_, options_.get_double("INTS_TOLERANCE")));
  const std::vector<std::pair<int, int> >& function_pairs = sieve->function_pairs();
  long int ntri = function_pairs.size();

  // read integrals that were written to disk in the scf
  long int nQ_scf = Process::environment.globals["NAUX (SCF)"];
  if ( options_.get_str("SCF_TYPE") == "DF" ) {

      std::shared_ptr<BasisSet> primary = get_basisset("ORBITAL");
      std::shared_ptr<BasisSet> auxiliary = get_basisset("DF_BASIS_SCF");

      nQ_scf = auxiliary->nbf();
      Process::environment.globals["NAUX (SCF)"] = nQ_scf;
  }

  std::shared_ptr<Matrix> Qmn = SharedMatrix(new Matrix("Qmn Integrals",nQ_scf,ntri));
  double** Qmnp = Qmn->pointer();
  std::shared_ptr<PSIO> psio(new PSIO());
  psio->open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DFSCF_BJ, "(Q|mn) Integrals", (char*) Qmnp[0], sizeof(double) * ntri * nQ_scf);
  psio->close(PSIF_DFSCF_BJ,1);

  // unpack and write again in my format
  std::shared_ptr<Matrix>L (new Matrix("3-index ERIs (SCF)", nQ_scf , nso*nso));
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

      std::shared_ptr<BasisSet> auxiliary = get_basisset("DF_BASIS_CC");
      std::shared_ptr<DFTensor> DF (new DFTensor(basisset(),auxiliary,Ca(),ndocc,nvirt+nfzv,ndoccact,nvirt,options_));
      nQ = auxiliary->nbf();
      std::shared_ptr<Matrix> tmp = DF->Qso();
      double ** Qso = tmp->pointer();

      // write Qso to disk
      psio->open(PSIF_DCC_QSO,PSIO_OPEN_OLD);
      psio->write_entry(PSIF_DCC_QSO,"Qso CC",(char*)&Qso[0][0],nQ*nso*nso*sizeof(double));
      psio->close(PSIF_DCC_QSO,1);
      outfile->Printf("    Number of auxiliary functions:       %5li\n",nQ);

      // stick nQ in process environment so ccsd can know it
      Process::environment.globals["NAUX (CC)"] = (double)nQ;
  }else{

      // Cholesky

      // read integrals from disk if they were generated in the SCF
      if ( options_.get_str("SCF_TYPE") == "CD" ) {
          outfile->Printf("        Reading Cholesky vectors from disk ...\n");
          nQ = Process::environment.globals["NAUX (SCF)"];
          outfile->Printf("        Cholesky decomposition threshold: %8.2le\n", options_.get_double("CHOLESKY_TOLERANCE"));
          outfile->Printf("        Number of Cholesky vectors:          %5li\n",nQ);

          // ntri comes from sieve above
          std::shared_ptr<Matrix> Qmn = SharedMatrix(new Matrix("Qmn Integrals",nQ,ntri));
          double** Qmnp = Qmn->pointer();
          // TODO: use my 3-index integral file in SCF for DFCC jobs
          psio->open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
          psio->read_entry(PSIF_DFSCF_BJ, "(Q|mn) Integrals", (char*) Qmnp[0], sizeof(double) * ntri * nQ);
          psio->close(PSIF_DFSCF_BJ,1);

          std::shared_ptr<Matrix>L (new Matrix("CD Integrals", nQ , nso*nso));
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
          outfile->Printf("        Generating Cholesky vectors ...\n");
          std::shared_ptr<BasisSet> primary = basisset();
          std::shared_ptr<IntegralFactory> integral (new IntegralFactory(primary,primary,primary,primary));
          double tol = options_.get_double("CHOLESKY_TOLERANCE");
          std::shared_ptr<CholeskyERI> Ch (new CholeskyERI(std::shared_ptr<TwoBodyAOInt>(integral->eri()),0.0,tol,Process::environment.get_memory()));
          Ch->choleskify();
          nQ  = Ch->Q();
          std::shared_ptr<Matrix> L = Ch->L();
          double ** Lp = L->pointer();

          // write Qso to disk
          psio->open(PSIF_DCC_QSO,PSIO_OPEN_OLD);
          psio->write_entry(PSIF_DCC_QSO,"Qso CC",(char*)&Lp[0][0],nQ*nso*nso*sizeof(double));
          psio->close(PSIF_DCC_QSO,1);
          outfile->Printf("        Cholesky decomposition threshold: %8.2le\n", tol);
          outfile->Printf("        Number of Cholesky vectors:          %5li\n",nQ);

      }

      // stick nQ in process environment so ccsd can know it
      Process::environment.globals["NAUX (CC)"] = (double)nQ;
  }
  outfile->Printf("\n");
}

/*
    build 4-index eri's from 3-index integrals
*/
void DFFrozenNO::FourIndexIntegrals() {

    outfile->Printf("  ==> Build 4-index ERI's from 3-index integrals <==\n");
    outfile->Printf("\n");

    long int o  = ndoccact;
    long int v  = nvirt;
    long int nQ = Process::environment.globals["NAUX (CC)"];

    double ** Cap = Ca()->pointer();

    // transform 3-index integrals to MO basis

    psio_address addr1 = PSIO_ZERO;
    psio_address addr2 = PSIO_ZERO;
    double * buf1 = (double*)malloc(nso*nso*sizeof(double));
    double * buf2 = (double*)malloc(nso*nso*sizeof(double));

    std::shared_ptr<PSIO> psio(new PSIO());
    psio->open(PSIF_DCC_QSO,PSIO_OPEN_OLD);
    for (int q = 0; q < nQ; q++) {
        psio->read(PSIF_DCC_QSO,"Qso CC",(char*)&buf1[0],nso*nso*sizeof(double),addr1,&addr1);
        F_DGEMM('n','n',nmo,nso,nso,1.0,&Cap[0][0],nmo,buf1,nso,0.0,buf2,nmo);
        F_DGEMM('n','t',nmo,nmo,nso,1.0,&Cap[0][0],nmo,buf2,nmo,0.0,buf1,nmo);
        for (int p = 0; p < nmo; p++) {
            for (int q = p; q < nmo; q++) {
                buf2[Position(p,q)] = buf1[p*nmo+q];
            }
        }
        psio->write(PSIF_DCC_QSO,"Qmo CC",(char*)&buf2[0],nmo*(nmo+1)/2*sizeof(double),addr2,&addr2);
    }
    free(buf2);
    free(buf1);

    // hopefully nQ*nmo*(nmo+1)/2 will fit in memory
    long int memory = Process::environment.get_memory();
    if ( memory < nmo*(nmo+1)/2*nQ*sizeof(double) ) {
        throw PsiException("Not enough memory (FourIndexIntegrals)",__FILE__,__LINE__);
    }
    double * Qmo = (double*)malloc(nmo*(nmo+1)/2*nQ*sizeof(double));

    psio->read_entry(PSIF_DCC_QSO,"Qmo CC",(char*)&Qmo[0],nmo*(nmo+1)/2*nQ*sizeof(double));
    psio->close(PSIF_DCC_QSO,1);

    IWL * iwl = new IWL(psio.get(), PSIF_MO_TEI, 1.0e-16, 0, 0);
    for (int p = nfzc; p < nmo; p++) {
        for (int q = p; q < nmo; q++) {
            int pq = Position(p,q);
            for (int r = nfzc; r < nmo; r++) {
                for (int s = r; s < nmo; s++) {
                    int rs = Position(r,s);
                    if ( rs > pq ) continue;
                    double val = C_DDOT(nQ,Qmo+pq,nmo*(nmo+1)/2,Qmo+rs,nmo*(nmo+1)/2);
                    iwl->write_value(p, q, r, s, val, false, "outfile", 0);
                }
            }
        }
    }
    iwl->flush(1);
    iwl->set_keep_flag(1);
    delete iwl;

    free(Qmo);
}

/*
 * build natural orbitals and transform TEIs
 */
void DFFrozenNO::ComputeNaturalOrbitals(){


  outfile->Printf( "  ==> Frozen Natural Orbitals <==\n");
  outfile->Printf("\n");


  long int o      = ndoccact;
  long int v      = nvirt;
  long int nQ     = Process::environment.globals["NAUX (CC)"];
  long int nQ_scf = Process::environment.globals["NAUX (SCF)"];
  long int memory = Process::environment.get_memory();

  if ( memory < 8L*(3L*nso*nso+nso*nso*nQ+o*v*nQ) ) {
      throw PsiException("not enough memory (fno)",__FILE__,__LINE__);
  }

  std::shared_ptr<PSIO> psio(new PSIO());

  // read in 3-index integrals specific to the CC method:
  double * tmp2 = (double*)malloc(nso*nso*nQ * sizeof(double));
  psio->open(PSIF_DCC_QSO,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_QSO,"Qso CC",(char*)&tmp2[0],nQ*nso*nso*sizeof(double));
  psio->close(PSIF_DCC_QSO,1);

  // transform Qso -> Qov:
  TransformQ(nQ,tmp2);
  double * Qov = (double*)malloc(o*v*nQ*sizeof(double));
  C_DCOPY(o*v*nQ,tmp2,1,Qov,1);
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

  std::shared_ptr<Vector> eps_test = epsilon_a();
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

  outfile->Printf("        Doubles contribution to MP2 energy in full space: %20.12lf\n",emp2);
  outfile->Printf("\n");

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
      C_DCOPY(v,Dab+(v-1-i)*v,1,temp+i*v,1);
  }

  // establish cutoff for frozen virtuals
  double cutoff = options_.get_double("OCC_TOLERANCE");
  nvirt_no = 0;

  if (!options_["ACTIVE_NAT_ORBS"].has_changed()) {

      if ( !options_["OCC_PERCENTAGE"].has_changed() ) {

          // use occupancy tolerance:
          for (long int i=0; i<v; i++) if (eigvalDab[i]>cutoff) nvirt_no++;
          outfile->Printf("        Cutoff for significant NO occupancy: %5.3le\n",cutoff);
          outfile->Printf("\n");
      }else {

          // retain orbitals with some percentage of total virtual occupation

          // calculate total trace of vv OPDM
          double occ_total = 0.0;
          nvirt_no = 0.0;
          for (int a = 0; a < v; a++) {
              occ_total += eigvalDab[a];
          }

          // choose virtuals from greatest occupation to least until
          // ratio of occupations matches occ_fraction
          double frac = options_.get_double("OCC_PERCENTAGE")/100.0;

          int * skip = (int*)malloc(v*sizeof(int));
          memset((void*)skip,'\0',v*sizeof(int));

          double my_occ_total = 0.0;
          for (int a = 0; a < v; a++) {

              int maxcount = -999;
              int maxb     = -999;
              double max   = -9e99;

              int count = 0;
              for (int b = 0; b < v; b++) {
                  if ( !skip[count] && eigvalDab[b] > max ) {
                      max      = eigvalDab[b];
                      maxb     = b;
                      maxcount = count;
                  }
                  count++;
              }
              if ( maxcount < 0 ) {
                  throw PsiException("having trouble sorting virtual NOs",__FILE__,__LINE__);
              }
              skip[maxcount] = 1;

              double my_occ = eigvalDab[maxb];
              if ( ( my_occ_total / occ_total ) < frac ) {
                  my_occ_total += my_occ;
                  nvirt_no++;
              }
          }
          free(skip);
          outfile->Printf("        Cutoff for retaining NOs is %5.2lf%% occupancy\n",frac*100.0);
          outfile->Printf("\n");
      }
  }else{
      // use user-specified number of virtuals
      nvirt_no = (long int)options_["ACTIVE_NAT_ORBS"][0].to_double();
  }

  outfile->Printf("        Number of virtual orbitals in original space:  %5li\n",v);
  outfile->Printf("        Number of virtual orbitals in truncated space: %5li\n",nvirt_no);
  outfile->Printf("\n");

  // transform Fock matrix to MP2 NO basis
  memset((void*)newFock,'\0',v*v*sizeof(double));
  C_DCOPY(v,F+ndoccact,1,newFock,v+1);
  F_DGEMM('n','n',v,nvirt_no,v,1.0,newFock,v,temp,v,0.0,Dab,v);
  F_DGEMM('t','n',nvirt_no,nvirt_no,v,1.0,temp,v,Dab,v,0.0,newFock,nvirt_no);

  // diagonalize new Fock matrix for semi-canonical orbitals
  Diagonalize(nvirt_no,newFock,neweps);

  // construct full mo -> no transformation matrix
  F_DGEMM('n','n',v,nvirt_no,nvirt_no,1.0,temp,v,newFock,nvirt_no,0.0,Dab,v);

  // put orbital energies back in F - doesn't matter in this implementation
  C_DCOPY(nvirt_no,neweps,1,F+ndoccact,1);

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

  std::shared_ptr<psi::Wavefunction> ref = reference_wavefunction_;

  std::shared_ptr<Matrix> Caomo = ref->Ca();

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

  std::shared_ptr<psi::Wavefunction> ref = reference_wavefunction_;

  std::shared_ptr<Matrix> Caomo = ref->Ca();

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
    C_DCOPY(nso*nso*nQ,Qso,1,tmp,1);
    F_DGEMM('n','n',nmo,nso*nQ,nso,1.0,&Cap[0][0],nmo,tmp,nso,0.0,Qso,nmo);
    #pragma omp parallel for schedule (static)
    for (long int q = 0; q < nQ; q++) {
        for (long int mu = 0; mu < nso; mu++) {
            C_DCOPY(nmo,Qso+q*nso*nmo+mu*nmo,1,tmp+q*nso*nmo+mu,nmo);
        }
    }
    F_DGEMM('n','n',nmo,nmo*nQ,nso,1.0,&Cap[0][0],nmo,tmp,nso,0.0,Qso,nmo);

    // build Fock matrix:

    // transform H
    // one-electron integrals
    std::shared_ptr<MintsHelper> mints(new MintsHelper(basisset_, options_, 0));
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
            dum += C_DDOT(nQ,temp2,1,Qso + i*nmo + j , nmo*nmo);
            F[i*nmo+j] = dum;
        }
    }

    free(h);
    free(tmp);
    free(temp2);
    free(temp3);
}

}} // end of namespaces

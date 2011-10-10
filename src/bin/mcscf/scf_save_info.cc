#include <utility>
#include <algorithm>
#include <cstdio>

#include <cmath>

#include <libchkpt/chkpt.hpp>
#include <libmoinfo/libmoinfo.h>
#include <liboptions/liboptions.h>
#include <libmints/molecule.h>
#include <libmints/matrix.h>
#include <libmints/factory.h>
#include <psifiles.h>

#include <psi4-dec.h>

#include "scf.h"

namespace psi{ namespace mcscf{

using namespace std;

void SCF::save_info()
{
    // No projection of MOs, just yet
    nmo_ = nso_;

    // figure out how many frozen orbitals per irrep
    int nfrzc = Process::environment.molecule()->nfrozen_core();
    intvec frz;
    for(int h = 0; h < nirreps; ++h) frz.push_back(0);
    vector<std::pair<double, int> > sorted_evals;
    for(int h = 0; h < nirreps; ++h)
      for(int i = 0; i < sopi[h]; ++i)
        sorted_evals.push_back( make_pair(epsilon->get(h,i),h) );
    sort(sorted_evals.begin(),sorted_evals.end());
    for(int i = 0; i < nfrzc; ++i)
      frz[sorted_evals[i].second]++;

    for(int h = 0; h < nirreps; ++h){
        doccpi_[h] = docc[h];
        soccpi_[h] = actv[h];
        nmopi_[h]  = nsopi_[h];
        frzcpi_[h] = frz[h];
        frzvpi_[h] = 0;
    }
    // Save the eigenvectors after rotating them
    if(options_.get_int("ROTATE_MO_ANGLE") != 0){
        int mo_rotate_angle = options_.get_int("ROTATE_MO_ANGLE");
        int p = options_.get_int("ROTATE_MO_P") -1;  // P, Q and IRREPS are one-based
        int q = options_.get_int("ROTATE_MO_Q") -1;
        int h = options_.get_int("ROTATE_MO_IRREP") - 1;

        fprintf(outfile,"\n\n  Rotating MOs %d and %d of irrep %d by %d degrees",
                        p,q,h,mo_rotate_angle);
        double angle = static_cast<double>(mo_rotate_angle) * acos(-1.0) / 180.0;
        for(int i = 0; i < sopi[h]; ++i){
            double Cp = cos(angle) * C->get(h,i,p) + sin(angle) * C->get(h,i,q);
            double Cq = cos(angle) * C->get(h,i,q) - sin(angle) * C->get(h,i,p);
            C->set(h,i,p,Cp);
            C->set(h,i,q,Cq);
        }
    }

    // Store the information in wavefunction's objects, for later access
    Ca_ = SharedMatrix(factory_->create_matrix("C"));
    Cb_ = Ca_;
    epsilon_a_ = SharedVector(factory_->create_vector());
    epsilon_b_ = epsilon_a_;
    for(int h = 0; h < nirreps; ++h){
        for(int so = 0; so < nsopi_[h]; ++so){
            epsilon_a_->set(h, so, epsilon->get(h, so));
            for(int mo = 0; mo < nmopi_[h]; ++mo){
                Ca_->set(h, so, mo, C->get(h, so, mo));
            }
        }
    }

    Process::environment.globals["MCSCF ENERGY"] = total_energy;
    Process::environment.globals["CURRENT ENERGY"] = total_energy;
    energy_ = total_energy;
    psio_->close(PSIF_CHKPT, 1);
    cleanup();

    return;

//    // Writes out the total number of irreps in the point group in which the molecule is being considered which have non-zero number of basis functions.
//      int n_so_typs = 0;
//      for(int h = 0; h < nirreps; ++h){
//        if( docc[h] + actv[h] > 0 ) n_so_typs++;
//      }
//      chkpt_->wt_nsymhf(n_so_typs);

//      // Write out the total number of molecular orbitals.
//      chkpt_->wt_nmo(nso); // TODO: find nmo

//      // Write out the dimensionality of ALPHA and BETA vectors of two-electron coupling coefficients for open shells.
//      int tmp_iopen = ioff[moinfo_scf->get_nactv()];
//      if(reference == tcscf) tmp_iopen = -tmp_iopen;
//      chkpt_->wt_iopen(tmp_iopen);

//      // Write open-shell coupling coefficients
//      if(moinfo_scf->get_nactv() > 0){
//        double** ccvecs;
//        allocate2(double,ccvecs,2,ioff[moinfo_scf->get_nactv()]);
//        for(int i=0; i < ioff[reference == tcscf]; i++) {
//          ccvecs[0][i] = 0.0;
//          ccvecs[1][i] = 0.0;
//        }
//        chkpt_->wt_ccvecs(ccvecs);
//        release2(ccvecs);
//      }

//      // Writes out the total energy.
//      chkpt_->wt_etot(total_energy);
//      chkpt_->wt_escf(total_energy);
//      chkpt_->wt_eref(total_energy);

//      // Write the orbitals per irreps arrays
//      chkpt_->wt_orbspi(&sopi[0]);
//      chkpt_->wt_clsdpi(&docc[0]);
//      chkpt_->wt_openpi(&actv[0]);


//      // Read the number of frozen MOs
//      int nfrzc = chkpt_->rd_nfzc();

//      int* frz = new int[nirreps];
//      for(int h = 0; h < nirreps; ++h) frz[h] = 0;
//      chkpt_->wt_frzvpi(frz);

//      vector<std::pair<double, int> > sorted_evals;

//      for(int h = 0; h < nirreps; ++h)
//        for(int i = 0; i < sopi[h]; ++i)
//          sorted_evals.push_back( make_pair(epsilon->get(h,i),h) );

//      // Sort the eigenvalues by energy
//      sort(sorted_evals.begin(),sorted_evals.end());

//      for(int i = 0; i < nfrzc; ++i){
//        frz[sorted_evals[i].second]++;
//      }

//      chkpt_->wt_frzcpi(frz);

//      delete[] frz;

//      // Save the eigenvectors after rotating them
//      if(options_.get_int("ROTATE_MO_ANGLE") != 0){
//        int mo_rotate_angle = options_.get_int("ROTATE_MO_ANGLE");
//        int p = options_.get_int("ROTATE_MO_P") -1;  // P, Q and IRREPS are one-based
//        int q = options_.get_int("ROTATE_MO_Q") -1;
//        int h = options_.get_int("ROTATE_MO_IRREP") - 1;

//        fprintf(outfile,"\n\n  Rotating MOs %d and %d of irrep %d by %d degrees",
//                        p,q,h,mo_rotate_angle);
//        double angle = static_cast<double>(mo_rotate_angle) * acos(-1.0) / 180.0;
//        for(int i = 0; i < sopi[h]; ++i){
//          double Cp = cos(angle) * C->get(h,i,p) + sin(angle) * C->get(h,i,q);
//          double Cq = cos(angle) * C->get(h,i,q) - sin(angle) * C->get(h,i,p);
//          C->set(h,i,p,Cp);
//          C->set(h,i,q,Cq);
//        }
//      }

//      double** C_save;
//      allocate2(double,C_save,nso,nso);

//      for(int h = 0; h < nirreps; ++h)
//        for(int i = 0; i < sopi[h]; ++i)
//          for(int j = 0; j < sopi[h]; ++j)
//            C_save[i + block_offset[h]][j + block_offset[h]] = C->get(h,i,j);
//      chkpt_->wt_scf(C_save);

//      release2(C_save);

//      int k = 0;
//      double* evals = new double[nso];
//      for(int h = 0; h < nirreps; ++h){
//        for(int i = 0; i < sopi[h]; ++i){
//          evals[k] = epsilon->get(h,i);
//          k++;
//        }
//      }
//      chkpt_->wt_evals(evals);
//      delete[] evals;
}

}} /* End Namespaces */

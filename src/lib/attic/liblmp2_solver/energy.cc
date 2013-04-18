/*! \defgroup LMP2 lmp2: LMP2 Evaluation of Energy */

/*!
 ** \file
 ** \ingroup LMP2
 ** \LMP2 evaluation of energy
 */

#include <liblmp2_solver/lmp2.h>
#include <libqt/qt.h>

namespace boost {
template<class T> class shared_ptr;
}

namespace psi{

class BasisSet;
class Options;
class PSIO;
class Chkpt;

namespace lmp2 {

#ifdef HAVE_MADNESS

madness::Future<double> LMP2::energy(const SharedMatrix T2, const int &ij, const int &iter) {

    // *** Compute the new MP2 energy ***

    int i = ij_i_map_[ij];
    int j = ij_j_map_[ij];

    double lmp2_energy = 0.0;

    if (i != j) {
        for(int a=0; a < pair_domain_len_[ij]; a++) {
            for(int b=0; b < pair_domain_len_[ij]; b++) {
                double Kab = eri_ij_[ij]->get(0, a, b);
                if(fabs(Kab) > 1e-14) {
                    double Tab = T2->get(0,a,b);
                    double Tba = T2->get(0,b,a);
                    //                    E_OS_ += 2 * Kab * ( Tab );
                    //                    E_SS_ += 2 * Kab * ( Tab - Tba );
                    lmp2_energy += 2 * Kab * (2 * Tab - Tba);
                }
            }
        }
    }
    else {
        for(int a=0; a < pair_domain_len_[ij]; a++) {
            for(int b=0; b < pair_domain_len_[ij]; b++) {
                double Kab = eri_ij_[ij]->get(0, a, b);
                if(fabs(Kab) > 1e-14) {
                    double Tab = T2->get(0,a,b);
                    double Tba = T2->get(0,b,a);

                    //                    E_OS_ += 2 * Kab * Tab ;
                    //                    E_SS_ += 2 * Kab * ( Tab - Tba );
                    lmp2_energy += Kab * (2 * Tab - Tba);
                }
            }
        }
    }

    return madness::Future<double>(lmp2_energy);
}

madness::Void LMP2::Elmp2_local_sum(const double &elmp2_local) {
    mutex_->lock();
    Elmp2_ += elmp2_local;
    mutex_->unlock();
    return madness::None;
}

madness::Void LMP2::Drms_local_sum(const double &drms_local) {
    mutex_->lock();
    Drms_T2_ += drms_local;
    mutex_->unlock();
    return madness::None;
}


void LMP2::print_results(const int &iter) const {

    if (me_ == 0) {
        if(iter != 0) {
            if(iter > diis_start_ & diis_) {
                fprintf(outfile, "  %3d %20.12f %20.12f %20.12f  DIIS\n", iter, Elmp2_, Delta_Elmp2_, Drms_T2_);
            }
            else  {
                fprintf(outfile, "  %3d %20.12f %20.12f %20.12f\n", iter, Elmp2_, Delta_Elmp2_, Drms_T2_);
            }
        }
        else {
            fprintf(outfile, "      \t     Correlation\t     Delta    \t\t   RMS\n");
            fprintf(outfile, "  Iter\t       Energy\t\t (Corr. Energy)\t     (T2 Amplitudes)\n");
            fprintf(outfile, "  ------------------------------------------------------------------------\n");

            fprintf(outfile, "  %3d %20.12f\n", iter, Elmp2_);
        }
    }

}

#endif // have_madness

}} // End of psi and lmp2 namespaces

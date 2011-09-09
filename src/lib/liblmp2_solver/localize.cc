/*! \defgroup LMP2 lmp2: LMP2 Evaluation of Energy */

/*!
 ** \file
 ** \ingroup LMP2
 ** \LMP2 evaluation of energy
 */

#include <liblmp2_solver/lmp2.h>
#include <libqt/qt.h>
#include <libciomr/libciomr.h>

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

void LMP2::localize_pipek_mezey() {


    // Get the orbital eigenvalues from the wavefunction object
    SharedVector evals = wfn_->epsilon_a();

    timer_on("Compute Density Matrix");
    C_DGEMM('n', 't', nso_, nso_, ndocc_, 1, C_->get_pointer(), nso_, C_->get_pointer(), nso_, 0, D_AO_->get_pointer(), nso_);
    timer_off("Compute Density Matrix");

    timer_on("Localization");
    std::vector< boost::shared_ptr<BasisSet> > ao_bas;
    for (int i=0; i < natom_; i++) {
        ao_bas.push_back(basisset_->atomic_basis_set(i));
    }


    for (int i=0; i < natom_; i++) {
        for (int j=0; j < nshell_; j++) {
            if (i == basisset_->shell(j)->ncenter()) {
                ao_start_.push_back(basisset_->shell(j)->function_index());
                ao_stop_.push_back(ao_start_[i]);
                for (int k=0; k < ao_bas[i]->nshell(); k++) {
                    ao_stop_[i] += ao_bas[i]->shell(k)->nfunction();
                }
                break;
            }
        }
    }
    ao_bas.clear();

    if (me_ == 0) {
        fprintf(outfile, "\n  ====> Orbital Localization <====\n\n");
        fprintf(outfile, "  Iter    Max Rotation    Convergence\n");
        fprintf(outfile, "  -----------------------------------\n");
    }

    SharedMatrix V(occ_occ_->create_matrix("V Matrix"));
    SharedMatrix U(occ_occ_->create_matrix("U Matrix"));
    SharedMatrix VV(occ_occ_->create_matrix("VV Matrix"));
    SharedMatrix F_oo(occ_occ_->create_matrix("F_occ"));


    V->identity();

    double **LCtmp;
    double alphalast, conv;

    timer_on("localization iterations");
    for(int iter=0; iter < 100; iter++) {

        double Ast,Bst;
        double AB, cos4a;
        double alpha, alphamax = 0.0;
        double Pst, Pss, Ptt;
        // Compute 2x2 rotations for Pipek-Mezey localization
        for (int s=nfocc_; s < ndocc_; s++) {
            for (int t=nfocc_; t < s; t++) {

                Ast = 0.0;
                Bst = 0.0;

                for (int A=0; A < natom_; A++) {
                    Pst = Pss = Ptt = 0.0;

                    for (int l=ao_start_[A]; l < ao_stop_[A]; l++)  {
                        for(int k=0; k < nso_; k++) {
                            Pst += 0.5 * (C_->get(0,k,s) * C_->get(0,l,t) +
                                          C_->get(0,l,s) * C_->get(0,k,t)) * S_->get(0,k,l);            // Eqn 31 (JCP 90, 4916)

                            Pss += C_->get(0,k,s) * C_->get(0,l,s) * S_->get(0,k,l);                    // Eqn 31 (JCP 90, 4916)

                            Ptt += C_->get(0,k,t) * C_->get(0,l,t) * S_->get(0,k,l);                    // Eqn 31 (JCP 90, 4916)

                        }
                        double Pss_test;
                        Pss_test = C_DDOT(nso_, &(C_->pointer(0)[0][s]), nso_, &(C_->pointer(0)[l][s]), nso_);
                    }

                    Bst += Pst * (Pss - Ptt);                                  // Eqn 29B (JCP 90, 4916)
                    Ast += Pst * Pst - 0.25 * (Pss - Ptt) * (Pss - Ptt);               // Eqn 29A (JCP 90, 4916)
                } // A-loop


                // Compute the rotation angle
                AB = (Ast * Ast) + (Bst * Bst);

                if(fabs(AB) > 0.0) {
                    cos4a = -Ast/sqrt(AB);                                     // Eqn 13b (JCP 90, 4916)
                    alpha = 0.25 * acos(cos4a) * (Bst > 0 ? 1 : -1);
                }
                else alpha = 0.0;

                // Keep up with the maximum 2x2 rotation angle
                alphamax = (fabs(alpha) > alphamax ? alpha : alphamax);


//                if (me_ == 0) {
//                    std::cout << "Ast[" << iter <<  "][" << s << "][" << t << "] = \t\t" << Ast << std::endl;
//                    std::cout << "Bst[" << iter <<  "][" << s << "][" << t << "] = \t\t" << Bst << std::endl;
//                    std::cout << "Alp[" << iter <<  "][" << s << "][" << t << "] = \t\t" << alpha << std::endl;
//                    std::cout << "Amax[" << iter <<  "][" << s << "][" << t << "] = \t\t" << alphamax << std::endl << std::endl;
//                }

                double Uss = cos(alpha);                                            // Eqn 10a/b (JCP 90, 4916)
                double Utt = cos(alpha);                                               // Eqn 10a/b (JCP 90, 4916)
                double Ust = sin(alpha);                                               // Eqn 10a/b (JCP 90, 4916)
                double Uts = -Ust;                                                  // Eqn 10a/b (JCP 90, 4916)

                // Now do the rotation
                for(int k=0; k < nso_; k++) {
                    double LCks = C_->get(0,k,s);
                    double LCkt = C_->get(0,k,t);
                    C_->set(0, k, s, (Uss * LCks + Ust * LCkt));
                    C_->set(0, k, t, (Uts * LCks + Utt * LCkt));
                }

                U->identity();

                U->set(0,s,s,Uss);
                U->set(0,t,t,Utt);
                U->set(0,s,t,Ust);
                U->set(0,t,s,Uts);

                VV->gemm('n', 't', 1.0, V, U, 0.0);

                V->copy(VV);

            } // t-loop
        } // s-loop

        conv = fabs(alphamax) - fabs(alphalast);
        if (me_ == 0) {
            fprintf(outfile, "  %2d  %16.10f  %12.3e\n", iter, alphamax, conv);
        }


        if((iter > 1) && ((fabs(conv) < 1e-12) || alphamax == 0.0))
            break;
        alphalast = alphamax;

    } // iter-loop

    timer_off("localization iterations");


    fflush(outfile);

    F_oo->zero();

    for(int i=0; i < ndocc_; i++) {
        for(int j=0; j < ndocc_; j++) {
            for(int k=0; k < ndocc_; k++) {
                F_oo->add(0,i, j, (V->get(0,k,i) * evals->get(0,k) * V->get(0,k,j)));
            }
        }
    }

    std::vector<int> orb_order(ndocc_,0);
    std::vector<int> orb_boolean(ndocc_,0);

    // First, find the overall maximum
    int i, max;
    for (i=0, max=0; i < ndocc_; i++) {
        if ( fabs(F_oo->get(0,i,i)) > fabs(F_oo->get(0,max,max)) ) {
            max = i;
        }
    }
    orb_order[0] = max;
    orb_boolean[max] = 1;

    // First, find the overall maximum
    for (i=0, max=0; i < ndocc_; i++) {
        if ( fabs(F_oo->get(0,i,i)) > fabs(F_oo->get(0,max,max)) ) {
            max = i;
        }
    }

    orb_order[0] = max;
    orb_boolean[max] = 1;

    for(i=1; i < ndocc_; i++) {
        max = 0;
        while(orb_boolean[max]) {
            max++; // Find an unused max
        }
        for(int j=0; j < ndocc_; j++) {
            if ((fabs( F_oo->get(0,j,j) ) >= fabs( F_oo->get(0,max,max) ) ) && !orb_boolean[j]) {
                max = j;
            }
        }
        orb_order[i] = max;
        orb_boolean[max] = 1;
    }

    // Now reorder the localized MO's according to F
    LCtmp = block_matrix(nso_,ndocc_);
    for(int i=0; i < ndocc_; i++) {
        C_DCOPY(nso_, &(C_->pointer(0)[0][i]), nso_, &(LCtmp[0][i]), ndocc_);
    }

    for(i=0; i < ndocc_; i++) {
        int iold = orb_order[i];
        for(int j=0; j < nso_; j++) {
            C_->set(0, j, i, LCtmp[j][iold]);
        }
        evals->set(0, i, F_oo->get(0, iold, iold));
    }

    C_->set_name("MO Coefficients (LO)");
    if(print_ >= 2) {
        C_->print(outfile);
        evals->print();
    }

    orb_order.clear();
    orb_boolean.clear();
    free_block(LCtmp);

    timer_off("Localization");

    if (me_ == 0)
        fprintf(outfile, "\n  ================================\n\n");


}

#endif // have_madness

}}

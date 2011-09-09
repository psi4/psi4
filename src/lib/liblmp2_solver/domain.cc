/*! \defgroup LMP2 lmp2: LMP2 Evaluation of Energy */

/*!
 ** \file
 ** \ingroup LMP2
 ** \LMP2 evaluation of energy
 */

#include <liblmp2_solver/lmp2.h>
#include <libqt/qt.h>
#include <physconst.h>
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


void LMP2::build_domains() {


    timer_on("Orbital Domains");

    int i,j,k,l, max;
    int xcoord, ycoord, zcoord;
    int next_atom;
    int row, col, errcod;
    double tmp;

    for (int i=0; i < ndocc_; i++) {
        domain_.push_back(std::vector<int>(natom_,0));
        domain_len_.push_back(0);
    }

    std::vector<double> charge(natom_,0.0);
    std::vector<int> rank(natom_, 0);
    std::vector<bool> boolean(natom_, false);
    std::vector<double> SR(nso_, 0.0);
    std::vector<double> Z(nso_,0.0);
    std::vector<int> ipiv(nso_,0);
    std::vector<double> fR(ndocc_,0.0);
    double **X = block_matrix(nso_,nso_);

    for (int i = 0; i < ndocc_; i++) {

        // Compute the contribution of each atom to this orbital's charge/population
        for (int A = 0; A < natom_; A++) {
            charge[A] = 0.0;
            for (int j = ao_start_[A]; j < ao_stop_[A]; j++) {
                double tmp = C_DDOT(nso_, &(S_->pointer(0)[j][0]), 1, &(C_->pointer(0)[0][i]), nso_);
                charge[A] += tmp * C_->get(0,j,i);
            }
        }

//            for(int j=0; j < glob.natom_; j++) {
//              fprintf(outfile, "charge[%d] = %20.12f\n", j, charge[j]);
//            }


        // Rank the atomic contributions to the orbital's charge
        for (j = 0; j < natom_; j++) {
            rank[j] = 0;
            boolean[j] = false;
        }
        // find the overall maximum
        for (j = 0, max = 0; j < natom_; j++) {
            if (fabs(charge[j]) >= fabs(charge[max]))
                max = j;
        }
        rank[0] = max;
        boolean[max] = 1;

        for (j = 1; j < natom_; j++) {
            max = 0;
            while (boolean[max]) max++; // find an unused max
            for (k = 0; k < natom_; k++) {
                if ((fabs(charge[k]) >= fabs(charge[max])) && !boolean[k])
                    max = k;
            }
            rank[j] = max;
            boolean[max] = 1;
        }

//        for(int j=0; j < glob.natom_; j++) {
//          fprintf(outfile, "rank[%d] = %d\n", j, rank[j]);
//        }


        // Build the orbital's domain starting in order of decreasing charge contribution
        for (j = 0; j < nso_; j++) {
            SR[j] = C_DDOT(nso_, &(S_->pointer(0)[j][0]), 1, &(C_->pointer(0)[0][i]), nso_);
        }

        domain_[i - nfocc_][rank[0]] = 1; // at least one atom must be in the domain
        domain_len_[i - nfocc_] = 1;

        fR[i - nfocc_] = 1.0;
        next_atom = 1;
        while (fabs(fR[i - nfocc_]) > cutoff_) {

            // Completeness check
            for (j = 0, row = 0; j < natom_; j++) {
                if (domain_[i - nfocc_][j]) {
                    for (k = ao_start_[j]; k < ao_stop_[j]; k++, row++) {

                        Z[row] = SR[k];

                        for (l = 0, col = 0; l < natom_; l++) {
                            if (domain_[i - nfocc_][l]) {

                                for (int m = ao_start_[l]; m < ao_stop_[l]; m++, col++)
                                    X[row][col] = S_->get(0,k,m);

                            }
                        } // l

                    } // k
                }
            } // j

            errcod = C_DGESV(row, 1, &(X[0][0]), nso_, &(ipiv[0]), &(Z[0]), nso_);
            if (errcod) {
                throw PsiException("DGESV error in LMP2 orbital domain construction", __FILE__, __LINE__);
            }
            fR[i - nfocc_] = 1.0;
            for (j = 0, row = 0; j < natom_; j++) {
                if (domain_[i - nfocc_][j]) {
                    for (k = ao_start_[j]; k < ao_stop_[j]; k++, row++) {
                        for (l = 0; l < nso_; l++)
                            fR[i - nfocc_] -= Z[row] * S_->get(0,k,l) * C_->get(0,l,i);
                    }
                }
            }

            // Augment the domain if necessary
            if (fabs(fR[i - nfocc_]) > cutoff_) {
                domain_[i - nfocc_][rank[next_atom++]] = 1;
                domain_len_[i - nfocc_]++;
            }
        } // glob.cutoff_ check
    } // i

    /* Print the orbital domains */
    if (me_ == 0) {
        fprintf(outfile, "\n  ====> Building the Domains <====\n\n");

        fprintf(outfile, "   ****** Boughton-Pulay Occupied Orbital Domains ******\n");
    }
    print_domains(fR);

//    std::vector<int> user_domain = options_.get_int_vector("DOMAINS");
//    for (int p=0; p < user_domain.size(); p++)
//        std::cout << "user_domain = " << user_domain[p] << std::endl;

//    if (glob.me_ == 0)
//        fprintf(outfile, "\n   ****** Boughton-Pulay Occupied Orbital Domains ******\n");
//    glob.print_domains(&fR[0]);


//    /* Allow user input of selected domains */
////    if (ip_exist("DOMAINS", 0)) {
////        ip_count("DOMAINS", &num_entries, 0);
////        for (int i = 0; i < num_entries; i++) {
////            ip_count("DOMAINS", &entry_len, 1, i);
////            ip_data("DOMAINS", "%d", &orbital, 2, i, 0);

////            /* Clear out the current domain for this orbital */
////            for (int j = 0; j < glob.natom_; j++) domain[orbital][j] = 0;
////            domain_len[orbital] = 0;

////            for (int j = 1; j < entry_len; j++) {
////                errcod = ip_data("DOMAINS", "%d", &atom, 2, i, j);
////                domain[orbital][atom] = 1;
////                domain_len[orbital]++;
////            }
////        }
////    }

    /* Recheck Completeness */
    for (i = 0; i < ndocc_; i++) {

        /* Build the orbital's domain starting in order of decreasing charge contribution */
        for (j = 0; j < nso_; j++) {
            SR[j] = C_DDOT(nso_, &(S_->pointer(0)[j][0]), 1, &(C_->pointer(0)[0][i]), nso_);
        }

        for (j = 0, row = 0; j < natom_; j++) {
            if (domain_[i - nfocc_][j]) {
                for (k = ao_start_[j]; k < ao_stop_[j]; k++, row++) {

                    Z[row] = SR[k];

                    for (l = 0, col = 0; l < natom_; l++) {
                        if (domain_[i - nfocc_][l]) {

                            for (int m = ao_start_[l]; m < ao_stop_[l]; m++, col++)
                                X[row][col] = S_->get(0,k,m);

                        }
                    } /* l */

                } /* k */
            }
        } /* j */

        errcod = C_DGESV(row, 1, &(X[0][0]), nso_, &(ipiv[0]), &(Z[0]), nso_);
        if (errcod) {
            //        if(Communicator::world->me() == 0)
            //          fprintf(outfile, "\nError in DGESV return in orbital domain construction.\n");
            //        exit(PSI_RETURN_FAILURE);
            throw PsiException("Error in DGESV in LMP2 orbital domain construction", __FILE__, __LINE__);
        }

    } /* i */

    /* Print the orbital domains */
    if (me_ == 0)
        fprintf(outfile, "\n   ****** Final Occupied Orbital Domains ******\n");
    print_domains(fR);

    int pairs = (ndocc_ * (ndocc_ + 1)) / 2;

    /* Determine if pairs are distant or not */
    pairdom_exist_ = std::vector<int>(pairs,0);

    /// The molecules xyz coordinates
    Matrix geometry = molecule_->geometry();

    ij_pairs_ = 0;
    for (i = 0; i < ndocc_; i++) {
        for (int m = 0; m < natom_; m++) {
            for (j = 0; j <= i; j++) {
                int ij = (i * (i + 1)) / 2 + j;
                for (int n = 0; n < natom_; n++) {
                    if (pairdom_exist_[ij] == 0) {
                        if(neglect_dp_) {
                            if (domain_[i][m] && domain_[j][n]) {
                                if(m == n || i == j) {
                                    pairdom_exist_[ij] = 1;
                                    ij_pairs_++;
                                }
                                else {
                                    xcoord = (geometry.get(0,n,0) - geometry.get(0,m,0)) * (geometry.get(0,n,0) - geometry.get(0,m,0));
                                    ycoord = (geometry.get(0,n,1) - geometry.get(0,m,1)) * (geometry.get(0,n,1) - geometry.get(0,m,1));
                                    zcoord = (geometry.get(0,n,2) - geometry.get(0,m,2)) * (geometry.get(0,n,2) - geometry.get(0,m,2));
                                    double dist = sqrt(xcoord + ycoord + zcoord) * _bohr2angstroms;
                                    if (dist < dp_cutoff_) {
                                        pairdom_exist_[ij] = 1;
                                        ij_pairs_++;
                                    }
                                }
                            }
                        }
                        else {
                            pairdom_exist_[ij] = 1;
                            ij_pairs_++;
                        }
                    }
                }
            }
        }
    }

    // This sets the maps of ij to i/j after distant pairs are negleted
    set_ij_maps();
    if(me_ == 0) {
        if(neglect_dp_) {
            fprintf(outfile, "\n   The following ij pairs are distant and will be negleted\n");
            for(int ij=0; ij < ij_pairs_; ij++) {
                if(pairdom_exist_[ij] == 0)
                    fprintf(outfile, "   i = %d \t j = %d\n", ij_i_map_[ij], ij_j_map_[ij]);
            }
        }
    }


    free_block(X);
    charge.clear();
    rank.clear();
    boolean.clear();
    SR.clear();
    Z.clear();
    ipiv.clear();
    fR.clear();

    if (only_print_domains_) {
        throw PsiException("Printing of orbital domains requested...exiting", __FILE__, __LINE__);
    }

    if (me_ == 0)
        fprintf(outfile, "\n  ================================\n\n");

    timer_off("Orbital Domains");


}

void LMP2::print_domains(const std::vector<double> &s){

  int j, cnt, max;
  double domain_tot, domain_ave;

  max = 0;
  for(int i=0; i < ndocc_; i++)
    if(domain_len_[i] > max) max = domain_len_[i];

  if(me_ == 0) {
    fprintf(outfile, "   Orbital  Domain");
    for(int i=0; i < max-2; i++) fprintf(outfile, "   "); /* formatting junk */
    fprintf(outfile, "  Completeness\n");
    fprintf(outfile, "   -------  ------");
    for(int i=0; i < max-2; i++) fprintf(outfile, "---"); /* more formatting junk */
    fprintf(outfile, "  ------------\n");
    for(int i=0; i < ndocc_; i++) {
      fprintf(outfile, "      %2d    ",i);
      for(j=0,cnt=0; j < natom_; j++) if(domain_[i][j]) { fprintf(outfile, " %2d", j); cnt++; }
      if(cnt < max) for(; cnt < max; cnt++) fprintf(outfile, "   ");
      fprintf(outfile, "     %7.5f\n", s[i]);
    }
    domain_tot = 0;
    for(int i=0; i < ndocc_; i++)
      domain_tot += domain_len_[i];
    domain_ave = domain_tot/ndocc_;
    fprintf(outfile, "\n   The average domain length is %4.2lf\n", domain_ave);    
    fflush(outfile);
  }
}

#endif // have_madness
}}


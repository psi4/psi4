#include "3index.h"
#include "libmints/mints.h"
#include <libqt/qt.h>

#include <string>
#include <sstream>
#include <iostream>
#include <cstdio>
#include <fstream>
#include <algorithm>
#include <utility>
#include <ctype.h>

using namespace boost;
using namespace std;
using namespace psi;

namespace psi {

DealiasBasisSet::DealiasBasisSet(boost::shared_ptr<BasisSet> primary, Options& options) :
    primary_(primary), options_(options)
{
    setDelta(options_.get_double("DEALIAS_DELTA"));
    setBeta(options_.get_double("DEALIAS_BETA"));
    setNCore(options_.get_int("DEALIAS_N_CORE"));
    setNIntercalater(options_.get_int("DEALIAS_N_INTERCALATER"));
    setNDiffuse(options_.get_int("DEALIAS_N_DIFFUSE"));
    setNCap(options_.get_int("DEALIAS_N_CAP"));
    setNL(options_.get_int("DEALIAS_N_L"));

    buildDealiasBasisSet();
}

DealiasBasisSet::~DealiasBasisSet()
{
}

void DealiasBasisSet::buildDealiasBasisSet()
{
    form_primary_alpha();
    form_core();
    form_intercalater();
    form_diffuse();
    form_cap();
    form_basis();
}

void DealiasBasisSet::form_primary_alpha()
{
    int natom = primary_->molecule()->natom();
    int max_am = primary_->max_am();

    primary_alpha_.clear();
    dealias_alpha_.clear();

    primary_alpha_.resize(natom);
    dealias_alpha_.resize(natom);

    // Little bigger than needed, but no big deal
    for (int A = 0; A < natom; A++) {
        primary_alpha_[A].resize(max_am + 1);
        dealias_alpha_[A].resize(max_am + 1 + nl_);
    }

    // Geometric mean of all contractions to extract a single alpha
    for (int M = 0; M < primary_->nshell(); M++) {
        const GaussianShell& shell = primary_->shell(M);
        int A = shell.ncenter();
        int l = shell.am();
        int nprim = shell.nprimitive();

        double numerator = 0.0;
        double denominator = 0.0;

        for (int K = 0; K < nprim; K++) {
            double c = shell.coef(K);
            double a = shell.exp(K);
            numerator   += c * log(a);
            denominator += c;
        }

        primary_alpha_[A][l].push_back(exp(numerator / denominator));
    }

    // Sort largest to smallest
    for (int A = 0; A < natom; A++) {
        for (int l = 0; l <= max_am; l++) {
            if (primary_alpha_[A][l].size() == 0) continue;
            std::sort(primary_alpha_[A][l].begin(), primary_alpha_[A][l].end(), std::greater<double>());
        }
    }
}

void DealiasBasisSet::form_core()
{
    int natom = primary_->molecule()->natom();
    int max_am = primary_->max_am();
    for (int A = 0; A < natom; A++) {
        for (int l = 0; l <= max_am; l++) {
            if (primary_alpha_[A][l].size() == 0) continue;
            double alpha_c = primary_alpha_[A][l][0];
            for (int i = 0; i < primary_alpha_[A][l].size(); i++) {
                alpha_c = (alpha_c < primary_alpha_[A][l][i] ? primary_alpha_[A][l][i] : alpha_c);
            }
            for (int j = ncore_; j >= 1; j--) {
                dealias_alpha_[A][l].push_back(alpha_c * pow(delta_,(double) j));
            }
        }
    }
}

void DealiasBasisSet::form_intercalater()
{
    int natom = primary_->molecule()->natom();
    int max_am = primary_->max_am();
    for (int A = 0; A < natom; A++) {
        for (int l = 0; l <= max_am; l++) {
            if (primary_alpha_[A][l].size() < 2) continue;
            for (int i = 0; i < primary_alpha_[A][l].size() - 1; i++) {
                double alpha_max = primary_alpha_[A][l][i];
                double alpha_min = primary_alpha_[A][l][i + 1];
                double e_max = log(alpha_max);
                double e_min = log(alpha_min);
                for (int j = nintercalater_; j >= 1; j--) {
                    double e = e_min + (e_max - e_min) * ((double) j) / ((double) nintercalater_ + 1);
                    dealias_alpha_[A][l].push_back(exp(e));
                }
            }
        }
    }
}

void DealiasBasisSet::form_diffuse()
{
    int natom = primary_->molecule()->natom();
    int max_am = primary_->max_am();
    for (int A = 0; A < natom; A++) {
        for (int l = 0; l <= max_am; l++) {
            if (primary_alpha_[A][l].size() == 0) continue;
            double alpha_c = primary_alpha_[A][l][0];
            for (int i = 0; i < primary_alpha_[A][l].size(); i++) {
                alpha_c = (alpha_c > primary_alpha_[A][l][i] ? primary_alpha_[A][l][i] : alpha_c);
            }
            for (int j = 1; j <= ndiffuse_; j++) {
                dealias_alpha_[A][l].push_back(alpha_c * pow(1.0 / delta_,(double) j));
            }
        }
    }
}

void DealiasBasisSet::form_cap()
{
    int natom = primary_->molecule()->natom();
    int max_am = primary_->max_am();

    for (int A = 0; A < natom; A++) {
        for (int l = max_am; l >= 0; l--) {
            if (primary_alpha_[A][l].size() > 0) {

                double numerator = 0.0;
                double denominator = 0.0;

                for (int i = 0; i < primary_alpha_[A][l].size(); i++) {
                    double a = primary_alpha_[A][l][i];
                    numerator   += log(a);
                    denominator += 1.0;
                }

                double alpha_c = exp(numerator / denominator);

                for (int dl = 1; dl <= nl_; dl++) {
                    int nfun = ncap_ + nl_ - dl;
                    for (int i = 1; i <= nfun; i++) {
                        dealias_alpha_[A][l+dl].push_back(alpha_c * pow(beta_, ((double) nfun + 1.0) / 2.0 - (double) i));
                    }
                }

                break;
            }
        }
    }
}

void DealiasBasisSet::form_basis()
{
    std::vector<GaussianShell> shells;
    int natom = primary_->molecule()->natom();
    int max_am = primary_->max_am();
    int max_l = max_am + nl_;
    std::vector<double> weight(1);
    weight.push_back(1.0);

    for (int A = 0; A < natom; A++) {
        for (int l = 0; l <= max_l; l++) {
            for (int i = 0; i < dealias_alpha_[A][l].size(); i++) {
                std::vector<double> e(1);
                e.push_back(dealias_alpha_[A][l][i]);
                Vector3 v = primary_->molecule()->xyz(A);
                shells.push_back(GaussianShell(l, weight, e, primary_->has_puream() ? Pure : Cartesian, A, v, 0, Unnormalized));
            }
        }
    }

    dealias_ = BasisSet::build(primary_->molecule(), shells);
}

}

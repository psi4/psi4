/*
 * frac.cc. So close to a great Battlestar Galactica joke. 
 * The way this works:
 * 1) form_C returns 1's normalized C matrices.
 * 2) find_occupation determines the occupations of said matrix via either Aufbau or MOM selection
 * 3) find_occupation then calls this, which renormalizes the C matrices with \sqrt(val) for each frac occ
 * 4) find_D is then computed with the renormalized C matrices, and everything is transparent until 1) on the next cycle
 *
 * Some executive decisions:
 *  -DIIS: Upon FRAC start, the old DIIS info is nuked, and DIIS_START is incremented by the current iteration count. 
 *         Thus, DIIS begins again on the next iteration. The exception is if you set FRAC_DIIS to false,
 *         in which case DIIS will cease for good upon FRAC start. 
 *  -MOM:  To use MOM with FRAC, set MOM_START to either FRAC_START or FRAC_START+1 (I think I prefer the latter).
 */

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>

#include <libmints/mints.h>
#include <libqt/qt.h>

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include "hf.h"

using namespace boost;
using namespace std;
using namespace psi;

namespace psi { namespace scf {

void HF::frac()
{
    // Perhaps no frac?
    if (iteration_ < options_.get_int("FRAC_START") || options_.get_int("FRAC_START") == 0)  return;
    frac_performed_ = true;    

    // First frac iteration, blow away the diis and print the frac task
    if (iteration_ == options_.get_int("FRAC_START")) {

        // Throw unless UHF/UKS
        if (!(options_.get_str("REFERENCE") == "UHF" || options_.get_str("REFERENCE") == "UKS"))
            throw PSIEXCEPTION("Fractional Occupation SCF is only implemented for UHF/UKS");

        // Throw if no frac tasks
        if (!options_["FRAC_OCC"].size())
            throw PSIEXCEPTION("Fractional Occupation SCF requested, but empty FRAC_OCC/FRAC_VAL vector");

        // Throw if inconsistent size
        if (options_["FRAC_OCC"].size() != options_["FRAC_VAL"].size()) 
            throw PSIEXCEPTION("Fractional Occupation SCF: FRAC_OCC/FRAC_VAL are of different dimensions");

        // Throw if the user is being an idiot with docc/socc
        if (input_docc_ || input_socc_) 
            throw PSIEXCEPTION("Fractional Occupation SCF: Turn off DOCC/SOCC");

        // Throw if the user is trying to start MOM before FRAC
        if (options_.get_int("MOM_START") <= options_.get_int("FRAC_START") && options_.get_int("MOM_START") != 0)
            throw PSIEXCEPTION("Fractional Occupation SCF: MOM must start after FRAC");

        // Throw if the use is just way too eager
        if (MOM_excited_) 
            throw PSIEXCEPTION("Fractional Occupation SCF: Don't try an excited-state MOM");
    
        // Close off a previous burn-in SCF
        fprintf(outfile, "\n");
        print_orbitals(); 
        
        // frac header
        fprintf(outfile, "\n  ==> Fractionally-Occupied SCF Iterations <==\n\n");
        for (int ind = 0; ind < options_["FRAC_OCC"].size(); ind++) {
            int i = options_["FRAC_OCC"][ind].to_integer();
            double val = options_["FRAC_VAL"][ind].to_double(); 
            
            // Throw if user requests frac occ above nalpha/nbeta 
            int max_i = (i > 0 ? nalpha_: nbeta_);
            if(abs(i) > max_i) {
                if (i > 0)
                    nalpha_++;
                else 
                    nbeta_++;
            }

            // Throw if the user is insane
            if (val < 0.0)
                throw PSIEXCEPTION("Fractional Occupation SCF: PSI4 is not configured for positrons. Please annihilate and start again");             

            fprintf(outfile, "    %-5s orbital %4d will contain %11.3E electron.\n", (i > 0 ? "Alpha" : "Beta"), abs(i), val);
        }
        fprintf(outfile, "\n");

        // Make sure diis restarts correctly/frac plays well with MOM
        if (initialized_diis_manager_) {
            diis_manager_->delete_diis_file();
            diis_manager_.reset();
            initialized_diis_manager_ = false;
            diis_start_ += iteration_ + 1;
        }

        // Turn yonder DIIS off if requested
        if (!options_.get_bool("FRAC_DIIS")) {
            diis_enabled_ = false;
        }

        // Load the old orbitals in if requested
        if (options_.get_bool("FRAC_LOAD")) {
            fprintf(outfile, "    Orbitals reloaded from file, your previous iterations are garbage.\n\n");
            load_orbitals();
        }

        // Keep the printing nice
        fprintf(outfile, "                        Total Energy        Delta E      Density RMS\n\n");
        fflush(outfile);

        // Prevent spurious convergence (technically this iteration comes from the N-electron system anyways)
        frac_performed_ = false;
    }
    
    // Every frac iteration: renormalize the Ca/Cb matrices

    // Sort the eigenvalues in the usual manner
    std::vector<boost::tuple<double,int,int> > pairs_a;
    std::vector<boost::tuple<double,int,int> > pairs_b;
    for (int h=0; h<epsilon_a_->nirrep(); ++h) {
        for (int i=0; i<epsilon_a_->dimpi()[h]; ++i)
            pairs_a.push_back(boost::tuple<double,int,int>(epsilon_a_->get(h, i), h, i));
    }
    for (int h=0; h<epsilon_b_->nirrep(); ++h) {
        for (int i=0; i<epsilon_b_->dimpi()[h]; ++i)
            pairs_b.push_back(boost::tuple<double,int,int>(epsilon_b_->get(h, i), h, i));
    }
    sort(pairs_a.begin(),pairs_a.end());
    sort(pairs_b.begin(),pairs_b.end());
    
    // Renormalize the C matrix entries
    for (int ind = 0; ind < options_["FRAC_OCC"].size(); ind++) {
        int i = options_["FRAC_OCC"][ind].to_integer();
        double val = options_["FRAC_VAL"][ind].to_double(); 
        bool is_alpha = (i > 0);
        i = abs(i) - 1; // Back to C ordering

        int i2  = ((is_alpha) ? get<2>(pairs_a[i]) : get<2>(pairs_b[i])); 
        int h   = ((is_alpha) ? get<1>(pairs_a[i]) : get<1>(pairs_b[i])); 

        int nso = Ca_->rowspi()[h];
        int nmo = Ca_->colspi()[h];

        double** Cp = ((is_alpha) ? Ca_->pointer(h) : Cb_->pointer(h));

        // And I say all that to say this
        C_DSCAL(nso, sqrt(val), &Cp[0][i], nmo); 
    } 
}
void HF::frac_renormalize()
{
    if (!options_.get_bool("FRAC_RENORMALIZE") || !frac_enabled_) return;

    // Renormalize the fractional occupations back to 1, if possible before storage 
    fprintf(outfile, "    FRAC: Renormalizing orbitals to 1.0 for storage.\n\n");

    // Sort the eigenvalues in the usual manner
    std::vector<boost::tuple<double,int,int> > pairs_a;
    std::vector<boost::tuple<double,int,int> > pairs_b;
    for (int h=0; h<epsilon_a_->nirrep(); ++h) {
        for (int i=0; i<epsilon_a_->dimpi()[h]; ++i)
            pairs_a.push_back(boost::tuple<double,int,int>(epsilon_a_->get(h, i), h, i));
    }
    for (int h=0; h<epsilon_b_->nirrep(); ++h) {
        for (int i=0; i<epsilon_b_->dimpi()[h]; ++i)
            pairs_b.push_back(boost::tuple<double,int,int>(epsilon_b_->get(h, i), h, i));
    }
    sort(pairs_a.begin(),pairs_a.end());
    sort(pairs_b.begin(),pairs_b.end());
    
    // Renormalize the C matrix entries
    for (int ind = 0; ind < options_["FRAC_OCC"].size(); ind++) {
        int i = options_["FRAC_OCC"][ind].to_integer();
        double val = options_["FRAC_VAL"][ind].to_double(); 
        bool is_alpha = (i > 0);
        i = abs(i) - 1; // Back to C ordering

        int i2  = ((is_alpha) ? get<2>(pairs_a[i]) : get<2>(pairs_b[i])); 
        int h   = ((is_alpha) ? get<1>(pairs_a[i]) : get<1>(pairs_b[i])); 

        int nso = Ca_->rowspi()[h];
        int nmo = Ca_->colspi()[h];

        double** Cp = ((is_alpha) ? Ca_->pointer(h) : Cb_->pointer(h));

        // And I say all that to say this
        if (val != 0.0) 
            C_DSCAL(nso, 1.0 / sqrt(val), &Cp[0][i], nmo); 
    } 
}

void HF::compute_spin_contamination()
{
    if (!(options_.get_str("REFERENCE") == "UHF" || options_.get_str("REFERENCE") == "UKS" || options_.get_str("REFERENCE") == "CUHF"))
        return;

    double nalpha = (double) nalpha_;
    double nbeta  = (double) nbeta_;

    // Adjust for fractional occupation
    if (frac_performed_) {
        for (int ind = 0; ind < options_["FRAC_OCC"].size(); ind++) {
            int i = options_["FRAC_OCC"][ind].to_integer();
            double val = options_["FRAC_VAL"][ind].to_double(); 
            bool is_alpha = (i > 0);
            if (is_alpha) {
                nalpha -= (1.0 - val);
            } else {
                nbeta  -= (1.0 - val);
            }
        }    
    }

    SharedMatrix S = SharedMatrix(factory_->create_matrix("S (Overlap)"));
    boost::shared_ptr<IntegralFactory> fact(new IntegralFactory(basisset_,basisset_, basisset_,basisset_));
    boost::shared_ptr<OneBodySOInt> so_overlap(fact->so_overlap());
    so_overlap->compute(S);

    double dN = 0.0;

    for (int h =0; h < S->nirrep(); h++) {
        int nbf = S->colspi()[h];
        int nmo = Ca_->colspi()[h];
        int na = nalphapi_[h];
        int nb = nbetapi_[h];
        if (na == 0 || nb == 0 || nbf == 0 || nmo == 0)
            continue;

        SharedMatrix Ht (new Matrix("H Temp", nbf, nb));
        SharedMatrix Ft (new Matrix("F Temp", na, nb));

        double** Sp = S->pointer(h);
        double** Cap = Ca_->pointer(h);
        double** Cbp = Cb_->pointer(h);
        double** Htp = Ht->pointer(0);
        double** Ftp = Ft->pointer(0);

        C_DGEMM('N','N',nbf,nb,nbf,1.0,Sp[0],nbf,Cbp[0],nmo,0.0,Htp[0],nb);
        C_DGEMM('T','N',na,nb,nbf,1.0,Cap[0],nmo,Htp[0],nb,0.0,Ftp[0],nb);

        dN += C_DDOT(na*(long int)nb, Ftp[0], 1, Ftp[0], 1);
    }

    double nmin = (nbeta < nalpha ? nbeta : nalpha);
    double dS = nmin - dN;
    double nm = (nalpha - nbeta) / 2.0;
    double S2 = fabs(nm) * (fabs(nm) + 1.0);

    fprintf(outfile, "   @Spin Contamination Metric: %17.9E\n", dS);
    fprintf(outfile, "   @S^2 Expected:              %17.9E\n", S2);
    fprintf(outfile, "   @S^2 Observed:              %17.9E\n", S2 + dS);
    fprintf(outfile, "   @S   Expected:              %17.9E\n", nm);
    fprintf(outfile, "   @S   Observed:              %17.9E\n", nm);

    if (frac_performed_) {
        fprintf(outfile, "   @Nalpha:                    %17.9E\n", nalpha);
        fprintf(outfile, "   @Nbeta:                     %17.9E\n", nbeta);
    }
    fprintf(outfile, "\n");
}

}}

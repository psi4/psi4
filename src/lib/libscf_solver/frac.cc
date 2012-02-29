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

}}

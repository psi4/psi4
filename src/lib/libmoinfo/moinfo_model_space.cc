#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>

#include <psifiles.h>
#include <exception.h>
#include <liboptions/liboptions.h>
#include <libutil/libutil.h>

#include "moinfo.h"

extern FILE *outfile;

using namespace std;

namespace psi {

std::string MOInfo::get_determinant_label(int i)
{
    return references[i].get_label();
}

void MOInfo::print_model_space()
{
    fprintf(outfile,"\n");
    fprintf(outfile,"\n  Model space");
    fprintf(outfile,"\n  ------------------------------------------------------------------------------");
    if(references.size() <= 20){
        for(int i = 0; i < references.size(); ++i){
            fprintf(outfile,"\n  %2d  %s",i,references[i].get_label().c_str());
        }
    }else{
        fprintf(outfile,"\n  There are %d determinants in the model space",static_cast<int>(references.size()));
    }
    fprintf(outfile,"\n  ==============================================================================\n");
}

void MOInfo::build_model_space()
{
    /********************************************************
    Generate all the Slater Determinants belonging to the
    model space using the following restrictions:
    -docc are doubly    occupied
    -actv are partially occupied
    -the generalized occupied orbital indexing (docc + actv)
     is assumed (see moinfo.cpp)
  ********************************************************/
    int index;
    MOInfo::SlaterDeterminant docc_det(this);

    /***********************************************
    Generate all combinations of active orbitals
  ***********************************************/
    if(nactv_docc == 0){
        /********************************************************
      Set up the doubly occupied part of all the
      determinants in the model space
    ********************************************************/
        index = 0;
        for(int h = 0; h < nirreps; ++h){
            for(int i = 0; i < docc[h]; ++i){
                docc_det.set(index);
                docc_det.set(index + nall);
                index++;
            }
            index += actv[h];
            index += extr[h];
        }

        /********************************************************
      Set up the a vectors containing the active orbitals and
      their symmetry
    ********************************************************/
        std::vector<int> alpha_active,alpha_active_sym,beta_active,beta_active_sym;
        index = 0;
        for(int h = 0; h < nirreps; ++h){
            index += docc[h];
            for(int i = 0; i < actv[h]; ++i){
                alpha_active.push_back(index);
                alpha_active_sym.push_back(h);
                beta_active.push_back(index + nall);
                beta_active_sym.push_back(h);
                index++;
            }
            index += extr[h];
        }

        std::vector<std::vector<int> > alpha_combinations,beta_combinations;
        generate_combinations(nactv,nactive_ael,alpha_combinations);
        generate_combinations(nactv,nactive_bel,beta_combinations);
        if(alpha_combinations.size()==0)
            alpha_combinations.push_back(vector<int>(0));
        if(beta_combinations.size()==0)
            beta_combinations.push_back(vector<int>(0));
        for(int a=0;a<alpha_combinations.size();a++){
            for(int b=0;b<beta_combinations.size();b++){
                int sym = 0; // Symmetry of the determinant
                // Create a copy of the docc_det
                SlaterDeterminant det(docc_det);
                // Fill the alpha active orbitals
                for(int i=0;i<nactive_ael;i++){
                    det.set(alpha_active[alpha_combinations[a][i]]);
                    sym = sym ^ alpha_active_sym[alpha_combinations[a][i]];
                }
                // Fill the beta active orbitals
                for(int i=0;i<nactive_bel;i++){
                    det.set(beta_active[beta_combinations[b][i]]);
                    sym = sym ^  beta_active_sym[beta_combinations[b][i]];
                }
                // Test the wfn symmetry
                if(sym==wfn_sym){
                    if(det.is_closed_shell()){
                        // Closed-shell determinant
                        closed_shell_refs.push_back(references.size());
                        unique_refs.push_back(references.size());
                        all_refs.push_back(references.size());
                    }else{
                        // Open-shell determinant
                        bool add_it = true;
                        int  spin_mirror = references.size();
                        if(options.get_bool("USE_SPIN_SYMMETRY")){
                            // Check if this is a spin-flipped determinant
                            for(int ref=0;ref<references.size();ref++){
                                if(references[ref].is_spin_flipped(det)){
                                    add_it      = false;
                                    spin_mirror = ref;
                                }
                            }
                            if(add_it){
                                unique_open_shell_refs.push_back(references.size());
                                unique_refs.push_back(references.size());
                                all_refs.push_back(references.size());
                            }else{
                                all_refs.push_back(spin_mirror);
                            }
                        }
                    }
                    references.push_back(det);
                }
            }
        }
    }else{
        /********************************************************
      Set up the doubly occupied part of all the
      determinants in the model space
    ********************************************************/
        index = 0;
        for(int h=0;h<nirreps;h++){
            for(int i=0;i<docc[h] + actv_docc[h];i++){
                docc_det.set(index);
                docc_det.set(index + nall);
                index++;
            }
            index += actv[h] - actv_docc[h];
            index += extr[h];
        }
        closed_shell_refs.push_back(references.size());
        unique_refs.push_back(references.size());
        all_refs.push_back(references.size());
        references.push_back(docc_det);
    }

    if(references.size() == 0){
        fprintf(outfile,"\n\n  MOInfo found no reference in the model space");
        fprintf(outfile,"\n  Please check the following:");
        fprintf(outfile,"\n  1) Definition of FROZEN_DOCC, RESTRICTED_DOCC, ACTIVE, and FROZEN_UOCC");
//        fprintf(outfile,"\n  1) Definition of FOCC, DOCC, ACTV, and FVIR");
        fprintf(outfile,"\n  2) Symmetry of the wavefunction");
        fprintf(outfile,"\n  3) Charge and multiplicity");
        fprintf(outfile,"\n\n  PSIMRCC will end the computation.\n");
        fflush(outfile);
        exit(PSI_RETURN_FAILURE);
    }
}

/*!
    \fn MOInfo::make_internal_excitations()
 */
void MOInfo::make_internal_excitations()
{
    /****************************************
    Build the mappings between references
    |m> = (+/-) ... b+ j a+ i |n>
  ****************************************/
    for(int m = 0; m < references.size(); ++m){
        vector<vector<pair<int,int> > > alpha_internals_ref_m;
        vector<vector<pair<int,int> > >  beta_internals_ref_m;
        vector<double>                   sign_internals_ref_m;
        //       fprintf(outfile,"\n\n\tReference %s",references[m].get_label().c_str());
        //       fprintf(outfile," gives:");
        for(int n=0;n<references.size();n++){
            double sign=1.0;
            std::vector<pair<int,int> > alpha_operators;
            std::vector<pair<int,int> > beta_operators;
            references[m].get_internal_excitations(references[n],sign,alpha_operators,beta_operators);
            alpha_internals_ref_m.push_back(alpha_operators);
            beta_internals_ref_m.push_back(beta_operators);
            sign_internals_ref_m.push_back(sign);
            //         fprintf(outfile,"\n\t  %s",references[n].get_label().c_str());
            //         fprintf(outfile," = %s{",sign > 0.0 ? "+" : (sign == 0.0 ? "0" : "-"));
            //         for(int i = 0; i<beta_operators.size();i++)
            //           fprintf(outfile," %db+ %db-",beta_operators[i].second,beta_operators[i].first);
            //         for(int i = 0; i<alpha_operators.size();i++)
            //           fprintf(outfile," %da+ %da-",alpha_operators[i].second,alpha_operators[i].first);
            //         fprintf(outfile," }");
            //         fprintf(outfile,"%s",references[m].get_label().c_str());
        }
        alpha_internal_excitations.push_back(alpha_internals_ref_m);
        beta_internal_excitations.push_back(beta_internals_ref_m);
        sign_internal_excitations.push_back(sign_internals_ref_m);
    }
}

vector<int> MOInfo::get_determinant(int i)
{
    vector<int> occupation(nall * 2,0);
    for(int p = 0; p < 2 * nall; ++p)
        if(references[i].test(p))
            occupation[p] = 1;
    return occupation;
}

vector<pair<int,int> > MOInfo::get_alpha_internal_excitation(int i,int j)
{
    return(alpha_internal_excitations[i][j]);
}

vector<pair<int,int> > MOInfo::get_beta_internal_excitation(int i,int j)
{
    return(beta_internal_excitations[i][j]);
}

double  MOInfo::get_sign_internal_excitation(int i,int j)
{
    return(sign_internal_excitations[i][j]);
}

/*!
 *  \fn MOInfo::get_ref_number(ReferenceType ref_type, int n)
 */
int MOInfo::get_ref_number(int n, ReferenceType ref_type)
{
    if(ref_type == AllRefs) // a
        return(all_refs[n]);
    if(ref_type == UniqueRefs) // u
        return(unique_refs[n]);
    if(ref_type == ClosedShellRefs) // c
        return(closed_shell_refs[n]);
    if(ref_type == UniqueOpenShellRefs) // o
        return(unique_open_shell_refs[n]);
    throw PSIEXCEPTION("MOInfo::get_ref_number(string str, int n) undefined space");
    return(NULL);
}

/*!
    \fn MOInfo::get_ref_size(string str)
 */
int MOInfo::get_ref_size(ReferenceType ref_type)
{
    if(ref_type == AllRefs) // a
        return(all_refs.size());
    if(ref_type == UniqueRefs) // u
        return(unique_refs.size());
    if(ref_type == ClosedShellRefs) // c
        return(closed_shell_refs.size());
    if(ref_type == UniqueOpenShellRefs) // o
        return(unique_open_shell_refs.size());
    throw PSIEXCEPTION("MOInfo::get_ref_size(string str) undefined space");
    return(NULL);
}

vector<string> MOInfo::get_matrix_names(std::string str)
{
    vector<string> names;
    if(str.find("{a}")!=string::npos){
        for(int n=0;n<all_refs.size();n++)
            names.push_back(find_and_replace(str,"{a}","{" + to_string(all_refs[n]) +"}"));
    }else if(str.find("{u}")!=string::npos){
        for(int n=0;n<unique_refs.size();n++)
            names.push_back(find_and_replace(str,"{u}","{" + to_string(unique_refs[n]) +"}"));
    }else if(str.find("{c}")!=string::npos){
        for(int n=0;n<closed_shell_refs.size();n++)
            names.push_back(find_and_replace(str,"{c}","{" + to_string(closed_shell_refs[n]) +"}"));
    }else if(str.find("{o}")!=string::npos){
        for(int n=0;n<unique_open_shell_refs.size();n++)
            names.push_back(find_and_replace(str,"{o}","{" + to_string(unique_open_shell_refs[n]) +"}"));
    }else
        names.push_back(str);
    return(names);
}

vector<int> MOInfo::get_aocc(int i,ReferenceType ref_type)
{
    int i_ref = get_ref_number(i,ref_type);
    return(references[i_ref].get_aocc());
}

vector<int> MOInfo::get_bocc(int i,ReferenceType ref_type)
{
    int i_ref = get_ref_number(i,ref_type);
    return(references[i_ref].get_bocc());
}

vector<int> MOInfo::get_avir(int i,ReferenceType ref_type)
{
    int i_ref = get_ref_number(i,ref_type);
    return(references[i_ref].get_avir());
}

vector<int> MOInfo::get_bvir(int i,ReferenceType ref_type)
{
    int i_ref = get_ref_number(i,ref_type);
    return(references[i_ref].get_bvir());
}

vector<bool> MOInfo::get_is_aocc(int i,ReferenceType ref_type)
{
    int i_ref = get_ref_number(i,ref_type);
    return(references[i_ref].get_is_aocc());
}

vector<bool> MOInfo::get_is_bocc(int i,ReferenceType ref_type)
{
    int i_ref = get_ref_number(i,ref_type);
    return(references[i_ref].get_is_bocc());
}

vector<bool> MOInfo::get_is_avir(int i,ReferenceType ref_type)
{
    int i_ref = get_ref_number(i,ref_type);
    return(references[i_ref].get_is_avir());
}

vector<bool> MOInfo::get_is_bvir(int i,ReferenceType ref_type)
{
    int i_ref = get_ref_number(i,ref_type);
    return(references[i_ref].get_is_bvir());
}

}

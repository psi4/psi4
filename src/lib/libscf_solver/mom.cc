/*
 *  hf.cpp
 *  matrix
 *
 *  Created by Justin Turney on 4/9/08.
 *  Copyright 2008 by Justin M. Turney, Ph.D.. All rights reserved.
 *
 */

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>

#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.hpp>
#include <libparallel/parallel.h>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>
#include <psifiles.h>
#include "integralfunctors.h"
#include "pseudospectral.h"
#include "pkintegrals.h"
#include "df.h"

#include <libmints/mints.h>

#include "hf.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace boost;
using namespace std;
using namespace psi;

namespace psi { namespace scf {

void HF::MOM_start()
{
    // Perhaps no MOM (or at least no MOM_start())? 
    if (iteration_ != options_.get_int("MOM_START") || iteration_ == 0 || MOM_started_)  return;
    
    // If we're here, we're doing MOM of some kind
    MOM_started_ = true;
    MOM_performed_ = true; // Gets printed next iteration   
    //
    // Build Ca_old_ matrices
    Ca_old_ = SharedMatrix(new Matrix("C Alpha Old (SO Basis)", nirrep_, nsopi_, nmopi_));
    if (!restricted()) {
        Cb_old_ = SharedMatrix(new Matrix("C Beta Old (SO Basis)", nirrep_, nsopi_, nmopi_));
    } else {
        Cb_old_ = Ca_old_;
    }
    
    Ca_old_->copy(Ca_);
    Cb_old_->copy(Cb_);

    // If no excitation requested, it's a stabilizer MOM, nothing fancy, don't print 
    if (!options_["MOM_OCC"].size()) return;      
   
    // If we're here, its an exciting MOM
    fprintf(outfile, "\n  ==> MOM Excited-State Iterations <==\n\n");
    
    // Reset iterations and DIIS (will automagically restart)
    iteration_ = 0;
    if (initialized_diis_manager_) {
        diis_manager_->delete_diis_file();
        diis_manager_.reset();
        initialized_diis_manager_ = false;
    }
 
    // Find out which orbitals are where
    std::vector<std::pair<double, std::pair<int, int> > > orbs_a;
    for (int h = 0; h < nirrep_; h++) {
        int nmo = nmopi_[h];
        if (nmo == 0) continue;
        double* eps = epsilon_a_->pointer(h);
        for (int a = 0; a < nmo; a++)
            orbs_a.push_back(make_pair(eps[a], make_pair(h, a)));
    }
    std::vector<std::pair<double, std::pair<int, int> > > orbs_b;
    for (int h = 0; h < nirrep_; h++) {
        int nmo = nmopi_[h];
        if (nmo == 0) continue;
        double* eps = epsilon_b_->pointer(h);
        for (int a = 0; a < nmo; a++)
            orbs_b.push_back(make_pair(eps[a], make_pair(h, a)));
    }
    sort(orbs_a.begin(),orbs_a.end());
    sort(orbs_b.begin(),orbs_b.end());

    if (options_["MOM_OCC"].size() != options_["MOM_VIR"].size())
        throw PSIEXCEPTION("SCF: MOM_OCC and MOM_VIR are not the same size");
    
    std::set<int> tot_check;

    for (int n = 0; n < options_["MOM_OCC"].size(); n++) {
        tot_check.insert(options_["MOM_OCC"][n].to_integer());
        tot_check.insert(options_["MOM_VIR"][n].to_integer());
    }

    if (tot_check.size() != 2*options_["MOM_OCC"].size())
        throw PSIEXCEPTION("SCF::MOM_start: Occupied/Virtual index array elements are not unique");

    CharacterTable ct = molecule_->point_group()->char_table();

    fprintf(outfile, "  Excitations:\n");

    for (int n = 0; n < options_["MOM_OCC"].size(); n++) {
        int i = options_["MOM_OCC"][n].to_integer(); 
        int a = options_["MOM_VIR"][n].to_integer();
       
        bool si = (i > 0);
        bool sa = (a > 0);
    
        i = abs(i) - 1;
        a = abs(a) - 1;

        if (options_.get_str("REFERENCE") == "RHF" || options_.get_str("REFERENCE") == "RKS") {

            if (!si || !sa) 
                throw PSIEXCEPTION("SCF::MOM_start: RHF/RKS requires + -> + in input, as only double excitations are performed");

            int hi = orbs_a[i].second.first;
            int ha = orbs_a[a].second.first;
 
            if (hi == ha) {
                // Same irrep
                int pi = orbs_a[i].second.second;
                int pa = orbs_a[a].second.second;
  
                int nso = nsopi_[hi];
                int nmo = nmopi_[hi];

                double** Ca = Ca_->pointer(hi);
                double* Ct = new double[nso]; 
                double* eps = epsilon_a_->pointer(hi);
                double  epst;
     
                // Swap eigvals 
                epst = eps[pi];
                eps[pi] = eps[pa];
                eps[pa] = epst; 

                // Swap eigvecs
                C_DCOPY(nso, &Ca[0][pi], nmo, Ct, 1);
                C_DCOPY(nso, &Ca[0][pa], nmo, &Ca[0][pi], nmo);
                C_DCOPY(nso, Ct, 1, &Ca[0][pa], nmo);
            
                delete[] Ct;

                fprintf(outfile, "   %8s: %4d%-4s -> %4d%-4s \n", "AB -> AB", pi + 1, ct.gamma(hi).symbol(), pa + 1, ct.gamma(ha).symbol()); 
            } else {
                // Different irrep
                // Occ -> Vir
                int pi = orbs_a[i].second.second;
                int pi2 = nalphapi_[hi] - 1;  
 
                int nso = nsopi_[hi];
                int nmo = nmopi_[hi];

                double** Ca = Ca_->pointer(hi);
                double* Ct = new double[nso]; 
                double* eps = epsilon_a_->pointer(hi);
                double  epst;
     
                // Swap eigvals 
                epst = eps[pi];
                eps[pi] = eps[pi2];
                eps[pi2] = epst; 

                // Swap eigvecs
                C_DCOPY(nso, &Ca[0][pi], nmo, Ct, 1);
                C_DCOPY(nso, &Ca[0][pi2], nmo, &Ca[0][pi], nmo);
                C_DCOPY(nso, Ct, 1, &Ca[0][pi2], nmo);
            
                delete[] Ct;
                
                // Redo indexing
                nalphapi_[hi]--;
                nbetapi_[hi]--;
                doccpi_[hi]--; 

                // Vir -> Occ 
                int pa = orbs_a[a].second.second;
                int pa2 = nalphapi_[ha];  
 
                nso = nsopi_[ha];
                nmo = nmopi_[ha];

                Ca = Ca_->pointer(ha);
                Ct = new double[nso]; 
                eps = epsilon_a_->pointer(ha);
     
                // Swap eigvals 
                epst = eps[pa];
                eps[pa] = eps[pa2];
                eps[pa2] = epst; 

                // Swap eigvecs
                C_DCOPY(nso, &Ca[0][pa], nmo, Ct, 1);
                C_DCOPY(nso, &Ca[0][pa2], nmo, &Ca[0][pa], nmo);
                C_DCOPY(nso, Ct, 1, &Ca[0][pa2], nmo);
            
                delete[] Ct;
                
                // Redo indexing
                nalphapi_[ha]++;
                nbetapi_[ha]++;
                doccpi_[ha]++; 

                fprintf(outfile, "   %8s: %4d%-4s -> %4d%-4s \n", "AB -> AB", pi + 1, ct.gamma(hi).symbol(), pa + 1, ct.gamma(ha).symbol()); 
            }

        } else if (options_.get_str("REFERENCE") == "UHF" || options_.get_str("REFERENCE") == "UKS") {

            if (si && sa) {

                // Alpha - alpha
                int hi = orbs_a[i].second.first;
                int ha = orbs_a[a].second.first;
 
                if (hi == ha) {
                    // Same irrep
                    int pi = orbs_a[i].second.second;
                    int pa = orbs_a[a].second.second;
  
                    int nso = nsopi_[hi];
                    int nmo = nmopi_[hi];

                    double** Ca = Ca_->pointer(hi);
                    double* Ct = new double[nso]; 
                    double* eps = epsilon_a_->pointer(hi);
                    double  epst;
     
                    // Swap eigvals 
                    epst = eps[pi];
                    eps[pi] = eps[pa];
                    eps[pa] = epst; 

                    // Swap eigvecs
                    C_DCOPY(nso, &Ca[0][pi], nmo, Ct, 1);
                    C_DCOPY(nso, &Ca[0][pa], nmo, &Ca[0][pi], nmo);
                    C_DCOPY(nso, Ct, 1, &Ca[0][pa], nmo);
                
                    delete[] Ct;

                    fprintf(outfile, "   %8s: %4d%-4s -> %4d%-4s \n", "A  -> A ", pi + 1, ct.gamma(hi).symbol(), pa + 1, ct.gamma(ha).symbol()); 
                } else {
                    // Different irrep
                    // Occ -> Vir
                    int pi = orbs_a[i].second.second;
                    int pi2 = nalphapi_[hi] - 1;  
 
                    int nso = nsopi_[hi];
                    int nmo = nmopi_[hi];

                    double** Ca = Ca_->pointer(hi);
                    double* Ct = new double[nso]; 
                    double* eps = epsilon_a_->pointer(hi);
                    double  epst;
     
                    // Swap eigvals 
                    epst = eps[pi];
                    eps[pi] = eps[pi2];
                    eps[pi2] = epst; 

                    // Swap eigvecs
                    C_DCOPY(nso, &Ca[0][pi], nmo, Ct, 1);
                    C_DCOPY(nso, &Ca[0][pi2], nmo, &Ca[0][pi], nmo);
                    C_DCOPY(nso, Ct, 1, &Ca[0][pi2], nmo);
                
                    delete[] Ct;
                    
                    // Redo indexing
                    nalphapi_[hi]--;

                    // Vir -> Occ 
                    int pa = orbs_a[a].second.second;
                    int pa2 = nalphapi_[ha];  
 
                    nso = nsopi_[ha];
                    nmo = nmopi_[ha];

                    Ca = Ca_->pointer(ha);
                    Ct = new double[nso]; 
                    eps = epsilon_a_->pointer(ha);
     
                    // Swap eigvals 
                    epst = eps[pa];
                    eps[pa] = eps[pa2];
                    eps[pa2] = epst; 

                    // Swap eigvecs
                    C_DCOPY(nso, &Ca[0][pa], nmo, Ct, 1);
                    C_DCOPY(nso, &Ca[0][pa2], nmo, &Ca[0][pa], nmo);
                    C_DCOPY(nso, Ct, 1, &Ca[0][pa2], nmo);
                
                    delete[] Ct;
                    
                    // Redo indexing
                    nalphapi_[ha]++;

                    fprintf(outfile, "   %8s: %4d%-4s -> %4d%-4s \n", "A  -> A ", pi + 1, ct.gamma(hi).symbol(), pa + 1, ct.gamma(ha).symbol()); 
                }
            } else if (!si && !sa) {
                // Beta->Beta
                int hi = orbs_b[i].second.first;
                int ha = orbs_b[a].second.first;
 
                if (hi == ha) {
                    // Same irrep
                    int pi = orbs_b[i].second.second;
                    int pa = orbs_b[a].second.second;
  
                    int nso = nsopi_[hi];
                    int nmo = nmopi_[hi];

                    double** Ca = Cb_->pointer(hi);
                    double* Ct = new double[nso]; 
                    double* eps = epsilon_b_->pointer(hi);
                    double  epst;
     
                    // Swap eigvals 
                    epst = eps[pi];
                    eps[pi] = eps[pa];
                    eps[pa] = epst; 

                    // Swap eigvecs
                    C_DCOPY(nso, &Ca[0][pi], nmo, Ct, 1);
                    C_DCOPY(nso, &Ca[0][pa], nmo, &Ca[0][pi], nmo);
                    C_DCOPY(nso, Ct, 1, &Ca[0][pa], nmo);
                
                    delete[] Ct;

                    fprintf(outfile, "   %8s: %4d%-4s -> %4d%-4s \n", "B  -> B ", pi + 1, ct.gamma(hi).symbol(), pa + 1, ct.gamma(ha).symbol()); 
                } else {
                    // Different irrep
                    // Occ -> Vir
                    int pi = orbs_b[i].second.second;
                    int pi2 = nbetapi_[hi] - 1;  
 
                    int nso = nsopi_[hi];
                    int nmo = nmopi_[hi];

                    double** Ca = Cb_->pointer(hi);
                    double* Ct = new double[nso]; 
                    double* eps = epsilon_b_->pointer(hi);
                    double  epst;
     
                    // Swap eigvals 
                    epst = eps[pi];
                    eps[pi] = eps[pi2];
                    eps[pi2] = epst; 

                    // Swap eigvecs
                    C_DCOPY(nso, &Ca[0][pi], nmo, Ct, 1);
                    C_DCOPY(nso, &Ca[0][pi2], nmo, &Ca[0][pi], nmo);
                    C_DCOPY(nso, Ct, 1, &Ca[0][pi2], nmo);
                
                    delete[] Ct;
                    
                    // Redo indexing
                    nbetapi_[hi]--;

                    // Vir -> Occ 
                    int pa = orbs_b[a].second.second;
                    int pa2 = nbetapi_[ha];  
 
                    nso = nsopi_[ha];
                    nmo = nmopi_[ha];

                    Ca = Cb_->pointer(ha);
                    Ct = new double[nso]; 
                    eps = epsilon_b_->pointer(ha);
     
                    // Swap eigvals 
                    epst = eps[pa];
                    eps[pa] = eps[pa2];
                    eps[pa2] = epst; 

                    // Swap eigvecs
                    C_DCOPY(nso, &Ca[0][pa], nmo, Ct, 1);
                    C_DCOPY(nso, &Ca[0][pa2], nmo, &Ca[0][pa], nmo);
                    C_DCOPY(nso, Ct, 1, &Ca[0][pa2], nmo);
                
                    delete[] Ct;
                    
                    // Redo indexing
                    nbetapi_[ha]++;

                    fprintf(outfile, "   %8s: %4d%-4s -> %4d%-4s \n", "B  -> B ", pi + 1, ct.gamma(hi).symbol(), pa + 1, ct.gamma(ha).symbol()); 
                }
            } else if (!si && sa) {
                // Beta->Alpha
                int hi = orbs_b[i].second.first;
                int ha = orbs_a[a].second.first;
 
                // Different irrep
                // Occ -> Vir
                int pi = orbs_b[i].second.second;
                int pi2 = nbetapi_[hi] - 1;  
 
                int nso = nsopi_[hi];
                int nmo = nmopi_[hi];

                double** Ca = Cb_->pointer(hi);
                double* Ct = new double[nso]; 
                double* eps = epsilon_b_->pointer(hi);
                double  epst;
     
                // Swap eigvals 
                epst = eps[pi];
                eps[pi] = eps[pi2];
                eps[pi2] = epst; 

                // Swap eigvecs
                C_DCOPY(nso, &Ca[0][pi], nmo, Ct, 1);
                C_DCOPY(nso, &Ca[0][pi2], nmo, &Ca[0][pi], nmo);
                C_DCOPY(nso, Ct, 1, &Ca[0][pi2], nmo);
                
                delete[] Ct;
                
                // Redo indexing
                nbetapi_[hi]--;
                nbeta_--;

                // Vir -> Occ 
                int pa = orbs_a[a].second.second;
                int pa2 = nalphapi_[ha];  
 
                nso = nsopi_[ha];
                nmo = nmopi_[ha];

                Ca = Ca_->pointer(ha);
                Ct = new double[nso]; 
                eps = epsilon_a_->pointer(ha);
     
                // Swap eigvals 
                epst = eps[pa];
                eps[pa] = eps[pa2];
                eps[pa2] = epst; 

                // Swap eigvecs
                C_DCOPY(nso, &Ca[0][pa], nmo, Ct, 1);
                C_DCOPY(nso, &Ca[0][pa2], nmo, &Ca[0][pa], nmo);
                C_DCOPY(nso, Ct, 1, &Ca[0][pa2], nmo);
                
                delete[] Ct;
                
                // Redo indexing
                nalphapi_[ha]++;
                nalpha_++;

                fprintf(outfile, "   %8s: %4d%-4s -> %4d%-4s \n", "B  -> A ", pi + 1, ct.gamma(hi).symbol(), pa + 1, ct.gamma(ha).symbol()); 
            } else if (si && !sa) {
                // Alpha->Beta
                int hi = orbs_a[i].second.first;
                int ha = orbs_b[a].second.first;
 
                // Different irrep
                // Occ -> Vir
                int pi = orbs_a[i].second.second;
                int pi2 = nalphapi_[hi] - 1;  
 
                int nso = nsopi_[hi];
                int nmo = nmopi_[hi];

                double** Ca = Ca_->pointer(hi);
                double* Ct = new double[nso]; 
                double* eps = epsilon_a_->pointer(hi);
                double  epst;
     
                // Swap eigvals 
                epst = eps[pi];
                eps[pi] = eps[pi2];
                eps[pi2] = epst; 

                // Swap eigvecs
                C_DCOPY(nso, &Ca[0][pi], nmo, Ct, 1);
                C_DCOPY(nso, &Ca[0][pi2], nmo, &Ca[0][pi], nmo);
                C_DCOPY(nso, Ct, 1, &Ca[0][pi2], nmo);
                
                delete[] Ct;
                
                // Redo indexing
                nalphapi_[hi]--;
                nalpha_--;

                // Vir -> Occ 
                int pa = orbs_b[a].second.second;
                int pa2 = nbetapi_[ha];  
 
                nso = nsopi_[ha];
                nmo = nmopi_[ha];

                Ca = Cb_->pointer(ha);
                Ct = new double[nso]; 
                eps = epsilon_b_->pointer(ha);
     
                // Swap eigvals 
                epst = eps[pa];
                eps[pa] = eps[pa2];
                eps[pa2] = epst; 

                // Swap eigvecs
                C_DCOPY(nso, &Ca[0][pa], nmo, Ct, 1);
                C_DCOPY(nso, &Ca[0][pa2], nmo, &Ca[0][pa], nmo);
                C_DCOPY(nso, Ct, 1, &Ca[0][pa2], nmo);
                
                delete[] Ct;
                
                // Redo indexing
                nbetapi_[ha]++;
                nbeta_++;

                fprintf(outfile, "   %8s: %4d%-4s -> %4d%-4s \n", "A  -> B ", pi + 1, ct.gamma(hi).symbol(), pa + 1, ct.gamma(ha).symbol()); 
            }
            if (nalpha_ < nbeta_) throw PSIEXCEPTION("PSI::MOM_start: Nbeta ends up being less than Nalpha, this is not supported");
           
            // Fix doccpi/soccpi. 
            for (int h = 0; h < nirrep_; h++) {
                std::vector<std::pair<double,std::pair<int, bool> > > alphas;    
                std::vector<std::pair<double,std::pair<int, bool> > > betas;

                for (int i = 0; i < nalphapi_[h]; i++) {
                    alphas.push_back(make_pair(epsilon_a_->get(h,i), make_pair(i, true)));                                
                }
                for (int i = 0; i < nbetapi_[h]; i++) {
                    betas.push_back(make_pair(epsilon_b_->get(h,i), make_pair(i, true)));                                
                }
                for (int i = nalphapi_[h]; i < nmopi_[h]; i++) {
                    alphas.push_back(make_pair(epsilon_a_->get(h,i), make_pair(i, false)));                                
                }
                for (int i = nbetapi_[h]; i < nmopi_[h]; i++) {
                    betas.push_back(make_pair(epsilon_b_->get(h,i), make_pair(i, false)));                                
                }
                sort(alphas.begin(),alphas.end());       
                sort(betas.begin(),betas.end());       

                doccpi_[h] = 0;
                soccpi_[h] = 0;

                for (int i = 0; i < nmopi_[h]; i++) {
                    bool alpha_occ = alphas[i].second.second;
                    bool beta_occ = betas[i].second.second;
                    if (alpha_occ && beta_occ)
                        doccpi_[h]++;
                    else if (alpha_occ || beta_occ) // Careful, could be beta occ
                        soccpi_[h]++;
                }
            }

        } else if (options_.get_str("REFERENCE") == "ROHF") {
            throw PSIEXCEPTION("SCF::MOM_start: MOM excited states are not implemented for ROHF");
        }
    }
    Ca_old_->copy(Ca_);
    Cb_old_->copy(Cb_);

    fprintf(outfile, "\n                        Total Energy        Delta E      Density RMS\n\n");
}
void HF::MOM()
{
    // Go MOM go!
    // Alpha
    for (int h = 0; h < nirrep_; h++) {
    
        // Indexing
        int nso = nsopi_[h];
        int nmo = nmopi_[h];
        int nalpha = nalphapi_[h];
    
        if (nso == 0 || nmo == 0 || nalpha == 0) continue;
        
        double** Cold = Ca_old_->pointer(h);
        double** Cnew = Ca_->pointer(h);
        double*  eps  = epsilon_a_->pointer(h);
        double** S = S_->pointer(h);

        double* c = new double[nso];
        double* d = new double[nso];
        double* p = new double[nmo]; 
        
        memset(static_cast<void*>(c), '\0', sizeof(double)*nso);
        
        for (int i = 0; i < nalpha; i++)
            C_DAXPY(nso,1.0,&Cold[0][i],nmo,c,1);

        C_DGEMV('N',nso,nso,1.0,S[0],nso,c,1,0.0,d,1);
        C_DGEMV('T',nso,nmo,1.0,Cnew[0],nmo,d,1,0.0,p,1);

        //fprintf(outfile,"  P_a:\n");
        //for (int a = 0; a < nmo; a++)
        //    fprintf(outfile,"   a = %3d: %14.10f\n", a + 1, p[a]);

        // Find the largest contributions
        std::vector<std::pair<double, int> > pvec;
        pvec.resize(nmo);
        for (int a = 0; a < nmo; a++)
            pvec[a] = make_pair(fabs(p[a]), a);
        sort(pvec.begin(),pvec.end(), greater<std::pair<double, int> >());

        //fprintf(outfile,"  P_a sorted:\n");
        //for (int a = 0; a < nmo; a++)
        //    fprintf(outfile,"   a = %3d: Index = %3d, %14.10f\n", a + 1, pvec[a].second, pvec[a].first);

        // Now order the mos in each group
        std::vector<std::pair<double, int> > occvec;
        occvec.resize(nalpha);
        for (int a = 0; a < nalpha; a++)
            occvec[a] = make_pair(eps[pvec[a].second], pvec[a].second);
        sort(occvec.begin(),occvec.end());

        //fprintf(outfile,"  P_a_occ sorted:\n");
        //for (int a = 0; a < nalpha; a++)
        //    fprintf(outfile,"   a = %3d: Index = %3d, %14.10f\n", a + 1, occvec[a].second, occvec[a].first);

        std::vector<std::pair<double, int> > virvec;
        virvec.resize(nmo - nalpha);
        for (int a = 0; a < nmo - nalpha; a++)
            virvec[a] = make_pair(eps[pvec[a + nalpha].second], pvec[a + nalpha].second);
        sort(virvec.begin(),virvec.end());

        //fprintf(outfile,"  P_a_vir sorted:\n");
        //for (int a = 0; a < nmo - nalpha; a++)
        //    fprintf(outfile,"   a = %3d: Index = %3d, %14.10f\n", a + 1, virvec[a].second, virvec[a].first);

        double** Ct = block_matrix(nso,nmo);

        // Use Cold and p as a buffer
        memcpy(static_cast<void*>(Ct[0]),static_cast<void*>(Cnew[0]),sizeof(double)*nso*nmo);
        memcpy(static_cast<void*>(p),      static_cast<void*>(eps),    sizeof(double)*nmo);

        for (int a = 0; a < nalpha; a++) {
            eps[a] = occvec[a].first;
            C_DCOPY(nso,&Ct[0][occvec[a].second],nmo,&Cnew[0][a],nmo);    
        } 

        for (int a = 0; a < nmo - nalpha; a++) {
            eps[a + nalpha] = virvec[a].first;
            C_DCOPY(nso,&Ct[0][virvec[a].second],nmo,&Cnew[0][a + nalpha],nmo);    
        } 

        free_block(Ct);

        delete[] c;
        delete[] d;
        delete[] p;
    }   
   
    // Beta
    if (restricted()) return;
  
    for (int h = 0; h < nirrep_; h++) {
    
        // Indexing
        int nso = nsopi_[h];
        int nmo = nmopi_[h];
        int nbeta = nbetapi_[h];
    
        if (nso == 0 || nmo == 0 || nbeta == 0) continue;
        
        double** Cold = Cb_old_->pointer(h);
        double** Cnew = Cb_->pointer(h);
        double*  eps  = epsilon_b_->pointer(h);
        double** S = S_->pointer(h);

        double* c = new double[nso];
        double* d = new double[nso];
        double* p = new double[nmo]; 
        
        memset(static_cast<void*>(c), '\0', sizeof(double)*nso);
        
        for (int i = 0; i < nbeta; i++)
            C_DAXPY(nso,1.0,&Cold[0][i],nmo,c,1);

        C_DGEMV('N',nso,nso,1.0,S[0],nso,c,1,0.0,d,1);
        C_DGEMV('T',nso,nmo,1.0,Cnew[0],nmo,d,1,0.0,p,1);

        //fprintf(outfile,"  P_a:\n");
        //for (int a = 0; a < nmo; a++)
        //    fprintf(outfile,"   a = %3d: %14.10f\n", a + 1, p[a]);

        // Find the largest contributions
        std::vector<std::pair<double, int> > pvec;
        pvec.resize(nmo);
        for (int a = 0; a < nmo; a++)
            pvec[a] = make_pair(fabs(p[a]), a);
        sort(pvec.begin(),pvec.end(), greater<std::pair<double, int> >());

        //fprintf(outfile,"  P_a sorted:\n");
        //for (int a = 0; a < nmo; a++)
        //    fprintf(outfile,"   a = %3d: Index = %3d, %14.10f\n", a + 1, pvec[a].second, pvec[a].first);

        // Now order the mos in each group
        std::vector<std::pair<double, int> > occvec;
        occvec.resize(nbeta);
        for (int a = 0; a < nbeta; a++)
            occvec[a] = make_pair(eps[pvec[a].second], pvec[a].second);
        sort(occvec.begin(),occvec.end());

        //fprintf(outfile,"  P_a_occ sorted:\n");
        //for (int a = 0; a < nbeta; a++)
        //    fprintf(outfile,"   a = %3d: Index = %3d, %14.10f\n", a + 1, occvec[a].second, occvec[a].first);

        std::vector<std::pair<double, int> > virvec;
        virvec.resize(nmo - nbeta);
        for (int a = 0; a < nmo - nbeta; a++)
            virvec[a] = make_pair(eps[pvec[a + nbeta].second], pvec[a + nbeta].second);
        sort(virvec.begin(),virvec.end());

        //fprintf(outfile,"  P_a_vir sorted:\n");
        //for (int a = 0; a < nmo - nalpha; a++)
        //    fprintf(outfile,"   a = %3d: Index = %3d, %14.10f\n", a + 1, virvec[a].second, virvec[a].first);

        double** Ct = block_matrix(nso,nmo);
        
        // Use Cold and p as a buffer
        memcpy(static_cast<void*>(Ct[0]),static_cast<void*>(Cnew[0]),sizeof(double)*nso*nmo);
        memcpy(static_cast<void*>(p),      static_cast<void*>(eps),    sizeof(double)*nmo);

        for (int a = 0; a < nbeta; a++) {
            eps[a] = occvec[a].first;
            C_DCOPY(nso,&Ct[0][occvec[a].second],nmo,&Cnew[0][a],nmo);    
        } 

        for (int a = 0; a < nmo - nbeta; a++) {
            eps[a + nbeta] = virvec[a].first;
            C_DCOPY(nso,&Ct[0][virvec[a].second],nmo,&Cnew[0][a + nbeta],nmo);    
        } 

        free_block(Ct);

        delete[] c;
        delete[] d;
        delete[] p;
    }   
    Cb_old_->copy(Cb_);
}

}}

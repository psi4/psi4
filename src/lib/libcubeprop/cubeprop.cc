/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

#include <boost/filesystem.hpp>

#include <psi4-dec.h>

#include <libpsi4util/libpsi4util.h>
#include <libmints/mints.h>

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include "cubeprop.h"
#include "csg.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace psi;
using namespace boost;
using namespace std;

namespace psi {

CubeProperties::CubeProperties() :
    options_(Process::environment.options)
{
    boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
    basisset_ = wfn->basisset();

    Ca_ = wfn->Ca_subset("AO", "ALL");
    Da_ = wfn->Da_subset("AO");

    if (wfn->same_a_b_orbs()) { 
        Cb_ = Ca_;
    } else { 
        Cb_ = wfn->Cb_subset("AO", "ALL");
    }

    if (wfn->same_a_b_dens()) { 
        Db_ = Da_;
    } else { 
        Db_ = wfn->Db_subset("AO");
    }

    int nirrep = wfn->nirrep();
    Dimension nmopi = wfn->nmopi();
    // Gather orbital information
    for (int h = 0; h < nirrep; h++) {
        for (int i = 0; i < (int)nmopi[h]; i++) {
            info_a_.push_back(boost::tuple<double,int,int>(wfn->epsilon_a()->get(h,i),i,h));
        }
    }
    std::sort(info_a_.begin(), info_a_.end(), std::less<boost::tuple<double,int,int> >()); // Sort as in wfn
    for (int h = 0; h < nirrep; h++) {
        for (int i = 0; i < (int)nmopi[h]; i++) {
            info_b_.push_back(boost::tuple<double,int,int>(wfn->epsilon_b()->get(h,i),i,h));
        }
    }
    std::sort(info_b_.begin(), info_b_.end(), std::less<boost::tuple<double,int,int> >()); // Sort as in wfn

    common_init();
}
CubeProperties::~CubeProperties()
{
}
void CubeProperties::common_init()
{
    grid_ = boost::shared_ptr<CubicScalarGrid>(new CubicScalarGrid(basisset_));
    grid_->set_filepath(options_.get_str("CUBEPROP_FILEPATH"));
}
void CubeProperties::print_header()
{
    outfile->Printf( "  ==> One Electron Grid Properties (v2.0) <==\n\n");
    grid_->print_header();
    outfile->Flush();
}
void CubeProperties::compute_properties()
{
    print_header();

    std::string filepath = options_.get_str("CUBEPROP_FILEPATH");
    std::stringstream ss;
    ss << filepath << "/" << "geom.xyz";

    // Is filepath a valid directory?
    boost::filesystem::path data_dir(filepath);
    if(not boost::filesystem::is_directory(data_dir)){
        printf("Filepath \"%s\" is not valid.  Please create this directory.\n",filepath.c_str());
        outfile->Printf("Filepath \"%s\" is not valid.  Please create this directory.\n",filepath.c_str());
        outfile->Flush();
        exit(Failure);
    }

    basisset_->molecule()->save_xyz_file(ss.str());

    for (size_t ind = 0; ind < options_["CUBEPROP_TASKS"].size(); ind++) {
        std::string task = options_["CUBEPROP_TASKS"][ind].to_string();

        if (task == "DENSITY") {
            boost::shared_ptr<Matrix> Dt(Da_->clone());        
            boost::shared_ptr<Matrix> Ds(Da_->clone());        
            Dt->copy(Da_);
            Ds->copy(Da_);
            Dt->add(Db_);
            Ds->subtract(Db_);
            compute_density(Dt, "Dt");
            compute_density(Ds, "Ds");
            compute_density(Da_, "Da");
            compute_density(Db_, "Db");
        } else if (task == "ESP") {
            boost::shared_ptr<Matrix> Dt(Da_->clone());        
            Dt->copy(Da_);
            Dt->add(Db_);
            compute_esp(Dt);
        } else if (task == "ORBITALS") {
            std::vector<int> indsa0;
            std::vector<int> indsb0;

            if (options_["CUBEPROP_ORBITALS"].size() == 0) {
                for (int ind = 0; ind < Ca_->colspi()[0]; ind++) {
                    indsa0.push_back(ind);
                } 
                for (int ind = 0; ind < Cb_->colspi()[0]; ind++) {
                    indsb0.push_back(ind);
                } 
            } else {
                for (size_t ind = 0; ind < options_["CUBEPROP_ORBITALS"].size(); ind++) {
                    int val = options_["CUBEPROP_ORBITALS"][ind].to_integer();
                    if (val > 0) {
                        indsa0.push_back(abs(val) - 1); 
                    } else {
                        indsb0.push_back(abs(val) - 1); 
                    }
                }
            }
            std::vector<string> labelsa;
            std::vector<string> labelsb;
            CharacterTable ct = basisset_->molecule()->point_group()->char_table();
            for (size_t ind = 0; ind < indsa0.size(); ++ind){
                int i = get<1>(info_a_[indsa0[ind]]);
                int h = get<2>(info_a_[indsa0[ind]]);
                labelsa.push_back(to_string(i + 1) + "-" + ct.gamma(h).symbol());
            }
            for (size_t ind = 0; ind < indsb0.size(); ++ind){
                int i = get<1>(info_b_[indsb0[ind]]);
                int h = get<2>(info_b_[indsb0[ind]]);
                labelsb.push_back(to_string(i + 1) + "-" + ct.gamma(h).symbol());
            }
            if (indsa0.size()) compute_orbitals(Ca_, indsa0,labelsa, "Psi_a");
            if (indsb0.size()) compute_orbitals(Cb_, indsb0,labelsb, "Psi_b");
        } else if (task == "BASIS_FUNCTIONS") {
            std::vector<int> inds0;
            if (options_["CUBEPROP_BASIS_FUNCTIONS"].size() == 0) {
                for (int ind = 0; ind < basisset_->nbf(); ind++) {
                    inds0.push_back(ind);
                } 
            } else {
                for (size_t ind = 0; ind < options_["CUBEPROP_BASIS_FUNCTIONS"].size(); ind++) {
                    inds0.push_back(options_["CUBEPROP_BASIS_FUNCTIONS"][ind].to_integer() - 1);
                }
            }
            compute_basis_functions(inds0, "Phi");
        } else if (task == "LOL") {
            compute_LOL(Da_, "LOLa");
            compute_LOL(Db_, "LOLb");
        } else if (task == "ELF") {
            compute_ELF(Da_, "ELFa");
            compute_ELF(Db_, "ELFb");
        } else {
            throw PSIEXCEPTION(task + "is an unrecognized PROPERTY_TASKS value");
        }
    }
}
void CubeProperties::compute_density(boost::shared_ptr<Matrix> D, const std::string& key)
{
    grid_->compute_density(D, key);
}
void CubeProperties::compute_esp(boost::shared_ptr<Matrix> Dt)
{
    grid_->compute_density(Dt, "Dt");
    grid_->compute_esp(Dt, "ESP"); 
}
void CubeProperties::compute_orbitals(boost::shared_ptr<Matrix> C, const std::vector<int>& indices, const std::vector<std::string>& labels, const std::string& key)
{
    grid_->compute_orbitals(C, indices, labels, key);
}
void CubeProperties::compute_basis_functions(const std::vector<int>& indices, const std::string& key)
{
    grid_->compute_basis_functions(indices, key);
}
void CubeProperties::compute_LOL(boost::shared_ptr<Matrix> D, const std::string& key)
{
    grid_->compute_LOL(D, key);
}
void CubeProperties::compute_ELF(boost::shared_ptr<Matrix> D, const std::string& key)
{
    grid_->compute_ELF(D, key);
}

}


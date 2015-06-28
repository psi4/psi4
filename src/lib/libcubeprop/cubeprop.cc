#include <boost/filesystem.hpp>
#include "prop.h"
#include "csg.h"
#include <libmints/mints.h>
#include <psi4-dec.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace psi;
using namespace boost;
using namespace std;

namespace psi {

Properties::Properties(boost::shared_ptr<Wavefunction> wfn) :
    options_(Process::environment.options)
{
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

    common_init();
}
Properties::~Properties()
{
}
void Properties::common_init()
{
    grid_ = boost::shared_ptr<CubicScalarGrid>(new CubicScalarGrid(basisset_));
    grid_->set_filepath(options_.get_str("PROPERTY_FILEPATH"));
}
void Properties::print_header()
{
    outfile->Printf( "  ==> One Electron Grid Properties (v2.0) <==\n\n");
    grid_->print_header();
    outfile->Flush();
}
void Properties::compute_properties()
{
    std::string filepath = options_.get_str("PROPERTY_FILEPATH");
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

    for (int ind = 0; ind < options_["PROPERTY_TASKS"].size(); ind++) {
        std::string task = options_["PROPERTY_TASKS"][ind].to_string();

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
            if (options_["PROPERTY_ORBITALS"].size() == 0) {
                for (int ind = 0; ind < Ca_->colspi()[0]; ind++) {
                    indsa0.push_back(ind);
                } 
                for (int ind = 0; ind < Cb_->colspi()[0]; ind++) {
                    indsb0.push_back(ind);
                } 
            } else {
                for (int ind = 0; ind < options_["PROPERTY_ORBITALS"].size(); ind++) {
                    int val = options_["PROPERTY_ORBITALS"][ind].to_integer();
                    if (val > 0) {
                        indsa0.push_back(abs(val) - 1); 
                    } else {
                        indsb0.push_back(abs(val) - 1); 
                    }
                }
            }
            if (indsa0.size()) compute_orbitals(Ca_, indsa0, "Psi_a");
            if (indsb0.size()) compute_orbitals(Cb_, indsb0, "Psi_b");
        } else if (task == "BASIS_FUNCTIONS") {
            std::vector<int> inds0;
            if (options_["PROPERTY_BASIS_FUNCTIONS"].size() == 0) {
                for (int ind = 0; ind < basisset_->nbf(); ind++) {
                    inds0.push_back(ind);
                } 
            } else {
                for (int ind = 0; ind < options_["PROPERTY_BASIS_FUNCTIONS"].size(); ind++) {
                    inds0.push_back(options_["PROPERTY_BASIS_FUNCTIONS"][ind].to_integer() - 1);
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
            throw PSIEXCEPTION("Unrecognized PROPERTY_TASKS value");
        }
    }
}
void Properties::compute_density(boost::shared_ptr<Matrix> D, const std::string& key)
{
    grid_->compute_density(D, key);
}
void Properties::compute_esp(boost::shared_ptr<Matrix> Dt)
{
    grid_->compute_density(Dt, "Dt");
    grid_->compute_esp(Dt, "ESP"); 
}
void Properties::compute_orbitals(boost::shared_ptr<Matrix> C, const std::vector<int>& indices, const std::string& key)
{
    grid_->compute_orbitals(C, indices, key);
}
void Properties::compute_basis_functions(const std::vector<int>& indices, const std::string& key)
{
    grid_->compute_basis_functions(indices, key);
}
void Properties::compute_LOL(boost::shared_ptr<Matrix> D, const std::string& key)
{
    grid_->compute_LOL(D, key);
}
void Properties::compute_ELF(boost::shared_ptr<Matrix> D, const std::string& key)
{
    grid_->compute_ELF(D, key);
}

}


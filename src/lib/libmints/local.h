#ifndef _psi_src_lib_libmints_local_h_
#define _psi_src_lib_libmints_local_h_

#include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>

namespace psi {

class Matrix;
class BasisSet;
class OrbitalDomain;

class Local {

protected:
    /// Print flag, defaults to 0
    int print_;
    /// Debug flag, defaults to 0
    int debug_;

    /// Basisset corresponding to C
    boost::shared_ptr<BasisSet> basisset_;

    /// Reference to externally provided C in USO basis, active part only
    boost::shared_ptr<Matrix> C_USO_;
    /// Current C in C1 basis
    boost::shared_ptr<Matrix> C_AO_;
    /// Localized C in C1 basis
    boost::shared_ptr<Matrix> L_AO_;
    /// AO2USO transformer, computed from basisset
    boost::shared_ptr<Matrix> AO2USO_;

    /// Orbital domains, computed/updated by compute_domains
    std::vector<boost::shared_ptr<OrbitalDomain> > domains_;

    /// Set up USO to AO and inital matrices
    void common_init();

    /// Localize via pivoted-Cholesky
    boost::shared_ptr<Matrix> localize_cholesky(double conv);
    /// Localize via Pipek-Mezey
    boost::shared_ptr<Matrix> localize_pipek_mezey(double conv);
    /// Localize via Boys 
    boost::shared_ptr<Matrix> localize_boys(double conv);
    /// Localize via ER 
    boost::shared_ptr<Matrix> localize_er(double conv);

public: 
    /// Constructor, builds USO->AO info, holds a reference to C 
    Local(boost::shared_ptr<BasisSet> basisset, boost::shared_ptr<Matrix> C);
    virtual ~Local();    

    void set_print(int print) { print_ = print; }
    void set_debug(int debug) { debug_ = debug; }

    boost::shared_ptr<Matrix> C_USO();
    boost::shared_ptr<Matrix> C_AO();
    boost::shared_ptr<Matrix> L_AO();
    boost::shared_ptr<Matrix> AO2USO();
    const std::vector<boost::shared_ptr<OrbitalDomain> >& domains() const {return domains_; } 

    boost::shared_ptr<Matrix> localize(const std::string& algorithm, double conv = 1.0E-10);
    const std::vector<boost::shared_ptr<OrbitalDomain> >& compute_domains(double Qcutoff = 0.02, double Rcutoff = 0.0); 
};

class OrbitalDomain {
    protected:
        std::vector<int>  functions_local_to_global_;
        std::map<int,int> functions_global_to_local_;
        std::vector<int>  shells_local_to_global_;
        std::map<int,int> shells_global_to_local_;
    public:
        OrbitalDomain();
        virtual ~OrbitalDomain();

        const std::vector<int>& functions_local_to_global() const { return functions_local_to_global_; }       
        const std::map<int,int>& functions_global_to_local() const { return functions_global_to_local_; }       
        const std::vector<int>& shells_local_to_global() const { return shells_local_to_global_; }       
        const std::map<int,int>& shells_global_to_local() const { return shells_global_to_local_; }       

    friend class Local;
};

} //Namespace psi

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
    /// Auxiliary Basisset 
    boost::shared_ptr<BasisSet> auxiliary_;
    /// AO2USO transformer, computed from basisset
    boost::shared_ptr<Matrix> AO2USO_;
    /// Inverse square root of the S matrix
    boost::shared_ptr<Matrix> X_; 

    /// Reference to externally provided C in USO basis, active part only
    boost::shared_ptr<Matrix> C_USO_;
    /// Current C in C1 basis
    boost::shared_ptr<Matrix> C_AO_;
    /// Localized C in C1 basis
    boost::shared_ptr<Matrix> L_AO_;

    /// gross Lowdin charges (nmo x natom)
    boost::shared_ptr<Matrix> Q_;
    /// Orbital domains, computed/updated by compute_domains
    std::vector<boost::shared_ptr<OrbitalDomain> > domains_;
    /// Orbital domains, auxiliary, computed/updated by compute_domains
    std::vector<boost::shared_ptr<OrbitalDomain> > auxiliary_domains_;

    /// Set up USO to AO and inital matrices
    void common_init();

    /// Localize via pivoted-Cholesky
    void localize_cholesky(double conv);
    /// Localize via Pipek-Mezey
    void localize_pipek_mezey(double conv);
    /// Localize via Boys 
    void localize_boys(double conv);
    /// Localize via ER 
    void localize_er(double conv);

    /// Compute gross Lowdin charges for coef matrix C (nmo x natom)
    boost::shared_ptr<Matrix> lowdin_charges(boost::shared_ptr<Matrix> C);

public: 
    /*!
    * Constructor, builds a Local object with a reference to basisset and C
    * \param basisset: the AO basis object
    * \param C: the C matrix corresponding the the basisset. This class
    * holds a reference to C, to use Local iteratively, you externally modify C,
    * then call Local::localize() 
    */
    Local(boost::shared_ptr<BasisSet> basisset, boost::shared_ptr<Matrix> C);
    /*!
    * Constructor, builds a Local object with a reference to basisset and C
    * \param basisset: the AO primary basis object
    * \param auxiliary: the AO auxiliary basis object
    * \param C: the C matrix corresponding the the basisset. This class
    * holds a reference to C, to use Local iteratively, you externally modify C,
    * then call Local::localize() 
    */
    Local(boost::shared_ptr<BasisSet> basisset, boost::shared_ptr<BasisSet> auxiliary, boost::shared_ptr<Matrix> C);
    virtual ~Local();    

    /// Set the print flag (default 0)
    void set_print(int print) { print_ = print; }
    /// Set the debug flag (default 0)
    void set_debug(int debug) { debug_ = debug; }

    /*!
    * Reference to the USO C matrix 
    * \return boost::shared_ptr<Matrix> 
    */
    boost::shared_ptr<Matrix> C_USO();
    /*!
    * Reference to the C1 C matrix 
    * \return boost::shared_ptr<Matrix> 
    */
    boost::shared_ptr<Matrix> C_AO();
    /*!
    * Reference to the C1 localized orbitals matrix 
    * \return boost::shared_ptr<Matrix> 
    */
    boost::shared_ptr<Matrix> L_AO();
    /*!
    * Reference to the AO2USO matrix 
    * \return boost::shared_ptr<Matrix> 
    */
    boost::shared_ptr<Matrix> AO2USO();
    /*!
    * The vector of orbital domains, as determined by compute_X_domains()
    * \return vector of boost::shared_ptr<OrbitalDomain>
    */
    const std::vector<boost::shared_ptr<OrbitalDomain> >& domains() const {return domains_; } 
    /*!
    * The vector of orbital auxiliary domains, as determined by compute_X_domains()
    * \return vector of boost::shared_ptr<OrbitalDomain>
    */
    const std::vector<boost::shared_ptr<OrbitalDomain> >& auxiliary_domains() const {return auxiliary_domains_; } 

    /*!
    * Localize the current C matrix (nondestructive) storing the result in L_USO (C1)
    * \param algorithm: upper-case string containing algorithm name, which is one of:
    *   CHOLESKY
    *   PIPEK_MEZEY
    *   BOYS
    *   ER
    * \param conv: double containing convergence criteria for iterative methods
    */
    void localize(const std::string& algorithm, double conv = 1.0E-10);
    /*!
    * Setup the domains according to the Boughton-Pulay algorithm
    * \param Qcutoff: residual of approximate span, default 0.02
    */
    void compute_boughton_pulay_domains(double Qcutoff = 0.02);
    /*!
    * Setup the domains according to the Polly algorithm
    * \param Qcutoff: Minimum Lowdin charge to add initially, default 0.05
    * \param Rcutoff: Minimum spatial extent outside Lowdin domains to add,
    *   in angstrom, default 3.0 
    * \param Qcheck: Required sum of Lowdin charges in domain, defaults to 0.98
    */
    void compute_polly_domains(double Qcutoff = 0.05, double Rcutoff = 3.0, double Qcheck = 0.98);
};

class OrbitalDomain {
    protected:
        std::set<int>     atoms_;
        std::vector<int>  atoms_local_to_global_;
        std::map<int,int> atoms_global_to_local_;
        std::vector<int>  functions_local_to_global_;
        std::map<int,int> functions_global_to_local_;
        std::vector<int>  shells_local_to_global_;
        std::map<int,int> shells_global_to_local_;
    public:
        OrbitalDomain();
        virtual ~OrbitalDomain();

        void print(FILE* = outfile, int label = 0);

        static boost::shared_ptr<OrbitalDomain> buildOrbitalDomain(boost::shared_ptr<BasisSet> basis, const std::set<int>& atoms);

        const std::set<int>& atoms() const { return atoms_; }
        const std::vector<int>& atoms_local_to_global() const { return atoms_local_to_global_; }       
        const std::map<int,int>& atoms_global_to_local() const { return atoms_global_to_local_; }       
        const std::vector<int>& functions_local_to_global() const { return functions_local_to_global_; }       
        const std::map<int,int>& functions_global_to_local() const { return functions_global_to_local_; }       
        const std::vector<int>& shells_local_to_global() const { return shells_local_to_global_; }       
        const std::map<int,int>& shells_global_to_local() const { return shells_global_to_local_; }       

};

} //Namespace psi

#endif

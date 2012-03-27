#ifndef LIBFOCK_DFT_H
#define LIBFOCK_DFT_H

#include <libmints/typedefs.h>
#include <vector>
#include <map>

namespace psi {

class Options;
class DFTGrid;
class PointFunctions;
class SuperFunctional;

// => BASE CLASS <= //

/**
 * Class VBase
 *
 * Class to compute KS-V matrices and
 * K-matrix-vector products
 **/

class VBase {

protected:
    /// Debug flag
    int debug_;
    /// Print flag
    int print_;
    /// Options object, used to build grid
    Options& options_;
    /// Basis set used in the integration
    boost::shared_ptr<BasisSet> primary_;
    /// Desired superfunctional kernal
    boost::shared_ptr<SuperFunctional> functional_;
    /// Point function computer (densities, gammas, basis values)
    boost::shared_ptr<PointFunctions> properties_;
    /// Integration grid, built by KSPotential
    boost::shared_ptr<DFTGrid> grid_;
    /// Quadrature values obtained during integration 
    std::map<std::string, double> quad_values_;

    /// AO2USO matrix (if not C1)
    SharedMatrix AO2USO_;

    /// Vector of V matrices (built by form_D) 
    std::vector<SharedMatrix> V_;
    /// Vector of C1 V matrices (built by USO2AO)
    std::vector<SharedMatrix> V_AO_;

    /// Vector of occupied C matrices (used for D and KE density)
    std::vector<SharedMatrix> C_;
    /// Vector of D matrices (built by form_D) 
    std::vector<SharedMatrix> D_;
    /// Vector of C1 C matrices (built by USO2AO)
    std::vector<SharedMatrix> C_AO_;
    /// Vector of C1 D matrices (built by USO2AO)
    std::vector<SharedMatrix> D_AO_;

    /// Vector of Caocc matrices (TDDFT)
    std::vector<SharedMatrix> Caocc_;
    /// Vector of Cavir matrices (TDDFT)
    std::vector<SharedMatrix> Cavir_;
    /// Vector of Perturbation matrices (TDDFT, ia) 
    std::vector<SharedMatrix> P_;
    /// Vector of Perturbation matrices (TDDFT, SO) 
    std::vector<SharedMatrix> P_SO_;
    /// Vector of Perturbation matrices (TDDFT, AO) 
    std::vector<SharedMatrix> P_AO_;

    virtual void compute_D();
    virtual void USO2AO();
    virtual void AO2USO();

    /// Actually build V_AO
    virtual void compute_V() = 0;
    /// Set things up
    void common_init();
public:
    VBase(boost::shared_ptr<SuperFunctional> functional,
        boost::shared_ptr<BasisSet> primary,
        Options& options);
    virtual ~VBase();
    
    static boost::shared_ptr<VBase> build_V(Options& options, const std::string& type = "RV");

    boost::shared_ptr<BasisSet> basis() const { return primary_; }
    boost::shared_ptr<SuperFunctional> functional() const { return functional_; }
    boost::shared_ptr<PointFunctions> properties() const { return properties_; }
    boost::shared_ptr<DFTGrid> grid() const { return grid_; }
    std::map<std::string, double>& quadrature_values() { return quad_values_; }

    /// Grab this, clear, and push Cocc matrices (with symmetry) to change GS density
    std::vector<SharedMatrix>& C() { return C_; }
    std::vector<SharedMatrix>& Caocc() { return Caocc_; }
    std::vector<SharedMatrix>& Cavir() { return Cavir_; }
    std::vector<SharedMatrix>& P() { return P_; }
    const std::vector<SharedMatrix>& V() const { return V_; }
    const std::vector<SharedMatrix>& D() const { return D_; }

    void set_print(int print) { print_ = print; }
    void set_debug(int debug) { debug_ = debug; }

    virtual void initialize();
    virtual void compute();
    virtual void finalize();

    virtual void print_header() const;
};

// => APPLIED CLASSES <= //

class RV : public VBase {

protected:

    // Actually build V_AO
    virtual void compute_V();

public:
    RV(boost::shared_ptr<SuperFunctional> functional,
        boost::shared_ptr<BasisSet> primary,
        Options& options);
    virtual ~RV();

    virtual void initialize();
    virtual void finalize();

    virtual void print_header() const;
};

class UV : public VBase {

protected:

    // Actually build V_AO
    virtual void compute_V();

public:
    UV(boost::shared_ptr<SuperFunctional> functional,
        boost::shared_ptr<BasisSet> primary,
        Options& options);
    virtual ~UV();

    virtual void initialize();
    virtual void finalize();

    virtual void print_header() const;
};

class RK : public RV {

protected:

    // Actually build V_AO
    virtual void compute_V();

public:
    RK(boost::shared_ptr<SuperFunctional> functional,
        boost::shared_ptr<BasisSet> primary,
        Options& options);
    virtual ~RK();

    virtual void print_header() const;
};

class UK : public UV {

protected:

    // Actually build V_AO
    virtual void compute_V();

public:
    UK(boost::shared_ptr<SuperFunctional> functional,
        boost::shared_ptr<BasisSet> primary,
        Options& options);
    virtual ~UK();

    virtual void print_header() const;
};

}
#endif

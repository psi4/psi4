#ifndef libscf_solver_dft_H
#define libscf_solver_dft_H

#include <psi4-dec.h>
#include <psiconfig.h>

namespace psi {

namespace functional {
class SuperFunctional;
}

class PSIO;
class BasisSet;
class Matrix;
class Vector;
class IntVector;
class Vector3;
class DFTGrid;
class RKSFunctions;
class UKSFunctions;
class Dimension;

namespace scf {

class KSPotential {

protected:
    /// Debug flag
    int debug_;
    /// Print flag
    int print_;
    /// Options object, used to build grid
    Options& options_;
    /// Molecule grid is built around
    boost::shared_ptr<Molecule> molecule_;
    /// Basis set used in the integration
    boost::shared_ptr<BasisSet> primary_;
    /// Desired superfunctional kernal
    boost::shared_ptr<functional::SuperFunctional> functional_;
    /// Integration grid, built by KSPotential
    boost::shared_ptr<DFTGrid> grid_;
    /// Quadrature values obtained during integration 
    std::map<std::string, double> quad_values_;
    /// AO2USO matrix (if not C1)
    SharedMatrix AO2USO_;

    /// Common setup
    void common_init();
    /// Setup the integration grid
    void buildGrid();
    /// Setup AO2USO (if not C1)
    void buildAO2USO();

public:
    KSPotential(boost::shared_ptr<functional::SuperFunctional> functional,
        boost::shared_ptr<Molecule> molecule,
        boost::shared_ptr<BasisSet> primary,
        Options& options);
    virtual ~KSPotential();

    boost::shared_ptr<Molecule> molecule() const { return molecule_; }
    boost::shared_ptr<BasisSet> basis() const { return primary_; }
    boost::shared_ptr<functional::SuperFunctional> functional() const { return functional_; }
    boost::shared_ptr<DFTGrid> grid() const { return grid_; }
    double quadrature_value(const std::string& key);

    virtual void print(FILE* out = outfile, int print = 2);
};

class RKSPotential : public KSPotential {

protected:
    /// Functions object
    boost::shared_ptr<RKSFunctions> properties_;
   
    /// AO D matrix (built or assigned internally)
    SharedMatrix D_AO_;
    /// AO C matrix (built or assigned internally)
    SharedMatrix C_AO_;
 
    /// Target V matrix, AO
    SharedMatrix V_AO_;
    /// Target V matrix, USO
    SharedMatrix V_USO_;

    /// Setup properties
    void buildProperties();
    /// Build D_AO_, C_AO_, allocate V_AO_ if need be
    void USO2AO(SharedMatrix D_USO, SharedMatrix C_USO,
        boost::shared_ptr<Dimension> noccpi);
    /// V_AO_->V_USO_
    void AO2USO();

public:
    RKSPotential(boost::shared_ptr<functional::SuperFunctional> functional,
        boost::shared_ptr<Molecule> molecule,
        boost::shared_ptr<BasisSet> primary,
        Options& options);
    virtual ~RKSPotential();

    void buildPotential(SharedMatrix D_USO, SharedMatrix C_USO,
        boost::shared_ptr<Dimension> noccpi);

    SharedMatrix V_AO() const { return V_AO_; }
    SharedMatrix V_USO() const { return V_USO_; }

    virtual void print(FILE* out = outfile, int print = 2);
};

class UKSPotential : public KSPotential {

protected:
    /// Functions object
    boost::shared_ptr<UKSFunctions> properties_;
   
    /// AO D matrix (built or assigned internally)
    SharedMatrix Da_AO_;
    /// AO C matrix (built or assigned internally)
    SharedMatrix Ca_AO_;
    /// AO D matrix (built or assigned internally)
    SharedMatrix Db_AO_;
    /// AO C matrix (built or assigned internally)
    SharedMatrix Cb_AO_;
 
    /// Target V matrix, AO
    SharedMatrix Va_AO_;
    /// Target V matrix, USO
    SharedMatrix Va_USO_;
    /// Target V matrix, AO
    SharedMatrix Vb_AO_;
    /// Target V matrix, USO
    SharedMatrix Vb_USO_;

    /// Setup properties
    void buildProperties();
    /// Build D_AO_, C_AO_, allocate V_AO_ if need be
    void USO2AO(SharedMatrix Da_USO, SharedMatrix Ca_USO, boost::shared_ptr<Dimension> napi,
                SharedMatrix Db_USO, SharedMatrix Cb_USO, boost::shared_ptr<Dimension> nbpi);
    /// V_AO_->V_USO_
    void AO2USO();

public:
    UKSPotential(boost::shared_ptr<functional::SuperFunctional> functional,
        boost::shared_ptr<Molecule> molecule,
        boost::shared_ptr<BasisSet> primary,
        Options& options);
    virtual ~UKSPotential();

    void buildPotential(SharedMatrix Da_USO, SharedMatrix Ca_USO, boost::shared_ptr<Dimension> napi,
                        SharedMatrix Db_USO, SharedMatrix Cb_USO, boost::shared_ptr<Dimension> nbpi);

    SharedMatrix Va_AO() const { return Va_AO_; }
    SharedMatrix Va_USO() const { return Va_USO_; }
    SharedMatrix Vb_AO() const { return Vb_AO_; }
    SharedMatrix Vb_USO() const { return Vb_USO_; }

    virtual void print(FILE* out = outfile, int print = 2);
};

}} // Namespace psi::scf
#endif

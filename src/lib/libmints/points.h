#ifndef libmints_points_H
#define libmints_points_H

#include <psi4-dec.h>
#include <psiconfig.h>
#include <boost/tuple/tuple.hpp>

namespace psi {

class PSIO;
class BasisSet;
class Matrix;
class Vector;
class IntVector;
class Vector3;
class BlockOPoints;

class PointFunctions {

protected:
    /// Basis set for this PointFunctions
    boost::shared_ptr<BasisSet> primary_;
    /// Current number of points
    int npoints_;
    /// Maximum number of points in a block
    int max_points_;
    /// Maximum number of functions in a block
    int max_functions_;
    /// Maximum derivative to compute
    int deriv_;
    /// Map of value names to Matrices containing values
    std::map<std::string, SharedMatrix > basis_values_;
    /// Map of temp names to Matrices containing temps
    std::map<std::string, SharedMatrix > basis_temps_;
    /// [L]: pure_index, cart_index, coef
    std::vector<std::vector<boost::tuple<int,int,double> > > spherical_transforms_;

    /// Setup spherical_transforms_
    void build_spherical();

public:
    PointFunctions(boost::shared_ptr<BasisSet> primary, int max_points, int max_functions);
    virtual ~PointFunctions();

    void computePoints(boost::shared_ptr<BlockOPoints> block);

    void set_derivative(int deriv);
    SharedMatrix basis_value(const std::string& key);

    virtual void print(FILE* out = outfile, int print = 2) const;

    int max_functions() const { return max_functions_; }
    int max_points() const { return max_points_; }
    int npoints() const { return npoints_; }
 
};

class RKSFunctions : public PointFunctions {

protected:
    /// Density matrix, AO
    SharedMatrix D_AO_;
    /// Occupied C matrix, AO 
    SharedMatrix Cocc_AO_;

    // => Temps <= //
    /// Buffer for half-transform
    SharedMatrix temp_;
    /// Buffer for KE density
    SharedMatrix meta_temp_;
    /// Local D matrix
    SharedMatrix D_local_;
    /// Local Cocc matrix
    SharedMatrix C_local_;

    /// RKS Ansatz (0 - LSDA, 1 - GGA, 2 - Meta-GGA)
    int ansatz_; 
    /// Map of value names to Vectors containing values
    std::map<std::string, boost::shared_ptr<Vector> > property_values_;

    ///Build temporary work arrays
    void build_temps();

public:
    RKSFunctions(boost::shared_ptr<BasisSet> primary, int max_points, int max_functions);
    virtual ~RKSFunctions();

    /// Set RKS Ansatz (0 - LSDA, 1 - GGA, 2 - Meta-GGA)
    void set_ansatz(int ansatz);
    /// Reset pointers, in case the D/Cocc matrices change (nocc *may* change) 
    void  reset_pointers(SharedMatrix D_AO, SharedMatrix Cocc_AO);

    /// Scratch array of size max_points() x max_functions()
    SharedMatrix scratch() const { return temp_; }

    void computeProperties(boost::shared_ptr<BlockOPoints> block);
    boost::shared_ptr<Vector> property_value(const std::string& key);

    virtual void print(FILE* out = outfile, int print = 2) const;
}; 

class UKSFunctions : public PointFunctions {

protected:
    /// Density matrix, AO
    SharedMatrix Da_AO_;
    /// Occupied C matrix, AO 
    SharedMatrix Caocc_AO_;
    /// Density matrix, AO
    SharedMatrix Db_AO_;
    /// Occupied C matrix, AO 
    SharedMatrix Cbocc_AO_;

    // => Temps <= //
    /// Buffer for half-transform
    SharedMatrix tempa_;
    /// Buffer for half-transform
    SharedMatrix tempb_;
    /// Buffer for KE density
    SharedMatrix meta_temp_;
    /// Local D matrix
    SharedMatrix Da_local_;
    /// Local Cocc matrix
    SharedMatrix Ca_local_;
    /// Local D matrix
    SharedMatrix Db_local_;
    /// Local Cocc matrix
    SharedMatrix Cb_local_;

    /// RKS Ansatz (0 - LSDA, 1 - GGA, 2 - Meta-GGA)
    int ansatz_; 
    /// Map of value names to Vectors containing values
    std::map<std::string, boost::shared_ptr<Vector> > property_values_;

    ///Build temporary work arrays
    void build_temps();

public:
    UKSFunctions(boost::shared_ptr<BasisSet> primary, int max_points, int max_functions);
    virtual ~UKSFunctions();

    /// Set RKS Ansatz (0 - LSDA, 1 - GGA, 2 - Meta-GGA)
    void set_ansatz(int ansatz);
    /// Reset pointers, in case the D/Cocc matrices change (nocc *may* change) 
    void  reset_pointers(SharedMatrix Da_AO, SharedMatrix Caocc_AO,
                         SharedMatrix Db_AO, SharedMatrix Cbocc_AO);

    /// Scratch array of size max_points() x max_functions()
    SharedMatrix scratchA() const { return tempa_; }
    /// Scratch array of size max_points() x max_functions()
    SharedMatrix scratchB() const { return tempb_; }

    void computeProperties(boost::shared_ptr<BlockOPoints> block);
    boost::shared_ptr<Vector> property_value(const std::string& key);

    virtual void print(FILE* out = outfile, int print = 2) const;
}; 


}
#endif

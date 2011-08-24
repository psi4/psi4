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
    /// Maximum derivative to compute
    int deriv_;
    /// Map of value names to Matrices containing values
    std::map<std::string, boost::shared_ptr<Matrix> > basis_values_;
    /// Map of temp names to Matrices containing temps
    std::map<std::string, boost::shared_ptr<Matrix> > basis_temps_;
    /// [L]: pure_index, cart_index, coef
    std::vector<std::vector<boost::tuple<int,int,double> > > spherical_transforms_;

    /// Setup spherical_transforms_
    void build_spherical();

public:
    PointFunctions(boost::shared_ptr<BasisSet> primary, int max_points);
    virtual ~PointFunctions();

    void computePoints(boost::shared_ptr<BlockOPoints> block);

    void set_derivative(int deriv);
    boost::shared_ptr<Matrix> basis_value(const std::string& key);

    virtual void print(FILE* out = outfile, int print = 2) const;

    int max_points() const { return max_points_; }
    int npoints() const { return npoints_; }
 
};

class RKSFunctions : public PointFunctions {

protected:
    /// Density matrix, AO
    boost::shared_ptr<Matrix> D_AO_;
    /// Occupied C matrix, AO 
    boost::shared_ptr<Matrix> Cocc_AO_;

    // => Temps <= //
    /// Buffer for half-transform
    boost::shared_ptr<Matrix> temp_;
    /// Local D matrix
    boost::shared_ptr<Matrix> D_local_;
    /// Local Cocc matrix
    boost::shared_ptr<Matrix> C_local_;

    /// RKS Ansatz (0 - LSDA, 1 - GGA, 2 - Meta-GGA)
    int ansatz_; 
    /// Map of value names to Vectors containing values
    std::map<std::string, boost::shared_ptr<Vector> > property_values_;

    ///Build temporary work arrays
    void build_temps();

public:
    RKSFunctions(boost::shared_ptr<BasisSet> primary, int max_points);
    virtual ~RKSFunctions();

    /// Set RKS Ansatz (0 - LSDA, 1 - GGA, 2 - Meta-GGA)
    void set_ansatz(int ansatz);
    /// Reset pointers, in case the D/Cocc matrices change (nocc *may* change) 
    void  reset_pointers(boost::shared_ptr<Matrix> D_AO, boost::shared_ptr<Matrix> Cocc_AO);

    void computeProperties(boost::shared_ptr<BlockOPoints> block);

    boost::shared_ptr<Vector> property_value(const std::string& key);

    virtual void print(FILE* out = outfile, int print = 2) const;
}; 

class UKSFunctions : public PointFunctions {

protected:
    /// Density matrix, AO
    boost::shared_ptr<Matrix> Da_AO_;
    /// Occupied C matrix, AO 
    boost::shared_ptr<Matrix> Caocc_AO_;
    /// Density matrix, AO
    boost::shared_ptr<Matrix> Db_AO_;
    /// Occupied C matrix, AO 
    boost::shared_ptr<Matrix> Cbocc_AO_;

    // => Temps <= //
    /// Buffer for half-transform
    boost::shared_ptr<Matrix> tempa_;
    /// Buffer for half-transform
    boost::shared_ptr<Matrix> tempb_;
    /// Local D matrix
    boost::shared_ptr<Matrix> Da_local_;
    /// Local Cocc matrix
    boost::shared_ptr<Matrix> Ca_local_;
    /// Local D matrix
    boost::shared_ptr<Matrix> Db_local_;
    /// Local Cocc matrix
    boost::shared_ptr<Matrix> Cb_local_;

    /// RKS Ansatz (0 - LSDA, 1 - GGA, 2 - Meta-GGA)
    int ansatz_; 
    /// Map of value names to Vectors containing values
    std::map<std::string, boost::shared_ptr<Vector> > property_values_;

    ///Build temporary work arrays
    void build_temps();

public:
    UKSFunctions(boost::shared_ptr<BasisSet> primary, int max_points);
    virtual ~UKSFunctions();

    /// Set RKS Ansatz (0 - LSDA, 1 - GGA, 2 - Meta-GGA)
    void set_ansatz(int ansatz);
    /// Reset pointers, in case the D/Cocc matrices change (nocc *may* change) 
    void  reset_pointers(boost::shared_ptr<Matrix> Da_AO, boost::shared_ptr<Matrix> Caocc_AO,
                         boost::shared_ptr<Matrix> Db_AO, boost::shared_ptr<Matrix> Cbocc_AO);

    void computeProperties(boost::shared_ptr<BlockOPoints> block);

    boost::shared_ptr<Vector> property_value(const std::string& key);

    virtual void print(FILE* out = outfile, int print = 2) const;
}; 


}
#endif

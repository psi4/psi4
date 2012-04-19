#ifndef libfock_points_H
#define libfock_points_H

#include <cstdio>
#include <map>

#include <libmints/typedefs.h>
#include <boost/tuple/tuple.hpp>

namespace psi {

class PSIO;
class BasisSet;
class Matrix;
class Vector;
class IntVector;
class Vector3;
class BlockOPoints;

extern FILE* outfile;

class BasisFunctions {

protected:
    /// Basis set for this BasisFunctions
    boost::shared_ptr<BasisSet> primary_;
    /// Pure AM or not
    bool puream_;
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
    /// Allocate registers
    virtual void allocate();

public:
    // => Constructors <= //

    BasisFunctions(boost::shared_ptr<BasisSet> primary, int max_points, int max_functions);
    virtual ~BasisFunctions();

    // => Computers <= //

    void compute_functions(boost::shared_ptr<BlockOPoints> block);

    // => Accessors <= //

    SharedMatrix basis_value(const std::string& key);
    std::map<std::string, SharedMatrix>& basis_values() { return basis_values_; }

    int max_functions() const { return max_functions_; }
    int max_points() const { return max_points_; }
    int deriv() const { return deriv_; }

    virtual void print(FILE* out = outfile, int print = 2) const;
    
    // => Setters <= //

    void set_deriv(int deriv) { deriv_ = deriv; allocate(); }
    void set_max_functions(int max_functions) { max_functions_ = max_functions; allocate(); }
    void set_max_points(int max_points) { max_points_ = max_points; allocate(); }
};

class PointFunctions : public BasisFunctions {

protected:
    /// Ansatz (0 - LSDA, 1 - GGA, 2 - Meta-GGA)
    int ansatz_;
    /// Map of value names to Vectors containing values
    std::map<std::string, boost::shared_ptr<Vector> > point_values_;

public:
    // => Constructors <= //

    PointFunctions(boost::shared_ptr<BasisSet> primary, int max_points, int max_functions);
    virtual ~PointFunctions();
    
    // => Computers <= //

    virtual void compute_points(boost::shared_ptr<BlockOPoints> block) = 0;

    // => Accessors <= //

    boost::shared_ptr<Vector> point_value(const std::string& key);
    std::map<std::string, SharedVector>& point_values() { return point_values_; }

    virtual std::vector<SharedMatrix> scratch() = 0;
    virtual std::vector<SharedMatrix> D_scratch() = 0;

    int ansatz() const { return ansatz_; }

    // => Setters <= //

    void set_ansatz(int ansatz) { ansatz_ = ansatz; deriv_ = ansatz; allocate(); } 
    virtual void set_pointers(SharedMatrix Da_occ_AO) = 0; 
    virtual void set_pointers(SharedMatrix Da_occ_AO, SharedMatrix Db_occ_AO) = 0;

}; 

class RKSFunctions : public PointFunctions {

protected:
    // => Pointers <= //
    
    /// Density matrix, AO
    SharedMatrix D_AO_;

    // => Temps <= //

    /// Buffer for half-transform
    SharedMatrix temp_;
    /// Local D matrix
    SharedMatrix D_local_;

    /// Build temporary work arrays
    void build_temps();
    /// Allocate registers
    void allocate();

public:
    RKSFunctions(boost::shared_ptr<BasisSet> primary, int max_points, int max_functions);
    virtual ~RKSFunctions();

    void set_pointers(SharedMatrix Da_occ_AO);
    void set_pointers(SharedMatrix Da_occ_AO, SharedMatrix Db_occ_AO);

    void compute_points(boost::shared_ptr<BlockOPoints> block);

    std::vector<SharedMatrix> scratch();
    std::vector<SharedMatrix> D_scratch();

    void print(FILE* out = outfile, int print = 2) const;
};

class UKSFunctions : public PointFunctions {

protected:
    // => Pointers <= //
    
    /// Density matrix, AO
    SharedMatrix Da_AO_;
    /// Density matrix, AO
    SharedMatrix Db_AO_;

    // => Temps <= //
    
    /// Buffer for half-transform
    SharedMatrix tempa_;
    /// Buffer for half-transform
    SharedMatrix tempb_;
    /// Local D matrix
    SharedMatrix Da_local_;
    /// Local D matrix
    SharedMatrix Db_local_;

    /// Build temporary work arrays
    void build_temps();
    /// Allocate registers
    void allocate();

public:
    UKSFunctions(boost::shared_ptr<BasisSet> primary, int max_points, int max_functions);
    virtual ~UKSFunctions();

    void set_pointers(SharedMatrix Da_occ_AO);
    void set_pointers(SharedMatrix Da_occ_AO, SharedMatrix Db_occ_AO);

    void compute_points(boost::shared_ptr<BlockOPoints> block);

    std::vector<SharedMatrix> scratch();
    std::vector<SharedMatrix> D_scratch();

    void print(FILE* out = outfile, int print = 2) const;
};


}
#endif

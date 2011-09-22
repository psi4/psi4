#ifndef CDSALCLIST_H
#define CDSALCLIST_H

#include <cstdio>
#include <vector>

namespace boost {
template<class T>
class shared_ptr;
}

namespace psi {

extern FILE *outfile;

class Molecule;
class Matrix;
class MatrixFactory;

class CdSalc {
public:
    class Component {
    public:
        double coef;
        int atom;
        int xyz;

        Component(double coef_, int atom_, int xyz_) :
            coef(coef_), atom(atom_), xyz(xyz_) {}
    };

private:
    /// All the components needed for this transformation
    std::vector<Component> components_;

    /// Symmetry of this SALC compatible for ^ (A1 = 0)
    char irrep_;
public:

    CdSalc(char irrep) : irrep_(irrep) {}

    void add(double coef, int atom, int xyz) {
        components_.push_back(Component(coef,atom, xyz));
    }

    size_t ncomponent() const { return components_.size(); }

    // made const function to access coef,atom,xyz through CdSalc[i]
    const Component& component(int com) const { return components_[com]; }

    char irrep() const { return irrep_; }
    void print() const;
};

class CdSalcWRTAtom {
public:
    class Component {
    public:
        double coef;
        int irrep;
        int salc;

        Component(double coef_, int irrep_, int salc_) :
            coef(coef_), irrep(irrep_), salc(salc_) {}
    };

private:

    /// We split the components up for simplicity
    std::vector<Component> x_;
    std::vector<Component> y_;
    std::vector<Component> z_;

public:
    CdSalcWRTAtom() {}

    void add(int xyz, double coef, int irrep, int salc) {
        if (xyz == 0)
            x_.push_back(Component(coef, irrep, salc));
        else if (xyz == 1)
            y_.push_back(Component(coef, irrep, salc));
        else if (xyz == 2)
            z_.push_back(Component(coef, irrep, salc));
    }

    int nx() const { return x_.size(); }
    int ny() const { return y_.size(); }
    int nz() const { return z_.size(); }

    const Component& x(int n) const { return x_[n]; }
    const Component& y(int n) const { return y_[n]; }
    const Component& z(int n) const { return z_[n]; }

    void print() const;
};

class CdSalcList
{
    boost::shared_ptr<Molecule> molecule_;
    boost::shared_ptr<MatrixFactory> factory_;

    char needed_irreps_;
    bool project_out_translations_;
    bool project_out_rotations_;

    int ncd_;
    int cdsalcpi_[8];
    char *atom_irreps_;
    double **cdsalc2cd_;
    int nirrep_;

    /// Vector of all requested SALCs
    std::vector<CdSalc> salcs_;
    std::vector<CdSalcWRTAtom> atom_salcs_;

public:
    CdSalcList() { throw PSIEXCEPTION("CdSalcList(): U R STOOOPID!!!"); }

    /*! Constructor for generating Cartesian displacement symmetry adapted
     *  linear combinations.
     *
     *  For a gradient calculation, you only need totally symmetric SALCs
     *  (needed_irreps=1).
     *
     *  The project_* ability is not coded, yet.
     *
     *  \param mol Molecule to form SALCs for
     *  \param needed_irreps bit representation for desired irreps (default all)
     *  \param project_out_translations Project out translational SALCs
     *  \param project_out_rotations Project out rotational SALCs
     */
    CdSalcList(boost::shared_ptr<Molecule> mol,
               boost::shared_ptr<MatrixFactory> fact,
               int needed_irreps=0xF,
               bool project_out_translations=true,
               bool project_out_rotations=true);
    virtual ~CdSalcList();

    /*! Returns the number of combintations. It may not be 3n-5 or 3n-6.
     *  The value returned depends on needed_irreps and the project_out*
     *  settings.
     */
    size_t ncd() const { return salcs_.size(); }

    std::vector<boost::shared_ptr<Matrix> > create_matrices(const std::string& basename);
    std::string name_of_component(int component);

    char needed_irreps() const { return needed_irreps_; }
    int nirrep(void) const { return nirrep_; }
    bool project_out_translations() const { return project_out_translations_; }
    bool project_out_rotations() const { return project_out_rotations_; }

    const CdSalc& operator[](int i) const { return salcs_[i]; }

    const CdSalcWRTAtom& atom_salc(int i) const { return atom_salcs_[i]; }

    boost::shared_ptr<Matrix> matrix();

    void print() const;
};

} // namespace psi

#endif // CDSALCLIST_H

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

class CdSalc {
public:
    class Element {
    public:
        double coef;
        int atom;
        int xyz;

        Element(double coef_, int atom_, int xyz_) :
            coef(coef_), atom(atom_), xyz(xyz_) {}
    };

private:
    /// All the components needed for this transformation
    std::vector<Element> components_;

    /// Symmetry of this SALC compatible for ^ (A1 = 0)
    char irrep_;
public:

    CdSalc(char irrep) : irrep_(irrep) {}

    void add(double coef, int atom, int xyz) {
        components_.push_back(Element(coef,atom, xyz));
    }

    size_t ncomponent() const { return components_.size(); }

    char irrep() const { return irrep_; }
    void print() const;
};

class CdSalcList
{
    const boost::shared_ptr<Molecule>& molecule_;

    char needed_irreps_;
    bool project_out_translations_;
    bool project_out_rotations_;

    int ncd_;
    double ***salc_symblock_;
    int cdsalcpi_[8];
    char *atom_irreps_;
    double **cdsalc2cd_;
    int nirrep_;

    /// Vector of all requested SALCs
    std::vector<CdSalc> salcs_;

public:
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
    CdSalcList(const boost::shared_ptr<Molecule>& mol,
               char needed_irreps=0xF,
               bool project_out_translations=true,
               bool project_out_rotations=true);
    virtual ~CdSalcList();

    size_t ncd() const { return salcs_.size(); }

    char needed_irreps() const { return needed_irreps_; }
    bool project_out_translations() const { return project_out_translations_; }
    bool project_out_rotations() const { return project_out_rotations_; }

    const CdSalc& operator[](int i) const { return salcs_[i]; }

    void print() const;
};

} // namespace psi

#endif // CDSALCLIST_H

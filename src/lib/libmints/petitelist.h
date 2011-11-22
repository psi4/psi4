#ifndef _psi_src_lib_libmints_petitelist_h_
#define _psi_src_lib_libmints_petitelist_h_

#include "typedefs.h"

#include <map>
#include <cstdio>
#include <stdint.h>
#include "pointgrp.h"

namespace psi {

extern FILE *outfile;

class BasisSet;
class Molecule;
class IntegralFactory;
class Matrix;
class Dimension;

inline int64_t ij_offset64(int64_t i, int64_t j)
{
    return (i>j) ? (((i*(i+1)) >> 1) + j) : (((j*(j+1)) >> 1) + i);
}

inline int64_t i_offset64(int64_t i)
{
    return ((i*(i+1)) >> 1);
}

/////////////////////////////////////////////////////////////////////////////
// These are helper functions for PetiteList and GenericPetiteList4

/*! @{
 *  Computes atom mappings during symmetry operations. Useful in generating
 *  SO information and Cartesian displacement SALCs.
 *
 *  \param mol Molecule to form mapping matrix from.
 *  \returns Integer matrix of dimension natoms X nirreps.
 */
int **compute_atom_map(const boost::shared_ptr<Molecule> &mol);
int **compute_atom_map(const Molecule* mol);
/// @}

/*! @{
 * Frees atom mapping for created for a molecule by compute_atom_map.
 *
 *  \param atom_map Map to free.
 *  \param mol Molecule used to create atom_map.
 */
void delete_atom_map(int **atom_map, const boost::shared_ptr<Molecule> &mol);
void delete_atom_map(int **atom_map, const Molecule* mol);
/// @}

int **compute_shell_map(int **atom_map, const boost::shared_ptr<BasisSet> &);
void delete_shell_map(int **shell_map, const boost::shared_ptr<BasisSet> &);

/////////////////////////////////////////////////////////////////////////////

struct contribution {
    int bfn;
    double coef;

    contribution();
    contribution(int b, double c);
    ~contribution();
};

struct SO {
    int len;
    int length;
    contribution *cont;

    SO();
    SO(int);
    ~SO();

    SO& operator=(const SO&);

    void set_length(int);
    void reset_length(int);

    // is this equal to so to within a sign
    int equiv(const SO& so);
};

struct SO_block {
    int len;
    SO *so;

    SO_block();
    SO_block(int);
    ~SO_block();

    void set_length(int);
    void reset_length(int);

    int add(SO& s, int i);
    void print(const char *title);
};

struct SOCoefficients{
    std::map<int, double> coefficients;
    int irrep;
    //        Contribution(std::map<int, double> coefficients_, int irrep_):
    //            coefficients(coefficients_), irrep(irrep_){}
    SOCoefficients(): irrep(-1){}
    void add_contribution(int bf, double coeff, int symm){
        if(irrep != -1 && irrep != symm)
            throw PSIEXCEPTION("Contribution::symmetry changed!");
        irrep = symm;
        coefficients[bf] += coeff;
    }

    void print() const {
        fprintf(outfile, "\nSOCoefficients, irrep = %d\n", irrep);
        std::map<int, double>::const_iterator iter;
        for(iter = coefficients.begin(); iter != coefficients.end(); ++iter){
            fprintf(outfile, "Basis function:%d Coefficient: %.5f\n", iter->first, iter->second);
        }
    }

    size_t size() const {return(coefficients.size());}

    void scale_coefficients(double factor){
        std::map<int, double>::iterator iter;
        for(iter = coefficients.begin(); iter != coefficients.end(); ++iter){
            iter->second *= factor;
        }
    }

    void delete_zeros(){
        std::map<int, double>::iterator iter;
        for(iter = coefficients.begin(); iter != coefficients.end();){
            std::map<int, double>::iterator erase_iter = iter++;

            if(fabs(erase_iter->second) < 1E-10) coefficients.erase(erase_iter);
        }

    }
};
/////////////////////////////////////////////////////////////////////////////

class PetiteList
{
    int natom_;
    int nshell_;
    int nunique_shell_;
    int ng_;
    int nirrep_;
    int nblocks_;
    bool c1_;

    boost::shared_ptr<BasisSet> basis_;
    const IntegralFactory* integral_;

    bool include_pure_transform_;

    char *p1_;
    int **atom_map_;
    int **shell_map_;
    int **unique_shell_map_;
    char *lamij_;
    int *nbf_in_ir_;
    unsigned short group_;

    unsigned short *stablizer_;
    int max_stablizer_;

    void init();

public:
    PetiteList(const boost::shared_ptr<BasisSet>&, const boost::shared_ptr<IntegralFactory>&, bool include_pure_transform = false);
    PetiteList(const boost::shared_ptr<BasisSet>&, const IntegralFactory*, bool include_pure_transform = false);
    ~PetiteList();

    bool include_pure_transform() const {return include_pure_transform_;}

    /// The AO basis set used to create this petite list
    boost::shared_ptr<BasisSet> basis() { return basis_; }
    /// The integral factory used to create this petite list
    const IntegralFactory* integral() { return integral_; }
    /// Create a clone of this petite list
    boost::shared_ptr<PetiteList> clone();

    /// Number of irreps
    int nirrep() const { return nirrep_; }
    /// The order of the group
    int order() const { return ng_; }
    /** How a full atom maps to another with the specified symmetry operation
     *  \param n Full atom index
     *  \param g Operation number for mapping
     *  \returns the atom n maps into when g is applied.
     */
    int atom_map(int n, int g) const { return (c1_) ? n : atom_map_[n][g]; }
    /** How a full shell index maps to another with the specified symmetry operation
     *  \param n Full shell index
     *  \param g Operation number for mapping
     *  \returns the shell n maps into when g is applied.
     */
    int shell_map(int n, int g) const { return (c1_) ? n : shell_map_[n][g]; }

    /** How a unique shell index maps to an absolute with the specified symmetry operation
     *  \param n Unique shell index
     *  \param g Operation number for mapping
     *  \returns the absolute shell n maps into when g is applied.
     */
    int unique_shell_map(int n, int g) const { return (c1_) ? n : unique_shell_map_[n][g]; }

    /** Number of functions in irrep.
     *  \param h Irrep of interest.
     */
    int nfunction(int h) const;

    /** Number of blocks in symmetry information. Should be same as nirrep().
     */
    int nblocks() const { return nblocks_; }

    void print(FILE *out=outfile);

    /** The symmetry operations that keep the atom unchanged in bit representation.
     *  \param atom The atom of interest.
     */
    unsigned short stablizer(int atom) const { return (c1_) ? group_ : stablizer_[atom]; }

    int max_stablizer() const { return (c1_) ? 1 : max_stablizer_;}

    /** The bit representation of the symmetry operation in the point group.
     */
    unsigned short group() const { return group_; }

    /** Returns the bit representation of the double coset representation.
     */
    unsigned short dcr(unsigned short subgroup1, unsigned short subgroup2) const {
        std::map<unsigned short,bool> uniqueCosets;
        for(int g = 0; g < 8; ++g){
            int coset = 0;
            if(SKIP_THIS_OPERATOR(group_, g)) continue;
            for(int mu = 0; mu < 8; ++mu){
                if(SKIP_THIS_OPERATOR(subgroup1, mu)) continue;
                for(int nu = 0; nu < 8; ++nu){
                    if(SKIP_THIS_OPERATOR(subgroup2, nu)) continue;
                    coset |= NUM_TO_OPERATOR_ID(mu^g^nu);
                    if(!NUM_TO_OPERATOR_ID(mu^g^nu)) coset |= SymmOps::ID;
                }
            }
            uniqueCosets[coset] = 1;
        }
        std::map<unsigned short, bool>::const_iterator iter = uniqueCosets.begin();
        std::map<unsigned short, bool>::const_iterator stop = uniqueCosets.end();
        int rOperators = 0;
        for(; iter != stop; ++iter){
            int coset = iter->first;
            for(int op = 1; op < 9; ++op){
                if(SKIP_THIS_OPERATOR(coset, op)) continue;
                rOperators |= (coset & SymmOps::ID ? SymmOps::E : NUM_TO_OPERATOR_ID(op));
                break;
            }
        }
        return rOperators;
    }

    /** Number of operators in the point group.
     *  @param group Get this from dcr()
     */
    int dcr_degeneracy(unsigned short group) const {
        int degeneracy = 0;
        for(int op = 0; op < 8; ++op){
            if(SKIP_THIS_OPERATOR(group, op)) continue;
            ++degeneracy;
        }
        return degeneracy;
    }
    unsigned short GnG(unsigned short group1, unsigned short group2) const {
         return group1 & group2;
    }

    std::vector<int> bits_to_operator_list(unsigned short list) {
        std::vector<int> positions;
        int position = 1;
        unsigned short g = group_;
        positions.push_back(0);
        for (int n=1; n<9; n++) {
            if (g & 1) {
                if (g & (list & 1))
                    positions.push_back(position);
                position += 1;
            }

            g >>= 1;
            list >>= 1;
        }
        return positions;
    }

    void print_group(unsigned short group) const;

    /// Returns the number of atomic orbitals in a convenient Dimension object.
    Dimension AO_basisdim();
    /// Returns the number of symmetry orbitals per irrep in a convenient Dimension object.
    Dimension SO_basisdim();

    /// Return the basis function rotation matrix R(g)
    /// @param g index of the group operation
    Matrix* r(int g);

    SO_block* compute_aotoso_info();

    /** @return the AO->SO coefficient matrix. The columns correspond to SOs (see SO_basisdim() )
        and rows to AOs (see AO_basisdim() ).

        This matrix can be used to transform operators from
        AO to SO basis and functions from SO to AO basis.
        An operator in the SO basis is obtained by \f$ X^T O_ao
        X\f$, where \f$X\f$ is the return value of this function and \f$ O_ao \f$
        is the operator in the AO basis.
        A function in the AO basis is obtained by \f$ X F_so \f$, where
        \f$ F_so \f$ is the function in the SO basis.
        */
    SharedMatrix aotoso();

    /** @return the SO->AO coefficient matrix (the inverse of AO->SO; for Abelian point groups it
        is a transpose of AO->SO matrix). The columns correspond to AOs (see AO_basisdim() )
        and rows to SOs (see SO_basisdim() ).

        This matrix can be used to transform operators from
        SO to AO basis and functions from AO to SO basis.
        An operator in the AO basis is obtained by \f$ X^T O_so
        X\f$, where \f$X\f$ is the return value of this function and \f$ O_so \f$
        is the operator in the SO basis.
        A function in the SO basis is obtained by \f$ X F_ao \f$, where
        \f$ F_ao \f$ is the function in the AO basis.
        */
    SharedMatrix sotoao();
};

}

#endif // _psi_src_lib_libmints_petitelist_h_

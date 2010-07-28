#ifndef _psi_src_lib_libmints_petitelist_h_
#define _psi_src_lib_libmints_petitelist_h_

#include <boost/shared_ptr.hpp>

#include <libmints/basisset.h>
#include <libmints/integral.h>

#include <stdint.h>

namespace psi {

class BasisSet;
class IntegralFactory;
class Matrix;

inline int64_t
ij_offset64(int64_t i, int64_t j)
{
  return (i>j) ? (((i*(i+1)) >> 1) + j) : (((j*(j+1)) >> 1) + i);
}

inline int64_t
i_offset64(int64_t i)
{
  return ((i*(i+1)) >> 1);
}

/////////////////////////////////////////////////////////////////////////////
// These are helper functions for PetiteList and GenericPetiteList4

int **compute_atom_map(const boost::shared_ptr<BasisSet> &);
void delete_atom_map(int **atom_map, const boost::shared_ptr<BasisSet> &);

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

/////////////////////////////////////////////////////////////////////////////

class PetiteList
{
    int natom_;
    int nshell_;
    int ng_;
    int nirrep_;
    int nblocks_;
    bool c1_;

    boost::shared_ptr<BasisSet> basis_;
    boost::shared_ptr<IntegralFactory> integral_;

    char *p1_;
    int **atom_map_;
    int **shell_map_;
    char *lamij_;
    int *nbf_in_ir_;

    void init();

public:
    PetiteList(const boost::shared_ptr<BasisSet>&, const boost::shared_ptr<IntegralFactory>&);
    ~PetiteList();

    boost::shared_ptr<BasisSet> basis() { return basis_; }
    boost::shared_ptr<IntegralFactory> integral() { return integral_; }
    boost::shared_ptr<PetiteList> clone() { return boost::shared_ptr<PetiteList>(new PetiteList(basis_, integral_)); }

    int nirrep() const { return nirrep_; }
    int order() const { return ng_; }
    int atom_map(int n, int g) const { return (c1_) ? n : atom_map_[n][g]; }
    int shell_map(int n, int g) const { return (c1_) ? n : shell_map_[n][g]; }
    int lambda(int ij) const { return (c1_) ? 1 : lamij_[ij]; }
    int lambda(int i, int j) const { return (c1_) ? 1 : lamij_[ij_offset64(i, j)]; }
    int in_p1(int n) const { return (c1_) ? 1 : p1_[n]; }
    int in_p2(int ij) const { return (c1_) ? 1 : lamij_[ij]; }
    int in_p2(int i, int j) const { return (c1_) ? 1 : lamij_[ij_offset64(i, j)]; }
    int in_p4(int ij, int kl, int i, int j, int k, int l) const;
    int in_p4(int i, int j, int k, int l) const;

    int nfunction(int i) const { return (c1_) ? basis_->nbf() : nbf_in_ir_[i]; }

    int nblocks() const { return nblocks_; }

    void print();

    int* AO_basisdim();
    int* SO_basisdim();

    /// Return the basis function rotation matrix R(g)
    /// @param g index of the group operation
    Matrix* r(int g);

    /// @return information about the transformation from AOs to SOs
    SO_block* aotoso_info();

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
    Matrix* aotoso();

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
    Matrix* sotoao();
};

inline int
PetiteList::in_p4(int ij, int kl, int i, int j, int k, int l) const
{
    if (c1_)
        return 1;

    int64_t ijkl = i_offset64(ij)+kl;
    int nijkl=1;

    for (int g=1; g < ng_; g++) {
        int64_t gij = ij_offset64(shell_map_[i][g],shell_map_[j][g]);
        int64_t gkl = ij_offset64(shell_map_[k][g],shell_map_[l][g]);
        int64_t gijkl = ij_offset64(gij,gkl);

        if (gijkl > ijkl)
            return 0;
        else if (gijkl == ijkl)
            nijkl++;
    }

    return ng_/nijkl;
}

inline int
PetiteList::in_p4(int i, int j, int k, int l) const
{
    if (c1_)
        return 1;

    int64_t ij = ij_offset64(i,j);
    int64_t kl = ij_offset64(k,l);
    int64_t ijkl = ij_offset64(ij,kl);
    int nijkl=1;

    for (int g=1; g < ng_; g++) {
        int64_t gij = ij_offset64(shell_map_[i][g],shell_map_[j][g]);
        int64_t gkl = ij_offset64(shell_map_[k][g],shell_map_[l][g]);
        int64_t gijkl = ij_offset64(gij,gkl);

        if (gijkl > ijkl)
            return 0;
        else if (gijkl == ijkl)
            nijkl++;
    }

    return ng_/nijkl;
}

}

#endif // _psi_src_lib_libmints_petitelist_h_

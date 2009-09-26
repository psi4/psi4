#ifndef _psi_src_lib_libmints_symmetry_h_
#define _psi_src_lib_libmints_symmetry_h_

/*!
    \file libmints/symmetry.h
    \ingroup MINTS
*/

#include <cstdio>
#include <libchkpt/chkpt.hpp>

#include <libutil/ref.h>
#include <libmints/molecule.h>

namespace psi {
    
/// Not used. Use at your own risk.
class Symmetry
{
    int nirreps_;               // number of irreps
    int max_stabilizer_index_;  // maximum stabilizer index
    int nunique_atoms_;         // number of symmetry unique atoms
    int nunique_shells_;        // number of symmetry unique shells
    int nso_;                   // number of SO's
    int *atom_position_;       // symmetry positions/stabilizers of atoms
    int **ict_;                 // transformation properties of nuclei under symmetry operations
    int *unique_atom_2_atom_;   // unique atom number to full atom number mapping
    int *unique_shell_2_shell_; // unique shell number to full shell number mapping array
    int *sopi_;                 // number of SO per irrep
    int *sym_operation_;        // mapping array between "canonical" and symmetry.h-defined ordering
    int *so2symblk_;            // SO number to symmetry block mapping array
    char *symlabel_;            // symmetry label
    char **irr_labels_;         // labels of irreps
    int **trans_vec_;           // a matrix of nshell*nirreps integers that contains symmetry information
    int nshells_;               // number of shells
public:
    Symmetry(psi::Chkpt* chkpt);
    Symmetry(const Symmetry&);
    ~Symmetry();
    
    int nirreps() const { return nirreps_; }
    int max_stabilizer_index() const { return max_stabilizer_index_; }
    int nunique_atoms() const { return nunique_atoms_; }
    int nunique_shells() const { return nunique_shells_; }
    int nso() const { return nso_; }
    int atom_position(int i) const { return atom_position_[i]; }
    int ict(int i, int j) const { return ict_[i][j]; }
    int unique_atom_2_atom(int i) const { return unique_atom_2_atom_[i]; }
    int unique_shell_2_shell(int i) const { return unique_shell_2_shell_[i]; }
    int sopi(int i) const { return sopi_[i]; }
    int sym_operation(int i) const { return sym_operation_[i]; }
    int so2symblk(int i) const { return so2symblk_[i]; }
    char* symlabel() const { return symlabel_; }
    char* irr_label(int i) const { return irr_labels_[i]; }
    int trans_vec(int i, int j) const { return trans_vec_[i][j]; }
    int nshells() const { return nshells_; }
};

}

#endif // _basis_symmetry_h_

#include <psi4-dec.h>
#include <libciomr/libciomr.h>
#include "molecule.h"
#include "pointgrp.h"
#include "petitelist.h" // For compute_atom_map
#include "cdsalclist.h"

using namespace psi;
using namespace boost;

CdSalcList::CdSalcList(const boost::shared_ptr<Molecule>& mol)
    : molecule_(mol)
{
    // Ensure point group has been set.
    if (!molecule_->point_group()) {
        throw PSIEXCEPTION("CdSalcList::CdSalcList: Molecule point group has not been set.");
    }

    // Obtain handy reference to point group
    PointGroup& pg = *molecule_->point_group().get();
    CharacterTable char_table = pg.char_table();
    int nirrep = char_table.nirrep();

    // Obtain atom mapping of atom * symm op to atom
    int **atom_map = compute_atom_map(molecule_);

    print_int_mat(atom_map, molecule_->natom(), nirrep, outfile);

    for (int uatom=0; uatom < molecule_->nunique(); ++uatom) {
        int atom = molecule_->unique(uatom);

        // Project each displacement
        for (int xyz=0; xyz<3; ++xyz) {

        }
    }

    // Free memory.
    delete_atom_map(atom_map, molecule_);
}

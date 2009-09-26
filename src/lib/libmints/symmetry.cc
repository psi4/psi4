#include <stdexcept>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.hpp>
#include <psifiles.h>

#include <libmints/symmetry.h>

using namespace psi;

Symmetry::Symmetry(Chkpt* chkpt) 
{
    symlabel_ = chkpt->rd_sym_label();
    nirreps_ = chkpt->rd_nirreps();
    nso_ = chkpt->rd_nso();
    nunique_atoms_ = chkpt->rd_num_unique_atom();
    nunique_shells_ = chkpt->rd_num_unique_shell();
    atom_position_ = chkpt->rd_atom_position();
    unique_atom_2_atom_ = chkpt->rd_ua2a();
    unique_shell_2_shell_ = chkpt->rd_us2s();
    sopi_ = chkpt->rd_sopi();
    sym_operation_ = chkpt->rd_symoper();
    irr_labels_ = chkpt->rd_irr_labs();
    ict_ = chkpt->rd_ict();
    nshells_ = chkpt->rd_nshell();
    trans_vec_ = chkpt->rd_shell_transm();

    int count = 0;
    so2symblk_ = new int[nso_];
    for (int i=0; i<nirreps_; ++i)
        for (int j=0; j<sopi_[i]; ++j)
            so2symblk_[count++] = i;
}

// Symmetry(const Symmetry&);
Symmetry::~Symmetry() 
{
    Chkpt::free(sopi_);
    Chkpt::free(atom_position_);
    Chkpt::free(unique_atom_2_atom_);
    Chkpt::free(unique_shell_2_shell_);
    Chkpt::free(sym_operation_);
    delete[] so2symblk_;
    Chkpt::free(symlabel_);
    Chkpt::free<int>(ict_);
    Chkpt::free<int>(trans_vec_);
    Chkpt::free(irr_labels_);    
}

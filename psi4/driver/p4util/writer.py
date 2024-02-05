#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2024 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#
"""Module for writing input files to external codes (e.g., Molden, NBO).
Only contains extensions to the :py:class:`psi4.core.Wavefunction` class.

"""

__all__ = []

from typing import Optional

from psi4 import core

from ..constants import constants


def _write_nbo(self, name: str):
    """Write wavefunction information in *wfn* to *name* in NBO format.

    Parameters
    ----------
    name
        Destination file name for NBO file.

    """
    basisset = self.basisset()
    mints = core.MintsHelper(basisset)
    mol = self.molecule()

    # Populate header and coordinates.
    NBO_file = f" $GENNBO NATOMS = {mol.natom()} NBAS = {basisset.nbf()} BODM "
    if self.nalpha() != self.nbeta():
        NBO_file += f" OPEN"
    NBO_file += " $END\n $NBO       $END\n $COORD\n"
    NBO_file += " GENNBO expects one comment line here. So, here's a comment line.\n"
    for atom in range(mol.natom()):
        NBO_file += f"{mol.true_atomic_number(atom):2d}  {int(mol.Z(atom)):2d}  {constants.bohr2angstroms * mol.x(atom):20.12f} {constants.bohr2angstroms * mol.y(atom):20.12f} {constants.bohr2angstroms * mol.z(atom):20.12f}\n"
    NBO_file += " $END\n"

    # Populate basis function information.
    pure_order = [
        [1],  # s
        [103, 101, 102],  # p
        [255, 252, 253, 254, 251],  # d: z2 xz yz x2-y2 xy
        [351, 352, 353, 354, 355, 356, 357],  # f
        [451, 452, 453, 454, 455, 456, 457, 458, 459],  #g
        [551, 552, 553, 554, 555, 556, 557, 558, 559, 560, 561]  #h
    ]
    # For historical reasons, the code loops over shells first and then basis functions within the shell.
    # This turns out to not give the same ordering as looping over basis functions directly.
    NBO_file += " $BASIS\n"
    center_string = ""
    label_string = ""
    count = 0
    for i in range(basisset.nshell()):
        shell = basisset.shell(i)
        am = shell.am
        for j in range(shell.nfunction):
            if not (count % 10):
                center_string += "\n  CENTER =" if not i else "\n          "
                label_string += "\n   LABEL =" if not i else "\n          "
            center_string += f" {shell.ncenter + 1:6d}"
            if basisset.has_puream():
                label = pure_order[am][j]
            else:
                label = 100 * am + j + 1
            label_string += f" {label:6d}"
            count += 1
    NBO_file += center_string + label_string + "\n $END\n"

    # Populate contraction information.Start with exponents.
    NBO_file += f" $CONTRACT\n  NSHELL = {basisset.nshell():6d}\n    NEXP = {basisset.nprimitive():6d}\n"
    function_nums = ""
    prim_nums = ""
    prim_indices = ""
    exponents = ""
    prim_index = 0
    coefficients = []  # [(AM int, coefficient), ...]
    for i in range(basisset.nshell()):
        if not (i % 10):
            function_nums += "\n   NCOMP =" if not i else "\n          "
            prim_nums += "\n   NPRIM =" if not i else "\n          "
            prim_indices += "\n    NPTR =" if not i else "\n          "
        shell = basisset.shell(i)
        nprim = shell.nprimitive
        function_nums += f" {shell.nfunction:6d}"
        prim_nums += f" {nprim:6d}"
        prim_indices += f" {prim_index + 1:6d}"
        for j in range(nprim):
            if not (prim_index % 4):
                exponents += "\n     EXP =" if not prim_index else "\n          "
            exponents += f"{shell.exp(j):15.6E}"
            prim_index += 1
            coefficients.append((shell.am, shell.coef(j)))
    NBO_file += function_nums + prim_nums + prim_indices + exponents

    # Populate contraction coefficients.Because some basis sets(Poples with S and P) use the same
    # coefficients for multiple angular momenta, we must supply coefficients for all primitives, for all
    # angular momenta.This leads to many zero elements.
    am_labels = ["S", "P", "D", "F", "G", "H"]
    for current_nbo_section_am in range(basisset.max_am() + 1):
        for i, (shell_am, coefficient) in enumerate(coefficients):
            if not (i % 4):
                NBO_file += f"\n      C{am_labels[current_nbo_section_am]} =" if not i else "\n          "
            if shell_am != current_nbo_section_am:
                coefficient = 0
            NBO_file += f"{coefficient:15.6E}"
    NBO_file += "\n $END"

    # That finishes most of the basis information. Next is the overlap. It would be great if we could just dump Psi's AO
    # overlap matrix, but we can 't. Per CCA guidelines, Psi' s Cartesian d and higher AM AOs aren't normalized to 1.
    # While NBO can "fix" this itself, it changes other AO quantities to match and gets the Fock matrix wrong.
    # Let's normalize ourselves instead.
    ao_overlap = mints.ao_overlap().np
    nbf = ao_overlap.shape[0]
    ao_normalizer = ao_overlap.diagonal()**(-1 / 2)

    def normalize(matrix, normalizer):
        return ((matrix * normalizer).T * normalizer).T

    normalized_ao_overlap = normalize(ao_overlap, ao_normalizer)

    def write_ao_quantity(*args):
        string = ""
        count = 0
        for quantity in args:
            for i in range(nbf):
                for j in range(nbf):
                    if not (count % 5):
                        string += "\n  "
                    string += f"{quantity[i][j]:15.6E}"
                    count += 1
        return string

    NBO_file += "\n $OVERLAP"
    NBO_file += write_ao_quantity(normalized_ao_overlap)
    NBO_file += "\n $END"

    normalized_alpha_density = normalize(self.Da_subset("AO"), 1 / ao_normalizer)
    normalized_beta_density = normalize(self.Db_subset("AO"), 1 / ao_normalizer)
    normalized_alpha_fock = normalize(self.Fa_subset("AO"), ao_normalizer)

    NBO_file += "\n $DENSITY"
    if self.same_a_b_dens():
        density = normalized_alpha_density + normalized_beta_density
        NBO_file += write_ao_quantity(density)
    else:
        NBO_file += write_ao_quantity(normalized_alpha_density, normalized_beta_density)
    NBO_file += "\n $END"

    NBO_file += "\n $FOCK"
    if not self.same_a_b_dens():
        normalized_beta_fock = normalize(self.Fb_subset("AO"), ao_normalizer)
        NBO_file += write_ao_quantity(normalized_alpha_fock, normalized_beta_fock)
    else:
        NBO_file += write_ao_quantity(normalized_alpha_fock)
    NBO_file += "\n $END"

    # The last step is to write the MO coefficients.
    NBO_file += "\n $LCAOMO"

    def write_C_matrix(C, count):
        # The C coefficients supplied the missing multiplication by the ao_normalizer in the overlap matrix before.
        # For NBO, we need that multiplication gone.
        C = (C.np.T / ao_normalizer).T
        string = ""
        for i in range(self.nmo()):
            for mu in range(nbf):
                count += 1
                if (count % 5 == 1):
                    string += ("\n  ")
                string += f"{C[mu][i]:15.6E}"
        # Pad linear dependencies
        for i in range((nbf - self.nmo()) * nbf):
            count += 1
            if (count % 5 == 1):
                string += ("\n  ")
            string += f"{0:15.6E}"
        return count, string

    count, alpha_LCAOMO = write_C_matrix(self.Ca_subset("AO", "ALL"), 0)
    NBO_file += alpha_LCAOMO
    if not self.same_a_b_orbs():
        NBO_file += write_C_matrix(self.Cb_subset("AO", "ALL"), count)[1]
    NBO_file += "\n $END\n"

    #Now time to write !
    with open(name, 'w') as f:
        f.write(NBO_file)

core.Wavefunction.write_nbo = _write_nbo

def _write_molden(
    self: core.Wavefunction,
    filename: Optional[str] = None,
    do_virtual: Optional[bool] = None,
    use_natural: bool = False,
):
    """Writes wavefunction information in *wfn* to *filename* in
    molden format. Will write natural orbitals from *density* (MO basis) if supplied.
    Warning! most post-SCF wavefunctions do not build the density as this is often
    much more costly than the energy. In addition, the wavefunction density attributes
    (Da and Db) return the SO density and must be transformed to the MO basis
    to use with this function.

    .. versionadded:: 0.5
       *wfn* parameter passed explicitly

    :returns: None

    :type filename:
    :param filename:

        Destination file name for MOLDEN file. If unspecified (None), a file
        name will be generated from the molecule name.

    :type do_virtual:
    :param do_virtual:

        Do write all the MOs to the MOLDEN file (True) or discard the unoccupied
        MOs (False). Not valid for NO's. If unspecified (None), value taken from
        :term:`MOLDEN_WITH_VIRTUAL <MOLDEN_WITH_VIRTUAL (GLOBALS)>`.

    :type use_natural:
    :param use_natural:

        Write natural orbitals determined from density on wavefunction.

    :examples:

    1. Molden file with the Kohn-Sham orbitals of a DFT calculation.

       >>> E, wfn = energy('b3lyp', return_wfn=True)
       >>> wfn.molden('mycalc.molden')

    2. Molden file with the natural orbitals of a CCSD computation. For correlated methods, an energy call will not compute the density.
       "properties" or "gradient" must be called.

       >>> E, wfn = properties('ccsd', return_wfn=True)
       >>> wfn.molden('ccsd_no.molden', use_natural=True)

    3. To supply a custom density matrix, manually set the Da and Db of the wavefunction.
       This is used, for example, to write natural orbitals coming from a root computed
       by a ``CIWavefunction`` computation, e.g., ``detci``, ``fci``, ``casscf``.
       The first two arguments of :py:meth:`~psi4.core.CIWavefunction.get_opdm`
       can be set to ``n, n`` where n => 0 selects the root to
       write out, provided these roots were computed, see :term:`NUM_ROOTS <NUM_ROOTS (DETCI)>`. The
       third argument controls the spin (``"A"``, ``"B"`` or ``"SUM"``) and the final
       boolean option determines whether inactive orbitals are included.

       >>> E, wfn = energy('detci', return_wfn=True)
       >>> wfn.Da() = wfn.get_opdm(0, 0, "A", True)
       >>> wfn.Db() = wfn.get_opdm(0, 0, "B", True)
       >>> molden(wfn, 'no_root1.molden', use_natural=True)

    """

    if filename is None:
        filename = core.get_writer_file_prefix(self.molecule().name()) + ".molden"

    if do_virtual is None:
        do_virtual = bool(core.get_option("SCF", "MOLDEN_WITH_VIRTUAL"))
    
    basisset = self.basisset()
    mol = self.molecule()
    # Header and geometry (Atom, Atom #, Z, x, y, z)
    mol_string = '[Molden Format]\n[Atoms] (AU)\n'
    for atom in range(mol.natom()):
        mol_string += f"{mol.symbol(atom):2s}  {atom+1:2d}  {int(mol.Z(atom)):3d}   {mol.x(atom):20.10f} {mol.y(atom):20.10f} {mol.z(atom):20.10f}\n"

    # Dump basis set
    mol_string += '[GTO]\n'
    for atom in range(mol.natom()):
        mol_string += f"  {atom+1:d} 0\n"
        for rel_shell_idx in range(basisset.nshell_on_center(atom)):
            abs_shell_idx = basisset.shell_on_center(atom, rel_shell_idx)
            shell = basisset.shell(abs_shell_idx)
            mol_string += f" {shell.amchar:s}{shell.nprimitive:5d}  1.00\n"
            for prim in range(shell.nprimitive):
                mol_string += f"{shell.exp(prim):20.10f} {shell.original_coef(prim):20.10f}\n"
        mol_string += '\n'
    
    # 
    if use_natural:
        # Alphas
        nmopi = self.nmopi()
        #MO_Da = core.Matrix("MO Alpha Density Matrix", nmopi, nmopi)
        #MO_Da.transform(self.Da(), self.Ca().transpose())
        MO_Da = self.Da_subset("MO") #MO_Da.transform(self.Da(), self.Ca())
        NO_Ra = core.Matrix("NO Alpha Rotation Matrix", nmopi, nmopi)
        occupation_a = core.Vector(nmopi)
        MO_Da.diagonalize(NO_Ra, occupation_a, core.DiagonalizeOrder.Descending)
        Ca = core.doublet(self.Ca(), NO_Ra, False, False)
        epsilon_a = occupation_a
        # Betas
        #MO_Db = core.Matrix("MO Beta Density Matrix", nmopi, nmopi)
        #MO_Db.transform(self.Db(), self.Cb().transpose())
        MO_Db = self.Db_subset("MO")
        NO_Rb = core.Matrix("NO Beta Rotation Matrix", nmopi, nmopi)
        occupation_b = core.Vector(nmopi)
        MO_Db.diagonalize(NO_Rb, occupation_b, core.DiagonalizeOrder.Descending)
        Cb = core.doublet(self.Cb(), NO_Rb, False, False)
        epsilon_b = occupation_b

    else:
        Ca = self.Ca()
        Cb = self.Cb()
        occupation_a = self.occupation_a()
        occupation_b = self.occupation_b()
        epsilon_a = self.epsilon_a()
        epsilon_b = self.epsilon_b()


    # Convert C matrices to AO MO basis. Ca_subset costs information about which symmetry an orbital originally had, which is why we can't use it.
    aotoso = self.aotoso()
    Ca_ao_mo = core.doublet(aotoso, Ca, False, False).nph
    Cb_ao_mo = core.doublet(aotoso, Cb, False, False).nph
    ao_overlap = self.mintshelper().ao_overlap().np
    # Convert from Psi4 internal normalization to the unit normalization expected by Molden
    ao_normalizer = ao_overlap.diagonal()**(-1 / 2)
    Ca_ao_mo = core.Matrix.from_array([(i.T / ao_normalizer).T for i in Ca_ao_mo])
    Cb_ao_mo = core.Matrix.from_array([(i.T / ao_normalizer).T for i in Cb_ao_mo])
    
    # Reorder AO x MO matrix to fit Molden conventions
    '''
    Reordering expected by Molden
    P: x, y, z
    5D: D 0, D+1, D-1, D+2, D-2
    6D: xx, yy, zz, xy, xz, yz
    7F: F 0, F+1, F-1, F+2, F-2, F+3, F-3
    10F: xxx, yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz
    9G: G 0, G+1, G-1, G+2, G-2, G+3, G-3, G+4, G-4
    15G: xxxx, yyyy, zzzz, xxxy, xxxz, yyyz, zzzx, zzzy, xxyy, xxzz, yyzz, xxyz, yyxz, zzxy
    
    Molden does not handle angular momenta higher than G
    '''
    molden_cartesian_order = [
        [2,0,1,0,0,0,0,0,0,0,0,0,0,0,0], # p
        [0,3,4,1,5,2,0,0,0,0,0,0,0,0,0], # d
        [0,4,5,3,9,6,1,8,7,2,0,0,0,0,0], # f
        [0,3,4,9,12,10,5,13,14,7,1,6,11,8,2] # g
    ]
    nirrep = self.nirrep()
    count = 0 # Keeps track of count for reordering
    temp_a = Ca_ao_mo.clone() # Placeholders for original AO x MO matrices
    temp_b = Cb_ao_mo.clone()

    for i in range(basisset.nshell()):
        am = basisset.shell(i).am
        if (am == 1 and basisset.has_puream()) or (am > 1 and am < 5 and basisset.shell(i).is_cartesian()):
            for j in range(basisset.shell(i).nfunction):
                for h in range(nirrep):        
                    for k in range(Ca_ao_mo.coldim()[h]):
                        Ca_ao_mo.set(h,count + molden_cartesian_order[am-1][j],k,temp_a.get(h,count+j,k))
                        Cb_ao_mo.set(h,count + molden_cartesian_order[am-1][j],k,temp_b.get(h,count+j,k))
        count += basisset.shell(i).nfunction
        
    # Dump MO information
    if basisset.has_puream():
        # For historical reasons, D and F can go on the same line, but setting D without F implicitly sets F. G must be on its own.
        mol_string += '[5D7F]\n[9G]\n\n'
    ct = mol.point_group().char_table()
    mol_string += '[MO]\n'
    mo_dim = self.nmopi() if do_virtual else (self.doccpi() + self.soccpi())
    
    # Alphas. If Alphas and Betas are the same, then only Alphas with double occupation will be written (see line marked "***")
    mos = []
    for h in range(nirrep):
        for n in range(mo_dim[h]):
            mos.append((epsilon_a.get(h, n), (h, n)))

    # Sort mos based on energy
    def mosSort(element):
        return element[0]
    mos.sort(key=mosSort)

    for i in range(len(mos)):
        h, n = mos[i][1]
        mol_string += f" Sym= {ct.gamma(h).symbol():s}\n Ene= {epsilon_a.get(h, n):24.10e}\n Spin= Alpha\n"
        if self.same_a_b_orbs() and self.epsilon_a() == self.epsilon_b() and self.same_a_b_dens():
            mol_string += f" Occup= {occupation_a.get(h, n) + occupation_b.get(h, n):24.10e}\n"
        else:
            mol_string += f" Occup= {occupation_a.get(h, n):24.10e}\n"
        for so in range(self.nso()):
            mol_string += f"{so+1:3d} {Ca_ao_mo.get(h, so, n):24.10e}\n"

    # Betas
    mos = []
    if not self.same_a_b_orbs() or self.epsilon_a() != self.epsilon_b() or not self.same_a_b_dens():
        for h in range(nirrep):
            for n in range(mo_dim[h]):
                mos.append((self.epsilon_b().get(h, n), (h, n)))
        mos.sort(key=mosSort)
        for i in range(len(mos)):
            h, n = mos[i][1]
            mol_string += f" Sym= {ct.gamma(h).symbol():s}\n Ene= {epsilon_b.get(h, n):24.10e}\n Spin= Beta\n " \
                          f"Occup= {occupation_b.get(h, n):24.10e}\n"
            for so in range(self.nso()):
                mol_string += f"{so+1:3d} {Cb_ao_mo.get(h, so, n):24.10e}\n"

    # Write Molden string to file
    with open(filename,'w') as fn:
        fn.write(mol_string)

core.Wavefunction.write_molden = _write_molden

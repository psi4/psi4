.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2017 The Psi4 Developers.
.. #
.. # The copyrights for code used from other parties are included in
.. # the corresponding files.
.. #
.. # This program is free software; you can redistribute it and/or modify
.. # it under the terms of the GNU General Public License as published by
.. # the Free Software Foundation; either version 2 of the License, or
.. # (at your option) any later version.
.. #
.. # This program is distributed in the hope that it will be useful,
.. # but WITHOUT ANY WARRANTY; without even the implied warranty of
.. # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
.. # GNU General Public License for more details.
.. #
.. # You should have received a copy of the GNU General Public License along
.. # with this program; if not, write to the Free Software Foundation, Inc.,
.. # 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
.. #
.. # @END LICENSE
.. #

.. include:: autodoc_abbr_options_c.rst

.. _`sec:handlingBasisSet`:

BasisSet
========

Circa November 2014, basis set handling in |PSIfour| has been revamped.
Fear not that the beloved libmints BasisSet object has changed- rather,
the user specification of basis sets, programmer's API to BasisSet
constructors, and the construction of BasisSet objects has been changed.

Advantages to New Scheme (aka *Why*)
------------------------------------

- Defaults for fitting basis sets set on a per-atom basis (*e.g.*, DF-SCF
  on metal-organic with cc-pVDZ uses cc-pVDZ-JKFIT for the organic and
  Def2-tzvpp (or something) for the metal) so that the user shouldn't
  experience a failed job on account of incomplete fitting basis sets.

- All default info for auxiliary basis sets in one place. Programmer when
  calling for a new auxiliary BasisSet gives the fitting role if defaults
  need to be computed (e.g., JKFIT) and the orbital basis to compute
  defaults off of (e.g., get_option(BASIS)). This eliminates all the
  "corresponding_jkfit" boilerplate in ``proc.py`` and also means defaults
  can be assigned for non-uniform orbital basis sets.

- Assignment of basis sets to atoms proceeds through "all", "by_symbol"
  (e.g., "Co"), or "by_label" (e.g., H1 or Co_mine). There is *no*
  assignment to atoms by number (except a bit internally where it's safe)
  which can be ambiguous when the Molecule has been fragmented as for SAPT.

- Users don't need to "set basis basisname" after every `molecule {...}`
  definition or activation because basis sets are not attached to the
  molecule at time option is set but at time BasisSet is built. Similarly,
  once can define a `basis basisname {...}` block and use it for multiple
  molecules.



BasisSet gives the option name where any user intentions as to proper
value may be found (DF_BASIS_SCF), the name by which the new basis can be
recalled (get_str('DF_BASIS_SCF')), the fitting role if defaults need to
be computed (JKFIT), the

*How* for Programmers
---------------------

To get a BasisSet object into your module, just call `pyconstruct` where
formerly you called `construct`. There are two flavors, one for orbital
basis sets and one for auxiliary basis sets. There's no difference in the
BasisSet objects they return or even the code used to assemble them- the
two flavors are just for sane argument naming and to establish different
signatures for Boost binding.

Orbital Basis
*************

Give the function a Molecule object for which to build basis, a label for
the basis (generally, BASIS), and a hint for finding the basis. This last
argument gets used to find a python function by that name camoflaged
(that's what ``basis {...}`` blocks in the input file get translated into)
or failing that a string to find a gbs file defining the basis. ::

    // simple
    boost::shared_ptr<BasisSet> primary = BasisSet::pyconstruct_orbital(molecule, 
        "BASIS", "CC-PVDZ");

    // self-contained
    boost::shared_ptr<BasisSet> primary = BasisSet::pyconstruct_orbital(Process::environment.molecule(),
        "BASIS", Process::environment.options.get_str("BASIS"));

Auxiliary Basis
***************

Give the function a Molecule object for which to build basis, a label for
the basis, a hint for finding the basis, a fitting role to apply if
defaults need to be generated, and a hint for finding the orbital basis to
build defaults against. ::

    // simple
    boost::shared_ptr<BasisSet> auxiliary = BasisSet::pyconstruct_auxiliary(molecule,
        "DF_BASIS_SCF", "", 
        "JKFIT", "CC-PVDZ");

    // self-contained and force Spherical
    boost::shared_ptr<BasisSet> auxiliary = BasisSet::pyconstruct_auxiliary(Process::environment.molecule(),
        "DF_BASIS_SCF", Process::environment.options.get_str("DF_BASIS_SCF"), 
        "JKFIT", Process:environment.options.get_str("BASIS"), 1);

Adding Basis Option to Code
***************************

- Register new basis keyword with :source:`src/bin/psi4/read_options.cc`
  (of course). The default should be the empty string. ::

    options.add_str("DF_BASIS_ELST", "");

- Register new basis keyword with the input parser
  :source:`share/python/inputparser.py`. In the main function
  `process_input`, add it to the regex below. This ensures that users can
  define ``basis_keyword basis_name {...}`` blocks where the contents of
  the block get associated with basis_name and assigned to your
  basis_keyword. ::

    basis_block = re.compile(r'^(\s*?)(basis|df_basis_scf|df_basis_mp2|df_basis_cc|df_basis_sapt)[=\s]*(\w*?)\s*\{(.*?)\}',
                             re.MULTILINE | re.DOTALL | re.IGNORECASE)

Deprecated Steps: Don't do these anymore!
-----------------------------------------

- Deprecated Step: registering non-basis keywords that contain the word
  BASIS in `check_for_basis` function in :source:`src/bin/psi4/python.cc`.
  Don't do this anymore!

- Deprecated Step: adding `corresponding_rifit` and surrounding
  boilerplate to `run_{method}` function :source:`share/python/proc.py`.
  Don't do this anymore!

- Building a parser object in module code as preparation to building a
  BasisSet. ::

    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser(old_forced_puream));
    boost::shared_ptr<BasisSet> dual_basis = BasisSet::construct(parser, molecule_, "DUAL_BASIS_SCF");





key: label by which basis name gets attached to mol's CoordEntry-s



TODO    check that can just pass stirng instead of options.get_str("BASIS") etc.
        make form fn on the fly. do this in pyconstruct instead of inputparser for ordinary (non-block) basis sets?

//  BasisSet::pyconstruct(mol, key, ?, ?, target)
//  mol, key, target --> smol, targetfunc --> (s O)
//  primary_   = BasisSet::pyconstruct(molecule_, "BASIS", options_.get_str("BASIS"));
//  auxiliary_ = BasisSet::pyconstruct(molecule_, "DF_BASIS_MP2", options.get_str("DF_BASIS_MP2"),
//      "RIFIT", options_get_str("BASIS"));
//  mol, key, target, role, orbital --> smol, orbitalfunc, targetfunc, role --> (s O O s)
//  BasisSet::pyconstruct(mol, key, aux, role, orb)

* Note that the basis set specification in psi4 does not permit the assignment of basis sets to
an atom *number*. This is because multi-fragment methods (e.g., SAPT, efp) can involve the internal
chopping up and reinstantiation of molecules, which coule make the user's instructions ambiguous.
Thus, basis set specifictation is
    molecule mymol {
        # water dimer where
        O  -2  0 0
        H_hb  -1  0 0
        H  -1 -1 0
        --
        O   2  0 0
        H_hb  1  0 0
        H   2  1 0
    }
    * per molecule
        set basis cc-pv(d+d)z
            --or--
        basis mydz {
            assign cc-pv(d+d)z
        }
    
    * per element
        basis mydz {
            assign cc-pv(d+d)z
            
    * per 


<<< Q for Jet/Andy/Rob >>>

* Shouldn't lock_frame_ be reset to False for set_basis_all_atoms/by_symbol/by_label?
    Need to trigger reeval of symm upon geometry_update(). Doing this with set_shell....

* Ok that maybe can't form a basisset name using a key that's not a keyword

* Ok to remove parser from arg list

* since set_basis by number being removed from user domain, switching it to 0-indexing (more natural
    for c-side prog code) and to not indluce dummies (why give a dummy a basis set)

* ok that symm lowering won't show up until basis built?

* order of searching for basissets

* get approval for bas search order: strings, here, PSIPATH, library (I think)

<<< todo >>>

* make sure PSIDATDIR getting searched right, esp for installed copy

* transfer load_basis printing to output file

* check puream handling

* get correct full PT basis aux sets

* empty mol before adding basis sets in basis {} block
* establish that a basis spec must cover the whole molecule

<<    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
<<    basisset_ = BasisSet::construct(parser, molecule_, "BASIS");

>>    basisset_ = BasisSet::pyconstruct(molecule_, "BASIS", options_.get_str("BASIS"));


<<    boost::shared_ptr<BasisSet> auxiliary = BasisSet::construct(parser, primary->molecule(), "DF_BASIS_SCF");

>>    boost::shared_ptr<BasisSet> auxiliary = BasisSet::pyconstruct(primary->molecule(),
            "DF_BASIS_SCF", options.get_str("DF_BASIS_SCF"), "JKFIT", options.get_str("BASIS"));








     boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser(old_forced_puream));
     molecule_->set_basis_all_atoms(basisname, "DUAL_BASIS_SCF");
     boost::shared_ptr<BasisSet> dual_basis = BasisSet::construct(parser, molecule_, "DUAL_BASIS_SCF");
+    // TODO: oh my, forced_puream!
+    // TODO: oh my, a basis for which a fn hasn't been set in the input translation
+    // TODO: oh my, a non-fitting basis to be looked up (in Mol) not under BASIS
+    //boost::shared_ptr<BasisSet> dual_basis = BasisSet::pyconstruct(molecule_, basisname,
+    //            "DUAL_BASIS_SCF");
+    // TODO: I think Rob was planning to rework this projection bit anyways


* check with everyone about order in which directories searched











.. 
.. See `Best Practices <http://sirius.chem.vt.edu/trac/wiki/BestPractices#point1>`_ 
.. 
.. .. comment options["AO_BASIS"].has_changed()
.. .. comment will return false if the default value is being used, and true if the user specified this keyword in the input.
.. 
.. 
.. .. warning:: |globals__puream| is an exception in that its value and
..    ``has_changed()`` value only reflect what the user has explicitly set.
..    This keyword should not be queried to find out the current
..    |globals__puream| state for the active basis; use instead,
..    ``psi4.MintsHelper().basisset().has_puream()``.
.. 
.. - get 
.. 
..   - :py:func:`~psi4.get_global_option()`
..   - :py:func:`~psi4.get_local_option()`
..   - :py:func:`~psi4.get_option()`
.. 
.. 
..   .. note:: Some options (BASIS, BASIS-like, and PUREAM) should always
..      be used globally (no module argument) with the OptionsState objects.
..      Similarly, within the body of the function, they should always be
..      queried and set globally. Same for FREEZE_CORE.
.. 
.. - **Setting-Up Calculations**
.. 
..   The other types of options calls in python driver functions are (a)


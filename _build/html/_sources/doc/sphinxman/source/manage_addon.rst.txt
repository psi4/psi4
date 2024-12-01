.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2024 The Psi4 Developers.
.. #
.. # The copyrights for code used from other parties are included in
.. # the corresponding files.
.. #
.. # This file is part of Psi4.
.. #
.. # Psi4 is free software; you can redistribute it and/or modify
.. # it under the terms of the GNU Lesser General Public License as published by
.. # the Free Software Foundation, version 3.
.. #
.. # Psi4 is distributed in the hope that it will be useful,
.. # but WITHOUT ANY WARRANTY; without even the implied warranty of
.. # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
.. # GNU Lesser General Public License for more details.
.. #
.. # You should have received a copy of the GNU Lesser General Public License along
.. # with Psi4; if not, write to the Free Software Foundation, Inc.,
.. # 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
.. #
.. # @END LICENSE
.. #

.. include:: autodoc_abbr_options_c.rst

.. _`sec:addAddOns`:

Adding Add-Ons
==============

.. _`faq:addonname`:

How to use an Add-On's name in directory structure, build, and distribution
---------------------------------------------------------------------------

* Select a name. May be mixed case with numerals and underscores
  (*e.g.*, CheMPS2, libefp, PCMSolver, v2rdm_casscf). Shouldn't start with a
  numeral. Needn't start with "lib", even if a library.

* GitHub repository name should be :samp:`{AddOn_name}` or
  :samp:`{AddOn_name}.lower()` (hereafter, :samp:`{addon_name}`. For example: CheMPS2, libefp, pcmsolver,
  v2rdm_casscf.

* CMake project name should be :samp:`{AddOn_name}`. For example:
  ``project(libefp)``, ``project(CheMPS2)``, ``project(PCMSolver)``,
  ``project(v2rdm_casscf)``. Namespacing in the directory structure used
  to detect the addon should have this name (*e.g.*,
  ``share/cmake/CheMPS2``).

* Restricted by the CMake project name, add-ons return CMake variables
  and compile definitions of :samp:`FOUND_{AddOn_name}` and
  :samp:`USING_{AddOn_name}`. For example: ``FOUND_libefp``,
  ``USING_CheMPS2``, ``PCMSolver_LIBRARIES``, ``USING_v2rdm_casscf``.

* The CMake target(s) formed use the full add-on name as the namespace,
  :samp:`{AddOn_name}::{lib_name_without_lib}.lower()`. For example:
  ``libefp::efp``, ``CheMPS2::chemps2``, ``PCMSolver::pcm``,
  ``v2rdm_casscf::v2rdm_casscf``.

* Following the CMake project name (though not restricted to it --
  |PSIfour| managment could change the pattern), the user flag to enable
  an add-on is :samp:`ENABLE_{AddOn_name}`. Note that runtime-only
  add-ons don't go through this enabling process.

* Internally, the ExternalProject_Add and dummy libraries as well as any
  tests/ and external/ subdirectories should all be lowercase,
  :samp:`{addon_name}`.

* The `conda package <https://anaconda.org/psi4/repo>`_ and internal to
  |PSIfour| (that is, the ExternalProject_Add, dummy libraries, and any
  tests/ and external/ subdirectories) should all be lowercase,
  :samp:`{addon_name}`.

* Alternatively, you can do everything mentioned here lowercase and just
  have a different capitalization for an advertising name. After all,
  that's what |PSIfour| does.


.. _`faq:addoncmake`:

How to integrate an Add-On into build, testing, and docs
--------------------------------------------------------

* In all cases, put Add-Ons in alphabetic order, ignoring any "lib" in the name.

* :source:`CMakeLists.txt`

  * Add the :samp:`ENABLE_{AddOn_name}` line

  * Add the :samp:`external_{addon_name}` dependency to the ``psi4-core`` external project

  * Add the :samp:`{AddOn_name}_DIR` variable passing to the ``psi4-core`` external project

* :source:`psi4/CMakeLists.txt`

  * Add a block imitating Libint if Add-On required or CheMPS2 if not
    required

  * If there are shared resources to the external that need
    to be found by |PSIfour| in PSIDATADIR, follow the ``efpfrag``
    pattern of libefp to symlink them in.

* :source:`psi4/src/CMakeLists.txt`

  * No changes should be required unless both (1) code in export_*
    or core.cc needs the :samp:`USING_{AddOn_name}` definition or
    AddOn header includes and (2) no binary |PSIfour| module (as
    opposed to library |PSIfour| module with the AddOn target linked
    is itself a direct dependency of target ``core``. Basically,
    try to leave this file alone, but if there are compile errors,
    add the definitions/headers as needed.

* :source:`psi4/src/psi4/`

  * If a module is needed to interface the AddOn to |PSIfour|, try to
    put "interface" in the name. Follow the pattern of CheMPS2 or gdma.
    If non-required, be sure to conditionalize it with ``if(TARGET
    AddOn::addon)`` in CMake files or ``#ifdef USING_AddOn`` in
    source files.

  * If a separate module is not required, follow the patter of dkh
    or simint with respect to libmints. Again, conditionalize as in
    preceding bullet.

* :source:`external/upstream/`

  * Add a CMakeLists.txt that imitates another AddOn of similar
    language and dependencies. Try to keep the format, messaging,
    and variables passed as similar as possible so that differences
    mean something. If BLAS/LAPACK or other common dependencies in
    :source:`external/common` are needed, be sure to add them to the
    ``DEPENDS`` argument.

  * The usual practice to to get everything cohesive between
    the CMake for the AddOn repository and |PSIfour| and then as a
    last step, mint a tag in the former and add it to two places in
    :samp:`external/upstream/{addon_name}/CMakeLists.txt` and one
    place in :source:`psi4/CMakeLists.txt` so that only that version
    and later are acceptable to |PSIfour| for detecting pre-built.

* :source:`tests/`

  * In :source:`tests/CMakeLists.txt`, add a block adding a tests subdirectory if Add-On enabled

  * Create new subdirectory :samp:`tests/{addon_name}` with a
    CMakeLists.txt. In that add a few tests. Imitate the pattern in
    other subdirs of including the addon prefix to the test name in the
    CMakeLists but not in the test dir name. Make sure the tests get the
    addon CTest label and that at least one of them gets the smoke label.

* :source:`doc/sphinxman/`

  * Create a new `.rst` page, copying one of the Add-Ons with similar
    language and dependency requirements. Edit it
    as appropriate. Add this page to the list in
    :source:`doc/sphinxman/source/interfacing.rst`.

  * Add a bullet to :source:`doc/sphinxman/source/build_planning.rst`

  * Add the new page to the long list in
    :source:`doc/sphinxman/CMakeLists.txt`. If there are any files or
    images referred to, add them to the file, too, following precedent.

else
----

* Build conda packages

  * Recipes in https://github.com/psi4/psi4meta/tree/master/conda-recipes

* |PSIfour| and Add-On Projects Working Together

  * Obligations of the External Project owners are to:

    (1) allow us to contribute some CMake files to your build system
        so that compile flags and dependencies (e.g., BLAS/LAPACK) can be
        consistent with the |PSIfour| build and so the installed project can
        be readily detected by |PSIfour| or any interested party (through a
        CMake imported target).

    (2) provide us a tag at a tested commit/version number so their
        development may be ongoing.

    (3) communicate with us when they've made improvements and minted
        a new tag.

  * In return, for Add-Ons the |PSIfour| project will:

    (1) leave control of their code under your purview.

    (2) maintain any interfacing code needed.

    (3) regularly run integration tests between |PSIfour| and your code.

    (4) build a mostly statically linked conda package so that any
        of your users can obtain a pre-built binary distribution through
        ``conda install addon --channel psi4``.

    (5) provide a development sandbox for your code through |PSIfour| plugins.

    (6) provide conda download counts independent of |PSIfour|.


.. _`faq:readoptions`:

How to name keywords in ``psi4/src/read_options.cc``
----------------------------------------------------

A few guidelines for standardizing option names among modules.

* ``TRIPLES`` (not trip), ``TRIPLETS`` (not trip), ``SINGLES`` (not sing),
  ``SINGLETS`` (not sing)

* ``CONVERGENCE`` (not conv, not converge) and ``TOLERANCE`` (not tol)

* Convergence of a method should be governed by an ``E_CONVERGENCE`` for
  energy and either a ``D_CONVERGENCE`` for density or a ``R_CONVERGENCE``
  for residual/amplitudes. All of these should be doubles- let the input
  parser handle the flexible input format.

* Diis should have a boolean ``DIIS`` (not do_diis, not use_diis) to turn
  on/off diis extrapolation, a ``DIIS_MIN_VECS`` and ``DIIS_MAX_VECS`` for
  minimum and maximum number of diis vectors to use, and a ``DIIS_START``
  which is the iteration at which to start saving vectors for diis. Not all
  modules conform to all these at present, but they're as standardized as
  they can be without changing code.

* ``AMPS`` (not amplitude, not amp) for amplitudes

* ``NUM_`` (not n) for number (e.g., ``NUM_AMPS_PRINT``, ``MAX_NUM_VECS``,
  ``NUM_THREADS``)

* Some names that could be split into multiple words are staying as one.
  Use ``MAXITER``, ``CACHELEVEL``, ``PUREAM``, ``DERTYPE``.

* ``INTS`` (not integrals), also ``OEI`` (not oe_integrals) for
  one-electron integrals and ``TEI`` (not te_integrals) for two-electron
  integrals

* ``PERTURB`` (not pert) for perturbation

* Use ``PRINT`` options to indicate printing to output file. Use ``WRITE``
  options to indicate printing to another file. This probably isn't
  entirely valid now but should be observed in future. The complement to
  ``WRITE`` is ``READ``. ``PRINT``, ``READ``, and ``WRITE`` will usually
  be the last words in an option name.

* Use ``FOLLOW_ROOT`` for the state to be followed in geometry optimizations

* ``WFN`` (not wavefunction)

* You're welcome to use ``WFN`` and ``DERTYPE`` as internal options, but
  plan to have these set by the python driver and mark them as ``!expert``
  options. Really avoid using ``JOBTYPE``.

* You're not welcome to add ``CHARGE`` or ``MULTP`` options. Plan to get
  these quantities from the molecule object. Since we frequently use subsets
  of systems (with their own charge and multiplicity), this is safer.

* Conform. Just grep ``'add' psi4/src/read_options.cc`` to get a list of
  all the option names in |PSIfour| and try to match any conventions you
  find.

* If you have a quantity you'd like to call a cutoff, a threshold, a
  tolerance, or a convergence, consider the following guidelines in naming
  it.

  * If its value is typically greater than ~0.001, give it a name with ``CUTOFF``.

  * If its value is typically less than ~0.001 and quantities being tested
    against the option are more valuable with larger values (e.g.,
    integrals, occupations, eigenvectors), give it a name with ``TOLERANCE``.

  * If its value is typically less than ~0.001 and quantities being tested
    against the option are more valuable with smaller values (e.g., energy
    changes, residual errors, gradients), give it a name with
    ``CONVERGENCE``.

* In deciding how to arrange words in an option name, place the context
  first (e.g., ``MP2_AMPS_PRINT``, ``TRIPLES_DIIS``). This means ``PRINT``
  will generally be at the end of an option name.

* Use ``INTS_TOLERANCE`` (not schwarz_cutoff)

* ``H`` in an option name is reserved for Hamiltonian (or hydrogen).
  Hessian should be ``HESS``.

* All option names should be all caps and separated by underscores.

* If you have an option that instructs your module to do something not too
  computationally intensive and then quit, append ``_EXIT`` to the option
  name.

* Scaling terms (like for scs) should follow the pattern ``MP2_SS_SCALE``
  and ``SAPT_OS_SCALE``.

* ``FRAG`` for fragment.

* ``AVG`` for average.

* For level-shifting, let's try to have it governed by (double)
  ``LEVEL_SHIFT`` only and not a boolean/double combo since the procedure
  can be turned on (role of boolean) if the value (role of double) has
  changed.

* For Tikhonow regularization, use ``TIKONOW_OMEGA``, not regularizer.

* ``SYM`` for symmetry.

* ``OCC`` for occupied/occupation (e.g., ``DOCC``, ``LOCK_OCC``, ``OCC_TOLERANCE``).

* ``COND`` for condition and ``CONDITIONER`` for conditioner.

* ``LOCAL`` (not localize).

* Use ``AO`` and ``MO`` for atomic and molecular orbitals. When 'O' for
  orbitals is too obsure or would make for too short a keyword, as in
  "bool NO" for "Do use natural orbitals", use ``ORBS`` for orbitals. So
  natural orbitals are ``NAT_ORBS`` and Brueckner orbitals are
  ``BRUECKNER_ORBS``.

* ``LEVEL`` (not ``LVL``, not ``LEV``).

* ``EX`` for excitation.

* ``VAL`` for valence.

* ``GEOM`` (not geo, not geometry).

* ``SYM`` (not symm, not symmetry).

* ``FILE`` (unless truly multiple FILES).

* ``WRITE``/``READ`` for info transfer across jobs. ``SAVE``/``RESTART``
  for same in context of restart.

* Damping should interface through option (double) ``DAMPING_PERCENTAGE``,
  where a value of 0.0 indicates no damping.

* Try to avoid ``COMPUTE`` or ``CALC`` in an option name. If it's a
  boolean like "opdm_compute" for "Do compute the one-particle density
  matrix", just use ``OPDM``.

* Properties should be governed by a ``PROPERTIES`` array for the root of
  interest or by a ``PROPERTIES_ALL`` array for all roots in a multi-root
  calc.  Since no module conforms to this right now, use ``PROPERTY``
  alone and ``PROP`` in multi-part option as ``PROP_ROOT``, ``PROP_ALL``,
  ``PROP_SYM`` to conform.

* Use ``DF`` (not ri) for density-fitting and resolution-of-the-identity
  option names. Only the basis sets are staying as -RI since that's what
  EMSL uses.


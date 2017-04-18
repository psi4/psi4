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

.. _`sec:addAddOns`:

Adding Add-Ons
==============

* Project Name

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

* CMake Within |PSIfour| Source

  * In all cases, put Add-Ons in alphabetic order, ignoring any "lib" in the name.

  * :source:`CMakeLists.txt`

    * Add the :samp:`ENABLE_{AddOn_name}` line

    * Add the :samp:`external_{addon_name}` dependency to the ``psi4-core`` external project

    * Add the :samp:`{AddOn_name}_DIR` variable passing to the ``psi4-core`` external project

  * :source:`psi4/CMakeLists.txt`

    * Add a block imitating libint if Add-On required or CheMPS2 if not
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

  * :source:`psi4/external/upstream/`

    * Add a CMakeLists.txt that imitates another AddOn of similar
      language and dependencies. Try to keep the format, messaging,
      and variables passed as similar as possible so that differences
      mean something. If BLAS/LAPACK or other common dependencies in
      :source:`psi4/common/` are needed, be sure to add them to the
      ``DEPENDS`` argument.

    * The usual practice to to get everything cohesive between
      the CMake for the AddOn repository and |PSIfour| and then as a
      last step, mint a tag in the former and add it to two places in
      :samp:`psi4/external/upstream/{addon_name}/CMakeLists.txt` and one
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

* Build conda packages

  * Recipes in https://github.com/psi4/psi4meta/tree/master/conda-recipes

* |PSIfour| and Add-On Projects Working Together

  * Obligations of the External Project owners are to:

    * (1) allow us to contribute some CMake files to your build system
      so that compile flags and dependencies (e.g., BLAS/LAPACK) can be
      consistent with the |PSIfour| build and so the installed project can
      be readily detected by |PSIfour| or any interested party (through a
      CMake imported target).

    * (2) provide us a tag at a tested commit/version number so their
      development may be ongoing.

    * (3) communicate with us when they've made improvements and minted
      a new tag.

  * In return, for Add-Ons the |PSIfour| project will:

    * (1) leave control of their code under your purview.

    * (2) maintain any interfacing code needed.

    * (3) regularly run integration tests between |PSIfour| and your code.

    * (4) build a mostly statically linked conda package so that any
      of your users can obtain a pre-built binary distribution through
      ``conda install addon --channel psi4``.

    * (5) provide a development sandbox for your code through |PSIfour| plugins.

    * (6) provide conda download counts independent of |PSIfour|.



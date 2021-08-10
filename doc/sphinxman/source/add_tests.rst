.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2021 The Psi4 Developers.
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

.. _`faq:add_tests`:

Adding PSIthon Test Cases
=========================

To create a new test case, first make a folder in :source:`tests`. Use hyphens, not spaces or underscores, in the directory name. This directory will need two files. The first is ``CMakeLists.txt``, which is necessary to add the test case to the suite. This file should have the following lines::

    include(TestingMacros)
    
    add_regression_test(directory_name "psi;semicolon_separated-list-of-applicable-test-labels")
    

The labels specify which groups of tests include the test case. The ``psi`` label should always be added, but the other labels are test-specific. The method tested should always be included, and this is often sufficient. If adding a test for an already existing module, the labels for other tests of the module will suggest other labels to add.

A test requiring over 15 minutes should be labeled ``longtests``. A short test under 30 seconds used for general bug checking should be labeled ``quicktests``. A test that confirms |PSIfour| is operational should be labeled ``smoketests``.

The other necessary file is the input file itself, ``input.dat``. The input file should be just a simple input file to run the test, with small additions. ::

    #! RI-SCF cc-pVTZ energy of water, with Z-matrix input and cc-pVTZ-RI auxiliary basis.
    #! Also a bit more to force a second line.
    
    nucenergy = 8.801466202085710  #TEST
    refenergy = -76.05098402733282  #TEST
    
    molecule h2o {
       symmetry c1
       O
       H 1 1.0
       H 1 1.0 2 104.5
    }
    
    set {
       basis cc-pVTZ
       scf_type df
       df_basis_scf cc-pVTZ-RI
       e_convergence 10
    }
    
    thisenergy = energy("hf")
    
    compare_values(nucenergy, h2o.nuclear_repulsion_energy(), 9, "Nuclear repulsion energy")  #TEST
    compare_values(refenergy, thisenergy, 9, "Reference energy")  #TEST
    compare_values(refenergy, get_variable('scf total energy'), 9, "Reference energy")  #TEST

Of those small modifications, first, note the special comment at the top (starting with the ``#!`` comment marker). This should be descriptive since it is inlined into the manual (unless ``!nosample`` in the comment) as a sample input.

The reference values are assigned to variables for later use. The compare_values function (along with several relatives in :source:`psi4/driver/p4util/testing.py` for comparing strings, matrices, etc.) checks that the computed values match these reference values to suitable precision. This function prints an error message and signals that the test failed to the make system, if the values don't match. Any lines of the input associated with the validation process should be flagged with ``#TEST`` at the end of each line, so that they can be removed when copying from the tests to the samples directory.

Finally, add the directory name to the list of tests in :source:`tests/CMakeLists.txt`.

In preparing the test case, turn energy, density, amplitude, and
geometry convergence criteria to very tight levels, and use these
results for reference energies, reference geometries, reference cube
files, *etc.*. Then, either remove or relax the convergence settings,
if these are not a vital part of the test. In choosing the number of
digits for :py:func:`~psi4.compare_values` and other compare_* functions,
select a number looser than the convergence set in the test or the
default convergence for the calculation type (energy, gradient, *etc.*).


.. _`faq:add_psiapi_tests`:

Adding PsiAPI Test Cases
========================

Sometimes you want to add tests that check several variations of a
template job or that test error handling or that are PsiAPI rather than
PSIthon focused. In these cases, you'll want to add to the second test
suite that lives at :source:`tests/pytests`. Presently, the "normal"
(everything in the ``tests/`` directory that isn't in ``tests/pytests/``)
are run through :program:`ctest`, while the pytests are run through :program:`pytest`. In
future, all will be run through Pytest, but the former will still be
run as PSIthon (``psi4 input.dat``) while the latter will still be
run as PsiAPI (``import psi4``). In other words, in designing a test,
choose its mode based on whether PSIthon or PsiAPI suits it better and
whether it's a simple model for users (probably PSIthon) or for expert
users (probably PsiAPI). Both will continue to work in future.

In developing a Pytest test, you probably want to edit it in place,
rather than running :program:`make` after each change. Easiest is from
<objdir>, run ``pytest ../tests/pytests``. Add any filters (``-k
test_name_fragment``) or parallelism (``-n <N>`` if ``pytest-xdist``
installed) or print test names (``-v``) or print warnings (``-rws``). To
see stdout output from an otherwise passing test, easiest to add ``assert
0`` at its end to trigger failure. An important point is that because
they're PsiAPI, ``import psi4`` is happening, so the <objdir> |PSIfour|
module must be in :envvar:`PYTHONPATH`. Also, any call to QCEngine is
using ``which psi4``, so the <objdir> |PSIfour| executable must be in
:envvar:`PATH`. The easiest way to prepare your local environment is to
execute the printout of ``<objdir>/stage/bin/psi4 --psiapi``.

* Test must be in the :source:`tests/pytests/` directory.
* Test file name must start with ``test_``. This is how pytest knows to collect it.
* Test file may contain many tests. To be recognized as a test, the Python function must start with ``test_``.
* No registration required to bring a test to pytest's attention.

There are individual "marks" that can be added to whole tests or parts
of parameterized tests so that they can be run by category (``pytest -m
<mark>`` vs. ``ctest -L <mark>``) rather than just by name (``pytest -k
<name_fragment>`` vs. ``ctest -R <name_fragment>``). Most important are
"quick" and "long" that opt tests into the quick CI suite or out of
the normal full suite. Mark with a decorator for the full test or the
marks argument in a parameterized test. Search "mark" in the test suite
for examples. Use "quick" freely for tests that cover functionality and
are under 15s. Use "long" sparingly to winnow out the longest examples,
particularly those over a minute.


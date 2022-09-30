.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2022 The Psi4 Developers.
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

==========
Test Suite
==========

A test suite plays a vital role in open-source software use and development.

* For a |PSIfour| user, tests provide models of inputs that should work
  "as-is" and a searchable collection of syntax and capabilities.
  The test suite also allows high-quality development snapshots of the
  codebase to be built automatically for users.

* For a user who has |PSIfour| as part of a complex computational molecular software environment, a test suite alongside installed |PSIfour| can be used to show that the |PSIfour| piece is working.

* For a feature developer, adding tests provides confidence that you
  can leave your code untouched and still advertise that the feature works
  years later. With tests, proposed changes to |PSIfour| that break your
  code fall upon the change proposer to fix, rather than being merged
  silently and lying in wait for a concientious user to detect and report
  and then likely falling upon *you* to fix.

* For a general developer, the test suite allows confidence in refactoring, switching out underlying libraries, maintenance, and upgrading.


CTest and pytest, PSIthon and PsiAPI
====================================

In designing a test, sometimes you want it to be a model input for the user in a single file or you don't want a lot of ``psi4.`` or Python syntax cluttering the input.
In this case, follow :ref:`faq:add_psithon_tests` to prepare as PSIthon (``psi4 input.dat``) for, roughly speaking, running through :program:`ctest`.
The PSIthon/CTest test suite occupies the whole of :source:`tests` *except* :source:`tests/pytests`.

At other times you want the test to check several variations of a template job or you want to test error handling or you want to focus on PsiAPI rather than PSIthon or you want to control the compute conditions with environment variables.
In this case, follow :ref:`faq:add_psiapi_tests` to prepare as PsiAPI (``import psi4``) for, roughly speaking, running through :program:`pytest`.
The PsiAPI/pytest test suite occupies :source:`tests/pytests`.

The above description sounds as if there are two disjoint test suites, and you have to run both ``ctest`` and ``pytest`` to fully test |PSIfour|.
This has indeed been the case until March 2022.
The difficulty has been that (1) two test suites is unexpected so some developers don't know to run both; and (2) there are important tests in the PSIthon suite that can't be run on a |PSIfour| installation since CTest only works in a build directory.
Now, by adding an extra file to the test directory (:ref:`faq:psithon_through_pytest`), PSIthon tests can also be run through :program:`pytest`.
This hasn't rolled out to all ~500 PSIthon tests (help wanted), but eventually |PSIfour| can be tested with a single command from a build or from an installation.
Therefore, in designing a test, choose its mode based on whether PSIthon or PsiAPI suits it better and whether it's a simple model for users (probably PSIthon) or for expert users (probably PsiAPI).
Both will continue to work in future, and neither have limitations.


Test Contents
=============

* Most |PSIfour| tests will be integration tests focusing on non-regression of user input to answers, and we insist on having these.
  But if you find unit tests helpful, by all means add them to the test suite.

* Most tests should store reference results (from literature or another implementation or a carefully run |PSIfour| calculation),
  run quantum chemistry, then apply one or more of the :ref:`faq:comparison_functions` so that the test will fail if the answer is unexpected.
  The functions are the same in CTest and pytest, but in the former they are, for example, ``compare_matrices(refmat, mat, ...)`` while in the latter it's asserted, like ``assert compare_matrices(refmat, mat, ...)``.
  The main advantage of the testing functions is that they provide helpful error printing upon failure. Deep down, they're NumPy functions.

* In preparing the test case reference values, aim for the converged value rather than many digits from your computer under default convergence conditions.
  This will make the test more robust for different OS, different BLAS libraries, and variations in SCF cycles.
  Turn energy, density, amplitude, and geometry convergence criteria to very tight levels, and use these results for reference energies, reference geometries, reference cube files, *etc.*.
  Then, either remove or relax the convergence settings, if these are not a vital part of the test.
  In choosing the number of digits for :py:func:`~psi4.compare_values` and other compare_* functions, select a number looser than the convergence set in the test
  or the default convergence for the calculation type (energy, gradient, *etc.*).

* Keep tests as short as possible without sacrificing coverage and variety. Under 30 seconds is a good aim.


.. _`faq:add_psithon_tests`:

Adding PSIthon Test Cases
=========================

To create a new test case, first make a folder in :source:`tests` or, for an addon, a subfolder under the addon folder.
Use hyphens, not spaces or underscores, in the directory name.
Add the directory name to the list of tests in :source:`tests/CMakeLists.txt` or, for an addon, ``tests/<addon>/CMakeLists.txt``.
The test directory will need at least two files, ``CMakeLists.txt`` and ``input.dat``.

``CMakeLists.txt``
------------------

This file adds the test case to the suite. It should have at least the following two uncommented lines::

    include(TestingMacros)

    # if extra files needed
    # file(COPY grid.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

    add_regression_test(directory_name "psi;semicolon_separated-list-of-applicable-test-labels")

    # if minutes long
    # set_tests_properties(isapt1 PROPERTIES COST 300)

The labels specify which groups of tests include the test case for ``ctest -L label`` purposes. The ``psi`` label should always be added, but the other labels are test-specific. The method tested should always be included, and this is often sufficient. If adding a test for an already existing module, the labels for other tests of the module will suggest other labels to add.
Labels have been added as developers needed, so they are not systematic or thorough. If you see labels to add or rename, please do.

A test requiring over 15 minutes should be labeled ``longtests``. A short test under 30 seconds used for general bug checking should be labeled ``quicktests``. A test that confirms |PSIfour| is operational should be labeled ``smoketests``.

If a test needs extra input files like ``grid.dat`` or extra reference files for checking against, like ``fchk``, specify these in the ``CMakeLists.txt`` as shown above. Such tests must be run through ``ctest`` and don't usually work when run "by hand" from the objdir via ``stage/bin/psi4 ../tests/directory_name/input.dat``.

If a test is multiple minutes long, load-balancing a parallel CTest run requires the test to be started early. Use the ``COST`` line as shown above to set a weighting to about the number of seconds the test takes.

``input.dat``
-------------

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
    compare_values(refenergy, variable('scf total energy'), 9, "Reference energy")  #TEST

Of those small modifications, first, note the special comment at the top (starting with the ``#!`` comment marker). This should be descriptive since it is inlined into the manual (unless ``!nosample`` in the comment) as a sample input.

Reference values are often assigned to variables for later use.
The compare_values function (along with several relatives in :source:`psi4/driver/p4util/testing.py` for comparing strings, matrices, etc.) checks that the computed values match these reference values to suitable precision. This function prints an error message and signals that the test failed to the make system, if the values don't match. Any lines of the input associated with the validation process should be flagged with ``#TEST`` at the end of each line, so that they can be removed when copying from the tests to the samples directory.

``output.ref``
--------------

When your test case is in final form, run it locally, rename the output to ``output.ref``, and check it into the repository alongside ``input.dat``.
While this isn't used for any testing machinery, it can be handy for users or developers to consult.

.. _`faq:psithon_through_pytest`:

``test_input.py``
-----------------

Starting March 2022, one can also run tests designed as above for CTest through pytest.
To bring the test to pytest's notice, add a file to the directory named ``test_input.py``.
Below is an example for the :source:`tests/ci-property/test_input.py` ::

    from addons import *

    @ctest_labeler("quick;ci;cas;properties;cart;noc1")
    def test_ci_property():
        ctest_runner(__file__, ["grid.dat"])

This file contains much the same information as the ``CMakeLists.txt``.
The ``def test_ci_property`` contains the name of the test, now with underscores rather than hyphens.
``def test_`` identifies it to pytest as a test.
That part of the function name and the name of the file, ``test_input.py`` are required, but no further registration with CMake is necessary.
Most tests need only the simple form of the runner line ``ctest_runner(__file__)``.
This uses QCEngine machinery to execute ``python psi4 input.dat``.
If additional input files are needed from the test directory, their names can be added to the the second argument list as shown above.
Those additional input files *do* need to be registered in :source:`psi4/CMakeLists.txt`.

Finally, the label string passed to CTest is here handed to pytest, with a few changes:

* ``psi`` added automatically, so exclude it when copying from CTest ``CMakeLists.txt``
* ``cli`` added automatically to distinguish CTest origin from deliberate pytest origin, which have ``api`` added
* ``smoke`` used instead of CTest ``smoketests``
* ``quick`` used instead of CTest ``quicktests``
* ``long`` used instead of CTest ``longtests``
* ``addon`` and ``<name-of-addon>`` added automatically when ``@uusing("<name-of-addon>")`` decorates the test or ``marks=using("<name-of-addon>")`` marks the test

CTest "labels" are called "marks" in pytest.
Any new marks should be added to :source:`pytest.ini`.

Running for Debugging
---------------------

* PSIthon tests that don't need extra files to run are easily run from ``<objdir>`` via ``stage/bin/psi4 ../tests/<test-name>/input.dat, with the output appearing in ``../tests/<test-name>/input.out``.
* All PSIthon tests are runable through CTest, and output files appear in ``<objdir>/tests/<test-name>/output.dat`` and stdout results appear in ``<objdir>/Testing/Temporary/LastTest.log*``.


.. _`faq:add_psiapi_tests`:

Adding PsiAPI Test Cases
========================

To create a new test case, either create a new file or add to an existing file under :source:`tests/pytests`.

* Test must be in the :source:`tests/pytests/` directory.
* Test file name must start with ``test_``. This is how pytest knows to collect it.
* A test file may contain many tests, each of which is an ordinary Python function with name starting ``test_``. This is how pytest knows to collect it.
* No registration required to bring a test to pytest's attention.
* No registration required to bring a test to CMake's attention. If a test needs additional files, register them in :source:`psi4/CMakeLists.txt`.

A few notes on test contents:

* Import testing functions from ``utils`` and use Python assert: ``assert compare_values(expected, computed, ...)``.
* Don't worry about cleaning up files or resetting options. A function in :source:`tests/pytests/conftest.py` does this automatically between every test.
* Especially if using data or functions from outside a test, run a variety of tests at different ``pytest -n <N>`` levels to mix up test ordering. If tests fail that pass when run alone, you've got a function of the same name changing state or some similar correctable phenomenon.

A few notes on test labels:

* For every new test file, add ``pytestmark = [pytest.mark.psi, pytest.mark.api]`` at the top.
  This ensures that every test has the ``psi`` mark and every PsiAPI test has the ``api`` mark to contrast with PSIthon tests with ``cli`` mark.

* There are individual "marks" that can be added to whole tests or parts
  of parameterized tests so that they can be run by category (``pytest -m <mark>``
  vs. ``ctest -L <mark>``) rather than just by name (``pytest -k <name_fragment>``
  vs. ``ctest -R <name_fragment>``). Far more complicated logic is allowed than for
  CTest: ``pytest -m "dftd3 and not api and not long"``.

* The most important marks are "quick" and "long" that opt tests into the quick CI suite or out of
  the normal full suite. Mark with a decorator for the full test or the
  marks argument in a parameterized test. Search "mark" in the test suite
  for examples. Use "quick" freely for tests that cover functionality and
  are under 15s. Use "long" sparingly to winnow out the longest examples,
  particularly those over a minute.

Running for Debugging
---------------------

There are many ways to run pytest, :ref:`faq:pytest`, and three different copies of the test file
(i.e., :source:`tests/pytests/test_mp2.py`, ``<objdir>/stage/lib/PYMOD_INSTALL_LIBDIR/psi4/tests/test_mp2.py``,
``CMAKE_INSTALL_PREFIX/lib/PYMOD_INSTALL_LIBDIR/psi4/tests/test_mp2.py``).
But for developing a pytest test, you probably want to use the first so you can edit it in place rather than running ``cmake --build`` after each change.

* Easiest is from <objdir>, run ``pytest ../tests``. Add any filters (``-k
  test_name_fragment``) or parallelism (``-n <N>``  or ``-n auto`` if ``pytest-xdist``
  installed) or print test names (``-v``) or print warnings (``-rws``).
* An important point is that because they're PsiAPI, ``import psi4`` is happening,
  so the <objdir> |PSIfour| module must be in :envvar:`PYTHONPATH`. Also, any call
  to QCEngine is using ``which psi4``, so the <objdir> |PSIfour| executable must be in
  :envvar:`PATH`. The easiest way to prepare your local environment is to
  execute the printout of ``<objdir>/stage/bin/psi4 --psiapi``.
* To see stdout output from an otherwise passing test, easiest to add ``assert 0``
  at its end to trigger failure.
* If stdout printing is insufficient, and you really need to see ``output.dat`` or other files,
  comment out their deletion in :source:`tests/pytests/conftest.py` and run the single test, deleting
  the file each time (since it appends).


.. _`faq:comparison_functions`:

Comparison Functions
====================

Plain Old Data
--------------

.. function:: psi4.compare_values(expected, computed, atol_exponent, label [, *, **kwargs])
   :noindex:

.. autofunction:: psi4.compare_values(expected, computed [, label, *, **kwargs])

.. comment compare_arrays covered by compare_values

.. autofunction:: qcelemental.testing.compare_values
   :noindex:

.. autofunction:: psi4.compare_integers(expected, computed [, label, *, **kwargs])

.. autofunction:: qcelemental.testing.compare
   :noindex:


Objects
-------

.. function:: psi4.compare_matrices(expected, computed, atol_exponent, label [, *, check_name=False, **kwargs])
   :noindex:

.. autofunction:: psi4.compare_matrices(expected, computed [, label, *, check_name=False, **kwargs])

.. autofunction:: qcelemental.testing.compare_recursive
   :noindex:

.. function:: psi4.compare_vectors(expected, computed, atol_exponent, label [, *, check_name=False, **kwargs])
   :noindex:

.. autofunction:: psi4.compare_vectors(expected, computed [, label, *, check_name=False, **kwargs])

.. function:: psi4.compare_wavefunctions(expected, computed, atol_exponent, label [, *, check_name=False, **kwargs])
   :noindex:

.. autofunction:: psi4.compare_wavefunctions(expected, computed [, label, *, check_name=False, **kwargs])

.. autofunction:: psi4.compare_molrecs(expected, computed [, label, *, check_name=False, **kwargs])

.. autofunction:: qcelemental.testing.compare_molrecs
   :noindex:


File Formats
------------

.. autofunction:: psi4.compare_cubes(expected, computed[, label, *, check_name=False, **kwargs])

.. autofunction:: psi4.compare_fchkfiles

.. autofunction:: psi4.compare_fcidumps

.. autofunction:: psi4.compare_moldenfiles

.. autofunction:: qcdb.compare_vibinfos


Extra QCA Functions
-------------------

.. autofunction:: psi4.compare

.. autofunction:: psi4.compare_recursive


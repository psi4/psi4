.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2019 The Psi4 Developers.
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

Adding Test Cases
=================

To create a new test case, first make a folder in :source:`tests`. The directory name may not contain an underscore; to indicate spaces, use a hyphen instead. This directory will need two files. The first is ``CMakeLists.txt``, which is necessary to add the test case to the suite. This file should have the following lines::

    include(TestingMacros)
    
    add_regression_test(directory_name "psi;semicolon_separated-list-of-applicable-test-labels")
    

The labels specify which groups of tests include the test case. The ``psi`` label should always be added, but the other labels are test-specific. The method tested should always be included, and this is often sufficient. If adding a test for an already existing module, the labels for other tests of the module will suggest other labels to add.

A test requiring over 15 minutes should be labeled ``longtests``. A short test used for general bug checking should be labeled ``quicktests``. A test that confirms |PSIfour| is operational should be labeled ``smoketests``.

The other necessary file is the input file itself, ``input.dat``. The input file should be just a simple input file to run the test, with small modifications. ::

    #! RI-SCF cc-pVTZ energy of water, with Z-matrix input and cc-pVTZ-RI auxilliary basis.
    #! Also a bit more to force a second line.
    
    nucenergy = 8.801466202085710  #TEST
    refenergy = -76.05098402733282 #TEST
    
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
    
    thisenergy = energy('rhf')
    
    compare_values(nucenergy, h2o.nuclear_repulsion_energy(), 9, "Nuclear repulsion energy") #TEST
    compare_values(refenergy, thisenergy, 9, "Reference energy") #TEST
    compare_values(refenergy, get_variable('scf total energy'), 9, "Reference energy") #TEST

Of those small modifications, first, note the special comment at the top (starting with the #! comment marker). This should be very descriptive since it is inlined into the manual (unless !nosample is present in this comment) as a sample input.

The reference values are assigned to variables for later use. The compare_values function (along with several relatives in :source:`psi4/driver/p4util/util.py` for comparing strings, matrices, etc.) checks that the computed values match these reference values to suitable precision. This function prints an error message and signals that the test failed to the make system, if the values don't match. Any lines of the input associated with the validation process should be flagged with #TEST at the end of each line, so that they can be removed when copying from the tests to the samples directory.

Finally, add the directory name to the list of tests in :source:`tests/CMakeLists.txt`.

In preparing the test case, turn energy, density, amplitude, and
geometry convergence criteria to very tight levels, and use these
results for reference energies, reference geometries, reference cube
files, *etc.*. Then, either remove or relax the convergence settings,
if these are not a vital part of the test. In choosing the number of
digits for :py:class:`compare_values` and other compare_* functions,
select a number looser than the convergence set in the test or the
default convergence for the calculation type (energy, gradient, *etc.*).

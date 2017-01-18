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

.. index::
   triple: setting; keywords; database()
   see: db(); database()
   single: database()

.. _`sec:db()`:

Database |w---w| :py:func:`~psi4.driver.wrapper_database.database`
==================================================================

.. codeauthor:: Lori A. Burns
.. sectionauthor:: Lori A. Burns

.. autofunction:: psi4.driver.wrapper_database.database(name, db_name [, func, mode, cp, rlxd, symm, zpe, benchmark, tabulate, subset])

.. index:: 
   pair: database(); output

Output
^^^^^^

At the beginning of a database job is printed a listing of the individual system 
calculations which will be performed. The output snippet below is from the example job [1]
above. It shows each reagent required for the subset of database reactions requested.
Note that this is an un-counterpoise-corrected example, and the wrapper is smart enough
to compute only once the monomer whose energy will be subtracted from each of the three dimers. ::

                    RGC1-HeHe-0.85-dimer
                    RGC1-He-mono-unCP
                    RGC1-HeHe-1.0-dimer
                    RGC1-HeHe-1.5-dimer

At the end of the job, the Requested Energy table is printed that gives the total
energies for the requested model chemistry for each reagent and each reaction, as
well as the stoichiometric weights by which the reagent energies are transformed
into the reaction energy. In this case, the dimer is +1 and the monomer is -2,
indicating the interaction energy is computed from dimer less first monomer
less second (identical) monomer. Error statistics are computed with respect to the reference
energies stored in the database. One of these, the mean absolute deviation, is 
returned by the wrapper as an ordinary Python variable. (For databases
without a stored reference energy, *e.g.*, BASIC, large and meaningless numbers are
printed for error.) The other two tables tabulate the PSI variables requested
through keyword ``tabulate``, in this case the total SCF energy and the number
of atoms in each reagent. ::

   ==> Mp2 Total Energy <==

   ------------------------------------------------------------------------------------------------------
               Reaction          Reaction Value                             Reagent 1           Reagent 2
                                                                             Value Wt            Value Wt
   ------------------------------------------------------------------------------------------------------
         RGC1-HeHe-0.85              0.00004802                        -5.69430495  1      -2.84717649 -2
          RGC1-HeHe-1.0             -0.00000064                        -5.69435362  1      -2.84717649 -2
   ------------------------------------------------------------------------------------------------------

   ==> Requested Energy <==

   ------------------------------------------------------------------------------------------------------
               Reaction     Reaction Energy        Reaction Error           Reagent 1           Reagent 2
                               Ref     Calc [kcal/mol]   [kJ/mol]             [Eh] Wt             [Eh] Wt
   ------------------------------------------------------------------------------------------------------
         RGC1-HeHe-0.85     0.0376   0.0301    -0.0075    -0.0312      -5.69430495  1      -2.84717649 -2
          RGC1-HeHe-1.0    -0.0219  -0.0004     0.0215     0.0899      -5.69435362  1      -2.84717649 -2
   ------------------------------------------------------------------------------------------------------
            Minimal Dev                        -0.0075    -0.0312
            Maximal Dev                         0.0215     0.0899
        Mean Signed Dev                         0.0070     0.0293
      Mean Absolute Dev                         0.0145     0.0605
                RMS Dev                         0.0161     0.0673
   ------------------------------------------------------------------------------------------------------

As well as being printed in the output file, database results from the
``tabulate`` option are available in the input file as ordinary Python
dictionaries ``DB_RGT`` and ``DB_RXN``, indexed firstly by reagent or reaction
name and secondly by the requested PSI variable name. See the first
paragraph of :ref:`sec:createDatabase` for the distinction between
reagents and reactions. For example, an input file like the following
requests a couple variables through ``tabulate`` and then makes use of the
resulting data structures, here, only to print. ::
   
   set basis 6-31g*
   db('mp2','s22',subset='small',tabulate=['CURRENT ENERGY','MP2 CORRELATION ENERGY'])

   from pprint import pprint

   print_stdout('\nDB_RGT')
   pprint(DB_RGT)

   print_stdout('\nDB_RXN')
   pprint(DB_RXN)

   print_stdout('\nmp2 interaction energy of water dimer (S22-2)')
   print_stdout(DB_RXN['S22-2']['CURRENT ENERGY'])

The output to the screen is as follows. ::

   DB_RGT
   {'S22-16-dimer': {'CURRENT ENERGY': -155.37373581838636,
                     'DF-MP2 CORRELATION ENERGY': -0.523870772178089},
    'S22-16-monoA-unCP': {'CURRENT ENERGY': -78.29412053242164,
                          'DF-MP2 CORRELATION ENERGY': -0.2629759351596186},
    'S22-16-monoB-unCP': {'CURRENT ENERGY': -77.07606823017188,
                          'DF-MP2 CORRELATION ENERGY': -0.2594122526144091},
    'S22-2-dimer': {'CURRENT ENERGY': -152.40958884746667,
                    'DF-MP2 CORRELATION ENERGY': -0.3797598812113561},
    'S22-2-monoA-unCP': {'CURRENT ENERGY': -76.19905879745446,
                         'DF-MP2 CORRELATION ENERGY': -0.1887118848315123},
    'S22-2-monoB-unCP': {'CURRENT ENERGY': -76.19902978067739,
                         'DF-MP2 CORRELATION ENERGY': -0.18857384937354635},
    'S22-8-dimer': {'CURRENT ENERGY': -80.67416758080654,
                    'DF-MP2 CORRELATION ENERGY': -0.2844102558783027},
    'S22-8-monoA-unCP': {'CURRENT ENERGY': -40.336952636980364,
                         'DF-MP2 CORRELATION ENERGY': -0.14185962536715307},
    'S22-8-monoB-unCP': {'CURRENT ENERGY': -40.336952636980506,
                         'DF-MP2 CORRELATION ENERGY': -0.14185962536715097}}
   
   DB_RXN
   {'S22-16': {'CURRENT ENERGY': -0.0035470557928363178,
               'DF-MP2 CORRELATION ENERGY': -0.0014825844040612934},
    'S22-2': {'CURRENT ENERGY': -0.011500269334817403,
              'DF-MP2 CORRELATION ENERGY': -0.0024741470062974724},
    'S22-8': {'CURRENT ENERGY': -0.0002623068456699684,
              'DF-MP2 CORRELATION ENERGY': -0.0006910051439986686}}
   
   mp2 interaction energy of water dimer (S22-2)
   -0.0115002693348


.. index:: database(); available

.. _`sec:availableDatabases`:

Available Databases
^^^^^^^^^^^^^^^^^^^

Below are documented for particular databases the availability of the generic
database function options **cp**, **rlxd**, **benchmark**, and the string
options for **subset**. The full reagent member list, which can also be used
in conjunction with **subset**, is not included here for consideration of space
and may be found in the database file. The database Python files are very
readable and should be consulted for more particular questions.

----

.. comment This toctree directive only here to suppress warning at build time.
   include line below is doing the work.

.. toctree::
   :hidden:

   autodoc_available_databases

.. include:: autodoc_available_databases.rst


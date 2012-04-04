
.. include:: autodoc_abbr_options_c.rst

.. index::
   triple: setting; keywords; database()
   see: db(); database()
   single: database()

.. _`sec:db()`:

Database
========

.. autofunction:: wrappers.database(name, db_name [, func, mode, cp, rlxd, symm, zpe, benchmark, tabulate, subset])

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
well as the stoichoimetric weights by which the reagent energies are transfromed
into the reaction energy. In this case, the dimer is +1 and the monomer is -2,
indicating the the interaction energy is computed from dimer less first monomer
less second (identical) monomer. Error statistics are computed with respect to the reference
energies stored in the database. One of these, the mean absolute deviation, is 
returned by the wrapper as an ordinary Python variable. (For databases
without a stored reference energy, e.g., BASIC, large and meaningless numbers are
printed for error.) The other two tables tabulate the PSI variables requested
through keyword ``tabulate``, in this case the total SCF energy and the number
of atoms in each reagent. ::

    ==> Scf Total Energy <==
    
    -----------------------------------------------------------------------------------
             Reaction          Reaction Value              Reagent 1       Reagent 2
                                                            Value Wt        Value Wt
    -----------------------------------------------------------------------------------
       RGC1-HeHe-0.85              0.00011520         -5.71020576  1  -2.85516048 -2
        RGC1-HeHe-1.0              0.00000153         -5.71031943  1  -2.85516048 -2
        RGC1-HeHe-1.5             -0.00000000         -5.71032096  1  -2.85516048 -2
    -----------------------------------------------------------------------------------
    
    ==> Natom <==
    
    -----------------------------------------------------------------------------------
             Reaction          Reaction Value              Reagent 1       Reagent 2
                                                            Value Wt        Value Wt
    -----------------------------------------------------------------------------------
       RGC1-HeHe-0.85              0.00000000          2.00000000  1   1.00000000 -2
        RGC1-HeHe-1.0              0.00000000          2.00000000  1   1.00000000 -2
        RGC1-HeHe-1.5              0.00000000          2.00000000  1   1.00000000 -2
    -----------------------------------------------------------------------------------
    
    ==> Requested Energy <==
    
    -----------------------------------------------------------------------------------
             Reaction     Reaction Energy      Error       Reagent 1       Reagent 2
                             Ref     Calc [kcal/mol]          [H] Wt          [H] Wt
    -----------------------------------------------------------------------------------
       RGC1-HeHe-0.85     0.0376   0.0723     0.0347  -5.71020576  1  -2.85516048 -2
        RGC1-HeHe-1.0    -0.0219   0.0010     0.0228  -5.71031943  1  -2.85516048 -2
        RGC1-HeHe-1.5    -0.0029  -0.0000     0.0029  -5.71032096  1  -2.85516048 -2
    -----------------------------------------------------------------------------------
          Minimal Dev                         0.0029
          Maximal Dev                         0.0347
      Mean Signed Dev                         0.0201
    Mean Absolute Dev                         0.0201
              RMS Dev                         0.0240
    -----------------------------------------------------------------------------------

.. index:: database(); available

.. _`sec:availableDatabases`:

Available Databases
^^^^^^^^^^^^^^^^^^^

Below are documented for particular databases the availibility of the generic
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


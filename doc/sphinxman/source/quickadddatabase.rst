
.. index::
   pair: database(); adding new

.. note:: No recompile of the PSI program is necessary for changes made to
    files in ``$PSIDATADIR``, including those described below.

.. _`sec:createDatabase`:

Creating a Database
===================

A necessary consideration in constructing a database is the distinction
between reagents and reactions. A reagent is a single molecular system
(may be a dimer) whose geometry you are possession of and whose electronic
energy may be of interest. A reaction is a combination of one or more
reagent energies whose value you are interested in and a reference value
for which you may or may not be in possession of. A few examples follow.
In a database of interaction energies, the reagents are dimers and their
component monomers (usually derived from the dimer geometry), and the
reactions are the dimer less monomers energies. In a database of barrier
heights, the reagents are reactants, products, and transition-state
structures, and the reactions are the transition-states less
minimum-energy structures. Possibly you may have a collection of
structures to simply be acted upon in parallel, in which case the
structures are both the reagents and the reactions. The role of the
database.py file is to collect arrays and dictionaries that define the
geometries of reagents (GEOS), their combination into reactions (RXNM &
ACTV), available reference values for reactions (BIND), and brief comments
for reagents and reactions (TAGL). The journey from reagent geometries to
functional database.py file is largely automated, in a process described
below.

* Prepare geometry files
    Assemble xyz files for all intended reagent systems in a directory.
    Follow the rules below for best results. The filename for each xyz
    file should be the name of the system. lowercase or MixedCase is
    preferable (according to Sherrill lab convention). Avoid dashes and
    dots in the name as python won't allow them. If you're determined to
    have dashes and dots, they must be replaced by other characters in the
    process_input line, then translated back in the GEOS section; see
    NBC10.py for an example.

    - The first line for each xyz file should be the number of atoms in the system.

    - The second line for each xyz file can be blank (interpreted as no comment), anything (interpreted as a comment), or two integers and anything (interpreted as charge, multiplicity, and remainder as comment).

    - The third and subsequent lines have four fields: the element symbol and the three cartesian coordinates in angstroms. The atom lines should not contain any dummy atoms (what's the use in cartesian form).  For dimer systems, an algorithm is used to apportion the atoms into two fragments; thus the atoms need not be arranged with all fragmentA atoms before all fragmentB atoms. The algorithm will fail for very closely arranged fragments. For dimers, any charge and multiplicity from the second line will be applied to fragmentA (python); charge and multiplicity may need to be redistributed later in the editing step.

* Run script :source:`lib/scripts/ixyz2database.pl`

    Move into the directory where all your xyz files are located. Run the
    script in place, probably as
    ``$PSIDATADIR/scripts/ixyz2database.pl``. It will ask a number of
    questions about your intended database and generate a python file
    named for your database. Uppercase is preferable for database names
    (according to Sherrill lab convention). Note your choice for the route
    variable for the next step.

* Edit file database.py

    According to your responses in to questions in the ixyz2database.pl script,
    several bullets will be printed of edits you necessarily or optionally
    should make. Copy your new database into :source:`lib/databases`.


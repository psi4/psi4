
.. note:: No recompile of the PSI program is necessary for changes made to
    files in ``$PSIDATADIR``, including those described below.

Defining a Method Alias
=======================

Since quantum chemical methods in PSI4 are accessed through Python functions, and
most important quantities are available as PSI variables, it is straightforward
to create aliases to commonly run calculations or to define hybrid methods. The
``$PSIDATADIR/python/aliases.py`` file is intended for editing by the user for
this purpose.

As an example, the MP2.5 method is the average of MP2 and MP3. The latter is
available through the arbitrary order MPn code and returns all lower energies
along with it in PSI variables. The following is basic code that will compute
and return the MP2.5 energy. ::

    def run_mp2_5(name, **kwargs):
    
        energy('mp3', **kwargs)
        e_scf = PsiMod.get_variable('SCF TOTAL ENERGY')
        ce_mp2 = PsiMod.get_variable('MP2 CORRELATION ENERGY')
        ce_mp3 = PsiMod.get_variable('MP3 CORRELATION ENERGY')
    
        ce_mp25 = 0.5 * (ce_mp2 + ce_mp3)
        e_mp25 = e_scf + ce_mp25
    
        print """  MP2.5 total energy:                      %16.8f\n""" % (e_mp25)
        print """  MP2.5 correlation energy:                %16.8f\n""" % (ce_mp25)
    
        return e_mp25

Compare the above to the method that resides in ``aliases.py``.
The rationale for the changes is indicated in the comments below. ::

    def run_mp2_5(name, **kwargs):
        lowername = name.lower()  # handy variable with name keyword in lowercase
        kwargs = kwargs_lower(kwargs)  # removes case sensitivity in keyword names
    
        # Run detci calculation and collect conventional quantities
        energy('mp3', **kwargs)
        e_scf = PsiMod.get_variable('SCF TOTAL ENERGY')
        ce_mp2 = PsiMod.get_variable('MP2 CORRELATION ENERGY')
        ce_mp3 = PsiMod.get_variable('MP3 CORRELATION ENERGY')
        e_mp2 = e_scf + ce_mp2  # reform mp2 and mp3 total energies for printing
        e_mp3 = e_scf + ce_mp3
    
        # Compute quantities particular to MP2.5
        ce_mp25 = 0.5 * (ce_mp2 + ce_mp3)
        e_mp25 = e_scf + ce_mp25
        PsiMod.set_variable('MP2.5 CORRELATION ENERGY', ce_mp25)  # add new method's important results
        PsiMod.set_variable('MP2.5 TOTAL ENERGY', e_mp25)         #     to PSI variable repository
        PsiMod.set_variable('CURRENT CORRELATION ENERGY', ce_mp25)
        PsiMod.set_variable('CURRENT ENERGY', e_mp25)  # geometry optimizer tracks this variable, permits
                                                       #     MP2.5 finite difference optimizations 
        # build string of title banner and print results
        banners = ''
        banners += """PsiMod.print_out('\\n')\n"""
        banners += """banner(' MP2.5 ')\n"""
        banners += """PsiMod.print_out('\\n')\n\n"""
        exec banners
    
        tables  = ''
        tables += """  SCF total energy:                        %16.8f\n""" % (e_scf)
        tables += """  MP2 total energy:                        %16.8f\n""" % (e_mp2)
        tables += """  MP2.5 total energy:                      %16.8f\n""" % (e_mp25)
        tables += """  MP3 total energy:                        %16.8f\n\n""" % (e_mp3)
        tables += """  MP2 correlation energy:                  %16.8f\n""" % (ce_mp2)
        tables += """  MP2.5 correlation energy:                %16.8f\n""" % (ce_mp25)
        tables += """  MP3 correlation energy:                  %16.8f\n""" % (ce_mp3)
        PsiMod.print_out(tables)  # prints nice header and table of all involved quantities to output file
    
        return e_mp25 

One final step is necessary. At the end of the ``aliases.py`` file, add 
the following line. ::

    procedures['energy']['mp2.5'] = run_mp2_5

This permits the newly defined MP2.5 method to be called in the input file
with the following command. ::

    energy('mp2.5')


Creating a New Database
=======================

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

* Run script ixyz2database.pl

    Move into the directory where all your xyz files are located. Run the
    script in place, probably as
    ``$PSIDATADIR/databases/ixyz2database.pl``. It will ask a number of
    questions about your intended database and generate a python file
    named for your database. Uppercase is preferable for database names
    (according to Sherrill lab convention). Note your choice for the route
    variable for the next step.

* Edit file database.py

    - All routes. Optionally, rearrange the order of reactions in HRXN as this will define the order for the database.

    - All routes. Optionally, add sources for geometries and any reference data to commented lines above the dbse variable.

    - All routes. Optionally, the comment lines of TAGL may be edited in plain text.

    - All routes. Optionally, fill in reference values (in kcal/mol) into BIND.

    - All routes. Optionally, define alternate sets of reference values in the array BIND_ALTREF in the database.py file to be called in a psi4 input file as benchmark='ALTREF' . See S22.py as an example.

    - All routes. Optionally, fill in the least computationally expensive 2-3 reactions into HRXN_SM and the most expensive into HRXN_LG 
    - All routes. Optionally, define subsets of reactions in the array SUBSETARRAY = ['reaction', 'reaction'] in the database.py file to be called in a psi4 input file as subset='SUBSETARRAY'. See NBC10.py as an example.

    - All routes. Necessarily (if charge and multiplicity not read in through line2 = cgmp and nor all neutral singlets), assign correct charge and multiplicity to all reagents.

    - Route 3. Necessarily (if any charge and multiplicity specified for the whole reagent is not intended for fragmentA with neutral singlet fragmentB), apportion charge and multiplicity properly between fragmentA and fragmentB. This is not likely necessary if all subsystems are neutral singlets.

    - Route 2. Necessarily, define the reagents that contribute to each reaction by filling in the empty single-quotes in ACTV. Add or delete lines as necessary if for each reaction more or fewer than three reagents contribute. See NHTBH.py as an example.  

    - Route 2. Necessarily, define the mathematical contribution of reagents to reactions by adding a number (most often +1 or -1) for each reagent to the RXNM of each reaction. See NHTBH.py as an example.

    - All routes. Necessarily, copy your new database into ``$PSIDATADIR/databases``.




.. _`sec:methodAlias`:

.. note:: No recompile of the PSI program is necessary for changes made to
    files in ``$PSIDATADIR``, including those described below.

Defining a Method Alias
=======================

Since quantum chemical methods in PSI4 are accessed through Python functions, and
most important quantities are available as PSI variables, it is straightforward
to create aliases to commonly run calculations or to define hybrid methods. The
:source:`lib/python/aliases.py` file is intended for editing by the user for
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

Compare the above to the method that resides in :source:`lib/python/aliases.py`.
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


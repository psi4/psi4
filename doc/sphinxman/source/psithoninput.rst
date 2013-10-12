
.. include:: autodoc_abbr_options_c.rst

.. _`sec:psithonInput`:

==================================
Psithon: Structuring an Input File
==================================

To allow arbitrarily complex computations to be performed, |PSIfour| was built
upon the Python interpreter. However, to make the input syntax simpler, some
pre-processing of the input file is performed before it is interpreted,
resulting in Python syntax that is customized for PSI, termed Psithon.  In
this section we will describe the essential features of the Psithon language.
|PSIfour| is distributed with an extensive test suite, described in section
:ref:`apdx:testSuite`; the input files for these test cases can be found in the
samples subdirectory of the top-level |PSIfour| source directory, and should
serve as useful examples.

.. index:: physical constants
.. _`sec:physicalConstants`:

Physical Constants
==================

For convenience, the Python interpreter will execute the contents of the
|psirc| file in the current user's home area (if present) before performing any
tasks in the input file.  This allows frequently used python variables to be
automatically defined in all input files.  For example, if we repeatedly make
use of the universal gravitational constant, the following line could be placed
in the |psirc| file ::

    UGC = 6.67384E-11 # m^3 / kg^-1 s^-2

which would make the variable ``UGC`` available in all |PSIfour| input files.
For convenience, the physical constants used within the |PSIfour| code (which
are obtained from the 3rd edition of the IUPAC Green
book [Cohen:GreenBook:2008]_) are also automatically loaded as Psithon
variables (before |psirc| is loaded, so that |psirc| values can be overridden by
the user).

.. _`table:physconst`:

The physical constants used within |PSIfour|, which are automatically
made available within all |PSIfour| input files.

.. literalinclude:: @SFNX_INCLUDE@lib/python/p4const/physconst.py
   :lines: 25-

The ``psi_`` prefix is to prevent clashes with user-defined variables in
|PSIfour| input files.

.. index:: memory
.. _`sec:memory`:

Memory Specification
====================

By default, |PSIfour| assumes that 256 Mb of memory are available. While this is
enough for many computations, many of the algorithms will perform better if
more is available. To specify memory, the ``memory`` keyword should be used. The following
lines are all equivalent methods for specifying that 2 Gb of RAM is available
to |PSIfour|::

    # all equivalent

    memory 2 Gb
    
    memory 2000 Mb
    
    memory 2000000 Kb

One convenient way to override the |PSIfour| default memory is to place a
memory command in the |psirc| file (Sec. :ref:`sec:psirc`). For example,
the following makes the default memory 2 Gb. ::

    set_memory(2000000000)

However, unless you're assured of having only one job running on a node at
a time (and all nodes on the filesystem with |psirc| have similar memory
capacities), it is advised to set memory in the input file on a
per-calculation basis.

.. note:: For parallel jobs, the ``memory`` keyword represents the total memory
   available to the job, *not* the memory per thread.

Molecule and Geometry Specification
===================================

.. toctree::
   :maxdepth: 2

   psithonmol

.. comment To add EFP fragments to a molecule, see :ref:`sec:usingEFPFragments`.

.. index::
   triple: setting; keywords; general
.. _`sec:jobControl`:

Job Control Keywords
====================

|PSIfour| comprises a number of modules, written in C++, that each perform
specific tasks and are callable directly from the Python front end. Each module
recognizes specific keywords in the input file which control its function.
These keywords are detailed in Appendix :ref:`apdx:options_c_module`.
The keywords can be made global, or scoped to apply to
certain specific modules. The following examples demonstrate some of the ways
that global keywords can be specified::

    # all equivalent

    set globals basis cc-pVDZ

    set basis cc-pVDZ

    set globals basis = cc-pVDZ

    set basis = cc-pVDZ

    set globals{
      basis cc-pVDZ
    }
    
    set {
      basis cc-pVDZ
    }
    
    set {
      basis = cc-pVDZ
    }

Note the lack of quotes around ``cc-pVDZ``, even though it is a string. The
Psithon preprocessor automatically wraps any string values in ``set`` commands in
strings. The last three examples provide a more convenient way for specifying
multiple keywords::

    set {
      basis = cc-pVDZ
      print = 1
      reference = rhf
    }

For arguments that require an array input, standard Python list syntax should
be used, *viz.*::

    set {
      docc = [3, 0, 1, 1]
    }

List/matrix inputs may span multiple lines, as long as the opening ``[`` is
on the same line as the name of the keyword.

Any of the above keyword specifications can be scoped to individual modules,
by adding the name of the module after the ``set`` keyword. Omitting the module
name, or using the name ``global`` or ``globals`` will result in the keyword being
applied to all modules. For example, in the following input ::

    molecule{
      o
      h 1 roh
      h 1 roh 2 ahoh
    
      roh = 0.957
      ahoh = 104.5
    }
    
    set basis cc-pVDZ
    set ccenergy print 3
    set scf print 1
    energy('ccsd')

the basis set is set to cc-pVDZ throughout, the SCF code will have a print
level of 1 and the ccenergy code, which performs coupled cluster computations,
will use a print level of 3. In this example a full CCSD computation is
performed by running the SCF code first, then the coupled cluster modules;
the ``energy()`` Python helper function ensures that this is performed correctly.
Note that the Python interpreter executes commands in the order they appear in
the input file, so if the last four commands in the above example were to read ::

    set basis cc-pVDZ
    energy('ccsd')
    set ccenergy print 3
    set scf print 1

the commands that set the print level would be ineffective, as they would be
processed after the CCSD computation completes.

Basis Sets
==========

.. toctree::
   :maxdepth: 2

   quickaddbasis

.. _`sec:psiVariables`:

PSI Variables & Return Values
=============================

To harness the power of Python, |PSIfour| makes the most pertinent results
of each computation available to the Python interpreter for
post-processing. To demonstrate, we can embellish the previous example of
H\ :sub:`2` and H atom::

    molecule h2 {
      H
      H 1 0.9
    }
    
    set basis cc-pvdz
    set reference rhf
    h2_energy = energy('scf')
    
    molecule h {
      H
    }
    
    set basis cc-pvdz
    set reference uhf
    h_energy = energy('scf')
    
    D_e = psi_hartree2kcalmol * (2*h_energy - h2_energy)
    print "De=%f" % D_e

The :py:func:`~driver.energy` function returns the final result of the
computation, the requested total energy in Hartrees, which we assign to a
Python variable. The two energies are then converted to a dissociation
energy and printed to the output file using standard Python notation.

Generally, there are multiple quantities of interest. Appendix
:ref:`apdx:psivariables_module` lists PSI variables variables set by each
module, and :ref:`apdx:psivariables_alpha` defines them.  These can be
accessed through the ``get_variable()`` function.  For example, after
performing a density fitted MP2 computation, both the spin component
scaled energy and the unscaled MP2 energy are made available::

    e_mp2 = get_variable('MP2 TOTAL ENERGY')
    e_scs_mp2 = get_variable('SCS-MP2 TOTAL ENERGY')

Each module and the Python driver set PSI variables over the course of a
calculation.  The values for all can be printed in the output file with
the input file command ``print_variables()``. Note that PSI variables
accumulate over a |PSIfour| instance unless cleared by ``clean_variables()``.
So if you run in a single input file a STO-3G FCI followed by a
aug-cc-pVQZ SCF followed by a ``print_variables()`` command, the last will
include both :psivar:`SCF TOTAL ENERGY <SCFTOTALENERGY>` and :psivar:`FCI
TOTAL ENERGY <FCITOTALENERGY>`. Don't get excited that you got a
high-quality calculation cheaply.

Most of the usual user computation functions (*i.e.*,
:py:func:`~driver.energy`, :py:func:`~driver.optimize`, and
:py:func:`~driver.frequency`) return simply the current total energy.
Consult the descriptions of other functions in :ref:`sec:psithonFunc` for
what quantities they return and for what data structures they make
available for post-processing.

.. _`sec:loops`:

Loops
=====

Python provides many control structures, any of which can be used within |PSIfour|
input files. For example, to loop over three basis sets, the following code can
be used::

    basis_sets = ["cc-pVDZ", "cc-pVTZ", "cc-pVQZ"]
    for basis_set in basis_sets:
        set basis = $basis_set
        energy('scf')

The declaration of ``basis_sets`` is completely standard Python, as is the next
line, which iterates over the list. However, because the Psithon preprocessor
wraps strings in quotes by default, we have to tell it that ``basis_set`` is a
Python variable, not a string, by prefixing it with a dollar sign. 

The geometry specification supports delayed initialization of variable,
which permits potential energy scans. As an example, we can scan both the
angle and bond length in water::

    molecule h2o{
      O
      H1 R
      H1 R2 A
    }
    
    Rvals=[0.9,1.0,1.1]
    Avals=range(102,106,2)
    
    set basis cc-pvdz
    set scf e_convergence=11
    for R in Rvals:
        h2o.R = R
        for A in Avals:
            h2o.A = A
            energy('scf')

The declarations of ``Rvals`` and ``Avals`` are both completely standard Python syntax.
Having named our molecule ``h2o`` we can then set the values of ``R`` and ``A`` within
the loops. Note that we do not need the dollar sign to access the Python
variable in this example; that is required only when using Python variables
with the ``set`` keyword.

Cartesian geometries, because of details of the geometry update process,
need to be specified within the loop(s) along with their basis set when
geometry scans are performed. See :srcsample:`scf4` for analogous Z-matrix
and Cartiesian scans.

.. _`sec:resultsTables`:

Tables of Results
=================

The results of computations can be compactly tabulated with the :py:func:`~text.Table` Psithon
function. For example, in the following potential energy surface scan for water ::

    molecule h2o {
      O
      H 1 R
      H 1 R 2 A
    }
    
    Rvals=[0.9,1.0,1.1]
    Avals=range(100,102,2)
    
    table=Table(rows=["R","A"], cols=["E(SCF)","E(SCS)","E(DFMP2)"])
    
    set basis cc-pvdz
    
    for R in Rvals:
        h2o.R = R
        for A in Avals:
            h2o.A = A
            energy('df-mp2')
            escf = get_variable('SCF TOTAL ENERGY')
            edfmp2 = get_variable('DF-MP2 TOTAL ENERGY')
            escsmp2 = get_variable('SCS-DF-MP2 TOTAL ENERGY')
            table[R][A] = [escf, escsmp2, edfmp2]
    
    print table
    relative=table.copy()
    relative.absolute_to_relative()
    print relative

we first define a table (on line 10) with two row indices and three column
indices. As the potential energy scan is performed, the results are stored
(line 22) and the final table is printed to the output file (line 24). The
table is converted from absolute energies to relative energies (in |kcalpermol|)
on line 26, before being printed again. The relative energies are reported with
respect to the lowest value in each column. More examples of how to control the
formatting of the tables can be found in the sample input files provided; see
Appendix :ref:`apdx:testSuite` for a complete listing.

.. _`sec:wrappers`:

Python Wrappers
===============

The Python foundations of the |PSIfour| driver and Psithon syntax permit
many commonly performed post-processing procedures to be integrated into
the |PSIfour| suite.  Among these are automated computations of
interaction energies through :py:func:`~wrappers.cp`, of a model chemistry
applied to a database of systems through :py:func:`~wrappers.database`,
and of several model chemistries together approximating greater accuracy
through :py:func:`~wrappers.complete_basis_set`.  These are discussed
separately in :ref:`sec:psithonFunc`.  Note that the options documented
for Python functions are placed as arguments in the command that calls the
function, not in the ``set {...}`` block or with any other ``set``
command.


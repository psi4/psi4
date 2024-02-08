.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2023 The Psi4 Developers.
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
samples subdirectory of the top-level |PSIfour| source directory and should
serve as useful examples.

The equivalent Python PsiAPI syntax is shown alongside the Psithon code snippets.
When using the Python API, one must import the |PSIfour| module with::

  import psi4

No such directive is neccesary when using Psithon, which is run using the ``psi4``
executable.

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

    UGC = 6.67384E-11  # m^3 / kg^-1 s^-2

which would make the variable ``UGC`` available in all |PSIfour| input files.
For convenience, the physical constants used within the |PSIfour| code (which
are obtained from `NIST CODATA 2014
<https://physics.nist.gov/cuu/Constants/archive2014.html>`_
are also automatically loaded as Psithon
variables (before |psirc| is loaded, so that the user's |psirc| values can
override the builtins (in the input file, not in the C++ code).

The physical constants used within |PSIfour|, which are automatically
made available within all |PSIfour| input files are in :ref:`table:physconst`.

.. .. literalinclude:: @SFNX_INCLUDE@psi4/driver/constants/physconst.py
..    :lines: 28-

In Psithon input files, prepend physical constants with ``psi_`` to
prevent clashes with user-defined variables (*e.g.*, ``psi_h``). In
PsiAPI mode, access as, *e.g.*, ``psi4.constants.h``.

.. index:: memory
.. _`sec:memory`:

Memory Specification
====================

By default, |PSIfour| assumes that 500 MiB of memory are available. While this is
enough for many computations, many of the algorithms will perform better if
more is available. To specify memory, the ``memory`` keyword should be used. The following
lines are all equivalent methods for specifying that 2 GB of RAM is available
to |PSIfour|::

    # all equivalent

    memory 2 GB
    
    memory 2000 MB
    
    memory 2000000 kB

Please note that memory can be specified both in IEC binary units (1 KiB = 1024 bytes) and SI units (1 kB = 1000 bytes). |PSIfour| recognizes and obeys both of them correctly. The units are not case sensitive (Kb and KB are equivalent to kB).

By default, |PSIfour| performs a "sanity check" when parsing Psithon input files, enforcing a minimum memory requirement of 250 MiB. While it is generally not recomennded to do so, expert users can bypass this check by directly setting the number of bytes availble to |PSIfour|::

    # setting available memory to 2 MB
    set_memory_bytes(2000000)

Please note that this memory setting only governs the maximal memory
usage of the major data structures, and actual total memory usage
is slightly higher. This is usually a negligible amount, except when
setting tiny memory allowances.

One convenient way to override the |PSIfour| default memory is to place a
memory command in the |psirc| file (Sec. :ref:`sec:psirc`). For example,
the following makes the default memory 2 GB. ::

    set_memory(2000000000)

However, unless you're assured of having only one job running on a node at
a time (and all nodes on the filesystem with |psirc| have similar memory
capacities), it is advised to set memory in the input file on a
per-calculation basis.

That same command can be used for PsiAPI mode::

    psi4.set_memory(int(5e8))

.. tabs::

   .. code-tab:: py PSIthon

        set_memory(2000000000)

   .. code-tab:: py PsiAPI

        psi4.set_memory(int(5e8))


.. note:: For parallel jobs, the ``memory`` keyword represents the total memory
   available to the job, *not* the memory per thread.

Molecule and Geometry Specification
===================================

.. toctree::
   :maxdepth: 2

   psithonmol

.. comment To add EFP fragments to a molecule, see :ref:`sec:usingEFPFragments`.

To add EFP fragments to a molecule, see :ref:`sec:usingEFPFragments`.

.. index::
   triple: setting; keywords; general
.. _`sec:jobControl`:

Job Control Keywords
====================

|PSIfour| comprises a number of C++ modules that each perform
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
the :py:func:`~psi4.driver.energy` Python helper function ensures that this is performed correctly.
Note that the Python interpreter executes commands in the order they appear in
the input file, so if the last four commands in the above example were to read ::

    set basis cc-pVDZ
    energy('ccsd')
    set ccenergy print 3
    set scf print 1

the commands that set the print level would be ineffective, as they would be
processed after the CCSD computation completes.

In PsiAPI mode, one can use the command :py:func:`~psi4.driver.set_options`
like below for general and module-specific options. Note that these values
should be of correct type, strings for strings, floats for floats like
convergences. The function :py:func:`~psi4.core.clean_options` that reinitializes
all options may also be useful to separate calculations in a PsiAPI
session. 

.. tabs::

   .. code-tab:: py PSIthon

        set {
            scf_type = pk
            e_convergence =  1.e-5
            soscf = True
            geom_maxiter = 50
        }

   .. code-tab:: py PsiAPI

        psi4.set_options({
            'scf_type': 'pk',
            'e_convergence': 1.e-5,
            'soscf': True,
            'geom_maxiter': 50
        })

Basis Sets
==========

.. toctree::
   :maxdepth: 2

   basissets

.. _`sec:psiVariables`:

PSI Variables
=============

To harness the power of Python, |PSIfour| makes the most pertinent results
of each computation available to the Python interpreter for
post-processing. To demonstrate, we can embellish the previous example of
H\ :sub:`2` and H atom

.. tabs::

   .. code-tab:: py PSIthon

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
        print ("De = %f" % D_e)

   .. code-tab:: py PsiAPI

        h2 = psi4.geometry("""
          H
          H 1 0.9
          """)
        
        psi4.set_options({"basis": "cc-pvdz",
                          "reference": "rhf"})
        h2_energy = psi4.energy('scf')
        
        h = psi4.geometry("""
          H
          """)
        
        psi4.set_options({"basis": "cc-pvdz",
                          "reference": "uhf"})
        h_energy = psi4.energy('scf')
        
        D_e = psi4.constants.hartree2kcalmol * (2*h_energy - h2_energy)
        print("De = %f" % D_e)

The :py:func:`~psi4.driver.energy` function returns the final result of the
computation, the requested total energy in Hartrees, which we assign to a
Python variable. The two energies are then converted to a dissociation
energy and printed to the output file using standard Python notation.

Generally, there are multiple quantities of interest. Appendix
:ref:`apdx:psivariables_module` lists PSI variables variables set by each
module, and :ref:`apdx:psivariables_alpha` defines them.  These can be
accessed through the :py:func:`~psi4.core.variable` function. For example, after
performing a density fitted MP2 computation, both the spin component
scaled energy and the unscaled MP2 energy are made available

.. tabs::

   .. code-tab:: py PSIthon

    e_mp2 = variable('MP2 TOTAL ENERGY')
    e_scs_mp2 = variable('SCS-MP2 TOTAL ENERGY')

   .. code-tab:: py PsiAPI

    e_mp2 = psi4.variable('MP2 TOTAL ENERGY')
    e_scs_mp2 = psi4.variable('SCS-MP2 TOTAL ENERGY')

Each module and the Python driver set PSI variables over the course of a
calculation.  The values for all can be printed in the output file with
the input file command :py:func:`~psi4.core.print_variables`. Note that PSI variables
are cleared at the start of each :py:func:`~psi4.driver.energy`, etc. in an input
file by :py:func:`~psi4.core.clean_variables()`.
So if you run in a single input file a STO-3G FCI followed by a
aug-cc-pVQZ SCF followed by a :py:func:`~psi4.core.print_variables` command, the
last will include :psivar:`SCF TOTAL ENERGY` but not
:psivar:`FCI TOTAL ENERGY`.
The entire dictionary of PSI variables can be obtained through
:py:func:`~psi4.core.get_variables`.

.. _`sec:returnvals`:

Return Values
=============

Most of the usual user computation functions (*i.e.*,
:py:func:`~psi4.driver.energy`, :py:func:`~psi4.driver.optimize`, and
:py:func:`~psi4.driver.frequency`) return simply the current total energy.
Consult the descriptions of other functions in :ref:`sec:psithonFunc` for
what quantities they return and for what data structures they make
available for post-processing. Many users need only deal with the simple return
form for the computation functions. 

.. tabs::

   .. code-tab:: py PSIthon

    # E is total energy float
    # G is gradient array
    # H is hessian array
    # wfn is class instance with many computational details

    # simple returns
    E = energy(...)
    E = optimize(...)
    E = frequency(...)
    G = gradient(...)  # used by optimize()
    H = hessian(...)  # used by frequency()

   .. code-tab:: py PsiAPI

    # E is total energy float
    # G is gradient array
    # H is hessian array
    # wfn is class instance with many computational details

    # simple returns
    E = psi4.energy(...)
    E = psi4.optimize(...)
    E = psi4.frequency(...)
    G = psi4.gradient(...)  # used by optimize()
    H = psi4.hessian(...)  # used by frequency()

For more elaborate post-processing of computations, adding
``return_wfn=True`` keyword argument additionally returns
:py:class:`~psi4.core.Wavefunction`. 

.. tabs::

   .. code-tab:: py PSIthon

    # power user returns
    E, wfn = energy(..., return_wfn=True)
    E, wfn = optimize(..., return_wfn=True)
    E, wfn = frequency(..., return_wfn=True)
    G, wfn = gradient(..., return_wfn=True)  # used by optimize()
    H, wfn = hessian(..., return_wfn=True)  # used by frequency()

    # print gradient array and its rms
    wfn.gradient.print_out()
    print wfn.gradient().rms()

    # format output for other programs
    molden(wfn, 'mycalc.molden')

    # access array in another format
    np.array(wfn.hessian())

   .. code-tab:: py PsiAPI

    # power user returns
    E, wfn = psi4.energy(..., return_wfn=True)
    E, wfn = psi4.optimize(..., return_wfn=True)
    E, wfn = psi4.frequency(..., return_wfn=True)
    G, wfn = psi4.gradient(..., return_wfn=True)  # used by optimize()
    H, wfn = psi4.hessian(..., return_wfn=True)  # used by frequency()

    # print gradient array and its rms
    wfn.gradient.print_out()
    print(wfn.gradient().rms())

    # format output for other programs
    psi4.molden(wfn, 'mycalc.molden')

    # access array in another format
    np.array(wfn.hessian())

.. _`sec:loops`:

Loops
=====

Python provides many control structures, any of which can be used within |PSIfour|
input files. For example, to loop over three basis sets, the following code can
be used:

.. tabs::

   .. code-tab:: py PSIthon
        :force:

        basis_sets = ["cc-pVDZ", "cc-pVTZ", "cc-pVQZ"]
        for basis_set in basis_sets:
            set basis = $basis_set
            energy('scf')

   .. code-tab:: py PsiAPI

        basis_sets = ["cc-pVDZ", "cc-pVTZ", "cc-pVQZ"]
        for basis_set in basis_sets:
            psi4.set_options({"basis": basis_set})
            psi4.energy('scf')

The declaration of ``basis_sets`` is completely standard Python, as is the next
line, which iterates over the list. However, because the Psithon preprocessor
wraps strings in quotes by default, we have to tell it that ``basis_set`` is a
Python variable, not a string, by prefixing it with a dollar sign. 

The geometry specification supports delayed initialization of variable,
which permits potential energy scans. As an example, we can scan both the
angle and bond length in water

.. tabs::

   .. code-tab:: py PSIthon

        molecule h2o{
          O
          H 1 R
          H 1 R 2 A
        }
        
        Rvals = [0.9, 1.0, 1.1]
        Avals = range(102, 106, 2)
        
        set basis cc-pvdz
        set scf e_convergence=11
        for R in Rvals:
            h2o.R = R
            for A in Avals:
                h2o.A = A
                energy('scf')

   .. code-tab:: py PsiAPI

        h2o = psi4.geometry("""
          O
          H 1 R
          H 1 R 2 A
          """)
        
        Rvals = [0.9, 1.0, 1.1]
        Avals = range(102, 106, 2)
        
        psi4.set_options({"basis": "cc-pvdz",
                          "e_convergence": 11})
        for R in Rvals:
            h2o.R = R
            for A in Avals:
                h2o.A = A
                psi4.energy('scf')

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

The Psithon function ``psi4.driver.p4util.Table`` has been removed,
as the Python ecosystem provides many more flexible alternatives. An
example tabulating a potential energy surface scan for water with Pandas
is shown below

.. tabs::

   .. code-tab:: py PSIthon

        molecule h2o {
          O
          H 1 R
          H 1 R 2 A
        }
        
        Rvals = [0.9, 1.0, 1.1]
        Avals = range(100, 103, 2)
        
        rows = []
        table = []
        
        set basis cc-pvdz
        
        for R in Rvals:
            h2o.R = R
            for A in Avals:
                h2o.A = A
                energy('mp2')
                escf = variable('SCF TOTAL ENERGY')
                edfmp2 = variable('MP2 TOTAL ENERGY')
                escsmp2 = variable('SCS-MP2 TOTAL ENERGY')
                rows.append((R, A))
                table.append([escf, escsmp2, edfmp2])
        
        import pandas as pd
        df = pd.DataFrame(table, columns = ["E(SCF)", "E(SCS)", "E(DFMP2)"], index=rows)
        print(df)

        #                E(SCF)     E(SCS)   E(DFMP2)
        # (0.9, 100) -76.020680 -76.217006 -76.221189
        # (0.9, 102) -76.021305 -76.217439 -76.221605
        # (1.0, 100) -76.021264 -76.224987 -76.228727
        # (1.0, 102) -76.021460 -76.224946 -76.228668
        # (1.1, 100) -75.990195 -76.201891 -76.205087
        # (1.1, 102) -75.990085 -76.201498 -76.204676

   .. code-tab:: py PsiAPI

        h2o = psi4.geometry("""
          O
          H 1 R
          H 1 R 2 A
          """)
        
        Rvals = [0.9, 1.0, 1.1]
        Avals = range(100, 103, 2)
        
        rows = []
        table = []
        
        psi4.set_options({"basis": "cc-pvdz"})
        
        for R in Rvals:
            h2o.R = R
            for A in Avals:
                h2o.A = A
                psi4.energy('mp2')
                escf = psi4.variable('SCF TOTAL ENERGY')
                edfmp2 = psi4.variable('MP2 TOTAL ENERGY')
                escsmp2 = psi4.variable('SCS-MP2 TOTAL ENERGY')
                rows.append((R, A))
                table.append([escf, escsmp2, edfmp2])
        
        import pandas as pd
        df = pd.DataFrame(table, columns = ["E(SCF)", "E(SCS)", "E(DFMP2)"], index=rows)
        print(df)

        #                E(SCF)     E(SCS)   E(DFMP2)
        # (0.9, 100) -76.020680 -76.217006 -76.221189
        # (0.9, 102) -76.021305 -76.217439 -76.221605
        # (1.0, 100) -76.021264 -76.224987 -76.228727
        # (1.0, 102) -76.021460 -76.224946 -76.228668
        # (1.1, 100) -75.990195 -76.201891 -76.205087
        # (1.1, 102) -75.990085 -76.201498 -76.204676


.. _`sec:wrappers`:

Python Wrappers
===============

The Python foundations of the |PSIfour| driver and Psithon syntax permit
many commonly performed post-processing procedures to be integrated into
the |PSIfour| suite.  

As seen in the neon dimer example from the :ref:`tutorial <sec:tutorial>` section,
the :py:func:`~psi4.driver.driver_nbody.nbody` wrapper provides automatic computation of
counterpoise-corrected interaction energies between two molecules.  For
example, 

.. tabs::

   .. code-tab:: py PSIthon

      energy('mp2', bsse_type='cp')

   .. code-tab:: py PsiAPI

      psi4.energy('mp2', bsse_type='cp')

will compute the counterpoise-corrected density-fitted MP2 interaction energy
between two molecules.

|PSIfour| also provides the :py:func:`~psi4.driver.cbs` wrapper,
which automatically computes a complete-basis-set extrapolation (and
automatically sets up the computations with different basis sets required to
do the extrapolation).  For example,

.. tabs::

   .. code-tab:: py PSIthon

      # all equivalent

      energy('mp2', corl_basis='cc-pv[dt]z', corl_scheme='corl_xtpl_helgaker_2')

      energy('mp2/cc-pv[dt]z')

   .. code-tab:: py PsiAPI

      # all equivalent

      psi4.energy('mp2', corl_basis='cc-pv[dt]z', corl_scheme='corl_xtpl_helgaker_2')

      psi4.energy('mp2/cc-pv[dt]z')

will compute a 2-point Helgaker extrapolation of the correlation energy
using the cc-pVDZ and cc-pVTZ basis sets (with method MP2) and add this
extrapolated correlation energy to the Hartree--Fock energy in the
largest basis (cc-pVTZ). :py:func:`~psi4.driver.cbs` can
be configured behind-the-scenes with explicit arguments, as in the
first example, or the convenience syntax of the equivalent second
example can be used.

Another very useful and powerful feature of |PSIfour| is the ability
to compute results on entire databases of molecules at a time,
as provided by the :py:func:`~psi4.driver.wrapper_database.database` wrapper.  For example,

.. tabs::

   .. code-tab:: py PSIthon

      database('mp2', 'S22', cp=1, benchmark='S22B')

   .. code-tab:: py PsiAPI
    
      psi4.wrapper_database.database('mp2', 'S22', cp=1, benchmark='S22B')

will perform DF-MP2 counterpoise-corrected interaction energies
(``cp=1``) on all members of Hobza's S22 database set of van der Waals
dimers, and then compare the results against the S22B benchmark energies.
Built-in databases include S22, A24, HTBH, HBC6, HSG, S22by5, S66, JSCH,
NCB31, S66by8, and NBC10, among others.

These wrapper functions are discussed separately in
:ref:`sec:psithonFunc`.  Note that the options documented for Python
functions are placed as arguments in the command that calls the function,
not in the ``set {...}`` block or with any other ``set`` command.


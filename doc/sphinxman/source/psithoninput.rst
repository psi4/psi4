
.. _`sec:psithoninput`:

==================================
Psithon: Structuring an Input File
==================================

To allow arbitrarily complex computations to be performed, |PSIfour| was built
upon the Python interpreter. However, to make the input syntax simpler, some
pre-processing of the input file is performed before it is interpreted,
resulting in Python syntax that is customized for \PSI, termed Psithon.  In
this section we will describe the essential features of the Psithon language.
|PSIfour| is distributed with an extensive test suite, described in section
:ref:`apdx:testSuite`; the input files for these test cases can be found in the
samples subdirectory of the top-level |PSIfour| source directory, and should
serve as useful examples.

.. _`sec:psirc`:

Scratch Files, the |psirc| File, and Physical Constants
=======================================================

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

.. literalinclude:: ../../../lib/python/physconst.py
   :lines: 3-

The ``psi_`` prefix is to prevent clashes with user-defined variables in
|PSIfour| input files.

Another important use of the |psirc| file is to control the
handling of scratch files.  |PSIfour| has a number of utilities that manage
input and output (I/O) of quantities to and from the hard disk.  Most
quantities, such as molecular integrals, are intermediates that are not of
interest to the user and can be deleted after the computation finishes, but
pertinent details of computations are also written to a checkpoint file and
might be useful in subsequent computations.  All files are sequentially
numbered and are written to ``/tmp``, then deleted at the end of the computation,
unless otherwise instructed by the user.

A Python callable handle to the |PSIfour| I/O management routines is available,
and is called ``psi4_io``.  To instruct the I/O manager to send all files to
another location, say ``/scratch/user``, add the following command to the |psirc|
file (note the trailing "/")::

    psi4_io.set_default_path('/scratch/user/')

For batch jobs running through a queue, it might be more convenient to use an
environmental variable (in this case ``$MYSCRATCH``) to set the scratch directory;
the following code will do that::

    scratch_dir = os.environ.get('MYSCRATCH')
    if scratch_dir:
        psi4_io.set_default_path(scratch_dir + '/')

Individual files can be send to specific locations.  For example, file 32 is
the checkpoint file that the user might want to retain in the working directory
(*i.e.*, where |PSIfour| was launched from) for restart purposes.  This is
accomplished by the commands below::

    psi4_io.set_specific_path(32, './')
    psi4_io.set_specific_retention(32, True)

To circumvent difficulties with running multiple jobs in the same scratch, the
process ID (PID) of the |PSIfour| instance is incorporated into the full file
name; therefore, it is safe to use the same scratch directory for calculations
running simultaneously.

To override any of these defaults for selected jobs, simply place the
appropriate commands from the snippets above in the input file itself.  During
excecution, the |psirc| defaults will be loaded in first, but then the commands
in the input file will be executed.  Executing |PSIfour| with the ``-m`` (for
messy) flag will prevent files being deleted at the end of the run::

    psi4 -m

.. _`sec:moleculespecification`:

Molecule Specification
======================

|PSIfour| has a very flexible input parser that allows the user to provide
geometries as Cartesian coordinates, Z-matrix variables, or a combination of
both. The use of fixed values and variables are supported for both. For
example, the geometry for H\ :sub:`2` can be specified a number of ways, using the
``molecule`` keyword::

    molecule{
      H
      H 1 0.9
    }
    
     or
    
    molecule{
      H
      H 1 r
      r = 0.9
    }
    
     or
    
    molecule{
      H1
      H2 H1 0.9
    }
    
     or
    
    molecule{
      H 0.0 0.0 0.0
      H 0.0 0.0 0.9
    }
    
     or
    
    molecule{
      H 0.0 0.0 0.0
      H 0.0 0.0 r
      r = 0.9
    }
    
     or
    
    molecule{
      H 0.0 0.0 -r
      H 0.0 0.0 r
      r = 0.45
    }

Blank lines are ignored and, unlike regular Python syntax, indentation within
the ``molecule`` block does not matter, although the ``molecule`` keyword itself must
be aligned within the input according to standard Python syntax. For more
examples of geometry specification, see the :source:`samples/mints1/input.dat` input file in the samples
folder. It is also possible to mix Cartesian and Z-matrix geometry
specifications, as demonstrated in the :source:`samples/mints4/input.dat` sample input file.

.. _`sec:multiplemolecules`:

Multiple Molecules
^^^^^^^^^^^^^^^^^^

To facilitate more elaborate computations, it is possible to provide a name for
each molecule, and tell |PSIfour| which one should be used in a given
calculation. For example, consider the following input file::

    molecule h2{
      H
      H 1 0.9
    }
    
    
    set basis cc-pvdz
    set reference rhf
    energy('scf')
    
    molecule h{
      H
    }
    
    set basis cc-pvdz
    set reference uhf
    energy('scf')

Here, two separate jobs are performed on two different molecules; the first is
performed on H\ :sub:`2`, while the second is for H atom. The last molecule to be
specified is the "active" molecule by default. To explicitly activate a named
molecule, the activate keyword is provided. Using this keyword, the above input
file can be equivalently written as follows::

    molecule h2{
      H
      H 1 0.9
    }
    
    molecule h{
      H
    }
    
    activate(h2)
    set basis cc-pvdz
    set reference rhf
    energy('scf')
    
    activate(h)
    set basis cc-pvdz
    set reference uhf
    energy('scf')

Note that whenever the molecule is changed, the basis set must be specified
again. The following section provides more details about the job control
keywords used in the above examples.

.. _`sec:moleculekeywords`:

Molecule Keywords
^^^^^^^^^^^^^^^^^

In addition to specifying the geometry, additional information can be provided
in the ``molecule`` block. If two integers are encountered on any line of the
``molecule`` block, they are interpreted as the molecular charge and multiplicity
(:math:`2 \times M_s + 1`), respectively.  The symmetry can be specified by a line reading
``symmetry symbol``, where ``symbol`` is
the Sch\"onflies symbol of the (Abelian) point group to use for the
computation; see Sec. :ref:`sec:symmetry` for more details. This need not be
specified, as the molecular symmetry is automatically detected by |PSIfour|.
Certain computations require that the molecule is not reoriented; this can be
achieved by adding either ``no_reorient`` or ``noreorient``. By default,
\AA ngstr\"om units are used; this is changed by adding a line that reads
``units spec``, where ``spec`` is one of ``ang``,
``angstrom``, ``a.u.``, ``au``, or ``bohr``.

.. _`sec:pubchem`:

Geometries from the `PubChem <http://pubchem.ncbi.nlm.nih.gov/>`_ Database
==========================================================================

Obtaining rough starting guess geometries can be burdensome.  The Z-matrix
coordinate system was designed to provide chemists with an intuitive method for
guessing structures in terms of bond lengths and angles.  While Z-matrix input is
intuitive for small molecules with few degrees of freedom, it quickly becomes
laborious as the system size grows.  To obtain a reasonable starting guess
geometry, |PSIfour| can take a chemical name as input; this is then used
to attempt to retrieve Cartesian coordinates from the [PubChem]_ database.

For example, to run a computation on benzene, we can use the following molecule specification::

    molecule benzene {
        pubchem:benzene
    }

If the computer is connected to the internet, the above code will instruct
|PSIfour| to search PubChem for a starting structure.  The search is actually
performed for compounds whose name *contains* "benzene", so multiple
entries will be returned.  If the name provided ("benzene" in the above
example) exactly matches one of the results, that entry will be used.  If no
exact match is found the results, along with a unique chemical identifier
(CID), are printed to the output file, prompting the user to provide a more
specific name.  For example, if we know that we want to run a computation on a
compound whose name(s) contain "benzene", but we're not sure of the exact IUPAC
name, the following input can be used::

    molecule benzene {
        pubchem:benzene*
    }

Appending the "*" prevents an exact match from being found and, at the time
of writing, the following results are displayed in the output file::

     Chemical ID     IUPAC Name
              241   benzene
             7371   benzenesulfonic acid
            91526   benzenesulfonate
              244   phenylmethanol
              727   1,2,3,4,5,6-hexachlorocyclohexane
              240   benzaldehyde
            65723   benzenesulfonohydrazide
            74296   N-phenylbenzenesulfonamide
              289   benzene-1,2-diol
              243   benzoic acid
             7370   benzenesulfonamide
           636822   1,2,4-trimethoxy-5-[(E)-prop-1-enyl]benzene
             7369   benzenesulfonyl chloride
            12932   N-[2-di(propan-2-yloxy)phosphinothioylsulfanylethyl]benzenesulfonamide
             7505   benzonitrile
            78438   N-[anilino(phenyl)phosphoryl]aniline
            12581   3-phenylpropanenitrile
           517327   sodium benzenesulfonate
           637563   1-methoxy-4-[(E)-prop-1-enyl]benzene
           252325   [(E)-prop-1-enyl]benzene

Note that some of these results do not contain the string "benzene"; these
compounds have synonyms containing that text.  We can now replace the
"benzene*" in the input file with one of the above compounds using either the
IUPAC name or the CID provided in the list, *viz*::

    molecule benzene {
        pubchem:637563
    }
    
     or
    
    molecule benzene {
        pubchem:1-methoxy-4-[(E)-prop-1-enyl]benzene
    }

Some of the structures in the database are quite loosely optimized and do not
have the correct symmetry.  Before starting the computation, |PSIfour| will
check to see if the molecule is close to having each of the possible
symmetries, and will adjust the structure accordingly so that the maximum
symmetry is utilized.

The standard keywords, described in Sec. :ref:`sec:moleculeKeywords`, can be
used in conjuction to specify charge, multiplicity, symmetry to use, *etc.* .

.. _`sec:symmetry`:

Symmetry
========

For efficiency, |PSIfour| can utilize the largest Abelian subgroup of the full
point group of the molecule.  Concomitantly a number of quantities, such as
:term:`SOCC` and :term:`DOCC`, are arrays whose entries pertain to irreducible
representations (irreps) of the molecular point group.  Ordering of irreps
follows the convention used in Cotton's :title:`Chemical Applications of Group
Theory`, as detailed in Table :ref:`Irreps <table:irrepOrdering>`.  We refer to this
convention as "Cotton Ordering" hereafter.

.. _`table:irrepOrdering`:

.. table:: Ordering of irreducible representations (irreps) used in |PSIfour|

    +----------------+-------------+----------------+----------------+----------------+-------------+----------------+----------------+----------------+
    | Point Group    | Irrep Order                                                                                                                     |
    +================+=============+================+================+================+=============+================+================+================+
    | :math:`C_1`    | :math:`A`   |                |                |                |             |                |                |                |  
    +----------------+-------------+----------------+----------------+----------------+-------------+----------------+----------------+----------------+
    | :math:`C_i`    | :math:`A_g` | :math:`A_u`    |                |                |             |                |                |                |
    +----------------+-------------+----------------+----------------+----------------+-------------+----------------+----------------+----------------+
    | :math:`C_2`    | :math:`A`   | :math:`B`      |                |                |             |                |                |                |
    +----------------+-------------+----------------+----------------+----------------+-------------+----------------+----------------+----------------+
    | :math:`C_s`    | :math:`A'`  | :math:`A''`    |                |                |             |                |                |                |
    +----------------+-------------+----------------+----------------+----------------+-------------+----------------+----------------+----------------+
    | :math:`D_2`    | :math:`A`   | :math:`B_1`    | :math:`B_2`    | :math:`B_3`    |             |                |                |                |
    +----------------+-------------+----------------+----------------+----------------+-------------+----------------+----------------+----------------+
    | :math:`C_{2v}` | :math:`A_1` | :math:`A_2`    | :math:`B_1`    | :math:`B_2`    |             |                |                |                |
    +----------------+-------------+----------------+----------------+----------------+-------------+----------------+----------------+----------------+
    | :math:`C_{2h}` | :math:`A_g` | :math:`B_g`    | :math:`A_u`    | :math:`B_u`    |             |                |                |                |
    +----------------+-------------+----------------+----------------+----------------+-------------+----------------+----------------+----------------+
    | :math:`D_{2h}` | :math:`A_g` | :math:`B_{1g}` | :math:`B_{2g}` | :math:`B_{3g}` | :math:`A_u` | :math:`B_{1u}` | :math:`B_{2u}` | :math:`B_{3u}` |
    +----------------+-------------+----------------+----------------+----------------+-------------+----------------+----------------+----------------+

For example, water (:math:`C_{2v}` symmetry) has 3 doubly occupied :math:`A_1`
orbitals, as well as 1 each of :math:`B_1` and :math:`B_2` symmetry; the
corresponding :term:`DOCC` array is therefore::

    DOCC = [3, 0, 1, 1]

Although |PSIfour| will detect the symmetry automatically, and use the largest
possible Abelian subgroup, the user might want to run in a lower point group.
To do this the ``symmetry`` keyword can be used when inputting the molecule
(see Sec. :ref:`sec:moleculeSpecification`).  In most cases the standard
Sch\"onflies symbol (one of ``c1``, ``c2``, ``ci``, ``cs``, ``d2``,
``c2h``, ``c2v``, ``d2h`` will suffice.
For certain computations, the user might want to specify which particular
subgroup is to be used by appending a unique axis specifier.  For example when
running a computation on a molecule with :math:`D_{2h}` symmetry in :math:`C_{2v}`, the
:math:`C_2` axis can be chosen as either the :math:`x`, the :math:`y`, or the :math:`z`; these can
be specified by requesing the symmetry as ``c2vx``, ``c2vy``, or ``c2vz``, respectively.
Likewise the ``c2x``, ``c2y``, ``c2z``, ``c2hx``, ``c2hy``, and ``c2hz``
labels are valid.  For :math:`C_s` symmetry the labels ``csx``, ``csy``, and
``csz`` request the :math:`yz`, :math:`xz`, and :math:`xy` planes be used as the mirror plane,
respectively.  If no unique axis is specified, |PSIfour| will choose an appropriate
subgroup.

.. _`sec:fragments`:

Non-Covalently Bonded Molecule Fragments
========================================

|PSIfour| has an extensive range of tools for treating non-covalent
intermolecular forces, including counterpoise corrections and symmetry adapted
perturbation theory methods. These require the definition of which fragments
are interacting within the complex. |PSIfour| provides a very simple mechanism
for doing so; simply define the complex's geometry using the standard
Cartesian, Z-matrix, or mixture thereof, specifications and then place two
dashes between nonbonded fragements. For example, to study the interaction
energy of ethane and ethyne molecules, we can use the following molecule
block::

    molecule{
      0 1
      C  0.000000 -0.667578  -2.124659
      C  0.000000  0.667578  -2.124659
      H  0.923621 -1.232253  -2.126185
      H -0.923621 -1.232253  -2.126185
      H -0.923621  1.232253  -2.126185
      H  0.923621  1.232253  -2.126185
      --
      0 1
      C 0.000000 0.000000 2.900503
      C 0.000000 0.000000 1.693240
      H 0.000000 0.000000 0.627352
      H 0.000000 0.000000 3.963929
    }

In this case, the charge and multiplicity of each interacting fragment is
explicitly specified. If the charge and multiplicity are specified for the
first fragment, it is assumed to be the same for all fragments. When
considering interacting fragments, the overall charge is simply the sum of all
fragment charges, and any unpaired electrons are assumed to be coupled to
yield the highest possible :math:`M_s` value.

.. _`sec:jobcontrol`:

Job Control
===========

|PSIfour| comprises a number of modules, written in C++, that each perform
specific tasks and are callable directly from the Python front end. Each module
recognizes specific keywords in the input file, detailed in Appendix :ref:`apdx:options_c_module`, which
control its function. The keywords can be made global, or scoped to apply to
certain specific modules. The following examples demonstrate some of the ways
that global keywords can be specified::

    set globals basis cc-pVDZ
    
     or
    
    set basis cc-pVDZ
    
     or
    
    set globals basis = cc-pVDZ
    
     or
    
    set basis = cc-pVDZ
    
     or
    
    set globals{
      basis cc-pVDZ
    }
    
     or
    
    set{
      basis cc-pVDZ
    }
    
     or
    
    set{
      basis = cc-pVDZ
    }

Note the lack of quotes around ``cc-pVDZ``, even though it is a string. The
Psithon preprocessor automatically wraps any string values in ``set`` commands in
strings. The last three examples provide a more convenient way for specifying
multiple keywords::

    set{
      basis = cc-pVDZ
      print = 1
      reference = rhf
    }

For arguments that require an array input, standard Python list syntax should
be used, *viz.*::

    set{
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

.. _`sec:psithonbasissets`:

Assigning Basis Sets
====================

While the above syntax will suffice for specifying basis sets in most cases,
the user may need to assign basis sets to specific atoms.  To achieve this, a
``basis`` block can be used.  We use a snippet from the :source:`samples/mints2/input.dat` sample
input file, which performs a benzene SCF computation, to demonstrate this
feature. ::

    basis {
       assign DZ
       assign C 3-21G
       assign H1 sto-3g
       assign C1 sto-3g
    }

The first line in this block assigns the DZ basis set to all atoms.  The next
line then assigns 3-21G to all carbon atoms, leaving the hydrogens with the DZ
basis set.  On the third line, the hydrogen atoms which have been specifically
labelled as ``H1`` are given the STO-3G basis set, leaving the unlabelled hydrogen
atoms with the DZ basis set.  Likewise, the fourth line assigns the STO-3G
basis set to just the carbon atoms labelled ``C1``.  This bizzare example was
constructed to demonstrate the syntax, but the flexibility of the basis set
specification is advantageous, for example, when selectivily omitting diffuse
functions to make computations more tractable.

In the above example the basis sets have been assigned asymmetrically, reducing
the effective symmetry from :math:`D_{6h}` to :math:`C_{2v}`; |PSIfour| will detect this
automatically and run in the appropriate point group.  The same syntax can be
used to specify basis sets other than that used to define orbitals.  For
example, ::

    set df_basis_mp2 cc-pvdz-ri
    
     or
    
    basis {
       assign cc-pVDZ-RI df_basis_mp2
    }

are both equivalent ways to set the auxiliary basis set for density fitted MP2
computations.  To assign the aug-cc-pVDZ-RI to carbon atoms, the following
command is used::

    basis {
       assign C aug-cc-pVDZ-RI df_basis_mp2
    }

When Dunning's correlation consistent basis sets (cc-pV*X*Z), and core-valence
and diffuse variants thereof, are being used the SCF and DF-MP2 codes will
chose the appropriate auxilliary basis set automatically, unless instructed
otherwise by setting theauxiliary basis set in the input.  Finally, we note
that the ``basis`` block may also be used for defining basis sets, as
detailed in Sec. :ref:`sec:basissetSpecification`.

.. _`sec:memory`:

Memory Specification
====================

By default, |PSIfour| assumes that 256 Mb of memory are available. While this is
enough for many computations, many of the algorithms will perform better if
more is available. To specify memory, the ``memory`` keyword should be used. The following
lines are all equivalent methods for specifying that 2 Gb of RAM is available
to |PSIfour|::

    memory 2 Gb
    
     or
    
    memory 2000 Mb
    
     or
    
    memory 2000000 Kb

One convenient way to override the |PSIfour| default memory is to place a memory
command in the |psirc| file, as detailed in Sec. :ref:`sec:psirc`.

.. _`sec:threading`:

Threading
=========

Most new modules in |PSIfour| are designed to run efficiently on SMP architectures
via application of several thread models. The de facto standard for |PSIfour|
involves using threaded BLAS/LAPACK (particularly Intel's excellent MKL package)
for most tensor-like operations, OpenMP for more general operations, and Boost
Threads for some special-case operations. Note: Using OpenMP alone is a really
bad idea. The developers make little to no effort to explicitly parallelize
operations which are already easily threaded by MKL or other threaded BLAS. Less
than 20% of the threaded code in |PSIfour| uses OpenMP, the rest is handled by
parallel DGEMM and other library routines. From this point forward, it is
assumed that you have compiled PSI4 with OpenMP and MKL (Note that it is
possible to use g++ or another compiler and yet still link against MKL).

Control of threading in |PSIfour| can be accomplished at a variety of levels,
ranging from global environment variables to direct control of thread count in
the input file, to even directives specific to each model. This hierarchy is
explained below. Note that each deeper level trumps all previous levels.

.. rubric:: (1) OpenMP/MKL Environment Variables

The easiest/least visible way to thread |PSIfour| is to set the standard OpenMP/MKL
environment variables. For instance, in tcsh::

    setenv OMP_NUM_THREADS 4
    setenv MKL_NUM_THREADS 4

|PSIfour| then detects these value via the API routines in ``<omp.h>`` and
``<mkl.h>``, and runs all applicable code with 4 threads. These environment
variables are typically defined in a ``.tcshrc`` or ``.bashrc``.

.. rubric:: (2) The -n Command Line Flag

To change the number of threads at runtime, the ``-n`` flag may be used. An
example is::

    psi4 -i input.dat -o output.dat -n 4

which will run on four threads.

.. rubric:: (3) Setting Thread Numbers in an Input

For more explicit control, the Process::environment class in |PSIfour| can
override the number of threads set by environment variables. This functionality
is accessed via the :py:func:`util.set_num_threads` Psithon function, which controls
both MKL and OpenMP thread numbers. The number of threads may be changed
multiple times in a |PSIfour| input file. An example input for this feature is::

    # A bit small-ish, but you get the idea
    molecule h2o {
    0 1
    O
    H 1 1.0
    H 1 1.0 2 90.0
    }
    
    set scf {
    basis cc-pvdz
    scf_type df
    }

    # Run from 1 to 4 threads, for instance, to record timings
    for nthread in range(1,5):
        set_num_threads(nthread)
        energy('scf')

.. rubric:: (4) Method-Specific Control

Even more control is possible in certain circumstances. For instance, the
threaded generation of AO density-fitted integrals involves a memory requirement
proportional to the number of threads. This requirement may exceed the total
memory of a small-memory node if all threads are involved in the generation of
these integrals. For general DF algorithms, the user may specify::

    set MODULE_NAME df_ints_num_threads n

to explicitly control the number of threads used for integral formation. Setting
this variable to 0 (the default) uses the number of threads specified by the
:py:func:`util.set_num_threads` Psithon method or the default environmental variables.

.. _`sec:psivariables`:

Return Values and PSI Variables
===============================

To harness the power of Python, |PSIfour| makes the most pertinent results of
each computation are made available to the Python interpreter for
post-processing. To demonstrate, we can embellish the previous example of H\ :sub:`2`
and H atom::

    molecule h2{
      H
      H 1 0.9
    }
    
    set basis cc-pvdz
    set reference rhf
    h2_energy = energy('scf')
    
    molecule h{
      H
    }
    
    set basis cc-pvdz
    set reference uhf
    h_energy = energy('scf')
    
    D_e = psi_hartree2kcalmol*(2*h_energy - h2_energy)
    print"De=%f"%D_e

The :py:func:`driver.energy` function returns the final result of the computation, which we
assign to a Python variable. The two energies are then converted to a
dissociation energy and printed to the output file using standard Python
notation. Sometimes there are multiple quantities of interest; these can be
accessed through the ``get_variable()`` function. For example, after performing a
density fitted MP2 computation, both the spin component scaled energy and the
unscaled MP2 energy are made available::

    e_mp2=get_variable('DF-MP2 TOTAL ENERGY')
    e_scs_mp2=get_variable('SCS-DF-MP2 TOTAL ENERGY')

Each module and the Python driver set PSI variables over the course of
a calculation.  The values for all can be printed in the output file
with the input file command ``print_variables()``. Note that
PSI variables accumulate over a |PSIfour| instance and are not cleared by
``clean()``. So if you run in a single input file a STO-3G FCI
followed by a aug-cc-pVQZ SCF followed by a ``print_variables()``
command, the last will include both {\tt SCF TOTAL ENERGY} and
{\tt FCI TOTAL ENERGY}. Don't get excited that you got a high-quality calculation
cheaply. Refer to Appendix :ref:`apdx:psiVariables` for a listing of the
variables set by each module.

.. _`sec:loops`:

Loops
=====

Python provides many control structures, which can be used within |PSIfour|
input files. For example, to loop over three basis sets, the following code can
be used::

    basis_sets=["cc-pVDZ","cc-pVTZ","cc-pVQZ"]
    for basis_set in basis_sets:
        set basis = $basis_set
        energy('scf')

The declaration of ``basis_set`` is completely standard Python, as is the next
line, which iterates over the list. However, because the Psithon preprocessor
wraps strings in quotes by default, we have to tell it that ``basis_set`` is a
Python variable, not a string, by prefixing it with a dollar sign. The geometry
specification supports delayed initialization of variable, which permits
potential energy scans. As an example, we can scan both the angle and bond
length in water::

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

.. _`sec:resultstables`:

Tables of Results
=================

The results of computations can be compactly tabulated with the :py:func:`text.Table()` Psithon
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
            energy('dfmp2')
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
the |PSIfour| suite.
Among these are automated computations of interaction energies through
:py:func:`wrappers.cp`, of a model chemistry applied to a database of systems through
:py:func:`wrappers.database`, and of several model chemistries together approximating greater
accuracy through :py:func:`wrappers.cbs`. 
These are discussed separately in Sec. :ref:`sec:psithonfunc`.
Note that the options documented for Python functions are placed as arguments
in the command that calls the function 
not in the ``set globals`` block or with any other ``set`` command.



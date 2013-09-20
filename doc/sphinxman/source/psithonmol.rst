
.. include:: autodoc_abbr_options_c.rst

.. index:: molecule; specification
.. _`sec:moleculeSpecification`:

===================================
Molecule and Geometry Specification
===================================

Coordinates
===========

|PSIfour| has a very flexible input parser that allows the user to provide
geometries as Cartesian coordinates, Z-matrix variables, or a combination of
both. The use of fixed values and variables are supported for both. For
example, the geometry for H\ :sub:`2` can be specified a number of ways, using the
:samp:`molecule {optional_molecule_name} \\{...\\}` block. ::

    molecule {
      H
      H 1 0.9
    }
    
or ::
    
    molecule {
      H
      H 1 r
      r = 0.9
    }
    
or ::
    
    molecule {
      H1
      H2 H1 0.9
    }
    
or ::
    
    molecule {
      H 0.0 0.0 0.0
      H 0.0 0.0 0.9
    }
    
or ::
    
    molecule {
      H 0.0 0.0 0.0
      H 0.0 0.0 r
      r = 0.9
    }
    
or ::
    
    molecule {
      H 0.0 0.0 -r
      H 0.0 0.0 r
      r = 0.45
    }

Blank lines are ignored and, unlike regular Python syntax, indentation within
the molecule block does not matter, although the ``molecule`` keyword itself must
be aligned within the input according to standard Python syntax. For more
examples of geometry specification, see the :srcsample:`mints1` input file in the samples
folder. It is also possible to mix Cartesian and Z-matrix geometry
specifications, as demonstrated in the :srcsample:`mints4` and
:srcsample:`mints6` sample input files.  For example, consider the following
geometry specification, taken from the :srcsample:`mints6` input::

    molecule alanine {
        N           -1.527107413251     0.745960643462     0.766603000356
        C           -0.075844098953     0.811790225041     0.711418672248
        C            0.503195220163    -0.247849447550    -0.215671574613
        O           -0.351261319421    -0.748978309671    -1.089590304723
        O            1.639498336738    -0.571249748886    -0.174705953194
        H           -1.207655674855    -0.365913941094    -0.918035522052
        # First, remove the H from the alpha carbon.  This line could be deleted
        # and is only included for completeness
        #H            0.429560656538     0.717651915252     1.673774709694
        # Now patch it, using a Z Matrix specification.  This patch can be applied
        # anywhere in the coord specification, as long as it appears lower than
        # the atoms referenced, as is usual for Z-Matrices
        C  2  rCC   3  aCCC   1  dCCCN
        H  7  rCH1  2  aHCC1  3  dHCCC1
        H  7  rCH2  2  aHCC2  3  dHCCC2
        H  7  rCH3  2  aHCC3  3  dHCCC3
        H            0.221781602033     1.772570540211     0.286988509018
        H           -1.833601608592     0.108401996052     1.481873213172
        H           -1.925572581453     1.640882152784     0.986471814808
    
        aCCC = 108.0
        rCC = 1.4
        dCCCN = 120
        rCH1 = 1.08
        rCH2 = 1.08
        rCH3 = 1.08
        aHCC1 = 109.0
        aHCC2 = 109.0
        aHCC3 = 109.0
        dHCCC1 = 0.0
        dHCCC2 = 120.0
        dHCCC3 = 240.0
    }

Here, we remove the hydrogen from the alpha carbon of glycine and replace it
with a methyl group.  Applying this patch using Cartesian coordinates is
difficult, because it depends on the orientation of the existing glycine unit.
In this example, we use Z-Matrix coordinates to define the methyl group, and
define the orientation in terms of the existing glycine Cartesian coordinates
which is much easier to visualize than the corresponding Cartesian-only
approach.

.. index:: molecule; multiple in input file
.. _`sec:multipleMolecules`:

.. index::
   triple: setting; keywords; molecule
   pair: molecule; charge
   pair: molecule; multiplicity
   pair: molecule; symmetry
   pair: molecule; no_reorient
   pair: molecule; units
.. _`sec:moleculeKeywords`:

Molecule Keywords
=================

In addition to specifying the geometry, additional information can be
provided in the molecule block :samp:`molecule {optional_molecule_name} \\{...\\}`.

**Charge & Multiplicity**
   If two integers :samp:`{charge} {multiplicity}` are encountered on any
   line of the molecule block, they are interpreted as the molecular charge
   and multiplicity (:math:`2 M_s + 1`), respectively. For multi-fragment 
   complexes, each fragment can have a :samp:`{charge} {multiplicity}` line.

**Units**
   By default, |Angstrom| units are used; this is changed by adding
   a line that reads :samp:`units {spec}`, where :samp:`{spec}` is one
   of ``ang``, ``angstrom``, ``a.u.``, ``au``, or ``bohr``.

**Orientation**
   Certain computations require that the molecule is not reoriented. This 
   can be achieved by adding either ``no_reorient`` or ``noreorient``. 
   To prevent even recentering of the molecule, add ``no_com`` or ``nocom``.

**PubChem**
   A line reading :samp:`pubchem:{mol}` fetches the geometry for molecule
   :samp:`{mol}` from the PubChem database, where :samp:`{mol}` is either
   the IUPAC molecule name or the CID number. See :ref:`sec:pubchem` for
   details.

**Symmetry**
   The symmetry can be specified by a line reading :samp:`symmetry
   {symbol}`, where :samp:`{symbol}` is the Sch\ |o_dots|\ nflies symbol
   of the (Abelian) point group to use for the computation, one of one of
   ``c1``, ``c2``, ``ci``, ``cs``, ``d2``, ``c2h``, ``c2v``, or ``d2h``.
   This need not be specified, as the molecular symmetry is automatically
   detected by |PSIfour|. See :ref:`sec:symmetry` for details.

**Fragments**
   A line reading ``--`` is interpreted as the separator between two non-covalently 
   bound molecular fragments. See :ref:`sec:fragments` for details.

Multiple Molecules
==================

To facilitate more elaborate computations, it is possible to provide a name for
each molecule and tell |PSIfour| which one should be used in a given
calculation. For example, consider the following input file::

    molecule h2 {
      H
      H 1 0.9
    }
    
    set basis cc-pvdz
    set reference rhf
    energy('scf')
    
    clean()

    molecule h {
      H
    }
    
    set basis cc-pvdz
    set reference uhf
    energy('scf')

Here, two separate jobs are performed on two different molecules; the first is
performed on H\ :sub:`2`, while the second is for H atom. The last molecule to be
specified is the "active" molecule by default. To explicitly activate a named
molecule, the activate command is provided. With it, the above input
file can be equivalently written as follows::

    molecule h2 {
      H
      H 1 0.9
    }
    
    molecule h {
      H
    }
    
    activate(h2)
    set basis cc-pvdz
    set reference rhf
    energy('scf')
    
    clean()

    activate(h)
    set basis cc-pvdz
    set reference uhf
    energy('scf')

Note that whenever the molecule is changed, the basis set must be specified
again. :ref:`sec:jobControl` provides more details about the job control
and calculation keywords used in the above examples.

.. index:: 
   single: Ghost Atoms
   single: molecule; ghost
.. _`sec:ghosts`:

Ghost Atoms
===========

While many common computations, particularly SAPT and counterpoise corrections, can
be greatly simplified using the notation described in :ref:`sec:fragments`,
manual specification of ghost atoms is sometimes required.  Either ::

    molecule he2 {
        He
        Gh(He) 1 2.0
    }

or ::

    molecule he2 {
        He
        @He 1 2.0
    }

will generate a helium dimer with the second atom ghosted, *i.e.*, possessing
basis functions but no electrons or nuclear charge.  See :srcsample:`dfmp2_1`
and :srcsample:`ghosts` for a demonstration of both mechanisms for specifying
ghost atoms.

.. index:: 
   single: PubChem
   single: molecule; PubChem
.. _`sec:pubchem`:

`PubChem <http://pubchem.ncbi.nlm.nih.gov/>`_ Database
======================================================

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
    
or ::
    
    molecule benzene {
        pubchem:1-methoxy-4-[(E)-prop-1-enyl]benzene
    }

Some of the structures in the database are quite loosely optimized and do not
have the correct symmetry.  Before starting the computation, |PSIfour| will
check to see if the molecule is close to having each of the possible
symmetries, and will adjust the structure accordingly so that the maximum
symmetry is utilized.

The standard keywords, described in :ref:`sec:moleculeKeywords`, can be
used in conjuction to specify charge, multiplicity, symmetry to use, *etc.* .

.. index:: symmetry, Cotton-ordering
.. _`sec:symmetry`:

Symmetry
========

For efficiency, |PSIfour| can utilize the largest Abelian subgroup of the full
point group of the molecule. Concomitantly, a number of quantities, such as
|globals__socc| and |globals__docc|, are arrays whose entries pertain to irreducible
representations (irreps) of the molecular point group.  Ordering of irreps
follows the convention used in Cotton's :title:`Chemical Applications of Group
Theory`, as detailed in Table :ref:`Irreps <table:irrepOrdering>`.  We refer to this
convention as "Cotton Ordering" hereafter.

.. _`table:irrepOrdering`:

.. table:: Ordering of irreducible representations (irreps) used in |PSIfour|

    +----------------+-------------+----------------+----------------+----------------+-------------+----------------+----------------+----------------+
    | Point Group    |      1      |       2        |       3        |      4         |     5       |        6       |       7        |       8        |  
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

For example, water (:math:`C_{2v}` symmetry) has three doubly occupied :math:`A_1`
orbitals, as well as one each of :math:`B_1` and :math:`B_2` symmetry; the
corresponding |globals__docc| array is therefore::

    DOCC = [3, 0, 1, 1]

Although |PSIfour| will detect the symmetry automatically, and use the largest
possible Abelian subgroup, the user might want to run in a lower point group.
To do this the molecule keyword :samp:`symmetry {symbol}` can be used 
(see :ref:`sec:moleculeKeywords`).  In most cases the standard
Schönflies symbol (one of ``c1``, ``c2``, ``ci``, ``cs``, ``d2``,
``c2h``, ``c2v``, ``d2h`` will suffice for :samp:`{symbol}`.
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

Certain types of finite difference computations, such as numerical vibrational
frequencies, might lower the symmetry of the molecule.  When this happens
symmetry-dependent arrays, such as |globals__socc|, are automatically remapped
to the lower symmetry.  For example, if we were to investigate the :math:`^2B_1`
state of water cation, we can specify

    SOCC = [0, 0, 1, 0]

in the input file.  If any ensuing computations lower the symmetry, the above
array will be appropriately remapped.  For example, reducing the symmetry to
:math:`C_s` (with the molecular plane defining the mirror plane), the above
array will be automatically interpreted as:

    SOCC = [0, 1]

Some caution is required, however.  The :math:`^2A_1` state can be obtained with
the

    SOCC = [1, 0, 0, 0]

specification, which would become

    SOCC = [1, 0]

under the above-mentioned reduction in symmetry.  The :math:`^2B_2` state,
whose singly-occupied orbitals are

    SOCC = [0, 0, 0, 1]

would be mapped to 

    SOCC = [1, 0]

which is the same occupation as the :math:`^2A_1` state.  In this case, the
:math:`^2A_1` state is lower in energy, and is not problematic.  The distorted
geometries for the :math:`^2B_2` state are excited states that are subject to
variational collapse.  One way to obtain reliable energies for these states is
to use a multi-state method; in this case it's easier to run the entire
computation in the lowest symmetry needed during the finite difference
procedure.

.. index:: molecule; multiple fragments
.. _`sec:fragments`:

Non-Covalently Bonded Molecule Fragments
========================================

|PSIfour| has an extensive range of tools for treating non-covalent
intermolecular forces, including counterpoise corrections and symmetry adapted
perturbation theory methods. These require the definition of which fragments
are interacting within the complex. |PSIfour| provides a very simple mechanism
for doing so: simply define the complex's geometry using the standard
Cartesian, Z-matrix, or mixture thereof, specifications and then place two
dashes between nonbonded fragements. For example, to study the interaction
energy of ethane and ethyne molecules, we can use the following molecule
block::

    molecule eneyne {
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

Having defined a molecule containing fragments like ``eneyne`` above, it
is a simple matter to perform calculations on only a subset of the
fragments. For instance, the commands below run a scf first on the ethene
fragment alone (``extract_subsets(1)`` pulls out fragment 1 as Real atoms
and discards remaining fragments) and next on the ethene fragment with the
ethyne fragment ghosted (``extract_subsets(1,2)`` pulls out fragment 1 as
Real atoms and sets fragment 2 as Ghost atoms). For beyond bimolecular
complexes, arrays can be used, e.g. ``extract_subsets(2,[1,3])``::

   mA = eneyne.extract_subsets(1)
   energy('scf')
   
   clean()
   
   mAcp = eneyne.extract_subsets(1,2)
   energy('scf')

If the molecule contains fragments but is not conveniently ordered for the
``--`` marker, the auto_fragment function can be applied, as shoown in
:srcsample:`pywrap-basis`, to return as active molecule the previous
active molecule, only fragmented.


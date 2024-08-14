.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2024 The Psi4 Developers.
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

.. index::
   single: geometry optimization, optimization

.. _`sec:optking`:

Geometry Optimization
=====================

.. codeauthor:: Rollin A. King and Alexander G. Heide
.. sectionauthor:: Rollin A. King, Alexander G. Heide, and Lori A. Burns

*Module:* :ref:`Keywords <apdx:optking>`, `OPTKING <https://github.com/psi-rking/optking>`_

|PSIfour| carries out molecular optimizations using a Python module called
optking.  The optking program takes as input nuclear gradients and,
optionally, nuclear second derivatives |w---w| both in Cartesian coordinates.
The default minimization algorithm employs an empirical model Hessian,
redundant internal coordinates, an RFO step with trust radius scaling, and the BFGS Hessian update.

The principal literature references include the introduction of redundant
internal coordinates by Peng et al. [Peng:1996:49]_.
The general approach employed in this code
is similar to the "model Hessian plus RF method" described and tested by Bakken and
Helgaker [Bakken:2002:9160]_. However, for separated
fragments, we have chosen not to employ their "extra-redundant" coordinates.

The internal coordinates are generated automatically based on an assumed bond
connectivity.  The connectivity is determined by testing if the interatomic
distance is less than the sum of atomic radii times the value of
|optking__covalent_connect|. If the user finds that some
connectivity is lacking by default, then this value may be increased.

.. warning:: The selection of a Z-matrix input, and in particular the inclusion
   of dummy atoms, has no effect on the behavior of the optimizer, which begins
   from a Cartesian representation of the system.

.. _DimerIntro:

Presently, by default, separate fragments are bonded by the
nearest atoms, and the whole system is treated as if it were part of one
molecule. However, with the option |optking__frag_mode|, fragments
may instead be related by a minimal set of interfragment coordinates defined by
reference points within each fragment.  The reference points can be atomic
positions (current default) or linear combinations of atomic positions
(automatic use of principal axes is under development).
These `dimer coordinates` can be directly specified through |optking__interfrag_coords|)
See `here <DimerSection_>` for two examples of their use.

Basic Keywords
^^^^^^^^^^^^^^

.. include:: autodir_options_c/optking__opt_type.rst
.. include:: autodir_options_c/optking__step_type.rst
.. include:: autodir_options_c/optking__geom_maxiter.rst
.. include:: autodir_options_c/optking__g_convergence.rst
.. include:: autodir_options_c/optking__full_hess_every.rst

.. index:: geometry optimization; minima
.. _`sec:optkingExamples`:

Optimizing Minima
^^^^^^^^^^^^^^^^^

First, define the molecule and basis in the input. ::

   molecule h2o {
     O
     H 1 1.0
     H 1 1.0 2 105.0
   }
    
   set basis dz

Then the following are examples of various types of calculations that can be completed.

* Optimize a geometry using default methods (RFO step)::

   optimize('scf')

* Optimize using Newton-Raphson steps instead of RFO steps::

   set step_type nr
   optimize('scf')

* Optimize using finite differences of energies instead of gradients::

   optimize('scf', dertype='energy')

* Optimize while limiting the initial step size to 0.1 au::

   set intrafrag_step_limit 0.1
   optimize('scf')

* Optimize while always limiting the step size to 0.1 au:

.. code-block:: none

   set {
     intrafrag_step_limit     0.1
     intrafrag_step_limit_min 0.1
     intrafrag_step_limit_max 0.1
   }

   optimize('scf')

* Optimize while calculating the Hessian at every step:

.. code-block:: none

   set full_hess_every 1
   optimize('scf')

.. code-block:: none
    
    import optking
    

Hessian
^^^^^^^

If Cartesian second derivatives are available, optking can read them
and transform them into internal coordinates to make an initial Hessian in
internal coordinates.  Otherwise, several empirical Hessians are available,
including those of Schlegel [Schlegel:1984:333]_ and Fischer and Almlof
[Fischer:1992:9770]_.
Either of these or a simple diagonal Hessian may be selected using the 
|optking__intrafrag_hess| keyword.

All the common Hessian update schemes are available.  For formulas, see
Schlegel [Schlegel:1987:AIMQC]_ and Bofill [Bofill:1994:1]_.

The Hessian may be computed during an optimization using the 
|optking__full_hess_every| keyword.

.. index:: 
   pair: geometry optimization; transition state
   pair: geometry optimization; IRC
   single: geometry optimization; constrained

Transition States and Reaction Paths
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Optking currently has two transition state algorithms. The current default is the
newer RS_I_RFO algorithm [Besalu:1998:265]_ . The old algorithm can be used by setting
`STEP_TYPE P_RFO` for `OPT_TYPE TS`

* Calculate a starting Hessian and optimize the "transition state" of
  linear water (note that without a reasonable starting geometry and
  Hessian, such a straightforward search often fails)::

   molecule h2o {
      O
      H 1 1.0
      H 1 1.0 2 160.0
   }

   set {
      basis dz
      full_hess_every 0
      opt_type ts
   }

   optimize('scf')

* At a transition state (planar HOOH), compute the second derivative, and
  then follow the intrinsic reaction path to the minimum::

   molecule hooh {
      symmetry c1
      H
      O 1 0.946347
      O 2 1.397780 1  107.243777
      H 3 0.946347 2  107.243777   1 0.0
   }

   set {
      basis dzp
      opt_type irc
      geom_maxiter 50
   }

   frequencies('scf')
   optimize('scf')

Constrained Optimizations
^^^^^^^^^^^^^^^^^^^^^^^^^
* Optimize a geometry (HOOH) at a frozen dihedral angle of 90 degrees. ::

   molecule {
     H
     O 1 0.90
     O 2 1.40 1 100.0
     H 3 0.90 2 100.0 1 90.0
   }

   set optking {
     frozen_dihedral = ("
       1 2 3 4
     ")
   }
   optimize('scf')

* To instead freeze the two O-H bond distances ::

   set optking {
     frozen_distance = ("
       1  2
       3  4
     ")
   }

For bends, the corresponding keyword is "frozen_bend".

* To freeze the cartesian coordinates of atom 2

.. code-block:: none

   freeze_list = """
     2 xyz
   """
   set optking frozen_cartesian $freeze_list

* To freeze only the y coordinates of atoms 2 and 3

.. code-block:: none

   freeze_list = """
     2 y
     3 y
   """
   set optking frozen_cartesian $freeze_list

* To optimize toward a value of 0.95 Angstroms for the distance between 
  atoms 1 and 3, as well as that between 2 and 4

.. code-block:: none

   set optking {
     ranged_distance = ("
       1  3 0.949 0.95
       2  4 0.949 0.95
     ")
   }

.. note:: 
    The effect of the frozen and ranged keywords is generally independent of
    how the geometry of the molecule was input (whether Z-matrix or Cartesian, etc.)..
    At this time; however, enforcing Cartesian constraints when using a zmatrix for
    molecular input is not supported. Freezing or constraining Cartesian coordinates
    requires Cartesian molecule input. If numerical errors results in symmetry 
    breaking, while Cartesian constraints are active, symmetrization cannot occur and
    an error will be raised, prompting you to restart the job.

* To scan the potential energy surface by optimizing at several fixed values
  of the dihedral angle of HOOH.

.. code-block:: none

   molecule hooh {
     0 1
     H  0.850718   0.772960    0.563468
     O  0.120432   0.684669   -0.035503
     O -0.120432  -0.684669   -0.035503
     H -0.850718  -0.772960    0.563468
   }
   
   set {
     basis cc-pvdz
     intrafrag_step_limit 0.1
   }

   lower_bound = [99.99, 109.99, 119.99, 129.99, 149.99]
   upper_bound = [100, 110, 120, 130, 140, 150]
   PES = []

   for lower, upper in zip(lower_bound, upper_bound):
   my_string = f"1 2 3 4 {lower} {upper}"
   set optking ranged_dihedral = $my_string
   E = optimize('scf')
   PES.append((upper, E))

   print("\n\tcc-pVDZ SCF energy as a function of phi\n")
   for point in PES:
     print("\t%5.1f%20.10f" % (point[0], point[1]))

* To scan the potential energy surface without the |optking__ranged_dihedral| keyword, a zmatrix
  can be used.

.. warning:: 
    Rotating dihedrals in large increments without allowing the molecule to relax
    in between increments can lead to unphysical geometries with overlapping functional groups in larger molecules,
    which may prevent successful constrained optimzations. Furthermore, such a relaxed scan of the PES does
    not always procude a result close to an IRC, or even a reaction path along which the energy changes in a
    continuous way.

.. code-block:: none

   molecule hooh {
     0 1 
     H   
     O 1 0.95
     O 2 1.39 1 103 
     H 3 0.95 2 103 1 D 

     D = 99

     units ang 
   }

   set {
     basis cc-pvdz
     intrafrag_step_limit 0.1 
     frozen_dihedral (" 1 2 3 4 ")
   }

   dihedrals = [100, 110, 120, 130, 140, 150]
   PES = []

   for phi in dihedrals:
     hooh.D = phi 
     E = optimize('scf')
     PES.append((phi, E))

   print("\n\tcc-pVDZ SCF energy as a function of phi\n")
   for point in PES:
     print("\t%5.1f%20.10f" % (point[0], point[1]))

Multi-Fragment Optimizations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _DimerSection:

In previous versions of optking, the metric for connecting atoms was increased until all atoms
were connected. This is the current behavior for |optking__frag_mode| `single`.
Setting |optking__frag_mode| to `multi` will now add a special
set of intermolecular coordinates between fragments - internally referred to as DimerFrag
coordinates (see `here <DimerIntro_>` for the brief description). 
For each pair of molecular fragments, a set of up to 3 reference points
are chosen on each fragment. Each reference point will be either an atom or a linear combination
of positions of atoms within that fragment. Stretches, bends, and dihedral angles between the two 
fragments will be created using these reference points. See 
:ref:`Dimer coordinate table <table:DimerFrag>` for how reference points are created.
For a set of three dimers A, B, and C, sets of coordinates are created between each pair:
AB, AC, and BC. Each `DimerFrag` may use different reference points. 
Creation of the intermolecular coordinates can be controlled through |optking__frag_ref_atoms| 
and |optking__interfrag_coords|. |optking__frag_ref_atoms| specifies which atoms 
(or linear combination of atoms) to use for the reference points and |optking__interfrag_coords|,
which encompasses |optking__frag_ref_atoms|, allows for constraints and labels to be added to the
intermolecular coordinates.

.. note:: Manual specification of the interfragment coordinates is supported for power users,
    and provides complete control of fragments' relative orientations.
    Setting |optking__interfrag_mode| to `multi` should suffice in almost all cases.
    :ref:`Dimer coordinate table <table:DimerFrag>`. provides the name and ordering
    convention for the coordinates.

* Basic multi-fragment optimization. Use automatically generated reference points.

.. code-block:: none

   memory 4GB 
   molecule mol {
       0 1 
       O   -0.5026452583       -0.9681078610       -0.4772692868
       H   -2.3292990446       -1.1611084524       -0.4772692868
       H   -0.8887241813        0.8340933116       -0.4772692868
       --  
       0 1 
       C    0.8853463281       -5.2235996493        0.5504918473
       C    1.8139169342       -2.1992967152        3.8040686146
       C    2.8624456357       -4.1143863257        0.5409035710
       C   -0.6240195463       -4.8153482424        2.1904642137
       C   -0.1646305764       -3.3031992532        3.8141619690
       C    3.3271056135       -2.6064153737        2.1669340785
       H    0.5244823836       -6.4459192939       -0.7478283184
       H    4.0823309159       -4.4449979205       -0.7680411190
       H   -2.2074914566       -5.7109913627        2.2110247636
       H   -1.3768100495       -2.9846751653        5.1327625515
       H    4.9209603634       -1.7288723155        2.1638694922
       H    2.1923374156       -0.9964630692        5.1155773223
       nocom
       units au
   }
   
   set {
       basis 6-31+G 
       frag_mode MULTI
   }
   
   optimize("mp2")

.. Warning:: The molecule input for psi4 has no effect upon optking, expect to provide Cartesian
    coordinates. Specifying independent fragments with the `--` seperator, will not trigger 
    optking to add specific interfragment coordinates. Use |optking__frag_mode| `multi`.

* Specify the reference points to use for coordinates via |optking__frag_ref_atoms|.
  Each list corresponds to a fragment. A list of indices denotes a linear combination
  of the atoms. In this case, the first reference point for the second fragment is the center
  of the benzene ring. Indexing starts at 1, so the second fragment in this example starts at index 4.

.. code-block:: none

   memory 4GB 
   molecule mol {
       0 1 
       O   -0.5026452583       -0.9681078610       -0.4772692868
       H   -2.3292990446       -1.1611084524       -0.4772692868
       H   -0.8887241813        0.8340933116       -0.4772692868
       --  
       0 1 
       C    0.8853463281       -5.2235996493        0.5504918473
       C    1.8139169342       -2.1992967152        3.8040686146
       C    2.8624456357       -4.1143863257        0.5409035710
       C   -0.6240195463       -4.8153482424        2.1904642137
       C   -0.1646305764       -3.3031992532        3.8141619690
       C    3.3271056135       -2.6064153737        2.1669340785
       H    0.5244823836       -6.4459192939       -0.7478283184
       H    4.0823309159       -4.4449979205       -0.7680411190
       H   -2.2074914566       -5.7109913627        2.2110247636
       H   -1.3768100495       -2.9846751653        5.1327625515
       H    4.9209603634       -1.7288723155        2.1638694922
       H    2.1923374156       -0.9964630692        5.1155773223
       nocom
       units au
   }
   
   set {
       basis 6-31+G 
       frag_mode MULTI

       # The line below specifies the reference points that will be used to construct the
       # interfragment coordinates between the two fragments (called A and B).
       # The format is the following:
       # [[A-1], [A-2], [A-3]], [[B-1], [B-2], [B-3]]
       #
       # In terms of atoms within each fragment, the line below chooses, for water:
       # H3 of water for the first reference point, O1 of water for the second reference point, and
       # H2 of water for the third reference point.
       # For benzene: the mean of the positions of all the C atoms, C2, one of the Carbon atoms,
       # and C6, another one of the carbon atoms.

       frag_ref_atoms [
           [[3], [1], [2]], [[4, 5, 6, 7, 8, 9], [5], [9]]
       ]   
   }
   
   optimize("mp2")

For even greater control, a dictionary can be passed to |optking__interfrag_coords|

The coordinates that are created between two dimers depend upon the number of atoms present
The fragments `A` and `B` have up to 3 reference atoms each as shown in
:ref:`Dimer coordinate table <table:DimerFrag>`.
The interfragment coordinates are named and can be frozen according to their names as show in 
example below. For specifying reference points, use 1 based indexing. 

.. _`table:DimerFrag`:

.. table:: Dimer coordinates

    +---------+----------+-------------+---------------------------------+
    | name    | type     | atom-labels | present, if                     |
    +=========+==========+=============+=================================+
    | RAB     | distance | A0-B0       | always                          | 
    +---------+----------+-------------+---------------------------------+
    | theta_A | angle    | A1-A0-B0    | A has > 1 atom                  |
    +---------+----------+-------------+---------------------------------+
    | theta_B | angle    | A0-B0-B1    | B has > 1 atom                  |
    +---------+----------+-------------+---------------------------------+
    | tau     | dihedral | A1-A0-B0-B1 | A and B have > 1 atom           |
    +---------+----------+-------------+---------------------------------+
    | phi_A   | dihedral | A2-A1-A0-B0 | A has > 2 atoms. Is not linear  |
    +---------+----------+-------------+---------------------------------+
    | phi_B   | dihedral | A0-B0-B1-B2 | B has > 2 atoms. Is not linear  |
    +---------+----------+-------------+---------------------------------+

* A constrained optimization is performed where the orientation of the two fragments is fixed but
  the distance between the fragments and all intrafragment coordinates are allowed to relax. In this
  example, the centers of the benzene and thiophene rings are selected for the first reference points.
  The methyl groups carbon and one hydrogen are selected for the other two reference points on the
  first fragments. For fragment two, two carbons of the benzene ring are chosen for the other reference points.

.. code-block:: none

   memory 4GB 
   molecule mol {
     C       -1.258686      0.546935      0.436840
     H       -0.683650      1.200389      1.102833
     C       -0.699036     -0.349093     -0.396608
     C       -2.693370      0.550414      0.355311
     H       -3.336987      1.206824      0.952052
     C       -3.159324     -0.343127     -0.536418
     H       -4.199699     -0.558111     -0.805894
     S       -1.883829     -1.212288     -1.301525
     C        0.786082     -0.656530     -0.606057
     H        1.387673     -0.016033      0.048976
     H        1.054892     -0.465272     -1.651226
     H        0.978834     -1.708370     -0.365860
     --
     C       -6.955593     -0.119764     -1.395442
     C       -6.977905     -0.135060      1.376787
     C       -7.111625      1.067403     -0.697024
     C       -6.810717     -1.314577     -0.707746
     C       -6.821873     -1.322226      0.678369
     C       -7.122781      1.059754      0.689090
     H       -7.226173      2.012097     -1.240759
     H       -6.687348     -2.253224     -1.259958
     H       -6.707325     -2.266920      1.222105
     H       -7.246150      1.998400      1.241304
     O       -6.944245     -0.111984     -2.805375
     H       -7.058224      0.807436     -3.049180
     C       -6.990227     -0.143507      2.907714
     H       -8.018305     -0.274985      3.264065
     H       -6.592753      0.807024      3.281508
     H       -6.368443     -0.968607      3.273516
     nocom
     unit angstrom
   }
   
   # Create a python dictionary and convert to string for pass through to optking
   MTdimer = """{
      "Natoms per frag": [12, 16],
      "A Frag": 1,
      "A Ref Atoms": [[1, 3, 4, 6, 8], [8], [11]],
      "A Label": "methylthiophene",
      "B Frag": 2,
      "B Ref Atoms": [[13, 14, 15, 16, 17, 18], [13], [15]],
      "B Label": "tyrosine",
      "Frozen": ["theta_A", "theta_B", "tau", "phi_A", "phi_B"],
   }"""
   
   set {
       basis 6-31+G 
       frag_mode MULTI
       interfrag_coords $MTdimer   
   }

   optimize("mp2")


Dealing with problematic optimizations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Although optking is continuously improved with robustness in mind, some
attempted optimizations will inevitably fail to converge to the desired minima.
For difficult cases, the following suggestions are made.

* As for any optimizer, computing the Hessian and limiting the step size will
  successfully converge a higher percentage of cases.  The default settings have
  been chosen because they perform efficiently for common, representative test sets.
  More restrictive, cautious steps are sometimes necessary.

* |optking__dynamic_level| allows optking to change the method of optimization
  toward algorithms that, while often less efficient, may help to converge difficult
  cases.  If this is initially set to 1, then optking, as poor steps are detected,
  will increase the dynamic level through several forms of more robust and cautious algorithms.
  The changes will reduce the trust radius, allow backward steps (partial line
  searching), add cartesian coordinates, switch to cartesian coordinates, and take
  steepest-descent steps.

* The developers have found the |optking__opt_coordinates| set to "BOTH" which
  includes both the redundant internal coordinate set, as well as cartesian coordinates,
  works well for systems with long 'arms' or floppy portions of a molecule poorly
  described by local internals.

* Optking does support the specification of ghost atoms. Certain internal coordinates such 
  as torsions become poorly defined when they contain near-linear bends. 
  An internal error `AlgError` may be raised in such cases. Optking will avoid such
  coordinates when choosing an initial coordinate system; however, they may arise in the course
  of an optimization. In such cases, try restarting from the most recent geometry.
  Alternatively, setting |optking__opt_coordinates| to cartesian will avoid any internal
  coordinate difficulties altogether. These coordinate changes can be automatically
  performed by turning |optking__dynamic_level| to 1.

.. warning:: In some cases, such as the coordinate issues described above, optking will reset to maintain
  a consistent history. If an error occurs in Psi4 due to |optking__geom_maxiter| being exceeded but
  the final step report indicates that optking has not taken |optking__geom_maxiter| steps, such a 
  reset has occured. Inspection will show that the step counter was reset to 1 somewhere in the
  optimization.

.. index:: 
   pair: geometry optimization; convergence criteria

Convergence Criteria
^^^^^^^^^^^^^^^^^^^^

Optking monitors five quantities to evaluate the progress of a geometry 
optimization. These are (with their keywords) the change in energy 
(|optking__max_energy_g_convergence|), the maximum element of 
the gradient (|optking__max_force_g_convergence|), the root-mean-square 
of the gradient (|optking__rms_force_g_convergence|), the maximum element
of displacement (|optking__max_disp_g_convergence|), and the 
root-mean-square of displacement (|optking__rms_disp_g_convergence|), 
all in internal coordinates and atomic units. Usually, these options will not 
be set directly. Primary control for geometry convergence lies with the keyword 
|optking__g_convergence| which sets the aforementioned in accordance 
with Table :ref:`Geometry Convergence <table:optkingconv>`.

|
|

.. _`table:optkingconv`:

.. table:: Summary of sets of geometry optimization criteria available in |PSIfour|

    +-----------------------------+----------------------------+----------------------------+----------------------------+----------------------------+----------------------------+
    | |optking__g_convergence|    | Max Energy                 | Max Force                  | RMS Force                  | Max Disp                   | RMS Disp                   |
    +=============================+============================+============================+============================+============================+============================+
    | NWCHEM_LOOSE [#fd]_         |                            | :math:`4.5 \times 10^{-3}` | :math:`3.0 \times 10^{-3}` | :math:`5.4 \times 10^{-3}` | :math:`3.6 \times 10^{-3}` |
    +-----------------------------+----------------------------+----------------------------+----------------------------+----------------------------+----------------------------+
    | GAU_LOOSE [#ff]_            |                            | :math:`2.5 \times 10^{-3}` | :math:`1.7 \times 10^{-3}` | :math:`1.0 \times 10^{-2}` | :math:`6.7 \times 10^{-3}` |
    +-----------------------------+----------------------------+----------------------------+----------------------------+----------------------------+----------------------------+
    | TURBOMOLE [#fd]_            | :math:`1.0 \times 10^{-6}` | :math:`1.0 \times 10^{-3}` | :math:`5.0 \times 10^{-4}` | :math:`1.0 \times 10^{-3}` | :math:`5.0 \times 10^{-4}` |
    +-----------------------------+----------------------------+----------------------------+----------------------------+----------------------------+----------------------------+
    | GAU [#fc]_ [#ff]_           |                            | :math:`4.5 \times 10^{-4}` | :math:`3.0 \times 10^{-4}` | :math:`1.8 \times 10^{-3}` | :math:`1.2 \times 10^{-3}` |
    +-----------------------------+----------------------------+----------------------------+----------------------------+----------------------------+----------------------------+
    | CFOUR [#fd]_                |                            |                            | :math:`1.0 \times 10^{-4}` |                            |                            |
    +-----------------------------+----------------------------+----------------------------+----------------------------+----------------------------+----------------------------+
    | QCHEM [#fa]_ [#fe]_         | :math:`1.0 \times 10^{-6}` | :math:`3.0 \times 10^{-4}` |                            | :math:`1.2 \times 10^{-3}` |                            |
    +-----------------------------+----------------------------+----------------------------+----------------------------+----------------------------+----------------------------+
    | MOLPRO [#fb]_ [#fe]_        | :math:`1.0 \times 10^{-6}` | :math:`3.0 \times 10^{-4}` |                            | :math:`3.0 \times 10^{-4}` |                            |
    +-----------------------------+----------------------------+----------------------------+----------------------------+----------------------------+----------------------------+
    | INTERFRAG_TIGHT [#fg]_      | :math:`1.0 \times 10^{-6}` | :math:`1.5 \times 10^{-5}` | :math:`1.0 \times 10^{-5}` | :math:`6.0 \times 10^{-4}` | :math:`4.0 \times 10^{-4}` |
    +-----------------------------+----------------------------+----------------------------+----------------------------+----------------------------+----------------------------+
    | GAU_TIGHT [#fc]_ [#ff]_     |                            | :math:`1.5 \times 10^{-5}` | :math:`1.0 \times 10^{-5}` | :math:`6.0 \times 10^{-5}` | :math:`4.0 \times 10^{-5}` |
    +-----------------------------+----------------------------+----------------------------+----------------------------+----------------------------+----------------------------+
    | GAU_VERYTIGHT [#ff]_        |                            | :math:`2.0 \times 10^{-6}` | :math:`1.0 \times 10^{-6}` | :math:`6.0 \times 10^{-6}` | :math:`4.0 \times 10^{-6}` | 
    +-----------------------------+----------------------------+----------------------------+----------------------------+----------------------------+----------------------------+

.. rubric:: Footnotes

.. [#fa] Default
.. [#fb] Baker convergence criteria are the same.
.. [#fc] Counterpart NWCHEM convergence criteria are the same.
.. [#fd] Convergence achieved when all active criteria are fulfilled.
.. [#fe] Convergence achieved when **Max Force** and one of **Max Energy** or **Max Disp** are fulfilled.
.. [#ff] Normal convergence achieved when all four criteria (**Max Force**, **RMS Force**,
         **Max Disp**, and **RMS Disp**) are fulfilled. To help with flat 
         potential surfaces, alternate convergence achieved when 100\ :math:`\times`\ *rms force* is less 
         than **RMS Force** criterion.
.. [#fg] Compensates for difficulties in converging geometry optmizations of supermolecular complexes 
         tightly, where large *rms disp* and *max disp* may result from flat potential surfaces even when
         *max force* and/or *rms force* are small.

For ultimate control, specifying a value for any of the five monitored options activates that
criterium and overwrites/appends it to the criteria set by |optking__g_convergence|.
Note that this revokes the special convergence arrangements detailed in notes [#fe]_ and [#ff]_ 
and instead requires all active criteria to be fulfilled to 
achieve convergence. To avoid this revokation, turn on keyword |optking__flexible_g_convergence|.

.. index::
   pair: geometry optimization; output

Interface to GeomeTRIC
^^^^^^^^^^^^^^^^^^^^^^

The GeomeTRIC optimizer developed by Wang and Song [Wang:2016:214108]_ may be used in place of
Psi4's native Optking optimizer. GeomeTRIC uses a translation-rotation-internal coordinate (TRIC)
system that works well for optimizing geometries of systems containing noncovalent interactions.

Use of the GeomeTRIC optimizer is specified with the ``engine`` argument to
:py:func:`~psi4.driver.optimize`. The optimization will respect the keywords |optking__g_convergence|
and |optking__geom_maxiter|. Any other GeomeTRIC-specific options (including constraints)
may be specified with the ``optimizer_keywords`` argument to :py:func:`~psi4.driver.optimize`.
Constraints may be placed on cartesian coordinates, bonds, angles, and dihedrals, and they can be
used to either freeze a coordinate or set it to a specific value. See the `GeomeTRIC github
<https://github.com/leeping/geomeTRIC>`_ 
for more information on keywords and JSON specification of constraints.

* Optimize the water molecule using GeomeTRIC::

   molecule h2o {
      O
      H 1 1.0
      H 1 1.0 2 160.0
   }

   set {
      maxiter 100
      g_convergence gau
   }

   optimize('hf/cc-pvdz', engine='geometric')

* Optimize the water molecule using GeomeTRIC, with one of the two OH bonds constrained to 2.0 au
  and the HOH angle constrained to 104.5 degrees::

   molecule h2o {
      O
      H 1 1.0
      H 1 1.0 2 160.0
   }

   set {
      maxiter 100
      g_convergence gau
   }

   geometric_keywords = { 
     'coordsys' : 'tric',
     'constraints' : { 
     'set' : [{'type'    : 'distance',
               'indices' : [0, 1], 
               'value'   : 2.0 },
              {'type'    : 'angle',
               'indices' : [1, 0, 2], 
               'value'   : 104.5 }]
      }   
   }   

   optimize('hf/cc-pvdz', engine='geometric', optimizer_keywords=geometric_keywords)

* Optimize the benzene/water dimer using GeomeTRIC, with the 6 carbon atoms of benzene frozen in 
  place::

   molecule h2o {
     C            0.833     1.221    -0.504
     H            1.482     2.086    -0.518
     C            1.379    -0.055    -0.486
     H            2.453    -0.184    -0.483
     C            0.546    -1.167    -0.474
     H            0.971    -2.162    -0.466
     C           -0.833    -1.001    -0.475
     H           -1.482    -1.867    -0.468
     C           -1.379     0.275    -0.490
     H           -2.453     0.404    -0.491
     C           -0.546     1.386    -0.506
     H           -0.971     2.381    -0.524
     --
     O            0.000     0.147     3.265
     H            0.000    -0.505     2.581
     H            0.000     0.965     2.790
     no_com
     no_reorient
   }

   set {
      maxiter 100
      g_convergence gau
   }

   geometric_keywords = { 
     'coordsys' : 'tric',
     'constraints' : { 
     'freeze' : [{'type'    : 'xyz',
                  'indices' : [0, 2, 4, 6, 8, 10]}]
      }   
   }   

   optimize('hf/cc-pvdz', engine='geometric', optimizer_keywords=geometric_keywords)


Output
^^^^^^

The progress of a geometry optimization can be monitored by grepping the output file for the
tilde character (``~``). This produces a table like the one below that shows
for each iteration the value for each of the five quantities and whether the criterion
is active and fulfilled (``*``), active and unfulfilled ( ),  or inactive (``o``). ::

   --------------------------------------------------------------------------------------------- ~
    Step     Total Energy     Delta E     MAX Force     RMS Force      MAX Disp      RMS Disp    ~
   --------------------------------------------------------------------------------------------- ~
     Convergence Criteria    1.00e-06 *    3.00e-04 *             o    1.20e-03 *             o  ~
   --------------------------------------------------------------------------------------------- ~
       1     -38.91591820   -3.89e+01      6.91e-02      5.72e-02 o    1.42e-01      1.19e-01 o  ~
       2     -38.92529543   -9.38e-03      6.21e-03      3.91e-03 o    2.00e-02      1.18e-02 o  ~
       3     -38.92540669   -1.11e-04      4.04e-03      2.46e-03 o    3.63e-02      2.12e-02 o  ~
       4     -38.92548668   -8.00e-05      2.30e-04 *    1.92e-04 o    1.99e-03      1.17e-03 o  ~
       5     -38.92548698   -2.98e-07 *    3.95e-05 *    3.35e-05 o    1.37e-04 *    1.05e-04 o  ~

The full list of keywords for optking is provided in Appendix :ref:`apdx:optking`.

Information on the Psithon function that drives geometry optimizations is provided
at :py:func:`~psi4.driver.optimize`.

Important User Changes from cpp-optking
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* `FIXED_COORD` keywords have been generalized to `RANGED_COORD` e.g. |optking__ranged_distance|

* Detailed optimization is now printed through the python logging system. If more information about
  the optimization is needed. Please see `\<output_name\>.log`

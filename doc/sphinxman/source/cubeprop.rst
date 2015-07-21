
.. include:: autodoc_abbr_options_c.rst

.. index::
   Cube
   cubeprop
   visualization

.. _`sec:cubeprop`:

Generation of cube files
========================

.. codeauthor:: Robert M. Parrish and Francesco A. Evangelista
.. sectionauthor:: Francesco A. Evangelista

Introduction
------------

|PSIfour| has the ability to export cube files that store information about
basis functions, molecular orbitals, the electron density, and
the electrostatic potential (ESP).  Cube files store the value of a scalar
quantity on a regular Cartesian grid, and can be visualized with several
visualization programs, some of which are free, like VMD
(http://www.ks.uiuc.edu/Research/vmd/).


An example utilization of the code is::

   molecule h2o {
   0 1
   O
   H 1 1.0
   H 1 1.0 2 104.5
   }

   set basis cc-pvdz
   set scf_type df
   set freeze_core True
   set cubeprop_tasks ['orbitals']
   set cubeprop_orbitals [5,6,-5,-6]

   energy('scf')
   cubeprop()

In this example, the ``cubeprop()`` call after the ``energy('scf')`` command
executes the cubeprop code.  The array ``CUBEPROP_TASKS`` specifies which
tasks should be executed.  In this case the task ``'orbitals'`` generates cube
files for orbitals.  The ``CUBEPROP_ORBITALS`` option specifies that cube files
should be generated only for alpha orbitals 5 (HOMO) and 6 (LUMO) and
beta orbitals 5 (indicated as -5) and 6.
If the option ``CUBEPROP_ORBITALS`` is not provided, then cube files are
generated for all orbitals.
After running, the above input will generate four files: `Psi_a_5.cube`,
`Psi_a_6.cube`, `Psi_b_5.cube`, and `Psi_b_6.cube`.

.. tip:: If your cube plots are too coarse, try to decrease the grid spacing via
    the option ``CUBIC_GRID_SPACING``.  If the edges of your plot are cut then
    increase the size of the grid via the option ``CUBIC_GRID_OVERAGE``.

Cubeprop Tasks
--------------

The cubeprop utility can be probided a list of tasks to perform.
Tasks are specified by the ``CUBEPROP_TASKS`` option, which is a list of strings
that identify the tasks.  Several tasks are available. These include:

ORBITALS [Default if  ``CUBEPROP_TASKS`` is not specified]
    Produces cube representations of the molecular orbitals
    :math:`\psi_q(\mathbf{r})`.  Orbitals are sorted according to increasing
    orbital energy ignoring symmetry.
DENSITY
    This task can be used to obtain the alpha and beta electron densities,
    :math:`\rho_\alpha(\mathbf{r})` and :math:`\rho_\beta(\mathbf{r})`, together
    with the total density
    :math:`\rho(\mathbf{r}) = \rho_\alpha(\mathbf{r}) + \rho_\beta(\mathbf{r})`,
    and the spin density
    :math:`\rho(\mathbf{r}) = \rho_\alpha(\mathbf{r}) - \rho_\beta(\mathbf{r})`.
BASIS_FUNCTIONS
    This task is useful to produce cube files of the atomic orbital basis
    functions :math:`\chi_\mu(\mathbf{r})`.
ESP
    Calculates the total (nuclear + electronic) electrostatic potential
    :math:`V(\mathbf{r})`.

.. note:: The ``ESP`` task requires the user to specify a density-fitting basis
    via the |scf__df_basis_scf| keyword.

.. warning:: It is important to specify the ``CUBEPROP_ORBITALS`` option when
   dealing with large molecules to avoid running out of disk space.
   For example, using the default grid spacing of
   0.2 |AA|\ ngstroms, the size of a single cube file for a molecule like water
   is of the order of 1.4 MB.  For a molecule with 200 basis functions, the cube
   files for all the orbitals occupy more than half a GB.

Cubeprop Options
----------------

CUBEPROP_FILEPATH
    The directory in which to place the cube files.  The default is the
    directory that contains the input file.

CUBEPROP_ORBITALS
    A list of the orbitals for which cube files are generated.  Alpha orbitals
    are specified by positive integers, where the lowest orbital is 1.
    Beta orbitals are indicated by prefixing the number with a minus sign.
    Orbitals are ordered according to their energy.
    Notice that if the ``CUBEPROP_ORBITALS`` list is empty then
    cubeprop will produce cube files for all alpha and beta orbitals.

CUBEPROP_BASIS_FUNCTIONS
    A list of basis functions for which cube files are generated.  Basis
    functions are numbered starting from 1.

CUBIC_GRID_SPACING
    A vector that specifies the grid spacing in the X, Y, and Z directions.
    The default value is [0.2,0.2,0.2] |AA|\ ngstroms.

CUBIC_GRID_OVERAGE
    A vector that controls the spatial extent of the grid in the X, Y, and Z
    directions.
    The default value is [4.0,4.0,4.0] |AA|\ ngstroms.



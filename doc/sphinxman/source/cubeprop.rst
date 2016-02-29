
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

Orbital Visualization with VMD
==============================

Included in |PSIfour| is functionality to automatically render specified surfaces, including molecular orbitals,
densities, and basis functions, using VMD. The vmd_cube.py script takes the .cube files generated
in a calculation and generates images alinged with user-input specifications. The script is located
in :source:`/share/scripts/vmd_cube.py`.

Script Prerequisites
--------------------

1. VMD must be installed, and it can be downloaded for free at (http://www.ks.uiuc.edu/Research/vmd/). Additionally,
   the script needs to know where to find the VMD executable, and this is defined as VMDPATH. VMDPATH must be defined as
   an environment variable.

2. To generate images with multiple surfaces, ImageMagick must also be installed. ImageMagick is a free program which
   can be installed using homebrew/pip or from http://www.imagemagick.org/script/binary-releases.php .

3. With ImageMagick installed, an environment variable called montage needs to be created which points to the montage executable.
   This executable can be found in the /bin/ sub-directory wherever ImageMagick was installed.


Running the Script
------------------

1. Run a |PSIfour| calculation, generating .cube files as detailed in the above documentation.

2. Copy vmd_script.py into a directory where the image files are desired, and pass the directory
   pointing to the .cube files as an argument to run::
	
	python vmd_cube.py /path/to/cube/files/

   Alternatively, the script can be run in the same directory as the cube files with no need to pass the
   directory as an argument.

3. For an additional image containing all surfaces in an array (very useful for hand-picking orbital spaces), set the montage
   flag to True::

	python vmd_cube.py /path/to/cube/files/ --montage=True

3. As an example, take the cube files generated from the water calculation from the above input file. Using the script and montage,
   the alpha molecular orbitals, for example, can be rendered and output in one image with::

	python vmd_cube.py /path/to/cube/files/ --montage=True --opacity=0.5 --rx=90 --ry=60

The desired image in this case is called "AlphaMOs.tga", and looks like	this:

.. image:: /AlphaMOs.png
    :align: center
    :scale: 100%
    :alt: Alpha MOs	


Script Options
--------------

    >>> ./vmd_cube.py --help
    usage: vmd_cube.py [-h] [--color1 [<integer>]] [--color2 [<integer>]]
                       [--iso [<isovalue>]] [--rx [<angle>]] [--ry [<angle>]]
                       [--rz [<angle>]] [--tx [<angle>]] [--ty [<angle>]]
                       [--tz [<angle>]] [--opacity [<opacity>]]
                       [--scale [<factor>]] [--montage [MONTAGE]]
                       [--imagesize [<integer>]] [--fontsize [<integer>]]
                       [<cubefile dir>]
                                                                    .
    vmd_cube is a script to render cube files with vmd. To generate cube files
    with Psi4 add the command cubeprop() at the end of your input file.
                                                                    .
    positional arguments:
      <cubefile dir>        The directory containing the cube files.
                                                                    .
    optional arguments:
      -h, --help            show this help message and exit
      --color1 [<integer>]  the color ID of surface 1 (integer, default = 3)
      --color2 [<integer>]  the color ID of surface 2 (integer, default = 23)
      --iso [<isovalue>]    the isosurface value (float, default = 0.05)
      --rx [<angle>]        the x-axis rotation angle (float, default = 30.0)
      --ry [<angle>]        the y-axis rotation angle (float, default = 40.0)
      --rz [<angle>]        the z-axis rotation angle (float, default = 15.0)
      --tx [<angle>]        the x-axis translation (float, default = 0.0)
      --ty [<angle>]        the y-axis translation (float, default = 0.0)
      --tz [<angle>]        the z-axis translation (float, default = 0.0)
      --opacity [<opacity>]
                            opacity of the isosurface (float, default = 1.0)
      --scale [<factor>]    the scaling factor (float, default = 1.0)
      --montage [MONTAGE]   call montage to combine images. (string, default =
                            false)
      --imagesize [<integer>]
                            the size of each image (integer, default = 250)
      --fontsize [<integer>]
                            the font size (integer, default = 20)


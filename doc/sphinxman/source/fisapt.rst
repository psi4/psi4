.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2017 The Psi4 Developers.
.. #
.. # The copyrights for code used from other parties are included in
.. # the corresponding files.
.. #
.. # This program is free software; you can redistribute it and/or modify
.. # it under the terms of the GNU General Public License as published by
.. # the Free Software Foundation; either version 2 of the License, or
.. # (at your option) any later version.
.. #
.. # This program is distributed in the hope that it will be useful,
.. # but WITHOUT ANY WARRANTY; without even the implied warranty of
.. # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
.. # GNU General Public License for more details.
.. #
.. # You should have received a copy of the GNU General Public License along
.. # with this program; if not, write to the Free Software Foundation, Inc.,
.. # 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
.. #
.. # @END LICENSE
.. #

.. include:: autodoc_abbr_options_c.rst

.. index::
   single: FISAPT
   pair: FISAPT; theory

.. _`sec:fisapt`:

F/I-SAPT: Functional Group and/or Intramolecular SAPT
=====================================================

.. codeauthor:: Robert M. Parrish
.. sectionauthor:: Robert M. Parrish

*Module:* :ref:`Keywords <apdx:fisapt>`, :ref:`PSI Variables
<apdx:fisapt_psivar>`, :source:`FISAPT <psi4/src/psi4/fisapt>`

The FISAPT module provides two extensions to standard SAPT theory to allow for
(1) an effective two-body partition of the various SAPT terms to localized
chemical functional groups (F-SAPT) and (2) a means to compute the SAPT
interaction between two moieties within the embedding field of a third body
(I-SAPT). F-SAPT is designed to provide additional insight into the chemical
origins of a noncovalent interaction, while I-SAPT allows for one to perform
a SAPT analysis for intramolecular interactions. F-SAPT and I-SAPT can be
deployed together in this module, yielding "F/I-SAPT." All F/I-SAPT computations
in |PSIfour| use density-fitted SAPT0 as the underlying SAPT methodology. Interested
users should consult the manual page for Ed Hohenstein's :ref:`SAPT0 <sec:sapt>` code
and the SAPT literature to understand the specifics of SAPT0 before beginning
with F/I-SAPT0.

F-SAPT is detailed over two papers: [Parrish:2014:044115]_ on our much-earlier
"atomic" SAPT (A-SAPT) and [Parrish:2014:4417]_ on the finished "functional
group" SAPT (F-SAPT). An additional paper describes how to use F-SAPT to analyze
differences under functional group substitutions [Parrish:2014:17386]_.  I-SAPT
is explained in [Parrish:2015:051103]_. There is also a reasonably-detailed
review of the aims of A/F/I-SAPT and the existing state-of-the-art in the field
in the introduction chapter on partitioned SAPT methods in `Parrish's thesis
<https://smartech.gatech.edu/handle/1853/53850>`_.

A video tutorial series for the use of the FISAPT module is available `here
<https://www.youtube.com/playlist?list=PLg_zUQpVYlA1Tc1X_HgAbqnFcHNydqN7W>`_.
Specific videos in the series include:

- `F-SAPT#1
  <https://www.youtube.com/watch?v=J22J0wh4mVo&index=1&list=PLg_zUQpVYlA1Tc1X_HgAbqnFcHNydqN7W>`_.
  Describes the use of F-SAPT to analyze the
  distribution of the intermolecular interaction energy components between the
  various hydroxyl and phenyl moieties of the phenol dimer.
- `F-SAPT#2
  <https://www.youtube.com/watch?v=fqlzXsayec0&index=2&list=PLg_zUQpVYlA1Tc1X_HgAbqnFcHNydqN7W>`_.
  Discusses how to plot the order-1 F-SAPT analysis with PyMol and perform a
  "difference F-SAPT" analysis
- `I-SAPT#1
  <https://www.youtube.com/watch?v=fD6mu_tTG_c&index=3&list=PLg_zUQpVYlA1Tc1X_HgAbqnFcHNydqN7W>`_.
  Describes the use of I-SAPT to analyze the interaction between the two phenol
  groups in a 2,4-pentanediol molecule.
- `I-SAPT#2
  <https://www.youtube.com/watch?v=hDbonAOD5dY&index=4&list=PLg_zUQpVYlA1Tc1X_HgAbqnFcHNydqN7W>`_.
  Discusses how to plot the density fields and ESPs of the various moieties of
  the I-SAPT embedding scheme with PyMol
- `F/I-SAPT Options
  <https://www.youtube.com/watch?v=KFkPKSUZVfI&index=5&list=PLg_zUQpVYlA1Tc1X_HgAbqnFcHNydqN7W>`_.
  Details all of the more-advanced options in the F/I-SAPT code (rarely needed).

The scripts discussed below are located in :source:`psi4/share/psi4/fsapt`.

F-SAPT: A Representative Example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Below, we show an example of using F-SAPT0/jun-cc-pVDZ to analyze the
distribution of the intermolecular interaction energy components between the
various hydroxyl and phenyl moieties of the phenol dimer. This example is
explicitly included in :srcsample:`fsapt1`. A video
lecture explaining this example is available `F-SAPT#1
<https://www.youtube.com/watch?v=J22J0wh4mVo&index=1&list=PLg_zUQpVYlA1Tc1X_HgAbqnFcHNydqN7W>`_,
while an additional video describing how to plot the order-1 F-SAPT analysis
with PyMol and perform a "difference F-SAPT"
analysis is available `F-SAPT#2
<https://www.youtube.com/watch?v=fqlzXsayec0&index=2&list=PLg_zUQpVYlA1Tc1X_HgAbqnFcHNydqN7W>`_::

    memory 1 GB
    
    molecule mol {
    0 1
    O    -1.3885044    1.9298523   -0.4431206
    H    -0.5238121    1.9646519   -0.0064609
    C    -2.0071056    0.7638459   -0.1083509
    C    -1.4630807   -0.1519120    0.7949930
    C    -2.1475789   -1.3295094    1.0883677
    C    -3.3743208   -1.6031427    0.4895864
    C    -3.9143727   -0.6838545   -0.4091028
    C    -3.2370496    0.4929609   -0.7096126
    H    -0.5106510    0.0566569    1.2642563
    H    -1.7151135   -2.0321452    1.7878417
    H    -3.9024664   -2.5173865    0.7197947
    H    -4.8670730   -0.8822939   -0.8811319
    H    -3.6431662    1.2134345   -1.4057590
    --
    0 1
    O     1.3531168    1.9382724    0.4723133
    H     1.7842846    2.3487495    1.2297110
    C     2.0369747    0.7865043    0.1495491
    C     1.5904026    0.0696860   -0.9574153
    C     2.2417367   -1.1069765   -1.3128110
    C     3.3315674   -1.5665603   -0.5748636
    C     3.7696838   -0.8396901    0.5286439
    C     3.1224836    0.3383498    0.8960491
    H     0.7445512    0.4367983   -1.5218583
    H     1.8921463   -1.6649726   -2.1701843
    H     3.8330227   -2.4811537   -0.8566666
    H     4.6137632   -1.1850101    1.1092635
    H     3.4598854    0.9030376    1.7569489
    symmetry c1
    no_reorient
    no_com
    }
    
    set {
    basis         jun-cc-pvdz
    scf_type df
    guess sad
    freeze_core true
    }
    
    energy('fisapt0')

This file runs a DF-HF computation on the full dimer using |PSIfours| existing
SCF code. The monomer SCF computations are performed inside the FISAPT module,
following which a complete DF-SAPT0 computation is performed. Additional bits of
analysis are performed to generate the order-2 partition of the SAPT terms to
the level of nuclei and localized occupied orbitals |w--w| this generally does not
incur much additional overhead beyond a standard SAPT0 computations. The
nuclear/orbital partition data is written to the folder :file:`fsapt/` in the same
directory as the input file (this can be changed by |fisapt__fisapt_fsapt_filepath|).

One obtains the desired F-SAPT partition by post-processing the data in
:file:`fsapt/`. Within this dir, the user is expected to provide the ASCII files
:file:`fA.dat` and :file:`fB.dat`, which describe the assignment of atoms to chemical
functional groups using 1-based ordering. *E.g.*, for the problem at hand,
:file:`fA.dat` contains::

    OH 1 2
    PH 3 4 5 6 7 8 9 10 11 12 13

while :file:`fB.dat` contains::
    
    OH 14 15
    PH 16 17 18 19 20 21 22 23 24 25 26

At this point, the user should run the ``fsapt.py`` post-processing script in
the ``fsapt`` directory as::

    >>> fsapt.py

This will generate, among other files, the desired functional-group partition in
``fsapt.dat``. For our problem, the bottom of this file contains the finished
partition::

    Frag1     Frag2          Elst      Exch     IndAB     IndBA      Disp     Total
    OH        OH           -8.425     6.216    -0.583    -1.512    -1.249    -5.553
    OH        PH            1.392     0.716     0.222    -0.348    -0.792     1.189
    PH        OH           -2.742     0.749    -0.147    -0.227    -0.674    -3.040
    PH        PH            0.680     2.187     0.007    -0.208    -2.400     0.266
    OH        All          -7.033     6.931    -0.362    -1.860    -2.040    -4.364
    PH        All          -2.062     2.936    -0.140    -0.435    -3.074    -2.774
    All       OH          -11.167     6.965    -0.730    -1.739    -1.923    -8.594
    All       PH            2.072     2.903     0.229    -0.556    -3.191     1.456
    All       All          -9.095     9.867    -0.501    -2.295    -5.114    -7.138

Note that the assignment of linking sigma bond contributions is a small point of
ambiguity in F-SAPT. The ``fsapt.dat`` file presents the "links-by-charge"
assignment at the top and the "links by 50-50" assignment at the bottom. We
generally prefer the latter, but both generally give qualitatively identical
energetic partitions.

Users should check the files ``fragA.dat`` and ``fragB.dat`` to ensure that
there is not too much charge delocalization from one fragment to another. This
is presented in the "Orbital Check" section in these files |w--w| a value larger than
0.1 docc is an indication that the picture of localizable functional groups may
be breaking down. We also *strongly discourage* the cutting of double,
triple, or aromatic bonding motifs when partitioning the molecule into fragments
|w--w| cuts across only simple sigma bonds are encouraged.

Order-1 Visualization with PyMol
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``fsapt.py`` script above also generates a number of order-1 ``.pdb`` files
that can be used to get a quick qualitative picture of the F-SAPT partition. The
preferred way to do this is to use PyMol to make plots of the molecular geometry
with the atoms colored according to their order-1 F-SAPT contributions. We have
a set of template ``.pymol`` scripts to help with this process. These can be
obtained by running::

    >>> copy_pymol.py

and then in PyMol::

    >>> @run.pymol

This last command runs all of the individual ``.pymol`` files (*e.g.*,
``Elst.pymol``), which in turn load in the molecule and order-1 analysis
(contained in the ``.pdb`` file), set up the visualization, and render a
``.png`` image of the scene. Generally the view orientation and some specific
details of the ``.pymol`` files require some small tweaks to permit
publication-quality renderings.

.. image:: /Total.png
    :align: center
    :scale: 30%
    :alt: Total Order-1 F-SAPT0

Difference F-SAPT Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^

For those interested in taking the differences between two F-SAPT partitions
(*e.g.*, to see how a substituent modulates a noncovalent interaction), we have
the ``fsapt-diff.py`` script to help with this. This is invoked as::

    >>> fsapt-diff.py source-fsapt-dir1 source-fsapt-dir2 target-diff-fsapt-dir

Where the use has already performed ``fsapt.py`` analysis using the same
functional group names in ``source-fsapt-dir-1`` and ``source-fsapt-dir-2``. The
difference F-SAPT partition entries are computed as :math:`E^{\Delta} = E^{1} -
E^{2}`, and the geometries for order-1 ``.pdb`` visualization files are taken
from system 1.

I-SAPT: A Representative Example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Below, we show an example of using I-SAPT0/jun-cc-pVDZ to analyze the
interaction between the two phenol groups in a 2,4-pentanediol molecule.
This example is
explicitly included in :srcsample:`isapt1`. A video
lecture explaining this example is available
`I-SAPT#1 <https://www.youtube.com/watch?v=fD6mu_tTG_c&index=3&list=PLg_zUQpVYlA1Tc1X_HgAbqnFcHNydqN7W>`_,
while an additional video describing how to plot the density and ESP fields from
the I-SAPT embedding procedure is available `I-SAPT#2 <https://www.youtube.com/watch?v=hDbonAOD5dY&index=4&list=PLg_zUQpVYlA1Tc1X_HgAbqnFcHNydqN7W>`_::

    memory 1 GB
    
    molecule mol {
    0 1
    O          0.39987        2.94222       -0.26535
    H          0.05893        2.05436       -0.50962
    --
    0 1
    O          0.48122        0.30277       -0.77763
    H          0.26106       -0.50005       -1.28451
    --
    0 1
    C          2.33048       -1.00269        0.03771
    C          1.89725        0.31533       -0.59009
    C          2.28232        1.50669        0.29709
    C          1.82204        2.84608       -0.29432
    C          2.37905        4.02099        0.49639
    H          3.41246       -1.03030        0.19825
    H          2.05362       -1.84372       -0.60709
    H          1.82714       -1.16382        0.99734
    H          2.36243        0.42333       -1.57636
    H          3.36962        1.51414        0.43813
    H          1.81251        1.38060        1.28140
    H          2.14344        2.92967       -1.33843
    H          3.47320        4.02400        0.48819
    H          2.03535        3.99216        1.53635
    H          2.02481        4.96785        0.07455
    symmetry c1
    no_reorient
    no_com
    }
    
    # => Standard Options <= #
    
    set {
    basis jun-cc-pvdz
    scf_type df
    guess sad
    freeze_core true
    fisapt_do_plot true  # For extra analysis
    }
    
    energy('fisapt0')

This is essentially the same input as for F-SAPT, except that the molecular
system is now divided into *three* moieties |w--w| subsystems A and B whose
intramolecular interaction we wish to compute, and a linking unit C.  This file
runs a DF-HF computation on the full system using |PSIfours| existing SCF code.
At the start of the FISAPT code, the occupied orbitals are localized and divided
by charge considerations into A, B, C, and link sets. By default, linking sigma
bonds are assigned to C (this can be changed by the |fisapt__fisapt_link_assignment|
options). Then, non-interacting Hartree--Fock solutions for A and B are optimized
in the embedding field of the linking moiety C. At this point, A and B are not
interacting with each other, but have any potential covalent links or other
interactions with C built in by the embedding.  A standard F-SAPT0 computation
is then performed between A and B, yielding the I-SAPT interaction energy. Any
F-SAPT considerations are also possible when I-SAPT is performed |w--w| F and I are
completely direct-product-separable considerations. 

Cube File Visualization with PyMol
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Setting |fisapt__fisapt_do_plot| ``true`` above generates a set of ``.cube`` files
containing the densities and ESPs of the various subsystems in the I-SAPT
embedding procedure. These can be used to gain a detailed understanding of the
intermolecular partition and the polarization between non-interacting and
Hartree--Fock-interacting moieties. We have developed a set of template
``.pymol`` scripts to help with this process. These can be obtained by running::

    >>> copy_pymol2.py

and then in PyMol::

    >>> @run.pymol

This last command runs all of the individual ``.pymol`` files (*e.g.*,
``DA.pymol``), which in turn load in the molecule and cube file data
(contained in the ``.cube`` file), set up the visualization, and render a
``.png`` image of the scene. Generally the view orientation and some specific
details of the ``.pymol`` files require some small tweaks to permit
publication-quality renderings.

.. image:: /VA.png
    :align: center
    :scale: 50%
    :alt: ESP of monomer A

F/I-SAPT Keywords
^^^^^^^^^^^^^^^^^

The input files described above cover roughly 90% of all F/I-SAPT analyses. For
more delicate or involved problems, there are a large number of user options
that permit the customization of the I-SAPT subsystem partition, the convergence
of the IBO localization procedure, numerical thresholds, etc. We have an entire
video tutorial devoted to these options `F/I-SAPT Options <https://www.youtube.com/watch?v=KFkPKSUZVfI&index=5&list=PLg_zUQpVYlA1Tc1X_HgAbqnFcHNydqN7W>`_.
Direct source-code documentation on these options is available :ref:`here
<apdx:fisapt_psivar>`_.

Additional Notes
^^^^^^^^^^^^^^^^

.. caution:: In constrast to Ed Hohenstein's SAPT0 code, FISAPT uses the -JKFIT
  auxiliary basis sets for all Fock-type terms (*e.g.*, electrostatics, exchange,
  induction, and core Fock matrix elements in exchange-dispersion), and the -RI
  auxiliary basis sets *only* for the dispersion term. Ed's code uses the -RI
  basis sets for all SAPT terms, which can be problematic for heavy elements.
  As such, Ed's SAPT0 code will yield slightly different results than FISAPT. The
  differences should be very minor for up to and including second-row elements,
  after which point one needs to use the |sapt__df_basis_elst| option in Ed's code to
  provide an accurate result. 


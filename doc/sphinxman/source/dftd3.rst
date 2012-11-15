
.. include:: autodoc_abbr_options_c.rst

.. index:: DFTD3
.. _`sec:dftd3`:

Interface to DFTD3 by S. Grimme
===============================

.. codeauthor:: Lori A. Burns
.. sectionauthor:: Lori A. Burns

.. *Module:* :ref:`Keywords <apdx:mrcc>`, :ref:`PSI Variables <apdx:mrcc_psivar>`, :source:`MRCC <src/bin/mrcc>`

|PSIfour| contains code to interface to the DFTD3 program of S. Grimme, which is freely
downloadable from `http://toc.uni-muenster.de/DFTD3/ <http://toc.uni-muenster.de/DFTD3/>`_).

Installation
~~~~~~~~~~~~

Follow the instructions provided with the source to build the DFTD3
program (essentially, unpack the source, edit the Makefile to select a
Fortran compiler, and run make). For the moment, you must apply a patch,
:source:`lib/scripts/patch_grimme_dftd3`, to the dftd3 source before building.
To be used by |PSIfour|, the program binary (``dftd3``) must be found in
your :envvar:`PATH`. If |PSIfour| is unable to execute the binary, an
error will be reported.

    >>> cd dftd3
    >>> ls 
    dftd3.tar
    patch_grimme_dftd3
    >>> tar -xvf dftd3.tar
    copyc6.f
    dftd3.f
    Makefile
    man.pdf
    pars.f
    param
    >>> patch < patch_grimme_dftd3 
    patching file Makefile
    patching file dftd3.f
    >>> make
    making dftd3.o from dftd3.f
    gfortran   -O  -c dftd3.f -o dftd3.o
    making copyc6.o from copyc6.f
    gfortran   -O  -c copyc6.f -o copyc6.o
    gfortran dftd3.o copyc6.o    -o ./dftd3 
    >>> ls
    Makefile           copyc6.o           dftd3.f            dftd3.tar          param              patch_grimme_dftd3
    copyc6.f           dftd3              dftd3.o            man.pdf            pars.f

Theory
~~~~~~

The local or semilocal character of conventional density functionals
necessarily leads to neglect of the long-range correlation interactions
which capture attractive van der Waals forces. Initially proposed by Yang
[Wu:2002:515]_ and assiduously developed by Grimme, [Grimme:2004:1463]_
[Grimme:2006:1787]_ [Grimme:2010:154104]_ the DFT+Dispersion method
appends to the base functional a scaled, damped, and fitted leading term
to the well-known dispersion energy series, :math:`E_{disp} = -C_6/R^6
-C_8/R^8 -C_{10}/R^{10}-\cdots`. The DFT-D2 [Grimme:2006:1787]_ variant
takes the explicit form below. Here, dispersion coefficients,
:math:`C_6^{ij}`, obtained from the geometric mean of tabulated elemental
values, are summed over interatomic distances, :math:`R_{ij}`, modulated
by a damping function, :math:`f_{damp}(R_{ij})`, that gradually activates
the dispersion correction (at a rate characterized by :math:`\alpha_6`)
over a distance characterized by the sum of the two atomic vdW radii,
:math:`R_{vdW}`, while an overall scaling term, :math:`s_6`, is optimized
to be unique to each :math:`E_{xc}` functional. (:math:`\alpha_6` is
sometimes allowed to vary as well.)

.. math:: E_{disp}^{\text{D2}}=-s_6 \sum_{i,j>i}^{N_{at}} \frac{C_6^{ij}}{(R_{ij})^6} f_{damp}(R_{ij})
   :label: DFTD2

.. math:: f_{damp}(R_{ij}) = \frac{1}{1 + e^{- \alpha_6 (R_{ij}/R_{vdW} - 1)}}

Grimme recently presented a refined method, DFT-D3, [Grimme:2010:154104]_
which incorporates an additional :math:`R^{-8}` term in the dispersion
series and adjusts the :math:`C_{6}^{ij}` combination formula and damping
function. The individual atomic :math:`C_6^i` are interpolated from
several reference values based upon coordination numbers extracted from
the molecular structure, rather than assigned solely by atomic identity as
in DFT-D2, and thereby incorporate some awareness of the chemical
environment into an otherwise largely heuristic correction.  The -D3
dispersion has the following form, where :math:`s_{r,6}` and :math:`s_8`
are the customary nonunity parameters fitted for individual functionals.

.. math:: E_{disp}^{\text{D3ZERO}}=-\sum_{n=6,8} s_n \sum_{i,j>i}^{N_{at}} 
   \frac{C_n^{ij}}{(R_{ij})^n} f_{damp}(R_{ij})
   :label: DFTD3ZERO

.. math:: f_{damp}(R_{ij}) = \frac{1}{1 + 6 (R_{ij}/(s_{r,n} R_0^{ij}))^{- \alpha_n}}


A modified damping scheme for DFT-D3 using the rational damping form of
Becke and Johnson was introduced in [Grimme:2011:1456]_.  The parameters
fit for individual functionals are now :math:`s_6`, :math:`s_8`,
:math:`a_1`, and :math:`a_2`.

.. math:: E_{disp}^{\text{D3BJ}}=-\sum_{n=6,8} s_n \sum_{i,j>i}^{N_{at}} 
   \frac{C_n^{ij}}{(R_{ij})^n + (f_{damp})^n}

.. math:: f_{damp} = a_1 \sqrt{\frac{C_8^{ij}}{C_6^{ij}}} + a_2

All parameters characterizing the dispersion correction are taken from
`http://toc.uni-muenster.de/DFTD3/ <http://toc.uni-muenster.de/DFTD3/>`_
or else from the literature.

Running DFTD3
~~~~~~~~~~~~~

A number of *a posteriori* dispersion corrections are available in
|PSIfour|.  While most are computed within |PSIfours| codebase (-D1, -D2,
-CHG, -DAS2009, -DAS2010), the -D3 correction and its variants are
available only through the ``DFTD3`` program.  Once installed, the
``dftd3``/|PSIfour| interface is transparent, and all corrections are
interfaced exactly alike.

Dispersion corrections are built into DFT functionals, so appending an *a
posteriori* correction to a computation is as simple as
``energy('b2plyp-d')`` *vs.* ``energy('b2plyp')``. For example, the
following input file computes (with much redundant work) for water a
B3LYP, a B3LYP-D2, and a B3LYP-D3 (zero-damping) energy. ::

   molecule h2o {
        O
        H 1 1.0
        H 1 1.0 2 104.5
    }
    set {
        basis cc-pVDZ
    }
    energy('b3lyp')
    energy('b3lyp-d')
    energy('b3lyp-d3')

Consult the table :ref:`-D Functionals <table:dft_disp>` to see for each
functional what corrections are available and what default parameters
define them. The dispersion correction is available after a calculation in
the PSI variable :psivar:`DISPERSION CORRECTION ENERGY <DISPERSIONCORRECTIONENERGY>`. 
By default, the output from the ``dftd3``
program is suppressed; to see it in the output file, set print > 2.


.. _`table:dashd`:

.. table:: Variants of S. Grimme's -D correction

   +-------------------------------------+---------------------------------------------------------------------------------------------------------+---------------------------------------------------------------+
   | Extension [#f0]_                    | Variant and Computing Program                                                                           | |scf__dft_dispersion_parameters|                              |
   +=====================================+=========================================================================================================+===============================================================+
   | -D                                  | alias to -D2P4                                                                                          |                                                               |
   +-------------------------------------+---------------------------------------------------------------------------------------------------------+---------------------------------------------------------------+
   | -D1                                 | -D1 [#f1]_ within |PSIfour|                                                                             |                                                               |
   +-------------------------------------+---------------------------------------------------------------------------------------------------------+---------------------------------------------------------------+
   | -D2                                 | alias to -D2P4                                                                                          |                                                               |
   +-------------------------------------+---------------------------------------------------------------------------------------------------------+---------------------------------------------------------------+
   | -D2P4                               | -D2 [#f2]_ within |PSIfour|                                                                             | [:math:`s_6`]                                                 |
   +-------------------------------------+---------------------------------------------------------------------------------------------------------+---------------------------------------------------------------+
   | -D2GR                               | -D2 [#f2]_ through ``dftd3``                                                                            | [:math:`s_6`, :math:`\alpha_6`]                               |
   +-------------------------------------+---------------------------------------------------------------------------------------------------------+---------------------------------------------------------------+
   | -D3                                 | alias to -D3ZERO                                                                                        |                                                               |
   +-------------------------------------+---------------------------------------------------------------------------------------------------------+---------------------------------------------------------------+
   | -D3ZERO                             | -D3 [#f3]_ w/ original zero-damping through ``dftd3``                                                   | [:math:`s_6`, :math:`s_8`, :math:`s_{r,6}`, :math:`\alpha_6`] |
   +-------------------------------------+---------------------------------------------------------------------------------------------------------+---------------------------------------------------------------+
   | -D3BJ                               | -D3 [#f4]_ w/ newer Becke-Johnson rational damping through ``dftd3``                                    | [:math:`s_6`, :math:`s_8`, :math:`a_1`, :math:`a_2`]          |
   +-------------------------------------+---------------------------------------------------------------------------------------------------------+---------------------------------------------------------------+

.. rubric:: Footnotes

.. [#f0] Note that there are functionals with these extensions (*e.g.*, wB97X-D) that, 
   not being Grimme corrections, have nothing to do with this table.
   
.. [#f1] [Grimme:2004:1463]_
.. [#f2] [Grimme:2006:1787]_
.. [#f3] [Grimme:2010:154104]_
.. [#f4] [Grimme:2011:1456]_

A few practical examples:

* DFT-D2 single point with default parameters (``dftd3`` not called) ::

   energy('bp86-d')

* DFT-D3BJ optimization with default parameters ::

   optimize('pbe-d3bj')

* DFT-D2 optimization with custom s6 parameter ::

   set dft_dispersion_parameters [1.20]
   optimize('b3lyp-d2')

* DFT-D3ZERO single point (b3lyp) with custom s8 parameter (reset all four values) ::

   set dft_dispersion_parameters [1.0, 2.0, 1.261, 14.0]
   energy('b3lyp-d3')

If only dispersion corrections (rather than total energies) are of
interest, the ``dftd3`` program can be run independently of the scf
through the python function :py:func:`~molutil.run_dftd3`. (This function
is the same |PSIfour|/``dftd3`` interface that is called during an scf job.)
This route is much faster than running a DFT-D energy.

* Some set-up::

   molecule nene {
   Ne
   Ne 1 2.0
   }
   
   nene.update_geometry()

* The same four dispersion corrections/gradients as the section above::

   >>> print nene.run_dftd3('bp86', 'd', dertype=0)
   -7.735e-05
   
   >>> E, G = nene.run_dftd3('pbe', 'd3bj')
   >>> print G
   [[0.0, 0.0, -1.1809087569358e-05], [0.0, 0.0, 1.1809087569358e-05]]
   
   >>> E, G = nene.run_dftd3('b3lyp', 'd2', {'s6': 1.20})
   >>> print E
   -8.84e-05
   
   >>> E, G = nene.run_dftd3(dashlvl='d3', dashparam={'s8': 2.0, 'alpha6': 14.0, 'sr6': 1.261, 's6': 1.0})
   >>> print E
   -0.00024762

.. autofunction:: molutil.run_dftd3


.. comment print_stdout('  -D correction from Py-side')
.. comment eneyne.update_geometry()
.. comment E, G = eneyne.run_dftd3('b3lyp', 'd2gr')
.. comment compare_values(ref_d2[0], E, 7, 'Ethene-Ethyne -D2')
.. comment mA = eneyne.extract_subsets(1)
.. comment E, G = mA.run_dftd3('b3lyp', 'd2gr')
.. comment compare_values(ref_d2[1], E, 7, 'Ethene -D2')
.. comment mB = eneyne.extract_subsets(2)
.. comment E, G = mB.run_dftd3('b3lyp', 'd2gr')
.. comment compare_values(ref_d2[2], E, 7, 'Ethyne -D2')
.. comment #mBcp = eneyne.extract_subsets(2,1)
.. comment #E, G = mBcp.run_dftd3('b3lyp', 'd2gr')
.. comment #compare_values(ref_d2[2], E, 7, 'Ethyne(CP) -D2')
.. comment 
.. comment E, G = eneyne.run_dftd3('b3lyp', 'd3zero')
.. comment compare_values(ref_d3zero[0], E, 7, 'Ethene-Ethyne -D3 (zero)')
.. comment mA = eneyne.extract_subsets(1)
.. comment E, G = mA.run_dftd3('b3lyp', 'd3zero')
.. comment compare_values(ref_d3zero[1], E, 7, 'Ethene -D3 (zero)')
.. comment mB = eneyne.extract_subsets(2)
.. comment E, G = mB.run_dftd3('b3lyp', 'd3zero')
.. comment compare_values(ref_d3zero[2], E, 7, 'Ethyne -D3 (zero)')
.. comment 
.. comment E, G = eneyne.run_dftd3('b3lyp', 'd3bj')
.. comment compare_values(ref_d3bj[0], E, 7, 'Ethene-Ethyne -D3 (bj)')
.. comment mA = eneyne.extract_subsets(1)
.. comment E, G = mA.run_dftd3('b3lyp', 'd3bj')
.. comment compare_values(ref_d3bj[1], E, 7, 'Ethene -D3 (bj)')
.. comment mB = eneyne.extract_subsets(2)
.. comment E, G = mB.run_dftd3('b3lyp', 'd3bj')
.. comment compare_values(ref_d3bj[2], E, 7, 'Ethyne -D3 (bj)')
.. comment 
.. comment E, G = eneyne.run_dftd3('b3lyp', 'd3')
.. comment compare_values(ref_d3zero[0], E, 7, 'Ethene-Ethyne -D3 (alias)')
.. comment E, G = eneyne.run_dftd3('b3lyp', 'd')
.. comment compare_values(ref_d2[0], E, 7, 'Ethene-Ethyne -D (alias)')
.. comment E, G = eneyne.run_dftd3('b3lyp', 'd2')
.. comment compare_values(ref_d2[0], E, 7, 'Ethene-Ethyne -D2 (alias)')
.. comment 
.. comment set basis sto-3g
.. comment set scf_type df
.. comment set dft_radial_points 50  # use really bad grid for speed since all we want is the -D value
.. comment set dft_spherical_points 110
.. comment #set scf print 3  # will print dftd3 program output to psi4 output file
.. comment 
.. comment 
.. comment print_stdout('  -D correction from C-side')
.. comment activate(mA)
.. comment energy('b3lyp-d2p4')
.. comment compare_values(ref_d2[1], get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D2 (calling psi4 Disp class)')
.. comment energy('b3lyp-d2gr')
.. comment compare_values(ref_d2[1], get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D2 (calling dftd3 -old)')
.. comment energy('b3lyp-d3zero')
.. comment compare_values(ref_d3zero[1], get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D3 (calling dftd3 -zero)')
.. comment energy('b3lyp-d3bj')
.. comment compare_values(ref_d3bj[1], get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D3 (calling dftd3 -bj)')
.. comment 
.. comment energy('b3lyp-d2')
.. comment compare_values(ref_d2[1], get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D2 (alias)')
.. comment energy('b3lyp-d3')
.. comment compare_values(ref_d3zero[1], get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D3 (alias)')
.. comment energy('b3lyp-d')
.. comment compare_values(ref_d2[1], get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D (alias)')
.. comment energy('wb97x-d')
.. comment compare_values(-0.000834247063, get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene wb97x-d (chg)')
.. comment 
.. comment print_stdout('  non-default -D correction from C-side')
.. comment activate(mB)
.. comment set dft_dispersion_parameters [0.75]
.. comment energy('b3lyp-d2p4')
.. comment compare_values(ref_pbe_d2[2], get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D2 (calling psi4 Disp class)')
.. comment set dft_dispersion_parameters [0.75, 20.0]
.. comment energy('b3lyp-d2gr')
.. comment compare_values(ref_pbe_d2[2], get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D2 (calling dftd3 -old)')
.. comment set dft_dispersion_parameters [1.0,  0.722, 1.217, 14.0]
.. comment energy('b3lyp-d3zero')
.. comment compare_values(ref_pbe_d3zero[2], get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D3 (calling dftd3 -zero)')
.. comment set dft_dispersion_parameters [1.000, 0.7875, 0.4289, 4.4407]
.. comment energy('b3lyp-d3bj')
.. comment compare_values(ref_pbe_d3bj[2], get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D3 (calling dftd3 -bj)')
.. comment 
.. comment set dft_dispersion_parameters [0.75]
.. comment energy('b3lyp-d2')
.. comment compare_values(ref_pbe_d2[2], get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D2 (alias)')
.. comment set dft_dispersion_parameters [1.0,  0.722, 1.217, 14.0]
.. comment energy('b3lyp-d3')
.. comment compare_values(ref_pbe_d3zero[2], get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D3 (alias)')
.. comment set dft_dispersion_parameters [0.75]
.. comment energy('b3lyp-d')
.. comment compare_values(ref_pbe_d2[2], get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D (alias)')
.. comment activate(mA)
.. comment set dft_dispersion_parameters [1.0]
.. comment energy('wb97x-d')
.. comment compare_values(-0.000834247063, get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene wb97x-d (chg)')
.. comment 
.. comment print_stdout('  non-default -D correction from Py-side')
.. comment eneyne.update_geometry()
.. comment eneyne.run_dftd3('b3lyp', 'd2gr', {'s6': 0.75})
.. comment compare_values(ref_pbe_d2[0], get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene-Ethyne -D2')
.. comment mA = eneyne.extract_subsets(1)
.. comment mA.run_dftd3('b3lyp', 'd2gr', {'s6': 0.75})
.. comment compare_values(ref_pbe_d2[1], get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D2')
.. comment mB = eneyne.extract_subsets(2)
.. comment mB.run_dftd3('b3lyp', 'd2gr', {'s6': 0.75})
.. comment compare_values(ref_pbe_d2[2], get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethyne -D2')
.. comment 
.. comment eneyne.run_dftd3('b3lyp', 'd3zero', {'s6': 1.0,  's8': 0.722, 'sr6': 1.217, 'alpha6': 14.0})
.. comment compare_values(ref_pbe_d3zero[0], get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene-Ethyne -D3 (zero)')
.. comment mA = eneyne.extract_subsets(1)
.. comment mA.run_dftd3('b3lyp', 'd3zero', {'s6': 1.0,  's8': 0.722, 'sr6': 1.217, 'alpha6': 14.0})
.. comment compare_values(ref_pbe_d3zero[1], get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D3 (zero)')
.. comment mB = eneyne.extract_subsets(2)
.. comment mB.run_dftd3('b3lyp', 'd3zero', {'s6': 1.0,  's8': 0.722, 'sr6': 1.217, 'alpha6': 14.0})
.. comment compare_values(ref_pbe_d3zero[2], get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethyne -D3 (zero)')
.. comment 
.. comment eneyne.run_dftd3('b3lyp', 'd3bj', {'s6': 1.000, 's8':  0.7875, 'a1':  0.4289, 'a2': 4.4407})
.. comment compare_values(ref_pbe_d3bj[0], get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene-Ethyne -D3 (bj)')
.. comment mA = eneyne.extract_subsets(1)
.. comment mA.run_dftd3('b3lyp', 'd3bj', {'s6': 1.000, 's8':  0.7875, 'a1':  0.4289, 'a2': 4.4407})
.. comment compare_values(ref_pbe_d3bj[1], get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D3 (bj)')
.. comment mB = eneyne.extract_subsets(2)
.. comment mB.run_dftd3('b3lyp', 'd3bj', {'s6': 1.000, 's8':  0.7875, 'a1':  0.4289, 'a2': 4.4407})
.. comment compare_values(ref_pbe_d3bj[2], get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethyne -D3 (bj)')
.. comment 
.. comment eneyne.run_dftd3('b3lyp', 'd3', {'s6': 1.0,  's8': 0.722, 'sr6': 1.217, 'alpha6': 14.0})
.. comment compare_values(ref_pbe_d3zero[0], get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene-Ethyne -D3 (alias)')
.. comment eneyne.run_dftd3('b3lyp', 'd', {'s6': 0.75})
.. comment compare_values(ref_pbe_d2[0], get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene-Ethyne -D (alias)')
.. comment eneyne.run_dftd3('b3lyp', 'd2', {'s6': 0.75})
.. comment compare_values(ref_pbe_d2[0], get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene-Ethyne -D2 (alias)')


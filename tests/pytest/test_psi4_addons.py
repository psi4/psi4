import psi4
import pytest


using_ambit = pytest.mark.skipif(psi4.addons("ambit") is False,
                                reason="Psi4 not compiled with ambit. Rebuild with -DENABLE_ambit")
using_cfour = pytest.mark.skipif(psi4.addons("cfour") is False,
                                reason="Psi4 not detecting CFOUR. Add `xcfour` to envvar PSIPATH or PATH")
using_chemps2 = pytest.mark.skipif(psi4.addons("chemps2") is False,
                                reason="Psi4 not compiled with CheMPS2. Rebuild with -DENABLE_CheMPS2")
using_dftd3 = pytest.mark.skipif(psi4.addons("dftd3") is False,
                                reason="Psi4 not detecting DFTD3. Add `dftd3` to envvar PSIPATH or PATH")
using_dkh = pytest.mark.skipif(psi4.addons("dkh") is False,
                                reason="Psi4 not compiled with dkh. Rebuild with -DENABLE_dkh")
using_libefp = pytest.mark.skipif(psi4.addons("libefp") is False,
                                reason="Psi4 not compiled with libefp. Rebuild with -DENABLE_libefp")
using_erd = pytest.mark.skipif(psi4.addons("erd") is False,
                                reason="Psi4 not compiled with erd. Rebuild with -DENABLE_erd")
using_gdma = pytest.mark.skipif(psi4.addons("gdma") is False,
                                reason="Psi4 not compiled with gdma. Rebuild with -DENABLE_gdma")
using_mrcc = pytest.mark.skipif(psi4.addons("mrcc") is False,
                                reason="Psi4 not detecting MRCC. Add `dmrcc` to envvar PSIPATH or PATH")
using_pcmsolver = pytest.mark.skipif(psi4.addons("pcmsolver") is False,
                                reason="Psi4 not compiled with PCMSolver. Rebuild with -DENABLE_PCMSolver")
using_simint = pytest.mark.skipif(psi4.addons("simint") is False,
                                reason="Psi4 not compiled with simint. Rebuild with -DENABLE_simint")


@using_gdma
def test_gdma():
    """gdma1"""
    #! Water RHF/cc-pVTZ distributed multipole analysis

    ref_energy = -76.0571685433842219
    ref_dma_mat = psi4.core.Matrix(3, 9)
    ref_dma_mat.name = 'Reference DMA values'
    ref_dma_arr = [
      [ -0.43406697290168, -0.18762673939633,  0.00000000000000,  0.00000000000000,  0.03206686487531,
         0.00000000000000, -0.00000000000000, -0.53123477172696,  0.00000000000000 ],
      [  0.21703348903257, -0.06422316619952,  0.00000000000000, -0.11648289410022,  0.01844320206227,
         0.00000000000000,  0.07409226544133, -0.07115302332866,  0.00000000000000 ],
      [  0.21703348903257, -0.06422316619952,  0.00000000000000,  0.11648289410022,  0.01844320206227,
         0.00000000000000, -0.07409226544133, -0.07115302332866,  0.00000000000000 ]
    ]
    for i in range(3):
        for j in range(9):
            ref_dma_mat.set(i, j, ref_dma_arr[i][j])
    ref_tot_mat = psi4.core.Matrix(1, 9)
    ref_tot_mat.name = "Reference total values"
    ref_tot_arr = [
         0.00000000516346, -0.79665315928128,  0.00000000000000,  0.00000000000000,  0.10813259329390,
         0.00000000000000,  0.00000000000000, -2.01989585894142,  0.00000000000000
    ]
    for i in range(9):
        ref_tot_mat.set(0, i, ref_tot_arr[i])

    # noreorient/nocom are not needed, but are used here to guarantee that the
    #   GDMA origin placement defined below is at the O atom.
    water = psi4.geometry("""
        O  0.000000  0.000000  0.117176
        H -0.000000 -0.756950 -0.468706
        H -0.000000  0.756950 -0.468706
     noreorient
     nocom
    """)

    psi4.set_options({"scf_type": "pk",
                       "basis": "cc-pvtz",
                       "d_convergence": 10,
                       "gdma_switch": 0,
                       "gdma_radius": [ "H", 0.65 ],
                       "gdma_limit": 2,
                       "gdma_origin": [ 0.000000,  0.000000,  0.117176 ]})

    energy, wfn = psi4.energy('scf', return_wfn=True)

    psi4.gdma(wfn)
    dmavals = psi4.core.get_array_variable("DMA DISTRIBUTED MULTIPOLES")
    totvals = psi4.core.get_array_variable("DMA TOTAL MULTIPOLES")
    ans = psi4.compare_values(ref_energy, energy, 8, "SCF Energy")
    assert ans
    ans = psi4.compare_matrices(dmavals, ref_dma_mat, 6, "DMA Distributed Multipoles")
    assert ans
    ans = psi4.compare_matrices(totvals, ref_tot_mat, 6, "DMA Total Multipoles")
    assert ans


@using_mrcc
def test_mrcc():
    """mrcc/ccsdt"""
    #! CCSDT cc-pVDZ energy for the H2O molecule using MRCC

    h2o = psi4.geometry("""
        o
        h 1 1.0
        h 1 1.0 2 104.5
    """)

    psi4.set_options({
        'basis': 'cc-pvdz',
        'freeze_core': 'true'})

    psi4.energy('mrccsdt')

    ans = psi4.compare_values(  8.801465529972, psi4.get_variable("NUCLEAR REPULSION ENERGY"), 6, 'NRE')
    assert ans
    ans = psi4.compare_values(-76.021418445155, psi4.get_variable("SCF TOTAL ENERGY"), 6, 'SCF')
    assert ans
    ans = psi4.compare_values( -0.204692406830, psi4.get_variable("MP2 CORRELATION ENERGY") , 6, 'MP2 correlation')
    assert ans
    ans = psi4.compare_values( -0.217715210258, psi4.get_variable("CCSDT CORRELATION ENERGY"), 6, 'CCSDT correlation')
    assert ans
    ans = psi4.compare_values(-76.239133655413, psi4.get_variable("CURRENT ENERGY"), 6, 'CCSDT')
    assert ans


@using_chemps2
def test_chemps2():
    """chemps2/scf-n2"""
    #! dmrg-scf on N2

    N2 = psi4.geometry("""
      N       0.0000   0.0000   0.0000
      N       0.0000   0.0000   2.1180
    units au
    """)

    psi4.set_options({
    'basis': 'cc-pVDZ',
    'reference': 'rhf',
    'e_convergence': 1e-12,
    'd_convergence': 1e-12,

    'dmrg_irrep': 0,
    'dmrg_multiplicity': 1,
    'restricted_docc': [ 1 , 0 , 0 , 0 , 0 , 1 , 0 , 0 ],
    'active': [ 2 , 0 , 1 , 1 , 0 , 2 , 1 , 1 ],

    'dmrg_sweep_states': [   500,  1000,  1000 ],
    'dmrg_sweep_energy_conv': [ 1e-10, 1e-10, 1e-10 ],
    'dmrg_sweep_dvdson_rtol': [  1e-4,  1e-6,  1e-8 ],
    'dmrg_sweep_max_sweeps': [     5,     5,    10 ],
    'dmrg_sweep_noise_prefac': [  0.05,  0.05,   0.0 ],
    'dmrg_print_corr': True,
    'dmrg_mps_write': False,

    'dmrg_unitary_write': True,
    'dmrg_diis': True,
    'dmrg_scf_diis_thr': 1e-2,
    'dmrg_diis_write': True,

    'dmrg_excitation': 0,   # Ground state
    'dmrg_scf_state_avg': False,
    'dmrg_scf_active_space': 'NO',  # INPUT; NO; LOC
    'dmrg_local_init': False,
    })

    psi4.energy("dmrg-scf")

    ref_energy = -109.1035023353
    ans = psi4.compare_values(ref_energy, psi4.get_variable("CURRENT ENERGY"), 6, "DMRG Energy")
    assert ans


@using_dftd3
def test_dftd3():
    """dftd3/energy"""
    #! Exercises the various DFT-D corrections, both through python directly and through c++

    ref_d2         = [-0.00390110, -0.00165271, -0.00058118]
    ref_d3zero     = [-0.00285088, -0.00084340, -0.00031923]
    ref_d3bj       = [-0.00784595, -0.00394347, -0.00226683]

    ref_pbe_d2     = [-0.00278650, -0.00118051, -0.00041513]
    ref_pbe_d3zero = [-0.00175474, -0.00045421, -0.00016839]
    ref_pbe_d3bj   = [-0.00475937, -0.00235265, -0.00131239]

    eneyne = psi4.geometry("""
    C   0.000000  -0.667578  -2.124659
    C   0.000000   0.667578  -2.124659
    H   0.923621  -1.232253  -2.126185
    H  -0.923621  -1.232253  -2.126185
    H  -0.923621   1.232253  -2.126185
    H   0.923621   1.232253  -2.126185
    --
    C   0.000000   0.000000   2.900503
    C   0.000000   0.000000   1.693240
    H   0.000000   0.000000   0.627352
    H   0.000000   0.000000   3.963929
    """)

    psi4.print_stdout('  -D correction from Py-side')
    eneyne.update_geometry()
    E, G = eneyne.run_dftd3('b3lyp', 'd2gr')
    ans = psi4.compare_values(ref_d2[0], E, 7, 'Ethene-Ethyne -D2')
    assert ans
    mA = eneyne.extract_subsets(1)
    E, G = mA.run_dftd3('b3lyp', 'd2gr')
    ans = psi4.compare_values(ref_d2[1], E, 7, 'Ethene -D2')
    assert ans
    mB = eneyne.extract_subsets(2)
    E, G = mB.run_dftd3('b3lyp', 'd2gr')
    ans = psi4.compare_values(ref_d2[2], E, 7, 'Ethyne -D2')
    assert ans
    #mBcp = eneyne.extract_subsets(2,1)
    #E, G = mBcp.run_dftd3('b3lyp', 'd2gr')
    #compare_values(ref_d2[2], E, 7, 'Ethyne(CP) -D2')

    E, G = eneyne.run_dftd3('b3lyp', 'd3zero')
    ans = psi4.compare_values(ref_d3zero[0], E, 7, 'Ethene-Ethyne -D3 (zero)')
    assert ans
    mA = eneyne.extract_subsets(1)
    E, G = mA.run_dftd3('b3lyp', 'd3zero')
    ans = psi4.compare_values(ref_d3zero[1], E, 7, 'Ethene -D3 (zero)')
    assert ans
    mB = eneyne.extract_subsets(2)
    E, G = mB.run_dftd3('b3lyp', 'd3zero')
    ans = psi4.compare_values(ref_d3zero[2], E, 7, 'Ethyne -D3 (zero)')
    assert ans

    E, G = eneyne.run_dftd3('b3lyp', 'd3bj')
    ans = psi4.compare_values(ref_d3bj[0], E, 7, 'Ethene-Ethyne -D3 (bj)')
    assert ans
    mA = eneyne.extract_subsets(1)
    E, G = mA.run_dftd3('b3lyp', 'd3bj')
    ans = psi4.compare_values(ref_d3bj[1], E, 7, 'Ethene -D3 (bj)')
    assert ans
    mB = eneyne.extract_subsets(2)
    E, G = mB.run_dftd3('b3lyp', 'd3bj')
    ans = psi4.compare_values(ref_d3bj[2], E, 7, 'Ethyne -D3 (bj)')
    assert ans

    E, G = eneyne.run_dftd3('b3lyp', 'd3')
    ans = psi4.compare_values(ref_d3zero[0], E, 7, 'Ethene-Ethyne -D3 (alias)')
    assert ans
    E, G = eneyne.run_dftd3('b3lyp', 'd')
    ans = psi4.compare_values(ref_d2[0], E, 7, 'Ethene-Ethyne -D (alias)')
    assert ans
    E, G = eneyne.run_dftd3('b3lyp', 'd2')
    ans = psi4.compare_values(ref_d2[0], E, 7, 'Ethene-Ethyne -D2 (alias)')
    assert ans

    psi4.set_options({'basis': 'sto-3g',
                      'scf_type': 'df',
                      'dft_radial_points': 50,  # use really bad grid for speed since all we want is the -D value
                      'dft_spherical_points': 110,
                      #'scf print': 3,  # will print dftd3 program output to psi4 output file
                    })

    psi4.print_stdout('  -D correction from C-side')
    psi4.activate(mA)
    psi4.energy('b3lyp-d2p4')
    ans = psi4.compare_values(ref_d2[1], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D2 (calling psi4 Disp class)')
    assert ans
    psi4.energy('b3lyp-d2gr')
    ans = psi4.compare_values(ref_d2[1], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D2 (calling dftd3 -old)')
    assert ans
    psi4.energy('b3lyp-d3zero')
    ans = psi4.compare_values(ref_d3zero[1], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D3 (calling dftd3 -zero)')
    assert ans
    psi4.energy('b3lyp-d3bj')
    ans = psi4.compare_values(ref_d3bj[1], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D3 (calling dftd3 -bj)')
    assert ans

    psi4.energy('b3lyp-d2')
    ans = psi4.compare_values(ref_d2[1], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D2 (alias)')
    assert ans
    psi4.energy('b3lyp-d3')
    ans = psi4.compare_values(ref_d3zero[1], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D3 (alias)')
    assert ans
    psi4.energy('b3lyp-d')
    ans = psi4.compare_values(ref_d2[1], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D (alias)')
    assert ans
    psi4.energy('wb97x-d')
    ans = psi4.compare_values(-0.000834247063, psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene wb97x-d (chg)')
    assert ans

    psi4.print_stdout('  non-default -D correction from C-side')
    psi4.activate(mB)
    psi4.set_options({'dft_dispersion_parameters': [0.75]})
    psi4.energy('b3lyp-d2p4')
    ans = psi4.compare_values(ref_pbe_d2[2], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D2 (calling psi4 Disp class)')
    assert ans
    psi4.set_options({'dft_dispersion_parameters': [0.75, 20.0]})
    psi4.energy('b3lyp-d2gr')
    ans = psi4.compare_values(ref_pbe_d2[2], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D2 (calling dftd3 -old)')
    assert ans
    psi4.set_options({'dft_dispersion_parameters': [1.0,  0.722, 1.217, 14.0]})
    psi4.energy('b3lyp-d3zero')
    ans = psi4.compare_values(ref_pbe_d3zero[2], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D3 (calling dftd3 -zero)')
    assert ans
    psi4.set_options({'dft_dispersion_parameters': [1.000, 0.7875, 0.4289, 4.4407]})
    psi4.energy('b3lyp-d3bj')
    ans = psi4.compare_values(ref_pbe_d3bj[2], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D3 (calling dftd3 -bj)')
    assert ans

    psi4.set_options({'dft_dispersion_parameters': [0.75]})
    psi4.energy('b3lyp-d2')
    ans = psi4.compare_values(ref_pbe_d2[2], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D2 (alias)')
    assert ans
    psi4.set_options({'dft_dispersion_parameters': [1.0,  0.722, 1.217, 14.0]})
    psi4.energy('b3lyp-d3')
    ans = psi4.compare_values(ref_pbe_d3zero[2], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D3 (alias)')
    assert ans
    psi4.set_options({'dft_dispersion_parameters': [0.75]})
    psi4.energy('b3lyp-d')
    ans = psi4.compare_values(ref_pbe_d2[2], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D (alias)')
    assert ans
    psi4.activate(mA)
    psi4.set_options({'dft_dispersion_parameters': [1.0]})
    psi4.energy('wb97x-d')
    ans = psi4.compare_values(-0.000834247063, psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene wb97x-d (chg)')
    assert ans

    psi4.print_stdout('  non-default -D correction from Py-side')
    eneyne.update_geometry()
    eneyne.run_dftd3('b3lyp', 'd2gr', {'s6': 0.75})
    ans = psi4.compare_values(ref_pbe_d2[0], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene-Ethyne -D2')
    assert ans
    mA = eneyne.extract_subsets(1)
    mA.run_dftd3('b3lyp', 'd2gr', {'s6': 0.75})
    ans = psi4.compare_values(ref_pbe_d2[1], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D2')
    assert ans
    mB = eneyne.extract_subsets(2)
    mB.run_dftd3('b3lyp', 'd2gr', {'s6': 0.75})
    ans = psi4.compare_values(ref_pbe_d2[2], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethyne -D2')
    assert ans

    eneyne.run_dftd3('b3lyp', 'd3zero', {'s6': 1.0,  's8': 0.722, 'sr6': 1.217, 'alpha6': 14.0})
    ans = psi4.compare_values(ref_pbe_d3zero[0], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene-Ethyne -D3 (zero)')
    assert ans
    mA = eneyne.extract_subsets(1)
    mA.run_dftd3('b3lyp', 'd3zero', {'s6': 1.0,  's8': 0.722, 'sr6': 1.217, 'alpha6': 14.0})
    ans = psi4.compare_values(ref_pbe_d3zero[1], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D3 (zero)')
    assert ans
    mB = eneyne.extract_subsets(2)
    mB.run_dftd3('b3lyp', 'd3zero', {'s6': 1.0,  's8': 0.722, 'sr6': 1.217, 'alpha6': 14.0})
    ans = psi4.compare_values(ref_pbe_d3zero[2], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethyne -D3 (zero)')
    assert ans

    eneyne.run_dftd3('b3lyp', 'd3bj', {'s6': 1.000, 's8':  0.7875, 'a1':  0.4289, 'a2': 4.4407})
    ans = psi4.compare_values(ref_pbe_d3bj[0], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene-Ethyne -D3 (bj)')
    assert ans
    mA = eneyne.extract_subsets(1)
    mA.run_dftd3('b3lyp', 'd3bj', {'s6': 1.000, 's8':  0.7875, 'a1':  0.4289, 'a2': 4.4407})
    ans = psi4.compare_values(ref_pbe_d3bj[1], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D3 (bj)')
    assert ans
    mB = eneyne.extract_subsets(2)
    mB.run_dftd3('b3lyp', 'd3bj', {'s6': 1.000, 's8':  0.7875, 'a1':  0.4289, 'a2': 4.4407})
    ans = psi4.compare_values(ref_pbe_d3bj[2], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethyne -D3 (bj)')
    assert ans
    eneyne.run_dftd3('b3lyp', 'd3', {'s6': 1.0,  's8': 0.722, 'sr6': 1.217, 'alpha6': 14.0})

    ans = psi4.compare_values(ref_pbe_d3zero[0], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene-Ethyne -D3 (alias)')
    assert ans
    eneyne.run_dftd3('b3lyp', 'd', {'s6': 0.75})
    ans = psi4.compare_values(ref_pbe_d2[0], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene-Ethyne -D (alias)')
    assert ans
    eneyne.run_dftd3('b3lyp', 'd2', {'s6': 0.75})
    ans = psi4.compare_values(ref_pbe_d2[0], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene-Ethyne -D2 (alias)')
    assert ans


@using_libefp
def test_libefp():
    """libefp/qchem-qmefp-sp"""
    #! EFP on mixed QM (water) and EFP (water + 2 * ammonia) system.
    #! An EFP-only calc performed first to test vales against q-chem.

    qmefp = psi4.geometry("""
    # QM fragment
    0 1
    units bohr
    O1     0.000000000000     0.000000000000     0.224348285559
    H2    -1.423528800232     0.000000000000    -0.897393142237
    H3     1.423528800232     0.000000000000    -0.897393142237
    # EFP as EFP fragments
    --
    efp h2o -4.014110144291     2.316749370493    -1.801514729931 -2.902133 1.734999 -1.953647
    --
    efp NH3,1.972094713645,,3.599497221584 ,    5.447701074734 -1.105309 2.033306 -1.488582
    --
    efp NH3 -7.876296399270    -1.854372164887    -2.414804197762  2.526442 1.658262 -2.742084
    """)

    #  <<<  EFP calc  >>>
    psi4.set_options({
        'basis': '6-31g*',
        'scf_type': 'pk',
        'guess': 'core',
        'df_scf_guess': False})

    psi4.energy('efp')
    ans = psi4.compare_values( 9.1793879214, qmefp.nuclear_repulsion_energy(), 6, 'QM NRE')
    assert ans
    ans = psi4.compare_values(-0.0004901368, psi4.get_variable('efp elst energy'), 6, 'EFP-EFP Elst')  # from q-chem
    assert ans
    ans = psi4.compare_values(-0.0003168768, psi4.get_variable('efp ind energy'), 6, 'EFP-EFP Indc')
    assert ans
    ans = psi4.compare_values(-0.0021985285, psi4.get_variable('efp disp energy'), 6, 'EFP-EFP Disp')  # from q-chem
    assert ans
    ans = psi4.compare_values( 0.0056859871, psi4.get_variable('efp exch energy'), 6, 'EFP-EFP Exch')  # from q-chem
    assert ans
    ans = psi4.compare_values( 0.0026804450, psi4.get_variable('efp total energy'), 6, 'EFP-EFP Totl')
    assert ans
    ans = psi4.compare_values( 0.0026804450, psi4.get_variable('current energy'), 6, 'Current')
    assert ans
    psi4.core.print_variables()

    psi4.core.clean()
    psi4.core.clean_variables()

    #  <<<  QM + EFP calc  >>>
    psi4.set_options({
        'e_convergence': 12,
        'd_convergence': 12})
    psi4.energy('scf')

    ans = psi4.compare_values( 9.1793879214, qmefp.nuclear_repulsion_energy(), 6, 'QM NRE')
    assert ans
    ans = psi4.compare_values( 0.2622598847, psi4.get_variable('efp total energy') - psi4.get_variable('efp ind energy'), 6, 'EFP corr to SCF')  # from q-chem
    assert ans
    ans = psi4.compare_values(-0.0117694790, psi4.get_variable('efp ind energy'), 6, 'QM-EFP Indc')  # from q-chem
    assert ans
    ans = psi4.compare_values(-0.0021985285, psi4.get_variable('efp disp energy'), 6, 'EFP-EFP Disp')  # from q-chem
    assert ans
    ans = psi4.compare_values( 0.0056859871, psi4.get_variable('efp exch energy'), 6, 'EFP-EFP Exch')  # from q-chem
    assert ans
    ans = psi4.compare_values( 0.2504904057, psi4.get_variable('efp total energy'), 6, 'EFP-EFP Totl')  # from q-chem
    assert ans
    ans = psi4.compare_values(-76.0139362744, psi4.get_variable('scf total energy'), 6, 'SCF')  # from q-chem
    assert ans
    psi4.core.print_variables()


@using_pcmsolver
def test_pcmsolver():
    """pcmsolver/scf"""
    #! pcm

    nucenergy   =  12.0367196636183458
    polenergy   =  -0.0053060443528559
    totalenergy = -55.4559426361734040

    NH3 = psi4.geometry("""
    symmetry c1
    N     -0.0000000001    -0.1040380466      0.0000000000
    H     -0.9015844116     0.4818470201     -1.5615900098
    H     -0.9015844116     0.4818470201      1.5615900098
    H      1.8031688251     0.4818470204      0.0000000000
    units bohr
    no_reorient
    no_com
    """)

    psi4.set_options({
      'basis': 'STO-3G',
      'scf_type': 'pk',
      'pcm': True,
      'pcm_scf_type': 'total',
    })

    psi4.pcm_helper("""
       Units = Angstrom
       Medium {
       SolverType = IEFPCM
       Solvent = Water
       }

       Cavity {
       RadiiSet = UFF
       Type = GePol
       Scaling = False
       Area = 0.3
       Mode = Implicit
       }
    """)

    print('RHF-PCM, total algorithm')
    energy_scf1 = psi4.energy('scf')
    ans = psi4.compare_values(NH3.nuclear_repulsion_energy(), nucenergy, 10, "Nuclear repulsion energy (PCM, total algorithm)")
    assert ans
    ans = psi4.compare_values(energy_scf1, totalenergy, 10, "Total energy (PCM, total algorithm)")
    assert ans
    ans = psi4.compare_values(psi4.get_variable("PCM POLARIZATION ENERGY"), polenergy, 6, "Polarization energy (PCM, total algorithm)")
    assert ans

    psi4.set_options({'pcm_scf_type': 'separate'})
    print('RHF-PCM, separate algorithm')
    energy_scf2 = psi4.energy('scf')
    ans = psi4.compare_values(energy_scf2, totalenergy, 10, "Total energy (PCM, separate algorithm)")
    assert ans
    ans = psi4.compare_values(psi4.get_variable("PCM POLARIZATION ENERGY"), polenergy, 6, "Polarization energy (PCM, separate algorithm)")
    assert ans

    # Now force use of UHF on NH3 to check sanity of the algorithm with PCM
    psi4.set_options({'pcm_scf_type': 'total', 'reference': 'uhf'})
    print('UHF-PCM, total algorithm')
    energy_scf3 = psi4.energy('scf')
    ans = psi4.compare_values(energy_scf3, totalenergy, 10, "Total energy (PCM, separate algorithm)")
    assert ans
    ans = psi4.compare_values(psi4.get_variable("PCM POLARIZATION ENERGY"), polenergy, 6, "Polarization energy (PCM, separate algorithm)")
    assert ans


@using_erd
def test_erd():
    """erd/scf5"""
    #! Test of all different algorithms and reference types for SCF, on singlet and triplet O2, using the cc-pVTZ basis set and using ERD integrals.

    psi4.set_options({'integral_package': 'ERD'})

    psi4.print_stdout(' Case Study Test of all SCF algorithms/spin-degeneracies: Singlet-Triplet O2')

    #Ensure that the checkpoint file is always nuked
    psi4.core.IOManager.shared_object().set_specific_retention(32,False)

    Eref_nuc      =   30.78849213614545
    Eref_sing_can = -149.58723684929720
    Eref_sing_df  = -149.58715054487624
    Eref_uhf_can  = -149.67135517240553
    Eref_uhf_df   = -149.67125624291961
    Eref_rohf_can = -149.65170765757173
    Eref_rohf_df  = -149.65160796208073

    singlet_o2 = psi4.geometry("""
        0 1
        O
        O 1 1.1
        units    angstrom
    """)

    triplet_o2 = psi4.geometry("""
        0 3
        O
        O 1 1.1
        units    angstrom
    """)
    singlet_o2.update_geometry()
    triplet_o2.update_geometry()

    psi4.print_stdout('   -Nuclear Repulsion:')
    ans = psi4.compare_values(Eref_nuc, triplet_o2.nuclear_repulsion_energy(), 9, "Triplet nuclear repulsion energy")
    assert ans
    ans = psi4.compare_values(Eref_nuc, singlet_o2.nuclear_repulsion_energy(), 9, "Singlet nuclear repulsion energy")
    assert ans

    psi4.set_options({
        'basis': 'cc-pvtz',
        'df_basis_scf': 'cc-pvtz-jkfit',
        'print': 2})

    print('    -Singlet RHF:')
    psi4.set_module_options('scf', {'reference': 'rhf'})

    psi4.set_module_options('scf', {'scf_type': 'pk'})
    E = psi4.energy('scf', molecule=singlet_o2)
    ans = psi4.compare_values(Eref_sing_can, E, 6, 'Singlet PK RHF energy')
    assert ans

    psi4.set_module_options('scf', {'scf_type': 'direct'})
    E = psi4.energy('scf', molecule=singlet_o2)
    ans = psi4.compare_values(Eref_sing_can, E, 6, 'Singlet Direct RHF energy')
    assert ans

    psi4.set_module_options('scf', {'scf_type': 'out_of_core'})
    E = psi4.energy('scf', molecule=singlet_o2)
    ans = psi4.compare_values(Eref_sing_can, E, 6, 'Singlet Disk RHF energy')
    assert ans

    psi4.set_module_options('scf', {'scf_type': 'df'})
    E = psi4.energy('scf', molecule=singlet_o2)
    ans = psi4.compare_values(Eref_sing_df, E, 6, 'Singlet DF RHF energy')
    assert ans

    print('    -Singlet UHF:')
    psi4.set_module_options('scf', {'reference': 'uhf'})

    psi4.set_module_options('scf', {'scf_type': 'pk'})
    E = psi4.energy('scf', molecule=singlet_o2)
    ans = psi4.compare_values(Eref_sing_can, E, 6, 'Singlet PK UHF energy')
    assert ans

    psi4.set_module_options('scf', {'scf_type': 'direct'})
    E = psi4.energy('scf', molecule=singlet_o2)
    ans = psi4.compare_values(Eref_sing_can, E, 6, 'Singlet Direct UHF energy')
    assert ans

    psi4.set_module_options('scf', {'scf_type': 'out_of_core'})
    E = psi4.energy('scf', molecule=singlet_o2)
    ans = psi4.compare_values(Eref_sing_can, E, 6, 'Singlet Disk UHF energy')
    assert ans

    psi4.set_module_options('scf', {'scf_type': 'df'})
    E = psi4.energy('scf', molecule=singlet_o2)
    ans = psi4.compare_values(Eref_sing_df, E, 6, 'Singlet DF UHF energy')
    assert ans

    print('    -Singlet CUHF:')
    psi4.set_module_options('scf', {'reference': 'cuhf'})

    psi4.set_module_options('scf', {'scf_type': 'pk'})
    E = psi4.energy('scf', molecule=singlet_o2)
    ans = psi4.compare_values(Eref_sing_can, E, 6, 'Singlet PK CUHF energy')
    assert ans

    psi4.set_module_options('scf', {'scf_type': 'direct'})
    E = psi4.energy('scf', molecule=singlet_o2)
    ans = psi4.compare_values(Eref_sing_can, E, 6, 'Singlet Direct CUHF energy')
    assert ans

    psi4.set_module_options('scf', {'scf_type': 'out_of_core'})
    E = psi4.energy('scf', molecule=singlet_o2)
    ans = psi4.compare_values(Eref_sing_can, E, 6, 'Singlet Disk CUHF energy')
    assert ans

    psi4.set_module_options('scf', {'scf_type': 'df'})
    E = psi4.energy('scf', molecule=singlet_o2)
    ans = psi4.compare_values(Eref_sing_df, E, 6, 'Singlet DF CUHF energy')
    assert ans

    psi4.set_options({
        'basis': 'cc-pvtz',
        'df_basis_scf': 'cc-pvtz-jkfit',
        'guess': 'core',
        'print': 2})

    print('    -Triplet UHF:')
    psi4.set_module_options('scf', {'reference': 'uhf'})

    psi4.set_module_options('scf', {'scf_type': 'pk'})
    E = psi4.energy('scf', molecule=triplet_o2)
    ans = psi4.compare_values(Eref_uhf_can, E, 6, 'Triplet PK UHF energy')
    assert ans

    psi4.set_module_options('scf', {'scf_type': 'direct'})
    E = psi4.energy('scf', molecule=triplet_o2)
    ans = psi4.compare_values(Eref_uhf_can, E, 6, 'Triplet Direct UHF energy')
    assert ans

    psi4.set_module_options('scf', {'scf_type': 'out_of_core'})
    E = psi4.energy('scf', molecule=triplet_o2)
    ans = psi4.compare_values(Eref_uhf_can, E, 6, 'Triplet Disk UHF energy')
    assert ans

    psi4.set_module_options('scf', {'scf_type': 'df'})
    E = psi4.energy('scf', molecule=triplet_o2)
    ans = psi4.compare_values(Eref_uhf_df, E, 6, 'Triplet DF UHF energy')
    assert ans
    psi4.core.clean()

    print('    -Triplet ROHF:')
    psi4.set_module_options('scf', {'reference': 'rohf'})

    psi4.set_module_options('scf', {'scf_type': 'pk'})
    E = psi4.energy('scf', molecule=triplet_o2)
    ans = psi4.compare_values(Eref_rohf_can, E, 6, 'Triplet PK ROHF energy')
    assert ans
    psi4.core.clean()

    psi4.set_module_options('scf', {'scf_type': 'direct'})
    E = psi4.energy('scf', molecule=triplet_o2)
    ans = psi4.compare_values(Eref_rohf_can, E, 6, 'Triplet Direct ROHF energy')
    assert ans
    psi4.core.clean()

    psi4.set_module_options('scf', {'scf_type': 'out_of_core'})
    E = psi4.energy('scf', molecule=triplet_o2)
    ans = psi4.compare_values(Eref_rohf_can, E, 6, 'Triplet Disk ROHF energy')
    assert ans
    psi4.core.clean()

    psi4.set_module_options('scf', {'scf_type': 'df'})
    E = psi4.energy('scf', molecule=triplet_o2)
    ans = psi4.compare_values(Eref_rohf_df, E, 6, 'Triplet DF ROHF energy')
    assert ans
    psi4.core.clean()

    print('    -Triplet CUHF:')
    psi4.set_module_options('scf', {'reference': 'cuhf'})

    psi4.set_module_options('scf', {'scf_type': 'pk'})
    E = psi4.energy('scf', molecule=triplet_o2)
    ans = psi4.compare_values(Eref_rohf_can, E, 6, 'Triplet PK CUHF energy')
    assert ans
    psi4.core.clean()

    psi4.set_module_options('scf', {'scf_type': 'direct'})
    E = psi4.energy('scf', molecule=triplet_o2)
    ans = psi4.compare_values(Eref_rohf_can, E, 6, 'Triplet Direct CUHF energy')
    assert ans
    psi4.core.clean()

    psi4.set_module_options('scf', {'scf_type': 'out_of_core'})
    E = psi4.energy('scf', molecule=triplet_o2)
    ans = psi4.compare_values(Eref_rohf_can, E, 6, 'Triplet Disk CUHF energy')
    assert ans
    psi4.core.clean()

    psi4.set_module_options('scf', {'scf_type': 'df'})
    E = psi4.energy('scf', molecule=triplet_o2)
    ans = psi4.compare_values(Eref_rohf_df, E, 6, 'Triplet DF CUHF energy')
    assert ans


def test_json():
    """json/energy"""

    import numpy as np
    #import psi4
    #from psi4.driver import p4util

    # Generate JSON data
    json_data = {}
    json_data["molecule"] = """He 0 0 0\n--\nHe 0 0 1"""
    json_data["driver"] = "gradient"
    json_data["args"] = 'SCF'
    json_data["kwargs"] = {}
    json_data["options"] = {"BASIS": "STO-3G"}
    json_data["return_output"] = True

    psi4.json_wrapper.run_json(json_data)

    ans = psi4.compare_strings("STO-3G", json_data["options"]["BASIS"], "Options test")
    assert ans
    ans = psi4.compare_integers(True, json_data["success"], "Success")
    assert ans

    bench_energy = -5.433191881443323
    cenergy = json_data["variables"]["CURRENT ENERGY"]

    assert ans

    bench_gradient = np.array([[  0.0 , 0.0 ,   0.4206844],
                               [  0.0 , 0.0 ,  -0.4206844]])
    cgradient = psi4.core.Matrix.from_serial(json_data["return_value"])
    ans = psi4.compare_arrays(bench_gradient, cgradient.np, 4, "SCF RETURN_VALUE")
    assert ans

    return_wfn = "return_wfn" not in json_data["kwargs"]
    ans = psi4.compare_integers(True, return_wfn, "Immutable input")
    assert ans

    with open("output.dat", "w") as f:
        f.write(json_data["output"])


@using_cfour
def test_cfour():
    """cfour/sp-rhf-ccsd_t_"""
    #! single-point CCSD(T)/qz2p on water

    print('        <<< Literal ZMAT to Cfour >>>')

    psi4.cfour_helper("""
    CCSD(T) single-point energy calculation
    O
    H 1 R
    H 1 R 2 A

    R=0.958
    A=104.5

    *ACES2(CALC=CCSD(T),BASIS=qz2p
    SCF_CONV=12
    CC_CONV=12
    MEMORY=20000000)
    """)

    psi4.energy('cfour')

    ans = psi4.compare_values(-76.062748460117, psi4.get_variable('scf total energy'), 6, 'SCF')
    assert ans
    ans = psi4.compare_values(-76.332940127333, psi4.get_variable('mp2 total energy'), 6, 'MP2')
    assert ans
    ans = psi4.compare_values(-76.338453951890, psi4.get_variable('ccsd total energy'), 6, 'CCSD')
    assert ans
    ans = psi4.compare_values(-0.275705491773, psi4.get_variable('ccsd correlation energy'), 6, 'CCSD corl')
    assert ans
    ans = psi4.compare_values(-76.345717549886, psi4.get_variable('ccsd(t) total energy'), 6, 'CCSD(T)')
    assert ans
    ans = psi4.compare_values(-0.282969089769, psi4.get_variable('ccsd(t) correlation energy'), 6, 'CCSD(T) corl')
    assert ans

    psi4.core.clean()
    psi4.core.clean_variables()
    psi4.cfour_helper('')

    print('        <<< Translation of ZMAT to Psi4 format to Cfour >>>')

    psi4.geometry("""
    O
    H 1 R
    H 1 R 2 A

    R=0.958
    A=104.5
    """)

    psi4.set_options({
    'cfour_CALC_level': 'CCSD(T)',
    'cfour_BASIS': 'qz2p',
    'cfour_SCF_CONV': 12,
    'cfour_CC_CONV': 12,
    })

    psi4.energy('cfour')

    ans = psi4.compare_values(-76.062748460117, psi4.get_variable('scf total energy'), 6, 'SCF')
    assert ans
    ans = psi4.compare_values(-76.332940127333, psi4.get_variable('mp2 total energy'), 6, 'MP2')
    assert ans
    ans = psi4.compare_values(-76.338453951890, psi4.get_variable('ccsd total energy'), 6, 'CCSD')
    assert ans
    ans = psi4.compare_values(-0.275705491773, psi4.get_variable('ccsd correlation energy'), 6, 'CCSD corl')
    assert ans
    ans = psi4.compare_values(-76.345717549886, psi4.get_variable('ccsd(t) total energy'), 6, 'CCSD(T)')
    assert ans
    ans = psi4.compare_values(-0.282969089769, psi4.get_variable('ccsd(t) correlation energy'), 6, 'CCSD(T) corl')
    assert ans

    psi4.core.clean()
    psi4.core.clean_variables()
    psi4.cfour_helper('')

    print('        <<< Translation of ZMAT to Psi4 format to Cfour >>>')

    psi4.geometry("""
    O
    H 1 R
    H 1 R 2 A

    R=0.958
    A=104.5
    """)

    psi4.set_options({
    'cfour_CALC_level': 'CCSD(T)',
    'cfour_BASIS': 'qz2p',
    'cfour_SCF_CONV': 12,
    'cfour_CC_CONV': 12,
    })

    psi4.energy('cfour')

    ans = psi4.compare_values(-76.062748460117, psi4.get_variable('scf total energy'), 6, 'SCF')
    assert ans
    ans = psi4.compare_values(-76.332940127333, psi4.get_variable('mp2 total energy'), 6, 'MP2')
    assert ans
    ans = psi4.compare_values(-76.338453951890, psi4.get_variable('ccsd total energy'), 6, 'CCSD')
    assert ans
    ans = psi4.compare_values(-0.275705491773, psi4.get_variable('ccsd correlation energy'), 6, 'CCSD corl')
    assert ans
    ans = psi4.compare_values(-76.345717549886, psi4.get_variable('ccsd(t) total energy'), 6, 'CCSD(T)')
    assert ans
    ans = psi4.compare_values(-0.282969089769, psi4.get_variable('ccsd(t) correlation energy'), 6, 'CCSD(T) corl')
    assert ans


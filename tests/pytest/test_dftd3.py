import psi4

def test_dftd3():
    """dftd3/energy"""
    #! Exercises the various DFT-D corrections, both through python directly and through c++

    ref_d2         = [-0.00390110, -0.00165271, -0.00058118]  #TEST
    ref_d3zero     = [-0.00285088, -0.00084340, -0.00031923]  #TEST
    ref_d3bj       = [-0.00784595, -0.00394347, -0.00226683]  #TEST

    ref_pbe_d2     = [-0.00278650, -0.00118051, -0.00041513]  #TEST
    ref_pbe_d3zero = [-0.00175474, -0.00045421, -0.00016839]  #TEST
    ref_pbe_d3bj   = [-0.00475937, -0.00235265, -0.00131239]  #TEST

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

    psi4.print_stdout('  -D correction from Py-side')                      #TEST
    eneyne.update_geometry()
    E, G = eneyne.run_dftd3('b3lyp', 'd2gr')
    ans = psi4.compare_values(ref_d2[0], E, 7, 'Ethene-Ethyne -D2')              #TEST
    assert ans
    mA = eneyne.extract_subsets(1)
    E, G = mA.run_dftd3('b3lyp', 'd2gr')
    ans = psi4.compare_values(ref_d2[1], E, 7, 'Ethene -D2')                     #TEST
    assert ans
    mB = eneyne.extract_subsets(2)
    E, G = mB.run_dftd3('b3lyp', 'd2gr')
    ans = psi4.compare_values(ref_d2[2], E, 7, 'Ethyne -D2')                     #TEST
    assert ans
    #mBcp = eneyne.extract_subsets(2,1)                               #TEST
    #E, G = mBcp.run_dftd3('b3lyp', 'd2gr')                           #TEST
    #compare_values(ref_d2[2], E, 7, 'Ethyne(CP) -D2')                #TEST

    E, G = eneyne.run_dftd3('b3lyp', 'd3zero')
    ans = psi4.compare_values(ref_d3zero[0], E, 7, 'Ethene-Ethyne -D3 (zero)')   #TEST
    assert ans
    mA = eneyne.extract_subsets(1)
    E, G = mA.run_dftd3('b3lyp', 'd3zero')
    ans = psi4.compare_values(ref_d3zero[1], E, 7, 'Ethene -D3 (zero)')          #TEST
    assert ans
    mB = eneyne.extract_subsets(2)
    E, G = mB.run_dftd3('b3lyp', 'd3zero')
    ans = psi4.compare_values(ref_d3zero[2], E, 7, 'Ethyne -D3 (zero)')          #TEST
    assert ans

    E, G = eneyne.run_dftd3('b3lyp', 'd3bj')
    ans = psi4.compare_values(ref_d3bj[0], E, 7, 'Ethene-Ethyne -D3 (bj)')       #TEST
    assert ans
    mA = eneyne.extract_subsets(1)
    E, G = mA.run_dftd3('b3lyp', 'd3bj')
    ans = psi4.compare_values(ref_d3bj[1], E, 7, 'Ethene -D3 (bj)')              #TEST
    assert ans
    mB = eneyne.extract_subsets(2)
    E, G = mB.run_dftd3('b3lyp', 'd3bj')
    ans = psi4.compare_values(ref_d3bj[2], E, 7, 'Ethyne -D3 (bj)')              #TEST
    assert ans

    E, G = eneyne.run_dftd3('b3lyp', 'd3')
    ans = psi4.compare_values(ref_d3zero[0], E, 7, 'Ethene-Ethyne -D3 (alias)')  #TEST
    assert ans
    E, G = eneyne.run_dftd3('b3lyp', 'd')
    ans = psi4.compare_values(ref_d2[0], E, 7, 'Ethene-Ethyne -D (alias)')       #TEST
    assert ans
    E, G = eneyne.run_dftd3('b3lyp', 'd2')
    ans = psi4.compare_values(ref_d2[0], E, 7, 'Ethene-Ethyne -D2 (alias)')      #TEST
    assert ans

    psi4.set_options({'basis': 'sto-3g',
                      'scf_type': 'df',
                      'dft_radial_points': 50,  # use really bad grid for speed since all we want is the -D value
                      'dft_spherical_points': 110,
                      #'scf print': 3,  # will print dftd3 program output to psi4 output file
                    })

    psi4.print_stdout('  -D correction from C-side')                                                                         #TEST
    psi4.activate(mA)
    psi4.energy('b3lyp-d2p4')
    ans = psi4.compare_values(ref_d2[1], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D2 (calling psi4 Disp class)')  #TEST
    assert ans
    psi4.energy('b3lyp-d2gr')
    ans = psi4.compare_values(ref_d2[1], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D2 (calling dftd3 -old)')       #TEST
    assert ans
    psi4.energy('b3lyp-d3zero')
    ans = psi4.compare_values(ref_d3zero[1], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D3 (calling dftd3 -zero)')  #TEST
    assert ans
    psi4.energy('b3lyp-d3bj')
    ans = psi4.compare_values(ref_d3bj[1], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D3 (calling dftd3 -bj)')      #TEST
    assert ans

    psi4.energy('b3lyp-d2')
    ans = psi4.compare_values(ref_d2[1], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D2 (alias)')                    #TEST
    assert ans
    psi4.energy('b3lyp-d3')
    ans = psi4.compare_values(ref_d3zero[1], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D3 (alias)')                #TEST
    assert ans
    psi4.energy('b3lyp-d')
    ans = psi4.compare_values(ref_d2[1], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D (alias)')                     #TEST
    assert ans
    psi4.energy('wb97x-d')
    ans = psi4.compare_values(-0.000834247063, psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene wb97x-d (chg)')            #TEST
    assert ans

    psi4.print_stdout('  non-default -D correction from C-side')                                                                 #TEST
    psi4.activate(mB)
    psi4.set_options({'dft_dispersion_parameters': [0.75]})
    psi4.energy('b3lyp-d2p4')
    ans = psi4.compare_values(ref_pbe_d2[2], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D2 (calling psi4 Disp class)')  #TEST
    assert ans
    psi4.set_options({'dft_dispersion_parameters': [0.75, 20.0]})
    psi4.energy('b3lyp-d2gr')
    ans = psi4.compare_values(ref_pbe_d2[2], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D2 (calling dftd3 -old)')       #TEST
    assert ans
    psi4.set_options({'dft_dispersion_parameters': [1.0,  0.722, 1.217, 14.0]})
    psi4.energy('b3lyp-d3zero')
    ans = psi4.compare_values(ref_pbe_d3zero[2], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D3 (calling dftd3 -zero)')  #TEST
    assert ans
    psi4.set_options({'dft_dispersion_parameters': [1.000, 0.7875, 0.4289, 4.4407]})
    psi4.energy('b3lyp-d3bj')
    ans = psi4.compare_values(ref_pbe_d3bj[2], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D3 (calling dftd3 -bj)')      #TEST
    assert ans

    psi4.set_options({'dft_dispersion_parameters': [0.75]})
    psi4.energy('b3lyp-d2')
    ans = psi4.compare_values(ref_pbe_d2[2], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D2 (alias)')                    #TEST
    assert ans
    psi4.set_options({'dft_dispersion_parameters': [1.0,  0.722, 1.217, 14.0]})
    psi4.energy('b3lyp-d3')
    ans = psi4.compare_values(ref_pbe_d3zero[2], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D3 (alias)')                #TEST
    assert ans
    psi4.set_options({'dft_dispersion_parameters': [0.75]})
    psi4.energy('b3lyp-d')
    ans = psi4.compare_values(ref_pbe_d2[2], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D (alias)')                     #TEST
    assert ans
    psi4.activate(mA)
    psi4.set_options({'dft_dispersion_parameters': [1.0]})
    psi4.energy('wb97x-d')
    ans = psi4.compare_values(-0.000834247063, psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene wb97x-d (chg)')                #TEST
    assert ans

    psi4.print_stdout('  non-default -D correction from Py-side')                                                         #TEST
    eneyne.update_geometry()
    eneyne.run_dftd3('b3lyp', 'd2gr', {'s6': 0.75})
    ans = psi4.compare_values(ref_pbe_d2[0], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene-Ethyne -D2')              #TEST
    assert ans
    mA = eneyne.extract_subsets(1)
    mA.run_dftd3('b3lyp', 'd2gr', {'s6': 0.75})
    ans = psi4.compare_values(ref_pbe_d2[1], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D2')                     #TEST
    assert ans
    mB = eneyne.extract_subsets(2)
    mB.run_dftd3('b3lyp', 'd2gr', {'s6': 0.75})
    ans = psi4.compare_values(ref_pbe_d2[2], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethyne -D2')                     #TEST
    assert ans

    eneyne.run_dftd3('b3lyp', 'd3zero', {'s6': 1.0,  's8': 0.722, 'sr6': 1.217, 'alpha6': 14.0})
    ans = psi4.compare_values(ref_pbe_d3zero[0], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene-Ethyne -D3 (zero)')   #TEST
    assert ans
    mA = eneyne.extract_subsets(1)
    mA.run_dftd3('b3lyp', 'd3zero', {'s6': 1.0,  's8': 0.722, 'sr6': 1.217, 'alpha6': 14.0})
    ans = psi4.compare_values(ref_pbe_d3zero[1], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D3 (zero)')          #TEST
    assert ans
    mB = eneyne.extract_subsets(2)
    mB.run_dftd3('b3lyp', 'd3zero', {'s6': 1.0,  's8': 0.722, 'sr6': 1.217, 'alpha6': 14.0})
    ans = psi4.compare_values(ref_pbe_d3zero[2], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethyne -D3 (zero)')          #TEST
    assert ans

    eneyne.run_dftd3('b3lyp', 'd3bj', {'s6': 1.000, 's8':  0.7875, 'a1':  0.4289, 'a2': 4.4407})
    ans = psi4.compare_values(ref_pbe_d3bj[0], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene-Ethyne -D3 (bj)')       #TEST
    assert ans
    mA = eneyne.extract_subsets(1)
    mA.run_dftd3('b3lyp', 'd3bj', {'s6': 1.000, 's8':  0.7875, 'a1':  0.4289, 'a2': 4.4407})
    ans = psi4.compare_values(ref_pbe_d3bj[1], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D3 (bj)')              #TEST
    assert ans
    mB = eneyne.extract_subsets(2)
    mB.run_dftd3('b3lyp', 'd3bj', {'s6': 1.000, 's8':  0.7875, 'a1':  0.4289, 'a2': 4.4407})
    ans = psi4.compare_values(ref_pbe_d3bj[2], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethyne -D3 (bj)')              #TEST
    assert ans
    eneyne.run_dftd3('b3lyp', 'd3', {'s6': 1.0,  's8': 0.722, 'sr6': 1.217, 'alpha6': 14.0})

    ans = psi4.compare_values(ref_pbe_d3zero[0], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene-Ethyne -D3 (alias)')  #TEST
    assert ans
    eneyne.run_dftd3('b3lyp', 'd', {'s6': 0.75})
    ans = psi4.compare_values(ref_pbe_d2[0], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene-Ethyne -D (alias)')       #TEST
    assert ans
    eneyne.run_dftd3('b3lyp', 'd2', {'s6': 0.75})
    ans = psi4.compare_values(ref_pbe_d2[0], psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene-Ethyne -D2 (alias)')      #TEST
    assert ans


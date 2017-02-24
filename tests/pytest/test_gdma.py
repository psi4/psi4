import psi4


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


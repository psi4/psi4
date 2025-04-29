import pytest

from utils import compare_values, compare
from addons import using

import psi4

pytestmark = [pytest.mark.psi, pytest.mark.api]

@pytest.fixture
def mols():
    return {
        "h2o" : psi4.geometry("""
0 1
O    0.000000000000     0.000000000000    -0.124038860300
H    0.000000000000    -1.431430901356     0.984293362719
H    0.000000000000     1.431430901356     0.984293362719
units au
"""),
        "nh2" : psi4.geometry("""
0 2
N    0.000000000000000   0.000000000000000  -0.145912918634892
H    0.000000000000000  -1.511214298139000   1.013682596946108
H    0.000000000000000   1.511214298139000   1.013682596946108
units au
"""),
        "h2o_nap1" : psi4.geometry("""
0 1
O    0.000000000000     0.000000000000    -0.124038860300
H    0.000000000000    -1.431430901356     0.984293362719
H    0.000000000000     1.431430901356     0.984293362719
--
1 1
Na   0.000000000000     0.000000000000    -4.124038860300
units au
""")
        }

@pytest.mark.smoke
@pytest.mark.parametrize("j_algo", [ 
        pytest.param("DFDIRJ") 
    ]
) #to be extended in the future
@pytest.mark.parametrize("k_algo", [ 
        pytest.param("LINK"),
        pytest.param("COSX"),
        pytest.param("SNLINK", marks=using('gauxc')),
    ]
) #to be extended in the future
def test_composite_call(j_algo, k_algo, mols, request):
    """Test all SCF_TYPE={J}+{K} combinations for an HF calculation.
    The correct algorithm pair should be called in each case.""" 

    # initial option setup 
    test_id = request.node.callspec.id
    
    molecule = mols["h2o"]
    
    scf_type = f'{j_algo}+{k_algo}'
    psi4.set_options({ "scf_type" : f'{j_algo}+{k_algo}', "incfock": True, "save_jk": True })
  
    if any([ _ in scf_type for _ in ["COSX", "SNLINK"] ]): 
        psi4.set_options({ "screening" : "schwarz" })
    else:
        psi4.set_options({ "screening" : "density" })

    # run composite JK algorithm
    E, wfn = psi4.energy("hf/6-31g", molecule=molecule, return_wfn=True) 

    clean_j_name, clean_k_name = wfn.jk().name().split("+")

    # check that correct J algo has been called
    clean_j_name = clean_j_name.replace("-", "") # replace DF-DirJ with DFDirJ
    assert clean_j_name.lower() == j_algo.lower(), f'{test_id} has correct J build method'

    # check that correct K algo has been called
    clean_k_name = clean_k_name.replace("-", "") # replace sn-LinK with snLinK
    assert clean_k_name.lower() == k_algo.lower(), f'{test_id} has correct K build method'

@pytest.mark.parametrize(
    "inp",
    [
        pytest.param({"method" : "hf",
                      "options": {"reference" : "rhf"},
                      "molecule" : "h2o",
                      "bsse_type" : None,
                      },
                      marks=pytest.mark.quick,
                      id="h2o (rhf)"),
        pytest.param({"method" : "b3lyp",
                      "options": {"reference" : "rhf"},
                      "molecule" : "h2o",
                      "bsse_type" : None,
                      },
                      marks=pytest.mark.quick,
                      id="h2o (rks)"),
        pytest.param({"method" : "hf",
                      "options": {"reference" : "uhf"},
                      "molecule" : "nh2",
                      "bsse_type" : None,
                      },
                      id="nh2 (uhf)"),
        pytest.param({"method" : "hf",
                      "options": {"reference" : "rohf"},
                      "molecule" : "nh2",
                      "bsse_type" : None,
                      },
                      id="nh2 (rohf)"),
        pytest.param({"method" : "hf",
                      "options": {"reference" : "rhf"},
                      "molecule" : "h2o_nap1",
                      "bsse_type" : "CP",
                      },
                      marks=pytest.mark.nbody,
                      id="h2o/na+ (rhf ie)"),
    ],
)
@pytest.mark.parametrize(
    "scf",
    [
        pytest.param({"scf_type" : "dfdirj+cosx",
                      "ref" : { 
                          "h2o (rhf)" : -76.026780223322,
                          "h2o (rks)" : -76.420402720419,
                          "nh2 (uhf)" : -55.566890252551,
                          "nh2 (rohf)" : -55.562689948780,
                          "h2o/na+ (rhf ie)" : -0.040121884077,
                      },
                      },
                      id="cosx"),
        pytest.param({"scf_type" : "dfdirj+snlink",
                      "snlink_force_cartesian": False,
                      "ref" : { 
                          "h2o (rhf)" : -76.026788692185, 
                          "h2o (rks)" : -76.420403557357,
                          "nh2 (uhf)" : -55.566911357539,
                          "nh2 (rohf)" : -55.562710424257,
                          "h2o/na+ (rhf ie)" : -0.040118757043,
                      },
                      },
                      id="snlink (spherical)", marks=using("gauxc")),
        pytest.param({"scf_type" : "dfdirj+snlink",
                      "snlink_force_cartesian": True,
                      "ref" : { 
                          "h2o (rhf)" : -76.026788692185, 
                          "h2o (rks)" : -76.420403557357,
                          "nh2 (uhf)" : -55.566911357539,
                          "nh2 (rohf)" : -55.562710424257,
                          "h2o/na+ (rhf ie)" : -0.040118757043,
                      },
                      },
                      id="snlink (cartesian)", marks=using("gauxc")),
 
    ]
)
def test_seminum(inp, scf, mols, request):
    """Test the DF-DirJ + COSX/sn-LinK JK objects via SCF calculations"""

    test_id = request.node.callspec.id
    
    molecule = mols[inp["molecule"]]
    psi4.set_options({"scf_type" : scf["scf_type"], "basis": "cc-pvdz"})
    psi4.set_options(inp["options"])
    if "snlink_force_cartesian" in scf.keys():
        psi4.set_options({"snlink_force_cartesian": scf["snlink_force_cartesian"]})

        #SNLINK_FORCE_CARTESIAN doesnt work with symmetry currently
        molecule.reset_point_group("C1")

    # does the SCF energy match a pre-computed reference?
    energy_seminum = psi4.energy(inp["method"], molecule=molecule, bsse_type=inp["bsse_type"])
    assert compare_values(scf["ref"][test_id.split("-")[1]], energy_seminum, 6, f'{test_id} accurate to reference (1e-6 threshold)')

    # is the SCF energy reasonably close to a conventional SCF?
    psi4.set_options({"scf_type" : "pk"})
    energy_pk = psi4.energy(inp["method"], molecule=molecule, bsse_type=inp["bsse_type"])
    assert compare_values(energy_pk, energy_seminum, 4, f'{test_id} DFDIRJ+COSX accurate to PK (1e-4 threshold)')

@pytest.mark.parametrize(
    "inp",
    [
        pytest.param({"method" : "hf",
                      "options": {"reference" : "rhf"},
                      "molecule" : "h2o",
                      "bsse_type" : None,
                      },
                      id="h2o (rhf)"),
        pytest.param({"method" : "b3lyp",
                      "options": {"reference" : "rhf"},
                      "molecule" : "h2o",
                      "bsse_type" : None,
                      },
                      id="h2o (rks)"),
        pytest.param({"method" : "hf",
                      "options": {"reference" : "uhf"},
                      "molecule" : "nh2",
                      "bsse_type" : None,
                      },
                      id="nh2 (uhf)"),
        pytest.param({"method" : "hf",
                      "options": {"reference" : "rohf"},
                      "molecule" : "nh2",
                      "bsse_type" : None,
                      },
                      id="nh2 (rohf)"),
        pytest.param({"method" : "hf",
                      "options": {"reference" : "rhf"},
                      "molecule" : "h2o_nap1",
                      "bsse_type" : "CP",
                      },
                      marks=pytest.mark.nbody,
                      id="h2o/na+ (rhf ie)"),
    ],
)
@pytest.mark.parametrize(
    "scf",
    [
        pytest.param({"scf_type" : "dfdirj+cosx"},
                      id="cosx"),
        pytest.param({"scf_type" : "dfdirj+snlink"},
                      id="snlink", marks=using("gauxc")),
    ]
)
def test_seminum_incfock(inp, scf, mols, request):
    """Test the efficiency of IncFock in DF-DirJ + COSX/sn-LinK JK objects via SCF calculations"""

    test_id = request.node.callspec.id

    molecule = mols[inp["molecule"]]
    psi4.set_options({"scf_type" : scf["scf_type"], "basis": "cc-pvdz", "incfock": False})
    psi4.set_options(inp["options"])

    # compute energy+wfn without IncFock 
    energy_seminum_noinc, wfn_seminum_noinc = psi4.energy(inp["method"], molecule=molecule, bsse_type=inp["bsse_type"], return_wfn=True)
    #assert compare_values(inp["ref"], energy_dfjcosk, atol=1e-6)

    # compute energy+wfn with Incfock 
    psi4.set_options({"incfock" : True})
    energy_seminum_inc, wfn_seminum_inc = psi4.energy(inp["method"], molecule=molecule, bsse_type=inp["bsse_type"], return_wfn=True)

    # how do energies compare?
    assert compare_values(energy_seminum_noinc, energy_seminum_inc, 6, f'{test_id} IncFock accurate (1e-6 threshold)')
    
    # how do SCF iteration counts compare?
    niter_noinc = int(wfn_seminum_noinc.variable("SCF ITERATIONS")) if wfn_seminum_noinc.has_scalar_variable("SCF ITERATIONS") else 0
    niter_inc = int(wfn_seminum_inc.variable("SCF ITERATIONS")) if wfn_seminum_inc.has_scalar_variable("SCF ITERATIONS") else 0
    
    assert compare(True, abs(niter_inc - niter_noinc) <= 3, f'{test_id} IncFock efficient')

@pytest.mark.parametrize("functional", [ "bp86", "b3lyp" ])
@pytest.mark.parametrize(
    "scf_type", 
    [ 
        pytest.param("DFDIRJ"), 
        pytest.param("LINK"),
        pytest.param("COSX"),
        pytest.param("SNLINK", marks=using('gauxc')),
        pytest.param("DFDIRJ+LINK"),
        pytest.param("DFDIRJ+COSX"),
        pytest.param("DFDIRJ+SNLINK", marks=using('gauxc')),
    ]
)
def test_dfdirj(functional, scf_type, mols):
    """Test the functionality of the SCF_TYPE keyword for CompositeJK methods under varying situations:
      - Using hybrid DFT functionals without specifying a K algorithm should cause a RuntimeError to be thrown.
      - Not specifying a J algorithm should cause a ValidationError to be thrown."""

    composite_algo_to_matrix = {
        "DFDIRJ": "J",
        "LINK" : "K",
        "COSX": "K",
        "SNLINK": "K",
    }

    molecule = mols["h2o"]
    screening = "CSAM" if any([ _ in scf_type for _ in [ "COSX", "SNLINK" ] ]) else "DENSITY"
    
    # if J algorithm isn't specified, code should throw here...
    if not any([ algo in scf_type for algo, matrix in composite_algo_to_matrix.items() if matrix == "J" ]):
        with pytest.raises(psi4.ValidationError) as e_info:
            psi4.set_options({"scf_type": scf_type, "reference": "rhf", "basis": "cc-pvdz", "screening": screening}) 

        # we keep this line just for printout purposes; should always pass if done correctly 
        assert compare(type(e_info), pytest.ExceptionInfo, f'{scf_type}+{functional} throws ValidationError')
 
    # ... else we continue as normal 
    else:  
        psi4.set_options({"scf_type": scf_type, "reference": "rhf", "basis": "cc-pvdz", "screening": screening}) 
    
        is_hybrid = True if functional == "b3lyp" else False
        k_algo_specified = True if any([ algo in scf_type for algo, matrix in composite_algo_to_matrix.items() if matrix == "K" ]) else False

        # if K algorithm isn't specified, but hybrid functional is used, code should throw...
        if is_hybrid and not k_algo_specified: 
            with pytest.raises(RuntimeError) as e_info:
                E = psi4.energy(functional, molecule=molecule) 

            # we keep this line just for printout purposes; should always pass if done correctly 
            assert compare(type(e_info), pytest.ExceptionInfo, f'{scf_type}+{functional} throws RuntimeError')
    
        # ... else code will run fine
        else:
            E = psi4.energy(functional, molecule=molecule) 

            # we keep this line just for printout purposes; should always pass if done correctly 
            assert compare(type(E), float, f'{scf_type}+{functional} executes')

@pytest.mark.parametrize("j_algo", [ "DFDIRJ" ]) #to be extended in the future
@pytest.mark.parametrize(
    "k_algo", 
    [ 
        pytest.param("LINK"),
        pytest.param("COSX"),
        pytest.param("SNLINK", marks=using('gauxc')),
    ]
) 
@pytest.mark.parametrize("df_basis_scf", [ "CC-PVDZ-JKFIT", "DEF2-UNIVERSAL-JFIT" ]) #to be extended in the future
def test_j_algo_bp86(j_algo, k_algo, df_basis_scf, mols):
    """Test SCF_TYPE={J} and all SCF_TYPE={J}+{K} combinations for a BP86 calculation.
    They should all give the exact same answer (within tolerance).""" 

    molecule = mols["h2o"]
    
    # run base composite J algorithm 
    psi4.set_options({"scf_type" : j_algo, "basis": "cc-pvdz", "df_basis_scf": df_basis_scf})
    energy_dfdirj = psi4.energy("bp86", molecule=molecule) 
    
    # compare composite combinations to base J algorithm
    scf_type = j_algo + "+" + k_algo    
    screening = "CSAM" if any([ _ in scf_type for _ in [ "COSX", "SNLINK" ] ]) else "DENSITY"

    psi4.set_options({"scf_type" : scf_type, "reference": "rhf", "basis": "cc-pvdz", "df_basis_scf": df_basis_scf, "screening": screening})
    energy_composite = psi4.energy("bp86", molecule=molecule) 
 
    assert compare_values(energy_dfdirj, energy_composite, 6, f'BP86/{df_basis_scf} {scf_type} accurate to {j_algo} (1e-6 threshold)')

@pytest.mark.parametrize("inp", [
    pytest.param({"molecule" : "h2o",
                  "options": {"reference" : "rhf"},
                  "ref":   [[ 0.0, -0.0,        -0.01496538],
                            [ 0.0, -0.01039195,  0.007477  ],
                            [ 0.0,  0.01039195,  0.007477  ]]
                  },
                  id="h2o (rhf)"),
    pytest.param({"molecule" : "nh2",
                  "options": {"reference" : "uhf"},
                  "ref" :  [[ 0.0,  0.0,         0.01517592],
                            [ 0.0,  0.00610143, -0.00739917],
                            [ 0.0, -0.00610143, -0.00739917]]
                  },
                  id="nh2 (uhf)"),
])
def test_cosx_grad(inp, mols):
    """ Test the functionality of DFJ+COSX gradient calculations."""

    molecule = mols[inp["molecule"]]
    
    psi4.set_options({ "scf_type" : "DFDIRJ+COSX"})
    psi4.set_options(inp["options"])
    grad_cosx = psi4.gradient("hf/cc-pvdz", molecule=molecule) 

    assert compare_values(inp["ref"], grad_cosx, 6, "COSX gradients match reference to 1e-6")

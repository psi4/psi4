import pytest

from utils import compare_values, compare

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

@pytest.mark.parametrize(
    "inp",
    [
        pytest.param({"method" : "hf",
                      "options": {"reference" : "rhf"},
                      "molecule" : "h2o",
                      "bsse_type" : None,
                      "ref" : -76.026780223322},
                      id="h2o (rhf)"),
        pytest.param({"method" : "b3lyp",
                      "options": {"reference" : "rhf"},
                      "molecule" : "h2o",
                      "bsse_type" : None,
                      "ref" : -76.420402720419},
                      id="h2o (rks)"),
        pytest.param({"method" : "hf",
                      "options": {"reference" : "uhf"},
                      "molecule" : "nh2",
                      "bsse_type" : None,
                      "ref" : -55.566890252551},
                      id="nh2 (uhf)"),
        pytest.param({"method" : "hf",
                      "options": {"reference" : "rohf"},
                      "molecule" : "nh2",
                      "bsse_type" : None,
                      "ref" : -55.562689948780},
                      id="nh2 (rohf)"),
        pytest.param({"method" : "hf",
                      "options": {"reference" : "rhf"},
                      "molecule" : "h2o_nap1",
                      "bsse_type" : "CP",
                      "ref" :  -0.040121884077},
                      marks=pytest.mark.nbody,
                      id="h2o/na+ (rhf ie)"),
    ],
)
def test_dfjcosk(inp, mols, request):
    """Test the DFJCOSK JK object via SCF calculations"""

    test_id = request.node.callspec.id
    
    molecule = mols[inp["molecule"]]
    psi4.set_options({"scf_type" : "dfdirj+cosx", "basis": "cc-pvdz"})
    psi4.set_options(inp["options"])

    # does the DFJCOSK SCF energy match a pre-computed reference?
    energy_dfjcosk = psi4.energy(inp["method"], molecule=molecule, bsse_type=inp["bsse_type"])
    assert compare_values(inp["ref"], energy_dfjcosk, 6, f'{test_id} DFDIRJ+COSX accurate to reference (1e-6 threshold)')

    # is the DFJCOSK SCF energy reasonably close to a conventional SCF?
    psi4.set_options({"scf_type" : "pk"})
    energy_pk = psi4.energy(inp["method"], molecule=molecule, bsse_type=inp["bsse_type"])
    assert compare_values(energy_pk, energy_dfjcosk, 4, f'{test_id} DFDIRJ+COSX accurate to PK (1e-4 threshold)')

@pytest.mark.parametrize(
    "inp",
    [
        pytest.param({"method" : "hf",
                      "options": {"reference" : "rhf"},
                      "molecule" : "h2o",
                      "bsse_type" : None,
                      "ref" : -76.026780223322},
                      id="h2o (rhf)"),
        pytest.param({"method" : "b3lyp",
                      "options": {"reference" : "rhf"},
                      "molecule" : "h2o",
                      "bsse_type" : None,
                      "ref" : -76.420402720419},
                      id="h2o (rks)"),
        pytest.param({"method" : "hf",
                      "options": {"reference" : "uhf"},
                      "molecule" : "nh2",
                      "bsse_type" : None,
                      "ref" : -55.566890252551},
                      id="nh2 (uhf)"),
        pytest.param({"method" : "hf",
                      "options": {"reference" : "rohf"},
                      "molecule" : "nh2",
                      "bsse_type" : None,
                      "ref" : -55.562689948780},
                      id="nh2 (rohf)"),
        pytest.param({"method" : "hf",
                      "options": {"reference" : "rhf"},
                      "molecule" : "h2o_nap1",
                      "bsse_type" : "CP",
                      "ref" :  -0.040121884077},
                      marks=pytest.mark.nbody,
                      id="h2o/na+ (rhf ie)"),
    ],
)
def test_dfjcosk_incfock(inp, mols, request):
    """Test the efficiency of IncFock in DFJCOSK JK object via SCF calculations"""

    test_id = request.node.callspec.id

    molecule = mols[inp["molecule"]]
    psi4.set_options({"scf_type" : "dfdirj+cosx", "basis": "cc-pvdz", "incfock": False})
    psi4.set_options(inp["options"])

    # compute DFJCOSK energy+wfn without IncFock 
    energy_dfjcosk_noinc, wfn_dfjcosk_noinc = psi4.energy(inp["method"], molecule=molecule, bsse_type=inp["bsse_type"], return_wfn=True)
    #assert compare_values(inp["ref"], energy_dfjcosk, atol=1e-6)

    # compute DFJCOSK energy+wfn with Incfock 
    psi4.set_options({"incfock" : True})
    energy_dfjcosk_inc, wfn_dfjcosk_inc = psi4.energy(inp["method"], molecule=molecule, bsse_type=inp["bsse_type"], return_wfn=True)

    # how do energies compare?
    assert compare_values(energy_dfjcosk_noinc, energy_dfjcosk_inc, 6, f'{test_id} IncFock accurate (1e-6 threshold)')
    
    # how do SCF iteration counts compare?
    niter_noinc = int(wfn_dfjcosk_noinc.variable("SCF ITERATIONS")) if wfn_dfjcosk_noinc.has_scalar_variable("SCF ITERATIONS") else 0
    niter_inc = int(wfn_dfjcosk_inc.variable("SCF ITERATIONS")) if wfn_dfjcosk_inc.has_scalar_variable("SCF ITERATIONS") else 0
    
    assert compare(True, abs(niter_inc - niter_noinc) <= 3, f'{test_id} IncFock efficient')

@pytest.mark.parametrize("functional", [ "bp86", "b3lyp" ])
@pytest.mark.parametrize("scf_type", [ "DFDIRJ", "LINK", "COSX", "DFDIRJ+COSX", "DFDIRJ+LINK" ])
def test_dfdirj(functional, scf_type, mols):
    """Test the functionality of the SCF_TYPE keyword for CompositeJK methods under varying situations:
      - Using hybrid DFT functionals without specifying a K algorithm should cause a RuntimeError to be thrown.
      - Not specifying a J algorithm should cause a ValidationError to be thrown."""

    composite_J_algo = [ "DFDIRJ" ]
    composite_K_algo = [ "LINK", "COSX" ]

    molecule = mols["h2o"]
    screening = "CSAM" if "COSX" in scf_type else "DENSITY"
    
    # if J algorithm isn't specified, code should throw here...
    if not any([ piece in scf_type for piece in composite_J_algo ]): 
        with pytest.raises(psi4.ValidationError) as e_info:
            psi4.set_options({"scf_type": scf_type, "reference": "rhf", "basis": "cc-pvdz", "screening": screening}) 

        # we keep this line just for printout purposes; should always pass if done correctly 
        assert compare(type(e_info), pytest.ExceptionInfo, f'{scf_type}+{functional} throws ValidationError')
    
    # ... else options will set normally
    else:  
        psi4.set_options({"scf_type": scf_type, "reference": "rhf", "basis": "cc-pvdz", "screening": screening}) 
    
    is_hybrid = True if functional == "b3lyp" else False
    k_algo_specified = True if any([ piece in scf_type for piece in composite_K_algo ]) else False

    # if K algorithm isn't specified, but hybrid functional is used, code should throw...
    if is_hybrid and not k_algo_specified: 
        with pytest.raises(RuntimeError) as e_info:
            E = psi4.energy(functional, molecule=molecule) 

        # we keep this line just for printout purposes; should always pass if done correctly 
        assert compare(type(e_info), pytest.ExceptionInfo, f'{scf_type}+{functional} throws RuntimeError')
    
    # ... else code will run fine
    else:
        E = psi4.energy(functional, molecule=molecule) 

@pytest.mark.parametrize("j_algo", [ "DFDIRJ" ]) #to be extended in the future
def test_j_algo_bp86(j_algo, mols):
    """Test SCF_TYPE={J} and all SCF_TYPE={J}+{K} combinations for a BP86 calculation.
    They should all give the exact same answer (within tolerance).""" 

    composite_K_algo = [ "LINK", "COSX" ]
  
    molecule = mols["h2o"]
    
    # run base composite J algorithm 
    psi4.set_options({"scf_type" : j_algo, "basis": "cc-pvdz"})
    energy_dfdirj = psi4.energy("bp86", molecule=molecule) 
    
    # compare composite combinations to base J algorithm
    for k_algo in composite_K_algo:
        scf_type = j_algo + "+" + k_algo    
        screening = "CSAM" if "COSX" in scf_type else "DENSITY"

        psi4.set_options({"scf_type" : scf_type, "reference": "rhf", "basis": "cc-pvdz", "screening": screening})
        energy_composite = psi4.energy("bp86", molecule=molecule) 
 
        assert compare_values(energy_dfdirj, energy_composite, 6, f'BP86 {scf_type} accurate to {j_algo} (1e-6 threshold)')

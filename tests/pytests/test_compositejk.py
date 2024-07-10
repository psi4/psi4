import pytest

from utils import compare_integers, compare_strings, compare_values, compare
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

@pytest.fixture
def tests():
    return {
        "h2o (rhf)": {
                      "method" : "hf",
                      "options": {"reference" : "rhf"},
                      "molecule" : "h2o",
                      "bsse_type" : None,
                      "ref" : -76.026780223322},
                      marks=pytest.mark.quick,
        "h2o (rks)": {
                      "method" : "b3lyp",
                      "options": {"reference" : "rhf"},
                      "molecule" : "h2o",
                      "bsse_type" : None,
                      "ref" : -76.420402720419},
                      marks=pytest.mark.quick,
        "nh2 (uhf)": {
                      "method" : "hf",
                      "options": {"reference" : "uhf"},
                      "molecule" : "nh2",
                      "bsse_type" : None,
                      "ref" : -55.566890252551},
        "nh2 (rohf)": {
                      "method" : "hf",
                      "options": {"reference" : "rohf"},
                      "molecule" : "nh2",
                      "bsse_type" : None,
                      "ref" : -55.562689948780},
        "h2o/na+ (rhf ie)": {
                      "method" : "hf",
                      "options": {"reference" : "rhf"},
                      "molecule" : "h2o_nap1",
                      "bsse_type" : "CP",
                      "ref" :  -0.040121884077},
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

@pytest.mark.parametrize("inp", 
    [
        pytest.param("h2o (rhf)", marks=pytest.mark.quick),
        pytest.param("h2o (rks)", marks=pytest.mark.quick),
    ], 
)
@pytest.mark.parametrize("opts", 
    [
        pytest.param({ "options": {
                           "scf_type": "DIRECT",
                           "screening": "density",
                           "maxiter": 30,
                       },
                       "ref_iter": {
                           "base": 8,
                           "df_scf_guess": 11,   
                           "scf_cosx_guess": 13, # assumes COSX_MAXITER_FINAL = -1   
                       },
                     }, id="DirectJK (MAXITER=30)"),
        pytest.param({ "options": {
                           "scf_type": "DIRECT",
                           "screening": "density",
                           "maxiter": 10,
                       },
                       "ref_iter": {
                           "base": 8,
                           "df_scf_guess": 11,   
                           "scf_cosx_guess": 13, # assumes COSX_MAXITER_FINAL = -1   
                       },
                     }, id="DirectJK (MAXITER=10)"),
        pytest.param({ "options": { 
                           "scf_type": "DFDIRJ+LINK",
                           "screening": "density",
                           "maxiter": 30,
                       },
                       "ref_iter": {
                           "base": 8,
                           "df_scf_guess": 11,   
                           "scf_cosx_guess": 13, # assumes COSX_MAXITER_FINAL = -1   
                       },
                     }, id="DF-DirJ+LinK"),
        pytest.param({ "options": { 
                           "scf_type": "DFDIRJ+COSX",
                           "screening": "schwarz",
                           "maxiter": 30,
                       },
                       "ref_iter": {
                           "base": 8, # for single initial COSX grid, i.e., COSX_MAXITER_FINAL = 0
                           "scf_cosx_guess": 20, # for full convergence on both COSX grids, i.e., COSX_MAXITER_FINAL = -1   
                       },
                     }, id="DF-DirJ+COSX (MAXITER=30)"),
        pytest.param({ "options": { 
                           "scf_type": "DFDIRJ+COSX",
                           "screening": "schwarz",
                           "maxiter": 10,
                       },
                       "ref_iter": {
                           "base": 8, # for single initial COSX grid, i.e., COSX_MAXITER_FINAL = 0
                           "scf_cosx_guess": 20, # for full convergence on both COSX grids, i.e., COSX_MAXITER_FINAL = -1   
                       },
                     }, id="DF-DirJ+COSX (MAXITER=10)"),
 
    ], 
)
@pytest.mark.parametrize("cosx_maxiter_final", [ -1, 0, 1 ])
@pytest.mark.parametrize("scf_cosx_guess", [ True, False ])
@pytest.mark.parametrize("df_scf_guess", [ True, False ])
def test_cosx_maxiter_final(inp, opts, cosx_maxiter_final, scf_cosx_guess, df_scf_guess, tests, mols, request):
    """Test the COSX_MAXITER_FINAL keyword in various situations via SCF calculations"""

    test_id = request.node.callspec.id
    
    # basic settings configuration
    molecule = mols[tests[inp]["molecule"]]
    psi4.set_options({
        "basis": "cc-pvdz",
        "cosx_maxiter_final": cosx_maxiter_final, 
        "df_scf_guess": df_scf_guess, 
        "scf_cosx_guess": scf_cosx_guess, 
        "save_jk": True
    })
    psi4.set_options(tests[inp]["options"])
    psi4.set_options(opts["options"])

    # DFDIRJ+COSX doesn't support either DF_SCF_GUESS or SCF_COSX_GUESS, so skip these
    if opts["options"]["scf_type"] == "DFDIRJ+COSX" and (scf_cosx_guess or df_scf_guess): 
        pytest.skip(f'{test_id} skipped: DFDIRJ+COSX doesnt support JK guesses') 

    # if both guesses are set to True, run should throw exception... 
    elif scf_cosx_guess and df_scf_guess:
        with pytest.raises(Exception) as e_info:
            E = psi4.energy(tests[inp]["method"], molecule=molecule)

        # we keep this line just for printout purposes; should always pass if done correctly
        assert compare(type(e_info), pytest.ExceptionInfo, f'{test_id} throws exception (both guesses are True)')
   
    # ... and enabling SCF_COSX_GUESS with COSX_MAXITER_FINAL set to 0 should also throw exception...
    elif scf_cosx_guess and cosx_maxiter_final == 0: 
        with pytest.raises(Exception) as e_info:
            E = psi4.energy(tests[inp]["method"], molecule=molecule)

        # we keep this line just for printout purposes; should always pass if done correctly
        assert compare(type(e_info), pytest.ExceptionInfo, f'{test_id} throws exception (invalid SCF_COSX_GUESS/COSX_MAXITER_FINAL combo)')

    # ... else proceed as normal
    else:
        # determine number of SCF iterations we expect based on settings
        guess_type = "base" # baseline number of iterations 
        if (opts["options"]["scf_type"] == "DFDIRJ+COSX" or scf_cosx_guess) and (cosx_maxiter_final == -1):
            guess_type = "scf_cosx_guess" # fully-converged SCF_COSX_GUESS (i.e., COSX_MAXITER_FINAL set to -1)
        elif df_scf_guess:
            guess_type = "df_scf_guess" # DF_SCF_GUESS enabled
        reference_iter = opts["ref_iter"][guess_type]

        # COSX_MAXITER_FINAL != -1 will add some number of iterations to baseline guess 
        if (opts["options"]["scf_type"] == "DFDIRJ+COSX" or scf_cosx_guess) and cosx_maxiter_final != -1:
            reference_iter += cosx_maxiter_final 
       
        # if maxiter < expected number of SCF iterations, run should throw exception... 
        if opts["options"]["maxiter"] < reference_iter: 
            with pytest.raises(Exception) as e_info:
                E = psi4.energy(tests[inp]["method"], molecule=molecule)
    
            # we keep this line just for printout purposes; should always pass if done correctly
            assert compare(type(e_info), pytest.ExceptionInfo, f'{test_id} throws exception')
        
        # ... otherwise, test should run with correct number of SCF iterations and post-guess method
        else:
            E, wfn = psi4.energy(tests[inp]["method"], molecule=molecule, return_wfn=True)
   
            # correct number of SCF iterations?            
            niter = wfn.variable("SCF ITERATIONS")
            assert compare_integers(niter, reference_iter, '{test_id} has correct number of SCF iterations')

            # correct post-guess method?
            clean_jk_name = wfn.jk().name().replace("-", "") # replace DF-DirJ with DFDirJ
            clean_jk_name = clean_jk_name.replace("DirectJK", "Direct") # DirectJK should be Direct instead
            assert compare_strings(opts["options"]["scf_type"].lower(), clean_jk_name.lower(), '{test_id} has correct end method')

@pytest.mark.parametrize("inp", 
    [
        pytest.param("h2o (rhf)"),
        pytest.param("h2o (rks)"),
        pytest.param("nh2 (uhf)"),
        pytest.param("nh2 (rohf)", marks=pytest.mark.nbody),
    ]
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
def test_seminum(inp, scf, tests, mols, request):
    """Test the DF-DirJ + COSX/sn-LinK JK objects via SCF calculations"""

    test_id = request.node.callspec.id
    
    molecule = mols[tests[inp]["molecule"]]
    psi4.set_options({"scf_type" : scf["scf_type"], "basis": "cc-pvdz"})
    psi4.set_options(tests[inp]["options"])
    if "snlink_force_cartesian" in scf.keys():
        psi4.set_options({"snlink_force_cartesian": scf["snlink_force_cartesian"]})

        #SNLINK_FORCE_CARTESIAN doesnt work with symmetry currently
        molecule.reset_point_group("C1")

    # does the SCF energy match a pre-computed reference?
    energy_seminum = psi4.energy(tests[inp]["method"], molecule=molecule, bsse_type=tests[inp]["bsse_type"])
    assert compare_values(scf["ref"][test_id.split("-")[1]], energy_seminum, 6, f'{test_id} accurate to reference (1e-6 threshold)')
)
    
    # is the SCF energy reasonably close to a conventional SCF?
    psi4.set_options({"scf_type" : "pk"})
    energy_pk = psi4.energy(tests[inp]["method"], molecule=molecule, bsse_type=tests[inp]["bsse_type"])
    assert compare_values(energy_pk, energy_seminum, 4, f'{test_id} DFDIRJ+COSX accurate to PK (1e-4 threshold)')

@pytest.mark.parametrize("inp", 
    [
        pytest.param("h2o (rhf)"),
        pytest.param("h2o (rks)"),
        pytest.param("nh2 (uhf)"),
        pytest.param("nh2 (rohf)"),
        pytest.param("nh2 (rohf)", marks=pytest.mark.nbody),
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
def test_seminum_incfock(inp, scf, tests, mols, request):
    """Test the efficiency of IncFock in DF-DirJ + COSX/sn-LinK JK objects via SCF calculations"""

    test_id = request.node.callspec.id

    molecule = mols[tests[inp]["molecule"]]
    psi4.set_options({"scf_type" : scf["scf_type"], "basis": "cc-pvdz", "incfock": False})
    psi4.set_options(tests[inp]["options"])

    # compute energy+wfn without IncFock 
    energy_seminum_noinc, wfn_seminum_noinc = psi4.energy(tests[inp]["method"], molecule=molecule, bsse_type=tests[inp]["bsse_type"], return_wfn=True)
    
    #assert compare_values(inp["ref"], energy_dfjcosk, atol=1e-6)

    # compute energy+wfn with Incfock 
    psi4.set_options({"incfock" : True})
    energy_seminum_inc, wfn_seminum_inc = psi4.energy(tests[inp]["method"], molecule=molecule, bsse_type=tests[inp]["bsse_type"], return_wfn=True)

    # how do energies compare?
    assert compare_values(energy_seminum_noinc, energy_seminum_inc, 6, f'{test_id} IncFock accurate (1e-6 threshold)')
    
    # how do SCF iteration counts compare?
    niter_noinc = int(wfn_seminum_noinc.variable("SCF ITERATIONS")) if wfn_seminum_noinc.has_scalar_variable("SCF ITERATIONS") else 0
    niter_inc = int(wfn_seminum_inc.variable("SCF ITERATIONS")) if wfn_seminum_inc.has_scalar_variable("SCF ITERATIONS") else 0
    
    assert compare(True, abs(niter_inc - niter_noinc) <= 3, f'{test_id} IncFock efficient')

@pytest.mark.parametrize("inp", 
    [
        pytest.param("h2o (rhf)"),
        pytest.param("h2o (rks)"),
        pytest.param("nh2 (uhf)"),
        pytest.param("nh2 (rohf)"),
        pytest.param("nh2 (rohf)", marks=pytest.mark.nbody),
    ], 
)
@pytest.mark.parametrize("scf_type", 
    [
        pytest.param("DIRECT"), 
        pytest.param("DFDIRJ+LINK"), 
    ], 
)
def test_scf_cos_guess(inp, scf_type, tests, mols, request):
    """Test the accuracy of the SCF_COSX_GUESS keyword via SCF calculations"""

    test_id = request.node.callspec.id

    molecule = mols[tests[inp]["molecule"]]
    psi4.set_options({"scf_type": scf_type, 
                      "screening": "density", 
                      "basis": "cc-pvdz", 
                      "incfock": True, 
                      "df_scf_guess": False, 
                      "scf_cosx_guess": False
    })
    psi4.set_options(tests[inp]["options"])

    # compute energy+wfn without SCF_COSX_GUESS 
    energy_noguess, wfn_noguess = psi4.energy(tests[inp]["method"], molecule=molecule, bsse_type=tests[inp]["bsse_type"], return_wfn=True)

    # compute DFJCOSK energy+wfn with SCF_COSX_GUESS
    psi4.set_options({"scf_cosx_guess" : True})
    energy_guess, wfn_guess = psi4.energy(tests[inp]["method"], molecule=molecule, bsse_type=tests[inp]["bsse_type"], return_wfn=True)

    # how do energies compare?
    assert compare_values(energy_noguess, energy_guess, 6, f'{test_id} SCF_COSX_GUESS accurate (1e-6 threshold)')

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

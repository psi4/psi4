import pytest

from utils import compare_integers, compare_values, compare

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
        "h2o (rks)": {
                      "method" : "b3lyp",
                      "options": {"reference" : "rhf"},
                      "molecule" : "h2o",
                      "bsse_type" : None,
                      "ref" : -76.420402720419},
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

@pytest.mark.parametrize("inp", 
    [
        pytest.param("h2o (rhf)"),
        pytest.param("h2o (rks)"),
        pytest.param("nh2 (uhf)"),
        pytest.param("nh2 (rohf)"),
        pytest.param("nh2 (rohf)", marks=pytest.mark.nbody),
    ], 
)
def test_dfjcosk(inp, tests, mols, request):
    """Test the DFJCOSK JK object via SCF calculations"""

    test_id = request.node.callspec.id
    
    molecule = mols[tests[inp]["molecule"]]
    psi4.set_options({"scf_type" : "dfdirj+cosx", "basis": "cc-pvdz"})
    psi4.set_options(tests[inp]["options"])

    # does the DFJCOSK SCF energy match a pre-computed reference?
    energy_dfjcosk = psi4.energy(tests[inp]["method"], molecule=molecule, bsse_type=tests[inp]["bsse_type"])
    assert compare_values(tests[inp]["ref"], energy_dfjcosk, 6, f'{test_id} DFDIRJ+COSX accurate to reference (1e-6 threshold)')

    # is the DFJCOSK SCF energy reasonably close to a conventional SCF?
    psi4.set_options({"scf_type" : "pk"})
    energy_pk = psi4.energy(tests[inp]["method"], molecule=molecule, bsse_type=tests[inp]["bsse_type"])
    assert compare_values(energy_pk, energy_dfjcosk, 4, f'{test_id} DFDIRJ+COSX accurate to PK (1e-4 threshold)')

@pytest.mark.parametrize("inp", 
    [
        pytest.param("h2o (rhf)"),
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
        
        # ... otherwise, test should run with correct number of iterations
        else:
            E, wfn = psi4.energy(tests[inp]["method"], molecule=molecule, return_wfn=True)
            
            niter = wfn.variable("SCF ITERATIONS")
            assert compare_integers(niter, reference_iter, '{test_id} has correct number of SCF iterations')

@pytest.mark.parametrize("inp", 
    [
        pytest.param("h2o (rhf)"),
        pytest.param("h2o (rks)"),
        pytest.param("nh2 (uhf)"),
        pytest.param("nh2 (rohf)"),
        pytest.param("nh2 (rohf)", marks=pytest.mark.nbody),
    ], 
)
def test_dfjcosk_incfock(inp, tests, mols, request):
    """Test the efficiency of IncFock in DFJCOSK JK object via SCF calculations"""

    test_id = request.node.callspec.id

    molecule = mols[tests[inp]["molecule"]]
    psi4.set_options({"scf_type" : "dfdirj+cosx", "basis": "cc-pvdz", "incfock": False})
    psi4.set_options(tests[inp]["options"])

    # compute DFJCOSK energy+wfn without IncFock 
    energy_dfjcosk_noinc, wfn_dfjcosk_noinc = psi4.energy(tests[inp]["method"], molecule=molecule, bsse_type=tests[inp]["bsse_type"], return_wfn=True)
    #assert compare_values(inp["ref"], energy_dfjcosk, atol=1e-6)

    # compute DFJCOSK energy+wfn with Incfock 
    psi4.set_options({"incfock" : True})
    energy_dfjcosk_inc, wfn_dfjcosk_inc = psi4.energy(tests[inp]["method"], molecule=molecule, bsse_type=tests[inp]["bsse_type"], return_wfn=True)

    # how do energies compare?
    assert compare_values(energy_dfjcosk_noinc, energy_dfjcosk_inc, 6, f'{test_id} IncFock accurate (1e-6 threshold)')
    
    # how do SCF iteration counts compare?
    niter_noinc = int(wfn_dfjcosk_noinc.variable("SCF ITERATIONS")) if wfn_dfjcosk_noinc.has_scalar_variable("SCF ITERATIONS") else 0
    niter_inc = int(wfn_dfjcosk_inc.variable("SCF ITERATIONS")) if wfn_dfjcosk_inc.has_scalar_variable("SCF ITERATIONS") else 0
    
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
@pytest.mark.parametrize("scf_type", [ "DFDIRJ", "LINK", "COSX", "DFDIRJ+COSX", "DFDIRJ+LINK" ])
def test_dfdirj(functional, scf_type, mols):
    """Test the functionality of the SCF_TYPE keyword for CompositeJK methods under varying situations:
      - Using hybrid DFT functionals without specifying a K algorithm should cause a RuntimeError to be thrown.
      - Not specifying a J algorithm should cause a ValidationError to be thrown."""

    composite_algo_to_matrix = {
        "DFDIRJ": "J",
        "LINK" : "K",
        "COSX": "K"
    }

    molecule = mols["h2o"]
    screening = "CSAM" if "COSX" in scf_type else "DENSITY"
    
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
@pytest.mark.parametrize("df_basis_scf", [ "CC-PVDZ-JKFIT", "DEF2-UNIVERSAL-JFIT" ]) #to be extended in the future
def test_j_algo_bp86(j_algo, df_basis_scf, mols):
    """Test SCF_TYPE={J} and all SCF_TYPE={J}+{K} combinations for a BP86 calculation.
    They should all give the exact same answer (within tolerance).""" 

    composite_K_algo = [ "LINK", "COSX" ]
  
    molecule = mols["h2o"]
    
    # run base composite J algorithm 
    psi4.set_options({"scf_type" : j_algo, "basis": "cc-pvdz", "df_basis_scf": df_basis_scf})
    energy_dfdirj = psi4.energy("bp86", molecule=molecule) 
    
    # compare composite combinations to base J algorithm
    for k_algo in composite_K_algo:
        scf_type = j_algo + "+" + k_algo    
        screening = "CSAM" if "COSX" in scf_type else "DENSITY"

        psi4.set_options({"scf_type" : scf_type, "reference": "rhf", "basis": "cc-pvdz", "df_basis_scf": df_basis_scf, "screening": screening})
        energy_composite = psi4.energy("bp86", molecule=molecule) 
 
        assert compare_values(energy_dfdirj, energy_composite, 6, f'BP86/{df_basis_scf} {scf_type} accurate to {j_algo} (1e-6 threshold)')

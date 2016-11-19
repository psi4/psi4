import numpy as np
np.set_printoptions(precision=14)

npoints = 5
keys = ["V", "V_RHO_A", "V_RHO_B", "V_GAMMA_AA", "V_GAMMA_AB", "V_GAMMA_BB"]
funcs = [
["S_X", "XC_LDA_X"],
["B88_X", "XC_GGA_X_B88"],
["B86B_X", "XC_GGA_X_B86_MGC"], 
["PW86_X", "XC_GGA_X_PW86"],
["PBE_X", "XC_GGA_X_PBE"],
#["RPBE_X", "XC_GGA_X_RPBE"],
#["SOGGA_X", "XC_GGA_X_SOGGA"],
["PW91_X", "XC_GGA_X_PW91"],
["FT97B_X", "XC_GGA_X_FT97_B"],
["PW92_C", "XC_LDA_C_PW"],
#["B_C"
#["M_C"
["LYP_C", "XC_GGA_C_LYP"],
["PZ81_C", "XC_LDA_C_PZ"],
["P86_C", "XC_GGA_C_P86"],
["PW91_C", "XC_GGA_C_PW91"],
["PBE_C", "XC_GGA_C_PBE"],
["FT97_C", "XC_GGA_C_FT97"],
["VWN3_C", "XC_LDA_C_VWN_3"],
["VWN5_C", "XC_LDA_C_VWN"],
["PW92A_C", "XC_LDA_C_PW_MOD"]]

#keys = ["V", "V_RHO_A", "V_GAMMA_AA", "V_GAMMA_AB", "V_GAMMA_BB", "V_RHO_B"]

rho_a = psi4.core.Vector.from_array(np.linspace(0.5, 0.99, npoints))
rho_b = psi4.core.Vector.from_array(np.linspace(0.5, 0.99, npoints))
sigma = psi4.core.Vector.from_array(np.ones((npoints)) * 0.3)
zeros = psi4.core.Vector.from_array(np.zeros((npoints)))

def build_in():
    inp = {
        'RHO_A' : rho_a,
        'RHO_B' : rho_b,
        'GAMMA_AA' : sigma,
        'GAMMA_AB' : sigma,
        'GAMMA_BB' : sigma,
    }
    return inp

def build_out():
    ret = {}
    for k in keys:
        ret[k] = psi4.core.Vector(npoints)
    return ret


for psi_func_name, xc_func_name in funcs:

    print("Building functional %s/%s" % (psi_func_name, xc_func_name))
    psi_fun = core.Functional.build_base(psi_func_name)
    if "GGA" in xc_func_name:
        psi_fun.set_gga(True)
    print_out("Psi4 functional\n")
    psi_fun.print_out()
    
    xc_fun = core.Functional.build_base(xc_func_name)
    print_out("XC functional\n")
    xc_fun.print_out()
    
    psi_out = build_out()
    psi_inp = build_in()
    psi_fun.compute_functional(psi_inp, psi_out, npoints, 1, 1.0)
    #psi_out["V_RHO_A"].np[:] *= 2
#    print("Called psi fun\n")
    
    xc_out = build_out()
    xc_inp = build_in()
    #xc_inp["GAMMA_AA"].np[:] *= 3
#    print("Called XC fun\n")
    xc_fun.compute_functional(xc_inp, xc_out, npoints, 1, 1.0)
    #xc_out["V"].np[:] *= (rho_b.np[:] + rho_a.np[:])
#    print('Here')
    
    for k in keys:
       if not np.allclose(psi_out[k].np, xc_out[k].np, atol=1.e-5):
            print("    Allclose failed for key %s" % k)
            print psi_out[k].np
            print xc_out[k].np
            print np.linalg.norm(psi_out[k].np - xc_out[k].np)
            raise Exception("Test failed")


import numpy as np
import psi4
np.set_printoptions(precision=14)

npoints = 5
keys = ["V", "V_RHO_A", "V_RHO_B", "V_GAMMA_AA", "V_GAMMA_AB", "V_GAMMA_BB", "V_TAU_A", "V_TAU_B"]
# xc_func_name = "XC_LDA_X"
# xc_func_name = "XC_GGA_X_B88"
# xc_func_name = "XC_GGA_X_B86_MGC" 
# xc_func_name = "XC_GGA_X_PW86"
# xc_func_name = "XC_GGA_X_PBE"
# xc_func_name = "XC_GGA_X_PW91"
# xc_func_name = "XC_GGA_X_FT97_B"
# xc_func_name = "XC_LDA_C_PW"
# xc_func_name = "XC_GGA_C_LYP"
# xc_func_name = "XC_LDA_C_PZ"
# xc_func_name = "XC_GGA_C_P86"
# xc_func_name = "XC_GGA_C_PW91"
# xc_func_name = "XC_GGA_C_PBE"
# xc_func_name = "XC_GGA_C_FT97"
# xc_func_name = "XC_LDA_C_VWN_3"
# xc_func_name = "XC_LDA_C_VWN"
# xc_func_name = "XC_LDA_C_PW_MOD"
xc_func_name = "XC_MGGA_X_MBEEF"

#keys = ["V", "V_RHO_A", "V_GAMMA_AA", "V_GAMMA_AB", "V_GAMMA_BB", "V_RHO_B"]

rho_a = psi4.core.Vector.from_array(np.linspace(0.5, 0.99, npoints))
rho_b = psi4.core.Vector.from_array(np.linspace(0.5, 0.99, npoints))
sigma = psi4.core.Vector.from_array(np.ones((npoints)) * 0.3)
sigma2 = psi4.core.Vector.from_array(np.ones((npoints)) * 0.3)
tau = psi4.core.Vector.from_array(np.ones((npoints)) * 5.0)
zeros = psi4.core.Vector.from_array(np.zeros((npoints)))

def build_in():
    inp = {
        'RHO_A' : rho_a,
        'RHO_B' : rho_b,
        'GAMMA_AA' : sigma,
        'GAMMA_AB' : sigma2,
        'GAMMA_BB' : sigma,
        'TAU_A' : tau,
        'TAU_B' : tau,
    }
    return inp

def build_out():
    ret = {}
    for k in keys:
        ret[k] = psi4.core.Vector(npoints)
    return ret



unp_fun = psi4.core.LibXCFunctional(xc_func_name, False)
unp_fun.print_out()

p_fun = psi4.core.LibXCFunctional(xc_func_name, True)
p_fun.print_out()

unp_out = build_out()
unp_inp = build_in()
print("Called psi fun\n")
unp_fun.compute_functional(unp_inp, unp_out, npoints, 1, 1.0)
#unp_out["V_RHO_A"].np[:] *= 2
print("Called psi fun\n")

p_out = build_out()
p_inp = build_in()
#p_inp["RHO_A"].np[:] *= 2.0
#p_inp["GAMMA_AA"].np[:] *= 4.0
# p_inp["TAU_A"].np[:] *= 0.5
#p_inp["GAMMA_AA"].np[:] += 2.0 * p_inp["GAMMA_AB"].np[:]
#p_inp["GAMMA_AA"].np[:] *= 3.0
# print("Called XC fun\n")
p_fun.compute_functional(p_inp, p_out, npoints, 1, 1.0)
#p_out["V_GAMMA_AA"].np[:] *= 4.0
#p_out["V"].np[:] *= (rho_b.np[:] + rho_a.np[:])
print('Here')

for k in keys:
    if unp_out[k].np[0] == 0: continue
    print(k)
    print('unpolar' + str(unp_out[k].np))
    print('polar  ' + str(p_out[k].np))
    print(np.linalg.norm(unp_out[k].np - p_out[k].np))
    print(" ")


unp_gamma = 2.0 * unp_out["V_GAMMA_AA"].np[:] + unp_out["V_GAMMA_AB"].np[:]
p_gamma =  2.0 * p_out["V_GAMMA_AA"].np[:]
print('Gamma (full)')
print('unpolar' + str(unp_gamma))
print('polar  ' + str(p_gamma))
print (np.linalg.norm(unp_gamma - p_gamma))
print(' ')

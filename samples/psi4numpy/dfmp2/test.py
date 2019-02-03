#! A simple Psi 4 input script to compute MP2 from a RHF reference

import psi4
import time
import numpy as np

psi4.set_output_file("output.dat", False)

mol = psi4.geometry("""
  H      0.5288      0.1610      0.9359
  C      0.0000      0.0000      0.0000
  H      0.2051      0.8240     -0.6786
  H      0.3345     -0.9314     -0.4496
  H     -1.0685     -0.0537      0.1921
symmetry c1
""")

psi4.set_options({'basis': 'cc-pVDZ'})

# Compute RHF refernce
psi4.core.print_out('\nStarting RHF...\n')

t = time.time()
RHF_E, wfn = psi4.energy('SCF', return_wfn=True)

psi4.core.print_out('...RHF finished in %.3f seconds:   %16.10f\n' % (time.time() - t, RHF_E))

# Split eigenvectors and eigenvalues into o and v
eps_occ = np.asarray(wfn.epsilon_a_subset("AO", "ACTIVE_OCC"))
eps_vir = np.asarray(wfn.epsilon_a_subset("AO", "ACTIVE_VIR"))
ndocc = eps_occ.shape[0]
nvirt = eps_vir.shape[0]

psi4.core.print_out('\nBuilding DF ERI tensor Qov...\n')
t = time.time()

# Build aux basis
aux = psi4.core.BasisSet.build(wfn.molecule(), "DF_BASIS_MP2",
                         psi4.core.get_option("DFMP2", "DF_BASIS_MP2"),
                         "RIFIT", psi4.core.get_global_option('BASIS'))

# Build and transform a OV basis
print("Initializing DF_Helper object.")
dfobj = psi4.core.DFHelper(wfn.basisset(), aux)
dfobj.set_method("DIRECT")
dfobj.initialize()


print("Transforming the AO->MO integrals.\n")
dfobj.add_space("O", wfn.Ca_subset("AO", "OCC"))
dfobj.add_space("V", wfn.Ca_subset("AO", "VIR"))
dfobj.add_transformation("Qov", "O", "V")

dfobj.transform()
Qov = dfobj.get_tensor("Qov")

psi4.core.print_out('...Qov build in %.3f seconds with a shape of %s, %.3f GB.\n'
                % (time.time() - t, str(Qov.shape), np.prod(Qov.shape) * 8.e-9))

psi4.core.print_out('\nComputing MP2 energy...\n')
t = time.time()

# This part of the denominator is identical for all i,j pairs
vv_denom = - eps_vir.reshape(-1, 1) - eps_vir

MP2corr_OS = 0.0
MP2corr_SS = 0.0
for i in range(ndocc):
    eps_i = eps_occ[i]
    i_Qv = Qov.np[:, i, :].copy()
    for j in range(i, ndocc):

        eps_j = eps_occ[j]
        j_Qv = Qov.np[:, j, :]
        tmp = np.dot(i_Qv.T, j_Qv)

        # Diagonal elements
        if i == j:
            div = 1.0 / (eps_i + eps_j + vv_denom)
        # Off-diagonal elements
        else:
            div = 2.0 / (eps_i + eps_j + vv_denom)

        MP2corr_OS += np.einsum('ab,ab,ab->', tmp, tmp, div)
        MP2corr_SS += np.einsum('ab,ab,ab->', tmp - tmp.T, tmp, div)

psi4.core.print_out('...finished computing MP2 energy in %.3f seconds.\n' % (time.time() - t))

MP2corr_E = MP2corr_SS + MP2corr_OS
MP2_E = RHF_E + MP2corr_E

SCS_MP2corr_E = MP2corr_SS / 3 + MP2corr_OS * (6. / 5)
SCS_MP2_E = RHF_E + SCS_MP2corr_E

psi4.core.print_out('\nMP2 SS correlation energy:         %16.10f\n' % MP2corr_SS)
psi4.core.print_out('MP2 OS correlation energy:         %16.10f\n' % MP2corr_OS)

psi4.core.print_out('\nMP2 correlation energy:            %16.10f\n' % MP2corr_E)
psi4.core.print_out('MP2 total energy:                  %16.10f\n' % MP2_E)

psi4.core.print_out('\nSCS-MP2 correlation energy:        %16.10f\n' % MP2corr_SS)
psi4.core.print_out('SCS-MP2 total energy:              %16.10f\n\n' % SCS_MP2_E)

psi4.compare_values(-0.0312016949, MP2corr_SS, 6, 'MP2 SS correlation energy')
psi4.compare_values(-0.1327414875, MP2corr_OS, 6, 'MP2 OS correlation energy')
psi4.compare_values(-0.1639431825, MP2corr_E, 6, 'MP2 correlation energy')
psi4.compare_values(-40.3626203385, MP2_E, 6, 'MP2 total energy')
psi4.compare_values(-40.3683675060, SCS_MP2_E, 6, 'SCS-MP2 total energy')

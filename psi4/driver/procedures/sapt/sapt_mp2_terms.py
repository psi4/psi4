#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2016 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

import numpy as np
import time

from psi4 import core
from psi4.driver.p4util.exceptions import *

from .sapt_util import print_sapt_var

def df_fdds_dispersion(primary, auxiliary, cache, do_print=True):

    if do_print:
        core.print_out("\n  ==> E20 FDDS Coupled-Dispersion <== \n\n")

    df_matrix_keys = ["Cocc_A", "Cvir_A", "Cocc_B", "Cvir_B"]
    fdds_matrix_cache = {key : cache[key] for key in df_matrix_keys}

    df_vector_keys = ["eps_occ_A", "eps_vir_A", "eps_occ_B", "eps_vir_B"]
    fdds_vector_cache = {key : cache[key] for key in df_vector_keys}

    fdds_obj = core.FDDS_Dispersion(primary, auxiliary, fdds_matrix_cache, fdds_vector_cache)

    # Densities
    D = fdds_obj.project_densities([cache["D_A"], cache["D_B"]])
    X_A = fdds_obj.form_unc_amplitude("A", 1.0)
    print(X_A.np[:5, :5])
    print(np.linalg.norm(X_A))

    print("Done with D!")
    # Check me!
    zero = core.BasisSet.zero_ao_basis_set()

    mints = core.MintsHelper(primary)
    metric = mints.ao_eri(zero, auxiliary, zero, auxiliary)
    auxiliary.print_detail_out()
    metric.power(-1.0, 1.e-12)
    metric = np.squeeze(metric)

    Saux = mints.ao_overlap(auxiliary, auxiliary)

    print("Metric close %s" % np.allclose(metric, fdds_obj.metric()))
    print("Saux close   %s" % np.allclose(Saux, fdds_obj.aux_overlap()))

    Qpq = np.squeeze(mints.ao_eri(zero, auxiliary, primary, primary))
    RA = np.einsum('Qpq,pq->Q', Qpq, cache["D_A"])
    RB = np.einsum('Qpq,pq->Q', Qpq, cache["D_B"])

    SA = np.dot(metric, RA)
    SB = np.dot(metric, RB)

    tc_aux = np.array(mints.ao_3coverlap(auxiliary, auxiliary, auxiliary))

    print("Dens close %s" % np.allclose(np.einsum("PQS,S->PQ", tc_aux, SA), D[0]))

    # Amps
    Ebs = - cache["eps_occ_A"].np.reshape(-1, 1) + cache["eps_vir_A"].np
    tmp = (4.0 * Ebs / (Ebs **2 + 1.0))


    # print(Ebs)
    print(tmp ** 0.5)
    trans = np.einsum('Qpq,qi,pa->Qia', Qpq, cache["Cocc_A"], cache["Cvir_A"])
    print(np.squeeze(trans.T)[:, 10:])
    X_Anp = np.einsum('Pia,ia,Qia->PQ', trans, tmp, trans)
    print(X_Anp[:5, :5])
    print(np.linalg.norm(X_Anp - X_A.np))





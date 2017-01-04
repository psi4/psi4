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

import psi4
import numpy as np
from psi4 import core

from psi4 import extras
from psi4.driver import p4util
from psi4.driver import qcdb
from psi4.driver.p4util.exceptions import *
from psi4.driver.molutil import *
from psi4.driver.procedures import scf_helper

def prepare_sapt_cache(wfn):
    jk = wfn.jk()
    jk.C_clear()
    jk.C_left_add(wfn.Ca_subset("SO", "OCC"))
    jk.compute()

    wfn.set_array("J mat", jk.J()[0].clone())
    wfn.set_array("K mat", jk.K()[0].clone())

    mints = core.MintsHelper(wfn.basisset())
    pot = mints.ao_potential()
    wfn.set_array("V mat", pot)

class DFT_SAPT(object):
    def __init__(self, **kwargs):
        optstash = p4util.OptionsState(
            ['SCF', 'SCF_TYPE'],
            ['SCF', 'SAVE_JK'])
    
        core.set_local_option('SCF', 'SAVE_JK', True)
    
        # Alter default algorithm
        if not core.has_option_changed('SCF', 'SCF_TYPE'):
            core.set_local_option('SCF', 'SCF_TYPE', 'DF')
    
        # Get the molecule of interest
        ref_wfn = kwargs.get('ref_wfn', None)
        if ref_wfn is None:
            sapt_dimer = kwargs.pop('molecule', core.get_active_molecule())
        else:
            core.print_out('Warning! SAPT argument "ref_wfn" is only able to use molecule information.')
            sapt_dimer = ref_wfn.molecule()
        sapt_dimer.update_geometry()  # make sure since mol from wfn, kwarg, or P::e
    
        # Shifting to C1 so we need to copy the active molecule
        if sapt_dimer.schoenflies_symbol() != 'c1':
            core.print_out('  SAPT does not make use of molecular symmetry, further calculations in C1 point group.\n')
            sapt_dimer = sapt_dimer.clone()
            sapt_dimer.reset_point_group('c1')
            sapt_dimer.fix_orientation(True)
            sapt_dimer.fix_com(True)
            sapt_dimer.update_geometry()
    
        sapt_dimer.reset_point_group('c1')
        sapt_dimer.fix_orientation(True)
        sapt_dimer.fix_com(True)
        sapt_dimer.update_geometry()

        functional = kwargs.pop('func', 'HF')
        core.set_global_option("DFT_FUNCTIONAL", functional.upper())
        core.set_local_option('SCF', 'REFERENCE', 'RKS')
        functional = "RKS"
    
        # if (core.get_option('SCF', 'REFERENCE') != 'RHF') and (name.upper() != "SAPT0"):
        #     raise ValidationError('Only SAPT0 supports a reference different from \"reference rhf\".')
    
        nfrag = sapt_dimer.nfragments()
        if nfrag != 2:
            raise ValidationError('SAPT requires active molecule to have 2 fragments, not %s.' % (nfrag))
    
        sapt_basis = 'dimer'
        if 'sapt_basis' in kwargs:
            sapt_basis = kwargs.pop('sapt_basis')
        sapt_basis = sapt_basis.lower()
    
        if sapt_basis == 'dimer':
            monomerA = sapt_dimer.extract_subsets(1, 2)
            monomerA.set_name('monomerA')
            monomerB = sapt_dimer.extract_subsets(2, 1)
            monomerB.set_name('monomerB')
        elif sapt_basis == 'monomer':
            monomerA = sapt_dimer.extract_subsets(1)
            monomerA.set_name('monomerA')
            monomerB = sapt_dimer.extract_subsets(2)
            monomerB.set_name('monomerB')
    
        ri = core.get_option('SCF', 'SCF_TYPE')
        df_ints_io = core.get_option('SCF', 'DF_INTS_IO')
        # inquire if above at all applies to dfmp2
    
        core.IO.set_default_namespace('dimer')
        core.print_out('\n')
        p4util.banner('Dimer HF')
        core.print_out('\n')
    
        # Compute dimer wavefunction
        if (sapt_basis == 'dimer') and (ri == 'DF'):
            core.set_global_option('DF_INTS_IO', 'SAVE')
    
        data = {}
        dimer_wfn = scf_helper(functional, molecule=sapt_dimer, **kwargs)
        data["SCF DIMER"] = psi4.core.get_variable("CURRENT ENERGY")
    
        if (sapt_basis == 'dimer') and (ri == 'DF'):
            core.set_global_option('DF_INTS_IO', 'LOAD')
    
        # Compute Monomer A wavefunction
        if (sapt_basis == 'dimer') and (ri == 'DF'):
            core.IO.change_file_namespace(97, 'dimer', 'monomerA')
  
        mon_a_shift = kwargs.pop('monomer_a_shift', False) 
        mon_b_shift = kwargs.pop('monomer_b_shift', False) 

        if mon_a_shift:
            core.set_global_option("DFT_GRAC_SHIFT", mon_a_shift) 

        core.IO.set_default_namespace('monomerA')
        core.print_out('\n')
        p4util.banner('Monomer A HF')
        core.print_out('\n')
        wfn_A = scf_helper(functional, molecule=monomerA, **kwargs)
        data["SCF MONOMERA"] = psi4.core.get_variable("CURRENT ENERGY")

        core.set_global_option("DFT_GRAC_SHIFT", 0.0) 
    
        # Compute Monomer B wavefunction
        if (sapt_basis == 'dimer') and (ri == 'DF'):
            core.IO.change_file_namespace(97, 'monomerA', 'monomerB')

        if mon_b_shift:
            core.set_global_option("DFT_GRAC_SHIFT", mon_b_shift) 

        core.IO.set_default_namespace('monomerB')
        core.print_out('\n')
        p4util.banner('Monomer B HF')
        core.print_out('\n')
        wfn_B = scf_helper(functional, molecule=monomerB, **kwargs)
        data["SCF MONOMERB"] = psi4.core.get_variable("CURRENT ENERGY")

        core.set_global_option("DFT_GRAC_SHIFT", 0.0) 
    
        # Delta MP2
        core.set_global_option('DF_INTS_IO', df_ints_io)
    
        if core.get_option('SCF', 'REFERENCE') == 'RHF':
            core.IO.change_file_namespace(p4const.PSIF_SAPT_MONOMERA, 'monomerA', 'dimer')
            core.IO.change_file_namespace(p4const.PSIF_SAPT_MONOMERB, 'monomerB', 'dimer') 
    
    
        # Prepare caches
        prepare_sapt_cache(wfn_A)
        prepare_sapt_cache(wfn_B)

        self.dimer_wfn = dimer_wfn
        self.wfn_A = wfn_A
        self.wfn_B = wfn_B
    
        self.Cocc_A = wfn_A.Ca_subset("SO", "OCC")
        self.Cvir_A = wfn_A.Ca_subset("SO", "VIR")
    
        self.Cocc_B = wfn_B.Ca_subset("SO", "OCC")
        self.Cvir_B = wfn_B.Ca_subset("SO", "VIR")

        self.eps_occ_A = wfn_A.epsilon_a_subset("SO", "OCC")
        self.eps_vir_A = wfn_A.epsilon_a_subset("SO", "VIR")

        self.eps_occ_B = wfn_B.epsilon_a_subset("SO", "OCC")
        self.eps_vir_B = wfn_B.epsilon_a_subset("SO", "VIR")
    
        self.nuc_rep = sapt_dimer.nuclear_repulsion_energy() - monomerA.nuclear_repulsion_energy() - monomerB.nuclear_repulsion_energy()
   

    def elst(self): 
        # ELST
        Elst10  = 4.0 * self.wfn_B.Da().vector_dot(self.wfn_A.get_array("J mat"))
        Elst10 += 2.0 * self.wfn_A.Da().vector_dot(self.wfn_B.get_array("V mat"))
        Elst10 += 2.0 * self.wfn_B.Da().vector_dot(self.wfn_A.get_array("V mat"))
        Elst10 += self.nuc_rep
        print('Electostatics  = % 16.8f' % (Elst10 * 1000))
        return Elst10
   
    def exch(self): 
        Cocc_A = self.Cocc_A
        Cvir_A = self.Cvir_A
        Cocc_B = self.Cocc_B
        Cvir_B = self.Cvir_B

        wfn_A = self.wfn_A
        wfn_B = self.wfn_B
        dimer_wfn = self.dimer_wfn

        SAB = core.Matrix.triplet(Cocc_A, wfn_A.S(), Cocc_B, True, False, False)
        num_occ = wfn_A.nalpha() + wfn_B.nalpha()
    
        h_A = wfn_A.get_array("V mat").clone()
        h_A.axpy(2.0, wfn_A.get_array("J mat"))
        h_A.axpy(-1.0, wfn_A.get_array("K mat"))
    
        h_B = wfn_B.get_array("V mat").clone()
        h_B.axpy(2.0, wfn_B.get_array("J mat"))
        h_B.axpy(-1.0, wfn_B.get_array("K mat"))
    
        Sab = core.Matrix(num_occ, num_occ)
        Sab.np[:wfn_A.nalpha(), wfn_A.nalpha():] = SAB.np
        Sab.np[wfn_A.nalpha():, :wfn_A.nalpha()] = SAB.np.T
        Sab.np[np.diag_indices_from(Sab.np)] += 1
        Sab.power(-1.0, 1.e-12)
        Sab.np[np.diag_indices_from(Sab.np)] -= 1.0

        Tmo_AA = core.Matrix.from_array(Sab.np[:wfn_A.nalpha(), :wfn_A.nalpha()])
        Tmo_BB = core.Matrix.from_array(Sab.np[wfn_A.nalpha():, wfn_A.nalpha():])
        Tmo_AB = core.Matrix.from_array(Sab.np[:wfn_A.nalpha(), wfn_A.nalpha():])
    
        T_A  = np.dot(Cocc_A, Tmo_AA).dot(Cocc_A.np.T)
        T_B  = np.dot(Cocc_B, Tmo_BB).dot(Cocc_B.np.T)
        T_AB = np.dot(Cocc_A, Tmo_AB).dot(Cocc_B.np.T)
    
        D_A = wfn_A.Da().np
        D_B = wfn_B.Da().np
        S = wfn_A.S().np
    
        P_A = np.dot(wfn_A.Ca(), wfn_A.Ca().np.T) - D_A
        P_B = np.dot(wfn_B.Ca(), wfn_B.Ca().np.T) - D_B
    
        w_A = wfn_A.get_array("V mat").clone()
        w_A.axpy(2.0, wfn_A.get_array("J mat"))
    
        w_B = wfn_B.get_array("V mat").clone()
        w_B.axpy(2.0, wfn_B.get_array("J mat"))
    
    
        jk = dimer_wfn.jk()   
        jk.C_clear()
    
        jk.C_left_add(Cocc_A)
        jk.C_right_add(core.Matrix.doublet(Cocc_A, Tmo_AA, False, False))
    
        jk.C_left_add(Cocc_B)
        jk.C_right_add(core.Matrix.doublet(Cocc_A, Tmo_AB, False, False))
    
        jk.C_left_add(psi4.core.Matrix.from_array( np.dot(P_B, wfn_B.S()).dot(Cocc_A)))
        jk.C_right_add(Cocc_A)
    
        jk.compute()
    
        JT_A, JT_AB, Jij = jk.J()
        KT_A, KT_AB, Kij = jk.K()
    
        Exch_s2 = 0.0
        Exch_s2 -= 2.0 * np.vdot(D_A.dot(S).dot(D_B).dot(S).dot(P_A), w_B)
        Exch_s2 -= 2.0 * np.vdot(D_B.dot(S).dot(D_A).dot(S).dot(P_B), w_A)
        Exch_s2 -= 2.0 * np.vdot(P_A.dot(S).dot(D_B), Kij.np.T)
        print('Exchange (S^2) = % 16.8f' % (Exch_s2 * 1000))
    
    
        # Start Sinf
        Exch10  = 0.0
        Exch10 -= 2.0 * np.vdot(wfn_A.Da(), wfn_B.get_array("K mat"))
        Exch10 += 2.0 * np.vdot(T_A, h_B.np)
        Exch10 += 2.0 * np.vdot(T_B, h_A.np)
        Exch10 += 2.0 * np.vdot(T_AB, h_A.np + h_B.np)
        Exch10 += 4.0 * np.vdot(T_B, JT_AB.np - 0.5 * KT_AB.np)
        Exch10 += 4.0 * np.vdot(T_A, JT_AB.np - 0.5 * KT_AB.np)
        Exch10 += 4.0 * np.vdot(T_B, JT_A.np - 0.5 *  KT_A.np)
        Exch10 += 4.0 * np.vdot(T_AB, JT_AB.np - 0.5 * KT_AB.np.T)
        print('Exchange       = % 16.8f' % (Exch10 * 1000))

    def ind(self):
        Cocc_A = self.Cocc_A
        Cvir_A = self.Cvir_A
        Cocc_B = self.Cocc_B
        Cvir_B = self.Cvir_B

        wfn_A = self.wfn_A
        wfn_B = self.wfn_B

        D_A = wfn_A.Da().np
        D_B = wfn_B.Da().np
        S = wfn_A.S().np


        dimer_wfn = self.dimer_wfn
        w_A = wfn_A.get_array("V mat").clone()
        w_A.axpy(2.0, wfn_A.get_array("J mat"))
    
        w_B = wfn_B.get_array("V mat").clone()
        w_B.axpy(2.0, wfn_B.get_array("J mat"))
    
        w_B_MOA = core.Matrix.triplet(Cocc_A, w_B, Cvir_A, True, False, False)
        w_A_MOB = core.Matrix.triplet(Cocc_B, w_A, Cvir_B, True, False, False)
    
        x_B_MOA = wfn_A.cphf_solve([w_B_MOA], 1.e-8, 20, 2)
        x_A_MOB = wfn_B.cphf_solve([w_A_MOB], 1.e-8, 20, 2)
    
        ind_ab = 2.0 * x_B_MOA[0].vector_dot(w_B_MOA)
        ind_ba = 2.0 * x_A_MOB[0].vector_dot(w_A_MOB)
        print('Induction B->A = % 16.8f' % (ind_ab * 1000))
        print('Induction A->B = % 16.8f' % (ind_ba * 1000))
    
    
        x_B = np.dot(Cocc_A, x_B_MOA[0]).dot(Cvir_A.np.T)
        x_A = np.dot(Cocc_B, x_A_MOB[0]).dot(Cvir_B.np.T)
    
        jk = dimer_wfn.jk()   
        jk.C_clear()
        
        C_O_A = psi4.core.Matrix.from_array(D_B.dot(S).dot(Cocc_A))
        C_P_A = psi4.core.Matrix.from_array(D_B.dot(S).dot(D_A).dot(S).dot(Cocc_B))
        C_P_B = psi4.core.Matrix.from_array(D_A.dot(S).dot(D_B).dot(S).dot(Cocc_A))

        jk.C_left_add(C_O_A)
        jk.C_right_add(Cocc_A)

        jk.C_left_add(C_P_A)
        jk.C_right_add(Cocc_B)

        jk.C_left_add(C_P_B)
        jk.C_right_add(Cocc_A)
    
        jk.compute()
    
        J_O, J_P_B, J_P_A = jk.J()
        K_O, K_P_B, K_P_A = jk.K()
        D_A = wfn_A.Da().np
        D_B = wfn_B.Da().np

        V_A = self.wfn_A.get_array("V mat").np
        J_A = self.wfn_A.get_array("J mat").np
        K_A = self.wfn_A.get_array("K mat").np

        V_B = self.wfn_B.get_array("V mat").np
        J_B = self.wfn_B.get_array("J mat").np
        K_B = self.wfn_B.get_array("K mat").np

        K_O = psi4.core.Matrix.from_array(K_O.np.T)
        W  = -1.0 * K_B.copy() 
        W -= 2.0 * np.dot(S, D_B).dot(J_A) 
        W += 1.0 * K_O.np
        W -= 2.0 * J_O.np

        W += 1.0 * np.dot(S, D_B).dot(K_A)
        W -= 2.0 * np.dot(J_B, D_B).dot(S)
        W += 1.0 * np.dot(K_B, D_B).dot(S)

        W += 2.0 * np.dot(S, D_B).dot(J_A).dot(D_B).dot(S)
        W += 2.0 * np.dot(J_B, D_A).dot(S).dot(D_B).dot(S)
        W -= 1.0 * np.dot(K_O, D_B).dot(S) 
        W += 2.0 * J_P_B.np

        W += 2.0 * np.dot(S, D_B).dot(S).dot(D_A).dot(J_B)
        W -= 1.0 * np.dot(S, D_B).dot(K_O.np.T)
        W -= 1.0 * np.dot(S, D_B).dot(V_A)
        W -= 1.0 * np.dot(V_B, D_B).dot(S)
        W += 1.0 * np.dot(S, D_B).dot(V_A).dot(D_B).dot(S)
        W += 1.0 * np.dot(V_B, D_A).dot(S).dot(D_B).dot(S)
        W += 1.0 * np.dot(S, D_B).dot(S).dot(D_A).dot(V_B)
        W = np.dot(Cocc_A.np.T, W).dot(Cvir_A)

        indexch_ab = 2.0 * np.vdot(W, x_B_MOA[0])
        print('IndExch   B->A = % 16.8f' % (indexch_ab * 1000))

        K_O = psi4.core.Matrix.from_array(K_O.np.T)
        W  = -1.0 * K_A.copy() 
        W -= 2.0 * np.dot(S, D_A).dot(J_B) 
        W += 1.0 * K_O.np
        W -= 2.0 * J_O.np
        W += 1.0 * np.dot(S, D_A).dot(K_B)
        W -= 2.0 * np.dot(J_A, D_A).dot(S)
        W += 1.0 * np.dot(K_A, D_A).dot(S)
        W += 2.0 * np.dot(S, D_A).dot(J_B).dot(D_A).dot(S)
        W += 2.0 * np.dot(J_A, D_B).dot(S).dot(D_A).dot(S)
        W -= 1.0 * np.dot(K_O, D_A).dot(S) 
        W += 2.0 * J_P_A.np
        W += 2.0 * np.dot(S, D_A).dot(S).dot(D_B).dot(J_A)
        W -= 1.0 * np.dot(S, D_A).dot(K_O.np.T)
        W -= 1.0 * np.dot(S, D_A).dot(V_B)
        W -= 1.0 * np.dot(V_A, D_A).dot(S)
        W += 1.0 * np.dot(S, D_A).dot(V_B).dot(D_A).dot(S)
        W += 1.0 * np.dot(V_A, D_B).dot(S).dot(D_A).dot(S)
        W += 1.0 * np.dot(S, D_A).dot(S).dot(D_B).dot(V_A)
        W = np.dot(Cocc_B.np.T, W).dot(Cvir_B)

        indexch_ba = 2.0 * np.vdot(W, x_A_MOB[0])
        print('IndExch   A->B = % 16.8f' % (indexch_ba * 1000))
    
    
    
        #return data

def sapt_printer(line, value):
    spacer = ' ' * (20 - len(line))
    print(line + spacer + '% 16.8f mH  % 16.8f kcal/mol  % 16.8f kelvin' % (value * 1000, value * 627.509, value * 315773))

#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2025 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

from psi4 import core

from ...constants import constants


def print_sapt_var(name, value, short=False, start_spacer="    "):
    """
    Converts the incoming value as hartree to a correctly formatted Psi print format.
    """

    vals = (name, value * 1000, value * constants.hartree2kcalmol, value * constants.hartree2kJmol)
    if short:
        return start_spacer + "%-28s % 15.8f [mEh]" % vals[:2]
    else:
        return start_spacer + "%-28s % 15.8f [mEh] % 15.8f [kcal/mol] % 15.8f [kJ/mol]" % vals


def print_sapt_hf_summary(data, name, short=False, delta_hf=False):

    ret = "   Partial %s Results, to compute Delta HF (dHF)\n" % name
    ret += "  " + "-" * 105 + "\n"

    # Elst
    ret += print_sapt_var("Electrostatics", data["Elst10,r"]) + "\n"
    ret += print_sapt_var("  Elst10,r", data["Elst10,r"]) + "\n"
    ret += "\n"
    core.set_variable("SAPT ELST ENERGY", data["Elst10,r"])

    # Exchange
    ret += print_sapt_var("Exchange", data["Exch10"]) + "\n"
    ret += print_sapt_var("  Exch10", data["Exch10"]) + "\n"
    ret += print_sapt_var("  Exch10(S^2)", data["Exch10(S^2)"]) + "\n"
    ret += "\n"
    core.set_variable("SAPT EXCH ENERGY", data["Exch10"])

    # Induction (no dHF)
    ind = data["Ind20,r"] + data["Exch-Ind20,r"]
    ind_ab = data["Ind20,r (A<-B)"] + data["Exch-Ind20,r (A<-B)"]
    ind_ba = data["Ind20,r (A->B)"] + data["Exch-Ind20,r (A->B)"]

    ret += print_sapt_var("Induction (no dHF)", ind) + "\n"
    ret += print_sapt_var("  Ind20,r", data["Ind20,r"]) + "\n"
    ret += print_sapt_var("  Exch-Ind20,r", data["Exch-Ind20,r"]) + "\n"
    ret += print_sapt_var("  Induction (A<-B) (no dHF)", ind_ab) + "\n"
    ret += print_sapt_var("  Induction (A->B) (no dHF)", ind_ba) + "\n"
    ret += "\n"
    core.set_variable("SAPT IND ENERGY", ind)

    if delta_hf:
        total_sapt = (data["Elst10,r"] + data["Exch10"] + ind)
        sapt_hf_delta = delta_hf - total_sapt

        core.set_variable("SAPT(DFT) Delta HF", sapt_hf_delta)

        ret += print_sapt_var("Subtotal SAPT(HF)", total_sapt, start_spacer="    ") + "\n"
        ret += print_sapt_var("Total HF", delta_hf, start_spacer="    ") + "\n"
        ret += print_sapt_var("Delta HF", sapt_hf_delta, start_spacer="    ") + "\n"

        ret += "  " + "-" * 105 + "\n"

        # Induction with dHF and scaling
        Sdelta = (ind + sapt_hf_delta) / ind

        ret += "   Partial %s Results, alternate SAPT0-like display\n" % name
        ret += "  " + "-" * 105 + "\n"

        ret += print_sapt_var("Induction", ind + sapt_hf_delta) + "\n"
        ret += print_sapt_var("  Ind20,r", data["Ind20,r"]) + "\n"
        ret += print_sapt_var("  Exch-Ind20,r", data["Exch-Ind20,r"]) + "\n"
        ret += print_sapt_var("  delta HF,r (2)", sapt_hf_delta) + "\n"
        ret += print_sapt_var("  Induction (A<-B)", ind_ab * Sdelta) + "\n"
        ret += print_sapt_var("  Induction (A->B)", ind_ba * Sdelta) + "\n"

        ret += "  " + "-" * 105 + "\n"
        return ret

    else:
        # Dispersion
        disp = data["Disp20"] + data["Exch-Disp20,u"]
        ret += print_sapt_var("Dispersion", disp) + "\n"
        ret += print_sapt_var("  Disp20", data["Disp20,u"]) + "\n"
        ret += print_sapt_var("  Exch-Disp20", data["Exch-Disp20,u"]) + "\n"
        ret += "\n"
        core.set_variable("SAPT DISP ENERGY", disp)

        # Total energy
        total = data["Elst10,r"] + data["Exch10"] + ind + disp
        ret += print_sapt_var("Total %-15s" % name, total, start_spacer="   ") + "\n"
        core.set_variable("SAPT0 TOTAL ENERGY", total)
        core.set_variable("SAPT TOTAL ENERGY", total)
        core.set_variable("CURRENT ENERGY", total)

        ret += "  " + "-" * 105 + "\n"
        return ret


def print_sapt_dft_summary(data, name, dimer_wfn, do_dft=True, short=False):
    ret = "   %s Results\n" % name
    ret += "  " + "-" * 105 + "\n"

    # Elst
    ret += print_sapt_var("Electrostatics", data["Elst10,r"]) + "\n"
    ret += print_sapt_var("  Elst1,r", data["Elst10,r"]) + "\n"
    ret += "\n"
    core.set_variable("SAPT ELST ENERGY", data["Elst10,r"])
    dimer_wfn.set_variable("SAPT ELST ENERGY", data["Elst10,r"])

    # Exchange
    ret += print_sapt_var("Exchange", data["Exch10"]) + "\n"
    ret += print_sapt_var("  Exch1", data["Exch10"]) + "\n"
    ret += print_sapt_var("  Exch1(S^2)", data["Exch10(S^2)"]) + "\n"
    ret += "\n"
    core.set_variable("SAPT EXCH ENERGY", data["Exch10"])
    dimer_wfn.set_variable("SAPT EXCH ENERGY", data["Exch10"])

    # Induction
    ind = data["Ind20,r"] + data["Exch-Ind20,r"]
    ind_ab = data["Ind20,r (A<-B)"] + data["Exch-Ind20,r (A<-B)"]
    ind_ba = data["Ind20,r (A->B)"] + data["Exch-Ind20,r (A->B)"]

    if "Delta HF Correction" in list(data):
        ind += data["Delta HF Correction"]

    ret += print_sapt_var("Induction", ind) + "\n"

    ret += print_sapt_var("  Ind2,r", data["Ind20,r"]) + "\n"
    ret += print_sapt_var("  Exch-Ind2,r", data["Exch-Ind20,r"]) + "\n"
    ret += print_sapt_var("  Induction (A<-B)", ind_ab) + "\n"
    ret += print_sapt_var("  Induction (A->B)", ind_ba) + "\n"

    if "Delta HF Correction" in list(data):
        ret += print_sapt_var("  delta HF,r (2)", data["Delta HF Correction"]) + "\n"

    ret += "\n"
    core.set_variable("SAPT IND ENERGY", ind)

    # Dispersion
    if do_dft:
        disp = data["Disp20"] + data["Exch-Disp20,r"]
        ret += print_sapt_var("Dispersion", disp) + "\n"
        ret += print_sapt_var("  Disp2,r", data["Disp20"]) + "\n"
        ret += print_sapt_var("  Disp2,u", data["Disp20,u"]) + "\n"
        if core.get_option("SAPT", "SAPT_DFT_EXCH_DISP_SCALE_SCHEME") != "NONE":
            ret += print_sapt_var("  Est. Exch-Disp2,r", data["Exch-Disp20,r"]) + "\n"
        ret += print_sapt_var("  Exch-Disp2,u", data["Exch-Disp20,u"]) + "\n"
        ret += "\n"
        core.set_variable("SAPT DISP ENERGY", disp)
        dimer_wfn.set_variable("SAPT DISP ENERGY", disp)
    else:
        disp = data["Disp20,u"] + data["Exch-Disp20,u"]
        ret += print_sapt_var("Dispersion", disp) + "\n"
        ret += print_sapt_var("  Disp20", data["Disp20,u"]) + "\n"
        ret += print_sapt_var("  Exch-Disp20", data["Exch-Disp20,u"]) + "\n"
        ret += "\n"
        core.set_variable("SAPT DISP ENERGY", disp)
        dimer_wfn.set_variable("SAPT DISP ENERGY", disp)
    
    # Total energy
    total = data["Elst10,r"] + data["Exch10"] + ind + disp
    ret += print_sapt_var("Total %-17s" % name, total, start_spacer="    ") + "\n"
    core.set_variable("SAPT(DFT) TOTAL ENERGY", total)
    core.set_variable("SAPT TOTAL ENERGY", total)
    core.set_variable("CURRENT ENERGY", total)
    dimer_wfn.set_variable("SAPT(DFT) TOTAL ENERGY", total)
    dimer_wfn.set_variable("SAPT TOTAL ENERGY", total)
    dimer_wfn.set_variable("CURRENT ENERGY", total)

    ret += "  " + "-" * 105 + "\n"
    return ret

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

from psi4 import core
from psi4.driver import p4const

def print_sapt_var(name, value, short=False):
    """
    Converts the incoming value as hartree to a correctly formatted Psi print format.
    """

    vals = (name, value * 1000, value * p4const.psi_hartree2kcalmol, value * p4const.psi_hartree2kJmol)
    if short:
        return "    %-20s %15.8f [mEh]" % vals[:2]
    else:
        return "    %-20s %15.8f [mEh] %15.8f [kcal/mol] %15.8f [kJ/mol]" % vals

def print_sapt_summary(data, name, short=False):

    ret = "    %s Results\n" % name
    ret += " " + "-" * 100 + "\n"

    # Elst
    ret += print_sapt_var(" Electrostatics", data["Elst10,r"]) + "\n"
    ret += print_sapt_var("   Elst10,r", data["Elst10,r"]) + "\n"
    ret += "\n"

    # Exchange
    ret += print_sapt_var(" Exchange", data["Exch10"]) + "\n"
    ret += print_sapt_var("   Exch10", data["Exch10"]) + "\n"
    ret += print_sapt_var("   Exch10(S^2)", data["Exch10(S^2)"]) + "\n"
    ret += "\n"

    ind = data["Ind20,r"] + data["Ind-Exch20,r"]
    ind_ab = data["Ind20,r (A<-B)"] + data["Ind-Exch20,r (A<-B)"]
    ind_ba = data["Ind20,r (A->B)"] + data["Ind-Exch20,r (A->B)"]

    ret += print_sapt_var(" Induction", ind) + "\n"
    ret += print_sapt_var("   Ind20,r", data["Ind20,r"]) + "\n"
    ret += print_sapt_var("   Ind-Exch20,r", data["Ind-Exch20,r"]) + "\n"
    ret += print_sapt_var("   Induction (A<-B)", ind_ab) + "\n"
    ret += print_sapt_var("   Induction (A->B)", ind_ba) + "\n"
    ret += "\n"

    total = data["Elst10,r"] + data["Exch10"] + ind
    ret += print_sapt_var("Total %s" % name, ind_ba) + "\n"
    ret += "\n"

    ret += " " + "-" * 100 + "\n"
    return ret


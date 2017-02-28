#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2017 The Psi4 Developers.
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

"""
Module to provide lightweight definitions of emperical dispersion terms.
"""
from psi4 import core
from psi4.driver.qcdb import interface_dftd3 as dftd3
from psi4.driver.qcdb import interface_gcp as gcp


class EmpericalDispersion(object):
    def __init__(self, alias, dtype, **kwargs):
        self.alias = alias.upper() # Functional Alias

        # Figure out dispersion type
        if dtype[0] != "-":
            dtype = "-" + dtype


        tuple_params = kwargs.pop('tuple_params', None)

        if dtype.lower() in dftd3.dash_alias.keys():
            self.dtype = dftd3.dash_alias[dtype.lower()].upper()
        else:
            self.dtype = dtype.upper() # Dispersion type

        if dtype.replace('-', '') in dftd3.dashcoeff.keys():
            self.dash_params = dftd3.dash_server(alias, dtype.replace('-', ''))
        else:
            self.dash_params = {'s6': 1.0}

        # Build coefficients
        if self.dtype in ["-D2GR", "-D3ZERO", "-D3BJ", "-D3MZERO", "-D3MBJ"]:
            self.dtype = self.dtype.replace('-D2GR', '-D2')
            self.disp_type = 'gr'

            # Odd tuple syntax favored by psi
            if (tuple_params is not None):
                self.tuple_params = None
                self.dash_params['s6'] = tuple_params[0]

                if len(tuple_params) > 1:
                    if "D2" in self.dtype:
                        self.dash_params["alpha6"] = tuple_params[1]
                    elif ("ZERO" in self.dtype) or ("BJ" in self.dtype):
                        self.dash_params["s8"] = tuple_params[1]

                if len(tuple_params) > 2:
                    if "ZERO" in self.dtype:
                        self.dash_params["sr6"] = tuple_params[2]
                    elif "BJ" in self.dtype:
                        self.dash_params["a1"] = tuple_params[2]

                if len(tuple_params) > 3:
                    if "ZERO" in self.dtype:
                        self.dash_params["alpha6"] = tuple_params[3]
                    elif "BJ" in self.dtype:
                        self.dash_params["a2"] = tuple_params[3]

                if len(tuple_params) > 4:
                    raise Exception("Too many parameter in input tuple param.")

        else: # Only other case at the moment
            self.dtype = self.dtype.replace('-D2P4', '-D2')
            self.disp_type = 'p4'
            if tuple_params is not None:
                self.dash_params = {}
                for k, v in zip(['s6', 'p1', 'p2', 'p3'], tuple_params):
                    self.dash_params[k] = v


        # Build the C++ dispersion class
        if self.disp_type == 'p4':
            self.disp = core.Dispersion.build(self.dtype, **self.dash_params)
        else:
            self.disp = None


        # Set user input
        for k, v in kwargs.keys():
            if k in self.dash_params.keys():
                self.dash_params[k] = kwargs.pop(k)

        if len(kwargs):
            raise Exception("The following DFTD3 parameters were not understood for %s dispersion type: %s" %
                            (dtype, ', '.join(kwargs.keys())))


        if self.dtype == "-D1":
            self.description = "    Grimme's -D1 Dispersion Correction"
            self.citation = "    Grimme, S. (2004), J. Comp. Chem., 25: 1463-1473"
            self.bibtex = "Grimme:2004:1463"

        elif self.dtype == "-D2":
            self.description = "    Grimme's -D2 Dispersion Correction"
            self.citation = "    Grimme, S. (2006),  J. Comp. Chem., 27: 1787-1799"
            self.bibtex = "Grimme:2006:1787"

        elif self.dtype == "-CHG":
            self.description = "    Chai and Head-Gordon Dispersion Correction"
            self.citation = "    Chai, J.-D.; Head-Gordon, M. (2010), J. Chem. Phys., 132: 6615-6620"
            self.bibtex = "Chai:2010:6615"

        elif self.dtype == "-DAS2009":
            self.description = "    Podeszwa and Szalewicz Dispersion Correction"
            self.citation = "    Pernal, K.; Podeszwa, R.; Patkowski, K.; Szalewicz, K. (2009), Phys. Rev. Lett., 103: 263201"
            self.bibtex = "Pernal:2009:263201"

        elif self.dtype == "-DAS2010":
            self.description = "    Podeszwa and Szalewicz Dispersion Correction"
            self.citation = "    Podeszwa, R.; Pernal, K.; Patkowski, K.; Szalewicz, K. (2010), J. Phys. Chem. Lett., 1: 550"
            self.bibtex = "Podeszwa:2010:550"

        elif self.dtype == "-D2GR":
            self.description = "    Grimme's -D2 Dispersion Correction"
            self.citation = "    Grimme, S. (2006),  J. Comp. Chem., 27: 1787-1799"
            self.bibtex = "Grimme:2006:1787"

        elif self.dtype == "-D3ZERO":
            self.description = "    Grimme's -D3 (zero-damping) Dispersion Correction"
            self.citation = "    Grimme S.; Antony J.; Ehrlich S.; Krieg H. (2010), J. Chem. Phys., 132: 154104"
            self.bibtex = "Grimme:2010:154104"

        elif self.dtype == "-D3BJ":
            self.description = "    Grimme's -D3 (BJ-damping) Dispersion Correction"
            self.citation = "    Grimme S.; Ehrlich S.; Goerigk L. (2011), J. Comput. Chem., 32: 1456"
            self.bibtex = "Grimme:2011:1456"

        elif self.dtype == "-D3MZERO":
            self.description = "    Grimme's -D3 (zero-damping, short-range refitted) Dispersion Correction"
            self.citation  = "    Grimme S.; Antony J.; Ehrlich S.; Krieg H. (2010), J. Chem. Phys., 132: 154104\n"
            self.citation += "    Smith, D. G. A.; Burns, L. A.; Patkowski, K.; Sherrill, C. D. (2016), J. Phys. Chem. Lett.; 7: 2197"
            self.bibtex = "Grimme:2010:154104"

        elif self.dtype == "-D3MBJ":
            self.description = "    Grimme's -D3 (BJ-damping, short-range refitted) Dispersion Correction"
            self.citation  = "    Grimme S.; Ehrlich S.; Goerigk L. (2011), J. Comput. Chem., 32: 1456\n"
            self.citation += "    Smith, D. G. A.; Burns, L. A.; Patkowski, K.; Sherrill, C. D. (2016), J. Phys. Chem. Lett.; 7: 2197"
            self.bibtex = "Grimme:2011:1456"

        else:
            raise Exception("Emperical Dispersion type %s not understood." % self.dtype)

    def print_out(self, level=1):

        core.print_out("   => %s: Empirical Dispersion <=\n\n" % self.dtype)
        core.print_out(self.description + "\n")

        core.print_out(self.citation + "\n\n")

        core.print_out("        S6 = %14.6E\n" % self.dash_params["s6"]);
        if "s8" in self.dash_params.keys():
            core.print_out("        S8 = %14.6E\n" % self.dash_params["s8"]);

        for k, v in self.dash_params.items():
            if k in ["s6", "s8"]: continue

            core.print_out("    %6s = %14.6E\n" % (k.upper(), v));

        # Psi auto sets this with no change of deviation
        if (self.disp_type == 'p4') and (self.dtype == "-D2"):
            core.print_out("    %6s = %14.6E\n" % ("A6", 20.0))

        core.print_out("\n");

    def compute_energy(self, molecule):
        if self.disp_type == 'gr':
            if self.alias in ['HF3C', 'PBEH3C']:
                dashd_part = dftd3.run_dftd3(molecule, dashlvl=self.dtype.lower().replace('-', ''),
                                             dashparam=self.dash_params, verbose=False, dertype=0)
                gcp_part = gcp.run_gcp(molecule, self.alias.lower(), verbose=False, dertype=0)
                return dashd_part + gcp_part
            else:
                return dftd3.run_dftd3(molecule, dashlvl=self.dtype.lower().replace('-', ''),
                                       dashparam=self.dash_params, verbose=False, dertype=0)
        else:
            return self.disp.compute_energy(molecule)

    def compute_gradient(self, molecule):
        if self.disp_type == 'gr':
            if self.alias in ['HF3C', 'PBEH3C']:
                dashd_part = dftd3.run_dftd3(molecule, dashlvl=self.dtype.lower().replace('-', ''),
                                             dashparam=self.dash_params, verbose=False, dertype=1)
                gcp_part = gcp.run_gcp(molecule, self.alias.lower(), verbose=False, dertype=1)
                dashd_part.add(gcp_part)
                return dashd_part
            else:
                return dftd3.run_dftd3(molecule, dashlvl=self.dtype.lower().replace('-', ''),
                                       dashparam=self.dash_params, verbose=False, dertype=1)
        else:
            return self.disp.compute_gradient(molecule)

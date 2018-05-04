#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2018 The Psi4 Developers.
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
"""
Module to provide lightweight definitions of emperical dispersion terms.
"""
from psi4 import core
from psi4.driver.qcdb import interface_dftd3 as dftd3
from psi4.driver.qcdb import interface_gcp as gcp
from psi4.driver.qcdb.dashparam import get_dispersion_aliases
from psi4.driver.qcdb.dashparam import get_default_dashparams
from psi4.driver import p4util
#import numpy as np

class EmpericalDispersion(object):
    def __init__(self, alias, dtype, **kwargs):
        # 1) Functional name processing:
        # 1a) Cleave out base functional from alias:
        for dash in ["-" + name.lower() for name in get_dispersion_aliases()]:
            if dash == alias.lower()[-len(dash):]:
                alias = alias[:-len(dash)]

        # 1b) Alias must be lowercase
        self.alias = alias.lower()

        # 2) Figure out dispersion type:
        # 2a) Strip "-" from dtype
        if dtype[0] == "-":
            dtype = dtype[1:]

        # 2b) Un-alias and capitalise dtype for printing
        if dtype.lower() in get_dispersion_aliases():
            self.dtype = "-" + get_dispersion_aliases()[dtype.lower()]
        else:
            self.dtype = "-" + dtype.lower()

        # 3) Get dispersion parameters:
        # 3a) Set defaults
        self.dash_params = get_default_dashparams(dtype)

        # 3b) Load passed variables from dictionary or from functional type
        tuple_params = kwargs.pop('tuple_params', None)
        if "dashparams" in kwargs:
            self.dash_params.update(kwargs.pop("dashparams"))
        elif dtype in dftd3.dashcoeff:
            self.dash_params.update(dftd3.dash_server(alias, dtype))
        else:
            self.dash_params = {'s6': 1.0}

        # 4) Dispersion class build process:
        # 4a) Build coefficients for dftd3
        if self.dtype in ["-d2gr", "-d3zero", "-d3bj", "-d3mzero", "-d3mbj"]:
            self.dtype = self.dtype.replace('-d2gr', '-d2')
            self.disp_type = 'gr'

            # Odd tuple syntax favored by psi
            if (tuple_params is not None):
                self.tuple_params = None
                self.dash_params['s6'] = tuple_params[0]

                if len(tuple_params) > 1:
                    if "d2" in self.dtype:
                        self.dash_params["alpha6"] = tuple_params[1]
                    elif ("zero" in self.dtype) or ("bj" in self.dtype):
                        self.dash_params["s8"] = tuple_params[1]

                if len(tuple_params) > 2:
                    if "zero" in self.dtype:
                        self.dash_params["sr6"] = tuple_params[2]
                    elif "bj" in self.dtype:
                        self.dash_params["a1"] = tuple_params[2]

                if len(tuple_params) > 3:
                    if "zero" in self.dtype:
                        self.dash_params["alpha6"] = tuple_params[3]
                    elif "bj" in self.dtype:
                        self.dash_params["a2"] = tuple_params[3]

                if len(tuple_params) > 4:
                    raise Exception("Too many parameter in input tuple param.")

        # DFT-NL dispersion. 
        elif self.dtype in ["-nl"]:
             self.disp_type = 'nl'


        # 4b) Build coefficients for psi4
        else:
            self.dtype = self.dtype.replace('-d2p4', '-d2')
            self.disp_type = 'p4'
            if tuple_params is not None:
                self.dash_params = {}
                for k, v in zip(['s6', 'p1', 'p2', 'p3'], tuple_params):
                    self.dash_params[k] = v

        # 4c) Build the C++ dispersion class for psi4
        if self.disp_type == 'p4':
            self.disp = core.Dispersion.build(self.dtype, **self.dash_params)
        else:
            self.disp = None

        # 5) Override parameters from user input
        # 5a) pop citation if present
        if "citation" in kwargs:
            custom_citation = kwargs.pop("citation")
        else:
            custom_citation = False

        # 5b) process other kwargs
        for k, v in kwargs.keys():
            if k in self.dash_params.keys():
                self.dash_params[k] = kwargs.pop(k)

        if len(kwargs):
            raise Exception("The following parameters in empirical_dispersion.py were not understood for %s dispersion type: %s" %
                            (dtype, ', '.join(kwargs.keys())))

        # 6) Process citations
        # 6a) Set default citations for method
        if self.dtype == "-d1":
            self.description = "    Grimme's -D1 Dispersion Correction"
            self.citation = "    Grimme, S. (2004), J. Comp. Chem., 25: 1463-1473"
            self.bibtex = "Grimme:2004:1463"

        elif self.dtype == "-d2":
            self.description = "    Grimme's -D2 Dispersion Correction"
            self.citation = "    Grimme, S. (2006),  J. Comp. Chem., 27: 1787-1799"
            self.bibtex = "Grimme:2006:1787"

        elif self.dtype == "-chg":
            self.description = "    Chai and Head-Gordon Dispersion Correction"
            self.citation = "    Chai, J.-D.; Head-Gordon, M. (2010), J. Chem. Phys., 132: 6615-6620"
            self.bibtex = "Chai:2010:6615"

        elif self.dtype == "-das2009":
            self.description = "    Podeszwa and Szalewicz Dispersion Correction"
            self.citation = "    Pernal, K.; Podeszwa, R.; Patkowski, K.; Szalewicz, K. (2009), Phys. Rev. Lett., 103: 263201"
            self.bibtex = "Pernal:2009:263201"

        elif self.dtype == "-das2010":
            self.description = "    Podeszwa and Szalewicz Dispersion Correction"
            self.citation = "    Podeszwa, R.; Pernal, K.; Patkowski, K.; Szalewicz, K. (2010), J. Phys. Chem. Lett., 1: 550"
            self.bibtex = "Podeszwa:2010:550"

        elif self.dtype == "-d2gr":
            self.description = "    Grimme's -D2 Dispersion Correction"
            self.citation = "    Grimme, S. (2006),  J. Comp. Chem., 27: 1787-1799"
            self.bibtex = "Grimme:2006:1787"

        elif self.dtype == "-d3zero":
            self.description = "    Grimme's -D3 (zero-damping) Dispersion Correction"
            self.citation = "    Grimme S.; Antony J.; Ehrlich S.; Krieg H. (2010), J. Chem. Phys., 132: 154104"
            self.bibtex = "Grimme:2010:154104"

        elif self.dtype == "-d3bj":
            self.description = "    Grimme's -D3 (BJ-damping) Dispersion Correction"
            self.citation = "    Grimme S.; Ehrlich S.; Goerigk L. (2011), J. Comput. Chem., 32: 1456"
            self.bibtex = "Grimme:2011:1456"

        elif self.dtype == "-d3mzero":
            self.description = "    Grimme's -D3 (zero-damping, short-range refitted) Dispersion Correction"
            self.citation = "    Grimme S.; Antony J.; Ehrlich S.; Krieg H. (2010), J. Chem. Phys., 132: 154104\n"
            self.citation += "    Smith, D. G. A.; Burns, L. A.; Patkowski, K.; Sherrill, C. D. (2016), J. Phys. Chem. Lett.; 7: 2197"
            self.bibtex = "Grimme:2010:154104"

        elif self.dtype == "-d3mbj":
            self.description = "    Grimme's -D3 (BJ-damping, short-range refitted) Dispersion Correction"
            self.citation = "    Grimme S.; Ehrlich S.; Goerigk L. (2011), J. Comput. Chem., 32: 1456\n"
            self.citation += "    Smith, D. G. A.; Burns, L. A.; Patkowski, K.; Sherrill, C. D. (2016), J. Phys. Chem. Lett.; 7: 2197"
            self.bibtex = "Grimme:2011:1456"

        elif self.dtype == "-nl":
            self.description = "    Grimme's -NL (DFT plus  VV10 correlation) "
            self.citation = "    Hujo, W.; Grimme, S; (2011), J. Chem. Theory Comput.; 7:3866 \n"
            self.bibtex = "Grimme:2011:3866"

        else:
            raise Exception("Empirical Dispersion type %s not understood." % self.dtype)

        # 6b) add custom citations if available
        if custom_citation:
            self.citation += "\n    Parametrisation from: \n" + custom_citation

    def print_out(self, level=1):

        core.print_out("   => %s: Empirical Dispersion <=\n\n" % self.dtype.upper())
        core.print_out(self.description + "\n")

        core.print_out(self.citation + "\n\n")
       
        if self.disp_type=='nl':
           return

        core.print_out("        S6 = %14.6E\n" % self.dash_params["s6"])
        if "s8" in self.dash_params.keys():
            core.print_out("        S8 = %14.6E\n" % self.dash_params["s8"])

        for k, v in self.dash_params.items():
            if k in ["s6", "s8"]: continue

            core.print_out("    %6s = %14.6E\n" % (k.upper(), v))

        # Psi auto sets this with no change of deviation
        if (self.disp_type == 'p4') and (self.dtype == "-D2"):
            core.print_out("    %6s = %14.6E\n" % ("A6", 20.0))

        core.print_out("\n")

    def compute_energy(self, molecule):
        if self.disp_type == 'gr':
            if self.alias in ['hf3c', 'pbeh3c']:
                dashd_part = dftd3.run_dftd3(
                    molecule,
                    dashlvl=self.dtype.lower().replace('-', ''),
                    dashparam=self.dash_params,
                    verbose=False,
                    dertype=0)
                gcp_part = gcp.run_gcp(molecule, self.alias.lower(), verbose=False, dertype=0)
                return dashd_part + gcp_part
            else:
                return dftd3.run_dftd3(
                    molecule,
                    dashlvl=self.dtype.lower().replace('-', ''),
                    dashparam=self.dash_params,
                    verbose=False,
                    dertype=0)
        else:
            return self.disp.compute_energy(molecule)

    def compute_gradient(self, molecule):
        if self.disp_type == 'gr':
            if self.alias in ['hf3c', 'pbeh3c']:
                dashd_part = dftd3.run_dftd3(
                    molecule,
                    dashlvl=self.dtype.lower().replace('-', ''),
                    dashparam=self.dash_params,
                    verbose=False,
                    dertype=1)
                gcp_part = gcp.run_gcp(molecule, self.alias.lower(), verbose=False, dertype=1)
                dashd_part.add(gcp_part)
                return dashd_part
            else:
                return dftd3.run_dftd3(
                    molecule,
                    dashlvl=self.dtype.lower().replace('-', ''),
                    dashparam=self.dash_params,
                    verbose=False,
                    dertype=1)
        else:
            return self.disp.compute_gradient(molecule)

    def compute_hessian(self, molecule):
        """
        #magic (if magic was easy)
        """

        optstash = p4util.OptionsState(['PRINT'])
        core.set_global_option('PRINT', 0)

        core.print_out("\n\n   Analytical Dispersion Hessians are not supported by dftd3 or gcp.\n")
        core.print_out("       Computing the Hessian through finite difference of gradients.\n\n")

        # Setup the molecule
        molclone = molecule.clone()
        molclone.reinterpret_coordentry(False)
        molclone.fix_orientation(True)
        molclone.fix_com(True)

        # Record undisplaced symmetry for projection of diplaced point groups
        core.set_parent_symmetry(molecule.schoenflies_symbol())

        gradients = []
        for geom in core.fd_geoms_freq_1(molecule, -1):
            molclone.set_geometry(geom)
            molclone.update_geometry()
            gradients.append(self.compute_gradient(molclone))

        H = core.fd_freq_1(molecule, gradients, -1)
        # H.print_out()
        optstash.restore()
        return H

#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2024 The Psi4 Developers.
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

import collections
from typing import Dict, List, Union

import numpy as np
import qcengine as qcng
from qcelemental.models import AtomicInput

from psi4 import core

from .. import p4util
from ..p4util.exceptions import ValidationError

_engine_can_do = collections.OrderedDict([
    # engine order establishes default for each disp
    ("libdisp",  ["d1", "d2",                                                                                                               "chg", "das2009", "das2010",]),
    ("s-dftd3",  [            "d3zero2b", "d3bj2b", "d3mzero2b", "d3mbj2b", "d3zeroatm", "d3bjatm", "d3mzeroatm", "d3mbjatm",                                           ]),
    ("dftd3",    [      "d2", "d3zero2b", "d3bj2b", "d3mzero2b", "d3mbj2b",                                                                                             ]),
    ("nl",       [                                                                                                                          "nl",                       ]),
    ("mp2d",     [                                                                                                                          "dmp2",                     ]),
    ("dftd4",    [                                                                                                            "d4bjeeqatm",                             ]),
    ("mctc-gcp", [                                                                                                                          "3c",                       ]),
    ("gcp",      [                                                                                                                          "3c",                       ]),
]) # yapf: disable


def _capable_engines_for_disp()-> Dict[str, List[str]]:
    """Invert _engine_can_do dictionary and check program detection.

    Returns a dictionary with keys all dispersion levels and values a list of all
    capable engines, where the engine in the first element is available, if any are.

    """
    try:
        from qcengine.testing import _programs as _programs_qcng
    except ModuleNotFoundError:
        # _programs_qcng is up-to-date with current harnesses but it requires pytest present, so let's provide a workaround
        from qcelemental.util import which, which_import
        _programs_qcng = {
            "dftd3": which("dftd3", return_bool=True),
            "dftd4": which_import("dftd4", return_bool=True),
            "s-dftd3": which_import("dftd3", return_bool=True),
            "mctc-gcp": which("mctc-gcp", return_bool=True),
            "gcp": which("gcp", return_bool=True),
            "mp2d": which("mp2d", return_bool=True),
        }

    programs_disp = {k: v for k, v in _programs_qcng.items() if k in _engine_can_do}
    programs_disp["libdisp"] = True
    programs_disp["nl"] = True

    capable = collections.defaultdict(list)
    capable_sorted_by_available = collections.defaultdict(list)
    for eng, disps in _engine_can_do.items():
        for disp in disps:
            capable[disp].append(eng)
    for disp, engines in capable.items():
        capable_sorted_by_available[disp] = sorted(engines, key=lambda x: programs_disp[x], reverse=True)

    return capable_sorted_by_available


class EmpiricalDispersion():
    """Lightweight unification of empirical dispersion calculation modes.

    Attributes
    ----------
    dashlevel : str
        {"d1", "d2", "chg", "das2009", "das2010", "nl", "dmp2", "d3zero2b", "d3bj2b", "d3mzero2b", "d3mbj2b", "d3zeroatm", "d3bjatm", "d3mzeroatm", "d3mbjatm", "d4bjeeqatm"}
        Name of dispersion correction to be applied. Resolved
        from `name_hint` and/or `level_hint` into a key of
        `empirical_dispersion_resources.dashcoeff`.
    dashparams : dict
        Complete set of parameter values defining the flexible parts
        of :py:attr:`dashlevel`. Number and parameter names vary by
        :py:attr:`dashlevel`. Resolved into a complete set (keys of
        dashcoeff[dashlevel]['default']) from `name_hint` and/or
        `dashcoeff_supplement` and/or user `param_tweaks`.
    fctldash : str
        If :py:attr:`dashparams` for :py:attr:`dashlevel` corresponds to a defined,
        named, untweaked "functional-dashlevel" set, then that
        functional. Otherwise, empty string.
    description : str
        Tagline for dispersion :py:attr:`dashlevel`.
    dashlevel_citation : str
        Literature reference for dispersion :py:attr:`dashlevel` in general,
        *not necessarily* for :py:attr:`dashparams`.
    dashparams_citation : str
        Literature reference for dispersion parameters, if :py:attr:`dashparams`
        corresponds to a defined, named, untweaked "functional-dashlevel"
        set with a citation. Otherwise, empty string.
    dashcoeff_supplement : dict
        See description in `qcengine.programs.empirical_dispersion_resources.from_arrays`. Used
        here to "bless" the dispersion definitions attached to
        the procedures/dft/<rung>_functionals-defined dictionaries
        as legit, non-custom, and of equal validity to
        `qcengine.programs.empirical_dispersion_resources.dashcoeff` itself for purposes of
        validating :py:attr:`fctldash`.
    engine : str
        {'libdisp', "s-dftd3", 'dftd3', 'nl', 'mp2d', "dftd4"}
        Compute engine for dispersion. One of Psi4's internal libdisp
        library, external Grimme or Beran projects, or nl.
    disp : Dispersion
        Only present for :py:attr:`engine` `=libdisp`. Psi4 class instance prepared
        to compute dispersion.
    ordered_params : list
        Fixed-order list of relevant parameters for :py:attr:`dashlevel`. Matches
        :rst:psivar:`DISPERSION CORRECTION ENERGY` ordering. Used for printing.

    Parameters
    ----------
    name_hint
        Name of functional (func only, func & disp, or disp only) for
        which to compute dispersion (e.g., blyp, BLYP-D2, blyp-d3bj,
        blyp-d3(bj), hf+d). Any or all parameters initialized from
        ``dashcoeff[dashlevel][functional-without-dashlevel]`` or
        ``dashcoeff_supplement[dashlevel][functional-with-dashlevel]``
        can be overwritten via `param_tweaks`.
    level_hint
        Name of dispersion correction to be applied (e.g., d, D2,
        d3(bj), das2010). Must be key in `dashcoeff` or "alias" or
        "formal" to one.
    param_tweaks
        Values for the same keys as `dashcoeff[dashlevel]['default']`
        (and same order if list) used to override any or all values
        initialized by `name_hint`.  Extra parameters will error.
    engine
        Override which code computes dispersion. See above for allowed
        values. Formerly (pre Nov 2022) only relevant for -D2, which can be computed by
        libdisp or dftd3. Now (post Nov 2022) also relevant for -D3 variants,
        which can be computed by dftd3 executable or simple-dftd3 Python module.
    gcp_engine
        Override which code computes the gcp correction. Now can use
        classic gcp or mctc-gcp executables.
    save_pairwise_disp
        Whether to request atomic pairwise analysis.

    """
    def __init__(self, *, name_hint: str = None, level_hint: str = None, param_tweaks: Union[Dict, List] = None, engine: str = None, gcp_engine: str = None, save_pairwise_disp: bool = False):
        from .dft import dashcoeff_supplement
        self.dashcoeff_supplement = dashcoeff_supplement
        self.save_pairwise_disp = save_pairwise_disp

        resolved = qcng.programs.empirical_dispersion_resources.from_arrays(
            name_hint=name_hint,
            level_hint=level_hint,
            param_tweaks=param_tweaks,
            dashcoeff_supplement=self.dashcoeff_supplement)
        self.fctldash = resolved['fctldash']
        self.dashlevel = resolved['dashlevel']
        self.dashparams = resolved['dashparams']
        self.description = qcng.programs.empirical_dispersion_resources.dashcoeff[self.dashlevel]['description']
        self.ordered_params = qcng.programs.empirical_dispersion_resources.dashcoeff[self.dashlevel]['default'].keys()
        self.dashlevel_citation = qcng.programs.empirical_dispersion_resources.dashcoeff[self.dashlevel]['citation']
        self.dashparams_citation = resolved['dashparams_citation']

        capable_engines_for_disp = _capable_engines_for_disp()
        if engine is None:
            self.engine = capable_engines_for_disp[self.dashlevel][0]
        else:
            if self.dashlevel in _engine_can_do[engine]:
                self.engine = engine
            else:
                raise ValidationError(f"This little engine ({engine}) can't ({self.dashlevel})")

        if self.engine == 'libdisp':
            self.disp = core.Dispersion.build(self.dashlevel, **resolved['dashparams'])

        if gcp_engine is None:
            self.gcp_engine = capable_engines_for_disp["3c"][0]
        else:
            if "3c" in _engine_can_do[gcp_engine]:
                self.gcp_engine = gcp_engine
            else:
                raise ValidationError(f"This little engine ({engine}) can't (3c)")

    def print_out(self):
        """Format dispersion parameters of `self` for output file."""

        text = []
        text.append("   => {}: Empirical Dispersion <=".format(
            (self.fctldash.upper() if self.fctldash.upper() else 'Custom')))
        text.append('')
        text.append(self.description)
        text.append(self.dashlevel_citation.rstrip())
        if self.dashparams_citation:
            text.append("    Parametrisation from:{}".format(self.dashparams_citation.rstrip()))
        text.append('')
        for op in self.ordered_params:
            text.append("    %6s = %14.6f" % (op, self.dashparams[op]))
        text.append('\n')

        core.print_out('\n'.join(text))

    def compute_energy(self, molecule: core.Molecule, wfn: core.Wavefunction = None) -> float:
        """Compute dispersion energy based on engine, dispersion level, and parameters in `self`.

        Parameters
        ----------
        molecule
            System for which to compute empirical dispersion correction.
        wfn
            Location to set QCVariables

        Returns
        -------
        float
            Dispersion energy [Eh].

        Notes
        -----
        :psivar:`DISPERSION CORRECTION ENERGY`
            Disp always set. Overridden in SCF finalization, but that only changes for "-3C" methods.
        :psivar:`fctl DISPERSION CORRECTION ENERGY`
            Set if :py:attr:`fctldash` nonempty.

        """
        if self.engine in ["s-dftd3", 'dftd3', 'mp2d', "dftd4"]:
            resi = AtomicInput(
                **{
                    'driver': 'energy',
                    'model': {
                        'method': self.fctldash,
                        'basis': '(auto)',
                    },
                    'keywords': {
                        'level_hint': self.dashlevel,
                        'params_tweaks': self.dashparams,
                        'dashcoeff_supplement': self.dashcoeff_supplement,
                        'pair_resolved': self.save_pairwise_disp,
                        'apply_qcengine_aliases': True,  # for s-dftd3
                        'verbose': 1,
                    },
                    'molecule': molecule.to_schema(dtype=2),
                    'provenance': p4util.provenance_stamp(__name__),
                })
            jobrec = qcng.compute(
                resi,
                self.engine,
                raise_error=True,
                task_config={"scratch_directory": core.IOManager.shared_object().get_default_path(), "ncores": core.get_num_threads()})

            dashd_part = float(jobrec.extras['qcvars']['DISPERSION CORRECTION ENERGY'])
            if wfn is not None:
                for k, qca in jobrec.extras['qcvars'].items():
                    if ("CURRENT" not in k) and ("PAIRWISE" not in k):
                        wfn.set_variable(k, float(qca) if isinstance(qca, str) else qca)

                # Pass along the pairwise dispersion decomposition if we need it
                if self.save_pairwise_disp is True:
                    wfn.set_variable("PAIRWISE DISPERSION CORRECTION ANALYSIS",
                                     jobrec.extras['qcvars']["2-BODY PAIRWISE DISPERSION CORRECTION ANALYSIS"])

            if self.fctldash in ['hf3c', 'pbeh3c', 'r2scan3c', 'b973c']:
                jobrec = qcng.compute(
                    resi,
                    self.gcp_engine,
                    raise_error=True,
                    task_config={"scratch_directory": core.IOManager.shared_object().get_default_path(), "ncores": core.get_num_threads()})
                gcp_part = jobrec.return_result
                dashd_part += gcp_part

            return dashd_part

        else:
            ene = self.disp.compute_energy(molecule)
            core.set_variable('DISPERSION CORRECTION ENERGY', ene)
            if self.fctldash:
                core.set_variable(f"{self.fctldash} DISPERSION CORRECTION ENERGY", ene)
            return ene

    def compute_gradient(self,
                         molecule: core.Molecule,
                         wfn: core.Wavefunction = None) -> core.Matrix:
        """Compute dispersion gradient based on engine, dispersion level, and parameters in `self`.

        Parameters
        ----------
        molecule
            System for which to compute empirical dispersion correction.
        wfn
            Location to set QCVariables

        Returns
        -------
        Matrix
            (nat, 3) dispersion gradient [Eh/a0].

        """
        if self.engine in ["s-dftd3", 'dftd3', 'mp2d', "dftd4"]:
            resi = AtomicInput(
                **{
                    'driver': 'gradient',
                    'model': {
                        'method': self.fctldash,
                        'basis': '(auto)',
                    },
                    'keywords': {
                        'level_hint': self.dashlevel,
                        'params_tweaks': self.dashparams,
                        'dashcoeff_supplement': self.dashcoeff_supplement,
                        'apply_qcengine_aliases': True,  # for s-dftd3
                        'verbose': 1,
                    },
                    'molecule': molecule.to_schema(dtype=2),
                    'provenance': p4util.provenance_stamp(__name__),
                })
            jobrec = qcng.compute(
                resi,
                self.engine,
                raise_error=True,
                task_config={"scratch_directory": core.IOManager.shared_object().get_default_path(), "ncores": core.get_num_threads()})

            dashd_part = core.Matrix.from_array(jobrec.extras['qcvars']['DISPERSION CORRECTION GRADIENT'])
            if wfn is not None:
                for k, qca in jobrec.extras['qcvars'].items():
                    if "CURRENT" not in k:
                        wfn.set_variable(k, float(qca) if isinstance(qca, str) else qca)

            if self.fctldash in ['hf3c', 'pbeh3c', 'r2scan3c', 'b973c']:
                jobrec = qcng.compute(
                    resi,
                    self.gcp_engine,
                    raise_error=True,
                    task_config={"scratch_directory": core.IOManager.shared_object().get_default_path(), "ncores": core.get_num_threads()})
                gcp_part = core.Matrix.from_array(jobrec.return_result)
                dashd_part.add(gcp_part)

            return dashd_part
        else:
            return self.disp.compute_gradient(molecule)

    def compute_hessian(self,
                        molecule: core.Molecule,
                        wfn: core.Wavefunction = None) -> core.Matrix:
        """Compute dispersion Hessian based on engine, dispersion level, and parameters in `self`.
        Uses finite difference, as no dispersion engine has analytic second derivatives.

        Parameters
        ----------
        molecule
            System for which to compute empirical dispersion correction.
        wfn
            Location to set QCVariables

        Returns
        -------
        Matrix
            (3*nat, 3*nat) dispersion Hessian [Eh/a0/a0].

        """
        from psi4.driver.driver_findif import assemble_hessian_from_gradients, hessian_from_gradients_geometries

        optstash = p4util.OptionsState(['PRINT'], ['PARENT_SYMMETRY'])
        core.set_global_option('PRINT', 0)

        core.print_out("\n\n   Analytical Dispersion Hessians are not supported by any engine.\n")
        core.print_out("       Computing the Hessian through finite difference of gradients.\n\n")

        # Setup the molecule
        molclone = molecule.clone()
        molclone.reinterpret_coordentry(False)
        molclone.fix_orientation(True)
        molclone.fix_com(True)

        # Record undisplaced symmetry for projection of diplaced point groups
        core.set_global_option("PARENT_SYMMETRY", molecule.schoenflies_symbol())

        findif_meta_dict = hessian_from_gradients_geometries(molclone, -1)
        for displacement in findif_meta_dict["displacements"].values():
            geom_array = np.reshape(displacement["geometry"], (-1, 3))
            molclone.set_geometry(core.Matrix.from_array(geom_array))
            molclone.update_geometry()
            displacement["gradient"] = self.compute_gradient(molclone).np.ravel().tolist()

        H = assemble_hessian_from_gradients(findif_meta_dict, -1)
        if wfn is not None:
            wfn.set_variable('DISPERSION CORRECTION HESSIAN', H)
        optstash.restore()
        return core.Matrix.from_array(H)

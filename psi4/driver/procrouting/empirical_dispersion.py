#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2022 The Psi4 Developers.
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
from qcelemental.models import AtomicInput
import qcengine as qcng

from psi4 import core
from psi4.driver import p4util
from psi4.driver import driver_findif
from psi4.driver.p4util.exceptions import ValidationError

_engine_can_do = collections.OrderedDict([('libdisp', ['d1', 'd2', 'chg', 'das2009', 'das2010']),
                                          ('dftd3', ['d2', 'd3zero', 'd3bj', 'd3mzero', 'd3mbj']),
                                          ('nl', ['nl']),
                                          ('mp2d', ['dmp2']),
                                          ("dftd4", ["d4bjeeqatm"]),
                                        ]) # yapf: disable

_capable_engines_for_disp = collections.defaultdict(list)
for eng, disps in _engine_can_do.items():
    for disp in disps:
        _capable_engines_for_disp[disp].append(eng)


class EmpiricalDispersion():
    """Lightweight unification of empirical dispersion calculation modes.

    Attributes
    ----------
    dashlevel : str
        {'d1', 'd2', 'd3zero', 'd3bj', 'd3mzero', 'd3mbj', 'chg', 'das2009', 'das2010', 'nl', 'dmp2', "d4bjeeqatm"}
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
        {'libdisp', 'dftd3', 'nl', 'mp2d', "dftd4"}
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
        values. Really only relevant for -D2, which can be computed by
        libdisp or dftd3.

    """
    def __init__(self, *, name_hint: str = None, level_hint: str = None, param_tweaks: Union[Dict, List] = None, engine: str = None, save_pairwise_disp=False):
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

        if engine is None:
            self.engine = _capable_engines_for_disp[self.dashlevel][0]
        else:
            if self.dashlevel in _engine_can_do[engine]:
                self.engine = engine
            else:
                raise ValidationError("""This little engine ({}) can't ({})""".format(engine, self.dashlevel))

        if self.engine == 'libdisp':
            self.disp = core.Dispersion.build(self.dashlevel, **resolved['dashparams'])

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
        if self.engine in ['dftd3', 'mp2d', "dftd4"]:
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

            if self.fctldash in ['hf3c', 'pbeh3c']:
                jobrec = qcng.compute(
                    resi,
                    "gcp",
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
        if self.engine in ['dftd3', 'mp2d', "dftd4"]:
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

            if self.fctldash in ['hf3c', 'pbeh3c']:
                jobrec = qcng.compute(
                    resi,
                    "gcp",
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

        findif_meta_dict = driver_findif.hessian_from_gradients_geometries(molclone, -1)
        for displacement in findif_meta_dict["displacements"].values():
            geom_array = np.reshape(displacement["geometry"], (-1, 3))
            molclone.set_geometry(core.Matrix.from_array(geom_array))
            molclone.update_geometry()
            displacement["gradient"] = self.compute_gradient(molclone).np.ravel().tolist()

        H = driver_findif.assemble_hessian_from_gradients(findif_meta_dict, -1)
        if wfn is not None:
            wfn.set_variable('DISPERSION CORRECTION HESSIAN', H)
        optstash.restore()
        return core.Matrix.from_array(H)

#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2019 The Psi4 Developers.
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
"""Module with functions that interface with Grimme's DFTD3 code."""

import re
import sys
import copy
import json
import pprint
pp = pprint.PrettyPrinter(width=120)
from decimal import Decimal

import numpy as np

import qcelemental as qcel

from .. import __version__
#from .. import qcvars
from ..util import parse_dertype
from ..pdict import PreservingDict
from ..exceptions import *
from ..datastructures import QCAspect, print_variables
from . import dashparam
from .worker import dftd3_subprocess
"""Compute dispersion correction using Grimme's DFTD3 executable."""


def run_dftd3(name, molecule, options, **kwargs):
    """QCDriver signature for computing `name` on `molecule` with `options` with engine `DFTD3`."""

    # * ONLY takes self-sufficient fctl-dash for name
    # * tweakparams are only valid options

    opts = {}

    jobrec = {}
    jobrec['error'] = ''
    jobrec['success'] = None
    jobrec['raw_output'] = None
    prov = {}
    prov['creator'] = 'QCDB'
    prov['version'] = __version__
    prov['routine'] = sys._getframe().f_code.co_name
    jobrec['provenance'] = prov

    # strip engine hint
    if name.startswith('d3-'):
        name = name[3:]

    jobrec['molecule'] = molecule.to_dict(np_out=False)
    jobrec['method'] = name
    _, jobrec['driver'] = parse_dertype(kwargs['ptype'], max_derivative=1)
    jobrec['options'] = opts
    #jobrec['options'] = copy.deepcopy(options)

    try:
        dftd3_driver(jobrec)
    except Exception as err:
        jobrec['success'] = False
        jobrec['error'] += repr(err)
    else:
        jobrec['success'] = True
        jobrec['qcvars']['CURRENT ENERGY'] = copy.deepcopy(jobrec['qcvars']['DISPERSION CORRECTION ENERGY'])
        
    return jobrec


def run_dftd3_from_arrays(molrec,
                          name_hint=None,
                          level_hint=None,
                          param_tweaks=None,
                          ptype='energy',
                          dashcoeff_supplement=None,
                          verbose=1):
    """Specialized signature disentangling dispersion level and
    parameters for computing on `molecule` with engine `DFTD3`. See
    `dashparam.from_array` for parameter details.

    """
    jobrec = {}
    jobrec['error'] = ''
    jobrec['success'] = None
    jobrec['return_output'] = True
    prov = {}
    prov['creator'] = 'QCDB'
    prov['version'] = __version__
    prov['routine'] = sys._getframe().f_code.co_name
    jobrec['provenance'] = prov

    # strip engine hint
    if name_hint.startswith('d3-'):
        name_hint = name_hint[3:]

    opts = {}
    opts['level_hint'] = level_hint
    opts['params_tweaks'] = param_tweaks
    opts['dashcoeff_supplement'] = dashcoeff_supplement

    jobrec['molecule'] = molrec
    jobrec['method'] = name_hint
    _, jobrec['driver'] = parse_dertype(ptype, max_derivative=1)
    jobrec['options'] = opts
    #jobrec['options'] = copy.deepcopy(options)

    try:
        dftd3_driver(jobrec)
    except Exception as err:
        jobrec['success'] = False
        jobrec['error'] += repr(err)
        raise RuntimeError(err) from err
    else:
        jobrec['success'] = True
        if ptype == "energy":
            jobrec['qcvars']['CURRENT ENERGY'] = copy.deepcopy(jobrec['qcvars']['DISPERSION CORRECTION ENERGY'])
        
    return jobrec


def dftd3_driver(jobrec, verbose=1):
    """Drive the jobrec@i (input) -> dftd3rec@i -> dftd3rec@io -> jobrec@io (returned) process.

    Input Fields
    ------------

    Optional Input Fields
    ---------------------

    Output Fields
    -------------

    Optional Output Fields
    ----------------------

    """
    if verbose > 2:
        print('[1] {} JOBREC PRE-PLANT (j@i) <<<'.format('DFTD3'))
        pp.pprint(jobrec)
        print('>>>')

    dftd3rec = dftd3_plant(jobrec)

    # test json roundtrip
    jdftd3rec = json.dumps(dftd3rec)
    dftd3rec = json.loads(jdftd3rec)

    if verbose > 3:
        print('[2] {}REC PRE-SUBPROCESS (m@i) <<<'.format('DFTD3'))
        pp.pprint(dftd3rec)
        print('>>>\n')

    dftd3_subprocess(dftd3rec)  # updates dftd3rec

    if verbose > 3:
        print('[3] {}REC POST-SUBPROCESS (m@io) <<<'.format('DFTD3'))
        pp.pprint(dftd3rec)
        print('>>>\n')

    dftd3_harvest(jobrec, dftd3rec)  # updates jobrec

    if verbose > 1:
        print('[4] {} JOBREC POST-HARVEST (j@io) <<<'.format('DFTD3'))
        pp.pprint(jobrec)
        print('>>>')

    return jobrec


def dftd3_plant(jobrec):
    """Transform the QC input specifications `jobrec` into the command
    and files for DFTD3: jobrec@i -> dftd3rec@i.

    Parameters
    ----------
    jobrec : dict
        Nested dictionary with input specifications for DFTD3 in generic
        QC terms.

    Returns
    -------
    dftd3rec : dict
        Nested dictionary with input specification for DFTD3 in
        program-specific commands and files.

    """
    try:
        jobrec['driver']
        jobrec['method']
        jobrec['options']
        jobrec['molecule']
    except KeyError as err:
        raise KeyError('Required field ({}) missing from ({})'.format(str(err), list(jobrec.keys()))) from err

    # temp until actual options object
    dftd3rec = dashparam.from_arrays(
        name_hint=jobrec['method'],
        level_hint=jobrec['options'].get('level_hint', None),
        param_tweaks=jobrec['options'].get('params_tweaks', None),
        dashcoeff_supplement=jobrec['options'].get('dashcoeff_supplement', None))
    # sketchy: adding to options during planting season
    jobrec['options'].update(dftd3rec)

    # this is what the dftd3 program needs, not what the job needs
    # * form dftd3_parameters string that governs dispersion calc
    # * form dftd3_geometry string that supplies geometry to dispersion calc
    # * form command and arguments

    dftd3rec['dftd3par'] = dftd3_coeff_formatter(dftd3rec['dashlevel'], dftd3rec['dashparams'])

    dftd3rec['dftd3_geometry'] = qcel.molparse.to_string(jobrec['molecule'], dtype='xyz', units='Angstrom', ghost_format='')

    command = ['dftd3', 'dftd3_geometry.xyz']
    if jobrec['driver'] == 'gradient':
        command.append('-grad')
    dftd3rec['command'] = command

    return dftd3rec


def dftd3_harvest(jobrec, dftd3rec):
    """Process raw results from read-only `dftd3rec` into QCAspect
    fields in returned `jobrec`: jobrec@i, dftd3rec@io -> jobrec@io.

    Parameters
    ----------
    jobrec : dict
        Nested dictionary with input specifications for DFTD3 in generic
        QC terms.
    dftd3rec : dict
        Nested dictionary with input specification and output collection
        from DFTD3 in program-specific commands, files, & output capture.

    Returns
    -------
    jobrec : dict
        Nested dictionary with input specification and output collection
        from DFTD3 in generic QC terms.

    """
    try:
        jobrec['molecule']['real']
        jobrec['driver']
        jobrec['provenance']
        jobrec['options']['fctldash']
    except KeyError as err:
        raise KeyError('Required field ({}) missing from ({})'.format(str(err), list(jobrec.keys()))) from err

    try:
        dftd3rec['stdout']
    except KeyError as err:
        raise KeyError('Required field ({}) missing from ({})'.format(str(err), list(dftd3rec.keys()))) from err

    # amalgamate output
    text = dftd3rec['stdout']
    text += '\n  <<<  DFTD3 Results  >>>\n'

    for fl in ['dftd3_gradient']:
        field = 'output_' + fl.lower()
        if field in dftd3rec:
            text += '\n  DFTD3 scratch file {} has been read.\n'.format(fl)
            text += dftd3rec[field]

    # parse energy output (could go further and break into E6, E8, E10 and Cn coeff)
    real = np.array(jobrec['molecule']['real'])
    full_nat = real.shape[0]
    real_nat = np.sum(real)

    for ln in dftd3rec['stdout'].splitlines():
        if re.search('DFTD3 V', ln):
            version = ln.replace('DFTD3', '').replace('|', '').strip().lower()
        elif re.match(' Edisp /kcal,au', ln):
            ene = Decimal(ln.split()[3])
        elif re.match(' normal termination of dftd3', ln):
            break
    else:
        if not ((real_nat == 1) and (jobrec['driver'] == 'gradient')):
            raise Dftd3Error('Unsuccessful run. Possibly -D variant not available in dftd3 version.')

    # parse gradient output
    # * DFTD3 crashes on one-atom gradients. Avoid the error (above) and just force the correct result (below).
    if 'output_dftd3_gradient' in dftd3rec:
        srealgrad = dftd3rec['output_dftd3_gradient'].replace('D', 'E')
        realgrad = np.fromstring(srealgrad, count=3 * real_nat, sep=' ').reshape((-1, 3))
    elif real_nat == 1:
        realgrad = np.zeros((1, 3))

    if jobrec['driver'] == 'gradient':
        ireal = np.argwhere(real).reshape((-1))
        fullgrad = np.zeros((full_nat, 3))
        try:
            fullgrad[ireal, :] = realgrad
        except NameError as err:
            raise Dftd3Error('Unsuccessful gradient collection.') from err

    qcvkey = jobrec['options']['fctldash'].upper()

    # OLD WAY
    calcinfo = []
    calcinfo.append(QCAspect('DISPERSION CORRECTION ENERGY', 'Eh', ene, ''))
    if qcvkey:
        calcinfo.append(QCAspect('{} DISPERSION CORRECTION ENERGY'.format(qcvkey), 'Eh', ene, ''))
    if jobrec['driver'] == 'gradient':
        calcinfo.append(QCAspect('DISPERSION CORRECTION GRADIENT', 'Eh/a0', fullgrad, ''))
        if qcvkey:
            calcinfo.append(QCAspect('{} DISPERSION CORRECTION GRADIENT'.format(qcvkey), 'Eh/a0', fullgrad, ''))
    calcinfo = {info.lbl: info for info in calcinfo}
    text += print_variables(calcinfo)

    # NEW WAY
    #module_vars = {}
    #module_vars['DISPERSION CORRECTION ENERGY'] =  ene
    #module_vars['{} DISPERSION CORRECTION ENERGY'.format(qcvkey)] =  ene
    #if jobrec['driver'] == 'gradient':
    #    module_vars['DISPERSION CORRECTION GRADIENT'] = fullgrad
    #    module_vars['{} DISPERSION CORRECTION GRADIENT'.format(qcvkey)] = fullgrad
    #
    #module_vars = PreservingDict(module_vars)
    #qcvars.build_out(module_vars)
    #calcinfo = qcvars.certify(module_vars)
    #text += print_variables(calcinfo)

    jobrec['raw_output'] = text
    jobrec['qcvars'] = calcinfo

    prov = {}
    prov['creator'] = 'DFTD3'
    prov['routine'] = sys._getframe().f_code.co_name
    prov['version'] = version
    jobrec['provenance'] = prov

    return jobrec


def dftd3_coeff_formatter(dashlvl, dashcoeff):
    """Return strings for DFTD3 program parameter file.

             s6 rs6 s18    rs8     alpha6      version
             -------------------------------------------
    d2:      s6 sr6 s8=0.0 a2=None alpha6      version=2
    d3zero:  s6 sr6 s8     a2=sr8  alpha6      version=3
    d3bj:    s6 a1  s8     a2      alpha6=None version=4
    d3mzero: s6 sr6 s8     beta    alpha6=14.0 version=5
    d3mbj:   s6 a1  s8     a2      alpha6=None version=6

    Parameters
    ----------
    dashlvl : {'d2', 'd3zero', d3bj', 'd3mzero', 'd3mbj'}
        Level of dispersion correction.
    dashcoeff : dict
        Dictionary fully specifying non-fixed parameters (table above) for `dashlvl` to drive DFTD3.

    Returns
    -------
    str
        Suitable for `.dftd3par` file.

    """
    dashformatter = """{:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:6}\n"""

    dashlvl = dashlvl.lower()
    if dashlvl == 'd2':
        return dashformatter.format(dashcoeff['s6'], dashcoeff['sr6'], 0.0, 0.0, dashcoeff['alpha6'], 2)
    elif dashlvl == 'd3zero':
        return dashformatter.format(dashcoeff['s6'], dashcoeff['sr6'], dashcoeff['s8'], dashcoeff['sr8'],
                                    dashcoeff['alpha6'], 3)
    elif dashlvl == 'd3bj':
        return dashformatter.format(dashcoeff['s6'], dashcoeff['a1'], dashcoeff['s8'], dashcoeff['a2'], 0.0, 4)
    elif dashlvl == 'd3mzero':
        return dashformatter.format(dashcoeff['s6'], dashcoeff['sr6'], dashcoeff['s8'], dashcoeff['beta'], 14.0, 5)
    elif dashlvl == 'd3mbj':
        return dashformatter.format(dashcoeff['s6'], dashcoeff['a1'], dashcoeff['s8'], dashcoeff['a2'], 0.0, 6)
    else:
        raise ValidationError(
            """-D correction level %s is not available. Choose among %s.""" % (dashlvl, dashcoeff.keys()))


"""
Notes
-----
The DFTD3 executable must be independently compiled and found in :envvar:`PATH` or :envvar:`PSIPATH`.
research site: https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/dft-d3
Psi4 mode: When `psi4` the python module is importable at `import qcdb`
           time, Psi4 mode is activated, with the following alterations:
           * output goes to output file
           * gradient returned as psi4.core.Matrix, not list o'lists
           * scratch is written to randomly named subdirectory of psi scratch
           * psivar "DISPERSION CORRECTION ENERGY" is set
           * `verbose` triggered when PRINT keywork of SCF module >=3
"""

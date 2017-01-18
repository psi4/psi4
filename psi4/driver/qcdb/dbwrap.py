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

from __future__ import absolute_import
from __future__ import print_function
import os
import sys
import math

try:
    import cPickle as pickle
except ImportError:
    import pickle
import itertools
# from collections import defaultdict
try:
    from collections import OrderedDict
except ImportError:
    from oldpymodules import OrderedDict
from .exceptions import *
from .molecule import Molecule
from .modelchems import Method, BasisSet, Error, methods, bases, errors, pubs
from . import psiutil
from . import textables


def initialize_errors():
    """Form OrderedDict of all possible statistical measures set to None"""
    error = OrderedDict()
    for e in ['e', 'pe', 'pbe', 'pce']:
        for m in ['pex', 'nex', 'max', 'min', 'm', 'ma', 'rms', 'std']:
            error[m + e] = None
    return error


def initialize_errors_elaborate(e=None, pe=None, pbe=None, pce=None, extrema=True):
    error = OrderedDict()
    error['maxe'] = None if (e is None or not extrema) else e  # LD_XA
    error['mine'] = None if (e is None or not extrema) else e  # LD_XI
    error['me'] = None if e is None else 0.0  # LD_MS
    error['mae'] = None if e is None else 0.0  # LD_MA
    error['rmse'] = None if e is None else 0.0  # LD_RA
    error['stde'] = None if e is None else 0.0
    error['maxpe'] = None if (pe is None or not extrema) else pe  # FD_XA
    error['minpe'] = None if (pe is None or not extrema) else pe  # FD_XI
    error['mpe'] = None if pe is None else 0.0  # FD_MS
    error['mape'] = None if pe is None else 0.0  # FD_MA
    error['rmspe'] = None if pe is None else 0.0  # FD_RA
    error['stdpe'] = None if pe is None else 0.0
    error['maxpbe'] = None if (pbe is None or not extrema) else pbe  # BD_XA
    error['minpbe'] = None if (pbe is None or not extrema) else pbe  # BD_XI
    error['mpbe'] = None if pbe is None else 0.0  # BD_MS
    error['mapbe'] = None if pbe is None else 0.0  # BD_MA
    error['rmspbe'] = None if pbe is None else 0.0  # BD_RA
    error['stdpbe'] = None if pbe is None else 0.0
    error['maxpce'] = None if (pce is None or not extrema) else pce  # BD_XA
    error['minpce'] = None if (pce is None or not extrema) else pce  # BD_XI
    error['mpce'] = None if pce is None else 0.0  # BD_MS
    error['mapce'] = None if pce is None else 0.0  # BD_MA
    error['rmspce'] = None if pce is None else 0.0  # BD_RA
    error['stdpce'] = None if pce is None else 0.0
    return error


def average_errors(*args):
    """Each item in *args* should be an error dictionary. Performs
    average-like operation over all items, which should be error
    dictionaries, in *args*. Defined for ME, MAE, STDE, and their
    relative-error variants. None returned for undefined statistics or
    when an item is missing.

    """
    Ndb = float(len(args))
    avgerror = initialize_errors()
    try:
        avgerror['pexe'] = max([x['pexe'] for x in args])
        avgerror['nexe'] = min([x['nexe'] for x in args])
        avgerror['maxe'] = max([x['maxe'] for x in args], key=lambda x: abs(x))
        avgerror['mine'] = min([x['mine'] for x in args], key=lambda x: abs(x))
        avgerror['me'] = sum([x['me'] for x in args]) / Ndb
        avgerror['mae'] = sum([x['mae'] for x in args]) / Ndb
        avgerror['rmse'] = sum([x['rmse'] for x in args]) / Ndb  # TODO: unsure of op validity
        avgerror['stde'] = math.sqrt(sum([x['stde'] ** 2 for x in args]) / Ndb)

        avgerror['pexpe'] = max([x['pexpe'] for x in args])
        avgerror['nexpe'] = min([x['nexpe'] for x in args])
        avgerror['maxpe'] = max([x['maxpe'] for x in args], key=lambda x: abs(x))
        avgerror['minpe'] = min([x['minpe'] for x in args], key=lambda x: abs(x))
        avgerror['mpe'] = sum([x['mpe'] for x in args]) / Ndb
        avgerror['mape'] = sum([x['mape'] for x in args]) / Ndb
        avgerror['rmspe'] = sum([x['rmspe'] for x in args]) / Ndb  # TODO: unsure of op validity
        avgerror['stdpe'] = math.sqrt(sum([x['stdpe'] * x['stdpe'] for x in args]) / Ndb)

        avgerror['pexpbe'] = max([x['pexpbe'] for x in args])
        avgerror['nexpbe'] = min([x['nexpbe'] for x in args])
        avgerror['maxpbe'] = max([x['maxpbe'] for x in args], key=lambda x: abs(x))
        avgerror['minpbe'] = min([x['minpbe'] for x in args], key=lambda x: abs(x))
        avgerror['mpbe'] = sum([x['mpbe'] for x in args]) / Ndb
        avgerror['mapbe'] = sum([x['mapbe'] for x in args]) / Ndb
        avgerror['rmspbe'] = sum([x['rmspbe'] for x in args]) / Ndb  # TODO: unsure of op validity
        avgerror['stdpbe'] = math.sqrt(sum([x['stdpbe'] * x['stdpbe'] for x in args]) / Ndb)

        avgerror['pexpce'] = max([x['pexpce'] for x in args])
        avgerror['nexpce'] = min([x['nexpce'] for x in args])
        avgerror['maxpce'] = max([x['maxpce'] for x in args], key=lambda x: abs(x))
        avgerror['minpce'] = min([x['minpce'] for x in args], key=lambda x: abs(x))
        avgerror['mpce'] = sum([x['mpce'] for x in args]) / Ndb
        avgerror['mapce'] = sum([x['mapce'] for x in args]) / Ndb
        avgerror['rmspce'] = sum([x['rmspce'] for x in args]) / Ndb  # TODO: unsure of op validity
        avgerror['stdpce'] = math.sqrt(sum([x['stdpce'] * x['stdpce'] for x in args]) / Ndb)
    except TypeError:
        pass
    return avgerror


def format_errors(err, mode=1):
    """From error dictionary *err*, returns a LaTeX-formatted string,
    after handling None entries.

    """
    onedecimal = r"""{0:8.1f}"""
    twodecimal = r"""{0:8.2f}"""
    threedecimal = r"""{0:12.3f}"""
    fourdecimal = r"""{0:12.4f}"""
    shortblank = r"""{0:8s}""".format('')
    longblank = r"""{0:12s}""".format('')

    if mode == 1:
        me = ' ----' if err['me'] is None else '%+.2f' % (err['me'])
        stde = '----' if err['stde'] is None else '%.2f' % (err['stde'])
        mae = '  ----' if err['mae'] is None else '%6.2f' % (err['mae'])
        mape = '  ----  ' if err['mape'] is None else '%6.1f\%%' % (100 * err['mape'])
        mapbe = '  ----  ' if err['mapbe'] is None else '%6.1f\%%' % (100 * err['mapbe'])
        mapce = '  ----  ' if err['mapce'] is None else '%6.1f\%%' % (100 * err['mapce'])
        text = """$\{%s; %s\}$ %s %s %s""" % \
               (me, stde, mae, mape, mapce)
        return text

    if mode == 2:
        sdict = OrderedDict()
        for lbl in ['pexe', 'nexe', 'maxe', 'mine', 'me', 'mae', 'rmse', 'stde']:
            sdict[lbl] = '        ----' if err[lbl] is None else fourdecimal.format(err[lbl])
        for lbl in ['pexpe', 'nexpe', 'maxpe', 'minpe', 'mpe', 'mape', 'rmspe', 'stdpe',
                    'pexpbe', 'nexpbe', 'maxpbe', 'minpbe', 'mpbe', 'mapbe', 'rmspbe', 'stdpbe',
                    'pexpce', 'nexpce', 'maxpce', 'minpce', 'mpce', 'mapce', 'rmspce', 'stdpce']:
            sdict[lbl] = '        ----' if err[lbl] is None else threedecimal.format(100 * err[lbl])
        text = """nex: {nexe}{nexpe}{nexpbe}{nexpce}\n""" \
               """pex: {pexe}{pexpe}{pexpbe}{pexpce}\n""" \
               """min: {mine}{minpe}{minpbe}{minpce}\n""" \
               """max: {maxe}{maxpe}{maxpbe}{maxpce}\n""" \
               """m:   {me}{mpe}{mpbe}{mpce}\n""" \
               """ma:  {mae}{mape}{mapbe}{mapce}\n""" \
               """rms: {rmse}{rmspe}{rmspbe}{rmspce}\n""" \
               """std: {stde}{stdpe}{stdpbe}{stdpce}\n""".format(**sdict)
        return text

    if mode == 3:
        sdict = OrderedDict()
        # shortblanks changed from empty strings Aug 2015
        for lbl in ['pexe', 'nexe', 'maxe', 'mine', 'me', 'mae', 'rmse', 'stde']:
            sdict[lbl] = shortblank if err[lbl] is None else twodecimal.format(err[lbl])
        for lbl in ['pexpe', 'nexpe', 'maxpe', 'minpe', 'mpe', 'mape', 'rmspe', 'stdpe',
                    'pexpbe', 'nexpbe', 'maxpbe', 'minpbe', 'mpbe', 'mapbe', 'rmspbe', 'stdpbe',
                    'pexpce', 'nexpce', 'maxpce', 'minpce', 'mpce', 'mapce', 'rmspce', 'stdpce']:
            sdict[lbl] = shortblank if err[lbl] is None else onedecimal.format(100 * err[lbl])
        return sdict


def string_contrast(ss):
    """From an array of strings, *ss*, returns maximum common prefix
    string, maximum common suffix string, and array of middles.

    """
    s = [item + 'q' for item in ss if item is not None]
    short = min(s, key=len)
    for ib in range(len(short)):
        if not all([mc[ib] == short[ib] for mc in s]):
            preidx = ib
            break
    else:
        preidx = 0
    for ib in range(len(short)):
        ie = -1 * (ib + 1)
        if not all([mc[ie] == short[ie] for mc in s]):
            sufidx = ie + 1
            break
    else:
        sufidx = -1 * (len(short))

    miditer = iter([mc[preidx:sufidx] for mc in s])
    prefix = short[:preidx]
    suffix = short[sufidx:-1]
    middle = ['' if mc is None else next(miditer) for mc in ss]

    return prefix, suffix, middle


def oxcom(lst):
    """Returns gramatical comma separated string of *lst*."""
    lst = [str(l) for l in lst]

    if not lst:
        return ''
    elif len(lst) == 1:
        return lst[0]
    elif len(lst) == 2:
        return ' and '.join(lst)
    else:
        return ', and '.join([', '.join(lst[:-1]), lst[-1]])


def cure_weight(refrxn, refeq, rrat, xi=0.2):
    """
    :param refeq: value of benchmark for equilibrium Reaction
    :param rrat: ratio of intermonomer separation for Reaction to equilibrium Reaction
    :param xi: parameter
    :return: weight for CURE

    """
    sigma = xi * abs(refeq) / (rrat ** 3)
    weight = max(abs(refrxn), sigma)
    return weight


def balanced_error(refrxn, refeq, rrat, m=0.03, p=10.0):
    """
    :param refrxn:
    :param refeq:
    :param rrat:
    :param m: minimum permitted weight for a point
    :param p: multiples of abs(refeq) above refeq to which zero-line in head is displaced
    :return:

    """
    one = float(1)
    q = one if rrat >= one else p
    qm1perat = q - 1 + refrxn / refeq
    weight = max(m, qm1perat / q)
    mask = weight * q / abs(qm1perat)
    return mask, weight


def fancify_mc_tag(mc, latex=False):
    """From the usual MTD-opt1_opt2-bas model chemistry identifier, return
    string based on fullname, if *latex* is False or latex if *latex* is True.

    """
    try:
        mtd, mod, bas = mc.split('-')
    except ValueError:
        text = mc
    else:
        if latex:
            text = r"""%20s / %-20s %s""" % (methods[mtd].latex, bases[bas].latex, mod)
        else:
            text = r"""%20s / %s, %s""" % (methods[mtd].fullname, bases[bas].fullname, mod)
    return text


class ReactionDatum(object):
    """Piece of quantum chemical information that describes a qcdb.Reaction object.

    """

    def __init__(self, dbse, rxn, method, mode, basis, value, units='kcal/mol', citation=None, doi=None, comment=None):
        # geometry
        self.dbrxn = dbse + '-' + str(rxn)
        # qcdb.Method
        self.method = method
        # mode, e.g., unCP, CP, RLX, etc.
        self.mode = mode
        # qcdb.BasisSet
        self.basis = basis
        # numerical value for reaction
        self.value = float(value)
        # energy unit attached to value, defaults to kcal/mol
        self.units = units
        # publication citation of value
        self.citation = citation
        # digital object identifier for publication (maybe this should be doi of datum, not of pub?)
        self.doi = doi
        # addl comments
        self.comment = comment

    @classmethod
    def library_modelchem(cls, dbse, rxn, method, mode, basis, value, units='kcal/mol', citation=None, doi=None,
                          comment=None):
        """Constructor when method and basis are strings corresponding to
        qcdb.Method and qcdb.BasisSet already defined in methods and bases.

        """
        # computational method
        try:
            tmp_method = methods[method.upper()]
        except KeyError as e:
            raise ValidationError("""Invalid ReactionDatum method %s: %s""" % (method, e))
        # computational basis set
        try:
            tmp_basis = bases[basis.lower()]
        except KeyError as e:
            raise ValidationError("""Invalid ReactionDatum basis %s: %s""" % (basis, e))
        # publication
        if citation is None:
            tmp_pub = citation
        else:
            try:
                tmp_pub = pubs[citation.lower()]
            except KeyError as e:
                raise ValidationError("""Invalid ReactionDatum publication %s: %s""" % (citation, e))
        return cls(dbse, rxn, tmp_method, mode, tmp_basis, value, units, citation=tmp_pub, doi=doi, comment=comment)

    def __str__(self):
        text = ''
        text += """  ==> ReactionDatum <==\n\n"""
        text += """  Database reaction:    %s\n""" % (self.dbrxn)
        text += """  Method:               %s\n""" % (self.method.fullname)
        text += """  Mode:                 %s\n""" % (self.mode)
        text += """  Basis:                %s\n""" % (self.basis.fullname)
        text += """  Value:                %f [%s]\n""" % (self.value, self.units)
        text += """  Citation:             %s %s\n""" % (self.citation.name, self.citation.doi)
        text += """  DOI:                  %s\n""" % (self.doi)
        text += """  Comment:              %s\n""" % (self.comment)
        text += """\n"""
        return text


class Subset(object):
    """Affiliated qcdb.Reaction-s

    """

    def __init__(self, name, hrxn, tagl=None, axis=None):
        # identifier
        self.name = name
        # array of reactions names
        self.hrxn = hrxn
        # description line
        self.tagl = tagl
        # mathematical relationships of reactions
        self.axis = OrderedDict()

    def __str__(self):
        text = ''
        text += """  ==> %s Subset <==\n\n""" % (self.name)
        text += """  Tagline:              %s\n""" % (self.tagl)
        text += """  %20s""" % ('Reactions')
        for ax in self.axis.keys():
            text += """  %8s""" % (ax)
        text += """\n"""
        for ix in range(len(self.hrxn)):
            text += """  %20s""" % (str(self.hrxn[ix]))
            for ax in self.axis.values():
                text += """  %8.3f""" % (ax[ix])
            text += """\n"""
        text += """\n"""
        return text


class Reagent(object):
    """Chemical entity only slightly dresed up from qcdb.Molecule.

    """

    def __init__(self, name, mol, tagl=None, comment=None):
        # full name, e.g., 'S22-2-dimer' or 'NBC1-BzMe-8.0-monoA-CP' or 'HTBH-HCl-reagent'
        self.name = name
        # qcdb.Molecule
        try:
            self.NRE = mol.nuclear_repulsion_energy()
        except AttributeError:
            raise ValidationError("""Reagent must be instantiated with qcdb.Molecule object.""")
        else:
            self.mol = mol.create_psi4_string_from_molecule()
        # description line
        self.tagl = tagl
        # # addl comments
        # self.comment = comment
        # # fragmentation
        # self.fragments = mol.fragments
        # # frag activation
        # self.frtype = mol.fragment_types
        # # frag charge
        # self.frchg = mol.fragment_charges
        # # frag multiplicity
        # self.frmult = mol.fragment_multiplicities
        self.charge = mol.molecular_charge()

    def __str__(self):
        text = ''
        text += """  ==> %s Reagent <==\n\n""" % (self.name)
        text += """  Tagline:              %s\n""" % (self.tagl)
        # text += """  Comment:              %s\n""" % (self.comment)
        text += """  NRE:                  %f\n""" % (self.NRE)
        # text += """  Charge:               %+d\n"""
        # text += """  Fragments:            %d\n""" % (len(self.fragments))
        # text += """    FrgNo  Actv  Chg  Mult  AtomRange\n"""
        # for fr in range(len(self.fragments)):
        #    text += """    %-4d   %1s     %+2d  %2d     %s\n""" % (fr + 1,
        #        '*' if self.frtype[fr] == 'Real' else '',
        #        self.frchg[fr], self.frmult[fr], self.fragments[fr])
        text += """  Molecule:             \n%s""" % (self.mol)
        text += """\n"""
        return text


class Reaction(object):
    """

    """

    def __init__(self, name, dbse, indx, tagl=None, latex=None, color='black', comment=None):
        # name, e.g., '2' or 'BzMe-8.0'
        self.name = name
        # database reaction name, e.g., 'S22-2' or 'NBC1-BzMe-8.0'
        self.dbrxn = dbse + '-' + str(name)
        # numerical index of reaction
        self.indx = indx
        # description line
        self.tagl = tagl
        # latex description
        self.latex = latex
        # addl comments
        self.comment = comment
        # reaction matrices, specifying reagent contributions per reaction
        self.rxnm = {}
        # qcdb.ReactionDatum objects of quantum chemical data pertaining to reaction
        self.data = {}
        # benchmark qcdb.ReactionDatum
        self.benchmark = None
        # color for plotting
        self.color = color

    def __str__(self):
        text = ''
        text += """  ==> %s Reaction <==\n\n""" % (self.name)
        text += """  Database reaction:    %s\n""" % (self.dbrxn)
        text += """  Index:                %s\n""" % (self.indx)
        text += """  LaTeX representation: %s\n""" % (self.latex)
        text += """  Tagline:              %s\n""" % (self.tagl)
        text += """  Comment:              %s\n""" % (self.comment)
        if self.benchmark is None:
            text += """  Benchmark:            %s\n""" % ('UNDEFINED')
        else:
            text += """  Benchmark:            %f\n""" % (self.data[self.benchmark].value)
        text += """  Color:                %s\n""" % (str(self.color))
        text += """  Reaction matrix:\n"""
        for mode, rxnm in self.rxnm.items():
            text += """      %s\n""" % (mode)
            for rgt, coeff in rxnm.items():
                text += """       %3d  %s\n""" % (coeff, rgt.name)
        text += """  Data:\n"""
        for label, datum in sorted(self.data.items()):
            text += """      %8.2f  %s\n""" % (datum.value, label)
        text += """\n"""
        return text

    def compute_errors(self, benchmark='default', mcset='default', failoninc=True, verbose=False):
        """For all data or modelchem subset *mcset*, computes raw reaction
        errors between *modelchem* and *benchmark* model chemistries.
        Returns error if model chemistries are missing for any reaction in
        subset unless *failoninc* set to False, whereupon returns partial.
        Returns dictionary of reaction labels and error forms.

        """
        if mcset == 'default':
            lsslist = self.data.keys()
        elif callable(mcset):
            # mcset is function that will generate subset of HRXN from sset(self)
            lsslist = [mc for mc in self.data.keys() if mc in mcset(self)]  # untested
        else:
            # mcset is array containing modelchemistries
            lsslist = [mc for mc in self.data.keys() if mc in mcset]
        # assemble dict of qcdb.Reaction objects from array of reaction names
        lsset = OrderedDict()
        for mc in lsslist:
            lsset[mc] = self.data[mc]

        lbench = self.benchmark if benchmark == 'default' else benchmark
        try:
            mcGreater = self.data[lbench].value
        except KeyError as e:
            raise ValidationError("""Reaction %s missing benchmark datum %s.""" % (self.name, str(e)))

        err = {}
        for label, datum in lsset.items():
            try:
                mcLesser = datum.value
            except KeyError as e:
                if failoninc:
                    raise ValidationError("""Reaction %s missing datum %s.""" % (label, str(e)))
                else:
                    continue

            err[label] = [mcLesser - mcGreater,
                          (mcLesser - mcGreater) / abs(mcGreater),
                          (mcLesser - mcGreater) / abs(mcGreater)]  # TODO define BER
            if verbose:
                print("""p = %6.2f, pe = %6.1f%%, bpe = %6.1f%% modelchem %s.""" %
                      (err[label][0], 100 * err[label][1], 100 * err[label][2], label))

        return err

    def plot(self, benchmark='default', mcset='default',
             failoninc=True, verbose=False, color='sapt',
             xlimit=4.0, labeled=True, view=True,
             mousetext=None, mouselink=None, mouseimag=None, mousetitle=None, mousediv=None,
             saveas=None, relpath=False, graphicsformat=['pdf']):
        """Computes individual errors over model chemistries in *mcset* (which
        may be default or an array or a function generating an array) versus
        *benchmark*. Thread *color* can be 'rgb' for old coloring, a color
        name or 'sapt' for spectrum coloring.

        *saveas* conveys directory ('/') and/or filename for saving the
        resulting plot. File extension is not accessible, but *graphicsformat*
        array requests among 'png', 'pdf', and 'eps' formats. *relpath*
        forces paths to saved files to be relative to current directory,
        rather than absolute paths for returned code and file dictionary.

        Prepares thread diagram instructions and either executes them if
        matplotlib available (Canopy or Anaconda) or prints them. Returns a
        dictionary of all saved plot filenames. If any of *mousetext*, *mouselink*,
        or *mouseimag* is specified, htmlcode will be returned with an image map of
        slats to any of text, link, or image, respectively.

        """
        # compute errors
        dbse = self.dbrxn.split('-')[0]
        indiv = self.compute_errors(benchmark=benchmark, mcset=mcset,
                                    failoninc=failoninc, verbose=verbose)

        # repackage
        dbdat = []
        for mc in indiv.keys():
            dbdat.append({'db': dbse,
                          'show': fancify_mc_tag(mc),
                          'sys': mc,
                          'color': self.color,
                          'data': [indiv[mc][0]]})
        mae = None  # [errors[ix][self.dbse]['mae'] for ix in index]
        mape = None  # [100 * errors[ix][self.dbse]['mape'] for ix in index]
        # form unique filename
        # ixpre, ixsuf, ixmid = string_contrast(index)
        # title = self.dbse + ' ' + ixpre + '[]' + ixsuf
        title = self.dbrxn
        labels = ['']
        # generate matplotlib instructions and call or print
        try:
            from . import mpl
            import matplotlib.pyplot as plt
        except ImportError:
            # if not running from Canopy, print line to execute from Canopy
            print("""filedict, htmlcode = mpl.threads(%s,\n    color='%s',\n    title='%s',\n    labels=%s,\n    mae=%s,\n    mape=%s\n    xlimit=%s\n    labeled=%s\n    saveas=%s\n    mousetext=%s\n    mouselink=%s\n    mouseimag=%s\n    mousetitle=%s,\n    mousediv=%s,\n    relpath=%s\n    graphicsformat=%s)\n\n""" %
                  (dbdat, color, title, labels, mae, mape, str(xlimit),
                   repr(labeled), repr(saveas), repr(mousetext), repr(mouselink), repr(mouseimag),
                   repr(mousetitle), repr(mousediv), repr(relpath), repr(graphicsformat)))
        else:
            # if running from Canopy, call mpl directly
            filedict, htmlcode = mpl.threads(dbdat, color=color, title=title, labels=labels, mae=mae, mape=mape,
                                             xlimit=xlimit, labeled=labeled, view=view,
                                             mousetext=mousetext, mouselink=mouselink,
                                             mouseimag=mouseimag, mousetitle=mousetitle, mousediv=mousediv,
                                             saveas=saveas, relpath=relpath, graphicsformat=graphicsformat)
            return filedict, htmlcode


class WrappedDatabase(object):
    """Wrapper class for raw Psi4 database modules that does some validation
    of contents, creates member data and accessors for database structures,
    defines error computation, and handles database subsets. Not to be used
    directly-- see qcdb.Database for handling single or multiple
    qdcb.WrappedDatabase objects and defining nice statistics, plotting, and
    table functionalities.

    >>> asdf = qcdb.WrappedDatabase('Nbc10')
    """

    def __init__(self, dbname, pythonpath=None):
        """Instantiate class with case insensitive name *dbname*. Module
        search path can be prepended with *pythonpath*.

        """
        #: internal name of database
        #:
        #: >>> print asdf.dbse
        #: 'NBC1'
        self.dbse = None

        #: description line
        #:
        #: >>> print asdf.tagl
        #: 'interaction energies of dissociation curves for non-bonded systems'
        self.tagl = None

        #: OrderedDict of reactions/members
        #:
        #: >>> print asdf.hrxn.keys()
        #: ['BzBz_S-3.2', 'BzBz_S-3.3', ...  'BzBz_PD36-2.8', 'BzBz_PD36-3.0']
        self.hrxn = None

        #: dict of reagents/geometries
        #:
        #: >>> print asdf.hrgt.keys()
        #: ['NBC1-BzBz_PD32-0.8-monoA-CP', 'NBC1-BzBz_PD34-0.6-dimer', ... 'NBC1-BzBz_PD34-1.7-dimer']
        self.hrgt = None

        #: dict of defined reaction subsets.
        #: Note that self.sset['default'] contains all the nonredundant information.
        #:
        #: >>> print asdf.sset.keys()
        #: ['meme', 'mxddpp', '5min', ... 'small']
        self.sset = None

        # Removing hrxn, hrgt etc. do not reduce the size of the object.
        # These attributes are stored for ease of access for adding qc info, etc.

        #: object of defined reaction subsets.
        self.oss = None

        # load database
        if pythonpath is not None:
            sys.path.insert(1, pythonpath)
        else:
            sys.path.append(os.path.dirname(__file__) + '/../databases')
        database = psiutil.import_ignorecase(dbname)
        if not database:
            print('\nPython module for database %s failed to load\n\n' % (dbname))
            print('\nSearch path that was tried:\n')
            print(", ".join(map(str, sys.path)))
            raise ValidationError("Python module loading problem for database " + str(dbname))

        # gross validation of database
        for item in ['dbse', 'GEOS', 'HRXN', 'ACTV', 'RXNM']:
            try:
                getattr(database, item)
            except AttributeError:
                raise ValidationError("""Database %s severely deformed with %s missing.""" % (database.__name__, item))
        for item in ['TAGL', 'BIND']:
            try:
                getattr(database, item)
            except AttributeError:
                print("""Warning: Database %s possibly deformed with %s missing.\n""" % (database.__name__, item))

        # form database name
        self.dbse = database.dbse
        try:
            self.tagl = database.TAGL['dbse']
        except KeyError:
            print("""Warning: TAGL missing for database %s""" % (self.dbse))

        # form array of database contents to process through
        pieces = []
        for item in dir(database):
            if item in ['qcdb', 'rxn', 'dbse', 'TAGL']:
                pass
            elif item.startswith('__'):
                pass
            else:
                pieces.append(item)

        # form qcdb.Reagent objects from all defined geometries, GEOS
        oHRGT = {}
        for rgt, mol in database.GEOS.items():
            mol.update_geometry()
            try:
                tagl = database.TAGL[rgt]
            except KeyError:
                tagl = None
                print("""Warning: TAGL missing for reagent %s""" % (rgt))
            oHRGT[rgt] = Reagent(name=rgt, mol=mol, tagl=tagl)
        pieces.remove('GEOS')
        self.hrgt = oHRGT

        # form qcdb.Reaction objects from comprehensive reaction list, HRXN
        oHRXN = OrderedDict()
        for rxn in database.HRXN:
            try:
                tagl = database.TAGL[database.dbse + '-' + str(rxn)]
            except KeyError:
                tagl = None
                print("""Warning: TAGL missing for reaction %s""" % (rxn))
            try:
                elst = database.DATA['SAPT ELST ENERGY'][database.dbse + '-' + str(rxn)]
                disp = database.DATA['SAPT DISP ENERGY'][database.dbse + '-' + str(rxn)]
                color = abs(elst) / (abs(elst) + abs(disp))
            except (KeyError, AttributeError):
                color = 'black'
                print("""Warning: DATA['SAPT * ENERGY'] missing for reaction %s""" % (rxn))

            oHRXN[rxn] = Reaction(name=rxn,
                                  dbse=database.dbse,
                                  indx=database.HRXN.index(rxn) + 1,
                                  color=color,
                                  tagl=tagl)
        pieces.remove('HRXN')
        self.hrxn = oHRXN

        # list and align database stoichiometry modes, ACTV* and RXNM*
        oACTV = {}
        for modactv in [item for item in pieces if item.startswith('ACTV')]:
            modrxnm = modactv.replace('ACTV', 'RXNM')
            mode = 'default' if modactv == 'ACTV' else modactv.replace('ACTV_', '')
            try:
                getattr(database, modrxnm)
            except AttributeError:
                modrxnm = 'RXNM'
            oACTV[mode] = [modactv, modrxnm]
        for item in [tmp for tmp in pieces if tmp.startswith('ACTV') or tmp.startswith('RXNM')]:
            pieces.remove(item)

        # populate reaction matrices in qcdb.Reaction objects
        for rxn in database.HRXN:
            dbrxn = database.dbse + '-' + str(rxn)
            for mode, actvrxnm in oACTV.items():
                tdict = OrderedDict()
                for rgt in getattr(database, actvrxnm[0])[dbrxn]:
                    tdict[oHRGT[rgt]] = getattr(database, actvrxnm[1])[dbrxn][rgt]
                oHRXN[rxn].rxnm[mode] = tdict

        # list embedded quantum chem info per rxn, incl. BIND*
        arrsbind = [item for item in pieces if item.startswith('BIND_')]
        if len(arrsbind) == 0:
            if 'BIND' in pieces:
                arrsbind = ['BIND']
            else:
                arrsbind = []
                print("""Warning: No BIND array with reference values.""")
        else:
            for arrbind in arrsbind:
                if getattr(database, arrbind) is database.BIND:
                    break
            else:
                print("""Warning: No BIND_* array assigned to be master BIND.""")

        oBIND = {}
        for arrbind in arrsbind:
            ref = database.dbse + 'REF' if arrbind == 'BIND' else arrbind.replace('BIND_', '')
            methods[ref] = Method(name=ref)
            bases[ref] = BasisSet(name=ref)
            try:
                getattr(database, 'BINDINFO_' + ref)
            except AttributeError:
                arrbindinfo = None
                print("""Warning: No BINDINFO dict with BIND attribution and modelchem for %s.""" % (ref))
            else:
                arrbindinfo = 'BINDINFO_' + ref
            oBIND[ref] = [methods[ref], 'default', bases[ref], arrbind,
                          (getattr(database, arrbind) is database.BIND),
                          arrbindinfo]
        for item in [tmp for tmp in pieces if tmp.startswith('BIND')]:
            pieces.remove(item)

        # populate data with reference values in qcdb.Reaction objects
        for rxn in database.HRXN:
            dbrxn = database.dbse + '-' + str(rxn)
            for ref, info in oBIND.items():
                bindval = getattr(database, info[3])[dbrxn]
                if info[5] is None:
                    methodfeed = info[0]
                    modefeed = info[1]
                    basisfeed = info[2]
                    citationkey = 'anon'
                else:
                    bindinforxn = getattr(database, info[5])[dbrxn]
                    methodfeed = methods[bindinforxn['method'].upper()] if 'method' in bindinforxn else info[0]
                    modefeed = bindinforxn['mode'] if 'mode' in bindinforxn else info[1]
                    basisfeed = bases[bindinforxn['basis'].lower()] if 'basis' in bindinforxn else info[2]
                    citationkey = bindinforxn['citation'].lower() if 'citation' in bindinforxn else 'anon'
                citationfeed = pubs[citationkey]

                if bindval is not None:
                    oHRXN[rxn].data[ref] = ReactionDatum(dbse=database.dbse, rxn=rxn,
                                                         method=methodfeed, mode=modefeed,
                                                         basis=basisfeed, citation=citationfeed,
                                                         value=bindval)
                    # oHRXN[rxn].data[ref] = ReactionDatum(dbse=database.dbse,
                    #                                     rxn=rxn,
                    #                                     method=info[0],
                    #                                     mode=info[1],
                    #                                     basis=info[2],
                    #                                     value=bindval)
                    #                                     #value=getattr(database, info[3])[dbrxn])
                    if info[4]:
                        oHRXN[rxn].benchmark = ref

        # Process subsets
        oSSET = {}
        fsHRXN = frozenset(database.HRXN)
        for sset in pieces:
            if not sset.startswith('AXIS_'):
                try:
                    fssset = frozenset(getattr(database, sset))
                except TypeError:
                    continue
                if fssset.issubset(fsHRXN):
                    oSSET[sset] = getattr(database, sset)
        for item in oSSET.keys():
            pieces.remove(item)
        oSSET['HRXN'] = database.HRXN

        self.sset = OrderedDict()
        self.oss = OrderedDict()  # just in case oss replaces sset someday
        for item in oSSET.keys():
            if item == 'HRXN_SM':
                label = 'small'
            elif item == 'HRXN_LG':
                label = 'large'
            elif item == 'HRXN_EQ':
                label = 'equilibrium'
            elif item == 'HRXN':
                label = 'default'
            elif item.startswith('HRXN_'):
                label = item.replace('HRXN_', '').lower()
            else:
                label = item.lower()

            # subsets may have different ordering from HRXN
            self.sset[label] = OrderedDict()
            for rxn in oSSET[item]:
                self.sset[label][rxn] = oHRXN[rxn]

            # initialize subset objects with light info
            try:
                sstagl = database.TAGL[item]
            except KeyError:
                try:
                    sstagl = database.TAGL[label]
                except KeyError:
                    sstagl = None
                    print("""Warning: TAGL missing for subset %s""" % (label))
            self.oss[label] = Subset(name=label,
                                     hrxn=self.sset[label].keys(),
                                     tagl=sstagl)

        # Process axes
        for axis in [item for item in pieces if item.startswith('AXIS_')]:
            label = axis.replace('AXIS_', '')
            try:
                defn = getattr(database, axis)
            except AttributeError:
                raise ValidationError("""Axis %s not importable.""" % (label))
            axisrxns = frozenset(defn.keys())
            attached = False
            for ss, rxns in self.sset.items():
                if frozenset(rxns).issubset(axisrxns):
                    ordered_floats = []
                    for rx in self.oss[ss].hrxn:
                        ordered_floats.append(defn[rx])
                    self.oss[ss].axis[label] = ordered_floats
                    attached = True
            if not attached:
                print("""Warning: AXIS %s not affiliated with a subset""" % (label))
            pieces.remove(axis)

        print("""WrappedDatabase %s: Unparsed attributes""" % (self.dbse), pieces)

    def __str__(self):
        text = ''
        text += """  ==> %s WrappedDatabase <==\n\n""" % (self.dbse)
        text += """  Reagents:             %s\n""" % (self.hrgt.keys())
        text += """  Reactions:            %s\n""" % (self.hrxn.keys())
        text += """  Subsets:              %s\n""" % (self.sset.keys())
        text += """  Reference:            %s\n""" % (self.benchmark())
        text += """\n"""
        return text

    def add_ReactionDatum(self, dbse, rxn, method, mode, basis, value, units='kcal/mol', citation=None, comment=None,
                          overwrite=False):
        """Add a new quantum chemical value to *rxn* by creating a
        qcdb.ReactionDatum from same arguments as that class's
        object-less constructor. *rxn* may be actual Reaction.name
        or Reaction.indx.

        """
        if (self.dbse == dbse):
            if rxn in self.hrxn:
                rxnname = rxn  # rxn is proper reaction name
            else:
                try:
                    if (rxn + 1 > 0) and (rxn == self.hrxn.items()[rxn - 1][1].indx):
                        rxnname = self.hrxn.items()[rxn - 1][1].name  # rxn is reaction index (maybe dangerous?)
                except (TypeError, IndexError):
                    raise ValidationError(
                        """Inconsistent to add ReactionDatum for %s to database %s with reactions %s.""" %
                        (dbse + '-' + str(rxn), self.dbse, self.hrxn.keys()))
            label = '-'.join([method, mode, basis])
            if overwrite or (label not in self.hrxn[rxnname].data):
                self.hrxn[rxnname].data[label] = ReactionDatum.library_modelchem(dbse=dbse, rxn=rxnname,
                                                                                 method=method, mode=mode, basis=basis,
                                                                                 value=value, units=units,
                                                                                 comment=comment, citation=citation)
            else:
                raise ValidationError("""ReactionDatum %s already present in Database.""" % (label))
        else:
            raise ValidationError("""Inconsistent to add ReactionDatum for %s to database %s.""" %
                                  (dbse + '-' + str(rxn), self.dbse))

    def add_Subset(self, name, func):
        """Define a new subset labeled *name* by providing a function
        *func* that filters *self.hrxn*.

        """
        sname = name.lower().split('\n')
        label = sname.pop(0)
        tagl = sname[0].strip() if sname else None
        try:
            filtered = func(self)
            lsslist = [rxn for rxn in self.sset['default'].keys() if rxn in filtered]
        except TypeError as e:
            raise ValidationError("""Function %s did not return list: %s.""" % (func.__name__, str(e)))
        if len(lsslist) == 0:
            print("""WrappedDatabase %s: Subset %s NOT formed: empty""" % (self.dbse, label))
            return

        self.sset[label] = OrderedDict()
        for rxn in lsslist:
            self.sset[label][rxn] = self.hrxn[rxn]
        self.oss[label] = Subset(name=label,
                                 hrxn=self.sset[label].keys(),
                                 tagl=tagl)
        print("""WrappedDatabase %s: Subset %s formed: %d""" % (self.dbse, label, len(self.sset[label].keys())))

    def compute_errors(self, modelchem, benchmark='default', sset='default', failoninc=True, verbose=False):
        """For full database or subset *sset*, computes raw reaction
        errors between *modelchem* and *benchmark* model chemistries.
        Returns error if model chemistries are missing for any reaction in
        subset unless *failoninc* set to False, whereupon returns partial.
        Returns dictionary of reaction labels and error forms.

        """
        if isinstance(sset, basestring):
            # sset is normal subset name 'MX' corresponding to HRXN_MX or MX array in database module
            try:
                lsset = self.sset[sset.lower()]
            except KeyError as e:
                # raise ValidationError("""Subset named %s not available""" % (str(e)))
                lsset = OrderedDict()
        else:
            if callable(sset):
                # sset is function that will generate subset of HRXN from sset(self)
                lsslist = [rxn for rxn in self.sset['default'].keys() if rxn in sset(self)]
            else:
                # sset is array containing reactions
                lsslist = [rxn for rxn in self.sset['default'].keys() if rxn in sset]
            # assemble dict of qcdb.Reaction objects from array of reaction names
            lsset = OrderedDict()
            for rxn in lsslist:
                lsset[rxn] = self.hrxn[rxn]

#        cureinfo = self.get_pec_weightinfo()
        err = {}
        for rxn, oRxn in lsset.items():
            lbench = oRxn.benchmark if benchmark == 'default' else benchmark
            try:
                mcLesser = oRxn.data[modelchem].value
            except KeyError as e:
                if failoninc:
                    raise ValidationError("""Reaction %s missing datum %s.""" % (str(rxn), str(e)))
                else:
                    continue
            try:
                mcGreater = oRxn.data[lbench].value
            except KeyError as e:
                if lbench == 'ZEROS':
                    pass
                else:
                    print("""Reaction %s missing benchmark""" % (str(rxn)))
                    continue
            # handle particulars of PEC error measures
#            rxncureinfo = cureinfo[rxn]
#            try:
#                mcGreaterCrvmin = self.hrxn[rxncureinfo['eq']].data[lbench].value
#            except KeyError as e:
#                print """Reaction %s missing benchmark""" % (str(eqrxn))

#            cure_denom = cure_weight(refrxn=mcGreater, refeq=mcGreaterCrvmin, rrat=rxncureinfo['Rrat'])
#            balanced_mask, balwt = balanced_error(refrxn=mcGreater, refeq=mcGreaterCrvmin, rrat=rxncureinfo['Rrat'])

            if lbench == 'ZEROS':
                err[rxn] = [mcLesser,
                            0.0, 0.0, 0.0, 1.0]  # FAKE
            else:
                err[rxn] = [mcLesser - mcGreater,
                        (mcLesser - mcGreater) / abs(mcGreater),
                        (mcLesser - mcGreater) / abs(mcGreater),  # FAKE
                        (mcLesser - mcGreater) / abs(mcGreater),  # FKAE
                        1.0  # FAKE
                        ]
#                        (mcLesser - mcGreater) / abs(cure_denom),
#                        (mcLesser - mcGreater) * balanced_mask / abs(mcGreaterCrvmin),
#                        balwt]
            if verbose:
                print("""p = %8.4f, pe = %8.3f%%, pbe = %8.3f%% pce = %8.3f%% reaction %s.""" %
                 (err[rxn][0], 100 * err[rxn][1], 100 * err[rxn][3], 100 * err[rxn][2], str(rxn)))
        return err

    def compute_statistics(self, modelchem, benchmark='default', sset='default',
                           failoninc=True, verbose=False, returnindiv=False):
        """For full database or subset *sset*, computes many error
        statistics between single *modelchem* and *benchmark* model
        chemistries. Returns error if model chemistries are missing
        for any reaction in subset unless *failoninc* set to False,
        whereupon returns partial statistics. Returns dictionary of
        statistics labels and values.

        """
        err = self.compute_errors(modelchem, benchmark=benchmark, sset=sset, failoninc=failoninc, verbose=verbose)
        if len(err) == 0:
            error = initialize_errors()
            if verbose:
                print("""Warning: nothing to compute.""")
        else:
            Nrxn = float(len(err))
            error = OrderedDict()
            # linear (absolute) error
            linear = [val[0] for val in err.values()]
            error['pexe'] = max(linear)
            error['nexe'] = min(linear)
            error['maxe'] = max(linear, key=lambda x: abs(x))
            error['mine'] = min(linear, key=lambda x: abs(x))
            error['me'] = sum(linear) / Nrxn
            error['mae'] = sum(map(abs, linear)) / Nrxn
            error['rmse'] = math.sqrt(sum(map(lambda x: x ** 2, linear)) / Nrxn)
            error['stde'] = math.sqrt((sum(map(lambda x: x ** 2, linear)) - (sum(linear) ** 2) / Nrxn) / Nrxn)
            # fractional (relative) error
            relative = [val[1] for val in err.values()]
            error['pexpe'] = max(relative)
            error['nexpe'] = min(relative)
            error['maxpe'] = max(relative, key=lambda x: abs(x))
            error['minpe'] = min(relative, key=lambda x: abs(x))
            error['mpe'] = sum(relative) / Nrxn
            error['mape'] = sum(map(abs, relative)) / Nrxn
            error['rmspe'] = math.sqrt(sum(map(lambda x: x ** 2, relative)) / Nrxn)
            error['stdpe'] = math.sqrt((sum(map(lambda x: x ** 2, relative)) - (sum(relative) ** 2) / Nrxn) / Nrxn)
            # balanced (relative) error
            balanced = [val[3] for val in err.values()]
            balwt = sum([val[4] for val in err.values()])  # get the wt fn. highly irregular TODO
            error['pexpbe'] = max(balanced)
            error['nexpbe'] = min(balanced)
            error['maxpbe'] = max(balanced, key=lambda x: abs(x))
            error['minpbe'] = min(balanced, key=lambda x: abs(x))
            error['mpbe'] = sum(balanced) / balwt #Nrxn
            error['mapbe'] = sum(map(abs, balanced)) / balwt #Nrxn
            error['rmspbe'] = math.sqrt(sum(map(lambda x: x ** 2, balanced)) / balwt) #Nrxn)
            error['stdpbe'] = None  # get math domain errors w/wt in denom math.sqrt((sum(map(lambda x: x ** 2, balanced)) - (sum(balanced) ** 2) / balwt) / balwt) #/ Nrxn) / Nrxn)
            # capped (relative) error
            capped = [val[2] for val in err.values()]
            error['pexpce'] = max(capped)
            error['nexpce'] = min(capped)
            error['maxpce'] = max(capped, key=lambda x: abs(x))
            error['minpce'] = min(capped, key=lambda x: abs(x))
            error['mpce'] = sum(capped) / Nrxn
            error['mapce'] = sum(map(abs, capped)) / Nrxn
            error['rmspce'] = math.sqrt(sum(map(lambda x: x ** 2, capped)) / Nrxn)
            error['stdpce'] = math.sqrt((sum(map(lambda x: x ** 2, capped)) - (sum(capped) ** 2) / Nrxn) / Nrxn)
            if verbose:
                print("""%d systems in %s for %s vs. %s, subset %s.\n%s""" %
                      (len(err), self.dbse, modelchem, benchmark, sset, format_errors(error, mode=2)))
        if returnindiv:
            return error, err
        else:
            return error

    def load_qcdata(self, modname, funcname, pythonpath=None, failoninc=True):
        """Loads qcdb.ReactionDatums from module *modname* function
        *funcname*. Module search path can be prepended with *pythonpath*.

        """
        if pythonpath is not None:
            sys.path.insert(1, pythonpath)
        else:
            sys.path.append(os.path.dirname(__file__) + '/../data')
        try:
            datamodule = __import__(modname)
        except ImportError:
            if not failoninc:
                print("""%s data unavailable for database %s.\n""" % (modname, self.dbse))
                return
            else:
                print("""\nPython module for database data %s failed to load\n\n""" % (modname))
                print("""\nSearch path that was tried:\n""")
                print(', '.join(map(str, sys.path)))
                raise ValidationError("""Python module loading problem for database data """ + str(modname))
        try:
            getattr(datamodule, funcname)(self)
        except AttributeError:
            if not failoninc:
                print("""%s %s data unavailable for database %s.\n""" % (modname, funcname, self.dbse))
                return
            else:
                raise ValidationError("Python module missing function %s for loading data " % (str(funcname)))

        print("""WrappedDatabase %s: %s %s results loaded""" % (self.dbse, modname, funcname))

    def load_qcdata_byproject(self, project, pythonpath=None):
        """Loads qcdb.ReactionDatums from standard location for *project*
        :module dbse_project and function load_project. Module search path
        can be prepended with *pythonpath*.

        """
        mod = self.dbse + '_' + project
        func = 'load_' + project
        self.load_qcdata(modname=mod, funcname=func, pythonpath=pythonpath)

    def load_qcdata_hrxn_byproject(self, project, path=None):
        """"""
        if path is None:
            path = os.path.dirname(__file__) + '/../data'
        pklfile = os.path.abspath(path) + os.sep + self.dbse + '_hrxn_' + project + '.pickle'
        if not os.path.isfile(pklfile):
            raise ValidationError(
                "Reactions pickle file for loading database data from file %s does not exist" % (pklfile))

        with open(pklfile, 'rb') as handle:
            hrxns = pickle.load(handle)
        # no error checking for speed
        for rxn, data in hrxns.items():
            self.hrxn[rxn].data.update(data)

    def load_qcdata_hdf5_trusted(self, project, path=None):
        """Loads qcdb.ReactionDatums from HDF5 file at path/dbse_project.h5 .
        If path not given, looks in qcdb/data. This file is written by
        reap-DB and so has been largely validated.

        """
        if path is None:
            path = os.path.dirname(__file__) + '/../data'
        hdf5file = os.path.abspath(path) + os.sep + self.dbse + '_' + project + '.h5'
        if not os.path.isfile(hdf5file):
            raise ValidationError("HDF5 file for loading database data from file %s does not exist" % (hdf5file))
        try:
            import pandas as pd
        except ImportError:
            raise ValidationError("Pandas data managment module must be available for import")

        try:
            next(self.hrxn.iterkeys()) + 1
        except TypeError:
            intrxn = False
        else:
            intrxn = True

        with pd.get_store(hdf5file) as handle:
            for mc in handle['pdie'].keys():
                lmc = mc.split('-')  # TODO could be done better
                method = lmc[0]
                bsse = '_'.join(lmc[1:-1])
                basis = lmc[-1]

                df = handle['pdie'][mc]
                for dbrxn in df.index[df.notnull()].values:
                    [dbse, rxn] = dbrxn.split('-', 1)
                    if intrxn:
                        rxn = int(rxn)
                    self.hrxn[rxn].data[mc] = ReactionDatum.library_modelchem(dbse=dbse, rxn=rxn,
                                                                              method=method, mode=bsse, basis=basis,
                                                                              value=df[dbrxn])

    def integer_reactions(self):
        """Returns boolean of whether reaction names need to be cast to integer"""
        try:
            next(self.hrxn.iterkeys()) + 1
        except TypeError:
            return False
        else:
            return True

    @staticmethod
    def load_pickled(dbname, path=None):
        """

        """
        if path is None:
            path = os.path.dirname(__file__) + '/../data'
        picklefile = psiutil.findfile_ignorecase(dbname,
                                                 pre=os.path.abspath(path) + os.sep, post='_WDb.pickle')
        if not picklefile:
            raise ValidationError("Pickle file for loading database data from file %s does not exist" % (
                os.path.abspath(path) + os.sep + dbname + '.pickle'))
        # with open('/var/www/html/bfdb_devel/bfdb/scratch/ASDFlogfile.txt', 'a') as handle:
        #    handle.write('<!-- PICKLE %s\n' % (picklefile))
        with open(picklefile, 'rb') as handle:
            instance = pickle.load(handle)
        return instance

    def available_modelchems(self, union=True):
        """Returns all the labels of model chemistries that have been
        loaded. Either all modelchems that have data for any reaction if
        *union* is True or all modelchems that have data for all reactions
        if *union* is False.

        """
        mcs = [set(v.data) for v in self.hrxn.itervalues()]
        if union:
            return sorted(set.union(*mcs))
        else:
            return sorted(set.intersection(*mcs))

    def benchmark(self):
        """Returns the model chemistry label for the database's benchmark."""
        bm = None
        rxns = self.hrxn.itervalues()
        while bm is None:
            try:
                bm = next(rxns).benchmark
            except StopIteration:
                break
        return bm
        # return next(self.hrxn.itervalues()).benchmark
        # TODO all rxns have same bench in db module so all have same here in obj
        #   but the way things stored in Reactions, this doesn't have to be so

    def load_subsets(self, modname='subsetgenerator', pythonpath=None):
        """Loads subsets from all functions in module *modname*.

        """
        if pythonpath is not None:
            sys.path.insert(1, pythonpath)
        else:
            sys.path.append(os.path.dirname(__file__))
        try:
            ssmod = __import__(modname)
        except ImportError:
            print("""\nPython module for database data %s failed to load\n\n""" % (modname))
            print("""\nSearch path that was tried:\n""")
            print(', '.join(map(str, sys.path)))
            raise ValidationError("Python module loading problem for database subset generator " + str(modname))

        for func in dir(ssmod):
            if callable(getattr(ssmod, func)):
                self.add_Subset(getattr(ssmod, func).__doc__, getattr(ssmod, func))

        print("""WrappedDatabase %s: Defined subsets loaded""" % (self.dbse))

    def get_pec_weightinfo(self):
        """

        """
        def closest(u, options):
            return max(options, key=lambda v: len(os.path.commonprefix([u, v])))

        dbdat = {}
        oss = self.oss['default']
        eqrxns = [rxn for rxn, rr in zip(oss.hrxn, oss.axis['Rrat']) if rr == 1.0]
        for rxnix, rxn in enumerate(oss.hrxn):
            dbdat[rxn] = {'eq': closest(rxn, eqrxns),
                          'Rrat': oss.axis['Rrat'][rxnix]}
        return dbdat

    # def table_simple1(self, mtd, bas, opt=['CP'], err=['mae'], benchmark='default', failoninc=True,
    #                   plotpath='analysis/flats/flat_', theme='smmerge'):
    #     rowplan = ['bas', 'mtd']
    #     columnplan = [
    #         ['l', r"""Method \& Basis Set""", '', textables.label, {}],
    #         ['d', r'S22', 'HB', textables.val, {'sset': 'hb'}],
    #         ['d', r'S22', 'MX', textables.val, {'sset': 'mx'}],
    #         ['d', r'S22', 'DD', textables.val, {'sset': 'dd'}],
    #         ['d', r'S22', 'TT', textables.val, {'sset': 'default'}],
    #     ]
    #
    # def table_simple2(self, mtd, bas, opt=['CP'], err=['mae'], benchmark='default', failoninc=True,
    #                   plotpath='analysis/flats/flat_', theme='smmerge'):
    #     rowplan = ['bas', 'mtd']
    #     columnplan = [
    #         ['l', r"""Method \& Basis Set""", '', textables.label, {}],
    #         ['d', r'MAE', 'HB', textables.val, {'sset': 'hb'}],
    #         ['d', r'MAE', 'MX', textables.val, {'sset': 'mx'}],
    #         ['d', r'MAE', 'DD', textables.val, {'sset': 'dd'}],
    #         ['d', r'MAE', 'TT', textables.val, {'sset': 'default'}],
    #         ['d', r'MA\%E', 'HB', textables.val, {'sset': 'hb', 'err': 'mape'}],
    #         ['d', r'MA\%E', 'MX', textables.val, {'sset': 'mx', 'err': 'mape'}],
    #         ['d', r'MA\%E', 'DD', textables.val, {'sset': 'dd', 'err': 'mape'}],
    #         ['d', r'MA\%E', 'TT', textables.val, {'sset': 'default', 'err': 'mape'}],
    #         ['d', r'maxE', 'TT ', textables.val, {'sset': 'default', 'err': 'maxe'}],
    #         ['d', r'min\%E', ' TT', textables.val, {'sset': 'default', 'err': 'minpe'}],
    #         ['d', r'rmsE', 'TT ', textables.val, {'sset': 'default', 'err': 'rmse'}],
    #         ['d', r'devE', ' TT', textables.val, {'sset': 'default', 'err': 'stde'}],
    #     ]
    #
    # def table_simple3(self, mtd, bas, opt=['CP'], err=['mae'], benchmark='default', failoninc=True,
    #                   plotpath='analysis/flats/flat_', theme='smmerge'):
    #     rowplan = ['err', 'bas', 'mtd']
    #     columnplan = [
    #         ['l', r"""Method \& Basis Set""", '', textables.label, {}],
    #         ['d', r'MAE', 'HB', textables.val, {'sset': 'hb'}],
    #         ['d', r'MAE', 'MX', textables.val, {'sset': 'mx'}],
    #         ['d', r'MAE', 'DD', textables.val, {'sset': 'dd'}],
    #         ['d', r'MAE', 'TT', textables.val, {'sset': 'default'}],
    #     ]
    #
    # def table_simple4(self, mtd, bas, opt=['CP'], err=['mae'], benchmark='default', failoninc=True,
    #                   plotpath='analysis/flats/flat_', theme='smmerge'):
    #     plotpath = 'autogen'  # TODO handle better
    #     rowplan = ['bas', 'mtd']
    #     columnplan = [
    #         ['l', r"""Method \& Basis Set""", '', textables.label, {}],
    #         ['d', r'S22', 'HB', textables.val, {'sset': 'hb'}],
    #         ['d', r'S22', 'MX', textables.val, {'sset': 'mx'}],
    #         ['d', r'S22', 'DD', textables.val, {'sset': 'dd'}],
    #         ['d', r'S22', 'TT', textables.val, {'sset': 'default'}],
    #         # ['l', r"""Error Distribution\footnotemark[1]""", r"""\includegraphics[width=6.67cm,height=3.5mm]{%s%s.pdf}""" % (plotpath, 'blank'), textables.graphics, {}],
    #         ['l', r"""Error Distribution\footnotemark[1]""", r"""""", textables.graphics, {}],
    #     ]


class Database(object):
    """Collection for handling single or multiple qcdb.WrappedDatabase objects.
    Particularly, unifying modelchem and subset names that when inconsistent
    across component databases. Also, defining statistics across databases.

    >>> asdf = qcdb.Database(['s22', 'Nbc10', 'hbc6', 'HSG'], 'DB4')
    >>> qwer = qcdb.Database('s22')
    """

    def __init__(self, dbnamelist, dbse=None, pythonpath=None, loadfrompickle=False, path=None):
        #: internal name of database collection
        #:
        #: >>> print asdf.dbse
        #: 'DB4'
        self.dbse = None

        #: ordered component Database objects
        #:
        #: >>> print asdf.dbdict
        #: XXXX
        self.dbdict = OrderedDict()

        #: subset assembly pattern
        #:
        #: >>> print asdf.sset.keys()
        #: XXXX
        self.sset = OrderedDict()

        #: assembly pattern for transspecies modelchems
        #:
        #: >>> print asdf.mcs.keys()
        #: XXXX
        self.mcs = {}

        self.benchmark = None

        # slight validation, repackaging into dbnamelist
        if isinstance(dbnamelist, basestring):
            dbnamelist = [dbnamelist]
        elif all(isinstance(item, basestring) for item in dbnamelist):
            pass
        else:
            raise ValidationError('Database::constructor: Inappropriate configuration of constructor arguments')

        # load databases
        for db in dbnamelist:
            if loadfrompickle:
                tmp = WrappedDatabase.load_pickled(db, path=path)
            else:
                tmp = WrappedDatabase(db, pythonpath=pythonpath)
            self.dbdict[tmp.dbse] = tmp

        # slurp up the obvious overlaps
        consolidated_bench = [odb.benchmark() for odb in self.dbdict.values()]
        if len(set(consolidated_bench)) == 1:
            self.benchmark = consolidated_bench[0]
        else:
            self.benchmark = ''.join(consolidated_bench)
        self.mcs[self.benchmark] = consolidated_bench

        # methods[ref] = Method(name=ref)
        # bases[ref] = BasisSet(name=ref)

        self.mcs['default'] = consolidated_bench
        # self.mcs['default'] = [odb.benchmark() for odb in self.dbdict.values()]
        self._intersect_subsets()
        self._intersect_modelchems()

        # complex subsets
        self.load_subsets()

        # collection name
        self.dbse = ''.join(self.dbdict.keys()) if dbse is None else dbse

        # merge Reaction-s
        self.hrxn = OrderedDict()
        for db, odb in self.dbdict.items():
            for rxn, orxn in odb.hrxn.items():
                self.hrxn[orxn.dbrxn] = orxn

        # merge Reagent-s
        self.hrgt = OrderedDict()
        for db, odb in self.dbdict.items():
            for rgt, orgt in odb.hrgt.items():
                self.hrgt[orgt.name] = orgt

        print("""Database %s: %s""" % (self.dbse, ', '.join(self.dbdict.keys())))

    def __str__(self):
        text = ''
        text += """  ===> %s Database <===\n\n""" % (self.dbse)
        # text += """  Reagents:             %s\n""" % (self.hrgt.keys())
        # text += """  Reactions:            %s\n""" % (self.hrxn.keys())
        text += """  Subsets:              %s\n""" % (self.sset.keys())
        # text += """  Reference:            %s\n""" % ('default: ' + ' + '.join(self.mcs['default']))
        try:
            text += """  Reference:            %s\n""" % (self.benchmark + ': ' + ' + '.join(self.mcs[self.benchmark]))
        except TypeError:
            text += """  Reference:            %s\n""" % ('UNDEFINED')
        text += """  Model Chemistries:    %s\n""" % (
            ', '.join(sorted([mc for mc in self.mcs.keys() if mc is not None])))
        text += """\n"""
        for db in self.dbdict.keys():
            text += self.dbdict[db].__str__()
        return text

    # def benchmark(self):
    #     """Returns the model chemistry label for the database's benchmark."""
    #     return self.benchmark  #TODO not sure if right way to go about this self.mcs['default']

    def fancy_mcs(self, latex=False):
        """

        """
        fmcs = {}
        for mc in self.mcs.keys():
            try:
                mtd, mod, bas = mc.split('-')
            except ValueError:
                fmcs[mc] = mc
            else:
                if latex:
                    tmp = """%s/%s, %s""" % \
                          (methods[mtd].latex, bases[bas].latex, mod.replace('_', '\\_'))
                    fmcs[mc] = """%45s""" % (tmp)
                else:
                    fmcs[mc] = """%20s / %-20s, %s""" % \
                               (methods[mtd].fullname, bases[bas].fullname, mod)
        return fmcs

    # def fancy_mcs_nested(self):
    #    """

    #    """
    #    fmcs = defaultdict(lambda: defaultdict(dict))
    #    for mc in self.mcs.keys():
    #        try:
    #            mtd, mod, bas = mc.split('-')
    #        except ValueError:
    #            fmcs['All']['All'][mc] = mc
    #            fmcs['Method']['Others'][mc] = mc
    #            fmcs['Options']['Others'][mc] = mc
    #            fmcs['Basis Treatment']['Others'][mc] = mc
    #        else:
    #            fancyrepr = """%20s / %-20s %s""" % (methods[mtd].latex, bases[bas].latex, mod)
    #            fmcs['All']['All'][mc] = fancyrepr
    #            fmcs['Method'][methods[mtd].latex][mc] = fancyrepr
    #            fmcs['Options'][mod][mc] = fancyrepr
    #            fmcs['Basis Treatment'][bases[bas].latex][mc] = fancyrepr
    #    return fmcs

    def integer_reactions(self):
        """Returns boolean of whether reaction names need to be cast to integer"""
        return {db: odb.integer_reactions() for db, odb in self.dbdict.items()}

    def load_qcdata_byproject(self, project, pythonpath=None):
        """For each component database, loads qcdb.ReactionDatums from
        standard location for *project* :module dbse_project and function
        load_project. Module search path can be prepended with *pythonpath*.

        """
        for db, odb in self.dbdict.items():
            odb.load_qcdata_byproject(project, pythonpath=pythonpath)
        self._intersect_modelchems()

    def load_qcdata_hdf5_trusted(self, project, path=None):
        """For each component database, loads qcdb.ReactionDatums from
        HDF5 file at path/dbse_project.h5 . If path not given, looks in
        qcdb/data. This file is written by reap-DB and so has been largely
        validated.

        """
        for db, odb in self.dbdict.items():
            odb.load_qcdata_hdf5_trusted(project, path=path)
        self._intersect_modelchems()

    def load_qcdata_hrxn_byproject(self, project, path=None):
        for db, odb in self.dbdict.items():
            odb.load_qcdata_hrxn_byproject(project, path=path)
        self._intersect_modelchems()

    def available_projects(self, path=None):
        """"""
        import glob

        if path is None:
            path = os.path.dirname(__file__) + '/../data'

        projects = []
        for pjfn in glob.glob(path + '/*_hrxn_*.pickle'):
            pj = pjfn[:-7].split('_')[-1]
            projects.append(pj)

        complete_projects = []
        for pj in set(projects):
            if all([os.path.isfile(path + '/' + db + '_hrxn_' + pj + '.pickle') for db in self.dbdict.keys()]):
                complete_projects.append(pj)

        return complete_projects

    def load_subsets(self, modname='subsetgenerator', pythonpath=None):
        """For each component database, loads subsets from all functions
        in module *modname*. Default *modname* usues standard generators.

        """
        for db, odb in self.dbdict.items():
            odb.load_subsets(modname=modname, pythonpath=pythonpath)
        self._intersect_subsets()

    def add_Subset(self, name, func):
        """Define a new subset labeled *name* by providing a database
        *func* whose keys are the keys of dbdict and whose values are a
        function that filters each WrappedDatabase's *self.hrxn*.

        """
        label = name.lower()
        merged = []
        for db, odb in self.dbdict.items():
            if callable(func[db]):
                ssfunc = func[db]
            else:
                ssfunc = lambda x: func[db]
            odb.add_Subset(name=name, func=ssfunc)
            if name in odb.sset:
                merged.append(name)
            else:
                merged.append(None)
        if any(merged):
            self.sset[label] = merged
            print("""Database %s: Subset %s formed: %s""" % (self.dbse, label, self.sset[label]))
        else:
            print("""Database %s: Subset %s NOT formed: empty""" % (self.dbse, label))

    def add_Subset_union(self, name, sslist):
        """
        Define a new subset labeled *name* (note that there's nothing to
        prevent overwriting an existing subset name) from the union of
        existing named subsets in *sslist*.

        """
        funcdb = {}
        for db, odb in self.dbdict.items():
            dbix = self.dbdict.keys().index(db)
            overlapping_dbrxns = []
            for ss in sslist:
                lss = self.sset[ss][dbix]
                if lss is not None:
                    overlapping_dbrxns.append(self.dbdict[db].sset[lss].keys())
            rxnlist = set().union(*overlapping_dbrxns)
            funcdb[db] = rxnlist
        self.add_Subset(name, funcdb)

    def add_sampled_Subset(self, sset='default', number_of_samples=1, sample_size=5, prefix='rand'):
        """Generate and register *number_of_samples* new subsets of size
        *sample_size* and name built from *prefix*. Reactions chosen from *sset*.

        """
        import random

        intrxn = self.integer_reactions()
        rxns = self.get_hrxn(sset=sset).keys()

        def random_sample(ssname):
            """Generate and register a single new subset of size *sample_size* and
            name *ssname*.

            """
            sample = {db: [] for db in self.dbdict.keys()}
            for dbrxn in random.sample(rxns, sample_size):
                db, rxn = dbrxn.split('-', 1)
                typed_rxn = int(rxn) if intrxn[db] else rxn
                sample[db].append(typed_rxn)
            self.add_Subset(ssname, sample)

        for sidx in range(number_of_samples):
            if number_of_samples == 1:
                ssname = prefix
            else:
                ssname = prefix + '_' + str(sidx)
            random_sample(ssname)

    def promote_Subset(self, name=None):
        """Examine component databases and elevate subset *name* not necessarily
        present for all component databases to a subset for the *self*. When *name*
        is None, promotes all subsets found for component databases. Also promotes
        entirety of each component database as a subset with name of component
        database dbse in lowercase.

        """
        if name is None:
            sss = [set(odb.sset.keys()) for db, odb in self.dbdict.items()]
            new = sorted(set.union(*sss))
        else:
            new = [name]
        for ss in new:
            if ss not in self.sset:
                self.sset[ss] = [ss if ss in odb.sset else None for db, odb in self.dbdict.items()]
                print("""Database %s: Subset %s promoted: %s""" % (self.dbse, ss, self.sset[ss]))
        if name is None and len(self.dbdict) > 1:
            for db, odb in self.dbdict.items():
                dbix = self.dbdict.keys().index(db)
                ss = odb.dbse.lower()
                if ss not in self.sset:
                    self.sset[ss] = ['default' if ix == dbix else None for ix in range(len(self.dbdict))]
                    print("""Database %s: Subset %s promoted: %s""" % (self.dbse, ss, self.sset[ss]))

    def _intersect_subsets(self):
        """Examine component database subsets and collect common names as
        Database subset.

        """
        sss = [set(odb.sset.keys()) for db, odb in self.dbdict.items()]
        new = sorted(set.intersection(*sss))
        for ss in new:
            self.sset[ss] = [ss] * len(self.dbdict.keys())

    def _intersect_modelchems(self):
        """Examine component database qcdata and collect common names as
        Database modelchem.

        """
        mcs = [set(odb.available_modelchems()) for odb in self.dbdict.itervalues()]
        new = sorted(set.intersection(*mcs))
        for mc in new:
            self.mcs[mc] = [mc] * len(self.dbdict.keys())

    # def reaction_generator(self):
    #    """

    #    """
    #    for db, odb in self.dbdict.items():
    #        for rxn, orxn in odb.hrxn.items():
    #            yield orxn

    def compute_statistics(self, modelchem, benchmark='default', sset='default',
                           failoninc=True, verbose=False, returnindiv=False):
        """Computes summary statistics and, if *returnindiv* True,
        individual errors for single model chemistry *modelchem* versus
        *benchmark* over subset *sset* over all component databases.
        Particularly, imposes cross-database definitions for sset and
        modelchem.
        #Returns error if model chemistries are missing
        #for any reaction in subset unless *failoninc* set to False,
        #whereupon returns partial statistics. Returns dictionary of
        #statistics labels and values.

        """
        errors = OrderedDict()
        indiv = OrderedDict()
        actvdb = []
        for db, odb in self.dbdict.items():
            dbix = self.dbdict.keys().index(db)
            if self.sset[sset][dbix] is None:
                errors[db], indiv[db] = (None, None)
            else:
                errors[db], indiv[db] = odb.compute_statistics(self.mcs[modelchem][dbix],
                                                               sset=self.sset[sset][dbix],
                                                               benchmark='ZEROS' if benchmark == 'ZEROS' else self.mcs[benchmark][dbix],
                                                               failoninc=failoninc, verbose=verbose, returnindiv=True)
                actvdb.append(errors[db])
        errors[self.dbse] = average_errors(*actvdb)

        if returnindiv:
            return errors, indiv
        else:
            return errors

    def analyze_modelchems(self, modelchem, benchmark='default', failoninc=True, verbose=False):
        """For each component database, compute and print nicely formatted
        summary error statistics for each model chemistry in array
        *modelchem* versus *benchmark* for all available subsets.

        """
        # compute errors
        errors = {}
        for mc in modelchem:
            errors[mc] = {}
            for ss in self.sset.keys():
                errors[mc][ss] = self.compute_statistics(mc, benchmark=benchmark, sset=ss,
                                                         failoninc=failoninc, verbose=verbose, returnindiv=False)
        # present errors
        pre, suf, mid = string_contrast(modelchem)
        text = """\n  ==> %s %s[]%s Errors <==\n""" % (self.dbse, pre, suf)
        text += """%20s    %44s""" % ('', '==> ' + self.dbse + ' <==')
        for db, odb in self.dbdict.items():
            text += """%44s""" % ('=> ' + odb.dbse + ' <=')
        text += '\n'

        collabel = """      {:5}   {:4}   {:6} {:6}    {:6}""".format(
                'ME', 'STDE', 'MAE', 'MA%E', 'MA%BE')

        text += """{:20}    """.format('') + collabel
        for db in self.dbdict.keys():
            text += collabel
        text += '\n'

        text += """{:20}    {}""".format('', '=' * 44)
        ul = False
        for db in self.dbdict.keys():
            text += """{}""".format('_' * 44 if ul else ' ' * 44)
            ul = not ul
        text += '\n'

        for ss in self.sset.keys():
            text += """   => %s <=\n""" % (ss)
            for mc in modelchem:
                perr = errors[mc][ss]
                text += """%20s    %44s""" % (mid[modelchem.index(mc)],
                                              format_errors(perr[self.dbse]))
                for db in self.dbdict.keys():
                    text += """%44s""" % ('' if perr[db] is None else format_errors(perr[db]))
                text += '\n'
        print(text)

    def plot_bars(self, modelchem, benchmark='default', sset=['default', 'hb', 'mx', 'dd'],
                  failoninc=True, verbose=False, view=True,
                  saveas=None, relpath=False, graphicsformat=['pdf']):
        """Prepares 'grey bars' diagram for each model chemistry in array
        *modelchem* versus *benchmark* over all component databases. A wide bar
        is plotted with three smaller bars, corresponding to the 'mae'
        summary statistic of the four subsets in *sset*.

        *saveas* conveys directory ('/') and/or filename for saving the
        resulting plot. File extension is not accessible, but *graphicsformat*
        array requests among 'png', 'pdf', and 'eps' formats. *relpath*
        forces paths to saved files to be relative to current directory,
        rather than absolute paths for returned code and file dictionary.

        Prepares bars diagram instructions and either executes them if
        matplotlib available (Canopy or Anaconda) or prints them. Returns a
        dictionary of all saved plot filenames.

        >>> asdf.plot_bars(['MP2-CP-adz', 'MP2-CP-adtz'], sset=['tt-5min', 'hb-5min', 'mx-5min', 'dd-5min'])
        """
        # compute errors
        errors = {}
        for mc in modelchem:
            if mc is not None:
                errors[mc] = {}
                for ss in sset:
                    errors[mc][ss] = self.compute_statistics(mc, benchmark=benchmark, sset=ss,
                                                             failoninc=failoninc, verbose=verbose, returnindiv=False)
        # repackage
        pre, suf, mid = string_contrast(modelchem)
        dbdat = []
        for mc in modelchem:
            if mc is None:
                dbdat.append(None)
            else:
                dbdat.append({'mc': mid[modelchem.index(mc)],
                              'data': [errors[mc][ss][self.dbse]['mae'] for ss in sset]})
        title = self.dbse + ' ' + pre + '[]' + suf + ' ' + ','.join(sset)
        # generate matplotlib instructions and call or print
        try:
            from . import mpl
            import matplotlib.pyplot as plt
        except ImportError:
            # if not running from Canopy, print line to execute from Canopy
            print("""filedict = mpl.bars(%s,\n    title='%s'\n    saveas=%s\n    relpath=%s\n    graphicsformat=%s)\n\n""" %
                  (dbdat, title, repr(saveas), repr(relpath), repr(graphicsformat)))
        else:
            # if running from Canopy, call mpl directly
            filedict = mpl.bars(dbdat, title=title,
                                view=view,
                                saveas=saveas, relpath=relpath, graphicsformat=graphicsformat)
            return filedict

    # def get_pec_weightinfo(self):
    #     """
    #
    #     """
    #     def closest(u, options):
    #         return max(options, key=lambda v: len(os.path.commonprefix([u, v])))
    #
    #     dbdat = {}
    #     for db, odb in self.dbdict.items():
    #         #dbix = self.dbdict.keys().index(db)
    #         oss = odb.oss['default']
    #         eqrxns = [rxn for rxn, rr in zip(oss.hrxn, oss.axis['Rrat']) if rr == 1.0]
    #         for rxnix, rxn in enumerate(oss.hrxn):
    #             dbrxn = '-'.join([db, rxn])
    #             rrat = oss.axis['Rrat'][rxnix]
    #             eq = closest(rxn, eqrxns)
    #             print rxn, rxnix, eq, rrat, dbrxn
    #             dbdat[dbrxn] = {'eq': eq, 'Rrat': rrat}
    #     return dbdat

    def plot_axis(self, axis, modelchem, benchmark='default', sset='default',
                  failoninc=True, verbose=False, color='sapt', view=True,
                  saveas=None, relpath=False, graphicsformat=['pdf']):
        """

        """
        dbdatdict = OrderedDict()
        for mc in modelchem:
            # compute errors
            errors, indiv = self.compute_statistics(mc, benchmark=benchmark, sset=sset,
                                                    failoninc=failoninc, verbose=verbose, returnindiv=True)
            # repackage
            dbdat = []
            for db, odb in self.dbdict.items():
                dbix = self.dbdict.keys().index(db)
                oss = odb.oss[self.sset[sset][dbix]]
                # TODO may need to make axis name distributable across wrappeddbs
                # TODO not handling mc present bm absent
                if indiv[db] is not None:
                    for rxn in oss.hrxn:
                        rxnix = oss.hrxn.index(rxn)
                        bm = self.mcs[benchmark][dbix]
                        bmpresent = False if (bm is None or bm not in odb.hrxn[rxn].data) else True
                        mcpresent = False if (self.mcs[mc][dbix] not in odb.hrxn[rxn].data) else True
                        entry = {'db': db,
                                 'sys': str(rxn),
                                 'color': odb.hrxn[rxn].color,
                                 'axis': oss.axis[axis][rxnix]}

                        if bmpresent:
                            entry['bmdata'] = odb.hrxn[rxn].data[self.mcs[benchmark][dbix]].value
                        else:
                            entry['bmdata'] = None

                        if mcpresent:
                            entry['mcdata'] = odb.hrxn[rxn].data[self.mcs[mc][dbix]].value
                        else:
                            continue

                        if bmpresent and mcpresent:
                            entry['error'] = [indiv[db][rxn][0]]
                        else:
                            entry['error'] = [None]
                        dbdat.append(entry)
            dbdatdict[fancify_mc_tag(mc).strip()] = dbdat

        pre, suf, mid = string_contrast(modelchem)
        title = """%s[%s]%s vs %s axis %s for %s subset %s""" % (pre, str(len(mid)), suf, benchmark, axis, self.dbse, sset)
        print(title)
        #for mc, dbdat in dbdatdict.items():
        #    print mc
        #    for d in dbdat:
        #        print '{:20s} {:8.2f}    {:8.2f} {:8.2f}'.format(d['sys'], d['axis'],
        #            0.0 if d['bmdata'] is None else d['bmdata'],
        #            0.0 if d['mcdata'] is None else d['mcdata'])
        # generate matplotlib instructions and call or print
        try:
            from . import mpl
            import matplotlib.pyplot as plt
        except ImportError:
            # if not running from Canopy, print line to execute from Canopy
            print("""filedict = mpl.valerr(%s,\n    color='%s',\n    title='%s',\n    xtitle='%s',\n    view=%s\n    saveas=%s\n    relpath=%s\n    graphicsformat=%s)\n\n""" %
                  (dbdat, color, title, axis, view, repr(saveas), repr(relpath), repr(graphicsformat)))
        else:
            # if running from Canopy, call mpl directly
            filedict = mpl.valerr(dbdatdict, color=color, title=title, xtitle=axis,
                                  view=view,
                                  saveas=saveas, relpath=relpath, graphicsformat=graphicsformat)
            return filedict

    def load_saptdata_frombfdb(self, sset='default',
        pythonpath='/Users/loriab/linux/bfdb/sapt_punt', failoninc=True):  # pythonpath=None
        """This is a stopgap function that loads sapt component data from
        sapt_punt in bfdb repo.

        """
        saptpackage = OrderedDict()
        for db, odb in self.dbdict.items():
            modname = 'sapt_' + odb.dbse
            if pythonpath is not None:
                sys.path.insert(1, pythonpath)
            else:
                sys.path.append(os.path.dirname(__file__) + '/../data')
            try:
                datamodule = __import__(modname)
            except ImportError:
                print("""\nPython module for database data %s failed to load\n\n""" % (modname))
                print("""\nSearch path that was tried:\n""")
                print(', '.join(map(str, sys.path)))
                raise ValidationError("Python module loading problem for database subset generator " + str(modname))

            try:
                saptdata = getattr(datamodule, 'DATA')
            except AttributeError:
                raise ValidationError("SAPT punt module does not contain DATA" + str(modname))
            saptmc = saptdata['SAPT MODELCHEM']

            dbix = self.dbdict.keys().index(db)
            for rxn, orxn in odb.hrxn.items():
                lss = self.sset[sset][dbix]
                if lss is not None:
                    if rxn in odb.sset[lss]:
                        dbrxn = orxn.dbrxn
                        try:
                            elst = saptdata['SAPT ELST ENERGY'][dbrxn]
                            exch = saptdata['SAPT EXCH ENERGY'][dbrxn]
                            ind = saptdata['SAPT IND ENERGY'][dbrxn]
                            disp = saptdata['SAPT DISP ENERGY'][dbrxn]
                        except (KeyError, AttributeError):
                            print("""Warning: DATA['SAPT * ENERGY'] missing for reaction %s""" % (dbrxn))
                            if failoninc:
                                break
                        else:
                            if not all([elst, ind, disp]):  # exch sometimes physically zero
                                print("""Warning: DATA['SAPT * ENERGY'] missing piece for reaction %s: %s""" % (dbrxn, [elst, exch, ind, disp]))
                                if failoninc:
                                    break
                        saptpackage[dbrxn] = {'mc': saptmc,
                                              'elst': elst,
                                              'exch': exch,
                                              'ind': ind,
                                              'disp': disp}
        return saptpackage

    def plot_ternary(self, sset='default', labeled=True,
        pythonpath='/Users/loriab/linux/bfdb/sapt_punt', failoninc=True,  # pythonpath=None
        view=True,
        saveas=None, relpath=False, graphicsformat=['pdf']):
        """This is a stopgap function that loads sapt component data from
        sapt_punt in bfdb repo, then formats it to plot a ternary diagram.

        """
        saptdata = self.load_saptdata_frombfdb(sset=sset, pythonpath=pythonpath,
            failoninc=failoninc)

        dbdat = []
        mcs = []
        for dat in saptdata.values():
            dbdat.append([dat['elst'], dat['ind'], dat['disp']])
            if dat['mc'] not in mcs:
                mcs.append(dat['mc'])

        title = ' '.join([self.dbse, sset, ' '.join(mcs)])

        # generate matplotlib instructions and call or print
        try:
            from . import mpl
            import matplotlib.pyplot as plt
        except ImportError:
            pass
            # if not running from Canopy, print line to execute from Canopy
        else:
            # if running from Canopy, call mpl directly
            filedict = mpl.ternary(dbdat, title=title, labeled=labeled,
                                   view=view,
                                   saveas=saveas, relpath=relpath, graphicsformat=graphicsformat)
            return filedict

    def plot_flat(self, modelchem, benchmark='default', sset='default',
                  failoninc=True, verbose=False, color='sapt', xlimit=4.0, xlines=[0.0, 0.3, 1.0],
                  view=True,
                  saveas=None, relpath=False, graphicsformat=['pdf']):
        """Computes individual errors and summary statistics for single
        model chemistry *modelchem* versus *benchmark* over
        subset *sset* over all component databases. Thread *color* can be
        'rgb' for old coloring, a color name or 'sapt' for spectrum coloring.

        *saveas* conveys directory ('/') and/or filename for saving the
        resulting plot. File extension is not accessible, but *graphicsformat*
        array requests among 'png', 'pdf', and 'eps' formats. *relpath*
        forces paths to saved files to be relative to current directory,
        rather than absolute paths for returned code and file dictionary.

        Prepares flat diagram instructions and either executes them if
        matplotlib available (Canopy or Anaconda) or prints them. Returns a
        dictionary of all saved plot filenames.

        asdf.plot_flat('CCSD-CP-atqzadz', failoninc=False)
        """
        # compute errors
        mc = modelchem
        errors, indiv = self.compute_statistics(mc, benchmark=benchmark, sset=sset,
                                                failoninc=failoninc, verbose=verbose, returnindiv=True)
        # repackage
        dbdat = []
        for db, odb in self.dbdict.items():
            if indiv[db] is not None:
                for rxn in indiv[db].keys():
                    dbdat.append({'db': db,
                                  'sys': str(rxn),
                                  'color': odb.hrxn[rxn].color,
                                  'data': [indiv[db][rxn][0]]})
        pre, suf, mid = string_contrast(mc)
        title = self.dbse + '-' + sset + ' ' + pre + '[]' + suf
        mae = errors[self.dbse]['mae']
        mape = None
        # mape = 100 * errors[self.dbse]['mape']
        mapbe = None
        # generate matplotlib instructions and call or print
        try:
            from . import mpl
            import matplotlib.pyplot as plt
        except ImportError:
            # if not running from Canopy, print line to execute from Canopy
            print("""filedict = mpl.flat(%s,\n    color='%s',\n    title='%s',\n    mae=%s,\n    mape=%s,\n    xlimit=%s,\n    xlines=%s,\n    view=%s\n    saveas=%s\n    relpath=%s\n    graphicsformat=%s)\n\n""" %
                  (dbdat, color, mc, mae, mape, xlimit, repr(xlines), view, repr(saveas), repr(relpath), repr(graphicsformat)))
        else:
            # if running from Canopy, call mpl directly
            filedict = mpl.flat(dbdat, color=color, title=mc, mae=mae, mape=mape,
                                xlimit=xlimit, xlines=xlines, view=view,
                                saveas=saveas, relpath=relpath, graphicsformat=graphicsformat)
            return filedict

    def write_xyz_files(self, path=None):
        """Writes xyz files for every reagent in the Database to directory
        in *path* or to directory dbse_xyzfiles that it createsin cwd if
        *path* is None. Additionally, writes a script to that directory
        that will generate transparent-background ray-traced png files for
        every reagent with PyMol.

        """
        if path is None:
            xyzdir = os.getcwd() + os.sep + self.dbse + '_xyzfiles' + os.sep
        else:
            xyzdir = os.path.abspath(path) + os.sep
        if not os.path.exists(xyzdir):
            os.mkdir(xyzdir)

        for rgt, orgt in self.hrgt.items():
            omol = Molecule(orgt.mol)
            omol.update_geometry()
            omol.save_xyz(xyzdir + rgt + '.xyz')

        with open(xyzdir + 'pymol_xyz2png_script.pml', 'w') as handle:
            handle.write("""
# Launch PyMOL and run from its command line:
# PyMOL> cd {}
# PyMOL> @{}
""".format(xyzdir, 'pymol_xyz2png_script.pml'))
            for rgt in self.hrgt.keys():
                handle.write("""
load {xyzfile}
hide lines
show sticks
color grey, name c
cmd.set('''opaque_background''','''0''',quiet=0)
reset
orient
cmd.zoom(buffer=0.3, complete=1)
ray
png {pngfile}
reinitialize
""".format(
                    xyzfile=xyzdir + rgt + '.xyz',
                    pngfile=xyzdir + rgt + '.png'))

    def plot_all_flats(self, modelchem=None, sset='default', xlimit=4.0,
                       failoninc=True,
                       saveas=None, relpath=False, graphicsformat=['pdf']):
        """Generate pieces for inclusion into tables. Supply list of
        modelchemistries to plot from *modelchem*, otherwise defaults to
        all those available. Can modify subset *sset* and plotting
        range *xlimit*.

        >>> asdf.plot_all_flats(sset='tt-5min', xlimit=4.0)
        """
        mcs = self.mcs.keys() if modelchem is None else modelchem
        filedict = OrderedDict()
        for mc in sorted(mcs):
            minifiledict = self.plot_flat(mc, sset=sset, xlimit=xlimit, view=False,
                                          failoninc=failoninc,
                                          saveas=saveas, relpath=relpath, graphicsformat=graphicsformat)
            filedict[mc] = minifiledict
        return filedict

    def get_hrxn(self, sset='default'):
        """

        """
        rhrxn = OrderedDict()
        for db, odb in self.dbdict.items():
            dbix = self.dbdict.keys().index(db)
            lss = self.sset[sset][dbix]
            if lss is not None:
                for rxn in odb.hrxn:
                    if rxn in odb.sset[lss]:
                        orxn = odb.hrxn[rxn]
                        rhrxn[orxn.dbrxn] = orxn  # this is a change and conflict with vergil version
        return rhrxn

    def get_hrgt(self, sset='default', actv='default'):
        """

        """
        rhrxn = self.get_hrxn(sset=sset)
        rhrgt = OrderedDict()
        for rxn, orxn in rhrxn.items():
            for orgt in orxn.rxnm[actv].keys():
                rhrgt[orgt.name] = orgt
        # TODO prob need to avoid duplicates or pass

        return rhrgt

    def get_reactions(self, modelchem, sset='default', benchmark='default',
                      failoninc=True):
        """Collects the reactions present in *sset* from each WrappedDatabase,
        checks that *modelchem* and *benchmark* ReactionDatum are present
        (fails if *failoninc* True), then returns in an array a tuple for
        each reaction containing the modelchem key needed to access
        *modelchem*, the modelchem key needed to access *benchmark*, and
        the Reaction object.

        """
        dbdat = []
        rhrxn = self.get_hrxn(sset=sset)
        for orxn in rhrxn.itervalues():
            dbix = self.dbdict.keys().index(orxn.dbrxn.split('-')[0])
            lmc = self.mcs[modelchem][dbix]
            lbm = self.mcs[benchmark][dbix]
            try:
                orxn.data[lbm]
            except KeyError as e:
                # not sure if should treat bm differently
                lbm = None
            try:
                orxn.data[lmc]
            except KeyError as e:
                if failoninc:
                    raise e
                else:
                    lmc = None
            dbdat.append((lmc, lbm, orxn))
            # this is diff in that returning empties not just pass over- may break bfdb
            # try:
            #     orxn.data[lmc]
            #     orxn.data[lbm]
            # except KeyError as e:
            #     if failoninc:
            #         raise e
            #     else:
            #         # not sure yet if should return empties or just pass over
            #         pass
            # else:
            #     dbdat.append((lmc, lbm, orxn))
        return dbdat

    def get_missing_reactions(self, modelchem, sset='default'):
        """Returns a dictionary (keys self.dbse and all component
        WrappedDatabase.dbse) of two elements, the first being the number
        of reactions *sset* should contain and the second being a list of
        the reaction names (dbrxn) not available for *modelchem*. Absence
        of benchmark not considered.

        """
        counts = OrderedDict()
        counts[self.dbse] = [0, []]
        soledb = True if (len(self.dbdict) == 1 and self.dbdict.items()[0][0] == self.dbse) else False
        if not soledb:
            for db in self.dbdict.keys():
                counts[db] = [0, []]
        for (lmc, lbm, orxn) in self.get_reactions(modelchem, benchmark='default',
                                                   sset=sset, failoninc=False):
            db, rxn = orxn.dbrxn.split('-', 1)
            mcdatum = orxn.data[lmc].value if lmc else None
            counts[self.dbse][0] += 1
            if not soledb:
                counts[db][0] += 1
            if mcdatum is None:
                counts[self.dbse][1].append(orxn.dbrxn)
                if not soledb:
                    counts[db][1].append(orxn.dbrxn)
        return counts

    def plot_disthist(self, modelchem, benchmark='default', sset='default',
                      failoninc=True, verbose=False, xtitle='', view=True,
                      saveas=None, relpath=False, graphicsformat=['pdf']):
        """Computes individual errors and summary statistics for single
        model chemistry *modelchem* versus *benchmark* over
        subset *sset* over all component databases. Computes histogram
        of errors and gaussian distribution.

        *saveas* conveys directory ('/') and/or filename for saving the
        resulting plot. File extension is not accessible, but *graphicsformat*
        array requests among 'png', 'pdf', and 'eps' formats. *relpath*
        forces paths to saved files to be relative to current directory,
        rather than absolute paths for returned code and file dictionary.

        Prepares disthist diagram instructions and either executes them if
        matplotlib available (Canopy or Anaconda) or prints them. Returns a
        dictionary of all saved plot filenames.

        >>>
        """
        # compute errors
        mc = modelchem
        errors, indiv = self.compute_statistics(mc, benchmark=benchmark, sset=sset,
                                                failoninc=failoninc, verbose=verbose, returnindiv=True)
        # repackage
        dbdat = []
        for db in self.dbdict.keys():
            if indiv[db] is not None:
                for rxn in indiv[db].keys():
                    dbdat.append(indiv[db][rxn][0])
        title = """%s vs %s for %s subset %s""" % (mc, benchmark, self.dbse, sset)
        me = errors[self.dbse]['me']
        stde = errors[self.dbse]['stde']
        # generate matplotlib instructions and call or print
        try:
            from . import mpl
            import matplotlib.pyplot as plt
        except ImportError:
            # if not running from Canopy, print line to execute from Canopy
            print("""filedict = mpl.disthist(%s,\n    title='%s',\n    xtitle='%s'\n    me=%s,\n    stde=%s,\n    saveas=%s,\n    relpath=%s\n    graphicsformat=%s)\n\n""" %
                  (dbdat, title, xtitle, me, stde, repr(saveas), repr(relpath), repr(graphicsformat)))
        else:
            # if running from Canopy, call mpl directly
            filedict = mpl.disthist(dbdat, title=title, xtitle=xtitle, me=me, stde=stde,
                                    view=view,
                                    saveas=saveas, relpath=relpath, graphicsformat=graphicsformat)
            return filedict

    def plot_modelchems(self, modelchem, benchmark='default', mbenchmark=None,
                        sset='default', msset=None, failoninc=True, verbose=False, color='sapt',
                        xlimit=4.0, labeled=True, view=True,
                        mousetext=None, mouselink=None, mouseimag=None, mousetitle=None, mousediv=None,
                        saveas=None, relpath=False, graphicsformat=['pdf']):
        """Computes individual errors and summary statistics over all component
        databases for each model chemistry in array *modelchem* versus *benchmark*
        over subset *sset*. *mbenchmark* and *msset* are array options (same
        length as *modelchem*) that override *benchmark* and *sset*, respectively,
        for non-uniform specification. Thread *color* can be 'rgb' for old
        coloring, a color name or 'sapt' for spectrum coloring.

        *saveas* conveys directory ('/') and/or filename for saving the
        resulting plot. File extension is not accessible, but *graphicsformat*
        array requests among 'png', 'pdf', and 'eps' formats. *relpath*
        forces paths to saved files to be relative to current directory,
        rather than absolute paths for returned code and file dictionary.

        Prepares thread diagram instructions and either executes them if
        matplotlib available (Canopy or Anaconda) or prints them. Returns a
        dictionary of all saved plot filenames. If any of *mousetext*, *mouselink*,
        or *mouseimag* is specified, htmlcode will be returned with an image map of
        slats to any of text, link, or image, respectively.

        """
        # distribute benchmark
        if mbenchmark is None:
            lbenchmark = [benchmark] * len(modelchem)  # normal bm modelchem name
        else:
            if isinstance(mbenchmark, basestring) or len(mbenchmark) != len(modelchem):
                raise ValidationError(
                    """mbenchmark must be array of length distributable among modelchem""" % (str(mbenchmark)))
            else:
                lbenchmark = mbenchmark  # array of bm for each modelchem
        # distribute sset
        if msset is None:
            lsset = [sset] * len(modelchem)  # normal ss name like 'MX'
        else:
            if isinstance(msset, basestring) or len(msset) != len(modelchem):
                raise ValidationError("""msset must be array of length distributable among modelchem""" % (str(msset)))
            else:
                lsset = msset  # array of ss for each modelchem
        # compute errors
        index = []
        errors = {}
        indiv = {}
        for mc, bm, ss in zip(modelchem, lbenchmark, lsset):
            ix = '%s_%s_%s' % (ss, mc, bm)
            index.append(ix)
            errors[ix], indiv[ix] = self.compute_statistics(mc, benchmark=bm, sset=ss,
                                                            failoninc=failoninc, verbose=verbose, returnindiv=True)
        # repackage
        dbdat = []
        for db, odb in self.dbdict.items():
            dbix = self.dbdict.keys().index(db)
            for rxn in odb.hrxn:
                data = []
                for ix in index:
                    if indiv[ix][db] is not None:
                        if rxn in odb.sset[self.sset[lsset[index.index(ix)]][dbix]]:
                            try:
                                data.append(indiv[ix][db][rxn][0])
                            except KeyError as e:
                                if failoninc:
                                    raise e
                                else:
                                    data.append(None)
                        else:
                            data.append(None)
                    else:
                        data.append(None)
                if not data or all(item is None for item in data):
                    pass  # filter out empty reactions
                else:
                    dbdat.append({'db': db,
                                  'sys': str(rxn),
                                  'show': str(rxn),
                                  'color': odb.hrxn[rxn].color,
                                  'data': data})
        mae = [errors[ix][self.dbse]['mae'] for ix in index]
        mape = [100 * errors[ix][self.dbse]['mape'] for ix in index]
        # form unique filename
        ixpre, ixsuf, ixmid = string_contrast(index)
        title = self.dbse + ' ' + ixpre + '[]' + ixsuf
        # generate matplotlib instructions and call or print
        try:
            from . import mpl
            import matplotlib.pyplot as plt
        except ImportError:
            # if not running from Canopy, print line to execute from Canopy
            print("""filedict, htmlcode = mpl.threads(%s,\n    color='%s',\n    title='%s',\n    labels=%s,\n    mae=%s,\n    mape=%s\n    xlimit=%s\n    labeled=%s\n    saveas=%s\n    mousetext=%s\n    mouselink=%s\n    mouseimag=%s\n    mousetitle=%s,\n    mousediv=%s,\n    relpath=%s\n    graphicsformat=%s)\n\n""" %
                  (dbdat, color, title, ixmid, mae, mape, str(xlimit),
                  repr(labeled), repr(saveas), repr(mousetext), repr(mouselink), repr(mouseimag),
                  repr(mousetitle), repr(mousediv), repr(relpath), repr(graphicsformat)))
        else:
            # if running from Canopy, call mpl directly
            filedict, htmlcode = mpl.threads(dbdat, color=color, title=title, labels=ixmid, mae=mae, mape=mape,
                                             xlimit=xlimit, labeled=labeled, view=view,
                                             mousetext=mousetext, mouselink=mouselink,
                                             mouseimag=mouseimag, mousetitle=mousetitle, mousediv=mousediv,
                                             saveas=saveas, relpath=relpath, graphicsformat=graphicsformat)
            return filedict, htmlcode

    def plot_liliowa(self, modelchem, benchmark='default',
                     failoninc=True, xlimit=2.0, view=True,
                     saveas=None, relpath=False, graphicsformat=['pdf']):
        """

        Note that not possible to access sset of component databases. That is, for Database SSIBBI, SSI-only arylaryl is accessible b/c not defined in BBI, but SSI-only neutral is not accessible.
        """
        # compute errors
        mc = modelchem
        errors = {}
        for ss in self.sset.keys():
            errors[ss] = self.compute_statistics(mc, benchmark=benchmark, sset=ss,
                                                 failoninc=failoninc, verbose=False, returnindiv=False)

        # repackage
        dbdat = []
        ssarray = ['pospos', 'posneg', 'pospolar', 'posaliph', 'posaryl',
                   None, 'negneg', 'negpolar', 'negaliph', 'negaryl',
                   None, None, 'polarpolar', 'polaraliph', 'polararyl',
                   None, None, None, 'aliphaliph', 'alipharyl',
                   None, None, None, None, 'arylaryl']
        for ss in ssarray:
            dbdat.append(0.0 if ss is None else errors[ss][self.dbse]['mae'])

        # generate matplotlib instructions and call or print
        try:
            from . import mpl
            import matplotlib.pyplot as plt
        except ImportError:
            print('Matplotlib not avail')
        else:
            filedict = mpl.liliowa(dbdat, xlimit=xlimit, view=view,
                                saveas=saveas, relpath=relpath, graphicsformat=graphicsformat)
            return filedict

    def plot_iowa(self, modelchem, benchmark='default', sset='default',
                  failoninc=True, verbose=False,
                  title='', xtitle='', xlimit=2.0,
                  view=True,
                  saveas=None, relpath=False, graphicsformat=['pdf']):
        """Computes individual errors for single *modelchem* versus
        *benchmark* over subset *sset*. Coloring green-to-purple with
        maximum intensity at *xlimit*. Prepares Iowa plot instructions and
        either executes them if matplotlib available (Canopy) or prints them.

        """
        title = self.dbse + ' ' + modelchem
        # compute errors
        mc = modelchem
        errors, indiv = self.compute_statistics(mc, benchmark=benchmark, sset=sset,
                                                failoninc=failoninc, verbose=verbose, returnindiv=True)
        # repackage
        dbdat = []
        dblbl = []
        for db in self.dbdict.keys():
            if indiv[db] is not None:
                for rxn in indiv[db].keys():
                    dbdat.append(indiv[db][rxn][0])
                    dblbl.append(str(rxn))
        title = """%s vs %s for %s subset %s""" % (mc, benchmark, self.dbse, sset)
        me = errors[self.dbse]['me']
        # generate matplotlib instructions and call or print
        try:
            from . import mpl
            import matplotlib.pyplot as plt
        except ImportError:
            # if not running from Canopy, print line to execute from Canopy
            print("""mpl.iowa(%s,\n    %s,\n    title='%s',\n    xtitle='%s'\n    xlimit=%s,\n    saveas=%s,\n    relpath=%s\n    graphicsformat=%s)\n\n""" %
                  (dbdat, dblbl, title, xtitle, xlimit, repr(saveas), repr(relpath), repr(graphicsformat)))
        else:
            # if running from Canopy, call mpl directly
            filedict = mpl.iowa(dbdat, dblbl, title=title, xtitle=xtitle, xlimit=xlimit,
                                view=view,
                                saveas=saveas, relpath=relpath, graphicsformat=graphicsformat)
            return filedict

    def export_pandas(self, modelchem=[], benchmark='default', sset='default', modelchemlabels=None,
                      failoninc=True):
        """
        *modelchem* is array of model chemistries, if modelchem is empty, get only benchmark
        is benchmark needed?
        """
        import pandas as pd
        import numpy as np

        if self.dbse not in ['ACONF', 'SCONF', 'PCONF', 'CYCONF']:
            saptdata = self.load_saptdata_frombfdb(sset=sset, pythonpath='/Users/loriab/linux/bfdb/sapt_punt',
                failoninc=failoninc)

        listodicts = []
        rhrxn = self.get_hrxn(sset=sset)
        for dbrxn, orxn in rhrxn.items():
            wdb = dbrxn.split('-')[0]
            dbix = self.dbdict.keys().index(wdb)
            wbm = self.mcs[benchmark][dbix]
            wss = self.sset[sset][dbix]
            woss = self.dbdict[wdb].oss[wss]
            try:
                Rrat = woss.axis['Rrat'][woss.hrxn.index(orxn.name)]
            except KeyError:
                Rrat = 1.0  # TODO generic soln?

            dictorxn = {}
            dictorxn['DB'] = wdb
            dictorxn['System'] = orxn.tagl
            dictorxn['Name'] = orxn.name
            dictorxn['R'] = Rrat
            dictorxn['System #'] = orxn.indx
            dictorxn['Benchmark'] = np.NaN if orxn.benchmark is None else orxn.data[
                wbm].value  # this NaN exception is new and experimental
            dictorxn['QcdbSys'] = orxn.dbrxn

            if self.dbse not in ['ACONF', 'SCONF', 'PCONF', 'CYCONF']:
                dictorxn['SAPT ELST ENERGY'] = saptdata[dbrxn]['elst']
                dictorxn['SAPT EXCH ENERGY'] = saptdata[dbrxn]['exch']
                dictorxn['SAPT IND ENERGY'] = saptdata[dbrxn]['ind']
                dictorxn['SAPT DISP ENERGY'] = saptdata[dbrxn]['disp']
                dictorxn['SAPT TOTAL ENERGY'] = dictorxn['SAPT ELST ENERGY'] + dictorxn['SAPT EXCH ENERGY'] + \
                                                dictorxn['SAPT IND ENERGY'] + dictorxn['SAPT DISP ENERGY']

            orgts = orxn.rxnm['default'].keys()
            omolD = Molecule(orgts[0].mol)  # TODO this is only going to work with Reaction ~= Reagent databases
            npmolD = omolD.format_molecule_for_numpy()
            omolA = Molecule(orgts[1].mol)  # TODO this is only going to work with Reaction ~= Reagent databases
            omolA.update_geometry()
            dictorxn['MonA'] = omolA.natom()

            # this whole member fn not well defined for db of varying stoichiometry
            if self.dbse in ['ACONF', 'SCONF', 'PCONF', 'CYCONF']:
                npmolD = omolD.format_molecule_for_numpy()
                npmolA = omolA.format_molecule_for_numpy()
                dictorxn['Geometry'] = np.vstack([npmolD, npmolA])
            else:
                dictorxn['Geometry'] = omolD.format_molecule_for_numpy()
            # print '\nD', npmolD.shape[0], npmolA.shape[0], dictorxn['MonA'], npmolD, npmolA, dictorxn['Geometry']

            for mc in modelchem:
                try:
                    wmc = self.mcs[mc][dbix]
                except KeyError:
                    # modelchem not in Database at all
                    print(mc, 'not found')
                    continue
                key = mc if modelchemlabels is None else modelchemlabels[modelchem.index(mc)]
                try:
                    dictorxn[key] = orxn.data[wmc].value
                except KeyError as e:
                    # reaction not in modelchem
                    if failoninc:
                        raise ValidationError("""Reaction %s missing datum %s.""" % (key, str(e)))
                    else:
                        print(mc, str(e), 'not found')
                        continue
            listodicts.append(dictorxn)

        df = pd.DataFrame(listodicts)
        pd.set_option('display.width', 500)
        print(df.head(5))
        print(df.tail(5))
        return df

    def table_reactions(self, modelchem, benchmark='default', sset='default',
                        failoninc=True,
                        columnplan=['indx', 'tagl', 'bm', 'mc', 'e', 'pe'],
                        title="""Reaction energies [kcal/mol] for {sset} $\subset$ {dbse} with {mc}""",
                        indextitle="""Detailed results for {sset} $\subset$ {dbse} with {mc}""",
                        plotpath='analysis/mols/',
                        standalone=True, theme='rxns', filename=None):
        r"""Prepare single LaTeX table to *filename* or return lines if None showing
        the per-reaction results for reactions in *sset* for single or array
        or 'all' *modelchem*, where the last uses self.mcs(), model chemistries
        versus *benchmark*. Use *failoninc* to toggle between command failing
        or blank lines in table. Use *standalone* to toggle between full
        compilable document and suitable for inclusion in another LaTeX document.
        Use *columnplan* to customize column (from among columnreservoir, below)
        layout. Use *title* and *indextitle* to customize table caption and
        table-of-contents caption, respectively; variables in curly braces will
        be substituted. Use *theme* to customize the \ref{tbl:} code.

        """
        # define eligible columns for inclusion
        columnreservoir = {
            'dbrxn': ['l', r"""\textbf{Reaction}""", """{0:25s}"""],
            'indx': ['r', '', """{0:14s}"""],
            'tagl': ['l', r"""\textbf{Reaction}""", """{0:50s}"""],
            'bm': ['d', r"""\multicolumn{1}{c}{\textbf{Benchmark}}""", """{0:8.2f}"""],
            'mc': ['d', r"""\multicolumn{1}{c}{\textbf{ModelChem}}""", """{0:8.2f}"""],
            'e': ['d', r"""\multicolumn{1}{c}{\textbf{Error}}""", """{0:8.2f}"""],
            'pe': ['d', r"""\multicolumn{1}{c}{\textbf{\% Err.}}""", """{0:8.1f}"""],
            'imag': ['l', '', r"""\includegraphics[width=1.0cm,height=3.5mm]{%s%%ss.png}""" % (plotpath)],  # untested
        }
        for col in columnplan:
            if col not in columnreservoir.keys():
                raise ValidationError('Column {0} not recognized. Register with columnreservoir.'.format(col))

        if isinstance(modelchem, basestring):
            if modelchem.lower() == 'all':
                mcs = sorted(self.mcs.keys())
            else:
                mcs = [modelchem]
        else:
            mcs = modelchem

        # commence to generate LaTeX code
        tablelines = []
        indexlines = []

        if standalone:
            tablelines += textables.begin_latex_document()

        # iterate to produce one LaTeX table per modelchem
        for mc in mcs:
            # prepare summary statistics
            perr = self.compute_statistics(mc, benchmark=benchmark, sset=sset,
                                           failoninc=failoninc, verbose=False,
                                           returnindiv=False)
            serrors = OrderedDict()
            for db in self.dbdict.keys():
                serrors[db] = None if perr[db] is None else format_errors(perr[db], mode=3)
            serrors[self.dbse] = format_errors(perr[self.dbse], mode=3)

            # prepare individual reactions and errors
            terrors = OrderedDict()
            isComplete = True
            for (lmc, lbm, orxn) in self.get_reactions(mc, benchmark=benchmark,
                                                       sset=sset, failoninc=failoninc):
                tmp = {}
                dbrxn = orxn.dbrxn
                tmp['dbrxn'] = dbrxn.replace('_', '\\_')
                tmp['indx'] = r"""\textit{""" + str(orxn.indx) + """}"""
                tmp['tagl'] = dbrxn.split('-')[0] + ' ' + \
                              (orxn.latex if orxn.latex else orxn.tagl.replace('_', '\\_'))
                tmp['imag'] = None  # name of primary rgt
                bmdatum = orxn.data[lbm].value if lbm else None
                mcdatum = orxn.data[lmc].value if lmc else None
                tmp['bm'] = bmdatum
                tmp['mc'] = mcdatum
                if lmc and lbm:
                    tmp['e'] = mcdatum - bmdatum
                    tmp['pe'] = 100 * (mcdatum - bmdatum) / abs(bmdatum)
                    # TODO redefining errors not good practice
                else:
                    isComplete = False
                    tmp['e'] = None
                    tmp['pe'] = None

                terrors[dbrxn] = {}
                for c in columnreservoir.keys():
                    terrors[dbrxn][c] = '' if tmp[c] is None else \
                        columnreservoir[c][2].format(tmp[c])

            fancymodelchem = self.fancy_mcs(latex=True)[mc]
            thistitle = title.format(dbse=self.dbse, mc=fancymodelchem,
                                     sset='All' if sset == 'default' else sset.upper())
            lref = [r"""tbl:qcdb"""]
            if theme:
                lref.append(theme)
            lref.append(self.dbse)
            if sset != 'default':
                lref.append(sset)
            lref.append(mc)
            ref = '-'.join(lref)

            # table intro
            tablelines.append(r"""\begingroup""")
            tablelines.append(r"""\squeezetable""")
            tablelines.append(r"""\LTcapwidth=\textwidth""")
            tablelines.append(r"""\begin{longtable}{%s}""" % (''.join([columnreservoir[col][0] for col in columnplan])))
            tablelines.append(r"""\caption{%s""" % (thistitle))
            tablelines.append(r"""\label{%s}} \\ """ % (ref))
            tablelines.append(r"""\hline\hline""")

            columntitles = [columnreservoir[col][1] for col in columnplan]
            # initial header
            tablelines.append(' & '.join(columntitles) + r""" \\ """)
            tablelines.append(r"""\hline""")
            tablelines.append(r"""\endfirsthead""")
            # to be continued header
            tablelines.append(r"""\multicolumn{%d}{@{}l}{\textit{\ldots continued} %s} \\ """ %
                              (len(columnplan), fancymodelchem))
            tablelines.append(r"""\hline\hline""")
            tablelines.append(' & '.join(columntitles) + r""" \\ """)
            tablelines.append(r"""\hline""")
            tablelines.append(r"""\endhead""")
            # to be continued footer
            tablelines.append(r"""\hline\hline""")
            tablelines.append(r"""\multicolumn{%d}{r@{}}{\textit{continued \ldots}} \\ """ %
                              (len(columnplan)))
            tablelines.append(r"""\endfoot""")
            # final footer
            tablelines.append(r"""\hline\hline""")
            tablelines.append(r"""\endlastfoot""")

            # table body
            for dbrxn, stuff in terrors.items():
                tablelines.append(' & '.join([stuff[col] for col in columnplan]) + r""" \\ """)

            # table body summary
            if any(col in ['e', 'pe'] for col in columnplan):
                field_to_put_labels = [col for col in ['tagl', 'dbrxn', 'indx'] if col in columnplan]
                if field_to_put_labels:

                    for block, blkerrors in serrors.items():
                        if blkerrors:  # skip e.g., NBC block in HB of DB4
                            tablelines.append(r"""\hline""")
                            summlines = [[] for i in range(8)]
                            for col in columnplan:
                                if col == field_to_put_labels[0]:
                                    summlines[0].append(
                                        r"""\textbf{Summary Statistics: %s%s}%s""" % \
                                        ('' if sset == 'default' else sset + r""" $\subset$ """,
                                         block,
                                         '' if isComplete else r""", \textit{partial}"""))
                                    summlines[1].append(r"""\textit{Minimal Signed Error}   """)
                                    summlines[2].append(r"""\textit{Minimal Absolute Error} """)
                                    summlines[3].append(r"""\textit{Maximal Signed Error}   """)
                                    summlines[4].append(r"""\textit{Maximal Absolute Error} """)
                                    summlines[5].append(r"""\textit{Mean Signed Error}      """)
                                    summlines[6].append(r"""\textit{Mean Absolute Error}    """)
                                    summlines[7].append(r"""\textit{Root-Mean-Square Error} """)
                                elif col in ['e', 'pe']:
                                    summlines[0].append('')
                                    summlines[1].append(blkerrors['nex' + col])
                                    summlines[2].append(blkerrors['min' + col])
                                    summlines[3].append(blkerrors['pex' + col])
                                    summlines[4].append(blkerrors['max' + col])
                                    summlines[5].append(blkerrors['m' + col])
                                    summlines[6].append(blkerrors['ma' + col])
                                    summlines[7].append(blkerrors['rms' + col])
                                else:
                                    for ln in range(len(summlines)):
                                        summlines[ln].append('')
                            for ln in range(len(summlines)):
                                tablelines.append(' & '.join(summlines[ln]) + r""" \\ """)

            # table conclusion
            tablelines.append(r"""\end{longtable}""")
            tablelines.append(r"""\endgroup""")
            tablelines.append(r"""\clearpage""")
            tablelines.append('\n\n')

            # form table index
            thisindextitle = indextitle.format(dbse=self.dbse, mc=fancymodelchem.strip(),
                                               sset='All' if sset == 'default' else sset.upper())
            indexlines.append(r"""\scriptsize \ref{%s} & \scriptsize %s \\ """ % \
                              (ref, thisindextitle))

        if standalone:
            tablelines += textables.end_latex_document()

        # form table and index return structures
        if filename is None:
            return tablelines, indexlines
        else:
            if filename.endswith('.tex'):
                filename = filename[:-4]
            with open(filename + '.tex', 'w') as handle:
                handle.write('\n'.join(tablelines))
            with open(filename + '_index.tex', 'w') as handle:
                handle.write('\n'.join(indexlines) + '\n')
            print("""\n  LaTeX index written to {filename}_index.tex\n"""
                  """  LaTeX table written to {filename}.tex\n"""
                  """  >>> pdflatex {filename}\n"""
                  """  >>> open /Applications/Preview.app {filename}.pdf\n""".format(filename=filename))
            filedict = {'data': os.path.abspath(filename) + '.tex',
                        'index': os.path.abspath(filename + '_index.tex')}
            return filedict

    def table_wrapper(self, mtd, bas, tableplan, benchmark='default',
                      opt=['CP'], err=['mae'], sset=['default'], dbse=None,
                      opttarget=None,
                      failoninc=True,
                      xlimit=4.0, xlines=[0.0, 0.3, 1.0],
                      ialimit=2.0,
                      plotpath='autogen',
                      subjoin=True,
                      title=None, indextitle=None,
                      suppressblanks=False,
                      standalone=True, theme=None, filename=None):
        """Prepares dictionary of errors for all combinations of *mtd*, *opt*,
        *bas* with respect to model chemistry *benchmark*, mindful of *failoninc*.
        The general plan for the table, as well as defaults for landscape,
        footnotes, *title*, *indextitle, and *theme* are got from function
        *tableplan*. Once error dictionary is ready, it and all other arguments
        are passed along to textables.table_generic. Two arrays, one of table
        lines and one of index lines are returned unless *filename* is given,
        in which case they're written to file and a filedict returned.

        """
        # get plan for table from *tableplan* and some default values
        kwargs = {'plotpath': plotpath,
                  'subjoin': subjoin,
                  'xlines': xlines,
                  'xlimit': xlimit,
                  'ialimit': ialimit}
        rowplan, columnplan, landscape, footnotes, \
        suggestedtitle, suggestedtheme = tableplan(**kwargs)
        #suggestedtitle, suggestedtheme = tableplan(plotpath=plotpath, subjoin=subjoin)

        # make figure files write themselves
        autothread = {}
        autoliliowa = {}
        if plotpath == 'autogen':
            for col in columnplan:
                if col[3].__name__ == 'flat':
                    if col[4] and autothread:
                        print('TODO: merge not handled')
                    elif col[4] or autothread:
                        autothread.update(col[4])
                    else:
                        autothread = {'dummy': True}
                elif col[3].__name__ == 'liliowa':
                    autoliliowa = {'dummy': True}

        # negotiate some defaults
        dbse = [self.dbse] if dbse is None else dbse
        theme = suggestedtheme if theme is None else theme
        title = suggestedtitle if title is None else title
        indextitle = title if indextitle is None else indextitle
        opttarget = {'default': ['']} if opttarget is None else opttarget

        def unify_options(orequired, opossible):
            """Perform a merge of options tags in *orequired* and *opossible* so
            that the result is free of duplication and has the mode at the end.

            """
            opt_combos = []
            for oreq in orequired:
                for opos in opossible:
                    pieces = sorted(set(oreq.split('_') + opos.split('_')))
                    if '' in pieces:
                        pieces.remove('')
                    for mode in ['CP', 'unCP', 'SA']:
                        if mode in pieces:
                            pieces.remove(mode)
                            pieces.append(mode)
                    pieces = '_'.join(pieces)
                    opt_combos.append(pieces)
            return opt_combos

        # gather list of model chemistries for table
        mcs = ['-'.join(prod) for prod in itertools.product(mtd, opt, bas)]
        mc_translator = {}
        for m, o, b in itertools.product(mtd, opt, bas):
            nominal_mc = '-'.join([m, o, b])
            for oo in unify_options([o], opttarget['default']):
                trial_mc = '-'.join([m, oo, b])
                try:
                    perr = self.compute_statistics(trial_mc, benchmark=benchmark, sset='default',  # prob. too restrictive by choosing subset
                                                   failoninc=False, verbose=False, returnindiv=False)
                except KeyError as e:
                    continue
                else:
                    mc_translator[nominal_mc] = trial_mc
                    break
            else:
                mc_translator[nominal_mc] = None

        # compute errors
        serrors = {}
        for mc in mcs:
            serrors[mc] = {}
            for ss in self.sset.keys():
                serrors[mc][ss] = {}
                if mc_translator[mc] in self.mcs:
                    # Note: not handling when one component Wdb has one translated pattern and another another
                    perr = self.compute_statistics(mc_translator[mc], benchmark=benchmark, sset=ss,
                                                   failoninc=failoninc, verbose=False, returnindiv=False)
                    serrors[mc][ss][self.dbse] = format_errors(perr[self.dbse], mode=3)
                    if not failoninc:
                        mcsscounts = self.get_missing_reactions(mc_translator[mc], sset=ss)
                        serrors[mc][ss][self.dbse]['tgtcnt'] = mcsscounts[self.dbse][0]
                        serrors[mc][ss][self.dbse]['misscnt'] = len(mcsscounts[self.dbse][1])
                    if autothread:
                        if ('sset' in autothread and ss in autothread['sset']) or ('sset' not in autothread):
                            mcssplots = self.plot_flat(mc_translator[mc], benchmark=benchmark, sset=ss,
                                failoninc=failoninc, color='sapt', xlimit=xlimit, xlines=xlines, view=False,
                                saveas='flat_' + '-'.join([self.dbse, ss, mc]), relpath=True, graphicsformat=['pdf'])
                            serrors[mc][ss][self.dbse]['plotflat'] = mcssplots['pdf']
                    if autoliliowa and ss == 'default':
                            mcssplots = self.plot_liliowa(mc_translator[mc], benchmark=benchmark,
                                failoninc=failoninc, xlimit=ialimit, view=False,
                                saveas='liliowa_' + '-'.join([self.dbse, ss, mc]), relpath=True, graphicsformat=['pdf'])
                            serrors[mc][ss][self.dbse]['plotliliowa'] = mcssplots['pdf']
                    for db in self.dbdict.keys():
                        if perr[db] is None:
                            serrors[mc][ss][db] = None
                        else:
                            serrors[mc][ss][db] = format_errors(perr[db], mode=3)
                            if not failoninc:
                                serrors[mc][ss][db]['tgtcnt'] = mcsscounts[db][0]
                                serrors[mc][ss][db]['misscnt'] = len(mcsscounts[db][1])
                else:
                    serrors[mc][ss][self.dbse] = format_errors(initialize_errors(), mode=3)
                    for db in self.dbdict.keys():
                        serrors[mc][ss][db] = format_errors(initialize_errors(), mode=3)
        for key in serrors.keys():
            print("""{:>35}{:>35}{}""".format(key, mc_translator[key], serrors[key]['default'][self.dbse]['mae']))

        # find indices that would be neglected in a single sweep over table_generic
        keysinplan = set(sum([col[-1].keys() for col in columnplan], rowplan))
        obvious = {'dbse': dbse, 'sset': sset, 'mtd': mtd, 'opt': opt, 'bas': bas, 'err': err}
        for key, vari in obvious.items():
            if len(vari) == 1 or key in keysinplan:
                del obvious[key]
        iteroers = [(prod) for prod in itertools.product(*obvious.values())]

        # commence to generate LaTeX code
        tablelines = []
        indexlines = []

        if standalone:
            tablelines += textables.begin_latex_document()

        for io in iteroers:
            actvargs = dict(zip(obvious.keys(), [[k] for k in io]))
            nudbse = actvargs['dbse'] if 'dbse' in actvargs else dbse
            nusset = actvargs['sset'] if 'sset' in actvargs else sset
            numtd = actvargs['mtd'] if 'mtd' in actvargs else mtd
            nuopt = actvargs['opt'] if 'opt' in actvargs else opt
            nubas = actvargs['bas'] if 'bas' in actvargs else bas
            nuerr = actvargs['err'] if 'err' in actvargs else err

            table, index = textables.table_generic(
                mtd=numtd, bas=nubas, opt=nuopt, err=nuerr, sset=nusset, dbse=nudbse,
                rowplan=rowplan, columnplan=columnplan, serrors=serrors,
                plotpath='' if plotpath == 'autogen' else plotpath,
                subjoin=subjoin,
                title=title, indextitle=indextitle,
                suppressblanks=suppressblanks,
                landscape=landscape, footnotes=footnotes,
                standalone=False, theme=theme)

            tablelines += table
            tablelines.append('\n\n')
            indexlines += index

        if standalone:
            tablelines += textables.end_latex_document()

        # form table and index return structures
        if filename is None:
            return tablelines, indexlines
        else:
            if filename.endswith('.tex'):
                filename = filename[:-4]
            with open(filename + '.tex', 'w') as handle:
                handle.write('\n'.join(tablelines))
            with open(filename + '_index.tex', 'w') as handle:
                handle.write('\n'.join(indexlines))
            print("""\n  LaTeX index written to {filename}_index.tex\n"""
                  """  LaTeX table written to {filename}.tex\n"""
                  """  >>> pdflatex {filename}\n"""
                  """  >>> open /Applications/Preview.app {filename}.pdf\n""".format(filename=filename))
            filedict = {'data': os.path.abspath(filename) + '.tex',
                        'index': os.path.abspath(filename + '_index.tex')}
            return filedict

    def table_scrunch(self, plotpath, subjoin):
        rowplan = ['mtd']
        columnplan = [
            ['l', r'Method', '', textables.label, {}],
            ['c', r'Description', '', textables.empty, {}],
            ['d', r'aug-cc-pVDZ', 'unCP', textables.val, {'bas': 'adz', 'opt': 'unCP'}],
            ['d', r'aug-cc-pVDZ', 'CP', textables.val, {'bas': 'adz', 'opt': 'CP'}],
            ['d', r'aug-cc-pVTZ', 'unCP', textables.val, {'bas': 'atz', 'opt': 'unCP'}],
            ['d', r'aug-cc-pVTZ', 'CP', textables.val, {'bas': 'atz', 'opt': 'CP'}]]

        footnotes = []
        landscape = False
        theme = 'summavg'
        title = r"""Classification and Performance of model chemistries. Interaction energy [kcal/mol] {{err}} statistics.""".format()
        return rowplan, columnplan, landscape, footnotes, title, theme

    def table_merge_abbr(self, plotpath, subjoin):
        """Specialization of table_generic into table with minimal statistics
        (three S22 and three overall) plus embedded slat diagram as suitable
        for main paper. A single table is formed in sections by *bas* with
        lines *mtd* within each section.

        """
        rowplan = ['bas', 'mtd']
        columnplan = [
            ['l', r"""Method \& Basis Set""", '', textables.label, {}],
            ['d', r'S22', 'HB', textables.val, {'sset': 'hb', 'dbse': 'S22'}],
            ['d', r'S22', 'MX/DD', textables.val, {'sset': 'mxdd', 'dbse': 'S22'}],
            ['d', r'S22', 'TT', textables.val, {'sset': 'tt', 'dbse': 'S22'}],
            ['d', r'Overall', 'HB', textables.val, {'sset': 'hb', 'dbse': 'DB4'}],
            ['d', r'Overall', 'MX/DD', textables.val, {'sset': 'mxdd', 'dbse': 'DB4'}],
            ['d', r'Overall', 'TT', textables.val, {'sset': 'tt', 'dbse': 'DB4'}],
            ['l', r"""Error Distribution\footnotemark[1]""",
             r"""\includegraphics[width=6.67cm,height=3.5mm]{%s%s.pdf}""" % (plotpath, 'blank'),
             textables.graphics, {}],
            ['d', r'Time', '', textables.empty, {}]]
        # TODO Time column not right at all

        footnotes = [fnreservoir['blankslat']]
        landscape = False
        theme = 'smmerge'
        title = r"""Interaction energy [kcal/mol] {{err}} subset statistics with computed with {{opt}}{0}.""".format(
            '' if subjoin else r""" and {bas}""")
        return rowplan, columnplan, landscape, footnotes, title, theme

    def table_merge_suppmat(self, plotpath, subjoin):
        """Specialization of table_generic into table with as many statistics
        as will fit (mostly fullcurve and a few 5min) plus embedded slat
        diagram as suitable for supplementary material. Multiple tables are
        formed, one for each in *bas* with lines *mtd* within each table.

        """
        rowplan = ['bas', 'mtd']
        columnplan = [
            ['l', r"""Method \& Basis Set""", '', textables.label, {}],
            ['d', 'S22', 'HB', textables.val, {'sset': 'hb', 'dbse': 'S22'}],
            ['d', 'S22', 'MX', textables.val, {'sset': 'mx', 'dbse': 'S22'}],
            ['d', 'S22', 'DD', textables.val, {'sset': 'dd', 'dbse': 'S22'}],
            ['d', 'S22', 'TT', textables.val, {'sset': 'tt', 'dbse': 'S22'}],
            ['d', 'NBC10', 'MX', textables.val, {'sset': 'mx', 'dbse': 'NBC1'}],
            ['d', 'NBC10', 'DD', textables.val, {'sset': 'dd', 'dbse': 'NBC1'}],
            ['d', 'NBC10', 'TT', textables.val, {'sset': 'tt', 'dbse': 'NBC1'}],
            ['d', 'HBC6', 'HB/TT', textables.val, {'sset': 'tt', 'dbse': 'HBC1'}],
            ['d', 'HSG', 'HB', textables.val, {'sset': 'hb', 'dbse': 'HSG'}],
            ['d', 'HSG', 'MX', textables.val, {'sset': 'mx', 'dbse': 'HSG'}],
            ['d', 'HSG', 'DD', textables.val, {'sset': 'dd', 'dbse': 'HSG'}],
            ['d', 'HSG', 'TT', textables.val, {'sset': 'tt', 'dbse': 'HSG'}],
            ['d', 'Avg', 'TT ', textables.val, {'sset': 'tt', 'dbse': 'DB4'}],
            ['l', r"""Error Distribution\footnotemark[1]""",
             r"""\includegraphics[width=6.67cm,height=3.5mm]{%s%s.pdf}""" % (plotpath, 'blank'),
             textables.graphics, {}],
            ['d', 'NBC10', r"""TT\footnotemark[2]""", textables.val, {'sset': 'tt-5min', 'dbse': 'NBC1'}],
            ['d', 'HBC6', r"""TT\footnotemark[2] """, textables.val, {'sset': 'tt-5min', 'dbse': 'HBC1'}],
            ['d', 'Avg', r"""TT\footnotemark[2]""", textables.val, {'sset': 'tt-5min', 'dbse': 'DB4'}]]

        footnotes = [fnreservoir['blankslat'], fnreservoir['5min']]
        landscape = True
        theme = 'lgmerge'
        title = r"""Interaction energy [kcal/mol] {{err}} subset statistics with computed with {{opt}}{0}.""".format(
            '' if subjoin else r""" and {bas}""")
        return rowplan, columnplan, landscape, footnotes, title, theme


class DB4(Database):
    def __init__(self, pythonpath=None, loadfrompickle=False, path=None):
        """Initialize FourDatabases object from SuperDatabase"""
        Database.__init__(self, ['s22', 'nbc10', 'hbc6', 'hsg'], dbse='DB4',
                          pythonpath=pythonpath, loadfrompickle=loadfrompickle, path=path)

        # # load up data and definitions
        # self.load_qcdata_byproject('dft')
        # self.load_qcdata_byproject('pt2')
        # #self.load_qcdata_byproject('dhdft')
        # self.load_subsets()
        self.define_supersubsets()
        self.define_supermodelchems()

    def define_supersubsets(self):
        """

        """
        self.sset['tt'] = ['default', 'default', 'default', 'default']
        self.sset['hb'] = ['hb', None, 'default', 'hb']
        self.sset['mx'] = ['mx', 'mx', None, 'mx']
        self.sset['dd'] = ['dd', 'dd', None, 'dd']
        self.sset['mxdd'] = ['mxdd', 'default', None, 'mxdd']
        self.sset['pp'] = ['mxddpp', 'mxddpp', None, None]
        self.sset['np'] = ['mxddnp', 'mxddnp', None, 'mxdd']
        self.sset['tt-5min'] = ['default', '5min', '5min', 'default']
        self.sset['hb-5min'] = ['hb', None, '5min', 'hb']
        self.sset['mx-5min'] = ['mx', 'mx-5min', None, 'mx']
        self.sset['dd-5min'] = ['dd', 'dd-5min', None, 'dd']
        self.sset['mxdd-5min'] = ['mxdd', '5min', None, 'mxdd']
        self.sset['pp-5min'] = ['mxddpp', 'mxddpp-5min', None, None]
        self.sset['np-5min'] = ['mxddnp', 'mxddnp-5min', None, 'mxdd']

    # def benchmark(self):
    #     """Returns the model chemistry label for the database's benchmark."""
    #     return 'C2001BENCH'

    def define_supermodelchems(self):
        """

        """
        self.benchmark = 'C2011BENCH'
        self.mcs['C2010BENCH'] = ['S22A', 'NBC100', 'HBC60', 'HSG0']
        self.mcs['C2011BENCH'] = ['S22B', 'NBC10A', 'HBC6A', 'HSGA']

        self.mcs['CCSD-CP-adz'] = ['CCSD-CP-adz', 'CCSD-CP-hadz', 'CCSD-CP-adz', 'CCSD-CP-hadz']
        self.mcs['CCSD-CP-atz'] = ['CCSD-CP-atz', 'CCSD-CP-hatz', 'CCSD-CP-atz', 'CCSD-CP-hatz']
        self.mcs['CCSD-CP-adtz'] = ['CCSD-CP-adtz', 'CCSD-CP-hadtz', 'CCSD-CP-adtz', 'CCSD-CP-hadtz']
        self.mcs['CCSD-CP-adtzadz'] = ['CCSD-CP-adtzadz', 'CCSD-CP-adtzhadz', 'CCSD-CP-adtzadz', 'CCSD-CP-adtzhadz']
        self.mcs['CCSD-CP-atzadz'] = ['CCSD-CP-atzadz', 'CCSD-CP-atzhadz', 'CCSD-CP-atzadz', 'CCSD-CP-atzhadz']
        self.mcs['CCSD-CP-atqzadz'] = ['CCSD-CP-atqzadz', 'CCSD-CP-atqzhadz', 'CCSD-CP-atqzadz', 'CCSD-CP-atqzhadz']
        self.mcs['CCSD-CP-atzadtz'] = ['CCSD-CP-atzadtz', 'CCSD-CP-atzhadtz', 'CCSD-CP-atzadtz', 'CCSD-CP-atzhadtz']
        self.mcs['CCSD-CP-atqzadtz'] = ['CCSD-CP-atqzadtz', 'CCSD-CP-atqzhadtz', 'CCSD-CP-atqzadtz',
                                        'CCSD-CP-atqzhadtz']
        self.mcs['CCSD-CP-atqzatz'] = ['CCSD-CP-atqzatz', 'CCSD-CP-atqzhatz', 'CCSD-CP-atqzatz', 'CCSD-CP-atqzhatz']

        self.mcs['SCSCCSD-CP-adz'] = ['SCSCCSD-CP-adz', 'SCSCCSD-CP-hadz', 'SCSCCSD-CP-adz', 'SCSCCSD-CP-hadz']
        self.mcs['SCSCCSD-CP-atz'] = ['SCSCCSD-CP-atz', 'SCSCCSD-CP-hatz', 'SCSCCSD-CP-atz', 'SCSCCSD-CP-hatz']
        self.mcs['SCSCCSD-CP-adtz'] = ['SCSCCSD-CP-adtz', 'SCSCCSD-CP-hadtz', 'SCSCCSD-CP-adtz', 'SCSCCSD-CP-hadtz']
        self.mcs['SCSCCSD-CP-adtzadz'] = ['SCSCCSD-CP-adtzadz', 'SCSCCSD-CP-adtzhadz', 'SCSCCSD-CP-adtzadz',
                                          'SCSCCSD-CP-adtzhadz']
        self.mcs['SCSCCSD-CP-atzadz'] = ['SCSCCSD-CP-atzadz', 'SCSCCSD-CP-atzhadz', 'SCSCCSD-CP-atzadz',
                                         'SCSCCSD-CP-atzhadz']
        self.mcs['SCSCCSD-CP-atqzadz'] = ['SCSCCSD-CP-atqzadz', 'SCSCCSD-CP-atqzhadz', 'SCSCCSD-CP-atqzadz',
                                          'SCSCCSD-CP-atqzhadz']
        self.mcs['SCSCCSD-CP-atzadtz'] = ['SCSCCSD-CP-atzadtz', 'SCSCCSD-CP-atzhadtz', 'SCSCCSD-CP-atzadtz',
                                          'SCSCCSD-CP-atzhadtz']
        self.mcs['SCSCCSD-CP-atqzadtz'] = ['SCSCCSD-CP-atqzadtz', 'SCSCCSD-CP-atqzhadtz', 'SCSCCSD-CP-atqzadtz',
                                           'SCSCCSD-CP-atqzhadtz']
        self.mcs['SCSCCSD-CP-atqzatz'] = ['SCSCCSD-CP-atqzatz', 'SCSCCSD-CP-atqzhatz', 'SCSCCSD-CP-atqzatz',
                                          'SCSCCSD-CP-atqzhatz']

        self.mcs['SCSMICCSD-CP-adz'] = ['SCSMICCSD-CP-adz', 'SCSMICCSD-CP-hadz', 'SCSMICCSD-CP-adz',
                                        'SCSMICCSD-CP-hadz']
        self.mcs['SCSMICCSD-CP-atz'] = ['SCSMICCSD-CP-atz', 'SCSMICCSD-CP-hatz', 'SCSMICCSD-CP-atz',
                                        'SCSMICCSD-CP-hatz']
        self.mcs['SCSMICCSD-CP-adtz'] = ['SCSMICCSD-CP-adtz', 'SCSMICCSD-CP-hadtz', 'SCSMICCSD-CP-adtz',
                                         'SCSMICCSD-CP-hadtz']
        self.mcs['SCSMICCSD-CP-adtzadz'] = ['SCSMICCSD-CP-adtzadz', 'SCSMICCSD-CP-adtzhadz', 'SCSMICCSD-CP-adtzadz',
                                            'SCSMICCSD-CP-adtzhadz']
        self.mcs['SCSMICCSD-CP-atzadz'] = ['SCSMICCSD-CP-atzadz', 'SCSMICCSD-CP-atzhadz', 'SCSMICCSD-CP-atzadz',
                                           'SCSMICCSD-CP-atzhadz']
        self.mcs['SCSMICCSD-CP-atqzadz'] = ['SCSMICCSD-CP-atqzadz', 'SCSMICCSD-CP-atqzhadz', 'SCSMICCSD-CP-atqzadz',
                                            'SCSMICCSD-CP-atqzhadz']
        self.mcs['SCSMICCSD-CP-atzadtz'] = ['SCSMICCSD-CP-atzadtz', 'SCSMICCSD-CP-atzhadtz', 'SCSMICCSD-CP-atzadtz',
                                            'SCSMICCSD-CP-atzhadtz']
        self.mcs['SCSMICCSD-CP-atqzadtz'] = ['SCSMICCSD-CP-atqzadtz', 'SCSMICCSD-CP-atqzhadtz', 'SCSMICCSD-CP-atqzadtz',
                                             'SCSMICCSD-CP-atqzhadtz']
        self.mcs['SCSMICCSD-CP-atqzatz'] = ['SCSMICCSD-CP-atqzatz', 'SCSMICCSD-CP-atqzhatz', 'SCSMICCSD-CP-atqzatz',
                                            'SCSMICCSD-CP-atqzhatz']

        self.mcs['CCSDT-CP-adz'] = ['CCSDT-CP-adz', 'CCSDT-CP-hadz', 'CCSDT-CP-adz', 'CCSDT-CP-hadz']
        self.mcs['CCSDT-CP-atz'] = ['CCSDT-CP-atz', 'CCSDT-CP-hatz', 'CCSDT-CP-atz', 'CCSDT-CP-hatz']
        self.mcs['CCSDT-CP-adtz'] = ['CCSDT-CP-adtz', 'CCSDT-CP-hadtz', 'CCSDT-CP-adtz', 'CCSDT-CP-hadtz']
        self.mcs['CCSDT-CP-adtzadz'] = ['CCSDT-CP-adtzadz', 'CCSDT-CP-adtzhadz', 'CCSDT-CP-adtzadz',
                                        'CCSDT-CP-adtzhadz']
        self.mcs['CCSDT-CP-atzadz'] = ['CCSDT-CP-atzadz', 'CCSDT-CP-atzhadz', 'CCSDT-CP-atzadz', 'CCSDT-CP-atzhadz']
        self.mcs['CCSDT-CP-atqzadz'] = ['CCSDT-CP-atqzadz', 'CCSDT-CP-atqzhadz', 'CCSDT-CP-atqzadz',
                                        'CCSDT-CP-atqzhadz']
        self.mcs['CCSDT-CP-atzadtz'] = ['CCSDT-CP-atzadtz', 'CCSDT-CP-atzhadtz', 'CCSDT-CP-atzadtz',
                                        'CCSDT-CP-atzhadtz']
        self.mcs['CCSDT-CP-atqzadtz'] = ['CCSDT-CP-atqzadtz', 'CCSDT-CP-atqzhadtz', 'CCSDT-CP-atqzadtz',
                                         'CCSDT-CP-atqzhadtz']
        self.mcs['CCSDT-CP-atqzatz'] = ['CCSDT-CP-atqzatz', 'CCSDT-CP-atqzhatz', 'CCSDT-CP-atqzatz',
                                        'CCSDT-CP-atqzhatz']

        # def make_pt2_flats(self):
        # def plot_all_flats(self):
        #    """Generate pieces for inclusion into tables for PT2 paper.
        #    Note that DB4 flats use near-equilibrium subset.
        #
        #   """
        # Database.plot_all_flats(self, modelchem=None, sset='tt-5min', xlimit=4.0,
        #    graphicsformat=['pdf'])

    def make_pt2_Figure_3(self):
        """Plot all the graphics needed for the calendar grey bars plot
        in Fig. 3 of PT2.

        Note that in the modern implementation of class DB4, would need to
        pass ``sset=['tt-5min', 'hb-5min', 'mx-5min', 'dd-5min']`` to get
        published figure.

        """
        # Fig. bars (a)
        self.plot_bars(['MP2-CP-dz', 'MP2-CP-jadz', 'MP2-CP-hadz', 'MP2-CP-adz',
                        'MP2-CP-tz', 'MP2-CP-matz', 'MP2-CP-jatz', 'MP2-CP-hatz', 'MP2-CP-atz',
                        'MP2-CP-dtz', 'MP2-CP-jadtz', 'MP2-CP-hadtz', 'MP2-CP-adtz',
                        'MP2-CP-qz', 'MP2-CP-aaqz', 'MP2-CP-maqz', 'MP2-CP-jaqz', 'MP2-CP-haqz', 'MP2-CP-aqz',
                        'MP2-CP-tqz', 'MP2-CP-matqz', 'MP2-CP-jatqz', 'MP2-CP-hatqz', 'MP2-CP-atqz',
                        'MP2-CP-a5z', 'MP2-CP-aq5z'])
        self.plot_bars(['SCSMP2-CP-dz', 'SCSMP2-CP-jadz', 'SCSMP2-CP-hadz', 'SCSMP2-CP-adz',
                        'SCSMP2-CP-tz', 'SCSMP2-CP-matz', 'SCSMP2-CP-jatz', 'SCSMP2-CP-hatz', 'SCSMP2-CP-atz',
                        'SCSMP2-CP-dtz', 'SCSMP2-CP-jadtz', 'SCSMP2-CP-hadtz', 'SCSMP2-CP-adtz',
                        'SCSMP2-CP-qz', 'SCSMP2-CP-aaqz', 'SCSMP2-CP-maqz', 'SCSMP2-CP-jaqz', 'SCSMP2-CP-haqz',
                        'SCSMP2-CP-aqz',
                        'SCSMP2-CP-tqz', 'SCSMP2-CP-matqz', 'SCSMP2-CP-jatqz', 'SCSMP2-CP-hatqz', 'SCSMP2-CP-atqz',
                        'SCSMP2-CP-a5z', 'SCSMP2-CP-aq5z'])
        self.plot_bars(['SCSNMP2-CP-dz', 'SCSNMP2-CP-jadz', 'SCSNMP2-CP-hadz', 'SCSNMP2-CP-adz',
                        'SCSNMP2-CP-tz', 'SCSNMP2-CP-matz', 'SCSNMP2-CP-jatz', 'SCSNMP2-CP-hatz', 'SCSNMP2-CP-atz',
                        'SCSNMP2-CP-dtz', 'SCSNMP2-CP-jadtz', 'SCSNMP2-CP-hadtz', 'SCSNMP2-CP-adtz',
                        'SCSNMP2-CP-qz', 'SCSNMP2-CP-aaqz', 'SCSNMP2-CP-maqz', 'SCSNMP2-CP-jaqz', 'SCSNMP2-CP-haqz',
                        'SCSNMP2-CP-aqz',
                        'SCSNMP2-CP-tqz', 'SCSNMP2-CP-matqz', 'SCSNMP2-CP-jatqz', 'SCSNMP2-CP-hatqz', 'SCSNMP2-CP-atqz',
                        'SCSNMP2-CP-a5z', 'SCSNMP2-CP-aq5z'])
        self.plot_bars([None, None, None, None,
                        'SCSMIMP2-CP-tz', 'SCSMIMP2-CP-matz', 'SCSMIMP2-CP-jatz', 'SCSMIMP2-CP-hatz', 'SCSMIMP2-CP-atz',
                        'SCSMIMP2-CP-dtz', 'SCSMIMP2-CP-jadtz', 'SCSMIMP2-CP-hadtz', 'SCSMIMP2-CP-adtz',
                        'SCSMIMP2-CP-qz', 'SCSMIMP2-CP-aaqz', 'SCSMIMP2-CP-maqz', 'SCSMIMP2-CP-jaqz',
                        'SCSMIMP2-CP-haqz', 'SCSMIMP2-CP-aqz',
                        'SCSMIMP2-CP-tqz', 'SCSMIMP2-CP-matqz', 'SCSMIMP2-CP-jatqz', 'SCSMIMP2-CP-hatqz',
                        'SCSMIMP2-CP-atqz',
                        None, None])
        self.plot_bars(['DWMP2-CP-dz', 'DWMP2-CP-jadz', 'DWMP2-CP-hadz', 'DWMP2-CP-adz',
                        'DWMP2-CP-tz', 'DWMP2-CP-matz', 'DWMP2-CP-jatz', 'DWMP2-CP-hatz', 'DWMP2-CP-atz',
                        'DWMP2-CP-dtz', 'DWMP2-CP-jadtz', 'DWMP2-CP-hadtz', 'DWMP2-CP-adtz',
                        'DWMP2-CP-qz', 'DWMP2-CP-aaqz', 'DWMP2-CP-maqz', 'DWMP2-CP-jaqz', 'DWMP2-CP-haqz',
                        'DWMP2-CP-aqz',
                        'DWMP2-CP-tqz', 'DWMP2-CP-matqz', 'DWMP2-CP-jatqz', 'DWMP2-CP-hatqz', 'DWMP2-CP-atqz',
                        'DWMP2-CP-a5z', 'DWMP2-CP-aq5z'])
        self.plot_bars(['MP2C-CP-dz', 'MP2C-CP-jadz', 'MP2C-CP-hadz', 'MP2C-CP-adz',
                        'MP2C-CP-tz', 'MP2C-CP-matz', 'MP2C-CP-jatz', 'MP2C-CP-hatz', 'MP2C-CP-atz',
                        'MP2C-CP-dtz', 'MP2C-CP-jadtz', 'MP2C-CP-hadtz', 'MP2C-CP-adtz',
                        None, None, None, None, None, 'MP2C-CP-aqz',
                        None, None, None, None, 'MP2C-CP-atqz',
                        None, None])
        self.plot_bars(['MP2C-CP-atqzdz', 'MP2C-CP-atqzjadz', 'MP2C-CP-atqzhadz', 'MP2C-CP-atqzadz',
                        'MP2C-CP-atqztz', 'MP2C-CP-atqzmatz', 'MP2C-CP-atqzjatz', 'MP2C-CP-atqzhatz', 'MP2C-CP-atqzatz',
                        'MP2C-CP-atqzdtz', 'MP2C-CP-atqzjadtz', 'MP2C-CP-atqzhadtz', 'MP2C-CP-atqzadtz'])

        # Fig. bars (c)
        self.plot_bars(['MP2F12-CP-dz', 'MP2F12-CP-jadz', 'MP2F12-CP-hadz', 'MP2F12-CP-adz',
                        'MP2F12-CP-tz', 'MP2F12-CP-matz', 'MP2F12-CP-jatz', 'MP2F12-CP-hatz', 'MP2F12-CP-atz',
                        'MP2F12-CP-dtz', 'MP2F12-CP-jadtz', 'MP2F12-CP-hadtz', 'MP2F12-CP-adtz',
                        'MP2F12-CP-aqz', 'MP2F12-CP-atqz'])
        self.plot_bars(['SCSMP2F12-CP-dz', 'SCSMP2F12-CP-jadz', 'SCSMP2F12-CP-hadz', 'SCSMP2F12-CP-adz',
                        'SCSMP2F12-CP-tz', 'SCSMP2F12-CP-matz', 'SCSMP2F12-CP-jatz', 'SCSMP2F12-CP-hatz',
                        'SCSMP2F12-CP-atz',
                        'SCSMP2F12-CP-dtz', 'SCSMP2F12-CP-jadtz', 'SCSMP2F12-CP-hadtz', 'SCSMP2F12-CP-adtz',
                        'SCSMP2F12-CP-aqz', 'SCSMP2F12-CP-atqz'])
        self.plot_bars(['SCSNMP2F12-CP-dz', 'SCSNMP2F12-CP-jadz', 'SCSNMP2F12-CP-hadz', 'SCSNMP2F12-CP-adz',
                        'SCSNMP2F12-CP-tz', 'SCSNMP2F12-CP-matz', 'SCSNMP2F12-CP-jatz', 'SCSNMP2F12-CP-hatz',
                        'SCSNMP2F12-CP-atz',
                        'SCSNMP2F12-CP-dtz', 'SCSNMP2F12-CP-jadtz', 'SCSNMP2F12-CP-adtz', 'SCSNMP2F12-CP-adtz',
                        'SCSNMP2F12-CP-aqz', 'SCSNMP2F12-CP-atqz'])
        self.plot_bars([None, None, None, None,
                        'SCSMIMP2F12-CP-tz', 'SCSMIMP2F12-CP-matz', 'SCSMIMP2F12-CP-jatz', 'SCSMIMP2F12-CP-hatz',
                        'SCSMIMP2F12-CP-atz',
                        'SCSMIMP2F12-CP-dtz', 'SCSMIMP2F12-CP-jadtz', 'SCSMIMP2F12-CP-hadtz', 'SCSMIMP2F12-CP-adtz',
                        'SCSMIMP2F12-CP-aqz', 'SCSMIMP2F12-CP-atqz'])
        self.plot_bars(['DWMP2F12-CP-dz', 'DWMP2F12-CP-jadz', 'DWMP2F12-CP-hadz', 'DWMP2F12-CP-adz',
                        'DWMP2F12-CP-tz', 'DWMP2F12-CP-matz', 'DWMP2F12-CP-jatz', 'DWMP2F12-CP-hatz', 'DWMP2F12-CP-atz',
                        'DWMP2F12-CP-dtz', 'DWMP2F12-CP-jadtz', 'DWMP2F12-CP-hadtz', 'DWMP2F12-CP-adtz',
                        'DWMP2F12-CP-aqz', 'DWMP2F12-CP-atqz'])
        self.plot_bars(['MP2CF12-CP-dz', 'MP2CF12-CP-jadz', 'MP2CF12-CP-hadz', 'MP2CF12-CP-adz',
                        'MP2CF12-CP-tz', 'MP2CF12-CP-matz', 'MP2CF12-CP-jatz', 'MP2CF12-CP-hatz', 'MP2CF12-CP-atz',
                        'MP2CF12-CP-dtz', 'MP2CF12-CP-jadtz', 'MP2CF12-CP-hadtz', 'MP2CF12-CP-adtz',
                        'MP2CF12-CP-aqz', 'MP2CF12-CP-atqz'])
        self.plot_bars(['MP2CF12-CP-atqzdz', 'MP2CF12-CP-atqzjadz', 'MP2CF12-CP-atqzhadz', 'MP2CF12-CP-atqzadz',
                        'MP2CF12-CP-atqztz', 'MP2CF12-CP-atqzmatz', 'MP2CF12-CP-atqzjatz', 'MP2CF12-CP-atqzhatz',
                        'MP2CF12-CP-atqzatz',
                        'MP2CF12-CP-atqzdtz', 'MP2CF12-CP-atqzjadtz', 'MP2CF12-CP-atqzhadtz', 'MP2CF12-CP-atqzadtz'])

    def make_pt2_Figure_2(self):
        """Plot all the graphics needed for the diffuse augmented grey
        bars plot in Fig. 2 of PT2.

        Note that in the modern implementation of class DB4, would need to
        pass ``sset=['tt-5min', 'hb-5min', 'mx-5min', 'dd-5min']`` to get
        published figure.

        """
        # Fig. bars (a)
        self.plot_bars(['MP2-CP-adz', 'MP2-CP-atz', 'MP2-CP-adtz',
                        'MP2-CP-aqz', 'MP2-CP-atqz', 'MP2-CP-a5z', 'MP2-CP-aq5z'])
        self.plot_bars(['SCSMP2-CP-adz', 'SCSMP2-CP-atz',
                        'SCSMP2-CP-adtz', 'SCSMP2-CP-aqz', 'SCSMP2-CP-atqz',
                        'SCSMP2-CP-a5z', 'SCSMP2-CP-aq5z'])
        self.plot_bars(['SCSNMP2-CP-adz', 'SCSNMP2-CP-atz',
                        'SCSNMP2-CP-adtz', 'SCSNMP2-CP-aqz', 'SCSNMP2-CP-atqz',
                        'SCSNMP2-CP-a5z', 'SCSNMP2-CP-aq5z'])
        self.plot_bars(['SCSMIMP2-CP-atz', 'SCSMIMP2-CP-atz',
                        'SCSMIMP2-CP-adtz', 'SCSMIMP2-CP-aqz', 'SCSMIMP2-CP-atqz'])
        self.plot_bars(['SCSMIMP2-CP-tz', 'SCSMIMP2-CP-tz',
                        'SCSMIMP2-CP-dtz', 'SCSMIMP2-CP-qz', 'SCSMIMP2-CP-tqz'])
        self.plot_bars(['DWMP2-CP-adz', 'DWMP2-CP-atz', 'DWMP2-CP-adtz',
                        'DWMP2-CP-aqz', 'DWMP2-CP-atqz', 'DWMP2-CP-a5z', 'DWMP2-CP-aq5z'])
        self.plot_bars(['MP2C-CP-adz', 'MP2C-CP-adtzadz',
                        'MP2C-CP-atqzadz', 'MP2C-CP-aq5zadz', 'MP2C-CP-atz',
                        'MP2C-CP-atqzatz', 'MP2C-CP-aq5zatz', 'MP2C-CP-adtz',
                        'MP2C-CP-atqzadtz', 'MP2C-CP-aqz', 'MP2C-CP-atqz'])

        # Fig. bars (b)
        self.plot_bars(['MP3-CP-adz', 'MP3-CP-adtzadz', 'MP3-CP-atqzadz',
                        'MP3-CP-atz', 'MP3-CP-atqzatz', 'MP3-CP-adtz', 'MP3-CP-atqzadtz'])
        self.plot_bars(['MP25-CP-adz', 'MP25-CP-adtzadz', 'MP25-CP-atqzadz',
                        'MP25-CP-atz', 'MP25-CP-atqzatz', 'MP25-CP-adtz', 'MP25-CP-atqzadtz'])
        self.plot_bars(['CCSD-CP-adz', 'CCSD-CP-adtzadz', 'CCSD-CP-atqzadz',
                        'CCSD-CP-atz', 'CCSD-CP-atqzatz', 'CCSD-CP-adtz', 'CCSD-CP-atqzadtz'])
        self.plot_bars(['SCSCCSD-CP-adz', 'SCSCCSD-CP-adtzadz',
                        'SCSCCSD-CP-atqzadz', 'SCSCCSD-CP-atz', 'SCSCCSD-CP-atqzatz',
                        'SCSCCSD-CP-adtz', 'SCSCCSD-CP-atqzadtz'])
        self.plot_bars(['SCSMICCSD-CP-adz', 'SCSMICCSD-CP-adtzadz',
                        'SCSMICCSD-CP-atqzadz', 'SCSMICCSD-CP-atz', 'SCSMICCSD-CP-atqzatz',
                        'SCSMICCSD-CP-adtz', 'SCSMICCSD-CP-atqzadtz'])
        self.plot_bars(['CCSDT-CP-adz', 'CCSDT-CP-adtzadz',
                        'CCSDT-CP-atqzadz', 'CCSDT-CP-atz', 'CCSDT-CP-atqzatz',
                        'CCSDT-CP-adtz', 'CCSDT-CP-atqzadtz'])

        # Fig. bars (c)
        self.plot_bars(['MP2F12-CP-adz', 'MP2F12-CP-atz', 'MP2F12-CP-adtz',
                        'MP2F12-CP-aqz', 'MP2F12-CP-atqz'])
        self.plot_bars(['SCSMP2F12-CP-adz', 'SCSMP2F12-CP-atz',
                        'SCSMP2F12-CP-adtz', 'SCSMP2F12-CP-aqz', 'SCSMP2F12-CP-atqz'])
        self.plot_bars(['SCSNMP2F12-CP-adz', 'SCSNMP2F12-CP-atz',
                        'SCSNMP2F12-CP-adtz', 'SCSNMP2F12-CP-aqz',
                        'SCSNMP2F12-CP-atqz'])
        self.plot_bars(['SCSMIMP2F12-CP-atz', 'SCSMIMP2F12-CP-atz',
                        'SCSMIMP2F12-CP-adtz', 'SCSMIMP2F12-CP-aqz',
                        'SCSMIMP2F12-CP-atqz'])
        self.plot_bars(['SCSMIMP2F12-CP-tz', 'SCSMIMP2F12-CP-tz', 'SCSMIMP2F12-CP-dtz'])
        self.plot_bars(['DWMP2F12-CP-adz', 'DWMP2F12-CP-atz',
                        'DWMP2F12-CP-adtz', 'DWMP2F12-CP-aqz', 'DWMP2F12-CP-atqz'])
        self.plot_bars(['MP2CF12-CP-adz', 'MP2CF12-CP-adtzadz',
                        'MP2CF12-CP-atqzadz', 'MP2CF12-CP-atz', 'MP2CF12-CP-atqzatz',
                        'MP2CF12-CP-adtz', 'MP2CF12-CP-atqzadtz', 'MP2CF12-CP-aqz',
                        'MP2CF12-CP-atqz'])

        # Fig. bars (d)
        self.plot_bars(['CCSDAF12-CP-adz', 'CCSDAF12-CP-adtzadz', 'CCSDAF12-CP-atqzadz'])
        self.plot_bars(['CCSDBF12-CP-adz', 'CCSDBF12-CP-adtzadz', 'CCSDBF12-CP-atqzadz'])
        self.plot_bars(['SCSCCSDAF12-CP-adz', 'SCSCCSDAF12-CP-adtzadz', 'SCSCCSDAF12-CP-atqzadz'])
        self.plot_bars(['SCSCCSDBF12-CP-adz', 'SCSCCSDBF12-CP-adtzadz', 'SCSCCSDBF12-CP-atqzadz'])
        self.plot_bars(['SCMICCSDAF12-CP-adz', 'SCMICCSDAF12-CP-adtzadz', 'SCMICCSDAF12-CP-atqzadz'])
        self.plot_bars(['SCMICCSDBF12-CP-adz', 'SCMICCSDBF12-CP-adtzadz', 'SCMICCSDBF12-CP-atqzadz'])
        self.plot_bars(['CCSDTAF12-CP-adz', 'CCSDTAF12-CP-adtzadz', 'CCSDTAF12-CP-atqzadz'])
        self.plot_bars(['CCSDTBF12-CP-adz', 'CCSDTBF12-CP-adtzadz', 'CCSDTBF12-CP-atqzadz'])
        self.plot_bars(['DWCCSDTF12-CP-adz', 'DWCCSDTF12-CP-adtzadz', 'DWCCSDTF12-CP-atqzadz'])

    def plot_dhdft_flats(self):
        """Generate pieces for grey bars figure for DH-DFT paper."""

        self.plot_all_flats(
            ['B97D3-CP-adz', 'PBED3-CP-adz', 'M11L-CP-adz', 'DLDFD-CP-adz', 'B3LYPD3-CP-adz', 'PBE0D3-CP-adz',
             'WB97XD-CP-adz', 'M052X-CP-adz', 'M062X-CP-adz', 'M08HX-CP-adz', 'M08SO-CP-adz', 'M11-CP-adz',
             'VV10-CP-adz',
             'LCVV10-CP-adz', 'WB97XV-CP-adz', 'PBE02-CP-adz', 'WB97X2-CP-adz', 'B2PLYPD3-CP-adz',
             'DSDPBEP86D2OPT-CP-adz', 'MP2-CP-adz'], sset='tt-5min')
        self.plot_all_flats(['B97D3-unCP-adz', 'PBED3-unCP-adz', 'M11L-unCP-adz', 'DLDFD-unCP-adz', 'B3LYPD3-unCP-adz',
                             'PBE0D3-unCP-adz',
                             'WB97XD-unCP-adz', 'M052X-unCP-adz', 'M062X-unCP-adz', 'M08HX-unCP-adz', 'M08SO-unCP-adz',
                             'M11-unCP-adz', 'VV10-unCP-adz',
                             'LCVV10-unCP-adz', 'WB97XV-unCP-adz', 'PBE02-unCP-adz', 'WB97X2-unCP-adz',
                             'B2PLYPD3-unCP-adz', 'DSDPBEP86D2OPT-unCP-adz', 'MP2-unCP-adz'], sset='tt-5min')
        self.plot_all_flats(
            ['B97D3-CP-atz', 'PBED3-CP-atz', 'M11L-CP-atz', 'DLDFD-CP-atz', 'B3LYPD3-CP-atz', 'PBE0D3-CP-atz',
             'WB97XD-CP-atz', 'M052X-CP-atz', 'M062X-CP-atz', 'M08HX-CP-atz', 'M08SO-CP-atz', 'M11-CP-atz',
             'VV10-CP-atz',
             'LCVV10-CP-atz', 'WB97XV-CP-atz', 'PBE02-CP-atz', 'WB97X2-CP-atz', 'B2PLYPD3-CP-atz',
             'DSDPBEP86D2OPT-CP-atz', 'MP2-CP-atz'], sset='tt-5min')
        self.plot_all_flats(['B97D3-unCP-atz', 'PBED3-unCP-atz', 'M11L-unCP-atz', 'DLDFD-unCP-atz', 'B3LYPD3-unCP-atz',
                             'PBE0D3-unCP-atz',
                             'WB97XD-unCP-atz', 'M052X-unCP-atz', 'M062X-unCP-atz', 'M08HX-unCP-atz', 'M08SO-unCP-atz',
                             'M11-unCP-atz', 'VV10-unCP-atz',
                             'LCVV10-unCP-atz', 'WB97XV-unCP-atz', 'PBE02-unCP-atz', 'WB97X2-unCP-atz',
                             'B2PLYPD3-unCP-atz', 'DSDPBEP86D2OPT-unCP-atz', 'MP2-unCP-atz'], sset='tt-5min')

    def make_dhdft_Figure_1(self):
        """Plot all the graphics needed for the grey bars plot
        in Fig. 1 of DHDFT.

        """
        # Fig. bars (a)
        self.plot_bars([
            'M052X-unCP-adz', 'M052X-CP-adz', 'M052X-unCP-atz', 'M052X-CP-atz', None,
            'M062X-unCP-adz', 'M062X-CP-adz', 'M062X-unCP-atz', 'M062X-CP-atz', None,
            'M08SO-unCP-adz', 'M08SO-CP-adz', 'M08SO-unCP-atz', 'M08SO-CP-atz', None,
            'M08HX-unCP-adz', 'M08HX-CP-adz', 'M08HX-unCP-atz', 'M08HX-CP-atz', None,
            'M11-unCP-adz', 'M11-CP-adz', 'M11-unCP-atz', 'M11-CP-atz', None,
            'M11L-unCP-adz', 'M11L-CP-adz', 'M11L-unCP-atz', 'M11L-CP-atz'],
            sset=['tt-5min', 'hb-5min', 'mx-5min', 'dd-5min'])

        # Fig. bars (b)
        self.plot_bars([
            'PBED3-unCP-adz', 'PBED3-CP-adz', 'PBED3-unCP-atz', 'PBED3-CP-atz', None,
            'B97D3-unCP-adz', 'B97D3-CP-adz', 'B97D3-unCP-atz', 'B97D3-CP-atz', None,
            'PBE0D3-unCP-adz', 'PBE0D3-CP-adz', 'PBE0D3-unCP-atz', 'PBE0D3-CP-atz', None,
            'B3LYPD3-unCP-adz', 'B3LYPD3-CP-adz', 'B3LYPD3-unCP-atz', 'B3LYPD3-CP-atz', None,
            'DLDFD-unCP-adz', 'DLDFD-CP-adz', 'DLDFD-unCP-atz', 'DLDFD-CP-atz', None,
            'WB97XD-unCP-adz', 'WB97XD-CP-adz', 'WB97XD-unCP-atz', 'WB97XD-CP-atz'],
            sset=['tt-5min', 'hb-5min', 'mx-5min', 'dd-5min'])

        # Fig. bars (c)
        self.plot_bars([
            'VV10-unCP-adz', 'VV10-CP-adz', 'VV10-unCP-atz', 'VV10-CP-atz', None, None,
            'LCVV10-unCP-adz', 'LCVV10-CP-adz', 'LCVV10-unCP-atz', 'LCVV10-CP-atz', None, None,
            'WB97XV-unCP-adz', 'WB97XV-CP-adz', 'WB97XV-unCP-atz', 'WB97XV-CP-atz'],
            sset=['tt-5min', 'hb-5min', 'mx-5min', 'dd-5min'])

        # Fig. bars (d)
        self.plot_bars([
            'PBE02-unCP-adz', 'PBE02-CP-adz', 'PBE02-unCP-atz', 'PBE02-CP-atz', None,
            'WB97X2-unCP-adz', 'WB97X2-CP-adz', 'WB97X2-unCP-atz', 'WB97X2-CP-atz', None,
            'B2PLYPD3-unCP-adz', 'B2PLYPD3-CP-adz', 'B2PLYPD3-unCP-atz', 'B2PLYPD3-CP-atz', None,
            'DSDPBEP86D2OPT-unCP-adz', 'DSDPBEP86D2OPT-CP-adz', 'DSDPBEP86D2OPT-unCP-atz', 'DSDPBEP86D2OPT-CP-atz'],
            sset=['tt-5min', 'hb-5min', 'mx-5min', 'dd-5min'])

        # Fig. bars (e)
        self.plot_bars([
            'MP2-unCP-adz', 'MP2-CP-adz', 'MP2-unCP-atz', 'MP2-CP-atz'],
            sset=['tt-5min', 'hb-5min', 'mx-5min', 'dd-5min'])

    def make_dhdft_Figure_2(self):
        """Plot all the graphics needed for the SAPT/DFT/WFN
        comparison plot in Fig. 2 of DHDFT.

        Note that benchmark set as reminder, not necessity, since default.

        """
        self.plot_bars([
            'SAPT0S-CP-jadz', 'SAPTDFT-CP-atz', 'SAPT2P-CP-adz', 'SAPT3M-CP-atz',
            'SAPT2PCM-CP-atz', None, 'B97D3-unCP-atz', 'B3LYPD3-CP-adz',
            'M052X-unCP-adz', 'WB97XD-CP-atz', 'WB97XV-CP-adz', 'WB97X2-CP-atz',
            'DSDPBEP86D2OPT-CP-atz', 'B2PLYPD3-CP-atz', None, 'MP2-CP-atz',
            'SCSMP2-CP-atz', 'SCSMIMP2-CP-qz', 'MP2C-CP-atqzadz',
            'MP2CF12-CP-adz', 'SCMICCSDAF12-CP-adz', 'CCSDT-CP-atz',
            'CCSDT-CP-atqzatz', 'DWCCSDTF12-CP-adz'],
            sset=['tt-5min', 'hb-5min', 'mx-5min', 'dd-5min'],
            benchmark='C2011BENCH')

    def plot_dhdft_modelchems(self):
        self.plot_modelchems(
            ['B97D3-CP-adz', 'PBED3-CP-adz', 'M11L-CP-adz', 'DLDFD-CP-adz', 'B3LYPD3-CP-adz', 'PBE0D3-CP-adz',
             'WB97XD-CP-adz', 'M052X-CP-adz', 'M062X-CP-adz', 'M08HX-CP-adz', 'M08SO-CP-adz', 'M11-CP-adz',
             'VV10-CP-adz',
             'LCVV10-CP-adz', 'WB97XV-CP-adz', 'PBE02-CP-adz', 'WB97X2-CP-adz', 'B2PLYPD3-CP-adz',
             'DSDPBEP86D2OPT-CP-adz', 'MP2-CP-adz'], sset='tt-5min')
        self.plot_modelchems(['B97D3-unCP-adz', 'PBED3-unCP-adz', 'M11L-unCP-adz', 'DLDFD-unCP-adz', 'B3LYPD3-unCP-adz',
                              'PBE0D3-unCP-adz',
                              'WB97XD-unCP-adz', 'M052X-unCP-adz', 'M062X-unCP-adz', 'M08HX-unCP-adz', 'M08SO-unCP-adz',
                              'M11-unCP-adz', 'VV10-unCP-adz',
                              'LCVV10-unCP-adz', 'WB97XV-unCP-adz', 'PBE02-unCP-adz', 'WB97X2-unCP-adz',
                              'B2PLYPD3-unCP-adz', 'DSDPBEP86D2OPT-unCP-adz', 'MP2-unCP-adz'], sset='tt-5min')
        self.plot_modelchems(
            ['B97D3-CP-atz', 'PBED3-CP-atz', 'M11L-CP-atz', 'DLDFD-CP-atz', 'B3LYPD3-CP-atz', 'PBE0D3-CP-atz',
             'WB97XD-CP-atz', 'M052X-CP-atz', 'M062X-CP-atz', 'M08HX-CP-atz', 'M08SO-CP-atz', 'M11-CP-atz',
             'VV10-CP-atz',
             'LCVV10-CP-atz', 'WB97XV-CP-atz', 'PBE02-CP-atz', 'WB97X2-CP-atz', 'B2PLYPD3-CP-atz',
             'DSDPBEP86D2OPT-CP-atz', 'MP2-CP-atz'], sset='tt-5min')
        self.plot_modelchems(['B97D3-unCP-atz', 'PBED3-unCP-atz', 'M11L-unCP-atz', 'DLDFD-unCP-atz', 'B3LYPD3-unCP-atz',
                              'PBE0D3-unCP-atz',
                              'WB97XD-unCP-atz', 'M052X-unCP-atz', 'M062X-unCP-atz', 'M08HX-unCP-atz', 'M08SO-unCP-atz',
                              'M11-unCP-atz', 'VV10-unCP-atz',
                              'LCVV10-unCP-atz', 'WB97XV-unCP-atz', 'PBE02-unCP-atz', 'WB97X2-unCP-atz',
                              'B2PLYPD3-unCP-atz', 'DSDPBEP86D2OPT-unCP-atz', 'MP2-unCP-atz'], sset='tt-5min')

    def plot_minn_modelchems(self):
        self.plot_modelchems(
            ['DLDFD-unCP-adz', 'M052X-unCP-adz', 'M062X-unCP-adz', 'M08HX-unCP-adz', 'M08SO-unCP-adz', 'M11-unCP-adz',
             'M11L-unCP-adz',
             'DLDFD-CP-adz', 'M052X-CP-adz', 'M062X-CP-adz', 'M08HX-CP-adz', 'M08SO-CP-adz', 'M11-CP-adz',
             'M11L-CP-adz'])
        self.plot_modelchems(
            ['DlDFD-unCP-atz', 'M052X-unCP-atz', 'M062X-unCP-atz', 'M08HX-unCP-atz', 'M08SO-unCP-atz', 'M11-unCP-atz',
             'M11L-unCP-atz',
             'DLDFD-CP-atz', 'M052X-CP-atz', 'M062X-CP-atz', 'M08HX-CP-atz', 'M08SO-CP-atz', 'M11-CP-atz',
             'M11L-CP-atz'])

    def make_dhdft_Table_I(self):
        """Generate the in-manuscript summary slat table for DHDFT.

        """
        self.table_wrapper(mtd=['B97D3', 'PBED3', 'M11L', 'DLDFD', 'B3LYPD3',
                           'PBE0D3', 'WB97XD', 'M052X', 'M062X', 'M08HX',
                           'M08SO', 'M11', 'VV10', 'LCVV10', 'WB97XV',
                           'PBE02', 'WB97X2', 'DSDPBEP86D2OPT', 'B2PLYPD3',
                           'MP2', 'SCSNMP2', 'SCSMIMP2', 'MP2CF12', 'SCMICCSDAF12',
                           'SAPTDFT', 'SAPT0S', 'SAPT2P', 'SAPT3M', 'SAPT2PCM'],
                           bas=['adz', 'atz'],
                           tableplan=self.table_scrunch,
                           opt=['CP', 'unCP'], err=['mae'],
                           subjoin=None,
                           plotpath=None,
                           standalone=False, filename='tblssets_ex1')

    def make_dhdft_Table_II(self):
        """Generate the in-manuscript CP slat table for DHDFT.

        """
        self.table_wrapper(mtd=['B97D3', 'PBED3', 'M11L', 'DLDFD', 'B3LYPD3',
                           'PBE0D3', 'WB97XD', 'M052X', 'M062X', 'M08HX',
                           'M08SO', 'M11', 'VV10', 'LCVV10', 'WB97XV',
                           'PBE02', 'WB97X2', 'DSDPBEP86D2OPT', 'B2PLYPD3', 'MP2'],
                           bas=['adz', 'atz'],
                           tableplan=self.table_merge_abbr,
                           opt=['CP'], err=['mae'],
                           subjoin=True,
                           plotpath='analysis/flats/mplflat_',  # proj still has 'mpl' prefix
                           standalone=False, filename='tblssets_ex2')

    def make_dhdft_Table_III(self):
        """Generate the in-manuscript unCP slat table for DHDFT.

        """
        self.table_wrapper(mtd=['B97D3', 'PBED3', 'M11L', 'DLDFD', 'B3LYPD3',
                           'PBE0D3', 'WB97XD', 'M052X', 'M062X', 'M08HX',
                           'M08SO', 'M11', 'VV10', 'LCVV10', 'WB97XV',
                           'PBE02', 'WB97X2', 'DSDPBEP86D2OPT', 'B2PLYPD3', 'MP2'],
                           bas=['adz', 'atz'],
                           tableplan=self.table_merge_abbr,
                           opt=['unCP'], err=['mae'],
                           subjoin=True,
                           plotpath='analysis/flats/mplflat_',  # proj still has 'mpl' prefix
                           standalone=False, filename='tblssets_ex3')

    def make_dhdft_Tables_SII(self):
        """Generate the subset details suppmat Part II tables and their indices for DHDFT.

        """
        self.table_wrapper(mtd=['B97D3', 'PBED3', 'M11L', 'DLDFD', 'B3LYPD3',
                                'PBE0D3', 'WB97XD', 'M052X', 'M062X', 'M08HX',
                                'M08SO', 'M11', 'VV10', 'LCVV10', 'WB97XV',
                                'PBE02', 'WB97X2', 'DSDPBEP86D2OPT', 'B2PLYPD3'],  # 'MP2']
                           bas=['adz', 'atz'],
                           tableplan=self.table_merge_suppmat,
                           opt=['CP', 'unCP'], err=['mae', 'mape'],
                           subjoin=False,
                           plotpath='analysis/flats/mplflat_',  # proj still has 'mpl' prefix
                           standalone=False, filename='tblssets')

    def make_dhdft_Tables_SIII(self):
        """Generate the per-reaction suppmat Part III tables and their indices for DHDFT.

        """
        self.table_reactions(
            ['B97D3-unCP-adz', 'B97D3-CP-adz', 'B97D3-unCP-atz', 'B97D3-CP-atz',
             'PBED3-unCP-adz', 'PBED3-CP-adz', 'PBED3-unCP-atz', 'PBED3-CP-atz',
             'M11L-unCP-adz', 'M11L-CP-adz', 'M11L-unCP-atz', 'M11L-CP-atz',
             'DLDFD-unCP-adz', 'DLDFD-CP-adz', 'DLDFD-unCP-atz', 'DLDFD-CP-atz',
             'B3LYPD3-unCP-adz', 'B3LYPD3-CP-adz', 'B3LYPD3-unCP-atz', 'B3LYPD3-CP-atz',
             'PBE0D3-unCP-adz', 'PBE0D3-CP-adz', 'PBE0D3-unCP-atz', 'PBE0D3-CP-atz',
             'WB97XD-unCP-adz', 'WB97XD-CP-adz', 'WB97XD-unCP-atz', 'WB97XD-CP-atz',
             'M052X-unCP-adz', 'M052X-CP-adz', 'M052X-unCP-atz', 'M052X-CP-atz',
             'M062X-unCP-adz', 'M062X-CP-adz', 'M062X-unCP-atz', 'M062X-CP-atz',
             'M08HX-unCP-adz', 'M08HX-CP-adz', 'M08HX-unCP-atz', 'M08HX-CP-atz',
             'M08SO-unCP-adz', 'M08SO-CP-adz', 'M08SO-unCP-atz', 'M08SO-CP-atz',
             'M11-unCP-adz', 'M11-CP-adz', 'M11-unCP-atz', 'M11-CP-atz',
             'VV10-unCP-adz', 'VV10-CP-adz', 'VV10-unCP-atz', 'VV10-CP-atz',
             'LCVV10-unCP-adz', 'LCVV10-CP-adz', 'LCVV10-unCP-atz', 'LCVV10-CP-atz',
             'WB97XV-unCP-adz', 'WB97XV-CP-adz', 'WB97XV-unCP-atz', 'WB97XV-CP-atz',
             'PBE02-unCP-adz', 'PBE02-CP-adz', 'PBE02-unCP-atz', 'PBE02-CP-atz',
             'WB97X2-unCP-adz', 'WB97X2-CP-adz', 'WB97X2-unCP-atz', 'WB97X2-CP-atz',
             'DSDPBEP86D2OPT-unCP-adz', 'DSDPBEP86D2OPT-CP-adz', 'DSDPBEP86D2OPT-unCP-atz', 'DSDPBEP86D2OPT-CP-atz',
             'B2PLYPD3-unCP-adz', 'B2PLYPD3-CP-adz', 'B2PLYPD3-unCP-atz', 'B2PLYPD3-CP-atz'],
            # 'MP2-unCP-adz', 'MP2-CP-adz', 'MP2-unCP-atz', 'MP2-CP-atz'],
            standalone=False, filename='tblrxn_all')


class ThreeDatabases(Database):
    """

    """

    def __init__(self, pythonpath=None):
        """Initialize ThreeDatabases object from Database"""
        Database.__init__(self, ['s22', 'a24', 'hsg'], dbse='DB3', pythonpath=None)

        # load up data and definitions
        self.load_qcdata_byproject('pt2')
        self.load_qcdata_byproject('dilabio')
        self.load_qcdata_byproject('f12dilabio')
        self.load_subsets()
        self.define_supersubsets()
        self.define_supermodelchems()

    def define_supersubsets(self):
        """

        """
        self.sset['tt'] = ['default', 'default', 'default']
        self.sset['hb'] = ['hb', 'hb', 'hb']
        self.sset['mx'] = ['mx', 'mx', 'mx']
        self.sset['dd'] = ['dd', 'dd', 'dd']
        self.sset['mxdd'] = ['mxdd', 'mxdd', 'mxdd']
        self.sset['pp'] = ['mxddpp', 'mxddpp', 'mxddpp']
        self.sset['np'] = ['mxddnp', 'mxddnp', 'mxddnp']
        self.sset['tt-5min'] = ['default', 'default', 'default']
        self.sset['hb-5min'] = ['hb', 'hb', 'hb']
        self.sset['mx-5min'] = ['mx', 'mx', 'mx']
        self.sset['dd-5min'] = ['dd', 'dd', 'dd']
        self.sset['mxdd-5min'] = ['mxdd', 'mxdd', 'mxdd']
        self.sset['pp-5min'] = ['mxddpp', 'mxddpp', 'mxddpp']
        self.sset['np-5min'] = ['mxddnp', 'mxddnp', 'mxddnp']
        self.sset['weak'] = ['weak', 'weak', 'weak']
        self.sset['weak_hb'] = ['weak_hb', None, 'weak_hb']
        self.sset['weak_mx'] = ['weak_mx', 'weak_mx', 'weak_mx']
        self.sset['weak_dd'] = ['weak_dd', 'weak_dd', 'weak_dd']

    def define_supermodelchems(self):
        """

        """
        self.mc['CCSD-CP-adz'] = ['CCSD-CP-adz', 'CCSD-CP-hadz', 'CCSD-CP-adz']
        self.mc['CCSD-CP-atz'] = ['CCSD-CP-atz', 'CCSD-CP-hatz', 'CCSD-CP-atz']
        self.mc['CCSD-CP-adtz'] = ['CCSD-CP-adtz', 'CCSD-CP-hadtz', 'CCSD-CP-adtz']
        self.mc['CCSD-CP-adtzadz'] = ['CCSD-CP-adtzadz', 'CCSD-CP-adtzhadz', 'CCSD-CP-adtzadz']
        self.mc['CCSD-CP-atzadz'] = ['CCSD-CP-atzadz', 'CCSD-CP-atzhadz', 'CCSD-CP-atzadz']
        self.mc['CCSD-CP-atqzadz'] = ['CCSD-CP-atqzadz', 'CCSD-CP-atqzhadz', 'CCSD-CP-atqzadz']
        self.mc['CCSD-CP-atzadtz'] = ['CCSD-CP-atzadtz', 'CCSD-CP-atzhadtz', 'CCSD-CP-atzadtz']
        self.mc['CCSD-CP-atqzadtz'] = ['CCSD-CP-atqzadtz', 'CCSD-CP-atqzhadtz', 'CCSD-CP-atqzadtz']
        self.mc['CCSD-CP-atqzatz'] = ['CCSD-CP-atqzatz', 'CCSD-CP-atqzhatz', 'CCSD-CP-atqzatz']

        self.mc['CCSDT-CP-adz'] = ['CCSDT-CP-adz', 'CCSDT-CP-hadz', 'CCSDT-CP-adz']
        self.mc['CCSDT-CP-atz'] = ['CCSDT-CP-atz', 'CCSDT-CP-hatz', 'CCSDT-CP-atz']
        self.mc['CCSDT-CP-adtz'] = ['CCSDT-CP-adtz', 'CCSDT-CP-hadtz', 'CCSDT-CP-adtz']
        self.mc['CCSDT-CP-adtzadz'] = ['CCSDT-CP-adtzadz', 'CCSDT-CP-adtzhadz', 'CCSDT-CP-adtzadz']
        self.mc['CCSDT-CP-atzadz'] = ['CCSDT-CP-atzadz', 'CCSDT-CP-atzhadz', 'CCSDT-CP-atzadz']
        self.mc['CCSDT-CP-atqzadz'] = ['CCSDT-CP-atqzadz', 'CCSDT-CP-atqzhadz', 'CCSDT-CP-atqzadz']
        self.mc['CCSDT-CP-atzadtz'] = ['CCSDT-CP-atzadtz', 'CCSDT-CP-atzhadtz', 'CCSDT-CP-atzadtz']
        self.mc['CCSDT-CP-atqzadtz'] = ['CCSDT-CP-atqzadtz', 'CCSDT-CP-atqzhadtz', 'CCSDT-CP-atqzadtz']
        self.mc['CCSDT-CP-atqzatz'] = ['CCSDT-CP-atqzatz', 'CCSDT-CP-atqzhatz', 'CCSDT-CP-atqzatz']

# print certain statistic for all 4 db and summary and indiv sys if min or max

fnreservoir = {}
fnreservoir['blankslat'] = r"""Errors with respect to Benchmark. Guide lines are at 0, 0.3, and 1.0 kcal/mol overbound ($-$) and underbound ($+$)."""
fnreservoir['5min'] = r"""Only equilibrium and near-equilibrium systems included. (All S22 and HSG, 50/194 NBC10, 28/118 HBC6.)"""
fnreservoir['liliowa'] = r"""{0}MAE (dark by {1} kcal/mol) for subsets in residue classes cation, anion, polar, aliphatic, \& aromatic (L to R)."""
fnreservoir['flat'] = r"""{0}Errors with respect to benchmark within $\pm${1} kcal/mol. Guide lines are at {2} overbound ($-$) and underbound ($+$)."""

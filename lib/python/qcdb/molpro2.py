#
#@BEGIN LICENSE
#
# PSI4: an ab initio quantum chemistry software package
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
#@END LICENSE
#

import math

from exceptions import *
import qcformat
import molpro_basissets
import options


class Infile(qcformat.InputFormat2):

    def __init__(self, mem, mol, mtd, der, opt):
        qcformat.InputFormat2.__init__(self, mem, mol, mtd, der, opt)

        #print self.method, self.molecule.nactive_fragments()
        if ('sapt' in self.method or 'mp2c' in self.method) and self.molecule.nactive_fragments() != 2:
            raise FragmentCountError("""Requested molecule has %d, not 2, fragments.""" % (self.molecule.nactive_fragments()))

#        # memory in MB --> MW
#        self.memory = int(math.ceil(mem / 8.0))
        # auxiliary basis sets
        [self.unaugbasis, self.augbasis, self.auxbasis] = self.corresponding_aux_basis()

    def muster_basis_options(self):
        text = ''
        lowername = self.method.lower()
        options = defaultdict(lambda: defaultdict(dict))
    
        options['BASIS']['ORBITAL']['value'] = self.basis
    
        if self.method in ['ccsd(t)-f12-optri']:
            if self.basis == 'cc-pvdz-f12':
                options['BASIS']['JKFIT']['value'] = 'aug-cc-pvtz/jkfit'
                options['BASIS']['JKFITC']['value'] = self.basis + '/optri'
                options['BASIS']['MP2FIT']['value'] = 'aug-cc-pvtz/mp2fit'
        elif ('df-' in self.method) or ('f12' in self.method) or (self.method in ['mp2c', 'dft-sapt', 'dft-sapt-pbe0acalda']):
            if self.unaugbasis and self.auxbasis:
                options['BASIS']['JKFIT']['value'] = self.auxbasis + '/jkfit'
                options['BASIS']['JKFITB']['value'] = self.unaugbasis + '/jkfit'
                options['BASIS']['MP2FIT']['value'] = self.auxbasis + '/mp2fit'
                options['BASIS']['DFLHF']['value'] = self.auxbasis + '/jkfit'
            else:
                raise ValidationError("""Auxiliary basis not predictable from orbital basis '%s'""" % (self.basis))
        return text, options

    def prepare_basis_for_molpro(self):
        text = ''
    
        for opt, val in self.options['BASIS'].items():
                #print opt, val['value']
                #print molpro_basissets.altbasis.keys()
                if not text:
                    text += """basis={\n"""
                try:
                    # jaxz, maxz, etc.
                    for line in molpro_basissets.altbasis[val['value']]:
                        text += """%s\n""" % (line)
                    text += '\n'
                except KeyError:
                    # haxz
                    if val['value'].startswith('heavy-aug-'):
                        text += """set,%s; default,%s,H=%s\n""" % (opt.lower(), self.augbasis, self.unaugbasis)
                    # xz, axz, 6-31g*
                    else:
                        text += """set,%s; default,%s\n""" % (opt.lower(), val['value'])
    
        if text:
            text += """}\n\n"""

        return text

    def format_infile_string(self):
        """

        """
        # Handle memory and comment
        memcmd, _memkw = """***, %s\nmemory,%d,m\n""" % (self.molecule.tagline, int(math.ceil(self.memory / 8.0))), {}

        # Handle molecule and basis set
        molcmd, _molkw = self.molecule.format_molecule_for_molpro(), {}



        # format global convergence directions
#        text += self.format_global_parameters()
        _cdscmd, cdskw = muster_cdsgroup_options(self.method)

        # Handle calc type and quantum chemical method
        mdccmd, mdckw, mdcls = procedures['energy'][self.method](self.method, self.dertype, self.molecule)
        _bascmd, baskw = self.muster_basis_options()


#        # format options
#        optcmd = qcdb.options.prepare_options_for_psi4(mdckw)

# make options from imdb only user options (currently non-existent). set basis and castup from here.
        # Handle driver vs input/default keyword reconciliation
        userkw = self.options
#        userkw = p4util.prepare_options_for_modules()
        #userkw = qcdb.options.reconcile_options(userkw, memkw)
        #userkw = qcdb.options.reconcile_options(userkw, molkw)
        userkw = options.reconcile_options2(userkw, cdskw)
        userkw = options.reconcile_options2(userkw, baskw)
        #userkw = qcdb.options.reconcile_options(userkw, psikw)
        userkw = options.reconcile_options2(userkw, mdckw)

        # Handle conversion of psi4 keyword structure into cfour format
        #optcmdB = options.prepare_options_for_psi4(userkw)
        optcmd = prepare_options_for_molpro(userkw, mdcls)
        bascmd, _baskw = self.prepare_basis_for_molpro(), {} #self.options['BASIS']), {}

        # Handle text to be passed untouched
        litcmd = """\nshow[1,20f20.12],ee*,ce*,te*\nshow[1,60f20.12],_E*\n\n"""


        # Assemble infile pieces
        return memcmd + molcmd + bascmd + optcmd + mdccmd + litcmd


def muster_cdsgroup_options(name):
    text = ''
    lowername = name.lower()
    options = defaultdict(lambda: defaultdict(dict))

    options['GTHRESH']['ZERO']['value'] = 1.0e-14
    options['GTHRESH']['ONEINT']['value'] = 1.0e-14
    options['GTHRESH']['TWOINT']['value'] = 1.0e-14
    options['GTHRESH']['ENERGY']['value'] = 1.0e-9

    if name in ['mp2c', 'dft-sapt-shift', 'dft-sapt', 'dft-sapt-pbe0ac', 'dft-sapt-pbe0acalda']:
        options['GTHRESH']['ENERGY']['value'] = 1.0e-8
        options['GTHRESH']['ORBITAL']['value'] = 1.0e-8
        options['GTHRESH']['GRID']['value'] = 1.0e-8
    elif name in ['b3lyp', 'b3lyp-d', 'df-b3lyp', 'df-b3lyp-d']:
        options['GTHRESH']['ENERGY']['value'] = 1.0e-8
        options['GTHRESH']['ORBITAL']['value'] = 1.0e-7
        options['GTHRESH']['GRID']['value'] = 1.0e-8
    else:
        pass

    return text, options





def prepare_options_for_molpro(options, proc):
    """Function to take the full snapshot of the liboptions object
    encoded in dictionary *options*, find the options directable toward
    Cfour (options['CFOUR']['CFOUR_**']) that aren't default, then write
    a CFOUR deck with those options.
    Note that unlike the cfour version, this uses complete options deck.

    """
    text = ''

    if len(options['GTHRESH']) > 0:
        text += 'gthresh'
        for opt, val in options['GTHRESH'].items():
            text += """,%s=%s""" % (opt, val['value'])
        text += '\n\n'

    for item in proc:
        if len(options[item.upper()]) > 0:
            text += """{%s%s}\n""" % (item, options[item.upper()]['OPTIONS']['value'])
        else:
            text += """%s\n""" % (item)

    if text:
        text += '\n'

    return text


def muster_modelchem(name, dertype, mol):
    """Transform calculation method *name* and derivative level *dertype*
    into options for cfour. While deliberately requested pieces,
    generally |cfour__cfour_deriv_level| and |cfour__cfour_calc_level|,
    are set to complain if contradicted ('clobber' set to True), other
    'recommended' settings, like |cfour__cfour_cc_program|, can be
    countermanded by keywords in input file ('clobber' set to False).
    Occasionally, want these pieces to actually overcome keywords in
    input file ('superclobber' set to True).

    """
    text = ''
    lowername = name.lower()
    options = defaultdict(lambda: defaultdict(dict))
    proc = []

    if dertype == 0:
        pass
    else:
        raise ValidationError("""Requested Psi4 dertype %d is not available.""" % (dertype))

    if lowername == 'mp2':
        pass
        options['GLOBALS']['FREEZE_CORE']['value'] = True
        options['SCF']['SCF_TYPE']['value'] = 'direct'
        options['MP2']['MP2_TYPE']['value'] = 'conv'
        text += """mp2')\n\n"""

    elif lowername == 'ccsd(t)-f12':
        proc.append('rhf')
        proc.append('ccsd(t)-f12')
        options['CCSD(T)-F12']['OPTIONS']['value'] = ',df_basis=mp2fit,df_basis_exch=jkfitb,ri_basis=jkfitb'

    elif lowername == 'ccsd(t)-f12-optri':
        proc.append('rhf')
        proc.append('ccsd(t)-f12')
        options['CCSD(T)-F12']['OPTIONS']['value'] = ',df_basis=mp2fit,df_basis_exch=jkfit,ri_basis=jkfitc'

    elif lowername == 'mp2c':
        proc.append('gdirect')
        proc.append(mol.extract_fragments(1, 2).format_molecule_for_molpro())
        proc.append('df-hf,')
        proc.append('df-ks,')
        proc.append('sapt; monomerA')
        options['DF-HF,']['OPTIONS']['value'] = """basis=jkfit,locorb=0; start,atdens; save,1101.2"""
        options['DF-KS,']['OPTIONS']['value'] = """lhf,df_basis=dflhf,basis_coul=jkfitb,basis_exch=jkfitb; dftfac,1.0; start,1101.2; save,2101.2"""

        proc.append(mol.extract_fragments(2, 1).format_molecule_for_molpro())
        proc.append('df-hf')
        proc.append('df-ks')
        proc.append('sapt; monomerB')
        options['DF-HF']['OPTIONS']['value'] = """,basis=jkfit,locorb=0; start,atdens; save,1102.2"""
        options['DF-KS']['OPTIONS']['value'] = """,lhf,df_basis=dflhf,basis_coul=jkfitb,basis_exch=jkfitb; dftfac,1.0; start,1102.2; save,2102.2"""

        proc.append(mol.format_molecule_for_molpro())
        proc.append('sapt; intermol')
        options['SAPT; INTERMOL']['OPTIONS']['value'] = """,saptlevel=3,ca=2101.2,cb=2102.2,icpks=0,fitlevel=3,nlexfac=0.0,cfac=0.0; dfit,basis_coul=jkfit,basis_exch=jkfit,cfit_scf=3"""

    else:
        raise ValidationError("""Requested Cfour computational methods %d is not available.""" % (lowername))

#    # Set clobbering
#    if 'CFOUR_DERIV_LEVEL' in options['CFOUR']:
#        options['CFOUR']['CFOUR_DERIV_LEVEL']['clobber'] = True
#        options['CFOUR']['CFOUR_DERIV_LEVEL']['superclobber'] = True
#    if 'CFOUR_CALC_LEVEL' in options['CFOUR']:
#        options['CFOUR']['CFOUR_CALC_LEVEL']['clobber'] = True
#        options['CFOUR']['CFOUR_CALC_LEVEL']['superclobber'] = True
#    if 'CFOUR_CC_PROGRAM' in options['CFOUR']:
#        options['CFOUR']['CFOUR_CC_PROGRAM']['clobber'] = False

    return text, options, proc

procedures = {
    'energy': {
        'mp2c'           : muster_modelchem,
        'ccsd(t)-f12'    : muster_modelchem,
        'ccsd(t)-f12-optri' : muster_modelchem,
        #'sapt0'         : muster_modelchem,
        #'sapt2+'        : muster_modelchem,
        #'sapt2+(3)'     : muster_modelchem,
        #'sapt2+3(ccd)'  : muster_modelchem,
    }
}

qcmtdIN = procedures['energy']


def psi4_list():
    """Return an array of Psi4 methods with energies.

    """
    return procedures['energy'].keys()

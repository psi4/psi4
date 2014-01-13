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

"""Module with functions for Psi4/Cfour interface. Portions that require
calls to Boost Python psi4 module are here, otherwise in qcdb module.
Also calls to qcdb module are here and not elsewhere in driver.
Organizationally, this module isolates qcdb code from psi4 code.

"""
from __future__ import print_function
import shutil
import os
import subprocess
import re
import inspect
import glob
import psi4
import p4const
import p4util
import qcdb
from p4regex import *
#from extend_Molecule import *
from molutil import *
from functional import *
from driver import *
# never import driver, wrappers, or aliases into this file


def run_cfour_module(xmod):
    # Find environment by merging PSIPATH and PATH environment variables
    lenv = os.environ
    lenv['PATH'] = ':'.join([os.path.abspath(x) for x in os.environ.get('PSIPATH', '').split(':')]) + ':' + lenv.get('PATH') + ':' + psi4.Process.environment["PSIDATADIR"] + '/basis' + ':' + psi4.psi_top_srcdir() + '/lib/basis' 

    # Call executable xcfour, directing cfour output to the psi4 output file
    try:
        retcode = subprocess.Popen([xmod], bufsize=0, stdout=subprocess.PIPE, env=lenv)
    except OSError as e:
        sys.stderr.write('Program %s not found in path or execution failed: %s\n' % (cfour_executable, e.strerror))
        #p4out.write('Program %s not found in path or execution failed: %s\n' % (cfour_executable, e.strerror))
        sys.exit(1)

    c4out = ''
    while True:
        data = retcode.stdout.readline()
        if not data:
            break
        #if psi4.outfile_name() == 'stdout':
        #    sys.stdout.write(data)
        #else:
        #    p4out.write(data)
        #    p4out.flush()
        c4out += data
    #internal_p4c4_info['output'] = c4out
    return c4out


def vpt2(name, **kwargs):
    """Perform vibrational second-order perturbation computation through
    Cfour to get anharmonic frequencies.

    .. caution:: Some features are not yet implemented. Buy a developer a coffee.

       - Presently uses all gradients. Could mix in analytic 2nd-derivs.

       - Collect resutls.

       - Manage scratch / subdir better.

       - Fix FD_PROJECT t/f issue so can be specified from input

       - Untangle CCSD(T) vs CCSD[T] and FJOBARC issue

       - Allow CFOUR_BASIS

       - Consider forcing some tighter convcrit, c4 and p4

       - Expand CURRENT DIPOLE XYZ beyond SCF

       - Remember additional FJOBARC record TOTENER2 if EXCITE .ne. NONE

    """
    lowername = name.lower()
    dertype = 1

    psi4.set_local_option('CFOUR', 'CFOUR_VIBRATION', 'FINDIF')
    psi4.set_local_option('CFOUR', 'CFOUR_FREQ_ALGORITHM', 'PARALLEL')
    psi4.set_local_option('CFOUR', 'CFOUR_ANH_ALGORITHM', 'PARALLEL')
    psi4.set_local_option('CFOUR', 'CFOUR_ANHARMONIC', 'VPT2')
    psi4.set_local_option('CFOUR', 'CFOUR_FD_PROJECT', True)  # equals OFF  # correct
    #psi4.set_local_option('CFOUR', 'CFOUR_FD_PROJECT', False)  # equals ON

    #*CFOUR(CALC=CCSD[T],BASIS=ANO0,CC_CONV=9,LINEQ_CONV=9
    #SCF_CONV=9,DROPMO=1,COORD=CARTESIAN,UNITS=BOHR
    #ANHARM=VPT2,FD_IRREPS=1
    #ANH_ALGORITHM=PARALLEL,FREQ_ALGORITHM=PARALLEL,FD_PROJECT=OFF
    #VIB=FINDIF)

    # Save submission directory
    current_directory = os.getcwd()

    # Move into job scratch directory
    psioh = psi4.IOManager.shared_object()
    psio = psi4.IO.shared_object()
    os.chdir(psioh.get_default_path())

    # Construct and move into cfour subdirectory of job scratch directory
    cfour_tmpdir = kwargs['path'] if 'path' in kwargs else \
        'psi.' + str(os.getpid()) + '.' + psio.get_default_namespace() + \
        '.cfour.' + str(random.randint(0, 99999))
    if not os.path.exists(cfour_tmpdir):
        os.mkdir(cfour_tmpdir)
    os.chdir(cfour_tmpdir)

    # Generate the ZMAT input file in scratch
    with open('ZMAT', 'w') as cfour_infile:
        cfour_infile.write(write_zmat(lowername, dertype))
    print('\n====== Begin ZMAT input for CFOUR ======')
    print(open('ZMAT', 'r').read())
    print('======= End ZMAT input for CFOUR =======\n')

    # <<<  SC2: generate the displacements that will form the harmonic freq  >>>
    run_cfour_module('xjoda')
    run_cfour_module('xsymcor')

    fdfiles = {}
    harmzmats = glob.glob('zmat*')
    for item in harmzmats:
        try:
            with open(psioh.get_default_path() + cfour_tmpdir + '/' + item, 'r') as handle:
                fdfiles[item] = handle.read()
                psi4.print_out('  CFOUR scratch file %s has been read\n' % (item))
                psi4.print_out('%s\n' % fdfiles[item])
         #       internal_p4c4_info[item.lower()] = fdfiles[item]
        except IOError:
            pass

    # <<<  SC3: run the gradients that will form the harmonic freq  >>>
    for item in harmzmats:
        harmnum = 'scr.000.' + item[-3:]
        if not os.path.exists(harmnum):
            os.mkdir(harmnum)
            os.chdir(harmnum)
        shutil.copy2('../' + item, 'ZMAT')
        shutil.copy2('../GENBAS', 'GENBAS')
        print(os.getcwd())
        
        outfile = run_cfour_module('xcfour')
        run_cfour_module('xja2fja')
        with open('../out.' + item[-3:], 'w') as cfour_outfile:
            cfour_outfile.write(outfile)
        shutil.copy2('FJOBARC', '../fja.' + item[-3:])

        os.chdir('..')

    # <<<  SC4: collect the harmonic gradients into harm.out  >>>
    harmout = run_cfour_module('xjoda')
    harmout += run_cfour_module('xsymcor')
    for part in sorted(glob.glob('fja.*')):
        shutil.copy2('fja.' + part[-3:], 'FJOBARC')
        harmout += run_cfour_module('xja2fja')
        harmout += run_cfour_module('xsymcor')
    harmout += run_cfour_module('xjoda')
    for item in harmzmats:
        os.remove(os.getcwd() + '/' + item)
    harmout += run_cfour_module('xcubic')

    with open('harm.out', 'w') as cfour_outfile:
        cfour_outfile.write(harmout)

    # <<<  SC5 & SC6 setup: generate displacements along normal modes and prep dirs  >>>
    anharmzmats = glob.glob('zmat*')
    for item in anharmzmats:
        anhidx = item[-3:]
        anhdir = 'scr.' + anhidx
        if not os.path.exists(anhdir):
            os.mkdir(anhdir)
            os.chdir(anhdir)
        shutil.copy2('../' + item, 'ZMAT')
        shutil.copy2('../GENBAS', 'GENBAS')
        #! ln -s $ecpdir/ECPDATA $j/ECPDATA
        run_cfour_module('xjoda')
        run_cfour_module('xsymcor')

        anharm2zmats = glob.glob('zmat*')
        for item2 in anharm2zmats:
            anh2idx = item2[-3:]
            anh2dir = 'scr.' + anh2idx
            if not os.path.exists(anh2dir):
                os.mkdir(anh2dir)
                os.chdir(anh2dir)
            shutil.copy2('../' + item2, 'ZMAT')
            shutil.copy2('../GENBAS', 'GENBAS')
            #! ln -s $ecpdir/ECPDATA $j/ECPDATA
            os.chdir('..')

        os.chdir('..')

    # <<<  SC6 run: run disp of disp  >>>
    for item in glob.glob('zmat*'):
        os.chdir('scr.' + item[-3:])

        for item2 in glob.glob('zmat*'):
            os.chdir('scr.' + item2[-3:])
            print(os.getcwd())

            outfile = run_cfour_module('xcfour')
            run_cfour_module('xja2fja')
            with open('../out.' + item2[-3:], 'w') as cfour_outfile:
                cfour_outfile.write(outfile)
            shutil.copy2('FJOBARC', '../fja.' + item2[-3:])

            os.chdir('..')
        os.chdir('..')

    # <<<  SC7: assemble gradients into anharm.out  >>>
    procdir = 'proc'
    if os.path.exists(procdir):
        shutil.rmtree(procdir)
    os.mkdir(procdir)
    #os.chdir(procdir)
    print("Processing initial harmonic frequency calculation")
    run_cfour_module('xclean')
    anharmout = run_cfour_module('xjoda')
    anharmout += run_cfour_module('xsymcor')

    for part in sorted(glob.glob('fja.*')):
        shutil.copy2('fja.' + part[-3:], 'FJOBARC')
        anharmout += run_cfour_module('xja2fja')
        anharmout += run_cfour_module('xsymcor')
    anharmout += run_cfour_module('xjoda')
    anharmout += run_cfour_module('xcubic')
    shutil.copy2('JOBARC', 'jobarc.0')
    shutil.copy2('JAINDX', 'jaindx.0')

    for item in sorted(glob.glob('scr.[0-9][0-9][0-9]')):
        print("Processing output from displaced point %s" % item)
        anharmout += """
----------------------------------------------------------------------
Processing output from displaced point %s
----------------------------------------------------------------------
""" % (item)
        os.chdir(item)

        run_cfour_module('xclean')
        anharmout += run_cfour_module('xjoda')
        anharmout += run_cfour_module('xsymcor')

        for part in sorted(glob.glob('fja.*')):
            shutil.copy2('fja.' + part[-3:], 'FJOBARC')
            anharmout += run_cfour_module('xja2fja')
            anharmout += run_cfour_module('xsymcor')
        anharmout += run_cfour_module('xjoda')
        os.remove(os.getcwd() + '/FJOBARC')
        anharmout += run_cfour_module('xja2fja')
        shutil.copy2('FJOBARC', '../' + procdir + '/fja.' + item[-3:])
        os.chdir('..')

    print("Processing final output")
    anharmout += """
----------------------------------------------------------------------
Processing final output
----------------------------------------------------------------------
"""

    os.chdir(procdir)
    shutil.copy2('../jobarc.0', 'JOBARC')
    shutil.copy2('../jaindx.0', 'JAINDX')

    for part in sorted(glob.glob('fja.*')):
        shutil.copy2('fja.' + part[-3:], 'FJOBARC')
        anharmout += run_cfour_module('xja2fja')
        anharmout += run_cfour_module('xcubic')
    os.chdir('..')

    with open('anharm.out', 'w') as cfour_outfile:
        cfour_outfile.write(anharmout)
    psi4.print_out(anharmout)



#   mkdir proc
#   rm proc\* 
#   echo "Processing initial harmonic frequency calculation"
#   xclean
#   xjoda > anharm.out
#   xsymcor >> anharm.out
#
#   foreach i (fja*)
#       cp $i FJOBARC
#       xja2fja >> anharm.out
#       xsymcor >> anharm.out
#
#   xjoda >> anharm.out
#   xcubic >> anharm.out
#   cp JOBARC jobarc.0
#   cp JAINDX jaindx.0
#   foreach i (0*):
#       echo "----------------------------------------------------------------------" >> anharm.out
#       echo "Processing output from displaceed point $i" >> anharm.out
#       echo "Processing output from displaceed point $i" 
#       echo "----------------------------------------------------------------------" >> anharm.out
#       cd $i
#       xclean
#       xjoda >> ../anharm.out
#       xsymcor >> ../anharm.out
#
#       foreach j (fja*):
#          cp $j FJOBARC
#          xja2fja >> ../anharm.out
#          xsymcor >> ../anharm.out
#
#       xjoda >> ../anharm.out
#       rm FJOBARC
#       xja2fja >> ../anharm.out
#       cp FJOBARC ../proc/fja_all.$i
#       cd ..
#
#   echo "----------------------------------------------------------------------" >> anharm.out
#   echo "Processing final output" >> anharm.out
#   echo "Processing final output" 
#   echo "----------------------------------------------------------------------" >> anharm.out
#   cd proc
#   cp ../jobarc.0 JOBARC
#   cp ../jaindx.0 JAINDX
#
#   foreach i (fja*):
#       cp $i FJOBARC
#       xja2fja >> ../anharm.out
#       xcubic >> ../anharm.out


def p4vpt2(name, **kwargs):
    """Perform vibrational second-order perturbation computation through
    Cfour to get anharmonic frequencies. This version uses c4 for the disp
    and pt2 but gets gradients from p4.

    .. caution:: Some features are not yet implemented. Buy a developer a coffee.

       - Presently uses all gradients. Could mix in analytic 2nd-derivs.

       - Collect resutls.

       - Manage scratch / subdir better.

       - Fix FD_PROJECT t/f issue so can be specified from input

       - Untangle CCSD(T) vs CCSD[T] and FJOBARC issue

       - Allow CFOUR_BASIS

       - Consider forcing some tighter convcrit, c4 and p4

    """
    optstash = p4util.OptionsState(
        ['BASIS'])
    user_basis = psi4.get_global_option('BASIS')
    lowername = name.lower()
    dertype = 1

    psi4.set_local_option('CFOUR', 'CFOUR_VIBRATION', 'FINDIF')
    psi4.set_local_option('CFOUR', 'CFOUR_FREQ_ALGORITHM', 'PARALLEL')
    psi4.set_local_option('CFOUR', 'CFOUR_ANH_ALGORITHM', 'PARALLEL')
    psi4.set_local_option('CFOUR', 'CFOUR_ANHARMONIC', 'VPT2')
    psi4.set_local_option('CFOUR', 'CFOUR_FD_PROJECT', True)  # equals OFF  # correct
    #psi4.set_local_option('CFOUR', 'CFOUR_FD_PROJECT', False)  # equals ON

    # A cheap (compchem-wise) Cfour skeleton is run underneath the
    #   Psi4 requested modelchem is order to have something to hang
    #   the P4 gradients upon. These are options for that calc.
    psi4.set_global_option('BASIS', '3-21G')
    psi4.set_local_option('CFOUR', 'CFOUR_CALC_LEVEL', 'SCF')
    psi4.set_local_option('CFOUR', 'CFOUR_SCF_CONV', 5)

    # Save submission directory
    current_directory = os.getcwd()

    # Move into job scratch directory
    psioh = psi4.IOManager.shared_object()
    psio = psi4.IO.shared_object()
    os.chdir(psioh.get_default_path())

    psioh.set_specific_retention(32, True)

    # Construct and move into cfour subdirectory of job scratch directory
    cfour_tmpdir = kwargs['path'] if 'path' in kwargs else \
        'psi.' + str(os.getpid()) + '.' + psio.get_default_namespace() + \
        '.cfour.' + str(random.randint(0, 99999))
    if not os.path.exists(cfour_tmpdir):
        os.mkdir(cfour_tmpdir)
    os.chdir(cfour_tmpdir)

    # Generate the ZMAT input file in scratch
    with open('ZMAT', 'w') as cfour_infile:
        cfour_infile.write(write_zmat('c4-scf', dertype))  # duplicate of CFOUR_CALC_LEVEL setting TODO
    print('\n====== Begin ZMAT input for CFOUR ======')
    print(open('ZMAT', 'r').read())
    print('======= End ZMAT input for CFOUR =======\n')

    # <<<  SC2: generate the displacements that will form the harmonic freq  >>>
    run_cfour_module('xjoda')
    run_cfour_module('xsymcor')

    fdfiles = {}
    harmzmats = sorted(glob.glob('zmat*'))
    for item in harmzmats:
        try:
            with open(psioh.get_default_path() + cfour_tmpdir + '/' + item, 'r') as handle:
                fdfiles[item] = handle.read()
                psi4.print_out('  CFOUR scratch file %s has been read\n' % (item))
                psi4.print_out('%s\n' % fdfiles[item])
         #       internal_p4c4_info[item.lower()] = fdfiles[item]
        except IOError:
            pass
    print(fdfiles['zmat005'])

    grdfiles = {}
    c4mols = {}  # at the c4 std orientation that we need
    print('<<<  SC3: (C4 skel) run the gradients that will form the harmonic freq  >>>')
    # <<<  SC3: (C4 skel) run the gradients that will form the harmonic freq  >>>
    for item in harmzmats:
        harmnum = 'scr.000.' + item[-3:]
        if not os.path.exists(harmnum):
            os.mkdir(harmnum)
            os.chdir(harmnum)
        shutil.copy2('../' + item, 'ZMAT')
        shutil.copy2('../GENBAS', 'GENBAS')
        print(os.getcwd())
        
        outfile = run_cfour_module('xcfour')
        run_cfour_module('xja2fja')
        with open('../outC.' + item[-3:], 'w') as cfour_outfile:
            cfour_outfile.write(outfile)
        shutil.copy2('FJOBARC', '../fjaC.' + item[-3:])
        try:
            with open(psioh.get_default_path() + cfour_tmpdir + '/' + harmnum + '/GRD', 'r') as handle:
                grdfiles[item] = handle.read()
                psi4.print_out('  CFOUR scratch file %s has been read\n' % (item + 'GRD'))
                psi4.print_out('%s\n' % grdfiles[item])
         #       internal_p4c4_info[item.lower()] = fdfiles[item]
        except IOError:
            pass
        print(grdfiles[item])
        c4mols[item], temp = qcdb.cfour.harvest_GRD(grdfiles[item])
        c4mols[item].print_out()


        os.chdir('..')

    # mirror file handle?
    os.chdir(current_directory)
    print('presc3', os.getcwd())
    print('<<<  SC3: (P4 reqd) run the gradients that will form the harmonic freq  >>>')
    # <<<  SC3: (P4 reqd) run the gradients that will form the harmonic freq  >>>
    psi4.set_global_option('BASIS', user_basis)
    for item in harmzmats:
        harmnum = 'scr.000.' + item[-3:]
        #print(fdfiles[item])

        psi4.IO.set_default_namespace('000-' + item[-3:])
        p4util.banner('VPT2: Displacement %s' % ('000-' + item[-3:]))
        qcdbmolecule = qcdb.cfour.harvest_zmat(fdfiles[item])
        psimolecule = geometry(qcdbmolecule.create_psi4_string_from_molecule(), 'disp-000-' + item[-3:])
        #print(harmnum, qcdbmolecule.create_psi4_string_from_molecule())
    #    psimolecule.reset_point_group('c1')
#        psimolecule.fix_orientation(False)
#        psimolecule.fix_com(False)
        psimolecule.update_geometry()
        psimolecule.print_out()
        pcrd = psimolecule.geometry()
        
        #shutil.copy2('../' + item, 'ZMAT')
        #shutil.copy2('../GENBAS', 'GENBAS')
        #print(os.getcwd())
        
        gradient(lowername, **kwargs)
        pgrd = psi4.get_gradient()
        pgrd.print_out()
        fjgrd = []
        fjcoord = []
        fjelem = []
        fjdip = []  # TODO this is going to need rotation to c4Native frame
        for at in range(psimolecule.natom()):
            fjgrd.append([pgrd.get(at, 0), pgrd.get(at, 1), pgrd.get(at, 2)])
            fjcoord.append([pcrd.get(at, 0), pcrd.get(at, 1), pcrd.get(at, 2)])
            fjelem.append(psimolecule.Z(at))
        print(harmnum)
        fje = psi4.get_variable('CURRENT ENERGY')
        print(fje)
        print(fjcoord)
        print(fjgrd)
        fjdip.append(psi4.get_variable('CURRENT DIPOLE X') / p4const.psi_dipmom_au2debye)
        fjdip.append(psi4.get_variable('CURRENT DIPOLE Y') / p4const.psi_dipmom_au2debye)
        fjdip.append(psi4.get_variable('CURRENT DIPOLE Z') / p4const.psi_dipmom_au2debye)

        #c4NativeElem, c4NativeCoord, c4NativeGrad, ElemMap = qcdb.cfour.backtransform_grad(qcdbmolecule, c4mols[item], fjgrd)
        c4NativeElem, c4NativeCoord, c4NativeGrad, ElemMap, c4NativeDip = qcdb.cfour.backtransform_grad(qcdbmolecule, c4mols[item], fjgrd, fjdip)

        fjobarc = qcdb.cfour.format_fjobarc(fje, c4NativeElem, c4NativeCoord, c4NativeGrad, ElemMap, c4NativeDip)
        #fjobarc = qcdb.cfour.format_fjobarc(fje, c4NativeElem, c4NativeCoord, c4NativeGrad, ElemMap, fjdip)
        #fjobarc = qcdb.cfour.format_fjobarc(fje, fjelem, fjcoord, fjgrd)
        print('fjobarc', item, '\n', fjobarc)
        psi4.clean()
        psi4.clean_variables()

        with open(psioh.get_default_path() + cfour_tmpdir + '/' + 'fjaP.' + item[-3:], 'w') as cfour_fjobarc:
            cfour_fjobarc.write(fjobarc)
        print('moving psi jobarc to ', psioh.get_default_path() + cfour_tmpdir + '/' + 'fjaP.' + item[-3:])
        #os.chdir('..')

    os.chdir(psioh.get_default_path())
    os.chdir(cfour_tmpdir)

    print('<<<  SC4: collect the harmonic gradients into harm.out  >>>')
    # <<<  SC4: collect the harmonic gradients into harm.out  >>>
    print(os.getcwd())
    harmout = run_cfour_module('xjoda')
    harmout += run_cfour_module('xsymcor')
    for part in sorted(glob.glob('fjaP.*')):
        shutil.copy2('fjaP.' + part[-3:], 'FJOBARC')
        harmout += run_cfour_module('xja2fja')
        harmout += run_cfour_module('xsymcor')
    harmout += run_cfour_module('xjoda')
    for item in harmzmats:
        os.remove(os.getcwd() + '/' + item)
    harmout += run_cfour_module('xcubic')

    with open('harm.out', 'w') as cfour_outfile:
        cfour_outfile.write(harmout)

    fdfiles2 = {}
    print('<<<  SC5 & SC6 setup: generate displacements along normal modes and prep dirs  >>>')
    # <<<  SC5 & SC6 setup: generate displacements along normal modes and prep dirs  >>>
    anharmzmats = sorted(glob.glob('zmat*'))
    for item in anharmzmats:
        fdfiles2[item] = {}
        anhidx = item[-3:]
        anhdir = 'scr.' + anhidx
        if not os.path.exists(anhdir):
            os.mkdir(anhdir)
            os.chdir(anhdir)
        shutil.copy2('../' + item, 'ZMAT')
        shutil.copy2('../GENBAS', 'GENBAS')
        #! ln -s $ecpdir/ECPDATA $j/ECPDATA
        run_cfour_module('xjoda')
        run_cfour_module('xsymcor')

        anharm2zmats = sorted(glob.glob('zmat*'))
        for item2 in anharm2zmats:
            anh2idx = item2[-3:]
            anh2dir = 'scr.' + anh2idx
            if not os.path.exists(anh2dir):
                os.mkdir(anh2dir)
                os.chdir(anh2dir)
            shutil.copy2('../' + item2, 'ZMAT')
            shutil.copy2('../GENBAS', 'GENBAS')
            #! ln -s $ecpdir/ECPDATA $j/ECPDATA

            print('getting zmat from ', psioh.get_default_path() + cfour_tmpdir + '/' + anhdir + '/' + item2)
            try:
                with open(psioh.get_default_path() + cfour_tmpdir + '/' + anhdir + '/' + item2, 'r') as handle:
                    fdfiles2[item][item2] = handle.read()
                    psi4.print_out('  CFOUR scratch file %s has been read\n' % (item + '-' + item2))
                    #psi4.print_out('%s\n' % fdfiles2[item][item2])
             #       internal_p4c4_info[item.lower()] = fdfiles2[item]
            except IOError:
                pass

            os.chdir('..')
        os.chdir('..')
    print('sample zmat\n', fdfiles2['zmat005']['zmat006'])

    grdfiles2 = {}
    c4mols2 = {}
    print('<<<  SC6 run: (c4 skel) run disp of disp  >>>')
    # <<<  SC6 run: run disp of disp  >>>
    for item in glob.glob('zmat*'):
        os.chdir('scr.' + item[-3:])
        grdfiles2[item] = {}
        c4mols2[item] = {}

        for item2 in glob.glob('zmat*'):
            os.chdir('scr.' + item2[-3:])
            print(os.getcwd())

            outfile = run_cfour_module('xcfour')
            run_cfour_module('xja2fja')
            with open('../out.' + item2[-3:], 'w') as cfour_outfile:
                cfour_outfile.write(outfile)
            shutil.copy2('FJOBARC', '../fjaC.' + item2[-3:])

            try:
                with open(psioh.get_default_path() + cfour_tmpdir + '/scr.' + item[-3:] + '/scr.' + item2[-3:] + '/GRD', 'r') as handle:
                    grdfiles2[item][item2] = handle.read()
                    psi4.print_out('  CFOUR scratch file %s has been read\n' % (item + '-' + item2 + 'GRD'))
                    psi4.print_out('%s\n' % grdfiles2[item][item2])
            except IOError:
                pass
            print(grdfiles2[item][item2])
            c4mols2[item][item2], temp = qcdb.cfour.harvest_GRD(grdfiles2[item][item2])
            c4mols2[item][item2].print_out()

            os.chdir('..')
        os.chdir('..')

    os.chdir(current_directory)
    print('<<<  SC6 run: (p4 reqd) run disp of disp  >>>')
    # <<<  SC6 run: run disp of disp  >>>
    for key, val in sorted(fdfiles2.items()):
        for key2 in sorted(val):
            print('\n<<<  ' + key + ' ' + key2 + '  >>>\n\n')
            #if key == 'zmat002' and key2 == 'zmat005':
            #    print(fdfiles2[key][key2])
    
            psi4.IO.set_default_namespace(key[-3:] + '-' + key2[-3:])
            p4util.banner('VPT2: Displacement %s-%s' % (key[-3:], key2[-3:]))
            qcdbmolecule = qcdb.cfour.harvest_zmat(fdfiles2[key][key2])
            psimolecule = geometry(qcdbmolecule.create_psi4_string_from_molecule(), 'disp-' + key[-3:] + '-' + key2[-3:])
            psimolecule.update_geometry()
            psimolecule.print_out()
            pcrd = psimolecule.geometry()

            gradient(lowername, **kwargs)
            pgrd = psi4.get_gradient()
            pgrd.print_out()
            fjgrd = []
            fjcoord = []
            fjelem = []
            fjdip = []
            for at in range(psimolecule.natom()):
                fjgrd.append([pgrd.get(at, 0), pgrd.get(at, 1), pgrd.get(at, 2)])
                fjcoord.append([pcrd.get(at, 0), pcrd.get(at, 1), pcrd.get(at, 2)])
                fjelem.append(psimolecule.Z(at))
            #print(harmnum)
            fje = psi4.get_variable('CURRENT ENERGY')
            fjdip.append(psi4.get_variable('CURRENT DIPOLE X') / p4const.psi_dipmom_au2debye)
            fjdip.append(psi4.get_variable('CURRENT DIPOLE Y') / p4const.psi_dipmom_au2debye)
            fjdip.append(psi4.get_variable('CURRENT DIPOLE Z') / p4const.psi_dipmom_au2debye)
            print(fje)
            print(fjcoord)
            print(fjgrd)

            #c4NativeElem, c4NativeCoord, c4NativeGrad, ElemMap = qcdb.cfour.backtransform_grad(qcdbmolecule, c4mols2[key][key2], fjgrd)
            c4NativeElem, c4NativeCoord, c4NativeGrad, ElemMap, c4NativeDip = qcdb.cfour.backtransform_grad(qcdbmolecule, c4mols2[key][key2], fjgrd, fjdip)

            #fjobarc = qcdb.cfour.format_fjobarc(fje, c4NativeElem, c4NativeCoord, c4NativeGrad, ElemMap, fjdip)
            fjobarc = qcdb.cfour.format_fjobarc(fje, c4NativeElem, c4NativeCoord, c4NativeGrad, ElemMap, c4NativeDip)
            print('fjobarc', key, key2, '\n', fjobarc)
            psi4.clean()
            psi4.clean_variables()

            with open(psioh.get_default_path() + cfour_tmpdir + '/scr.' + key[-3:] + '/fjaP.' + key2[-3:], 'w') as cfour_fjobarc:
                cfour_fjobarc.write(fjobarc)
            print('moving psi jobarc to ', psioh.get_default_path() + cfour_tmpdir + '/scr.' + key[-3:] + '/fjaP.' + key2[-3:])

    os.chdir(psioh.get_default_path())
    os.chdir(cfour_tmpdir)
    print('<<<  SC7: assemble gradients into anharm.out  >>>')
    # <<<  SC7: assemble gradients into anharm.out  >>>
    procdir = 'proc'
    if os.path.exists(procdir):
        shutil.rmtree(procdir)
    os.mkdir(procdir)
#    #os.chdir(procdir)
    print("Processing initial harmonic frequency calculation")
    run_cfour_module('xclean')
    anharmout = run_cfour_module('xjoda')
    anharmout += run_cfour_module('xsymcor')

    for part in sorted(glob.glob('fjaP.*')):
        shutil.copy2('fjaP.' + part[-3:], 'FJOBARC')
        anharmout += run_cfour_module('xja2fja')
        anharmout += run_cfour_module('xsymcor')
    anharmout += run_cfour_module('xjoda')
    anharmout += run_cfour_module('xcubic')
    shutil.copy2('JOBARC', 'jobarc.0')
    shutil.copy2('JAINDX', 'jaindx.0')

    for item in sorted(glob.glob('scr.[0-9][0-9][0-9]')):
        print("Processing output from displaced point %s" % item)
        anharmout += """
----------------------------------------------------------------------
Processing output from displaced point %s
----------------------------------------------------------------------
""" % (item)
        os.chdir(item)

        run_cfour_module('xclean')
        anharmout += run_cfour_module('xjoda')
        anharmout += run_cfour_module('xsymcor')

        for part in sorted(glob.glob('fjaP.*')):
            shutil.copy2('fjaP.' + part[-3:], 'FJOBARC')
            anharmout += run_cfour_module('xja2fja')
            anharmout += run_cfour_module('xsymcor')
        anharmout += run_cfour_module('xjoda')
        os.remove(os.getcwd() + '/FJOBARC')
        anharmout += run_cfour_module('xja2fja')
        shutil.copy2('FJOBARC', '../' + procdir + '/fja.' + item[-3:])
        os.chdir('..')

    print("Processing final output")
    anharmout += """
----------------------------------------------------------------------
Processing final output
----------------------------------------------------------------------
"""

    os.chdir(procdir)
    shutil.copy2('../jobarc.0', 'JOBARC')
    shutil.copy2('../jaindx.0', 'JAINDX')

    for part in sorted(glob.glob('fja.*')):
        shutil.copy2('fja.' + part[-3:], 'FJOBARC')
        anharmout += run_cfour_module('xja2fja')
        anharmout += run_cfour_module('xcubic')
    os.chdir('..')

    with open('anharm.out', 'w') as cfour_outfile:
        cfour_outfile.write(anharmout)
    psi4.print_out(anharmout)

    optstash.restore()

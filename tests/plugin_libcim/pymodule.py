import PsiMod
import re
import os
import input
import math
import warnings
from driver import *
from wrappers import *
from molutil import *
from text import *
from procutil import *
# for the parallel version:
#from mpi4py import MPI

def run_plugin_libcim(name, **kwargs):
    """Function encoding sequence of PSI module and plugin calls so that
    plugin_libcim can be called via :py:func:`driver.energy`.

    >>> energy('cim-ccsd(t)')

    """
    lowername = name.lower()
    kwargs = kwargs_lower(kwargs)

    # override symmetry:
    molecule = PsiMod.get_active_molecule()
    molecule.update_geometry()
    molecule.reset_point_group('c1')
    molecule.fix_orientation(1)
    molecule.update_geometry()

    # what type of cim?
    if (lowername == 'cim-ccsd'):
        PsiMod.set_global_option('compute_triples', False)
        PsiMod.set_global_option('cim_correlated_method', 'CCSD')
        method = 0
    if (lowername == 'cim-ccsd(t)'):
        PsiMod.set_global_option('compute_triples', True)
        PsiMod.set_global_option('cim_correlated_method', 'CCSD(T)')
        method = 0
    if (lowername == 'cim-cepa(0)'):
        PsiMod.set_global_option('cepa_level', 'cepa0')
        PsiMod.set_global_option('cim_correlated_method', 'CEPA')
        method = 1
    if (lowername == 'cim-cepa(1)'):
        PsiMod.set_global_option('cepa_level', 'cepa1')
        PsiMod.set_global_option('cim_correlated_method', 'CEPA')
        method = 1
    if (lowername == 'cim-cepa(2)'):
        PsiMod.set_global_option('cepa_level', 'cepa2')
        PsiMod.set_global_option('cim_correlated_method', 'CEPA')
        method = 1
    if (lowername == 'cim-cepa(3)'):
        PsiMod.set_global_option('cepa_level', 'cepa3')
        PsiMod.set_global_option('cim_correlated_method', 'CEPA')
        method = 1
    if (lowername == 'cim-cisd'):
        PsiMod.set_global_option('cepa_level', 'cisd')
        PsiMod.set_global_option('cim_correlated_method', 'CEPA')
        method = 1
    if (lowername == 'cim-acpf'):
        PsiMod.set_global_option('cepa_level', 'acpf')
        PsiMod.set_global_option('cim_correlated_method', 'CEPA')
        method = 1
    if (lowername == 'cim-aqcc'):
        PsiMod.set_global_option('cepa_level', 'aqcc')
        PsiMod.set_global_option('cim_correlated_method', 'CEPA')
        method = 1

    # the options that are not unique to the plugins are not set 
    # properly when i load more than one plugin.  set them here
    if PsiMod.has_option_changed('r_convergence') == False:
       PsiMod.set_global_option('r_convergence', 1e-7)
    if PsiMod.has_option_changed('maxiter') == False:
       PsiMod.set_global_option('maxiter', 100)
    if PsiMod.has_option_changed('occ_tolerance') == False:
       PsiMod.set_global_option('occ_tolerance', 5e-5)

    # if df basis not set, pick one. there is no easy way to guarantee 
    # the ri basis corresponding to the primary basis exists, so default
    # to the primary basis
    if (PsiMod.get_global_option('df_basis_mp2').lower()==''):
       PsiMod.set_global_option('df_basis_mp2', PsiMod.get_global_option('basis'))

    # run scf
    energy('scf', **kwargs)

    # initialize cim (determine clusters)
    PsiMod.set_global_option('cim_initialize', True)
    PsiMod.plugin('plugin_libcim.so')

    cluster_ccsd = []
    cluster_ccsdt = []
    cluster_t = []
    built_ccsd = 0.0
    built_t = 0.0
    built_ccsdt = 0.0
    built_energy = 0.0
    escf = PsiMod.get_variable('SCF TOTAL ENERGY')
    PsiMod.set_global_option('cim_initialize', False)
    cim_n = 0
    while cim_n < PsiMod.get_variable('CIM CLUSTERS'):
        # run plugin_libcim:
        PsiMod.set_global_option('CIM_CLUSTER_NUM', cim_n)
        PsiMod.plugin('plugin_libcim.so')

        # accumulate correlation energies
        if (method == 0):
           cluster_ccsd.append(PsiMod.get_variable('CURRENT CLUSTER CCSD CORRELATION ENERGY'))
           built_ccsd += PsiMod.get_variable('CURRENT CLUSTER CCSD CORRELATION ENERGY')
        if (method == 1):
           cluster_ccsd.append(PsiMod.get_variable('CURRENT CLUSTER CEPA CORRELATION ENERGY'))
           built_ccsd += PsiMod.get_variable('CURRENT CLUSTER CEPA CORRELATION ENERGY')
        built_energy = built_ccsd
        if (lowername == 'cim-ccsd(t)'):
            cluster_ccsdt.append(PsiMod.get_variable('CURRENT CLUSTER CCSD(T) CORRELATION ENERGY'))
            cluster_t.append(PsiMod.get_variable('CURRENT CLUSTER (T) CORRELATION ENERGY'))
            built_ccsdt += PsiMod.get_variable('CURRENT CLUSTER CCSD(T) CORRELATION ENERGY')
            built_t += PsiMod.get_variable('CURRENT CLUSTER (T) CORRELATION ENERGY')
            built_energy = built_ccsdt

        cim_n += 1

    PsiMod.set_variable('CURRENT ENERGY', built_energy + escf)
    if (method == 0):
       PsiMod.set_variable('CIM-CCSD CORRELATION ENERGY', built_ccsd)
       PsiMod.set_variable('CIM-CCSD TOTAL ENERGY', built_ccsd + escf)
       PsiMod.print_out('\n')
       PsiMod.print_out('        CIM-CCSD correlation energy    %20.12lf\n' % PsiMod.get_variable('CIM-CCSD CORRELATION ENERGY'))
       PsiMod.print_out('      * CIM-CCSD total energy          %20.12lf\n' % PsiMod.get_variable('CIM-CCSD TOTAL ENERGY'))
       PsiMod.print_out('\n')
       if (lowername == 'cim-ccsd(t)'):
           PsiMod.set_variable('CIM-CCSD(T) CORRELATION ENERGY', built_ccsdt)
           PsiMod.set_variable('CIM-CCSD(T) TOTAL ENERGY', built_ccsdt + escf)
           PsiMod.set_variable('CIM-(T) CORRELATION ENERGY', built_t)
           PsiMod.print_out('        CIM-(T) correlation energy     %20.12lf\n' % PsiMod.get_variable('CIM-(T) CORRELATION ENERGY'))
           PsiMod.print_out('        CIM-CCSD(T) correlation energy %20.12lf\n' % PsiMod.get_variable('CIM-CCSD(T) CORRELATION ENERGY'))
           PsiMod.print_out('      * CIM-CCSD(T) total energy       %20.12lf\n' % PsiMod.get_variable('CIM-CCSD(T) TOTAL ENERGY'))
           PsiMod.print_out('\n')
    if (method == 1):
       PsiMod.set_variable('CIM-CEPA CORRELATION ENERGY', built_ccsd)
       PsiMod.set_variable('CIM-CEPA TOTAL ENERGY', built_ccsd + escf)
       PsiMod.print_out('\n')
       if (lowername == 'cim-cepa(0)'):
          PsiMod.print_out('        CIM-CEPA(0) correlation energy    %20.12lf\n' % PsiMod.get_variable('CIM-CEPA CORRELATION ENERGY'))
          PsiMod.print_out('      * CIM-CEPA(0) total energy          %20.12lf\n' % PsiMod.get_variable('CIM-CEPA TOTAL ENERGY'))
       if (lowername == 'cim-cepa(1)'):
          PsiMod.print_out('        CIM-CEPA(1) correlation energy    %20.12lf\n' % PsiMod.get_variable('CIM-CEPA CORRELATION ENERGY'))
          PsiMod.print_out('      * CIM-CEPA(1) total energy          %20.12lf\n' % PsiMod.get_variable('CIM-CEPA TOTAL ENERGY'))
       if (lowername == 'cim-cepa(2)'):
          PsiMod.print_out('        CIM-CEPA(2) correlation energy    %20.12lf\n' % PsiMod.get_variable('CIM-CEPA CORRELATION ENERGY'))
          PsiMod.print_out('      * CIM-CEPA(2) total energy          %20.12lf\n' % PsiMod.get_variable('CIM-CEPA TOTAL ENERGY'))
       if (lowername == 'cim-cepa(3)'):
          PsiMod.print_out('        CIM-CEPA(3) correlation energy    %20.12lf\n' % PsiMod.get_variable('CIM-CEPA CORRELATION ENERGY'))
          PsiMod.print_out('      * CIM-CEPA(3) total energy          %20.12lf\n' % PsiMod.get_variable('CIM-CEPA TOTAL ENERGY'))
       if (lowername == 'cim-cisd'):
          PsiMod.print_out('        CIM-CISD correlation energy    %20.12lf\n' % PsiMod.get_variable('CIM-CEPA CORRELATION ENERGY'))
          PsiMod.print_out('      * CIM-CISD total energy          %20.12lf\n' % PsiMod.get_variable('CIM-CEPA TOTAL ENERGY'))
       if (lowername == 'cim-aqcc'):
          PsiMod.print_out('        CIM-AQCC correlation energy    %20.12lf\n' % PsiMod.get_variable('CIM-CEPA CORRELATION ENERGY'))
          PsiMod.print_out('      * CIM-AQCC total energy          %20.12lf\n' % PsiMod.get_variable('CIM-CEPA TOTAL ENERGY'))
       if (lowername == 'cim-acpf'):
          PsiMod.print_out('        CIM-ACPF correlation energy    %20.12lf\n' % PsiMod.get_variable('CIM-CEPA CORRELATION ENERGY'))
          PsiMod.print_out('      * CIM-ACPF total energy          %20.12lf\n' % PsiMod.get_variable('CIM-CEPA TOTAL ENERGY'))
       PsiMod.print_out('\n')

    return built_energy + escf

# the parallel version:
#def run_plugin_libcim_parallel(name, **kwargs):
#    """Function encoding sequence of PSI module and plugin calls so that
#    Eugene's plugin_libcim can be called via :py:func:`driver.energy`.
#
#    >>> energy('cim-ccsd(t)')
#
#    """
#    lowername = name.lower()
#    kwargs = kwargs_lower(kwargs)
#
#    # create a MPI communicator
#    comm = MPI.COMM_WORLD
#    # get the rank
#    rank = comm.Get_rank()
#    # get the number of processes
#    nproc = comm.Get_size()
#
#    # override symmetry:
#    molecule = PsiMod.get_active_molecule()
#    molecule.update_geometry()
#    molecule.reset_point_group('c1')
#    molecule.fix_orientation(1)
#    molecule.update_geometry()
#
#    # what type of cim?
#    if (lowername == 'cim-ccsd'):
#        PsiMod.set_global_option('compute_triples', False)
#    if (lowername == 'cim-ccsd(t)'):
#        PsiMod.set_global_option('compute_triples', True)
#
#    # some options are not correct when i load two plugins ... set them here
#    PsiMod.set_global_option('r_convergence', 1e-7)
#    PsiMod.set_global_option('maxiter', 100)
#
#    energy('scf', **kwargs)
#    PsiMod.set_global_option('cim_initialize', True)
#    PsiMod.plugin('plugin_libcim.so')
#
#    cluster_ccsd = []
#    cluster_ccsdt = []
#    cluster_t = []
#    for i in range(nproc):
#        cluster_ccsd.append(0.0)
#        cluster_ccsdt.append(0.0)
#        cluster_t.append(0.0)
#
#    built_ccsd = 0.0
#    built_t = 0.0
#    built_ccsdt = 0.0
#    built_energy = 0.0
#    escf = PsiMod.get_variable('SCF TOTAL ENERGY')
#    PsiMod.set_global_option('cim_initialize', False)
#    cim_n = 0
#    while cim_n < PsiMod.get_variable('CIM CLUSTERS'):
#        # if rank owns this work
#        if (rank == cim_n%nproc):
#           # run plugin_libcim
#           PsiMod.set_global_option('CIM_CLUSTER_NUM', cim_n)
#           PsiMod.plugin('plugin_libcim.so')
#           # accumulate correlation energies
#           cluster_ccsd[rank] += PsiMod.get_variable('CURRENT CLUSTER CCSD CORRELATION ENERGY')
#           if (lowername == 'cim-ccsd(t)'):
#               cluster_ccsdt[rank] += PsiMod.get_variable('CURRENT CLUSTER CCSD(T) CORRELATION ENERGY')
#               cluster_t[rank] += PsiMod.get_variable('CURRENT CLUSTER (T) CORRELATION ENERGY')
#        cim_n += 1
#
#    PsiMod.flush_outfile()
#
#    # grab energies from other all procs:
#    for i in range(nproc):
#        if (i>0):
#           if (rank == i):
#              comm.send(cluster_ccsd[i], dest=0, tag=i+1)
#           if (rank == 0):
#              cluster_ccsd[i] = comm.recv(source=i, tag=i+1)
#   
#           if (rank == i):
#              comm.send(cluster_ccsdt[i], dest=0, tag=i+1)
#           if (rank == 0):
#              cluster_ccsdt[i] = comm.recv(source=i, tag=i+1)
#   
#           if (rank == i):
#              comm.send(cluster_t[i], dest=0, tag=i+1)
#           if (rank == 0):
#              cluster_t[i] = comm.recv(source=i, tag=i+1)
#
#        built_ccsd += cluster_ccsd[i]
#        built_ccsdt += cluster_ccsdt[i]
#        built_t += cluster_t[i]
#
#    built_energy = built_ccsd
#    if (lowername == 'cim-ccsd(t)'):
#       built_energy = built_ccsdt
#
#    PsiMod.set_variable('CURRENT ENERGY', built_energy + escf)
#    PsiMod.set_variable('CIM-CCSD CORRELATION ENERGY', built_ccsd)
#    PsiMod.set_variable('CIM-CCSD TOTAL ENERGY', built_ccsd + escf)
#    if (lowername == 'cim-ccsd(t)'):
#        PsiMod.set_variable('CIM-CCSD(T) CORRELATION ENERGY', built_ccsdt)
#        PsiMod.set_variable('CIM-CCSD(T) TOTAL ENERGY', built_ccsdt + escf)
#        PsiMod.set_variable('CIM-(T) CORRELATION ENERGY', built_t)
#
#    if (rank == 0):
#       PsiMod.print_out('\n')
#       PsiMod.print_out('        CIM-CCSD correlation energy    %20.12lf\n' % PsiMod.get_variable('CIM-CCSD CORRELATION ENERGY'))
#       PsiMod.print_out('      * CIM-CCSD total energy          %20.12lf\n' % PsiMod.get_variable('CIM-CCSD TOTAL ENERGY'))
#       PsiMod.print_out('\n')
#       if (lowername == 'cim-ccsd(t)'):
#          PsiMod.print_out('        CIM-(T) correlation energy     %20.12lf\n' % PsiMod.get_variable('CIM-(T) CORRELATION ENERGY'))
#          PsiMod.print_out('        CIM-CCSD(T) correlation energy %20.12lf\n' % PsiMod.get_variable('CIM-CCSD(T) CORRELATION ENERGY'))
#          PsiMod.print_out('      * CIM-CCSD(T) total energy       %20.12lf\n' % PsiMod.get_variable('CIM-CCSD(T) TOTAL ENERGY'))
#          PsiMod.print_out('\n')
#
#    return built_energy + escf

# Integration with driver routines
procedures['energy']['cim-ccsd(t)'] = run_plugin_libcim
procedures['energy']['cim-ccsd']    = run_plugin_libcim
procedures['energy']['cim-cepa(0)'] = run_plugin_libcim
procedures['energy']['cim-cepa(1)'] = run_plugin_libcim
procedures['energy']['cim-cepa(2)'] = run_plugin_libcim
procedures['energy']['cim-cepa(3)'] = run_plugin_libcim
procedures['energy']['cim-cisd']    = run_plugin_libcim
procedures['energy']['cim-acpf']    = run_plugin_libcim
procedures['energy']['cim-aqcc']    = run_plugin_libcim

def exampleFN():
    # Your Python code goes here
    pass

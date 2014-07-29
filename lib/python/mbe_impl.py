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
from parallel import *
from qmmm import *
import sys

#Our instance of the MBE c++ API
lfrag=psi4.LibFragHelper()

       

def RunCalc(name,molecule,atoms,caps,charges,ghosts,Egys,**kwargs):
    old_geom=molecule.save_string_xyz()
    temp_geom=old_geom.split('\n')
    new_geom=[]
    qmmm=QMMM()
    for i in range(0,len(charges)):
        qmm.addChargeBohr(charges[i*4],charges[i*4+1],
                          charges[i*4+2],charges[i*4+3])
    if(len(charges)>=1):
        qmmm.populateExtern()
        psi4.set_global_option_python('EXTERN',qmmm.extern)
    for i in range(0,len(atoms)):
        new_geom.append(temp_geom[atoms[i]+1])
    new_geom.append(caps)
    for i in range(0,len(ghosts)):
        new_geom.append('@'+temp_geom[ghosts[i]+1].lstrip(' '))
    new_geom.append('symmetry c1')
    new_geom.append('no_com')
    new_geom.append('no_reorient')
    fragstr='\n'.join(new_geom)
    frag=geometry(fragstr)
    activate(frag)
    frag.update_geometry()
    Egys.append(energy(name,**kwargs))
    activate(molecule)
    molecule.update_geometry

#This call is incharge of setting up basic MPI stuff
def DoMPI(PMan,size): 
    current_comm=PMan.CurrentComm()
    MyProcN=PMan.me()
    NProc=PMan.nproc()
    PMan.make_comm("NewComm",MyProcN)
    if NProc > 1:
        batchsize_remainder=size%NProc
        batchsize=(size-batchsize_remainder)/NProc
        mystart=MyProcN*batchsize
        myend=(MyProcN+1)*batchsize
    else:
        batchsize_remainder=0
        mystart=0
        myend=size
    return mystart,myend,batchsize_remainder        



#This is where the magic happens
def Magic(frag,N,Egys,name,molecule,**kwargs):
    atoms=lfrag.GetNMerN(N,frag)
    caps=lfrag.CapHelper(N,frag)
    charges=[]#lfrag.EmbedHelper("NNNNO")
    ghosts=lfrag.GetGhostsNMerN(N,frag)
    oldguess=psi4.get_global_option("GUESS")    
    if(N!=0):
        lfrag.WriteMOs(N,frag)
        psi4.set_global_option("GUESS","READ")
    RunCalc(name,molecule,atoms,caps,charges,ghosts,Egys,**kwargs)
    if(N==0):
        lfrag.ReadMOs()
    psi4.set_global_option("GUESS",oldguess)
    psi4.clean()
    

#Both the fragments and NMers go through the same sequence of calls,
#we term them the basecall
def BaseCall(name,molecule,size,Egys,i,**kwargs):
    PMan=Parallel() 
    old_comm=PMan.CurrentComm()
    PMan.sync(old_comm)
    psi4.be_quiet()
    [mystart,myend,remainder]=DoMPI(PMan,size)
    for x in range(mystart,myend):
        Magic(x,i,Egys,name,molecule,**kwargs)
    PMan.all_gather(Egys,old_comm)
    offset=size-remainder
    for x in range(0,remainder):
        Egys.append(0)
    #Last element in the list is the value
    for x in range(0,remainder):
        if PMan.me(old_comm)==x:
            Magic(x+offset,i,Egys,name,molecule,**kwargs)
    PMan.sync(old_comm)
    for x in range(0,remainder):
        #Has to be a list...
        value=[]
        if PMan.me(old_comm)==x:
            value.append(Egys.pop())
        else:
            value.append(0)
        PMan.bcast(value,x,old_comm)
        Egys[offset+x]=value[0]
    lfrag.Sync(old_comm,i)
    PMan.free_comm(PMan.CurrentComm())
    psi4.reopen_outfile()
    PMan.sync(old_comm)

def setup(frag_method,N,embed_method,cap_method,bsse_method):
    lcfrag_method=frag_method.lower()
    lcbsse_method=bsse_method.lower()
    lcembed_method=embed_method.lower()
    lccap_method=cap_method.lower()
    lfrag.FragHelper(lcfrag_method,N,lcembed_method,lccap_method,
                     lcbsse_method)

def fragment(name,molecule,Egys,**kwargs):
    size=lfrag.GetNFrags()
    if(lfrag.RunFrags()!=0):
        psi4.print_out('The system contains %d monomers\n'% (size))
        BaseCall(name,molecule,size,Egys,0,**kwargs)
        psi4.print_out('Finished running monomer calculations\n')
    
def nmers(name,molecule,N,Egys,**kwargs):
    lfrag.NMerHelper(N)
    i=1
    while lfrag.GetNNMers(i)!=0:
        NNMers=lfrag.GetNNMers(i)
        psi4.print_out('The system contains %d %d-mers\n'% (NNMers,i+1))
        Egys.append([])
        BaseCall(name,molecule,NNMers,Egys[i],i,**kwargs)
        psi4.print_out('Finished running %d-body calculations\n'% (i+1))
        i+=1
       
def SystemEnergy(Egys):
    return lfrag.CalcEnergy(Egys)

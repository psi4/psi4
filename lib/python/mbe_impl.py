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

def Done():
    global lfrag
    del lfrag


def RunCalc(name,molecule,atoms,caps,charges,ghosts,Egys,**kwargs):
    old_geom=molecule.save_string_xyz()
    temp_geom=old_geom.split('\n')
    new_geom=[]
    qmmm=QMMM()
    for i in range(0,len(charges)/4):
        qmmm.addChargeBohr(charges[i*4],charges[i*4+1],
                              charges[i*4+2],charges[i*4+3])
    if(len(charges)>=4):
        qmmm.populateExtern()
        psi4.set_global_option_python('EXTERN',qmmm.extern)
    for i in range(0,len(atoms)):
        new_geom.append(temp_geom[atoms[i]+1])
    new_geom.append(caps)
    for i in range(0,len(ghosts)):
        new_geom.append('@'+temp_geom[ghosts[i]+1].lstrip(' '))
    new_geom.append('symmetry c1')
    new_geom.append('no_reorient')
    new_geom.append('no_com')
    fragstr='\n'.join(new_geom)
    frag=geometry(fragstr)
    activate(frag)
    frag.update_geometry()
    Egys.append(energy(name,**kwargs))
    if(len(charges)>=4):
        psi4.set_global_option_python("EXTERN", None)
    oep=psi4.OEProp()
    oep.add("MULLIKEN_CHARGES")
    oep.compute()
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
    charges=lfrag.EmbedHelper(N,frag)
    ghosts=lfrag.GetGhostsNMerN(N,frag)
    oldguess=psi4.get_global_option("GUESS")
    if((N!=0 and lfrag.IsGMBE()==0) or len(charges)>=4):
        lfrag.WriteMOs(N,frag)
        psi4.set_global_option("GUESS","READ")
    RunCalc(name,molecule,atoms,caps,charges,ghosts,Egys,**kwargs)
    if(N==0):
        lfrag.GatherData()
        psi4.set_global_option("GUESS",oldguess)
    psi4.clean()
        

#Both the fragments and NMers go through the same sequence of calls,
#we term them the basecall
def BaseCall(name,molecule,size,Egys,i,SuppressPrint,itr,**kwargs):
    PMan=Parallel() 
    old_comm=PMan.CurrentComm()
    PMan.sync(old_comm)
    if(SuppressPrint):
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
    lfrag.Sync(old_comm,i,itr)
    PMan.free_comm(PMan.CurrentComm())
    if(SuppressPrint):
        psi4.reopen_outfile()
    PMan.sync(old_comm)

def setup(frag_method,N,embed_method,cap_method,bsse_method):
    lcfrag_method=frag_method.lower()
    lcbsse_method=bsse_method.lower()
    lcembed_method=embed_method.lower()
    lccap_method=cap_method.lower()
    lfrag.FragHelper(lcfrag_method,N,lcembed_method,lccap_method,
                     lcbsse_method)

def fragment(name,molecule,Egys,SuppressPrint,**kwargs):
    size=lfrag.GetNFrags()
    if(lfrag.RunFrags()!=0):
        psi4.print_out('The system contains %d monomers\n'% (size))
        itr=0;
        done =False
        while (not done): 
            if(itr==1):
                psi4.print_out('Rerunning monomer\'s with embedding\n')
            if(itr>1):
                psi4.print_out('Embedding Iteration %d \n'%itr)
            if(itr>=1):
                #Reset energies
                del Egys[:]
            BaseCall(name,molecule,size,Egys,0,SuppressPrint,itr,**kwargs)
            if(lfrag.Iterate(itr)==0):
                done=True
            itr=itr+1
        psi4.print_out('Finished running monomer calculations\n')
    
def nmers(name,molecule,N,Egys,SuppressPrint,**kwargs):
    lfrag.NMerHelper(N)
    i=1
    while lfrag.GetNNMers(i)!=0:
        NNMers=lfrag.GetNNMers(i)
        psi4.print_out('The system contains %d %d-mers\n'% (NNMers,i+1))
        Egys.append([])
        BaseCall(name,molecule,NNMers,Egys[i],i,SuppressPrint,0,**kwargs)
        psi4.print_out('Finished running %d-body calculations\n'% (i+1))
        i+=1
       
def SystemEnergy(Egys):
    return lfrag.CalcEnergy(Egys)

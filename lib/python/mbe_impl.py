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
import sys

#Our instance of the MBE c++ API
lfrag=psi4.LibFragHelper()

       
def RunCalc(name,molecule,atoms,Egys,**kwargs):
    old_geom=molecule.save_string_xyz()
    temp_geom=old_geom.split('\n')
    new_geom=[]
    for i in range(0,len(atoms)):
        new_geom.append(temp_geom[atoms[i]+1])
    fragstr='\n'.join(new_geom)
    frag=geometry(fragstr)
    activate(frag)
    frag.update_geometry()
    Egys.append(energy(name,**kwargs))
    psi4.clean()
    activate(molecule)
    molecule.update_geometry

#This call is incharge of setting up basic MPI stuff
def DoMPI(PMan,size): 
    current_comm=PMan.current_comm
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

#Both the fragments and NMers go through the same sequence of calls,
#we term them the basecall
def BaseCall(name,molecule,size,Egys,i,**kwargs):
    PMan=Parallel() 
    old_comm=PMan.current_comm
    [mystart,myend,remainder]=DoMPI(PMan,size)
    for x in range(mystart,myend):
        atoms=lfrag.GetNMerN(i,x)
        psi4.be_quiet()
        RunCalc(name,molecule,atoms,Egys,**kwargs)
        psi4.reopen_outfile()
    PMan.all_gather(Egys,old_comm)
    offset=size-remainder
    for x in range(0,remainder):
        Egys.append(0)
    #Last element in the list is the value
    for x in range(0,remainder):
        if PMan.me(old_comm)==x:
            atoms=lfrag.GetNMerN(i,offset+x)
            psi4.be_quiet()
            RunCalc(name,molecule,atoms,Egys,**kwargs)
            psi4.reopen_outfile()
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
    PMan.sync(old_comm)
    PMan.free_comm(PMan.current_comm)

def fragment(name,molecule,frag_method,Egys,**kwargs):
    lcfrag_method=frag_method.lower()
    lfrag.FragHelper(lcfrag_method)
    size=lfrag.GetNFrags()
    BaseCall(name,molecule,size,Egys,0,**kwargs)
    
def nmers(name,molecule,N,Egys,**kwargs):
    lfrag.NMerHelper(N)
    i=1
    while lfrag.GetNNMers(i)!=0:
        NNMers=lfrag.GetNNMers(i)
        Egys.append([])
        BaseCall(name,molecule,NNMers,Egys[i],i,**kwargs)
        i+=1
       
def SystemEnergy(Egys):
    return lfrag.CalcEnergy(Egys)
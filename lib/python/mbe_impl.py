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
    
def MakePointChargeEmbed(N,frag):
    charges=lfrag.EmbedHelper(N,frag)
    qmmm=QMMM()
    for i in range(0,len(charges)/4):
        qmmm.addChargeBohr(charges[i*4],charges[i*4+1],charges[i*4+2],charges[i*4+3])
    if(len(charges)>=4):
        qmmm.populateExtern()
        psi4.set_global_option_python('EXTERN',qmmm.extern)
    return (len(charges)>=4)

def NewMolecule(N,frag,molecule):
    atoms=lfrag.GetNMerN(N,frag)
    caps=lfrag.CapHelper(N,frag)
    ghosts=lfrag.GetGhostsNMerN(N,frag)
    old_geom=molecule.save_string_xyz()
    temp_geom=old_geom.split('\n')
    new_geom=[]
    for i in range(0,len(atoms)):
        new_geom.append(temp_geom[atoms[i]+1])
    new_geom.append(caps)
    for i in range(0,len(ghosts)):
        new_geom.append('@'+temp_geom[ghosts[i]+1].lstrip(' '))
    new_geom.append('symmetry c1')
    new_geom.append('no_reorient')
    new_geom.append('no_com')
    fragstr='\n'.join(new_geom)
    return geometry(fragstr)

def GetOptions():
    DaOptions={}
    DaOptions['DF_BASIS_SCF']=psi4.get_global_option('DF_BASIS_SCF')
    DaOptions['DF_BASIS_MP2']=psi4.get_global_option('DF_BASIS_MP2')
    DaOptions['DF_BASIS_CC']=psi4.get_global_option('DF_BASIS_CC')
    return DaOptions

def UpdateOptions(DaOptions):
    for k,v in DaOptions.iteritems():
        if(not ((DaOptions[k] == '') or (DaOptions[k] == 'NONE'))):
            psi4.set_global_option(k,DaOptions[k])

def RunCalc(frag,name,molecule,Egys,CEgys,N,**kwargs):
    OldOptions=GetOptions()
    oldguess=psi4.get_global_option("GUESS")
    if(lfrag.WriteMOs(N,frag)==1):        
        psi4.set_global_option("GUESS","READ")
    AreEmbedding=MakePointChargeEmbed(N,frag)
    frag=NewMolecule(N,frag,molecule)
    activate(frag)
    frag.update_geometry()
    UpdateOptions(OldOptions)
    Egys.append(energy(name,**kwargs))
    for k,v in CEgys.iteritems():
        CEgys[k].append(psi4.get_variable(k))
    if(AreEmbedding):
        psi4.set_global_option_python("EXTERN", None)
    oep=psi4.OEProp()
    oep.add("MULLIKEN_CHARGES")
    oep.compute()
    activate(molecule)
    molecule.update_geometry()
    if(N==0):
        lfrag.GatherData()
    psi4.set_global_option("GUESS",oldguess)
    psi4.clean()    

#Both the fragments and NMers go through the same sequence of calls,
#we term them the basecall

def BaseCall(name,molecule,Egys,CEgys,NStart,NEnd,SuppressPrint,itr,**kwargs):
    if(SuppressPrint):
        psi4.be_quiet()
    PMan=Parallel()
    sizes=[] 
    for N in range(NStart,NEnd):
        sizes.append(lfrag.GetNNMers(N))
        for frag in range(0,sizes[N-NStart]):
            PMan.AddTask(str(N)+' '+str(frag),(NEnd-1)-N)
    PMan.MakeJob()
    task=PMan.Begin()
    tempEgys=[]
    tempCEgys={}
    tempCEgys.fromkeys(CEgys.keys(),[])
    while(not PMan.Done()):
        frag=int(task.split(' ')[1])
        order=int(task.split(' ')[0])
        RunCalc(frag,name,molecule,tempEgys,tempCEgys,order,**kwargs)
        task=PMan.Next()
    R1=PMan.Synch(tempEgys,1)
    R2={}
    R2.fromkeys(CEgys.keys(),[])
    for k,v in tempCEgys.iteritems():
        R2[k]=PMan.Synch(tempCEgys[k],1)
    total=0
    print sizes
    for N in range(NStart,NEnd):
        offset=sizes[N-NStart]
        lfrag.Sync("COMM_WORLD",N,itr)
        Egys[N]=R1[total:total+offset]
        for k,v in tempCEgys.iteritems():
            CEgys[k][N]=R2[k][total:total+offset]
        total+=offset
    del PMan
    if(SuppressPrint):
        psi4.reopen_outfile()

def setup(frag_method,N,embed_method,cap_method,bsse_method,IsSymm):   
    lcfrag_method=frag_method.lower()
    lcbsse_method=bsse_method.lower()
    lcembed_method=embed_method.lower()
    lccap_method=cap_method.lower()
    lfrag.FragHelper(lcfrag_method,N,lcembed_method,lccap_method,
                     lcbsse_method,IsSymm)

def fragment(name,molecule,Egys,CEgys,SuppressPrint,**kwargs):
    if(lfrag.RunFrags()!=0):
        itr=0;
        done =False
        while (not done): 
            if(itr>=1):
                #Reset energies
                del Egys[0][:]
                for k,v in CEgys.iteritems():
                    del CEgys[v][0][:]
            BaseCall(name,molecule,Egys,CEgys,0,1,SuppressPrint,itr,**kwargs)
            if(lfrag.Iterate(itr)==0):
                done=True
            itr=itr+1
        lfrag.PrintEnergies(Egys,1,'Total Energy')

    
def nmers(name,molecule,N,Egys,CEgys,SuppressPrint,**kwargs):
    lfrag.NMerHelper(N)
    PMan=Parallel()
    NStart=1
    NEnd=NStart
    while lfrag.GetNNMers(NEnd)!=0:
        NEnd+=1  
        Egys.append([])
        for k,v in CEgys.iteritems():
            CEgys[k].append([])
    BaseCall(name,molecule,Egys,CEgys,NStart,NEnd,SuppressPrint,0,**kwargs)
    for i in range(NStart,NEnd):
        lfrag.PrintEnergies(Egys,i+1,'Total Energy')
       
def SystemEnergy(Egys,IsCorr):
    return lfrag.CalcEnergy(Egys,IsCorr)

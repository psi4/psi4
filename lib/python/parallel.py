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
from driver import *

class Parallel:
    """Parallel is the base class for managing MPI communication"""
    def CurrentComm(self):
        return psi4.get_comm()
    
    def GetComm(self,Comm):
        """ Returns a communicator based on input.
        
        If you pass NONE, this will return current_comm.
        Otherwise it returns whatever you passed.
        Basically I use this all over the place and don't want
        to type it multiple times.  Oh, and this function should
        be considered private, so don't use it outside the class
        """
        if Comm=="NONE":
            return self.CurrentComm()
        return Comm
    
    def __init__(self):
        """ Makes a Parallel class, which amounts to figuring out the current comm"""
        #self.current_comm[0]=psi4.get_comm()
    
    def make_comm(self,NewComm,Color,Comm2Split="NONE"):
        """This splits Comm2Split, base on color, and calls the new comm NewComm."""
        psi4.make_comm(NewComm,Color,self.GetComm(Comm2Split))
        #self.current_comm[0]=psi4.get_comm()
    
    def free_comm(self,Comm="NONE"):
        """Frees the communicator passed to it.  No check to ensure COMM_WORLD isn't freed"""
        psi4.free_comm(self.GetComm(Comm))
        #self.current_comm[0]=psi4.get_comm()
    
    def me(self,Comm="NONE"):
        """Returns the current process number for the communicator Comm"""
        return psi4.me(self.GetComm(Comm))
    
    def nproc(self,Comm="NONE"):
        """Returns the number of processes in comm Comm"""
        return psi4.nproc(self.GetComm(Comm))
    
    def all_gather(self,data,Comm="NONE"):
        """Gathers lists of doubles.  Only lists and only doubles.
        
        Assume you have N processes, each with a list called data.
        data, needs to be the same length on each process, call it
        n (for those of you versed in MPI, this is not MPI all 
        gatherv, sorry). In the end is N*n elements long where
        elements i*n through (i+1)*n came from process i (0..i..N-1)
        """ 
        target=psi4.all_gather_double(data,self.GetComm(Comm))
        del data[0:len(data)]
        for i in range(0,len(target)):
            data.append(target[i])
    
    def bcast(self,data,broadcaster,Comm="NONE"):
        """Given a list of doubles, broadcasts the list to everyone.
        
        Everyones list should be the proper size beforehand.
        
        Process number broadcaster, takes the data contained in data,
        and sends it out to everyone on comm Comm.  Afterwards, everyones
        data is the same.  I am unsure of the boost::python list symmenatics
        so I do this copy thing.
        """  
        value=[]
        tempcomm=self.GetComm(Comm)
        for i in range(0,len(data)):
            if self.me(tempcomm)==broadcaster:
                value.append(data[i])
            else:
                value.append(0)
        temp=psi4.bcast_double(value,broadcaster,tempcomm)
        for i in range(0,len(data)):
            if self.me(tempcomm)!=broadcaster:
                data[i]=temp[i]
        
    
    def sync(self,Comm="NONE"):
        """Makes all processes in comm Comm wait at this point until they are all here"""
        psi4.sync(self.GetComm(Comm))
    
    def debug_print(self,thing2print,Comm="NONE"):
        """Useful printing function that tells you what each process has"""
        nproc_var=self.nproc(self.GetComm(Comm))
        for x in range(0,nproc_var):
            me_var=self.me(self.GetComm(Comm))
            if me_var==x:
                print("I am process: ",me_var," and I contain:")
                print(thing2print)
                sys.stdout.flush()
            self.sync(self.GetComm(Comm))
 
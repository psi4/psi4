from __future__ import absolute_import
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
import sys


class Parallel:
    
    """Parallel is the base class for managing MPI communication"""
    
    def __init__(self):
        self.PGuts=psi4.LibParallelHelper()
    
    def AddTask(self,name,priority):
        self.PGuts.AddTask(name,priority)

    def MakeJob(self):
        self.PGuts.MakeJob()

    def Begin(self):
        return self.PGuts.Begin()

    def Next(self):
        return self.PGuts.Next()

    def Done(self):
        return self.PGuts.Done()

    def Synch(self,local,N):
        return self.PGuts.Synch(local,N)

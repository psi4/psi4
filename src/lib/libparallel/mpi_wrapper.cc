/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

#include "mpi_wrapper.h"

#include <omp.h>
#if HAVE_MPI
namespace psi{

MPICommunicator::MPICommunicator(int &argc,char **argv):
        Env(new boost::mpi::environment(argc,argv)){
        boost::mpi::communicator world;
        Communicators[CurrentComm.back()]=world;
        //The next three lines are what the old local comm did
        //the code breaks if I do not include them
        //The way I understand this is that the number of openmp threads
        //is getting hard-coded to 1, and somewhere in the code people are
        //counting on this behavior...
        omp_set_nested(0);
        if (Process::environment("OMP_NUM_THREADS") == "")
            Process::environment.set_n_threads(1);
}

void MPICommunicator::MakeComm(const std::string& Name,const int Color,const std::string& Comm2Split){
    boost::mpi::communicator comm=GetComm(Comm2Split);
    Communicators[Name]=comm.split(Color);
    CurrentComm.push_back(Name);
}

void MPICommunicator::FreeComm(const std::string& Name){
    //Refuse to free mpi_comm_world
    std::cout<<"C++: "<<Name<<std::endl;
    if(CurrentComm.size()>1){
        Communicators.erase(Name);
        CurrentComm.erase(CurrentComm.end()-1);
        for(int i=0;i<CurrentComm.size();i++)
           std::cout<<"CurrentComm "<<i<<" : "<<CurrentComm[i]<<std::endl;
    }
}

}//End namespace psi


#endif






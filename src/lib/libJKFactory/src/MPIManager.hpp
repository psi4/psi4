/*
 * JKFactory: Interface and code for highly parallel J and K
 *             builds.
 *
 *  Copyright (c) 2014 Ryan M. Richard
 *
 *  This file is part of JKFactory.
 *
 *  JKFactory is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef MPIMANAGER_HPP_
#define MPIMANAGER_HPP_

#include <boost/mpi.hpp>
#include <boost/shared_ptr.hpp>
#include <map>
/*** \brief The class that manages our MPI session
 *
 *  Perhaps the most important aspect of this class is making sure you set
 *  the communicator appropriately so that it reflects the state of your
 *  program prior to the JKFactory call.  Also note that it is the only class
 *  that does not have to derive from JKFactoryBase, because it is part of it.
 *
 *  This class is very similar to the one I wrote for Psi4, but is slimmed down.
 */
namespace JKFactory{
class MPIManager{
	private:
		///This is the actual MPI environment
		//boost::shared_ptr<boost::mpi::environment> Env;
		///This is a mapping between Comm names and Comms
		std::map<std::string,boost::mpi::communicator> Comms;
		///This is the "que" of comms, i.e. the order they were derived in
		std::vector<std::string> CommList;
	public:
        MPIManager(){
           CommList.push_back("COMM_WORLD");
        }
        void Initialize(const MPI_Comm& BaseComm){
           boost::mpi::communicator temp(BaseComm,boost::mpi::comm_attach);
           Comms[CommList.back()]=temp;
        }

		std::string Comm(){return CommList[CommList.size()-1];}
		///Wait for all processes to catch up
		void Sync(const std::string& Comm){Comms[Comm].barrier();}
		///The number of processes on this communicator
		int NProc(const std::string& Comm){return Comms[Comm].size();}
		///My identity
		int Me(const std::string& Comm){return Comms[Comm].rank();}
		///Does an All_Gather
		template<class T>
		void AllGather(T* Data,int NElem,T* Target,const std::string Comm){
				boost::mpi::all_gather(Comms[Comm],Data,NElem,Target);
		}
		///Makes NewComm by splitting Comm2Split by color, and ordering by rank
		void MakeComm(const std::string& Comm2Split,int color,int rank,
				const std::string& NewComm){
			Comms[NewComm]=Comms[Comm2Split].split(color,rank);
		}
		template<class T>
		void Bcast(T* Data,int NElem,int Broadcaster,const std::string& Comm){
		      boost::mpi::broadcast(Comms[Comm],Data,NElem,Broadcaster);
		}

};

}

#endif /* MPIMANAGER_HPP_ */

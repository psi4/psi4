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

#include "ParallelComm.h"
#include <boost/mpi/communicator.hpp>

namespace psi{
namespace LibParallel{
typedef boost::shared_ptr<ParallelComm> SharedThis;
typedef boost::mpi::communicator BComm;
ParallelComm::ParallelComm():BaseType_(this),Comm_(new BComm()){}

SharedThis ParallelComm::MakeComm(const int Color)const{
   SharedThis temp(new ParallelComm());
   temp->Comm_=
         boost::shared_ptr<BComm>(new BComm(this->Comm_->split(Color)));
   return temp;
}

void ParallelComm::Barrier()const{
   Comm_->barrier();
}

int ParallelComm::Me()const{
   return Comm_->rank();
}

int ParallelComm::NProc()const{
   return Comm_->size();
}

int ParallelComm::Probe(const int Sender,const int MessageTag,
                        const bool Block)const{
   int Sendcopy=(Sender<0?boost::mpi::any_source:Sender);
   int Messagecopy=(MessageTag<0?boost::mpi::any_tag:MessageTag);
   boost::optional<boost::mpi::status> status=
         (Block?Comm_->probe(Sendcopy,Messagecopy):
                Comm_->iprobe(Sendcopy,Messagecopy));
   return (status?status.get().source():-1);
}

}}//End namespaces


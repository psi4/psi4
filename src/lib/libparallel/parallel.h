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

/*
 * File:   parallel.h
 * Author: jturney, jjwilke
 *
 * Created on December 11, 2009, 3:34 PM
 */

#ifndef _psi_src_lib_libparallel_parallel_h_
#define	_psi_src_lib_libparallel_parallel_h_

#include <boost/shared_ptr.hpp>
#include <psiconfig.h>
#include <cstdio>

// Define the default parallel environment. In general, developers should just use this type.
#if defined(HAVE_MPI)
#   include <libparallel/mpi_wrapper.h>
typedef psi::MPICommunicator            worldcomm;
#else
#   include <libparallel/local.h>
typedef psi::LocalCommWrapper           worldcomm;
#endif

namespace psi {

    extern FILE *outfile;


    // A templated version of init comunicator.
    template <typename comm_type>
    static comm_type* initialize_specific_communicator(const int &argc, char **argv) {
        return new comm_type(argc, argv);
    }

    // Create a communicator from Comm typedef'ed above.
    boost::shared_ptr<worldcomm> initialize_communicator(const int &argc, char **argv);
}
/*
template <typename DerivedClass>
class Parallel{
	public:
	   ///This typedef is the type of this class
	   typedef Parallel<DerivedType> ThisType;
	protected:
	   ///Array of our communicator names, in the order they are derived
	   std::vector<std::string> CurrentComm;
	public:
	   Parallel(){CurrentComm.push_back("COMM_WORLD");}
	   ~Parallel(){}
	   virtual void sync(const std::string& CommName="None") const{}

	   template<class T>
	   void all_gather(const T* localdata,const int nelem,T* target,
	                   const std::string& CommName="None") const{
	      static_cast<const DerivedType*>(this)->all_gatherImpl(data,
	          nelem,target,CommName);
	   }

	   template<class T>
	   void bcast(T* data,const int nelem,const int broadcaster,
	              const std::string& CommName="NONE")const{
			static_cast<const DerivedType*>(this)->bcastImpl(data,nelem,
			     broadcaster,CommName);
		}

		virtual int me(const std::string& CommName="None") const{
		    return 0
		};

		virtual int nproc(const std::string& CommName="None) const{}

private:
    int me_;
    int nproc_;
    int nthread_;
    std::string communicator_;

public:

    LocalCommWrapper(const int &argc, char **argv)
    {
        init_openmp();
        me_ = 0;
        nproc_ = 1;
        nthread_ = Process::environment.get_n_threads();
        communicator_ = "LOCAL";
    }

    ~LocalCommWrapper() { }


    inline int nproc() const { return nproc_; }


    inline int nthread() const { return nthread_; }


    inline std::string communicator() const {return communicator_;}


    inline int thread_id(const pthread_t &thread) { return 0; }



    template<typename type>
    inline void send(const type data, int nelem, int target) { }


    template<typename type>
    inline void recv(type data, int nelem, int sender) { }



    template<typename type>
    inline void bcast_serializable(type &data, int broadcaster=0) { }


    template<typename type>
    inline void sum(type data, int nelem, type *receive_buffer=0, int target=-1)
    {
        if (receive_buffer != 0)
            std::memcpy(receive_buffer, data, sizeof(type) * nelem);
    }


    inline void print(FILE *out) const
    {
        fprintf(out, "\n    Using LocalCommunicator (Number of processes = 1)\n\n");
    }


    inline void finalize()
    { }

//#ifdef HAVE_MADNESS
//    SharedLock get_spinlock() { }
//    SharedMutex get_mutex() { }
//    SharedMadWorld get_madworld() { }
//#endif

}; // End of LocalCommWrapper class
*/
//#include "threaded_storage.h"

#endif  /* _psi_src_lib_libparallel_parallel_h_ */

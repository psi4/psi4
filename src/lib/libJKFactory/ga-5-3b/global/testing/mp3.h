#if defined(PVM)
#   include <pvm3.h>
    extern void pvm_init(int argc, char *argv[]);
    extern double armci_timer();
#   ifdef CRAY
#       define MPGROUP              (char *)NULL
#       define MP_INIT(argc,argv)
#   else
#       define MPGROUP              "mp_working_group"
#       define MP_INIT(argc,argv)   pvm_init(argc, argv)
#   endif
#   define MP_FINALIZE()        pvm_exit()
#   define MP_BARRIER()         pvm_barrier(MPGROUP,-1)
#   define MP_MYID(pid)         *(pid)   = pvm_getinst(MPGROUP,pvm_mytid())
#   define MP_PROCS(pproc)      *(pproc) = (int)pvm_gsize(MPGROUP)
#   define MP_TIMER             armci_timer
#   define MP_ASSERT(code)      code
#elif defined(MSG_COMMS_TCGMSG) || defined(MSG_COMMS_TCGMSG5) || defined(MSG_COMMS_TCGMSGMPI)
#   include <tcgmsg.h>
#   define MP_BARRIER()         tcg_synch(30000)
#   define MP_INIT(argc,argv)   tcg_pbegin((argc),(argv))
#   define MP_FINALIZE()        tcg_pend()
#   define MP_MYID(pid)         *(pid)   = (int)tcg_nodeid()
#   define MP_PROCS(pproc)      *(pproc) = (int)tcg_nnodes()
#   define MP_TIMER             tcg_time
#   define MP_ASSERT(code)      code
#elif defined(BGML)
    extern double armci_timer();
#   define MP_BARRIER()         armci_msg_barrier()
#   define MP_FINALIZE()     
#   define MP_INIT(argc,argv) 
#   define MP_MYID(pid)         *(pid)=armci_msg_me()
#   define MP_PROCS(pproc)      *(pproc)=armci_msg_nproc()
#   define MP_TIMER             armci_timer
#   define MP_ASSERT(code)      code
#else
#   include <mpi.h>
#   define MP_BARRIER()         MPI_Barrier(MPI_COMM_WORLD)
#   define MP_FINALIZE()        MPI_Finalize()
#   if defined(DCMF) || defined(MPI_MT)
    static inline int MPI_INIT_THREAD(int *argc, char ***argv) {
        int status;
        int provided;
        status = MPI_Init_thread(argc, argv, MPI_THREAD_MULTIPLE, &provided);
        return status;
    }
#       define MP_INIT(argc,argv)   MPI_INIT_THREAD(&(argc),&(argv))
#   else
#       define MP_INIT(argc,argv)   MPI_Init(&(argc),&(argv))
#   endif
#   define MP_MYID(pid)         MPI_Comm_rank(MPI_COMM_WORLD, (pid))
#   define MP_PROCS(pproc)      MPI_Comm_size(MPI_COMM_WORLD, (pproc))
#   define MP_TIMER             MPI_Wtime
#   define MP_ASSERT(code) do { \
        if (MPI_SUCCESS != (code)) { \
            MPI_Abort(MPI_COMM_WORLD, (code)); \
        } \
    } while (0)
#endif
#ifdef MPI_SPAWN 
#   define GA_INIT(argc,argv) GA_Initialize_args(&(argc),&(argv)) 
#   define ARMCI_INIT(argc,argv) ARMCI_Init_args(&(argc),&(argv)) 
#else 
#   define GA_INIT(argc,argv) GA_Initialize() 
#   define ARMCI_INIT(argc,argv) ARMCI_Init() 
#endif 

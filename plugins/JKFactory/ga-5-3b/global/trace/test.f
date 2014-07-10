       real a(10000)
       integer i
       call trace_init(1000)
       do k = 1,10
         call trace_stime()
         do i = 1,10000
            a(i) = sin(real(i+k))
         enddo
         call trace_etime()
         call trace_genrec(k,k,i,k,i,999)
       enddo
       call trace_end(99)
       end

       program async
       integer nodeid, nnodes
       integer me, proc
       integer req, ack, len,node, to
c
       call pbeginf
       proc = nnodes()
       me = nodeid()
       print *,'proc=',proc,' node=',me
       if(proc.lt.2) then
         call pend()
         return 
       endif
c
       call setdbg(1)
       if(me.eq.0)then
         print *, 'Checking non-blocking receive'
         call rcv(33, req,1, len, 1, node, 1) !blocking
         print *, 'received request ',me
         call snd(34, ack,1, 1, 1)
       endif
c
       if(me.eq.1)then
         call rcv(34, ack, 1, len, 0, node, 0)!nonblocking
         print *,'after nonblocking receive =',me
         call snd(33, req,1, 0, 0)
         print *,'after nonblocking send =',me
         call waitcom(0)
         print *, 'OK'
       endif
       call pend()
       end
         
       

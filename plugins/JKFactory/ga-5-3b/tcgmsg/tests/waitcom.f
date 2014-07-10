      program main
      implicit none
      integer maxloop
      parameter (maxloop = 67 )
      integer buf,me,nproc,loop,lenmes,node, received
      integer nnodes, nodeid 
c
      call pbeginf
c
      nproc = nnodes()
      me = nodeid()
      if(nproc.lt.3)then
          print *,'min 3 processes required ',nproc
          call parerr(0)
      endif
      received =0
      do loop = 1, maxloop
          node = Mod(loop,2)+1
          if(me.eq.0) then
            call snd(loop, buf, 1, node, 0)
          endif
          if(me.eq.node) then
             received = received +1
             call rcv(loop, buf, 1, lenmes, 0, node, 1)
          endif
      enddo
      if(me.eq.0)print *,'0: waiting for coms to node 1 to complete'
      call waitcom(1)
      if(me.eq.0)print *,'0: waiting for remaining coms to complete'
      call waitcom(-1)

      if(me.eq.0) then
          print *,'node=',me, maxloop,' messages sent asynchronously'
      else
          print *, 'node=',me,  received,' messages received'
      endif

c
      call pend
      end

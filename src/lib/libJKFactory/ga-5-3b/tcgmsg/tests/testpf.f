c     $Header: /tmp/hpctools/ga/tcgmsg/ipcv4.0/testpf.f,v 1.3 1995-02-24 02:18:01 d3h325 Exp $
      character*60 fname
      call pbeginf
      fname = ' '
      write(fname,'(a,i3.3)') '/tmp/pfcopy.test',nodeid()
      call pfcopy(5, 0, fname)
      call pend
      end

        program matrix_gen
        integer*4 i, k
        double precision amat(8,8), xvec(8), bvec(8)
        integer *4 irow(8), icol(8*8)
        integer*4 timeArray(3), ei,sf;
        open (unit = 8, file = "matrix.bin")
        ei = 8
        sf = 64
        do i = 1, 8
          do k = 1, 8
            do n=1,100000
              call itime(timeArray)
              amat(i,k) = amat(i,k) + 1212.12*
     &              rand(timeArray(1)/(timeArray(2)+timeArray(3)))
            end do
            call itime(timeArray)
            amat(i,k) = amat(i,k)+rand ( timeArray(1)+timeArray(2)+
     &                      timeArray(3))
          end do
          call itime(timeArray)
          xvec(i) = rand ( timeArray(1)+timeArray(2)+timeArray(3) )
        end do
        do i = 1,8
          bvec(i)=0D0
          irow(i)=(i-1)*8+1
          do k =1, 8
            bvec(i) = bvec(i)+xvec(k)*amat(i,k)
            icol((i-1)*8+k)=i
          end do
        end do
        write(8,*) ei
        write(8,*) sf
        write(8,*) amat 
        write(8,*) irow
        write(8,*) icol
        write(8,*) xvec
        close(8)
        write(6,*) 'amat =', amat 
        write(6,*) 'xvec =', xvec
        write(6,*) 'bvec=', bvec
        write(6,*) 'irow=', irow
        write(6,*) 'icol=', icol
        stop
        end

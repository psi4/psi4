*
* $Id: fci_hsort.F,v 1.2 1999/07/28 00:23:40 d3e129 Exp $
*
      subroutine fci_hsort(A, ionel, n)
      implicit none
      integer n
      integer A(n), temp, ionel(5,*)
      integer i, hsize
c     
c     RJH modified to sort on ionel(4,A(*))
c
c    build heap
c
      do i = n/2, 1, -1
         call fci_heapify(A, ionel, n, i)
      enddo
c     
c    main part of sort algorithm
c
      hsize = n
      do i    = n, 2, -1
c     
c    swap element i and 1
c
         temp = A(1)
         A(1) = A(i)
         A(i) = temp
c     
c    maintain heap property from element 1 down
c
         hsize= hsize - 1
         call fci_heapify(A, ionel, hsize, 1)
c     
      enddo
      end 
      
c     
c     establish heap property for a tree branch rooted at elem 
c     
      subroutine fci_heapify(A, ionel, n, elem)
      implicit none
      integer n, elem 
      integer A(n), ionel(5,*)
      integer left, right, largest, i, temp
c     
      i = elem
c     
c     Main Loop
c     
 100  continue 
      left  = 2*i
      right = 2*i + 1
      if (left. gt. n .and. right .gt. n) return !we traversed entire branch
c     
c     check heap property among element i and its children
c     
      if (left .le. n .and. ionel(4,A(left)) .gt. ionel(4,A(i))) then
         largest = left
      else
         largest = i
      endif
      if (right .le. n .and. 
     $     ionel(4,A(right)) .gt. ionel(4,A(largest))) largest = right
c     
      if (largest .ne. i) then
c     
c     swap array elements if largest is not i
c     
         temp = A(i)
         A(i) = A(largest)
         A(largest) = temp
c     
c     traverse down the tree
c     
         i = largest 
         goto 100
      endif
      end

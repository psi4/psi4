#if HAVE_CONFIG_H
#   include "config.fh"
#endif
      subroutine clear_hash
#include "common.fh"
      integer i
c
c   Clear all information from hash tables prior to setting new tables
c
      do i = 1, MAXAT
        link_a(i) = 0
        top_a(i) = 0
        hash_key_a(i) = 0
        hash_value_a(i) = 0
      end do
      a_cnt = 0
      return
      end
c
      subroutine add_hash_a(igx, idx)
#include "common.fh"
      integer igx, idx, ifunc, itmp
c
c   Add an element to hash table A, which tracks local particles.
c   Start by computing hash function value of idx.
c
      ifunc = mod(igx,MAXAT) + 1
c
c   Store value in a linked list
c
      a_cnt = a_cnt+1
      itmp = top_a(ifunc)
      top_a(ifunc) = a_cnt
      link_a(a_cnt) = itmp
      hash_key_a(a_cnt) = igx
      hash_value_a(a_cnt) = idx
      return
      end
c
      integer function get_hash_a(igx)
#include "common.fh"
      integer igx, ifunc, itmp, jtmp
c
c   Return the local index of a locally held particle from the global
c   index. Start by computing the hash function value of idx.
c
      ifunc = mod(igx,MAXAT) + 1
      get_hash_a = 0
      itmp = top_a(ifunc)
      do while (hash_key_a(itmp).ne.igx.and.itmp.gt.0)
        itmp = link_a(itmp)
      end do
      if (itmp.gt.0) then
        get_hash_a = hash_value_a(itmp)
      endif
      return
      end


      subroutine dbdqr(n, D, E, c1, c2)
      implicit none
      integer n
      double precision D(*),E(*),c1,c2
      
      integer i
      double precision cs,sn,r

      if (n.lt.2) return
      do i=1,n-1
         call dlartg(d(i),e(i),cs,sn,r)
         d(i) = r
         e(i) = sn*d(i+1)
         d(i+1) = cs*d(i+1)
      enddo
      call dlartg(d(n),e(n),cs,sn,r)
      d(n) = r
      e(n) = 0.0
      c1 = sn
      c2 = cs
      end


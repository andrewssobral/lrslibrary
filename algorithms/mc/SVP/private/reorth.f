      subroutine reorth(n,k,V,ldv,vnew,normv,index,alpha,work,
     c     iflag,nre)
c
c FORTRAN 77 version of MATLAB routine REORTH:
c
c REORTH   Reorthogonalize a vector using iterated Gram-Schmidt
c
c   [R_NEW,NORMR_NEW,NRE] = reorth(Q,R,NORMR,INDEX,ALPHA,METHOD)
c   reorthogonalizes R against the subset of columns of Q given by INDEX. 
c   If INDEX==[] then R is reorthogonalized all columns of Q.
c   If the result R_NEW has a small norm, i.e. if norm(R_NEW) < ALPHA*NORMR,
c   then a second reorthogonalization is performed. If the norm of R_NEW
c   is once more decreased by  more than a factor of ALPHA then R is 
c   numerically in span(Q(:,INDEX)) and a zero-vector is returned for R_NEW.
c
c   If method==0 then iterated modified Gram-Schmidt is used.
c   If method==1 then iterated classical Gram-Schmidt is used.
c
c   The default value for ALPHA is 0.5. 
c   NRE is the number of reorthogonalizations performed (1 or 2).

c References: 
c  Aake Bjorck, "Numerical Methods for Least Squares Problems",
c  SIAM, Philadelphia, 1996, pp. 68-69.
c
c  J.~W. Daniel, W.~B. Gragg, L. Kaufman and G.~W. Stewart, 
c  ``Reorthogonalization and Stable Algorithms Updating the
c  Gram-Schmidt QR Factorization'', Math. Comp.,  30 (1976), no.
c  136, pp. 772-795.
c
c  B. N. Parlett, ``The Symmetric Eigenvalue Problem'', 
c  Prentice-Hall, Englewood Cliffs, NJ, 1980. pp. 105-109

c  Rasmus Munk Larsen, DAIMI, 1998.
      implicit none
      integer n,k,ldv,i,iflag,nre
      double precision V(ldv,*),vnew(*),normv,index(*),work(*)
      double precision alpha,normv_old,dnrm2
      integer MAXTRY
      parameter (MAXTRY=4)
      external dgemv,dnrm2

c     Hack: If index .ne. 1:k we do MGS to avoid reshuffling.
      if (iflag.eq.1) then     
         do i=1,k
            if (int(index(i)).ne.i) then
               iflag=0
               goto 100
            endif
         enddo
      endif
 100  normv_old = 0
      nre = 0
      normv = dnrm2(n,vnew,1)
      do while ((normv.lt.alpha*normv_old .or. nre.eq.0))
         if (iflag.eq.1) then
c           CGS:
            call dgemv('T',n,k,1D0,V,ldv,vnew,1,0D0,work,1)
            call dgemv('N',n,k,-1D0,V,ldv,work,1,1D0,vnew,1)
         else
c           MGS:
            call MGS(n,k,V,ldv,vnew,index)
         endif
         normv_old = normv
         normv = dnrm2(n,vnew,1)         
         nre = nre + 1
         
         if  ( nre.gt.MAXTRY )  then
c
c     vnew is numerically in span(V) => return vnew = (0,0,...,0)^T
            normv = 0d0
            do i=1,n
               vnew(i) = 0d0
            enddo
            return
         endif
      enddo
      end
c
c****************************************************************************
c

      subroutine MGS(n,k,V,ldv,vnew,index)
      implicit none
      integer n,k,ldv
      double precision V(ldv,*),vnew(*),index(*)
      integer i,j,idx
      double precision s

c     
c     Modified Gram-Schmidt orthogonalization:
c     Orthogalizes vnew against the k vectors in V by the
c     iterative process     
c     
c     FOR i=1...k DO          
c       vnew = vnew - DOT( V(:,i), vnew ) * V(:,i) 
c

c     This simple version is faster on Pentium machines.
c     Compile with "g77 -O6 -funroll-all-loops -fomit-frame-pointer"
     
      do i=1,k
         idx = int(index(i))
         s = 0
         do j=1,n
            s = s + V(j,idx)*vnew(j)
         enddo
         do j=1,n
            vnew(j) = vnew(j) - s*V(j,idx)
         enddo
      enddo
      end
c
c****************************************************************************
c


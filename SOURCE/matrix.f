      subroutine matml(a,b,c,ndim)
c
c*** matrix multiplication c = a*b, with a, b and c square matrices
c
      implicit real*8 (a-h,o-z),integer(i-n)
      dimension a(ndim*ndim),b(ndim*ndim),c(ndim*ndim)
c
      c(1:ndim*ndim)=0.d0
      do i=1,ndim
        ii=(i-1)*ndim
        do j=1,ndim
          ij=ii+j
          do k=1,ndim
            ik=(i-1)*ndim+k
            kj=(k-1)*ndim+j
            c(ij)=c(ij)+a(ik)*b(kj)
          enddo
        enddo
      enddo
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine vecmat(vec1,xmat,ndim,vec2)
c
c**  multiplication of vector * matrix
c**    vec2(i) = vec1(j) * xmat(j,i) sum over j
c
      implicit real*8 (a-h,o-z),integer(i-n)
      dimension vec1(ndim),xmat(ndim,ndim),vec2(ndim)
      do i=1,ndim
        vec2(i) = 0.d0
        do j = 1,ndim
          vec2(i) = vec2(i)+vec1(j)*xmat(j,i)
        enddo
      enddo
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine sqrmat(a,n)
c
c** expand matrix in triangular form to full  n*n matrix.
c
c note: if the off-diagonal element in the triangular matrix is the
c       sum of the two off-diagonal elements (i,j) and (j,i) then
c       use routine fulmat (also see comment in routine fulmat).
c
      implicit real*8 (a-h,o-z),integer(i-n)
      dimension a(n*n)
      do i=n,1,-1
        do j=i,1,-1
          a(n*(i-1)+j)=a(i*(i-1)/2+j)
        enddo
      enddo
      do i=1,n
        do j=1,i-1
          a(n*(j-1)+i)=a(n*(i-1)+j)
        enddo
      enddo
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine diaglz(a,r,n)
c
c Diagonalise a matrix by the jacobi method
c
c    arguments
c          a  the matrix to be diagonalised
c                 this matrix is destroyed - the eigenvalues are
c                 put onto the main diagonal, zeroes elsewhere
c          r   the matrix of eigenvectors produced
c          n   the overall dimension of the matrix  a
c
      implicit real*8 (a-h,o-z),integer(i-n)
      parameter(small = 1.d-24)
c
      dimension a(n,n),r(n,n)
c
      root2 = sqrt(2.d0)
c
c** generate identity matrix
c
      r(1:n,1:n)=0.d0
      do i = 1,n
        r(i,i)=1.d0
      enddo
c
c** compute initial and final norms (anorm and anormx)
c
      anorm=0.d0
      do i=1,n-1
        do j=i+1,n
          anorm=anorm+abs(a(i,j))
        enddo
      enddo
      if (anorm.gt.0.d0) then
        anorm=root2*anorm
        anrmx=anorm*small/n
c
c** initialize indicators and compute threshold, thr 
c
        ind=0
        thr=anorm
        thr=thr/n
        l=1
        m=l+1
c
c** compute the sin and cos of the rotation which transforms
c** element a(l,m) to zero
c
  260   if(abs(a(l,m)).gt.thr) then
          ind=1
          x=(a(l,l)-a(m,m))/2.d0
          y=-a(l,m)/sqrt(a(l,m)*a(l,m)+x*x)
          if (x.lt.0.d0) y=-y
          sinx=y/sqrt(2.d0*(1.d0+(sqrt(1.d0-y*y))))
          sinx2=sinx*sinx
          cosx=sqrt(1.d0-sinx2)
          cosx2=cosx*cosx
          sincs =sinx*cosx
c
c** rotate l and m columns
c
          do i=1,n
            if ((i.ne.l).and.(i.ne.m)) then
              x=a(i,l)*cosx-a(i,m)*sinx
              a(i,m)=a(i,l)*sinx+a(i,m)*cosx
              a(m,i)=a(i,m)
              a(i,l)=x
              a(l,i)=a(i,l)
            endif
            x=r(i,l)*cosx-r(i,m)*sinx
            r(i,m)=r(i,l)*sinx+r(i,m)*cosx
            r(i,l)=x
          enddo
          x=2.d0*a(l,m)*sincs
          y=a(l,l)*cosx2+a(m,m)*sinx2-x
          x=a(l,l)*sinx2+a(m,m)*cosx2+x
          a(l,m)=(a(l,l)-a(m,m))*sincs+a(l,m)*(cosx2-sinx2)
          a(m,l)=a(l,m)
          a(l,l)=y
          a(m,m)=x
        endif
c
c** tests for completion; test for m = last column
c
        if (m.ne.n) then
          m=m+1
          goto 260
        endif
c
c** test for l = second from last column
c
        if (l.ne.(n-1)) then
          l=l+1
          m=l+1
          goto 260
        endif
c
        if (ind.eq.1) then
          ind=0
          l=1
          m=l+1
          goto 260
        endif
c
c** compare threshold with final norm
c
        if (thr.gt.anrmx) then
          thr=thr/n
          l=1
          m=l+1
          goto 260
        endif
      endif
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine eigsrt(v,d,ism,n)
c
c** sort eigenvectors v, eigenvalues d, and orbital symmetry ism  ***
c
      implicit real*8 (a-h,o-z),integer(i-n)
      dimension v(n,n),d(n),ism(n)
      dimension vtmp(n)
c
      do i=1,n-1
        k=i
        p=d(i)
        do j=i+1,n
          if (d(j).ge.p) then
            k=j
            p=d(j)
          endif
        enddo
        if (k.ne.i) then
          d(k)=d(i)
          d(i)=p
          itmp=ism(k)
          ism(k)=ism(i)
          ism(i)=itmp
          vtmp(1:n)=v(1:n,k)
          v(1:n,k)=v(1:n,i)
          v(1:n,i)=vtmp(1:n)
        endif
      enddo
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine tmtdag(a,na,b,nb,q,c)
c
c subroutine for matrix transformation : b = qaq+
c      b(kl) = q(ki)*a(ij)*q(lj)  sum over i,j
c
c a  : symmetric matrix with dimension na*na.
c      element a(i,j) stored in position i*(i-1)/2+j (i>=j)
c b  : symmetric matrix with dimension nb*nb.
c      element b(i,j) stored in position i*(i-1)/2+j (i>=j)
c q  : transformation matrix with dimension nb*na
c
c before transformation multiply off-diagonal elements of
c matrix a with constant c.  after transformation divide
c off-diagonal elements of matrix a and b by constant c.
c if a is a triangular matrix with element i*(i-1)/2+j equal to the
c sum of the off-diagonal elements (i,j) and (j,i) then use c=0.5.
c
c
      implicit real*8 (a-h,o-z),integer(i-n)
c
      dimension a(na*(na+1)/2),b(nb*(nb+1)/2),q(nb,na)
      dimension t(na)
c
      ci = 1.d0/c
      b(1:nb*(nb+1)/2)=0.d0
c
c** multiply off diagonal elements with factor c 
c
	open(99,file='tmtdag.txt')
        do j=1,nb
        write(99,'(i3,14(f11.6))') j,(a(j*(j-1)/2+k),
     +   k=1,j)
        enddo
	do j=1,nb
        write(99,'(i3,14(f11.6))') j,(q(j,k),
     +   k=1,na)
        enddo
        close(99)
      a(1:na*(na+1)/2)=a(1:na*(na+1)/2)*c
      do i=1,na
        ii=i*(i+1)/2
        a(ii)=a(ii)*ci
      enddo
      do k=1,nb
        do j=1,na
          t(j)=0.d0
          jj=j*(j-1)/2
          do i=1,j
            t(j)=t(j)+a(jj+i)*q(k,i)
c	if (j.eq.1) print*,k,j,i,a(jj+i),q(k,i),t(j)
          enddo
          do i=j+1,na
            t(j)=t(j)+a(i*(i-1)/2+j)*q(k,i)
c	if (j.eq.1) print*,i,j,k,
c     +  a(i*(i-1)/2+j),q(k,i),t(j)
          enddo
        enddo
        do l=1,k
          kl=k*(k-1)/2+l
          do j=1,na
            b(kl)=b(kl)+t(j)*q(l,j)
c	if (kl.eq.1) print*,k,l,j,t(j),q(l,j)
          enddo
        enddo
      enddo
c
      b(1:nb*(nb+1)/2)=b(1:nb*(nb+1)/2)*ci
      do i=1,nb
        ii=i*(i+1)/2
        b(ii)=b(ii)*c
      enddo
      a(1:na*(na+1)/2)=a(1:na*(na+1)/2)*ci
      do i=1,na
        ii=i*(i+1)/2
        a(ii)=a(ii)*c
      enddo
c
      return
      end
c
c-----------------------------------------------------------------------
c
      function valmat(xmat,norb,val,qfull)
c
c valmat = xmat(i,j) * val(i) * val(j)
c
c if qfull=true then a full norb*norb matrix is used, else a matrix
c in triangular form
c
      implicit real*8 (a-h,o-z),integer(i-n)
      logical qfull
      dimension xmat(norb*norb),val(norb)
      dimension a(norb)
      px = 0.0d0
      if (qfull) then
        do 702 j = 1,norb
          ijbas = (j-1)*norb
          do 705 i = 1,norb
            a(i) = xmat(ijbas+i)*val(j)
 705      continue
          do 706 i = 1,norb
            px = px+a(i)*val(i)
 706      continue
 702    continue
        valmat = px
      else
        do 708 i = 1,norb
          ijbas = i*(i-1)/2
          b = 0.0d0
          do 709 j = 1,i
            b = b+xmat(ijbas+j)*val(j)
 709      continue
          px = px+b*val(i)
 708    continue
        valmat = px
      end if
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine gtriad(n,scalar,r,a,b)
c..    r = a + scalar*b
      implicit double precision (a-h,o-z), integer (i-n)
c
      dimension a(*),r(*),b(*)
      do 1 i=1,n
 1    r(i)=a(i)+scalar*b(i)
      return
c
      end

c
      subroutine minvr(a,tol,det,ier,n)
c
c in-place inversion of a square matrix by means of Gauss-Jordan 
c elimination with partial pivoting
c
c input  : n   - order of the (square) matrix a, to be inverted.
c          tol - criterium for singularity: if any pivot has
c                absolute value less than tol, the matrix is
c                regarded as singular, and no inverse is calculated.
c
c output : det - determinant of input matrix.
c          ier - error parameter. 0: inverse computed, 1: singular
c
c in/out : a   - in  : matrix to be inverted
c                out : inverse (if a-input was not singular)
c
      implicit real*8 (a-h,o-z),integer(i-n)
      dimension a(n,n)
      dimension ind(n),w(n)
c
c** initialization
c
      det = 1.d0
      ier = 0
      do i=1,n
        ind(i)=i
      enddo
c
c** loop over rows
c
      do j=1,n
c
c**determination of maximum element in remaining columns
c
        c=tol
        k=0
        do i=j,n
          if (abs(a(j,i)).gt.c) then
            k=i
            c=abs(a(j,i))
          endif
        enddo
        if (k.eq.0) then
          ier=1
          return
        endif
c
c** calculate determinant and pivot
c
        det = det*a(j,k)
        if (j.ne.k) det = -det
        pivot = 1.d0/a(j,k)
c
c** update index array for final row order restoring 
c
        jtemp = ind(j)
        ind(j) = ind(k)
        ind(k) = jtemp
c
c** multiply column k with pivot and store in w
c
        call scaler(n,pivot,w,a(1,k))
c
c** move column j to k
c
        call fmove(a(1,j),a(1,k),n)
c
c** elimination of row j
c
        do k=1,n
          rowj=-a(j,k)
          call gtriad(n,rowj,a(1,k),a(1,k),w)
          a(j,k)=pivot*rowj
        enddo
c
c** restore column j 
c
        call fmove(w,a(1,j),n)
        a(j,j)=pivot
      enddo
c
c** restore row order of last column and store temporarely in w
c
      do j=1,n
        w(ind(j))=a(j,n)
      enddo
c
c** restore row order and move all columns one to the right
c
      do j=n,2,-1
        do i=1,n
          a(ind(i),j)=a(i,j-1)
        enddo
      enddo
c
c** copy and shift columns
c
      do j=1,n-1
        call fmove(a(1,j+1),a(1,j),n)
      enddo
c
c** restore last column
c
      call fmove(w,a(1,n),n)
c
      return
      end

c** generic fmove
c  
      subroutine fmove(a,b,n)
      implicit real*8  (a-h,p-w),integer(i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension a(*),b(*)
      do 3 i=1,n
 3    b(i)=a(i)
      return
      end

      
      subroutine scaler(n,scalar,r,b)
c...   ***note***  more general then dscal
      implicit double precision (a-h,o-z), integer (i-n)
      dimension r(n),b(n)
c     
      do 1 i=1,n
 1    r(i)=scalar*b(i)
c     
      return         
      end            



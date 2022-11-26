c*******************************************************
      program test_bezier
c*******************************************************
      implicit double precision (a-h,o-z)
      dimension val(1000),x(1000),y(1000)
      dimension xb(25),yb(25)
      common /bez1/ bezm(25,25)
c
      open(1,file='solution.txt')
	  do i=1,26
	    read(1,*,end=10)xb(i),yb(i)
	  enddo
	  stop  'Increase dimension (xb.yb)'
  10  nb=i-1
      close(1)
c
      call inibezier(nb)
c
      write(*,*)' Type the number of points along the fival curve'
	  read(*,*)ideg
      call usebezier(xb,yb,nb,x,y,ideg-1)
c
      open(1,file='curve.txt')
	  do i=1,ideg
	  write(1,'(2(3x,f12.7))')x(i),y(i)
	  enddo
      close(1)
c
      stop
      end
c
c
c *********************************************************************
      subroutine inibezier(nco)
c *********************************************************************
c     NCO control points
c
      implicit double precision (a-h,o-z)
      common /bez1/ bezm(25,25)
c
      do 1 mi=0,nco-1      !  for the control points
        b=0.d0
        c=0.d0
          do 1 i=0,nco-1
            call paragon (nco-1,i ,kres1)
            call paragon (i    ,mi,kres2)
            kres3=(-1)**(i-mi)
            coeffi = dfloat(kres1*kres2*kres3)
            if(mi.gt.i) coeffi=0.d0
            bezm(mi+1,i+1) = coeffi
  1   continue
c
      return
      end
c
c
c *********************************************************************
      subroutine paragon (n,i,k)  ! n=UP, i=LOW, k=RESULT
c *********************************************************************
      implicit double precision (a-h,o-z)
c
      ks=max(i,n-i)+1
      kp=min(i,n-i)
      k=1
      do iii=ks,n
        k=k*iii
      enddo
      do iii=1,kp
        k=k/iii
      enddo
      return
      end
c
c
c *********************************************************************
      subroutine usebezier   (xco,yco,nco,x,y,ideg)
c *********************************************************************
c     NCO control points (xco,yco)
c     IDEG+1 final points with coordinates (X,Y)
c
      implicit double precision (a-h,o-z)
      common /bez1/ bezm(25,25)
      dimension x(1),y(1),xco(1),yco(1)
c
      aa1  =  0.5d0
      dd   =  0.1
	  do k=1,ideg+1
	  y(k)=dfloat(k-1)/dfloat(ideg)
	  enddo
c
      do 10 kpoi=0,ideg    !  for the outcome
      kpoi1 = kpoi+1
      tlocal = y(kpoi1)    !  current parametric value
      x  (kpoi1) = 0.d0
      y  (kpoi1) = 0.d0
        do mi=0,nco-1      !  for the control points
        b=0.d0
          do  i=0,nco-1
            b = b + bezm(mi+1,i+1) * tlocal**i
          enddo  ! i
        x(kpoi1)   = x  (kpoi1) + b*xco(mi+1)
        y(kpoi1)   = y  (kpoi1) + b*yco(mi+1)
      enddo  !  mi
  10  continue
c
      return
      end
c
c

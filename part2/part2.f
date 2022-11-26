C*******************************************************
      PROGRAM AMAKSI2
C*******************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (KDIM=30)
      PARAMETER (LDIM=15)
      DIMENSION X(KDIM),Y(KDIM), T(KDIM), T0(KDIM), FR(KDIM),aux(kdim),
     &sol(kdim), E(KDIM)
      DIMENSION XB(LDIM),YB(LDIM), B(LDIM)
      DIMENSION A(LDIM,LDIM)
      COMMON /BEZ1/ NB/BEZ2/ BEZM(LDIM,LDIM)

C
      NB=10
      CALL INIBEZIER(NB)
      
      OPEN(1,FILE='C.txt')     !!!! MONO GIA MENA
	  DO I=1,NB
          WRITE(1,*)(BEZM(I,J),j=1,nb)
	  ENDDO
       CLOSE(1)

      OPEN(1,FILE='syntetagmenes.txt')
      DO I=1,KDIM
      READ(1,*,END=33)X(I),Y(I)
      ENDDO
 33   CLOSE(1)
      N=I-1
C
      XB(1)=0.D0
      XB(NB)=40.D0
      XB(2)=4.D0
      XB(3)=6.D0
      XB(4)=10.D0
      XB(5)=17.D0
      XB(6)=20.D0
      XB(7)=23.D0
      XB(8)=29.D0
      XB(9)=37.D0
C
      T0(1)=0.D0
      DO I=2,N
      T0(I)=T0(I-1)+1.D0/(N-1)
      ENDDO

      ER=0.0001D0
      NMAX=10
      DO I=1,N
      CALL NEWTON (X(I),T0(I),er,nmax,T(I),FR(I))
      ENDDO


C      DO I=1,N
C      T0(I)=X(I)/(XB(2)*BEZM(2,2))
C      ENDDO


      DO I=1,(NB)
         DO J=1,(NB)
            A(I,J)=0.d0
            IF (I.NE.1) THEN
            IF (J.NE.1) THEN
            IF (I.NE.NB) THEN
            IF (J.NE.NB) THEN

           DO K=1,N
               C1=0.d0
               C2=0.d0
                DO M=1,(NB)
                C1=C1+BEZM(I,M)*(T(K)**(M-1))
                C2=C2+BEZM(J,M)*(T(K)**(M-1))
                ENDDO
                C=C1*C2
              A(I,J)=A(I,J)+C
            ENDDO
            ENDIF
            ENDIF
            ENDIF
            ENDIF
          ENDDO
      ENDDO

      DO I=1,(NB)
      B(I)=0.d0
      IF (I.NE.1) THEN
      IF (I.NE.NB) THEN
            DO K=1,N
            C=0.d0
               DO M=1,(NB)
               C=C+BEZM(I,M)*(T(K)**(M-1))
               ENDDO
             B(I)=B(I)+C*Y(K)
        ENDDO
        ENDIF
        ENDIF
      ENDDO

      do i=1,(nb-1)
       do j=1,(nb-1)
       A(I,J)=A(I+1,J+1)
       ENDDO
       B(I)=B(I+1)
      ENDDO
       
      CALL lu((NB-2),A,B,SOL)
      
      DO I=(NB-1),2,-1
      SOL(I)=SOL(I-1)
      ENDDO
      SOL(1)=0
      SOL(NB)=0


      open(1,file='solution.txt')
      do i=1,NB
      write(1,*)XB(I), sol(i)
      enddo
      close(1)

c
      STOP
      END



C *********************************************************************
      SUBROUTINE INIBEZIER(NCO)
C *********************************************************************
C     NCO CONTROL POINTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (LDIM=15)
      COMMON /BEZ1/ NB/BEZ2/ BEZM(LDIM,LDIM)

C
      DO 1 MI=0,NCO-1      !  FOR THE CONTROL POINTS
        B=0.D0
        C=0.D0
          DO 1 I=0,NCO-1
            CALL PARAGON (NCO-1,I ,KRES1)
            CALL PARAGON (I    ,MI,KRES2)
            KRES3=(-1)**(I-MI)
            COEFFI = DFLOAT(KRES1*KRES2*KRES3)
            IF(MI.GT.I) COEFFI=0.D0
            BEZM(MI+1,I+1) = COEFFI
  1   CONTINUE
C
      RETURN
      END
C
C
C *********************************************************************
      SUBROUTINE PARAGON (N,I,K)  ! N=UP, I=LOW, K=RESULT
C *********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      KS=MAX(I,N-I)+1
      KP=MIN(I,N-I)
      K=1
      DO III=KS,N
        K=K*III
      ENDDO
      DO III=1,KP
        K=K/III
      ENDDO
      RETURN
      END

C***********************************************************************
      subroutine lu(N,A,B,SOL)
C***********************************************************************
      implicit double precision (a-h,o-z)
      PARAMETER (LDIM=15)
      dimension A(ldim,ldim),b(ldim),aux(ldim),sol(ldim)


c     L and U matrices will be written over A
c     Crout : diagonal(U)=UNIT
      call crout(ldim,n,a)
c
      open(2,file='res_lu.txt')
        write(2,*)' L-matrix: '
            do k1=1,n
            write(2,'(6(1x,f10.5))')(a(k1,k2),k2=1,k1)
            enddo
        write(2,*)' U-matrix: (UNIT diagonal omitted)'
            do k1=1,n
            write(2,'(6(1x,f10.5))')(a(k1,k2),k2=k1+1,n)
            enddo
      close(2)
c
c     Solution - Step (1) :  L.AUX=b  -->  find AUX
      aux(1)=b(1)/a(1,1)
      do i=2,n
        sum=0.d0
        do j=1,i-1
        sum=sum+a(i,j)*aux(j)
        enddo
      aux(i)=(b(i)-sum)/a(i,i)
      enddo
c
c     Solution - Step (2) :  U.SOL=AUX  -->  find SOL
      sol(n)=aux(n)
      do i=n-1,1,-1
        sum=0.d0
        do j=n,i+1,-1
        sum=sum+a(i,j)*sol(j)
        enddo
      sol(i)=aux(i)-sum
      enddo
c
      RETURN
      END
C *********************************************************************
      subroutine crout(kdim,n,a)
C *********************************************************************
      implicit double precision (a-h,o-z)
      dimension a(kdim,kdim)
      eps=1.d-8
c
      do j=2,n
       a(1,j) = a(1,j) / a(1,1)
      enddo
c
      do k=2,n
        do i=k,n
         sum=0.d0
         do j=1,k-1
           sum=sum+a(i,j)*a(j,k)
         enddo
         a(i,k)=a(i,k)-sum
        enddo
c
        do j=k+1,n
         sum=0.d0
         do i=1,k-1
           sum=sum+a(k,i)*a(i,j)
         enddo
         a(k,j)=(a(k,j)-sum)/a(k,k)
        enddo
      enddo
c
      return
      end


c**************************************************
      FUNCTION FF(T)
c**************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (LDIM=15)
      DIMENSION XB(LDIM)
      COMMON /BEZ1/ NB/BEZ2/ BEZM(LDIM,LDIM)

      XB(1)=0.D0
      XB(NB)=40.D0
      XB(2)=4.D0
      XB(3)=6.D0
      XB(4)=10.D0
      XB(5)=17.D0
      XB(6)=20.D0
      XB(7)=23.D0
      XB(8)=29.D0
      XB(9)=37.D0

      G=0.d0
      DO I=1,NB
      C=0.d0
      DO J=1,NB
      C=C+BEZM(I,J)*T**(J-1)*XB(I)
      ENDDO
      G=G+C
      ENDDO
      FF=G
      RETURN
      END


c**************************************************
      FUNCTION FD1(X)
c**************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (LDIM=15)
      DIMENSION XB(LDIM)
      COMMON /BEZ1/ NB/BEZ2/ BEZM(LDIM,LDIM)
      XB(1)=0.D0
      XB(NB)=40.D0
      XB(2)=4.D0
      XB(3)=6.D0
      XB(4)=10.D0
      XB(5)=17.D0
      XB(6)=20.D0
      XB(7)=23.D0
      XB(8)=29.D0
      XB(9)=37.D0

      F=0.d0
      DO I=1,NB
      C=0.d0
      DO J=2,NB
      C=C+((J-1)*BEZM(I,J)*X**(J-2))*XB(I)
      ENDDO
      F=F+C
      ENDDO

      FD1=F
      RETURN
      END


c**************************************************
      SUBROUTINE NEWTON (C,x0,er,nmax,xr,FR)
c**************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (LDIM=15)
      COMMON /BEZ1/ NB/BEZ2/ BEZM(LDIM,LDIM)

c
c-- initial values
      niter=0
      xold=x0
      fr=FF(xold)-C
c
c___ start of iterations _______________
c
      write(*,*) '  niter   xr   ea   fr '
      DO
         niter=niter+1
         fgr=FD1(xold)
         IF(fgr .NE. 0.) THEN
            xr=xold-fr/fgr
         ELSE
            WRITE(*,*) 'error!  zero derivative'
            RETURN
         ENDIF
c-- xr is a roote?
         fr=FF(xr)-C
         if(fr. EQ. 0.) THEN
            ea=0.
            EXIT
         ENDIF
c-- relative (or absolute) error
         IF (xr .NE. 0.) THEN
            ea=ABS((xr-xold)/xr)
         ELSE
            ea=ABS(xr-xold)
         ENDIF
c-- monitoring (optional)
         WRITE(*,901) niter,xr,ea,fr
c-- exit checks
         IF (ea .LE. er)  EXIT
         IF (niter .GE. nmax)  THEN
            WRITE(*,*) 'warning!  iterations limit'
            EXIT
         ENDIF
         xold=xr
      END DO

      RETURN
 901  FORMAT(2X,I4,E20.10,2E20.6)
      END




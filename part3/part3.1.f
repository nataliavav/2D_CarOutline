      PROGRAM KINHSH1
C*******************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (LDIM=15)
      REAL *8 T
      DIMENSION YB(LDIM), B(LDIM)
      DIMENSION aux(2),ak(2,4)
      DIMENSION A(LDIM,LDIM)
      COMMON /BEZ1/ NB/BEZ2/ BEZM(LDIM,LDIM)

      OPEN(1,FILE='yb.txt')
      DO I=1,LDIM
      READ(1,*,END=33)YB(I)
      ENDDO
 33   CLOSE(1)
      NB=I-1
      CALL INIBEZIER(NB)

C
      x0=0.7d0
      ER=0.0001D0
      NMAX=10
      CALL NEWTON (YB,x0,er,nmax,T,FR)
      
      
      Y=0.d0
      DO I=1,NB
      C=0.d0
      DO J=1,NB
      C=C+BEZM(I,J)*T**(J-1)*YB(I)
      ENDDO
      Y=Y+C
      ENDDO
      
      OPEN(1,FILE='Ty.txt')     !
        WRITE(1,*) y, ’
      CLOSE(1)

      stop
      end


c**************************************************
      SUBROUTINE FF(YB,X,F)
c**************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (LDIM=15)
      DIMENSION YB(LDIM)
      COMMON /BEZ1/ NB/BEZ2/ BEZM(LDIM,LDIM)

      F=0.d0
      DO I=1,NB
      C=0.d0
      DO J=2,NB
      C=C+((J-1)*BEZM(I,J)*X**(J-2))*YB(I)
      ENDDO
      F=F+C
      ENDDO

      RETURN
      END


c**************************************************
      SUBROUTINE FD1(YB,X,F)
c**************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (LDIM=15)
      DIMENSION YB(LDIM)
      COMMON /BEZ1/ NB/BEZ2/ BEZM(LDIM,LDIM)
      
      F=0.d0
      DO I=1,NB
      A=0.d0
      DO J=3,NB
      A=A+((J-1)*(J-2)*BEZM(I,J)*X**(J-3))*YB(I)
      ENDDO
      F=F+A
      ENDDO
      
      RETURN
      END


c**************************************************
      SUBROUTINE NEWTON (YB,x0,er,nmax,xr,FR)
c**************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (LDIM=15)
      DIMENSION YB(LDIM)
      COMMON /BEZ1/ NB/BEZ2/ BEZM(LDIM,LDIM)

c
c-- initial values
      niter=0
      xold=x0
      CALL FF(YB,xold,fr)
c
c___ start of iterations _______________
c
      write(*,*) '  niter   xr   ea   fr '
      DO
         niter=niter+1
         CALL Fd1(YB,xold,fgr)
         IF(fgr .NE. 0.) THEN
            xr=xold-fr/fgr
         ELSE
            WRITE(*,*) 'error!  zero derivative'
            RETURN
         ENDIF
c-- xr is a roote?
       CALL FF(YB,xr,fr)
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

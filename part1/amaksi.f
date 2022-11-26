	PROGRAM CAR
	IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	PARAMETER (KDIM=30)
	PARAMETER (LDIM=60)
	DIMENSION BEF(KDIM),AFT(KDIM),DIAG(KDIM),X(KDIM),Y(KDIM),M(KDIM)
        DIMENSION XNEW(LDIM), YNEW(LDIM), RHSOL(KDIM)
	CHARACTER F
	REAL *8 M
C
	OPEN(1,FILE='syntetagmenes.txt')
	DO I=1,KDIM
	READ(1,*,END=33)X(I),Y(I)
	ENDDO
 33     CLOSE(1)
	N=I-1
C
        DO J=1,2
 	RHSOL(1)=0.d0
	RHSOL(N)=0.d0
	DIAG(1)=1.d0
	DIAG(N)=1.d0
	AFT(1)=0.d0
	BEF(N)=0.d0
C
        IF (J.EQ.1) THEN
	DO I=2,(N-1)
	DIAG(I)=4.d0
	BEF(I)=1.d0
	AFT(I)=1.d0
	RHSOL(I)=6*(X(I-1)-2*X(I)+X(I+1))
	ENDDO
C
	CALL trdiag (kdim,n,bef,diag,aft,rhsol)
	DO I=1,N
	M(I)=RHSOL(I)
	ENDDO
	ENDIF
C
        IF (J.EQ.2) THEN
        DO I=2,(N-1)
	DIAG(I)=4.d0
	BEF(I)=1.d0
	AFT(I)=1.d0
	RHSOL(I)=6*(Y(I-1)-2*Y(I)+Y(I+1))
	ENDDO
	CALL trdiag (kdim,n,bef,diag,aft,rhsol)
	ENDIF
        ENDDO
C
	J=1
	DO I=1,N
	XNEW(J)=X(I)
	YNEW(J)=Y(I)
	J=J+2
	ENDDO
C
	J=1
	DO I=2,(2*N),2
	XNEW(I)=XNEW(I-1)+((XNEW(I+1)-XNEW(I-1))-(1.d0/6.d0)*M(J+1)-
     &  (1.d0/3.d0)*M(J))*(1.d0/2.d0)+(1.d0/2.d0)**3*M(J)+
     &  (1.d0/6.d0)*(M(J+1)-M(J))*(1.d0/2.d0)**3
	YNEW(I)=YNEW(I-1)+((YNEW(I+1)-YNEW(I-1))-(1.d0/6.d0)*RHSOL(J+1)-
     &  (1.d0/3.d0)*RHSOL(J))*(1.d0/2.d0)+(1.d0/2.d0)**3*RHSOL(J)+
     &  (1.d0/6.d0)*(RHSOL(J+1)-RHSOL(J))*(1.d0/2.d0)**3
C
	J=J+1
	ENDDO
C
	OPEN(2,FILE='news.txt')
	DO I=1,(2*N-1)
	WRITE(2,*)XNEW(I),YNEW(I)
	ENDDO
	CLOSE (2)
C
	STOP
	END

C###############################################################
      subroutine trdiag (kdim,n,bef,diag,aft,rhsol)
      implicit double precision (a-h,o-z)
      dimension bef(kdim),aft(kdim),rhsol(kdim),diag(kdim)
C
      eps=1.d-8
      do 10 i=2,n
      if(dabs(diag(i-1)).lt.eps) stop  'Divide by zero'
      r = bef(i)/diag(i-1)
      diag(i)  = diag(i)  - r* aft(i-1)
      rhsol(i) = rhsol(i) - r* rhsol(i-1)
  10  continue
C
      rhsol(n) = rhsol(n)/diag(n)
      do 20 i=n-1,1,-1
      rhsol(i) = (rhsol(i)-aft(i)*rhsol(i+1)) / diag(i)
  20  continue
C
      return
      end
C

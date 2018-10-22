	SUBROUTINE MARUPV(TVX,TX,GB,H,NI,NM)
	IMPLICIT NONE
	INTEGER NI,NM,I,J,K
	REAL*8 TVX(NI-1,NI,NM),TX(NI,NI,NM),GB(NI),H(NI)
C
	DO K=1,NM
	  DO J=1,NI
	    DO I=1,NI-1
	      TVX(I,J,K)=GB(I)*(TX(I,J,K)-TX(I+1,J,K))+H(I)*TVX(I,J,K)
	    END DO
	  END DO
	END DO
C
	RETURN
	END
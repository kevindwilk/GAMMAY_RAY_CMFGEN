	SUBROUTINE EDGEWR(EDGE,N,DESC,LU)
	IMPLICIT NONE
	INTEGER N,LU,I
	REAL*8 EDGE(N)
	CHARACTER*(*) DESC
C
	IF(EDGE(1) .NE. 0)THEN
	  WRITE(LU,*)' '
	  WRITE(LU,*)'Ionization frequencies for'//DESC
	  DO I=1,N
	    WRITE(LU,*)EDGE(I),I
	  END DO
	END IF
C
	RETURN
	END

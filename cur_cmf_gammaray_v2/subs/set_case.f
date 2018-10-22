C
C Altered 09-Aug-1991 : If either of STREND and STRBEG is zero, the full
C                       STRING is converted.

C Altered 20-Nov-1989 : If STREND is ZERO, STREND is determined by the
C                       string length.
C
	SUBROUTINE SET_CASE_LOW(STRING,STRBEG,STREND)
	IMPLICIT NONE
	CHARACTER*(*) STRING
	INTEGER VAL,J,STRBEG,STREND,END_STR,ST_STR
C
	IF(STREND .EQ. 0 .OR. STRBEG .EQ. 0)THEN
	  ST_STR=1
	  END_STR=LEN(STRING)
	ELSE
	  ST_STR=STRBEG
	  END_STR=STREND
	END IF
C
	DO J=ST_STR,END_STR
	  IF( STRING(J:J) .GE. 'A' .AND. STRING(J:J) .LE. 'Z' )THEN
	    VAL=ICHAR(STRING(J:J))
	    VAL=VAL+32
	    STRING(J:J)=CHAR(VAL)
	  END IF
	END DO
C
	RETURN
	END
C
	SUBROUTINE SET_CASE_UP(STRING,STRBEG,STREND)
	IMPLICIT NONE
	CHARACTER*(*) STRING
	INTEGER VAL,J,STRBEG,STREND,END_STR,ST_STR
C
	IF(STREND .EQ. 0 .OR. STRBEG .EQ. 0)THEN
	  ST_STR=1
	  END_STR=LEN(STRING)
	ELSE
	  ST_STR=STRBEG
	  END_STR=STREND
	END IF
C
	DO J=ST_STR,END_STR
	  IF( STRING(J:J) .GE. 'a' .AND. STRING(J:J) .LE. 'z' )THEN
	    VAL=ICHAR(STRING(J:J))
	    VAL=VAL-32
	    STRING(J:J)=CHAR(VAL)
	  END IF
	END DO
C
	RETURN
	END
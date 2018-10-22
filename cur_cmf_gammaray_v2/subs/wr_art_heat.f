!
! Routine to write out the artificial heating term.
!
	SUBROUTINE WR_ART_HEAT(ARTIFICIAL_HEAT_TERM,NETCR,TOTCR,COUNTER,ND,LU)
	IMPLICIT NONE
!
! Created 28-Sep-2004
!
	INTEGER COUNTER,ND,LU
	REAL*8 ARTIFICIAL_HEAT_TERM(ND)
	REAL*8 NETCR(ND)
	REAL*8 TOTCR(ND)
!
	REAL*8 T1
	INTEGER MS,MF,I
!
	MS=(COUNTER-1)*10+1	
	MF=COUNTER*10
	IF(MF .GT. ND)MF=ND
!
	T1=SUM(ARTIFICIAL_HEAT_TERM)
	IF(T1 .EQ. 0.0D0)RETURN
!
	DO I=MS,MF
	  NETCR(I)=NETCR(I)-ARTIFICIAL_HEAT_TERM(I)
	  TOTCR(I)=TOTCR(I)+ABS(ARTIFICIAL_HEAT_TERM(I))
	END DO
	WRITE(LU,'(/3X,A)')'Artificial heating term [ergs/cm**3/sec]'
	WRITE(LU,999)(ARTIFICIAL_HEAT_TERM(I),I=MS,MF)
C
999	FORMAT(1X,1P,10E12.4)
	RETURN
	END
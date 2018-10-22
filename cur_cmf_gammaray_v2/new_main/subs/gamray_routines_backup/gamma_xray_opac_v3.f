! Subroutine to compute the xray opacities
! My goal is calculate the xray opacity as function of frequency
! I have the compton opacities set up this way
!
!
	SUBROUTINE GAMMA_XCROSS_V3(NU,ND,CHI,E_SCAT_CHI,ED_TOT)
	USE XRAY_DATA_MOD
	USE MOD_CMFGEN
	IMPLICIT NONE
!
	REAL*8 :: CHI(ND)
	REAL*8 :: E_SCAT_CHI(ND)
	REAL*8 :: ED_TOT(ND)
	REAL*8 :: T1,T2,T3,T4
	REAL*8 :: NU 		!frequency in Hz (which will be converted for XCROSS_V2)
	REAL*8 :: XCROSS_V2
	REAL*8 :: CROSS_KN
	REAL*8, PARAMETER :: PLANCK = 4.135668E-021
	EXTERNAL XCROSS_V2
	EXTERNAL CROSS_KN
!
	INTEGER :: ND
	INTEGER :: I,J
	INTEGER :: ISPEC,ID
	INTEGER :: LUNIT
	INTEGER, PARAMETER :: IZERO=0
	INTEGER :: IOS
!
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
!
	CHI=0.0D0
	E_SCAT_CHI=0.0D0
	T4=CROSS_KN(NU)
	E_SCAT_CHI(1:ND)=E_SCAT_CHI(1:ND)+ED_TOT(1:ND)*T4
	DO ISPEC=1,NUM_SPECIES
!	  E_SCAT_CHI(1:ND)=E_SCAT_CHI(1:ND)+AT_NO(ISPEC)*POP_SPECIES(1:ND,ISPEC)*T4
	  DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)
	    IF (ATM(ID)%XzV_PRES) THEN
	      T2=AT_NO(SPECIES_LNK(ID))+1-ATM(ID)%ZXzV
	      T3=NU/1.0D+15
	      T1=XCROSS_V2(T3,AT_NO(SPECIES_LNK(ID)),T2,IZERO,IZERO,L_FALSE,L_FALSE)
	      IF(T1 .NE. 0.0D0 .AND. ID .NE. SPECIES_END_ID(ISPEC))THEN
	        DO I=1,ND
	          T2=0.0D0			!Temporary CHI
	          DO J=1,ATM(ID)%NXzV_F
		    T2=T2+ATM(ID)%XzV_F(J,I)
	          END DO
	          CHI(I)=CHI(I)+T1*T2
!		  E_SCAT_CHI(I)=E_SCAT_CHI(I)+AT_NO(ISPEC)*POP_SPECIES(I,ISPEC)*CROSS_KN(NU)
	        END DO
	      ELSE IF (T1 .NE. 0.0D0 .AND. ID .EQ. SPECIES_END_ID(ISPEC)) THEN
	        DO I=1,ND
	          T2=ATM(ID)%DXzV(I) 		!Temporary CHI
	          CHI(I)=CHI(I)+T1*T2
!		  E_SCAT_CHI(I)=E_SCAT_CHI(I)+AT_NO(ISPEC)*POP_SPECIES(I,ISPEC)*CROSS_KN(NU)
	        END DO
	      END IF
	    END IF
	  END DO
	END DO
!
	RETURN
	END SUBROUTINE
!
!*********************************************************************************************************
!
	FUNCTION CROSS_KN(NU)
	IMPLICIT NONE
!
	REAL*8 :: NU
	REAL*8 :: CROSS_KN 	! will return in units of cm^2
	REAL*8 :: EPSI
	REAL*8 :: Pi
	REAL*8 :: r_e
	REAL*8 :: r_e_squared
	REAL*8, PARAMETER :: r_c=2.42630694715D-10 ! h/(m_e*c) in units of cm
	REAL*8, PARAMETER :: ONE=1.0D0
	REAL*8, PARAMETER :: PLANCK = 4.135668E-021	! UNITS OF MeV*s SINCE PHOTON ENERGIES IN MeV
	REAL, PARAMETER :: M_elect = 0.5109989		! UNITS of MeV after multiply by c^2
	REAL*8, PARAMETER :: fine_struc=0.007297355257 	! Fine structure constant
!
	Pi=ACOS(-ONE)
	r_e=(fine_struc*r_c)/(2.0D0*Pi) ! units of cm
	r_e_squared=r_e*r_e!*(1.0D+18) ! units of cm^2 -x (old) megabarns x-
!
	EPSI = ( NU*PLANCK )/M_elect
!
	CROSS_KN = ( 2.0D0*Pi*r_e_squared )*( ((ONE+EPSI)/(EPSI*EPSI))*(2.0D0*(ONE+EPSI)/(ONE+2.0D0*EPSI)-LOG(ONE+2.0D0*EPSI)/EPSI)
	1		+LOG(ONE+2.0D0*EPSI)/(2.0D0*EPSI)-(ONE+3.0D0*EPSI)/(ONE+2.0D0*EPSI)/(ONE+2.0D0*EPSI) )
!
	RETURN
	END FUNCTION CROSS_KN

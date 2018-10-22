! Subroutine to evaluate the downscattered emissivity given an intensity
! at a specific frequency. This will use interpolation to sample the
! angular dependence in order to use Gauss-Chebyshev quadrature.
! 
! This routine is different than gamma_eta_sub_v*.f95 because it
! will calculate the emissivity of all the downscattered photons given 
! the specific intensity at a particular frequency.
!
! One thing to note is that for the integration, the angle grid
! goes from 1 to -1, so the first index always starts at 1.
! R_MU(DI)%MU_VECTOR(1)=1, R_MU(DI)%MU_VECTOR(last index)=-1
!
! Subject to modifications
!
!
! Created June 10, 2014
!
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
	SUBROUTINE ETA_SCAT_V6(NF_GRID_PTS,NU_GRID_VEC,IN_NU,NU_END, &
		ETA,INTENSITY,ND,NA,MON_NA,KLEIN_ARRAY,E_VEC,CHEB_ORDER,ED, &
		DO_NORM_ETA)
	USE GAM_MU_MOD
	IMPLICIT NONE
!
	INTEGER :: I,J,K,L,MS,MK,MON_MS
	INTEGER :: CHEB_ORDER
	INTEGER :: NF_GRID_PTS
	INTEGER :: ND 
        INTEGER :: LOC
        INTEGER :: IN_NU,FI,MI,DI ! FI=frequency index, MI=mu index, DI=depth index
	INTEGER :: IONE=1,ITWO=2
	INTEGER :: IOS
	INTEGER :: NU_END
	INTEGER :: NA,MON_NA
!
	LOGICAL :: DO_NORM_ETA
!
	REAL*8 :: NU_GRID_VEC(NF_GRID_PTS)
        REAL*8 :: ETA(NA,ND,NF_GRID_PTS) ! (MU,DEPTH,FREQ)
        REAL*8 :: INTENSITY(NA,ND,NF_GRID_PTS) ! (MU,ND,FREQ)
	REAL*8 :: MON_I(MON_NA,ND)
	REAL*8 :: KLEIN_ARRAY(NF_GRID_PTS,NF_GRID_PTS)
	REAL*8 :: ED(ND)
!
        REAL*8 :: C1,C2
        REAL*8 :: R_ELEC_SQRD
        REAL*8 :: A1
	REAL*8 :: A2
	REAL*8 :: GAM_DIFF
	REAL*8 :: ROOT1
	REAL*8 :: ROOT2
	REAL*8 :: B1,B2
        REAL*8 :: DISC
        REAL*8 :: CHEB_W
        REAL*8 :: X(CHEB_ORDER),X2(CHEB_ORDER)
        REAL*8 :: E_VEC(NF_GRID_PTS) ! will be h*nu/mec^2 thus unitless
	REAL*8 :: TEMPINT
	REAL*8 :: MU
	REAL*8 :: PI
	REAL*8 :: NU_PRIME_ON_NU
	REAL*8 :: FAC
	REAL*8 :: PRE_FAC(ND)
	REAL*8 :: T1,T2,T3,T4
	REAL*8 :: DELTA_NU
	REAL*8, DIMENSION(:), ALLOCATABLE :: INTERP
	REAL*8, DIMENSION(:), ALLOCATABLE :: INTERP2
	REAL*8, DIMENSION(:,:), ALLOCATABLE :: COEF
	REAL*8, DIMENSION(:,:), ALLOCATABLE :: COEF2
!
! Variables to norm for photon number.
!
	INTEGER :: NLEN
	REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: ETA_NORM
	REAL*8, DIMENSION(:), ALLOCATABLE :: I_NORM
	REAL*8, DIMENSION(:), ALLOCATABLE :: WRK
	REAL*8, DIMENSION(:,:), ALLOCATABLE :: WRK2
	REAL*8, DIMENSION(:), ALLOCATABLE :: NORM
	REAL*8 :: COMPTON_SIG
	EXTERNAL COMPTON_SIG
!
! In order to make the code faster, we can interpolate
! outside of the main loop and create a structure variable
! of length ND with components that will be the coefficients
!
	TYPE INTERP_COEF
	   REAL*8, DIMENSION(:,:), ALLOCATABLE :: C
	END TYPE
	TYPE(INTERP_COEF) :: Y(ND)
!
	REAL*8, PARAMETER :: ONE=1.0D0
	REAL*8, PARAMETER :: TWO=2.0D0
	REAL*8, PARAMETER :: PLANCK = 4.135668E-021 ! UNITS OF MeV*s SINCE PHOTON ENERGIES IN MeV
	REAL*8, PARAMETER :: R_ELEC = 2.8179403227D-13 ! Units of cm
	REAL*8, PARAMETER :: HoMC2 = 8.09330118D-21     ! HoMC2 = h/mc^2 in units of seconds
!
	MON_I=0.0D0
	IF(DO_NORM_ETA)THEN
	  NLEN=NU_END-IN_NU
	  ALLOCATE (ETA_NORM(NA,ND,NF_GRID_PTS)) ;  ETA_NORM=0.0D0
	  ALLOCATE (I_NORM(ND)); I_NORM=0.0D0
	  ALLOCATE (WRK(ND));    WRK=0.0D0
	  ALLOCATE (WRK2(NF_GRID_PTS,ND)); WRK2=0.0D0
	  ALLOCATE (NORM(ND));	NORM=0.0D0
	END IF
!
	PI=ACOS(-ONE)
	CHEB_W=PI/CHEB_ORDER
	R_ELEC_SQRD=R_ELEC*R_ELEC !*0.5D0
!
! ^^ 0.5D0 is removed from R_ELEC_SQRD because it is cancelled by a 2
! from having two roots for Dirac delta which each contribute the same
! term to the integral. See paper or notes. 
!
        DO K=1,CHEB_ORDER
          X(K)=COS(((TWO*K-ONE)*PI)/(TWO*CHEB_ORDER))
        END DO
!
	IF(IN_NU .LT. NF_GRID_PTS .AND. IN_NU .GT. 1)THEN
	  DELTA_NU=0.5D0*(NU_GRID_VEC(IN_NU-1)-NU_GRID_VEC(IN_NU+1))
	ELSE IF(IN_NU .EQ. NF_GRID_PTS)THEN
	  DELTA_NU=0.5D0*(NU_GRID_VEC(IN_NU-1)-NU_GRID_VEC(IN_NU))
	ELSE IF(IN_NU .EQ. 1)THEN
	  DELTA_NU=0.5D0*(NU_GRID_VEC(IN_NU)-NU_GRID_VEC(IN_NU+1))
	END IF
!
	DO DI=1,ND
	  MS=R_MU(DI)%MU_PTS
	  MON_MS=R_MU(DI)%MON_MU_PTS
	  ALLOCATE(INTERP(MS),COEF(MS,4))
	  ALLOCATE(INTERP2(MON_MS),COEF2(MON_MS,4))
!
	  INTERP(1:MS)=INTENSITY(1:MS,DI,IN_NU)
	  CALL MON_INT_FUNS_V2(COEF,INTERP,R_MU(DI)%MU_VECTOR,MS)
	  J=1
	  DO I=2,MON_MS-1
	    IF(R_MU(DI)%MON_MU(I) .LE. R_MU(DI)%MU_VECTOR(J) .AND. &
		R_MU(DI)%MON_MU(I) .GT. R_MU(DI)%MU_VECTOR(J+1))THEN
	      T2=R_MU(DI)%MON_MU(I)-R_MU(DI)%MU_VECTOR(J)
	      MON_I(I,DI)=COEF(J,4)+T2*(COEF(J,3)+T2*(COEF(J,2)+&
		COEF(J,1)*T2))
	    ELSE
	      DO WHILE(R_MU(DI)%MON_MU(I) .LT. R_MU(DI)%MU_VECTOR(J))
		J=J+1
		IF (J .EQ. MS) EXIT
	      END DO
	      J=J-1
	      T2=R_MU(DI)%MON_MU(I)-R_MU(DI)%MU_VECTOR(J)
	      MON_I(I,DI)=COEF(J,4)+T2*(COEF(J,3)+T2*(COEF(J,2)+&
		COEF(J,1)*T2))
	    END IF
	  END DO
	  MON_I(1,DI)=INTENSITY(1,DI,IN_NU)
	  MON_I(MON_MS,DI)=INTENSITY(MS,DI,IN_NU)
	  INTERP2(1:MON_MS)=MON_I(1:MON_MS,DI)
	  CALL MON_INT_FUNS_V2(COEF2,INTERP2,R_MU(DI)%MON_MU,MON_MS)
	  ALLOCATE(Y(DI)%C(MON_MS,4))
	  Y(DI)%C(1:MON_MS,1)=COEF2(1:MON_MS,1)
	  Y(DI)%C(1:MON_MS,2)=COEF2(1:MON_MS,2)
	  Y(DI)%C(1:MON_MS,3)=COEF2(1:MON_MS,3)
	  Y(DI)%C(1:MON_MS,4)=COEF2(1:MON_MS,4)
	  DEALLOCATE(INTERP,COEF)
	  DEALLOCATE(INTERP2,COEF2)
	END DO
!
!
!-------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------
! Starting the main loops for calculating ETA. Outermost loop will be frequency so I can update all downgraded photons 
!-------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------
!
	CALL TUNE(IONE,'Gamma scat full')
	DO FI=IN_NU+1,NU_END !Loop over outgoing frequency
	  NU_PRIME_ON_NU=NU_GRID_VEC(FI)/NU_GRID_VEC(IN_NU)
	  FAC=R_ELEC_SQRD*KLEIN_ARRAY(IN_NU,FI)*DELTA_NU*NU_PRIME_ON_NU
	  PRE_FAC(1:ND)=FAC*ED(1:ND)
	  GAM_DIFF=ONE+(E_VEC(FI)-E_VEC(IN_NU))/E_VEC(FI)/E_VEC(IN_NU)
!
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(DI,MS,MON_MS,I,T1,MU,& 
!$OMP& TEMPINT,A1,A2,ROOT1,ROOT2,DISC,C1,C2,X2,T2,J,LOC,B1,B2,K)
	  DO DI=1,ND
	    MS=R_MU(DI)%MU_PTS
	    MON_MS=R_MU(DI)%MON_MU_PTS
!
!-------------------------------------------------------------------------------------------------------------------------------
! This piece will handle the cases when MU'=1 and -1
!-------------------------------------------------------------------------------------------------------------------------------
!
! MU'=1 first
! FI is frequency index and the IN prefix means incoming, the MI means MU index
!
!	    MU=ONE-( E_VEC(IN_NU)-E_VEC(FI) )/(E_VEC(IN_NU)*E_VEC(FI))
	    MU=GAM_DIFF
	    IF (MU .LE. ONE .AND. MU .GE. -ONE) THEN
	      LOC=1+(1.0D0-MU)/R_MU(DI)%dMU
	      IF(LOC .GE. 1 .AND. LOC .LE. MON_MS)THEN
		T1=MU-R_MU(DI)%MON_MU(LOC)
		TEMPINT=( Y(DI)%C(LOC,4)+T1*(Y(DI)%C(LOC,3)&
			+T1*(Y(DI)%C(LOC,2)+Y(DI)%C(LOC,1)*T1)) )
	      ELSE
		WRITE(6,'(A,1X,I3)')'LOC:',LOC
		WRITE(6,'(A,1X)')'INDEX RANGE: 1',R_MU(DI)%MON_MU_PTS
		WRITE(6,'(A6,1X,I4,4X,A3,1X,I4)')'IN_NU:',IN_NU,'FI:',FI
		STOP "ERROR Failing with Y(DI)%C in ETA_SCAT_V3 out of range"
	      END IF
	    END IF ! NU is in range of scattering for emissivity at outgoing frequency
	    IF(DO_NORM_ETA)THEN
	      ETA_NORM(1,DI,FI)=PI*PRE_FAC(DI)*TEMPINT
	    ELSE
	      ETA(1,DI,FI)=ETA(1,DI,FI)+PI*PRE_FAC(DI)*TEMPINT
	    END IF
!
! MU'=-1
! FI is frequency index and the IN prefix means incoming, the MI means MU index
!
!	    MU=-(ONE-( E_VEC(IN_NU)-E_VEC(FI) )/(E_VEC(IN_NU)*E_VEC(FI)))
	    MU=-GAM_DIFF
	    IF (MU .LE. ONE .AND. MU .GE. -ONE) THEN
	      LOC=1+(1.0D0-MU)/R_MU(DI)%dMU
	      IF(LOC .GE. 1 .AND. LOC .LE. MON_MS)THEN
		T1=MU-R_MU(DI)%MON_MU(LOC)
		TEMPINT=( Y(DI)%C(LOC,4)+T1*(Y(DI)%C(LOC,3)&
			+T1*(Y(DI)%C(LOC,2)+Y(DI)%C(LOC,1)*T1)) )
	      ELSE
		WRITE(6,'(A,1X,I3)')'LOC:',LOC
		WRITE(6,'(A,1X)')'INDEX RANGE: 1',R_MU(DI)%MON_MU_PTS
		WRITE(6,'(A6,1X,I4,4X,A3,1X,I4)')'IN_NU:',IN_NU,'FI:',FI
		STOP "ERROR Failing with Y(DI)%C in ETA_SCAT_V3 out of range"
	      END IF
	    END IF ! NU is in range of scattering for emissivity at outgoing frequency
	    IF(DO_NORM_ETA)THEN
	      ETA_NORM(MS,DI,FI)=PI*PRE_FAC(DI)*TEMPINT
	    ELSE
	      ETA(MS,DI,FI)=ETA(MS,DI,FI)+PI*PRE_FAC(DI)*TEMPINT
	    END IF
!
!-------------------------------------------------------------------------------------------------------------------------------
! This section does all other cases
!-------------------------------------------------------------------------------------------------------------------------------
!
! FI is frequency index and the IN prefix means incoming, the MI means MU index
!
	    DO MI=2,MS-1
!
! This next line is to say that the incoming scattering frequency has to be less than or equal to the
! maximum frequency that can scatter to the outgoing frequency
!
!
! This next little piece is to to figure out the limits on the incoming mu integration
! which is a result of removing the phi dependence of the dirac delta 
!
! Original form was (-A2 +- SQRT(A2*A2 - 4*A1*A3))/(2*A1) However, A1 is
! -1 (see notes) or paper, so I can just ignore it and relabel A2->A1
! and A3->A2 for convenience.
!
!	      B1=TWO*R_MU(DI)%MU_VECTOR(MI)
!	      B2=R_MU(DI)%MU_VECTOR(MI)*R_MU(DI)%MU_VECTOR(MI)
!	      A1=B1*(ONE+ONE/E_VEC(IN_NU)-ONE/E_VEC(FI))
!	      A2=(ONE/E_VEC(IN_NU))*(TWO/E_VEC(FI)-ONE/E_VEC(IN_NU)-TWO)&
!			+(ONE/E_VEC(FI))*(TWO-ONE/E_VEC(FI))-B2
!	      DISC=A1*A1+4.0D0*A2 
	      B1=R_MU(DI)%MU_VECTOR(MI)
	      B2=R_MU(DI)%MU_VECTOR(MI)*R_MU(DI)%MU_VECTOR(MI)
	      DISC=(ONE-B2)*(ONE-GAM_DIFF*GAM_DIFF)
	      IF (DISC .LT. 0.0D0) THEN
		WRITE(6,*) "*** ERROR! NAN PROBLEM WITH DISCRIMINANT IN CALCULATING ETA SCAT ***"
		WRITE(6,'(A,ES16.6)')'DISCRIMINANT:',DISC
		WRITE(6,'(2(A3,1X,ES16.6,2X))')'A1:',A1,'A2:',A2
		WRITE(6,'(A3,1X,I3,2X,A,F10.6)')'MI:',MI,"MU':",R_MU(DI)%MU_VECTOR(MI)
		WRITE(6,'(A6,1X,I4,3X,A3,1X,I4)')'IN_NU:',IN_NU,'FI:',FI
		WRITE(6,'(A,1X,ES16.6,2X,ES16.6)')'Incoming/outgoing photon E [MeV]',&
			PLANCK*NU_GRID_VEC(IN_NU),PLANCK*NU_GRID_VEC(FI)
		STOP "***ERROR! NAN PROBLEM WITH DISCRIMINANT IN CALCULATING ETA SCAT ***"
	      ELSE
!
! ROOT1 to be the upper integration limit and ROOT2 will be the lower integration limit
! therefore, ROOT1 must be larger than ROOT2
!
		A1=B1*GAM_DIFF
		A2=-GAM_DIFF*GAM_DIFF+ONE-B2
		DISC=SQRT(DISC)
		IF(A1 .GT. 0.0D0)THEN
		  ROOT1=(A1+DISC)
		  ROOT2=-A2/ROOT1
		ELSE
		  ROOT2=(A1+DISC)
		  ROOT1=-A2/ROOT2
		END IF
!
! The quadratic formula guarantees that ROOT1 will be larger. If B1 is
! positive or negative, the discriminant is still positive, so we need
! to take the positive solution to make sure it is larger.
!
		C1=TWO/(ROOT1-ROOT2)
		C2=-(ROOT1+ROOT2)/(ROOT1-ROOT2)
		TEMPINT=0.0D0
		X2=(X-C2)/C1
!
! I need to find the locations for X2 to use with the interpolating coefficients
!
		DO K=1,CHEB_ORDER
		  J=1+(1.0D0-X2(K))/R_MU(DI)%dMU
		  T2=X2(K)-R_MU(DI)%MON_MU(J)
		  TEMPINT=TEMPINT+( Y(DI)%C(J,4)+T2*(Y(DI)%C(J,3)&
			+T2*(Y(DI)%C(J,2)+Y(DI)%C(J,1)*T2)) )
		END DO
!
		IF(DO_NORM_ETA)THEN
	          ETA_NORM(MI,DI,FI)=PRE_FAC(DI)*CHEB_W*TEMPINT
	        ELSE
	          ETA(MI,DI,FI)=ETA(MI,DI,FI)+PRE_FAC(DI)*CHEB_W*TEMPINT
		END IF
!
	      END IF
	    END DO ! (MI) out mu' loop
!	
	  END DO ! Over depth index
!$OMP END PARALLEL DO
	END DO ! Over outgoing frequencies
	CALL TUNE(ITWO,'Gamma scat full')
!
!*******************************************************************************************
! This next piece is meant to normalize and conserve the photons before
! and after scattering. (i.e. I*ED*SIGMA*dMU*dNU/NU = Eta/NU*dMU*dNU)
!*******************************************************************************************
!
      IF(DO_NORM_ETA .AND. NU_END .LT. NF_GRID_PTS)THEN
	DO DI=1,ND
	  MS=R_MU(DI)%MU_PTS
	  DO FI=IN_NU+1,NU_END
	    DO MI=1,MS-1
	      T1=R_MU(DI)%MU_VECTOR(MI)-R_MU(DI)%MU_VECTOR(MI+1)
	      WRK2(FI,DI)=WRK2(FI,DI)+0.5D0*T1*(ETA_NORM(MI,DI,FI)+&
				ETA_NORM(MI+1,DI,FI))
	    END DO ! (MI) out mu' loop
	  END DO ! Over FI index
	  DO FI=IN_NU+1,NU_END-1
	    T1=NU_GRID_VEC(FI)-NU_GRID_VEC(FI+1)
	    T2=NU_GRID_VEC(FI)
	    T3=NU_GRID_VEC(FI+1)
	    WRK(DI)=WRK(DI)+0.5D0*T1*(WRK2(FI,DI)/T2 + &
			WRK2(FI+1,DI)/T3)
	  END DO
	  DO MI=1,MS-1
	    T1=R_MU(DI)%MU_VECTOR(MI)-R_MU(DI)%MU_VECTOR(MI+1)
	    I_NORM(DI)=I_NORM(DI)+0.5D0*T1*(INTENSITY(MI,DI,IN_NU)&
			+INTENSITY(MI+1,DI,IN_NU) )
	  END DO
	  T1=NU_GRID_VEC(IN_NU)
	  T2=COMPTON_SIG(T1)
	  I_NORM(DI)=I_NORM(DI)*ED(DI)*T2*DELTA_NU/T1
	  IF(WRK(DI) .GT. 0.0D0)THEN
	    NORM(DI)=I_NORM(DI)/WRK(DI)
	  ELSE IF(NORM(DI) .EQ. 0.0D0)THEN
	    NORM(DI)=1.0D0
	  ELSE
	    WRITE(6,*)'NORM(DI):',NORM(DI)
	    WRITE(6,*)'WRK(DI):',WRK(DI)
	    WRITE(6,*)'I_NORM(DI):',I_NORM(DI)
	    STOP
	  END IF
	END DO
	DO FI=IN_NU+1,NU_END
	  DO DI=1,ND
	    MS=R_MU(DI)%MU_PTS
	    ETA(1:MS,DI,FI)=ETA(1:MS,DI,FI)+ETA_NORM(1:MS,DI,FI)*NORM(DI)
	  END DO
	END DO
	DEALLOCATE (ETA_NORM)
	DEALLOCATE (I_NORM)
	DEALLOCATE (WRK)
	DEALLOCATE (WRK2)
	DEALLOCATE (NORM)
      END IF
!
	RETURN
	END SUBROUTINE 
!
	FUNCTION COMPTON_SIG(NU)
	IMPLICIT NONE
!
	REAL*8 :: NU
	REAL*8 :: COMPTON_SIG
!
	REAL*8 :: EPSI
	REAL*8 :: X
	REAL*8 :: PI
	REAL*8 :: R_E2
!
	REAL*8, PARAMETER :: R_EL=2.8179403227D-13   ! cm
	REAL*8, PARAMETER :: ONE=1.0D0
	REAL*8, PARAMETER :: TWO=2.0D0
	REAL*8, PARAMETER :: H  =4.135668D-21	     ! MeV*sec
        REAL*8, PARAMETER :: M_E=5.109989D-1         ! Units of MeV
!
	PI=ACOS(-ONE)
	R_E2=R_EL*R_EL
!
	X=NU*H/M_E
	EPSI=ONE+TWO*X
!
	COMPTON_SIG=TWO*PI*R_E2*( ((ONE+X)/X/X/X)*( (TWO*X*(ONE+X))/EPSI-LOG(EPSI) ) &
			+ LOG(EPSI)/TWO/X - (EPSI+X)/EPSI/EPSI )
!
	RETURN	
	END FUNCTION

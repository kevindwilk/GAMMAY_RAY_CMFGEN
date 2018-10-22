! Subroutine to call the routines necessary to set up the
! gamma-ray routine.
!
!
	SUBROUTINE GAMRAY_SUB_V3(GAM_LU,ND,NP,NA,NA_MON,R,V,INST_DECAY)
	USE MOD_GAMMA_V3
        USE GAMMA_NUC_DECAY_V2
        USE GAM_MU_MOD
	USE MOD_RD_GAMRAY_CNTRL_VARIABLES
	IMPLICIT NONE
!
	INTEGER :: I,J,K
	INTEGER :: FI,ID
	INTEGER :: GAM_LU,IOS
	INTEGER :: ND,NP
	INTEGER :: NA,NA_MON
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: ITWO=2
	REAL*8 :: T1,T2,T3,T4,T5
	REAL*8 :: R(ND)
	REAL*8 :: V(ND)
	REAL*8 :: PI
	REAL*8 :: FOURPI
	LOGICAL :: INST_DECAY
	CHARACTER(LEN=30) :: FILENAME
	CHARACTER(LEN=30) :: OPTION
!
	REAL*8, PARAMETER :: PLANCK = 4.135668E-18 ! keV*s
	REAL*8, PARAMETER :: SOL = 299792 ! Units of km/s
	REAL*8, PARAMETER :: ERGS_TO_MEV = 624150.9
	REAL*8, PARAMETER :: HoMC2 = 8.09330118D-21  ! HoMC2 = h/mc^2 in units of seconds
!
	INTEGER :: LUER
	INTEGER :: ERROR_LU
	EXTERNAL :: ERROR_LU
!
	PI=ACOS(-1.0D0)
	FOURPI=4.0D0*PI
	LUER=ERROR_LU()
!
! RD_GAMRAY_CNTRL to read in the control parameters (GAMRAY_PARAMS) for the gamma-ray code
!
	FILENAME='GAMRAY_PARAMS'
	CALL RD_GAMRAY_CNTRL_V2(FILENAME)
!
	ALLOCATE (NU_GRID_VEC(NU_GRID_MAX))
!
	CALL GAM_NUC_DECAY_DATA_SUB_V3()
!
	ALLOCATE(NU_VEC(N_GAMMA))
!
	NU_VEC=0.0D0
	CALL NU_GRID_V15(N_GAMMA,NU_VEC,NU_GRID_VEC,NF_GRID_PTS)
!
! These next pieces are just diagnosing and supplemental pieces. They
! can be removed if needed.
!
	OPEN(UNIT=7,FILE='./data/scattering_diff.dat',STATUS='UNKNOWN',ACTION='WRITE',IOSTAT=IOS)
	  DO I=1,NF_GRID_PTS
	    WRITE(7,'(ES18.8,2X,ES18.8,2X,ES18.8)')NU_GRID_VEC(I),NU_GRID_VEC(I)*4.135668D-18,&
		(NU_GRID_VEC(I)-NU_GRID_VEC(I)/(1.0D0+2.0D0*GAM_E_ME(I)))/NU_GRID_VEC(I)*299792
	  END DO
	CLOSE(7)
!
	OPEN(UNIT=7,FILE='./data/velocity_step.dat',STATUS='UNKNOWN',ACTION='WRITE',IOSTAT=IOS)
	  DO I=1,NF_GRID_PTS-1
	    WRITE(7,'(ES18.8,2X,ES18.8,2X,ES18.8)')NU_GRID_VEC(I),NU_GRID_VEC(I)*4.135668D-18,&
		(NU_GRID_VEC(I)-NU_GRID_VEC(I+1))/NU_GRID_VEC(I)*299792
	  END DO
	CLOSE(7)
!
! ^
!
	ALLOCATE(NU_VEC_15(NF_GRID_PTS))
	NU_VEC_15(1:NF_GRID_PTS)=NU_GRID_VEC(1:NF_GRID_PTS)/1.0D15
	ALLOCATE(GAMMA_TAU(NF_GRID_PTS),XRAY_TAU(NF_GRID_PTS))
	ALLOCATE(NU_END(NF_GRID_PTS))
!	ALLOCATE(GAM_NU_MAX(NF_GRID_PTS))
!
	OPEN(UNIT=7,FILE='./data/gamma_nu_grid.dat',ACTION='WRITE',STATUS='UNKNOWN')
	CALL WRITV_V2(NU_GRID_VEC,NF_GRID_PTS,6,'Gamma ray freq grid',7)
	CLOSE(UNIT=7)
!
	OPEN(UNIT=GAM_LU,FILE='kevin_testing',STATUS='UNKNOWN',ACTION='WRITE',IOSTAT=IOS)
	  IF (IOS .NE. 0 ) STOP '*****Cannot open kevin_testing*****'
	CALL SET_LINE_BUFFERING(GAM_LU)
!
	ALLOCATE(GAM_E_ME(NF_GRID_PTS))
	ALLOCATE(KLEIN_ARRAY(NF_GRID_PTS,NF_GRID_PTS))
	ALLOCATE(ETA_ISO(ND,NF_GRID_PTS))
	ALLOCATE(GAMMA_FLUX(NF_GRID_PTS))
	ALLOCATE(DECAY_KIN_E(ND))
	ALLOCATE(ED_TOT(ND))
!
	GAMMA_FLUX=0.0D0
!
! KLEIN_SCAT is a subroutine to set up an array for the calculation of the scattering emissivity
!
	GAM_E_ME=home*NU_GRID_VEC
	CALL KLEIN_SCAT(KLEIN_ARRAY,NF_GRID_PTS,GAM_E_ME)
!
! This part is necessary for the scattering routine.
! For this next piece, I have to break it up into parts. The limit of
! backscattered photon energy is equal to half of the electron rest mass energy
! because of compton scattering. Therefore, I will take the upper limit for an outgoing
! frequency higher than half of the rest mass to be the first frequency in my vector
! which will represent infinity.
!
! Edit: this is no longer relevant since scattering emissivity is solved
! for all downscattered frequencies. 26-Oct-2017 KDW25
!
!	DO I=1,NF_GRID_PTS
!	  IF (GAM_E_ME(I) .LT. 0.5D0) THEN
!	    GAM_NU_MAX(I)=GAM_E_ME(I)/(1.0D0-2.0D0*GAM_E_ME(I))
!	  ELSE
!	    GAM_NU_MAX(I)=GAM_E_ME(1)
!	  END IF
!	END DO
!
! Calculating the intrinsic local emission from decays
!
	CALL TUNE(IONE,'ISO EMISSION')
	CALL GAMMA_INT_EMISS_V7(ETA_ISO,NU_GRID_VEC,NF_GRID_PTS,R,ND,V_GAUSS,&
				DECAY_KIN_E,NORM_GAM_LINES,INST_DECAY)
	CALL TUNE(ITWO,'ISO EMISSION')
!
	ALLOCATE(GAMRAY_EMISS(ND))
	ALLOCATE(TOTAL_DECAY_LUM(ND))
!
	GAMRAY_EMISS=0.0D0
	DO I=1,ND
	  DO K=1,NF_GRID_PTS-1
	    T1=NU_GRID_VEC(K)
	    T2=NU_GRID_VEC(K+1)
	    GAMRAY_EMISS(I)=GAMRAY_EMISS(I)+0.5D0*(T1-T2)*(ETA_ISO(I,K)+ETA_ISO(I,K+1))
	  END DO
	  TOTAL_DECAY_LUM(I)=GAMRAY_EMISS(I)*FOURPI+DECAY_KIN_E(I)
	END DO
	T3=0.0D0
	DO I=1,ND-1
	  T1=R(I)
	  T2=R(I+1)
	  T3=T3+0.5D0*(T1-T2)*(TOTAL_DECAY_LUM(I)*T1*T1+TOTAL_DECAY_LUM(I+1)*T2*T2)
	END DO
	T3=T3*FOURPI*1.0D+30
	GAMRAY_EMISS=GAMRAY_EMISS*FOURPI+DECAY_KIN_E
	OPEN(UNIT=7,FILE='./data/gamma_ray_local_emission.dat',STATUS='UNKNOWN',ACTION='WRITE',IOSTAT=IOS)
	  WRITE(7,'(A,ES16.6)')'Total Nuclear Luminosity in ergs/s:',T3
	  WRITE(7,'(A,ES16.6)')'Total Nuclear Luminosity in Lsun:',T3/3.826D+33
	  WRITE(7,'(A)')'!Gamma ray energies from the different lines and K.E. for each depth in ergs/cm^3/s'
	  WRITE(7,'(A)')'!Local emission is both line energies and positron/electron K.E.'
	  WRITE(7,'(A16,4X,A16,4X,A16)')'!Velocity km/s','Local Emiss','Decay K.E.'
	  DO I=1,ND
	    WRITE(7,'(ES18.8,2X,ES16.6,4X,ES16.6,4X,ES16.6)')R(I),V(I),GAMRAY_EMISS(I),DECAY_KIN_E(I)
	  END DO
	CLOSE(UNIT=7)
!
	FILENAME='./data/ETA_ISO_'
	OPTION='ETA_ISO'
	CALL WRITE_ARRAY_ISO(ETA_ISO,ND,NF_GRID_PTS,NU_GRID_VEC,FILENAME)
!
	CALL ELECTRON_DENSITY_CALC_V2(ED_TOT,ND)
	OPEN(UNIT=7,FILE='./data/electron_density.dat',ACTION='WRITE',STATUS='UNKNOWN')
	WRITE(7,'(A6,4X,A18)') 'Depth:','Electron density:'
	DO I=1,ND
	  WRITE(7,'(2X,I3,5X,ES16.6)')I,ED_TOT(I)
	END DO
	CLOSE(UNIT=7)
!
	NA=2*NP-1  ! Max number of mu angles
	NA_MON=NA*ANG_MULT ! Max number of linear mu grid angles
!
!
! NU_END is used to find the NU_GRID_VEC index of the back-scattered frequecny of all the frequencies in the grid
!
	OPEN(UNIT=7,FILE='./data/nu_end.dat',STATUS='UNKNOWN',ACTION='WRITE',IOSTAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Trouble opening ./data/nu_end.dat'
	  WRITE(LUER,*)'IOSTAT:',IOS
	  STOP
	END IF
	WRITE(7,'(A4,2X,A4,2X,5(A16,2X)))')'!  I','IEND','NU','NU_BACK_SCAT','NU(IEND-1)','NU(IEND)','NU(IEND+1)'
	DO I=1,NF_GRID_PTS
	  T1=NU_GRID_VEC(I)/(1.0D0+2.0D0*home*NU_GRID_VEC(I))
	  J=1
	  DO WHILE(NU_GRID_VEC(J) .GT. T1)
	    J=J+1
	  END DO
	  IF(I .GE. J-1)THEN
	    NU_END(I)=I
	  ELSE IF(I .LT. J-1)THEN
	    NU_END(I)=J-1
	  END IF
	  IF(NU_END(I) .EQ.1)THEN
	    WRITE(7,'(I4,2X,I4,2X,2(ES16.6,2X),A18,2(ES16.6,2X))')I,NU_END(I),NU_GRID_VEC(I),T1,' ',NU_GRID_VEC(NU_END(I)),NU_GRID_VEC(NU_END(I)+1)
	  ELSE IF(NU_END(I) .NE. 1 .AND. NU_END(I) .NE. NF_GRID_PTS)THEN
	    WRITE(7,'(I4,2X,I4,2X,5(ES16.6,2X))')I,NU_END(I),NU_GRID_VEC(I),T1,NU_GRID_VEC(NU_END(I)-1),NU_GRID_VEC(NU_END(I)),NU_GRID_VEC(NU_END(I)+1)
	  ELSE IF(I .EQ. NF_GRID_PTS .OR. NU_END(I) .EQ. NF_GRID_PTS)THEN
	    WRITE(7,'(I4,2X,I4,2X,4(ES16.6,2X),A18)')I,NU_END(I),NU_GRID_VEC(I),T1,NU_GRID_VEC(NU_END(I)-1),NU_GRID_VEC(NU_END(I)),' '
	  END IF
	END DO
	CLOSE(UNIT=7)
!
	ALLOCATE(GAM_I(NA,ND,NF_GRID_PTS),GAM_ETA_SCAT(NA,ND,NF_GRID_PTS))
	ALLOCATE(GAM_I_T(ND,NA,NF_GRID_PTS),GAM_ETA_SCAT_T(ND,NA,NF_GRID_PTS))
	ALLOCATE(GAMMA_J(NF_GRID_PTS,ND))
	ALLOCATE(GAM_ETA_MUAVG(NF_GRID_PTS,ND))
! ^^^^  "_T" means transpose on arrays to use them in the transfer routine  ^^^^
!
	GAM_I=0.0D0 !gamma-ray intensity
	GAM_I_T=0.0D0 !gamma-ray intensity
	GAM_ETA_SCAT=0.0D0 !gamma-ray scattering array
	GAM_ETA_SCAT_T=0.0D0 !gamma-ray scattering array
!
	ALLOCATE(GAM_ETA(ND,NA),GAM_INT(ND,NA))
	GAM_ETA=0.0D0
	GAM_INT=0.0D0
!
	ALLOCATE(XCHI(ND),E_SCAT_CHI(ND),GAM_OPAC(ND),GAM_OPAC_COPY(ND))
	ALLOCATE(GAM_OPAC_CLUMP(ND))
	ALLOCATE(ETA_NORM_CONST(ND))
	ALLOCATE(E_SCAT_ARRAY(NF_GRID_PTS,ND))
	ALLOCATE(CHI_ARRAY_ABS(NF_GRID_PTS,ND))
	ALLOCATE(GAM_ENERGY_DEP(ND))
	ALLOCATE(GAMRAY_LUM_VEC(ND))
!	ALLOCATE(NU_AVG_I(NA,ND))
	ALLOCATE(GAM_ETA_NUAVG(NA,ND))
!
	GAM_ENERGY_DEP=0.0D0
	ETA_NORM_CONST=1.0D0
	XCHI=0.0D0
	E_SCAT_CHI=0.0D0
	GAM_OPAC=0.0D0
	GAM_OPAC_COPY=0.0D0
	GAM_OPAC_CLUMP=0.0D0
!
	OPEN(UNIT=7,FILE='./data/E_SCAT_ARRAY',STATUS='UNKNOWN',ACTION='WRITE')
        WRITE(7,*)"!Opened E_SCAT_ARRAY"
        CLOSE(UNIT=7)
        OPEN(UNIT=7,FILE='./data/GAMRAY_E_DEP',STATUS='UNKNOWN',ACTION='WRITE')
        WRITE(7,*)"!Opened GAMRAY_E_DEP"
        CLOSE(UNIT=7)
        OPEN(UNIT=7,FILE='./data/gamma_ray_lum.dat',STATUS='UNKNOWN',ACTION='WRITE')
        WRITE(7,*)"!Opened gamma_ray_lum.dat"
        CLOSE(UNIT=7)
        OPEN(UNIT=7,FILE='./data/gamma_ray_lum_J.dat',STATUS='UNKNOWN',ACTION='WRITE')
        WRITE(7,*)"!Opened gamma_ray_lum.dat"
        CLOSE(UNIT=7)
        OPEN(UNIT=7,FILE='./data/ray1_intensity.dat',STATUS='UNKNOWN',ACTION='WRITE')
        WRITE(7,*)'!Opened ray1_intensity.dat'
        CLOSE(UNIT=7)
        OPEN(UNIT=7,FILE='./data/DIAGN_EDEP',STATUS='UNKNOWN',ACTION='WRITE')
        WRITE(7,*)'!Opened DIAGN_EDEP'
        CLOSE(UNIT=7)
        OPEN(UNIT=7,FILE='./data/photons.dat',STATUS='UNKNOWN',ACTION='WRITE')
        WRITE(7,*)'!Opened photons.dat'
        CLOSE(UNIT=7)
        OPEN(UNIT=7,FILE='./data/TAU_RAY.dat',STATUS='UNKNOWN',ACTION='WRITE')
        WRITE(7,*)'!Opened TAU_RAY.dat'
        WRITE(7,*)'ND:',ND
        CLOSE(UNIT=7)
	WRITE(LUER,*)'Set up gamma-ray transfer routine variables...'
!
!************************************************************************************
!************************************************************************************
!
! Now to solve the radiative transfer equation
!
!************************************************************************************
!************************************************************************************
!
	FIRST_FREQ=.FALSE.
	dLOG_NU=0.0D0
	DO FI=1,NF_GRID_PTS
	  WRITE(KW_LU,'(A5,I7)')"FI:",FI
	  INNER_BAND_METHOD='ZERO_FLUX'
	  BNUE=0.0D0
	  GAM_INT=0.0D0
	  GAM_ETA=0.0D0
	  FL=NU_GRID_VEC(FI)/1.0D15
	  IF(FI .EQ. 1)THEN
	    FIRST_FREQ=.TRUE.
	  ELSE
	    FIRST_FREQ=.FALSE.
	    dLOG_NU=DLOG(NU_GRID_VEC(FI-1)/NU_GRID_VEC(FI))
	  END IF
	  NEW_FREQ=.TRUE.
!
! Setting up the opacities each frequencies
!
!         CALL GAMMA_XCROSS_V3(NU_GRID_VEC(FI),ND,XCHI,E_SCAT_CHI,ED_TOT)
	  CALL GAMMA_XCROSS_V4_TEST(NU_GRID_VEC(FI),ND,XCHI,E_SCAT_CHI,ED_TOT)
!
! E_SCAT_CHI is not in units CMFGEN uses so I multiply by the factor of 10^10
!
	  GAM_OPAC(1:ND)=XCHI(1:ND) + 1.0D+10*E_SCAT_CHI(1:ND)
	  GAM_OPAC_CLUMP(1:ND)=GAM_OPAC(1:ND)*CLUMP_FAC(1:ND)
	  GAM_OPAC_COPY=GAM_OPAC
!
! I store the separate opacities as a function of frequency and depth to
! use later on after Intensity is solved. Note these arrays are in cm^-1
!
	  E_SCAT_ARRAY(FI,1:ND)=E_SCAT_CHI(1:ND)
	  CHI_ARRAY_ABS(FI,1:ND)=XCHI(1:ND)*1D-10
!
!
! Using the array GAM_ETA to be the 2D array to pass to the subroutines
! along with GAM_INT (for gamma intensity). Initialized the total
! emissivity as the scattered plus always the isotropic emission. 
!
	  DO J=1,NA ! NA=2*NP-1
	    DO I=1,ND
	      GAM_ETA(I,J)=ETA_ISO(I,FI)+GAM_ETA_SCAT(J,I,FI)
	    END DO
	  END DO
	  GAM_ETA=GAM_ETA*1.0D+10
!
! Now we solve for the intensity for all rays for the current frequency
! We are solving it from blue to red since Compton scattering downgrades
! photon energy.
!
	  CALL TUNE(IONE,'CMF_FORMAL SOLVER')
	  CALL CMF_FORMAL_REL_V4_GAM( & 
		GAM_OPAC_COPY,GAM_OPAC_CLUMP,CHI_SCAT_CLUMP,V,SIGMA,R,P, &
		GAM_OPAC_COPY,FEDD,HFLUX_AT_IB,HFLUX_AT_OB,IPLUS, &
		FL,dLOG_NU,BNUE,DBB, &
		INNER_BAND_METHOD,THK_CONT, &
		VDOP_VEC,DELV_FRAC_FG,REXT_FAC, &
		METHOD,FIRST_FREQ,NEW_FREQ,NC,NP,ND, &
		NA,FI,GAM_ETA,GAM_INT,ANG_MULT)
	  CALL TUNE(ITWO,'CMF_FORMAL SOLVER')
!
! We need to save the intensity from GAM_INT(ND,NA)
!
	  DO I=1,ND
	    DO J=1,NA
	      GAM_I(J,I,FI)=GAM_INT(I,J)
	    END DO
	  END DO
!
! Next we calculate the scattering emissivity from all downscattered
! photons from the current frequency (assuming no emission at the
! current frequency)
!
	  CALL ETA_SCAT_V4_BETA(NF_GRID_PTS,NU_GRID_VEC,FI,NU_END(FI),&
		GAM_ETA_SCAT,GAM_I,ND,NA,NA_MON,KLEIN_ARRAY,&
		GAM_E_ME,GAM_NU_MAX,CHEB_ORDER,ED_TOT)
!
	END DO ! FI (over frequency)
!
!************************************************************************************
!************************************************************************************
!
! We have now solved for the specific intensity, so we can now get an
! observers frame flux. Note this code currently uses the same CMF
! frequencies as it does the OBS frequencies
!
	ALLOCATE(GAMMA_NU_STORE(NF_GRID_PTS))
	ALLOCATE(GAM_I_STORE(NF_GRID_PTS,NP))
	GAMMA_NU_STORE=0.0D0
	GAM_I_STORE=0.0D0
	FIRST_OBS_COMP=.TRUE.
	OBS_REL_FULL=.TRUE.
	DO I=1,NF_GRID_PTS
	  FL=NU_VEC_15(I)
	  CALL COMP_OBS_V2(GAM_I(1:NP,1,I),FL, &
		GAM_I_STORE,GAMMA_NU_STORE,NF_GRID_PTS, &
		MU_AT_RMAX,HQW_AT_RMAX,NU_VEC_15,GAMMA_FLUX,NF_GRID_PTS, &
		V(1),R(1),'IPLUS','LIN_INT',OBS_REL_FULL,
		FIRST_OBS_COMP,NP)
	END DO
!
! The gamma-ray flux is now in Janskies for 1 kpc. I will output this
! the same as OBSFLUX does for CMFGEN
!
! Later, I will try add an option in plt_spec and like files to change the units to
! photons/cm^2/s/keV when using rd_mod in plt_spec. For now use
! cnvrt_gamflux.exe 
!
	CALL GEN_ASCI_OPEN(LU_GAM,'GAMFLUX','UNKNOWN',' ',' ',IZERO,IOS)
	WRITE(STRING,'(I10)')NF_GRID_PTS
	DO WHILE(STRING(1:1) .EQ. ' ')
	  STRING(1:)=STRING(2:)
	END DO
	STRING='Continuum Frequencies ( '//TRIM(STRING)//' )'
	CALL WRITV_V2(NU_VEC_15,NF_GRID_PTS,ISIX,TRIM(STRING),LU_GAM)
	CALL WRITV_V2(GAMMA_FLUX,NF_GRID_PTS,IFOUR,&
		'Observed intensity (Janskys)',LU_GAM)
	CLOSE(UNIT=LU_GAM)
!
!
!------------------------------------------------------------------------------
! Now I write some variables to check the outputs
!------------------------------------------------------------------------------
!
! Checking the optical depth to both x-rays and gamma-rays 
!
	GAMMA_TAU=0.0D0
	XRAY_TAU=0.0D0
	DO K=1,NF_GRID_PTS
	  DO I=1,ND-1
	    T1=E_SCAT_ARRAY(K,I)
	    T2=E_SCAT_ARRAY(K,I+1)
	    T3=CHI_ARRAY_ABS(K,I)
	    T4=CHI_ARRAY_ABS(K,I+1)
	    GAMMA_TAU(K)=GAMMA_TAU(K)+0.5D0*(R(I)-R(I+1))*(T1+T2)
	    XRAY_TAU(K)=XRAY_TAU(K)+0.5D0*(R(I)-R(I+1))*(T3+T4)
	  END DO
	END DO
	GAMMA_TAU=GAMMA_TAU*1.0D10
	XRAY_TAU=XRAY_TAU*1.0D10
	OPEN(UNIT=7,FILE='./data/TAU_gam_xray.dat',STATUS='UNKNOWN',ACTION='WRITE')
	WRITE(7,'(A1,1X,A12,2X,A14,2X,A14,)')'!','NU (keV)','TAU_GAM','TAU_XRAY'
	DO I=1,NF_GRID_PTS
	  WRITE(7,'(2X,F14.6,ES16.6,ES16.6)')NU_GRID_VEC(I)*PLANCK,GAMMA_TAU(I),XRAY_TAU(I)
	END DO
	CLOSE(UNIT=7)
!
! Now we calculate the energy deposited from gamma-ray scattering.
! This calculates the difference between integral (chi*I-eta) dnu*dOmega
!
	CALL GAMRAY_ENERGY_DEPOSIT_V5(GAM_I,GAM_ETA_SCAT,E_SCAT_ARRAY,CHI_ARRAY_ABS, &
		GAM_ENERGY_DEP,DECAY_KIN_E,NU_GRID_VEC,NF_GRID_PTS,V,R,NA,ND,GAMRAY_EMISS)
!
	GAMMA_J=0D0
	DO K=1,NF_GRID_PTS
	  DO I=1,ND
	    ID=R_MU(I)%MU_PTS
	    DO J=1,ID-1
	      T1=R_MU(I)%MU_VECTOR(J)
	      T2=R_MU(I)%MU_VECTOR(J+1)
	      GAMMA_J(K,I)=GAMMA_J(K,I)+5.0D-1*(T1-T2)*(GAM_I(J,I,K)+GAM_I(J+1,I,K))
	    END DO
	  END DO
	END DO
	GAMMA_J=GAMMA_J/2.0D0
	FILE_NAME='./data/GAMMA_J_'
	OPTION='FREQ'
	CALL WRITE_ARRAY_V2(GAMMA_J,ND,NF_GRID_PTS,NU_GRID_VEC,FILE_NAME,COUNTER,OPTION)
!
	GAM_ETA_MUAVG=0D0
	DO K=1,NF_GRID_PTS
	  DO I=1,ND
	    ID=R_MU(I)%MU_PTS
	    DO J=1,ID-1
	      T1=R_MU(I)%MU_VECTOR(J)
	      T2=R_MU(I)%MU_VECTOR(J+1)
	      GAM_ETA_MUAVG(K,I)=GAM_ETA_MUAVG(K,I)+5.0D-1*(T1-T2)*( &
			GAM_ETA_SCAT(J,I,K)+GAM_ETA_SCAT(J+1,I,K))
	    END DO
	  END DO
	END DO
	FILE_NAME='./data/ETA_MUAVG_'
	OPTION='FREQ'
	CALL WRITE_ARRAY_V2(GAM_ETA_MUAVG,ND,NF_GRID_PTS,NU_GRID_VEC,FILE_NAME,COUNTER,OPTION)
!
	CALL TUNE(3,'')
	CLOSE(KW_LU)
	WRITE(LUER,'(A)')'Finished gamma-ray scattering and deposition'
	WRITE(LUER,'(A)')'Leaving GAMRAY_SUB_V3'
!
	RETURN
	END SUBROUTINE

!
!======================================================================
!
!*******3D NON-SPHERICAL PARTICLE TRACKER SUBROUTINES BEGIN************
!
!======================================================================
!
SUBROUTINE PART_MAIN
!
USE  FDIM_MODULE
USE  FGRID_MODULE
USE  FVAR_MODULE
USE  FTURB_MODULE
USE  FTIME_MODULE
USE  FFILE_MODULE
USE  FLOGICAL_MODULE
USE  FCPUTIME_MODULE
USE  FPROP_MODULE
USE  PART_NUM_MODULE
USE  PART_COUNT_MODULE
USE  PART_INITIAL_MODULE
USE  PART_LOGICAL_MODULE
USE  PART_STAT_MODULE
USE  PART_VECTOR_MODULE
USE  PART_TIME_MODULE
USE  PART_PROP_MODULE
USE  PART_EXTRA_MODULE
USE  PART_SGS_MODULE
USE  GEO_SIZE_MODULE
USE  FROTAT_MODULE
USE  PART_PWI_MODULE
USE  PART_PDF_MODULE
USE  PART_CONCEN_MODULE
USE  FUTILITY_MODULE
USE  PART_TWOWAY_MODULE
USE  PART_DEPOSIT_MODULE
USE  PART_PSTATE_MODULE
USE  PART_RSTATE_MODULE
USE  PART_FLAG_MODULE
USE  PART_COLLISION_MODULE

!
IMPLICIT NONE
INTEGER ::  I,LOOP
REAL,EXTERNAL:: GASDEV,RAN1
!
!EXECUTE
!:::::::::::::::::::::::::::::::::::::::::
!ILPT COUNTS NUMBER TIMES BOFFIN CALLS LPT
!:::::::::::::::::::::::::::::::::::::::::

ILPT = ILPT + 1 ! UPDATE NUMBER OF TIMES LPT HAS BEEN CALLED

IF(TWOWAY)THEN ! DO WE HAVE TWO WAY COUPLING ON?
   CALL PART_ALLOCATE(21)
   FX = 0.0
   FY = 0.0
   FZ = 0.0
END IF

DO LOOP = 1, LPT_STEP
   !:::::::::::::::::
   !INTEGERATION TIME
   !:::::::::::::::::
   PTIM = DTIM
   IF(FIRST)THEN
      WRITE(SCRN,*)'LPT INSIDE FIRST'
      !:::::
      !ZEROS
      !:::::
      LPTSTEP = 0
      LPTTIME = 0.0
      SPRAYTIME = 0.0
      PART_STEP  = 0
      INJECTION 	= 	.TRUE.
      FIRST 		= 	.FALSE.
      REPLACE 		= 	.TRUE.
      !
      !*********************
      CALL PART_ALLOCATE(20)
      !*********************
      !INITIALISE ROTATION AND ORIENTATION
      !--------------
      CALL PART_INPUT
      !--------------
      !::::::::::::::
      !INITIAL VALUES
      !::::::::::::::
      NPCT = LPTTIME*WTIM/FFT

      IMIN_CELL = LP1
      IMAX_CELL = 1
      KMIN_CELL = NP1
      KMAX_CELL = 1
   ELSE
      INDEX = 0
      PART_STEP = PART_STEP + 1
   END IF
   !
   IF(MOD(LPTSTEP,INJECT_STEP) == 0 .OR. INJECTION)THEN
    !****************************
   	!INJECTION OF PARTICLE BEGINS
   	!****************************
   	!-----------------
   	CALL PART_INJECTOR
   	!-----------------
   END IF
      LPTTIME = LPTTIME + PTIM
   IF(PART_INTERACT)THEN
   	!*********************************
   	!CALL PARTICLE-PARTICE INTERACTION
   	!*********************************
   	!----------------------
      CALL PART_COLLISION
   	!----------------------
      MAXNUMCOL = 0
      IPOPULATION = 0
   END IF
   !
   !*********************

   !*********************
   DO ISET = 1,NPSET
      IMIN_CELL(ISET) = LP1
      KMIN_CELL(ISET) = NP1
      PTIM = DTIM/TIMEFAC(ISET)
      DO IP = 1,NPAR
         JPART = IP + (ISET - 1)*NPART
         IF(INDEX_OUT(JPART) == 1)THEN
            CYCLE
         END IF
         INDEX(ISET) = INDEX(ISET) + 1
         IF(PTIM*WTIM > PST_SET(ISET))THEN
            !---------------------
            CALL PART_MESSAGE(400)
            !---------------------
         END IF
         DO I = 1,ND
            PART0(I,   JPART) = PART(I,   JPART)
            PART0(I+ND,JPART) = PART(I+ND,JPART)
            !FLUID(I+ND,JPART) = FLUID(I,  JPART)
         END DO
         !
         !SOLVE FOR PARTICLE TRANSLATION
         !----------------
         CALL PART_BBO_YIN
         !----------------
         !
         IF(LAMBDA_SET(ISET) /= 1)THEN
            !SOLVE FOR PARTICLE ORIENTATION AND ROTATION
            !----------------
            CALL PART_ROT_YIN
            !----------------
         END IF
         !----------------------
         CALL PART_PWI_ELLIPSOID
         !----------------------
         !:::::::::::::::::::::::::::::::::::::::::::::::::::
         !UPDATE THE TOTAL INTEGRATION TIME FOR EACH PARTICLE
         !:::::::::::::::::::::::::::::::::::::::::::::::::::
         PART(7,JPART) = PART(7,JPART) + PART0(10,JPART)
         !-------------
         CALL PART_FIND
         !-------------
         IF(DEPOSITION)THEN
            DO I = 1, DSLAB
               IF(PART0(3,JPART) <= ZSLAB(I) .AND. PART(3,JPART) > ZSLAB(I))THEN
                  DEPOSITNUM(ISET,I) = DEPOSITNUM(ISET,I) + 1
               END IF
            END DO
            IF(PART0(3,JPART) <= ZSLAB(DSLAB+1) .AND. PART(3,JPART) > ZSLAB(DSLAB+1))THEN
               DEPOSITNUM(ISET,DSLAB+1) = DEPOSITNUM(ISET,DSLAB+1) + 1
            END IF
            IF(IMIN_CELL(ISET) > IPS(1,JPART))THEN
               IMIN_CELL(ISET) = IPS(1,JPART)
            END IF
            IF(IMAX_CELL(ISET) < IPS(1,JPART))THEN
               IMAX_CELL(ISET) = IPS(1,JPART)
            END IF
            IF(KMIN_CELL(ISET) > IPS(3,JPART))THEN
               KMIN_CELL(ISET) = IPS(3,JPART)
            END IF
            IF(KMAX_CELL(ISET) < IPS(3,JPART))THEN
               KMAX_CELL(ISET) = IPS(3,JPART)
            END IF
         END IF
         !
         IF(PART_INTERACT)THEN
            !COMPUTE MAXIMUM NUMBER OF TIMES PARTICLE-PARTICLE INTERACTION HAS OCCURED
         	IF(MAXNUMCOL < INDEX_COL_COUNT(JPART))THEN
            	MAXNUMCOL = INDEX_COL_COUNT(JPART)
         	END IF
            !COUNT THE NUMBER OF "SINGLET", "DOUBLET", "TRIPLET", "QUADRUPLETS", ETC
            IPOPULATION(INDEX_COL_COUNT(JPART)) = IPOPULATION(INDEX_COL_COUNT(JPART)) + 1
         END IF
      END DO
      IF(DEPOSITION)DEPOSITNUM(ISET,1) = NPAR
      SPRAYTIME(ISET) = SPRAYTIME(ISET) + PTIM
   END DO
   !INCREASE TOTAL PARTICLE INTEGRATION TIME STEP
   !LPTTIME = LPTTIME + DTIM

   !INCREASE TIME STEP
   LPTSTEP = LPTSTEP + 1
   !NPCT = NUMBER OF PARTICLE CROSSING TIME (IN WALL UNIT)
   NPCT = LPTTIME*WTIM/FFT
   !-----------------------
   !ANIMATION AND DEBURGING
   !-----------------------
   IF( MOD(LPTSTEP,DEBUGSKIP) == 0)THEN
      CALL PART_INFORM
   END IF
   !---------------------------------------
   !SAMPLE PARTICLE FOR  STATISTICS AND PDF
   !---------------------------------------
   IF(SAMPLE)THEN
      CALL PART_SAMPLING
      IF(MOD(LPTSTEP,ISAVE) == 0)THEN
        !OUTPUT SAMPLED STATISTICS
         CALL SAVE_OUTPUT(21)
      END IF
   END IF
   !----------------------
   !PARTICLE CONCENTRATION
   !----------------------
   IF(PARTCONC .AND. MOD(LPTSTEP,CONCSTEP) == 0)THEN
		!OUTPUT PARTICLE PREFERENCIAL CONCENTRATION
   	CALL PART_CONC
	END IF
   !
   CALL CPU_TIME(USERTIM2)
   ELAPSEDTIME = USERTIM2 - USERTIM0


   IF(TEXEC > 0.0)THEN
      IF(ELAPSEDTIME >= 47.0*60.0*60.0 .OR. ELAPSEDTIME >= TEXEC)THEN
         WRITE(SCRN,'(3(A,1PE20.12))')'USERTIM1=',USERTIM1,' USERTIM2=',USERTIM2,' ELAPSEDTIME =',ELAPSEDTIME
     		CALL SAVE_OUTPUT(0)
     		WRITE(SCRN,'(2(A,F12.4))') &
     		'BOFFIN STOPS ALOCATED CPU RUN, TEXEC = ',TEXEC, ' USER TTIME = ',ELAPSEDTIME
     		CALL DHMS(ELAPSEDTIME,SCRN)
     		STOP
      END IF
   ELSE IF(LPTSTEP > MAXSTEP) THEN
      !FINISHED LPT
      CALL SAVE_OUTPUT(0)
      WRITE(SCRN,*)'*** LAGRANGIAN PARTICLE TRACKING IS COMPLETED ***'
      WRITE(SCRN,*)'MAXIMIUM TIME STEP IS REACHED: LPTSTEP= ',LPTSTEP
      STOP
   END IF
   IF(ALL(INDEX(1:NPSET) == 0))THEN
      WRITE(SCRN,*)'*** DEPOSITION IS COMPLETED ***'
      DO ISET = 1, NPSET
      	WRITE(SCRN,'(2(A,I8))')'ISET=',ISET,'   INDEX(ISET)=',INDEX(ISET)
      END DO
      DEPACTIVE = 1
      CALL SAVE_OUTPUT(0)
      WRITE(SCRN,'(2(A,F12.4))') &
     'BOFFIN STOPS ALOCATED CPU RUN, TEXEC = ',TEXEC, ' USER TTIME = ',ELAPSEDTIME
     CALL DHMS(ELAPSEDTIME,SCRN)
     STOP
   END IF
END DO
IF(TWOWAY .AND. ALLOCATED(FORCE))CALL DEALLOCATE_TWOWAY_FORCE
RETURN
END
!-------------------
!END OF MAIN PROGRAM
!-------------------
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!							DEFINITION OF VARIABLES
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!
!-----------------------------------------------------------------------
!		PART(1,JPART)			X				NEW PARTICLE POSITION (m)
!		PART(2,JPART)			Y				NEW PARTICLE POSITION (m)
!		PART(3,JPART)			Z				NEW PARTICLE POSITION (m)
!		PART(4,JPART)			U				NEW PARTICLE VELOCITY (m/s)
!		PART(5,JPART)			V				NEW PARTICLE VELOCITY
!		PART(6,JPART)			W				NEW PARTICLE VELOCITY
!		part(7,jpart)			TDT			TOTAL Dt for each particle (s)
!		part(8,jpart)			dp				PARTICLE DIAMETER (m)
!		part(9,jpart)			rad			PARTICLE RADIUS, a (m); used for collision
!		part(10,jpart)			mp				PARTICLE MASS (kg)
!
!		............................................................
!
!		PART0(1,JPART)			X				OLD PARTICLE POSITION (m)
!		PART0(2,JPART)			Y				OLD PARTICLE POSITION (m)
!		PART0(3,JPART)			Z				OLD PARTICLE POSITION (m)
!		PART0(4,JPART)			U				OLD PARTICLE VELOCITY (m/s)
!		PART0(5,JPART)			V				OLD PARTICLE VELOCITY
!		PART0(6,JPART)			W				OLD PARTICLE VELOCITY
!		PART0(7,JPART)			a				OLD PARTICLE ACCELERATION IN X DIRECTION (m/s2)
!		PART0(8,JPART)			a				OLD PARTICLE ACCELERATION IN y DIRECTION (m/s2)
!		PART0(9,JPART)			a				OLD PARTICLE ACCELERATION IN z DIRECTION (m/s2)
!		PART0(10,JPART)		DT				PARTICLE INSTANT INTEGRATION TIME (s)
!
!     ............................................................
!
!		PIN(01,JPART)			USLIP			PARTICLE SLIP VELOCITY (m/s)
!		PIN(02,JPART)			Rep			PARTICLE REYNOLDS NUMBER
!		PIN(03,JPART)			CD				PARTICLE DRAG COEFFICIENT
!		PIN(04,JPART)			TAUP			PARTICLE RELAXATION TIME (S)
!		PIN(05,JPART)			ALPHA			PARTICLE INCIDENT ANGLE (rad)
!		PIN(06,JPART)			SEFFD			PARTICLE PROJECTED AREA (m2)
!		PIN(07,JPART)			Cos(Tx)		DIRECTION COSINE (COSINE(TitaX)
!		PIN(08,JPART)			Cos(Ty)		DIRECTION COSINE (COSINE(TitaY)
!		PIN(09,JPART)			Cos(Tz)		DIRECTION COSINE (COSINE(TitaZ)
!		PIN(10,JPART)			AFA			ANGLE BTW Z' AND PLANE YZ-AXIS (WALL) (rad)
!		PIN(11,JPART)			DIAMETER		PARTICLE DYNAMIC DIAMETER DUE TO COLLISION AND AGGLOMERATION (m)
!
!		.............................................................
!
!		ROT(1,JPART)			q1				Euler Parameter 1
!		ROT(2,JPART)			q2				Euler Parameter 2
!		ROT(3,JPART)			q3				Euler Parameter 3
!		ROT(4,JPART)			wx				ANGULAR VELOCITY X
!		ROT(5,JPART)			wy				ANGULAR VELOCITY y
!		ROT(6,JPART)			wz				ANGULAR VELOCITY z
!
!		............................................................
!
!		PSTATE(01,ISET)		LAMBDA_SET(ISET)		PARTICLE ASPECT RATIO
!		PSTATE(02,ISET)		RADA						PARTICLE SEMI-MINOR AXIS (m)
!		PSTATE(03,ISET)		ESDV						PARTICLE EQUIVALENT VOLUME OF SPHERE DIAMETER (m)
!		PSTATE(04,ISET)		PST_SET(ISET)			PARTICLE STOKES NUMBER
!		PSTATE(05,ISET)		MASS						PARTICLE MASS (kg)
!		PSTATE(06,ISET)		DENP_SET(ISET)			PARTICLE DENSITY (kg/m3)
!		PSTATE(07,ISET)		PS							PARTICLE SPHERICITY
!		PSTATE(08,ISET)		CCF						PARTICLE CUMMIL
!		PSTATE(09,ISET)		OFFSHERE					PARTICLE OFFSET FROM SPHERE FOR PARTICLE RELAXATION TIME
!
!		.............................................................
!
!		RSTATE(01,ISET)		KXX
!		RSTATE(02,ISET)		KYY
!		RSTATE(03,ISET)		KZZ
!		RSTATE(04,ISET)		THX
!		RSTATE(05,ISET)		THY
!		RSTATE(06,ISET)		THZ
!		RSTATE(07,ISET)		MOMX						PARTICLE MOMENT OF INERTIA IN X AXIS ( )
!		RSTATE(08,ISET)		MOMY						PARTICLE MOMENT OF INERTIA IN Y AXIS ( )
!		RSTATE(09,ISET)		MOMZ						PARTICLE MOMENT OF INERTIA IN Z AXIS ( )
!
!		..............................................................
!
!		FLUID(1,JPART)			Ug							FLUID VELOCITY AT PARTICLE POSITION (m/s)
!		FLUID(2,JPART)			Vg							FLUID VELOCITY AT PARTICLE POSITION (m/s)
!		FLUID(3,JPART)			Wg							FLUID VELOCITY AT PARTICLE POSITION (m/s)
!		FLUID(4-12,JPART)		dU/dx						FLUID ROTATION AT PARTICLE POSITION ( /s)
!		FLUID(13,JPART)		GAMsgs					FLUID SGS KINETIC ENERGY AT PARTICLE POSITION (??)

!

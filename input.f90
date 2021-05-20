SUBROUTINE INPUT(NITER)
!--------------------------------------------------------------------------!
!   Authors: W.P. Jones w.jones@ic.ac.uk - Thu Jul 25 11:30:32 BST 2002!   !
!            F.di Mare f.di.mare@ic.ac.uk                                  !
!   No part of this code can be modified without                           !
!   written permission from the Authors. Thu Jul 25 11:30:32 BST 2002      !
!                                                                          !
!   MODIFICATION:                                                          !
!   DR D.O. NJOBUENWU d.o.njobuenwu@leeds.ac.uk Mon Sept 10 13:10 BST 2012 !
!   CONVERTED TO .f90 and inclusion of MODULE                              !
!--------------------------------------------------------------------------!
!   SUBROUTINE INPUT READS IN SEVERAL VARIABLES AND PROBLEM CONSTANTS      !
!   THAT MAY BE REQUIRED TO BE CHANGED FREQUENTLY BY THE USER              !
!--------------------------------------------------------------------------!
!---------!
! MODULES !
!---------!
USE FDIM_MODULE
USE FBCFLAG_MODULE
USE FGRID_MODULE
USE FVAR_MODULE
USE FV02_MODULE
USE FV03_MODULE
USE FV04_MODULE
USE FV05_MODULE
!USE FV06_MODULE
USE FSTAT_MODULE
USE FPROP_MODULE
USE FCHEM_MODULE
USE FTURB_MODULE
USE FMODEL_MODULE
USE FFILE_MODULE
USE FUTILITY_MODULE
USE FTIME_MODULE
USE FLOGICAL_MODULE
USE FCPUTIME_MODULE
USE FINFLOW_MODULE
USE FANIMATE_MODULE
USE FCOURANT_MODULE
USE FSCALARS_MODULE
USE FBEND_MODULE
USE PART_SGS_MODULE
!USE PARTICLE
USE PART_NUM_MODULE
USE PART_COUNT_MODULE
USE PART_INITIAL_MODULE
USE PART_LOGICAL_MODULE
USE PART_STAT_MODULE
USE PART_VECTOR_MODULE
USE PART_TIME_MODULE
USE PART_PROP_MODULE
USE PART_EXTRA_MODULE
USE PART_SGS_MODULE
USE GEO_SIZE_MODULE
USE FROTAT_MODULE
USE PART_PWI_MODULE
USE PART_PDF_MODULE
USE PART_CONCEN_MODULE
USE PART_FILENAME_MODULE
USE PART_TWOWAY_MODULE
USE PART_DEPOSIT_MODULE
USE PART_FLAG_MODULE
USE PART_PSTATE_MODULE
USE PART_RSTATE_MODULE
USE PART_COLLISION_MODULE
!
IMPLICIT NONE
!-----------------------!
! VARIABLE DECLARATIONS !
!-----------------------!
INTEGER :: NITER
INTEGER :: INB,NOV
INTEGER :: I,J,K
REAL FAC
LOGICAL	::	EXIST
CHARACTER (LEN=120) :: INSTFILE
!------------------------!
! FORMATS FOR READING IN !
!------------------------!
1000 FORMAT (A80)
2100 FORMAT(1P8E15.7)
!-----------!
! OPEN FILE !
!-----------!
OPEN (IMIN,FILE='./INPUT/input.dat',STATUS='OLD')
!-----------!
!  READ IN  !
!-----------!
READ (IMIN,*)
READ (IMIN,*)
READ (IMIN,*)
READ (IMIN,1000) GEOFILE
READ (IMIN,1000) RESFILE
READ (IMIN,1000) STATFILE
READ (IMIN,1000) WALLFILE
READ (IMIN,1000) INFOFILE
READ (IMIN,*)
READ (IMIN,*)
READ (IMIN,*)
READ (IMIN,*) FORM
READ (IMIN,*) READ,WRITE
READ (IMIN,*) PSCRN,SCRN
READ (IMIN,*) NEWCOURANT
READ (IMIN,*)
READ (IMIN,*)
READ (IMIN,*)
READ (IMIN,*) LSTEP,SAVESTEP,INFO
READ (IMIN,*) TLAST
READ (IMIN,*) TEXEC
READ (IMIN,*)
READ (IMIN,*)
READ (IMIN,*)
DO NOV=1,NMAX
  READ (IMIN,*) TVD(NOV),ICYCL(NOV),TOL(NOV),NVORDR(NOV)
END DO
READ (IMIN,*)
READ (IMIN,*)
READ (IMIN,*)
READ (IMIN,*) DYN_STRESS,DYN_SFLUX,SIMIL
READ (IMIN,*) DYN_RESTART_STRESS,DYN_RESTART_SFLUX
READ (IMIN,*) CS,CXI
READ (IMIN,*)
READ (IMIN,*)
READ (IMIN,*)
READ (IMIN,*) TURBSTAT,TURBREAD
READ (IMIN,*) FSAMPLESTEP
READ (IMIN,*) FLOWSTAT,FLOWSTATREAD
READ (IMIN,*)
READ (IMIN,*)
READ (IMIN,*)
READ (IMIN,*) ALPHA
READ (IMIN,*) ALPHAT,NITER
READ (IMIN,*) IREF,JREF,KREF
READ (IMIN,*)
READ (IMIN,*)
READ (IMIN,*)
READ (IMIN,*) PERIODIC,CHANNELTYPE
READ (IMIN,*) CONTOUR_MOVIE
READ (IMIN,*)
READ (IMIN,*)
READ (IMIN,*)
READ (IMIN,*) BURN,NFCHEM
READ (IMIN,*)
READ (IMIN,*)
READ (IMIN,*)
READ (IMIN,*) CFLMIN,CFLMAX,CFL_OK
READ (IMIN,*) RHORELAX
READ (IMIN,*) CONV_OUT
READ (IMIN,*)
READ (IMIN,*)
READ (IMIN,*)
READ (IMIN,*) DIGIT_TURB,START_TURB,FORM_TURB
READ (IMIN,*) BEND
READ (IMIN,'(A120)') PATH_TURB
READ (IMIN,*)
READ (IMIN,*)
READ (IMIN,*)
!--------------------------!
! END OF LES INPUT READING !
!--------------------------!
!--------------------------!
!   CALL FLUID_ALLOCATE    !
!--------------------------!
CALL FLUID_ALLOCATE(10)
!------------------------------------------------------------------------------------------------------------------!
!    	PARTICLE (1) LOGICALS                                                                                      !
!------------------------------------------------------------------------------------------------------------------!
READ (IMIN,*) LPT
IF(LPT)THEN
   READ (IMIN,*) LPTREAD
   READ (IMIN,*) POINTX,POINTY,POINTZ
   READ (IMIN,*) IDEBURG,W2S
   READ (IMIN,*)!---------------------------------------------------------------------------------------------------
   READ (IMIN,*)!									PARTICLE  (2) INPUT PARAMETERS (NUMBERS)                       !
   READ (IMIN,*)!---------------------------------------------------------------------------------------------------
   READ (IMIN,*) NPSET, NPART
   READ (IMIN,*) INJECT_TYPE
   READ (IMIN,*) NPAR_STEP,INJECT_STEP
   READ (IMIN,*) MAXSTEP
   READ (IMIN,*) IBEGIN
   READ (IMIN,*)!---------------------------------------------------------------------------------------------------
   READ (IMIN,*)!									PARTICLE  (3) INITAL POSITION (REAL) IF POINT SOURCE IS TRUE   !
   READ (IMIN,*)!---------------------------------------------------------------------------------------------------
   !********************
   CALL PART_ALLOCATE(10)
   !********************
   READ(IMIN,*)(PPOINTX(I),I=1,NPSET)
   READ(IMIN,*)(PPOINTY(I),I=1,NPSET)
   READ(IMIN,*)(PPOINTZ(I),I=1,NPSET)
   READ (IMIN,*)!---------------------------------------------------------------------------------------------------
   READ (IMIN,*)!									PARTICLE  (4) INPUT PARAMETERS (FACTORS OR MULTIPLES OF STEPS) !
   READ (IMIN,*)!---------------------------------------------------------------------------------------------------
   READ (IMIN,*) LPT_STEP
   READ (IMIN,*) ROTTIMFAC
   READ (IMIN,*) ISAVE
   READ (IMIN,*) DEBUGSKIP
   READ (IMIN,*)!---------------------------------------------------------------------------------------------------
   READ (IMIN,*)!									PARTICLE  (5) FLAGS                                            !
   READ (IMIN,*)!---------------------------------------------------------------------------------------------------
   READ (IMIN,*) PARTSIZING !(SPECIFY PST OR DIAM, TO COMPUTE THE OTHER) [1 = SPECIFY DIAMETER, 2 =  SPECIFY STOKES NUMBER,PST]
   READ (IMIN,*) VEPWHICH
   READ (IMIN,*) OMEGAWHICH
   READ (IMIN,*) TYPE_SGS_VEL
   READ (IMIN,*) DRAGFLAG
   READ (IMIN,*) PARTCCF
   READ (IMIN,*)!----------------------------------------------------------------------------------------------------------------------------!
   READ (IMIN,*)!	PARTICLE  (6) FLAGS SET PARTICLE FORCES (DRAG, GRAVITY, SHEAR LIFT, BROWNIAN, PROFILE, ADDED MASS and PRESSURE GRADIENT )!
   READ (IMIN,*)!----------------------------------------------------------------------------------------------------------------------------!
   DO I = 1, FORCECOUNT
      READ(IMIN,*) LFORCE(I)
   END DO
   READ (IMIN,*)!---------------------------------------------------------------------------------------------------
   READ (IMIN,*)!									PARTICLE  (7) INPUT PARAMETERS                                 !
   READ (IMIN,*)!---------------------------------------------------------------------------------------------------
   READ (IMIN,*) (LAMBDA_SET(I),		I = 1, NPSET)
   READ (IMIN,*) (PST_SET(I),			I = 1, NPSET)
   READ (IMIN,*) (DIAP_SET(I),		I = 1, NPSET)
   READ (IMIN,*) (DENP_SET(I),		I = 1, NPSET)
   READ (IMIN,*) (TIMEFAC(I), 		I = 1, NPSET)
   READ (IMIN,*)!---------------------------------------------------------------------------------------------------
   READ (IMIN,*)!									PARTICLE  (6) PROPERTIES                                       !
   READ (IMIN,*)!---------------------------------------------------------------------------------------------------
   READ (IMIN,*) EA10, EA20, EA30
   READ (IMIN,*) RAN_ANGLE
   READ (IMIN,*) TEMPERATURE
   READ (IMIN,*) GBX,GBY,GBZ
   READ (IMIN,*)!---------------------------------------------------------------------------------------------------
   READ (IMIN,*)!									PARTICLE  (7) STATISTICS                                       !
   READ (IMIN,*)!---------------------------------------------------------------------------------------------------
   READ (IMIN,*) SAMPLE
   READ (IMIN,*) PART_STATIST,PART_STAT_READ
   READ (IMIN,*) PSTAT_STEP, PSTAT_TOTAL
   READ (IMIN,*) PARTCONC,CONCSTEP
   READ (IMIN,*) PDF
   READ (IMIN,*)!---------------------------------------------------------------------------------------------------
   READ (IMIN,*)!									PARTICLE  (8)  FILENAME                                        !
   READ (IMIN,*)!---------------------------------------------------------------------------------------------------
   READ (IMIN,10) PART_REST_TRAJ_FILE
   READ (IMIN,10) PART_REST_STAT_FILE
   READ (IMIN,10) PART_STAT_FILE
   READ (IMIN,10) PART_POSIT_FILE
   READ (IMIN,10) PART_CONC_FILE
   READ (IMIN,10) PART_PDF_FILE
   READ (IMIN,*)!---------------------------------------------------------------------------------------------------
   READ (IMIN,*)!									PARTICLE  (9)	 MOVIE                                         !
   READ (IMIN,*)!---------------------------------------------------------------------------------------------------
   READ (IMIN,*) PART_MOVIE
   READ (IMIN,*) MOVIE_WHERE
   READ (IMIN,*) MOVIE_STEP
   READ (IMIN,*) MOVIE_END
   READ (IMIN,*)!---------------------------------------------------------------------------------------------------
   READ (IMIN,*)!									PARTICLE  (14)	TWO WAY COUPLING                               !
   READ (IMIN,*)!---------------------------------------------------------------------------------------------------
   READ (IMIN,*) TWOWAY,SPREADER
   READ (IMIN,*) IPARCEL
   READ (IMIN,*) INJECT
   READ (IMIN,*) VOLFRAC
   READ (IMIN,*) MASFRAC
   READ (IMIN,*)!---------------------------------------------------------------------------------------------------
   READ (IMIN,*)!									PARTICLE  (15)	DEPOSITION                                     !
   READ (IMIN,*)!---------------------------------------------------------------------------------------------------
   READ (IMIN,*) DEPOSITION,DEPOSIT_TRAJ
   READ (IMIN,*) FORCE_RATIO,FORCE_RATIO_READ
   READ (IMIN,*) DSLAB
   READ (IMIN,*)!---------------------------------------------------------------------------------------------------
   READ (IMIN,*)!									PARTICLE  (16)	PARTICLE INTERACTION                           !
   READ (IMIN,*)!---------------------------------------------------------------------------------------------------
   READ (IMIN,*) PART_INTERACT,PART_INTERACT_READ
   READ (IMIN,*) LPX,LPY,LPZ
   READ (IMIN,*) POSTCOLLISION
   !--------------------------!
   !   CALL PART_ALLOCATE     !
   !--------------------------!
   CALL PART_ALLOCATE(15)
   !=== FORMAT 10 ===!
   10  FORMAT(A80)
END IF
CLOSE (IMIN,STATUS='KEEP')
!--------------------------!
!   END OF READING IN      !
!--------------------------!
!
!=== CHECK TO SEE IF WE'RE PRINTING TO THE SCREEN ===!
!
IF(PSCRN.AND.SCRN.NE.6) THEN
   OPEN  (SCRN,FILE='SCREEN',STATUS='UNKNOWN',FORM='FORMATTED',ACTION='WRITE')
   REWIND SCRN ! REWIND TO START
END IF
!
OPEN (MOUT,FILE=INFOFILE,  STATUS='UNKNOWN')
IF(FORM)THEN
   OPEN (LIN,FILE=GEOFILE,STATUS='OLD')
   READ (LIN,*,END=700) LP1,MP1,NP1
ELSE
   OPEN (LIN,FILE=GEOFILE,FORM='UNFORMATTED',STATUS='OLD')
   READ (LIN,END=700) LP1,MP1,NP1
END IF
!
!=== CHECK THAT LP1 MP1 AND NP1 HAVE ALLOWABLE VALUES ===!
!
IF(LP1 > IMAX .OR. MP1 > JMAX .OR. NP1 > KMAX) THEN
   WRITE(6,*) '!!! INSUFFICIENT STORAGE : INCREASE IMAX,JMAX AND/OR KMAX !!!'
   STOP
END IF
!----------------------------------!
!    INITIALISE BOUNDARY FLAGS     !
!----------------------------------!
DO K=1,NP1
  DO J=1,MP1
    IBN(J,K) = -1000
    IBS(J,K) = -1000
  END DO
END DO
DO K=1,NP1
  DO I=1,LP1
    IBE(I,K) = -1000
    IBW(I,K) = -1000
  END DO
END DO
DO J=1,MP1
  DO I=1,LP1
    IBR(I,J) = -1000
    IBL(I,J) = -1000
  END DO
END DO

!---------------------!
!    READ IN GRID     !
!---------------------!
IF (FORM) THEN
  READ (LIN,*,END=700) DH,UBULK,UTAU,REB,RET,RHOAIR,RHOJET,VISCO
  READ (LIN,*,END=700) RLEN,RTIME
  READ (LIN,*,END=600) (((X(I+(J-1)*LP1+(K-1)*LP1*MP1),  &
      I=1,LP1),J=1,MP1),K=1,NP1)
  READ (LIN,*,END=600) (((Y(I+(J-1)*LP1+(K-1)*LP1*MP1),  &
      I=1,LP1),J=1,MP1),K=1,NP1)
  READ (LIN,*,END=600) (((Z(I+(J-1)*LP1+(K-1)*LP1*MP1),  &
      I=1,LP1),J=1,MP1),K=1,NP1)
  READ(LIN,*) ((IBN(J,K),J=1,MP1),K=1,NP1)
  READ(LIN,*) ((IBS(J,K),J=1,MP1),K=1,NP1)
  READ(LIN,*) ((IBE(I,K),I=1,LP1),K=1,NP1)
  READ(LIN,*) ((IBW(I,K),I=1,LP1),K=1,NP1)
  READ(LIN,*) ((IBR(I,J),I=1,LP1),J=1,MP1)
  READ(LIN,*) ((IBL(I,J),I=1,LP1),J=1,MP1)

  DO NOV=1,NMAX
    READ(LIN,*) ((FNTH(NOV,J,K,1),J=1,MP1),K=1,NP1)
    READ(LIN,*) ((FSTH(NOV,J,K,1),J=1,MP1),K=1,NP1)
    READ(LIN,*) ((FEST(NOV,I,K,1),I=1,LP1),K=1,NP1)
    READ(LIN,*) ((FWST(NOV,I,K,1),I=1,LP1),K=1,NP1)
    READ(LIN,*) ((FRHT(NOV,I,J,1),I=1,LP1),J=1,MP1)
    READ(LIN,*) ((FLFT(NOV,I,J,1),I=1,LP1),J=1,MP1)
    READ(LIN,*) FINIT(NOV,1)
  END DO
  CLOSE (LIN,STATUS='KEEP')
ELSE
  READ (LIN,END=700) DH,UBULK,UTAU,REB,RET,RHOAIR,RHOJET,VISCO
  READ (LIN,END=700) RLEN,RTIME
  READ (LIN,END=600) (((X(I+(J-1)*LP1+(K-1)*LP1*MP1),I=1,LP1),  &
      J=1,MP1),K=1,NP1)
  READ (LIN,END=600) (((Y(I+(J-1)*LP1+(K-1)*LP1*MP1),I=1,LP1),  &
      J=1,MP1),K=1,NP1)
  READ (LIN,END=600) (((Z(I+(J-1)*LP1+(K-1)*LP1*MP1),I=1,LP1),  &
      J=1,MP1),K=1,NP1)
  READ (LIN,END=700) ( (IBN(J,K),J=1,MP1), K=1,NP1)
  READ (LIN,END=700) ( (IBS(J,K),J=1,MP1), K=1,NP1)
  READ (LIN,END=700) ( (IBE(I,K),I=1,LP1), K=1,NP1)
  READ (LIN,END=700) ( (IBW(I,K),I=1,LP1), K=1,NP1)
  READ (LIN,END=700) ( (IBR(I,J),I=1,LP1), J=1,MP1)
  READ (LIN,END=700) ( (IBL(I,J),I=1,LP1), J=1,MP1)
  READ (LIN,END=700) ((FINIT(NOV,INB),INB=MEAN,DEVIAT),NOV=1,NF)
  DO INB=1,2
    READ (LIN,END=700) (((FNTH(NOV,J,K,INB),NOV=1,NF),J=1,MP1),K=1,NP1)
    READ (LIN,END=700) (((FSTH(NOV,J,K,INB),NOV=1,NF),J=1,MP1),K=1,NP1)
    READ (LIN,END=700) (((FEST(NOV,I,K,INB),NOV=1,NF),I=1,LP1),K=1,NP1)
    READ (LIN,END=700) (((FWST(NOV,I,K,INB),NOV=1,NF),I=1,LP1),K=1,NP1)
    READ (LIN,END=700) (((FRHT(NOV,I,J,INB),NOV=1,NF),I=1,LP1),J=1,MP1)
    READ (LIN,END=700) (((FLFT(NOV,I,J,INB),NOV=1,NF),I=1,LP1),J=1,MP1)
  END DO
  CLOSE (LIN,STATUS='KEEP')
END IF
!--------------------------------------------!
! USE AVERAGE WALL STRESS IF STATS COLLECTED !
!--------------------------------------------!
IF(TURBREAD) THEN
   DO J = 1, MP1
      DO K = 1, NP1
         IBN(J,K) = -42
         IBS(J,K) = -42
      END DO
   END DO
ELSE
   DO J = 1, MP1
      DO K = 1, NP1
         IBN(J,K) = -41
         IBS(J,K) = -41
      END DO
   END DO
END IF
!
WRITE(SCRN,*)'DH FROM MESH = ',DH ! PRINT DH FROM MESH
!
!--------------------------------------------!
!   CHECK UU_RMS EXISTS, IF SO READ IT IN    !
!--------------------------------------------!
!
EXIST=.TRUE.
WRITE(INSTFILE,'(A)')'./INPUT/UU_RMS.dat'
INQUIRE(FILE=INSTFILE,EXIST=EXIST)
IF(EXIST) THEN
    WRITE(SCRN,'(A)')TRIM(INSTFILE)// 'FILE EXIST'
    OPEN(1222,FILE=INSTFILE,ACTION='READ')
    READ(1222,*) ((FLFT(3,I,J,2), I = 1,LP1),J = 1,MP1)
    CLOSE(1222)
ELSE
    WRITE(SCRN,'(A)')'UU_RMS.dat DOES NOT EXIT'
END IF
!--------------------------------------------!
!   CHECK UU_RMS EXISTS, IF SO READ IT IN    !
!--------------------------------------------!
TNPART = NPSET*NPART ! TOTAL NUMBER OF PARTICLES
FAC=UTAU*RHOAIR/VISCO ! FACTOR FOR NON DIMENSIONALISING
!--------------------------------------------!
!              CHECK FOR BENDING             !
!--------------------------------------------!
IF(BEND)THEN
  OPEN(101,FILE="./INPUT/BENDINPUT")
  REWIND(101)
  READ(101,*)BENDTYPE
  READ(101,*)NODEH,NODEV,NODEB,NBIN,NBOUT
  READ(101,*)RDIV,RBEND,ROVAD
  CLOSE(101)
END IF
!--------------------------------------------!
!          CHECK FOR CONTOUR MOVIE           !
!--------------------------------------------!
IF(CONTOUR_MOVIE)THEN
   OPEN(584,FILE='./INPUT/animation_input.dat',ACTION='READ')
   READ(584,*)!HEADER
   READ(584,*)FLUID_MOVIE
   READ(584,*)SIZE_ANIMATE
   READ(584,*)STEP_ANIMATE,END_ANIMATE
   READ(584,*)!HEADER FOR FIXED CONTOUR
   READ(584,*)LCONTOUR
   READ(584,*)SCONTOUR
   !*********************
   CALL FLUID_ALLOCATE(11)
   !*********************
   READ(584,*)(POST_ANIMATE(I),I=1,SIZE_ANIMATE)
   READ(584,*)(XCONTOUR(I),I=1,SCONTOUR)
   CLOSE(584)
   !WRITE(SCRN,*)'LCONTOUR=',LCONTOUR,' SCONTOUR=',SCONTOUR
   !DO I=1,SCONTOUR
      !WRITE(SCRN,*)'I=',I,' XCONTOUR(I)=',XCONTOUR(I)
   !END DO
   IF(LCONTOUR)THEN
      !*********************
      CALL FLUID_ALLOCATE(12)
      !*********************
      DO K=1,SCONTOUR
         DO I=1,LP1-1
            IF(XCONTOUR(K).GE.X(I)*FAC .AND. XCONTOUR(K).LT. X(I+1)*FAC)THEN
               IJK_CONTOUR(K)=I
               !WRITE(SCRN,*)'K=',K,' I=',I,' IJK_CONTOUR(K)=',IJK_CONTOUR(K)
               EXIT
            END IF
         END DO
      END DO
      !WRITE(SCRN,*)'LCONTOUR=',LCONTOUR,' SCONTOUR=',SCONTOUR
      IJK_CONTOUR(SCONTOUR)=LP1/2+1
      XCONTOUR(SCONTOUR)=X(IJK_CONTOUR(SCONTOUR))*FAC
      !WRITE(SCRN,*)'IJK_CONTOUR(SCONTOUR)=',IJK_CONTOUR(SCONTOUR),' XCONTOUR(SCONTOUR)=',XCONTOUR(SCONTOUR)
   END IF
END IF
!------------------------------------------------------------------!
!          CHECK BOUNDARY CONDITION FLAGS FOR CONSISTENCY          !
!------------------------------------------------------------------!
DO  K = 2,NP1-1
  DO  J = 2,MP1-1
    IF ( IBS(J,K) == 0 .AND. (IBS(J,K) /= IBN(J,K)) ) THEN
      WRITE(*,*)' --- INCONSISTENT BOUNDARY FLAGS --- '
      WRITE(*,*) 'J=',J,' K=',K
      WRITE(*,*) 'IBS(J,K)=',IBS(J,K),' IBN(J,K)=',IBN(J,K)
      STOP
    END IF
  END DO
END DO

DO  K = 2,NP1-1
  DO  I = 2,LP1-1
    IF ( IBW(I,K) == 0 .AND. (IBW(I,K) /= IBE(I,K)) ) THEN
      WRITE(*,*)' --- INCONSISTENT BOUNDARY FLAGS --- '
      WRITE(*,*) 'I=',I,' K=',K
      WRITE(*,*) 'IBW(I,K)=',IBW(I,K),' IBE(I,K)=',IBE(I,K)
      STOP
    END IF
  END DO
END DO

DO  J = 2,MP1-1
  DO  I = 2,LP1-1
    IF ( IBL(I,J) == 0 .AND. (IBL(I,J) /= IBR(I,J)) ) THEN
      WRITE(*,*)' --- INCONSISTENT BOUNDARY FLAGS --- '
      WRITE(*,*) 'I=',I,' J=',J
      WRITE(*,*) 'IBL(I,J)=',IBL(I,J),' IBR(I,J)=',IBR(I,J)
      STOP
    END IF
  END DO
END DO

RETURN

600 CONTINUE
WRITE(*,*) ' --- END-OF-FILE REACHED : GEOMETRY DATA FILE --- '
STOP

700 CONTINUE
WRITE(*,*) ' --- END-OF-FILE REACHED : BOUNDARY CONDITION DATA FILE --- '
STOP

800 CONTINUE
WRITE(*,*) ' --- END-OF-FILE REACHED IN COMBUSTION DATA FILE --- '
STOP

END SUBROUTINE INPUT

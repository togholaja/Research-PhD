!
!=======================================================================
!
!   NUMERICAL INTEGRATION OF PARTICLE EQUATION OF MOTION
!
!=======================================================================
!
SUBROUTINE PART_BBO_YIN

USE  FDIM_MODULE
USE  FFILE_MODULE
USE  FPROP_MODULE
USE  FGRID_MODULE
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
USE  FUTILITY_MODULE
USE  PART_TWOWAY_MODULE
USE  FTURB_MODULE
use  PART_FLAG_MODULE
USE  PART_PSTATE_MODULE
USE  PART_RSTATE_MODULE

!
IMPLICIT NONE
!
INTEGER	::	I,IJK
REAL	::	LAMBDA,MASS
DOUBLE PRECISION	::	H1,H2,H6
DOUBLE PRECISION	::	XX(N6),DYDX(N6),DYM(N6),DYT(N6)
!
!::::::::::::::::::::::::::::::::::
MASS = PART(10,JPART)
LAMBDA = PSTATE(01,ISET)
!::::::::::::::::::::::::::::::::::
!
H1 = PTIM
H2 = H1*0.5D0
H6 = H1/6.0D0
!
IF(.NOT.ALLOCATED (YT))CALL ALLOCATE_PART_VECTOR_RK4

DO I = 1, N6
   XX(I) = PART0(I,JPART)
   YT(I) = XX(I)
END DO
!
IF(LAMBDA /= 1.0)THEN
   CALL PART_MATRIX
END IF
!
!APPLY FOURTH ORDER RUNGE KUTTA RK4 FOR POSITION X AND VELOCITY U
!---------------------------
CALL PART_FORCE_YIN(YT,DYDX)
!---------------------------
!
UREL0 = PIN(1,JPART)
CD0  	= PIN(3,JPART)
TAUP0 = PIN(4,JPART)
DO I = 1, N6
   YT(I) = XX(I) + H2*DYDX(I)
END DO
!--------------------------
CALL PART_FORCE_YIN(YT,DYT)
!--------------------------
DO I = 1, N6
   YT(I) = XX(I) + H2*DYT(I)
END DO
!--------------------------
CALL PART_FORCE_YIN(YT,DYM)
!--------------------------
DO I = 1, N6
   YT(I) = XX(I) + H1*DYM(I)
   DYM(I) = DYT(I) + DYM(I)
END DO
!--------------------------
CALL PART_FORCE_YIN(YT,DYT)
!--------------------------
DO I = 1, N6
   YT(I) = XX(I) + H6*(DYDX(I) + DYT(I) + 2.0D0*DYM(I))
END DO
!
!UPDATE PARTICLE VELOCITIES AND POSITIONS
DO I = 1, N6
   PART(I,JPART) = YT(I)
END DO
!
!INCLUDE THE EFFECT OF SGS ENERGY ON PARTICLE DISPERSION
IF(TYPE_SGS_VEL == 1)THEN
   !----------------
   CALL PART_SGS_VEL
   !----------------
   PART(1,JPART) = PART(1,JPART) + STOCA1*PTIM
   PART(2,JPART) = PART(2,JPART) + STOCA2*PTIM
   PART(3,JPART) = PART(3,JPART) + STOCA3*PTIM

   PART(4,JPART) = PART(4,JPART) + STOCA1
   PART(5,JPART) = PART(5,JPART) + STOCA2
   PART(6,JPART) = PART(6,JPART) + STOCA3
   !
END IF
IF(ALLOCATED(YY))DEALLOCATE(YY)
!
IF(TWOWAY .AND. LPTSTEP > 1)THEN
   FORCEU = MASS*(PART(4,JPART) - PART0(4,JPART))/PTIM
   FORCEV = MASS*(PART(5,JPART) - PART0(5,JPART))/PTIM
   FORCEW = MASS*(PART(6,JPART) - PART0(6,JPART))/PTIM

   FORCE(1,JPART) = FORCEU*IPARCEL
   FORCE(2,JPART) = FORCEV*IPARCEL
   FORCE(3,JPART) = FORCEW*IPARCEL

   IJK = IPS(1,JPART) + JO(IPS(2,JPART)) + KO(IPS(3,JPART))
   FX(IJK)  = FX (IJK) - FORCE(1,JPART)/ABS(AJC(IJK))
   FY(IJK)  = FY (IJK) - FORCE(2,JPART)/ABS(AJC(IJK))
   FZ(IJK)  = FZ (IJK) - FORCE(3,JPART)/ABS(AJC(IJK))

END IF

PART0(10,JPART) = PTIM

RETURN
END SUBROUTINE PART_BBO_YIN
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!HYDRODYNAMIC FORCES:
!DRAG, LIFT, VIRTUAL MASS, PRESSURE GRADIENT, CGRAVITY, BROWNIAN MOTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
SUBROUTINE PART_FORCE_YIN(SOLU,DERIVAT)
!
USE  FDIM_MODULE
USE  FPROP_MODULE
USE  PART_COUNT_MODULE
USE  PART_LOGICAL_MODULE
USE  PART_VECTOR_MODULE
USE  PART_TIME_MODULE
USE  PART_PROP_MODULE
USE  PART_EXTRA_MODULE
USE  PART_NUM_MODULE
USE  FTURB_MODULE
USE  FFILE_MODULE
USE  FUTILITY_MODULE
USE  PART_INITIAL_MODULE
use  PART_FLAG_MODULE
USE  PART_PSTATE_MODULE
USE  PART_RSTATE_MODULE
!
IMPLICIT NONE
!
INTEGER I
REAL,DIMENSION(N6)::SOLU,DERIVAT
REAL,DIMENSION(ND)::VP0,VG0,USLIP
REAL,DIMENSION(ND)::TRANSF,CON_PG,CON_VM
REAL RHOFP,FRE,TAOP,GRAV
REAL RK1,RK2
REAL FD,FL,CL0,CL1,CL2,AFAI
REAL EPACD
REAL CC1,CC2,CC3,ZC1,ZC2,ZC3
REAL BB1,BB2,BB3,BB4
REAL CD,REP,DENOM
REAL MFROT
REAL::SAFC
REAL RANDX,RANDY,RANDZ
real TEMP1,SIGMA1,DXX
REAL,EXTERNAL:: GASDEV,RAN1
REAL FAC
REAL	::	FORCE1M,FORCE2M,FORCE3M,FORCE4M,FORCE5M,FORCE6M,FORCE7M
REAL	::	RADA,ESDV,LAMBDA,MASS,PSN,CCF,RHOD,PST,OFFSPHERE,RHOPF
DOUBLE PRECISION	::	EPS1,EPS2,EPS3,EPS0
REAL, allocatable :: saffman(:)


!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!:::::::::::::::		VARIABLE DEFINITIONS		::::::::::::::::::::::::::
!
!	NAME			UNIT			DEFINITION
! 	AFAI			rad			INCIDENCE ANGLE BETWEEN Z'-AXIS AND RELATIVE VELOCITY VECTOR:
!  CD				-				DRAG COEFFICIENT
!	EPACD 		m				EQUAL PROJECTED AREA CIRCLE DIAMETER
!  FRE			-				SCHILLER AND NAUMANN DRAG FACTOR
!	FLFD			-				RATIO OF LIFT TO DRAG FORCES
!  PS_NORM		-				PARTICLE SPHERICITY (NORMAL) =s/S = SPHVOL/PASA
!  REP			-				PARTICLE REYNOLDS NUMBER BASED ON EQUIVALENT SPHERE DIAMETER
!	RHOFP		-					FLUID-TO-PARTICLE DENSITY RATIO=RHOAIR/RHOP
!	s=SPHVOL		m2				SURFACE AREA OF A SPHERE HAVING THE SAME VOLUME AS THE NON-SPHERICAL PARTICLE
! 	S=PASA		m2				ACTUAL SURFACE AREA OF THE NON-SPHERICAL PARTICLE, PASA
!	SEFFD			m2				PARTICLE AREA NORMAL TO THE DIRECTION OF THE DRAG FORCE
!	SEFFL			m2				PARTICLE AREA NORMAL TO THE DIRECTION OF THE PROFILE LIFT FORCE
!
!:::::::::::::::::::::	FLAGS	::::::::::::::::::::::::::::::::::::::::::::::::::::::
!	FLAG				VALUE				DESCRIPTION
!  DRAGFLAG(0)		0, DEFAULT		GANSER DRAG COEFFICIENT
!	DRAGFLAG(1)		1					SCHILLER AND NAUMANN DRAG COEFFICIENT

!	LFORCE(1)		1 OR 0	1 = DRAG;						0 = NO DRAG						FORCE
!	LFORCE(2)		1 OR 0	1 = GRAVITY & BOUYANCY;		0 = NO GRAVITY & BOUYANCY		FORCE
!	LFORCE(3)		1 OR 0	1 = SHEAR LIFT;				0 = NO SHEAR LIFT 				FORCE
!	LFORCE(4)		1 OR 0	1 = BROWNIAN MOTION;			0 = NO BROWNIAN MOTION			FORCC
!	LFORCE(5)		1 OR 0	1 = PROFILE LIFT;				0 = NO PROFILE LIFT				FORCE (FOR NON-SPHERES ONLY)
!	LFORCE(6)		1 OR 0	1 = ADDED MASS;				0 = NO PRESSURE GRADIENT		FORCE
!	LFORCE(7)		1 OR 0	1 = PRESSURE GRADIENT;		0 = NO ADDED MASS					FORCE
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!:::::::::::::::::::::::::::::::
LAMBDA 	=	PSTATE(01,ISET)
PST		= 	PSTATE(04,ISET)
RHOD		= 	PSTATE(06,ISET)
PSN		= 	PSTATE(07,ISET)
CCF		= 	PSTATE(08,ISET)
OFFSPHERE= 	PSTATE(09,ISET)
RHOPF		=  RHOD/RHOAIR
RHOFP 	= 	1.0D0/RHOPF
!
ESDV     =  PART(08,JPART)
RADA 		= 	PART(09,JPART)
MASS		= 	PART(10,JPART)
!
!:::::::::::::::::::::::::::::::
!
DO I = 1, ND

   VG0(I) = FLUID(I,JPART)
   VP0(I) = SOLU(I+ND)

   !IF(JPART == 1) THEN
   !   WRITE(*,*) 'VG0: ',FLUID(I,JPART)
   !   WRITE(*,*) 'VP0: ',SOLU(I+ND)
   !END IF

END DO
!
FORCE1M = 0.0
FORCE2M = 0.0
FORCE3M = 0.0
FORCE4M = 0.0
FORCE5M = 0.0
FORCE6M = 0.0
FORCE7M = 0.0
!
!SLIP VELOCITY BETWEEN PARTICLE AND FLUID
USLIP(1) = VG0(1) - VP0(1)
USLIP(2) = VG0(2) - VP0(2)
USLIP(3) = VG0(3) - VP0(3)

!IF(JPART == 1) THEN
!  WRITE(*,*) 'USLIP1: ',USLIP(1)
!  WRITE(*,*) 'USLIP2: ',USLIP(2)
!  WRITE(*,*) 'USLIP3: ',USLIP(3)
!END IF

IF(USLIP(1) < SMALL) USLIP(1) = SMALL
IF(USLIP(2) < SMALL) USLIP(2) = SMALL
IF(USLIP(3) < SMALL) USLIP(3) = SMALL

!MAGNETITUDE OF SLIP VELOCITY
UREL = SQRT(USLIP(1)**2 + USLIP(2)**2 + USLIP(3)**2)
PIN(1,JPART) = UREL
!
!REP IS PARTICLE REYNOLDS NUMBER BASED ON RELATIVE VELOCITY AND
!EQUAL-VOLUME SPHERE DIAMETER???
REP = UREL*ESDV*RHOAIR/VISCO
PIN(2,JPART) = REP
!
IF(LAMBDA == 1.0)THEN
   AFAI = RAD90
   SEFFD = PIE*RADA**2
   SEFFL = SEFFD
ELSE
   !ORIENTATION OF THE MAJOR PARTICLE-AXIS Z' IN THE INERTIAL FRAME: EQN(4) 3RD ROW
   EPS1 = ROT(1,JPART)
   EPS2 = ROT(2,JPART)
   EPS3 = ROT(3,JPART)
   EPS0 = ROT(0,JPART)

   ZC1 = 2.0D0*(EPS3*EPS1 + EPS2*EPS0)
   ZC2 = 2.0D0*(EPS3*EPS2 - EPS1*EPS0)
   ZC3 = 1.0 - 2.0D0*(EPS1**2 + EPS2**2)
   BB1 = SQRT(ZC1**2 + ZC2**2 + ZC3**2)
   IF(BB1 /= 1.0)THEN
      ZC1 = ZC1/BB1
      ZC2 = ZC2/BB1
      ZC3 = ZC3/BB1
   END IF
   !INCIDENCE ANGLE BETWEEN Z'-AXIS AND RELATIVE VELOCITY VECTOR:
   !SEE FIG.3 ON PAGE 3492 CHEM ENG SCI 58, 2003, 3489-3498
   BB1 = ZC1*USLIP(1)
   BB2 = ZC2*USLIP(2)
   BB3 = ZC3*USLIP(3)

   BB4 = (BB1 + BB2 + BB3)/UREL

   IF(BB4 < -1.0)THEN
      BB4 = -1.0
   ELSE IF(BB4 > 1.0)THEN
      BB4 = 1.0
   END IF
   AFAI = ACOS(BB4)
   !USED BY YIN ET AL CHEM ENG SCI 2013 FOR NON-SPHERE
   SEFFD = PIE*RADA**2*SQRT(COS(AFAI)**2 + (4*LAMBDA/PIE)**2*SIN(AFAI)**2)
   SEFFL = PIE*RADA**2*SQRT(SIN(AFAI)**2 + (4*LAMBDA/PIE)**2*COS(AFAI)**2)
END IF
PIN(5,JPART) = AFAI
!----------
!DRAG FORCE
!----------
IF(LFORCE(1) == 1)THEN
   !SELECT DRAG COEFFICIENT CORELECTION
   SELECT CASE(DRAGFLAG)
      CASE(1)
         !SCHILLER AND NAUMANN DRAG COEFFICIENT
         IF(REP > 1000)THEN
            FRE = 0.44D0*REP/24.0D0
         ELSE IF(REP < 1)THEN
            FRE = 1.0
         ELSE
            FRE = 1.0 + 0.15D0*REP**0.687D0
         END IF
         !DRAG COEFFICIENT
         CD = FRE *24.0/REP
         !DYNAMIC PARTICLE RELAXATION TIME
         IF(LFORCE(7) == 1)THEN
            !INCLUDE THE EFFECT OF ADDED MASS FORCE ON PARTICLE RELAXATION TIME
            TAOP = FOURTHIRD*ESDV*(RHOPF + 0.5)*CCF/(CD*UREL)
         ELSE
            !STANDARD TAO WITHOUT ADDED MASS FORCE
            TAOP = FOURTHIRD*ESDV*(RHOPF + 0.0)*CCF/(CD*UREL)
         END IF

         FD = MASS/TAOP

         IF(JPART == 1) THEN
            !WRITE(*,*) 'FD = ', FD
         END IF

         FORCEX(1) = FD*USLIP(1)
         FORCEY(1) = FD*USLIP(2)
         FORCEZ(1) = FD*USLIP(3)
         FORCE5M = 0.0
      CASE DEFAULT
         !GANSER'S DRAG COEFFICIENT
         !THE BEST CD TILL NOW THAT IS BOTH SHAPE-FACTOR AND  ...
         !ORIENTATION-DEPENDENT: EQN(18)
         EPACD = SQRT(4.0D0*SEFFD/PIE)
         PSN = PSTATE(07,ISET)

         RK1 = EPACD/(3.0D0*ESDV) + 2.0D0/(3.0D0*SQRT(PSN))
         RK1 = 1.0/RK1
         RK2 = 1.8148D0*(-1.0*LOG10(PSN))**0.5743D0
         RK2 = 10.0**RK2

         BB1 = REP*RK1*RK2
         !DRAG COFFICIENT
         CD = RK2*((24.0/BB1)*(1.0 + 0.1118D0*(BB1**0.6567)) +  0.4305D0/(1.0 + 3305D0/BB1))
         CD = CD/CCF
         FAC = CD*UREL
         IF(LFORCE(7) == 1)THEN
            !INCLUDE THE EFFECT OF ADDED MASS FORCE ON PARTICLE RELAXATION TIME
            TAOP = FOURTHIRD*ESDV*(RHOPF + 0.5)/(CD*UREL)
         ELSE
            !STANDARD TAO WITHOUT ADDED MASS FORCE
            TAOP = FOURTHIRD*ESDV*(RHOPF + 0.0)/(CD*UREL)
         END IF
         !Drag Force calculation: see Eqn(16) of Yin et al 2003
         FD = 0.5D0*CD*RHOAIR*SEFFD*UREL
         FORCEX(1) = FD*USLIP(1)
         FORCEY(1) = FD*USLIP(2)
         FORCEZ(1) = FD*USLIP(3)
         FORCE1M = SQRT(FORCEX(1)**2 + FORCEY(1)**2 + FORCEZ(1)**2)
         IF(FORCE1M < SMALL)THEN
            FORCE1M = SMALL
         END IF
   END SELECT
ELSE IF(LFORCE(1) == 0)THEN
   !NO DRAG FORCE
    FORCEX(1) = 0.0
    FORCEY(1) = 0.0
    FORCEZ(1) = 0.0
    FORCE1M = 0.0
    TAOP = SMALL
    CD = SMALL
END IF
PIN(3,JPART) = CD
PIN(4,JPART) = TAOP
PIN(6,JPART) = SEFFD
!---------------------------
! GRAVITY AND BOUYANCY FORCE
!---------------------------
IF(LFORCE(2) == 1)THEN
   !FORCE = ma (mass * acceleration)
   GRAV = 9.81D0*(1.0 - RHOFP)*MASS
   FORCEX(2) = GBX*GRAV
   FORCEY(2) = GBY*GRAV
   FORCEZ(2) = GBZ*GRAV
   FORCE2M = SQRT(FORCEX(2)**2+FORCEY(2)**2+FORCEZ(2)**2)
ELSE IF(LFORCE(2) == 0)THEN
   !NO DRAG FORCE
   FORCEX(2) = 0.0
   FORCEY(2) = 0.0
   FORCEZ(2) = 0.0
END IF
!
!----------------
!SHEAR LIFT FORCE
!----------------
IF(LFORCE(3) == 1)THEN
   !---------------------------------------
   !OBTAIN FLUID ROTATION IN INERTIAL FRAME
   !---------------------------------------
   CALL PART_FLUID_TRANS(0)

   !MAGNITUDE OF THE FLUID ROTATION
   MFROT = SQRT(FROTX(1)*FROTX(1) + FROTX(2)*FROTX(2) + FROTX(3)*FROTX(3))
   !
   IF(.NOT.ALLOCATED (saffman))ALLOCATE (saffman(3))

	SAFC = 6.46D0*RHOAIR*RADA**2*SQRT(NU)/SQRT(MFROT)
	SAFFMAN(1) = USLIP(2)*FROTX(3) - USLIP(3)*FROTX(2)
	SAFFMAN(2) = USLIP(3)*FROTX(1) - USLIP(1)*FROTX(3)
	SAFFMAN(3) = USLIP(1)*FROTX(2) - USLIP(2)*FROTX(1)

   FORCEX(3) = (SAFC*SAFFMAN(1))
   FORCEY(3) = (SAFC*SAFFMAN(2))
   FORCEZ(3) = (SAFC*SAFFMAN(3))

   FORCE3M = SQRT(FORCEX(3)**2+FORCEY(3)**2+FORCEZ(3)**2)

ELSE IF(LFORCE(3) == 0)THEN
   !NO SHEAR-SLIP FORCE
   FORCEX(3) = 0.0
   FORCEY(3) = 0.0
   FORCEZ(3) = 0.0
END IF
!
!---------------
!BROWNIAN MOTION
!---------------
!
IF(LFORCE(4) == 1)THEN
   RANDX = GASDEV(NSEEDX3(ISET))
   RANDY = GASDEV(NSEEDY3(ISET))
   RANDZ = GASDEV(NSEEDZ3(ISET))

   !BASED ON AHMADI FOR SPHERICAL PARTICLES, Li and Ahmadi (1993)
   !SEE EQ 25 Kleinstreuer and Feng JBioEng 2013,
   !BM ADDED TO FORCE TERM IN BBO

   TEMP1 = TEMPERATURE + 273.0

   DXX = 216.0*NU*KB*TEMP1/(PIE**2*RHOAIR*ESDV**5*RHOPF**2*CCF)
   SIGMA1 = SQRT(PIE*DXX/PTIM)
   !FORCE = MASS * ACCELERATION
   FORCEX(4) = SIGMA1*RANDX*MASS
   FORCEY(4) = SIGMA1*RANDY*MASS
   FORCEZ(4) = SIGMA1*RANDZ*MASS

   FORCE4M = SQRT(FORCEX(4)**2 + FORCEY(4)**2 + FORCEZ(4)**2)
ELSE IF(LFORCE(4) == 0)THEN
   !NO BROWNIAN MOTION
   FORCEX(4) = 0.0
   FORCEY(4) = 0.0
   FORCEZ(4) = 0.0
END IF
!
!-------------------------------------------
!PROFILE LIFT FORCE FOR NON-SPHERICAL PARTICLES ONLY
!-------------------------------------------
!
IF(LFORCE(5) == 1 .AND. LAMBDA /= 1.0)THEN
   !Lift Force calculation: see Eqn(19) & (21) of Yin et al 2003
   CL0 = CD*SIN(AFAI)
   FL = 0.5D0*CL0*RHOAIR*SEFFL*COS(AFAI)
   BB1 = ZC2*USLIP(3) - ZC3*USLIP(2)
   BB2 = ZC3*USLIP(1) - ZC1*USLIP(3)
   BB3 = ZC1*USLIP(2) - ZC2*USLIP(1)
   !
   FORCEX(5) = FL*(BB2*USLIP(3) - BB3*USLIP(2))
   FORCEY(5) = FL*(BB3*USLIP(1) - BB1*USLIP(3))
   FORCEZ(5) = FL*(BB1*USLIP(2) - BB2*USLIP(1))
   !
   FORCE5M = SQRT(FORCEX(5)**2 + FORCEY(5)**2 + FORCEZ(5)**2)
   !
   IF(FORCE1M <= SMALL)THEN
      CL1 = SMALL
   ELSE
      CL1 = FORCE5M/FORCE1M
   END IF

   CL2 = SIN(AFAI)**2*ABS(COS(AFAI))/(0.65 + 40*REP**0.72)

   IF(CL1 /= 0)THEN
      IF((CL1 - CL2) /= 0)THEN
         FORCEX(5) = FORCEX(5)*(CL2/CL1)
         FORCEY(5) = FORCEY(5)*(CL2/CL1)
         FORCEZ(5) = FORCEZ(5)*(CL2/CL1)
      END IF
   END IF
   FORCE5M = SQRT(FORCEX(5)**2 + FORCEY(5)**2 + FORCEZ(5)**2)
   IF(FORCE1M < SMALL)THEN
      FLFD = SMALL
   ELSE
      FLFD = FORCE5M/FORCE1M
   END IF
ELSE IF(LFORCE(5) == 1)THEN
   !NO PROFILE LIFT FORCE
   FORCEX(5) = 0.0
   FORCEY(5) = 0.0
   FORCEZ(5) = 0.0
END IF
!
!-----------------------------------------------------
!DEFINITION OF VIRTUAL-MASS & PRESSURE-GRADIENT FORCES
!EQUATIONS (12), (13) & (14) IN YIN ET AL CHEM ENG SCI 2003
!-----------------------------------------------------
!
TransF(1)=0.0
TransF(2)=0.0
TransF(3)=0.0
!
!--------------------------------
!VIRTUAL-MASS OR ADDED-MASS FORCE
!--------------------------------
IF(LFORCE(6) == 1)THEN
   !VIRTUAL MASS FORCE IN X DIRECTION USING FLUID ROTATION IN INERTIAL FRAME
   CC1 = VP0(1)*FLUID(04,JPART)
   CC2 = VP0(2)*FLUID(05,JPART)
   CC3 = VP0(3)*FLUID(06,JPART)
   CON_VM(1) = CC1 + CC2 + CC3
   FORCEX(6) = 0.5*RHOFP*MASS*(TRANSF(1) + CON_VM(1))
   !VIRTUAL MASS FORCE IN Y DIRECTION
   CC1 = VP0(1)*FLUID(07,JPART)
   CC2 = VP0(2)*FLUID(08,JPART)
   CC3 = VP0(3)*FLUID(09,JPART)
   CON_VM(2) = CC1 + CC2 + CC3
   FORCEY(6) = 0.5*RHOFP*MASS*(TRANSF(2) + CON_VM(2))
   !VIRTUAL MASS FORCE IN Z DIRECTION
   CC1 = VP0(1)*FLUID(10,JPART)
   CC2 = VP0(2)*FLUID(11,JPART)
   CC3 = VP0(3)*FLUID(12,JPART)
   CON_VM(3) = CC1 + CC2 + CC3
   FORCEZ(6) = 0.5*RHOFP*MASS*(TRANSF(3) + CON_VM(3))
ELSE IF(LFORCE(6) == 0)THEN
   !NO VIRTUAL-MASS OR ADDED-MASS FORCE
   FORCEX(6) = 0.0
   FORCEY(6) = 0.0
   FORCEZ(6) = 0.0
END IF
!
!-----------------------
!PRESSURE-GRADIENT FORCE
!-----------------------
!
IF(LFORCE(7) == 1)THEN
   !PRESSURE GRADIENT FORCE IN X DIRECTION USING FLUID ROTATION IN INERTIAL FRAME
   CC1 = VG0(1)*FLUID(04,JPART)
   CC2 = VG0(2)*FLUID(05,JPART)
   CC3 = VG0(3)*FLUID(06,JPART)
   CON_PG(1) = CC1 + CC2 + CC3
   FORCEX(7) = RHOFP*MASS*(TRANSF(1) + CON_PG(1))
   !PRESSURE GRADIENT FORCE IN Y DIRECTION
   CC1 = VG0(1)*FLUID(07,JPART)
   CC2 = VG0(2)*FLUID(08,JPART)
   CC3 = VG0(3)*FLUID(09,JPART)
   CON_PG(2) = CC1 + CC2 + CC3
   FORCEY(7) = RHOFP*MASS*(TRANSF(2) + CON_PG(2))
   !PRESSURE GRADIENT FORCE IN Z DIRECTION
   CC1 = VG0(1)*FLUID(10,JPART)
   CC2 = VG0(2)*FLUID(11,JPART)
   CC3 = VG0(3)*FLUID(12,JPART)
   CON_PG(3) = CC1 + CC2 + CC3
   FORCEZ(7) = RHOFP*MASS*(TRANSF(3) + CON_PG(3))
   FORCE7M = SQRT(FORCEX(7)**2+FORCEY(7)**2+FORCEZ(7)**2)
ELSE
   !NO PRESSURE-GRADIENT FORCE
   FORCEX(7) = 0.0
   FORCEY(7) = 0.0
   FORCEZ(7) = 0.0
END IF


!NOTE GRAVITY IS NOT INCLUDED IN HYDRODYNAMICE FORCES HFORCE FOR PITCHING TORQUE
IF(LFORCE(6))THEN
   !HYDRODYNAMIC FORCE WITH VIRTUAL-MASS FORCE
   DENOM = (1 + 0.5D0*RHOFP)
ELSE
   DENOM = 1.0
END IF

!IF(JPART == 1) THEN
!    WRITE(*,*) 'FY1: ',FORCEY(1)
!    WRITE(*,*) 'FY2: ',FORCEY(2)
!    WRITE(*,*) 'FY3: ',FORCEY(3)
!    WRITE(*,*) 'FY4: ',FORCEY(4)
!    WRITE(*,*) 'FY5: ',FORCEY(5)
!    WRITE(*,*) 'FY6: ',FORCEY(6)
!    WRITE(*,*) 'FY7: ',FORCEY(7)
!END IF

HFORCE(1) = (FORCEX(1) + FORCEX(3) + FORCEX(5) + FORCEX(6) + FORCEX(7))/DENOM
HFORCE(2) = (FORCEY(1) + FORCEY(3) + FORCEY(5) + FORCEY(6) + FORCEY(7))/DENOM
HFORCE(3) = (FORCEZ(1) + FORCEZ(3) + FORCEZ(5) + FORCEZ(6) + FORCEZ(7))/DENOM

DERIVAT(1+ND) = (FORCEX(1) + FORCEX(2) + FORCEX(3) + FORCEX(4) + FORCEX(5) + FORCEX(6) + FORCEX(7))/(DENOM*MASS)
DERIVAT(2+ND) = (FORCEY(1) + FORCEY(2) + FORCEY(3) + FORCEY(4) + FORCEY(5) + FORCEY(6) + FORCEY(7))/(DENOM*MASS)
DERIVAT(3+ND) = (FORCEZ(1) + FORCEZ(2) + FORCEZ(3) + FORCEZ(4) + FORCEZ(5) + FORCEZ(6) + FORCEZ(7))/(DENOM*MASS)

!DERIVAT(1+ND) = (HFORCE(1) + FORCEX(2) + FORCEX(4))/MASS
!DERIVAT(2+ND) = (HFORCE(2) + FORCEY(2) + FORCEY(4))/MASS
!DERIVAT(3+ND) = (HFORCE(3) + FORCEZ(2) + FORCEZ(4))/MASS

DERIVAT(1) = SOLU(1+ND)
DERIVAT(2) = SOLU(2+ND)
DERIVAT(3) = SOLU(3+ND)

DO I = 1, FORCECOUNT
   !SELECT CASE(I)
      !CASE(1,3,5:7)
      	!HYDRODYNAMIC FORCES
         FORCEX(I) = FORCEX(I)/(DENOM*MASS)
         FORCEY(I) = FORCEY(I)/(DENOM*MASS)
         FORCEZ(I) = FORCEZ(I)/(DENOM*MASS)
!      CASE(2,4)
!      	!OTHER FORCES
!         FORCEX(I) = FORCEX(I)/(MASS)
!         FORCEY(I) = FORCEY(I)/(MASS)
!         FORCEZ(I) = FORCEZ(I)/(MASS)
!   END SELECT
END DO

IF(FORCE_RATIO .AND. LPTSTEP >= 100)THEN
   CALL PART_FORCESTAT(1)
END IF

RETURN
END SUBROUTINE PART_FORCE_YIN

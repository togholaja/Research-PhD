SUBROUTINE PART_COLLISION
!
USE FTURB_MODULE
USE FUTILITY_MODULE
USE FFILE_MODULE
USE FTIME_MODULE
USE PART_NUM_MODULE
USE PART_VECTOR_MODULE
USE PART_COLLISION_MODULE
USE PART_COUNT_MODULE
USE PART_TIME_MODULE
USE PART_PSTATE_MODULE
! COMPUTES PARTICLE-PARTICLE COLLISIONS
!
IMPLICIT NONE
!
INTEGER	::	I,II,J,K,IJKP
INTEGER	:: IP1,IP2
INTEGER	::	IPX,IPY,IPZ
INTEGER	::	JALLOC,NN,N1,N2
INTEGER :: COLLISION_SLAB(32)
REAL		:: DPAR1,DPAR2,DPAR3,DPAR13,DPAR23,DPAR33
REAL		:: XMAXCOL,YMAXCOL,ZMAXCOL,DXCOL,DYCOL,DZCOL
REAL		:: RNX,RNY,RNZ,WNX,WNY,WNZ,RMAG,WMAG
REAL		::	DOTPRODUCT,THETA,D12,CONST
REAL		::	TIMEMIN,SMIN
REAL		::	TEMP,ITEMP1,ITEMP2,EMASS,RHOP
REAL		::	PX,PY,PZ
REAL		::	RADA,PVOL,ESDV,MASS
! 
!***********************************************************************************
!LPX,LPY,LPZ = NUMBER OF DIVSIONS IN THE X,Y,Z DIRECTIONS TO FORM PARTICLE LISTS
!NPARSLICE(IJKP) = NUMBER OF PARTICLES IN SLICE I,J,K
!IPSLICE(IJKP,NN) = GLOBAL PARTICLE NUMBER OF THE NTH PARTICLE IN THE I,J,KTH SLICE
!
!***********************************************************************************
!
DO i=1,32
	COLLISION_SLAB(i) = 0
END DO

MIJK = LPX*LPY*LPZ
!
! *****PROBLEM SPECIFIC INFORMATION ON XMAXCOL,YMAXCOL,ZMAXCOL*********
XMAXCOL = DIAX
YMAXCOL = DIAY
ZMAXCOL = DIAZ

DXCOL = XMAXCOL/DFLOAT(LPX)
DYCOL = YMAXCOL/DFLOAT(LPY)
DZCOL = ZMAXCOL/DFLOAT(LPZ)

IF(.NOT. ALLOCATED(CFLAG))ALLOCATE(CFLAG(NPART))
CFLAG = 0
!
!********************************************************************************
!DIVIDE THE DOMAIN INTO SLICES AND DUMP THE APPROPRIATE PARTICLES INTO EACH SLICE
!********************************************************************************
!
!ZERO THE NUMBER OF PARTICLES IN EACH SLICE
JALLOC = 0
!
IF(.NOT. ALLOCATED(NPARSLICE))ALLOCATE(NPARSLICE(MIJK))
NPARSLICE = 0
!
!*******************************************
!COUNT THE NUMBER OF PARTICLES IN EACH SLICE
!*******************************************
!
DO JPART = 1, NPAR
   IF(INDEX_OUT(JPART) == 1)THEN
      CYCLE
   END IF
  	IPX = INT(PART(1,JPART)/DXCOL) + 1
  	IPY = INT(PART(2,JPART)/DYCOL) + 1
  	IPZ = INT(PART(3,JPART)/DZCOL) + 1
	! MODIFY IPX,IPY,IPZ INCASE PARTICLE IS LYING EXACTLY ON THE LAST BOUNDARY
  	IF (IPX > LPX) IPX = LPX
  	IF (IPY > LPY) IPY = LPY
  	IF (IPZ > LPZ) IPZ = LPZ
   IJKP = (IPX + (IPY-1)*LPX + (IPZ-1)*LPX*LPY)
   NPARSLICE(IJKP) = NPARSLICE(IJKP) + 1
  	IF(NPARSLICE(IJKP) > JALLOC) JALLOC = NPARSLICE(IJKP)
END DO
!
!**************************************
!DYNAMICALLY ALLOCATE ARRAY FOR IPSLICE
!**************************************
!
IF(.NOT. ALLOCATED(IPSLICE))ALLOCATE(IPSLICE(MIJK,JALLOC))
!
!*************************************************
!ZERO THE NUMBER OF PARTICLES IN EACH SLICE AGAIN
!SO IT MAY BE USED AS AN INDEX IN THE NEXT LOOP
!*************************************************
!
NPARSLICE = 0
!
!DUMP GLOBAL PARTICLE NUMBER INTO LOCAL BIN NUMBER AND RECOUNT NUMBER OF PARTICLES IN EACH BIN
DO JPART = 1, NPAR
   IF(INDEX_OUT(JPART) == 1)THEN
      CYCLE
   END IF
  	IPX = INT(PART(1,JPART)/DXCOL) + 1
  	IPY = INT(PART(2,JPART)/DYCOL) + 1
  	IPZ = INT(PART(3,JPART)/DZCOL) + 1
	! MODIFY IPX,IPY,IPZ INCASE PARTICLE IS LYING EXACTLY ON THE LAST BOUNDARY
  	IF (IPX > LPX) IPX = LPX
  	IF (IPY > LPY) IPY = LPY
  	IF (IPZ > LPZ) IPZ = LPZ
   IJKP = (IPX + (IPY-1)*LPX + (IPZ-1)*LPX*LPY)
   NPARSLICE(IJKP) = NPARSLICE(IJKP) + 1
   K = NPARSLICE(IJKP)
  	IPSLICE(IJKP,K) = JPART
END DO
!
!*************************************************
!LOOK FOR COLLISIONS BETWEEN PARTILCES IP1 AND IP2
!*************************************************
!
JALLOC = JALLOC*(JALLOC - 1)/2

IF(.NOT. ALLOCATED(IPAR))ALLOCATE(IPAR(JALLOC,2))
IF(.NOT. ALLOCATED(TIMECONTACT))ALLOCATE(TIMECONTACT(JALLOC))
!
DO IPX = 1, LPX
  	DO IPY = 1, LPY
    	DO IPZ = 1, LPZ
      	I= 0
         !
         IJKP = (IPX + (IPY-1)*LPX + (IPZ-1)*LPX*LPY)
         DO IP1 = 1, NPARSLICE(IJKP) - 1
        		DO IP2 = IP1 + 1, NPARSLICE(IJKP)
          		N1 = IPSLICE(IJKP,IP1)
          		N2 = IPSLICE(IJKP,IP2)
          		WNX = PART(4,N1) - PART(4,N2)
          		WNY = PART(5,N1) - PART(5,N2)
          		WNZ = PART(6,N1) - PART(6,N2)
          		RNX = PART(1,N1) - PART(1,N2)
          		RNY = PART(2,N1) - PART(2,N2)
          		RNZ = PART(3,N1) - PART(3,N2)
          		RMAG = DSQRT(RNX**2 + RNY**2 + RNZ**2)
          		WMAG = DSQRT(WNX**2 + WNY**2 + WNZ**2)
               DOTPRODUCT = WNX*RNX + WNY*RNY + WNZ*RNZ

          		IF(DOTPRODUCT > 0.0)CYCLE

               THETA = PIE - DACOS((DOTPRODUCT/RMAG/WMAG))
               TIMEMIN = LPTTIME + RMAG*DCOS(THETA)/WMAG
               SMIN = RMAG*DSIN(THETA)
               !DIAMETER OF PARTICLES #1 ADND #2
               DPAR1 = PART(8,N1)
               DPAR2 = PART(8,N2)
               D12 = (DPAR1 + DPAR2)/2.0

               IF(SMIN > D12)CYCLE
               IF(TIMEMIN < LPTTIME) CYCLE
               IF(TIMEMIN > LPTTIME + PTIM)CYCLE
          		!
          		!**************************************************
					!AT THIS POINT, IT HAS BEEN DETERMINED THAT THE TWO
               !PARTICLES COLLIDE. NOW FIND THE TIME OF IMPACT.
               !**************************************************
               !
          		I = I + 1
          		CONST = RMAG**2*WMAG**2*(1.0 - D12**2/RMAG**2)/DOTPRODUCT**2
          		TIMECONTACT(I) = TIMEMIN + DOTPRODUCT*DSQRT(1.0 - CONST)/WMAG**2
          		IPAR(I,1) = N1
          		IPAR(I,2) = N2
               NPARTPARTCOL = NPARTPARTCOL + 1
        		END DO
      	END DO
         !
         NN = I
      	!
         !************************************************
			!NOW SORT COLLISIONS ACCORDING TO TIME OF CONTACT
         !************************************************
         !
        	DO I = 1, NN - 1
         	DO J = I + 1, NN
             	IF(TIMECONTACT(J) < TIMECONTACT(I))THEN
               	TEMP = TIMECONTACT(I)
               	ITEMP1 = IPAR(I,1)
               	ITEMP2 = IPAR(I,2)
               	TIMECONTACT(I) = TIMECONTACT(J)
               	IPAR(I,1) = IPAR(J,1)
               	IPAR(I,2) = IPAR(J,2)
               	TIMECONTACT(J) = TEMP
               	IPAR(J,1) = ITEMP1
               	IPAR(J,2) = ITEMP2
             	END IF
           	END DO
         END DO
         !
         !********************************************************
			!CARRY OUT COLLISIONS FOR THE PARTICLE LIST IN THIS SLICE
         !********************************************************
      	DO I = 1, NN
            IF(CFLAG(IPAR(I,1)) == 0 .AND. CFLAG(IPAR(I,2)) == 0)THEN
               !ADVANCE THIS PARTICLE PAIR TO THE POINT OF IMPACT
               !AND FIND THE NEW RELATIVE POSITION VECTOR
               !
               PART0(1,IPAR(I,1)) = PART(1,IPAR(I,1))
               PART0(2,IPAR(I,1)) = PART(2,IPAR(I,1))
               PART0(3,IPAR(I,1)) = PART(3,IPAR(I,1))

               !THIS IS THE COORDINATES OF THE CENTRES OF PARTICLE 1 A THE INSTANT OF COLLISION
               PART(1,IPAR(I,1)) = PART(1,IPAR(I,1)) + (TIMECONTACT(I) - LPTTIME)*PART(4,IPAR(I,1))
               PART(2,IPAR(I,1)) = PART(2,IPAR(I,1)) + (TIMECONTACT(I) - LPTTIME)*PART(5,IPAR(I,1))
               PART(3,IPAR(I,1)) = PART(3,IPAR(I,1)) + (TIMECONTACT(I) - LPTTIME)*PART(6,IPAR(I,1))
               !
               !*****************************************************************************
               !CHECK IF PARTICLE 1 HAS CROSSED A BOUNDARY APPLY PARTICLE BOUNDARY CONDITIONS
               !CHECK FOR DEPOSITION HERE
               !*****************************************************************************
               !
               ISET = 1
               JPART = IPAR(I,1)
               !----------------------
               CALL PART_PWI_ELLIPSOID
               !----------------------
               !
               PART0(1,IPAR(I,2)) = PART(1,IPAR(I,2))
               PART0(2,IPAR(I,2)) = PART(2,IPAR(I,2))
               PART0(3,IPAR(I,2)) = PART(3,IPAR(I,2))

               !THIS IS THE COORDINATES OF THE CENTRES OF PARTICLE 2 A THE INSTANT OF COLLISION
               PART(1,IPAR(I,2)) = PART(1,IPAR(I,2)) + (TIMECONTACT(I) - LPTTIME)*PART(4,IPAR(I,2))
               PART(2,IPAR(I,2)) = PART(2,IPAR(I,2)) + (TIMECONTACT(I) - LPTTIME)*PART(5,IPAR(I,2))
               PART(3,IPAR(I,2)) = PART(3,IPAR(I,2)) + (TIMECONTACT(I) - LPTTIME)*PART(6,IPAR(I,2))               !
               !
               !*****************************************************************************
               !CHECK IF PARTICLE 2 HAS CROSSED A BOUNDARY APPLY PARTICLE BOUNDARY CONDITIONS
               !CHECK FOR DEPOSITION HERE
               !*****************************************************************************
               !
               ISET = 1
               JPART = IPAR(I,2)
               !----------------------
               CALL PART_PWI_ELLIPSOID
               !----------------------
               !
               !
               !IF(TEST IF PARTICLE HAS DEPOSITED)CYCLE
               !
               !********************************************************
               !POST COLLISION PROPERTIES OF THE TWO COLLIDING PARTICLES
               !m1V1i + m2V2i = miV1f + m2V2f
               !********************************************************
               !
               CFLAG(IPAR(I,1)) = 1
               CFLAG(IPAR(I,2)) = 1

               !COUNT AND RECORD THE NUMBER OF COLLISION
               COLLISION_COUNT(0) = COLLISION_COUNT(0) + 1
			   
			   if(mod(istep,500) .EQ. 0)
			   
				do ii=1,32
						if(PART(1,IPAR(I,2)) .GE. (ii-1)*(0.04/32) .AND. PART(1,IPAR(I,2)) < (ii)*(0.04/32)
							COLLISION_SLAB(ii)++ 
						end if
				end do
				end if	
			   
			   WRITE(*,*) '================== COLLISION COUNTING ==================='
			   do ii=1,32
					WRITE(*,*) ((ii-1)*(0.04/32) + (0.02/32)),COLLISION_SLAB(ii) 
			   end do
			   WRITE(*,*) '================== END OF COLLISION COUNTING ==================='
               !
              IF(POSTCOLLISION == 1)THEN
                  !
                  !***************
                  !PARTICLE BOUNCE
                  !***************
                  !
                  WNX = PART(4,IPAR(I,1)) - PART(4,IPAR(I,2))
                  WNY = PART(5,IPAR(I,1)) - PART(5,IPAR(I,2))
                  WNZ = PART(6,IPAR(I,1)) - PART(6,IPAR(I,2))
                  RNX = PART(1,IPAR(I,1)) - PART(1,IPAR(I,2))
                  RNY = PART(2,IPAR(I,1)) - PART(2,IPAR(I,2))
                  RNZ = PART(3,IPAR(I,1)) - PART(3,IPAR(I,2))
                  RMAG = DSQRT(RNX**2 + RNY**2 + RNZ**2)
                  RNX = RNX/RMAG
                  RNY = RNY/RMAG
                  RNZ = RNZ/RMAG

                  DOTPRODUCT =WNX*RNX + WNY*RNY + WNZ*RNZ
                  DPAR1 = PART(8,IPAR(I,1))
                  RHOP = PSTATE(06,ISET)
                  EMASS = RHOP*PIE*DPAR1**3/6.0

                  PX = -EMASS*DOTPRODUCT*RNX
                  PY = -EMASS*DOTPRODUCT*RNY
                  PZ = -EMASS*DOTPRODUCT*RNZ

                  PART(4,IPAR(I,1)) =  PX/EMASS + PART(4,IPAR(I,1))
                  PART(5,IPAR(I,1)) =  PY/EMASS + PART(5,IPAR(I,1))
                  PART(6,IPAR(I,1)) =  PZ/EMASS + PART(6,IPAR(I,1))
                  PART(4,IPAR(I,2)) = -PX/EMASS + PART(4,IPAR(I,2))
                  PART(5,IPAR(I,2)) = -PY/EMASS + PART(5,IPAR(I,2))
                  PART(6,IPAR(I,2)) = -PZ/EMASS + PART(6,IPAR(I,2))


                  PART0(1,IPAR(I,1)) = PART(1,IPAR(I,1))
                  PART0(2,IPAR(I,1)) = PART(2,IPAR(I,1))
                  PART0(3,IPAR(I,1)) = PART(3,IPAR(I,1))

                  PART(1,IPAR(I,1)) = PART(1,IPAR(I,1)) + (LPTTIME + PTIM - TIMECONTACT(I))*PART(4,IPAR(I,1))
                  PART(2,IPAR(I,1)) = PART(2,IPAR(I,1)) + (LPTTIME + PTIM - TIMECONTACT(I))*PART(5,IPAR(I,1))
                  PART(3,IPAR(I,1)) = PART(3,IPAR(I,1)) + (LPTTIME + PTIM - TIMECONTACT(I))*PART(6,IPAR(I,1))

                  !
                  !*****************************************************************************
                  !CHECK IF PARTICLE 1 HAS CROSSED A BOUNDARY APPLY PARTICLE BOUNDARY CONDITIONS
                  !CHECK FOR DEPOSITION HERE
                  !*****************************************************************************
                  !
                  ISET = 1
                  JPART = IPAR(I,1)
                  !----------------------
                  CALL PART_PWI_ELLIPSOID
                  !----------------------
                  !
                  PART0(1,IPAR(I,2)) = PART(1,IPAR(I,2))
                  PART0(2,IPAR(I,2)) = PART(2,IPAR(I,2))
                  PART0(3,IPAR(I,2)) = PART(3,IPAR(I,2))

                  PART(1,IPAR(I,2)) = PART(1,IPAR(I,2)) + (LPTTIME + PTIM - TIMECONTACT(I))*PART(4,IPAR(I,2))
                  PART(2,IPAR(I,2)) = PART(2,IPAR(I,2)) + (LPTTIME + PTIM - TIMECONTACT(I))*PART(5,IPAR(I,2))
                  PART(3,IPAR(I,2)) = PART(3,IPAR(I,2)) + (LPTTIME + PTIM - TIMECONTACT(I))*PART(6,IPAR(I,2))
                  !
                  !*****************************************************************************
                  !CHECK IF PARTICLE 2 HAS CROSSED A BOUNDARY APPLY PARTICLE BOUNDARY CONDITIONS
                  !CHECK FOR DEPOSITION HERE
                  !*****************************************************************************
                  !
                  ISET = 1
                  JPART = IPAR(I,2)
                  !----------------------
                  CALL PART_PWI_ELLIPSOID
                  !----------------------
                  !
                  !********************************************************
                  !COUNT THE NUMBER OF COLLISIONS EACH PARTICLE IS INVOLVED
                  !********************************************************
                  INDEX_COL_COUNT(IPAR(I,1)) = INDEX_COL_COUNT(IPAR(I,1)) + 1
                  INDEX_COL_COUNT(IPAR(I,2)) = INDEX_COL_COUNT(IPAR(I,2)) + 1
                  !
                  COLLISION_COUNT(POSTCOLLISION) = COLLISION_COUNT(POSTCOLLISION) + 1
                  !
                  !OUTPUT TIMESTEP, TIME[s], TIME[+], TOTAL COLLISION, TOTAL BOUNCE, TOTAL COALECENCE, POSITION OF DAUGHTER, DIAMETER, MASS
						WRITE(COLFILE(POSTCOLLISION),'(I10,2X,2(1PE15.6),2X,3(I10,2X),4(1PE15.6,2X))')LPTSTEP,LPTTIME,LPTTIME*WTIM,(COLLISION_COUNT(II), II = 0, 2),PART(01,IPAR(I,1))*WLEN,PART(8,IPAR(I,1))*1E6,PART(10,IPAR(I,1))
                  !
               ELSE IF(POSTCOLLISION == 2)THEN
                  !
                  !*************************************************************************************************
                  !PARTICLE COALECENCE:
                  !TWO PARTICLES WOULD COMBINE INTANTANEOUSLY UPON CONTACT TO FORM LARGER PARTICLE OF THE SAME SHAPE
                  !*************************************************************************************************
                  !
                  DPAR1 = PART(8,IPAR(I,1))
                  DPAR2 = PART(8,IPAR(I,2))
                  DPAR13 = DPAR1**3
                  DPAR23 = DPAR2**3
                  !COMPUTE THE CHILD PARTICLE DIAMETER USING CONSERVATION OF MASS, ASSUMING DENSITY IS CONSTANT
                  DPAR3 = (DPAR13 + DPAR23)**ONETHIRD
                  DPAR33 = DPAR3**3
                  !ASSIGN THE CHILD PARTICLE DIAMETER TO PARTICLE 1
                  PART(8,IPAR(I,1)) = DPAR3
                  !UPDATE THE GLOBAL PARTICLE PROPERTIES
                  ESDV = PART(8,IPAR(I,1))
                  !NEW PARTICLE VOLUME
      				PVOL = PIE*(1.0D0/6.0D0)*ESDV**3
                  !NEW PARTICLE RADIUS
      				RADA = (PVOL/(FOURTHIRD*LAMBDA_SET(ISET)*PIE))**ONETHIRD
						!NEW PARTICLE ACTUAL AREA OF NON-SPHERICAL PARTICLE
   					!PASA = 4.0D0*PIE*RADA**2
						!!SURFACE, PROJECTED AREAS
   					!SPHVOL=PIE*ESDV**2
						!!PARTICLE SPHERICITY, PSN
   					!PSN = SPHVOL/PASA
                  !MASS OF THE NEW PARTICLE
                  MASS = DENP_SET(ISET)*PVOL
                  !SAVE TO ARRAY
                  PART(8,JPART) = ESDV
      				PART(9,JPART) = RADA
      				PART(10,JPART) = MASS

                  !COMPUTE THE VELOCITY OF THE NEW PARTICLE
                  PART(4,IPAR(I,1)) =  (PART(4,IPAR(I,1))*DPAR13 + PART(4,IPAR(I,2))*DPAR23)/DPAR33
                  PART(5,IPAR(I,1)) =  (PART(5,IPAR(I,1))*DPAR13 + PART(5,IPAR(I,2))*DPAR23)/DPAR33
                  PART(6,IPAR(I,1)) =  (PART(6,IPAR(I,1))*DPAR13 + PART(6,IPAR(I,2))*DPAR23)/DPAR33

                  PART0(1,IPAR(I,1)) = PART(1,IPAR(I,1))
                  PART0(2,IPAR(I,1)) = PART(2,IPAR(I,1))
                  PART0(3,IPAR(I,1)) = PART(3,IPAR(I,1))

                  PART(1,IPAR(I,1)) = 0.5*(PART(1,IPAR(I,1)) + PART(1,IPAR(I,2))) + (LPTTIME + PTIM - TIMECONTACT(I))*PART(4,IPAR(I,1))
                  PART(2,IPAR(I,1)) = 0.5*(PART(2,IPAR(I,1)) + PART(2,IPAR(I,2))) + (LPTTIME + PTIM - TIMECONTACT(I))*PART(5,IPAR(I,1))
                  PART(3,IPAR(I,1)) = 0.5*(PART(3,IPAR(I,1)) + PART(3,IPAR(I,2))) + (LPTTIME + PTIM - TIMECONTACT(I))*PART(6,IPAR(I,1))
                  !
                  !*****************************************************************************
                  !CHECK IF PARTICLE CHILD HAS CROSSED A BOUNDARY APPLY PARTICLE BOUNDARY CONDITIONS
                  !CHECK FOR DEPOSITION HERE
                  !*****************************************************************************
                  !
                  ISET = 1
                  JPART = IPAR(I,1)
                  !PARTICLE 2 WILL BE OUT
                  INDEX_OUT(IPAR(I,2)) = 1
                  !COUNT THE NUMBER OF COLLISIONS FOR CHILD PARTICLE
                  INDEX_COL_COUNT(IPAR(I,1)) = INDEX_COL_COUNT(IPAR(I,1)) + INDEX_COL_COUNT(IPAR(I,2))
                  !
                  COLLISION_COUNT(POSTCOLLISION) = COLLISION_COUNT(POSTCOLLISION) + 1
                  !
                  !OUTPUT TIMESTEP, TIME[s], TIME[+], TOTAL COLLISION, TOTAL BOUNCE, TOTAL COALECENCE, POSITION OF DAUGHTER, DIAMETER, MASS
						WRITE(COLFILE(POSTCOLLISION),'(I10,2X,2(1PE15.6),2X,3(I10,2X),4(1PE15.6,2X))')LPTSTEP,LPTTIME,LPTTIME*WTIM,(COLLISION_COUNT(II), II = 0, 2),PART(01,IPAR(I,1))*WLEN,PART(8,IPAR(I,1))*1E6,PART(10,IPAR(I,1))
                  !
                  !----------------------
                  CALL PART_PWI_ELLIPSOID
                  !----------------------
                  !
               end if
               !
					!OUTPUT TIMESTEP, TIME[s], TIME[+], TOTAL COLLISION, TOTAL BOUNCE, TOTAL COALECENCE, POSITION OF DAUGHTER, DIAMETER, MASS
					WRITE(COLFILE(0),'(I10,2X,2(1PE15.6),2X,3(I10,2X),4(1PE15.6,2X))')LPTSTEP,LPTTIME,LPTTIME*WTIM,(COLLISION_COUNT(II), II = 0, 2),PART(01,IPAR(I,1))*WLEN,PART(8,IPAR(I,1))*1E6,PART(10,IPAR(I,1))
               !
         	END IF
      	END DO
   	END DO
	END DO
END DO

IF(ALLOCATED(IPSLICE))		DEALLOCATE(IPSLICE)
IF(ALLOCATED(IPAR))			DEALLOCATE(IPAR)
IF(ALLOCATED(TIMECONTACT))	DEALLOCATE(TIMECONTACT)

RETURN
END SUBROUTINE PART_COLLISION










! -*- f90 -*-
MODULE DES_TIME_MARCH

      
  use discretelement
  use functions
  use machine
  
  use output, only: DLB,DLB_TIME
  
  use output_man, only: OUTPUT_MANAGER
  
  
  use run, only: NSTEP
  use run, only: TIME, TSTOP, DT

  
!---------------------------------------------------------------------//
      ! Total number of particles
      INTEGER, SAVE :: NP=0, P

      ! loop counter index for any initial particle settling incoupled cases
      INTEGER :: FACTOR
      ! Temporary variables when des_continuum_coupled is T to track
      ! changes in solid time step
      DOUBLE PRECISION :: DTSOLID_TMP

      LOGICAL :: EXIT_LOOP

      DOUBLE PRECISION ::  RADIUS, EN

      INTEGER :: NLGS, C, CC, CC_START, CC_END, I, IERR, NB_CONTACTS

      DOUBLE PRECISION :: DX, DY, DZ
      DOUBLE PRECISION :: NX, NY, NZ

      DOUBLE PRECISION :: S
      
      DOUBLE PRECISION :: UN_NEW,UT_NEW
      DOUBLE PRECISION :: UN_OLD,UT_OLD
      DOUBLE PRECISION :: VM_N,VM_T
      DOUBLE PRECISION :: TAU_N, GAMMA, MEQ
      
!......................................................................!

    CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!     Subroutine: DES_TIME_INIT                                        !
!     Author: SOUFIANE                                Date: 21-Jun-04  !
!                                                                      !
!     Purpose: Main DEM driver routine                                 !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DES_TIME_INIT
        USE DISCRETELEMENT
        IMPLICIT NONE
       
        PRINT*,"DES_TIME_INIT_ME"
        
        RADIUS = HUGE(0D0)
        DO P = 1, PARTICLES
           RADIUS = MIN(RADIUS,DES_RADIUS(P))
        END DO
        
        EXIT_LOOP = .FALSE.
        
        !> Initialize time stepping variable for pure granular simulations.
        FACTOR = CEILING(real((TSTOP-TIME)/DTSOLID))
        
        DT = 0.1 * RADIUS
        EN = 0.9
        
        !> PARTICLES INITIALIZATION
        DO P = 1, PARTICLES
           
           !> MFIX VARIABLES --> NEW ME VARIABLES
           PPOS_ME(P,1:3) = PPOS(P,1:3)
           
           DES_POS_NEW_ME(P,1:3) = DES_POS_NEW(P,1:3)
           
           DES_VEL_NEW_ME(P,1:3) = DES_VEL_NEW(P,1:3)
           
           RO_ME(P) = RO_Sol(P)
           DES_RADIUS_ME(P) = DES_RADIUS(P)
           PMASS_ME(P) = RO_ME(P) * ACOS(-1.0) *  DES_RADIUS_ME(P)**2
           
           FC_OLD_ME(P,1) = GRAV(1) * PMASS_ME(P)
           FC_OLD_ME(P,2) = GRAV(2) * PMASS_ME(P)
           FC_OLD_ME(P,3) = GRAV(3) * PMASS_ME(P)
           
           FC_NEW_ME(P,1) = GRAV(1) * PMASS_ME(P)
           FC_NEW_ME(P,2) = GRAV(2) * PMASS_ME(P)
           FC_NEW_ME(P,3) = GRAV(3) * PMASS_ME(P)
           
           !> MY CFUPDATEOLD
           DES_POS_OLD_ME(P,1:3) = DES_POS_NEW_ME(P,1:3)
           
           DES_VEL_OLD_ME(P,1:3) = DES_VEL_NEW_ME(P,1:3)
           
           !> OLD ME VARIABLES --> FREE ME VARIABLES
           DES_POS_FREE_ME(P,1:3) = DES_POS_OLD_ME(P,1:3)
           
           DES_VEL_FREE_ME(P,1:3) = DES_VEL_OLD_ME(P,1:3)
           
           !> OLD ME VARIABLES --> DEMI ME VARIABLES
           DES_POS_DEMI_ME(P,1:3) = DES_POS_OLD_ME(P,1:3)
           
           DES_VEL_DEMI_ME(P,1:3) = DES_VEL_OLD_ME(P,1:3)
            
        END DO
        
        CALL EXPORT_ME('OUT.DAT',DES_POS_OLD_ME,DES_VEL_OLD_ME,DES_RADIUS_ME,&
             FC_OLD_ME,NP)
        
        CALL OUTPUT_MANAGER(.FALSE., .FALSE.)
        
       
      END SUBROUTINE DES_TIME_INIT


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!     Subroutine: DES_TIME_STEP                                        !
!     Author: Jay Boyalakuntla                        Date: 21-Jun-04  !
!                                                                      !
!     Purpose: Main DEM driver routine                                 !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DES_TIME_STEP(NN)
         use compar, only: ADJUST_PARTITION

! Modules
!---------------------------------------------------------------------//
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NN
         
         PRINT*,"DES_TIME_STEP_ME: ",NN
                  
         DO P =1, PARTICLES
            PRINT*,"======================================================"
            PRINT*,"POSITION OLD PARTICLE:",P,": ",DES_POS_OLD_ME(P,1:3) 
            PRINT*,"VELOCITY OLD PARTICLE:",P,": ",DES_VEL_OLD_ME(P,1:3)
            PRINT*,"======================================================"
         END DO
         
         !CALL MPI_FINALIZE(IERR)
         !STOP
                  
         DO P = 1, PARTICLES
            
            !> Q_{K+1/2} =  Q_{K} + 0.5 DT \DOT{Q}_{K} 
            DES_POS_DEMI_ME(P,1) = DES_POS_OLD_ME(P,1) + 0.5*DT*DES_VEL_OLD_ME(P,1)
            DES_POS_DEMI_ME(P,2) = DES_POS_OLD_ME(P,2) + 0.5*DT*DES_VEL_OLD_ME(P,2)
            DES_POS_DEMI_ME(P,3) = DES_POS_OLD_ME(P,3) + 0.5*DT*DES_VEL_OLD_ME(P,3)
            
            !> THETA-SCHEMA 0.5
            
            !PRINT*,DES_VEL_FREE_ME(P,1)
            !PRINT*,DES_VEL_FREE_ME(P,2)
            
            DES_VEL_FREE_ME(P,1) = DES_VEL_OLD_ME(P,1) + 0.5*DT*(FC_OLD_ME(P,1)+&
                 FC_NEW_ME(P,1))/PMASS_ME(P)
            DES_VEL_FREE_ME(P,2) = DES_VEL_OLD_ME(P,2) + 0.5*DT*(FC_OLD_ME(P,2)+&
                 FC_NEW_ME(P,2))/PMASS_ME(P)
            DES_VEL_FREE_ME(P,3) = DES_VEL_OLD_ME(P,3) + 0.5*DT*(FC_OLD_ME(P,3)+&
                 FC_NEW_ME(P,3))/PMASS_ME(P)
                        
            DES_POS_FREE_ME(P,1) = DES_POS_OLD_ME(P,1) + 0.5*DT*(DES_VEL_FREE_ME(P,1)&
                 +DES_VEL_OLD_ME(P,1))
            DES_POS_FREE_ME(P,2) = DES_POS_OLD_ME(P,2) + 0.5*DT*(DES_VEL_FREE_ME(P,2)&
                 +DES_VEL_OLD_ME(P,2))
            DES_POS_FREE_ME(P,3) = DES_POS_OLD_ME(P,3) + 0.5*DT*(DES_VEL_FREE_ME(P,3)&
                 +DES_VEL_OLD_ME(P,3))

            !PRINT*,DES_VEL_FREE_ME(P,1)
            !PRINT*,DES_VEL_FREE_ME(P,2)
            
            !> DOT{Q}_{0}{K+1} =  DOT{Q}_{FREE}{K}
            DES_VEL_NEW_ME(P,1) = DES_VEL_FREE_ME(P,1)
            DES_VEL_NEW_ME(P,2) = DES_VEL_FREE_ME(P,2)
            DES_VEL_NEW_ME(P,3) = DES_VEL_FREE_ME(P,3)

            !PRINT*,DES_VEL_NEW_ME(P,1)
            !PRINT*,DES_VEL_NEW_ME(P,2)
            !PRINT*,DES_VEL_NEW_ME(P,3)
            
         END DO

         !CALL MPI_FINALIZE(IERR)
         !STOP
         
         !> PARTICLES IN MOTION WITHOUT CONTACT FORCES
         DO P = 1, PARTICLES
            DES_POS_NEW(P,1:3) = DES_POS_FREE_ME(P,1:3)
            DES_VEL_NEW(P,1:3) = DES_VEL_NEW_ME(P,1:3)
            !PRINT*,"======================================================"
            !PRINT*,"POSITION FREE PARTICLE:",P,": ",DES_POS_NEW(P,1:3)
            !PRINT*,"VELOCITY FREE PARTICLE:",P,": ",DES_VEL_NEW(P,1:3)
            !PRINT*,"======================================================"            
         END DO

         !CALL MPI_FINALIZE(IERR)
         !STOP

         NB_CONTACTS = 0

         !> NEIGHBOR SEARCH
         CALL NEIGHBOUR_ME

         !CALL MPI_FINALIZE(IERR)
         !STOP
         
         !> INITIALIZATION OF ACTIVE SET PARAMETERS
         PN = 0
         VM = 0
         DPN = 0
         ACTIVE = .FALSE.
         
         !> NLGS LOOP
         DO NLGS = 1, 3
         
            DO P = 1, PARTICLES
               
               NB_CONTACTS = 0
               CC_START = 1
               IF (P.GT.1) CC_START = NEIGHBOR_INDEX_ME(P-1)
               CC_END = NEIGHBOR_INDEX_ME(P)
               NB_CONTACTS = NB_CONTACTS + CC_END - CC_START
               
               !> CONTACTS LOOP
               DO CC = CC_START, CC_END-1
                  I = NEIGHBORS_ME(CC)
                  IF(IS_NONEXISTENT(I)) CYCLE
                  ! IF (I<=P) CYCLE !> À TESTER
                  PRINT*,"! ==================================================== !"
                  PRINT*,"! NOMBRE DE CONTACTS: ",NB_CONTACTS,"                    !"
                  PRINT*,"! CONNECTIVITE P -- I: ",P,"<==>",I,"  !"
                  PRINT*,"! ==================================================== !"
                  
! ==================================================================== !               
! ================== CONTACT CONDITION + ACTIVE SET ================== !      
! ==================================================================== !                                       
                  DX = DES_POS_DEMI_ME(P,1) - DES_POS_DEMI_ME(I,1)
                  DY = DES_POS_DEMI_ME(P,2) - DES_POS_DEMI_ME(I,2)
                  DZ = DES_POS_DEMI_ME(P,3) - DES_POS_DEMI_ME(I,3)
                  
                  !CALL MPI_FINALIZE(IERR)
                  !STOP
                  
                  PRINT*,"--> DX: ",DX
                  
                  !NEW VARIABLE S = SQRT(DX**2 + DY**2 + DZ**2)
                  NX = DX / SQRT(DX**2 + DY**2 + DZ**2)
                  NY = DY / SQRT(DX**2 + DY**2 + DZ**2)
                  NZ = DZ / SQRT(DX**2 + DY**2 + DZ**2)
                  
                  PRINT*,"--> NX: ",NX
                  
                  S = SQRT(DX**2+DY**2+DZ**2)-(DES_RADIUS_ME(P)+DES_RADIUS_ME(I))
                  
                  PRINT*,"--> S: ",S

                  !CALL MPI_FINALIZE(IERR)
                  !STOP
                  
                  IF (S<0) THEN     ! CONTACT CONDITION
                     
                     UN_OLD = (DES_VEL_OLD_ME(P,1)-DES_VEL_OLD_ME(I,1))*NX + &
                          (DES_VEL_OLD_ME(P,2)-DES_VEL_OLD_ME(I,2))*NY + &
                          (DES_VEL_OLD_ME(P,3)-DES_VEL_OLD_ME(I,3))*NZ
                        
                     PRINT*,"--> UN_OLD: ",UN_OLD
                     
                     UN_NEW = (DES_VEL_NEW_ME(P,1)-DES_VEL_NEW_ME(I,1))*NX + &
                          (DES_VEL_NEW_ME(P,2)-DES_VEL_NEW_ME(I,2))*NY + &
                          (DES_VEL_NEW_ME(P,3)-DES_VEL_NEW_ME(I,3))*NZ
                     
                     PRINT*,"--> UN_NEW: ",UN_NEW
                     !PRINT*,PMASS_ME(P),PMASS_ME(I)
                     !CALL MPI_FINALIZE(IERR)
                     !STOP
                     
                     VM(CC) = (UN_NEW + UN_OLD * EN) / (1 + EN)
                     
                     PRINT*,"--> VM(CC): ",VM(CC)
                     
                     MEQ = 2*PMASS_ME(P)*PMASS_ME(I) / (PMASS_ME(P) + PMASS_ME(I))
                     
                     PRINT*,"--> MEQ: ",MEQ
                     
                     GAMMA = 1E3                                          
                     TAU_N = PN(CC) - GAMMA * VM(CC) * MEQ / DT
                     
                     PRINT*,"--> TAU_N: ",TAU_N
                     
                     !CALL MPI_FINALIZE(IERR)
                     !STOP
                     
                     !> ACTIVE SET
                     IF (TAU_N>0) THEN
                        DPN(CC) = - VM(CC) * MEQ / DT * EN
                        ACTIVE(CC) = .TRUE.
                     ELSE
                        DPN(CC) = 0
                     END IF
                     
                     PRINT*,"--> DPN(CC): ",DPN(CC)
                     
                     !CALL MPI_FINALIZE(IERR)
                     !STOP
                     
                     PRINT*,"1 VEL NEW P: ",DES_VEL_NEW_ME(P,1)
                     PRINT*,"1 VEL NEW I: ",DES_VEL_NEW_ME(I,1)
                     PRINT*,((-DPN(CC))*DT/PMASS_ME(P)*NX)
                     
                     DES_VEL_NEW_ME(P,1)=DES_VEL_NEW_ME(P,1)+(DPN(CC))*DT/PMASS_ME(P)*NX
                     DES_VEL_NEW_ME(P,2)=DES_VEL_NEW_ME(P,2)+(DPN(CC))*DT/PMASS_ME(P)*NY
                     DES_VEL_NEW_ME(P,3)=DES_VEL_NEW_ME(P,3)+(DPN(CC))*DT/PMASS_ME(P)*NZ
                     
                     DES_VEL_NEW_ME(I,1)=DES_VEL_NEW_ME(I,1)+(-DPN(CC))*DT/PMASS_ME(I)*NX
                     DES_VEL_NEW_ME(I,2)=DES_VEL_NEW_ME(I,2)+(-DPN(CC))*DT/PMASS_ME(I)*NY
                     DES_VEL_NEW_ME(I,3)=DES_VEL_NEW_ME(I,3)+(-DPN(CC))*DT/PMASS_ME(I)*NZ
                     
                     PRINT*,"2 VEL NEW P: ",DES_VEL_NEW_ME(P,1)
                     PRINT*,"2 VEL NEW I: ",DES_VEL_NEW_ME(I,1)
                     
                     !CALL MPI_FINALIZE(IERR)
                     !STOP
                     
                     PN(CC) = PN(CC) + DPN(CC)
                     PRINT*,"--> PN(CC): ",PN(CC)
                        
                     !CALL MPI_FINALIZE(IERR)
                     !STOP
                     
                     !PRINT*,"1 POS FREE P: ",DES_POS_FREE_ME(P,1:3)
                     !PRINT*,"1 POS FREE I: ",DES_POS_FREE_ME(I,1:3)
                     
                     DES_POS_FREE_ME(P,1)=DES_POS_OLD_ME(P,1) + DT*DES_VEL_NEW_ME(P,1)
                     DES_POS_FREE_ME(P,2)=DES_POS_OLD_ME(P,2) + DT*DES_VEL_NEW_ME(P,2)
                     DES_POS_FREE_ME(P,3)=DES_POS_OLD_ME(P,3) + DT*DES_VEL_NEW_ME(P,3)
                     
                     DES_POS_FREE_ME(I,1)=DES_POS_OLD_ME(I,1) + DT*DES_VEL_NEW_ME(I,1)
                     DES_POS_FREE_ME(I,2)=DES_POS_OLD_ME(I,2) + DT*DES_VEL_NEW_ME(I,2)
                     DES_POS_FREE_ME(I,3)=DES_POS_OLD_ME(I,3) + DT*DES_VEL_NEW_ME(I,3)
                     
                     !PRINT*,"2 POS FREE P: ",DES_POS_FREE_ME(P,1:3)
                     !PRINT*,"2 POS FREE I: ",DES_POS_FREE_ME(I,1:3)
                     
                     !CALL MPI_FINALIZE(IERR)
                     !STOP
                     
                  ELSE
                     PN(CC) = 0
                  END IF
                     
               END DO

               IF ((CC_END - CC_START) == 0) THEN
                  PRINT*,"! ==================================================== !"
                  PRINT*,"! CONNECTIVITE P -- P: ",P,"<==>",P,"  !"
                  PRINT*,"! ==================================================== !"
                  EXIT
               END IF

            END DO
            
            !CALL MPI_FINALIZE(IERR)
            !STOP
                       
            !> CONDITION TO EXIT NLGS LOOP
            !> PROBLEME VENAIT DE LA !!!!!
            IF (MAXVAL(ABS(DT*DPN(1:NEIGHBOR_INDEX_ME(PARTICLES)-1))) .LT.1E-16) THEN
               PRINT*,"EXIT NLGS LOOP",MAXVAL(ABS(DT*DPN(1:NEIGHBOR_INDEX_ME(PARTICLES)-1))) 
               EXIT 
            ELSE
               PRINT*,"NLGS LOOP AGAIN"
            END IF
            
         END DO
                  
         !CALL MPI_FINALIZE(IERR)
         !STOP
               
         DO P = 1, PARTICLES
            
            DES_POS_NEW_ME(P,1) = DES_POS_DEMI_ME(P,1) + 0.5*DT*DES_VEL_NEW_ME(P,1)
            DES_POS_NEW_ME(P,2) = DES_POS_DEMI_ME(P,2) + 0.5*DT*DES_VEL_NEW_ME(P,2)
            DES_POS_NEW_ME(P,3) = DES_POS_DEMI_ME(P,3) + 0.5*DT*DES_VEL_NEW_ME(P,3)

            !> NEW ME VARIABLES --> OLD ME VARIABLES 
            DES_POS_OLD_ME(P,1:3) = DES_POS_NEW_ME(P,1:3)  
            DES_VEL_OLD_ME(P,1:3) = DES_VEL_NEW_ME(P,1:3)
            
            !PRINT*,DES_POS_OLD_ME(P,1)
            !PRINT*,DES_VEL_OLD_ME(P,1)
            
            !> NEW ME VARIABLES --> MFIX VARIABLES 
            DES_POS_NEW(P,1:3) = DES_POS_NEW_ME(P,1:3)
            DES_VEL_NEW(P,1:3) = DES_VEL_NEW_ME(P,1:3)
            
            !PRINT*,DES_POS_NEW(P,1)
            !PRINT*,DES_VEL_NEW(P,1)
            
         END DO
         
         !CALL MPI_FINALIZE(IERR)
         !STOP
             
! ==================================================================== !                 
! ==================================================================== !               
! ==================================================================== !                        
         ! Update time to reflect changes
         S_TIME = S_TIME + DTSOLID
         
         ! Keep track of TIME and number of steps for DEM simulations
         TIME = S_TIME
         NSTEP = NSTEP + 1
      
         DLB = .TRUE.
         
         ! Call the output manager to write RES and SPx data.
         
         CALL OUTPUT_MANAGER(.FALSE., .FALSE.)
         
         IF(ADJUST_PARTITION) THEN
            EXIT_LOOP = .TRUE.
            RETURN
         ENDIF
         
         IF(NN == 5) THEN
            CALL MPI_FINALIZE(IERR)
            STOP
         END IF
         
       END SUBROUTINE DES_TIME_STEP
    
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!     Subroutine: DES_TIME_END                                         !
!     Author: Jay Boyalakuntla                        Date: 21-Jun-04  !
!                                                                      !
!     Purpose: Main DEM driver routine                                 !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
    SUBROUTINE DES_TIME_END
        IMPLICIT NONE
         
         print*,"DES_TIME_END_ME"
         
         ! Reset the discrete time step to original value.
         DTSOLID = DTSOLID_TMP
         
         RETURN
       END SUBROUTINE DES_TIME_END
     END MODULE DES_TIME_MARCH
     





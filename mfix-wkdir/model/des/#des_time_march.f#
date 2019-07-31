! -*- f90 -*-
MODULE DES_TIME_MARCH

      
  use discretelement
  use des_allocate
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

      DOUBLE PRECISION :: EX(2)
      
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
        use geometry, only: X_MIN, X_MAX, Y_MIN, Y_MAX, Z_MIN, Z_MAX
        IMPLICIT NONE
        
        PRINT*,"! ========================================= ! "
        PRINT*,"! =========> INIT DEM SIMULATIONS <======== ! "
        PRINT*,"! ========================================= ! "
        
        
        RADIUS = HUGE(0D0)
        DO P = 1, PARTICLES
           RADIUS = MIN(RADIUS,DES_RADIUS(P))
        END DO
        
        EXIT_LOOP = .FALSE.
        
        !> Initialize time stepping variable for pure granular simulations.
        FACTOR = CEILING(real((TSTOP-TIME)/DTSOLID))
        
        DT = 0.1 * RADIUS
        !DT = 1E-2
        EN = 0.5
        
        DT = SQRT(2*(0.5-RADIUS)/10.)/100.
        DT = DT/2

        !PRINT*,"DT",DT,EN,RADIUS
        !STOP
        
        !> PARTICLES INITIALIZATION
        DO P = 1, PARTICLES
           
           !> MFIX VARIABLES --> NEW ME VARIABLES
           PPOS_ME(P,1:3) = PPOS(P,1:3)
           DES_POS_OLD(P,1:3) = DES_POS_NEW(P,1:3)
           DES_VEL_OLD(P,1:3) = DES_VEL_NEW(P,1:3)
           RO_ME(P) = RO_Sol(P)
           DES_RADIUS_ME(P) = DES_RADIUS(P)
           PMASS_ME(P) = RO_ME(P) * ACOS(-1.0) * DES_RADIUS_ME(P)**2
           !PRINT*,PMASS_ME(P),RO_ME(P),ACOS(-1.0) * DES_RADIUS_ME(P)**2
           !STOP
           FC_OLD_ME(P,1) = GRAV(1) * PMASS_ME(P)
           FC_OLD_ME(P,2) = -10 * PMASS_ME(P)
           FC_OLD_ME(P,3) = GRAV(3) * PMASS_ME(P)
           
           FC_NEW_ME(P,1) = GRAV(1) * PMASS_ME(P)
           FC_NEW_ME(P,2) = -10 * PMASS_ME(P)
           FC_NEW_ME(P,3) = GRAV(3) * PMASS_ME(P)

           ! PRINT*,P,DES_POS_OLD(P,1:3)
           ! PRINT*,RO_ME(P)
           ! PRINT*,DES_RADIUS_ME(P)
           ! PRINT*,PMASS_ME(P)
           ! PRINT*,FC_OLD_ME(P,2)
                       
        END DO
                                        
        !> 4 particules fictives (2D)

        ! ========================== !        
        !               4            !
        !              ___           !
        !          1  |__|  3        !
        !                            !
        !              2             !
        ! ========================== !
        
        !> MUR 1
        PARTICLES = PARTICLES + 1
        PIP = PIP + 1
        DES_RADIUS(PARTICLES) = 1e8
        DES_POS_NEW(PARTICLES,1) = X_MIN - DES_RADIUS(PARTICLES)
        DES_POS_NEW(PARTICLES,2) = (Y_MIN + Y_MAX)*0.5
        DES_POS_NEW(PARTICLES,3) = (Z_MIN + Z_MAX)*0.5
        DES_VEL_NEW(PARTICLES,1:3) = 0.0
        PMASS_ME(PARTICLES) = 1e12
        RO_Sol(PARTICLES) = 1e12
        PARTICLE_STATE(PARTICLES) = NORMAL_PARTICLE
        
        !> MUR 2
        PARTICLES = PARTICLES + 1
        PIP = PIP + 1
        DES_RADIUS(PARTICLES) = 1e8
        DES_POS_NEW(PARTICLES,1) = (X_MIN + X_MAX)*0.5
        DES_POS_NEW(PARTICLES,2) = Y_MIN - DES_RADIUS(PARTICLES)
        DES_POS_NEW(PARTICLES,3) = (Z_MIN + Z_MAX)*0.5
        DES_VEL_NEW(PARTICLES,1:3) = 0.0
        PMASS_ME(PARTICLES) = 1e12
        RO_Sol(PARTICLES) = 1e12
        PARTICLE_STATE(PARTICLES) = NORMAL_PARTICLE
        
        !> MUR 3
        PARTICLES = PARTICLES + 1
        PIP = PIP + 1
        DES_RADIUS(PARTICLES) = 1e8
        DES_POS_NEW(PARTICLES,1) = X_MAX + DES_RADIUS(PARTICLES)
        DES_POS_NEW(PARTICLES,2) = (Y_MIN + Y_MAX)*0.5
        DES_POS_NEW(PARTICLES,3) = (Z_MIN + Z_MAX)*0.5
        DES_VEL_NEW(PARTICLES,1:3) = 0.0
        PMASS_ME(PARTICLES) = 1e12
        RO_Sol(PARTICLES) = 1e12
        PARTICLE_STATE(PARTICLES) = NORMAL_PARTICLE
        
        !> MUR 4
        PARTICLES = PARTICLES + 1
        PIP = PIP + 1
        DES_RADIUS(PARTICLES) = 1e8
        DES_POS_NEW(PARTICLES,1) = (X_MIN + X_MAX)*0.5
        DES_POS_NEW(PARTICLES,2) = Y_MAX +  DES_RADIUS(PARTICLES) 
        DES_POS_NEW(PARTICLES,3) = (Z_MIN + Z_MAX)*0.5
        DES_VEL_NEW(PARTICLES,1:3) = 0.0
        PMASS_ME(PARTICLES) = 1e12
        RO_Sol(PARTICLES) = 1e12
        PARTICLE_STATE(PARTICLES) = NORMAL_PARTICLE

        !> 6 en 3D
        
        PARTICLES = PARTICLES + 1
        PIP = PIP + 1
        DES_RADIUS(PARTICLES) = 1e8
        DES_POS_NEW(PARTICLES,1) = (X_MIN + X_MAX)*0.5
        DES_POS_NEW(PARTICLES,2) = (Y_MIN + Y_MAX)*0.5
        DES_POS_NEW(PARTICLES,3) = Z_MIN - DES_RADIUS(PARTICLES)
        DES_VEL_NEW(PARTICLES,1:3) = 0.0
        PMASS_ME(PARTICLES) = 1e12
        RO_Sol(PARTICLES) = 1e12
        PARTICLE_STATE(PARTICLES) = NORMAL_PARTICLE
        
        
        PARTICLES = PARTICLES + 1
        PIP = PIP + 1
        DES_RADIUS(PARTICLES) = 1e8
        DES_POS_NEW(PARTICLES,1) = (X_MIN + X_MAX)*0.5
        DES_POS_NEW(PARTICLES,2) = (Y_MIN + Y_MAX)*0.5
        DES_POS_NEW(PARTICLES,3) = Z_MAX + DES_RADIUS(PARTICLES)
        DES_VEL_NEW(PARTICLES,1:3) = 0.0
        PMASS_ME(PARTICLES) = 1e12
        RO_Sol(PARTICLES) = 1e12
        PARTICLE_STATE(PARTICLES) = NORMAL_PARTICLE
        
        
        DO P = 1, PARTICLES            
           DES_POS_OLD(P,1:3) = DES_POS_NEW(P,1:3)
           DES_VEL_OLD(P,1:3) = DES_VEL_NEW(P,1:3)
           DES_RADIUS_ME(P) = DES_RADIUS(P) 
           !PMASS_ME(P) =  PMASS(P)
           
        END DO
      
        CALL OUTPUT_MANAGER(.FALSE., .FALSE.)

        OPEN (UNIT=20, FILE='exact_free_falling_particle.dat')
        OPEN (UNIT=1024, FILE='pos.dat')
       
      END SUBROUTINE DES_TIME_INIT


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!     Subroutine: DES_TIME_STEP                                        !
!     Author: SOUFIANE                                Date: 21-Jun-04  !
!                                                                      !
!     Purpose: Main DEM driver routine                                 !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DES_TIME_STEP(NN)

        use compar, only: ADJUST_PARTITION
        use geometry, only: X_MIN, X_MAX, Y_MIN, Y_MAX, Z_MIN, Z_MAX
        
        IMPLICIT NONE

         INTEGER, INTENT(IN) :: NN
         DOUBLE PRECISION :: X, Y, EC
         

         DO P = 1, PARTICLES
            DES_POS_OLD_ME(P,1:3) = DES_POS_OLD(P,1:3) 
            DES_VEL_OLD_ME(P,1:3) = DES_VEL_OLD(P,1:3) 
            
            DES_POS_NEW_ME(P,1:3) = DES_POS_NEW(P,1:3)
            DES_VEL_NEW_ME(P,1:3) = DES_VEL_NEW(P,1:3)
            
         end DO

         DO P = 1, PARTICLES            
            !> Q_{K+1/2} =  Q_{K} + 0.5 DT \DOT{Q}_{K} 
            DES_POS_DEMI_ME(P,1:3) = DES_POS_OLD_ME(P,1:3) + 0.5*DT*DES_VEL_OLD_ME(P,1:3)
            
            !> THETA-SCHEMA 0.5                      
            DES_VEL_FREE_ME(P,1:3) = DES_VEL_OLD_ME(P,1:3) + 0.5*DT*(FC_OLD_ME(P,1:3)+FC_NEW_ME(P,1:3))/PMASS_ME(P)            
            DES_POS_FREE_ME(P,1:3) = DES_POS_OLD_ME(P,1:3) + 0.5*DT*(DES_VEL_FREE_ME(P,1:3)+DES_VEL_OLD_ME(P,1:3))
            
            !> DOT{Q}_{0}{K+1} =  DOT{Q}_{FREE}{K}
            DES_VEL_NEW_ME(P,1:3) = DES_VEL_FREE_ME(P,1:3)

         END DO
         
         !> NEIGHBOR SEARCH
         CALL NEIGHBOUR_ME
         
         !> INITIALIZATION OF ACTIVE SET PARAMETERS
         PN = 0
         VM = 0
         DPN = 0
         ACTIVE = .FALSE.
         
         !> NLGS LOOP
         DO NLGS = 1, 10

            NB_CONTACTS = 0
            DO P = 1, PARTICLES - 6 !> because of walls
               CC_START = 1
               IF (P.GT.1) CC_START = NEIGHBOR_INDEX_ME(P-1)
               CC_END = NEIGHBOR_INDEX_ME(P)
               NB_CONTACTS = NB_CONTACTS + CC_END - CC_START
               
               !> CONTACTS LOOP
               DO CC = CC_START, CC_END-1
                  I = NEIGHBORS_ME(CC)
                  IF(IS_NONEXISTENT(I)) CYCLE
                  
                  ! ==================================================================== !               
                  ! ================== CONTACT CONDITION + ACTIVE SET ================== !      
                  ! ==================================================================== !                        
                  DX = DES_POS_DEMI_ME(I,1) - DES_POS_DEMI_ME(P,1)
                  DY = DES_POS_DEMI_ME(I,2) - DES_POS_DEMI_ME(P,2)
                  DZ = DES_POS_DEMI_ME(I,3) - DES_POS_DEMI_ME(P,3)
                                    
                  !NEW VARIABLE S = SQRT(DX**2 + DY**2 + DZ**2)
                  NX = DX / SQRT(DX**2 + DY**2 + DZ**2)
                  NY = DY / SQRT(DX**2 + DY**2 + DZ**2)
                  NZ = DZ / SQRT(DX**2 + DY**2 + DZ**2)
                                    
                  S = SQRT(DX**2+DY**2+DZ**2)-(DES_RADIUS_ME(P)+DES_RADIUS_ME(I))
                  
                  IF (S<0) THEN     ! CONTACT CONDITION
                     
                     UN_OLD = &
                          (DES_VEL_OLD_ME(I,1)-DES_VEL_OLD_ME(P,1))*NX + &
                          (DES_VEL_OLD_ME(I,2)-DES_VEL_OLD_ME(P,2))*NY + &
                          (DES_VEL_OLD_ME(I,3)-DES_VEL_OLD_ME(P,3))*NZ
                        
                     UN_NEW = &
                          (DES_VEL_NEW_ME(I,1)-DES_VEL_NEW_ME(P,1))*NX + &
                          (DES_VEL_NEW_ME(I,2)-DES_VEL_NEW_ME(P,2))*NY + &
                          (DES_VEL_NEW_ME(I,3)-DES_VEL_NEW_ME(P,3))*NZ
                     
                     VM(CC) = (UN_NEW + UN_OLD * EN) / (1 + EN)
                     
                     MEQ = PMASS_ME(P) * PMASS_ME(I) / (PMASS_ME(P) + PMASS_ME(I))
                     
                     GAMMA = 1E3                                          
                     TAU_N = PN(CC) - GAMMA * VM(CC) * MEQ / DT
                     
                     !> ACTIVE SET
                     IF (TAU_N>0) THEN
                        DPN(CC) = - VM(CC) * MEQ / DT !!* EN
                        ACTIVE(CC) = .TRUE.
                        
                     ELSE
                        DPN(CC) = 0
                     END IF
                     
                     DES_VEL_NEW_ME(P,1)=DES_VEL_NEW_ME(P,1)+(-DPN(CC))*DT/PMASS_ME(P)*NX!/(1+EN)
                     DES_VEL_NEW_ME(P,2)=DES_VEL_NEW_ME(P,2)+(-DPN(CC))*DT/PMASS_ME(P)*NY!/(1+EN)
                     DES_VEL_NEW_ME(P,3)=DES_VEL_NEW_ME(P,3)+(-DPN(CC))*DT/PMASS_ME(P)*NZ!/(1+EN)
                     
                     DES_VEL_NEW_ME(I,1)=DES_VEL_NEW_ME(I,1)+(DPN(CC))*DT/PMASS_ME(I)*NX!/(1+EN)
                     DES_VEL_NEW_ME(I,2)=DES_VEL_NEW_ME(I,2)+(DPN(CC))*DT/PMASS_ME(I)*NY!/(1+EN)
                     DES_VEL_NEW_ME(I,3)=DES_VEL_NEW_ME(I,3)+(DPN(CC))*DT/PMASS_ME(I)*NZ!/(1+EN)
                     
                     PN(CC) = PN(CC) + DPN(CC)
                     !PRINT*,PN(CC)
                                          
                     DES_POS_FREE_ME(P,1)=DES_POS_OLD_ME(P,1) + DT*DES_VEL_NEW_ME(P,1)
                     DES_POS_FREE_ME(P,2)=DES_POS_OLD_ME(P,2) + DT*DES_VEL_NEW_ME(P,2)
                     DES_POS_FREE_ME(P,3)=DES_POS_OLD_ME(P,3) + DT*DES_VEL_NEW_ME(P,3)
                     
                     DES_POS_FREE_ME(I,1)=DES_POS_OLD_ME(I,1) + DT*DES_VEL_NEW_ME(I,1)
                     DES_POS_FREE_ME(I,2)=DES_POS_OLD_ME(I,2) + DT*DES_VEL_NEW_ME(I,2)
                     DES_POS_FREE_ME(I,3)=DES_POS_OLD_ME(I,3) + DT*DES_VEL_NEW_ME(I,3)
                                          
                  ELSE
                     PN(CC) = 0
                  END IF
                     
               END DO
               
               IF ((CC_END - CC_START) == 0) THEN
                  !PRINT*,"! ==================================================== !"
                  !PRINT*,"! CONNECTIVITE P -- P: ",P,"<==>",P,"  !"
                  !PRINT*,"! ==================================================== !"
                  EXIT
               END IF
               
            END DO
            
            !> CONDITION TO EXIT NLGS LOOP
            IF (MAXVAL(ABS(DT*DPN(1:NEIGHBOR_INDEX_ME(PARTICLES)-1))) .LT. 1E-16) EXIT             
         END DO
                  
         DO P = 1, PARTICLES
            
            DES_POS_NEW_ME(P,1:3) = DES_POS_DEMI_ME(P,1:3) + 0.5*DT*DES_VEL_NEW_ME(P,1:3)
            
            !> NEW ME VARIABLES --> OLD ME VARIABLES 
            DES_POS_OLD_ME(P,1:3) = DES_POS_NEW_ME(P,1:3)  
            DES_VEL_OLD_ME(P,1:3) = DES_VEL_NEW_ME(P,1:3)
         END DO
         
         
         
         !> OLD ME VARIABLES --> MFIX VARIABLES 
         DO P = 1, PARTICLES
            
            DES_POS_OLD(P,1:3) = DES_POS_OLD_ME(P,1:3)
            DES_VEL_OLD(P,1:3) = DES_VEL_OLD_ME(P,1:3)
            
            DES_POS_NEW(P,1:3) = DES_POS_NEW_ME(P,1:3)
            DES_VEL_NEW(P,1:3) = DES_VEL_NEW_ME(P,1:3)
            
         END DO

         !> KINETIC ENERGY
         EC = 0
         DO P = 1, PARTICLES
            EC = EC + 0.5*PMASS_ME(P)*(DES_VEL_NEW(P,1)**2+DES_VEL_NEW(P,2)**2)
         END DO
      
         WRITE(1024,*)"ZONE"
         DO P = 1, PARTICLES - 6
            PRINT*,TIME,NSTEP,DES_POS_NEW(1,2),DES_POS_NEW(2,2),DES_POS_NEW(3,2)
            WRITE(1024,'(10(e15.8,1x))'),DES_POS_NEW(P,1:3),DES_VEL_NEW(P,1:3),DES_RADIUS(P)
         END DO

         ! ==================================================================== !        
         ! ==================================================================== !        
         ! ==================================================================== !    

         DO P = 1, PARTICLES - 6
            EX = EXACT_FREEFALLINGPARTICLE(S_TIME+DT)
            WRITE(20,'(10(E15.8,1X))')S_TIME+DT,DES_POS_NEW(P,2),EX(1),ABS(DES_POS_NEW(P,2)-EX(1))
            IF (S_TIME+DT>=0.34*10) EXIT
         END DO
         
         DTSOLID = DT
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
         
         ! IF(NN == 500) THEN
         !    CALL DES_TIME_END
         ! END IF         
         
       END SUBROUTINE DES_TIME_STEP
    
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!     Subroutine: DES_TIME_END                                         !
!     Author: SOUFIANE                                Date: 21-Jun-04  !
!                                                                      !
!     Purpose: Main DEM driver routine                                 !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
       SUBROUTINE DES_TIME_END
         USE DISCRETELEMENT
         IMPLICIT NONE
         
         PRINT*,""
         PRINT*,""
         PRINT*,"! ======================================== ! "
         PRINT*,"! =========> END DEM SIMULATIONS <======== ! "
         PRINT*,"! ======================================== ! "
                  
         ! Reset the discrete time step to original value.
         DTSOLID = DTSOLID_TMP
         
         !CALL MPI_FINALIZE(IERR)
         !STOP
         
       END SUBROUTINE DES_TIME_END
       
     END MODULE DES_TIME_MARCH
     





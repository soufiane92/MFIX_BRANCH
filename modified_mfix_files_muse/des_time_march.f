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
          RADIUS = MIN(RADIUS,DES_RADIUS_OLD_ME(P))
       END DO
       
       EXIT_LOOP = .FALSE.
       
       
! Initialize time stepping variable for pure granular simulations.

       FACTOR = CEILING(real((TSTOP-TIME)/DTSOLID))
       
       DT = 0.1 * RADIUS
       EN = 0.9
       
       !PRINT*,DT,EN
       
       !CALL MPI_FINALIZE(IERR)
       !STOP
       
       CALL EXPORT_ME('OUT.DAT',DES_POS_OLD_ME,DES_VEL_OLD_ME,DES_RADIUS_OLD_ME,NP)
       
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
         
         PRINT*,"DES_TIME_STEP_ME", NN         
         
         !PRINT*,PARTICLES
         
         DO P = 1, PARTICLES
            CALL PARTICLE_COPY_NEW_OLD_ME
            CALL PARTICLE_COPY_FREE_OLD_ME
            CALL PARTICLE_COPY_DEMI_OLD_ME
         END DO
         
         DO P = 1, PARTICLES
            
            !> Q_{K+1/2} =  Q_{K} + 0.5 DT \DOT{Q}_{K} 
            DES_POS_DEMI_ME(P,1) = DES_POS_OLD_ME(P,1) + 0.5*DT*DES_VEL_OLD_ME(P,1)
            DES_POS_DEMI_ME(P,2) = DES_POS_OLD_ME(P,2) + 0.5*DT*DES_VEL_OLD_ME(P,2)
            DES_POS_DEMI_ME(P,3) = DES_POS_OLD_ME(P,3) + 0.5*DT*DES_VEL_OLD_ME(P,3)
            
            !> THETA-SCHEMA 0.5
            DES_VEL_FREE_ME(P,1) = DES_VEL_OLD_ME(P,1) + 0.5*DT*(FC_OLD_ME(P,1)+&
                 FC_NEW_ME(P,1)/PMASS_NEW_ME(P))
            DES_VEL_FREE_ME(P,2) = DES_VEL_OLD_ME(P,2) + 0.5*DT*(FC_OLD_ME(P,2)+&
                 FC_NEW_ME(P,2)/PMASS_NEW_ME(P))
            DES_VEL_FREE_ME(P,3) = DES_VEL_OLD_ME(P,3) + 0.5*DT*(FC_OLD_ME(P,3)+&
                 FC_NEW_ME(P,3)/PMASS_NEW_ME(P))
            
            DES_POS_FREE_ME(P,1) = DES_POS_OLD_ME(P,1) + 0.5*DT*(DES_VEL_FREE_ME(P,1)&
                 +DES_VEL_OLD_ME(P,1))
            DES_POS_FREE_ME(P,2) = DES_POS_OLD_ME(P,2) + 0.5*DT*(DES_VEL_FREE_ME(P,2)&
                 +DES_VEL_OLD_ME(P,2))
            DES_POS_FREE_ME(P,3) = DES_POS_OLD_ME(P,3) + 0.5*DT*(DES_VEL_FREE_ME(P,3)&
                 +DES_VEL_OLD_ME(P,3))
            
            !> DOT{Q}_{0}{K+1} =  DOT{Q}_{FREE}{K}
            DES_VEL_NEW_ME(P,1) = DES_VEL_FREE_ME(P,1)
            DES_VEL_NEW_ME(P,2) = DES_VEL_FREE_ME(P,2)
            DES_VEL_NEW_ME(P,3) = DES_VEL_FREE_ME(P,3)
            
         END DO
         
! ========================================================================= !         
! ========================================================================= !         
! ========================================================================= !         

         ! ! Calculate inter particle forces acting (collisional, cohesion)
         ! CALL CALC_FORCE_DEM
         
         ! ! Update the old values of particle position and velocity with the new
         ! ! values computed
         ! IF (DO_OLD) CALL CFUPDATEOLD
         
         ! ! Call user functions.
         
         ! !IF(CALL_USR) CALL USR1_DES
         
         ! ! Update position and velocities
         ! CALL CFNEWVALUES

! ========================================================================= !         
! ========================================================================= !         
! ========================================================================= !         
         
         NB_CONTACTS = 0

         !PRINT*,PARTICLES

         !PRINT*,"NB CONTACTS 1",NB_CONTACTS
         !PRINT*,"NEIGHBORS INDEX 1",NEIGHBOR_INDEX(1),NEIGHBOR_INDEX(2)
         
         !CALL MPI_FINALIZE(IERR)
         !STOP
         
         ! DO P = 1, PARTICLES
         !    CC_START = 1
         !    IF (P.GT.1) CC_START = NEIGHBOR_INDEX_ME(P-1)
         !    CC_END = NEIGHBOR_INDEX_ME(P)
         !    NB_CONTACTS = NB_CONTACTS + CC_END - CC_START
         ! END DO

         !> NEIGHBOR SEARCH
         CALL NEIGHBOUR_ME

         !PRINT*,"NB CONTACTS 2",NB_CONTACTS
         !PRINT*,"NEIGHBORS INDEX 2",NEIGHBOR_INDEX(1),NEIGHBOR_INDEX(2)
         
         !CALL MPI_FINALIZE(IERR)
         !STOP
         
         PN = 0
         VM = 0
         DPN = 0
         ACTIVE = .FALSE.

         !PRINT*,ACTIVE

         !> CONNECTIVITY
         NB_CONTACTS = 0
         
         DO P = 1, PARTICLES
            !DO NLGS = 1, 10
               CC_START = 1
               IF (P.GT.1) CC_START = NEIGHBOR_INDEX_ME(P-1)
               CC_END = NEIGHBOR_INDEX_ME(P)
               NB_CONTACTS = NB_CONTACTS + CC_END - CC_START

               DO CC = CC_START, CC_END-1
                  I = NEIGHBORS_ME(CC)
                  IF(IS_NONEXISTENT(I)) CYCLE

                  PRINT*,"NOMBRE DE CONTACTS: ",NB_CONTACTS
                  PRINT*,"CONNECTIVITE P -- I: ",P," <==> ",I
               
! ==================================================================== !               
! ================== CONTACT CONDITION + ACTIVE SET ================== !      
! ==================================================================== !               
                  !DO C = 1, NB_CONTACTS
                  
                  !END DO
! ==================================================================== !                 
! ==================================================================== !               
! ==================================================================== !               
               END DO
               IF ((CC_END-CC_START)==0) PRINT*,"CONNECTIVITE P -- P: ",P," <==> ",P
               
            !END DO
         END DO

         CALL MPI_FINALIZE(IERR)
         STOP
         
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
     

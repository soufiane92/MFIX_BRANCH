
PROGRAM MAIN
  USE M_PARTICLE_BASE
  IMPLICIT NONE
  
  INTEGER :: NP
  REAL(KIND=8),DIMENSION(2) :: XMIN,XMAX
  REAL(KIND=8) :: RADIUS
  REAL(KIND=8) :: DT
  REAL(KIND=8) :: EN
  INTEGER :: IT
  
  TYPE(PARTICLE),DIMENSION(:),ALLOCATABLE:: P_NEW,P_OLD
  TYPE(PARTICLE),DIMENSION(:),ALLOCATABLE:: P_FREE,P_DEMI
  TYPE(PARTICLE),DIMENSION(:),ALLOCATABLE:: P_NLGS
  
  INTEGER :: NB_CONTACTS
  INTEGER ,DIMENSION(:,:),ALLOCATABLE :: LIST_OF_CONTACTS
  INTEGER :: C,P,Q
  
  REAL(KIND=8) :: UN_NEW,UT_NEW
  REAL(KIND=8) :: UN_OLD,UT_OLD
  REAL(KIND=8) :: VM_N,VM_T

  LOGICAL ,DIMENSION(:) , ALLOCATABLE :: ACTIVE
  REAL(KIND=8) ,DIMENSION(:) , ALLOCATABLE :: PN,VM,DPN
  
  REAL(KIND=8) :: RX,RY,TAU_N,GAMMA,WNN,MEQ
  REAL(KIND=8) :: S,SP,SQ,PI

  REAL(KIND=8) :: XC,YC
  REAL(KIND=8) :: NX,NY
  REAL(KIND=8) :: DX,DY
  REAL(KIND=8) :: G
  
  REAL(KIND=8) :: XCP,XCQ
  REAL(KIND=8) :: YCP,YCQ
  REAL(KIND=8) :: C_DOT_N,CP_DOT_N,CQ_DOT_N

  INTEGER :: NL_IT

  REAL(KIND=8) :: DMOY,RES,DMAX,TC,EX(2)


  INTEGER           :: TEST_CASE
  INTEGER,PARAMETER :: FREE_FALLING_PARTICLE = 101
  INTEGER,PARAMETER :: PARTICLE_COLLISION    = 102
  INTEGER,PARAMETER :: DAM_BREAK             = 103
  INTEGER,PARAMETER :: TWO_PARTICLES           = 104
  
  XMIN=0
  XMAX=1

  !PRINT*,XMIN,XMAX
  !STOP
  
  TEST_CASE = TWO_PARTICLES
  !TEST_CASE = FREE_FALLING_PARTICLE
  !TEST_CASE = PARTICLE_COLLISION
  !TEST_CASE = DAM_BREAK

  !CALL INITIALIZE_FREEFALLINGPARTICLE(P_OLD,NP)
  
  IF      (TEST_CASE==FREE_FALLING_PARTICLE) THEN
     CALL INITIALIZE_FREEFALLINGPARTICLE(P_OLD,NP)
     
  ELSE IF (TEST_CASE==TWO_PARTICLES) THEN
     CALL INITIALIZE_2_PARTICLE(P_OLD,NP)   
     
  ELSE IF (TEST_CASE==DAM_BREAK) THEN
     CALL INITIALIZE_DAM(P_OLD,NP)
  END IF
  
  
  ALLOCATE( P_NEW (NP) )
  ALLOCATE( P_FREE(NP) )
  ALLOCATE( P_DEMI(NP) )
  
  ALLOCATE( LIST_OF_CONTACTS(100*NP,2) )
  ALLOCATE( ACTIVE( 100*NP ))

  ALLOCATE( PN( 100*NP ))
  ALLOCATE(DPN( 100*NP ))
  ALLOCATE( VM( 100*NP ))

  !> JE CHERCHE LE RAYON LE PLUS PETIT
  RADIUS = HUGE(0D0)
  DO P=1,NP-4
     RADIUS=MIN(RADIUS,P_OLD(P)%R)
  END DO
  
  !PRINT*,RADIUS,DT
  !STOP
  
  !> ON FIXE LE PAS DE TEMPS ET LE COEFFICIENT 
  !> DE RESTITUTION EN
  DT = 0.1*RADIUS
  EN = 0.9
  
  !PRINT*,DT,EN
  !STOP

  IF  (TEST_CASE==FREE_FALLING_PARTICLE) THEN
     ! 100 PAS DE TEMPS AVANT DE TOMBER EXACTEMENT SUR LA PAROI
     DT = SQRT(2*(0.5-RADIUS)/10.)/100.
     DT = DT/2
     
  END IF
  !PRINT*,NP
  CALL EXPORT('OUT.DAT',P_OLD,NP-4)
  !STOP
  
  !> INITIALISATION DES LISTES P_NEW,P_FREE ET P_DEMI
  !> À PARTIR DE P_OLD
  DO P=1,NP
     CALL PARTICLE_COPY(P_NEW(P),P_OLD(P))
     !PRINT*,P_NEW(1),P_OLD(1)
     !STOP
     CALL PARTICLE_COPY(P_FREE(P),P_OLD(P))
     CALL PARTICLE_COPY(P_DEMI(P),P_OLD(P))
  END DO
  
  TC=0D0
  !> INITIAL KINETIC ENERGY
  DO IT=1,2
     
     PRINT*,"ITÉRATION: ",IT,"TIME: ",TC
     
     TC = TC + DT
     DO P=1,NP

        !> VOIR AVEC SLIDE 47 DE DUBOIS !
        !> Q_{K+1/2} =  Q_{K} + 0.5 DT \DOT{Q}_{K} 
        P_DEMI(P)%X = P_OLD(P)%X + 0.5*DT*P_OLD(P)%U
        P_DEMI(P)%Y = P_OLD(P)%Y + 0.5*DT*P_OLD(P)%V
        
        !PRINT*,P_FREE(P)%U
        !PRINT*,P_FREE(P)%V
        
        !> THETA-SCHEMA 0.5
        P_FREE(P)%U = P_OLD(P)%U + DT*0.5*(P_OLD(P)%FX+P_NEW(P)%FX)/P_NEW(P)%M
        P_FREE(P)%V = P_OLD(P)%V + DT*0.5*(P_OLD(P)%FY+P_NEW(P)%FY)/P_NEW(P)%M
        P_FREE(P)%X = P_OLD(P)%X + DT*0.5*(P_FREE(P)%U+P_OLD(P)%U)
        P_FREE(P)%Y = P_OLD(P)%Y + DT*0.5*(P_FREE(P)%V+P_OLD(P)%V)
        
        !PRINT*,P_FREE(P)%U
        !PRINT*,P_FREE(P)%V

        !> THETA-SCHEMA 1.0
        !P_FREE(P)%U = P_OLD(P)%U + DT*1.0*(P_NEW(P)%FX)/P_NEW(P)%M
        !P_FREE(P)%V = P_OLD(P)%V + DT*1.0*(P_NEW(P)%FY)/P_NEW(P)%M

        !> DOT{Q}_{0}{K+1} =  DOT{Q}_{FREE}{K}
        P_NEW(P)%U = P_FREE(P)%U
        P_NEW(P)%V = P_FREE(P)%V

        !PRINT*,P_NEW(P)%U
        !PRINT*,P_NEW(P)%V
        
     END DO

     !STOP

     !> RECHERCHE DES CONTACTS À PARTIR DE Q_{K+1/2} 
     CALL UPDATE_CONTACT_LIST(P_DEMI,NP,LIST_OF_CONTACTS,NB_CONTACTS&
          &,RADIUS)
!!$
     !STOP
     
     PN = 0
     VM = 0
     DPN = 0
     ACTIVE = .FALSE.
     
     !> ITERATION NLGS
     DO NL_IT=1,3
        RES = 0
        DO C=1,NB_CONTACTS
           P = LIST_OF_CONTACTS(C,1)
           Q = LIST_OF_CONTACTS(C,2)
           !>
           DX = (P_DEMI(Q)%X-P_DEMI(P)%X)
           DY = (P_DEMI(Q)%Y-P_DEMI(P)%Y)
           
           NX = DX/SQRT(DX**2+DY**2)
           NY = DY/SQRT(DX**2+DY**2)
                      

           ! ATTENTION L'ESTIMATION DISTANCE N'EST PAS
           ! SUR UN CRANK-NICHOLSON
           !DX = (P_FREE(Q)%X-P_FREE(P)%X)
           !DY = (P_FREE(Q)%Y-P_FREE(P)%Y)
           DX = (P_DEMI(Q)%X-P_DEMI(P)%X)
           DY = (P_DEMI(Q)%Y-P_DEMI(P)%Y)

           S = SQRT(DX**2+DY**2)-(P_DEMI(P)%R+P_DEMI(Q)%R)
           
           PRINT*,DX
           !STOP

           PRINT*,NX
           
           PRINT*,S
           
           IF ( S<0 ) THEN   ! CONTACT CONDITION 

              UN_OLD = (P_OLD(Q)%U-P_OLD(P)%U)*NX + (P_OLD(Q)%V-P_OLD(P)%V)*NY
              UN_NEW = (P_NEW(Q)%U-P_NEW(P)%U)*NX + (P_NEW(Q)%V-P_NEW(P)%V)*NY

              PRINT*,UN_OLD
              PRINT*,UN_NEW
              
              VM(C) = (UN_NEW + UN_OLD*EN)/(1+EN)

              PRINT*,VM(C)
           
              MEQ = 2*P_NEW(P)%M*P_NEW(Q)%M/(P_NEW(P)%M + P_NEW(Q)%M) 

              PRINT*,MEQ
           
              GAMMA=1E3
              TAU_N = PN(C)-GAMMA*VM(C)*MEQ/DT

              PRINT*,TAU_N
           
              !> ACTIVE SET
              IF (TAU_N>0) THEN
                 DPN(C) = - VM(C)*MEQ/DT*EN
                 ACTIVE(C) = .TRUE.
              ELSE
                 DPN(C) = 0
                 !ACTIVE(C) = .FALSE.
              END IF

              PRINT*,DPN(C)
           
              PRINT*,"1 VEL NEW P: ",P_NEW(P)%U
              PRINT*,"1 VEL NEW I: ",P_NEW(Q)%U
              PRINT*,((-DPN(C))*DT/P_NEW(P)%M*NX)
              P_NEW(P)%U = P_NEW(P)%U + (-DPN(C))*DT/P_NEW(P)%M*NX!/(1+EN)
              P_NEW(P)%V = P_NEW(P)%V + (-DPN(C))*DT/P_NEW(P)%M*NY!/(1+EN)
              
              P_NEW(Q)%U = P_NEW(Q)%U + (+DPN(C))*DT/P_NEW(Q)%M*NX!/(1+EN)
              P_NEW(Q)%V = P_NEW(Q)%V + (+DPN(C))*DT/P_NEW(Q)%M*NY!/(1+EN)

              PRINT*,"2 VEL NEW P: ",P_NEW(P)%U
              PRINT*,"2 VEL NEW I: ",P_NEW(Q)%U
              !STOP
              PN(C)=PN(C)+DPN(C)

              !P_FREE(P)%X = P_OLD(P)%X + DT*0.5*(P_NEW(P)%U + P_OLD(P)%U )
              !P_FREE(P)%Y = P_OLD(P)%Y + DT*0.5*(P_NEW(P)%V + P_OLD(P)%V )
              
              P_FREE(P)%X = P_OLD(P)%X + DT*P_NEW(P)%U 
              P_FREE(P)%Y = P_OLD(P)%Y + DT*P_NEW(P)%V

              P_FREE(Q)%X = P_OLD(Q)%X + DT*P_NEW(Q)%U 
              P_FREE(Q)%Y = P_OLD(Q)%Y + DT*P_NEW(Q)%V 

              PRINT*,PN(C)
              
           ELSE   ! S > 0
              PN(C) = 0
           END IF
           !PRINT*,">",S,DT*DPN(1),P_NEW(P)%V,P_OLD(P)%V
           !PRINT*,">",S,DT*DPN(1),P_NEW(P)%V,P_OLD(P)%V
        END DO
        
        RES= MAXVAL(ABS(DPN(1:NB_CONTACTS)))
        IF (MAXVAL(ABS(DT*DPN(1:NB_CONTACTS))).LT.1E-16) EXIT ! CONDITION TO EXIT NLGS
        !P_NEW(1)%V=0.30991902E+01
     END DO
     PRINT*,NL_IT,MAXVAL(ABS(DT*DPN(1:NB_CONTACTS)))    
     
     DO P=1,NP

        !> COMME DANS LE PAPIER DE SERGE ET MIKAEL
        P_NEW(P)%X = P_DEMI(P)%X + 0.5*DT*P_NEW(P)%U
        P_NEW(P)%Y = P_DEMI(P)%Y + 0.5*DT*P_NEW(P)%V 

        !> JE DECENTRE EN TEMP POUR LE CHOC
!!$        IF (NB_CONTACTS==0) THEN
!!$           P_NEW(P)%X = P_OLD(P)%X + DT*0.5*(P_NEW(P)%U + P_OLD(P)%U )
!!$           P_NEW(P)%Y = P_OLD(P)%Y + DT*0.5*(P_NEW(P)%V + P_OLD(P)%V )
!!$        ELSE
!!$           P_NEW(P)%X = P_OLD(P)%X + DT*1.0*(P_NEW(P)%U )
!!$           P_NEW(P)%Y = P_OLD(P)%Y + DT*1.0*(P_NEW(P)%V)
!!$        END IF

        !> LE THETA-SCHEMA
!!$        P_NEW(P)%X = P_OLD(P)%X + DT*0.5*(P_NEW(P)%U + P_OLD(P)%U )
!!$        P_NEW(P)%Y = P_OLD(P)%Y + DT*0.5*(P_NEW(P)%V + P_OLD(P)%V )
        
        CALL PARTICLE_COPY(P_OLD(P),P_NEW(P))
     END DO
     !<> CALCUL DE DMAX SUR LES CONNECTIONS DE TYPE ACTIVE
     
     
     !PRINT'(I5,1X,I5,1X,10(E15.8,1X))',IT,NL_IT,TC,P_NEW(1)%V,EX(2),P_NEW(1)%Y,EX(1)
     !PRINT'(I5,1X,I5,1X,10(E15.8,1X))',IT,NL_IT,P_NEW(1)%U,P_NEW(2)%U
     
     
     
     IF  (TEST_CASE==FREE_FALLING_PARTICLE) THEN
        
        EX = EXACT_FREEFALLINGPARTICLE(TC)
        PRINT'(I5,1X,I2,1X,10(E15.8,1X))',IT,NB_CONTACTS,P_NEW(1)%Y,EX(1),P_NEW(1)%V,EX(2)
        WRITE(20,'(10(E15.8,1X))')TC,P_NEW(1)%Y,EX(1),ABS(P_NEW(1)%Y-EX(1))
        IF (TC>=0.34*10) EXIT

     ELSE IF  (TEST_CASE==TWO_PARTICLES) THEN
        WRITE(30,'(10(E15.8,1X))')TC,P_NEW(1)%X,P_NEW(1)%Y,P_NEW(1)%U,P_NEW(1)%V
        
     ELSE IF  (TEST_CASE==DAM_BREAK) THEN
        
        PRINT'(I5,1X,I2,1X,10(E15.8,1X))',IT,NB_CONTACTS
        IF (TC>=5.) EXIT
        
     END IF
        
     
     
     
     IF (MOD(IT,10)==0) THEN
        !PN = DT*PN
        CALL EXPORT_REACTION(P_NEW,NP,LIST_OF_CONTACTS,NB_CONTACTS,PN)
        CALL EXPORT('OUT.DAT',P_NEW,NP-4,IT)
     END IF
     
  END DO

! ==================================== CONTAINS ======================================= !
  
CONTAINS
  
! =================================================================================== ! 
  SUBROUTINE SET_BOUNDARY_WALL(PCL,NP)
    INTEGER :: NP
    TYPE(PARTICLE),DIMENSION(:),ALLOCATABLE:: PCL
    
    ! ========================================== !
    ! PRINT*,"2",NP ====== > "ICI NP = 1"
    ! ========================================== !
   
    
    !> FOUR ACTIVE SET WALL => PARTICLE WITH HUGE RADIUS
    PCL(NP+1)%X = XMIN(1)-1E8
    PCL(NP+1)%Y = (XMIN(2)+XMAX(2))*0.5
    PCL(NP+1)%U = 0
    PCL(NP+1)%V = 0
    PCL(NP+1)%R = 1E8
    PCL(NP+1)%RHO = 1E12
    PCL(NP+1)%FX = 0.0
    PCL(NP+1)%FY = 0.0

    
    PCL(NP+2)%X = (XMIN(1)+XMAX(1))*0.5
    PCL(NP+2)%Y = XMIN(2)-1E8
    PCL(NP+2)%U = 0
    PCL(NP+2)%V = 0
    PCL(NP+2)%R = 1E8
    PCL(NP+2)%RHO = 1E12
    PCL(NP+2)%M = 1E12
    PCL(NP+2)%FX = 0.0
    PCL(NP+2)%FY = 0.0
    
    PCL(NP+3)%X = XMAX(1)+1E8
    PCL(NP+3)%Y = (XMIN(2)+XMAX(2))*0.5
    PCL(NP+3)%U = 0
    PCL(NP+3)%V = 0
    PCL(NP+3)%R = 1E8
    PCL(NP+3)%RHO = 1E12
    PCL(NP+3)%FX = 0.0
    PCL(NP+3)%FY = 0.0
    
    PCL(NP+4)%X = (XMIN(1)+XMAX(1))*0.5
    PCL(NP+4)%Y = XMAX(2)+1E8
    PCL(NP+4)%U = 0
    PCL(NP+4)%V = 0
    PCL(NP+4)%R = 1E8
    PCL(NP+4)%RHO = 1E12
    PCL(NP+4)%FX = 0.0
    PCL(NP+4)%FY = 0.0
    !>
    
  END SUBROUTINE SET_BOUNDARY_WALL
! =================================================================================== !  

! =================================================================================== !  
  SUBROUTINE INITIALIZE_STACKEDPARTICLES(PCL,NP)
    IMPLICIT NONE
    INTEGER :: NP
    TYPE(PARTICLE),DIMENSION(:),ALLOCATABLE:: PCL

    INTEGER :: N
    INTEGER :: I,J,P

    REAL(KIND=8) :: XC,YC,U,V,R1,R2,NX,NY
    REAL(KIND=8) :: RADIUS
    
    NP = 3
    ALLOCATE( PCL (NP+4) )
    
    !> 1 SECOND OF SIMULATION
    RADIUS = 1E-2
    XC = 0.5
    YC = 0.5
    U = 1.00
    V = 0.00
    R1 = RADIUS
    R2 = RADIUS

    NX = U/SQRT(U**2+V**2)
    NY = V/SQRT(U**2+V**2)
    
    !> PARTICULE DE GAUCHE 
    PCL(1)%X = XC-NX*0.5
    PCL(1)%Y = YC-NY*0.5
    PCL(1)%U = U
    PCL(1)%V = V
    PCL(1)%R = RADIUS
    PCL(1)%RHO = 2600
    PCL(1)%FX = 10.0*NX
    PCL(1)%FY = 10.0*NY
    
    !> PARTICULE DE DROITE
    PCL(2)%X = XC+NX*0.5
    PCL(2)%Y = YC+NY*0.5
    PCL(2)%U = -U
    PCL(2)%V = -V
    PCL(2)%R = RADIUS
    PCL(2)%RHO = 2600
    PCL(2)%FX = -10.0*NX
    PCL(2)%FY = -10.0*NY
    

    
    
    DO P=1,(NP-2)
       PCL(2+P)%X = XC - (P-1)*2*RADIUS*NX
       PCL(2+P)%Y = YC - (P-1)*2*RADIUS*NY
       PCL(2+P)%U = 0
       PCL(2+P)%V = 0
       PCL(2+P)%R = RADIUS
       PCL(2+P)%RHO = 2600
       PCL(2+P)%FX = 0.0
       PCL(2+P)%FY = 0.0
    END DO
    
    
    CALL SET_BOUNDARY_WALL(PCL,NP)
    
    NP = NP+4
    PI = ACOS(-1.0)
    DO P=1,NP-4
       PCL(P)%M= PCL(P)%RHO * PI*PCL(P)%R**2
       PCL(P)%FX = PCL(P)%FX*PCL(P)%M
       PCL(P)%FY = PCL(P)%FY*PCL(P)%M
    END DO
    
  END SUBROUTINE INITIALIZE_STACKEDPARTICLES
! =================================================================================== !  

! =================================================================================== !  
  SUBROUTINE INITIALIZE_FREEFALLINGPARTICLE(PCL,NP)
    IMPLICIT NONE
    INTEGER :: NP
    TYPE(PARTICLE),DIMENSION(:),ALLOCATABLE:: PCL

    INTEGER :: N
    INTEGER :: I,J,P

    REAL(KIND=8) :: XC,YC,U,V,R1,R2,NX,NY
    REAL(KIND=8) :: RADIUS,TC,VC
    
    NP = 1
    ALLOCATE( PCL (NP+4) )
    
    ! ========================================== !
    ! PRINT*,"1",NP ====== > "ICI NP = 1"
    ! ========================================== !
    
    !> 1 SECOND OF SIMULATION
    RADIUS = 1E-2
    XC = 0.5
    YC = 0.5
    PCL(1)%X = XC
    PCL(1)%Y = YC
    PCL(1)%U = 0
    PCL(1)%V = 0
    PCL(1)%R = RADIUS
    PCL(1)%RHO = 2600
    PCL(1)%FX = 0.0
    PCL(1)%FY = -10.0

    !PRINT*,PCL(1)%X
    
    CALL SET_BOUNDARY_WALL(PCL,NP)
    
    NP = NP+4
    
    ! ========================================== !
    ! PRINT*,"3",NP ====== > "ICI NP = 5"
    ! ========================================== !
   
    PI = ACOS(-1.0)
    DO P=1,NP-4
       PCL(P)%M= PCL(P)%RHO * PI*PCL(P)%R**2
       PCL(P)%FX = PCL(P)%FX*PCL(P)%M
       PCL(P)%FY = PCL(P)%FY*PCL(P)%M
    END DO
    
  END SUBROUTINE INITIALIZE_FREEFALLINGPARTICLE
! =================================================================================== !  

! =================================================================================== !  
  SUBROUTINE INITIALIZE_2_PARTICLE(PCL,NP)
    IMPLICIT NONE
    INTEGER :: NP
    TYPE(PARTICLE),DIMENSION(:),ALLOCATABLE:: PCL

    INTEGER :: N
    INTEGER :: I,J,P

    REAL(KIND=8) :: XC,YC,U,V,R1,R2,NX,NY
    REAL(KIND=8) :: RADIUS,TC,VC
    
    NP = 2
    ALLOCATE( PCL (NP+4) )
    
    !> 1 SECOND OF SIMULATION
    RADIUS = 1E-1
    XC = 0.5
    YC = 0.5
    U = 1.00
    V = 0.00
    R1 = RADIUS
    R2 = RADIUS

    NX = U/SQRT(U**2+V**2)
    NY = V/SQRT(U**2+V**2)
    
    !> PARTICULE DE GAUCHE 
    PCL(1)%X = 0.51
    PCL(1)%Y = YC
    PCL(1)%U = U
    PCL(1)%V = V
    PCL(1)%R = RADIUS
    PCL(1)%RHO = 2600
    PCL(1)%FX = 0.0
    PCL(1)%FY = -9.8
    
    !> PARTICULE DE DROITE
    PCL(2)%X = 0.59
    PCL(2)%Y = 0.5
    PCL(2)%U = -U
    PCL(2)%V = V
    PCL(2)%R = RADIUS
    PCL(2)%RHO = 2600
    PCL(2)%FX = 0.0
    PCL(2)%FY = -9.8
    

    
    
    ! DO P=1,(NP-2)
    !    PCL(2+P)%X = XC - (P-1)*2*RADIUS*NX
    !    PCL(2+P)%Y = YC - (P-1)*2*RADIUS*NY
    !    PCL(2+P)%U = 0
    !    PCL(2+P)%V = 0
    !    PCL(2+P)%R = RADIUS
    !    PCL(2+P)%RHO = 2600
    !    PCL(2+P)%FX = 0.0
    !    PCL(2+P)%FY = 0.0
    ! END DO
    
    
    CALL SET_BOUNDARY_WALL(PCL,NP)
    
    NP = NP+4
    PI = ACOS(-1.0)
    DO P=1,NP-4
       PCL(P)%M= PCL(P)%RHO * PI*PCL(P)%R**2
       PCL(P)%FX = PCL(P)%FX*PCL(P)%M
       PCL(P)%FY = PCL(P)%FY*PCL(P)%M
    END DO
    
  END SUBROUTINE INITIALIZE_2_PARTICLE
! =================================================================================== !  

  
! =================================================================================== !  
  SUBROUTINE INITIALIZE_DAM(PCL,NP)
    IMPLICIT NONE
    INTEGER :: NP
    TYPE(PARTICLE),DIMENSION(:),ALLOCATABLE:: PCL

    INTEGER :: N
    INTEGER :: I,J,P

    REAL(KIND=8) :: XC,YC,U,V,R1,R2,NX,NY
    REAL(KIND=8) :: RADIUS,TC,VC
    
    N=20
    
    NP = N*N    
    ALLOCATE( PCL (NP+4) )
    
    !> DISTRIBUTION DES PARTICULES 
    !> SUR LE COTE GAUCHE
    RADIUS = (XMAX(1)/4)/N*0.5
    P=0
    DO I=1,N
       DO J=1,N
          P=P+1
          P_OLD(P)%X= (I-0.5)*(2*RADIUS)
          P_OLD(P)%Y= (J-0.5)*(2*RADIUS)
          P_OLD(P)%U= 0.0
          P_OLD(P)%V= 0.0
          CALL RANDOM_NUMBER(S)

          P_OLD(P)%R= RADIUS*(1-S*0.5)
          P_OLD(P)%RHO = 2600.0
          P_OLD(P)%FX=-0.0
          P_OLD(P)%FY=-10.
       END DO
    END DO

    
    CALL SET_BOUNDARY_WALL(PCL,NP)
    
    NP = NP+4
    PI = ACOS(-1.0)
    DO P=1,NP
       PCL(P)%M= PCL(P)%RHO * PI*PCL(P)%R**2
       PCL(P)%FX = PCL(P)%FX*PCL(P)%M
       PCL(P)%FY = PCL(P)%FY*PCL(P)%M
    END DO
    
    
  END SUBROUTINE INITIALIZE_DAM
! =================================================================================== !  

! =================================================================================== !  
  FUNCTION EXACT_FREEFALLINGPARTICLE(T)
    IMPLICIT NONE
    REAL(KIND=8) EXACT_FREEFALLINGPARTICLE(2)
    REAL(KIND=8) :: T,R,YC,H0,TN,G,E,Y,V,TT,T0,V0,HN
    INTEGER :: NB,N

    
    
    R=1E-2
    G=10
    YC=0.5
    E = EN
    H0 = YC - R
    T0 = -SQRT(2*G*H0)/G 
    
    TT = T0
    DO N=1,100
       TT = TT + SQRT(8*H0/G)*E**(N-1)
       IF (TT>=T) EXIT
    END DO
    TT = TT - SQRT(8*H0/G)*E**(N-1) 
    
    HN = H0*E**(2*(N-1))
    V0 = SQRT(2*G*HN)
    
    TT = T-TT
    Y = TT*(-0.5*G*TT+V0) + R
    V = -1.0*G*TT+V0

!!$    IF (N==2) THEN
!!$       PRINT'(2(E15.8,1X))',SQRT(2*G*HN)/SQRT(2*G*H0)
!!$    END IF

    
    EXACT_FREEFALLINGPARTICLE(1) = Y
    EXACT_FREEFALLINGPARTICLE(2) = V
    
  END FUNCTION EXACT_FREEFALLINGPARTICLE
! =================================================================================== !  
  
  


END PROGRAM MAIN

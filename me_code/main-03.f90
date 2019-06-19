
program main
  use m_particle_base
  implicit none
  
  integer :: np
  real(kind=8),dimension(2) :: xmin,xmax
  real(kind=8) :: radius
  real(kind=8) :: dt
  real(kind=8) :: en
  integer :: it
  
  type(particle),dimension(:),allocatable:: p_new,p_old
  type(particle),dimension(:),allocatable:: p_free,p_demi
  type(particle),dimension(:),allocatable:: p_nlgs
  
  integer :: nb_contacts
  integer ,dimension(:,:),allocatable :: list_of_contacts
  integer :: c,p,q
  
  real(kind=8) :: un_new,ut_new
  real(kind=8) :: un_old,ut_old
  real(kind=8) :: vm_n,vm_t

  logical ,dimension(:) , allocatable :: active
  real(kind=8) ,dimension(:) , allocatable :: pn,vm,dpn
  
  real(kind=8) :: rx,ry,tau_n,gamma,wnn,meq
  real(kind=8) :: s,sp,sq,pi

  real(kind=8) :: xc,yc
  real(kind=8) :: nx,ny
  real(kind=8) :: dx,dy
  real(kind=8) :: g
  
  real(kind=8) :: xcp,xcq
  real(kind=8) :: ycp,ycq
  real(kind=8) :: c_dot_n,cp_dot_n,cq_dot_n

  integer :: nl_it

  real(kind=8) :: dmoy,res,dmax,tc,ex(2)


  integer           :: test_case
  integer,parameter :: free_falling_particle = 101
  integer,parameter :: particle_collision    = 102
  integer,parameter :: dam_break             = 103
  
  xmin=0
  xmax=1

  !print*,xmin,xmax
  !stop
  
  !test_case = free_falling_particle
  !test_case = particle_collision
  test_case = dam_break

  !print*,test_case
  !stop

  
  !> 1:np-4
  !> np-5:np
  !call Initialize_FreeFallingParticle(p_old,np)
  
  ! ========================================== !
  ! print*,"0",np ====== > "ICI np = 16777216"
  ! ========================================== !
  
  if      (test_case==free_falling_particle) then
     call Initialize_FreeFallingParticle(p_old,np)

     ! ========================================== !
     ! print*,"4",np ====== > "ICI np = 5"
     ! ========================================== !

     
  else if (test_case==dam_break) then
     call Initialize_DAM(p_old,np)
  end if
  
  
  allocate( p_new (np) )
  allocate( p_free(np) )
  allocate( p_demi(np) )
  
  allocate( list_of_contacts(100*np,2) )
  allocate( active( 100*np ))

  allocate( pn( 100*np ))
  allocate(dpn( 100*np ))
  allocate( vm( 100*np ))

  !print*,np
  !print*,p_old(1)
  !stop
  
  !> je cherche le rayon le plus petit
  radius = huge(0d0)
  do p=1,np-4
     radius=min(radius,p_old(p)%r)
  end do
  
  !print*,radius
  !stop
  
  !> on fixe le pas de temps et le coefficient 
  !> de restitution en
  dt = 0.1*radius
  en = 0.9
  
  !print*,dt,en
  !stop

  if  (test_case==free_falling_particle) then
     ! 100 pas de temps avant de tomber exactement sur la paroi
     dt = sqrt(2*(0.5-radius)/10.)/100.
     dt = dt/2
     
  end if
  !print*,np
  call export('out.dat',p_old,np-4)
  !stop
  
  !> initialisation des listes p_new,p_free et p_demi
  !> à partir de p_old
  do p=1,np
     call particle_copy(p_new(p),p_old(p))
     !print*,p_new(1),p_old(1)
     !stop
     call particle_copy(p_free(p),p_old(p))
     call particle_copy(p_demi(p),p_old(p))
  end do
  
  tc=0d0
  !> initial kinetic energy
  do it=1,1000
     !print*,dt,np
     !stop
     tc = tc + dt
     do p=1,np

        !> voir avec slide 47 de Dubois !
        !> q_{k+1/2} =  q_{k} + 0.5 dt \dot{q}_{k} 
        p_demi(p)%x = p_old(p)%x + 0.5*dt*p_old(p)%u
        p_demi(p)%y = p_old(p)%y + 0.5*dt*p_old(p)%v
        
        
        !> theta-schema 0.5
        p_free(p)%u = p_old(p)%u + dt*0.5*(p_old(p)%fx+p_new(p)%fx)/p_new(p)%m
        p_free(p)%v = p_old(p)%v + dt*0.5*(p_old(p)%fy+p_new(p)%fy)/p_new(p)%m
        p_free(p)%x = p_old(p)%x + dt*0.5*(p_free(p)%u+p_old(p)%u)
        p_free(p)%y = p_old(p)%y + dt*0.5*(p_free(p)%v+p_old(p)%v)
        
        !> theta-schema 1.0
        !p_free(p)%u = p_old(p)%u + dt*1.0*(p_new(p)%fx)/p_new(p)%m
        !p_free(p)%v = p_old(p)%v + dt*1.0*(p_new(p)%fy)/p_new(p)%m

        !> dot{q}_{0}{k+1} =  dot{q}_{free}{k}
        p_new(p)%u = p_free(p)%u
        p_new(p)%v = p_free(p)%v
        
     end do
     
     !> recherche des contacts à partir de q_{k+1/2} 
     call update_contact_list(p_demi,np,list_of_contacts,nb_contacts&
          &,radius)
!!$
     !stop
     
     pn = 0
     vm = 0
     dpn = 0
     active = .false.
     
     !> iteration nlgs
     do nl_it=1,1000
        res = 0
        do c=1,nb_contacts
           p = list_of_contacts(c,1)
           q = list_of_contacts(c,2)
           !>
           dx = (p_demi(q)%x-p_demi(p)%x)
           dy = (p_demi(q)%y-p_demi(p)%y)
           
           nx = dx/sqrt(dx**2+dy**2)
           ny = dy/sqrt(dx**2+dy**2)

           ! attention l'estimation distance n'est pas
           ! sur un crank-Nicholson
           !dx = (p_free(q)%x-p_free(p)%x)
           !dy = (p_free(q)%y-p_free(p)%y)
           dx = (p_demi(q)%x-p_demi(p)%x)
           dy = (p_demi(q)%y-p_demi(p)%y)

           s = sqrt(dx**2+dy**2)-(p_demi(p)%r+p_demi(q)%r)
           
           !print*,"s ",nl_it,s

           if ( s<0 ) then   ! Contact condition 

              un_old = (p_old(q)%u-p_old(p)%u)*nx + (p_old(q)%v-p_old(p)%v)*ny
              un_new = (p_new(q)%u-p_new(p)%u)*nx + (p_new(q)%v-p_new(p)%v)*ny
              
              vm(c) = (un_new + un_old*en)/(1+en)
              
              meq = 2*p_new(p)%m*p_new(q)%m/(p_new(p)%m + p_new(q)%m) 
              
              gamma=1e3
              tau_n = pn(c)-gamma*vm(c)*meq/dt

              
              !> active set
              if (tau_n>0) then
                 dpn(c) = - vm(c)*meq/dt*en
                 active(c) = .true.
              else
                 dpn(c) = 0
                 !active(c) = .false.
              end if
              
              p_new(p)%u = p_new(p)%u + (-dpn(c))*dt/p_new(p)%m*nx!/(1+en)
              p_new(p)%v = p_new(p)%v + (-dpn(c))*dt/p_new(p)%m*ny!/(1+en)
              
              p_new(q)%u = p_new(q)%u + (+dpn(c))*dt/p_new(q)%m*nx!/(1+en)
              p_new(q)%v = p_new(q)%v + (+dpn(c))*dt/p_new(q)%m*ny!/(1+en)
              
              pn(c)=pn(c)+dpn(c)

              !p_free(p)%x = p_old(p)%x + dt*0.5*(p_new(p)%u + p_old(p)%u )
              !p_free(p)%y = p_old(p)%y + dt*0.5*(p_new(p)%v + p_old(p)%v )
              
              p_free(p)%x = p_old(p)%x + dt*p_new(p)%u 
              p_free(p)%y = p_old(p)%y + dt*p_new(p)%v

              p_free(q)%x = p_old(q)%x + dt*p_new(q)%u 
              p_free(q)%y = p_old(q)%y + dt*p_new(q)%v 

              
           else   ! s > 0
              pn(c) = 0
           end if
           !print*,">",s,dt*dpn(1),p_new(p)%v,p_old(p)%v
           !print*,">",s,dt*dpn(1),p_new(p)%v,p_old(p)%v
        end do
        
        res= maxval(abs(dpn(1:nb_contacts)))
        if (maxval(abs(dt*dpn(1:nb_contacts))).lt.1e-16) exit ! Condition to exit nlgs
        !p_new(1)%v=0.30991902E+01
     end do
     
     
     do p=1,np

        !> comme dans le papier de Serge et Mikael
        p_new(p)%x = p_demi(p)%x + 0.5*dt*p_new(p)%u
        p_new(p)%y = p_demi(p)%y + 0.5*dt*p_new(p)%v 

        !> je decentre en temp pour le choc
!!$        if (nb_contacts==0) then
!!$           p_new(p)%x = p_old(p)%x + dt*0.5*(p_new(p)%u + p_old(p)%u )
!!$           p_new(p)%y = p_old(p)%y + dt*0.5*(p_new(p)%v + p_old(p)%v )
!!$        else
!!$           p_new(p)%x = p_old(p)%x + dt*1.0*(p_new(p)%u )
!!$           p_new(p)%y = p_old(p)%y + dt*1.0*(p_new(p)%v)
!!$        end if

        !> le theta-schema
!!$        p_new(p)%x = p_old(p)%x + dt*0.5*(p_new(p)%u + p_old(p)%u )
!!$        p_new(p)%y = p_old(p)%y + dt*0.5*(p_new(p)%v + p_old(p)%v )
        
        call particle_copy(p_old(p),p_new(p))
     end do
     !<> calcul de dmax sur les connections de type active
     
     
     !print'(i5,1x,i5,1x,10(e15.8,1x))',it,nl_it,tc,p_new(1)%v,ex(2),p_new(1)%y,ex(1)
     !print'(i5,1x,i5,1x,10(e15.8,1x))',it,nl_it,p_new(1)%u,p_new(2)%u
     
     
     
     if  (test_case==free_falling_particle) then
        
        ex = exact_FreeFallingParticle(tc)
        print'(i5,1x,i2,1x,10(e15.8,1x))',it,nb_contacts,p_new(1)%y,ex(1),p_new(1)%v,ex(2)
        write(20,'(10(e15.8,1x))')tc,p_new(1)%y,ex(1),abs(p_new(1)%y-ex(1))
        if (tc>=0.34*10) exit
        
     else if  (test_case==dam_break) then
        
        print'(i5,1x,i2,1x,10(e15.8,1x))',it,nb_contacts
        if (tc>=5.) exit
        
     end if
        
     
     
     
     if (mod(it,10)==0) then
        !pn = dt*pn
        call export_reaction(p_new,np,list_of_contacts,nb_contacts,pn)
        call export('out.dat',p_new,np-4,it)
     end if
     
  end do

! ==================================== CONTAINS ======================================= !
  
contains
  
! =================================================================================== ! 
  subroutine set_boundary_wall(pcl,np)
    integer :: np
    type(particle),dimension(:),allocatable:: pcl
    
    ! ========================================== !
    ! print*,"2",np ====== > "ICI np = 1"
    ! ========================================== !
   
    
    !> four active set wall => particle with huge radius
    pcl(np+1)%x = xmin(1)-1e8
    pcl(np+1)%y = (xmin(2)+xmax(2))*0.5
    pcl(np+1)%u = 0
    pcl(np+1)%v = 0
    pcl(np+1)%r = 1e8
    pcl(np+1)%rho = 1e12
    pcl(np+1)%fx = 0.0
    pcl(np+1)%fy = 0.0

    
    pcl(np+2)%x = (xmin(1)+xmax(1))*0.5
    pcl(np+2)%y = xmin(2)-1e8
    pcl(np+2)%u = 0
    pcl(np+2)%v = 0
    pcl(np+2)%r = 1e8
    pcl(np+2)%rho = 1e12
    pcl(np+2)%m = 1e12
    pcl(np+2)%fx = 0.0
    pcl(np+2)%fy = 0.0
    
    pcl(np+3)%x = xmax(1)+1e8
    pcl(np+3)%y = (xmin(2)+xmax(2))*0.5
    pcl(np+3)%u = 0
    pcl(np+3)%v = 0
    pcl(np+3)%r = 1e8
    pcl(np+3)%rho = 1e12
    pcl(np+3)%fx = 0.0
    pcl(np+3)%fy = 0.0
    
    pcl(np+4)%x = (xmin(1)+xmax(1))*0.5
    pcl(np+4)%y = xmax(2)+1e8
    pcl(np+4)%u = 0
    pcl(np+4)%v = 0
    pcl(np+4)%r = 1e8
    pcl(np+4)%rho = 1e12
    pcl(np+4)%fx = 0.0
    pcl(np+4)%fy = 0.0
    !>
    
  end subroutine set_boundary_wall
! =================================================================================== !  

! =================================================================================== !  
  subroutine Initialize_StackedParticles(pcl,np)
    implicit none
    integer :: np
    type(particle),dimension(:),allocatable:: pcl

    integer :: n
    integer :: i,j,p

    real(kind=8) :: xc,yc,u,v,r1,r2,nx,ny
    real(kind=8) :: radius
    
    np = 3
    allocate( pcl (np+4) )
    
    !> 1 second of simulation
    radius = 1e-2
    xc = 0.5
    yc = 0.5
    u = 1.00
    v = 0.00
    r1 = radius
    r2 = radius

    nx = u/sqrt(u**2+v**2)
    ny = v/sqrt(u**2+v**2)
    
    !> particule de gauche 
    pcl(1)%x = xc-nx*0.5
    pcl(1)%y = yc-ny*0.5
    pcl(1)%u = u
    pcl(1)%v = v
    pcl(1)%r = radius
    pcl(1)%rho = 2600
    pcl(1)%fx = 10.0*nx
    pcl(1)%fy = 10.0*ny
    
    !> particule de droite
    pcl(2)%x = xc+nx*0.5
    pcl(2)%y = yc+ny*0.5
    pcl(2)%u = -u
    pcl(2)%v = -v
    pcl(2)%r = radius
    pcl(2)%rho = 2600
    pcl(2)%fx = -10.0*nx
    pcl(2)%fy = -10.0*ny
    

    
    
    do p=1,(np-2)
       pcl(2+p)%x = xc - (p-1)*2*radius*nx
       pcl(2+p)%y = yc - (p-1)*2*radius*ny
       pcl(2+p)%u = 0
       pcl(2+p)%v = 0
       pcl(2+p)%r = radius
       pcl(2+p)%rho = 2600
       pcl(2+p)%fx = 0.0
       pcl(2+p)%fy = 0.0
    end do
    
    
    call set_boundary_wall(pcl,np)
    
    np = np+4
    pi = acos(-1.0)
    do p=1,np-4
       pcl(p)%m= pcl(p)%rho * pi*pcl(p)%r**2
       pcl(p)%fx = pcl(p)%fx*pcl(p)%m
       pcl(p)%fy = pcl(p)%fy*pcl(p)%m
    end do
    
  end subroutine Initialize_StackedParticles
! =================================================================================== !  

! =================================================================================== !  
  subroutine Initialize_FreeFallingParticle(pcl,np)
    implicit none
    integer :: np
    type(particle),dimension(:),allocatable:: pcl

    integer :: n
    integer :: i,j,p

    real(kind=8) :: xc,yc,u,v,r1,r2,nx,ny
    real(kind=8) :: radius,tc,vc
    
    np = 1
    allocate( pcl (np+4) )
    
    ! ========================================== !
    ! print*,"1",np ====== > "ICI np = 1"
    ! ========================================== !
    
    !> 1 second of simulation
    radius = 1e-2
    xc = 0.5
    yc = 0.5
    pcl(1)%x = xc
    pcl(1)%y = yc
    pcl(1)%u = 0
    pcl(1)%v = 0
    pcl(1)%r = radius
    pcl(1)%rho = 2600
    pcl(1)%fx = 0.0
    pcl(1)%fy = -10.0

    !print*,pcl(1)%x
    
    call set_boundary_wall(pcl,np)
    
    np = np+4
    
    ! ========================================== !
    ! print*,"3",np ====== > "ICI np = 5"
    ! ========================================== !
   
    pi = acos(-1.0)
    do p=1,np-4
       pcl(p)%m= pcl(p)%rho * pi*pcl(p)%r**2
       pcl(p)%fx = pcl(p)%fx*pcl(p)%m
       pcl(p)%fy = pcl(p)%fy*pcl(p)%m
    end do
    
  end subroutine Initialize_FreeFallingParticle
! =================================================================================== !  

! =================================================================================== !  
  subroutine Initialize_DAM(pcl,np)
    implicit none
    integer :: np
    type(particle),dimension(:),allocatable:: pcl

    integer :: n
    integer :: i,j,p

    real(kind=8) :: xc,yc,u,v,r1,r2,nx,ny
    real(kind=8) :: radius,tc,vc
    
    n=20
    
    np = n*n    
    allocate( pcl (np+4) )
    
    !> distribution des particules 
    !> sur le cote gauche
    radius = (xmax(1)/4)/n*0.5
    p=0
    do i=1,n
       do j=1,n
          p=p+1
          p_old(p)%x= (i-0.5)*(2*radius)
          p_old(p)%y= (j-0.5)*(2*radius)
          p_old(p)%u= 0.0
          p_old(p)%v= 0.0
          call random_number(s)

          p_old(p)%r= radius*(1-s*0.5)
          p_old(p)%rho = 2600.0
          p_old(p)%fx=-0.0
          p_old(p)%fy=-10.
       end do
    end do

    
    call set_boundary_wall(pcl,np)
    
    np = np+4
    pi = acos(-1.0)
    do p=1,np
       pcl(p)%m= pcl(p)%rho * pi*pcl(p)%r**2
       pcl(p)%fx = pcl(p)%fx*pcl(p)%m
       pcl(p)%fy = pcl(p)%fy*pcl(p)%m
    end do
    
    
  end subroutine Initialize_DAM
! =================================================================================== !  

! =================================================================================== !  
  function exact_FreeFallingParticle(t)
    implicit none
    real(kind=8) exact_FreeFallingParticle(2)
    real(kind=8) :: t,r,yc,h0,tn,g,e,y,v,tt,t0,v0,hn
    integer :: nb,n

    
    
    r=1e-2
    g=10
    yc=0.5
    e = en
    h0 = yc - r
    t0 = -sqrt(2*g*h0)/g 
    
    tt = t0
    do n=1,100
       tt = tt + sqrt(8*h0/g)*e**(n-1)
       if (tt>=t) exit
    end do
    tt = tt - sqrt(8*h0/g)*e**(n-1) 
    
    hn = h0*e**(2*(n-1))
    v0 = sqrt(2*g*hn)
    
    tt = t-tt
    y = tt*(-0.5*g*tt+v0) + r
    v = -1.0*g*tt+v0

!!$    if (n==2) then
!!$       print'(2(e15.8,1x))',sqrt(2*g*hn)/sqrt(2*g*h0)
!!$    end if

    
    exact_FreeFallingParticle(1) = y
    exact_FreeFallingParticle(2) = v
    
  end function exact_FreeFallingParticle
! =================================================================================== !  
  
  


end program main

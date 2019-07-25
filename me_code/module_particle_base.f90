module m_particle_base
  implicit none
  type particle
     real(kind=8) :: x=0,y=0
     real(kind=8) :: u=0,v=0
     real(kind=8) :: fx=0,fy=0
     real(kind=8) :: m=1
     real(kind=8) :: rho=1
     real(kind=8) :: r=1
     
     integer bc_x_flag 
     integer bc_y_flag 

     real(kind=8) :: data
     real(kind=8) :: rx=0,ry=0
  end type particle

contains
  subroutine particle_copy(p_old,p_new)
    implicit none
    type(particle)::p_new,p_old
    real(kind=8) :: fx,fy,dt

    !> velocity update
    p_old%u = p_new%u 
    p_old%v = p_new%v
    p_old%x = p_new%x 
    p_old%y = p_new%y

    p_old%fx = p_new%fx 
    p_old%fy = p_new%fy

    p_old%m = p_new%m
    p_old%r = p_new%r
    p_old%rho = p_new%rho

  end subroutine particle_copy

  subroutine particle_zero(p)
    implicit none
    type(particle)::p
    !> velocity update
    p%u = 0.0
    p%v = 0.0
    p%x = 0.0
    p%y = 0.0
    p%fx = 0.0
    p%fy = 0.0
  end subroutine particle_zero
  
  subroutine make_a_time_step(p_old,p_new,dt,fx,fy)
    implicit none
    type(particle)::p_new,p_old
    real(kind=8) :: fx,fy,dt

    !> velocity update
    p_new%u = p_old%u + dt*p_old%fx/p_old%m
    p_new%v = p_old%v + dt*p_old%fy/p_old%m
    !> displacement update
    p_new%x = p_old%x + dt*p_new%u
    p_new%y = p_old%y + dt*p_new%v
    
  end subroutine make_a_time_step

  subroutine export(filename,plist,N,it)
    implicit none
    character(len=*) :: filename
    type(particle),dimension(N) :: plist
    integer,optional :: it
    integer :: N,k
    
    
    if (present(it)) then
       open(unit=24,file=trim(filename),Access = 'append',Status='old')
       write(24,'("ZONE T=""",i5 ,""" ")'),it
    else
       open(unit=24,file=trim(filename))
       write(24,'("ZONE T=""",i5 ,""" ")'),0
    end if
    do k=1,N
       write(24,'(10(e15.8,1x))')plist(k)%x,plist(k)%y,plist(k)%u&
            &,plist(k)%v,plist(k)%r,plist(k)%rx,plist(k)%ry
    end do
    close(24)
    
  end subroutine export

  subroutine update_contact_list(list_of_particles,nb_particles,list_of_contacts,nb_contacts,distance_contact)
    implicit none
 
    integer :: nb_particles,nb_contacts
    type(particle),dimension(:),allocatable::list_of_particles
    integer ,dimension(:,:),allocatable:: list_of_contacts
    real(kind=8) :: distance_contact

    integer :: p,q
    real(kind=8) :: d2,dc2
    real(kind=8) :: xp,xq
    real(kind=8) :: yp,yq
    real(kind=8) :: rq,rp

    real(kind=8) :: xc,yc
    real(kind=8) :: nx,ny
    real(kind=8) :: dx,dy
    real(kind=8) :: xcp,xcq
    real(kind=8) :: ycp,ycq
    real(kind=8) :: c_dot_n,cp_dot_n,cq_dot_n,g
   
    dc2 = (0.1*distance_contact)**2
    !print*,dc2

    list_of_particles(1:nb_particles)%data=0
    nb_contacts=0
    do p = 1 , nb_particles-4
       do q = p+1 , nb_particles
          
          dx = (list_of_particles(q)%x-list_of_particles(p)%x)
          dy = (list_of_particles(q)%y-list_of_particles(p)%y)

!!$          nx = dx/sqrt(dx**2+dy**2)
!!$          ny = dy/sqrt(dx**2+dy**2)
!!$          
!!$          xcp= (list_of_particles(p)%x + list_of_particles(p)%r*nx)
!!$          ycp= (list_of_particles(p)%y + list_of_particles(p)%r*ny)
!!$
!!$          xcq= (list_of_particles(q)%x - list_of_particles(q)%r*nx)
!!$          ycq= (list_of_particles(q)%y - list_of_particles(q)%r*ny)
!!$           
!!$           xc =  0.5*(xcp+xcq)
!!$           yc =  0.5*(ycp+ycq)
!!$           
!!$           c_dot_n  = xc*nx+yc*ny !> projection du centre en local
!!$           cp_dot_n = (xcp*nx+ycp*ny) - c_dot_n
!!$           cq_dot_n = (xcq*nx+ycq*ny) - c_dot_n
!!$          
!!$!           print*,p,q,cp_dot_n**2
!!$           if ((cp_dot_n**2+cq_dot_n**2)<dc2) then
!!$              nb_contacts = nb_contacts + 1
!!$              list_of_contacts(nb_contacts,1) = p
!!$              list_of_contacts(nb_contacts,2) = q
!!$              list_of_particles(p)%data =  p
!!$           end if
           

           g = sqrt(dx**2+dy**2)-(list_of_particles(p)%r+list_of_particles(q)%r)
           if ( g<0 ) then
              nb_contacts = nb_contacts + 1
              list_of_contacts(nb_contacts,1) = p
              list_of_contacts(nb_contacts,2) = q
              list_of_particles(p)%data =  p
           end if



        end do
       
     end do
     
   end subroutine update_contact_list
   
   subroutine get_distance(list_of_particles,nb_particles,list_of_contacts,nb_contacts, active,dmax)
    implicit none
    integer :: nb_particles,nb_contacts
    type(particle),dimension(:),allocatable::list_of_particles
    integer ,dimension(:,:),allocatable:: list_of_contacts
    logical ,dimension(:) , allocatable :: active
    real(kind=8) :: dmax
    

    integer :: c,p,q
    real(kind=8) :: d2,dc2
    real(kind=8) :: xp,xq
    real(kind=8) :: yp,yq
    real(kind=8) :: rq,rp

    real(kind=8) :: xc,yc
    real(kind=8) :: nx,ny
    real(kind=8) :: dx,dy
    real(kind=8) :: xcp,xcq
    real(kind=8) :: ycp,ycq
    
   
    dmax=0
    do c = 1 , nb_contacts

       if (active(c)) then
          p = list_of_contacts(c,1)
          q = list_of_contacts(c,2)
          
          !>
          dx = (list_of_particles(q)%x-list_of_particles(p)%x)
          dy = (list_of_particles(q)%y-list_of_particles(p)%y)
          
          nx = dx/sqrt(dx**2+dy**2)
          ny = dy/sqrt(dx**2+dy**2)
          
          xcp= (list_of_particles(p)%x + list_of_particles(p)%r*nx)
          ycp= (list_of_particles(p)%y + list_of_particles(p)%r*ny)
          
          xcq= (list_of_particles(q)%x - list_of_particles(q)%r*nx)
          ycq= (list_of_particles(q)%y - list_of_particles(q)%r*ny)
          
          dmax = sqrt((xcq-xcp)**2+(ycq-ycp)**2)
          dmax = max(0d0,dmax)
       end if
       
    end do
    
    
  end subroutine get_distance
   

  subroutine export_reaction(particles,nb_p,contact_particles,nb_c,pn)
    implicit none
    integer :: nb_p,nb_c
    type(particle),dimension(:),allocatable::particles
    integer ,dimension(:,:),allocatable:: contact_particles
!    logical ,dimension(:) , allocatable :: active
    real(kind=8) ,dimension(:) , allocatable :: pn
    real(kind=8) :: dmax
    

    integer :: c,p,q
    real(kind=8) :: nx,ny
    real(kind=8) :: dx,dy
    
    dmax=0
    particles(:)%rx=0
    particles(:)%ry=0
    do c = 1 , nb_c
       p = contact_particles(c,1)
       q = contact_particles(c,2)
       !>
       dx = (particles(q)%x-particles(p)%x)
       dy = (particles(q)%y-particles(p)%y)
       
       nx = dx/sqrt(dx**2+dy**2)
       ny = dy/sqrt(dx**2+dy**2)

!       particles(p)%rx = particles(p)%rx - pn(c)*nx
!       particles(p)%ry = particles(p)%ry - pn(c)*ny 
       particles(p)%rx = particles(p)%rx + abs(pn(c))
       particles(q)%ry = 0
       
    end do
    
    
  end subroutine export_reaction

  


   subroutine get_local_vel(pi,pj,pti,ptj,un,ut)
     implicit none
     type(particle) :: pi,pj   !> particule dont on 
     type(particle) :: pti,ptj !> pour le repere local
     real(kind=8)   :: un,ut   !> vitesse relative tangentielle
 
    real(kind=8) :: du,dv
    real(kind=8) :: dx,dy
    real(kind=8) :: nx,ny
    real(kind=8) :: tx,ty
    

    !> vitesse relative repere global
    du = pj%u - pi%u
    dv = pj%v - pi%v

    !> vitesse relative repere local

    dx = ptj%x - pti%x
    dy = ptj%y - pti%y
    !> n(nx,ny) normale a partir de ptj et pti 
    nx = dx/(sqrt(dx**2+dy**2))
    ny = dy/(sqrt(dx**2+dy**2))
    !> tau(tx,ty) tq (n vec tau) > 0 en 2d 
    tx = -ny 
    ty = +nx
    
    !> projection sur la tenge
    un = (du*nx+dv*ny)
    ut = du*tx+dv*ty
    
  end subroutine get_local_vel
  
  
end module m_particle_base


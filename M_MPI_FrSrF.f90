module M_MPI_FrSrf

  use mpi
  use M_MPI_General, only:u,v,x,y,nx,ny,pi,Parallel
  use M_MPI_Exch

  implicit none
  double precision,dimension (:,:),allocatable       ::   phi



  contains


  subroutine Levelset_Ini(yfree,s,e,comm1d, nbrbottom, nbrtop)
    implicit none

    integer                                     ,intent(in)          ::  s,e
    double precision                            ,intent(in)          ::  yfree
    integer                                     ,intent(in)          ::  comm1d, nbrbottom, nbrtop


    double precision                                                 ::  phiges
    integer                                                          ::  i,j,ipr

    allocate( Phi(-1:nx+1,s-2:e+2) )
    do i=-1,nx+1
      do j=s,e
        phi(i,j)=y(j)-yfree
      end do 
    end do

    !! for validation purpose 
    !!!!!!!!!!!!!!!!!! for free surface testing with Wu et al. (2001)
!    do j=s,e 
!      do i=1,nx-1

!        phi(i,j)=10000
!
!        do ipr=1,nx-1
!
!          phiges=sqrt(   ( x(i)-x(ipr) )**2  +  (  y(j)-( 2.57+0.05*cos(pi*(x(ipr)) ) )  )**2   )

!          if (phiges.lt.phi(i,j) ) then 

!            phi(i,j)=phiges

!          end if 

!        end do 

!        if ( y(j).lt.( 2.57+0.05*cos(pi*(x(i)) ) ) ) then 

!          phi(i,j)=-phi(i,j)

!        end if 
!      end do 
!    end do 
!!!!!!!!!!!!!!!!!!!!!!
 call Levelset_Bound(s,e)  !! for y=-1,0,ny, and ny+1 


 if (Parallel) then
   call UVPhi_Exch(Phi, nx, s, e, comm1d, nbrbottom, nbrtop)      
 end if



return 
end subroutine  


  subroutine levelset(dt,s,e,comm1d, nbrbottom, nbrtop)

    implicit none

    integer                                     ,intent(in)          :: s,e
    double precision                            ,intent(in)          :: dt
    integer                                     ,intent(in)          :: comm1d, nbrbottom, nbrtop

    Double precision,dimension (-1:nx+1,s-2:e+2)                     ::   phioldc
    double precision,dimension (0:nx,s-1:e+1)                        ::   dmphix,dpphix,dmphiy,dpphiy
    double precision,dimension (1:nx-1,s:e)                          ::   phix,phiy,Lphin   !! Maybe phix and phiy can be real instead of matrix

    Double precision Lphis


    integer i,j,kk

 !  print*, "I am in level set"

 
  phioldc(:,:)=phi(:,:)
 
!! Second ordeer ENO convective terms for Advection terms of Level set equation !! 
!! Second order in time !! 

  do kk=1,2  !!prediction correction method!!


    do i=0,nx 
      do j=s,e                                
        Dpphix(i,j)=( phi(i+1,j)-phi(i,j)   )/(  x(i+1)-x(i)  )
        Dmphix(i,j)=( phi(i,j)  -phi(i-1,j) )/(  x(i)-x(i-1)  )
        Dpphiy(i,j)=( phi(i,j+1)-phi(i,j)   )/(  y(j+1)-y(j)  )
        Dmphiy(i,j)=( phi(i,j)  -phi(i,j-1) )/(  y(j)-y(j-1)  )
      end do 
    end do

    if (s.eq.1) then
      do i=0,nx
        j=0
        Dpphix(i,j)=( phi(i+1,j)-phi(i,j)   )/(  x(i+1)-x(i)  )
        Dmphix(i,j)=( phi(i,j)  -phi(i-1,j) )/(  x(i)-x(i-1)  )
        Dpphiy(i,j)=( phi(i,j+1)-phi(i,j)   )/(  y(j+1)-y(j)  )
        Dmphiy(i,j)=( phi(i,j)  -phi(i,j-1) )/(  y(j)-y(j-1)  )
      end do 
    end if

    if (e.eq.(ny-1)) then
      do i=0,nx
        j=ny
        Dpphix(i,j)=( phi(i+1,j)-phi(i,j)   )/(  x(i+1)-x(i)  )
        Dmphix(i,j)=( phi(i,j)  -phi(i-1,j) )/(  x(i)-x(i-1)  )
        Dpphiy(i,j)=( phi(i,j+1)-phi(i,j)   )/(  y(j+1)-y(j)  )
        Dmphiy(i,j)=( phi(i,j)  -phi(i,j-1) )/(  y(j)-y(j-1)  )
      end do 
    end if 
 
    if (Parallel) then
      call Dy_A_Exch(Dpphiy,nx, s, e, comm1d, nbrbottom, nbrtop)
      call Dy_A_Exch(Dmphiy,nx, s, e, comm1d, nbrbottom, nbrtop)
    end if
  

    do i=1,nx-1
      do j=s,e    !!A!!

        if (0.5*( u(i,j)+u(i+1,j) ).gt.0.d0) then


          if (  abs(  Dmphix(i,j)-Dmphix(i-1,j) ).lt.abs(  Dpphix(i,j)-Dpphix(i-1,j) )   ) then
            phix(i,j)=Dmphix(i,j)+  0.5*(  Dmphix(i,j)-Dmphix(i-1,j)  ) 
          else
            phix(i,j)=Dmphix(i,j)+  0.5*(  Dpphix(i,j)-Dpphix(i-1,j)  )
          end if 
    
       else

         if (  abs(  Dmphix(i+1,j)-Dmphix(i,j) ).lt.abs(  Dpphix(i+1,j)-Dpphix(i,j) )   ) then
           phix(i,j)=Dpphix(i,j)-  0.5*(  Dmphix(i+1,j)-Dmphix(i,j)  ) 
         else
           phix(i,j)=Dpphix(i,j)-  0.5*(  Dpphix(i+1,j)-Dpphix(i,j)  )
         end if 

       end if 



       if (0.5*( V(i,j)+V(i,j+1) ).gt.0.d0) then


         if (  abs(  DmphiY(i,j)-DmphiY(i,j-1) ).lt.abs(  DpphiY(i,j)-DpphiY(i,j-1) )   ) then
           phiY(i,j)=DmphiY(i,j)+  0.5*(  DmphiY(i,j)-DmphiY(i,j-1)  ) 
         else
           phiY(i,j)=DmphiY(i,j)+  0.5*(  DpphiY(i,j)-DpphiY(i,j-1)  )
         end if 
    
      else

        if (  abs(  DmphiY(i,j+1)-DmphiY(i,j) ).lt.abs(  DpphiY(i,j+1)-DpphiY(i,j) )   ) then
          phiY(i,j)=DpphiY(i,j)-  0.5*(  DmphiY(i,j+1)-DmphiY(i,j)  ) 
        else
          phiY(i,j)=DpphiY(i,j)-  0.5*(  DpphiY(i,j+1)-DpphiY(i,j)  )
        end if 

      end if 

 
    end do 
  end do  !!A!!

  if (kk.eq.1) then

    do i=1,nx-1 
      do j=s,e
 
        Lphin(i,j)=  -0.5*( u(i,j)+u(i+1,j) )* phix(i,j)+&
        &   -0.5*( V(i,j)+V(i,j+1) )* phiy(i,j) 
           
        phi(i,j)=phioldc(i,j)+dt*( -0.5*( u(i,j)+u(i+1,j) )* phix(i,j)+&
        & -0.5*( V(i,j)+V(i,j+1) )* phiy(i,j) )
                
      end do 
    end do 

                  
  end if                        

  if (kk.eq.2) then
   
    do i=1,nx-1 
      do j=s,e 
 
        Lphis=   -0.5*( u(i,j)+u(i+1,j) )* phix(i,j)&
        & -0.5*( V(i,j)+V(i,j+1) )* phiy(i,j) 
            
        phi(i,j)=phioldc(i,j)+0.5*dt*(Lphis+Lphin(i,j) )                       
      end do
    end do 

  end if

  call Levelset_Bound(s,e)   

  if (Parallel) then
    call UVPhi_Exch(Phi, nx, s, e, comm1d, nbrbottom, nbrtop)      
  end if 


end do  !!prediction corection method!! 
                 

RETURN
END subroutine


  subroutine Reinitialize(eps,dt,s,e,comm1d, nbrbottom, nbrtop)
    implicit none
    integer                                     ,intent(in)          ::  s,e
    double precision                            ,intent(in)          ::  dt,eps
    integer                                     ,intent(in)          ::  comm1d, nbrbottom, nbrtop

    Double precision,dimension(-1:nx+1,s-2:e+2)                      ::  phioldcR
    double precision,dimension(0:nx,s-1:e+1)                         ::  dmphix,dpphix,dmphiy,dpphiy
    double precision,Dimension (1:nx-1,s:e)                          ::  phin,lphin
    Double precision                                                 ::  phixm,phiym,phixp,phiyp,phixr,phiyr
    Double precision                                                 ::  sphi,ss,dtau
    Double precision                                                 ::  lphis,bb,avgh,avgh_cpu
    integer                                                          ::  i,j,kk,m,ierr



!   print*, "I am in Reinitialzation"
! Ref; Journal of computational physics vol. 152 pp-493-516 (1999)

    phioldCR(:,:)=phi(:,:)

!!! for non-uniform mesh !! it should be modified for parallel !!

    avgh_cpu=1000 !! large number!!
    do i=1,nx-1 
      do j=s, e 
        bb=min ((x(i+1)-x(i)),(y(j+1)-y(j)))
        if (bb.lt.avgh_cpu) then
          avgh_cpu=bb
        end if 
      end do 
    end do
   
    if (Parallel) then 
      call MPI_Allreduce( avgh_cpu  ,avgh  ,1 ,MPI_Double_Precision, MPI_Min, comm1d, ierr ) 
    endif

    dtau=0.5*avgh  !! fictious time step !
  
!    print*,"avgh",avgh

    do kk=1,3

      do m=1,2   !!prediction correction method for time step the same as Level set for advection terms  !! 


        do i=0,nx 
          do j=s,e                                
            Dpphix(i,j)=( phi(i+1,j)-phi(i,j)   )/(  x(i+1)-x(i)  )
            Dmphix(i,j)=( phi(i,j)  -phi(i-1,j) )/(  x(i)-x(i-1)  )
            Dpphiy(i,j)=( phi(i,j+1)-phi(i,j)   )/(  y(j+1)-y(j)  )
            Dmphiy(i,j)=( phi(i,j)  -phi(i,j-1) )/(  y(j)-y(j-1)  )
          end do 
        end do

        if (s.eq.1) then
          do i=0,nx
            j=0
            Dpphix(i,j)=( phi(i+1,j)-phi(i,j)   )/(  x(i+1)-x(i)  )
            Dmphix(i,j)=( phi(i,j)  -phi(i-1,j) )/(  x(i)-x(i-1)  )
            Dpphiy(i,j)=( phi(i,j+1)-phi(i,j)   )/(  y(j+1)-y(j)  )
            Dmphiy(i,j)=( phi(i,j)  -phi(i,j-1) )/(  y(j)-y(j-1)  )
          end do 
        end if

        if (e.eq.(ny-1)) then
          do i=0,nx
            j=ny
            Dpphix(i,j)=( phi(i+1,j)-phi(i,j)   )/(  x(i+1)-x(i)  )
            Dmphix(i,j)=( phi(i,j)  -phi(i-1,j) )/(  x(i)-x(i-1)  )
            Dpphiy(i,j)=( phi(i,j+1)-phi(i,j)   )/(  y(j+1)-y(j)  )
            Dmphiy(i,j)=( phi(i,j)  -phi(i,j-1) )/(  y(j)-y(j-1)  )
          end do 
        end if 
 
        if (Parallel) then
          call Dy_A_Exch(Dpphiy,nx, s, e, comm1d, nbrbottom, nbrtop)
          call Dy_A_Exch(Dmphiy,nx, s, e, comm1d, nbrbottom, nbrtop)
        end if


        do i=1,nx-1 
          do j=s,e !!A!!  !!A!!


!sphi=phiold(i,j)/(  dsqrt(phiold(i,j)*phiold(i,j)+h*h)  )
!sphi=sgn(phiold)



        if (phioldCR(i,j).gt.eps) then
          sphi=1.d0
        else if (phioldCR(i,j).lt.-eps) then
          sphi=-1.d0
        else
          sphi=phioldCR(i,j)/eps -(1/pi)*dsin(pi*phioldCR(i,j)/eps) 
        end if 

 

        if (  abs(  Dmphix(i,j)-Dmphix(i-1,j) ).lt.abs(  Dpphix(i,j)-Dpphix(i-1,j) )   ) then
          phixm=Dmphix(i,j)+  0.5*(  Dmphix(i,j)-Dmphix(i-1,j)  ) 
        else
          phixm=Dmphix(i,j)+  0.5*(  Dpphix(i,j)-Dpphix(i-1,j)  )
        end if 
    

        if (  abs(  Dmphix(i+1,j)-Dmphix(i,j) ).lt.abs(  Dpphix(i+1,j)-Dpphix(i,j) )   ) then
          phixp=Dpphix(i,j)-  0.5*(  Dmphix(i+1,j)-Dmphix(i,j)  ) 
        else
          phixp=Dpphix(i,j)-  0.5*(  Dpphix(i+1,j)-Dpphix(i,j)  )
        end if 



        if (  abs(  DmphiY(i,j)-DmphiY(i,j-1) ).lt.abs(  DpphiY(i,j)-DpphiY(i,j-1) )   ) then
          phiYm=DmphiY(i,j)+  0.5*(  DmphiY(i,j)-DmphiY(i,j-1)  ) 
        else
          phiYm=DmphiY(i,j)+  0.5*(  DpphiY(i,j)-DpphiY(i,j-1)  )
        end if 
    


        if (  abs(  DmphiY(i,j+1)-DmphiY(i,j) ).lt.abs(  DpphiY(i,j+1)-DpphiY(i,j) )   ) then
          phiYp=DpphiY(i,j)-  0.5*(  DmphiY(i,j+1)-DmphiY(i,j)  ) 
        else
          phiYp=DpphiY(i,j)-  0.5*(  DpphiY(i,j+1)-DpphiY(i,j)  )
        end if
  
  
   

   
        if (sphi*phixp.ge.0.d0.AND.sphi*phixm.ge.0.d0) then
           phixr=phixm
        else if (sphi*phixp.le.0.d0.AND.sphi*phixm.le.0.d0) then 
           phixr=phixp
        else if (sphi*phixp.gt.0.d0.AND.sphi*phixm.lt.0.d0) then
           phixr=0.0
        else if (sphi*phixp.lt.0.d0.AND.sphi*phixm.gt.0.d0) then
           ss=sphi*( abs(phixp)-abs(phixm) )/(phixp-phixm)
           if (ss.gt.0.0) then
             phixr=phixm
           else
             phixr=phixp 
           end if
        end if
   
   
       if (sphi*phiyp.ge.0.d0.AND.sphi*phiym.ge.0.d0) then
         phiyr=phiym
       else if (sphi*phiyp.le.0.d0.AND.sphi*phiym.le.0.d0) then 
         phiyr=phiyp
       else if (sphi*phiyp.gt.0.d0.AND.sphi*phiym.lt.0.d0) then
         phiyr=0.0
       else if (sphi*phiyp.lt.0.d0.AND.sphi*phiym.gt.0.d0) then
         ss=sphi*( abs(phiyp)-abs(phiym) )/(phiyp-phiym)
         if (ss.gt.0.0) then
           phiyr=phiym
         else
           phiyr=phiyp 
         end if
       end if
    
        
   
   
       if (m.eq.1) then
         lphin(i,j)=sphi*(  1.d0-dsqrt(phiyr*phiyr+phixr*phixr)  )
         phin(i,j)=phi(i,j)
         phi(i,j)=phi(i,j)+dtau*sphi*(  1.d0-dsqrt(phiyr*phiyr+phixr*phixr)  )
       end if
   
   
       if (m.eq.2) then
         lphis   =sphi*(  1.d0-dsqrt(phiyr*phiyr+phixr*phixr)  )
         phi(i,j)=phin(i,j)+0.5*dtau*(  lphis+lphin(i,j)  )                   
       end if

       

     end do 
   end do 

   call Levelset_Bound(s,e) 

   if (Parallel) then
     call UVPhi_Exch(Phi, nx, s, e, comm1d, nbrbottom, nbrtop)      
   end if 

 end do 



end do !!fictious time step!!

return 
END subroutine


subroutine Levelset_Bound(s,e) 
implicit none

integer,intent(in)            :: s,e
integer                       :: i,j


 do j=s,e 
   phi(0,j)=2*phi(1,j)-phi(2,j)
   phi(-1,j)=3*phi(1,j)-2*phi(2,j)

   phi(nx,j)=2*phi(nx-1,j)-phi(nx-2,j)
   phi(nx+1,j)=3*phi(nx-1,j)-2*phi(nx-2,j)
 end do 

 if (s.eq.1) then 
   do i=-1,nx+1 
     phi(i,0)=2*phi(i,1)-phi(i,2)
     phi(i,-1)=3*phi(i,1)-2*phi(i,2)
   end do
 end if 


 if (e.eq.(ny-1)) then 
   do i=-1,nx+1
     phi(i,ny)=2*phi(i,ny-1)-phi(i,ny-2)
     phi(i,ny+1)=3*phi(i,ny-1)-2*phi(i,ny-2)
   end do
 end if 
 

 return 
end subroutine 





end module 



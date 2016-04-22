

program Code_3Phase_2D

!Parallel
use mpi
use M_MPI_General
use M_MPI_Exch 
use M_MPI_Mesh
use M_MPI_Solver  
use M_MPI_FrSrf
use M_MPI_Solid

!Parallel
   
  implicit none
 
  double precision                                   ::    dt,div,divmax,divmax2,beta,maxp,pdif
  double precision                                   ::    ro_w,miu_w,ro_A,miu_A
  double precision                                   ::    ro_s,miu_s,ro_s2,miu_s2
  double precision                                   ::    eps_cpu,eps,Gap,XRef,YRef

  double precision                                   ::    StartTime,EndTime

  double precision,dimension (:,:),allocatable       ::   P
  double precision,dimension (:,:),allocatable       ::   Tx,Ty
  double precision,dimension (:,:),allocatable       ::   advectu,advectv
  double precision,dimension (:,:),allocatable       ::   Amx,Apx,Amy,Apy,Ap,Q

  double precision,dimension (:,:),allocatable       ::   dpvx,dmvx,dpvy,dmvy,dpux,dmux,dpuy,dmuy

  double precision                                   ::   YFree
  double precision                                   ::   Drag_cpu,Lift_cpu
  double precision                                   ::   Drag,Lift,U_Inflow


  double precision,dimension (2)                     ::   ChV
  integer                                            ::   i,j,k,tp,count,ADVECT,plot,viscose,tstep,number2
  integer                                            ::   comm1d,nbrbottom,nbrtop,myid,MyRank,numprocs,ierr 
  integer                                            ::    s,e    !,svs,svl,evs,evl
  integer,dimension(1)                               ::   DIMS
  logical,dimension(1)                               ::   PERIODS

 
  if (Parallel) then 

    call MPI_INIT( ierr )
    call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )

    print *, "my id is",myid,"of",numprocs
    PERIODS=.FALSE.
    DIMS=numprocs

    call MPI_CART_CREATE(MPI_COMM_WORLD,1,DIMS, PERIODS,.true., comm1d, ierr)
    call MPI_COMM_RANK( comm1d, MyRank, ierr )
    call MPI_Cart_shift(comm1d,0,1,nbrbottom,nbrtop,ierr)

  else

    numprocs=1
    MyRank=0
    myid=0

  endif

  if (myid.eq.0) then 
    call CPU_TIME(StartTime)
    print*,"Start time is", StartTime
  end if 

  

  s=1+MyRank*((ny-1)/numprocs)
  e=s+((ny-1)/numprocs)-1
  print*,s,e,"MyRank",MyRank
  call flush (6)

  call General_Ini(s,e)
  count=100000000  !! big number 
  include 'Par_General.txt'


  call MeshGenerator(Chv,s,e,MyRank,comm1d,nbrbottom,nbrtop)
  !print*,"I am out of MeshGen"


  !CellSize(1)=Lx/real(nx-1)
  !CellSize(2)=Ly/real(ny-1)
  !hx(-1:nx+1)    =CellSize(1)
  !hy(s-2:e+2)    =CellSize(2)   

  !x(-1)=-1.5d0*hx(1)
  !x(0) =-0.5*hx(1)
  !x(1)=0.5*hx(1)  
  !Do i=2,nx+1   !! other x(-) and x(0) 
  !  x(i)=x(i-1)+0.5d0*hx(i-1)+0.5d0*hx(i)
  !end do 

    
    
    !y(1)=0.5d0*hy(1)
    !Do j=2,ny

    !  y(j)=y(j-1)+0.5d0*hy(j-1)+0.5d0*hy(j)
    !end do
  !! Uniform grid no Parallelization is required!! 
  ! do j=s-2,e+2
  !   y(j)=real(j)*CellSize(2)-0.5*CellSize(2)
  ! end do 

 

!! end of grid generation !!

!! intial values !! 


  !ro(0:nx,s-1:e+1)=ro_W
  !miuv(0:nx,s-1:e+1)=Miu_W    


  allocate(  p(0:nx,s-1:e+1) )
  allocate(  Tx(1:nx-1,s:e),Ty(1:nx-1,s:e) )
  allocate(  dpux(0:nx,s-1:e+1),dmux(0:nx,s-1:e+1),dpuy(0:nx,s-1:e+1),dmuy(0:nx,s-1:e+1)  )    
  allocate(  dpvx(0:nx,s-1:e+1),dmvx(0:nx,s-1:e+1),dpvy(0:nx,s-1:e+1),dmvy(0:nx,s-1:e+1)  )
  allocate(  advectu(1:nx-1,  s:e),advectv(1:nx-1,s:e)  )
  allocate(  Amx(1:nx-1,s:e),Apx(1:nx-1,s:e),Amy(1:nx-1,s:e),Apy(1:nx-1,s:e),Ap(1:nx-1,s:e),Q(1:nx-1,s:e) )

  p( 0:nx,s-1:e+1)=0


  if (MyRank.eq.0) then  
    OPEN(35,file='result.plt')
    OPEN(75,file='FloatLocation.plt')
    OPEN(85,file='WaveGage.plt')
    OPEN(95,file='Force_Cyl.plt')

    write(35,*) 'variables="x","y","u","v","p","ro","phi","Ib","Vorticity"'
    write(75,*) 'variables="Time (s)","Xbar","YBar","Theta (deg)"'
    write(85,'(A150)') 'variables="time","Gage1","Gage2","Gage3","Gage4","Gage5","Gage6","Gage7","Gage8","Gage9","Gage10","Gage11","Gage12","Gage13","Gage14"'
    write(95,*) 'variables="time","Drag Cof.","Lift Cof."'

  end if 


   do i=1,nx-1

     if (x(i).le.XRef.AND.x(i+1).ge.XRef) then
       gap=0.75*(x(i+1)-x(i))
       exit
      end if
   end do
  
    !! If the free surface is off, it still needs some values, otherwise we encoutre errors in code.
   eps_cpu=0 
   if (int(ny/2).ge.s.AND.int(ny/2).le.e) then 
     eps_cpu=3*( y(int(ny/2))-y(int(ny/2)-1) )
     print*,eps_cpu, "eps equal to"    
   end if
    
   if (parallel) then
     call MPI_Allreduce( eps_CPU ,eps  ,1 ,MPI_Double_Precision, MPI_Sum, comm1d, ierr )
   endif    
   print*, "total eps",eps
   call Levelset_Ini(yfree,s,e,comm1d, nbrbottom, nbrtop)
          
   call Solid_Ini(XRef,YRef,Gap,s,e,comm1d,nbrbottom,nbrtop)
  
   call Boundarycond(0,dt,s,e,comm1d)



   do tp=1,tstep !time step



     print*,tp


     call Prop(ro_W ,ro_A ,ro_s,ro_s2,eps,Gap,ro,s,e)
     call Prop(miu_W,miu_A,miu_s,miu_s2,eps,Gap,miuv,s,e)
       
    
     if (Parallel) then
       call Dy_A_Exch(ro  ,nx, s, e, comm1d, nbrbottom, nbrtop)
       call Dy_A_Exch(miuv,nx, s, e, comm1d, nbrbottom, nbrtop)
     end if


     call Advection_Ini(dpux,dmux,dpuy,dmuy,dpvx,dmvx,dpvy,dmvy,s,e)    
    
     if (Parallel) then
       call Dy_A_Exch(Dpuy,nx, s, e, comm1d, nbrbottom, nbrtop)
       call Dy_A_Exch(Dmuy,nx, s, e, comm1d, nbrbottom, nbrtop)
       call Dy_A_Exch(Dpvy,nx, s, e, comm1d, nbrbottom, nbrtop)
       call Dy_A_Exch(Dmvy,nx, s, e, comm1d, nbrbottom, nbrtop)
     end if

     call Advection (dpux,dmux,dpuy,dmuy,dpvx,dmvx,dpvy,dmvy,advectU,advectV,s,e)   
 
     call Viscosity (Tx,Ty,s,e)

     Do i=1,nx-1 
       Do j=s,e 
         U(i,j)=U(i,j)+advectu(i,j)*dt+Tx(i,j)*dt
       end do 
     end do

     Do j=s,e 
       Do i=1,nx-1
         V(i,j)=V(i,j)+advectv(i,j)*dt+Ty(i,j)*dt+gy*dt
       end do 
     end do 
     
     call Boundarycond_Poi(tp,dt,s,e,comm1d)

     if (Parallel) then
       call UVPhi_Exch(U, nx, s, e, comm1d, nbrbottom, nbrtop)
       call UVPhi_Exch(V, nx, s, e, comm1d, nbrbottom, nbrtop)
     end if
    !call PrintData(u,v,p,ro,x,y,count,plot,MyRank,nx,ny,s,e,comm1d,Parallel,tp,tstep,dt)
     call Poisson_COF(dt,Apx,Amx,Apy,Amy,Ap,Q,s,e)

     if ( Parallel ) then 
       call Poisson_Solver_Parallel(Apx,Amx,Apy,Amy,Ap,Q,pdif,beta,P,s,e,comm1d,nbrbottom,nbrtop)
     else
       call Poisson_Solver(Apx,Amx,Apy,Amy,Ap,Q,pdif,beta,P,s,e)
     end if !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!solve pressure SOR end !!!!!!!

 !!! start update velocity !!!!!!!!!!!!!!!!!!!
     do i=1,nx-1 
       do j=s,e
         U(i,j)=U(i,j)-(  1/(  0.5*(ro(i,j)+ro(i-1,j))  )   )*( p(i,j)-p(i-1,j) )*dt/( 0.5*(hx(i)+hx(i-1))  )
       end do 
     end do 

     do i=1,nx-1 
       do j=s,e 
         V(i,j)=V(i,j)-(  1/(  0.5*(ro(i,j)+ro(i,j-1))  )   )*( p(i,j)-p(i,j-1) )*dt/( 0.5*(hy(j)+hy(j-1))  )
       end do 
     end do 

     call Boundarycond(tp,dt,s,e,comm1d)
! 
     if (Parallel) then
       call UVPhi_Exch(U, nx, s, e, comm1d, nbrbottom, nbrtop)
       call UVPhi_Exch(V, nx, s, e, comm1d, nbrbottom, nbrtop)
     end if

  
     divmax2=0
     do i=1,nx-1 
       do j=s,e     !! need to work !!
         div=(U(i+1,j)-U(i,j))/hx(i)+(v(i,j+1)-v(i,j))/hy(j)

         if(abs(div).gt.divmax2)then
           divmax2=abs(div)
         else
         end if
     end do 
    end do 

    WRITE(*,*)'DIVERGENCE Second order in time =', divmax2,"MyRank",MyRank,"e=",e

    Drag_cpu=0 ; Lift_cpu=0
    do i=2,nx-1
      do j=s,e
        Drag_cpu=Drag_cpu+ 0.5*(ro(i,j)+ro(i-1,j))*0.5*(Ib_1(i,j)+Ib_1(i-1,j))*(x(i)-x(i-1))*( 0.5*(y(j)+y(j+1))-0.5*(y(j)+y(j-1)) )* &
        &          ( +advectu(i,j)-(  1/(  0.5*(ro(i,j)+ro(i-1,j))  )   )*( p(i,j)-p(i-1,j) )/( 0.5*(hx(i)+hx(i-1))  )+Tx(i,j) )

        Lift_cpu=Lift_cpu+ 0.5*(ro(i,j)+ro(i,j-1))*0.5*(Ib_1(i,j)+Ib_1(i,j-1))*(y(j)-y(j-1))*( 0.5*(x(i)+x(i+1))-0.5*(x(i)+x(i-1)) )* &
        &          ( +advectv(i,j)-(  1/(  0.5*(ro(i,j)+ro(i,j-1))  )   )*( p(i,j)-p(i,j-1) )/( 0.5*(hy(j)+hy(j-1))  )+Ty(i,j)+gy )

      end do
    end do

    if (Parallel) then

       call MPI_Allreduce( Drag_cpu   ,Drag  ,1 ,MPI_Double_Precision, MPI_Sum, comm1d, ierr )
       call MPI_Allreduce( Lift_cpu   ,Lift  ,1 ,MPI_Double_Precision, MPI_Sum, comm1d, ierr )

     end if

     call Solid_Update(XRef,YRef,Gap,dt,s,e,comm1d,nbrbottom,nbrtop)



     if (Free_Surface) then 


       call Levelset(dt,s,e,comm1d,nbrbottom,nbrtop)
       call reinitialize(eps,dt,s,e,comm1d,nbrbottom,nbrtop)


     end if 



     if (PrintD) then
       call PrintData(p,count,plot,MyRank,s,e,comm1d,tp,tstep,dt)
     end if

     
     include "Par_UniFlow.txt"
     if (MyRank.eq.0) then
       write (75,'(4(1x,e15.7))') tp*dt,XBar(1),XBar(2),180/pi*Atan( ( points_1(5,2)-points_1(6,2) )/(points_1(5,1)-points_1(6,1))  )
       write (95,'(3(1x,e15.7))') tp*dt,Drag/(0.5*ro_w*(U_inflow**2)*2*r_cyl),Lift/(0.5*ro_w*(U_inFlow**2)*2*r_cyl)
     end if




   end do !time step !! end of the main loop of the code !! 




   if (myid.eq.0) then 
      write(*,*)'end'

      call CPU_TIME(EndTime)
      print*,"Simulation time is",EndTime-StartTime,"seconds"
   end if 
  
   if (Parallel) then
      call MPI_FINALIZE(ierr)
    end if 

   !read (*,*)
   stop





end program 






   subroutine Boundarycond_Poi(tp,dt,s,e,comm1d)
    
    use mpi
    use M_MPI_General, only: u,v,x,y,Landa,gy,pi,nx,ny,Parallel,Uni_Flow,Wave_Gen     
    use M_MPI_FrSrf,   only: phi 
    implicit none
 
    integer,intent(in)                                              ::   s,e,tp
    double precision,intent(in)                                     ::   dt
    integer,intent(in)                                              ::   comm1d

 
    double precision                                                ::   period,hait,HH,HH1,Wwave,Kwave
    double precision                                                ::   sumvv_CPU,sumvv,Yfree_t_CPU,YFree_t,U_Inflow
    integer                                                         ::   i,j,ierr

    

    if (Wave_Gen) then 
      include "Par_Wave_2D.txt"
    else if (Uni_Flow ) then
      include 'Par_UniFlow.txt'
    end if 


    if (Wave_Gen) then

 
      period=landa/(  sqrt(-gy*hait)*sqrt( tanh(2*pi/landa)/(2*pi/landa*hait) )   )
      kwave=2*pi/landa
      wwave=2*pi/period

      !print*, "BC",s,e
      yfree_t_CPU=0
      do j=s,e
        if (0.5*( phi(1,j)+phi(0,j) ).le.0.AND.0.5*( phi(1,j+1)+phi(0,j+1) ).ge.0 ) then             !! may be more precsiness can be done !! 
          yfree_t_CPU=y(j)+ ( y(j+1)-y(j) )/ ( 0.5*( phi(1,j+1)+phi(0,j+1) )-0.5*( phi(1,j)+phi(0,j) ) ) *( -0.5*( phi(1,j)+phi(0,j) ) ) 
          exit 
        end if 
      end do  
      if (Parallel) then
        call MPI_Allreduce( yfree_t_CPU ,YFree_t  ,1 ,MPI_Double_Precision, MPI_Sum, comm1d, ierr )
      end if

      if (YFree_t.eq.0)then 
        print*,"error" 
      end if

      if ( (tp*dt).lt.(1*period) ) then 
        HH=(tp*dt)/(1*period)*HH1
      else 
        HH=HH1
      end if
 
      sumvv_cpu=0
      do j=s ,e 

        if (phi(0,j).lt.0 ) then
          u(1 ,j)=HH/2* wwave *cos(kwave*(0)        -wwave*(tp*dt))*cosH(kwave* ( (y(j)-yfree_t)+hait) )/sinH(kwave *hait)
        else
          u(1 ,j)=HH/2* wwave *cos(kwave*(0)        -wwave*(tp*dt))*cosH(kwave* (-(y(j)-yfree_t)+hait) )/sinH(kwave *hait)
        end if 
        sumvv_cpu=sumvv_cpu+u(1,j)*0.5*(y(j+1)-y(j-1))*1.0d0 

        u(nx,j)=0

      end do 

      if (Parallel) then
        call MPI_Allreduce( sumvv_CPU ,sumvv  ,1 ,MPI_Double_Precision, MPI_Sum, comm1d, ierr )
      end if

      if (s.eq.1) then 
        Do i=1, nx-1 
          V(i,1)=0
        end do 
      end if 
 
      if (e.eq.(ny-1)) then 
        do i=1, nx-1
          V(i,ny)=sumvv/( real(nx-1) *0.5*(x(i+1)-x(i-1)) )   !0
        end do
      end if


    else if (Uni_Flow) then

      
      do j=s,e
        u(1,j)=U_Inflow
        u(nx,j)=u(nx-1,j)
      end do
      
      if      (s==1) then
        do i=1,nx-1
          V(i,1)=0
        end do 
      else if (e==(ny-1)) then
        do i=1,nx-1
          v(i,ny)=0
        end do
      end if  


    end if 

    return 
  end subroutine 







  subroutine Boundarycond(tp,dt,s,e,comm1d)
    
    use mpi
    use M_MPI_General, only: u,v,x,y,Landa,gy,pi,nx,ny,Parallel,Uni_Flow,Wave_Gen     
    use M_MPI_FrSrf,   only: phi 
    implicit none
 
    integer,intent(in)                                              ::   s,e,tp
    double precision,intent(in)                                     ::   dt
    integer,intent(in)                                              ::   comm1d

 
    double precision                                                ::   period,hait,HH,HH1,Wwave,Kwave
    double precision                                                ::   sumvv_CPU,sumvv,Yfree_t_CPU,YFree_t,U_Inflow
    integer                                                         ::   i,j,ierr

    
    if (Wave_Gen) then
      include "Par_Wave_2D.txt"
    elseif (Uni_Flow) then
      include "Par_UniFlow.txt"
    end if 


    if (Wave_Gen) then


      period=landa/(  sqrt(-gy*hait)*sqrt( tanh(2*pi/landa)/(2*pi/landa*hait) )   )
      kwave=2*pi/landa
      wwave=2*pi/period

      !print*, "BC",s,e
      yfree_t_CPU=0
      do j=s,e
        if (0.5*( phi(1,j)+phi(0,j) ).le.0.AND.0.5*( phi(1,j+1)+phi(0,j+1) ).ge.0 ) then             !! may be more precsiness can be done !! 
          yfree_t_CPU=y(j)+ ( y(j+1)-y(j) )/ ( 0.5*( phi(1,j+1)+phi(0,j+1) )-0.5*( phi(1,j)+phi(0,j) ) ) *( -0.5*( phi(1,j)+phi(0,j) ) ) 
          exit 
        end if 
      end do  
      if (Parallel) then
        call MPI_Allreduce( yfree_t_CPU ,YFree_t  ,1 ,MPI_Double_Precision, MPI_Sum, comm1d, ierr )
      end if

      if (YFree_t.eq.0)then 
        print*,"error" 
      end if

      if ( (tp*dt).lt.(1*period) ) then 
        HH=(tp*dt)/(1*period)*HH1
      else 
        HH=HH1
      end if
 
      sumvv_cpu=0
      do j=s ,e 

        if (phi(0,j).lt.0 ) then

          u(1 ,j)=HH/2* wwave *cos(kwave*(0)        -wwave*(tp*dt))*cosH(kwave* ( (y(j)-yfree_t)+hait) )/sinH(kwave *hait)
          u(0 ,j)=HH/2* wwave *cos(kwave*(-2)*x(1)  -wwave*(tp*dt))*cosH(kwave* ( (y(j)-yfree_t)+hait) )/sinH(kwave *hait)
          u(-1,j)=HH/2* wwave *cos(kwave*(-4)*x(1)  -wwave*(tp*dt))*cosH(kwave* ( (y(j)-yfree_t)+hait) )/sinH(kwave *hait)

        else

          u(1 ,j)=HH/2* wwave *cos(kwave*(0)        -wwave*(tp*dt))*cosH(kwave* (-(y(j)-yfree_t)+hait) )/sinH(kwave *hait)
          u(0 ,j)=HH/2* wwave *cos(kwave*(-2)*x(1)  -wwave*(tp*dt))*cosH(kwave* (-(y(j)-yfree_t)+hait) )/sinH(kwave *hait)
          u(-1,j)=HH/2* wwave *cos(kwave*(-4)*x(1)  -wwave*(tp*dt))*cosH(kwave* (-(y(j)-yfree_t)+hait) )/sinH(kwave *hait)

        end if 
        sumvv_cpu=sumvv_cpu+u(1,j)*0.5*(y(j+1)-y(j-1))*1.0d0 

        u(nx,j)=0
        u(nx+1,j)=-u(nx-1,j) 

      end do 
      if (Parallel) then
        call MPI_Allreduce( sumvv_CPU ,sumvv  ,1 ,MPI_Double_Precision, MPI_Sum, comm1d, ierr )
      end if


    elseif (Uni_Flow) then
 

      do j=s ,e 
        u(1,j)   =U_inflow    
        u(0,j)   =U_inflow
        u(-1,j)  =U_inflow

        u(nx,j)  =u(nx-1,j) 
        u(nx+1,j)=u(nx-2,j) 
      end do


    end if



!    do j=s ,e
!      U(1,j)=0      
!      U(0,j)=-U(2,j)
!      U(-1,j)=-U(3,J)   !! New for Parallel 
!      u(nx,j)=0 
!      u(nx+1,j)=-u(nx-1,j) 
!    end do


!! at first smaller and then use them for boundary of bigger !!!
      
    

    if (Uni_Flow.OR.Wave_Gen) then

  
      if (s.eq.1) then
        do i=-1,nx+1 
          U(i,0)=U(i,1)  !! NO slip for cavity 
          U(i,-1)=U(i,2) !! NO slip for cavity 
        end do 
      endif

      if (e.eq.(ny-1)) then
        do i=-1,nx+1
          U(i,ny)=U(i,ny-1)  !! 1.0 
          U(i,ny+1)=U(i,ny-2) ! 1.0
        end do
      end if

 
    end if 


    if (Uni_Flow.OR.Wave_Gen) then


      if (s.eq.1) then 
        Do i=1, nx-1 
          V(i,1)=0
          V(i,0)=-V(i,2)
          V(i,-1)=-V(i,3) 
        end do 
      end if

 
    end if


    if (Wave_Gen) then 


      if (e.eq.(ny-1)) then 
        do i=1, nx-1
          V(i,ny)=sumvv/( real(nx-1) *0.5*(x(i+1)-x(i-1)) )   !0
          V(i,ny+1)=V(i,ny)                                   !-V(i,ny-1)
        end do
      end if


    else if (Uni_Flow) then


      if (e.eq.(ny-1)) then 
        do i=1, nx-1
          V(i,ny)=0         
          V(i,ny+1)=-V(i,ny-1)                                   !-V(i,ny-1)
        end do
      end if

    
    end if 


    if (Uni_Flow.OR.Wave_Gen) then


      do j=s,e 
        V(0,j)=V(1,j)       !! NO slip for cavity 
        V(-1,j)=V(2,j)      !! NO slip for cavity 
        V(nx,j)=V(nx-1,j)   !! NO slip for cavity 
        V(nx+1,j)=V(nx-2,j) !! NO slip for cavity 
      end do

      if (s.eq.1) then 
        do j=-1,0
          V(0,j)=V(1,j)       !! NO slip for cavity 
          V(-1,j)=V(2,j)      !! NO slip for cavity 
          V(nx,j)=V(nx-1,j)   !! NO slip for cavity 
          V(nx+1,j)=V(nx-2,j) !! NO slip for cavity 
        end do 
      end if
   
      if (e.eq.(ny-1)) then    
        do j=ny,ny+1
          V(0,j)=V(1,j)       !! NO slip for cavity 
          V(-1,j)=V(2,j)      !! NO slip for cavity 
          V(nx,j)=V(nx-1,j)   !! NO slip for cavity 
          V(nx+1,j)=V(nx-2,j) !! NO slip for cavity 
        end do 
      end if

 
    end if 
 


    return 
  end subroutine 

 

!!!!!!!!!!!!! smoothed heaviside function !!!!!!!!!!!!!
!!!!!!subroutine for finding weather a point is inside the solid or not !!!


  subroutine Prop(Pro_W,Pro_A,Pro_s,Pro_s2,eps,gap,Pro,s,e)

    use M_MPI_General,   only: nx,ny,Free_Surface
    use M_MPI_FrSrf, only: phi 
    use M_MPI_Solid, only: Ib_1,Ib_2
    implicit none 

    integer,intent(in)                                               ::   s,e
    dOUBLE PRECISION,intent(in)                                      ::   Pro_W,Pro_A,Pro_s,Pro_s2
    dOUBLE PRECISION,intent(in)                                      ::   eps,gap

    double precision,dimension (0:nx,s-1:e+1)   ,intent(out)         ::   pro

    integer i,j
    REAL                                                             ::   HV

!! it should be written more general !!

    if (Free_Surface) then 


      do i=0,nx 
        do j=s,e
          Pro(i,j)=Pro_A    +HV( -phi(i,j),eps )*(Pro_W  -Pro_A    )
          Pro(i,j)=Pro(i,j)+Ib_1(i,j)          *(Pro_s-Pro(i,j) )
          Pro(i,j)=Pro(i,j)+Ib_2(i,j)*(Pro_s2 )
        end do 
      end do
      !! Written for special case. It should be written more general !!

      if (s.eq.1) then
        do i=0,nx
          Pro(i,0)=Pro_W
        end do 
      end if

      if (e.eq.(ny-1)) then
        do i=0,nx
          Pro(i,ny)=Pro_a
        end do 
      end if
   
      
    else


      do i=0,nx 
        do j=s,e
          Pro(i,j)=Pro_W
          Pro(i,j)=Pro(i,j)+Ib_1(i,j)*(Pro_s-Pro(i,j) )
          Pro(i,j)=Pro(i,j)+Ib_2(i,j)*(Pro_s2 )
        end do 
      end do

      if (s.eq.1) then
        do i=0,nx
          Pro(i,0)=Pro_W
        end do 
      end if

      if (e.eq.(ny-1)) then
        do i=0,nx
          Pro(i,ny)=Pro_w
        end do 
      end if

   
    end if


    return 

  end subroutine  

 
function HV(phi_D,eps)

double precision phi_D,pi,eps


pi=3.141592654D0
if (phi_D.gt.eps) then
HV=1.d0
else if (phi_D.lt.-eps) then
HV=0.d0
else
HV=0.5d0*(  1.d0+ phi_D/eps +(1/pi)*dsin(pi*phi_D/eps)  )
end if 

return
end function 



   
  subroutine Advection_Ini(dpux,dmux,dpuy,dmuy,dpvx,dmvx,dpvy,dmvy,s,e)    
    use M_MPI_General, only:  u,v,hx,hy,nx,ny
    implicit none 
    integer                                      ,intent(in)              ::  s,e

    double precision,dimension (0:nx,s-1:e+1),intent(out)            ::   dpux,dmux,dpuy,dmuy  !,advectu
    double precision,dimension (0:nx,s-1:e+1) ,intent(out)           ::   dpvx,dmvx,dpvy,dmvy  !,advectv

    integer i,j,k
    !!!!!!!!!!!! second order ENO method for advection terms !!!!!!!!!!!!
 

    !  advectu(0:nx,s-1,e+1)=0
    !  advectv(0:nx,s-1:e+1)=0 
    do i=0,nx 
      do j=s,e 
        Dpux(i,j)=( u(i+1,j)-u(i,j)   )/hx(i)
        Dmux(i,j)=( u(i,j)  -u(i-1,j) )/hx(i-1)
        Dpuy(i,j)=( u(i,j+1)-u(i,j)   )/( 0.5d0*( hy(j+1)+hy(j)) )
        Dmuy(i,j)=( u(i,j)  -u(i,j-1) )/( 0.5d0*( hy(j)+hy(j-1)) )
      end do 
    end do 

    do i=0,nx 
      do j=s,e 
        Dpvx(i,j)=( v(i+1,j)-v(i,j)   )/( 0.5d0*( hx(i+1)+hx(i)) )
        Dmvx(i,j)=( v(i,j)  -v(i-1,j) )/( 0.5d0*( hx(i)+hx(i-1)) )
        Dpvy(i,j)=( v(i,j+1)-v(i,j)   )/hy(j)
        Dmvy(i,j)=( v(i,j)  -v(i,j-1) )/hy(j-1)
      end do 
    end do 

    if (s.eq.1) then 

      do i=0,nx 
        j=0
        Dpux(i,j)=( u(i+1,j)-u(i,j)   )/hx(i)
        Dmux(i,j)=( u(i,j)  -u(i-1,j) )/hx(i-1)
        Dpuy(i,j)=( u(i,j+1)-u(i,j)   )/( 0.5d0*( hy(j+1)+hy(j)) )
        Dmuy(i,j)=( u(i,j)  -u(i,j-1) )/( 0.5d0*( hy(j)+hy(j-1)) )
      end do 

      do i=0,nx 
        j=0
        Dpvx(i,j)=( v(i+1,j)-v(i,j)   )/( 0.5d0*( hx(i+1)+hx(i)) )
        Dmvx(i,j)=( v(i,j)  -v(i-1,j) )/( 0.5d0*( hx(i)+hx(i-1)) )
        Dpvy(i,j)=( v(i,j+1)-v(i,j)   )/hy(j)
        Dmvy(i,j)=( v(i,j)  -v(i,j-1) )/hy(j-1)
      end do 

    end if 

    if (e.eq.(ny-1)) then 

      do i=0,nx 
        j=ny
        Dpux(i,j)=( u(i+1,j)-u(i,j)   )/hx(i)
        Dmux(i,j)=( u(i,j)  -u(i-1,j) )/hx(i-1)
        Dpuy(i,j)=( u(i,j+1)-u(i,j)   )/( 0.5d0*( hy(j+1)+hy(j)) )
        Dmuy(i,j)=( u(i,j)  -u(i,j-1) )/( 0.5d0*( hy(j)+hy(j-1)) )
      end do 

      do i=0,nx 
        j=ny
        Dpvx(i,j)=( v(i+1,j)-v(i,j)   )/( 0.5d0*( hx(i+1)+hx(i)) )
        Dmvx(i,j)=( v(i,j)  -v(i-1,j) )/( 0.5d0*( hx(i)+hx(i-1)) )
        Dpvy(i,j)=( v(i,j+1)-v(i,j)   )/hy(j)
        Dmvy(i,j)=( v(i,j)  -v(i,j-1) )/hy(j-1)
      end do 

    end if 

    return

  end subroutine 


 !subroutine advection_unifrom(u,v,nx,ny,advectU,advectV,s,e,svS,svl,evl,cellsize)    
  subroutine Advection (dpux,dmux,dpuy,dmuy,dpvx,dmvx,dpvy,dmvy,advectU,advectV,s,e)    
    use M_MPI_General, only: u,v,nx
    implicit none
    integer                                       ,intent(in)             ::  s,e
    double precision,dimension (0:nx,s-1:e+1) ,intent(in)                 ::  dpux,dmux,dpuy,dmuy  !,advectu
    double precision,dimension (0:nx,s-1:e+1) ,intent(in)                 ::  dpvx,dmvx,dpvy,dmvy  !,advectv

    double precision,dimension (1:nx-1,  s:e)     ,intent(out)            ::  advectu
    double precision,dimension (1:nx-1,  S:e)     ,intent(out)            ::  advectv

    double precision                                                      ::  ux,uy,vx,vy
    integer i,j,k
     

    
    do i=1,nx-1 
      do j=s,e    !!A!!
 
 
 !!U1    
        if (u(i,j).gt.0.d0) then

          if (  abs(  Dmux(i,j)-Dmux(i-1,j) ).lt.abs(  Dpux(i,j)-Dpux(i-1,j) )   ) then
            ux=Dmux(i,j)+  0.5*(  Dmux(i,j)-Dmux(i-1,j)  ) 
          else
            ux=Dmux(i,j)+  0.5*(  Dpux(i,j)-Dpux(i-1,j)  )
          end if 
    
        else

          if (  abs(  Dmux(i+1,j)-Dmux(i,j) ).lt.abs(  Dpux(i+1,j)-Dpux(i,j) )   ) then
            ux=Dpux(i,j)-  0.5*(  Dmux(i+1,j)-Dmux(i,j)  ) 
          else
            ux=Dpux(i,j)-  0.5*(  Dpux(i+1,j)-Dpux(i,j)  )
          end if 

        end if 


!!U2
        if (  (0.25*( V(i,j)+V(i,j+1)+v(i-1,j)+v(i-1,j+1) )).gt.0.d0) then

          if (  abs(  DmuY(i,j)-DmuY(i,j-1) ).lt.abs(  DpuY(i,j)-DpuY(i,j-1) )   ) then
            uY=DmuY(i,j)+  0.5*(  DmuY(i,j)-DmuY(i,j-1)  ) 
          else
            uY=DmuY(i,j)+  0.5*(  DpuY(i,j)-DpuY(i,j-1)  )
          end if 
    
        else

          if (  abs(  DmuY(i,j+1)-DmuY(i,j) ).lt.abs(  DpuY(i,j+1)-DpuY(i,j) )   ) then
            uY=DpuY(i,j)-  0.5*(  DmuY(i,j+1)-DmuY(i,j)  ) 
          else
            uY=DpuY(i,j)-  0.5*(  DpuY(i,j+1)-DpuY(i,j)  )
          end if 

       end if 
 
 
 
       advectu(i,j)=-(  u(i,j)*uX+0.25*( V(i,j)+V(i,j+1)+v(i-1,j)+v(i-1,j+1) )*uY  )  !!i-1 j+1 in 2 dimensiona l Parallelization !!


     end do
   end do 
  

 do i=1,nx-1 
   do j=S,e  !!A!!
 
 
!!V1
     if (  (0.25*( u(i,j)+u(i+1,j)+u(i,j-1)+u(i+1,j-1) )).gt.0.d0) then

       if (  abs(  Dmvx(i,j)-Dmvx(i-1,j) ).lt.abs(  Dpvx(i,j)-Dpvx(i-1,j) )   ) then
         vx=Dmvx(i,j)+  0.5*(  Dmvx(i,j)-Dmvx(i-1,j)  ) 
       else
         vx=Dmvx(i,j)+  0.5*(  Dpvx(i,j)-Dpvx(i-1,j)  )
       end if 
    
     else

       if (  abs(  Dmvx(i+1,j)-Dmvx(i,j) ).lt.abs(  Dpvx(i+1,j)-Dpvx(i,j) )   ) then
         vx=Dpvx(i,j)-  0.5*(  Dmvx(i+1,j)-Dmvx(i,j)  ) 
       else
         vx=Dpvx(i,j)-  0.5*(  Dpvx(i+1,j)-Dpvx(i,j)  )
       end if 

     end if 


!!V2
     if ( V(i,j).gt.0.d0) then

       if (  abs(  DmvY(i,j)-DmvY(i,j-1) ).lt.abs(  DpvY(i,j)-DpvY(i,j-1) )   ) then
         vY=DmvY(i,j)+  0.5*(  DmvY(i,j)-DmvY(i,j-1)  ) 
       else
         vY=DmvY(i,j)+  0.5*(  DpvY(i,j)-DpvY(i,j-1)  )
       end if 
    
     else

       if (  abs(  DmvY(i,j+1)-DmvY(i,j) ).lt.abs(  DpvY(i,j+1)-DpvY(i,j) )   ) then
         vY=DpvY(i,j)-  0.5*(  DmvY(i,j+1)-DmvY(i,j)  ) 
       else
         vY=DpvY(i,j)-  0.5*(  DpvY(i,j+1)-DpvY(i,j)  )
       end if 

     end if 
 
 
     advectv(i,j)=-(0.25*(  u(i,j)+u(i+1,j)+u(i,j-1)+u(i+1,j-1)  )*vX+ V(i,j)*vY)   !! i+1,j-1 Parallel 

   end do 
 end do 
 
 
 return 
 
end subroutine 


  



  Subroutine Viscosity (Tx,Ty,s,e)

    use M_MPI_General, only:u,v,hx,hy,nx,ro,miuv
    implicit none 
    integer,intent(in)                                             ::   s,e
    double precision,dimension(1:nx-1,s:e),intent(out)             ::   Tx,Ty 

    dOUBLE PRECISION                                               ::   Txxr,Txxl,Tyxd,Tyxu,Tyyu,Tyyd,Txyr,Txyl
    integer                                                        ::   i,j

  !!!!!!!!!!!! simple centeral difference for viscose terms !!!!!!!!!!!!!!!!!!


  Do i=1,nx-1 
   Do j=s, e 

     Txxr=2*miuv(i,j)  *( U(i+1,j)-U(i,j)   )/hx(i)
 
     Txxl=2*miuv(i-1,j)*( U(i,j)  -U(i-1,j) )/hx(i-1)

     Tyxu=0.25*( miuv(i,j)+miuv(i,j+1)+miuv(i-1,j)+miuv(i-1,j+1) )*&
     &(      ( U(i,j+1)-U(i,j)     )/( 0.5*(hy(j)+hy(j+1)) )+ ( V(i,j+1)-V(i-1,j+1) )/( 0.5*(hx(i)+hx(i-1)) )     )
                 
     Tyxd=0.25*( miuv(i,j-1)+miuv(i,j)+miuv(i-1,j-1)+miuv(i-1,j) )*&
     &(      ( U(i,j)-U(i,j-1)     )/( 0.5*(hy(j)+hy(j-1)) )+ ( V(i,j)  -V(i-1,j)   )/( 0.5*(hx(i)+hx(i-1)) )     )


     Tx(i,j)= (  (Txxr-Txxl)/( 0.5*( hx(i)+hx(i-1)) ) + (Tyxu-Tyxd)/hy(j)  )/(  0.5*(ro(i,j)+ro(i-1,j))  )
 
   
    end do 
  end do   




    Do j=s,e 
      Do i=1,nx-1 


        Tyyu=2*miuv(i,j)  *( V(i,j+1)-V(i,j)   )/hy(j)
 
        Tyyd=2*miuv(i,j-1)*( V(i,j)  -V(i,j-1) )/hy(j-1)
 
 
        Txyr=0.25*( miuv(i,j)+miuv(i,j-1)+miuv(i+1,j)+miuv(i+1,j-1) )*&
        & (     ( V(i+1,j)-V(i,j)     )/( 0.5*(hx(i)+hx(i+1)) )+ ( U(i+1,j)  -U(i+1,j-1) )/( 0.5*(hy(j)+hy(j-1))  )   )
  
        Txyl=0.25*( miuv(i-1,j)+miuv(i-1,j-1)+miuv(i,j)+miuv(i,j-1) )*&
        & (     ( V(i,j)-V(i-1,j)     )/( 0.5*(hx(i)+hx(i-1)) )+ ( U(i,j)    -U(i,j-1)   )/( 0.5*(hy(j)+hy(j-1))  )   )
  
 

       Ty(i,j)= (  (Tyyu-Tyyd)/( 0.5*( hy(j)+hy(j-1)) ) + (Txyr-Txyl)/hx(i) )/(  0.5*(ro(i,j)+ro(i,j-1))  )

 
     end do 
   end do 
 

return 
end 


!!!!!!!!!!!!!! simple iterative method for solving poisson equation !!!!!!!!!!!!
  subroutine Poisson_COF(dt,Apx,Amx,Apy,Amy,Ap,Q,s,e)

    use M_MPI_General, only: u,v,hx,hy,ro,nx,ny
    implicit none

    integer,intent(in)                                                 ::    s,e
    double precision,intent(in)                                        ::    dt

    double precision,dimension (1:nx-1,s:e),intent(out)                ::     Apx,Amx,Apy,Amy,AP,Q
    integer i,j
 


    do i=1,nx-1 
      do j=s,e
                         
 
 if (i==nx-1) then
 Apx(i,j)=0
 else
 
 Apx(i,j)=1/( 0.5d0*(hx(i)+hx(i+1))*hx(i) )/(  0.5*(ro(i,j)+ro(i+1,j))  )
 end if 
 
 if (i==1)then 
 Amx(i,j)=0       
 else 
 
 Amx(i,j)=1/( 0.5d0*(hx(i)+hx(i-1))*hx(i) )/(  0.5*(ro(i,j)+ro(i-1,j))  )
 end if 
 
 if (j==ny-1) then
 Apy(i,j)=0
 else
 
 Apy(i,j)=1/( 0.5d0*(hy(j)+hy(j+1))*hy(j) )/(  0.5*(ro(i,j)+ro(i,j+1))  )

 end if 
 if (j==1) then
 Amy(i,j)=0
 else
 
 Amy(i,j)=1/( 0.5d0*(hy(j)+hy(j-1))*hy(j) )/(  0.5*(ro(i,j)+ro(i,j-1))  )
 
 end if 
 
 
 
 AP(i,j)=-( Apx(i,j)+Amx(i,j)+Apy(i,j)+Amy(i,j) )   
 
 Q(I,J)=(  (U(I+1,J)-U(I,J))/hx(i)+(V(I,J+1)-V(I,J))/hy(j) )/dt
 

 end do ; end do

    return

  end 


  subroutine Poisson_Solver(Apx,Amx,Apy,Amy,Ap,Q,pdif,beta,P,s,e)

    use M_MPI_General, only: nx,ny
    implicit none  
    integer,intent(in)                                                 ::    s,e
    double precision,intent(in)                                        ::    pdif,beta
    double precision,dimension (1:nx-1,s:e)  ,intent(in)               ::    Apx,Amx,Apy,Amy,AP,Q
    double precision,dimension (0:nx,s-1:e+1),intent(inout)            ::    p

    double precision                                                   ::   maxp
    double precision,dimension (0:nx,s-1:e+1)                          ::   pold
    integer i,j,number2

 




      number2=0
      maxp=1.5
      do while (maxp.gt.pdif.AND.number2.lt.100000)

      maxp=0 
      number2=number2+1
 
      do i=1,nx-1
        do j=s,e 

          pold(i,j)=p(i,j)
          p(i,j)=beta*( ( Apx(i,j)*P(I+1,j)+Amx(i,j)*P(I-1,j)+Apy(i,j)*P(I,j+1)+Amy(i,j)*P(I,j-1)-Q(I,j)  )/(-AP(i,j))  ) &
          & +(1-beta)*p(i,j)

          if (abs( pold(i,j)-p(i,j) ).gt.maxp) then 
            maxp=abs( pold(i,j)-p(i,j) )
          end if 
 
        end do 
      end do
 
    end do  !!while !
                        
 
    print *,"pdif=",maxp,"number it=",number2

   
    
    
    do j=s,e
      p(0,j)=p(1,j)
      p(nx,j)=p(nx-1,j)
    end do


    if (s.eq.1) then 
      do i=0,nx
        p(i,0)=p(i,1)
      end do
    end if

    if (e.eq.(ny-1)) then 
      do i=0,nx 
        p(i,ny)=p(i,ny-1)
      end do 
    end if 

  
    return 
  end     
 
 
  


subroutine PrintData(p,count,plot,MyRank,s,e,comm1d,tp,tstep,dt)

!Parallel
  use mpi
  use M_MPI_General, only: u,v,ro,x,y,Landa,Lx,nx,ny,Parallel,Wave_Gen 
  use M_MPI_FrSrf,   only: phi 
  use M_MPI_Solid,   only: Ib_1
!Parallel
  implicit none

  integer,intent(in)                                             :: s,e
  double precision,dimension (0:nx   ,s-1:e+1),intent(in)        :: p
  dOUBLE PRECISION,intent(in)                                    :: dt
  integer,intent(in)                                             :: MyRank,comm1d 
  integer,intent(in)                                             :: tp,tstep,plot
  integer,intent(inout)                                          :: count
  
  double precision,dimension(:,:),allocatable                    :: Vort
  double precision,dimension ((nx-1)*(e-s+1))                    :: Upak,VPak,PPak,roPak,PhiPak,Ib_1Pak,VortPak
  double precision,dimension ((nx-1)*(ny-1))                     :: UTo,VTo,PTo,roTo,PhiTo,Ib_1To,VortTo
  double precision,dimension (1:ny-1)                            :: YTo
  double precision,dimension(14)                                 :: XGage,YFreeG,YFreeG_cpu
  integer                                                        :: i,j,kk,ii,igage
  integer                                                        :: status(MPI_STATUS_SIZE), ierr
  
    include "Par_Print.txt"
    allocate (vort(1:nx-1,s:e))

    count=count+1
    if (count.ge.plot)then


      count=0
      print*,"data is written"
      call flush(6)

      do i=1,nx-1 
        do j=s,e
          vort(i,j)= 0.25*(( v(i,j)    -v(i-1,j)   )/(x(i)-x(i-1))   &
          &               +( v(i+1,j)  -v(i,j)     )/(x(i+1)-x(i))   &
          &               +( v(i,j+1)  -v(i-1,j+1) )/(x(i)-x(i-1))   &
          &               +( v(i+1,j+1)-v(i,j+1)   )/(x(i+1)-x(i))  )&
          &         -0.25*(( u(i,j)    -u(i,j-1)   )/(y(j)-y(j-1))   &
          &               +( u(i+1,j)  -u(i+1,j-1) )/(y(j)-y(j-1))   &
          &               +( u(i,j+1)  -u(i,j)     )/(y(j+1)-y(j))   &
          &               +( u(i+1,j+1)-u(i+1,j)   )/(y(j+1)-y(j))  )
        end do
      end do 

      kk=0
      do j=s,e
        do i=1,nx-1
          kk=kk+1
          roPak(kk)=ro(i,j)
          PPak(kk) =p(i,j)
          UPak(kk) =u(i,j)
          VPak(kk) =v(i,j)
          PhiPak(kk) =Phi(i,j)
          Ib_1Pak(kk)=Ib_1(i,j)
          VortPak(kk)=Vort(i,j)
        end do
      end do

      if (Parallel) then 

        call MPI_GATHER (Y(s)       ,e-s+1         ,MPI_dOUBLE_PRECISION ,YTo(1)   ,e-s+1         ,MPI_dOUBLE_PRECISION,0,comm1d,ierr )  !! I don't know if the receive bufer should be one or not
        call MPI_GATHER (roPak(1)   ,(nx-1)*(e-s+1),MPI_dOUBLE_PRECISION ,roTo(1)  ,(nx-1)*(e-s+1),MPI_dOUBLE_PRECISION,0,comm1d,ierr )
        call MPI_GATHER (PPak(1)    ,(nx-1)*(e-s+1),MPI_dOUBLE_PRECISION ,PTo(1)   ,(nx-1)*(e-s+1),MPI_dOUBLE_PRECISION,0,comm1d,ierr )
        call MPI_GATHER (UPak(1)    ,(nx-1)*(e-s+1),MPI_dOUBLE_PRECISION ,UTo(1)   ,(nx-1)*(e-s+1),MPI_dOUBLE_PRECISION,0,comm1d,ierr )
        call MPI_GATHER (VPak(1)    ,(nx-1)*(e-s+1),MPI_dOUBLE_PRECISION ,VTo(1)   ,(nx-1)*(e-s+1),MPI_dOUBLE_PRECISION,0,comm1d,ierr )
        call MPI_GATHER (PhiPak(1)  ,(nx-1)*(e-s+1),MPI_dOUBLE_PRECISION ,PhiTo(1) ,(nx-1)*(e-s+1),MPI_dOUBLE_PRECISION,0,comm1d,ierr )
        call MPI_GATHER (Ib_1Pak(1) ,(nx-1)*(e-s+1),MPI_dOUBLE_PRECISION ,Ib_1To(1),(nx-1)*(e-s+1),MPI_dOUBLE_PRECISION,0,comm1d,ierr )
        call MPI_GATHER (VortPak(1) ,(nx-1)*(e-s+1),MPI_dOUBLE_PRECISION ,VortTo(1),(nx-1)*(e-s+1),MPI_dOUBLE_PRECISION,0,comm1d,ierr )
      end if 
      
      if (MyRank.eq.0) then
        kk=0
        write(35,*) 'zone i=',nx-1,' j=',ny-1
        call flush(35)
        Do j=1,ny-1
          Do i=1,nx-1
            kk=kk+1
            write(35,135) x(i),YTo(j),UTo(kk),VTo(kk),pTo(kk),roTo(kk),PhiTo(kk),Ib_1To(kk),VortTo(kk)
            call flush(35)
            135  format (9(1x,e15.7))
          end do
        end do
      end if

 
    end if


    if (Wave_Gen) then


      do ii=1,14

        igage=1 !! for problems with out free surface 
        do  i=1,nx-1
          if (x(i).le.xgage(ii).AND.x(i+1).ge.xgage(ii)) then 
            igage=i
            exit 
          end if 
        end do
 
        yfreeG_cpu(ii)=0
        do j=s,e
          if (0.5*( phi(igage+1,j)+phi(igage,j) ).le.0.AND.0.5*( phi(igage+1,j+1)+phi(igage,j+1) ).ge.0 ) then             !! may be more precsiness can be done !! 
            yfreeg_cpu(ii)=y(j)+ ( y(j+1)-y(j) )/ ( 0.5*( phi(igage+1,j+1)+phi(igage,j+1) )-0.5*( phi(igage+1,j)+phi(igage,j) ) ) & 
            & *( -0.5*( phi(igage+1,j)+phi(igage,j) ) ) 
            exit 
          end if 
        end do
 
      end do 

      if (Parallel) then
        call MPI_Allreduce( yfreeG_CPU ,YFreeG  ,14,MPI_Double_Precision, MPI_Sum, comm1d, ierr )
      end if
  
      if (MyRank.eq.0) then  
        write(85,1200) tp*dt,YFreeG(1),YFreeG(2),YFreeG(3),YFreeG(4),YFreeG(5),YFreeG(6), &
                &            YFreeG(7),YFreeG(8),YFreeG(9),YFreeG(10),YFreeG(11),YFreeG(12),YFreeG(13),YFreeG(14)
      end if 
      1200  format (15(1x,e15.7))


    end if 




    !if (tp.eq.tstep) then 
    !  Do j=1,ny-1
    !    write(45,145) UTo((nx-1)*j-int(nx/2) ) 
    !    call flush(45)
    !  end do

    !  Do i=1,nx-1
    !    write(55,155) VTo( (nx-1)*int((ny-1)/2)+i)  
    !    call flush(55)
    !  end do
    !end if 
    ! 145  format (1(1x,e15.5))
    ! 155  format (1(1x,e15.5))
  !end if


  !YFreeG=0
  !do j=s,e
  !  if ( (Phi(int((nx-1)/2),j)*Phi(int((nx-1)/2),j+1)).lt.0 ) then 
  !    YFreeG=y(j)
  !    exit
  !  end if
  !end do 

  !if (Parallel) then
  !  call MPI_Allreduce( YFreeG, YFreeG_T, 1, MPI_DOUBLE_PRECISION,MPI_SUM, comm1d, ierr ) 
  !end if 
 
  !165  format (2(1x,e15.5))

  !if (MyRank.eq.0) then
  !  write (65,165) tp*dt, YFreeG_T
  !end if


  return
end subroutine 








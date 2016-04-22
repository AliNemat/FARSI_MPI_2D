module M_MPI_Solid

  use mpi
  use M_MPI_General, only: U,V,ro,x,y,hx,hy,Landa,Lx,Ly,nx,Parallel,Solid_Finder_Num_R
  use M_MPI_Exch
  use M_MPI_SolidFinder_Num

  implicit none
  integer,parameter                                  ::   NPoints_1=1
  integer,parameter                                  ::   NPoints_2=8
  
  double precision,dimension (:,:),allocatable       ::   Ib_1,Ib_2
  double precision,dimension (2)                     ::   XBar,UBar
  double precision                                   ::   Omega,IniH

  double precision,dimension (NPoints_1+1,2)         ::   Points_1
  double precision,dimension (NPoints_2+1,2)         ::   Points_2
  double precision                                   ::   LenX_1,LenY_1
  double precision                                   ::   LenX_2,LenY_2
  double precision                                   ::   R_Cyl
  logical                                            ::   smooth 

  contains




  subroutine Solid_Ini(XRef,YRef,Gap,s,e,comm1d,nbrbottom,nbrtop)

    implicit none
    integer,intent(in)                                             ::   s,e

    double precision,intent(in)                                    ::   XRef,YRef,Gap
    integer,intent(in)                                             ::   comm1d,nbrbottom,nbrtop


    allocate ( Ib_1(0:nx,s-1:e+1),Ib_2(0:nx,s-1:e+1) )

    include 'Par_Solid.txt'

    if (Solid_Finder_Num_R) then 
      include 'Par_Solid_1.txt'
      include 'Par_Solid_2.txt'
    else
      include 'Par_Solid_Cyl.txt'
    end if    
    
    if (Solid_Finder_Num_R) then

 
      points_1(npointS_1+1,1)=points_1(1,1)
      points_1(npointS_1+1,2)=points_1(1,2)

      points_2(npointS_2+1,1)=points_2(1,1)
      points_2(npointS_2+1,2)=points_2(1,2)
      !call Solid_Finder_Anal(Ib_1,Points_1(1,:),r_Cyl,Gap,x,y,nx,s,e)

      !! For now npoints_1 should be equal to npoints_2 for allocating !!
      call Solid_Finder_Num_Ini(npoints_1,s,e)
   
      call Solid_Finder_Num(s,e,Lenx_1,LenY_1,npoints_1,points_1,Xbar,Ib_1)    
      call Solid_Finder_Num(s,e,Lenx_2,LenY_2,npoints_2,points_2,Xbar,Ib_2)    

      if (Parallel) then 
        call Dy_A_Exch(Ib_1  ,nx, s, e, comm1d, nbrbottom, nbrtop)
        call Dy_A_Exch(Ib_2  ,nx, s, e, comm1d, nbrbottom, nbrtop)
      end if

      if (smooth) then 

        call Solid_Smoother(Ib_1,s,e)
        call Solid_Smoother(Ib_2,s,e)

        if (Parallel) then 
          call Dy_A_Exch(Ib_1  ,nx, s, e, comm1d, nbrbottom, nbrtop)
          call Dy_A_Exch(Ib_2  ,nx, s, e, comm1d, nbrbottom, nbrtop)
        end if

      end if 


    else
      

      points_1(npointS_1+1,1)=points_1(1,1)
      points_1(npointS_1+1,2)=points_1(1,2)
      call Solid_Finder_Anal(Ib_1,XBar,R_Cyl,Gap,s,e)
      Ib_2(:,:)=0

      if (Parallel) then 
        call Dy_A_Exch(Ib_1  ,nx, s, e, comm1d, nbrbottom, nbrtop)
      end if


    end if 

    return
  end subroutine



  subroutine Solid_Update(XRef,YRef,Gap,dt,s,e,comm1d,nbrbottom,nbrtop)
    
    implicit none

    integer,intent(in)                                             ::   s,e
    double precision,intent(in)                                    ::   XRef,YRef,dt,Gap
    integer,intent(in)                                             ::   comm1d,nbrbottom,nbrtop


    call Solid_Cor_Velocity(XRef,YRef,dt,s,e,comm1d,nbrbottom,nbrtop)

    if (Parallel) then
      call UVPhi_Exch(U  , nx,s, e, comm1d, nbrbottom, nbrtop)
      call UVPhi_Exch(V  , nx,s, e, comm1d, nbrbottom, nbrtop)
    end if

    if (Solid_Finder_Num_R) then 


      call Solid_Update_Position(npoints_1,points_1,dt)
      call Solid_Update_Position(npoints_2,points_2,dt)

      call Solid_Finder_Num(s,e,Lenx_1,LenY_1,npoints_1,points_1,Xbar,Ib_1)
      call Solid_Finder_Num(s,e,Lenx_2,LenY_2,npoints_2,points_2,Xbar,Ib_2)

      if (Parallel) then 
        call Dy_A_Exch(Ib_1  ,nx, s, e, comm1d, nbrbottom, nbrtop)
        call Dy_A_Exch(Ib_2  ,nx, s, e, comm1d, nbrbottom, nbrtop)
      end if

      if (smooth) then 

        call Solid_Smoother(Ib_1,s,e)
        call Solid_Smoother(Ib_2,s,e)

        if (Parallel) then 
          call Dy_A_Exch(Ib_1  ,nx, s, e, comm1d, nbrbottom, nbrtop)
          call Dy_A_Exch(Ib_2  ,nx, s, e, comm1d, nbrbottom, nbrtop)
        end if

      end if


    else

        
      call Solid_Update_Position(npoints_1,points_1,dt)
      call Solid_Finder_Anal(Ib_1,Xbar,r_Cyl,Gap,s,e)

      if (Parallel) then 
        call Dy_A_Exch(Ib_1  ,nx, s, e, comm1d, nbrbottom, nbrtop)
      end if


    end if    

    return
  end subroutine 





  subroutine Solid_Cor_Velocity(XRef,YRef,dt,s,e,comm1d,nbrbottom,nbrtop)
    implicit none

    integer,intent(in)                                            ::    s,e
    integer,intent(in)                                            ::  comm1d,nbrbottom,nbrtop
    double precision,intent(in)                                    ::   XRef,YRef,dt



    double precision                                              ::   Tmass,Txm,Tym,sumu,sumv,sumIzb,Izzb
    double precision                                              ::   Tmass_cpu,Txm_cpu,Tym_cpu,sumu_cpu,sumv_cpu,sumIzb_cpu,Izzb_cpu
    double precision                                              ::   AAA,BBB,dm
    double precision,dimension(1:nx-1,s:e)                        ::   mom
    integer                                                       ::   i,j
    integer                                                       ::   ierr


    Tmass_cpu=0 
    Txm_cpu=0 
    Tym_cpu=0
    do i=1,nx-1 
      do j=s,e 

        dm=0.25*(x(i+1)-x(i-1))*(y(j+1)-y(j-1))*ro(i,j)*Ib_1(i,j)
        Tmass_cpu=Tmass_cpu+dm
        Txm_cpu=Txm_cpu+x(i)*dm
        Tym_cpu=Tym_cpu+y(j)*dm

       end do 
     end do

     if (Parallel) then
      
       call MPI_Allreduce( Tmass_cpu  ,Tmass ,1 ,MPI_Double_Precision, MPI_Sum, comm1d, ierr ) 
       call MPI_Allreduce( Txm_cpu    ,Txm   ,1 ,MPI_Double_Precision, MPI_Sum, comm1d, ierr ) 
       call MPI_Allreduce( Tym_cpu    ,Tym   ,1 ,MPI_Double_Precision, MPI_Sum, comm1d, ierr ) 

     end if  

     XBar(1)=XRef !Txm/Tmass  
     XBar(2)=YRef !Tym/Tmass
    ! print*, XBar(2),"Xbar(2)"
    

     Tmass_cpu=0 
     sumu_cpu=0 
     sumv_cpu=0 
     sumIzb_cpu=0     !!!sumIx=Hx , SumIy=Hy, SumIz=Hz !!
     Izzb_cpu=0   
     mom(:,:)=0 
     do i=1,nx-1 
       do j=s,e
 
         dm=0.25*(x(i+1)-x(i-1))*(y(j+1)-y(j-1))*ro(i,j)*Ib_1(i,j)
         Tmass_cpu=Tmass_cpu+dm
         Izzb_cpu=Izzb_cpu+ ( (x(i)-xbar(1))**2+ (y(j)-Xbar(2))**2 )*dm   

         AAA=0.5*(u(i,j)+u(i+1,j))*dm 
         BBB=0.5*(v(i,j)+v(i,j+1))*dm  
 
         sumu_cpu=sumu_cpu+AAA 
         sumv_cpu=sumv_cpu+BBB 
         mom(i,j)=(x(i)-XBar(1))*BBB -(y(j)-XBar(2))*AAA 
         sumIzb_cpu=sumIzb_cpu+ mom(i,j)  

       end do 
     end do
 

     if (Parallel) then
      
       call MPI_Allreduce( Tmass_cpu  ,Tmass  ,1 ,MPI_Double_Precision, MPI_Sum, comm1d, ierr ) 
       call MPI_Allreduce( Izzb_cpu   ,Izzb   ,1 ,MPI_Double_Precision, MPI_Sum, comm1d, ierr ) 
       call MPI_Allreduce( sumu_cpu   ,sumu   ,1 ,MPI_Double_Precision, MPI_Sum, comm1d, ierr ) 
       call MPI_Allreduce( sumv_cpu   ,sumv   ,1 ,MPI_Double_Precision, MPI_Sum, comm1d, ierr ) 
       call MPI_Allreduce( sumIzb_cpu ,sumIzb ,1 ,MPI_Double_Precision, MPI_Sum, comm1d, ierr ) 

     end if  

     !! calculating linear and angular velocity of the solid !! 
     Omega  =0 ! (sumIzb)/(Izzb)
     UBar(1) =0 !(sumu)/Tmass +dt/Tmass*(-197.58*(xbar(1)-XRef)-19.8*(UBar(1)))
     UBar(2) =0 !(sumv)/Tmass

     do i=1,nx-1 
       do j=s,e 

         u(i,j)=u(i,j)+0.5d0*( Ib_1(i-1,j)+Ib_1(i,j) )*( (  UBar(1)  -Omega*(y(j)-Xbar(2)) )-u(i,j) ) 
         v(i,j)=v(i,j)+0.5d0*( Ib_1(i,j-1)+Ib_1(i,j) )*( (  UBar(2)  +Omega*(x(i)-XBar(1)) )-v(i,j) )
 
       end do 
     end do 

     
     return
   end subroutine


  subroutine Solid_Finder_Anal(Ib_G,X_G,r_G,gap,s,e)   
 
    implicit none
    integer,intent(in)                                             ::   s,e

    double precision,dimension(2),intent(in)                       ::   X_G
    double precision,intent(in)                                    ::   r_G,Gap

    double precision,dimension(0:nx,s-1:e+1),intent(out)           ::   Ib_G

    double precision                                               ::   d1
    integer                                                        ::   i,j

    Ib_G(0:nx,s-1:e+1)=0
    do i=1,nx-1
      do j=s,e
        d1=sqrt( (x(i)-x_G(1))**2 + (y(j)-X_G(2))**2)
        Ib_G(i,j)=0.5d0+0.5d0*( (r_G-d1)**3 + 1.5d0* gap*gap *(r_G-d1) )  / ( (r_G-d1)*(r_G-d1)+gap*gap )**1.5d0 
 
        if (Ib_G(i,j).le.0.001) then 
          Ib_G(i,j)=0.d0
        else if (Ib_G(i,j).ge.0.999) then 
          Ib_G(i,j)=1.0d0 
        end if

      end do
    end do   
    
    print*, "SolidFinder",X_g(1),X_g(2)
    return 
  end subroutine

  subroutine Solid_Update_Position(npoints,point,dt) 

    implicit none

    integer         ,intent(in)                               :: npoints
    double precision,intent(inout), dimension(npoints+1,2)    :: point
    double precision,intent(in)                               :: dt

    double precision,dimension(npoints+1,2)    :: pointTmp
    double precision,dimension(2)              :: RTmp
    integer                                    :: i

    pointTmp(:,:)=point(:,:)

    do i=1,NPoints+1
      RTmp(:)=point(i,:)-XBar(:)
      point(i,1)=pointTmp(i,1)+dt*(Ubar(1)-Omega*(pointTmp(i,2)-xBar(2)))
      point(i,2)=pointTmp(i,2)+dt*(Ubar(2)+Omega*(pointTmp(i,1)-XBar(1)))
    end do 

    return

  end subroutine
 
  Subroutine Solid_Smoother(Ib_x,s,e)
    implicit none
    integer,intent(in)                                             ::   s,e
    double precision,dimension(0:nx,s-1:e+1),intent(inout)         ::   Ib_x

    double precision,dimension(0:nx,s-1:e+1)                       ::   Ib_x_Temp
    integer                                                        ::   i,j

    Ib_x_Temp(:,:)=Ib_x(:,:)  !! for boundary points !!


    do i=1,nx-1 
      do j=s,e 
        Ib_x(i,j) = 0.25*(Ib_x_Temp(i,j))+ &
        &   0.125*(Ib_x_Temp(i+1,j)  +Ib_x_Temp(i-1,j)  +Ib_x_Temp(i,j-1)  +Ib_x_Temp(i,j+1)) + &
        &  0.0625*(Ib_x_Temp(i+1,j+1)+Ib_x_Temp(i+1,j-1)+Ib_x_Temp(i-1,j+1)+Ib_x_Temp(i-1,j-1))
      end do 
    end do 
  
    return 
  end subroutine   

  
end module



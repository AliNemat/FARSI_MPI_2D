module M_MPI_Mesh


  contains
 
  subroutine MeshGenerator(Chv,s,e,MyRank,comm1d,nbrbottom,nbrtop)
  use mpi
  use M_MPI_General,   only : x,y,hx,hy,lx,ly,Landa,nx,ny,Parallel
  use M_MPI_Exch
    implicit none
    double precision,dimension(2),intent(in)                       ::    Chv
    integer,intent(in)                                             ::    s,e
    integer,intent(in)                                             ::    MyRank,comm1d,nbrtop,nbrbottom


    integer                                                        :: meshgen
    dOUBLE PRECISION                                               :: XE(0:nx-1)
    double precision,allocatable,dimension(:)                      :: YE
    dOUBLE PRECISION                                               :: Lx1,Lx2,Lx3,Lx4,Ly1,Ly2,Ly3
    dOUBLE PRECISION                                               :: hxx1,hxx2,hxx3,hxx4,hyy1,hyy2,hyy3
    dOUBLE PRECISION                                               :: Amesh12,Amesh14,Amesh21,Amesh23,betam12,betam14,betam21,betam23
    dOUBLE PRECISION                                               :: Dmesh12,Dmesh14,Dmesh21,Dmesh23
    dOUBLE PRECISION                                               :: SPoint_1,SPoint_2
    dOUBLE PRECISION                                               :: SPoint_1_cpu,SPoint_2_cpu
    integer                                                        :: i,j,nx1,nx2,nx3,ny1,ny2
    integer                                                        :: ierr

    allocate(YE(s-1:e))
    SPoint_1_cpu=0
    SPoint_2_cpu=0

    include 'Par_Mesh_2D.txt'
!if (meshgen.eq.1) then 
!    x(1)=0.5*Lx/(nx-1)
!    x(-1)=-3*x(1) ; x(0)=-x(1)
!    do i=1,nx+1
!      x(i)=x(i-1)+Lx/(nx-1)
!    end do 
!    hx(:)=Lx/(nx-1)


!    y(1)=0.5*Ly/(ny-1)
!    y(-1)=-3*y(1) ; y(0)=-y(1)
!    do j=1,ny+1
!      y(j)=y(j-1)+Ly/(ny-1)
!    end do 
!    hy(:)=Ly/(ny-1)



!else

    print*,"domain size",lx,ly

    nx2=nx2+nx1
    nx3=nx3+nx2  
    ny2=ny2+ny1


    
    hxx1=1.0/dble(nx1-1) ; hxx2=1.0/dble(nx2-nx1) ;  hxx3=1.0/dble(nx3-nx2) ; hxx4=1.0/dble(nx-nx3)
    hyy1=1.0/dble(ny1-1) ; hyy2=1.0/dble(ny2-ny1) ;  hyy3=1.0/dble(ny-ny2)  
    

    Amesh12=1/(2*betam12)*log(  (  1+( exp(betam12)-1 )*Dmesh12/Lx2 )/(  1+( exp(-betam12)-1 )*Dmesh12 /Lx2   )    ) 
    Amesh14=1/(2*betam14)*log(  (  1+( exp(betam14)-1 )*Dmesh14/Lx4 )/(  1+( exp(-betam14)-1 )*Dmesh14 /Lx4   )    ) 


    Amesh21=1/(2*betam21)*log(  (  1+( exp(betam21)-1 )*Dmesh21/Ly1 )/(  1+( exp(-betam21)-1 )*Dmesh21 /Ly1   )    ) 
    Amesh23=1/(2*betam23)*log(  (  1+( exp(betam23)-1 )*Dmesh23/Ly3 )/(  1+( exp(-betam23)-1 )*Dmesh23 /Ly3   )    )

  
     
    do i=0,nx1-1   !! like the whole domain disceretization, nx1-1 number of grid points are used for this section. nx1 is used for the next section.
      XE(i)=hxx1*Lx1*dble(i) 
    end do
    
    do i=nx1,nx2-1  !nx-1
      XE(i)=XE(nx1-1)+Dmesh12* (  1+ (  sinh ( betam12*(dble(i-(nx1-1))*hxx2-Amesh12) )  )/sinh(betam12*Amesh12)  )
    end do 
     
    do i=nx2,nx3-1
      XE(i)=XE(nx2-1)+hxx3*Lx3*dble( i-(nx2-1) )
    end do 

    do i=nx3,nx-1
      XE(i)=XE(nx3-1)+Dmesh14* (  1+ (  sinh ( betam14*(dble(i-(nx3-1))*hxx4-Amesh14) )  )/sinh(betam14*Amesh14)  )
    end do
     

    if (s.eq.1) then 
      YE(0)=0
    endif 

     if ( (s.ge.1.AND.s.le.ny1).OR.(e.ge.1.AND.e.le.ny1 ).OR.(s.le.1.AND.e.ge.ny1) )then
       do j=max(s,1),min(ny1-1,e)  
         YE(j)=Dmesh21* (  1+ (  sinh ( betam21*(dble(j)*hyy1-Amesh21) )  )/sinh(betam21*Amesh21)  )
       end do
     end if
 
     if ((ny1-1).ge.s.AND.(ny1-1).le.e) then
       SPoint_1_CPU=YE(ny1-1)
     end if 

     if (parallel) then 
       call MPI_Allreduce( SPoint_1_CPU ,SPoint_1  ,1 ,MPI_Double_Precision, MPI_Sum, comm1d, ierr )      
     endif 
     
          
     if ( (s.ge.ny1.AND.s.le.(ny2-1)).OR.(e.ge.ny1.AND.e.le.(ny2-1)).OR.(s.le.ny1.AND.e.ge.(ny2-1)) )then
       do j=max(ny1,s),min(ny2-1,e)  
         YE(j)=SPoint_1+hyy2*Ly2*dble(j-(ny1-1))
       end do 
     end if 

     if ((ny2-1).ge.s.AND.(ny2-1).le.e) then
       SPoint_2_CPU=YE(ny2-1)
     end if
        
     if (parallel) then 
       call MPI_Allreduce( SPoint_2_CPU ,SPoint_2  ,1 ,MPI_Double_Precision, MPI_Sum, comm1d, ierr )      
     end if 



     if ( (s.ge.ny2.AND.s.le.(ny-1)).OR.(e.ge.ny2.AND.e.le.(ny-1)).OR.(s.le.ny2.AND.e.ge.(ny-1)) )then
       do j=max(ny2,s),min(ny-1,e) 
         YE(j)=SPoint_2+Dmesh23* (  1+ (  sinh ( betam23*(dble(j-(ny2-1))*hyy3-Amesh23) )  )/sinh(betam23*Amesh23)  )
       end do 
     end if 
     
     if (Parallel) then 
       call YE_Exch(YE, s, e, comm1d, nbrbottom, nbrtop)
     end if
     !!! it should be searched how important is the cell size on the boundary!!

     
     do i=1,nx-1
       x(i)=XE(i-1)+0.5*(XE(i)-XE(i-1))
       hx(i)=XE(i)-XE(i-1)
     end do   
     x(-1)=-3*x(1) 
     x(0)=-x(1) 

     hx(0)=hx(1)
     hx(-1)=hx(0)    

     x(nx)=x(nx-1)+ ( x(nx-1)-x(nx-2) ) 
     x(nx+1)=x(nx)+ ( x(nx)-x(nx-1) )
     
     hx(nx)=hx(nx-1)
     hx(nx+1)=hx(nx) 

     
     do j=s,e
       y(j)=YE(j-1)+0.5*(YE(j)-YE(j-1))
       hy(j)=YE(j)-YE(j-1)
       !print*, j,y(j)
     end do
     
     if (s.eq.1) then 
       y(-1)=-3*y(1) 
       y(0)=-y(1)

       hy(0)=hy(1)
       hy(-1)=hy(0)    
     end if 
       
     
     if (e.eq.(ny-1)) then 
       y(ny)=y(ny-1)+ ( y(ny-1)-y(ny-2) ) 
       y(ny+1)=y(ny)+ ( y(ny)-y(ny-1) )

       hy(ny)=hy(ny-1)
       hy(ny+1)=hy(ny) 
     end if 

     
     if (Parallel) then 
       call Yhy_Exch(Y,  s, e, comm1d, nbrbottom, nbrtop)
       call Yhy_Exch(hy, s, e, comm1d, nbrbottom, nbrtop)
     endif 

! end if 



!OPEN(205,file='meshx.plt')
!OPEN(215,file='meshy.plt')
!OPEN(235,file='grideval2DXY.plt')

!do i=1,nx-1
!write(205,*) i,x(i),hx(i)
!end do 
!call flush (205)
!do j=1,ny-1
!write(215,*) j,y(j),hy(j)
!end do 
!call flush (215)


!write(235,*) 'zone i=',nx-1,' k=',ny-1
!do j=1,ny-1 ;  Do i=1,nx-1

!write(235,350) x(i),y(j)
!end do ; end do
!350 format (2(1x,e15.7))
!call flush (235)

    return 
  end subroutine 
  
  


     
end module 



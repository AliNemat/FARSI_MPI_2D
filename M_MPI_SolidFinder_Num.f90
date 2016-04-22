Module M_MPI_SolidFinder_Num
  use M_MPI_General,  only: nx,x,y,hx,hy

  implicit none

  integer,dimension(:,:),allocatable                 :: ccL       
  double precision,dimension(:,:),allocatable        :: cc,CI,CF,CFS
  double precision,dimension(:),allocatable          :: mm,mm1,mm2,alpha



  contains


  subroutine  Solid_Finder_Num_Ini(npoints,s,e)
    implicit none
    integer,intent(in)                                            :: npoints,s,e


    allocate (ccL(1:nx-1,s:e))
    allocate (cc(1:nx-1,s:e),CI(1:nx-1,s:e),CF(1:nx-1,s:e),CFS(1:nx-1,s:e))
    allocate (mm(npoints),mm1(npoints),mm2(npoints),alpha(npoints))

    return 
  end subroutine


  subroutine Solid_Finder_Num(s,e,Lenx,LenY,npoints,point,XBarV,CFSOut)
    implicit none 

    integer,intent(in)                                            :: s,e,npoints
    double precision,dimension(npoints+1,2),intent(in)            :: point
    double precision,dimension(2),intent(in)                      :: XBarV
    double precision,intent(in)                                   :: LenX,LenY

    double precision,intent(out),dimension(0:nx,s-1:e+1)          :: CFSOut
    
    integer                                                       :: i,j

    call BorderFinder  (s,e,Lenx,LenY,npoints,point,XBarV)
    call insideFinder  (s,e,Lenx,lenY,npoints,point,XBarV)
    call FractionFinder(s,e          ,npoints,point      )
    CFS(:,:)=CF(:,:)
    CFSOut(:,:)=0
!! Assume no contact with boundary   
    do i=1,nx-1
      do j=s,e
        CFSOUT(i,j)=CFS(i,j)
      end do 
    end do 
 
    return

  end subroutine



     

  subroutine BorderFinder(s,e,Lenx,LenY,npoints,point,XBarV)
    implicit none

    integer,intent(in)                                            :: s,e,npoints
    
    double precision,dimension(npoints+1,2),intent(in)            :: point
    double precision,dimension(2),intent(in)                      :: XBarV
    double precision,intent(in)                                   :: LenX,LenY
    
    integer kk,i,j
 
    do kk=1,npointS

      if (point(kk+1,1).ne.point(kk,1)) then 

        mm(kk)=(point(kk+1,2)-point(kk,2))/(point(kk+1,1)-point(kk,1))

        mm1(kk)=-mm(kk)/sqrt(1+mm(kk)*mm(kk))
        mm2(kk)=  1.0/sqrt(1+mm(kk)*mm(kk))
       
        alpha(kk)=(point(kk,2)-mm(kk)*point(kk,1))/sqrt(1+mm(kk)*mm(kk))

      else

        mm1(kk)=1
        mm2(kk)=0
        alpha(kk)=point(kk,1)

      end if

    end do 
    !print*, "mm1",mm1(:)
    !print*, "mm2",mm2(:)
    !print*, "alpha",alpha(:)

    ccl(:,:)=0
    do i=1, nx-1
      if ( abs (x(i)-XBarV(1)).gt.(LenX) )  then
        cc(i,:)=0
      else
     
        do j=s,e
          if ( abs (y(j)-XBarV(2)).gt.(lenY) ) then
            cc(i,j)=0
          else 
            do kk=1,npoints

              if (mm2(kk).eq.0) then
    !         print*, " I am here1 " 
                if ( (x(i)-0.5*hx(i)).lt.alpha(kk).AND.(x(i)+0.5*hx(i)).gt.Alpha(kk) ) then
                  if ( ( (x(i)-point(kk,1  ))**2 + (y(j)-point(kk,2  ))**2 )  .lt.( (point(kk,1)-point(kk+1,1))**2 + (point(kk,2)-point(kk+1,2))**2 ).AND. &
                     (   (x(i)-point(kk+1,1))**2 + (y(j)-point(kk+1,2))**2 )  .lt.( (point(kk,1)-point(kk+1,1))**2 + (point(kk,2)-point(kk+1,2))**2 ) ) then

                   ! print*, " I am here2 " 
                    cc(i,j)=1.0
                    ccl(i,j)=kk
                    exit
                  end if  
                end if 
              else
                if (      max(  (Alpha(kk)-mm1(kk)*(x(i)-0.5*hx(i)))/mm2(kk),(Alpha(kk)-mm1(kk)*(x(i)+0.5*hx(i)))/mm2(kk) ).lt. (y(j)-0.5*hy(j)) ) then
                  cc(i,j)=0
                else if ( min ( (Alpha(kk)-mm1(kk)*(x(i)-0.5*hx(i)))/mm2(kk),(Alpha(kk)-mm1(kk)*(x(i)+0.5*hx(i)))/mm2(kk) ).gt. (y(j)+0.5*hy(j)) ) then 
                  cc(i,j)=0
                else
                 
                  if ( ( (x(i)-point(kk,1  ))**2 + (y(j)-point(kk,2  ))**2 )  .lt.( (point(kk,1)-point(kk+1,1))**2 + (point(kk,2)-point(kk+1,2))**2 ).AND. &
                     (   (x(i)-point(kk+1,1))**2 + (y(j)-point(kk+1,2))**2 )  .lt.( (point(kk,1)-point(kk+1,1))**2 + (point(kk,2)-point(kk+1,2))**2 ) ) then

                       !print*, " I am here2 " 
                    cc(i,j)=1.0
                    ccl(i,j)=kk
                    exit
                  end if  
                end if 
              end if
 
            end do
       
          end if 
        end do

      end if
    end do

!OPEN(2350,file='BorderFinder.plt')

!write(2350,*) 'zone i=',nx-1,' k=',ny-1
!do j=1,ny-1 ;  Do i=1,nx-1

!write(2350,3500) x(i),y(j),cc(i,j)
!end do ; end do
!3500  format (3(1x,e15.7))
!call flush (2350)


    return

  end subroutine


  subroutine InsideFinder(s,e,Lenx,lenY,npoints,point,XBarV)
    implicit none

    integer,intent(in)                                           :: s,e,npoints
    double precision,dimension(2),intent(in)                     :: XBarV
    double precision,dimension(npoints+1,2),intent(in)           :: point
    double precision,intent(in)                                  :: LenX,LenY
  

    double precision                                             :: YCross
    integer                                                      :: i,j,kk

    CI(:,:)=0

    do i=1, nx-1
      if ( abs (x(i)-XBarV(1)).gt.(LenX) )  then
        cI(i,:)=0
      else
     
        do j=s,e
          if ( abs (y(j)-XBarV(2)).gt.(lenY) ) then
            cI(i,j)=0
          else
            do kk=1,npoints
            
              if (mm2(kk).ne.0) then
                YCross=(Alpha(kk)-mm1(kk)*x(i))/mm2(kk)
                if (((YCross-point(kk  ,2))**2 + (X(i)-point(kk  ,1))**2).le.( (point(kk,1)-point(kk+1,1))**2 + (point(kk,2)-point(kk+1,2))**2 ).AND. &
                &   ((YCross-point(kk+1,2))**2 + (X(i)-point(kk+1,1))**2).le.( (point(kk,1)-point(kk+1,1))**2 + (point(kk,2)-point(kk+1,2))**2 ).AND. &
                &  ((Alpha(kk)-mm1(kk)*x(i))/mm2(kk)).gt.y(j)                                                                                      ) then
                    
                  cI(i,j)=cI(i,j)+1
                end if  
              end if
          
            end do
       
            if (cI(i,j).eq.1.AND.cc(i,j).eq.0) then
              CI(i,j)=1.0
            else 
              CI(i,j)=0
            end if 
          end if 
        end do

      end if
    end do
 
         
!OPEN(2351,file='InsideFinder.plt')

!write(2351,*) 'zone i=',nx-1,' k=',ny-1
!do j=1,ny-1 ;  Do i=1,nx-1

!write(2351,3501) x(i),y(j),cI(i,j),cc(i,j)
!end do ; end do
!3501  format (4(1x,e15.7))
!call flush (2351)


    return
  end subroutine 


  subroutine FractionFinder(s,e,npoints,point)
    implicit none

    integer,intent(in)                                           :: s,e,npoints
    double precision,dimension(npoints+1,2),intent(in)           :: point


    double precision,dimension(npoints)                          :: mm1N,mm2N,alphaN  
    integer                                                      :: kk,i,j


    CF(:,:)=0 
    do i=1,nx-1
      do j=s,e
        if (cc(i,j).ne.0) then
          kk=ccl(i,j) 
          if       ( (point(kk+1,1)-point(kk,1)).gt.0.AND.( point(kk+1,2)-point(kk,2) ).gt.0 ) then
            mm1N(kk)= mm1(kk)
            mm2N(kk)=-mm2(kk)
            AlphaN(kk)=Alpha(kk) -mm1(kk)*(x(i)-0.5*hx(i)) -mm2(kk)*(y(j)+0.5*hy(j))
            !print*, "hx,hy",hx(i),hy(j)
            !print*,"old value",mm1(1),mm2(1), Alpha(1)
            !print*,"new value",mm1N(1),mm2N(1), AlphaN(1)
            !print*, "x(i),y(j)", x(i),y(j) 
 
            !print*,  "I am here 1"
            if (AlphaN(kk).le.0.0) then
              mm1N(kk)=-mm1N(kk)
              mm2N(kk)=-mm2N(kk)
              AlphaN(kk)=-AlphaN(kk)
            end if 
          else if  ( (point(kk+1,1)-point(kk,1)).le.0.AND.( point(kk+1,2)-point(kk,2) ).ge.0 ) then
            mm1N(kk)= mm1(kk)
            mm2N(kk)= mm2(kk)
            AlphaN(kk)=Alpha(kk) -mm1(kk)*(x(i)-0.5*hx(i)) -mm2(kk)*(y(j)-0.5*hy(j))
          ! print *, " I am here 2"
            if (AlphaN(kk).le.0.0) then
              mm1N(kk)=-mm1N(kk)
              mm2N(kk)=-mm2N(kk)
              AlphaN(kk)=-AlphaN(kk)
            end if 
          else if ( (point(kk+1,1)-point(kk,1)).le.0.AND.( point(kk+1,2)-point(kk,2) ).le.0 ) then
          ! print*, "I am here 3 "
            mm1N(kk)=-mm1(kk)
            mm2N(kk)= mm2(kk)
            AlphaN(kk)=Alpha(kk) -mm1(kk)*(x(i)+0.5*hx(i)) -mm2(kk)*(y(j)-0.5*hy(j))
            if (AlphaN(kk).le.0.0) then
              mm1N(kk)=-mm1N(kk)
              mm2N(kk)=-mm2N(kk)
              AlphaN(kk)=-AlphaN(kk)
            end if 
          else
            !print*,  "I am here 4"
            mm1N(kk)=-mm1(kk)
            mm2N(kk)=-mm2(kk)
            AlphaN(kk)=Alpha(kk) -mm1(kk)*(x(i)+0.5*hx(i)) -mm2(kk)*(y(j)+0.5*hy(j))
            if (AlphaN(kk).le.0.0) then
              mm1N(kk)=-mm1N(kk)
              mm2N(kk)=-mm2N(kk)
              AlphaN(kk)=-AlphaN(kk)
            end if 
          end if
        
          if     (mm2N(kk).eq.0) then
            CF(i,j)=(AlphaN(kk)/mm1N(kk))/hx(i)
            !print*," I am here 5 "
          else if(mm1N(kk).eq.0) then
            CF(i,j)=(AlphaN(kk)/mm2N(kk))/hy(j)
            !print *,  "I am here 6" 
          else
            ! print*, "I am here 7"
            CF(i,j)=(AlphaN(kk)**2)/(2*mm1N(kk)*mm2N(kk))/(hx(i)*hy(j))*(1 &
                                                           & -HFF( AlphaN(kk)-mm1N(kk)*hx(i) )*( (AlphaN(kk)-mm1N(kk)*hx(i))/AlphaN(kk) )**2 &
                                                           & -HFF( AlphaN(kk)-mm2N(kk)*hy(j) )*( (AlphaN(kk)-mm2N(kk)*hy(j))/AlphaN(kk) )**2  )                                                        
            !print*,"CF",CF(i,j) 
          end if
      
        end if
        if (cI(i,j).eq.1.0) then
          CF(i,j)=1.0

             !print*, "I am here inside 2 "
        end if 
      end do
    end do 
  

  !if (Test) then   
!write(2352,*) 'zone i=',nx-1,' k=',ny-1
!  do j=1,ny-1 ;  Do i=1,nx-1

!    write(2352,3502) x(i),y(j),cI(i,j),cc(i,j),CF(i,j),dble (ccl(i,j))
!  end do ; end do
!  3502  format (6(1x,e15.7))
!  call flush (2352)

!end if 

    return

  end subroutine 





  DOUBLE PRECISION FUNCTION HFF(a)

    double precision :: a 
    if (a.gt.0.0) then 
      HFF=1.0
    else 
      HFF=0 
    end if 

  end function HFF
          
      

  end  module 



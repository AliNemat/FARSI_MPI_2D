Module M_MPI_Solver
  implicit none 
  private                                                  ::    DIFF
  contains

subroutine Poisson_Solver_Parallel(Apx,Amx,Apy,Amy,Ap,Q,pdif,beta,P,s,e,comm1d,nbrbottom,nbrtop)

!Parallel
    use mpi
!Parallel
    use M_MPI_General, only: nx,ny
    implicit none  
    integer,intent(in)                                                 ::    s,e
    integer,intent(in)                                                 ::    comm1d,nbrbottom,nbrtop
    double precision,intent(in)                                        ::    pdif,beta
    double precision,dimension (1:nx-1,s:e)  ,intent(in)               ::    Apx,Amx,Apy,Amy,AP,Q
    double precision,dimension (0:nx,s-1:e+1),intent(inout)            ::    p

    double precision                                                   ::    diffnorm,dwork
    double precision,dimension (0:nx,s-1:e+1)                          ::    p1,p2
    integer i,j,it
    integer                                                            ::    status(MPI_STATUS_SIZE), ierr


    p1(0:nx,s-1:e+1)=p(0:nx,s-1:e+1)
    p2(0:nx,s-1:e+1)=p(0:nx,s-1:e+1)
    do it=1,10000
       call Sweep1D   (Apx,Amx,Apy,Amy,Ap,Q,beta,p1,p2,nx,s,e)  
       call PTemp_Exch( p2, nx, s, e, comm1d, nbrbottom, nbrtop )     
       call Sweep1D   (Apx,Amx,Apy,Amy,Ap,Q,beta,p2,p1,nx,s,e)  
       call PTemp_Exch( p1, nx, s, e, comm1d, nbrbottom, nbrtop )   
  
       dwork = DIFF   ( p1, p2, nx, s, e )
       call MPI_Allreduce( dwork, diffnorm, 1, MPI_DOUBLE_PRECISION,MPI_SUM, comm1d, ierr ) 
       !diffnorm=dwork

       if (diffnorm .lt. pdif) then
            exit
       end if
     end do 

     if (s.eq.1) then
      print *, 'Converged after ', 2*it, 'Iterations','diff',diffnorm
     end if

     p(0:nx,s-1:e+1)=p1(0:nx,s-1:e+1)


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
   end subroutine 




  Subroutine Sweep1D (Apx,Amx,Apy,Amy,Ap,Q,beta,pold,pnew,nx1,s,e)  
    implicit none

    integer,intent(in)                                                 ::    nx1,s,e
    double precision,intent(in)                                        ::    beta
    double precision,dimension (1:nx1-1,s:e)  ,intent(in)               ::    Apx,Amx,Apy,Amy,AP,Q
    double precision,dimension (0:nx1,s-1:e+1),intent(in)            ::       pold
    double precision,dimension (0:nx1,s-1:e+1),intent(out)            ::      pnew


    integer i,j

    do i=1,nx1-1
      do j=s,e 
                              
        pnew(i,j)=beta*( ( Apx(i,j)*Pold(I+1,j)+Amx(i,j)*Pold(I-1,j)+Apy(i,j)*Pold(I,j+1)+Amy(i,j)*Pold(I,j-1)-Q(I,j)  )/(-AP(i,j))  ) &
        & +(1-beta)*pold(i,j)
 
      end do 
    end do



    return

  end subroutine


subroutine PTemp_Exch( pex, nx1, s, e, comm1d, nbrbottom, nbrtop )
!Parallel
  use mpi
!Parallel
  implicit none
  integer,intent(in)                                            :: nx1, s, e
  double precision,dimension(0:nx1,s-1:e+1),intent(inout)        :: Pex
  integer,intent(in)                                            :: comm1d, nbrbottom, nbrtop
  integer                                                       :: status(MPI_STATUS_SIZE), ierr


        call MPI_SENDRECV(  &
     &            pex(0,e)  , nx1+1, MPI_DOUBLE_PRECISION, nbrtop   , 0, &
     &            pex(0,s-1), nx1+1, MPI_DOUBLE_PRECISION, nbrbottom, 0, &
     &            comm1d, status, ierr )
        call MPI_SENDRECV( &
     &            pex(0,s),   nx1+1, MPI_DOUBLE_PRECISION, nbrbottom, 1, &
     &            pex(0,e+1), nx1+1, MPI_DOUBLE_PRECISION, nbrtop,    1, &
     &            comm1d, status, ierr )
  return
end subroutine 





  double precision function DIFF( a, b, nx1, s, e )
      integer nx1, s, e
      double precision a(0:nx1, s-1:e+1), b(0:nx1, s-1:e+1)
      double precision sump
      integer i, j

      sump = 0.0d0
      do  j=s,e
         do  i=1,nx1-1
            sump = sump + (a(i,j) - b(i,j)) ** 2
         end do
      end do

      diff = sump
      
      !print*,diff,"diff"
    return

  end function


end module 


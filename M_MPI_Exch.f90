Module M_MPI_Exch

implicit none

contains




subroutine Dy_A_Exch(Dy_A,nx1, s, e, comm1d, nbrbottom, nbrtop)

!Parallel
  use mpi
!Parallel
  implicit none
  integer,intent(in)                                         ::  nx1,s,e
  double precision,dimension (0:nx1,s-1:e+1),intent(inout)    ::  Dy_A
  integer,intent(in)                                         ::  comm1d, nbrbottom, nbrtop
  integer                                                    ::  status(MPI_STATUS_SIZE), ierr



     call MPI_SENDRECV( &
     &            Dy_A(0,e)  ,nx1+1, MPI_DOUBLE_PRECISION, nbrtop,    0, &
     &            Dy_A(0,s-1),nx1+1, MPI_DOUBLE_PRECISION, nbrbottom, 0, &
     &            comm1d, status, ierr )

     call MPI_SENDRECV( &
     &            Dy_A(0,s)  ,nx1+1, MPI_DOUBLE_PRECISION, nbrbottom, 1, &
     &            Dy_A(0,e+1),nx1+1, MPI_DOUBLE_PRECISION,nbrtop    , 1, &
     &            comm1d, status, ierr )


  return
end subroutine 



subroutine UVPhi_Exch(UV, nx1, s, e, comm1d, nbrbottom, nbrtop)
!Parallel
  use mpi
!Parallel
  implicit none
  integer,intent(in)                                             :: s, e,nx1
  double precision,dimension (-1:nx1+1,s-2:e+2),intent(inout)     :: UV
  integer,intent(in)                                             :: comm1d, nbrbottom, nbrtop
  integer                                                        :: status(MPI_STATUS_SIZE), ierr


     call MPI_SENDRECV( &
     &            UV(-1,e)  ,nx1+3, MPI_DOUBLE_PRECISION, nbrtop,    0, &
     &            UV(-1,s-1),nx1+3, MPI_DOUBLE_PRECISION, nbrbottom, 0, &
     &            comm1d, status, ierr )

     call MPI_SENDRECV( &
     &            UV(-1,s)  ,nx1+3, MPI_DOUBLE_PRECISION, nbrbottom, 1, &
     &            UV(-1,e+1),nx1+3, MPI_DOUBLE_PRECISION, nbrtop   , 1, &
     &            comm1d, status, ierr )


!    print*, "I am here UV "
  return
end subroutine

subroutine YE_Exch(YE, s, e, comm1d, nbrbottom, nbrtop)
!Parallel
  use mpi
!Parallel
  implicit none
  integer,intent(in)                                             :: s, e
  double precision,dimension (s-1:e),intent(inout)               :: YE
  integer,intent(in)                                             :: comm1d, nbrbottom, nbrtop
  integer                                                        :: status(MPI_STATUS_SIZE), ierr


     call MPI_SENDRECV( &
     &            YE(e)  ,1, MPI_DOUBLE_PRECISION, nbrtop,    0, &
     &            YE(s-1),1, MPI_DOUBLE_PRECISION, nbrbottom, 0, &
     &            comm1d, status, ierr )



  return
end subroutine

subroutine Yhy_Exch(Yhy, s, e, comm1d, nbrbottom, nbrtop)
!Parallel
  use mpi
!Parallel
  implicit none
  integer,intent(in)                                             :: s, e
  double precision,dimension (s-2:e+2),intent(inout)             :: Yhy
  integer,intent(in)                                             :: comm1d, nbrbottom, nbrtop
  integer                                                        :: status(MPI_STATUS_SIZE), ierr


     call MPI_SENDRECV( &
     &            Yhy(e)  ,1, MPI_DOUBLE_PRECISION, nbrtop,    0, &
     &            Yhy(s-1),1, MPI_DOUBLE_PRECISION, nbrbottom, 0, &
     &            comm1d, status, ierr )

     call MPI_SENDRECV( &
     &            Yhy(s)  ,1, MPI_DOUBLE_PRECISION, nbrbottom, 1, &
     &            Yhy(e+1),1, MPI_DOUBLE_PRECISION, nbrtop   , 1, &
     &            comm1d, status, ierr )


  return
end subroutine

end module 



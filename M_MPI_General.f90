module M_MPI_General

  implicit none 
  save
  logical, parameter                                  :: Parallel           =.True.
  logical, parameter                                  :: PrintD             =.True.
  logical, parameter                                  :: Uni_Flow           =.True.
  logical, parameter                                  :: Wave_Gen           =.False.
  logical, parameter                                  :: Solid_Finder_Num_R =.False.
  logical, parameter                                  :: Free_Surface       =.False.
  integer, parameter                                  :: nx=285, ny=161  !nx =400, ny=120 !nx=460, ny=186 !nx=820, ny=373 !nx=512, ny=224 !nx=390,ny=224  !nx=355, ny=204 
   
  double precision, parameter                         :: pi=3.141592653, gy=0 , froude=1.0, landa=2.57 
  double precision, parameter                         :: Lx=10.0, Ly=2.5 !!lx=22.0*Landa,  Ly=2.0*Landa

  double precision,dimension(:,:), allocatable        ::  u
  double precision,dimension(:,:), allocatable        ::  v
  double precision,dimension(:,:), allocatable        ::  ro,miuv
  double precision,dimension(:)  , allocatable        ::  y,hy
  double precision, dimension   (-1:nx+1)             ::  x,hx
       
  contains 
  subroutine General_Ini(s,e)
    implicit none 
    integer,intent(in)                      :: s,e

    allocate (y(s-2:e+2),hy(s-2:e+2) )
    allocate ( u(-1:nx+1,s-2:e+2),v(-1:nx+1,s-2:e+2) )
    allocate ( ro(0:nx,s-1:e+1),miuv(0:nx,s-1:e+1) )

    u(-1:nx+1,s-2:e+2)=0
    v(-1:nx+1,s-2:e+2)=0

    
    return
  end subroutine 

end module




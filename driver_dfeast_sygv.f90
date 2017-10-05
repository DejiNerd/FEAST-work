!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! FEAST Driver dense example !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! solving Ax=eBx with A real and B spd --- A and B dense matrix!!!!
!!!!!!! by Eric Polizzi- 2009-2011!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program driver

  implicit none

!!!!!!!!!!!!!!!!! Feast declaration variable
  integer,dimension(64) :: fpm 
  double precision :: epsout
  double precision, parameter :: pi=3.14159265358979d0
  integer :: loop
  character(len=1) :: UPLO='F' ! 'L' or 'U' also fine

!!!!!!!!!!!!!!!!! Matrix declaration variable
  character(len=100) :: name
  integer :: n,nnz
  double precision,dimension(:,:),allocatable :: A,B

!!!!!!!!!!!!!!!!! Others
  integer :: t1,t2,tim
  integer :: i,j,k
  integer :: M0,M,info
  double precision :: Emin,Emax

  !! :,: represents [][]
  double precision,dimension(:,:),allocatable :: X ! eigenvectors
  double precision,dimension(:),allocatable :: E,res,theta ! eigenvalue+residual

  !!!!
  complex(kind=kind(1.0d0)),dimension(:),allocatable :: Zne,Wne

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!Read Coordinate format and convert to dense format
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  name='./system1'

  open(10,file=trim(name),status='old')
  read(10,*) n,nnz
  allocate(A(1:n,1:n))
  allocate(B(1:n,1:n))
  A=0.0d0
  B=0.0d0
  do k=1,nnz
     read(10,*) i,j,A(i,j),B(i,j)
  enddo
  close(10)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! INFORMATION ABOUT MATRIX !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  print *,'dense matrix -system1- size',n


  call system_clock(t1,tim)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! FEAST in dense format !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! search interval [Emin,Emax] including M eigenpairs
  Emin=0.00d0
  Emax=1.00d0 
  M0=30!400!25 !! M0>=M

!!!!!!!!!!!!! ALLOCATE VARIABLE 
  allocate(e(1:M0))     ! Eigenvalue
  allocate(X(1:n,1:M0)) ! Eigenvectors
  allocate(res(1:M0))   ! Residual (if needed)

!!!!!!!!!!!!  FEAST
  call feastinit(fpm)
  !! fpm1 just prints
  fpm(1)=1
  !! fpm2 no of contour pts
  fpm(2)=200
  !! fpm4 does loops
  fpm(4)=0
  !! fpm16 equidistant trapezoid
  fpm(16)=1

  allocate(Zne(fpm(2)),Wne(fpm(2)))
  allocate(theta(fpm(2)))
  ! call zfeast_contour(Emin, Emax,fpm(2), fpm(16), fpm(18), Zne, Wne)
  ! print *, Wne
  call dlarnv(1, (/100,200,40,7/), fpm(2), theta)
  ! print *, theta
  Zne = (cos(theta * pi) * (1.0d0,0.0d0) + sin(theta * pi) * (0.0d0, 1.0d0)) * (Emax-Emin)/2.0d0
  Zne = Zne + (Emax+Emin)/2.0d0
  ! print *, Zne
  Wne = (0.0d0,1.0d0)
  ! Wne = (1.0d0, 0.0d0)/fpm(2)
  
! !!! x for expert routine
  call dfeast_sygvx('L',N,A(1,1),N,B(1,1),N,fpm,epsout,loop,Emin,Emax,M0,E,X,M,res,info,Zne,Wne)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! POST-PROCESSING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 call system_clock(t2,tim)
  print *,'FEAST OUTPUT INFO',info
  if (info==0) then
     print *,'*************************************************'
     print *,'************** REPORT ***************************'
     print *,'*************************************************'
     print *,'SIMULATION TIME',(t2-t1)*1.0d0/tim
     print *,'# Search interval [Emin,Emax]',Emin,Emax
     print *,'# mode found/subspace',M,M0
     print *,'# iterations',loop
     print *,'TRACE',sum(E(1:M))
     print *,'Relative error on the Trace',epsout
     print *,'Eigenvalues/Residuals'
     do i=1,M
        print *,i,E(i),res(i)
     enddo
     print *,"OUTSIDE INTERVAL"	
     do i=M+1,M0
        print *,i,E(i),res(i)
     enddo
  endif


end program driver




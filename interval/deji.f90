!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! FEAST general Driver sparse !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! solving Ax=eBx or Ax=eX      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! symmetric, hermitian or general matrices !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! single or double precision !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! where A (and B if any), provided by user in coordinate format !!!!                         
!!!!!!! by Eric Polizzi- 2009-2015  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program driver_feast_sparse

  implicit none

!!!!!!!!!!!!!!!!! Feast declaration variable
  integer,dimension(64) :: fpm
  double precision :: depsout
  real :: sepsout
  integer :: loop

!!!!!!!!!!!!!!!!! Matrix declaration variable
  character(len=100) :: name
  character(len=1) :: UPLO,PRE,SHG,EG 
  integer :: n,nnza,nnzb
  real,dimension(:),allocatable :: ssa,ssb,sca,scb
  double precision,dimension(:),allocatable :: dsa,dsb,dca,dcb
  integer,dimension(:),allocatable :: isa,jsa,isb,jsb,ica,jca,icb,jcb

!!!!!!!!!!!!!!!!! Others
  integer :: t1,t2,tim
  integer :: i,j,k,pc
  integer :: M0,M,info, count, M1
  character(len=1) :: cc

  double precision :: dEmin,dEmax,dr,minE,ddiff,dEmax0
  real :: sEmin,sEmax,sr, sdiff
  double precision:: drea,dimg
  real:: srea,simg

  ! eigenvectors
  double precision,dimension(:,:),allocatable :: dX
  real,dimension(:,:),allocatable :: sX
  complex,dimension(:,:),allocatable :: cX
  complex(kind=kind(1.0d0)),dimension(:,:),allocatable :: zX

  ! eigenvalue + residual
  double precision,dimension(:),allocatable :: dres
  real,dimension(:),allocatable :: sres
  double precision,dimension(:),allocatable :: dE
  real,dimension(:),allocatable :: sE
  complex,dimension(:),allocatable :: cE
  complex(kind=kind(1.0d0)),dimension(:),allocatable :: zE



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!! read main input file !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call getarg(1,name)



  call feastinit(fpm)


!!!!!!!!!!!! DRIVER_FEAST_SPARSE input file  
  open(10,file=trim(name)//'.in',status='old')
  read(10,*) SHG ! type of eigenvalue problem "General, Hermitian, Symmetric" 
  read(10,*) EG ! type of eigenvalue probl  g== sparse generalized, e== sparse standard
  read(10,*) PRE  ! "PRE"==(s,d,c,z) resp. (single real,double real,complex,double complex) 
  read(10,*) UPLO ! UPLO==(F,L,U) reps. (Full csr, Lower csr, Upper csr) 

  if (SHG=='s') then

     if (PRE=='d') then
        read(10,*) dEmin
        read(10,*) dEmax
     elseif (PRE=='s') then
        read(10,*) sEmin
        read(10,*) sEmax
     end if
  end if


  read(10,*) M0   ! size subspace

  read(10,*) pc ! Some changes from default for fpm

  do i=1,pc
     read(10,*) j,fpm(j)
  enddo

  close(10)

  dEmin = huge(dEmin)
  dEmax = -huge(dEmax)
  M = huge(M)
  fpm(14)=2 !! stochastic approach
  ! fpm(14) = 0 !! non stochastic



!!!!!!!!!!!read matrix A
  open(10,file=trim(name)//'_A.mtx',status='old')
  k=0
  cc='%'
  do while(cc=='%')
     k=k+1 
     read(10,'(A1)') cc
  end do
  close(10)

  open(10,file=trim(name)//'_A.mtx',status='old')
  do i=1,k-1
     read(10,'(A1)') cc
  enddo
  read(10,*) n,n,nnza
  allocate(ica(nnza))
  allocate(jca(nnza))
  if (PRE=='s') then
     allocate(sca(nnza))
     do i=1,nnza
        read(10,*) ica(i),jca(i),sca(i)
     end do
  elseif (PRE=='d') then
     allocate(dca(nnza))
     do i=1,nnza
        read(10,*) ica(i),jca(i),dca(i)
     end do
  end if
  close(10)

  !! create csr format
  allocate(isa(1:n+1))
  allocate(jsa(1:nnza))


  if (PRE=='s') then
     allocate(ssa(1:nnza))
     call scoo2csr(n,nnza,ica,jca,sca,isa,jsa,ssa)
  elseif (PRE=='d') then
     allocate(dsa(1:nnza))
     call dcoo2csr(n,nnza,ica,jca,dca,isa,jsa,dsa)
     call dgershgorin(n,nnza,isa,jsa,dsa,dEmin,dEmax)
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! INFORMATION ABOUT MATRIX !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  print *,'matrix name ',trim(name)
  print *,'matrix -coordinate format- size',n
  print *,'sparse matrix A- nnz',nnza
  if (EG=='g') print *,'sparse matrix B- nnz',nnzb
  print *,''

  info=-1
  call system_clock(t1,tim)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! FEAST  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!  FEAST SYMMETRIC

if (SHG=='s') then 

   if  ((PRE=='d').and.(EG=='e')) then 
      print *,'Routine  ','dfeast_scsrev'
      allocate(dE(1:M0))     ! Eigenvalue
      allocate(dres(1:M0))   ! Residual 
      allocate(dX(1:n,1:M0))
      dEmax0= dEmax
      M1 = M0/2 ! number of eigenvalues
      ddiff = (dEmax - dEmin) * M0/n
      count = 0
      fpm(1) = 0
      i = 3
      dEmax = dEmax0/(2**i)
      do while ((M < M1).or.(M > M0))
        count = count + 1
        i = i + 1
        ddiff = (dEmax0 -dEmin)/(2**i)
        if (M > M0) then
          dEmax = dEmax - ddiff
        else if (M < M1) then
          dEmax = dEmax + ddiff
        endif
        fpm(4) = 0
        call dfeast_scsrev(UPLO,N,dsa,isa,jsa,fpm,depsout,loop,dEmin,dEmax,10,dE,dX,M,dres,info)
        print *, 'Iteration count', count, 'Estimate', M
      enddo
      fpm(1)=0
      fpm(4)=20
      fpm(14)=0
      call dfeast_scsrev(UPLO,N,dsa,isa,jsa,fpm,depsout,loop,dEmin,dEmax,M0,dE,dX,M,dres,info)
      print *, 'Number of iterations' , count

   elseif ((PRE=='s').and.(EG=='e')) then 
      print *,'Routine  ','sfeast_scsrev'
      allocate(sE(1:M0))     ! Eigenvalue
      allocate(sres(1:M0))   ! Residual 
      allocate(sX(1:n,1:M0))
      do while (M>M0)
       !!! M is N, but figure out what N is
        sdiff = (sEmax - sEmin) * M0/M
        sEmin = sEmin + sdiff
        call sfeast_scsrev(UPLO,N,ssa,isa,jsa,fpm,sepsout,loop,sEmin,sEmax,M0,sE,sX,M,sres,info)
      enddo
      fpm(1)=1
      fpm(14)=0
      call sfeast_scsrev(UPLO,N,ssa,isa,jsa,fpm,sepsout,loop,sEmin,sEmax,M0,sE,sX,M,sres,info)

   end if

end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! POST-PROCESSING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call system_clock(t2,tim)
  print *,'FEAST OUTPUT INFO',info

  !!!IF ((info==0).or.(info==5).or.(info==6)) then
IF (.true.) then
     print *,'*************************************************'
     print *,'************** REPORT ***************************'
     print *,'*************************************************'

     if (fpm(14)==0) then 
        print *,'Eigenvalues/Residuals'

        if (SHG/='g') then
           print *,'inside interval'
           do i=1,M0

              if ((SHG=='s').and.(PRE=='d')) then
                 print *,i,dE(i),dres(i)
              elseif ((SHG=='s').and.(PRE=='s')) then
                 print *,i,sE(i),sres(i)
              end if
              if (i==M) then 
                print *,''
                print *,'outside interval'
              endif
           enddo
        end if

     elseif (fpm(14)==2) then
        print *,'Running average for estimation (1 -> M0)'
        do i=1,M0
           if ((PRE=='s').or.(PRE=='c')) then
              print *,i,sres(i)
           else
              print *,i,dres(i)
           end if
        enddo
     endif

     print *,'------------------------------------'
     print *,'SIMULATION TIME',(t2-t1)*1.0d0/tim



     
      print *,'# Search interval [Emin,Emax]'
      if ((PRE=='s').or.(PRE=='c')) then
         print *,sEmin,sEmax
      else
         print *,dEmin,dEmax
      end if

     if (fpm(14)==2) then
        print *,'Subspace size M0        ',M0
        print *,'# estimated eigenvalue M',M

     else

        print *,'Subspace size M0    ',M0
        print *,'# eigenvalue found M',M
        print *,'# FEAST iterations  ',loop
        if ((PRE=='s').or.(PRE=='c')) then
           print *,'Relative error on the Trace  ',sepsout
           if (SHG/='g') then
              print *,'Maximum eigenvector residual ',maxval(sres(1:M))
           else
              print *,'Maximum eigenvector residuals',maxval(sres(1:M)),maxval(sres(M0+1:M0+M))
           end if
        else
           print *,'Relative error on the Trace  ',depsout
           if (SHG/='g') then
              print *,'Maximum eigenvector residual ',maxval(dres(1:M))
           else
              print *,'Maximum eigenvector residuals',maxval(dres(1:M)),maxval(dres(M0+1:M0+M))
           end if
        end if


     end if

  end IF

end program driver_feast_sparse






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine dcoo2csr(n,nnz,ic,jc,c,isa,jsa,sa)
  implicit none
  integer :: n,nnz
  integer,dimension(*) :: ic,jc,isa,jsa
  double precision,dimension(*) :: c,sa
!!!
  integer :: k,k1,i,j,idum
  integer,dimension(n) :: iloc
  double precision :: adum


  isa(1:N+1) = 0  
  !find how many elements in row
  do  k=1, nnz
     isa(ic(k)) = isa(ic(k))+1
  end do



  !build isa
  k = 1
  do  i=1,N+1
     k1 = isa(i)
     isa(i) = k
     k = k+k1
  end do


  iloc(1:n)=isa(1:n)
  !Build jsa, sa - increment local row counter
  do  k=1, nnz
     sa(iloc(ic(k))) =  c(k)
     jsa(iloc(ic(k))) = jc(k)
     iloc(ic(k)) = iloc(ic(k))+1
  end do
  ! Reorder by increasing column
  do i=1,n
     do k=isa(i),isa(i+1)-1
        do k1=k,isa(i+1)-1
           if (jsa(k1)<jsa(k)) then
              idum=jsa(k)
              jsa(k)=jsa(k1)
              jsa(k1)=idum
              adum=sa(k)
              sa(k)=sa(k1)
              sa(k1)=adum
           endif
        enddo
     enddo
  enddo


end subroutine dcoo2csr




subroutine scoo2csr(n,nnz,ic,jc,c,isa,jsa,sa)
  implicit none
  integer :: n,nnz
  integer,dimension(*) :: ic,jc,isa,jsa
  real,dimension(*) :: c,sa
!!!
  integer :: k,k1,i,j,idum
  integer,dimension(n) :: iloc
  real :: adum


  isa(1:N+1) = 0  
  !find how many elements in row
  do  k=1, nnz
     isa(ic(k)) = isa(ic(k))+1
  end do

  !build isa
  k = 1
  do  i=1,N+1
     k1 = isa(i)
     isa(i) = k
     k = k+k1
  end do

  iloc(1:n)=isa(1:n)
  !Build jsa, sa - increment local row counter
  do  k=1, nnz
     sa(iloc(ic(k))) =  c(k)
     jsa(iloc(ic(k))) = jc(k)
     iloc(ic(k)) = iloc(ic(k))+1
  end do
  ! Reorder by increasing column
  do i=1,n
     do k=isa(i),isa(i+1)-1
        do k1=k,isa(i+1)-1
           if (jsa(k1)<jsa(k)) then
              idum=jsa(k)
              jsa(k)=jsa(k1)
              jsa(k1)=idum
              adum=sa(k)
              sa(k)=sa(k1)
              sa(k1)=adum
           endif
        enddo
     enddo
  enddo


end subroutine scoo2csr


subroutine dgershgorin(n,nnza,isa,jsa,dsa,dEmin,dEmax)
    !!! The array dsa is of length NNZ and holds all- 
    !!! -the nonzero entries of the matrix(M) in left-to-right top-to-bottom ("row-major") order.
    !!! The array isa is of length n + 1. It's what enables us to iterate the CSR
    !!! jsa, contains the column index in M of each element of dsa and hence is of length NNZ as well.
    !!! dEmin and dEmax represent the variables that store the bound when gersh is complete

    implicit none
    integer :: n,nnza
    integer,dimension(*) :: isa,jsa
    double precision :: dEmin, dEmax
    double precision,dimension(*) :: dsa
    !!!
    integer :: x,z,i,index
    double precision :: radius, center

    z = 1
    x = 1
    radius = 0
    do i = 1, n
        do index = isa(i),isa(i+1)-1
            if (x.EQ.jsa(z)) then
                center = dsa(z)
            else
                radius = radius + ABS(dsa(z))
            endif
            z = z + 1
        enddo
        x = x + 1
        if ((center - radius) < dEmin) then
          dEmin = center - radius
        endif
        if ((center + radius) > dEmax) then
          dEmax = center + radius
        endif
        center = 0
        radius = 0
    enddo
    print *,'Min/Max:    ',dEmin,dEmax
end subroutine dgershgorin

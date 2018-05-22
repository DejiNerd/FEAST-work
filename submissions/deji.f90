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
  double precision,dimension(:),allocatable :: dsa,dsb,dca,dcb,dsa1,radius
  integer,dimension(:),allocatable :: isa,jsa,isb,jsb,ica,jca,icb,jcb,track

  !!!!!!!!!!!!!!!!! Others
  integer :: t1,t2,tim
  integer :: i,j,k,pc
  integer :: M0,M,info,iterations,M1,ms
  character(len=1) :: cc

  double precision :: dEmin,dEmax,dr,minE,ddiff,dEmax0,dEmin0
  real :: sEmin,sEmax,sr,sdiff
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
  read(10,*) dEmin
  read(10,*) dEmax
  read(10,*) M0   ! size subspace
  read(10,*) pc ! Some changes from default for fpm

  do i=1,pc
   read(10,*) j,fpm(j)
  enddo

  close(10)
  M = huge(M)

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
    print *, M0
    call dgershgorin1(UPLO,n,nnza,isa,jsa,dsa,dEmin,dEmax)
    ! call rgersh(UPLO,n,nnza,isa,jsa,dsa,dEmin,dEmax,M0)
  end if

  !!!!!!!!!!!read matrix B if any
  if (EG=='g') then

    open(10,file=trim(name)//'_B.mtx',status='old')
    k=0
    cc='%'
    do while(cc=='%')
      k=k+1 
      read(10,'(A1)') cc
    end do
    close(10)

    open(10,file=trim(name)//'_B.mtx',status='old')
    do i=1,k-1
      read(10,'(A1)') cc
    enddo
    read(10,*) n,n,nnzb
    allocate(icb(nnzb))
    allocate(jcb(nnzb))
    if (PRE=='s') then
      allocate(scb(nnzb))
      do i=1,nnzb
        read(10,*) icb(i),jcb(i),scb(i)
      end do
      elseif (PRE=='d') then
        allocate(dcb(nnzb))
      do i=1,nnzb
        read(10,*) icb(i),jcb(i),dcb(i)
      end do
    end if
    close(10)

    !! create csr format
    allocate(isb(1:n+1))
    allocate(jsb(1:nnzb))


    if (PRE=='s') then
      allocate(ssb(1:nnzb))
      call scoo2csr(n,nnzb,icb,jcb,scb,isb,jsb,ssb)
    elseif (PRE=='d') then
      allocate(dsb(1:nnzb))
      allocate(dsa1(1:nnza))
      allocate(radius(1:n))
      call dcoo2csr(n,nnzb,icb,jcb,dcb,isb,jsb,dsb)
      call condense(UPLO,n,nnzb,isb,jsb,dsb,radius)
    end if
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
      allocate(track(1:3))
      open(10,file="spectrum",status="replace")
      write(10,*) dEmin, dEmax, n
      dEmax0 = dEmax
      dEmin0 = dEmin 
      M1 = M0/2 ! number of eigenvalues
      i = 3
      ms = 10
      if (fpm(56) == 0) then
        dEmax = (dEmax - dEmin)/(2d0**i) + dEmin
      else
        dEmin = dEmax - (dEmax - dEmin)/(2d0**i)
      endif
      call dfeast_scsrev(UPLO,N,dsa,isa,jsa,fpm,depsout,loop,dEmin,dEmax,ms,dE,dX,M,dres,info)
      iterations = 1
      track(1) = M
        if (M == 1) then
          write(10,*) dEmin, dEmax, 0
        else
          write(10,*) dEmin, dEmax, M
        endif
      print *, 'min', dEmin, 'max', dEmax
      print *, 'Iteration count', iterations, 'Estimate', M

      do while ((M < (M0/2 + 0.05*M0)).or.(M > (M0)))
        i = i+1
        ddiff = (dEmax0-dEmin0)/(2d0**i)
        if (M == 1) then
          i =3
          if (fpm(56) == 0) then
            dEmin = dEmax
            dEmax = (dEmax0 - dEmin)/(2d0**i) + dEmin
            print *, 'No eigenvalues found in this range. min = max, max = min + interval/8'
            print*,''
            print*,''
          else
            dEmax = dEmin
            dEmin = dEmax - (dEmax - dEmin0)/(2d0**i)
            print *, 'No eigenvalues found in this range. Max = min, min = max - interval/8'
                        print*,''
            print*,''
          endif
        else if (M > M0) then
          if(fpm(56) == 0) then
            dEmax = dEmax - ddiff
            print *,'Too many eigenvalues found. MAX is reduced by', ddiff
                        print*,''
            print*,''
          else
            dEmin = dEmin + ddiff
            print *,'Too many eigenvalues found. Min shrinks by' ,ddiff
                        print*,''
            print*,''
          endif
        else if (M < (M0/2 + 0.05*M0)) then
          if (M < 0.2*M0/2) then
          endif
          if (fpm(56) == 0) then
            dEmax = dEmax + ddiff
            print *, 'Few eigenvalues found. MAX is increased by', ddiff
                        print*,''
            print*,''
          else
            dEmin = dEmin - ddiff
            print *,' Few eigenvalues found. Min widens by', ddiff
            print*,''
            print*,''
          endif
        endif
        call dfeast_scsrev(UPLO,N,dsa,isa,jsa,fpm,depsout,loop,dEmin,dEmax,ms,dE,dX,M,dres,info)
        iterations = iterations + 1
        if (M == 1) then
          write(10,*) dEmin, dEmax, 0
        else
          write(10,*) dEmin, dEmax, M
        endif
        print *, 'min', dEmin, 'max', dEmax
        print *, 'Iteration count', iterations, 'Estimate', M
      enddo

      write(10,*)
      write(10,*)
      close(10)
      ! stop
      fpm(1)=1
      fpm(4)=20
      fpm(14)=0
      call dfeast_scsrev(UPLO,N,dsa,isa,jsa,fpm,depsout,loop,dEmin,dEmax,M0,dE,dX,M,dres,info)
      print *, 'Number of iterations' , iterations

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


subroutine dgershgorin1(UPLO,n,nnza,isa,jsa,dsa,dEmin,dEmax)
  !!! The array dsa is of length NNZ and holds all- 
    !!! UPLO represents the shape of the matrix (whether full or triangular)
  !!! -the nonzero entries of the matrix(M) in left-to-right top-to-bottom ("row-major") order.
  !!! The array isa is of length n + 1. It's what enables us to iterate the CSR
  !!! jsa, contains the column index in M of each element of dsa and hence is of length NNZ as well.
  !!! dEmin and dEmax represent the variables that store the bound when gersh is complete

  implicit none
  character(len=1) :: UPLO
  integer :: n,nnza
  integer,dimension(*) :: isa,jsa
  double precision :: dEmin, dEmax, theta
  double precision,dimension(*) :: dsa
  double precision, dimension(:), allocatable :: radius, center
  !!!
  integer :: i,index,k
  open(10,file="gershcircle1_g1",status="replace")
  allocate(radius(n))
  allocate(center(n))

  radius = 0d0
  print *, n
  do i = 1, n
    do index = isa(i),isa(i+1)-1
      if (UPLO == 'F') then
        if (i/=jsa(index)) radius(i) = radius(i) + ABS(dsa(index))
        if (i == jsa(index)) center(i) = dsa(index)
      else
        if (i/=jsa(index)) then
          radius(i) = radius(i) + ABS(dsa(index))
          radius(jsa(index)) = radius(jsa(index)) + ABS(dsa(index))
        endif
        if (i == jsa(index)) center(i) = dsa(index)
      endif
    enddo
  enddo
  dEmin = huge(dEmin)
  dEmax = -dEmin
  do i = 1,n  
    if ((center(i) - radius(i)) < dEmin) then
      dEmin = center(i) - radius(i)
    endif
    if ((center(i) + radius(i)) > dEmax) then
      dEmax = center(i) + radius(i)
    endif

    !         ! print *, center,radius
    do k=1,256
      theta=(k-1)*2.0d0*3.14159265359d0/256.0d0
      ! write(10,*) center+sin(theta)*radius,center+cos(theta)*radius
      write(10,*) radius(i)*cos(theta)+center(i), radius(i)*sin(theta)
    enddo
    theta=(1-1)*2.0d0*3.14159265359d0/256.0d0
    ! write(10,*) center+sin(theta)*radius,center+cos(theta)*radius
    write(10,*) radius(i)*cos(theta)+center(i), radius(i)*sin(theta)
    write(10,*)
    write(10,*)
    ! print*, i, center(i), radius(i)
  enddo
  close(10)
  print *,'Min/Max:    ',dEmin,dEmax
  deallocate(radius,center)
end subroutine dgershgorin1

subroutine rgersh(UPLO,n,nnza,isa,jsa,dsa,dEmin,dEmax,M0)
  !!! The array dsa is of length NNZ and holds all- 
  !!! -the nonzero entries of the matrix(M) in left-to-right top-to-bottom ("row-major") order.
  !!! The array isa is of length n + 1. It's what enables us to iterate the CSR
  !!! jsa, contains the column index in M of each element of dsa and hence is of length NNZ as well.
  !!! dEmin and dEmax represent the variables that store the bound when gersh is complete

  implicit none
  character(len=1) :: UPLO
  integer :: n,nnza,M0
  integer,dimension(*) :: isa,jsa
  double precision :: dEmin, dEmax
  double precision,dimension(*) :: dsa
  double precision, dimension(:,:), allocatable :: x,y
  !!!
  integer :: i,j
  double precision, dimension(:),allocatable :: center,radius
  integer,dimension(:),allocatable :: rand_i
  REAL :: num, rand

  allocate(x(n,m0))
  allocate(y(n,m0))
  allocate(rand_i(m0))
  allocate(radius(m0))
  allocate(center(m0))

  print *, 'Performing the reduce gersh subroutine'
  ! CALL srand(86456)
  dEmin = huge(dEmax)
  dEmax = -dEmin
  x = 0d0
  radius = 0d0
  center = 0d0
  call RANDOM_SEED
  do j = 1,m0 
    call RANDOM_NUMBER(num)
    ! print*, j, num
    rand_i(j) = (num * n) + 1
    ! print*, rand_i(j)
    if (rand_i(j) < 1 .or. rand_i(j)>n) stop
    x(rand_i(j),j) = 1.0d0
  enddo
  call dcsrmm(UPLO,'N',n,n,m0,1.0d0,dsa,isa,jsa,x,0.0d0,y)

  do j =1,m0
    do i =1,n
      if (rand_i(j)/=i) then
        radius(j) = radius(j) + ABS(y(i,j))
      else
        center(j) = y(i,j)
      endif
    enddo
  enddo

  do j = 1,m0  
    if ((center(j) - radius(j)) < dEmin) then
      dEmin = center(j) - radius(j)
    endif
    if ((center(j) + radius(j)) > dEmax) then
      dEmax = center(j) + radius(j)
    endif
  enddo


  print *, 'min',dEmin,'max',dEmax
  deallocate(y,x,rand_i)
  stop
end subroutine rgersh

subroutine condense(UPLO,n,nnza,isa,jsa,dsa,radius)
  !!! The array dsa is of length NNZ and holds all- 
  !!! -the nonzero entries of the matrix(M) in left-to-right top-to-bottom ("row-major") order.
  !!! The array isa is of length n + 1. It's what enables us to iterate the CSR
  !!! jsa, contains the column index in M of each element of dsa and hence is of length NNZ as well.
  !!! dEmin and dEmax represent the variables that store the bound when gersh is complete

  implicit none
  character(len=1) :: UPLO
  integer :: n,nnza
  integer,dimension(*) :: isa,jsa
  double precision,dimension(*) :: dsa,radius
  ! double precision, dimension(:), allocatable :: radius
  !!!
  integer :: i,index
  ! allocate(radius(n))
  open(10,file="condensed_mtx",status="replace")

  radius(1:n) = 0d0
  do i = 1, n
    do index = isa(i),isa(i+1)-1
      if (UPLO == 'F') then
        radius(i) = radius(i) + ABS(dsa(index))
      else
        radius(i) = radius(i) + ABS(dsa(index))
        if (i/=jsa(index)) radius(jsa(index)) = radius(jsa(index)) + ABS(dsa(index))
      endif
    enddo
  enddo
  do i = 1,n  
    write(10,*) radius(i)
  enddo
  close(10)
  ! deallocate(radius)
end subroutine condense

subroutine rquotientsev(UPLO,n,nnza,isa,jsa,dsa,dEmin,dEmax)
  !!! The array dsa is of length NNZ and holds all- 
  !!! -the nonzero entries of the matrix(M) in left-to-right top-to-bottom ("row-major") order.
  !!! The array isa is of length n + 1. It's what enables us to iterate the CSR
  !!! jsa, contains the column index in M of each element of dsa and hence is of length NNZ as well.
  !!! dEmin and dEmax represent the variables that store the bound when gersh is complete

  implicit none
  character(len=1) :: UPLO
  integer :: n,nnza
  integer,dimension(*) :: isa,jsa
  double precision :: dEmin, dEmax
  double precision,dimension(*) :: dsa
  double precision, dimension(:), allocatable :: sigma
  double precision, dimension(:,:), allocatable :: x,y
  !!!
  integer :: i,j,m0

  m0 = 10
  allocate(sigma(m0))
  allocate(x(n,m0))
  allocate(y(n,m0))

  call dlarnv(1, (/100,200,40,7/), n*m0, x)
  do j=1,m0
    do i=1,n
      if (abs(x(i,j)) > 0.1d0) then
       x(i,j)=0d0
     else
      x(i,j)=1d0
    endif
  enddo
enddo

call dcsrmm(UPLO,'N',n,n,m0,1.0d0,dsa,isa,jsa,x,0.0d0,y)
do i=1,m0
  sigma(i) = dot_product(x(:,i),y(:,i))
enddo

do i=1,m0
  sigma(i) =sigma(i)/dot_product(x(:,i),x(:,i))
  print *, i, sigma(i)
enddo

dEmin = minval(sigma)
dEmax = maxval(sigma)
print *, 'min',dEmin,'max',dEmax
deallocate(sigma,y,x)

end subroutine rquotientsev

subroutine rquotients(UPLO,n,nnza,isa,jsa,dsa,isb,jsb,dsb,dEmin,dEmax)
  !!! The array dsa is of length NNZ and holds all- 
  !!! -the nonzero entries of the matrix(M) in left-to-right top-to-bottom ("row-major") order.
  !!! The array isa is of length n + 1. It's what enables us to iterate the CSR
  !!! jsa, contains the column index in M of each element of dsa and hence is of length NNZ as well.
  !!! dEmin and dEmax represent the variables that store the bound when gersh is complete

  implicit none
  character(len=1) :: UPLO
  integer :: n,nnza
  integer,dimension(*) :: isa,jsa,isb,jsb
  double precision :: dEmin, dEmax
  double precision,dimension(*) :: dsa, dsb
  double precision, dimension(:), allocatable :: sigma
  double precision, dimension(:,:), allocatable :: x,y
  !!!
  integer :: i,j,m0

  m0 = 10
  allocate(sigma(m0))
  allocate(x(n,m0))
  allocate(y(n,m0))

  call dlarnv(1, (/100,200,40,7/), n*m0, x)
  do j=1,m0
    do i=1,n
      if (abs(x(i,j)) > 0.1d0) then
       x(i,j)=0d0
     else
      x(i,j)=1d0
    endif
  enddo
enddo
call dcsrmm(UPLO,'N',n,n,m0,1.0d0,dsa,isa,jsa,x,0.0d0,y)
do i=1,m0
  sigma(i) = dot_product(x(:,i),y(:,i))
enddo

call dcsrmm(UPLO,'N',n,n,m0,1.0d0,dsb,isb,jsb,x,0.0d0,y)
do i=1,m0
  sigma(i) =sigma(i)/dot_product(x(:,i),y(:,i))
enddo

dEmin = minval(sigma)
dEmax = maxval(sigma)

print *, 'min',dEmin,'max',dEmax
deallocate(sigma,y,x)

end subroutine rquotients


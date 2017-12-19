program gershgorin

    implicit none

    integer :: N,temp,center,x,z,i,index
    double precision :: sum
    double precision,dimension(:),allocatable :: val, gc_ctr, gc_sum
    integer,dimension(:),allocatable :: row, row_ptr, col, index_ptr

    ! allocate(val(1:8))
    ! allocate(index_ptr(1:6))
    ! allocate(col(1:8))

    !!!! test case 1 !!!!
    ! val = (/10,20,30,40,50,60,70,80/)
    ! col = (/1,2,2,4,3,4,5,6/)
    ! index_ptr = (/0,2,4,7,8/)

    !!!! test case 2 !!!!
    ! val = (/10,0,30,10,50,40,10,0,10,50,10,60/)
    ! col = (/1,2,2,4,5,3,5,1,2,4,4,5/)
    ! index_ptr = (/0,2,5,7,10,12/)


    ! allocata(gc_ctr(1:N_rows))
    ! allocata(gc_sum(1:N_rows))


    ! temp = 1
    ! sum = 0
    ! do i=1,N
    !     if(row(i).NE.temp) then
    !         print *, center, sum
    !         gc_ctr(i-1) = center
    !         gc_sum(i-1) = sum
    !         temp = row(i)
    !         center = 0
    !         sum = 0
    !     end if
    !     if (row(i) == col(i)) then
    !         center = val(i)
    !     else
    !         sum = sum + val(i)
    !     end if
    ! end do

    N = 6 ! Size of index_ptr
    z = 1
    x = 1
    sum = 0
    do i = 1,N-1
        do index = index_ptr(i),index_ptr(i+1)-1
            if (x.EQ.col(z)) then
                center = val(z)
            else
                sum = sum + val(z)
            endif
            z = z + 1
        enddo
        x = x + 1
        print *, center, sum
        center = 0
        sum = 0
    enddo
end program gershgorin
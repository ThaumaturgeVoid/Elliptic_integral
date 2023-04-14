program main
    use module_calc_ellint
    use moduleMELPsubroutine
    use moduleMELPfunction
    use, intrinsic :: iso_c_binding
    implicit none 
    real(c_double) :: k
    real(c_double) :: fc,ec,ff,ef,fs,es
    integer :: i,n
    real(8) ::CPU1,CPU2,time
    call random_seed
    ! 用于测试准确性
    do i=0,10
        k=0.1*i
        call ellint(k,fc,ec)
        call melp(k,PI/2.d0,fs,es)
        write(*,"(5X,A,6(A20))") "k","F(k)(C++)","F(k)(F90func)","F(k)(F90subr)","E(k)(C++)","E(k)(F90func)","E(k)(F90subr)"
        write(*,"(F8.1,6(E20.10))") k,Fc,melp1(k,PI/2.d0),fs,Ec,melp2(k,PI/2.d0),es
        ! write(*,*) "F(k)(C++)=",F,"E(k)(C++)=",E
        ! write(*,*) "F(k)(F90)=",melp1(k,PI/2.d0),"E(k)(F90)=",melp2(k,PI/2.d0)
    end do

    ! 用于测试速度
    n=10000000
    write(*,"(A30, I)") "test loop:",n
    call cpu_time(CPU1)
    do i=1,n
        call random_number(k)
        call ellint(k,fc,ec)
    end do
    call cpu_time(CPU2)
    time=cpu2-CPU1
    write(*,"(A30, F)") "C++ Time:", time

    call cpu_time(CPU1)
    do i=1,n
        call random_number(k)
        ff=melp1(k,PI/2.d0)
        ef=melp2(k,PI/2.d0)
    end do
    call cpu_time(CPU2)
    time=cpu2-CPU1
    write(*,"(A30, F)") "F90 Function Time:", time

    call cpu_time(CPU1)
    do i=1,n
        call random_number(k)
        call melp(k,PI/2.d0,fs,es)
    end do
    call cpu_time(CPU2)
    time=cpu2-CPU1
    write(*,"(A30, F)") "F90 Subroutine Time:", time


end program main

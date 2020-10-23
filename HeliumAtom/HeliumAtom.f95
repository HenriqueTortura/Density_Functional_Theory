program Teste

    use util

    implicit none

    ! Setting main parameters
    real (kind = 8), parameter, dimension(2) :: r_range = (/ 0., 30./) , Eigenvalue_Range = (/ -5., 0./)

    integer, parameter :: KS_int_max = 100, SelfCons_int_max = 100

    real (kind = 8), parameter :: Eigenvalue_tol = 0.001, u0_tol = 0.001, SelfCons_tol = 0.001

    logical, parameter :: Uniform = .FALSE.

    real (kind = 8), parameter :: h = 0.0001

    integer, parameter :: j_max = 200000
    real (kind = 8), parameter :: delta = 0.0001

    real :: start, finish

    !---------------------------------------------------------------------------------------!

    call cpu_time(start)

    call HeliumAtom(r_range, Eigenvalue_range, SelfCons_int_max, KS_int_max,&
    &Eigenvalue_tol, u0_tol, SelfCons_tol, Uniform, h, j_max, delta)

    call cpu_time(finish)
    print *,"CPU time = ",finish-start,"s"

end program


program Helium

    use dft

    implicit none

    ! Setting main parameters
    real (kind = 8), parameter, dimension(2) :: r_range = (/ 0., 20./), Eigenvalue_Range = (/ -5., 0./)

    integer, parameter :: KS_int_max = 100, SelfCons_int_max = 100

    real (kind = 8), parameter :: Eigenvalue_tol = 0.0001, u0_tol = 0.0001, SelfCons_tol = 0.0001

    logical, dimension(2), parameter :: Uniform_Numerov = (/ .FALSE., .FALSE./)

    real (kind = 8), parameter :: h = 0.0001

    integer, parameter :: j_max = 200000
    real (kind = 8), parameter :: delta = 0.0001

    real :: start, finish

    !---------------------------------------------------------------------------------------!

    call cpu_time(start)

    call HeliumAtom(r_range, Eigenvalue_range, SelfCons_int_max, KS_int_max,&
    &Eigenvalue_tol, u0_tol, SelfCons_tol, Uniform_Numerov, h, j_max, delta)

    call cpu_time(finish)
    print *,"CPU time = ",finish-start,"s"

end program

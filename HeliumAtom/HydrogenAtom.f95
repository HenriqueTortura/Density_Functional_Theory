program Hydrogen

    use dft

    implicit none

    ! Setting main parameters
    real (kind = 8), dimension(2) :: r_range = (/ 0., 50./), Eigenvalue_Range = (/-5., 0./)

    integer, parameter :: KS_int_max = 100

    real (kind = 8) , parameter :: eigenvalue_tol = 0.00001, u0_tol= 0.001

    logical, dimension(2), parameter :: Uniform_Numerov = (/ .FALSE., .TRUE./)

    real (kind = 8), parameter :: h = 0.0001

    integer, parameter :: j_max = 1000
    real (kind = 8), parameter :: delta = 0.0001

    logical, parameter :: write_data = .FALSE.

    real :: start, finish

    !---------------------------------------------------------------------------------------!

    call cpu_time(start)

    call HydrogenAtom(r_range, Eigenvalue_range, KS_int_max, Eigenvalue_tol, u0_tol,&
    &Uniform_Numerov, h, j_max, delta, write_data)

    call cpu_time(finish)
    print *,"CPU time = ",finish-start,"s"

end program

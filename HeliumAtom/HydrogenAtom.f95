program Hydrogen

    use util

    implicit none

    ! Setting main parameters
    real (kind = 8), dimension(2) :: r_range = (/ 0., 30./), Eigenvalue_Range = (/-2., -0./)

    integer, parameter :: KS_int_max = 100

    real (kind = 8) , parameter :: eigenvalue_tol = 0.00001, u0_tol= 0.001

    logical, dimension(2), parameter :: Uniform_Numerov = (/ .FALSE., .TRUE./)

    real (kind = 8), parameter :: h = 0.0001

    integer, parameter :: j_max = 200000
    real (kind = 8), parameter :: delta = 0.0001


    call HydrogenAtom(r_range, Eigenvalue_range, KS_int_max, Eigenvalue_tol, u0_tol, Uniform_Numerov, h, j_max, delta)


end program

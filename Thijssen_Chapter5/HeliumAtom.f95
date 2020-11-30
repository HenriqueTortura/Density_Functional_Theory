program Helium

    use dft

    implicit none

    ! Setting main parameters
    real (kind = 8), parameter, dimension(2) :: r_range = (/ 0., 20./), Eigenvalue_Range = (/ -5., 0./)

    integer, parameter :: KS_int_max = 100, SelfCons_int_max = 100,&
    &Z = 2, N_electrons = 2, Correlation_Method = 1

    real (kind = 8), parameter :: Eigenvalue_tol = 0.0001, u0_tol = 0.0001, SelfCons_tol = 0.0001

    logical, dimension(2), parameter :: Uniform_Numerov = (/ .FALSE., .TRUE./)

    real (kind = 8), parameter :: h = 0.0001

    integer, parameter :: j_max = 200000
    real (kind = 8), parameter :: delta = 0.0001

    real (kind = 8) :: Energy, Hartree_Correction, Exchange_Correction, Correlation_Correction

    real :: start, finish

    !---------------------------------------------------------------------------------------!

    call cpu_time(start)

    call HeliumAtom(r_range, Eigenvalue_range, SelfCons_int_max, KS_int_max,&
    &Eigenvalue_tol, u0_tol, SelfCons_tol, Uniform_Numerov, h, j_max, delta,&
    &Correlation_Method, Z, N_electrons, Energy, Hartree_Correction,&
    &Exchange_Correction, Correlation_Correction)

    call cpu_time(finish)
    print *,"CPU time = ",finish-start,"s"

end program


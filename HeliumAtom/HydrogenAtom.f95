program HydrogenAtom

    use util

    implicit none

    ! Setting main parameters
    integer, parameter :: r0 = 0, r_max = 50, int_max = 100 ! Radial coordinates range
    real (kind = 8) , parameter :: h = 0.0001, tolerance = 0.00001 ! Discretization

    real (kind = 8), dimension(1:2) :: Eigenvalue_Range = (/-1., -0.2/) ! Eigenvalues range

    integer, parameter :: n = int((r_max-r0)/h) ! Number of steps in the radial coordinate

    integer i
    real (kind = 8), dimension(1:n) :: r, Potential, u, Potential_U
    real (kind = 8) :: Eigenvalue

    ! Initializing radial coordinates and effective potential
    do i=1, n
        r(i) = r0 + h*i
        Potential(i) = -1/r(i)
    end do

    ! Solve Kohn-Sham
    call KohnSham1D(r, u, Potential, Eigenvalue_Range, EigenValue, h, int_max, tolerance)

    call Poisson(r, u, Potential_U, h)

    ! To plot u(r) and Potential U
    open(1, file='Hydrogen_u.dat', status='replace')
        do i=1, n
            write(1,*) r(i), u(i)
        end do
    close(1)
    open(2, file='Hydrogen_Potential_U.dat', status='replace')
        do i=1, n
            write(2,*) r(i), Potential_U(i)
        end do
    close(2)

end program

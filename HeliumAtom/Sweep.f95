program Sweep

    use util

    implicit none

    ! Setting main parameters
    integer, parameter :: r0 = 0, r_max = 50, int_max = 100 ! Radial coordinates range
    real (kind = 8) , parameter :: h = 0.0001, tolerance = 0.00001 ! Discretization

    real (kind = 8), dimension(1:2) :: Eigenvalue_Range = (/-1., 5./)
    integer, parameter :: n_eigenvalues = 10001

    integer, parameter :: n = int((r_max-r0)/h) ! Number of steps

    integer i
    real (kind = 8), dimension(1:n) :: r, Potential, u, Potential_U

    real (kind = 8), dimension(1:n_eigenvalues) :: Eigenvalues, u0s

    ! Initializing radial coordinates and effective potential
    do i=1, n
        r(i) = r0 + h*i
        Potential(i) = -1/r(i)
    end do

    ! Initializing eigenvalues
    do i=1, n_eigenvalues
        Eigenvalues(i) = Eigenvalue_Range(1) + (i-1)*(Eigenvalue_Range(2)-Eigenvalue_Range(1))/(n_eigenvalues-1)
    end do

    call KohnShamSweep(r, u, Potential, h, Eigenvalues, u0s)

    ! To plot u0s
    open(1, file='Hydrogen_u0s.dat', status='replace')
        do i=1, n_eigenvalues
            write(1,*) Eigenvalues(i), u0s(i)
        end do
    close(1)

end program


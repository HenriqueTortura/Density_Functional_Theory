program Sweep

    use util

    implicit none

    integer i
    integer, parameter :: r0 = 0, r_max = 50, int_max = 100
    real (kind = 8) , parameter :: h = 0.0001 ! Discretization
    integer, parameter :: n = int((r_max-r0)/h) ! Number of steps
    integer, parameter :: n_eigenvalues = 10001

    real (kind = 8) , dimension(1:n) :: r, Potential, u, Potential_U

    real (kind = 8), dimension(1:n_eigenvalues) :: EigenValues, u0s
    real (kind = 8)  :: EigenValue_min, EigenValue_max

    do i=1, n
        r(i) = r0 + h*i
        Potential(i) = -1/r(i)
    end do

    EigenValue_min = -1
    EigenValue_max = 5

    do i=1, n_eigenvalues
        EigenValues(i) = EigenValue_min + (i-1)*(EigenValue_max-EigenValue_min)/(n_eigenvalues-1)
    end do

    call KohnShamSweep(r, u, Potential, h, EigenValues, u0s)

    ! To plot u0s
    open(1, file='u0s.dat', status='replace')
        do i=1, n_eigenvalues
            write(1,*) EigenValues(i), u0s(i)
        end do
    close(1)

end program


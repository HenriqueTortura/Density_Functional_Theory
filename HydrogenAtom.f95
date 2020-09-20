program HydrogenAtom

    use util

    implicit none

    integer i
    integer, parameter :: r0 = 0, r_max = 20, int_max = 100
    real (kind = 8) , parameter :: h = 0.0004 ! Discretization
    integer, parameter :: n = int((r_max-r0)/h) ! Number of steps

    real (kind = 8) , dimension(1:n) :: r, Potential, u, Potential_U

    real (kind = 8)  :: EigenValue_min, EigenValue_max, E_aux, tolerance

    do i=1, n
        r(i) = r0 + h*i
        Potential(i) = -1/r(i)
    end do

    tolerance = 0.00001

    u(n) = r(n)*exp(-r(n))
    u(n-1) = r(n-1)*exp(-r(n-1))

    EigenValue_min = -10
    EigenValue_max = 0

    call KohnSham1D(u, Potential, h, int_max, EigenValue_min, EigenValue_max, tolerance)

    call Poisson(r, u, Potential_U, h)

    ! To plot Potential U
    open(1, file='Potential_U.dat', status='replace')
        do i=1, n
            write(1,*) r(i), Potential_U(i)
        end do
    close(1)

end program

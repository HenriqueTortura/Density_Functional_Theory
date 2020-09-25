program HydrogenAtom

    use util

    implicit none

    integer i
    integer, parameter :: r0 = 0, r_max = 50, int_max = 100
    real (kind = 8) , parameter :: h = 0.0001 ! Discretization
    integer, parameter :: n = int((r_max-r0)/h) ! Number of steps

    real (kind = 8) , dimension(1:n) :: r, Potential, u, Potential_U

    real (kind = 8)  :: EigenValue_min, EigenValue_max, tolerance

    do i=1, n
        r(i) = r0 + h*i
        Potential(i) = -1/r(i)
    end do

    tolerance = 0.00001

    EigenValue_min = -20
    EigenValue_max = 0

    call KohnSham1D(r, u, Potential, h, int_max, EigenValue_min, EigenValue_max, tolerance)

    call Poisson(r, u, Potential_U, h)

    ! To plot u(r) and Potential U
    open(1, file='u.dat', status='replace')
        do i=1, n
            write(1,*) r(i), u(i)
        end do
    close(1)
    open(2, file='Potential_U.dat', status='replace')
        do i=1, n
            write(2,*) r(i), Potential_U(i)
        end do
    close(2)

end program

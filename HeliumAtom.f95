program HeliumAtom

    use util

    implicit none

    integer i
    integer, parameter :: r0 = 0, r_max = 20, int_max = 100
    real (kind = 8) , parameter :: h = 0.00004 ! Parâmetro da discretização
    integer, parameter :: n = int((r_max-r0)/h) ! Número de intervalos
    real (kind = 8) , dimension(1:n) :: r, Potential, u
    real (kind = 8)  :: EigenValue_min, EigenValue_max, E_aux, tolerance

    ! Inicializando as posições e o Potencial
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

end program

program HeliumAtom

    use util

    implicit none


    integer i
    integer, parameter :: r0 = 0, r_max = 20, int_max = 10
    real (kind = 8) , parameter :: h = 0.0004 ! Parâmetro da discretização
    integer, parameter :: n = int((r_max-r0)/h) ! Número de intervalos
    real (kind = 8) , dimension(1:n) :: r, f, Potential, u
    real (kind = 8)  :: EigenValue0, EigenValue, u0, E_aux

    ! Inicializando as posições e o Potencial
    do i=1, n
        r(i) = r0 + h*i
        Potential(i) = -1/r(i)
    end do

    ! Chutes iniciais

    u0 = h
    u(n) = r(n)*exp(-r(n))
    u(n-1) = r(n-1)*exp(-r(n-1))
    EigenValue0 = -0.6
    EigenValue = -0.4

    call RadialSchrodinger(u, Potential, h, int_max, EigenValue0, EigenValue)


end program

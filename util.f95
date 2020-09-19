module util
    implicit none


contains

    subroutine RadialSchrodinger(u, Potential, h, int_max, EigenValue0, EigenValue)


        integer i, j, int_max
        real (kind = 8) , dimension(:) :: Potential, u
        integer :: n
        real (kind = 8)  :: EigenValue0, EigenValue, u0, E_aux, h

        n = size(u)

        u = Numerov(h, n, u, -2*(EigenValue0-Potential))
        u0 = u(1)
        print *,EigenValue0
        print *,u(1)

        do i=1, int_max

            u = Numerov(h, n, u, -2*(EigenValue0-Potential))

            E_aux = EigenValue
            print *,EigenValue
            EigenValue = (EigenValue0*u(1) - EigenValue*u0)/(u(1)-u0)

            u0 = u(1)
            EigenValue0 = E_aux
        end do

    contains

        function Numerov(h, n, u, f)
            integer :: i, n
            real (kind = 8)  :: h, a, b, c
            real (kind = 8) , dimension(1:n) :: u, f, Numerov

            do i=2, n-1
                a = 2*(1 + 5*(h**2)*f(n-i+1)/12)
                b = (1 - (h**2)*f(n-i+2)/12)
                c = (1 - (h**2)*f(n-i)/12)
                u(n-i) = (a*u(n-i+1) - b*u(n-i+2)) / c
            end do

            Numerov = u
        end function Numerov

    end subroutine RadialSchrodinger

end module util

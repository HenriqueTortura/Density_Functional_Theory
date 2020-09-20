module util
    implicit none

contains

    subroutine KohnSham1D(u, Potential, h, int_max, EigenValue_min, EigenValue_max, tolerance)

        integer i, j, int_max
        real (kind = 8) , dimension(:) :: Potential, u
        integer :: n
        real (kind = 8) :: EigenValue_min, EigenValue_max, EigenValue, E_aux, h, tolerance

        n = size(u)
        i = 1
        EigenValue = (EigenValue_min + EigenValue_max)/2

        ! Main loop for finding the Kohn-Sham eigenvalue
        do while (abs(EigenValue-EigenValue_min)>=tolerance .and. i<int_max)

            ! Integrating wave function via Numerov
            u = Numerov(h, n, u, -2*(EigenValue-Potential))

            print *,"Eigen Value tested:", EigenValue
            print *,"Boundary term:", u(1)
            print *,"*-----------------------*"

            ! Bisection method approaching boundary condition
            if (u(1)<0) then
                EigenValue_max = EigenValue
            else
                EigenValue_min = EigenValue
            end if
            EigenValue = (EigenValue_min + EigenValue_max)/2

            i = i + 1

        end do

        ! Satisfying tolerance for final result
        if (abs(EigenValue-EigenValue_min)<tolerance) then
            print *,"Converged for eigenvalue", EigenValue
        else
            print *, "Did not conveged within ", int_max, " iterations"
        end if

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

    end subroutine KohnSham1D

end module util

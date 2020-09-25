module util
    implicit none

contains

    subroutine KohnSham1D(r, u, Potential, h, int_max, EigenValue_min, EigenValue_max, tolerance)
    ! Routine to solve a one dimentional Kohn-Sham (Schrödinger) equation.
    ! Parameters
    !   r: 1D array of coordinates to run by;
    !   Potential: Effective potential of Kohn-Sham equation;
    !   h: Discretization of array r;
    !   int_max: Maximum number of iterations;
    !   EigenValue_min & EigenValue_max: boundaries for possible eigenvalue;
    !   tolerance: Maximum difference between eigenvalue in each iteration to determinate convergency.
    ! Returns
    !   u: Radial wave function u(r) = r*\psi(r);
    !   EigenValue: Eigenvalue of function u(r) satisfying boundary conditions.

        integer i, int_max
        real (kind = 8), dimension(:) :: Potential, u, r
        integer :: n
        real (kind = 8) :: EigenValue_min, EigenValue_max, EigenValue, E_aux, h, tolerance

        ! Arrays size and interation counter
        n = size(u)
        i = 1

        ! Initial wave function guess
        u(n) = r(n)*exp(-r(n))
        u(n-1) = r(n-1)*exp(-r(n-1))

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

        ! Normalizing u
        u = u/sqrt(sum(u*u*h))

    contains

        function Numerov(h, n, u, f)
        ! Applies Numerov's method for the wave function.
        ! Parameters
        !   h: Discretization;
        !   n: Number of steps, i.e., size of array;
        !   (part of) u: Two guesses at boundary of wave function on the form u(r) = r*\psi(r);
        !   f: f function on Numerov's method.
        ! Returns
        !   u: Integrated wave function.

            integer :: i, n
            real (kind = 8) :: h, a, b, c
            real (kind = 8), dimension(1:n) :: u, f, Numerov

            do i=2, n-1
                a = 2*(1 + 5*(h**2)*f(n-i+1)/12)
                b = (1 - (h**2)*f(n-i+2)/12)
                c = (1 - (h**2)*f(n-i)/12)
                u(n-i) = (a*u(n-i+1) - b*u(n-i+2)) / c
            end do

            Numerov = u

        end function Numerov

    end subroutine KohnSham1D

    subroutine Poisson(r, u, Potential_U, h)

        real (kind = 8), dimension(:) :: r, u, Potential_U
        real (kind = 8) :: h, a
        integer :: n, i

        n = size(r)

        ! Boundary condition at r=0 and first guess
        Potential_U(1) = 0
        Potential_U(2) = h

        ! Integrating Potential U via Verlet
        Potential_U = Verlet(Potential_U, -u**2/r, h, n)

        ! Satisfaying boundary condition at r_max by adding term a*r
        a = (1 - Potential_U(size(Potential_U))) / r(size(r))
        Potential_U = Potential_U + a*r

    contains
        function Verlet(Potential_U, F, h, n)
            integer :: i, n
            real (kind = 8), dimension(1:n) :: F, Potential_U, Verlet
            real (kind = 8) :: h

            do i=3, n
                Potential_U(i) = 2*Potential_U(i-1) - Potential_U(i-2) +h**2*F(i-1)
            end do

            Verlet = Potential_U

        end function Verlet

    end subroutine Poisson

!! Fooling around
    subroutine KohnShamSweep(r, u, Potential, h, EigenValues, u0s)

        integer i
        real (kind = 8), dimension(:) :: Potential, u, r
        real (kind = 8), dimension(:) :: EigenValues, u0s
        integer :: n, n_eigenvalues
        real (kind = 8) :: h

        ! Arrays size and interation counter
        n = size(u)
        n_eigenvalues = size(EigenValues)
        i = 1

        ! Initial wave function guess
        u(n) = r(n)*exp(-r(n))
        u(n-1) = r(n-1)*exp(-r(n-1))

        ! Main loop for finding the Kohn-Sham eigenvalue
        do while (i <= n_eigenvalues)

            ! Integrating wave function via Numerov
            u = Numerov(h, n, u, -2*(EigenValues(i)-Potential))

            ! Normalizing u
            u = u/sqrt(sum(u*u*h))
            u0s(i) = u(1)

            print *,"Eigen Value tested:", EigenValues(i)
            print *,"Boundary term:", u(1)
            print *,"*-----------------------*"

            i = i + 1

        end do

    contains

        function Numerov(h, n, u, f)
        ! Applies Numerov's method for the wave function.
        ! Parameters
        !   h: Discretization;
        !   n: Number of steps, i.e., size of array;
        !   (part of) u: Two guesses at boundary of wave function on the form u(r) = r*\psi(r);
        !   f: f function on Numerov's method.
        ! Returns
        !   u: Integrated wave function.

            integer :: i, n
            real (kind = 8) :: h, a, b, c
            real (kind = 8), dimension(1:n) :: u, f, Numerov

            do i=2, n-1
                a = 2*(1 + 5*(h**2)*f(n-i+1)/12)
                b = (1 - (h**2)*f(n-i+2)/12)
                c = (1 - (h**2)*f(n-i)/12)
                u(n-i) = (a*u(n-i+1) - b*u(n-i+2)) / c
            end do

            Numerov = u

        end function Numerov

    end subroutine KohnShamSweep

end module util

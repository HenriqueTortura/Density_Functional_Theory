module util
    implicit none

contains

    subroutine KohnSham1D(r, u, Potential, Eigenvalue_Range, EigenValue, h, int_max, eigenvalue_tol, u0_tol)
    ! Routine to solve a one dimentional Kohn-Sham (Schrödinger) equation.
    ! Parameters
    !   r: 1D array of coordinates to run by;
    !   Potential: Effective potential of Kohn-Sham equation;
    !   Eigenvalue_Range: boundaries for possible eigenvalue;
    !   h: Discretization of array r;
    !   int_max: Maximum number of iterations;
    !   eigenvalue_tol: Maximum difference between eigenvalue in each iteration to determinate convergency.
    !   u0_tol: Maximum absolute value for u at 'origin', expecting to converge for u(1) = 0
    ! Returns
    !   u: Radial wave function u(r) = r*\psi(r);
    !   EigenValue: Eigenvalue of function u(r) satisfying boundary conditions.

        integer i, int_max
        real (kind = 8), dimension(:) :: Potential, u, r
        real (kind =8), dimension(2) :: Eigenvalue_Range
        integer :: n !Dimension of Potential, u and r
        real (kind = 8) :: EigenValue, E_aux, h, eigenvalue_tol, u0_tol

        ! Arrays size and interation counter
        n = size(u)

        ! Initial wave function guess
        u(n) = r(n)*exp(-r(n))
        u(n-1) = r(n-1)*exp(-r(n-1))
        u(1) = abs(2*u0_tol)

        ! Main loop for finding the Kohn-Sham eigenvalue
        do i=1, int_max, 1

            EigenValue = (Eigenvalue_Range(1) + Eigenvalue_Range(2))/2

            ! Integrating wave function via Numerov
            u = Numerov(h, n, u, -2*(EigenValue-Potential))

            print *,"Eigen Value tested:", EigenValue
            print *,"Boundary term:", u(1)
            print *,"*-----------------------*"

            ! Bisection method approaching boundary condition
            if ((abs(u(1))>=u0_tol .or. abs(EigenValue-Eigenvalue_Range(1))>=eigenvalue_tol)) then
                if (u(1)<0) then
                    Eigenvalue_Range(2) = EigenValue
                else
                    Eigenvalue_Range(1) = EigenValue
                end if
            else
                exit
            end if

        end do

        ! Satisfying tolerance for final result
        if ((abs(u(1))<u0_tol .and. abs(EigenValue-Eigenvalue_Range(1))<eigenvalue_tol)) then
            print *,"Converged for eigenvalue: ", EigenValue
        else
            print *, "Did not conveged within ", (i-1), " iterations"
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
    ! Routine to solve a Poisson equation of radial potential Hartree potential.
    ! Parameters
    !   r: 1D array of coordinates to run by;
    !   u: Radial wave function u(r) = r*\psi(r);
    !   h: Discretization of array r;
    ! Returns
    !   Potential_U: U radial function, related with Hartree Potential by Hartree = r*U(r).

        real (kind = 8), dimension(:) :: r, u, Potential_U
        real (kind = 8) :: h, a
        integer :: n, i

        n = size(r)

        ! Boundary condition at r=0 and first guess
        Potential_U(1) = r(1)
        Potential_U(2) = r(2)

        ! Integrating Potential U via Verlet
        Potential_U = Verlet(Potential_U, -u**2/r, h, n)

        ! Satisfaying boundary condition at r_max by adding term a*r
        a = (1 - Potential_U(size(Potential_U))) / r(size(r))
        Potential_U = Potential_U + a*r

    contains
        function Verlet(Potential_U, F, h, n)
        ! Applies Verlet integration for Potential_U.
        ! Parameters
        !   h: Discretization;
        !   n: Number of steps, i.e., size of array;
        !   (part of) Potential_U: Two guesses at boundary on the form Potential_U(r) = r;
        !   F: F function on Verlet algorithm.
        ! Returns
        !   Potential_U: Integrated potential.

            integer :: i, n
            real (kind = 8), dimension(1:n) :: F, Potential_U, Verlet
            real (kind = 8) :: h

            do i=3, n
                Potential_U(i) = 2*Potential_U(i-1) - Potential_U(i-2) + (h**2)*F(i-1)
            end do

            Verlet = Potential_U

        end function Verlet

    end subroutine Poisson

!! Fooling around
    subroutine KohnShamSweep(r, u, Potential, h, Eigenvalues, u0s)

        integer i
        real (kind = 8), dimension(:) :: Potential, u, r
        real (kind = 8), dimension(:) :: Eigenvalues, u0s
        integer :: n, n_eigenvalues
        real (kind = 8) :: h

        ! Arrays size and interation counter
        n = size(u)
        n_eigenvalues = size(Eigenvalues)
        i = 1

        ! Initial wave function guess
        u(n) = r(n)*exp(-r(n))
        u(n-1) = r(n-1)*exp(-r(n-1))

        ! Main loop for finding the Kohn-Sham eigenvalue
        do while (i <= n_eigenvalues)

            ! Integrating wave function via Numerov
            u = Numerov(h, n, u, -2*(Eigenvalues(i)-Potential))

            ! Normalizing u
            u = u/sqrt(sum(u*u*h))
            u0s(i) = u(1)

            !print *,"Eigen Value tested:", Eigenvalues(i)
            !print *,"Boundary term:", u(1)
            !print *,"*-----------------------*"

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

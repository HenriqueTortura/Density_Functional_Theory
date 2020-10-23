module util

    implicit none

contains

    subroutine KohnSham1D(r, u, Potential, Eigenvalue_Range, EigenValue,&
    &int_max, eigenvalue_tol, u0_tol, n, h, rp, delta, Uniform)
    ! Routine to solve a one dimentional Kohn-Sham (Schrödinger) equation.
    ! Parameters
    !   r: 1D array of coordinates to run by;
    !   Potential: Effective potential of Kohn-Sham equation;
    !   Eigenvalue_Range: boundaries for possible eigenvalue;
    !   int_max: Maximum number of iterations;
    !   eigenvalue_tol: Maximum difference between eigenvalue in each iteration to determinate convergency.
    !   u0_tol: Maximum absolute value for u at 'origin', expecting to converge for u(1) = 0;
    !   h: Discretization of array r (in case of r being uniform);
    !
    ! Returns
    !   u: Radial wave function u(r) = r*\psi(r);
    !   EigenValue: Eigenvalue of function u(r) satisfying boundary conditions.

        integer :: i, int_max, n !Dimension of Potential, u and r
        real (kind = 8), dimension(n) :: Potential, u, r, du, Delta_r
        real (kind =8), dimension(2) :: Eigenvalue_Range
        real (kind = 8) :: EigenValue, E_aux, eigenvalue_tol, u0_tol, h, rp, delta
        logical :: Uniform

        ! Initial wave function guess
        u(n) = r(n)*exp(-r(n))
        u(n-1) = r(n-1)*exp(-r(n-1)) !Uniform guess, using Numerov
        du(n) = (1 - r(n))*exp(-r(n))*rp*delta*exp(n*delta) !Non-uniform guess, using Runge-Kutta
        u(1) = abs(2*u0_tol)

        ! Main loop for finding the Kohn-Sham eigenvalue
        do i=1, int_max

            EigenValue = (Eigenvalue_Range(1) + Eigenvalue_Range(2))/2

            if (Uniform) then
                ! Integrating wave function via Numerov
                u = Numerov(u, -2*(EigenValue-Potential), h, n)
            else
                ! Integrating wave function via Runge-Kutta
                u =  RungeKutta_KohnSham(u, du, -2*(EigenValue-Potential), rp**2*delta**2, delta, n)
            end if

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
        if (Uniform) then
            u = u/sqrt(sum(u*u*h))
        else
            do i=1, n-1
                Delta_r(i) = r(i+1)-r(i)
            end do
            Delta_r(n) = r(n)-r(n-1)
            u = u/sqrt(sum(u*u*Delta_r))
        end if

    end subroutine KohnSham1D

    subroutine Poisson(r, u, Potential_U, n, h, rp, delta, Uniform)
    ! Routine to solve a Poisson equation of radial potential Hartree potential.
    ! Parameters
    !   r: 1D array of coordinates to run by;
    !   u: Radial wave function u(r) = r*\psi(r);
    !   h: Discretization of array r;
    ! Returns
    !   Potential_U: U radial function, related with Hartree Potential by Hartree = r*U(r).

        integer :: n, i
        real (kind = 8), dimension(n) :: r, u, Potential_U, dPotential_U
        real (kind = 8) :: h, rp, delta, a
        logical :: Uniform

        ! Boundary condition at r=0 and first guess
        Potential_U(1) = r(1)

        if (Uniform) then
            Potential_U(2) = r(2)
            Potential_U = Verlet(Potential_U, -u**2/r, h, n) ! Integrating Potential U via Verlet
        else
            dPotential_U(1) = (2*r(1)+1)*exp(-2*r(1))*rp*delta*exp(n*delta)
            Potential_U = RungeKutta_Poisson(Potential_U, dPotential_U,&
            &-u**2, rp*delta**2, delta, n) ! Integrating Potential U via Rung-Kutta
        end if

        ! Satisfaying boundary condition at r_max by adding term a*r
        a = (1 - Potential_U(size(Potential_U))) / r(size(r))
        Potential_U = Potential_U + a*r

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
            u = Numerov(u, -2*(Eigenvalues(i)-Potential), h, n)

            ! Normalizing u
            u = u/sqrt(sum(u*u*h))
            u0s(i) = u(1)

            !print *,"Eigen Value tested:", Eigenvalues(i)
            !print *,"Boundary term:", u(1)
            !print *,"*-----------------------*"

            i = i + 1

        end do

    end subroutine KohnShamSweep

    function Numerov(u, f, h, n)
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
        real (kind = 8), dimension(n) :: u, f, Numerov

        do i=2, n-1
            a = 2*(1 + 5*(h**2)*f(n-i+1)/12)
            b = (1 - (h**2)*f(n-i+2)/12)
            c = (1 - (h**2)*f(n-i)/12)
            u(n-i) = (a*u(n-i+1) - b*u(n-i+2)) / c
        end do

        Numerov = u

    end function Numerov

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
        real (kind = 8), dimension(n) :: F, Potential_U, Verlet
        real (kind = 8) :: h

        do i=3, n
            Potential_U(i) = 2*Potential_U(i-1) - Potential_U(i-2) + (h**2)*F(i-1)
        end do

        Verlet = Potential_U

    end function Verlet

    function  RungeKutta_KohnSham(u, du, f, cte, delta, n)

        integer :: i, j, n
        real (kind = 8) :: k0, k1, k2, k3, l0, l1, l2, l3, cte, delta
        real (kind = 8), dimension(n) :: u, du, f, RungeKutta_KohnSham

        do j=1, n-1
            i = n-j+1

            k0 = -du(i)
            l0 = -(cte*exp(2*i*delta)*f(i)*u(i) + delta*du(i))

            k1 = -(du(i) + l0/2)
            l1 = -(cte*exp(2*(i-1/2)*delta)*((f(i)+f(i-1))/2)*(u(i)+k0/2) + delta*(du(i)+l0/2))

            k2 = -(du(i) + l1/2)
            l2 = -(cte*exp(2*(i-1/2)*delta)*((f(i)+f(i-1))/2)*(u(i)+k1/2) + delta*(du(i)+l1/2))

            k3 = -(du(i) + l2)
            l3 = -(cte*exp(2*(i-1)*delta)*f(i-1)*(u(i)+k2) + delta*(du(i)+l2))

            u(i-1) = u(i) + (k0+2*k1+2*k2+k3)/6
            du(i-1) = du(i) + (l0+2*l1+2*l2+l3)/6
        end do

        RungeKutta_KohnSham = u

    end function RungeKutta_KohnSham

    function  RungeKutta_Poisson(u, du, f, cte, delta, n)

        integer :: i, j, n
        real (kind = 8) :: k0, k1, k2, k3, l0, l1, l2, l3, cte, delta
        real (kind = 8), dimension(n) :: u, du, f, RungeKutta_Poisson

        do i=1, n-1

            k0 = du(i)
            l0 = cte*f(i)/(exp(-i*delta)-exp(-2*i*delta)) + delta*du(i)

            k1 = du(i) + l0/2
            l1 = cte*((f(i)+f(i+1))/2)/(exp(-(i+1/2)*delta)-exp(-2*(i+1/2)*delta)) + delta*(du(i)+l0/2)

            k2 = du(i) + l1/2
            l2 = cte*((f(i)+f(i+1))/2)/(exp(-(i+1/2)*delta)-exp(-2*(i+1/2)*delta)) + delta*(du(i)+l1/2)

            k3 = du(i) + l2
            l3 = cte*f(i+1)/(exp(-(i+1)*delta)-exp(-2*(i+1)*delta)) + delta*(du(i)+l2)

            u(i+1) = u(i) + (k0+2*k1+2*k2+k3)/6
            du(i+1) = du(i) + (l0+2*l1+2*l2+l3)/6
        end do

        RungeKutta_Poisson = u

    end function RungeKutta_Poisson

end module util

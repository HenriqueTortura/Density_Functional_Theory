subroutine HydrogenAtom(r_range, Eigenvalue_range, KS_int_max,&
&Eigenvalue_tol, u0_tol, Uniform_Numerov, h, j_max, delta, write_data, u0)

    integer, intent(in) :: j_max, KS_int_max
    integer :: n,  i, AllocateStatus
    real (kind = 8) , parameter :: pi = 3.141592653589793

    real (kind = 8), intent(in), dimension(2) :: r_range, Eigenvalue_Range
    real (kind = 8), intent(in) ::  Eigenvalue_tol, u0_tol, h, delta
    real (kind = 8) :: Eigenvalue, rp

    real (kind = 8), intent(out) :: u0

    real (kind = 8), dimension(:), allocatable :: r, u, Potential, Potential_U

    logical, intent(in), dimension(2) :: Uniform_Numerov
    logical, intent(in):: write_data

    if (Uniform_Numerov(1)) then
        n = int((r_range(2)-r_range(1))/h)
    else
        rp = r_range(2)/(exp(j_max*delta)-1)
        n = j_max
    end if

    allocate(r(n), u(n), Potential(n), Potential_U(n), stat = AllocateStatus)
    if (AllocateStatus /= 0) stop "*** Not enough memory ***"

    ! Initializing radial coordinates and effective potential
    do i=1, n

        if (Uniform_Numerov(1)) then
            r(i) = r_range(1) + h*i
        else
            r(i) = rp*(exp(i*delta)-1)
        end if

        Potential(i) = -1/r(i)
    end do

    ! Solve Kohn-Sham
    call KohnSham1D(r, u, Potential, Eigenvalue_Range, Eigenvalue,&
    &KS_int_max, Eigenvalue_tol, u0_tol, n, h, rp, delta, Uniform_Numerov)

    u0 = u(1)

    call Poisson(r, u, Potential_U, n, h, rp, delta, Uniform_Numerov)

    ! To plot u(r) and Potential U
    if (write_data) then
        open(1, file='Hydrogen_u.dat', status='replace')
            do i=1, n
                write(1,*) r(i), u(i)
            end do
        close(1)
        open(2, file='Hydrogen_Potential_U.dat', status='replace')
            do i=1, n
                write(2,*) r(i), Potential_U(i)
            end do
        close(2)
    end if

end subroutine HydrogenAtom

subroutine HeliumAtom(r_range, Eigenvalue_range, SelfCons_int_max, KS_int_max,&
&Eigenvalue_tol, u0_tol, SelfCons_tol, Uniform_Numerov, h, j_max, delta)

    integer :: n, j_max, i, KS_int_max, SelfCons_int_max, AllocateStatus
    real (kind = 8) , parameter :: pi = 3.141592653589793

    real (kind = 8), dimension(2) :: r_range, Eigenvalue_Range, E_Range_aux
    real (kind = 8)  :: Eigenvalue, Eigenvalue_aux, Energy,&
    &Eigenvalue_tol, u0_tol, SelfCons_tol, h, delta, rp

    real (kind = 8), dimension(:), allocatable :: r, u, Ext_Potential, Hartree, Exchange, Potential_U, j_array

    logical, dimension(2) :: Uniform_Numerov

    if (Uniform_Numerov(1)) then
        n = int((r_range(2)-r_range(1))/h)
    else
        rp = r_range(2)/(exp(j_max*delta)-1)
        n = j_max
    end if

    allocate(r(n), u(n), Ext_Potential(n), Hartree(n), Exchange(n), Potential_U(n), j_array(N), stat = AllocateStatus)
    if (AllocateStatus /= 0) stop "*** Not enough memory ***"

    ! Initializing radial coordinates and potentials
    do i=1, n
        j_array(i) = i

        if (Uniform_Numerov(1)) then
            r(i) = r_range(1) + h*i
        else
            r(i) = rp*(exp(i*delta)-1)
        end if

        Ext_Potential(i) = -2/r(i)
        Hartree(i) = 0
        Exchange(i) = 0
    end do

    Eigenvalue_aux = 0
    Eigenvalue = Eigenvalue_aux + 2*SelfCons_tol

    do i=1, SelfCons_int_max

        print *,'**************************'
        print *,'**************************'
        print *,"Iteration: ",i

        E_Range_aux = Eigenvalue_Range
        Eigenvalue_aux = Eigenvalue

        call KohnSham1D(r, u, Ext_Potential+Hartree+Exchange, E_Range_aux, Eigenvalue,&
        &KS_int_max, Eigenvalue_tol, u0_tol, n, h, rp, delta, Uniform_Numerov)

        if (abs(Eigenvalue-Eigenvalue_aux)>=SelfCons_tol) then
            call Poisson(r, u, Potential_U, n, h, rp, delta, Uniform_Numerov)

            Hartree = 2 * Potential_U / r

            Exchange = -((3./2.)*(u/(pi * r))**2.)**(1.0/3.0)
        else
            exit
        end if

    end do

    if (Uniform_Numerov(1)) then
        Energy = 2. * Eigenvalue - sum(Hartree*(u**2.)*h) - (1./2.)*sum(Exchange*(u**2.)*h)
    else
        Energy = 2. * Eigenvalue - rp*delta*sum(Hartree*(u**2.)*exp(j_array*delta))&
        & - rp*delta*(1./2.)*sum(Exchange*(u**2.)*exp(j_array*delta))
    end if

    print *,"Eigenvalue absolute true error: ", abs(-0.52-Eigenvalue)
    print *,"**********"
    print *,"Energy: ",Energy
    print *,"Energy absolute true error: ", abs(-2.72-Energy)
    print *,"**********"

end subroutine HeliumAtom

subroutine KohnSham1D(r, u, Potential, Eigenvalue_Range, EigenValue,&
&KS_int_max, eigenvalue_tol, u0_tol, n, h, rp, delta, Uniform_Numerov)
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

    integer :: i, j, KS_int_max, n !Dimension of Potential, u and r
    real (kind = 8), dimension(n) :: Potential, u, r, du, Delta_r, v, f
    real (kind =8), dimension(2) :: Eigenvalue_Range
    real (kind = 8) :: EigenValue, eigenvalue_tol, u0_tol, h, rp, delta
    logical, dimension(2) :: Uniform_Numerov

    ! Initial wave function guess
    !Uniform guess, using Numerov
    u(n) = r(n)*exp(-r(n))
    u(n-1) = r(n-1)*exp(-r(n-1))

    !Non-uniform guess, using Runge-Kutta
    du(n) = (1 - r(n))*exp(-r(n))*rp*delta*exp(n*delta)

    !Non-uniform guess, using Numerov
    v(n) = u(n)*exp(-n*delta/2)
    v(n-1) = u(n-1)*exp(-(n-1)*delta/2)

    ! Main loop for finding the Kohn-Sham eigenvalue
    do i=1, KS_int_max

        EigenValue = (Eigenvalue_Range(1) + Eigenvalue_Range(2))/2

        if (Uniform_Numerov(1)) then ! Integrating wave function via Numerov
            u = Numerov(u, -2*(EigenValue-Potential), h, n)

        else if (Uniform_Numerov(2)) then ! Non-uniform, integrating wave function via Numerov

            do j=1, n
                f(j) = (delta**2)/4 - 2*(rp**2)*(delta**2)*exp(2*j*delta)*(EigenValue-Potential(j))
            end do

            v = Numerov(v, f, 1._8, n)

            do j=2, n-1
                u(n-j) = v(n-j)*exp((n-j)*delta/2)
            end do

        else ! Non-uniform, integrating wave function via Runge-Kutta
            u =  RungeKutta_KohnSham(u, du, -2*(EigenValue-Potential), rp**2*delta**2, delta, n)

        end if


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

    ! Normalizing u
    if (Uniform_Numerov(1)) then
        u = u/sqrt(sum(u*u*h))
    else
        do i=1, n-1
            Delta_r(i) = r(i+1)-r(i)
        end do
        Delta_r(n) = r(n)-r(n-1)
        u = u/sqrt(sum(u*u*Delta_r))
    end if

contains

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

end subroutine KohnSham1D

subroutine Poisson(r, u, Potential_U, n, h, rp, delta, Uniform_Numerov)
! Routine to solve a Poisson equation of radial potential Hartree potential.
! Parameters
!   r: 1D array of coordinates to run by;
!   u: Radial wave function u(r) = r*\psi(r);
!   h: Discretization of array r;
! Returns
!   Potential_U: U radial function, related with Hartree Potential by Hartree = r*U(r).

    integer :: n
    real (kind = 8), dimension(n) :: r, u, Potential_U, dPotential_U
    real (kind = 8) :: h, rp, delta, a
    logical, dimension(2) :: Uniform_Numerov

    ! Boundary condition at r=0 and first guess
    Potential_U(1) = r(1)

    if (Uniform_Numerov(1)) then
        Potential_U(2) = r(2)
        Potential_U = Verlet(Potential_U, -u**2/r, h, n) ! Integrating Potential U via Verlet
    else
        dPotential_U(1) = (2*r(1)+1)*exp(-2*r(1))*rp*delta*exp(n*delta)
        Potential_U = RungeKutta_Poisson(Potential_U, dPotential_U,&
        &-u**2, rp*delta**2, delta, n) ! Integrating Potential U via Rung-Kutta
    end if

    ! Satisfaying boundary condition at r_max by adding term a*r
    a = (1 - Potential_U(n)) / r(n)
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
        real (kind = 8), dimension(n) :: F, Potential_U, Verlet
        real (kind = 8) :: h

        do i=3, n
            Potential_U(i) = 2*Potential_U(i-1) - Potential_U(i-2) + (h**2)*F(i-1)
        end do

        Verlet = Potential_U

    end function Verlet

    function  RungeKutta_Poisson(u, du, f, cte, delta, n)

        integer :: i, n
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

end subroutine Poisson

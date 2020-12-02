module dft

    implicit none

contains

    subroutine HydrogenAtom(r_range, Eigenvalue_range, KS_int_max,&
    &Eigenvalue_tol, u0_tol, Uniform_Numerov, h, j_max, delta, verbose,&
    &write_data, path, Eigenvalue, u0, Hartree_Energy, Exchange_Energy,&
    &Correlation_Energy)

        integer, intent(in) :: j_max, KS_int_max
        integer :: n,  i, AllocateStatus
        real (kind = 8) , parameter :: pi = 3.141592653589793

        real (kind = 8), intent(in), dimension(2) :: r_range, Eigenvalue_Range
        real (kind = 8), intent(in) ::  Eigenvalue_tol, u0_tol, h, delta
        real (kind = 8) :: rp
        real (kind = 8), intent(out) :: Hartree_Energy, Exchange_Energy, Correlation_Energy

        real (kind = 8), intent(out) :: Eigenvalue, u0

        real (kind = 8), dimension(:), allocatable :: r, u, Potential_U, Potential,&
        &j_array, Hartree, Exchange, r_s, e_c

        logical, intent(in), dimension(2) :: Uniform_Numerov
        logical, intent(in):: verbose, write_data
        character(len=256), intent(in) :: path
        character(len=:), allocatable :: realpath

        realpath = trim(path)

        if (Uniform_Numerov(1)) then
            n = int((r_range(2)-r_range(1))/h)
        else
            rp = r_range(2)/(exp(j_max*delta)-1)
            n = j_max
        end if

        allocate(r(n), u(n), Potential(n), Potential_U(n), j_array(n), Hartree(n),&
        &Exchange(n), r_s(n), e_c(n), stat = AllocateStatus)
        if (AllocateStatus /= 0) stop '*** Not enough memory ***'

        ! Initializing radial coordinates and effective potential
        do i=1, n
            j_array(i) = i

            if (Uniform_Numerov(1)) then
                r(i) = r_range(1) + h*i
            else
                r(i) = rp*(exp(i*delta)-1)
            end if

            Potential(i) = -1/r(i)
        end do

        ! Solve Kohn-Sham
        call KohnSham1D(r, u, Potential, Eigenvalue_Range, Eigenvalue,&
        &KS_int_max, Eigenvalue_tol, u0_tol, n, h, rp, delta, Uniform_Numerov, verbose)

        u0 = u(1)

        call Poisson(r, u, Potential_U, n, h, rp, delta, Uniform_Numerov)

        Hartree = Potential_U / r
        Exchange = -((3./4.)*(u/(pi * r))**2.)**(1.0/3.0)
        r_s = ((3.*r**2.)/u**2.)**(1./3.)
        e_c = -(0.045/2)*( (1+(r_s/21)**3)*log(1+(21/r_s)) + r_s/42 - (r_s/21)**2. - 1/3 )

        if (Uniform_Numerov(1)) then
            Hartree_Energy = (1./2.)*sum(Hartree*(u**2.)*h)
            Exchange_Energy = -( (9./(16.*pi))**(2./3.) )*sum(( ((u**4.)/r)**(2./3.) )*h)
            Correlation_Energy = sum((u**2.)*e_c*h)
        else
            Hartree_Energy = rp*delta*(1./2.)*sum(Hartree*(u**2.)*exp(j_array*delta))
            Exchange_Energy = -rp*delta*( (9./(16.*pi))**(2./3.) )*sum(( ((u**4.)/r)**(2./3.) )*exp(j_array*delta))
            Correlation_Energy = rp*delta*sum((u**2.)*e_c*exp(j_array*delta))
        end if

        if (verbose) then
            print *,'**********'
            print *,'Hartree energy (eV): ', 27.2113845*Hartree_Energy
            print *,'Hartree energy: ', Hartree_Energy

            print *,'**********'
            print *,'Exchange energy: ', Exchange_Energy
            print *,'Ratio: ', Exchange_Energy/Hartree_Energy
            print *,'**********'
            print *,'Correlation energy (Hedin-Lundqvist): ', Correlation_Energy
            print *,'Ratio: ', Correlation_Energy/Hartree_Energy
            print *,'**********'
            print *,'Combined ratio: ', (Exchange_Energy+Correlation_Energy)/Hartree_Energy
            print *,'**********'
        end if

        if (write_data) then
            realpath = trim(path)
            print *,'Writing data'
            open(1, file=realpath//'Hydrogen_data.dat', status='replace')
                do i=1, size(r)
                    write(1,*) r(i), u(i), Potential_U(i)
                end do
            close(1)
        end if

    end subroutine HydrogenAtom

    subroutine HeliumAtom(r_range, Eigenvalue_range, SelfCons_int_max, KS_int_max,&
    &Eigenvalue_tol, u0_tol, SelfCons_tol, Uniform_Numerov, h, j_max, delta, verbose,&
    &SelfInteractionCorrection, Exchange_Method, Correlation_Method, Z, N_electrons,&
    &Energy, Eigenvalue, Hartree_Correction, Exchange_Correction, Correlation_Correction)

        integer, intent(in) :: SelfCons_int_max, KS_int_max, j_max,&
        &Z, N_electrons, Exchange_Method, Correlation_Method
        integer :: n, i, j, AllocateStatus

        integer, dimension(:), allocatable :: indices1, indices2
        real (kind = 8) , parameter :: pi = 3.141592653589793

        real (kind = 8), dimension(2), intent(in) :: r_range, Eigenvalue_Range
        real (kind = 8), dimension(2) :: E_Range_aux

        real (kind = 8), intent(in) :: Eigenvalue_tol, u0_tol, SelfCons_tol,&
        &h, delta
        real (kind = 8), intent(out) :: Energy, Eigenvalue, Hartree_Correction,&
        &Exchange_Correction, Correlation_Correction
        real (kind = 8)  :: Eigenvalue_aux, rp, ExHartree_ratio

        real (kind = 8), dimension(:), allocatable :: r, u, Ext_Potential,&
        &Hartree, Exchange, Correlation, Potential_U, j_array, r_s, e_c

        logical, dimension(2), intent(in) :: Uniform_Numerov
        logical, intent(in) :: verbose, SelfInteractionCorrection

        if (Uniform_Numerov(1)) then
            n = int((r_range(2)-r_range(1))/h)
        else
            rp = r_range(2)/(exp(j_max*delta)-1)
            n = j_max
        end if

        allocate(r(n), u(n), Ext_Potential(n), Hartree(n), Exchange(n),&
        &Correlation(n), Potential_U(n), j_array(n), r_s(n), e_c(n), stat = AllocateStatus)
        if (AllocateStatus /= 0) stop '*** Not enough memory ***'

        ! Initializing radial coordinates and potentials
        do i=1, n
            j_array(i) = i

            if (Uniform_Numerov(1)) then
                r(i) = r_range(1) + h*i
            else
                r(i) = rp*(exp(i*delta)-1)
            end if

            Ext_Potential(i) = -Z/r(i)
            Hartree(i) = 0
            Exchange(i) = 0
            Correlation(i) = 0
        end do

        Eigenvalue_aux = 0
        Eigenvalue = Eigenvalue_aux + 2*SelfCons_tol

        do i=1, SelfCons_int_max

            if (verbose) then
                print *,'**************************'
                print *,'**************************'
                print *,'Iteration: ',i
            end if

            E_Range_aux = Eigenvalue_Range
            Eigenvalue_aux = Eigenvalue

            call KohnSham1D(r, u, Ext_Potential+Hartree+Exchange+Correlation, E_Range_aux, Eigenvalue,&
            &KS_int_max, Eigenvalue_tol, u0_tol, n, h, rp, delta, Uniform_Numerov, verbose)

            if (abs(Eigenvalue-Eigenvalue_aux)>=SelfCons_tol) then
                call Poisson(r, u, Potential_U, n, h, rp, delta, Uniform_Numerov)

                if (SelfInteractionCorrection) then
                    Hartree = (N_electrons/2)*Potential_U / r
                else
                    Hartree = N_electrons*Potential_U / r
                end if

                if (Exchange_Method==0) then
                    Exchange = 0
                else if (Exchange_Method==1) then
                    Exchange = -((3.*N_electrons/4.)*(u/(pi * r))**2.)**(1.0/3.0)
                end if

                if (Correlation_Method == 0) then !Only exchange correction
                    Correlation = 0

                else if (Correlation_Method == 1) then !Hedin-Lundqvist
                    r_s = ((3.*r**2.)/u**2.)**(1./3.)
                    e_c = -(0.045/2)*( (1+(r_s/21)**3)*log(1+(21/r_s)) + r_s/42 - (r_s/21)**2. - 1/3 )
                    Correlation = -(0.045/2)*log(1+(21/r_s))

                else if (Correlation_Method == 2) then !Perdew-Zunger
                    r_s = ((3.*r**2.)/u**2.)**(1./3.)
                    indices1 = pack([(j, j=1,size(r_s))], r_s .GE. 1)
                    indices2 = pack([(j, j=1,size(r_s))], r_s .LT. 1)

                    e_c(indices1) = -0.1423/( 1. + 1.0529*sqrt(r_s(indices1)) + 0.3334*r_s(indices1) )

                    e_c(indices2) = 0.0311*log(r_s(indices2)) -0.048&
                    &+0.002*r_s(indices2)*log(r_s(indices2)) -0.0116*r_s(indices2)

                    Correlation(indices1) = e_c(indices1)*&
                    &((1. + (7./6.)*1.0529*sqrt(r_s(indices1)) + 0.3334*r_s(indices1))&
                    & / (1. + 1.0529*sqrt(r_s(indices1)) + 0.3334*r_s(indices1)))

                    Correlation(indices2) = 0.0311*log(r_s(indices2)) -0.048 -(0.0311/3.)&
                    &+(2./3.)*0.002*r_s(indices2)*log(r_s(indices2))&
                    &+(2.*(-0.0116)-0.002)*r_s(indices2)/3
                end if
            else
                exit
            end if

        end do

        if (Uniform_Numerov(1)) then
            Hartree_Correction = -(N_electrons/2.)*sum(Hartree*(u**2.)*h)
            Exchange_Correction = -(N_electrons/4.)*sum(Exchange*(u**2.)*h)
            Correlation_Correction = -sum(Correlation*(u**2.)*h) + sum((u**2.)*e_c*h)
        else
            Hartree_Correction = -rp*delta*(N_electrons/2.)*sum(Hartree*(u**2.)*exp(j_array*delta))
            Exchange_Correction = -rp*delta*(N_electrons/4.)*sum(Exchange*(u**2.)*exp(j_array*delta))
            Correlation_Correction = -rp*delta*sum(Correlation*(u**2.)*exp(j_array*delta))&
            &+ rp*delta*sum((u**2.)*e_c*exp(j_array*delta))
        end if

        Energy = N_electrons*Eigenvalue + Hartree_Correction + Exchange_Correction + Correlation_Correction

        if (verbose) then
            !print *,'Eigenvalue absolute true error: ', abs(-0.52-Eigenvalue)
            print *,'**********'
            print *,'Energy: ', Energy
            !print *,'Energy absolute true error: ', abs(-2.72-Energy)
            print *,'**********'
            print *,'Correction terms'
            print *, 'Hartree: ', Hartree_Correction
            print *, 'Exchange: ', Exchange_Correction
            print *, 'Correlation: ', Correlation_Correction
            print *,'**********'
        end if

    end subroutine HeliumAtom

    subroutine KohnSham1D(r, u, Potential, Eigenvalue_Range, EigenValue,&
    &KS_int_max, eigenvalue_tol, u0_tol, n, h, rp, delta, Uniform_Numerov, verbose)
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
        logical :: verbose

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

            if (verbose) then
                print *,'Eigen Value tested:', EigenValue
                print *,'Boundary term:', u(1)
                print *,'*-----------------------*'
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

        ! Satisfying tolerance for final result
        if (verbose) then
            if ((abs(u(1))<u0_tol .and. abs(EigenValue-Eigenvalue_Range(1))<eigenvalue_tol)) then
                print *,'Converged for eigenvalue: ', EigenValue
            else
                print *,'Did not conveged within ', (i-1), ' iterations'
            end if
        end if


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

    end subroutine Poisson

    subroutine KohnShamSweep(r_range, Eigenvalue_Range, h, n_eigenvalues, path)

        real (kind = 8), dimension(2), intent(in) :: r_range, Eigenvalue_Range
        real (kind = 8), intent(in) :: h
        integer, intent(in) :: n_eigenvalues
        character(len=256), intent(in) :: path

        real (kind = 8), dimension(:), allocatable :: Potential, u, r
        real (kind = 8) :: Eigenvalue
        integer :: i, n, AllocateStatus
        character(len=:), allocatable :: realpath

        realpath = trim(path)
        n = int((r_range(2)-r_range(1))/h)

        allocate(r(n), u(n), Potential(n), stat = AllocateStatus)
        if (AllocateStatus /= 0) stop '*** Not enough memory ***'

        ! Initializing radial coordinates and effective potential
        do i=1, n
            r(i) = r_range(1) + h*i
            Potential(i) = -1/r(i)
        end do

        ! Initial wave function guess
        u(n) = r(n)*exp(-r(n))
        u(n-1) = r(n-1)*exp(-r(n-1))

        open(1, file=realpath//'Hydrogen_u0s.dat', status='replace')
            do i=1, n_eigenvalues ! Main loop for finding the Kohn-Sham eigenvalue

                !print *,i,'/',n_eigenvalues
                ! Integrating wave function via Numerov
                Eigenvalue = Eigenvalue_Range(1) + (i-1)*(Eigenvalue_Range(2)-Eigenvalue_Range(1))/(n_eigenvalues-1)
                u = Numerov(u, -2*(Eigenvalue-Potential), h, n)

                ! Normalizing u
                u = u/sqrt(sum(u*u*h))

                !print *,'Eigen Value tested:', Eigenvalue
                !print *,'Boundary term:', u(1)
                !print *,'*-----------------------*'
                write(1,*) Eigenvalue, u(1)

            end do
        close(1)

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

    function RungeKutta_KohnSham(u, du, f, cte, delta, n)

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

    function RungeKutta_Poisson(u, du, f, cte, delta, n)

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

end module dft

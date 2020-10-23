program HeliumAtom

    use util

    implicit none

    ! Setting main parameters
    integer, parameter :: r0 = 0, r_max = 30, int_max_KS = 100, int_max = 100 ! Radial and iterations limits
    real (kind = 8) , parameter :: pi = 3.141592653589793, h = 0.0001, & ! Discretization
    &Eigenvalue_tol = 0.001, u0_tol= 0.001, self_tol=0.001 ! Tolerances

    real (kind = 8), dimension(2) :: E_Range_aux, Eigenvalue_Range = (/ -5., 0./) ! Eigenvalue range

    integer, parameter :: n = int((r_max-r0)/h) ! Number of steps

    integer i, j
    real (kind = 8) , dimension(1:n) :: r, u, Ext_Potential, Hartree, Exchange, Potential_U
    real (kind = 8)  :: Eigenvalue, Eigenvalue_aux, Energy

    real :: start, finish

    call cpu_time(start)

    ! Initializing radial coordinates and potentials
    do i=1, n
        r(i) = r0 + h*i
        Ext_Potential(i) = -2/r(i)
        Hartree(i) = 0
        Exchange(i) = 0
    end do

    Eigenvalue_aux = 0
    Eigenvalue = Eigenvalue_aux + 2*self_tol

    i = 1

    do i=1, int_max

        print *,'**************************'
        print *,'**************************'
        print *,"Iteration: ",i

        E_Range_aux = Eigenvalue_Range
        Eigenvalue_aux = Eigenvalue

        call KohnSham1D(r, u, Ext_Potential+Hartree+Exchange, E_Range_aux, Eigenvalue,&
        &int_max_KS, Eigenvalue_tol, u0_tol, n, h, 0d8, 0d8, .TRUE.)

        if (abs(Eigenvalue-Eigenvalue_aux)>=self_tol) then
            call Poisson(r, u, Potential_U, n, h, 0d8, 0d8, .TRUE.)

            Hartree = 2 * Potential_U / r

            Exchange = -((3./2.)*(u/(pi * r))**2.)**(1.0/3.0)
        else
            exit
        end if

    end do

    Energy = 2. * Eigenvalue - sum(Hartree*(u**2.)*h) - (1./2.)*sum(Exchange*(u**2.)*h)
    print *,"Energy: ",Energy

    call cpu_time(finish)
    print *,"CPU time = ",finish-start,"s"

end program

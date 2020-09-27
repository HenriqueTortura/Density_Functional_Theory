program HeliumAtom

    use util

    implicit none

    ! Setting main parameters
    integer, parameter :: r0 = 0, r_max = 30, int_max_KS = 100, int_max = 100
    real (kind = 8) , parameter :: pi = 3.141592653589793, h = 0.0001, Eigenvalue_tol = 0.001, u0_tol= 0.001, self_tol=0.001 ! Discretization

    real (kind = 8), dimension(1:2) :: Eigenvalue_Range

    integer, parameter :: n = int((r_max-r0)/h) ! Number of steps

    integer i, j
    real (kind = 8) , dimension(1:n) :: r, u, Ext_Potential, Hartree, Exchange, Potential_U
    real (kind = 8)  :: Eigenvalue_min, Eigenvalue_max, Eigenvalue, Eigenvalue_aux, Energy

    integer, parameter :: n_Eigenvalues = 101
    real (kind = 8), dimension(1:n_Eigenvalues) :: Eigenvalues, u0s


    ! Initializing radial coordinates and potentials
    do i=1, n
        r(i) = r0 + h*i
        Ext_Potential(i) = -2/r(i)
        Hartree(i) = 0
        Exchange(i) = 0
    end do

    Eigenvalue = 42
    Eigenvalue_aux = 0

    i = 1

    do while(abs(Eigenvalue-Eigenvalue_aux)>=self_tol .and. i<=int_max)
        Eigenvalue_Range = (/ -10., 0. /)
        Eigenvalue_aux = Eigenvalue

        print *,'**************************'
        print *,'**************************'
        print *,"Iteration: ",i

        call KohnSham1D(r, u, Ext_Potential+Hartree+Exchange, Eigenvalue_Range, Eigenvalue, h, int_max_KS, Eigenvalue_tol, u0_tol)

        call Poisson(r, u, Potential_U, h)

        Hartree = 2 * Potential_U / r

        Exchange = -((3./2.)*(u/(pi * r))**2.)**(1.0/3.0)

        i = i + 1

    end do

    Energy = 2. * Eigenvalue - sum(Hartree*(u**2.)*h) - (1./2.)*sum(Exchange*(u**2.)*h)
    print *,"Energy: ",Energy

end program

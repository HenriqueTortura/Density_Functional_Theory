program Sweep

    use dft

    implicit none

    ! Setting main parameters
    real (kind = 8) , parameter :: h = 0.0001

    real (kind = 8), dimension(2) :: r_range = (/ 0., 50./), Eigenvalue_Range = (/-1., 5./)

    integer, parameter :: n_eigenvalues = 501

    character(len=256) :: path = 'data/'

    !---------------------------------------------------------------------------------------!

    call KohnShamSweep(r_range, Eigenvalue_Range, h, n_eigenvalues, path)

end program


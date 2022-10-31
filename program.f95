
program finalproject

!Make actual variable and parameters!

    Lx = 6e+06 !domain size in x direction
    Ly = 2e+06 !domain size in y direction

    d = 5e+05       !Resolution (changes for each run)
    !d = 2.5e+05
    !d = 1.25e+05

    Nx = Lx/d + 1 !number of grid points in the x direction (13, 25, 49)
    Ny = Ly/d + 1 !number of grid points in the y direction (5, 9, 17)

    hs_t = 2e+03 !height of the topography

    delt = 10*60 !time step (10 min) in seconds (changes for each run)
    !delt = 5*60
    !delt = 2.5*60

    an = (55*delt)/24   !alphas for 3rd order adams-bashforth schemes
    an_1 = (59*delt)/24
    an_2 = (37*delt)/24
    an_3 = (9*delt)/24


!Brian and Lauren : topography variable, time discretization by Monday (11/07)
!Michael and Taylor: begin momentum equations after Brain and Lauren finish







end program finalproject

program finalproject

!Make actual variable and parameters!

    Lx = 6e+06 !domain size in x direction, real numbers
    Ly = 2e+06 !domain size in y direction, real numbers

    d = 5e+05       !Resolution (changes for each run)
    !d = 2.5e+05
    !d = 1.25e+05

    Nx = Lx/d + 1 !number of grid points in the x direction (13, 25, 49), integers
    Ny = Ly/d + 1 !number of grid points in the y direction (5, 9, 17), integers

    hs_t = 2e+03 !height of the topography
    ! need to initialize h0, cannot equal/exceed top of the model or 5000! 

    delt = 10*60 !time step (10 min) in seconds (changes for each run)
    !delt = 5*60
    !delt = 2.5*60

    an = (55*delt)/24   !alphas for 3rd order adams-bashforth schemes, double check these 
    an_1 = (59*delt)/24
    an_2 = (37*delt)/24
    an_3 = (9*delt)/24


!Brian and Lauren : topography variable, time discretization by Monday (11/07)

    if (d == 5e+05) then       !topography variable for first resolution
        hs(Nx/2) = hs_t
    end if

    if (d == 2.5e+05) then.    !topography variables for second resolution
        hs(Nx/2-1) = 1e+03
        hs(Nx/2) = hs_t
        hs(Nx/2+1) = 1e+03
    end if

    if (d == 1.25e+05) then.   !topography variable for third resolution
        hs(Nx/2-2) = 0.5e+03
        hs(Nx/2-1) = 1.5e+03
        hs(Nx/2) = hs_t
        hs(Nx/2+1) = 0.5e+03
        hs(Nx/2+2) = 0.5e+03
    end if

!Time discretization using 3rd order adams-bashforth scheme 

!variables 

hu0(Nx,Ny,3)
hv0(Nx,Ny,3)
hq0(Nx,Ny,3)
us0(Nx,Ny,3)
vs0(Nx,Ny,3)
h0 = 5000                     !top of fluid

h(Nx,Ny)
u(Nx,Ny)
v(Nx,Ny)
z(Nx,Ny)                      !add remaining new variables and initial conditions for q etc

!initial conditions, t=0 or n=1

u(1:Nx,1:Ny) = 0.5            !value of u
v(1:Nx,1) = 0                 !value of v
v(1:Nx,Ny) = 0                !value of v
V(1:Nx,2:Ny-1) = 0.1          !horizontal velocity 


do i = 1, Nx
    h(i,:) = h0-hs(i)   !we only have the height of the topography so we account for fluid above
end do 

hu0(1,:,1) = (h(Nx,:) + h(2,:))/2.0
hu0(Nx,:,1) = (h(Nx-1,:) + h(1,:))/2.0

do i = 2, Nx-1
    hu0(i,:,1) = (h(i-1,:) + h(i+1,:))/2.0
end do

hv0(:,1,1) = (h(:,Ny) +h(:,2))/2.0
hv0(:,Ny,1) = (h(:,Ny-1) + h(:,1))/2.0

do j = 2, Ny-1
    hv0(:,j,1) = (h(:,j-1) + h(:,j+1))/2.0
end do 

us0(:,:,1) = hu0(:,:,1)*u(:,:)
vs0(:,:,1) = hv0(:,:,1)*v(:,:) 

do n = 2,3

do i = 2, Nx-1
    do j = 2, Ny-1 
        h(i,j) = h(i,j)-delt*(us0(i+1,j+1,n-1)-us0(i,j+1,n-1)+vs0(i+1,j+1,n-1)-vs0(i+1,j,n-1))/d
    end do
end do

hu0(1,:,n) = (h(Nx,:) + h(2,:))/2.0
hu0(Nx,:,n) = (h(Nx-1,:) + h(1,:))/2.0

do i = 2, Nx-1
    hu0(i,:,n) = (h(i-1,:) + h(i+1,:))/2.0
end do 

hv0(:,1,n) = (h(:,Ny) +h(:,2))/2.0
hv0(:,Ny,n) = (h(:,Ny-1) + h(:,1))/2.0

do j = 2, Ny-1
    hv0(:,j,n) = (h(:,j-1) + h(:,j+1))/2.0
end do 

us0(:,:,n) = hu0(:,:,n)*u(:,:)
vs0(:,:,n) = hv0(:,:,n)*v(:,:)

end do

us1 = us0(:,:,1)
us2 = us0(:,:,2)
us3 = us0(:,:,3)
vs1 = vs0(:,:,1)
vs2 = vs0(:,:,2)
vs3 = vs0(:,:,3)
hu1 = hu0(:,:,1)
hu2 = hu0(:,:,2)
hu3 = hu0(:,:,3)
hv1 = hv0(:,:,1) 
hv2 = hv0(:,:,2)
hv3 = hv0(:,:,3)

nstep = 4

do n = 4, ntime
    nstep = nstep + 1
    do i = 2, Nx-1
        do j = 2, Ny-1
            h(i,j) = h(i,j)-an*(us1(i+1,j+1)-us1(i,j+1)+vs1(i+1,j+1)-vs1(i+1,j))+an_1*(us2(i+1,j+1)-us2(i,j+1)+vs2(i+1,j+1)-vs2(i+1,j))-an_2*(us3(i+1,j+1)-us3(i,j+1)+vs3(i+1,j+1)-vs3(i+1,j))
        end do
    end do

    do i = 2, Nx-1
        hu(i,:) = (h(i-1,:) + h(i+1,:))/2.0
    end do

    do j = 2, Ny-1
        hv(:,j) = (h(:,j-1) + h(:,j+1))/2.0
    end do

    us(:,:) = hu(:,:)*u(:,:)
    vs(:,:) = hv(:,:)*v(:,:)

    us1 = us2
    us2 = us3
    us3 = us
    vs1 = vs2
    vs2 = vs3
    vs3 = vs
    hu1 = hu2
    hu2 = hu3
    hu3 = hu
    hv1 = hv2 
    hv2 = hv3
    hv3 = hv

    if (nstep == 1440) then
        write(10,*) h, u, V
        nstep = 0
    end if 
end do !time loop 


!Michael and Taylor: begin momentum equations after Brian and Lauren finish







end program finalproject
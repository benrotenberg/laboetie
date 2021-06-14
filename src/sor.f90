! This subroutine is a Successive Over Relaxation (SOR) method implementation.
! see: Horbach and Frenkel, Phys. Rev. E 64, 061507 (2001)

subroutine sor

    USE precision_kinds, ONLY: dp, i2b
    USE system, ONLY: bjl, tot_sol_charge, pbc, phi, LaplacianOfPhi, kbt, c_plus, c_minus, supercell, node, time, t
    use module_convergence
    USE constants, ONLY: pi, x, y, z
    USE mod_lbmodel, ONLY: lbm
    USE myallocations
    use module_input, only: getinput, verbose

    IMPLICIT NONE

    integer, parameter :: maxiterations = 500000
    real(dp), parameter :: threshold = 1.0e-6 ! threshold on value of phi to be considered or not -> to be improved?
    real(dp), parameter :: omega = 1.4_dp     ! the over-ralaxation factor. 1.45 proposed by Horbach & Frenkel, PRE64 Eq.17
    real(dp) :: factor, h
    real(dp) :: dphi 
    integer :: iter ! number of iterations to achieve the tolerance
    integer :: i, j, k, l, imin, jmin, kmin, imax, jmax, kmax, Constant_Potential
    integer :: pmin, qmin, rmin, timeTMP, n1
    integer :: ipass,isw,jsw,ksw 
    integer :: isw2,jsw2,ksw2
    real(dp) :: phiold, PoissonDiff, LaplOfPhi
    real(dp), dimension(:,:,:), allocatable :: phitmp, phiR, phi_old
    real(dp), parameter :: zerodp = 0._dp

    
    Constant_Potential = getinput%int("Constant_Potential", defaultvalue=0) ! Ade : 1=true, 0=false
    if( tot_sol_charge==0.0_dp .and. Constant_Potential==0) then
        if(.not.allocated(phi)) call allocateReal3D(phi)
        phi = 0.0_dp
        return 
    end if
    
    call allocateReal3D( phitmp )
    phitmp = 0.0_dp
    if (.not. allocated(LaplacianOfPhi)) call allocateReal3D( LaplacianOfPhi )
    LaplacianOfPhi = 0.0_dp
    if (.not. allocated(phi_old)) call allocateReal3D(phi_old)
    phi_old = phi

    ! BR: corresponds to (4*pi*bjl)*(cs^2 /2) in Eq (17) of PRE64, 061507) ok for cs^2 = 1/2 for D3Q18 only
    factor = 4.0_dp*pi*bjl*kbt/2.0_dp

    ! Solving the Poisson equation Laplacian(beta*e*psi) = Laplacian(phi) = - (c_plus-c_minus) * e^2/(epsilon_0*epsilon_r)
    !     is equivalent to solving (2/(c^2*dt^2)*sum_i w_i*[phi(r+ci*dt)-phi(r)] = - [cplus(r)-cminus(r)] * 4*pi*bjl
    !     ie PoissonDiff = sum_i w_i[phi(r+ci*dt)-phi(r)] + factor*[cplus(r)-cminus(r)] = 0
    !     this is achieved iteratively with SOR
    !
    ! In the presence of metallic walls, we introduce slight modifications, see Asta et al. JCP 151, 114104 (2019)
    !     phi is not updated on nodes where the potential is imposed
    !     the contribution of interfacial links is modified to improve accuracy (location of the interface between nodes)
    !     the charge induced on the metallic nodes is computed from the Poisson equation (see LaplacianOfPhi)

    convergenceloop: do iter=1, maxiterations 

        ksw = 0
        ksw2 = 0
        imin = supercell%geometry%dimensions%indiceMin(x)
        imax = supercell%geometry%dimensions%indiceMax(x)
        jmin = supercell%geometry%dimensions%indiceMin(y)
        jmax = supercell%geometry%dimensions%indiceMax(y)
        kmin = supercell%geometry%dimensions%indiceMin(z)
        kmax = supercell%geometry%dimensions%indiceMax(z)
        do k = kmin, kmax
            do j = jmin, jmax
                do i = imin, imax
                    if( node(i,j,k)%isFixedPotential ) then
                        cycle ! don't use SOR on nodes where you impose the potential
                    else 
                        phitmp(i,j,k) = 0.0_dp
                        PoissonDiff = 0.0_dp
                        phiold = phi(i,j,k)                     ! Ade : used to compute the Laplacian of Phi
                        do l= lbm%lmin, lbm%lmax 
                            pmin = pbc(i-lbm%vel(l)%coo(x),x)  
                            qmin = pbc(j-lbm%vel(l)%coo(y),y)  
                            rmin = pbc(k-lbm%vel(l)%coo(z),z)
                            if( node(pmin,qmin,rmin)%isFixedPotential ) then
                                 PoissonDiff = PoissonDiff + lbm%vel(l)%a0 * (phi(pmin,qmin,rmin)-phiold)*2
                            else 
                                 PoissonDiff = PoissonDiff + lbm%vel(l)%a0 * (phi(pmin,qmin,rmin)-phiold)
                            endif 
                        end do
                        PoissonDiff = PoissonDiff + factor*(c_plus(i,j,k)-c_minus(i,j,k)) ! see PRE64, Horbach: Eq. 17
                        phitmp(i,j,k) = phiold + omega*PoissonDiff
                    endif
                end do
            end do
         end do
         where( .not. node%isFixedPotential ) phi = phitmp

        ! dphi is the difference between the electric potential between the beginning and end of the iteration.
        dphi = 0.0_dp 
 
        ! count the number of times the array is not zero
        n1 = count(abs(phi_old) > threshold)

        ! at each node, if phiold is different from 0, compute the normalized relative difference 
        ! between new and old potential
        if(n1/=0) dphi = sum( abs(  (phi - phi_old)/phi_old ), mask= abs(phi_old)>threshold) / real(n1,kind=dp) 
        phi_old = phi

        ! Have we reached convergence?
        if(iter>1 .and. dphi<target_error%target_error_sor ) then
            exit convergenceloop
        end if

        ! Report every 100 steps
        if( verbose ) then
          if( modulo(iter, 100)==0 ) then
            print*, iter,dphi,n1
          end if
        end if

    end do convergenceloop

    ! tell user if maximum convergence steps is reached, ie if no convergence is found
    if( iter >= maxiterations ) stop 'maximum iterations 500 000 reached without convergence in sor'

    ! compute charge corresponding to this distribution of potential
    !     this is used to compute the charge on the electrodes (when relevant)
    do k = kmin, kmax
        do j = jmin, jmax
            do i = imin, imax
                LaplOfPhi = 0.0_dp
                do l= lbm%lmin, lbm%lmax               
                    pmin = pbc(i-lbm%vel(l)%coo(x),x) 
                    qmin = pbc(j-lbm%vel(l)%coo(y),y) 
                    rmin = pbc(k-lbm%vel(l)%coo(z),z)
                    LaplOfPhi = LaplOfPhi + lbm%vel(l)%a0 * phi(pmin,qmin,rmin)   
                    if( .not. node(pmin,qmin,rmin)%isFixedPotential ) then
                        LaplOfPhi = LaplOfPhi + lbm%vel(l)%a0 * (phi(pmin,qmin,rmin)-phiold)*2
                    else 
                        LaplOfPhi = LaplOfPhi + lbm%vel(l)%a0 * (phi(pmin,qmin,rmin)-phiold)
                    endif 
                end do
                LaplacianOfPhi(i,j,k) = -LaplOfPhi/factor   ! this is in fact the charge on the node
            end do
        end do
    end do

    if( verbose ) then
        if(iter>1 .and. dphi<target_error%target_error_sor ) then
            print*,'SOR converged in',iter-1,'steps'
        end if
    end if

    !if( t==1 ) then 
    !    open(109, file='output/PHIsor.dat')
    !    do k = kmin, kmax
    !        write(109,*) k, SUM(phi(:,:,k))
    !    end do
    !    close(109)
    !endif


end subroutine sor

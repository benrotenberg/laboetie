! This subroutine is a Successive Over Relaxation (SOR) method implementation.
! It is now used for the computation of the phi electrostatic potential but it
! is a much more generic method, so it should be (TODO) an object called
! using a poisson solver.
! see: Horbach and Frenkel, Phys. Rev. E 64, 061507 (2001)

! Note de Max du 28 juin : je suis très intrigué par SOR qui fonctionne pour le potentiel constant.
! J'aimerais réécrire cette routine de façon plus générique avec des entrées et sorties explicites,
! pour pouvoir tester ce solveur dans MDFT, où finalement il pourrait très bien faire un job propre.

subroutine sor

    USE precision_kinds, ONLY: dp, i2b
    USE system, ONLY: bjl, tot_sol_charge, pbc, anormf0, phi, LaplacianOfPhi, kbt, c_plus, c_minus, supercell, node, time, t
    use module_convergence
    USE constants, ONLY: pi, x, y, z
    USE mod_lbmodel, ONLY: lbm
    USE myallocations
    use module_input, only: getinput, verbose

    IMPLICIT NONE

    integer, parameter :: maxiterations = 500000
    real(dp), parameter :: threshold = 1.0e-6 ! convergence tolerance
    real(dp), parameter :: omega = 1.4_dp ! the over-ralaxation factor. 1.45 proposed by Horbach & Frenkel, PRE64 Eq.17
    real(dp) :: factor, h
    real(dp) :: anorm ! what we want to minimize, ie diff between phi and phi-old
    real(dp) :: anormf, dphi
    integer :: iter ! number of iterations to achieve the tolerance
    integer :: i, j, k, l, imin, jmin, kmin, imax, jmax, kmax, Constant_Potential
    integer :: pmin, qmin, rmin, timeTMP, n1
    integer :: ipass,isw,jsw,ksw 
    integer :: isw2,jsw2,ksw2
    real(dp) :: phistar, phiold
    real(dp), dimension(:,:,:), allocatable :: phitmp, phiR, phi_old
    real(dp), parameter :: zerodp = 0._dp

    
    Constant_Potential = getinput%int("Constant_Potential", defaultvalue=0) ! Ade : 1=true, 0=false
    if( tot_sol_charge==0.0_dp .and. Constant_Potential==0) then
        if(.not.allocated(phi)) call allocateReal3D(phi)
        phi = 0.0_dp
        return ! phi has been computed, go on !
    end if
    
    call allocateReal3D( phitmp )
    phitmp = 0.0_dp
    if (.not. allocated(LaplacianOfPhi)) call allocateReal3D( LaplacianOfPhi )
    LaplacianOfPhi = 0.0_dp
    if (.not. allocated(phi_old)) call allocateReal3D(phi_old)
    phi_old = phi

    anormf = sum(abs(phi)) 
    if(time==0 .or. t==1) anormf = anormf0

    ! Ben: corresponds to (4*pi*bjl)*(cs^2 /2) in Eq (17) of PRE64, 061507) ok for cs^2 = 1/2 for D3Q18 only
    factor = 4.0_dp*pi*bjl*kbt/2.0_dp

    open(105, file='output/anorm.dat')

    convergenceloop: do iter=1, maxiterations 

        anorm = 0.0_dp ! cumulative diff between new and old phi
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
                        phistar = 0.0_dp
                        phiold = phi(i,j,k)                     ! Ade : used to compute the Laplacian of Phi
                        do l= lbm%lmin, lbm%lmax 
                            pmin = pbc(i-lbm%vel(l)%coo(x),x)  
                            qmin = pbc(j-lbm%vel(l)%coo(y),y)  
                            rmin = pbc(k-lbm%vel(l)%coo(z),z)
                            phistar = phistar + lbm%vel(l)%a0 * phi(pmin,qmin,rmin)   
                        end do
                        phistar = phistar + factor*(c_plus(i,j,k)-c_minus(i,j,k)) ! see PRE64, Horbach: Eq. 17
                        phitmp(i,j,k) = omega*phistar +(1.0_dp-omega)*phiold
                        anorm = anorm + abs(phistar-phiold)
                    endif
                end do
            end do
         end do
         where( .not. node%isFixedPotential ) phi = phitmp

    
        dphi = 0.0_dp ! dphi is the difference between the electric potential between the beginning and end of the iteration.
        ! count the number of times the array is not zero
        n1 = count(abs(phi_old) > threshold)

        ! at each node, if phiold is different from 0, compute the normalized relative difference 
        ! between new and old potential
        if(n1/=0) dphi = sum( abs(  (phi - phi_old)/phi_old ), mask= abs(phi_old)>threshold) / real(n1,kind=dp) 
        phi_old = phi

        ! BR : the convergence criteria should be thought more carefully
        if(anorm <= threshold*anormf) then
            exit convergenceloop
        !else if(iter>1 .and. dphi<1.0d-8) then
        else if(iter>1 .and. dphi<target_error%target_error_sor ) then
            exit convergenceloop
        end if

        ! report every 100 steps
        if( verbose ) then
          if( modulo(iter, 100)==0 ) then
            print*, iter,anorm,threshold*anormf
            write(105,*) '----------------------------------------------'
            write(105,*) 'anormf = ', anormf
            write(105,*) '----------------------------------------------'
          end if
        end if

    end do convergenceloop
    close(105)

    ! tell user if maximum convergence steps is reached, ie if no convergence is found
    if( iter >= maxiterations ) stop 'maximum iterations 500 000 reached without convergence in sor'

    ! compute charge corresponding to this distribution of potential
    !     this is used to compute the charge on the electrodes (when relevant)
    do k = kmin, kmax
        do j = jmin, jmax
            do i = imin, imax
                phistar = 0.0_dp
                do l= lbm%lmin, lbm%lmax               
                    pmin = pbc(i-lbm%vel(l)%coo(x),x) 
                    qmin = pbc(j-lbm%vel(l)%coo(y),y) 
                    rmin = pbc(k-lbm%vel(l)%coo(z),z)
                    phistar = phistar + lbm%vel(l)%a0 * phi(pmin,qmin,rmin)   
                end do
                LaplacianOfPhi(i,j,k) = -(phistar-phi(i,j,k))/factor
            end do
        end do
    end do

    !if( iter > 1 ) then
    if( iter > 1 .and. verbose ) then
        if(anorm <= threshold*anormf) then
            !print*,'SOR converged in',iter-1,'steps with anormf0 =', anormf0,' because anorm <= threshold*anormf'
            print*,'SOR converged in',iter-1,'steps'
        else if(iter>1 .and. dphi<target_error%target_error_sor ) then
            !print*,'SOR converged in',iter-1,'steps with anormf0 =', anormf0,' because dphi < 1.0d-8'
            print*,'SOR converged in',iter-1,'steps'
        end if
    end if

    if( t==1 ) then 
        open(109, file='output/PHIsor.dat')
        do k = kmin, kmax
            write(109,*) k, SUM(phi(:,:,k))
        end do
        close(109)
    endif


end subroutine sor

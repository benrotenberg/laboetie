module module_transient_regime
  implicit none
  private
  public transient_regime

  contains


  !##########################################################################
  SUBROUTINE transient_regime( jx, jy, jz)

    USE precision_kinds, only: dp
    USE system, only: fluid, supercell, node, lbm, n, pbc, solute_force, phi, c_plus, c_minus, LaplacianOfPhi, &
                      el_curr_x, el_curr_y, el_curr_z, t, ion_curr_x, ion_curr_y, ion_curr_z, elec_slope 
    use module_convergence
    use module_collision, only: collide
    use module_input, only: getinput, verbose
    USE constants, only: x, y, z
    use module_bounceback, only: bounceback
    use module_propagation, only: propagation
    use module_advect, only: advect
    use module_external_forces, only: apply_external_forces

    implicit none
    real(dp), intent(inout), dimension(:,:,:) :: jx, jy, jz
    integer :: i,j,k,l, lmin, lmax, ios, lx, ly, lz, Constant_Potential ! Ade : I have declared t in system
    integer :: supercellgeometrylabel, tfext, GL
    integer :: print_frequency, input_print_frequency, print_files_frequency
    integer(kind(fluid)), allocatable, dimension(:,:,:) :: nature
    real(dp) :: fext_tmp(3)
    real(dp) :: l2err, target_error_mass_flux, target_error_charge, max_error_charge
    real(dp) :: Jxx, Jyy, Jzz, Somma4, Somma5, Somma6, FixedPotentialUP, FixedPotentialDOWN
    real(dp), allocatable, dimension(:,:,:) :: density, jx_old, jy_old, jz_old, fextx, fexty, fextz, F1, F2, F3
    real(dp), allocatable, dimension(:) :: a0, a1
    integer, allocatable, dimension(:) :: cx, cy, cz
    logical :: compensate_f_ext
    logical :: mass_flux_convergence_IsReached, mass_flux_convergence_IsReached_without_fext, mass_flux_convergence_IsReached_with_fext
    logical :: charge_convergence_IsReached
    REAL(dp), PARAMETER :: eps=EPSILON(1._dp)
    LOGICAL :: write_total_mass_flux
    integer, allocatable :: il(:,:), jl(:,:), kl(:,:)
    character(200) :: ifile
    LOGICAL :: RestartPNP, i_exist
    integer :: maxEquilibrationTimestep, LW, UW, LB1, UB1, LB2, UB2 !, Starte 
    real(dp), parameter :: zerodp = 0._dp



    open(325, file = "output/phi.dat")
    open(395, file='output/elec_ion_curr.dat')
    open(1324, file = "output/SFZ.dat")
    if (verbose) then
        open(1316, file = "output/SFX.dat")
        open(1323, file = "output/SFY.dat")
    endif

    Constant_Potential = getinput%int("Constant_Potential", 0) ! Ade : 1=true, 0=false
    if (Constant_Potential.eq.1) then
      open(391, file = "output/electrode_charge.dat")
    end if


    RestartPNP = getinput%log("RestartPNP", defaultValue=.TRUE.)
    if( .not. RestartPNP) then
        inquire( file='output/phi.bin'    , exist=i_exist); if( .not. i_exist) error stop "output/phi.bin does not exist"
        inquire( file='output/c_plus.bin' , exist=i_exist); if( .not. i_exist) error stop "output/c_plus.bin does not exist"
        inquire( file='output/c_minus.bin', exist=i_exist); if( .not. i_exist) error stop "output/c_minus.bin does not exist"

        open(731, file='output/phi.bin'    , form="unformatted"); read(731) phi    ; close(731)
        open(731, file='output/c_plus.bin' , form="unformatted"); read(731) c_plus ; close(731)
        open(731, file='output/c_minus.bin', form="unformatted"); read(731) c_minus; close(731)
    end if


    lmin = lbm%lmin
    lmax = lbm%lmax

    supercellgeometrylabel = supercell%geometry%label ! -1 for solid free cell

    lx = supercell%geometry%dimensions%indiceMax(x)
    ly = supercell%geometry%dimensions%indiceMax(y)
    lz = supercell%geometry%dimensions%indiceMax(z)

    !
    ! Print info to terminal every that number of steps
    !
    ! start with 1, then increase to 10, 100, 1000 -- unless the input frequency is smaller than 1000
    input_print_frequency = getinput%int('print_frequency', defaultvalue=1000, assert=">0" )
    print_frequency = 1


    !
    ! WRITE velocity profiles to terminal every that number of steps
    !
    print_files_frequency = getinput%int( "print_files_frequency", defaultvalue = huge(1) )


    !
    ! Target values to monitor convergence
    !

    !   on mass flux (note: this on L1 norm for difference of absolute fluxes)
    target_error_mass_flux =target_error%target_error_mass_flux 

    !   on potential and ion concentrations (note: this on relative values)
    target_error_charge = target_error%target_error_charge


    allocate( density(lx,ly,lz), source=sum(n,4), stat=ios)
    if (ios /= 0) stop "density: Allocation request denied"

    if(.not.allocated(solute_force)) allocate(solute_force(lx,ly,lz,x:z),source=0.0_dp)
    allocate( jx_old (lx,ly,lz), source=0._dp )
    allocate( jy_old (lx,ly,lz), source=0._dp )
    allocate( jz_old (lx,ly,lz), source=0._dp )

    OPEN(66, FILE="output/mass-flux_profile_along_z.dat")
    WRITE(66,*) "# z, <ρ.v_x>_{x,y}, <ρ.v_y>_{x,y}, <ρ.v_z>_{x,y}"

    if ( verbose ) then
      OPEN(67, FILE="output/mass-flux_profile_along_y.dat")
      OPEN(68, FILE="output/mass-flux_profile_along_x.dat")
      WRITE(67,*) "# y, <ρ.v_x>_{x,z}, <ρ.v_y>_{x,z}, <ρ.v_z>_{x,z}"
      WRITE(68,*) "# x, <ρ.v_x>_{y,z}, <ρ.v_y>_{y,z}, <ρ.v_z>_{y,z}"

      OPEN(56, FILE="output/mean-density_profile_along_z.dat")
      OPEN(57, FILE="output/mean-density_profile_along_y.dat")
      OPEN(58, FILE="output/mean-density_profile_along_x.dat")
    end if

    allocate( nature (lx,ly,lz), source=node%nature)
    allocate( fextx(lx,ly,lz), source=zerodp)
    allocate( fexty(lx,ly,lz), source=zerodp)
    allocate( fextz(lx,ly,lz), source=zerodp)
    !--------------- Ade ----------------------
    allocate( F1(lx,ly,lz), source=zerodp)
    allocate( F2(lx,ly,lz), source=zerodp)
    allocate( F3(lx,ly,lz), source=zerodp)
    !--------------- Ade ----------------------

    fext_tmp = zerodp ! this is important and spagetty like... please read carefully before modifying this line
    allocate( cx(lmax), source=lbm%vel(:)%coo(1))
    allocate( cy(lmax), source=lbm%vel(:)%coo(2))
    allocate( cz(lmax), source=lbm%vel(:)%coo(3))
    allocate( a0(lmax), source=lbm%vel(:)%a0)
    allocate( a1(lmax), source=lbm%vel(:)%a1)

    !
    ! Tabulate the index of the node one finishes if one starts from a node and a velocity index l per direction
    !
    allocate( il(lbm%lmin:lbm%lmax, 1:lx), stat=ios)
    if (ios /= 0) stop "il: Allocation request denied"
    allocate( jl(lbm%lmin:lbm%lmax, 1:ly), stat=ios)
    if (ios /= 0) stop "jl: Allocation request denied"
    allocate( kl(lbm%lmin:lbm%lmax, 1:lz), stat=ios)
    if (ios /= 0) stop "kl: Allocation request denied"
    do l= lmin, lmax
        il(l,:) = [( pbc(i+cx(l),x) ,i=1,lx )]
        jl(l,:) = [( pbc(j+cy(l),y) ,j=1,ly )]
        kl(l,:) = [( pbc(k+cz(l),z) ,k=1,lz )]
    end do

    mass_flux_convergence_IsReached_without_fext = .false.
    mass_flux_convergence_IsReached_with_fext = .false.

    mass_flux_convergence_IsReached = .false.
    charge_convergence_IsReached = .false.
    
    compensate_f_ext = getinput%log( "compensate_f_ext", defaultvalue = .false.)
    if(compensate_f_ext) then
        open(79, file = "./output/v_centralnode.dat")
    end if

    write_total_mass_flux = getinput%log( "write_total_mass_flux", defaultvalue = .false.)
        if( write_total_mass_flux ) open(65, file = "./output/total_mass_flux.dat" )



    PRINT*
    PRINT*,'Lattice Boltzmann'
    PRINT*,'================='
    PRINT*,'       step'
    PRINT*,'       ----'


    ! ADE : We initialise tfext, which is the time from when the force f_ext is applied
    !       upon a certain number of nodes. 
    tfext = HUGE(tfext)

    ! Initialize the measures of convergence
    l2err = -1
    max_error_charge = -1

    ! ADE : Max had originally read this variable in PNP, but in case we restart a simulation after PNP, 
    !       it is better to read it again in this module
    elec_slope = getinput%dp3("elec_slope")


    ! ADE: constants needed only for constant Potential simulations
    if (Constant_Potential.eq.1) then

        ! index of interfacial electrode planes
        LW = getinput%int("SolidSheetsEachSideOfSlit", defaultvalue=1, assert=">0") ! defines the index of the last lower wall in slit geometry
        UW = lz + 1 - LW                                                            ! defines the index of the first upper wall in slit geometry

        LB1= getinput%int("LB1", defaultvalue=1)  ! Lower Boundary at z = LW
        UB1= getinput%int("UB1", defaultvalue=ly) ! Upper Boundary at z = LW
        LB2= getinput%int("LB2", defaultvalue=1)  ! Lower Boundary at z = UW
        UB2= getinput%int("UB2", defaultvalue=ly) ! Upper Boundary at z = UW

        FixedPotentialUP = getinput%dp('FixedPotentialUP',0._dp)
        FixedPotentialDOWN = getinput%dp('FixedPotentialDOWN',0._dp)

    end if

    !
    ! Main temporal loop
    !
    do t = 1, huge(t)

        ! 
        ! Print some info to stdout: time step and current/target values for convergence.
        !
        if( modulo(t, print_frequency) == 0) then
            if( t==10    )   print_frequency = min( input_print_frequency, 10    )
            if( t==100   )   print_frequency = min( input_print_frequency, 100   )
            if( t==1000  )   print_frequency = min( input_print_frequency, 1000  )

            WRITE(*,"(4X,I8,A,2(E12.6,A,E12.6,A))") t, " crit. on mass flux ",l2err," ( target ",target_error_mass_flux," ) and on charge ", max_error_charge, " ( target ", target_error_charge, " )"
            !PRINT*, t, "crit. on mass flux ",real(l2err),"(target",real(target_error_mass_flux,4),")"
            !PRINT*, "                  on charge    ",real(max_error_charge),"(target",real(target_error_charge,4),")"
        end if

        !
        ! Print velocity profiles, etc.
        !
        if( modulo(t, print_files_frequency) == 0 ) then

            ! most frequent case
            WRITE(66,*)"# timestep",t
            DO k=1,lz
                WRITE(66,*) k, SUM(jx(:,:,k)), SUM(jy(:,:,k)), SUM(jz(:,:,k))
            END DO
 
            ! only if we really want this
            if (verbose) then
               WRITE(67,*)"# timestep",t
               WRITE(68,*)"# timestep",t
               WRITE(56,*)"# timestep",t
               WRITE(57,*)"# timestep",t
               WRITE(58,*)"# timestep",t
               DO k=1,lz
                   WRITE(56,*) k, SUM(density(:,:,k))/ MAX( COUNT(density(:,:,k)>eps)  ,1)
               END DO
               DO k=1,ly
                   WRITE(57,*) k, SUM(density(:,k,:))/ MAX( COUNT(density(:,k,:)>eps)  ,1)
                   WRITE(67,*) k, SUM(jx(:,k,:)), SUM(jy(:,k,:)), SUM(jz(:,k,:))
               END DO
               DO k=1,lx
                   WRITE(58,*) k, SUM(density(k,:,:))/ MAX( COUNT(density(k,:,:)>eps)  ,1)
                   WRITE(68,*) k, SUM(jx(k,:,:)), SUM(jy(k,:,:)), SUM(jz(k,:,:))
               END DO
               WRITE(66,*)
               WRITE(67,*)
               WRITE(68,*)
            end if

        end if

        !
        ! Compensate_f_ext is an option to have a local force one could apply on a given node (only the central node for now)
        ! compensated by a continuum background force, a little bit like a compensating electric field in a charged supercell.
        !
        ! Ade: By doing so, we can apply a force fx or fy in order to analyse and observe the velocity streamlines around
        ! the particle in the output file v_centralnode.dat.
        !
        
        ! f_ext is obtained reading input file lb.in (=> pressure gradient)
        ! solute_force is computed in smolu.f90

        F1(:,:,:)  = fextx(:,:,:) + solute_force(:,:,:,1)
        F2(:,:,:)  = fexty(:,:,:) + solute_force(:,:,:,2)
        F3(:,:,:)  = fextz(:,:,:) + solute_force(:,:,:,3) 

        !##################
        !# Collision step #
        !##################
        call collide(n, density, jx, jy, jz, F1, F2, F3)


        !###############
        !# Bounce Back # to simplify propagation
        !###############
        call bounceback(n, nature)


        !###############
        !# PROPAGATION #
        !###############
        call propagation(n, lmin, lmax, lx, ly, lz, il, jl, kl)


        !###############
        !# CHECK POPULATIONS ARE POSITIVE
        !###############
        IF( ANY(n<0) ) ERROR STOP "In transient_regime, some population n(x,y,z,vel) gets negative !"


        !###############
        !# UPDATE DENSITIES
        !# Remember: Density(x,y,z) is the first moment of the populations(x,y,z,v): it is the weighted sum over the velocities
        !################
        density = SUM(n,4)

        !
        ! backup moment density (velocities) to test convergence at the end of the timestep
        !
        jx_old = jx
        jy_old = jy
        jz_old = jz

        !################################
        !# UPDATE OF IONS AND POTENTIAL #
        !################################

        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! This modification is for constant potential simulations. One of the
        ! idea is that we can restart from equilibration
        ! using the phi.bin and c_plus.bin files in order to simulate a
        ! discharge, setting Delta V = 0 for this part.
        ! N.B. this part could possibly be put before the temporal loop and
        ! would save up compuational time. -> I need to check this
        if (Constant_Potential.eq.1) then
         do k=1,lz
            do j=1,ly
                do i=1,lx
                    if( k<=LW .and. j>=LB1 .and. j<=UB1) then
                        node(i,j,k)%isFixedPotential = .true.
                        phi(i,j,k) = FixedPotentialDOWN
                    else if( k>=UW .and. j>=LB2 .and. j<=UB2 ) then
                        node(i,j,k)%isFixedPotential = .true.
                        phi(i,j,k) = FixedPotentialUP
                    else
                        node(i,j,k)%isFixedPotential = .false.
                    end if
                end do
            end do
         end do
        end if
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        ! backup potential and solute concentrations from last step
        !     to check convergence at the end of this iteration

        call backup_phi_c_plus_c_minus          


        ! compute the electric potential phi from charge distribution
        !     with the Successive OverRelation method (SOR)

        call sor  


        ! correction to the mass fluxes jx, jy, jz due to the local force F1, F2, F3
        !     + output mass fluxes if requested

        call update_solventCurrent( jx, jy, jz, n, cx, cy, cz, F1, F2, F3, t, write_total_mass_flux)


        ! update ion concentrations due to advection, based on local fluxes jx, jy, jz

        call advect( density, jx, jy, jz )      


        ! add external electric field to that arising from charge distribution (computed in SOR)
        !     to update concentrations via migration in smolu

        call electrostatic_pot           


        ! update ion concentrations due to diffusion an migration 

        call smolu


        ! check that charge is conserved

        call check_charge_conservation 



        !#####################
        !# check convergence #
        !#####################

        ! check the convergence of the mass fluxes with respect to previous step

        call check_mass_flux_convergence(t, target_error_mass_flux, l2err, jx, jy, jz, jx_old, jy_old, jz_old, mass_flux_convergence_IsReached )

        call check_charge_distribution_equilibrium( t, target_error_charge, max_error_charge, charge_convergence_IsReached)


        ! First, we reach convergence without external forces.
        ! Then, we turn on the external forces and iterate again until convergence
        !
        if(mass_flux_convergence_IsReached) then

            ! If this was only without the external force, now we apply it
            if( .not.mass_flux_convergence_IsReached_without_fext ) then
                mass_flux_convergence_IsReached_without_fext = .true.
                mass_flux_convergence_IsReached = .false.
                print*,'Mass flux convergence reached without fext after ',t,' steps... Now applying it'

       !############################################
       !# Apply external force if it wasn't already there
       !############################################
                tfext=t+1
                call apply_external_forces( fextx, fexty, fextz, nature)

            ! If we now also have converged with the external force
            else if( mass_flux_convergence_IsReached_without_fext ) then
                mass_flux_convergence_IsReached_with_fext = .true.

                ! We also check whether convergence has also been reached for
                !   electrostatic potential and ionic concentrations
                ! If yes, then we are done and we can exit the loop over time steps
                !
                if ( charge_convergence_IsReached .and. t>2 ) exit

            end if

        end if



        !############################################
        !# Apply external force
        !############################################
        !if( mass_flux_convergence_IsReached ) then

        !    ! if you are already converged without then with f_ext then exit time loop. Stationary state is found.
        !    if( mass_flux_convergence_IsReached_without_fext .and. mass_flux_convergence_IsReached_with_fext .and. t>2) then
        !        if( charge_convergence_IsReached ) exit    ! we are done: exit loop over time steps

        !    ! if you have already converged without fext, but not yet with fext, then enable fext
        !    else if(mass_flux_convergence_IsReached_without_fext .and. .not.mass_flux_convergence_IsReached_with_fext) then
        !        tfext=t+1
        !        call apply_external_forces( fextx, fexty, fextz, nature)
        !        print*,'Applying external force: ', fextx, fexty, fextz
        !    end if  

        !end if 



        !############################################
        !# Compute (and output) the charge on the electrodes for parallel plane capacitor
        !############################################

        ! ADE : here the Somma* indicate the charge accumulated on different layers of the plate geometry. 
        ! Only the layer corresponding to the interfacial solid node is in fact needed. However, since my 
        ! python scripts are implemented in order to reads all of these entries I shall leave till the next
        ! Laboetie user

        ! note BR : indeed this should be rationalized 
        !       (only interfacial planes, and more transparent variable name)  
        if (Constant_Potential.eq.1) then
          DO k=lz-2, lz
            IF (k==lz-2) Somma4 = sum( LaplacianOfPhi(:,:,k) )
            IF (k==lz-1) Somma5 = sum( LaplacianOfPhi(:,:,k) )
            IF (k==lz)   Somma6 = sum( LaplacianOfPhi(:,:,k) )
          ENDDO
          write(391,*) t, Somma4, Somma5, Somma6
        end if

        write(395,*) t, el_curr_x, el_curr_y, el_curr_z, ion_curr_x, ion_curr_y, ion_curr_z
   

     
   end do ! end of temporal loop


   !##########################################################
   !# Store phi, c_plus, c_minus in binary file for restarts #
   !##########################################################

   open(731, file='output/phi.bin'    , form="unformatted"); write(731) phi    ; close(731)
   open(731, file='output/c_plus.bin' , form="unformatted"); write(731) c_plus ; close(731)
   open(731, file='output/c_minus.bin', form="unformatted"); write(731) c_minus; close(731)


   !########################
   !# ADE: post-processing #
   !########################

   ! 1. Print velocity field
   WRITE(66,*)"# Steady state with convergence threshold on mass flux", REAL(target_error_mass_flux)

   if ( verbose ) then
      WRITE(67,*)"# Steady state with convergence threshold on mass flux", REAL(target_error_mass_flux)
      WRITE(68,*)"# Steady state with convergence threshold on mass flux", REAL(target_error_mass_flux)
      WRITE(56,*)"# Steady state with convergence threshold on mass flux", REAL(target_error_mass_flux)
      WRITE(57,*)"# Steady state with convergence threshold on mass flux", REAL(target_error_mass_flux)
      WRITE(58,*)"# Steady state with convergence threshold on mass flux", REAL(target_error_mass_flux)
   end if  

   GL = getinput%int("geometryLabel", defaultvalue=0) ! if GL=-1 =>bulk case

   if (GL==2) then                      ! Cylindrical geometry

      if ( verbose ) then
         Jxx = 0
         Jyy = 0
         Jzz = 0
         DO i=1,lx
           Jxx = jx(i,max(ly/2,1),max(lz/2,1))
           Jyy = jy(i,max(ly/2,1),max(lz/2,1))
           Jzz = jz(i,max(ly/2,1),max(lz/2,1))
           WRITE(67,*) i, Jxx, Jyy, Jzz
         END DO
     end if

   else

     DO k=1,lz
        WRITE(66,*) k, SUM(jx(:,:,k)), SUM(jy(:,:,k)), SUM(jz(:,:,k))
     END DO

     if ( verbose ) then
        DO k=1,lz
           WRITE(56,*) k, SUM(density(:,:,k))/ MAX( COUNT(density(:,:,k)>eps) ,1)
        END DO
        DO k=1,ly
           WRITE(67,*) k, SUM(jx(:,k,:)), SUM(jy(:,k,:)), SUM(jz(:,k,:))
           WRITE(57,*) k, SUM(density(:,k,:))/ MAX( COUNT(density(:,k,:)>eps) ,1)
        END DO
        DO k=1,lx
           WRITE(68,*) k, SUM(jx(k,:,:)), SUM(jy(k,:,:)), SUM(jz(k,:,:))
           WRITE(58,*) k, SUM(density(k,:,:))/ MAX( COUNT(density(k,:,:)>eps) ,1)
        END DO
     end if

   end if

 
   ! 2. Solute Force (averages over xy slabs, as a function of z)
   DO k=1,lz
     WRITE(1324,*) k, SUM(solute_force(:,:,k,3))/(lx*ly) 
   ENDDO 
   if ( verbose ) then
     DO k=1,lz
        WRITE(1316,*) k, SUM(solute_force(:,:,k,1))/(lx*ly)  
        WRITE(1323,*) k, SUM(solute_force(:,:,k,2))/(lx*ly) 
     ENDDO 
   end if

   ! 3. Potential PHI (averages over xy slabs, as a function of z)
   DO k=1,lz
      write(325,*) k, SUM(phi(:,:,k))/(lx*ly) 
   END DO


   !
   ! Print velocity 2D profiles
   !
   OPEN(69, FILE="output/mass-flux_field_2d_at_x.eq.1.dat")
   DO concurrent( j=1:ly, k=1:lz); WRITE(69,*) j, k, jy(1,j,k), jz(1,j,k); END DO
   DO j=1,ly
      DO k=1,lz
          WRITE(69,*) j, k, jy(1,j,k), jz(1,j,k)
      END DO
   END DO
   CLOSE(69)

   ! close files that may still be open
   close(325)
   close(395)
   close(66)
   close(1324)

   if ( verbose ) then
      close(68)
      close(58)
      close(67)
      close(57)
      close(56)
      close(1316)
      close(1323)
   end if


   if( write_total_mass_flux ) close(65)
   if(compensate_f_ext) close(79)
   if (Constant_Potential.eq.1) close(391)

end subroutine transient_regime




!##########################################################################
subroutine update_solventCurrent( jx, jy, jz, n, cx, cy, cz, F1, F2, F3, timestep, write_total_mass_flux)
    use precision_kinds, only: dp
    implicit none
    real(dp), intent(inout), dimension(:,:,:) :: jx, jy, jz
    real(dp), intent(in) :: n(:,:,:,:), F1(:,:,:), F2(:,:,:), F3(:,:,:)
    integer, intent(in) :: cx(:), cy(:), cz(:)
    integer, intent(in) :: timestep
    integer :: lmin, lmax, l
    logical, intent(in) :: write_total_mass_flux
    lmin = lbound(cx,1)
    lmax = ubound(cx,1)
    ! update momentum densities after the propagation
    ! this is completely local in space and my be parallelized very well
    jx = F1/2._dp
    jy = F2/2._dp
    jz = F3/2._dp
    !$OMP PARALLEL DO DEFAULT(NONE)&
    !$OMP PRIVATE(l)&
    !$OMP SHARED(lmin,lmax,n,cx,cy,cz)&
    !$OMP REDUCTION(+:jx)&
    !$OMP REDUCTION(+:jy)&
    !$OMP REDUCTION(+:jz)
    do l=lmin,lmax
        jx = jx +n(:,:,:,l)*cx(l)
        jy = jy +n(:,:,:,l)*cy(l)
        jz = jz +n(:,:,:,l)*cz(l)
    end do
    !$OMP END PARALLEL DO
    if( write_total_mass_flux ) write(65,*) timestep, real([  sum(jx), sum(jy), sum(jz)  ])

  end subroutine update_solventCurrent
  !##########################################################################


  !##########################################################################
  subroutine check_mass_flux_convergence( timestep, target_error_mass_flux, l2err, jx, jy, jz, jx_old, jy_old, jz_old, mass_flux_convergence_IsReached )
    use precision_kinds, only: dp
    use module_input, only: getinput
    implicit none
    logical :: itIsOpen
    real(dp) :: maxjx, maxjy, maxjz
    real(dp), intent(inout) :: l2err
    real(dp), intent(in), dimension(:,:,:) :: jx, jy, jz, jx_old, jy_old, jz_old
    logical, intent(out) :: mass_flux_convergence_IsReached
    character(18), parameter :: filename = "./output/l2err.dat"
    integer, intent(in) :: timestep
    real(dp), intent(in) :: target_error_mass_flux
    integer :: maxEquilibrationTimestep

    maxEquilibrationTimestep = getinput%int( 'maxEquilibrationTimestep' , defaultvalue = huge(1))
    if( timestep > maxEquilibrationTimestep ) then
        error stop "You reached the maximum number of iterations allowed before convergence, as defined in input"
    end if

    ! We start by equilibrating the densities, fluxes and other moments without the external forces.
    ! This equilibration step is done up to equilibration 
    ! or if the timestep t gets higher than a maximum value called maxEquilibrationTimestep
    inquire(file=filename, opened=itIsOpen)
    if(.not. itIsOpen) open(13,file=filename)

    maxjx = maxval(abs(jx-jx_old))
    maxjy = maxval(abs(jy-jy_old))
    maxjz = maxval(abs(jz-jz_old))
    l2err = max(maxjx, maxjy, maxjz)            ! note BR: this is not l2err but l1err
    write(13,*) timestep, l2err

    if( l2err <= target_error_mass_flux .and. timestep > 1 ) then
        mass_flux_convergence_IsReached = .true.
    else
        mass_flux_convergence_IsReached = .false.
    end if

  end subroutine check_mass_flux_convergence
  !##########################################################################

end module module_transient_regime

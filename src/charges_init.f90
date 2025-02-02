! Here we initiate the salt local densities (solutes + and -)

subroutine charges_init

    use precision_kinds, only: i2b, dp
    use system, only: lambda_D, c_plus, c_minus, phi, charge_distrib, tot_sol_charge, node,&
                       fluid, solid, bjl, kBT, rho_0, D_plus, D_minus, supercell, A_zero 
    use constants, only: x, y, z
    use module_input, only: getinput
    use myallocations

    implicit none
    integer  :: count_solid, count_fluid, count_solid_int, MyGeometry, i,j,k, lx, ly, lz, ENNE,&
                capacitor, m, slit_sol, LB1, UB1, LB2, UB2, LW, UW, width1, width2, InnerRadius, OuterRadius
    real(dp) :: tot_sol_charge_solid, tot_sol_charge_fluid, SF, PrefactLP1, PrefactLP2, PrefactLP3, Alpha
    real(dp) :: in_c_plus_solid, in_c_plus_fluid, in_c_minus_solid, in_c_minus_fluid, FixedPotentialUP, &
                FixedPotentialDOWN
    real(dp), dimension(2) :: rnode, rorigin ! coordinates of each node and center of cylinder in x,y coordinates
    REAL(dp), PARAMETER :: eps=EPSILON(1._dp)
    integer :: Constant_Potential
    real(dp), parameter :: pi = acos(-1._dp)

    !open(279, file='output/c_ions_from_charges_init.dat')

    !
    ! Is there any charge into the solute ?
    !
    bjl = getinput%dp('bjl',0._dp)
    capacitor = getinput%int('capacitor', defaultvalue=0) ! 1 := simulate a capacitor. 0:=we're not simulating a capacitor

    IF( bjl <= epsilon(1.0_dp) ) return

    lx = supercell%geometry%dimensions%indiceMax(x)
    ly = supercell%geometry%dimensions%indiceMax(y)
    lz = supercell%geometry%dimensions%indiceMax(z)
    
    MyGeometry = getinput%int('geometryLabel',-1) ! Ade: 15/03/2017
    slit_sol = getinput%int('slit_sol',0) ! 1= true 0 = false
    tot_sol_charge = getinput%dp('tot_sol_charge',0._dp)
    lambda_d = getinput%dp('lambda_D',0._dp)

    if( abs(lambda_D) <= epsilon(1._dp) ) then
      print*, 'salt free fluid'
      rho_0 = 0.0_dp
    else
      rho_0 = 1.0_dp / (4.0_dp*pi*bjl*lambda_D**2)
    end if

    ! count of solid, fluid and interfacial nodes
    count_solid = count(node%nature==solid)
    count_solid_int = count( node%nature == solid .and. node%isInterfacial )
    count_fluid = count(node%nature==fluid)

    ! read where are distributed the charges
    ! call read_charge_distrib

    !Constant_Potential = getinput%log("Constant_Potential", .FALSE.)
    Constant_Potential = getinput%int("Constant_Potential", 0) ! Ade : 1=true, 0=false
    print*, 'Constant_Potential = ', Constant_Potential


    IF( ABS(tot_sol_charge) > EPSILON(1._dp) .OR. ABS(lambda_d)> EPSILON(1._dp) .OR. Constant_Potential==1) THEN ! Ade : this if statement should be removed 27/01/2017

    ! init ion (solute) concentrations
    if (.not. allocated(c_plus)) call allocateReal3D(c_plus)
    if (.not. allocated(c_minus)) call allocateReal3D( c_minus)

    ! init potential
    ! Ade : 19/05/17
    ! Constant Potential Simulation
    LW = getinput%int("SolidSheetsEachSideOfSlit", defaultvalue=1, assert=">0") ! defines the index of the last lower wall in slit geometry
    UW = lz + 1 - LW ! defines the index of the first upper wall in slit geometry
    LB1= getinput%int("LB1", defaultvalue=1) ! Lower Boundary at z = LW
    UB1= getinput%int("UB1", defaultvalue=ly) ! Upper Boundary at z = LW
    LB2= getinput%int("LB2", defaultvalue=1) ! Lower Boundary at z = UW
    UB2= getinput%int("UB2", defaultvalue=ly) ! Upper Boundary at z = UW
    width1 = getinput%int("width1",3)
    width2 = getinput%int("width2",1)
    InnerRadius = int(real(lx-width1,dp)/2.0_dp)
    OuterRadius = int(real(lx-width2,dp)/2.0_dp)
    FixedPotentialUP = getinput%dp('FixedPotentialUP',0._dp)
    FixedPotentialDOWN = getinput%dp('FixedPotentialDOWN',0._dp)
    call allocateReal3D( phi) !allocate( phi(lx,ly,lz), source=0.0_dp )
    phi = 0._dp
    !if(.NOT.Constant_Potential) then ! Two distinct regions where. Constant_Potential = true here
    if(Constant_Potential==1) then ! Two distinct regions where. Constant_Potential = true here
      if(MyGeometry==2 .or. MyGeometry==15) then ! Cylindrical Geometry
        rorigin = [ real(lx+1,dp)/2.0_dp, real(ly+1,dp)/2.0_dp ]
        print*, 'Coaxial capacitor'
        print*, 'rorigin in charges_init is =', rorigin
        do i = 1, lx
          do j = 1, ly
          rnode = [real(i,dp),real(j,dp)] - rorigin
            do k = 1, lz
              if(int(norm2(rnode)) <= (InnerRadius-1)) then
                  phi(i,j,k) = FixedPotentialDOWN
                  node(i,j,k)%isFixedPotential = .true. ! Ade : we defined the following nodes as nodes
              else if(int(norm2(rnode)) >= OuterRadius) then
                  phi(i,j,k) = FixedPotentialUP
                  node(i,j,k)%isFixedPotential = .true.
              endif
            enddo
          enddo
        enddo
      else ! Slit geometry
        print*, '=================================================================================================='
        print*, 'Parallel plate capacitor'
        do i = 1, lx
          do j = 1, ly
            do k = 1, lz
              if(k<=LW) then
                if(j>=LB1 .and. j<=UB1) then
                    phi(i,j,k) = FixedPotentialDOWN
                    node(i,j,k)%isFixedPotential = .true. ! Ade : we defined the following nodes as nodes
                                                          ! which have a fixed potential. This is nature of nodes
                                                          ! is then used in SOR
                endif
              else if(k>=UW) then
                if(j>=LB2 .and. j<=UB2) then
                    phi(i,j,k) = FixedPotentialUP
                    node(i,j,k)%isFixedPotential = .true.
                endif
              endif
            enddo
          enddo
        enddo
      endif
    endif
    !print*,'==========================================================='
    !print*, 'phi is initialised as'
    !print*,'==========================================================='
    charge_distrib = getinput%char("charge_distrib")
    if( charge_distrib(1:3) /= 'int' .and. charge_distrib(1:3)/='sol') stop 'charge_distrib can only be int or sol for now'

    ! distribute charge, depending on where user asked
    if( charge_distrib(1:3) == 'sol') then
      if(count_solid/=0) then
        tot_sol_charge_solid = tot_sol_charge / count_solid ! charges distributed in all solid nodes
      else
        tot_sol_charge_solid = 0
      end if

    else if( charge_distrib(1:3) == 'int') then
      if(count_solid_int/=0) then
        tot_sol_charge_solid = tot_sol_charge / count_solid_int ! charges distributed in interfacial solid nodes only
      else
        tot_sol_charge_solid = 0
      end if
    end if

    if(count_fluid/=0) then
      tot_sol_charge_fluid = -1.0_dp * tot_sol_charge / count_fluid ! charges distributed in all fluid nodes
    else
      tot_sol_charge_fluid = 0
    end if

    in_c_plus_solid  = +0.5_dp*tot_sol_charge_solid;
    in_c_minus_solid = -0.5_dp*tot_sol_charge_solid;

    if( tot_sol_charge_fluid > 0.0_dp ) then
          in_c_plus_fluid  = 0.5_dp*rho_0 + tot_sol_charge_fluid
          in_c_minus_fluid = 0.5_dp*rho_0
    else
      in_c_plus_fluid  = 0.5_dp*rho_0
      in_c_minus_fluid = 0.5_dp*rho_0 - tot_sol_charge_fluid
    end if

    m = getinput%int("SolidSheetsEachSideOfSlit", defaultvalue=1, assert=">0")
    
    if(abs(tot_sol_charge)>epsilon(1._dp)) then

      if(charge_distrib(1:3)=='int') then
        print*,'The total charge is set ONLY on the solid nodes at the interface'
        print*,' # of interfacial solid sites = ', count_solid_int, ' (charge per node =',tot_sol_charge_solid,')'
        print*,' # of fluid sites             = ', count_fluid,     ' (charge per node =',tot_sol_charge_fluid,')'
      else if(charge_distrib(1:3)=='sol') then
        print*,' # of interfacial solid sites = ',count_solid,' (charge per node =',tot_sol_charge_solid,')'
        print*,' # of fluid sites             = ',count_fluid,' (charge per node =',tot_sol_charge_fluid,')'
        stop 'only surface charge is implemented for now in charge_init.f90'
      else
        stop 'pb in charges_init.f90'
      end if

      print*,'Salt concentration ',0.5*rho_0
      print*,'Initial density values:'
      print*,'  p_solid =',in_c_plus_solid
      print*,'  p_fluid =',in_c_plus_fluid
      print*,'  m_solid =',in_c_minus_solid
      print*,'  m_fluid =',in_c_minus_fluid
      print*,'*********************************************************************'

    end if

    where(node%nature==solid .and. node%isInterfacial )
      c_plus = in_c_plus_solid
      c_minus = in_c_minus_solid
    else where(node%nature==solid .and. .not. node%isInterfacial )
      c_plus = 0.0_dp
      c_minus = 0.0_dp
    else where(node%nature==fluid)
      c_plus = in_c_plus_fluid
      c_minus = in_c_minus_fluid
    end where


    ! This modication was done in order to simulate a parallel plate capacitor with fixed charges.
    !     The charge is positive on the left wall, negative on the right wall
    if (capacitor.EQ.1) then
          c_plus(:,:,m) = in_c_plus_solid
          c_minus(:,:,m) = -in_c_plus_solid
        print*, ' ------------------------------------ '
        print*, 'c_plus wall1 = ', c_plus
          c_plus(:,:,lz-(m-1)) = -in_c_plus_solid
          c_minus(:,:,lz-(m-1)) = in_c_plus_solid
        print*, ' ------------------------------------ '
        print*, 'c_plus wall2 = ', c_plus
    end if 


    ! Attention!!! you should write also the corresponding c_minus
    if (capacitor==1) then
        where(node%nature==solid .and. .not. node%isInterfacial )
          c_plus = 0.0_dp
          c_minus = 0.0_dp
        ! Ade : I added the part below because if I insert two sheets of solid nodes
        ! I need to have the non interfacial nodes to be free of charge
        else where(node%nature==fluid)
          c_plus = 0.0_dp
          c_minus = 0.0_dp
        end where
    end if
    

    if(capacitor.NE.1 .and. slit_sol==1) then
      if( MyGeometry == 1 .or. MyGeometry==12 .or. MyGeometry==13 .or. MyGeometry==14 ) then   ! slit pore geometry
          ! Question BR->ADE: how can it be 12, 13 or 14?
          Alpha = getinput%dp('Alpha',0._dp)  ! Attention!!!!!! This is dangerous. We should probably  
                                              ! compute its value in another subroutine
          ENNE = lz-(2*m) ! Nb of fluid nodes
          if(mod(lz,2) == 0) then
              SF = real(lz)/2
          else
              SF = real(lz)/2 + 0.5
          endif
          if( lambda_d > EPSILON(1._dp)) then ! Low potential condition - salt added
              !SurfArea = (lx*ly)_dp
              PrefactLP1 = 1/(8*PI*bjl*lambda_d**2) 
              PrefactLP2 = 4*PI*tot_sol_charge*bjl*lambda_d/(2*lx*ly)
              PrefactLP3 = sinh(real(ENNE) / ( 2*lambda_d) )
              do i = 1, lx    
                do j = 1, ly   
                  do k = m+1, lz-m ! first and last node are solid
                    c_plus(i,j,k) = PrefactLP1 * ( 1 - PrefactLP2 * cosh( (k-SF)/lambda_d )/(PrefactLP3) ) !* SurfArea
                    c_minus(i,j,k) = PrefactLP1 * ( 1 + PrefactLP2 * cosh( (k-SF)/lambda_d )/(PrefactLP3) )
                  end do
                end do
              end do
          else      ! Normal c_plus density - no salt
              do i = 1, lx      
                do j = 1, ly  
                  do k = m+1, lz-m ! first and last node are solid
                      c_plus(i,j,k) = (Alpha**2/( 2*PI*bjl * ( cos(Alpha*(k-SF)) )**2 )) !* SurfArea
                      c_minus(i,j,k) = 0._dp ! Ade : this is redundant as it has already this value
                  end do
                end do
              end do
          endif
      endif
    end if

    
    !DO k=supercell%geometry%dimensions%indiceMin(z), supercell%geometry%dimensions%indiceMax(z)
    !    WRITE(279,*) k, SUM(c_plus(:,:,k))/(lx*ly), SUM(c_minus(:,:,k))/(lx*ly),&
    !                c_plus(:,:,k), c_minus(:,:,k)
    !ENDDO
    !close(279)


    ! read diffusion coefficients of solute + and solute -
    d_plus = getinput%dp( 'D_plus', defaultvalue=0._dp, assert=">0")
    d_plus = d_plus/A_zero
    
    d_minus = getinput%dp('D_minus', defaultvalue=0._dp, assert=">0")
    d_minus = d_minus/A_zero

END IF
end subroutine charges_init

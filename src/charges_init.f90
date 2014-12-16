! Here we initiate the charged densities (solutes + and -)

subroutine charges_init

  use precision_kinds, only: i2b, dp
  use system, only: lambda_D, c_plus, c_minus, phi, charge_distrib, sigma, node,&
                     fluid, solid, anormf0, bjl, kBT, rho_ch, D_plus, D_minus, supercell
  use constants, only: pi, x, y, z
  use input
  use myallocations
  implicit none
  integer(kind=i2b) :: count_solid, count_fluid, count_solid_int
  real(kind=dp) :: sigma_solid, sigma_fluid
  real(kind=dp) :: in_c_plus_solid, in_c_plus_fluid, in_c_minus_solid, in_c_minus_fluid
  real(kind=dp) :: rho_0 ! don't understand why here it is different from the value read in input file.
  ! init ion (solute) concentrations
  if (.not. allocated(c_plus)) call allocateReal3D(c_plus)
  if (.not. allocated(c_minus)) call allocateReal3D( c_minus)

  ! init potential
  call allocateReal3D( phi) !allocate( phi(lx,ly,lz), source=0.0_dp )
  phi = 0._dp
  ! read
  bjl = input_dp('bjl')
  sigma = input_dp('sigma')
  lambda_d = input_dp('lambda_D')

  if( lambda_D <= epsilon(1._dp) ) then
    print*, 'salt free fluid'
    rho_0 = 0.0_dp
  else
    rho_0 = 1.0_dp / (4.0_dp*pi*bjl*lambda_D**2)
  end if
  rho_ch = 1.0_dp / (4.0_dp*pi*bjl*lambda_D**2) ! DONT UNDERSTAND AWFULL THING WHAT IS THE F* DIFFERENCE BETWEEN RHO_0 AND CH.RHO_0 IN C CODE?

  ! count of solid, fluid and interfacial nodes
  count_solid = count(node%nature==solid)
  count_solid_int = count( node%nature == solid .and. node%isInterfacial )
  print*,'nb of solid interfacial nodes ',count_solid_int
  count_fluid = count(node%nature==fluid)
  print*,'nb of fluid nodes ',count_fluid

  ! read where are distributed the charges
!  call read_charge_distrib



  charge_distrib = input_char("charge_distrib")
  if( charge_distrib(1:3) /= 'int' .and. charge_distrib(1:3)/='sol') stop 'charge_distrib can only be int or sol for now'

  ! distribute charge, depending on where user asked
  if( charge_distrib(1:3) == 'sol') then
    if(count_solid/=0) then
      sigma_solid = sigma / count_solid ! charges distributed in all solid nodes
    else
      sigma_solid = 0
    end if

  else if( charge_distrib(1:3) == 'int') then
    if(count_solid_int/=0) then
      sigma_solid = sigma / count_solid_int ! charges distributed in interfacial solid nodes only
    else
      sigma_solid = 0
    end if
  end if

  if(count_fluid/=0) then
    sigma_fluid = -1.0_dp * sigma / count_fluid ! charges distributed in all fluid nodes
  else
    sigma_fluid = 0
  end if

  in_c_plus_solid  = +0.5_dp*sigma_solid;
  in_c_minus_solid = -0.5_dp*sigma_solid;

  if( sigma_fluid > 0.0_dp ) then
    in_c_plus_fluid  = 0.5_dp*rho_0 + sigma_fluid
    in_c_minus_fluid = 0.5_dp*rho_0
  else
    in_c_plus_fluid  = 0.5_dp*rho_0
    in_c_minus_fluid = 0.5_dp*rho_0 - sigma_fluid
  end if

  if(abs(sigma)>epsilon(1._dp)) then
    if(charge_distrib(1:3)=='int') then
      print*,'The total charge is set ONLY on the solid nodes at the interface (',count_solid_int,'/',count_solid,')'
      print*,'Internal sites (at the interface) = ',count_solid_int,', charge per link =',sigma_solid
      print*,'External sites = ',count_fluid,' charge per link =',sigma_fluid
    else if(charge_distrib(1:3)=='sol') then
      print*,'Internal sites =',count_solid,' charge per link =',sigma_solid
      print*,'External sites =',count_fluid,' charge per link =',sigma_fluid
      stop 'only surface charge is implemented for now in charge_init.f90'
    else
      stop 'pb in charges_init.f90'
    end if

    print*,'Salt concentration ',0.5*rho_0
    print*,'Init density values :'
    print*,'p_solid =',in_c_plus_solid
    print*,'p_fluid =',in_c_plus_fluid
    print*,'m_solid =',in_c_minus_solid
    print*,'m_fluid =',in_c_minus_fluid
    print*,'*********************************************************************'

    print*,'ATTENTION ONLY SURFACE CHARGE IS OK FOR NOW'
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

  anormf0 = 4.0_dp*pi*bjl*kBT/2.0_dp *sum(c_plus+c_minus)

  ! TODO call charge_test

  ! read diffusion coefficients of solute + and solute -
  d_plus = input_dp('D_plus')
  if( D_plus < 0.0_dp ) stop 'D_plus <0. critical.'
  d_minus = input_dp('D_minus')
  if( D_minus < 0.0_dp ) stop 'D_minus <0. critical.'
end subroutine charges_init

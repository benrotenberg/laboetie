PROGRAM main

  use precision_kinds, only: dp
  use mod_time, only: tick, tock
  use io, ONLY: print_tail
  use module_input, only: getinput
  use system, only: supercell
  use module_transient_regime, only: transient_regime
  use module_tracers, only: drop_tracers

  implicit none
  
  integer :: t, lx, ly, lz
  character(8)  :: date
  character(10) :: time
  real(dp), allocatable, dimension(:,:,:) :: jx, jy, jz

  call tick(t) ! this sets up the wall clock time 0 of the simulation

  !
  ! init the system: read input, geometry, composition, periodicity, ...
  !
  call init_simu
  
  call charges_init ! init charge distribution

  !
  ! One solves coupled Poisson and Nernst-Planck equations without solvent flux nor external forces.
  ! i.e., the solvent does not move. This leads to the Poisson-Boltzmann distribution.
  ! To solve the electrostatics equations, we use the Successive Over Relaxation method.
  ! To solve Nernst-Planck, we use Link-Flux without advection.
  !
  if (getinput%log("RestartPNP", defaultvalue=.TRUE.)) call poisson_nernst_planck

  lx = supercell%geometry%dimensions%indiceMax(1)
  ly = supercell%geometry%dimensions%indiceMax(2)
  lz = supercell%geometry%dimensions%indiceMax(3)
  allocate( jx(lx, ly, lz), source=0._dp)
  allocate( jy(lx, ly, lz), source=0._dp)
  allocate( jz(lx, ly, lz), source=0._dp)

  call transient_regime(jx, jy, jz)
  call drop_tracers( jx, jy, jz)
  call print_tail

  print*,"Execution time:", real(tock(t),4),"sec"
END PROGRAM main

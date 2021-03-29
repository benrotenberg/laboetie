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

  !
  ! Setup the wall clock time 0 of the simulation
  !
  call tick(t) 

  !
  ! init the system: read input, geometry, composition, periodicity, ...
  !
  call init_simu
 
  !
  ! initialize the charge distribution
  ! 
  call charges_init 

  !
  ! Here we solve the PNP equations without flow / external perturbation,
  ! i.e., the solvent does not move. This should converge toward the Poisson-Boltzmann distribution.
  !    The Poisson equation is solved by the Successive Over Relaxation (SOR) method.
  !    For the Nernst-Planck equations (diffusion/migration of ions), we use Link-Flux without advection.
  !
  ! This step can be bypassed using a restart from a previous simuation.
  !    It is only performed if "RestartPNP = T" in lb.in
  !
  if (getinput%log("RestartPNP", defaultvalue=.TRUE.)) call poisson_nernst_planck

  !
  ! Isn't there a way to put this somewhere else than the main?
  !
  lx = supercell%geometry%dimensions%indiceMax(1)
  ly = supercell%geometry%dimensions%indiceMax(2)
  lz = supercell%geometry%dimensions%indiceMax(3)
  allocate( jx(lx, ly, lz), source=0._dp)
  allocate( jy(lx, ly, lz), source=0._dp)
  allocate( jz(lx, ly, lz), source=0._dp)

  !
  ! Here we do the actual LBE:
  !     Starting from the PB equilibrium,
  !     we can apply an external force (mimicking a pressure gradient) or an
  !     electric field and iterate PNP coupled with Navier-Stokes until convergence
  !
  !     Note that using a restart we may start from a PB equilibrium for a
  !     different situation (e.g. voltage for capacitors) and investigate the
  !     transient regime to reach a new steady-state
  !
  call transient_regime(jx, jy, jz)

  !
  ! Once the steady state is reached, we perform Moment Propagation for (charged) tracers
  !
  call drop_tracers( jx, jy, jz)

  !
  ! Finalize simulation
  !
  call print_tail

  print*,"Execution time:", real(tock(t),4),"sec"

END PROGRAM main

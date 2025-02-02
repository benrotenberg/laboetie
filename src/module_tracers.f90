! Here tracers are dropped in the fluid. They may be charged or not. They evolve in the
! equilibrated fluid and its solutes. They do not change the potential even if they
! have a charge. The idea is to follow them in order to get the velocity auto correlation
! function, while making them not to change the equilibrium properties of the system.
module module_tracers
  implicit none
  contains

  !
  ! This performs the Moment Propagation
  !     by calling routines in module moment_propagation
  ! 
  SUBROUTINE drop_tracers( solventCurrentx, solventCurrenty, solventCurrentz)

    use precision_kinds, only: dp
    USE system, only: n, elec_slope
    USE moment_propagation, only: init, propagate, deallocate_propagated_quantity
    use module_input, ONLY: getinput
    use myallocations, only: allocateReal3D

    IMPLICIT NONE
    real(dp), intent(in), dimension(:,:,:) :: solventCurrentx, solventCurrenty, solventCurrentz
    integer :: it, maximum_moment_propagation_steps
    logical :: is_converged
    real(dp), allocatable :: solventDensity(:,:,:)

    maximum_moment_propagation_steps = getinput%int("maximum_moment_propagation_steps", 0) ! negative value means make it converge
    IF( maximum_moment_propagation_steps == 0) RETURN

    PRINT*
    PRINT*,'Moment propagation'
    PRINT*,'=================='
    PRINT*,'       step           VACF(x)                   VACF(y)                   VACF(z)'
    PRINT*,'       ----------------------------------------------------------------------------------'

    allocate( solventDensity(size(n,1),size(n,2),size(n,3)) , source=sum(n,4) )
    print*, 'STEP1'

    CALL update_tracer_population( solventDensity, solventCurrentx, solventCurrenty, solventCurrentz ) ! include elec_slope in population n
    PRINT*, 'STEP2'
    elec_slope = 0.0_dp ! turn the electric field off for the rest of mom_prop (included in n)

    ! add electrostatic potential computed by the SOR routine an external contribution
    ! elec_pot(singlx,ly,lz, ch, phi, phi2, t, t_equil);
    ! call elec_pot
    CALL init( solventDensity) ! init moment propagation
    print*, '+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+'
    print*, 'init called'
    print*, '+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+'

    !
    ! Propagation in time
    !
    if( maximum_moment_propagation_steps < 0) maximum_moment_propagation_steps = HUGE(1)
    DO it= 1, maximum_moment_propagation_steps
      !  it=0
      !  do while (not_yet_converged(it))
      !   call elec_pot
      CALL propagate (it,is_converged, solventDensity) ! propagate the propagated quantity
      !    if( modulo(it,50000)==0 ) print_vacf
      !    it = it + 1
      IF( is_converged ) exit
    END DO

    PRINT*,

    !  CALL integrate_vacf ! compute the integral over time of vacf in each direction
    !  CALL print_vacf ! print vacf to file vacf.dat
    CALL deallocate_propagated_quantity

!contains
  end SUBROUTINE drop_tracers


  !
  ! Modify populations to account for external force and electric field
  !     this latter contribution depends on the tracer charge
  !
  SUBROUTINE update_tracer_population( solventDensity, solventCurrentx, solventCurrenty, solventCurrentz )

    use precision_kinds, only: dp, i2b
    use system, only: f_ext, fluid, elec_slope, supercell, lbm, x, y, z, node, n
    use module_collision, only: check_population
    use module_input, only: input_dp, input_dp3

    implicit none

    real(dp), intent(in) :: solventDensity(:,:,:), solventCurrentx(:,:,:), solventCurrenty(:,:,:), solventCurrentz(:,:,:)
    integer(i2b) :: l, ll, lu, lx, ly, lz, i, j, k
    type tracer
      real(dp) :: D ! diffusion coefficient
      real(dp) :: q ! charge
    end type
    type(tracer) :: tr

    ll = lbm%lmin
    lu = lbm%lmax
    lx = supercell%geometry%dimensions%indiceMax(x)
    ly = supercell%geometry%dimensions%indiceMax(y)
    lz = supercell%geometry%dimensions%indiceMax(z)

    deallocate(n) ! the purpose here is to prepare n with the good index arrangement to fit correctly in memory with inner loop over l in moment propagation.
    allocate( n(ll:lu,lx,ly,lz) ,source=0._dp)

    tr%D = input_dp('tracer_Db',0._dp)
    if (tr%D<=epsilon(1._dp)) ERROR STOP 'The diffusion coefficient (tracer_Db in input file) is invalid'

    tr%q = input_dp('tracer_z', 0._dp)
    f_ext(:) = input_dp3('f_ext', [0._dp,0._dp,0._dp] )

    !
    ! Apply force on all fluid nodes and update populations
    !
    do concurrent (l=ll:lu, i=1:lx, j=1:ly, k=1:lz)
      if( node(i,j,k)%nature==fluid ) then
        n(l,i,j,k) = lbm%vel(l)%a0 *solventdensity(i,j,k) + lbm%vel(l)%a1 &
          *sum( lbm%vel(l)%coo(:)*( [solventCurrentx(i,j,k), solventCurrenty(i,j,k), solventCurrentz(i,j,k)] + f_ext(:) - solventDensity(i,j,k)*tr%q*tr%D*elec_slope(:) ) ) 
      else
        n(l,i,j,k) = lbm%vel(l)%a0 *solventdensity(i,j,k) + lbm%vel(l)%a1 &
          *( lbm%vel(l)%coo(x)*solventCurrentx(i,j,k) + lbm%vel(l)%coo(y)*solventCurrenty(i,j,k) + lbm%vel(l)%coo(z)*solventCurrentz(i,j,k) )
      end if
    end do

    !
    ! do concurrent( l= lbm%lmin: lbm%lmax )
    !   where( node%nature ==fluid )
    !     n(l,:,:,:) = lbm%vel(l)%a0*solventDensity &
    !                + lbm%vel(l)%a1*(&
    !        lbm%vel(l)%coo(x)*(node%solventFlux(x) + f_ext(x) - solventDensity*tr%q *tr%D *elec_slope(x)) &
    !      + lbm%vel(l)%coo(y)*(node%solventFlux(y) + f_ext(y) - solventDensity*tr%q *tr%D *elec_slope(y)) &
    !      + lbm%vel(l)%coo(z)*(node%solventFlux(z) + f_ext(z) - solventDensity*tr%q *tr%D *elec_slope(z)) )
    !   elsewhere
    !     n(l,:,:,:) = lbm%vel(l)%a0*node%solventDensity + lbm%vel(l)%a1*( &
    !                       lbm%vel(l)%coo(x)*node%solventFlux(x) + &
    !                       lbm%vel(l)%coo(y)*node%solventFlux(y) + &
    !                       lbm%vel(l)%coo(z)*node%solventFlux(z) )
    !   end where
    ! end do

  end SUBROUTINE update_tracer_population

end module module_tracers

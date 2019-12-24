! This is the first part of the LaBo code. Here, no forces are applied on solutes,
! and no flux is allowed. The objective here is to find the initial equilibrium
! distribution of the charges.

SUBROUTINE poisson_nernst_planck

  USE precision_kinds, ONLY: dp
  USE system, ONLY: D_equil, tot_sol_charge, time, elec_slope, lncb_slope, node, c_plus, c_minus, supercell, phi, & 
                    fluid, bjl, lambda_D, solid, LaplacianOfPhi, el_curr_x, el_curr_y, el_curr_z
  use module_input, only: getinput
  USE io, ONLY: print_everything_related_to_charge_equil
  use constants, only: x,y,z, pi
  use myallocations
  use module_convergence


  IMPLICIT NONE

  INTEGER :: timestep, timestepmax, i, j, k, MyGeometry, lx, ly, lz, capacitor,&
             slit_sol, m, ENNE
  LOGICAL :: is_converged = .FALSE.
  !real(dp), allocatable, dimension(:,:,:) :: phiTMP ! Ade
  real(dp) :: target_error_charge, max_error
  real(dp) :: SF, Alpha, FACT1, FACT2, FLAT, Somma1, Somma2, Somma3, Somma4, Somma5, Somma6
  integer :: print_freq
  character(len=200) :: ifile

  open(314, file='output/c_plus_alongZ.dat')
  open(315, file='output/c_minus_alongZ.dat')
  open(316, file='output/c_plus_beforeLoop.dat')
  open(280, file='output/PHI_PNP1.dat')
  open(281, file='output/PHI_PNP2.dat')
  ! If tot_sol_charge is 0, neutral system, and thus no need to go continue here
  tot_sol_charge = getinput%dp('tot_sol_charge',0._dp)
  !IF ( ABS(tot_sol_charge) <= EPSILON(1._dp) .and. ABS(lambda_D) <= EPSILON(1._dp)) RETURN ! Ade: 22/05/17 I added the and

  PRINT*
  PRINT*,'Poisson + Nernst-Planck'
  PRINT*,'======================='

  ! read the number of iterations one does for the first step of equilibration (D_equil)
  CALL get_timestepmax( timestepmax )
  D_equil = timestepmax ! Ade: 30/03/2017 D_equil was never given a value. Thus it was always
                        ! equal to zero. This explains why the profiles would never iterate. 

  !
  ! BR: We should probably get rid of D_equil, not so useful anymore
  !
  IF (D_equil<0) STOP "D_equil in input file should be >= 0"
  !print*, 'D_equil = ', D_equil


  ! Ade: modification 5/04/2018
  el_curr_x = 0._dp
  el_curr_y = 0._dp
  el_curr_z = 0._dp

  ! read electrostatic related stuff
  lncb_slope = getinput%dp3("lncb_slope")
  elec_slope = getinput%dp3("elec_slope")
  
  ! BR: added target for convergence of charge distribution
  target_error_charge = target_error%target_error_charge

  ! Ade : modification 20/03/2017
  lx = supercell%geometry%dimensions%indiceMax(x)
  ly = supercell%geometry%dimensions%indiceMax(y)
  lz = supercell%geometry%dimensions%indiceMax(z)
  MyGeometry = getinput%int('geometryLabel',-1)
  
  !IF( .NOT. ALLOCATED(phiTMP) ) CALL allocateReal3D(phiTMP)
  IF( .NOT. ALLOCATED(phi) ) CALL allocateReal3D(phi)


  !phi = 0.0_dp ! Ade : initialise phi
  capacitor = getinput%int('capacitor',0) ! 1 = true 0 = false
  slit_sol = getinput%int('slit_sol',0) ! 1= true 0 = false
  ! Ade : m is used for the initial solution input for c_plus        
  m = getinput%int("SolidSheetsEachSideOfSlit", defaultvalue=1, assert=">0")

  if(capacitor.NE.1 .and. slit_sol==1) then

     if( MyGeometry == 1 ) then   ! slit pore geometry
        Alpha = getinput%dp('Alpha',0._dp)  ! Attention!!!!!! This is dangerous. We should probably  
                                            ! compute its value in another subroutine
        ENNE = lz-(2*m) ! Nb of fluid nodes
        if(mod(lz,2) == 0) then
           SF = real(lz)/2
        else
           SF = real(lz)/2 + 0.5
        endif

        !TODO : BOUGER LE TRUC CI DESSOUS POUR QUE LA SOLUTION ANALYTIQUE SOIT APPELEE SEULVEMENT SI ON VEUT ET A LEXTERIEUR
        if( lambda_d == 0.0) then ! no salt added
            do i = 1, lx      
                do j = 1, ly  
                  do k = m+1, lz-m ! first and last nodes are solid
                    ! Ade : below is the analytical solution for the potential (called phi) for a slit geometry.
                    !phiTMP(i,j,k) = 2.0*log(cos(Alpha*(k-SF))/(cos(Alpha*(real(ENNE)*0.5)) ))
                    phi(i,j,k) = 2.0*log( cos(Alpha*(k-SF)) / cos(Alpha*(real(ENNE)*0.5)) )
                  end do
                end do
            end do
            FLAT =  2.0*log(cos(Alpha*(m-SF))/(cos(Alpha*(real(ENNE)*0.5)) ))
        else ! salt added
            FACT1 = 2.0*pi*tot_sol_charge*bjl*lambda_D/(lx*ly)
            FACT2 = ( sinh(real(ENNE)/(2.0*lambda_D)) )
            do i = 1, lx      
                do j = 1, ly  
                  do k = m+1, lz-m ! first and last node are solid
                        !phiTMP(i,j,k) = FACT1 * cosh( (k-SF)/lambda_D )/( FACT2 )
                    phi(i,j,k) = FACT1 * cosh( (k-SF)/lambda_D )/( FACT2 )
                  end do
                end do
            end do
            FLAT =  FACT1 * cosh( (m-SF)/lambda_D )/( FACT2 )
        endif
        !phi = phiTMP 
        where(node%nature==solid)
            phi = FLAT
        end where

     endif ! end if( MyGeometry == 1 )

  end if ! end if(capacitor.NE.1 .and. slit_sol==1) then


  DO k=1,supercell%geometry%dimensions%indiceMax(z)
    write(280,*) k, SUM(phi(:,:,k)) !, sum(phiTMP(:,:,k))
  ENDDO
  close(280)
 ! Ade : end modification 20/03/2017
  
  DO k=supercell%geometry%dimensions%indiceMin(z), supercell%geometry%dimensions%indiceMax(z)
      WRITE(316,*) k, c_plus(:,:,k)
  ENDDO
  close(316)

  !##########################################################################
  !#  Main loop to converge towards Poisson-Boltzmann equilibrium 
  !##########################################################################

  ! iterate until charges are equilibrated
  ! Ade : modifications on 23/03/2017
  time = 0
  print_freq = 1
  DO WHILE ((.NOT. is_converged) .AND. (time<=timestepmax))

      CALL backup_phi_c_plus_c_minus  ! backup potential and solute concentrations from last step

      CALL sor                        ! compute phi from concentrations with the Successive Overrelation Routine (SOR)

      CALL just_eq_smolu              ! solve Smoluchowski (diffusion + electrostatic part)

      ! monitor evolution of phi, c_plus, c_minus w.r.t. the backup at the beginning of the iteration
      ! this is done only every 10 loops in order not to waste too much time. arbitrary number.
      IF ((time/= timestepmax .and. modulo(time,10)==0) .or. tot_sol_charge==0.0_dp ) THEN 
          CALL check_charge_distribution_equilibrium (time, target_error_charge, max_error, is_converged)
      END IF

      time = time + 1

      ! inform of current state with regular outputs 
      if( modulo(time, print_freq) == 0) then
            if( time==10    )   print_freq = 10
            if( time==100   )   print_freq = 100
            if( time==1000  )   print_freq = 1000

            WRITE(*,"(4X,I8,A,E12.6,A,E12.6,A)") time, " crit. on charge ", max_error, " ( target ", target_error_charge, " )"
        end if


  END DO
  ! end of main loop to converge towards Poisson-Boltzmann equilibrium


  !##########################################################################
  !# Postprocessing once Poisson-Boltzmann equilibrium is reached... or not
  !##########################################################################

  !
  ! First, backup the arrays phi, c_plus and c_minus for restarts
  !
  open(731, file='output/phi.bin'    , form="unformatted"); write(731) phi    ; close(731)
  open(731, file='output/c_plus.bin' , form="unformatted"); write(731) c_plus ; close(731)
  open(731, file='output/c_minus.bin', form="unformatted"); write(731) c_minus; close(731)

  !
  ! Some output (to be cleaned?)
  !
  if (MyGeometry==2 .or. MyGeometry==15) then ! Cylindrical geometry
      DO j=1,ly
          WRITE(281,*) j, phi(j,ly/2,lz/2)
          WRITE(314,*) j, c_plus(j,ly/2,lz/2)
          WRITE(315,*) j, c_minus(j,ly/2,lz/2)
      END DO
  else
      DO k=supercell%geometry%dimensions%indiceMin(z), supercell%geometry%dimensions%indiceMax(z)
          write(281,*) k, SUM(phi(:,:,k)), phi(1,1,k)
          WRITE(314,*) k, SUM(c_plus(:,:,k))
          WRITE(315,*) k, SUM(c_minus(:,:,k))
      END DO
  end if

  if( tot_sol_charge == 0.0_dp ) return
  call print_everything_related_to_charge_equil

  IF (is_converged) THEN
      PRINT*,'PNP converged after',time,' steps'
  ELSE
      STOP 'Poisson-Nernst-Planck did not converge to Poisson-Boltzmann equilibrium'
  END IF

  close(314)
  close(315)
  close(281)


  !#########################################################################
  CONTAINS ! in a subroutine?

      SUBROUTINE get_timestepmax( a )
          IMPLICIT NONE
          INTEGER :: a
          a = getinput%int("timestepmax_for_PoissonNernstPlanck",1000000)
          IF (a<=0) STOP "timestepmax_for_PoissonNernstPlanck must be >=1"
      END SUBROUTINE



end subroutine poisson_nernst_planck

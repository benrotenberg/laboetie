! This is the first part of the LaBo code. Here, no forces are applied on solutes,
! and no flux is allowed. The objective here is to find the initial equilibrium
! distribution of the charges.

SUBROUTINE poisson_nernst_planck

  USE precision_kinds, ONLY: dp
  USE system, ONLY: D_equil, tot_sol_charge, time, elec_slope, lncb_slope, node, c_plus, c_minus, supercell, phi, & ! ADE: I added phi, c_plus and c_minus
                    fluid, bjl, lambda_D, solid, LaplacianOfPhi, el_curr_x, el_curr_y, el_curr_z
  use module_input, only: getinput
  USE io, ONLY: print_everything_related_to_charge_equil
  use constants, only: x,y,z, pi
  use myallocations


  IMPLICIT NONE

  INTEGER :: timestep, timestepmax, i, j, k, geometrie, lx, ly, lz, capacitor,&
             slit_sol, m, ENNE
  LOGICAL :: is_converged = .FALSE.
  real(dp), allocatable, dimension(:,:,:) :: phiTMP ! Ade
  real(dp) :: SF, Alpha, FACT1, FACT2, FLAT, Somma1, Somma2, Somma3, Somma4, Somma5, Somma6
  character*200 :: ifile

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

  ! read the number of iterations one does for the first step of equilibration (D_iter)
  CALL get_timestepmax (timestepmax)
  D_equil = timestepmax ! Ade: 30/03/2017 D_equil was never given a value. Thus it was always
                        ! equal to zero. This explains why the profiles would never iterate. 

  IF (D_equil<0) STOP "D_equil in input file should be >= 0"
  print*, 'D_equil = ', D_equil


  ! Ade: modification 5/04/2018
  el_curr_x = 0._dp
  el_curr_y = 0._dp
  el_curr_z = 0._dp

  ! read electrostatic related stuff
  lncb_slope = getinput%dp3("lncb_slope")
  elec_slope = getinput%dp3("elec_slope")
  

  ! Ade : modification 20/03/2017
  lx = supercell%geometry%dimensions%indiceMax(x)
  ly = supercell%geometry%dimensions%indiceMax(y)
  lz = supercell%geometry%dimensions%indiceMax(z)
  geometrie = getinput%int('geometryLabel',-1)
  ! Ade : attention phiTMP here is not the same as in SOR
  IF( .NOT. ALLOCATED(phiTMP) ) CALL allocateReal3D(phiTMP)
  IF( .NOT. ALLOCATED(phi) ) CALL allocateReal3D(phi)
  !phi = 0.0_dp ! Ade : initialise phi
  capacitor = getinput%int('capacitor',0) ! 1 = true 0 = false

  slit_sol = getinput%int('slit_sol',0) ! 1= true 0 = false
  ! Ade : m is used for the initial solution input for c_plus        
    m = getinput%int("SolidSheetsEachSideOfSlit", defaultvalue=1, assert=">0")
    if(capacitor.NE.1 .and. slit_sol==1) then
    if( geometrie == 1 ) then   ! slit pore geometry
        Alpha = getinput%dp('Alpha',0._dp)  ! Attention!!!!!! This is dangerous. We should probably  
                                            ! compute its value in another subroutine
        ENNE = lz-(2*m) ! Nb of fluid nodes
        if(mod(lz,2) == 0) then
           SF = real(lz)/2
        else
           SF = real(lz)/2 + 0.5
        endif
!TODO : BOUGER LE TRUC CI DESSOUS POUR QUE LA SOLUTION ANALYTIQUE SOIT APPELEE SEULVEMENT SI ON VEUT ET A LEXTERIEUR
        ! Below : cas sans sel
        if( lambda_d == 0.0) then ! no salt added
            do i = 1, lx      
                do j = 1, ly  
                  do k = m+1, lz-m ! first and last node are solid
                    ! Ade : below is the analytical solution for the potential (called phi) for a slit geometry.
                    phiTMP(i,j,k) = 2.0*log(cos(Alpha*(k-SF))/(cos(Alpha*(real(ENNE)*0.5)) ))! Potential phi analytical solution for a slit
                  end do
                end do
            end do
            FLAT =  2.0*log(cos(Alpha*(m-SF))/(cos(Alpha*(real(ENNE)*0.5)) ))
        ! below : cas avec sel
        else ! salt added
            FACT1 = 2.0*pi*tot_sol_charge*bjl*lambda_D/(lx*ly)
            FACT2 = ( sinh(real(ENNE)/(2.0*lambda_D)) )
            do i = 1, lx      
                do j = 1, ly  
                  do k = m+1, lz-m ! first and last node are solid
                        phiTMP(i,j,k) = FACT1 * cosh( (k-SF)/lambda_D )/( FACT2 )
                  end do
                end do
            end do
            FLAT =  FACT1 * cosh( (m-SF)/lambda_D )/( FACT2 )
        endif
        phi = phiTMP 
        where(node%nature==solid)
            phi = FLAT
        end where
    endif
  end if

    DO k=1,supercell%geometry%dimensions%indiceMax(z)
      write(280,*) k, SUM(phi(:,:,k)), sum(phiTMP(:,:,k))
    ENDDO
    close(280)
  ! Ade : end modification 20/03/2017
    
    DO k=supercell%geometry%dimensions%indiceMin(z), supercell%geometry%dimensions%indiceMax(z)
        WRITE(316,*) k, c_plus(:,:,k)
    ENDDO
    close(316)

  ! iterate until charges are equilibrated
  ! Ade : modifications on 23/03/2017
    time = 0
    DO WHILE ((.NOT. is_converged) .AND. (time<=timestepmax))
        CALL backup_phi_c_plus_c_minus ! backup potential and solute concentrations from last step
        CALL sor ! TODO    ! compute phi with the Successive Overrelation Routine (SOR)
        CALL just_eq_smolu ! solve smoluchowski (diffusion + electrostatic part) ie not a full smolu
        ! monitor evolution of phi, c_plus, c_minus w.r.t. the backup at the beginning of the iteration
        ! this is done only every 10 loops in order not to waste too much time. arbitrary number.

        IF ((time/= timestepmax .and. modulo(time,10)==0) .or. tot_sol_charge==0.0_dp ) THEN ! Ade : removed minus sign in front of timestepmax
            CALL check_charge_distribution_equilibrium (time, is_converged)
        END IF
        time = time + 1 ! Ade : missing incrementing number
    END DO

    ! ---------------------------------------------------
    ! POSTPROCESSING AT EQUILIBRIUM 
    ! ---------------------------------------------------

    !
    ! First, backup the arrays phi, c_plus and c_minus for restarts
    !
    open(731, file='output/phi.bin'    , form="unformatted"); write(731) phi    ; close(731)
    open(731, file='output/c_plus.bin' , form="unformatted"); write(731) c_plus ; close(731)
    open(731, file='output/c_minus.bin', form="unformatted"); write(731) c_minus; close(731)

    if (geometrie==2 .or. geometrie==15) then ! Cylindrical geometry
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
    !call charge_test ! check charge conservation ! TODO rename to call check_charge_conservation ! check charge conservation every 1000 steps (arbitrary number)
    if( tot_sol_charge == 0.0_dp ) return
    call print_everything_related_to_charge_equil

    IF (is_converged) THEN
        PRINT*,'Convergence found at step ',time,' after',time,' steps'
    ELSE
        STOP 'Equilibrium distribution of salts not found'
    END IF

    close(314)
    close(315)
    close(281)


    CONTAINS

        SUBROUTINE get_timestepmax (a)
            IMPLICIT NONE
            INTEGER :: a
            a = getinput%int("timestepmax_for_PoissonNernstPlanck",1000000)
            IF (a<=0) STOP "timestepmax_for_PoissonNernstPlanck must be >=1"
        END SUBROUTINE



end subroutine poisson_nernst_planck

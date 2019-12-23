! Here we define the criteria used to define convergence of the
! iterative processes: SOR, PNP->PB, main loop, MP

module module_convergence

        use precision_kinds
        implicit none

        type type_convergence
                real(dp) :: target_error_mass_flux      ! for mass flux in module_transient_regime
                real(dp) :: target_error_charge         ! for potential and ion concentrations in module_transient_regime
                real(dp) :: target_error_sor            ! for potential in SOR
        end type

        type (type_convergence), public :: target_error


contains

        subroutine init_convergence_targets

          use module_input, only: getinput
          implicit none        

          target_error%target_error_mass_flux = getinput%dp("target_error_mass_flux", 1.D-10 )

          target_error%target_error_charge    = getinput%dp("target_error_charge",    1.D-12 )

          target_error%target_error_sor       = getinput%dp("target_error_sor",       1.D-8  )

        end subroutine

end module module_convergence

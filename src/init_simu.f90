SUBROUTINE init_simu

    USE precision_kinds, only: dp
    USE mod_lbmodel, only: init_everything_related_to_lb_model => initialize, lbm
    USE system, only: node, n, solid
    USE io, only: print_header, print_input_in_output_folder, inquireNecessaryFilesExistence
    USE myallocations
    use module_input, only: getinput

    IMPLICIT NONE

    REAL(dp) :: initialSolventDensity
    integer :: l

    CALL print_header
    CALL inquireNecessaryFilesExistence  ! check that input, output folder and file ./lb.in exist
    CALL init_everything_related_to_lb_model ! init everything related to D3Q15 or D3Q19 etc ie LB models
    CALL supercell_definition ! prepare supercell geometry
    CALL scheduler ! tmom, tmax! schedule simulation

    !
    ! Initialize solvent populations, noted n(x,y,z,l) in laboetie.
    ! It should be 0 in solid nodes and 1 by default in the liquid phase.
    !
    if( .not. allocated(n) ) call allocatereal4D(n)
    initialSolventDensity = getinput%dp("initialSolventDensity", defaultvalue=1._dp, assert=">0")
    do l = lbm%lmin, lbm%lmax
        where( node%nature /= solid )
            n(:,:,:,l) = initialSolventDensity * lbm%vel(l)%a0
        elsewhere
            n(:,:,:,l) = 0._dp
        end where
    end do

CONTAINS
    !
    !
    !
SUBROUTINE scheduler
        use system, only: tmom, tmax, D_iter, time
        use module_input, only: getinput
        ! 4 times are important :
        ! - 0 at which simulation starts
        ! - tmom at which we're looking at tracer moment propagation
        ! - tmax at which simulation stops
        ! init to non-physical value catched later in order to be sure they are modified
        time = 0
        D_iter = getinput%int('D_iter',1) ! ADE : the default value should be 1 (check)
        tmax = getinput%int('tmax',-1)
        tmom = getinput%int('tmom',-1)
end subroutine
    !
    !
    !
end subroutine

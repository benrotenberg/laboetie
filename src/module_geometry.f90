module geometry

    use precision_kinds ! all precision kinds defined in the dedicated module
    use system, only: fluid, solid, supercell, node
    use constants, only: x, y, z
    use module_input, only: getinput

    implicit none

contains

SUBROUTINE CONSTRUCT_TUBE_WITH_VARYING_DIAMETER
! construct the system described in Makhnovskii et al., Chem. Phys. 367, 110 (2010)
! 4 effective variables: length of big pore lw, narrow pore ln, radius of each, R and a, respectively.
! the tube is along x, with infinity sym along z (ie lz = 1 ok)
  real(dp) :: R, a, ln, lw, lw_o_lx, a_o_R
  real(dp), dimension(2) :: r0, rjk
  integer(i2b) :: i, j, k
  logical :: is_here
            integer(i2b) :: lx, ly, lz
            lx = supercell%geometry%dimensions%indiceMax(x)
            ly = supercell%geometry%dimensions%indiceMax(y)
            lz = supercell%geometry%dimensions%indiceMax(z)
  if(ly/=lz) stop 'ly should be equal to lz cause it is a tube'
  node%nature = solid
  inquire(file='tube.in',exist=is_here)
  if(.not.is_here) stop 'tube.in cannot be read'
  open(unit=12,file='tube.in')
  read(12,*) lw_o_lx, a_o_R
  R = real(ly-1)/2.0
  a = a_o_R*R
  r0 = [real(ly+1)/2.,real(lz+1)/2.]
  do concurrent(i=1:lx, j=1:ly, k=1:lz)
    rjk = [real(j),real(k)]
    if(i>lw_o_lx*lx) then
      if( norm2( rjk-r0 ) < a ) then
        node(i,j,k)%nature = fluid
      else
        node(i,j,k)%nature = solid
      end if
    else if(i<=lw_o_lx*lx) then
      if( norm2( rjk-r0 ) < R ) then
        node(i,j,k)%nature = fluid
      else
        node(i,j,k)%nature = solid
      end if
    end if
  end do
END SUBROUTINE CONSTRUCT_TUBE_WITH_VARYING_DIAMETER






SUBROUTINE CONSTRUCT_PLANES_WITH_VARIOUS_RADIUS_2D
! construct the system described in Makhnovskii et al., Chem. Phys. 367, 110 (2010)
! 4 effective variables: length of big pore lw, narrow pore ln, radius of each, R and a, respectively.
! the tube is along x, with infinity sym along z (ie lz = 1 ok)
  real(dp) :: R, a, ln, lw
  integer(i2b) :: i
            integer(i2b) :: lx, ly, lz
            lx = supercell%geometry%dimensions%indiceMax(x)
            ly = supercell%geometry%dimensions%indiceMax(y)
            lz = supercell%geometry%dimensions%indiceMax(z)
!  R = real(ly-2,dp)/2._dp
!  a = R/2._dp
!  lw = real(lx,dp)/2._dp
!  ln = lw
  node%nature = fluid
  node(:,1,:)%nature = solid
  node(:,ly,:)%nature = solid
  do i=1,lx
    if( i>lx/2 ) then
      node(i,2:1+(ly-2)/4,:)%nature = solid
      node(i,ly-(ly-2)/4:ly-1,:)%nature = solid
    end if
  end do
END SUBROUTINE CONSTRUCT_PLANES_WITH_VARIOUS_RADIUS_2D


SUBROUTINE CONSTRUCT_SINUSOIDAL_WALLS_2D
  ! tunnel along lx, ly is the width of the tunnel
  real(dp) :: xi, yi
  real(dp) :: a, b
  real(dp), parameter :: pi = acos(-1.0_dp)
  integer(i2b) :: i, j, ninty, mirror, midly
              integer(i2b) :: lx, ly, lz
            lx = supercell%geometry%dimensions%indiceMax(x)
            ly = supercell%geometry%dimensions%indiceMax(y)
            lz = supercell%geometry%dimensions%indiceMax(z)
  ! ly should be odd, so that the middle of ly is on a node (w(x)=0)
  if( mod(ly,2)==0 ) stop 'ly should be odd'
  midly = (ly+1)/2
  b = 2./3.*midly
  a = midly - b - 1
  if( a<=0 .or. b<=0 ) stop 'pb in def geometry function'
  do i = 1, lx
    xi = real(i-1,dp)
    yi = a*sin(2._dp*pi*xi/Lx) + b + midly
    ninty = nint(yi)
    mirror = midly
    do j = midly, ly
      if( j < ninty ) then
        node(i,j,:)%nature = fluid
      else
        node(i,j,:)%nature = solid
      end if
      node(i, mirror, :)%nature = node(i, j, :)%nature
      mirror = mirror - 1
    end do
  end do
END SUBROUTINE CONSTRUCT_SINUSOIDAL_WALLS_2D


SUBROUTINE construct_slit_Nsheets(n)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    INTEGER :: i, j, k
    i = LBOUND( node%nature, 3)
    j = UBOUND( node%nature, 3)
    node%nature = fluid
    do k = 0, n-1
        node(:,:,i+k)%nature = solid
        node(:,:,j-k)%nature = solid
    end do
    ! node(:,:,i)%nature = solid   ! The lower bound of the third dimension is the first solid sheet
    ! node(:,:,i+1)%nature = solid ! the second sheet, etc.
    ! ...
    ! node(:,:,j)%nature = solid   ! Same for the upper bound
    ! node(:,:,j-1)%nature = solid ! the sheet before the upper bound is solid
    ! ...
END SUBROUTINE




SUBROUTINE CONSTRUCT_SPHERICAL_CAVITY
  implicit none
  real(dp) :: radius ! radius of cylinder
  real(dp), dimension(3) :: rnode, rorigin ! coordinates of each node and center of cylinder in x,y coordinates
  integer(i2b) :: i, j, k ! dummy
              integer(i2b) :: lx, ly, lz
            lx = supercell%geometry%dimensions%indiceMax(x)
            ly = supercell%geometry%dimensions%indiceMax(y)
            lz = supercell%geometry%dimensions%indiceMax(z)
  if( lx /= ly .or. lx/=lz) stop 'wall=8 for a spherical cavity so lx=ly=lz. check input file.'
  if( lx<3 ) stop 'the diameter of the cylinder (lx) should be greater than 3'
  rorigin = [ real(lx+1,dp)/2.0_dp, real(ly+1,dp)/2.0_dp, real(lz+1,dp)/2.0_dp ]
  radius = real(lx-1,dp)/2.0_dp
  do i = 1, lx
    do j = 1, ly
      do k = 1, lz
        rnode = [real(i,dp),real(j,dp),real(k,dp)] - rorigin
        if( norm2(rnode) >= radius ) then ! = radius is important because without it one has exists
          node(i,j,k)%nature = solid
        else
          node(i,j,k)%nature = fluid
        end if
      end do
    end do
  end do
END SUBROUTINE CONSTRUCT_SPHERICAL_CAVITY







SUBROUTINE CONSTRUCT_CC
  implicit none
  integer(i2b) :: i,j,k ! dummy
  integer(i2b) :: lx, ly, lz
            lx = supercell%geometry%dimensions%indiceMax(x)
            ly = supercell%geometry%dimensions%indiceMax(y)
            lz = supercell%geometry%dimensions%indiceMax(z)
  if( lx /= ly .or. lx /= lz ) stop 'with wall = 3, i.e. cfc cell, the supercell should be cubic with lx=ly=lz'
  node%nature = fluid
  do concurrent( i=1:lx, j=1:ly, k=1:lz )
    if( is_in_solid_sphere(i,j,k) ) then
      node(i,j,k)%nature = solid
    else
      node(i,j,k)%nature = fluid
    end if
  end do
  contains
    PURE LOGICAL FUNCTION IS_IN_SOLID_SPHERE(i,j,k)
      implicit none
      integer(i2b), intent(in) :: i,j,k
      real(dp), dimension(9) :: distances
      real(dp), dimension(3) :: r
      r = real([i,j,k],dp)
      distances(1) = norm2( r - real([1,1,1]) ) ! origin is not 0,0,0 but 1,1,1 in our referential 1:lx
      distances(2) = norm2( r - real([lx,1,1]) )
      distances(3) = norm2( r - real([1,ly,1]) )
      distances(4) = norm2( r - real([1,1,lz]) )
      distances(5) = norm2( r - real([lx,ly,1]) )
      distances(6) = norm2( r - real([lx,1,lz]) )
      distances(7) = norm2( r - real([1,ly,lz]) )
      distances(8) = norm2( r - real([lx,ly,lz]) )
      distances(9) = norm2( r - real([lx+1,ly+1,lz+1])/2._dp )
      ! put distances to all of corners of the cube and to its center (=9 points) in an array
      if ( any(distances <= real(lx-1,dp)*sqrt(3.0_dp)/4.0_dp) ) then
        is_in_solid_sphere = .true.
      else
        is_in_solid_sphere = .false.
      end if
    END FUNCTION IS_IN_SOLID_SPHERE
END SUBROUTINE CONSTRUCT_CC







SUBROUTINE CONSTRUCT_CYLINDER
  implicit none
  real(dp) :: radius ! radius of cylinder
  real(dp), dimension(2) :: rnode, rorigin ! coordinates of each node and center of cylinder in x,y coordinates
  integer(i2b) :: i, j, width ! dummy
            integer(i2b) :: lx, ly, lz
            lx = supercell%geometry%dimensions%indiceMax(x)
            ly = supercell%geometry%dimensions%indiceMax(y)
            lz = supercell%geometry%dimensions%indiceMax(z)
  if( lx /= ly) stop 'wall=2 is for cylinders, which should have same lx and ly'
  if( lx<3 ) stop 'the diameter of the cylinder (lx) should be greater than 3'
  rorigin = [ real(lx+1,dp)/2.0_dp, real(ly+1,dp)/2.0_dp ]
  width = getinput%int("width", defaultvalue=2) ! Ade : solid width. This parameter is used to make smaller tubes
                                                ! within our box of simulation. This is usefull as error increases
                                                ! if w/2 = delta_x. Error is small if w/2 = 50*delta_x. Please read
                                                ! Amael Obliger thesis p 46.
  radius = real(lx-width,dp)/2.0_dp

  do i = 1, lx
    do j = 1, ly
      rnode = [real(i,dp),real(j,dp)] - rorigin
      if( norm2(rnode) >= radius ) then ! = radius is important because without it one has exists
        node(i,j,:)%nature = solid
      else
        node(i,j,:)%nature = fluid
      end if
    end do
  end do
END SUBROUTINE CONSTRUCT_CYLINDER

SUBROUTINE CONSTRUCT_COAXIAL_CYLINDER
  implicit none
  real(dp) :: radius1, radius2 ! radius of cylinder
  real(dp), dimension(2) :: rnode, rorigin ! coordinates of each node and center of cylinder in x,y coordinates
  integer(i2b) :: i, j, width1, width2 ! dummy
            integer(i2b) :: lx, ly, lz
            lx = supercell%geometry%dimensions%indiceMax(x)
            ly = supercell%geometry%dimensions%indiceMax(y)
            lz = supercell%geometry%dimensions%indiceMax(z)
  if( lx /= ly) stop 'wall=2 is for cylinders, which should have same lx and ly'
  if( lx<3 ) stop 'the diameter of the cylinder (lx) should be greater than 3'
  rorigin = [ real(lx+1,dp)/2.0_dp, real(ly+1,dp)/2.0_dp ]
  print*, 'rorigin in module is =', rorigin


  width1 = getinput%int("width1", defaultvalue=2) ! Number of nodes outside the cylinder
  width2 = getinput%int("width2", defaultvalue=2)

  radius1 = real(lx-width1,dp)/2.0_dp
  radius2 = real(lx-width2,dp)/2.0_dp

  print*, '*****************************************************************'
  print*, 'The inner cylinder has a radius of R1 = ', radius1
  print*, 'The outer cylinder has a radius of R2 = ', radius2
  print*, '*****************************************************************'

  do i = 1, lx
    do j = 1, ly
      rnode = [real(i,dp),real(j,dp)] - rorigin
      if( norm2(rnode) <= radius1 ) then ! = radius is important because without it one has exists
        node(i,j,:)%nature = solid
      elseif (norm2(rnode) > radius2) then
        node(i,j,:)%nature = solid
      else
        node(i,j,:)%nature = fluid
      end if
    end do
  end do
END SUBROUTINE CONSTRUCT_COAXIAL_CYLINDER



    subroutine construct_custom
    ! reads each point from input. Solid nodes only are indicated
        character(len=len("geom.in")), parameter :: filename="geom.in"
        integer(i1b), parameter :: u=77
        integer(i2b) :: i,j,k,lx,ly,lz,stat
        logical :: file_exists
        lx = supercell%geometry%dimensions%indiceMax(x)
        ly = supercell%geometry%dimensions%indiceMax(y)
        lz = supercell%geometry%dimensions%indiceMax(z)
        inquire(file="geom.in", exist=file_exists)
        if( .not. file_exists) then
            print*,"Cant find file containing custom geometry: "//filename//". Check lb.in if you really wanted custom geometry"
            stop "stop"
        end if
        node%nature = fluid ! everything but what is precised in geom.in is fluid
        open(unit=u,file=filename)
        stat= 0
        do while (.not. is_iostat_end(stat) )
            read(u,*,iostat=stat)i,j,k
            if(i<=0) then
                print*,"Index ",i," in 1st column (x column) of geom.in is negative or null. It must be between 1 and lx"
                stop
            end if
            if(j<=0) then
                print*,"Index ",j," in 2nd column (y column) of geom.in is negative or null. It must be between 1 and ly"
                stop
            end if
            if(k<=0) then
                print*,"Index ",k," in 3rd column (z column) of geom.in is negative or null. It must be between 1 and lz"
                stop
            end if
            if(i>lx) then
                print*,"Index ",i," in 1st column (x column) of geom.in is sup to nb of nodes. It must be between 1 and lx"
                stop
            end if
            if(j>ly) then
                print*,"Index ",j," in 2nd column (y column) of geom.in is sup to nb of nodes. It must be between 1 and ly"
                stop
            end if
            if(k>lz) then
                print*,"Index ",k," in 3rd column (z column) of geom.in is sup to nb of nodes. It must be between 1 and lz"
                stop
            end if
            node(i,j,k)%nature = solid
        end do
        close(u)
    end subroutine construct_custom
!
!
!
SUBROUTINE construct_pbm

    use system, only: node, supercell

    IMPLICIT NONE

    CHARACTER(2) :: magic_number
    CHARACTER(1) :: pix
    integer, parameter :: ncolumnmax=780
    CHARACTER(ncolumnmax) :: line
    INTEGER :: nline, ncolumn, i, j, nx, ny, nz
    nx = supercell%geometry%dimensions%indiceMax(x)
    ny = supercell%geometry%dimensions%indiceMax(y)
    nz = supercell%geometry%dimensions%indiceMax(z)
    if( nx/=1 ) error stop "lx must be 1 if geom.pbm is used"
    !
    ! pbm images are 2D.
    ! I chose arbitrarily that nx=1 (in fact, openmp parallelization is sometimes over z so bad idea to have nz=1)
    !
    OPEN(36, FILE="geom.pbm")
    READ(36,*) magic_number
    IF( magic_number /= "P1" ) THEN
        PRINT*, "geom.pbm doesnt seem to be valid. It's magic number (first line) is not P1 but",magic_number
        ERROR STOP
    END IF
    READ(36,*) ncolumn, nline

    if( ncolumn > ncolumnmax) then
        print*,"ncolumn from geom.pbm =",ncolumn
        print*,"maximum ncolumn implemented in construct_pbm =",ncolumnmax
        print*,"it is easy to implement higher numbers !"
        error stop "nevertheless it is not now. stop"
    end if

    if( ncolumn /= ny) then
        print*,"ncolumn from geom.pbm =",ncolumn
        print*,"ny in lb.in =",ny
        error stop "ncolumn /= ny in geom.pbm"
    end if

    if( nline /= nz) then
        print*,"nline from geom.pbm = ",nline
        print*,"nz in lb.in =",nz
        error stop "nline /= nz in geom.pbm"
    end if

    !
    ! init to fluid everywhere
    !
    node%nature = fluid

    !
    ! then find the solid nodes
    !
    do j=1,nline
        line=""
        read(36,*) line(1  :min(ncolumn, 70))
        if(ncolumn>70)  read(36,*) line(71 :min(ncolumn,140))
        if(ncolumn>140) read(36,*) line(141:min(ncolumn,210))
        if(ncolumn>210) read(36,*) line(211:min(ncolumn,280))
        if(ncolumn>280) read(36,*) line(281:min(ncolumn,350))
        if(ncolumn>350) read(36,*) line(351:min(ncolumn,420))
        if(ncolumn>420) read(36,*) line(421:min(ncolumn,490))
        if(ncolumn>490) read(36,*) line(491:min(ncolumn,560))
        if(ncolumn>560) read(36,*) line(561:min(ncolumn,630))
        if(ncolumn>630) read(36,*) line(631:min(ncolumn,700))
        if(ncolumn>700) read(36,*) line(701:min(ncolumn,770))
        if(ncolumn>770) error stop "look at module_geometry.f90 ncolumn>770 not implemented"
        do i=1,ncolumn
            select case (line(i:i))
            case("1")
                node(1,i,j)%nature = solid
            case("0")
            case default
                print*,"module_geometry. Pbm format allow 0 or 1 only"
                print*,"error with j,i=",j,i
                error stop
            end select
        end do
    end do
    close(36)
END SUBROUTINE construct_pbm


end module geometry

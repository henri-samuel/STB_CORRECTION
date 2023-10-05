! --------------------------------------!
! Copyright (c) 2023 Henri Samuel       !
! All rights reserved.                  !
! --------------------------------------!

! ================================================================================================= !

PURE SUBROUTINE lininterp3D(x0,x1,y0,y1,z0,z1,xi,yi,zi,S,si)
!$ USE OMP_LIB 
IMPLICIT NONE
REAL                        , INTENT(IN) :: x0,y0,z0,x1,y1,z1,xi,yi, zi
REAL, DIMENSION(0:1,0:1,0:1), INTENT(IN) :: S       ! Cell center values
REAL                        , INTENT(OUT):: si      ! Interpolated value at point xi,yi,zi
REAL                                     :: x, y, z
! Performs a trilinear interpolation 

x = (xi-x0) / (x1-x0)
y = (yi-y0) / (y1-y0)
z=  (zi-z0) / (z1-z0)

si = S(0,0,0) * (1.0-x) * (1.0-y) * (1.0-z) &
   + S(1,0,0) * x       * (1.0-y) * (1.0-z) &
   + S(0,1,0) * (1.0-x) * y       * (1.0-z) &
   + S(0,0,1) * (1.0-x) * (1.0-y) * z       &
   + S(1,0,1) * x       * (1.0-y) * z       &
   + S(0,1,1) * (1.0-x) * y       * z       &
   + S(1,1,0) * x       * y       * (1.0-z) &
   + S(1,1,1) * x       * y       * z

END SUBROUTINE lininterp3D

! ============================================================================ !

SUBROUTINE interp_vel3D(nx,ny,nz,x,y,z,dx,dy,dz,vel,vi)
IMPLICIT NONE
INTEGER                                  , INTENT(IN) :: nx,ny,nz
REAL                                     , INTENT(IN) :: x,y,z,dx,dy,dz
REAL, DIMENSION(1:3,0:nx+2,0:ny+2,0:nz+2), INTENT(IN) :: vel
REAL, DIMENSION(1:3)                     , INTENT(OUT):: vi
INTEGER                                               :: ix,iy,iz
REAL                                                  :: x0,y0,z0,x1,y1,z1

ix = INT(x/dx) +1   ;  x0 = FLOAT(ix-1)*dx ; x1 = x0+dx
iy = INT(y/dy) +1   ;  y0 = FLOAT(iy-1)*dy ; y1 = y0+dy
iz = INT(z/dz) +1   ;  z0 = FLOAT(iz-1)*dz ; z1 = z0+dz

CALL lininterp3D(x0,x1,y0,y1,z0,z1,x,y,z,vel(1,ix:ix+1,iy:iy+1,iz:iz+1),vi(1))
CALL lininterp3D(x0,x1,y0,y1,z0,z1,x,y,z,vel(2,ix:ix+1,iy:iy+1,iz:iz+1),vi(2))
CALL lininterp3D(x0,x1,y0,y1,z0,z1,x,y,z,vel(3,ix:ix+1,iy:iy+1,iz:iz+1),vi(3))

END SUBROUTINE interp_vel3D

! ======================================================================================= !

SUBROUTINE vxi3D(nx,ny,nz,x,y,z,dx,dy,dz,vel,vx_i)
IMPLICIT NONE
INTEGER                                  , INTENT(IN)  :: nx,ny,nz
REAL                                     , INTENT(IN)  :: x,y,z,dx,dy,dz
REAL, DIMENSION(1:3,0:nx+2,0:ny+2,0:nz+2), INTENT(IN)  :: vel
REAL                                     , INTENT(OUT) :: vx_i
INTEGER                                                :: ix,iy, iz
REAL                                                   :: x0,y0,z0,x1,y1,z1

ix = INT(x/dx) +1 ; x0 = FLOAT(ix-1)*dx ; x1 = x0+dx
iy = INT(y/dy) +1 ; y0 = FLOAT(iy-1)*dy ; y1 = y0+dy
iz = INT(z/dz) +1 ; z0 = FLOAT(iz-1)*dz ; z1 = z0+dz

CALL lininterp3D(x0,x1,y0,y1,z0,z1,x,y,z,vel(1,ix:ix+1,iy:iy+1,iz:iz+1),vx_i)

END SUBROUTINE vxi3D

! ======================================================================================= !

SUBROUTINE vyi3D(nx,ny,nz,x,y,z,dx,dy,dz,vel,vy_i)
IMPLICIT NONE
INTEGER                                  , INTENT(IN)  :: nx,ny,nz
REAL                                     , INTENT(IN)  :: x,y,z,dx,dy,dz
REAL, DIMENSION(1:3,0:nx+2,0:ny+2,0:nz+2), INTENT(IN)  :: vel
REAL                                     , INTENT(OUT) :: vy_i
INTEGER                                                :: ix,iy,iz
REAL                                                   :: x0,y0,z0,x1,y1,z1

ix = INT(x/dx) +1 ; x0 = FLOAT(ix-1)*dx ; x1 = x0+dx
iy = INT(y/dy) +1 ; y0 = FLOAT(iy-1)*dy ; y1 = y0+dy
iz = INT(z/dz) +1 ; z0 = FLOAT(iz-1)*dz ; z1 = z0+dz

CALL lininterp3D(x0,x1,y0,y1,z0,z1,x,y,z,vel(2,ix:ix+1,iy:iy+1,iz:iz+1),vy_i)

END SUBROUTINE vyi3D

! ======================================================================================= !

SUBROUTINE vzi3D(nx,ny,nz,x,y,z,dx,dy,dz,vel,vz_i)
IMPLICIT NONE
INTEGER                                  , INTENT(IN)  :: nx,ny,nz
REAL                                     , INTENT(IN)  :: x,y,z,dx,dy,dz
REAL, DIMENSION(1:3,0:nx+2,0:ny+2,0:nz+2), INTENT(IN)  :: vel
REAL                                     , INTENT(OUT) :: vz_i
INTEGER                                                :: ix,iy, iz
REAL                                                   :: x0,y0,z0,x1,y1,z1

ix = INT(x/dx) +1 ; x0 = FLOAT(ix-1)*dx ; x1 = x0+dx
iy = INT(y/dy) +1 ; y0 = FLOAT(iy-1)*dy ; y1 = y0+dy
iz = INT(z/dz) +1 ; z0 = FLOAT(iz-1)*dz ; z1 = z0+dz

CALL lininterp3D(x0,x1,y0,y1,z0,z1,x,y,z,vel(3,ix:ix+1,iy:iy+1,iz:iz+1),vz_i)

END SUBROUTINE vzi3D
! ============================================================================================ !

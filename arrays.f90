! --------------------------------------!
! Copyright (c) 2023 Henri Samuel       !
! All rights reserved.                  !
! --------------------------------------!

! ========================================================================================================================= !

SUBROUTINE assign_array_value(nx,ny,nz,val,array)
!$ USE OMP_LIB
IMPLICIT NONE
INTEGER                        , INTENT(IN)    :: nx,ny,nz
REAL                           , INTENT(IN)    :: val
REAL, DIMENSION(1:nx,1:ny,1:nz), INTENT(INOUT) :: array
INTEGER                                        :: ix,iy,iz
! This routine serves mostly as an OpenMP wrapper parallel assignment
! Assumes the array is 3D 

!$OMP PARALLEL IF(nx*ny*nz >1000) PRIVATE(ix,iy,iz)
!$OMP DO SCHEDULE(RUNTIME)
DO iz=1,nz
  DO iy=1,ny
    DO ix=1,nx
      array(ix,iy,iz) = val   
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE assign_array_value

! ================================================================================================================= !
SUBROUTINE assign_array2array(nx,ny,nz,array_in,array_out)
!$ USE OMP_LIB
IMPLICIT NONE
INTEGER                        , INTENT(IN)  :: nx,ny,nz
REAL, DIMENSION(1:nx,1:ny,1:nz), INTENT(IN)  :: array_in
REAL, DIMENSION(1:nx,1:ny,1:nz), INTENT(OUT) :: array_out
INTEGER                                      :: ix,iy,iz
! This routine serves mostly as an OpenMP wrapper parallel assignment
! Assumes the array is 3D 

!$OMP PARALLEL IF(nx*ny*nz >1000) PRIVATE(ix,iy,iz)
!$OMP DO SCHEDULE(RUNTIME)
DO iz=1,nz
  DO iy=1,ny
    DO ix=1,nx
      array_out(ix,iy,iz) = array_in(ix,iy,iz) 
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE assign_array2array
! ================================================================================================================= !

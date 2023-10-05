! --------------------------------------!
! Copyright (c) 2023 Henri Samuel       !
! All rights reserved.                  !
! --------------------------------------!

! ======================================================================================================================== !

SUBROUTINE apply_allVBCS(nx,ny,nz,Vx,Vy,Vz)
!$USE OMP_LIB
IMPLICIT NONE
INTEGER                              , INTENT(IN)    :: nx,ny,nz
REAL, DIMENSION(0:nx+1,0:ny+1,0:nz+1), INTENT(INOUT) :: Vx, Vy,Vz 
INTEGER                                              :: ix,iy,iz
LOGICAL                                              :: boundary

CALL apply_VxBCS(nx,ny,nz,Vx)
CALL apply_VyBCS(nx,ny,nz,Vy)
CALL apply_VzBCS(nx,ny,nz,Vz)


END SUBROUTINE apply_allVBCS

! ================================================================================================================ !
SUBROUTINE apply_VxBCS(nx,ny,nz,Vx)
!$USE OMP_LIB
IMPLICIT NONE
INTEGER                              , INTENT(IN)    :: nx,ny,nz
REAL, DIMENSION(0:nx+1,0:ny+1,0:nz+1), INTENT(INOUT) :: Vx
INTEGER                                              :: ix,iy,iz
LOGICAL                                              :: boundary

!$OMP PARALLEL DO  PRIVATE(ix,iy,iz,boundary)
DO iz=0, nz+1
  DO iy=0, ny+1
    DO ix=0, nx+1
      boundary = ix <= 1 .OR. ix >= nx+1 .OR. iy <= 0 .OR. iy >= ny+1 .OR. iz <= 0 .OR. iz >= nz+1
      IF (boundary) CALL xmom_bcs_pt(ix,iy,iz,nx,ny,nz,Vx)
    ENDDO
  ENDDO
ENDDO
!$OMP END PARALLEL DO

END SUBROUTINE apply_VxBCS

! ================================================================================================================ !
SUBROUTINE apply_VyBCS(nx,ny,nz,Vy)
IMPLICIT NONE
INTEGER                              , INTENT(IN)    :: nx,ny,nz
REAL, DIMENSION(0:nx+1,0:ny+1,0:nz+1), INTENT(INOUT) :: Vy
INTEGER                                              :: ix,iy,iz
LOGICAL                                              :: boundary

!$OMP PARALLEL DO  PRIVATE(ix,iy,iz,boundary)
DO iz=0, nz+1
  DO iy=0, ny+1
    DO ix=0, nx+1
      boundary = ix <= 0 .OR. ix >= nx+1 .OR. iy <= 1 .OR. iy >= ny+1 .OR. iz <= 0 .OR. iz >= nz+1
      IF (boundary) CALL ymom_bcs_pt(ix,iy,iz,nx,ny,nz,Vy)
    ENDDO
  ENDDO
ENDDO
!$OMP END PARALLEL DO

END SUBROUTINE apply_VyBCS

! ================================================================================================================ !
SUBROUTINE apply_VzBCS(nx,ny,nz,Vz)
!$USE OMP_LIB
IMPLICIT NONE
INTEGER                              , INTENT(IN)    :: nx,ny,nz
REAL, DIMENSION(0:nx+1,0:ny+1,0:nz+1), INTENT(INOUT) :: Vz
INTEGER                                              :: ix,iy,iz
LOGICAL                                              :: boundary

!$OMP PARALLEL DO  PRIVATE(ix,iy,iz,boundary)
DO iz=0, nz+1   
  DO iy=0, ny+1
    DO ix=0, nx+1
      boundary = ix <= 0 .OR. ix >= nx+1 .OR. iy <= 0 .OR. iy >= ny+1 .OR. iz <= 1 .OR. iz >= nz+1
      IF (boundary) CALL zmom_bcs_pt(ix,iy,iz,nx,ny,nz,Vz)
    ENDDO
  ENDDO
ENDDO
!$OMP END PARALLEL DO

END SUBROUTINE apply_VzBCS

! ================================================================================================================================== !

SUBROUTINE xmom_bcs_pt(ix,iy,iz,nx,ny,nz,Vx)
USE   bcs_mod, ONLY : VxBC
IMPLICIT NONE
INTEGER                                 , INTENT(IN)    :: ix,iy,iz,nx,ny,nz
REAL   , DIMENSION(0:nx+1,0:ny+1,0:nz+1), INTENT(INOUT) :: Vx

! -------------------------------------------z boundaries first --------------------------------- !
IF (iz==0) THEN
  Vx(ix,iy,iz) = 2.0*VxBC%Botval-Vx(ix,iy,iz+1) ; RETURN
ELSEIF (iz ==nz+1) THEN
  Vx(ix,iy,iz) = 2.0*VxBC%Topval-Vx(ix,iy,iz-1) ; RETURN
ENDIF
! ---------------------------------------- y boundaries next ------------------------------------ !

IF (iy ==0 ) THEN
    Vx(ix,iy,iz) = 2.0*VxBC%Southval-Vx(ix,iy+1,iz)  ; RETURN
ELSEIF (iy == ny+1 ) THEN
    Vx(ix,iy,iz) = 2.0*VxBC%Northval-Vx(ix,iy-1,iz)  ; RETURN
ENDIF
! ---------------------------------------- x boundaries next ------------------------------------ !
IF (ix ==1 ) THEN
    Vx(ix,iy,iz) = VxBC%Westval ; RETURN
ELSEIF (ix == nx+1) THEN
    vx(ix,iy,iz) = VxBC%Eastval  ; RETURN
ENDIF

END SUBROUTINE xmom_bcs_pt

SUBROUTINE ymom_bcs_pt(ix,iy,iz,nx,ny,nz,Vy)
USE   bcs_mod, ONLY : VyBC
IMPLICIT NONE
INTEGER                                 , INTENT(IN)    :: ix,iy,iz,nx,ny,nz
REAL   , DIMENSION(0:nx+1,0:ny+1,0:nz+1), INTENT(INOUT) :: Vy

! -------------------------------------------z boundaries first --------------------------------- !

IF (iz==0) THEN
    Vy(ix,iy,iz) = 2.0*VyBC%Botval-Vy(ix,iy,iz+1) ; RETURN
ELSEIF (iz ==nz+1) THEN
    Vy(ix,iy,iz) = 2.0*VyBC%Topval-Vy(ix,iy,iz-1) ; RETURN
ENDIF
! ---------------------------------------- y boundaries next ------------------------------------ !

IF (iy ==1 ) THEN
  Vy(ix,iy,iz) = VyBC%Southval  ; RETURN
ELSEIF (iy == ny+1 ) THEN
  Vy(ix,iy,iz) = VyBC%Northval  ; RETURN
ENDIF

! ---------------------------------------- x boundaries next ------------------------------------ !

IF (ix ==0 ) THEN
  Vy(ix,iy,iz) = 2.0*VyBC%Westval-Vy(ix+1,iy,iz)
ELSEIF (ix == nx+1 ) THEN
  Vy(ix,iy,iz) = 2.0*VyBC%Westval-Vy(ix-1,iy,iz)
ENDIF

END SUBROUTINE ymom_bcs_pt

! ================================================================================================================================== !
SUBROUTINE zmom_bcs_pt(ix,iy,iz,nx,ny,nz,Vz)
USE   bcs_mod, ONLY: VzBC
IMPLICIT NONE
INTEGER                                 , INTENT(IN)    :: ix,iy,iz,nx,ny,nz
REAL   , DIMENSION(0:nx+1,0:ny+1,0:nz+1), INTENT(INOUT) :: Vz
REAL                                                    :: x,y

! -------------------------------------------z boundaries first --------------------------------- !

IF (iz==1) THEN
  Vz(ix,iy,iz) = VzBC%Botval ; RETURN
ELSEIF (iz ==nz+1) THEN
  Vz(ix,iy,iz) = VzBC%Topval ; RETURN
ENDIF

! -------------------------------------------y boundaries next --------------------------------- !

IF (iy ==0  ) THEN
  Vz(ix,iy,iz) = 2.0*VzBC%Southval-Vz(ix,iy+1,iz)  ; RETURN
ELSEIF (iy ==ny+1  ) THEN
  Vz(ix,iy,iz) = 2.0*VzBC%Northval-Vz(ix,iy-1,iz)  ; RETURN
ENDIF

! ---------------------------------------- x boundaries next ------------------------------------ !

IF (ix ==0 ) THEN
  Vz(ix,iy,iz) = 2.0*VzBC%Westval-Vz(ix+1,iy,iz) ;  RETURN
ELSEIF (ix == nx+1 ) THEN
  Vz(ix,iy,iz) = 2.0*VzBC%Eastval-Vz(ix-1,iy,iz) ;  RETURN
ENDIF

END SUBROUTINE zmom_bcs_pt

! ================================================================================================================================== !

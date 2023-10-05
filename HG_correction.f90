! --------------------------------------!
! Copyright (c) 2023 Henri Samuel       !
! All rights reserved.                  !
! --------------------------------------!

! ================================================================================================================================== !
SUBROUTINE HG_correction(nx,ny,nz,dx,dy,dz,Vstag,P)
!$ USE OMP_LIB 
IMPLICIT NONE
INTEGER                                      , INTENT(IN)    :: nx,ny,nz
REAL                                         , INTENT(IN)    :: dx,dy,dz
REAL    , DIMENSION(1:3,0:nx+1,0:ny+1,0:nz+1), INTENT(INOUT) :: Vstag 
REAL    , DIMENSION(    0:nx+1,0:ny+1,0:nz+1), INTENT(INOUT) :: P
REAL    , DIMENSION(    0:nx+1,0:ny+1,0:nz+1)                :: div
INTEGER                                                      :: ix,iy,iz
REAL                                                         :: dxi,dyi,dzi
! Correct velocities using Helmholtz-Hodge decomp in order to obtain a divergence free velocity field
! The initial general field U can be decomposed as: U = V + grad(p) (1)
! With V solenoidal: div(V) = 0 (2)
! And p a scalar field such that grad(p) is irrotational
! Taking the divergence of (1) ==>  div(U) = div(V) + div(grad(p) 
!                                   div(U) =  0     + laplacian(p)  (3)
! (3) defines implicily a projection operator P such that V = P(U)       (NB: P /= p here)
! Therefore P has the following properties: P(V) = V  and P(grad(p) = 0  (NB: P /= p here)
! ==> solving for (3) using a muligrid method with newmann BCs we can compute V using (1):
! V = U - grad(p)
! ---------------------------------------------------------------------------------------------- !

dxi = 1.0/dx  ; dyi = 1.0/dy ; dzi = 1.0/dz
!------------------------------------------- Solve for pressure ------------------------------------------ !

div(0   ,:,:) = 0.0 ; div(:,   0,:) = 0.0 ; div(:,:,   0) = 0.0
div(nx+1,:,:) = 0.0 ; div(:,ny+1,:) = 0.0 ; div(:,:,nz+1) = 0.0

CALL calc_div(nx,ny,nz,dx,dy,dz,Vstag,div(1:nx,1:ny,1:nz))

CALL poisson_pressure_solver(nx,ny,nz,dx,dy,dz,Vstag,div,P)

! ------------------------------------------- Correct velocities ----------------------------------------- !

CALL correct_vel(nx,ny,nz,dx,dy,dz,P,Vstag,Vstag)

! -----------------------------------------------------------------------------------------------!

END SUBROUTINE HG_correction 
! ================================================================================================================================== !
SUBROUTINE correct_vel(nx,ny,nz,dx,dy,dz,P,Vstag,Vstag_corr)
!$ USE OMP_LIB 
IMPLICIT NONE
INTEGER                                  , INTENT(IN)  :: nx,ny,nz
REAL                                     , INTENT(IN)  :: dx,dy,dz
REAL, DIMENSION(    0:nx+1,0:ny+1,0:nz+1), INTENT(IN)  :: P
REAL, DIMENSION(1:3,0:nx+1,0:ny+1,0:nz+1), INTENT(IN)  :: Vstag
REAL, DIMENSION(1:3,0:nx+1,0:ny+1,0:nz+1), INTENT(OUT) :: Vstag_corr
INTEGER                                                :: ix,iy,iz
REAL                                                   :: dxi,dyi,dzi

dxi = 1.0/dx  ; dyi = 1.0/dy ; dzi = 1.0/dz

! ------------------------------------------- Correct velocities ----------------------------------------- !

!$OMP PARALLEL PRIVATE(ix,iy,iz) 
!$OMP DO SCHEDULE(RUNTIME)
DO iz=1,nz
  DO iy=1,ny
    DO ix=2,nx
      Vstag_corr(1,ix,iy,iz)=Vstag(1,ix,iy,iz)-(P(ix,iy,iz)-P(ix-1,iy  ,iz))*dxi
    ENDDO
  ENDDO
ENDDO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(RUNTIME)
DO iz=1,nz
  DO iy=2,ny
    DO ix=1,nx
      Vstag_corr(2,ix,iy,iz)=Vstag(2,ix,iy,iz)-(P(ix,iy,iz)-P(ix  ,iy-1,iz))*dyi
    ENDDO
  ENDDO
ENDDO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(RUNTIME)
DO iz=2,nz
  DO iy=1,ny
    DO ix=1,nx
      Vstag_corr(3,ix,iy,iz)=Vstag(3,ix,iy,iz)-(P(ix,iy,iz)-P(ix  ,iy,iz-1))*dzi
    ENDDO
  ENDDO
ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

CALL apply_allVBCS(nx,ny,nz,Vstag_corr(1,:,:,:),Vstag_corr(2,:,:,:),Vstag_corr(3,:,:,:))

! -----------------------------------------------------------------------------------------------!
END SUBROUTINE correct_vel

! ================================================================================================================================== !
SUBROUTINE poisson_pressure_solver(nx,ny,nz,dx,dy,dz,Vstag,rhs,sol)
IMPLICIT NONE
INTEGER                                       , INTENT(IN)    :: nx,ny,nz
REAL                                          , INTENT(IN)    :: dx,dy,dz
REAL     , DIMENSION(1:3,0:nx+1,0:ny+1,0:nz+1), INTENT(INOUT) :: Vstag
REAL     , DIMENSION(    0:nx+1,0:ny+1,0:nz+1), INTENT(INOUT) :: rhs
REAL     , DIMENSION(    0:nx+1,0:ny+1,0:nz+1), INTENT(OUT)   :: sol 

CALL multigrid(nx,ny,nz,dx,dy,dz,rhs,Vstag,sol)

END SUBROUTINE poisson_pressure_solver

! ======================================================================================================= !

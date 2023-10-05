! --------------------------------------!
! Copyright (c) 2023 Henri Samuel       !
! All rights reserved.                  !
! --------------------------------------!

! ================================================================================================================================== !
SUBROUTINE multigrid(nx,ny,nz,dx,dy,dz,rhs,Vstag,P)
USE    bcs_mod
USE solver_mod, ONLY: niter_finest,rincr,w_cont,ngmax,nvee,nmaxit,errmax_multi
USE  types_mod, ONLY: mgrid_all
IMPLICIT NONE
INTEGER                                              , INTENT(IN)    :: nx,ny,nz
REAL                                                 , INTENT(IN)    :: dx,dy,dz
REAL           , DIMENSION(    0:nx+1,0:ny+1,0:nz+1) , INTENT(IN)    :: rhs
REAL           , DIMENSION(1:3,0:nx+1,0:ny+1,0:nz+1) , INTENT(IN)    :: Vstag 
REAL           , DIMENSION(    0:nx+1,0:ny+1,0:nz+1) , INTENT(INOUT) :: P
REAL           , DIMENSION(1:3,0:nx+1,0:ny+1,0:nz+1)                 :: Vstag_new
REAL           , DIMENSION(1:3,1:nx+1,1:ny+1,1:nz+1)                 :: Vnode_new
REAL           , DIMENSION(    0:nx+1,0:ny+1,0:nz+1)                 :: div
INTEGER, SAVE                                                        :: ng
INTEGER                                                              :: ivee,ig,npost,npre,nxg,nyg,nzg
REAL                                                                 :: dxg,dyg,dzg
TYPE(mgrid_all), DIMENSION(:), ALLOCATABLE,SAVE                      :: MG
LOGICAL, SAVE                                                        :: first
REAL                                                                 :: rmsP,vrms,mean_div
REAL, DIMENSION(:), ALLOCATABLE,SAVE                                 :: res_cont

DATA first /.TRUE./

! Solves A u = f  (1)
! u: solution
! A: operator stencil (ie from FD or FV scheme)
! f: source/rhs 
! Let v be an approximation of u  (ie an unconverged solution)
! The error e is : e = u - v (2)
! Let r be the residue or defect such that : r = f - A v (3)
! Recast : (1) using (2) as A (v+e) = f 
! ==> A e = f - Av
! Then using (3): 
! ==>  A e = r  (4)  is the Residual equation solved recursively
! The MG method solves (1) on a hierarchy of grids 
! (a) by relaxing (1) on a fine grid 
! (b) Tranfering (restricting) the residue r to coarser grids
! (c) Solving for the error e using (4) with an initial guess of e=0
! (d) prolongating the corrections (i.e. the errors e) from coarse to fine grids
! (e) using the error to correct the approximate solution: v=v+e  (~uses (2))
! and keep on cycling until convergence...

! -------------------------------------------------------------------------------- !
IF (first) THEN 
  ng =MIN(INT(LOG(FLOAT(MIN(nx,ny,nz)))/LOG(2.0)),ngmax)
  ALLOCATE(MG(1:ng))  ; ALLOCATE(res_cont(1:ng))
  WRITE(*,*) '======================================================================'
  WRITE(*,*) 'Geometric multigrid ==> number of grid levels:', ng

  DO ig =1, ng  ! Finest (ig=1) to coarsest (ig=ng) grid level
    nxg = nx /(2**(ig-1))  ; MG(ig)%grid%nx = nxg
    nyg = ny /(2**(ig-1))  ; MG(ig)%grid%ny = nyg
    nzg = nz /(2**(ig-1))  ; MG(ig)%grid%nz = nzg 

    dxg = dx *FLOAT(2**(ig-1))  ; MG(ig)%grid%dx = dxg
    dyg = dy *FLOAT(2**(ig-1))  ; MG(ig)%grid%dy = dyg
    dzg = dz *FLOAT(2**(ig-1))  ; MG(ig)%grid%dz = dzg

    WRITE(*,*) 'Grid :', ig, 'nx=', nxg, 'ny=', nyg, 'nz=', nzg; CALL FLUSH(6)

    ALLOCATE(MG(ig)%val%cont(0:nxg+1,0:nyg+1,0:nzg+1))
    ALLOCATE(MG(ig)%rhs%cont(0:nxg+1,0:nyg+1,0:nzg+1))
    ALLOCATE(MG(ig)%res%cont(0:nxg+1,0:nyg+1,0:nzg+1)) 

    CALL assign_array_value(MG(ig)%grid%nx+2,MG(ig)%grid%ny+2,MG(ig)%grid%nz+2,0.0,MG(ig)%rhs%cont)
    CALL assign_array_value(MG(ig)%grid%nx+2,MG(ig)%grid%ny+2,MG(ig)%grid%nz+2,0.0,MG(ig)%res%cont)

  ENDDO
  first = .FALSE.

  WRITE(*,*) '======================================================================'
ENDIF

! ------------------------------------------------------------------------------------------------------------------- !

CALL assign_array2array(nx+2,ny+2,nz+2,  P(0:nx+1,0:ny+1,0:nz+1),MG(1)%val%cont(0:nx+1,0:ny+1,0:nz+1))
CALL assign_array2array(nx  ,ny  ,nz  ,rhs(1:nx  ,1:ny  ,1:nz  ),MG(1)%rhs%cont(1:nx  ,1:ny  ,1:nz  ))

! --------------------------------------------------------------------------------------------------------------- ! 

! Just standard relaxation
IF (ng ==1) THEN
  CALL relax(MG(ng)%grid,MG(ng)%rhs,MG(ng)%val,ng,nmaxit,w_cont,errmax_multi,res_cont(ng))
  CALL assign_array2array(nx+2,ny+2,nz+2,MG(ng)%val%cont(0:nx+1,0:ny+1,0:nz+1),P(0:nx+1,0:ny+1,0:nz+1))
  RETURN
ENDIF

! -------------------------------------------------------------------------------------------------------------- ! 

DO ivee =1,nvee
  
  DO ig=2,ng
    CALL assign_array_value(MG(ig)%grid%nx+2,MG(ig)%grid%ny+2,MG(ig)%grid%nz+2,0.0,MG(ig)%val%cont) ! important: initial guess for e (ig> 1) 
    CALL assign_array_value(MG(ig)%grid%nx+2,MG(ig)%grid%ny+2,MG(ig)%grid%nz+2,0.0,MG(ig)%rhs%cont) ! important: rhs of ig > 1 used to carry the residuals (see below)
  ENDDO

! ----------------------------------------- Downward V-cycle ------------------------------------------------------- !

  DO ig =1,ng-1
    npre = niter_finest*INT(rincr**FLOAT(ig-1))
    CALL    relax(MG(ig  )%grid,MG(ig)%rhs,MG(ig  )%val,ig,npre,w_cont,1e-15,res_cont(ig))
    CALL    resid(MG(ig  )%grid,MG(ig)%rhs,MG(ig  )%val,MG(ig)%res) ! Compute defect on fine grid
    CALL restrict(MG(ig+1)%grid,MG(ig)%res,MG(ig+1)%rhs)            ! Restrict defect to coarse grid 
  ENDDO

! ------------------------------------- Coarsest grid (bottom V) -------------------------------------------------- !

  CALL relax(MG(ng)%grid,MG(ng)%rhs,MG(ng)%val,ng,nmaxit,w_cont,1e-13,res_cont(ng))

! ------------------------------------------ Upward V-cycle ------------------------------------------------------- !

  DO ig =ng-1,1,-1
    npost = niter_finest*INT(rincr**FLOAT(ig-1))
    CALL    prolongate(      MG(ig)%grid,MG(ig+1)%val,MG(ig)%res)
    CALL       correct(ig==1,MG(ig)%grid,MG(ig  )%res,MG(ig)%val)
    CALL relax(              MG(ig)%grid,MG(ig  )%rhs,MG(ig)%val,ig,npost,w_cont,1e-15,res_cont(ig))
  ENDDO

! --------------------- Check divergence ---------------------------------- ! 
  CALL correct_vel(nx,ny,nz,dx,dy,dz,MG(1)%val%cont,Vstag,Vstag_new)
  CALL calc_div(nx,ny,nz,dx,dy,dz,Vstag_new,div(1:nx,1:ny,1:nz))
  CALL stag2node(nx,ny,nz,Vstag_new,Vnode_new)
  CALL calc_vrms(nx+1,ny+1,nz+1,Vnode_new,vrms) 
  CALL calc_mean(nx,ny,nz,ABS(div(1:nx,1:ny,1:nz)),mean_div)
! ------------------------------------------------------------------------- !

  IF (ivee ==1) WRITE(*,*)  'Poisson MG V-cycle       RMS residual           <div(v)>/Vrms           Vrms'
  WRITE(*,2) ivee,res_cont(1), mean_div/vrms, vrms
  IF (res_cont(1) < errmax_multi .AND. mean_div/vrms < 1E-6) EXIT
  IF (ivee ==nvee) STOP 'Convergence not reached, aborting'
ENDDO  ! ivee

WRITE(*,*)

CALL assign_array2array(nx+2,ny+2,nz+2,MG(1)%val%cont(0:nx+1,0:ny+1,0:nz+1),P(0:nx+1,0:ny+1,0:nz+1))

1 FORMAT('  Poisson MG V-cycle:',i5,',     RMS residual: ',1(E11.4,2x), '     <div(v)>/Vrms=',1(E11.4,2x), &
                  'Vrms=',E12.6,2x )
2 FORMAT(i5,"              |     ",2(E11.4,'          |  '),E12.6)

! -----------------------------------------------------------------------------------------------!

END SUBROUTINE multigrid

! ================================================================================================================================== !

SUBROUTINE relax(grid,rhs,val,ig,niter,w_cont,resmax,res_cont)
!$ USE OMP_LIB
USE solver_mod, ONLY: nmaxit,ngmax
USE  types_mod, ONLY: all_sol_array,grid_type
IMPLICIT NONE
TYPE(grid_type)      , INTENT(IN)                       :: grid
TYPE(all_sol_array)  , INTENT(IN)                       :: rhs
TYPE(all_sol_array)  , INTENT(INOUT)                    :: val
INTEGER              , INTENT(IN)                       :: ig,niter
REAL                 , INTENT(IN)                       :: w_cont,resmax
REAL                 , INTENT(OUT)                      :: res_cont
INTEGER                                                 :: nx,ny,nz,i,ix,iy,iz
REAL                                                    :: dx,dy,dz,dxi,dyi,dzi,dx2i,dy2i,dz2i
REAL                                                    :: rmsP,Pref 
REAL   , DIMENSION(1:grid%nx  ,1:grid%ny  ,1:grid%nz  ) :: res  ! NB: only for interior cells
REAL   , DIMENSION(0:grid%nx+1,0:grid%ny+1,0:grid%nz+1) :: P,cci,c_w,c_e,c_s,c_n,c_b,c_t,c_c
LOGICAL, DIMENSION(0:grid%nx+1,0:grid%ny+1,0:grid%nz+1) :: boundary

! ----------------------------------------------------------------------------------------------------- !

dx   = grid%dx ; dy   = grid%dy ; dz   = grid%dz
nx   = grid%nx ; ny   = grid%ny ; nz   = grid%nz
dxi  = 1.0/dx  ; dyi  = 1.0/dy  ; dzi  = 1.0/dz
dx2i = dxi**2  ; dy2i = dyi**2  ; dz2i = dzi**2

! --------------------------------------------------------------------------------------------- !
!$OMP PARALLEL IF (nx*ny*nz > 1000) PRIVATE(i,ix,iy,iz) 
!$OMP DO SCHEDULE (RUNTIME) 
  DO iz=0, nz+1
    DO iy=0, ny+1
      DO ix=0, nx+1
        boundary(ix,iy,iz) = ix==0 .OR. ix==nx+1 .OR. iy==0 .OR. iy==ny+1 .OR. iz==0 .OR. iz==nz+1

        P(ix,iy,iz) = val%cont(ix,iy,iz) 

        c_w(ix,iy,iz) = dx2i
        c_e(ix,iy,iz) = dx2i
        c_s(ix,iy,iz) = dy2i
        c_n(ix,iy,iz) = dy2i
        c_b(ix,iy,iz) = dz2i
        c_t(ix,iy,iz) = dz2i

        c_c(ix,iy,iz) = -(c_w(ix,iy,iz)+c_e(ix,iy,iz)+c_s(ix,iy,iz)+c_n(ix,iy,iz)+c_b(ix,iy,iz)+c_t(ix,iy,iz))

        cci(ix,iy,iz) = 1.0/c_c(ix,iy,iz)
        IF (.NOT. boundary(ix,iy,iz) ) res(ix,iy,iz)  = 0.0
      ENDDO
    ENDDO
  ENDDO
!$OMP END DO 

! --------------------------------------------------------------------------------------------- !

DO i=1, niter
  !$OMP DO SCHEDULE(RUNTIME) 
  DO iz=1, nz
    DO iy=1, ny
      DO ix=1, nx
        res(ix,iy,iz)= rhs%cont(ix,iy,iz)             -( c_C(ix,iy,iz)*P(ix  ,iy  ,iz  ) &
                     + c_W(ix,iy,iz)*P(ix-1,iy  ,iz  ) + c_E(ix,iy,iz)*P(ix+1,iy  ,iz  ) &
                     + c_S(ix,iy,iz)*P(ix  ,iy-1,iz  ) + c_N(ix,iy,iz)*P(ix  ,iy+1,iz  ) &
                     + c_B(ix,iy,iz)*P(ix  ,iy  ,iz-1) + c_T(ix,iy,iz)*P(ix  ,iy  ,iz+1) )
      ENDDO
    ENDDO
  ENDDO
  !$OMP END DO 

  !$OMP DO SCHEDULE(RUNTIME)
  DO iz=1, nz
    DO iy=1, ny
      DO ix=1 ,nx
        P(ix,iy,iz) = P(ix,iy,iz) + w_cont * res(ix,iy,iz) *cci(ix,iy,iz)
      ENDDO
    ENDDO
  ENDDO
  !$OMP END DO 

  !$OMP DO SCHEDULE(RUNTIME)
  DO iz=0, nz+1
    DO iy=0, ny+1
      DO ix=0, nx+1
        IF (boundary(ix,iy,iz)) CALL cont_bcs_pt(ig==1,ix,iy,iz,nx,ny,nz,P)
      ENDDO
    ENDDO
  ENDDO
  !$OMP END DO

! ---------------------------------- evaluate exit condition ------------------------ !
res_cont = 0.0 ; rmsP = 0.0
!$OMP BARRIER        ! Explicit synchronisation

!$OMP DO REDUCTION (+: rmsP)
DO iz=1,nz
  DO iy=1,ny
    DO ix=1,nx
      rmsP = rmsP + (res(ix,iy,iz))**2
    ENDDO
  ENDDO
ENDDO
!$OMP END DO 
res_cont = SQRT(rmsP/FLOAT(nx*ny*nz)) 

  IF (res_cont < resmax) EXIT
  IF (ngmax==1) WRITE(*,*)'Relax solver:',i,'resmax=',res_cont
ENDDO


!IF (ig==1) THEN 
!  Pref = P(1,1,1)  ! Reduces pressure values (Does NOT affect gradients, so OK)
!ELSE
!  Pref= 0.0
!ENDIF

Pref = 0.0
!$OMP WORKSHARE
val%cont(:,:,:) = P(:,:,:) - Pref 
!$OMP END WORKSHARE

!$OMP END PARALLEL

END SUBROUTINE relax
! ================================================================================================================================== !

SUBROUTINE resid(grid,rhs,val,res)
!$ USE OMP_lib
USE   types_mod, ONLY: all_sol_array, grid_type
IMPLICIT NONE
TYPE(grid_type)     , INTENT(IN)  :: grid
TYPE(all_sol_array) , INTENT(IN)  :: rhs,val
TYPE(all_sol_array) , INTENT(OUT) :: res
INTEGER                           :: ix,iy, iz
REAL                              :: dx,dy,dz,cc,cw,ce,cs,cn,cb,ct,dxi,dyi,dzi,dx2i,dy2i,dz2i
! -------------------------------------------------------------------------------------------- !
!$OMP PARALLEL WORKSHARE
  res%cont(:,:,:) = 0.0
!$OMP END PARALLEL WORKSHARE

dx   = grid%dx ; dy   = grid%dy ; dz   = grid%dz
dxi  = 1.0/dx  ; dyi  = 1.0/dy  ; dzi  = 1.0/dz
dx2i = dxi**2  ; dy2i = dyi**2  ; dz2i = dzi**2

!$OMP PARALLEL DO IF (grid%nx*grid%ny*grid%nz> 1000) SCHEDULE(RUNTIME)  PRIVATE(ix,iy,iz,cC,cW,cE,cS,cN,cB,cT) 
DO iz=1, grid%nz
  DO iy=1, grid%ny
    DO ix=1, grid%nx
      cw = dx2i
      ce = dx2i
      cs = dy2i
      cn = dy2i
      cb = dz2i
      ct = dz2i

      cc = -(cw+ce+cs+cn+cb+ct)
  
      res%cont(ix,iy,iz)= rhs%cont(ix,iy,iz) -( cC*val%cont(ix  ,iy  ,iz)                                 &
                                              + cW*val%cont(ix-1,iy  ,iz  ) + cE*val%cont(ix+1,iy  ,iz  ) &
                                              + cS*val%cont(ix  ,iy-1,iz  ) + cN*val%cont(ix  ,iy+1,iz  ) &
                                              + cB*val%cont(ix  ,iy  ,iz-1) + cT*val%cont(ix  ,iy  ,iz+1) )
    ENDDO
  ENDDO
ENDDO
!$OMP END PARALLEL DO 

! ------------------------------------------------------------------------------------------------------- !

END SUBROUTINE resid

! ================================================================================================================================== !

SUBROUTINE restrict(coarse_grid,fine_field,coarse_field)
USE types_mod, ONLY: all_sol_array,grid_type
IMPLICIT NONE
TYPE(grid_type)    , INTENT(IN)  :: coarse_grid
TYPE(all_sol_array), INTENT(IN)  :: fine_field
TYPE(all_sol_array), INTENT(OUT) :: coarse_field
INTEGER                          :: nxc,nyc,nzc
! ------------------------------------------------------------------------------------------------------------ !

nxc = coarse_grid%nx  ; nyc = coarse_grid%ny  ;  nzc = coarse_grid%nz

CALL restrict_P(nxc,nyc,nzc,fine_field%cont,coarse_field%cont)

! ----------------------------------------------------------------------------------------------------------- !

END SUBROUTINE restrict

! ================================================================================================================================= !
SUBROUTINE prolongate(fine_grid,coarse_field,fine_field)
USE types_mod, ONLY: all_sol_array,grid_type
IMPLICIT NONE
TYPE(grid_type)     , INTENT(IN)  :: fine_grid
TYPE(all_sol_array) , INTENT(IN)  :: coarse_field
TYPE(all_sol_array) , INTENT(OUT) :: fine_field
INTEGER                           :: nxf,nyf,nzf
! ------------------------------------------------------------------------------------------------------------ !
nxf = fine_grid%nx  ; nyf = fine_grid%ny  ;  nzf = fine_grid%nz

CALL prolongate_P(nxf,nyf,nzf,coarse_field%cont,fine_field%cont)

! ------------------------------------------------------------------------------------------------------------ !

END SUBROUTINE prolongate
! ================================================================================================================================= !
SUBROUTINE correct(finest,grid,res,val)
!$ USE OMP_LIB
USE types_mod, ONLY: all_sol_array,grid_type
IMPLICIT NONE
LOGICAL            , INTENT(IN)    :: finest
TYPE(grid_type)    , INTENT(IN)    :: grid
TYPE(all_sol_array), INTENT(IN)    :: res 
TYPE(all_sol_array), INTENT(INOUT) :: val 
INTEGER                            :: ix,iy,iz
INTEGER                            :: nx,ny,nz

! ------------------------------------------------------------------------------------------------------------ !
nx = grid%nx  ; ny = grid%ny  ;  nz = grid%nz
!$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(ix,iy,iz)
DO iz=1,nz
  DO iy=1,ny
    DO ix=1,nx
     val%cont(ix,iy,iz) = val%cont(ix,iy,iz) + res%cont(ix,iy,iz)
    ENDDO
  ENDDO
ENDDO
!$OMP END PARALLEL DO

IF (finest) THEN 
  !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(ix,iy,iz)
  DO iz=0,nz+1
    DO iy=0,ny+1
      DO ix=0,nx+1
        IF (ix==0 .OR. ix==nx+1 .OR. iy==0 .OR. iy==ny+1 .OR. iz==0 .OR. iz==nz+1) THEN 
          CALL cont_bcs_pt(finest,ix,iy,iz,nx,ny,nz,val%cont)
        ENDIF 
      ENDDO
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO
ENDIF
! ------------------------------------------------------------------------------------------------------------ !

END SUBROUTINE correct

! =================================================================================================================================== !
PURE SUBROUTINE cont_bcs_pt(finest,ix,iy,iz,nx,ny,nz,P)
IMPLICIT NONE
LOGICAL                                 , INTENT(IN)    :: finest
INTEGER                                 , INTENT(IN)    :: ix,iy,iz,nx,ny,nz
REAL   , DIMENSION(0:nx+1,0:ny+1,0:nz+1), INTENT(INOUT) :: P 

! ----------------------------------- z boundaries first ----------------------------------------- !

IF (iz ==0) THEN
  P(ix,iy,iz) =  P(ix,iy,iz+1) 
  RETURN
ELSEIF (iz ==nz+1) THEN
  P(ix,iy,iz) =  P(ix,iy,iz-1)
  RETURN
ENDIF

! ----------------------------------- y boundaries next ----------------------------------------- !
IF (iy ==0) THEN
  P(ix,iy,iz) =  P(ix,iy+1,iz)
  RETURN
ELSEIF (iy ==ny+1) THEN
  P(ix,iy,iz) =  P(ix,iy-1,iz)
  RETURN
ENDIF

! ---------------------------------- x boundaries next -------------------------------------- !

IF (ix ==0) THEN
  P(ix,iy,iz) = P(ix+1,iy,iz)
  RETURN
ELSEIF (ix ==nx+1) THEN
  P(ix,iy,iz) = P(ix-1,iy,iz)
  RETURN
ENDIF

END SUBROUTINE cont_bcs_pt

! ================================================================================================================================== !

SUBROUTINE prolongate_P(nxf,nyf,nzf,Pcoarse,Pfine)
!$USE OMP_LIB
IMPLICIT NONE
INTEGER                                          , INTENT(IN)  :: nxf,nyf,nzf
REAL   , DIMENSION(0:nxf/2+1,0:nyf/2+1,0:nzf/2+1), INTENT(IN)  :: Pcoarse
REAL   , DIMENSION(0:nxf+1  ,0:nyf+1  ,0:nzf+1  ), INTENT(OUT) :: Pfine
REAL                                                           :: facx,facy,facz,dPdx,dPdy,dPdz
INTEGER                                                        :: ixc,iyc,izc,ixf,iyf,izf,ii_x,ii_y,ii_z
INTEGER                                                        :: ifacxp,ifacxm,ifacyp,ifacym,ifaczp,ifaczm
LOGICAL, PARAMETER                                             :: injection = .FALSE. !.TRUE.
! Simple averaging as in [Albers, 2000]
! index of the lower left corner of the finer grid:
! ixf =2*ixc-1
! iyf =2*iyc-1
! izf =2*izc-1

!$OMP PARALLEL IF(nxf*nyf*nzf/8 > 1000 )                  &
!$OMP PRIVATE(ixc,iyc,izc,ixf,iyf,izf,ii_x,ii_y,ii_z    ) &
!$OMP PRIVATE(ifacxp,ifacxm,ifacyp,ifacym,ifaczp,ifaczm ) &
!$OMP PRIVATE(facx,facy,facz,dPdx,dPdy,dPdz             )       

dPdx = 0.0
dPdy = 0.0
dPdz = 0.0

!$OMP DO SCHEDULE(RUNTIME)
DO izf=0,nzf+1
  DO iyf=0,nyf+1
    DO ixf =0,nxf+1
      IF (ixf==0 .OR. ixf==nxf+1 .OR. iyf==0 .OR. iyf==nyf+1 .OR. izf==0 .OR. izf==nzf+1) Pfine(ixf,iyf,izf) = 0.0
    ENDDO
  ENDDO
ENDDO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(RUNTIME)
DO ixc=1,nxf/2
  DO iyc=1,nyf/2
    DO izc = 1,nzf/2
      DO ii_x=-1,0
        DO ii_y=-1,0
          DO ii_z=-1,0
            ixf = 2*ixc+ii_x
            iyf = 2*iyc+ii_y
            izf = 2*izc+ii_z
            IF (.NOT. injection)  THEN
              facx = 1.0+2*FLOAT(ixf-2*ixc) ; ifacxp = MAX(0,INT(facx))  ; ifacxm = MIN(0,INT(facx))
              facy = 1.0+2*FLOAT(iyf-2*iyc) ; ifacyp = MAX(0,INT(facy))  ; ifacym = MIN(0,INT(facy))
              facz = 1.0+2*FLOAT(izf-2*izc) ; ifaczp = MAX(0,INT(facz))  ; ifaczm = MIN(0,INT(facz))
              dPdx =facx*0.25*(Pcoarse(ixc+ifacxp,iyc       ,izc       )-Pcoarse(ixc+ifacxm,iyc       ,izc        ))
              dPdy =facy*0.25*(Pcoarse(ixc       ,iyc+ifacyp,izc       )-Pcoarse(ixc       ,iyc+ifacym,izc        ))
              dPdz =facz*0.25*(Pcoarse(ixc       ,iyc       ,izc+ifaczp)-Pcoarse(ixc       ,iyc        ,izc+ifaczm))
            ENDIF
            Pfine(ixf,iyf,izf) = Pcoarse(ixc,iyc,izc) + dPdx + dPdy + dPdz
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

END SUBROUTINE prolongate_P

! ================================================================================================================================= !
SUBROUTINE restrict_P(nxc,nyc,nzc,Pfine,Pcoarse)
!$ USE OMP_LIB
IMPLICIT NONE
INTEGER                                          , INTENT(IN)  :: nxc,nyc,nzc
REAL   , DIMENSION(0:2*nxc+1,0:2*nyc+1,0:2*nzc+1), INTENT(IN)  :: Pfine
REAL   , DIMENSION(0:  nxc+1,0:  nyc+1,0:  nzc+1), INTENT(OUT) :: Pcoarse
INTEGER                                                        :: ixc,iyc,izc,ixf,iyf,izf,iix,iiy,iiz

! Simple averaging as in [Albers, 2000]
! index of the lower left corner of the finer grid:
! ixf =2*ixc-1
! izf =2*izc-1

!$OMP PARALLEL IF(nxc*nxc*nxc> 1000) PRIVATE(ixc,iyc,izc,ixf,iyf,izf,iix,iiy,iiz)
!$OMP DO SCHEDULE(RUNTIME)
DO izc=1, nzc
izf =2*izc-1
  DO iyc=1, nyc
    iyf =2*iyc-1
    DO ixc=1,nxc
      ixf = 2*ixc-1
      Pcoarse(ixc,iyc,izc) = 0.0
      DO iiz =0,1
        DO iiy=0,1
          DO iix=0,1
            Pcoarse(ixc,iyc,izc) = Pcoarse(ixc,iyc,izc) + 0.125*Pfine(ixf+iix,iyf+iiy,izf+iiz)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE restrict_P

! ================================================================================================================================= !

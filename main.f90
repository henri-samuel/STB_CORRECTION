! -----------------------------------------!
! Author: Henri Samuel, 05/10/2023         !
! Version: 1.0.0                           !
! Distributed under the GNU GPLv3 Licence  !
! -----------------------------------------!
!
! STB_correction is a software that performs corrections on a 3D velocity field in a Cartesian domain 
! to apply rigid boundary conditions on all surfaces under a divergence-free constraint.
! STB_correction was developped to process velocity fields extracted from Particle Tracking Velocimetry in a convection tank using a Shake-The-Box approach
! implemented in the LaVision DaVis software.
! The resulting velocity field was found inconsistent with the existence of rigid boundaries near the edges of the domains where the PTV has trouble
! detecting and tracking particles. 
! The  STB algorithm used by DaVis did not have the capability of handling these issues, which motivated the implementation of STB_correction. 
!
! The STB_correction code proceeds as follows: 
! (1) Reads input (ASCII) file(s) to be corrected for velocities at the boundaries.
!  The corresponding format for each file is:
!  READ(11,*) iframe_recorded, time
!  READ(11,*) nx_orig,ny_orig,nz_orig
!  READ(11,*) Lx,Ly,Lz
!  READ(11,*) dx_orig,dy_orig,dz_orig
!  READ(11,*) x0_orig,y0_orig,z0_orig
!  READ(11,*) rhomin,rhomax,rho0
!  READ(11,*) etamin,etamax,eta0
!  READ(11,*) alphamin,alphamax,alpha0
!  READ(11,*) g,DeltaT,L,KappaT
 ! DO ix=1,nx_orig+1
!    DO iy=1,ny_orig+1
!      DO iz=1,nz_orig+1
!        READ(11,*) Vnode_orig(1:3,ix,iy,iz)
!      ENDDO
!    ENDDO
!  ENDDO
!
! (2) Linearly projects/interpolates velocity field components onto another regular Cartesian grid suitable for the multigrid solver 
!     (ideally using a number of grid cells that is a power of 2 in each direction).
! (3) Performs velocity field corrections using a Helmoltz-Hodge decomposition to ensure a solenoidal velocity field. 
!     This involves a geometric multrigrid solver on a Staggered FV-grid using Jacobi V-cycles.
! (4) Re-interpolates the corrected velocity field onto the original grid.
! (5) Outputs the corrected velocity field in ASCII (and optionally in vtk that can be visualized with ParaView for example) format.
!
! A single frame, or a series of frames can be processed.
! The code is written in Fortran and uses modules. 
! It could also calls for Linux commands to copy, remove, uncompress/compress input/output files (athough this is commented in the current version).
!
! Tests performded on more than 15 PTV Shake-The-Box convection experiments in a Cartesian tank have reliably and sucessfuly shown that the correction 
! essentially affects the velocity field near the boundaries. 
!
! This sofware is released to allow reproducing the results published in: 
! A. Walbecq, A. Limare, H. Samuel, Fully determined velocity field in a convection experiment with boundary conditions and zero divergence constraints (2023). 
! 
! STB_correction It is also intended to help academic researchers concerned with similar/related issues.
! The GNU GPLv3 licence requires in particular that any copy or altered version of the present code to be distributed (with its sources) 
! under the same licence.
! The author requires that the original source and the publication listed above to be cited upon each instance of STB_correction use.
!
! Results obtained with STB_correction are not subject to any guarantee from the author.
!
! --------------------------------------------------- Variables and Files ----------------------------------------------------- !
!
! Input parameters are specified in the "par" file:
!
! -------------------------------  IO parameters (from io_nl in par file)  ------------------------------------- !
!  dir_in          : Input directory 
!  stem_in         : Input stem file 
!  dir_out         : Output directory 
!  stem_out        : Output stem file  
!  iframe_start    : Starting frame for i/o files 
!  iframe_end      : Ending frame for i/o files
!  nwrite          : Frequency of frame output  (1 ==> all frames are output) 
!  write_vtk       : Switches on vtk output files (initial and final velocity, their difference, divergence and "Pressure") in addition to ASCII output 
!  binary_vtk      : Switches on binary mode for vtk output
!
! --------------------  (Multigrid) solver parameters (from solver_nl in par file)  ----------------------------- !
!  ngmax           : Maximum number of grid levels. 
!  nmaxit          : Maximum number of iteration on the coarser level.
!  w_cont          : Relaxation parameter (<1). Decrease if convergence is problematic. 
!  nvee            : Maximum number of V-cycles. Increase if convergence is problematic.
!  errmax_multi    : Maximum RMS residual allowed to declare convergence for the corrected velocity field (in addition to <div(V)> /V_RMS < 1E-6 required)  
!  niter_finest    : Maximum number of iteration on the finest grid. Increase if convergence is problematic.
!  nx              : Number of grid cells along the x-direction in the grid used for velocity correction. Best to use values close to nx_orig but that are power of 2.
!  ny              : Number of grid cells along the y-direction in the grid used for velocity correction. Best to use values close to ny_orig but that are power of 2. 
!  nz              : Number of grid cells along the z-direction in the grid used for velocity correction. Best to use values close to nz_orig but that are power of 2. 
!  non_dim         : Switches on non-dimensionalization of distances and velocities before applying the velocity correction 
!
! -----------------  Other parameters read in file header (see Unit 11 below for format)  ---------------------------------- !
! Each frame read corresponds to an ASCII file (Unit 11) named as:  TRIM(stem_in)//'_'//TRIM(str(iframe))//'.dat'
! The output ASCII file (corrected velocity field, Unit 12) has the same format and is named as TRIM(stem_out)//'_'//TRIM(str(iframe))//'.dat'
! The ASCII output/input files could compressed/uncompressed using the bzip2/bunzip2 command line.
!
! iframe_recorded : Frame number (should correspond to file name as well otherwise the code stops)
! time            : Time                                                            [s]
! nx_orig         : Number of grid cells along the x-direction in the original grid              ==> Influences the velocity correction
! ny_orig         : Number of grid cells along the y-direction in the original grid              ==> Influences the velocity correction
! nz_orig         : Number of grid cells along the z-direction in the original grid              ==> Influences the velocity correction
! Lx              : Tank/Domain length                                               [m]         ==> Influences the velocity correction
! Ly              : Tank/Domain depth                                                [m]         ==> Influences the velocity correction
! Lz              : Tank/Domain height                                               [m]         ==> Influences the velocity correction
! dx_orig         : Grid spacing along the x-direction in the original grid          [m]         ==> Influences the velocity correction
! dy_orig         : Grid spacing along the y-direction in the original grid          [m]         ==> Influences the velocity correction
! dz_orig         : Grid spacing along the z-direction in the original grid          [m]         ==> Influences the velocity correction
! x0_orig         : x-coordinate of the front,lower,left corner in the original grid [m]         ==> Influences the velocity correction
! y0_orig         : y-coordinate of the front,lower,left corner in the original grid [m]         ==> Influences the velocity correction
! z0_orig         : z-coordinate of the front,lower,left corner in the original grid [m]         ==> Influences the velocity correction
! rhomin          : Minimum value for the fluid density in the tank                  [kg/m^3]   
! rhomax          : Maximum value for the fluid density in the tank                  [kg/m^3]
! rho0            : Reference value for the fluid density in the tank                [kg/m^3]
! etamin          : Minimum value for the fluid dynamic viscosity in the tank        [Pa s]
! etamax          : Maximum value for the fluid dynamic viscosity in the tank        [Pa s]
! eta0            : Reference value for the fluid dynamic viscosity in the tank      [Pa s]
! alphamin        : Minimum value for the fluid thermal expansion in the tank        [K^(-1)]
! alphamax        : Maximum value for the fluid dynamic viscosity in the tank        [Pa s]
! alpha0          : Reference value for the fluid dynamic viscosity in the tank      [Pa s]
! g               : Gravity acceleration                                             [m/s^2]
! DeltaT          : Temperature difference between the horizontal surfaces           [K]
! L               : Length scale chosen here to be height of the tank                [m]  Typically L=Lz here since gravity is parallel to the z-axis 
! KappaT          : Thermal diffusivity for the fluid in the tank                    [m^2/s]
! -------------------------------------------------- Compiling ------------------------------------------------- !
!
! Compiling requires make, a 90+ Fortran compiler (gfortran,ifort...) and OpenMP to speed up calculations on shared memory machines.
! Compilation rules are contained in the "makefile" file.
! (1) Edit the makefile configuration file if necessary (two examples are in the MAKE_INC local directory) 
! (2) Simply type "make"
! Upon sucessful compilation, the executable "stb_correction.exe" will be generated.
! To clean up eveything type "make cleanall".
!
! --------------------------------------------------- Running -------------------------------------------------- !
!
! (1) [Optional] Set the number of OpenMP threads (e.g., "export OMP_NUM_THREADS=4" for 4 threads using bash shell)
! (2) Type "./stb_correction.exe"
!
! ---------------------------------------------------- Example ------------------------------------------------- !
!
! An example is set up in the "par" file contained in the current directory. 
! It uses 10 compressed input files contained in the local directory INPUT_DIR.
!
! -------------------------------------------------------------------------------------------------------------- !
PROGRAM SBT_corr
!$ USE OMP_LIB 
USE    bcs_mod
USE solver_mod
IMPLICIT NONE
INTEGER                               :: nx,ny,nz,nx_orig,ny_orig,nz_orig,ix,iy,iz
INTEGER                               :: i,ixi,iyi,izi,nwrite
INTEGER                               :: iframe_recorded,iframe, iframe_start, iframe_end
REAL                                  :: arx,ary,arz,dx,dy,dz,dx_orig,dy_orig,dz_orig,Ra,Pr
REAL                                  :: x0,x1,y0,y1,z0,z1,xi,yi,zi
REAL                                  :: arx_orig, ary_orig, arz_orig,x0_orig,y0_orig,z0_orig,L,KappaT,DeltaT
REAL                                  :: Lx,Ly,Lz,alphamin,alphamax,alpha0
REAL                                  :: rhomin,rhomax,etamin,etamax,vrms,rho0,eta0,g,div_mean
REAL, DIMENSION(  :,:,:), ALLOCATABLE :: P,div,div_orig
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: Vstag_orig, Vnode_orig, Vstag, Vnode,deltaVnode
LOGICAL                               :: non_dim, binary_vtk,write_vtk
CHARACTER(LEN=50)                     :: dir_in,dir_out,namefile_in,namefile_out,stem_in,stem_out
CHARACTER(LEN=20)  , EXTERNAL         :: str

NAMELIST /    io_nl/dir_in,stem_in,dir_out,stem_out,iframe_start,iframe_end,nwrite,write_vtk,binary_vtk
NAMELIST /solver_nl/ngmax,nmaxit,w_cont,nvee,errmax_multi,niter_finest,nx,ny,nz,non_dim
                                    
! ------------------------------- Initialization ------------------------------------------ !
DO i=1,50
  dir_in(i:i)      =' '
  dir_out(i:i)     =' '
  namefile_in(i:i) =' '
  namefile_out(i:i)=' '
  stem_in(i:i)     =' '
  stem_out(i:i)    =' '
ENDDO

VxBC%Botval   = 0.0  ; VyBC%Botval   = 0.0  ; VzBC%Botval   = 0.0
VxBC%Topval   = 0.0  ; VyBC%Topval   = 0.0  ; VyBC%Topval   = 0.0
VxBC%Westval  = 0.0  ; VyBC%Westval  = 0.0  ; VzBC%Westval  = 0.0
VxBC%Eastval  = 0.0  ; VyBC%Eastval  = 0.0  ; VzBC%Eastval  = 0.0
VxBC%Southval = 0.0  ; VyBC%Southval = 0.0  ; VzBC%Southval = 0.0
VxBC%Northval = 0.0  ; VyBC%Northval = 0.0  ; VzBC%Northval = 0.0

! ---------------------- Default values for input parameters ---------------------------- !
ngmax             =  6
nmaxit            = 2000
w_cont            = 0.95 
nvee              = 50
errmax_multi      = 1e-5
rincr             = 2.0
niter_finest      = 8 
nx                = 64 
ny                = 64
nz                = 64
non_dim           = .FALSE.
! -------------------------- Default values for input parameters   -------------- !
dir_in            = '/gpfs/scratch/walbecq/calibration_corrigee/'
stem_in           = 'STB_in'
dir_out           = 'OUTPUT/'
stem_out          = 'STB_Vbox_out'
iframe_start      = 1
iframe_end        = 2000 
nwrite            = 1
write_vtk         = .TRUE.
binary_vtk        = .TRUE.
! ------------------------------------------------------------------------------------------- ! 

OPEN(1,FILE='par',STATUS='OLD')
READ(1,io_nl)
READ(1,solver_nl)
CLOSE(1)

! -------------------------------------------------- Main time/frame-loop -------------------------- !
DO iframe=iframe_start, iframe_end

! ------------------------------------- SBT uncorrected file  -----------------------
  namefile_in =TRIM(dir_in)//TRIM(stem_in)//'_'//TRIM(str(iframe))//'.dat'
!  CALL SYSTEM('bunzip2 '//TRIM(namefile_in)//'.bz2 ')
  WRITE(*,*) '------------------------------------------------------------------- '
  WRITE(*,*) 'Reading file: '//TRIM(namefile_in)
  OPEN(11,FILE=TRIM(namefile_in))
  READ(11,*) iframe_recorded, time
  READ(11,*) nx_orig,ny_orig,nz_orig
  READ(11,*) Lx,Ly,Lz
  READ(11,*) dx_orig,dy_orig,dz_orig
  READ(11,*) x0_orig,y0_orig,z0_orig
  READ(11,*) rhomin,rhomax,rho0
  READ(11,*) etamin,etamax,eta0
  READ(11,*) alphamin,alphamax,alpha0
  READ(11,*) g,DeltaT,L,KappaT  

  IF (iframe_recorded /= iframe) THEN 
    WRITE(*,*) iframe_recorded , iframe
    STOP 'iframe_recorded /= iframe '
  ENDIF

  IF (non_dim) THEN
    dx_orig= dx_orig/L
    dy_orig= dy_orig/L
    dz_orig= dz_orig/L
    x0_orig= x0_orig/L
    y0_orig= y0_orig/L
    z0_orig= z0_orig/L
  ENDIF

  arx_orig = nx_orig * dx_orig
  ary_orig = ny_orig * dy_orig
  arz_orig = nz_orig * dz_orig

  WRITE(*,*) 'Processing frame:', iframe
  WRITE(*,*) 'arx_orig=',arx_orig, 'ary_orig=', ary_orig, 'arz_orig=',arz_orig
  WRITE(*,*) 'x orig range:',x0_orig,x0_orig+arx_orig; IF (non_dim) WRITE(*,*) 'dimensional:',x0_orig*L,(x0_orig+arx_orig)*L
  WRITE(*,*) 'y orig range:',y0_orig,y0_orig+ary_orig; IF (non_dim) WRITE(*,*) 'dimensional:',y0_orig*L,(y0_orig+ary_orig)*L
  WRITE(*,*) 'z orig range:',z0_orig,z0_orig+arz_orig; IF (non_dim) WRITE(*,*) 'dimensional:',z0_orig*L,(z0_orig+arz_orig)*L

  IF (iframe == iframe_start) THEN
    ALLOCATE(div_orig(      1:nx_orig  ,1:ny_orig  ,1:nz_orig  ))
    ALLOCATE(Vstag_orig(1:3,0:nx_orig+1,0:ny_orig+1,0:nz_orig+1))
    ALLOCATE(Vnode_orig(1:3,1:nx_orig+1,1:ny_orig+1,1:nz_orig+1))
    !$OMP PARALLEL WORKSHARE
    Vstag_orig(:,:,:,:) = 0.0
    Vnode_orig(:,:,:,:) = 0.0
    !$OMP END PARALLEL WORKSHARE
  ENDIF

  DO ix=1,nx_orig+1
    DO iy=1,ny_orig+1
      DO iz=1,nz_orig+1
      READ(11,*) Vnode_orig(:,ix,iy,iz)
      ENDDO
    ENDDO
  ENDDO

  CLOSE(11)

  ! ----------------------------------------------------------------------------------------------------- !

  CALL calc_vrms(nx_orig+1,ny_orig+1,nz_orig+1,Vnode_orig,vrms)
  WRITE(*,*) 'vrms dim (m/s)=', SNGL(vrms), 'vms ND=',vrms*L/KappaT,'Re= rho0  vrms L /eta0=  ', SNGL(vrms*rho0*L/eta0)
  WRITE(*,*) 'Pr = eta/(rho*Kappa) =',SNGL(eta0  /(rho0*KappaT)), 'Kappa re-estimated=',SNGL(eta0/(rho0*g))

  CALL Vnode2Vstag(nx_orig,ny_orig,nz_orig,Vnode_orig,Vstag_orig)

  CALL calc_div(nx_orig,ny_orig,nz_orig,dx_orig,dy_orig,dz_orig,Vstag_orig,div_orig)

  IF (MOD(iframe,nwrite) ==0 .AND. write_vtk) THEN
    CALL vtk_scalar_asciibin(nx_orig,ny_orig,nz_orig,dx_orig,dy_orig,dz_orig,iframe,'div_initial',binary_vtk,div_orig)
    CALL vtk_vector_asciibin(nx_orig+1,ny_orig+1,nz_orig+1,dx_orig,dy_orig,dz_orig,iframe,'V_initial_dim',binary_vtk,Vnode_orig)
  ENDIF


  IF (non_dim) THEN
    !$OMP PARALLEL WORKSHARE 
    Vnode_orig(:,:,:,:) = Vnode_orig(:,:,:,:) / (KappaT/L)
    !$OMP END PARALLEL WORKSHARE 
    IF (MOD(iframe,nwrite) ==0 .AND. write_vtk) THEN 
    CALL vtk_vector_asciibin(nx_orig+1,ny_orig+1,nz_orig+1,dx_orig,dy_orig,dz_orig,iframe,'V_initial',binary_vtk,Vnode_orig)
    ENDIF
    CALL calc_vrms(nx_orig+1,ny_orig+1,nz_orig+1,Vnode_orig,vrms)
    WRITE(*,*) 'Vrms ND=', vrms
  ENDIF

  CALL FLUSH(6)
! ---------------------------------------------------------------------------------------------------- !

! Interpolate orig grid onto another grid
  IF (iframe== iframe_start) THEN
    arx = Lx; IF (non_dim) arx=arx/L; dx =arx/FLOAT(nx)
    ary = Ly; IF (non_dim) ary=ary/L; dy =ary/FLOAT(ny)
    arz = Lz; IF (non_dim) arz=arz/L; dz =arz/FLOAT(nz)

    ALLOCATE(div(           1:nx  ,1:ny  ,1:nz  ))
    ALLOCATE(P(             0:nx+1,0:ny+1,0:nz+1))  
    ALLOCATE(Vstag(     1:3,0:nx+1,0:ny+1,0:nz+1)) 
    ALLOCATE(Vnode(     1:3,1:nx+1,1:ny+1,1:nz+1)) 
    ALLOCATE(deltaVnode(1:3,1:nx+1,1:ny+1,1:nz+1))  
    !$OMP PARALLEL WORKSHARE 
    div(    :,:,:) = 0.0
    Vstag(:,:,:,:) = 0.0
    Vnode(:,:,:,:) = 0.0
    !$OMP END PARALLEL WORKSHARE  
  ENDIF

  !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(xi,yi,zi,ix,iy,iz,ixi,iyi,izi,x0,y0,z0,x1,y1,z1)
  DO iz=2,nz
    zi = FLOAT(iz-1)*dz ; izi = INT((zi-z0_orig)/dz_orig)+1 ; z0=z0_orig + FLOAT(izi-1)*dz_orig ; z1=z0+dz_orig
    DO iy=2,ny
      yi = FLOAT(iy-1)*dy  ; iyi = INT((yi-y0_orig)/dy_orig)+1 ; y0=y0_orig + FLOAT(iyi-1)*dy_orig ; y1=y0+dy_orig
      DO ix=2,nx
        xi = FLOAT(ix-1)*dx  ; ixi = INT((xi-x0_orig)/dx_orig)+1 ; x0=x0_orig + FLOAT(ixi-1)*dx_orig ; x1=x0+dx_orig
        CALL lininterp3D(x0,x1,y0,y1,z0,z1,xi,yi,zi,Vnode_orig(1,ixi:ixi+1,iyi:iyi+1,izi:izi+1),Vnode(1,ix,iy,iz)) 
        CALL lininterp3D(x0,x1,y0,y1,z0,z1,xi,yi,zi,Vnode_orig(2,ixi:ixi+1,iyi:iyi+1,izi:izi+1),Vnode(2,ix,iy,iz)) 
        CALL lininterp3D(x0,x1,y0,y1,z0,z1,xi,yi,zi,Vnode_orig(3,ixi:ixi+1,iyi:iyi+1,izi:izi+1),Vnode(3,ix,iy,iz))  
      ENDDO
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO 

! ------------------------------------------------------------------------------------------- !

  CALL Vnode2Vstag(nx,ny,nz,Vnode,Vstag) ! Includes BCs
  CALL calc_div(nx,ny,nz,dx,dy,dz,Vstag,div)

  IF (MOD(iframe,nwrite) ==0 .AND. write_vtk) THEN
    CALL vtk_scalar_asciibin(nx  ,ny  ,nz  ,dx,dy,dz,iframe,'div_initial_interp',binary_vtk,div)
    CALL vtk_vector_asciibin(nx+1,ny+1,nz+1,dx,dy,dz,iframe,'V_initial_interp'  ,binary_vtk,Vnode)
  ENDIF

! ---------------------------------------------------------------------------------------------- !

  IF (MOD(iframe,nwrite) ==0 .AND. write_vtk) THEN
    !$OMP PARALLEL WORKSHARE 
    deltaVnode(:,:,:,:) =Vnode(:,:,:,:)
    !$OMP END PARALLEL WORKSHARE 
  ENDIF

  !$OMP PARALLEL WORKSHARE 
    P(:,:,:) = 0.0  ! better to reinitialize for each frame (maybe because P here is not really a pressure)
  !$OMP END PARALLEL WORKSHARE 

  CALL HG_correction(nx,ny,nz,dx,dy,dz,Vstag,p)

  CALL calc_div(nx,ny,nz,dx,dy,dz,Vstag,div)
  CALL calc_mean(nx,ny,nz,ABS(div),div_mean)
  WRITE(*,*) 'div(v): max=', maxval(abs(div)),'mean=',div_mean,'mean/vrms=',div_mean/vrms

  CALL stag2node(nx,ny,nz,Vstag,Vnode)

  WRITE(*,*) 'vrms initial=', SNGL(vrms)
  CALL calc_vrms(nx+1,ny+1,nz+1,Vnode,vrms)
  WRITE(*,*) 'vrms final=', SNGL(vrms)

  IF (MOD(iframe,nwrite) ==0 .AND. write_vtk) THEN
    !$OMP PARALLEL WORKSHARE 
    deltaVnode(:,:,:,:) =deltaVnode(:,:,:,:) - Vnode(:,:,:,:)
    !$OMP END PARALLEL WORKSHARE 
    CALL vtk_scalar_asciibin(nx  ,ny  ,nz  ,dx,dy,dz,iframe,'div_final',binary_vtk,div       )    
    CALL vtk_scalar_asciibin(nx+2,ny+2,nz+2,dx,dy,dz,iframe,'P_final'  ,binary_vtk,p         )
    CALL vtk_vector_asciibin(nx+1,ny+1,nz+1,dx,dy,dz,iframe,'V_final'  ,binary_vtk,Vnode     )
    CALL vtk_vector_asciibin(nx+1,ny+1,nz+1,dx,dy,dz,iframe,'deltaV'   ,binary_vtk,deltaVnode)
  ENDIF

! ------------------------------------------------------------------------------------------------- !
  namefile_out = TRIM(dir_out)//TRIM(stem_out)//'_'//TRIM(str(iframe))//'.dat' 
  WRITE(*,*) 'Output results in :',TRIM(namefile_out)
  OPEN(12,FILE=TRIM(namefile_out))
  WRITE(12,*) iframe_recorded, time
  WRITE(12,*) nx,ny,nz
  WRITE(12,*) Lx,Ly,Lz
  WRITE(12,*) dx,dy,dz
  WRITE(12,*) rhomin,rhomax,rho0
  WRITE(12,*) etamin,etamax,eta0 
  WRITE(12,*) alphamin,alphamax,alpha0
  WRITE(12,*) g,DeltaT,L,KappaT

   DO ix=1,nx+1
     DO iy=1,ny+1
       DO iz=1,nz+1
         WRITE(12,12)  Vnode(:,ix,iy,iz)
       ENDDO
     ENDDO
   ENDDO

12 FORMAT(3(E10.4,1x))

  CLOSE(12)
!  CALL SYSTEM('bzip2 '//TRIM(namefile_out))
ENDDO

END PROGRAM SBT_corr
! ================================================================================================== !
CHARACTER(LEN=20) FUNCTION str(k) !   Convert an integer to string 
IMPLICIT NONE
INTEGER, INTENT(IN) :: k

WRITE (str, *) k ; str = ADJUSTL(str)

END FUNCTION str
! ================================================================================================== !

SUBROUTINE Vnode2Vstag(nx,ny,nz,Vnode,Vstag)
!$ USE  OMP_LIB 
IMPLICIT NONE
INTEGER                                  , INTENT(IN) :: nx,ny,nz
REAL, DIMENSION(1:3,1:nx+1,1:ny+1,1:nz+1), INTENT(IN) :: Vnode
REAL, DIMENSION(1:3,0:nx+1,0:ny+1,0:nz+1), INTENT(OUT):: Vstag
INTEGER                                               :: ix,iy,iz

!$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(ix,iy,iz)
DO ix=1,nx+1
  DO iy=1,ny
    DO iz=1,nz
      Vstag(1,ix,iy,iz) = 0.25*(Vnode(1,ix,iy,iz)+Vnode(1,ix,iy+1,iz)+Vnode(1,ix,iy,iz+1)+Vnode(1,ix,iy+1,iz+1))
    ENDDO
  ENDDO
ENDDO
!$OMP END PARALLEL DO 


!$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(ix,iy,iz) 
DO ix=1,nx
  DO iy=1,ny+1
    DO iz=1,nz
      Vstag(2,ix,iy,iz) = 0.25*(Vnode(2,ix,iy,iz)+Vnode(2,ix+1,iy,iz)+Vnode(2,ix,iy,iz+1)+Vnode(2,ix+1,iy,iz+1))
    ENDDO
  ENDDO
ENDDO
!$OMP END PARALLEL DO 


!$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(ix,iy,iz)
DO ix=1,nx
  DO iy=1,ny
    DO iz=1,nz+1
      Vstag(3,ix,iy,iz) = 0.25*(Vnode(3,ix,iy,iz)+Vnode(3,ix+1,iy,iz)+Vnode(3,ix,iy+1,iz)+Vnode(3,ix+1,iy+1,iz))
    ENDDO
  ENDDO
ENDDO
!$OMP END PARALLEL DO 

CALL apply_allVBCS(nx,ny,nz,Vstag(1,:,:,:),Vstag(2,:,:,:),Vstag(3,:,:,:))

END SUBROUTINE Vnode2Vstag

! ================================================================================================== !

SUBROUTINE stag2node(nx,ny,nz,Ustag,Unode)
!$ USE  OMP_LIB 
IMPLICIT NONE
INTEGER                                  , INTENT(IN)  :: nx,ny,nz
REAL, DIMENSION(1:3,0:nx+1,0:ny+1,0:nz+1), INTENT(IN)  :: Ustag
REAL, DIMENSION(1:3,1:nx+1,1:ny+1,1:nz+1), INTENT(OUT) :: Unode
INTEGER                                                :: ix,iy,iz

!$OMP PARALLEL DO IF (nx*ny*nz > 1000) SCHEDULE(STATIC) PRIVATE(ix,iy,iz)
DO iz=1,nz+1
  DO iy=1,ny+1
    DO ix=1,nx+1
    Unode(1,ix,iy,iz) =0.25*( Ustag(1,ix,iy-1,iz-1) + Ustag(1,ix,iy-1,iz  )       &
                            + Ustag(1,ix,iy  ,iz  ) + Ustag(1,ix,iy  ,iz-1) )

    Unode(2,ix,iy,iz) =0.25*( Ustag(2,ix-1,iy  ,iz-1) + Ustag(2,ix-1,iy  ,iz  )   &
                            + Ustag(2,ix  ,iy  ,iz  ) + Ustag(2,ix  ,iy  ,iz-1) )

    Unode(3,ix,iy,iz) =0.25*( Ustag(3,ix-1,iy-1,iz) + Ustag(3,ix-1,iy  ,iz)       &
                            + Ustag(3,ix  ,iy  ,iz) + Ustag(3,ix  ,iy-1,iz) )
    ENDDO
  ENDDO
ENDDO
!$OMP END PARALLEL DO

END SUBROUTINE stag2node
! ========================================================================================================================= !
SUBROUTINE calc_div(nx,ny,nz,dx,dy,dz,Vstag,div)
!$ USE  OMP_LIB 
IMPLICIT NONE
INTEGER                                  , INTENT(IN)  :: nx,ny,nz
REAL                                     , INTENT(IN)  :: dx,dy,dz
REAL, DIMENSION(1:3,0:nx+1,0:ny+1,0:nz+1), INTENT(IN)  :: Vstag
REAL, DIMENSION(    1:nx  ,1:ny  ,1:nz  ), INTENT(OUT) :: div 
INTEGER                                                :: ix,iy,iz

!$OMP PARALLEL DO IF (nx*ny*nz > 1000) SCHEDULE(STATIC) PRIVATE(ix,iy,iz)
DO iz=1,nz
  DO iy=1,ny
    DO ix=1,nx
      div(ix,iy,iz) = (Vstag(1,ix+1,iy  ,iz  ) - Vstag(1,ix,iy,iz))/dx &
                    + (Vstag(2,ix  ,iy+1,iz  ) - Vstag(2,ix,iy,iz))/dy &
                    + (Vstag(3,ix  ,iy  ,iz+1) - Vstag(3,ix,iy,iz))/dz
    ENDDO
  ENDDO
ENDDO
!$OMP END PARALLEL DO

END SUBROUTINE calc_div

! ========================================================================================================================= !
SUBROUTINE calc_vrms(nx,ny,nz,Vnode,vrms)
!$ USE  OMP_LIB 
IMPLICIT NONE
INTEGER                            , INTENT(IN) :: nx,ny,nz
REAL, DIMENSION(1:3,1:nx,1:ny,1:nz), INTENT(IN) :: Vnode
REAL                               , INTENT(OUT):: vrms
INTEGER                                         :: ix,iy,iz

vrms = 0.0
!$OMP PARALLEL PRIVATE(ix,iy,iz)
!$OMP DO REDUCTION (+: vrms)
DO iz=1,nz
  DO iy=1,ny
    DO ix=1,nx
      vrms = vrms + Vnode(1,ix,iy,iz)**2+Vnode(2,ix,iy,iz)**2+Vnode(3,ix,iy,iz)**2
   ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL 

vrms = SQRT(vrms/FLOAT(nx*ny*nz))

END SUBROUTINE calc_vrms

! ========================================================================================================================= !
SUBROUTINE calc_mean(nx,ny,nz,a,amean)
!$ USE  OMP_LIB 
IMPLICIT NONE
INTEGER                        , INTENT(IN) :: nx,ny,nz
REAL, DIMENSION(1:nx,1:ny,1:nz), INTENT(IN) :: a 
REAL                           , INTENT(OUT):: amean 
INTEGER                                     :: ix,iy,iz

amean = 0.0
!$OMP PARALLEL PRIVATE(ix,iy,iz)
!$OMP DO REDUCTION (+: amean)
DO iz=1,nz
  DO iy=1,ny
    DO ix=1,nx
      amean = amean + a(ix,iy,iz)
   ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL 

amean = amean/FLOAT(nx*ny*nz)

END SUBROUTINE calc_mean
! ========================================================================================================================= !

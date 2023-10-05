# STB_CORRECTION
STB_correction is a software that performs corrections on a 3D velocity field in a Cartesian domain  
to apply rigid boundary conditions on all surfaces under a divergence-free constraint.

! -----------------------------------------!
! Author: Henri Samuel, 31/08/2023         !
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
! The STB algorithm used by DaVis did not have the capability of handling these issues, which motivated the implementation of STB_correction. 
!
! The STB_correction code proceeds as follows: 
! (1) Reads input file(s) to be corrected for velocities at the boundaries.
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
!  DO ix=1,nx_orig+1
!    DO iy=1,ny_orig+1
!      DO iz=1,nz_orig+1
!        READ(11,*) Vnode_orig(1:3,ix,iy,iz)
!      ENDDO
!    ENDDO
!  ENDDO
!
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

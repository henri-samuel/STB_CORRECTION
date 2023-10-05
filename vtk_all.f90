! --------------------------------------!
! Copyright (c) 2023 Henri Samuel       !
! All rights reserved.                  !
! --------------------------------------!
! ======================================================================= !

SUBROUTINE vtk_header3D(ifile,nx,ny,nz,dx,dy,dz,binary)
IMPLICIT NONE
INTEGER, INTENT(IN)               :: ifile,nz,ny,nx
REAL   , INTENT(IN)               :: dx, dy,dz
LOGICAL, INTENT(IN)               :: binary
LOGICAL                           ::  struct_grid
CHARACTER(LEN=50)                 :: cs
REAL*4, DIMENSION(:), ALLOCATABLE :: x,y,z
INTEGER                           :: ix,iy,iz
struct_grid =.true.
! The total number of characters in the header is 343 (just add up the numbers in the format)

! This number is needed when reading the file in binary format

IF (binary) THEN 
  WRITE(ifile) "# vtk DataFile Version 3.0"//char(10)     ! 26 characters
  WRITE(ifile) "STREAMV"//char(10)                        !  6 characters
  WRITE(ifile) "BINARY"//char(10)                         !  5 characters
  IF (struct_grid) THEN
    WRITE(ifile) "DATASET STRUCTURED_POINTS"//char(10)    ! 
    WRITE(cs,FMT='(A10,3(i5,2x))')   "DIMENSIONS",nx, ny, nz   ; WRITE(ifile) cs//char(10)
    WRITE(cs,FMT='(A10,3(f8.3,2x))') "ORIGIN", 0.,0.,0.        ; WRITE(ifile) cs//char(10)
    WRITE(cs,FMT='(A10,3(f8.3,2x))') "SPACING", dx,dy,dz       ; WRITE(ifile) cs//char(10)
    WRITE(cs,FMT='(A10,i10)')     "POINT_DATA", (nx)*(ny)*(nz) ; WRITE(ifile) cs//char(10)
  ELSE 
    ALLOCATE(x(1:nx),y(1:ny),z(1:nz))
    FORALL (ix=1:nx) x(ix) = 0.0+ FLOAT(ix-1)*dx
    FORALL (iy=1:ny) y(iy) = 0.0+ FLOAT(iy-1)*dy
    FORALL (iz=1:nz) z(iz) = 0.0+ FLOAT(iz-1)*dz
    WRITE(ifile) "DATASET RECTILINEAR_GRID "//char(10)
    WRITE(cs,fmt='(A10,3I5)') "DIMENSIONS",nx,ny,nz          ; WRITE(ifile) cs//char(10)
    WRITE(cs,fmt='(A13,I6,A6)') "X_COORDINATES",nx," float"  ; WRITE(ifile)           cs//char(10) ; WRITE(ifile) x(:)
    WRITE(cs,fmt='(A13,I6,A6)') "Y_COORDINATES",ny," float"  ; WRITE(ifile) char(10)//cs//char(10) ; WRITE(ifile) y(:)
    WRITE(cs,fmt='(A13,I6,A6)') "Z_COORDINATES",nz," float"  ; WRITE(ifile) char(10)//cs//char(10) ; WRITE(ifile) z(:)
    DEALLOCATE(x,y,z)
    WRITE(ifile) char(10)
  ENDIF
ELSE
  CALL vtk_header_ascii3D(ifile,nx,ny,nz,dx,dy,dz)
ENDIF

CALL FLUSH(ifile)  ! Seems important to avoid overwriting on needed blank spaces 
END SUBROUTINE vtk_header3D

! ============================================================================================================ !
SUBROUTINE vtk_header_ascii3D(ifile,nx,ny,nz,dx,dy,dz)
IMPLICIT NONE
INTEGER, INTENT(IN) :: ifile,nz,ny,nx
REAL   , INTENT(IN) :: dx, dy,dz

WRITE(ifile,1)
WRITE(ifile,2) nx, ny, nz 
WRITE(ifile,3) 0.,0.,0.
WRITE(ifile,4) dx,dy,dz    
WRITE(ifile,5) (nx)*(ny)*(nz)

1       FORMAT('# vtk DataFile Version 2.0',/,'STREAMV',/,'ASCII',/,'DATASET STRUCTURED_POINTS')
2       FORMAT('DIMENSIONS',1x,3(i7,1x))
3       FORMAT('ORIGIN',1x,3(f5.3,1x))
4       FORMAT('SPACING',1x,3(f12.8,1x))
5       FORMAT('POINT_DATA',1x,i15)

END SUBROUTINE vtk_header_ascii3D

! ======================================================================= !

SUBROUTINE vtk_header(ifile,nx,ny,nz,arx,ary,arz)
IMPLICIT NONE
INTEGER, INTENT(IN) :: ifile,nz,ny,nx
REAL   , INTENT(IN) :: arx, ary,arz
REAL                :: dx,dy,dz

dx=  arx/FLOAT(nx-1)
dy = ary/FLOAT(ny-1); IF (ny==1) dy=0.0001 
dz = arz/FLOAT(nz-1)

WRITE(ifile,1)
WRITE(ifile,2) nx, nz, ny ! For fake 2D
WRITE(ifile,3) INT(arx), INT(ary), int(arz) ! FAKE 2D
WRITE(ifile,4) 0.,0.,0.
WRITE(ifile,5) dx,dz,dy    ! For fake 2D
WRITE(ifile,6) (nx)*(ny)*(nz)

1       FORMAT('# vtk DataFile Version 2.0',/,'STREAMV',/,'ASCII',/,'DATASET STRUCTURED_POINTS')
2       FORMAT('DIMENSIONS',1x,3(i7,1x))
3       FORMAT('ASPECT_RATIO',1x,3(i2,1x))
4       FORMAT('ORIGIN',1x,3(f5.3,1x))
5       FORMAT('SPACING',1x,3(f12.8,1x))
6       FORMAT('POINT_DATA',1x,i15)

END subroutine vtk_header

! ======================================================================= !

SUBROUTINE vtk_scalar(nx,ny,nz,nstep,arx,ary,arz,a,fieldname)
IMPLICIT NONE
INTEGER                        , INTENT(IN) :: nx,ny,nz,nstep
REAL                           , INTENT(IN) :: arx, ary,arz
REAL, DIMENSION(1:nx,1:ny,1:nz), INTENT(IN) :: a
CHARACTER(LEN=*)               , INTENT(IN) :: fieldname
CHARACTER(LEN=50)                           :: namefile
INTEGER                                     :: iz,iy,ix
REAL                                        :: dx,dy,dz
CHARACTER(LEN=20), EXTERNAL                 :: str

namefile='VTK/'//TRIM(fieldname)//'_'//TRIM(str(nstep))//'.vtk'

WRITE(*,*) '==> Writing in vtk file: ',TRIM(namefile)
OPEN(41,FILE=TRIM(namefile))

dx=  arx/FLOAT(nx-1)
dy = ary/FLOAT(ny-1); IF (ny==1) dy=0.0001 
dz = arz/FLOAT(nz-1)

WRITE(41,1)
WRITE(41,2) nx, ny, nz ! For fake 2D
WRITE(41,3) INT(arx), INT(ary), INT(arz) 
WRITE(41,4) 0.,0.,0.
WRITE(41,5) dx,dz,dy    ! For fake 2D
WRITE(41,6) (nx)*(ny)*(nz)
WRITE(41,7) 'Test'
WRITE(41,8)
     
DO iz=1,nz
  DO iy=1,ny
    DO ix=1,nx
      WRITE(41,*)  a(ix,iy,iz)
!     if (abs(a(ix,iy,iz)) > 1e-4)  then
!     write(41,*) a(ix,iy,iz) 
!     else
!     write(41,*) 0.0 
!     endif
    ENDDO
  ENDDO
ENDDO

CLOSE(41)

1       FORMAT('# vtk DataFile Version 2.0',/,'TEST',/,'ASCII',/,'DATASET STRUCTURED_POINTS')
2       FORMAT('DIMENSIONS',1x,3(i7,1x))
3       FORMAT('ASPECT_RATIO',1x,3(i2,1x))
4       FORMAT('ORIGIN',1x,3(f5.3,1x))
5       FORMAT('SPACING',1x,3(f12.8,1x))
6       FORMAT('POINT_DATA',1x,i15)
7       FORMAT('SCALARS ',a5,' float 1')
8       FORMAT('LOOKUP_TABLE default')
9       FORMAT(f10.2)

END SUBROUTINE vtk_scalar
! ====================================================================================================== !

SUBROUTINE vtk_vector(nx,ny,nz,arx,ary,arz,a,nstep,fieldname)
IMPLICIT NONE
INTEGER, INTENT(IN)                            :: nz,ny,nx,nstep
REAL, DIMENSION(3,1:nx,1:ny,1:nz), INTENT(IN)  :: a
REAL, INTENT(IN)                               :: arx,ary,arz
CHARACTER(LEN=50)                              :: namefile
INTEGER                                        :: iz, iy, ix, l
REAL                                           :: dx,dy,dz
CHARACTER(LEN=*)                               :: fieldname
CHARACTER(LEN=20), EXTERNAL                   :: str

namefile='VTK/'//TRIM(fieldname)//'_'//TRIM(str(nstep))//'.vtk'

WRITE(*,*) '==> Writing in vtk file: ',TRIM(namefile)
OPEN(41,file=TRIM(namefile))

dx=  arx/FLOAT(nx-1)
dy = ary/FLOAT(ny-1); IF (ny==1) dy=0.0001
dz = arz/FLOAT(nz-1)

WRITE(41,1)
WRITE(41,2) nx, ny, nz
WRITE(41,3) INT(arx), INT(ary), INT(arz) 
WRITE(41,4) 0.,0.,0.
WRITE(41,5) dx,dy,dz
WRITE(41,6) (nx)*(ny)*(nz)
WRITE(41,7)

DO iz=1,nz
  DO iy=1,ny
    DO ix=1,nx
    WRITE(41,*) (a(l,ix,iy,iz),l=1,3)!, 0.0
    ENDDO
  ENDDO
ENDDO
CLOSE(41)

1       FORMAT('# vtk DataFile Version 2.0',/, &
              'TEST',/,'ASCII',/,'DATASET STRUCTURED_POINTS')
2       FORMAT('DIMENSIONS',1x,3(i7,1x))
3       FORMAT('ASPECT_RATIO',1x,3(i2,1x))
4       FORMAT('ORIGIN',1x,3(f5.3,1x))
5       FORMAT('SPACING',1x,3(f12.8,1x))
6       FORMAT('POINT_DATA',1x,i15)
7       FORMAT('VECTORS Velocity float')
8       FORMAT(1pe12.5,2x,1pe12.5,2x,1pe12.5)

END SUBROUTINE vtk_vector
! ========================================================================================= !

SUBROUTINE add_vtk_scalar(ifile,nx,ny,nz,a,fieldname,binary)
IMPLICIT NONE
INTEGER                        , INTENT(IN) :: ifile,nz,ny,nx
REAL, DIMENSION(1:nx,1:ny,1:nz), INTENT(IN) :: a
LOGICAL                        , INTENT(IN) :: binary
INTEGER                                     :: iz, iy, ix
CHARACTER(LEN=*)                            :: fieldname
CHARACTER(LEN=50)                           :: cs
REAL*4, DIMENSION(1:nx,1:ny,1:nz)           ::asngl 

IF (binary) THEN
  FORALL (ix=1:nx,iy=1:ny,iz=1:nz) asngl(ix,iy,iz) = SNGL(a(ix,iy,iz))   ! NB : This is mandatory 
  WRITE(cs,FMT='(A10,1x,A15,1x,A10)') "SCALARS",fieldname, "float 1"; WRITE(ifile) cs//char(10)
  WRITE(ifile) "LOOKUP_TABLE default"//char(10)
  WRITE(ifile) asngl 
ELSE
  CALL add_vtk_scalar_ascii(ifile,nx,ny,nz,a,fieldname)
ENDIF

END SUBROUTINE add_vtk_scalar

! ========================================================================================= !

SUBROUTINE add_vtk_scalar_ascii(ifile,nx,ny,nz,a,fieldname)
IMPLICIT NONE
INTEGER                        , INTENT(IN)  :: ifile,nz,ny,nx
REAL, DIMENSION(1:nx,1:ny,1:nz), INTENT(IN)  :: a
INTEGER                                      :: iz, iy, ix
CHARACTER(LEN=5)                             :: fieldname


WRITE(ifile,'(A10,1x,A5,1x,A10)') "SCALARS",fieldname, "float 1"
WRITE(ifile,'(A30)') "LOOKUP_TABLE default"
DO iz=1,nz
  DO iy=1,ny
    DO ix=1,nx
    WRITE(ifile,'(f10.4)')  a(ix,iy,iz)
    ENDDO
  ENDDO
ENDDO

END SUBROUTINE add_vtk_scalar_ascii

! ========================================================================================== !

SUBROUTINE add_vtk_vector(ifile,nx,ny,nz,a,fieldname,binary)
IMPLICIT NONE
INTEGER                          , INTENT(IN) :: ifile,nx,ny,nz
REAL, DIMENSION(3,1:nx,1:ny,1:nz), INTENT(IN) :: a
CHARACTER(LEN=*)                 , INTENT(IN) :: fieldname
LOGICAL                          , INTENT(IN) :: binary
CHARACTER(LEN=50)                             :: cs
INTEGER                                       :: ix, iy, iz, iv
REAL*4,DIMENSION(1:3,1:nx,1:ny,1:nz)          :: asngl


IF (binary) THEN
  FORALL (iv=1:3,ix=1:nx,iy=1:ny,iz=1:nz) asngl(iv,ix,iy,iz) = SNGL(a(iv,ix,iy,iz))
  WRITE(cs,FMT='(A10,1x,A15,1x,A10)') "VECTORS",fieldname, "float 1"; WRITE(ifile) cs//char(10)
!  WRITE(cs,FMT='(A10,1x,A15,1x,A10)') "VECTORS",fieldname, "float "; WRITE(ifile) cs//char(10)
  WRITE(ifile) asngl
ELSE
  CALL add_vtk_vector_ascii_3D(ifile,nx,ny,nz,a,fieldname)
ENDIF

END subroutine add_vtk_vector

! ======================================================================== !

SUBROUTINE add_vtk_vector_ascii_3D(ifile,nx,ny,nz,a,fieldname)
IMPLICIT NONE
INTEGER                          , INTENT(IN)  :: ifile,nz,ny,nx
REAL, DIMENSION(3,1:nx,1:ny,1:nz), INTENT(IN)  :: a
INTEGER                                        :: iz, iy, ix, l
CHARACTER(LEN=*)                               :: fieldname

WRITE(ifile,7)

DO iz=1,nz
  DO iy=1,ny
    DO ix=1,nx
    WRITE(ifile,*) (a(l,ix,iy,iz),l=1,3)
    ENDDO
  ENDDO
ENDDO

7       FORMAT('VECTORS Velocity float')
8       FORMAT(1pe12.5,2x,1pe12.5,2x,1pe12.5)

END subroutine add_vtk_vector_ascii_3D

! ======================================================================== !

SUBROUTINE add_vtk_vector_ascii(ifile,nx,ny,nz,a,fieldname)
IMPLICIT NONE
INTEGER, INTENT(IN)                            :: ifile,nz,ny,nx
REAL, DIMENSION(2,1:nx,1:ny,1:nz), INTENT(IN)  :: a
INTEGER                                        :: iz, iy, ix, l
CHARACTER(LEN=5)                               :: fieldname

WRITE(ifile,7) 
      
DO iz=1,nz
  DO iy=1,ny
    DO ix=1,nx
    WRITE(ifile,*) (a(l,ix,iy,iz),l=1,2), 0.0
    ENDDO
  ENDDO
ENDDO

7       FORMAT('VECTORS Velocity float')
8       FORMAT(1pe12.5,2x,1pe12.5,2x,1pe12.5)

END subroutine add_vtk_vector_ascii
! ============================================================================================================ !

SUBROUTINE vtk_scalar_asciibin(nx,ny,nz,dx,dy,dz,iframe,fieldname,binary,a)

IMPLICIT NONE
INTEGER                              , INTENT(IN) :: nx,ny,nz
REAL                                 , INTENT(IN) :: dx,dy,dz
INTEGER                              , INTENT(IN) :: iframe
CHARACTER(LEN=*)                     , INTENT(IN) :: fieldname
LOGICAL                              , INTENT(IN) :: binary
REAL, DIMENSION(1:nx,1:ny,1:nz), INTENT(IN) :: a
INTEGER, PARAMETER                                :: ID = 777

CALL open_file(     ID     ,iframe,binary,TRIM(fieldname)//'_node')   
CALL vtk_header3D(  ID,nx,ny,nz,dx,dy,dz,binary)
CALL add_vtk_scalar(ID,nx,ny,nz,a,fieldname,binary)

CLOSE(ID)

END SUBROUTINE vtk_scalar_asciibin

! =============================================================================================== !
SUBROUTINE vtk_vector_asciibin(nx,ny,nz,dx,dy,dz,iframe,fieldname,binary,v)

IMPLICIT NONE
INTEGER                              , INTENT(IN) :: nx,ny,nz
REAL                                 , INTENT(IN) :: dx,dy,dz
INTEGER                              , INTENT(IN) :: iframe
CHARACTER(LEN=*)                     , INTENT(IN) :: fieldname
LOGICAL                              , INTENT(IN) :: binary
REAL, DIMENSION(1:3,1:nx,1:ny,1:nz)  , INTENT(IN) :: v
INTEGER, PARAMETER                                :: ID = 888 

CALL open_file(     ID     ,iframe,binary,TRIM(fieldname)//'_node')
CALL vtk_header3D(  ID,nx,ny,nz,dx,dy,dz,binary)
CALL add_vtk_vector(ID,nx,ny,nz,v,fieldname,binary)

CLOSE(ID)

END SUBROUTINE vtk_vector_asciibin


! ======================================================================= !

SUBROUTINE open_file(ifile,iplot,binary,point_type)
INTEGER         , INTENT(IN) :: ifile, iplot
LOGICAL         , INTENT(IN) :: binary
CHARACTER(LEN=*), INTENT(IN) :: point_type
CHARACTER(LEN=50)            :: filename
CHARACTER(LEN=20), EXTERNAL  :: str

IF (iplot >=0) THEN
  filename='VTK/All_'//TRIM(point_type)//'_'//TRIM(str(iplot))//'.vtk'
ELSE
  filename='VTK/All_'//TRIM(point_type)//'.vtk'
ENDIF


WRITE(*,*) '==> Writing in vtk file: ',TRIM(filename)

IF (binary) THEN
 OPEN(ifile,file=TRIM(filename),access='stream',status='unknown',convert='big_endian')
ELSE
  OPEN(ifile,file=TRIM(filename))
ENDIF

END SUBROUTINE open_file

! ======================================================================================================================== !


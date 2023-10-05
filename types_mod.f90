! --------------------------------------!
! Copyright (c) 2023 Henri Samuel       !
! All rights reserved.                  !
! --------------------------------------!

MODULE types_mod
IMPLICIT NONE

! ------------------------------------ Boundary conditions --------------------------------------------------------- !

TYPE scalbcs
  REAL              :: Botval,Topval,Westval,Eastval,Southval,Northval
  CHARACTER(LEN=15) :: West,East,South,North,Top,Bot 
END TYPE scalbcs

! ---------------------------------------- field arrays ------------------------------------------------------------- !

TYPE all_sol_array
  REAL, DIMENSION(:,:,:), POINTER :: cont
END TYPE all_sol_array

TYPE grid_type
  INTEGER  :: nx,ny,nz
  REAL     :: dx,dy,dz
END TYPE grid_type

TYPE mgrid_all
  TYPE(grid_type)      :: grid
  TYPE(all_sol_array)  :: val,val_prev,rhs,res
END TYPE mgrid_all

END MODULE types_mod

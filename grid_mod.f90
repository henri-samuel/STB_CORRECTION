! --------------------------------------!
! Copyright (c) 2023 Henri Samuel       !
! All rights reserved.                  !
! --------------------------------------!

MODULE grid_mod
USE types_mod, ONLY: grid_specs
IMPLICIT NONE
INTEGER             , SAVE :: nx,ny,nz
REAL                , SAVE :: dx,dy,dz,arx,ary,arz,xmin,xmax,ymin,ymax,zmin,zmax
TYPE(grid_specs)           :: grid

END MODULE grid_mod

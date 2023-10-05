! --------------------------------------!
! Copyright (c) 2023 Henri Samuel       !
! All rights reserved.                  !
! --------------------------------------!

MODULE bcs_mod
USE types_mod
IMPLICIT NONE
CHARACTER(LEN=32), SAVE :: WestVBC_mode,EastVBC_mode,SouthVBC_mode,NorthVBC_mode,BotVBC_mode,TopVBC_mode
TYPE(scalbcs)    , SAVE :: Vx_BC,Vz_BC,VxBC,VyBC,VzBC

END MODULE bcs_mod

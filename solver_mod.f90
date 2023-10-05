! --------------------------------------!
! Copyright (c) 2023 Henri Samuel       !
! All rights reserved.                  !
! --------------------------------------!
MODULE solver_mod
IMPLICIT NONE
INTEGER, SAVE :: nsteps,iframestep
INTEGER, SAVE :: ngmax,nvee,npre,npost,nmaxit
REAL   , SAVE :: errmax_multi
REAL   , SAVE :: time
REAL   , SAVE :: w_cont,rincr,niter_finest
END MODULE solver_mod

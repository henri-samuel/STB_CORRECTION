###########################################################################
##############              Author: Henri Samuel               ############
##############                      CNRS                       ############
###########################################################################

########################## Object and module files  ###################################

ODIR = ./OBJ
MDIR = $(ODIR)

########################################################################################
################################ COMPILERS #############################################
########################################################################################

#NB mpifor --showme   shows explicitly all the paths and options
#F90    =  mpif90   
#F90    =  mpifort   
F90  = gfortran
FOPT   = -O3 -fopenmp -funroll-loops  -finit-local-zero -finit-real=zero -finit-logical=false  #-framework veclib # -ffast-math  -Ofast
FCHECK = -fbounds-check    #    -fcheck=all -finit-real=snan -Wuninitialized -Wsurprising   -Werror
FFLAGS =      -fdefault-real-8 -fdefault-double-8 -I$(MDIR) $(FOPT)  $(FCHECK)
FFLAGSLIGHT = -fdefault-real-8 -fdefault-double-8 -I$(MDIR) $(FOPT)

#########################################################################################

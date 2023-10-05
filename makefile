.SUFFIXES: .f90 .o .mod .f .c .f90.o .f.o .c.o
##########################################################################################
##############           Makefile for SBT_correction CODE, Version 1.0.0      ############      
##############                      Author: Henri Samuel                      ############
##############                             CNRS                               ############ 
##############                             2023                               ############
############## This freeware is distributed under the GNU GPLv3 Licence terms ############
##########################################################################################

#include MAKE_INC/Makefile_LINUX_STELLA_PAR_INTEL.inc
include MAKE_INC/Makefile_OSX_PAR_GFORTRAN.inc  

FFLAGS = -fdefault-real-8 -fdefault-double-8 $(FOPT)   -I$(MDIR)
#################################### Subroutines ##################################

SUBS =  $(addprefix OBJ/, $(addsuffix .o, apply_bcs arrays interpolations main \
         multi  HG_correction vtk_all  ))

##########################################################################################

MODS = $(addprefix $(ODIR)/, $(addsuffix .o, bcs_mod \
         solver_mod types_mod  ))

MMODS = $(MODS:%.o=%.mod)

##########################################################################################
EXEC    =stb_correction.exe

stb   : $(SUBS) 
	$(F90) $(FFLAGS)  -o $(EXEC)   $(SUBS) $(MODS)  $(LIBS)   

$(ODIR)/apply_bcs.o : apply_bcs.f90 $(addsuffix .mod, $(addprefix $(MDIR)/, bcs_mod  solver_mod types_mod ))
	$(F90) $(FFLAGS) -c -o $@  $<  
$(ODIR)/arrays.o : arrays.f90 $(addsuffix .mod, $(addprefix $(MDIR)/,    ))
	$(F90) $(FFLAGS) -c -o $@  $<  
$(ODIR)/interpolations.o : interpolations.f90  $(addsuffix .mod, $(addprefix $(MDIR)/, types_mod ))
	$(F90) $(FFLAGS) -c -o $@  $<  
$(ODIR)/main.o: main.f90  $(MMODS) 
	$(F90) $(FFLAGS) -c -o $@  $<  
$(ODIR)/multi.o: multi.f90  $(addsuffix .mod, $(addprefix $(MDIR)/, bcs_mod solver_mod  types_mod ))
	$(F90) $(FFLAGS) -c -o $@  $<
$(ODIR)/HG_correction.o : HG_correction.f90 $(addsuffix .mod, $(addprefix $(MDIR)/, bcs_mod types_mod ))
	$(F90) $(FFLAGS) -c -o $@  $<
$(ODIR)/vtk_all.o: vtk_all.f90
	$(F90) $(FFLAGS) -c -o $@  $<  

####################### MODULES ##############################
$(MDIR)/bcs_mod.mod :: bcs_mod.f90 $(addprefix $(MDIR)/, types_mod.mod )
	$(F90) $(FFLAGS) -c bcs_mod.f90 -o $(addsuffix .o ,$(basename $@)) -J $(MDIR) 
$(MDIR)/solver_mod.mod :: solver_mod.f90
	$(F90) $(FFLAGS) -c solver_mod.f90  -o $(addsuffix .o ,$(basename $@)) -J $(MDIR)  
$(MDIR)/types_mod.mod :: types_mod.f90
	$(F90) $(FFLAGS) $(INC) -c  types_mod.f90  -o $(addsuffix .o ,$(basename $@)) -J $(MDIR)  
#################### clean-up commands ########################

clean :
	rm $(ODIR)/*.o $(MDIR)/*.mod

cleanexe :
	rm *.exe	

cleandat :
	rm ./DAT/*.dat

cleanvtk :
	rm ./VTK/*.vtk

cleanall : clean cleanexe  cleandat cleanvtk


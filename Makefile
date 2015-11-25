#*************************** EULER ***********************************
#
#
#
#
#
### intel fortran compiler flags 
HOSTNAME=$(shell hostname)
COMP=mpif90
ifeq ($(strip $(HOSTNAME)),himura)

      COMP=/usr/local/mpiintel/bin/mpif90
endif
FLAGS   = 
#FLAGS  =-autodouble
LIBS   = 
#EXTRA_FLAGS=-fpe0       
OBJSM  = vpm.o pmlib.o pmproject.o yaps.o
OBJSHP   = testmod.o mkl_dfti.o main_pm.o pmlib.o pmproject.o libpm.o yaps.o vpm.o vpm_mpi.o vpm_time.o vpm_gcalc.o vpm_remesh.o test_applic.o test.o
OBJS     = $(OBJSHP) 
#OBJSPS   = yaps.o
EXENAMEIN = vpm
MODPATH = ./modules
SRCHP   = ./source
     
     
DEBUG ?=0
ifeq ($(debug), 1)
     FLAGS = -C -fpe0 -traceback -mkl -O3  #-openmp
     EXTRA_FLAGS= -warn interfaces -warn general                   
     EXENAME=$(EXENAMEIN)_debug
     PATHOBJHP = ./debug
     VPATH = ./debug
else
     PATHOBJHP = ./release
     EXENAME=$(EXENAMEIN)
     FLAGS =-mkl -O3 -openmp -warn all -diag-disable id,7712 
     EXTRA_FLAGS= #-warn interfaces -warn general
     VPATH = ./release
endif
#
#
$(EXENAME): $(OBJSHP) 
	 $(COMP) $(patsubst %.o,$(PATHOBJHP)/%.o, $(OBJSHP))  -module $(MODPATH)  $(FLAGS) $(EXTRA_FLAGS) $(LIBS)  -o $(EXENAME)

#

main_pm.o: $(SRCHP)/main_pm.f90               Makefile
	$(COMP)  $(FLAGS) $(EXTRA_FLAGS) -c -module  $(MODPATH) $< -o $(PATHOBJHP)/$@
#
test.o: $(SRCHP)/test.f90 	
	$(COMP)  $(FLAGS) $(EXTRA_FLAGS)  -c -module $(MODPATH) $< -o $(PATHOBJHP)/$@
#
testmod.o: $(SRCHP)/testmod.f90 	
	$(COMP)  $(FLAGS) $(EXTRA_FLAGS)  -c -module $(MODPATH) $< -o $(PATHOBJHP)/$@
#
test_applic.o: $(SRCHP)/test_applic.f90 	
	$(COMP)  $(FLAGS) $(EXTRA_FLAGS)  -c -module $(MODPATH) $< -o $(PATHOBJHP)/$@
#
vpm.o: $(SRCHP)/vpm.f90 $(SRCHP)/yaps.f90 $(SRCHP)/vpm_mpi.f90 $(SRCHP)/vpm_gcalc.f90 Makefile
	$(COMP)  $(FLAGS) $(EXTRA_FLAGS)  -c -module $(MODPATH) $< -o $(PATHOBJHP)/$@
#
vpm_mpi.o: $(SRCHP)/vpm_mpi.f90
	$(COMP)  $(FLAGS) $(EXTRA_FLAGS)  -c -module $(MODPATH) $< -o $(PATHOBJHP)/$@
#
vpm_remesh.o: $(SRCHP)/vpm_remesh.f90
	$(COMP)  $(FLAGS) $(EXTRA_FLAGS)  -c -module $(MODPATH) $< -o $(PATHOBJHP)/$@    
#
vpm_gcalc.o: $(SRCHP)/vpm_gcalc.f90
	$(COMP)  $(FLAGS) $(EXTRA_FLAGS)  -c -module $(MODPATH) $< -o $(PATHOBJHP)/$@    
#
vpm_time.o: $(SRCHP)/vpm_time.f90
	$(COMP)  $(FLAGS) $(EXTRA_FLAGS)  -c -module $(MODPATH) $< -o $(PATHOBJHP)/$@    
#
yaps.o: $(SRCHP)/yaps.f90  $(SRCHP)/yaps2d.f90 $(SRCHP)/yaps3d.f90 Makefile
	$(COMP)  $(FLAGS) $(EXTRA_FLAGS)  -c -module $(MODPATH) $< -o $(PATHOBJHP)/$@
#
pmlib.o: $(SRCHP)/pmlib.f90  $(SRCHP)/pmbound.f90 $(SRCHP)/pinfdomain.f90 $(SRCHP)/pmsolve.f90  Makefile
	$(COMP)  $(FLAGS) $(EXTRA_FLAGS) -c -module  $(MODPATH) $< -o $(PATHOBJHP)/$@
#
pmproject.o: $(SRCHP)/pmproject.f90 $(SRCHP)/main_pm.f90  Makefile
	$(COMP)  $(FLAGS) $(EXTRA_FLAGS) -c -module  $(MODPATH) $< -o $(PATHOBJHP)/$@
#
mkl_dfti.o: $(SRCHP)/mkl_dfti.f90 #              Makefile
	$(COMP)  $(FLAGS) $(EXTRA_FLAGS) -c -module  $(MODPATH) $< -o $(PATHOBJHP)/$@

#
libpm.o: $(SRCHP)/libpm.f90                            
	$(COMP)  -c -module  $(MODPATH) $< -o $(PATHOBJHP)/$@

clean :
	rm -rf $(EXENAME) $(PATHOBJNS)/*.o $(PATHOBJNS)/*.mod *.o $(MODPATH)/*.mod $(PATHOBJPS)/*.mod $(PATHOBJHP)/*.mod $(PATHOBJHP)/*.o
install :
	cp $(EXENAME) /home/papis/bin

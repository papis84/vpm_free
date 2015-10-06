#*************************** EULER ***********************************
#
#
#
#
#
### intel fortran compiler flags 
COMP= /usr/local/mpiintel/bin/mpif90
FLAGS   = 
#FLAGS  =-autodouble
LIBS   = 
#EXTRA_FLAGS=-fpe0       
OBJSMOD  = vpm.o pmlib.o pmproject.o yaps.o
OBJSHP   =  main_pm.o testmod.o mkl_dfti.o pmlib.o pmproject.o libpm.o yaps.o vpm.o vpm_mpi.o vpm_time.o vpm_gcalc.o vpm_remesh.o test.o
OBJS     = $(OBJSHP) 
#OBJSPS   = yaps.o
EXENAMEIN = vpm
MODPATH = ./modules
SRCHP   = ./source
     
     
DEBUG ?=0
ifeq ($(DEBUG), 1)
     FLAGS = -C -fpe0 -traceback -mkl -O3  #-openmp
     EXTRA_FLAGS= -warn interfaces -warn general                   
     EXENAME=$(EXENAMEIN)_debug
     PATHOBJHP = ./debug
else
     PATHOBJHP = ./release
     EXENAME=$(EXENAMEIN)
     FLAGS =-mkl -O3 -openmp -warn all -diag-disable id,7712 
     EXTRA_FLAGS= #-warn interfaces -warn general
endif
#
#
$(EXENAME): $(OBJS) 
	    $(COMP) -module $(MODPATH) $(patsubst %.o,$(PATHOBJHP)/%.o, $(OBJSHP)) $(FLAGS) $(EXTRA_FLAGS)  $(LIBS)  -o $(EXENAME)

#

#
test.o: $(SRCHP)/test.f90 	
	$(COMP)  $(FLAGS) $(EXTRA_FLAGS)  -c -module $(MODPATH) $< -o $(PATHOBJHP)/$@
#
testmod.o: $(SRCHP)/testmod.f90 	
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
main_pm.o: $(SRCHP)/main_pm.f90               Makefile
	$(COMP)  $(FLAGS) $(EXTRA_FLAGS) -c -module  $(MODPATH) $< -o $(PATHOBJHP)/$@
#
pmlib.o: $(SRCHP)/pmlib.f90  $(SRCHP)/pmbound.f90 $(SRCHP)/pinfdomain.f90 $(SRCHP)/pmsolve.f90  Makefile
	$(COMP)  $(FLAGS) $(EXTRA_FLAGS) -c -module  $(MODPATH) $< -o $(PATHOBJHP)/$@
#
pmproject.o: $(SRCHP)/pmproject.f90 $(SRCHP)/main_pm.f90  Makefile
	$(COMP)  $(FLAGS) $(EXTRA_FLAGS) -c  -module  $(MODPATH) $< -o $(PATHOBJHP)/$@
#
pmesh.o: $(SRCHP)/pmesh.f90                  $(SRCHP)/main_pm.f90  Makefile
	$(COMP)  $(FLAGS) $(EXTRA_FLAGS) -c -module  $(MODPATH) $< -o $(PATHOBJHP)/$@
#
pmbound.o: $(SRCHP)/pmbound.f90 $(SRCHP)/main_pm.f90  Makefile
	$(COMP)  $(FLAGS) $(EXTRA_FLAGS) -c -module  $(MODPATH) $< -o $(PATHOBJHP)/$@
#
pmsolve.o: $(SRCHP)/pmsolve.f90 $(SRCHP)/main_pm.f90  Makefile
	$(COMP)  $(FLAGS) $(EXTRA_FLAGS) -c  -module  $(MODPATH) $< -o $(PATHOBJHP)/$@
#
pmgcalc.o: $(SRCHP)/pmgcalc.f90    $(SRCHP)/main_pm.f90  Makefile
	$(COMP)  $(FLAGS) $(EXTRA_FLAGS) -c  -module  $(MODPATH) $< -o $(PATHOBJHP)/$@
#
pinfdomain.o: $(SRCHP)/pinfdomain.f90  $(SRCHP)/main_pm.f90  Makefile
	$(COMP)  $(FLAGS) $(EXTRA_FLAGS) -c -module  $(MODPATH) $< -o $(PATHOBJHP)/$@
#
libpm.o: $(SRCHP)/libpm.f90                            
	$(COMP)  -c -module  $(MODPATH) $< -o $(PATHOBJHP)/$@
#
mkl_dfti.o: $(SRCHP)/mkl_dfti.f90               Makefile
	$(COMP)  $(FLAGS) $(EXTRA_FLAGS) -c -module  $(MODPATH) $< -o $(PATHOBJHP)/$@


clean :
	rm -rf $(EXENAME) $(PATHOBJNS)/*.o $(PATHOBJNS)/*.mod *.o $(MODPATH)/*.mod $(PATHOBJPS)/*.mod $(PATHOBJHP)/*.mod $(PATHOBJHP)/*.o
install :
	cp $(EXENAME) /home/papis/bin

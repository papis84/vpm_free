#*************************** EULER ***********************************
#
#
#
#
#
### intel fortran compiler flags 
HOSTNAME  = $(shell hostname)
COMP      = mpif90
OBJSM     = vpm.o pmlib.o pmproject.o yaps.o
OBJSHP    = testmod.o mkl_dfti.o main_pm.o pmlib.o pmproject.o libpm.o yaps.o vpm.o vpm_mpi.o vpm_time.o vpm_gcalc.o vpm_remesh.o test_applic.o test.o
OBJS      = $(OBJSHP) 
EXENAMEIN = vpm
SRCHP     = ./source
     
     
dbg ?=0
ifeq ($(dbg), 1)
     FLAGS       = -C -fpe0 -traceback -mkl -O3  #-openmp
     EXTRA_FLAGS = -warn interfaces -warn general
     EXENAME     = $(EXENAMEIN)_debug
     PATHOBJHP   = $(SRCHP)/debug
else
     FLAGS       = -mkl -O3  #-warn all -diag-disable id,7712 
     EXTRA_FLAGS = #-warn interfaces -warn general
     EXENAME     = $(EXENAMEIN)
     PATHOBJHP   = $(SRCHP)/release
endif
     VPATH       = $(PATHOBJHP)
     MODPATH     = $(PATHOBJHP)


.PHONY        : all clean cleanall complete install prepare remake uninstall

complete      : prepare $(EXENAME) install

$(EXENAME)    : $(OBJSHP) 
		$(COMP) $(patsubst %.o,$(PATHOBJHP)/%.o, $(OBJSHP))  -module $(MODPATH)  $(FLAGS) $(EXTRA_FLAGS) -o $(EXENAME)



main_pm.o     : $(SRCHP)/main_pm.f90               Makefile
		$(COMP)  $(FLAGS) $(EXTRA_FLAGS) -c -module $(MODPATH) $< -o $(PATHOBJHP)/$@

test.o        : $(SRCHP)/test.f90 	
		$(COMP)  $(FLAGS) $(EXTRA_FLAGS) -c -module $(MODPATH) $< -o $(PATHOBJHP)/$@

testmod.o     : $(SRCHP)/testmod.f90 	
		$(COMP)  $(FLAGS) $(EXTRA_FLAGS) -c -module $(MODPATH) $< -o $(PATHOBJHP)/$@

test_applic.o : $(SRCHP)/test_applic.f90 	
		$(COMP)  $(FLAGS) $(EXTRA_FLAGS) -c -module $(MODPATH) $< -o $(PATHOBJHP)/$@

vpm.o         : $(SRCHP)/vpm.f90 $(SRCHP)/yaps.f90 $(SRCHP)/vpm_mpi.f90 $(SRCHP)/vpm_gcalc.f90 Makefile
		$(COMP)  $(FLAGS) $(EXTRA_FLAGS) -c -module $(MODPATH) $< -o $(PATHOBJHP)/$@

vpm_mpi.o     : $(SRCHP)/vpm_mpi.f90
		$(COMP)  $(FLAGS) $(EXTRA_FLAGS) -c -module $(MODPATH) $< -o $(PATHOBJHP)/$@

vpm_remesh.o  : $(SRCHP)/vpm_remesh.f90
		$(COMP)  $(FLAGS) $(EXTRA_FLAGS) -c -module $(MODPATH) $< -o $(PATHOBJHP)/$@    

vpm_gcalc.o   : $(SRCHP)/vpm_gcalc.f90
		$(COMP)  $(FLAGS) $(EXTRA_FLAGS) -c -module $(MODPATH) $< -o $(PATHOBJHP)/$@    

vpm_time.o    : $(SRCHP)/vpm_time.f90
		$(COMP)  $(FLAGS) $(EXTRA_FLAGS) -c -module $(MODPATH) $< -o $(PATHOBJHP)/$@    

yaps.o        : $(SRCHP)/yaps.f90  $(SRCHP)/yaps2d.f90 $(SRCHP)/yaps3d.f90 Makefile
		$(COMP)  $(FLAGS) $(EXTRA_FLAGS) -c -module $(MODPATH) $< -o $(PATHOBJHP)/$@

pmlib.o       : $(SRCHP)/pmlib.f90  $(SRCHP)/pmbound.f90 $(SRCHP)/pinfdomain.f90 $(SRCHP)/pmsolve.f90  Makefile
		$(COMP)  $(FLAGS) $(EXTRA_FLAGS) -c -module $(MODPATH) $< -o $(PATHOBJHP)/$@

pmproject.o   : $(SRCHP)/pmproject.f90 $(SRCHP)/main_pm.f90  Makefile
		$(COMP)  $(FLAGS) $(EXTRA_FLAGS) -c -module $(MODPATH) $< -o $(PATHOBJHP)/$@

mkl_dfti.o    : $(SRCHP)/mkl_dfti.f90 #              Makefile
		$(COMP)  $(FLAGS) $(EXTRA_FLAGS) -c -module $(MODPATH) $< -o $(PATHOBJHP)/$@

libpm.o       : $(SRCHP)/libpm.f90                            
		$(COMP)  -c -module  $(MODPATH) $< -o $(PATHOBJHP)/$@


clean         :
		@rm -r -f $(PATHOBJHP)                      $(EXENAME)
		@reset
cleanall      :
		@rm -r -f $(SRCHP)/debug $(SRCHP)/release* $(EXENAMEIN)*
		@reset
install       :
		@ln -s -f $(shell pwd)/$(EXENAME) ~/bin
prepare       :
		@mkdir -p $(PATHOBJHP)
		@mkdir -p ~/bin
remake        : clean complete
uninstall     :
		@rm ~/bin/$(EXENAMEIN)*
all           :
		make dbg=0
		make dbg=1


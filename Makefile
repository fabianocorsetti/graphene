include make.inc

.SUFFIXES:
.SUFFIXES: .x .f90 .F90

all : routines_lib MatrixSwitch_lib screening.x TB.x adatom_post.x

screening : clear_MS routines_lib screening.x

TB : routines_lib MatrixSwitch_lib TB.x

adatom_post : clear_MS clear_ROUTINES adatom_post.x

clear_MS :
	$(eval MS_INC := )
	$(eval MS_LIB := )

clear_ROUTINES :
	$(eval ROUTINES_INC := )
	$(eval ROUTINES_LIB := )

clean : clean_routines clean_MatrixSwitch
	rm -f *.x
	rm -f *.o
	rm -f *.mod

routines_lib :
	cd routines; \
	make

clean_routines :
	cd routines; \
	make clean

MatrixSwitch_lib :
	cd MatrixSwitch/src; \
	make

clean_MatrixSwitch :
	cd MatrixSwitch/src; \
	make clean

.F90.x : 
	$(FORTRAN) $(OPTS) $(FPPFLAGS) $(ROUTINES_INC) $(SPGLIB_INC) $(FFTW_INC) $(MS_INC) $(LINALG_INC) $< $(ROUTINES_LIB) $(SPGLIB_LIB) $(FFTW_LIB) $(MS_LIB) $(LINALG_LIB) -o $@

.f90.x : 
	$(FORTRAN) $(OPTS) $(ROUTINES_INC) $(SPGLIB_INC) $(FFTW_INC) $(MS_INC) $(LINALG_INC) $< $(ROUTINES_LIB) $(SPGLIB_LIB) $(FFTW_LIB) $(MS_LIB) $(LINALG_LIB) -o $@

include make.inc

.SUFFIXES:
.SUFFIXES: .o .a .f90 .F90

all : $(LIBNAME)

ALLOBJ = MatrixSwitch_ops.o MatrixSwitch_mm_multiply.o MatrixSwitch_m_add.o \
         MatrixSwitch_m_set.o MatrixSwitch_m_copy.o MatrixSwitch_m_register.o \
         MatrixSwitch.o MatrixSwitch_wrapper_params.o MatrixSwitch_wrapper.o

$(LIBNAME) : $(ALLOBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

install :
	if [ ! -d $(BUILDPATH) ]; then mkdir $(BUILDPATH); fi; \
	if [ ! -d $(BUILDPATH)/include ]; then mkdir $(BUILDPATH)/include; fi; \
	if [ ! -d $(BUILDPATH)/lib ]; then mkdir $(BUILDPATH)/lib; fi; \
	cp *.mod $(BUILDPATH)/include; \
	cp *.a $(BUILDPATH)/lib;

clean : 
	rm -f *.a
	rm -f *.o
	rm -f *.mod

.F90.o : 
	$(FORTRAN) $(OPTS) $(FPPFLAGS) -c $< -I$(PSPINCLUDEPATH) -o $@

.f90.o : 
	$(FORTRAN) $(OPTS) -c $< -I$(PSPINCLUDEPATH) -o $@

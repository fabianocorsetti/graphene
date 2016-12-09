SHELL     = /bin/sh
FORTRAN   = mpif90
OPTS      = -O3
FPPFLAGS  = -DCONV -DMPI -DLAP -DSLAP # do not modify
ARCH      = ar
ARCHFLAGS = cr
RANLIB    = ranlib

# Spglib
SPGLIB_INC = 
SPGLIB_LIB = -lsymspg

# FFTW3
FFTW_INC = 
FFTW_LIB = -lfftw3

# BLAS+LAPACK+ScaLAPACK
LINALG_INC = -I${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include
LINALG_LIB = ${MKLROOT}/lib/intel64/libmkl_blas95_lp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a ${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group ${MKLROOT}/lib/intel64/libmkl_blacs_openmpi_lp64.a -lpthread -lm

#do not modify below this line

MSLIBNAME = MatrixSwitch.a
MS_INC    = -I./MatrixSwitch/src/
MS_LIB    = ./MatrixSwitch/src/${MSLIBNAME}

ROUTINESLIBNAME = routines.a
ROUTINES_INC    = -I./routines/
ROUTINES_LIB    = ./routines/${ROUTINESLIBNAME}

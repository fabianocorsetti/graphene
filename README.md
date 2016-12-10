Graphene model screening + tight-binding
========================================

Author
------

Fabiano Corsetti, Imperial College London, UK

*   Email: fabiano.corsetti08 {at} imperial.ac.uk
*   Homepage: <http://www.cmth.ph.ic.ac.uk/people/f.corsetti/>

Description
-----------

This code compiles three separate executables, `screening.x`, `TB.x`, and
`adatom_post.x`. `screening.x` is a serial program which calculates the model
screened potential for Coulomb charges on graphene and the symmetries of the
supercell. `TB.x` is a parallel program for performing a tight-binding
simulation of the system. `adatom_post.x` is a serial program for post-
processing the output of the tight-binding simulation for the case of a single
adatom at the origin of an n x n supercell.

Installation
------------

1.  Modify `make.inc` to suit your needs. In particular, you will need to
    provide a Fortran compiler (MPI-enabled for `TB.x`) and the following
    libraries:

    *   Spglib for crystal symmetries (only needed for `screening.x`)
    *   FFTW3 for FFTs (only needed for `screening.x`)
    *   BLAS+LAPACK+ScaLAPACK for linear algebra (only needed for `TB.x`)

    The example provided will work on the cmth cluster after loading the
    following modules:

    *   intel-suite/2015.3.187
    *   openmpi/1.10.2
    *   mkl/2015.3.187
    *   fftw/3.3.4
    *   spglib/1.9.8-20-gb63ea50

2.  Type `make` to compile both codes. Alternatively use `make screening` to
    only compile `screening.x`, `make TB` to only compile `TB.x`, or
    `make adatom_post` to only compile `adatom_post.x`.

Usage
-----

*   `screening.x`:

    *   Run as: `./screening.x`
    *   Input: The program reads a single text file called `screening.in`. An
        example is provided in the `example` directory. All variables must be
        present in the given order. An explanation for each variable is given
        in the example. The input file specifies the dimensions of the
        supercell, the impurity Coulomb charges providing the bare external
        potential, the specifications of the screening model to use, and the
        k-point mesh to calculate.
    *   Output: The screened potential is written to `potential.dat`.
        Additional info for the tight-binding code is written to `weights*`
        (atomic site info), `kweights*` (k-point info), and `sym_ops*`
        (symmetry operations info, only present if k-point symmetry is
        calculated).

*   `TB.x`:

    *   Run as: `mpirun -np * ./TB.x`
    *   Input: The program requires all the output files produced by
        `screening.x`, and reads a text file called `TB.in`. An example is
        provided in the `example` directory. All variables must be
        present in the given order; note however that the `sym_ops*` filename
        line is only present if the file is present.  An explanation for each
        variable is given in the example. The input file specifies the
        dimensions of the supercell, physical parameters used in the
        simulation, the names of the additional input files produced by
        `screening.x`, and some information necessary for the parallel linear
        algebra routines.
    *   Output: The total DOS is written to `dos.dat`, and the LDOS of each
        site to `ldos.*.dat`. The electronic density is written to
        `density.dat`. Note that both the LDOS and density are only written for
        symmetrically non-equivalent sites (as specified in the `weights*`
        file).

*   `adatom_post.x`:

    *   Run as: `./adatom_post.x`
    *   Input: The program reads both the `TB.in` input file used by `TB.x` and
        the `ldos.*.dat` and `density.dat` output files it has produced. It
        also reads the `weights*` file from `screening.x`. No further input is
        needed.
    *   Output: The radially-averaged DOS is written to `rdos.*.dat`, and the
        radially-averaged electronic density is written to `rdensity.dat`.

License
-------

Copyright &copy; 2016, Fabiano Corsetti

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1.   Redistributions of source code must retain the above copyright notice,
     this list of conditions and the following disclaimer.

2.   Redistributions in binary form must reproduce the above copyright notice,
     this list of conditions and the following disclaimer in the documentation
     and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

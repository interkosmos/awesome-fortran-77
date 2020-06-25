# Awesome FORTRAN 77 [![Awesome](https://awesome.re/badge.svg)](https://awesome.re)
A curated list of awesome FORTRAN 77 frameworks, libraries, and software.

 * [Compilers](#compilers)
    + [Free Compilers](#free-compilers)
    + [Commercial Compilers](#commercial-compilers)
    + [Translators](#translators)
  * [Software](#software)
    + [Repositories & Directories](#repositories---directories)
    + [Numerical Software](#numerical-software)
    + [Pseudo-Random Number Generators](#pseudo-random-number-generators)
    + [Graphics, Plotting](#graphics--plotting)
    + [Message Passing](#message-passing)
    + [Serialisiation](#serialisiation)
    + [System](#system)
    + [Utilities](#utilities)
    + [Static Code Analysis](#static-code-analysis)
  * [Standards](#standards)
  * [Documentation, FAQs, Tutorials](#documentation--faqs--tutorials)
    + [Introduction](#introduction)
    + [Coding Style](#coding-style)
    + [Transition to Modern Fortran](#transition-to-modern-fortran)
  * [Free Books](#free-books)
  * [Games](#games)
  * [Humor](#humor)

## Compilers
### Free Compilers

  * [BC-FORTRAN 77](https://ftp.sunet.se/mirror/archive/ftp.sunet.se/pub/simtelnet/msdos/fortran/bcf7713b.zip) (DOS): Freeware compiler/linker, supporting the ANSI/ISO FORTRAN 77 language standard only. Written by Andre Koestli.
  * [f2c](http://www.netlib.org/f2c/) (DOS, Win32, Linux, UNIX): An open-source FORTRAN 77-to-C transpiler. The output can be compiled with any K&R or ANSI C compiler.
  * [g77](http://www.kilmnj.com/g77/) (DOS, Win32, Linux, UNIX): The open-source FORTRAN 77 compiler of the GNU Compiler Collection (GCC) prior version 4.0. Has been replaced by GNU Fortran since then. Available for Win32 through MinGW or Cygwin.
  * [lcc-win32](http://mermaja.act.uji.es/docencia/is37/data/LCC/LCC-Win32.html) (Win32): C and FORTRAN development environment, by Jacob Navia.
  * [Open Watcom FORTRAN 77](http://www.openwatcom.com/) (DOS, OS/2, Win32, Linux): Multi-platform open-source compiler, based on the commercial Watcom FORTRAN 77. Includes the DOS/4GW 32-bit DOS extender with run-time and virtual memory support up to 32 MiB.

### Commercial Compilers

  * [Absoft ProFortran F77](https://web.archive.org/web/20000706211705/http://www.absoft.com/) (DOS, Win32, MacOS 7.5, Linux): Full ANSI/ISO FORTRAN 77 compiler with VAX/VMS language extensions; Cray pointers; Sun, IBM, and HP extensions; and LS Fortran language extensions.
  * [Apogee-FORTRAN 77](https://web.archive.org/web/19990221061617/http://www.apogee.com/page11.html) (UNIX): Supports almost all language extensions found in Sun, VAX/VMS, Cray, and MIL-STD 1753 FORTRAN 77.
  * [Fujitsu Fortran](https://web.archive.org/web/20000305112103/http://www.tools.fujitsu.com/linux/index.shtml) (Linux): Fortran 90/95 compiler with FORTRAN 77 support.
  * [Lahey F77L3-EM/32](https://web.archive.org/web/19980523094643/http://www.lahey.com/em32.htm) (DOS): Highly optimised 32-bit FORTRAN 77 compiler for DOS that supports the full ANSI standard and the most popular DEC VAX, IBM VS, and Fortran 90 features. Includes a basic graphics library. The more advanced Graphoria library was sold separately.
  * [Microsoft FORTRAN 5.1](https://winworldpc.com/product/microsoft-fortran/5x) (DOS, OS/2, Win16): Microsoft’s implementation of a FORTRAN 77 compiler. The compiler provides its own `HIMEM.SYS` to access up to 16 MiB of real memory.
  * [Microway NDP Fortran](https://web.archive.org/web/19980120153749/http://www.microway.com/ndpfort.html) (DOS, Win32, OS/2, Linux, UNIX): Globally optimising VMS-compatible compiler that supports Intel, Cyrix, and Weitek co-processors. Full FORTRAN 77 support with DOD, FORTRAN 66, BSD4.2, and VMS extensions. Linkable with NDP C/C++/Pascal. Includes linker, DOS extender, and librarian.
  * [Prospero Pro Fortran-77](https://web.archive.org/web/20040415132626/http://www.prosperosoftware.com/f2iix.html) (DOS, OS/2): FORTRAN 77 development environment that includes compiler, linker, librarian, source level debugger and IDE. The F77PC Graphics BIOS Library even provides “turtle graphics”, better known from the Logo programming language. Prospero also published FORTRAN 77 compilers for Atari ST, OS/2, and Sinclair QL. An optional DOS Extender Kit was available. The PC version was supported until 2003.
  * [Ryan-McFarland RM/FORTRAN](https://vetusware.com/download/RMFORTRAN%202.43/?id=9864) (DOS): Supports the full ANSI FORTRAN 77 standard and includes popular extensions from VAX, VS, and FORTRAN 66. Featured code optimisation and floating-point co-processor support/emulation.
  * [Salford FTN77](https://www.silverfrost.com/53/ftn77_personal_edition.aspx) (DOS, Win32, UNIX): The FTN77 compiler uses the Salford-developed DBOS extender. Later versions of FTN77 were ported to Win32 and UNIX. A personal edition is available free of charge.
  * [Siemens Nixdorf SNI F77](https://web.archive.org/web/19970210142028/http://www.linuxland.de/fortran/fortran.html) (Linux): optimised compiler for Intel 486/586/686, developed by Siemens Nixdorf Informationssysteme.
  * [Sun WorkShop Compiler FORTRAN 77](https://web.archive.org/web/19970803185812/http://www.sun.com/workshop/literature/ProductGuide/SWP.SW/PWSF_tech.html) (UNIX): ANSI FORTRAN 77 compiler, with VMS, Cray, and MIL-STD 1753 extensions; POSIX binding; and interlanguage calling (C, C++, Pascal, Ada).
  * [The Portland Group PGF77 Workstation](https://web.archive.org/web/20000815200338/http://www.pgroup.com/prodworkpgf77.htm) (Win32, Linux, UNIX): SMP-parallel FORTRAN 77 compiler for Intel Pentium Pro processors.
  * [Watcom FORTRAN 77](https://web.archive.org/web/19970719190132/http://www.powersoft.com/products/languages/fort77.html) (DOS, OS/2, Win16, Win32): 32-bit protected-mode compiler that supports up to 4 GiB memory. Run-time compatible with Watcom C. Includes protected-mode version of compiler, source-level debugger, full ANSI FORTRAN 77 language support plus extensions, profiler, IBM SAA compatibility, linker, librarian, and make utility. Requires 386/486 and DOS extender. In 2003, the compiler was made available free of charage as Open Watcom.

### Translators

  * [C2F](https://web.archive.org/web/20100328094906/http://home.earthlink.net/~dave_gemini/): David Frank’s C to FORTRAN 77 translator.
  * [F2J](https://web.archive.org/web/20040220231611/http://www.npac.syr.edu/projects/pcrc/f2j.html): FORTRAN 77-to-Java converter.
  * [Mac F2C](https://web.archive.org/web/20040404090614/http://www.alumni.caltech.edu/~igormt/Mac_F2C.html): Free FORTRAN-to-C translator for Mac OS.

## Software

### Repositories & Directories

  * [Computer Physics Communications Program Library](http://cpc01.cs.qub.ac.uk/): Almost 2000 programs in computational physics and chemistry.
  * [GAMS](https://gams.nist.gov/): Repository of mathematical and statistical software components of use in computational science and engineering.
  * [John Burkardt](https://people.sc.fsu.edu/~jburkardt/f77_src/f77_src.html): Collection of FORTRAN 77 software.
  * [NCAR’s Mathematical and Statistical Libraries](https://web.archive.org/web/20050211013748/http://www.scd.ucar.edu/softlib/mathlib.html): Collection of mathematical and statistical software in FORTRAN 77.
  * [Netlib](http://www.netlib.org/): Collection of mathematical software.

### Numerical Software

  * [ADIFOR](https://www.mcs.anl.gov/research/projects/adifor/): Automatic differentiation of FORTRAN 77 programs.
  * [ALFPACK](https://web.archive.org/web/20050207230333/http://www.scd.ucar.edu/softlib/ALFPACK.html): Legendre functions of first kind.
  * [ARPACK](https://www.caam.rice.edu/software/ARPACK/): Collection of FORTRAN 77 subroutines designed to solve large scale eigenvalue problems.
  * [ASAD](https://web.archive.org/web/20060113132300/www.atm.ch.cam.ac.uk/acmsu/asad/index.html): Software package developed for creating and integrating chemistry schemes in atmospheric models without the need to write any FORTRAN code to solve the chemical rate equations.
  * [Aztec](http://www.cs.sandia.gov/CRF/aztec1.html): An iterative sparse linear solver package.
  * [BLAS](http://www.netlib.org/blas/): Reference implementation of the Basic Linear Algebra Subprograms.
  * [BiM](http://web.math.unifi.it/users/brugnano/BiM/index.html): Implements a variable order-variable stepsize method for (stiff) initial value problems for ODEs.
  * [CERNLIB](https://cernlib.web.cern.ch/cernlib/): CERN general-purpose program library.
  * [CMLIB](https://gams.nist.gov/cgi-bin/serve.cgi/Package/CMLIB/): NIST core math library.
  * [DAEPAK](http://www.netlib.org/contin/manpak/): Differential algebraic equations.
  * [DASPK](http://www.netlib.org/ode/): Differential-algebraic system solver (BDF/Krylov method).
  * [DCDFLIB](https://biostatistics.mdanderson.org/SoftwareDownload/SingleSoftware/Index/21): Cumulative distribution functions, inverses, and parameters for common statistical distributions.
  * [EIGENTRI](http://calgo.acm.org/): Set of Fortran programs for reducing a nonsymmetric matrix to tridiagonal form, computing the eigenvalues of the tridiagonal matrix, improving the accuracy of an eigenvalue, and computing the corresponding eigenvector
  * [EISPACK](http://www.netlib.org/eispack/): Collection of FORTRAN 77 subroutines that compute the eigenvalues and eigenvectors. Superseded by LAPACK.
  * [FFTPACK](https://www2.cisl.ucar.edu/resources/legacy/fft5): FORTRAN 77 library of fast Fourier transforms.
  * [FISHPAK](https://web.archive.org/web/20040303192401/http://www.scd.ucar.edu/css/software/fishpack/): Subprograms for the solution of separable elliptic PDEs.
  * [GSLIB](http://www.statios.com/GSLIB/index.html): Geostatistical Software Library. Collection of geostatistical programs developed at Stanford University.
  * [HSL](http://www.hsl.rl.ac.uk/archive/index.html): Collection of FORTRAN codes for large scale scientific computation.
  * [HSL](https://web.archive.org/web/20020820093441/http://www.cse.clrc.ac.uk/Activity/HSL): Harwell Subroutine Library. Ccollection of ISO FORTRAN 77 subprograms for large scale scientific computation.
  * [IMSL](https://web.archive.org/web/20010405123018/http://www.vni.com/products/imsl/imslfort.html): More than 900 FORTRAN 77 subroutines for use in general applied mathematics, statistical data analysis and presentation in scientific and business applications.
  * [ITPACK](http://www.netlib.org/itpack/): Subroutines for solving large sparse linear systems by iterative methods.
  * [LAIPE](http://www.equation.com/servlet/equation.cmd?fa=laipe): High-performance package of parallel direct solvers.
  * [LANCELOT](http://www.numerical.rl.ac.uk/lancelot/lancelot.html): Standard FORTRAN 77 package for solving large-scale nonlinearly constrained optimisation problems.
  * [LAPACK](http://www.netlib.org/lapack/): Subroutines for solving systems of simultaneous linear equations, least-squares solutions of linear systems of equations, eigenvalue problems, and singular value problems.
  * [LINPACK](http://www.netlib.org/linpack/): Collection of Fortran subroutines that analyse and solve linear equations and linear least-squares problems. Superseded by LAPACK.
  * [MGMPI](https://web.archive.org/web/20050405232502/http://cosmos.ucsd.edu/mgmpi/index-home.html): Parallel multigrid solver for 3D elliptic problems.
  * [MINPACK](http://www.netlib.org/minpack/): Subroutines for solving nonlinear equations and nonlinear least squares problems.
  * [MINUIT](https://github.com/ramos/minuit): Numerical minimisation program, by CERN.
  * [MPFUN77](https://www.netlib.org/mpfun/): Multiple precision arithmetic in FORTRAN 77.
  * [MUDPACK](https://web.archive.org/web/20070607150336/http://www.cisl.ucar.edu/css/software/mudpack/): Multigrid software for elliptic PDEs.
  * [NCARM](https://web.archive.org/web/20050204103302/http://www.scd.ucar.edu/softlib/NCARM.html): Non-portable math library that is very Cray parallel-vector specific, and dates back to the mid-1970s.
  * [ODEPACK](https://people.sc.fsu.edu/~jburkardt/f77_src/odepack/odepack.html): FORTRAN 77 library which implements a variety of solvers for ODEs.
  * [PIM](https://web.archive.org/web/19991001072853/http://www.mat.ufrgs.br/%7Erudnei/pim/pim-i.html): Collection of FORTRAN 77 routines designed to solve systems of linear equations (SLEs) on parallel computers using a variety of iterative methods.
  * [PSIDE](https://web.archive.org/web/19991006225752/https://www.cwi.nl/cwi/projects/PSIDE/): Parallel software for implicit differential equations.
  * [REGRIDPACK](https://web.archive.org/web/20041214084306/http://www.scd.ucar.edu/softlib/REGRIDPACK.html): FORTRAN 77 routines for interpolating values between one-, two-, three-, and four-dimensional arrays defined on uniform or nonuniform orthogonal grids.
  * [SCILIB](http://www.netlib.org/scilib/): Portable FORTRAN emulation of Cray SCILIB.
  * [SLATEC](http://www.netlib.org/slatec/): Comprehensive software library containing over 1400 general purpose mathematical and statistical routines written in FORTRAN 77.
  * [SLEIGN2](http://www.math.niu.edu/SL2/): Sturm-Liouville problems.
  * [SPARSEKIT](https://www-users.cs.umn.edu/~saad/software/SPARSKIT/): Basic tool-kit for sparse matrix computations.
  * [SPBLAS](https://math.nist.gov/spblas/original.html): Sparse matrix computational kernels.
  * [SPHEREPACK](https://web.archive.org/web/20051107222328/http://www.scd.ucar.edu/css/software/spherepack/): Collection of FORTRAN programs that facilitates computer modeling of geophysical processes.
  * [STARPAC](https://water.usgs.gov/software/OTIS/addl/starpac/nls.html): Time-series and regressions package.
  * [TLCPACK](https://web.archive.org/web/20050309013545/http://www.scd.ucar.edu/softlib/TLCPACK.html): Suite of FORTRAN 77 routines for interpolating values between one-, two-, three-, and four-dimensional arrays defined on uniform and nonuniform orthogonal grids.
  * [TOMS](http://www.netlib.org/toms/): Collected algorithms of the ACM.
  * [UMFPACK](https://web.archive.org/web/20041209062517/http://www.cise.ufl.edu/research/sparse/umfpack/): Set of routines for solving unsymmetric sparse linear systems.
  * [WSMP](https://researcher.watson.ibm.com/researcher/view_group.php?id=1426): A collection of algorithms for efficiently solving large sparse systems of linear equations.

### Pseudo-Random Number Generators

  * [ACORN](http://acorn.wikramaratna.org/download.html): Original FORTRAN 77 version of the *Additive Congruential Random Number* (ACORN) generator, by R. S. Wikramaratna.
  * [Mersenne Twister](http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/VERSIONS/FORTRAN/fortran.html): Various adaptions of the Mersenne Twister (MT19937), collected by Makoto Matsumoto.
  * [RANLIB](https://people.sc.fsu.edu/~jburkardt/f77_src/ranlib/ranlib.html): General random number generators.
  * [RANARRAY](https://www-cs-faculty.stanford.edu/~knuth/programs.html): Portable random number generator.
  * [Using random number generators on UNIX systems](https://web.archive.org/web/20060303093233/https://www.cisl.ucar.edu/zine/96/spring/articles/3.random-6.html): FORTRAN 77 adaption of the ”Interger Version 2” PRNG, originally written by Steven K. Park and Keith W. Miller in Pascal.

### Graphics, Plotting

  * [DISLIN](http://www.dislin.de/): High-level scientific plotting library. Free for non-commercial use.
  * [DRAWCGM](https://web.archive.org/web/20050418002941/http://www.psc.edu/general/software/packages/drawcgm/drawcgm.html): Package of FORTRAN and C routines which can be used to create CGM metafiles.
  * [GGG](http://www.unige.ch/~hairer/GGGraphics.html): Geneva Group Graphics library for FORTRAN and LaTeX.
  * [GRAFIC](https://web.archive.org/web/20021204062152/http://www.physics.ohio-state.edu/~koepf/grafic.html): Library of interactive graphics routines for UNIX.
  * [moriplot](https://web.archive.org/web/19970727210155/http://phase.etl.go.jp/contrib/moriplot/): Plotting program for Sun OS 1.4.2, by Masatake Mori.
  * [PGPLOT](https://www.astro.caltech.edu/~tjp/pgplot/): Scientific graphics subroutine library for FORTRAN 77 and C.
  * [PLplot](http://plplot.sourceforge.net/): Cross-platform software package for creating scientific plots.
  * [PSPLOT](https://web.archive.org/web/20050831063710/http://www.nova.edu/ocean/psplot.html): Free Fortran-callable PostScript plotting library.

### Message Passing

  * [F77_ZMQ](https://github.com/zeromq/f77_zmq): ZeroMQ bindings for FORTRAN 77.

### Serialisiation

  * [Harwell-Boeing](https://people.sc.fsu.edu/~jburkardt/data/hb/hb.html): File format used to store sparse matrices.

### System

  * [ANSI.SYS](https://web.archive.org/web/20020604045509/http://www.ccwr.ac.za/~lynch2/ansi.sys.html): Using the `ANSI.SYS` driver to enhance FORTRAN 77 programs (MS-DOS).

### Utilities

  * [Perl for Fortran](https://web.archive.org/web/20011004055401/https://marine.rutgers.edu/po/perl.html): Several Perl scripts useful when working with Fortran programs.

### Static Code Analysis

  * [ftnchek](https://www.dsm.fordham.edu/~ftnchek/): Static analyser for FORTRAN 77 programs that is designed to detect certain errors in a FORTRAN program that a compiler usually does not.
  * [FLOPPY](http://www.netlib.org/floppy/): Software tool that takes as input a file of FORTRAN 77 code and checks it according to various “coding conventions”. FLOW produces reports based on the output of FLOPPY.

## Standards

  * [ANSI X3J3/90.4](https://web.archive.org/web/20070702014846/https://www.fortran.com/fortran/F77_std/f77_std.html): Official FORTRAN 77 standard.
  * [NIST FORTRAN 77 Test Suite](http://www.fortran-2000.com/ArnaudRecipes/fcvs21_f95.html): FCVS78 validation programs by NIST Information Technology Laboratory (ITL) to check compiler conformance to the FORTRAN 77 language standard ([documentation](https://sourceryinstitute.github.io/RefactorF4Acc-test-docs/tests/NIST_F78_test_suite/fcvs21_f95/doc/index.html)).

## Documentation, FAQs, Tutorials
### Introduction

  * [FORTRAN FAQ](http://www.faqs.org/faqs/fortran-faq/): Frequently Asked Questions, compiled by Keith H. Bierman.
  * [Professional Programmer’s Guide to FORTRAN 77](https://www.tat.physik.uni-tuebingen.de/~kley/lehre/ftn77/f77prof.pdf): Introduction by Clive G. Page.
  * [User Notes on Fortran Programming](http://www.ibiblio.org/pub/languages/fortran/) (UNFP): An open cooperative practical guide.

### Coding Style

  * [FORTRAN 77 Coding Guidelines](https://web.archive.org/web/20000530084632/http://ics.uci.edu/pub/levine/F77_Style_Guide): By David L. Levine.
  * [Fortran Coding Style for the Modern Programmer](https://web.archive.org/web/20000421165600/http://studbolt.physast.uga.edu/templon/fortran/fortran_style): By Glen Reesor.

### Transition to Modern Fortran

  * [Fortran 90 for the FORTRAN 77 Programmer](https://www.nsc.liu.se/~boein/f77to90/f77to90.html)
  * [Fortran 95 for the FORTRAN 77 Programmer](https://web.archive.org/web/20080123082959/http://www.soks.org/view/Fortran95ForFortran77Programmers)

## Free Books

  * [Numerical Recipes](http://numerical.recipes/oldverswitcher.html): Free download of the book series on numerical analysis and algorithms in FORTRAN 77 and Fortran 90.

## Games

  * [COMEL](https://core.ac.uk/download/pdf/36712909.pdf): A military communications-oriented war game, developed by the Joint Telecommunications Staff Officers’ Course at Keesler Air Force Base in 1983.
  * [Colossal Cave Adventure](https://ifarchive.org/indexes/if-archive/games/source/): Several FORTRAN 77 ports of the famous text adventure, originally written by Will Crowther.
  * [Conway’s Game of Life](https://github.com/owainkenwayucl/fortlife): Simple implementation in FORTRAN 77 for DOS.
  * [LIFE_SERIAL](https://people.sc.fsu.edu/~jburkardt/f77_src/life_serial/life_serial.html): Another implementation of Conway’s Game of Life, by John Burkardt.
  * [Empire](https://www.classicempire.com/): Strategy and tactics war game, ported to various platforms.
  * [RFK-77](http://cyber.dabamos.de/programming/fortran/rfk/): Implementation of [robotfindskitten](http://robotfindskitten.org/) in FORTRAN 77 for 32-bit DOS.
  * [Rock, Paper, Scissors](https://craftofcoding.files.wordpress.com/2018/06/prac_fortrps.pdf): FORTRAN 77 and Fortran 90/95 implementations.
  * [University of Toronto FORTRAN Games](http://freshmeat.sourceforge.net/projects/fortran-games): Chess, Tic Tac Toe, Minefield, and other games, some dating back to the late 70’s.

## Humor

  * [COME FROM](http://www.fortranlib.com/gotoless.htm): A linguistic contribution of GOTO-less programming.
  * [Fortran Songs](http://www.poppyfields.net/filks/f.html). Singin’ and dancin’ to the joy of programming in FORTRAN.
  * [Real Programmers Don’t Use PASCAL](https://web.mit.edu/humor/Computers/real.programmers): Essay about computer programming, written by Ed Post.


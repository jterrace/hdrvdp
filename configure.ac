AC_INIT(vdp, 1.6)
AM_INIT_AUTOMAKE

AC_MSG_CHECKING([for pfs library])
PKG_CHECK_MODULES(PFS, pfs >= 1.2,,
  AC_MSG_ERROR(pfstools package required. Download it from http://pfstools.sourceforge.net/. If it is installed, make sure that PKG_CONFIG_PATH points to the directory, in which pfs.pc is located.))
AC_SUBST(PFS_CFLAGS)
AC_SUBST(PFS_LIBS)

dnl
dnl Internal Debug
dnl
CXXFLAGS="-O3"
AC_ARG_ENABLE(debug,
		AC_HELP_STRING([--enable-debug],
				 [Enable debug symbols and asserts]),
		 [if test "$enable_debug" = "yes" ; then
		     CXXFLAGS="-g"
		     AC_DEFINE([DEBUG])
		  fi],[
                  enable_debug="no"
                  AC_DEFINE(NDEBUG, 1, [Disable assertion checks])
                  ])
        
AC_PROG_CXX
AC_PROG_CC

PKG_CHECK_MODULES(FFTW, fftw3f >= 3.0,,                             
AC_CHECK_LIB(fftw3f, fftwf_plan_dft_c2r_2d,,
             [ AC_MSG_ERROR(Can not find fftw3f library. Download it from http://www.fftw.org/. Remember to add --enable-float switch! Make also sure that PKG_CONFIG_PATH points to the directory, in which fftw3f.pc is located.) ],-lm)
)

AC_CHECK_LIB(fftw3f_threads, fftwf_init_threads,
             [ FFTW_THREADS_LIBS="-lfftw3f_threads -lpthread"; AC_DEFINE(HAVE_LIBFFTW_THREADS, 1, [Support for FFTW with multiple threads]) ],
             [ AC_MSG_WARN(Can not find fftw3f_threads library. Thread support will be disabled.) ],$FFTW_LIBS -lpthread)
             
AC_SUBST(FFTW_CFLAGS)
AC_SUBST(FFTW_LIBS)
AC_SUBST(FFTW_THREADS_LIBS)

        
AC_CONFIG_FILES([Makefile src/Makefile distort/Makefile])

AC_MSG_RESULT([
vdp is now configured

  Source directory:     $srcdir
  Installation prefix:  $prefix
  C++ compiler:         $CXX $XTRA_CXXFLAGS $WARN_CXXFLAGS $CXXFLAGS
  LIBS:                 $LIBS

  pfs cflags:           $PFS_CFLAGS
  pfs libs:             $PFS_LIBS

  fftw cflags:          $FFTW_CFLAGS
  fftw libs:            $FFTW_LIBS
  fftw_threads libs:    $FFTW_THREADS_LIBS
  
  Debug mode            $enable_debug
])
                                           
AC_OUTPUT
        


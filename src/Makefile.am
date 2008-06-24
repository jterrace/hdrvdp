bin_PROGRAMS = vdpcmp vdpvis
bin_SCRIPTS  = vdp
man_MANS = vdp.1 vdpcmp.1 vdpvis.1

EXTRA_DIST = $(man_MANS) $(bin_SCRIPTS)

vdpcmp_SOURCES = vdpcmp.cpp \
	csf_psi_nonlinearity.cpp csf_psi_nonlinearity.h \
	CSF.cpp CSF.h \
	array2d_algorithm.cpp array2d_algorithm.h \
	otf.cpp otf.h \
	nonlinearity.cpp nonlinearity.h csf_psi_mapping.h \
	csf_filter.cpp csf_filter.h \
	cortex_transform.cpp cortex_transform.h \
	dump_image_aspect.cpp dump_image_aspect.h \
	fftutils.cpp fftutils.h fftw_array2d.h

vdpvis_SOURCES = vdpvis.cpp

LIBS += $(PFS_LIBS) $(FFTW_THREADS_LIBS) $(FFTW_LIBS) 
INCLUDES = $(PFS_CFLAGS) $(FFTW_CFLAGS)

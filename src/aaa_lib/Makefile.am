
## ###########################################################################
## File:    Makefile.am
##
## Purpose: Automatically generate Makefile.in for autoconf using automake.
##
## Notes:   The final makefile produced by autoconf will use the following:
##
##     COMPILE = $(CC) $(DEFS) $(INCLUDES) $(CPPFLAGS) $(CFLAGS)
##     LINK = $(CC) $(CFLAGS) $(LDFLAGS) -o $@
##
## Author:  Michael Holst
## ###########################################################################

AUTOMAKE_OPTIONS = no-dependencies
SUFFIXES = .y .l .H .h .m .C .cc .c .f .java .o
always_built_SUBDIRS =
SUBDIRS = $(always_built_SUBDIRS) .

BASE_LIBS   = ../base/libbase.la

BIOM_LIBS   = ../biom/libbiom.la

TRIANGLE_LIBS = ../triangle/libtriangle.la
#if BUILD_TRIANGLE
LIBTRIANGLE   = libtriangle.la
#else
#LIBTRIANGLE   = 
#endif

TETGEN_LIBS = ../tetgen/libtetgen.la
#if BUILD_TETGEN
LIBTETGEN   = libtetgen.la
#else
#LIBTETGEN   = 
#endif

libdir = ${prefix}/lib
lib_LTLIBRARIES = libgamer.la ${LIBTRIANGLE} ${LIBTETGEN}

libgamer_la_SOURCES = 
libgamer_la_LIBADD  = ${BASE_LIBS} ${BIOM_LIBS} ${TETGEN_LIBS} ${TRIANGLE_LIBS}
libgamer_la_LDFLAGS = -version-info ${FETK_VERSION}

libtriangle_la_SOURCES = 
libtriangle_la_LIBADD  = ${TRIANGLE_LIBS}
libtriangle_la_LDFLAGS = -version-info ${FETK_VERSION}

libtetgen_la_SOURCES = 
libtetgen_la_LIBADD  = ${TETGEN_LIBS}
libtetgen_la_LDFLAGS = -version-info ${FETK_VERSION}

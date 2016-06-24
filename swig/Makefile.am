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
ACLOCAL_AMFLAGS  = -I m4
SUFFIXES = .y .l .H .h .m .C .cc .c .f .java .o
always_built_SUBDIRS =
SUBDIRS = $(always_built_SUBDIRS) .

INCLUDES               = @gamer_inc@ @maloc_inc@ 
AM_CFLAGS              = @profile@ @wall_pedantic_ansi@ 

BUILT_SOURCES          = cgamer_wrap.cc
SWIG_SOURCES           = gamer.i gamer_typemaps.i biom.i post.i pre.i

pkgpython_PYTHON       = cgamer.py  ${top_srcdir}/src/__init__.py  ${top_srcdir}/src/upy_gui.py
pkgpyexec_LTLIBRARIES  = _cgamer.la
_cgamer_la_SOURCES     = cgamer_wrap.cc ${SWIG_SOURCES}
_cgamer_la_CPPFLAGS    = ${SWIG_PYTHON_CPPFLAGS} -I. 
_cgamer_la_LDFLAGS     = -module
_cgamer_la_LIBADD      = @gamer_lib@

cgamer_wrap.cc : ${top_srcdir}/src/gamer.i
	${SWIG} ${SWIG_PYTHON_OPT} -I. ${maloc_inc} ${gamer_inc} -o cgamer_wrap.cc -oh cgamer_wrap.h ${top_srcdir}/src/gamer.i

clean-local :
	rm -f cgamer.py cgamer_wrap.cc


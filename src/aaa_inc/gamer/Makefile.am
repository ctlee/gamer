
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

TRIANGLE_HEADS = triangle.h
TETGEN_HEADS   = tetgen.h
BIOM_HEADS     = biom.h

GAMER_HEADS    = gamer.h gamer_base.h
OTHER_HEADS    = ${TRIANGLE_HEADS} ${TETGEN_HEADS} ${BIOM_HEADS}

includedir = ${prefix}/include/gamer
include_HEADERS = ${GAMER_HEADS} ${OTHER_HEADS}

distclean-local: clean-local

clean-local:
	-rm -f ${OTHER_HEADS}


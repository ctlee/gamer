
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
SUBDIRS = $(always_built_SUBDIRS) . @vtst@

bin_PROGRAMS     = ImproveSurfMesh
ImproveSurfMesh_SOURCES    = ImproveSurfMesh.C
ImproveSurfMesh_LDADD      = @gamer_lib@ @maloc_lib@ @mpi_lib@ @rl_lib@ @liberty_lib@ @socket_lib@ @nsl_lib@ @thread_lib@ @math_lib@ @zfence_lib@ @profile_lib@
INCLUDES         = @gamer_inc@ @maloc_inc@
AM_CFLAGS        = @profile@ @wall_pedantic_ansi@


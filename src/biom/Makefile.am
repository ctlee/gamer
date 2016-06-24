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
THISLIB = biom

noinst_LTLIBRARIES = libbiom.la
libbiom_la_SOURCES = Triangulate.C GemMesh.C ReadLattice.C CorrectNormals.C ReadWriteMesh.C SurfSmooth.C MarchingCube.C PDB2Volume.C GenerateSphereMesh.C SurfExtract.C SurfRefine.C ReadPDB.C PDB2Mesh.C RegionGrowing.C SurfaceMesh.C
INCLUDES           = @maloc_inc@
AM_CFLAGS          = @profile@ @wall_pedantic_ansi@

all-local:
	-cp -p ${top_srcdir}/src/${THISLIB}/gamer/*.h ../aaa_inc/gamer/.


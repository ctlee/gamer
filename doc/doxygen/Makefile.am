
## ###########################################################################
## File:    Makefile
##
## Purpose: Automatically generate documentation using doxygen.
##
## Notes:   Both HTML and LaTeX (dvi/ps/pdf) documentation is created.
##
## Author:  Michael Holst
## ###########################################################################

docdir = @prefix@/doc/gamer

doc_DATA = gamer_doc

gamer_doc: 
	$(doxygen_path) $(top_srcdir)/doc/doxygen/gamer.dox

install-docDATA:
	mkdir -p $(docdir)
	cp -r ../api $(docdir)/

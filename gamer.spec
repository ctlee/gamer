## ###########################################################################
## File:    gamer.spec
##
## Purpose: Spec file for building RPMS
##
## Notes:   If this is installed in the top directory the user can build a
##          full set of src and arch rpms with one command:
##
##          rpm -ta gamer.tar.gz
##
## Author:  Michael Holst
## ###########################################################################

Summary: Geometric Mesh generator
Name: gamer
Version: 0.1
Release: 1
Copyright: GPL
Group: Applications/Science
Prefix: /usr/local
Buildroot: %{_topdir}/buildroot
Source: gamer.tar.gz
URL: http://cam.ucsd.edu/~mholst
Packager: Michael Holst <mholst@math.ucsd.edu>
%description
GAMER (Geometry-preserving Adaptive Mesher) is a small self-contained simplex 
mesh-generation package for generating high-quality tetrahedral meshes with 
complex boundary and imbedded surface triangulations.  It uses various 
surface triangulation and improvement algorithms, together with Tetgen for 
volume mesh generation.  It is designed to work with collaboratively with 
the FETK family of research software (MALOC, SG, MC, APBS).  GAMER can be 
used as a stand-alone mesh generation package that takes input from files 
and produces output meshes in files in one of several formats.  Since it 
uses the MALOC platform portability layer from FETK, it can also be started 
as a mesh generation server, taking input from UNIX and/or INET sockets, 
and producing ouput to UNIX and/or INET sockets.

%prep
%setup -n gamer

%build

./configure --prefix=${RPM_BUILD_ROOT}/%{prefix}
make 

%install
mkdir -p ${RPM_BUILD_ROOT}/%{prefix}
make install

%clean
rm -rf ${RPM_BUILD_ROOT}

%post

%postun

%files
%defattr(-,root,root)
%{prefix}/lib
%{prefix}/include
%doc AUTHORS COPYING INSTALL NEWS ChangeLog doc

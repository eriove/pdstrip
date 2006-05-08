# You can get a compiler binary from http://gcc.gnu.org/wiki/GFortranBinaries
# or http://ftp.g95.org/ or ftp://download.intel.com/software/products/compilers/downloads/

FC = gfortran

#FC = g95

#FC = ifort
FFLAGS = -g

all: pdstrip

pdstrip: pdstrip.f90
	nice $(FC) -o $@ $(FFLAGS) $^

#%.o : %.f90
#	nice $(FC) -c $(FFLAGS) $^

clean:
	$(RM) *.mod

dist-clean: clean
	$(RM) *~ svn-commit.tmp pdstrip
	cd examples; $(MAKE) $@
	cd windows; $(MAKE) $@

# NOTE do a svn update first, to kick svn info into the right revision.
# Archives are created in the parent directory
dist:
	REVISION=`svn info|egrep "^Revision"|sed "s/Revision: //"`;\
	$(MAKE) "pdstrip-bin_rev"$$REVISION".tar.gz";\
	$(MAKE) "pdstrip-src_rev"$$REVISION".tar.gz"

pdstrip-src_rev%.tar.gz: dist-clean
	cd ..; tar -cvzf $@ --exclude *.svn* pdstrip

pdstrip-bin_rev%.tar.gz: pdstrip
	cd ..; tar -cvzf $@ --exclude *.svn*  --exclude *Makefile pdstrip/doc \
		pdstrip/examples pdstrip/README pdstrip/LICENSE pdstrip/pdstrip

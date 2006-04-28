# You can get a compiler binary from http://gcc.gnu.org/wiki/GFortranBinaries
# or http://ftp.g95.org/ or ftp://download.intel.com/software/products/compilers/downloads/

#FC = /home/bastiaan/programming/fortran/gfortran/irun/bin/gfortran

#FC = g95

FC = ifort
FFLAGS = -g

all: pdstrip

pdstrip: pdstrip.f90
	nice $(FC) -o $@ $(FFLAGS) $^

#%.o : %.f90
#	nice $(FC) -c $(FFLAGS) $^

clean:
	rm -rf *.mod

dist-clean: clean
	rm -rf *~ svn-commit.tmp pdstrip
	cd examples; $(MAKE) $@
	cd windows; $(MAKE) $@

pdstrip-%.tar.gz: dist-clean
	cd ..; tar -cvzf $@ --exclude *.svn* pdstrip
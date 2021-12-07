# ================================================
# flatPERM makefile

# Compiler and Linker
# fc          := gfortran
# Use HDF5 compile script
# run .../h5fc -show in terminal to see compile options used
# fc			:= /usr/local/hdf5/bin/h5fc
fc	:= h5fc

OMP_FLAG =
# OMP_FLAG += -openmp	# Intel
OMP_FLAG += -fopenmp	# gfortran
# OMP_FLAG += -homp		# Cray

# The target Binary Program
target      := main

# The Directories, Source, Includes, objects, Binary and Resources
srcdir    := src
incdir    := inc
blddir    := build
tgtdir		:= .
datadir		:= savedata
srcext    := f90
modext    := mod
objext    := o

# flags, Libraries and Includes
ldflags		:= -J$(blddir) -I$(incdir)
fcflags   := -O2 -std=f2008 -cpp -fbackslash
fcflags		+= -Wall -Wextra -Warray-bounds #-Wrealloc-lhs-all -Wrealloc-lhs
fcflags		+= -g -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow -fdump-core
# fcflags		+= -fprofile-arcs -ftest-coverage
fcflags		+= -fall-intrinsics -fcheck=bounds,array-temps, -Wno-unused-function
# hdf5libs		+= -I$(hdf5)/include $(hdf5)/lib/libhdf5_fortran.a $(hdf5)/lib/libhdf5.a

# ================================================
sources     := $(shell find $(srcdir) -type f -name *.$(srcext))
objects		:= $(patsubst $(srcdir)/%,$(blddir)/%,$(sources:.$(srcext)=.$(objext)))
modules		:= $(objects:.$(objext)=.$(modext))
mtobj     := $(shell find $(incdir) -type f -name *.$(objext))

# Default make
all: directories $(target)

# Dependencies -- need to list manually because Fortran
$(target): $(objects) $(mtobj)

%.$(objext): %.$(modext)		# prevent some bullshit

$(blddir)/globals.o: $(blddir)/fhash_modules.o

$(blddir)/hash_wrappers.o: $(blddir)/globals.o $(blddir)/lattice.o $(blddir)/fhash_modules.o

$(blddir)/save_data.o: $(blddir)/globals.o $(blddir)/lattice.o

$(blddir)/frontend.o: $(blddir)/globals.o

$(blddir)/buildchain.o: $(blddir)/globals.o $(blddir)/lattice.o $(blddir)/hash_wrappers.o

$(blddir)/microcanonical.o: $(blddir)/globals.o $(blddir)/lattice.o $(blddir)/hash_wrappers.o

# $(blddir)/radius.o: $(blddir)/lattice.o

$(blddir)/flatperm.o: $(blddir)/globals.o $(blddir)/lattice.o $(blddir)/microcanonical.o $(blddir)/radius.o $(blddir)/buildchain.o $(blddir)/hash_wrappers.o

$(blddir)/$(target).o: $(blddir)/globals.o $(blddir)/frontend.o $(blddir)/flatperm.o $(blddir)/save_data.o $(blddir)/hash_wrappers.o

# ================================================
# Link
%: $(blddir)/%.$(objext)
	$(fc) $(fcflags) $(OMP_FLAG) -o $@ $^
#$(ldflags)

# Compile
$(blddir)/%.$(objext): $(srcdir)/%.$(srcext)
	$(fc) $(fcflags) $(ldflags) $(OMP_FLAG) -c $<  -o $@

# ================================================
# Utility targets
.PHONY: clean cleaner directories

# Remake
remake: cleaner all

# Setup directories
directories:
	@mkdir -p $(blddir)
	@mkdir -p $(datadir)
	@cp -vu ~/punim0106/mt_stream_f90-1.11/mt_stream.o $(incdir)		# use prebuilt MT19937 parallel RNG
	@cp -vu ~/punim0106/mt_stream_f90-1.11/mt_stream.mod $(incdir)
	@cp -vu ~/punim0106/mt_stream_f90-1.11/f_jump_ahead_coeff/f_get_coeff.o $(incdir)
	@cp -vu ~/punim0106/mt_stream_f90-1.11/f_jump_ahead_coeff/gf2xe.o $(incdir)

list:
	@echo "target:	"$(target)
	@echo "sources:	"$(sources)
	@echo "objects:	"$(objects)
	@echo "modules:	"$(modules)

# Clean objects
clean:
	rm -f $(blddir)/*.o $(blddir)/*.mod $(blddir)/*.MOD

# Clean objects and binaries
cleaner: clean
	rm -f *~ $(target)

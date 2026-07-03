#------------------------------------------------------------------------------
# Makefile for ExoRT
#

# Set NETCDF Library Directory
# make sure to add the library to your LD_LIBRARY_PATH (Linux)
# or DYLD_LIBRARY_PATH (macOS) to point to your NetCDF lib directory
LIB_NETCDF_STRING := $(shell nf-config --flibs)

# Use the NetCDF library's own netcdf.inc (matches the linked library version)
# rather than a stale local copy. nf-config --includedir gives the directory
# containing the installed netcdf.inc.
INC_NETCDF_STRING := -I$(shell nf-config --includedir)

# directories for up to date builds
SRC_DIR  = ../source
MAIN_DIR = ../source/src.main
MISC_DIR = ../source/src.misc
N68EQUIV_DIR = ../source/src.n68equiv
N84EQUIV_DIR = ../source/src.n84equiv
EXORT_DIR = ../source/src.exort

# build directories
MAIN_BLD = bld.main
N68EQUIV_BLD = bld.n68equiv
N84EQUIV_BLD = bld.n84equiv
EXORT_BLD = bld.exort
EXORT_LIB_BLD = bld.libexort

# shared-library suffix (Stage D: libexort)
ifeq ($(shell uname -s),Darwin)
  SHLIB_EXT = .dylib
else
  SHLIB_EXT = .so
endif

OBJDIR = obj

# Compiler selection: default is ifort (Discover/Linux).
# On Apple Silicon Macs, ifort is unavailable (x86_64 only, no arm64 port).
# Use: make USER_FC=gfortran <target>
#
# gfortran compatibility flags for CESM-heritage code:
#   -fallow-invalid-boz      : BOZ literals in infnan.F90 (ifort accepts silently)
#   -fno-range-check         : integer overflow in mcica_random_numbers.F90
#   -fallow-argument-mismatch: assumed-size arr(*) in wrap_nf.F90
#   -ffree-line-length-none  : long lines in calc_opd_mod.F90
#   -w                       : suppress remaining legacy warnings from CESM infrastructure
#                              (remove when refactoring src.misc)

# set ifort as default
ifeq ($(USER_FC),$(null))
  FC := ifort
else
  FC := $(USER_FC)
endif

# set gfortran flags
ifeq ($(FC),gfortran)
  FCFLAGS += -fallow-invalid-boz \
             -fno-range-check \
             -fallow-argument-mismatch \
             -ffree-line-length-none \
             -w
endif

# OpenMP (Stage E2): the column loops in main.F90 and exort_run_columns are
# !$omp parallel do. Threads only spawn for multi-column runs (if(ncol>1)
# clause), so single-column behavior is unchanged. The large solver work
# arrays were converted from automatics to allocatable (heap) in E2, so
# default thread stacks suffice; if a future large automatic reappears, set
# OMP_STACKSIZE (e.g. 512M). Build without OpenMP with: make OMPFLAGS= <target>
ifeq ($(FC),gfortran)
  OMPFLAGS = -fopenmp
endif
ifeq ($(FC),ifort)
  OMPFLAGS = -qopenmp
endif
FCFLAGS += $(OMPFLAGS)


##  
##
## ==============================================================================
## 2021, 1D-3D common code builds
## Up to date builds
## ==============================================================================

##
## exort  (v2 single bundle: n84 grid + NH3/CO; merges n68equiv + n84equiv)
##
OBJS_EXORT := $(MISC_DIR)/shr_kind_mod.F90 $(SRC_DIR)/exoplanet_mod.F90 $(MAIN_DIR)/sys_rootdir.F90 $(MISC_DIR)/shr_const_mod.F90 $(MISC_DIR)/spmd_utils.F90 $(MISC_DIR)/time_manager.F90 $(MISC_DIR)/wrap_nf.F90 $(MISC_DIR)/infnan.F90 $(MISC_DIR)/ppgrid.F90 $(MISC_DIR)/pmgrid.F90 $(MISC_DIR)/ioFileMod.F90 $(EXORT_DIR)/radgrid.F90 $(MISC_DIR)/physconst.F90 $(MAIN_DIR)/exort_column_mod.F90 $(MAIN_DIR)/io_1D.F90 $(EXORT_DIR)/cloud.F90 $(EXORT_DIR)/kabs.F90 $(MAIN_DIR)/mcica_random_numbers.F90 $(MAIN_DIR)/mcica.F90 $(EXORT_DIR)/rad_interp_mod.F90 $(MAIN_DIR)/rayleigh_data.F90 $(EXORT_DIR)/calc_opd_mod.F90 $(EXORT_DIR)/initialize_rad_mod_1D.F90 $(MAIN_DIR)/planck_mod.F90 $(EXORT_DIR)/model_specific.F90 $(MAIN_DIR)/exo_init_ref.F90 $(MAIN_DIR)/exo_radiation_mod.F90 $(MAIN_DIR)/exort_column_run.F90 $(MAIN_DIR)/main.F90

exort.exe : $(OBJS_EXORT)
	$(FC) -o $@  $(OBJS_EXORT) $(FCFLAGS) $(LIB_NETCDF_STRING) $(INC_NETCDF_STRING) -I$(MISC_DIR)
	$(RM) -r -f $(EXORT_BLD)
	mkdir $(EXORT_BLD)
	mv -f *.mod $(EXORT_BLD)
	mv -f exort.exe ../run


##
## libexort  (Stage D: shared library with the C API from exort_lib_mod.F90;
##            same physics objects as exort.exe, main.F90 replaced by the
##            column/lib modules. Output goes to ../run alongside the .exe)
##
OBJS_EXORT_LIB := $(MISC_DIR)/shr_kind_mod.F90 $(SRC_DIR)/exoplanet_mod.F90 $(MAIN_DIR)/sys_rootdir.F90 $(MISC_DIR)/shr_const_mod.F90 $(MISC_DIR)/spmd_utils.F90 $(MISC_DIR)/time_manager.F90 $(MISC_DIR)/wrap_nf.F90 $(MISC_DIR)/infnan.F90 $(MISC_DIR)/ppgrid.F90 $(MISC_DIR)/pmgrid.F90 $(MISC_DIR)/ioFileMod.F90 $(EXORT_DIR)/radgrid.F90 $(MISC_DIR)/physconst.F90 $(MAIN_DIR)/exort_column_mod.F90 $(MAIN_DIR)/io_1D.F90 $(EXORT_DIR)/cloud.F90 $(EXORT_DIR)/kabs.F90 $(MAIN_DIR)/mcica_random_numbers.F90 $(MAIN_DIR)/mcica.F90 $(EXORT_DIR)/rad_interp_mod.F90 $(MAIN_DIR)/rayleigh_data.F90 $(EXORT_DIR)/calc_opd_mod.F90 $(EXORT_DIR)/initialize_rad_mod_1D.F90 $(MAIN_DIR)/planck_mod.F90 $(EXORT_DIR)/model_specific.F90 $(MAIN_DIR)/exo_init_ref.F90 $(MAIN_DIR)/exo_radiation_mod.F90 $(MAIN_DIR)/exort_column_run.F90 $(MAIN_DIR)/exort_lib_mod.F90

libexort$(SHLIB_EXT) : $(OBJS_EXORT_LIB)
	$(FC) -shared -fPIC -o $@  $(OBJS_EXORT_LIB) $(FCFLAGS) $(LIB_NETCDF_STRING) $(INC_NETCDF_STRING) -I$(MISC_DIR)
	$(RM) -r -f $(EXORT_LIB_BLD)
	mkdir $(EXORT_LIB_BLD)
	mv -f *.mod $(EXORT_LIB_BLD)
	mv -f libexort$(SHLIB_EXT) ../run


##
## n68equiv
##
OBJS_N68EQUIV := $(MISC_DIR)/shr_kind_mod.F90 $(SRC_DIR)/exoplanet_mod.F90 $(MAIN_DIR)/sys_rootdir.F90 $(MISC_DIR)/shr_const_mod.F90 $(MISC_DIR)/spmd_utils.F90 $(MISC_DIR)/time_manager.F90 $(MISC_DIR)/wrap_nf.F90 $(MISC_DIR)/infnan.F90 $(MISC_DIR)/ppgrid.F90 $(MISC_DIR)/pmgrid.F90 $(MISC_DIR)/ioFileMod.F90 $(N68EQUIV_DIR)/radgrid.F90 $(MISC_DIR)/physconst.F90 $(MAIN_DIR)/exort_column_mod.F90 $(MAIN_DIR)/io_1D.F90 $(N68EQUIV_DIR)/cloud.F90 $(N68EQUIV_DIR)/kabs.F90 $(MAIN_DIR)/mcica_random_numbers.F90 $(MAIN_DIR)/mcica.F90 $(N68EQUIV_DIR)/rad_interp_mod.F90 $(MAIN_DIR)/rayleigh_data.F90 $(N68EQUIV_DIR)/calc_opd_mod.F90 $(N68EQUIV_DIR)/initialize_rad_mod_1D.F90 $(MAIN_DIR)/planck_mod.F90 $(N68EQUIV_DIR)/model_specific.F90 $(MAIN_DIR)/exo_init_ref.F90 $(MAIN_DIR)/exo_radiation_mod.F90 $(MAIN_DIR)/exort_column_run.F90 $(MAIN_DIR)/main.F90

n68equiv.exe : $(OBJS_N68EQUIV)
	$(FC) -o $@  $(OBJS_N68EQUIV) $(FCFLAGS) $(LIB_NETCDF_STRING) $(INC_NETCDF_STRING) -I$(MISC_DIR)
	$(RM) -r -f $(N68EQUIV_BLD)
	mkdir $(N68EQUIV_BLD)
	mv -f *.mod $(N68EQUIV_BLD)
	mv -f n68equiv.exe ../run


##
## n84equiv
##
OBJS_N84EQUIV := $(MISC_DIR)/shr_kind_mod.F90 $(SRC_DIR)/exoplanet_mod.F90 $(MAIN_DIR)/sys_rootdir.F90 $(MISC_DIR)/shr_const_mod.F90 $(MISC_DIR)/spmd_utils.F90 $(MISC_DIR)/time_manager.F90 $(MISC_DIR)/wrap_nf.F90 $(MISC_DIR)/infnan.F90 $(MISC_DIR)/ppgrid.F90 $(MISC_DIR)/pmgrid.F90 $(MISC_DIR)/ioFileMod.F90 $(N84EQUIV_DIR)/radgrid.F90 $(MISC_DIR)/physconst.F90 $(MAIN_DIR)/exort_column_mod.F90 $(MAIN_DIR)/io_1D.F90 $(N84EQUIV_DIR)/cloud.F90 $(N84EQUIV_DIR)/kabs.F90 $(MAIN_DIR)/mcica_random_numbers.F90 $(MAIN_DIR)/mcica.F90 $(N84EQUIV_DIR)/rad_interp_mod.F90 $(MAIN_DIR)/rayleigh_data.F90 $(N84EQUIV_DIR)/calc_opd_mod.F90 $(N84EQUIV_DIR)/initialize_rad_mod_1D.F90 $(MAIN_DIR)/planck_mod.F90 $(N84EQUIV_DIR)/model_specific.F90 $(MAIN_DIR)/exo_init_ref.F90 $(MAIN_DIR)/exo_radiation_mod.F90 $(MAIN_DIR)/exort_column_run.F90 $(MAIN_DIR)/main.F90

n84equiv.exe : $(OBJS_N84EQUIV)
	$(FC) -o $@  $(OBJS_N84EQUIV) $(FCFLAGS) $(LIB_NETCDF_STRING) $(INC_NETCDF_STRING) -I$(MISC_DIR)
	$(RM) -r -f $(N84EQUIV_BLD)
	mkdir $(N84EQUIV_BLD)
	mv -f *.mod $(N84EQUIV_BLD)
	mv -f n84equiv.exe ../run





## ==============================================================================

#all
exort            : exort.exe
libexort         : libexort$(SHLIB_EXT)
n68equiv         : n68equiv.exe
n84equiv         : n84equiv.exe

clean:
	$(RM) -f *.o *.mod *.exe
	$(RM) -f ../run/*.exe
	$(RM) -f ../run/libexort$(SHLIB_EXT)
	$(RM) -r -f $(EXORT_BLD)
	$(RM) -r -f $(EXORT_LIB_BLD)
	$(RM) -r -f $(N68EQUIV_BLD)
	$(RM) -r -f $(N84EQUIV_BLD)

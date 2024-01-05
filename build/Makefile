#------------------------------------------------------------------------------
# Makefile for ExoRT
#

# Set NETCDF Library Directory
# make sure to add the library to your LD_LIBRARY_PATH
#LIB_NETCDF := /gpfsm/dulocal/sles11/other/netcdf/3.6.2_intel-14.0.3.174/lib
LIB_NETCDF := /gpfsm/dulocal/sles11/other/netcdf/4.1.2_intel-14.0.3/lib

# use "-lnetcdf" for NETCDF3
# use "-lnetcdf -lnetcdff" for NETCDF4
ARGS := -lnetcdf -lnetcdff


# directories for up to date builds
SRC_DIR  = ../source
MAIN_DIR = ../source/src.main
MISC_DIR = ../source/src.misc
N28ARCHEAN_DIR = ../source/src.n28archean
N42H2O_DIR = ../source/src.n42h2o
N68H2O_DIR = ../source/src.n68h2o
N68EQUIV_DIR = ../source/src.n68equiv
N84EQUIV_DIR = ../source/src.n84equiv
N68EQUIV_EXP_DIR = ../source/experimental/src.n68equiv_exp

# build directories
MAIN_BLD = bld.main
N28ARCHEAN_BLD = bld.n28archean
N42H2O_BLD = bld.n42h2o
N68H2O_BLD = bld.n68h2o
N68EQUIV_BLD = bld.n68equiv
N84EQUIV_BLD = bld.n84equiv
N68EQUIV_EXP_BLD = bld.n68equiv_exp


OBJDIR = obj

ifeq ($(USER_FC),$(null))
  FC := ifort
else
  FC := $(USER_FC)
endif

##  
##
## ==============================================================================
## 2021, 1D-3D common code builds
## Up to date builds
## ==============================================================================

##
## n68equiv
## 
OBJS_N68EQUIV := $(MISC_DIR)/shr_kind_mod.F90 $(SRC_DIR)/exoplanet_mod.F90 $(MAIN_DIR)/sys_rootdir.F90 $(MISC_DIR)/shr_const_mod.F90 $(MISC_DIR)/shr_sys_mod.F90 $(MISC_DIR)/spmd_utils.F90 $(MISC_DIR)/time_manager.F90 $(MISC_DIR)/wrap_nf.F90 $(MISC_DIR)/infnan.F90 $(MISC_DIR)/ppgrid.F90 $(MISC_DIR)/pmgrid.F90 $(MISC_DIR)/ioFileMod.F90 $(N68EQUIV_DIR)/radgrid.F90 $(MAIN_DIR)/input.F90 $(MISC_DIR)/physconst.F90 $(MAIN_DIR)/output.F90 $(N68EQUIV_DIR)/cloud.F90 $(N68EQUIV_DIR)/kabs.F90 $(MAIN_DIR)/mcica_random_numbers.F90 $(MAIN_DIR)/mcica.F90 $(N68EQUIV_DIR)/rad_interp_mod.F90 $(MAIN_DIR)/rayleigh_data.F90 $(N68EQUIV_DIR)/calc_opd_mod.F90 $(N68EQUIV_DIR)/initialize_rad_mod_1D.F90 $(MAIN_DIR)/planck_mod.F90 $(N68EQUIV_DIR)/model_specific.F90 $(MAIN_DIR)/exo_init_ref.F90 $(MAIN_DIR)/exo_radiation_mod.F90 $(MAIN_DIR)/main.F90

n68equiv.exe : $(OBJS_N68EQUIV)
	$(FC) -o $@  $(OBJS_N68EQUIV) -L$(LIB_NETCDF) $(ARGS) -I$(MISC_DIR) 
	$(RM) -r -f $(N68EQUIV_BLD)
	mkdir $(N68EQUIV_BLD)
	mv -f *.mod $(N68EQUIV_BLD)
	mv -f n68equiv.exe ../run


##
## n84equiv
##
OBJS_N84EQUIV := $(MISC_DIR)/shr_kind_mod.F90 $(SRC_DIR)/exoplanet_mod.F90 $(MAIN_DIR)/sys_rootdir.F90 $(MISC_DIR)/shr_const_mod.F90 $(MISC_DIR)/shr_sys_mod.F90 $(MISC_DIR)/spmd_utils.F90 $(MISC_DIR)/time_manager.F90 $(MISC_DIR)/wrap_nf.F90 $(MISC_DIR)/infnan.F90 $(MISC_DIR)/ppgrid.F90 $(MISC_DIR)/pmgrid.F90 $(MISC_DIR)/ioFileMod.F90 $(N84EQUIV_DIR)/radgrid.F90 $(MAIN_DIR)/input.F90 $(MISC_DIR)/physconst.F90 $(MAIN_DIR)/output.F90 $(N84EQUIV_DIR)/cloud.F90 $(N84EQUIV_DIR)/kabs.F90 $(MAIN_DIR)/mcica_random_numbers.F90 $(MAIN_DIR)/mcica.F90 $(N84EQUIV_DIR)/rad_interp_mod.F90 $(MAIN_DIR)/rayleigh_data.F90 $(N84EQUIV_DIR)/calc_opd_mod.F90 $(N84EQUIV_DIR)/initialize_rad_mod_1D.F90 $(MAIN_DIR)/planck_mod.F90 $(N84EQUIV_DIR)/model_specific.F90 $(MAIN_DIR)/exo_init_ref.F90 $(MAIN_DIR)/exo_radiation_mod.F90 $(MAIN_DIR)/main.F90

n84equiv.exe : $(OBJS_N84EQUIV)
	$(FC) -o $@  $(OBJS_N84EQUIV) -L$(LIB_NETCDF) $(ARGS) -I$(MISC_DIR)
	$(RM) -r -f $(N84EQUIV_BLD)
	mkdir $(N84EQUIV_BLD)
	mv -f *.mod $(N84EQUIV_BLD)
	mv -f n84equiv.exe ../run


##
## n42h2o
##
OBJS_N42H2O := $(MISC_DIR)/shr_kind_mod.F90 $(SRC_DIR)/exoplanet_mod.F90 $(MAIN_DIR)/sys_rootdir.F90 $(MISC_DIR)/shr_const_mod.F90 $(MISC_DIR)/shr_sys_mod.F90 $(MISC_DIR)/spmd_utils.F90 $(MISC_DIR)/time_manager.F90 $(MISC_DIR)/wrap_nf.F90 $(MISC_DIR)/infnan.F90 $(MISC_DIR)/ppgrid.F90 $(MISC_DIR)/pmgrid.F90 $(MISC_DIR)/ioFileMod.F90 $(N42H2O_DIR)/radgrid.F90 $(MAIN_DIR)/input.F90 $(MISC_DIR)/physconst.F90 $(MAIN_DIR)/output.F90 $(N42H2O_DIR)/cloud.F90 $(N42H2O_DIR)/kabs.F90 $(MAIN_DIR)/mcica_random_numbers.F90 $(MAIN_DIR)/mcica.F90 $(N42H2O_DIR)/rad_interp_mod.F90 $(MAIN_DIR)/rayleigh_data.F90 $(N42H2O_DIR)/calc_opd_mod.F90 $(N42H2O_DIR)/initialize_rad_mod_1D.F90 $(MAIN_DIR)/planck_mod.F90 $(N28ARCHEAN_DIR)/model_specific.F90 $(MAIN_DIR)/exo_init_ref.F90 $(MAIN_DIR)/exo_radiation_mod.F90 $(MAIN_DIR)/main.F90

n42h2o.exe : $(OBJS_N42H2O) 
	$(FC) -o $@  $(OBJS_N42H2O) -L$(LIB_NETCDF) $(ARGS) -I$(MISC_DIR)
	$(RM) -r -f $(N42H2O_BLD)
	mkdir $(N42H2O_BLD)
	mv -f *.mod $(N42H2O_BLD)
	mv -f n42h2o.exe ../run


##
## n68h2o
##
OBJS_N68H2O := $(MISC_DIR)/shr_kind_mod.F90 $(SRC_DIR)/exoplanet_mod.F90 $(MAIN_DIR)/sys_rootdir.F90 $(MISC_DIR)/shr_const_mod.F90 $(MISC_DIR)/shr_sys_mod.F90 $(MISC_DIR)/spmd_utils.F90 $(MISC_DIR)/time_manager.F90 $(MISC_DIR)/wrap_nf.F90 $(MISC_DIR)/infnan.F90 $(MISC_DIR)/ppgrid.F90 $(MISC_DIR)/pmgrid.F90 $(MISC_DIR)/ioFileMod.F90 $(N68H2O_DIR)/radgrid.F90 $(MAIN_DIR)/input.F90 $(MISC_DIR)/physconst.F90 $(MAIN_DIR)/output.F90 $(N68H2O_DIR)/cloud.F90 $(N68H2O_DIR)/kabs.F90 $(MAIN_DIR)/mcica_random_numbers.F90 $(MAIN_DIR)/mcica.F90 $(N68H2O_DIR)/rad_interp_mod.F90 $(MAIN_DIR)/rayleigh_data.F90 $(N68H2O_DIR)/calc_opd_mod.F90 $(N68H2O_DIR)/initialize_rad_mod_1D.F90 $(MAIN_DIR)/planck_mod.F90 $(N28ARCHEAN_DIR)/model_specific.F90 $(MAIN_DIR)/exo_init_ref.F90 $(MAIN_DIR)/exo_radiation_mod.F90 $(MAIN_DIR)/main.F90

n68h2o.exe : $(OBJS_N68H2O) 
	$(FC) -o $@  $(OBJS_N68H2O) -L$(LIB_NETCDF) $(ARGS) -I$(MISC_DIR)
	$(RM) -r -f $(N68H2O_BLD)
	mkdir $(N68H2O_BLD)
	mv -f *.mod $(N68H2O_BLD)
	mv -f n68h2o.exe ../run


##
## n28archean, should be identical to legacy highco2 build
##
OBJS_N28ARCHEAN := $(MISC_DIR)/shr_kind_mod.F90 $(SRC_DIR)/exoplanet_mod.F90 $(MAIN_DIR)/sys_rootdir.F90 $(MISC_DIR)/shr_const_mod.F90 $(MISC_DIR)/shr_sys_mod.F90 $(MISC_DIR)/spmd_utils.F90 $(MISC_DIR)/time_manager.F90 $(MISC_DIR)/wrap_nf.F90 $(MISC_DIR)/infnan.F90 $(MISC_DIR)/ppgrid.F90 $(MISC_DIR)/pmgrid.F90 $(MISC_DIR)/ioFileMod.F90 $(N28ARCHEAN_DIR)/radgrid.F90 $(MAIN_DIR)/input.F90 $(MISC_DIR)/physconst.F90 $(MAIN_DIR)/output.F90 $(N28ARCHEAN_DIR)/cloud.F90 $(N28ARCHEAN_DIR)/kabs.F90 $(MAIN_DIR)/mcica_random_numbers.F90 $(MAIN_DIR)/mcica.F90 $(N28ARCHEAN_DIR)/rad_interp_mod.F90 $(MAIN_DIR)/rayleigh_data.F90 $(N28ARCHEAN_DIR)/calc_opd_mod.F90 $(N28ARCHEAN_DIR)/initialize_rad_mod_1D.F90 $(MAIN_DIR)/planck_mod.F90 $(N28ARCHEAN_DIR)/model_specific.F90 $(MAIN_DIR)/exo_init_ref.F90 $(MAIN_DIR)/exo_radiation_mod.F90 $(MAIN_DIR)/main.F90

n28archean.exe : $(OBJS_N28ARCHEAN) 
	$(FC) -o $@  $(OBJS_N28ARCHEAN) -L$(LIB_NETCDF) $(ARGS) -I$(MISC_DIR)
	$(RM) -r -f $(N28ARCHEAN_BLD)
	mkdir $(N28ARCHEAN_BLD)
	mv -f *.mod $(N28ARCHEAN_BLD)
	mv -f n28archean.exe ../run


#================================================#
#======= Experimental Builds ====================#
#==== Used for model development and testing ====#

##
## n68equiv_exp
##
OBJS_N68EQUIV_EXP := $(MISC_DIR)/shr_kind_mod.F90 $(N68EQUIV_EXP_DIR)/exoplanet_mod.F90 $(N68EQUIV_EXP_DIR)/shr/sys_rootdir.F90 $(MISC_DIR)/shr_const_mod.F90 $(MISC_DIR)/shr_sys_mod.F90 $(MISC_DIR)/spmd_utils.F90 $(MISC_DIR)/time_manager.F90 $(MISC_DIR)/wrap_nf.F90 $(MISC_DIR)/infnan.F90 $(MISC_DIR)/ppgrid.F90 $(MISC_DIR)/pmgrid.F90 $(MISC_DIR)/ioFileMod.F90 $(N68EQUIV_EXP_DIR)/shr/radgrid.F90 $(N68EQUIV_EXP_DIR)/input.F90 $(MISC_DIR)/physconst.F90 $(N68EQUIV_EXP_DIR)/output.F90 $(N68EQUIV_EXP_DIR)/shr/cloud.F90 $(N68EQUIV_EXP_DIR)/shr/kabs.F90 $(N68EQUIV_EXP_DIR)/shr/mcica_random_numbers.F90 $(N68EQUIV_EXP_DIR)/shr/mcica.F90 $(N68EQUIV_EXP_DIR)/shr/rad_interp_mod.F90 $(N68EQUIV_EXP_DIR)/shr/rayleigh_data.F90 $(N68EQUIV_EXP_DIR)/shr/calc_opd_mod.F90 $(N68EQUIV_EXP_DIR)/initialize_rad_mod_1D.F90 $(N68EQUIV_EXP_DIR)/shr/exo_init_ref.F90 $(N68EQUIV_EXP_DIR)/shr/exo_radiation_mod.F90 $(N68EQUIV_EXP_DIR)/main.F90

n68equiv_exp.exe : $(OBJS_N68EQUIV_EXP)
	$(FC) -o $@  $(OBJS_N68EQUIV_EXP) -L$(LIB_NETCDF) $(ARGS) -I$(MISC_DIR)
	$(RM) -r -f $(N68EQUIV_EXP_BLD)
	mkdir $(N68EQUIV_EXP_BLD)
	mv -f *.mod $(N68EQUIV_EXP_BLD)
	mv -f n68equiv_exp.exe ../run



## ==============================================================================

#all        
n28archean       : n28archean.exe
n42h2o           : n42h2o.exe
n68h2o           : n68h2o.exe
n68equiv         : n68equiv.exe
n84equiv         : n84equiv.exe
n68equiv_exp     : n68equiv_exp.exe

clean: 
	$(RM) -f *.o *.mod *.exe
	$(RM) -f ../run/*.exe 
	$(RM) -r -f $(N28ARCHEAN_BLD)
	$(RM) -r -f $(N42H2O_BLD)
	$(RM) -r -f $(N68H2O_BLD)
	$(RM) -r -f $(N68EQUIV_BLD)
	$(RM) -r -f $(N84EQUIV_BLD)
	$(RM) -r -f $(N68EQUIV_EXP_BLD)

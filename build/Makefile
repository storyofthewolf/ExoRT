
#ifeq ($(LIB_NETCDF),$(null))

# netcdf directory for Hyak
#LIB_NETCDF := /sw/netcdf-fortran+c-4.4.1_icc-15.0/lib
# NOTE for hyak add "-lnetcdff"

# netcdf directory for Summit
#LIB_NETCDF := /projects/btoon/Libraries/netcdf-3.6.2-ifort12.0/lib       

# netcdf directory for discover
LIB_NETCDF := /gpfsm/dulocal/sles11/other/netcdf/3.6.2_intel-14.0.3.174/lib

#endif

# directories for up to date builds
MISC_DIR = ../source/src.misc
MAIN_DIR = ../source/src.main
N28ARCHEAN_DIR = ../source/src.n28archean
N42H2O_DIR = ../source/src.n42h2o
N42H2O_16GPT_DIR = ../source/experimental/src.n42h2o_16gpt
N42H2O_32GPT_DIR = ../source/experimental/src.n42h2o_32gpt
N68H2O_DIR = ../source/src.n68h2o
N68H2O_16GPT_DIR = ../source/experimental/src.n68h2o_16gpt
N68EQUIV_DIR = ../source/src.n68equiv
N68EQUIV_KV_DIR = ../source/experimental/src.n68equiv_kv
N68EQUIV_EXP_DIR = ../source/experimental/src.n68equiv_exp
N73EQUIV_DIR = ../source/experimental/src.n73equiv
N84EQUIV_DIR = ../source/src.n84equiv
N84EQUIV_EXP_DIR = ../source/experimental/src.n84equiv_exp

# legacy directories, old stuff, stuff not really used
MAIN_DIR_LEGACY = ../source/legacy/src.main_old
HIGHCO2_DIR  = ../source/legacy/src.highco2
H2ON28DIR = ../source/src.h2on28
H2ON35DIR = ../source/src.h2on35
H2ON68DIR = ../source/src.h2on68

# build directories
MAIN_BLD = bld.main
N28ARCHEAN_BLD = bld.n28archean
N42H2O_BLD = bld.n42h2o
N42H2O_16GPT_BLD = bld.n42h2o_16gpt
N42H2O_32GPT_BLD = bld.n42h2o_32gpt
N68H2O_BLD = bld.n68h2o
N68H2O_16GPT_BLD = bld.n68h2o_16gpt
N63CO2H2O_BLD = bld.n63co2h2o
N68EQUIV_BLD = bld.n68equiv
N68EQUIV_EXP_BLD = bld.n68equiv_exp
N68EQUIV_KV_BLD = bld.n68equiv_kv
N73EQUIV_BLD = bld.n73equiv
N84EQUIV_BLD = bld.n84equiv
N84EQUIV_EXP_BLD = bld.n84equiv_exp
HIGHCO2_BLD = bld.highco2


#extras
H2ON42_BLD = bld.h2on42
H2ON28_BLD = bld.h2on28
H2ON35_BLD = bld.h2on35



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
## n28archean, should be identical to legacy highco2 build
##
OBJS_N28ARCHEAN := $(MISC_DIR)/shr_kind_mod.F90 $(MAIN_DIR)/exoplanet_mod_1D.F90 $(MAIN_DIR)/sys_rootdir.F90 $(MISC_DIR)/shr_const_mod.F90 $(MISC_DIR)/shr_sys_mod.F90 $(MISC_DIR)/spmd_utils.F90 $(MISC_DIR)/time_manager.F90 $(MISC_DIR)/wrap_nf.F90 $(MISC_DIR)/infnan.F90 $(MISC_DIR)/ppgrid.F90 $(MISC_DIR)/pmgrid.F90 $(MISC_DIR)/ioFileMod.F90 $(N28ARCHEAN_DIR)/radgrid.F90 $(MAIN_DIR)/input.F90 $(MISC_DIR)/physconst.F90 $(MAIN_DIR)/output.F90 $(N28ARCHEAN_DIR)/cloud.F90 $(N28ARCHEAN_DIR)/kabs.F90 $(MAIN_DIR)/mcica_random_numbers.F90 $(MAIN_DIR)/mcica.F90 $(N28ARCHEAN_DIR)/rad_interp_mod.F90 $(MAIN_DIR)/rayleigh_data.F90 $(N28ARCHEAN_DIR)/calc_opd_mod.F90 $(N28ARCHEAN_DIR)/initialize_rad_mod_1D.F90 $(N28ARCHEAN_DIR)/exo_init_ref.F90 $(MAIN_DIR)/exo_radiation_mod.F90 $(MAIN_DIR)/main.F90

n28archean.exe : $(OBJS_N28ARCHEAN) 
	$(FC) -o $@  $(OBJS_N28ARCHEAN) -L$(LIB_NETCDF) -lnetcdf -I$(MISC_DIR)
	$(RM) -r -f $(N28ARCHEAN_BLD)
	mkdir $(N28ARCHEAN_BLD)
	mv -f *.mod $(N28ARCHEAN_BLD)
	mv -f n28archean.exe ../run


##
## n42h2o
##
OBJS_N42H2O := $(MISC_DIR)/shr_kind_mod.F90 $(MAIN_DIR)/exoplanet_mod_1D.F90 $(MAIN_DIR)/sys_rootdir.F90 $(MISC_DIR)/shr_const_mod.F90 $(MISC_DIR)/shr_sys_mod.F90 $(MISC_DIR)/spmd_utils.F90 $(MISC_DIR)/time_manager.F90 $(MISC_DIR)/wrap_nf.F90 $(MISC_DIR)/infnan.F90 $(MISC_DIR)/ppgrid.F90 $(MISC_DIR)/pmgrid.F90 $(MISC_DIR)/ioFileMod.F90 $(N42H2O_DIR)/radgrid.F90 $(MAIN_DIR)/input.F90 $(MISC_DIR)/physconst.F90 $(MAIN_DIR)/output.F90 $(N42H2O_DIR)/cloud.F90 $(N42H2O_DIR)/kabs.F90 $(MAIN_DIR)/mcica_random_numbers.F90 $(MAIN_DIR)/mcica.F90 $(N42H2O_DIR)/rad_interp_mod.F90 $(MAIN_DIR)/rayleigh_data.F90 $(N42H2O_DIR)/calc_opd_mod.F90 $(N42H2O_DIR)/initialize_rad_mod_1D.F90 $(N42H2O_DIR)/exo_init_ref.F90 $(MAIN_DIR)/exo_radiation_mod.F90 $(MAIN_DIR)/main.F90

n42h2o.exe : $(OBJS_N42H2O) 
	$(FC) -o $@  $(OBJS_N42H2O) -L$(LIB_NETCDF) -lnetcdf -I$(MISC_DIR)
	$(RM) -r -f $(N42H2O_BLD)
	mkdir $(N42H2O_BLD)
	mv -f *.mod $(N42H2O_BLD)
	mv -f n42h2o.exe ../run


##
## n68h2o
##
OBJS_N68H2O := $(MISC_DIR)/shr_kind_mod.F90 $(MAIN_DIR)/exoplanet_mod_1D.F90 $(MAIN_DIR)/sys_rootdir.F90 $(MISC_DIR)/shr_const_mod.F90 $(MISC_DIR)/shr_sys_mod.F90 $(MISC_DIR)/spmd_utils.F90 $(MISC_DIR)/time_manager.F90 $(MISC_DIR)/wrap_nf.F90 $(MISC_DIR)/infnan.F90 $(MISC_DIR)/ppgrid.F90 $(MISC_DIR)/pmgrid.F90 $(MISC_DIR)/ioFileMod.F90 $(N68H2O_DIR)/radgrid.F90 $(MAIN_DIR)/input.F90 $(MISC_DIR)/physconst.F90 $(MAIN_DIR)/output.F90 $(N68H2O_DIR)/cloud.F90 $(N68H2O_DIR)/kabs.F90 $(MAIN_DIR)/mcica_random_numbers.F90 $(MAIN_DIR)/mcica.F90 $(N68H2O_DIR)/rad_interp_mod.F90 $(MAIN_DIR)/rayleigh_data.F90 $(N68H2O_DIR)/calc_opd_mod.F90 $(N68H2O_DIR)/initialize_rad_mod_1D.F90 $(N68H2O_DIR)/exo_init_ref.F90 $(MAIN_DIR)/exo_radiation_mod.F90 $(MAIN_DIR)/main.F90

n68h2o.exe : $(OBJS_N68H2O) 
	$(FC) -o $@  $(OBJS_N68H2O) -L$(LIB_NETCDF) -lnetcdf -I$(MISC_DIR)
	$(RM) -r -f $(N68H2O_BLD)
	mkdir $(N68H2O_BLD)
	mv -f *.mod $(N68H2O_BLD)
	mv -f n68h2o.exe ../run


##
## n68equiv
## 
OBJS_N68EQUIV := $(MISC_DIR)/shr_kind_mod.F90 $(MAIN_DIR)/exoplanet_mod_1D.F90 $(MAIN_DIR)/sys_rootdir.F90 $(MISC_DIR)/shr_const_mod.F90 $(MISC_DIR)/shr_sys_mod.F90 $(MISC_DIR)/spmd_utils.F90 $(MISC_DIR)/time_manager.F90 $(MISC_DIR)/wrap_nf.F90 $(MISC_DIR)/infnan.F90 $(MISC_DIR)/ppgrid.F90 $(MISC_DIR)/pmgrid.F90 $(MISC_DIR)/ioFileMod.F90 $(N68EQUIV_DIR)/radgrid.F90 $(MAIN_DIR)/input.F90 $(MISC_DIR)/physconst.F90 $(MAIN_DIR)/output.F90 $(N68EQUIV_DIR)/cloud.F90 $(N68EQUIV_DIR)/kabs.F90 $(MAIN_DIR)/mcica_random_numbers.F90 $(MAIN_DIR)/mcica.F90 $(N68EQUIV_DIR)/rad_interp_mod.F90 $(MAIN_DIR)/rayleigh_data.F90 $(N68EQUIV_DIR)/calc_opd_mod.F90 $(N68EQUIV_DIR)/initialize_rad_mod_1D.F90 $(N68EQUIV_DIR)/exo_init_ref.F90 $(MAIN_DIR)/exo_radiation_mod.F90 $(MAIN_DIR)/main.F90

n68equiv.exe : $(OBJS_N68EQUIV)
	$(FC) -o $@  $(OBJS_N68EQUIV) -L$(LIB_NETCDF) -lnetcdf -I$(MISC_DIR)
	$(RM) -r -f $(N68EQUIV_BLD)
	mkdir $(N68EQUIV_BLD)
	mv -f *.mod $(N68EQUIV_BLD)
	mv -f n68equiv.exe ../run


##
## n84equiv
##
OBJS_N84EQUIV := $(MISC_DIR)/shr_kind_mod.F90 $(MAIN_DIR)/exoplanet_mod_1D.F90 $(MAIN_DIR)/sys_rootdir.F90 $(MISC_DIR)/shr_const_mod.F90 $(MISC_DIR)/shr_sys_mod.F90 $(MISC_DIR)/spmd_utils.F90 $(MISC_DIR)/time_manager.F90 $(MISC_DIR)/wrap_nf.F90 $(MISC_DIR)/infnan.F90 $(MISC_DIR)/ppgrid.F90 $(MISC_DIR)/pmgrid.F90 $(MISC_DIR)/ioFileMod.F90 $(N84EQUIV_DIR)/radgrid.F90 $(MAIN_DIR)/input.F90 $(MISC_DIR)/physconst.F90 $(MAIN_DIR)/output.F90 $(N84EQUIV_DIR)/cloud.F90 $(N84EQUIV_DIR)/kabs.F90 $(MAIN_DIR)/mcica_random_numbers.F90 $(MAIN_DIR)/mcica.F90 $(N84EQUIV_DIR)/rad_interp_mod.F90 $(MAIN_DIR)/rayleigh_data.F90 $(N84EQUIV_DIR)/calc_opd_mod.F90 $(N84EQUIV_DIR)/initialize_rad_mod_1D.F90 $(N84EQUIV_DIR)/exo_init_ref.F90 $(MAIN_DIR)/exo_radiation_mod.F90 $(MAIN_DIR)/main.F90

n84equiv.exe : $(OBJS_N84EQUIV)
	$(FC) -o $@  $(OBJS_N84EQUIV) -L$(LIB_NETCDF) -lnetcdf -I$(MISC_DIR)
	$(RM) -r -f $(N84EQUIV_BLD)
	mkdir $(N84EQUIV_BLD)
	mv -f *.mod $(N84EQUIV_BLD)
	mv -f n84equiv.exe ../run


#================================================#
#======= Experimental Builds ====================#
#==== Used for model development and testing ====#

##
## n73equiv
## 
OBJS_N73EQUIV := $(MISC_DIR)/shr_kind_mod.F90 $(MAIN_DIR)/exoplanet_mod_1D.F90 $(MAIN_DIR)/sys_rootdir.F90 $(MISC_DIR)/shr_const_mod.F90 $(MISC_DIR)/shr_sys_mod.F90 $(MISC_DIR)/spmd_utils.F90 $(MISC_DIR)/time_manager.F90 $(MISC_DIR)/wrap_nf.F90 $(MISC_DIR)/infnan.F90 $(MISC_DIR)/ppgrid.F90 $(MISC_DIR)/pmgrid.F90 $(MISC_DIR)/ioFileMod.F90 $(N73EQUIV_DIR)/radgrid.F90 $(MAIN_DIR)/input.F90 $(MISC_DIR)/physconst.F90 $(MAIN_DIR)/output.F90 $(N73EQUIV_DIR)/cloud.F90 $(N73EQUIV_DIR)/kabs.F90 $(MAIN_DIR)/mcica_random_numbers.F90 $(MAIN_DIR)/mcica.F90 $(N73EQUIV_DIR)/rad_interp_mod.F90 $(MAIN_DIR)/rayleigh_data.F90 $(N73EQUIV_DIR)/calc_opd_mod.F90 $(N73EQUIV_DIR)/initialize_rad_mod_1D.F90 $(N73EQUIV_DIR)/exo_init_ref.F90 $(MAIN_DIR)/exo_radiation_mod.F90 $(MAIN_DIR)/main.F90

n73equiv.exe : $(OBJS_N73EQUIV)
	$(FC) -o $@  $(OBJS_N73EQUIV) -L$(LIB_NETCDF) -lnetcdf -I$(MISC_DIR)
	$(RM) -r -f $(N73EQUIV_BLD)
	mkdir $(N73EQUIV_BLD)
	mv -f *.mod $(N73EQUIV_BLD)
	mv -f n73equiv.exe ../run


##
## n68equiv_kv
##
OBJS_N68EQUIV_KV := $(MISC_DIR)/shr_kind_mod.F90 $(MAIN_DIR)/exoplanet_mod_1D.F90 $(MAIN_DIR)/sys_rootdir.F90 $(MISC_DIR)/shr_const_mod.F90 $(MISC_DIR)/shr_sys_mod.F90 $(MISC_DIR)/spmd_utils.F90 $(MISC_DIR)/time_manager.F90 $(MISC_DIR)/wrap_nf.F90 $(MISC_DIR)/infnan.F90 $(MISC_DIR)/ppgrid.F90 $(MISC_DIR)/pmgrid.F90 $(MISC_DIR)/ioFileMod.F90 $(N68EQUIV_KV_DIR)/radgrid.F90 $(MAIN_DIR)/input.F90 $(MISC_DIR)/physconst.F90 $(MAIN_DIR)/output.F90 $(N68EQUIV_KV_DIR)/cloud.F90 $(N68EQUIV_KV_DIR)/kabs.F90 $(MAIN_DIR)/mcica_random_numbers.F90 $(MAIN_DIR)/mcica.F90 $(N68EQUIV_KV_DIR)/rad_interp_mod.F90 $(MAIN_DIR)/rayleigh_data.F90 $(N68EQUIV_KV_DIR)/calc_opd_mod.F90 $(N68EQUIV_KV_DIR)/initialize_rad_mod_1D.F90 $(N68EQUIV_KV_DIR)/exo_init_ref.F90 $(MAIN_DIR)/exo_radiation_mod.F90 $(MAIN_DIR)/main.F90

n68equiv_kv.exe : $(OBJS_N68EQUIV_KV)
	$(FC) -o $@  $(OBJS_N68EQUIV_KV) -L$(LIB_NETCDF) -lnetcdf -I$(MISC_DIR)
	$(RM) -r -f $(N68EQUIV_KV_BLD)
	mkdir $(N68EQUIV_KV_BLD)
	mv -f *.mod $(N68EQUIV_KV_BLD)
	mv -f n68equiv_kv.exe ../run


##
## n68equiv_exp
##
OBJS_N68EQUIV_EXP := $(MISC_DIR)/shr_kind_mod.F90 $(MAIN_DIR)/exoplanet_mod_1D.F90 $(N68EQUIV_EXP_DIR)/sys_rootdir.F90 $(MISC_DIR)/shr_const_mod.F90 $(MISC_DIR)/shr_sys_mod.F90 $(MISC_DIR)/spmd_utils.F90 $(MISC_DIR)/time_manager.F90 $(MISC_DIR)/wrap_nf.F90 $(MISC_DIR)/infnan.F90 $(MISC_DIR)/ppgrid.F90 $(MISC_DIR)/pmgrid.F90 $(MISC_DIR)/ioFileMod.F90 $(N68EQUIV_EXP_DIR)/radgrid.F90 $(N68EQUIV_EXP_DIR)/input.F90 $(MISC_DIR)/physconst.F90 $(N68EQUIV_EXP_DIR)/output.F90 $(N68EQUIV_EXP_DIR)/cloud.F90 $(N68EQUIV_EXP_DIR)/kabs.F90 $(N68EQUIV_EXP_DIR)/mcica_random_numbers.F90 $(N68EQUIV_EXP_DIR)/mcica.F90 $(N68EQUIV_EXP_DIR)/rad_interp_mod.F90 $(N68EQUIV_EXP_DIR)/rayleigh_data.F90 $(N68EQUIV_EXP_DIR)/calc_opd_mod.F90 $(N68EQUIV_EXP_DIR)/initialize_rad_mod_1D.F90 $(N68EQUIV_EXP_DIR)/exo_init_ref.F90 $(N68EQUIV_EXP_DIR)/exo_radiation_mod.F90 $(N68EQUIV_EXP_DIR)/main.F90

n68equiv_exp.exe : $(OBJS_N68EQUIV_EXP)
	$(FC) -o $@  $(OBJS_N68EQUIV_EXP) -L$(LIB_NETCDF) -lnetcdf -I$(MISC_DIR)
	$(RM) -r -f $(N68EQUIV_EXP_BLD)
	mkdir $(N68EQUIV_EXP_BLD)
	mv -f *.mod $(N68EQUIV_EXP_BLD)
	mv -f n68equiv_exp.exe ../run


##
## n84equiv_exp
##
OBJS_N84EQUIV_EXP := $(MISC_DIR)/shr_kind_mod.F90 $(MAIN_DIR)/exoplanet_mod_1D.F90 $(MAIN_DIR)/sys_rootdir.F90 $(MISC_DIR)/shr_const_mod.F90 $(MISC_DIR)/shr_sys_mod.F90 $(MISC_DIR)/spmd_utils.F90 $(MISC_DIR)/time_manager.F90 $(MISC_DIR)/wrap_nf.F90 $(MISC_DIR)/infnan.F90 $(MISC_DIR)/ppgrid.F90 $(MISC_DIR)/pmgrid.F90 $(MISC_DIR)/ioFileMod.F90 $(N84EQUIV_EXP_DIR)/radgrid.F90 $(MAIN_DIR)/input.F90 $(MISC_DIR)/physconst.F90 $(MAIN_DIR)/output.F90 $(N84EQUIV_EXP_DIR)/cloud.F90 $(N84EQUIV_EXP_DIR)/kabs.F90 $(MAIN_DIR)/mcica_random_numbers.F90 $(MAIN_DIR)/mcica.F90 $(N84EQUIV_EXP_DIR)/rad_interp_mod.F90 $(MAIN_DIR)/rayleigh_data.F90 $(N84EQUIV_EXP_DIR)/calc_opd_mod.F90 $(N84EQUIV_EXP_DIR)/initialize_rad_mod_1D.F90 $(N84EQUIV_EXP_DIR)/exo_init_ref.F90 $(N84EQUIV_EXP_DIR)/exo_radiation_mod.F90 $(MAIN_DIR)/main.F90

n84equiv_exp.exe : $(OBJS_N84EQUIV_EXP)
	$(FC) -o $@  $(OBJS_N84EQUIV_EXP) -L$(LIB_NETCDF) -lnetcdf -I$(MISC_DIR)
	$(RM) -r -f $(N84EQUIV_EXP_BLD)
	mkdir $(N84EQUIV_EXP_BLD)
	mv -f *.mod $(N84EQUIV_EXP_BLD)
	mv -f n84equiv_exp.exe ../run


##
## n42h2o_16gpt
##
OBJS_N42H2O_16GPT := $(MISC_DIR)/shr_kind_mod.F90 $(MAIN_DIR)/exoplanet_mod_1D.F90 $(MAIN_DIR)/sys_rootdir.F90 $(MISC_DIR)/shr_const_mod.F90 $(MISC_DIR)/shr_sys_mod.F90 $(MISC_DIR)/spmd_utils.F90 $(MISC_DIR)/time_manager.F90 $(MISC_DIR)/wrap_nf.F90 $(MISC_DIR)/infnan.F90 $(MISC_DIR)/ppgrid.F90 $(MISC_DIR)/pmgrid.F90 $(MISC_DIR)/ioFileMod.F90 $(N42H2O_16GPT_DIR)/radgrid.F90 $(MAIN_DIR)/input.F90 $(MISC_DIR)/physconst.F90 $(MAIN_DIR)/output.F90 $(N42H2O_16GPT_DIR)/cloud.F90 $(N42H2O_16GPT_DIR)/kabs.F90 $(MAIN_DIR)/mcica_random_numbers.F90 $(MAIN_DIR)/mcica.F90 $(N42H2O_16GPT_DIR)/rad_interp_mod.F90 $(MAIN_DIR)/rayleigh_data.F90 $(N42H2O_16GPT_DIR)/calc_opd_mod.F90 $(N42H2O_16GPT_DIR)/initialize_rad_mod_1D.F90 $(N42H2O_16GPT_DIR)/exo_init_ref.F90 $(MAIN_DIR)/exo_radiation_mod.F90 $(MAIN_DIR)/main.F90

n42h2o_16gpt.exe : $(OBJS_N42H2O_16GPT) 
	$(FC) -o $@  $(OBJS_N42H2O_16GPT) -L$(LIB_NETCDF) -lnetcdf -I$(MISC_DIR)
	$(RM) -r -f $(N42H2O_16GPT_BLD)
	mkdir $(N42H2O_16GPT_BLD)
	mv -f *.mod $(N42H2O_16GPT_BLD)
	mv -f n42h2o_16gpt.exe ../run


##
## n42h2o_32gpt
##
OBJS_N42H2O_32GPT := $(MISC_DIR)/shr_kind_mod.F90 $(MAIN_DIR)/exoplanet_mod_1D.F90 $(MAIN_DIR)/sys_rootdir.F90 $(MISC_DIR)/shr_const_mod.F90 $(MISC_DIR)/shr_sys_mod.F90 $(MISC_DIR)/spmd_utils.F90 $(MISC_DIR)/time_manager.F90 $(MISC_DIR)/wrap_nf.F90 $(MISC_DIR)/infnan.F90 $(MISC_DIR)/ppgrid.F90 $(MISC_DIR)/pmgrid.F90 $(MISC_DIR)/ioFileMod.F90 $(N42H2O_32GPT_DIR)/radgrid.F90 $(MAIN_DIR)/input.F90 $(MISC_DIR)/physconst.F90 $(MAIN_DIR)/output.F90 $(N42H2O_32GPT_DIR)/cloud.F90 $(N42H2O_32GPT_DIR)/kabs.F90 $(MAIN_DIR)/mcica_random_numbers.F90 $(MAIN_DIR)/mcica.F90 $(N42H2O_32GPT_DIR)/rad_interp_mod.F90 $(MAIN_DIR)/rayleigh_data.F90 $(N42H2O_32GPT_DIR)/calc_opd_mod.F90 $(N42H2O_32GPT_DIR)/initialize_rad_mod_1D.F90 $(N42H2O_32GPT_DIR)/exo_init_ref.F90 $(MAIN_DIR)/exo_radiation_mod.F90 $(MAIN_DIR)/main.F90

n42h2o_32gpt.exe : $(OBJS_N42H2O_32GPT) 
	$(FC) -o $@  $(OBJS_N42H2O_32GPT) -L$(LIB_NETCDF) -lnetcdf -I$(MISC_DIR)
	$(RM) -r -f $(N42H2O_32GPT_BLD)
	mkdir $(N42H2O_32GPT_BLD)
	mv -f *.mod $(N42H2O_32GPT_BLD)
	mv -f n42h2o_32gpt.exe ../run


##
## n68h2o_16gpt
##
OBJS_N68H2O_16GPT := $(MISC_DIR)/shr_kind_mod.F90 $(MAIN_DIR)/exoplanet_mod_1D.F90 $(MAIN_DIR)/sys_rootdir.F90 $(MISC_DIR)/shr_const_mod.F90 $(MISC_DIR)/shr_sys_mod.F90 $(MISC_DIR)/spmd_utils.F90 $(MISC_DIR)/time_manager.F90 $(MISC_DIR)/wrap_nf.F90 $(MISC_DIR)/infnan.F90 $(MISC_DIR)/ppgrid.F90 $(MISC_DIR)/pmgrid.F90 $(MISC_DIR)/ioFileMod.F90 $(N68H2O_16GPT_DIR)/radgrid.F90 $(MAIN_DIR)/input.F90 $(MISC_DIR)/physconst.F90 $(MAIN_DIR)/output.F90 $(N68H2O_16GPT_DIR)/cloud.F90 $(N68H2O_16GPT_DIR)/kabs.F90 $(MAIN_DIR)/mcica_random_numbers.F90 $(MAIN_DIR)/mcica.F90 $(N68H2O_16GPT_DIR)/rad_interp_mod.F90 $(MAIN_DIR)/rayleigh_data.F90 $(N68H2O_16GPT_DIR)/calc_opd_mod.F90 $(N68H2O_16GPT_DIR)/initialize_rad_mod_1D.F90 $(N68H2O_16GPT_DIR)/exo_init_ref.F90 $(MAIN_DIR)/exo_radiation_mod.F90 $(MAIN_DIR)/main.F90

n68h2o_16gpt.exe : $(OBJS_N68H2O_16GPT) 
	$(FC) -o $@  $(OBJS_N68H2O_16GPT) -L$(LIB_NETCDF) -lnetcdf -I$(MISC_DIR)
	$(RM) -r -f $(N68H2O_16GPT_BLD)
	mkdir $(N68H2O_16GPT_BLD)
	mv -f *.mod $(N68H2O_16GPT_BLD)
	mv -f n68h2o_16gpt.exe ../run



#==============================
#==== Legacy ====
# older builds

#Rules for each executable that could be built
##
## "highco2" legacy build, should give the same results as "archean"
##
OBJS_HIGHCO2 := $(MAIN_DIR_LEGACY)/shr_kind_mod.F90 $(MAIN_DIR_LEGACY)/shr_const_mod.F90 $(MAIN_DIR_LEGACY)/shr_sys_mod.F90 $(MAIN_DIR_LEGACY)/wrap_nf.F90 $(MAIN_DIR_LEGACY)/infnan.F90 $(MAIN_DIR_LEGACY)/ppgrid.F90 $(MAIN_DIR_LEGACY)/ioFileMod.F90 $(HIGHCO2_DIR)/radgrid.F90 $(MAIN_DIR_LEGACY)/input.F90 $(MAIN_DIR_LEGACY)/physconst.F90 $(MAIN_DIR_LEGACY)/output.F90 $(HIGHCO2_DIR)/kabs.F90 $(HIGHCO2_DIR)/solar.F90 $(MAIN_DIR_LEGACY)/mcica_random_numbers.F90 $(MAIN_DIR_LEGACY)/mcica.F90 $(HIGHCO2_DIR)/rad_interp_mod.F90 $(MAIN_DIR_LEGACY)/rayleigh_data.F90 $(HIGHCO2_DIR)/calc_opd_mod.F90 $(HIGHCO2_DIR)/radiation.F90 $(HIGHCO2_DIR)/initialize_rad_mod.F90 $(MAIN_DIR_LEGACY)/main.F90

highco2.exe : $(OBJS_HIGHCO2) 
	$(FC) -o $@  $(OBJS_HIGHCO2) -L$(LIB_NETCDF) -lnetcdf -I$(MAIN_DIR_LEGACY)
	$(RM) -r -f $(HIGHCO2_BLD)
	mkdir $(HIGHCO2_BLD)
	mv -f *.mod $(HIGHCO2_BLD)
	mv -f highco2.exe ../run


##
## h2o n28
##
OBJSH2ON28 := $(MAINDIR)/shr_kind_mod.F90 $(MAINDIR)/shr_const_mod.F90 $(MAINDIR)/shr_sys_mod.F90 $(MAINDIR)/wrap_nf.F90 $(MAINDIR)/infnan.F90 $(MAINDIR)/ppgrid.F90 $(MAINDIR)/ioFileMod.F90 $(H2ON28DIR)/radgrid.F90 $(MAINDIR)/input.F90 $(MAINDIR)/physconst.F90 $(MAINDIR)/output.F90 $(H2ON28DIR)/kabs.F90 $(H2ON28DIR)/solar.F90 $(MAINDIR)/mcica_random_numbers.F90 $(MAINDIR)/mcica.F90 $(H2ON28DIR)/rad_interp_mod.F90 $(MAINDIR)/rayleigh_data.F90 $(H2ON28DIR)/calc_opd_mod.F90 $(H2ON28DIR)/radiation.F90 $(H2ON28DIR)/initialize_rad_mod.F90 $(MAINDIR)/main.F90

h2on28.exe : $(OBJSH2ON28) 
	$(FC) -o $@  $(OBJSH2ON28) -L$(LIB_NETCDF) -lnetcdf -I$(MAINDIR)
	$(RM) -r -f $(H2ON28DIR)
	mkdir $(H2ON28BLD)
	mv -f *.mod $(H2ON28BLD)
	mv -f h2on28.exe run


##
## h2o n35
##
OBJSH2ON35 := $(MAINDIR)/shr_kind_mod.F90 $(MAINDIR)/shr_const_mod.F90 $(MAINDIR)/shr_sys_mod.F90 $(MAINDIR)/wrap_nf.F90 $(MAINDIR)/infnan.F90 $(MAINDIR)/ppgrid.F90 $(MAINDIR)/ioFileMod.F90 $(H2ON35DIR)/radgrid.F90 $(MAINDIR)/input.F90 $(MAINDIR)/physconst.F90 $(MAINDIR)/output.F90 $(H2ON35DIR)/kabs.F90 $(H2ON35DIR)/solar.F90 $(MAINDIR)/mcica_random_numbers.F90 $(MAINDIR)/mcica.F90 $(H2ON35DIR)/rad_interp_mod.F90 $(MAINDIR)/rayleigh_data.F90 $(H2ON35DIR)/calc_opd_mod.F90 $(H2ON35DIR)/radiation.F90 $(H2ON35DIR)/initialize_rad_mod.F90 $(MAINDIR)/main.F90

h2on35.exe : $(OBJSH2ON35) 
	$(FC) -o $@  $(OBJSH2ON35) -L$(LIB_NETCDF) -lnetcdf -I$(MAINDIR)
	$(RM) -r -f $(blddir)
	mkdir $(blddir)
	mkdir $(blddir)/$(H2ON35BLD)
	mv -f *.mod $(blddir)/$(H2ON35BLD)
	mv -f h2on35.exe run


##
## h2o n42
##
OBJSH2ON42 := $(MAINDIR)/shr_kind_mod.F90 $(MAINDIR)/shr_const_mod.F90 $(MAINDIR)/shr_sys_mod.F90 $(MAINDIR)/wrap_nf.F90 $(MAINDIR)/infnan.F90 $(MAINDIR)/ppgrid.F90 $(MAINDIR)/ioFileMod.F90 $(H2ON42DIR)/radgrid.F90 $(MAINDIR)/input.F90 $(MAINDIR)/physconst.F90 $(MAINDIR)/output.F90 $(H2ON42DIR)/kabs.F90 $(H2ON42DIR)/solar.F90 $(MAINDIR)/mcica_random_numbers.F90 $(MAINDIR)/mcica.F90 $(H2ON42DIR)/rad_interp_mod.F90 $(MAINDIR)/rayleigh_data.F90 $(H2ON42DIR)/calc_opd_mod.F90 $(H2ON42DIR)/radiation.F90 $(H2ON42DIR)/initialize_rad_mod.F90 $(MAINDIR)/main.F90

h2on42.exe : $(OBJSH2ON42) 
	$(FC) -o $@  $(OBJSH2ON42) -L$(LIB_NETCDF) -lnetcdf -I$(MAINDIR)
	$(RM) -r -f $(blddir)
	mkdir $(blddir)
	mkdir $(blddir)/$(H2ON42BLD)
	mv -f *.mod $(blddir)/$(H2ON42BLD)
	mv -f h2on42.exe run


##
## h2o n68
##
OBJSH2ON68 := $(MAINDIR)/shr_kind_mod.F90 $(MAINDIR)/shr_const_mod.F90 $(MAINDIR)/shr_sys_mod.F90 $(MAINDIR)/wrap_nf.F90 $(MAINDIR)/infnan.F90 $(MAINDIR)/ppgrid.F90 $(MAINDIR)/ioFileMod.F90 $(H2ON68DIR)/radgrid.F90 $(MAINDIR)/input.F90 $(MAINDIR)/physconst.F90 $(MAINDIR)/output.F90 $(H2ON68DIR)/kabs.F90 $(H2ON68DIR)/solar.F90 $(MAINDIR)/mcica_random_numbers.F90 $(MAINDIR)/mcica.F90 $(H2ON68DIR)/rad_interp_mod.F90 $(MAINDIR)/rayleigh_data.F90 $(H2ON68DIR)/calc_opd_mod.F90 $(H2ON68DIR)/radiation.F90 $(H2ON68DIR)/initialize_rad_mod.F90 $(MAINDIR)/main.F90

h2on68.exe : $(OBJSH2ON68) 
	$(FC) -o $@  $(OBJSH2ON68) -L$(LIB_NETCDF) -lnetcdf -I$(MAINDIR)
	$(RM) -r -f $(blddir)
	mkdir $(blddir)
	mkdir $(blddir)/$(H2ON68BLD)
	mv -f *.mod $(blddir)/$(H2ON68BLD)
	mv -f h2on68.exe


## ==============================================================================

#all        
n28archean       : n28archean.exe
n42h2o           : n42h2o.exe
n42h2o_16gpt     : n42h2o_16gpt.exe
n42h2o_32gpt     : n42h2o_32gpt.exe
n68h2o           : n68h2o.exe
n68h2o_16gpt     : n68h2o_16gpt.exe
n63co2h2o        : n63co2h2o.exe
n68equiv         : n68equiv.exe
n68equiv_kv      : n68equiv_kv.exe
n68equiv_exp     : n68equiv_exp.exe
n73equiv         : n73equiv.exe
n84equiv         : n84equiv.exe
n84equiv_exp     : n84equiv_exp.exe

#legacy
highco2      :  highco2.exe
h2on28       :  h2on28.exe
h2on35       :  h2on35.exe
h2on42       :  h2on42.exe
h2on68       :  h2on68.exe


#ox : oxygen.exe
clean: 
	$(RM) -f *.o *.mod *.exe
	$(RM) -f ../run/*.exe 
	$(RM) -r -f $(N63CO2H2O_BLD)
	$(RM) -r -f $(H2ON42_BLD)
	$(RM) -r -f $(H2ON35_BLD)
	$(RM) -r -f $(H2ON28_BLD)
	$(RM) -r -f $(HIGHCO2_BLD)
	$(RM) -r -f $(N28ARCHEAN_BLD)
	$(RM) -r -f $(N42H2O_BLD)
	$(RM) -r -f $(N42H2O_16GPT_BLD)
	$(RM) -r -f $(N42H2O_32GPT_BLD)
	$(RM) -r -f $(N68H2O_BLD)
	$(RM) -r -f $(N68H2O_16GPT_BLD)
	$(RM) -r -f $(N68EQUIV_BLD)
	$(RM) -r -f $(N68EQUIV_KV_BLD)
	$(RM) -r -f $(N68EQUIV_EXP_BLD)
	$(RM) -r -f $(N73EQUIV_BLD)
	$(RM) -r -f $(N84EQUIV_BLD)
	$(RM) -r -f $(N84EQUIV_EXP_BLD)
	

# Build/run the standalone C test for libexort (Stage D).
#
#   make        build test_exort_c
#   make run    build, then run against run/libexort + the TS300K_G2V
#               fixture/baseline (paths default relative to this directory)
#
# Requires the netCDF C library (nc-config on PATH) and a built
# run/libexort.dylib|.so (cd build && USER_FC=gfortran make libexort).

CC ?= cc
NC_CFLAGS := $(shell nc-config --cflags)
# rpath: anaconda/homebrew netCDF dylibs carry @rpath install names
NC_LIBS   := $(shell nc-config --libs) -Wl,-rpath,$(shell nc-config --libdir)

ifeq ($(shell uname -s),Darwin)
  SHLIB_EXT = .dylib
else
  SHLIB_EXT = .so
endif

test_exort_c : test_exort_c.c
	$(CC) -o $@ $< -DSHLIB_EXT='"$(SHLIB_EXT)"' $(NC_CFLAGS) $(NC_LIBS) -lm

run : test_exort_c
	./test_exort_c

clean:
	$(RM) test_exort_c

.PHONY: run clean

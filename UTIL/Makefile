# Legacy UTIL Makefile, kept for reference.
#
# Machine: "rotor", the Linux x86_64 server the original SAM-ASP-TUV was
# built on (the top-level Makefile tests $(HOSTNAME) for rotor). Evidence:
# the paths (/rotor/data/p1808/CMAQv5.0.1/... for NetCDF built with ifort;
# Intel-built HDF5 1.8.6 and zlib 1.2.5 under /usr/local), and the old
# prebuilt UTIL binaries (removed; most likely built with this file): SUSE
# Linux (glibc 2.11.1, gcc 4.3.4, i.e. SLES 11 era), Intel Fortran 15.0
# run-time library (2015), NetCDF built with Intel icc 12.0. Its project
# data later sat under /nas/data_cellar/rotor/data/ (see
# SRC_TUV/ASP/Makefile.out).
# These paths do not exist elsewhere: build the tools with build_util.sh.
#
# Makefile to compile conversion utilities. You need to set include and library paths for NETCDF
#

PLATFORM := $(shell uname -s)

# uncomment line for your platform or add your own:

#---------------------------------------------------------


  INC_NETCDF = /rotor/data/p1808/CMAQv5.0.1/lib/x86_64/ifort/netcdf/include
  LIB_NETCDF = /rotor/data/p1808/CMAQv5.0.1/lib/x86_64/ifort/netcdf/lib -lnetcdff -lnetcdf -lnetcdf
  HDF5_HL =  /usr/local/hdf5/1.8.6/intel/lib -lhdf5_hl -lhdf5
  Z =  /usr/local/zlib/1.2.5/intel/lib -lz -lm -lhdf5_hl -lhdf5 -lz

 # FF77 = ifort -c -fixed -extend_source
  #FF90 = ifort -c -free
  #FFLAGS = -I$(INC_NETCDF) -I$(INC_MPI) -g -O2 -pad 
  #LD = ifort
  #LDFLAGS = -L$(LIB_MPI) -L$(LIB_NETCDF) -L$(HDF5_HL) -L$(Z)

  #CC = icc -c -DLINUX
# AIX 
#ifeq ($(PLATFORM),AIX)

#LIB_NETCDF = /usr/local/xlf/lib
#INC_NETCDF = /usr/local/xlf/include

#FF = xlf90_r -qfixed -I$(INC_NETCDF) 
#LDFLAGS = -L$(LIB_NETCDF) -lnetcdf

#endif
#---------------------------------------------------------
# SGI IRIX64


#ifeq ($(PLATFORM),IRIX64)


#LIB_NETCDF = /usr/local/lib

#FF = f90 -fixedform -extend_source -I$(INC_NETCDF) 
#LDFLAGS = -L$(LIB_NETCDF) -lnetcdf

#endif
#---------------------------------------------------------
# Mac OS X 


#ifeq ($(PLATFORM),Darwin)


#INC_NETCDF      := /usr/local/netcdf/include
#LIB_NETCDF       := /usr/local/netcdf/lib

#FF = mpif90 -fixed -O3 -pad -extend_source -I$(INC_NETCDF)

#FF = mpif90 -fixed -extend_source -I$(INC_NETCDF)

#FF = xlf90 -qfixed -I$(INC_NETCDF) 

#FF = f95 -qfixed -kind=byte -gline -I$(INC_NETCDF) 

#FF = gfortran  -O3 -ffixed-form -ffixed-line-length-0 -I$(INC_NETCDF)

#LDFLAGS = -L${LIB_NETCDF} -lnetcdf


#endif
#---------------------------------------------------------
#
#ifeq ($(PLATFORM),Linux)



#  LIB_MPI = /usr/local/mpich/mpich2-1.5-intel11.1/lib
#  INC_MPI = /usr/local/mpich/mpich2-1.5-intel11.1/include

  #INC_NETCDF = /usr/local/netcdf/3.6.3-intel11.1/include
  #LIB_NETCDF = /usr/local/netcdf/3.6.3-intel11.1/lib  

  #FFLAGS = -O3 

  #FFLAGS += -I${INC_NETCDF}
  #LDFLAGS = -L${LIB_NETCDF} -lnetcdf
  #FFLAGS += -I${INC_MPI}
  #LDFLAGS += -L${LIB_MPI} -lmpich -lmpl
  #FF77= /opt/intel/Compiler/11.1/080/bin/intel64/ifort -c -fixed -extend_source
  #FF= /opt/intel/Compiler/11.1/080/bin/intel64/ifort 
  #FF90= /opt/intel/Compiler/11.1/080/bin/intel64/ifort -c
  #CC= /opt/intel/Compiler/11.1/080/bin/intel64/icc -c -DLINUX
  #LD= /opt/intel/Compiler/11.1/080/bin/intel64/ifort 



#---------------------------------------------------------
#
# Linux Cluster (seawulf@sunysb)
#

#INC_NETCDF = /nfs/user04/marat/local/include
#LIB_NETCDF = /nfs/user04/marat/local/lib


FF = ifort -O3 -fixed -extend_source -I${INC_NETCDF}
LDFLAGS = -L${LIB_NETCDF} -lnetcdf -L$(HDF5_HL) -L$(Z)

#endif
#---------------------------------------------------------

VPATH = ./SRC

all: bin2D2nc bin3D2nc 2Dbin2nc 2Dbin2nc_mean bin3D2nc_mean com3D2bin 2Dcom2nc 2Dcom2nc_mean com3D2nc com3D2nc_mean com2D2nc stat2nc isccp2nc modis2nc misr2nc com3D2nc_sep 2Dbin2nc_sep 2Dcom_sep2one 2Dbin_sep2one com3D_sep2one bin3D_sep2one glue_movie_raw
#all: bin2D2nc bin3D2nc 2Dbin2nc 2Dbin2nc_mean bin3D2nc_mean stat2nc isccp2nc modis2nc misr2nc com3D2nc_sep 2Dbin2nc_sep 2Dcom_sep2one 2Dbin_sep2one com3D_sep2one bin3D_sep2one glue_movie_raw

.f:   
	$(FF) -o $@ -I./SRC $< ./SRC/hbuf_lib.f ./SRC/cape.f ./SRC/cin.f $(LDFLAGS) 
	#$(FF) $@ -I./SRC $< ./SRC/hbuf_lib.f ./SRC/cape.f ./SRC/cin.f $(LDFLAGS) 

clean: 
	rm bin* com* stat* 2* isccp* modis* misr* *.o glue*

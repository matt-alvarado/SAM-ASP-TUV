c     gfortran only (UTIL/build_util.sh with SAMASP_FC=gfortran): tools that
c     declare "external iargc" need a real routine; gfortran's iargc is intrinsic.
      integer function iargc()
      iargc = command_argument_count()
      end

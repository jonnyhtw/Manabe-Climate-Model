      subroutine check_err(iret)

      integer iret

#include 'netcdf.inc'

      if (iret .ne. NF_NOERR) then
ccc     print *, nc_strerror(iret)
        print *, iret, NC_NOERR
        call aquit (20)
        stop
      endif
      
      return
      end


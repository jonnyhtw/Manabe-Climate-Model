      subroutine check_err(iret)
      integer iret
      include 'netcdf.inc'
      if (iret .ne. NF_NOERR) then
      print *, nf_strerror(iret)
      call aquit (20)
      stop
      endif
      end


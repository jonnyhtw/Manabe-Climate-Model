      subroutine odam (lu, nwrs, nfrst, a)
ccc   implicit real*8 (a-h,o-z)

      include 'netcdf.inc'

#include "ncstuf.h"
      integer var_id

#ifdef multitasking
cfpp$ noconcur r
#endif
#ifdef diskless
c
c=======================================================================
c
c  odam (ocean direct access manager) is a set of routines which
c       simulates memory/disc data transfer management, using the
c       central memory array, "big" as a virtual disc. this results
c       in a totally memory contained system. the arguments, in terms
c       of normal disc i/o, are:
c         lu   =i/o unit number
c         ntot =length of unit, in words
c         nblk =length of each block on the unit
c         nbuf =number of buffers supplied to the unit
c         nwrs =number of words to transfer
c         nfrst=unit address of the first word to be transferred
c         a    =origination/destination array in memory
c
c=======================================================================
c
#include "param.h"
#include "iounit.h"
c
c---------------------------------------------------------------------
c     "sdisk" is the simulated disk using memory
c---------------------------------------------------------------------
c
      parameter (n11=20, n12=imt*jmt*nkflds
     $,           nsl=nslab*jmt, ntb=n11+n12+ntlev*nsl)
# ifdef t90_1st
      parameter (nta=n11+(imt*jmt*5))
      dimension adisk(nta), kkmt(imt,jmt)
# endif
      dimension sdisk(ntb)
      dimension a(nwrs)
      save sdisk
c
c-----------------------------------------------------------------------
c     begin executable code
c-----------------------------------------------------------------------
c
      entry ostart (lu, ntot, nblk, nbuf)
      return
      entry ofind (lu, nwrs, nfrst)
      return
c
c-----------------------------------------------------------------------
c     virtual disc-to-memory transfer. used to move contiguous portions
c     of the virtual disk to contiguous portions of memory
c-----------------------------------------------------------------------
c
      entry oget (lu, nwrs, nfrst, a)
      nskp = 0
      if (lu .gt. kontrl) nskp = nskp + n11
      if (lu .gt. kflds) nskp = nskp + n12
      ns = nfrst + nskp
      if (lu .ge. labs(1)) ns = nskp + (ntlev*((nfrst-1)/nwrs)+
     $                                 (lu-labs(1)))*nwrs +1
      ne = ns + nwrs - 1
      do 100 n=ns,ne
        a(n-ns+1) = sdisk(n)
100   continue
      return
c
c---------------------------------------------------------------------
c     memory-to-virtual disc transfer. used to move contiguous portions
c     of memory to contiguous portions of the virtual disc
c---------------------------------------------------------------------
c
      entry oput (lu, nwrs, nfrst, a)
      nskp = 0
      if (lu .gt. kontrl) nskp = nskp + n11
      if (lu .gt. kflds) nskp = nskp + n12
      ns =nfrst + nskp
      if (lu .ge. labs(1)) ns = nskp+(ntlev*((nfrst-1)/nwrs)+
     $                          (lu-labs(1)))*nwrs+1
      ne = ns + nwrs - 1
      do 200 n=ns,ne
        sdisk(n) = a(n-ns+1)
200   continue
      return
c
      entry oclose (lu)
      return
c
c---------------------------------------------------------------------
c     initialize virtual disc from tape
c---------------------------------------------------------------------
c
      entry ord (lu)
      if (lu .ne. 21) then
        read (lu) sdisk
      else
cread in netcdf data
        iret = nf_open('ocean_r30_in.nc',nf_read, ncid)
        call check_err(iret)
        iret = nf_inq_varid(ncid, "sdisk", var_id)
        call check_err(iret)
        iret = nf_get_var_double (ncid, var_id, sdisk)
        call check_err(iret)
      endif
cnetcdf

      return
c
# ifdef t90_1st
c---------------------------------------------------------------------
c     when starting from a C90 or YMP restart file, ordkmt allows
c     one to read in the integer kmt field
c---------------------------------------------------------------------
c
      entry ordkmt (lu,kkmt)
      rewind lu
      read (lu) adisk, kkmt
      rewind (lu)
      return
# endif
c
c---------------------------------------------------------------------
c     save virtual disc to tape
c---------------------------------------------------------------------
c
      entry owrt (lu)
      if (lu .ne. 21) then
        rewind lu
        write (lu) sdisk
      else
cnetcdf
        iret = nf_put_var_double(ncid, sdsk_id, sdisk)
        call check_err(iret)
        iret = nf_close(ncid)
        call check_err(iret)
cnetcdf
      endif
      return
#define ioalreadyset
#endif
#if defined crayio && !defined ioalreadyset
#include "param.h"
c
c=======================================================================
c
c     odam (ocean direct access manager) is a set of cray routines for
c     doing disk i/o. The arguments are:
c
c     lu   = i/o unit number
c     ntot = length of unit, in words
c     nblk = length of each block on the unit
c     nbuf = number of buffers supplied to the unit
c     nwrs = number of words to transfer
c     nfrst= unit address of the first word to be transferred
c     a    = origination/destination array in memory
c
c=======================================================================
c
      dimension a(nwrs)
c
      entry ostart (lu, ntot, nblk, nbuf)
      istats=1
      call wopen (lu, nblk/512+1, istats)
      return
      entry ofind (lu, nwrs, nfrst)
cmic$ guard
      call seek (lu, nfrst, nwrs)
cmic$ end guard
      return
      entry oget (lu, nwrs, nfrst, a)
cmic$ guard
      call getwa (lu, a, nfrst, nwrs)
cmic$ end guard
      return
      entry oput (lu, nwrs, nfrst, a)
cmic$ guard
#ifdef multitasking
      call aputwa (lu, a, nfrst, nwrs)
#else
      call putwa (lu, a, nfrst, nwrs)
#endif
cmic$ end guard
      return
      entry oclose (lu)
      call wclose (lu)
      return
#define ioalreadyset
#endif
#if defined fio && !defined ioalreadyset
c
c=======================================================================
c
c     odam (ocean direct access manager) is a set of fio routines for
c     doing disk i/o. The arguments are:
c
c     lu   = i/o unit number
c     ntot = length of unit, in words
c     nblk = length of each block on the unit
c     nbuf = number of buffers supplied to the unit
c     nwrs = number of words to transfer
c     nfrst= unit address of the first word to be transferred
c     a    = origination/destination array in memory
c
c=======================================================================
c
      dimension a(nwrs)
c
      entry ostart (lu, ntot, nblk, nbuf)
      open (unit=lu, access='DIRECT', form='UNFORMATTED',
     $           recl=8*nblk)
      return
      entry ofind (lu, nwrs, nfrst)
      return
      entry oget (lu, nwrs, nfrst, a)
      read (lu, rec=(nfrst-1)/nwrs+1) a
      return
      entry oput (lu, nwrs, nfrst, a)
      write (lu, rec=(nfrst-1)/nwrs+1) a
      return
      entry oclose (lu)
      close (lu)
      return
#define ioalreadyset
#endif
#ifndef ioalreadyset
      dimension a(nwrs)
c
      entry ostart (lu, ntot, nblk, nbuf)
      open (unit=lu, access='DIRECT', form='UNFORMATTED',
     $           recl=8*nblk)
      return
      entry ofind (lu, nwrs, nfrst)
      return
      entry oget (lu, nwrs, nfrst, a)
      read (lu, rec=(nfrst-1)/nwrs+1) a
      return
      entry oput (lu, nwrs, nfrst, a)
      write (lu, rec=(nfrst-1)/nwrs+1) a
      return
      entry oclose (lu)
      close (lu)
      return
#endif
#undef ioalreadyset
      end
      subroutine xfer (buffer, slab)
ccc   implicit real*8 (a-h,o-z)
#ifdef multitasking
cfpp$ noconcur r
#endif
c
c     transfer contiguous data from the buffer area to non-contiguous
c     locations (variables) in the memory slab window
c
#include "param.h"
      dimension slab(imt,km,nslabs,ntau,nvar), buffer(imt,km,nvar)
c
      do 100 n=1,nvar
        do 90 k=1,km
          do 80 i=1,imt
            slab(i,k,1,1,n) = buffer(i,k,n)
80        continue
90      continue
100   continue
c
      return
      end
c23456789a123456789b123456789c123456789d123456789e123456789f123456789g12
# ifdef f90type
c
c=======================================================================
c
      subroutine odami (lu, nwrs, nfrst, inta)
ccc   implicit real*8 (a-h,o-z)
c
c=======================================================================
c
c   The f90type ifdef option is designed to adapt parts of the MOM 1
c   odam i/o routines to the Fortran 90 requirement of having the
c   data types of all arguments in calling statements match the dummy
c   arguments of the corresponding subroutine or entry statements.
c   All data actually written to disk or virtual disk are of type
c   real when this option in used.  Conversions to or from integer
c   type data are performed as needed by the code associated with the
c   new ogeti and oputi entry points.
c   Designed to work with all 3 i/o options (diskless, crayio and fio)
c   (Follows interpretation of CF90 Fortran Language Reference Manual,
c   Vol 1, SR-3902 2.0,  section 2.5 - Argument Association)
c                                                     (KD - Sept 1996)
c
c=======================================================================
c
cSGI
ccc   integer*8 inta
      integer*4 inta
cSGI
      dimension    a(nwrs)
      dimension inta(nwrs)
c
c---------------------------------------------------------------------
c
      entry ogeti (lu, nwrs, nfrst, inta)
c
c  introduce ogeti entry point to allow conversion from
c  real to integer after retrieving real data from disk
c
      call oget (lu, nwrs, nfrst, a)
c
      do n=1,nwrs
        inta(n) = nint(a(n))
      end do
c
      return
c
c---------------------------------------------------------------------
c
      entry oputi (lu, nwrs, nfrst, inta)
c
c  introduce oputi entry point to allow integer type data
c  to be loaded into real arrays before being put to disk
c
      do n=1,nwrs
        a(n) = float(inta(n))
      end do
c
      call oput (lu, nwrs, nfrst, a)
c
c---------------------------------------------------------------------
      return
c
      end
# endif

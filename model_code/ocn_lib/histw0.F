      subroutine histw0
      implicit real*8 (a-h,o-z)
c
#if defined worldocn && defined selfdoc
c     NOTE: the user should consult with keith before modifying the
c           version number - the convention will be that the version
c           number is yyyymmdd (year month day)
c
c     NOTE: this code is specifically designed for GFDL global ocean
c           climate studies and produces self documenting output for
c           use with the smmaps / latxht family of graphical analysis
c           programs    ( Keith Dixon - 1993 )
#ifdef cyclic
c     NOTE: this code assumes cyclic boundary conditions
c           a few minor changes would be needed to i loops to remove
c           this restriction
c
c=======================================================================
c
c     "histw0" writes out the self documenting header information
c     of history files.
c
c=======================================================================
c
#include "param.h"
#include "accel.h"
#include "coord.h"
#include "ctmngr.h"
#include "docnam.h"
#include "grdvar.h"
#include "histjr.h"
#include "iounit.h"
#include "scalar.h"
c
      logical season
c
# ifdef seasonal
      data season /.true./
# else
      data season /.false./
# endif
c
c
c     variables:
c       idvers = selfdoc version date, the convention will be that 
c           the version represents the day of modification, the form of
c           the integer is yyyymmdd (year month day)
c
c   19930930:  the original release of selfdoc
c   19931104:  the release that switch from pointwise to slab-wise
c              history output (the header info did not change)
ccccc data idvers / 19930930 / 
      data idvers / 19931104 /  
c
c=======================================================================
c
c---------------------------------------------------------------------
c     write self documenting info to "iohist"
c---------------------------------------------------------------------
c
c---------------------------------------------------------------------
c     1st, write the version number date
c---------------------------------------------------------------------

      write(iohist) idvers
      print *,'> selfdoc version date number: ',idvers
c
c---------------------------------------------------------------------
c     2nd, write the experiment name
c---------------------------------------------------------------------
c
      write(iohist) expnam
c
c---------------------------------------------------------------------
c     3rd, write in info concerning the dimensions of the grid
c     the number of tracers carried in the model, and the amount of
c     data to be written out for each i,j point
c---------------------------------------------------------------------
c
      idim = imt - 2
      jdim = jmt - 2
      kdim = km
      ntr = nt
      nvtot = iblk
c
      write(iohist) idim, jdim, kdim, ntr, nvtot
c
c---------------------------------------------------------------------
c     4th, write info on x,y,z grid locations (from these values
c     all grid locations can be easy computed
c---------------------------------------------------------------------
c
      write(iohist) (xt(its),its=2,imt-1),
     $              (xu(iuv),iuv=2,imt-1),
     $              (yt(jts),jts=2,jmt-1),
     $              (yu(juv),juv=2,jmt-1), dzt
c
c---------------------------------------------------------------------
c     5th, write out info on time levels and the timesteps and
c     acceleration techniques used (4 writes)
c---------------------------------------------------------------------
c
#ifdef gtime
      totday = ctotdays
      tyear = icyear + year0
      tday = icdayofmonth
#endif
      write(iohist) totday, tday, tyear, relyr, itt
      write(iohist) stampb, stampe, cmonth
      write(iohist) knthst
      write(iohist) dtts, dtuv, dtsf, dtxcel
c
ctestself
cSGI debug
      arelyr = relyr 
cc    print *,'> totday: ',totday
cc    print *,'> tday  : ',tday
cc    print *,'> tyear : ',tyear
cc    print *,'> arelyr : ',arelyr
cc    print *,'> itt   : ',itt
      print *,'> totday,tday,tyear,arelyr,itt: ',totday, tday,
     $           tyear, arelyr, itt
cSGI debug
ccc   print *,'> stampb stampe knthst: ',stampb, stampe, knthst
      print *,'> stampb stampe cmonth: ',stampb, stampe, cmonth
      print *,'> knthst: ',knthst
      print *,'> dtts dtuv dtsf: ',dtts, dtuv, dtsf
      print *,'> dxcel: ',dtxcel
ctestself
c
c---------------------------------------------------------------------
c     6th, write logical for (seasonal/annual mean) forcing (T/F)
c---------------------------------------------------------------------
c
      write(iohist) season
      print *, '> season: ',season
c
c---------------------------------------------------------------------
c     7th, write variable names and units and grid type info
c     (cheat_sheet) to "iohist"
c---------------------------------------------------------------------
c
      write(iohist) varnam, varunt, varhgr, varvgr, klevar
c
ctestself
ccc   do 199 nv=1,iblk
ccc     print 9199, nv, varnam(nv), varunt(nv), varhgr(nv),
ccc  $              varvgr(nv), klevar(nv)
199   continue
ctestself
#endif
c
#endif
      return
c
9199  format(' >',1x,i4,2x,a16,2x,a16,2x,a2,2x,a2,2x,i3)
c
      end

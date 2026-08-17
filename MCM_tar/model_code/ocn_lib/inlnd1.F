      subroutine inlnd1
      implicit real*8 (a-h,o-z)
#ifdef inland
c
c     this subroutine handles the initial set up for inland
c     sea communications.
c
c     using data provided in blkdta, checks are performed to ensure that
c     the locations and rates of inland mixing requested are valid, and
c     the requested mixing rates of cm**3/sec are converted into the
c     fraction of the model grid boxes to be mixed per second, and
c     checked.
c
#include "param.h"
#include "accel.h"
#include "cinlnd.h"
#include "grdvar.h"
#include "levind.h"
#include "scalar.h"
#if defined cpldleg3 || defined cpldleg2
#include "switch.h"
#endif
c
c     local variables:
c
c     depth   = depth thickness over which inland mixing is done [cm]
c     colvol  = column volume over which inland mixing is done [cm**3]
c     errmix  = logical flag for checking validity of mixing rates
c               level by level
c
      logical errmix
c
c
c=======================================================================
c
c     for each pair of i,j points to be mixed in this manner
c     the following is done:
c
c=======================================================================
c 
      do 100 nin=1,ninlnd
c
c       check for invalid i,j grid locations
c
        if (jinlnd(nin)   .lt. 1 .or. jinlnd(nin)   .gt. jmt .or.
     $      iinlnd(nin,1) .lt. 2 .or. iinlnd(nin,2) .gt. imtm2) then
          write(stdout,996) nin, jinlnd(nin), (iinlnd(nin,ip),ip=1,2)
          write(stderr,996) nin, jinlnd(nin), (iinlnd(nin,ip),ip=1,2)
          stop 'inlnd1'
        endif
        jj = jinlnd(nin)
c
c       check for i positions that are not in ascending order
c
         if (iinlnd(nin,1) .gt. iinlnd(nin,2)) then
             write (stdout,997) nin, iinlnd(nin,1), iinlnd(nin,2)
             write (stderr,997) nin, iinlnd(nin,1), iinlnd(nin,2)
             stop 'inlnd2'
         endif
c
c       check for attempts to mix land rather than sea
c
         ktop = kinlnd(nin,1)
         kbot = kinlnd(nin,2) 
666   format(/'nin, iinlnd(nin,1), jj, kmt(iinlnd(nin,1),jj)'/4i5)
         if (ktop.gt.kbot .or. kbot.gt.kmt(iinlnd(nin,1),jj) .or. 
     $       kbot.gt.kmt(iinlnd(nin,2),jj) .or. ktop.lt.1) then 
             write (stdout,998) nin, ktop, kbot, 
     $         kmt(iinlnd(nin,1),jj), kmt(iinlnd(nin,2),jj) 
             write (stderr,998) nin, ktop, kbot, 
     $         kmt(iinlnd(nin,1),jj), kmt(iinlnd(nin,2),jj) 
             stop 'inlnd3' 
         endif 
c
c        calculate depth range over which inland mixing will be done
c
          depth = c0 
          do 80 k=ktop,kbot
            depth = depth + dzt(k) 
80        continue 
c
c         for each of the two points to be mixed, convert from volume to
c         be mixed per second "vinlnd(nin)" to model box fraction to be
c         mixed per second "binlnd(nin,l)"      (see "cinlnd.h")
c
          do 90 l=1,2 
            boxar  = cst(jj)*dxt(iinlnd(nin,l))*dyt(jj)
            colvol = boxar*depth
            binlnd(nin,l) = vinlnd(nin) / colvol
c
c           check for attempts to mix at a rate so fast that it can not
c           be achieved in a single, leapfrog timestep
c           (for the two boxes to be thoroughly mixed, no more than
c           one-half of the volume of a given grid box can be
c           in effect transported into the other gridbox)
c
            errmix = .false.
            do 85 k=ktop,kbot
              factor = binlnd(nin,l)*c2*dtts*dtxcel(k)
              if (factor .gt. p5) then
                write (stdout,999) nin, factor, iinlnd(nin,l), 
     $                jinlnd(nin), k ,binlnd(nin,l), vinlnd(nin), 
     $                ktop, kbot, colvol
                write (stderr,999) nin, factor, iinlnd(nin,l), 
     $                jinlnd(nin), k ,binlnd(nin,l), vinlnd(nin), 
     $                ktop, kbot, colvol
                stop 'inlnd4'               
              endif
85          continue
c 
90        continue
c
#if defined cpldleg3 || defined cpldleg2
          if (pintro) then
#endif
          write (stdout,9101) nin, jinlnd(nin), ktop, 
     $          kbot, (iinlnd(nin,ll),ll=1,2), vinlnd(nin), 
     $          (binlnd(nin,ll), iinlnd(nin,ll),ll=1,2) 
#if defined cpldleg3 || defined cpldleg2
          endif
#endif
c
100   continue
c
#endif
      return
c
996   format(/' Error => problem with inland sea mixing #',i3,
     $': out of bounds i,j grid location requested'/' asked for j-row ',
     $i5,' i points ',i5,' and',i5,/)
997   format(/' Error => problem with inland sea mixing #',i3,
     $': improper i-values requested'/' iinlnd(nin,1)=',i5,' and is not'
     $' less than iinlnd(nin,2)=',i5,/)
998   format(/' Error => problem with inland sea mixing #',i3,
     $': improper k-levels requested'/' ktop, kbot =',2i5,', kmt at the'
     $,' pair of points to be mixed =',2i5,/)
999   format(/' Error => problem with inland sea mixing (#',i3,')'/
     $' attempting to mix at to fast a rate: asking for more than half'
     $' the volume of a'/' grid box (',e12.6,') to be transported into'
     $,' its neighbor on leapfrog timesteps'/' at i, j, k: ',3i5,/,
     $' - as specified now, fraction of box mixed per sec = ',e12.6,/,
     $' potential fixes: shorten dtts*dtxcel(k), reduce mixing rate ',
     $'vinlnd (now ',e12.6,' cm**3/sec)'/
     $,' or mix over a larger depth range (now from k=',i3,' to ',i3
     $,' => water column volume = ',e12.6,' cm**3)'/)
9101  format(' for inland sea communication pair number',i2,/
     $' j row=',i3,' from level',i3,' to',i3,' at i=',i4,',',i4/
     $' simulated flow in = flow out = ',e12.6,' cm**3/sec'/' so mix ',
     $e12.6,' fraction of i=',i4,' with ',e12.6,' of i=',i4,' per sec'/)
c
      end


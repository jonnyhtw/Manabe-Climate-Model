      subroutine icedyn (j)
      implicit real*8 (a-h,o-z)
c
#ifdef seaice
#include "param.h"
#include "ctmngr.h"
#include "cvbc.h"
#include "grdvar.h"
#include "scalar.h"
#include "seasnl.h"
#include "seayce.h"
#include "slabs.h"
#include "switch.h"
c
c   calculate movement and diffusion of sea ice 
c
c========================================================================
c
c  variables:
c         for the next 4 variables, the n, s, and w suffixes represent
c         the northern, southern, and western edges of the gridboxes.
c    ycjam[n,s,w] = sea ice jam masks; set to zero if conditions are
c                   to prevent sea ice advection, otherwise set to one.
c    ycthk[n,s,w] = mean sea ice thickness at gridbox edges (average of
c                   the central gridbox and an adjacent gridbox).
c    ycdel[n,s,w] = difference between the sea ice thickness at an
c                   adjacent gridbox and the central gridbox (positive
c                   if thicker sea ice is to the north or east).
c    ycspd[n,s,w] = average speed of advective surface currents flowing
c                   normal to gridbox edges (computed at the gridbox
c                   edges; positive if northward or eastward).
c
c    snowfx       = rate of snowfall over sea ice gridbox [cm/sec] to
c                   be added to sea ice thickness.
c
c
      dimension ycjamn(imt), ycjams(imt), ycjamw(imt)
      dimension ycthkn(imt), ycthks(imt), ycthkw(imt)
      dimension ycdeln(imt), ycdels(imt), ycdelw(imt)
      dimension ycspdn(imt), ycspds(imt), ycspdw(imt)
c
c--------------------------------------------------------------------- 
c     construct physical operators using statement functions 
c--------------------------------------------------------------------- 
c 
c     the zonal advection of sea ice
c
      uycx(i) = ( ycjamw(i  )*ycthkw(i  )*ycspdw(i  )
     $          - ycjamw(i+1)*ycthkw(i+1)*ycspdw(i+1) )
     $          *cstr(j)*dxtr(i)
c
c     the meridional advection of sea ice
c
      vycy(i) = ( ycjams(i)*ycthks(i)*ycspds(i)*csu(j-1)
     $          - ycjamn(i)*ycthkn(i)*ycspdn(i)*csu(j  ) )
     $          *dytr(j)*cstr(j)
c
c     the zonal diffusion of sea ice
c
      ycxx(i) = ahyc*cstr(j)*cstr(j)*
     $        ( fm(i+1,1,jc)*ycdelw(i+1)*dxur(i+1)*dxur(i+1) -
     $          fm(i-1,1,jc)*ycdelw(i  )*dxur(i  )*dxur(i  ) )
c
c     the meridional diffusion of sea ice
c
      ycyy(i) = ahyc*dytr(j)*cstr(j)*
     $        ( fm(i,1,jp)*csu(j  )*dyur(j  )*ycdeln(i) -
     $          fm(i,1,jm)*csu(j-1)*dyur(j-1)*ycdels(i) )
c
c========================================================================
c
#ifdef cyclic
      u (1,1,jc,nc) = u (imtm1,1,jc,nc)
      u (1,1,jm,nc) = u (imtm1,1,jm,nc)
      v (1,1,jc,nc) = v (imtm1,1,jc,nc)
      v (1,1,jm,nc) = v (imtm1,1,jm,nc)
      fm(1,1,jc)    = fm(imtm1,1,jc)
      yc(1,j,nc)    = yc(imtm1,j,nc)
      dxu(1)        = dxu(imtm1)
      yc(imt,j+1,nc) = yc(2,j+1,nc)
      yc(imt,j  ,nc) = yc(2,j  ,nc)
      yc(imt,j-1,nc) = yc(2,j-1,nc)
      fm(imt,1,jc)  = fm(2,1,jc)
      dxur(imt)     = dxur(2)
#endif
#ifdef skipland
      do 1300 lstrip=1,lseg
        istart = ist(j,lstrip)
        iend   = iet(j,lstrip)
        if (istart .ne. 0) then
#endif
c
c     prepare to compute effect of horizontal advection and horizontal
c     diffusion on the rate of change of sea ice thicknesses
c
      do 100 i=istart,iend+1
        ycthkn(i) = p5*(yc(i,j+1,nc) + yc(i,  j  ,nc))
        ycthks(i) = p5*(yc(i,j  ,nc) + yc(i,  j-1,nc))
        ycthkw(i) = p5*(yc(i,j  ,nc) + yc(i-1,j  ,nc))
        ycdeln(i) =     yc(i,j+1,nc) - yc(i  ,j  ,nc)
        ycdels(i) =     yc(i,j  ,nc) - yc(i  ,j-1,nc)
        ycdelw(i) =     yc(i,j  ,nc) - yc(i-1,j  ,nc)
#ifdef newice
c       calculate simple mean of current speeds at gridbox edges
        ycspdn(i) = p5*(v(i-1,1,jc,nc) + v(i  ,1,jc,nc))
        ycspds(i) = p5*(v(i-1,1,jm,nc) + v(i  ,1,jm,nc))
        ycspdw(i) = p5*(u(i-1,1,jm,nc) + u(i-1,1,jc,nc))
#else
c       calculate weighted mean of current speeds at gridbox edges
        ycspdn(i) = p5*dxtr(i)*
     $              (v(i-1,1,jc,nc)*dxu(i-1) + v(i,1,jc,nc)*dxu(i))
        ycspds(i) = p5*dxtr(i)*
     $              (v(i-1,1,jm,nc)*dxu(i-1) + v(i,1,jm,nc)*dxu(i))
        ycspdw(i) = p5*dytr(j)*
     $              (u(i-1,1,jm,nc)*dyu(j-1) + u(i-1,1,jc,nc)*dyu(j))
#endif
100   continue
#ifdef newice
c     if the sea ice at a box edge is thicker than "ycejam", then
c     do not allow advection to add ice to that gridbox by setting
c     the ice jam mask for that gridbox edge to zero; otherwise set
c     the mask to one.
#else
c     if the sea ice at a box edge is thicker than "ycejam", then
c     do not allow advection to remove ice from that box by setting
c     the ice jam mask for that gridbox edge to zero; otherwise set
c     the mask to one.
#endif
c
      do 200 i=istart,iend
#ifdef newice
        if (ycthkn(i).gt.ycejam .and. (ycspdn(i)*ycdeln(i)).lt.c0) then
          ycjamn(i) = c0
        else
          ycjamn(i) = c1
        endif
        if (ycthks(i).gt.ycejam .and. (ycspds(i)*ycdels(i)).lt.c0) then
          ycjams(i) = c0
        else
          ycjams(i) = c1
        endif
        if (ycthkw(i).gt.ycejam .and. (ycspdw(i)*ycdelw(i)).lt.c0) then
          ycjamw(i) = c0
        else
          ycjamw(i) = c1
        endif
#else
        if (ycthkn(i).gt.ycejam .and. (ycspdn(i)*ycdeln(i)).gt.c0) then
          ycjamn(i) = c0
        else
          ycjamn(i) = c1
        endif
        if (ycthks(i).gt.ycejam .and. (ycspds(i)*ycdels(i)).gt.c0) then
          ycjams(i) = c0
        else
          ycjams(i) = c1
        endif
        if (ycthkw(i).gt.ycejam .and. (ycspdw(i)*ycdelw(i)).gt.c0) then
          ycjamw(i) = c0
        else
          ycjamw(i) = c1
        endif
#endif
200   continue
c  
#ifdef cyclic
      ycjamw(imt) = ycjamw(2)
#endif
c  
      do 999 i=istart,iend
        yca(i,j) = yc(i,j,nm) + c2dtts*fm(i,1,jc)*
     $           (uycx(i)+vycy(i)+ycxx(i)+ycyy(i)+snowfx(i))
999   continue
c  
#ifdef ycetest
cc   Artifically set yca > 400cm to test jammimg:
cc    if (j.gt.38) then
cc      do 9030 i=10,20
cc        yca(i,j) = ycejam + float(10*i) + float(j)
9030    continue
cc    endif
      if ( j.eq.74 ) then
        do 9990 i=11,15
          write(6,9190)
9190      format(1x,'in subr ice: j, i, yca, ycnc, ycnm, ',
     $              'fm, uycx, vycy, ycxx, ycyy, snowfx=')
          if ( fm(i,1,jc).ne.c0 ) then
            facw = cstr(j) * 2.0
            facn = cstr(j) * dytr(j) * csu(j) * 0.5
            facs = cstr(j) * dytr(j) * csu(j-1) * 0.5
            write(stdout,9191) j,i,
     $       yca(i,j), yc(i,j,nc), yc(i,j,nm), 
     $        fm(i,1,jc), uycx(i), vycy(i),
     $      ycxx(i), ycyy(i), snowfx(i)
cc   $,  2.*ycthkw(i), ycspdw(i)*facw, ycdelw(i), ycjamw(i),
cc   $   2.*ycthkn(i), ycspdn(i)*facn, ycdeln(i), ycjamn(i),
cc   $   2.*ycthks(i), ycspds(i)*facs, ycdels(i), ycjams(i)
cc 9191        format(1x,2i5,9e13.5,/,(11x,4e13.5))
9191        format(1x,2i5,9e13.5)
          endif
9990    continue
      endif
#endif
c
#ifdef skipland
1300  continue
c
#endif
      return
      end
c
      subroutine icethm (j, ycconv)
      implicit real*8 (a-h,o-z)
c
#include "param.h"
#include "accel.h"
#include "ctask.h"
#include "ctmngr.h"
#include "cvbc.h"
#include "grdvar.h"
#include "levind.h"
cms
#include "seasnl.h"
cms
#include "scalar.h"
#include "seayce.h"
#include "slabs.h"
#include "switch.h"
c
      logical ycconv
c
      dimension aheat(imt), factor(imt), ll(imt), zdep(imt)
      dimension tastor(imt,2), ycstor(imt)
      dimension tastok(imt,2)
c
      data sfact /0.035/
c
#ifdef skipland
      do 1300 lstrip=1,lseg
        istart = ist(j,lstrip)
        iend   = iet(j,lstrip)
        if (istart .ne. 0) then
#endif
c
c   compute depth of mixed layer
c
c    "ycconv" = indicator of whether or not convective adjustment has
c           been performed in tracer yet this timestep 
c           (false=no, true=yes)
c    "ll" = number of levels with temperatures identical to level 1
c           (the number of model levels that are mixed)
c    "zdep" = mixed layer depth [cm] (dtxcel factor takes into effect
c           the diminished heat capacity implied by accelerating
c           dtts with increasing depth)
c
        if (.not. ycconv) then
c    before convective adjustment, set mixed layer depth "zdep" to
c    dzt(1)
          do 100 i=istart,iend
            ll(i) = 1
            zdep(i) = dzt(1)
100       continue
c    before convective adjustment, initialize arrays for history data
          if (histts) then
            if (.not. mxpas2) then
              do 150 i=istart,iend
                dtyce(i,j,1) = c0
                dtyce(i,j,2) = c0
                stfyc(i,j,1) = c0
                stfyc(i,j,2) = c0
150           continue
              ycknt = c0
            endif
          endif
c
        else
c    after convective adjustment, attempt to calculate mixed layer 
c    depth "zdep" by determining extent of isothermal conditions
c    extending from the surface to lower levels
          do 200 i=istart,iend
            ll(i) = 1
            zdep(i) = c0
            kz = kmt(i,j)
            do 170  k=1,kz
              if (ta(i,1,1) .ne. ta(i,k+1,1)) go to 180
              ll(i) = ll(i) + 1
170         continue
180         continue
            llmax = ll(i)
            do 190 lp=1,llmax
              zdep(i) = zdep(i) + (dzt(lp) / dtxcel(lp))
190         continue
200       continue
        endif
c
c    store T & S before ice thermodynamics computation for use 
c     in computing changes caused by ice thickness changes.
        do 210 i=istart,iend
          tastor(i,1) = ta(i,1,1)
          tastor(i,2) = ta(i,1,2)
          ycstor(i  ) = yca(i,j)
210     continue
c
c     set ocean temp to freezing and convert heat lost/gained
c     into ice. if ice<0.0, freeze negative ice and add heat to ocean.
c
c   "factor" is a dirty trick, inserted by computer nerds seeking
c   vectorization at the expense of code understandability.
c   They have been sacked.  sorry for the inconvenience.  (kd)
c
        althtr = c1 / altht
c
c   this code will produce negative ice thicknesses over any 
c   gridbox having a level 1 temperature above freezing, and 
c   where yca(i,j) is zero
c   "factor" stores the magnitude of this silliness, to
c   allow things to be corrected below
c
        do 300 i=istart,iend
          factor(i) = (freeze - ta(i,1,1)) * althtr
          yca(i,j) = yca(i,j) + (factor(i) * zdep(i))
300     continue
c
c   freeze the entire surface of the world ocean & make it realllllll
c   fresh so that we can melt it and salt it later is a nicely 
c   vectorized manner (what a rush!)  aheat, like factor, helps
c   recover from these sins later
c
        do 400 i=istart,iend
          ta    (i,1,1) = freeze
          ta    (i,1,2) = ta(i,1,2) + (sfact * factor(i))
          aheat (i    ) = yca(i,j)/ zdep(i)
400   continue
c
c     for negative sea ice (most of the world at this point) adjust
c     surface temperature (warming effect) and salinity (more saline)
c
      do 500 i=istart,iend
        if (yca(i,j) .lt. c0) then
          ta(i,1,1) = ta(i,1,1) - (altht * aheat(i))
          ta(i,1,2) = ta(i,1,2) - (sfact * aheat(i))
        endif
500   continue
c
c     after convection, set temperatures to freezing throughout mixed
c     layer "zdep" anywhere unfortunate enough to have had convection
c     extending down from the surface
c
      if (ycconv) then
        do 600 i=istart,iend
          if (ll(i) .ge. 2) then
            do 590 lp = 2, ll(i)
              tastok(i,1) = ta(i,lp,1)
              tastok(i,2) = ta(i,lp,2)
              ta(i,lp,1) = freeze
              ta(i,lp,2) = ta(i,lp,2) + (sfact * factor(i))
c     (rush #2) for negative sea ice (remember your last rush?), set 
c     ice thickness to zero, and adjust surface temperature (warming 
c     effect) and salinity (more saline) to undo your sins
              if (yca(i,j) .lt. c0) then
                ta(i,lp,1) = ta(i,lp,1) - (aheat(i) * altht)
                ta(i,lp,2) = ta(i,lp,2) - (aheat(i) * sfact)
              endif
              stfyc(i,j,1) = stfyc(i,j,1) + (ta(i,lp,1) - tastok(i,1)) /
     $                       (dzt(lp) * c2dtts*dtxcel(lp))
              stfyc(i,j,2) = stfyc(i,j,2) + (ta(i,lp,2) - tastok(i,2)) /
     $                       (dzt(lp) * c2dtts*dtxcel(lp))
590         continue
          endif
600     continue
c
      endif
c
c     one last rush ... time to come down ... set all those
c     silly negative sea ice values back to zero, and all is well 
c     again in the nicely vectorized world of the computer nerds!
c
      do 700 i=istart,iend
        yca(i,j) = max(yca(i,j),c0)
700   continue
c
      if (histts) then
        if (j .eq. jerow) ycknt = ycknt + c1
        do 710 i=istart,iend
          if (yca(i,j) .ne. ycstor(i)) then
            stfyc(i,j,1) = stfyc(i,j,1) + (ta(i,1,1) - tastor(i,1)) /
     $                                    (dzt(1) * c2dtts*dtxcel(1))
            stfyc(i,j,2) = stfyc(i,j,2) + (ta(i,1,2) - tastor(i,2)) /
     $                                    (dzt(1) * c2dtts*dtxcel(1))
          endif
#if defined cpldleg2 || defined cpld_restore
          dtyce(i,j,1) = dtyce(i,j,1) + adjyc(i) * altht
          dtyce(i,j,2) = dtyce(i,j,2) + adjyc(i) * sfact
#endif
710     continue
      endif
c
#if defined cpldleg2 || defined cpld_restore
c-----------------------------------------------------------------------
c  compute serr and serr, reseting yca toward obs
c  store FLUXES of terr (cal/cm2/sec) and serr (gm/cm2/sec); see setvbc
c-----------------------------------------------------------------------
      if (.not. ycconv) then
        if (.not. mxpas2) then
        do 800 i=istart,iend
          terr(i,j) = terr(i,j) + adjyc(i) * altht
          serr(i,j) = serr(i,j) + adjyc(i) * sfact
c   Store terr and serr for history tape, current timestep only
c    (but accumulated with terrj, serrj from setvbc.F)
          terrj(i)  = terrj(i)  + adjyc(i) * altht
          serrj(i)  = serrj(i)  + adjyc(i) * sfact
          yca (i,j) = yca (i,j) - adjyc(i) * c2dtts*fm(i,1,jc)
800     continue
#ifdef testterr
c   also print terr, serr in setvbc.F before adding ice effects.
        if (j .eq. 12) then
          write(stdout,9112) itt, icntts,
     $      terr(12,12), serr(12,12), terrj(12), serrj(12), adjyc(12)
9112      format(1x,'in ice:    at (12,12), itt,icntts,terr,serr=',
     $                        2i10,5e13.5)
        endif
#endif
        endif
      endif
#endif
c
#ifdef ycetest2
      if (j.eq.74) then
        write(6,9192)
9192    format(1x,'in subr ice: j, i, yca(i,j), snowfx, ta(..1), ',
     $            'ta(..2), yceobs(..mahd), yceobs(..mbhd), adjyc, fm=')
        do 1110 i=11,15
          if(fm(i,1,jc) .ne. 0.)
     $      write(6,9193) j, i, yca(i,j), yc(i,j,nm), snowfx(i),
     $                    ta(i,1,1), ta(i,1,2), yceobs(i,j,mahd),
     $                    yceobs(i,j,mbhd), adjyc(i), fm(i,1,jc)
9193        format(1x,'j,i=',2i3,9e13.5)
1110    continue
      endif
#endif
#endif
      return
      end

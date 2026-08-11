      subroutine albsfc(ipt,seaalb)
c
c  This subroutine computes the surface albedos.
c  This calculation is based on temperature, snow depth, and
c  sea ice thickness.
c     
c       Called from: impt
c       Calls      : nothing
c
c  ipt is the i-loop index in impt (do 800 loop)
c
#include "parm.h"
c
#include "parmkp.h"
#include "parmjx.h"
#include "parmrad.h"
c
#include "albedo.h"
c
#include "cld1.h"
c
#include "ctlblk.h"
c
#include "grdfld.h"
c
#include "param.h"
c
#include "surfc.h"
c
#ifdef ice_dynamics
c
c variables of ice model used in climate model
#include "varice.h"
#endif
c
c
c compute critical temperature for use below
      tcrit = tfra - 10.0
c
c determine surface albedo for ocean grid points
c =================================================
c
      if (zin(ipt) .eq. 0.0) then
        awater = seaalb
c
c compute albedo of thick sea ice based on surface temperature
#ifdef ice_dynamics
c
c use ice surface temperature (not average) for ice points (cfi)
c
        if ( abs(tsice(ipt,jrow)) .le. tcrit) then
            aice = 0.80
          elseif ( abs(tsice(ipt,jrow)) .eq. tfra) then
            aice = 0.55
          else
            aice = ((tfra-abs(tsice(ipt,jrow)))*0.80 +
     $         (abs(ts(ipt,jrow))-tcrit)*0.55) / (tfra-tcrit)
        endif
#else
c
        if ( abs(ts(ipt,jrow)) .le. tcrit) then
            aice = 0.80
          elseif ( abs(ts(ipt,jrow)) .eq. tfra) then
            aice = 0.55
          else
            aice = ((tfra-abs(ts(ipt,jrow)))*0.80 +
     $         (abs(ts(ipt,jrow))-tcrit)*0.55) / (tfra-tcrit)
        endif
#endif
c
c assign final albedo based on sea ice thickness
        if (pkice(ipt) .eq. 0.0) then
            ca(nc) = awater
          elseif (pkice(ipt) .ge. 100.0) then
            ca(nc) = aice
          else
            ca(nc) = awater + sqrt(pkice(ipt)/100.0)*(aice-awater)
        endif
c
#ifdef ice_dynamics
c
c determine average albedo for grid point based on fractional
c ice coverage (cfi)
c
	if (pkice(ipt) .gt. 0.) then
        ca(nc) = ca(nc)*area(ipt,jrow,2) + awater*(1.-area(ipt,jrow,2))
        end if
#endif
      endif
c
c determine surface albedo for land grid points (including glaciers)
c =================================================
c
      if (zin(ipt) .eq. 1.0) then
c
c albedo determination for unglaciated points (albedo > 0.0)
        if (alb(ipt,jrow) .ge. 0.0) then
            abase = ralb(ipt)
c
c compute albedo for thick snow based on surface temperature
            if (abs(ts(ipt,jrow)) .le. tcrit) then
                asnow = 0.60
              elseif (abs(ts(ipt,jrow)) .eq. tfra) then
                asnow = 0.45
              else
                asnow = ((tfra-abs(ts(ipt,jrow)))*0.60 +
     $            (abs(ts(ipt,jrow))-tcrit)*0.45) / (tfra-tcrit)
            endif
c
c albedo determination for glaciated points (albedo < 0.0)
          elseif (alb(ipt,jrow) .lt. 0.0) then
c
c compute base albedo of glacier based on surface temperature
            if (abs(ts(ipt,jrow)) .le. tcrit) then
                abase = 0.65
              elseif (abs(ts(ipt,jrow)) .eq. tfra) then
                abase = 0.55
              else
                abase = ((tfra-abs(ts(ipt,jrow)))*0.65 +
     $            (abs(ts(ipt,jrow))-tcrit)*0.55) / (tfra-tcrit)
            endif
c
c compute albedo for thick snow based on surface temperature
            if (abs(ts(ipt,jrow)) .le. tcrit) then
                asnow = 0.80
              elseif (abs(ts(ipt,jrow)) .eq. tfra) then
                asnow = 0.65
              else
                asnow = ((tfra-abs(ts(ipt,jrow)))*0.80 +
     $            (abs(ts(ipt,jrow))-tcrit)*0.65) / (tfra-tcrit)
            endif
c
        endif
c
c assign final albedo based on snow depth
        if (snwdpt(ipt,jrow) .le. 0.0) then
            ca(nc) = abase
          elseif (snwdpt(ipt,jrow) .ge. 2.0) then
            ca(nc) = asnow
          else
            ca(nc) = abase + sqrt(snwdpt(ipt,jrow)/2.0)*(asnow-abase)
        endif
c
      endif
c
      return
      end

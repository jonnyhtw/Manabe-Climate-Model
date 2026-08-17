      subroutine bciol ( dayrad )
#if defined fixed_sst || defined q_flux || defined restore
c
c          called from : airsea
c         entry points : none
c   subroutines called : none
c
c  this routine reads monthly mean values for prescribed
c  boundary conditions such as to (mixed layer temperature), 
c  yce (sea ice thickness), conc (sea ice concentration) and
c  qflx (heat flux adjustment) as appropriate.  these are 
c  interpolated linearly to obtain values for each day (dayrad).
c
#include "parm.h"
c
#include "bcs.h"
c
#include "calndr.h"
c
#include "grdfld.h"
c
#include "param.h"
c
      dimension dmidmo(12)
#if defined fixed_sst || defined restore
      dimension sst(ix,12), ssice(ix,12)
#ifdef ice_dynamics
      dimension sconc(ix,12)
#endif
#endif
#ifdef q_flux
      dimension qflux(ix,12)
#endif
c
c disk file for boundary condition input (41)
      data iunitr / 41 /
c
c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c
c compute the day number at the middle of each calendar month
      do 500 mo=1,12
        dmidmo(mo) = 0.5*(dy(mo+1)+dy(mo))
  500 continue
c
c determine the calendar months whose mid-month dates lie just behind
c and ahead of the current date
      if (dayrad.le.dmidmo(1) .or. dayrad.gt.dmidmo(12)) then
          mbhd = 12
          mahd = 1
        else
          do 510 mo=1,11
            if (dayrad.gt.dmidmo(mo) .and. dayrad.le.dmidmo(mo+1)) then
              mbhd = mo
              mahd = mo+1
            endif
  510     continue
      endif
c
c load the value of the day number at the middle of the month behind
c and month ahead,
      if (dayrad.lt.dmidmo(mbhd)) then
c
c if the current date is before mid-january, then adjust so
c dmbhd < dayrad
          dmbhd = dmidmo(mbhd) - 365.0
          dmahd = dmidmo(mahd)
        elseif (dayrad.gt.dmidmo(mahd)) then
c
c if the current date is after mid-december, then adjust so
c dmahd > dayrad
          dmbhd = dmidmo(mbhd)
          dmahd = dmidmo(mahd) + 365.0
        else
          dmbhd = dmidmo(mbhd)
          dmahd = dmidmo(mahd)
      endif
c
c compute weights used in linear interpolation between dmbhd & dmbhd
      wmahd = (dayrad-dmbhd)/(dmahd-dmbhd)
      wmbhd = 1.0 - wmahd
c
c read the boundary condition data from disk (unit 41)
      rewind iunitr
      do 900 j=1,il
#if defined fixed_sst || defined restore
        read ( iunitr ) sst , ssice
#ifdef ice_dynamics
     $                , sconc
#endif
#endif
c
#ifdef q_flux
        read ( iunitr ) qflux
#endif
c
c load bc's with new boundary data
        do 4021 i=1,ix
#if defined fixed_sst || defined restore
          bcto(i,j) = sst(i,mahd)*wmahd + sst(i,mbhd)*wmbhd
          bcyce(i,j) = ssice(i,mahd)*wmahd + ssice(i,mbhd)*wmbhd
#ifdef ice_dynamics
          bcconc(i,j) = sconc(i,mahd)*wmahd + sconc(i,mbhd)*wmbhd
#endif
#endif
c
#ifdef q_flux
          bcqflx(i,j) = qflux(i,mahd)*wmahd + qflux(i,mbhd)*wmbhd
#endif
 4021   continue
c
 900  continue
c
#if defined fixed_sst || defined restore
c
c make certain sst and sea ice are physically consistent
c if inconsistency exists, adjust sst
      do 520 j=1,il
        do 520 i=1,ix
c
c test for sea ice thickness < 0.
         bcyce(i,j) = max(bcyce(i,j), 0.0)
c
c test for mixed layer temperature < tfro
         bcto(i,j) = max(bcto(i,j), tfro)
c
c test for sea ice thickness > 0 and mixed layer temperature > tfro
         if (bcyce(i,j).gt.0.0) then
           bcto(i,j) = tfro
         endif
#ifdef ice_dynamics
c
c test for sea ice concentration < 0.
         bcconc(i,j) = max(bcconc(i,j), 0.0)
c
c allow non-zero sea ice concentration only where thickness > 0
c and vice versa
         if (bcyce(i,j).eq.0.0) then
           bcconc(i,j) = 0.0
         endif
c
         if (bcconc(i,j).eq.0.0) then
           bcyce(i,j) = 0.0
         endif
#endif
  520 continue
#endif
c
c set bcto, bcyce and bcconc for land points to 290, 0 & 0
      do 350 j=1,il
        do 350 i=1,ix
          if (ts(i,j) .gt. 0.0) then
#if defined fixed_sst || defined restore
            bcto(i,j)   = 290.0
            bcyce(i,j)  =   0.0
#ifdef ice_dynamics
            bcconc(i,j) =   0.0
#endif
#endif
c
#ifdef q_flux
            bcqflx(i,j) =   0.0
#endif
          endif
  350 continue
c
#endif
      return
      end

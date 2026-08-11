      subroutine zenith(fjd,r,dlt,slag,xlat,ha,dhr,nlng,cosz,frac)
c
c  zenith is called by impt in diurnal models (see notes of may 17,'65).
c   zenith computes effective mean cosine of zenith angle and daylight
c   fraction from latitude and parameters given by subroutine solar
c   by means of the following equation -
c
c      cosz = sin(xlat)*sin(dlt) + cos(xlat)*cos(dlt)*sin(arg)/arg
c
c   definition of variables appearing in this subroutine
c
c     fjd   = day fraction beyond integral julian day(begins at noon ut)
c     r     = radius vector (distance to sun in astronomical units)
c     dlt   = declination of sun
c     slag  = lag in ra of apparent sun behind mean sun
c     xlat  = latitude
c     ha    = hour angle of sun at sunset
c     dhr   = integration (or averaging) period in hours
c     nlng  = number of gridpoints in a row
c     cosz  = mean cosine of zenith angle for all longitudes
c     frac  = daylight fraction at all longitudes
c     cvpr  = difference in longitude between gridpoints
c     gha   = hour angle of sun at greenwich (west of meridian is plus)
c     arg   = half of integration period in radians
c     xlng  = longitude of gridpoint
c     hloc  = hour angle of sun plus arg at longitude xlng
c     delsh = truncated averaging period resulting from late sunrise or
c               early sunset (dele and delw are similarly truncated)
c
c   input arguments to circular functions are in radians
c
      logical rise, set
c
#include "param.h"
c
      dimension cosz(nlng), frac(nlng)
c
c
c  compute basic constants used by this subroutine
c
      cvpr   = twopi / float(nlng)
      gha    = fjd * twopi + slag
      arg    = dhr * pid24
      sinfac = sin(arg) / arg
      ss     = sin(xlat) * sin(dlt)
      cc     = cos(xlat) * cos(dlt)
      if (ha .gt. 0.0) cons = ss + cc * sin(ha) / ha
c
c  compute cosz and frac at all longitudes
c
      do 20 i = 1,nlng
        xlng = cvpr * (float(i) - 1.0)
c
c  local hour angle shifted by half of the averaging period
c   to make average representative of period.
c
        hloc = gha + xlng + arg
        rise = .false.
        set  = .false.
c
c  confine hloc to within plus and minus pi
c
        hloc = amod(hloc,twopi)
        if (hloc .gt. pi) hloc = hloc - twopi
c
c  determine if sun rises or sets during averaging period
c
        hlpar = hloc + arg
        armhl = arg - hloc
        if (hlpar .gt. ha) set  = .true.
        if (armhl .gt. ha) rise = .true.
        if (rise .and. set) then
          if (ha .le. 0.0) then
            cosz(i) = 0.0
            frac(i) = 0.0
          else
c
c  averaging period covers the entire duration of daylight
c
            cosz(i) = cons
            frac(i) = ha / arg
          endif
        else
          if (hlpar .gt. pi) then
            dele = 0.5 * amax1(hlpar + ha - twopi, 0.0)
            delw = 0.5 * amax1(ha + armhl, 0.0)
          elseif (armhl .gt. pi) then
            dele = 0.5 * amax1(ha + hlpar, 0.0)
            delw = 0.5 * amax1(armhl + ha - twopi, 0.0)
          endif
          if (hlpar.gt.pi .or. armhl.gt.pi) then
c
c  averaging interrupted by both sunset and sunrise - done in 2 parts
c
            frac(i) = (dele + delw) / arg
            if (frac(i) .eq. 0.0) then
              cosz(i) = 0.0
            else
              cosz(i) = ss + cc * (cos(ha - dele) * sin(dele)
     *                             + cos(ha - delw) * sin(delw))
     *                          / (dele + delw)
            endif
          endif
        endif
        if (.not.(rise.and.set) .and. (hlpar.le.pi)
     *      .and. (armhl.le.pi)) then
          if (set) delsh = 0.5 * (ha + armhl)
          if (rise) delsh = 0.5 * (ha + hlpar)
c
c   either a sunrise or sunset occurs during the averaging period
c
          if (rise .or. set) then
            if (delsh .le. 0.0) then
              cosz(i) = 0.0
              frac(i) = 0.0
            else
              cosz(i) = ss + cc * cos(ha -delsh) * sin(delsh) / delsh
              frac(i) = delsh / arg
            endif
          else
c
c  normal case where sun is up during entire averaging period
c
            cosz(i) = ss + cc * cos(hloc) * sinfac
            frac(i) = 1.0
          endif
        endif
c
c  limit cosz and fraction to within range of zero to one
c
        cosz(i) = amin1(1.0, cosz(i))
        cosz(i) = amax1(0.0, cosz(i))
        frac(i) = amin1(1.0, frac(i))
   20 continue
c
      return
      end

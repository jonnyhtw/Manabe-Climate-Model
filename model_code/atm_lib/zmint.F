      subroutine zmint
c
#include "parm.h"
#include "parmkp.h"
#include "parmjx.h"
c
      parameter( ilt = 2 * iy + 1 )
c
#include "ctlblk.h"
c
#include "gauss.h"
c
#include "lats.h"
c
#include "param.h"
c
#include "intgrl.h"
c
c
c   <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
c
      ha   = 0.5 * float(3 - nhem) * erad
      drag = ha * drvor
      hmf  = float(3 - nhem)
      write(6,8) ktau
      rain = 0.0
      snow = 0.0
      t9   = 0.0
      evap = 0.0
      u3   = 0.0
c
      do 3 l=1,il
        ls   = min0(l, ilt - l)
        xlat = ylat(l)
        cnvrl = -0.5 * cp * zmcnvr(l) / grav
        dragl = ha * zmdrag(l)
        write(6,7) xlat,cnvrl,dragl,zmt9(l),zmpc(l),zmrn(l),
     *             zmsn(l),zmevap(l),zmu3(l)
        wt   = 0.5 * hmf * cosqr(ls) * weight(ls)
        rain = rain + wt * zmrn  (l)
        snow = snow + wt * zmsn  (l)
        t9   = t9   + wt * zmt9  (l)
        evap = evap + wt * zmevap(l)
        u3   = u3   + wt * zmu3  (l)
    3 continue
c
      prec = rain + snow
      write(6,9) cnvr,drag,t9,prec,rain,snow,evap,u3
c
c
    7 format(1x,f8.3,8(1x,e14.7))
    8 format(//,44x,'zonal mean quantities at tau=',i8,//,
     * '  latitude  conversion    surface drag  temperature(kx) '
     * 'precipitation    rainfall      snowfall      evaporation  '
     * 'zonal wind k=3',/)
    9 format(/' global  ',8(1x,e14.7),//)
c
      return
      end

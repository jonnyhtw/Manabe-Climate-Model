      subroutine cldalb (zangle, zza)
      implicit real*8 (a-h,o-z)
c
c this subroutine is called by impt and introduces a zenith angle
c dependence to cloud albedos.  the cloud albedos are interpolated
c using data adapted from fritz (1954).  the solar zenith angle
c is the only input required.
c
#include "parm.h"
#include "parmncc.h"
#ifdef fixcld 
      parameter ( lx=ncc+1 )
#endif
#ifdef prdcld
      parameter ( lx=kx+2 )
#endif
      parameter ( nbands=2, nangs=17 )
c
      dimension zza(lx, nbands)
CSGI
      dimension albch1(nangs,2), albcm1(nangs,2), albcl1(nangs,2)
CSGI
c
#include "cabd.h"
c
      common / clalbs / albch(nangs,2), albcm(nangs,2), albcl(nangs,2)
c
c high cloud albedos for zenith angles from 0-80 deg. by 5 degs.
CSGI
cc    data albch
      data albch1
     $  /.04,.05,.05,.05,.06,.06,.07,.07,.08,.11,.13,.16,.21,
     $   .28,.39,.48,.61,.04,.05,.05,.05,.06,.06,.07,.07,.08,.10,.11,
     $   .14,.19,.26,.35,.44,.55/
c
c middle cloud albedos (same format as above)
cSGI
c     data albcm
      data albcm1
     $  /.18,.18,.19,.20,.21,.23,.24,.26,.29,.33,.37,.42,.47,
     $   .55,.64,.71,.79,.14,.14,.15,.16,.17,.18,.18,.20,.23,.25,.29,
     $   .32,.37,.43,.50,.55,.61/
c
c low cloud albedos (same format as above)
cSGI
cc    data albcl
      data albcl1
     $      /.50,.50,.51,.51,.52,.53,.54,.56,.58,.62,.65,.67,.69,
     $       .73,.78,.82,.86,.42,.42,.43,.43,.44,.45,.46,.48,.50,.52,
     $       .55,.57,.59,.63,.66,.70,.74/
c
      albch = albch1
      albcm = albcm1
      albcl = albcl1
c
c if zenith angle > 80 degrees, use albedos for zenith angle = 80
      if (zangle .ge. 80.0) then
        do 1 nband=1,nbands
          zza(2,nband) = albch(nangs,nband)
          zza(3,nband) = albcm(nangs,nband)
          zza(4,nband) = albcl(nangs,nband)
    1   continue
        return
      endif
c
c otherwise, interpolate albedos from tables for each cloud level
      indx   = ifix (zangle/5.0) + 1
      remain = amod (zangle,5.0)
      do 2 nband=1,nbands
        zza(2,nband) = albch(indx,nband) +
     $    (remain/5.0) * (albch(indx+1,nband)-albch(indx,nband))
        zza(3,nband) = albcm(indx,nband) +
     $    (remain/5.0) * (albcm(indx+1,nband)-albcm(indx,nband))
        zza(4,nband) = albcl(indx,nband) +
     $    (remain/5.0) * (albcl(indx+1,nband)-albcl(indx,nband))
    2 continue
c
      return
      end

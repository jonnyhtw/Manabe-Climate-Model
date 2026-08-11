      subroutine swnew
c
#include "parm.h"
#include "parmkp.h"
#include "parmrad.h"
      parameter (nbsw=9)
c
c
#include "cld.h"
#include "cld1.h"
#include "maimc.h"
#include "swnuin.h"
#include "swnuot.h"
      dimension         f(kp),   h(kp)
      equivalence       (f,flx), (h,heat)
c
#include "swtabl0.h"
c
      dimension dfn(nbsw,kp),  ufn(nbsw,kp),  ttd(nbsw,kp)
      dimension ttu(nbsw,kp),  t1(nbsw,kp),   t2(nbsw,kp)
      dimension cca(nbsw,k2),  ccb(nbsw,k2),  ct(nbsw,k2)
      dimension cr(nbsw,k2),   tclu(nbsw,kp), tcld(nbsw,kp)
      dimension alfat(nbsw,kp),tdh2o(nbsw,kp),tuh2o(nbsw,kp)
      dimension udrv(nbsw,kx), ddrv(nbsw,kx), alfau(nbsw,kp)
c
      dimension press(kp), ppress(kp), rh2o(kp), ro3(kp)
      equivalence (press,p(1)), (ppress,pp(0))
      equivalence (rh2o,r(1,1)), (ro3,r(1,3))
      equivalence (solrsw,solar)
c
      dimension abcff(nbsw),pwts(nbsw)
      dimension df(kp),uf(kp)
      dimension dpsw(kp),ur(kp),ud(kp),dusw(kp)
      dimension ff(kp),alfa(nbsw)
      dimension vv(nbsw),vvu(nbsw),vvd(nbsw)
      dimension udo3(kp),uro3(kp),ffo3(kp), absdo3(kp),
     1    absuo3(kp),duo3(kp)
      dimension cro3(k2),cto3(k2)
      dimension udco2(kp),urco2(kp),ffco2(kp),adco2(kp),auco2(kp),
     1     duco2(kp)
      dimension axx(kpx2),ayy(kpx2),azz(kpx2)
      dimension in(kpx2),ixx(kpx2)
      dimension aco2(kpx2),ao3(kpx2),uco2(kpx2),uo3(kpx2)
      dimension dpcld(kx)
c
      equivalence (adco2(1),aco2(1)),(auco2(1),aco2(k2)),
     1            (udco2(1),uco2(1)),(urco2(1),uco2(k2)),
     2            (absdo3(1),ao3(1)),(absuo3(1),ao3(k2)),
     3            (udo3(1),uo3(1)),(uro3(1),uo3(k2))
c
      data abcff/2*4.0e-5,.002,.035,.377,1.95,9.40,44.6,190./
      data pwts/.5000,.1470,.0698,.1443,.0584,.0335,.0225,.0158,.0087/
      data cfco2,cfo3/508.96,466.64/
      data reflo3/1.9/
      data rrayav/0.144/
      data g1/1.020408e-3/
c
c     ..................................................................
c
      do 351 i=1,kp
      ff(i)=1.66
      ffco2(i)=1.66
      ffo3(i)=1.90
351   continue
      kthsw(1)=1
      kbhsw(1)=1
      kthsw(nc2)=kp
      kbhsw(nc2)=kp
      do 339 k=1,nc2
      cro3(k)=coca(k)*cloudy(k)
339   continue
      do 334 kk=1,nc2
      cto3(kk)=1.-cro3(kk)
334   continue
      do 335 n=2,nbsw
      do 335 k=1,nc2
      cca(n,k)=cwca(k)
      ccb(n,k)=cwcb(k)
335   continue
      do 336 n=2,nbsw
      do 333 kk=1,nc2
      cr(n,kk)=cca(n,kk)*cloudy(kk)
      ct(n,kk)=cloudy(kk)*(1.-(cca(n,kk)+ccb(n,kk)))+
     1     1.-cloudy(kk)
333   continue
336   continue
      do 337 kk=1,nc2
      cr(1,kk)=cro3(kk)
      ct(1,kk)=cto3(kk)
337   continue
c     kbhsw,kthsw are cloud level indices for bottom and top of cloud
c     cr,ct are cloud reflectivity,transmissivity with index 2
c     representing the topmost cloud
c     ***********************************************************
c     calculate initial ozone reflectivity
      rray=0.219/(1.+0.816*cosz)
      rg=cro3(nc2)
      refl=rray+(1.-rray)*(1.-rrayav)*rg/(1.-rg*rrayav)
      cr(1,nc2)=refl
      ltop=kthsw(2)
      secz=1./cosz
      do 2 i=1,kx
2     dpsw(i)=ppress(i+1)-ppress(i)
      do 3 i=1,kx
      duo3(i)=ro3(i)*dpsw(i)*g1
      duco2(i)=rco2*dpsw(i)*g1
3     dusw(i)=rh2o(i)*dpsw(i)*g1
      do 131 i=1,ltop
      ffo3(i)=secz
      ffco2(i)=secz
      ff(i)=secz
131   continue
c     calculate pressure-weighted optical paths in units of g/cm2.
c     pressure weighting is by p**0.5
c     ud is the downward path,ur the upward path,
c     and the calculation is made by taking a path with an angle
c     of (secz) from the top of the atmosphere to the topmost cloud,
c     then using the diffusivity factor (1.66) to the surface and for
c     reflected radiation. the code below reflects this.
      ud(1)=0.
      udco2(1)=0.
      udo3(1)=0.
      do 4 i=2,kp
      ud(i)=ud(i-1)+dusw(i-1)*pr2(i-1)*ff(i)
      udco2(i)=udco2(i-1)+duco2(i-1)*pr2(i-1)*ffco2(i)*cfco2
      udo3(i)=udo3(i-1)+duo3(i-1)*ffo3(i)*cfo3
4     continue
c     udo3,uro3 are in units of cm. cfco2,cfo3 is the conversion
c     factor from gm/cm2 to cm.
      ur(kp)=ud(kp)
      urco2(kp)=udco2(kp)
      uro3(kp)=udo3(kp)
      do 5 i=1,kx
      ur(i)=ud(kp)+1.66*(ud(kp)/ff(kp)-ud(i)/ff(i+1))
     $          +ud(ltop)*(ff(kp)-ff(i+1))/ff(i+1)
      urco2(i)=urco2(kp)+1.66*(udco2(kp)/ffco2(kp)-udco2(i)/ffco2(i+1))
     $   +udco2(ltop)*(ffco2(kp)-ffco2(i+1))/ffco2(i+1)
      uro3(i)=uro3(kp)+reflo3*(udo3(kp)/ffo3(kp)-udo3(i)/ffo3(i+1))
     $     +udo3(ltop)*(ffo3(kp)-ffo3(i+1))/ffo3(i+1)
5     continue
c
c     maximize the size of o3 and co2 path lengths to avoid going
c     off tables
      do 621 k=1,kpx2
       uco2(k) =amin1(uco2(k),998.9)
       uo3(k) =amin1(uo3(k),3.998)
621   continue
c
c     calculate entering flux at the top
      do 6 n=1,nbsw
      dfn(n,1)=solar*6.97667e5*cosz*tauda*pwts(n)
6     continue
c     calculate water vapor transmission functions for bands 2-9;
c     t.f. for band 1= t.f for band 2
      do 7 k=1,kp
      do 7 n=2,nbsw
      t1(n,k)=amin1(abcff(n)*ud(k),50.)
      t2(n,k)=amin1(abcff(n)*ur(k),50.)
7     continue
      do 8 n=2,nbsw
      tdh2o(n,1)=1.
8     continue
      do 9 n=2,nbsw
      do 9 k=2,kp
      tdh2o(n,k)=exp(-t1(n,k))
9     continue
      do 10 n=2,nbsw
      tuh2o(n,kp)=tdh2o(n,kp)
10    continue
      do 11 n=2,nbsw
      do 11 k=1,kx
      tuh2o(n,k)=exp(-t2(n,k))
11    continue
      do 12 k=1,kp
      tdh2o(1,k)=tdh2o(2,k)
      tuh2o(1,k)=tuh2o(2,k)
12    continue
c     calculate co2 absorptions . they will be used in bands 2-9.
c     since these occupy 50 percent of the solar spectrum the
c     absorptions will be multiplied by 2.
c     the absorptions are obtained by table lookup in array aab,
c     common block swtabl,and then interpolation.
      do 614 k=1,kpx2
      ixx(k)=uco2(k)+1.
614   continue
      do 615 k=1,kpx2
      axx(k)=uco2(k)-ifix(uco2(k))
615   continue
      call sif1d ( ixx, kpx2m, aab, ayy, azz )
      do 617 k=1,kpx2
      aco2(k)=ayy(k)+axx(k)*(azz(k)-ayy(k))
617   continue
      do 26 k=1,kp
      adco2(k)=2.*adco2(k)
      auco2(k)=2.*auco2(k)
26    continue
c     now calculate ozone absorptions. these will be used in
c     band 1. as this occupies 50 percent of the solar spectrum
c     the ozone absorptions will be multiplied by 2.
c     the ozone absorptions are obtained by table lookup from
c     array aaa, common block swtabl. no interpolation is done.
      do 603 k=1,kpx2
      in(k)=500.*uo3(k)+1
603   continue
      call sif1 ( in, kpx2m, aaa, ao3 )
      do 33 k=1,kp
      absdo3(k)=2.*absdo3(k)
      absuo3(k)=2.*absuo3(k)
33    continue
c     combine absorptions and transmissions to obtain a
c     transmission function for each of the 9 bands.
      do 41 n=1,nbsw
      ttd(n,1)=1.
41    continue
      do 42 k=2,kp
      ttd(1,k)=tdh2o(1,k)*(1.-absdo3(k))
42    continue
      do 43 k=1,kx
      ttu(1,k)=tuh2o(1,k)*(1.-absuo3(k))
43    continue
      do 44 n=2,nbsw
      do 44 k=2,kp
      ttd(n,k)=tdh2o(n,k)*(1.-adco2(k))
44    continue
      do 45 n=1,nbsw
      ttu(n,kp)=ttd(n,kp)
45    continue
      do 46 n=2,nbsw
      do 46 k=1,kx
      ttu(n,k)=tuh2o(n,k)*(1.-auco2(k))
46    continue
c     the following calculation is for alfat: the ratio between
c     the downward flux at the top of the highest cloud (if
c     present) or the ground to the upward flux at the same
c     level, taking into account multiple reflections from
c     clouds, if present
      do 53 n=1,nbsw
      do 53 nn=1,nc1
      alfat(n,nn)=cr(n,nc3-nn)
53    continue
      do 55 nn=1,nc1
      do 55 n=1,nbsw
      alfau(n,nn)=0.
55    continue
      if ( nc0 .eq. 0 )  go to 58
      do 51 ll=1,nc1
      kt1=kthsw(nc3-ll)
      kt2=kthsw(nc2-ll)
      kt3=kbhsw(nc2-ll)
c     tclu(ll) is transmission function from top of next lower
c     cloud to top of upper cloud. tcld(ll) is t.f. from top
c     of next lower cloud to bottom of upper cloud. ll=1 is
c     the lowest bottom cloud (the ground) ; ll=nc1 is the
c     highest upper cloud. this calculation is used for the
c     alfat calculation only if cloud is present.
      do 52 n=1,nbsw
      tclu(n,ll)=ttd(n,kt1)/ttd(n,kt2)
      tcld(n,ll)=ttd(n,kt1)/ttd(n,kt3)
52    continue
51    continue
      do 56 nn=1,nc0
      do 57 n=1,nbsw
      alfau(n,nn+1)=(cr(n,nc3-nn)+alfau(n,nn))*ct(n,nc2-nn)*
     1     ct(n,nc2-nn)*tclu(n,nn)*tclu(n,nn)/
     2     (1.-(cr(n,nc3-nn)+alfau(n,nn))*cr(n,nc2-nn)*
     3     tcld(n,nn)*tcld(n,nn))
57    continue
56    continue
58    continue
      do 59 n=1,nbsw
      alfa(n)=alfat(n,nc1)+alfau(n,nc1)
59    continue
c     downward flux above topmost cloud
      do 61 n=1,nbsw
      vv(n)=dfn(n,1)
61    continue
      do 62 n=1,nbsw
      do 62 k=2,ltop
      dfn(n,k)=vv(n)*ttd(n,k)
62    continue
c     upward flux above topmost cloud
      ltopm=ltop-1
      do 63 n=1,nbsw
      ufn(n,ltop)=alfa(n)*dfn(n,ltop)
63    continue
      do 64 n=1,nbsw
      vv(n)=ufn(n,ltop)/ttu(n,ltop)
64    continue
      do 65 n=1,nbsw
      do 65 k=1,ltopm
      ufn(n,k)=vv(n)*ttu(n,k)
65    continue
      if ( nc0 .eq. 0 )  go to 91
c     calculate ufn at cloud tops and dfn at cloud bottoms
      do 71 ll=2,nc1
      do 72 n=1,nbsw
      ufn(n,kthsw(ll+1))=alfau(n,nc3-ll)*dfn(n,kbhsw(ll-1))*
     1     tcld(n,nc3-ll)/(tclu(n,nc2-ll)*ct(n,ll))
72    continue
      do 73 n=1,nbsw
      dfn(n,kbhsw(ll))=dfn(n,kbhsw(ll-1))*tcld(n,nc3-ll)*tclu(n,nc2-ll)*
     1     ct(n,ll)/tcld(n,nc2-ll)+ufn(n,kthsw(ll+1))*tcld(n,nc2-ll)*
     2     cr(n,ll)
73    continue
71    continue
c     now obtain dfn and ufn for levels between the clouds
      do 74 ll=1,nc0
      ltop=kbhsw(ll+1)
      ltopp=ltop+1
      lbot=kthsw(ll+2)
      lbotm=lbot-1
      if (ltop.eq.lbot) go to 74
      do 75 n=1,nbsw
      vvu(n)=ufn(n,lbot)/ttu(n,lbot)
      vvd(n)=dfn(n,ltop)/ttd(n,ltop)
75    continue
      do 76 k=ltop,lbotm
      do 76 n=1,nbsw
      ufn(n,k)=vvu(n)*ttu(n,k)
76    continue
      do 77 k=ltopp,lbot
      do 77 n=1,nbsw
      dfn(n,k)=vvd(n)*ttd(n,k)
77    continue
74    continue
c     now obtain downward and upward fluxes for levels,if any,
c     between the tops and bottoms of clouds. the assumption of
c     constant heating rate in these regions is used.
      do 78 ll=1,nc0
      ltop=kthsw(ll+1)
      lbot=kbhsw(ll+1)
      if((lbot-ltop).le.1) go to 78
      ltopp=ltop+1
      lbotm=lbot-1
      dpcld(ll)=ppress(ltop)-ppress(lbot)
      do 79 n=1,nbsw
      udrv(n,ll)=(ufn(n,ltop)-ufn(n,lbot))/dpcld(ll)
      ddrv(n,ll)=(dfn(n,ltop)-dfn(n,lbot))/dpcld(ll)
79    continue
      do 80 n=1,nbsw
      vvu(n)=ufn(n,ltop)
      vvd(n)=dfn(n,ltop)
80    continue
      do 81 k=ltopp,lbotm
      do 81 n=1,nbsw
      ufn(n,k)=vvu(n)+udrv(n,ll)*(ppress(k)-ppress(ltop))
      dfn(n,k)=vvd(n)+ddrv(n,ll)*(ppress(k)-ppress(ltop))
81    continue
78    continue
91    continue
c     sum over bands
c
c
      do 14 k=1,kp
      df(k)=0.
      uf(k)=0.
      do 15 n=1,nbsw
      df(k)=df(k)+dfn(n,k)
      uf(k)=uf(k)+ufn(n,k)
15    continue
14    continue
      do 16 k=1,kp
      f(k)=uf(k)-df(k)
16    continue
      do 17 k=1,kx
      h(k)=8.42668*(f(k+1)-f(k))/dpsw(k)
17    continue
c     8.42668=g(cgs units)*(no.sec/da)/cp(cgs units)
      grdflx=(1.-refl)*dfn(1,kp)
      do 19 n=2,nbsw
      grdflx=grdflx+ct(1,nc2)*dfn(n,kp)
19    continue
      return
      end

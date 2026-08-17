c====================== include file "grdvar.h" ========================
c
c     variables which are functions of the grid defined by "coord.h"
c
c     dxt     = longitudinal width of "t" grid box at the
c               equator (in cm)
c     dxtr    = reciprocal of "dxt"
c     dxt2r   = reciprocal of "2*dxt"
c     dxt4r   = reciprocal of "4*dxt"
c     dxu     = longitudinal width of "u,v" grid box at the
c               equator (in cm)
c     dxur    = reciprocal of "dxu"
c     dxu2r   = reciprocal of "2*dxu"
c     dxu4r   = reciprocal of "4*dxu"
c
c     dyt     = latitudinal height of "t" grid box (in cm)
c     dytr    = reciprocal of "dyt"
c     dyt2r   = reciprocal of "2*dyt"
c     dyt4r   = reciprocal of "4*dyt"
c     dyu     = latitudinal height of "u,v" grid box (in cm)
c     dyur    = reciprocal of "dyu"
c     dyu2r   = reciprocal of "2*dyu"
c     dyu4r   = reciprocal of "4*dyu"
c
c     csu     = cosine of "u,v" grid point latitude
c     csur    = reciprocal of "csu"
c     cst     = cosine of "t" grid point latitude
c     cstr    = reciprocal of "cst"
c     phi     = latitude of "u,v" grid point in radians
c     phit    = latitude of "t" grid point in radians
c     sine    = sine of "u,v" grid point latitude
c     tng     = tan of "u,v" grid point latitude
c     fcor    = 2*omega*sine(j)
c
c     dzt(k)  = level thickness of "t" and "u,v" grid boxes (in cm)
c               (vertical separation between "w" velocity points)
c     c2dzt(k)= "2*dzt"
c     dztr(k) = reciprocal of dzt
c     dzt2r(k)= reciprocal of "2*dzt"
c     dzw(k)  = vertical separation between centers of levels k & k+1
c     dzwr(k) = reciprocal of dzw
c     dzw2r(k)= reciprocal of "2*dzw"
c     dztur(k)= upper diffusion grid factor = 1.0/(dzw(k-1)*dzt(k))
c     dztlr(k)= lower diffusion grid factor = 1.0/(dzw(k)*dzt(k))
c     dzwur(k)= upper diffusion grid factor = 1.0/(dzt(k)*dzw(k))
c     dzwlr(k)= lower diffusion grid factor = 1.0/(dzt(k+1)*dzw(k))
c     tanra(j)= tan(j)/radius on the "u,v" grid
c
      common /grdvar/ dxt(imt),  dxtr(imt),  dxt2r(imt), dxu(imt)
     $,               dxur(imt), dxu2r(imt), dxu4r(imt), dxt4r(imt)
     $,               dyt(jmt),  dytr(jmt),  dyt2r(jmt), dyu(jmt)
     $,               dyur(jmt), dyu2r(jmt), dyu4r(jmt), dyt4r(jmt)
     $,               csu(jmt),  csur(jmt),  cst(jmt),   cstr(jmt)
     $,               phi(jmt),  phit(jmt),  sine(jmt),  tng(jmt)
     $,               fcor(jmt), c2dzt(km),  dzt(km),    dztr(km)
     $,               dzt2r(km), dzw(0:km),  dzwr(0:km), dzw2r(0:km)
#ifdef tcvmix
     $,               dzwur(km), dzwlr(km)
#endif
#if defined implicitvmix || defined isopycmix
     $,               dztur(km), dztlr(km)
#endif
#ifndef nohilats
     $,               tanra(jmt)
#endif
c

c====================== include file "fdift.h" =========================
c
c     finite difference numerics for tracers
c
c
c     zonal advection of tracer
c
      UTx(i,k) = (aux1(i,k) - aux1(i+1,k))*dxt4r(i)
c
c     meridional advection of tracer
c
      VTy(i,k) = fvst(i,k)*(t(i,k,jc,nc,m) + t(i,k,jm,nc,m))
     $          -fvn (i,k)*(t(i,k,jp,nc,m) + t(i,k,jc,nc,m))
c
c     vertical advection of tracer
c
      WTz(i,k) = (fwb1(i,k) - fwb1(i,k-1))*dzt2r(k)
#ifndef isopycmix
# ifdef biharmonic
c
c     zonal diffusion (del**4) of tracer
c
      Txx(i,k) = bbt(j)*dxt4r(i)*
     $           (fm(i-1,k,jc)*aux2(i,k) - fm(i+1,k,jc)*aux2(i+1,k))
c
c     meridional diffusion (del**4)  of tracer
c
      Tyy(i,k) = cct(j)*fm(i,k,jp)*(del2(i,k,m,jpt2)-del2(i,k,m,jpt3))
     $         + ddt(j)*fm(i,k,jm)*(del2(i,k,m,jpt2)-del2(i,k,m,jpt1))
#  define hmixalreadyset
# endif
# if defined nlhmix && !defined hmixalreadyset
c
c     zonal diffusion (nonlinear) of tracer
c
      Txx(i,k) = bbtj(i,k)*dxt4r(i)*
     $           (fm(i+1,k,jc)*aux2(i+1,k) - fm(i-1,k,jc)*aux2(i,k))
c
c     meridional diffusion (nonlinear) of tracer
c
      Tyy(i,k) = cctj(i,k)*fm(i,k,jp)*(t(i,k,jp,nm,m)-t(i,k,jc,nm,m))
     $         + ddtj(i,k)*fm(i,k,jm)*(t(i,k,jm,nm,m)-t(i,k,jc,nm,m))
#  define hmixalreadyset
# endif
# if defined consthmix && !defined hmixalreadyset
c
c     zonal diffusion (del**2) of tracer
c
      Txx(i,k) = hdcsck(k)*bbt(j)*dxt4r(i)*
     $           (fm(i+1,k,jc)*aux2(i+1,k) - fm(i-1,k,jc)*aux2(i,k))
c
c     meridional diffusion (del**2) of tracer
c
      Tyy(i,k) = (cct(j)*fm(i,k,jp)*(t(i,k,jp,nm,m)-t(i,k,jc,nm,m))
     $         +  ddt(j)*fm(i,k,jm)*(t(i,k,jm,nm,m)-t(i,k,jc,nm,m)))
     $         *hdcsck(k)
#  define hmixalreadyset
# endif
#endif
#ifdef isopycmix
c
c     zonal diffusion of tracer through E & W side walls
c     (xx & zx terms included ... yx is neglected)
c
      Tisox(i,k) =  (e(i,k,1) - e(i-1,k,1))*dxtr(i)*cstr(j)
c
c     meridional diffusion of tracer through N & S side walls
c     (xy & zy terms included ... xy is neglected)
c
      Tisoy(i,k) =  (e(i,k,2) - esav(i,k,m))*dytr(j)*cstr(j)
c
c     diffusion of tracer oriented thru top & bot along isopycnal slopes
c     (xz & yz terms included ... zz term done implicitly) 
c
      Tisoz(i,k) = (e(i,k,3) - e(i,k+1,3))*dztr(k)
# define hmixalreadyset
#endif
#ifndef hmixalreadyset
      Txx(i,k) = c0
      Tyy(i,k) = c0
#endif
c
c     vertical diffusion of tracer
c
      Tzz(i,k) = (vtf(i,k-1) - vtf(i,k))*dztr(k)
#undef hmixalreadyset
c

c====================== include file "fdifm.h" =========================
c
c     finite difference numerics for momentum
c
c
c     zonal advection of momentum
c
      UUx(i,k) = (aux1(i,k) - aux1(i+1,k))*dxu2r(i)
      UVx(i,k) = (aux2(i,k) - aux2(i+1,k))*dxu2r(i)
c
c     meridional advection of momentum
c
      VUy(i,k) = fvsu(i,k)*(u(i,k,jc,nc) + u(i,k,jm,nc))
     $          -fvn (i,k)*(u(i,k,jp,nc) + u(i,k,jc,nc))
      VVy(i,k) = fvsu(i,k)*(v(i,k,jc,nc) + v(i,k,jm,nc))
     $          -fvn (i,k)*(v(i,k,jp,nc) + v(i,k,jc,nc))
c
c     vertical advection of momentum
c
      WUz(i,k) = (fwb1(i,k) - fwb1(i,k-1))*dzt2r(k)
      WVz(i,k) = (fwb2(i,k) - fwb2(i,k-1))*dzt2r(k)
#ifndef nohilats
c
c     metric advective terms for high latitudes
c
      UVtan(i,k) = u(i,k,jc,nc)*v(i,k,jc,nc)*tanra(j)
      UUtan(i,k) = u(i,k,jc,nc)**2*tanra(j)
#endif
#ifdef biharmonic
c
c     zonal mixing (del**4) of momentum
c
      Uxx(i,k) = bbu(j)*(dxu2r(i)*(aux3(i,k) - aux3(i+1,k)))
      Vxx(i,k) = bbu(j)*(dxu2r(i)*(aux4(i,k) - aux4(i+1,k)))
c
c     meridional mixing (del**4) of momentum
c
      Uyy(i,k) = ccu(j)*(del2(i,k,m1,jpt2) - del2(i,k,m1,jpt3)) +
     $           ddu(j)*(del2(i,k,m1,jpt2) - del2(i,k,m1,jpt1))
      Vyy(i,k) = ccu(j)*(del2(i,k,m2,jpt2) - del2(i,k,m2,jpt3)) +
     $           ddu(j)*(del2(i,k,m2,jpt2) - del2(i,k,m2,jpt1))
c
c     metric term for mixing (del**4) of momentum
c
      Umet(i,k) = - ggu(j)*del2(i,k,m1,jpt2) +
     $            hhu(j)*dxu2r(i)*(del2(i+1,k,m2,jpt2)
     $                           - del2(i-1,k,m2,jpt2))
      Vmet(i,k) = - ggu(j)*del2(i,k,m2,jpt2) -
     $            hhu(j)*dxu2r(i)*(del2(i+1,k,m1,jpt2)
     $                           - del2(i-1,k,m1,jpt2))
# define hmixalreadyset
#endif
#if defined nlhmix && !defined hmixalreadyset
c
c     zonal mixing (nonlinear) of momentum
c
      Uxx(i,k) = ((dpure(i+1,k)*dten(i+1,k) 
     $            + dpures(i+1,k)*dtens(i+1,k))
     $           - (dpure(i-1,k)*dten(i-1,k)
     $            + dpures(i-1,k)*dtens(i-1,k)))
     $            *dxu2r(i)*effk(i)*bbu(j)
      Vxx(i,k) = ((dpure(i+1,k)*dshear(i+1,k) 
     $            + dpures(i+1,k)*dshers(i+1,k))
     $           - (dpure(i-1,k)*dshear(i-1,k)
     $            + dpures(i-1,k)*dshers(i-1,k)))
     $            *dxu2r(i)*effk(i)*bbu(j)
c
c     meridional mixing (nonlinear)  of momentum
c
      Uyy(i,k) = (dpure(i,k)*dshear(i,k)*cst(j+1)**2
     $           - dpures(i,k)*dshers(i,k)*cst(j)**2)
     $           *ccu(j)
      Vyy(i,k) = -(dpure(i,k)*dten(i,k)*cst(j+1)**2
     $           -dpures(i,k)*dtens(i,k)*cst(j)**2)
     $           *ccu(j)
c
c     metric mixing (nonlinear) of momentum
c
      Umet(i,k) = c0
      Vmet(i,k) = (tng(j)/radius)*(effky(j)*p5
     $            *(dpure(i,k)*dten(i,k)
     $             + dpures(i,k)*dtens(i,k)))
# define hmixalreadyset
#endif
#if defined consthmix && !defined hmixalreadyset
# define hmixalreadyset
c
c     zonal mixing (del**2) of momentum
c
      Uxx(i,k) = hvcsck(k)*bbu(j)*(dxu2r(i)*(aux3(i+1,k) - aux3(i,k)))
      Vxx(i,k) = hvcsck(k)*bbu(j)*(dxu2r(i)*(aux4(i+1,k) - aux4(i,k)))
c
c     meridional mixing (del**2) of momentum
c
      Uyy(i,k) = (ccu(j)*(u(i,k,jp,nm) - u(i,k,jc,nm)) +
     $            ddu(j)*(u(i,k,jm,nm) - u(i,k,jc,nm)))*hvcsck(k)
      Vyy(i,k) = (ccu(j)*(v(i,k,jp,nm) - v(i,k,jc,nm)) +
     $            ddu(j)*(v(i,k,jm,nm) - v(i,k,jc,nm)))*hvcsck(k)
c
c     metric term for mixing (del**2) of momentum
c
      Umet(i,k) = hvcsck(k)*(ggu(j)*u(i,k,jc,nm) -
     $            hhu(j)*dxu2r(i)*(v(i+1,k,jc,nm) - v(i-1,k,jc,nm)))
      Vmet(i,k) = hvcsck(k)*(ggu(j)*v(i,k,jc,nm) +
     $            hhu(j)*dxu2r(i)*(u(i+1,k,jc,nm) - u(i-1,k,jc,nm)))
#endif
#ifndef hmixalreadyset
      Uxx(i,k)  = c0
      Vxx(i,k)  = c0
      Uyy(i,k)  = c0
      Vyy(i,k)  = c0
      Umet(i,k) = c0
      Vmet(i,k) = C0
#endif
c
c     vertical mixing of momentum
c
      Uzz(i,k) = (vmf(i,k-1,1) - vmf(i,k,1))*dztr(k)
      Vzz(i,k) = (vmf(i,k-1,2) - vmf(i,k,2))*dztr(k)
c
c     coriolis term
c
      fV(i,k) = fcor(j)*v(i,k,jc,lcor)
      fU(i,k) = fcor(j)*u(i,k,jc,lcor)
#undef hmixalreadyset
c

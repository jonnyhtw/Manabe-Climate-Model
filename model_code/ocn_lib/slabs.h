c====================== include file "slabs.h" =========================
c
c                         slab layout on disk
c
c     on disk, each model latitude is represented by a slab (or latitude
c     row of data). there are three disks and each contains a different
c     time level (one for past "n-1", present "n", & predicted "n+1" 
c     data )
c
c                         memory slab window
c
c     the memory slab window is "nslabs" wide and "ntau" time levels
c     deep. in subroutine "step", as row j cycles from 1 through jmt-1,
c     (south to north) slab data is read into the memory slab window
c     from disk "n" & disk "n-1". this data can be referenced in the
c     memory slab window using the following pointers (indices) :
c
c     jm     = points to the row j-1 data
c     jc     = points to the the central row j data
c     jp     = points to the the row j+1 data
c     jpp    = points to the row j+2 data (biharmonic option only)
c
c     jpt1   = cyclic pointer (index) from 1 to numjpt 
c              (only for internal mode & biharmonic use)
c              points to data at jrow - 1
c     jpt2   = cyclic pointer (index) from 1 to numjpt
c              (only for internal mode & biharmonic use)
c              points to data at jrow
c     jpt3   = cyclic pointer (index) from 1 to numjpt
c              (only for internal mode & biharmonic use)
c              points to data at jrow + 1
c
c     nm     = points to time level n-1 (in "timelv.h")
c     nc     = points to current time level n (in "timelv.h")
c     np     = points to time level n+1 (in "timelv.h")
c
c     pointers (indices) jm,jc,jp & jpp take on values from 1 to
c     "nslabs"
c     pointers (indices) jpt1,jpt2,jpt3 take on values from 1 to
c     "numjpt"
c     pointers (indices) nm,nc & np take on values from 1 to "ntau"
c
c     lcor   = points to time level for coriolis term(implicit/explicit)
c     m1     = points to location of "u"  within slab (for biharmonic)
c     m2     = points to location of "v"  within slab (for biharmonic)
c
c     the slab data is :
c
c     t = tracer quantites. the 1st is temperature & the 2nd is salinity
c         if nt > 2 then other tracers are allowed.
c     u = zonal component of velocity.
c     v = meridional component of velocity.
c
c         (note: only the baroclinic component of the velocity is stored
c                on the slab disks. in memory, the barotropic component
c                is added on to give the full velocity)
c
c
#ifndef multitasking
# define task
#endif
#include "timelv.h"
#ifdef biharmonic
cSGI
cc    task common /slabs/ 
           common /slabs/ jm, jc, jp, jpp, jpt1, jpt2, jpt3, lcor
     $,                   m1, m2
#else
cSGI
cc    task common /slabs/ 
           common /slabs/ jm, jc, jp, jpt1, jpt2, lcor
#endif
     $,                   t(imt,km,nslabs,ntau,nt)
     $,                   u(imt,km,nslabs,ntau)
     $,                   v(imt,km,nslabs,ntau)
c
c
c     normally, "np" is not used since predicted quantities at level
c     "np" are stored into the buffer slab before being written to disk
c
c     ta = buffer area to hold the updated "n+1" tracer slab data
c     ua = buffer area to hold the updated "n+1" u slab data
c     va = buffer area to hold the updated "n+1" v slab data
c     
c     when all prognostic slab variables have been updated, the 
c     common /bufout/ area is written to disk 
c
cSGI
cc    task common /bufout/ 
           common /bufout/ ta(imt,km,nt), ua(imt,km), va(imt,km)
c     
c     all prognostic slab variables are read from disk using the 
c     common /bufin/ area and transferred to the appropriate locations  
c     in the memory slab window via "xfer"
c
cSGI
cc    task common /bufin/ 
           common /bufin/ bufsl(imt,km,nvar)
c
c
c     uclin  = internal mode u (used for advective fluxes & diagnostics)
c              (,,1) is for "jrow"; (,,2) is for "jrow+1"
c     vclin  = internal mode v (used for advective fluxes & diagnostics)
c              (,,1) is for "jrow"; (,,2) is for "jrow+1"
c     rhon   = density at row j+1 on the "t" grid
c     rhos   = density at row j on the "t" grid
c     dpdx   = zonal gradient of pressure on "u,v" grid
c     dpdy   = meridional gradient of pressure on "u,v" grid
c     fuw    = advective coeff for western face of "u,v" grid box
c              in "clinic" & "t" grid box in "tracer"
c     fvn    = advective coeff for northern face of "u,v" grid box
c              in "clinic" & "t" grid box in "tracer"
c     fvst   = advective coeff for southern face of "t" grid box
c     fvsu   = advective coeff for southern face of "u,v" grid box
c     fm     = (0,1) over "t" grid (land,ocean) points
c     gm     = (0,1) over "u,v" grid (land,ocean) points
c     vmf    = array (1) used for vertical differences of u
c     vmf    = array (2) used for vertical differences of v
c     vtf    = array used for vertical tracer flux
c     w      = vertical velocity defined at the bottom of "u,v" boxes
c              in "clinic" and "t" boxes in "tracer"
c     fwb1   = w * (quantity 1) defined at the bottoms of the boxes
c     fwb2   = w * (quantity 2) defined at the bottoms of the boxes
c     restr  = term for restoring surface tracers to prescribed values
c              via newtonain damping
c     rests  = time scale for restoring surface tracers (days)
c
cSGI
cc    task common /work/
           common /work/
     $ uclin(imt,km,numjpt),   vclin(imt,km,numjpt)
     $,dpdx(imt,km),      dpdy(imt,km),     fuw(imt,km),   fvn(imt,km)
     $,fvsu (imt,km),     fvst (imt,km),    rhon(imt,km),  rhos(imt,km)
     $,fm(imt,km,nslabs), gm(imt,km,nslabs)
     $,vmf(imt,0:km,2),   vtf(imt,0:km)
     $,w(imt,0:km),       fwb1(imt,0:km),   fwb2(imt,0:km)
#if defined restorst || defined cpldleg2 || defined cpld_restore
     $,restr(imt,nt)
      common /workr/ rests
#endif
c
c     aux = auxiliary arrays for storing intermediate results to
c           prevent redundant calculation
c
cSGI
cc    task common /aux/ aux1(imt,km), aux2(imt,km), aux3(imt,km)
           common /aux/ aux1(imt,km), aux2(imt,km), aux3(imt,km)
     $,            aux4(imt,km)

c
c     sourcu = source term for zonal component of momentum
c     sourcv = source term for meridional component of momentum
c     sourct = source term for tracer
c
cSGI
cc    task common /source/ sourcu(imt,km), sourcv(imt,km)
           common /source/ sourcu(imt,km), sourcv(imt,km)
     $,                    sourct(imt,km,nt)
c
c     if momentum & tracer terms are to be kept (as opposed to being
c     statement functions) then they are dimensioned here.
c
#ifdef keepterms
cSGI
cc    task common /terms/
           common /terms/
     $               UTx(imt,km), VTy(imt,km), WTz(imt,km)
#ifdef isopycmix
     $,              Tisox(imt,km), Tisoy(imt,km), Tisoz(imt,km)
#else
     $,              Txx(imt,km), Tyy(imt,km)
#endif
     $,              Tzz(imt,km)
     $,              UUx(imt,km), UVx(imt,km), VUy(imt,km)
     $,              VVy(imt,km), WUz(imt,km), WVz(imt,km)
#ifndef nohilats
     $,              UVtan(imt,km), UUtan(imt,km)
#endif
     $,              Uxx(imt,km), Uyy(imt,km), Umet(imt,km)
     $,              Vxx(imt,km), Vyy(imt,km), Vmet(imt,km)
     $,              Uzz(imt,km), Vzz(imt,km), fU(imt,km), fV(imt,km)
#endif
c
c     include further extensions to the "slabs.h" here. care must
c     be exercised when adding additional variables to the slab block.
c     increase "nvar" & "nvarbh" appropriately. also be careful to
c     add the variables to "delsq"
c
#ifdef tcvmix
#include "tcslab.h"
#endif

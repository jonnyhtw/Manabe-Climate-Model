c====================== include file "emode.h" =========================
c
c
c     variables for rigid lid calculation of external mode
c
c     p     = stream function (,,1) is for tau; (,,2) is for tau-1
c     ztd   = time change of vorticity
c     ptd   = time change of stream function
c     h     = depth over "u,v" points
c     hr    = reciprocal depth over "u,v" points
c     zu    = time change of vertically averaged zonal forcing on row j
c     zv    = time change of vertically averaged merid forcing on row j
#ifdef neptune
c     pnept = streamfunction given by -fLH**2, used to "correct" the
c             transport streamfunction tendency when parameterizing
c             topographic stress effects (the Neptune effect)
c     unept = depth independent e-w velocity arising from Neptune effect
c     vnept = depth independent n-s velocity arising from Neptune effect
#endif
#ifdef sfctopog
c
c     these five 2-D arrays are associated with the computation of
c     diagnostic surface pressure fields
c
c     psx    = zonal gradient of sea surface pressure
c     psy    = meridional gradient of sea surface pressure
c     divslh = divergence of sea level times depth
c     sealev = sea level anomalies [in centimeters]
c
      common /dsprsr/ psx(imt,jmt), psy(imt,jmt)
      common /dsprsr/ divslh(imt,jmt), sealev(imt,jmt)
#endif
c
c     cfn   = relaxation coefficient to the north of point (i,j)
c     cfs   = relaxation coefficient to the south of point (i,j)
c     cfe   = relaxation coefficient to the east of point (i,j)
c     cfw   = relaxation coefficient to the west of point (i,j)
c     cpf   = relaxation coefficient for point (i,j)
c     cfne  = relaxation coefficient to the northeast of point (i,j)
c     cfse  = relaxation coefficient to the southeast of point (i,j)
c     cfnw  = relaxation coefficient to the northwest of point (i,j)
c     cfsw  = relaxation coefficient to the southwest of point (i,j)
c     rmsk  = relaxation mask of (0,1) on (land & perimeter, open ocean)
c     res   = relaxation residual
c     cof   = normalizing coefficient for islands
c     cofis = flow around islands
c
      common /emode/
     $ cfn(imt,jmt), cfs(imt,jmt), cfe(imt,jmt), cfw(imt,jmt)
     $,ptd(imt,jmt), res(imt,jmt), cpf(imt,jmt)
     $,p(imt,jmt,2), ztd(imt,jmt), hr (imt,jmt), h(imt,jmt)
#ifdef neptune
     $,pnept(imt,jmt), unept(imt,jmt), vnept(imt,jmt)
#endif
     $,zu(imt,jmt),  zv(imt,jmt)
#ifdef hypergrid
     $, rmsk(imt,jmt)
#endif
#if defined congrad5pt || defined congrad9pt
     $,gdir(imt,jmt), hdir(imt,jmt), fdir(imt,jmt), asum(jmt)
#endif
#ifdef congrad9pt
     $,cfne(imt,jmt), cfnw(imt,jmt), cfse(imt,jmt), cfsw(imt,jmt)
#endif
#ifdef islands
     $,cof(imt,jmt), cofis(nisle)
#endif
#ifdef oldcongrad
congrad
     $,ismask(imt,jmt)
#endif
c

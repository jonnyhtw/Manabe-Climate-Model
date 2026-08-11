c====================== include file "convt.h" =========================
#ifdef worldocn
c
c     convt common block stores a convection flags.
c
c     flgcnv = i,j,k points affected by the the convective adjustment
c            code are given flags of 1...those with no convection are
c            given flags of zero.  
c
      common / convt / flgcnv(imt,km)
c
#endif
c====================== end file "convt.h" =============================

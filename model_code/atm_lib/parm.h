#ifdef r15
      parameter ( jtrun=15, mtrun=15, ix=48, iy=20 )
#endif
#ifdef r21
      parameter ( jtrun=21, mtrun=21, ix=64, iy=27 )
#endif
#ifdef r30
      parameter ( jtrun=30, mtrun=30, ix=96, iy=40 )
#endif
#ifdef r63
      parameter ( jtrun=63, mtrun=63, ix=192, iy=80 )
#endif
#ifdef t21
      parameter ( jtrun=21, mtrun=21, ix=64, iy=16 )
#endif
#ifdef t42
      parameter ( jtrun=42, mtrun=42, ix=128, iy=32 )
#endif
#ifdef t84
      parameter ( jtrun=84, mtrun=84, ix=256, iy=64 )
#endif
#ifdef k9
      parameter ( kx= 9 )
#endif
#ifdef k13
      parameter ( kx=13 )
#endif
#ifdef k14
      parameter ( kx=14 )
#endif
#ifdef k20
      parameter ( kx=20 )
#endif
#ifdef k30
      parameter ( kx=30 )
#endif
#ifdef dry
      parameter (lmx=3 )
#endif
#ifdef moist
      parameter (lmx=4 )
#endif
#ifdef dirnal
      parameter ( ltm = 24 )
#endif
      parameter ( ixkx = ix*kx )
      parameter ( iclt = 100 )
      parameter ( nhem = 2 )
      parameter ( isc = 1 )
      parameter ( il = nhem*iy )
      parameter ( ilp = il + 1 , ilm = il - 1 )
      parameter ( jlp = ilp-(2-nhem)*il )
      parameter ( jtrun1 = jtrun + 1 )
      parameter ( maxwve = 3 )
      parameter ( ixm = ix - 1, ixp = ix + 1, ixp2 = ix + 2 )
      parameter ( kxp2 = kx + 2 )
      parameter ( jin = 3 - nhem, ji = jin + 1 )
      parameter ( mnth = 12 )
c
c   iclt =      length of control block in words
c  jtrun =      the # of zonally truncated waves kept
c  mtrun =      the # of meridionally truncated waves kept
c     ix =      the zonal grid dimension
c     iy =      meridional grid dimension ( one hemisphere )
c     kx =      # of vertical levels
c    isc =      the sector index
c     il =      total meridional grid dimension ( both hemispheres )
c    lmx =      # of vertically diffused variables
c   nhem =      # of hemispheres
c     

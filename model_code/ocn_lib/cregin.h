c====================== include file "cregin.h" ========================
c
c     variables used for computing regional tracer averages (see 
c     "reg1st.F" & "region.F") and for computing term balances for 
c     tracer and momentum equations (see "clinic.F", "tracer.F"
c      and "diag2.F")
c
c
c     mskhr = mask field defining regions in the horizontal
c              (ie: mskhr(i,j) = n indicates point (i,j) is in the
c               "nth" horizontal region   where n=1..nhreg)
c               The "mskhr" masks are used in "region.F" when
c               computing volume weighted tracer averages and in 
c               "clinic.F", "tracer.F" and "diag2.F" when computing
c               term balances for tracers and momentum.
c     mskvr  =  mask field defining regions in the vertical
c              (ie: mskvr(k) = m indicates all points under a horizontal
c               mask at level "k" are in the "mth" vertical region 
c               where m=1..nvreg)     
c               The "mskvr" masks are used in "diag.F", but not
c               in "region.F", where tracer averages are calculated
c               for each k-level.
c
c     hregnm = horizontal region name
c     vregnm = vertical region name
c     volbt  = total volume under a given horizontal region
c     volbk  = volume contained in a horizontal region at level "k"
c     volgt  = total ocean volume
c     volgk  = total ocean volume at level "k"
c     areab  = total ocean surface area for a given horizontal region 
c     areag  = total ocean surface area
c
c     volt   = ocean volume within a particular horizontal & vertical 
c              region (on the "t" grid) for tracer term balances
c     areat  = horizontal ocean area corresponding to "volt"
c     volu   = ocean volume within a particular horizontal & vertical 
c              region (on the "u" grid") for momentum term balances
c     areau  = horizontal ocean area corresponding to "volu"
c     llvreg = level limits for defining vertical regions in term 
c              balance calculations (not used in computing volume 
c              weighted tracer averages)
c              (ie: llvreg(3,1) = 4... means that starting level for 
c                  the third region in the vertical is 4. similarly,
c                  llvreg(3,2) = 6 means the ending level is 6 for that
c                  region. note regions should not overlap.)
c     sumbk  = volume weighted tracer sum for a region at a k level
c     sumbf  = area weighted tracer flux sum for a given region
c     avgbf  = average area weighted tracer flux for a region
c  sumbt = volume weighted tracer sum for a given region
c  avgbt = volume weighted tracer average for a given region
c  sumbk = volume weighted tracer sum for a region at a k level
c  avgbk = volume weighted tracer average for a region at a k level
c  sumgt = total global volume wtd tracer sum
c  avggt = total global volume wtd tracer average
c  sumgk = total global volume wtd tracer sum at a k level
c  avggk = total global volume wtd tracer average at a k level
c  sumbf = area weighted tracer flux sum for a given region
c  avgbf = average area weighted tracer flux for a region
c  sumgf = global area weighted tracer flux sum
c  avggf = global average area weighted tracer flux
c  iceflg  = (false,true) =  when the ifdef seaice option is
c                            (off, on); written on unit 87 in region.F.
#ifdef seaice
c  sumbi = area weighted sea ice sum for a given region
c  avgbi = average area weighted sea ice for a region
c  sumgi = global area weighted sea ice sum
c  avggi = global average area weighted sea ice
#endif
c
      logical  iceflg
c
      character*40 hregnm
      character*20 vregnm
c
      common /cregni/ mskhr(imt,jmt), mskvr(km), llvreg(numreg,2)
      common /cregnr/ volbk(nhreg,km), volbt(nhreg), volgk(km)
     $,               volgt, areab(nhreg), areag, volt(0:numreg)
     $,               volu(0:numreg), areat(0:numreg), areau(0:numreg)
     $,               avgbk(nhreg,km,nt), avggk(km,nt), avggt(nt)
     $,               avgbt(nhreg,nt), avgbf(nhreg,nt), avggf(nt)
     $,               sumbk(nhreg,km,nt), sumbf(nhreg,nt)
     $,               sumbt(nhreg,nt), sumgk(km,nt), sumgt(nt)
     $,               sumgf(nt)
#ifdef seaice
      common /cregnr/ sumbi(nhreg), avgbi(nhreg), avggi, sumgi
#endif
#ifdef cpldleg3
ccccMEDtest
ccccc common /cregnr/ delmed(nt), dtrdt(km,nt), areamed
ccccMEDtest
#endif
      common /cregnc/ hregnm(nhreg), vregnm(nvreg)
      common /cregnl/ iceflg
c
c====================== end of  file "cregin.h" ========================

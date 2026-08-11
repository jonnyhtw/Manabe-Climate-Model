#ifdef prdcld
cSGI 
      integer*8 kflags, iflg, ipower, ishift
cSGI
      common / cloudp / kflags(ix,il), iflg  (ix,kx), cldjrw(il)
      common / cloudp / ipower(kx)   , ishift(kx)
# ifdef varhy
      common / cloudp / hyk(kx)
# endif
# ifdef covcld
      common / cloudp / crx
# endif
#endif
c
c  hyk is the critical humidity value for cloud/no cloud at each level.

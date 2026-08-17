c====================== include file "index.h" =========================
c
c    various starting & ending indices for controlling where quantities
c    are calculated:
c
c    isz   = starting longitudinal index for stream function calculation
c    iez   = ending longitudinal index for stream function calculation
c    ist   = starting longitudinal index for tracer calculations
c    iet   = ending longitudinal index for tracer calculations
c    isu   = starting longitudinal index for velocity calculations
c    ieu   = ending longitudinal index for velocity calculations
c
c    alonis= longiude of arbitrary point within island. 
c    alatis= latitude of arbitrary point within island.
c            coordinate "alonis", "alatis" must point to a land "t" 
c            grid box within the island.
c    nippts= number of island perimeter points
c    iperm = "i" coordinate for the island perimeter point
c    jperm = "j" coordinate for the island perimeter point
c    iofs  = offset for indexing into the island perimeter coordinates
c
c    istf  = starting index for filtering "t" grid points
c    ietf  = ending index for filtering "t" grid points
c    isuf  = starting index for filtering "u,v" grid points
c    ieuf  = ending index for filtering "u,v" grid points
c    iszf  = starting index for filtering "vorticity" grid points
c    iezf  = ending index for filtering "vorticity" grid points
c
      common /index0/ ist(jmt,lseg), iet(jmt,lseg)
     $,               isu(jmt,lseg), ieu(jmt,lseg)
#ifdef rigidlid
      common /index1/ isz(jmt,lseg), iez(jmt,lseg)
# ifdef islands
#  ifdef oldcongrad
congrad
c    aisi  = starting longitudinal for island (used to calculate "isis")
c    aiei  = ending longitudinal for island   (used to calculate "ieis")
c    ajsi  = starting latitudinal for island  (used to calculate "jsis")
c    ajei  = ending latitudinal for island    (used to calculate "jeis")
c    isis  = starting longitudinal model index for island
c    ieis  = ending longitudinal model index for island
c    jsis  = starting latitudinal model index for island
c    jeis  = ending latitudinal model index for island
c
      common /index2/ isis(nisle), ieis(nisle), jsis(nisle), jeis(nisle)
      common /index3/ aisi(nisle), aiei(nisle), ajsi(nisle), ajei(nisle)
#  endif
      common /index2/ alonis(nisle), alatis(nisle), nippts(nisle)
     $,               iofs(nisle), iperm(maxipp), jperm(maxipp) 
# endif
#endif
#if defined fourfil || defined firfil
c
c    bring in filter parameters
c
#include "pfil.h"
      common /index4/ istf(jmtfil,lsegf,km), ietf(jmtfil,lsegf,km)
     $,               isuf(jmtfil,lsegf,km), ieuf(jmtfil,lsegf,km)
     $,               iszf(jmtfil,lsegf),    iezf(jmtfil,lsegf)
#endif
c

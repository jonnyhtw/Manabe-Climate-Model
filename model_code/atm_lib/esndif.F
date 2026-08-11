      subroutine esndif(tmp,es,esdif,ixxx)
c
c  this subroutine is called by tsc or tscs and by vertdf or vrtdfb.
c
#include "parm.h"
c
c  esndif returns saturation vapor pressure and its derivative given
c   temperatures in celsius between -153 deg. c to 101 deg. c.
c
c   subroutine arguments are as follows:
c        tmp    = input temperatures in degrees celsius
c        es     = output saturation vapor pressure in microbars
c       esdif   = output derivative of es in microbars per deg. c
c    iabs(ixxx) = number of input temperatures (maximum = ix*kx)
c    sign of ixxx gives options.  if ixxx is negative, only es returned.
c  odd numbered entries in 5102-element vector table in common block,
c   /establ/, are es values and even numbered ones are esdif for every
c   tenth of a degree from -153 deg. c to 102 deg. c.
c
      dimension tmp(1), es(1), esdif(1)
      dimension t(ixkx), it(ixkx), f(ixkx), esd(ixkx), esm(ixkx)
      dimension esp(ixkx), esc(ixkx)
c
#include "establ.h"
c
      equivalence (esm(1),t(1))
c
c  determine vector length ixx from input parameter ixxx
c
      ixx = iabs(ixxx)
c
c  compute temperature above base (-153c) in tenths of degree c
c
      do 10 i=1,ixx
       t(i) = 10.0 * (tmp(i)+153.0)
   10 continue
c
c  determine integer value of nearest whole tenth of a degree
c
      do 20 i=1,ixx
       it(i) = t(i) + 0.5
   20 continue
c
c  determine algebraic fraction of degree = (t minus it)/10
c
      do 30 i=1,ixx
       f(i) = 0.1 * (t(i)-float(it(i)))
   30 continue
c
c  determine index for fetching esc from table given whole tenth of a
c   degree in vector it (now overwritten by this result) making sure
c   no value is fetched outside the table in subsequent code.
c
      do 40 i=1,ixx
       it(i) = 2*it(i) + 1
       it(i) = max(it(i),1)
       it(i) = min(it(i),5102-3)
   40 continue
c
c  table look-up to fetch esc given index in vector (it)
c
      do 50 i=1,ixx
       esc(i) = etabl(it(i))
   50 continue
c
c  increment index it by one and fetch es derivative, esd
c
      do 60 i=1,ixx
       it(i) = it(i) + 1
   60 continue
      do 70 i=1,ixx
       esd(i) = etabl(it(i))
   70 continue
c
c  apply linear equation to obtain es corresponding to temperature=t
c
      do 80 i=1,ixx
       es(i) = esc(i) + f(i)*esd(i)
   80 continue
c
c  return without computing esdif if ixxx is negative
c
      if (ixxx .lt. 0) return
c
c  increment index it by one and fetch es at (t + 0.1) deg. c, esp
c
      do 90 i=1,ixx
       esp(i) = etabl(it(i)+1)
   90 continue
c
c  decrement index it by 4 and fetch es at (t - 0.1) deg. c, esm,
c    making sure value of index does not drop below one
c
      do 100 i=1,ixx
       it(i) = it(i) - 4
       it(i) = max(it(i),1)
  100 continue
      do 110 i=1,ixx
       esm(i) = etabl(it(i)+1)
  110 continue
c
c  apply linear equation to obtain esdif corresponding to temperature=t
c
      do 120 i=1,ixx
       esdif(i) = esd(i) + 100.0*f(i)*(esp(i)+esm(i) - 2.0*esc(i))
  120 continue
c
c  note: 100*(esp-2*esc+esm)=second derivative of es or d(esd)/dt
c
      return
      end

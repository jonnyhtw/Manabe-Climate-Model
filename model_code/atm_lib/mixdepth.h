c
c mixdepth.h file
c
c the first five level depths must be the same as dzt in ocean model
c (a modified version of thick.h in ocean model) 
c
c     specifications for box thickness (cm) in the vertical
c
c the default is depths for km12 model
c
#if defined km12 || defined icecpld
c
c     (values taken from cyber 205 12-level runs)
c
      data dzti   / 50.9E2,  68.4E2, 100.4E2, 151.1E2, 224.0E2 /
c This was computed using an atan equation and then modified by
c  hand.
#endif
#ifdef km18
c  In this 18-level setup, more of the extra resolution
c  is spent in the upper ocean... the bottom
c  of the 10th level is at exactly 1000 meters.
c
c     (values taken from /net/kd/gauss/setdepths/depths program)
c
c     specifications for box thickness (cm) in the vertical
c     (see kd mail "another 18 levels" in ~/mail/MOM, 6/22/92)
      data dzti/ 32.00e2,   36.98e2,   44.26e2,   54.62e2,   68.94e2 /
c
c     dzt = gaussian depth thicknesses computed as
c     dzt(k)=0.700684E+03exp(-((k-1)-19)**2/0.824360E+02)+0.232162E+02
c     (beta =0.723900E+03  nchop=  2)
c
#endif
#ifdef km25
c
c     (values taken /net/kd/gauss/setdepths/depths program)
c
      data dzti
     $ /    25.00e2,    27.28e2,    30.23e2,    33.98e2,    38.74e2 /
c
c     deltaz = dz = gaussian depth thicknesses computed as
c     dz(k)=0.992964E+03exp(-((k-1)-34)**2/0.232117E+03)+0.181760E+02
#endif
#ifdef km15
c
      data dzti/   30.00e2,   46.15e2,   68.93e2,   99.93e2,  140.63e2 /
c
c     dzt = model level thicknesses computed as
c     dzt(km+1-k)=exp(-(k/sigma)**2) + offset) for k=1,km
c      km =  15
c      sigma=  8.204686
c      offset= 0.0000000
c      1st level thickness =  30.00000
c      # of levels with this thickness =  1
c      depth at bottom of last level =  5700.000
c
#endif
c

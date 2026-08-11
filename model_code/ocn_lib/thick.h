c====================== include file "thick.h" ========================
c
c     specifications for box thickness (cm) in the vertical
#ifdef km12
c
c     (values taken from cyber 205 12-level runs)
c
      data dzt    / 50.9E2,  68.4E2, 100.4E2, 151.1E2, 224.0E2,
     $             319.6E2, 432.5E2, 551.0E2, 660.9E2, 751.9E2,
     $             820.4E2, 868.9E2 /
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
# ifdef old_dzt
c     specifications for box thickness (cm) in the vertical
c     (see kd mail "another 18 levels" in ~/mail/MOM, 6/22/92)
      data dzt / 32.00e2,   36.98e2,   44.26e2,   54.62e2,   68.94e2
     $       ,   88.22e2,  113.40e2,  145.38e2,  184.68e2,  231.52e2
     $       ,  285.50e2,  345.60e2,  409.92e2,  475.96e2,  540.60e2
     $       ,  600.30e2,  651.42e2,  690.70e2
     $ / 
c
c     dzt = gaussian depth thicknesses computed as
c     dzt(k)=0.700684E+03exp(-((k-1)-19)**2/0.824360E+02)+0.232162E+02
c     (beta =0.723900E+03  nchop=  2)
# else
c     specifications for box thickness (cm) in the vertical
c     (see kd mail "another 18 levels" in ~/mail/MOM, 7/31/95)
      data dzt / 40.00e2,   43.42e2,   48.80e2,   56.97e2,   68.98e2
     $       ,   86.02e2,  109.34e2,  140.12e2,  179.26e2,  227.09e2
     $       ,  283.19e2,  346.16e2,  413.53e2,  481.81e2,  546.75e2
     $       ,  603.71e2,  648.23e2,  676.61e2
     $ /
c
c     dzt = gaussian depth thicknesses computed as
c     dzt(k)=0.651276E+03exp(-((k-1)-18)**2/0.663177E+02)+0.350798E+02
c     (beta =0.686356E+03  nchop=  1)
# endif
c
#endif
#ifdef km25
c
c     (values taken /net/kd/gauss/setdepths/depths program)
c
      data dzt
     $ /    25.00e2,    27.28e2,    30.23e2,    33.98e2,    38.74e2
     $ ,    44.69e2,    52.07e2,    61.13e2,    72.14e2,    85.40e2
     $ ,   101.21e2,   119.84e2,   141.59e2,   166.70e2,   195.40e2
     $ ,   227.83e2,   264.06e2,   304.07e2,   347.75e2,   394.84e2
     $ ,   444.97e2,   497.61e2,   552.13e2,   607.75e2,   663.59e2
     $ /
c
c     deltaz = dz = gaussian depth thicknesses computed as
c     dz(k)=0.992964E+03exp(-((k-1)-34)**2/0.232117E+03)+0.181760E+02
#endif
#ifdef km15
c
      data dzt /   30.00e2,   46.15e2,   68.93e2,   99.93e2,  140.63e2
     $      ,  192.11e2,  254.76e2,  327.95e2,  409.81e2,  497.11e2
     $      ,  585.36e2,  669.09e2,  742.41e2,  799.65e2,  836.10e2
     $   /
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

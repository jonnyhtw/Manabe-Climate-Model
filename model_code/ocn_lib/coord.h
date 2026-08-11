c====================== include file "coord.h" =========================
c
c                    model grid point coordinates
c
c     grid definition:
c
c     the model uses a staggard arakawa "b" grid. horizontally, tracer
c     quantities are defined at the centers of "t" grid boxes and
c     velocities are defined at the centers of "u,v" grid boxes.
c     the centers of "u,v" grid boxes are located at the northeast
c     corners of "t" grid boxes. the 1st "t" grid box is located in
c     the southwest corner of the "t" grid. this grid system is
c     replicated and stacked vertically one on top of another from the
c     surface of the ocean downward. vertically, tracers & velocities
c     are defined at the centers of their corresponding boxes and are
c     at the same depths.
c
c
c     inputs:  (data initialized in "blkdta.F")
c
c     stlon = starting longitude of the "t" grid in degrees. It defines
c             the west face of the 1st "t" grid box.
c     stlat = starting latitude of the "t" grid in degrees. It defines
c             the south face of the 1st "t" grid  box
c
c
c
c     in the longitudinal direction, the "t" grid may be divided into
c     one or more sub domains with the following characteristics:
c
c     nxres = number of sub domains of either constant or smoothly
c             varying grid box width (set in "param.h")
c     xmax  = max grid box width on side of the sub domain (degrees)
c     xmin  = min grid box width on side of the sub domain (degrees)
c     xwid  = longitudinal width of the sub domain (degrees)
c     idir  = (1,-1) implies that the longitudinal box width
c             (increases, decreases) going from from west to east
c
c     note:   2*xwid/(xmax+xmin) is the number of boxes which exactly
c             fit within the sub domain. this number must be an integer.
c
c     example:  parameter (imt=140, nxres=2)
c               data xmax /0.5, 1.0/
c               data xmin /0.5, 0.5/
c               data xwid /50.0, 30.0/
c               data idir /1, 1/
c               data stlon /10.0/
c
c     describes a grid with 2 distinct sub domains of resolution. the
c     1st sub domain contains 100 "t" grid boxes
c     { 100 = xwid(1)/(0.5*(xmax(1)+xmin(1))) } between longitude
c     10.0 deg east { stlon=10.0 } & 60.0 deg east {xwid(1) = 50.0}.
c     the grid box width { delta x = 0.5*(xmax(1)+xmin(1)) } is a
c     constant 0.5 deg. the 1st "t" grid point is at 10.25 deg east
c     and the 100th is at 59.75 deg east. the 2nd sub domain contains
c     40 "t" grid boxes 
c     { 40 = xwid(2)/(0.5*(xmax(2)+xmin(2))) } between longitude
c     60.0 deg east { where the first domain ends } & 90.0 deg east.
c     the box width (delta x) gradually increases from 0.5 deg at
c     60.0 deg east to 1.0 deg at 90 deg east. the variation is
c     fitted to a cosine function. idir = 1 says the grid delta x 
c     increases with increasing longitude. ( 0.5 deg in the west and
c     1.0 deg in the east. idir = -1 would indicate the reverse)
c     Note that 100 + 40 must equal imt.
c
c
c     in the latitudinal direction, the "t" grid may be divided into
c     one or more sub domains with the following characteristics:
c
c     nyres = number of sub domains of either constant or smoothly
c             varying grid box height (set in "param.h")
c     ymax  = max grid box height on side of the sub domain (degrees)
c     ymin  = min grid box height on side of the sub domain (degrees)
c     ywid  = latitudinal height of the sub domain (degrees)
c     jdir  = (1,-1) implies that the latitudinal box height
c             (increases, decreases) going from from south to north
c
c     note:   2*ywid/(ymax+ymin) is the number of boxes which exactly
c             fit within the sub domain. it must be an integer.
c
c     see above example:
c
c     primary grid coordinates:
c
c     xt(i) = longitude of the ith "t" point in degrees.   i=1..imt
c     xu(i) = longitude of the ith "u,v" point in degrees. i=1..imt
c     yt(j) = latitude of the jth "t" point in degrees.    j=1..jmt
c     yu(j) = latitude of the jth "u,v" point in degrees.  j=1..jmt
c     zt(k)   = distance from surface down to center of level k (in cm)
c               (depth of "t" and "u,v" grid points)
c     zw(k)   = distance from surface down to bottom of level k (in cm)
c
c
      common /coord/ stlon, stlat, xt(imt), yt(jmt), xu(imt), yu(jmt)
     $,              zw(km), zt(kmp1)
c
      common /cdmesh/ idir(nxres), jdir(nyres)
     $,               xmax(nxres), xmin(nxres), xwid(nxres)
     $,               ymax(nyres), ymin(nyres), ywid(nyres)
c

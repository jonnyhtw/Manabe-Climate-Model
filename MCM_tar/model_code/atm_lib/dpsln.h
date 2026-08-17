      double precision apsp, apsln, apsdt, apsm
      double precision adtx2
c
      common / dpsln  / apsp, apsln, apsdt, apsm
      common / dpsln  / adtx2
c
c     apsp  is the real part of the tau+1 psln
c     apsm  is the real part of the tau-1 psln
c     apsln is the real part of the tau   psln
c     apsdt is the real part newly computed tendency
c
c     adtx2 is dtx2 in double precision (time step length * 2)
c
c  These variables are in double precision for accuracy
c  Without this the model gained .3mb/200yrs in global mean ps.
c

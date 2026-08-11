      block data s1gwgt
c
#include "parm.h"
#include "parmkp.h"
c
#include "sgwgt1.h"
c
#include "albsea.h"
c
#ifdef old_levels
c
c   load array qmh with the sigma-coordinates of the kx+1 atmospheric
c   half-levels.  The coordinates of the full levels will be computed
c   in parmtr
c
# ifdef k9
      data qmh/ 0.0,.050,.140,.270,.430,.600,.760,.900,.980,1.000/
# endif
c
# ifdef k13
      data qmh/
     1     0.0000000e 00, 0.2701777e-01, 0.7307184e-01, 0.1487764e 00,
     2     0.2646376e 00, 0.4199635e 00, 0.5897341e 00, 0.7369612e 00,
     3     0.8430778e 00, 0.9114627e 00, 0.9530230e 00, 0.9775428e 00,
     4     0.9917893e 00, 0.1000000e 01/
# endif
c
# ifdef k14
      data qmh /0.0000,0.0300,0.0707,0.1311,0.2102,0.3036,0.4062,
     1          0.5138,0.6226,0.7284,0.8255,0.9066,0.9640,0.9933,
     2          1.0000/
# endif
c
# ifdef k20
      data qmh /0.0000,0.0207,0.0442,0.0758,0.1174,0.1686,0.2279,
     1          0.2938,0.3644,0.4381,0.5138,0.5901,0.6656,0.7386,
     2          0.8071,0.8686,0.9202,0.9595,0.9849,0.9970,1.0000/
# endif
c
# ifdef k30
      data qmh /0.0000,0.0139,0.0278,0.0442,0.0642,0.0885,0.1174,
     1          0.1505,0.1875,0.2279,0.2712,0.3169,0.3644,0.4133,
     2          0.4632,0.5138,0.5647,0.6154,0.6656,0.7147,0.7621,
     3          0.8071,0.8490,0.8870,0.9202,0.9479,0.9695,0.9849,
     4          0.9943,0.9986,1.0000/
# endif
#endif
c
c   load array "oas" with zonal means of the ocean albedo for each ten
c   degrees of latitude. for models with two hemispheres, the data are
c   used symmetrically. the sea surface albedoes were derived from two
c   data sources. from 70 degrees to the equator, the data are from
c   posey & clapp. from 70 degrees to the north pole, the data were
c   extrapolated by eye from Budyko table.
c
      data oas / 0.206, 0.161, 0.110, 0.097, 0.089, 0.076, 0.068, 0.063,
     1 3*0.060,  0.063, 0.068, 0.076, 0.089, 0.097, 0.110, 0.161, 0.206/
c
      end

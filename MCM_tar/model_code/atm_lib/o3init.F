      subroutine o3init 
c
c  this subroutine is called from parmtr. it is called only once at the
c  beginning of a run. it does not call any subroutine. it can be run
c  independently with a driver by loading arrays q and qmh with their
c  appropriate values, declaring "seasnl" a logical variable and
c  assigning it a value of ".true." for seasonal values of ozone to be
c  calculated or ".false." for annual mean values to be calculated, and
c  adding the following excutable statements :
c             call o3init
c             stop
c
c  this subroutine computes ozone for annual mean models or seasonal
c  models. the observed ozone data in arrays ro31, spring, fall,
c  winter, and summer are for the northern hemisphere only. the units
c  are micrograms/cm**2 and the data are given for every 10 degrees of
c  latitude for 41 linear levels for every 2 kilometers from 0 to 80
c  kilometers. seasonal data are given only for the bottom 16 levels.
c  the data were obtained by richard wetherald from jerry mahlman.
c
c  the units of the resultant arrays dduo3n, xo3sf, and xo3ws are
c  grams/gram. dduo3n contains annual mean data, xo3sf contains spring
c  and fall data, and xo3ws contains winter and summer data.
c
#include "parm.h"
#include "parmncc.h"
#include "parmkp.h"
      parameter ( nk=81 )
      parameter ( nkk=41 )
      parameter ( nkp=nk+1 )
      parameter ( nkkm=nkk-1 )
      parameter ( nkm=nk-1 )
      parameter ( nkmm=nk-3 )
c
      character*16 seasns(2)
c
#include "ctlblk.h"
c
c  by setting the logical variables aprto3, bprto3, and cprto3 to true,
c  one can verify that the input and results are correct. aprt03 is for
c  printing the input arrays and the intermediate result of array rdata
c  which contains ozone data for a 1 kilometer resolution obtained by
c  means of the bessel interpolation formula. bprto3 can be used to
c  print the total integrated values of the observed and resultant
c  arrays. cprto3 can be used to print the resultant arrays dduo3n or
c  xo3sf and xo3ws. if pcho3 is true, the resultant array(s) will be
c  punched.
c
c  variable seasnl is true for seasonal models and false for annual mean
c  models.
c
      logical aprto3, bprto3, cprto3, pcho3
c
      dimension annmn1(190)  , annmn2(190), annmn3(30)
      dimension dump1 (nkkm) , dump2 (nkk)
      dimension dumph1(nkk)  , dumph2(nkk)
      dimension latn(10)
cSGI
      real*8 spring8(10,16), summer8(10,16)
      real*8 fall8  (10,16), winter8(10,16)
cSGI
c
#include "latdat.h"
c
#include "sgwgt1.h"
c
      common / pdata  /  p(nk), ph(nkp)
c
      common / ro3dat /  ro31(10,nkk)
c
      common / season /  spring(10,16), summer(10,16)
      common / season /  fall  (10,16), winter(10,16)
c
      dimension  lat  (19)
      dimension  duo3n(19,nkk), rdata(nk), ro32(10,nkk)
      dimension  rstd (   nkk)
      dimension  phalf(kp),     pstd(kx) , rbar(kx)
c
      equivalence ( annmn1(1) , ro31(1, 1) )
      equivalence ( annmn2(1) , ro31(1,20) )
      equivalence ( annmn3(1) , ro31(1,39) )
c
      equivalence ( dumph1(1) , ph(1) ) , ( dumph2(1) , ph(42) )
      equivalence ( dump1 (1) , p (1) ) , ( dump2 (1) , p (41) )
c
      equivalence ( latn  (1) , lat(10) )
c
      data aprto3 / .false. /
      data bprto3 / .false. /
      data cprto3 / .false. /
cpunchdata pcho3  / .false. /
c
      data seasns / ' winter & summer' , '  spring & fall' /
c
      data half   /  0.50e+00 /
      data oconst /  0.98e+03 /
      data quartr /  0.25e+00 /
      data sxteen /  0.16e+02 /
      data zero   /  0.00e+00 /
c
c  load array ph with the pressure values of the half-levels which will
c  be used in the Bessel interpolation.
c
      data dumph1 /       0.0000000e+00,
     $     0.1027246e-04, 0.1239831e-04, 0.1491845e-04, 0.1788053e-04,
     $     0.2135032e-04, 0.2540162e-04, 0.3011718e-04, 0.3558949e-04,
     $     0.4192172e-04, 0.4922875e-04, 0.5763817e-04, 0.6729146e-04,
     $     0.7834518e-04, 0.9097232e-04, 0.1053635e-03, 0.1217288e-03,
     $     0.1402989e-03, 0.1613270e-03, 0.1850904e-03, 0.2119495e-03,
     $     0.2423836e-03, 0.2768980e-03, 0.3160017e-03, 0.3602623e-03,
     $     0.4103126e-03, 0.4668569e-03, 0.5306792e-03, 0.6026516e-03,
     $     0.6839018e-03, 0.7759249e-03, 0.8803303e-03, 0.9987843e-03,
     $     0.1133178e-02, 0.1285955e-02, 0.1460360e-02, 0.1660001e-02,
     $     0.1888764e-02, 0.2151165e-02, 0.2452466e-02, 0.2798806e-02/
c
      data dumph2/
     $     0.3197345e-02, 0.3656456e-02, 0.4185934e-02, 0.4797257e-02,
     $     0.5503893e-02, 0.6321654e-02, 0.7269144e-02, 0.8368272e-02,
     $     0.9644873e-02, 0.1112946e-01, 0.1285810e-01, 0.1487354e-01,
     $     0.1722643e-01, 0.1997696e-01, 0.2319670e-01, 0.2697093e-01,
     $     0.3140135e-01, 0.3660952e-01, 0.4274090e-01, 0.4996992e-01,
     $     0.5848471e-01, 0.6847525e-01, 0.8017242e-01, 0.9386772e-01,
     $     0.1099026e 00, 0.1286765e 00, 0.1506574e 00, 0.1763932e 00,
     $     0.2065253e 00, 0.2415209e 00, 0.2814823e 00, 0.3266369e 00,
     $     0.3774861e 00, 0.4345638e 00, 0.4984375e 00, 0.5697097e 00,
     $     0.6490189e 00, 0.7370409e 00, 0.8344896e 00, 0.9421190e 00,
     $     0.1000000e 01/
c
c  load array p with the pressure values of the full-levels which will
c  be used in the Bessel interpolation.
c
      data dump1/
     $     0.9300000e-05, 0.1129521e-04, 0.1360915e-04, 0.1635370e-04,
     $     0.1954990e-04, 0.2331653e-04, 0.2767314e-04, 0.3277707e-04,
     $     0.3864321e-04, 0.4547839e-04, 0.5328839e-04, 0.6234301e-04,
     $     0.7263268e-04, 0.8450696e-04, 0.9793231e-04, 0.1133587e-03,
     $     0.1307170e-03, 0.1505832e-03, 0.1728373e-03, 0.1982122e-03,
     $     0.2266389e-03, 0.2592220e-03, 0.2957792e-03, 0.3376068e-03,
     $     0.3844381e-03, 0.4379281e-03, 0.4976965e-03, 0.5658476e-03,
     $     0.6418494e-03, 0.7287094e-03, 0.8261995e-03, 0.9380076e-03,
     $     0.1063498e-02, 0.1207423e-02, 0.1369594e-02, 0.1557141e-02,
     $     0.1769657e-02, 0.2015887e-02, 0.2295520e-02, 0.2620143e-02/
c
      data dump2/
     $     0.2989651e-02, 0.3419469e-02, 0.3909867e-02, 0.4481491e-02,
     $     0.5135272e-02, 0.5898971e-02, 0.6774619e-02, 0.7799763e-02,
     $     0.8978218e-02, 0.1036103e-01, 0.1195488e-01, 0.1382957e-01,
     $     0.1599631e-01, 0.1855114e-01, 0.2151235e-01, 0.2501293e-01,
     $     0.2908220e-01, 0.3390544e-01, 0.3952926e-01, 0.4621349e-01,
     $     0.5403168e-01, 0.6330472e-01, 0.7406807e-01, 0.8677983e-01,
     $     0.1015345e 00, 0.1189603e 00, 0.1391863e 00, 0.1630739e 00,
     $     0.1908004e 00, 0.2235461e 00, 0.2609410e 00, 0.3036404e 00,
     $     0.3513750e 00, 0.4055375e 00, 0.4656677e 00, 0.5335132e 00,
     $     0.6083618e 00, 0.6923932e 00, 0.7845676e 00, 0.8875882e 00,
     $     0.1000000e 01/
c
c  load array ro31 with the observed annual mean data from equator to
c  pole and from the top of the atmosphere to the bottom.
c
      data annmn1/
     *    0.55, 0.50, 0.45, 0.45, 0.40, 0.35, 0.35, 0.30, 0.30, 0.30,
     *    0.55, 0.51, 0.46, 0.47, 0.42, 0.38, 0.37, 0.36, 0.35, 0.35,
     *    0.55, 0.53, 0.48, 0.49, 0.44, 0.42, 0.41, 0.40, 0.38, 0.38,
     *    0.60, 0.55, 0.52, 0.52, 0.50, 0.47, 0.46, 0.44, 0.42, 0.41,
     *    0.65, 0.60, 0.55, 0.56, 0.53, 0.52, 0.50, 0.48, 0.45, 0.45,
     *    0.75, 0.65, 0.60, 0.60, 0.55, 0.55, 0.55, 0.50, 0.48, 0.47,
     *    0.80, 0.75, 0.75, 0.75, 0.70, 0.70, 0.65, 0.63, 0.60, 0.60,
     *    0.90, 0.85, 0.85, 0.80, 0.80, 0.75, 0.75, 0.74, 0.72, 0.71,
     *    1.10, 1.05, 1.00, 0.90, 0.90, 0.90, 0.85, 0.83, 0.80, 0.80,
     *    1.40, 1.30, 1.25, 1.25, 1.25, 1.20, 1.15, 1.10, 1.05, 1.00,
     *    1.70, 1.70, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.50, 1.50,
     *    2.10, 2.00, 1.90, 1.90, 1.90, 1.80, 1.80, 1.80, 1.70, 1.70,
     *    2.40, 2.30, 2.20, 2.20, 2.20, 2.10, 2.10, 2.10, 2.00, 2.00,
     *    2.70, 2.50, 2.50, 2.50, 2.50, 2.50, 2.40, 2.40, 2.30, 2.30,
     *    2.90, 2.80, 2.70, 2.70, 2.70, 2.70, 2.70, 2.70, 2.60, 2.60,
     *    3.10, 3.10, 3.00, 3.00, 3.00, 3.00, 3.00, 3.00, 2.90, 2.80,
     *    3.30, 3.40, 3.40, 3.60, 3.70, 3.90, 4.00, 4.10, 4.00, 3.80,
     *    3.60, 3.80, 3.90, 4.20, 4.70, 5.30, 5.60, 5.70, 5.50, 5.20,
     *    4.10, 4.30, 4.70, 5.20, 6.00, 6.70, 7.00, 6.80, 6.40, 6.20/
c
      data annmn2/
     *    5.40, 5.70, 6.00, 6.60, 7.30, 8.00, 8.40, 7.70, 7.10, 6.70,
     *    6.70, 6.80, 7.00, 7.60, 8.30,10.00, 9.60, 8.20, 7.50, 7.20,
     *    9.20, 9.30, 9.40, 9.60,10.30,10.60,10.00, 8.50, 7.70, 7.30,
     *   12.60,12.10,12.00,12.10,11.70,11.00,10.00, 8.60, 7.80, 7.40,
     *   14.20,13.50,13.10,12.80,11.90,10.90, 9.80, 8.50, 7.80, 7.50,
     *   14.30,14.00,13.40,12.70,11.60,10.60, 9.30, 8.40, 7.60, 7.30,
     *   14.60,14.00,13.40,12.10,11.20,10.20, 8.90, 8.10, 7.50, 7.10,
     *   13.40,13.20,12.60,11.50,10.50, 9.60, 8.30, 7.60, 7.10, 6.80,
     *   11.00,10.90,10.80,10.30, 9.70, 9.20, 7.90, 7.20, 6.70, 6.50,
     *    7.20, 7.50, 7.80, 8.10, 8.10, 7.90, 7.10, 6.80, 6.50, 6.30,
     *    4.10, 4.40, 5.10, 5.90, 6.30, 6.80, 6.60, 6.40, 6.20, 6.00,
     *    1.80, 2.00, 2.40, 3.30, 4.00, 5.00, 5.50, 5.60, 5.80, 5.60,
     *    0.50, 0.58, 0.83, 1.40, 2.40, 3.20, 4.00, 4.20, 4.60, 4.70,
     *    0.10, 0.14, 0.24, 0.29, 1.00, 1.80, 2.60, 2.80, 3.10, 3.20,
     *    0.07, 0.09, 0.13, 0.25, 0.55, 1.00, 1.50, 1.70, 2.10, 2.30,
     *    0.06, 0.07, 0.08, 0.14, 0.35, 0.55, 0.98, 1.10, 1.20, 1.30,
     *    0.05, 0.05, 0.06, 0.09, 0.16, 0.26, 0.48, 0.55, 0.66, 0.68,
     *    0.05, 0.05, 0.06, 0.07, 0.08, 0.10, 0.14, 0.16, 0.20, 0.25,
     *    0.04, 0.04, 0.05, 0.06, 0.06, 0.07, 0.08, 0.08, 0.08, 0.09/
c
      data annmn3/
     *    0.04, 0.04, 0.04, 0.05, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06,
     *    0.03, 0.03, 0.04, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
     *    0.02, 0.02, 0.03, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04/
c
      data spring8/
     *   14.90,14.20,13.30,12.50,11.20,10.30, 9.50, 8.60, 7.50, 7.40,
     *   14.50,14.10,13.00,11.80,10.50, 9.80, 9.20, 7.90, 7.40, 7.40,
     *   11.80,11.50,10.90,10.50, 9.90, 9.60, 8.90, 7.50, 7.20, 7.20,
     *    7.30, 7.70, 7.80, 8.40, 8.40, 8.50, 7.90, 7.40, 7.10, 7.10,
     *    4.10, 4.40, 5.30, 6.60, 6.90, 7.50, 7.40, 7.20, 7.00, 6.90,
     *    1.80, 1.90, 2.50, 3.30, 4.50, 5.80, 6.30, 6.30, 6.40, 6.10,
     *    0.40, 0.50, 0.80, 1.20, 2.70, 3.60, 4.60, 4.70, 5.00, 5.20,
     *    0.10, 0.15, 0.20, 0.50, 1.40, 2.10, 3.00, 3.20, 3.50, 3.90,
     *    0.07, 0.10, 0.12, 0.30, 1.00, 1.40, 1.80, 1.90, 2.30, 2.50,
     *    0.06, 0.08, 0.10, 0.15, 0.60, 0.80, 1.40, 1.50, 1.50, 1.60,
     *    0.05, 0.05, 0.06, 0.09, 0.20, 0.40, 0.70, 0.80, 0.90, 0.90,
     *    0.05, 0.05, 0.06, 0.08, 0.10, 0.13, 0.20, 0.25, 0.30, 0.40,
     *    0.05, 0.05, 0.05, 0.06, 0.07, 0.07, 0.08, 0.09, 0.10, 0.13,
     *    0.05, 0.05, 0.05, 0.05, 0.06, 0.06, 0.06, 0.06, 0.07, 0.07,
     *    0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.06, 0.06, 0.06,
     *    0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.05, 0.05, 0.05/
c
      data summer8/
     *   14.20,13.80,13.20,12.50,11.70,10.50, 8.60, 7.80, 7.50, 6.60,
     *   12.50,12.40,12.20,11.70,10.80, 9.80, 7.80, 7.20, 6.50, 6.10,
     *   10.60,10.50,10.40,10.10, 9.60, 9.00, 7.10, 6.80, 6.10, 5.90,
     *    7.00, 7.40, 7.90, 7.80, 7.60, 7.30, 6.20, 6.10, 5.80, 5.60,
     *    4.20, 4.60, 5.10, 5.60, 5.90, 5.90, 5.90, 5.80, 5.60, 5.30,
     *    2.10, 2.30, 2.60, 2.90, 3.50, 4.30, 4.80, 4.90, 5.10, 5.10,
     *    0.70, 0.80, 1.00, 1.50, 2.00, 2.80, 3.50, 3.60, 3.70, 4.00,
     *    0.15, 0.20, 0.40, 0.50, 0.60, 1.40, 2.10, 2.20, 2.30, 2.50,
     *    0.08, 0.10, 0.15, 0.25, 0.30, 0.90, 1.20, 1.30, 1.40, 1.60,
     *    0.07, 0.08, 0.10, 0.14, 0.20, 0.50, 0.70, 0.90, 0.90, 0.80,
     *    0.05, 0.06, 0.08, 0.12, 0.14, 0.20, 0.35, 0.40, 0.60, 0.50,
     *    0.05, 0.05, 0.08, 0.09, 0.09, 0.09, 0.11, 0.12, 0.15, 0.18,
     *    0.04, 0.05, 0.06, 0.07, 0.07, 0.08, 0.08, 0.08, 0.08, 0.08,
     *    0.04, 0.04, 0.05, 0.07, 0.07, 0.07, 0.07, 0.07, 0.06, 0.05,
     *    0.02, 0.02, 0.04, 0.05, 0.05, 0.05, 0.05, 0.05, 0.04, 0.04,
     *    0.02, 0.02, 0.03, 0.04, 0.04, 0.04, 0.04, 0.04, 0.03, 0.03/
c
      data fall8/
     *   14.80,14.20,13.80,12.20,11.00, 9.80, 8.50, 7.80, 7.40, 6.90,
     *   13.20,13.00,12.50,11.30,10.40, 9.00, 7.80, 7.50, 7.00, 6.60,
     *   10.60,10.60,10.70,10.10, 9.40, 8.60, 7.50, 7.00, 6.50, 6.10,
     *    7.00, 7.30, 7.50, 7.50, 7.50, 7.30, 6.70, 6.40, 6.00, 5.80,
     *    3.80, 4.00, 4.70, 5.00, 5.20, 5.90, 5.80, 5.60, 5.50, 5.50,
     *    1.40, 1.60, 2.40, 3.00, 3.70, 4.10, 4.60, 4.80, 5.10, 5.00,
     *    0.40, 0.50, 0.90, 1.20, 2.00, 2.70, 3.20, 3.60, 4.30, 4.10,
     *    0.07, 0.10, 0.20, 0.30, 0.80, 1.40, 2.10, 2.40, 2.70, 3.00,
     *    0.06, 0.07, 0.09, 0.15, 0.30, 0.70, 1.20, 1.40, 1.60, 2.00,
     *    0.05, 0.05, 0.06, 0.12, 0.15, 0.30, 0.60, 0.70, 0.80, 0.80,
     *    0.04, 0.05, 0.06, 0.08, 0.09, 0.15, 0.30, 0.40, 0.40, 0.40,
     *    0.04, 0.04, 0.05, 0.05, 0.06, 0.09, 0.12, 0.13, 0.15, 0.15,
     *    0.03, 0.03, 0.04, 0.05, 0.05, 0.06, 0.07, 0.07, 0.06, 0.07,
     *    0.03, 0.03, 0.04, 0.05, 0.05, 0.05, 0.06, 0.06, 0.05, 0.05,
     *    0.02, 0.02, 0.03, 0.05, 0.05, 0.05, 0.04, 0.04, 0.04, 0.04,
     *    0.02, 0.02, 0.02, 0.04, 0.04, 0.04, 0.03, 0.03, 0.03, 0.03/
c
      data winter8/
     *   14.50,14.00,13.50,11.30,11.00,10.00, 9.00, 8.30, 7.50, 7.30,
     *   13.50,13.20,12.50,11.10,10.40, 9.70, 8.20, 7.80, 7.40, 6.80,
     *   10.80,10.90,11.00,10.40,10.00, 9.60, 7.90, 7.50, 7.00, 6.70,
     *    7.30, 7.50, 7.80, 8.50, 9.00, 8.50, 7.70, 7.40, 6.90, 6.50,
     *    4.10, 4.50, 5.30, 6.20, 7.30, 7.70, 7.30, 7.00, 6.60, 6.40,
     *    1.80, 2.00, 2.20, 3.80, 4.30, 5.60, 6.20, 6.20, 6.40, 6.20,
     *    0.30, 0.50, 0.60, 1.50, 2.80, 3.70, 4.50, 4.70, 5.50, 5.60,
     *    0.09, 0.10, 0.15, 0.60, 1.20, 2.10, 3.00, 3.50, 4.00, 4.30,
     *    0.06, 0.08, 0.10, 0.30, 0.60, 1.10, 1.90, 2.20, 2.90, 3.00,
     *    0.04, 0.05, 0.06, 0.15, 0.45, 0.60, 1.10, 1.30, 1.60, 1.80,
     *    0.04, 0.04, 0.04, 0.08, 0.20, 0.30, 0.55, 0.60, 0.75, 0.90,
     *    0.04, 0.04, 0.04, 0.05, 0.06, 0.10, 0.12, 0.15, 0.20, 0.25,
     *    0.04, 0.04, 0.03, 0.04, 0.05, 0.06, 0.07, 0.07, 0.07, 0.08,
     *    0.03, 0.03, 0.04, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
     *    0.03, 0.03, 0.03, 0.04, 0.04, 0.04, 0.05, 0.05, 0.04, 0.04,
     *    0.02, 0.02, 0.02, 0.04, 0.04, 0.04, 0.04, 0.04, 0.03, 0.03/
c
 200  format (10e12.5)
 201  format (19f6.2)
 202  format (20f6.2)
 203  format (/9(i5,'s'),i5,1x,9(i5,'n')/)
 205  format ('- ph data ')
 206  format ('- p data')
 207  format ('- phalf data to be computed')
 208  format ('- pstd data to be computed')
 209  format ('- duo3n data computed')
 210  format ('- rdata computed for latitude =',i4/)
 211  format ('-latitude =',i4,'  o3rd =',e12.5,'  o3tot =',e12.5,
     1                         '  o3du =',e12.5)
 212  format ('-',a16,' data scaled by 1.0e+04')
 214  format (10f6.2)
 215  format ('- ro31 data'//10(i5,'n')/)
 216  format ('- ro32 data'//10(i5,'n')/)
 217  format ('- annual mean data')
c
 300  format ((5x,9(1h*,4pf6.5,1h,)))
 301  format (1h1)
c
c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c
cSGI
      spring = spring8
      summer = summer8
      fall   = fall8
      winter = winter8
cSGI
      qp0 = 1.0e+06
      p0  = 1013250.0e+00
c
c  compute the observed latitudes
c
      if ( aprto3 .or. bprto3 )  then
           do 3300 l = 1,19
           lat(l) = 10 * ( l - 10 )
3300       continue
      end if
c
c  convert the model sigma full-levels to dynes/cm**2.
c
      do 3400 k = 1,kx
      pstd(k) = q(k)
      pstd(k) = pstd(k) * qp0
3400  continue
c
c  convert the model sigma half-levels to dynes/cm**2.
c
      do 3500 k = 1,kp
      phalf(k) = qmh(k-1)
      phalf(k) = phalf(k) * qp0
3500  continue
c
c  convert the observed sigma half-levels to dynes/cm**2.
c
      do 3600 k = 1,nkp
      ph(k) = ph(k) * p0
3600  continue
c
c  convert the observed sigma full-levels to dynes/cm**2.
c
      do 3700 k = 1,nk
      p(k) = p(k) * p0
3700  continue
c
c  for seasonal models copy the first 25 levels of ozone from array
c  ro31 to array r032. the annual mean ozone data for the bottom 16
c  levels can be computed from the arrays for the four seasons.
c
      if ( seasnl )  then
           do 3800 j = 1,25
           do 3800 i = 1,10
3800       ro32(i,j) = ro31(i,j)
ctemp else
ctemp      do 3900 j = 26,nkk
ctemp      do 3900 i = 1,10
ctemp      ro31(i,j) = ( winter(i,j-25) +
ctemp*                   summer(i,j-25) +
ctemp*                   spring(i,j-25) +
ctemp*                   fall  (i,j-25) ) * quartr
3900       continue
      end if
c
      if ( aprto3 )  then
           write ( 6 , 205 )
           write ( 6 , 200 )  ph
           write ( 6 , 206 )
           write ( 6 , 200 )  p
           write ( 6 , 207 )
           write ( 6 , 200 )  phalf
           write ( 6 , 208 )
           write ( 6 , 200 )  pstd
      end if
c
c  whenever seasnl is true, the winter & summer ozone values will
c  be computed in the first iteration of the folowing do-loop and
c  the spring & fall values during the second. when seasnl is false,
c  there will be only one iteration of the do-loop to calculate the
c  annual mean ozone values.
c
#ifdef seanal
      iters = 2
#endif
#ifdef annlmn
      if ( .not. seasnl )  iters=1
#endif
c
      do 2500 iseasn = 1,iters
c
c  for seasonal computation fill the bottom 16 levels of ro31 with the
c  winter or spring values and the bottom 16 levels of r032 with summer
c  or fall values. ro32 contains the data for the "opposite" season of
c  that in array ro31.
c
      if ( seasnl )  then
      if ( iseasn .eq. 1 )  then
           do 800 j = 26,nkk
           do 800 i = 1,10
           ro31(i,j) = winter(i,j-25)
           ro32(i,j) = summer(i,j-25)
 800       continue
      else
           do 900 j = 26,nkk
           do 900 i = 1,10
           ro31(i,j) = spring(i,j-25)
           ro32(i,j) = fall  (i,j-25)
 900       continue
      end if
      end if
c
c  move the observed, northern hemisphere ozone values to array duo3n.
c  since no ozone data is available for the southern hemisphere, fill
c  array du03n with the data for the corresponding latitude of the
c  northern hemisphere. for seasonal models, duo3n will have summer
c  (fall) values for the southern hemisphere when it has winter (spring)
c  values for the northern hemisphere. for seasonal models, the data
c  for the equator is an average of the two opposite seasons.
c
      if ( .not. seasnl )  then
           do 1000 kk=1,nkk
           do 1000 l=1,10
           duo3n(l  ,kk) = ro31(11-l,kk)
1000       duo3n(l+9,kk) = ro31(   l,kk)
      else
           do 1200 kk=1,nkk
           do 1100 l=1,10
           duo3n(l  ,kk) = ro31(11-l,kk)
1100       duo3n(l+9,kk) = ro32(   l,kk)
           duo3n( 10,kk) = half * ( ro31(1,kk) + ro32(1,kk) )
1200       continue
      end if
c
      if ( aprto3 )  then
           write ( 6 , 301 )
           write ( 6 , 215 )  latn
           write ( 6 , 214 )  ro31
           write ( 6 , 216 )  latn
           write ( 6 , 214 )  ro32
           write ( 6 , 209 )
           write ( 6 , 203 )  lat
           write ( 6 , 201 )  duo3n
           write ( 6 , 203 )  lat
      end if
c
      do 2200 l=1,19
c
      do 1300 kk=1,nkk
1300  rstd(kk) = duo3n(l,kk)
c
c   use bessels half-point interpolation formula to fill in 1 kilometer
c   resolution.
c
      do 1400 k = 4,nkmm,2
      ki = k / 2
1400  rdata(  k) = ( rstd(ki)   + rstd(ki+1) ) * half
     *           - ( rstd(ki+2) - rstd(ki+1) - rstd(ki) + rstd(ki-1)  )
     *           / sxteen
      rdata(  2) = half * ( rstd(  2) + rstd(   1) )
      rdata(nkm) = half * ( rstd(nkk) + rstd(nkkm) )
c
c     put unchanged data into new array
c
      do 1500 k = 1,nk,2
      kq = ( k + 1 ) / 2
1500  rdata(k) = rstd(kq)
c
      if ( aprto3 )  then
           write ( 6 , 210 )  lat(l)
           write ( 6 , 202 )  rdata
      end if
c
c     calculate layer-mean ozone mixing ratio for each model level
c
      do 1800 kk = 1,kx
      rbar(kk) = zero
c
c     loop to calculate sums to get layer ozone mean
c
      do 1600 k = 1,nk
      if ( ph(k+1) .lt. phalf(kk)   )  go to 1600
      if ( ph(k)   .gt. phalf(kk+1) )  go to 1600
      if ( ph(k+1) .lt. phalf(kk+1) .and. ph(k) .lt. phalf(kk) )
     * rbar(kk) = rbar(kk) + rdata(k) * ( ph(k+1) - phalf(kk) )
      if ( ph(k+1) .lt. phalf(kk+1) .and. ph(k) .ge. phalf(kk) )
     * rbar(kk) = rbar(kk) + rdata(k) * ( ph(k+1) - ph(k) )
      if ( ph(k+1) .gt. phalf(kk+1) .and. ph(k) .gt. phalf(kk) )
     * rbar(kk) = rbar(kk) + rdata(k) * ( phalf(kk+1) - ph(k) )
1600  continue
      rbar(kk) = rbar(kk) / ( phalf(kk+1) - phalf(kk) )
c
c     code to cover case when model resolution is so fine that no value
c     of p(k) in the ozone data array falls betwen phalf(kk+1) and
c     phalf(kk).   procedure is to simply grab the nearest value from
c     rdata
c
      if ( rbar(kk) .le. zero )  then
           do 1700 k = 1,nk
           if ( ph(k) .lt. phalf(kk) .and. ph(k+1) .ge. phalf(kk+1) )
     *          rbar(kk) = rdata(k)
1700       continue
      end if
c
1800  continue
c
c     calculate total ozone
c
      if ( bprto3 )  then
           o3rd = zero
           do 1900 kk = 1,nkm
1900       o3rd = o3rd + rdata(kk) * ( ph(kk+1) - ph(kk) )
           o3rd = o3rd + rdata(nk) * ( p (nk)   - ph(nk) )
           o3rd = o3rd / oconst
c
           o3tot = zero
           do 2000 kk = 1,kx
2000       o3tot = o3tot + rbar(kk) * ( phalf(kk+1) - phalf(kk) )
           o3tot = o3tot / oconst
c
c     o3tot units are micrograms/cm**2
c
           o3du = o3tot / 2.144e+00
c
c     o3du units are dobson units (10**-3 atm-cm)
c
           write ( 6 , 211 )  lat(l), o3rd, o3tot, o3du
      end if
c
c  convert the resultant ozone from micrograms/cm**2 to grams/gram.
c
      do 2100 kk = 1,kx
2100  dduo3n(l,kk) = ( rbar(kk) * 1.00e-06 )
c
2200  continue
c
#ifdef seanal
      if ( seasnl )  then
      if ( iseasn .eq. 1 )  then
           do 2300 k=1,kx
           do 2300 l=1,19
           xo3ws(l,k) = dduo3n(l,k)
 2300 continue
      else
           do 2400 k=1,kx
           do 2400 l=1,19
           xo3sf(l,k) = dduo3n(l,k)
 2400 continue
      end if
      end if
#endif
c
      if ( cprto3 )  then
           write ( 6 , 301 )
           if ( .not. seasnl )  write ( 6 , 217 )
           if (    seasnl    )  write ( 6 , 212 )  seasns(iseasn)
           write ( 6 , 300 )    dduo3n
c          write ( 6 , 300 )  ((dduo3n(20-j,k),j=1,19),k=1,kx)
      end if
c
cpunchif ( pcho3 )  then
cpunch     if ( .not. seasnl )  write ( 7 , 217 )
cpunch     if (    seasnl    )  write ( 7 , 212 )  seasns(iseasn)
cpunch     write ( 7 , 300 )    dduo3n
c          write ( 7 , 300 )  ((dduo3n(20-j,k),j=1,19),k=1,kx)
cpunchend if
2500  continue
c
      return
      end

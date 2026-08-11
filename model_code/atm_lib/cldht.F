      subroutine cldht
#ifdef fixcld
c
c  This subroutine computes the heights of the cloud tops
c  and bottoms for the fixed cloud model.  The observed data
c  are from London (1954, 1957).  This data is a function of 10 deg.
c  latitude bands (0-10, 10-20 and etc.), season and height 
c  in the orginal paper and only for
c  the Northern Hemisphere for various cloud types.
c  Dick and Suki averaged the four seasons together to get annual 
c  mean cloud heights for three type of clouds (hi, middle and low).
c  Somebody also interpolated the data from the 10 deg latitude
c  bands to 5 deg bands.  At the equator, this interpolation
c  was more like an extrapolation.
c  These heights were then put in pressure coordinates using
c  a Skew-T diagram which assumes a "standard atmosphere".
c
c  The Dick and Suki (Ron pressures) data follow:
c
c TYPE:      hi      |   middle    |            low
c     |              |             |     ctop    |   cbase
c Lat |   ht   press |  ht   press |  ht   press |  ht   press
c 0   |  9.80  272   | 4.35  590   | 3.00  700   | 1.40  855
c 5   |  9.82  271   | 4.40  586   | 3.04  698   | 1.47  848
c 10  | 10.13  259   | 4.45  581   | 3.08  693   | 1.61  836
c 15  | 10.35  250   | 4.50  578   | 3.08  693   | 1.70  828 
c 20  | 10.50  244   | 4.50  578   | 3.01  699   | 1.72  825
c 25  | 10.50  244   | 4.41  584   | 2.91  710   | 1.71  826
c 30  | 10.38  248   | 4.26  595   | 2.80  719   | 1.70  828
c 35  | 10.03  263   | 4.10  614   | 2.70  729   | 1.65  830
c 40  |  9.44  285   | 3.92  621   | 2.60  735   | 1.58  839
c 45  |  8.65  322   | 3.79  633   | 2.47  750   | 1.50  846
c 50  |  7.97  357   | 3.67  647   | 2.35  760   | 1.40  859
c 55  |  7.55  379   | 3.56  651   | 2.24  770   | 1.31  867
c 60  |  7.29  392   | 3.51  657   | 2.17  780   | 1.25  871
c 65  |  7.13  401   | 3.50  658   | 2.10  783   | 1.20  875
c 70  |  7.03  406   | 3.48  659   | 2.03  788   | 1.12  881
c 75  |  7.01  409   | 3.44  660   | 1.98  795   | 1.05  890
c 80  |  6.99  410   | 3.43  661   | 1.91  800   | 1.02  891
c 85  |  6.98  411   | 3.43  661   | 1.88  803   | 1.00  896
c 90  |  6.98  411   | 3.43  661   | 1.87  804   | 1.00  896
c
c  Note that the heights are in kilometers and the pressures
c  in millibars.
c
c
c
c
#include "parm.h"
#include "parmncc.h"
#include "parmkp.h"
      parameter ( lx=ncc+1 )
      parameter (latobs=19)
c
#include "sgwgt1.h"
c
#include "latdat.h"
c
c  The data in the following arrays run pole to equator,
c  starting at 90, 80, 70....10, 0.
      dimension ciht(10), asht(10), cltop(10), clbase(10)
c
c  Observed high cloud heights running pole to equator (mb)
      data ciht / 411, 410, 406, 392, 357, 285, 248, 244, 259, 272 /
c
c  Observed middle cloud heights running pole to equator (mb)
      data asht / 661, 661, 659, 657, 647, 621, 595, 578, 581, 590 /
c
c  Observed low cloud top heights running pole to equator (mb)
ccc   data cltop / 804, 800, 788, 780, 760, 735, 719, 699, 693, 700 /
      data cltop / 804, 800, 788, 780, 766, 735, 719, 699, 693, 700 /
c  The above data statement was changed so that this code would
c  reproduce exactly the indexes used in the 9 level model.  The 6 mb
c  error should not affect the results very much....if at all.
c
c  Observed low cloud bottom heights running pole to equator (mb)
      data clbase / 896, 891, 881, 871, 859, 839, 828, 825, 836, 855 /
c
c
c  Initialize cloud heights to -1 for use in radiation code
c
      do 5 j=1,latobs
       do 5 nl=1,5
        kkth(j,nl) = -1
        kkbh(j,nl) = -1
    5 continue
c
c  First compute the cloud top indexes for the high clouds
c  nl is the index for cloud type nl=1 => high, nl=3 => low
c  
      nl = 1
      call cldint(ciht, kkth, nl)
c
c  Set cloud bottom height equal to cloud top height.
c  This assumes cirrus clouds are one level thick.
      do 10 j=1,latobs
       kkbh(j,nl) = kkth(j,nl)
   10 continue
c
c  Second compute the cloud top indexes for the middle clouds
c  nl is the index for cloud type nl=2 => middle
c 
      nl = 2
      call cldint(asht, kkth, nl)
c
c  Set cloud bottom height equal to cloud top height.
c  This assumes middle clouds are one level thick.
      do 20 j=1,latobs
       kkbh(j,nl) = kkth(j,nl)
   20 continue
c
c  Third compute the cloud top indexes for the low clouds
c  nl is the index for cloud type nl=3 => low
c 
      nl = 3
      call cldint(cltop, kkth, nl)
c
c  Lastly compute the cloud bottom indexes for the low clouds
c  This assumes that low clouds can be thicker than one level.
c  nl is the index for cloud type nl=3 => low
c 
      nl = 3
      call cldint(clbase, kkbh, nl)
c
#endif
      return
      end

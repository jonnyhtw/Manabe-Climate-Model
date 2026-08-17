c====================== include file "seasnl.h" ========================
c
#ifdef worldocn
# ifdef seasonal
c     variables used to cycle through seasonal surface boundary 
c     conditions:
c     calculated in "tmngr.F" and used in "surfbc.F"
c
#  ifdef oldtime
c     replicate the calendar of the old Cyber world ocean model runs.
c     assumes 12 equally spaced months when cycling through the
c     boundary conditions.  Month name and day of the month variables
c     for a 365 day julian calendar are kept for labeling purposes, 
c     but are not used when cycling through sets of 12 monthly averages
c     to simulate the seasonal cycle of surface boundary conditions
#  else
c     when cycling through sets of 12 monthly averages to simulate the
c     seasonal cycle of surface boundary conditions we assume that the
c     monthly means are representative of the mid-point of the month.
c     to determine the current boundary conditions we need to determine
c     the two monthly midpoints closest to the current time (the month
c     ahead "mahd" and month behind "mbhd").  additionally, we need to
c     determine the appropriate weights ("wmahd" & "wmbhd") for the 2
c     monthly values that bracket the current time.  so, at midnight
c     on the morning of September 1st, "mahd"=9 and "mbhd"=8, and the
c     weightings are equal ("wmahd" = "wmbhd" = 0.50).  at midnight
c     on the morning of September 16th, when 15.0 days are completed
c     and 15.0 days remain in the month (the midpoint of the month), 
c     the wieght for month 9 boundary conditions is 1.00.  between
c     midnight September 1st and midnight September 16th the weighting
c     for month 9 would have linearly increased from 0.50 to 1.00.
c     the September monthly means will have some influence on the
c     boundary conditions from just after the middle of August to
c     just before the middle of October.
#  endif
c     
c     
c     mahd & mbhd   = month ahead & month behind
c     wmahd & wmbhd = weighting factors for month ahead & month behind
c
      common /seasn/ mahd, mbhd, wmahd, wmbhd
# endif
#endif
c

c====================== include file "cinlnd.h" ==========================
c
c     establish communication with inland seas via mixing of tracers
c     in water columns, separated by land.  the rate of mixing of the
c     water columns, the vertical extent of the mixing, and the location
c     of the pairs of points to be mixed, are specifiec in "blkdta.F".
c     each of the pairs of water columns to be mixed must lie in
c     the same j-row.  location values are extremely model dependent.
c
c     "iinlnd" specifies the 2 i column locations of the water columns
c              to be mixed for each of the "ninlnd" pairs of points  
c     "jinlnd" specifies the j row of each of the "ninlnd" pairs of pts
c     "kinlnd" specifies the vertical extent of the mixing by setting 
c              the top & bottom model k-levels to be mixed   
c     "vinlnd" specifies the rate at which of the water columns are to
c              be mixed [cubic cm per second] equivalent to 
c              the flow to the east = the flow to the west
c     "binlnd" this calculated variable is the fraction of the model
c              grid boxes specified that must be mixed per second
c              in order to meet the requested "vinlnd" mixing rate
c
c        Note: for "iinlnd", and "kinlnd" pairs the (n,1) element's 
c              value should be < the corresponding (n,2) value
c              similarly, the "jinlnd"(n) values should be in 
c              ascending order if ninlnd > 1
c
#if defined worldocn
# if defined inland
# if defined ocnr15 
c     for the world ocean model that is coupled to the 
c     atmosphere model, inland mixing is done between the Mediterranean
c     and the Atlantic
c
      parameter (ninlnd = 1)
# endif
# if defined ocnr30 
c     for the world ocean model that is coupled to the 
c     atmosphere model, inland mixing is done between the Mediterranean
c     and the Atlantic and between Hudson Bay and the Atlantic
c
      parameter (ninlnd = 2)
# endif
# else
      parameter (ninlnd = ?)
# endif
#endif
c
      common /inlnd/ iinlnd(ninlnd,2), jinlnd(ninlnd), kinlnd(ninlnd,2)
     $              ,vinlnd(ninlnd), binlnd(ninlnd,2)
c

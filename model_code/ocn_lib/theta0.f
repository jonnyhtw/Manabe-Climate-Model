      function theta0 (ydeg, depth)
c
c=======================================================================
c     this subroutine returns estimates of global mean potential
c     temperature for model initialization as a function of depth.
c     it is used to produce a reference thermal stratification for the
c     upper 2000m of the MOM's test case.  below 2000m, the
c     potential temperature returned is 2.0 degrees C.  surface
c     values are set slightly above 18.4 degrees C at the reference
c     latitude "reflat".
c     the estimates are produced from a 7th order ploynomial fit to
c     the annual mean world ocean potential temperature observations
c     of Levitus (1982).
c
c     input [units]:
c       a latitdue (ydeg): [degrees]
c       a zt value (depth): [centimeters]
c     output [units]:
c       potential temperature estimate (est): [degrees centigrade]
c
c     variables:
c       coeft     = coefficients for the polynomial fit of potential
c                   temperature vs. depth
c       reflat    = reference latitude at which observed surface
c                   temperatures approximately equal coeft(1)
c       factor    = the ratio of the cosine of the latitude requested
c                   ("ydeg") to the reference latitude ("reflat")
c                   used to scale the upper 2000 meters of the vertical
c                   temperature profile
c       tmin,tmax = the minumum and maximum potential temperatures
c                   allowed at the time of model initialization
c
c     reference:
c       Levitus, S., Climatological atlas of the world ocean, NOAA
c     Prof. Paper 13, US Gov't printing Office, Washington, DC, 1982.
c
      parameter (ndeg=7)
      dimension coeft(ndeg+1)
      save coeft, tmin, tmax, reflat
      data coeft / 0.184231944E+02,-0.430306621E-01, 0.607121504E-04
     $           ,-0.523806281E-07, 0.272989082E-10,-0.833224666E-14
     $           , 0.136974583E-17,-0.935923382E-22/
      data tmin, tmax, reflat /2.0, 25.0, 34.0/
c
c=======================================================================
c
      c0 = 0.0
      pi = atan(1.0) * 4.0
      refcos = abs(cos(pi*reflat/180.))
c
      coslat = abs(cos(pi*ydeg/180.))
      factor = coslat/refcos
      z = depth * 0.01
c
      if (z .gt. 2000.) then
        est = 2.0
      else
        est = c0
        do 100 nn=1,ndeg+1
          bb = 1.0
          if (nn.gt.1) bb = z**(nn-1)
          est = est + coeft(nn)*bb
100     continue
        est = est * factor
      endif
c
      if (est .gt. tmax) est = tmax
      if (est .lt. tmin) est = tmin
c
      theta0 = est
      return
      end

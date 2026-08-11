      subroutine depths
c      program depths
c
c     due to the simple UNIX linker, only one main program may
c     be in the directory at one time. To run this program, wipe out
c     the subroutine call & uncomment the program line. all other
c     main programs must be either removed from the directory or
c     commented out also!
c     note: this is no problem if "makefiles" are used for compiling
c           & linking.
c
c=======================================================================
c      interactive program to allow generation of model level thickness
c      file "thick.h"
c=======================================================================
c
      parameter (kmax=1000)
cSGI
cc    implicit double precision (a-h,o-z)
      implicit real*16 (a-h,o-z)
cSGI
      character *1 ch
      dimension z(kmax), dzt(kmax), zw(kmax), dzw(0:kmax), thick(kmax)
cSGI
cc    double precision fun
      real*16 fun
cc    fun (k,sigma,c) = dexp(-(k/sigma)**2) + c
      fun (k,sigma,c) =  exp(-(k/sigma)**2) + c
cSGI
c
      print *,' '
      print *, 'fit model level thickness "dzt(k)" to a gaussian with'
      print *, 'the following form:'
      print *,' '
      print *, 'dzt(km+1-k) = exp(-(k/sigma)**2) + offset'
      print *,' '
      print *, 'for levels k=1..km by shooting for "sigma" using a'
      print *, 'newton or secant method. the distribution can be '
      print *, 'changed by altering "offset" and the initial "sigma".'
      print *, 'for a start, try km=15 with a 30m thick 1st level'
      print *, '(only one surface level with this thickness) and'
      print *, 'a 5700m deep ocean. now try the following:'
      print *, 'offset=0.0  & sigma=5'
      print *, 'offset=0.01 & sigma=8'
      print *, 'offset=0.01 & sigma=3'
      print *,' '
      print *,' '
      print *, 'enter the number of model levels "km" ='
      read *, kmt
      print *,'enter thickness of the first level (meters) ='
      read *, dzt1
      print *,'enter number of levels with this thickness='
      read *, ktop
      km  = kmt - ktop + 1
      top = (ktop-1)*dzt1
      print *,'enter depth at the bottom of the last level (meters) ='
      read *, depth
      depth = depth - top
      print *,'enter "offset" ='
      read *, offset
      print *,'enter initial "sigma" ='
      read *, sigma
      if (depth/km .lt. dzt1) then
        print *,'no solution!  try another configuration.'
        stop
      endif
c
c     shoot for sigma using the newton method
c
      maxl = 200
      c0   = 0.0
      l    = 0
      deps = 1.e-3
c
      do 300 l=1,maxl
c
c       first guess
c
        sum = c0
        do 100 k=1,km
          dzt(km+1-k) = fun(k,sigma,offset)
          sum = sum + dzt(km+1-k)
100     continue
c
        a  = dzt1/dzt(1)
        en = depth - a*sum
c
c       2nd guess
c
        dsigma =  sigma*1.e-10
        sigmap = sigma + dsigma
        sum = c0
        do 200 k=1,km
          dzt(km+1-k) = fun(k,sigmap,offset)
          sum = sum + dzt(km+1-k)
200     continue
c
        a   = dzt1/dzt(1)
        enp = depth - a*sum
        print *,'trial shoot #=',l,' sigma=',sigmap,' error=',enp
c
c       construct new 1st guess (newton or secant method)
c
        sigma = sigmap - enp * dsigma / (enp - en)
c
        if (abs(enp) .le. deps) go to 301
300   continue
      print *,' Either there is no convergence or it is slow.'
      print *,' Try another set of inputs.'
      stop '=>done'
301   continue
c
      do 400 k=1,km
        dzt(k) = a*dzt(k)
400   continue
c
      zw(1)  = top + dzt(1)
      z(1)   = top + 0.5*dzt(1)
      dzw(0) = 0.5*dzt(1)
      do 500 k=2,km
        z(k) = z(k-1) + 0.5*(dzt(k)+dzt(k-1))
        zw(k) = zw(k-1) + dzt(k)
        dzw(k-1) = z(k) - z(k-1)
500   continue
c
      print *,' '
      print *,'gaussian fit for dzt(km+1-k) k=',ktop,'..',kmt
      print *,'dzt(km+1-k) = exp(-(k/',sigmap,')**2)+',offset
      print *,' '
      write (6,'(1x,a3,a10,3a15)') 'k','z(k)','zw(k)','dzt(k)','dzw(k)'
      do 600 k=1,kmt
        if (k .lt. ktop) then
          dzw1 = dzt1
          thick(k) = dzt1
          if (k .eq. 1) dzw1 = 0.5*dzt1
          write (6,'(1x,i3,4(1x,g14.7))') k,
     1      (k-0.5)*dzt1, k*dzt1, dzt1, dzw1
        else
          m = k - ktop + 1
          thick(k) = dzt(m)
          if (m .eq. 1) dzw(m-1) = 2.0*dzw(m-1)
          write (6,'(1x,i3,4(1x,g14.7))') k, z(m), zw(m), dzt(m),
     $    dzw(m-1)
        endif
600   continue
c
      print *,' do you want to produce a "thick.h" file ? (y/n)'
      print *,' (note: this will write over any old "thick.h" file)'
      read (*,'(a1)') ch
      if (ch .eq. 'Y' .or. ch .eq. 'y') then
        open (unit=40, file='thick.h')
c
c       write sample thick.h file "dzt" data statements
c
        ks = 1
        ke = ks + 4
        ke = min(kmt,ke)
        write(40,1958) (thick(k),k=1,ke)
        do 700 ii=1,99
          ks = ks + 5
          if (ks .gt. kmt) goto 701
          ke = ks + 4
          ke = min(kmt,ke)
          write(40,1959) (thick(k),k=ks,ke)
700     continue
701     continue
        write(40,1960) km, sigmap, offset, dzt1, ktop, zw(km)
        write(40,1961)
c
1958  format('c====================== include file "thick.h" ',
     $   '========================',/
     $   'c     specifications for box thickness (cm) in the vertical'/
     $   ,'c     produced by program "depths.f"'/
     $   ,'c'/'      data dzt /'f8.2,'e2',4(',',f8.2,'e2'))
1959  format('     $      ',5(',',f8.2,'e2'))
1960  format('     $   /',/,'c'/
     $   'c     dzt = model level thicknesses computed as',/
     $   'c     dzt(km+1-k)=exp(-(k/sigma)**2) + offset) for k=1,km',/,
     $   'c      km =',i4,/
     $   'c      sigma=',g14.7,/
     $   'c      offset=',g14.7,/
     $   'c      1st level thickness =',g14.7,/
     $   'c      # of levels with this thickness =',i3,/
     $   'c      depth at bottom of last level =',g14.7)
1961  format('c'/)
c
        close(40)
      endif
      stop
      end

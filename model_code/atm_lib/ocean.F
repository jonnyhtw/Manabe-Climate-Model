#if !defined cpldleg3 && !defined cpldleg1
       subroutine ocean
c
c           called from: airsea
c          entry points: none
c   subroutines called : ascmm
c                        fixsst (fixed sst models only)
c
c   ocean updates 'to' (mixed layer ocn temp) and 'yce' (sea ice).  for
c    a fixed sst model, both varialbles are updated in fixsst.  for a
c    mixed layer model, this is done by computing the mixed layer heat
c    content from the old values of 'to' and 'yce', then adding the
c    contributions from the ocean heat flux 'qao' and snowfall 'sno'.
c    the new values of 'to' and 'yce' are then computed from the
c    updated mixed layer heat content.
c
#include "parm.h"
c
      parameter ( intglo= 87  )
c
c
c
#include "anlocn.h"
c
#include "bcs.h"
c
#include "ctlblk.h"
c
#include "gauss.h"
c
#include "grdfld.h"
c
#include "qa.h"
c
#include "param.h"
c
#include "surfo.h"
c
      dimension yceold(ix,il),sea(ix,il)
c
#ifdef ice_dynamics
c
c variables in ice model used in ocean model
c
#include "varice.h"
c
#include "variceco.h"
#endif
c
#include "zmanl.h"
c
c
c  the below assumes a mixed layer depth of 50 meters (5000 cm)
      data    co  / 5000.0 /,   alhfrz /    80.0 /
      data yrday  /  365.0 /, yrmon  /   12.0 /, daysec / 86400.0 /
      data alandtm /  290.0 /, daymin / 1440.0 /
c
#ifdef ice_dynamics
c
c constant from the dynamic ice model (cfi)
c
      data rhoic /1./, alhfice /80./
#endif
c
      data cyear / 365.0 /
c
c     <<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      xi   = 1.0/float(ix)
      do 50 j=1,il
      do 50 i=1,ix
       if (ts(i,j) .gt. 0.0) then
         sea(i,j)=0.0
       else
         sea(i,j)=1.0
       endif
 50   continue
c
      write(6,480) ktauo,ktau
c
c   units:  idto=sec/tmstep, co=cal/(cm*cm*deg), min=60 sec
c           qao=cal/(cm*cm*min)      so total units are deg/tmstep
c   over land to = alandtm  (290 deg k).
c
c
c save old ice thickness value (after adding snowfall) for use in
c  computing melt
      do 285 j=1,il
       do 285 i=1,ix
#ifdef ice_dynamics
c
c used stored liquid equivalent ice in ice model (cfi)
c
        yceold(i,j)=yceol(i,j)+sno(i,j)*float(idto)
#else
        yceold(i,j)=yce(i,j)+sno(i,j)*float(idto)
#endif
  285 continue
c
      if (.not.swamp) then
#ifdef mixed_layer
#ifdef old_qflux_fixedsst
          astday = amod(day, cyear)
c
c for q-flux mixed layer models get heat flux
c  at the bottom of the mixed layer; otherwise,
c  set it to zero.
#ifdef q_flux
            call getqf (astday)
#else
            do 60 j=1,il
              do 60 i=1,ix
                qbot(i,j) = 0.
   60       continue
#endif
#endif
#ifdef ice_dynamics
c
c calculate ocean mixed-layer temperature from heat flux from ice (cfi)
c
        do 19 j = 1, il
        do 19 i = 1, ix
c          htcont = (to(i,j)-tfro)*co + qocean(i,j)*daymin
c          to(i,j) = htcont/co + tfro
           to(i,j) = to(i,j) + qocean(i,j)*daymin/co
   19   continue
#else
c
c compute to and yce using qao and sno
        do 300 j=1,il
          do 300 i=1,ix
c
c compute heat content of mixed layer (relative to tfro)
            htcont = (to(i,j)-tfro)*co - yce(i,j)*alhfrz
c
c update heat content using ocean heat fluxes and snowfall
            htcont = htcont + qao(i,j)*daymin + qbot(i,j)*daymin
     $               - sno(i,j)*float(idto)*alhfrz
c
c determine new values of to and yce from new heat content
            if (htcont .le. 0.0) then
              to (i,j) = tfro
              yce(i,j) = -htcont/alhfrz
            else
              to (i,j) = tfro + htcont/co
              yce(i,j) = 0.0
            endif
  300   continue
#endif
#endif
#ifdef fixed_sst
c
#ifdef old_qflux_fixedsst
c
c specify to and yce using observed data
          astday = amod (day,yrday)
          call fixsst (astday)
c
#endif
c
c specify to using observed data
	do 330 j=1,il
	  do 330 i=1,ix
	    to(i,j) = bcto(i,j)
  330   continue
#ifndef ice_dynamics
c
c for version without ice dynamics, specify yce using observed data
c (otherwise, sea ice thickness and concentration are specified in
c subroutine growths)
	do 335 j=1,il
	  do 335 i=1,ix
            yce(i,j) = bcyce(i,j)
  335   continue
#endif
#endif
c
#if defined restore && defined mixed_layer
c
c restoring to observed sst and sea ice with a time constant (cfi)
c
	call heatres
#endif
      endif
c
c  set to and yce for swamp models
      if (swamp) then
        do 333 j= 1,il
         do 333 i =1,ix
          yce(i,j) = 0.0
          to (i,j) = abs(ts(i,j)) * sea(i,j)
  333   continue
      endif
c
c compute ice melt and mask land values for 'to' and 'yce'
       do 315 j=1,il
        do 315 i=1,ix
         ycemlt(i,j) = ycemlt(i,j) + sea(i,j)*(yce(i,j)-yceold(i,j))
         to    (i,j) = sea(i,j)*to (i,j) + (1.0-sea(i,j))*alandtm
         yce   (i,j) = sea(i,j)*yce(i,j)
 315   continue
c
      do 200 j=1,il
       do 200 i=1,ix
        zmto (j) = zmto (j) +  to(i,j)*sea(i,j)
        zmqo (j) = zmqo (j) + qao(i,j)*sea(i,j)
        zmyce(j) = zmyce(j) + yce(i,j)*sea(i,j)
        zmsno(j) = zmsno(j) + sno(i,j)*sea(i,j)
200   continue
c
c    if not the end  of  the month return
c
c
c    ktauo  = the timestep number in the ocean
c    ifcsto = the forecast time in seconds in the ocean
c    idto   = the timestep in seconds in the ocean
c    dayo   = the ocean day
c    oyear  = the ocean year
c
      ktauo  = ktauo+1
      kpero  = kpero+1
      ifcsto = ifcsto+idto
      dayo   = float(ifcsto)/daysec
      oyear  = dayo/yrday
c
c   test for the end of the month
      if( float(ifcsto+idto)/daysec .gt. apo ) then
c
        write (intglo,102) expnr,day,iday,month,iyear
c
c       compute monthly means in ascmm
        call ascmm
        mono   = mono+1
        apo1   = yrday/yrmon
        apo    = apo + apo1
        kpero  = 0
c
c  the test below indicates the end of the year
        if (mono .eq. (yrmon+1)) then
          mono = 1
          its  = its+1
        endif
      endif
c
  102 format (2a8,f8.0,i4,a4,i4)
  480 format (1x,'in ocean at ktauo=',i8,'atmos ktau=',i8)
c
c
c  write ocean history tape
c
      if (rsrt) call rstrto
c
      return
      end
#endif

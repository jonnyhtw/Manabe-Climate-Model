c
c constants in ice model
c
      common / varcona / deltat, dticet, dlam, dlamd, r, piice
      common / varconb / cliqsol, rhoice, fcor, pstari
c
c ice velocity at both hemispheres
c
      common / varvelall / uic(ix,il,3), vic(ix,il,3)
c
c wind and ocean drag turning angle
c
      common / varangle / phi, sinphi, cosphi, theta, sinth, costh
c
c varying latitude points
c
      common / varlat / arlat(iy), arlatc(iy+1)
c
c varying Coriolis parameters
c
      common / varcor / fcorp(iy), fcorc(iy+1)
c
c center point mask for ice model for each hemisphere
c
      common / varmaskh / out(ix+1,iy+1), heffm(ix+1,iy+1)
c
c velocity mask for ice model for each hemisphere
c
      common / varmaskv / uvmx(ix,iy), uvmy(ix,iy), uvmb(ix,iy)
c
c multiplicative factors and forces in dynamic calculation
c
      common / varmul / xx(ix,iy), yy(ix,iy),
     >           axi(ix,iy), ayi(ix,iy), bxi(ix,iy), byi(ix,iy),
     >           dwatn(ix,iy), dwatx(ix,iy), dwaty(ix,iy),
     >           dairn(ix,iy), dairx(ix,iy), dairy(ix,iy)
c
c geostrophic wind and ocean velocities in dynamic calculation
c
      common / vargeo / gairx(ix,iy), gairy(ix,iy),
     >           gwatx(ix,iy), gwaty(ix,iy)
c
c ice velocity in dynamic calculation
c
      common / varvel / uice(ix,iy,3), vice(ix,iy,3)
c
c ice thickness and compactness in dynamic calculation
c
      common / varthi / hef(ix+1,iy+1,3), are(ix+1,iy+1,3)
c
c wind stress vectors in dynamic calculation
c
      common / varwind / tauxo(ix+1,iy+1), tauyo(ix+1,iy+1),
     >           tauxin(ix,iy), tauyin(ix,iy)
c
c internal ice force vectors in dynamic calculation
c
      common / varficxy / fxi(ix+1,iy+1), fyi(ix+1,iy+1)
c
c ice strength variables in dynamic calculation
c
      common / varicest / xmpress(ix+1,iy+1), presur(ix+1,iy+1)
c
c indicator in dynamic calculations
c
      common / varindyn / np, nl, ihem, jlimit

#ifdef restore
c
c restoring heat flux (flux adjustment term)
c
      common / varestore / qres(ix,il)
#endif
c
#ifdef icecpld
c
c depth of first few layers in ocean model
c
      common / vardepth / dzti(5)
#endif
c
c wind stress and ocean currents (at corner of grid)
c
      common / varwinstr / taxic(ix,il), tayic(ix,il)
      common / varocecur / gwxic(ix,il), gwyic(ix,il)
c
c ice thickness, ice compactness and ice internal pressure
c
      common / varthick / heff(ix,il,3), area(ix,il,3), ar(ix,il)
      common / varpres / presu(ix,il)
c
c heat flux and moisture flux over open water and ice
c
      common / varqfluxes / qfow(ix,il), qfice(ix,il) 
      common / varqfluxatm / qfowatm(ix), qficeatm(ix)
      common / varvflux / vapori(ix), vaporo(ix)
c
c heat flux and water flux from ice to ocean
c
      common / varoflux / qocean(ix,il), pocean(ix,il)
      common / varoiflux / qoceani(ix,il), poceani(ix,il)
c
c internal ice force in the horizontal directions
c
      common / varfice / fxice(ix,il), fyice(ix,il)
      common / varifice / fxicei(ix,il), fyicei(ix,il)
c
c wind stress vector (surface wind) over open water and ice
c
      common / varwstress / tauxic(ix,il), tauyic(ix,il)
c     
c open water latent heat, mixing ratio, and drag parameter when ice is present
c
      common / varhow / hltow(ix), rstsow(ix), rhokow(ix)
c
c drag parameters over open water and grid point average
c
      common / varcd / cdow(ix,lmx), cdio(ix,lmx)
c 
c land, ice compactness
c
      common / varlaicc / arlaic(ix)
c
c sensible and latent heat flux for open water
c
      common / varqfow / shfow(ix,il), altntow(ix,il)
c
c total growth rate for ice thickness and compactness
c
      common / varrate / fice(ix,il), fow(ix,il)
      common / vargrow / sarea(ix,il), sice(ix,il)
c
c center point mask
c
      common / varmask / heffmg(ix,il)
c
c ocean mixed layer heat and ice thickness normalized by compactness
c
      common / varheat / heatml(ix,il), hice(ix,il)
c
c intermediate ocean temperature and heat fluxes 
c
      common / varocei / toint(ix,il), qfowi(ix,il), qficei(ix,il)
c
c gaussian weighting function
c
      common / varweight / weightin(il)
c
c ice velocities at grid center for history
c
      common / varicevel / uicec(ix,il), vicec(ix,il) 
      common / variceadd / uicea(ix,il), vicea(ix,il)
c
c heat fluxes and growth rates at grid center for history
c
      common / varqfluxc / qficec(ix,il), qfowc(ix,il)
      common / vargrrate / siceic(ix,il), siceow(ix,il), hmlc(ix,il)
      common / vargrowc / sicec(ix,il), sareac(ix,il), qflc(ix,il)
c
c ice growth from advection, thermodynamics, and from ice model
c
      common / variceadth / giceadv(ix,il), gicetherm(ix,il) 
      common / varareadth / gareadv(ix,il), garetherm(ix,il) 
      common / varicet / gicet(ix,il), garet(ix,il) 
c
c surface temperatures and ice temperature 
c
      common / varsurft / tts(ix,il), tsice(ix,il), tscd(ix,il)
c
c testing variable
c
      common / ttice / pkiceo(ix), cdfi(ix,lmx), cdfo(ix,lmx)
      common / ttic / yceol(ix,il), ths(ix,il)
      common / hcice / hswaic(ix), hlwaic(ix), hsheic(ix), hlheic(ix)
      common / hcowa / hswaow(ix), hlwaow(ix), hsheow(ix), hlheow(ix)
c
c daily average wind velocities and stress
c
      common / wiav / windvx(ix,il),windvy(ix,il),tax(ix,il),tay(ix,il)
c
c depth of oceanic mixed layer
c
      common / varmixl / depoml(ix,il)
c
c indicator variables for restart
c
c     common / inice / irestra, irestic

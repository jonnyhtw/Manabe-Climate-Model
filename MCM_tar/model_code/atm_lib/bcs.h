c
#if defined fixed_sst || defined restore
      common / bcs / bcto(ix,il), bcyce(ix,il)
#ifdef ice_dynamics
      common / bcs / bcconc(ix,il)
#endif
#endif
c
#ifdef q_flux
      common / bcs / bcqflx(ix,il)
#endif
c
#if defined fixed_sst || defined restore
c bcto   = current day's observed mixed layer temperature
c bcyce  = current day's observed sea ice thickness
#ifdef ice_dynamics
c bcconc = current day's observed sea ice concentration
#endif
#endif
#ifdef q_flux
c bcqflx = current day's heat flux adjustment
#endif

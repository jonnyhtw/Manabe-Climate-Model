      common /  qa    / qao(ix,il), sno(ix,il), qbot(ix,il)
#if defined cpldleg3 || defined cpldleg1
      common /  qa    / taux (ix,il), tauy(ix,il)
      common /  qa    / per  (ix,il)
      common /  qa    / ijrnf(ix,il)
# ifdef windbc
      common /  qa    / winds(ix,il), wind3(ix,il)
# endif
#endif
#ifdef cpldleg1
      common /  qa    / qaot (ix,il), snot (ix,il)
      common /  qa    / tauxt(ix,il), tauyt(ix,il)
      common /  qa    / pert (ix,il)
# ifdef windbc
      common /  qa    / windst(ix,il), wind3t(ix,il)
# endif
#endif

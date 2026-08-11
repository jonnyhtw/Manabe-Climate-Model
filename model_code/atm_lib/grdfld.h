      common / grdfld / ts     (ix,il), shfx   (ix,il)
#ifdef rad
      common / grdfld / slwdp  (ix,il)
#endif
#ifdef moist
      common / grdfld / altnt  (ix,il), precip (ix,il)
      common / grdfld / evapor (ix,il), snowfl (ix,il)
#endif
#ifdef hydrol
      common / grdfld / soilm  (ix,il), runoff (ix,il)
      common / grdfld / snwdpt (ix,il), snwmlt (ix,il)
      common / grdfld / smelts (ix,il), calb   (ix,il)
#endif
#ifdef moist
      common / grdfld / sublim (ix,il)
#endif
      common / grdfld / ycemlt (ix,il)
#if defined cpldleg3 || defined cpldleg1
      common / grdfld / snwrnf (ix,il)
#endif

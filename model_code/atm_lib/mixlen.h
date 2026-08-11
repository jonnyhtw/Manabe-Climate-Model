      common / mixlen / elmix(kx), zstd(kx), zhstd(kp)
#ifdef garratt
      common / mixlen / fgarr(kx)
#endif
c
c elmix(kx)  =  mixing lengths [cm] for vertical diffusion
c zstd(kx)   =  heights [cm] of full levels assuming standard lapse rate
c zhstd(kp)  =  heights [cm] of half levels assuming standard lapse rate
c fgarr(kx)  =  
c
c Note: elmix, zstd, zhstd are used in PARMTR to compute ks and
c the vertical diffusion coefficients.

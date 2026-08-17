#ifdef fixcld
      common / latdat / ccd(19,nc),dduo3n(19,kx),kkbh(19,nc)
      common / latdat / kkth(19,nc)
#endif
#ifdef prdcld
      common / latdat / dduo3n(19,kx)
#endif
#ifdef seanal
      common / latdat / xo3sf(19,kx), xo3ws(19,kx)
      common / latdat / ccosz(il),ttauda(il)
#endif
#ifdef annlmn
      common / latdat / ccosz(19),ttauda(19)
#endif
#ifdef dry
      common / latdat / rclim(19,0:kp)
#endif

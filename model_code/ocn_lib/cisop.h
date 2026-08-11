c====================== include file "cisop.h" =========================
c
c           variables for isopycnal diffusion scheme
c
c    "ahisop" = along isopycnal diffusion coefficient for mixing of
c               tracers  [cm**2/sec]
c               ("ah" = in "chmix.h is used as background horizontal 
c               diffusion to prevent checkerboarding)
c    "slmxr"  = maximum slope of isopycnal allowed in mixing scheme
c               limited to prevent excessively large vertical mixing
c               that could create numerical instabilities
c               a value of 200 translates to a slope of 1:200
c    "rx", "ry" & "rz" terms = store density differences used to compute 
c               isopycnal diffusion tensor components "fk#"
c    "rxp" = east - central density difference at the row j+1
c    "rx"  = east - central density difference at the current j row
c    "ryn" = north - central density difference (row j+1  -  row j)
c    "rys" = central - south density difference (row j  -  row j-1)
c    "rzp" = above - central vertical density difference at the row j+1
c    "rz"  = above - central vertical density difference at current row
c
c    "e" = different things in different places...
c      in isop0 = density gradients (1=zonal; 2=meridional; 3=vertical)
c      in isop1 = who knows what
c    "fk#" = isopycnal diffusion tensor components (compute five)
c          "fk1(3)" is resolved at center of eastern face of t,s box
c          "fk2(3)" is resolved at center of northern face of t,s box
c          "fk3(1-3)" are resolved at center of top face of t,s box
c     components (1,1) & (2,2) of tensor are assumed = 1.0
c     components (1,2) & (2,1) of tensor are assumed = 0.0
c
c
      common /cisop/ rxp(imt,km), ryn(imt,km), rzp(imt,kmp1),
     $               rx (imt,km), rys(imt,km), rz (imt,kmp1),
     $               fk1(imt,km,3:3), fk2(imt,km,3:3), fk3(imt,km,3),
     $               e(imt,kmp1,3), esav(imt,km,nt),
     $               ahisop, slmxr
c

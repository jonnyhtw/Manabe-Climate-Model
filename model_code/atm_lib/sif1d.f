      subroutine sif1d(ndx,len,src,tgt1,tgt2)
      dimension ndx(1),src(1),tgt1(1),tgt2(1)
      leng=len+1
      do 10 j=1,leng
        jj=ndx(j)
        tgt1(j)=src(jj)
        tgt2(j)=src(jj+1)
  10  continue
      return
      end

      subroutine sif1(ndx,len,src,tgt)
      dimension ndx(1),src(1),tgt(1)
      leng=len+1
      do 10 j=1,leng
        jj=ndx(j)
        tgt(j)=src(jj)
  10  continue
      return
      end

      subroutine qqb_QQb_ew2(p,msq)
      implicit none
      include 'types.f'
      include 'nf.f'
      include 'constants.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'zcouple.f'
      include 'qcdcouple.f'
      include 'breit.f'
      include 'sprods_com.f'
      include 'masses.f'
      include 'first.f'
      integer j,k
      logical first
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),ss,tt,uu,dM,DL,SL,
     .     aem0,alpha,aemmz,m,sw2,cw2,c1,c2,is,it,iu,r13,r23,
     .     qqb(-nf:nf,-1:1,-1:1),ttb(-1:1,-1:1)
      common/em/aemmz

      if(first) then
         first=.false.
         write(6,*) 'Heavy Quark mass:', mass2
      endif

      m=mass2
      sw2=xw
      cw2=1._dp-sw2
      dM=log(zmass**2/wmass**2)

c--- alpha_em(0)
c      aem0=1._dp/126.3_dp
c--- alpha -> alpha_em(0)
c      alpha=aem0
c--- alpha -> standard MCFM value
      alpha=aemmz


      call dotem(4,p,s)

      ss=s(1,2)
      r13=dabs(s(1,3))
      r23=dabs(s(2,3))
      is=1._dp/ss


      DL=alpha/4._dp/pi*log(ss/wmass**2)**2
      SL=alpha/4._dp/pi*log(ss/wmass**2)

      msq=0._dp
      qqb=0._dp
      ttb=0._dp

c--- qqb(-nf:nf,-1:1,-1:1) is the EW Sudakov logs w.r.t helicity
c--- process: qqb -> ttb
c--- arguments: 1st = (anti)quark; 2nd = helicity of q; 3rd = helicity of t

      qqb(-5,-1,-1) = 
     &    ( + SL * ( 3._dp*sw2**(-1) + 3._dp/2._dp*sw2**(-1)*cw2 + 
     &    sw2**(-1)*cw2*dM + 1._dp/6._dp*sw2*cw2**(-1) + 1._dp/9._dp*sw2*
     &    cw2**(-1)*dM - 1._dp/2._dp*m**2*wmass**(-2)*sw2**(-1) - log(
     &    r13*r23**(-1))*sw2**(-1)*cw2**(-1) + 2._dp*log(r13*r23**(-1))
     &    *cw2**(-1) - 8._dp/9._dp*log(r13*r23**(-1))*sw2*cw2**(-1) )
     &    + DL * (  - sw2**(-1) - 1._dp/2._dp*sw2**(-1)*cw2 - 1._dp/18._dp
     &    *sw2*cw2**(-1) ))

      qqb(-5,-1,1) = 
     &    ( + SL * ( 1._dp/2._dp + 3._dp/2._dp*sw2**(-1) + 3._dp/4._dp*
     &    sw2**(-1)*cw2 + 1._dp/2._dp*sw2**(-1)*cw2*dM + 1._dp/3._dp*dM + 
     &    17._dp/12._dp*sw2*cw2**(-1) + 17._dp/18._dp*sw2*cw2**(-1)*dM - 
     &    3._dp/4._dp*m**2*wmass**(-2)*sw2**(-1) + 4._dp/3._dp*log(r13*
     &    r23**(-1))*cw2**(-1) - 8._dp/9._dp*log(r13*r23**(-1))*sw2*
     &    cw2**(-1) )
     &    + DL * (  - 1._dp/6._dp - 1._dp/2._dp*sw2**(-1) - 1._dp/4._dp*
     &    sw2**(-1)*cw2 - 17._dp/36._dp*sw2*cw2**(-1) ))

      qqb(-5,1,-1) = 
     &    ( + SL * (  - 1._dp/2._dp + 3._dp/2._dp*sw2**(-1) + 3._dp/4._dp
     &    *sw2**(-1)*cw2 + 1._dp/2._dp*sw2**(-1)*cw2*dM - 1._dp/3._dp*dM + 
     &    5._dp/12._dp*sw2*cw2**(-1) + 5._dp/18._dp*sw2*cw2**(-1)*dM - 1._dp/
     &    4._dp*m**2*wmass**(-2)*sw2**(-1) + 2._dp/3._dp*log(r13*
     &    r23**(-1))*cw2**(-1) - 8._dp/9._dp*log(r13*r23**(-1))*sw2*
     &    cw2**(-1) )
     &    + DL * ( 1._dp/6._dp - 1._dp/2._dp*sw2**(-1) - 1._dp/4._dp*
     &    sw2**(-1)*cw2 - 5._dp/36._dp*sw2*cw2**(-1) ))

      qqb(-5,1,1) = 
     &    ( + SL * ( 5._dp/3._dp*sw2*cw2**(-1) + 10._dp/9._dp*sw2*
     &    cw2**(-1)*dM - 1._dp/2._dp*m**2*wmass**(-2)*sw2**(-1) - 
     &    8._dp/9._dp*log(r13*r23**(-1))*sw2*cw2**(-1) )
     &    + DL * (  - 5._dp/9._dp*sw2*cw2**(-1) ))

      qqb(-2,-1,-1) = 
     &    ( + SL * (  - 1._dp + 3._dp*sw2**(-1) + 3._dp/2._dp*
     &    sw2**(-1)*cw2 + sw2**(-1)*cw2*dM - 2._dp/3._dp*dM + 1._dp/6._dp*
     &    sw2*cw2**(-1) + 1._dp/9._dp*sw2*cw2**(-1)*dM - 1._dp/4._dp*m**2*
     &    wmass**(-2)*sw2**(-1) + log(r13*r23**(-1))*sw2**(-1)*
     &    cw2**(-1) - 8._dp/3._dp*log(r13*r23**(-1))*cw2**(-1) + 
     &    16._dp/9._dp*log(r13*r23**(-1))*sw2*cw2**(-1) )
     &    + DL * ( 1._dp/3._dp - sw2**(-1) - 1._dp/2._dp*sw2**(-1)*cw2 - 
     &    1._dp/18._dp*sw2*cw2**(-1) ))

      qqb(-2,-1,1) = 
     &    ( + SL * (  - 1._dp/2._dp + 3._dp/2._dp*sw2**(-1) + 3._dp/4._dp
     &    *sw2**(-1)*cw2 + 1._dp/2._dp*sw2**(-1)*cw2*dM - 1._dp/3._dp*dM + 
     &    17._dp/12._dp*sw2*cw2**(-1) + 17._dp/18._dp*sw2*cw2**(-1)*dM - 
     &    1._dp/2._dp*m**2*wmass**(-2)*sw2**(-1) - 4._dp/3._dp*log(r13*
     &    r23**(-1))*cw2**(-1) + 16._dp/9._dp*log(r13*r23**(-1))*sw2*
     &    cw2**(-1) )
     &    + DL * ( 1._dp/6._dp - 1._dp/2._dp*sw2**(-1) - 1._dp/4._dp*
     &    sw2**(-1)*cw2 - 17._dp/36._dp*sw2*cw2**(-1) ))

      qqb(-2,1,-1) = 
     &    ( + SL * (  - 1._dp/2._dp + 3._dp/2._dp*sw2**(-1) + 3._dp/4._dp
     &    *sw2**(-1)*cw2 + 1._dp/2._dp*sw2**(-1)*cw2*dM - 1._dp/3._dp*dM + 
     &    17._dp/12._dp*sw2*cw2**(-1) + 17._dp/18._dp*sw2*cw2**(-1)*dM - 
     &    1._dp/4._dp*m**2*wmass**(-2)*sw2**(-1) - 4._dp/3._dp*log(r13*
     &    r23**(-1))*cw2**(-1) + 16._dp/9._dp*log(r13*r23**(-1))*sw2*
     &    cw2**(-1) )
     &    + DL * ( 1._dp/6._dp - 1._dp/2._dp*sw2**(-1) - 1._dp/4._dp*
     &    sw2**(-1)*cw2 - 17._dp/36._dp*sw2*cw2**(-1) ))

      qqb(-2,1,1) = 
     &    ( + SL * ( 8._dp/3._dp*sw2*cw2**(-1) + 16._dp/9._dp*sw2*
     &    cw2**(-1)*dM - 1._dp/2._dp*m**2*wmass**(-2)*sw2**(-1) + 
     &    16._dp/9._dp*log(r13*r23**(-1))*sw2*cw2**(-1) )
     &    + DL * (  - 8._dp/9._dp*sw2*cw2**(-1) ))

      qqb(-1,-1,-1) = 
     &    ( + SL * ( 3._dp*sw2**(-1) + 3._dp/2._dp*sw2**(-1)*cw2 + 
     &    sw2**(-1)*cw2*dM + 1._dp/6._dp*sw2*cw2**(-1) + 1._dp/9._dp*sw2*
     &    cw2**(-1)*dM - 1._dp/4._dp*m**2*wmass**(-2)*sw2**(-1) - log(
     &    r13*r23**(-1))*sw2**(-1)*cw2**(-1) + 2._dp*log(r13*r23**(-1))
     &    *cw2**(-1) - 8._dp/9._dp*log(r13*r23**(-1))*sw2*cw2**(-1) )
     &    + DL * (  - sw2**(-1) - 1._dp/2._dp*sw2**(-1)*cw2 - 1._dp/18._dp
     &    *sw2*cw2**(-1) ))

      qqb(-1,-1,1) = 
     &    ( + SL * ( 1._dp/2._dp + 3._dp/2._dp*sw2**(-1) + 3._dp/4._dp*
     &    sw2**(-1)*cw2 + 1._dp/2._dp*sw2**(-1)*cw2*dM + 1._dp/3._dp*dM + 
     &    17._dp/12._dp*sw2*cw2**(-1) + 17._dp/18._dp*sw2*cw2**(-1)*dM - 
     &    1._dp/2._dp*m**2*wmass**(-2)*sw2**(-1) + 4._dp/3._dp*log(r13*
     &    r23**(-1))*cw2**(-1) - 8._dp/9._dp*log(r13*r23**(-1))*sw2*
     &    cw2**(-1) )
     &    + DL * (  - 1._dp/6._dp - 1._dp/2._dp*sw2**(-1) - 1._dp/4._dp*
     &    sw2**(-1)*cw2 - 17._dp/36._dp*sw2*cw2**(-1) ))

      qqb(-1,1,-1) = 
     &    ( + SL * (  - 1._dp/2._dp + 3._dp/2._dp*sw2**(-1) + 3._dp/4._dp
     &    *sw2**(-1)*cw2 + 1._dp/2._dp*sw2**(-1)*cw2*dM - 1._dp/3._dp*dM + 
     &    5._dp/12._dp*sw2*cw2**(-1) + 5._dp/18._dp*sw2*cw2**(-1)*dM - 1._dp/
     &    4._dp*m**2*wmass**(-2)*sw2**(-1) + 2._dp/3._dp*log(r13*
     &    r23**(-1))*cw2**(-1) - 8._dp/9._dp*log(r13*r23**(-1))*sw2*
     &    cw2**(-1) )
     &    + DL * ( 1._dp/6._dp - 1._dp/2._dp*sw2**(-1) - 1._dp/4._dp*
     &    sw2**(-1)*cw2 - 5._dp/36._dp*sw2*cw2**(-1) ))

      qqb(-1,1,1) = 
     &    ( + SL * ( 5._dp/3._dp*sw2*cw2**(-1) + 10._dp/9._dp*sw2*
     &    cw2**(-1)*dM - 1._dp/2._dp*m**2*wmass**(-2)*sw2**(-1) - 
     &    8._dp/9._dp*log(r13*r23**(-1))*sw2*cw2**(-1) )
     &    + DL * (  - 5._dp/9._dp*sw2*cw2**(-1) ))

      qqb(1,-1,-1) = 
     &    ( + SL * ( 3._dp*sw2**(-1) + 3._dp/2._dp*sw2**(-1)*cw2 + 
     &    sw2**(-1)*cw2*dM + 1._dp/6._dp*sw2*cw2**(-1) + 1._dp/9._dp*sw2*
     &    cw2**(-1)*dM - 1._dp/4._dp*m**2*wmass**(-2)*sw2**(-1) + log(
     &    r13*r23**(-1))*sw2**(-1)*cw2**(-1) - 2._dp*log(r13*r23**(-1))
     &    *cw2**(-1) + 8._dp/9._dp*log(r13*r23**(-1))*sw2*cw2**(-1) )
     &    + DL * (  - sw2**(-1) - 1._dp/2._dp*sw2**(-1)*cw2 - 1._dp/18._dp
     &    *sw2*cw2**(-1) ))

      qqb(1,-1,1) = 
     &    ( + SL * ( 1._dp/2._dp + 3._dp/2._dp*sw2**(-1) + 3._dp/4._dp*
     &    sw2**(-1)*cw2 + 1._dp/2._dp*sw2**(-1)*cw2*dM + 1._dp/3._dp*dM + 
     &    17._dp/12._dp*sw2*cw2**(-1) + 17._dp/18._dp*sw2*cw2**(-1)*dM - 
     &    1._dp/2._dp*m**2*wmass**(-2)*sw2**(-1) - 4._dp/3._dp*log(r13*
     &    r23**(-1))*cw2**(-1) + 8._dp/9._dp*log(r13*r23**(-1))*sw2*
     &    cw2**(-1) )
     &    + DL * (  - 1._dp/6._dp - 1._dp/2._dp*sw2**(-1) - 1._dp/4._dp*
     &    sw2**(-1)*cw2 - 17._dp/36._dp*sw2*cw2**(-1) ))

      qqb(1,1,-1) = 
     &    ( + SL * (  - 1._dp/2._dp + 3._dp/2._dp*sw2**(-1) + 3._dp/4._dp
     &    *sw2**(-1)*cw2 + 1._dp/2._dp*sw2**(-1)*cw2*dM - 1._dp/3._dp*dM + 
     &    5._dp/12._dp*sw2*cw2**(-1) + 5._dp/18._dp*sw2*cw2**(-1)*dM - 1._dp/
     &    4._dp*m**2*wmass**(-2)*sw2**(-1) - 2._dp/3._dp*log(r13*
     &    r23**(-1))*cw2**(-1) + 8._dp/9._dp*log(r13*r23**(-1))*sw2*
     &    cw2**(-1) )
     &    + DL * ( 1._dp/6._dp - 1._dp/2._dp*sw2**(-1) - 1._dp/4._dp*
     &    sw2**(-1)*cw2 - 5._dp/36._dp*sw2*cw2**(-1) ))

      qqb(1,1,1) = 
     &    ( + SL * ( 5._dp/3._dp*sw2*cw2**(-1) + 10._dp/9._dp*sw2*
     &    cw2**(-1)*dM - 1._dp/2._dp*m**2*wmass**(-2)*sw2**(-1) + 
     &    8._dp/9._dp*log(r13*r23**(-1))*sw2*cw2**(-1) )
     &    + DL * (  - 5._dp/9._dp*sw2*cw2**(-1) ))

      qqb(2,-1,-1) = 
     &    ( + SL * (  - 1._dp + 3._dp*sw2**(-1) + 3._dp/2._dp*
     &    sw2**(-1)*cw2 + sw2**(-1)*cw2*dM - 2._dp/3._dp*dM + 1._dp/6._dp*
     &    sw2*cw2**(-1) + 1._dp/9._dp*sw2*cw2**(-1)*dM - 1._dp/4._dp*m**2*
     &    wmass**(-2)*sw2**(-1) - log(r13*r23**(-1))*sw2**(-1)*
     &    cw2**(-1) + 8._dp/3._dp*log(r13*r23**(-1))*cw2**(-1) - 
     &    16._dp/9._dp*log(r13*r23**(-1))*sw2*cw2**(-1) )
     &    + DL * ( 1._dp/3._dp - sw2**(-1) - 1._dp/2._dp*sw2**(-1)*cw2 - 
     &    1._dp/18._dp*sw2*cw2**(-1) ))

      qqb(2,-1,1) = 
     &    ( + SL * (  - 1._dp/2._dp + 3._dp/2._dp*sw2**(-1) + 3._dp/4._dp
     &    *sw2**(-1)*cw2 + 1._dp/2._dp*sw2**(-1)*cw2*dM - 1._dp/3._dp*dM + 
     &    17._dp/12._dp*sw2*cw2**(-1) + 17._dp/18._dp*sw2*cw2**(-1)*dM - 
     &    1._dp/2._dp*m**2*wmass**(-2)*sw2**(-1) + 4._dp/3._dp*log(r13*
     &    r23**(-1))*cw2**(-1) - 16._dp/9._dp*log(r13*r23**(-1))*sw2*
     &    cw2**(-1) )
     &    + DL * ( 1._dp/6._dp - 1._dp/2._dp*sw2**(-1) - 1._dp/4._dp*
     &    sw2**(-1)*cw2 - 17._dp/36._dp*sw2*cw2**(-1) ))

      qqb(2,1,-1) = 
     &    ( + SL * (  - 1._dp/2._dp + 3._dp/2._dp*sw2**(-1) + 3._dp/4._dp
     &    *sw2**(-1)*cw2 + 1._dp/2._dp*sw2**(-1)*cw2*dM - 1._dp/3._dp*dM + 
     &    17._dp/12._dp*sw2*cw2**(-1) + 17._dp/18._dp*sw2*cw2**(-1)*dM - 
     &    1._dp/4._dp*m**2*wmass**(-2)*sw2**(-1) + 4._dp/3._dp*log(r13*
     &    r23**(-1))*cw2**(-1) - 16._dp/9._dp*log(r13*r23**(-1))*sw2*
     &    cw2**(-1) )
     &    + DL * ( 1._dp/6._dp - 1._dp/2._dp*sw2**(-1) - 1._dp/4._dp*
     &    sw2**(-1)*cw2 - 17._dp/36._dp*sw2*cw2**(-1) ))

      qqb(2,1,1) = 
     &    ( + SL * ( 8._dp/3._dp*sw2*cw2**(-1) + 16._dp/9._dp*sw2*
     &    cw2**(-1)*dM - 1._dp/2._dp*m**2*wmass**(-2)*sw2**(-1) - 
     &    16._dp/9._dp*log(r13*r23**(-1))*sw2*cw2**(-1) )
     &    + DL * (  - 8._dp/9._dp*sw2*cw2**(-1) ))

      qqb(5,-1,-1) = 
     &    ( + SL * ( 3._dp*sw2**(-1) + 3._dp/2._dp*sw2**(-1)*cw2 + 
     &    sw2**(-1)*cw2*dM + 1._dp/6._dp*sw2*cw2**(-1) + 1._dp/9._dp*sw2*
     &    cw2**(-1)*dM - 1._dp/2._dp*m**2*wmass**(-2)*sw2**(-1) + log(
     &    r13*r23**(-1))*sw2**(-1)*cw2**(-1) - 2._dp*log(r13*r23**(-1))
     &    *cw2**(-1) + 8._dp/9._dp*log(r13*r23**(-1))*sw2*cw2**(-1) )
     &    + DL * (  - sw2**(-1) - 1._dp/2._dp*sw2**(-1)*cw2 - 1._dp/18._dp
     &    *sw2*cw2**(-1) ))

      qqb(5,-1,1) = 
     &    ( + SL * ( 1._dp/2._dp + 3._dp/2._dp*sw2**(-1) + 3._dp/4._dp*
     &    sw2**(-1)*cw2 + 1._dp/2._dp*sw2**(-1)*cw2*dM + 1._dp/3._dp*dM + 
     &    17._dp/12._dp*sw2*cw2**(-1) + 17._dp/18._dp*sw2*cw2**(-1)*dM - 
     &    3._dp/4._dp*m**2*wmass**(-2)*sw2**(-1) - 4._dp/3._dp*log(r13*
     &    r23**(-1))*cw2**(-1) + 8._dp/9._dp*log(r13*r23**(-1))*sw2*
     &    cw2**(-1) )
     &    + DL * (  - 1._dp/6._dp - 1._dp/2._dp*sw2**(-1) - 1._dp/4._dp*
     &    sw2**(-1)*cw2 - 17._dp/36._dp*sw2*cw2**(-1) ))

      qqb(5,1,-1) = 
     &    ( + SL * (  - 1._dp/2._dp + 3._dp/2._dp*sw2**(-1) + 3._dp/4._dp
     &    *sw2**(-1)*cw2 + 1._dp/2._dp*sw2**(-1)*cw2*dM - 1._dp/3._dp*dM + 
     &    5._dp/12._dp*sw2*cw2**(-1) + 5._dp/18._dp*sw2*cw2**(-1)*dM - 1._dp/
     &    4._dp*m**2*wmass**(-2)*sw2**(-1) - 2._dp/3._dp*log(r13*
     &    r23**(-1))*cw2**(-1) + 8._dp/9._dp*log(r13*r23**(-1))*sw2*
     &    cw2**(-1) )
     &    + DL * ( 1._dp/6._dp - 1._dp/2._dp*sw2**(-1) - 1._dp/4._dp*
     &    sw2**(-1)*cw2 - 5._dp/36._dp*sw2*cw2**(-1) ))

      qqb(5,1,1) = 
     &    ( + SL * ( 5._dp/3._dp*sw2*cw2**(-1) + 10._dp/9._dp*sw2*
     &    cw2**(-1)*dM - 1._dp/2._dp*m**2*wmass**(-2)*sw2**(-1) + 
     &    8._dp/9._dp*log(r13*r23**(-1))*sw2*cw2**(-1) )
     &    + DL * (  - 5._dp/9._dp*sw2*cw2**(-1) ))

c--- qb + q: exchange of t and u
      uu=s(1,3)
      tt=s(2,3)

      msq(-5,5) = 
     & (4._dp*qqb(-5,-1,-1)*uu**2 + 4._dp*qqb(-5,-1,-1)*m**2*ss + 4._dp
     & *qqb(-5,-1,1)*tt**2 + 4._dp*qqb(-5,-1,1)*m**2*ss + 4._dp*qqb(-5,1,
     & -1)*tt**2 + 4._dp*qqb(-5,1,-1)*m**2*ss + 4._dp*qqb(-5,1,1)*uu**2
     &  + 4._dp*qqb(-5,1,1)*m**2*ss)

      msq(-2,2) = 
     & (4._dp*qqb(-2,-1,-1)*uu**2 + 4._dp*qqb(-2,-1,-1)*m**2*ss + 4._dp
     & *qqb(-2,-1,1)*tt**2 + 4._dp*qqb(-2,-1,1)*m**2*ss + 4._dp*qqb(-2,1,
     & -1)*tt**2 + 4._dp*qqb(-2,1,-1)*m**2*ss + 4._dp*qqb(-2,1,1)*uu**2
     &  + 4._dp*qqb(-2,1,1)*m**2*ss)

      msq(-1,1) = 
     & (4._dp*qqb(-1,-1,-1)*uu**2 + 4._dp*qqb(-1,-1,-1)*m**2*ss + 4._dp
     & *qqb(-1,-1,1)*tt**2 + 4._dp*qqb(-1,-1,1)*m**2*ss + 4._dp*qqb(-1,1,
     & -1)*tt**2 + 4._dp*qqb(-1,1,-1)*m**2*ss + 4._dp*qqb(-1,1,1)*uu**2
     &  + 4._dp*qqb(-1,1,1)*m**2*ss)

c--- q + qb
      tt=s(1,3)
      uu=s(2,3)
      it=1._dp/tt
      iu=1._dp/uu

      msq(1,-1) = 
     & (4._dp*qqb(1,-1,-1)*uu**2 + 4._dp*qqb(1,-1,-1)*m**2*ss + 4._dp*
     & qqb(1,-1,1)*tt**2 + 4._dp*qqb(1,-1,1)*m**2*ss + 4._dp*qqb(1,1,-1)*
     & tt**2 + 4._dp*qqb(1,1,-1)*m**2*ss + 4._dp*qqb(1,1,1)*uu**2 + 4._dp*
     & qqb(1,1,1)*m**2*ss)

      msq(2,-2) = 
     & (4._dp*qqb(2,-1,-1)*uu**2 + 4._dp*qqb(2,-1,-1)*m**2*ss + 4._dp*
     & qqb(2,-1,1)*tt**2 + 4._dp*qqb(2,-1,1)*m**2*ss + 4._dp*qqb(2,1,-1)*
     & tt**2 + 4._dp*qqb(2,1,-1)*m**2*ss + 4._dp*qqb(2,1,1)*uu**2 + 4._dp*
     & qqb(2,1,1)*m**2*ss)

      msq(5,-5) = 
     & (4._dp*qqb(5,-1,-1)*uu**2 + 4._dp*qqb(5,-1,-1)*m**2*ss + 4._dp*
     & qqb(5,-1,1)*tt**2 + 4._dp*qqb(5,-1,1)*m**2*ss + 4._dp*qqb(5,1,-1)*
     & tt**2 + 4._dp*qqb(5,1,-1)*m**2*ss + 4._dp*qqb(5,1,1)*uu**2 + 4._dp*
     & qqb(5,1,1)*m**2*ss)

      msq(-3,3)=msq(-1,1)
      msq(3,-3)=msq(1,-1)
      msq(-4,4)=msq(-2,2)
      msq(4,-4)=msq(2,-2)

c--- 1st factor of 2 corresponds to 2*(dMs*Mb)
c--- 2/9/4 color structure and spin ave
      msq=2._dp*msq*gsq**2*2._dp/9._dp/4._dp/ss**2

c--- ttb(-1:1,-1:1) is the EW Sudakov logs to process gg -> ttb
c--- arguments: 1st = helicity of top; 2nd = helicity of anti-top
      ttb(-1,-1) = 
     &    ( + SL * (  - 1._dp/2._dp - 1._dp/4._dp*sw2**(-1)*wmass**(-2)*
     &    m**2 + 3._dp/2._dp*sw2**(-1) + 3._dp/4._dp*sw2**(-1)*cw2 + 
     &    1._dp/2._dp*sw2**(-1)*cw2*dM - 1._dp/3._dp*dM + 
     &    1._dp/12._dp*sw2*cw2**(-1) + 1._dp/18._dp*sw2*cw2**(-1)*dM )
     &    + DL * ( 1._dp/6._dp - 1._dp/2._dp*sw2**(-1) - 1._dp/4._dp*
     &    sw2**(-1)*cw2 - 1._dp/36._dp*sw2*cw2**(-1) ))

      ttb(-1,1) = 
     &    ( + SL * (  - 1._dp/4._dp - 3._dp/8._dp*sw2**(-1)*wmass**(-2)*
     &    m**2 + 3._dp/4._dp*sw2**(-1) + 3._dp/8._dp*sw2**(-1)*cw2 + 
     &    1._dp/4._dp*sw2**(-1)*cw2*dM - 1._dp/6._dp*dM + 17._dp/24._dp*sw2*
     &    cw2**(-1) + 17._dp/36._dp*sw2*cw2**(-1)*dM )
     &    + DL * ( 1._dp/12._dp - 1._dp/4._dp*sw2**(-1) - 1._dp/8._dp*
     &    sw2**(-1)*cw2 - 17._dp/72._dp*sw2*cw2**(-1) ))

      ttb(1,-1) = 
     &    ( + SL * (  - 1._dp/4._dp - 3._dp/8._dp*sw2**(-1)*wmass**(-2)*
     &    m**2 + 3._dp/4._dp*sw2**(-1) + 3._dp/8._dp*sw2**(-1)*cw2 + 
     &    1._dp/4._dp*sw2**(-1)*cw2*dM - 1._dp/6._dp*dM + 17._dp/24._dp*sw2*
     &    cw2**(-1) + 17._dp/36._dp*sw2*cw2**(-1)*dM )
     &    + DL * ( 1._dp/12._dp - 1._dp/4._dp*sw2**(-1) - 1._dp/8._dp*
     &    sw2**(-1)*cw2 - 17._dp/72._dp*sw2*cw2**(-1) ))

      ttb(1,1) = 
     &    ( + SL * (  - 1._dp/2._dp*sw2**(-1)*wmass**(-2)*m**2 + 
     &    4._dp/3._dp*sw2*cw2**(-1) + 8._dp/9._dp*sw2*cw2**(-1)*dM )
     &    + DL * (  - 4._dp/9._dp*sw2*cw2**(-1) ))

c--- broken down (dMs*Mb) w.r.t. color structures of c1 and c2
c--- c1 = (-iTaTb)*(-iTaTb)^*; c2 = (-iTaTb)*(iTbTa)^*
      c1 = 
     &    ( + ttb(-1,-1) * (  - 4._dp + 6._dp*iu*tt + 6._dp*it*uu - 
     &    2._dp*is*uu- 2._dp*is*tt + 2._dp*is*iu*tt**2 + 2._dp*is*it*uu**2
     &    - 8._dp*is**2*uu**2 - 8._dp*is**2*tt**2 - 12._dp*m**2*iu + 
     &    4._dp*m**2*iu**2*tt - 12._dp*m**2*it + 4._dp*m**2*it**2*uu - 
     &    32._dp*m**2*is + 8._dp*m**2*is*iu*tt + 8._dp*m**2*is*it*uu - 
     &    8._dp*m**4*iu**2 - 8._dp*m**4*it**2 )
     &    + ttb(-1,1) * ( 12._dp*m**2*iu - 4._dp*m**2*iu**2*tt + 12._dp*
     &    m**2*it - 4._dp*m**2*it**2*uu + 32._dp*m**2*is + 8._dp*m**2*is*
     &    iu*tt + 8._dp*m**2*is*it*uu - 8._dp*m**4*iu**2 - 8._dp*m**4*
     &    it**2 )
     &    + ttb(1,-1) * ( 12._dp*m**2*iu - 4._dp*m**2*iu**2*tt + 12._dp*
     &    m**2*it - 4._dp*m**2*it**2*uu + 32._dp*m**2*is + 8._dp*m**2*is*
     &    iu*tt + 8._dp*m**2*is*it*uu - 8._dp*m**4*iu**2 - 8._dp*m**4*
     &    it**2 )
     &    + ttb(1,1) * (  - 4._dp + 6._dp*iu*tt + 6._dp*it*uu - 2._dp*is*
     &    uu - 2._dp*is*tt + 2._dp*is*iu*tt**2 + 2._dp*is*it*uu**2 - 8._dp*
     &    is**2*uu**2 - 8._dp*is**2*tt**2 - 12._dp*m**2*iu + 4._dp*m**2*
     &    iu**2*tt - 12._dp*m**2*it + 4._dp*m**2*it**2*uu - 32._dp*m**2*is
     &    + 8._dp*m**2*is*iu*tt + 8._dp*m**2*is*it*uu - 8._dp*m**4*iu**2
     &    - 8._dp*m**4*it**2 ))

      c2 = 
     &    ( + ttb(-1,-1) * (  - 4._dp + 2._dp*iu*tt + 2._dp*it*uu - 
     &    2._dp*is*uu - 2._dp*is*tt + 2._dp*is*iu*tt**2 + 2._dp*is*it*uu**2 
     &    - 8._dp*is**2*uu**2 - 8._dp*is**2*tt**2 - 32._dp*m**2*is + 
     &    8._dp*m**2*is*iu*tt + 8._dp*m**2*is*it*uu + 16._dp*m**4*it*iu )
     &     + ttb(-1,1) * ( 16._dp*m**2*iu + 16._dp*m**2*it + 32._dp*m**2*
     &    is + 8._dp*m**2*is*iu*tt + 8._dp*m**2*is*it*uu + 16._dp*m**4*it*
     &    iu )
     &     + ttb(1,-1) * ( 16._dp*m**2*iu + 16._dp*m**2*it + 32._dp*m**2*
     &    is + 8._dp*m**2*is*iu*tt + 8._dp*m**2*is*it*uu + 16._dp*m**4*it*
     &    iu )
     &    + ttb(1,1) * (  - 4._dp + 2._dp*iu*tt + 2._dp*it*uu - 2._dp*is*
     &    uu - 2._dp*is*tt + 2._dp*is*iu*tt**2 + 2._dp*is*it*uu**2 - 8._dp*
     &    is**2*uu**2 - 8._dp*is**2*tt**2 - 32._dp*m**2*is + 8._dp*m**2*is
     &    *iu*tt + 8._dp*m**2*is*it*uu + 16._dp*m**4*it*iu ))

      msq(0,0)=1.6d1/3._dp*c1+2._dp/3._dp*c2
      msq(0,0)=2._dp*msq(0,0)*avegg*gsq**2

      end

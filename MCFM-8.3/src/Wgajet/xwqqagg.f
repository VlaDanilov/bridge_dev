      subroutine xwqqagg_qq(j1,j2,j3,j4,j5,j6,j7,xa,xb,a70qq)
      implicit none
      include 'types.f'
      
************************************************************************
*     return the helicity amplitudes for                               *
*     0 -> q(p1)+qb(p2)+l(p3)+vb(p4)+gam(p5)+g(p6)+g(p7)               *
*     Qu and Qd parts of the amplitude                                 *
************************************************************************

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zerowidth.f'
      integer:: st1,st2,st3,st4,st5,st6,st7,st8
      integer:: j1,j2,j3,j4,j5,j6,j7,ic,h5,h6,h7
      complex(dp):: a70qq(2,2,2,2)
      complex(dp):: prp34,amp_wqqagg
      complex(dp):: xa(mxpart,mxpart),xb(mxpart,mxpart)

c-----initialize helicity amplitudes to zero
      do h7=1,2
      do h6=1,2
      do h5=1,2
      do ic=1,2
      a70qq(ic,h5,h6,h7)=czip
      enddo
      enddo
      enddo
      enddo

c---- notation is amp_wqqagg(st,pq,pqb,pl,pvb,pga,pg,pg,xa,xb)

c---- helicity stamps (st)
c---- Qd = photon on leg 2 = [a] in Dixon,Signer [hep-ph/9803250] (e.g. see (4.9)-(4.12) for NLO)
c---- Qu = photon on leg 1 = [b] in Dixon,Signer [hep-ph/9803250] (e.g. see (4.9)-(4.12) for NLO)

c---- 'q-,qb+,l-,vb+,ga+,g+,g+,Qd'=1
c---- 'q-,qb+,l-,vb+,ga+,g+,g-,Qd'=2
c---- 'q-,qb+,l-,vb+,ga+,g-,g+,Qd'=3
c---- 'q-,qb+,l-,vb+,ga-,g+,g+,Qd'=4

c---- 'q-,qb+,l-,vb+,ga+,g+,g+,Qu'=5
c---- 'q-,qb+,l-,vb+,ga+,g+,g-,Qu'=6
c---- 'q-,qb+,l-,vb+,ga+,g-,g+,Qu'=7
c---- 'q-,qb+,l-,vb+,ga-,g+,g+,Qu'=8

c---- helicity stamps
      st1=1
      st2=2
      st3=3
      st4=4
      st5=5
      st6=6
      st7=7
      st8=8

      prp34=s(j3,j4)/cplx2(s(j3,j4)-wmass**2,wmass*wwidth)

      a70qq(1,2,2,2) = (Qd*amp_wqqagg(st1,j1,j2,j3,j4,j5,j6,j7,xa,xb)
     &                 +Qu*amp_wqqagg(st5,j1,j2,j3,j4,j5,j6,j7,xa,xb))*prp34
      a70qq(1,2,2,1) = (Qd*amp_wqqagg(st2,j1,j2,j3,j4,j5,j6,j7,xa,xb)
     &                 +Qu*amp_wqqagg(st6,j1,j2,j3,j4,j5,j6,j7,xa,xb))*prp34
      a70qq(1,2,1,2) = (Qd*amp_wqqagg(st3,j1,j2,j3,j4,j5,j6,j7,xa,xb)
     &                 +Qu*amp_wqqagg(st7,j1,j2,j3,j4,j5,j6,j7,xa,xb))*prp34
      a70qq(1,1,2,2) = (Qd*amp_wqqagg(st4,j1,j2,j3,j4,j5,j6,j7,xa,xb)
     &                 +Qu*amp_wqqagg(st8,j1,j2,j3,j4,j5,j6,j7,xa,xb))*prp34

      a70qq(1,1,1,1) = (Qd*amp_wqqagg(st5,j2,j1,j4,j3,j5,j7,j6,xb,xa)
     &                 +Qu*amp_wqqagg(st1,j2,j1,j4,j3,j5,j7,j6,xb,xa))*prp34
      a70qq(1,1,1,2) = (Qd*amp_wqqagg(st7,j2,j1,j4,j3,j5,j7,j6,xb,xa)
     &                 +Qu*amp_wqqagg(st3,j2,j1,j4,j3,j5,j7,j6,xb,xa))*prp34
      a70qq(1,1,2,1) = (Qd*amp_wqqagg(st6,j2,j1,j4,j3,j5,j7,j6,xb,xa)
     &                 +Qu*amp_wqqagg(st2,j2,j1,j4,j3,j5,j7,j6,xb,xa))*prp34
      a70qq(1,2,1,1) = (Qd*amp_wqqagg(st8,j2,j1,j4,j3,j5,j7,j6,xb,xa)
     &                 +Qu*amp_wqqagg(st4,j2,j1,j4,j3,j5,j7,j6,xb,xa))*prp34

      a70qq(2,2,2,2) = (Qd*amp_wqqagg(st1,j1,j2,j3,j4,j5,j7,j6,xa,xb)
     &                 +Qu*amp_wqqagg(st5,j1,j2,j3,j4,j5,j7,j6,xa,xb))*prp34
      a70qq(2,2,2,1) = (Qd*amp_wqqagg(st3,j1,j2,j3,j4,j5,j7,j6,xa,xb)
     &                 +Qu*amp_wqqagg(st7,j1,j2,j3,j4,j5,j7,j6,xa,xb))*prp34
      a70qq(2,2,1,2) = (Qd*amp_wqqagg(st2,j1,j2,j3,j4,j5,j7,j6,xa,xb)
     &                 +Qu*amp_wqqagg(st6,j1,j2,j3,j4,j5,j7,j6,xa,xb))*prp34
      a70qq(2,1,2,2) = (Qd*amp_wqqagg(st4,j1,j2,j3,j4,j5,j7,j6,xa,xb)
     &                 +Qu*amp_wqqagg(st8,j1,j2,j3,j4,j5,j7,j6,xa,xb))*prp34

      a70qq(2,1,1,1) = (Qd*amp_wqqagg(st5,j2,j1,j4,j3,j5,j6,j7,xb,xa)
     &                 +Qu*amp_wqqagg(st1,j2,j1,j4,j3,j5,j6,j7,xb,xa))*prp34
      a70qq(2,1,1,2) = (Qd*amp_wqqagg(st6,j2,j1,j4,j3,j5,j6,j7,xb,xa)
     &                 +Qu*amp_wqqagg(st2,j2,j1,j4,j3,j5,j6,j7,xb,xa))*prp34
      a70qq(2,1,2,1) = (Qd*amp_wqqagg(st7,j2,j1,j4,j3,j5,j6,j7,xb,xa)
     &                 +Qu*amp_wqqagg(st3,j2,j1,j4,j3,j5,j6,j7,xb,xa))*prp34
      a70qq(2,2,1,1) = (Qd*amp_wqqagg(st8,j2,j1,j4,j3,j5,j6,j7,xb,xa)
     &                 +Qu*amp_wqqagg(st4,j2,j1,j4,j3,j5,j6,j7,xb,xa))*prp34

      return
      end


      subroutine xwqqagg_ql(j1,j2,j3,j4,j5,j6,j7,xa,xb,a70ql)
      implicit none
      include 'types.f'

************************************************************************
*     return the helicity amplitudes for                               *
*     0 -> q(p1)+qb(p2)+l(p3)+vb(p4)+gam(p5)+g(p6)+g(p7)               *
*     Ql part of the amplitude                                         *
************************************************************************

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zerowidth.f'
      integer:: st9,st10,st11,st12,st13,st14,st15,st16
      integer:: j1,j2,j3,j4,j5,j6,j7,ic,h5,h6,h7
      real(dp):: s345
      complex(dp):: a70ql(2,2,2,2)
      complex(dp):: prp345,amp_wqqagg
      complex(dp):: xa(mxpart,mxpart),xb(mxpart,mxpart)

c-----initialize helicity amplitudes to zero
      do h7=1,2
      do h6=1,2
      do h5=1,2
      do ic=1,2
      a70ql(ic,h5,h6,h7)=czip
      enddo
      enddo
      enddo
      enddo

c---- notation is amp_wqqagg(st,pq,pqb,pl,pvb,pga,pg,pg,xa,xb)

c---- helicity stamps (st)
c---- Ql = photon emitted by lepton (amplitudes include part of the Wga vertex diagrams)

c---- 'q-,qb+,l-,vb+,ga+,g+,g+,Ql'=9
c---- 'q-,qb+,l-,vb+,ga+,g+,g-,Ql'=10
c---- 'q-,qb+,l-,vb+,ga+,g-,g+,Ql'=11
c---- 'q-,qb+,l-,vb+,ga-,g+,g+,Ql'=12
c---- 'q-,qb+,l-,vb+,ga-,g-,g-,Ql'=13
c---- 'q-,qb+,l-,vb+,ga-,g-,g+,Ql'=14
c---- 'q-,qb+,l-,vb+,ga-,g+,g-,Ql'=15
c---- 'q-,qb+,l-,vb+,ga+,g-,g-,Ql'=16

c---- helicity stamps
      st9=9
      st10=10
      st11=11
      st12=12
      st13=13
      st14=14
      st15=15
      st16=16

      s345=s(j3,j4)+s(j3,j5)+s(j4,j5)

      if (zerowidth) then
c---- zerowidth: no final state radiation, so we can set prp345 to zero
        prp345=czip
      else
c---- otherwise, usual Breit-Wigner form
        prp345=s345/cplx2(s345-wmass**2,wmass*wwidth)
      endif

      a70ql(1,2,2,2) =  (Qd-Qu)*prp345*amp_wqqagg(st9,j1,j2,j3,j4,j5,j6,j7,xa,xb)
      a70ql(1,2,2,1) = (Qd-Qu)*prp345*amp_wqqagg(st10,j1,j2,j3,j4,j5,j6,j7,xa,xb)
      a70ql(1,2,1,2) = (Qd-Qu)*prp345*amp_wqqagg(st11,j1,j2,j3,j4,j5,j6,j7,xa,xb)
      a70ql(1,1,2,2) = (Qd-Qu)*prp345*amp_wqqagg(st12,j1,j2,j3,j4,j5,j6,j7,xa,xb)
      a70ql(1,1,1,1) = (Qd-Qu)*prp345*amp_wqqagg(st13,j1,j2,j3,j4,j5,j6,j7,xa,xb)
      a70ql(1,1,1,2) = (Qd-Qu)*prp345*amp_wqqagg(st14,j1,j2,j3,j4,j5,j6,j7,xa,xb)
      a70ql(1,1,2,1) = (Qd-Qu)*prp345*amp_wqqagg(st15,j1,j2,j3,j4,j5,j6,j7,xa,xb)
      a70ql(1,2,1,1) = (Qd-Qu)*prp345*amp_wqqagg(st16,j1,j2,j3,j4,j5,j6,j7,xa,xb)

      a70ql(2,2,2,2) =  (Qd-Qu)*prp345*amp_wqqagg(st9,j1,j2,j3,j4,j5,j7,j6,xa,xb)
      a70ql(2,2,2,1) = (Qd-Qu)*prp345*amp_wqqagg(st11,j1,j2,j3,j4,j5,j7,j6,xa,xb)
      a70ql(2,2,1,2) = (Qd-Qu)*prp345*amp_wqqagg(st10,j1,j2,j3,j4,j5,j7,j6,xa,xb)
      a70ql(2,1,2,2) = (Qd-Qu)*prp345*amp_wqqagg(st12,j1,j2,j3,j4,j5,j7,j6,xa,xb)
      a70ql(2,1,1,1) = (Qd-Qu)*prp345*amp_wqqagg(st13,j1,j2,j3,j4,j5,j7,j6,xa,xb)
      a70ql(2,1,1,2) = (Qd-Qu)*prp345*amp_wqqagg(st15,j1,j2,j3,j4,j5,j7,j6,xa,xb)
      a70ql(2,1,2,1) = (Qd-Qu)*prp345*amp_wqqagg(st14,j1,j2,j3,j4,j5,j7,j6,xa,xb)
      a70ql(2,2,1,1) = (Qd-Qu)*prp345*amp_wqqagg(st16,j1,j2,j3,j4,j5,j7,j6,xa,xb)

      return
      end
 

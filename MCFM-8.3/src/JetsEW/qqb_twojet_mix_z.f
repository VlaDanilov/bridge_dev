      subroutine qqb_twojet_mix_z(p,z)
c--- Integrated subtractions, averaged over initial colors and spins,
c--- for mixed QCD/Z matrix element for dijet production,
c     q(-p1)+qbar(-p2) --> q(p3)+q~(p4) + g(p5)
c
c--- all momenta are incoming
c
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'constants.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'agq.f'
      include 'PR_mix.f'
      include 'breit.f'
      real(dp):: z,p(mxpart,4),dot,metric,Q34sq,
     & xl12,xl13,xl14,xl23,xl24,xl34,tempgq,tempqg,
     . ii_qq,ii_gg,if_qq,if_gg,fi_qq,ff_qq,ii_qg,ii_gq
      integer is,nu,icol

      xl12=log(+2._dp*dot(p,1,2)/musq)
      xl13=log(-2._dp*dot(p,1,3)/musq)
      xl14=log(-2._dp*dot(p,1,4)/musq)
      xl23=log(-2._dp*dot(p,2,3)/musq)
      xl24=log(-2._dp*dot(p,2,4)/musq)
      xl34=log(+2._dp*dot(p,3,4)/musq)

c--- The variables R and P provide the Regular and Plus pieces associated
c--- with radiation from leg 1 (M1(a,b,c,is) and leg 2 (M2(a,b,c,is)
c--- In each case the parton labelling is Using the normal QM notation of putting
c--- everything backward
c---       emitted line after emission =   a
c---       emitter before emission     =   b
c---       spectator                   =   c
c--- There is no label for he or she who is emitted.
c--- Note that in general each piece will be composed of many different
c--- dipole contributions

c--- note that, in order to distinguish identical and non-identical
c--- contributions in quark-antiquark processes, the integer isame (=4)
c--- is used to fill higher elements of arrays


      do is=1,3
c--- quark-quark
c---   color-structure = 0
      M1(q,q,q,0,is)=ason4pi*(
     & +if_qq(z,xl14,is)+fi_qq(z,xl14,is)
     & -ii_qq(z,xl12,is)-ff_qq(z,xl34,is))
     
      M2(q,q,q,0,is)=ason4pi*(
     & +if_qq(z,xl23,is)+fi_qq(z,xl23,is)
     & -ii_qq(z,xl12,is)-ff_qq(z,xl34,is))

c---   color-structure = 1
      M1(q,q,q,1,is)=ason4pi/xn*(
     & +2._dp*Cf*xn*(if_qq(z,xl13,is)+fi_qq(z,xl13,is))
     & -if_qq(z,xl14,is)-fi_qq(z,xl14,is)
     & +ii_qq(z,xl12,is)+ff_qq(z,xl34,is))
     
      M2(q,q,q,1,is)=ason4pi/xn*(
     & +2._dp*Cf*xn*(if_qq(z,xl24,is)+fi_qq(z,xl24,is))
     & -if_qq(z,xl23,is)-fi_qq(z,xl23,is)
     & +ii_qq(z,xl12,is)+ff_qq(z,xl34,is))

c---   color-structure = 2
      M1(q,q,q,2,is)=ason4pi/xn*(
     & +2._dp*Cf*xn*(if_qq(z,xl14,is)+fi_qq(z,xl14,is))
     & -if_qq(z,xl13,is)-fi_qq(z,xl13,is)
     & +ii_qq(z,xl12,is)+ff_qq(z,xl34,is))
     
      M2(q,q,q,2,is)=ason4pi/xn*(
     & +2._dp*Cf*xn*(if_qq(z,xl23,is)+fi_qq(z,xl23,is))
     & -if_qq(z,xl24,is)-fi_qq(z,xl24,is)
     & +ii_qq(z,xl12,is)+ff_qq(z,xl34,is))

c---   color-structure = 3
      M1(q,q,q,3,is)=ason4pi*(
     & +if_qq(z,xl13,is)+fi_qq(z,xl13,is)
     & -ii_qq(z,xl12,is)-ff_qq(z,xl34,is))
     
      M2(q,q,q,3,is)=ason4pi*(
     & +if_qq(z,xl24,is)+fi_qq(z,xl24,is)
     & -ii_qq(z,xl12,is)-ff_qq(z,xl34,is))

c--- quark-antiquark: non-identical quarks
c---   color-structure = 0
      M1(q,q,a,0,is)=ason4pi*(
     & -if_qq(z,xl14,is)-fi_qq(z,xl14,is)
     & +ii_qq(z,xl12,is)+ff_qq(z,xl34,is))
     
      M2(a,a,q,0,is)=ason4pi*(
     & -if_qq(z,xl23,is)-fi_qq(z,xl23,is)
     & +ii_qq(z,xl12,is)+ff_qq(z,xl34,is))

c---   color-structure = 1
      M1(q,q,a,1,is)=ason4pi/xn*(
     & +2._dp*Cf*xn*(if_qq(z,xl13,is)+fi_qq(z,xl13,is))
     & +if_qq(z,xl14,is)+fi_qq(z,xl14,is)
     & -ii_qq(z,xl12,is)-ff_qq(z,xl34,is))
     
      M2(a,a,q,1,is)=ason4pi/xn*(
     & +2._dp*Cf*xn*(if_qq(z,xl24,is)+fi_qq(z,xl24,is))
     & +if_qq(z,xl23,is)+fi_qq(z,xl23,is)
     & -ii_qq(z,xl12,is)-ff_qq(z,xl34,is))

c---   color-structure = 2
      M1(q,q,a,2,is)=ason4pi/xn*(
     & +2._dp*Cf*xn*(ii_qq(z,xl12,is)+ff_qq(z,xl34,is))
     & -if_qq(z,xl13,is)-fi_qq(z,xl13,is)
     & +if_qq(z,xl14,is)+fi_qq(z,xl14,is))
     
      M2(a,a,q,2,is)=ason4pi/xn*(
     & +2._dp*Cf*xn*(ii_qq(z,xl12,is)+ff_qq(z,xl34,is))
     & -if_qq(z,xl24,is)-fi_qq(z,xl24,is)
     & +if_qq(z,xl23,is)+fi_qq(z,xl23,is))

c---   color-structure = 3
      M1(q,q,a,3,is)=ason4pi*(
     & +if_qq(z,xl13,is)+fi_qq(z,xl13,is)
     & -if_qq(z,xl14,is)-fi_qq(z,xl14,is))
     
      M2(a,a,q,3,is)=ason4pi*(
     & +if_qq(z,xl24,is)+fi_qq(z,xl24,is)
     & -if_qq(z,xl23,is)-fi_qq(z,xl23,is))

c--- quark-antiquark: identical quarks
c---   color-structure = 0
      M1(q,q,a,0+isame,is)=ason4pi*(
     & +if_qq(z,xl13,is)+fi_qq(z,xl13,is)
     & -if_qq(z,xl14,is)-fi_qq(z,xl14,is))
     
      M2(a,a,q,0+isame,is)=ason4pi*(
     & +if_qq(z,xl24,is)+fi_qq(z,xl24,is)
     & -if_qq(z,xl23,is)-fi_qq(z,xl23,is))

c---   color-structure = 1
      M1(q,q,a,1+isame,is)=ason4pi/xn*(
     & +2._dp*Cf*xn*(ii_qq(z,xl12,is)+ff_qq(z,xl34,is))
     & -if_qq(z,xl13,is)-fi_qq(z,xl13,is)
     & +if_qq(z,xl14,is)+fi_qq(z,xl14,is))
     
      M2(a,a,q,1+isame,is)=ason4pi/xn*(
     & +2._dp*Cf*xn*(ii_qq(z,xl12,is)+ff_qq(z,xl34,is))
     & -if_qq(z,xl24,is)-fi_qq(z,xl24,is)
     & +if_qq(z,xl23,is)+fi_qq(z,xl23,is))

c---   color-structure = 2
      M1(q,q,a,2+isame,is)=ason4pi/xn*(
     & +2._dp*Cf*xn*(if_qq(z,xl13,is)+fi_qq(z,xl13,is))
     & +if_qq(z,xl14,is)+fi_qq(z,xl14,is)
     & -ii_qq(z,xl12,is)-ff_qq(z,xl34,is))
     
      M2(a,a,q,2+isame,is)=ason4pi/xn*(
     & +2._dp*Cf*xn*(if_qq(z,xl24,is)+fi_qq(z,xl24,is))
     & +if_qq(z,xl23,is)+fi_qq(z,xl23,is)
     & -ii_qq(z,xl12,is)-ff_qq(z,xl34,is))

c---   color-structure = 3
      M1(q,q,a,3+isame,is)=ason4pi*(
     & -if_qq(z,xl14,is)-fi_qq(z,xl14,is)
     & +ii_qq(z,xl12,is)+ff_qq(z,xl34,is))
     
      M2(a,a,q,3+isame,is)=ason4pi*(
     & -if_qq(z,xl23,is)-fi_qq(z,xl23,is)
     & +ii_qq(z,xl12,is)+ff_qq(z,xl34,is))

c-- copy to charge-conjugate states
      M1(a,a,a,:,is)=M1(q,q,q,:,is)
      M2(a,a,a,:,is)=M2(q,q,q,:,is)
      
      M1(a,a,q,:,is)=M1(q,q,a,:,is)
      M2(q,q,a,:,is)=M2(a,a,q,:,is)

c--- quark-gluon pieces
      tempqg=+ason2pi*tr*ii_qg(z,xl12,is)

      M1(a,g,q,0:3,is)=0._dp
      M1(q,g,a,0:3,is)=0._dp
      M1(a,g,a,0:3,is)=0._dp
      M1(q,g,q,0:3,is)=0._dp
      M1(a,g,q,1:2,is)=tempqg
      M1(q,g,a,1:2,is)=tempqg
      M1(a,g,a,1:2,is)=tempqg
      M1(q,g,q,1:2,is)=tempqg

      M2(a,g,q,0:3,is)=0._dp
      M2(q,g,a,0:3,is)=0._dp
      M2(a,g,a,0:3,is)=0._dp
      M2(q,g,q,0:3,is)=0._dp
      M2(a,g,q,1:2,is)=tempqg
      M2(q,g,a,1:2,is)=tempqg
      M2(a,g,a,1:2,is)=tempqg
      M2(q,g,q,1:2,is)=tempqg
      enddo

      return
      end


      subroutine qqb_QQb_mix_z(p,z)
c--- Integrated subtractions, averaged over initial colors and spins,
c--- for mixed QCD/Z matrix element for top pair production,
c     q(-p1)+qbar(-p2) -->  Q(p3)+Q~(P4)

c--- all momenta are incoming
c
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'constants.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'agq.f'
      include 'PR_new.f'
      include 'breit.f'
      real(dp):: z,p(mxpart,4),dot,metric,Q34sq,
     & xl12,xl13,xl14,xl23,xl24,xl34,
     & mbar12,mbar13,mbar14,mbar23,mbar24,mbar34,tempgq,tempqg,
     .                 ii_mqq,ii_mgg,
     .                 if_mqq,if_mgg,
     .                 fi_mqq,
     .                 ff_mqq,
     .                 ii_qg,ii_gq
      integer is,nu,icol

CDTS (5.45,5.77)
      mbar13=mass2/sqrt(-two*Dot(p,1,3))
      mbar23=mass2/sqrt(-two*Dot(p,2,3))
      mbar14=mass2/sqrt(-two*Dot(p,1,4))
      mbar24=mass2/sqrt(-two*Dot(p,2,4))
        
      xl13=log(-two*Dot(p,1,3)/musq)
      xl14=log(-two*Dot(p,1,4)/musq)
      xl23=log(-two*Dot(p,2,3)/musq)
      xl24=log(-two*Dot(p,2,4)/musq)

c--- The variables R and P provide the Regular and Plus pieces associated
c--- with radiation from leg 1 (Q1(a,b,c,is) and leg 2 (Q2(a,b,c,is)
c--- In each case the parton labelling is Using the normal QM notation of putting
c--- everything backward
c---       emitted line after emission =   a
c---       emitter before emission     =   b
c---       spectator                   =   c
c--- There is no label for he or she who is emitted.
c--- Note that in general each piece will be composed of many different
c--- dipole contributions

      do is=1,3
      Q1(q,q,a,is)=ason4pi*4._dp*(
     & +if_mqq(z,xl13,mbar13,is)+fi_mqq(z,xl13,mbar13,is)
     & -if_mqq(z,xl14,mbar14,is)-fi_mqq(z,xl14,mbar14,is))
      Q2(a,a,q,is)=ason4pi*4._dp*(
     & -if_mqq(z,xl23,mbar23,is)-fi_mqq(z,xl23,mbar23,is)
     & +if_mqq(z,xl24,mbar24,is)+fi_mqq(z,xl24,mbar24,is))

c--- interchanging 1 and 2 in above just flips the sign
      Q1(a,a,q,is)=-Q1(q,q,a,is)
      Q2(q,q,a,is)=-Q2(a,a,q,is)

      enddo

      return
      end


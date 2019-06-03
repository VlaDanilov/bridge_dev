      subroutine gg_hg_v(p,msq)
        implicit none
        include 'types.f'
        include 'mxpart.f'
        include 'masses.f'
        include 'hdecaymode.f'
        include 'nf.f'

        real(dp), intent(in) :: p(mxpart,4)
        real(dp), intent(out) :: msq(-nf:nf,-nf:nf)

        real(dp) :: dotvec, s34
        real(dp) :: hdecay,msqhgamgam
        integer, parameter :: iglue = 5

        s34 = dotvec(p(3,:)+p(4,:),p(3,:)+p(4,:))

        hdecay = 0._dp
        msq = 0._dp

        if (hdecaymode == 'tlta') then
            call htautaudecay(p,3,4,hdecay)
        elseif (hdecaymode == 'bqba') then
            call hbbdecay(p,3,4,hdecay)
        elseif (hdecaymode == 'gaga') then
            hdecay=msqhgamgam(s34)
        else
          write(6,*) 'Unimplemented process in gg_hg_v'
          call abort
        endif

        hdecay=hdecay/((s34-hmass**2)**2+(hmass*hwidth)**2)

        call gg_hg_v_nodecay(p,iglue,msq)

        msq = msq*hdecay

      end subroutine

      subroutine gg_hg_zgam_v(p,msq)
        implicit none
        include 'types.f'
        include 'mxpart.f'
        include 'masses.f'
        include 'hdecaymode.f'
        include 'nf.f'

        real(dp), intent(in) :: p(mxpart,4)
        real(dp), intent(out) :: msq(-nf:nf,-nf:nf)

        real(dp) :: dotvec, hdecay
        real(dp) :: shsq, HZgamMSQ
        integer, parameter :: iglue = 6

        hdecay = 0._dp
        
        shsq = dotvec(p(3,:)+p(4,:)+p(5,:), p(3,:)+p(4,:)+p(5,:))
        hdecay = HZgamMSQ(3,4,5)
        hdecay = hdecay/((shsq-hmass**2)**2+(hmass*hwidth)**2)

        call gg_hg_v_nodecay(p,iglue,msq)

        msq = msq*hdecay

      end subroutine
      
      subroutine gg_hg_v_nodecay(p,iglue,msq)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'scheme.f'
C     (Taken from Ravindran, Smith, van Neerven hep-ph/0201114)
C     Modified by overall factors

      real(dp), intent(in) :: p(mxpart,4)
      integer, intent(in) :: iglue
      real(dp), intent(out) :: msq(-nf:nf,-nf:nf)

      integer:: j,k
      real(dp):: ss,tt,uu,
     & virtgg,virtqa,virtaq,virtqg,virtgq,Asq,fac

      scheme='tH-V'

      call dotem(iglue,p,s)
      ss=s(1,2)
      tt=s(1,iglue)
      uu=s(2,iglue)

      Asq=(as/(3._dp*pi))**2/vevsq

      fac=ason2pi*Asq*gsq
      call hjetfill(ss,tt,uu,virtgg,virtqa,virtaq,virtqg,virtgq)

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      if ((j==0).and.(k==0)) msq(j,k)=avegg*fac*virtgg
      if ((j>0).and.(k==-j)) msq(j,k)=aveqq*fac*virtqa
      if ((j<0).and.(k==-j)) msq(j,k)=aveqq*fac*virtaq
      if ((j==0).and.(k.ne.0)) msq(j,k)=aveqg*fac*virtgq
      if ((j.ne.0).and.(k==0)) msq(j,k)=aveqg*fac*virtqg
      enddo
      enddo

      return
      end

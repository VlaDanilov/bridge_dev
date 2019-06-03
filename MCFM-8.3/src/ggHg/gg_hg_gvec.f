      
      subroutine gg_hg_zgam_gvec(p,n,in,msq)
        implicit none
        include 'types.f'
        include 'mxpart.f'
        include 'nf.f'
        include 'masses.f'

        real(dp), intent(in) :: p(mxpart,4), n(4)
        integer, intent(in) :: in
        real(dp), intent(out) :: msq(-nf:nf,-nf:nf)

        integer, parameter :: iglue = 6

        real(dp) :: shsq, hdecay, HZgamMSQ, dotvec

        hdecay = 0._dp

        shsq = dotvec(p(3,:)+p(4,:)+p(5,:), p(3,:)+p(4,:)+p(5,:))
        hdecay = HZgamMSQ(3,4,5)
        hdecay = hdecay/((shsq-hmass**2)**2+(hmass*hwidth)**2)

        call gg_hg_gvec_nodecay(p,n,in,iglue,msq)

        msq = msq*hdecay

      end subroutine

      subroutine gg_hg_gvec(p,n,in,msq)
        implicit none
        include 'types.f'
        include 'mxpart.f'
        include 'nf.f'
        include 'masses.f'
        include 'hdecaymode.f'

        real(dp), intent(in) :: p(mxpart,4), n(4)
        integer, intent(in) :: in
        real(dp), intent(out) :: msq(-nf:nf,-nf:nf)

        integer, parameter :: iglue = 5
        real(dp) :: hdecay, dotvec, s34, msqhgamgam

        s34 = dotvec(p(3,:)+p(4,:),p(3,:)+p(4,:))

        hdecay = 0._dp

C   De  al with Higgs decay
        if (hdecaymode == 'tlta') then
            call htautaudecay(p,3,4,hdecay)
        elseif (hdecaymode == 'bqba') then
            call hbbdecay(p,3,4,hdecay)
        elseif (hdecaymode == 'gaga') then
             hdecay=msqhgamgam(s34)
        else
          write(6,*) 'Unimplemented process in gg_hgg_gvec'
          call abort
        endif

        hdecay=hdecay/((s34-hmass**2)**2+(hmass*hwidth)**2)

        call gg_hg_gvec_nodecay(p,n,in,iglue,msq)

        msq = msq*hdecay

      end subroutine
      
      subroutine gg_hg_gvec_nodecay(p,n,in,iglue,msq)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'

      real(dp), intent(in) :: p(mxpart,4), n(4)
C  in is the label of the momentum contracted with n
      integer, intent(in) :: in, iglue
      real(dp), intent(out) :: msq(-nf:nf,-nf:nf)

      integer :: j,k
      real(dp) :: qqghn,ggghn,p1p2(-1:1,-1:1)

      msq(:,:)=zip

      p1p2(:,:)=zip

      if (in == 1) then
      p1p2(0,-1)=-aveqg*qqghn(2,iglue,1,p,n)
      p1p2(0,+1)=-aveqg*qqghn(2,iglue,1,p,n)
      p1p2(0,0)=+avegg*ggghn(iglue,2,1,p,n)
      elseif (in == 2) then
      p1p2(+1,0)=-aveqg*qqghn(1,iglue,2,p,n)
      p1p2(-1,0)=-aveqg*qqghn(iglue,1,2,p,n)
      p1p2(0,0)=+avegg*ggghn(1,iglue,2,p,n)
      elseif (in == iglue) then     
      p1p2(1,-1)=+aveqq*qqghn(1,2,iglue,p,n)
      p1p2(-1,1)=+aveqq*qqghn(2,1,iglue,p,n)
      p1p2(0,0)=+avegg*ggghn(1,2,iglue,p,n)
      else
        write(6,*) 'Illegal value of in in gg_hg_gvec: ',in
        stop
      endif

      do j=-nf,nf
      do k=-nf,nf
      if     ((j > 0) .and. (k == -j)) then
          msq(j,k)=p1p2(1,-1)
      elseif ((j < 0) .and. (k == -j)) then
          msq(j,k)=p1p2(-1,1)
      elseif ((j == 0) .and. (k == 0)) then
          msq(j,k)=p1p2(0,0)
      elseif ((j > 0) .and. (k == 0)) then
          msq(j,k)=p1p2(+1,0)
      elseif ((j < 0) .and. (k == 0)) then
          msq(j,k)=p1p2(-1,0)
      elseif ((j == 0) .and. (k > 0)) then
          msq(j,k)=p1p2(0,+1)
      elseif ((j == 0) .and. (k < 0)) then
          msq(j,k)=p1p2(0,-1)
      endif
      enddo
      enddo
 
      return
      end

      function qqghn(j1,j2,j5,p,n)
      implicit none
      include 'types.f'
      real(dp):: qqghn
      
C---calculates the amplitude squared for the process
c   q(p1)+qbar(p2) --> H((p3+p4)+g(p5)
c   contracted with the vector n(mu)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      integer:: j1,j2,j5
      real(dp):: Asq,p(mxpart,4),n(4),nDn,nDp1,nDp2,dot,s,t,u

      nDp1=n(4)*p(j1,4)-n(3)*p(j1,3)-n(2)*p(j1,2)-n(1)*p(j1,1)
      nDp2=n(4)*p(j2,4)-n(3)*p(j2,3)-n(2)*p(j2,2)-n(1)*p(j2,1)
      nDn=n(4)**2-n(3)**2-n(2)**2-n(1)**2

      call checkndotp(p,n,j5)

      Asq=(as/(three*pi))**2/vevsq

      s=two*Dot(p,j1,j2)
      t=two*Dot(p,j1,j5)
      u=two*Dot(p,j2,j5)

c--- RKE Schoonship answer, ncalc.out
c     Id,bit=g^2*A^2*V/2*(2*(nDp1*u-nDp2*t)^2/s^2
c           + 0.5*nDn*(u+t)^2/s)

      qqghn=-Asq*gsq*V/two*(two*(nDp1*u-nDp2*t)**2/s**2
     &                    +0.5_dp*nDn*(u+t)**2/s)

      return
      end

      function ggghn(j1,j2,j5,p,n)
      implicit none
      include 'types.f'
      real(dp):: ggghn
      
C---calculates the amplitude squared for the process
c   g(p1)+g(p2) --> H((p3+p4)+g(p5)
c   contracted with the vector n(mu)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      integer:: j1,j2,j5
      real(dp):: Asq,p(mxpart,4),n(4),nDn,nDp1,nDp2,dot,s,t,u,sh

      nDp1=n(4)*p(j1,4)-n(3)*p(j1,3)-n(2)*p(j1,2)-n(1)*p(j1,1)
      nDp2=n(4)*p(j2,4)-n(3)*p(j2,3)-n(2)*p(j2,2)-n(1)*p(j2,1)
      nDn=n(4)**2-n(3)**2-n(2)**2-n(1)**2

      call checkndotp(p,n,j5)

      Asq=(as/(three*pi))**2/vevsq

      s=two*Dot(p,j1,j2)
      t=two*Dot(p,j1,j5)
      u=two*Dot(p,j2,j5)
      sh=s+t+u

c--- JMC answer, gggH.frm
c -f(b,a,c)^2/s12/s13/s23*(
c -(n.n)/2*(s12^4+s13^4+s23^4+mHsq^4-2*(s13^2*s23^2+s12^2*mHsq^2))
c +2*(p1.n*s23-p2.n*s13)^2/s13/s23*(s13^2*s23^2+s12^2*mHsq^2)/s12);

      ggghn=Asq*gsq*V*xn*(
     & -nDn/two*(s**4+t**4+u**4+sh**4-two*(t**2*u**2+s**2*sh**2))
     & +two*(nDp1*u-nDp2*t)**2/t/u*(t**2*u**2+s**2*sh**2)/s)/s/t/u

      return
      end
      

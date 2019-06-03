      subroutine qqb_zaj_gs_vdecay(p,msq)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ptilde.f'
      include 'nflav.f'
      include 'qqgg.f'
      include 'frag.f'
      include 'ewcharge.f'
      include 'phot_dip.f'
      include 'ipsgen.f'        
      integer:: j,k,nd,f
      real(dp):: p(mxpart,4),msq(maxd,-nf:nf,-nf:nf),nup,ndo
      real(dp):: 
     & msq16_2(-nf:nf,-nf:nf),msq26_1(-nf:nf,-nf:nf),
     & msq17_2(-nf:nf,-nf:nf),msq27_1(-nf:nf,-nf:nf),
     & msq16_7(-nf:nf,-nf:nf),msq27_6(-nf:nf,-nf:nf),
     & msq17_6(-nf:nf,-nf:nf),msq26_7(-nf:nf,-nf:nf),
     & msq67_1v(-nf:nf,-nf:nf),msq67_2v(-nf:nf,-nf:nf),
     & msq27_6v(-nf:nf,-nf:nf),msq27_1v(-nf:nf,-nf:nf),
     & msq16_7v(-nf:nf,-nf:nf),msq17_2v(-nf:nf,-nf:nf),
     & msq17_6v(-nf:nf,-nf:nf),msq26_7v(-nf:nf,-nf:nf),
     & msq26_1v(-nf:nf,-nf:nf),
     & msq16_2v(-nf:nf,-nf:nf),
     & dummy(-nf:nf,-nf:nf),
     & sub16_2(4),sub26_1(4),sub17_2(4),sub27_1(4),
     & sub16_7(4),sub17_6(4),sub26_7(4),sub27_6(4),
     & sub67_1(4),sub67_2(4),sub67_1v,sub67_2v,
     & sub27_6v,sub26_1v,sub27_1v,sub17_6v,sub17_2v,sub16_2v,sub16_7v,
     & sub26_7v
      real(dp):: sub56_1,sub56_2,sub56_7,sub57_6,sub57_1,sub57_2
      real(dp):: msq56_1(-nf:nf,-nf:nf),msq56_2(-nf:nf,-nf:nf),
     &     msq56_7(-nf:nf,-nf:nf),msq57_1(-nf:nf,-nf:nf), 
     &     msq57_2(-nf:nf,-nf:nf),msq57_6(-nf:nf,-nf:nf)  
      real(dp)::
     &  msq56_1x(0:2,-nf:nf,-nf:nf,-nf:nf,-nf:nf),
     &  msq57_1x(0:2,-nf:nf,-nf:nf,-nf:nf,-nf:nf),
     &  msq56_7x(0:2,-nf:nf,-nf:nf,-nf:nf,-nf:nf),
     &  msq57_6x(0:2,-nf:nf,-nf:nf,-nf:nf,-nf:nf)
      real(dp):: 
     &  msq57_1x_swap(0:2,-nf:nf,-nf:nf,-nf:nf,-nf:nf),
     &  msq57_1_swap(-nf:nf,-nf:nf)
      real(dp):: 
     &  msq56_1x_swap(0:2,-nf:nf,-nf:nf,-nf:nf,-nf:nf),
     &  msq56_1_swap(-nf:nf,-nf:nf)
      external qqb_zaj_vdecay
      external qqb_zaj_vdecay_cross
      external qqb_zaj_gvec_vdecay
      external donothing_gvec


      real(dp) :: sub37_4(4), sub37_4v
      real(dp) :: sub37_6(4), sub37_6v

      real(dp) :: sub47_3(4), sub47_3v
      real(dp) :: sub47_6(4), sub47_6v

      real(dp) :: sub47_3_cross(4), sub47_3v_cross
      real(dp) :: sub47_6_cross(4), sub47_6v_cross

      real(dp) :: sub67_3(4), sub67_3v
      real(dp) :: sub67_4(4), sub67_4v
      real(dp) :: msq37_4(-nf:nf,-nf:nf), msq37_4v(-nf:nf,-nf:nf)
      real(dp) :: msq37_6(-nf:nf,-nf:nf), msq37_6v(-nf:nf,-nf:nf)

      real(dp) :: msq47_3(-nf:nf,-nf:nf), msq47_3v(-nf:nf,-nf:nf)
      real(dp) :: msq47_6(-nf:nf,-nf:nf), msq47_6v(-nf:nf,-nf:nf)

      real(dp) :: msq47_3_cross(-nf:nf,-nf:nf), msq47_3v_cross(-nf:nf,-nf:nf)
      real(dp) :: msq47_6_cross(-nf:nf,-nf:nf), msq47_6v_cross(-nf:nf,-nf:nf)

      real(dp) :: msq67_3(-nf:nf,-nf:nf), msq67_3v(-nf:nf,-nf:nf)
      real(dp) :: msq67_4(-nf:nf,-nf:nf), msq67_4v(-nf:nf,-nf:nf)

      ! dip 7
      real(dp) :: sub36_4(4), sub36_4v
      real(dp) :: msq36_4(-nf:nf,-nf:nf), msq36_4v(-nf:nf,-nf:nf)

      ! dip 8
      real(dp) :: sub36_7(4), sub36_7v
      real(dp) :: msq36_7(-nf:nf,-nf:nf), msq36_7v(-nf:nf,-nf:nf)

      ! dip 9
      real(dp) :: sub46_3(4), sub46_3v
      real(dp) :: msq46_3(-nf:nf,-nf:nf), msq46_3v(-nf:nf,-nf:nf)

      ! dip 10
      real(dp) :: sub46_7(4), sub46_7v
      real(dp) :: msq46_7(-nf:nf,-nf:nf), msq46_7v(-nf:nf,-nf:nf)

      msq = 0._dp
      ndmax = 6

c --- vdecay only fifi dipoles

      call dips(1,p,3,7,4,sub37_4,sub37_4v,msq37_4,msq37_4v,
     &  qqb_zaj_vdecay,donothing_gvec)
      call dips(2,p,3,7,6,sub37_6,sub37_6v,msq37_6,msq37_6v,
     &  qqb_zaj_vdecay,donothing_gvec)

      ! for qa gg
      call dips(3,p,4,7,3,sub47_3,sub47_3v,msq47_3,msq47_3v,
     &  qqb_zaj_vdecay,donothing_gvec)
      call dips(4,p,4,7,6,sub47_6,sub47_6v,msq47_6,msq47_6v,
     &  qqb_zaj_vdecay,donothing_gvec)

      ! for qa qa identical quarks
      call dips(3,p,4,7,3,sub47_3_cross,sub47_3v_cross,msq47_3_cross,msq47_3v_cross,
     &  qqb_zaj_vdecay_cross,qqb_zaj_gvec_vdecay)
      call dips(4,p,4,7,6,sub47_6_cross,sub47_6v_cross,msq47_6_cross,msq47_6v_cross,
     &  qqb_zaj_vdecay_cross,qqb_zaj_gvec_vdecay)

      call dips(5,p,6,7,3,sub67_3,sub67_3v,msq67_3,msq67_3v,
     &  qqb_zaj_vdecay,qqb_zaj_gvec_vdecay)
      call dips(6,p,6,7,4,sub67_4,sub67_4v,msq67_4,msq67_4v,
     &  qqb_zaj_vdecay,qqb_zaj_gvec_vdecay)

      do j=-nf,nf
      do k=-nf,nf
      
      if ( (j /= 0) .and. (k == -j) ) then
C-----half=statistical factor

      msq(1,j,k)=half*cf*msq37_4(j,k)*sub37_4(qq)
      msq(2,j,k)=half*cf*msq37_6(j,k)*sub37_6(qq)

      msq(3,j,k)=half*cf*msq47_3(j,k)*sub47_3(qq)
      msq(4,j,k)=half*cf*msq47_6(j,k)*sub47_6(qq)

      ! splitting for g->gg
      msq(5,j,k)=half*ca*(
     &     +msq67_3(j,k)*sub67_3(gg)
     &     +msq67_3v(j,k)*sub67_3v
     &  )
      msq(6,j,k)=half*ca*(
     &     +msq67_4(j,k)*sub67_4(gg)
     &     +msq67_4v(j,k)*sub67_4v
     &  )

      ! splitting for g -> qqb, identical flavors
      msq(3,j,k)= msq(3,j,k) + (half)*tr*(
     &     +msq47_3_cross(j,k)*sub47_3_cross(gq)
     &     -msq47_3v_cross(j,k)*sub47_3v_cross
     &  )
      msq(4,j,k)= msq(4,j,k) + (half)*tr*(
     &     +msq47_6_cross(j,k)*sub47_6_cross(gq)
     &     -msq47_6v_cross(j,k)*sub47_6v_cross
     &  )

      ! splitting for g -> qqb
      msq(5,j,k)= msq(5,j,k) + (nf-1 + half)*tr*(
     &     +msq67_3(j,k)*sub67_3(gq)
     &     -msq67_3v(j,k)*sub67_3v
     &  )
      msq(6,j,k)= msq(6,j,k) + (nf-1 + half)*tr*(
     &     +msq67_4(j,k)*sub67_4(gq)
     &     -msq67_4v(j,k)*sub67_4v
     &  )


      endif

      enddo
      enddo


      end


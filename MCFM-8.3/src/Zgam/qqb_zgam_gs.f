      subroutine qqb_zgam_gs(p,msq)
        use VVconfig_m
      implicit none
      include 'types.f'
c---Matrix element SUBTRACTION squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  (Z decaying into (p3,p4))+a(p5)+g(p6)
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ptilde.f'
      include 'qqgg.f'
      include 'frag.f'
      include 'ewcharge.f'
      include 'phot_dip.f'
      real(dp), intent(in) :: p(mxpart,4)
      real(dp), intent(out) :: msq(maxd,-nf:nf,-nf:nf)

      integer :: j,k,nd
      real(dp) :: msq16_2(-nf:nf,-nf:nf),msq26_1(-nf:nf,-nf:nf),
     & sub16_2(4),sub26_1(4),dummyv(-nf:nf,-nf:nf),dsubv
      real(dp) :: msq36_4(-nf:nf,-nf:nf),msq46_3(-nf:nf,-nf:nf),
     & sub36_4(4),sub46_3(4)
      real(dp) :: sub56_1,sub56_2,msq56_1(-nf:nf,-nf:nf),
     & msq56_2(-nf:nf,-nf:nf)
      external qqb_zgam_new,donothing_gvec
      external qqb_z1jet

      phot_dip(:) = .false.
      msq = 0._dp

      if (decayChannel() == decayQuarks) then
        ndmax = 2

        call dips(1,p,3,6,4,sub36_4,dsubv,msq36_4,dummyv,
     &      qqb_zgam_new,donothing_gvec)
        call dips(2,p,4,6,3,sub46_3,dsubv,msq46_3,dummyv,
     &      qqb_zgam_new,donothing_gvec)

        do j=-nf,nf;  do k=-nf,nf
          if ((j > 0 .and. k == - j) .or. (j < 0 .and. k == - j)) then
            msq(1,j,k) = 2._dp*cf*sub36_4(qq)*msq36_4(j,k)
            msq(2,j,k) = 2._dp*cf*sub46_3(qq)*msq46_3(j,k)
          endif
        enddo; enddo

      else

        if (frag) then
           ndmax = 3
        else
           ndmax = 2
        endif

c---- calculate both initial-initial dipoles
c---- note that we do not require the gg dipoles, so the v-type
c---- entries are left as dummies
        call dips(1,p,1,6,2,sub16_2,dsubv,msq16_2,dummyv,
     &   qqb_zgam_new,donothing_gvec)
        call dips(2,p,2,6,1,sub26_1,dsubv,msq26_1,dummyv,
     &   qqb_zgam_new,donothing_gvec)

        if (frag) then
           call dipsfrag(3,p,5,6,2,sub56_2,msq56_2,qqb_z1jet)
           phot_dip(3)=.true.
        endif
    
        do j=-nf,nf
        do k=-nf,nf

        if  ((j == 0) .and. (k == 0)) then
          cycle
        elseif  ((j > 0) .and. (k == -j)
     &          .or.(j < 0) .and. (k == -j)) then
           msq(1,j,k)=2._dp*cf*sub16_2(qq)*msq16_2(j,k)
           msq(2,j,k)=2._dp*cf*sub26_1(qq)*msq26_1(j,k)
        elseif ((j .ne. 0) .and. (k == 0)) then
           msq(2,j,k)=2._dp*tr*sub26_1(qg)*msq26_1(j,-j)
           if(frag) then
              msq(3,j,k)=Q(j)**2*sub56_2*msq56_2(j,k)
           endif
        elseif ((j == 0) .and. (k .ne. 0)) then
           msq(1,j,k)=2._dp*tr*sub16_2(qg)*msq16_2(-k,k)
           if(frag) then
              msq(3,j,k)=Q(k)**2*sub56_2*msq56_2(j,k)
           endif

        endif

        enddo
        enddo

      endif
   
      end

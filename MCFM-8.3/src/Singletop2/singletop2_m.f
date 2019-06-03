      module singletop2_m
        implicit none

        public :: singletop2_tree, singletop2_virt
        public :: singletop2_gs, singletop2_z
        public :: singletop2_gs_light, singletop2_gs_heavy
        private

        contains

      subroutine singletop2_tree(p,msq)
        use types
        use singletop2_realamps_m
        use anomcoup_tbW
      implicit none

c     Matrix element for t-bbar production
C     (nwz=+1)
c      b(-p1)+u(-p2)-->n(p3)+e^+(p4)+b(p5)+d(p6)
C     or for
C     (nwz=-1)
c      b~(-p1)+d(-p2)-->e^-(p3)+n~(p4)+b~(p5)+u(p6)
C     averaged(summed) over initial(final) colours and spins
c--NB average over spins only -- colour factors cancel
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'ewcouple.f'
      include 'ckm.f'
      include 'nwz.f'
      include 'zprods_com.f'
      include 'masses.f'

      real(dp), intent(in) :: p(mxpart,4)
      real(dp), intent(out) :: msq(-nf:nf,-nf:nf)

      abstract interface
        function sttree(p, ju,jb,jn,je,jc,jd)
            use types
            implicit none
            include 'mxpart.f'
            real(dp) :: sttree
            real(dp), intent(in) :: p(mxpart,4)
            integer, intent(in) :: ju,jb,jn,je,jc,jd
        end function
      end interface

      procedure (sttree), pointer :: tree => null ()

      if (mode_anomcoup) then
          tree => singletop2_amp_protos
      else
          tree => singletop2_amp_tree
      endif

      call spinoru(6,p,za,zb)

      associate ( fac => aveqq*xn**2 )

      msq(:,:) = 0._dp

      if (nwz == +1) then
          msq(+2,5) = fac*tree(p,1,2,3,4,5,6)
          msq(+4,5) = msq(+2,5)

          msq(5,+2) = fac*tree(p,2,1,3,4,5,6)
          msq(5,+4) = msq(5,+2)

          msq(-1,5) = fac*tree(p,6,2,3,4,5,1)
          msq(-3,5) = msq(-1,5)

          msq(5,-1) = fac*tree(p,6,1,3,4,5,2)
          msq(5,-3) = msq(5,-1)
      else
          ! use global CP transformation for anti-top amplitudes
          ! and then cross 1 <-> 6 and reidentify 3 <-> 4
          msq(1,-5) = fac*tree(p,6,2,4,3,5,1)
          msq(3,-5) = msq(1,-5)

          msq(-5,+1) = fac*tree(p,6,1,4,3,5,2)
          msq(-5,+3) = msq(-5,+1)

          msq(-2,-5) = fac*tree(p,1,2,4,3,5,6)
          msq(-4,-5) = msq(-2,-5)

          msq(-5,-2) = fac*tree(p,2,1,4,3,5,6)
          msq(-5,-4) = msq(-5,-2)
      endif

      end associate

      end

      subroutine singletop2_gs_light(p,msq)
        use types
        implicit none
        include 'nf.f'
        include 'mxpart.f'
        include 'maxd.f'

        real(dp), intent(in) :: p(mxpart,4)
        real(dp), intent(out) :: msq(maxd, -nf:nf, -nf:nf)

        call singletop2_gs(p,msq,.true.,.false.)
      end subroutine

      subroutine singletop2_gs_heavy(p,msq)
        use types
        implicit none
        include 'nf.f'
        include 'mxpart.f'
        include 'maxd.f'

        real(dp), intent(in) :: p(mxpart,4)
        real(dp), intent(out) :: msq(maxd, -nf:nf, -nf:nf)

        call singletop2_gs(p,msq,.false.,.true.)
      end subroutine

      subroutine singletop2_gs(p,msq,light,heavy)
        use ieee_arithmetic
        use types
        use singletop2_scale_m
        implicit none

        include 'constants.f'
        include 'nf.f'
        include 'mxpart.f'
        include 'ptilde.f'
        include 'nwz.f'

        include 'stopbmass.f'
        include 'qqgg.f'
        include 'breit.f'

        include 'zprods_decl.f'

        real(dp), intent(in) :: p(mxpart,4)
        real(dp), intent(out) :: msq(maxd, -nf:nf, -nf:nf)
        logical, intent(in) :: light,heavy

        real(dp)::
     &  msq17_6(-nf:nf,-nf:nf),msq17_5(-nf:nf,-nf:nf),
     &  msq27_5(-nf:nf,-nf:nf),msq27_6(-nf:nf,-nf:nf),
     &  msq16_2(-nf:nf,-nf:nf),msq26_1(-nf:nf,-nf:nf),
     &  msq67_1(-nf:nf,-nf:nf),msq57_2(-nf:nf,-nf:nf),msq67_2(-nf:nf,-nf:nf),msq57_1(-nf:nf,-nf:nf),
     &  sub17_6(4),sub17_2_light(4),sub17_2_heavy(4),sub17_5(4),
     &  sub27_5(4),sub27_6(4),sub27_1_light(4),sub27_1_heavy(4),
     &  sub57_2(4),sub57_1(4),sub57_6(4),
     &  sub67_2(4),sub67_1(4),sub67_5(4),
     &  sub16_2(4),sub26_1(4),
     &  dummyv(-nf:nf,-nf:nf),dummy(-nf:nf,-nf:nf),dsubv,
     &  msq17_2_light(-nf:nf,-nf:nf), msq17_2_heavy(-nf:nf,-nf:nf),
     &  msq27_1_light(-nf:nf,-nf:nf), msq27_1_heavy(-nf:nf,-nf:nf)

        integer :: j,k,ib,nd,isub

        external donothing_gvec

c       call spinoru(7,p,za,zb)

        msq17_6 = 0._dp
        msq17_5 = 0._dp
        msq27_5 = 0._dp
        msq27_6 = 0._dp
        msq16_2 = 0._dp
        msq26_1 = 0._dp
        msq67_1 = 0._dp
        msq57_2 = 0._dp
        msq67_2 = 0._dp
        msq57_1 = 0._dp
        sub17_6 = 0._dp
        sub17_2_light = 0._dp
        sub17_2_heavy = 0._dp
        sub17_5 = 0._dp
        sub27_5 = 0._dp
        sub27_6 = 0._dp
        sub27_1_light = 0._dp
        sub27_1_heavy = 0._dp
        sub57_2 = 0._dp
        sub57_1 = 0._dp
        sub57_6 = 0._dp
        sub67_2 = 0._dp
        sub67_1 = 0._dp
        sub67_5 = 0._dp
        sub16_2 = 0._dp
        sub26_1 = 0._dp
        dummyv = 0._dp
        dummy = 0._dp
        dsubv = 0._dp
        msq17_2_light = 0._dp
        msq17_2_heavy = 0._dp
        msq27_1_light = 0._dp
        msq27_1_heavy = 0._dp

        ndmax = 14

        if (light) then
          corr_islight = .true.
          corr_beam1 = .true.
          call dips(1,p,  1,7,6,sub17_6,dsubv,msq17_6,dummyv,singletop2_tree,donothing_gvec)
          call dips(2,p,  6,7,1,sub67_1,dsubv,msq67_1,dummyv,singletop2_tree,donothing_gvec)
        endif

        if (light) then
          corr_islight = .true.
          corr_beam1 = .false.
          call dips(3,p,  2,7,6,sub27_6,dsubv,msq27_6,dummyv,singletop2_tree,donothing_gvec)
          call dips(4,p,  6,7,2,sub67_2,dsubv,msq67_2,dummyv,singletop2_tree,donothing_gvec)
        endif

        if (light) then
          corr_islight = .true.
          corr_beam1 = .true.
          call dips(5,p, 1,7,2,sub17_2_light,dsubv,msq17_2_light,dummyv,singletop2_tree,donothing_gvec)
          call dips(6,p, 1,6,2,sub16_2,dsubv,msq16_2,dummyv,singletop2_tree,donothing_gvec)
        endif

        if (light) then
          corr_islight = .true.
          corr_beam1 = .false.
          call dips(7,p, 2,7,1,sub27_1_light,dsubv,msq27_1_light,dummyv,singletop2_tree,donothing_gvec)
          call dips(8,p, 2,6,1,sub26_1,dsubv,msq26_1,dummyv,singletop2_tree,donothing_gvec)
        endif

        if (heavy) then
          corr_islight = .false.
          corr_beam1 = .false.
          call dips(9,p,  2,7,5,sub27_5,dsubv,msq27_5,dummyv,singletop2_tree,donothing_gvec)
          call dips(10,p,  5,7,2,sub57_2,dsubv,msq57_2,dummyv,singletop2_tree,donothing_gvec)
        endif

        if (heavy) then
          corr_islight = .false.
          corr_beam1 = .true.
          call dips(11,p,  1,7,5,sub17_5,dsubv,msq17_5,dummyv,singletop2_tree,donothing_gvec)
          call dips(12,p,  5,7,1,sub57_1,dsubv,msq57_1,dummyv,singletop2_tree,donothing_gvec)
        endif

        ! g splits into b pair
        if (heavy) then
          corr_islight = .false.
          corr_beam1 = .false.
          call dips(14,p, 2,7,1,sub27_1_heavy,dsubv,msq27_1_heavy,dummyv,singletop2_tree,donothing_gvec)

          corr_islight = .false.
          corr_beam1 = .true.
          call dips(13,p, 1,7,2,sub17_2_heavy,dsubv,msq17_2_heavy,dummyv,singletop2_tree,donothing_gvec)
        endif

        ib=5*nwz

        msq = 0._dp

      do j=-nf,nf
      do k=-nf,nf

      if ((j /= 0) .and. (k == ib)) then
        msq(1,j,k) = 2._dp*cf*sub17_6(qq)*msq17_6(j,k) ! light
        msq(2,j,k) = 2._dp*cf*sub67_1(qq)*msq67_1(j,k) ! light

        msq(9,j,k) = 2._dp*cf*sub27_5(qq)*msq27_5(j,k) ! heavy
        msq(10,j,k) = 2._dp*cf*sub57_2(qq)*msq57_2(j,k) ! heavy
      elseif ((j == ib) .and. (k /= 0)) then
        msq(3,j,k) = 2._dp*cf*sub27_6(qq)*msq27_6(j,k) ! light
        msq(4,j,k) = 2._dp*cf*sub67_2(qq)*msq67_2(j,k) ! light

        msq(11,j,k) = 2._dp*cf*sub17_5(qq)*msq17_5(j,k) ! heavy
        msq(12,j,k) = 2._dp*cf*sub57_1(qq)*msq57_1(j,k) ! heavy
      elseif ((j == 0) .and. (k == ib)) then
        msq(5,j,k) = 2._dp*tr*sub17_2_light(qg) * sum(msq17_2_light(1:5,k)) ! light
        msq(6,j,k) = 2._dp*tr*sub16_2(qg) * sum(msq16_2(-5:-1,k)) ! light
      elseif ((j == ib) .and. (k==0)) then
        msq(7,j,k) = 2._dp*tr*sub27_1_light(qg) * sum(msq27_1_light(j,1:5)) ! light
        msq(8,j,k) = 2._dp*tr*sub26_1(qg) * sum(msq26_1(j,-5:-1)) ! light
      elseif ((j == 0) .and. (k /= 0)) then
        ! g splits into b pair, but crossed
        msq(13,j,k) = 2._dp*tr*sub17_2_heavy(qg)*msq17_2_heavy(ib,k)
      elseif ((j /= 0) .and. (k == 0)) then
        ! g splits into b pair
        msq(14,j,k) = 2._dp*tr*sub27_1_heavy(qg)*msq27_1_heavy(j,ib)
      endif

      enddo
      enddo

      end subroutine

      subroutine singletop2_z(p,z)
        use types
        use singletop2_scale_m
        implicit none

        include 'constants.f'
        include 'nf.f'
        include 'mxpart.f'
        include 'PR_new.f'
        include 'PR_stop.f'
        include 'agq.f'
        include 'scale.f'
        include 'qcdcouple.f'
        include 'nwz.f'

        real(dp), intent(in) :: p(mxpart,4), z
        real(dp) :: dot, if_qq, fi_qq, ii_qg

        integer :: is
        real(dp) :: xl16,xl25,xl26,xl15,xl12

        xl16 = log(-2*dot(p,1,6)/renscale_beam1_islight_onlight**2)
        xl25 = log(-2*dot(p,2,5)/renscale_beam2_isheavy_onheavy**2)

        xl15 = log(-2*dot(p,1,5)/renscale_beam1_isheavy_onheavy**2)
        xl26 = log(-2*dot(p,2,6)/renscale_beam2_islight_onlight**2)

        Q1 = zip
        Q2 = zip
        B1 = zip
        B2 = zip

        do is=1,3
          B1(q,q,b,is) = as_light_beam1/2/pi * cf*(if_qq(z,xl16,is) + fi_qq(z,xl16,is))
          B2(b,b,q,is) = as_heavy_beam2/2/pi * cf*(if_qq(z,xl25,is) + fi_qq(z,xl25,is)) ! massive line

          B1(b,b,q,is) = as_heavy_beam1/2/pi * cf*(if_qq(z,xl15,is) + fi_qq(z,xl15,is)) ! massive line
          B2(q,q,b,is) = as_light_beam2/2/pi * cf*(if_qq(z,xl26,is) + fi_qq(z,xl26,is))

          Q1(q,g,q,is) = as_light_beam1/2/pi * tr * ii_qg(z, log(2*dot(p,1,2)/renscale_beam1_islight_onlight**2),is)
          Q2(q,g,q,is) = as_light_beam2/2/pi * tr * ii_qg(z, log(2*dot(p,1,2)/renscale_beam2_islight_onlight**2),is)

          B1(q,g,q,is) = as_heavy_beam1/2/pi * tr * ii_qg(z, log(2*dot(p,1,2)/renscale_beam1_isheavy_onheavy**2),is)
          B2(q,g,q,is) = as_heavy_beam2/2/pi * tr * ii_qg(z, log(2*dot(p,1,2)/renscale_beam2_isheavy_onheavy**2),is)
        enddo

      end subroutine

      subroutine singletop2_virt(p,msqv,light,heavy)
        use types
        use singletop2_virtamps_m
        use singletop2_scale_m
        implicit none
        include 'constants.f'
        include 'nf.f'
        include 'mxpart.f'
        include 'qcdcouple.f'
        include 'ewcouple.f'
        include 'ckm.f'
        include 'nwz.f'
        include 'zprods_com.f'
        include 'scheme.f'

        real(dp), intent(in) :: p(mxpart,4)
        real(dp), intent(out) :: msqv(-nf:nf,-nf:nf)
        logical, intent(in) :: light,heavy

        if (light .and. heavy) then
          write (*,*) "light and heavy must be called separately"
          call abort
        endif

        scheme='tH-V'

        call spinoru(6,p,za,zb)

        msqv = 0._dp

        associate (fac => cf*aveqq*gw**4*xn**2)

        if (nwz == +1) then
            msqv(+2,5) = singletop2_amp_virt(p,1,2,3,4,5,6, renscale_beam1_islight_onlight**2,
     &            renscale_beam2_isheavy_onheavy**2, light,heavy)

            msqv(5,+2) = singletop2_amp_virt(p,2,1,3,4,5,6, renscale_beam2_islight_onlight**2,
     &            renscale_beam1_isheavy_onheavy**2, light,heavy)

            msqv(-1,5) = singletop2_amp_virt(p,6,2,3,4,5,1, renscale_beam1_islight_onlight**2,
     &            renscale_beam2_isheavy_onheavy**2, light,heavy)

            msqv(5,-1) = singletop2_amp_virt(p,6,1,3,4,5,2, renscale_beam2_islight_onlight**2,
     &            renscale_beam1_isheavy_onheavy**2, light,heavy)
        else 
          ! use global CP transformation for anti-top amplitudes
          ! and then cross 1 <-> 6 and reidentify 3 <-> 4
            msqv(+1,-5) = singletop2_amp_virt(p,6,2,4,3,5,1, renscale_beam1_islight_onlight**2,
     &            renscale_beam2_isheavy_onheavy**2, light,heavy)

            msqv(-5,+1) = singletop2_amp_virt(p,6,1,4,3,5,2, renscale_beam2_islight_onlight**2,
     &            renscale_beam1_isheavy_onheavy**2, light,heavy)

            msqv(-2,-5) = singletop2_amp_virt(p,1,2,4,3,5,6, renscale_beam1_islight_onlight**2,
     &            renscale_beam2_isheavy_onheavy**2, light,heavy)

            msqv(-5,-2) = singletop2_amp_virt(p,2,1,4,3,5,6, renscale_beam2_islight_onlight**2,
     &            renscale_beam1_isheavy_onheavy**2, light,heavy)
        endif

        ! for single top
        if (nwz == +1) then
            if (light) then
                msqv(+2,5) = msqv(+2,5) * fac * as_light_beam1/2/pi
                msqv(5,+2) = msqv(5,+2) * fac * as_light_beam2/2/pi

                msqv(-1,5) = msqv(-1,5) * fac * as_light_beam1/2/pi
                msqv(5,-1) = msqv(5,-1) * fac * as_light_beam2/2/pi
            elseif (heavy) then
                msqv(+2,5) = msqv(+2,5) * fac * as_heavy_beam2/2/pi
                msqv(5,+2) = msqv(5,+2) * fac * as_heavy_beam1/2/pi

                msqv(-1,5) = msqv(-1,5) * fac * as_heavy_beam2/2/pi
                msqv(5,-1) = msqv(5,-1) * fac * as_heavy_beam1/2/pi
            endif

            msqv(+4,5) = msqv(+2,5)
            msqv(5,+4) = msqv(5,+2)
            msqv(-3,5) = msqv(-1,5)
            msqv(5,-3) = msqv(5,-1)

        else ! anti top
            if (light) then
                msqv(+1,-5) = msqv(+1,-5) * fac * as_light_beam1/2/pi
                msqv(-5,+1) = msqv(-5,+1) * fac * as_light_beam2/2/pi

                msqv(-2,-5) = msqv(-2,-5) * fac * as_light_beam1/2/pi
                msqv(-5,-2) = msqv(-5,-2) * fac * as_light_beam2/2/pi
            elseif (heavy) then
                msqv(+1,-5) = msqv(+1,-5) * fac * as_heavy_beam2/2/pi
                msqv(-5,+1) = msqv(-5,+1) * fac * as_heavy_beam1/2/pi

                msqv(-2,-5) = msqv(-2,-5) * fac * as_heavy_beam2/2/pi
                msqv(-5,-2) = msqv(-5,-2) * fac * as_heavy_beam1/2/pi
            endif

            msqv(+3,-5) = msqv(+1,-5)
            msqv(-5,+3) = msqv(-5,+1)
            msqv(-4,-5) = msqv(-2,-5)
            msqv(-5,-4) = msqv(-5,-2)
        endif

        end associate

      end subroutine

      end module

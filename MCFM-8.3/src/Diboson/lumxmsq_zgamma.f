      subroutine lumxmsq_zgamma(p,xx,z1,z2,QB,order,xmsq)
          implicit none
          include 'types.f'
          include 'nf.f'
          include 'constants.f'
          include 'mxpart.f'
          include 'facscale.f'
          include 'qcdcouple.f'
          include 'tiny.f'
          include 'taucut.f'

          real(dp), intent(in) :: p(mxpart,4), xx(2), z1, z2, QB(2)
          integer, intent(in) :: order
          real(dp), intent(out) :: xmsq

          real(dp) :: soft1(-1:1), soft2(-1:3)
          real(dp) :: beama0(-5:5), beamb0(-5:5),
     &                beama1(-5:5, -1:1), beamb1(-5:5, -1:1),
     &                beama2(-5:5, -1:3), beamb2(-5:5, -1:3)
          real(dp) :: hard(2), bit, assemble

          real(dp) :: msq(-nf:nf, -nf:nf, 0:2), ggcontrib
          real(dp) :: Q, dot, msqpow(-5:5,-5:5)

          integer :: ih1, ih2
          common/density/ih1,ih2

          integer j,k

          call softqqbis(order,soft1,soft2)

          if (order >= 0) then
              call fdist(ih1,xx(1),facscale,beama0)
              call fdist(ih2,xx(2),facscale,beamb0)
          endif
          if (order >= 1) then
              call xbeam1bis(ih1,z1,xx(1),QB(1),beama1)
              call xbeam1bis(ih1,z2,xx(2),QB(2),beamb1)
          endif
          if (order >= 2) then
              call xbeam2bis(ih1,z1,xx(1),QB(1),beama2)
              call xbeam2bis(ih2,z2,xx(2),QB(2),beamb2)
          endif

          xmsq = 0

          call zgam_mat(p,msq)

! compute power corrections if required
          if ((incpowcorr) .or. (onlypowcorr)) then
              Q=sqrt(two*dot(p,1,2))
              call tau0_powcorr_qa(order,xx(1),xx(2),Q,beama0,beamb0,msqpow)
          endif
      
          do j=-nf,nf
              k=-j
              if (j*k > 0 .or. j*k == 0) cycle
              if (msq(j,k,0) < tiny) cycle

              hard(1) = msq(j,k,1)/msq(j,k,0)/ason2pi
              hard(2) = msq(j,k,2)/msq(j,k,0)/ason2pi**2

              bit = assemble(order,beama0(j),beamb0(k),
     &                         beama1(j,:),beamb1(k,:),
     &                         beama2(j,:),beamb2(k,:),soft1,soft2,hard)
     
              if (incpowcorr) then
                bit=bit+msqpow(j,k)+msqpow(j,0)+msqpow(0,k)
              endif
              if (onlypowcorr) then
                bit=msqpow(j,k)+msqpow(j,0)+msqpow(0,k)
              endif

              bit = bit * msq(j,k,0)

              xmsq = xmsq + bit
        
          enddo 

          if ((order >= 2) .and. (onlypowcorr .eqv. .false.)) then
            call gg_zgam(p,ggcontrib)
            xmsq = xmsq + ggcontrib*beama0(0)*beamb0(0)
          endif

      end subroutine

      subroutine set_anomcoup(p)
        implicit none
        include 'types.f'
        include 'mxpart.f'
        include 'anomcoup.f'
        include 'masses.f'

        real(dp), intent(in) :: p(mxpart,4)
        real(dp) :: shat, dotvec, xfac

        ! by using s345 (the lepton pair + the photon) we don't have to
        ! recompute the form factor for different parts
        shat = dotvec(p(3,:)+p(4,:)+p(5,:),p(3,:)+p(4,:)+p(5,:))

        if (tevscale > 0._dp) then
          xfac = 1._dp/(1 + shat/(tevscale*1.e3_dp)**2)
        else
          xfac = 1._dp
        endif

        hitZ(1) = xfac**3*h1Z/zmass**2
        hitZ(2) = xfac**4*h2Z/zmass**4
        hitZ(3) = xfac**3*h3Z/zmass**2
        hitZ(4) = xfac**4*h4Z/zmass**4
        hitgam(1) = xfac**3*h1gam/zmass**2
        hitgam(2) = xfac**4*h2gam/zmass**4
        hitgam(3) = xfac**3*h3gam/zmass**2
        hitgam(4) = xfac**4*h4gam/zmass**4


      end subroutine

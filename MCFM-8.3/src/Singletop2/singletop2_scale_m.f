      module singletop2_scale_m
          use types
        implicit none
        include 'maxd.f'
        include 'nf.f'

        public :: singletop2_scale_init
        public :: singletop2_scale_setup
        public :: singletop2_set_dipscale
        public :: singletop2_scale_reset
        public :: singletop2_fillAP

        ! we need to distinguish the following 8 cases:
        ! beam number (1,2) is (heavy,light) one and corrections are on (heavy,light) line
        ! both for renormalization and factorization scale
        ! let's just be verbose and stupid!

        real(dp), save, public, protected ::
     &          renscale_beam1_isheavy_onheavy,
     &          renscale_beam1_isheavy_onlight,
     &          renscale_beam1_islight_onheavy,
     &          renscale_beam1_islight_onlight,
     &          renscale_beam2_isheavy_onheavy,
     &          renscale_beam2_isheavy_onlight,
     &          renscale_beam2_islight_onheavy,
     &          renscale_beam2_islight_onlight,
     &          facscale_beam1_isheavy_onheavy,
     &          facscale_beam1_isheavy_onlight,
     &          facscale_beam1_islight_onheavy,
     &          facscale_beam1_islight_onlight,
     &          facscale_beam2_isheavy_onheavy,
     &          facscale_beam2_isheavy_onlight,
     &          facscale_beam2_islight_onheavy,
     &          facscale_beam2_islight_onlight

!$omp threadprivate(renscale_beam1_isheavy_onheavy)
!$omp threadprivate(renscale_beam1_isheavy_onlight)
!$omp threadprivate(renscale_beam1_islight_onheavy)
!$omp threadprivate(renscale_beam1_islight_onlight)
!$omp threadprivate(renscale_beam2_isheavy_onheavy)
!$omp threadprivate(renscale_beam2_isheavy_onlight)
!$omp threadprivate(renscale_beam2_islight_onheavy)
!$omp threadprivate(renscale_beam2_islight_onlight)
!$omp threadprivate(facscale_beam1_isheavy_onheavy)
!$omp threadprivate(facscale_beam1_isheavy_onlight)
!$omp threadprivate(facscale_beam1_islight_onheavy)
!$omp threadprivate(facscale_beam1_islight_onlight)
!$omp threadprivate(facscale_beam2_isheavy_onheavy)
!$omp threadprivate(facscale_beam2_isheavy_onlight)
!$omp threadprivate(facscale_beam2_islight_onheavy)
!$omp threadprivate(facscale_beam2_islight_onlight)

        ! PDF array couting for real emission
        integer, parameter, public :: i_beam1_islight_onlight = 1
        integer, parameter, public :: i_beam2_isheavy_onlight = 2

        integer, parameter, public :: i_beam1_islight_onheavy = 3
        integer, parameter, public :: i_beam2_isheavy_onheavy = 4

        integer, parameter, public :: i_beam1_isheavy_onlight = 5
        integer, parameter, public :: i_beam2_islight_onlight = 6

        integer, parameter, public :: i_beam1_isheavy_onheavy = 7
        integer, parameter, public :: i_beam2_islight_onheavy = 8

        ! PDF array couting for virtual, also argument for fillAP function
        integer, parameter, public :: i_beam1_light = 1
        integer, parameter, public :: i_beam1_light_z = 2

        integer, parameter, public :: i_beam2_light = 3
        integer, parameter, public :: i_beam2_light_z = 4

        integer, parameter, public :: i_beam1_heavy = 5
        integer, parameter, public :: i_beam1_heavy_z = 6

        integer, parameter, public :: i_beam2_heavy = 7
        integer, parameter, public :: i_beam2_heavy_z = 8

        real(dp), save, public :: singletop2_dipscale(maxd,2)
!$omp threadprivate(singletop2_dipscale)

        real(dp), save, public :: singletop2_dipole_pdfs(maxd,2,-nf:nf)
!$omp threadprivate(singletop2_dipole_pdfs)

        ! this works for the real emission
        real(dp), save, public :: singletop2_pdfs(8,-nf:nf)
!$omp threadprivate(singletop2_pdfs)

        real(dp), save, public, protected ::
     &      as_light_beam1, as_light_beam2, as_heavy_beam1, as_heavy_beam2
!$omp threadprivate(as_light_beam1, as_light_beam2, as_heavy_beam1, as_heavy_beam2)

        logical, save, public :: corr_islight
        logical, save, public :: corr_beam1
!$omp threadprivate(corr_islight, corr_beam1)


        logical, public, save :: use_DDIS
        private

        contains

      subroutine singletop2_scale_init()
        implicit none
        include 'dynamicscale.f'

        if (dynstring == 'DDIS') then
          use_DDIS = .true.
        else
          use_DDIS = .false.
        endif
      end subroutine

      subroutine singletop2_scale_reset()
        implicit none
        singletop2_dipscale = 0._dp
        singletop2_dipole_pdfs = 0._dp
        singletop2_pdfs = 0._dp
      end subroutine

      subroutine singletop2_fillAP(z, i, AP)
        implicit none
        include 'mxpart.f'
        include 'epinv.f'
        include 'constants.f'
        include 'b0.f'
        include 'agq.f'

        real(dp), intent(in) :: z
        integer, intent(in) :: i
        real(dp), intent(out) :: AP(-1:1, -1:1, 3)

        real(dp) :: epcorr
        real(dp) :: omz, ason2pi
        real(dp) :: facscale,renscale

        if (i==i_beam1_light) then
            ason2pi = as_light_beam1/2/pi
            facscale = facscale_beam1_islight_onlight
            renscale = renscale_beam1_islight_onlight
        elseif (i==i_beam2_heavy) then
            ason2pi = as_heavy_beam2/2/pi
            facscale = facscale_beam2_isheavy_onheavy
            renscale = renscale_beam2_isheavy_onheavy
        elseif (i==i_beam2_light) then
            ason2pi = as_light_beam2/2/pi
            facscale = facscale_beam2_islight_onlight
            renscale = renscale_beam2_islight_onlight
        elseif (i==i_beam1_heavy) then
            ason2pi = as_heavy_beam1/2/pi
            facscale = facscale_beam1_isheavy_onheavy
            renscale = renscale_beam1_isheavy_onheavy
        else
          write (*,*) "undefined case in singletop2_fillAP!"
          call abort
        endif

        omz = 1._dp-z
        epcorr=epinv+2._dp*log(renscale/facscale)

        AP(:,:,:) = 0._dp

        AP(q,q,1)=+ason2pi*Cf*1.5_dp*epcorr
        AP(q,q,2)=+ason2pi*Cf*(-1._dp-z)*epcorr
        AP(q,q,3)=+ason2pi*Cf*2._dp/omz*epcorr
c       AP(a,a,1)=+ason2pi*Cf*1.5_dp*epcorr
c       AP(a,a,2)=+ason2pi*Cf*(-1._dp-z)*epcorr
c       AP(a,a,3)=+ason2pi*Cf*2._dp/omz*epcorr

        AP(q,g,1)=0._dp
        AP(q,g,2)=ason2pi*Tr*(z**2+omz**2)*epcorr
        AP(q,g,3)=0._dp
c       AP(a,g,1)=0._dp
c       AP(a,g,2)=ason2pi*Tr*(z**2+omz**2)*epcorr
c       AP(a,g,3)=0._dp

c       AP(g,q,1)=0._dp
c       AP(g,q,2)=ason2pi*Cf*(1._dp+omz**2)/z*epcorr
c       AP(g,q,3)=0._dp
c       AP(g,a,1)=0._dp
c       AP(g,a,2)=ason2pi*Cf*(1._dp+omz**2)/z*epcorr
c       AP(g,a,3)=0._dp

c       AP(g,g,1)=+ason2pi*b0*epcorr
c       AP(g,g,2)=+ason2pi*xn*2._dp*(1._dp/z+z*omz-2._dp)*epcorr
c       AP(g,g,3)=+ason2pi*xn*2._dp/omz*epcorr

      end subroutine

      subroutine singletop2_set_dipscale(nd, ptrans, gsq)
        use types
        implicit none

        include 'maxd.f'
        include 'mxpart.f'
        include 'couple.f'
        include 'nlooprun.f'
        include 'constants.f'

        integer, intent(in) :: nd
        real(dp), intent(in) :: ptrans(mxpart,4)
        real(dp), intent(out) :: gsq

        real(dp) :: alphas

c       write (*,*) "setting dipole scale, nd = ", nd

        call singletop2_scale_setup(ptrans)

        if (corr_islight .and. corr_beam1) then
            gsq = 4._dp*pi*alphas(renscale_beam1_islight_onlight,amz,nlooprun)
            singletop2_dipscale(nd, 1) = facscale_beam1_islight_onlight
            singletop2_dipscale(nd, 2) = facscale_beam2_isheavy_onlight
        elseif (corr_islight .and. (corr_beam1 .eqv. .false.)) then
            gsq = 4._dp*pi*alphas(renscale_beam2_islight_onlight,amz,nlooprun)
            singletop2_dipscale(nd, 1) = facscale_beam1_isheavy_onlight
            singletop2_dipscale(nd, 2) = facscale_beam2_islight_onlight
        elseif ((corr_islight .eqv. .false.) .and. corr_beam1) then
            gsq = 4._dp*pi*alphas(renscale_beam1_isheavy_onheavy,amz,nlooprun)
            singletop2_dipscale(nd, 1) = facscale_beam1_isheavy_onheavy
            singletop2_dipscale(nd, 2) = facscale_beam2_islight_onheavy
        elseif ((corr_islight .eqv. .false.) .and. (corr_beam1 .eqv. .false.)) then
            gsq = 4._dp*pi*alphas(renscale_beam2_isheavy_onheavy,amz,nlooprun)
            singletop2_dipscale(nd, 1) = facscale_beam1_islight_onheavy
            singletop2_dipscale(nd, 2) = facscale_beam2_isheavy_onheavy
        else
            call abort
        endif

      end subroutine

      subroutine singletop2_scale_setup(p)
         use types
       implicit none
        include 'nf.f'
        include 'mxpart.f'
        include 'constants.f'
        include 'masses.f'
        include 'couple.f' ! amz
        include 'nlooprun.f' ! nlooprun

        include 'initialscales.f'
        include 'facscale.f'
        include 'scale.f'
        include 'dynamicscale.f'

        real(dp), intent(in) :: p(mxpart,4)
        real(dp) :: dotvec, alphas

        real(dp) :: renscale
        real(dp) :: minscale

        if (use_DDIS) then

        ! read as:
        ! scale for beam1 when this beam is connected to the light line and corrections are on the light line
        renscale_beam1_islight_onlight = -dotvec(p(2,:)+p(3,:)+p(4,:)+p(5,:), p(2,:)+p(3,:)+p(4,:)+p(5,:))
        renscale_beam2_isheavy_onlight = renscale_beam1_islight_onlight + mt**2

        renscale_beam1_islight_onheavy = -dotvec(-p(1,:) - p(6,:), -p(1,:) - p(6,:))
        renscale_beam2_isheavy_onheavy = renscale_beam1_islight_onheavy + mt**2

        renscale_beam1_isheavy_onlight = -dotvec(p(1,:)+p(3,:)+p(4,:)+p(5,:), p(1,:)+p(3,:)+p(4,:)+p(5,:)) + mt**2
        renscale_beam2_islight_onlight = renscale_beam1_isheavy_onlight - mt**2

        renscale_beam1_isheavy_onheavy = -dotvec(-p(2,:) - p(6,:), -p(2,:) - p(6,:)) + mt**2
        renscale_beam2_islight_onheavy = renscale_beam1_isheavy_onheavy - mt**2

        minscale = 1._dp

        renscale_beam1_islight_onlight = max(sqrt(renscale_beam1_islight_onlight),minscale)
        renscale_beam2_isheavy_onlight = max(sqrt(renscale_beam2_isheavy_onlight),minscale)

        renscale_beam1_islight_onheavy = max(sqrt(renscale_beam1_islight_onheavy),minscale)
        renscale_beam2_isheavy_onheavy = max(sqrt(renscale_beam2_isheavy_onheavy),minscale)

        renscale_beam1_isheavy_onlight = max(sqrt(renscale_beam1_isheavy_onlight),minscale)
        renscale_beam2_islight_onlight = max(sqrt(renscale_beam2_islight_onlight),minscale)

        renscale_beam1_isheavy_onheavy = max(sqrt(renscale_beam1_isheavy_onheavy),minscale)
        renscale_beam2_islight_onheavy = max(sqrt(renscale_beam2_islight_onheavy),minscale)

        ! debug fixed scale
c         renscale_beam1_islight_onlight = 200d0
c         renscale_beam2_isheavy_onlight = 200d0

c         renscale_beam1_islight_onheavy = 200d0
c         renscale_beam2_isheavy_onheavy = 200d0

c         renscale_beam1_isheavy_onlight = 200d0
c         renscale_beam2_islight_onlight = 200d0

c         renscale_beam1_isheavy_onheavy = 200d0
c         renscale_beam2_islight_onheavy = 200d0

          ! fixed scale for light, other for heavy
c         renscale_beam1_islight_onlight = 100d0
c         renscale_beam1_islight_onheavy = 100d0
c         renscale_beam2_islight_onlight = 100d0
c         renscale_beam2_islight_onheavy = 100d0

c         renscale_beam2_isheavy_onlight = 350d0
c         renscale_beam2_isheavy_onheavy = 350d0
c         renscale_beam1_isheavy_onlight = 350d0
c         renscale_beam1_isheavy_onheavy = 350d0

        ! equal factorization scales
        facscale_beam1_islight_onlight = renscale_beam1_islight_onlight
        facscale_beam2_isheavy_onlight = renscale_beam2_isheavy_onlight
        facscale_beam1_islight_onheavy = renscale_beam1_islight_onheavy
        facscale_beam2_isheavy_onheavy = renscale_beam2_isheavy_onheavy
        facscale_beam1_isheavy_onlight = renscale_beam1_isheavy_onlight
        facscale_beam2_islight_onlight = renscale_beam2_islight_onlight
        facscale_beam1_isheavy_onheavy = renscale_beam1_isheavy_onheavy
        facscale_beam2_islight_onheavy = renscale_beam2_islight_onheavy

        ! as for light and heavy line corrections, calculated to formula depending on
        ! whether beam1 or beam2 holds the corrections

        as_light_beam1 = alphas(renscale_beam1_islight_onlight,amz,nlooprun)
        as_light_beam2 = alphas(renscale_beam2_islight_onlight,amz,nlooprun)

        as_heavy_beam1 = alphas(renscale_beam1_isheavy_onheavy,amz,nlooprun)
        as_heavy_beam2 = alphas(renscale_beam2_isheavy_onheavy,amz,nlooprun)

        else ! use other scale settings

        if (dynamicscale) then
          call scaleset(initscale,initfacscale,p)
        endif

        renscale = sqrt(musq)

        renscale_beam1_islight_onlight = renscale
        renscale_beam2_isheavy_onlight = renscale

        renscale_beam1_islight_onheavy = renscale
        renscale_beam2_isheavy_onheavy = renscale

        renscale_beam1_isheavy_onlight = renscale
        renscale_beam2_islight_onlight = renscale

        renscale_beam1_isheavy_onheavy = renscale
        renscale_beam2_islight_onheavy = renscale


        facscale_beam1_islight_onlight = facscale
        facscale_beam2_isheavy_onlight = facscale

        facscale_beam1_islight_onheavy = facscale
        facscale_beam2_isheavy_onheavy = facscale

        facscale_beam1_isheavy_onlight = facscale
        facscale_beam2_islight_onlight = facscale

        facscale_beam1_isheavy_onheavy = facscale
        facscale_beam2_islight_onheavy = facscale

        as_light_beam1 = alphas(renscale_beam1_islight_onlight,amz,nlooprun)
        as_light_beam2 = alphas(renscale_beam2_islight_onlight,amz,nlooprun)

        as_heavy_beam1 = alphas(renscale_beam1_isheavy_onheavy,amz,nlooprun)
        as_heavy_beam2 = alphas(renscale_beam2_isheavy_onheavy,amz,nlooprun)

        endif

      end subroutine

      end module


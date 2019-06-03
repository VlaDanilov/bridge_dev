      module VVconfig_m
          implicit none
      ! use these as bitflags, bosW should not be used directly,
      ! it's only for internal comparison purposes
      enum, bind(c)
          enumerator :: bosZ=1, bosW=2, bosWPlus=4, bosWMinus=8,
     &                  bosGamma=16
      endenum

      enum, bind(c)
          enumerator :: helLeft = 1, helRight = 2
          enumerator :: helMinus = 1, helPlus = 2
      endenum

c     enum, bind(c)
c         enumerator :: diagClassA = 1, diagClassB = 2,
c    &                  diagClassC = 3, diagClassF = 4
c     endenum

      enum, bind(c)
          enumerator :: decayUnset, decayElAntiEl, decayNuAntiNu, decayElAntiNu,
     &                  decayAntiElNu, decayQuarks, decayIgnore
      endenum

c     * schemeCatani returns finite amplitudes with Catani IR subtractions
c     * schemeMSBAR returns finite amplitudes with just IR pole
c           subtractions (finite Catani pieces added back)
c     * schemeMSFM returns fin. amplitudes with Catani IR-subtractions
c           (epinv) added back, and converted to DRED scheme
c           (see hep-ph/9610553 for DRED conversion factor, and a good
c            overview of different schemes)
      enum, bind(c)
          enumerator :: schemeCatani, schemeMSBAR, schemeMCFM
      endenum

      integer(kind(schemeMCFM)), save :: zgam_scheme = schemeMSBAR
!$omp threadprivate(zgam_scheme)

      integer (kind(decayElAntiEl)), private, save :: vDecay = decayUnset

      contains

      subroutine setDecayChannel(chan)
          implicit none
          integer (kind(decayElAntiEl)), intent(in) :: chan
          vDecay = chan
      end subroutine

      function decayChannel()
          implicit none
          include 'types.f'
          include 'mxpart.f'
          include 'is_functions_com.f'
          include 'plabel.f'

          integer (kind(decayElAntiEl)) :: decayChannel
          logical :: is_electron, is_neutrino, is_hadronic

          if (vDecay /= decayUnset) then
            decayChannel = vDecay
            return
          endif

          if (is_neutrino(3) .and. is_neutrino(4)) then
              decayChannel = decayNuAntiNu 
          elseif (is_electron(3) .and. is_electron(4)) then
              decayChannel =  decayElAntiEl
          elseif (is_hadronic(3) .and. is_hadronic(4)) then
              decayChannel = decayQuarks
          else
            call abort
            write (*,*) "undefined decay channel"
          endif

          vDecay = decayChannel

      end function

      end module

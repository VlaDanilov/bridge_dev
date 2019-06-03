      
      module helicities 
          implicit none
          enum, bind(c)
              enumerator :: helPlus = 1, helMinus = 2
              enumerator :: helRight = 1, helLeft = 2
          end enum
      end module

c      computes all helicity amplitudes. couplings stripped,
c      label ordering,  q qb g a lb l
c      from (g,a) = (+,+), (+,-), (-,+) and q=+, l=+ basis amps
c      this routine is only for q/isr-type amplitudes
        subroutine zaj_crossings(i1,i2,i3,i4,i5,i6,za,zb,
     &       amp_pp, amp_pm, amp_mp, amps)
        use helicities 
        implicit none
        include 'types.f'
        include 'constants.f'
        include 'mxpart.f'

        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        ! helicities for particles q, l, g, gam
        complex(dp), intent(out) :: amps(2,2,2,2)

        interface 
          function amp_pp(i1,i2,i3,i4,i5,i6,za,zb)
            integer, intent(in) :: i1,i2,i3,i4,i5,i6
            include 'types.f'
            include 'mxpart.f'
            complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
            complex(dp) :: amp_pp
          end function

          function amp_pm(i1,i2,i3,i4,i5,i6,za,zb)
            integer, intent(in) :: i1,i2,i3,i4,i5,i6
            include 'types.f'
            include 'mxpart.f'
            complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
            complex(dp) :: amp_pm
          end function

          function amp_mp(i1,i2,i3,i4,i5,i6,za,zb)
            integer, intent(in) :: i1,i2,i3,i4,i5,i6
            include 'types.f'
            include 'mxpart.f'
            complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
            complex(dp) :: amp_mp
          end function

        end interface  


        ! first block q=+, l=+
        amps(helPlus, helPlus, helPlus, helPlus) = 
     &     amp_pp(i1,i2,i3,i4,i5,i6,za,zb)
        amps(helPlus, helPlus, helPlus, helMinus) =
     &     amp_pm(i1,i2,i3,i4,i5,i6,za,zb)
        amps(helPlus, helPlus, helMinus, helPlus) =
     &     amp_mp(i1,i2,i3,i4,i5,i6,za,zb)
        amps(helPlus, helPlus, helMinus, helMinus) =
     &     amp_pp(i2,i1,i3,i4,i6,i5,zb,za)

        ! second block, q=+, l=-
        amps(helPlus, helMinus, helPlus, helPlus) = 
     &     amp_pp(i1,i2,i3,i4,i6,i5,za,zb)
        amps(helPlus, helMinus, helPlus, helMinus) =
     &     amp_pm(i1,i2,i3,i4,i6,i5,za,zb)
        amps(helPlus, helMinus, helMinus, helPlus) =
     &     amp_mp(i1,i2,i3,i4,i6,i5,za,zb)
        amps(helPlus, helMinus, helMinus, helMinus) = 
     &     amp_pp(i2,i1,i3,i4,i5,i6,zb,za)

        ! third block q=-, l=-
        amps(helMinus, helMinus, helPlus, helPlus) = 
     &     amp_pp(i2,i1,i3,i4,i6,i5,za,zb)
        amps(helMinus, helMinus, helPlus, helMinus) =
     &     amp_pm(i2,i1,i3,i4,i6,i5,za,zb)
        amps(helMinus, helMinus, helMinus, helPlus) =
     &     amp_mp(i2,i1,i3,i4,i6,i5,za,zb)
        amps(helMinus, helMinus, helMinus, helMinus) =
     &     amp_pp(i1,i2,i3,i4,i5,i6,zb,za)

        ! fourth block q=-, l=+
        amps(helMinus, helPlus, helPlus, helPlus) =
     &     amp_pp(i2,i1,i3,i4,i5,i6,za,zb)
        amps(helMinus, helPlus, helPlus, helMinus) = 
     &     amp_pm(i2,i1,i3,i4,i5,i6,za,zb)
        amps(helMinus, helPlus, helMinus, helPlus) =
     &     amp_mp(i2,i1,i3,i4,i5,i6,za,zb)
        amps(helMinus, helPlus, helMinus, helMinus) =
     &     amp_pp(i1,i2,i3,i4,i6,i5,zb,za)
     
      end subroutine


      ! this routine is identical to zaj_crossings, except for the added
      ! flip argument for the amplitudes. I prefer this for now instead of
      ! giving all other amps a dummy flip argument
      subroutine zaj_crossings_anom(i1,i2,i3,i4,i5,i6,za,zb,
     &       amp_pp, amp_pm, amp_mp, amps)
        use helicities 
        implicit none
        include 'types.f'
        include 'constants.f'
        include 'mxpart.f'

        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        ! helicities for particles q, l, g, gam
        complex(dp), intent(out) :: amps(2,2,2,2)

        interface 
          function amp_pp(i1,i2,i3,i4,i5,i6,za,zb,flip)
            integer, intent(in) :: i1,i2,i3,i4,i5,i6
            include 'types.f'
            include 'mxpart.f'
            complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
            logical, intent(in) :: flip
            complex(dp) :: amp_pp
          end function

          function amp_pm(i1,i2,i3,i4,i5,i6,za,zb,flip)
            integer, intent(in) :: i1,i2,i3,i4,i5,i6
            include 'types.f'
            include 'mxpart.f'
            complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
            logical, intent(in) :: flip
            complex(dp) :: amp_pm
          end function

          function amp_mp(i1,i2,i3,i4,i5,i6,za,zb,flip)
            integer, intent(in) :: i1,i2,i3,i4,i5,i6
            include 'types.f'
            include 'mxpart.f'
            complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
            logical, intent(in) :: flip
            complex(dp) :: amp_mp
          end function

        end interface  


        ! first block q=+, l=+
        amps(helPlus, helPlus, helPlus, helPlus) = 
     &     amp_pp(i1,i2,i3,i4,i5,i6,za,zb,.false.)
        amps(helPlus, helPlus, helPlus, helMinus) =
     &     amp_pm(i1,i2,i3,i4,i5,i6,za,zb,.false.)
        amps(helPlus, helPlus, helMinus, helPlus) =
     &     amp_mp(i1,i2,i3,i4,i5,i6,za,zb,.false.)
        amps(helPlus, helPlus, helMinus, helMinus) =
     &     amp_pp(i2,i1,i3,i4,i6,i5,zb,za,.true.)

        ! second block, q=+, l=-
        amps(helPlus, helMinus, helPlus, helPlus) = 
     &     amp_pp(i1,i2,i3,i4,i6,i5,za,zb,.false.)
        amps(helPlus, helMinus, helPlus, helMinus) =
     &     amp_pm(i1,i2,i3,i4,i6,i5,za,zb,.false.)
        amps(helPlus, helMinus, helMinus, helPlus) =
     &     amp_mp(i1,i2,i3,i4,i6,i5,za,zb,.false.)
        amps(helPlus, helMinus, helMinus, helMinus) = 
     &     amp_pp(i2,i1,i3,i4,i5,i6,zb,za,.true.)

        ! third block q=-, l=-
        amps(helMinus, helMinus, helPlus, helPlus) = 
     &     amp_pp(i2,i1,i3,i4,i6,i5,za,zb,.false.)
        amps(helMinus, helMinus, helPlus, helMinus) =
     &     amp_pm(i2,i1,i3,i4,i6,i5,za,zb,.false.)
        amps(helMinus, helMinus, helMinus, helPlus) =
     &     amp_mp(i2,i1,i3,i4,i6,i5,za,zb,.false.)
        amps(helMinus, helMinus, helMinus, helMinus) =
     &     amp_pp(i1,i2,i3,i4,i5,i6,zb,za,.true.)

        ! fourth block q=-, l=+
        amps(helMinus, helPlus, helPlus, helPlus) =
     &     amp_pp(i2,i1,i3,i4,i5,i6,za,zb, .false.)
        amps(helMinus, helPlus, helPlus, helMinus) = 
     &     amp_pm(i2,i1,i3,i4,i5,i6,za,zb, .false.)
        amps(helMinus, helPlus, helMinus, helPlus) =
     &     amp_mp(i2,i1,i3,i4,i5,i6,za,zb, .false.)
        amps(helMinus, helPlus, helMinus, helMinus) =
     &     amp_pp(i1,i2,i3,i4,i6,i5,zb,za, .true.)
     
      end subroutine

c      computes all helicity amplitudes. couplings stripped,
c      label ordering, q qb g a lb l
c      from (g,a) = (+,+), (+,-) and q=+, l=+ basis amps
c      this routine is only for l/fsr-type amplitudes
        subroutine zaj_crossings_l(i1,i2,i3,i4,i5,i6,za,zb,amp_pp,amp_pm,amps)
        use helicities 
        implicit none
        include 'types.f'
        include 'constants.f'
        include 'mxpart.f'

        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        ! helicities for particles q, l, g, gam
        complex(dp), intent(out) :: amps(2,2,2,2)

        interface 
          function amp_pp(i1,i2,i3,i4,i5,i6,za,zb)
            integer, intent(in) :: i1,i2,i3,i4,i5,i6
            include 'types.f'
            include 'mxpart.f'
            complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
            complex(dp) :: amp_pp
          end function

          function amp_pm(i1,i2,i3,i4,i5,i6,za,zb)
            integer, intent(in) :: i1,i2,i3,i4,i5,i6
            include 'types.f'
            include 'mxpart.f'
            complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
            complex(dp) :: amp_pm
          end function

        end interface  

        ! first block q=+, l=+
        amps(helPlus, helPlus, helPlus, helPlus) = 
     &     amp_pp(i1,i2,i3,i4,i5,i6,za,zb)
        amps(helPlus, helPlus, helPlus, helMinus) =
     &     amp_pm(i1,i2,i3,i4,i5,i6,za,zb)
        amps(helPlus, helPlus, helMinus, helPlus) =
     &     amp_pm(i2,i1,i3,i4,i6,i5,zb,za)
        amps(helPlus, helPlus, helMinus, helMinus) =
     &     amp_pp(i2,i1,i3,i4,i6,i5,zb,za)

        ! second block, q=+, l=-, in principle with a minus sign
        amps(helPlus, helMinus, helPlus, helPlus) = 
     &     amp_pp(i1,i2,i3,i4,i6,i5,za,zb)
        amps(helPlus, helMinus, helPlus, helMinus) =
     &     amp_pm(i1,i2,i3,i4,i6,i5,za,zb)
        amps(helPlus, helMinus, helMinus, helPlus) =
     &     amp_pm(i2,i1,i3,i4,i5,i6,zb,za)
        amps(helPlus, helMinus, helMinus, helMinus) = 
     &     amp_pp(i2,i1,i3,i4,i5,i6,zb,za)

        ! third block q=-, l=-
        amps(helMinus, helMinus, helPlus, helPlus) = 
     &     amp_pp(i2,i1,i3,i4,i6,i5,za,zb)
        amps(helMinus, helMinus, helPlus, helMinus) =
     &     amp_pm(i2,i1,i3,i4,i6,i5,za,zb)
        amps(helMinus, helMinus, helMinus, helPlus) =
     &     amp_pm(i1,i2,i3,i4,i5,i6,zb,za)
        amps(helMinus, helMinus, helMinus, helMinus) =
     &     amp_pp(i1,i2,i3,i4,i5,i6,zb,za)

        ! fourth block q=-, l=+, in principle with a minus sign
        amps(helMinus, helPlus, helPlus, helPlus) =
     &     amp_pp(i2,i1,i3,i4,i5,i6,za,zb)
        amps(helMinus, helPlus, helPlus, helMinus) = 
     &     amp_pm(i2,i1,i3,i4,i5,i6,za,zb)
        amps(helMinus, helPlus, helMinus, helPlus) =
     &     amp_pm(i1,i2,i3,i4,i6,i5,zb,za)
        amps(helMinus, helPlus, helMinus, helMinus) =
     &     amp_pp(i1,i2,i3,i4,i6,i5,zb,za)
     
      end subroutine


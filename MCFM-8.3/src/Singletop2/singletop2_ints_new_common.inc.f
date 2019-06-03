        use types

       public :: initamp

       public :: ampHeavyMM
       public :: ampHeavyMMtest
       public :: ampHeavyMMrealc8
       public :: ampHeavyMMrealc3
       public :: ampHeavyMMimagc3
       public :: ampHeavyMMrealc6
       public :: ampHeavyMMimagc6
       public :: ampHeavyMPPPrealc9

       public :: ampHeavyPMrealc2
       public :: ampHeavyPMimagc2
       public :: ampHeavyPMrealc4
       public :: ampHeavyPMimagc4
       public :: ampHeavyPMrealc7
       public :: ampHeavyPMimagc7

       public :: ampHeavyMPrealc2
       public :: ampHeavyMPimagc2
       public :: ampHeavyMPrealc4
       public :: ampHeavyMPimagc4
       public :: ampHeavyMPrealc7
       public :: ampHeavyMPimagc7

       logical, public, save :: initialized = .false.
!$omp threadprivate(initialized)

      private

      complex(dp), parameter :: im = (0._dp, 1._dp)

      complex(dp) :: M1000(0:1), M1110(0:1), M1111(0:2), M110S34(0:1)
      complex(dp) :: M0111(0:2), M1011(0:1), M1101(0:1), M1100(0:1)
      complex(dp) :: M11S126(0:1), M1111D2(0:0)
!$omp threadprivate(M1000, M1110, M1111, M110S34, M0111, M1011, M1101, M1100, M11S126, M1111D2)
      complex(dp) :: mtsq
!$omp threadprivate(mtsq)
      real(dp) :: s126, s12, s13, s14, s15, s16, s23, s24, s25, s26, s34, s35, s36, s45, s46, s56
!$omp threadprivate(s126, s12, s13, s14, s15, s16, s23, s24, s25, s26, s34, s35, s36, s45, s46, s56)
      real(dp) :: musq
!$omp threadprivate(musq)
      real(dp) :: epinv,epinv2
!$omp threadprivate(epinv,epinv2)


      complex(dp), public :: mytree, struc11step4, struc4step5, struc8step5, struc9step5
!$omp threadprivate(mytree, struc11step4, struc4step5, struc8step5, struc9step5)

      complex(dp) :: st10MPPP, st11MPPP, st12MPPP, st13MPPP, st14MPPP, st15MPPP, st16MPPP
      complex(dp) :: st17MPPP, st1MPPP, st27MPPP, st2MPPP, st37MPPP, st39MPPP, st3MPPP
      complex(dp) :: st55MPPP, st56MPPP, st57MPPP, st59MPPP, st5MPPP, st7MPPP, st8MPPP
      complex(dp) :: st9MPPP
!$omp threadprivate(st10MPPP, st11MPPP, st12MPPP, st13MPPP, st14MPPP, st15MPPP, st16MPPP)
!$omp threadprivate(st17MPPP, st1MPPP, st27MPPP, st2MPPP, st37MPPP, st39MPPP, st3MPPP)
!$omp threadprivate(st55MPPP, st56MPPP, st57MPPP, st59MPPP, st5MPPP, st7MPPP, st8MPPP)
!$omp threadprivate(st9MPPP)

      complex(dp), public :: struc10MP, struc13MP, struc14MP, struc15MP, struc17MP, struc18MP, struc19MP
      complex(dp), public :: struc1MP, struc22MP, struc24MP, struc27MP, struc29MP, struc2MP
      complex(dp), public :: struc31MP, struc33MP, struc34MP, struc35MP, struc36MP, struc37MP
      complex(dp), public :: struc38MP, struc39MP, struc3MP, struc42MP, struc43MP, struc44MP
      complex(dp), public :: struc45MP, struc47MP, struc50MP, struc51MP, struc53MP, struc54MP
      complex(dp), public :: struc55MP, struc56MP, struc57MP, struc58MP, struc59MP, struc5MP
      complex(dp), public :: struc6MP, struc7MP, struc8MP, struc9MP, struc41MP
!$omp threadprivate(struc10MP, struc13MP, struc14MP, struc15MP, struc17MP, struc18MP, struc19MP)
!$omp threadprivate(struc1MP, struc22MP, struc24MP, struc27MP, struc29MP, struc2MP)
!$omp threadprivate(struc31MP, struc33MP, struc34MP, struc35MP, struc36MP, struc37MP)
!$omp threadprivate(struc38MP, struc39MP, struc3MP, struc42MP, struc43MP, struc44MP)
!$omp threadprivate(struc45MP, struc47MP, struc50MP, struc51MP, struc53MP, struc54MP)
!$omp threadprivate(struc55MP, struc56MP, struc57MP, struc58MP, struc59MP, struc5MP)
!$omp threadprivate(struc6MP, struc7MP, struc8MP, struc9MP, struc41MP)

      complex(dp), public :: struc10PM, struc11PM, struc12PM, struc13PM, struc14PM, struc15PM, struc16PM
      complex(dp), public :: struc17PM, struc18PM, struc19PM, struc1PM, struc21PM, struc22PM
      complex(dp), public :: struc23PM, struc24PM, struc25PM, struc26PM, struc27PM, struc28PM
      complex(dp), public :: struc29PM, struc2PM, struc31PM, struc32PM, struc33PM, struc34PM
      complex(dp), public :: struc35PM, struc36PM, struc37PM, struc38PM, struc3PM, struc42PM
      complex(dp), public :: struc43PM, struc44PM, struc45PM, struc47PM, struc48PM, struc49PM
      complex(dp), public :: struc50PM, struc51PM, struc52PM, struc55PM, struc56PM, struc57PM
      complex(dp), public :: struc59PM, struc5PM, struc6PM, struc7PM, struc8PM, struc9PM
!$omp threadprivate(struc10PM, struc11PM, struc12PM, struc13PM, struc14PM, struc15PM, struc16PM)
!$omp threadprivate(struc17PM, struc18PM, struc19PM, struc1PM, struc21PM, struc22PM)
!$omp threadprivate(struc23PM, struc24PM, struc25PM, struc26PM, struc27PM, struc28PM)
!$omp threadprivate(struc29PM, struc2PM, struc31PM, struc32PM, struc33PM, struc34PM)
!$omp threadprivate(struc35PM, struc36PM, struc37PM, struc38PM, struc3PM, struc42PM)
!$omp threadprivate(struc43PM, struc44PM, struc45PM, struc47PM, struc48PM, struc49PM)
!$omp threadprivate(struc50PM, struc51PM, struc52PM, struc55PM, struc56PM, struc57PM)
!$omp threadprivate(struc59PM, struc5PM, struc6PM, struc7PM, struc8PM, struc9PM)


      contains

      subroutine initstrucs_mm(za,zb,ju,jb,jn,je,jc,jd)
          implicit none
          include 'mxpart.f'

          complex(dp), intent(in) :: za(mxpart,4), zb(mxpart,4)
          integer, intent(in) :: ju,jb,jn,je,jc,jd

          mytree = za(jc,jn)*zb(ju,jb)*(zb(je,jc)*za(jc,jd)+zb(je,jn)*za(jn,jd))
          struc9Step5 = 2._dp*za(jb,jd)*za(jn,jc)*zb(jb,ju)*zb(je,jb)
          struc4Step5 = -2._dp*za(jn,jc)*za(jc,jd)*zb(jb,ju)*zb(jc,je)
          struc8Step5 = 2._dp*za(jn,jc)*za(jn,jd)*zb(jn,ju)*zb(je,jb)
          struc11step4 = -(za(jc,jd)*za(jn,jc)*za(jn,jd)*zb(jc,je)*zb(jd,jb)*zb(jn,ju))

      end subroutine

      subroutine initstrucs_mppp(za,zb,ju,jb,jn,je,jc,jd)
          implicit none
          include 'mxpart.f'

          complex(dp), intent(in) :: za(mxpart,4), zb(mxpart,4)
          integer, intent(in) :: ju,jb,jn,je,jc,jd

       st10MPPP = za(jb,jc)*za(jn,jc)*za(ju,je)*zb(jc,je)*zb(jd,je)
       st11MPPP = 2*za(ju,jb)*za(ju,jc)*za(ju,jn)*zb(jd,ju)*zb(je,ju)
       st12MPPP = 2*za(jb,jn)*za(ju,jb)*za(ju,jc)*zb(jd,ju)*zb(je,jb)
       st13MPPP = 2*za(jb,jn)*za(ju,jb)*za(ju,jc)*zb(jd,jb)*zb(je,ju)
       st14MPPP = 2*za(jb,jn)*za(ju,jc)*za(ju,jn)*zb(jd,jn)*zb(je,ju)
       st15MPPP = 2*za(jb,jn)*za(ju,jc)*za(ju,je)*zb(jd,je)*zb(je,ju)
       st16MPPP = 2*za(jn,jc)*za(ju,jb)*za(ju,jc)*zb(jc,je)*zb(jd,ju)
       st17MPPP = 2*za(jb,jn)*za(ju,jc)*za(ju,jn)*zb(jd,je)*zb(jn,ju)
       st1MPPP = 2*za(jb,jc)*za(ju,jn)*zb(jd,je)
       st27MPPP = 2*za(jb,jd)*za(ju,jc)*za(ju,jn)*zb(jd,je)*zb(jd,ju)
       st2MPPP = za(jb,jc)*za(ju,jb)*za(ju,jn)*zb(jd,jb)*zb(je,ju)
       st37MPPP = 2*za(jn,jc)*za(ju,jb)*za(ju,jn)*zb(jd,jn)*zb(je,ju)
       st39MPPP = 2*za(jb,jn)*za(jn,jc)*za(ju,jb)*zb(jd,jb)*zb(je,jn)
       st3MPPP = za(jb,jc)*za(jb,jn)*za(ju,jb)*zb(jd,jb)*zb(je,jb)
       st55MPPP = 2*za(jb,jn)*za(jc,jd)*za(ju,jb)*zb(jd,jb)*zb(jd,je)
       st56MPPP = 2*za(jb,jn)*za(jc,jd)*za(ju,jn)*zb(jd,je)*zb(jd,jn)
       st57MPPP = 2*za(jb,jn)*za(jc,jd)*za(ju,je)*zb(jd,je)**2
       st59MPPP = 4*za(jn,jc)*za(ju,jb)*zb(jd,je)
       st5MPPP = za(jb,jc)*za(jb,jn)*za(ju,jn)*zb(jd,jn)*zb(je,jb)
       st7MPPP = za(jb,jc)*za(jb,jn)*za(ju,je)*zb(jd,je)*zb(je,jb)
       st8MPPP = za(jb,jc)*za(jn,jc)*za(ju,jb)*zb(jc,je)*zb(jd,jb)
       st9MPPP = za(jb,jc)*za(jn,jc)*za(ju,jn)*zb(jc,je)*zb(jd,jn)

      end subroutine

      subroutine initstrucs_pm(za,zb,ju,jb,jn,je,jc,jd)
          implicit none
          include 'mxpart.f'

          complex(dp), intent(in) :: za(mxpart,4), zb(mxpart,4)
          integer, intent(in) :: ju,jb,jn,je,jc,jd

       struc10PM = -(za(je,jd)*za(jn,jc)*zb(jc,jb)*zb(jc,je)*zb(je,ju))
       struc11PM = -2*za(ju,jd)*za(ju,jn)*zb(jb,ju)*zb(jc,ju)*zb(je,ju)
       struc12PM = -2*za(jb,jn)*za(ju,jd)*zb(jb,ju)*zb(jc,ju)*zb(je,jb)
       struc13PM = -2*za(jb,jd)*za(ju,jn)*zb(jb,ju)*zb(jc,ju)*zb(je,jb)
       struc14PM = -2*za(jn,jd)*za(ju,jn)*zb(jc,ju)*zb(je,jb)*zb(jn,ju)
       struc15PM = -2*za(je,jd)*za(ju,jn)*zb(jc,ju)*zb(je,jb)*zb(je,ju)
       struc16PM = -2*za(jn,jc)*za(ju,jd)*zb(jb,ju)*zb(jc,je)*zb(jc,ju)
       struc17PM = -2*za(jn,jd)*za(ju,jn)*zb(jc,ju)*zb(je,ju)*zb(jn,jb)
       struc18PM = -(za(jb,jd)*za(ju,jn)**2*zb(jb,ju)*zb(jc,ju)*zb(je,ju)*zb(jn,jb))
       struc19PM = -(za(jb,jd)*za(jb,jn)*za(ju,jn)*zb(jb,ju)*zb(jc,ju)*zb(je,jb)*zb(jn,jb))
       struc1PM = -2*za(jn,jd)*zb(jc,jb)*zb(je,ju)
       struc21PM = -(za(jb,jn)*za(jn,jd)*za(ju,jn)*zb(jc,ju)*zb(je,jb)*zb(jn,jb)*zb(jn,ju))
       struc22PM = -(za(je,jd)*za(ju,jn)**2*zb(jc,ju)*zb(je,ju)**2*zb(jn,jb))
       struc23PM = -(za(jb,jn)*za(je,jd)*za(ju,jn)*zb(jc,ju)*zb(je,jb)*zb(je,ju)*zb(jn,jb))
       struc24PM = -(za(jb,jd)*za(jn,jc)*za(ju,jn)*zb(jb,ju)*zb(jc,je)*zb(jc,ju)*zb(jn,jb))
       struc25PM = -(za(jn,jc)*za(jn,jd)*za(ju,jn)*zb(jc,je)*zb(jc,ju)*zb(jn,jb)*zb(jn,ju))
       struc26PM = -(za(je,jd)*za(jn,jc)*za(ju,jn)*zb(jc,je)*zb(jc,ju)*zb(je,ju)*zb(jn,jb))
       struc27PM = -2*za(jn,jd)*za(ju,jd)*zb(jc,ju)*zb(jd,jb)*zb(je,ju)
       struc28PM = -(za(jb,jd)*za(ju,jd)*za(ju,jn)*zb(jb,ju)*zb(jc,ju)*zb(jd,jb)*zb(je,ju))
       struc29PM = -(za(jb,jd)*za(jb,jn)*za(ju,jd)*zb(jb,ju)*zb(jc,ju)*zb(jd,jb)*zb(je,jb))
       struc2PM = -(za(jb,jd)*za(ju,jn)*zb(jb,ju)*zb(jc,jb)*zb(je,ju))
       struc31PM = -(za(jb,jn)*za(jn,jd)*za(ju,jd)*zb(jc,ju)*zb(jd,jb)*zb(je,jb)*zb(jn,ju))
       struc32PM = -(za(je,jd)*za(ju,jd)*za(ju,jn)*zb(jc,ju)*zb(jd,jb)*zb(je,ju)**2)
       struc33PM = -(za(jb,jn)*za(je,jd)*za(ju,jd)*zb(jc,ju)*zb(jd,jb)*zb(je,jb)*zb(je,ju))
       struc34PM = -(za(jb,jd)*za(jn,jc)*za(ju,jd)*zb(jb,ju)*zb(jc,je)*zb(jc,ju)*zb(jd,jb))
       struc35PM = -(za(jn,jc)*za(jn,jd)*za(ju,jd)*zb(jc,je)*zb(jc,ju)*zb(jd,jb)*zb(jn,ju))
       struc36PM = -(za(je,jd)*za(jn,jc)*za(ju,jd)*zb(jc,je)*zb(jc,ju)*zb(jd,jb)*zb(je,ju))
       struc37PM = -2*za(jn,jd)*za(ju,jn)*zb(jb,ju)*zb(jc,jn)*zb(je,ju)
       struc38PM = -2*za(jb,jn)*za(jn,jd)*zb(jb,ju)*zb(jc,jn)*zb(je,jb)
       struc3PM = -(za(jb,jd)*za(jb,jn)*zb(jb,ju)*zb(jc,jb)*zb(je,jb))
       struc42PM = -2*za(jn,jc)*za(jn,jd)*zb(jb,ju)*zb(jc,je)*zb(jc,jn)
       struc43PM = -2*za(jn,jd)**2*zb(jc,jn)*zb(jd,jb)*zb(je,ju)
       struc44PM = -(za(jb,jd)*za(jn,jd)*za(ju,jn)*zb(jb,ju)*zb(jc,jn)*zb(jd,jb)*zb(je,ju))
       struc45PM = -(za(jb,jd)*za(jb,jn)*za(jn,jd)*zb(jb,ju)*zb(jc,jn)*zb(jd,jb)*zb(je,jb))
       struc47PM = -(za(jb,jn)*za(jn,jd)**2*zb(jc,jn)*zb(jd,jb)*zb(je,jb)*zb(jn,ju))
       struc48PM = -(za(je,jd)*za(jn,jd)*za(ju,jn)*zb(jc,jn)*zb(jd,jb)*zb(je,ju)**2)
       struc49PM = -(za(jb,jn)*za(je,jd)*za(jn,jd)*zb(jc,jn)*zb(jd,jb)*zb(je,jb)*zb(je,ju))
       struc50PM = -(za(jb,jd)*za(jn,jc)*za(jn,jd)*zb(jb,ju)*zb(jc,je)*zb(jc,jn)*zb(jd,jb))
       struc51PM = -(za(jn,jc)*za(jn,jd)**2*zb(jc,je)*zb(jc,jn)*zb(jd,jb)*zb(jn,ju))
       struc52PM = -(za(je,jd)*za(jn,jc)*za(jn,jd)*zb(jc,je)*zb(jc,jn)*zb(jd,jb)*zb(je,ju))
       struc55PM = -2*za(jb,jd)*za(jn,jd)*zb(jb,ju)*zb(jd,jc)*zb(je,jb)
       struc56PM = -2*za(jn,jd)**2*zb(jd,jc)*zb(je,jb)*zb(jn,ju)
       struc57PM = -2*za(je,jd)*za(jn,jd)*zb(jd,jc)*zb(je,jb)*zb(je,ju)
       struc59PM = -4*za(jn,jd)*zb(jb,ju)*zb(jc,je)
       struc5PM = -(za(jb,jn)*za(jn,jd)*zb(jc,jb)*zb(je,jb)*zb(jn,ju))
       struc6PM = -(za(je,jd)*za(ju,jn)*zb(jc,jb)*zb(je,ju)**2)
       struc7PM = -(za(jb,jn)*za(je,jd)*zb(jc,jb)*zb(je,jb)*zb(je,ju))
       struc8PM = -(za(jb,jd)*za(jn,jc)*zb(jb,ju)*zb(jc,jb)*zb(jc,je))
       struc9PM = -(za(jn,jc)*za(jn,jd)*zb(jc,jb)*zb(jc,je)*zb(jn,ju))

      end subroutine

      subroutine initstrucs_mp(za,zb,ju,jb,jn,je,jc,jd)
          implicit none
          include 'mxpart.f'

          complex(dp), intent(in) :: za(mxpart,4), zb(mxpart,4)
          integer, intent(in) :: ju,jb,jn,je,jc,jd


           struc10MP = za(jb,jc)*za(je,jd)*za(jn,jc)*zb(jc,je)*zb(je,ju)
           struc13MP = 2*za(jb,jd)*za(jb,jn)*za(ju,jc)*zb(jb,ju)*zb(je,ju)
           struc14MP = 2*za(jb,jn)*za(jn,jd)*za(ju,jc)*zb(je,ju)*zb(jn,ju)
           struc15MP = 2*za(jb,jn)*za(je,jd)*za(ju,jc)*zb(je,ju)**2
           struc17MP = 2*za(jb,jn)*za(jn,jd)*za(ju,jc)*zb(je,ju)*zb(jn,ju)
           struc18MP = za(jb,jd)*za(jb,jn)*za(ju,jc)*za(ju,jn)*zb(jb,ju)*zb(je,ju)*zb(jn,ju)
           struc19MP = za(jb,jd)*za(jb,jn)**2*za(ju,jc)*zb(jb,ju)*zb(je,jb)*zb(jn,ju)
           struc1MP = 2*za(jb,jc)*za(jn,jd)*zb(je,ju)
           struc22MP = za(jb,jn)*za(je,jd)*za(ju,jc)*za(ju,jn)*zb(je,ju)**2*zb(jn,ju)
           struc24MP = za(jb,jd)*za(jb,jn)*za(jn,jc)*za(ju,jc)*zb(jb,ju)*zb(jc,je)*zb(jn,ju)
           struc27MP = 2*za(jb,jd)*za(jn,jd)*za(ju,jc)*zb(jd,ju)*zb(je,ju)
           struc29MP = za(jb,jd)**2*za(jb,jn)*za(ju,jc)*zb(jb,ju)*zb(jd,ju)*zb(je,jb)
           struc2MP = za(jb,jc)*za(jb,jd)*za(ju,jn)*zb(jb,ju)*zb(je,ju)
           struc31MP = za(jb,jd)*za(jb,jn)*za(jn,jd)*za(ju,jc)*zb(jd,ju)*zb(je,jb)*zb(jn,ju)
           struc33MP = za(jb,jd)*za(jb,jn)*za(je,jd)*za(ju,jc)*zb(jd,ju)*zb(je,jb)*zb(je,ju)
           struc34MP = za(jb,jd)**2*za(jn,jc)*za(ju,jc)*zb(jb,ju)*zb(jc,je)*zb(jd,ju)
           struc35MP = za(jb,jd)*za(jn,jc)*za(jn,jd)*za(ju,jc)*zb(jc,je)*zb(jd,ju)*zb(jn,ju)
           struc36MP = za(jb,jd)*za(je,jd)*za(jn,jc)*za(ju,jc)*zb(jc,je)*zb(jd,ju)*zb(je,ju)
           struc37MP = 2*za(jb,jd)*za(jn,jc)*za(ju,jn)*zb(je,ju)*zb(jn,ju)
           struc38MP = 2*za(jb,jd)*za(jb,jn)*za(jn,jc)*zb(je,jb)*zb(jn,ju)
           struc39MP = 2*za(jb,jd)*za(jb,jn)*za(jn,jc)*zb(jb,ju)*zb(je,jn)
           struc3MP = za(jb,jc)*za(jb,jd)*za(jb,jn)*zb(jb,ju)*zb(je,jb)
           struc42MP = 2*za(jb,jd)*za(jn,jc)**2*zb(jc,je)*zb(jn,ju)
           struc43MP = 2*za(jb,jd)*za(jn,jc)*za(jn,jd)*zb(jd,jn)*zb(je,ju)
           struc44MP = za(jb,jd)**2*za(jn,jc)*za(ju,jn)*zb(jb,ju)*zb(jd,jn)*zb(je,ju)
           struc45MP = za(jb,jd)**2*za(jb,jn)*za(jn,jc)*zb(jb,ju)*zb(jd,jn)*zb(je,jb)
           struc47MP = za(jb,jd)*za(jb,jn)*za(jn,jc)*za(jn,jd)*zb(jd,jn)*zb(je,jb)*zb(jn,ju)
           struc50MP = za(jb,jd)**2*za(jn,jc)**2*zb(jb,ju)*zb(jc,je)*zb(jd,jn)
           struc51MP = za(jb,jd)*za(jn,jc)**2*za(jn,jd)*zb(jc,je)*zb(jd,jn)*zb(jn,ju)
           struc53MP = -2*za(jb,jd)*za(jc,jd)*za(ju,jn)*zb(jd,ju)*zb(je,ju)
           struc54MP = -2*za(jb,jd)*za(jb,jn)*za(jc,jd)*zb(jd,ju)*zb(je,jb)
           struc55MP = 2*za(jb,jd)*za(jb,jn)*za(jc,jd)*zb(jb,ju)*zb(jd,je)
           struc56MP = 2*za(jb,jn)*za(jc,jd)*za(jn,jd)*zb(jd,je)*zb(jn,ju)
           struc57MP = 2*za(jb,jn)*za(jc,jd)*za(je,jd)*zb(jd,je)*zb(je,ju)
           struc58MP = -2*za(jb,jd)*za(jc,jd)*za(jn,jc)*zb(jc,je)*zb(jd,ju)
           struc59MP = 4*za(jb,jd)*za(jn,jc)*zb(je,ju)
           struc5MP = za(jb,jc)*za(jb,jn)*za(jn,jd)*zb(je,jb)*zb(jn,ju)
           struc6MP = za(jb,jc)*za(je,jd)*za(ju,jn)*zb(je,ju)**2
           struc7MP = za(jb,jc)*za(jb,jn)*za(je,jd)*zb(je,jb)*zb(je,ju)
           struc8MP = za(jb,jc)*za(jb,jd)*za(jn,jc)*zb(jb,ju)*zb(jc,je)
           struc9MP = za(jb,jc)*za(jn,jc)*za(jn,jd)*zb(jc,je)*zb(jn,ju)
           struc41MP = 2*za(jb,jn)*za(je,jd)*za(jn,jc)*zb(je,jn)*zb(je,ju)

      end subroutine

      subroutine initamp(musq_in, mtsq_in, za,zb,ju,jb,jn,je,jc,jd, epinv_in, epinv2_in)
          use mod_qcdloop_c
          implicit none
          include 'mxpart.f'

          real(dp), intent(in) :: musq_in
          complex(dp), intent(in) :: mtsq_in
          complex(dp), intent(in) :: za(mxpart,4), zb(mxpart,4)
          integer, intent(in) :: ju,jb,jn,je,jc,jd
          real(dp), intent(in) :: epinv_in, epinv2_in

          integer :: j,k,l

          real(dp) :: s,t
          complex(dp), parameter :: mzero = (0d0, 0d0)

          s(j,k) = real(za(j,k)*zb(k,j),dp)
          t(j,k,l) = s(j,k) + s(j,l) + s(k,l)

          epinv = epinv_in
          epinv2 = epinv2_in
          musq = musq_in
          mtsq = mtsq_in

          s12 = s(ju,jb)
          s13 = s(ju,jn)
          s14 = s(ju,je)
          s15 = s(ju,jc)
          s16 = s(ju,jd)
          s23 = s(jb,jn)
          s24 = s(jb,je)
          s25 = s(jb,jc)
          s26 = s(jb,jd)
          s34 = s(jn,je)
          s35 = s(jn,jc)
          s36 = s(jn,jd)
          s45 = s(je,jc)
          s46 = s(je,jd)
          s24 = s(jb,je)
          s56 = s(jc,jd)
          s126 = t(ju,jb,jd)

          call initstrucs_mm(za,zb,ju,jb,jn,je,jc,jd)
          call initstrucs_mp(za,zb,ju,jb,jn,je,jc,jd)
          call initstrucs_pm(za,zb,ju,jb,jn,je,jc,jd)
          call initstrucs_mppp(za,zb,ju,jb,jn,je,jc,jd)

          M1000(0) = qli1c(mtsq_in,musq_in,0)
          M1000(1) = qli1c(mtsq_in,musq_in,1)

          M11S126(0) = qli2c(dble(s126), mzero, mtsq_in, musq_in, 0)
          M11S126(1) = qli2c(dble(s126), mzero, mtsq_in, musq_in, 1)

          M1100(0) = qli2c(dble(s16), mzero, mtsq_in, musq_in, 0)
          M1100(1) = qli2c(dble(s16), mzero, mtsq_in, musq_in, 1)

          M110s34(0) = qli2c(dble(s34), mzero, mtsq_in, musq_in, 0)
          M110s34(1) = qli2c(dble(s34), mzero, mtsq_in, musq_in, 1)

          M1110(0) = qli3c(0d0, dble(s16), dble(s126), mzero, mzero, mtsq_in, musq_in, 0)
          M1110(1) = qli3c(0d0, dble(s16), dble(s126), mzero, mzero, mtsq_in, musq_in, 1)

          M1111(0) = qli4c(0d0, 0d0, dble(s34), dble(s16), dble(s25), dble(s126), mzero, mzero, mzero, mtsq_in, musq_in, 0)
          M1111(1) = qli4c(0d0, 0d0, dble(s34), dble(s16), dble(s25), dble(s126), mzero, mzero, mzero, mtsq_in, musq_in, 1)
          M1111(2) = qli4c(0d0, 0d0, dble(s34), dble(s16), dble(s25), dble(s126), mzero, mzero, mzero, mtsq_in, musq_in, 2)

          M0111(0) = qli3c(0d0, 0d0, dble(s25), mzero, mzero, mzero, musq_in, 0)
          M0111(1) = qli3c(0d0, 0d0, dble(s25), mzero, mzero, mzero, musq_in, 1)
          M0111(2) = qli3c(0d0, 0d0, dble(s25), mzero, mzero, mzero, musq_in, 2)

          M1101(0) = qli3c(dble(s25), dble(s34), dble(s16), mzero, mzero, mtsq_in, musq_in, 0)
          M1101(1) = qli3c(dble(s25), dble(s34), dble(s16), mzero, mzero, mtsq_in, musq_in, 1)

          M1011(0) = qli3c(0d0, dble(s126), dble(s34), mzero, mzero, mtsq_in, musq_in, 0)
          M1011(1) = qli3c(0d0, dble(s126), dble(s34), mzero, mzero, mtsq_in, musq_in, 1)

          M1111D2(0) = ints126s16s25s34mtx1111D2eps0()

          initialized = .true.

      end subroutine


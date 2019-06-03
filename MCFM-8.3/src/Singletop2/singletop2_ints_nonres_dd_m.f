      module singletop2_ints_nonres_dd_m
        use ddmodule
        use eftcouple
        use types

       public :: initamp_nonres

       public :: ampNonresonantHeavyMM
       public :: ampNonresonantLightFullMM

       public :: ampNonresonantHeavyPP
       public :: ampNonresonantLightFullPP

       public :: ampNonresonantHeavyReC4PM
       public :: ampNonresonantHeavyReC4MP
       public :: ampNonresonantHeavyImC4PM
       public :: ampNonresonantHeavyImC4MP

       public :: ampNonresonantHeavyReC7PM
       public :: ampNonresonantHeavyReC7MP
       public :: ampNonresonantHeavyImC7PM
       public :: ampNonresonantHeavyImC7MP

       public :: ampNonresonantLightFullReC4PM
       public :: ampNonresonantLightFullReC4MP
       public :: ampNonresonantLightFullImC4PM
       public :: ampNonresonantLightFullImC4MP

      private

      type(dd_complex), parameter :: im =  dd_complex([0.0d0, 0.0d0, 1.0d0, 0.0d0])

      ! W radiation on right side
      type(dd_complex) :: MB1001(0:1), MB1110(0:1), MB1101(0:0), MB1011(0:1)
      type(dd_complex) :: MB1111(0:2), MB1111D2(0:0)
!$omp threadprivate(MB1001, MB1110, MB1101, MB1011, MB1111,  MB1111D2)

      ! W radiation on left side
      type(dd_complex) :: MBL1001(0:1), MBL1110(0:1), MBL1101(0:0), MBL1011(0:1)
      type(dd_complex) :: MBL1111(0:2), MBL1111D2(0:0)
!$omp threadprivate(MBL1001, MBL1110, MBL1101, MBL1011, MBL1111, MBL1111D2)

      type(dd_complex) :: I300s25(0:2), I300s16(0:2)
!$omp threadprivate(I300s25, I300s16)

      type(dd_complex) :: propW16, propW34, propZ25
!$omp threadprivate (propW16, propW34, propZ25)


      type(dd_real) :: s126, s12, s13, s14, s15, s16, s23, s24, s25, s26, s34, s35, s36, s45, s46, s56
!$omp threadprivate(s126, s12, s13, s14, s15, s16, s23, s24, s25, s26, s34, s35, s36, s45, s46, s56)
      type(dd_real) :: musqLight, musqHeavy
!$omp threadprivate(musqLight, musqHeavy)
      type(dd_real) :: epinv,epinv2
!$omp threadprivate(epinv,epinv2)


      type(dd_complex), private :: mytree, struc11step4, struc4step5, struc8step5, struc9step5
!$omp threadprivate(mytree, struc11step4, struc4step5, struc8step5, struc9step5)

      type(dd_complex), private :: struc10MP, struc13MP, struc14MP, struc15MP, struc17MP, struc18MP, struc19MP
      type(dd_complex), private :: struc1MP, struc22MP, struc24MP, struc27MP, struc29MP, struc2MP
      type(dd_complex), private :: struc31MP, struc33MP, struc34MP, struc35MP, struc36MP, struc37MP
      type(dd_complex), private :: struc38MP, struc39MP, struc3MP, struc42MP, struc43MP, struc44MP
      type(dd_complex), private :: struc45MP, struc47MP, struc50MP, struc51MP, struc53MP, struc54MP
      type(dd_complex), private :: struc55MP, struc56MP, struc57MP, struc58MP, struc59MP, struc5MP
      type(dd_complex), private :: struc6MP, struc7MP, struc8MP, struc9MP, struc41MP
!$omp threadprivate(struc10MP, struc13MP, struc14MP, struc15MP, struc17MP, struc18MP, struc19MP)
!$omp threadprivate(struc1MP, struc22MP, struc24MP, struc27MP, struc29MP, struc2MP)
!$omp threadprivate(struc31MP, struc33MP, struc34MP, struc35MP, struc36MP, struc37MP)
!$omp threadprivate(struc38MP, struc39MP, struc3MP, struc42MP, struc43MP, struc44MP)
!$omp threadprivate(struc45MP, struc47MP, struc50MP, struc51MP, struc53MP, struc54MP)
!$omp threadprivate(struc55MP, struc56MP, struc57MP, struc58MP, struc59MP, struc5MP)
!$omp threadprivate(struc6MP, struc7MP, struc8MP, struc9MP, struc41MP)

      type(dd_complex) :: private, struc40MP, struc60MP, struc61MP
!$omp threadprivate(struc40MP, struc60MP, struc61MP)

      type(dd_complex), private :: struc10PM, struc11PM, struc12PM, struc13PM, struc14PM, struc15PM, struc16PM
      type(dd_complex), private :: struc17PM, struc18PM, struc19PM, struc1PM, struc21PM, struc22PM
      type(dd_complex), private :: struc23PM, struc24PM, struc25PM, struc26PM, struc27PM, struc28PM
      type(dd_complex), private :: struc29PM, struc2PM, struc31PM, struc32PM, struc33PM, struc34PM
      type(dd_complex), private :: struc35PM, struc36PM, struc37PM, struc38PM, struc3PM, struc42PM
      type(dd_complex), private :: struc43PM, struc44PM, struc45PM, struc47PM, struc48PM, struc49PM
      type(dd_complex), private :: struc50PM, struc51PM, struc52PM, struc55PM, struc56PM, struc57PM
      type(dd_complex), private :: struc59PM, struc5PM, struc6PM, struc7PM, struc8PM, struc9PM
!$omp threadprivate(struc10PM, struc11PM, struc12PM, struc13PM, struc14PM, struc15PM, struc16PM)
!$omp threadprivate(struc17PM, struc18PM, struc19PM, struc1PM, struc21PM, struc22PM)
!$omp threadprivate(struc23PM, struc24PM, struc25PM, struc26PM, struc27PM, struc28PM)
!$omp threadprivate(struc29PM, struc2PM, struc31PM, struc32PM, struc33PM, struc34PM)
!$omp threadprivate(struc35PM, struc36PM, struc37PM, struc38PM, struc3PM, struc42PM)
!$omp threadprivate(struc43PM, struc44PM, struc45PM, struc47PM, struc48PM, struc49PM)
!$omp threadprivate(struc50PM, struc51PM, struc52PM, struc55PM, struc56PM, struc57PM)
!$omp threadprivate(struc59PM, struc5PM, struc6PM, struc7PM, struc8PM, struc9PM)

      type(dd_complex), private :: struc60PM, struc61PM
!$omp threadprivate(struc60PM, struc61PM)

      type(dd_complex), private :: struc10PP, struc11PP, struc12PP, struc13PP, struc14PP, struc15PP, struc16PP
      type(dd_complex), private :: struc17PP, struc19PP, struc20PP, struc21PP, struc22PP, struc24PP
      type(dd_complex), private :: struc25PP, struc27PP, struc30PP, struc31PP, struc32PP, struc33PP
      type(dd_complex), private :: struc34PP, struc35PP, struc36PP, struc4PP, struc5PP, struc6PP, struc7PP
      type(dd_complex), private :: struc9PP

!$omp threadprivate(struc10PP, struc11PP, struc12PP, struc13PP, struc14PP, struc15PP, struc16PP)
!$omp threadprivate(struc17PP, struc19PP, struc20PP, struc21PP, struc22PP, struc24PP)
!$omp threadprivate(struc25PP, struc27PP, struc30PP, struc31PP, struc32PP, struc33PP)
!$omp threadprivate(struc34PP, struc35PP, struc36PP, struc4PP, struc5PP, struc6PP, struc7PP)
!$omp threadprivate(struc9PP)


      contains

      subroutine initstrucs_pp_nonres(za,zb,ju,jb,jn,je,jc,jd)
          implicit none
          include 'mxpart.f'

          complex(dp), intent(in) :: za(mxpart,4), zb(mxpart,4)
          integer, intent(in) :: ju,jb,jn,je,jc,jd

       struc10PP = 2*za(jb,jn)*za(jn,jd)*zb(jc,jn)*zb(je,ju)
       struc11PP = -(za(jb,jd)*za(je,jd)*za(jn,jc)*zb(jc,je)*zb(jd,jc)*zb(je,ju))
       struc12PP = -(za(jb,jd)*za(jn,jc)*za(jn,jd)*zb(jc,je)*zb(jd,jc)*zb(jn,ju))
       struc13PP = -(za(jb,jd)**2*za(jn,jc)*zb(jb,ju)*zb(jc,je)*zb(jd,jc))
       struc14PP = -(za(jb,jd)*za(jb,jn)*za(je,jd)*zb(jd,jc)*zb(je,jb)*zb(je,ju))
       struc15PP = -(za(jb,jd)*za(jb,jn)*za(jn,jd)*zb(jd,jc)*zb(je,jb)*zb(jn,ju))
       struc16PP = -(za(jb,jd)**2*za(jb,jn)*zb(jb,ju)*zb(jd,jc)*zb(je,jb))
       struc17PP = -(za(jb,jd)*za(je,jd)*za(ju,jn)*zb(jd,jc)*zb(je,ju)**2)
       struc19PP = -(za(jb,jd)**2*za(ju,jn)*zb(jb,ju)*zb(jd,jc)*zb(je,ju))
       struc20PP = -2*za(jb,jd)*za(jn,jd)*zb(jd,jc)*zb(je,ju)
       struc21PP = -(za(je,jd)*za(jn,jc)*za(ju,jb)*zb(jc,je)*zb(jc,ju)*zb(je,ju))
       struc22PP = -(za(jn,jc)*za(jn,jd)*za(ju,jb)*zb(jc,je)*zb(jc,ju)*zb(jn,ju))
       struc24PP = -(za(jb,jn)*za(je,jd)*za(ju,jb)*zb(jc,ju)*zb(je,jb)*zb(je,ju))
       struc25PP = -(za(jb,jn)*za(jn,jd)*za(ju,jb)*zb(jc,ju)*zb(je,jb)*zb(jn,ju))
       struc27PP = -(za(je,jd)*za(ju,jb)*za(ju,jn)*zb(jc,ju)*zb(je,ju)**2)
       struc30PP = -2*za(jn,jd)*za(ju,jb)*zb(jc,ju)*zb(je,ju)
       struc31PP = -2*za(jb,jd)*za(jn,jc)*zb(jc,je)*zb(jc,ju)
       struc32PP = -2*za(jb,jd)*za(jb,jn)*zb(jc,ju)*zb(je,jb)
       struc33PP = -2*za(jb,jd)*za(ju,jn)*zb(jc,ju)*zb(je,ju)
       struc34PP = 2*za(jb,jn)*za(je,jd)*zb(jc,je)*zb(je,ju)
       struc35PP = 2*za(jb,jn)*za(jn,jd)*zb(jc,je)*zb(jn,ju)
       struc36PP = 2*za(jb,jd)*za(jb,jn)*zb(jb,ju)*zb(jc,je)
       struc4PP = za(jb,jn)**2*za(je,jd)*zb(jc,jn)*zb(je,jb)*zb(je,ju)
       struc5PP = za(jb,jn)**2*za(jn,jd)*zb(jc,jn)*zb(je,jb)*zb(jn,ju)
       struc6PP = za(jb,jd)*za(jb,jn)**2*zb(jb,ju)*zb(jc,jn)*zb(je,jb)
       struc7PP = za(jb,jn)*za(je,jd)*za(ju,jn)*zb(jc,jn)*zb(je,ju)**2
       struc9PP = za(jb,jd)*za(jb,jn)*za(ju,jn)*zb(jb,ju)*zb(jc,jn)*zb(je,ju)

      end subroutine

      subroutine initstrucs_mm_nonres(za,zb,ju,jb,jn,je,jc,jd)
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

      subroutine initstrucs_pm_nonres(za,zb,ju,jb,jn,je,jc,jd)
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

       struc60PM = za(jb,jd)*za(jb,jn)*za(jn,je)*zb(jc,jb)*zb(je,jb)*zb(je,ju)*zb(jn,jb)
       struc61PM = za(jb,jd)*za(jn,jc)*za(jn,je)*zb(jc,jb)*zb(jc,je)*zb(je,ju)*zb(jn,jb)

      end subroutine

      subroutine initstrucs_mp_nonres(za,zb,ju,jb,jn,je,jc,jd)
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

           struc40MP = 2*za(jb,jn)*za(jn,jc)*za(jn,jd)*zb(je,jn)*zb(jn,ju)
           struc60MP = -(za(jb,jc)*za(jb,jd)*za(jb,jn)*za(jn,je)*zb(je,jb)*zb(je,ju)*zb(jn,jb))
           struc61MP = -(za(jb,jc)*za(jb,jd)*za(jn,jd)*za(jn,je)*zb(jd,je)*zb(je,ju)*zb(jn,jb))


      end subroutine

      subroutine initamp_nonres(musqHeavy_in, musqLight_in, za,zb,ju,jb,jn,je,jc,jd, epinv_in, epinv2_in, sc_in)
          use mod_qcdloop_c
          implicit none
          include 'mxpart.f'

          real(dp), intent(in) :: musqHeavy_in, musqLight_in
          complex(dp), intent(in) :: za(mxpart,4), zb(mxpart,4)
          integer, intent(in) :: ju,jb,jn,je,jc,jd
          real(dp), intent(in) :: epinv_in, epinv2_in
          real(dp), intent(in) :: sc_in

          integer :: j,k,l

          real(dp) :: s,t
          complex(dp), parameter :: mzero = (0d0, 0d0)

          real(dp) :: s134, s125

          s(j,k) = real(za(j,k)*zb(k,j),dp)
          t(j,k,l) = s(j,k) + s(j,l) + s(k,l)

          epinv = epinv_in
          epinv2 = epinv2_in
          musqHeavy = musqHeavy_in
          musqLight = musqLight_in

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

          call initstrucs_mm_nonres(za,zb,ju,jb,jn,je,jc,jd)
          call initstrucs_mp_nonres(za,zb,ju,jb,jn,je,jc,jd)
          call initstrucs_pm_nonres(za,zb,ju,jb,jn,je,jc,jd)
          call initstrucs_pp_nonres(za,zb,ju,jb,jn,je,jc,jd)

          propW16 = ddreal(1) / (s16 - wmass**2*sc_in)
          propW34 = ddreal(1) / (s34 - wmass**2*sc_in + im*wmass*wwidth*sc_in)
          propZ25 = ddreal(1) / (s25 - zmass**2*sc_in)

          s125 = s12 + s15 + s25
          s134 = s13 + s14 + s34

          ! W radiation on right side
          MB1001(0) = qli2(dble(s34), 0d0, 0d0, musqLight_in, 0)
          MB1001(1) = qli2(dble(s34), 0d0, 0d0, musqLight_in, 1)

          MB1110(0) = qli3(0d0,dble(s25),dble(s125), 0d0, 0d0, 0d0, musqLight_in, 0)
          MB1110(1) = qli3(0d0,dble(s25),dble(s125), 0d0, 0d0, 0d0, musqLight_in, 1)

          MB1101(0) = qli3(dble(s16), dble(s25), dble(s34), 0d0, 0d0, 0d0, musqLight_in, 0)

          MB1011(0) = qli3(0d0, dble(s125), dble(s34), 0d0, 0d0, 0d0, musqLight_in, 0)
          MB1011(1) = qli3(0d0, dble(s125), dble(s34), 0d0, 0d0, 0d0, musqLight_in, 1)

          MB1111(0) = qli4(0d0, 0d0, dble(s25), dble(s34),  dble(s16), dble(s125), 0d0, 0d0, 0d0, 0d0, musqLight_in, 0)
          MB1111(1) = qli4(0d0, 0d0, dble(s25), dble(s34),  dble(s16), dble(s125), 0d0, 0d0, 0d0, 0d0, musqLight_in, 1)
          MB1111(2) = qli4(0d0, 0d0, dble(s25), dble(s34),  dble(s16), dble(s125), 0d0, 0d0, 0d0, 0d0, musqLight_in, 2)

          ! light line corrections, W radiation on left side
          MBL1001(0) = qli2(dble(s25), 0d0, 0d0, musqLight_in, 0)
          MBL1001(1) = qli2(dble(s25), 0d0, 0d0, musqLight_in, 1)

          MBL1110(0) = qli3(0d0,dble(s34),dble(s134), 0d0, 0d0, 0d0, musqLight_in, 0)
          MBL1110(1) = qli3(0d0,dble(s34),dble(s134), 0d0, 0d0, 0d0, musqLight_in, 1)

          MBL1101(0) = qli3(dble(s16), dble(s34), dble(s25), 0d0, 0d0, 0d0, musqLight_in, 0)

          MBL1011(0) = qli3(0d0, dble(s134), dble(s25), 0d0, 0d0, 0d0, musqLight_in, 0)
          MBL1011(1) = qli3(0d0, dble(s134), dble(s25), 0d0, 0d0, 0d0, musqLight_in, 1)

          MBL1111(0) = qli4(0d0, 0d0, dble(s34), dble(s25),  dble(s16), dble(s134), 0d0, 0d0, 0d0, 0d0, musqLight_in, 0)
          MBL1111(1) = qli4(0d0, 0d0, dble(s34), dble(s25),  dble(s16), dble(s134), 0d0, 0d0, 0d0, 0d0, musqLight_in, 1)
          MBL1111(2) = qli4(0d0, 0d0, dble(s34), dble(s25),  dble(s16), dble(s134), 0d0, 0d0, 0d0, 0d0, musqLight_in, 2)

          ! heavy line vertex corrections
          I300s25(0) = qli3(0d0, 0d0, dble(s25), 0d0, 0d0, 0d0, musqHeavy_in, 0)
          I300s25(1) = qli3(0d0, 0d0, dble(s25), 0d0, 0d0, 0d0, musqHeavy_in, 1)
          I300s25(2) = qli3(0d0, 0d0, dble(s25), 0d0, 0d0, 0d0, musqHeavy_in, 2)

          ! light line vertex corrections
          I300s16(0) = qli3(0d0, 0d0, dble(s16), 0d0, 0d0, 0d0, musqLight_in, 0)
          I300s16(1) = qli3(0d0, 0d0, dble(s16), 0d0, 0d0, 0d0, musqLight_in, 1)
          I300s16(2) = qli3(0d0, 0d0, dble(s16), 0d0, 0d0, 0d0, musqLight_in, 2)

          MB1111D2(0) = intHs16s25s26s34s56x1111D2eps0()
          MBL1111D2(0) = intHLs16s25s26s34s56x1111D2eps0()

      end subroutine

      include 'newints_nonres_dd.inc.f'

      end module

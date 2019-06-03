       function intHL0s25s260s56x1010D0eps0()
           implicit none
           complex(dp) :: intHL0s25s260s56x1010D0eps0

           complex(dp) :: result

      result = MBL1001(0) + (s26 + s56) * MBL1011(1)

           intHL0s25s260s56x1010D0eps0 = result
       end function intHL0s25s260s56x1010D0eps0

       function intHL0s25s260s56x1010D0eps1()
           implicit none
           complex(dp) :: intHL0s25s260s56x1010D0eps1

           complex(dp) :: result

      result = 1.0_dp

           intHL0s25s260s56x1010D0eps1 = result
       end function intHL0s25s260s56x1010D0eps1

       function intHL0s25s260s56x1011D2eps0()
           implicit none
           complex(dp) :: intHL0s25s260s56x1011D2eps0
           complex(dp) ::  t1

           complex(dp) :: result

      t1 = 0.1e1_dp / 0.2e1_dp
      result = t1 * ((s25 + s26 + s56) * MBL1011(1) + MBL1001(0) + 1.0_d
     #p)

           intHL0s25s260s56x1011D2eps0 = result
       end function intHL0s25s260s56x1011D2eps0

       function intHL0s25s260s56x1011D2eps1()
           implicit none
           complex(dp) :: intHL0s25s260s56x1011D2eps1

           complex(dp) :: result

      result = 0.1e1_dp / 0.2e1_dp

           intHL0s25s260s56x1011D2eps1 = result
       end function intHL0s25s260s56x1011D2eps1

       function intHL0s25s260s56x1012D2eps0()
           implicit none
           complex(dp) :: intHL0s25s260s56x1012D2eps0
           complex(dp) ::  t1,t2

           complex(dp) :: result

      t1 = s26 + s56
      t2 = s25 + s26 + s56
      t1 = 0.1e1_dp / t1
      result = 2.0_dp * t2 * MBL1011(1) * t1 + t1 * (t2 * MBL1011(0) + M
     #BL1001(0))

           intHL0s25s260s56x1012D2eps0 = result
       end function intHL0s25s260s56x1012D2eps0

       function intHL0s25s260s56x1012D2eps1()
           implicit none
           complex(dp) :: intHL0s25s260s56x1012D2eps1
           complex(dp) ::  t1

           complex(dp) :: result

      t1 = s26 + s56
      t1 = 0.1e1_dp / t1
      result = t1 * ((s25 + s26 + s56) * MBL1011(1) + 1.0_dp)

           intHL0s25s260s56x1012D2eps1 = result
       end function intHL0s25s260s56x1012D2eps1

       function intHL0s25s260s56x1013D4eps0()
           implicit none
           complex(dp) :: intHL0s25s260s56x1013D4eps0
           complex(dp) ::  t1,t2,t3

           complex(dp) :: result

      t1 = s25 + s26 + s56
      t2 = s26 + s56
      t3 = 0.1e1_dp / t2 ** 2.0_dp
      result = t3 * ((2.0_dp * s25 + 3.0_dp * t2) * MBL1001(0) / 4.0_dp 
     #+ 0.3e1_dp / 0.2e1_dp * t1 ** 2.0_dp * MBL1011(1)) + t1 * t3 * (t1
     # * MBL1011(0) + 1.0_dp) / 2.0_dp

           intHL0s25s260s56x1013D4eps0 = result
       end function intHL0s25s260s56x1013D4eps0

       function intHL0s25s260s56x1013D4eps1()
           implicit none
           complex(dp) :: intHL0s25s260s56x1013D4eps1
           complex(dp) ::  t1,t2,t3

           complex(dp) :: result

      t1 = s26 + s56
      t2 = s25 + s26 + s56
      t3 = 0.1e1_dp / t1 ** 2.0_dp
      result = t3 * (s25 / 2.0_dp + 0.3e1_dp / 0.4e1_dp * t1 + t2 ** 2.0
     #_dp * MBL1011(1) / 2.0_dp)

           intHL0s25s260s56x1013D4eps1 = result
       end function intHL0s25s260s56x1013D4eps1

       function intHL0s25s260s56x1020D2eps0()
           implicit none
           complex(dp) :: intHL0s25s260s56x1020D2eps0
           complex(dp) ::  t1

           complex(dp) :: result

      t1 = 0.1e1_dp / 0.2e1_dp
      result = t1 * ((s26 + s56) * MBL1011(1) + MBL1001(0))

           intHL0s25s260s56x1020D2eps0 = result
       end function intHL0s25s260s56x1020D2eps0

       function intHL0s25s260s56x1020D2eps1()
           implicit none
           complex(dp) :: intHL0s25s260s56x1020D2eps1

           complex(dp) :: result

      result = 0.1e1_dp / 0.2e1_dp

           intHL0s25s260s56x1020D2eps1 = result
       end function intHL0s25s260s56x1020D2eps1

       function intHL0s25s260s56x1021D2eps0()
           implicit none
           complex(dp) :: intHL0s25s260s56x1021D2eps0
           complex(dp) ::  t1,t2

           complex(dp) :: result

      t1 = s26 + s56
      t2 = 0.1e1_dp / t1
      result = -(s25 * MBL1011(0) + (2.0_dp * s25 + t1) * MBL1011(1) + M
     #BL1001(0)) * t2

           intHL0s25s260s56x1021D2eps0 = result
       end function intHL0s25s260s56x1021D2eps0

       function intHL0s25s260s56x1021D2eps1()
           implicit none
           complex(dp) :: intHL0s25s260s56x1021D2eps1
           complex(dp) ::  t1

           complex(dp) :: result

      t1 = s26 + s56
      t1 = 0.1e1_dp / t1
      result = -t1 * (s25 * MBL1011(1) + 1.0_dp)

           intHL0s25s260s56x1021D2eps1 = result
       end function intHL0s25s260s56x1021D2eps1

       function intHL0s25s260s56x1022D4eps0()
           implicit none
           complex(dp) :: intHL0s25s260s56x1022D4eps0
           complex(dp) ::  t1,t2,t3

           complex(dp) :: result

      t1 = s26 + s56
      t2 = s25 + s26 + s56
      t3 = 0.1e1_dp / t1 ** 2.0_dp
      result = -t3 * ((MBL1001(0) + 1.0_dp) * (2.0_dp * s25 + t1) + t2 *
     # (6.0_dp * s25 + t1) * MBL1011(1)) / 2.0_dp - s25 * t2 * MBL1011(0
     #) * t3

           intHL0s25s260s56x1022D4eps0 = result
       end function intHL0s25s260s56x1022D4eps0

       function intHL0s25s260s56x1022D4eps1()
           implicit none
           complex(dp) :: intHL0s25s260s56x1022D4eps1
           complex(dp) ::  t1,t2

           complex(dp) :: result

      t1 = s26 + s56
      t2 = 0.1e1_dp / t1 ** 2.0_dp
      result = t2 * (-s25 - t1 / 2.0_dp - s25 * (s25 + s26 + s56) * MBL1
     #011(1))

           intHL0s25s260s56x1022D4eps1 = result
       end function intHL0s25s260s56x1022D4eps1

       function intHL0s25s260s56x1031D4eps0()
           implicit none
           complex(dp) :: intHL0s25s260s56x1031D4eps0
           complex(dp) ::  t1,t2,t3

           complex(dp) :: result

      t1 = s26 + s56
      t2 = 2.0_dp
      t3 = 0.1e1_dp / t1 ** 2.0_dp
      result = s25 * t3 * (s25 * MBL1011(0) + 1.0_dp) / 2.0_dp - t3 * (-
     #(s25 * t2 - t1) * MBL1001(0) + (-t2 * ((s25 - s26) * s56 + s25 * s
     #26) - 6.0_dp * s25 ** 2.0_dp + s26 ** 2.0_dp + s56 ** 2.0_dp) * MB
     #L1011(1)) / 4.0_dp

           intHL0s25s260s56x1031D4eps0 = result
       end function intHL0s25s260s56x1031D4eps0

       function intHL0s25s260s56x1031D4eps1()
           implicit none
           complex(dp) :: intHL0s25s260s56x1031D4eps1
           complex(dp) ::  t1,t2,t3,t4

           complex(dp) :: result

      t1 = s26 + s56
      t2 = 2.0_dp
      t3 = 0.1e1_dp / 0.4e1_dp
      t4 = 0.1e1_dp / t1 ** 2.0_dp
      result = t4 * (t3 * (s25 * t2 - t1) + s25 ** 2.0_dp * MBL1011(1) /
     # 2.0_dp)

           intHL0s25s260s56x1031D4eps1 = result
       end function intHL0s25s260s56x1031D4eps1

       function intHL0s25s26s34s56x1110D2eps0()
           implicit none
           complex(dp) :: intHL0s25s26s34s56x1110D2eps0
           complex(dp) ::  t1

           complex(dp) :: result

      t1 = 0.1e1_dp / 0.2e1_dp
      result = t1 * (s34 * MBL1110(1) + (s26 + s56) * MBL1011(1) + MBL10
     #01(0) + 1.0_dp)

           intHL0s25s26s34s56x1110D2eps0 = result
       end function intHL0s25s26s34s56x1110D2eps0

       function intHL0s25s26s34s56x1110D2eps1()
           implicit none
           complex(dp) :: intHL0s25s26s34s56x1110D2eps1

           complex(dp) :: result

      result = 0.1e1_dp / 0.2e1_dp

           intHL0s25s26s34s56x1110D2eps1 = result
       end function intHL0s25s26s34s56x1110D2eps1

       function intHL0s25s26s34s56x1120D2eps0()
           implicit none
           complex(dp) :: intHL0s25s26s34s56x1120D2eps0
           complex(dp) ::  t1

           complex(dp) :: result

      t1 = s25 + s26 - s34 + s56
      t1 = 0.1e1_dp / t1
      result = -2.0_dp * s34 * MBL1110(1) * t1 - (s34 * MBL1110(0) + (s2
     #6 + s56) * MBL1011(1) + MBL1001(0)) * t1

           intHL0s25s26s34s56x1120D2eps0 = result
       end function intHL0s25s26s34s56x1120D2eps0

       function intHL0s25s26s34s56x1120D2eps1()
           implicit none
           complex(dp) :: intHL0s25s26s34s56x1120D2eps1
           complex(dp) ::  t1

           complex(dp) :: result

      t1 = s25 + s26 - s34 + s56
      t1 = 0.1e1_dp / t1
      result = -t1 * (s34 * MBL1110(1) + 1.0_dp)

           intHL0s25s26s34s56x1120D2eps1 = result
       end function intHL0s25s26s34s56x1120D2eps1

       function intHL0s25s26s34s56x1130D4eps0()
           implicit none
           complex(dp) :: intHL0s25s26s34s56x1130D4eps0
           complex(dp) ::  t1

           complex(dp) :: result

      t1 = s34 - s25 - s26 - s56
      t1 = 0.1e1_dp / t1 ** 2.0_dp
      result = -(s25 + s26 - 3.0_dp * s34 + s56) * t1 * ((s26 + s56) * M
     #BL1011(1) + MBL1001(0)) / 4.0_dp + s34 * t1 * (s34 * MBL1110(0) + 
     #1.0_dp) / 2.0_dp + 0.3e1_dp / 0.2e1_dp * s34 ** 2.0_dp * MBL1110(1
     #) * t1

           intHL0s25s26s34s56x1130D4eps0 = result
       end function intHL0s25s26s34s56x1130D4eps0

       function intHL0s25s26s34s56x1130D4eps1()
           implicit none
           complex(dp) :: intHL0s25s26s34s56x1130D4eps1
           complex(dp) ::  t1,t2

           complex(dp) :: result

      t1 = s25 + s26 - s34 + s56
      t2 = -0.1e1_dp / 0.4e1_dp
      t1 = 0.1e1_dp / t1 ** 2.0_dp
      result = t1 * (t2 * (s25 + s26 - 3.0_dp * s34 + s56) + s34 ** 2.0_
     #dp * MBL1110(1) / 2.0_dp)

           intHL0s25s26s34s56x1130D4eps1 = result
       end function intHL0s25s26s34s56x1130D4eps1

       function intHL0s25s26s34s56x1210D2eps0()
           implicit none
           complex(dp) :: intHL0s25s26s34s56x1210D2eps0
           complex(dp) ::  t1

           complex(dp) :: result

      t1 = s25 + s26 - s34 + s56
      t1 = 0.1e1_dp / t1
      result = ((s26 + s56) * MBL1011(1) + (s25 + s26 + s56) * MBL1110(0
     #) + (s25 + s26 + s34 + s56) * MBL1110(1) + MBL1001(0)) * t1

           intHL0s25s26s34s56x1210D2eps0 = result
       end function intHL0s25s26s34s56x1210D2eps0

       function intHL0s25s26s34s56x1210D2eps1()
           implicit none
           complex(dp) :: intHL0s25s26s34s56x1210D2eps1
           complex(dp) ::  t1

           complex(dp) :: result

      t1 = s25 + s26 - s34 + s56
      t1 = 0.1e1_dp / t1
      result = t1 * ((s25 + s26 + s56) * MBL1110(1) + 1.0_dp)

           intHL0s25s26s34s56x1210D2eps1 = result
       end function intHL0s25s26s34s56x1210D2eps1

       function intHL0s25s26s34s56x1220D4eps0()
           implicit none
           complex(dp) :: intHL0s25s26s34s56x1220D4eps0
           complex(dp) ::  t1,t2

           complex(dp) :: result

      t1 = s25 + s26 - s34 + s56
      t2 = s25 + s26 + s56
      t1 = 0.1e1_dp / t1 ** 2.0_dp
      result = -t1 * (((s26 + s56) * MBL1011(1) + MBL1001(0) + 1.0_dp) *
     # (s25 + s26 + s34 + s56) + s34 * (s34 + 5.0_dp * t2) * MBL1110(1))
     # / 2.0_dp - s34 * t2 * MBL1110(0) * t1

           intHL0s25s26s34s56x1220D4eps0 = result
       end function intHL0s25s26s34s56x1220D4eps0

       function intHL0s25s26s34s56x1220D4eps1()
           implicit none
           complex(dp) :: intHL0s25s26s34s56x1220D4eps1
           complex(dp) ::  t1,t2

           complex(dp) :: result

      t1 = s25 + s26 - s34 + s56
      t2 = -0.1e1_dp / 0.2e1_dp
      t1 = 0.1e1_dp / t1 ** 2.0_dp
      result = t1 * (t2 * (s25 + s26 + s34 + s56) - s34 * (s25 + s26 + s
     #56) * MBL1110(1))

           intHL0s25s26s34s56x1220D4eps1 = result
       end function intHL0s25s26s34s56x1220D4eps1

       function intHL0s25s26s34s56x1310D4eps0()
           implicit none
           complex(dp) :: intHL0s25s26s34s56x1310D4eps0
           complex(dp) ::  t1,t2,t3,t4

           complex(dp) :: result

      t1 = s25 + s26 + s56
      t2 = s25 + s26 - s34 + s56
      t3 = 3.0_dp
      t4 = 6.0_dp
      t2 = 0.1e1_dp / t2 ** 2.0_dp
      result = t2 * (((s26 + s56) * MBL1011(1) + MBL1001(0)) * (t1 * t3 
     #- s34) + (t3 * (s25 ** 2.0_dp + s26 ** 2.0_dp + s56 ** 2.0_dp) + t
     #4 * ((s25 + s26) * s56 + s25 * s26) + 4.0_dp * t1 * s34 - s34 ** 2
     #.0_dp) * MBL1110(1)) / 4.0_dp + t1 * t2 * (t1 * MBL1110(0) + 1.0_d
     #p) / 2.0_dp

           intHL0s25s26s34s56x1310D4eps0 = result
       end function intHL0s25s26s34s56x1310D4eps0

       function intHL0s25s26s34s56x1310D4eps1()
           implicit none
           complex(dp) :: intHL0s25s26s34s56x1310D4eps1
           complex(dp) ::  t1,t2,t3,t4

           complex(dp) :: result

      t1 = s25 + s26 + s56
      t2 = 3.0_dp
      t3 = s25 + s26 - s34 + s56
      t4 = 0.1e1_dp / 0.4e1_dp
      t3 = 0.1e1_dp / t3 ** 2.0_dp
      result = t3 * (t4 * (t1 * t2 - s34) + t1 ** 2.0_dp * MBL1110(1) / 
     #2.0_dp)

           intHL0s25s26s34s56x1310D4eps1 = result
       end function intHL0s25s26s34s56x1310D4eps1

       function intHLs160000x0111D0eps0()
           implicit none
           complex(dp) :: intHLs160000x0111D0eps0

           complex(dp) :: result

      result = I300s16(0)

           intHLs160000x0111D0eps0 = result
       end function intHLs160000x0111D0eps0

       function intHLs160000x0111D0eps1()
           implicit none
           complex(dp) :: intHLs160000x0111D0eps1

           complex(dp) :: result

      result = I300s16(1)

           intHLs160000x0111D0eps1 = result
       end function intHLs160000x0111D0eps1

       function intHLs160000x0111D2eps0()
           implicit none
           complex(dp) :: intHLs160000x0111D2eps0

           complex(dp) :: result

      result = 0.3e1_dp / 0.2e1_dp + s16 * I300s16(1) / 2.0_dp

           intHLs160000x0111D2eps0 = result
       end function intHLs160000x0111D2eps0

       function intHLs160000x0111D2eps1()
           implicit none
           complex(dp) :: intHLs160000x0111D2eps1

           complex(dp) :: result

      result = 0.1e1_dp / 0.2e1_dp

           intHLs160000x0111D2eps1 = result
       end function intHLs160000x0111D2eps1

       function intHLs160000x0112D2eps0()
           implicit none
           complex(dp) :: intHLs160000x0112D2eps0

           complex(dp) :: result

      result = -2.0_dp / s16 - I300s16(1)

           intHLs160000x0112D2eps0 = result
       end function intHLs160000x0112D2eps0

       function intHLs160000x0112D2eps1()
           implicit none
           complex(dp) :: intHLs160000x0112D2eps1

           complex(dp) :: result

      result = -0.1e1_dp / s16

           intHLs160000x0112D2eps1 = result
       end function intHLs160000x0112D2eps1

       function intHLs160000x0113D4eps0()
           implicit none
           complex(dp) :: intHLs160000x0113D4eps0

           complex(dp) :: result

      result = -0.1e1_dp / s16 / 2.0_dp - I300s16(1) / 4.0_dp

           intHLs160000x0113D4eps0 = result
       end function intHLs160000x0113D4eps0

       function intHLs160000x0113D4eps1()
           implicit none
           complex(dp) :: intHLs160000x0113D4eps1

           complex(dp) :: result

      result = -0.1e1_dp / s16 / 4.0_dp

           intHLs160000x0113D4eps1 = result
       end function intHLs160000x0113D4eps1

       function intHLs160000x0121D2eps0()
           implicit none
           complex(dp) :: intHLs160000x0121D2eps0

           complex(dp) :: result

      result = 4.0_dp / s16 + I300s16(0) + 2.0_dp * I300s16(1)

           intHLs160000x0121D2eps0 = result
       end function intHLs160000x0121D2eps0

       function intHLs160000x0121D2eps1()
           implicit none
           complex(dp) :: intHLs160000x0121D2eps1

           complex(dp) :: result

      result = 2.0_dp / s16 + I300s16(1)

           intHLs160000x0121D2eps1 = result
       end function intHLs160000x0121D2eps1

       function intHLs160000x0122D4eps0()
           implicit none
           complex(dp) :: intHLs160000x0122D4eps0

           complex(dp) :: result

      result = -0.3e1_dp / 0.2e1_dp / s16 - I300s16(1) / 2.0_dp

           intHLs160000x0122D4eps0 = result
       end function intHLs160000x0122D4eps0

       function intHLs160000x0122D4eps1()
           implicit none
           complex(dp) :: intHLs160000x0122D4eps1

           complex(dp) :: result

      result = -0.1e1_dp / s16 / 2.0_dp

           intHLs160000x0122D4eps1 = result
       end function intHLs160000x0122D4eps1

       function intHLs160000x0211D2eps0()
           implicit none
           complex(dp) :: intHLs160000x0211D2eps0

           complex(dp) :: result

      result = -2.0_dp / s16 - I300s16(1)

           intHLs160000x0211D2eps0 = result
       end function intHLs160000x0211D2eps0

       function intHLs160000x0211D2eps1()
           implicit none
           complex(dp) :: intHLs160000x0211D2eps1

           complex(dp) :: result

      result = -0.1e1_dp / s16

           intHLs160000x0211D2eps1 = result
       end function intHLs160000x0211D2eps1

       function intHLs16s25s26s34s56x1111D2eps0()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1111D2eps0
           complex(dp) ::  t1,t2,t3,t4,t5,t6,t7,t8

           complex(dp) :: result

      t1 = s25 + s26 + s56
      t2 = s13 + s14 + s34
      t3 = 2.0_dp
      t4 = s16 * s25 + (-s34 + s16 + s25 + s26) * s26 + (s26 * t3 + s16 
     #+ s25 - s34 + s56) * s56
      t5 = s26 + s56
      t6 = s25 - s16
      t7 = s25 + s26 - s34 + s56
      t8 = 0.1e1_dp / t4
      t4 = 0.1e1_dp / t4
      t4 = t4 * t1
      t2 = 0.1e1_dp / t2
      result = t8 * (t1 * ((-t1 * MBL1111(0) + I300s16(0)) * s16 + t5 * 
     #MBL1011(0) + t7 * MBL1110(0)) + (-s25 * t6 - s26 * t6 + (s25 - s26
     #) * s34 - (s25 + s34 - s16) * s56) * MBL1101(0)) / 2.0_dp + t4 * (
     #(-t1 * MBL1111(1) + I300s16(1)) * s16 + t5 * MBL1011(1) + t7 * MBL
     #1110(1)) + t4 * t3 * (s13 + s14 - s25 - s26 + s34 - s56) * t2

           intHLs16s25s26s34s56x1111D2eps0 = result
       end function intHLs16s25s26s34s56x1111D2eps0

       function intHLs16s25s26s34s56x1111D2eps1()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1111D2eps1
           complex(dp) ::  t1,t2,t3,t4

           complex(dp) :: result

      t1 = s25 + s26 + s56
      t2 = s13 + s14 + s34
      t3 = 2.0_dp
      t3 = s16 * s25 + (-s34 + s16 + s25 + s26) * s26 + (s26 * t3 + s16 
     #+ s25 - s34 + s56) * s56
      t4 = 0.1e1_dp / t3
      t2 = 0.1e1_dp / t2
      t3 = 0.1e1_dp / t3
      result = t4 * t1 * ((-t1 * MBL1111(1) + I300s16(1)) * s16 + (s26 +
     # s56) * MBL1011(1) + (s25 + s26 - s34 + s56) * MBL1110(1)) / 2.0_d
     #p - t1 * (-s13 - s14 + s25 + s26 - s34 + s56) * t2 * t3

           intHLs16s25s26s34s56x1111D2eps1 = result
       end function intHLs16s25s26s34s56x1111D2eps1

       function intHLs16s25s26s34s56x1112D2eps0()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1112D2eps0
           complex(dp) ::  t1

           complex(dp) :: result

      t1 = s25 + s26 + s56
      t1 = 0.1e1_dp / t1
      result = -0.1e1_dp / s16 * ((s34 * MBL1101(0) - (s34 - s25 - s26 -
     # s56) * MBL1111D2(0)) * t1 - MBL1011(0))

           intHLs16s25s26s34s56x1112D2eps0 = result
       end function intHLs16s25s26s34s56x1112D2eps0

       function intHLs16s25s26s34s56x1112D2eps1()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1112D2eps1

           complex(dp) :: result

      result = MBL1011(1) / s16

           intHLs16s25s26s34s56x1112D2eps1 = result
       end function intHLs16s25s26s34s56x1112D2eps1

       function intHLs16s25s26s34s56x1112D4eps0()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1112D4eps0
           complex(dp) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t2,t3,t4
           complex(dp) ::  t5,t6,t7,t8,t9

           complex(dp) :: result

      t1 = s16 - s25
      t2 = s16 + s25
      t3 = s25 ** 2.0_dp
      t4 = s16 ** 2.0_dp
      t5 = 2.0_dp * s16
      t6 = s34 ** 2.0_dp
      t7 = t4 + t3 + t6
      t8 = t1 ** 2.0_dp
      t9 = t8 * s26
      t10 = -t5 * t3 + t7 * s25 - (s26 * t2 + 2.0_dp * t3) * s34 + (-s16
     # * s34 - s25 * s34 - t5 * s25 + t3 + t4) * s56 + t9
      t11 = s16 * s25
      t7 = -2.0_dp * s34 * t2 - 2.0_dp * t11 + t7
      t11 = (-s34 + t2 + s26) * s26 + (s16 + s25 + 2.0_dp * s26 - s34 + 
     #s56) * s56 + t11
      t12 = 3.0_dp * s25
      t4 = (t4 + t3) * s25
      t13 = s16 * t3
      t14 = t1 * s25
      t1 = t1 * s26
      t15 = s25 + s26 - s34 + s56
      t7 = 0.1e1_dp / t7
      t16 = 0.1e1_dp / t11
      t11 = 0.1e1_dp / t11
      t17 = 0.1e1_dp / 0.2e1_dp
      result = t17 * ((s25 + s26 + s56) * t15 * MBL1111D2(0) * t11 + ((-
     #s25 * ((s16 + 2.0_dp * s25 + s26 - s34) * s34 + (s16 - s25 + s34) 
     #* s56 + t14 + t1) * MBL1101(0) - t15 * ((s25 - s26) * s34 + (s16 -
     # s25 - s34) * s56 + t1 + t14) * MBL1110(1)) * s34 + t10 * (-s16 * 
     #I300s16(1) + MBL1001(0)) + (-2.0_dp * s25 * (t13 - t9) + t8 * s26 
     #** 2.0_dp + (t6 * (s25 + 3.0_dp * s26) + t4) * s25 + ((-(3.0_dp * 
     #s34 + t5) * s25 - (s34 - s16) * s16 + t3) * s56 + (t12 * s34 - (3.
     #0_dp * s16 + 5.0_dp * s25 + 6.0_dp * s26) * s25 - t5 * s26) * s34 
     #+ 2.0_dp * t4 + 2.0_dp * t9 - 4.0_dp * t13) * s56 - (s25 + s26) * 
     #((s16 + t12) * s26 + 2.0_dp * s25 * t2) * s34) * MBL1011(1)) * t16
     # * t7) - t10 * t7 * t16

           intHLs16s25s26s34s56x1112D4eps0 = result
       end function intHLs16s25s26s34s56x1112D4eps0

       function intHLs16s25s26s34s56x1112D4eps1()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1112D4eps1

           complex(dp) :: result

      result = 0.0_dp

           intHLs16s25s26s34s56x1112D4eps1 = result
       end function intHLs16s25s26s34s56x1112D4eps1

       function intHLs16s25s26s34s56x1113D4eps0()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1113D4eps0
           complex(dp) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           complex(dp) ::  t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30
           complex(dp) ::  t31,t32,t33,t34,t35,t36,t37,t38,t4,t5,t6,t7
           complex(dp) ::  t8,t9

           complex(dp) :: result

      t1 = s16 - s25
      t2 = 2.0_dp
      t3 = t2 * s26
      t4 = 3.0_dp * s25
      t5 = s16 + s26
      t6 = s25 ** 2.0_dp
      t7 = s25 * t6
      t8 = s16 ** 2.0_dp
      t9 = s16 * t8
      t10 = s34 ** 2.0_dp
      t11 = t10 + t8 + t6
      t12 = t8 + t6
      t13 = t1 ** 2.0_dp
      t14 = t13 * s26
      t15 = t12 * s25
      t16 = t2 * s16
      t17 = t16 * t6
      t18 = -t17 + (s25 * (s16 - t4 - t3) + (-s34 + t4 + t5) * s34) * s3
     #4 + (-t2 * s25 * (s16 + s34) + t11) * s56 + t14 + t15
      t19 = s16 + s25
      t20 = s16 * s25
      t11 = -t2 * (s34 * t19 + t20) + t11
      t21 = (-s34 + t19 + s26) * s26 + (t3 + s16 + s25 - s34 + s56) * s5
     #6 + t20
      t22 = s26 ** 2.0_dp
      t23 = 4.0_dp * s26
      t24 = s16 * t6
      t25 = s26 + s56
      t26 = s16 * t25
      t27 = 4.0_dp * s16
      t28 = s25 + s26
      t29 = 5.0_dp * t8
      t30 = 3.0_dp * s16
      t31 = s25 + t30
      t32 = t2 * s25
      t30 = t30 + t32
      t33 = 7.0_dp * s16
      t34 = 4.0_dp * s25
      t35 = t1 * s25
      t36 = s16 + t32
      t1 = t1 * s26
      t37 = s25 + s26 - s34 + s56
      t26 = 0.1e1_dp / t26
      t38 = 0.1e1_dp / s16
      t25 = 0.1e1_dp / t25
      t21 = 0.1e1_dp / t21
      t11 = 0.1e1_dp / t11
      t7 = (t2 * (t13 * t19 * t22 + t24 * t12) - ((-(-(3.0_dp * s34 + t3
     #2) * s16 + (-s34 + s25) * s25 + t8) * s56 + t2 * ((t35 + t8) * s25
     # - t9) + s34 * (-s34 * t30 + s16 * (5.0_dp * s16 + 9.0_dp * s26) +
     # (3.0_dp * s26 + t33 + t34) * s25) - 3.0_dp * t14) * s56 - ((-s25 
     #* s34 + 6.0_dp * s16 * t28 + (t4 + t23) * s25) * s34 - 3.0_dp * t2
     #2 * t31 - 3.0_dp * t7 - t3 * ((t33 + t34) * s25 + t29) - t20 * (9.
     #0_dp * s16 + 8.0_dp * s25)) * s34 - t6 ** 2.0_dp - 4.0_dp * s25 * 
     #t9 - 4.0_dp * t14 * t19 - 3.0_dp * t13 * t22 - t24 * (t32 - t33)) 
     #* s56 + t10 * ((s26 * t30 + t4 * (s25 + t16)) * s26 + t17) - s25 *
     # s26 * s34 * t10 + t14 * (s25 * (s25 + t27) + t22) - t28 * (4.0_dp
     # * t20 * t19 + s26 * (t31 * s26 + (t4 + t27) * s25 + t29)) * s34 -
     # 4.0_dp * t8 * t7) * MBL1011(1) * t38
      t9 = t38 * t11
      t12 = 0.1e1_dp / 0.2e1_dp
      result = t12 * (t21 * (-t9 * (((-t1 - t35) * s25 + (-s16 * t19 - s
     #26 * t36 + (s26 + t4 + t16 - s34) * s34 - 3.0_dp * t6) * s34 + ((-
     #t36 + s34) * s34 - t35) * s56) * MBL1101(0) + t37 * ((-t32 - t5 + 
     #s34) * s34 - (s16 - s25 + s34) * s56 - t1 - t35) * MBL1110(1)) * s
     #34 - t37 ** 2.0_dp * MBL1111D2(0) * t38 + t11 * ((-t2 * (s34 * s56
     # ** 2.0_dp + t24) + ((s25 + t3) * s34 - t2 * (s25 * t19 + t22) - 3
     #.0_dp * s26 * t19) * s34 + (t2 * (-t20 + t10) - (3.0_dp * t19 + t2
     #3) * s34 + t6 + t8) * s56 + t14 + t15) * MBL1001(0) + t7) * t25 + 
     #t18 * I300s16(1) * t11) + (s25 + s26 + s56) * MBL1011(0) * t26) + 
     #t9 * t18 * t21

           intHLs16s25s26s34s56x1113D4eps0 = result
       end function intHLs16s25s26s34s56x1113D4eps0

       function intHLs16s25s26s34s56x1113D4eps1()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1113D4eps1
           complex(dp) ::  t1

           complex(dp) :: result

      t1 = s16 * (s26 + s56)
      t1 = 0.1e1_dp / t1
      result = t1 * ((s25 + s26 + s56) * MBL1011(1) + 1.0_dp) / 2.0_dp

           intHLs16s25s26s34s56x1113D4eps1 = result
       end function intHLs16s25s26s34s56x1113D4eps1

       function intHLs16s25s26s34s56x1114D6eps0()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1114D6eps0
           complex(dp) ::  t1,t10,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109
           complex(dp) ::  t11,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t12
           complex(dp) ::  t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t13,t130
           complex(dp) ::  t131,t132,t133,t134,t135,t136,t137,t138,t139,t14,t140,t141
           complex(dp) ::  t142,t143,t144,t145,t146,t147,t148,t149,t15,t150,t151,t152
           complex(dp) ::  t153,t154,t155,t156,t157,t158,t159,t16,t160,t161,t162,t163
           complex(dp) ::  t164,t165,t166,t167,t168,t169,t17,t170,t171,t172,t173,t174
           complex(dp) ::  t175,t176,t177,t178,t179,t18,t180,t181,t182,t183,t184,t185
           complex(dp) ::  t186,t187,t188,t189,t19,t190,t191,t192,t193,t194,t195,t196
           complex(dp) ::  t197,t198,t199,t2,t20,t200,t21,t22,t23,t24,t25,t26
           complex(dp) ::  t27,t28,t29,t3,t30,t31,t32,t33,t34,t35,t36,t37
           complex(dp) ::  t38,t39,t4,t40,t41,t42,t43,t44,t45,t46,t47,t48
           complex(dp) ::  t49,t5,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59
           complex(dp) ::  t6,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t7
           complex(dp) ::  t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t8,t80
           complex(dp) ::  t81,t82,t83,t84,t85,t86,t87,t88,t89,t9,t90,t91
           complex(dp) ::  t92,t93,t94,t95,t96,t97,t98,t99

           complex(dp) :: result

      t1 = s16 - s25
      t2 = 2.0_dp * s16
      t3 = s16 + s25
      t4 = 6.0_dp * s16
      t5 = 5.0_dp * s25
      t6 = t5 + t4
      t7 = s16 ** 2.0_dp
      t8 = t7 ** 2.0_dp
      t9 = t7 * t8
      t10 = s16 * t8
      t11 = s16 * t7
      t12 = 13.0_dp * s16
      t13 = 17.0_dp * s25
      t14 = (t12 + t13) * s25 + 2.0_dp * t7
      t15 = 19.0_dp * s25
      t16 = 5.0_dp * s16
      t17 = t16 + t15
      t18 = s25 + s26
      t19 = 22.0_dp * s16
      t20 = 15.0_dp * s25
      t21 = 19.0_dp * t7
      t22 = 30.0_dp * s25
      t23 = 32.0_dp * s16
      t24 = 11.0_dp * t11
      t25 = 31.0_dp * t7
      t26 = (61.0_dp * s25 + t23) * s25
      t27 = 23.0_dp * s16
      t28 = 27.0_dp * s25
      t29 = t28 + t27
      t30 = s16 * s25
      t31 = s25 ** 2.0_dp
      t32 = t31 ** 2.0_dp
      t33 = t32 ** 2.0_dp
      t34 = s25 * t31
      t35 = t34 * t32
      t36 = t31 * t32
      t37 = s25 * t32
      t38 = t7 * t31
      t39 = 8.0_dp
      t40 = t31 + t7
      t41 = s25 * t40
      t42 = t39 * t7
      t43 = 38.0_dp * s16
      t44 = 75.0_dp * s25
      t45 = 71.0_dp * t7
      t46 = (t44 + t43) * s25 + t45
      t47 = 60.0_dp * t11 * t3 + t31 * t46
      t48 = 3.0_dp
      t49 = 200.0_dp
      t50 = t48 * s16
      t51 = 24.0_dp * t8
      t52 = (((s25 * t49 - t50) * s25 + 58.0_dp * t7) * s25 + 89.0_dp * 
     #t11) * s25
      t53 = t52 + t51
      t54 = 53.0_dp * t11
      t55 = ((s16 + 162.0_dp * s25) * s25 + 24.0_dp * t7) * s25
      t56 = t55 + t54
      t57 = (39.0_dp * s25 + t4) * s25 + t21
      t58 = s26 ** 2.0_dp
      t59 = t58 ** 2.0_dp
      t60 = s26 * t59
      t61 = s26 * t58
      t62 = t57 * s26
      t63 = t7 * t34
      t64 = 12.0_dp
      t65 = 6.0_dp * s26
      t66 = s16 * t31
      t67 = t39 * s16
      t68 = 7.0_dp * t7
      t69 = 100.0_dp * s25
      t70 = ((102.0_dp * s16 + t69) * s25 + 113.0_dp * t7) * s25 + 67.0_
     #dp * t11
      t71 = 107.0_dp * t7
      t72 = 32.0_dp * t11
      t73 = ((149.0_dp * s16 + 230.0_dp * s25) * s25 + t71) * s25
      t74 = t73 + t72
      t75 = 81.0_dp * s16
      t76 = 148.0_dp * s25
      t77 = (t76 + t75) * s25 + 45.0_dp * t7
      t78 = 9.0_dp * s16
      t79 = 25.0_dp * s25
      t80 = -t79 - t78
      t81 = s25 * t70
      t82 = t63 * t3
      t44 = (83.0_dp * s16 + t44) * s25 + 49.0_dp * t7
      t83 = 32.0_dp * t7
      t84 = (110.0_dp * s16 + 145.0_dp * s25) * s25
      t85 = t84 + t83
      t86 = 30.0_dp * s16
      t87 = 67.0_dp * s25
      t88 = t87 + t86
      t89 = s25 * t44
      t90 = t48 * s26
      t91 = t90 * t66
      t92 = t22 + t27
      t93 = 21.0_dp * s16
      t94 = 47.0_dp * s25
      t95 = t94 + t93
      t96 = s25 * t92
      t97 = t66 * t48
      t98 = t48 * s25
      t99 = 10.0_dp * s25
      t100 = 64.0_dp * s16
      t101 = 72.0_dp * t7
      t102 = ((91.0_dp * s25 + t100) * s25 + t101) * s25
      t103 = 55.0_dp * s16
      t104 = 88.0_dp * s25
      t105 = (t104 + t103) * s25 + t25
      t106 = 4.0_dp
      t107 = 2.0_dp * s25
      t108 = t105 * t106
      t109 = (t102 + t24) * t48
      t110 = 15.0_dp * s16
      t111 = t11 * t31
      t112 = s16 * t40
      t113 = t112 * t34
      t114 = t56 * t106
      t115 = 5.0_dp * t57
      t116 = t53 * t48
      t117 = t8 * t31
      t118 = 36.0_dp * t7
      t119 = t48 * t74
      t120 = t77 * t106
      t121 = t48 * t85
      t122 = 18.0_dp
      t123 = t38 * t122
      t124 = 30.0_dp * t59
      t125 = 48.0_dp * s26
      t126 = 9.0_dp * s26
      t127 = t1 ** 2.0_dp
      t128 = t127 ** 2.0_dp
      t129 = t1 * t127
      t130 = t128 * t14
      t131 = s34 ** 2.0_dp
      t132 = t131 ** 2.0_dp
      t133 = t131 * t132
      t134 = s34 * t132
      t135 = s34 * t131
      t136 = t128 * s25
      t137 = 35.0_dp * s26
      t138 = t128 * t61
      t139 = t34 + t11
      t140 = 10.0_dp * s26 * t29
      t141 = 63.0_dp * s16
      t142 = t11 * s25
      t143 = 10.0_dp * t62
      t144 = -10.0_dp * t80
      t145 = 60.0_dp * s26
      t146 = 72.0_dp * s26
      t147 = t122 * s26
      t148 = t8 * t34
      t149 = t128 * t17
      t150 = 5.0_dp * t10
      t151 = t106 * s26
      t152 = 6.0_dp * t133
      t153 = 31.0_dp * s16
      t154 = s16 + s26
      t155 = s16 * t154
      t156 = t32 + t8
      t157 = t30 * t40
      t158 = t127 * t18
      t13 = (-(-28.0_dp * t157 + 6.0_dp * t132 - ((-s34 * t80 - t57) * s
     #34 + t127 * t29) * s34 + 7.0_dp * t156 + 42.0_dp * t38) * s56 + t1
     #34 * t64 - t31 * (((-71.0_dp * s16 + t15) * s25 + 94.0_dp * t7) * 
     #s25 - 46.0_dp * t11) - t150 + s25 * t8 + ((-t115 * s26 + (-(30.0_d
     #p * t154 + t87) * s34 + (125.0_dp * s26 + t75 + t76) * s25 + 45.0_
     #dp * t155) * s34 - t54 - t55) * s34 + t127 * (s16 * (115.0_dp * s2
     #6 + t153) + (135.0_dp * s26 + t103 + t104) * s25)) * s34 - t137 * 
     #t128) * s56 - s25 * ((t31 * ((-t103 + t13) * s25 + 52.0_dp * t7) -
     # 23.0_dp * t8) * s25 + t150) + (((((t125 + t94 + t93) * s34 - (120
     #.0_dp * s16 + 268.0_dp * s25 + t145) * s26 - t83 - t84) * s34 + (s
     #26 * t144 + t120) * s26 + t72 + t73) * s34 - (t143 + t114) * s26 -
     # t51 - t52) * s34 + t127 * ((t140 + t108) * s26 + t102 + t24)) * s
     #34 - 2.0_dp * t11 * t139 - t151 * t149 - 70.0_dp * t128 * t58 - t1
     #52
      t13 = t13 * s56 - t31 * (t31 * (((t5 - t78) * s25 - t42) * s25 + 2
     #2.0_dp * t11) - 13.0_dp * t10) - t48 * (t130 * s26 - t148) - 6.0_d
     #p * s25 * t9 - 6.0_dp * t149 * t58 + (((((-(t5 + t147) * s34 + (14
     #1.0_dp * s25 + t141 + t146) * s26 + t96) * s34 - ((6.0_dp * t88 + 
     #t145) * s26 + t121) * s26 - t89) * s34 + t81 + s26 * (t144 * t58 +
     # t119) + 6.0_dp * t77 * t58) * s34 - ((6.0_dp * t56 + t143) * s26 
     #+ t116) * s26 - 60.0_dp * t142 * t3 - t34 * t46) * s34 + t127 * (3
     #0.0_dp * s25 * t139 + ((6.0_dp * t105 + t140) * s26 + t109) * s26 
     #+ t66 * (t94 + t141))) * s34 - 70.0_dp * t138
      t13 = t13 * s56 - t48 * (s16 * t35 + t130 * t58) + ((-t107 * t47 *
     # s26 + (t107 * t70 * s26 + (-t107 * t44 * s26 - t106 * t88 * t61 +
     # (((t48 * t95 + t125) * s26 + t107 * t92) * s26 + t97) * s34 - t11
     #0 * t34 - t121 * t58 - t123 - t124) * s34 + t119 * t58 + t120 * t6
     #1 + t86 * t32 - 5.0_dp * t80 * t59 + 42.0_dp * t111 + 48.0_dp * t6
     #3) * s34 - t114 * t61 - t115 * t59 - t116 * t58 - t118 * t32 - 30.
     #0_dp * t113 - 48.0_dp * t117) * s34 + t127 * (27.0_dp * t111 + t10
     #7 * (t98 + t2) * ((t99 + t78) * s25 + 15.0_dp * t7) * s26 + t109 *
     # t58 + t108 * t61 + 5.0_dp * t29 * t59 + 30.0_dp * t63 + t110 * t3
     #2)) * s34 - 2.0_dp * s26 * ((t5 + t126) * t133 + t136 * t3 * t6) -
     # 6.0_dp * t38 * (t34 * (-s25 - s16) + t8) + 21.0_dp * t10 * t34 - 
     #t51 * t32 - t138 * (t106 * t17 + t137)
      t15 = s25 * t3
      t44 = t128 * t58
      t46 = t11 * t32
      t51 = 7.0_dp * t128 * t60
      t6 = -t39 * t46 * t40 - (((((s26 * (t61 * t64 + s26 * (s26 * t95 +
     # t96) + t97) - t58 * (t5 + t65) * s34) * s34 - ((s26 * t88 + t85) 
     #* s26 + t89) * t58 - 6.0_dp * t60 - t91 * t6 - 2.0_dp * t63) * s34
     # + t39 * t82 + (((-s26 * t80 + t77) * s26 + t74) * s26 + t81) * t5
     #8 + t66 * t65 * ((t5 + t67) * s25 + t68)) * s34 - t64 * t63 * t40 
     #- (s25 * t47 + ((t56 + t62) * s26 + t53) * s26) * t58 - t66 * (t65
     # * ((6.0_dp * t31 + t42) * s16 + 5.0_dp * t41) + t38 * t39)) * s34
     # + t158 * (t39 * t38 * t3 + (((s26 * t29 + t25 + t26) * s26 + ((t2
     #3 + t22) * s25 + 41.0_dp * t7) * s25 + t24) * s26 + t30 * ((t20 + 
     #t19) * s25 + t21)) * s26)) * s34 - t13 * s56 - t44 * ((-s26 * t17 
     #- t14) * s26 - t15 * t6) + 2.0_dp * t63 * t156 + t91 * t128 * (t2 
     #+ s25) + t51 + t64 * t8 * t37
      t13 = s34 * t3
      t14 = t31 + t7 + t131 - 2.0_dp * t30 - 2.0_dp * t13
      t17 = s26 + s56
      t21 = 2.0_dp * s26
      t24 = (s34 - t3 - s26) * s26 + (-t21 - s16 - s25 + s34 - s56) * s5
     #6 - t30
      t25 = 7.0_dp * s25
      t29 = t39 * s25
      t47 = t29 + s16
      t51 = s16 * t3
      t52 = t64 * s25
      t53 = t52 + t78
      t54 = 6.0_dp * t7
      t55 = 85.0_dp
      t56 = 5.0_dp * t7
      t57 = t64 * t11
      t62 = -t31 + t7
      t70 = 72.0_dp * t34
      t72 = 13.0_dp * s16 * t62
      t73 = t72 + t70
      t74 = t31 * t48 + t7
      t75 = s25 * (((s25 * t55 - t86) * s25 + t56) * s25 + t57)
      t76 = 14.0_dp * s16
      t77 = 40.0_dp * s16
      t80 = t122 * t7
      t81 = t48 * t11
      t83 = ((t77 + t69) * s25 + t80) * s25
      t84 = 27.0_dp * s16
      t85 = 68.0_dp * s25
      t86 = t64 * t7
      t88 = (t84 + t85) * s25 + t86
      t89 = t106 * s25
      t91 = t89 + s16
      t92 = 40.0_dp * s25
      t94 = 65.0_dp * s25
      t95 = 9.0_dp * t7
      t96 = 22.0_dp * s25
      t97 = 7.0_dp * s16
      t102 = t90 + t107
      t104 = t53 / 3.0_dp
      t105 = 38.0_dp * s25
      t108 = 11.0_dp * s16
      t109 = 14.0_dp * t7
      t114 = 9.0_dp * t58 + 9.0_dp * t7
      t115 = 24.0_dp * s16
      t116 = t64 * s26
      t119 = t128 * s26
      t120 = 26.0_dp * s16
      t121 = 36.0_dp * s26
      t130 = 17.0_dp * s16
      t137 = 27.0_dp * s26
      t139 = t64 * s16
      t140 = 32.0_dp * s25
      t70 = (s25 * ((((t29 - t153) * s25 + 44.0_dp * t7) * s25 - 26.0_dp
     # * t11) * s25 + t106 * t8) - 6.0_dp * t134 - (((-(t140 + t139 + t1
     #26) * s34 + s16 * (t139 + t126) + (t84 + t85 + t121) * s25) * s34 
     #- (54.0_dp * t31 + t80) * s26 - t70 - t72) * s34 + t127 * (s16 * (
     #t137 + t67) + (t130 + t105 + t121) * s25)) * s34 + t10 + t126 * t1
     #28) * s56 + 2.0_dp * s25 * t10 + 2.0_dp * t119 * t47 + t133 * t48 
     #+ t34 * (((t25 - t120) * s25 + 34.0_dp * t7) * s25 - 16.0_dp * t11
     #) - (((((t96 + t116 + t78) * s34 - (64.0_dp * s26 + t77 + t94) * s
     #25 - t114 - t115 * s26) * s34 + (t126 * t91 + 2.0_dp * t88) * s26 
     #+ t81 + t83) * s34 - (t147 * t74 + 2.0_dp * t73) * s26 - t75) * s3
     #4 + t127 * (s25 * ((t108 + t105) * s25 + t109) + ((34.0_dp * s16 +
     # 76.0_dp * s25) * s25 + 16.0_dp * t7 + t126 * t104) * s26)) * s34 
     #- t117 + 9.0_dp * t44
      t72 = t10 * t31
      t84 = s56 ** 2.0_dp
      t85 = s16 * t37
      t47 = t48 * ((-t106 * t157 - (t104 * t127 + t131 * (t91 - s34)) * 
     #s34 + t32 + t8 + 2.0_dp * t74 * t131 + 6.0_dp * t38) * s56 * t84 +
     # t138) + t70 * s56 + t131 * (t34 * ((t22 - t76) * s25 + t42) + s26
     # * (s26 * t73 + t75) + 6.0_dp * t74 * t61) + t132 * (s25 * ((t22 +
     # t93) * s25 + t42) + ((t106 * (t29 + t50) + t90) * s26 + (t77 + t9
     #4) * s25 + t95) * s26) + t133 * t102 - t134 * (s25 * (t52 + t97) +
     # (t96 + t78 + t65) * s26) - t135 * (t31 * ((t76 + t92) * s25 + t86
     #) + t48 * (t61 * t91 + t142) + s26 * (s26 * t88 + t81 + t83)) - t1
     #58 * (s25 * ((t52 + t50) * s25 + t54) + s26 * (s26 * t53 + t39 * t
     #51 + 26.0_dp * t31)) * s34 + t136 * (t25 + t2) * s26 + t44 * t47 +
     # t72 + 2.0_dp * t34 * ((t34 - t11) * s25 - t8) + t85 * (-t25 + t67
     #)
      t53 = 13.0_dp * s25
      t67 = t53 + t67
      t70 = 14.0_dp * s25
      t73 = t70 + t50
      t5 = t5 + t97
      t74 = 24.0_dp * s25
      t75 = t19 + t74
      t77 = t48 * t61
      t83 = 20.0_dp * s25
      t88 = 45.0_dp * t11
      t26 = (t26 + 42.0_dp * t7) * s25 + t88
      t91 = 100.0_dp * t7
      t94 = ((73.0_dp * s16 + 126.0_dp * s25) * s25 + t91) * s25
      t105 = t94 + t88
      t121 = 57.0_dp * s25
      t141 = 48.0_dp * s16
      t142 = 47.0_dp * t7
      t143 = (t141 + t121) * s25
      t144 = t143 + t142
      t149 = s25 * s26
      t150 = (132.0_dp * s16 + 118.0_dp * s25) * s25 + 97.0_dp * t7
      t153 = 57.0_dp * t7
      t159 = (130.0_dp * s16 + 116.0_dp * s25) * s25
      t160 = t159 + t153
      t161 = 53.0_dp * s25
      t162 = t141 + t161
      t163 = s25 * t150
      t164 = t106 * t31
      t165 = 52.0_dp * s16
      t166 = t165 + t121
      t167 = 52.0_dp * s25
      t168 = 36.0_dp * s16
      t169 = t168 + t167
      t170 = 11.0_dp * s25
      t171 = 60.0_dp * s16
      t87 = (t87 + t171) * s25
      t172 = t87 + t80
      t173 = 88.0_dp * s16
      t174 = 108.0_dp * s25
      t175 = t174 + t173
      t176 = t48 * t58
      t177 = s25 * ((21.0_dp * s25 + t19) * s25 + t109)
      t178 = 15.0_dp * t59
      t179 = 58.0_dp * t11
      t180 = t105 * t48
      t181 = t48 * t160
      t182 = 45.0_dp * t59
      t183 = 84.0_dp * t7
      t184 = t3 * s26
      t185 = 120.0_dp * t184
      t186 = 90.0_dp * s26
      t187 = 44.0_dp * s16
      t188 = t106 * t11
      t189 = 30.0_dp * s26
      t190 = 2.0_dp * t30
      t87 = (t122 * t132 + 5.0_dp * t156 - 20.0_dp * t157 + ((-(45.0_dp 
     #* s26 + t161 + t141) * s34 + t142 + t143 + 60.0_dp * t184) * s34 -
     # t127 * (15.0_dp * s26 + t28 + t19)) * s34 + 30.0_dp * t38) * s56 
     #+ ((((-53.0_dp * s16 + t70) * s25 + t101) * s25 - 38.0_dp * t11) *
     # s25 + 2.0_dp * t8) * s25 + t10 * t48 - 9.0_dp * t134 - (((-(t168 
     #+ t167 + t146) * s34 + (192.0_dp * s16 + 212.0_dp * s25 + t186) * 
     #s26 + t153 + t159) * s34 - (t106 * t144 + t185) * s26 - t88 - t94)
     # * s34 + t127 * ((t174 + t189 + t173) * s26 + t80 + t87)) * s34 + 
     #20.0_dp * t119
      t53 = t87 * s56 + s25 * (((((-t187 + t53) * s25 + 46.0_dp * t7) * 
     #s25 - t188) * s25 - 19.0_dp * t8) * s25 + t10 * t39) + t48 * t127 
     #* (-(10.0_dp * t61 + s26 * (t175 * s26 / 2.0_dp + t172) + t177) * 
     #s34 + t127 * t73 * s26) + ((((-t137 - t170) * s34 + (156.0_dp * s2
     #6 + t121 + t165) * s25 + 108.0_dp * s26 * t154) * s34 - ((6.0_dp *
     # t162 + t186) * s26 + t181) * s26 - t163) * s34 + s25 * (((122.0_d
     #p * s25 + t100) * s25 + t183) * s25 + 90.0_dp * t11) + ((6.0_dp * 
     #t144 + t185) * s26 + t180) * s26) * t131 + 30.0_dp * t44
      t53 = t53 * s56 - t106 * t37 * t62 - (((t107 * t150 * s26 + t106 *
     # t162 * t61 + (((t96 + t137) * s26 + t164) * s34 - ((t64 * t169 / 
     #4.0_dp + t146) * s26 + t107 * t166) * s26 - t31 * (t83 + t27)) * s
     #34 + t171 * t34 + t181 * t58 + t182 + 40.0_dp * t32 + 52.0_dp * t3
     #8) * s34 - t106 * s26 * (s25 * t26 + t144 * t58) - t31 * (((42.0_d
     #p * s16 + t92) * s25 + t118) * s25 + t179) - t58 * (60.0_dp * t3 *
     # t58 + t180)) * s34 + t127 * (t168 * t34 + t176 * t172 + t175 * t6
     #1 + t177 * t65 + t178 + 20.0_dp * t32 + 32.0_dp * t38)) * s34 + t1
     #19 * t107 * t67 + t44 * t48 * t73 + 20.0_dp * t138 + t66 * (((-9.0
     #_dp * t31 + 26.0_dp * t7) * s25 - 24.0_dp * t11) * s25 + 7.0_dp * 
     #t8)
      t62 = s16 * t34
      t87 = t7 * t32
      t88 = t48 * s34
      t94 = 5.0_dp * t128 * t59
      t13 = -t39 * t87 * t40 + t53 * s56 + t131 * (t64 * (t3 * t60 + t11
     #3) + 2.0_dp * t149 * (s25 * (((t83 + t93) * s25 + t80) * s25 + 29.
     #0_dp * t11) + s26 * t26) + t61 * (s26 * t144 + t105) + t42 * t32) 
     #+ t132 * (t122 * t59 + s26 * ((s25 * t166 + s26 * t169) * s26 + t3
     #1 * (t83 + t27)) + t2 * t34) - t135 * ((9.0_dp * t59 + t164 * ((t1
     #10 + t99) * s25 + 13.0_dp * t7)) * s26 + t39 * t51 * t34 + ((s26 *
     # t162 + t160) * s26 + t163) * t58) + t136 * t67 * t58 - t158 * ((t
     #77 + t89 * (s25 * t5 + t54)) * s26 + t39 * t66 * t3 + t58 * (t75 *
     # s26 + (43.0_dp * s25 + t43) * s25 + t80)) * s34 - s26 * ((t126 + 
     #t170) * s26 + t164) * t134 + t138 * t73 + t119 * t31 * (t97 + t89)
     # + 2.0_dp * t62 * t156 - t88 * (-t106 * t13 + t131 * t48 - t190 + 
     #t40) * s56 * t84 ** 2.0_dp + t94 + t57 * t37
      t26 = s25 + s26 + s56
      t27 = (t187 + t98) * s25 + t118
      t43 = t106 * t7
      t51 = (t108 + t107) * s25 + t43
      t53 = 9.0_dp * s25
      t57 = (t19 + t53) * s25 + t48 * t7
      t67 = t75 / 2.0_dp
      t73 = 55.0_dp * s25
      t75 = 35.0_dp * t7
      t80 = t122 * s25
      t84 = 50.0_dp * t7
      t93 = 45.0_dp * s25
      t96 = 83.0_dp * t8
      t100 = 25.0_dp * t7
      t101 = (t139 + t170) * s25 + t100
      t105 = t101 * s26
      t113 = 99.0_dp * t11
      t118 = t76 + t89
      t121 = 66.0_dp * s25
      t137 = 59.0_dp * t7
      t141 = (140.0_dp * s16 + t121) * s25 + t137
      t142 = t97 + t98
      t143 = s25 * ((279.0_dp * s16 + t69) * s25 + 184.0_dp * t7)
      t144 = t31 * ((212.0_dp * s16 + t93) * s25 + 186.0_dp * t7)
      t146 = t106 * t142
      t150 = t16 + t107
      t153 = 16.0_dp * s25
      t159 = 25.0_dp * s16
      t160 = t159 + t153
      t161 = 6.0_dp * s25
      t162 = t97 + t161
      t163 = s26 * t162
      t164 = s25 * t160
      t165 = 35.0_dp * t11
      t166 = 6.0_dp * t101 * s26
      t169 = 177.0_dp * t11
      t171 = ((355.0_dp * s16 + 144.0_dp * s25) * s25 + 336.0_dp * t7) *
     # s25
      t172 = t171 + t113
      t109 = (t130 + t29) * s25 + t109
      t174 = t118 / 2.0_dp
      t175 = t172 * t58
      t177 = s25 * (((555.0_dp * s16 + 160.0_dp * s25) * s25 + 568.0_dp 
     #* t7) * s25 + 315.0_dp * t11)
      t180 = t141 * t106
      t181 = t18 ** 2.0_dp
      t45 = t127 * t181 * (t58 * (t64 * t5 * s26 + (262.0_dp * s16 + 120
     #.0_dp * s25) * s25 + 250.0_dp * t7) + t151 * (((77.0_dp * s16 + t5
     #3) * s25 + t45) * s25 + t165) + t30 * ((118.0_dp * s16 + t73) * s2
     #5 + t71))
      t22 = s34 * (((((((t163 * t39 + 6.0_dp * t164) * s26 + t31 * t150 
     #* t122) * s26 + 11.0_dp * t62) * s34 - t103 * t32 - 20.0_dp * t142
     # * t59 - t143 * t176 - t144 * t21 - t180 * t61 - 68.0_dp * t63) * 
     #s34 + t106 * s26 * (t31 * (((199.0_dp * s16 + t22) * s25 + 238.0_d
     #p * t7) * s25 + t169) + t175) + t177 * t176 + t124 * t109 + t64 * 
     #t174 * t60 + 110.0_dp * t85 + 162.0_dp * t11 * t34 + 184.0_dp * t8
     #7) * s34 - t18 * (s25 * (t21 * ((((317.0_dp * s16 + t93) * s25 + t
     #49 * t7) * s25 + 229.0_dp * t11) * s25 + 257.0_dp * t8) + t30 * ((
     #(144.0_dp * s16 + 110.0_dp * s25) * s25 + 118.0_dp * t7) * s25 + 1
     #88.0_dp * t11)) + ((t166 + ((433.0_dp * s16 + 294.0_dp * s25) * s2
     #5 + 350.0_dp * t7) * s25 + 475.0_dp * t11) * s26 + (((839.0_dp * s
     #16 + 330.0_dp * s25) * s25 + 518.0_dp * t7) * s25 + 781.0_dp * t11
     #) * s25 + 332.0_dp * t8) * t58)) * s34 + t45)
      t45 = 20.0_dp * t7
      t49 = 24.0_dp * t9
      t53 = 100.0_dp * t8
      t71 = 101.0_dp * s16
      t103 = t189 * t5
      t185 = (((318.0_dp * s16 + 156.0_dp * s25) * s25 + 217.0_dp * t7) 
     #* s25 + 314.0_dp * t11) * s25
      t186 = t185 + t96
      t191 = 72.0_dp * s25
      t192 = 95.0_dp * t11
      t91 = ((t191 + t71) * s25 + t91) * s25
      t193 = t91 + t192
      t194 = 10.0_dp * t193
      t195 = s25 * ((((491.0_dp * s16 + 140.0_dp * s25) * s25 + 306.0_dp
     # * t7) * s25 + 413.0_dp * t11) * s25 + 282.0_dp * t8)
      t142 = 40.0_dp * t142
      t196 = ((229.0_dp * s16 + 84.0_dp * s25) * s25 + 196.0_dp * t7) * 
     #s25
      t197 = 43.0_dp * s16
      t100 = (t197 + t74) * s25 + t100
      t198 = t174 * s26
      t199 = 24.0_dp * s26
      t200 = 28.0_dp * s16
      t9 = ((t101 * t131 + t156 * t48 - t157 * t64 + t123 - 2.0_dp * s34
     # * (t127 * t5 + t131 * t174)) * s56 + (t31 * ((-37.0_dp * s16 + t5
     #2) * s25 + 28.0_dp * t7) - 32.0_dp * t8) * s25 + (((-6.0_dp * (t13
     #0 + t29 + t151) * s25 - 84.0_dp * t155 + t146 * s34) * s34 + t166 
     #+ t192 + t91) * s34 - 2.0_dp * t127 * (s16 * (42.0_dp * s26 + t159
     #) + (t197 + t74 + t189) * s25)) * s34 + t122 * (t119 + t111) + 11.
     #0_dp * t10) * s56 + (((((-t200 + t80) * s25 - 62.0_dp * t7) * s25 
     #+ 168.0_dp * t11) * s25 - 122.0_dp * t8) * s25 + 20.0_dp * t10) * 
     #s25 + 6.0_dp * t9 + ((((140.0_dp * s16 * t18 + (t121 + t145) * s25
     # + t137) * s34 - 30.0_dp * s26 * (t109 + t198) - t113 - t171) * s3
     #4 + (5.0_dp * t193 + 15.0_dp * t105) * s26 + t185 + t96) * s34 - t
     #127 * ((10.0_dp * t100 + t103) * s26 + t165 + t196)) * s34 + 5.0_d
     #p * t119 * t67 - 2.0_dp * t162 * t134 + 45.0_dp * t44
      t9 = t9 * s56 + t35 * t64 + (((((t142 * s26 + t180) * s26 + t143) 
     #* s34 - ((60.0_dp * t109 + 40.0_dp * t198) * s26 + t106 * t172) * 
     #s26 - t177) * s34 + ((t194 + 20.0_dp * t105) * s26 + t106 * t186) 
     #* s26 + t195) * s34 - t127 * (s25 * (((311.0_dp * s16 + 64.0_dp * 
     #s25) * s25 + 312.0_dp * t7) * s25 + 129.0_dp * t11) + ((40.0_dp * 
     #t5 * s26 + 20.0_dp * t100) * s26 + t106 * (t196 + t165)) * s26)) *
     # s34 + 10.0_dp * t44 * t67 + t119 * t39 * t57 - (2.0_dp * (t159 + 
     #t153 + t199) * s25 + 2.0_dp * t200 * s26) * t134 + t89 * t133 + s2
     #5 * t49 + 60.0_dp * t138 + t122 * s16 * t36 + 252.0_dp * t46 - 108
     #.0_dp * t148 - 168.0_dp * t7 * t37 - 30.0_dp * t72
      t9 = t9 * s56 + t48 * (t134 * (-((t110 + t161) * s25 + t21 * t160)
     # * s25 - t106 * t162 * t58 + s25 * (t151 + s25) * s34) + t33) - ((
     #(398.0_dp * s16 * t32 - (6.0_dp * t141 * t58 + t142 * t61 + t143 *
     # t90 + 45.0_dp * t32 + 186.0_dp * t38 + 212.0_dp * t62) * s34 + 60
     #.0_dp * t109 * t61 + t124 * t174 + t177 * t90 + 354.0_dp * t111 + 
     #6.0_dp * t175 + 60.0_dp * t37 + 476.0_dp * t63) * s34 - t178 * t10
     #1 - 288.0_dp * t11 * t34 - 6.0_dp * t186 * t58 - t194 * t61 - t195
     # * t90 - 351.0_dp * t117 - 45.0_dp * t36 - 372.0_dp * t85 - 272.0_
     #dp * t87) * s34 + t158 * (s25 * (((209.0_dp * s16 + t80) * s25 + 2
     #60.0_dp * t7) * s25 + t169) + ((10.0_dp * (65.0_dp * s16 + 33.0_dp
     # * s25) * s25 + 10.0_dp * t84 + t103) * s26 + 2.0_dp * (t50 + s25)
     # * ((87.0_dp * s25 + t71) * s25 + t75)) * s26)) * s34 + t147 * t13
     #6 * t51 + 10.0_dp * t138 * t67 + t44 * t64 * t57 + t182 * t128 + t
     #66 * (((((-122.0_dp * s16 + t140) * s25 + 108.0_dp * t7) * s25 + 4
     #3.0_dp * t11) * s25 - t53) * s25 + 36.0_dp * t10)
      t29 = t18 * (24.0_dp * t82 + ((t118 * s26 + (44.0_dp * s25 + t173)
     # * s25 + t183) * s26 + ((267.0_dp * s16 + t69) * s25 + 252.0_dp * 
     #t7) * s25 + t113) * t61 + t149 * (t190 * ((80.0_dp * s16 + t73) * 
     #s25 + 69.0_dp * t7) + t151 * (((72.0_dp * s16 + t20) * s25 + 79.0_
     #dp * t7) * s25 + 54.0_dp * t11))) * t135
      t38 = t181 * ((((77.0_dp * t31 + 95.0_dp * t7) * s16 + 50.0_dp * t
     #41 + t105) * s26 + (((152.0_dp * s16 + t93) * s25 + 92.0_dp * t7) 
     #* s25 + 124.0_dp * t11) * s25 + t96) * t58 + 36.0_dp * t38 * t40 +
     # t30 * (t21 * (((t168 + t73) * s25 + t75) * s25 + t179) + 24.0_dp 
     #* t38)) * t131
      t5 = (-t134 * (t31 * (t126 * t150 + 11.0_dp * t30) + 2.0_dp * t58 
     #* (t164 + t163)) + t62 * t128 * (t170 + t115)) * s26 + (t9 * s56 +
     # t122 * t44 * (s25 * t51 + t61) - t22 + t94 * t67 + t138 * t39 * t
     #57 + t108 * t33 - t45 * t35 + t49 * t34 - t55 * t10 * t32 + t53 * 
     #t37 - 30.0_dp * t11 * t36 + t149 * (t107 * t128 * t27 + t152 * (t2
     #1 + s25))) * s56 + t132 * (s26 * (((s26 * t141 + t143) * s26 + t14
     #4) * s26 + t62 * (68.0_dp * s16 + t73)) + t146 * t60 + 6.0_dp * t8
     #7) - t127 * t18 * t181 * (s26 * (s16 * ((46.0_dp * s16 + t73) * s2
     #5 + t75) + s26 * ((t187 + t80) * s25 + t84)) + 24.0_dp * t15 * t7 
     #+ 2.0_dp * t5 * t61) * s34 - t29 + t128 * t67 * t60 + t38 + s25 * 
     #(t133 * (t151 + t98) + t136 * t27) * t58 + 6.0_dp * s25 * (t138 * 
     #t51 + t63 * (t32 + t8)) - 24.0_dp * t11 * t37 * t40 + 36.0_dp * t8
     # * t36 + t128 * t59 * (2.0_dp * t57 + t176)
      t9 = t48 * t155
      t10 = 2.0_dp * t11
      t15 = s16 * t106
      t21 = (t15 + t161) * s25 + t7
      t22 = ((33.0_dp * s16 + t92) * s25 + t45) * s25
      t27 = (20.0_dp * s16 + t74) * s25 + t54
      t29 = t107 + s16
      t33 = 2.0_dp * t29
      t35 = s25 * t1
      t36 = -t34 + t11
      t37 = t129 * t58
      t38 = t31 * t36 + t37
      t39 = t129 * s26
      t40 = t36 * s25
      t41 = t62 * t1
      t44 = t66 * t1
      t11 = ((t44 * t48 + ((t131 + t21) * s34 + t35 * t104) * s34 - t40 
     #- t33 * t135) * s56 - t48 * (t40 + t39) * s25 + ((t131 * (t170 + t
     #4 + t90) + ((26.0_dp * s25 + t12) * s25 + t54) * s25 + t10 + t90 *
     # t21) * s34 + t35 * (s16 * (t126 + t50) + (t15 + t116 + t70) * s25
     #)) * s34 - 2.0_dp * t135 * (t131 + (10.0_dp * s16 + t65 + t52) * s
     #25 + t9) + 9.0_dp * t41) * s56 - t48 * t38 * s25 + ((((t131 + 22.0
     #_dp * s25 * t154 + (t139 + t90) * s26 + 25.0_dp * t31 + t54) * s34
     # - (t65 * t29 + 2.0_dp * t27) * s26 - t188 - t22) * s34 + s25 * ((
     #(35.0_dp * s25 + t97) * s25 + t86) * s25 + 6.0_dp * t11) + (((t167
     # + t120) * s25 + t86) * s25 + t188) * s26 + t8 + t176 * t21) * s34
     # + t35 * t18 * (s16 * (t126 + t4) + (t153 + t116 - s16) * s25)) * 
     #s34 + 9.0_dp * s16 * t32 * t1 - t106 * (t107 + t154) * t134 - t65 
     #* t129 * t31
      t29 = t36 * t34
      t32 = t129 * t61
      t34 = t149 * t129
      t37 = t62 * (-s16 + s25) + t37 + t34
      t45 = s25 + s26 - s34 + s56
      t46 = -t89 + s16
      t49 = -2.0_dp * t46
      t51 = 42.0_dp * s25
      t52 = 21.0_dp * t184
      t3 = ((((-t88 + t49) * s34 + 7.0_dp * t1 * t3) * s34 + 6.0_dp * t3
     #0 * t1 - 2.0_dp * t36) * s56 + t122 * t44 + 6.0_dp * t132 - 6.0_dp
     # * t39 - 6.0_dp * t40 + ((-(t74 + t126 + t97) * s34 - s16 * (t15 +
     # t65) + (36.0_dp * s25 + t199 + t76) * s25) * s34 + t1 * ((t12 + t
     #74) * s25 + t52 + t56)) * s34) * s56 + t122 * t41 + t48 * t132 * (
     #-s34 + t161 + t151 + t50) + ((-((t23 + t51 + t125) * s25 + t114 + 
     #t76 * s26) * s34 + ((19.0_dp * s16 + 48.0_dp * s25) * s25 + t43) *
     # s25 + (-t65 * t46 + (t200 + t191) * s25 - t42) * s26 + t81) * s34
     # + t1 * t18 * ((t16 + t28) * s25 + 10.0_dp * t7 + t52)) * s34 - 6.
     #0_dp * t38 - t34 * t64
      t16 = 0.1e1_dp / s16
      t14 = 0.1e1_dp / t14 ** 2.0_dp
      t17 = 0.1e1_dp / t17 ** 2.0_dp
      t24 = 0.1e1_dp / t24 ** 2.0_dp
      t28 = t45 ** 2.0_dp
      t2 = t16 * (t24 * (-t26 * t45 * t28 * MBL1111D2(0) - s34 * (-t48 *
     # t31 * t37 + (-t32 - t29) * s25 + ((((t131 * t18 + s25 * ((t20 + t
     #12) * s25 + t54) + ((t170 + t4 + s26) * s26 + (t19 + t79) * s25 + 
     #t54) * s26) * s34 - s25 * (((t83 + t139) * s25 + 11.0_dp * t7) * s
     #25 + t188) - s26 * (s26 * t27 + t188 + t22) - t33 * t61) * s34 + t
     #18 * (s25 * (((t20 - t2) * s25 + t68) * s25 + t188) + (t21 * s26 +
     # ((t83 + t78) * s25 + t56) * s25 + t10) * s26 + t8)) * s34 + t35 *
     # t181 * ((t161 + t151 - t2) * s25 + t9)) * s34 + t11 * s56 - 2.0_d
     #p * t18 * (t98 + t2 + s26) * t134) * MBL1101(0) * t14) + t26 ** 2.
     #0_dp * MBL1011(0) * t17)
      result = t14 * t24 * (t47 * I300s16(1) + (t13 * MBL1001(0) + (t5 *
     # MBL1011(1) + t6) * t16) * t17 - s34 * t45 * (s34 * ((((-s34 * t10
     #2 + s25 * (t99 + t97) + (t80 + t78 + t65) * s26) * s34 - s25 * ((t
     #110 + t83) * s25 + t42) - t77 - s26 * ((t74 + t97) * s26 + (t23 + 
     #t51) * s25 + t95)) * s34 + t18 * ((20.0_dp * t31 + t42) * s25 + (t
     #49 * s26 - t106 * ((-t15 - t25) * s25 + t7)) * s26 + t112 * t48)) 
     #* s34 + t1 * t181 * ((t99 - s16) * s25 + 7.0_dp * t184 + t56)) + t
     #3 * s56 - 2.0_dp * t32 - 2.0_dp * t29 - 6.0_dp * s25 * t37) * MBL1
     #110(1) * t16) / 12.0_dp + t2 / 6.0_dp

           intHLs16s25s26s34s56x1114D6eps0 = result
       end function intHLs16s25s26s34s56x1114D6eps0

       function intHLs16s25s26s34s56x1114D6eps1()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1114D6eps1
           complex(dp) ::  t1,t2,t3

           complex(dp) :: result

      t1 = s26 + s56
      t2 = s25 + s26 + s56
      t3 = 0.1e1_dp / t1 ** 2.0_dp
      result = 0.1e1_dp / s16 * t3 * (s25 / 6.0_dp + t1 / 4.0_dp + t2 **
     # 2.0_dp * MBL1011(1) / 6.0_dp)

           intHLs16s25s26s34s56x1114D6eps1 = result
       end function intHLs16s25s26s34s56x1114D6eps1

       function intHLs16s25s26s34s56x1121D2eps0()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1121D2eps0
           complex(dp) ::  t1,t2,t3

           complex(dp) :: result

      t1 = s25 + s26 + s56
      t2 = -s16 + s25
      t3 = 0.1e1_dp / s16
      t1 = 0.1e1_dp / t1
      result = 2.0_dp * s25 * s34 * MBL1101(0) * t3 * t1 ** 2.0_dp - t1 
     #* ((s25 * MBL1011(0) + s34 * MBL1110(0) + (-s25 * t2 - s26 * t2 + 
     #(s25 - s26) * s34 - (-s16 + s25 + s34) * s56) * MBL1111D2(0) * t1)
     # * t3 - I300s16(0))

           intHLs16s25s26s34s56x1121D2eps0 = result
       end function intHLs16s25s26s34s56x1121D2eps0

       function intHLs16s25s26s34s56x1121D2eps1()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1121D2eps1
           complex(dp) ::  t1

           complex(dp) :: result

      t1 = s25 + s26 + s56
      t1 = 0.1e1_dp / t1
      result = -t1 * ((s25 * MBL1011(1) + s34 * MBL1110(1)) / s16 - I300
     #s16(1))

           intHLs16s25s26s34s56x1121D2eps1 = result
       end function intHLs16s25s26s34s56x1121D2eps1

       function intHLs16s25s26s34s56x1121D4eps0()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1121D4eps0
           complex(dp) ::  t1,t2,t3,t4

           complex(dp) :: result

      t1 = s25 + s26 + s56
      t2 = 2.0_dp
      t2 = s16 * s25 + (-s34 + s16 + s25 + s26) * s26 + (s26 * t2 + s16 
     #+ s25 - s34 + s56) * s56
      t3 = s25 - s16
      t4 = 0.1e1_dp / t2
      t2 = 0.1e1_dp / t2
      result = t1 * t2 + t4 * (t1 * (s16 * I300s16(1) - s34 * MBL1110(1)
     # - t1 * MBL1011(1) - MBL1001(0)) + (-s25 * t3 - s26 * t3 + (s25 - 
     #s26) * s34 - (s25 - s16 + s34) * s56) * MBL1111D2(0) - s25 * s34 *
     # MBL1101(0)) / 2.0_dp

           intHLs16s25s26s34s56x1121D4eps0 = result
       end function intHLs16s25s26s34s56x1121D4eps0

       function intHLs16s25s26s34s56x1121D4eps1()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1121D4eps1

           complex(dp) :: result

      result = 0.0_dp

           intHLs16s25s26s34s56x1121D4eps1 = result
       end function intHLs16s25s26s34s56x1121D4eps1

       function intHLs16s25s26s34s56x1122D4eps0()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1122D4eps0
           complex(dp) ::  t1,t10,t11,t12,t13,t14,t2,t3,t4,t5,t6,t7
           complex(dp) ::  t8,t9

           complex(dp) :: result

      t1 = s25 + s26 - s34 + s56
      t2 = s16 + s25
      t3 = 2.0_dp
      t4 = s16 * s25 + (-s34 + t2 + s26) * s26 + (s26 * t3 + s16 + s25 -
     # s34 + s56) * s56
      t5 = s25 + s26 + s56
      t6 = s26 + s56
      t7 = s16 * t6
      t8 = 3.0_dp * s16
      t9 = s56 ** 2.0_dp
      t10 = s26 ** 2.0_dp
      t11 = -s16 + s25
      t12 = t11 * s26
      t4 = 0.1e1_dp / t4
      t13 = 0.1e1_dp / t5
      t7 = 0.1e1_dp / t7
      t6 = 0.1e1_dp / t6
      t14 = 0.1e1_dp / s16
      result = -t3 * t1 * t14 * t4 + (t3 * ((-s34 * (s25 + s26) + t12) *
     # s25 + ((-s16 - s34 + s25) * s25 + t12) * s56) + (s25 * t11 + s34 
     #** 2.0_dp) * s25 + t11 * (t9 + t10)) * MBL1111D2(0) * t14 * t13 * 
     #t4 - t4 * (-s34 * t14 * (s25 * MBL1101(0) * t13 + MBL1110(1)) + I3
     #00s16(1)) * t1 - t4 * (t5 * MBL1001(0) + (t3 * s16 * s25 ** 2.0_dp
     # + (t3 * t2 * s26 + s25 * (s25 - s34) + t8 * s25) * s56 + t2 * t9 
     #+ t2 * t10 - s25 * (s34 - s25 - t8) * s26) * MBL1011(1) * t14) * t
     #6 - s25 * MBL1011(0) * t7

           intHLs16s25s26s34s56x1122D4eps0 = result
       end function intHLs16s25s26s34s56x1122D4eps0

       function intHLs16s25s26s34s56x1122D4eps1()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1122D4eps1
           complex(dp) ::  t1

           complex(dp) :: result

      t1 = s16 * (s26 + s56)
      t1 = 0.1e1_dp / t1
      result = -t1 * (s25 * MBL1011(1) + 1.0_dp)

           intHLs16s25s26s34s56x1122D4eps1 = result
       end function intHLs16s25s26s34s56x1122D4eps1

       function intHLs16s25s26s34s56x1123D6eps0()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1123D6eps0
           complex(dp) ::  t1,t10,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109
           complex(dp) ::  t11,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t12
           complex(dp) ::  t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t13,t130
           complex(dp) ::  t131,t132,t133,t134,t135,t136,t137,t138,t139,t14,t140,t141
           complex(dp) ::  t142,t143,t144,t145,t146,t147,t148,t15,t16,t17,t18,t19
           complex(dp) ::  t2,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3
           complex(dp) ::  t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t4,t40
           complex(dp) ::  t41,t42,t43,t44,t45,t46,t47,t48,t49,t5,t50,t51
           complex(dp) ::  t52,t53,t54,t55,t56,t57,t58,t59,t6,t60,t61,t62
           complex(dp) ::  t63,t64,t65,t66,t67,t68,t69,t7,t70,t71,t72,t73
           complex(dp) ::  t74,t75,t76,t77,t78,t79,t8,t80,t81,t82,t83,t84
           complex(dp) ::  t85,t86,t87,t88,t89,t9,t90,t91,t92,t93,t94,t95
           complex(dp) ::  t96,t97,t98,t99

           complex(dp) :: result

      t1 = 2.0_dp
      t2 = t1 * s25
      t3 = s16 - s25
      t4 = 3.0_dp
      t5 = t4 * s25
      t6 = s16 - t5
      t7 = s25 + s26
      t8 = s16 ** 2.0_dp
      t9 = t8 ** 2.0_dp
      t10 = s16 * t8
      t11 = t4 * s16
      t12 = 18.0_dp * s25
      t13 = s16 + s25
      t14 = 19.0_dp
      t15 = 4.0_dp
      t16 = 20.0_dp * s25
      t17 = t14 * s16
      t18 = 7.0_dp * t8
      t19 = t13 * s26
      t20 = 21.0_dp
      t21 = 6.0_dp
      t22 = s16 * t21
      t23 = (s25 * t20 + t22) * s25
      t24 = 8.0_dp
      t25 = 17.0_dp * s25
      t26 = t24 * s16
      t27 = s26 ** 2.0_dp
      t28 = t27 ** 2.0_dp
      t29 = s26 * t27
      t30 = t8 * s25
      t31 = s25 ** 2.0_dp
      t32 = t31 ** 2.0_dp
      t33 = s25 * t32
      t34 = s25 * t31
      t35 = s16 * t31
      t36 = 30.0_dp
      t37 = s16 + s26
      t38 = s26 * t37
      t39 = 7.0_dp * s16
      t40 = t1 * s26
      t41 = 5.0_dp * s25
      t42 = s16 * t13
      t43 = 22.0_dp * t31
      t44 = 11.0_dp
      t45 = t44 * s16
      t46 = t15 * t8
      t47 = t4 * s26
      t48 = t15 * s16
      t49 = t24 * s26
      t50 = t3 ** 2.0_dp
      t51 = t10 * s25
      t52 = t50 * s26
      t53 = s34 ** 2.0_dp
      t54 = t53 ** 2.0_dp
      t55 = s34 * t53
      t56 = 12.0_dp * s25
      t57 = t21 * s26
      t58 = 9.0_dp * s25
      t59 = t52 * t21
      t60 = s16 * s25
      t61 = s34 * t13
      t62 = -t1 * (t60 + t61) + t53 + t31 + t8
      t63 = t10 * t31
      t64 = t50 * t27
      t65 = s56 ** 2.0_dp
      t66 = t50 * t29
      t67 = 13.0_dp * s16
      t68 = t31 * ((t67 - t41) * s25 - t44 * t8) + t59 * (-t2 + s16)
      t69 = (-s34 + t13 + s26) * s26 + (s16 + s25 - s34 + t40 + s56) * s
     #56 + t60
      t70 = 5.0_dp * s16
      t71 = t2 + s16
      t72 = 7.0_dp * s25
      t73 = s16 + t72
      t74 = t31 + t8
      t75 = t4 * t74
      t76 = t1 * s16
      t77 = t15 * s25
      t78 = t11 + t41
      t79 = t15 * t35 * t13
      t80 = t78 * t29
      t81 = t44 * s25
      t82 = t11 + t81
      t83 = s25 * (t67 + t12)
      t84 = t70 + t72
      t85 = s25 * t84 + t46
      t86 = 36.0_dp * s16
      t87 = 9.0_dp * t8
      t88 = (51.0_dp * s25 + t86) * s25 + t87
      t89 = s16 * t34
      t90 = 10.0_dp * t89
      t91 = t21 * s25
      t92 = t8 * t31
      t93 = t15 * s26
      t94 = s25 + t47
      t95 = s16 * t32
      t96 = 33.0_dp * s25
      t97 = 9.0_dp * s16
      t98 = t89 - t64
      t99 = t15 * t10
      t100 = t4 * t8
      t101 = 5.0_dp * t8
      t102 = t48 * s25
      t103 = t1 * t74
      t78 = ((((s34 - t11 - t41) * s34 - t102 + t103) * s56 - s25 * ((t6
     #7 - t72) * s25 - t101) - ((-t11 - t93 - t81 + s34) * s34 + 12.0_dp
     # * s16 * t7 + (20.0_dp * s26 + t25) * s25 + t100) * s34 + t10 + t5
     #2 * t24) * s56 - (-t24 * t34 - t99) * s25 + t4 * (-(s25 * t85 + s2
     #6 * (t1 * t78 * s26 + t88 / 3.0_dp)) * s34 + t52 * t73) - 12.0_dp 
     #* t98 + t53 * (-(t47 + t41) * s34 + (t96 + t57 + t97) * s26 + t83)
     #) * s56 + t24 * t52 * (s25 * t71 + t27) + t4 * (t64 * t73 + t33) -
     # (t91 * t85 * s26 + (-((t4 * t82 + t93) * s26 + t83 * t1) * s26 - 
     #t31 * (t58 + t45) + (s26 + t5) * t94 * s34) * s34 + t80 * t15 + t2
     #7 * t88 + 9.0_dp * t32 + t90 + 13.0_dp * t92) * s34 - t95 - t92 * 
     #(-t70 + t72)
      t85 = s34 * t7
      t88 = t8 * t34
      t93 = t89 * t74
      t104 = t50 * t28
      t105 = s26 + s56
      t106 = s16 - t77
      t107 = -s16 * t3
      t108 = t24 * s25
      t109 = (t22 + t12) * s25
      t110 = t71 * t27
      t111 = s25 * t7
      t112 = t60 * t4
      t113 = t2 + s26
      t114 = 10.0_dp * t31
      t115 = 16.0_dp * s25
      t116 = 12.0_dp * s16
      t117 = t19 * t21
      t118 = t21 * t8
      t119 = t113 * s34
      t120 = 10.0_dp * s25
      t121 = s16 + t120
      t122 = t84 * s26
      t123 = 25.0_dp * s25
      t124 = 17.0_dp * s16
      t125 = t123 + t124
      t126 = t115 + t70
      t67 = t31 * (t67 + t56)
      t127 = t15 * t31
      t128 = 29.0_dp * s25
      t129 = (t124 + t128) * s25 + 14.0_dp * t8
      t130 = (t48 + t91) * s25 + t8
      t131 = s25 * (t48 + t81)
      t132 = 48.0_dp
      t133 = 5.0_dp * s26
      t134 = t8 * t24
      t135 = t1 * t53
      t84 = ((-((-s34 * t1 + t70 + t72) * s34 - t75 + t91 * s16) * s56 -
     # s25 * ((-t120 + t17) * s25 - t134) - (t15 * (s16 * (s16 + t133) +
     # (7.0_dp * s26 + t48 + t91) * s25) + t135) * s34 + t53 * (t115 + t
     #70 + t49) + t10 + 12.0_dp * t52) * s56 - (-t34 * t44 - t99) * s25 
     #+ t4 * (t52 * t121 + t92) - 18.0_dp * t98 - (((t72 + t57) * s34 - 
     #s25 * (t123 + t124) - (s25 * t132 + 15.0_dp * s16 + 12.0_dp * s26)
     # * s26) * s34 + s25 * t129 + (t122 * t21 + 12.0_dp * t130) * s26) 
     #* s34) * s56 - t1 * (t113 * t94 * t55 - t131 * t52) - t31 * ((-t12
     #7 + t118) * s25 - 5.0_dp * t10) - t4 * (-t64 * t121 + t95) - s34 *
     # (t2 * t129 * s26 + t84 * t15 * t29 - (((t4 * t126 + t49) * s26 + 
     #t2 * t125) * s26 + t67) * s34 + 12.0_dp * t130 * t27 + 12.0_dp * t
     #32 + t90 + 14.0_dp * t92) + 12.0_dp * t66
      t90 = (t116 + s25) * s25 + t18
      t94 = 14.0_dp * t60
      t98 = t94 + t75
      t113 = 20.0_dp * s16
      t120 = (t113 + t58) * s25 + t8
      t94 = 5.0_dp * t74 + t94
      t123 = 15.0_dp * s25
      t124 = t13 * t1
      t129 = t2 + t70
      t130 = t5 + t48
      t136 = t71 * s26
      t137 = t21 * t31
      t44 = t89 * t44
      t138 = 5.0_dp * t94
      t139 = 42.0_dp * t8
      t140 = (91.0_dp * s16 + t96) * s25 + t139
      t86 = (t25 + t86) * s25
      t141 = t101 + t86
      t142 = t5 * t140
      t143 = t141 * t15
      t144 = t15 * t120
      t128 = ((95.0_dp * s16 + t128) * s25 + 83.0_dp * t8) * s25 + 33.0_
     #dp * t10
      t145 = t8 * t20
      t146 = ((26.0_dp * s16 + t81) * s25 + t145) * s25 + t1 * t10
      t147 = 10.0_dp * t94
      t148 = t3 * s25
      t9 = (-((-t124 * s34 + t94) * s34 + t4 * ((t8 + t148) * s25 - t10)
     #) * s56 - s25 * ((t36 * t8 - 9.0_dp * t31) * s25 - 18.0_dp * t10) 
     #+ t1 * (-t55 * t71 + t89) + ((10.0_dp * t19 + t101 + t86) * s34 - 
     #((52.0_dp * s16 + 22.0_dp * s25) * s25 + t139) * s25 - t99 - t138 
     #* s26) * s34 + t9 + 15.0_dp * t19 * t50) * s56 - (-t1 * t54 + 66.0
     #_dp * t92) * s25 + 9.0_dp * s25 * (t32 + t9) + 24.0_dp * t35 * t74
     # - ((((16.0_dp * s26 + t113 + t123) * s25 + t26 * s26) * s34 - s25
     # * t140 - (t143 + 20.0_dp * t19) * s26) * s34 + s25 * t128 + (t147
     # * s26 + t146 * t24) * s26) * s34 + t52 * (t19 * t36 + t144)
      t9 = t9 * s56 - s25 * (-9.0_dp * t52 * t98 + t30 * (t132 * t31 - t
     #145)) - t21 * (t10 * t34 - t64 * t120) + t36 * (s16 * t33 + t66 * 
     #t13) + t4 * (t31 * t32 + t55 * (-(t133 * t130 + 10.0_dp * t60) * s
     #25 - t15 * (t34 + t110) + s25 * (s25 + t40) * s34)) - s34 * (t5 * 
     #t128 * s26 - (t141 * t21 * t27 + t142 * s26 + 20.0_dp * t13 * t29 
     #+ 18.0_dp * t32 + 90.0_dp * t89 + 72.0_dp * t92) * s34 + 12.0_dp *
     # t146 * t27 + t147 * t29 + 12.0_dp * t33 + 66.0_dp * t63 + 72.0_dp
     # * t88 + 90.0_dp * t95)
      t22 = t7 ** 2.0_dp * (12.0_dp * t92 * t13 + s26 * (t27 * t94 + t60
     # * ((14.0_dp * s16 + t96) * s25 + 25.0_dp * t8)) + t15 * (((t22 + 
     #t5) * s25 + t134) * s25 + t10) * t27) * s34
      t86 = t7 * (s26 * ((s25 * ((57.0_dp * s16 + t12) * s25 + 37.0_dp *
     # t8) + ((34.0_dp * s16 + t123) * s25 + t101) * s26) * s26 + t35 * 
     #(35.0_dp * s16 + t96)) + t124 * t28 + t118 * t34) * t53
      t92 = t64 * ((s25 * t98 + t13 * t27) * s26 + t31 * t90)
      t9 = -t21 * t8 * t32 * t74 - t4 * t92 - (-t55 * (((t136 * t1 + t41
     # * t130) * s26 + t137 * t129) * s26 + t44) + t89 * t50 * (t81 + t1
     #7)) * s26 - (t9 * s56 + t21 * s25 * s26 * (s25 * t50 * t90 + t54 *
     # t7) - (((((t123 * t130 + t136 * t24) * s26 + 12.0_dp * t31 * t129
     #) * s26 + t44) * s34 - 36.0_dp * t19 * t31 * (s25 + t48) - t142 * 
     #t27 - 10.0_dp * t13 * t28 - t143 * t29 - 41.0_dp * t88 - 33.0_dp *
     # t95) * s34 + t7 * (s26 * (s25 * (((147.0_dp * s16 + 24.0_dp * s25
     #) * s25 + 106.0_dp * t8) * s25 + 83.0_dp * t10) + s26 * (((138.0_d
     #p * s16 + 63.0_dp * s25) * s25 + 143.0_dp * t8) * s25 + 16.0_dp * 
     #t10)) + t138 * t29 + t35 * ((38.0_dp * s16 + t96) * s25 + 49.0_dp 
     #* t8))) * s34 + t58 * t64 * t98 + t144 * t66 + 15.0_dp * t104 * t1
     #3 - t89 * (((t11 - t81) * s25 + 27.0_dp * t8) * s25 - t10 * t14)) 
     #* s56 - s25 * t27 * (t40 + t5) * t54 + t22 - t104 * t120 - t86 + 1
     #2.0_dp * t10 * t33
      t22 = t21 * t111
      t29 = t74 * s25
      t33 = t74 * t31
      t44 = s25 * (t35 - t52)
      t50 = s25 + s26 - s34 + s56
      t54 = t3 * s26
      t74 = 0.1e1_dp / s16
      t81 = 0.1e1_dp / t62
      t86 = 0.1e1_dp / t105 ** 2.0_dp
      t69 = 0.1e1_dp / t69 ** 2.0_dp
      t9 = (-t1 * t93 - t84 * s56 - t53 * (t1 * (t28 + t89) + s26 * ((s2
     #5 * t125 + s26 * t126) * s26 + t67)) + t85 * (t27 * ((t25 + t45) *
     # s25 + t46 + t122) + t79 + t2 * ((t11 + t91) * s25 + t101) * s26) 
     #+ s26 * ((t40 + t72) * s26 + t127) * t55 - t52 * ((s26 * t121 + t1
     #31) * s26 + t31 * (t77 + t70)) - t104 * t4 + t46 * t32) * MBL1001(
     #0) - t1 * (t93 + t104) - t15 * s25 * (t64 * t71 - t88) - t78 * s56
     # - t53 * (s26 * (((s26 + t82) * s26 + t83) * s26 + t31 * (t58 + t4
     #5)) + t76 * t34) + t85 * (t4 * s26 * (s25 * (t76 * s25 + t75) + s2
     #6 * ((t11 + t77) * s25 + t8)) + t79 + t80) + s26 * ((s26 + t41) * 
     #s26 + t31 * t4) * t55 - t52 * (t27 * t73 + t31 * (t5 + t70)) + t9 
     #* MBL1011(1) * t74
      t5 = s34 * t50 * (t1 * (t33 + t64) - t15 * t44 - ((t119 - (s26 + t
     #76) * s26 - t112 - t22) * s34 + t7 * (t19 * t4 - t107 + t137)) * s
     #34 - (-(-t4 * t61 - t102 + t103 + t53) * s56 - t15 * (t52 + t29) +
     # s34 * (t53 + (t76 + t58) * s25 + t117 + t8) - t135 * (t37 + t5) +
     # t35 * t24) * s56) * MBL1110(1)
      t5 = t5 + s25 * t68 - t1 * (t62 * s56 * t65 + t66) + t4 * (t64 * t
     #6 + t63) + ((((-t40 - t41) * s34 + (18.0_dp * s26 + t16 + t39) * s
     #25 + t15 * t38) * s34 - t1 * ((t27 + t8 + t23) * s26 + t30) - t36 
     #* t34 - (t25 + t26) * t27 - t35) * s34 + t7 * (s25 * ((t16 - t17) 
     #* s25 + t18) + ((-t11 + t12) * s25 + t8 + t19 * t15) * s26)) * s34
     # + ((((s16 * t20 - t58) * s25 - 15.0_dp * t8) * s25 + t10 * t4 + t
     #15 * t55 + s34 * (-(t25 + t57 + t26) * s34 + 12.0_dp * t19 + t42 +
     # t43) - t59) * s56 - t1 * ((t53 + (t25 + t47 + t26) * s26 + t23 + 
     #t8) * s34 - s25 * ((s25 * t14 - t45) * s25 + t46) - (t19 * t21 + t
     #42 + t43) * s26) * s34 + t21 * (t51 + t52 * (-s26 + t6)) + t31 * (
     #(s16 * t36 - t56) * s25 - 24.0_dp * t8) + (t48 + t49 + t12) * t55)
     # * s56
      t4 = t69 * t81 * ((t34 * ((-t2 + t70) * s25 - t46) + (((t1 * t38 +
     # t111 * t24 + t112 - t119) * s34 - 12.0_dp * t34 - (t27 + t8 + t10
     #9) * s26 - t60 * t13 - t110 * t15) * s34 + t7 * (s25 * ((-t39 + t1
     #08) * s25 + t100) + (t19 * t1 + t24 * t31 - t107) * s26)) * s34 + 
     #((-s56 * t62 + s25 * ((-t77 + t97) * s25 - t118) + t1 * t55 + (-(t
     #47 + t48 + t108) * s34 + t114 + t117 + t42) * s34 + t10 - t52 * t4
     #) * s56 + t1 * ((t53 * (s16 + t40 + t77) + s25 * ((-t48 + t108) * 
     #s25 + t1 * t8) + (t19 * t4 + t114 + t42) * s26) * s34 + t51 + t52 
     #* t106) + t31 * ((t116 - t41) * s25 - t87) - t53 * (t53 + (t115 + 
     #t47 + t26) * s26 + t109 + t8) - t64 * t4) * s56 + t63 - t66 + t52 
     #* (t76 - t41) * s25 + t64 * t106) * I300s16(1) + t74 * t5 + t9 * t
     #86)
      result = t4 / 4.0_dp + t74 * (t69 * (-t50 * (t1 * ((t54 + t85) * s
     #25 + ((s16 + s34 - s25) * s25 + t54) * s56) + (-t53 + t148) * s25 
     #+ t3 * (t65 + t27)) * MBL1111D2(0) + s25 * s34 * (-t1 * (t55 * (t2
     # + t37) + t44) + ((t53 + (t2 + t47) * s16 + t27 + t8 + t22) * s34 
     #- t7 * ((-t76 + t77) * s25 + t136 + t8)) * s34 + ((-t1 * (s16 + s3
     #4) * s25 - (s16 - s34) * s34 + t31 + t8) * s56 + t1 * (-t55 + t52 
     #+ t29) - s34 * (-(t40 + t11 + t91) * s34 + (t76 + t77) * s26 - t10
     #7 + t137) - t35 * t15) * s56 + t64 + t33) * MBL1101(0) * t81) - s2
     #5 * (s25 + s26 + s56) * MBL1011(0) * t86) / 2.0_dp

           intHLs16s25s26s34s56x1123D6eps0 = result
       end function intHLs16s25s26s34s56x1123D6eps0

       function intHLs16s25s26s34s56x1123D6eps1()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1123D6eps1
           complex(dp) ::  t1,t2

           complex(dp) :: result

      t1 = s26 + s56
      t2 = 0.1e1_dp / t1 ** 2.0_dp
      result = 0.1e1_dp / s16 * t2 * (-s25 / 2.0_dp - t1 / 4.0_dp - s25 
     #* (s25 + s26 + s56) * MBL1011(1) / 2.0_dp)

           intHLs16s25s26s34s56x1123D6eps1 = result
       end function intHLs16s25s26s34s56x1123D6eps1

       function intHLs16s25s26s34s56x1131D4eps0()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1131D4eps0
           complex(dp) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           complex(dp) ::  t20,t21,t22,t23,t24,t25,t3,t4,t5,t6,t7,t8
           complex(dp) ::  t9

           complex(dp) :: result

      t1 = s16 + s25 - s34
      t2 = s26 + s56
      t3 = 2.0_dp
      t4 = t2 * t3 + t1
      t5 = s16 + s25
      t6 = s16 * s25
      t1 = (-s34 + t5 + s26) * s26 + (s26 * t3 + s56 + t1) * s56 + t6
      t7 = s25 + s26 + s56
      t8 = s56 ** 2.0_dp
      t9 = s26 ** 2.0_dp
      t10 = (s25 - s34) * s25
      t11 = s25 + s26
      t12 = s25 * s26
      t13 = s25 + s26 - s34 + s56
      t14 = s25 ** 2.0_dp
      t15 = s34 ** 2.0_dp
      t16 = s26 * s34
      t17 = t3 * s16
      t18 = s16 - s25
      t19 = t18 * s25
      t20 = t18 * s26
      t21 = s16 ** 2.0_dp
      t22 = t21 + t14
      t18 = t18 ** 2.0_dp
      t23 = t18 * s26
      t24 = 0.1e1_dp / t1
      t7 = 0.1e1_dp / t7
      t1 = 0.1e1_dp / t1
      t25 = 0.1e1_dp / s16
      t13 = 0.1e1_dp / t13
      t2 = 0.1e1_dp / t2
      t18 = (t3 * ((s16 * t14 - t23) * s25 - (-t17 * t14 + s25 * t22 + (
     #-s25 * t5 - t17 * s26 + t16) * s34 + t23) * s56 + t11 * (s16 * s26
     # + t14) * s34) - t14 * t22 - t15 * (t9 + t14) - t8 * (-t3 * s16 * 
     #(s25 + s34) + t14 + t15 + t21) - t9 * t18) * MBL1111D2(0) * t25 * 
     #t7 * t1 + I300s16(0)
      t21 = 0.1e1_dp / 0.2e1_dp
      result = t21 * (t13 * (t15 * MBL1110(0) * t7 * t25 + t1 * ((t3 * (
     #s56 * t11 + t12) + t10 + t8 + t9) * MBL1001(0) + (t17 * s25 * t14 
     #+ s16 * s26 * t9 + s16 * s56 * t8 + ((t16 - t3 * (s26 * t5 + t6) -
     # t9) * s34 + t9 * (3.0_dp * s16 + s25)) * s25 + (s16 * ((6.0_dp * 
     #s25 + 3.0_dp * s26) * s26 + 4.0_dp * t14) - t3 * s25 * ((s16 + s25
     # + s26) * s34 - t12) + s25 * (t15 + t14)) * s56 + t8 * (3.0_dp * s
     #16 * t11 + t10) + t14 * (4.0_dp * s16 + s25) * s26) * MBL1011(1) *
     # t25) * t2 - s34 * ((-s25 * t3 - s16 - s26 + s34) * s34 - (s16 - s
     #25 + s34) * s56 - t19 - t20) * MBL1110(1) * t1 * t25) + t7 * t18 +
     # t14 * MBL1011(0) * t25 * t2 * t7 - s25 * s34 * (-(s25 - s26) * s3
     #4 - (s16 - s25 - s34) * s56 - t19 - t20) * MBL1101(0) * t25 * t7 *
     #* 2.0_dp * t1 + t4 * I300s16(1) * t24) + t4 * t25 * t1

           intHLs16s25s26s34s56x1131D4eps0 = result
       end function intHLs16s25s26s34s56x1131D4eps0

       function intHLs16s25s26s34s56x1131D4eps1()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1131D4eps1
           complex(dp) ::  t1,t2,t3,t4,t5,t6

           complex(dp) :: result

      t1 = s26 + s56
      t2 = 2.0_dp
      t3 = s25 + s26 - s34 + s56
      t4 = s25 + s26 + s56
      t4 = 0.1e1_dp / t4
      t3 = 0.1e1_dp / t3
      t5 = 0.1e1_dp / t1
      t6 = 0.1e1_dp / 0.2e1_dp
      result = t6 * (I300s16(1) * t4 + (t5 * ((t1 * t2 + s25 - s34) * t3
     # + s25 ** 2.0_dp * MBL1011(1) * t4) + s34 ** 2.0_dp * MBL1110(1) *
     # t4 * t3) / s16)

           intHLs16s25s26s34s56x1131D4eps1 = result
       end function intHLs16s25s26s34s56x1131D4eps1

       function intHLs16s25s26s34s56x1132D6eps0()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1132D6eps0
           complex(dp) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           complex(dp) ::  t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30
           complex(dp) ::  t31,t32,t33,t34,t35,t36,t37,t38,t39,t4,t40,t41
           complex(dp) ::  t42,t43,t44,t45,t46,t47,t48,t49,t5,t50,t51,t52
           complex(dp) ::  t53,t54,t55,t56,t57,t58,t59,t6,t60,t61,t62,t63
           complex(dp) ::  t64,t65,t7,t8,t9

           complex(dp) :: result

      t1 = s25 + s26 + s56
      t2 = 3.0_dp
      t3 = 2.0_dp
      t4 = t3 * s16
      t5 = t2 * s25
      t6 = s26 ** 2.0_dp
      t7 = t6 ** 2.0_dp
      t8 = s26 * t6
      t9 = (-s25 + s34) * s25
      t10 = 4.0_dp
      t11 = t2 * s26
      t12 = t10 * s25
      t13 = t10 * s26
      t14 = t4 * s25
      t15 = s26 + s56
      t16 = t3 * s26
      t17 = s16 * s25
      t18 = (-s34 + s16 + s25 + s26) * s26 + (t16 + s25 + s16 - s34 + s5
     #6) * s56 + t17
      t19 = 8.0_dp * s16
      t20 = 5.0_dp
      t21 = 10.0_dp * s25
      t22 = s25 + s26
      t23 = t20 * s25
      t24 = s16 * t22
      t25 = -s16 + s34
      t26 = s25 ** 2.0_dp
      t27 = t26 ** 2.0_dp
      t28 = s25 * t26
      t29 = 10.0_dp * s26
      t30 = s56 ** 2.0_dp
      t31 = t30 ** 2.0_dp
      t32 = s56 * t30
      t33 = t3 * s25
      t34 = s26 + s16
      t35 = (-s26 + t33) * s34
      t36 = t4 * t26
      t37 = 9.0_dp * s16
      t38 = s25 + t37
      t39 = s16 + t5
      t40 = t10 * s16
      t41 = s25 + t40
      t42 = s16 ** 2.0_dp
      t43 = (t5 + t19) * s25
      t44 = t43 - t42
      t45 = 6.0_dp * s25
      t46 = s16 * s26
      t47 = s16 * t26
      t48 = t20 * s16
      t49 = t40 * s26
      t50 = s34 ** 2.0_dp
      t51 = 18.0_dp * t26 * t41
      t52 = -t13 * t44
      t53 = t3 * t26
      t54 = s16 * t6
      t55 = s16 * t2 + t45
      t56 = s16 * t25
      t57 = s34 * t22
      t58 = s25 - s16
      t59 = t58 * s26
      t60 = t58 * s25
      t61 = t58 ** 2.0_dp
      t62 = t61 * s26
      t63 = t61 * t6
      t64 = (t42 + t26) * s25
      t65 = 0.1e1_dp / s16
      t18 = 0.1e1_dp / t18 ** 2.0_dp
      t15 = 0.1e1_dp / t15 ** 2.0_dp
      t7 = (6.0_dp * t26 * (t26 * t42 + t41 * t8) + s16 * s26 * t7 + s16
     # * s56 * t31 + (s16 * t28 * (14.0_dp * s16 + 11.0_dp * s25) + t57 
     #* ((-s25 * (t40 + t45) + t46) * s26 - 11.0_dp * t47)) * s26 + ((t2
     #0 * t7 + 11.0_dp * t27) * s16 + s26 * ((-t52 + t51) * s26 + t53 * 
     #t38 * t39) - (((9.0_dp * s25 * (s16 + t33) - t49) * s26 + 6.0_dp *
     # t26 * (t48 + t33)) * s26 + 11.0_dp * s16 * t28) * s34 + 6.0_dp * 
     #t26 * s26 * t50 + 14.0_dp * t42 * t28) * s56 + t30 * (((t29 * s16 
     #+ 6.0_dp * t44) * s26 + t51) * s26 + t2 * ((s34 * t26 - s25 * (t23
     # * s16 + s26 * t55) + t3 * (t54 - t28)) * s34 + t27) + t47 * (28.0
     #_dp * s25 + t37)) + t31 * (t48 * s26 + t43 + t56) + t32 * ((-s25 *
     # t55 + t49) * s34 + 6.0_dp * t28 - t52 + 24.0_dp * t47 + 10.0_dp *
     # t54) + t44 * t7 + t26 * t38 * t39 * t6 + t2 * t26 * t6 * t50) * M
     #BL1011(1) * t15
      t9 = t15 * ((t36 + (s25 * (s16 + t12) + (-s16 + t23 + s26) * s26) 
     #* s26 + s26 * (s26 - t12) * s34 + ((t23 + t11 + t25 + s56) * s56 +
     # t3 * s26 * t25 + (t21 + t11) * s26 - t10 * t9 + t17) * s56) * MBL
     #1001(0) + s25 * ((-s34 * t2 + t13) * s26 + t14) + (s25 * (t5 + t4)
     # + t6) * s26 + ((t12 + t11 + s56) * s56 + (8.0_dp * s26 + t4) * s2
     #5 - t2 * (t9 - t6)) * s56) * t1
      t7 = t18 * ((-t2 * s16 * (t26 + t6) + t28 * t3 - (t6 - s25 * (-6.0
     #_dp * s16 + t5)) * s26 - (-t22 * (s16 - t12 + t16) - t35) * s34 - 
     #((-s34 * t3 + t2 * t34 + s56) * s56 - t2 * (t26 - t6) + 6.0_dp * t
     #24 + s34 * (s34 - s16 - t13 + t33)) * s56) * I300s16(1) + (t7 - t2
     # * (t32 + t8) - t26 * (-t23 + t19) - (s25 * (16.0_dp * s16 - 7.0_d
     #p * s25) + (s25 + t19) * s26) * s26 - ((-t23 + t16) * s34 + (-s26 
     #* t20 + t21 - t4) * t22) * s34 - ((-s34 * t20 + s25 + 9.0_dp * s26
     # + t19) * s56 + t3 * (s25 * s26 + s34 * t25) + 9.0_dp * t6 + 16.0_
     #dp * t24 - 7.0_dp * t26 + t23 * s34 - t29 * s34) * s56 - s34 * (-t
     #3 * (-t26 * t58 - t58 * t6) + (t22 * (-t12 + t34) + t35) * s34 + (
     #(t3 * t58 + s34) * s56 - t10 * (-t59 - t60) - s34 * (s34 - s16 + t
     #5 - t16)) * s56 + t60 * t13) * MBL1110(1)) * t65 + t9)
      t1 = 0.1e1_dp / t1
      result = t7 / 4.0_dp + t65 * (((-t2 * (s25 * (s25 * t22 * t50 + t6
     #2 * s25 + t63) + ((-t36 - (t46 + t26) * s34 + t62 + t64) * s56 - t
     #3 * s25 * (t47 - t62) + t26 * (t50 + t42 + t26) + t63 + t57 * (-t5
     #3 + s16 * (s25 - s26))) * s56) - t26 * (-s34 * t50 + t64) - t32 * 
     #(-t14 - t56 + t26) - t61 * t8 + t22 ** 2.0_dp * ((t5 - t4) * s25 +
     # t46) * s34 + t4 * t27) * MBL1111D2(0) + s25 * s34 * (-t3 * ((-t57
     # + t59) * s25 - ((s16 + s34 - s25) * s25 - t59) * s56) - (t60 + t5
     #0) * s25 - t58 * (t30 + t6)) * MBL1101(0)) * t1 * t18 + t26 * MBL1
     #011(0) * t15) / 2.0_dp

           intHLs16s25s26s34s56x1132D6eps0 = result
       end function intHLs16s25s26s34s56x1132D6eps0

       function intHLs16s25s26s34s56x1132D6eps1()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1132D6eps1
           complex(dp) ::  t1,t2,t3,t4

           complex(dp) :: result

      t1 = s26 + s56
      t2 = 2.0_dp
      t3 = 0.1e1_dp / 0.4e1_dp
      t4 = 0.1e1_dp / t1 ** 2.0_dp
      result = 0.1e1_dp / s16 * t4 * (t3 * (s25 * t2 - t1) + s25 ** 2.0_
     #dp * MBL1011(1) / 2.0_dp)

           intHLs16s25s26s34s56x1132D6eps1 = result
       end function intHLs16s25s26s34s56x1132D6eps1

       function intHLs16s25s26s34s56x1141D6eps0()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1141D6eps0
           complex(dp) ::  t1,t10,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109
           complex(dp) ::  t11,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t12
           complex(dp) ::  t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t13,t130
           complex(dp) ::  t131,t132,t133,t134,t14,t15,t16,t17,t18,t19,t2,t20
           complex(dp) ::  t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30,t31
           complex(dp) ::  t32,t33,t34,t35,t36,t37,t38,t39,t4,t40,t41,t42
           complex(dp) ::  t43,t44,t45,t46,t47,t48,t49,t5,t50,t51,t52,t53
           complex(dp) ::  t54,t55,t56,t57,t58,t59,t6,t60,t61,t62,t63,t64
           complex(dp) ::  t65,t66,t67,t68,t69,t7,t70,t71,t72,t73,t74,t75
           complex(dp) ::  t76,t77,t78,t79,t8,t80,t81,t82,t83,t84,t85,t86
           complex(dp) ::  t87,t88,t89,t9,t90,t91,t92,t93,t94,t95,t96,t97
           complex(dp) ::  t98,t99

           complex(dp) :: result

      t1 = 6.0_dp
      t2 = 2.0_dp
      t3 = 38.0_dp
      t4 = s16 ** 2.0_dp
      t5 = t2 * s25
      t6 = t1 * t4
      t7 = 17.0_dp * s25
      t8 = 12.0_dp * s16
      t9 = 4.0_dp * s25
      t10 = 5.0_dp
      t11 = t1 * s26
      t12 = t10 * s25
      t13 = s34 ** 2.0_dp
      t14 = t13 ** 2.0_dp
      t15 = s34 * t13
      t16 = s16 - s34
      t17 = s26 * t16
      t18 = 19.0_dp
      t19 = 3.0_dp
      t20 = s16 * t1
      t21 = s26 ** 2.0_dp
      t22 = t21 ** 2.0_dp
      t23 = s26 * t21
      t24 = s26 * t22
      t25 = s16 * t18
      t26 = (-t20 - t5) * s34
      t27 = t19 * (t13 + t4)
      t28 = s56 ** 2.0_dp
      t29 = t28 ** 2.0_dp
      t30 = s56 * t28
      t31 = s56 * t29
      t32 = t23 + t30
      t33 = 13.0_dp * s16
      t34 = s16 + s25
      t35 = t2 * s26
      t36 = s16 * s25
      t37 = (-s34 + t34 + s26) * s26 + (s16 + s25 - s34 + t35 + s56) * s
     #56 + t36
      t38 = 7.0_dp
      t39 = t38 * s16
      t40 = t19 * s25
      t41 = t40 + t39
      t42 = t2 * s16
      t43 = t40 + t42
      t44 = 9.0_dp
      t45 = t44 * s25
      t46 = s16 + t45
      t47 = t1 * s25
      t48 = s16 - t47
      t49 = s16 - t40
      t50 = s25 + s26
      t51 = t38 * s25
      t52 = t19 * s16
      t53 = s25 ** 2.0_dp
      t54 = t53 ** 2.0_dp
      t55 = s25 * t54
      t56 = s25 * t53
      t57 = 4.0_dp * s16
      t58 = t57 * t53
      t59 = t19 * s26
      t60 = s25 * t46
      t61 = t42 * t53
      t62 = 15.0_dp * s26
      t63 = 10.0_dp * s16
      t64 = s25 * s34
      t65 = t21 * t48
      t66 = s25 * s26
      t67 = t49 * t10
      t68 = 30.0_dp
      t69 = t68 * s26
      t70 = 12.0_dp * s25
      t71 = 25.0_dp * s25
      t72 = t10 * s16
      t73 = 18.0_dp * t21
      t74 = 8.0_dp * s16
      t75 = -s16 + s26
      t76 = s16 * t53
      t77 = t67 * s26
      t78 = t44 * t56
      t79 = t19 * s34
      t80 = s25 + s34
      t81 = t50 ** 2.0_dp
      t82 = t66 * (-t15 * (s26 + t40) + t41 * t56)
      t18 = ((t28 * (-(-t19 * t80 + s16) * s56 - (-t70 + t42) * s25 - (s
     #25 - t52 - t62 + t79) * s34 - t77) - (-t53 * t46 * t1 + 10.0_dp * 
     #t49 * t21) * s26 - 12.0_dp * s25 * (t65 - t56) - ((-(s16 + 14.0_dp
     # * s25 + t11) * s25 + t73 + t64) * s34 + ((-18.0_dp * s16 + t47 - 
     #t69) * s26 - t70 * (s16 - t12)) * s26 + t53 * (t72 + t71)) * s34 +
     # t74 * t56) * s56 - t53 * (-t21 * t46 * t1 - t41 * t53) + 8.0_dp *
     # t66 * (t43 * t53 - t65) - ((-((28.0_dp * s26 + t45 + t52) * s25 +
     # s26 * (t42 + t11)) * s25 + 12.0_dp * t23 + t64 * (t40 + t35)) * s
     #34 + t50 * (((s25 * t18 - t62 - t8) * s26 + 41.0_dp * t53) * s26 +
     # t53 * (t45 + t63))) * s34 - t67 * t22) * s56 - t24 * t49 - t81 * 
     #(-t19 * t23 + s26 * (s25 * (t42 + t45) + s26 * (-t52 + t51)) + t58
     #) * s34 + t50 * (((-t59 + t12) * s26 + t60) * s26 + t61) * t13 + t
     #82
      t41 = s16 * t54
      t46 = s26 + s56
      t49 = s25 + s26 - s34 + s56
      t65 = s25 + s26 + s56
      t67 = 16.0_dp
      t82 = t19 * t4
      t83 = (s16 * t67 - s25) * s25
      t84 = 11.0_dp * s16
      t85 = t59 - t5
      t86 = t9 + t72
      t87 = s16 - t5
      t88 = s16 - s25
      t89 = 23.0_dp
      t90 = t21 - t53
      t91 = t53 * s26
      t92 = 15.0_dp * s16
      t93 = s25 * t88
      t94 = -t9 + t52
      t95 = t68 * t21
      t96 = s16 * t21
      t97 = t22 - t54
      t98 = s16 + s26
      t99 = s26 * t98
      t63 = ((t28 * (-(s34 * t10 - t51 - t69 - t72) * s34 - 10.0_dp * t9
     #3 - 15.0_dp * t99 - t28) - t68 * s26 * (t53 * t87 + t96) + ((-(s16
     # - 13.0_dp * s25 + 12.0_dp * s26) * s25 - t95 + t64) * s34 + ((60.
     #0_dp * s26 + t1 * (t51 + t72)) * s26 + t70 * t94) * s26 + t53 * (-
     #29.0_dp * s25 + t72)) * s34 - 15.0_dp * t97 - 60.0_dp * t93 * t21)
     # * s56 - t13 * (12.0_dp * s25 * t90 + (s25 * (-26.0_dp * s25 + t42
     #) + 20.0_dp * t21) * s26 - t76) - t2 * s34 * (t50 * (t53 * (t52 + 
     #t47) - 15.0_dp * t23 + s26 * ((s25 - t63) * s26 + s25 * (s25 * t89
     # - t74))) + s25 * (-s26 + t5) * t13) - t23 * ((t92 + t11) * s26 + 
     #40.0_dp * t93) + t54 * t86 - t68 * t91 * (s26 * t87 - t53)) * s56 
     #- t21 * t22 - t81 * (((-t10 * t88 - t11) * s26 - t5 * t48) * s26 +
     # t58) * s34 + t50 * ((-t10 * t21 + 12.0_dp * t53) * s26 + t66 * t7
     #5 + t61) * t13 + t66 * (-t15 * (-s26 + t9) + t56 * t86)
      t69 = t88 * s26
      t74 = s16 * t55
      t86 = 11.0_dp * s25
      t89 = 21.0_dp * s16
      t100 = 14.0_dp * s16
      t101 = (s25 + t100) * s25 + 8.0_dp * t4
      t102 = (t40 + t89) * s25 + t4
      t103 = -46.0_dp * s16
      t104 = (t103 - t45) * s25 + 15.0_dp * t4
      t105 = (s25 + t72) * s25 - 4.0_dp * t4
      t106 = t57 + t5
      t34 = s25 * t34
      t107 = t10 * t4
      t108 = -t5 - t72
      t109 = s16 + t5
      t110 = t40 - t72
      t111 = 4.0_dp * s26
      t112 = t84 * t53
      t113 = -t104
      t114 = 8.0_dp * s25
      t115 = 33.0_dp * s16
      t116 = t115 * t54
      t117 = (21.0_dp * s25 + t100) * s25 - t82
      t118 = -s25 + t8
      t119 = -36.0_dp * s26
      t120 = t67 * s26
      t121 = s25 - s26
      t122 = 18.0_dp * t4
      t123 = s25 * t105
      t124 = t53 * t113
      t125 = t78 * t102
      t71 = (82.0_dp * s16 + t71) * s25 - 11.0_dp * t4
      t126 = s25 + t42
      t127 = (-s25 + t84) * s25 + t107
      t128 = (t57 - t51) * t126
      t129 = 10.0_dp * t127
      t130 = -10.0_dp * t108
      t131 = 8.0_dp * s26
      t132 = 45.0_dp * t4
      t55 = ((((t89 * s26 + t132) * s26 - t123 * t68) * s26 - t124 * t1)
     # * s26 - t125) * s26 + t19 * s25 * (t15 * (-(s25 - t35) * s34 + s2
     #5 * (t9 + t42) + (-t131 + t40 - t72) * s26) - t55) + s34 * (((s26 
     #* (t11 * t118 - t117 * t19) + t53 * (-54.0_dp * s16 - 18.0_dp * s2
     #5)) * s25 + t130 * t23) * s34 + (t53 * t71 * t19 - t129 * t21) * s
     #26 - 90.0_dp * s16 * t97 - 12.0_dp * s25 * (t128 * t21 - t54) + t1
     #22 * t56) + t41 * (-24.0_dp * s16 - 42.0_dp * s25)
      t97 = s16 * s26
      t133 = -63.0_dp
      t134 = s16 * t23
      t45 = ((((s16 * t133 - t45) * s25 - t82) * s25 - t111 * t113) * s2
     #5 + t14 * t2 - t95 * t105) * s25 + ((-t100 * t53 + 4.0_dp * t66 * 
     #t118 + t130 * t21 + t40 * t4 - 21.0_dp * t56 + t64 * (t40 - t120 -
     # t72)) * s34 + (s25 * t71 - t131 * t128) * s25 + t21 * (-120.0_dp 
     #* t97 - t129)) * s34 + t134 * (60.0_dp * s16 + 35.0_dp * s26)
      t71 = t4 * s25
      t67 = -t44 * t54 - 15.0_dp * s25 * (s26 * t105 - t71) + s34 * (-t1
     #0 * t127 * s26 + ((-s25 + 10.0_dp * s26 + t8) * s25 + 25.0_dp * t9
     #7) * s34 - t71 * t67 + 14.0_dp * t56 + 20.0_dp * t76 - 90.0_dp * t
     #96) + t103 * t56 - t9 * t15 + t132 * t21 + 35.0_dp * t134
      t71 = -t19 * t56 + s34 * (-t108 * s34 + (t119 - t86 - t72) * s16 +
     # t53) + t70 * t4 + t122 * s26 - t92 * t53 + t89 * t21
      t95 = t1 * s34
      t100 = t50 * s34
      t91 = t100 * (-(t19 * t66 * (s26 * ((t39 + t47) * s25 - t4) + t36 
     #* (s16 + t86)) + t23 * (t108 * s26 + s25 * (t40 - t39)) + t6 * t56
     #) * s34 + t50 * (t1 * t21 * (s25 * (s25 * t106 - t4) - t96) + t23 
     #* (t34 - t107) + 12.0_dp * t4 * t56 + t52 * t91 * (t42 + t86)))
      t96 = t66 * (t15 * (((t110 - t111) * s26 + t47 * t109) * s26 + t11
     #2) - t41 * (t89 + t86))
      t22 = s16 * t23 * t22 + s16 * t30 * t29 + (-t1 * t66 * (t101 * t56
     # + t121 * t14) + t21 * ((((t39 * s26 + t122) * s26 - 15.0_dp * t12
     #3) * s26 - 4.0_dp * t124) * s26 - t125) + ((-t40 * t117 * t21 + t9
     # * t118 * t23 - t116 + t119 * t56 * (s25 + t52) - t78 * t4 - t108 
     #* t10 * t22 + t64 * (((t19 * t110 - t120) * s26 + t70 * t109) * s2
     #6 + t112)) * s34 + t50 * (t19 * t66 * (s25 * ((49.0_dp * s16 + t11
     #4) * s25 + t2 * t4) + s26 * ((t115 + t7) * s25 - 13.0_dp * t4)) + 
     #t23 * ((t12 - t25) * s25 - 25.0_dp * t4) + t116 + t68 * t4 * t56 -
     # 36.0_dp * s16 * t22)) * s34 - t74 * (t89 + t86)) * s56 + t28 * t5
     #5 + t29 * t67 + t30 * t45 + t31 * t71 - s25 * (t40 - t35) * t21 * 
     #t14 + s16 * (s26 * t38 + t52 - t95) * t28 * t29 + t53 * t104 * t22
     # + t91 + t96
      t29 = t97 + t53
      t38 = t53 + t4
      t39 = t88 ** 2.0_dp
      t40 = t38 * s25
      t45 = t39 * s26
      t55 = t2 * s16 * t80
      t67 = s25 * (t45 - t76)
      t68 = t39 * t21
      t71 = t38 * t53
      t80 = t50 * s26
      t86 = 0.1e1_dp / s16
      t46 = 0.1e1_dp / t46 ** 2.0_dp
      t49 = 0.1e1_dp / t49 ** 2.0_dp
      t37 = 0.1e1_dp / t37 ** 2.0_dp
      t14 = ((t19 * t14 * (-s34 + s25 + t111 + t42) - t13 * (((t35 + t52
     # + t51) * s16 + t44 * t90 + 12.0_dp * t66) * s34 + t50 * ((-t33 - 
     #t59) * s25 + 15.0_dp * t29 + 10.0_dp * t4)) + t95 * t39 * t81) * s
     #56 + t14 * (-t85 * s34 + (s16 - t114 + t59) * s25 + t1 * t99) + t2
     #8 * (t1 * s34 * (t39 * t50 + t15) - t13 * ((s26 * t44 + s16 + t47)
     # * s34 + s16 * (t62 + t72) + (s16 - t59 + t47) * s25)) - s34 * ((-
     #s25 + t79 + t72) * s34 - t2 * t38 + t57 * s25) * t30 - t81 * (t10 
     #* s16 * t98 + (-s16 * t44 - s26 + t114) * s25) * t13 + t50 * ((-t2
     #0 + t70) * s25 - t19 * (t80 + t4) - t97) * t15 + t2 * t39 * t50 * 
     #t81 * s34) * MBL1110(1) * t86
      t10 = t46 * ((t19 * (s16 * t24 + t31 * (t2 * (s26 - s34) + s16)) +
     # t2 * (t30 * ((t92 * s26 + 20.0_dp * t93) * s26 - t2 * s34 * (-(s2
     #6 * t10 + s25) * s34 + s25 * t94 + (t62 + t51 + t72) * s26) + 10.0
     #_dp * t23 - 10.0_dp * t56 + t72 * t53) - t74) - t63 + 10.0_dp * s2
     #5 * t23 * (s25 * t87 + t69) - 15.0_dp * t54 * t21) * MBL1001(0) - 
     #t2 * (s25 * (t23 * (-s26 * t48 + t60) + t41) + t30 * ((-t9 * t48 -
     # t77) * s26 + ((s25 - t11) * s34 - t2 * s25 * t75 + (t20 + t62) * 
     #s26 - 10.0_dp * t53) * s34 + t76 + t78)) - t18 - 4.0_dp * t56 * t4
     #3 * t21 + (-t19 * t21 * ((t102 * t56 + t21 * (-s26 * t4 + t123)) *
     # s26 + t101 * t54) + t22 - t6 * t53 * t54) * MBL1011(1) * t86) + t
     #14
      t14 = 0.1e1_dp / t65
      t4 = t14 * (t86 * ((-t64 * (t2 * (((s26 * s34 - t42 * s26 - t34) *
     # s34 + t40 + t45 - t61) * s56 + t67 - t100 * t29) + t13 * (t21 + t
     #53) + t28 * (t13 + t53 + t4 - t55) + t68 + t71) * MBL1101(0) + (s3
     #4 * t121 + (s16 - s25 - s34) * s56 + t69 + t93) * (t2 * t67 - (-(t
     #80 + t53) * s34 + t50 * (s26 * t126 + t2 * t53 + t36)) * s34 - (-(
     #(s34 - s25) * s34 + t4 + t53 - t55) * s56 - t2 * (t45 + t40) + t58
     # + s34 * (s26 * t106 - (s25 + t35) * s34 + t19 * t34)) * s56 + t68
     # + t71) * MBL1111D2(0)) * t37 * t14 + t56 * MBL1011(0) * t46 + t15
     # * MBL1110(0) * t49) - I300s16(0))
      result = t37 * ((s25 * ((t33 - t12) * s25 + t6) + 14.0_dp * t32 + 
     #((25.0_dp * s16 + t7) * s26 + (s16 * t3 - t5) * s25 + t6) * s26 + 
     #((-t12 + t11) * s34 - s25 * (s16 - 10.0_dp * s25) - (25.0_dp * s26
     # + t8 + t9) * s26) * s34 + (42.0_dp * s26 + 25.0_dp * t16 + t7) * 
     #t28 + t2 * (s25 * (-s25 + t25) + 25.0_dp * t17 + 21.0_dp * t21 + t
     #26 + t27 + t7 * s26) * s56) * t86 + (s25 * ((-t5 + t72) * s25 + t8
     #2) + t1 * t32 + ((t84 + t51) * s26 + t82 + t83) * s26 + (t85 * s34
     # - s25 * (s16 - t9) - (11.0_dp * s26 + t20 + t5) * s26) * s34 + ((
     #18.0_dp * s26 + 11.0_dp * t16 + t51) * s56 + 22.0_dp * t17 + t26 +
     # t27 + 14.0_dp * t66 + t73 + t83) * s56) * I300s16(1) + t49 * t10)
     # / 12.0_dp - t4 / 6.0_dp

           intHLs16s25s26s34s56x1141D6eps0 = result
       end function intHLs16s25s26s34s56x1141D6eps0

       function intHLs16s25s26s34s56x1141D6eps1()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1141D6eps1
           complex(dp) ::  t1,t10,t11,t12,t13,t2,t3,t4,t5,t6,t7,t8
           complex(dp) ::  t9

           complex(dp) :: result

      t1 = -11.0_dp
      t2 = s25 ** 2.0_dp
      t3 = 2.0_dp * s25
      t4 = 3.0_dp
      t5 = -14.0_dp
      t6 = 18.0_dp * s26
      t7 = s26 ** 2.0_dp
      t8 = 6.0_dp
      t9 = s26 + s56
      t10 = s25 + s26 - s34 + s56
      t11 = s25 + s26 + s56
      t11 = 0.1e1_dp / t11
      t12 = 0.1e1_dp / s16
      t9 = 0.1e1_dp / t9 ** 2.0_dp
      t10 = 0.1e1_dp / t10 ** 2.0_dp
      t13 = 0.1e1_dp / 0.12e2_dp
      result = t11 * (-(s34 ** 3.0_dp * MBL1110(1) * t10 + s25 * t2 * MB
     #L1011(1) * t9) * t12 + I300s16(1)) / 6.0_dp + t13 * ((-2.0_dp * t2
     # + 7.0_dp * t7) * s25 + t8 * (s56 ** 3.0_dp + s26 * t7) - s26 * t2
     # + ((s26 * t4 - t3) * s34 + (s26 * t1 - t3) * s26 + 4.0_dp * t2) *
     # s34 + ((s34 * t1 + 7.0_dp * s25 + t6) * s56 - (s25 * t5 - t6) * s
     #26 - (-s34 * t4 + 22.0_dp * s26 + t3) * s34 - t2) * s56) * t12 * t
     #9 * t10

           intHLs16s25s26s34s56x1141D6eps1 = result
       end function intHLs16s25s26s34s56x1141D6eps1

       function intHLs16s25s26s34s56x1211D2eps0()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1211D2eps0
           complex(dp) ::  t1

           complex(dp) :: result

      t1 = s25 + s26 + s56
      t1 = 0.1e1_dp / t1
      result = -0.1e1_dp / s16 * ((s25 * MBL1101(0) + (s26 + s56) * MBL1
     #111D2(0)) * t1 - MBL1110(0))

           intHLs16s25s26s34s56x1211D2eps0 = result
       end function intHLs16s25s26s34s56x1211D2eps0

       function intHLs16s25s26s34s56x1211D2eps1()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1211D2eps1

           complex(dp) :: result

      result = MBL1110(1) / s16

           intHLs16s25s26s34s56x1211D2eps1 = result
       end function intHLs16s25s26s34s56x1211D2eps1

       function intHLs16s25s26s34s56x1211D4eps0()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1211D4eps0
           complex(dp) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           complex(dp) ::  t20,t21,t22,t23,t24,t3,t4,t5,t6,t7,t8,t9

           complex(dp) :: result

      t1 = s16 - s25
      t2 = 2.0_dp * s16
      t3 = s25 + t2
      t4 = s16 * s25
      t5 = s25 ** 2.0_dp
      t6 = s16 ** 2.0_dp
      t7 = s34 ** 2.0_dp
      t8 = (s26 * t1 - t5) * s16
      t9 = t6 * s25
      t10 = s26 * t7
      t11 = (s26 * t3 + t4) * s34
      t12 = t7 + t6 + t5
      t13 = s16 + s25
      t14 = -2.0_dp * s34 * t13 + t12 - 2.0_dp * t4
      t15 = s16 + s25 - s34
      t16 = (-s34 + t13 + s26) * s26 + (2.0_dp * s26 + t15 + s56) * s56 
     #+ t4
      t17 = s26 + s56
      t18 = s26 * s34
      t19 = -3.0_dp * s16
      t20 = 2.0_dp * s25
      t21 = s25 - t2
      t22 = -4.0_dp
      t14 = 0.1e1_dp / t14
      t23 = 0.1e1_dp / t16
      t16 = 0.1e1_dp / t16
      t24 = 0.1e1_dp / 0.2e1_dp
      result = t24 * (t17 * (s25 + s26 + s56) * MBL1111D2(0) * t23 + (((
     #((-t20 + t19) * s25 + t6) * s26 + (s26 * t21 + t18 - t4) * s34 + (
     #(s26 * t22 + t19 - t20) * s25 + (t21 + s34) * s34 + t6) * s56 + t9
     # - 2.0_dp * s25 * (s26 ** 2.0_dp + s56 ** 2.0_dp) + t19 * t5) * MB
     #L1110(1) - s25 * (t2 * s25 - (s34 - t13) * s26 + s56 * t15) * MBL1
     #101(0)) * s34 + (-s16 * I300s16(1) + MBL1001(0)) * ((-s25 * s34 - 
     #t2 * s34 - t4 + t6 + t7) * s56 + t10 - t11 + t8 + t9) + t17 * (-t2
     # * t5 + (t6 + t5) * s25 + t1 ** 2.0_dp * s26 + (-s25 * t13 - t2 * 
     #s26 + t18) * s34 + (t12 - 2.0_dp * s16 * (s25 + s34)) * s56) * MBL
     #1011(1)) * t16 * t14) - ((s16 * t1 + (-t3 + s34) * s34) * s56 + t8
     # + t9 + t10 - t11) * t14 * t16

           intHLs16s25s26s34s56x1211D4eps0 = result
       end function intHLs16s25s26s34s56x1211D4eps0

       function intHLs16s25s26s34s56x1211D4eps1()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1211D4eps1

           complex(dp) :: result

      result = 0.0_dp

           intHLs16s25s26s34s56x1211D4eps1 = result
       end function intHLs16s25s26s34s56x1211D4eps1

       function intHLs16s25s26s34s56x1212D4eps0()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1212D4eps0
           complex(dp) ::  t1,t10,t11,t12,t13,t14,t2,t3,t4,t5,t6,t7
           complex(dp) ::  t8,t9

           complex(dp) :: result

      t1 = s16 - s25
      t2 = s25 * t1 + s26 * t1 + (s25 - s26) * s34 + (s16 - s25 - s34) *
     # s56
      t3 = s16 + s25
      t4 = s16 * s25
      t5 = s16 ** 2.0_dp
      t6 = s25 ** 2.0_dp
      t7 = s34 ** 2.0_dp + t5 + t6
      t8 = 2.0_dp
      t9 = -t8 * (s34 * t3 + t4) + t7
      t10 = s16 + s25 - s34
      t11 = (-s34 + t3 + s26) * s26 + (s26 * t8 + s56 + t10) * s56 + t4
      t12 = s26 + s56
      t13 = t8 * t6
      t14 = s25 * s34
      t4 = t4 * t8
      t9 = 0.1e1_dp / t9
      t11 = 0.1e1_dp / t11
      result = t8 * t2 * t9 * t11 - t11 * (-(s25 + s26 + s56) * MBL1111D
     #2(0) + t9 * (t2 * MBL1001(0) + (-s34 * (s25 + s26 - s34 + s56) * (
     #t4 - (s34 - t3) * s26 + s56 * t10) * MBL1110(1) + t12 * (-t13 * s1
     #6 + s25 * t7 + t1 ** 2.0_dp * s26 - (s26 * t3 + t13) * s34 + (-s16
     # * s34 - t14 - t4 + t5 + t6) * s56) * MBL1011(1)) / s16 - s16 * t2
     # * I300s16(1) - t14 * (t12 * t8 + t10) * MBL1101(0)))

           intHLs16s25s26s34s56x1212D4eps0 = result
       end function intHLs16s25s26s34s56x1212D4eps0

       function intHLs16s25s26s34s56x1212D4eps1()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1212D4eps1

           complex(dp) :: result

      result = 0.0_dp

           intHLs16s25s26s34s56x1212D4eps1 = result
       end function intHLs16s25s26s34s56x1212D4eps1

       function intHLs16s25s26s34s56x1213D6eps0()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1213D6eps0
           complex(dp) ::  t1,t10,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109
           complex(dp) ::  t11,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t12
           complex(dp) ::  t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t13,t130
           complex(dp) ::  t131,t132,t133,t134,t135,t136,t137,t138,t139,t14,t140,t141
           complex(dp) ::  t142,t143,t144,t145,t146,t147,t148,t149,t15,t150,t151,t152
           complex(dp) ::  t153,t154,t155,t16,t17,t18,t19,t2,t20,t21,t22,t23
           complex(dp) ::  t24,t25,t26,t27,t28,t29,t3,t30,t31,t32,t33,t34
           complex(dp) ::  t35,t36,t37,t38,t39,t4,t40,t41,t42,t43,t44,t45
           complex(dp) ::  t46,t47,t48,t49,t5,t50,t51,t52,t53,t54,t55,t56
           complex(dp) ::  t57,t58,t59,t6,t60,t61,t62,t63,t64,t65,t66,t67
           complex(dp) ::  t68,t69,t7,t70,t71,t72,t73,t74,t75,t76,t77,t78
           complex(dp) ::  t79,t8,t80,t81,t82,t83,t84,t85,t86,t87,t88,t89
           complex(dp) ::  t9,t90,t91,t92,t93,t94,t95,t96,t97,t98,t99

           complex(dp) :: result

      t1 = s16 - s25
      t2 = 11.0_dp
      t3 = s16 ** 2.0_dp
      t4 = t3 ** 2.0_dp
      t5 = s16 * t4
      t6 = s16 * t3
      t7 = 6.0_dp * t3
      t8 = 3.0_dp
      t9 = 2.0_dp
      t10 = 7.0_dp * s16
      t11 = t9 * s25
      t12 = t8 * t3
      t13 = (-t10 - t11) * s25 + t12
      t14 = s16 + s25
      t15 = s25 + s26
      t16 = t8 * s16
      t17 = 5.0_dp
      t18 = t17 * s25
      t19 = t9 * t3
      t20 = 12.0_dp * s16
      t21 = t8 * s25
      t22 = t17 * t3
      t23 = (t20 + t21) * s25 + t22
      t24 = s25 ** 2.0_dp
      t25 = t24 ** 2.0_dp
      t26 = s25 * t24
      t27 = t26 * t25
      t28 = s25 * t25
      t29 = 20.0_dp * s16
      t30 = 10.0_dp * s25
      t31 = 63.0_dp * s16
      t32 = 30.0_dp * t3
      t33 = t17 * t6
      t34 = 6.0_dp * s25
      t35 = (23.0_dp * s16 + t34) * s25 + t3
      t36 = 8.0_dp * s16
      t37 = 7.0_dp * t3
      t38 = (t36 + t21) * s25 - t37
      t39 = s26 ** 2.0_dp
      t40 = s26 * t39
      t41 = 30.0_dp * s25
      t42 = s16 * s25
      t43 = 18.0_dp
      t44 = 8.0_dp * s25
      t45 = s16 * t43
      t46 = (-t44 - t45) * s25 + t7
      t47 = -s25 + t16
      t48 = s25 * ((47.0_dp * s16 + t30) * s25 + 17.0_dp * t3)
      t49 = s16 * t24
      t50 = 13.0_dp * s16
      t51 = -t16 + t11
      t52 = 4.0_dp
      t53 = 17.0_dp * s16
      t54 = t52 * s25
      t55 = 12.0_dp * t3
      t56 = t8 * t6
      t57 = 42.0_dp
      t58 = 6.0_dp * t6
      t59 = t8 * s26
      t60 = t59 * t23
      t61 = t8 * t38
      t62 = t52 * s26
      t63 = 6.0_dp * s16
      t64 = t63 * s26
      t65 = s25 * s34
      t66 = t1 ** 2.0_dp
      t67 = t66 ** 2.0_dp
      t68 = t1 * t67
      t69 = t1 * t66
      t70 = t59 * t14
      t71 = t69 * s26
      t72 = t4 * t24
      t73 = -34.0_dp
      t74 = 24.0_dp * t3
      t75 = 12.0_dp * s25
      t76 = 9.0_dp * s26
      t77 = s16 * t25
      t78 = t24 - t3
      t79 = t1 * t15
      t80 = t14 * s26
      t73 = -(-t9 * t42 * t78 + ((s34 * t47 + t38) * s34 + t1 * t23) * s
     #34 + t25 - t4) * s56 - s25 * (((t24 * t9 - t3 * t43) * s25 + 28.0_
     #dp * t6) * s25 - 16.0_dp * t4) - t8 * (-t71 * t14 + t5) - (((s34 *
     # t51 + s16 * (t63 + t76) + (-t59 - t44 - t45) * s25) * s34 + s25 *
     # ((46.0_dp * s16 + t75) * s25 + t19) + s26 * t61) * s34 - t1 * (((
     #s16 * t73 - t44) * s25 - t74) * s25 + t58 - t60)) * s34 - t77
      t60 = t73 * s56 - ((t54 * t35 * s26 + t33 * s25 + (((t50 + t18 + t
     #62) * s25 - t64 - t65) * s34 - s26 * (-t8 * t47 * s26 - t9 * t46) 
     #- t48) * s34 + t32 * t24 + t31 * t26 + t39 * t61 + 10.0_dp * t25) 
     #* s34 - t1 * (s25 * (((-s16 * t57 - t18) * s25 - 19.0_dp * t3) * s
     #25 + t58) + (t52 * (((-t53 - t54) * s25 - t55) * s25 + t56) - t60)
     # * s26)) * s34 - t24 * t25 - t71 * (t13 * t9 - t70) - t36 * t28 + 
     #29.0_dp * t72 + 36.0_dp * t3 * t25 - 50.0_dp * t6 * t26 - t34 * t5
      t61 = s34 ** 2.0_dp
      t73 = t61 ** 2.0_dp
      t81 = s34 * t73
      t82 = s34 * t61
      t83 = t24 + t3 + t61
      t84 = s34 * t14
      t85 = -t9 * (t84 + t42) + t83
      t86 = s16 + s25 - s34
      t87 = t9 * s26
      t88 = (-s34 + t14 + s26) * s26 + (t86 + t87 + s56) * s56 + t42
      t89 = t9 * s16
      t90 = s16 - t54
      t91 = 9.0_dp * s16
      t92 = t52 * s16
      t93 = t34 + t92
      t94 = t42 * t14
      t95 = 8.0_dp * t26 + (t93 * s26 + (t75 + t91) * s25 - t3) * s26 + 
     #t94
      t96 = 15.0_dp * s25
      t97 = s16 + t11
      t98 = 22.0_dp * s25
      t99 = (-s16 + t98) * s25 - t12
      t100 = t17 * s16
      t101 = -t100 + t21
      t102 = t24 * ((t10 + t75) * s25 + t12) + ((s26 * t101 + t99) * s26
     # - s25 * (s16 - t96) * t97) * s26
      t103 = t24 * (t10 + t44) + t39 * (t34 - t100) - t9 * s26 * (-s25 *
     # (t30 + t16) + t39)
      t104 = t9 * (t24 - t39) + t18 * s26
      t99 = t9 * t99
      t105 = t8 * t101
      t106 = t3 * s25
      t107 = (-t41 - t50) * t24
      t108 = -t18 + t62
      t109 = s26 * t90 * t66
      t110 = s34 * (t24 * (-t30 - t100) + s26 * (-0.3e1_dp / 0.2e1_dp * 
     #t93 * s26 + (-s25 * t43 - t50) * s25 + t3))
      t111 = ((t34 - t100) * s26 - t8 * (t39 - t42) + 10.0_dp * t24) * t
     #82
      t112 = t4 * s25
      t113 = (t53 - t18) * s25
      t114 = 21.0_dp * t3
      t115 = t2 * t6
      t116 = t69 * t39
      t117 = t116 * t8
      t118 = 12.0_dp * s26
      t119 = 13.0_dp * s25
      t120 = s25 * (((-t50 + t54) * s25 + 15.0_dp * t3) * s25 - 7.0_dp *
     # t6) - t73 * t9 - ((-(-6.0_dp * s25 + 6.0_dp * s26 + t100) * s34 +
     # s16 * (15.0_dp * s26 + t16) + (s16 - t98 - t76) * s25) * s34 + t1
     # * (-t43 * s25 * t15 + (s16 - t118 - t119) * s16)) * s34 + t4 - t5
     #9 * t69
      t121 = t8 * t42 * t1
      t93 = t9 * s34 * (t1 * t93 / 2.0_dp + t61)
      t122 = s56 ** 2.0_dp
      t123 = s56 * t122
      t124 = s25 * (-t18 + t89)
      t125 = (t10 - t11) * s25
      t126 = 9.0_dp * t3
      t127 = s16 * t26
      t128 = t26 * ((-t125 + t126) * s25 - t33) - ((-(t61 * t101 + t121 
     #+ t26 - t6 + t93) * s56 - t120) * s56 - t24 * ((-t113 + t114) * s2
     #5 - t115) + t9 * (t1 * (t110 - t109) + t111 - t112) + t117 - t61 *
     # (-t61 * t108 + (t105 * s26 + t99) * s26 - t106 - t107)) * s56 + t
     #102 * t61 - t103 * t82 + t104 * t73 - t40 * t69 + t72 + t79 * t95 
     #* s34 + t116 * t90 + t124 * t71
      t129 = s16 - t21
      t130 = t14 * t39
      t131 = t3 * t52
      t132 = t24 + t3
      t133 = 6.0_dp * t42
      t134 = -t132 * t17 + t133
      t135 = -t17 * s16 * t14 + 12.0_dp * t24
      t136 = (t96 - t63) * s25
      t137 = t136 + t22
      t138 = t2 * s25
      t139 = -t138 + t89
      t140 = s25 * (s16 - t18)
      t141 = t9 * t97
      t142 = s25 + t87
      t143 = 9.0_dp * t6
      t144 = ((51.0_dp * s16 - 63.0_dp * s25) * s25 + 33.0_dp * t3) * s2
     #5 - t143
      t143 = s25 * (((29.0_dp * s16 - 26.0_dp * s25) * s25 + t74) * s25 
     #+ t143)
      t145 = t6 * t43
      t146 = -t134 * t52
      t147 = t26 * t1
      t148 = t137 * t8
      t149 = t129 * t67
      t150 = 7.0_dp * s25
      t151 = 22.0_dp * t3
      t152 = t52 * t67
      t119 = (-t119 - t92) * s25
      t153 = -t26 + t6
      t154 = t147 * t3
      t155 = t15 ** 2.0_dp
      t53 = ((t134 * t61 + t52 * (t42 * t132 + t84 * t66) - t25 - t4 - t
     #7 * t24 + t9 * t14 * t82) * s56 + s25 * ((((t50 - t21) * s25 - t15
     #1) * s25 + t145) * s25 - 7.0_dp * t4) + (((t22 + 8.0_dp * t80 + t1
     #36) * s34 + ((-21.0_dp * s25 + t53) * s25 + t2 * t3) * s25 - t56 -
     # t146 * s26) * s34 - t66 * (t3 - 16.0_dp * t80 + t119)) * s34 + t5
     # - t152 * s26 - t141 * t73) * s56 - t8 * ((t28 - t5) * s25 - t149 
     #* s26) + (((((-t118 - t138 + t89) * s25 - t64) * s34 + (t148 + 12.
     #0_dp * t80) * s26 + 24.0_dp * t26 - 10.0_dp * t94) * s34 + (6.0_dp
     # * t134 * s26 + t144) * s26 + t143) * s34 - t66 * (s25 * ((-14.0_d
     #p * s25 + t36) * s25 + t7) + (t8 * (t3 + t119) - 24.0_dp * t80) * 
     #s26)) * s34 - 15.0_dp * t49 * t153 + 30.0_dp * t154 - 6.0_dp * t67
     # * t39 + t11 * t81
      t56 = t67 * s26 * (t24 * t47 + t39 * (t129 - s26))
      t64 = t66 * t155 * (-t52 * (-t94 + t130) + t1 * (s16 + t18) * s26)
     # * s34
      t67 = t15 * ((t87 * (((t91 - t18) * s25 + t131) * s25 + t58) + t13
     #1 * t24) * s25 + t39 * (t134 * s26 + (t2 * t24 - t12) * s16 - 16.0
     #_dp * s25 * t78) + 6.0_dp * t49 * t132) * t61
      t5 = s16 * t27 + (t53 * s56 + t8 * (t149 * t39 + t24 * t5) + (((((
     #t87 * t139 + t140) * s25 - 6.0_dp * t97 * t39 + t65 * (s25 + t62))
     # * s34 + ((t148 + 8.0_dp * t80) * s26 + t54 * t135) * s26 - 10.0_d
     #p * t147 - 14.0_dp * t3 * t24) * s34 + t143 * t87 + t144 * t39 + t
     #145 * t24 - t146 * t40 + t55 * t26 - 10.0_dp * t28 + 24.0_dp * t77
     #) * s34 - t66 * t15 * (s25 * ((t20 - t18) * s25 + t126) + ((-23.0_
     #dp * s25 + t92) * s25 + t12 - 16.0_dp * t80) * s26)) * s34 - t27 +
     # t127 * (((t150 - t45) * s25 + t151) * s25 - 13.0_dp * t6) - t152 
     #* t40 + t34 * t68 * s26) * s56 + t26 * t5 + t73 * (s25 * ((s26 * t
     #139 + t140) * s26 + t49) - t141 * t40) + t82 * (t137 * t40 - t52 *
     # t127 * t14 + t9 * s26 * ((s25 * t135 + t130) * s26 + t24 * (-t17 
     #* s25 * t1 - t37))) + t56 + s25 * s26 * t142 * t81 - t64 + t67
      t20 = (t100 + t54) * s25
      t27 = t17 * s26
      t28 = (t36 + t150) * s25
      t37 = s25 + s26 - s34 + s56
      t40 = 9.0_dp * s25
      t53 = 10.0_dp * t42
      t55 = t132 + t53
      t43 = t42 * t43 + t78 * t8
      t56 = -t26 + t6
      t57 = (t9 * s16 * t153 + (t61 * (t15 * t8 + s16 - s34) + ((19.0_dp
     # * s16 + s25) * s25 + t114) * s25 - t33 + t59 * t55) * s34 - 6.0_d
     #p * t106 * t1 + t8 * (-t80 * t9 - t133 - t78) * t61) * s56 + t52 *
     # s16 * (s25 * t56 + t71) - t9 * t61 * (s26 * t61 + (t43 + t70) * s
     #26 + t42 * (t40 + t100)) + s34 * (t61 * ((t59 + t34 + t89) * s26 +
     # 7.0_dp * t42) - s25 * t6 + (((38.0_dp * s16 + t11) * s25 + t3 * t
     #57) * s25 - 10.0_dp * t6) * s26 + t49 * (22.0_dp * s16 + t96) + t8
     # * t55 * t39) - 12.0_dp * t3 * t24 * t1
      t59 = s25 + s26 + s56
      t62 = 0.1e1_dp / t88 ** 2.0_dp
      t64 = 0.1e1_dp / t85 ** 2.0_dp
      t11 = s34 * t37 * (t9 * s16 * (t24 * t56 + t116) - 6.0_dp * t154 +
     # (((-s34 * t39 + t39 * (s26 + s16 + t21) + t42 * (7.0_dp * s26 + t
     #11)) * s34 - t39 * t43 - t9 * s26 * (t130 + t42 * (t40 + t100)) - 
     #t49 * (t36 + t34)) * s34 + t123 * (-t84 * t9 + t53 + t83) + t15 * 
     #(s26 * (t55 * s26 + s16 * ((t29 + t40) * s25 - t22)) + t42 * ((t34
     # + t89) * s25 + t131))) * s34 + s56 * t57 + t71 * t42 * t52) * MBL
     #1110(1)
      t2 = (-t52 * t3 * t25 * t132 + (t39 * t68 * t8 + t58 * t25) * s25 
     #+ t5) * MBL1011(1) + ((((-s26 * (s25 * (t50 + t18) + s26 * t51) - 
     #t49 * t17 + t65 * s26) * s34 - ((s26 * t47 + t46) * s26 - t48) * s
     #26 - t49 * (-16.0_dp * s16 - 20.0_dp * s25)) * s34 - s25 * (t42 * 
     #((14.0_dp * s16 + t41) * s25 + 8.0_dp * t3) + t9 * t35 * t39) - s2
     #6 * (s25 * (((t30 + t31) * s25 + t32) * s25 + t33) + t38 * t39)) *
     # s34 + t79 * (-s26 * (t23 * s26 - (s25 + t16) * ((-t10 - t18) * s2
     #5 + t19)) - t29 * t26)) * s34 + t60 * s56 + t71 * (-s25 * ((-s16 *
     # t2 - s25) * s25 + t7) + (-t13 + t80) * s26) - t49 * ((((-t45 + t1
     #8) * s25 + t74) * s25 - 14.0_dp * t6) * s25 + t4 * t8) - t11
      t3 = t62 * ((t65 * ((((s34 * t142 - t52 * s25 * t14 - (t44 + t63 +
     # t27) * s26) * s34 + t17 * t94 + t39 * (t138 + t91) + t8 * s26 * (
     #t20 + t19 + t39) + 6.0_dp * t26) * s34 - t15 * ((t20 + t12 + t70) 
     #* s26 + t26 * t52 + t9 * s16 * ((s16 - s25) * s25 + t3))) * s34 + 
     #((-t17 * t82 + t121 - (-(9.0_dp * s16 + 9.0_dp * s26 + t138) * s34
     # + t12 + t28 + 9.0_dp * t80) * s34 + t26 - t6) * s56 - t9 * (s25 *
     # t153 + t82 * (t27 + t54 + t16 - s34) + t71) - s34 * (-((15.0_dp *
     # s16 + t75) * s25 + (t45 + t98 + t76) * s26 + t7) * s34 + s25 * ((
     #t44 + t16) * s25 + t22) + t9 * ((t28 + t12) * s26 + t6) + 9.0_dp *
     # t130) + 6.0_dp * t49 * t1) * s56 - t155 * t69) - t21 * t86 * t61 
     #* t123) * MBL1101(0) * t64 + t59 ** 2.0_dp * t37 * MBL1111D2(0))
      result = t64 * t62 * ((-s16 * (t61 * ((s34 * t104 - t103) * s34 + 
     #t102) + t72) + (s16 * (t24 * ((t113 - t114) * s25 + t115) + t117) 
     #+ t9 * s16 * (t1 * (t110 - t109) + t111 - t112) + t61 * (s16 * ((-
     #t105 * s26 - t99) * s26 + t106 + t107) + t61 * s16 * t108)) * s56 
     #- t79 * s16 * t95 * s34 - s16 * t120 * t122 + s16 * (-t61 * t101 -
     # t121 - t26 + t6 - t93) * t123 - t71 * s16 * ((t90 - s26) * s26 + 
     #t124) + t127 * ((t125 - t126) * s25 + t33)) * I300s16(1) + t128 * 
     #MBL1001(0) + t2 / s16) / 4.0_dp + t3 / 2.0_dp

           intHLs16s25s26s34s56x1213D6eps0 = result
       end function intHLs16s25s26s34s56x1213D6eps0

       function intHLs16s25s26s34s56x1213D6eps1()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1213D6eps1

           complex(dp) :: result

      result = 0.0_dp

           intHLs16s25s26s34s56x1213D6eps1 = result
       end function intHLs16s25s26s34s56x1213D6eps1

       function intHLs16s25s26s34s56x1221D4eps0()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1221D4eps0
           complex(dp) ::  t1,t10,t11,t2,t3,t4,t5,t6,t7,t8,t9

           complex(dp) :: result

      t1 = s26 + s56
      t2 = 2.0_dp
      t3 = s16 * s25
      t4 = (-s34 + s16 + s25 + s26) * s26 + (s26 * t2 + s16 + s25 - s34 
     #+ s56) * s56 + t3
      t5 = s25 + s26 + s56
      t6 = s25 + s26 - s34 + s56
      t7 = t2 * s16
      t8 = s16 - s34
      t9 = 0.1e1_dp / t5
      t10 = 0.1e1_dp / s16
      t6 = 0.1e1_dp / t6
      t4 = 0.1e1_dp / t4
      t11 = t10 * t4
      result = -t11 * t2 * t1 - t6 * (t4 * (t5 * MBL1001(0) + (s34 * (t7
     # * s25 + (-s34 + s25 + t7 + s26) * s26 + (t2 * (s26 + s16) + s25 -
     # s34 + s56) * s56) * MBL1110(1) + t1 * ((s16 + s34 - s25) * s25 + 
     #t1 * (s16 - s25)) * MBL1011(1)) * t10) + s34 * MBL1110(0) * t10) -
     # t1 * I300s16(1) * t4 + t11 * ((-t2 * ((s26 * t8 + t3) * s56 + t3 
     #* s26) - s16 * s25 ** 2.0_dp - t8 * s26 ** 2.0_dp - t8 * s56 ** 2.
     #0_dp) * MBL1111D2(0) + s25 * s34 * t1 * MBL1101(0)) * t9

           intHLs16s25s26s34s56x1221D4eps0 = result
       end function intHLs16s25s26s34s56x1221D4eps0

       function intHLs16s25s26s34s56x1221D4eps1()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1221D4eps1
           complex(dp) ::  t1

           complex(dp) :: result

      t1 = s25 + s26 - s34 + s56
      t1 = 0.1e1_dp / t1
      result = -0.1e1_dp / s16 * t1 * (s34 * MBL1110(1) + 1.0_dp)

           intHLs16s25s26s34s56x1221D4eps1 = result
       end function intHLs16s25s26s34s56x1221D4eps1

       function intHLs16s25s26s34s56x1222D6eps0()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1222D6eps0
           complex(dp) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           complex(dp) ::  t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30
           complex(dp) ::  t31,t32,t33,t34,t35,t36,t37,t38,t39,t4,t40,t41
           complex(dp) ::  t42,t43,t44,t45,t46,t47,t48,t49,t5,t50,t51,t52
           complex(dp) ::  t53,t54,t55,t56,t57,t58,t59,t6,t60,t61,t62,t63
           complex(dp) ::  t64,t65,t66,t67,t7,t8,t9

           complex(dp) :: result

      t1 = s16 - s25
      t2 = 5.0_dp
      t3 = 2.0_dp
      t4 = t3 * s25
      t5 = t2 * s16
      t6 = t5 + t4
      t7 = s25 + s26
      t8 = 6.0_dp
      t9 = s16 * t8
      t10 = 3.0_dp * s25
      t11 = t10 + t9
      t12 = s16 ** 2.0_dp
      t13 = s16 * t12
      t14 = 11.0_dp * t12
      t15 = s16 + s25
      t16 = 10.0_dp * s25
      t17 = t15 * s26
      t18 = s16 * s25
      t19 = 8.0_dp * s16
      t20 = 4.0_dp
      t21 = t20 * s25
      t22 = 7.0_dp * s16
      t23 = s25 ** 2.0_dp
      t24 = t23 ** 2.0_dp
      t25 = s25 * t24
      t26 = s25 * t23
      t27 = t5 * t23
      t28 = t2 * s25
      t29 = (t19 + t28) * s25
      t30 = s26 ** 2.0_dp
      t31 = t3 * s26
      t32 = (s25 + t31) * s34
      t33 = t1 ** 2.0_dp
      t34 = 3.0_dp * s26
      t35 = t33 * s26
      t36 = s16 * t23
      t37 = t3 * t23
      t38 = t2 * t13
      t39 = 3.0_dp * t35
      t40 = s34 ** 2.0_dp
      t41 = s34 * t40
      t42 = -t3 * (s34 * t15 + t18) + t12 + t23 + t40
      t43 = s26 * t7
      t44 = t12 + t23
      t45 = (-s34 + t15 + s26) * s26 + (t31 + s16 + s25 - s34 + s56) * s
     #56 + t18
      t46 = t3 * s16
      t47 = s25 + t46
      t48 = t47 * s26
      t49 = s25 * t15
      t50 = t20 * s16
      t51 = t36 * t20
      t52 = t3 * (s25 * t44 + t35) - s34 * (-t32 + (t50 + t4) * s26 + 3.
     #0_dp * t49) - t51
      t53 = (-t3 * s16 * (s25 + s34) + (s34 - s25) * s34 + t12 + t23) * 
     #s56
      t54 = t33 * t30
      t55 = t44 * t23
      t56 = (-(t43 + t23) * s34 + t7 * (t37 + t48 + t18)) * s34
      t57 = t3 * s25 * (-t36 + t35)
      t58 = t5 - t4
      t59 = -s25 + t46
      t60 = t20 * t12
      t61 = -s25 * t1 - t60
      t62 = 3.0_dp * s16
      t63 = -t4 + t62
      t64 = t23 * t1
      t65 = s16 * t15
      t66 = 9.0_dp * s16
      t67 = t40 + t12
      t5 = -t3 * t12 * t24 + s16 * t25 + s34 * (((t46 * t30 + 3.0_dp * t
     #64) * s26 + s25 * (t30 * t63 + t36)) * s34 + t7 * (-t3 * t36 * t15
     # + s26 * (s25 * ((-t22 + t10) * s25 - t60) + s26 * t61))) + (((t3 
     #* s16 * t67 + s34 * t61 - t26 + t18 * (t21 - t5)) * s56 + s25 * ((
     #(-t4 + t66) * s25 - 12.0_dp * t12) * s25 + t38) + s34 * ((s25 * t6
     #3 + t9 * s26) * s34 - 8.0_dp * t18 * t15 + t20 * t26 + t34 * t61) 
     #+ t39 * t59) * s56 + t36 * ((s25 * t8 - t66) * s25 + t60) + 3.0_dp
     # * t54 * t59 + t4 * t35 * t58 - ((-t8 * t43 * s16 + t20 * t23 * s2
     #6 - s34 * t23 - 3.0_dp * t64) * s34 - s25 * ((-8.0_dp * t3 * t65 +
     # 8.0_dp * t23) * s26 - t18 * (9.0_dp * s25 + t9)) - 3.0_dp * t30 *
     # t61 - 3.0_dp * t24) * s34 - t25) * s56 + t13 * t26 + t54 * t58 * 
     #s25 + t33 * t59 * s26 * t30 + t23 * (t33 * (-s25 + t50) + t41) * s
     #26
      t15 = s16 - t10
      t25 = s16 + t4
      t33 = -s16 + s25 + s34
      t50 = t12 * s25
      t49 = t49 - t12
      t58 = s16 * t1
      t59 = s25 * s34
      t60 = s56 ** 2.0_dp
      t28 = -(t31 * t40 ** 2.0_dp - ((s26 * (t21 + t46 + t34) * s34 + t3
     # * (-s26 * t49 + t50) + t27 - 3.0_dp * t47 * t30) * s34 - t58 * t7
     # * (t46 - t28 - t34)) * s34) * s56 - t41 * (s34 * t30 + t36 * t3 -
     # t30 * (s26 + t25)) - ((s26 * t11 + (-s16 - t34 - t4 + s34) * s34 
     #+ t49) * s34 + t58 * (s16 - t21 - t34)) * s34 * t60 - t58 * (-t4 +
     # s16 - s26) * t7 ** 2.0_dp * s34 + t7 * (t51 + (-t48 + t65) * s26 
     #+ t50) * t40 + s34 * (-t46 * s34 + t12 - t18 + t40 - t59) * s56 * 
     #t60
      t46 = 0.1e1_dp / t42
      t45 = 0.1e1_dp / t45 ** 2.0_dp
      t1 = t59 * (-t23 * t63 + ((s26 + t4) * s34 + (-t21 - t62) * s25 + 
     #(-t25 * t3 + s26) * s26) * s34 + (t33 * s56 - t20 * s25 * (s16 + s
     #34) + t3 * (-s16 * s34 + t33 * s26) + 3.0_dp * t23 + t67) * s56 + 
     #t50 + t1 * s26 * (-s26 + t15)) * MBL1101(0)
      t4 = (s25 + s26 + s56) * (t37 + (-t15 + s26) * s26 + (s26 - t4) * 
     #s34 + (t31 + t10 - s16 + s34 + s56) * s56 - t18) * MBL1111D2(0)
      t15 = t46 * t45 * ((-t57 + (-t53 - t52) * s56 - t54 - t55 + t56) *
     # MBL1001(0) + s16 * (t57 - (-t53 - t52) * s56 + t54 + t55 - t56) *
     # I300s16(1))
      result = t45 * (t46 * ((-t5 * MBL1011(1) + t28 * MBL1110(1) + t2 *
     # t36 * t44 + ((t42 * s56 + (-8.0_dp * t12 + t37) * s25 - ((-t22 - 
     #t21 - t34 + s34) * s34 + t17 * t8 + t14 + t29) * s34 + t36 + t38 +
     # t39) * s56 + t35 * (t3 * t6 + t34) - 19.0_dp * t12 * t23 + t16 * 
     #t13 + t19 * t26 + ((-t32 + 3.0_dp * t30 + 3.0_dp * t23 + 8.0_dp * 
     #s25 * (s16 + s26) + 14.0_dp * s16 * s26) * s34 - s25 * ((16.0_dp *
     # s16 + t10) * s25 + 17.0_dp * t12) - (t3 * (t29 + t14) + t17 * t8)
     # * s26) * s34 + t24) * s56 + (t27 + (s25 * (t19 + t10) + (t22 + t2
     #1 + s26) * s26) * s26) * t40 - t43 * t41 + t35 * (s25 * (10.0_dp *
     # s16 + s25) + (t6 + s26) * s26) - t7 * ((s25 * t11 + t17 * t3 + t1
     #4) * s26 + t18 * (t16 + t9)) * s34 - 10.0_dp * t12 * t26) / s16 - 
     #t1) - t4) / 2.0_dp + t15

           intHLs16s25s26s34s56x1222D6eps0 = result
       end function intHLs16s25s26s34s56x1222D6eps0

       function intHLs16s25s26s34s56x1222D6eps1()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1222D6eps1

           complex(dp) :: result

      result = 0.0_dp

           intHLs16s25s26s34s56x1222D6eps1 = result
       end function intHLs16s25s26s34s56x1222D6eps1

       function intHLs16s25s26s34s56x1231D6eps0()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1231D6eps0
           complex(dp) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           complex(dp) ::  t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30
           complex(dp) ::  t31,t32,t33,t34,t35,t36,t37,t38,t39,t4,t40,t41
           complex(dp) ::  t42,t43,t44,t45,t46,t47,t48,t49,t5,t50,t51,t52
           complex(dp) ::  t53,t54,t55,t56,t57,t58,t59,t6,t60,t61,t62,t63
           complex(dp) ::  t64,t65,t66,t67,t68,t69,t7,t8,t9

           complex(dp) :: result

      t1 = s25 + s26 + s56
      t2 = s25 - s16
      t3 = 2.0_dp
      t4 = s25 + s26
      t5 = t3 * s26
      t6 = s25 - t5
      t7 = s26 + s16
      t8 = s26 ** 2.0_dp
      t9 = t8 ** 2.0_dp
      t10 = s26 * t8
      t11 = s25 ** 2.0_dp
      t12 = t11 ** 2.0_dp
      t13 = s25 * t11
      t14 = 3.0_dp
      t15 = t3 * s34
      t16 = 6.0_dp * s25
      t17 = s26 * t4
      t18 = t17 * s25
      t19 = t6 * s34
      t20 = s25 + s26 - s34 + s56
      t21 = s16 + s25
      t22 = s16 * s25
      t23 = (-s34 + t21 + s26) * s26 + (t5 + s25 + s16 - s34 + s56) * s5
     #6 + t22
      t24 = 14.0_dp * s16
      t25 = t3 * s16
      t26 = s16 - s34
      t27 = s26 * t26
      t28 = 8.0_dp * s25
      t29 = 9.0_dp * s26
      t30 = 4.0_dp
      t31 = t30 * s25
      t32 = s56 ** 2.0_dp
      t33 = t32 ** 2.0_dp
      t34 = s56 * t32
      t35 = 6.0_dp * s16
      t36 = t14 * s16
      t37 = t3 * s25
      t38 = t14 * s26
      t39 = 5.0_dp * s16
      t40 = (t38 + t37) * s26
      t41 = 5.0_dp * s26
      t42 = t14 * t8
      t43 = 5.0_dp * s25
      t44 = s16 * s26
      t45 = s16 * t11
      t46 = s16 ** 2.0_dp
      t47 = t46 + t11
      t48 = 7.0_dp
      t49 = t48 * s16
      t50 = -t49 * s25 + t47
      t51 = t14 * t47
      t52 = 16.0_dp * s16
      t53 = -t52 * s25 + t51
      t54 = t30 * s16
      t55 = s25 * t14
      t56 = -t31 + t36
      t57 = 9.0_dp * t46
      t58 = s25 * t47
      t59 = 12.0_dp
      t60 = s26 * t2
      t61 = t14 * t58
      t62 = s16 * t59
      t63 = t50 * t30
      t64 = 9.0_dp * s16
      t65 = s34 ** 2.0_dp
      t66 = s34 * t65
      t67 = 10.0_dp * s16
      t68 = t3 * s34 * t2 - (s25 * t48 + t41) * s16 + t47
      t69 = t4 * s34
      t51 = t69 * (((-t56 * s26 - s25 * (t55 - t54)) * s26 + t45) * s34 
     #+ t4 * ((-t28 * s16 + t51) * s26 - t3 * (t2 * t8 + t45)))
      t6 = s25 * (-t11 * t50 + t6 * t66) * s26
      t6 = s16 * s25 * t12 + s16 * s26 * t9 + s16 * s56 * t33 - ((-t12 *
     # t48 - 5.0_dp * t9) * s16 + ((-s25 * (-s26 * t30 + s25) * s34 + t8
     # * t14 * t56 + ((-t39 + t55) * s25 - t5 * t21) * s25) * s34 - t4 *
     # ((-8.0_dp * t60 + (-20.0_dp * s16 + t43) * s25 + t57) * s26 + t61
     # - t62 * t11)) * s34 + t13 * t47 + t40 * s25 * t53 + t63 * t10) * 
     #s56 - t32 * (t14 * (s26 * t53 + t58) * s25 + s34 * (-(s25 * t21 + 
     #(s25 * t59 - t64) * s26) * s34 + s25 * ((t52 - t31) * s25 - 6.0_dp
     # * t46) + (t60 * t59 + (t62 + t55) * s25 - t57) * s26) + t37 * t66
     # - t52 * t13 - t67 * t10 + 6.0_dp * t50 * t8) - t33 * t68 - t34 * 
     #((-t67 * s26 + t63) * s26 + t61 + s34 * (s34 * t56 + (s25 + t54) *
     # s25 - t14 * t46 + 8.0_dp * t60) - t52 * t11) - t50 * t9 - t11 * t
     #53 * t8 - s25 * t53 * t10 + t51 + t6
      t9 = 6.0_dp * s26
      t21 = t59 * s26
      t47 = t4 ** 2.0_dp
      t50 = s16 * t4
      t51 = s16 * t2
      t52 = -t51 + t65
      t53 = 0.1e1_dp / t1
      t54 = 0.1e1_dp / s16
      t20 = 0.1e1_dp / t20 ** 2.0_dp
      t23 = 0.1e1_dp / t23 ** 2.0_dp
      t6 = t6 * MBL1011(1) + ((t25 * t47 ** 2.0_dp + t42 * t66) * s34 + 
     #((8.0_dp * s16 * t4 * t47 + t9 * t66) * s34 - t65 * (((t21 + t49) 
     #* s25 + 18.0_dp * s26 * t7) * s34 - t4 * ((21.0_dp * s16 + t16) * 
     #s26 + t59 * (t46 + t8) - t22))) * s56 + t32 * ((((10.0_dp * s25 + 
     #21.0_dp * s26 + t35) * s16 + t11 * t14 + 18.0_dp * t17) * s34 + t6
     #2 * t47) * s34 + t14 * t66 * (s34 - t9 - t36 - t37)) + t4 * ((-t9 
     #- t64) * s26 + t37 * s16) * t66 + s34 * (s34 * t14 + t25) * t33 + 
     #s34 * ((-6.0_dp * s34 + t21 + t49 + t16) * s34 + 8.0_dp * t50) * t
     #34 - t47 * ((-s26 * t48 + t31 - t35) * s16 - t42) * t65) * MBL1110
     #(1)
      t1 = t20 * ((t13 * t3 - ((s16 - t31 - s26) * s26 - s25 * (-t25 + t
     #43)) * s26 - (-(-t38 + t37) * s34 - t4 * (-t31 + t36 + t5)) * s34 
     #- ((s16 - t31 - t38 - t15 - s56) * s56 - (s34 * t30 + t28) * s26 -
     # t14 * (s34 * t26 + t8) + t3 * ((s16 + s34) * s25 + t44) - 5.0_dp 
     #* t11) * s56 - t45) * MBL1001(0) + t14 * t18 + ((-t2 * t3 + s26) *
     # t4 + t19) * s34 + ((t14 * t4 + s34 + s56) * s56 + (t3 * t7 - s25 
     #- t15) * s34 + t14 * (t11 + t8) + t16 * s26) * s56 + t10 + t13) * 
     #t1
      t1 = t23 * (-((t41 * s25 + t11 * t3) * s16 + (t39 * s25 + (t14 * (
     #s26 + s16 - s34) + s25 + s56) * s56 + 6.0_dp * t27 + t40) * s56 + 
     #t8 * (s25 + t36 + s26) - t42 * s34) * I300s16(1) + (t6 * t20 - t14
     # * (t34 + t10) - (s25 * (s25 + t24) - s34 * (s25 + 8.0_dp * s26)) 
     #* s26 - ((8.0_dp * t26 + t31 + t29) * s56 + (t29 + t28) * s26 + 16
     #.0_dp * t27 + s25 * (s25 - s34) + t24 * s25) * s56 - t35 * t11 - t
     #30 * (s25 + t25) * t8) * t54 + t1)
      result = t54 * (((t14 * ((-(-t5 * s16 * s34 - ((s26 + s34 + s25) *
     # s25 - t50) * s16 + s26 * t65) * s56 + s16 * t13 - t11 * t46 + t3 
     #* t44 * (s25 * t2 + t69) - t52 * t8) * s56 + t51 * t18) + s16 * t1
     #2 - t10 * t65 - t13 * t46 - t34 * (-t25 * s34 + t52) - t19 * s16 *
     # t47 + t51 * t10) * MBL1111D2(0) + s25 * s34 * (t3 * ((t27 + t22) 
     #* s56 + t22 * s26) + t26 * t32 + t26 * t8 + t45) * MBL1101(0)) * t
     #23 * t53 + t65 * MBL1110(0) * t20) / 2.0_dp + t1 / 4.0_dp

           intHLs16s25s26s34s56x1231D6eps0 = result
       end function intHLs16s25s26s34s56x1231D6eps0

       function intHLs16s25s26s34s56x1231D6eps1()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1231D6eps1
           complex(dp) ::  t1,t2

           complex(dp) :: result

      t1 = s25 + s26 - s34 + s56
      t2 = -0.1e1_dp / 0.4e1_dp
      t1 = 0.1e1_dp / t1 ** 2.0_dp
      result = 0.1e1_dp / s16 * t1 * (t2 * (s25 + s26 - 3.0_dp * s34 + s
     #56) + s34 ** 2.0_dp * MBL1110(1) / 2.0_dp)

           intHLs16s25s26s34s56x1231D6eps1 = result
       end function intHLs16s25s26s34s56x1231D6eps1

       function intHLs16s25s26s34s56x1311D4eps0()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1311D4eps0
           complex(dp) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           complex(dp) ::  t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30
           complex(dp) ::  t31,t32,t33,t34,t35,t36,t37,t38,t39,t4,t40,t41
           complex(dp) ::  t42,t5,t6,t7,t8,t9

           complex(dp) :: result

      t1 = s16 ** 2.0_dp
      t2 = s16 * t1
      t3 = s25 ** 2.0_dp
      t4 = s16 + s25
      t5 = 2.0_dp
      t6 = s16 * s25
      t7 = s34 * t4
      t8 = s34 ** 2.0_dp
      t9 = t3 + t1 + t8
      t10 = s26 * s34
      t11 = s16 * t3
      t12 = t1 * s25
      t13 = (t3 + t1) * s26
      t14 = (-s26 * t4 * t5 + t10 - t6) * s34 + (-t5 * t7 + t9) * s56 + 
     #t11 + t12 + t13
      t7 = -t5 * (t7 + t6) + t9
      t9 = t5 * s26
      t15 = (-s34 + t4 + s26) * s26 + (t9 + s16 + s25 - s34 + s56) * s56
     # + t6
      t16 = 3.0_dp * s16
      t17 = t5 * s25
      t18 = -t17 - t16
      t19 = t5 * s16
      t20 = s25 - t19
      t21 = 4.0_dp
      t22 = s56 ** 2.0_dp
      t23 = s26 ** 2.0_dp
      t24 = s25 + s26 - s34 + s56
      t25 = s26 + s56
      t26 = t4 * t3
      t27 = 3.0_dp * t1
      t28 = t27 * s25
      t29 = s16 - s25
      t30 = s25 * t29
      t31 = s25 + s26
      t32 = t4 * t8
      t33 = s26 * t29
      t34 = s16 + t17
      t35 = s25 + t19
      t36 = 5.0_dp * s25
      t37 = t35 * s26
      t38 = s16 * (t36 - t19)
      t39 = (7.0_dp * s16 + s25) * s25
      t40 = -9.0_dp
      t41 = 6.0_dp * s16
      t24 = 0.1e1_dp / t24
      t7 = 0.1e1_dp / t7
      t42 = 0.1e1_dp / s16
      t15 = 0.1e1_dp / t15
      t13 = t25 * (t11 * (t17 - t16) + s25 * t2 + (t32 - t28 + t26 + t2)
     # * s26 - (t5 * (t13 + t11) + s25 * (t23 + t1)) * s34 + (-t5 * ((s2
     #5 * t31 + t1) * s34 + t33 * s25) + t2 + t26 - t28 + t32) * s56 - t
     #23 * t30 - s25 * (s16 - s25 + s34) * t22) * MBL1011(1)
      t26 = s25 * ((t33 + t1) * s25 - s26 * t8 + (s26 * t34 + t6) * s34 
     #+ ((t34 - s34) * s34 + t30) * s56 - t11) * MBL1101(0)
      t10 = (-t5 * s25 * (t23 + t22) + (s25 * t18 + t1) * s26 + (s26 * t
     #20 + t10 - t6) * s34 + ((-s26 * t21 - t16 - t17) * s25 + (t20 + s3
     #4) * s34 + t1) * s56 + t12 - t16 * t3) * MBL1001(0) * t24
      t20 = (s25 + s26 + s56) * MBL1110(0) * t42 * t24
      t1 = t15 * (t7 * (t14 * I300s16(1) + (t24 * ((-t1 * s25 * t3 + (((
     #-s25 * s34 - t19 * s34 + t1 - t6 + t8) * s56 + (-t3 * t5 + 3.0_dp 
     #* t33) * s16 + (-3.0_dp * s25 * s26 - t41 * s26 + (3.0_dp * s16 + 
     #3.0_dp * s26 + t17 - s34) * s34 - 3.0_dp * t1 - t39) * s34 + t1 * 
     #t4) * s56 + t5 * s16 * (t33 * t34 + t12) - ((-(s16 - t9) * s34 - (
     #s25 * t21 + 3.0_dp * s26 + t41) * s26 - t38) * s34 - s16 * ((s16 *
     # t40 - t36) * s25 + t1) - (-t5 * (t27 + t39) - 3.0_dp * t37) * s26
     #) * s34 - t11 * t4 + t16 * t29 * t23) * s56 + t2 * t3 - t8 * (((t1
     #8 - s26) * s26 - t38) * s26 + t12) - t31 * (6.0_dp * t12 - (-s16 *
     # (t36 + t16) - t37) * s26 - t2) * s34 + t33 * (s25 * t35 + (t34 + 
     #s26) * s26) * s16 + s26 * (s16 - s26) * s34 * t8) * MBL1110(1) + t
     #13) + t26) * t42 + t10) - t25 ** 2.0_dp * MBL1111D2(0) * t42) + t2
     #0
      result = t1 / 2.0_dp + t14 * t42 * t7 * t15

           intHLs16s25s26s34s56x1311D4eps0 = result
       end function intHLs16s25s26s34s56x1311D4eps0

       function intHLs16s25s26s34s56x1311D4eps1()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1311D4eps1
           complex(dp) ::  t1

           complex(dp) :: result

      t1 = s25 + s26 - s34 + s56
      t1 = 0.1e1_dp / t1
      result = 0.1e1_dp / s16 * t1 * ((s25 + s26 + s56) * MBL1110(1) + 1
     #.0_dp) / 2.0_dp

           intHLs16s25s26s34s56x1311D4eps1 = result
       end function intHLs16s25s26s34s56x1311D4eps1

       function intHLs16s25s26s34s56x1312D6eps0()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1312D6eps0
           complex(dp) ::  t1,t10,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109
           complex(dp) ::  t11,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t12
           complex(dp) ::  t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t13,t130
           complex(dp) ::  t14,t15,t16,t17,t18,t19,t2,t20,t21,t22,t23,t24
           complex(dp) ::  t25,t26,t27,t28,t29,t3,t30,t31,t32,t33,t34,t35
           complex(dp) ::  t36,t37,t38,t39,t4,t40,t41,t42,t43,t44,t45,t46
           complex(dp) ::  t47,t48,t49,t5,t50,t51,t52,t53,t54,t55,t56,t57
           complex(dp) ::  t58,t59,t6,t60,t61,t62,t63,t64,t65,t66,t67,t68
           complex(dp) ::  t69,t7,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79
           complex(dp) ::  t8,t80,t81,t82,t83,t84,t85,t86,t87,t88,t89,t9
           complex(dp) ::  t90,t91,t92,t93,t94,t95,t96,t97,t98,t99

           complex(dp) :: result

      t1 = s16 - s25
      t2 = 2.0_dp
      t3 = s16 ** 2.0_dp
      t4 = t3 ** 2.0_dp
      t5 = s16 * t4
      t6 = s16 * t3
      t7 = t2 * t3
      t8 = s25 ** 2.0_dp
      t9 = t8 ** 2.0_dp
      t10 = t8 * t9
      t11 = s25 * t9
      t12 = s25 * t8
      t13 = t2 * t8
      t14 = 3.0_dp
      t15 = t14 * s25
      t16 = s25 + s26
      t17 = -37.0_dp * s16
      t18 = 11.0_dp * s25
      t19 = 8.0_dp * s16
      t20 = 7.0_dp * t3
      t21 = t2 * t6
      t22 = ((s25 - t19) * s25 - t20) * s25 - t21
      t23 = t3 * t12
      t24 = 7.0_dp * s16
      t25 = t6 * s25
      t26 = 15.0_dp
      t27 = 17.0_dp * s16
      t28 = 4.0_dp
      t29 = 5.0_dp * s16
      t30 = t28 * s25
      t31 = (t29 - t30) * s25
      t32 = 11.0_dp * t3
      t33 = -22.0_dp * t6
      t34 = (t31 + t32) * s25 + t33
      t35 = t14 * s16
      t36 = -s25 + t35
      t37 = s26 ** 2.0_dp
      t38 = s26 * t37
      t39 = s16 * s25
      t40 = t28 * s16
      t41 = 18.0_dp
      t42 = 6.0_dp
      t43 = t42 * s25
      t44 = t41 * t3
      t45 = (t43 + t29) * s25 + t44
      t46 = t2 * s16
      t47 = t15 + t46
      t48 = t3 * t8
      t49 = t14 * s26
      t50 = t39 * t2
      t51 = -44.0_dp * t3
      t52 = 21.0_dp * t6
      t53 = 13.0_dp * t4
      t54 = 20.0_dp
      t55 = s16 * t26
      t56 = 9.0_dp * s25
      t57 = 12.0_dp * s16
      t58 = 8.0_dp * s25
      t59 = 14.0_dp * s16
      t60 = t1 ** 2.0_dp
      t61 = t60 ** 2.0_dp
      t62 = t4 * s25
      t63 = t60 * (s16 * t1 + t13)
      t64 = t6 * t8
      t65 = (s16 - t15) * t60
      t66 = s34 ** 2.0_dp
      t67 = t66 ** 2.0_dp
      t68 = s34 * t67
      t69 = s34 * t66
      t70 = t2 * s26
      t71 = t3 * t9
      t72 = s16 * t12
      t73 = 9.0_dp * s26
      t74 = 11.0_dp * s16
      t75 = s16 * s26
      t76 = 9.0_dp * t6
      t77 = s16 * t8
      t78 = 7.0_dp * s25
      t79 = t3 * s25
      t33 = (t14 * s25 * (t36 * t66 - t77) + (t66 * (t47 - s34) + t22) *
     # s34 + t4 + t79 * (-t29 + t78)) * s56 - t26 * t48 * t1 - t14 * s16
     # * (-t65 * s26 + t4) + ((((-t24 - t30 - t49 + s34) * s34 + (t73 + 
     #t29) * s25 + t42 * (t8 + t75) + t44) * s34 + (s16 * (27.0_dp * s26
     # + t74) + (-t73 + t29 - t30) * s25) * s25 + t33) * s34 + ((t8 + t5
     #1) * s25 - t52) * s25 + t14 * (t22 * s26 + t72) + t53) * s34 + t39
     # * (-t12 * t42 + t76)
      t44 = t16 * ((s16 * (((t17 + t18) * s25 - 19.0_dp * t3) * s25 + 13
     #.0_dp * t6) + t22 * s26) * s26 - t23 + t25 * (-22.0_dp * s25 + t24
     #)) * s34
      t18 = (t33 * s56 - t14 * s16 * (-t65 * t37 + t11) + t26 * t64 * t1
     # - t42 * s16 * (t63 * s26 + t62) + ((t2 * t34 * s26 + t56 * t36 * 
     #t37 + (-((t58 + t59 + t49) * s26 + t50) * s34 + (t14 * t47 * s26 +
     # t2 * t45) * s26 + t39 * (t57 + t56)) * s34 - t55 * t12 - 24.0_dp 
     #* t25 + 17.0_dp * t48) * s34 + (t2 * ((((s25 + t35) * s25 + t51) *
     # s25 - t52) * s25 + t53) + t49 * t22) * s26 + t39 * (((-38.0_dp * 
     #s16 + t18) * s25 - 41.0_dp * t3) * s25 + t54 * t6)) * s34 + t70 * 
     #t68 + 9.0_dp * t71) * s56 + t11 * t3 + t37 * t68 + ((t38 * t36 * t
     #14 - 5.0_dp * t25) * s25 + t34 * t37 - t23 - t39 * ((s25 * t26 - t
     #27) * s25 + 24.0_dp * t3) * s26) * t66 + (t37 * (s26 * t47 + t45) 
     #+ t48 + t49 * t39 * (t40 + t15)) * t69 - s26 * ((s26 + t24 + t30) 
     #* s26 + t50) * t67 + t65 * s16 * t38 + t44
      t22 = t4 * t8
      t33 = t60 * s25
      t34 = 5.0_dp * s25
      t44 = t6 * t12
      t45 = s16 + s25
      t47 = -t2 * (s34 * t45 + t39) + t3 + t8 + t66
      t51 = (-s34 + t45 + s26) * s26 + (t70 + s16 + s25 - s34 + s56) * s
     #56 + t39
      t52 = t8 + t3
      t53 = t2 * t52
      t65 = t2 * s25
      t68 = s16 - t65
      t80 = 8.0_dp * t3
      t81 = t28 * t6
      t50 = t14 * t52 + t50
      t82 = -10.0_dp * s25
      t83 = -10.0_dp * s16
      t84 = t42 * t3
      t85 = (t83 - t34) * s25 + t84
      t86 = s16 + t65
      t87 = t14 * t86
      t88 = t2 * s25 * t45
      t89 = (-s26 + t40 - t15) * s26 + t88
      t90 = 8.0_dp * t6
      t91 = ((-40.0_dp * s16 + t43) * s25 - 22.0_dp * t3) * s25 + t90
      t92 = t41 * t77 * t45
      t93 = t42 * s25 * (t6 + t12)
      t94 = s25 * t52
      t95 = (t2 * t85 - t73 * t86) * s26 + t42 * t94 - 11.0_dp * t77
      t88 = (-t43 + t19 - t49) * s26 + t88
      t96 = s16 * s34
      t97 = t60 * ((t30 - t46) * s25 + t3)
      t98 = t48 * t52
      t99 = t14 * t68 * t60 * t37
      t100 = t84 * t12
      t101 = t14 * t50 * s26
      t102 = s16 * t42
      t103 = s25 * (t6 - t12)
      t104 = 9.0_dp * s16
      t54 = -t103 * t28 + (((t14 * t16 + s34 - t40) * s34 + s16 * (-t73 
     #+ t102) + (-s26 * t41 - t34 + t83) * s25) * s34 + ((s16 * t54 - t1
     #5) * s25 + t32) * s25 - t81 + t101) * s34 + t4 + t77 * (t82 + t104
     #) - t49 * t68 * t60
      t105 = t28 * t3
      t106 = s25 * ((t29 - t65) * s25 - t105) + (-t66 - t50) * s34 + t6 
     #+ t87 * t66
      t107 = t97 + t67
      t108 = s56 ** 2.0_dp
      t109 = s56 * t108
      t82 = ((-t87 * t38 + (s25 * (t42 * t52 - 11.0_dp * t39) + t85 * s2
     #6) * s26 + t48) * s34 - t16 * (t79 * (t82 + t46) + s26 * (-t50 * s
     #26 + ((-s16 * t41 + t43) * s25 - t80) * s25 + t81))) * s34
      t85 = t33 * (-t39 + t53)
      t87 = t107 * t37
      t110 = s26 + s56
      t36 = t36 * t86
      t111 = 16.0_dp * t79
      t112 = t42 * (t1 * t8 + t6)
      t113 = -t111 + t112
      t114 = t14 * t3
      t115 = (t102 + s25) * s25 + t114
      t116 = (t29 + t65) * s25
      t117 = t116 - t105
      t118 = t8 - t75
      t119 = 12.0_dp * s25
      t120 = 26.0_dp * t3
      t121 = t6 * t42
      t122 = t36 * t49
      t123 = t49 * t115
      t124 = t5 * s25
      t125 = 16.0_dp * t3
      t126 = 13.0_dp * t3
      t127 = t4 * t42
      t13 = (t28 * t39 * t52 + ((-t115 + t96) * s34 + t36 * t1) * s34 - 
     #t4 - t9 - t48 * t42) * s56 + s25 * ((((t104 - t65) * s25 - t125) *
     # s25 + 14.0_dp * t6) * s25 - t127) - (((-(t34 - t40 + t49) * s16 -
     # t13 - t96) * s34 + t111 - t112 + t123) * s34 + t1 * (((t43 - t55)
     # * s25 - t126) * s25 - t122 + t81)) * s34 + t5 - t49 * t61
      t40 = t5 * t8
      t111 = t1 * t16
      t36 = t111 * (t2 * (t45 * (-7.0_dp * t39 + t53) * s26 + t25) - t77
     # * (t29 + t15) - t36 * t37) * s34
      t10 = s16 * t10 + (t13 * s56 + t2 * (t68 * t61 * s26 + t124) - (((
     #-t2 * t117 * s26 - (t70 * s16 - t8) * s34 - t12 * t28 + t79 * t2 -
     # t35 * t37 + 5.0_dp * t77) * s34 - (t28 * t113 / 2.0_dp - t123) * 
     #s26 + t42 * s25 * ((t8 + t3) * s25 - t6) - 16.0_dp * t72) * s34 + 
     #t1 * (s25 * (((-13.0_dp * s16 + t30) * s25 - t26 * t3) * s25 + t12
     #1) + (-t122 + ((-30.0_dp * s16 + t119) * s25 - t120) * s25 + t90) 
     #* s26)) * s34 - t10 + t77 * (((t43 - t59) * s25 + t125) * s25 - t7
     #6) - t14 * t61 * t37) * s56 + t66 * ((-s26 * t118 * s34 + (s25 * (
     #-t31 - t7) + (t117 + t75) * s26) * s26 - t72) * s34 + s25 * (t77 *
     # (t102 + t15) + t70 * (-t14 * (-t6 + t94) + 8.0_dp * t77)) + t37 *
     # (-s26 * t115 + t113) + t64) + t40 - t36 + (s25 * (-s25 + t46) + (
     #t68 - s26) * s26) * s26 * t61
      t13 = t39 * t28 + t52
      t31 = (s16 - t78) * s25 + t114
      t36 = t56 * t37
      t45 = t28 * s26
      t53 = t14 * t8
      t61 = s26 * s34
      t55 = (((t55 + t65) * s25 - t105) * s25 - t90) * s25
      t66 = (t34 - t46) * s25 + t114
      t72 = t66 * s26
      t76 = ((t43 + t24) * s25 - 5.0_dp * t3) * s25 + t14 * t6
      t90 = t3 + 5.0_dp * t8
      t94 = t30 + t35
      t112 = (t58 + s16) * s25
      t113 = t112 + t105
      t115 = t78 + t35
      t117 = ((13.0_dp * s25 + t104) * s25 - t32) * s25
      t41 = t41 * s25
      t122 = (t6 + t117) * t14
      t123 = t45 * t66
      t125 = t94 * t37
      t116 = t1 * (t3 - t116)
      t128 = t1 * t13
      t129 = -29.0_dp
      t130 = t116 * t37
      t26 = ((s25 * ((t35 + t65) * s25 - t84) + ((t94 - s34) * s34 - t66
     #) * s34 + t6) * s56 + t28 * (t116 * s26 - t103) - ((-16.0_dp * s25
     # * s26 + (t78 + t45 + t35 - s34) * s34 - t26 * t8 - t114 - 12.0_dp
     # * t75) * s34 + t117 + t123 + t6) * s34 - 12.0_dp * t77 * t1) * s5
     #6 + s25 * (-12.0_dp * t128 * s26 + t55) + t42 * (-t62 + t130) + ((
     #((s16 + t49 + t65) * s34 - (s26 * t42 + 21.0_dp * s25 + t104) * s2
     #6 - t105 - t112) * s34 + s25 * ((t119 + t59) * s25 - 10.0_dp * t3)
     # + t42 * (t6 + t125) + t73 * t90) * s34 + (((s16 * t129 - t58) * s
     #25 + t120) * s25 + t26 * t6) * s25 + (-t72 * t42 - t122) * s26 - t
     #28 * t4) * s34 + t5
      t11 = t26 * s56 + t2 * (t1 * ((((-t27 - t65) * s25 - t126) * s25 -
     # 5.0_dp * t6) * s25 + t4) * s26 + t124) + t28 * (t116 * t38 + t71)
     # + (((t61 * (t30 + t49 + t46) - t2 * (s26 * t113 + t79) - t37 * (t
     #14 * t115 + t45) - t77 * t42) * s34 + t28 * s26 * (t125 + t76) + t
     #39 * ((t83 + t41) * s25 + t84) + 9.0_dp * t90 * t37) * s34 + ((-t1
     #22 - t123) * s26 - t1 * (((-t58 + t17) * s25 - t32) * s25 + t81) *
     # t2) * s26 + t39 * (((-t41 + t19) * s25 + 28.0_dp * t3) * s25 - t1
     #21)) * s34 - 12.0_dp * s25 * (t128 * t37 + t62) + t102 * t11
      t17 = s25 + s26 + s56
      t19 = 0.1e1_dp / t47 ** 2.0_dp
      t26 = 0.1e1_dp / t51 ** 2.0_dp
      t5 = s34 * (s25 * (-t128 * t28 * t38 + t64 * (-t102 + t15)) + t2 *
     # t39 * (s16 * t9 + ((t8 * (t15 + t46) - t121) * s25 + t4) * s26) +
     # (((-s26 * (t39 * (t43 + t46) + ((t115 + s26) * s26 + t113) * s26)
     # + t37 * (t65 + s16 + s26) * s34) * s34 + t2 * (t37 * t76 + t23 + 
     #t39 * ((-t29 + t56) * s25 + t114) * s26) + t37 ** 2.0_dp * t94 + t
     #64 + t14 * t90 * t38) * s34 - t16 * (t2 * s16 * ((((-t102 + t56) *
     # s25 - t80) * s25 + t21) * s26 + t25) + t37 * (t72 + ((t58 + t74) 
     #* s25 - 14.0_dp * t3) * s25 + t6) + t48 * (-t57 + t30))) * s34 + t
     #11 * s56 + t37 * (t130 + (t55 - t127) * s25 + t5) + t40) * MBL1110
     #(1) - t14 * s16 * (t63 * t37 + t22 + t33 * (-s25 * t1 + t7) * s26)
     # + t18 + t44 * (-t34 + t24) + t110 * (-t28 * t23 * t52 + t121 * t9
     # + t10) * MBL1011(1)
      t1 = t26 * ((s25 * s34 * (s34 * (t61 * (-(s25 - s26) * s34 + t118 
     #* t14 - t2 * s25 * (s16 - s26)) + t14 * s25 * (-t38 + t79) - s26 *
     # (s25 * ((s16 + t15) * s25 - t20) - t31 * s26)) + (((-t73 * t1 + t
     #28 * t8 - t114) * s25 - ((t35 - t65 - s34) * s34 - (-t73 + s16 - t
     #78) * s25 - t114) * s34 - t6) * s56 + (((-s25 + t70) * s34 + (t45 
     #+ t15 - t46) * s25 - t75 * t42) * s34 - t12 * t14 + t70 * t31 - t3
     #6 - t77 + 7.0_dp * t79) * s34 + t9 + t39 * (-t105 + t53) - t36 * t
     #1 - t70 * t1 * t86 ** 2.0_dp) * s56 - t111 * (s26 * t13 + t14 * s2
     #5 * (t3 + t37))) - t53 * (s16 - s25 + s34) * s34 * t109) * MBL1101
     #(0) * t19 + t110 * t17 ** 2.0_dp * MBL1111D2(0))
      result = t19 * t26 * ((t2 * t4 * t12 - s16 * (t82 + t98) + ((t99 -
     # t100) * s16 - t2 * s16 * ((t9 + t4) * s25 + (t97 + t67) * s26) + 
     #t96 * ((s34 * t88 - t95) * s34 + (-t49 * t50 + t91) * s26 - t92 + 
     #t93) + 5.0_dp * t98) * s56 - t75 * (-t69 * t89 + t85) - t87 * s16 
     #+ s16 * t106 * t109 - s16 * t54 * t108 + s16 * t68 * t60 * t38) * 
     #I300s16(1) + (-t44 * t2 + (-t69 * t89 + t85) * s26 + ((-t106 * s56
     # + t54) * s56 + t2 * ((t9 + t4) * s25 + s26 * t107) + t100 - t99 -
     # 5.0_dp * t77 * t52 + ((-s34 * t88 + t95) * s34 + t92 - t93 + s26 
     #* (-t91 + t101)) * s34) * s56 + t22 + t82 + t87 + t71 - t68 * t60 
     #* t38) * MBL1001(0) + t5 / s16) / 4.0_dp + t1 / 2.0_dp

           intHLs16s25s26s34s56x1312D6eps0 = result
       end function intHLs16s25s26s34s56x1312D6eps0

       function intHLs16s25s26s34s56x1312D6eps1()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1312D6eps1

           complex(dp) :: result

      result = 0.0_dp

           intHLs16s25s26s34s56x1312D6eps1 = result
       end function intHLs16s25s26s34s56x1312D6eps1

       function intHLs16s25s26s34s56x1321D6eps0()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1321D6eps0
           complex(dp) ::  t1,t10,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109
           complex(dp) ::  t11,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t12
           complex(dp) ::  t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t13,t130
           complex(dp) ::  t131,t132,t133,t134,t135,t136,t137,t14,t15,t16,t17,t18
           complex(dp) ::  t19,t2,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29
           complex(dp) ::  t3,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t4
           complex(dp) ::  t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t5,t50
           complex(dp) ::  t51,t52,t53,t54,t55,t56,t57,t58,t59,t6,t60,t61
           complex(dp) ::  t62,t63,t64,t65,t66,t67,t68,t69,t7,t70,t71,t72
           complex(dp) ::  t73,t74,t75,t76,t77,t78,t79,t8,t80,t81,t82,t83
           complex(dp) ::  t84,t85,t86,t87,t88,t89,t9,t90,t91,t92,t93,t94
           complex(dp) ::  t95,t96,t97,t98,t99

           complex(dp) :: result

      t1 = s16 - s25
      t2 = s25 * t1
      t3 = 3.0_dp
      t4 = 2.0_dp
      t5 = s16 ** 2.0_dp
      t6 = t5 ** 2.0_dp
      t7 = s16 * t5
      t8 = t3 * t5
      t9 = t4 * t2
      t10 = -t9 + t8
      t11 = t4 * s16
      t12 = t11 + s25
      t13 = s16 + s25
      t14 = s25 * t13
      t15 = 9.0_dp * t5
      t16 = -t15 + t14
      t17 = 6.0_dp
      t18 = s26 ** 2.0_dp
      t19 = t18 ** 2.0_dp
      t20 = s26 * t18
      t21 = s25 ** 2.0_dp
      t22 = t21 ** 2.0_dp
      t23 = s25 * t22
      t24 = s25 * t21
      t25 = 4.0_dp * t13
      t26 = t17 * s16
      t27 = t26 * s25
      t28 = 9.0_dp * s16
      t29 = t4 * s26
      t30 = 4.0_dp * s25
      t31 = t5 * s25
      t32 = t13 * t18
      t33 = 12.0_dp
      t34 = t26 * t21
      t35 = s16 * s25
      t36 = s34 ** 2.0_dp
      t37 = t36 ** 2.0_dp
      t38 = s34 * t36
      t39 = s26 * t1
      t40 = t1 ** 2.0_dp
      t41 = t31 * t1
      t42 = t33 * s26
      t43 = t17 * s26
      t44 = t4 * s25
      t45 = 4.0_dp * s16
      t46 = 5.0_dp * t5
      t47 = t43 * t40
      t48 = s34 * t13
      t49 = -t4 * (t35 + t48) + t21 + t5 + t36
      t50 = s26 * s34
      t51 = t18 * t38
      t52 = t5 * t21
      t53 = s56 ** 2.0_dp
      t54 = t40 * t20
      t55 = t31 * t39
      t56 = (s34 - t13 - s26) * s26 + (-t29 - s16 - s25 + s34 - s56) * s
     #56 - t35
      t57 = t3 * s16
      t58 = t57 - s25
      t59 = t21 - t5
      t60 = t11 * s25
      t61 = t60 - t59
      t62 = t3 * s25
      t63 = (-t62 + t26) * s25 + t5
      t64 = t11 - s25
      t65 = s25 + s26
      t66 = t21 + t5
      t67 = 8.0_dp * s16
      t68 = 5.0_dp * s25
      t69 = t68 + t45
      t70 = t69 * s26
      t71 = t31 * t17
      t72 = 19.0_dp * s25
      t73 = 32.0_dp * s16
      t74 = t4 * t7
      t75 = -33.0_dp
      t76 = 4.0_dp * t5
      t77 = 8.0_dp * s26
      t78 = t3 * t18
      t79 = t1 * t63
      t80 = s16 * t21
      t81 = t17 * s25
      t82 = t18 * t64
      t83 = 13.0_dp * s25
      t84 = 11.0_dp * s25
      t85 = t5 * t33
      t86 = t5 + t36
      t87 = 5.0_dp * s34
      t88 = t65 ** 2.0_dp
      t63 = (((-t57 - t87 + s25) * s25 + t4 * t86 - t45 * s34) * s56 - s
     #25 * ((-t62 + t28) * s25 - t46) + ((t57 + t77 + t84 - s34) * s34 -
     # 16.0_dp * s16 * t65 - (20.0_dp * s26 + t83) * s25 - t8) * s34 + t
     #7 + 4.0_dp * t39 * t64) * s56 + t3 * (s25 * (s25 * t66 + t7) + t36
     # * ((s16 - s25 - s26) * s34 - (-t57 - t84) * s26 + 4.0_dp * t14 + 
     #4.0_dp * t18 - t5) + t39 * t63) - s34 * (t38 + ((20.0_dp * s16 + t
     #84) * s25 + t85) * s25 + (t43 * t69 + t13 * (t57 + t83) * t3) * s2
     #6 - t7) + t82 * t17 * t1 - t28 * t24
      t69 = t18 * t37
      t89 = t5 * t24
      t90 = t62 * t1
      t91 = s25 + s26 - s34 + s56
      t92 = t2 - t5
      t93 = -t8 + t14
      t94 = t57 + s25
      t95 = t43 * t13
      t96 = t3 * s26
      t97 = t1 * t92
      t98 = s16 - s26
      t99 = t7 * t21
      t93 = t55 * t4 + (t18 * t93 + t4 * s26 * (-t12 * t35 + t32) - t52)
     # * s34 + ((-s56 * t49 - t4 * t35 * t1 + ((t3 * t98 + s25 - s34) * 
     #s34 + t14 - t8 + t95) * s34 - t24 + t7 - t96 * t40) * s56 - t4 * (
     #(t97 + t38) * s26 - t41) + s34 * (((-t96 + t26) * s26 + t4 * s25 *
     # (s16 + s26)) * s34 + (t4 * t93 + t95) * s26 - t35 * (t45 + t44)) 
     #- t78 * t40) * s56 - t51 - t54 - t89 + t99 - t97 * t18 + s26 * (t6
     #0 - s26 * (s26 - t94)) * t36
      t97 = 5.0_dp * s16
      t100 = t3 * t7
      t101 = (t57 - t44) * s25 + t5
      t102 = (-t81 + t28) * s25 + t5
      t103 = t57 - t44
      t104 = 7.0_dp
      t105 = t104 * s25
      t106 = t105 + t26
      t107 = t106 * s26
      t108 = 11.0_dp * s16
      t109 = (t108 + t81) * s25
      t110 = 8.0_dp * s25
      t111 = 8.0_dp * t5
      t112 = 16.0_dp * s25
      t113 = t1 * t102
      t114 = 4.0_dp * t103
      t115 = 5.0_dp * t14
      t116 = 22.0_dp * s16
      t117 = 14.0_dp * t5
      t118 = -48.0_dp
      t119 = 15.0_dp * s16
      t120 = 18.0_dp * s26
      t121 = 16.0_dp * s16
      t122 = t7 * s25
      t123 = t103 * t18
      t124 = t26 * s34
      t86 = t3 * t86
      t102 = (((s34 * t104 - t44 + t97) * s25 + t124 - t86) * s56 + s25 
     #* ((-t81 + t119) * s25 - t111) + t38 * t4 + 4.0_dp * ((t105 + t26)
     # * s26 + t115 + t5) * s34 - 4.0_dp * t39 * t103 - t36 * (t97 + t42
     # + t112) - t7) * s56 - t17 * (t1 * t123 + t21 * t66) - t3 * (t102 
     #* t39 + t122) + (((t43 + t68 - t57 + s34) * s34 + (-t72 - t121) * 
     #s25 + (s25 * t118 - t119 - t120) * s26 + t8) * s34 + ((t116 + t72)
     # * s25 + t117) * s25 + (t33 * (t115 + t5) + t43 * t106) * s26 - t7
     #) * s34 + t119 * t24
      t105 = t7 * t24
      t9 = -t90 * t101 * t18 + (t38 * ((t68 - t57) * s26 - t4 * (t2 - t1
     #8)) + t21 * (t21 * (t97 - t44) - t100)) * s26 + (t102 * s56 - t3 *
     # (t113 * t18 + t99) + t4 * ((s26 * t38 + t65 * (((t45 + t62) * s25
     # + t111) * s25 + t4 * s26 * (t107 + (t28 + t110) * s25 + t8) - t7)
     #) * s34 - t23) + t36 * ((-s26 * t17 * t98 + 10.0_dp * s25 * s26 - 
     #t9) * s34 + s25 * (-t109 + t76) + ((-t3 * (t97 + t112) - t42) * s2
     #6 + t13 * (t57 - t72) * t4) * s26) - t81 * t39 * t101 - t114 * t1 
     #* t20 + t97 * t22) * s56 + t22 * t5 + t69 - t1 * t103 * t19 - t113
     # * t20 - t65 * (t20 * t3 - (-(t97 + t83) * s26 - t109 + t8) * s26 
     #- t31) * t36 + t88 * (t71 - (-t107 - (t67 + t81) * s25 - t76) * s2
     #6 - t7) * s34 - t105
      t83 = s26 + s56
      t101 = t45 * s25 - t59
      t102 = 10.0_dp * s16
      t106 = s16 * t24
      t107 = s16 * t13
      t109 = 9.0_dp * s25
      t113 = s25 * (t109 - s16) + t8
      t115 = -t44 + s16
      t118 = -t4 * t115
      t119 = t44 + s16
      t104 = t104 * s16
      t125 = 8.0_dp * t7
      t95 = t95 * t64
      t126 = t17 * t7
      t127 = s25 * (t67 - t84) + t46
      t128 = t17 * t21
      t129 = t40 * (s25 * (-t62 + t102) + t5)
      t130 = -t68 + t104
      t131 = (s25 * t130 + t111) * s25
      t132 = t40 * t58
      t133 = 18.0_dp * s25
      t84 = s56 * (s56 * (-t3 * (t13 * t64 * s34 - t7) + t103 * t36 - t2
     #4 - t35 * t130) + s25 * (s25 * (s25 * (-t62 + t121) - 22.0_dp * t5
     #) + t125) + t36 * (s34 * t118 + s16 * (t97 + t42) + 8.0_dp * s25 *
     # t98 - 11.0_dp * t21) + t6 - t4 * (t95 + t131 + t74) * s34 + 4.0_d
     #p * t132 * s26) + t17 * (t132 * t18 + t99) + t3 * (s25 * (-t22 + t
     #6) + t129 * s26) + (s34 * (s34 * (-s34 * t119 + s16 * (-t43 + t57)
     # + s25 * (t42 + t109 - s16)) + s25 * (-15.0_dp * t21 + t117) + t17
     # * (t123 + t80) - t3 * (-s26 * t127 + t7)) + s25 * (s25 * (s25 * (
     #t84 - t116) - t85) - 14.0_dp * t7) + s26 * (-t17 * (t131 + t74) - 
     #t120 * t13 * t64) + t6) * s34 + t106 * (-24.0_dp * s16 + t133)
      t85 = t65 * (-t57 * t24 + s26 * (s26 * (-s26 * t103 - 5.0_dp * t10
     #7 + 9.0_dp * t21) + s25 * (s25 * (t81 - s16) - t15) + t100) + t31 
     #* t13) * t36
      t98 = t88 * (-t3 * (t82 * t13 + t106) + t5 * (s25 * (t44 - t26) + 
     #t5) - (-4.0_dp * s25 * t59 + 4.0_dp * t11 * t21 + 4.0_dp * t7) * s
     #26) * s34
      t84 = t62 * t40 * t101 * t18 + s26 * (-t37 * (s26 * t119 + t21) + 
     #t21 * t40 * (s25 * (t26 - s25) + t8)) + s56 * (s56 * t84 + t3 * (t
     #129 * t18 + t21 * t6) + t4 * s34 * (t36 * (t24 * t4 + s26 * (-s26 
     #* t3 * t115 + t113) + t35 * t13) + t65 * (s25 * (t21 * (t44 - t104
     #) - t125) + s26 * (-t3 * (s25 * (s25 * (t97 - t62) + t76) + t74) -
     # t95) + t6)) - t21 * t22 - t36 * (t36 * (s26 * (t30 + t11) + t21) 
     #- t21 * (t111 - t128) - 4.0_dp * t35 * t59 - s26 * (t18 * t114 + s
     #25 * (s25 * (s16 * t33 - 30.0_dp * s25) + 28.0_dp * t5) - t126) - 
     #t78 * t127) + t47 * s25 * t101 + 4.0_dp * t54 * t58 + s16 * t22 * 
     #(-t102 + t110)) + t38 * (t44 * (t21 * t4 + t107) * s26 - t18 * (-s
     #26 * t118 - t113) - t106) - t85 + t132 * t19 + t129 * t20 + t98 + 
     #t106 * (s25 * (-t2 - t5) + t7)
      t85 = s26 * t13
      t95 = s25 * (t44 + t108) - 51.0_dp * t5
      t97 = 15.0_dp * t6
      t98 = s25 * (s25 * t95 + t100) + t97
      t100 = -49.0_dp * t5
      t101 = 21.0_dp * t7
      t102 = 5.0_dp * t6
      t57 = s25 * (s25 * (s25 * (-t57 + t81) + t100) + t101) + t102
      t103 = 10.0_dp * t5
      t107 = s25 * (-21.0_dp * s16 + t81) + t103
      t109 = 14.0_dp * s16
      t111 = 26.0_dp * t5
      t68 = t68 + t26
      t113 = 15.0_dp * s25
      t114 = 15.0_dp * s26
      t94 = t94 * (-t110 + s16)
      t118 = 47.0_dp
      t119 = s25 * (192.0_dp * s16 + 81.0_dp * s25) + 78.0_dp * t5
      t120 = 88.0_dp * s16
      t121 = 72.0_dp * s25
      t123 = t121 + t120
      t125 = 9.0_dp * s26
      t127 = t25 * t107
      t129 = 5.0_dp * t58 * t115
      t130 = 19.0_dp * t7
      t131 = (s25 * (46.0_dp * s16 + 21.0_dp * s25) + 66.0_dp * t5) * s2
     #5
      t132 = 23.0_dp * s16
      t134 = s25 * (t132 + t133) + t111
      t135 = 10.0_dp * t68
      t115 = t58 * t115
      t136 = 10.0_dp * t115 * s26
      t107 = t13 * t107
      t137 = 30.0_dp * s26
      t86 = s56 * (s56 * (s25 * (t44 - t87 - t104) - t124 + t86) + s16 *
     # (s25 * (-t108 - t113) + t103) + t17 * (t24 - t38) - s34 * (-s34 *
     # (t114 + t116 + t133) + s16 * (26.0_dp * s16 + t137) + s25 * (25.0
     #_dp * s26 + t132 + t133)) + t129 * s26) + s25 * (s25 * (t100 + t12
     #8) + t101) + s26 * (t127 + t136) - t3 * (t38 * (5.0_dp * t13 + t77
     # - s34) + t106) + t102 - s34 * (-s34 * (s25 * (64.0_dp * s16 + 27.
     #0_dp * s25) + s26 * (t121 + t120 + t137) + t111) + t131 + t130 + s
     #26 * (t135 * s26 + 4.0_dp * t134))
      t86 = s56 * t86 + s25 * (t21 * t95 + t97) + t18 * (t107 * t17 + t1
     #36) + t3 * (s26 * t57 + t99) + (s34 * (s34 * (s34 * (t44 + t125 - 
     #s16) + s25 * (-t132 - t110) - 36.0_dp * t18 + t8 - 45.0_dp * t85) 
     #+ s16 * (s25 * (58.0_dp * s16 + 60.0_dp * s25) - t8) + s26 * (30.0
     #_dp * t18 + t119) + t33 * (t18 * t123 / 8.0_dp + t24)) + s25 * (s2
     #5 * (s25 * (-s16 * t118 - t110) - 58.0_dp * t5) - 52.0_dp * t7) + 
     #s26 * (s26 * (-t135 * s26 - t134 * t17) - t3 * (t131 + t130)) + t6
     #) * s34
      t87 = s25 * t98 - t37 * (-t96 - t44 + s16)
      t95 = t65 * (t19 * t3 - s26 * (s26 * (-s26 * (t113 + t116) - s25 *
     # (s25 * t33 + 42.0_dp * s16) - t111) - s16 * (s25 * (t73 + t133) -
     # t8)) - t122 + t52 * t4) * t36
      t6 = t88 * (t31 * (t30 + t109) + s26 * (s16 * (s25 * (t109 + t133)
     # + 19.0_dp * t5) + s26 * (s26 * t68 + s25 * (t108 + t110) + t111))
     # - t6) * s34
      t96 = t89 * (s25 * (t44 - t28) + t46)
      t99 = -s16 + s34
      t100 = 0.1e1_dp / t49
      t101 = 0.1e1_dp / s16
      t56 = 0.1e1_dp / t56 ** 2.0_dp
      t91 = 0.1e1_dp / t91 ** 2.0_dp
      t6 = s34 * (t96 + s26 * (-t38 * (s26 * (t114 * t13 - t94) + t17 * 
     #(t80 + t20) - t31 * t4) + t80 * t58 * (s25 * (-t45 - t81) + t46)) 
     #+ s56 * (s56 * t86 - 17.0_dp * t105 + t97 * t21 - t117 * t22 + t78
     # * t57 + t26 * t23 + t129 * t19 + t127 * t20 + t44 * t98 * s26 - (
     #s34 * (s34 * (t50 * (-t30 + t11 - t125) - t18 * (-24.0_dp * s26 - 
     #45.0_dp * t13) - t4 * (t94 * s26 + t31) + t34) - 15.0_dp * t19 - (
     #s26 * (s26 * t123 + t119) + s25 * (s25 * (120.0_dp * s16 + 24.0_dp
     # * s25) + 116.0_dp * t5) - t126) * s26 - t35 * (s25 * (34.0_dp * s
     #16 + t133) - t76)) + t65 * (s16 * (s25 * (s25 * (t116 + t133) + t1
     #18 * t5) - t74) + s26 * (s26 * (s25 * (s25 * t118 + 62.0_dp * s16)
     # + 104.0_dp * t5) + s25 * (s25 * (76.0_dp * s16 + t112) + 94.0_dp 
     #* t5) + 57.0_dp * t7) + 5.0_dp * t68 * t20)) * s34) + t18 * t87 + 
     #t20 * t57 + t115 * s26 * t19 + t107 * t19 + t95 - t6)
      t6 = t9 * MBL1001(0) - (t6 * MBL1110(1) + t83 * t84 * MBL1011(1)) 
     #* t101 - t90 * t61 * t18 - (-t38 * (-t3 * t39 + t18 + t21 - t60) +
     # t21 * t1 * t13 * t58) * s26 - (t63 * s56 - s26 * (-t1 * (t61 * t8
     #1 + 4.0_dp * t82) + t37 * t4) + t3 * (t18 * t79 - t59 * t80) - (((
     #(-t11 + s25) * s25 - t17 * t39 + t78) * s34 + s25 * ((-t62 - t28) 
     #* s25 + t76) + ((s25 * t75 - t28 - t77) * s26 + t17 * (-4.0_dp * t
     #14 + t5)) * s26) * s34 + t65 * (((t62 + t67) * s25 + 15.0_dp * t5)
     # * s25 + (4.0_dp * t70 + (t73 + t72) * s25 + t15) * s26 - t74)) * 
     #s34 + t24 * t59) * s56 + t69 - t1 * t64 * t19 - t79 * t20 - t65 * 
     #(t18 * (t3 * (t62 + s16) + t29) + (s25 * t28 + t3 * t59) * s26 - t
     #31) * t36 - t89 * t1 + t88 * (t71 - (-s25 * t67 - t3 * t66 - t70) 
     #* s26 - t7) * s34
      t2 = t101 * (t56 * (s25 * s34 * (s26 * (s26 * t92 + s34 * (t4 * t8
     #5 + t35 - t50)) + s56 * (-s56 * (-t4 * t48 - t2 + t36 + t5) - t4 *
     # s26 * (s25 * (-s16 + s25) + t36 + t5) + t35 * (-s16 - s25 + s34) 
     #+ 4.0_dp * t50 * t13) - t52 - t35 * t13 * s26) * MBL1101(0) * t100
     # + t83 * (-t18 * t99 - t4 * (s56 * (s26 * t99 - t35) - t35 * s26) 
     #- t53 * t99 + t80) * MBL1111D2(0)) + s34 * (s25 + s26 + s56) * MBL
     #1110(0) * t91)
      result = t100 * t56 * ((t3 * (t1 * t52 - t51) - t4 * (s56 * t49 * 
     #t53 + t54) + (t50 * ((t30 - t29 + t28) * s26 + t27) + t18 * (s26 *
     # t25 + t16) - t27 * t12 * s26 - t8 * t21) * s34 + ((s25 * ((t45 - 
     #t44) * s25 - t46) + t3 * (t7 - t38) + s34 * ((-t43 + t30 + t28) * 
     #s34 + t42 * t13 + t14 - t15) - t47) * s56 + (t4 * t16 * s26 - t33 
     #* (-t32 + t31) - t34) * s34 - t17 * ((s26 * t40 + t38) * s26 - t41
     #) + t4 * (((t30 + t28) * s26 + t3 * (t35 - t18)) * t36 + t39 * t10
     #)) * s56 + t1 * t10 * t18 + t55 * t17) * t101 + t93 * I300s16(1) +
     # t6 * t91) / 4.0_dp - t2 / 2.0_dp

           intHLs16s25s26s34s56x1321D6eps0 = result
       end function intHLs16s25s26s34s56x1321D6eps0

       function intHLs16s25s26s34s56x1321D6eps1()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1321D6eps1
           complex(dp) ::  t1,t2

           complex(dp) :: result

      t1 = s25 + s26 - s34 + s56
      t2 = -0.1e1_dp / 0.4e1_dp
      t1 = 0.1e1_dp / t1 ** 2.0_dp
      result = 0.1e1_dp / s16 * t1 * (t2 * (s25 + s26 + s34 + s56) - s34
     # * (s25 + s26 + s56) * MBL1110(1) / 2.0_dp)

           intHLs16s25s26s34s56x1321D6eps1 = result
       end function intHLs16s25s26s34s56x1321D6eps1

       function intHLs16s25s26s34s56x1411D6eps0()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1411D6eps0
           complex(dp) ::  t1,t10,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109
           complex(dp) ::  t11,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t12
           complex(dp) ::  t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t13,t130
           complex(dp) ::  t131,t132,t133,t134,t135,t136,t137,t138,t139,t14,t140,t141
           complex(dp) ::  t142,t143,t144,t145,t146,t147,t148,t149,t15,t150,t151,t152
           complex(dp) ::  t153,t154,t155,t156,t157,t158,t159,t16,t160,t161,t162,t163
           complex(dp) ::  t164,t165,t166,t167,t168,t169,t17,t170,t171,t172,t173,t174
           complex(dp) ::  t175,t176,t177,t178,t179,t18,t180,t181,t182,t183,t184,t185
           complex(dp) ::  t186,t187,t188,t189,t19,t190,t191,t192,t193,t194,t195,t196
           complex(dp) ::  t197,t198,t199,t2,t20,t200,t201,t202,t203,t204,t205,t206
           complex(dp) ::  t207,t208,t209,t21,t210,t211,t212,t213,t214,t215,t216,t217
           complex(dp) ::  t218,t219,t22,t220,t221,t222,t223,t224,t225,t226,t227,t228
           complex(dp) ::  t229,t23,t230,t231,t232,t233,t234,t235,t236,t237,t238,t239
           complex(dp) ::  t24,t240,t241,t242,t243,t244,t245,t246,t247,t248,t249,t25
           complex(dp) ::  t250,t251,t252,t253,t254,t255,t256,t257,t258,t259,t26,t260
           complex(dp) ::  t261,t262,t263,t264,t265,t266,t267,t268,t269,t27,t270,t271
           complex(dp) ::  t272,t273,t274,t275,t276,t277,t278,t279,t28,t280,t281,t282
           complex(dp) ::  t283,t284,t29,t3,t30,t31,t32,t33,t34,t35,t36,t37
           complex(dp) ::  t38,t39,t4,t40,t41,t42,t43,t44,t45,t46,t47,t48
           complex(dp) ::  t49,t5,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59
           complex(dp) ::  t6,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t7
           complex(dp) ::  t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t8,t80
           complex(dp) ::  t81,t82,t83,t84,t85,t86,t87,t88,t89,t9,t90,t91
           complex(dp) ::  t92,t93,t94,t95,t96,t97,t98,t99

           complex(dp) :: result

      t1 = 3.0_dp
      t2 = t1 * s25
      t3 = -t2 + s16
      t4 = s16 - s25
      t5 = s16 + s25
      t6 = 6.0_dp
      t7 = 5.0_dp
      t8 = 4.0_dp
      t9 = s16 ** 2.0_dp
      t10 = t9 ** 2.0_dp
      t11 = t9 * t10
      t12 = s16 * t10
      t13 = s16 * t9
      t14 = t8 * s25
      t15 = t7 * s16
      t16 = t6 * t9
      t17 = s25 ** 2.0_dp
      t18 = t17 ** 2.0_dp
      t19 = s25 * t17
      t20 = t17 * t18
      t21 = s25 * t18
      t22 = t12 - t21
      t23 = 33.0_dp
      t24 = 49.0_dp * s16
      t25 = t1 * t13
      t26 = s16 * s25
      t27 = t22 * t6 + t26 * (((-s25 * t23 - t24) * s25 - 41.0_dp * t9) 
     #* s25 + t25)
      t28 = t10 + t18
      t29 = 42.0_dp
      t30 = 30.0_dp
      t31 = 2.0_dp
      t32 = t29 * t9
      t33 = t30 * t10
      t34 = t31 * t11
      t35 = (((18.0_dp * t17 + t32) * s25 - 50.0_dp * t13) * s25 - t33) 
     #* t17
      t36 = 9.0_dp * t26 * t28
      t37 = 18.0_dp * s25
      t38 = t1 * s16
      t39 = 53.0_dp * t9
      t40 = 9.0_dp * t13
      t41 = t7 * t10
      t42 = 16.0_dp
      t43 = s25 * t6
      t44 = t9 * t42
      t45 = 7.0_dp * t13
      t46 = ((t38 - t43) * s25 - t44) * s25 + t45
      t47 = s25 + s26
      t48 = 37.0_dp
      t49 = 15.0_dp
      t50 = 12.0_dp
      t51 = s16 * t50
      t52 = t31 * t13
      t53 = t13 - t19
      t54 = 9.0_dp * s25
      t55 = 7.0_dp * s16
      t56 = 18.0_dp * t13
      t57 = 28.0_dp
      t58 = 25.0_dp
      t59 = t6 * s16
      t60 = 19.0_dp * t9
      t61 = ((s25 * t58 - t59) * s25 - t60) * s25 + t13 * t57
      t62 = s26 ** 2.0_dp
      t63 = t62 ** 2.0_dp
      t64 = s26 * t63
      t65 = s26 * t62
      t66 = t9 - t17
      t67 = t9 * t17
      t68 = t50 * t12
      t69 = 98.0_dp * s16
      t70 = 140.0_dp * s25
      t71 = t49 * t9
      t72 = t31 * t10
      t73 = 102.0_dp * s25
      t74 = 51.0_dp * t10
      t75 = 95.0_dp * s16
      t76 = 64.0_dp * t9
      t77 = 94.0_dp * t13
      t78 = 31.0_dp
      t79 = 39.0_dp * s25
      t80 = (s16 * t78 + t79) * s25 + t32
      t81 = t13 * t17
      t82 = t9 * s25
      t83 = t42 * t13
      t84 = 75.0_dp * s25
      t85 = 84.0_dp * t13
      t86 = (((188.0_dp * s16 + t84) * s25 + 136.0_dp * t9) * s25 + t85)
     # * s25 + t41
      t87 = ((363.0_dp * s16 + 290.0_dp * s25) * s25 + 273.0_dp * t9) * 
     #s25
      t88 = t77 + t87
      t89 = 84.0_dp * s25
      t90 = (87.0_dp * s16 + t89) * s25 + 58.0_dp * t9
      t91 = t57 * s16
      t92 = 27.0_dp * s25
      t93 = -t92 - t91
      t94 = t26 * (((152.0_dp * s16 + 160.0_dp * s25) * s25 + 89.0_dp * 
     #t9) * s25 + t83)
      t95 = 36.0_dp
      t96 = t42 * s16
      t97 = s25 * t95
      t98 = 79.0_dp * s16
      t99 = 24.0_dp * t9
      t100 = 150.0_dp * t9
      t101 = ((239.0_dp * s16 + 130.0_dp * s25) * s25 + t100) * s25
      t102 = 20.0_dp * t13
      t103 = t102 + t101
      t104 = (55.0_dp * s16 + 50.0_dp * s25) * s25 + t99
      t105 = 82.0_dp * s25
      t106 = 69.0_dp * s16
      t107 = t106 + t105
      t108 = t26 * ((100.0_dp * s25 + t98) * s25 + t99)
      t109 = 14.0_dp * s25
      t110 = t31 * s25
      t111 = s16 + t110
      t112 = 66.0_dp * s25
      t113 = (78.0_dp * s16 + t112) * s25
      t114 = 20.0_dp * t9
      t115 = t114 + t113
      t116 = 73.0_dp * s25
      t117 = 51.0_dp * s16
      t118 = t117 + t116
      t119 = 10.0_dp
      t120 = t119 * s16
      t121 = t26 * t8
      t122 = 27.0_dp * s16
      t123 = 19.0_dp * t13
      t124 = 11.0_dp * t10
      t125 = t31 * t12
      t126 = ((((-t122 - t37) * s25 - 69.0_dp * t9) * s25 - t123) * s25 
     #+ t124) * s25 + t125
      t127 = 187.0_dp * t13
      t128 = 54.0_dp * s16
      t129 = t18 * t5
      t130 = 39.0_dp * t9
      t131 = t7 * t65
      t132 = t13 * s25
      t133 = 51.0_dp * s25
      t134 = (((837.0_dp * s16 + 705.0_dp * s25) * s25 + 468.0_dp * t9) 
     #* s25 + 561.0_dp * t13) * s25 + 153.0_dp * t10
      t135 = 63.0_dp * s16
      t136 = 47.0_dp * t13
      t137 = ((86.0_dp * s25 + t135) * s25 + t39) * s25 + t136
      t138 = t12 * s25
      t139 = 8.0_dp
      t140 = 27.0_dp * t9
      t141 = t7 * s26
      t142 = t141 * t80
      t143 = t137 * t139
      t144 = s16 * t17
      t145 = t90 * t139
      t146 = t88 * t1
      t147 = 35.0_dp * s26
      t148 = t104 * t50
      t149 = 64.0_dp * s26
      t150 = 32.0_dp * s25
      t151 = t23 * s26
      t152 = s34 ** 2.0_dp
      t153 = t152 ** 2.0_dp
      t154 = s34 * t152
      t155 = t154 * t153
      t156 = t152 * t153
      t157 = s34 * t153
      t39 = t4 * ((((-t38 - t37) * s25 - t39) * s25 + t40) * s25 + t41)
      t158 = t39 * t62
      t128 = ((((-((20.0_dp * s16 + t97 + t151) * s26 + t121) * s34 + ((
     #t1 * t118 + t149) * s26 + t8 * t115 / 2.0_dp) * s26 + t26 * (t96 +
     # t150)) * s34 - (((t107 * t8 + t147) * s26 + t148) * s26 + t31 * t
     #103) * s26 - t108) * s34 + (((-t141 * t93 + t145) * s26 + t146) * 
     #s26 + t8 * t86) * s26 + t94) * s34 - t62 * ((t143 + t142) * s26 + 
     #t134) - t8 * ((((((162.0_dp * s16 + t133) * s25 + 95.0_dp * t9) * 
     #s25 + 54.0_dp * t13) * s25 + 57.0_dp * t10) * s25 + t12) * s26 + t
     #138) - t144 * (((142.0_dp * s16 + t70) * s25 + t140) * s25 + 75.0_
     #dp * t13)) * s34 + t47 * (64.0_dp * t129 * s16 + s26 * (((((227.0_
     #dp * s25 + t128) * s25 - t9) * s25 - 68.0_dp * t13) * s25 + 144.0_
     #dp * t10) * s26 + ((((228.0_dp * s16 + 76.0_dp * s25) * s25 + t76)
     # * s25 - t127) * s25 + 66.0_dp * t10) * s25 + 45.0_dp * t12) + t13
     #1 * t61 + t132 * ((-s25 * t48 - t128) * s25 + t130))
      t159 = t50 * s25
      t160 = 11.0_dp * s16
      t161 = 17.0_dp * t13
      t162 = 13.0_dp * t10
      t163 = t12 * t6
      t164 = t1 * t4
      t165 = 21.0_dp
      t166 = (94.0_dp * s16 + 101.0_dp * s25) * s25
      t167 = t165 * t9
      t168 = 85.0_dp * t13
      t169 = 70.0_dp * t10
      t170 = t49 * t12
      t171 = 44.0_dp * s25
      t172 = t50 * t9
      t173 = 18.0_dp * t10
      t174 = (((t38 + t171) * s25 - t172) * s25 + t40) * s25 + t173
      t175 = 112.0_dp * t13
      t176 = t10 * t6
      t177 = t12 * t29
      t178 = s26 * t119
      t179 = t178 * t61
      t180 = 108.0_dp * t13
      t181 = t178 * t80
      t182 = -t119 * t93 * s26
      t183 = 70.0_dp * s26
      t184 = t12 * t17
      t185 = t139 * t9
      t186 = t119 * t4 * t46
      t187 = s16 * t19
      t188 = t187 * (((t92 + t96) * s25 - t185) * s25 - 44.0_dp * t13)
      t189 = -88.0_dp
      t190 = 168.0_dp * s25
      t191 = 116.0_dp * s16
      t192 = t4 * t46
      t193 = s16 * t49
      t194 = t10 + t153
      t195 = 23.0_dp
      t196 = 9.0_dp * s16
      t197 = (-t196 + t43) * s25
      t198 = t195 * t13
      t199 = t47 ** 2.0_dp
      t105 = (s25 * ((t197 + t60) * s25 - t198) + 7.0_dp * t194 + ((s34 
     #* t93 + t80) * s34 - t61) * s34) * s56 + s25 * ((((-t193 + t37) * 
     #s25 + 50.0_dp * t9) * s25 - 62.0_dp * t13) * s25 + t10 * t8) - t15
     #7 * t42 + t7 * (s26 * t192 + t12) + ((((t147 + t106 + t105) * s34 
     #- s16 * (140.0_dp * s26 + t191) - (174.0_dp * s16 + 135.0_dp * s26
     # + t190) * s25) * s34 + ((126.0_dp * s16 + 172.0_dp * s25) * s25 +
     # 106.0_dp * t9) * s25 + t142 + t77) * s34 + (((s25 * t189 - t59) *
     # s25 + t99) * s25 - t56) * s25 - t10 * t95 - t141 * t61) * s34
      t87 = t105 * s56 + t34 + t35 + t36 + ((((-(t117 + t116 + t149) * s
     #34 + (220.0_dp * s16 + 200.0_dp * s25) * s25 + (276.0_dp * s16 + 3
     #28.0_dp * s25 + t183) * s26 + 96.0_dp * t9) * s34 - (t145 + t182) 
     #* s26 - t77 - t87) * s34 + (((279.0_dp * s16 + 235.0_dp * s25) * s
     #25 + 156.0_dp * t9) * s25 + t127) * s25 + (t181 + t143) * s26 + t7
     #4) * s34 + (((-t166 - t167) * s25 + t168) * s25 - t169) * s25 + (-
     #t139 * t174 - t179) * s26 - t170) * s34 + 11.0_dp * t156 + t39 * t
     #8 * s26 + t186 * t62
      t87 = t87 * s56 - t1 * (-s26 * t126 * t4 + t184) - t155 * t31 + t6
     # * ((t11 + t20) * s25 + t158) + ((((((t120 + t151 + t37) * s34 - (
     #153.0_dp * s16 + 219.0_dp * s25 + 96.0_dp * s26) * s26 - t113 - t1
     #14) * s34 + ((t107 * t6 + t183) * s26 + t148) * s26 + t101 + t102)
     # * s34 - (((376.0_dp * s16 + 150.0_dp * s25) * s25 + 272.0_dp * t9
     #) * s25 + 168.0_dp * t13) * s25 - ((t50 * t90 + t182) * s26 + t146
     #) * s26 - t10 * t119) * s34 + ((((324.0_dp * s16 + t73) * s25 + 19
     #0.0_dp * t9) * s25 + t180) * s25 + 114.0_dp * t10) * s25 + ((t137 
     #* t50 + t181) * s26 + t134) * s26 + t125) * s34 + s25 * (((((-146.
     #0_dp * s16 - 38.0_dp * s25) * s25 - t76) * s25 + t175) * s25 - t17
     #6) * s25 - t177) + ((-t174 * t50 - t179) * s26 - t1 * ((((t166 + t
     #167) * s25 - t168) * s25 + t169) * s25 + t170)) * s26) * s34 + t18
     #6 * t65 + t188
      t101 = t4 * s25
      t70 = t47 * (t50 * t81 * t5 + (s16 * ((((t69 + t70) * s25 + t71) *
     # s25 + 63.0_dp * t13) * s25 + t72) + ((t80 * s26 + ((133.0_dp * s2
     #5 + t75) * s25 + t76) * s25 + t77) * s26 + (((184.0_dp * s16 + t73
     #) * s25 + 92.0_dp * t9) * s25 + 93.0_dp * t13) * s25 + t74) * s26)
     # * s26 + t82 * (44.0_dp * t19 + t52)) * t152
      t56 = t199 * (-26.0_dp * t67 * t66 + s26 * (s16 * ((((64.0_dp * s2
     #5 + t51) * s25 - t48 * t9) * s25 - t52) * s25 + t10 * t49) + t61 *
     # t62) + t68 * s25 + t31 * (s16 * (t17 * (t54 + t55) + t56) - 19.0_
     #dp * s25 * t53) * t62) * s34
      t61 = t13 * t19
      t73 = -t155 * t62 + t53 * t61
      t74 = t9 * t18
      t27 = t31 * t73 - t6 * t74 * t53 + (t156 * ((11.0_dp * s26 + t37 +
     # t120) * s26 + t121) + t144 * t3 * t4 * t5 * ((t14 + t15) * s25 + 
     #t16)) * s26 + (t87 * s56 + t8 * s26 * (t158 - t155) - s34 * t128 +
     # t13 * t21 + t110 * t4 * t27 * s26 + t164 * t126 * t62 + t7 * t4 *
     # t46 * t63 + t144 * (((t17 * (t159 + t160) - t161) * s25 - t162) *
     # s25 + t163)) * s56 + (t35 + t36 + t34) * t65 + ((-(t42 * s26 * (t
     #111 * t26 + t65) + t62 * (s26 * t118 + t115) + t67 * t31) * s34 + 
     #t65 * (t104 * t8 + 7.0_dp * t62) + ((t107 * t62 + t103) * s26 + t1
     #08) * s26 + t67 * (t109 + t59)) * s34 - t31 * t62 * (t62 * t90 + t
     #86) - (t62 * (-t62 * t93 + t88) + t94) * s26 + t67 * ((-t96 - t97)
     # * s25 - t16)) * t154 + t101 * t27 * t62 + t70 - t56 + t192 * t64 
     #+ t39 * t63
      t34 = -t31 * (s34 * t5 + t26) + t9 + t152 + t17
      t35 = s25 + s26 - s34 + s56
      t36 = t31 * s26
      t39 = (s34 - t5 - s26) * s26 + (-s16 - s25 + s34 - t36 - s56) * s5
     #6 - t26
      t46 = t7 * t9
      t56 = t10 * t5
      t70 = 14.0_dp * s16
      t73 = t67 * (-t70 + t43)
      t77 = t1 * t129
      t80 = t82 * t31
      t53 = -t80 + t53
      t86 = t119 * s25
      t87 = t59 + t86
      t88 = t8 * t13
      t90 = 13.0_dp * s25
      t93 = 18.0_dp * s16
      t94 = t119 * t13
      t97 = (t17 * (t90 + t93) + t94) * s25 + t41
      t103 = t13 + t19
      t104 = t31 * t103 - t82
      t105 = t6 * s26
      t106 = t1 * t10
      t107 = 17.0_dp * s16
      t108 = 24.0_dp * s25
      t113 = 32.0_dp * s16
      t115 = 22.0_dp * s25
      t116 = ((t113 + t115) * s25 + t99) * s25
      t118 = t116 + t94
      t125 = t9 + t17
      t126 = t125 * t6 + t26 * t7
      t127 = t26 * ((t107 + t108) * s25 + t172)
      t128 = t31 * s16
      t129 = (t54 + t160) * s25
      t134 = t46 + t129
      t137 = t50 * s26
      t142 = t1 * s26
      t143 = 7.0_dp * s25
      t145 = t26 * t31
      t146 = 9.0_dp * s26
      t148 = t146 * t126
      t149 = t5 * s26
      t151 = t31 * t154
      t158 = -t38 + t110
      t109 = -(t77 + t73 + ((t152 * (t143 + t146 + t15 - s34) + t116 + t
     #148 + t94) * s34 - (13.0_dp * t19 + t94) * s25 - 18.0_dp * s26 * t
     #104 - 18.0_dp * t187 - t41) * s34 + t56 + t146 * t4 * t53 - t151 *
     # (18.0_dp * t149 + t129 + t46)) * s56 + t31 * ((s26 * t153 + s26 *
     # t97 + 9.0_dp * t104 * t62 - t81 + t26 * (t17 * t87 + t88)) * s34 
     #- t138 - t4 * ((((-t2 - t59) * s25 - t172) * s25 + t52) * s25 + t1
     #0) * s26) + t8 * (t154 * ((t146 * t5 + t134) * s26 + t26 * (t2 + t
     #128)) + t61) - t152 * (t152 * ((t120 + t146 + t109) * s26 + t145) 
     #+ (t118 * t31 + t148) * s26 + t127) - 9.0_dp * t4 * t53 * t62 - t1
     #44 * (t17 * (t38 + t43) - t45)
      t116 = t26 * t4
      t129 = t13 * t18
      t148 = s56 ** 2.0_dp
      t160 = t148 ** 2.0_dp
      t166 = s56 * t160
      t168 = s56 * t148
      t169 = t9 * t21
      t174 = 22.0_dp * s16
      t179 = t8 * s16
      t181 = t1 * t28 + t26 * ((s25 * t195 - t179) * s25 - t46)
      t182 = t7 * t13
      t135 = (((t92 + t135) * s25 - 38.0_dp * t9) * s25 + t182) * s25 + 
     #t106
      t183 = 17.0_dp * t9
      t186 = ((t92 + t193) * s25 - t183) * s25 + t182
      t188 = t57 * s25
      t189 = (t188 + t59) * s25
      t192 = 32.0_dp * t9
      t200 = t111 * t65
      t201 = s25 * (t10 - t200)
      t202 = t165 * s25
      t203 = t29 * t13
      t204 = t49 * s25
      t205 = t30 * t9
      t206 = 13.0_dp * s16
      t207 = 89.0_dp * s16
      t208 = ((s25 * t42 + t207) * s25 + 90.0_dp * t9) * s25
      t209 = t94 + t208
      t210 = s16 * t48
      t211 = (t210 + t204) * s25
      t212 = t99 + t211
      t213 = t2 + t15
      t91 = t26 * ((t91 + t204) * s25 + t16)
      t214 = t213 * t8
      t215 = 11.0_dp * s25
      t117 = (29.0_dp * s25 + t117) * s25
      t216 = t119 * t9
      t217 = t216 + t117
      t218 = t78 * s25
      t219 = t218 + t122
      t220 = 27.0_dp * t10
      t221 = s16 * t165
      t222 = 11.0_dp * t9
      t223 = 14.0_dp * t13
      t224 = ((-t159 + t38) * s25 + t60) * s25 + t223
      t225 = (-t54 + s16) * s25 + t71
      t226 = t10 * s25
      t227 = 34.0_dp
      t228 = t174 * t18
      t229 = t9 * t227
      t230 = s16 * t139
      t231 = 47.0_dp * t9
      t232 = 9.0_dp * t10
      t233 = t8 * t4 * t186
      t234 = t13 * t139
      t235 = ((108.0_dp * s16 + 312.0_dp * s25) * s25 + 192.0_dp * t9) *
     # s25 - 120.0_dp * t13
      t236 = 60.0_dp * s25
      t237 = t9 * t19
      t238 = t29 * t10
      t239 = t30 * s25
      t240 = t30 * s26
      t241 = 81.0_dp * s16
      t242 = 93.0_dp * s25
      t243 = 18.0_dp * s26
      t244 = 24.0_dp * t13
      t245 = t3 * t65
      t246 = s25 * s26
      t247 = t9 + t246
      t248 = 40.0_dp * s26
      t249 = 20.0_dp * s26
      t250 = s16 * t95
      t251 = t6 * t157
      t252 = 60.0_dp * s26
      t253 = s26 * t49
      t92 = (s25 * (((-t92 + t51) * s25 + t192) * s25 - 22.0_dp * t13 - 
     #t253 * t3 * t4) + t194 * t7 + s34 * (((-t253 + t128 - t37) * s25 +
     # t205) * s34 + (s16 * (t240 + t113) + (52.0_dp * s25 + t93 + t252)
     # * s25) * s25 - t102) - t214 * t154) * s56 + t1 * t12 + t17 * (((-
     #t92 - t250) * s25 + 101.0_dp * t9) * s25 - 43.0_dp * t13) + t31 * 
     #(-t154 * (s16 * t248 + t211 + 24.0_dp * t247) + t226) + ((t139 * t
     #225 * s26 + t152 * (t218 + t122 + t249) - t19 * t95 - t239 * t62 +
     # 9.0_dp * t144 + t203 + 57.0_dp * t82) * s34 + 0.2e1_dp / 0.3e1_dp
     # * t235 * s26 - t173 + 68.0_dp * t18 + t236 * t111 * t62 + 40.0_dp
     # * t67 + 74.0_dp * t187 - t132 * t42) * s34 + t233 * s26 - t101 * 
     #t30 * t3 * t62 - t251
      t22 = t92 * s56 + 9.0_dp * s25 * t22 + t1 * (s26 * t135 * t4 + t61
     #) + ((t50 * t225 * t62 + (((-t159 - t15 - t243 + s34) * s34 + (t24
     #2 + t240 + t241) * s26 + t117 + t216) * s34 - (24.0_dp * s26 * t21
     #3 + t212 * t6) * s26 - t208 - t94) * s34 + t146 * t224 - t165 * t1
     #8 - t187 * t49 - t239 * t65 + 84.0_dp * t132 + t41 + 57.0_dp * t67
     #) * s34 + 118.0_dp * s16 * t18 - t238 * s25 + t200 * t236 + t21 * 
     #t57 + t235 * t62 + t237 * t31 - t12 + 27.0_dp * t81 - t105 * ((((-
     #s25 * t227 - t210) * s25 - t114) * s25 + t234) * s25 + t232)) * s3
     #4 - t245 * t101 * t30 + t6 * t4 * t186 * t62 + t144 * (t17 * (-t23
     #6 + t241) - t244)
      t22 = t22 * s56 + t31 * ((t152 * (t152 * (s26 * s34 - (t159 + t146
     # + t15) * s26 - t26) - t139 * t213 * t65 - t91 - s26 * (s26 * t1 *
     # t212 + t209)) + t47 * (-t201 * t49 - (-(((74.0_dp * s25 + t221) *
     # s25 + t76) * s25 - 40.0_dp * t13) * s26 - (((90.0_dp * s16 + t188
     #) * s25 - t8 * t9) * s25 + t83) * s25 + t220) * s26 - t12 + t144 *
     # (t222 + t189))) * s34 + t169) + t152 * (t152 * (20.0_dp * t65 + s
     #26 * (s26 * t1 * t219 + t217 * t31) + t26 * (t230 + t115)) + t139 
     #* (t225 * t65 + t226) - t228 - t204 * t63 + 9.0_dp * t224 * t62 + 
     #t229 * t19 + 55.0_dp * t81 + t36 * ((((-t193 - t202) * s25 + 57.0_
     #dp * t9) * s25 + t85) * s25 + t41)) - t101 * t49 * t3 * t63 + t164
     # * t135 * t62 + t233 * t65 + t105 * t101 * t181 + t144 * (((-24.0_
     #dp * t17 + t231) * s25 - t13 * t227) * s25 + t232)
      t22 = t22 * s56 + t153 * ((t145 * (t179 + t215) + t131) * s26 + t6
     #2 * (s26 * t219 + t217) + t67) + t154 * (-t31 * s26 * (t212 * t62 
     #+ t91) + t62 * (-t214 * t62 - t209) - t67 * (t38 + t86)) + t156 * 
     #t62 - s26 * ((t105 + t159 + t15) * s26 + t145) * t157 + t47 * (-t1
     # * s25 * (-t10 + t63) + (s16 * (((s16 * t30 - t115) * s25 + t32) *
     # s25 + t182) + (((t128 - t204) * s25 + t205) * s26 + ((t55 - t202)
     # * s25 + t140) * s25 + t203) * s26) * s26 + t67 * (t206 + t14)) * 
     #t152 + t199 * ((18.0_dp * t187 - t36 * (-t116 * t50 - t19 * t57 + 
     #t40)) * s16 - t201 * t6 - t125 * t13 + t62 * ((t192 + t189) * s25 
     #- t102)) * s34 + (t62 * (s26 * t186 + t135) + t144 * (((t174 + t10
     #8) * s25 - t58 * t9) * s25 + t40)) * s26 * t4
      t76 = s26 + s56
      t83 = 18.0_dp * t9
      t85 = t6 * t13
      t91 = (((t2 + t113) * s25 - t83) * s25 - t85) * s25 + t232
      t92 = -t30 * t144 * t4 + (9.0_dp * t19 + t234) * s25 + t106
      t93 = ((t54 + t230) * s25 - t172) * s25 + t182
      t101 = t13 * t50
      t108 = t9 * t57
      t117 = t8 * t17
      t131 = -t158
      t135 = t30 * t13
      t164 = t65 * t131
      t173 = t36 * s16
      t186 = 24.0_dp * s16
      t188 = 14.0_dp * t9
      t189 = t38 + t143
      t201 = (t204 + t250) * s25
      t208 = ((t201 + t231) * s25 + 45.0_dp * t13) * s25 + t41
      t209 = t103 * t6 + t26 * (t120 + t143)
      t211 = (t186 + t202) * s25 + t114
      t212 = t1 * t9
      t213 = t211 * s26
      t214 = 51.0_dp * t9
      t201 = (t201 + t214) * s25
      t217 = t201 + t94
      t219 = 17.0_dp * s25
      t224 = (t210 + t219) * s25 + t140
      t225 = t14 + t15
      t185 = t26 * ((t215 + t174) * s25 + t185)
      t226 = t38 + s25
      t231 = t226 * t62
      t233 = ((t2 + t51) * s25 + t46) * s26
      t234 = 123.0_dp * t9
      t235 = t29 * s25
      t241 = 105.0_dp * s16
      t254 = 58.0_dp * t13
      t255 = 84.0_dp * t10
      t256 = t7 * t12
      t257 = ((s16 + t79) * s25 + t188) * s25 + t135
      t258 = t138 * t139
      t259 = t7 * s25
      t99 = t1 * ((((62.0_dp * s25 + t210) * s25 + t99) * s25 + 59.0_dp 
     #* t13) * s25 + t238)
      t260 = t8 * t257
      t261 = 55.0_dp * t10
      t262 = t259 * t63
      t263 = t1 * t224
      t264 = t4 ** 2.0_dp
      t265 = s16 * s34
      t266 = t246 * t264
      t267 = t152 * (t152 * (t8 * t225 * t65 + t185 + s26 * (s26 * t263 
     #+ t217 * t31)) + t259 * t189 * t63 + 66.0_dp * s16 * t21 + t261 * 
     #t17 + t260 * t65 - t74 * t6 + t99 * t62 + t258 + 27.0_dp * t61 + t
     #36 * (((((t241 + t239) * s25 - t188) * s25 + t254) * s25 + t255) *
     # s25 + t256))
      t268 = 9.0_dp * t12
      t269 = t264 * t62
      t270 = t8 * t264
      t96 = t31 * ((t152 * (t152 * (s26 * t265 - t1 * t231 - t233 - t82)
     # - t176 * s25 - t81 * t57 - 17.0_dp * t237 - t50 * t209 * t62 - t2
     #62 - t228 - s26 * (t211 * t31 * t62 + t208 * t31)) - t47 * (s25 * 
     #((((22.0_dp * t26 - t83) * s25 - t123) * s25 - t124 - t164 * t119)
     # * s25 + t170) + (((((-26.0_dp * s16 + t235) * s25 - t205) * s25 -
     # 56.0_dp * t13) * s25 + 40.0_dp * t10) * s26 + ((((62.0_dp * s16 +
     # t204) * s25 - t234) * s25 + t182) * s25 - t10 * t42) * s25 + 27.0
     #_dp * t12) * s26 + t11)) * s34 + t266 * t91) + t267 - t262 * t3 * 
     #t264 + t269 * t1 * t92 + t270 * t93 * t65 + t144 * (s25 * (s25 * (
     #s25 * (s25 * (t215 - t96) - t183) + t136) - t10 * t227) + t268)
      t123 = 19.0_dp * s25
      t124 = t10 * t139
      t136 = t26 * t125
      t170 = 20.0_dp * t10
      t57 = t57 * t136
      t71 = t17 * (t17 * t78 - t71) + t170 - t57
      t78 = 84.0_dp * s16
      t176 = t86 * t189
      t182 = t264 * s26
      t183 = 24.0_dp * t10
      t205 = t269 * t6
      t28 = t1 * (t182 * t92 + t19 * t28) + (s34 * (t6 * t257 * t62 + t2
     #55 * s25 + t99 * s26 + t152 * (s34 * (t265 - s16 * (t15 + t243) - 
     #s25 * (t105 + t2 + t51)) + s26 * (t105 * t225 + t263) + t201 + t94
     #) + t176 * t65 + t241 * t18 + t21 * t30 - 14.0_dp * t237 + t256 + 
     #58.0_dp * t81) + s25 * (s25 * (s25 * (s25 * (s25 * (-t78 - t204) +
     # 141.0_dp * t9) + t223) + t220 + 40.0_dp * t164) - t177) - t6 * s2
     #6 * (s26 * t71 + s25 * (s25 * (s25 * (s25 * (t51 + t123) - t214) -
     # t161) + t124) + t268) - t11) * s34 - t86 * t245 * t264 + t205 * t
     #93 - t151 * (t62 * t1 * t211 + t137 * t209 + t18 * t49 + t187 * t9
     #5 + t86 * t65 + 45.0_dp * t132 + t41 + 47.0_dp * t67) + t26 * (s25
     # * (t17 * (s25 * (26.0_dp * s25 - t98) + 62.0_dp * t9) - t183) + t
     #268)
      t86 = 24.0_dp * t67
      t57 = s25 * (s25 * (s25 * (s25 * (t54 - t120) - t60) + t13 * t48) 
     #- 22.0_dp * t10) + t7 * (-t266 * t3 + t12) + (s34 * (s34 * (s34 * 
     #t225 - s25 * (t178 + t186 + t202) - t114) + s25 * (s16 * (t70 + t2
     #53) + s25 * (s16 + t147 + t79)) + t135) + t17 * (s16 * (t252 + t19
     #3) + s25 * (-t218 - t248)) - t170 + t57) * s34
      t37 = (s34 * (t62 * ((s26 * t225 + t224) * s26 + t217) + t81 + t18
     #5 * s26) - t31 * s26 * ((s25 * t65 + t208) * s26 + t26 * (((t107 +
     # t115) * s25 + t108) * s25 + t85)) + t65 * (-t139 * t209 - t213) +
     # t67 * ((-t120 - t43) * s25 - t212)) * s34 + t47 * (t82 * (((-t128
     # + t37) * s25 + 13.0_dp * t9) * s25 + t25) + (s16 * ((((-t186 + t1
     #12) * s25 + 29.0_dp * t9) * s25 + t203) * s25 + t41) + ((t246 * t1
     #89 + ((t150 - t128) * s25 + t188) * s25 + t135) * s26 + (((39.0_dp
     # * s16 + t239) * s25 + t216) * s25 + 29.0_dp * t13) * s25 + t238) 
     #* s26) * s26)
      t37 = t37 * s34 - t199 * ((((-t135 + 18.0_dp * t82) * s25 - t164 *
     # t8) * s25 + t163) * s25 + t10 * t125 + t62 * (-t49 * t17 * t66 + 
     #s16 * ((-t108 - t117) * s25 + t102)) + t173 * ((((t115 - t250) * s
     #25 + t222) * s25 - t101) * s25 + t232))
      t60 = t144 * t264
      t70 = t264 * s25
      t79 = s16 * t156
      t95 = s25 * (t14 * t131 * s34 + t1 * t19 + t152 * t189 + t82 * t7 
     #- t13 - 7.0_dp * t144 - t151) * t166
      t98 = 19.0_dp * s16
      t99 = t237 * (s25 * (s25 * (s25 * (-t98 + t43) + t195 * t9) - 13.0
     #_dp * t13) + t106)
      t102 = s25 + t128
      t107 = -t8 * t17 * t4 + t3 * t9
      t114 = s25 * (t179 + t43) + t9
      t115 = s25 * (s16 + t43) + t212
      t147 = s25 * t114
      t150 = t38 + t14
      t151 = t1 * t107
      t161 = t1 * t115
      t163 = s34 * (s34 * (t62 * (s26 * t115 + t147) + t237 - t36 * t144
     # * (-t2 + t128)) - t47 * (-s26 * (-s26 * t107 + t26 * (-t197 + t9)
     #) + t67 * (-t179 + t110))) + t237 * t125
      t164 = t17 * t264
      t170 = s25 * (t164 * t102 + t153 * t47) * t62
      t177 = t246 * (-t154 * (t150 * t62 + t31 * s25 * (s26 * t111 + t26
     #)) + t60 * t111)
      t179 = s25 + s26 + s56
      t38 = -t259 + t38
      t16 = s25 * (s25 + t230) + t16
      t185 = t26 * t1 + t125
      t188 = s25 * (t230 + t43) + t9
      t189 = t1 * t9 * t4 + t19
      t193 = t189 * t65
      t197 = -48.0_dp
      t201 = 117.0_dp * s25
      t203 = 20.0_dp * s25
      t204 = 60.0_dp * t13
      t208 = s25 * (t215 + t196) + t83
      t209 = t208 * s26
      t211 = 177.0_dp * s25
      t212 = 52.0_dp * t13
      t214 = 40.0_dp * s25
      t217 = 60.0_dp * t9
      t220 = t259 + t59
      t223 = t31 * t220
      t224 = 108.0_dp * t9
      t225 = (s25 * (165.0_dp * s16 + t214) + t229) * s25
      t227 = t225 - t254
      t112 = (106.0_dp * s16 + t112) * s25
      t228 = t112 + t167
      t230 = t10 * t17
      t232 = 11.0_dp * t61
      t216 = s25 * (58.0_dp * s16 + t2) + t216
      t241 = (-t128 - t203) * s25
      t248 = t192 + t241
      t250 = s16 * t216
      t256 = t90 + t15
      t257 = -t55 + t14
      t262 = -362.0_dp
      t100 = s25 * (s25 * (s25 * (s16 * t262 - t211) - t100) + 106.0_dp 
     #* t13) + t41
      t262 = -393.0_dp
      t89 = (s25 * (293.0_dp * s16 + t89) + 226.0_dp * t9) * s25
      t263 = t244 + t89
      t207 = (54.0_dp * s25 + t207) * s25
      t266 = t217 + t207
      t267 = (s25 * (s25 * (s25 * (s16 * t262 - t214) - 378.0_dp * t9) -
     # 121.0_dp * t13) + 52.0_dp * t10) * t1
      t271 = t263 * t8
      t272 = t266 * t7
      t273 = (-t98 - t242) * s25
      t274 = t1 * t227
      t275 = t8 * t228
      t276 = t1 * t248
      t277 = s16 * s26
      t278 = t269 * t38
      t279 = t38 * t264
      t280 = t111 * t63
      t72 = s34 * (t50 * t220 * t64 + s34 * (s34 * (-s34 * (-t31 * s16 *
     # (s26 * t256 + t26) + t62 * (t1 * t257 + t137)) - t132 * t139 - t1
     #73 * t216 - t200 * t197 - t276 * t62 + t30 * t63 - t86) + t50 * t8
     #2 * t66 - s26 * (t65 * (150.0_dp * t5 + t243) - t128 * (s25 * (t27
     #3 + t224) + t94)) - t62 * (s26 * t275 + t274) + 65.0_dp * t81) - t
     #173 * t100 - t267 * t62 + t271 * t65 + t272 * t63 - 70.0_dp * t230
     # - t258 + 69.0_dp * t61 + 150.0_dp * t74) - t47 * (t9 * (s25 * (s2
     #5 * (s25 * (131.0_dp * s16 + t190) + 142.0_dp * t9) - 41.0_dp * t1
     #3) - t72) + (s26 * (s26 * (s25 * (s25 * (436.0_dp * s16 + 144.0_dp
     # * s25) + 282.0_dp * t9) + 300.0_dp * t13) + s25 * (s25 * (s25 * (
     #708.0_dp * s16 + t236) + 514.0_dp * t9) + 672.0_dp * t13) + t255) 
     #+ s16 * (s25 * (s25 * (s25 * (698.0_dp * s16 + 234.0_dp * s25) + 3
     #03.0_dp * t9) + 336.0_dp * t13) - 69.0_dp * t10)) * s26 + t6 * t20
     #8 * t63)
      t216 = 38.0_dp * s16
      t239 = t239 - t216
      t255 = t270 * t38
      t258 = t255 * t188
      t270 = 18.0_dp * t12
      t281 = t252 * t189
      t133 = (s25 * (s25 * (286.0_dp * s16 + t133) + 199.0_dp * t9) + 24
     #3.0_dp * t13) * s25
      t282 = t165 * t10
      t283 = t282 + t133
      t24 = s25 * (s25 * (t24 + t202) + t130) + t135
      t130 = t253 * t208
      t240 = t240 * t220
      t253 = 300.0_dp * t5
      t284 = 72.0_dp
      t90 = s34 * (s34 * (s34 * (s34 * (t265 - s16 * (t90 + t15) + s26 *
     # (-t221 + t159 + t243)) + s26 * (s26 * (-t111 * t284 - t252) + t27
     #6) + t250) + s16 * (s25 * (-t273 - t224) - t94) + s26 * (s26 * (s2
     #6 * (45.0_dp * s26 + t253) + t228 * t6) + t274)) + s16 * t100 + s2
     #6 * (s26 * (s26 * (-t119 * t266 - t240) - t263 * t6) + t267)) + t6
     #2 * (s26 * (20.0_dp * t24 + t130) + t283 * t6) + s26 * (s25 * (s25
     # * (s25 * (s25 * (942.0_dp * s16 + t236) + 1212.0_dp * t9) + 975.0
     #_dp * t13) + 420.0_dp * t10) - 69.0_dp * t12) - t11 + t26 * (s25 *
     # (s25 * (s25 * (433.0_dp * s16 + t201) + 217.0_dp * t9) + 239.0_dp
     # * t13) - t261)
      t100 = 57.0_dp * s25
      t186 = (-t186 + t100) * s25
      t221 = 109.0_dp * t9
      t12 = t50 * (t12 + t21)
      t236 = (s25 * (s25 * (t98 + t159) - 46.0_dp * t9) + t40) * s25
      t243 = t236 + t33
      t261 = t182 * t38
      t206 = s25 * (-t259 + t206) - t222
      t84 = t47 * (-t42 * t81 * t4 - (s26 * (s26 * (-s26 * (s25 * (77.0_
     #dp * s16 + t171) + t217) - s25 * (s25 * (216.0_dp * s16 + t214) + 
     #166.0_dp * t9) - t244) - s16 * (s25 * (s25 * (212.0_dp * s16 + t21
     #1) + 97.0_dp * t9) - t212)) + t9 * (s16 - t110) * (s25 * (64.0_dp 
     #* s16 + t84) + t46)) * s26 + t223 * t64 - t138 * t1) * t152
      t32 = t199 * (-s26 * (s26 * (-s16 * (s25 * (s25 * (97.0_dp * s16 +
     # t201) + t234) + t13 * t165) + s26 * (-t209 - s25 * (s25 * (80.0_d
     #p * s16 + t203) + t32) - t204)) - t9 * (s25 * (s25 * (t190 + t98) 
     #+ 98.0_dp * t9) - t198)) - t11 + t132 * (s25 * (56.0_dp * s25 + t1
     #74) - 9.0_dp * t9)) * s34
      t23 = s56 * (s56 * (s34 * (s34 * (-t223 * s34 + t208) - t189 * t8)
     # + t1 * t194 + t26 * t206) + t6 * (s16 * (t261 + t10) + t153 * (t5
     # * t7 - s34 + t142)) - (24.0_dp * t189 * s26 + s34 * (s34 * (t277 
     #* t284 + t207 + 60.0_dp * t247) - s25 * (s25 * (t235 + t69) + 78.0
     #_dp * t9) - t204 - t105 * t208) + t236 + t33) * s34 + t26 * (s25 *
     # (s25 * (s16 * t29 - t203) - t83) - t94)) + s16 * (s25 * (s25 * (s
     #25 * (-s25 * t239 + t23 * t9) - 57.0_dp * t13) + t162) + t279 * t1
     #78 * t111 + t278 * t49) + t1 * (t11 + t156) + (s34 * (s34 * (s34 *
     # (45.0_dp * t62 + 150.0_dp * t149 + t112 + t167) - s26 * (t240 + t
     #272) - t244 - t89) + s26 * (t119 * t24 + t130) + t133 + t282) + s2
     #6 * (-t243 * t7 - t281) - t12 + t26 * (s25 * (-t186 + t221) - t175
     #)) * s34 - t251 * (t14 + t128 + t141)
      t12 = s56 * t23 - t139 * t67 * t103 + s16 * (t258 * s26 + t82 * (t
     #17 * (-68.0_dp * s16 + 92.0_dp * s25) + t101)) - 20.0_dp * s16 * (
     #-t278 * (s26 + t111) + t20) + (s34 * (s34 * (s34 * (s34 * (s34 * (
     #t137 - t55 + t14) + s26 * (t111 * t197 - t252) + t192 + t241) + s2
     #6 * (s26 * (t252 + t253) + t275) + t225 - t254) + s16 * (s25 * (s2
     #5 * (s25 * t262 - 378.0_dp * s16) - 121.0_dp * t9) + t212) + s26 *
     # (-t178 * t266 - t271) - 40.0_dp * t220 * t65 - 40.0_dp * t18) + s
     #16 * (s25 * (s25 * (s25 * (404.0_dp * s16 + 314.0_dp * s25) + 325.
     #0_dp * t9) + 140.0_dp * t13) - t10 * t195) + 20.0_dp * t62 * (t209
     # + t24) + 20.0_dp * t21 + t8 * t283 * s26) + t62 * (-80.0_dp * t18
     #9 * s26 - t119 * t243) + t8 * (-s26 * (t12 + t26 * (s25 * (t186 - 
     #t221) + t175)) + t11 - t20) + t26 * (s25 * (s25 * (s25 * (-t78 - t
     #100) + 236.0_dp * t9) - t180) - 67.0_dp * t10)) * s34
      t12 = s56 * t12 + s16 * (t17 * (s25 * (s25 * (s25 * (s25 * (-t259 
     #- t122) + 63.0_dp * t9) - t45) - t238) + t270) + t205 * t38 * t188
     # + 20.0_dp * t279 * t200 + t70 * t137 * t38 * t185 + t279 * t49 * 
     #t63) + (s34 * t90 - t47 * (s16 * (s25 * (s25 * (s25 * (s25 * (t191
     # + t123) - 144.0_dp * t9) + t13 * t197) + 129.0_dp * t10) - t68) +
     # s26 * (s26 * (t119 * (s25 * (s25 * (s25 * (t98 + t43) - t108) - t
     #40) + t33) + t281) + t8 * (s25 * (s25 * (s25 * (s25 * (t2 + t216) 
     #+ t229) - 141.0_dp * t13) + 93.0_dp * t10) + t270)))) * s34
      t12 = s56 * t12 + s16 * (t129 * t239 + t258 * t65) - t119 * s16 * 
     #(s16 * t18 * t19 - t280 * t279) + t31 * t277 * (t164 * t38 * t16 +
     # t155) + t50 * t26 * (t278 * t185 + t184) + t6 * s16 * (t20 * t9 +
     # t279 * t64) - s34 * (s34 * t72 + t199 * (-t11 * t50 + t62 * (s25 
     #* (s25 * (s25 * (t75 + t159) - 86.0_dp * t9) - 99.0_dp * t13) + 15
     #0.0_dp * t10) + t232 + t173 * (t17 * (19.0_dp * t17 - 155.0_dp * t
     #9) + 74.0_dp * t136 + t183) + 24.0_dp * t193 - 122.0_dp * t230 + t
     #74 * t29 + 105.0_dp * t138))
      t20 = t26 * (t255 * t185 * t65 + t237 * (s25 * t206 + t25))
      t23 = 0.1e1_dp / s16
      t24 = 0.1e1_dp / t39 ** 2.0_dp
      t25 = 0.1e1_dp / t34 ** 2.0_dp
      t29 = 0.1e1_dp / t35 ** 2.0_dp
      t3 = t76 * (s26 * (-t157 * (t31 * (t231 + t82) + t233) + t60 * (((
     #t215 + t59) * s25 - t44) * s25 + t40)) + t37 * s34 + s56 * t96 + t
     #148 * t28 + t160 * t57 + t168 * (s25 * (s25 * (s25 * (s25 * (s25 *
     # (t54 + t51) - 81.0_dp * t9) + 98.0_dp * t13) - 43.0_dp * t10) - t
     #269 * t3 * t119) + t1 * t11 + t31 * (s34 * (-s26 * t31 * t71 + 20.
     #0_dp * t131 * t17 * t62 - t124 * s25 - t153 * t226 - t51 * t18 - 1
     #9.0_dp * t21 + 51.0_dp * t237 - t268 + 17.0_dp * t81) + t138) + t8
     # * (-t154 * (t103 * t50 + t259 * t62 + 14.0_dp * t144 + t213 + 20.
     #0_dp * t82) + t182 * t93) + t152 * (t260 * s26 + t152 * (s16 * (t1
     #22 + t249) + s25 * (s26 * t42 + t210 + t219)) + t176 * t62 + t187 
     #* t48 + 59.0_dp * t132 + 62.0_dp * t18 + t238 + t86)) + t62 * (t70
     # * t91 + t79) - t70 * t3 * t64 + t95 + t264 * t93 * t63 + t264 * t
     #92 * t65 + t99) * MBL1011(1)
      t3 = (t31 * t261 * s16 * (t187 * t102 + t280) + s34 * (-t199 * t47
     # * (s16 * (t62 * (s25 * (s25 * (-t120 + t123) - t140) + t135) + t1
     #73 * (s25 * (s25 * (-t113 + t202) + t222) + t85) + t237 * t58 + t8
     #1 * t197 + t218 * t10) - t8 * (t11 - t193)) + t154 * (s34 * (s26 *
     # (t65 * (-t111 * t50 - t105) + t82 * (t2 + s16) * t139) + t62 * (s
     #26 * t248 + t250) + t81) - t1 * (-t62 * t63 + t230) - (s26 * (s26 
     #* (-s26 * t228 - t227) - s16 * (s25 * (s25 * (t242 + t98) - t224) 
     #- t94)) + t82 * (t50 * t66 + 65.0_dp * t26)) * s26 - t232 + t30 * 
     #t5 * t64) + s26 * (t1 * t65 - t80 - s26 * (s16 * t256 - s26 * t257
     #)) * t157 + t79 * t62 - t84 + t32) + s56 * t12 + t278 * s16 * (t16
     # * t17 + t62 * (t188 + t62)) + t20) * MBL1110(1) + t27 + t3
      t3 = (t1 * s25 * (-t166 * (t1 * t17 - t8 * s25 * (s16 + s34) + t9 
     #+ t152 - s34 * t128) + t184 + t4 * t62 * (-t245 + t181)) + t22 - 1
     #3.0_dp * t74 * t125 + t198 * t21) * MBL1001(0) + t3 * t23
      t11 = t76 ** 2.0_dp
      t6 = (-s25 * t163 + s56 * (t1 * t19 * (t67 - t269) - t31 * t19 * (
     #t182 * t102 + t187) - s34 * (s34 * (s25 * (t110 * t114 * s26 + t16
     #1 * t62 + t187 * t6 - t67 * t8) + s34 * (-s25 * (s26 * (s26 * t1 *
     # t150 + t14 * t111) + t144 * t31) + t246 * (t142 + t110) * s34)) -
     # s25 * (-t117 * t131 * t102 * s26 + t151 * t62 + t59 * t18 - t81 *
     # t7 - 7.0_dp * t237)) - t10 * t19) + t148 * (t1 * t19 * (-t182 + t
     #82) - s25 * (s34 * (s34 * (t161 * s26 + s34 * (s34 * (t142 + s25) 
     #- s25 * (t137 + t14 + t128) - t146 * s16) + t147) - t17 * (s25 * (
     #t14 + t128) - t172) - t151 * s26) + t21) - t61 * t31) - s25 * (-t1
     #87 * t31 + s34 * (s34 * (s34 * (-t150 + s34) + t115) - t107) + t12
     #5 * t17) * t168 - t264 * t19 * t65 - t170 - t177 + t52 * t21) * MB
     #L1101(0)
      result = t25 * t24 * ((t1 * ((-((-t126 - t152) * s34 + t31 * t104)
     # * s34 - t19 * t4 + t10 + t82 * t158 - t8 * t5 * t154) * t168 + t1
     #69 + t4 * t53 * t65) - ((((-s26 * ((t142 + t143 + t15) * s26 + t14
     #5) + s34 * t62) * s34 + t62 * (t134 * t31 + t137 * t5) + t67 + t12
     #1 * (t2 + t128) * s26) * s34 - t1 * (t126 * t65 + t81) - s26 * (s2
     #6 * t118 + t127) - t46 * t19) * s34 + t62 * (t104 * t105 + t97) + 
     #t36 * t26 * ((s25 * t87 - t9) * s25 + t88) + t106 * t17 + 7.0_dp *
     # t74) * s34 - t109 * s56 + t62 * (t77 + t73 + t56) + t184 + t116 *
     # (((-t196 - t43) * s25 - t46) * s25 + t52) * s26 + t129 - t41 * t1
     #9) * I300s16(1) + t29 * t3) / 12.0_dp + t23 * (t24 * (t6 * t25 - t
     #76 * t11 * t179 * MBL1111D2(0)) + t179 ** 2.0_dp * MBL1110(0) * t2
     #9) / 6.0_dp

           intHLs16s25s26s34s56x1411D6eps0 = result
       end function intHLs16s25s26s34s56x1411D6eps0

       function intHLs16s25s26s34s56x1411D6eps1()
           implicit none
           complex(dp) :: intHLs16s25s26s34s56x1411D6eps1
           complex(dp) ::  t1,t2,t3,t4

           complex(dp) :: result

      t1 = s25 + s26 + s56
      t2 = 3.0_dp
      t3 = s25 + s26 - s34 + s56
      t4 = 0.1e1_dp / 0.12e2_dp
      t3 = 0.1e1_dp / t3 ** 2.0_dp
      result = 0.1e1_dp / s16 * t3 * (t4 * (t1 * t2 - s34) + t1 ** 2.0_d
     #p * MBL1110(1) / 6.0_dp)

           intHLs16s25s26s34s56x1411D6eps1 = result
       end function intHLs16s25s26s34s56x1411D6eps1

       function intHLs250000x0111D2eps0()
           implicit none
           complex(dp) :: intHLs250000x0111D2eps0

           complex(dp) :: result

      result = 0.3e1_dp / 0.2e1_dp + s25 * I300s25(1) / 2.0_dp

           intHLs250000x0111D2eps0 = result
       end function intHLs250000x0111D2eps0

       function intHLs250000x0111D2eps1()
           implicit none
           complex(dp) :: intHLs250000x0111D2eps1

           complex(dp) :: result

      result = 0.1e1_dp / 0.2e1_dp

           intHLs250000x0111D2eps1 = result
       end function intHLs250000x0111D2eps1

       function intHLs250000x0112D2eps0()
           implicit none
           complex(dp) :: intHLs250000x0112D2eps0

           complex(dp) :: result

      result = -2.0_dp / s25 - I300s25(1)

           intHLs250000x0112D2eps0 = result
       end function intHLs250000x0112D2eps0

       function intHLs250000x0112D2eps1()
           implicit none
           complex(dp) :: intHLs250000x0112D2eps1

           complex(dp) :: result

      result = -0.1e1_dp / s25

           intHLs250000x0112D2eps1 = result
       end function intHLs250000x0112D2eps1

       function intHLs250000x0113D4eps0()
           implicit none
           complex(dp) :: intHLs250000x0113D4eps0

           complex(dp) :: result

      result = -0.1e1_dp / s25 / 2.0_dp - I300s25(1) / 4.0_dp

           intHLs250000x0113D4eps0 = result
       end function intHLs250000x0113D4eps0

       function intHLs250000x0113D4eps1()
           implicit none
           complex(dp) :: intHLs250000x0113D4eps1

           complex(dp) :: result

      result = -0.1e1_dp / s25 / 4.0_dp

           intHLs250000x0113D4eps1 = result
       end function intHLs250000x0113D4eps1

       function intHLs250000x0121D2eps0()
           implicit none
           complex(dp) :: intHLs250000x0121D2eps0

           complex(dp) :: result

      result = 4.0_dp / s25 + I300s25(0) + 2.0_dp * I300s25(1)

           intHLs250000x0121D2eps0 = result
       end function intHLs250000x0121D2eps0

       function intHLs250000x0121D2eps1()
           implicit none
           complex(dp) :: intHLs250000x0121D2eps1

           complex(dp) :: result

      result = 2.0_dp / s25 + I300s25(1)

           intHLs250000x0121D2eps1 = result
       end function intHLs250000x0121D2eps1

       function intHLs250000x0122D4eps0()
           implicit none
           complex(dp) :: intHLs250000x0122D4eps0

           complex(dp) :: result

      result = -0.3e1_dp / 0.2e1_dp / s25 - I300s25(1) / 2.0_dp

           intHLs250000x0122D4eps0 = result
       end function intHLs250000x0122D4eps0

       function intHLs250000x0122D4eps1()
           implicit none
           complex(dp) :: intHLs250000x0122D4eps1

           complex(dp) :: result

      result = -0.1e1_dp / s25 / 2.0_dp

           intHLs250000x0122D4eps1 = result
       end function intHLs250000x0122D4eps1

       function intHs160000x0111D0eps0()
           implicit none
           complex(dp) :: intHs160000x0111D0eps0

           complex(dp) :: result

      result = I300s16(0)

           intHs160000x0111D0eps0 = result
       end function intHs160000x0111D0eps0

       function intHs160000x0111D0eps1()
           implicit none
           complex(dp) :: intHs160000x0111D0eps1

           complex(dp) :: result

      result = I300s16(1)

           intHs160000x0111D0eps1 = result
       end function intHs160000x0111D0eps1

       function intHs160000x0112D2eps0()
           implicit none
           complex(dp) :: intHs160000x0112D2eps0

           complex(dp) :: result

      result = -2.0_dp / s16 - I300s16(1)

           intHs160000x0112D2eps0 = result
       end function intHs160000x0112D2eps0

       function intHs160000x0112D2eps1()
           implicit none
           complex(dp) :: intHs160000x0112D2eps1

           complex(dp) :: result

      result = -0.1e1_dp / s16

           intHs160000x0112D2eps1 = result
       end function intHs160000x0112D2eps1

       function intHs160000x0121D2eps0()
           implicit none
           complex(dp) :: intHs160000x0121D2eps0

           complex(dp) :: result

      result = 4.0_dp / s16 + I300s16(0) + 2.0_dp * I300s16(1)

           intHs160000x0121D2eps0 = result
       end function intHs160000x0121D2eps0

       function intHs160000x0121D2eps1()
           implicit none
           complex(dp) :: intHs160000x0121D2eps1

           complex(dp) :: result

      result = 2.0_dp / s16 + I300s16(1)

           intHs160000x0121D2eps1 = result
       end function intHs160000x0121D2eps1

       function intHs160000x0211D2eps0()
           implicit none
           complex(dp) :: intHs160000x0211D2eps0

           complex(dp) :: result

      result = -2.0_dp / s16 - I300s16(1)

           intHs160000x0211D2eps0 = result
       end function intHs160000x0211D2eps0

       function intHs160000x0211D2eps1()
           implicit none
           complex(dp) :: intHs160000x0211D2eps1

           complex(dp) :: result

      result = -0.1e1_dp / s16

           intHs160000x0211D2eps1 = result
       end function intHs160000x0211D2eps1

       function intHs160s26s34s56x1011D2eps0()
           implicit none
           complex(dp) :: intHs160s26s34s56x1011D2eps0
           complex(dp) ::  t1

           complex(dp) :: result

      t1 = 0.1e1_dp / 0.2e1_dp
      result = t1 * (-(s16 + s26 - s34 + s56) * MB1011(1) + MB1001(0) + 
     #1.0_dp)

           intHs160s26s34s56x1011D2eps0 = result
       end function intHs160s26s34s56x1011D2eps0

       function intHs160s26s34s56x1011D2eps1()
           implicit none
           complex(dp) :: intHs160s26s34s56x1011D2eps1

           complex(dp) :: result

      result = 0.1e1_dp / 0.2e1_dp

           intHs160s26s34s56x1011D2eps1 = result
       end function intHs160s26s34s56x1011D2eps1

       function intHs160s26s34s56x1012D2eps0()
           implicit none
           complex(dp) :: intHs160s26s34s56x1012D2eps0
           complex(dp) ::  t1,t2

           complex(dp) :: result

      t1 = s16 + s26 + s56
      t2 = s16 + s26 - s34 + s56
      t1 = 0.1e1_dp / t1
      result = 2.0_dp * t2 * MB1011(1) * t1 + t1 * (t2 * MB1011(0) - MB1
     #001(0))

           intHs160s26s34s56x1012D2eps0 = result
       end function intHs160s26s34s56x1012D2eps0

       function intHs160s26s34s56x1012D2eps1()
           implicit none
           complex(dp) :: intHs160s26s34s56x1012D2eps1
           complex(dp) ::  t1

           complex(dp) :: result

      t1 = s16 + s26 + s56
      t1 = 0.1e1_dp / t1
      result = t1 * ((s16 + s26 - s34 + s56) * MB1011(1) - 1.0_dp)

           intHs160s26s34s56x1012D2eps1 = result
       end function intHs160s26s34s56x1012D2eps1

       function intHs160s26s34s56x1013D4eps0()
           implicit none
           complex(dp) :: intHs160s26s34s56x1013D4eps0
           complex(dp) ::  t1,t2,t3

           complex(dp) :: result

      t1 = s16 + s26 - s34 + s56
      t2 = s16 + s26 + s56
      t3 = 0.1e1_dp / t2 ** 2.0_dp
      result = t3 * (-(-2.0_dp * s34 + 3.0_dp * t2) * MB1001(0) / 4.0_dp
     # + 0.3e1_dp / 0.2e1_dp * t1 ** 2.0_dp * MB1011(1)) - t1 * t3 * (-t
     #1 * MB1011(0) + 1.0_dp) / 2.0_dp

           intHs160s26s34s56x1013D4eps0 = result
       end function intHs160s26s34s56x1013D4eps0

       function intHs160s26s34s56x1013D4eps1()
           implicit none
           complex(dp) :: intHs160s26s34s56x1013D4eps1
           complex(dp) ::  t1,t2,t3

           complex(dp) :: result

      t1 = s16 + s26 + s56
      t2 = s16 + s26 - s34 + s56
      t3 = 0.1e1_dp / t1 ** 2.0_dp
      result = t3 * (s34 / 2.0_dp - 0.3e1_dp / 0.4e1_dp * t1 + t2 ** 2.0
     #_dp * MB1011(1) / 2.0_dp)

           intHs160s26s34s56x1013D4eps1 = result
       end function intHs160s26s34s56x1013D4eps1

       function intHs160s26s34s56x1020D2eps0()
           implicit none
           complex(dp) :: intHs160s26s34s56x1020D2eps0
           complex(dp) ::  t1

           complex(dp) :: result

      t1 = 0.1e1_dp / 0.2e1_dp
      result = t1 * (-(s16 + s26 + s56) * MB1011(1) + MB1001(0))

           intHs160s26s34s56x1020D2eps0 = result
       end function intHs160s26s34s56x1020D2eps0

       function intHs160s26s34s56x1020D2eps1()
           implicit none
           complex(dp) :: intHs160s26s34s56x1020D2eps1

           complex(dp) :: result

      result = 0.1e1_dp / 0.2e1_dp

           intHs160s26s34s56x1020D2eps1 = result
       end function intHs160s26s34s56x1020D2eps1

       function intHs160s26s34s56x1021D2eps0()
           implicit none
           complex(dp) :: intHs160s26s34s56x1021D2eps0
           complex(dp) ::  t1,t2

           complex(dp) :: result

      t1 = s16 + s26 + s56
      t2 = 0.1e1_dp / t1
      result = (s34 * MB1011(0) - (-2.0_dp * s34 + t1) * MB1011(1) + MB1
     #001(0)) * t2

           intHs160s26s34s56x1021D2eps0 = result
       end function intHs160s26s34s56x1021D2eps0

       function intHs160s26s34s56x1021D2eps1()
           implicit none
           complex(dp) :: intHs160s26s34s56x1021D2eps1
           complex(dp) ::  t1

           complex(dp) :: result

      t1 = s16 + s26 + s56
      t1 = 0.1e1_dp / t1
      result = t1 * (s34 * MB1011(1) + 1.0_dp)

           intHs160s26s34s56x1021D2eps1 = result
       end function intHs160s26s34s56x1021D2eps1

       function intHs160s26s34s56x1022D4eps0()
           implicit none
           complex(dp) :: intHs160s26s34s56x1022D4eps0
           complex(dp) ::  t1,t2,t3

           complex(dp) :: result

      t1 = s16 + s26 + s56
      t2 = s16 + s26 - s34 + s56
      t3 = 0.1e1_dp / t1 ** 2.0_dp
      result = -t3 * (-(MB1001(0) + 1.0_dp) * (-2.0_dp * s34 + t1) + (-6
     #.0_dp * s34 + t1) * t2 * MB1011(1)) / 2.0_dp + s34 * t2 * MB1011(0
     #) * t3

           intHs160s26s34s56x1022D4eps0 = result
       end function intHs160s26s34s56x1022D4eps0

       function intHs160s26s34s56x1022D4eps1()
           implicit none
           complex(dp) :: intHs160s26s34s56x1022D4eps1
           complex(dp) ::  t1,t2

           complex(dp) :: result

      t1 = s16 + s26 + s56
      t2 = 0.1e1_dp / t1 ** 2.0_dp
      result = t2 * (-s34 + t1 / 2.0_dp + s34 * (s16 + s26 - s34 + s56) 
     #* MB1011(1))

           intHs160s26s34s56x1022D4eps1 = result
       end function intHs160s26s34s56x1022D4eps1

       function intHs160s26s34s56x1031D4eps0()
           implicit none
           complex(dp) :: intHs160s26s34s56x1031D4eps0
           complex(dp) ::  t1,t2

           complex(dp) :: result

      t1 = s16 + s26 + s56
      t2 = 0.1e1_dp / t1 ** 2.0_dp
      result = s34 * t2 * (s34 * MB1011(0) + 1.0_dp) / 2.0_dp + t2 * ((2
     #.0_dp * s34 + t1) * MB1001(0) - (s16 ** 2.0_dp + s26 ** 2.0_dp + s
     #56 ** 2.0_dp - 6.0_dp * s34 ** 2.0_dp + 2.0_dp * (s16 + s26 + s34)
     # * s56 + 2.0_dp * (s16 + s26) * s34 + 2.0_dp * s16 * s26) * MB1011
     #(1)) / 4.0_dp

           intHs160s26s34s56x1031D4eps0 = result
       end function intHs160s26s34s56x1031D4eps0

       function intHs160s26s34s56x1031D4eps1()
           implicit none
           complex(dp) :: intHs160s26s34s56x1031D4eps1
           complex(dp) ::  t1,t2

           complex(dp) :: result

      t1 = s16 + s26 + s56
      t2 = 0.1e1_dp / t1 ** 2.0_dp
      result = t2 * (s34 / 2.0_dp + t1 / 4.0_dp + s34 ** 2.0_dp * MB1011
     #(1) / 2.0_dp)

           intHs160s26s34s56x1031D4eps1 = result
       end function intHs160s26s34s56x1031D4eps1

       function intHs16s25s26s34s56x1110D2eps0()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1110D2eps0
           complex(dp) ::  t1

           complex(dp) :: result

      t1 = 0.1e1_dp / 0.2e1_dp
      result = t1 * (s25 * MB1110(1) - (s16 + s26 + s56) * MB1011(1) + M
     #B1001(0) + 1.0_dp)

           intHs16s25s26s34s56x1110D2eps0 = result
       end function intHs16s25s26s34s56x1110D2eps0

       function intHs16s25s26s34s56x1110D2eps1()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1110D2eps1

           complex(dp) :: result

      result = 0.1e1_dp / 0.2e1_dp

           intHs16s25s26s34s56x1110D2eps1 = result
       end function intHs16s25s26s34s56x1110D2eps1

       function intHs16s25s26s34s56x1111D2eps0()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1111D2eps0
           complex(dp) ::  t1,t2,t3,t4,t5,t6,t7,t8

           complex(dp) :: result

      t1 = s16 + s26 - s34 + s56
      t2 = s12 + s15 + s25
      t3 = 2.0_dp
      t4 = -s16 * s25 + (s34 - s16 - s25 - s26) * s26 + (-s26 * t3 - s16
     # - s25 + s34 - s56) * s56
      t5 = s16 + s26 + s56
      t6 = s16 - s25
      t7 = s16 + s25 + s26 - s34 + s56
      t4 = 0.1e1_dp / t4
      t8 = t4 * t1
      t2 = 0.1e1_dp / t2
      result = t4 * (t1 * ((t1 * MB1111(0) + I300s16(0)) * s16 - t5 * MB
     #1011(0) - t7 * MB1110(0)) - (-s16 * t6 - s26 * t6 + (s16 * t3 + s2
     #5 + s26 - s34) * s34 - (s16 - s34 - s25) * s56) * MB1101(0)) / 2.0
     #_dp + t8 * ((t1 * MB1111(1) + I300s16(1)) * s16 - t5 * MB1011(1) -
     # t7 * MB1110(1)) + t8 * t3 * (s12 + s15 + s16 + s25 + s26 - s34 + 
     #s56) * t2

           intHs16s25s26s34s56x1111D2eps0 = result
       end function intHs16s25s26s34s56x1111D2eps0

       function intHs16s25s26s34s56x1111D2eps1()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1111D2eps1
           complex(dp) ::  t1,t2,t3,t4

           complex(dp) :: result

      t1 = s16 + s26 - s34 + s56
      t2 = s12 + s15 + s25
      t3 = 2.0_dp
      t3 = s16 * s25 + (-s34 + s16 + s25 + s26) * s26 + (s26 * t3 + s16 
     #+ s25 - s34 + s56) * s56
      t4 = 0.1e1_dp / t3
      t3 = 0.1e1_dp / t3
      t2 = 0.1e1_dp / t2
      result = -t4 * t1 * ((t1 * MB1111(1) + I300s16(1)) * s16 - (s16 + 
     #s26 + s56) * MB1011(1) - (s16 + s25 + s26 - s34 + s56) * MB1110(1)
     #) / 2.0_dp - t1 * (s12 + s15 + s16 + s25 + s26 - s34 + s56) * t2 *
     # t3

           intHs16s25s26s34s56x1111D2eps1 = result
       end function intHs16s25s26s34s56x1111D2eps1

       function intHs16s25s26s34s56x1112D2eps0()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1112D2eps0
           complex(dp) ::  t1

           complex(dp) :: result

      t1 = s16 + s26 - s34 + s56
      t1 = 0.1e1_dp / t1
      result = 0.1e1_dp / s16 * ((s25 * MB1101(0) - (s16 + s25 + s26 - s
     #34 + s56) * MB1111D2(0)) * t1 + MB1011(0))

           intHs16s25s26s34s56x1112D2eps0 = result
       end function intHs16s25s26s34s56x1112D2eps0

       function intHs16s25s26s34s56x1112D2eps1()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1112D2eps1

           complex(dp) :: result

      result = MB1011(1) / s16

           intHs16s25s26s34s56x1112D2eps1 = result
       end function intHs16s25s26s34s56x1112D2eps1

       function intHs16s25s26s34s56x1112D4eps0()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1112D4eps0
           complex(dp) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           complex(dp) ::  t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t4
           complex(dp) ::  t5,t6,t7,t8,t9

           complex(dp) :: result

      t1 = s16 - s25
      t2 = 2.0_dp
      t3 = t2 * s16
      t4 = s25 + t3
      t5 = s16 + s25
      t6 = s25 * t5
      t7 = s16 ** 2.0_dp
      t8 = 3.0_dp * s16
      t9 = s34 ** 2.0_dp
      t10 = s34 * t9
      t11 = s16 * s25
      t12 = t1 * s16
      t13 = t1 * s26
      t14 = (-t13 - t12) * s16
      t15 = (s26 * t4 + (-s25 * t2 - s26 + s34 - t8) * s34 + t6 + 3.0_dp
     # * t7) * s34
      t16 = -(-s25 * s34 - t3 * s34 - t11 + t7 + t9) * s56 + t14 + t15
      t17 = s25 ** 2.0_dp
      t18 = -t2 * (s34 * t5 + t11) + t17 + t7 + t9
      t19 = s16 + s25 - s34
      t20 = t2 * s26
      t11 = (-s34 + t5 + s26) * s26 + (t20 + t19 + s56) * s56 + t11
      t21 = s16 + s26
      t22 = s25 + s26
      t23 = 4.0_dp * s16
      t24 = s26 ** 2.0_dp
      t25 = t3 + t22
      t26 = t1 * t7
      t27 = s16 + s25 + s26 - s34 + s56
      t18 = 0.1e1_dp / t18
      t28 = 0.1e1_dp / t11
      t11 = 0.1e1_dp / t11
      t29 = 0.1e1_dp / 0.2e1_dp
      result = t29 * ((s16 + s26 - s34 + s56) * t27 * MB1111D2(0) * t28 
     #+ ((s34 * (s26 * t5 + (-t2 * t5 - s26 + s34) * s34 + s56 * t19 + t
     #17 + t7) * MB1101(0) + t27 * ((t25 - s34) * s34 - (s16 - s25 - s34
     #) * s56 - t12 - t13) * MB1110(1)) * s25 + (((t4 - s34) * s34 - t12
     #) * s56 + t14 + t15) * MB1001(0) + (-t2 * (-t26 * s26 + t10 * t25)
     # - (-t1 * t24 - t26) * s16 - ((-t9 - 6.0_dp * s16 * t21 - (5.0_dp 
     #* s26 + t8) * s25 - t17 - t24) * s34 + t21 * ((t23 + t20) * s16 + 
     #3.0_dp * s25 * t22)) * s34 - (-(-(3.0_dp * s25 + t3) * s34 + t12 +
     # t9) * s56 + t2 * (t14 + t10) + s34 * (6.0_dp * s25 * s26 + t23 * 
     #s26 - (6.0_dp * s16 + 5.0_dp * s25 + t20) * s34 + 3.0_dp * t6 + 6.
     #0_dp * t7)) * s56) * MB1011(1) - s16 * t16 * I300s16(1)) * t11 * t
     #18) - t16 * t18 * t11

           intHs16s25s26s34s56x1112D4eps0 = result
       end function intHs16s25s26s34s56x1112D4eps0

       function intHs16s25s26s34s56x1112D4eps1()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1112D4eps1

           complex(dp) :: result

      result = 0.0_dp

           intHs16s25s26s34s56x1112D4eps1 = result
       end function intHs16s25s26s34s56x1112D4eps1

       function intHs16s25s26s34s56x1113D4eps0()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1113D4eps0
           complex(dp) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           complex(dp) ::  t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30
           complex(dp) ::  t31,t32,t33,t34,t35,t36,t37,t38,t39,t4,t40,t41
           complex(dp) ::  t42,t43,t5,t6,t7,t8,t9

           complex(dp) :: result

      t1 = s25 ** 2.0_dp
      t2 = s25 * t1
      t3 = s16 ** 2.0_dp
      t4 = s16 * t3
      t5 = t3 + t1
      t6 = s16 + s25
      t7 = s25 * t6
      t8 = 3.0_dp
      t9 = 2.0_dp
      t10 = s34 * t6
      t11 = s34 ** 2.0_dp
      t12 = t3 + t1 + t11
      t13 = s26 * t5 + ((t6 * t8 + s26 - s34) * s34 - t8 * (t7 + t3) - t
     #9 * s26 * t6) * s34 + (-t10 * t9 + t12) * s56 + t2 + t4
      t14 = s16 * s25
      t10 = -t9 * (t10 + t14) + t12
      t12 = s16 + s25 - s34
      t15 = (-s34 + t6 + s26) * s26 + (s26 * t9 + s56 + t12) * s56 + t14
      t16 = s16 * t6
      t17 = t9 * t1
      t18 = t9 * s16
      t19 = t9 * s25
      t20 = 4.0_dp
      t21 = t20 * s26
      t22 = s26 ** 2.0_dp
      t23 = -s16 + s25
      t24 = s16 + s26 + s56
      t25 = s16 * t23
      t26 = 6.0_dp
      t27 = t8 * s25
      t28 = (-s16 * t20 - t27) * s25 + t3
      t29 = s16 + s26
      t30 = (s25 + t18) * s26
      t31 = 5.0_dp * s16
      t32 = t20 * s25
      t33 = t8 * t1
      t34 = t8 * s26
      t35 = t20 * t3
      t36 = (s16 - t27) * s26
      t37 = s16 * t26
      t37 = (-(-(s34 * t9 + t27) * s16 - (s25 - s34) * s34 + t3) * s56 +
     # t9 * s34 * t11 + t8 * s16 * (-t36 + t1) - ((t32 + t31 + t34) * s3
     #4 - (t37 + t27) * s26 + t20 * t25 - t17) * s34 - t4 + t35 * s25) *
     # s56 - s16 * ((t28 * t9 + t36 * t8) * s26 - t37 * t1) - (((-t21 - 
     #t18 - t27 + s34) * s34 + (10.0_dp * s16 + 8.0_dp * s25) * s26 + t8
     # * (t22 + t1)) * s34 + (-t20 * (-t25 * t9 + t1) - t30 * t8) * s26 
     #+ t9 * s16 * ((s16 + s25) * s25 + t3) - t2) * s34 - t23 * t4
      t38 = s16 + t19
      t39 = t23 * s25
      t40 = s16 + s25 + s26 - s34 + s56
      t15 = 0.1e1_dp / t15
      t24 = 0.1e1_dp / t24
      t10 = 0.1e1_dp / t10
      t41 = 0.1e1_dp / s16
      t42 = t15 * t10
      t10 = t13 * I300s16(1) * t10 + t40 ** 2.0_dp * MB1111D2(0) * t41
      t43 = t42 * t41
      t6 = t43 * (-((-s26 * t23 - t39) * s25 + (s26 * t38 + (-s26 - t18 
     #- t27 + s34) * s34 + t16 + t33) * s34 + ((t38 - s34) * s34 - t39) 
     #* s56) * MB1101(0) + t40 * (s26 * t6 + (-t6 * t9 - s26 + s34) * s3
     #4 + s56 * t12 + t1 + t3) * MB1110(1)) * s25
      t12 = 0.1e1_dp / 0.2e1_dp
      result = t12 * (-t15 * t10 + t24 * (t42 * ((t9 * (s56 ** 2.0_dp + 
     #t14 + t22) * s25 + (t17 + t16) * s26 - ((-s16 * t8 - s26 + s34 - t
     #19) * s34 + t8 * (s25 * s26 + t3) + t7 + t18 * s26) * s34 + ((-s34
     # * t8 + t21) * s25 - t9 * (s16 * s34 - t1) + t11 + t16) * s56 - t3
     # * t23) * MB1001(0) + (t4 * s25 * (t18 - t27) - s16 * ((s16 * (t1 
     #* t26 - t25) + (-t36 - t28) * s26) * s26 + t3 ** 2.0_dp) - ((((s16
     # - s26) * s34 + t9 * s26 * t29 - (s16 - t34) * s25 - t35) * s34 + 
     #t26 * t4 - s16 * t1 - t22 * (t32 + t31 + s26) - t33 * s26) * s34 -
     # t29 * (t20 * s16 * t5 - ((s16 - t19) * (-s25 + t18) + t30) * s26 
     #- t2 - t27 * t3)) * s34 - t37 * s56) * MB1011(1) * t41) + (s16 + s
     #26 - s34 + s56) * MB1011(0) * t41) + t6) - t43 * t13

           intHs16s25s26s34s56x1113D4eps0 = result
       end function intHs16s25s26s34s56x1113D4eps0

       function intHs16s25s26s34s56x1113D4eps1()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1113D4eps1
           complex(dp) ::  t1

           complex(dp) :: result

      t1 = s16 + s26 + s56
      t1 = 0.1e1_dp / t1
      result = 0.1e1_dp / s16 * t1 * ((s16 + s26 - s34 + s56) * MB1011(1
     #) - 1.0_dp) / 2.0_dp

           intHs16s25s26s34s56x1113D4eps1 = result
       end function intHs16s25s26s34s56x1113D4eps1

       function intHs16s25s26s34s56x1114D6eps0()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1114D6eps0
           complex(dp) ::  t1,t10,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109
           complex(dp) ::  t11,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t12
           complex(dp) ::  t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t13,t130
           complex(dp) ::  t131,t132,t133,t134,t135,t136,t137,t138,t139,t14,t140,t141
           complex(dp) ::  t142,t143,t144,t145,t146,t147,t148,t149,t15,t150,t151,t152
           complex(dp) ::  t153,t154,t155,t156,t157,t158,t159,t16,t160,t161,t162,t163
           complex(dp) ::  t164,t165,t166,t167,t168,t169,t17,t170,t171,t172,t173,t174
           complex(dp) ::  t175,t176,t177,t178,t179,t18,t180,t181,t182,t183,t184,t185
           complex(dp) ::  t186,t187,t188,t189,t19,t190,t191,t192,t193,t194,t195,t196
           complex(dp) ::  t197,t198,t199,t2,t20,t200,t201,t202,t203,t204,t205,t206
           complex(dp) ::  t207,t208,t209,t21,t210,t211,t212,t213,t214,t215,t216,t217
           complex(dp) ::  t218,t219,t22,t220,t221,t222,t223,t224,t225,t226,t227,t228
           complex(dp) ::  t229,t23,t230,t231,t232,t233,t234,t235,t236,t237,t238,t239
           complex(dp) ::  t24,t240,t241,t242,t243,t244,t245,t246,t247,t248,t249,t25
           complex(dp) ::  t250,t251,t252,t253,t254,t255,t256,t257,t258,t259,t26,t260
           complex(dp) ::  t261,t262,t263,t264,t265,t266,t267,t268,t269,t27,t270,t271
           complex(dp) ::  t272,t273,t274,t275,t276,t277,t278,t279,t28,t280,t281,t282
           complex(dp) ::  t283,t284,t285,t286,t287,t288,t289,t29,t290,t3,t30,t31
           complex(dp) ::  t32,t33,t34,t35,t36,t37,t38,t39,t4,t40,t41,t42
           complex(dp) ::  t43,t44,t45,t46,t47,t48,t49,t5,t50,t51,t52,t53
           complex(dp) ::  t54,t55,t56,t57,t58,t59,t6,t60,t61,t62,t63,t64
           complex(dp) ::  t65,t66,t67,t68,t69,t7,t70,t71,t72,t73,t74,t75
           complex(dp) ::  t76,t77,t78,t79,t8,t80,t81,t82,t83,t84,t85,t86
           complex(dp) ::  t87,t88,t89,t9,t90,t91,t92,t93,t94,t95,t96,t97
           complex(dp) ::  t98,t99

           complex(dp) :: result

      t1 = s16 - s25
      t2 = s16 + s25
      t3 = s16 ** 2.0_dp
      t4 = t3 ** 2.0_dp
      t5 = t4 ** 2.0_dp
      t6 = s16 * t5
      t7 = s16 * t3
      t8 = t7 * t4
      t9 = t3 * t4
      t10 = s16 * t4
      t11 = t3 * s25
      t12 = t11 * t2
      t13 = 8.0_dp
      t14 = 6.0_dp
      t15 = 5.0_dp
      t16 = s25 ** 2.0_dp
      t17 = t16 ** 2.0_dp
      t18 = s25 * t16
      t19 = t16 * t17
      t20 = s25 * t17
      t21 = s25 * t14
      t22 = t15 * s16
      t23 = 7.0_dp * t4
      t24 = 3.0_dp
      t25 = 18.0_dp * s25
      t26 = 27.0_dp * s16
      t27 = t24 * t3
      t28 = 46.0_dp * t10
      t29 = ((((-t26 - t25) * s25 - t27) * s25 - 87.0_dp * t7) * s25 - 3
     #1.0_dp * t4) * s25 + t28
      t30 = 52.0_dp
      t31 = 33.0_dp
      t32 = t3 * t15
      t33 = 65.0_dp * t4
      t34 = 2.0_dp
      t35 = 7.0_dp * t7
      t36 = ((t18 * t34 + t35) * s25 - 14.0_dp * t4) * s25 + t10 * t15
      t37 = t24 * s16
      t38 = t37 - t21
      t39 = 16.0_dp * t3
      t40 = (s25 * t38 - t39) * s25 + t35
      t41 = s16 + s26
      t42 = 47.0_dp * s25
      t43 = 24.0_dp * s16
      t44 = 20.0_dp * t3
      t45 = 59.0_dp * t7
      t46 = s16 * t16
      t47 = 50.0_dp
      t48 = 83.0_dp * t4
      t49 = 19.0_dp
      t50 = t14 * s16
      t51 = 25.0_dp * s25
      t52 = t3 * t49
      t53 = 28.0_dp * t7
      t54 = ((t51 - t50) * s25 - t52) * s25 + t53
      t55 = t54 * s26
      t56 = t15 * s25
      t57 = t14 * t3
      t58 = 16.0_dp * t7
      t59 = 11.0_dp * t4
      t60 = t13 * t10
      t61 = 30.0_dp * s25
      t62 = 145.0_dp * s25
      t63 = 148.0_dp * s25
      t64 = 31.0_dp * s16
      t65 = 42.0_dp * t3
      t66 = (39.0_dp * s25 + t64) * s25 + t65
      t67 = t3 + t16
      t68 = t66 * s26
      t69 = 80.0_dp * s16
      t70 = 27.0_dp * s25
      t71 = 85.0_dp * t4
      t72 = s16 * s25
      t73 = 135.0_dp * s25
      t74 = 75.0_dp * s25
      t75 = 835.0_dp * t4
      t76 = (((588.0_dp * s16 + t74) * s25 + 1052.0_dp * t3) * s25 + 115
     #0.0_dp * t7) * s25
      t77 = t75 + t76
      t78 = 654.0_dp * t7
      t79 = ((651.0_dp * s16 + 230.0_dp * s25) * s25 + 789.0_dp * t3) * 
     #s25
      t80 = t78 + t79
      t81 = (128.0_dp * s16 + 81.0_dp * s25) * s25 + 117.0_dp * t3
      t82 = 28.0_dp * s16
      t83 = -t82 - t70
      t84 = s26 ** 2.0_dp
      t85 = t84 ** 2.0_dp
      t86 = t84 * t85
      t87 = s26 * t85
      t88 = s26 * t84
      t89 = s16 * ((((534.0_dp * s16 + t73) * s25 + 769.0_dp * t3) * s25
     # + 798.0_dp * t7) * s25 + 525.0_dp * t4)
      t90 = 140.0_dp * t4
      t91 = 100.0_dp * s25
      t92 = 630.0_dp * t7
      t93 = ((562.0_dp * s16 + t91) * s25 + 881.0_dp * t3) * s25
      t94 = t92 + t93
      t95 = 406.0_dp * t3
      t96 = (470.0_dp * s16 + 200.0_dp * s25) * s25
      t97 = t95 + t96
      t98 = 106.0_dp
      t99 = 88.0_dp * s25
      t100 = t98 * s16
      t101 = t100 + t99
      t102 = s16 * (((560.0_dp * s16 + 170.0_dp * s25) * s25 + 722.0_dp 
     #* t3) * s25 + 455.0_dp * t7)
      t103 = 78.0_dp * s25
      t104 = 120.0_dp * s25
      t105 = 271.0_dp * t3
      t74 = (286.0_dp * s16 + t74) * s25
      t106 = t105 + t74
      t107 = 131.0_dp * s16
      t108 = 91.0_dp * s25
      t109 = t108 + t107
      t110 = s16 * ((299.0_dp * s16 + t104) * s25 + 231.0_dp * t3)
      t111 = 57.0_dp
      t112 = t111 * s25
      t113 = 7.0_dp * s16
      t114 = t34 * s16
      t115 = t114 + s25
      t116 = 9.0_dp
      t117 = 17.0_dp * s26
      t118 = s16 * t116
      t119 = 4.0_dp
      t120 = t15 * s26
      t121 = t119 * t3
      t122 = 48.0_dp * t3
      t123 = 79.0_dp * t4
      t124 = 12.0_dp
      t125 = 95.0_dp * s16
      t126 = 10.0_dp * s25
      t127 = 36.0_dp * t3
      t128 = t124 * t7
      t129 = t15 * t54
      t130 = 15.0_dp * s25
      t131 = 176.0_dp * s16
      t132 = (((443.0_dp * s16 + t62) * s25 + 444.0_dp * t3) * s25 + 415
     #.0_dp * t7) * s25
      t133 = 581.0_dp * t4
      t134 = 82.0_dp
      t135 = 128.0_dp * t7
      t136 = ((97.0_dp * s16 + 74.0_dp * s25) * s25 + t134 * t3) * s25 +
     # t135
      t137 = t120 * t66
      t138 = t136 * t13
      t139 = (t133 + t132) * t24
      t140 = t81 * t13
      t141 = t80 * t24
      t142 = 35.0_dp * s26
      t143 = t97 * t14 / 2.0_dp
      t144 = 20.0_dp * s25
      t145 = 40.0_dp * s16
      t146 = 21.0_dp * s16
      t147 = 10.0_dp * s26
      t148 = s34 ** 2.0_dp
      t149 = t148 ** 2.0_dp
      t150 = t149 ** 2.0_dp
      t151 = s34 * t148
      t152 = t151 * t149
      t153 = t148 * t149
      t154 = s34 * t149
      t155 = ((((-t147 - t113) * t148 - 76.0_dp * t88 - t110 - s26 * (t2
     #4 * t109 * s26 + t34 * t106)) * s34 + (((t13 * t101 / 2.0_dp + t14
     #2) * s26 + t143) * s26 + t34 * t94) * s26 + t102) * s34 - (((-t120
     # * t83 + t140) * s26 + t141) * s26 + t34 * t77) * s26 - t89) * s34
     # + s16 * (((((287.0_dp * s16 + t112) * s25 + 329.0_dp * t3) * s25 
     #+ 293.0_dp * t7) * s25 + 377.0_dp * t4) * s25 + 357.0_dp * t10) + 
     #(((t137 + t138) * s26 + t139) * s26 + t119 * (((((t131 + t130) * s
     #25 + 277.0_dp * t3) * s25 + 249.0_dp * t7) * s25 + 287.0_dp * t4) 
     #* s25 + 320.0_dp * t10)) * s26
      t156 = t1 * (((((-s16 * t31 - t21) * s25 - t32) * s25 - 63.0_dp * 
     #t7) * s25 - t33) * s25 + t10 * t30)
      t157 = t8 * s25
      t158 = 15.0_dp * s16
      t159 = t15 * t85
      t31 = t3 * t31
      t160 = (158.0_dp * s16 + t42) * s25
      t161 = 53.0_dp * t7
      t162 = 112.0_dp * t4
      t163 = 271.0_dp * t10
      t164 = 67.0_dp * s25
      t165 = 44.0_dp * s16
      t166 = (t165 + t164) * s25
      t167 = 88.0_dp * t7
      t168 = 139.0_dp * t4
      t169 = ((-t57 + t166) * s25 - t167) * s25 + t168
      t170 = 118.0_dp
      t171 = 123.0_dp
      t172 = 18.0_dp * t7
      t173 = 10.0_dp * t55
      t174 = 352.0_dp * s16
      t175 = 10.0_dp * t68
      t176 = -10.0_dp * t83 * s26
      t177 = 70.0_dp * s26
      t178 = 60.0_dp
      t179 = 84.0_dp * t4
      t180 = t3 * t116
      t181 = 18.0_dp * s16
      t182 = 10.0_dp * t1 * t40
      t183 = t24 * s26
      t184 = t9 * s25
      t185 = 58.0_dp
      t186 = t7 * t185
      t187 = 42.0_dp * t7
      t188 = t4 + t149
      t189 = 23.0_dp * t7
      t190 = t41 ** 2.0_dp
      t77 = (((-((t120 + t113) * s26 + t121) * s34 + ((30.0_dp * t115 + 
     #t117) * s26 + t118 * (t113 + t56)) * s26 + t3 * (23.0_dp * s25 + t
     #82)) * s34 - t3 * ((112.0_dp * s16 + t112) * s25 + 84.0_dp * t3) -
     # t49 * t85 - ((s26 * t109 + t106) * s26 + t110) * s26) * s34 + t3 
     #* (((190.0_dp * s16 + t103) * s25 + 215.0_dp * t3) * s25 + 140.0_d
     #p * t7) + 7.0_dp * t87 + (((s26 * t101 + t97) * s26 + t94) * s26 +
     # t102) * s26) * s34 + t3 * ((((-168.0_dp * s16 - 62.0_dp * s25) * 
     #s25 - 218.0_dp * t3) * s25 - 200.0_dp * t7) * s25 - t90) - (((-t83
     # * t84 + t80) * s26 + t77) * s26 + t89) * s26 - t34 * t81 * t85
      t52 = (s25 * (((-t118 + t21) * s25 + t52) * s25 - t189) + ((s34 * 
     #t83 + t66) * s34 - t54) * s34 + 7.0_dp * t188) * s56 + ((t124 * t1
     #8 + t187) * s25 - t179) * s25 - t154 * t49 + 30.0_dp * t10 + ((((t
     #142 + t100 + t99) * s34 - s16 * (234.0_dp * s16 + 140.0_dp * s26) 
     #- (256.0_dp * s16 + 162.0_dp * s25 + 135.0_dp * s26) * s25) * s34 
     #+ ((194.0_dp * s16 + t63) * s25 + 164.0_dp * t3) * s25 + t137 + 25
     #6.0_dp * t7) * s34 + ((t57 - t166) * s25 + t167) * s25 - t168 - t1
     #29 * s26) * s34 + t120 * t1 * t40
      t52 = t52 * s56 + (((((t26 + t21) * s25 - 28.0_dp * t3) * s25 + t1
     #86) * s25 + t34 * t4) * s25 - 117.0_dp * t10) * s25 + (t182 * s26 
     #+ 24.0_dp * t36) * s26 + t30 * t9 + 17.0_dp * t153 + ((((-(76.0_dp
     # * s26 + t107 + t108) * s34 + (424.0_dp * s16 + 352.0_dp * s25 + t
     #177) * s26 + t95 + t96) * s34 - (t140 + t176) * s26 - t78 - t79) *
     # s34 + (t175 + t138) * s26 + t132 + t133) * s34 + (((-t31 - t160) 
     #* s25 + t161) * s25 + t162) * s25 + (-t119 * t169 - t173) * s26 - 
     #t163) * s34
      t31 = t52 * s56 - t15 * t152 + 46.0_dp * t8 + ((((((s16 * t178 + 5
     #1.0_dp * s26 + t61) * s34 - (393.0_dp * s16 + 273.0_dp * s25 + 114
     #.0_dp * s26) * s26 - t105 - t74) * s34 + ((t14 * t101 + t177) * s2
     #6 + t143) * s26 + t92 + t93) * s34 - ((t124 * t81 + t176) * s26 + 
     #t141) * s26 - t75 - t76) * s34 + ((((t61 + t174) * s25 + 554.0_dp 
     #* t3) * s25 + 498.0_dp * t7) * s25 + 574.0_dp * t4) * s25 + ((t124
     # * t136 + t175) * s26 + t139) * s26 + 640.0_dp * t10) * s34 + ((((
     #(-s16 * t170 - t56) * s25 - 121.0_dp * t3) * s25 - t172) * s25 + t
     #171 * t4) * s25 + 28.0_dp * t10) * s25 + ((-t14 * t169 - t173) * s
     #26 - t24 * ((((t31 + t160) * s25 - t161) * s25 - t162) * s25 + t16
     #3)) * s26 - 265.0_dp * t9) * s34 - 56.0_dp * t10 * t16 + t179 * t1
     #8 - 24.0_dp * t7 * t17 + t180 * t20 + t181 * t19 + t182 * t88 + t1
     #83 * t156 + 36.0_dp * t36 * t84 - 77.0_dp * t184
      t31 = t31 * s56 + t24 * ((s16 * (t130 + t146) + (t145 + t117 + t14
     #4) * s26) * t153 + t156 * t84) + t3 * (t16 * ((((t25 - t158) * s25
     # + t57) * s25 + 39.0_dp * t7) * s25 - 45.0_dp * t4) + 21.0_dp * t9
     #) + (t155 * s34 - t41 * ((((((t126 + t125) * s25 + t127) * s25 - t
     #128) * s25 - 81.0_dp * t4) * s25 + 23.0_dp * t10) * s25 + 133.0_dp
     # * t9 + s26 * (((((51.0_dp * s16 + 268.0_dp * s25) * s25 + t57) * 
     #s25 - 257.0_dp * t7) * s25 + 416.0_dp * t4) * s26 + ((((206.0_dp *
     # s16 + 141.0_dp * s25) * s25 + t122) * s25 - 165.0_dp * t7) * s25 
     #- t123) * s25 + 397.0_dp * t10) + t129 * t88)) * s34 + 24.0_dp * t
     #36 * t88 - 24.0_dp * t157 + t114 * t1 * t29 * s26 + t159 * t1 * t4
     #0
      t52 = s16 * t1
      t54 = t2 * s26
      t66 = t3 * (t54 * t1 * (-t12 * t13 + t18 * (-t21 + t22) + t23) - t
     #184)
      t74 = t7 + t18
      t75 = t116 * s25
      t76 = s16 * t13
      t78 = 10.0_dp * t7
      t79 = t4 * t16
      t29 = t119 * t9 * t74 + t14 * (t19 * t7 + t36 * t85) + t24 * t66 +
     # s34 * ((t77 * s34 + t41 * (84.0_dp * t4 * t67 + (((t68 + ((155.0_
     #dp * s16 + t63) * s25 + 133.0_dp * t3) * s25 + 214.0_dp * t7) * s2
     #6 + (((295.0_dp * s16 + t62) * s25 + 289.0_dp * t3) * s25 + 282.0_
     #dp * t7) * s25 + 367.0_dp * t4) * s26 + ((((207.0_dp * s16 + t61) 
     #* s25 + 259.0_dp * t3) * s25 + 209.0_dp * t7) * s25 + 292.0_dp * t
     #4) * s25 + 273.0_dp * t10) * s26 + t72 * (t16 * ((t70 + t69) * s25
     # + 70.0_dp * t3) + t71))) * s34 - t190 * ((((((t43 + t56) * s25 + 
     #t57) * s25 - t58) * s25 - t59) * s25 + t60) * s25 + 28.0_dp * t9 +
     # ((t14 * t46 * t1 + (-t47 * t7 + 67.0_dp * t18) * s25 + t48 + t55)
     # * s26 + ((((t42 + t43) * s25 + t44) * s25 - t45) * s25 + t23) * s
     #25 + 77.0_dp * t10) * s26)) + t31 * s56 + t52 * t29 * t84 + t1 * t
     #40 * t87 + t156 * t88 + t79 * (t16 * (t76 - t75) - t78)
      t31 = -t34 * (s34 * t2 + t72) + t3 + t16 + t148
      t36 = s16 + s26 + s56
      t40 = t34 * s26
      t42 = (-s34 + t2 + s26) * s26 + (t40 + s16 + s25 - s34 + s56) * s5
     #6 + t72
      t55 = t24 * s25
      t62 = t119 * t4
      t63 = 7.0_dp * t10
      t66 = t15 * t7
      t68 = t4 * t49
      t77 = t24 * t52 * t18
      t80 = ((t14 * t18 + t66) * s25 - t68) * s25
      t81 = -t11 * t34 - t18 + t7
      t83 = 22.0_dp
      t89 = 42.0_dp * t10
      t92 = ((t16 * (-s25 * t83 - t43) + t128) * s25 + t14 * t4) * s25
      t93 = 32.0_dp * s25
      t94 = (-t118 - t93) * t16
      t95 = t7 * t49
      t96 = 40.0_dp * t4
      t97 = (t94 + t95) * s25 - t96
      t99 = t34 * t74 - t11
      t100 = t34 * s25
      t101 = t124 * t4
      t102 = 96.0_dp * s16
      t105 = 105.0_dp * t4
      t106 = (((65.0_dp * s25 + t102) * s25 + 85.0_dp * t3) * s25 + 86.0
     #_dp * t7) * s25
      t107 = t3 * t111
      t108 = 80.0_dp * t7
      t109 = ((67.0_dp * s16 + 68.0_dp * s25) * s25 + t107) * s25
      t110 = t108 + t109
      t117 = t14 * t67 + t72 * t15
      t120 = t124 * s25
      t129 = t4 * t47
      t132 = t24 * t88
      t133 = 70.0_dp * t7
      t136 = ((s25 * t47 + 87.0_dp * s16) * s25 + 92.0_dp * t3) * s25
      t137 = 72.0_dp * s25
      t138 = 80.0_dp * t3
      t139 = (t137 + t125) * s25
      t140 = t138 + t139
      t141 = 110.0_dp
      t142 = 86.0_dp * s16
      t143 = 70.0_dp * t4
      t155 = 105.0_dp * t3
      t156 = (141.0_dp * s16 + 85.0_dp * s25) * s25
      t160 = 38.0_dp * s25
      t161 = 115.0_dp
      t163 = 40.0_dp * s25
      t166 = t3 * t161
      t167 = s16 * t41
      t168 = 59.0_dp * s16
      t169 = 14.0_dp * s16
      t173 = t16 * t67
      t175 = t24 * t16
      t176 = 11.0_dp * t3
      t177 = t13 * t7
      t179 = 18.0_dp * s26
      t182 = t179 * t99
      t191 = t116 * s26
      t192 = t191 * t117
      t193 = t49 * s25
      t194 = t13 * s26
      t195 = t20 * t2
      t196 = t124 * s16
      t197 = t3 * t13
      t198 = -t10 + t154
      t199 = t1 * t18
      t200 = t119 * t2
      t94 = (-t13 * t198 + t77 + t80 - (((-(t145 + t160 + t191) * s34 + 
     #t138 + t139 + 36.0_dp * t54) * s34 - t108 - t109 - t192) * s34 + (
     #-t94 - t95) * s25 + t182 + t96) * s34 + t191 * t1 * t81) * s56 + t
     #24 * t195 - t34 * (t151 * (t148 * (t194 + t146 + t193) + (t140 + 1
     #8.0_dp * t54) * s26 + t133 + t136) - t1 * (s16 * ((-t176 - t175) *
     # s25 + t177) - t14 * t173) * s26) + ((t148 * ((76.0_dp * s25 + t69
     # + t191) * s26 + t155 + t156) + (t110 * t34 + t192) * s26 + t105 +
     # t106) * s34 + (t34 * t97 - t182) * s26 - t89 + t92) * s34 + 7.0_d
     #p * t9 + 7.0_dp * t153 + t11 * (((t196 - t21) * s25 - t197) * s25 
     #- 11.0_dp * t7) + t116 * t1 * t81 * t84
      t95 = s56 ** 2.0_dp
      t108 = s56 * t95
      t81 = (t11 * (t100 - t37) + ((t148 + t117) * s34 - t34 * t99) * s3
     #4 + t4 - t199 - t200 * t151) * t108 - t52 * t20 + t1 * t81 * t88
      t60 = t24 * t81 - t34 * (-t8 + t152) + ((t77 + t60 + t80) * s26 + 
     #t1 * (((t16 * (-t50 - t55) - t128) * s25 - t62) * s25 + t63)) * s2
     #6 + ((((((7.0_dp * s26 + t169 + t120) * s34 - (38.0_dp * s26 + t61
     # + t168) * s25 - t13 * t84 - 42.0_dp * t167) * s34 + ((99.0_dp * s
     #16 + t163) * s25 + t166) * s25 + t132 + t133 + s26 * ((t145 + t160
     #) * s26 + t155 + t156)) * s34 - ((t54 * t124 + t140) * s26 + t34 *
     # (t133 + t136)) * s26 - t12 * t141 - t18 * (t61 + t142) - t143) * 
     #s34 + ((((t165 + t120) * s25 + 40.0_dp * t3) * s25 + 36.0_dp * t7)
     # * s25 + t129) * s25 + t89 + s26 * (s26 * t110 + t105 + t106) + t1
     #32 * t117) * s34 + ((t18 * (-t100 - t158) + t101) * s25 - t63) * s
     #25 - 14.0_dp * t9 + s26 * (s26 * t97 - t89 + t92) - t14 * t99 * t8
     #8) * s34 + t94 * s56 - t184 + t79 * (-t113 + t21)
      t77 = -t116 * s25 * (t16 * (s16 + s25) + t7) + t23
      t80 = 45.0_dp * s16
      t81 = (t16 * (-t75 - t80) - t189) * s25 + 17.0_dp * t4
      t92 = t3 * t34
      t94 = (-t25 - t158) * s25
      t97 = -t55 + s16
      t99 = 53.0_dp * s25
      t105 = 20.0_dp * t7
      t100 = t100 + s16
      t106 = 11.0_dp * t16
      t109 = t88 * t100
      t110 = 147.0_dp
      t113 = 30.0_dp * t3
      t117 = t119 * s25
      t132 = 340.0_dp * t7
      t136 = ((s25 * t170 + 314.0_dp * s16) * s25 + 357.0_dp * t3) * s25
      t138 = t132 + t136
      t139 = 63.0_dp * s25
      t140 = 91.0_dp * s16
      t155 = (t140 + t139) * s25 + 76.0_dp * t3
      t156 = 20.0_dp * s16 + t70
      t160 = ((202.0_dp * s16 + t163) * s25 + 301.0_dp * t3) * s25
      t182 = 245.0_dp * t7
      t192 = (261.0_dp * s16 + 122.0_dp * s25) * s25
      t201 = 235.0_dp * t3
      t202 = t192 + t201
      t203 = 73.0_dp * s16
      t204 = t203 + t164
      t205 = 88.0_dp * s16
      t206 = 28.0_dp * s25
      t207 = 35.0_dp * s16
      t208 = 40.0_dp * t7
      t209 = 26.0_dp * s16
      t210 = 11.0_dp * s25
      t211 = t7 * t14
      t212 = 13.0_dp * t4
      t213 = 47.0_dp * s16
      t214 = (t213 + t112) * s25 + t113
      t215 = (((465.0_dp * s16 + 348.0_dp * s25) * s25 + 369.0_dp * t3) 
     #* s25 + 474.0_dp * t7) * s26
      t216 = 15.0_dp * t85
      t217 = 152.0_dp * t7
      t218 = t24 * t84
      t219 = 20.0_dp * s26
      t66 = t1 * ((-t92 + t94) * s25 + t66)
      t45 = (t148 * (t148 * (63.0_dp * s25 * s26 + t142 * s26 + 20.0_dp 
     #* t16 + 63.0_dp * t3 + 65.0_dp * t72 + 21.0_dp * t84) + (((t126 + 
     #t205) * s25 + 143.0_dp * t3) * s25 + t217) * s25 + t90 + s26 * (t8
     #4 * t34 * t156 + t138) + t218 * t155) + t41 * (((((t210 + t209) * 
     #s25 + t197) * s25 - t211) * s25 - t212 + 15.0_dp * t109) * s25 + 3
     #5.0_dp * t10 + s26 * ((((s25 * t98 - t50) * s25 + t176) * s25 + t2
     #08) * s26 + (((17.0_dp * s16 + t103) * s25 + t122) * s25 - t45) * 
     #s25 + t48))) * s34
      t48 = 36.0_dp * s16
      t197 = 24.0_dp * t3
      t220 = t14 * s26
      t221 = s16 * t74
      t222 = 26.0_dp * t16
      t223 = 14.0_dp * t3
      t224 = 41.0_dp * t221
      t225 = 13.0_dp * t3
      t226 = ((t196 + t99) * s25 + t225) * s25 + t105
      t227 = 25.0_dp * t3
      t228 = t170 * t10
      t229 = t3 * t16
      t230 = t61 * t88
      t231 = t14 * t84
      t232 = 201.0_dp * s25
      t233 = 30.0_dp * s26
      t234 = 48.0_dp * s25
      t235 = t183 * t1
      t236 = t58 * s25
      t237 = t61 * t84
      t238 = t119 * s26
      t239 = 40.0_dp * s26
      t240 = t178 * s26
      t241 = 13.0_dp * s16
      t242 = 15.0_dp * s26
      t65 = ((-(s16 * (t120 + t169) + (13.0_dp * s26 + t144 + t207) * s2
     #6) * s34 + s16 * ((t206 + t168) * s25 + t65) + ((14.0_dp * s26 + t
     #142 + t139) * s26 + t34 * ((65.0_dp * s16 + t144) * s25 + 63.0_dp 
     #* t3)) * s26) * s34 - s16 * (((t93 + t205) * s25 + t166) * s25 + t
     #133) - t159 - ((s26 * t204 + t202) * s26 + t160 + t182) * s26) * s
     #34 + s16 * ((((66.0_dp * s16 + t25) * s25 + 86.0_dp * t3) * s25 + 
     #t141 * t7) * s25 + t143) + ((t156 * t84 + t138) * s26 + (((t131 + 
     #t144) * s25 + 286.0_dp * t3) * s25 + 304.0_dp * t7) * s25 + 280.0_
     #dp * t4) * s26 + t34 * t155 * t88
      t61 = -(s25 * (((-t37 + t25) * s25 - t225) * s25 - t35 + t242 * t9
     #7 * t1) + t15 * t188 - ((t156 * s34 - (t213 + t112 + t242) * s25 -
     # t113) * s34 + (s16 * (t233 + t241) + (t196 + t240 + t99) * s25) *
     # s25 + t105) * s34) * s56 - ((((t75 + t48) * s25 - 45.0_dp * t3) *
     # s25 + t189) * s25 - t96) * s25 + t34 * s34 * (t34 * t226 * s26 + 
     #t61 * t100 * t84 + t148 * (s16 * (76.0_dp * s16 + t239) + (54.0_dp
     # * s26 + t139 + t140) * s25) + 26.0_dp * t17 + t224 + 14.0_dp * t2
     #29 - t236) - 17.0_dp * t10 + 14.0_dp * t154 - t148 * (t11 * t171 +
     # t148 * (t203 + t164 + t219) + t238 * t214 + 116.0_dp * t18 + t237
     # + 155.0_dp * t46 + 158.0_dp * t7) - t237 * t97 * t1 - t238 * t66
      t58 = t61 * s56 - t10 * (-t234 + t146) - 13.0_dp * t153 - 27.0_dp 
     #* t46 * ((t16 - t3) * s25 + t7) - ((216.0_dp * s16 * t18 + ((-(42.
     #0_dp * s26 + t139 + t142) * s34 + (219.0_dp * s16 + t232 + t233) *
     # s26 + t192 + t201) * s34 - t14 * s26 * (s26 * t156 + t155) - t132
     # - t136) * s34 + t111 * t17 + t231 * t214 + t73 * t7 + t215 + 207.
     #0_dp * t229 + t230 + 275.0_dp * t4) * s34 - ((((104.0_dp * s16 + t
     #210) * s25 + t227) * s25 + t187) * s25 - 72.0_dp * t4 + t109 * t17
     #8) * s25 - t14 * s26 * (t226 * s26 + s25 * ((t222 + t223) * s25 - 
     #t58) + t224) - t228) * s34 - t230 * t97 * t1 - t231 * t66 - t235 *
     # t81
      t45 = t58 * s56 + t119 * (-t66 * t88 + t152) - t3 * (((((t70 - t48
     #) * s25 + t197) * s25 - t35) * s25 - t68) * s25 + 11.0_dp * t10) +
     # t34 * t45 - (t148 * (t148 * (26.0_dp * s26 + t144 + t207) + ((t20
     #4 * t24 + t219) * s26 + t34 * t202) * s26 + t160 + t182) + (((136.
     #0_dp * t3 + 91.0_dp * t72) * s25 + 119.0_dp * t7) * s25 + t4 * t98
     # + t216) * s25 + 189.0_dp * t10 + t119 * (t214 * t88 + t20) + s26 
     #* (t215 + (((432.0_dp * s16 + 114.0_dp * s25) * s25 + 414.0_dp * t
     #3) * s25 + 270.0_dp * t7) * s25 + 550.0_dp * t4)) * t148 - t130 * 
     #t97 * t1 * t85 - t220 * t52 * t77 - t218 * t1 * t81
      t58 = t97 * t1
      t61 = s25 * (-t119 * s25 * (s16 + s34) + t175 + t3 + t148 - t114 *
     # s34) * s56 * t95 ** 2.0_dp + t58 * s25 * t87 + t52 * t77 * t84
      t68 = t7 * t16
      t73 = s16 + s26 - s34 + s56
      t77 = 14.0_dp * s25
      t96 = 11.0_dp * s16
      t98 = (-t96 + t77) * s25 + t27
      t130 = 90.0_dp
      t131 = 64.0_dp * t3
      t132 = t16 - t3
      t136 = -t14 * s16 * t132 + (t16 * t34 - t180) * s25
      t138 = 38.0_dp * s16
      t139 = 77.0_dp * t7
      t140 = t7 * t1
      t142 = 10.0_dp * s16
      t143 = 84.0_dp * s25
      t39 = (((-t125 - t77) * s25 + t39) * s25 - t7) * s25 + t62
      t125 = (((-s16 * t47 - t206) * s25 + 32.0_dp * t3) * s25 - 27.0_dp
     # * t7) * s25 + t212
      t155 = 211.0_dp
      t156 = 66.0_dp * s25
      t47 = t47 * t3
      t160 = 18.0_dp * t3
      t164 = (t210 + t50) * s25 + t160
      t166 = t164 * s26
      t168 = 30.0_dp * s16
      t169 = 38.0_dp * t3
      t180 = t50 + t56
      t187 = 7.0_dp * s25
      t189 = 28.0_dp * t16
      t192 = 55.0_dp * t3
      t201 = 34.0_dp * t7
      t202 = s16 * (((t192 - t189) * s25 - t201) * s25 + t23)
      t204 = t24 * t7
      t212 = (((-39.0_dp * s16 + t21) * s25 + t223) * s25 - 47.0_dp * t7
     #) * s25 + 20.0_dp * t4
      t213 = t10 * s25
      t214 = t136 * t88
      t215 = -264.0_dp * s25
      t219 = (((460.0_dp * s16 + t215) * s25 + 492.0_dp * t3) * s25 + 60
     #8.0_dp * t7) * s25 - 232.0_dp * t4
      t223 = t130 * t7
      t224 = ((t234 - t207) * s25 - t169) * s25
      t225 = t223 + t224
      t226 = 10.0_dp * t225
      t230 = 278.0_dp * t10
      t237 = ((((-t93 + t158) * s25 + 372.0_dp * t3) * s25 + 335.0_dp * 
     #t7) * s25 + 54.0_dp * t4) * s25
      t243 = 76.0_dp * t7
      t244 = ((t137 + t146) * s25 + t176) * s25
      t245 = t243 + t244
      t246 = (36.0_dp * s25 + t207) * s25 + t47
      t247 = 20.0_dp * t246
      t248 = t245 * t13
      t249 = 168.0_dp * t3
      t250 = (179.0_dp * s16 + 156.0_dp * s25) * s25
      t251 = t249 + t250
      t252 = 55.0_dp * s16
      t253 = t234 + t252
      t254 = t251 * t119
      t255 = 36.0_dp * s26
      t256 = t98 * t88
      t257 = 10.0_dp * t202
      t258 = 74.0_dp * t4
      t259 = t136 * s26
      t260 = 30.0_dp * t259
      t242 = t242 * t164
      t261 = 240.0_dp * s25
      t262 = 81.0_dp * s16
      t263 = 600.0_dp
      t264 = 78.0_dp * s16
      t265 = (((63.0_dp * s16 + 300.0_dp * s25) * s25 - 738.0_dp * t3) *
     # s25 - 684.0_dp * t7) * s25 - 780.0_dp * t4
      t266 = 325.0_dp * t3
      t267 = (t205 + t25) * s25
      t268 = 558.0_dp * t7
      t269 = 360.0_dp * t4
      t270 = 772.0_dp * t10
      t271 = ((-399.0_dp * s16 - 480.0_dp * s25) * s25 + t113) * s25 + 2
     #10.0_dp * t7
      t272 = (459.0_dp * s16 + 420.0_dp * s25) * s25 + 192.0_dp * t3
      t273 = 528.0_dp * s16 + 504.0_dp * s25
      t274 = t178 * s25
      t275 = 48.0_dp * s16
      t276 = 45.0_dp * s25
      t277 = 68.0_dp * t3
      t49 = (-s16 * t49 - t276) * s25
      t278 = t119 * s16
      t245 = (((((t255 + t278 + t25) * s34 + (-162.0_dp * s16 - 192.0_dp
     # * s25 - 108.0_dp * s26) * s26 + t277 + t49) * s34 + ((t275 + t274
     #) * s25 - 185.0_dp * t3) * s25 - 324.0_dp * t7 + 120.0_dp * t88 + 
     #s26 * (s26 * t273 + t272)) * s34 + s16 * (((-s25 * t134 + 265.0_dp
     # * s16) * s25 + 490.0_dp * t3) * s25 + 670.0_dp * t7) + ((-t14 * t
     #251 - t147 * t253) * s26 + t271) * s26 - 45.0_dp * t17 - 45.0_dp *
     # t85) * s34 + (((-t266 + t267) * s25 - t268) * s25 - t269) * s25 +
     # (((t233 * t180 + t247) * s26 + t124 * t245) * s26 + t265) * s26 -
     # t270) * s34 + s16 * (((((255.0_dp * s16 - 51.0_dp * s25) * s25 + 
     #300.0_dp * t3) * s25 + 431.0_dp * t7) * s25 - 103.0_dp * t4) * s25
     # + 516.0_dp * t10) - t24 * (-(t230 + t237) * s26 + t19) + t84 * ((
     #-t226 - t242) * s26 + 0.3e1_dp / 0.2e1_dp * t219)
      t279 = s16 * t125
      t280 = 438.0_dp
      t281 = 56.0_dp * t4
      t282 = 96.0_dp * t10
      t283 = 309.0_dp
      t284 = 96.0_dp * s25
      t285 = 382.0_dp * t7
      t286 = 241.0_dp
      t134 = (((-s16 * t134 - t276) * s25 + 265.0_dp * t3) * s25 + 490.0
     #_dp * t7) * s25 + 670.0_dp * t4
      t276 = 132.0_dp
      t287 = 35.0_dp * s25
      t276 = (((s16 * t276 - t287) * s25 + 423.0_dp * t3) * s25 + 570.0_
     #dp * t7) * s25 + 630.0_dp * t4
      t102 = ((t102 + t104) * s25 - 370.0_dp * t3) * s25 - 648.0_dp * t7
      t104 = 419.0_dp * t3
      t288 = (-180.0_dp * s16 + t126) * s25
      t289 = 448.0_dp * t7
      t49 = t277 + t49
      t277 = -t284 - t262
      t203 = (t126 + t203) * s25 + 98.0_dp * t3
      t290 = s25 + s26
      t51 = ((t148 * (s34 * (-t220 + t22) + (t194 - t275 - t193) * s16 +
     # 36.0_dp * s26 * t290) + s16 * ((-t104 + t288) * s25 - t289) + t17
     #8 * t85 + ((0.2e1_dp / 0.3e1_dp * s26 * t273 + t272) * s26 + t102)
     # * s26) * s34 + s16 * t276 + ((((-t15 * t253 - t179) * s26 - t254)
     # * s26 + t271) * s26 + t134 * t34) * s26) * s34 + s16 * (((((-56.0
     #_dp * s16 + t51) * s25 - 299.0_dp * t3) * s25 - t285) * s25 - 345.
     #0_dp * t4) * s25 - 560.0_dp * t10) + (-t34 * ((((t266 - t267) * s2
     #5 + t268) * s25 + t269) * s25 + t270) + t84 * (t147 * t246 + t248)
     #) * s26 + t84 * (t88 * t124 * t180 + t265)
      t28 = s34 * ((t51 * s34 + t41 * (-t14 * (t164 * t85 + t195) + t3 *
     # ((((184.0_dp * s16 + t232) * s25 + 108.0_dp * t3) * s25 + 14.0_dp
     # * t7) * s25 + 308.0_dp * t4) + (((((s16 * t286 - t261) * s25 + 22
     #6.0_dp * t3) * s25 - 342.0_dp * t7) * s26 + (((700.0_dp * s16 + t2
     #15) * s25 + 251.0_dp * t3) * s25 + t285) * s25 + t141 * t4) * s26 
     #+ ((((s16 * t283 - t284) * s25 + 416.0_dp * t3) * s25 + 754.0_dp *
     # t7) * s25 - 220.0_dp * t4) * s25 + 724.0_dp * t10) * s26)) * s34 
     #- t190 * (s16 * (((((-s16 * t30 + 102.0_dp * s25) * s25 + 87.0_dp 
     #* t3) * s25 - t172) * s25 - t123) * s25 + t282) + (((((s16 * t280 
     #- t274) * s25 - 284.0_dp * t3) * s25 + 254.0_dp * t7) * s25 - t281
     #) * s26 + t119 * (((((s16 * t141 - t120) * s25 - 99.0_dp * t3) * s
     #25 + t170 * t7) * s25 - 80.0_dp * t4) * s25 + t28)) * s26 - t124 *
     # (t214 + t19)))
      t51 = t140 * ((((-t143 - t142) * s25 + 26.0_dp * t3) * s25 - 51.0_
     #dp * t7) * s25 + 29.0_dp * t4)
      t123 = 26.0_dp * s25
      t82 = -(-t148 * t164 - t188 * t24 + t34 * s34 * (t148 * t180 + t13
     #6) + 14.0_dp * t72 * t67 - 25.0_dp * t229) * s56 - s16 * (((-t192 
     #+ t189) * s25 + t201) * s25 - t23 - t220 * t1 * t98) - t151 * (t12
     #4 * t148 + t34 * (36.0_dp * s16 * s26 + (t207 + t233) * s25 + 36.0
     #_dp * t16 + t47)) - ((-t148 * (t234 + t179 + t252) - t220 * t164 -
     # t223 - t224) * s34 + s16 * (((-t103 + t82) * s25 - 94.0_dp * t3) 
     #* s25 + t208) + t124 * (t17 + t259)) * s34
      t69 = t82 * s56 - s16 * (t10 * t119 + t16 * (((t262 + t77) * s25 -
     # 111.0_dp * t3) * s25 + 17.0_dp * t7) - 15.0_dp * t1 * t98 * t84) 
     #+ t15 * (s26 * t202 + t213) - t151 * (t34 * ((15.0_dp * t180 * s26
     # + t15 * t246) * s26 + t243 + t244) - 18.0_dp * t151) + ((t148 * (
     #-(t205 + t143 + t240) * s34 + (275.0_dp * s16 + 45.0_dp * s26 + t2
     #61) * s26 + t249 + t250) - (((s16 * t161 - t156) * s25 + t171 * t3
     #) * s25 + t217) * s25 - (-t15 * t225 - t242) * s26 + t185 * t4) * 
     #s34 - s25 * ((((t120 - t69) * s25 - 114.0_dp * t3) * s25 + t211) *
     # s25 - t258) - (10.0_dp * t212 + t260) * s26) * s34
      t69 = t69 * s56 - t153 * (s34 * t124 - t34 * (t93 + t255 + t26)) -
     # 26.0_dp * t8 - ((((((153.0_dp * s16 + 140.0_dp * s25) * s25 + (33
     #6.0_dp * s25 + 120.0_dp * s26 + t174) * s26 + t131) * s34 + (-t178
     # * t84 - t254) * s26 + t16 * (-133.0_dp * s16 - 160.0_dp * s25) - 
     #10.0_dp * t253 * t84 + 10.0_dp * t11 + t133) * s34 + (((t91 + t146
     #) * s25 - 246.0_dp * t3) * s25 - 228.0_dp * t7) * s25 + ((t239 * t
     #180 + t247) * s26 + t248) * s26 - 260.0_dp * t4) * s34 + ((-t226 -
     # 20.0_dp * t166) * s26 + t219) * s26 + t230 + t237) * s34 - t114 *
     # t20 + t119 * t19 - 256.0_dp * t3 * t17 + 89.0_dp * t7 * t18 + 155
     #.0_dp * t213 - 214.0_dp * t79 - 136.0_dp * t9 + t194 * s25 * ((((-
     #t145 + t21) * s25 - t107) * s25 + t204) * s25 - 37.0_dp * t4) + 40
     #.0_dp * t214 + 20.0_dp * t212 * t84) * s34 + 20.0_dp * t256 * t52 
     #+ 80.0_dp * t184 - 56.0_dp * t3 * t20 - 44.0_dp * t7 * t17 - t228 
     #* t16 + 164.0_dp * t4 * t18 - t238 * t52 * t39 + t257 * t84
      t8 = (t69 * s56 - t14 * t52 * s26 * (s26 * t39 + t279) + t150 * t2
     #4 - 29.0_dp * t5 - (t245 * s34 - t41 * ((((((t264 - t120) * s25 + 
     #t249) * s25 - 111.0_dp * t7) * s25 + 218.0_dp * t4) * s25 - 239.0_
     #dp * t10) * s25 + ((-10.0_dp * (((-84.0_dp * s16 + t120) * s25 + 4
     #6.0_dp * t3) * s25 - 67.0_dp * t7) * s25 - 10.0_dp * t4 * t83 - t2
     #60) * s26 + ((((s16 * t263 - t137) * s25 - 156.0_dp * t3) * s25 + 
     #424.0_dp * t7) * s25 - 226.0_dp * t4) * s25 + 220.0_dp * t10) * s2
     #6 + 188.0_dp * t9)) * s34 + t216 * t52 * t98 + 80.0_dp * t157 + 36
     #.0_dp * t10 * t18 + t258 * t17 - 77.0_dp * t9 * t16 - 84.0_dp * t7
     # * t20 + t257 * t88) * s56 - t14 * t52 * t84 * (t279 - t256) - t34
     # * (t153 * (s16 * t203 - 36.0_dp * t88 + s26 * (s26 * t277 + t49))
     # + t51 * s26) - t28 - 13.0_dp * t6 - t52 * t119 * t39 * t88 - t178
     # * t9 * t18 - t281 * t20 + t282 * t17 + 7.0_dp * t8 * t16 + t123 *
     # t5 + t159 * t202
      t28 = t41 * t190 * ((t16 * (t16 * (t138 - t117) + t139) + t63) * s
     #25 + t34 * s26 * (-t84 * t136 + ((((s16 * t185 - t21) * s25 - 72.0
     #_dp * t3) * s25 + 75.0_dp * t7) * s25 - t129) * s25 + 24.0_dp * t1
     #0) + t84 * ((((s16 * t130 - t120) * s25 - t131) * s25 + t208) * s2
     #5 - t62) + 16.0_dp * t9 - 56.0_dp * t229 * t67) * s34
      t33 = t190 * ((((((-t55 + t241) * s25 + 27.0_dp * t3) * s25 + 93.0
     #_dp * t7) * s25 - t33) * s25 + 63.0_dp * t10) * s25 + 56.0_dp * t9
     # + (((-t166 + ((s16 * t111 - t234) * s25 + t47) * s25 - 54.0_dp * 
     #t7) * s26 + (((s16 * t155 - t156) * s25 + t44) * s25 + t186) * s25
     # + 68.0_dp * t4) * s26 + ((((s16 * t110 - t93) * s25 - t92) * s25 
     #+ 238.0_dp * t7) * s25 - t162) * s25 + 196.0_dp * t10) * s26) * t1
     #48
      t47 = t41 * (s16 * (((((t165 - t187) * s25 + 97.0_dp * t3) * s25 +
     # t208) * s25 + 175.0_dp * t4) * s25 + 112.0_dp * t10) + t34 * s26 
     #* (-t85 * t180 + ((((-t75 + t50) * s25 + 101.0_dp * t3) * s25 + 17
     #1.0_dp * t7) * s25 + t71) * s25 + 224.0_dp * t10) + t84 * ((((-144
     #.0_dp * s25 + t168) * s25 + t169) * s25 - 64.0_dp * t7) * s26 + ((
     #(s16 * t171 - t91) * s25 + 216.0_dp * t3) * s25 + 190.0_dp * t7) *
     # s25 + 324.0_dp * t4) - t119 * (t3 * t83 - t94) * t85) * t151
      t6 = s25 * t6 + t8 * s56 - t149 * (-t24 * t86 + t3 * ((((t264 + t1
     #17) * s25 + 145.0_dp * t3) * s25 + t182) * s25 + t90) + (s16 * t27
     #6 + (((-t253 * s26 - t251) * s26 + t271 / 3.0_dp) * s26 + t134) * 
     #s26) * s26) - t152 * (t124 * t88 + t3 * (-16.0_dp * s16 - 13.0_dp 
     #* s25) + s26 * (-s16 * (t275 + t193) + s26 * (t278 + t25))) - t153
     # * (t84 * (0.2e1_dp / 0.3e1_dp * s26 * t277 + t49) - 18.0_dp * t85
     # + t281 + 29.0_dp * t229 + t139 * s25 + t114 * t203 * s26) - t154 
     #* (t124 * t87 - 112.0_dp * t7 * t67 + (-s16 * ((t104 - t288) * s25
     # + t289) + ((s26 * t273 / 6.0_dp + t272 / 3.0_dp) * s26 + t102 / 2
     #.0_dp) * s26) * s26 - t11 * (189.0_dp * t3 + t222)) + t87 * t202 +
     # t28 - t4 * t1 * (t16 * ((-56.0_dp * s25 + t145) * s25 - t44) + 13
     #.0_dp * t140) * s26 - t52 * t39 * t85 - t51 * t84 + t52 * t98 * t8
     #6 - t33 + t47 - (t114 - s26) * (t183 + s16) * t150
      t8 = t24 * t4
      t25 = (-t210 + t76) * t18
      t28 = -t119 * t16 * t1 + t3 * t97
      t33 = t28 * t84
      t39 = -t117 + t22
      t47 = 15.0_dp * t221
      t35 = s25 * ((25.0_dp * t16 + t121) * s25 + t35)
      t22 = ((24.0_dp * s25 + t22) * s25 + t92) * s25
      t49 = t128 + t22
      t21 = (t21 + s16) * s25 + t27
      t27 = 30.0_dp * t7
      t51 = ((t163 + t80) * s25 + 39.0_dp * t3) * s25
      t55 = s16 * t83
      t62 = (t123 + t55) * s25 + t160
      t37 = t37 + t117
      t63 = t17 + t4
      t69 = t187 + t50
      t71 = 16.0_dp * s25
      t12 = t14 * t12
      t79 = t183 * t21
      t82 = t1 ** 2.0_dp
      t83 = s26 * t82
      t86 = t7 * s25
      t87 = t24 * (t16 * (t82 * t84 - t229) + t153) * s25 + t34 * t18 * 
     #(t83 * t100 + t86) + ((s25 * ((t34 * t49 + t79) * s26 + t35 + t47)
     # + (-s25 * (t27 + t51 + s26 * (t24 * t37 * s26 + t34 * t62)) + (s2
     #5 * ((t80 + t287) * s25 + (t206 + t43 + t183) * s26 + t113) - s25 
     #* (t220 + t71 + t158) * s34) * s34) * s34) * s34 + (-t40 * (-t25 +
     # t8 - t12) + t68 * t14) * s25 - t13 * t17 * t132 - t24 * s25 * (t4
     # * (s16 - s25) + t33)) * s34 + t18 * t17
      t82 = t82 * t18
      t90 = t67 * s16
      t91 = t83 * t17 * t115
      t71 = s25 * ((t183 + t71) * s26 + 15.0_dp * t16 + 15.0_dp * t167 +
     # 24.0_dp * t72) * t154
      t48 = s25 * (((t34 * t69 + s26) * s26 + t15 * ((t118 + t187) * s25
     # + t57)) * s26 + 20.0_dp * t74 + t72 * (37.0_dp * s25 + t48)) * t1
     #49
      t27 = s25 * (15.0_dp * t63 + ((s26 * t37 + t62) * s26 + t27 + t51)
     # * s26 + t72 * ((t209 + t206) * s25 + t197)) * t151
      t51 = -t92 * t20 + t153 * (t2 * t34 + s26) * t24
      t8 = t51 * s25 - s25 * t152 + (-((t2 * (-t13 * t199 + t38 * t7) + 
     #t33) * s26 + t19 + t9) * s25 + t84 * (t14 * t229 * t2 + s25 * (t25
     # - t8)) + s16 * t18 * (t16 * t39 + t204)) * s34 + t87 * s56 + t148
     # * (t14 * s25 * ((t17 + t4) * s25 + t10) + s25 * (((s26 * t21 + t4
     #9) * s26 + t35 + t47) * s26 + t68) + s16 * t17 * (t120 + t114)) + 
     #t95 * (t19 * t34 - t24 * (t16 * (t46 - t83) + t154) * s25 + s25 * 
     #(((((t196 + t183 + t77) * s34 - s16 * (t181 + t191) - (s26 * t124 
     #+ t123 + t55) * s25) * s34 + t128 + t22 + t79) * s34 - t24 * (s26 
     #* t28 + t4) + t12 + t25) * s34 + t68)) + s25 * ((((-t37 + s34) * s
     #34 + t21) * s34 - t28) * s34 + t173 - t114 * t18) * t108 + t82 * t
     #88 + t90 * t20 + t91 + t82 * t100 * t84 - t71 + t48 - t27
      t9 = s16 + s25 + s26 - s34 + s56
      t12 = t17 * t2
      t19 = t46 * t2
      t20 = -t177 * s25 + t24 * t63 + t19
      t21 = 15.0_dp * t4
      t22 = (t16 * (t75 + t76) - t204) * s25
      t25 = -t278 + t56
      t27 = t124 * t74 + t72 * t25
      t28 = -t13 * t16 - t57
      t33 = t178 * t7
      t35 = ((53.0_dp * s16 + 42.0_dp * s25) * s25 + t122) * s25
      t37 = t116 * t67 + t72 * t13
      t38 = t58 * t115 * t84
      t13 = (((-7.0_dp * t3 + t175) * s25 + t34 * (t90 - t151) + s34 * (
     #s34 * t69 + t28)) * s56 + t14 * (t151 * (-t200 - s26 + s34) + t17 
     #+ t4) + t148 * ((21.0_dp * s25 + t181) * s26 + 36.0_dp * t67 + 32.
     #0_dp * t72) + t34 * (-(s16 * (s25 * t25 + t118 * s26) + t124 * (t1
     #6 * t290 + t7)) * s34 + t19) - t236 + t235 * t97 * t115) * s56 + t
     #119 * (s26 * t20 + t18 * t3) - t14 * t198 + t148 * ((t13 * t37 + t
     #183 * t69) * s26 + t33 + t35) + t24 * (t149 * (t142 + t75 + t238) 
     #+ t12 + t38) - t34 * s34 * (t148 * ((t168 + t287) * s16 + 24.0_dp 
     #* t16 + t218 + 24.0_dp * t54) + t21 + s26 * (-0.3e1_dp / 0.2e1_dp 
     #* s26 * t28 + t34 * t27) + t22) - t86 * (t96 + t56)
      t19 = 0.1e1_dp / s16
      t25 = 0.1e1_dp / t42 ** 2.0_dp
      t36 = 0.1e1_dp / t36 ** 2.0_dp
      t31 = 0.1e1_dp / t31 ** 2.0_dp
      t6 = (t34 * t3 * (t1 * t125 * t88 + t5) - t6 + t10 * t16 * (((-t64
     # + t77) * s25 + 29.0_dp * t3) * s25 - 13.0_dp * t7)) * MB1011(1) +
     # t29
      t5 = -(18.0_dp * t199 * t4 + t24 * t61 + t34 * (-t152 * (t40 + s16
     #) + t5) - s34 * ((t65 * s34 - t41 * (((((34.0_dp * s16 + t117) * s
     #25 + t127) * s25 + t78) * s25 + t24 * t85 + t129) * s25 + t89 + ((
     #((t112 + t165) * s25 + t113) * s26 + ((98.0_dp * s16 + 116.0_dp * 
     #s25) * s25 + 79.0_dp * t3) * s25 + t135) * s26 + (((100.0_dp * s16
     # + t112) * s25 + 109.0_dp * t3) * s25 + 56.0_dp * t7) * s25 + t110
     # * t4) * s26)) * s34 + t190 * ((((t106 + t44) * s25 - 26.0_dp * t7
     #) * s25 + t23 + t109 * t14) * s25 + 14.0_dp * t10 + s26 * ((((-t19
     #6 + t99) * s25 + t3) * s25 + t105) * s26 - 40.0_dp * t11 * t1 + t1
     #8 * (s25 * t30 - t43) + 42.0_dp * t4))) - t45 * s56 - t157 + t3 * 
     #t1 * ((((t118 - t70) * s25 - 15.0_dp * t3) * s25 - t177) * s25 + t
     #59) * s26 + t1 * t81 * t88 + t66 * t85 + t68 * (t116 * t18 - t78))
     # * MB1001(0) + t6 * t19
      t6 = t9 ** 2.0_dp
      t6 = t19 * (t25 * (-t8 * MB1101(0) * t31 + t73 * t9 * t6 * MB1111D
     #2(0)) - t73 ** 2.0_dp * MB1011(0) * t36)
      result = -t25 * t31 * (t60 * I300s16(1) + t5 * t36 - t9 * s25 * (-
     #t15 * t68 * t1 - t24 * t52 * t17 - t34 * (-t1 * t10 + t154 * (t50 
     #+ t183 + t56 - s34) - t20 * t84) + (t12 * t24 + t3 * ((-s25 * t39 
     #- t176) * s25 + t211) + t38) * s26 + (((((27.0_dp * s26 + t144 + t
     #138) * s25 + 30.0_dp * t167 + t231) * s34 - ((t144 + t80) * s25 + 
     #t3 * t30) * s25 - ((24.0_dp * t2 + t40) * s26 + (70.0_dp * s16 + t
     #234) * s25 + t178 * t3) * s26 - t208) * s34 + (((t126 + t26) * s25
     # + t227) * s25 + t53) * s25 + 30.0_dp * t4 + s26 * (t69 * t84 + t3
     #3 + t35) + t119 * t37 * t84) * s34 + s16 * (t16 * (-t106 + t32) - 
     #t101) + t28 * t88 - t34 * (s25 * t63 + (s26 * t27 + t21 + t22) * s
     #26)) * s34 + t13 * s56) * MB1110(1) * t19) / 12.0_dp - t6 / 6.0_dp

           intHs16s25s26s34s56x1114D6eps0 = result
       end function intHs16s25s26s34s56x1114D6eps0

       function intHs16s25s26s34s56x1114D6eps1()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1114D6eps1
           complex(dp) ::  t1,t2,t3

           complex(dp) :: result

      t1 = s16 + s26 + s56
      t2 = s16 + s26 - s34 + s56
      t3 = 0.1e1_dp / t1 ** 2.0_dp
      result = 0.1e1_dp / s16 * t3 * (s34 / 6.0_dp - t1 / 4.0_dp + t2 **
     # 2.0_dp * MB1011(1) / 6.0_dp)

           intHs16s25s26s34s56x1114D6eps1 = result
       end function intHs16s25s26s34s56x1114D6eps1

       function intHs16s25s26s34s56x1120D2eps0()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1120D2eps0
           complex(dp) ::  t1

           complex(dp) :: result

      t1 = s16 + s25 + s26 - s34 + s56
      t1 = 0.1e1_dp / t1
      result = 2.0_dp * s25 * MB1110(1) * t1 - (-s25 * MB1110(0) + (s16 
     #+ s26 + s56) * MB1011(1) - MB1001(0)) * t1

           intHs16s25s26s34s56x1120D2eps0 = result
       end function intHs16s25s26s34s56x1120D2eps0

       function intHs16s25s26s34s56x1120D2eps1()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1120D2eps1
           complex(dp) ::  t1

           complex(dp) :: result

      t1 = s16 + s25 + s26 - s34 + s56
      t1 = 0.1e1_dp / t1
      result = t1 * (s25 * MB1110(1) + 1.0_dp)

           intHs16s25s26s34s56x1120D2eps1 = result
       end function intHs16s25s26s34s56x1120D2eps1

       function intHs16s25s26s34s56x1121D2eps0()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1121D2eps0
           complex(dp) ::  t1,t2,t3

           complex(dp) :: result

      t1 = s16 + s26 - s34 + s56
      t2 = s16 - s25
      t1 = 0.1e1_dp / t1
      t3 = 0.1e1_dp / s16
      result = 2.0_dp * s25 * s34 * MB1101(0) * t3 * t1 ** 2.0_dp + t1 *
     # ((s25 * MB1110(0) + s34 * MB1011(0) - (-s16 * t2 - s26 * t2 + (2.
     #0_dp * s16 + s25 + s26 - s34) * s34 - (s16 - s25 - s34) * s56) * M
     #B1111D2(0) * t1) * t3 - I300s16(0))

           intHs16s25s26s34s56x1121D2eps0 = result
       end function intHs16s25s26s34s56x1121D2eps0

       function intHs16s25s26s34s56x1121D2eps1()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1121D2eps1
           complex(dp) ::  t1

           complex(dp) :: result

      t1 = s16 + s26 - s34 + s56
      t1 = 0.1e1_dp / t1
      result = t1 * ((s25 * MB1110(1) + s34 * MB1011(1)) / s16 - I300s16
     #(1))

           intHs16s25s26s34s56x1121D2eps1 = result
       end function intHs16s25s26s34s56x1121D2eps1

       function intHs16s25s26s34s56x1121D4eps0()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1121D4eps0
           complex(dp) ::  t1,t2,t3,t4,t5

           complex(dp) :: result

      t1 = s16 + s26 - s34 + s56
      t2 = 2.0_dp
      t3 = s16 * s25 + (-s34 + s16 + s25 + s26) * s26 + (s26 * t2 + s16 
     #+ s25 - s34 + s56) * s56
      t4 = -s16 + s25
      t5 = 0.1e1_dp / t3
      t3 = 0.1e1_dp / t3
      result = -t1 * t3 - t5 * (t1 * (s16 * I300s16(1) - s25 * MB1110(1)
     # + t1 * MB1011(1) - MB1001(0)) - (s16 * t4 + s26 * t4 + (s16 * t2 
     #+ s25 + s26 - s34) * s34 + (-s16 + s25 + s34) * s56) * MB1111D2(0)
     # + s25 * s34 * MB1101(0)) / 2.0_dp

           intHs16s25s26s34s56x1121D4eps0 = result
       end function intHs16s25s26s34s56x1121D4eps0

       function intHs16s25s26s34s56x1121D4eps1()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1121D4eps1

           complex(dp) :: result

      result = 0.0_dp

           intHs16s25s26s34s56x1121D4eps1 = result
       end function intHs16s25s26s34s56x1121D4eps1

       function intHs16s25s26s34s56x1122D4eps0()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1122D4eps0
           complex(dp) ::  t1,t10,t11,t12,t13,t14,t15,t16,t2,t3,t4,t5
           complex(dp) ::  t6,t7,t8,t9

           complex(dp) :: result

      t1 = s16 + s25 + s26 - s34 + s56
      t2 = 2.0_dp
      t3 = t2 * s26
      t4 = s16 * s25 + (-s34 + s16 + s25 + s26) * s26 + (t3 + s16 + s25 
     #- s34 + s56) * s56
      t5 = s16 + s26 - s34 + s56
      t6 = s16 + s26 + s56
      t7 = s25 + s26
      t8 = t2 * s16
      t9 = s16 + s26
      t10 = s16 * t9
      t11 = s26 ** 2.0_dp
      t12 = s16 ** 2.0_dp
      t13 = (t11 + t12) * s16
      t14 = 3.0_dp * s16
      t4 = 0.1e1_dp / t4
      t6 = 0.1e1_dp / t6
      t15 = 0.1e1_dp / s16
      t16 = 0.1e1_dp / t5
      result = t2 * t1 * t15 * t4 - t6 * (t15 * (-s34 * MB1011(0) - (t3 
     #* t12 + ((s16 - s26) * s34 - (t8 - t7) * t9) * s34 + ((s16 + s34) 
     #* s56 + t10 * t2 - s34 * (s34 - t3 + s16 - s25)) * s56 + t13) * MB
     #1011(1) * t4) + t5 * MB1001(0) * t4) + (t2 * (s26 * t12 + ((-t8 - 
     #t7 + s34) * s34 + t10) * s56) + ((t2 * t7 - s34 + t14) * s34 - t2 
     #* s25 * t9 - (4.0_dp * s26 + t14) * s16 - s25 ** 2.0_dp - t11) * s
     #34 + (s16 - s34) * s56 ** 2.0_dp + t13) * MB1111D2(0) * t15 * t16 
     #* t4 - t4 * (s25 * t15 * (-s34 * MB1101(0) * t16 + MB1110(1)) - I3
     #00s16(1)) * t1

           intHs16s25s26s34s56x1122D4eps0 = result
       end function intHs16s25s26s34s56x1122D4eps0

       function intHs16s25s26s34s56x1122D4eps1()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1122D4eps1
           complex(dp) ::  t1

           complex(dp) :: result

      t1 = s16 + s26 + s56
      t1 = 0.1e1_dp / t1
      result = 0.1e1_dp / s16 * t1 * (s34 * MB1011(1) + 1.0_dp)

           intHs16s25s26s34s56x1122D4eps1 = result
       end function intHs16s25s26s34s56x1122D4eps1

       function intHs16s25s26s34s56x1123D6eps0()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1123D6eps0
           complex(dp) ::  t1,t10,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109
           complex(dp) ::  t11,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t12
           complex(dp) ::  t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t13,t130
           complex(dp) ::  t131,t132,t133,t134,t135,t136,t137,t138,t139,t14,t140,t141
           complex(dp) ::  t142,t143,t144,t145,t146,t147,t148,t149,t15,t150,t151,t152
           complex(dp) ::  t153,t154,t155,t156,t157,t158,t159,t16,t160,t161,t162,t163
           complex(dp) ::  t164,t165,t166,t167,t168,t169,t17,t170,t171,t172,t173,t174
           complex(dp) ::  t175,t176,t177,t178,t179,t18,t19,t2,t20,t21,t22,t23
           complex(dp) ::  t24,t25,t26,t27,t28,t29,t3,t30,t31,t32,t33,t34
           complex(dp) ::  t35,t36,t37,t38,t39,t4,t40,t41,t42,t43,t44,t45
           complex(dp) ::  t46,t47,t48,t49,t5,t50,t51,t52,t53,t54,t55,t56
           complex(dp) ::  t57,t58,t59,t6,t60,t61,t62,t63,t64,t65,t66,t67
           complex(dp) ::  t68,t69,t7,t70,t71,t72,t73,t74,t75,t76,t77,t78
           complex(dp) ::  t79,t8,t80,t81,t82,t83,t84,t85,t86,t87,t88,t89
           complex(dp) ::  t9,t90,t91,t92,t93,t94,t95,t96,t97,t98,t99

           complex(dp) :: result

      t1 = s16 - s25
      t2 = 2.0_dp
      t3 = t2 * s16
      t4 = s25 + t3
      t5 = s16 + s25
      t6 = s25 * t5
      t7 = 3.0_dp
      t8 = s16 ** 2.0_dp
      t9 = t8 ** 2.0_dp
      t10 = s16 * t8
      t11 = t10 * t9
      t12 = t8 * t9
      t13 = s16 * t9
      t14 = t7 * t8
      t15 = 4.0_dp
      t16 = t15 * s25
      t17 = 9.0_dp * t8
      t18 = s25 ** 2.0_dp
      t19 = t18 ** 2.0_dp
      t20 = s25 * t18
      t21 = t18 + t8
      t22 = s25 * t21
      t23 = 9.0_dp * t22
      t24 = (11.0_dp * t18 + 24.0_dp * t8) * s16
      t25 = 11.0_dp * s16
      t26 = 17.0_dp * s25
      t27 = (t25 + t26) * s25
      t28 = 27.0_dp * t8
      t29 = t27 + t28
      t30 = s16 * s25
      t31 = t30 * t21
      t32 = 25.0_dp
      t33 = 7.0_dp
      t34 = 5.0_dp
      t35 = t15 * s26
      t36 = t18 * t34
      t37 = t33 * t8
      t38 = 42.0_dp * s25
      t39 = 66.0_dp * s16
      t40 = (t38 + t39) * s25
      t41 = 72.0_dp * t8
      t42 = 27.0_dp * s16
      t43 = 22.0_dp * s25
      t44 = s26 ** 2.0_dp
      t45 = t44 ** 2.0_dp
      t46 = s26 * t44
      t47 = 20.0_dp * s25
      t48 = t2 * t46
      t49 = 48.0_dp * s16
      t50 = 9.0_dp * s26
      t51 = 57.0_dp * s16
      t52 = 50.0_dp * t8
      t53 = 12.0_dp * s26
      t54 = t32 * s16
      t55 = 6.0_dp
      t56 = t55 * s26
      t57 = 44.0_dp
      t58 = 24.0_dp * s16
      t59 = s34 ** 2.0_dp
      t60 = t59 ** 2.0_dp
      t61 = s34 * t60
      t62 = s34 * t59
      t63 = t1 * ((-t3 - t16) * s25 + t17)
      t64 = t60 + t9
      t65 = 10.0_dp
      t66 = t1 ** 2.0_dp
      t67 = 8.0_dp * s16
      t68 = t65 * t8
      t69 = s26 * t5
      t70 = -t62 + t10
      t71 = 11.0_dp * t8
      t72 = s34 * t5
      t73 = -t2 * (t30 + t72) + t59 + t18 + t8
      t74 = s56 ** 2.0_dp
      t75 = t1 * t44
      t76 = t1 * s26
      t77 = t7 * s16
      t78 = (t77 + t16) * s25
      t23 = t2 * (t73 * s56 * t74 + s16 * t19 + t76 * (t4 * (-t6 + t14) 
     #+ t75)) + t34 * (-t61 + t13) + ((((t47 + t53 + t54) * s34 - (30.0_
     #dp * s25 + t51) * s25 - (38.0_dp * s25 + t49 + t50) * s26 - t52) *
     # s34 + ((41.0_dp * s16 + t47) * s25 + 51.0_dp * t8) * s25 + 50.0_d
     #p * t10 + s26 * ((t42 + t43) * s26 + t40 + t41) + t48) * s34 - ((t
     #35 * t5 + t29) * s26 + t2 * (t23 + t24)) * s26 - t18 * (t36 + t37)
     # - t32 * t9 - 11.0_dp * t31) * s34 + ((s25 * ((-t3 + t16) * s25 - 
     #t71) - s34 * (-(t42 + t43 + t56) * s34 + t27 + t28 + 12.0_dp * t69
     #) + 9.0_dp * t70 + t56 * t66) * s56 - t2 * (((19.0_dp * s25 + t50 
     #+ t58) * t59 + (t56 * t5 + t29) * s26 + t23 + t24) * s34 - t19 - t
     #63 * s26) + ((s25 * t57 + 54.0_dp * s16 + t56) * s26 + t40 + t41) 
     #* t59 + 12.0_dp * t64 + t30 * ((t16 - t67) * s25 - t68) + t55 * t6
     #6 * t44) * s56 + t63 * t44 - t78 * t10
      t24 = t2 * s26
      t27 = (-s34 + t5 + s26) * s26 + (s16 + s25 - s34 + t24 + s56) * s5
     #6 + t30
      t29 = t7 * s25
      t40 = t34 * s16
      t41 = -t18 + t14
      t63 = -s25 + t3
      t79 = s16 + s26
      t80 = 11.0_dp * s25
      t81 = 13.0_dp * t8
      t82 = t15 * s16
      t83 = t34 * s25
      t84 = t82 + t83
      t85 = 19.0_dp * s16
      t86 = 18.0_dp * s25
      t87 = t33 * s16
      t88 = 15.0_dp * t8
      t89 = t65 * t10
      t90 = 9.0_dp * s25
      t91 = 32.0_dp * s16
      t92 = (t90 + t91) * s25
      t93 = 38.0_dp * t8
      t94 = s16 * t18
      t95 = t46 + t94
      t96 = t65 * s16
      t97 = 15.0_dp * s25
      t98 = 8.0_dp * s26
      t99 = 17.0_dp * s16
      t100 = t7 * s26
      t101 = s16 + t100
      t102 = t18 * t5
      t103 = t8 * (t87 - t90) + t102
      t104 = 33.0_dp * s25
      t105 = t65 * s25
      t106 = 8.0_dp * t8
      t107 = t55 * s25
      t108 = 12.0_dp * s16
      t109 = 17.0_dp * t8
      t110 = 34.0_dp * t8
      t111 = 42.0_dp * t10
      t112 = 16.0_dp
      t113 = t112 * s26
      t114 = t18 - t8
      t115 = t8 * s25
      t116 = t10 * t114
      t117 = t15 * t44
      t118 = t55 * s16
      t119 = 20.0_dp * s16
      t120 = s25 * s26
      t121 = 21.0_dp * s16
      t122 = s25 + s26
      t123 = t94 * t1
      t124 = t10 * t1
      t125 = t35 * t1
      t126 = t59 + t8
      t127 = t82 * s34
      t128 = t2 * t126
      t129 = t20 + t10
      t130 = t79 ** 2.0_dp
      t131 = (((-s34 * t34 + s25 - t77) * s25 - t127 + t128) * s56 + t33
     # * t70 - (-(t98 + t121 + t26) * s34 + s16 * (t113 + t121) + (20.0_
     #dp * s26 + t80 + t67) * s25) * s34 + t102 - 9.0_dp * t115 + t125 *
     # t63) * s56 - (-8.0_dp * t62 + ((t58 + t104) * s25 + 63.0_dp * t8)
     # * s26 + 35.0_dp * t10 + t18 * (t83 + t119) + t55 * (t44 * t84 + t
     #115)) * s34 + t7 * (s26 * t103 + t59 * (t121 * s26 - (t122 * t33 +
     # t25) * s34 + t18 * t55 + t117 + 17.0_dp * t120 + 12.0_dp * t30 + 
     #17.0_dp * t8) - t123) + 9.0_dp * t124 + t75 * t55 * t63
      t38 = t131 * s56 - t34 * t116 + t7 * (t103 * t44 + t115 * t114 - t
     #61) + ((((t113 + t90 + t99) * s34 - (21.0_dp * s26 + t38 + t39) * 
     #s26 - t92 - t93) * s34 + ((t29 + t54) * s25 + t110) * s25 + ((63.0
     #_dp * s16 + 51.0_dp * s25 + t98) * s26 + t55 * ((t107 + t108) * s2
     #5 + t109)) * s26 + t111) * s34 - t79 * (((t87 + t105) * s25 + t106
     #) * s25 + (t35 * t84 + (t82 + t104) * s25 + 47.0_dp * t8) * s26 + 
     #23.0_dp * t10)) * s34 + t76 * (t117 * t63 + t118 * t41)
      t39 = t129 * t10
      t103 = t76 * t8
      t131 = t2 * t18
      t132 = s16 + s26 + s56
      t133 = t7 * t18
      t134 = t34 * t8
      t135 = t2 * s25
      t136 = (-s16 - t135) * s25
      t137 = t15 * t8 + t136
      t138 = 20.0_dp * t8
      t139 = t18 * t65
      t140 = 8.0_dp * s25
      t141 = (t40 + t140) * s25
      t142 = 12.0_dp * t8
      t143 = t142 + t141
      t144 = (28.0_dp * s16 + t86) * s25
      t145 = 30.0_dp * t8
      t146 = 21.0_dp * t8
      t147 = 20.0_dp * t10
      t148 = t2 * t44
      t149 = s16 * t79
      t150 = 12.0_dp * s25
      t151 = 23.0_dp * s16
      t152 = t34 * s26
      t153 = t7 * t5
      t154 = s16 * t20
      t155 = t7 * t44
      t70 = t15 * t70
      t156 = t76 * ((-t133 + t134) * s16 - s25 * t114)
      t19 = t2 * (-t61 - t116) + ((((t152 + t140 + t96) * s34 - (t113 + 
     #t150 + t151) * s25 - t117 - 20.0_dp * t149) * s34 + ((t140 + t99) 
     #* s25 + t146) * s25 + s26 * (t144 + t145 + t44) + t147 + t148 * (t
     #118 + t83)) * s34 - t2 * (t46 * t5 + t19) - t31 * t34 - t8 * (t133
     # + t68) - s26 * (t143 * s26 + (t139 + t138) * s16 + 8.0_dp * t22))
     # * s34 + ((s56 * t73 + (t131 - t134) * s25 - (-(t100 + t105 + t108
     #) * s34 + t55 * t69 + t141 + t142) * s34 - t94 + t70 + t100 * t66)
     # * s56 - t15 * (t62 * (t24 + t16 + t40) + t115 * t5) - t2 * (((t36
     # + t68) * s16 + (t153 * s26 + t143) * s26 + t15 * t22) * s34 - t15
     #4 - t76 * t137) + t34 * t64 + t59 * ((t58 + t47 + t100) * s26 + t1
     #44 + t145) + t19 + t155 * t66) * s56 + t46 * t66 + t30 * (t20 - t1
     #0) + t75 * t137 + t156
      t31 = (s16 - t107) * s25
      t36 = t136 + t134
      t47 = (-t135 - t40) * s25 + t71
      t64 = -t135 + t77
      t66 = t112 * s25
      t68 = t33 * s25
      t105 = t118 + t68
      t113 = 31.0_dp * s16
      t116 = t32 * s25
      t117 = s25 + t82
      t136 = 29.0_dp
      t137 = t136 * s25
      t141 = 55.0_dp * s16
      t143 = (t141 + t116) * s25
      t144 = 78.0_dp * t8
      t156 = 72.0_dp * s25
      t157 = 93.0_dp * s16
      t158 = 37.0_dp * s16
      t159 = 15.0_dp * s26
      t160 = t64 * t1
      t161 = t160 * t46
      t162 = 18.0_dp * s26
      t163 = 30.0_dp * s26
      t164 = t8 * t18
      t165 = t75 * t64
      t166 = 15.0_dp * s16
      t64 = t76 * t64
      t126 = (((s34 * t33 - t135 + t40) * s25 - t126 * t7 + t118 * s34) 
     #* s56 + ((-t135 - t77) * s25 + t112 * t8) * s25 - t15 * (-(s16 * (
     #t56 + t67) + (s26 * t33 + t16 + t77) * s25) * s34 + t64) + t59 * (
     #s34 * t65 - 24.0_dp * s25 - t113 - t53) - 11.0_dp * t10) * s56 + t
     #10 * (t86 - t166) - t55 * (t165 + t154) - t7 * (t76 * t47 - t164) 
     #- 11.0_dp * t60 + (((t49 + t163 + t137) * s34 - (t156 + t157 + t16
     #2) * s26 - t143 - t144) * s34 + ((t113 + t68) * s25 + t106) * s25 
     #+ (t56 * t105 + 12.0_dp * t106 + 12.0_dp * t78) * s26 + 56.0_dp * 
     #t10) * s34
      t53 = t126 * s56 + t10 * ((t68 + t67) * s25 - t17) - t15 * (-t61 +
     # t161) + t2 * s34 * (t59 * ((t58 + t107) * s25 + (t49 + t159 + t13
     #7) * s26 + t145) + t79 * (t102 * t33 + t15 * (((s16 + t107) * s25 
     #+ t17) * s26 + t115) + t147 + t148 * t105)) - t55 * s16 * (t76 * t
     #36 + t154) - t59 * (t59 * (22.0_dp * s26 + t150 + t54) + ((t16 + t
     #158) * s25 + 52.0_dp * t8) * s25 + ((t156 + t157 + t53) * s26 + t2
     # * (t143 + t144)) * s26 + 70.0_dp * t10) - t75 * t7 * t47
      t58 = t160 * t45
      t67 = t10 * t18
      t17 = t2 * (t61 * (s16 + t24) - t39) - t7 * s16 * (t75 * t36 - t67
     #) + s25 * t13 + s34 * (((-(s16 * (t140 + t96) + (11.0_dp * s26 + t
     #150 + t54) * s26) * s34 + (12.0_dp * s25 * t117 + (t49 + t137) * s
     #26 + 12.0_dp * t134) * s26 + t65 * t95 + t8 * (23.0_dp * s25 + t11
     #9)) * s34 - t79 * (s16 * ((t150 + t121) * s25 + t138) + ((t113 + t
     #116) * s25 + t52 + t155) * s26 + t15 * (t44 * (t107 + t87) + t20))
     #) * s34 + t130 * ((t18 * t33 + t134) * s25 + t89 + (t105 * s26 + (
     #-t3 + t66) * s25 + t138) * s26 - t94)) + t53 * s56 - t58 - t1 * t4
     #7 * t46 - t103 * (t31 + t17)
      t36 = -t29 + t82
      t47 = t135 + t77
      t49 = 26.0_dp
      t52 = 8.0_dp * t18
      t53 = t49 * t10
      t54 = t8 * ((-t110 + t52) * s25 + t53)
      t68 = 14.0_dp
      t102 = t68 * t8
      t105 = (-t135 - t87) * s25 + t102
      t113 = 24.0_dp * t10
      t14 = (s16 * t68 - t135) * s25 + t14
      t116 = t14 * s26
      t119 = t34 * t10
      t83 = -t77 - t83
      t121 = t2 * t8
      t25 = t107 + t25
      t126 = -t108 - t137
      t138 = s16 * ((-t29 + t141) * s25 + 105.0_dp * t8)
      t141 = t15 * t10
      t143 = 62.0_dp
      t144 = (s16 * t57 - t107) * s25
      t145 = t152 * t14
      t147 = -45.0_dp * s25
      t152 = 75.0_dp * t8
      t156 = -t83
      t157 = (t87 - t104) * s25
      t57 = t57 * t8
      t167 = 52.0_dp * s16
      t168 = 88.0_dp * s25
      t169 = t167 + t168
      t170 = (t57 + t157) * t7
      t171 = 36.0_dp
      t172 = t171 * s26
      t173 = t7 * (s16 * (t29 + t108) + (t3 - t140 - t50) * s26) * t61
      t174 = 13.0_dp * t10
      t175 = (t42 - t26) * s25 + t88
      t176 = t156 * t65 * s26
      t177 = s16 * s26
      t178 = ((s16 * t49 - t16) * s25 + 28.0_dp * t8) * s25
      t179 = 33.0_dp * t10
      t14 = ((t7 * (t62 + t10) + s34 * (s34 * t83 - t14) + t30 * (t135 -
     # t40)) * s56 + s16 * ((t131 - t146) * s25 + t10 * t68) + t34 * s16
     # * (t64 + t94) - 9.0_dp * t60 - ((-(13.0_dp * s16 + t43 + t159) * 
     #s34 - (-s26 * t32 - t26 + t42) * s25 - 15.0_dp * s16 * (s16 - s26)
     #) * s34 + t145 + t178 + t179) * s34) * s56 + s16 * (s16 * ((-t110 
     #+ t52) * s25 + t53) + t165 * t65 + t125 * t105) + 9.0_dp * t61 + (
     #((-(t108 + t137 + t172) * s34 + (t167 + t168 + t163) * s26 - t157 
     #- t57) * s34 + ((-t97 + t151) * s25 + t152) * s25 + (t15 * t175 - 
     #t176) * s26 + 114.0_dp * t10) * s34 + ((t131 - 82.0_dp * t8) * s25
     # - t119) * s25 - t15 * ((t178 + t179) * s26 + t154) - 93.0_dp * t9
     # - t14 * t65 * t44) * s34
      t14 = t14 * s56 + t10 * ((-t49 * t8 + 12.0_dp * t18) * s25 + t113)
     # + t55 * t177 * (s16 * ((t15 * t18 - t109) * s25 + t174) + t76 * t
     #105) + t65 * s16 * (t161 - t67) + t7 * ((-s16 + t16 + t50 - s34) *
     # s34 + (-t108 - t137 - t162) * s26 + t31 + 22.0_dp * t8) * t60 - (
     #(-(((9.0_dp * s16 + t150) * s25 - 81.0_dp * t8) * s25 + ((0.3e1_dp
     # / 0.2e1_dp * t169 + t163) * s26 - t170) * s26 - 174.0_dp * t10) *
     # s34 - (((-t29 - t166) * s25 + 99.0_dp * t8) * s25 + 79.0_dp * t10
     #) * s25 - ((t175 * t55 - t176) * s26 + ((69.0_dp * s16 + t147) * s
     #25 + 225.0_dp * t8) * s25 + 342.0_dp * t10) * s26 - 201.0_dp * t9)
     # * s34 + t79 * ((((s16 * t171 - t107) * s25 + 70.0_dp * t8) * s25 
     #- t174) * s25 + (t15 * ((t144 + t37) * s25 + t111) + t116 * t65) *
     # s26 + 111.0_dp * t9)) * s34
      t13 = t14 * s56 + t15 * t76 * s16 * (t8 * t36 * t47 + t105 * t44) 
     #+ 11.0_dp * t11 - s34 * ((((-t59 * (-t56 + t40) - ((t126 * t7 - t1
     #72) * s26 + t63 * t25 * t55) * s26 - t138) * s34 - s16 * (((-40.0_
     #dp * s16 + 13.0_dp * s25) * s25 - 102.0_dp * t8) * s25 - 160.0_dp 
     #* t10) - (((t169 + t159) * s26 - t170) * s26 - t55 * ((-t78 + t28)
     # * s25 + 58.0_dp * t10)) * s26) * s34 - t79 * ((((-t107 + t166) * 
     #s25 + 61.0_dp * t8) * s25 + 66.0_dp * t10) * s25 + 135.0_dp * t9 +
     # s26 * (((133.0_dp * s16 - 68.0_dp * s25) * s25 + t152) * s26 + ((
     #137.0_dp * s16 + t147) * s25 + 92.0_dp * t8) * s25 + 267.0_dp * t1
     #0) - t34 * t156 * t46)) * s34 + t130 * (((t144 + t106) * s25 + t10
     #) * s25 + (t2 * (((s16 * t143 - t140) * s25 - t102) * s25 + 51.0_d
     #p * t10) + t145) * s26 + 60.0_dp * t9)) - t139 * t13 + 8.0_dp * t9
     # * t20 - t90 * t12 + t55 * t54 * t44 / 2.0_dp - t173 + t40 * t58
      t14 = t160 * s16 * s26 * t45
      t11 = -s25 * t11 + t13 * s56 + t46 * t54 + t60 * (t10 * (30.0_dp *
     # s16 + 35.0_dp * s25) + s26 * (t126 * t44 + t138) - 9.0_dp * t45 +
     # 9.0_dp * t164 + t155 * t63 * t25) + t130 * (t18 * ((-t29 + t166) 
     #* s25 + t121) + ((t83 * s26 + (t158 - t26) * s25 + t146) * s26 + (
     #(-t97 + t51) * s25 + t55 * t8) * s25 + 75.0_dp * t10) * s26 + 30.0
     #_dp * t10 * t5) * t59 + (-s26 + t3) * t101 * t59 * t60 - t79 * ((t
     #15 * (((-t29 + t118) * s25 + t81) * s25 + 30.0_dp * t10) + t44 * (
     #-t2 * (t40 + t80) - t100)) * s26 + t8 * ((50.0_dp * s16 + t66) * s
     #25 + 40.0_dp * t8) + t44 * ((s16 * t136 - t104) * s25 + 54.0_dp * 
     #t8) - t154) * t62 - t79 * t130 * ((((s16 * t112 - t135) * s25 - t1
     #06) * s25 + t119) * s25 + 12.0_dp * t9 + s26 * (((-t16 + t91) * s2
     #5 - t102) * s25 + t113 + t116)) * s34 + t76 * t9 * ((t3 - t140) * 
     #s25 + t71) + t14 + s16 * t1 * t105 * t45
      t13 = t1 * s25
      t14 = t30 * t34
      t25 = t55 * t21
      t28 = t2 * t5
      t31 = (-t8 + t13) * s26
      t32 = s16 + s25 + s26 - s34 + s56
      t42 = 0.1e1_dp / t73
      t43 = 0.1e1_dp / s16
      t49 = 0.1e1_dp / t132 ** 2.0_dp
      t27 = 0.1e1_dp / t27 ** 2.0_dp
      t1 = t17 * MB1001(0) + (-t75 * t41 * t7 + t131 * t10) * s16 - s34 
     #* ((((-t101 * s34 + t34 * s16 * t5 + (t90 + t98 + t99) * s26) * s3
     #4 - t33 * t95 - t8 * (t96 + t97) - s26 * ((33.0_dp * s16 + 21.0_dp
     # * s25) * s26 + t92 + t93)) * s34 + t79 * (((t29 + t87) * s25 + t8
     #8) * s25 + t48 + s26 * ((t85 + t26) * s26 + (t85 + t86) * s25 + 32
     #.0_dp * t8) + t89)) * s34 - t130 * (t34 * (t10 + t22) + s26 * (t84
     # * s26 + (-t3 + t80) * s25 + t81) - t3 * t18)) - t38 * s56 - t39 -
     # t1 * t63 * t45 - t1 * (-s25 * t4 + t37) * t46 - t103 * (-t29 + t4
     #0) * t5
      t1 = t42 * t27 * (t19 * I300s16(1) + (-t32 * s25 * (-t2 * (t62 * (
     #t24 + t29 + t82 - s34) - t124) + (((s25 - t77) * s25 + t121) * s26
     # + (t18 - t134) * s25 + t141) * s26 + (((t107 + t50 + t96) * s25 +
     # t148 + 12.0_dp * t149) * s34 - (t133 + t106) * s16 - t2 * t22 - s
     #26 * ((t29 + t82) * s26 + (t82 + t107) * s25 + t142)) * s34 + ((-t
     #7 * s25 * (s16 + s34) - t127 + t128 + t18) * s56 - t2 * ((t7 * s25
     # * t122 + (t135 + t35 + t118) * s16) * s34 - t76 * t63) + t70 + t5
     #9 * (t35 + t90 + t108) + t20 - t115 * t34) * s56 - t123) * MB1110(
     #1) + t23 + (t2 * (t75 * t36 * t47 + t129 * t8) * t10 - t7 * (t12 *
     # t18 + t61 * (t44 * (s16 - t16) + t7 * ((s16 * t117 - t44) * s26 +
     # t115) + t141)) + t11) * MB1011(1) * t49) * t43 + t49 * t1)
      result = t1 / 4.0_dp + t43 * (t27 * (t32 * (-t2 * (-s26 * t8 + ((t
     #122 + t3 - s34) * s34 - t149) * s56) - (-t44 - t8) * s16 - ((-t122
     # * t2 + s34 - t77) * s34 + t2 * s25 * t79 + (t77 + t35) * s16 + t1
     #8 + t44) * s34 + t74 * (s16 - s34)) * MB1111D2(0) + s25 * s34 * ((
     #t5 * (t2 * t21 - t30 * t7) - t31) * s26 + ((t55 * (t18 + t8 + t120
     # + t177) + t44 + 8.0_dp * t30 + t59) * s34 - (t28 * s26 + t14 + t2
     #5) * s26 - t15 * ((t8 + t6) * s25 + t10)) * s34 + ((-t2 * t72 - t1
     #3 + t59 + t8) * s56 + t2 * (t59 * (s26 + t153 - s34) + t10 + t20 -
     # t31) - (t15 * t69 + t14 + t25) * s34 - t30 * t5) * s56 + t114 * t
     #18 + t9 - t2 * (s26 + t28) * t62) * MB1101(0) * t42) + s34 * (s16 
     #+ s26 - s34 + s56) * MB1011(0) * t49) / 2.0_dp

           intHs16s25s26s34s56x1123D6eps0 = result
       end function intHs16s25s26s34s56x1123D6eps0

       function intHs16s25s26s34s56x1123D6eps1()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1123D6eps1
           complex(dp) ::  t1,t2

           complex(dp) :: result

      t1 = s16 + s26 + s56
      t2 = 0.1e1_dp / t1 ** 2.0_dp
      result = 0.1e1_dp / s16 * t2 * (-s34 / 2.0_dp + t1 / 4.0_dp + s34 
     #* (s16 + s26 - s34 + s56) * MB1011(1) / 2.0_dp)

           intHs16s25s26s34s56x1123D6eps1 = result
       end function intHs16s25s26s34s56x1123D6eps1

       function intHs16s25s26s34s56x1130D4eps0()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1130D4eps0
           complex(dp) ::  t1,t2

           complex(dp) :: result

      t1 = s16 + s25 + s26 - s34 + s56
      t2 = 3.0_dp
      t1 = 0.1e1_dp / t1 ** 2.0_dp
      result = -(s25 * t2 + s16 + s26 - s34 + s56) * t1 * ((s16 + s26 + 
     #s56) * MB1011(1) - MB1001(0)) / 4.0_dp + s25 * t1 * (s25 * MB1110(
     #0) + 1.0_dp) / 2.0_dp + 0.3e1_dp / 0.2e1_dp * s25 ** 2.0_dp * MB11
     #10(1) * t1

           intHs16s25s26s34s56x1130D4eps0 = result
       end function intHs16s25s26s34s56x1130D4eps0

       function intHs16s25s26s34s56x1130D4eps1()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1130D4eps1
           complex(dp) ::  t1,t2,t3

           complex(dp) :: result

      t1 = 3.0_dp
      t2 = s16 + s25 + s26 - s34 + s56
      t3 = 0.1e1_dp / 0.4e1_dp
      t2 = 0.1e1_dp / t2 ** 2.0_dp
      result = t2 * (t3 * (s25 * t1 + s16 + s26 - s34 + s56) + s25 ** 2.
     #0_dp * MB1110(1) / 2.0_dp)

           intHs16s25s26s34s56x1130D4eps1 = result
       end function intHs16s25s26s34s56x1130D4eps1

       function intHs16s25s26s34s56x1131D4eps0()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1131D4eps0
           complex(dp) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           complex(dp) ::  t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30
           complex(dp) ::  t4,t5,t6,t7,t8,t9

           complex(dp) :: result

      t1 = s16 + s25 - s34
      t2 = 2.0_dp
      t3 = -t2 * (s26 + s56) - t1
      t4 = s16 + s25
      t5 = t2 * s26
      t6 = s16 * s25 + (-s34 + t4 + s26) * s26 + (t5 + t1 + s56) * s56
      t7 = s16 + s26 - s34 + s56
      t8 = s16 + s26
      t9 = s16 ** 2.0_dp
      t10 = s26 ** 2.0_dp
      t11 = s56 ** 2.0_dp
      t12 = s16 + s26 + s56
      t13 = s16 + s25 + s26 - s34 + s56
      t14 = s16 * t8
      t15 = -s16 + s26 + s25
      t16 = 3.0_dp
      t17 = t16 * t14
      t18 = s16 - s25
      t14 = 6.0_dp * t14
      t19 = s34 ** 2.0_dp
      t20 = t16 * s16
      t21 = s25 + s26
      t22 = t2 * s16
      t23 = t22 + t21
      t24 = s25 ** 2.0_dp
      t25 = t9 + t24
      t26 = t18 ** 2.0_dp
      t27 = t26 * s26
      t28 = t9 * s25
      t29 = 0.1e1_dp / s16
      t7 = 0.1e1_dp / t7
      t13 = 0.1e1_dp / t13
      t30 = 0.1e1_dp / t6
      t6 = 0.1e1_dp / t6
      t12 = 0.1e1_dp / t12
      t21 = -(t2 * (s34 * t19 * t23 + (s34 * t8 * (t22 - s25 + s26) + t2
     #8 - t27) * s16 - (-t28 * t2 + s16 * t25 + (-s16 * (-s25 + t20 + t5
     #) + (-s34 + t20 + t21) * s34) * s34 + t27) * s56) - t10 * t26 - t1
     #1 * (-t2 * s16 * (s25 + s34) + t19 + t24 + t9) - t19 * (t2 * s25 *
     # t8 + t10 + t14 + t19 + t24) - t25 * t9) * MB1111D2(0) * t29 * t7 
     #* t30 + I300s16(0)
      t6 = -t3 * I300s16(1) * t6
      t1 = t13 * (-t24 * MB1110(0) * t29 * t7 + t30 * (-(t2 * (s16 * s26
     # + (s16 + s26 - s34) * s56) + (-t2 * t8 - s25 + s34) * s34 + t10 +
     # t11 + t9) * MB1001(0) + (t16 * t9 * s26 * t8 + s16 * (s16 * t9 + 
     #s26 * t10) - ((-t2 * s26 * t18 + (s16 - s26) * s34 - t16 * t9 + t1
     #0) * s34 + t8 * (-t15 * s25 + t17)) * s34 - ((-s16 * s56 - t17 + s
     #34 * (s34 - s25 + t20)) * s56 - t16 * s16 * (t10 + t9) - s34 * (t1
     #9 + (s25 + t5) * s25 - t14) - 6.0_dp * t9 * s26 + t2 * t15 * t19) 
     #* s56) * MB1011(1) * t29) * t12 - s25 * (s26 * t4 + (-t2 * t4 - s2
     #6 + s34) * s34 + s56 * t1 + t24 + t9) * MB1110(1) * t29 * t30) + t
     #7 * t21 + s25 * s34 * (s16 * t18 + s26 * t18 + (-t23 + s34) * s34 
     #+ (s16 - s25 - s34) * s56) * MB1101(0) * t29 * t7 ** 2.0_dp * t30 
     #+ t6 - t19 * MB1011(0) * t29 * t12 * t7
      result = -t1 / 2.0_dp + t3 * t29 * t30

           intHs16s25s26s34s56x1131D4eps0 = result
       end function intHs16s25s26s34s56x1131D4eps0

       function intHs16s25s26s34s56x1131D4eps1()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1131D4eps1
           complex(dp) ::  t1,t2,t3,t4,t5,t6,t7

           complex(dp) :: result

      t1 = s16 + s26 + s56
      t2 = 2.0_dp
      t3 = s16 + s25 + s26 - s34 + s56
      t4 = s16 + s26 - s34 + s56
      t3 = 0.1e1_dp / t3
      t5 = 0.1e1_dp / t1
      t6 = 0.1e1_dp / s16
      t4 = 0.1e1_dp / t4
      t7 = -0.1e1_dp / 0.2e1_dp
      result = t7 * (t4 * (-(s34 ** 2.0_dp * MB1011(1) * t5 + s25 ** 2.0
     #_dp * MB1110(1) * t3) * t6 + I300s16(1)) + (t1 * t2 + s25 - s34) *
     # t6 * t5 * t3)

           intHs16s25s26s34s56x1131D4eps1 = result
       end function intHs16s25s26s34s56x1131D4eps1

       function intHs16s25s26s34s56x1132D6eps0()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1132D6eps0
           complex(dp) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           complex(dp) ::  t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30
           complex(dp) ::  t31,t32,t33,t34,t35,t36,t37,t38,t39,t4,t40,t41
           complex(dp) ::  t42,t43,t44,t45,t46,t47,t48,t49,t5,t50,t51,t52
           complex(dp) ::  t53,t54,t55,t56,t57,t58,t59,t6,t60,t61,t62,t63
           complex(dp) ::  t64,t65,t66,t67,t68,t69,t7,t70,t71,t72,t73,t74
           complex(dp) ::  t8,t9

           complex(dp) :: result

      t1 = s16 + s26 - s34 + s56
      t2 = s16 + s26
      t3 = 4.0_dp
      t4 = 2.0_dp
      t5 = t4 * s16
      t6 = 3.0_dp * s25
      t7 = t3 * s26
      t8 = 3.0_dp * s26
      t9 = s16 + t8
      t10 = 8.0_dp
      t11 = 6.0_dp * s16
      t12 = t10 * s26
      t13 = s16 ** 2.0_dp
      t14 = t13 ** 2.0_dp
      t15 = s16 * t13
      t16 = s26 ** 2.0_dp
      t17 = s26 * t16
      t18 = s34 ** 2.0_dp
      t19 = s34 * t18
      t20 = s16 * s26
      t21 = s16 + s26 + s56
      t22 = s16 + s25
      t23 = t4 * s26
      t24 = s16 * s25
      t25 = (-s34 + t22 + s26) * s26 + (t23 + s16 - s34 + s25 + s56) * s
     #56 + t24
      t26 = 7.0_dp
      t27 = t4 * s25
      t28 = t10 * s16
      t29 = t26 * t13
      t30 = 5.0_dp
      t31 = t30 * s25
      t32 = s26 + s25
      t33 = s25 * t32
      t34 = 14.0_dp
      t35 = 15.0_dp * s16
      t36 = s25 ** 2.0_dp
      t37 = (s16 + t31) * s26
      t38 = s16 - s34
      t39 = 9.0_dp * s26
      t40 = -t19 + t15
      t41 = s56 ** 2.0_dp
      t42 = t41 ** 2.0_dp
      t43 = s56 * t41
      t44 = t13 * s25
      t45 = s16 - s25
      t46 = s16 * t45
      t47 = s16 * t2
      t48 = 6.0_dp * t47
      t49 = 3.0_dp * s16
      t50 = t3 * s25
      t51 = s34 * t4
      t52 = t30 * s16
      t53 = t3 * s16
      t54 = t30 * s26
      t55 = t4 * s25 * t2
      t56 = t28 - t6
      t57 = s16 * t26 - t27
      t58 = t36 + t16
      t59 = 6.0_dp * s26
      t60 = t27 * s26
      t61 = s26 * t32
      t62 = -s25 + t11
      t63 = 32.0_dp * t2
      t64 = t6 * s26
      t65 = 24.0_dp
      t66 = t2 ** 2.0_dp
      t67 = s16 * t66
      t68 = s25 * t22
      t69 = t62 * s26
      t70 = 6.0_dp * s25
      t71 = 3.0_dp * t18
      t72 = t17 * (-t20 + s16 * (s25 - t11)) + t14 * (-9.0_dp * s16 + t5
     #0)
      t73 = s16 * t2 * t66 * (t10 * t2 + t6) * s34
      t74 = t2 * ((16.0_dp * s26 + t28 + t31) * s16 - 6.0_dp * t61) * t1
     #9
      t42 = s16 * t14 * s25 - s16 * s56 * t42 + t72 * s26 + ((-t16 ** 2.
     #0_dp * t30 - 9.0_dp * t14) * s16 + t3 * s16 * (s25 * t15 + (-t13 *
     # t56 - t16 * t62) * s26) + (t18 * ((t59 - t52) * s34 + s16 * (s16 
     #* t65 - s25) + (20.0_dp * s16 - 12.0_dp * s25 - 18.0_dp * s26) * s
     #26) + t67 * (9.0_dp * s25 + t63)) * s34 - 6.0_dp * t13 * t57 * t16
     # - 6.0_dp * t2 * (-t4 * (-t24 + t16) + t29 - t36 - t64 + t20 * t10
     #) * t18) * s56 + t41 * ((-16.0_dp * t15 - 10.0_dp * t17) * s16 + t
     #18 * (t51 * (-t39 + t52 - t6) + (-45.0_dp * s16 - 36.0_dp * s26) *
     # s16 + 18.0_dp * t61 + 3.0_dp * t68) + 3.0_dp * s34 * (t47 * (16.0
     #_dp * t2 + t6) + t19) + 6.0_dp * s16 * ((-s16 * t57 - t69) * s26 +
     # t44)) + t42 * ((s34 * t10 + s25 - t11 - t54) * s16 + t71) + t43 *
     # ((-t13 * t34 - 10.0_dp * t16) * s16 + t3 * s16 * (t24 - t69) - 6.
     #0_dp * t19 + s34 * (s16 * (t6 + t63) + s34 * (-12.0_dp * s16 + 12.
     #0_dp * s26 + t70))) + t73 - t9 * (-s26 + t5) * t18 ** 2.0_dp + t74
      t61 = t32 * t4 + t49
      t62 = t5 + t32
      t63 = t13 * s26
      t65 = s16 * t36
      t69 = 10.0_dp * s16
      t72 = t63 * t45
      t73 = 0.1e1_dp / s16
      t21 = 0.1e1_dp / t21 ** 2.0_dp
      t25 = 0.1e1_dp / t25 ** 2.0_dp
      t14 = (-t4 * t13 * (t16 * (s16 * t56 + s26 * t57) + t14) + t42 - t
     #71 * t66 * ((t59 + t53 + t6) * s16 - t58 - t60)) * MB1011(1) * t21
     # - t30 * t40 + 3.0_dp * t17 + 3.0_dp * t43 + 3.0_dp * t44 + (t37 +
     # (t28 + t27) * s25 - t29) * s26 + (-(s26 * t26 + 10.0_dp * s25 + t
     #35) * s34 + (s25 * t26 + s26 * t34 + t35) * s16 + t30 * t33 - t16)
     # * s34 + ((t39 + t31 + t38) * s56 + (s16 * t34 - t23 + t31) * s34 
     #+ 9.0_dp * t16 - t26 * (t18 + t13) + t4 * (t37 + t36) + t28 * s25)
     # * s56 + t5 * t36 + s25 * (t4 * t40 - (-(-s25 + t5) * s26 - t13 * 
     #t3 + t68) * s26 - (-t4 * t61 * s34 + (t12 + t11 + t50) * s16 + t4 
     #* t58 + t64) * s34 - (-(t38 * t4 - s25) * s56 - (-t28 - t6) * s34 
     #- t3 * (-s26 * s34 + t13 + t18 + t20) + t60 + t68) * s56 - t65) * 
     #MB1110(1)
      t26 = 0.1e1_dp / t1
      t10 = (((-t15 * t45 - t17 * t45) * s16 + ((((s34 - 3.0_dp * t32 - 
     #t52) * s34 + (s25 * t10 + 12.0_dp * s26 + t69) * s16 + 3.0_dp * t5
     #8 + t70 * s26) * s34 - s25 * t36 - t13 * (t70 + t69) - t17 - 3.0_d
     #p * ((s25 + t49) * s26 + (s25 + t49) * s25 + 6.0_dp * t13) * s26 -
     # 3.0_dp * t65) * s34 + t67 * (t52 + t23)) * s34 + (-t41 * (-t5 * s
     #34 + t18 + t46) + 6.0_dp * (t18 * t62 + t47 * (s26 + t5)) * s34 - 
     #6.0_dp * t72 - 3.0_dp * (t13 * t45 + t16 * t45) * s16 - 3.0_dp * t
     #18 * (t18 + (t49 + t23) * s25 + t48 + t58)) * s56 - 3.0_dp * t41 *
     # ((s26 * t45 - s34 * (t49 + t23) + t46) * s16 + t18 * (-s34 + t49 
     #+ t32)) - 3.0_dp * t72 * t2) * MB1111D2(0) - s25 * s34 * (t4 * (((
     #-t62 + s34) * s34 + t47) * s56 + t63) + (t16 + t13) * s16 + ((t61 
     #- s34) * s34 - (t7 + t49) * s16 - t55 - t58) * s34 + t38 * t41) * 
     #MB1101(0)) * t25 * t26
      t10 = t73 * (t10 + t18 * MB1011(0) * t21)
      result = t25 * ((-t4 * (-s25 * t16 + t15 - t19) + (t16 - 3.0_dp * 
     #t46 + t36) * s26 + (t49 * s25 - (t11 + t50 + t8) * s34 + t33 * t4 
     #+ t48) * s34 + ((t8 + t27 + s56) * s56 + (s25 + t51 + t7) * s25 + 
     #3.0_dp * t16 - 3.0_dp * t18 - 3.0_dp * t46 + t11 * s34) * s56 + t2
     #4 * t22) * I300s16(1) + t14 * t73 + t21 * ((t4 * (t18 * (s16 + t23
     #) + t15) + (s16 * (t52 - t27) + (s26 - s25 + t53) * s26) * s26 + (
     #(-s34 * t30 - s25 + s56 + t53 + t8) * s56 + s16 * (-9.0_dp * s34 +
     # t52) + (-10.0_dp * s34 + t28 + t8) * s26 + t3 * s34 * (-s25 + s34
     #) - t55) * s56 - t44 - t2 * (t22 * t3 + t54) * s34) * MB1001(0) + 
     #(t9 * s34 - t2 * (t7 + t6 + t5)) * s34 + ((-s34 * t3 + s56 + 3.0_d
     #p * t2) * s56 - (t12 + t11 + t6) * s34 + 3.0_dp * t18 + 3.0_dp * t
     #16 + 3.0_dp * t13 + t11 * s26) * s56 + t15 + t17 + 3.0_dp * t20 * 
     #t2) * t1) / 4.0_dp + t10 / 2.0_dp

           intHs16s25s26s34s56x1132D6eps0 = result
       end function intHs16s25s26s34s56x1132D6eps0

       function intHs16s25s26s34s56x1132D6eps1()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1132D6eps1
           complex(dp) ::  t1,t2

           complex(dp) :: result

      t1 = s16 + s26 + s56
      t2 = 0.1e1_dp / t1 ** 2.0_dp
      result = 0.1e1_dp / s16 * t2 * (s34 / 2.0_dp + t1 / 4.0_dp + s34 *
     #* 2.0_dp * MB1011(1) / 2.0_dp)

           intHs16s25s26s34s56x1132D6eps1 = result
       end function intHs16s25s26s34s56x1132D6eps1

       function intHs16s25s26s34s56x1141D6eps0()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1141D6eps0
           complex(dp) ::  t1,t10,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109
           complex(dp) ::  t11,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t12
           complex(dp) ::  t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t13,t130
           complex(dp) ::  t131,t132,t133,t134,t135,t136,t137,t138,t139,t14,t140,t141
           complex(dp) ::  t142,t143,t144,t145,t146,t147,t148,t149,t15,t150,t151,t152
           complex(dp) ::  t153,t154,t155,t156,t157,t158,t159,t16,t160,t161,t162,t163
           complex(dp) ::  t164,t165,t166,t167,t168,t169,t17,t170,t171,t172,t173,t18
           complex(dp) ::  t19,t2,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29
           complex(dp) ::  t3,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t4
           complex(dp) ::  t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t5,t50
           complex(dp) ::  t51,t52,t53,t54,t55,t56,t57,t58,t59,t6,t60,t61
           complex(dp) ::  t62,t63,t64,t65,t66,t67,t68,t69,t7,t70,t71,t72
           complex(dp) ::  t73,t74,t75,t76,t77,t78,t79,t8,t80,t81,t82,t83
           complex(dp) ::  t84,t85,t86,t87,t88,t89,t9,t90,t91,t92,t93,t94
           complex(dp) ::  t95,t96,t97,t98,t99

           complex(dp) :: result

      t1 = 3.0_dp
      t2 = s16 ** 2.0_dp
      t3 = t2 ** 2.0_dp
      t4 = s16 * t2
      t5 = t2 * t3
      t6 = s16 * t3
      t7 = 19.0_dp * s16
      t8 = t1 * s25
      t9 = 25.0_dp
      t10 = t9 * s25
      t11 = 17.0_dp * s16
      t12 = s16 - s25
      t13 = s26 * t12
      t14 = 5.0_dp
      t15 = 4.0_dp
      t16 = s26 ** 2.0_dp
      t17 = t16 ** 2.0_dp
      t18 = s26 * t17
      t19 = s26 * t16
      t20 = t14 * s25
      t21 = t1 * s16
      t22 = 15.0_dp * t2
      t23 = 2.0_dp
      t24 = 15.0_dp * s16
      t25 = 10.0_dp * s25
      t26 = t23 * s26
      t27 = s16 - s34
      t28 = s26 * t27
      t29 = s34 ** 2.0_dp
      t30 = t29 ** 2.0_dp
      t31 = s34 * t30
      t32 = s34 * t29
      t33 = t29 + t2
      t34 = t23 * s34 * t12
      t35 = 42.0_dp
      t36 = -t32 + t4
      t37 = s56 ** 2.0_dp
      t38 = t37 ** 2.0_dp
      t39 = s56 * t37
      t40 = s56 * t38
      t41 = t39 + t19
      t42 = 13.0_dp * s16
      t43 = 6.0_dp * s25
      t44 = s16 * s25
      t45 = s16 + s25
      t46 = (-s34 + t45 + s26) * s26 + (s16 + s25 - s34 + t26 + s56) * s
     #56 + t44
      t47 = 12.0_dp * s25
      t48 = t14 * t2
      t49 = 9.0_dp
      t50 = t49 * s25
      t51 = 6.0_dp * s16
      t52 = s16 * t49
      t53 = (-t43 - t52) * s25
      t54 = 12.0_dp * s16
      t55 = t8 + t54
      t56 = s25 * t55 - t48
      t57 = t8 - s16
      t58 = s16 + s26
      t59 = t23 * s25
      t60 = 13.0_dp * t2
      t61 = 11.0_dp * s16
      t62 = s25 ** 2.0_dp
      t63 = s25 * t62
      t64 = t1 * t19
      t65 = 35.0_dp
      t66 = 14.0_dp * s25
      t67 = 34.0_dp * s16
      t68 = 20.0_dp * s25
      t69 = 40.0_dp * s16
      t70 = (t50 + t69) * s25
      t71 = 38.0_dp * t2
      t72 = 46.0_dp * s16
      t73 = 7.0_dp
      t74 = t73 * s25
      t75 = 10.0_dp * t2
      t76 = 12.0_dp * s26
      t77 = t1 * s26
      t78 = t77 + s16
      t79 = t15 * t56
      t80 = -t57
      t81 = 20.0_dp * s16
      t82 = 53.0_dp * s16
      t83 = t15 * s25
      t84 = t15 * s16
      t85 = t59 + t84
      t86 = t85 * s25
      t87 = 15.0_dp * t19
      t88 = (55.0_dp * s16 + t66) * s25
      t89 = 66.0_dp * t2
      t90 = 138.0_dp * s16
      t91 = 60.0_dp * s25
      t92 = t9 * s16
      t93 = 50.0_dp * t2
      t94 = t35 * t4
      t95 = 54.0_dp * s26
      t96 = 24.0_dp
      t97 = t96 * s26
      t98 = t3 * s25
      t99 = s26 * (t48 + t53)
      t100 = 6.0_dp * t99
      t101 = s16 * ((-t51 - t50) * s25 + t48)
      t102 = 12.0_dp * t62
      t103 = (t84 + s25) * s25
      t104 = 19.0_dp * t2
      t105 = t14 * s16
      t106 = 18.0_dp * t2
      t107 = t4 * s25
      t108 = s16 * t62
      t109 = 15.0_dp * s26
      t110 = 23.0_dp * s16
      t111 = t80 * t16
      t112 = t2 * s25
      t113 = t1 * t62
      t114 = s25 + s34
      t115 = t58 ** 2.0_dp
      t104 = (((t1 * t114 - s16) * s56 + (t109 + t11 + s25) * s34 - t14 
     #* (s26 * t80 + t2) - 12.0_dp * t29 + 12.0_dp * t44 + t113) * s56 +
     # s26 * t79 + t23 * s34 * (-(t25 + t110 + t97) * s34 + (t59 + t109 
     #+ t67) * s26 + t103 + t104) - 10.0_dp * t111 - 10.0_dp * t4 + 18.0
     #_dp * t32 + 18.0_dp * t112 + t102 * s16) * s56 + (18.0_dp * t108 -
     # t100) * s16 + (((t10 + t95 + t72) * s34 - (72.0_dp * s26 + t90 + 
     #t91) * s26 - t88 - t89) * s34 + ((t105 + s25) * s25 + t106) * s25 
     #+ 30.0_dp * t19 + t94 + 6.0_dp * s26 * ((t11 + s25) * s26 + t103 +
     # t104)) * s34 + 6.0_dp * t16 * t56 - 10.0_dp * t19 * t80 - 10.0_dp
     # * t3 - 12.0_dp * t30 + 12.0_dp * t107
      t79 = t104 * s56 + s16 * ((-t101 * t15 - t100) * s26 + t102 * t2) 
     #+ t1 * (t31 + t98) - t14 * (t17 * t80 + t6) + (((-(t50 + t11 + t97
     #) * s34 + (92.0_dp * s16 + 50.0_dp * s25 + t95) * s26 + t70 + t71)
     # * s34 - ((t8 + t92) * s25 + t93) * s25 - ((48.0_dp * s26 + t90 + 
     #t91) * s26 + t23 * (t89 + t88)) * s26 - t94) * s34 + t58 * ((16.0_
     #dp * t2 + t86) * s25 + 23.0_dp * t4 + s26 * ((t83 + t82) * s26 + (
     #t43 + t81) * s25 + 61.0_dp * t2) + t87)) * s34 + t19 * t79
      t88 = s16 * t16
      t89 = s16 + s26 + s56
      t90 = s16 + s25 + s26 - s34 + s56
      t91 = s16 + s26 - s34 + s56
      t94 = 16.0_dp * s16
      t95 = 11.0_dp * s25
      t97 = -22.0_dp
      t100 = 18.0_dp * s26
      t104 = 10.0_dp * s16
      t116 = -s25 * t23 * t45 + t2
      t86 = t2 - t86
      t117 = t105 + s25
      t118 = t59 + s16
      t119 = t23 * s16
      t120 = t59 + t105
      t121 = -t74 + t104
      t122 = 13.0_dp * s25
      t123 = 10.0_dp * s26
      t124 = 20.0_dp * t4
      t125 = 60.0_dp * t2
      t126 = 29.0_dp * s25
      t127 = 60.0_dp * s16
      t128 = 40.0_dp * t2
      t129 = 20.0_dp * t2
      t82 = (t82 + t122) * s25
      t130 = 90.0_dp * t2
      t131 = 8.0_dp * s25
      t132 = 70.0_dp * t4
      t133 = 26.0_dp * s16
      t134 = 30.0_dp * s26
      t135 = s16 * s26
      t136 = 6.0_dp * s26
      t137 = s25 * t117
      t138 = s16 * t58
      t139 = t73 * s26
      t140 = 40.0_dp * s25
      t141 = s25 + s26
      t103 = (((-s34 * t73 + t134 + t20 + t92) * s25 + 15.0_dp * s26 * t
     #58 + 10.0_dp * s34 * t27 + t37) * s56 + (t15 * ((-t84 - t139 - s25
     #) * s25 + 10.0_dp * t138) - (20.0_dp * s26 + t92 + t131) * t23 * s
     #34) * s34 + t2 * (t140 - t104) + 20.0_dp * (t16 + t137) * s26 + 20
     #.0_dp * t108 + 20.0_dp * t32 + 30.0_dp * t118 * t16) * s56 + (((60
     #.0_dp * t58 + t126) * s34 - (150.0_dp * s16 + 48.0_dp * s25 + 60.0
     #_dp * s26) * s26 - t130 - t82) * s34 + s16 * ((-t51 - t95) * s25 +
     # t125) + 6.0_dp * t121 * t16 - t63 + t76 * (t75 - t103)) * s34 + 3
     #0.0_dp * (-s16 * t86 + (s26 * t118 + t137) * s26) * s26 + 30.0_dp 
     #* t112 * t45 - 15.0_dp * t30 - 15.0_dp * t3 + 15.0_dp * t17
      t52 = t103 * s56 + t15 * t31 + t17 * (15.0_dp * t118 + t136) + t23
     # * s34 * (t29 * ((29.0_dp * s26 + t43 + t133) * s25 + 30.0_dp * t1
     #6 + 30.0_dp * t2 + t127 * s26) + t58 * (s16 * ((t84 - t20) * s25 +
     # t129) + (t23 * t121 * s26 + (-t43 - t104) * s25 + t128) * s26 - t
     #63)) + t3 * (t25 - t52) - t29 * (t29 * (t47 + t92 + t134) + ((t83 
     #+ t92) * s25 + 58.0_dp * t2) * s25 + ((40.0_dp * s26 + 6.0_dp * t1
     #31 + 6.0_dp * t92) * s26 + t23 * (t130 + t82)) * s26 + t132) + 20.
     #0_dp * s25 * (t117 * t19 + t107) - 30.0_dp * t135 * (s16 * t116 + 
     #s26 * t86)
      t52 = t1 * (t118 * t18 + t40 * (t141 * t23 + s16)) + t14 * s25 * (
     #t117 * t17 + t98) + t23 * (t31 * (t26 + s16) - t5) + s25 * t6 + (t
     #4 * ((t104 + t68) * s25 - t2 * t49) + t18) * s26 + s34 * (((-(t119
     # + t77) * (t14 * t58 + t83) * s34 + s26 * ((t126 + t127) * s26 + (
     #52.0_dp * s16 + t47) * s25 + t125) + 20.0_dp * t4 + 20.0_dp * t19 
     #+ t44 * (t25 + t110)) * s34 - t58 * (((t83 + t54) * s25 + 21.0_dp 
     #* t2) * s25 + ((8.0_dp * t120 + t123) * s26 + (37.0_dp * s16 + t12
     #2) * s25 + t93) * s26 + t124)) * s34 + t115 * (s16 * ((-t8 + t105)
     # * s25 + t75) + t121 * t16 - t63 + t26 * (t119 - s25) * t120)) + t
     #52 * s56 - t88 * (t116 * t24 + t123 * t86)
      t82 = s16 * t12
      t86 = t119 - t20
      t98 = t20 + t51
      t103 = 14.0_dp * s16
      t107 = (-t25 - t103) * s25 + t48
      t108 = (-t20 - t42) * s25 + t2 * t23
      t110 = 31.0_dp * s16
      t116 = (t20 + t110) * s25
      t117 = t1 * t2 + t116
      t118 = t119 + t8
      t120 = 8.0_dp * s16
      t121 = 6.0_dp * t2
      t125 = t96 * t4
      t126 = t84 - s25
      t127 = (-t59 + t21) * t126
      t130 = t127 * s26
      t137 = s16 * t63
      t142 = t14 * t4
      t143 = 75.0_dp
      t144 = 48.0_dp * t2
      t145 = t21 + s25
      t146 = t24 - s25
      t147 = 30.0_dp * s16
      t148 = 21.0_dp * s25
      t149 = (t105 - t43) * s25
      t150 = 37.0_dp * t2
      t151 = t150 + t149
      t152 = 27.0_dp
      t153 = t152 * s16
      t154 = -t10 + t153
      t155 = t15 * t4
      t156 = 60.0_dp * t4
      t157 = t130 * t14
      t158 = 65.0_dp * s16
      t159 = ((-t21 - t83) * s25 + t144) * s25
      t160 = 79.0_dp * t4
      t161 = (-78.0_dp * s16 + 63.0_dp * s25) * s25 - 453.0_dp * t2
      t162 = -124.0_dp * s16
      t163 = 56.0_dp * s25
      t164 = t163 + t162
      t165 = 36.0_dp * s26
      t166 = 12.0_dp * t118
      t167 = t117 * t14
      t168 = s16 * t108
      t169 = t22 * t107
      t170 = t9 * t3
      t55 = t1 * (s16 * t55 + (-s26 * t49 + t104 - t131) * s26) * t31
      t55 = s16 * (s26 * (s26 * (s26 * (s26 * (s26 * (-t166 - t139) - t1
     #67) + 20.0_dp * t168) + t169) + t155 * t86 * t98) - t170 * t62) + 
     #s34 * ((((t29 * (t105 - t136) + s16 * ((s16 * t143 - t8) * s25 + 1
     #05.0_dp * t2) + ((t1 * t154 - t165) * s26 + 6.0_dp * t151) * s26) 
     #* s34 + s16 * (((-36.0_dp * s16 + t122) * s25 - 154.0_dp * t2) * s
     #25 - 160.0_dp * t4) + (-6.0_dp * t160 - 6.0_dp * t159 + t87) * s26
     # + t16 * (s26 * t164 + t161)) * s34 + t58 * ((((t21 - t43) * s25 +
     # t71) * s25 + 108.0_dp * t4) * s25 + 135.0_dp * t3 + ((t14 * t146 
     #* s26 + (-t83 + t158) * s25 + 285.0_dp * t2) * s26 + ((-t50 + t67)
     # * s25 + 178.0_dp * t2) * s25 + 345.0_dp * t4) * s26)) * s34 - t11
     #5 * ((((t43 + s16) * s25 - 29.0_dp * t2) * s25 + t4 * t49) * s25 +
     # (t23 * (((t131 - t103) * s25 - 23.0_dp * t2) * s25 + t156) + t157
     #) * s26 + 60.0_dp * t3)) + 11.0_dp * t5 * t12 - t55
      t87 = -43.0_dp
      t103 = t127 * t123
      t122 = 120.0_dp * t3
      t123 = 90.0_dp * s16
      t171 = s25 * (-t43 + t123) + 540.0_dp * t2
      t172 = 10.0_dp * t146 * s26
      t173 = 10.0_dp * t19
      t18 = s16 * (t169 * s26 - t173 * t117 - t140 * t3 - 50.0_dp * t4 *
     # t62 + t6 * t96 - 21.0_dp * t18) + t1 * (s34 * (s26 * (-t100 - t10
     # + t153) + t149 + t150 - t29) + t173 - s26 * (-s26 * t164 / 2.0_dp
     # - t161 / 3.0_dp) - t160 - t159) * t32 + (s34 * (t32 * (s26 * t152
     # - t24 + t47) - t1 * t63 * t45 + t2 * (s25 * (143.0_dp * s16 + 36.
     #0_dp * s25) + 240.0_dp * t2) + (s26 * (t171 + t172) + s25 * (s25 *
     # (-t50 + t147) + 243.0_dp * t2) + 630.0_dp * t4) * s26) - t58 * (s
     #25 * (s25 * (t2 * t87 - t53) - 14.0_dp * t4) + s26 * (t15 * (s25 *
     # (s25 * (t43 - t120) - 31.0_dp * t2) + t156) + t103) + t122)) * s3
     #4 + 30.0_dp * t88 * (-t118 * t16 + t168)
      t53 = (-t119 + t83) * s25
      t87 = 39.0_dp * t2
      t109 = s16 * (s16 * (s25 * (-t10 - t158) + t75) + s26 * (s26 * (-s
     #26 * t65 - 30.0_dp * t118) - t167)) - t30 * t49 + (s34 * (s34 * (t
     #109 + t66 - t110) + s16 * (s26 * t143 + t123) + s25 * (-s26 * t14 
     #- s25 + t24)) + s25 * (t87 - t53) - t156 - t157) * s34
      t110 = s16 * (s26 * (21.0_dp * s26 + t166) + t116) + t1 * t36 + s3
     #4 * (-s34 * t146 + t127)
      t116 = s25 * t5 + s26 * (t17 * (t117 + t16) - t3 * (-t62 * t9 + 11
     #.0_dp * t82))
      t123 = t58 * t115 * (((-t2 * t73 + t23 * t62) * s25 + t142) * s25 
     #+ 12.0_dp * t3 + (t130 + ((t83 - t120) * s25 - t121) * s25 + t125)
     # * s26 - t137) * s34
      t54 = t115 * (t1 * t63 * t12 + ((t146 * s26 + (-s25 + t81) * t145)
     # * s26 + ((-t8 + t54) * s25 + t144) * s25 + t143 * t4) * s26 + 30.
     #0_dp * t4 * t45 + t60 * t62) * t29
      t69 = t58 * ((t15 * (t1 * t62 * t12 + t2 * (26.0_dp * s25 + t147))
     # - t64) * s26 + t2 * ((s25 * t96 + 50.0_dp * s16) * s25 + t128) + 
     #t16 * ((-t66 + t67) * s26 + (-t148 + t69) * s25 + 117.0_dp * t2) -
     # t137) * t32
      t18 = s16 * t39 * t38 + s16 * t116 - t55 * s56 - t30 * (t1 * s26 *
     # (s16 * ((t92 - s25) * s25 + t2 * t65) + t151 * s26) + t154 * t19 
     #+ t4 * (s25 * t65 + t147) + t49 * (t2 * t62 - t17)) - t37 * t18 - 
     #t38 * t109 - t39 * (s16 * (-t132 * s25 + 20.0_dp * t135 * t108 - 1
     #0.0_dp * t117 * t16 - 40.0_dp * t118 * t19 - t17 * t65 - t93 * t62
     # + t170) + t31 * t49 + (s34 * (s34 * (s34 * (-t10 - t165 + t153) +
     # s25 * (-t133 + t148) + s26 * (t134 + t163 + t162) - 151.0_dp * t2
     #) + s25 * (s25 * (-t8 + t104) + 81.0_dp * t2) + 210.0_dp * t4 + s2
     #6 * (0.2e1_dp / 0.3e1_dp * t171 + t172)) + s25 * (s25 * (s25 * (-t
     #59 - t61) + t2 * t9) + 46.0_dp * t4) + s26 * (-t15 * (s25 * (-t87 
     #+ t53) + t156) - t103) - t122) * s34) + t40 * t110 + t123 - t54 + 
     #s16 * (t84 + t43 + t139) * t37 * t38 - t78 * (t119 - s26) * t29 * 
     #t30 + t69
      t38 = t16 + t62
      t39 = 6.0_dp * t138
      t40 = t141 + t119
      t49 = t2 + t62
      t53 = s34 - t141 - t21
      t54 = t12 ** 2.0_dp
      t55 = s16 * t49
      t69 = t54 * s26
      t81 = t23 * s16 * t114
      t87 = t40 * t32
      t49 = t2 * t49
      t54 = t54 * t16
      t93 = -t4 + t63
      t103 = s16 * t93
      t93 = t1 * t93 + t44 * (t83 + t105)
      t109 = t119 + s25
      t110 = t1 * t63
      t114 = t14 * t44 * t12
      t106 = s25 * (-t84 - t43) + t106
      t104 = t131 + t104
      t111 = t111 * t109
      t116 = -t111 + s25 * (s25 * (t113 + t75) + t142) + 6.0_dp * t103
      t80 = s25 * (s26 * t116 + s34 * (s34 * (s34 * (t96 * s16 * t141 + 
     #s26 * (15.0_dp * s25 + t136) + t102 + t129) - s25 * (s25 * (t131 +
     # t24) + t2 * t96) - s26 * (s26 * (t26 + 6.0_dp * t145) + s25 * (t5
     #0 + t92) + 36.0_dp * t2) - t124) - t23 * t62 * (s25 * (s16 - s25) 
     #+ t2) + t4 * t104 + (s26 * (s26 * t126 + t106) - t110 + t114 + t12
     #5) * s26) + s56 * (s56 * (s56 * (s25 * (t8 + t105 - s34) - t23 * t
     #33 + t84 * s34) + s34 * (s16 * (18.0_dp * s16 + t76) + s25 * (-t77
     # - t84 - t43)) + 6.0_dp * t29 * t53 - 6.0_dp * t4 + 6.0_dp * t63 -
     # t77 * t80 * t109 + t44 * t104) - t1 * (-t32 * (s26 * t15 + t120 +
     # t20) - t62 ** 2.0_dp + t111) + t15 * t93 * s26 + s34 * (-s34 * (s
     #25 * (t76 + t50 + t92) + 36.0_dp * t138 + 6.0_dp * t16) + s26 * (s
     #26 * t1 * t126 + t15 * t106 / 2.0_dp) - t110 + t114 + t125) - 6.0_
     #dp * t30 + 6.0_dp * t103 + t112 * (t105 + t25)))
      t92 = 0.1e1_dp / s16
      t90 = 0.1e1_dp / t90 ** 2.0_dp
      t46 = 0.1e1_dp / t46 ** 2.0_dp
      t89 = 0.1e1_dp / t89 ** 2.0_dp
      t1 = (-t14 * t2 * (t19 * (-s16 * t107 - s26 * t108) + t3 * t62) + 
     #t23 * s16 * (t16 * (t4 * t86 * t98 - t118 * t17) + t3 * t4) - t18 
     #- t1 * (t1 * ((s16 * (t84 + s25) - t16) * s26 + t112) + t16 * (-t8
     #3 + t105) + t155) * t31) * MB1011(1) * t92
      t1 = t89 * (t52 * MB1001(0) - t23 * t88 * (t99 + t101) + (t19 * (s
     #26 * t57 + t56) + t4 * ((t21 + t47) * s25 - t48)) * s26 + s34 * ((
     #((t78 * s34 - t14 * s16 * t45 - (t76 + t50 + t11) * s26) * s34 + s
     #16 * ((t74 + t24) * s25 + t75) + 18.0_dp * t19 + s26 * ((t10 + t72
     #) * s26 + t70 + t71)) * s34 - t58 * (((t8 + t61) * s25 + t22) * s2
     #5 + 12.0_dp * t19 + 10.0_dp * t4 + s26 * ((t67 + t68) * s26 + (s16
     # * t65 + t66) * s25 + 32.0_dp * t2))) * s34 + t115 * (t14 * t2 * t
     #45 + ((t61 + s25) * s26 + (t59 + t51) * s25 + t60) * s26 + t45 * t
     #62 + t64)) + t79 * s56 - t5 + t113 * t3 + t1) + (-t23 * s25 * (-t1
     #6 * t93 + t30 * (t77 + t83 + t105 - s34) + t6) + t80 + t137 * (t11
     #3 + t48)) * MB1110(1) * t92
      t3 = 0.1e1_dp / t91
      t3 = t3 * (t92 * (((s26 * t12 + s34 * (-t40 + s34) + (s16 - s25 - 
     #s34) * s56 + t82) * (t23 * (s16 * (-t112 + t69) - t87) - s34 * (s3
     #4 * (-s25 * (t119 + t77) - t38 - t39 - t29) + t58 * (s16 * (-t59 +
     # t84) + s26 * t109 + t62)) - s56 * (-s56 * (-s34 * (s25 - s34) + t
     #2 + t62 - t81) - t23 * (-t32 + t69 + t55) + s34 * (-s25 * t12 + s2
     #6 * t85 - s34 * (t26 + t8 + t51) + t121) + t112 * t15) + t49 + t54
     #) * MB1111D2(0) + s25 * s34 * (t23 * (s16 * (-s34 * t58 * (-s25 + 
     #s26 + t119) - t112 + t69) + s56 * (-t112 * t23 - s34 * (s16 * (t26
     # + t21 - s25) + s34 * t53) + t55 + t69) - t87) + t29 * (t23 * s25 
     #* t58 + t29 + t38 + t39) + t37 * (t29 + t2 + t62 - t81) + t49 + t5
     #4) * MB1101(0)) * t3 * t46 + t32 * MB1011(0) * t89 + t63 * MB1110(
     #0) * t90) - I300s16(0))
      result = t46 * ((t14 * t36 + t23 * (((-t8 - t7) * s25 + t2) * s26 
     #+ ((-s26 * t9 - t7 - t8) * s25 - 21.0_dp * t16 - 17.0_dp * t28 + t
     #33 - t34) * s56) - 14.0_dp * t41 + ((t26 + t25 + t24) * s34 + (t21
     # - t20) * s25 - t13 * t15 + 17.0_dp * t16 - t22) * s34 + (-t11 - t
     #10) * t16 + (-s26 * t35 - t10 - 17.0_dp * t27) * t37 - t44 * (t43 
     #+ t42)) * t92 + (t23 * t36 - 6.0_dp * t41 + ((-s16 * t73 - t95) * 
     #s26 + (-t8 - t94) * s25 + t2) * s26 + ((t83 + t51 + s26) * s34 + (
     #-t51 + s25) * s16 + t16 * t73 - t23 * (t62 + t13)) * s34 + ((-t27 
     #* t73 - t100 - t95) * s56 + (s26 * t97 - t8 - t94) * s25 - 18.0_dp
     # * t16 - 14.0_dp * t28 + t33 - t34) * s56 - t44 * (t8 + t105)) * I
     #300s16(1) + t1 * t90) / 12.0_dp + t3 / 6.0_dp

           intHs16s25s26s34s56x1141D6eps0 = result
       end function intHs16s25s26s34s56x1141D6eps0

       function intHs16s25s26s34s56x1141D6eps1()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1141D6eps1
           complex(dp) ::  t1,t10,t11,t12,t13,t14,t2,t3,t4,t5,t6,t7
           complex(dp) ::  t8,t9

           complex(dp) :: result

      t1 = s16 ** 2.0_dp
      t2 = 3.0_dp * s25
      t3 = 11.0_dp * s25
      t4 = s26 ** 2.0_dp
      t5 = t4 + t1
      t6 = s16 + s26
      t7 = s25 ** 2.0_dp
      t8 = 6.0_dp
      t9 = s34 ** 3.0_dp
      t10 = s16 + s26 + s56
      t11 = s16 + s25 + s26 - s34 + s56
      t12 = s16 + s26 - s34 + s56
      t13 = 0.1e1_dp / s16
      t12 = 0.1e1_dp / t12
      t11 = 0.1e1_dp / t11 ** 2.0_dp
      t10 = 0.1e1_dp / t10 ** 2.0_dp
      t14 = -0.1e1_dp / 0.12e2_dp
      result = t12 * ((s25 * t7 * MB1110(1) * t11 + t9 * MB1011(1) * t10
     #) * t13 - I300s16(1)) / 6.0_dp + t14 * (t8 * (s56 ** 3.0_dp + s16 
     #* t1 + s26 * t4) + ((18.0_dp * s16 + t3) * s26 + (22.0_dp * s16 + 
     #t2) * s25 + 18.0_dp * t1) * s26 + (-(4.0_dp * s25 + t6) * s34 - 7.
     #0_dp * t5 + 2.0_dp * (-s16 - s26 + s25) * s25 - 14.0_dp * s16 * s2
     #6) * s34 + ((-7.0_dp * s34 + t3 + 18.0_dp * t6) * s56 + 22.0_dp * 
     #s25 * t6 - 14.0_dp * s34 * t6 + (-2.0_dp * s25 - s34) * s34 + 18.0
     #_dp * t5 + 3.0_dp * t7 + 36.0_dp * s16 * s26) * s56 + 2.0_dp * t9 
     #+ s16 * s25 * (11.0_dp * s16 + t2)) * t13 * t10 * t11

           intHs16s25s26s34s56x1141D6eps1 = result
       end function intHs16s25s26s34s56x1141D6eps1

       function intHs16s25s26s34s56x1210D2eps0()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1210D2eps0
           complex(dp) ::  t1

           complex(dp) :: result

      t1 = s16 + s25 + s26 - s34 + s56
      t1 = 0.1e1_dp / t1
      result = -(-(s16 + s26 + s56) * MB1011(1) - (s16 + s26 - s34 + s56
     #) * MB1110(0) - (s16 - s25 + s26 - s34 + s56) * MB1110(1) + MB1001
     #(0)) * t1

           intHs16s25s26s34s56x1210D2eps0 = result
       end function intHs16s25s26s34s56x1210D2eps0

       function intHs16s25s26s34s56x1210D2eps1()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1210D2eps1
           complex(dp) ::  t1

           complex(dp) :: result

      t1 = s16 + s25 + s26 - s34 + s56
      t1 = 0.1e1_dp / t1
      result = t1 * ((s16 + s26 - s34 + s56) * MB1110(1) - 1.0_dp)

           intHs16s25s26s34s56x1210D2eps1 = result
       end function intHs16s25s26s34s56x1210D2eps1

       function intHs16s25s26s34s56x1211D2eps0()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1211D2eps0
           complex(dp) ::  t1

           complex(dp) :: result

      t1 = s34 - s16 - s26 - s56
      t1 = 0.1e1_dp / t1
      result = -0.1e1_dp / s16 * ((s34 * MB1101(0) - (s16 + s26 + s56) *
     # MB1111D2(0)) * t1 - MB1110(0))

           intHs16s25s26s34s56x1211D2eps0 = result
       end function intHs16s25s26s34s56x1211D2eps0

       function intHs16s25s26s34s56x1211D2eps1()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1211D2eps1

           complex(dp) :: result

      result = MB1110(1) / s16

           intHs16s25s26s34s56x1211D2eps1 = result
       end function intHs16s25s26s34s56x1211D2eps1

       function intHs16s25s26s34s56x1211D4eps0()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1211D4eps0
           complex(dp) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           complex(dp) ::  t20,t3,t4,t5,t6,t7,t8,t9

           complex(dp) :: result

      t1 = s16 - s25
      t2 = s16 + s25
      t3 = s16 ** 2.0_dp
      t4 = t2 * s26
      t5 = 2.0_dp * t3
      t6 = s25 ** 2.0_dp
      t7 = s16 * s25
      t8 = s34 ** 2.0_dp
      t9 = t6 + t3 + t8
      t10 = t1 ** 2.0_dp * s26
      t11 = t5 * s25
      t5 = -t11 + t9 * s16 - (t5 + t4) * s34 + (-s16 * s34 - s25 * s34 +
     # t3 + t6 - 2.0_dp * t7) * s56 + t10
      t12 = -2.0_dp * s34 * t2 - 2.0_dp * t7 + t9
      t13 = 2.0_dp * s26
      t14 = (-s34 + t2 + s26) * s26 + (t13 + s16 + s25 - s34 + s56) * s5
     #6 + t7
      t15 = s16 + s26 + s56
      t16 = 3.0_dp * s16
      t17 = (t6 + t3) * s16
      t18 = t1 * s16
      t12 = 0.1e1_dp / t12
      t19 = 0.1e1_dp / t14
      t14 = 0.1e1_dp / t14
      t20 = -0.1e1_dp / 0.2e1_dp
      result = t20 * (-t15 * (s16 + s26 - s34 + s56) * MB1111D2(0) * t19
     # + (((2.0_dp * s34 * s56 ** 2.0_dp - 2.0_dp * s25 * t3 + (2.0_dp *
     # s26 ** 2.0_dp + (s16 - t13) * s34 - 2.0_dp * t18 + t4) * s34 + (-
     #2.0_dp * t7 - 2.0_dp * t8 + (4.0_dp * s26 + t2) * s34 + t3 + t6) *
     # s56 + t10 + t17) * MB1110(1) + s34 * (-s26 * t1 + (s16 - s26) * s
     #34 - (s16 - s25 + s34) * s56 - t18) * MB1101(0)) * s25 + t5 * (-s1
     #6 * I300s16(1) + MB1001(0)) - t15 * (-t11 + (-s16 * (-s25 + t13 + 
     #t16) + (-s34 + t16 + s25 + s26) * s34) * s34 + (t9 - 2.0_dp * s16 
     #* (s25 + s34)) * s56 + t10 + t17) * MB1011(1)) * t14 * t12) + t5 *
     # t12 * t14

           intHs16s25s26s34s56x1211D4eps0 = result
       end function intHs16s25s26s34s56x1211D4eps0

       function intHs16s25s26s34s56x1211D4eps1()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1211D4eps1

           complex(dp) :: result

      result = 0.0_dp

           intHs16s25s26s34s56x1211D4eps1 = result
       end function intHs16s25s26s34s56x1211D4eps1

       function intHs16s25s26s34s56x1212D4eps0()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1212D4eps0
           complex(dp) ::  t1,t10,t11,t12,t2,t3,t4,t5,t6,t7,t8,t9

           complex(dp) :: result

      t1 = s16 - s25
      t2 = 2.0_dp * s16
      t3 = t1 * s16
      t1 = t1 * s26
      t4 = (t2 + s25 + s26 - s34) * s34 - (s16 - s25 - s34) * s56 - t3 -
     # t1
      t5 = s16 + s25
      t6 = s16 * s25
      t7 = s34 ** 2.0_dp
      t8 = s16 ** 2.0_dp
      t9 = s25 ** 2.0_dp - 2.0_dp * s34 * t5 - 2.0_dp * t6 + t7 + t8
      t10 = s16 + s25 - s34
      t11 = (-s34 + t5 + s26) * s26 + (2.0_dp * s26 + t10 + s56) * s56 +
     # t6
      t12 = s25 * s34
      t9 = 0.1e1_dp / t9
      t11 = 0.1e1_dp / t11
      result = 2.0_dp * t4 * t9 * t11 - t11 * ((s16 + s26 - s34 + s56) *
     # MB1111D2(0) + t9 * (t4 * (-s16 * I300s16(1) + MB1001(0)) + (s25 *
     # (s16 + s25 + s26 - s34 + s56) * ((s16 - s26) * s34 - (s16 - s25 +
     # s34) * s56 - t1 - t3) * MB1110(1) - (s16 + s26 + s56) * ((-t1 - t
     #3) * s16 + ((-3.0_dp * s16 - 2.0_dp * s25 - s26 + s34) * s34 + 3.0
     #_dp * t8 + s25 * t5 + (s25 + t2) * s26) * s34 - (-t2 * s34 - t12 -
     # t6 + t7 + t8) * s56) * MB1011(1)) / s16 + t12 * (t10 + 2.0_dp * s
     #26 + 2.0_dp * s56) * MB1101(0)))

           intHs16s25s26s34s56x1212D4eps0 = result
       end function intHs16s25s26s34s56x1212D4eps0

       function intHs16s25s26s34s56x1212D4eps1()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1212D4eps1

           complex(dp) :: result

      result = 0.0_dp

           intHs16s25s26s34s56x1212D4eps1 = result
       end function intHs16s25s26s34s56x1212D4eps1

       function intHs16s25s26s34s56x1213D6eps0()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1213D6eps0
           complex(dp) ::  t1,t10,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109
           complex(dp) ::  t11,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t12
           complex(dp) ::  t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t13,t130
           complex(dp) ::  t131,t132,t133,t134,t135,t136,t137,t138,t139,t14,t140,t141
           complex(dp) ::  t142,t143,t144,t145,t146,t147,t148,t149,t15,t150,t151,t152
           complex(dp) ::  t153,t154,t155,t156,t157,t158,t159,t16,t160,t161,t162,t163
           complex(dp) ::  t164,t165,t166,t167,t168,t169,t17,t170,t171,t172,t173,t174
           complex(dp) ::  t175,t176,t177,t18,t19,t2,t20,t21,t22,t23,t24,t25
           complex(dp) ::  t26,t27,t28,t29,t3,t30,t31,t32,t33,t34,t35,t36
           complex(dp) ::  t37,t38,t39,t4,t40,t41,t42,t43,t44,t45,t46,t47
           complex(dp) ::  t48,t49,t5,t50,t51,t52,t53,t54,t55,t56,t57,t58
           complex(dp) ::  t59,t6,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69
           complex(dp) ::  t7,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t8
           complex(dp) ::  t80,t81,t82,t83,t84,t85,t86,t87,t88,t89,t9,t90
           complex(dp) ::  t91,t92,t93,t94,t95,t96,t97,t98,t99

           complex(dp) :: result

      t1 = s16 - s25
      t2 = 2.0_dp
      t3 = s25 ** 2.0_dp
      t4 = t3 ** 2.0_dp
      t5 = s25 * t4
      t6 = s25 * t3
      t7 = t2 * t3
      t8 = 3.0_dp * s16 * t1 - t7
      t9 = 3.0_dp * s16
      t10 = -t9 - s25
      t11 = s16 ** 2.0_dp
      t12 = t11 ** 2.0_dp
      t13 = s16 * t11
      t14 = s16 * t12
      t15 = t2 * t11
      t16 = 3.0_dp * s25
      t17 = s16 + s26
      t18 = -t3 + t11
      t19 = s16 * s25
      t20 = 8.0_dp
      t21 = 7.0_dp * t11
      t22 = 10.0_dp * t12
      t23 = t20 * s16
      t24 = t2 * t13
      t25 = ((t23 - s25) * s25 + t21) * s25 + t24
      t26 = s26 ** 2.0_dp
      t27 = s26 * t26
      t28 = t11 + t3
      t29 = t28 * t3
      t30 = 24.0_dp * s16
      t31 = t2 * s26
      t32 = 5.0_dp * s25
      t33 = 31.0_dp * s16
      t34 = 53.0_dp * t11
      t35 = 22.0_dp * t13
      t36 = 85.0_dp * t12
      t37 = s16 + s25
      t38 = 37.0_dp * t19 * t37
      t39 = 28.0_dp * t13
      t40 = t20 * t6 - t38 - t39
      t41 = s25 * (-9.0_dp * s16 + t16)
      t42 = 23.0_dp
      t43 = 21.0_dp * s16
      t44 = 6.0_dp * t11
      t45 = 60.0_dp * t13
      t46 = 39.0_dp
      t47 = 10.0_dp * s25
      t48 = ((s16 * t46 - t47) * s25 + 66.0_dp * t11) * s25
      t49 = 80.0_dp * t13
      t50 = -12.0_dp * t18
      t51 = 17.0_dp * t19
      t52 = t51 - t50
      t53 = t2 * s16
      t54 = -t16 - t53
      t55 = 37.0_dp * s16
      t56 = (-t43 + t47) * s25
      t57 = -35.0_dp * t11
      t58 = 4.0_dp
      t59 = t58 * s25
      t60 = t59 + s16
      t61 = -49.0_dp
      t62 = 20.0_dp * s25
      t63 = 15.0_dp
      t64 = t58 * s16
      t65 = t63 * s25
      t66 = t11 * t28
      t67 = 30.0_dp * s16
      t68 = 9.0_dp * s25
      t69 = 16.0_dp * s16
      t70 = 25.0_dp * t11
      t71 = 3.0_dp * t26
      t72 = 3.0_dp * s26
      t73 = t1 ** 2.0_dp
      t74 = s16 * t4
      t75 = t73 * (s25 * t10 + t15)
      t76 = t11 * t6
      t77 = (-t16 + s16) * t73
      t78 = 27.0_dp * s16
      t79 = t12 * ((-21.0_dp * s25 + t78) * s25 - 9.0_dp * t11)
      t80 = t2 * s25
      t81 = 9.0_dp * t13
      t82 = t20 * s25
      t83 = 9.0_dp * s26
      t84 = 6.0_dp * s16
      t85 = t13 * s25
      t86 = t85 * t1
      t87 = s34 ** 2.0_dp
      t88 = t87 ** 2.0_dp
      t89 = t87 * t88
      t90 = s34 * t88
      t91 = s34 * t87
      t92 = 7.0_dp * s16
      t93 = 5.0_dp * t11
      t38 = ((((t54 + s34) * s34 + t41) * s34 + t25) * s34 - t12 + t19 *
     # ((-t92 + t16) * s25 + t93)) * s56 - t2 * t90 - 6.0_dp * t14 + (((
     #(t82 + t72 + t53) * s34 + (-t84 - t68) * s26 - t50 + t51) * s34 - 
     #((-t83 - t82) * s25 + t78 * s26) * s25 - t38 - t39) * s34 + s25 * 
     #(t3 * (s16 * t63 - t80) - t81) + 22.0_dp * t66 + t72 * t25) * s34 
     #- 3.0_dp * s16 * (t6 * (-s16 - s25) + t77 * s26) + 21.0_dp * t86
      t38 = t38 * s56 + (((((t58 * (s16 - s26) - t32 + s34) * s34 + (t58
     # * t60 + t72) * s26 + t56 + t57) * s34 + s26 * (3.0_dp * s26 * t54
     # + t2 * t52) + t49 + t48) * s34 + t2 * t40 * s26 - t35 * s25 - t34
     # * t3 - t33 * t6 - t36 + 5.0_dp * t4 - t68 * (t9 - s25) * t26) * s
     #34 + s16 * ((((t69 + t68) * s25 + t70) * s25 - 45.0_dp * t13) * s2
     #5 + 44.0_dp * t12) + (s25 * (t3 * (t67 - t59) - 18.0_dp * t13) + 4
     #4.0_dp * t66) * s26 - t5 + t71 * t25) * s34 - 6.0_dp * s16 * (t75 
     #* s26 - t74) - 3.0_dp * s16 * (t77 * t26 + t76) + t79
      t39 = s16 * t73
      t50 = 7.0_dp * s25
      t51 = 11.0_dp * s16
      t66 = t58 * t11
      t94 = t87 + t11 + t3
      t95 = -t2 * (s34 * t37 + t19) + t94
      t96 = (-s34 + t37 + s26) * s26 + (s16 + s25 - s34 + t31 + s56) * s
     #56 + t19
      t97 = s25 * t37
      t98 = 5.0_dp * s16
      t99 = (-t98 - t80) * s25 + t66
      t100 = -t80 + s16
      t101 = 6.0_dp * s25
      t102 = 10.0_dp * t11
      t103 = 13.0_dp * t13
      t104 = t19 * t2 + 3.0_dp * t28
      t105 = 3.0_dp * t11
      t106 = 5.0_dp * t3
      t107 = ((t33 + t62) * s25 + 18.0_dp * t11) * s25
      t108 = 50.0_dp * t13
      t109 = s16 * t42
      t110 = 22.0_dp * s25
      t111 = 24.0_dp * t11
      t112 = (t109 + t110) * s25 + t111
      t113 = t80 + s16
      t114 = t13 + t6
      t115 = 12.0_dp * s16
      t116 = 13.0_dp * s25
      t117 = 3.0_dp * t27
      t118 = 46.0_dp * s16
      t119 = (30.0_dp * s25 + t118) * s25
      t120 = 50.0_dp * t11
      t121 = 40.0_dp * t13
      t122 = 25.0_dp * s16
      t123 = t58 * s26
      t124 = 12.0_dp * s25
      t125 = 3.0_dp * t104 * s26
      t126 = 32.0_dp * s16
      t127 = t20 * s26
      t128 = s16 * t6
      t129 = t72 * t100
      t130 = t73 * s26
      t131 = 14.0_dp * t12
      t132 = 10.0_dp * s16
      t133 = 16.0_dp * t13
      t134 = t88 + t12
      t110 = ((s25 * ((t98 - t80) * s25 - t66) - (t104 + t87) * s34 + t1
     #3 + 3.0_dp * t113 * t87) * s56 + s25 * ((12.0_dp * t11 - t7) * s25
     # - t103) + t58 * t134 - (((18.0_dp * s25 + t72 + t69) * s34 - s16 
     #* (t83 + t30) - (18.0_dp * s26 + t109 + t110) * s25) * s34 + ((t13
     #2 + t101) * s25 - t11 * t20) * s25 + t125 + t133) * s34 - t128 + t
     #129 * t73) * s56 + t58 * t128 * t1 - (((-(t127 + t62 + t122) * s34
     # + (36.0_dp * s25 + t72 + t126) * s26 + t119 + t120) * s34 - (t2 *
     # t112 + t83 * t113) * s26 - t107 - t108) * s34 + (((t84 + t32) * s
     #25 + t102) * s25 - t35) * s25 + (t58 * (((t98 + t16) * s25 - t66) 
     #* s25 + t13 * t20) + t125) * s26 + 25.0_dp * t12) * s34 - 5.0_dp *
     # t90 + 5.0_dp * t14 + t130 * (t2 * t99 + t129) + t81 * t3 - t131 *
     # s25
      t125 = 3.0_dp * t3
      t129 = -t92 - t80
      t135 = s25 * t129 + t93
      t10 = t2 * t11 * t1 - t10 * t3
      t136 = t10 * t26
      t137 = 6.0_dp * t3
      t138 = 11.0_dp * s25
      t139 = 40.0_dp * s16
      t140 = 69.0_dp * t12
      t141 = 16.0_dp * t11
      t142 = (t32 + t64) * s25 + t44
      t143 = t142 * s26
      t144 = (((58.0_dp * s16 + t32) * s25 + 109.0_dp * t11) * s25 + 124
     #.0_dp * t13) * s25
      t145 = 155.0_dp * t12
      t146 = 57.0_dp * t11
      t147 = ((t118 + t124) * s25 + t146) * s25 + t45
      t148 = s16 * t37
      t149 = 38.0_dp * t148
      t150 = 21.0_dp * t3
      t151 = t150 + t149
      t152 = t58 * t37
      t153 = ((76.0_dp * s16 + t47) * s25 + 131.0_dp * t11) * s25
      t154 = 125.0_dp * t13
      t155 = 75.0_dp * t148
      t156 = 26.0_dp * t3
      t157 = t156 + t155
      t158 = 17.0_dp * s16
      t159 = t116 + t158
      t160 = 14.0_dp * s25
      t161 = (50.0_dp * s16 + t47) * s25
      t162 = 16.0_dp * s25
      t163 = 13.0_dp * s16
      t164 = 14.0_dp * s16
      t127 = t127 * t10
      t165 = (((55.0_dp * s16 + t138) * s25 + t146) * s25 + 32.0_dp * t1
     #3) * s25
      t166 = 105.0_dp * t12
      t167 = 20.0_dp * t11
      t168 = ((t65 + t55) * s25 + t167) * s25
      t169 = 42.0_dp * t13
      t170 = 11.0_dp * t12
      t171 = t123 * t142
      t172 = 3.0_dp * t169 + 3.0_dp * t168
      t173 = 6.0_dp * s26
      t174 = t73 * t26
      t175 = 13.0_dp * t11
      t176 = s26 * t37
      t177 = 5.0_dp * t12
      t67 = (((((t53 + s26) * s34 - s16 * (t164 + t68) - (t163 + t32 + t
     #72) * s26) * s34 + s16 * ((t139 + t162) * s25 + 42.0_dp * t11) + s
     #26 * ((t160 + t30) * s26 + t146 + t161) + t117) * s34 - s16 * (((t
     #118 + t160) * s25 + 65.0_dp * t11) * s25 + 70.0_dp * t13) - (((t15
     #9 + s26) * s26 + t157) * s26 + t153 + t154) * s26) * s34 + s16 * (
     #(((t67 + t101) * s25 + 44.0_dp * t11) * s25 + t121) * s25 + 70.0_d
     #p * t12) + t26 * (t147 * t2 + t152 * t26) + s26 * (t151 * t26 + t1
     #44 + t145)) * s34 - t17 * (s16 * ((((t51 + t124) * s25 + t141) * s
     #25 - 5.0_dp * t13) * s25 + 42.0_dp * t12) + ((((t126 + t65) * s25 
     #+ t141) * s25 + 36.0_dp * t13 + t143) * s26 + (((t139 + t138) * s2
     #5 + t70) * s25 + t133) * s25 + t140) * s26 + t5)
      t10 = (t2 * (s34 * t10 + t128) + t58 * (t37 * t91 + t85) - t87 * (
     #t142 + t87) - t12 - t93 * t3) * s56 + (-t130 * t100 * t58 + t2 * t
     #4 - t177) * s16 + 3.0_dp * t76 + 17.0_dp * t86 + 3.0_dp * t90 + ((
     #(-(t116 + t123 + t158) * s34 + t149 + t150 + 16.0_dp * t176) * s34
     # - t168 - t169 - t171) * s34 + (((t164 + t59) * s25 + t175) * s25 
     #- t35) * s25 + t12 * t42 + t127) * s34
      t10 = t10 * s56 - 3.0_dp * (t130 * t135 + t76) * s16 - 3.0_dp * t8
     #9 - 6.0_dp * s16 * (t174 * t100 - t74) + (((((t160 + t83 + t30) * 
     #s34 - (s25 * t46 + 51.0_dp * s16 + t173) * s26 - t155 - t156) * s3
     #4 + s16 * ((114.0_dp * s16 + 92.0_dp * s25) * s25 + 120.0_dp * t11
     #) + 24.0_dp * t26 * t37 + 24.0_dp * t6 + t72 * t151) * s34 - (6.0_
     #dp * t143 + t172) * s26 - t165 - t166) * s34 + ((((42.0_dp * s16 +
     # t124) * s25 + t11 * t46) * s25 - 66.0_dp * t13) * s25 + t140) * s
     #26 + t6 * ((t164 + t80) * s25 + t167) + 12.0_dp * t13 * t3 + 12.0_
     #dp * t136 + 48.0_dp * t12 * t1) * s34 + t79
      t10 = t10 * s56 + s16 * (19.0_dp * t14 * s25 - t170 * t3 - t174 * 
     #(t123 * t100 + 3.0_dp * t135)) + ((((((-t163 - t32 - t173 + s34) *
     # s34 + (48.0_dp * s16 + 28.0_dp * s25 + t83) * s26 + t146 + t161) 
     #* s34 - ((3.0_dp * t159 + t123) * s26 + t157 * t2) * s26 - t153 - 
     #t154) * s34 + ((16.0_dp * s26 * t37 + 3.0_dp * t151) * s26 + t147 
     #* t58) * s26 + t144 + t145) * s34 - ((((t109 + s25) * s25 + 51.0_d
     #p * t11) * s25 + 41.0_dp * t13) * s25 + t170) * s25 - ((t172 + t17
     #1) * s26 + t2 * (t166 + t165)) * s26 - 111.0_dp * t14) * s34 + t17
     # * (((((t69 + t59) * s25 + t44) * s25 + t81) * s25 - 46.0_dp * t12
     #) * s25 + ((((34.0_dp * s16 + t124) * s25 + t11 * t63) * s25 - t10
     #8) * s25 + 53.0_dp * t12 + t127) * s26 + 43.0_dp * t14)) * s34 - 6
     #.0_dp * t11 * (t130 * t8 - t74) - 7.0_dp * t12 * t114
      t10 = -t2 * t12 * (t3 * (-t3 + t11) + t12) + (-t73 * t100 * t26 **
     # 2.0_dp - t73 * t135 * t27) * s16 + s34 * (t67 * s34 + t17 ** 2.0_
     #dp * (s16 * (((t137 + t66) * s25 - t133) * s25 + t131) + ((((t132 
     #+ t59) * s25 + t11) * s25 - 14.0_dp * t13) * s25 + t12 * t63) * s2
     #6 + t2 * (-t18 * t6 + t136))) + t10 * s56 - t130 * t13 * ((-t101 -
     # t98) * s25 + t21) - 3.0_dp * t11 * (t13 * t6 + t174 * t8) + t32 *
     # t13 * t12
      t37 = t20 * t114
      t46 = (t138 + s16) * s25 + t105
      t28 = 12.0_dp * t28
      t63 = 11.0_dp * t19
      t67 = t50 + t9
      t70 = t125 * (s16 - s25 + s34) * s34 * s56 ** 3.0_dp
      t7 = t70 - s25 * s34 * (((((t98 + t59 + t31 - s34) * s34 - t176 * 
     #t20 - t137 - 10.0_dp * t148 - t26) * s34 + ((t92 + t59) * s25 + t4
     #4) * s25 + 10.0_dp * t13 + s26 * (s26 * t67 + t28 + t63)) * s34 - 
     #(t64 * t3 + t117 - t24) * s25 - t177 - (t46 * s26 + t19 * (t50 - t
     #53) + t37) * s26 - t29) * s34 + (((-t83 * t1 + t106 - t44) * s25 +
     # ((t67 - s34) * s34 - (t138 + t83 + s16) * s25 - t105) * s34 + t13
     #) * s56 + s25 * (s16 * ((t59 - t9) * s25 - t93) - 9.0_dp * t1 * t2
     #6) + t2 * (((t106 - t44) * s25 + t13) * s26 + t91 * (-t152 - s26 +
     # s34) + t12 + t4) - s34 * (t2 * t46 * s26 - t15 * s25 - ((t84 + t1
     #60) * s26 + t28 + t63) * s34 + t68 * t26 + t92 * t3 + t37)) * s56 
     #+ t1 * t17 * ((-t71 - t7) * s25 + s16 * (-t97 + t11) + ((-t32 - t5
     #3) * s25 + t11) * s26))
      t28 = s16 + s25 + s26 - s34 + s56
      t15 = (-t64 - s25) * s25 + t15
      t37 = t98 - s25
      t46 = (t69 - t16) * s25 + t44
      t63 = s25 * s26
      t53 = t53 - s25
      t6 = (-t2 * s16 * (t13 - t6) + 3.0_dp * (s25 * ((10.0_dp * s26 + t
     #9 + t59) * s16 - t3) + t2 * (-s26 * t3 + t13)) * s34 + 3.0_dp * t1
     #30 * s25 + t87 * (t53 * s34 - (t162 + t84) * s16 + 3.0_dp * s25 * 
     #(s25 + s26)) + t4 - t11 * s25 * (t68 - t23) + s25 * (-t2 * s25 * (
     #s16 + s34) + t94 + t132 * s34) * s56) * s56 - s25 * (t13 * (t124 -
     # t163) - t71 * t73) - t2 * (t130 * t15 - t74) - t58 * s16 * t134 +
     # ((-t2 * t46 * s26 - t175 * s25 + t71 * s25 + ((t69 + t116 + t123)
     # * s16 - t31 * s25) * s34 - t115 * t3 - 24.0_dp * t13) * s34 + s16
     # * (((20.0_dp * s16 + s25) * s25 - t175) * s25 + t133) + 6.0_dp * 
     #s26 * (t63 * t37 + ((t64 - s25) * s25 + t105) * s25 + t24)) * s34 
     #+ t76
      t9 = s16 + s26 - s34 + s56
      t24 = 0.1e1_dp / t95 ** 2.0_dp
      t59 = 0.1e1_dp / t96 ** 2.0_dp
      t4 = s25 * t28 * (t2 * s16 * (t88 * (t98 + t16 + t31 - s34) + t14)
     # - (((s16 * ((t101 + t115) * s25 + t167) + s26 * (s16 * (t69 + t11
     #6) + s26 * t53)) * s34 - s16 * (-t129 * t3 + 20.0_dp * t13) - (s16
     # * ((t124 + t163) * s25 + t111) + (-t63 + t46) * s26) * s26) * s34
     # + t17 * ((((t164 - t16) * s25 - t11) * s25 + 6.0_dp * t13) * s26 
     #+ t22 + t44 * t3 + t128 * t58 - 12.0_dp * t85 + t80 * t37 * t26)) 
     #* s34 - t6 * s56 - t11 * t4 - t73 * s25 * t27 + t174 * t15 + t39 *
     # t99 * s26 + t12 * s25 * (-t84 + t32)) * MB1110(1)
      t4 = t10 * MB1011(1) + t14 * ((t51 - t50) * s25 - t66) + ((((((-t6
     #4 + s26) * s34 + s16 * (t65 + t30) + (-t32 + t64 - t31) * s26) * s
     #34 + s16 * ((s16 * t61 - t62) * s25 - 60.0_dp * t11) + s26 * (t57 
     #+ t56 + t26) + t2 * t60 * t26) * s34 + s16 * (((t55 + t47) * s25 +
     # 46.0_dp * t11) * s25 + t49) + ((s26 * t54 + t52) * s26 + t48 + t4
     #9) * s26) * s34 + t11 * (((-s25 * t42 - t43) * s25 + t44) * s25 - 
     #t45) + ((s26 * t41 + t40) * s26 - (((t33 - t32) * s25 + t34) * s25
     # + t35) * s25 - t36) * s26) * s34 + t17 * (11.0_dp * t29 * s16 + t
     #12 * (-29.0_dp * s25 + t30) + t25 * t26 - t5 + t31 * (-t18 * t19 *
     # t20 + t3 * (t21 - t3) + t22))) * s34 + t38 * s56 - t77 * s16 * t2
     #7 - 3.0_dp * s16 * (t39 * t8 * s26 + t76 * (s16 - s25) + t75 * t26
     #) + t4
      result = t24 * t59 * ((-s16 * I300s16(1) + MB1001(0)) * (t2 * (t11
     # * (t18 * t3 + t12) + t89) + ((((-(5.0_dp * s26 + t82 + t115) * s3
     #4 + (t78 + t124) * s25 + (t62 + t123 + t122) * s26 + 30.0_dp * t11
     #) * s34 - ((t109 + t82) * s25 + 28.0_dp * t11) * s25 - s26 * (t120
     # + t119 + t26) - t121 - t2 * (t68 + t23) * t26) * s34 + s16 * (t3 
     #* (t116 + t115) + 30.0_dp * t13) + t2 * s25 * t114 + s26 * (s26 * 
     #t112 + t107 + t108) + t117 * t113) * s34 - t17 * (t3 * (t106 + t10
     #5) + s26 * (t104 * s26 + ((t92 + t101) * s25 - t102) * s25 + t103)
     # + 12.0_dp * t13 * t1)) * s34 + t110 * s56 + t130 * (s16 * (-t58 *
     # t97 + t93) + (s26 * t100 + t99) * s26) + t85 * (t125 - t93)) + t4
     # / s16) / 4.0_dp + t59 * (t7 * MB1101(0) * t24 - t9 ** 2.0_dp * t2
     #8 * MB1111D2(0)) / 2.0_dp

           intHs16s25s26s34s56x1213D6eps0 = result
       end function intHs16s25s26s34s56x1213D6eps0

       function intHs16s25s26s34s56x1213D6eps1()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1213D6eps1

           complex(dp) :: result

      result = 0.0_dp

           intHs16s25s26s34s56x1213D6eps1 = result
       end function intHs16s25s26s34s56x1213D6eps1

       function intHs16s25s26s34s56x1220D4eps0()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1220D4eps0
           complex(dp) ::  t1,t2

           complex(dp) :: result

      t1 = s16 + s25 + s26 - s34 + s56
      t2 = s16 + s26 - s34 + s56
      t1 = 0.1e1_dp / t1 ** 2.0_dp
      result = -t1 * (((s16 + s26 + s56) * MB1011(1) - MB1001(0) - 1.0_d
     #p) * (s16 - s25 + s26 - s34 + s56) + s25 * (s25 - 5.0_dp * t2) * M
     #B1110(1)) / 2.0_dp + s25 * t2 * MB1110(0) * t1

           intHs16s25s26s34s56x1220D4eps0 = result
       end function intHs16s25s26s34s56x1220D4eps0

       function intHs16s25s26s34s56x1220D4eps1()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1220D4eps1
           complex(dp) ::  t1,t2

           complex(dp) :: result

      t1 = s16 + s25 + s26 - s34 + s56
      t2 = 0.1e1_dp / 0.2e1_dp
      t1 = 0.1e1_dp / t1 ** 2.0_dp
      result = t1 * (t2 * (s16 - s25 + s26 - s34 + s56) + s25 * (s16 + s
     #26 - s34 + s56) * MB1110(1))

           intHs16s25s26s34s56x1220D4eps1 = result
       end function intHs16s25s26s34s56x1220D4eps1

       function intHs16s25s26s34s56x1221D4eps0()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1221D4eps0
           complex(dp) ::  t1,t10,t11,t12,t13,t2,t3,t4,t5,t6,t7,t8
           complex(dp) ::  t9

           complex(dp) :: result

      t1 = s16 + s26 + s56
      t2 = 2.0_dp * s26
      t3 = s16 * s25
      t4 = (-s34 + s16 + s25 + s26) * s26 + (t2 + s16 + s25 - s34 + s56)
     # * s56 + t3
      t5 = s16 + s26 - s34 + s56
      t6 = s16 + s25 + s26 - s34 + s56
      t7 = s16 + s26
      t8 = s26 ** 2.0_dp
      t9 = s16 - s25
      t10 = t9 * s26
      t11 = 0.1e1_dp / s16
      t12 = 0.1e1_dp / t5
      t6 = 0.1e1_dp / t6
      t4 = 0.1e1_dp / t4
      t13 = t11 * t4
      result = 2.0_dp * t13 * t1 - t6 * (t4 * (t5 * MB1001(0) - (s25 * (
     #s25 * s26 + (s16 - s26) * s34 + (t2 + s25 - s34 + s56) * s56 + t3 
     #- s16 ** 2.0_dp + t8) * MB1110(1) + t1 * (t7 * s16 + (-2.0_dp * s1
     #6 - s25 - s26 + s34) * s34 + (s16 - s34) * s56) * MB1011(1)) * t11
     #) - s25 * MB1110(0) * t11) + t1 * I300s16(1) * t4 + t13 * (((s16 *
     # t9 + s34 ** 2.0_dp) * s16 + t9 * (s56 ** 2.0_dp + t8) + 2.0_dp * 
     #(-s34 * t7 + t10) * s16 + 2.0_dp * ((-s25 - s34 + s16) * s16 + t10
     #) * s56) * MB1111D2(0) + s25 * s34 * t1 * MB1101(0)) * t12

           intHs16s25s26s34s56x1221D4eps0 = result
       end function intHs16s25s26s34s56x1221D4eps0

       function intHs16s25s26s34s56x1221D4eps1()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1221D4eps1
           complex(dp) ::  t1

           complex(dp) :: result

      t1 = s16 + s25 + s26 - s34 + s56
      t1 = 0.1e1_dp / t1
      result = 0.1e1_dp / s16 * t1 * (s25 * MB1110(1) + 1.0_dp)

           intHs16s25s26s34s56x1221D4eps1 = result
       end function intHs16s25s26s34s56x1221D4eps1

       function intHs16s25s26s34s56x1222D6eps0()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1222D6eps0
           complex(dp) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           complex(dp) ::  t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30
           complex(dp) ::  t31,t32,t33,t34,t35,t36,t37,t38,t39,t4,t40,t41
           complex(dp) ::  t42,t43,t44,t45,t46,t47,t48,t49,t5,t50,t51,t52
           complex(dp) ::  t53,t54,t55,t56,t57,t58,t59,t6,t60,t61,t62,t63
           complex(dp) ::  t64,t65,t66,t67,t68,t69,t7,t70,t71,t72,t73,t8
           complex(dp) ::  t9

           complex(dp) :: result

      t1 = 2.0_dp
      t2 = t1 * s25
      t3 = 7.0_dp * s16
      t4 = s16 - s25
      t5 = t1 * s16
      t6 = -s25 + t5
      t7 = s16 + s26
      t8 = s16 ** 2.0_dp
      t9 = t8 ** 2.0_dp
      t10 = s16 * t8
      t11 = s25 * t4
      t12 = s16 + s25
      t13 = s25 ** 2.0_dp
      t14 = s25 * t13
      t15 = s16 * t13
      t16 = 11.0_dp
      t17 = 4.0_dp
      t18 = s26 ** 2.0_dp
      t19 = 3.0_dp
      t20 = 18.0_dp
      t21 = 6.0_dp * s16
      t22 = t19 * s25
      t23 = 5.0_dp * s25
      t24 = 5.0_dp * s16
      t25 = (t24 + t2) * s25
      t26 = t1 * s26
      t27 = t17 * s16
      t28 = 7.0_dp * s25
      t29 = 16.0_dp * s16
      t30 = s16 + t2
      t31 = t12 * s26
      t32 = 7.0_dp * t8
      t33 = t13 + t18
      t34 = 6.0_dp * s25
      t35 = t17 * s26
      t36 = t19 * t33
      t37 = t4 ** 2.0_dp
      t38 = s16 * t14
      t39 = t37 * s26
      t40 = s25 * t16
      t41 = t37 * t18
      t42 = t41 * t19
      t43 = t19 * t31
      t44 = t19 * s26
      t45 = s34 ** 2.0_dp
      t46 = t45 ** 2.0_dp
      t47 = s34 * t45
      t48 = t17 * s25
      t49 = t39 * t19
      t50 = s16 * s25
      t51 = -t1 * (s34 * t12 + t50) + t8 + t13 + t45
      t52 = t8 * t14
      t37 = t37 * s26 * t18
      t53 = t39 * s16
      t54 = t17 * t8
      t6 = t10 * ((-9.0_dp * s16 + t34) * s25 + t54) - (((-(-s26 + t27) 
     #* s34 + s16 * (t28 + t29) - (-t27 + t22 + t26) * s26) * s34 - s16 
     #* (24.0_dp * t8 + t25) + ((t23 + t5 + s26) * s26 - (-t22 + t21) * 
     #s25 - t20 * t8) * s26) * s34 + t7 * ((-t16 * t8 - t13) * s25 + 16.
     #0_dp * t10 + t17 * (-(-t11 - t8) * s26 + t15) - t1 * t12 * t18)) *
     # s34 - ((t51 * s56 - t1 * ((-s16 * t12 + t1 * t13 + t43 + t45) * s
     #34 + t10) + (t23 + t44 + t5) * t45 + t14 + t49 - t50 * (t48 - t24)
     #) * s56 + t1 * (-t39 * t6 + t38) - t8 * ((-t29 + t40) * s25 + t32)
     # + (((t17 * (s16 - s26) - t22 + s34) * s34 - (s16 * t20 + t34 - t3
     #5) * s16 + t36 + 10.0_dp * s25 * s26) * s34 - (t32 + t13) * s25 - 
     #(-t4 * t30 * t17 + 6.0_dp * t31) * s26 + 20.0_dp * t10) * s34 + t4
     #2) * s56 - t52 - t37 + t53 * (t3 - t2) + t41 * t6
      t20 = (-s34 + t12 + s26) * s26 + (t26 + s16 + s25 - s34 + s56) * s
     #56 + t50
      t29 = s16 * t7
      t55 = 6.0_dp * t29
      t33 = (t44 + t5) * s25 + t33 + t55
      t56 = s25 + s26
      t57 = t8 + t13
      t58 = s16 * t57
      t59 = t54 * s25
      t11 = t1 * (t58 + t39 - t47) - s34 * (-(t22 + t26 + t21) * s34 + (
     #t27 + t2) * s26 - t11 + 6.0_dp * t8) - t59
      t60 = s25 * t8 - t39
      t61 = t60 * s16 + t47 * (t5 + t56)
      t62 = t7 * ((t27 - t2) * s16 + (s25 + t5) * s26 + t13)
      t63 = (-t1 * s16 * (s25 + s34) - (s25 - s34) * s34 + t13 + t8) * s
     #56
      t64 = t57 * t8
      t65 = s25 + t24
      t66 = 10.0_dp * t8
      t67 = s25 + t27
      t68 = t67 * t18
      t69 = t19 * s16
      t70 = 20.0_dp * s16
      t71 = 13.0_dp * s16
      t72 = 15.0_dp * s26
      t73 = 8.0_dp * s26
      t16 = ((t1 * t58 - t59 - (s16 * t65 + (-t67 + s34) * s34) * s34) *
     # s56 + t1 * t46 + 6.0_dp * (t58 + t39) * s16 - (s16 * (t19 * s25 *
     # t56 + (-t23 + t72 + t70) * s16) + ((12.0_dp * s16 + t44 + t48) * 
     #s34 - s16 * (24.0_dp * s16 + 12.0_dp * s26) - (s16 * t16 + t2 + t4
     #4) * s25) * s34) * s34 - 12.0_dp * t10 * s25) * s56 - 12.0_dp * t8
     # * t60 + 6.0_dp * (t41 + t64) * s16 - ((((-10.0_dp * s16 - t22 - t
     #35 + s34) * s34 + (30.0_dp * s16 + 17.0_dp * s25 + 24.0_dp * s26) 
     #* s16 + t36 + t73 * s25) * s34 - s16 * ((t28 + t71) * s25 + 40.0_d
     #p * t8) - ((22.0_dp * s16 + t48) * s25 + 48.0_dp * t8) * s26 - t14
     # - t68 * t19) * s34 + t29 * (s16 * (25.0_dp * s16 + t72) + (t44 + 
     #t34 - t71) * s25)) * s34
      t34 = s16 - s25 - s34
      t36 = t8 + t45
      t56 = (-s25 + t69 + s26) * s26
      t47 = t10 - t47
      t58 = 8.0_dp * t10
      t59 = 15.0_dp * s16
      t60 = s25 * t12
      t71 = s16 * s34
      t72 = s25 * s34
      t12 = t1 * t50 * (t9 + t46) + s25 * (((-t71 * (8.0_dp * s16 + 5.0_
     #dp * s26 + t48) + s16 * (t25 + 12.0_dp * t8) + s26 * (s16 * (t28 +
     # t59) + s26 * t67)) * s34 - t7 * (-t1 * t50 * t4 + s26 * (s25 * t3
     #0 + t31 + t32) + t58)) * s34 + (((-t50 * t1 + t13 - t71 - t72 + t8
     #) * s56 + s16 * ((-t3 + t2) * s25 + t54) + t49 + (s34 * t67 - t1 *
     # t60 - t43 - 8.0_dp * t8) * s34 + t14) * s56 + 5.0_dp * s16 * t47 
     #+ t1 * (t39 * t67 + t38) + t42 + (((t73 + t28 + t59) * s16 + t26 *
     # s25) * s34 - t17 * ((t54 + t60) * s26 + t15) - t8 * (-s25 + t59) 
     #- t19 * t12 * t18) * s34 + t13 * t8 - t58 * s25) * s56 + t37 + t52
     # + t53 * (t24 + t2) + t41 * t67) - t19 * t9 * t13
      t25 = 0.1e1_dp / t51
      t20 = 0.1e1_dp / t20 ** 2.0_dp
      t27 = (s16 + s26 - s34 + s56) * (t1 * t36 - (t27 + t44 + t2) * s34
     # + (t19 * (s16 - s34) - s25 + t26 + s56) * s56 - t50 + t56) * MB11
     #11D2(0)
      t11 = t25 * t20 * ((-t1 * t61 * s16 + s16 * (((t45 + t33) * s34 - 
     #t62) * s34 + (t63 + t11) * s56 + t8 * (t8 + t13) + t41)) * I300s16
     #(1) + (t1 * t61 - ((t45 + t33) * s34 - t62) * s34 - (t63 + t11) * 
     #s56 - t41 - t64) * MB1001(0))
      result = -t20 * (t25 * ((-(t1 * (t10 * t57 + t37) * s16 + s34 * ((
     #((-(s26 + t5) * s34 + (t22 + t26) * s26 + 10.0_dp * t29 + 5.0_dp *
     # t50) * s34 - t17 * (t18 * (s25 + t69) + t15) - t8 * (t40 + t70) -
     # s26 * (t18 + (17.0_dp * s16 + t22) * s25 + 30.0_dp * t8)) * s34 +
     # t7 * (s16 * ((t23 + t69) * s25 + 20.0_dp * t8) + t14 + t68 + t26 
     #* ((s25 + t24) * s25 + t66))) * s34 - s16 * t7 ** 2.0_dp * ((-t3 +
     # t22) * s25 + s26 * t65 + t66)) + t16 * s56 + 6.0_dp * t39 * t8 * 
     #t7 - t48 * s16 * t9) * MB1011(1) + t12 * MB1110(1) - t6) / s16 + t
     #72 * (-t1 * t47 - ((t48 + t44 + t21) * s34 - (t35 + t2) * s25 - t1
     #8 - t50 - t55) * s34 - (t34 * s56 + t1 * t34 * s26 - t17 * s25 * (
     #s16 + s34) + t19 * t36 + t13 - t21 * s34) * s56 - t15 - t56 * t4 +
     # t22 * t8) * MB1101(0)) - t27) / 2.0_dp + t11

           intHs16s25s26s34s56x1222D6eps0 = result
       end function intHs16s25s26s34s56x1222D6eps0

       function intHs16s25s26s34s56x1222D6eps1()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1222D6eps1

           complex(dp) :: result

      result = 0.0_dp

           intHs16s25s26s34s56x1222D6eps1 = result
       end function intHs16s25s26s34s56x1222D6eps1

       function intHs16s25s26s34s56x1231D6eps0()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1231D6eps0
           complex(dp) ::  t1,t10,t100,t101,t102,t103,t11,t12,t13,t14,t15,t16
           complex(dp) ::  t17,t18,t19,t2,t20,t21,t22,t23,t24,t25,t26,t27
           complex(dp) ::  t28,t29,t3,t30,t31,t32,t33,t34,t35,t36,t37,t38
           complex(dp) ::  t39,t4,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49
           complex(dp) ::  t5,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t6
           complex(dp) ::  t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t7,t70
           complex(dp) ::  t71,t72,t73,t74,t75,t76,t77,t78,t79,t8,t80,t81
           complex(dp) ::  t82,t83,t84,t85,t86,t87,t88,t89,t9,t90,t91,t92
           complex(dp) ::  t93,t94,t95,t96,t97,t98,t99

           complex(dp) :: result

      t1 = 9.0_dp
      t2 = s16 ** 2.0_dp
      t3 = t2 ** 2.0_dp
      t4 = s16 * t3
      t5 = s16 * t2
      t6 = 37.0_dp * s16
      t7 = 16.0_dp * s25
      t8 = t2 * t1
      t9 = 13.0_dp * t5
      t10 = 12.0_dp
      t11 = 8.0_dp * s25
      t12 = t10 * t5
      t13 = ((47.0_dp * s16 + t11) * s25 + 35.0_dp * t2) * s25
      t14 = t13 - t12
      t15 = 2.0_dp
      t16 = t15 * s16
      t17 = 19.0_dp
      t18 = t17 * s25
      t19 = 8.0_dp * s16
      t20 = 14.0_dp * s25
      t21 = t20 + t19
      t22 = s16 + s26
      t23 = s25 ** 2.0_dp
      t24 = s25 * t23
      t25 = -t2 + t23
      t26 = -22.0_dp
      t27 = 18.0_dp * s16
      t28 = 4.0_dp
      t29 = t28 * s16
      t30 = -31.0_dp * s25
      t31 = s26 ** 2.0_dp
      t32 = t31 ** 2.0_dp
      t33 = s26 * t31
      t34 = s16 * t1
      t35 = 7.0_dp * s25
      t36 = 16.0_dp * t2
      t37 = 3.0_dp
      t38 = t37 * s16
      t39 = t38 - s25
      t40 = t29 + s25
      t41 = t16 * t23
      t42 = t37 * s26
      t43 = t2 * s25
      t44 = t37 * s25
      t45 = 6.0_dp * s26
      t46 = 10.0_dp * s16
      t47 = 16.0_dp * s16
      t48 = -t29 + s26
      t49 = 42.0_dp * s16
      t50 = -93.0_dp * s25
      t51 = -29.0_dp
      t52 = 38.0_dp * t2
      t53 = s16 + s25
      t54 = s25 * t53
      t55 = 40.0_dp * s26
      t56 = t10 * s26
      t57 = 15.0_dp * s26
      t58 = (t16 + s25) * (t18 + s16)
      t59 = t58 * t37
      t60 = 60.0_dp
      t61 = 5.0_dp * s25
      t62 = t1 * s26
      t63 = s34 ** 2.0_dp
      t64 = t63 ** 2.0_dp
      t65 = s34 * t63
      t66 = 14.0_dp * s16
      t67 = 32.0_dp * s16
      t13 = (((-10.0_dp * s34 + t19 + t20 + t57) * s56 + (39.0_dp * s16 
     #+ t18) * s25 + (56.0_dp * s25 + 30.0_dp * s26 + t67) * s26 + (s34 
     #* t10 + t30 - t55 - t66) * s34 + t15 * t2) * s56 + (30.0_dp * t31 
     #+ t59) * s26 + (t28 * (t62 + t61) - 6.0_dp * s34) * t63 - t12 + t1
     #3 + ((s25 * t26 - 49.0_dp * s16) * s25 + (-s26 * t60 - t49 + t50) 
     #* s26 + 18.0_dp * t2) * s34 + 6.0_dp * t31 * t21) * s56 + s16 * ((
     #(t7 + t6) * s25 + t8) * s25 - t9) + t31 * ((4.0_dp * t21 + t57) * 
     #s26 + t59) + t15 * t14 * s26 + ((t55 * s25 + (-t56 + t46 - t44 + s
     #34) * s34 + t37 * t54 - 36.0_dp * t2 + 36.0_dp * t31) * s34 + s16 
     #* ((s25 * t51 - t34) * s25 + t52) - 40.0_dp * t33 + ((-t49 + t50) 
     #* s26 + (-98.0_dp * s16 - 44.0_dp * s25) * s25 + 36.0_dp * t2) * s
     #26 - t24) * s34
      t18 = t3 * s25
      t20 = s56 ** 2.0_dp
      t49 = t20 ** 2.0_dp
      t50 = s56 * t49
      t51 = s56 * t20
      t55 = s16 + s25 + s26 - s34 + s56
      t57 = t15 * s26
      t59 = s16 * s25
      t68 = (-s34 + t53 + s26) * s26 + (s16 + s25 - s34 + t57 + s56) * s
     #56 + t59
      t69 = t15 * s25
      t70 = -s16 + s34
      t71 = t70 * s16
      t72 = s25 + s26
      t73 = s34 * t22
      t74 = t2 + t63
      t75 = s16 + s26 - s34 + s56
      t76 = (-t44 - s16) * s25
      t77 = 10.0_dp * s26
      t78 = 6.0_dp * s16
      t79 = t28 * s25
      t80 = 5.0_dp * s26
      t81 = t28 * s26
      t82 = s34 * t15
      t83 = t5 - t65
      t84 = (-s25 * t1 - t38) * s25 + t36
      t85 = 5.0_dp * s16
      t86 = (-t85 - t44) * s25 + 14.0_dp * t2
      t87 = 7.0_dp * s16
      t88 = t87 - t69
      t89 = 8.0_dp * t2
      t90 = 20.0_dp * t2
      t91 = 20.0_dp * s16
      t92 = 42.0_dp * t2
      t93 = (t44 + t91) * s25
      t94 = 25.0_dp * s16
      t95 = t16 + s26
      t96 = 6.0_dp * s25
      t97 = 13.0_dp * s16
      t98 = 26.0_dp
      t60 = t60 * s16
      t99 = 8.0_dp * t39
      t100 = s16 * s26
      t101 = t1 * t23
      t102 = s16 * t23
      t103 = t22 ** 2.0_dp
      t54 = ((t15 * s25 * t70 + (-7.0_dp * s34 + t78 + t80) * s16 + t63)
     # * s56 + s16 * ((-t61 + t66) * s16 + (t77 + t99) * s26) - t37 * (t
     #102 + t65) + s34 * ((t81 + t91 - s25) * s34 + (-28.0_dp * s16 + t1
     #1) * s26 + t28 * t54 - 31.0_dp * t2)) * s56 + s16 * ((t36 - t101) 
     #* s16 + t31 * (t10 * t39 + t77)) + t37 * ((s26 * t86 - t43) * s16 
     #+ t64) + ((-(t62 + t79 + t94) * s34 + (t97 - s25) * s25 + (t60 + t
     #45 - t44) * s26 + 57.0_dp * t2) * s34 + ((11.0_dp * s16 + t69) * s
     #25 - 6.0_dp * t2) * s25 - 51.0_dp * t5 + s26 * (-6.0_dp * s26 * t8
     #8 + t10 * t54 - 93.0_dp * t2)) * s34
      t17 = t22 * (s16 * ((t27 + t79) * s25 + t90) + ((s16 * t17 - s25 +
     # s26) * s26 + (t66 - s25) * s25 + t52) * s26 + t24) * t63
      t52 = t103 * (((-t38 - t69) * s25 + t89) * s25 + 10.0_dp * t5 + s2
     #6 * (s26 * t88 - t23 * t28 + 17.0_dp * t2)) * s34
      t66 = s16 * (t32 * t39 + t4)
      t91 = t37 * t23
      t1 = t15 * t66 - t95 * s34 * t64 + (t33 * (t86 + t31) + t50) * s16
     # + s25 * t4 + (t54 * s56 - t1 * t5 * t25 + ((((t97 + t45 + t44 - s
     #34) * s34 - (50.0_dp * s16 + t11 + t62) * s26 - t92 - t93) * s34 +
     # s16 * ((t44 + t67) * s25 + 58.0_dp * t2) + t28 * t33 + ((t60 - t4
     #4) * s26 + (s16 * t98 - t69) * s25 + 114.0_dp * t2) * s26 + t24) *
     # s34 - t22 * (((-t46 - t79) * s25 + t36) * s25 + 37.0_dp * t5 + s2
     #6 * (t28 * t88 * s26 + (t97 - t96) * (t85 + t69)))) * s34 + t18 + 
     #t100 * (((t80 + t99) * s26 + t37 * t86) * s26 + t16 * t84)) * s56 
     #+ t64 * (t37 * s26 * t72 + (13.0_dp * s26 + t61 + t46) * s16) - t6
     #5 * (s16 * ((t79 + t47) * s25 + t90) + t33 * t37 + s26 * ((t79 + t
     #94) * s26 + t92 + t93)) + t5 * (-t1 * t25 + t59) * s26 + t2 * t84 
     #* t31 + t17 - t52 - t91 * t3
      t4 = s16 - s25
      t17 = t4 * s26
      t52 = t2 + t23
      t54 = s25 * t52
      t59 = (t2 * t28 + t91) * s16 + t54
      t60 = (t85 + t96) * s25 + t89
      t66 = t44 + t16
      t19 = t19 + t96
      t67 = t10 * t2
      t84 = s16 * t22
      t64 = (t5 * t53 + t64) * s16
      t65 = t29 * (t15 * t22 + s25) * t65
      t86 = t4 ** 2.0_dp
      t88 = s26 * t86
      t90 = 0.1e1_dp / s16
      t55 = 0.1e1_dp / t55 ** 2.0_dp
      t68 = 0.1e1_dp / t68 ** 2.0_dp
      t12 = s25 * (t15 * t64 + t37 * (-t2 * t23 * t4 + t31 * t59) + s16 
     #* ((t2 * t37 + 6.0_dp * t23) * s25 + 8.0_dp * t5) * s26 + ((t40 * 
     #t37 * t31 + 10.0_dp * t43 + t100 * (24.0_dp * s16 + t35) + t12 + t
     #41) * s34 - t22 * (t19 * t31 + 8.0_dp * t2 * t53 - t102 + t57 * (-
     #t76 + t89))) * s34 + (6.0_dp * s16 * t24 + 6.0_dp * s26 * t59 + t3
     #7 * (s25 * t5 + t31 * t60) + t33 * t28 * t66 + (((t87 + t45) * s25
     # + 24.0_dp * t84) * s34 - s16 * ((t46 + t61) * s25 + 24.0_dp * t2)
     # - (t28 * ((t44 + t29) * s25 + t67) + t45 * t19 / 2.0_dp) * s26) *
     # s34 + 8.0_dp * s16 * t83) * s56 + t20 * ((t101 + t67) * s16 + t37
     # * (s26 * t60 + t40 * t63 + t54) - t82 * ((t62 + t44 + t29) * s25 
     #+ t10 * t84) + 6.0_dp * t66 * t31) + t32 * t66 + t33 * t60 + t49 *
     # t66 + t51 * ((t56 + t85) * s25 + 8.0_dp * (s26 - s34 + s16) * s16
     # - 6.0_dp * s25 * (-s25 + s34)) - t65)
      t19 = t75 * (t15 * t83 + ((t29 - t69 + s26) * s26 + 5.0_dp * t2 + 
     #t76) * s26 + ((t80 + t78 + t79) * s34 - t15 * s25 * t72 - (t61 + t
     #78 + t77) * s16 - t28 * t31) * s34 + ((-t28 * t70 + s56 + t42 - t6
     #9) * s56 + (-t82 - t81 - s16) * s25 - 8.0_dp * s26 * t70 - t37 * (
     #t23 - t31) + 5.0_dp * t74 - t46 * s34) * s56 + t43 - t38 * t23) * 
     #MB1001(0)
      t1 = t68 * ((-t15 * s16 * t74 - t37 * ((s16 * (-t69 + s16) - s25 *
     # s26) * s26 - t43) + ((t37 * t72 - s34 + s56) * s56 + 6.0_dp * s25
     # * t22 + t37 * (t71 + t31) - t57 * s34) * s56 + t33 - t73 * t48) *
     # I300s16(1) + t55 * (t90 * (-t1 * MB1011(1) + t2 * (t23 * (t34 + t
     #11) - t28 * t5) + t37 * (s26 * t32 + t50) + (((s26 * t21 + t58) * 
     #s26 + t14) * s26 + s16 * (((t7 + t6) * s25 + t8) * s25 - t9)) * s2
     #6 + ((-t42 * t39 * t40 + 20.0_dp * s25 * t31 + (t48 * s34 + s16 * 
     #(t35 + t47) + (t46 - t45 - t44) * s26) * s34 + t10 * t33 - t41 - 1
     #5.0_dp * t43 - 24.0_dp * t5) * s34 + t22 * (s16 * ((-t35 + t34) * 
     #s25 + t36) - 10.0_dp * t33 + (-t27 * s25 + (t30 - t29) * s26 + t25
     # * t26) * s26 - t24)) * s34 + t13 * s56 - t18 + t12 * MB1110(1)) +
     # t19))
      t5 = 0.1e1_dp / t75
      result = t1 / 4.0_dp + t90 * (t23 * MB1110(0) * t55 + ((-t37 * (s1
     #6 * (t22 * t95 * t63 + t88 * t22) + ((-t43 * t15 + s16 * ((s25 - s
     #26 - t16 + s34) * s34 + t2 + t23) + t88) * s56 + t15 * s16 * (-t43
     # + t88) + s16 * (s16 * t52 + ((t38 + t57 - s34) * s34 - t22 * (t38
     # - t69 + s26)) * s34) + t31 * t86) * s56) - s16 * ((-t103 * (-t44 
     #+ t29 + s26) + t63 * (-t29 - t42 - s25 + s34)) * s34 + t2 * t52) -
     # t33 * t86 - t51 * (-t16 * s25 + t23 - t71) + t69 * t3) * MB1111D2
     #(0) - s25 * s34 * (t15 * ((-t73 + t17) * s16 + ((-s25 - s34 + s16)
     # * s16 + t17) * s56) + (s16 * t4 + t63) * s16 + t4 * (t31 + t20)) 
     #* MB1101(0)) * t68 * t5) / 2.0_dp

           intHs16s25s26s34s56x1231D6eps0 = result
       end function intHs16s25s26s34s56x1231D6eps0

       function intHs16s25s26s34s56x1231D6eps1()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1231D6eps1
           complex(dp) ::  t1,t2,t3

           complex(dp) :: result

      t1 = 3.0_dp
      t2 = s16 + s25 + s26 - s34 + s56
      t3 = 0.1e1_dp / 0.4e1_dp
      t2 = 0.1e1_dp / t2 ** 2.0_dp
      result = 0.1e1_dp / s16 * t2 * (t3 * (s25 * t1 + s16 + s26 - s34 +
     # s56) + s25 ** 2.0_dp * MB1110(1) / 2.0_dp)

           intHs16s25s26s34s56x1231D6eps1 = result
       end function intHs16s25s26s34s56x1231D6eps1

       function intHs16s25s26s34s56x1310D4eps0()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1310D4eps0
           complex(dp) ::  t1,t2,t3,t4,t5

           complex(dp) :: result

      t1 = s16 + s26 - s34 + s56
      t2 = s16 + s25 + s26 - s34 + s56
      t3 = 3.0_dp
      t4 = 6.0_dp
      t5 = 4.0_dp * s25
      t2 = 0.1e1_dp / t2 ** 2.0_dp
      result = t2 * (((s16 + s26 + s56) * MB1011(1) - MB1001(0)) * (t1 *
     # t3 + s25) + (t3 * (s16 ** 2.0_dp + s26 ** 2.0_dp + s34 ** 2.0_dp 
     #+ s56 ** 2.0_dp) + (s16 * t4 - t5) * s26 + (-t4 * (s16 + s26) + t5
     #) * s34 + (t4 * (s16 + s26 - s34) - t5) * s56 - s25 ** 2.0_dp - t5
     # * s16) * MB1110(1)) / 4.0_dp - t1 * t2 * (-t1 * MB1110(0) + 1.0_d
     #p) / 2.0_dp

           intHs16s25s26s34s56x1310D4eps0 = result
       end function intHs16s25s26s34s56x1310D4eps0

       function intHs16s25s26s34s56x1310D4eps1()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1310D4eps1
           complex(dp) ::  t1,t2

           complex(dp) :: result

      t1 = s16 + s26 - s34 + s56
      t2 = s16 + s25 + s26 - s34 + s56
      t2 = 0.1e1_dp / t2 ** 2.0_dp
      result = t2 * (-s25 / 4.0_dp - 0.3e1_dp / 0.4e1_dp * t1 + t1 ** 2.
     #0_dp * MB1110(1) / 2.0_dp)

           intHs16s25s26s34s56x1310D4eps1 = result
       end function intHs16s25s26s34s56x1310D4eps1

       function intHs16s25s26s34s56x1311D4eps0()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1311D4eps0
           complex(dp) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           complex(dp) ::  t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30
           complex(dp) ::  t31,t32,t33,t34,t35,t36,t37,t38,t39,t4,t40,t41
           complex(dp) ::  t42,t5,t6,t7,t8,t9

           complex(dp) :: result

      t1 = s16 - s25
      t2 = s16 + s25
      t3 = 2.0_dp
      t4 = t3 * s25
      t5 = t4 * s26
      t6 = s25 ** 2.0_dp
      t7 = s25 * t6
      t8 = s16 ** 2.0_dp
      t9 = s16 * t8
      t10 = s34 ** 2.0_dp
      t11 = t10 + t8 + t6
      t12 = t8 + t6
      t13 = t1 ** 2.0_dp
      t14 = t13 * s26
      t15 = t12 * s16
      t16 = t4 * t8
      t17 = -t16 + (-s16 * t2 + s26 * s34 - t5) * s34 + (-t3 * s25 * (s1
     #6 + s34) + t11) * s56 + t14 + t15
      t18 = s16 * s25
      t11 = -t3 * (s34 * t2 + t18) + t11
      t19 = s16 + s25 - s34
      t20 = t3 * s26
      t21 = (-s34 + t2 + s26) * s26 + (t20 + t19 + s56) * s56 + t18
      t22 = s26 ** 2.0_dp
      t23 = s16 * t1
      t24 = s26 * t2
      t25 = 4.0_dp
      t26 = t25 * s26
      t27 = s56 ** 2.0_dp
      t28 = t8 * s25
      t29 = s16 + s25 + s26 - s34 + s56
      t30 = s16 + s26 + s56
      t31 = t8 - t6
      t32 = t4 * s16
      t33 = s16 * t6
      t34 = 3.0_dp
      t35 = t34 * s16
      t36 = s16 + t4
      t37 = t3 * s16
      t38 = s25 + t37
      t39 = t3 * t23 + t6
      t40 = t2 * t22
      t41 = s16 * s34
      t21 = 0.1e1_dp / t21
      t42 = 0.1e1_dp / s16
      t11 = 0.1e1_dp / t11
      t29 = 0.1e1_dp / t29
      t5 = t30 * (-t4 * t9 + (t15 + t14) * s16 + ((-(s16 - s26) * s34 + 
     #t34 * t8 - t5 - s26 * (s16 + s26)) * s34 + t19 * t27 + t8 * (t4 - 
     #t35) + (t24 + t32 - t31) * s26 + t33) * s34 + (-t16 + ((-t3 * (s25
     # + s26) - s16 + s34) * s34 + t3 * (t24 + t18) - t31) * s34 + t15) 
     #* s56) * MB1011(1)
      t4 = -(t3 * t28 * t12 + (-t25 * t28 * t1 - s26 * (t23 * (s25 + t35
     #) + t40) - t3 * t39 * t22) * s34 - ((-(-s25 * s34 - t32 - t41 + t6
     # + t8) * s56 - t3 * t9 + t34 * (t28 - t14) - (-t34 * s26 * t2 + s3
     #4 * t38 - t25 * t23 - t3 * t6) * s34 - t7) * s56 - t34 * (s16 * t7
     # + t13 * t22) - ((-t41 + (s16 * t25 + t4) * s26 + t23 * t34) * s34
     # - t34 * (t40 + t9) - t26 * t39 + t18 * (s25 + t37)) * s34 - t2 * 
     #t9 - t20 * t13 * t38 + 5.0_dp * t8 * t6) * s56 + t10 * ((s26 * t38
     # + t23 * t34) * s26 + t16) - t14 * (-s16 * (s25 * t34 + s16) + (-t
     #38 - s26) * s26) - s16 * s26 * s34 * t10 - t25 * t9 * t6) * MB1110
     #(1) + t5
      t5 = -0.1e1_dp / 0.2e1_dp
      result = t5 * (t21 * (t11 * (t17 * I300s16(1) + (-(-t3 * (-s34 * t
     #27 + t28) + ((s16 - t20) * s34 + t3 * (-t23 + t22) + t24) * s34 + 
     #(-t3 * (t18 + t10) + (t26 + t2) * s34 + t6 + t8) * s56 + t14 + t15
     #) * MB1001(0) + t4 * t42) * t29 - s34 * ((-s26 * t1 - t8) * s25 + 
     #s26 * t10 - (s26 * t36 + t18) * s34 - (s25 * t1 + (t36 - s34) * s3
     #4) * s56 + t33) * MB1101(0) * t42) + t30 ** 2.0_dp * MB1111D2(0) *
     # t42) - (s16 + s26 - s34 + s56) * MB1110(0) * t42 * t29) - t17 * t
     #42 * t11 * t21

           intHs16s25s26s34s56x1311D4eps0 = result
       end function intHs16s25s26s34s56x1311D4eps0

       function intHs16s25s26s34s56x1311D4eps1()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1311D4eps1
           complex(dp) ::  t1

           complex(dp) :: result

      t1 = s16 + s25 + s26 - s34 + s56
      t1 = 0.1e1_dp / t1
      result = 0.1e1_dp / s16 * t1 * ((s16 + s26 - s34 + s56) * MB1110(1
     #) - 1.0_dp) / 2.0_dp

           intHs16s25s26s34s56x1311D4eps1 = result
       end function intHs16s25s26s34s56x1311D4eps1

       function intHs16s25s26s34s56x1312D6eps0()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1312D6eps0
           complex(dp) ::  t1,t10,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109
           complex(dp) ::  t11,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t12
           complex(dp) ::  t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t13,t130
           complex(dp) ::  t131,t132,t133,t134,t135,t136,t137,t138,t139,t14,t140,t141
           complex(dp) ::  t142,t143,t144,t145,t146,t147,t148,t149,t15,t150,t151,t152
           complex(dp) ::  t16,t17,t18,t19,t2,t20,t21,t22,t23,t24,t25,t26
           complex(dp) ::  t27,t28,t29,t3,t30,t31,t32,t33,t34,t35,t36,t37
           complex(dp) ::  t38,t39,t4,t40,t41,t42,t43,t44,t45,t46,t47,t48
           complex(dp) ::  t49,t5,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59
           complex(dp) ::  t6,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t7
           complex(dp) ::  t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t8,t80
           complex(dp) ::  t81,t82,t83,t84,t85,t86,t87,t88,t89,t9,t90,t91
           complex(dp) ::  t92,t93,t94,t95,t96,t97,t98,t99

           complex(dp) :: result

      t1 = s16 - s25
      t2 = 9.0_dp
      t3 = 2.0_dp
      t4 = s16 ** 2.0_dp
      t5 = t4 ** 2.0_dp
      t6 = s16 * t4
      t7 = t3 * s25
      t8 = 5.0_dp * s16
      t9 = t2 * t4
      t10 = s25 * t1
      t11 = 6.0_dp * t4
      t12 = -t10 + t11
      t13 = s16 + s25
      t14 = s16 + s26
      t15 = 4.0_dp
      t16 = t15 * s25
      t17 = 12.0_dp * t4
      t18 = 19.0_dp * t6
      t19 = 12.0_dp * s16
      t20 = 3.0_dp * s25
      t21 = (t20 + t19) * s25 + 5.0_dp * t4
      t22 = 20.0_dp
      t23 = 14.0_dp * t4
      t24 = t22 * t6
      t25 = t3 * s16
      t26 = 23.0_dp * t4
      t27 = -66.0_dp * t6
      t28 = 6.0_dp * s25
      t29 = 14.0_dp * s16
      t30 = t22 * t4
      t31 = ((t28 + t29) * s25 - t30) * s25
      t32 = 36.0_dp * t6
      t33 = t31 - t32
      t34 = 8.0_dp * s16
      t35 = 7.0_dp * t4
      t36 = (t20 + t34) * s25 - t35
      t29 = (-t29 - t7) * s25
      t37 = 40.0_dp * t6
      t38 = s25 + t25
      t39 = 27.0_dp * s16
      t40 = 7.0_dp * s25
      t41 = s25 ** 2.0_dp
      t42 = t41 ** 2.0_dp
      t43 = s25 * t42
      t44 = s25 * t41
      t45 = 24.0_dp * t4
      t46 = t15 * t41
      t47 = t45 - t46
      t48 = 3.0_dp * s16
      t49 = t48 - s25
      t50 = 10.0_dp
      t51 = t50 * s25
      t52 = 18.0_dp * s16
      t53 = 7.0_dp * s16
      t54 = t20 + t53
      t55 = 6.0_dp * s16
      t56 = s25 - t55
      t57 = s26 ** 2.0_dp
      t58 = s26 * t57
      t59 = 11.0_dp * s25
      t60 = t22 * s16
      t61 = 3.0_dp * s26
      t62 = t15 * s16
      t63 = t4 * t13
      t64 = 17.0_dp * s16
      t65 = 3.0_dp * t21 * s26
      t66 = 3.0_dp * s26 * t36
      t67 = s25 * t2
      t68 = 12.0_dp * s26
      t69 = 21.0_dp * s16
      t70 = t1 ** 2.0_dp
      t70 = t1 * t70
      t71 = t70 * s26
      t72 = s34 ** 2.0_dp
      t73 = t72 ** 2.0_dp
      t74 = s34 * t73
      t75 = s34 * t72
      t76 = t13 * t57
      t77 = s16 * t74
      t78 = 30.0_dp
      t79 = 32.0_dp * t6
      t80 = t2 * s16
      t81 = 22.0_dp * t4
      t82 = t61 * t70
      t83 = s16 * s25
      t84 = t3 * t83 * (-t41 + t4)
      t31 = (((s34 * t49 + t36) * s34 + t1 * t21) * s34 + t42 - t5 + t84
     #) * s56 + s16 * ((((-t16 + t19) * s25 - t81) * s25 + t18) * s25 - 
     #6.0_dp * t5) + (((-t61 * s25 + t80 * s26 + s34 * t56 + t45 - t46) 
     #* s34 + t31 - t32 + t66) * s34 + t1 * (t41 * (t16 + t34) + 24.0_dp
     # * t63 + t65)) * s34 + t43 - t82 * t13
      t31 = t31 * s56 - t3 * (-s16 * t43 + t71 * t12) + t4 * ((((s16 * t
     #78 - t59) * s25 - 44.0_dp * t4) * s25 + t79) * s25 - t2 * t5) + ((
     #(-((t68 + t69 + t67) * s16 - t7 * s26) * s34 + s16 * ((41.0_dp * s
     #16 + t40) * s25 + 54.0_dp * t4) + s26 * (3.0_dp * s26 * t49 + t3 *
     # t47)) * s34 + s16 * (((t20 + t25) * s25 - t26) * s25 + t27) + s26
     # * (t3 * t33 + t66)) * s34 + t1 * (s16 * (((t28 + t64) * s25 - t3 
     #* t4) * s25 + 39.0_dp * t6) + (8.0_dp * t38 * t41 + 48.0_dp * t63 
     #+ t65) * s26)) * s34 - 3.0_dp * t76 * t70 + 3.0_dp * t77
      t32 = t70 * t58
      t34 = t71 * s16
      t63 = t14 * s34
      t65 = t63 * t1
      t66 = 16.0_dp * s16
      t85 = 15.0_dp * t6
      t86 = -t3 * (s34 * t13 + t83) + t72 + t41 + t4
      t87 = s16 + s25 - s34
      t88 = t3 * s26
      t89 = (-s34 + t13 + s26) * s26 + (t88 + t87 + s56) * s56 + t83
      t90 = t8 - t7
      t91 = t62 - s25
      t92 = t62 + t28
      t93 = t92 * s26
      t94 = t83 * t1
      t95 = 11.0_dp * s16
      t96 = 5.0_dp * s25
      t97 = (-t96 + t95) * s25 + t45
      t98 = -t20 + t8
      t99 = 3.0_dp * t4
      t100 = 24.0_dp * s16
      t101 = t3 * t58
      t102 = s25 * t14
      t103 = t15 * s26
      t104 = 13.0_dp * s16
      t105 = 0.3e1_dp / 0.2e1_dp * t92 * s26
      t106 = 13.0_dp * t4
      t107 = s25 + s26
      t108 = 8.0_dp * s26
      t109 = s16 * t42
      t110 = t71 * t91
      t111 = 17.0_dp * t4
      t112 = t4 * (((-t59 + t69) * s25 - t111) * s25 + 5.0_dp * t6)
      t113 = 3.0_dp * t70 * t57
      t114 = 6.0_dp * s26
      t115 = 15.0_dp * s16
      t116 = (-t40 + t115) * s25 - t106
      t92 = t92 / 2.0_dp
      t117 = 3.0_dp * t94
      t46 = ((-t3 * s34 * (t1 * t92 + t72) + t72 * t98 - t44 + t6 - t117
     #) * s56 + t15 * (t73 + t5) + t82 + ((-(t20 + t114 + t66) * s34 + s
     #16 * (15.0_dp * s26 + t100) + (-s26 * t2 + t95 - t96) * s25) * s34
     # - t1 * (s16 * (t68 + t66) + (18.0_dp * s26 + t20 + t95) * s25)) *
     # s34 + t42 + t83 * t116) * s56 - t3 * ((t1 * ((11.0_dp * t4 + t46)
     # * s16 + s26 * ((t20 + t95) * s25 + 16.0_dp * t4 + t105)) + t72 * 
     #(16.0_dp * s16 * t14 + 3.0_dp * s26 * t107 + 3.0_dp * t41 + t72 + 
     #12.0_dp * t83)) * s34 - t109 - t110) + t72 * (t72 * (t28 + t108 + 
     #t104) + (t3 * t41 + t106) * s25 + 38.0_dp * t6 + s16 * t41 + (3.0_
     #dp * s26 * t98 + t3 * t97) * s26) + t112 + t113
      t118 = t4 * t42
      t91 = t70 * t91 * t57
      t119 = ((-t96 + t80) * s25 - t35) * s25
      t120 = s16 + s26 + s56
      t121 = t6 * t15
      t122 = ((t80 + t7) * s25 + t26) * s25
      t123 = t79 + t122
      t124 = s25 + t55
      t125 = t13 * t124
      t126 = t78 * t6
      t127 = 28.0_dp * s16
      t128 = 46.0_dp * t4
      t129 = 28.0_dp * t4
      t130 = (t28 + t39) * s25 + t129
      t69 = (t69 + t67) * s25
      t131 = s16 * t50
      t132 = 18.0_dp * t4
      t133 = (t96 + t53) * s25
      t134 = t4 * t15
      t135 = 12.0_dp * s25
      t136 = s16 * t1
      t137 = -t44 + t6
      t138 = t137 * s16
      t139 = t10 * t4
      t79 = ((-t3 * s34 * (t136 * t92 + t38 * t72) + t72 * (t72 + t125) 
     #+ t138 - 3.0_dp * t139) * s56 - t3 * t74 + t4 * (s25 * t116 + t121
     #) + ((-((27.0_dp * s25 + t68 + t127) * s16 + 6.0_dp * s25 * t107) 
     #* s34 + t125 * t61 + t122 + t79) * s34 - t136 * (18.0_dp * s25 * s
     #26 + t19 * s26 + t133 + 18.0_dp * t4)) * s34 + t109 + 3.0_dp * t73
     # * (t62 + s26 + t7) + 3.0_dp * t34) * s56 + s16 * (t113 + t112) - 
     #t3 * (-s16 * (t110 + t109) + (t72 * (t72 * (t107 * t3 + t8) + (-t2
     #9 + t26) * s25 + (t61 * t38 + t130) * s26 + t126) + t136 * (s16 * 
     #((-t62 + t28) * s25 + t106) + s26 * (t132 + t133 + t105))) * s34) 
     #+ (t72 * ((29.0_dp * s16 + t28) * s25 + (t61 + t100 + t135) * s26 
     #+ 35.0_dp * t4 + t72) + t3 * t123 * s26 + s16 * ((t69 + t134) * s2
     #5 + 55.0_dp * t6) + t42 + 3.0_dp * t76 * t124) * t72
      t92 = t71 * t4
      t105 = t4 * s25
      t106 = t41 + t4
      t107 = t57 + t4
      t110 = t53 + t7
      t112 = t62 + s25
      t113 = 19.0_dp * s16
      t55 = t16 + t55
      t116 = 31.0_dp
      t122 = 23.0_dp * s16
      t133 = 8.0_dp * s25
      t136 = 15.0_dp * s25
      t45 = (-t95 - t136) * s25 + t45
      t140 = t48 - t96
      t141 = t140 * s26
      t142 = t13 * ((-t133 - t122) * s25 + t128)
      t143 = s16 * (((-t53 - t51) * s25 + t116 * t4) * s25 + 44.0_dp * t
     #6)
      t144 = 13.0_dp * s25
      t145 = -t62 + t144
      t146 = (-t28 + t53) * t38
      t147 = s16 * (s25 + t53)
      t148 = (s16 - s26) * s26
      t149 = (69.0_dp * s16 + 21.0_dp * s25) * s25 + 60.0_dp * t4
      t150 = t55 / 2.0_dp
      t151 = 3.0_dp * t45
      t152 = 36.0_dp * t4
      t35 = ((-t3 * s34 * (t1 * t150 - t72) - t117 + t140 * t72 - t44 + 
     #t6) * s56 + t15 * (-t73 + t71 + t5) - ((-(-t62 + t108 + t144) * s3
     #4 - s16 * (t68 + t100) - (-s26 * t22 - t136 - t95) * s25) * s34 + 
     #t1 * (s16 * (24.0_dp * s26 + t60) + (16.0_dp * s26 + t40 + t122) *
     # s25)) * s34 - t13 * t44 + t105 * (t67 - t95)) * s56 - t3 * (-t74 
     #+ t109) + t6 * ((t136 - t113) * s25 + t35) + ((((t48 - t133 - t68)
     # * s34 + (39.0_dp * s26 + t131) * s25 - t129 - 12.0_dp * t148 + 12
     #.0_dp * t41) * s34 + ((-s16 * t116 - t133) * s25 + t26) * s25 + (t
     #151 + 6.0_dp * t141) * s26 + 46.0_dp * t6) * s34 - t1 * (((22.0_dp
     # * s16 + t7) * s25 + t152) * s25 + t126 + s26 * (6.0_dp * t55 * s2
     #6 + t149))) * s34 - t4 * t44 + 6.0_dp * t70 * t57 + t82 * t112
      t1 = t35 * s56 + t3 * (t73 * ((t48 - t114 - t133) * s26 + t147) + 
     #t34 * t110) + t5 * ((-t64 + t136) * s25 + t11) + (((t72 * (t103 - 
     #t25) + s16 * ((t28 - t39) * s25 - t152) + ((3.0_dp * t145 + t108) 
     #* s26 - t146 * t15) * s26) * s34 + ((t141 * t15 + t151) * s26 + t1
     #42 * t3) * s26 + t143) * s34 - t1 * (t15 * (s16 * t44 + (((s25 + t
     #95) * s25 + t132) * s25 + t85) * s26) + t4 * ((25.0_dp * s25 + t11
     #5) * s25 + 26.0_dp * t4) + t57 * (t108 * t150 + t149))) * s34 - t1
     #18 + 3.0_dp * t70 * t112 * t57 - 3.0_dp * t44 * t6 + t32 * t15
      t35 = s16 + s26 - s34 + s56
      t64 = 0.1e1_dp / t86 ** 2.0_dp
      t67 = 0.1e1_dp / t89 ** 2.0_dp
      t42 = t120 * (t5 * (t3 * t6 + t119) + (s26 + t25) * t72 * t73 + t7
     #9 * s56 + t42 * t6 + t72 * (t4 * (7.0_dp * t44 + t126) + ((t125 * 
     #s26 + t123) * s26 + (((s25 + t80) * s25 + 21.0_dp * t4) * s25 + t1
     #21) * s25 + 55.0_dp * t5) * s26 + t109 - 14.0_dp * t10 * t6) + t73
     # * (s16 * (t4 * t78 + t69) + s26 * (t54 * (t8 + t7) + t57) + 6.0_d
     #p * t38 * t57) - t74 * (s16 * (t40 + t19) + (t88 + t16 + t131) * s
     #26) - t75 * (s16 * (((t96 + t66) * s25 + t23) * s25 + t37) + s26 *
     # (t130 * s26 + ((t16 + t127) * s25 + t128) * s25 + 60.0_dp * t6) +
     # t101 * t38) + t32 * s16 + t91 * s16 + t92 * t90 - t65 * s16 * (s1
     #6 * ((t40 - t80) * s25 + t17) + s26 * ((s16 + t96) * s25 + t23 + t
     #93))) * MB1011(1)
      t1 = s25 * (-t3 * (-t137 * t5 + t74 * (t148 + t4)) + t57 ** 2.0_dp
     # * t70 + t1 * s56 + t72 * (t4 * (t41 * (t131 - t7) + t24) + (((t45
     # + t141) * s26 + t142) * s26 + t143) * s26) + t73 * (((-t103 + t48
     # - t133) * s26 + t147 * t3) * s26 + t4 * (t28 + t131)) + t75 * (t1
     #45 * t58 - t22 * t5 - t3 * t57 * (t146 - t57) - t134 * t41 - t135 
     #* t6 - t61 * s16 * ((t80 - t7) * s25 + t17)) + s16 * t70 * t110 * 
     #t57 + t92 * t124 + t32 * t112 - t65 * (t50 * t4 * t106 - t84 + ((t
     #55 * s26 + (t40 + t113) * s25 + t23) * s26 + ((t115 + t7) * s25 + 
     #t111) * s25 + 16.0_dp * t6) * s26) - 6.0_dp * t10 * s16 * t5) * MB
     #1110(1)
      t1 = -t1 - t6 * ((((-t28 + t66) * s25 - t81) * s25 + t85) * s25 - 
     #t15 * t5) - t31 * s56 - t4 * t43 - (t4 * ((t29 + t17) * s25 - t37)
     # + ((s26 * t36 + t33) * s26 + s16 * (((t20 + t25) * s25 - t26) * s
     #25 + t27)) * s26) * t72 - (t4 * ((t52 + t51) * s25 + 40.0_dp * t4)
     # + ((s26 * t49 + t47) * s26 + s16 * t38 * (t40 + t39)) * s26) * t7
     #5 - (((-t59 - t60) * s16 - t61 * t54) * s16 + t56 * t57) * t73 + t
     #32 * t13 + t70 * t12 * t57 - t77 * (t61 + t62) + t34 * ((-t8 + t7)
     # * s25 + t9) - t65 * (s16 * (((t19 + t7) * s25 - t23) * s25 + t24)
     # + s26 * (t21 * s26 + ((t16 + t8) * s25 + t17) * s25 + t18)) + t42
      t2 = t67 * (t120 * t35 ** 2.0_dp * MB1111D2(0) - (-s25 * s34 * (t7
     #0 * t14 ** 2.0_dp + (((-t8 * s26 - t107 * t15 - 3.0_dp * t102 + t6
     #3) * s34 + 6.0_dp * t6 + s26 * ((t40 + t80) * s26 + (t20 + t104) *
     # s25 + t9) + 3.0_dp * t83 * t13 + 3.0_dp * t58) * s34 - t14 * (t15
     # * s16 * t106 - 3.0_dp * t105 + 3.0_dp * t76 + (t48 + s25) * (s16 
     #+ t7) * s26 + t44)) * s34 + (-(t15 * t75 + t117 + (-(t14 * t2 + t4
     #0) * s34 + t2 * t13 * s26 + (t131 + t7) * s25 + t11) * s34 + t44 -
     # t6) * s56 + t3 * (t138 + t71) - 6.0_dp * t139 - (((t20 + t8 + t10
     #8 - s34) * s34 - (14.0_dp * s26 + t20 + t104) * s25 - t107 * t2 - 
     #t52 * s26) * s34 + t124 * t41 + t15 * (((s25 + t8) * s25 + t99) * 
     #s26 + t105) + 7.0_dp * t6 + t76 * t2) * s34) * s56) + t20 * t87 * 
     #t72 * s56 ** 3.0_dp) * MB1101(0) * t64)
      result = -t64 * t67 * ((s16 * I300s16(1) - MB1001(0)) * (-t3 * (t1
     #4 * t74 - t4 * t5) + t46 * s56 + t72 * (s16 * (((t53 + t7) * s25 -
     # t99) * s25 + t24) + ((s26 * t98 + t97) * s26 + t38 * ((-t48 + t7)
     # * s25 + 19.0_dp * t4)) * s26) + t73 * ((t103 + t104) * s26 + t4 *
     # t50 + 6.0_dp * t102) - t75 * (s16 * ((t28 + t95) * s25 + t30) + s
     #26 * ((t20 + t66) * s26 + (t28 + t100) * s25 + 32.0_dp * t4) + t10
     #1) + t32 + t118 + t91 + t34 * t90 - t65 * (t50 * t6 + s26 * ((t20 
     #+ t8) * s25 + t17 + t93) - 5.0_dp * t94) + t119 * t6) + t1 / s16) 
     #/ 4.0_dp - t2 / 2.0_dp

           intHs16s25s26s34s56x1312D6eps0 = result
       end function intHs16s25s26s34s56x1312D6eps0

       function intHs16s25s26s34s56x1312D6eps1()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1312D6eps1

           complex(dp) :: result

      result = 0.0_dp

           intHs16s25s26s34s56x1312D6eps1 = result
       end function intHs16s25s26s34s56x1312D6eps1

       function intHs16s25s26s34s56x1321D6eps0()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1321D6eps0
           complex(dp) ::  t1,t10,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109
           complex(dp) ::  t11,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t12
           complex(dp) ::  t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t13,t130
           complex(dp) ::  t131,t132,t133,t134,t135,t136,t137,t138,t139,t14,t140,t141
           complex(dp) ::  t142,t143,t144,t145,t146,t147,t148,t149,t15,t150,t151,t152
           complex(dp) ::  t153,t154,t155,t156,t157,t158,t159,t16,t160,t161,t162,t163
           complex(dp) ::  t164,t165,t166,t167,t168,t169,t17,t170,t171,t172,t173,t174
           complex(dp) ::  t175,t176,t177,t178,t179,t18,t180,t19,t2,t20,t21,t22
           complex(dp) ::  t23,t24,t25,t26,t27,t28,t29,t3,t30,t31,t32,t33
           complex(dp) ::  t34,t35,t36,t37,t38,t39,t4,t40,t41,t42,t43,t44
           complex(dp) ::  t45,t46,t47,t48,t49,t5,t50,t51,t52,t53,t54,t55
           complex(dp) ::  t56,t57,t58,t59,t6,t60,t61,t62,t63,t64,t65,t66
           complex(dp) ::  t67,t68,t69,t7,t70,t71,t72,t73,t74,t75,t76,t77
           complex(dp) ::  t78,t79,t8,t80,t81,t82,t83,t84,t85,t86,t87,t88
           complex(dp) ::  t89,t9,t90,t91,t92,t93,t94,t95,t96,t97,t98,t99

           complex(dp) :: result

      t1 = s16 - s25
      t2 = 2.0_dp
      t3 = t2 * s16
      t4 = 3.0_dp
      t5 = t4 * s16
      t6 = -s25 + t5
      t7 = s16 + s26
      t8 = 4.0_dp
      t9 = s16 ** 2.0_dp
      t10 = t9 ** 2.0_dp
      t11 = s16 * t9
      t12 = t9 * t10
      t13 = 7.0_dp * s16
      t14 = t8 * s25
      t15 = s16 + s25
      t16 = t2 * s25
      t17 = 9.0_dp * s16
      t18 = 15.0_dp * t9
      t19 = t15 * s26
      t20 = 5.0_dp
      t21 = t20 * s16
      t22 = 10.0_dp * s16
      t23 = t22 + s25
      t24 = s16 * s25
      t25 = s26 ** 2.0_dp
      t26 = t25 ** 2.0_dp
      t27 = s26 * t26
      t28 = s26 * t25
      t29 = t2 * s26
      t30 = 15.0_dp * t11
      t31 = t8 * s16
      t32 = t20 * t9
      t33 = 11.0_dp * s16
      t34 = (-t14 + t33) * s25
      t35 = 17.0_dp * t9
      t36 = 6.0_dp
      t37 = 14.0_dp * t11
      t38 = t36 * t25
      t39 = t4 * s26
      t40 = t1 ** 2.0_dp
      t41 = s25 ** 2.0_dp
      t42 = t41 ** 2.0_dp
      t43 = s25 * t41
      t44 = s16 * t43
      t45 = t40 * s26
      t46 = s34 ** 2.0_dp
      t47 = t46 ** 2.0_dp
      t48 = s34 * t47
      t49 = s34 * t46
      t50 = 12.0_dp * t9
      t51 = t36 * s26
      t52 = 21.0_dp
      t53 = 15.0_dp * s16
      t54 = t4 * s25
      t55 = t9 * t52
      t56 = 9.0_dp * t11
      t57 = t45 * t36
      t58 = -t2 * (s34 * t15 + t24) + t41 + t46 + t9
      t59 = t9 * t43
      t60 = t40 * t25
      t61 = s56 ** 2.0_dp
      t62 = t40 * t28
      t63 = 13.0_dp * s16
      t64 = 11.0_dp * s25
      t65 = (-s34 + t15 + s26) * s26 + (s16 + s25 - s34 + t29 + s56) * s
     #56 + t24
      t66 = s16 * t15
      t67 = t2 * t41
      t68 = t36 * s16
      t69 = t13 + s25
      t70 = 16.0_dp * t11
      t71 = 14.0_dp * t9
      t72 = t20 * s25
      t73 = t72 + t5
      t74 = t73 * s26
      t75 = t11 * t4
      t76 = 34.0_dp * t11
      t77 = 12.0_dp * s25
      t78 = (19.0_dp * s16 + t77) * s25
      t79 = 36.0_dp * t9
      t80 = t79 + t78
      t81 = 10.0_dp * t9
      t82 = 13.0_dp * t15
      t83 = 14.0_dp * s16
      t84 = (t83 + t54) * s25
      t85 = 26.0_dp * t9
      t86 = 18.0_dp * s16
      t87 = s16 * t41
      t88 = 9.0_dp * s25
      t89 = (s25 + s26) * s26 + t24
      t90 = 9.0_dp * s26
      t91 = 9.0_dp * t9
      t92 = t2 * t9
      t93 = 10.0_dp * t11
      t94 = 12.0_dp * s16
      t95 = (33.0_dp * s25 + t94) * s25 + 51.0_dp * t9
      t96 = 8.0_dp * s25
      t97 = 15.0_dp * s25
      t98 = t8 * s26
      t99 = 36.0_dp * s16
      t100 = 22.0_dp * s25
      t101 = t10 * t1
      t102 = s16 * t42
      t103 = t9 * t41
      t104 = 8.0_dp * t28
      t105 = 7.0_dp * s25
      t106 = t36 * t41
      t107 = 20.0_dp * s26
      t108 = 17.0_dp * s16
      t109 = 12.0_dp * s26
      t110 = 7.0_dp * t9
      t111 = t41 + t9
      t73 = (((-t72 + s34 - t5) * s34 + t111 * t2 - t24 * t8) * s56 - s1
     #6 * ((-t72 + t63) * s25 - t110) - t4 * t49 + ((t82 + t98) * s34 - 
     #s16 * (t109 + t108) - (t31 + t64 + t107) * s25) * s34 + t43 + 8.0_
     #dp * t45) * s56 - s16 * ((t71 - t106) * s25 - t56) + t4 * ((-s16 *
     # ((t105 - t3) * s25 + t81) - (t2 * t73 * s26 + t95 / 3.0_dp) * s26
     # - t43 + t49) * s34 + t45 * t69) - t46 * ((t90 + t64 + t86) * s34 
     #- 39.0_dp * t19 - t38 - t78 - t79) - t42 + 12.0_dp * t60
      t73 = t73 * s56 + t101 * t20 - t2 * (-t45 * ((t31 - s25) * s25 + t
     #91) + t102) + t4 * (t47 * (s25 + t5 + t29) + t60 * t69) - ((((t90 
     #+ t99 + t100) * s26 + t84 + t85 + t46) * s34 - t2 * t80 * s26 - s1
     #6 * ((t97 + t83) * s25 + 34.0_dp * t9) - t25 * (39.0_dp * t15 + t9
     #8) - t43) * s34 + s16 * (((t21 + t96) * s25 - t92) * s25 + t11 * t
     #52) + ((t74 * t8 + t95) * s26 + t36 * (((t13 + s25) * s25 - t92) *
     # s25 + t93)) * s26) * s34 - t103 * (t13 - t88) + t104 * t40
      t78 = t7 * t48
      t79 = s34 * t7
      t95 = t45 * s16
      t112 = t11 * t41
      t113 = t112 * t1
      t114 = t40 * t26
      t115 = s16 + s25 + s26 - s34 + s56
      t116 = t31 - s25
      t117 = t36 * t9
      t118 = s16 * (t14 + t17)
      t119 = (t21 - s25) * s25
      t120 = 8.0_dp * t9
      t121 = t36 * t11
      t122 = t60 * t4
      t123 = t36 * s25
      t124 = (t123 - t17) * s25
      t125 = t9 * t8
      t126 = t45 * t4
      t116 = t11 * ((t14 - t21) * s25 + t92) + ((-(t2 * s16 * t7 + t25) 
     #* s34 + t9 * (t68 + t16) + ((t21 + s25 + s26) * s26 + t118) * s26)
     # * s34 - t7 * (t9 * (t68 - t54) + ((-s25 + t5) * s25 + t117) * s26
     # + t87 + t2 * t15 * t25)) * s34 + ((s56 * t58 + s16 * (t125 + t124
     #) - ((-t21 - s25 - t39 + s34) * s34 + t19 * t36 + t119 + t120) * s
     #34 - t43 + t126) * s56 - t2 * ((t46 * t7 + (t19 * t4 + t119 + t120
     #) * s26 + t121) * s34 + t44 - t45 * t116) + t46 * ((t22 + t39 + t1
     #6) * s26 + t118) + t9 * ((-t94 + t88) * s25 + t32) + t122) * s56 -
     # t59 + t62 + t60 * t116 + t95 * (t21 - t16)
      t118 = 8.0_dp * s16
      t119 = 15.0_dp * t10
      t33 = t33 + t16
      t127 = t41 - t9
      t128 = s25 * t127
      t129 = 28.0_dp * t9
      t130 = 16.0_dp * s25
      t131 = 23.0_dp * t9
      t132 = t21 + t105
      t133 = t132 * s26
      t134 = t9 * t4
      t135 = 27.0_dp * s16
      t136 = ((t135 + t16) * s25 + 32.0_dp * t9) * s25
      t137 = 62.0_dp * t11
      t138 = 19.0_dp * s25
      t139 = (38.0_dp * s16 + t138) * s25
      t140 = 63.0_dp * t9
      t141 = t140 + t139
      t142 = 23.0_dp * s16
      t143 = 20.0_dp * s25
      t144 = t142 + t143
      t26 = t44 + t26
      t145 = 20.0_dp * t9
      t146 = 33.0_dp * s16
      t147 = 17.0_dp * s25
      t148 = 20.0_dp * s16
      t149 = (t13 - s25) * s25 + t18
      t150 = 32.0_dp * s16
      t125 = ((t72 + t150) * s25 - t125) * s25
      t151 = 51.0_dp * t11
      t152 = t14 + t3
      t153 = t152 * s25
      t154 = t110 + t153
      t155 = t11 * s25
      t156 = (t53 + t54) * s25
      t157 = 24.0_dp * t9
      t158 = t45 * t149
      t159 = t59 + t62
      t122 = t122 * t33
      t160 = 10.0_dp * s25
      t161 = t160 * t10
      t162 = 18.0_dp * s26
      t163 = (-9.0_dp * t41 + t131) * s25
      t164 = t36 * t47
      t165 = 18.0_dp * t60
      t166 = 8.0_dp * s26
      t167 = (-t13 - t16) * s25
      t168 = 12.0_dp * t45
      t169 = t2 * s34
      t170 = t4 * t111
      t171 = t24 * t36
      t104 = ((((t21 + t105 - t169) * s34 - t170 + t171) * s56 + (t145 +
     # t167) * s25 + ((s34 * t36 - t142 - t143 - t166) * s34 + t8 * (7.0
     #_dp * s25 * s26 + t21 * s26 + t153 + 7.0_dp * t9)) * s34 - 11.0_dp
     # * t11 - t168) * s56 + s16 * (-t30 + t163) + (((t162 + t138 + t146
     #) * s34 - (69.0_dp * s16 + 60.0_dp * s25 + t109) * s26 - t139 - t1
     #40) * s34 + (t133 * t36 + 12.0_dp * t154) * s26 + t125 + t151) * s
     #34 + t42 - t164 - t165 - t126 * t33) * s56 + 9.0_dp * t11 * t127 +
     # t2 * (t49 * ((t90 + t138 + t146) * s26 + t156 + t157 + t46) + t10
     #2 - t158) - 12.0_dp * t159 + ((-t46 * (t123 + t109 + t108) - t104 
     #- t137 - t136 - s26 * (t4 * t144 * s26 + t2 * t141)) * s34 - t2 * 
     #(-(t151 + t125) * s26 + t155) + 38.0_dp * t10 + 12.0_dp * t87 * t1
     #5 + 12.0_dp * t154 * t25 + t8 * t132 * t28) * s34 - t122 + t161
      t108 = t95 * ((t118 - t16) * s25 + t91)
      t125 = t62 * t33
      t132 = s16 * t10 * s25
      t139 = t114 * t4
      t104 = t113 * t20 - t2 * (-t78 + t12) + (((-((17.0_dp * s26 + t22)
     # * s16 + t36 * t89) * s34 + t25 * (t138 + t146) + t36 * (((t21 + s
     #25) * s25 + t120 + t25) * s26 + t87) + t9 * (t147 + t148)) * s34 -
     # t2 * t26 - t9 * ((t105 + t53) * s25 + t145) - ((s26 * t144 + t141
     #) * s26 + t136 + t137) * s26) * s34 + t7 * (s16 * (((-t31 + t105) 
     #* s25 + t134) * s25 + t93) + ((t133 + (s16 + t130) * s25 + t131) *
     # s26 + (16.0_dp * t41 + t129) * s16 + t128 * t20) * s26)) * s34 + 
     #t104 * s56 + t25 * ((t41 * (-t17 + s25) + 23.0_dp * t11) * s25 - t
     #119) + t42 * t9 - t108 - t125 + t132 - t139
      t131 = s16 + s26 + s56
      t133 = (t5 + t16) * s25 + t110
      t136 = t133 * s26
      t137 = 16.0_dp * s16
      t138 = ((39.0_dp * s16 + t88) * s25 + 48.0_dp * t9) * s25 + 102.0_
     #dp * t11
      t140 = 26.0_dp * s16
      t141 = (t64 + t140) * s25 + 42.0_dp * t9
      t142 = t21 + t54
      t143 = 22.0_dp * t9
      t144 = (52.0_dp * s16 + t97) * s25
      t146 = 78.0_dp * t9
      t151 = t146 + t144
      t153 = s25 + t3
      t154 = 40.0_dp * t11
      t99 = (t123 + t99) * s25
      t172 = 55.0_dp * t9
      t173 = (((t63 + t16) * s25 + t143) * s25 - t70) * s25
      t174 = 63.0_dp * t10
      t91 = (t41 * t8 + t91) * s16
      t175 = t91 + t128
      t176 = 37.0_dp * s16
      t177 = 16.0_dp * t9
      t178 = s16 * t111
      t179 = t43 + t11
      t180 = t9 * s25
      t108 = -t60 * t149 + t180 * t179 - t108 - t125
      t33 = (-(t178 * t4 - ((-t142 + s34) * s34 + t133) * s34 - t117 * s
     #25) * s56 - s16 * (s16 * ((t105 - t148) * s25 + 11.0_dp * t9) + t1
     #68) - ((t4 * t46 + s16 * (42.0_dp * s16 + t107) + (t109 + t64 + t1
     #40) * s25) * s34 - t8 * (t91 + t128 + t136)) * s34 - t2 * (-t49 * 
     #(t72 + t22 + t29) + t44)) * s56 - s16 * (s16 * (t30 - t163) + t165
     #) + t4 * (t46 * (-((t63 + t54) * s25 + t177) * s25 - t76 - s26 * (
     #t2 * t142 * s26 + t141) + t49) - t95 * t33) + (t46 * (-(t90 + t64 
     #+ t135) * s34 + (60.0_dp * s16 + 30.0_dp * s25 + t51) * s26 + t144
     # + t146) + (t136 * t36 + 12.0_dp * t175) * s26 + t173 + t174) * s3
     #4 + t102
      t33 = t33 * s56 - (-t161 + t122) * s16 + 9.0_dp * t10 * t127 + t2 
     #* ((t102 - t158) * s16 + (t46 * (t46 * (t13 + t39 + t16) + s16 * (
     #(t176 + t147) * s25 + 50.0_dp * t9) + (15.0_dp * t153 * s26 + t151
     #) * s26 + t2 * (t43 + t28)) + s16 * ((t41 * (t118 + t16) - t75) * 
     #s25 + 23.0_dp * t10) + s26 * (t25 * t2 * t133 + t173 + t174) + t38
     # * t175) * s34) - (t46 * ((54.0_dp * s16 + t100 + t90) * s26 + t17
     #2 + t99 + t46) + t25 * t4 * t141 + s16 * (((t176 + t130) * s25 + 4
     #6.0_dp * t9) * s25 + 95.0_dp * t11) + (t36 * t138 / 3.0_dp + t142 
     #* t8 * t25) * s26 + t42) * t46 - 12.0_dp * s16 * t159
      t33 = s16 * t108 + (((((-(s26 + t3) * s34 + s16 * (t105 + t94) + (
     #t14 + t83 + t39) * s26) * s34 - s16 * ((t135 + t88) * s25 + 30.0_d
     #p * t9) - t28 * t4 - s26 * ((t64 + t135) * s26 + t172 + t99)) * s3
     #4 + s16 * (((t72 + t148) * s25 + 38.0_dp * t9) * s25 + t154) + ((t
     #151 + t25) * s26 + ((34.0_dp * s16 + t14) * s25 + 74.0_dp * t9) * 
     #s25 + 100.0_dp * t11) * s26 + 10.0_dp * t153 * t28) * s34 - t9 * (
     #((t105 + t118) * s25 + t143) * s25 + 30.0_dp * t11) - (((s26 * t14
     #2 + t141) * s26 + t138) * s26 + (((s25 + t137) * s25 + 37.0_dp * t
     #9) * s25 + 46.0_dp * t11) * s25 + 95.0_dp * t10) * s26 - t102) * s
     #34 + t7 * (s16 * (((-t120 - t167) * s25 + t75) * s25 + 12.0_dp * t
     #10) + ((((t14 + t83) * s25 - t110) * s25 + 29.0_dp * t11 + t136) *
     # s26 + 9.0_dp * t24 * t127 + t41 * (t120 + t67) + 34.0_dp * t10) *
     # s26)) * s34 + t33 * s56
      t83 = (-t118 - t160) * s25 + t18
      t91 = t32 - t156
      t99 = t21 + t123
      t100 = 18.0_dp * s25
      t102 = -t72 - t13
      t52 = ((-s16 * t52 - t96) * s25 - t129) * s25
      t108 = t154 + t52
      t109 = (27.0_dp * s25 + t140) * s25
      t110 = 31.0_dp * t9
      t120 = t110 + t109
      t122 = 25.0_dp * s16
      t123 = -t100 - t122
      t125 = 52.0_dp * t9
      t77 = (t77 + t137) * s25
      t127 = -t117 + t77
      t128 = t105 + t17
      t129 = t20 * t11
      t118 = ((t118 + t97) * s25 + t157) * s25 + t129
      t130 = (s16 + t88) * s25 + t50
      t133 = -10.0_dp * t102
      t72 = (((t72 + t13 - t169) * s34 - t170 + t171) * s56 + s16 * ((t3
     #1 + t105) * s25 - t32) - t36 * (t43 - t49) + s34 * (-(10.0_dp * s2
     #6 + t122 + t100) * s34 + s16 * (24.0_dp * s16 + 35.0_dp * s26) + (
     #25.0_dp * s26 + t100 + t3) * s25) - 15.0_dp * t45) * s56 + s16 * (
     #((-t88 + t150) * s25 - 25.0_dp * t9) * s25 + t129) - t4 * (-t49 * 
     #(t105 + t166 + t17) + t42) - t164 + s34 * (-((100.0_dp * s16 + 72.
     #0_dp * s25 + t107) * s26 + t109 + t110) * s34 + (15.0_dp * t41 + t
     #157) * s25 + 8.0_dp * s26 * t130 + t133 * t25 + t129 + 8.0_dp * t8
     #7) - 30.0_dp * t60 - t45 * t8 * t99
      t52 = t72 * s56 + s16 * ((((-t160 + t94) * s25 + t55) * s25 - 38.0
     #_dp * t11) * s25 + t119) + t2 * t48 + (((-(t13 + t96 + t162) * s34
     # + (81.0_dp * s16 + 63.0_dp * s25 + 36.0_dp * s26) * s26 - t117 + 
     #t77) * s34 + ((t36 * t123 - t107) * s26 - t120 * t4) * s26 + t154 
     #+ t52) * s34 + (((22.0_dp * s16 + t16) * s25 - t85) * s25 + 58.0_d
     #p * t11) * s25 + ((t133 * s26 + 12.0_dp * t130) * s26 + t118 * t4)
     # * s26 - 44.0_dp * t10) * s34 - 30.0_dp * t62 + t126 * t91 - t38 *
     # t40 * t99
      t5 = (((((-t13 - t96) * s26 + t3 * (t17 + s25)) * s26 - t36 * (-t1
     #80 + t28) + t93) * s34 + t28 * (t128 * t4 + t51) + t9 * ((-t14 - t
     #86) * s25 - t145) - s26 * (s16 * (t125 - t124) - s26 * t127)) * s3
     #4 + t10 * (t100 + t148) - t2 * (t103 * t15 + t27) + (s16 * (((t21 
     #- t160) * s25 - t134) * s25 + 68.0_dp * t11) + ((s26 * t123 - t120
     #) * s26 + t108) * s26) * s26) * s34 - t7 * (10.0_dp * t9 * t179 + 
     #(((t102 * s26 + (-t100 + t5) * s25 - t35) * s26 + ((t22 - t97) * s
     #25 - 27.0_dp * t9) * s25 + 12.0_dp * t11) * s26 + ((t177 + t167) *
     # s25 - 31.0_dp * t11) * s25 + 32.0_dp * t10) * s26 + t24 * ((-t177
     # - t67) * s25 + t121))
      t5 = t5 * s34 + (t52 * s56 - 11.0_dp * t103 * t111 + t11 * (27.0_d
     #p * t43 + t93) + t2 * (t47 * ((-t13 - t96) * s26 + t24 - 9.0_dp * 
     #t25 + 9.0_dp * t9) + t95 * t83) - 15.0_dp * t114 - 15.0_dp * t132 
     #+ (((t46 * (t98 - t3) + (-t125 + t106) * s16 + (t2 * t127 + 24.0_d
     #p * t25) * s26 + 9.0_dp * t128 * t25 - 9.0_dp * t180) * s34 + (t12
     #3 * t8 * t25 + t2 * t108) * s26 - t4 * (t120 * t25 + t155) + t9 * 
     #(t20 * t41 + 68.0_dp * t9) - 10.0_dp * t26) * s34 + s16 * ((t43 * 
     #t8 + 25.0_dp * t11) * s25 - 42.0_dp * t10) + (t25 * (-t102 * t20 *
     # s26 + 8.0_dp * t130) + (((44.0_dp * s16 + t14) * s25 - t125) * s2
     #5 + 116.0_dp * t11) * s25 - 88.0_dp * t10) * s26 - t4 * (-t118 * t
     #25 + t59)) * s34 + t60 * (t4 * t91 - t98 * t99)) * s56 + t45 * ((s
     #16 * t83 + (-s26 * t99 + t91) * s26) * s26 + t9 * ((t21 - t64) * s
     #25 + t81))
      t26 = t1 * s26
      t43 = 0.1e1_dp / t58
      t52 = 0.1e1_dp / s16
      t72 = 0.1e1_dp / t115 ** 2.0_dp
      t65 = 0.1e1_dp / t65 ** 2.0_dp
      t5 = -s25 * (t2 * (t11 * t10 - t48 * ((s16 - s26) * s26 + t9)) + t
     #5 + t112 * ((-t14 + t94) * s25 - t81) - t4 * t40 * t27) * MB1110(1
     #) + t131 * (-t20 * t101 * t41 + (t12 * t2 + t139) * s16 - t33) * M
     #B1011(1)
      t5 = t104 * MB1001(0) + t5 * t52 + t8 * t113 - t73 * s56 + t25 * (
     #(t41 * (-t68 + s25) + t37) * s25 - 9.0_dp * t10) - t46 * (t9 * ((t
     #17 + t16) * s25 + t81) + ((t80 + t25) * s26 + ((s25 + t53) * s25 +
     # t71) * s25 + t76) * s26 + t44 + t82 * t28) - t47 * ((t21 + t90) *
     # s16 + t4 * t89) + t49 * (t4 * (t87 + t28) + t9 * (t22 + t88) + s2
     #6 * ((t64 + t86) * s26 + t84 + t85)) - t9 * (-t42 + t10) + t78 - t
     #62 * t69 + t79 * (t20 * s16 * (t41 * (s25 - s16) + t11) + ((t74 + 
     #(t64 - s16) * s25 + t71) * s26 + ((t22 + t54) * s25 - t32) * s25 +
     # t70) * s26 + t75 * s25) - t95 * (t20 * t66 - t67) - t114 * t2
      t3 = t65 * t43 * ((s16 * (t9 * ((t64 - t63) * s25 + t32) + t57 * (
     #-s25 + t3)) + t2 * (t58 * s56 * t61 + t62) - t4 * (-t60 * t6 + t59
     #) - ((((t31 + t29) * s26 + t32) * s34 - t8 * s16 * ((t21 + t16) * 
     #s26 + t24) - (t23 + t29) * t25 - t30) * s34 + t7 * (s16 * ((-t17 +
     # t16) * s25 + t18) + (t19 * t8 + (-t14 + t13) * s25 + 13.0_dp * t9
     #) * s26)) * s34 - (-(((t53 - t54) * s25 - t55) * s25 - t2 * t49 - 
     #s34 * (-(t51 + t22 + s25) * s34 + 12.0_dp * t19 + t34 + t35) + t56
     # + t57) * s56 + t2 * s34 * (-((t14 + t22) * s16 + (t23 + t39) * s2
     #6) * s34 + t38 * t15 + t37 + (t35 + t34) * s26 - t24 * t15) + t36 
     #* (-t45 * (s26 + t6) + t44) - t9 * ((-30.0_dp * s16 + 24.0_dp * s2
     #5) * s25 + t50) + t8 * t7 * t49) * s56) * t52 + t116 * I300s16(1) 
     #+ t5 * t72)
      result = t3 / 4.0_dp + t52 * (t65 * (t131 * (t2 * ((-t79 + t26) * 
     #s16 + ((-s25 - s34 + s16) * s16 + t26) * s56) + (s16 * t1 + t46) *
     # s16 + t1 * (t61 + t25)) * MB1111D2(0) + s25 * s34 * (-t2 * s16 * 
     #(t180 - t45) + ((s26 * t7 + t9) * s34 - t7 * (t92 + (s16 + t16) * 
     #s26 + t24)) * s34 + ((-t2 * s25 * (s16 + s34) + (s34 - s16) * s34 
     #+ t41 + t9) * s56 + t2 * (t45 + t178) - s34 * (s26 * t152 - (s16 +
     # t29) * s34 + t4 * t66) - t180 * t8) * s56 + t111 * t9 + t60) * MB
     #1101(0) * t43) + s25 * (s16 + s26 - s34 + s56) * MB1110(0) * t72) 
     #/ 2.0_dp

           intHs16s25s26s34s56x1321D6eps0 = result
       end function intHs16s25s26s34s56x1321D6eps0

       function intHs16s25s26s34s56x1321D6eps1()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1321D6eps1
           complex(dp) ::  t1,t2

           complex(dp) :: result

      t1 = s16 + s25 + s26 - s34 + s56
      t2 = 0.1e1_dp / 0.4e1_dp
      t1 = 0.1e1_dp / t1 ** 2.0_dp
      result = 0.1e1_dp / s16 * t1 * (t2 * (s16 - s25 + s26 - s34 + s56)
     # + s25 * (s16 + s26 - s34 + s56) * MB1110(1) / 2.0_dp)

           intHs16s25s26s34s56x1321D6eps1 = result
       end function intHs16s25s26s34s56x1321D6eps1

       function intHs16s25s26s34s56x1411D6eps0()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1411D6eps0
           complex(dp) ::  t1,t10,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109
           complex(dp) ::  t11,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t12
           complex(dp) ::  t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t13,t130
           complex(dp) ::  t131,t132,t133,t134,t135,t136,t137,t138,t139,t14,t140,t141
           complex(dp) ::  t142,t143,t144,t145,t146,t147,t148,t149,t15,t150,t151,t152
           complex(dp) ::  t153,t154,t155,t156,t157,t158,t159,t16,t160,t161,t162,t163
           complex(dp) ::  t164,t165,t166,t167,t168,t169,t17,t170,t171,t172,t173,t174
           complex(dp) ::  t175,t176,t177,t178,t179,t18,t180,t181,t182,t183,t184,t185
           complex(dp) ::  t186,t187,t188,t189,t19,t190,t191,t192,t193,t194,t195,t196
           complex(dp) ::  t197,t198,t199,t2,t20,t200,t201,t202,t203,t204,t205,t206
           complex(dp) ::  t207,t208,t209,t21,t210,t211,t212,t213,t214,t215,t216,t217
           complex(dp) ::  t218,t219,t22,t220,t221,t222,t223,t224,t225,t226,t227,t228
           complex(dp) ::  t229,t23,t230,t231,t232,t233,t234,t235,t236,t237,t238,t239
           complex(dp) ::  t24,t240,t241,t242,t243,t244,t245,t246,t247,t248,t249,t25
           complex(dp) ::  t250,t251,t252,t253,t254,t255,t256,t257,t258,t259,t26,t260
           complex(dp) ::  t261,t262,t27,t28,t29,t3,t30,t31,t32,t33,t34,t35
           complex(dp) ::  t36,t37,t38,t39,t4,t40,t41,t42,t43,t44,t45,t46
           complex(dp) ::  t47,t48,t49,t5,t50,t51,t52,t53,t54,t55,t56,t57
           complex(dp) ::  t58,t59,t6,t60,t61,t62,t63,t64,t65,t66,t67,t68
           complex(dp) ::  t69,t7,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79
           complex(dp) ::  t8,t80,t81,t82,t83,t84,t85,t86,t87,t88,t89,t9
           complex(dp) ::  t90,t91,t92,t93,t94,t95,t96,t97,t98,t99

           complex(dp) :: result

      t1 = s16 - s25
      t2 = 15.0_dp * s16
      t3 = 8.0_dp * s25
      t4 = t2 + t3
      t5 = 4.0_dp
      t6 = s16 ** 2.0_dp
      t7 = t6 ** 2.0_dp
      t8 = t7 ** 2.0_dp
      t9 = s16 * t6
      t10 = t9 * t7
      t11 = t6 * t7
      t12 = s16 * t7
      t13 = t5 * s25
      t14 = 29.0_dp * s16
      t15 = 43.0_dp * t6
      t16 = s16 + s25
      t17 = 2.0_dp
      t18 = t17 * s25
      t19 = 46.0_dp * t6
      t20 = t5 * s16
      t21 = s25 + t20
      t22 = 11.0_dp * s25
      t23 = 13.0_dp * s16 + t22
      t24 = s16 + s26
      t25 = 9.0_dp
      t26 = 43.0_dp * s16
      t27 = t25 * s25
      t28 = 30.0_dp * t6
      t29 = 73.0_dp * s25
      t30 = 141.0_dp * t9
      t31 = 98.0_dp * s16
      t32 = 82.0_dp * s25
      t33 = 96.0_dp * t6
      t34 = 23.0_dp * s16
      t35 = 27.0_dp * s25
      t36 = t35 + t34
      t37 = s26 ** 2.0_dp
      t38 = t37 ** 2.0_dp
      t39 = s26 * t38
      t40 = s26 * t37
      t41 = s25 ** 2.0_dp
      t42 = t41 ** 2.0_dp
      t43 = s25 * t41
      t44 = t43 * t42
      t45 = t41 * t42
      t46 = s25 * t42
      t47 = t41 + t6
      t48 = t6 * s25
      t49 = t48 * t47
      t50 = 52.0_dp * t6
      t51 = t17 * s26
      t52 = 762.0_dp * s16
      t53 = ((219.0_dp * s25 + t52) * s25 + 846.0_dp * t6) * s25 + 711.0
     #_dp * t9
      t54 = (750.0_dp * s16 + 492.0_dp * s25) * s25 + 714.0_dp * t6
      t55 = 10.0_dp
      t56 = t55 * t36 * s26
      t57 = 18.0_dp * s25
      t58 = 49.0_dp
      t59 = 42.0_dp * t6
      t60 = t58 * t9
      t61 = (((309.0_dp * s16 + 200.0_dp * s25) * s25 + t59) * s25 + t60
     #) * s25
      t62 = 408.0_dp * t7
      t63 = t62 + t61
      t64 = 6.0_dp
      t65 = t64 * t6
      t66 = ((69.0_dp * s16 + 168.0_dp * s25) * s25 + t65) * s25
      t67 = 157.0_dp * t9
      t68 = t67 + t66
      t69 = t64 * s16
      t70 = 39.0_dp * s25
      t71 = 19.0_dp * t6
      t72 = (t69 + t70) * s25 + t71
      t73 = t9 * t41
      t74 = 3.0_dp
      t75 = 66.0_dp * s25
      t76 = (361.0_dp * s16 + t75) * s25 + 191.0_dp * t6
      t77 = t63 * t74
      t78 = t55 * t72
      t79 = 177.0_dp * t7
      t80 = 475.0_dp * t12
      t81 = ((409.0_dp * s16 + 290.0_dp * s25) * s25 + 247.0_dp * t6) * 
     #s25
      t82 = 340.0_dp * t9
      t83 = t81 + t82
      t84 = 119.0_dp * s16
      t85 = 172.0_dp * s25
      t86 = 95.0_dp * t6
      t87 = (t84 + t85) * s25 + t86
      t88 = 90.0_dp * s16
      t89 = 250.0_dp * s25 + t88
      t90 = t89 * s26
      t91 = t83 * t74
      t92 = (774.0_dp * s16 + 705.0_dp * s25) * s25 + 486.0_dp * t6
      t93 = 88.0_dp * s25
      t94 = 45.0_dp * s16
      t95 = t94 + t93
      t96 = 60.0_dp * s26
      t97 = 303.0_dp * s25
      t98 = 189.0_dp * s16
      t99 = 108.0_dp * s26
      t100 = 102.0_dp * s25
      t101 = 27.0_dp * s16
      t102 = 54.0_dp * s26
      t103 = 38.0_dp * s25
      t104 = t1 ** 2.0_dp
      t105 = t104 ** 2.0_dp
      t106 = t1 * t104
      t107 = 76.0_dp
      t108 = 12.0_dp
      t109 = s34 ** 2.0_dp
      t110 = t109 ** 2.0_dp
      t111 = s34 * t109
      t112 = t111 * t110
      t113 = t109 * t110
      t114 = s34 * t110
      t115 = 37.0_dp * t6
      t116 = t108 * t105
      t117 = 70.0_dp * t105
      t118 = t74 * s26
      t119 = t64 * t112
      t120 = 0.2e1_dp / 0.3e1_dp * t54
      t121 = 254.0_dp * s16
      t122 = t68 * t5
      t123 = t5 * t87
      t124 = 110.0_dp * t6
      t125 = t105 * s26
      t126 = 5.0_dp * s26
      t127 = t126 * t72
      t128 = 45.0_dp * s26
      t129 = 30.0_dp * s26
      t130 = 34.0_dp
      t131 = 16.0_dp * s25
      t132 = s16 * t130
      t133 = (t131 - t132) * s25
      t134 = 24.0_dp * t6
      t89 = t89 / 10.0_dp
      t135 = t7 + t42
      t136 = s16 * s25
      t137 = t136 * t47
      t138 = 36.0_dp
      t139 = t138 * s25
      t140 = t9 * t64
      t141 = ((s25 * t76 - t140) * s25 + t79) * s25 + t80
      t142 = ((((201.0_dp * s16 + 100.0_dp * s25) * s25 + 38.0_dp * t6) 
     #* s25 + 18.0_dp * t9) * s25 + 162.0_dp * t7) * s25 + 265.0_dp * t1
     #2
      t143 = 260.0_dp * s25
      t144 = (((988.0_dp * s16 + t143) * s25 + 770.0_dp * t6) * s25 + 77
     #8.0_dp * t9) * s25 + 992.0_dp * t7
      t52 = ((300.0_dp * s25 + t52) * s25 + 620.0_dp * t6) * s25 + 568.0
     #_dp * t9
      t145 = 215.0_dp * t9
      t146 = ((199.0_dp * s16 + 140.0_dp * s25) * s25 + 235.0_dp * t6) *
     # s25 + t145
      t147 = 204.0_dp * s25
      t148 = (314.0_dp * s16 + t147) * s25 + 194.0_dp * t6
      t149 = t97 + t98
      t150 = 64.0_dp * s25
      t151 = (75.0_dp * s16 + t150) * s25 + t107 * t6
      t152 = s16 * t108
      t153 = 11.0_dp * t6
      t154 = (((((-((s26 * t107 + t152) * s25 + 54.0_dp * s26 * t24 + t1
     #53) * s34 + s16 * t151 + 72.0_dp * t40 + s26 * (s26 * t149 + t148)
     #) * s34 - s16 * t146 - t40 * (t5 * t95 + t129) - s26 * (s26 * t92 
     #+ t52)) * s34 + s16 * ((((270.0_dp * s16 + 160.0_dp * s25) * s25 +
     # 204.0_dp * t6) * s25 + 334.0_dp * t9) * s25 + 320.0_dp * t7) + t3
     #7 * (5.0_dp * t89 * t37 + t91) + s26 * (t123 * t37 + t144)) * s34 
     #- s16 * t142 - (((t122 + t127) * s26 + t77) * s26 + t141 * t17) * 
     #s26) * s34 + t104 * (s16 * ((((143.0_dp * s16 + 32.0_dp * s25) * s
     #25 + 198.0_dp * t6) * s25 + 191.0_dp * t9) * s25 + 116.0_dp * t7) 
     #+ ((s26 * t120 + t53) * s26 + (((354.0_dp * s16 + t139) * s25 + 63
     #6.0_dp * t6) * s25 + 672.0_dp * t9) * s25 + 466.0_dp * t7) * s26 +
     # 5.0_dp * t36 * t38)) * s34
      t155 = t105 * t16 * ((t14 + t18) * s25 + t19)
      t156 = t10 * s25
      t157 = 13.0_dp * t6
      t158 = t74 * t105
      t159 = t9 * s25
      t160 = -86.0_dp
      t161 = 68.0_dp
      t162 = 44.0_dp
      t163 = 80.0_dp * t7
      t164 = 26.0_dp * s25
      t165 = 41.0_dp * s16
      t166 = 69.0_dp * t6
      t167 = t64 * t39
      t168 = t64 * s25
      t148 = t6 * ((19.0_dp * s16 + t168) * s25 + t134) + 18.0_dp * t38 
     #+ (s16 * t151 + (s26 * t149 / 3.0_dp + t148 / 2.0_dp) * s26) * s26
      t149 = t108 * s25
      t151 = 11.0_dp * s16
      t169 = t5 * t9
      t66 = (t110 * t64 - 28.0_dp * t137 - ((t89 * s34 - t72) * s34 + t1
     #04 * t36) * s34 + 7.0_dp * t135 + t59 * t41) * s56 + (((t133 - t13
     #4) * s25 + 116.0_dp * t9) * s25 - 104.0_dp * t7) * s25 - 18.0_dp *
     # t114 + 30.0_dp * t12 - (((-(t94 + t129 + t93) * s34 + s16 * (95.0
     #_dp * s16 + t128) + (125.0_dp * s26 + t84 + t85) * s25) * s34 - t1
     #27 - t66 - t67) * s34 + t104 * (s16 * (115.0_dp * s26 + t84) + (12
     #5.0_dp * s16 + 135.0_dp * s26 + t32) * s25)) * s34 + 35.0_dp * t12
     #5
      t61 = t66 * s56 + ((((s25 * t23 - t124) * s25 + 90.0_dp * t9) * s2
     #5 + 95.0_dp * t7) * s25 - 151.0_dp * t12) * s25 + 52.0_dp * t11 + 
     #18.0_dp * t113 - (((((63.0_dp * s16 + 101.0_dp * s25 + 72.0_dp * s
     #26) * s34 - (258.0_dp * s16 + 235.0_dp * s25) * s25 - (180.0_dp * 
     #s16 + 352.0_dp * s25 + t96) * s26 - 162.0_dp * t6) * s34 + s26 * (
     #t123 + t90) + t81 + t82) * s34 - (t78 * s26 + t122) * s26 - t61 - 
     #t62) * s34 + t104 * (((t121 + t29) * s25 + 282.0_dp * t6) * s25 + 
     #237.0_dp * t9 + s26 * (t56 + t120))) * s34 + 8.0_dp * t125 * t4 + 
     #t117 * t37
      t53 = t61 * s56 + ((((((t18 + t34) * s25 - t115) * s25 - t107 * t9
     #) * s25 + 144.0_dp * t7) * s25 + 7.0_dp * t12) * s25 - 109.0_dp * 
     #t11) * s25 + 46.0_dp * t10 - (((((-(t102 + t103 + t101) * s34 + (1
     #57.0_dp * s16 + t100) * s25 + (t97 + t99 + t98) * s26 + 97.0_dp * 
     #t6) * s34 - ((381.0_dp * s16 + 150.0_dp * s25) * s25 + 310.0_dp * 
     #t6) * s25 - ((t64 * t95 + t96) * s26 + t92) * s26 - 284.0_dp * t9)
     # * s34 + (((494.0_dp * s16 + 130.0_dp * s25) * s25 + 385.0_dp * t6
     #) * s25 + 389.0_dp * t9) * s25 + ((t64 * t87 + t90) * s26 + t91) *
     # s26 + 496.0_dp * t7) * s34 - (t41 * t76 + t79) * s25 - (t78 * t37
     # + t77) * s26 - t64 * (t37 * t68 - t73) - t80) * s34 + t104 * ((((
     #177.0_dp * s16 + t57) * s25 + 318.0_dp * t6) * s25 + 336.0_dp * t9
     #) * s25 + 233.0_dp * t7 + ((t56 + t54) * s26 + t53) * s26)) * s34 
     #+ t116 * t4 * t37 + t117 * t40 + t118 * t105 * t21 * t23 - t119
      t53 = t53 * s56 + (-t108 * t112 + t155 * t17) * s26 - t154 + 21.0_
     #dp * t8 - 41.0_dp * t156 + 8.0_dp * t105 * t4 * t40 + 35.0_dp * t1
     #05 * t38 - 17.0_dp * t11 * t41 - t60 * t46 + 62.0_dp * t12 * t43 +
     # t157 * t45 + 7.0_dp * t7 * t42 + t20 * t44 + t158 * t21 * t23 * t
     #37
      t54 = t9 * t45
      t56 = t125 * s16
      t14 = t56 * (((t13 + t14) * s25 + t15) * s25 + 21.0_dp * t9)
      t61 = t104 * t24
      t29 = t61 * (s16 * (t41 * (14.0_dp * t41 + t50) + 24.0_dp * t7) + 
     #39.0_dp * t49 + ((t36 * s26 + (t31 + t32) * s25 + t33) * s26 + ((1
     #72.0_dp * s16 + t29) * s25 + 184.0_dp * t6) * s25 + t30) * t37 + t
     #51 * t16 * (((t26 + t27) * s25 + t28) * s25 + 46.0_dp * t9)) * s34
      t14 = t53 * s56 + t109 * (t130 * t159 * ((t41 + t6) * s25 + t9) + 
     #t6 * (t43 * (t138 * t41 - t65) + 60.0_dp * t12) + (s16 * t142 + ((
     #(s26 * t72 + t68) * s26 + t63) * s26 + t141) * s26) * s26) + t110 
     #* (t6 * (((t164 + t165) * s25 + t166) * s25 + 60.0_dp * t9) + (s16
     # * t146 + ((s26 * t95 + t92 / 3.0_dp) * s26 + t52 / 2.0_dp) * s26)
     # * s26 + t167) + t111 * (s16 * (s16 * ((((-s25 * t162 - 50.0_dp * 
     #s16) * s25 - t161 * t6) * s25 + t160 * t9) * s25 - t163) - t51 * (
     #(((135.0_dp * s16 + 80.0_dp * s25) * s25 + 102.0_dp * t6) * s25 + 
     #167.0_dp * t9) * s25 + 160.0_dp * t7)) - (((t89 * s26 + t87) * s26
     # + t83) * s26 + t144 / 2.0_dp) * t37) + t113 * (18.0_dp * t40 + s2
     #6 * (s16 * (t151 + t149) + s26 * (t103 + t101)) + t169) - t114 * t
     #148 + t54 + t105 * t21 * t23 * t40 + t155 * t37 + t14 - t29
      t21 = t8 * s25
      t23 = t105 * t39
      t29 = 14.0_dp * s16
      t32 = 15.0_dp * s25
      t36 = t74 * t9
      t52 = t7 * t41
      t53 = s34 * t16
      t62 = -t17 * (t136 + t53) + t6 + t41 + t109
      t63 = s16 + s25 + s26 - s34 + s56
      t66 = (s34 - t16 - s26) * s26 + (-t51 - s16 - s25 + s34 - s56) * s
     #56 - t136
      t67 = 7.0_dp * s16
      t68 = 8.0_dp * s16
      t72 = s25 + t68
      t76 = 16.0_dp * s16
      t77 = 7.0_dp * s25
      t78 = t25 * s16
      t80 = t78 + t149
      t81 = 5.0_dp * s16
      t82 = 24.0_dp * s25
      t83 = t17 * t6
      t85 = 22.0_dp
      t87 = s16 * t85
      t89 = t41 * (t57 + t87) + 32.0_dp * t9
      t90 = t74 * t41
      t91 = t6 + t90
      t92 = s16 * (((t67 + t82) * s25 + t83) * s25 + 39.0_dp * t9)
      t93 = 7.0_dp * t41
      t94 = t108 * t9
      t95 = 20.0_dp * s25
      t97 = 18.0_dp * s16
      t98 = t85 * s25
      t117 = (t97 + t98) * s25 + 15.0_dp * t6
      t120 = s16 + t13
      t122 = s16 * ((t95 + t152) * s25 + 25.0_dp * t6)
      t123 = t6 * t41
      t127 = 13.0_dp * s25
      t141 = s16 * t16
      t142 = t17 * s16
      t144 = s25 + t142
      t146 = t80 / 3.0_dp
      t148 = 24.0_dp * s16
      t154 = 27.0_dp * t6
      t155 = t25 * s26
      t170 = t6 * t42
      t171 = 26.0_dp * t6
      t172 = t125 * t25
      t173 = 27.0_dp * s26
      t174 = t138 * s26
      t175 = t6 * t162
      t176 = t65 * t41
      t177 = s56 ** 2.0_dp
      t178 = s56 * t177 ** 2.0_dp
      t133 = (s16 * ((((-26.0_dp * s16 + t13) * s25 + t175) * s25 - 31.0
     #_dp * t9) * s25 + 8.0_dp * t7) + t109 * (-t111 * t74 + t17 * (t25 
     #* (s26 * t6 + t43) + t41 * (t151 + t173) + 16.0_dp * t9)) + t172 +
     # (-t104 * ((t101 + t139) * s26 + 28.0_dp * t141 + t93) + t109 * ((
     #t155 + t69 + t127) * s34 - s16 * (t2 + t155) - (t97 + t174 + t98) 
     #* s25)) * s34 + t46) * s56 + (-t114 * t64 + t172) * s26 + t17 * (s
     #16 * t46 + t125 * t72) + t9 * ((-t133 - t171) * s25 + 7.0_dp * t9)
     # - (((-((t164 + t155 + t152) * s26 + t141 * t64) * s34 + (t17 * t1
     #17 + t155 * t120) * s26 + t122) * s34 - (18.0_dp * t91 * s26 + t17
     # * t89) * s26 - t92) * s34 + t104 * (s16 * ((t149 + t148) * s25 + 
     #t154) + (14.0_dp * t144 ** 2.0_dp + t155 * t146) * s26)) * s34 - t
     #170
      t164 = t6 * t46
      t172 = t43 * (-s25 - s16) + t7
      t179 = t12 * s25
      t180 = (t18 + t76) * s25 + t153
      t181 = 17.0_dp * s16
      t182 = (s25 + t181) * s25 + 21.0_dp * t6
      t183 = t181 + t168
      t184 = 31.0_dp * s25
      t185 = t55 * s25
      t186 = t6 * t55
      t187 = 78.0_dp * t6
      t188 = -105.0_dp * t9
      t189 = (((s16 + 29.0_dp * s25) * s25 - t187) * s25 + t188) * s25
      t190 = t189 + t79
      t191 = 51.0_dp
      t192 = 81.0_dp * s16
      t193 = 30.0_dp * s25
      t194 = ((t193 - t192) * s25 - 88.0_dp * t6) * s25
      t195 = t191 * t9
      t196 = t195 + t194
      t197 = 21.0_dp * s16
      t198 = (t27 + t197) * s25
      t199 = 14.0_dp * t6
      t200 = t199 + t198
      t201 = s25 * t1
      t202 = 16.0_dp * t6
      t203 = 93.0_dp
      t204 = ((s16 * t203 - t131) * s25 + 92.0_dp * t6) * s25
      t205 = 138.0_dp * t9
      t206 = -t205 + t204
      t207 = (134.0_dp * s16 + t139) * s25
      t208 = t207 + t115
      t209 = 52.0_dp * s25
      t210 = 57.0_dp * s16
      t211 = t210 + t209
      t212 = 48.0_dp
      t213 = t212 * s25
      t214 = t51 * s16
      t215 = 23.0_dp * t6
      t216 = 21.0_dp * s25
      t160 = (s16 * t160 - t216) * s25
      t217 = t160 + t154
      t218 = 72.0_dp * s16
      t219 = -s25 * t161 - t218
      t220 = s16 * t41
      t221 = 30.0_dp * s16
      t222 = 41.0_dp * t6
      t223 = ((t13 + t87) * s25 + 20.0_dp * t6) * s25
      t224 = 33.0_dp * t9
      t225 = (s16 * t138 + t184) * s25
      t226 = 60.0_dp * t6
      t227 = t226 + t225
      t228 = 28.0_dp * s16
      t229 = t213 + t228
      t230 = 15.0_dp * s26
      t231 = 32.0_dp * s16
      t232 = 40.0_dp * t6
      t233 = 24.0_dp * t9
      t234 = 8.0_dp * t200
      t235 = t16 * s26
      t236 = t196 * t74
      t237 = 32.0_dp * t6
      t238 = 170.0_dp * t9
      t239 = t208 * t74
      t240 = 81.0_dp * s26
      t241 = 42.0_dp * s16
      t242 = s16 * t45
      t243 = t125 * t182 + t242
      t169 = t6 * (((((-t165 + t3) * s25 + t175) * s25 + t169) * s25 - 2
     #8.0_dp * t7) * s25 + 11.0_dp * t12)
      t175 = 20.0_dp * t105 * t40
      t244 = t158 * t183 * t37
      t245 = 120.0_dp * t235
      t246 = 90.0_dp * s26
      t247 = -216.0_dp
      t248 = 28.0_dp * s25
      t249 = t183 * s26
      t250 = s16 * (((((-t165 + t127) * s25 + t199) * s25 + 59.0_dp * t9
     #) * s25 - 67.0_dp * t7) * s25 + 21.0_dp * t12)
      t251 = 30.0_dp * t105
      t252 = t251 * t37
      t253 = t25 * t113
      t254 = 17.0_dp * s25
      t255 = 78.0_dp * t9
      t256 = ((((-t67 + t168) * s25 - t237) * s25 + t255) * s25 - 62.0_d
     #p * t7) * s25
      t257 = 17.0_dp * t12
      t258 = 20.0_dp * t125
      t259 = t136 * t17
      t28 = ((t109 * (-t17 * (30.0_dp * t235 + t198 + t199) - 27.0_dp * 
     #t109) + 5.0_dp * t135 - 20.0_dp * t137 - s34 * (t104 * (-t230 + t6
     #7 + t149) - t109 * (t128 + t210 + t209)) + t28 * t41) * s56 + t110
     # * (-t5 * (t97 + t254 + t173) + 27.0_dp * s34) - ((-((228.0_dp * s
     #16 + 208.0_dp * s25 + t246) * s26 + t115 + t207) * s34 - (-t245 - 
     #t234) * s26 - t194 - t195) * s34 + t104 * ((t213 + t228 - t129) * 
     #s26 + t225 + t226)) * s34 + t256 + t257 + t258) * s56 + t74 * t104
     # * (-(-t40 * t55 + t223 + t224 + s26 * (t229 * s26 / 2.0_dp + t227
     #)) * s34 + t249 * t104) + t109 * ((((t240 + t248 + t197) * s34 + (
     #s16 * t247 - 162.0_dp * s26 - t147) * s26 + t154 + t160) * s34 + (
     #(t211 * t64 + t246) * s26 + t239) * s26 + t204 - t205) * s34 + ((-
     #t108 * t200 - t245) * s26 + t236) * s26 + t189 + t79) + t45 + t250
     # + t252 - t253
      t26 = t28 * s56 + t17 * t243 - ((((-((56.0_dp * s25 + t240 + t241)
     # * s26 - t222) * s34 - ((t108 * t219 / 4.0_dp - t99) * s26 + t17 *
     # t217) * s26 - 20.0_dp * t136 * t1 - 115.0_dp * t9) * s34 - s16 * 
     #(((s16 * t162 - t18) * s25 - t237) * s25 - t238) - (((t211 * t5 + 
     #t128) * s26 + t239) * s26 + t17 * t206) * s26) * s34 - s16 * ((((-
     #t231 + t139) * s25 - t232) * s25 - t233) * s25 + 140.0_dp * t7) - 
     #(((-60.0_dp * t235 - t234) * s26 + t236) * s26 + t17 * t190) * s26
     #) * s34 + t104 * (s16 * (((t26 + t98) * s25 + 58.0_dp * t6) * s25 
     #+ 61.0_dp * t9) + (((t229 - t230) * s26 + t227 * t74) * s26 + t64 
     #* (t224 + t223)) * s26)) * s34 + t169 + t175 + t244 + t64 * (s16 -
     # t118) * t113
      t28 = t61 * (s16 * (((t185 + t68) * s25 + t6 * t85) * s25 + t94) +
     # t40 * (t17 * (t81 + t168) - t118) + s26 * (((t184 + t148) * s25 +
     # 50.0_dp * t6) * s26 + ((35.0_dp * s16 + t149) * s25 + t138 * t6) 
     #* s25 + t60)) * s34
      t45 = 7.0_dp * t6
      t60 = (t97 - t149) * s25 - t45
      t79 = t17 * t7
      t97 = t41 * t60 + t79
      t98 = 5.0_dp * t105 * t38
      t26 = t7 * t97 + t74 * (s34 * (t109 * t74 - t5 * t53 - t259 + t47)
     # * t178 - t156) + t26 * s56 + t109 * (t85 * t201 * t7 + 30.0_dp * 
     #t11 + t37 * (s26 * t196 + t190) - t108 * t16 * t39 - t17 * t200 * 
     #t38 + 8.0_dp * t170 - t140 * t43 + t20 * ((((-t68 + t27) * s25 - t
     #186) * s25 - t140) * s25 + 35.0_dp * t7) * s26) + t110 * (t9 * (37
     #.0_dp * s25 + t221) - 27.0_dp * t38 + 5.0_dp * s16 * ((t5 * t201 +
     # t215) * s26 + t220) + t37 * (s26 * t219 + t217)) + t111 * (-40.0_
     #dp * t12 + t25 * t39 + ((s26 * t211 + t208) * s26 + t206) * t37 - 
     #t213 * t7 + t5 * t6 * t43 - t214 * (((s25 - t87) * s25 + t202) * s
     #25 + 85.0_dp * t9)) + t113 * ((-t155 + t69) * s26 + t83) + t114 * 
     #(((7.0_dp * t146 + t173) * s26 - t222) * s26 - t6 * (t185 + t152))
     # - t28 + t56 * t180 + t105 * t183 * t40 + t105 * t182 * t37 + t164
     # * t16 + t98
      t27 = s16 + s26 + s56
      t28 = 57.0_dp * t9
      t47 = 18.0_dp * t6
      t53 = s16 + t18
      t83 = t53 * s26
      t85 = t6 * t212
      t87 = 42.0_dp * t9
      t97 = t74 * s25
      t2 = ((((-t2 - t97) * s25 - t50) * s25 - t255) * s25 - 105.0_dp * 
     #t7) * s25 + 277.0_dp * t12
      t50 = 94.0_dp
      t99 = 90.0_dp * t6
      t115 = t50 * t9
      t128 = (((-t254 - t228) * s25 - t99) * s25 - t115) * s25
      t139 = 141.0_dp * t7
      t140 = t128 + t139
      t147 = 45.0_dp * t6
      t155 = ((-t132 - t216) * s25 - t147) * s25
      t156 = t94 + t155
      t160 = t81 + t77
      t189 = t235 * t160
      t190 = s16 * (((((-s16 - t168) * s25 - t85) * s25 - t87) * s25 - 1
     #8.0_dp * t7) * s25 + 195.0_dp * t12)
      t194 = 42.0_dp * t12
      t195 = 19.0_dp * s25
      t132 = (-t132 - t195) * s25
      t196 = 148.0_dp * t9
      t198 = 295.0_dp * t7
      t199 = 109.0_dp * t6
      t162 = (((t32 + t218) * s25 + t199) * s25 - t162 * t9) * s25
      t200 = -335.0_dp * t7
      t204 = t200 + t162
      t84 = ((t84 + t213) * s25 + t99) * s25
      t99 = 104.0_dp * t9
      t205 = -t99 + t84
      t161 = s16 * t161
      t206 = (t70 + t161) * s25 + t202
      t207 = t67 + t3
      t208 = 80.0_dp * t9
      t145 = ((-t193 - t192) * s25 + t154) * s25 + t145
      t154 = (-100.0_dp * s16 - 62.0_dp * s25) * s25 + t153
      t192 = 25.0_dp * s16
      t209 = -t184 - t192
      t131 = s16 * (((-t131 + t221) * s25 + 192.0_dp * t6) * s25 + 265.0
     #_dp * t9)
      t211 = 61.0_dp * t6
      t217 = 30.0_dp * s25 * t16
      t219 = -t211 + t217
      t222 = t103 + t148
      t223 = 42.0_dp * t7
      t225 = t47 * t41
      t165 = s16 * (14.0_dp * s25 + t165)
      t22 = -t25 * t40 + t6 * (t29 + t22) + s26 * (-s26 * (s16 + t32) + 
     #t165)
      t33 = ((s16 * t191 - t168) * s25 + t33) * s25 + 240.0_dp * t9
      t226 = (-t13 + t69) * s25 + 17.0_dp * t6
      t227 = t226 * t5
      t160 = t16 * t160
      t228 = t160 * t126
      t229 = t156 * t5
      t234 = t74 * t205
      t235 = t206 * t5
      t236 = t74 * t154
      t239 = t64 * s26
      t240 = 64.0_dp * s16
      t245 = t53 * t37
      t246 = 52.0_dp * s16
      t247 = t189 * t55
      t255 = t55 * t207
      t260 = -124.0_dp
      t261 = 40.0_dp * s26
      t262 = 30.0_dp * t9
      t49 = (t114 * t25 - 20.0_dp * t49 - ((((t230 + t184 + t192) * s34 
     #- s16 * (35.0_dp * s26 + t76) - (t261 + t70 + t161) * s25) * s34 -
     # t155 + t228 - t94) * s34 + t104 * (s16 * (-t126 + t181) + (-s26 *
     # t55 - t13 + t69) * s25)) * s34 + 5.0_dp * s16 * t135 + t262 * t41
     #) * s56 + s16 * (t256 + t257 + t258) - t253 - ((((-(t174 + t103 + 
     #t148) * s34 - (s25 * t260 - t129) * s26 - t153 + 62.0_dp * t41 + 1
     #00.0_dp * (s25 + s26) * s16) * s34 - (s26 * t255 + t235) * s26 - t
     #84 + t99) * s34 - (t229 - t247) * s26 - t128 - t139) * s34 + t104 
     #* (((t181 - t18) * s25 + t237) * s25 + (-t83 * t55 + t227) * s26 +
     # t208)) * s34
      t49 = t49 * s56 + s16 * (t250 + t252) + t74 * (t56 * t183 + t112) 
     #- ((((((s16 + t32 + t173) * s34 - (114.0_dp * s25 + t102 + t218) *
     # s26 + t211 - t217) * s34 - s16 * ((-81.0_dp * s25 + t101) * s25 +
     # 215.0_dp * t6) - (t239 * t209 + t236) * s26 + 30.0_dp * t40 + 30.
     #0_dp * t43) * s34 - s26 * (t255 * t37 + t234) - t64 * t206 * t37 -
     # t200 - t162) * s34 - s16 * ((((-t32 - t246) * s25 - t187) * s25 +
     # t188) * s25 + 277.0_dp * t7) - t37 * (t156 * t64 - t247) + t74 * 
     #(-s26 * t140 + t46)) * s34 + t104 * ((-t245 * t55 + t33) * s26 + t
     #6 * ((t213 + t240) * s25 + 119.0_dp * t6) + t64 * (s16 * t43 + t22
     #6 * t37))) * s34 + t242
      t33 = t49 * s56 + s16 * (t169 + t175 + t244) + t17 * s16 * t243 - 
     #((((((-(t239 - t81) * s34 - (-t193 - t173 - t142) * s26 - t165) * 
     #s34 - s16 * ((-t88 - t168) * s25 - 141.0_dp * t6) - t138 * t40 - s
     #26 * (t74 * t222 * s26 + t17 * t219)) * s34 - t40 * (t5 * t209 - t
     #230) - t131 - s26 * (s26 * t236 + t17 * t145)) * s34 - s16 * (((-t
     #132 - t186) * s25 - t196) * s25 - t198) - t40 * (t126 * t207 + t23
     #5) - s26 * (s26 * t234 + t17 * t204)) * s34 - (((t229 - t228) * s2
     #6 + t140 * t74) * s26 + t17 * t2) * s26 - t190) * s34 + t104 * (t1
     #50 * t7 + t202 * t43 + t224 * t41 + t227 * t40 + t33 * t37 - 5.0_d
     #p * t53 * t38 + 71.0_dp * t12 + t214 * (((s16 * t212 + t168) * s25
     # + 64.0_dp * t6) * s25 + 119.0_dp * t9))) * s34
      t49 = t178 * (t17 * t43 - t74 * (t111 + t220) - (-s34 * t207 + t16
     #0) * s34 + t9) - t61 * (t6 * (((t69 + t3) * s25 + t134) * s25 + 14
     #.0_dp * t9) + (s16 * (((t101 + t3) * s25 + t232) * s25 + t28) + ((
     #(-t13 + t68) * s25 + t47 - t83) * s26 + ((-t18 + t197) * s25 + t13
     #4) * s25 + 62.0_dp * t9) * s26) * s26)
      t70 = t105 * t37
      t83 = t7 * t46
      t2 = t49 * s34 + t33 * s56 + t109 * (t6 * (((((-t97 - t142) * s25 
     #- t65) * s25 - t233) * s25 + 25.0_dp * t7) * s25 + t194) + ((((t15
     #6 - t189) * s26 + t140) * s26 + t2) * s26 + t190) * s26) + t110 * 
     #(-t39 * t74 + t6 * (((47.0_dp * s16 + t3) * s25 + t86) * s25 + 70.
     #0_dp * t9) + (((s26 * t209 + t154) * s26 + t145) * s26 + t131) * s
     #26) + t111 * (t6 * (((t41 * t5 - t171) * s25 - t208) * s25 - 70.0_
     #dp * t7) - (s16 * (((t132 + t186) * s25 + t196) * s25 + t198) + ((
     #(-s26 * t207 - t206) * s26 - t205) * s26 - t204) * s26) * s26) + t
     #113 * t22 + t114 * (t25 * t38 + t37 * (s26 * t222 + t219) - t223 -
     # t225 - 52.0_dp * t159 - t118 * s16 * ((t221 + t18) * s25 + 47.0_d
     #p * t6)) + t54 + t125 * t6 * t180 + t70 * (t182 + t249) * s16 - (-
     #s26 + t142) * (s16 + t118) * t112 + t83
      t22 = s26 * t5
      t33 = -t6 + t41
      t20 = (t20 + t168) * s25 + t6
      t49 = s16 * t55
      t84 = t74 * t6
      t86 = (t49 + t168) * s25 + t84
      t53 = t17 * t53
      t88 = s16 * t74
      t99 = s26 * t108
      t101 = t25 * s16 * t24
      t102 = (t13 + t29) * t41
      t103 = t118 * t20
      t126 = t24 ** 2.0_dp
      t128 = t97 * t106
      t131 = s26 * s34
      t84 = ((-t53 * t110 - s34 * ((-s34 * t146 * t1 + t106) * s25 - t10
     #9 * (t109 + t20))) * s56 - t128 * t24 * s34 + ((((t24 * t74 - s34 
     #+ t13) * s34 - (t49 + t99) * s25 - t64 * (s16 * s26 + t41) - t84) 
     #* s34 + t102 + t103 + t9) * s34 + t201 * ((s25 + t99 + t151) * s25
     # + t101)) * t109) * s56 + (-t51 * t114 - t128 * t126) * s34 + ((((
     #(t69 + t118 + t3) * s26 + t259) * s34 - t17 * t86 * s26 - t64 * (t
     #245 + t220) - t48) * s34 + t64 * t141 * t41 + (t17 * (t9 + t102) +
     # t103) * s26 - 7.0_dp * t159) * s34 + t201 * t24 * ((t49 + t99 + t
     #18) * s25 + t101)) * t109
      t99 = s16 + s26 - s34 + s56
      t65 = (s25 + t49) * s25 + t65
      t101 = 31.0_dp * s16
      t102 = ((38.0_dp * s16 + s25) * s25 + t85) * s25 + t36
      t103 = 56.0_dp * s16
      t128 = (t103 + t195) * s25 + t108 * t6
      t132 = 78.0_dp * s16
      t139 = 23.0_dp * t7
      t140 = t185 + t151
      t145 = ((((-80.0_dp * s16 + 63.0_dp * s25) * s25 - 234.0_dp * t6) 
     #* s25 - 357.0_dp * t9) * s25 + 389.0_dp * t7) * s25 + 75.0_dp * t1
     #2
      t148 = 99.0_dp
      t150 = ((((s16 * t148 + t95) * s25 - 254.0_dp * t6) * s25 - 289.0_
     #dp * t9) * s25 + 84.0_dp * t7) * s25 + 220.0_dp * t12
      t34 = (((t75 - t34) * s25 - 126.0_dp * t6) * s25 - 121.0_dp * t9) 
     #* s25
      t75 = 228.0_dp * t7
      t151 = t34 + t75
      t153 = 54.0_dp * s25
      t154 = 96.0_dp * t9
      t155 = ((-t231 + t153) * s25 - t19) * s25 + t154
      t156 = 11.0_dp * t41 + t157
      t157 = s16 * t145
      t160 = (-313.0_dp * s16 + t35) * s25
      t161 = 301.0_dp * t6
      t162 = 386.0_dp * t9
      t165 = 135.0_dp * t7
      t169 = ((t160 - t161) * s25 + t162) * s25 + t165
      t171 = 40.0_dp * s25
      t173 = s25 * (t41 * (-t171 + t210) - t262)
      t175 = -320.0_dp * t6 * t33
      t178 = t173 - t175
      t180 = 84.0_dp * s25
      t107 = (s16 * t107 - t180) * s25
      t182 = 59.0_dp * t6
      t183 = 252.0_dp * t9
      t184 = (t182 + t107) * s25 - t183
      t185 = 42.0_dp * s25
      t186 = (t81 + t185) * s25 + 72.0_dp * t6
      t187 = -t13 - t81
      t188 = -203.0_dp
      t189 = 57.0_dp * s25
      t188 = ((s16 * t188 - t189) * s25 + 174.0_dp * t6) * s25 + 145.0_d
     #p * t9
      t190 = 110.0_dp * s16
      t134 = ((t171 - t190) * s25 - t134) * s25 + 260.0_dp * t9
      t171 = (s16 + t254) * s25 + t19
      t192 = t81 + t97
      t196 = s16 * t188
      t198 = t9 * t43
      t200 = t5 * t192
      t201 = t7 * s25
      t181 = (-t189 + t181) * s25 + t203 * t6
      t189 = (t78 - t95) * s25 - 112.0_dp * t6
      t202 = t81 + t18
      t203 = s16 * t181
      t204 = t64 * t202
      t58 = -t51 * t48 * (s16 * t58 - t97) - 15.0_dp * t73 - t201 * t108
     # - t204 * t38 - t37 * (-s26 * t189 + t203)
      t97 = 33.0_dp * s16
      t205 = s25 + t81
      t206 = s16 * (t97 - t77)
      t207 = ((t22 * t205 + t206) * s26 + 16.0_dp * t48) * s26 + t159 * 
     #t17
      t199 = (((130.0_dp * s16 + 25.0_dp * s25) * s25 + t199) * s25 + 23
     #5.0_dp * t9) * s25 + t139
      t209 = ((92.0_dp * s16 + t82) * s25 + 91.0_dp * t6) * s25
      t210 = 102.0_dp * t9
      t213 = t210 + t209
      t217 = (t193 + t101) * s25 + 56.0_dp * t6
      t124 = ((((102.0_dp * s16 + t13) * s25 + t124) * s25 + 211.0_dp * 
     #t9) * s25 + t163) * t74
      t163 = 5.0_dp * t217
      t219 = t5 * t151
      t222 = t155 / 2.0_dp
      t224 = t222 * t55
      t226 = t150 * t74
      t227 = -t187
      t228 = -t178 * t74
      t229 = -t184 * t5
      t232 = 5.0_dp * t186
      t233 = t134 * t64 / 2.0_dp
      t234 = 84.0_dp * t6
      t189 = t74 * t189
      t235 = s16 * t44
      t206 = (t125 * t102 + t235) * s16 + t113 * ((t239 * t205 + t206) *
     # s26 + 8.0_dp * t48)
      t217 = t217 * t55
      t236 = 15.0_dp * s26 * t140
      t237 = 15.0_dp * s26 * t156
      t55 = t186 * t55
      t242 = -55.0_dp
      t243 = 66.0_dp * s26
      t19 = (-t108 * t137 + t135 * t74 + t225 - ((-s34 * t187 - t156) * 
     #s34 + t104 * t140) * s34) * s56 + s25 * ((t212 * t9 + t43 * t64) *
     # s25 - t223) + t108 * s16 * t172 - (((s16 * (t218 + t129) + (24.0_
     #dp * s26 + t81 + t185) * s25) * s34 - ((-t231 + t153 + t243) * s25
     # - t19) * s25 - t6 * (96.0_dp * s16 + 78.0_dp * s26)) * s34 + t104
     # * (s16 * (t103 + t243) + (t96 + t193 + t101) * s25)) * s34 + 18.0
     #_dp * t125 + t200 * t110
      t19 = t19 * s56 + s25 * ((((-t234 + t90) * s25 + t154) * s25 - t25
     # * t7) * s25 - t194) - (((-(s16 * (138.0_dp * s16 + 100.0_dp * s26
     #) + (s25 * t191 + t88 + t96) * s25) * s34 + (-t182 - t107) * s25 +
     # (-15.0_dp * s26 * t187 + t232) * s26 + t183) * s34 - (t224 + t237
     #) * s26 - t34 - t75) * s34 + t104 * ((t163 + t236) * s26 + t209 + 
     #t210)) * s34 + 18.0_dp * s16 * (t12 + t46) + 30.0_dp * t125 * t144
     # + 45.0_dp * t70 - t204 * t114
      t19 = t19 * s56 + t108 * (t125 * t65 + t10) + t110 * (t17 * (s16 *
     # ((s25 * t242 - t152) * s25 + 130.0_dp * t6) + 20.0_dp * t192 * t3
     #7 + 20.0_dp * t43 + t239 * t171) + t205 * t5 * t109) + 60.0_dp * t
     #70 * (t144 + s26) - (((t109 * (s16 * (112.0_dp * s16 + 120.0_dp * 
     #s26) + (s26 * t212 - t78 + t95) * s25) - ((-20.0_dp * t227 * s26 -
     # t55) * s26 - t229) * s26 - t173 + t175) * s34 - s16 * ((((s25 * t
     #148 - t121) * s25 - 289.0_dp * t6) * s25 + 84.0_dp * t9) * s25 + 2
     #20.0_dp * t7) - 20.0_dp * t222 * t37 - 20.0_dp * t46 - s26 * (20.0
     #_dp * t37 * t156 + t219)) * s34 + t104 * (s16 * (((t190 + t100) * 
     #s25 + 211.0_dp * t6) * s25 + t208) + t37 * (20.0_dp * s26 * t140 +
     # t217) + t5 * (s26 * t213 + t42))) * s34 + t136 * (((((-20.0_dp * 
     #s16 + t195) * s25 - 98.0_dp * t6) * s25 + 212.0_dp * t9) * s25 - 1
     #33.0_dp * t7) * s25 + 8.0_dp * t12)
      t19 = t19 * s56 + s16 * (-5.0_dp * t112 + t136 * (((((s25 * t130 -
     # t31) * s25 + 35.0_dp * t6) * s25 + 125.0_dp * t9) * s25 - 136.0_d
     #p * t7) * s25 + t12 * t138)) + t74 * s16 * (t125 * t128 + t10) - (
     #((((-((t96 + t97 - t77) * s16 + s26 * t149) * s34 + (t174 * t202 -
     # t189) * s26 + t203) * s34 - ((t261 * t192 + 18.0_dp * t171) * s26
     # + t233) * s26 - t196) * s34 - s16 * (((-t160 + t161) * s25 - t162
     #) * s25 - t165) - (((-t230 * t227 - t55) * s26 + t64 * t184) * s26
     # - t228) * s26) * s34 - (((20.0_dp * t222 + t237) * s26 + t151 * t
     #64) * s26 + t226) * s26 - t157) * s34 + t104 * (s16 * t199 + (((t2
     #17 + t236) * s26 + t213 * t64) * s26 + t124) * s26)) * s34 + t235 
     #+ 18.0_dp * t70 * t65 + 60.0_dp * t105 * t144 * t40 + 45.0_dp * t1
     #05 * t38
      t10 = t19 * s56 + t17 * t206 - ((((-t108 * t171 * t40 - (-t214 * t
     #181 + t189 * t37 - 24.0_dp * t202 * t40 - 98.0_dp * t159 + t176) *
     # s34 - t214 * t188 - 20.0_dp * t192 * t38 - t233 * t37 + t234 * t4
     #3 + t73 * t74 - 250.0_dp * t201) * s34 + t167 * t227 + t214 * t169
     # + t228 * t37 + t229 * t40 + t232 * t38 - 54.0_dp * t170 + 340.0_d
     #p * t179 - 188.0_dp * t198 - 14.0_dp * t52) * s34 - t143 * t11 - t
     #214 * t145 + 84.0_dp * t12 * t41 + 138.0_dp * t7 * t43 + 128.0_dp 
     #* t9 * t42 - t226 * t37 - t59 * t46 - t224 * t38 - t40 * (t37 * t1
     #56 * t64 + t219)) * s34 + t104 * (t39 * t64 * t140 + t5 * t213 * t
     #40 + t124 * t37 + t170 * t138 + t163 * t38 + t214 * t199 + 106.0_d
     #p * t179 + 54.0_dp * t198 + 104.0_dp * t52)) * s34 + 18.0_dp * t21
     # + 18.0_dp * t23 - t49 * s26 * t112 + t116 * t65 * t40 + t251 * t1
     #44 * t38 - 14.0_dp * t11 * t43 - 41.0_dp * t10 * t41 + 23.0_dp * t
     #54 - t83 * t50 + 106.0_dp * t12 * t42 + t158 * s16 * t128 * t37
      t12 = t61 * (t48 * (((t32 + t142) * s25 + t15) * s25 + t94) + ((((
     #s26 * t140 + (t193 + t197) * s25 + t147) * s26 + ((62.0_dp * s16 +
     # t82) * s25 + 70.0_dp * t6) * s25 + t28) * s26 + (((t132 + t13) * 
     #s25 + t85) * s25 + t30) * s25 + t139) * s26 + t136 * (((t246 + t21
     #6) * s25 + t211) * s25 + t115)) * s26) * s34
      t15 = t125 * t48 * ((t101 + t18) * s25 + t47)
      t8 = t10 * s56 + t109 * (t48 * (t6 * (((-t132 - t57) * s25 + 66.0_
     #dp * t6) * s25 + t262) + t51 * ((((-t240 + t216) * s25 - t166) * s
     #25 - t87) * s25 + 130.0_dp * t7)) + ((((s26 * t156 + t155) * s26 +
     # t151) * s26 + t150) * s26 + t157) * t37) + t110 * (t74 * (t171 * 
     #t38 - t198) + s26 * ((s26 * t134 + t196) * s26 + t48 * ((-t180 - t
     #88) * s25 + 250.0_dp * t6)) + t200 * t39 + t201 * (61.0_dp * s25 +
     # t221)) + t111 * (t48 * (s16 * (((t231 + t193) * s25 - t50 * t6) *
     # s25 - 40.0_dp * t9) - t51 * (((-s16 * t50 - t35) * s25 - t45) * s
     #25 + t238)) - (s16 * t169 + (((-s26 * t187 + t186) * s26 - t184) *
     # s26 - t178) * s26) * t37) + t113 * t207 + t114 * t58 + t41 * t8 +
     # t44 * t9 - t12 + t70 * (s26 * t128 + t102) * s16 + t15
      t10 = 0.1e1_dp / s16
      t12 = 0.1e1_dp / t66 ** 2.0_dp
      t15 = 0.1e1_dp / t63 ** 2.0_dp
      t19 = 0.1e1_dp / t62 ** 2.0_dp
      t2 = t26 * MB1001(0) + ((t74 * t105 * t38 * (t37 + t65) - 5.0_dp *
     # s16 * (t112 * t37 - t54) + t8 + t167 * t105 * t144 + t179 * (t41 
     #* ((t241 - t248) * s25 - t215) + t79)) * MB1110(1) - t17 * (t105 *
     # t4 * t38 + t44 * t6) - t5 * t11 * (t9 + t43) - t14 + 7.0_dp * t21
     # - 7.0_dp * t23 - t52 * (t41 * (-t32 + t29) - t36) + t119 * t37 - 
     #t27 * (s16 * (t7 * ((s25 * t60 - t36) * s25 + t79) + t98) + t2) * 
     #MB1011(1)) * t10
      t4 = t27 ** 2.0_dp
      t1 = t10 * (t12 * (-(s34 * ((-t106 * t24 * t126 + s34 * t126 * ((t
     #22 + s25) * s25 + t74 * (s25 + s26 + s16) * s16) * t1) * s25 + t10
     #9 * (s34 * (t131 * (-t131 + (s26 + t88 + t13) * s26 + t259) - t53 
     #* t40 - t37 * t86 - t123 + t159 - t136 * (s16 + t168) * s26) - t24
     # * (t17 * t220 * t1 + (-s26 * t20 - t13 * (t33 + t259)) * s26 + t1
     #59 * t74))) + t84 * s56) * MB1101(0) * t19 + t27 * t4 * t99 * MB11
     #11D2(0)) - t99 ** 2.0_dp * MB1110(0) * t15)
      result = t12 * t19 * (-(t17 * t9 * t172 - t74 * (t37 * (t114 - t12
     #5) - (-t137 * t5 + (-t104 * t146 + t109 * (-t120 + s34)) * s34 + t
     #42 + t7 + t176 + t17 * t91 * t109) * s56 * t177) + t133 * s56 + t1
     #09 * (t6 * ((5.0_dp * t6 + t93) * s25 + t94) + s26 * (t89 * s26 + 
     #t92) + t64 * t91 * t40) + t110 * (((t69 + t127 + t118) * s26 + t14
     #1 * t64) * s26 + t17 * t9) - t111 * (t74 * (t120 * t40 + t123) + t
     #9 * (t168 + t68) + s26 * (s26 * t117 + t122)) - t61 * (8.0_dp * t6
     # * t16 + s26 * (t80 * s26 + (t77 + t76) * s25 + t71) + t81 * t41) 
     #* s34 + t56 * (t67 + t18) + t105 * t72 * t37 + t164 + t179 * (-t67
     # + t3)) * I300s16(1) + t15 * t2) / 12.0_dp - t1 / 6.0_dp

           intHs16s25s26s34s56x1411D6eps0 = result
       end function intHs16s25s26s34s56x1411D6eps0

       function intHs16s25s26s34s56x1411D6eps1()
           implicit none
           complex(dp) :: intHs16s25s26s34s56x1411D6eps1
           complex(dp) ::  t1,t2

           complex(dp) :: result

      t1 = s16 + s26 - s34 + s56
      t2 = s16 + s25 + s26 - s34 + s56
      t2 = 0.1e1_dp / t2 ** 2.0_dp
      result = 0.1e1_dp / s16 * t2 * (-s25 / 12.0_dp - t1 / 4.0_dp + t1 
     #** 2.0_dp * MB1110(1) / 6.0_dp)

           intHs16s25s26s34s56x1411D6eps1 = result
       end function intHs16s25s26s34s56x1411D6eps1

       function ampNonresonantHeavyImC4MM()
           implicit none
           complex(dp) :: ampNonresonantHeavyImC4MM

           complex(dp) :: result

      result = 0.0_dp

           ampNonresonantHeavyImC4MM = result
       end function ampNonresonantHeavyImC4MM

       function ampNonresonantHeavyImC4MP()
           implicit none
           complex(dp) :: ampNonresonantHeavyImC4MP
           complex(dp) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           complex(dp) ::  t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30
           complex(dp) ::  t31,t32,t33,t34,t35,t36,t37,t4,t5,t6,t7,t8
           complex(dp) ::  t9

           complex(dp) :: result

      t1 = intHLs250000x0121D2eps0()
      t2 = intHLs250000x0121D2eps1() * epinv
      t3 = (t2 + t1) * s25 + epinv2
      t4 = gb ** 2.0_dp
      t5 = propZ25 * s25
      t6 = t5 * gw ** 2.0_dp
      t7 = s25 + s26 + s56
      t8 = s16 + s26 - s34 + s56
      t9 = 2.0_dp * propW16
      t10 = 6.0_dp * propW16
      t11 = 3.0_dp * t6
      t12 = t4 * (t10 * t7 + t5 + 2.0_dp) + t11 * (t9 * t7 + 1.0_dp)
      t13 = intHLs250000x0112D2eps1()
      t14 = intHLs250000x0113D4eps1()
      t15 = intHLs250000x0122D4eps1()
      t16 = 2.0_dp * t14
      t17 = (-t16 + t2 + t1 + t13 - t15) * s25 + epinv2
      t9 = t9 * t8 - 1.0_dp
      t18 = s25 + s16 - s34
      t19 = s26 + s56
      t20 = t18 + 2.0_dp * t19
      t21 = s34 ** 2.0_dp
      t22 = s16 ** 2.0_dp
      t23 = s25 ** 2.0_dp + t21 + t22
      t24 = s56 ** 2.0_dp
      t21 = s25 * (s16 * s25 - s25 * s34 + t21 + t22)
      t22 = s12 * t7
      t25 = s26 * propW16
      t26 = t25 * ((s25 + s56) * s16 + (-s25 - s56) * s34 + s56 * (s25 +
     # s26))
      t27 = s26 ** 2.0_dp
      t28 = s16 * s34
      t28 = (((-s16 + s34) * s56 + t28) * s25 - s26 * t27 + (t28 - t24) 
     #* s26 + t28 * s56) * propW16 - s45 * t19 - s46 * t19
      t29 = t27 * propW16 * (s25 + s16 - s34)
      t10 = t10 * t8 + t5 - 4.0_dp
      t30 = -4.0_dp * s25 + 2.0_dp * s34 - 6.0_dp * t19 - s16 * (t5 + 2.
     #0_dp) + t5 * (s25 + s34)
      t31 = -1.0_dp + t5
      t32 = s25 * (-1.0_dp + (s25 + s45 + s46) * propZ25)
      t33 = t31 * s26
      t31 = t31 * s56
      t34 = -s45 - s46
      t35 = t4 * t10
      t36 = t11 * t9 + t35
      t37 = (t2 + t1 - t15) * t4
      t10 = t12 * t8 * t3 * (struc15MP + struc41MP + struc57MP) - t17 * 
     #(t4 * (12.0_dp * t26 + 9.0_dp * t29 - 6.0_dp * t28 + 3.0_dp * (pro
     #pW16 * t23 + s25) * s56 + 3.0_dp * t24 * (propW16 * t18 + 1.0_dp) 
     #+ 3.0_dp * t25 * t23 + 3.0_dp * t21 * propW16 + s13 * t30 + (t31 +
     # t32 + t33) * s16 + (s26 * (propZ25 * t19 - 1.0_dp) + t5 * (s26 - 
     #s45 - s46)) * s25 + (-t31 - t32 - t33) * s34 - t27 + t22 * t10 - 2
     #.0_dp * s16 * t34 + 2.0_dp * s26 * s56 + 2.0_dp * s34 * t34 - 4.0_
     #dp * s25 * t34) + t11 * (4.0_dp * t26 + (s26 * t23 + s56 * t23 + t
     #18 * t24 + t21) * propW16 - s13 * t20 + s45 * t18 + s46 * t18 + s5
     #6 * (s25 + s26) + t24 + t22 * t9 - 2.0_dp * t28 + 3.0_dp * t29)) *
     # struc1MP
      t18 = 0.1e1_dp / t8
      t19 = 0.1e1_dp / t7
      result = -0.2e1_dp / 0.3e1_dp * im * propW34 * (-2.0_dp * t17 * (-
     #struc2MP * (t11 * t20 - t4 * t30) + t7 * t36 * (struc3MP + struc8M
     #P)) + 2.0_dp * t7 * (epinv2 * (t4 * (-t5 + 4.0_dp) + t11) + s25 * 
     #(6.0_dp * t6 * t14 * t9 + 6.0_dp * (t4 * t8 + t6 * t8) * t15 * pro
     #pW16 - t13 * t36 - t37 * t5 + 3.0_dp * t6 * (t2 + t1 - t15) + 4.0_
     #dp * t37 + t35 * t16)) * struc5MP - 2.0_dp * (-t16 + t13 - t15) * 
     #s25 * t8 * t12 * struc7MP + t10 + 6.0_dp * t8 * t7 * (t6 + t4) * t
     #3 * propW16 * (struc14MP + struc40MP + struc56MP)) / ecossin ** 2.
     #0_dp / s25 * t19 * t18

           ampNonresonantHeavyImC4MP = result
       end function ampNonresonantHeavyImC4MP

       function ampNonresonantHeavyImC4PM()
           implicit none
           complex(dp) :: ampNonresonantHeavyImC4PM
           complex(dp) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           complex(dp) ::  t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30
           complex(dp) ::  t31,t32,t33,t34,t35,t36,t37,t38,t39,t4,t40,t41
           complex(dp) ::  t42,t43,t5,t6,t7,t8,t9

           complex(dp) :: result

      t1 = s25 + s26
      t2 = s25 + s26 + s56
      t3 = 2.0_dp
      t4 = s25 + s16 - s34
      t5 = s16 + s26 - s34 + s56
      t6 = s25 + s26 - s34 + s56
      t7 = s12 + s13 + s16
      t8 = t3 * s26
      t9 = t8 + s25 + s16 - s34 + s56
      t10 = propW16 * (t3 * (s26 + s12) + t4) * t2 * t5
      t11 = (-s16 - s25 - s26 + s34 - s56) * s13
      t12 = t1 * s46
      t13 = (s12 * t3 + s25 - t6) * s25
      t14 = (t3 * (s13 - s46 + s12 - s45) - t2 - s56) * s56
      t15 = t1 * s34
      t16 = intHLs250000x0121D2eps0()
      t17 = intHLs250000x0112D2eps1()
      t18 = s13 + s12
      t19 = s45 + s46
      t20 = -t3 * (t19 * t9 + t10 + t11) + (t18 * t3 + s25 + s34 - s56) 
     #* s26 + t13 + t14 - t15
      t21 = intHLs250000x0113D4eps1()
      t22 = intHLs250000x0122D4eps1()
      t18 = s16 + 4.0_dp * t18
      t23 = 6.0_dp
      t24 = 4.0_dp * s12
      t25 = t23 * s13
      t26 = t3 * s13
      t27 = s16 - s34
      t28 = 3.0_dp * s26
      t29 = s25 * t3 + t27 + t28
      t30 = 3.0_dp * s25
      t31 = s25 - s16 + s34
      t32 = s13 * s16
      t33 = t7 * s25
      t34 = s56 ** 2.0_dp
      t1 = t3 * (-t19 * t29 + t32) - 3.0_dp * t10 - 3.0_dp * t34 + (-pro
     #pZ25 * ((s12 + s16 + s25 + s26) * s26 - (s25 - s13 + s26) * s34 + 
     #(s12 + s16 + s26 - s34) * s56 - t19 * t31 - t32 + t33) + t18) * s2
     #5 + (t25 + t24 + s25 + s16 + s26) * s26 - (t26 + t1) * s34 + (t23 
     #* (s13 - s46 - s45) + t24 + t27 - t8 - t30) * s56
      t8 = intHLs250000x0121D2eps1() * epinv
      t19 = (t8 + t16) * s25 + epinv2
      t35 = s25 * (t17 - t22)
      t36 = propZ25 * s25
      t37 = t36 * gw ** 2.0_dp
      t38 = gb ** 2.0_dp
      t39 = t36 + t3
      t40 = 3.0_dp * t37
      t41 = s26 + s56
      t42 = t3 * propW16
      t43 = t23 * propW16
      t4 = ((-t21 * t3 + t16 + t17 - t22 + t8) * s25 + epinv2) * (-(t38 
     #* (-s34 * t3 + t23 * t41 + 4.0_dp * s25 + s16 * t39 - t36 * (s25 +
     # s34)) + t40 * (t3 * t41 + t4)) * struc2PM + t5 * (t38 * (t43 * t2
     # + t3 + t36) + t40 * (t42 * t2 + 1.0_dp)) * struc7PM + (t38 * (t43
     # * t5 + t36 - 4.0_dp) + t40 * (t42 * t5 - 1.0_dp)) * t2 * (struc3P
     #M + struc5PM + struc8PM))
      t17 = 0.1e1_dp / t5
      t2 = 0.1e1_dp / t2
      result = im * propW34 * (t3 * (t38 * t39 + t40) * ((t8 + t16) * s2
     #5 + epinv2) * t5 * (struc27PM + struc43PM - struc57PM) + (t3 * t38
     # * (t19 * (3.0_dp * s26 ** 2.0_dp + t3 * (-s45 * t29 + t32) - 3.0_
     #dp * t10 - 3.0_dp * t34 + (propZ25 * (-(s12 + s25) * s26 + (s25 - 
     #s13) * s34 + s45 * t31 - (s12 + s16 - s34 - s46) * s56 + t12 + t32
     # - t33) + t18) * s25 + (3.0_dp * s16 + s25 + t24 + t25) * s26 - (s
     #25 + t28 + t26) * s34 + (t23 * (s13 - s45) + t27 - t30 - 4.0_dp * 
     #s46 + 4.0_dp * s12) * s56 - 4.0_dp * t12) + t35 * t1) + 3.0_dp * t
     #37 * (t19 * (-t3 * (s45 * t9 + t10 + t11 + t12) + (t3 * t7 + s26 +
     # t6) * s26 + t13 + t14 - t15) + t35 * t20) - t21 * s25 * (t37 * t2
     #3 * t20 + 4.0_dp * t38 * t1)) * struc1PM - 4.0_dp * t4) / ecossin 
     #** 2.0_dp / s25 * t2 * t17 / 3.0_dp

           ampNonresonantHeavyImC4PM = result
       end function ampNonresonantHeavyImC4PM

       function ampNonresonantHeavyImC4PP()
           implicit none
           complex(dp) :: ampNonresonantHeavyImC4PP

           complex(dp) :: result

      result = 0.0_dp

           ampNonresonantHeavyImC4PP = result
       end function ampNonresonantHeavyImC4PP

       function ampNonresonantHeavyImC7MM()
           implicit none
           complex(dp) :: ampNonresonantHeavyImC7MM

           complex(dp) :: result

      result = 0.0_dp

           ampNonresonantHeavyImC7MM = result
       end function ampNonresonantHeavyImC7MM

       function ampNonresonantHeavyImC7MP()
           implicit none
           complex(dp) :: ampNonresonantHeavyImC7MP
           complex(dp) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           complex(dp) ::  t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30
           complex(dp) ::  t31,t32,t33,t34,t35,t36,t37,t38,t39,t4,t40,t41
           complex(dp) ::  t42,t43,t44,t45,t46,t47,t48,t49,t5,t50,t51,t52
           complex(dp) ::  t53,t54,t55,t56,t57,t58,t59,t6,t60,t61,t62,t63
           complex(dp) ::  t64,t65,t66,t67,t68,t69,t7,t70,t71,t72,t73,t74
           complex(dp) ::  t75,t76,t77,t78,t79,t8,t80,t81,t82,t83,t84,t85
           complex(dp) ::  t86,t87,t88,t89,t9,t90,t91,t92,t93

           complex(dp) :: result

      t1 = 2.0_dp
      t2 = t1 * s26
      t3 = -t2 - s25 + s13 - s16 + s34 - s56 + s12 - s45 - s46
      t4 = s26 + s12
      t5 = t1 * t4
      t6 = s25 + s26 + s56
      t7 = s16 + s26 - s34 + s56
      t8 = -s12 + s16 + s25 + s26 - s34 + s56
      t9 = s45 + s46
      t10 = s56 ** 2.0_dp
      t11 = s25 ** 2.0_dp
      t12 = t9 * (t2 + s25 + s16 - s34 + s56)
      t13 = propW16 * (t5 + s25) * t6
      t14 = t13 * t7
      t15 = s26 * (s13 - s16 - s26 + s34 - s56 + s12)
      t16 = (-s16 - s26 + s34) * s13 + (-s13 + t8) * s25 + (-s13 - t3) *
     # s56 + t10 - t11 + t12 + t14 - t15
      t17 = s13 + s12
      t18 = t1 * t17
      t19 = t18 - s16
      t20 = s16 + s25
      t21 = 3.0_dp
      t22 = t1 * s12
      t23 = t21 * s13
      t24 = -t23 - t22 + t20
      t25 = s16 - s34
      t26 = t1 * s25
      t27 = t21 * s26
      t28 = s25 - s13 + s45 + s46
      t29 = -t21 * t28 - 4.0_dp * s26 + t22 - t25
      t30 = s12 + s16 + s25
      t31 = -s16 + s25 + s34
      t32 = s12 + s16 + s26 - s34
      t33 = t9 * t31
      t34 = s13 * s16
      t35 = (-s13 + s25 + s26) * s34
      t36 = (s12 + s13 + s16) * s25
      t37 = (t30 + s26) * s26 + s56 * t32 - t33 - t34 - t35 + t36
      t38 = t9 * (t27 + t26 + t25)
      t39 = propZ25 * t37
      t40 = t14 * t21
      t41 = gw ** 2.0_dp
      t42 = t41 ** 2.0_dp
      t43 = gb ** 2.0_dp
      t44 = t43 ** 2.0_dp
      t45 = t21 * s25
      t46 = t45 * t41 * propZ25
      t47 = intHLs250000x0113D4eps0()
      t48 = s25 + s16 - s12
      t49 = t1 * s13
      t50 = s25 + s16 - s34
      t51 = s16 - s34 - s12
      t52 = t9 * (t2 + t50)
      t53 = (s12 + s13 - s16) * s25
      t28 = t1 * t10 + (-t49 + t48 + s26) * s26 + (t1 * t28 + t27 + t51)
     # * s56 + t14 - t34 - t35 + t52 - t53
      t17 = 4.0_dp * t17
      t54 = s25 + s26
      t55 = s13 - s45 - s46
      t56 = 9.0_dp
      t57 = 4.0_dp * s12
      t58 = 6.0_dp * t55
      t59 = t1 * t55
      t60 = propW16 * (t1 * t32 + s25) * t6
      t61 = t54 * s34
      t62 = 10.0_dp
      t63 = t45 * propZ25
      t64 = propW16 * t6
      t62 = (-8.0_dp * s12 - 12.0_dp * s13) * s26 - t1 * ((-5.0_dp * s16
     # + t17) * s25 + (-s25 * t56 - 14.0_dp * s26 - 5.0_dp * t25 + t57 +
     # t58) * s56) + 18.0_dp * t10 + t7 * (6.0_dp * t64 * (t5 + t31) + t
     #63 * (t60 - t59)) + t62 * ((t20 + s26) * s26 - t61) - 4.0_dp * s13
     # * t25 + 4.0_dp * t38
      t65 = propZ25 * s25
      t66 = t65 * t43
      t67 = t66 * t37
      t68 = t65 * t56 * t42
      t69 = 6.0_dp * s13
      t70 = 4.0_dp * t30 + t69
      t71 = 4.0_dp * t32 + t58
      t60 = t60 * t7
      t72 = (s13 + s25 + s26) * s34
      t30 = (-t49 - t30 - s26) * s26 - (t59 + t32) * s56 - t34 - t36 + t
     #52 + t60 + t72
      t32 = s26 ** 2.0_dp
      t52 = (t1 * t54 + s13) * s34
      t36 = 4.0_dp * t36 + 4.0_dp * t32
      t59 = intHLs250000x0111D2eps0()
      t55 = -t21 * t55 - t22 + t31
      t73 = 4.0_dp * t25
      t74 = t64 * t25 * t7
      t4 = t64 * (6.0_dp * t4 + t45 + t73) * t7
      t23 = (-t23 - t22 - t20 - s26) * s26
      t75 = (t18 + s16) * s25
      t76 = intHLs250000x0112D2eps1()
      t77 = t1 * (s25 - s13 + s26 + s56 + s45 + s46) + t13
      t78 = s16 - s34 + s56
      t79 = -s25 * t19 + (t24 + s26) * s26 - s56 * t29 - t34 - t35 + t38
      t80 = 1.0_dp - epinv
      t81 = t65 * t42
      t82 = intHLs250000x0113D4eps1()
      t83 = -t1 + epinv
      t84 = intHLs250000x0122D4eps1()
      t85 = s26 + s56
      t32 = (-s13 - s12 - s16) * s25 + (-2.0_dp * s13 - s12 - s16 - s25)
     # * s26 + (s13 + t54) * s34 + (-2.0_dp * s13 - s12 - s16 - s26 + s3
     #4 + s45 + s46) * s56 + t12 - t32 - t34 + t60
      t86 = s25 + s26 + s45 + s46 + s56
      t87 = intHLs250000x0111D2eps1()
      t88 = (-s16 - s25 - s26 + s34 - s56) * s13
      t5 = t64 * (t5 + t50) * t7
      t89 = t43 * intHLs250000x0122D4eps0()
      t9 = t76 * (t56 * t81 * t16 * t80 + t43 * (-t41 * (-t21 * (epinv *
     # (t1 * (s13 * t78 + (2.0_dp * s13 + t22 - s25 - s16 - s26 + s34 - 
     #s56) * s25 + (s13 + t18 - t7) * s26 + (-t1 * (s25 - s13 - s12 + s4
     #5 + s46) - t27 - t78) * s56 - t9 * (t27 + t26 + t78) + t11) - 4.0_
     #dp * t10 - t40 - t65 * t37) + t65 * t7 * t77) - 12.0_dp * t10 - 4.
     #0_dp * t79 - 6.0_dp * t14) + t67)) + t84 * (-t68 * t16 - t43 * (t4
     #1 * (epinv * (-s56 * t71 - t1 * (s26 * t70 / 2.0_dp + t34 - t38 - 
     #t52) + t21 * (-t65 * t32 + t60) - t36) + t7 * (t63 * t77 + 6.0_dp 
     #* t13) + 12.0_dp * s26 * t9 + 12.0_dp * (s25 + s45 + s46) * s56 + 
     #12.0_dp * t10 - 4.0_dp * s13 * (t21 * t85 + t25 + t26) + 4.0_dp * 
     #t86 * s16 + 4.0_dp * s26 * t54 - 4.0_dp * t86 * s34 - 8.0_dp * s12
     # * t6 + 8.0_dp * s25 * t9 + 16.0_dp * s26 * s56) + t67 * t80)) + t
     #89 * (t41 * (t1 * (-t52 + t34 - t38) - t21 * (-t65 * t30 + t60) + 
     #s26 * t70 + s56 * t71 + t36) + t67)
      t11 = s25 * intHLs250000x0112D2eps0()
      t14 = t11 * (-t43 * (t1 * ((t24 + s26) * s26 - t34 - t35 + t38) + 
     #6.0_dp * t10 - t1 * t29 * s56 + (-t1 * t19 + t39) * s25 + t40) - t
     #46 * t16)
      t13 = s25 * (t47 * (t68 * t28 + t43 * (t41 * t62 + t67)) + t82 * (
     #t81 * t28 * (epinv * t56 - 18.0_dp) + t43 * (t41 * (epinv * t62 - 
     #t7 * (6.0_dp * t65 * t77 + 12.0_dp * t13) - 24.0_dp * t10 - 8.0_dp
     # * t79) + t67 * t83)))
      t16 = t1 * t25
      t19 = 6.0_dp * t85
      t24 = (propZ25 * t31 - 4.0_dp) * s25 - t16 - t19
      t28 = t1 * t85 + t50
      t29 = t46 * t28
      t35 = -t83
      t36 = propZ25 * t7
      t37 = t1 * t43 * t41
      t40 = (-propZ25 * t28 * t21 + 4.0_dp) * s25 + t16 + t19
      t50 = t66 * t31
      t52 = (6.0_dp * t36 + 8.0_dp) * s25
      t60 = 12.0_dp * t85
      t62 = t21 * epinv
      t67 = t87 * t35
      t35 = t82 * t35
      t70 = t35 - t47
      t71 = 6.0_dp * propW16 * t7
      t77 = -t71 - t65 + 4.0_dp
      t78 = -t1 * propW16 * t7 + 1.0_dp
      t79 = -t78
      t82 = 1.0_dp - t65
      t83 = t71 * t82 + t63 - 4.0_dp
      t85 = t65 + t1
      t86 = -t21 * propW16 * t85 * t7 + 4.0_dp
      t71 = t71 * t85
      t90 = t66 * t80
      t91 = t64 * t1 + 1.0_dp
      t92 = 6.0_dp * t64
      t93 = t92 + t65 + t1
      t63 = t92 * t82 + t1 - t63
      t64 = t64 * t85 + t65
      t82 = 6.0_dp * t64
      t63 = t7 * (t1 * s25 * (t65 * (t42 * t91 * t56 + t44) + t37 * (t21
     # * t64 + t1)) * t70 + t21 * (t43 * t93 + t46 * t91) * t41 * (t67 -
     # t59 + t11) - (t76 * (t56 * t81 * t91 * t80 + t43 * (-t41 * (t62 *
     # t93 - t82 - 4.0_dp) + t66)) - t84 * (t68 * t91 + t43 * (t41 * (ep
     #inv * t63 + t82 + 4.0_dp) + t90)) + t89 * (-t41 * t63 + t66)) * s2
     #5) * struc7MP
      t35 = (-t1 * s25 * (-t47 * (-t37 * t86 + t65 * (t42 * t79 * t56 + 
     #t44)) + t35 * (-t37 * t86 + t65 * (-t42 * t78 * t56 + t44))) + t21
     # * t41 * (t59 * (-t43 * t77 + t46 * t79) + (t67 + t11) * (t43 * t7
     #7 + t46 * t78)) + (t76 * (-t56 * t81 * (epinv * t79 + t78) + t43 *
     # (t41 * (t62 * t77 + t71 - 8.0_dp) + t66)) + t84 * (-t68 * t79 - t
     #43 * (t41 * (epinv * t83 + t71 - 8.0_dp) + t90)) + t89 * (-t41 * t
     #83 + t66)) * s25) * t6 * (struc3MP + struc5MP + struc8MP)
      t11 = (-t1 * s25 * (t65 * (-t42 * t28 * t56 + t31 * t44) - t37 * (
     #(t36 * t21 + 4.0_dp) * s25 + t16 + t19)) * t70 - t21 * t41 * ((t67
     # - t59) * (t24 * t43 - t29) - t11 * (-t24 * t43 + t29)) + (t76 * (
     #-t56 * t81 * t28 * t80 + t43 * (t41 * (-t62 * t24 - t52 - t60 - t7
     #3) + t50)) + t84 * (t68 * t28 - t43 * (t41 * (-epinv * t40 - t52 -
     # t60 - t73) + t50 * t80)) + t89 * (t40 * t41 + t50)) * s25) * stru
     #c2MP - t63 + t35
      t1 = t1 * t11 + (s25 * t9 + t1 * t13 + t21 * (t59 * (t43 * (t1 * (
     #(t55 + s56) * s56 - t34 + t38 + t72 + t23 - t75) + t65 * (t1 * (t7
     #4 - t10) - (t48 + s26) * s26 + (-t27 - t26 - t51) * s56 - t33 - t3
     #4 + t53 + t72) + t4) + t46 * t30) + t87 * (t21 * t65 * t41 * (epin
     #v * t32 - t1 * (t5 + t12 + t88) + (t22 - s26 + s34 - s56) * s25 + 
     #(t18 + s25 + s34 - s56) * s26 + (t1 * (s13 + s12 - s45 - s46) - t6
     # - s56) * s56 - t61) - t43 * (epinv * (-t1 * ((t55 + s56) * s56 - 
     #t34 + t38 + t72 + t23 - t75) + t65 * (-t74 * t1 - (t49 - t8 + s25)
     # * s25 - (s45 + s46 + t3 - s56) * s56 - t15 + t33 - t88) - t4) - t
     #1 * ((-t39 + t17 + s16) * s25 + (t69 + t57 + t20 + s26) * s26 - (t
     #49 + t54) * s34 + (t58 + t57 - t45 - t2 + t25) * s56) + 6.0_dp * t
     #5 + 6.0_dp * t10 - 4.0_dp * t34 + 4.0_dp * t38)) + t14) * t41) * s
     #truc1MP
      t2 = 0.1e1_dp / t7
      t3 = 0.1e1_dp / t6
      result = im * propW34 * t1 / ecossin ** 2.0_dp / gw ** 2.0_dp / s2
     #5 * t3 * t2 / 9.0_dp

           ampNonresonantHeavyImC7MP = result
       end function ampNonresonantHeavyImC7MP

       function ampNonresonantHeavyImC7PM()
           implicit none
           complex(dp) :: ampNonresonantHeavyImC7PM
           complex(dp) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           complex(dp) ::  t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30
           complex(dp) ::  t31,t32,t33,t34,t35,t36,t37,t38,t39,t4,t40,t41
           complex(dp) ::  t42,t43,t44,t45,t46,t47,t48,t49,t5,t50,t51,t52
           complex(dp) ::  t53,t54,t55,t56,t57,t58,t59,t6,t60,t61,t62,t63
           complex(dp) ::  t64,t65,t66,t67,t68,t69,t7,t70,t71,t72,t73,t74
           complex(dp) ::  t75,t76,t77,t78,t79,t8,t80,t81,t82,t83,t84,t85
           complex(dp) ::  t86,t87,t88,t89,t9,t90,t91,t92

           complex(dp) :: result

      t1 = intHLs250000x0122D4eps0()
      t2 = s25 + s26 + s56
      t3 = s16 + s26 - s34 + s56
      t4 = s13 - s45 - s46
      t5 = gb ** 2.0_dp
      t6 = t5 ** 2.0_dp
      t7 = gw ** 2.0_dp
      t8 = t7 ** 2.0_dp
      t9 = t3 * s25
      t10 = s12 + s16 + s26 - s34
      t11 = t3 * t2
      t12 = t11 * propW16
      t13 = t12 * t7 * t5
      t14 = t9 * propW16
      t15 = propZ25 * s25
      t16 = t15 * t8
      t17 = (s12 + s13 + s16 + s26 - s34 - s45 - s46) * t2
      t18 = 4.0_dp
      t19 = 3.0_dp
      t20 = s25 ** 2.0_dp
      t21 = s12 + s16 + s25
      t22 = -s16 + s25 + s34
      t23 = s12 + s16 + s26 - s34
      t24 = -s45 - s46
      t25 = s13 * s16
      t26 = (s12 + s13 + s16) * s25
      t27 = (-s13 + s25 + s26) * s34
      t28 = (t21 + s26) * s26 + s56 * t23 + t22 * t24 - t25 + t26 - t27
      t29 = 2.0_dp
      t30 = t29 * s26
      t31 = s12 + s26
      t32 = t29 * t31
      t33 = s26 + s16 - s34
      t34 = s56 ** 2.0_dp
      t35 = t33 * s13
      t36 = t24 * (t30 + s25 + s56 + s16 - s34)
      t37 = t12 * (t32 + s25)
      t38 = (s13 + s12 - s16 - s25 - s26 + s34 - s56) * s25 + s26 * (s12
     # - s26 - s56 + s13 - s16 + s34) + (2.0_dp * s13 - t30 - s25 + s12 
     #- s56 - s16 + s34 - s45 - s46) * s56 + t20 - t34 + t35 + t36 - t37
      t39 = t15 - t29
      t40 = s56 + s16 - s34
      t41 = t19 * s26
      t42 = s12 + s13
      t43 = t29 * t42
      t44 = -t43 + t3
      t45 = t29 * s12
      t46 = t29 * s25
      t47 = t24 * (t41 + t46 + t40)
      t48 = (-t29 * (s25 - s12 - s13 + s45 + s46) - t40 - t41) * s56
      t40 = -t18 * t34 + t29 * (s13 * t40 + (2.0_dp * s13 + t45 - s25 - 
     #s26 - s56 - s16 + s34) * s25 + (s13 - t44) * s26 + t20 + t47 + t48
     #) - t37 * t19
      t49 = t15 * t5
      t50 = t49 * t28
      t51 = t7 * t39
      t52 = t16 * t19 * t38
      t53 = intHLs250000x0113D4eps0()
      t54 = s16 - s34
      t42 = t18 * t42
      t55 = s16 + s25
      t56 = s25 + s26
      t57 = s13 - s45 - s46
      t58 = t18 * s12
      t59 = -t43 + s16
      t60 = t29 * s13
      t61 = s25 + s16 - s34
      t62 = t19 * s25
      t63 = t59 * s25
      t64 = t12 * (t18 * t31 + t46 + t54)
      t30 = t24 * (t30 + t61)
      t65 = t56 * s34
      t24 = t24 * (t41 + t46 + t54)
      t66 = 10.0_dp
      t67 = 18.0_dp * t34
      t68 = (8.0_dp * s12 + 12.0_dp * s13) * s26
      t69 = 6.0_dp * t12 * (t32 + t22)
      t70 = t29 * ((-5.0_dp * s16 + t42) * s25 + (-9.0_dp * s25 - 14.0_d
     #p * s26 - 5.0_dp * t54 + 6.0_dp * t57 + t58) * s56)
      t71 = t18 * (-s13 * t54 - t24)
      t66 = t66 * ((t55 + s26) * s26 - t65)
      t72 = t15 * t19
      t28 = t6 * t28 * t29
      t11 = 9.0_dp * t11 * t8 * (-propW16 * t54 + 1.0_dp)
      t73 = t7 * t5
      t74 = intHLs250000x0112D2eps1()
      t75 = t19 * epinv
      t76 = t75 - t18
      t77 = t19 * s13
      t78 = t77 + t45 - t55
      t79 = s25 - s13 + s45 + s46
      t80 = -s26 * t18 - t79 * t19 + t45 - t54
      t81 = 24.0_dp * t34
      t59 = 8.0_dp * t59 * s25 + 8.0_dp * s26 * (-t78 + s26) - 8.0_dp * 
     #s56 * t80 - 8.0_dp * t24 - 8.0_dp * t25 - 8.0_dp * t27
      t82 = 12.0_dp * t15 * t38 + 12.0_dp * t37
      t75 = t75 * t39
      t38 = epinv * t38
      t83 = intHLs250000x0111D2eps0()
      t84 = t18 * t54
      t31 = 6.0_dp * t31
      t85 = t19 * t57
      t86 = t19 * s16
      t87 = s12 + s13 - s45 - s46
      t88 = t29 * t87
      t89 = s26 ** 2.0_dp
      t90 = (s13 + s25 + s26) * s34
      t31 = -t18 * (s25 * (t43 + s16) + s26 * (t77 + t45 + t55 + s26) + 
     #s56 * (t85 + t45 - t22 - s56) + t24 + t25 - t90) + t15 * (-t18 * (
     #t65 - t89) + t29 * (s25 * (2.0_dp * s13 + t86 + t45 + s25 + s26 + 
     #s56 - s34) + s26 * (t29 * (s25 + s12 + s13) + t86 + s26 + s56 - s3
     #4) + s56 * (t19 * t33 + s13 + s56 + t88) - t20 + t35 + t47) - t12 
     #* (8.0_dp * t54 + t62 + t31)) + t12 * t29 * (t31 + t84 + t62)
      t33 = intHLs250000x0122D4eps1()
      t35 = -t18 + epinv
      t57 = t29 * t57
      t86 = t29 * t23
      t91 = t86 + s25
      t92 = propW16 * t91 * t2
      t21 = -9.0_dp * t16 * epinv * (s26 * (t60 + t21 + s26) + s56 * (t5
     #7 + t23) + t25 + t26 + t30 - t90 - t12 * t91) + t5 * (t7 * (epinv 
     #* (-t18 * (s26 * (t21 * t29 + t77) - s34 * (t29 * t56 + s13) + s56
     # * (t85 + t86) + t24 + t25) + t3 * (6.0_dp * t92 + t72 * (t92 - t5
     #7)) - 8.0_dp * t26 - 8.0_dp * t89) + t59 + t81 - 12.0_dp * t15 * (
     #t29 * t34 - (s12 + s13 - s16) * s25 - s26 * (t60 - s25 + s12 - s16
     # - s26) + s56 * (t29 * t79 - s12 + s16 - s34 + t41) - t25 - t27 - 
     #t30 + t37) + 12.0_dp * t37) + t50 * t35)
      t23 = intHLs250000x0113D4eps1()
      t26 = -t29 + epinv
      t37 = s25 + s26 + s56 + s16 - s34
      t41 = t37 * s13
      t57 = intHLs250000x0111D2eps1()
      t12 = t12 * (t32 + t61)
      t32 = s25 * intHLs250000x0112D2eps0()
      t11 = s25 * (t53 * (t73 * (t72 * (t19 * t34 + t29 * (-t25 - t30) +
     # (-s13 * t18 + s26 - t45 + t55) * s26 + (t60 - t56) * s34 + (t18 *
     # (s26 - s13 + s45 + s46) - t45 + t54 + t62) * s56 + t63 + t64) - t
     #66 - t71 + t70 - t69 + t68 - t67) + t15 * (-t11 + t28)) + t23 * (-
     #t73 * (epinv * (t72 * (-t29 * (-t36 - t41 + t34) - s25 * (-t45 + t
     #37 - s25) - s26 * t44 + t48 - t64) + t66 + t71 - t70 + t69 - t68 +
     # t67) - t81 - t82 - 8.0_dp * s26 * (-t78 + s26) + 8.0_dp * s56 * t
     #80 + 8.0_dp * t24 + 8.0_dp * t25 + 8.0_dp * t27 - 8.0_dp * t63) + 
     #t15 * (-t11 * epinv + t28 * t26)))
      t4 = (-t1 * (6.0_dp * t9 * t7 * t5 * (propZ25 * t4 + (-1.0_dp + (-
     #s12 - s16 - s26 + s34) * propZ25) * t2 * propW16) - 12.0_dp * t13 
     #* t10 - 9.0_dp * t16 * ((-t2 - t3) * s13 + s45 * t3 + (t2 + t3) * 
     #s46 + t2 * (t14 - s12 - s16 - s26 + s34 + s45)) + 8.0_dp * t17 * t
     #7 * t5 + t18 * t3 * t7 * t5 * t4 - 18.0_dp * t14 * t2 * propZ25 * 
     #t8 * t10 - t15 * t6 * (t3 * (-s13 + s45 + s46) + t17) - t13 * t19 
     #* propZ25 * t20) - t74 * (9.0_dp * t38 * t16 + t5 * (t7 * (-t75 * 
     #t40 + t59 + t81 + t82) + t50 * t76)) + t33 * t21) * s25
      t9 = s26 + s56
      t10 = t29 * t9 + t61
      t13 = t19 * t9 + t46 + t54
      t14 = t6 * t22
      t17 = t8 * t10
      t20 = t51 * t29 * t5
      t21 = -t20 * t13 + t15 * (t17 * t19 + t14)
      t24 = t29 * t54
      t25 = 6.0_dp * t9
      t27 = s25 * (propZ25 * t3 * t19 + t18) + t24 + t25
      t17 = 9.0_dp * t17
      t28 = t73 * t29
      t10 = propZ25 * t10
      t30 = (6.0_dp * t10 - 8.0_dp) * s25
      t9 = 12.0_dp * t9
      t37 = t17 * epinv
      t44 = -t57 * t26 + t32
      t23 = t23 * t26 + t53
      t48 = s25 * t5
      t53 = t19 * propW16
      t54 = t53 * t3
      t55 = -t54 + t29
      t59 = t29 * propW16
      t61 = -t59 * t3 + 1.0_dp
      t62 = t8 * t61
      t63 = t15 * (t62 * t19 + t6)
      t64 = -t28 * t39
      t66 = -t55
      t67 = 6.0_dp * propW16 * (1.0_dp - t15)
      t61 = -t61
      t68 = t15 + t29
      t54 = t54 * t68 - t18
      t69 = 6.0_dp * t15
      t70 = t69 * t61
      t71 = 12.0_dp * propW16
      t73 = t71 * t3
      t77 = t6 * t76
      t62 = 9.0_dp * t62 * epinv
      t78 = t6 * t35
      t59 = t59 * t2 + 1.0_dp
      t53 = t53 * t2 + 1.0_dp
      t79 = t8 * t59
      t80 = t15 * (-t79 * t19 + t6)
      t68 = t19 * (propW16 * t68 * t2 + t15) + t29
      t79 = 9.0_dp * t79
      t59 = t69 * t59
      t69 = t71 * t2
      t71 = t79 * epinv
      t8 = (t19 * (-t44 * (t20 * t66 + t63) + t83 * (t64 * t55 + t63)) +
     # 8.0_dp * t48 * (t7 * (-t67 * t3 + t18 - t72) + t49) * t23 + (t1 *
     # (t28 * t54 + t15 * (9.0_dp * t8 * t61 + t6)) + t33 * (t15 * (-t62
     # + t78) - t28 * (-epinv * t54 + t70 - t73 + 8.0_dp)) - t74 * (t15 
     #* (t62 + t77) - t28 * (-t75 * t66 + t70 - t73 + 8.0_dp))) * s25) *
     # t2 * (-struc3PM - struc5PM - struc8PM)
      t1 = struc2PM * (t19 * (t21 * t44 - t83 * t21) - 8.0_dp * t48 * (t
     #7 * (s25 * (-t10 * t19 + t18) + t24 + t25) + t49 * t22) * t23 - (t
     #1 * (-t28 * t27 + t15 * (-t17 + t14)) + t33 * (t15 * (t14 * t35 - 
     #t37) - t28 * (epinv * t27 - t30 + t84 + t9)) - t74 * (t15 * (t14 *
     # t76 + t37) + t28 * (-t75 * t13 + t30 - t84 - t9))) * s25) - t3 * 
     #(t19 * (-t44 * (t20 * t53 + t80) + t83 * (-t64 * t53 + t80)) + 8.0
     #_dp * t48 * (t7 * (-t67 * t2 - t29 + t72) + t49) * t23 + (t1 * (t2
     #8 * t68 + t15 * (t79 + t6)) + t33 * (t28 * (epinv * t68 + t18 - t5
     #9 + t69) + t15 * (t71 + t78)) - t74 * (t15 * (-t71 + t77) + t28 * 
     #(t75 * t53 + t18 - t59 + t69))) * s25) * struc7PM + t8
      t3 = 0.1e1_dp / t3
      t6 = 0.1e1_dp / t2
      result = im * propW34 * (-t29 * t1 - struc1PM * (-t18 * t11 - t19 
     #* (t83 * (t52 + t5 * (-t31 * t7 + t50)) + t57 * (-t19 * t16 * (-t3
     #8 - t29 * (t12 - t36 - t41) - s25 * (-t45 + s26 + s56 - s34) - s26
     # * (-t43 - s25 + s56 - s34) - s56 * (-t88 + t2 + s56) - t65) + t5 
     #* (-t7 * (epinv * t31 - t29 * t39 * (-t29 * (-s13 * s26 + s13 * (-
     #s16 - s25 - s56) + s34 * (s13 + t56) + t34 - t47 - t89) + s25 * (t
     #58 + t60 - s26 - s56 + s16 + s34) - s26 * (-t42 - t46 + s26 + s56 
     #- s16 - s34) + (t18 * t87 + s16 - s26 - s34 - s56 - t46) * s56 - t
     #12 * t19)) + t50 * t26)) - t32 * (t52 - t5 * (t40 * t51 - t50))) -
     # t4)) / ecossin ** 2.0_dp / gw ** 2.0_dp / s25 * t6 * t3 / 18.0_dp

           ampNonresonantHeavyImC7PM = result
       end function ampNonresonantHeavyImC7PM

       function ampNonresonantHeavyImC7PP()
           implicit none
           complex(dp) :: ampNonresonantHeavyImC7PP

           complex(dp) :: result

      result = 0.0_dp

           ampNonresonantHeavyImC7PP = result
       end function ampNonresonantHeavyImC7PP

       function ampNonresonantHeavyMM()
           implicit none
           complex(dp) :: ampNonresonantHeavyMM
           complex(dp) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           complex(dp) ::  t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30
           complex(dp) ::  t31,t32,t33,t34,t35,t36,t37,t38,t39,t4,t40,t41
           complex(dp) ::  t42,t43,t44,t45,t46,t47,t48,t49,t5,t50,t51,t52
           complex(dp) ::  t53,t54,t55,t56,t57,t58,t59,t6,t60,t61,t62,t63
           complex(dp) ::  t64,t65,t66,t67,t68,t69,t7,t70,t71,t72,t73,t74
           complex(dp) ::  t75,t76,t77,t78,t79,t8,t80,t81,t82,t83,t84,t9

           complex(dp) :: result

      t1 = intHLs250000x0122D4eps0()
      t2 = s25 + s26 + s56
      t3 = s16 + s26 - s34 + s56
      t4 = s26 + s56
      t5 = s56 ** 2.0_dp
      t6 = s13 * t3
      t7 = 2.0_dp
      t8 = (s25 + s26 + s56 - s13 + s45 + s46) * (s34 + s45 + s46 + s56)
      t9 = (s12 * t2 - t7 * s45 * t4 - s16 * (s45 + s46 + s56) + (-s25 +
     # s34) * s45 + s46 * (-s26 + s34) + (-s26 + s34 - s46) * s56 - t5 +
     # t6) * s12 + t8 * t3
      t10 = t7 * s16
      t11 = s25 + s16
      t12 = s12 - s13
      t13 = -s13 + s25 + s26
      t14 = s16 + s26
      t15 = s25 + s26
      t16 = t13 * t14
      t17 = -s16 + s25 + s34
      t18 = s25 + s26 - s13 + s34
      t19 = s16 + s26 - s34
      t20 = t18 * t19
      t21 = s12 * t17 + t20
      t22 = -s12 - s13 + s25 + s26 + s34
      t23 = t7 * s45
      t24 = t22 + t23
      t25 = 3.0_dp
      t26 = t25 * s12
      t27 = t7 * s12
      t28 = s25 - s13 + s16
      t29 = 5.0_dp
      t30 = t25 * s26
      t31 = (-t7 * (s45 + s16) + s12 + s13 - s25 + s34 - t30 - s46) * s4
     #6
      t32 = t7 * (s26 + s45 + s46)
      t8 = (-s25 - s56 + s13 - s46 + s16 - s34 - t23 + s12) * s12 + t8
      t23 = s45 ** 2.0_dp
      t33 = s46 ** 2.0_dp
      t34 = s13 * s16
      t35 = t13 * s34
      t36 = t19 * t23
      t37 = t19 * t33
      t38 = t19 * t24 * s46
      t39 = propW16 * t2
      t40 = t39 * t3 * t8
      t41 = s12 * s26 ** 2.0_dp
      t42 = t40 - t41
      t43 = t7 * t42
      t12 = -(-s25 * (t10 + s12) - s26 * (t11 * t7 + t12) + t34) * s12 -
     # (-s12 * (-t15 * t7 + s13) + t16 - t35) * s34 - s45 * t21 - ((t28 
     #- t26 + t32 + s56) * s56 - s12 * t12 - (t26 + s13) * s16 - (t27 - 
     #s16) * s25 - (s12 * t29 - s26 - t28) * s26 - (-t14 - t26 + s34) * 
     #s34 - (-s25 + s13 + s34 - t30 - t10 - s45) * s45 - t31) * s56 - t3
     #6 - t37 - t38 - t43
      t29 = t25 * s16
      t16 = s12 * (-t15 * t25 + s13) - t16
      t44 = t7 * s25
      t45 = -s26 + s16 - s34 - t44
      t46 = 4.0_dp
      t47 = t46 * s12
      t48 = 7.0_dp
      t49 = s12 ** 2.0_dp
      t26 = t7 * t49 - s12 * s13 + (t47 + s13) * s16 + (t26 - s16) * s25
     # + (s12 * t48 - s26 - t28) * s26 + (-t14 - t47 + s34) * s34 + (-s2
     #5 - s12 + s13 + s34 - t30 - t10 - s45) * s45 + t31
      t28 = -t28 + t47 - t32
      t31 = t7 * (s12 - s45 - s46) - t18
      t32 = (t22 + s45) * s45
      t47 = (t24 + s46) * s46
      t8 = t39 * t8
      t48 = -s12 * t13 + (-t31 + s56) * s56 + t35 + t32 + t47 + t8
      t11 = (-s25 * (t29 + t27) - s26 * (t11 * t25 - s13 + t27) + t34) *
     # s12
      t34 = 12.0_dp * t42
      t42 = propZ25 * s25
      t50 = 6.0_dp * t42
      t51 = t50 * t3 * t48
      t52 = gw ** 2.0_dp
      t53 = t52 ** 2.0_dp
      t54 = gb ** 2.0_dp
      t55 = t54 ** 2.0_dp
      t56 = t42 * t54 * t9
      t57 = t52 * (-t46 * ((-t16 - t35) * s34 - (s12 * t45 - t20) * s45 
     #+ ((-t28 + s56) * s56 - t26) * s56 + t36 + t37 + t38 + t11) - t34 
     #- t51)
      t58 = 9.0_dp
      t59 = t42 * t58 * t53
      t60 = intHLs250000x0111D2eps0()
      t61 = 0.1e1_dp / s25
      t62 = t54 * propZ25 * t9
      t63 = t58 * t53
      t64 = t63 * propZ25
      t65 = 1.0_dp - epinv
      t66 = s12 - s16 - s26 + s34 - s56
      t67 = -t7 * (s45 + s46) - t18
      t68 = s25 - s45
      t69 = s16 - s34
      t70 = t25 * t69
      t71 = t46 * s26
      t72 = t67 * t3
      t73 = (-s12 - s26 - s56 - s13 + s34 - t10) * s25
      t74 = (-s25 - s12 - 2.0_dp * s26 - s56 + s34 - t10) * s26
      t40 = t40 * t7
      t69 = t7 * t69
      t19 = -t19 * t67
      t67 = -t3 + t27
      t75 = s25 + s12
      t76 = s13 * (s12 + s16 + s25 + s26 - s34 + s56)
      t49 = t49 * s13
      t29 = ((s26 + s56 + s13 - s34 + t29 + t27) * s25 + s26 * (t7 * t75
     # + s26 - s34 + s56 + t29) - t76) * s12
      t77 = (s13 * s34 - t23 - t47) * t3
      t27 = t15 * s34 * (t3 + t27)
      t44 = s45 * (s12 * (-s26 + s56 + s16 - s34 - t44) - t18 * t3)
      t50 = -t3 * (t50 * t48 + 12.0_dp * t8) + t46 * (((t31 + t67) * s56
     # + t7 * s12 * (s12 + t68 + t30 + t69) - t19) * s56 + t49 + t29 + t
     #77 - t27 + t44) + 8.0_dp * t41
      t78 = intHLs250000x0112D2eps1()
      t14 = s13 * t14 - t15 * (s12 + s16 + s26)
      t79 = t68 + t71 + t70
      t10 = t13 * s34 ** 2.0_dp - ((s34 - s56 - t75 - t10) * s26 + t73 +
     # t76) * s12 - (s56 * t13 - t14) * s34 - (s56 * t22 + t21) * s45 - 
     #((-t31 - t66) * s56 + (-t79 - s12) * s12 + t19) * s56 - t3 * (t23 
     #+ t47) + t49 - t43
      t13 = t42 * t53
      t22 = t7 * t52
      t23 = t22 * t50 / 2.0_dp
      t43 = t65 * t54 * (t23 - t56)
      t47 = intHLs250000x0122D4eps1()
      t75 = epinv * t12
      t80 = intHLs250000x0111D2eps1()
      t81 = -t7 + epinv
      t82 = 18.0_dp
      t9 = t55 * t9
      t83 = intHLs250000x0113D4eps0()
      t22 = t22 * t54
      t84 = intHLs250000x0112D2eps0()
      t11 = t52 * (t46 * ((-t16 - t35) * s34 + (-s12 * t45 + t20) * s45 
     #+ ((-t28 + s56) * s56 - t26) * s56 + t11 + t36 + t37 + t38) + t34 
     #+ t51)
      t16 = intHLs250000x0121D2eps0()
      t20 = intHLs250000x0121D2eps1() * epinv
      t26 = intHLs250000x0113D4eps1() * t65
      t6 = -t83 * (t22 * (-t3 * (t42 * t48 * t25 + 6.0_dp * t8) + t7 * (
     #(s56 * t67 + t7 * s12 * (s12 + s25 + s56 - s45 + t30 + t69) + t72)
     # * s56 - t27 + t29 + t44 + t49 + t77) + t41 * t46) + t42 * (t63 * 
     #t10 - t9)) + t26 * (-t54 * (-t50 * t52 + t56) + t59 * (-(s13 * (s1
     #6 + s25 + s26 - s34 + s56) + s45 * (s25 - s56 - s16 + s34) + t73 +
     # t74) * s12 - s34 * (t15 * (s12 + s16 + s26 - s34 + s56) - t6) - (
     #-s56 * t66 - (s56 * t7 + s12 + t68 + t70 + t71) * s12 - t72) * s56
     # - t3 * ((t18 + s45) * s45 + (t24 + s46) * s46) - t40))
      t8 = t47 * t65 + t61 * (-t80 * t81 + epinv2 - t60) - t65 * t78 - t
     #1 + t16 + t20 + t84
      t15 = t7 * (-t26 + t83)
      t18 = 0.1e1_dp / s12
      t26 = 0.1e1_dp / t3
      t1 = t7 * mytree * (-t7 * t6 + epinv2 * (-t64 * t10 + t54 * (t11 *
     # t61 + t62)) + t1 * (t59 * t12 - t54 * (t56 - t57)) + t60 * (t54 *
     # (t57 * t61 - t62) + t64 * t12) + t78 * (t58 * t13 * (-epinv * t10
     # + t12) + t43) + t47 * (-t58 * t13 * (t10 - t75) - t43) + t80 * (t
     #23 * t81 * t54 * t61 - t9 * t81 * propZ25 + t53 * propZ25 * (t82 *
     # (-(-t74 - t73 - t76) * s12 - (t14 + t35) * s34 + s45 * t21 + s56 
     #* t5 - (-(-t79 - s12) * s12 - t19 - t33 - t32 - t35) * s56 - t5 * 
     #(s12 - s26 - s16 + s34 + t31) + t36 + t37 - t49 + t24 * s46 * t3 +
     # t40) + t75 * t58)) + (t84 + t16 + t20) * (-t59 * t10 + t54 * (t11
     # + t56))) * t18 * t26
      t5 = s12 + s16 + s26 - s34 - s45
      t2 = 0.1e1_dp / t2
      t5 = 0.1e1_dp / t5
      result = propW34 * (t1 + (s12 - s34 - s45 - s46) * ((-t8 + t15) * 
     #(t42 * (t17 * t55 - t63 * (t4 * t7 + s16 + s25 - s34)) - t22 * ((p
     #ropZ25 * t3 * t25 + t46) * s25 + 6.0_dp * t4 + t69)) * struc4Step5
     # * t26 + (t8 - t15) * (t42 * (t63 * (t39 * t7 + 1.0_dp) + t55) + t
     #22 * (t25 * (t39 * (t42 + t7) + t42) + t7)) * ((s12 - s34 - s45 - 
     #s46 - s56) * struc9Step5 * t18 + struc8Step5))) / ecossin ** 2.0_d
     #p / gw ** 2.0_dp * t5 * t2 / 18.0_dp

           ampNonresonantHeavyMM = result
       end function ampNonresonantHeavyMM

       function ampNonresonantHeavyMP()
           implicit none
           complex(dp) :: ampNonresonantHeavyMP

           complex(dp) :: result

      result = 0.0_dp

           ampNonresonantHeavyMP = result
       end function ampNonresonantHeavyMP

       function ampNonresonantHeavyPM()
           implicit none
           complex(dp) :: ampNonresonantHeavyPM

           complex(dp) :: result

      result = 0.0_dp

           ampNonresonantHeavyPM = result
       end function ampNonresonantHeavyPM

       function ampNonresonantHeavyPP()
           implicit none
           complex(dp) :: ampNonresonantHeavyPP
           complex(dp) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           complex(dp) ::  t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30
           complex(dp) ::  t31,t32,t33,t34,t35,t36,t37,t38,t39,t4,t40,t41
           complex(dp) ::  t42,t43,t44,t45,t46,t47,t48,t49,t5,t50,t51,t52
           complex(dp) ::  t53,t54,t55,t56,t57,t58,t59,t6,t60,t61,t62,t63
           complex(dp) ::  t64,t65,t66,t67,t68,t69,t7,t70,t71,t72,t73,t74
           complex(dp) ::  t75,t76,t77,t78,t79,t8,t80,t81,t82,t83,t84,t85
           complex(dp) ::  t86,t87,t88,t89,t9,t90

           complex(dp) :: result

      t1 = intHLs250000x0112D2eps0()
      t2 = 3.0_dp
      t3 = t2 * s25
      t4 = t3 * propZ25
      t5 = t4 - 2.0_dp
      t6 = gw ** 2.0_dp
      t7 = gb ** 2.0_dp
      t8 = propZ25 * s25
      t9 = t8 * t7
      t10 = t5 * t6 + t9
      t11 = intHLs250000x0113D4eps0()
      t12 = intHLs250000x0121D2eps0()
      t13 = intHLs250000x0122D4eps0()
      t14 = intHLs250000x0112D2eps1()
      t15 = 1.0_dp - epinv
      t16 = intHLs250000x0113D4eps1()
      t17 = intHLs250000x0122D4eps1()
      t18 = intHLs250000x0111D2eps0()
      t19 = -s16 + s25 + s34
      t20 = s26 + s56
      t21 = s25 + s16 - s34 + 2.0_dp * t20
      t22 = s16 - s34
      t23 = s16 + s26 - s34 + s56
      t24 = intHLs250000x0111D2eps1()
      t25 = -2.0_dp + epinv
      t26 = 0.1e1_dp / t23
      t27 = 0.1e1_dp / s25
      t28 = (-t14 + t17) * t15
      t29 = epinv2 * t27
      t30 = intHLs250000x0121D2eps1() * epinv
      t31 = t24 * t25
      t32 = t16 * t15
      t33 = -t32 + t11
      t34 = t10 * t33
      t35 = s25 + s26 + s56
      t36 = 2.0_dp * s26
      t37 = t36 + s25 + s16 - s34 + s56 + s12 - s13 + s46
      t38 = -1.0_dp + t8
      t39 = t23 * t37
      t40 = t39 * t6
      t41 = t40 * propW16 * t35 * t38
      t42 = (s12 + s16 + s25 + s26 - s34 + s56) * t35
      t43 = t9 * (t39 + t42)
      t44 = t8 * t6
      t45 = t44 * (-t39 + t42)
      t46 = t42 * t6
      t47 = s25 ** 2.0_dp
      t48 = t35 * t47
      t49 = t40 * epinv
      t50 = 12.0_dp
      t51 = 8.0_dp
      t52 = -t38
      t53 = 6.0_dp * propW16
      t54 = t53 * t52
      t55 = t6 * (-t35 * t54 + t4 - 2.0_dp) + t9
      t56 = t7 * propW34
      t57 = t56 * t55
      t37 = 0.1e1_dp / t37
      t58 = 0.1e1_dp / t35
      t59 = s26 * t26 * t37
      t60 = 0.1e1_dp / gw ** 2.0_dp
      t61 = s26 * t35
      t62 = t61 - t39
      t63 = t61 + t39
      t64 = 2.0_dp * t23
      t65 = t9 * t21
      t66 = -t6 * ((propZ25 * (-propW16 * t35 * t64 + t19) * t2 - 4.0_dp
     #) * s25 + 2.0_dp * s16 - 2.0_dp * s26 - 2.0_dp * s34 - 2.0_dp * s5
     #6 + t53 * t35 * t23) + t65
      t67 = t2 * propZ25
      t68 = s34 * t67 - 4.0_dp
      t69 = t67 * t47
      t70 = 2.0_dp * s26 + 2.0_dp * s34 + 2.0_dp * s56
      t71 = (t6 * (s16 * t5 - s25 * t68 - t69 + t70) + t65) * t60 * t58 
     #* t26 - t54
      t72 = epinv2 * propW16
      t73 = -t30 - t1 - t12
      t74 = t23 * t54 + t4 - 4.0_dp
      t75 = t6 * t74 - t9
      t74 = -t6 * t74 + t9
      t76 = s25 - s34 + s56 - s13
      t77 = s16 + s12 + s46
      t78 = t2 * s26
      t53 = s26 * t50 - t53 * t39 + t51 * t76 + 4.0_dp * t77 - t4 * (-2.
     #0_dp * propW16 * t39 - 2.0_dp * s13 + 2.0_dp * s25 - 2.0_dp * s34 
     #+ 2.0_dp * s56 + t77 + t78)
      t76 = t9 * (2.0_dp * t76 + t77 + t78)
      t77 = t31 * (t53 * t6 + t76)
      t53 = t18 * t56 * (-t53 * t6 - t76)
      t52 = t56 * propW16 * t52
      t76 = s12 - s13
      t79 = 2.0_dp * s13
      t80 = 5.0_dp * s25
      t81 = 2.0_dp * s16
      t82 = 2.0_dp * s25
      t83 = s13 - s46
      t84 = (t76 + s16) * s16
      t85 = s56 ** 2.0_dp
      t86 = s26 ** 2.0_dp
      t87 = 2.0_dp * s12
      t88 = 5.0_dp * s16
      t89 = t2 * s12
      t90 = 5.0_dp * s26
      t3 = s25 * (s16 * t2 + s12 - t79) + s26 * (-s13 * t2 + 7.0_dp * s2
     #5 + t87 + t88) + s34 * (-t3 - t90 - t81 - t76 + s34) + s46 * (t82 
     #+ t78 + t22) + s56 * (-t2 * t83 + 6.0_dp * s25 + 10.0_dp * s26 + 4
     #.0_dp * t22 + t87) + t84 + 4.0_dp * t85 + 6.0_dp * t86 + 2.0_dp * 
     #t47
      t76 = s25 * (-4.0_dp * s13 + t87 + t88) + s26 * (-5.0_dp * s13 + 7
     #.0_dp * s16 + 13.0_dp * s25 + t89) + s34 * (-7.0_dp * s26 - t76 - 
     #t80 - t81 + s34) + s46 * (4.0_dp * s25 + t22 + t90) + s56 * (11.0_
     #dp * s25 + 17.0_dp * s26 + 6.0_dp * t22 - 5.0_dp * t83 + t89) + t8
     #4
      t78 = t9 * (s25 * (-t79 + s16 + s12) - s26 * (-t80 + s16 + s13) - 
     #s34 * (-t81 + s25 - s26 - s12 + s13 + s34) + s46 * (t82 - s16 + s3
     #4 + s26) + s56 * (4.0_dp * s25 + 4.0_dp * s26 - t83) - t84 + 2.0_d
     #p * t85 + 2.0_dp * t86 + 2.0_dp * t47)
      t79 = t30 + t1 + t12 - t13
      t80 = t23 + t35
      t81 = -t23 + t35
      t82 = 4.0_dp * t35
      t64 = 6.0_dp * t6 * propW16 * t35 * t23 * t38 - t2 * t44 * t81 + t
     #6 * (t82 - t64) + t9 * t80
      t83 = epinv2 * t6
      t84 = t17 * t15
      t87 = t14 * t15
      t16 = t26 * t27 * (struc35PP * t74 * (s25 * (t87 - t30 - t84 - t1 
     #- t12 + t13) - epinv2 + t18 + t31 + 2.0_dp * s25 * (t16 * (-1.0_dp
     # + epinv) + t11)) - struc36PP * (s25 * (-t87 + t30 + t84 + t1 + t1
     #2 - t13 - 2.0_dp * t33) + epinv2) * (t6 * (-t4 + 4.0_dp) + t9)) + 
     #struc10PP * (-t10 * (t30 + t29 + t28 + t1 + t12 - t13) - t27 * (t3
     #1 + t18) * t26 * (t6 * ((-propZ25 * t2 * t21 + 4.0_dp) * s25 + 6.0
     #_dp * t20 + 2.0_dp * t22) + t9 * t19) + 2.0_dp * t34) * t58
      t8 = struc34PP * (-t2 * t83 * t9 * propW34 * t81 + t8 * epinv2 * t
     #80 * propW34 * t7 ** 2.0_dp - t18 * t23 * t57 + t56 * (s25 * (t28 
     #* t66 + t64 * t79) - t23 * t24 * t55 * t25) + 2.0_dp * t56 * (s25 
     #* (-t11 * t64 + t32 * t66) - t83 * t23) + 6.0_dp * t72 * t56 * t6 
     #* t35 * t23 * t38 + t82 * t83 * t56) * t26 * t27
      t3 = -struc20PP * (t57 * (t1 + t12 - t13) + t56 * (t55 * (t30 + t2
     #9 + t28) + t27 * (-t18 * (-t2 * t45 - 2.0_dp * t40 + 6.0_dp * t41 
     #+ t43 + 4.0_dp * t46) + t24 * (t2 * t45 * epinv - t43 * epinv + t4
     #1 * t50 + t46 * t51 + 6.0_dp * t6 * (((t35 * (-s12 - s16 - s26 + s
     #34 - s56) + t39 * (-epinv * propW16 * t35 + 1.0_dp)) * s25 - t48) 
     #* propZ25 + t39 * epinv * propW16 * t35) + 2.0_dp * ((t35 * (s12 +
     # s16 + s26 - s34 + s56) + t39) * s25 + t48) * propZ25 * t7 + 2.0_d
     #p * t49 - 4.0_dp * t6 * (epinv * t42 + t39))) * t37 * t26) - 2.0_d
     #p * t57 * t33) + struc33PP * (t56 * (t10 * t79 + t27 * (epinv2 * t
     #10 + t31 * (t6 * (s25 * (s25 * t51 - t67 * t3) + 2.0_dp * t76 + 14
     #.0_dp * t85 + 20.0_dp * t86) + t78) * t26 * t37) + t28 * t10) + t1
     #8 * t7 * (-t6 * propW34 * (t4 * t3 - 8.0_dp * t47 - 2.0_dp * t76 -
     # 14.0_dp * t85 - 20.0_dp * t86) + t78 * propW34) * t27 * t26 * t37
     # - 2.0_dp * t34 * t56) - t8
      t2 = t56 * (t60 * t16 + struc30PP * (-t2 * epinv2 * propZ25 * (t58
     # - t26) - 6.0_dp * t72 * (propZ25 - t27) + t13 * t71 + t18 * ((t2 
     #* (t59 + t58) + t7 * (-t59 + t58) * t60) * propZ25 + t27 * (-2.0_d
     #p * t26 * t36 * t37 - t54 - 2.0_dp * t58)) - t60 * (t28 * t66 * t5
     #8 * t26 + (t58 + t26) * propZ25 * t7 * epinv2) - t73 * (-(-t6 * (-
     #s16 * t5 + s25 * t68 + t69 - t70) + t65) * t60 * t58 * t26 + t54) 
     #+ t24 * (-t54 * t25 + (t2 * t44 * epinv * t63 - t9 * epinv * t62 -
     # 4.0_dp * t6 * (epinv * t61 - t39) + 2.0_dp * t9 * t62 - 2.0_dp * 
     #t49 - 6.0_dp * t44 * t63 + t61 * t51 * t6) * t60 * t58 * t26 * t37
     #) * t27 + 2.0_dp * t11 * t71 + 2.0_dp * t58 * (-t26 * t32 * t60 * 
     #t66 + t29) - 4.0_dp * t29 * t26))
      t4 = 0.1e1_dp / 0.9e1_dp
      result = t4 * (t60 * (t3 * t58 + (-t53 * t27 * t37 + t56 * (-t13 *
     # t75 + t27 * (-epinv2 * t74 + t77 * t37) - t73 * t75 - t28 * t74) 
     #- 2.0_dp * t56 * (t11 * t75 + t32 * t74)) * struc31PP * t26) + str
     #uc32PP * (-t50 * t52 * t33 + 6.0_dp * t52 * ((-t14 + t17) * t15 + 
     #t1 + t12 - t13 + t29 + t30) - t60 * t27 * t26 * t37 * (-t77 * t56 
     #+ t53)) + t2) / ecossin ** 2.0_dp

           ampNonresonantHeavyPP = result
       end function ampNonresonantHeavyPP

       function ampNonresonantHeavyReC4MM()
           implicit none
           complex(dp) :: ampNonresonantHeavyReC4MM

           complex(dp) :: result

      result = 0.0_dp

           ampNonresonantHeavyReC4MM = result
       end function ampNonresonantHeavyReC4MM

       function ampNonresonantHeavyReC4MP()
           implicit none
           complex(dp) :: ampNonresonantHeavyReC4MP
           complex(dp) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           complex(dp) ::  t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30
           complex(dp) ::  t31,t32,t33,t34,t35,t36,t37,t4,t5,t6,t7,t8
           complex(dp) ::  t9

           complex(dp) :: result

      t1 = intHLs250000x0121D2eps0()
      t2 = intHLs250000x0121D2eps1() * epinv
      t3 = (t2 + t1) * s25 + epinv2
      t4 = gb ** 2.0_dp
      t5 = propZ25 * s25
      t6 = t5 * gw ** 2.0_dp
      t7 = s25 + s26 + s56
      t8 = s16 + s26 - s34 + s56
      t9 = 2.0_dp * propW16
      t10 = 6.0_dp * propW16
      t11 = 3.0_dp * t6
      t12 = t4 * (t10 * t7 + t5 + 2.0_dp) + t11 * (t9 * t7 + 1.0_dp)
      t13 = intHLs250000x0112D2eps1()
      t14 = intHLs250000x0113D4eps1()
      t15 = intHLs250000x0122D4eps1()
      t16 = 2.0_dp * t14
      t17 = (-t16 + t2 + t1 + t13 - t15) * s25 + epinv2
      t9 = t9 * t8 - 1.0_dp
      t18 = s25 + s16 - s34
      t19 = s26 + s56
      t20 = t18 + 2.0_dp * t19
      t21 = s34 ** 2.0_dp
      t22 = s16 ** 2.0_dp
      t23 = s25 ** 2.0_dp + t21 + t22
      t24 = s56 ** 2.0_dp
      t21 = s25 * (s16 * s25 - s25 * s34 + t21 + t22)
      t22 = s12 * t7
      t25 = s26 * propW16
      t26 = t25 * ((s25 + s56) * s16 + (-s25 - s56) * s34 + s56 * (s25 +
     # s26))
      t27 = s26 ** 2.0_dp
      t28 = s16 * s34
      t28 = (((-s16 + s34) * s56 + t28) * s25 - s26 * t27 + (t28 - t24) 
     #* s26 + t28 * s56) * propW16 - s45 * t19 - s46 * t19
      t29 = t27 * propW16 * (s25 + s16 - s34)
      t10 = t10 * t8 + t5 - 4.0_dp
      t30 = -4.0_dp * s25 + 2.0_dp * s34 - 6.0_dp * t19 - s16 * (t5 + 2.
     #0_dp) + t5 * (s25 + s34)
      t31 = -1.0_dp + t5
      t32 = s25 * (-1.0_dp + (s25 + s45 + s46) * propZ25)
      t33 = t31 * s26
      t31 = t31 * s56
      t34 = -s45 - s46
      t35 = t4 * t10
      t36 = t11 * t9 + t35
      t37 = (t2 + t1 - t15) * t4
      t10 = t12 * t8 * t3 * (struc15MP + struc41MP + struc57MP) - t17 * 
     #(t4 * (12.0_dp * t26 + 9.0_dp * t29 - 6.0_dp * t28 + 3.0_dp * (pro
     #pW16 * t23 + s25) * s56 + 3.0_dp * t24 * (propW16 * t18 + 1.0_dp) 
     #+ 3.0_dp * t25 * t23 + 3.0_dp * t21 * propW16 + s13 * t30 + (t33 +
     # t31 + t32) * s16 + (s26 * (propZ25 * t19 - 1.0_dp) + t5 * (s26 - 
     #s45 - s46)) * s25 + (-t33 - t31 - t32) * s34 - t27 + t22 * t10 - 2
     #.0_dp * s16 * t34 + 2.0_dp * s26 * s56 + 2.0_dp * s34 * t34 - 4.0_
     #dp * s25 * t34) + t11 * (4.0_dp * t26 + (s26 * t23 + s56 * t23 + t
     #18 * t24 + t21) * propW16 - s13 * t20 + s45 * t18 + s46 * t18 + s5
     #6 * (s25 + s26) + t24 + t22 * t9 - 2.0_dp * t28 + 3.0_dp * t29)) *
     # struc1MP
      t18 = 0.1e1_dp / t8
      t19 = 0.1e1_dp / t7
      result = -0.2e1_dp / 0.3e1_dp * propW34 * (-2.0_dp * t17 * (-struc
     #2MP * (t11 * t20 - t4 * t30) + t7 * t36 * (struc3MP + struc8MP)) +
     # 2.0_dp * t7 * (epinv2 * (t4 * (-t5 + 4.0_dp) + t11) + s25 * (6.0_
     #dp * t6 * t14 * t9 + 6.0_dp * (t4 * t8 + t6 * t8) * t15 * propW16 
     #- t13 * t36 - t37 * t5 + 3.0_dp * t6 * (t2 + t1 - t15) + 4.0_dp * 
     #t37 + t35 * t16)) * struc5MP - 2.0_dp * (-t16 + t13 - t15) * s25 *
     # t8 * t12 * struc7MP + t10 + 6.0_dp * t8 * t7 * (t6 + t4) * t3 * p
     #ropW16 * (struc14MP + struc40MP + struc56MP)) / ecossin ** 2.0_dp 
     #/ s25 * t19 * t18

           ampNonresonantHeavyReC4MP = result
       end function ampNonresonantHeavyReC4MP

       function ampNonresonantHeavyReC4PM()
           implicit none
           complex(dp) :: ampNonresonantHeavyReC4PM
           complex(dp) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           complex(dp) ::  t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30
           complex(dp) ::  t31,t32,t33,t34,t35,t36,t37,t38,t39,t4,t40,t41
           complex(dp) ::  t42,t43,t44,t45,t46,t5,t6,t7,t8,t9

           complex(dp) :: result

      t1 = intHLs250000x0113D4eps1()
      t2 = s12 + s13
      t3 = 2.0_dp
      t4 = t3 * s13
      t5 = s25 + s16 - s34
      t6 = t3 * s26
      t7 = s13 - s45 - s46
      t8 = s25 + s26 + s56
      t9 = s16 + s26 - s34 + s56
      t10 = s16 - s34
      t11 = -s45 - s46
      t12 = propW16 * (t3 * (s12 + s26) + t5) * t8 * t9
      t13 = s25 + s26
      t14 = s25 + s26 - s34 + s56
      t15 = s12 + s13 + s16
      t16 = t6 + s25 + s16 - s34 + s56
      t17 = (s16 + s25 + s26 - s34 + s56) * s13
      t18 = t13 * s46
      t19 = (s12 * t3 + s25 - t14) * s25
      t20 = (t3 * (s13 + s12 - s45 - s46) - t8 - s56) * s56
      t21 = t13 * s34
      t22 = intHLs250000x0121D2eps0()
      t23 = intHLs250000x0112D2eps1()
      t24 = intHLs250000x0122D4eps1()
      t25 = s16 + 4.0_dp * t2
      t26 = 6.0_dp
      t27 = t26 * s13
      t28 = 4.0_dp * s12
      t29 = 3.0_dp * s26
      t30 = s25 * t3 + t10 + t29
      t31 = 3.0_dp * s25
      t32 = s25 - s16 + s34
      t33 = s13 * s16
      t34 = t15 * s25
      t35 = (-s12 - s16 - s25 - s26) * s26 + (s25 - s13 + s26) * s34 - (
     #s12 + s16 + s26 - s34) * s56 - t11 * t32 + t33 - t34
      t36 = s56 ** 2.0_dp
      t37 = (t28 + t27 + s25 + s16 + s26) * s26
      t13 = (t4 + t13) * s34
      t38 = (t26 * t7 + t10 + t28 - t31 - t6) * s56
      t39 = t3 * (t11 * t30 + t33)
      t40 = -3.0_dp * t12 - 3.0_dp * t36 + (propZ25 * t35 + t25) * s25 +
     # t37 - t13 + t38 + t39
      t41 = propZ25 * s25
      t42 = intHLs250000x0121D2eps1() * epinv
      t43 = (t42 + t22) * s25 + epinv2
      t44 = gb ** 2.0_dp
      t45 = gw ** 2.0_dp
      t46 = t41 * t45
      t6 = -t24 * t44 * (-t3 * propW34 * (s25 * t25 + t41 * t35 - t13 - 
     #3.0_dp * t36 + t37 + t38 + t39) + t12 * t26 * propW34) + (4.0_dp *
     # t44 * t40 + 12.0_dp * t46 * (t10 * s13 + s25 * t2 + (s12 + t4) * 
     #s26 + (t3 * t7 + s12 - s25 - s26 - s56) * s56 + t11 * (t6 + t5) - 
     #t12)) * t1 * propW34
      t7 = t41 + t3
      t13 = 3.0_dp * t46
      t35 = s26 + s56
      t37 = t3 * propW16
      t38 = t26 * propW16
      t1 = ((-t1 * t3 + t22 + t23 - t24 + t42) * s25 + epinv2) * propW34
     # * ((t44 * (-s34 * t3 + t26 * t35 + 4.0_dp * s25 + s16 * t7 - t41 
     #* (s25 + s34)) + t13 * (t3 * t35 + t5)) * struc2PM - t9 * (t44 * (
     #t38 * t8 + t3 + t41) + t13 * (t37 * t8 + 1.0_dp)) * struc7PM + (t4
     #4 * (t38 * t9 + t41 - 4.0_dp) + t13 * (t37 * t9 - 1.0_dp)) * t8 * 
     #(-struc3PM - struc5PM - struc8PM))
      t5 = 0.1e1_dp / t9
      t8 = 0.1e1_dp / t8
      t35 = 0.1e1_dp / 0.3e1_dp
      result = t35 * (-t3 * (t44 * t7 + t13) * ((t42 + t22) * s25 + epin
     #v2) * t9 * propW34 * (struc27PM + struc43PM - struc57PM) + (s25 * 
     #t6 - t3 * propW34 * t44 * (t43 * (t3 * (-s45 * t30 + t33) + 3.0_dp
     # * s26 ** 2.0_dp - 3.0_dp * t12 - 3.0_dp * t36 + (propZ25 * (-(s12
     # + s25) * s26 + (s25 - s13) * s34 + s45 * t32 - (s12 + s16 - s34 -
     # s46) * s56 + t18 + t33 - t34) + t25) * s25 + (3.0_dp * s16 + s25 
     #+ t27 + t28) * s26 - (s25 + t4 + t29) * s34 + (t26 * (s13 - s45) +
     # t10 - t31 + 4.0_dp * s12 - 4.0_dp * s46) * s56 - 4.0_dp * t18) + 
     #s25 * t23 * t40) - 3.0_dp * t41 * propW34 * t45 * (t43 * (t3 * (-s
     #45 * t16 - t12 + t17 - t18) + (t15 * t3 + s26 + t14) * s26 + t19 +
     # t20 - t21) + s25 * (t23 - t24) * (t3 * (t11 * t16 - t12 + t17) + 
     #(t2 * t3 + s25 + s34 - s56) * s26 + t19 + t20 - t21))) * struc1PM 
     #- 4.0_dp * t1) / ecossin ** 2.0_dp / s25 * t8 * t5

           ampNonresonantHeavyReC4PM = result
       end function ampNonresonantHeavyReC4PM

       function ampNonresonantHeavyReC4PP()
           implicit none
           complex(dp) :: ampNonresonantHeavyReC4PP

           complex(dp) :: result

      result = 0.0_dp

           ampNonresonantHeavyReC4PP = result
       end function ampNonresonantHeavyReC4PP

       function ampNonresonantHeavyReC7MM()
           implicit none
           complex(dp) :: ampNonresonantHeavyReC7MM

           complex(dp) :: result

      result = 0.0_dp

           ampNonresonantHeavyReC7MM = result
       end function ampNonresonantHeavyReC7MM

       function ampNonresonantHeavyReC7MP()
           implicit none
           complex(dp) :: ampNonresonantHeavyReC7MP
           complex(dp) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           complex(dp) ::  t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30
           complex(dp) ::  t31,t32,t33,t34,t35,t36,t37,t38,t39,t4,t40,t41
           complex(dp) ::  t42,t43,t44,t45,t46,t47,t48,t49,t5,t50,t51,t52
           complex(dp) ::  t53,t54,t55,t56,t57,t58,t59,t6,t60,t61,t62,t63
           complex(dp) ::  t64,t65,t66,t67,t68,t69,t7,t70,t71,t72,t73,t74
           complex(dp) ::  t75,t76,t77,t78,t79,t8,t80,t81,t82,t83,t84,t85
           complex(dp) ::  t86,t87,t88,t89,t9,t90,t91,t92,t93,t94,t95,t96
           complex(dp) ::  t97,t98,t99

           complex(dp) :: result

      t1 = intHLs250000x0112D2eps0()
      t2 = 2.0_dp
      t3 = t2 * s26
      t4 = s26 + s12
      t5 = t2 * t4
      t6 = s25 + s26 + s56
      t7 = s16 + s26 - s34 + s56
      t8 = -s45 - s46
      t9 = s56 ** 2.0_dp
      t10 = s25 ** 2.0_dp
      t11 = propW16 * (t5 + s25) * t6
      t12 = t11 * t7
      t13 = t8 * (t3 + s25 + s16 - s34 + s56)
      t14 = (s16 + s26 - s34) * s13 + (s13 + s12 - s16 - s25 - s26 + s34
     # - s56) * s25 + s26 * (s13 - s16 - s26 + s34 - s56 + s12) + (2.0_d
     #p * s13 - t3 - s25 - s16 + s34 - s56 + s12 - s45 - s46) * s56 - t9
     # + t10 - t12 + t13
      t15 = s12 + s13
      t16 = t2 * t15
      t17 = t16 - s16
      t18 = s16 + s25
      t19 = 3.0_dp
      t20 = t19 * s13
      t21 = t2 * s12
      t22 = -t21 - t20 + t18
      t23 = s16 - s34
      t24 = t19 * s26
      t25 = t2 * s25
      t26 = s25 - s13 + s45 + s46
      t27 = -t19 * t26 - 4.0_dp * s26 + t21 - t23
      t28 = s12 + s16 + s25
      t29 = -s16 + s25 + s34
      t30 = s12 + s16 + s26 - s34
      t31 = s13 * s16
      t32 = t8 * t29
      t33 = (-s13 + s25 + s26) * s34
      t34 = (s12 + s13 + s16) * s25
      t35 = (-t28 - s26) * s26 - s56 * t30 + t31 - t32 + t33 - t34
      t36 = t8 * (t25 + t24 + t23)
      t37 = propZ25 * t35
      t38 = 6.0_dp
      t39 = t12 * t19
      t40 = gb ** 2.0_dp
      t41 = t40 ** 2.0_dp
      t42 = gw ** 2.0_dp
      t43 = t42 ** 2.0_dp
      t44 = t19 * s25
      t45 = t44 * t42
      t46 = t45 * propZ25
      t47 = intHLs250000x0113D4eps0()
      t48 = s25 + s16 - s12
      t49 = t2 * s13
      t50 = s25 + s16 - s34
      t51 = s16 - s34 - s12
      t52 = t8 * (t3 + t50)
      t53 = (s12 + s13 - s16) * s25
      t26 = t2 * t9 - (t49 - t48 - s26) * s26 + (t2 * t26 + t24 + t51) *
     # s56 + t12 - t31 - t33 - t52 - t53
      t54 = 4.0_dp * t15
      t55 = s25 + s26
      t56 = s13 - s45 - s46
      t57 = 9.0_dp
      t58 = t38 * t56
      t59 = 4.0_dp * s12
      t60 = t2 * t56
      t61 = propW16 * (t2 * t30 + s25) * t6
      t62 = t23 * s13
      t63 = 10.0_dp
      t64 = t44 * propZ25
      t65 = propW16 * t6
      t63 = (-8.0_dp * s12 - 12.0_dp * s13) * s26 - t2 * ((-5.0_dp * s16
     # + t54) * s25 + (-s25 * t57 - 14.0_dp * s26 - 5.0_dp * t23 + t58 +
     # t59) * s56) + t7 * (t65 * (t5 + t29) * t38 + t64 * (-t60 + t61)) 
     #+ 18.0_dp * t9 + t63 * ((t18 + s26) * s26 - s34 * t55) - 4.0_dp * 
     #t36 - 4.0_dp * t62
      t66 = propZ25 * s25
      t67 = t66 * t40
      t68 = t67 * t35
      t69 = t66 * t57 * t43
      t70 = intHLs250000x0122D4eps0()
      t71 = t38 * s13
      t72 = 4.0_dp * t28 + t71
      t73 = 4.0_dp * t30 + t58
      t61 = t61 * t7
      t74 = (s13 + s25 + s26) * s34
      t28 = (t49 + t28 + s26) * s26 + (t60 + t30) * s56 + t31 + t34 + t5
     #2 - t61 - t74
      t30 = (t2 * t55 + s13) * s34
      t75 = s26 ** 2.0_dp
      t34 = 4.0_dp * t75 + 4.0_dp * t34
      t76 = t40 * propW34
      t77 = intHLs250000x0111D2eps0()
      t78 = 4.0_dp * t23
      t4 = -t2 * ((t16 + s16) * s25 + (t21 + t20 + t18 + s26) * s26 + (t
     #19 * t56 - s56 + t21 - t29) * s56 + t31 + t36 - t74) + t66 * (t2 *
     # (t65 * t23 * t7 - t9) - (t48 + s26) * s26 + (-t25 - t24 - t51) * 
     #s56 - t31 + t32 + t53 + t74) + t65 * (t38 * t4 + t44 + t78) * t7
      t20 = intHLs250000x0112D2eps1()
      t32 = t2 * (s25 - s13 + s26 + s56 + s45 + s46) + t11
      t48 = s16 - s34 + s56
      t51 = s25 * t17 + (-t22 - s26) * s26 + s56 * t27 + t31 + t33 + t36
      t53 = 1.0_dp - epinv
      t56 = t66 * t43
      t74 = intHLs250000x0113D4eps1()
      t79 = -t2 + epinv
      t80 = intHLs250000x0111D2eps1()
      t5 = t65 * (t5 + t50) * t7
      t81 = t19 * epinv
      t82 = t66 * t42
      t83 = intHLs250000x0122D4eps1()
      t84 = s26 + s56
      t85 = s25 + s26 + s45 + s46 + s56
      t13 = t69 * t14 + t40 * (-t42 * (epinv * (-s56 * t73 + t19 * (t66 
     #* ((s12 + s13 + s16) * s25 + (s12 + 2.0_dp * s13 + s16 + s25) * s2
     #6 + (-s13 - t55) * s34 + (2.0_dp * s13 + s12 + s16 + s26 - s34 - s
     #45 - s46) * s56 + t13 + t31 - t61 + t75) + t61) - t2 * (s26 * t72 
     #/ 2.0_dp - t30 + t31 + t36) - t34) + t7 * (t11 * t38 + t64 * t32) 
     #- 4.0_dp * s13 * (t19 * t84 + t23 + t25) + 4.0_dp * t85 * s16 + 4.
     #0_dp * s26 * t55 - 4.0_dp * t85 * s34 - 8.0_dp * s12 * t6 - 8.0_dp
     # * s25 * t8 - 12.0_dp * s26 * t8 + 12.0_dp * (s25 + s45 + s46) * s
     #56 + 12.0_dp * t9 + 16.0_dp * s26 * s56) + t68 * t53)
      t75 = s25 * t1
      t8 = s25 * ((t20 * (-t57 * t56 * t14 * t53 - t40 * (t42 * (-t19 * 
     #(epinv * (t2 * (s13 * t48 + (2.0_dp * s13 + t21 - s25 - s16 - s26 
     #+ s34 - s56) * s25 + (s13 + t16 - t7) * s26 + (-t2 * (s25 - s13 - 
     #s12 + s45 + s46) - t24 - t48) * s56 + t8 * (t25 + t24 + t48) + t10
     #) - t39 - 4.0_dp * t9 + t66 * t35) + t66 * t7 * t32) + 4.0_dp * t5
     #1 - 12.0_dp * t9 - t12 * t38) + t68)) + t83 * t13) * propW34 + t70
     # * t76 * (t42 * (-t19 * (t66 * t28 + t61) + t2 * (-t30 + t36 + t31
     #) + s26 * t72 + s56 * t73 + t34) - t68))
      t10 = propW34 * s25
      t12 = t38 * t84
      t13 = t2 * t23
      t16 = (propZ25 * t29 - 4.0_dp) * s25 - t12 - t13
      t21 = t2 * t84 + t50
      t24 = t46 * t21
      t25 = -t79
      t30 = propZ25 * t7
      t34 = t2 * t40
      t35 = t34 * t42
      t48 = (-propZ25 * t21 * t19 + 4.0_dp) * s25 + t12 + t13
      t50 = t67 * t29
      t61 = (t30 * t38 + 8.0_dp) * s25
      t72 = 12.0_dp * t84
      t73 = t74 * t25
      t84 = t73 - t47
      t85 = t40 * t70
      t25 = t80 * t25
      t86 = propW16 * t38 * t7
      t87 = -t86 - t66 + 4.0_dp
      t88 = -t2 * propW16 * t7 + 1.0_dp
      t89 = -t88
      t90 = -t87
      t91 = 1.0_dp - t66
      t92 = t86 * t91 + t64 - 4.0_dp
      t93 = t2 + t66
      t94 = -t19 * propW16 * t93 * t7 + 4.0_dp
      t86 = t86 * t93
      t95 = t67 * t53
      t96 = t65 * t2 + 1.0_dp
      t97 = t65 * t38
      t98 = t97 + t2 + t66
      t99 = t40 * t98
      t64 = t97 * t91 + t2 - t64
      t65 = t65 * t93 + t66
      t91 = t65 * t38
      t93 = t69 * t96
      t64 = t7 * (-t19 * t42 * ((t75 + t25) * (t46 * t96 + t99) - t99 * 
     #t77) - t2 * s25 * (t66 * (t43 * t96 * t57 + t41) + t35 * (t19 * t6
     #5 + t2)) * t84 + (t20 * (t57 * t56 * t96 * t53 + t40 * (-t42 * (t8
     #1 * t98 - t91 - 4.0_dp) + t67)) - t83 * (t93 + t40 * (t42 * (epinv
     # * t64 + t91 + 4.0_dp) + t95)) + t85 * (-t42 * t64 + t67)) * s25 +
     # t93 * t77) * struc7MP
      t1 = (t2 * s25 * ((t20 * (-t57 * t56 * (epinv * t89 + t88) + t40 *
     # (t42 * (-t81 * t90 + t86 - 8.0_dp) + t67)) + t83 * (-t69 * t89 - 
     #t40 * (t42 * (epinv * t92 + t86 - 8.0_dp) + t95))) * propW34 + t70
     # * t76 * (-t42 * t92 + t67)) + t38 * t42 * propW34 * (t77 * (t40 *
     # t90 + t46 * t89) + t25 * (t40 * t87 + t46 * t88)) + 4.0_dp * t10 
     #* (t47 * (-t35 * t94 + t66 * (t43 * t89 * t57 + t41)) - t73 * (-t3
     #5 * t94 + t66 * (-t43 * t88 * t57 + t41))) - t45 * t1 * propW34 * 
     #(t82 * t38 * t89 + t34 * t90)) * t6 * (struc3MP + struc5MP + struc
     #8MP)
      t1 = (t19 * (t77 * propW34 * (-t46 * t28 + t4 * t40) + t75 * propW
     #34 * (-t40 * (-t2 * ((-t22 - s26) * s26 + t31 + t33 + t36) + t9 * 
     #t38 - t2 * t27 * s56 + (-t2 * t17 - t37) * s25 + t39) + t46 * t14)
     # + propW34 * t80 * (t40 * (epinv * t4 + t2 * ((t37 + t54 + s16) * 
     #s25 + (t59 + t71 + t18 + s26) * s26 - (t49 + t55) * s34 + (-t3 - t
     #44 + t59 + t58 + t23) * s56) - t38 * (t9 + t5) + 4.0_dp * t31 + 4.
     #0_dp * t36) + t82 * (t38 * (s25 * t15 + (t49 + s12) * s26 + (t60 -
     # s25 - s26 + s12 - s56) * s56 + t52 + t62 - t5) - t81 * t28))) * t
     #42 + t2 * t10 * (t47 * (t69 * t26 - t40 * (-t42 * t63 + t68)) + t7
     #4 * (t56 * t26 * (epinv * t57 - 18.0_dp) - t40 * (-t42 * (epinv * 
     #t63 - t7 * (t66 * t32 * t38 + 12.0_dp * t11) + 8.0_dp * t51 - 24.0
     #_dp * t9) + t68 * t79))) + t8) * struc1MP + t1
      t3 = 0.1e1_dp / t6
      t4 = 0.1e1_dp / t7
      t5 = 0.1e1_dp / 0.9e1_dp
      result = t5 * (t2 * propW34 * ((-t19 * t42 * ((t25 - t77) * (t16 *
     # t40 - t24) - t75 * (-t16 * t40 + t24)) - t2 * s25 * (t66 * (-t43 
     #* t21 * t57 + t29 * t41) - t35 * ((t30 * t19 + 4.0_dp) * s25 + t12
     # + t13)) * t84 + (t20 * (-t57 * t56 * t21 * t53 + t40 * (t42 * (-t
     #81 * t16 - t61 - t72 - t78) + t50)) + t83 * (t69 * t21 - t40 * (t4
     #2 * (-epinv * t48 - t61 - t72 - t78) + t50 * t53)) + t85 * (t42 * 
     #t48 + t50)) * s25) * struc2MP + t64) + t1) / ecossin ** 2.0_dp / g
     #w ** 2.0_dp / s25 * t3 * t4

           ampNonresonantHeavyReC7MP = result
       end function ampNonresonantHeavyReC7MP

       function ampNonresonantHeavyReC7PM()
           implicit none
           complex(dp) :: ampNonresonantHeavyReC7PM
           complex(dp) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           complex(dp) ::  t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30
           complex(dp) ::  t31,t32,t33,t34,t35,t36,t37,t38,t39,t4,t40,t41
           complex(dp) ::  t42,t43,t44,t45,t46,t47,t48,t49,t5,t50,t51,t52
           complex(dp) ::  t53,t54,t55,t56,t57,t58,t59,t6,t60,t61,t62,t63
           complex(dp) ::  t64,t65,t66,t67,t68,t69,t7,t70,t71,t72,t73,t74
           complex(dp) ::  t75,t76,t77,t78,t79,t8,t80,t81,t82,t83,t84,t85
           complex(dp) ::  t86,t87,t88,t89,t9,t90,t91,t92

           complex(dp) :: result

      t1 = intHLs250000x0122D4eps0()
      t2 = s25 + s26 + s56
      t3 = s16 + s26 - s34 + s56
      t4 = s13 - s45 - s46
      t5 = propW16 * t2
      t6 = gb ** 2.0_dp
      t7 = t6 ** 2.0_dp
      t8 = gw ** 2.0_dp
      t9 = t8 ** 2.0_dp
      t10 = t6 * t8
      t11 = t10 * s25 * t3 * (propZ25 * t4 + t5 * (-1.0_dp + (-s12 - s16
     # - s26 + s34) * propZ25))
      t12 = s12 + s16 + s26 - s34
      t13 = t10 * t5 * t3
      t14 = t2 + t3
      t15 = -t2 - t3
      t16 = propW16 * s25 * t3
      t17 = propZ25 * s25
      t18 = t17 * t9
      t19 = (s12 + s13 + s16 + s26 - s34 - s45 - s46) * t2
      t20 = (-s13 + s45 + s46) * t3
      t21 = t17 * t7 * (t20 + t19)
      t22 = t18 * t5 * t3 * t12
      t23 = 3.0_dp
      t24 = s25 ** 2.0_dp
      t25 = t13 * t23 * propZ25 * t24
      t26 = intHLs250000x0122D4eps1()
      t27 = propZ25 * t14
      t28 = t5 * t3
      t29 = t3 * t24
      t30 = t18 * epinv
      t20 = -24.0_dp * t13 * (t17 * (-s12 - s26) + s12 + s26) + 16.0_dp 
     #* t10 * t2 * (s12 + s13 - s45 - s46 - s56) - 12.0_dp * t10 * ((t3 
     #* (-propZ25 * s26 + t5) + propZ25 * s12 * t2 - t27 * s45 + t27 * s
     #13 - t27 * s56 - t27 * s46) * s25 + t28 * t12 * epinv - t29 * prop
     #Z25 * (t5 + 1.0_dp)) - t21 * epinv + 9.0_dp * t30 * (s13 * t3 + s4
     #5 * t15 + s46 * t15 + t2 * (-t16 + s12 + s13 + s16 + s26 - s34)) +
     # 8.0_dp * t10 * (t3 * (s13 - s25 - s26 - s45 - s46 - s56) + t19 * 
     #epinv) - 4.0_dp * t6 * (t20 * t8 * epinv + ((t2 * (-s12 - s13 + s4
     #5 + s46 + s56) + t3 * (s13 - s26 - s45 - s46 - s56)) * s25 - t29) 
     #* propZ25 * t6) - 18.0_dp * t22 * epinv + 6.0_dp * t11 * epinv - t
     #25 * epinv
      t27 = -s16 + s25 + s34
      t29 = -s45 - s46
      t31 = (-s13 + s25 + s26) * s34
      t32 = s13 * s16
      t33 = (s12 + s13 + s16) * s25 + (s12 + s16 + s25 + s26) * s26 + (s
     #12 + s16 + s26 - s34) * s56 + t29 * t27 - t31 - t32
      t34 = 2.0_dp
      t35 = t34 * s26
      t36 = s12 + s26
      t37 = t34 * t36
      t38 = s56 ** 2.0_dp
      t39 = t29 * (t35 + s25 + s56 + s16 - s34)
      t40 = t28 * (t37 + s25)
      t41 = (-s16 - s26 + s34) * s13 + (-s13 - s12 + s16 + s25 + s26 - s
     #34 + s56) * s25 - s26 * (s12 - s26 - s56 + s13 - s16 + s34) + (-2.
     #0_dp * s13 + t35 + s25 - s12 + s56 + s16 - s34 + s45 + s46) * s56 
     #- t24 + t38 - t39 + t40
      t42 = t17 - t34
      t43 = s56 + s16 - s34
      t44 = t23 * s26
      t45 = s12 + s13
      t46 = t34 * t45
      t47 = -t46 + t3
      t48 = t34 * s12
      t49 = t34 * s25
      t50 = (-t34 * (s25 - s12 - s13 + s45 + s46) - t43 - t44) * s56
      t24 = t34 * (s13 * t43 + (2.0_dp * s13 + t48 - s25 - s26 - s56 - s
     #16 + s34) * s25 + (s13 - t47) * s26 + t29 * (t44 + t49 + t43) + t2
     #4 + t50) - 4.0_dp * t38 - t40 * t23
      t43 = t17 * t6
      t51 = t43 * t33
      t52 = t18 * t23
      t53 = intHLs250000x0111D2eps0()
      t54 = t34 * s13
      t55 = s25 + s16 - s34
      t56 = s25 - s13 + s45 + s46
      t57 = t29 * (t35 + t55)
      t58 = t34 * t38 - (s12 + s13 - s16) * s25 - (t54 - s25 + s12 - s16
     # - s26) * s26 + (t34 * t56 - s12 + s16 - s34 + t44) * s56 - t31 - 
     #t32 + t40 - t57
      t59 = s16 - s34
      t60 = t23 * s25
      t61 = 4.0_dp * t59
      t62 = 6.0_dp * t36
      t63 = s16 + s25
      t64 = t23 * s13
      t65 = s13 - s45 - s46
      t66 = s25 + s26
      t67 = 4.0_dp * s26
      t68 = t59 * s13
      t29 = t29 * (t44 + t49 + t59)
      t44 = t66 * s34
      t69 = 4.0_dp * s12
      t62 = -4.0_dp * (t46 + s16) * s25 - 4.0_dp * (t64 + t48 + t63 + s2
     #6) * s26 + 4.0_dp * (s13 + s25 + s26) * s34 - 4.0_dp * (t23 * t65 
     #- s56 - t27 + t48) * s56 - 4.0_dp * t32 - 4.0_dp * t29 - t17 * (t3
     #4 * (-(s16 * t23 + t46) * s25 + (-t23 * (s13 + s16 - s34 - s45 - s
     #46) - s25 - t48 - t67 - s56) * s56 - t68 - t29) - 6.0_dp * (s13 + 
     #s16 + s25 + s26) * s26 + 6.0_dp * t44 + t28 * (8.0_dp * t59 + t60 
     #+ t62) - t69 * s26) + t28 * t34 * (t62 + t61 + t60)
      t70 = intHLs250000x0113D4eps0()
      t71 = 4.0_dp * t45
      t72 = 6.0_dp * t65
      t46 = -t46 + s16
      t73 = t46 * s25
      t36 = t28 * (4.0_dp * t36 + t59 + t49)
      t74 = 10.0_dp
      t75 = 18.0_dp * t38
      t76 = (8.0_dp * s12 + 12.0_dp * s13) * s26
      t77 = 6.0_dp * t28 * (t37 + t27)
      t78 = t34 * ((-5.0_dp * s16 + t71) * s25 + (-9.0_dp * s25 - 14.0_d
     #p * s26 - 5.0_dp * t59 + t69 + t72) * s56)
      t79 = -4.0_dp * t29 - 4.0_dp * t68
      t44 = t74 * ((-t63 - s26) * s26 + t44)
      t74 = t17 * t23
      t80 = 9.0_dp * t9 * t2 * t3 * (-propW16 * t59 + 1.0_dp)
      t33 = t7 * t33 * t34
      t81 = intHLs250000x0112D2eps1()
      t82 = t23 * epinv
      t83 = t82 - 4.0_dp
      t64 = t64 + t48 - t63
      t56 = -t23 * t56 + t48 - t59 - t67
      t67 = 24.0_dp * t38
      t40 = 12.0_dp * t17 * t41 - 12.0_dp * t40
      t84 = t82 * t42
      t85 = intHLs250000x0111D2eps1()
      t86 = -t34 + epinv
      t28 = t28 * (t37 + t55)
      t37 = t34 * t42
      t87 = intHLs250000x0113D4eps1()
      t88 = s25 + s26 + s56 + s16 - s34
      t89 = s25 * intHLs250000x0112D2eps0()
      t90 = propW34 * s25
      t33 = t90 * (t70 * (t17 * (t33 - t80) - t10 * (-t74 * (t23 * t38 +
     # t34 * (-t32 - t57) + (-4.0_dp * s13 + t63 - t48 + s26) * s26 + (t
     #54 - t66) * s34 + (t59 + 4.0_dp * s26 - 4.0_dp * s13 + 4.0_dp * s4
     #5 + 4.0_dp * s46 - t48 + t60) * s56 + t73 + t36) - t44 + t79 - t78
     # + t77 - t76 + t75)) + t87 * (t10 * (epinv * (-t74 * (-t34 * (-s13
     # * t88 + t38 - t39) - (-t48 + t88 - s25) * s25 - s26 * t47 - t36 +
     # t50) + t44 - t79 + t78 - t77 + t76 - t75) - t40 + t67 + 8.0_dp * 
     #(-t64 + s26) * s26 - 8.0_dp * s56 * t56 - 8.0_dp * t29 - 8.0_dp * 
     #t31 - 8.0_dp * t32 + 8.0_dp * t73) + t17 * (-epinv * t80 + t33 * t
     #86)))
      t4 = t90 * (t1 * (-12.0_dp * t13 * t12 + 6.0_dp * t11 - 9.0_dp * t
     #18 * (s13 * t15 + s45 * t3 + s46 * t14 + t2 * (t16 - s12 - s16 - s
     #26 + s34 + s45)) + 8.0_dp * t19 * t10 + 4.0_dp * t10 * t3 * t4 - 1
     #8.0_dp * t22 - t21 - t25) + t26 * t20 + t81 * (-9.0_dp * t30 * t41
     # + t6 * (-t8 * (-8.0_dp * t46 * s25 + 8.0_dp * (t64 - s26) * s26 +
     # 8.0_dp * s56 * t56 + 8.0_dp * t29 + 8.0_dp * t31 + 8.0_dp * t32 -
     # t67 + t40 + t84 * t24) + t51 * t83)))
      t11 = -t86
      t12 = s26 + s56
      t13 = t12 * t34 + t55
      t14 = t12 * t23 + t49 + t59
      t15 = t7 * t27
      t16 = t9 * t13
      t19 = t37 * t10
      t20 = -t19 * t14 + t17 * (t16 * t23 + t15)
      t21 = t34 * t59
      t22 = 6.0_dp * t12
      t25 = s25 * (propZ25 * t3 * t23 + 4.0_dp) + t21 + t22
      t16 = 9.0_dp * t16
      t10 = t10 * t34
      t13 = propZ25 * t13
      t30 = (6.0_dp * t13 - 8.0_dp) * s25
      t12 = 12.0_dp * t12
      t31 = t16 * epinv
      t36 = -4.0_dp + epinv
      t39 = t85 * t11
      t40 = t89 + t39
      t11 = t87 * t11
      t44 = t11 - t70
      t46 = t6 * s25
      t47 = t23 * propW16 * t3
      t48 = -t47 + t34
      t49 = t34 * propW16 * t3 - 1.0_dp
      t50 = t9 * t49
      t55 = -t48
      t56 = 1.0_dp - t17
      t64 = 6.0_dp * propW16 * t56 * t3 + t74 - 4.0_dp
      t67 = t17 + t34
      t47 = t47 * t67 - 4.0_dp
      t73 = 9.0_dp * t50
      t75 = 12.0_dp * propW16 * t3
      t76 = 6.0_dp * t17
      t77 = t76 * t49
      t78 = t73 * epinv
      t79 = t7 * t83
      t80 = t7 * t36
      t87 = t5 * t34 + 1.0_dp
      t88 = t5 * t23 + 1.0_dp
      t90 = t9 * t87
      t91 = t17 * (-t90 * t23 + t7)
      t67 = t23 * (t5 * t67 + t17) + t34
      t90 = 9.0_dp * t90
      t92 = 12.0_dp * t5
      t76 = t76 * t87
      t87 = t90 * epinv
      t9 = (-t23 * ((t19 * t55 + t17 * (-t9 * t49 * t23 + t7)) * (-t39 +
     # t53) - t89 * (-t19 * t48 + t17 * (-t50 * t23 + t7))) + 8.0_dp * t
     #46 * (-t70 * (-t64 * t8 + t43) + t11 * (-t64 * t8 + t43)) - (t1 * 
     #(t10 * t47 + t17 * (t73 + t7)) + t26 * (-t10 * (-epinv * t47 - t75
     # + t77 + 8.0_dp) + t17 * (t80 + t78)) - t81 * (-t10 * (-t84 * t55 
     #- t75 + t77 + 8.0_dp) + t17 * (t79 - t78))) * s25) * t2 * (struc3P
     #M + struc5PM + struc8PM)
      t1 = propW34 * (struc2PM * (t23 * (t20 * t40 - t53 * t20) + 8.0_dp
     # * t46 * (t8 * (s25 * (-t13 * t23 + 4.0_dp) + t21 + t22) + t43 * t
     #27) * t44 - (t1 * (-t10 * t25 + t17 * (-t16 + t15)) + t26 * (-t10 
     #* (t25 * epinv + t12 - t30 + t61) + t17 * (t15 * t36 - t31)) - t81
     # * (t10 * (-t84 * t14 - t12 + t30 - t61) + t17 * (t15 * t83 + t31)
     #)) * s25) - t3 * (t23 * (-t40 * (t19 * t88 + t91) + t53 * (t10 * t
     #42 * t88 + t91)) - 8.0_dp * t46 * (t8 * (-6.0_dp * t5 * t56 - t34 
     #+ t74) + t43) * t44 - (-t1 * (t10 * t67 + t17 * (t90 + t7)) - t26 
     #* (t17 * (t80 + t87) + t10 * (epinv * t67 - t76 + t92 + 4.0_dp)) +
     # t81 * (t17 * (t79 - t87) + t10 * (t84 * t88 - t76 + t92 + 4.0_dp)
     #)) * s25) * struc7PM + t9)
      t2 = 0.1e1_dp / t2
      t3 = 0.1e1_dp / t3
      t5 = 0.1e1_dp / 0.18e2_dp
      result = t5 * (t34 * t1 + (-t23 * (-(t53 * (t52 * t58 + t6 * (t62 
     #* t8 - t51)) + t89 * (-t52 * t41 + t6 * (-t24 * t42 * t8 + t51))) 
     #* propW34 + t85 * (-t18 * propW34 * (6.0_dp * s25 * t45 + 6.0_dp *
     # (t54 + s12) * s26 + 6.0_dp * (t34 * t65 + s12 - s25 - s26 - s56) 
     #* s56 + 6.0_dp * t57 + 6.0_dp * t68 - 6.0_dp * t28 + t82 * t58) + 
     #t6 * (-t8 * propW34 * (epinv * t62 + t37 * (t23 * (t28 + t38) + t3
     #4 * (-t32 - t29) - (t71 + s16) * s25 - (6.0_dp * s13 + t63 + t69 +
     # s26) * s26 + (t54 + t66) * s34 + (-t72 + t35 - t69 + t60 - t59) *
     # s56)) + t51 * t86 * propW34))) - 4.0_dp * t33 + t4) * struc1PM) /
     # ecossin ** 2.0_dp / gw ** 2.0_dp / s25 * t2 * t3

           ampNonresonantHeavyReC7PM = result
       end function ampNonresonantHeavyReC7PM

       function ampNonresonantHeavyReC7PP()
           implicit none
           complex(dp) :: ampNonresonantHeavyReC7PP

           complex(dp) :: result

      result = 0.0_dp

           ampNonresonantHeavyReC7PP = result
       end function ampNonresonantHeavyReC7PP

       function ampNonresonantLightFullImC4MM()
           implicit none
           complex(dp) :: ampNonresonantLightFullImC4MM

           complex(dp) :: result

      result = 0.0_dp

           ampNonresonantLightFullImC4MM = result
       end function ampNonresonantLightFullImC4MM

       function ampNonresonantLightFullImC4MP()
           implicit none
           complex(dp) :: ampNonresonantLightFullImC4MP
           complex(dp) ::  t1,t10,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109
           complex(dp) ::  t11,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t12
           complex(dp) ::  t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t13,t130
           complex(dp) ::  t131,t132,t133,t134,t135,t136,t137,t138,t139,t14,t140,t141
           complex(dp) ::  t142,t143,t144,t145,t146,t147,t148,t149,t15,t150,t151,t152
           complex(dp) ::  t153,t154,t155,t156,t157,t158,t159,t16,t160,t161,t162,t163
           complex(dp) ::  t164,t165,t166,t167,t168,t169,t17,t170,t171,t172,t173,t174
           complex(dp) ::  t175,t176,t177,t178,t179,t18,t180,t181,t182,t183,t184,t185
           complex(dp) ::  t186,t187,t188,t189,t19,t190,t191,t192,t193,t194,t195,t196
           complex(dp) ::  t197,t198,t199,t2,t20,t200,t201,t202,t203,t204,t205,t206
           complex(dp) ::  t207,t208,t209,t21,t210,t211,t212,t213,t214,t215,t216,t217
           complex(dp) ::  t218,t219,t22,t220,t221,t222,t223,t224,t225,t226,t227,t228
           complex(dp) ::  t229,t23,t230,t231,t232,t233,t234,t235,t236,t237,t238,t239
           complex(dp) ::  t24,t240,t241,t242,t243,t244,t245,t246,t247,t248,t249,t25
           complex(dp) ::  t250,t251,t252,t253,t254,t255,t256,t257,t258,t259,t26,t260
           complex(dp) ::  t261,t262,t263,t264,t265,t266,t267,t268,t269,t27,t270,t271
           complex(dp) ::  t272,t273,t274,t275,t276,t277,t278,t279,t28,t280,t281,t282
           complex(dp) ::  t283,t284,t285,t286,t287,t288,t289,t29,t290,t291,t292,t293
           complex(dp) ::  t294,t295,t296,t297,t298,t299,t3,t30,t300,t301,t302,t303
           complex(dp) ::  t304,t305,t306,t307,t308,t309,t31,t310,t311,t312,t313,t314
           complex(dp) ::  t315,t316,t317,t318,t319,t32,t320,t321,t322,t323,t324,t325
           complex(dp) ::  t326,t327,t328,t329,t33,t330,t331,t332,t333,t334,t335,t336
           complex(dp) ::  t337,t338,t339,t34,t340,t341,t342,t343,t344,t345,t346,t347
           complex(dp) ::  t348,t349,t35,t350,t351,t352,t353,t354,t355,t356,t357,t358
           complex(dp) ::  t359,t36,t360,t361,t362,t363,t364,t365,t366,t367,t368,t369
           complex(dp) ::  t37,t370,t371,t372,t373,t374,t375,t376,t377,t378,t379,t38
           complex(dp) ::  t380,t381,t382,t383,t384,t385,t386,t387,t388,t389,t39,t390
           complex(dp) ::  t391,t392,t393,t394,t395,t396,t397,t398,t399,t4,t40,t400
           complex(dp) ::  t401,t402,t403,t404,t405,t406,t407,t408,t409,t41,t410,t411
           complex(dp) ::  t412,t413,t414,t415,t416,t417,t418,t419,t42,t420,t421,t422
           complex(dp) ::  t423,t424,t425,t426,t427,t428,t429,t43,t430,t431,t432,t433
           complex(dp) ::  t434,t435,t436,t437,t438,t439,t44,t440,t441,t442,t443,t444
           complex(dp) ::  t445,t446,t447,t448,t449,t45,t450,t451,t452,t453,t454,t455
           complex(dp) ::  t456,t457,t458,t459,t46,t460,t461,t462,t463,t464,t465,t466
           complex(dp) ::  t467,t468,t469,t47,t470,t471,t472,t473,t474,t475,t476,t477
           complex(dp) ::  t478,t479,t48,t480,t481,t482,t483,t484,t485,t486,t487,t488
           complex(dp) ::  t489,t49,t490,t491,t492,t493,t494,t495,t496,t497,t498,t499
           complex(dp) ::  t5,t50,t500,t51,t52,t53,t54,t55,t56,t57,t58,t59
           complex(dp) ::  t6,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t7
           complex(dp) ::  t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t8,t80
           complex(dp) ::  t81,t82,t83,t84,t85,t86,t87,t88,t89,t9,t90,t91
           complex(dp) ::  t92,t93,t94,t95,t96,t97,t98,t99

           complex(dp) :: result

      t1 = intHLs160000x0112D2eps0()
      t2 = intHLs16s25s26s34s56x1111D2eps0()
      t3 = intHs160000x0112D2eps0()
      t4 = intHs160s26s34s56x1021D2eps1()
      t5 = intHs160s26s34s56x1022D4eps0()
      t6 = intHs16s25s26s34s56x1112D4eps0()
      t7 = intHs16s25s26s34s56x1121D4eps0()
      t8 = intHs16s25s26s34s56x1211D4eps0()
      t9 = intHs16s25s26s34s56x1131D4eps0()
      t10 = 4.0_dp
      t11 = propZ25 * s25
      t12 = t11 - t10
      t13 = intHs16s25s26s34s56x1221D4eps0()
      t14 = intHs16s25s26s34s56x1311D4eps0()
      t15 = intHs16s25s26s34s56x1311D4eps1()
      t16 = 2.0_dp
      t17 = -t16 + epinv
      t18 = 3.0_dp
      t19 = gb ** 2.0_dp
      t20 = gw ** 2.0_dp
      t21 = t18 * s25
      t22 = t21 * t20 * propZ25
      t23 = -t19 * t12
      t24 = t22 + t23
      t25 = 1.0_dp - epinv
      t26 = t19 * t12
      t27 = -t22 + t26
      t28 = intHLs16s25s26s34s56x1112D4eps0()
      t29 = t11 + t16
      t30 = t19 * t29
      t31 = t30 + t22
      t32 = intHs16s25s26s34s56x1310D4eps0()
      t33 = s16 + s26 - s34 + s56
      t34 = intHLs16s25s26s34s56x1222D6eps1()
      t35 = s12 + s25
      t36 = t16 * t35
      t37 = t10 * s26
      t38 = t37 + t36 + s16
      t39 = intHLs16s25s26s34s56x1211D2eps1()
      t40 = intHs16s25s26s34s56x1222D6eps1()
      t41 = 6.0_dp
      t42 = 5.0_dp
      t43 = t18 * s26
      t44 = t42 * s25
      t45 = t41 * s12
      t46 = t45 + s34 - s56 + t43 + t44
      t47 = intHLs16s25s26s34s56x1121D2eps1()
      t48 = s25 + s26 + s56
      t49 = intHLs16s25s26s34s56x1112D2eps1()
      t50 = -s13 + s26 + s45 + s46 + s56
      t51 = t16 * t50 + s25
      t52 = intHs16s25s26s34s56x1213D6eps1()
      t53 = s12 + s25 + s26
      t54 = t16 * t53 + s34
      t55 = s25 + s26 - s56
      t56 = intHLs16s25s26s34s56x1123D6eps1()
      t57 = t16 * s25
      t58 = t43 + t57 - s56
      t59 = intHLs16s25s26s34s56x1213D6eps1()
      t60 = t43 + t36 - s56
      t61 = intHs16s25s26s34s56x1132D6eps1()
      t62 = t10 * s12
      t63 = s26 + s34 - s56 + t57 + t62
      t64 = intHs16s25s26s34s56x1123D6eps1()
      t65 = t16 * (s12 + s34) + t55
      t66 = intHs16s25s26s34s56x1110D2eps1()
      t67 = intHs16s25s26s34s56x1220D4eps0()
      t68 = intHs16s25s26s34s56x1130D4eps0()
      t69 = t16 * s12
      t70 = t69 + s25
      t71 = intHs16s25s26s34s56x1120D2eps1()
      t72 = t25 * s25
      t73 = intHs160s26s34s56x1020D2eps1()
      t74 = intHs16s25s26s34s56x1321D6eps1()
      t75 = t69 + t21
      t76 = intHs16s25s26s34s56x1231D6eps1()
      t77 = t21 + t62
      t78 = intHLs16s25s26s34s56x1132D6eps1()
      t79 = t16 * s26
      t80 = t79 + s25
      t81 = intHLs16s25s26s34s56x1312D6eps1()
      t82 = s25 + s16
      t83 = t16 * (s12 + s26)
      t84 = t82 + t83
      t85 = intHs16s25s26s34s56x1312D6eps1()
      t86 = t21 + t83
      t87 = intHs160s26s34s56x1031D4eps0()
      t88 = 0.1e1_dp / s25
      t89 = t19 * (-propZ25 * t16 + 8.0_dp * t88) + t41 * t20 * propZ25
      t90 = intHs16s25s26s34s56x1211D2eps1()
      t91 = intHs16s25s26s34s56x1222D6eps0()
      t92 = intHLs16s25s26s34s56x1121D2eps0()
      t93 = intHs16s25s26s34s56x1120D2eps0()
      t94 = intHs16s25s26s34s56x1210D2eps0()
      t95 = intHs16s25s26s34s56x1210D2eps1()
      t96 = intHs16s25s26s34s56x1122D4eps1()
      t97 = t17 * s25
      t98 = t16 * epinv
      t99 = intHLs16s25s26s34s56x1112D2eps0()
      t100 = intHLs16s25s26s34s56x1112D4eps1()
      t101 = -t18 + t98
      t102 = intHs16s25s26s34s56x1122D4eps0()
      t103 = intHs16s25s26s34s56x1212D4eps0()
      t104 = intHLs16s25s26s34s56x1222D6eps0()
      t105 = intHs16s25s26s34s56x1212D4eps1()
      t106 = intHLs16s25s26s34s56x1211D2eps0()
      t107 = intHLs16s25s26s34s56x1212D4eps1()
      t108 = s13 - s45 - s46
      t109 = t16 * t108
      t110 = t18 * s56
      t111 = s25 + s26 - s34 + t109 + t62 - t110
      t112 = epinv * t111
      t113 = intHLs16s25s26s34s56x1212D4eps0()
      t114 = intHs160s26s34s56x1022D4eps1()
      t115 = intHs160000x0112D2eps1()
      t116 = intHs16s25s26s34s56x1121D2eps1()
      t117 = 1.0_dp + epinv
      t118 = t69 * t117
      t119 = intHs16s25s26s34s56x1121D4eps1()
      t120 = intHs16s25s26s34s56x1211D4eps1()
      t121 = intHs16s25s26s34s56x1112D4eps1()
      t122 = intHLs160000x0111D0eps0()
      t123 = intHLs160000x0111D0eps1()
      t124 = intHLs160000x0112D2eps1()
      t125 = intHs16s25s26s34s56x1220D4eps1()
      t126 = intHs16s25s26s34s56x1130D4eps1()
      t127 = intHs16s25s26s34s56x1411D6eps0()
      t128 = intHLs16s25s26s34s56x1122D4eps1()
      t129 = t16 * s56
      t130 = -t108 + t129
      t131 = intHs16s25s26s34s56x1121D2eps0()
      t132 = intHs16s25s26s34s56x1131D4eps1()
      t133 = t97 - t62
      t134 = intHs16s25s26s34s56x1221D4eps1()
      t135 = t97 - t69
      t136 = intHs16s25s26s34s56x1141D6eps0()
      t137 = intHs16s25s26s34s56x1321D6eps0()
      t138 = intHs16s25s26s34s56x1231D6eps0()
      t139 = intHLs16s25s26s34s56x1132D6eps0()
      t140 = intHLs16s25s26s34s56x1312D6eps0()
      t141 = intHs16s25s26s34s56x1312D6eps0()
      t142 = intHs16s25s26s34s56x1113D4eps1()
      t143 = intHs16s25s26s34s56x1114D6eps0()
      t144 = intHs16s25s26s34s56x1213D6eps0()
      t145 = intHLs16s25s26s34s56x1113D4eps0()
      t146 = -s26 + t109 - t110
      t147 = intHLs16s25s26s34s56x1122D4eps0()
      t148 = intHLs16s25s26s34s56x1114D6eps0()
      t149 = intHLs16s25s26s34s56x1123D6eps0()
      t150 = intHLs16s25s26s34s56x1213D6eps0()
      t151 = intHs16s25s26s34s56x1132D6eps0()
      t152 = intHs16s25s26s34s56x1123D6eps0()
      t153 = intHs16s25s26s34s56x1310D4eps1()
      t154 = intHs16s25s26s34s56x1110D2eps0()
      t155 = intHs160s26s34s56x1020D2eps0()
      t156 = intHLs16s25s26s34s56x1113D4eps1()
      t157 = 0.1e1_dp / t33
      t158 = t66 * t17
      t159 = t25 * t40
      t160 = -t159 + t91
      t161 = t25 * t34
      t162 = t161 - t104
      t163 = epinv * t49 + t99
      t164 = t96 * (t98 * s26 - s26 - s34 + s56 - t62 + t97)
      t165 = (t71 * (t72 + t69) + t158) * t157
      t166 = (t102 + t103) * t80
      t167 = t100 * t101
      t168 = (epinv * t47 + t92) * t48
      t169 = s34 * (epinv * t39 + t106)
      t170 = t25 * (t73 * t88 + t95) * t157
      t171 = t87 * t89
      t172 = t155 * t88 + t93 + t94
      t173 = t119 + t120
      t174 = t7 + t8
      t175 = t1 * t29
      t176 = epinv * t115
      t177 = t154 * t157
      t178 = t114 * t25
      t179 = t116 * (t118 + s25)
      t180 = t101 * t121
      t111 = t111 * t113
      t181 = (-t123 + t124) * epinv
      t182 = t25 * t125
      t183 = t25 * t126
      t184 = -t68 + t183
      t185 = t52 * t54
      t186 = t61 * t63
      t187 = t64 * t65
      t188 = t74 * t75
      t189 = t76 * t77
      t190 = t85 * t86
      t191 = t184 * t70
      t192 = (-t67 + t182) * t35
      t193 = t15 * t17
      t194 = t78 * t80
      t195 = t81 * t84
      t196 = t6 * t88
      t197 = t88 * t25
      t198 = t88 * t31
      t199 = t145 * t146
      t200 = -t147 * t130
      t201 = t80 * t139
      t84 = t84 * t140
      t202 = t54 * t144
      t203 = t63 * t151
      t204 = t77 * t138
      t205 = t86 * t141
      t206 = t152 * t65
      t75 = t137 * t75
      t207 = s34 * t142
      t208 = t132 * t133
      t209 = intHs160s26s34s56x1031D4eps1() * t25
      t210 = t131 * s12
      t211 = t2 * t29
      t212 = t211 * t19
      t213 = t25 * t153
      t58 = t88 * (t24 * (-t134 * t135 - t208 + t209 - t210) - t212) + t
     #24 * (t157 * (t88 * (t192 + t191) - t32) - t193) + t197 * (t188 + 
     #t189 + t190 + t185 + t186 + t187) * t27 + t198 * (-t25 * (t56 * t5
     #8 + t59 * t60 + t194 + t195) + t28) + t26 * (t9 + t13 + t14 - t196
     #) + t213 * t24 * t157 - (t207 + t206 + t75 + t205 + t202 + t203 + 
     #t204) * t88 * t27 + t198 * (t149 * t58 + t150 * t60 + t156 * (-epi
     #nv * t146 + s25 + s26 - s56) - t199 - t200 + t201 + t84)
      t60 = t25 * intHLs16s25s26s34s56x1114D6eps1()
      t214 = -t60 + t148
      t215 = t25 * intHs16s25s26s34s56x1114D6eps1()
      t216 = t25 * intHs16s25s26s34s56x1141D6eps1()
      t217 = t216 - t136
      t218 = intHs16s25s26s34s56x1411D6eps1() * t25
      t219 = t217 * t70
      t220 = (t215 - t143) * s34
      t221 = (t9 + t13 + t14) * s25
      t222 = t20 * propZ25
      t223 = intHL0s25s26s34s56x1130D4eps0()
      t224 = -s13 + s25 + s26 + s56
      t225 = s45 + s46
      t226 = t16 * t225
      t227 = t224 + t226
      t228 = s13 - s45 - s56
      t229 = -t16 * t228 + s25 + s46
      t230 = intHLs160000x0211D2eps0()
      t231 = -s16 + s25 + s34
      t232 = t16 * (s12 + s13 - s56)
      t233 = t18 * s46
      t234 = t10 * s45
      t235 = t231 + t233 - t232 + t234
      t236 = t10 * t225
      t232 = t231 - t232 + t236
      t237 = intHs160000x0121D2eps0()
      t238 = s25 - s26 - s16 + s34
      t239 = s13 - s56
      t240 = t18 * t239
      t241 = t238 - t69 - t240 + t236
      t242 = intHLs16s25s26s34s56x1211D4eps0()
      t243 = s25 - s26
      t244 = s16 - s34
      t245 = t16 * t244
      t246 = t42 * s46
      t247 = t41 * s45
      t248 = t247 + t243 - t240 + t246 - t245 - t62
      t249 = intHs16s25s26s34s56x1111D2eps0()
      t250 = t16 * (s12 - s45 - s46)
      t251 = s26 + s16 - s34 + s13 - s56 + t250
      t252 = intHL0s25s260s56x1012D2eps0()
      t253 = t48 - t109
      t254 = intHL0s25s260s56x1013D4eps0()
      t255 = intHL0s25s260s56x1021D2eps0()
      t256 = intHs160000x0111D0eps0()
      t257 = intHL0s25s26s34s56x1210D2eps1()
      t258 = intHL0s25s26s34s56x1310D4eps1()
      t259 = intHs160s26s34s56x1012D2eps1()
      t260 = s26 + s16 - s34 - s56
      t261 = t16 * (s12 + s13 - s45 - s46)
      t262 = -t260 - t261
      t263 = t222 * s25
      t264 = t19 + t263
      t265 = s25 + s16 - s34
      t266 = t265 + t83
      t267 = s13 + s16
      t268 = t267 + t69
      t269 = -s16 + s25
      t270 = t18 * s13
      t271 = t16 * s34
      t272 = -t269 - t271
      t273 = t18 * s34
      t274 = -t269 - t273
      t275 = intHLs160000x0211D2eps1()
      t276 = intHs160000x0121D2eps1()
      t277 = intHL0s25s260s56x1022D4eps0()
      t278 = s26 + s56
      t279 = t278 + t57
      t280 = s45 + s46 + s56
      t281 = t16 * t280
      t282 = s13 * t279
      t283 = t278 * t280
      t284 = (t281 + s26 + s25) * s25 - t282 + t283
      t285 = s12 + s13
      t286 = t16 * s13
      t287 = t265 + t79
      t288 = t285 * s25
      t289 = propW16 * t266
      t290 = t289 * t48 * t33
      t291 = t10 * t285
      t292 = s25 + s26
      t293 = s12 + s13 + s16
      t294 = s12 + s16 + s26 - s34
      t295 = t293 * s25
      t296 = (-s13 + s25 + s26) * s34
      t297 = s13 * s16
      t298 = t225 * t231
      t299 = (s12 + s16 + s25 + s26) * s26
      t300 = s56 ** 2.0_dp
      t301 = intHL0s25s260s56x1010D0eps0()
      t302 = intHL0s25s260s56x1011D2eps0()
      t303 = intHL0s25s260s56x1020D2eps0()
      t304 = intHL0s25s26s34s56x1110D2eps0()
      t305 = intHs160000x0111D0eps1()
      t306 = intHL0s25s26s34s56x1130D4eps1()
      t307 = intHs16s25s26s34s56x1113D4eps0()
      t308 = -t16 * (s12 + s26 + s13 - s45 - s46 - s56) - t244
      t309 = intHLs160000x0111D2eps0()
      t310 = intHLs16s25s26s34s56x1411D6eps0()
      t311 = intHs160s26s34s56x1011D2eps0()
      t312 = t16 * (s12 + s13 - s45) + s26 + s16 - s34 - s46 - s56
      t313 = intHLs16s25s26s34s56x1311D4eps0()
      t314 = t16 * s45
      t315 = s25 + s26 - s13 + s46 + s56 + t314
      t316 = intHL0s25s260s56x1012D2eps1()
      t317 = intHL0s25s260s56x1013D4eps1()
      t318 = t16 * t33
      t319 = t33 - t57
      t281 = t281 + s25
      t320 = s13 * (t318 - t21) - t319 * t281
      t321 = s25 + s56
      t322 = t16 * t278
      t323 = t281 * (t322 + t21) - t286 * (t16 * t321 + s26)
      t324 = s12 + s13 - s16
      t325 = s25 - s12
      t326 = -t286 + s25
      t327 = s16 - s34 - s13
      t328 = s12 - s26
      t329 = s13 - s45
      t330 = t285 * s13
      t331 = (-s34 + t293) * s25 + (s25 - s12 - s13) * s26 + (-s25 - s12
     # + s26 - t286 + s45) * s45 + (-t16 * t329 + s46 - t328) * s46 + (-
     #t328 - t109 + s56) * s56 + t330
      t332 = (t16 * t327 + s45 + t325 + t43) * s45
      t333 = t18 * t300
      t334 = s25 ** 2.0_dp
      t335 = s26 + s34
      t336 = s25 * t225
      t337 = t10 * t336
      t338 = t16 * (s34 * t280 + t334) + s13 * (-t21 + t245) - s16 * t28
     #1 + t335 * s25 + t337 + t21 * s56
      t339 = s13 * (t129 + s25)
      t340 = t48 * t281
      t341 = -t340 + t339
      t342 = t318 + s25
      t343 = t69 * s13
      t344 = t48 + t226
      t345 = s26 * s34
      t346 = s13 * t274 - s16 * t344 + s25 * t48 + t10 * s34 * t225 + t1
     #6 * (s34 + s45 + s46) * s25 + t345 + t273 * s56
      t347 = intHL0s25s26s34s56x1220D4eps0()
      t348 = t48 + t271
      t349 = s25 + s26 + s34 + s56
      t350 = intHL0s25s26s34s56x1120D2eps1()
      t351 = t98 * s34
      t352 = epinv * s34
      t353 = t352 - s25 - s26 - s56
      t354 = t16 * (s13 - s45 - s46 - s56)
      t355 = t10 * s56
      t356 = s46 ** 2.0_dp
      t357 = t10 * s56 * t225 + t16 * (s45 * t244 + (s25 + s16 - s34 + s
     #45) * s46 + s56 * t265) + t18 * (s26 * t280 + t300) - s12 * (s25 +
     # s26 - s13 + s45 - s46 + s56) + (-t16 * (s16 - s34 + s45 + s46) - 
     #s25 - t43 - t355 + s13) * s13 + (s26 + s16 - s34 + s45) * s25 + t3
     #56 + s45 ** 2.0_dp
      t358 = 7.0_dp
      t359 = t41 * t225
      t360 = t358 * s56
      t361 = s46 * s56
      t362 = s26 + s45 + s56
      t363 = s26 * s46
      t364 = t42 * s13
      t365 = t364 * s25
      t366 = t16 * (-s13 * t278 + s26 * s45 + s56 * t362) + t18 * (t361 
     #+ t334) + s25 * (t360 + t359 + s26) + t363 - t365
      t367 = intHL0s25s260s56x1021D2eps1()
      t368 = epinv * s56
      t369 = t25 * s26
      t370 = t25 * s34
      t371 = -t370 + s26 + s16 + s56
      t226 = t260 - t226
      t372 = t351 + t33
      t373 = intHs160000x0211D2eps1()
      t374 = -t354 + s25
      t375 = t18 * epinv
      t376 = -t10 + t375
      t377 = intHL0s25s260s56x1031D4eps1()
      t378 = s26 - s56
      t379 = s46 * t378
      t380 = intHLs16s25s26s34s56x1311D4eps1()
      t381 = -s13 + s16
      t382 = t364 - t57 - s16
      t383 = s26 + s16
      t384 = s26 + s45 + s46
      t385 = t18 * t381
      t386 = t10 * t384
      t387 = t42 * s34
      t388 = s12 + s16
      t389 = s13 * t388
      t390 = t57 * t381
      t391 = t10 * t389 - t333 - (-t382 + s26) * s26 + (-t364 + t21 + s2
     #6) * s34 + (-t57 - t385 + t387 - t386) * s56 + t225 * (s34 * t41 -
     # t10 * t383) - t390
      t392 = t18 * s16
      t393 = t10 * s25
      t394 = -s13 * t358 + s26 + t393
      t395 = t41 * s16
      t396 = t358 * s34
      t397 = s12 * s13
      t243 = t16 * (t397 + t334) - t333 - (t392 + s13) * s25 - (-t364 + 
     #t82 + s26) * s26 + t394 * s34 + (-s16 * t42 + s25 + t270 - t386 + 
     #t396) * s56 + t225 * (t10 * t243 + 8.0_dp * s34 - t395) + t364 * s
     #16
      t398 = t16 * s16
      t399 = s25 - s12 + s16 - s34
      t400 = t16 * (s26 - s13 + s45)
      t401 = s34 ** 2.0_dp
      t329 = -t16 * (-s16 * t285 + (s12 + s13 + s16 + s26) * s34 + t300)
     # + (t18 * t285 + s26 + t398) * s26 + (-t399 - t400 - s46) * s46 + 
     #(t10 * t329 + t18 * (s12 - s46) - t292) * s56 + t288 - t330 - t332
     # + s16 ** 2.0_dp + t401
      t402 = intHLs16s25s26s34s56x1321D6eps1()
      t403 = s25 + s16 - s13
      t404 = t16 * t325
      t405 = t10 * s34
      t406 = -t16 * (t225 * (t292 + t271) + t397) + (t69 + s16) * s25 - 
     #(-t403 - t69 - s26) * s26 - (-t270 + t404 - s16 + s34) * s34 - (s2
     #5 - s16 - s13 + t405 - t250 + s56) * s56 - t297
      t407 = t10 * s13
      t408 = t381 * s25
      t409 = (t79 + t21 - t407) * s34
      t384 = t16 * (-t225 * (t383 - t273) + t389) - (t82 - t270 + s26) *
     # s26 - (t16 * t384 + s56 + t403 - t405) * s56 - t408 + t409
      t403 = s25 - s13
      t410 = s26 - s13
      t411 = t16 * t410
      t412 = s16 * s25
      t413 = t403 * s34
      t399 = t16 * t300 + s13 * t324 + (t325 - t286) * s26 + (-s12 + s16
     # - s34 + t411 + s45) * s45 + (t265 + t400 + s46) * s46 + (-t108 * 
     #t18 + t399 + t79) * s56 + t412 - t413
      t414 = intHLs16s25s26s34s56x1231D6eps1()
      t415 = t16 * t292
      t416 = s25 + s45 + s46
      t328 = t16 * t328
      t417 = t10 * t416
      t381 = (-t381 - t69 - s26) * s26
      t418 = (-t268 + s25) * s25
      t419 = t16 * (t225 * (t415 + s34) + t397) + t333 + (t16 * t403 + s
     #26) * s34 + (t18 * (s34 - s13) - s16 - t328 + t417) * s56 + t297 +
     # t381 + t418
      t420 = t358 * s25
      t421 = (s25 + s16 + s13 + t69 + s26) * s26
      t335 = t10 * (t335 + t57) * t225 + t16 * (s16 * t108 + t334 + t397
     #) + t333 - (t392 + t291) * s25 + (t10 * t403 + s26) * s34 + (-t18 
     #* t267 + t236 - t328 + t387 + t420) * s56 - t421
      t422 = s16 + s26 + s56
      t423 = t422 - t273
      t424 = t422 - t271
      t425 = s13 * t423 + t226 * t424 + t69 * t424
      t426 = t16 * t388
      t413 = t16 * (t413 + t297) + (-s26 - t270 - t426 + s25) * s25 + (-
     #t245 + t393) * s45 - (s26 + t245 - t393) * s46 + (t21 + s46 - t245
     #) * s56
      t427 = intHs16s25s26s34s56x1112D2eps1()
      t428 = intHL0s25s26s34s56x1220D4eps1()
      t429 = t359 + s26 - s16 + s34 + t404 - t240
      t430 = intHLs16s25s26s34s56x1121D4eps0()
      t431 = s46 + s56
      t247 = t42 * t431 + s26 + t21 + t247 - t407
      t432 = s16 + s26 - s34
      t359 = s56 * t42 + t359 + t404 - t407 - t432
      t246 = -t228 * t41 + t21 + t246
      t404 = s26 - s34
      t293 = -t16 * (s13 * t293 + (s25 - s16 + s34 - s13) * s45 - t295 +
     # t296 - t299 - t356) + (t16 * (s16 + s45) + t404) * s46 + (t16 * (
     #s12 + s26 + s16 - s34 + s13) + t233) * s56
      t433 = t80 * s34
      t410 = -t410 - t57
      t434 = (t364 + t398 + t62) * s25
      t314 = (-t314 + s25) * s46
      t435 = s13 ** 2.0_dp
      t436 = t292 * s34
      t437 = (s12 + s13 + s16 + s25 + s26) * s26
      t438 = (-s25 + s12 - s45 + s46 + t286 - s56) * s56
      t439 = intHLs16s25s26s34s56x1321D6eps0()
      t440 = intHLs16s25s26s34s56x1131D4eps0()
      t285 = t16 * t285
      t441 = s12 - s45
      t442 = t16 * t441
      t443 = (t285 + s16) * s25
      t444 = intHLs16s25s26s34s56x1221D4eps0()
      t445 = -t43 - t57 + s13
      t446 = s25 + s26 + s34
      t447 = -s25 + s16 - s34 + s13 - s46 + t442
      t448 = t446 * s46
      t449 = t16 * (s45 * t446 - t397)
      t450 = (s12 - s13 + s16 + s25 + s26) * s26
      t411 = (-t411 - s25) * s34 + t16 * ((t16 * t383 + s12 - s13 + s25 
     #- s34 + s45 + s46 + s56) * s56 + t225 * (t404 + t398) - t397 + t41
     #2 + t450) - t407 * s16
      t412 = t16 * (s26 + s16 - s13 + s45 + s46)
      t271 = t225 * (t383 - t271)
      t401 = t16 * (t271 - t297 - t401 + t300) + (t267 + s26) * s25 + (t
     #43 + t398 + t291) * s34 + (-t273 + t412 + s25) * s56
      t451 = s25 - s34
      t452 = t326 * s26
      t453 = t225 * t432
      t267 = t16 * (t453 - t297 + t300) + s25 * t267 + (t383 + t285 - s3
     #4) * s34 + (t451 + t412) * s56 + t452
      t383 = t270 - t57 - s26
      t327 = (t16 * (s12 + s45 + s46) + t18 * t327 + s25 + t37) * s56
      t403 = (t18 * t403 + s16 + s26 + t69) * s26
      t238 = t225 * t238
      t412 = t16 * (-t295 + t238 + t397)
      t454 = t18 * (-t297 + t300)
      t455 = intHL0s25s260s56x1022D4eps1()
      t456 = intHs16s25s26s34s56x1211D2eps0()
      t457 = intHLs160000x0121D2eps0()
      t458 = t289 * s16
      t459 = t16 * (t458 - s25)
      t460 = t458 * t41
      t234 = t19 * (t29 * (-t18 * t431 - s26 - t234 + t270 - t57) + t460
     #) + t22 * (t18 * (s13 - s46 - s56) - s26 - t234 + t459)
      t459 = t10 * t228 - t233 + t459
      t233 = t29 * (-t10 * (s45 + s56) - t233 - t57 + t407) + t460
      t461 = intHs160s26s34s56x1013D4eps1()
      t462 = intHL0s25s26s34s56x1210D2eps0()
      t463 = intHL0s25s26s34s56x1310D4eps0()
      t464 = intHs160s26s34s56x1012D2eps0()
      t465 = intHs160s26s34s56x1021D2eps0()
      t466 = intHs160s26s34s56x1013D4eps0()
      t467 = intHLs160000x0122D4eps0()
      t468 = intHs160s26s34s56x1011D2eps1()
      t469 = intHL0s25s260s56x1011D2eps1()
      t470 = intHL0s25s26s34s56x1110D2eps1()
      t471 = intHL0s25s260s56x1010D0eps1()
      t472 = intHL0s25s260s56x1020D2eps1()
      t473 = intHL0s25s26s34s56x1120D2eps0()
      t285 = (-t285 - t392 + s25) * s25
      t450 = t16 * ((-s25 + s12 + s26 + t245) * s56 - t397 + t450)
      t474 = t10 * (t298 + t297)
      t475 = -t21 + t407
      t476 = s26 + s16 - s34 - s13 + s45 + s46
      t477 = -s34 * t383
      t478 = intHLs16s25s26s34s56x1221D4eps1()
      t479 = t16 * t476 + s25
      t453 = t453 - t389 + t300
      t480 = t16 * t326
      t441 = t441 * s25
      t481 = (t362 - s13) * s13
      t260 = t260 * s25
      t361 = -t16 * (s45 * s46 + t356 + t441 + t481) - t18 * (s13 * s25 
     #+ t361) - t260 - t363
      t363 = s25 * s26
      t260 = -t16 * (s46 * t362 + t356 + t441 + t481) - t260 - t21 * s13
      t362 = intHLs16s25s26s34s56x1141D6eps0()
      t441 = s26 ** 2.0_dp
      t481 = s12 * t48
      t482 = t416 * s16
      t483 = intHLs16s25s26s34s56x1131D4eps1()
      t484 = intHLs160000x0121D2eps1()
      t485 = intHLs16s25s26s34s56x1231D6eps0()
      t486 = t18 * t244
      t487 = intHs160000x0211D2eps0()
      t488 = intHL0s25s260s56x1031D4eps0()
      t489 = intHLs16s25s26s34s56x1211D4eps1()
      t490 = t10 * t239
      t491 = 8.0_dp * s45
      t492 = t358 * s46
      t493 = 8.0_dp * t225
      t494 = intHLs16s25s26s34s56x1121D4eps1()
      t495 = 0.1e1_dp / t48
      t496 = s34 * t39
      t497 = t25 * t316
      t498 = t155 * t157
      t499 = t120 * t376
      t500 = intHs16s25s26s34s56x1112D2eps0() * s34
      t226 = t88 * (t105 * (epinv * (-t16 * (-t298 - t397 + t441) - (t42
     #6 + s13) * s25 - (t16 * (s12 - s13) + s16 + t21) * s26 + (s25 - t8
     #3 - t486) * s56 + t477 + t270 * s16) + t285 + t409 - t450 + t474) 
     #+ t119 * (epinv * t359 + s16 + s26 - s34 - t21 - t360 + t364 - t49
     #3 + t69) + t121 * (epinv * t429 + s16 - s34 - t21 + t328 + t490 - 
     #t493) + t157 * (t4 * (s13 * t372 + t226 * t371 + t69 * t371) + t5 
     #* t425 + t498 * t329) - t159 * t243 + t499 * t374) * t27
      t360 = t88 * t24
      t324 = t360 * (t157 * (-t331 * (t93 + t94) - t71 * (t286 * s16 + e
     #pinv * t331 + s25 * t324 + (-t325 + t270) * s26 + t326 * s34 - (t1
     #6 * (s25 + s16 - s34 - s13 + s45) + s12 + t43 + s46) * s46 + (t10 
     #* t108 - t16 * t265 + s12 - t43) * s56 - t330 - t332 - t333)) + t5
     #00 * t374)
      t325 = 0.1e1_dp / ecossin ** 2.0_dp * propW34 * im
      t330 = epinv * t123
      t332 = s34 * t106
      t371 = epinv * t275
      t397 = t305 * t27
      t339 = t198 * (t104 * t335 + t229 * (-t332 + t330) - t47 * (epinv 
     #* (-t16 * (s13 * (s25 - s16 + s56) - s25 * s45 + (s25 + s16 - s34 
     #- s45) * s26 + (-s25 + s16 - s45) * s56 - t300 + t441 + t481 + t48
     #2) + s25 * s34 + s46 * t48) - t339 + t340) - t49 * (epinv * t260 +
     # s25 * t227) - t371 * t235)
      t115 = t325 * (t24 * (t227 * t456 + t88 * ((t115 * t232 + t241 * t
     #276) * epinv + t103 * (t16 * (s12 * t224 + s26 * t278 - s34 * t416
     # - t336 + t482) + t18 * (s56 * t244 + t363) + s13 * (-t79 + s25 - 
     #t486) - s25 * s56 + s26 * t244) - t116 * (epinv * t411 - t326 * s3
     #4 + s56 * t479 + t16 * t453 + t408 + t452) + t157 * (-t370 * t259 
     #* t262 + t154 * t312 - t311 * t262) + t96 * (t10 * (t271 - t389) +
     # t333 + t390 + epinv * t401 + (-t382 + s26) * s26 - t394 * s34 + (
     #-t396 + t57 + t385 + t386) * s56))) - t90 * t27 * (epinv * t227 - 
     #s13 - s16 - s26 + s34 + s56 - t250) + t339 + (t397 * t227 + t234 *
     # t484) * t88 * epinv + t198 * (-t370 * t257 + t301 - t302 - t303 -
     # t304) * t495 * t253)
      t224 = t31 * (t252 + t255) * t495
      t232 = t325 * (t88 * (t1 * (t19 * t233 + t22 * t459) + t124 * (t19
     # * (epinv * t233 - t460) + t263 * (t375 * t459 - t460)) + t234 * t
     #457 + t24 * (t232 * t3 + t237 * t241 + t249 * (s13 + s26 - s56)) +
     # t27 * (t131 * t411 + t227 * t256 + t243 * t91 + t359 * t7 + t429 
     #* t6) + t31 * (-t113 * t293 + t122 * t229 - t147 * (t16 * (s45 * t
     #410 + t356 - t435 - t436 + t437 + t438) - t314 + t434) - t2 * (s25
     # + s26 + s16 - s34 - s13 - s46 + s56 + t69) - t230 * t235 + t242 *
     # t248 + t28 * t246 + t92 * (-t16 * ((s25 + s26 - s16) * s45 - t295
     # + t297 - t299 + t300) + (-t292 + t398) * s46 + (-t16 * (s25 - s12
     # - s16 - s13 + s45) - s46) * s56 - t433) + t430 * t247 - t444 * (s
     #34 * t445 + (t447 - s56) * s56 + t421 + t443 - t448 - t449))) + t2
     #24 * t253)
      t233 = -t471 + t472
      t234 = t464 + t465
      t235 = t262 * t157
      t241 = t25 * t157
      t243 = (t462 + t473) * s34
      t250 = t107 * (-epinv * t293 + t285 + t409 - t450 + t474) + t128 *
     # (t10 * (t225 * (t57 + s26) + t300 + t334) - t365 + epinv * (-t16 
     #* t410 * s45 - t16 * (-t436 + t438 + t437 - t435 + t356) + t314 - 
     #t434) + (-t286 + t21) * s26 + (9.0_dp * s25 + t10 * (s26 - s13 + s
     #45 + s46)) * s56) + t478 * (t10 * t446 * t225 + t333 + t343 + epin
     #v * (-s34 * t445 + (-t447 + s56) * s56 - t421 - t443 + t448 + t449
     #) + (t21 + s26 - t407) * s34 + (-t270 - t328 + t387 - s16 + t417) 
     #* s56 + t297 + t381 + t418)
      t271 = t325 * t88
      t66 = t271 * (t24 * ((t117 * t374 * t427 + t235 * t234) * s34 + t1
     #02 * t401 + t13 * (s34 * t383 + t327 + t403 - t412 + t454) + t157 
     #* (-t262 * t468 + t312 * t66) * t17 + t241 * t95 * t331) + t27 * (
     #t134 * (t10 * (t238 + t297 - t300) + t343 + epinv * (t477 - t327 -
     # t403 + t412 - t454) - (t69 + t392 + s13 - s25) * s25 + s26 * t475
     # - t475 * s34 + (-t10 * t476 + s25) * s56) - t241 * (t73 * t329 * 
     #t157 + t114 * t425)) + t31 * t250 + t31 * (t17 * (-t469 - t470) + 
     #t233 * t25 + t243) * t495 * t253)
      t238 = t317 * t495
      t250 = epinv * t373
      t262 = t25 * t76
      t285 = t262 - t138
      t293 = t25 * t81
      t312 = t25 * t306
      t314 = t25 * t461
      t326 = t132 * t24
      t82 = epinv2 * (t19 * (t16 * (t225 * (t244 + t43 + t57) - t297) + 
     #t18 * (t290 + t300) - (-propZ25 * (s56 * t294 + t295 - t296 - t297
     # - t298 + t299) + s16 + t291) * s25 - (s13 * t41 + s26 + t62 + t82
     #) * s26 + (t292 + t286) * s34 - (t108 * t41 - t21 + t244 + t62 - t
     #79) * s56) + t22 * (-s13 * t244 - (t286 + s12) * s26 + (s25 - s12 
     #+ s26 - t109 + s56) * s56 + t225 * t287 - t288 + t290)) * t495 * t
     #157 + t326 * (t16 * t479 * s56 + epinv * t267 + s26 * t480 - t480 
     #* s34 + t10 * t453 + t390)
      t108 = t198 * (-t139 * t323 - t150 * t413 - t156 * (-t18 * s25 * t
     #239 + epinv * t361 + t16 * t334 + t337 + t363 - t379) + (-t293 + t
     #140) * (t16 * t272 * t225 - (t69 + s13 + s16) * s16 + s25 * t268 +
     # s26 * t269 - (t16 * (s25 - s12 - s16) - s26 - t270 + s34) * s34 +
     # s56 * t274) + t439 * t406 - t485 * t419 - t483 * (-epinv * (t16 *
     # (-s45 * t48 - t345 + t481) + (s16 - s34 - s46) * s25 + (s25 + s16
     # - s46) * s26 + (-s25 + s16 - s46) * s56 + t282 - t300 + t441) - t
     #16 * t341) + (-t370 * t253 * t258 + t277 * t284) * t495 - t312 * t
     #227)
      t82 = t325 * (t27 * (-t14 * t374 - t15 * (epinv * t374 - s16 + s25
     # - s26 + s34 + t236 - t240 - t69) + t88 * ((t314 * t235 - t307 * t
     #308) * s34 - t137 * t320 - t141 * t338 - t151 * t391 - t152 * t384
     # + t285 * ((t10 * t33 + t69) * s13 - t281 * (t318 - s25)) + t209 *
     # t251)) - t88 * t82 + t108)
      t108 = t325 * (t88 * (t31 * (t223 * t227 + t440 * (-s46 * t292 - (
     #s25 - s16 - s13 + s46 - t442 + s56) * s56 + t421 - t433 + t443 - t
     #415 * s45)) - t27 * t87 * t251 - t9 * t24 * t267) + t254 * t31 * t
     #253 * t495 + t197 * (t157 * (-t27 * (t125 * t399 + t153 * t331) + 
     #t24 * t126 * t357) - t31 * t284 * t455 * t495) + t88 * (t31 * t253
     # * t463 * t495 - t235 * t27 * t466) * s34)
      t82 = t108 + t82 + t325 * (t25 * (t31 * (t88 * (t323 * t78 + t366 
     #* t56 - t374 * t377 - t402 * t406 + t413 * t59 + t414 * t419) - t2
     #38 * t253) + (t320 * t74 + t338 * t85 + t346 * t52 + t384 * t64 + 
     #t391 * t61) * t88 * t27) + t88 * (-t27 * (t144 * t346 + t207 * (ep
     #inv * t308 + s26 + t236 - t240 + t57)) + t31 * (-(t313 * t315 + t3
     #80 * (epinv * t315 - s16 + s25 - s26 + s34 + t236 - t240 - t69)) *
     # s34 - t145 * t361 - t149 * t366 + t374 * t488) - t374 * (t250 + t
     #487) * t24) + t88 * (t27 * (t32 * t331 + t399 * t67) - t24 * t68 *
     # t357) * t157)
      t108 = t25 * intHLs160000x0122D4eps1()
      t225 = intHLs160000x0111D2eps1() * t17
      t235 = (-t108 + t467) * s16 + t309 + t225
      t236 = t25 * intHLs16s25s26s34s56x1411D6eps1()
      t239 = t25 * intHLs16s25s26s34s56x1141D6eps1()
      t227 = t325 * (t88 * (t27 * ((s13 * t342 - t281 * t33 + t343) * (t
     #216 - t136) + (t215 - t143) * t227 * s34) + t31 * (t341 * (-t239 +
     # t362) + (t60 - t148) * (t16 * s25 * (t280 - s13) + t334 - t379) +
     # (-t236 + t310) * t251 * s34)) + t27 * (t218 - t127) * t374 + t289
     # * t88 * t235 * t264)
      t240 = t25 * intHLs160000x0113D4eps1() - intHLs160000x0113D4eps0()
      t251 = 0.1e1_dp / 0.3e1_dp
      t267 = -s16 + s25 - s26 + s34 - s56
      t268 = epinv * s25
      t270 = epinv * t33
      t274 = s25 + s34
      t280 = t117 * s34
      t281 = t207 * t88
      t76 = (-t31 * t483 + t88 * (t124 * t31 + t326 * t267)) * epinv + t
     #25 * (t24 * (-t153 * t157 + t88 * (-(t157 * t461 + t52) * s34 - t2
     #31 * t85 - t267 * t74)) + t198 * ((-t258 * t495 + t81) * s34 + t23
     #1 * t59 + t279 * t56 + t48 * t78) - t76 * t27 * t33 * t88) + t24 *
     # (t88 * ((t144 + t307) * s34 + t141 * t231 + t267 * (t9 + t137)) +
     # t14) + t31 * (t25 * (-t377 * t88 - t238) - t88 * (t139 * t48 + t1
     #49 * t279 + t150 * t231 + t156 * (t352 + s25)) - t440 - (t380 + t1
     #40 + t145) * t88 * s34) + t88 * (t19 * (t12 * (t487 + t68) + t29 *
     # t488) + t27 * (t138 * t33 - t183 + t250)) + t175 * t19 * t88 + (s
     #34 * t463 * t88 + t254) * t495 * t31 + t24 * (t15 + t281) * t117 +
     # t24 * (s34 * t466 * t88 + t32) * t157
      t238 = t127 - t218
      t250 = t238 * t24
      t282 = t214 * t31
      t284 = t25 * t494
      t286 = t100 * t376
      t288 = epinv * t484
      t289 = t276 * epinv
      t290 = t311 * t157
      t259 = -s34 * t259 - t73
      t291 = t25 * t259 + t352 * t4
      t295 = -t25 * (s34 * t114 + t125 * t267) + t158
      t296 = t279 * t455
      t252 = t31 * (-t495 * (t252 + t255) + t88 * (-t104 * t349 - t107 *
     # (t280 + s25 - s16) - t128 * (epinv * t274 + s25 + s26 + s56) + t3
     #3 * t444 + t478 * (t270 - s34) + t495 * (t25 * ((-t428 + t257) * s
     #34 - t296) + t302 + t303 + t304 - t352 * t350) - t430)) + t88 * (t
     #91 * t27 * t424 - t134 * t24 * (epinv * t319 + s16 + s26 - s34 + s
     #56)) + t31 * (t495 * (-t301 * t88 + t497) + t88 * (-t147 * t274 + 
     #t269 * t92 + t284 - t286 + t288) + t88 * (t99 - t113) * s34) + t88
     # * (t19 * (t12 * t237 - (t249 + t7) * t12 + t29 * t457) + t27 * (t
     #131 * t269 + t289 + t500)) + t360 * (-t102 * t272 + t103 * t231 - 
     #t119 * t25 - t13 * t319 + t96 * (-epinv * t272 + s34) + t290 + t49
     #9) + t88 * (t24 * t295 + t27 * t291) * t157
      t255 = epinv * t116
      t272 = t25 * t95
      t274 = t105 * t117
      t297 = t71 * (t268 + s26 + s16 - s34 + s56)
      t298 = (t464 + t465) * s34
      t299 = s34 * t5
      t300 = t17 * t468
      t308 = t67 * t267
      t277 = t277 * t279
      t315 = t367 * (t72 + s26 + s56)
      t319 = t28 * t31
      t320 = t26 * t8
      t323 = 9.0_dp * t222
      t326 = s25 + s26 - s34 + s56
      t327 = t398 * propW16
      t328 = 1.0_dp + t327
      t329 = t395 * propW16
      t331 = t11 + t329 + t16
      t333 = t41 * propW16
      t314 = t466 - t314
      t293 = -t140 + t293
      t334 = t25 * t59
      t336 = t334 - t150
      t337 = t25 * t56
      t338 = t337 - t149
      t339 = t25 * t78
      t340 = (t32 - t213) * t157
      t341 = t102 + t103
      t343 = t7 + t8
      t345 = t1 * t331
      t346 = s34 * t427
      t356 = (t276 + t373) * epinv
      t318 = epinv2 * (t19 * (t333 * t33 - t10 + t11) + t22 * (t318 * pr
     #opW16 - 1.0_dp)) * t157
      t331 = t124 * ((t19 * t331 + t22 * t328) * epinv - t329 * t264)
      t234 = t157 * (t27 * (t88 * (s34 * t234 + t155 + t297) + t93 + t94
     #) + t360 * (t308 + t154 + t299 + t300)) + t157 * (t27 * (t291 * t8
     #8 - t272) + t360 * t295) + t88 * (t274 * t24 * t48 + t255 * t27 * 
     #t326 + t212 + t331) + t88 * (t19 * (-t12 * t343 + t12 * (t237 + t4
     #87) + t345) + t24 * (t25 * (-t119 - t120) + t48 * t341 + t96 * (ep
     #inv * t48 + s34) + t290) + t27 * (t131 * t326 + t346 + t356) + t31
     # * (s34 * t99 - t107 * (epinv * t422 + s25 + s26 + s56) - t113 * t
     #422 - t128 * (t72 + t352 + s26 + s56) + t147 * t451) + t318) + t19
     #8 * (((t39 + t49) * epinv + t106) * s34 + t162 * t349 - t122 - t16
     #7 + t168 - t330)
      t290 = t1 * t328
      t291 = propW16 * s16
      t295 = t291 * t264
      t326 = t295 * t88 * t240
      t108 = t16 * (t88 * (t31 * (t293 * s34 + t145 * t451 - t156 * (t72
     # + t352) + t231 * t336 + t279 * t338 + t48 * (t339 - t139) - t28) 
     #+ s34 * (t157 * t314 + t142) * t24 - t183 * t27) + t340 * t24 + t2
     #6 * (t68 * t88 - t13 - t14 - t9) + (t132 + t134 + t15) * t24 * epi
     #nv) + t18 * t222 * (t290 + t2 + t7 + t8 - t237 - t487) - t41 * (((
     #t467 - t457 - t108 - t288) * s16 + t225 + t309) * t88 * t264 * pro
     #pW16 + t222 * (-t221 + t68) + t282) + 12.0_dp * t326 + t234
      t225 = t16 + epinv
      t234 = 12.0_dp * t291
      t279 = t11 + t234 + t16
      t291 = t291 * t10
      t309 = 1.0_dp + t291
      t328 = t10 * t388
      t352 = t42 * s26
      t357 = s12 + s26 + s16 - s34 + s56
      t359 = t16 * t357
      t357 = t10 * t357
      t361 = -s13 + s25 + s26 + s45 + s56
      t363 = -t16 * (s16 - s34 + s13 - s45 - s56) + s25 - t62
      t365 = t400 + t21 + t355
      t64 = t25 * t64
      t61 = t25 * t61
      t52 = t25 * t52
      t366 = t25 * t85
      t194 = (-t194 - t195) * t25
      t195 = t15 * t225
      t370 = t27 * (t88 * (-t298 - t155 - t297) - t93 - t94) - t360 * (t
     #154 + t300 + t308)
      t95 = t157 * t370 + t157 * (t25 * (-t259 * t88 + t95) * t27 + t360
     # * (t182 * t267 - t158 - t311)) + t19 * (t12 * (-t88 * (t237 + t48
     #7) + t13 + t14 + t9) - t345 * t88) + t88 * (t24 * (-t131 * t388 - 
     #t132 * (t97 - t357) - t134 * (t268 - t359) + (t366 - t141) * (-t16
     # * t294 + s25) + (t25 * t74 - t137) * (-t359 + s25)) + t27 * (t285
     # * (t357 + s25) + (t52 - t144) * (t16 * (s12 + s26 + s16) - t273) 
     #+ (t64 - t152) * (t321 + t43 + t426 - t405) + (t61 - t151) * (t404
     # * t42 + t110 + t328 + t57) + t209 - t346) + t31 * ((epinv * t380 
     #+ t313) * s34 + t156 * (t98 * t361 + s25) + t48 * (epinv * t483 + 
     #t440)) - t318 - t26 * t87) - t88 * (t276 + t373) * t27 * epinv - t
     #195 * t24 - t331 * t88 + t198 * (t194 - t168 - t169 + t122 + t330)
     # + t198 * (t201 + t84)
      t241 = t24 * (t341 * t80 + t96 * (t16 * (epinv * s26 - s25) - t110
     # + t268 + t396 - t328 - t352) + t499) + t31 * (t38 * (t161 - t104)
     # - t242 + t284 - t286 - t430) + t212 - t241 * t114 * t27 * t423
      t142 = t222 * (-t290 + t487 - t221 + t237 + t87) + t88 * t217 * (t
     #359 + s25) * t27 + t360 * (-t142 - t143 + t215) * s34 - t250 + t28
     #2
      t184 = t88 * (-t184 * t27 + t31 * (-t145 * t361 + t292 * t338 + t3
     #36 * t53)) + (t314 * t88 * s34 - t213 + t32) * t157 * t24
      t217 = t88 * (-t26 * t174 - t319) + t222 * (t309 * t457 + t2 + t23
     #0 - t3)
      t95 = t10 * t184 - t16 * t95 + t18 * t217 + 24.0_dp * t326 - t41 *
     # t142 + t88 * t241 + t27 * (t88 * (t157 * (t372 * t4 + t423 * t5) 
     #+ t160 * (t358 * t404 + t388 * t41 + s25 + t110) + t176) + t90) + 
     #t6 * t89 + t88 * (t19 * (t12 * t3 + t230 * t29 + t279 * t457) + t2
     #4 * (t105 * (s25 * t117 - t10 * (s12 + s16 - s34) - t16 * t369) + 
     #t116 * (t16 * (t117 * t388 + s26 - s34 + s56) + s25) + t119 * t376
     # + t180) + t288 * (t19 * t279 + t22 * t309)) + t198 * (-t107 * (ep
     #inv * t363 + s25 + t79) - t113 * t363 - t128 * (epinv * t365 + s25
     # + t79) - t147 * t365 + t163 * (t16 * (s26 - s13 + s45 + s56) + t2
     #1) + t25 * t489 + t349 * (-epinv * t478 - t444) + t371) - 12.0_dp 
     #* t88 * t264 * propW16 * t235 + 9.0_dp * t222 * t174
      t73 = t73 * t25
      t142 = t117 * s12
      t184 = t16 * (t142 + s25 + s26)
      t217 = s26 - s16 + s34 - s56 + t69 + t21
      t241 = s25 - s26 + s56
      t259 = t43 + t69 + s56 + t393
      t267 = s16 + s56
      t273 = epinv * t51
      t129 = s25 - s13 + s45 + s46 + t129
      t276 = t292 - t109 + t110
      t52 = t52 - t144
      t144 = t25 * t414
      t61 = (t61 - t151) * t63
      t63 = t380 * s34
      t15 = t19 * (t12 * t9 + t88 * (-t12 * t7 - t211)) + t24 * (t88 * (
     #-t210 + t209 - t176 - t208) - t340 + t191 * t88 * t157 + t192 * t8
     #8 * t157) + t27 * (t25 * (t74 + t85) - t281 + t61 * t88 + t88 * (t
     #64 - t152) * t65) + t27 * (t231 * t88 * t52 - t137 - t141 - t15) +
     # t31 * (t88 * (t145 * t276 + t147 * t129 + t156 * (t276 * epinv - 
     #s25 + s26 - s56) + t231 * t293 + t241 * (t25 * (t56 + t59) - t149 
     #- t150) + t483 * (t97 - t37) + (t144 - t485) * (t267 + t69 + t21 +
     # t352) + (t25 * t402 - t439) * (s16 + s34 + s56 + t43 + t36) - t31
     #2) + t440) - t30 * t242 * t88 + t198 * (t223 - t63)
      t30 = t51 * t99
      t37 = t198 * (-t107 * (epinv * (t404 + t109 - t110) + t184) - t113
     # * (t404 + t261 - t110) + t128 * (t16 * (epinv * t129 - s26) - t21
     #) + t478 * (-t69 * t25 - s16 - s56 - t21 + t268 - t352) + t122 - t
     #169 + t230 + t28 + t347 + t350 + t371 - t30)
      t23 = t157 * (-t360 * t154 - t172 * t27) + t24 * (-t117 * t90 + t8
     #8 * (-t102 * t80 - t164 - t165)) + t27 * (-t91 * t217 * t88 + t170
     #) - t171 + t24 * (-t289 * t88 - t456) - t3 * t89 + (-epinv * t305 
     #- t101 * t119 - t103 * t266 - t105 * (epinv * t287 + t184) - t121 
     #* t376 - t13 * (t33 + t69) - t134 * (t270 + t118 + t21) + t159 * t
     #217) * t88 * t27 + t198 * (-t101 * t489 - t104 * t259 - t17 * t494
     # - t241 * t92 + t25 * (t259 * t34 - t100 - t428) + t444 * t70 + t4
     #7 * (s26 * t225 - t368 + t72) - t49 * (t273 - s25)) - (-t12 * t237
     # + t12 * t256 + t29 * t430) * t88 * t19 - t23 * t88 * (t249 + t4) 
     #+ t37 - t181 * t198 - t88 * (-t5 * t12 + t175) * t19 + t360 * (t17
     #8 - t179)
      t15 = -t10 * t88 * (t285 * t35 * t27 + t31 * t292 * (t339 - t139))
     # - t16 * t15 + t18 * (t196 * t26 + t222 * (t430 - t256 + t237 + t2
     #49 + t1 + t4 + t5)) - t41 * (t88 * (t27 * (t219 + t220) + t31 * ((
     #t236 - t310) * s34 + t214 * t378 + t80 * (t239 - t362))) + t222 * 
     #(-s25 * t9 - t2 - t242 + t7)) - t23 - t323 * t6
      t23 = -t392 + s25 + t405 - t83
      t29 = 1.0_dp + t98
      t34 = -s13 + s25 + s45 + s46 + s56
      t35 = t16 * (s25 - s12 - s16 + s34) + s56 - t43
      t37 = t48 + t398
      t43 = t10 * t404 + t16 * t267 + s25
      t64 = t16 * t48
      t65 = -s26 + t393 + t110
      t44 = t79 + t44 + t355
      t69 = s25 * t41 + t10 * (s34 + s56) - t392 - t69
      t74 = t19 * (t11 - t329 + t16) + t22 * (1.0_dp - t327)
      t83 = -t16 * t432 + s25
      t85 = -s12 - s13 + s25 + s45 + s46 + s56
      t89 = s12 - s13 + s25 + s45 + s46 + s56
      t97 = 1.0_dp - t291
      t109 = t11 - t234 + t16
      t118 = t157 * t24
      t9 = t271 * (t24 * (t13 * t33 + t132 * (epinv * t342 - t57 - t62) 
     #+ t342 * t9 - t209 + t356 + t500) + t27 * (t134 * (-t270 + t36) + 
     #t204 + t205) + t31 * (t101 * (t489 + t494) + t139 * t44 + t149 * t
     #65 + t150 * t35 + t156 * (epinv * t146 - s26 + s56 + t21) + t495 *
     # (-t302 - t303 - t304) - t371) + epinv2 * (t19 * (t333 * t48 + t11
     # + t16) + t22 * (t64 * propW16 + 1.0_dp)) * t495 - t288 * t74)
      t11 = t325 * (t88 * (t24 * (t87 + t237 + t487) + t256 * t27 + t31 
     #* (t2 + t28 - t230) - t457 * t74) + t224 + t195 * t27)
      t9 = t11 + t9 + t325 * (t25 * (t88 * (t27 * (-t190 - t186) + t118 
     #* (-t125 * t34 - t126 * t89 - t153 * t85)) + (t88 * (-s34 * t257 +
     # t349 * t428 + t296) - t316) * t495 * t31) + t88 * (t27 * (t207 * 
     #t29 + t202 + t203 + t206 + t75) + t31 * (t140 * t23 + t495 * (-t27
     #7 + t301) + t199 + t200) + t131 * t24 * (s12 - s16 - s26 + s34 - s
     #56)) + t88 * (t346 * t24 + t397) * epinv) + t325 * (t24 * (t88 * t
     #157 * (t32 * t85 + t34 * t67 + t68 * t89) + t456) + t27 * (t197 * 
     #(-t189 - t185 - t187 - t188) + t14) + t198 * (t25 * (-t23 * t81 + 
     #t233 * t495 - t35 * t59 - t44 * t78 - t56 * t65) + t495 * (-t347 *
     # t349 + t350 * t353 + t243 - t315) - t495 * (t469 + t470) * t17))
      t11 = t157 ** 2.0_dp
      t13 = t325 * (t24 * (t29 * t90 + t88 * (t102 * t43 + t116 * (-t16 
     #* (t270 - t142) + s25) + t96 * (epinv * t43 - s26 - s34 + s56 - t5
     #7 - t62) + t177 - t178 + t4 + t5)) + t88 * (-t1 * (t109 * t19 + t2
     #2 * t97) + t124 * ((-t109 * t19 - t22 * t97) * epinv - 12.0_dp * t
     #295)) + t198 * (t104 * t69 + t128 * (t16 * (-epinv * t130 + s26) +
     # t355 + t420) + t37 * t92 + t47 * (epinv * t37 - t64) + t49 * (t27
     #3 - t57) - t122 + t167 + t30 + t332))
      t14 = t349 * t478
      t23 = t325 * (t31 * (t495 * (-t25 * t317 + t254 + (-t25 * t258 + t
     #463) * t88 * s34) + t88 * (t25 * (t348 * t402 - t306 - t377) - t34
     #8 * t439 + t223 - t242 - t430 + t488 - t14) + t88 * (t144 - t485) 
     #* (s34 + t64)) - t88 * (s34 * t307 + t6) * t27)
      t29 = t48 * t483 + t63
      t9 = t10 * t325 * (t27 * (t88 * (-t215 + t143) * s34 + t127 - t218
     # + t88 * (-t216 + t136) * t70) + t198 * (-t60 + t148) * t241) + t1
     #6 * t271 * (t24 * (t176 + t3) + t27 * t343) + 0.4e1_dp / 0.3e1_dp 
     #* t9 - 0.8e1_dp / 0.3e1_dp * t23 + 16.0_dp * t295 * t271 * t240 + 
     #0.16e2_dp / 0.3e1_dp * t325 * t198 * t29 + 0.2e1_dp / 0.3e1_dp * t
     #13 + 0.2e1_dp / 0.3e1_dp * t271 * (t27 * (t103 * t83 + t105 * (epi
     #nv * t83 + t21 + t62 + t79) + t121 * (epinv * t10 - t42) + t376 * 
     #(t119 + t120) + t46 * (-t159 + t91) + (-t272 + t93 + t94) * t85 * 
     #t157 + (s12 + s13 + s16 + s26 - s34 - s45 - s46) * (-t73 + t155) *
     # t11) + t31 * (t107 * (-t10 * t244 + t112 + t21 - t79) + t111 - t3
     #30 - t161 * t69 + t496 * t17) + t118 * (t71 * (t25 * t34 + t142) +
     # t158)) - 8.0_dp * t271 * (t31 * ((t236 - t310) * s34 + t48 * (t23
     #9 - t362)) + t264 * t235 * propW16)
      t13 = t354 + s25
      t23 = epinv * t80
      t34 = t354 + s25 + t62
      t35 = -t87 + t6
      t29 = t24 * (t221 + t210 - t209 - t206 - t75) + t24 * (t132 * t133
     # + t134 * t135 + t25 * (t188 + t187) + t52 * t54 + t77 * (t262 - t
     #138) + t86 * (t366 - t141) + t61 + t193 * s25) + t31 * (epinv * t2
     #9 + s34 * t313 - t156 * (t98 * t50 + s25) - t194 - t201 - t84) + t
     #207 * t27 + t26 * t35 + t440 * t31 * t48
      t11 = struc60MP * t27 * (t33 * (-t16 * (t68 + t32 + t67 - t183 - t
     #213 - t182) - t25 * t71 - t272 + t93 + t94) - t155 + t73) * t11
      t11 = struc9MP * (t10 * t31 * (-t145 * t50 + t292 * (t337 - t149) 
     #+ t53 * (t334 - t150)) + t16 * t29 - t18 * t263 * t174 + t41 * (t2
     #4 * (-s25 * t238 + t219 + t220) - t263 * t35 - t282 * s25) + t24 *
     # (-t105 * (-t23 + t79 + t21 + t62) - t46 * t91 - t96 * (-t23 + s26
     # + s34 - s56 + t57 + t62) + t166 + t179) + t31 * ((-t275 - t484) *
     # epinv - t100 * t17 - t104 * t38 - t107 * (-epinv * t34 + s25 + t7
     #9) + t113 * t34 - t128 * (-epinv * t13 + s25 + t79) + t13 * t147 +
     # t25 * (-t489 - t494) + t349 * t444 + t30) + t24 * (s25 * t90 - t1
     #01 * t121 - t173 * t17 + t25 * (t40 * t46 - t114) + t176 + t3 + t4
     # + t5) + t31 * ((t49 * t51 + t14) * epinv + t2 - t230 + t242 - t28
     # + t430 - t457 + t161 * t38) + t320 + t26 * t7) - t11
      t13 = t198 * (-t16 * (epinv * t156 + t145) + (-t107 - t128 + t49) 
     #* epinv - t113 - t147 + t99) * struc61MP
      t13 = t325 * ((t16 * t76 + t18 * (t88 * (-t320 - t319) + t222 * (t
     #457 - t237 + t249 + t2 + t7)) + t41 * (t250 - t282 + t222 * (-t487
     # + t488 - t68 + t1)) + t157 * (t27 * (t88 * (t298 + t297) + t94) +
     # t360 * (t299 + t300 + t308) - epinv2 * (t19 * ((-propZ25 * t231 +
     # t10) * s25 + t278 * t41 + t245) + t22 * (t265 + t322)) * t495 * t
     #88) + t27 * (t157 * (t93 - t272) + t88 * (-t159 * t424 + t255 * t2
     #69 + t280 * t427 + t498)) + t88 * (t24 * (t274 * t231 + t177) + t2
     #12) + t252 + t198 * ((s34 * t49 + t269 * t47 - t123) * epinv + (t1
     #17 * t39 + t106) * s34 + t495 * ((t347 - t462 - t473) * s34 + t17 
     #* (t469 + t470) + t25 * (t471 - t472)) - t122 + t161 * t349) + t19
     #8 * t495 * (t315 + t277) + t323 * t8) * struc2MP + (struc3MP + str
     #uc8MP) * t108 + t13)
      t9 = (t10 * t271 * t458 * t264 * t240 - t16 * t227 + t251 * (t66 +
     # t232 + t115 + t325 * (t31 * (t88 * (t100 * (epinv * t246 + 8.0_dp
     # * t228 - t393 - t492) - t99 * t260 + t489 * (epinv * t248 - s25 +
     # t45 + t486 + t490 - t491 - t492 + t79) + t494 * (epinv * t247 - t
     #358 * t431 - s26 + t364 - t393 - t491) + ((-t25 * t428 + t347) * (
     #s13 * t348 - t344 * t349) + t350 * (s13 * (-t351 + t48) + t344 * t
     #353) + t367 * (-t16 * (s13 * (-t72 - s26 - s56) + t283) + s25 * (-
     #t72 - t16 * (s45 * t25 + s46 * t25) - t110 + t368 - t369))) * t495
     # - t161 * t335 - t496 * (epinv * t229 + s25 - t354)) - t497 * t253
     # * t495) + t226 + t324)) - 0.2e1_dp / 0.3e1_dp * t82 + t271 * t8 *
     # t27 * t374) * struc1MP + t9 * struc7MP
      result = 0.4e1_dp / 0.3e1_dp * t13 + 0.2e1_dp / 0.3e1_dp * t325 * 
     #(t88 * t11 + (-t16 * t58 + t18 * t222 * (t1 + t3 + t4 + t5 - t7 - 
     #t8) - t41 * (t27 * (t219 * t88 + t220 * t88 - t127 + t218) + t198 
     #* t214 * t55 + (-t221 - t2 + t6) * propZ25 * t20) + t27 * (t88 * t
     #173 * t17 + t157 * t172 + t180 * t88) + t88 * (t19 * (t12 * t174 -
     # (t3 + t4 + t5) * t12 + t175) + t24 * (-t178 + t179 + t177 + t176)
     #) + t24 * (t88 * (t105 * ((-t18 + epinv) * s25 - t62 - t79 * t25) 
     #+ t164 + t165 + t166) + t90) + t27 * (t88 * t160 * t46 - t170) + t
     #171 + t88 * (t162 * t38 + t163 * t51 - t122 - t167 + t168 + t169) 
     #* t31 + t88 * (t107 * (t112 - s25 - t79) - t128 * (t16 * (epinv * 
     #t130 + s26) + s25) + t111 + t181) * t31) * struc10MP + t95 * struc
     #5MP + t15 * struc6MP) + t9

           ampNonresonantLightFullImC4MP = result
       end function ampNonresonantLightFullImC4MP

       function ampNonresonantLightFullImC4PM()
           implicit none
           complex(dp) :: ampNonresonantLightFullImC4PM
           complex(dp) ::  t1,t10,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109
           complex(dp) ::  t11,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t12
           complex(dp) ::  t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t13,t130
           complex(dp) ::  t131,t132,t133,t134,t135,t136,t137,t138,t139,t14,t140,t141
           complex(dp) ::  t142,t143,t144,t145,t146,t147,t148,t149,t15,t150,t151,t152
           complex(dp) ::  t153,t154,t155,t156,t157,t158,t159,t16,t160,t161,t162,t163
           complex(dp) ::  t164,t165,t166,t167,t168,t169,t17,t170,t171,t172,t173,t174
           complex(dp) ::  t175,t176,t177,t178,t179,t18,t180,t181,t182,t183,t184,t185
           complex(dp) ::  t186,t187,t188,t189,t19,t190,t191,t192,t193,t194,t195,t196
           complex(dp) ::  t197,t198,t199,t2,t20,t200,t201,t202,t203,t204,t205,t206
           complex(dp) ::  t207,t208,t209,t21,t210,t211,t212,t213,t214,t215,t216,t217
           complex(dp) ::  t218,t219,t22,t220,t221,t222,t223,t224,t225,t226,t227,t228
           complex(dp) ::  t229,t23,t230,t231,t232,t233,t234,t235,t236,t237,t238,t239
           complex(dp) ::  t24,t240,t241,t242,t243,t244,t245,t246,t247,t248,t249,t25
           complex(dp) ::  t250,t251,t252,t253,t254,t255,t256,t257,t258,t259,t26,t260
           complex(dp) ::  t261,t262,t263,t264,t265,t266,t267,t268,t269,t27,t270,t271
           complex(dp) ::  t272,t273,t274,t275,t276,t277,t278,t279,t28,t280,t281,t282
           complex(dp) ::  t283,t284,t285,t286,t287,t288,t289,t29,t290,t291,t292,t293
           complex(dp) ::  t294,t295,t296,t297,t298,t299,t3,t30,t300,t301,t302,t303
           complex(dp) ::  t304,t305,t306,t307,t308,t309,t31,t310,t311,t312,t313,t314
           complex(dp) ::  t315,t316,t317,t318,t319,t32,t320,t321,t322,t323,t324,t325
           complex(dp) ::  t326,t327,t328,t329,t33,t330,t331,t332,t333,t334,t335,t336
           complex(dp) ::  t337,t338,t339,t34,t340,t341,t342,t343,t344,t345,t346,t347
           complex(dp) ::  t348,t349,t35,t350,t351,t352,t353,t354,t355,t356,t357,t358
           complex(dp) ::  t359,t36,t360,t361,t362,t363,t364,t365,t366,t367,t368,t369
           complex(dp) ::  t37,t370,t371,t372,t373,t374,t375,t376,t377,t378,t379,t38
           complex(dp) ::  t380,t381,t382,t383,t384,t385,t386,t387,t388,t389,t39,t390
           complex(dp) ::  t391,t392,t393,t394,t395,t396,t397,t398,t399,t4,t40,t400
           complex(dp) ::  t401,t402,t403,t404,t405,t406,t407,t408,t409,t41,t410,t411
           complex(dp) ::  t412,t413,t414,t415,t416,t417,t418,t419,t42,t420,t421,t422
           complex(dp) ::  t423,t424,t425,t426,t427,t428,t429,t43,t430,t431,t432,t433
           complex(dp) ::  t434,t435,t436,t437,t438,t439,t44,t440,t441,t442,t443,t444
           complex(dp) ::  t445,t446,t447,t448,t449,t45,t450,t451,t452,t453,t454,t455
           complex(dp) ::  t456,t457,t458,t459,t46,t460,t461,t462,t463,t464,t465,t466
           complex(dp) ::  t467,t468,t469,t47,t470,t471,t472,t473,t474,t475,t476,t477
           complex(dp) ::  t478,t479,t48,t480,t481,t482,t483,t484,t485,t486,t487,t488
           complex(dp) ::  t489,t49,t490,t491,t492,t493,t494,t495,t5,t50,t51,t52
           complex(dp) ::  t53,t54,t55,t56,t57,t58,t59,t6,t60,t61,t62,t63
           complex(dp) ::  t64,t65,t66,t67,t68,t69,t7,t70,t71,t72,t73,t74
           complex(dp) ::  t75,t76,t77,t78,t79,t8,t80,t81,t82,t83,t84,t85
           complex(dp) ::  t86,t87,t88,t89,t9,t90,t91,t92,t93,t94,t95,t96
           complex(dp) ::  t97,t98,t99

           complex(dp) :: result

      t1 = intHLs160000x0111D0eps0()
      t2 = intHLs16s25s26s34s56x1112D4eps0()
      t3 = intHs160000x0112D2eps0()
      t4 = intHs160s26s34s56x1021D2eps1()
      t5 = intHs160s26s34s56x1022D4eps0()
      t6 = intHs160s26s34s56x1031D4eps0()
      t7 = intHs16s25s26s34s56x1121D4eps0()
      t8 = intHs16s25s26s34s56x1211D4eps0()
      t9 = 1.0_dp - epinv
      t10 = 4.0_dp
      t11 = propZ25 * s25
      t12 = -t11 + t10
      t13 = 3.0_dp
      t14 = gb ** 2.0_dp
      t15 = gw ** 2.0_dp
      t16 = t13 * s25
      t17 = t16 * t15 * propZ25
      t18 = t12 * t14 + t17
      t19 = intHs16s25s26s34s56x1311D4eps1()
      t20 = 2.0_dp
      t21 = -t20 + epinv
      t22 = -t12
      t23 = t14 * t22
      t24 = t23 - t17
      t25 = intHLs16s25s26s34s56x1111D2eps0()
      t26 = t11 + t20
      t27 = t14 * t26
      t28 = t27 + t17
      t29 = intHs16s25s26s34s56x1310D4eps0()
      t30 = s16 + s26 - s34 + s56
      t31 = intHLs16s25s26s34s56x1222D6eps1()
      t32 = s12 + s25
      t33 = t20 * t32
      t34 = t10 * s26
      t35 = t34 + t33 + s16
      t36 = intHLs16s25s26s34s56x1211D2eps1()
      t37 = intHs16s25s26s34s56x1222D6eps1()
      t38 = 6.0_dp
      t39 = 5.0_dp
      t40 = t13 * s26
      t41 = t38 * s12
      t42 = t39 * s25
      t43 = t42 + t41 + t40 + s34 - s56
      t44 = intHLs16s25s26s34s56x1121D2eps1()
      t45 = s25 + s26 + s56
      t46 = intHLs16s25s26s34s56x1112D2eps1()
      t47 = -s13 + s26 + s45 + s46 + s56
      t48 = t20 * t47 + s25
      t49 = t20 * s12
      t50 = t49 + s25
      t51 = intHs16s25s26s34s56x1321D6eps1()
      t52 = t16 + t49
      t53 = intHs16s25s26s34s56x1231D6eps1()
      t54 = t10 * s12
      t55 = t54 + t16
      t56 = intHLs16s25s26s34s56x1132D6eps1()
      t57 = t20 * s26
      t58 = t57 + s25
      t59 = intHLs16s25s26s34s56x1312D6eps1()
      t60 = s16 + s25
      t61 = t20 * (s12 + s26)
      t62 = t61 + t60
      t63 = intHs16s25s26s34s56x1312D6eps1()
      t64 = t16 + t61
      t65 = intHs16s25s26s34s56x1213D6eps1()
      t66 = s12 + s25 + s26
      t67 = t20 * t66 + s34
      t68 = s25 + s26 - s56
      t69 = intHLs16s25s26s34s56x1123D6eps1()
      t70 = t20 * s25
      t71 = t40 + t70 - s56
      t72 = intHLs16s25s26s34s56x1213D6eps1()
      t73 = t40 + t33 - s56
      t74 = intHs16s25s26s34s56x1132D6eps1()
      t75 = t54 + t70 + s26 + s34 - s56
      t76 = intHs16s25s26s34s56x1123D6eps1()
      t77 = t20 * (s12 + s34) + t68
      t78 = intHs160s26s34s56x1020D2eps1()
      t79 = intHs16s25s26s34s56x1110D2eps1()
      t80 = intHs16s25s26s34s56x1220D4eps0()
      t81 = intHs16s25s26s34s56x1130D4eps0()
      t82 = intHs16s25s26s34s56x1120D2eps1()
      t83 = t9 * s25
      t84 = intHs16s25s26s34s56x1211D2eps1()
      t85 = intHs16s25s26s34s56x1112D4eps0()
      t86 = 0.1e1_dp / s25
      t87 = t14 * (-propZ25 * t20 + 8.0_dp * t86) + t38 * t15 * propZ25
      t88 = intHLs16s25s26s34s56x1212D4eps0()
      t89 = s26 - s34 + s25
      t90 = s13 - s45 - s46
      t91 = t20 * t90
      t92 = t13 * s56
      t93 = t54 - t92 + t91 + t89
      t94 = intHs16s25s26s34s56x1222D6eps0()
      t95 = intHLs16s25s26s34s56x1121D2eps0()
      t96 = intHs16s25s26s34s56x1120D2eps0()
      t97 = intHs16s25s26s34s56x1210D2eps0()
      t98 = intHs16s25s26s34s56x1210D2eps1()
      t99 = intHs160000x0112D2eps1()
      t100 = intHs16s25s26s34s56x1121D2eps1()
      t101 = 1.0_dp + epinv
      t102 = t49 * t101
      t103 = intHLs160000x0112D2eps0()
      t104 = intHLs160000x0111D0eps1()
      t105 = intHLs160000x0112D2eps1()
      t106 = intHLs16s25s26s34s56x1112D4eps1()
      t107 = t20 * epinv
      t108 = t107 - t13
      t109 = intHs16s25s26s34s56x1122D4eps0()
      t110 = intHs16s25s26s34s56x1212D4eps0()
      t111 = intHLs16s25s26s34s56x1222D6eps0()
      t112 = intHs16s25s26s34s56x1212D4eps1()
      t113 = intHLs16s25s26s34s56x1211D2eps0()
      t114 = intHLs16s25s26s34s56x1212D4eps1()
      t115 = intHs16s25s26s34s56x1122D4eps1()
      t116 = t21 * s25
      t117 = intHLs16s25s26s34s56x1112D2eps0()
      t118 = intHLs16s25s26s34s56x1122D4eps1()
      t119 = t20 * s56
      t120 = t119 - t90
      t121 = intHs16s25s26s34s56x1121D4eps1()
      t122 = intHs16s25s26s34s56x1211D4eps1()
      t123 = intHs16s25s26s34s56x1112D4eps1()
      t124 = intHs160s26s34s56x1022D4eps1()
      t125 = intHs16s25s26s34s56x1411D6eps0()
      t126 = intHs16s25s26s34s56x1131D4eps0()
      t127 = intHs16s25s26s34s56x1221D4eps0()
      t128 = intHs16s25s26s34s56x1311D4eps0()
      t129 = intHs16s25s26s34s56x1220D4eps1()
      t130 = intHs16s25s26s34s56x1130D4eps1()
      t131 = intHs16s25s26s34s56x1131D4eps1()
      t132 = -t54 + t116
      t133 = intHs16s25s26s34s56x1221D4eps1()
      t134 = t116 - t49
      t135 = intHLs16s25s26s34s56x1132D6eps0()
      t136 = intHLs16s25s26s34s56x1312D6eps0()
      t137 = intHs16s25s26s34s56x1312D6eps0()
      t138 = intHs16s25s26s34s56x1113D4eps1()
      t139 = intHs16s25s26s34s56x1114D6eps0()
      t140 = intHs16s25s26s34s56x1213D6eps0()
      t141 = intHLs16s25s26s34s56x1113D4eps0()
      t142 = -t92 + t91 - s26
      t143 = intHLs16s25s26s34s56x1122D4eps0()
      t144 = intHLs16s25s26s34s56x1114D6eps0()
      t145 = intHLs16s25s26s34s56x1123D6eps0()
      t146 = intHLs16s25s26s34s56x1213D6eps0()
      t147 = intHs16s25s26s34s56x1132D6eps0()
      t148 = intHs16s25s26s34s56x1123D6eps0()
      t149 = intHs160s26s34s56x1020D2eps0()
      t150 = intHs16s25s26s34s56x1310D4eps1()
      t151 = intHs16s25s26s34s56x1110D2eps0()
      t152 = intHLs16s25s26s34s56x1113D4eps1()
      t153 = intHs16s25s26s34s56x1141D6eps0()
      t154 = intHs16s25s26s34s56x1321D6eps0()
      t155 = intHs16s25s26s34s56x1231D6eps0()
      t156 = intHs160s26s34s56x1031D4eps1()
      t157 = intHs16s25s26s34s56x1121D2eps0()
      t158 = 0.1e1_dp / t30
      t159 = t98 * t9
      t160 = t9 * t37
      t161 = t108 * t123
      t162 = (t121 + t122) * t21
      t163 = s34 * t36
      t164 = t46 * t48
      t165 = t9 * t31
      t166 = t100 * (t102 + s25)
      t167 = t9 * t124
      t168 = t99 * epinv
      t169 = t1 * t26
      t170 = -t7 - t8
      t171 = t106 * t108
      t172 = s34 * t113
      t173 = t35 * t111
      t174 = t48 * t117
      t175 = t82 * (t49 + t83)
      t176 = t21 * t79
      t177 = (-t109 - t110) * t58
      t178 = t115 * (epinv * t57 - s26 - s34 + s56 + t116 - t54)
      t179 = t9 * t78
      t180 = t86 * t14
      t181 = t9 * intHLs16s25s26s34s56x1114D6eps1()
      t182 = -t144 + t181
      t183 = t9 * intHs16s25s26s34s56x1141D6eps1()
      t184 = t183 - t153
      t185 = t9 * intHs16s25s26s34s56x1114D6eps1()
      t186 = intHs16s25s26s34s56x1411D6eps1() * t9
      t187 = (t185 - t139) * s34
      t188 = t184 * t50
      t189 = t182 * t28
      t68 = t18 * (t187 * t86 + t188 * t86 - t125 + t186) + (-t2 + t6) *
     # propZ25 * t15 + t189 * t86 * t68
      t190 = t9 * t59
      t191 = -t136 + t190
      t192 = t56 * t9
      t193 = t192 - t135
      t194 = t131 * t132
      t195 = t9 * t158
      t196 = t19 * t21
      t197 = t51 * t52
      t198 = t53 * t55
      t199 = t138 * s34
      t200 = t137 * t64
      t201 = t140 * t67
      t202 = t86 * t18
      t64 = t63 * t64
      t67 = t65 * t67
      t203 = t74 * t75
      t204 = t76 * t77
      t205 = t204 + t64 + t67 + t203
      t206 = t52 * t154
      t207 = t55 * t155
      t208 = t75 * t147
      t209 = t77 * t148
      t210 = t50 * t81
      t211 = t157 * s12
      t212 = t2 * t26
      t71 = t86 * (t18 * (-t206 - t207 - t208 - t209) + t24 * (-t158 * (
     #t32 * t80 + t210) - t211) + t28 * (-t145 * t71 - t146 * t73 - t152
     # * (-epinv * t142 + s25 + s26 - s56))) + t9 * (t24 * (t150 * t158 
     #+ t156 * t86) + t86 * (t69 * t71 + t72 * t73) * t28 + t202 * t205)
     # + t24 * (-t158 * t29 + t86 * (-t133 * t134 - t194 + t195 * (t129 
     #* t32 + t130 * t50)) - t126 - t127 - t128 - t196) + t202 * (t9 * (
     #t197 + t198) - t199 - t200 - t201) + t86 * (-t120 * t143 + t141 * 
     #t142 + t191 * t62 + t193 * t58 + t25) * t28 - t180 * (t22 * t6 + t
     #212)
      t15 = t15 * propZ25
      t68 = -t13 * t15 * (-t1 + t3 + t4 + t5 - t7 - t8) - t20 * t71 - t3
     #8 * t68 - t18 * (t158 * (-t149 * t86 + t159) + t86 * (t160 - t94) 
     #* t43 - t161 * t86 - t162 * t86) - t24 * (t86 * (-t151 * t158 - t1
     #66 + t167 - t168) - t84) + t85 * t87 + t86 * (t22 * t5 + t169) * t
     #14 - t86 * ((t44 * t45 - t104 + t163 + t164) * epinv + t45 * t95 +
     # t88 * t93 + t165 * t35) * t28 - (t179 * t86 - t96 - t97) * t158 *
     # t18 - t180 * (t12 * t170 + t22 * (-t3 - t4)) - (epinv * t105 + t1
     #14 * (epinv * t93 - s25 - t57) - t118 * (t20 * (epinv * t120 + s26
     #) + s25) + t103 - t171 + t172 - t173 + t174) * t86 * t28 - t86 * (
     #-t112 * ((-t13 + epinv) * s25 - t54 - t57 * t9) + (-t175 - t176) *
     # t158 + t177 - t178) * t24
      t71 = intHL0s25s26s34s56x1130D4eps0()
      t73 = -s13 + s25 + s26 + s56
      t93 = s45 + s46
      t120 = t20 * t93
      t142 = t120 + t73
      t213 = s13 - s45 - s56
      t214 = -t20 * t213 + s25 + s46
      t215 = intHLs160000x0211D2eps0()
      t216 = s16 - s25 - s34
      t217 = t20 * (s12 + s13 - s56)
      t218 = t13 * s46
      t219 = t10 * s45
      t220 = t219 + t218 - t217 - t216
      t221 = t10 * t93
      t217 = t221 - t217 - t216
      t222 = intHs160000x0121D2eps0()
      t223 = s16 - s25 + s26 - s34
      t224 = s13 - s56
      t225 = t13 * t224
      t226 = t221 - t225 - t49 - t223
      t227 = intHLs16s25s26s34s56x1211D4eps0()
      t228 = s26 - s25
      t229 = s16 - s34
      t230 = t20 * t229
      t231 = t38 * s45
      t232 = t39 * s46
      t233 = t232 + t231 - t54 - t225 - t230 - t228
      t234 = intHs16s25s26s34s56x1111D2eps0()
      t235 = s26 + s16 - s34 + s13 - s56
      t236 = t20 * (s12 - s45 - s46)
      t237 = t236 + t235
      t238 = intHL0s25s260s56x1012D2eps0()
      t239 = -t91 + t45
      t240 = intHL0s25s260s56x1013D4eps0()
      t241 = intHL0s25s260s56x1021D2eps0()
      t242 = intHs160000x0111D0eps0()
      t243 = intHL0s25s26s34s56x1210D2eps1()
      t244 = intHL0s25s26s34s56x1310D4eps1()
      t245 = intHs160s26s34s56x1013D4eps1()
      t246 = s26 + s16 - s34 - s56
      t247 = t20 * (s12 + s13 - s45 - s46)
      t248 = -t247 - t246
      t249 = t15 * s25
      t250 = t249 + t14
      t251 = s16 - s34 + s25
      t252 = t61 + t251
      t253 = s12 + s13 - s16
      t254 = s12 - s25
      t255 = t13 * s13
      t256 = t20 * s13
      t257 = -t256 + s25
      t258 = s16 - s34 - s13
      t259 = s12 + s13 + s16
      t260 = s12 - s26
      t261 = s13 - s45
      t262 = t20 * t261
      t263 = s12 + s13
      t264 = t263 * s13
      t265 = (-s34 + t259) * s25 - (s12 + s13 - s25) * s26 + (-t256 - s1
     #2 + s26 - s25 + s45) * s45 + (-t262 - t260 + s46) * s46 + (-t91 - 
     #t260 + s56) * s56 + t264
      t266 = (t20 * t258 + s45 - t254 + t40) * s45
      t267 = t257 * s34
      t268 = s56 ** 2.0_dp
      t269 = t13 * t268
      t270 = -t20 * (s12 + s26 + s13 - s45 - s46 - s56) - t229
      t271 = t10 * s56
      t272 = s46 ** 2.0_dp
      t273 = s45 + s46 + s56
      t274 = t10 * s56 * t93 + t13 * (s26 * t273 + t268) + t20 * (s45 * 
     #t229 + (s16 - s34 + s45 + s25) * s46 + s56 * t251) - s12 * (s26 - 
     #s13 + s45 - s46 + s56 + s25) + (-t20 * (s16 - s34 + s45 + s46) - s
     #25 - t40 - t271 + s13) * s13 + (s26 + s16 - s34 + s45) * s25 + s45
     # ** 2.0_dp + t272
      t275 = s26 + s56
      t276 = 7.0_dp
      t277 = t38 * t93
      t278 = t276 * s56
      t279 = s25 ** 2.0_dp
      t280 = t39 * s13
      t281 = t280 * s25
      t282 = intHL0s25s260s56x1021D2eps1()
      t283 = t83 + s26 + s56
      t284 = epinv * s56
      t285 = t9 * s26
      t286 = t275 * t273
      t287 = t9 * s34
      t288 = s26 + s16 + s56 - t287
      t246 = -t120 + t246
      t289 = t107 * s34
      t290 = t289 + t30
      t291 = t13 * epinv
      t292 = t291 - t10
      t293 = t20 * (s13 - s45 - s46 - s56)
      t294 = -t293 + s25
      t295 = intHs160000x0211D2eps1()
      t296 = intHL0s25s260s56x1031D4eps1()
      t297 = intHs16s25s26s34s56x1112D2eps0()
      t298 = s26 - s56
      t299 = intHLs16s25s26s34s56x1311D4eps1()
      t300 = s45 * t20 - s13 + s25 + s26 + s46 + s56
      t301 = s16 - s25
      t302 = t13 * s34
      t303 = -t302 + t301
      t304 = t120 + t45
      t305 = s26 * s34
      t306 = s13 * t303 - s16 * t304 + s25 * t45 + t10 * s34 * t93 + t20
     # * (s34 + s45 + s46) * s25 + t305 + t302 * s56
      t307 = intHL0s25s26s34s56x1220D4eps0()
      t308 = t20 * s34
      t309 = t308 + t45
      t310 = s25 + s26 + s34 + s56
      t311 = s13 * t309 - t304 * t310
      t312 = intHL0s25s26s34s56x1120D2eps1()
      t313 = epinv * s34
      t314 = -t313 + s26 + s56 + s25
      t315 = s13 - s25
      t316 = t20 * (s26 - s13)
      t317 = -s13 + s26 + s45
      t318 = t20 * t317
      t319 = s12 - s16 + s34 - s25
      t320 = t315 * s34
      t321 = s16 * s25
      t322 = intHLs16s25s26s34s56x1231D6eps1()
      t323 = s13 + s16
      t324 = t49 + t323
      t325 = s16 - s13
      t326 = s25 + s26
      t327 = t20 * t326
      t328 = s45 + s46 + s25
      t260 = t20 * t260
      t329 = t10 * t328
      t330 = (-t49 - t325 - s26) * s26
      t331 = (-t324 + s25) * s25
      t332 = s13 * s16
      t333 = s12 * s13
      t334 = t20 * (t93 * (t327 + s34) + t333) + t269 + (-t20 * t315 + s
     #26) * s34 + (t13 * (s34 - s13) - s16 - t260 + t329) * s56 + t330 +
     # t331 + t332
      t335 = t13 * s16
      t336 = t10 * t263
      t337 = s26 + s34
      t338 = t39 * s34
      t339 = t276 * s25
      t340 = (t49 + s16 + s13 + s25 + s26) * s26
      t341 = t10 * (t70 + t337) * t93 + t20 * (s16 * t90 + t279 + t333) 
     #+ t269 - (t336 + t335) * s25 + (-t10 * t315 + s26) * s34 + (-t13 *
     # t323 + t221 - t260 + t338 + t339) * s56 - t340
      t342 = s16 + s26 + s56
      t343 = -t302 + t342
      t344 = -t308 + t342
      t345 = s13 * t343 + t246 * t344 + t344 * t49
      t346 = t10 * s25
      t347 = -s13 * t276 + s26 + t346
      t348 = t38 * s16
      t349 = s26 + s45 + s46
      t350 = t10 * t349
      t351 = t276 * s34
      t228 = t20 * (t333 + t279) - t269 - (t335 + s13) * s25 - (-t280 + 
     #t60 + s26) * s26 + t347 * s34 + (-s16 * t39 + s25 + t255 - t350 + 
     #t351) * s56 + t93 * (-t10 * t228 + 8.0_dp * s34 - t348) + t280 * s
     #16
      t352 = t20 * s16
      t353 = s34 ** 2.0_dp
      t354 = t263 * s25
      t355 = intHLs16s25s26s34s56x1321D6eps1()
      t356 = s16 - s13 + s25
      t357 = t20 * t254
      t358 = t10 * s34
      t359 = t10 * s13
      t360 = s26 + s16
      t361 = t325 * s25
      t362 = (-t359 + t16 + t57) * s34
      t363 = -s12 - s16
      t364 = t363 * s13
      t365 = intHLs160000x0211D2eps1()
      t366 = intHs160000x0121D2eps1()
      t367 = -t308 + t301
      t303 = t20 * t367 * t93 - (t49 + s13 + s16) * s16 + s25 * t324 - s
     #26 * t301 - (-t20 * (s12 + s16 - s25) - s26 - t255 + s34) * s34 + 
     #s56 * t303
      t368 = -t363
      t369 = t20 * t368
      t370 = t280 - t70 - s16
      t371 = t13 * t325
      t325 = t70 * t325
      t372 = intHL0s25s26s34s56x1130D4eps1()
      t373 = intHs16s25s26s34s56x1113D4eps0()
      t374 = intHLs16s25s26s34s56x1411D6eps0()
      t375 = intHs160s26s34s56x1011D2eps0()
      t376 = t20 * (s12 + s13 - s45) + s26 + s16 - s34 - s46 - s56
      t377 = intHL0s25s260s56x1012D2eps1()
      t378 = intHL0s25s260s56x1013D4eps1()
      t379 = intHL0s25s260s56x1010D0eps0()
      t380 = intHL0s25s260s56x1011D2eps0()
      t381 = intHL0s25s260s56x1020D2eps0()
      t382 = intHL0s25s26s34s56x1110D2eps0()
      t383 = intHLs160000x0111D2eps0()
      t384 = s56 + s25
      t385 = t20 * t275
      t386 = t20 * t273
      t387 = t386 + s25
      t388 = t20 * t30
      t389 = intHL0s25s260s56x1022D4eps0()
      t390 = t70 + t275
      t391 = s13 * t390
      t386 = (t386 + s26 + s25) * s25 + t286 - t391
      t392 = t57 + t251
      t393 = propW16 * t252
      t394 = t393 * t45 * t30
      t395 = s12 + s16 + s26 - s34
      t396 = (s12 + s16 + s25 + s26) * s26
      t397 = t93 * t216
      t259 = t259 * s25
      t398 = (-s13 + s25 + s26) * s34
      t399 = s25 * t93
      t337 = t10 * t399 + t20 * (s34 * t273 + t279) + s13 * (-t16 + t230
     #) - s16 * t387 + s25 * t337 + t16 * s56
      t400 = s13 * (t119 + s25)
      t401 = t45 * t387
      t402 = -t401 + t400
      t403 = t388 + s25
      t404 = t49 * s13
      t405 = -t70 + t30
      t406 = intHs160s26s34s56x1012D2eps1()
      t407 = intHL0s25s260s56x1022D4eps1()
      t408 = intHs16s25s26s34s56x1112D2eps1()
      t409 = intHL0s25s26s34s56x1220D4eps1()
      t410 = s26 - s13 + s45 + s56
      t411 = s12 - s45
      t412 = t20 * t263
      t413 = s26 - s34
      t414 = -t20 * (s16 - s13 + s45) - t413
      t415 = t20 * t395
      t416 = (s12 - s13 + s16 + s25 + s26) * s26
      t417 = (t218 + t415) * s56
      t418 = t20 * (s25 * t368 + s45 * t216 + t272 + t364 - t398 + t416)
      t419 = (-t335 - t412 + s25) * s25
      t420 = t20 * ((t230 + s12 + s26 - s25) * s56 - t333 + t416)
      t421 = t10 * (t332 - t397)
      t422 = t359 - t16
      t423 = s26 + s16 - s34 - s13 + s45 + s46
      t424 = -t255 + t70 + s26
      t315 = (-t13 * t315 + s16 + s26 + t49) * s26
      t425 = t424 * s34
      t258 = (t13 * t258 + t20 * (s12 + s45 + s46) + s25 + t34) * s56
      t223 = t93 * t223
      t426 = t20 * (-t333 + t259 + t223)
      t427 = t13 * (-t332 + t268)
      t428 = s16 + s26 - s34
      t429 = -t20 * t423 - s25
      t430 = t316 + s25
      t416 = t20 * ((t20 * t360 + s12 - s13 + s25 - s34 + s45 + s46 + s5
     #6) * s56 + t93 * (t352 + t413) + t321 - t333 + t416)
      t431 = t359 * s16
      t432 = t93 * t428
      t433 = t268 + t432 + t364
      t434 = intHLs16s25s26s34s56x1221D4eps1()
      t435 = s25 + s26 + s34
      t436 = t40 + t70 - s13
      t437 = t20 * t411
      t438 = -t437 - s16 + s34 - s13 + s46 + s25
      t439 = t435 * s46
      t440 = (t412 + s16) * s25
      t441 = t20 * (-s45 * t435 + t333)
      t442 = t70 + s26
      t443 = -t20 * (-s34 * t326 - s45 * t442 + s56 * (s12 + s13 - s25 -
     # s45 + s46 - s56) + t272 + t396) - (t54 + t255 + t352) * s25 + (t2
     #62 + s25) * s46
      t444 = t20 * (s26 + s16 - s13 + s45 + s46)
      t445 = t93 * (-t308 + t360)
      t446 = -t20 * (t332 - t268 + t353 - t445) + s25 * (s26 + t323) + (
     #t336 + t40 + t352) * s34 + (-t302 + t444 + s25) * s56
      t447 = t20 * t257
      t448 = s34 - s25
      t323 = t20 * (-t332 + t268 + t432) + s25 * t323 + t257 * s26 - s34
     # * (-t412 - t360 + s34) + (t444 - t448) * s56
      t412 = intHLs16s25s26s34s56x1141D6eps0()
      t432 = s26 ** 2.0_dp
      t444 = t328 * s16
      t449 = s25 * s45
      t450 = s12 * t45
      t451 = intHLs16s25s26s34s56x1131D4eps1()
      t452 = intHLs160000x0121D2eps1()
      t453 = t393 * s16
      t454 = t20 * (-t453 + s25)
      t455 = s46 + s56
      t456 = t453 * t38
      t457 = t14 * (t26 * (t13 * t455 + s26 + t219 - t255 + t70) - t456)
     # - t17 * (t13 * (s13 - s46 - s56) - s26 - t219 - t454)
      t458 = intHLs16s25s26s34s56x1231D6eps0()
      t459 = t13 * t229
      t460 = s25 * s26
      t461 = intHs160000x0211D2eps0()
      t462 = intHL0s25s260s56x1031D4eps0()
      t463 = intHLs16s25s26s34s56x1211D4eps1()
      t464 = t10 * t224
      t465 = 8.0_dp * s45
      t466 = t276 * s46
      t467 = s26 - s16 + s34 + t277 - t225 - t357
      t468 = 8.0_dp * t93
      t469 = intHLs16s25s26s34s56x1121D4eps1()
      t231 = t39 * t455 + s26 + t16 + t231 - t359
      t470 = s56 * t39 + t277 - t357 - t359 - t428
      t232 = -t213 * t38 + t16 + t232
      t324 = (s26 - s34 + t324) * s25
      t262 = t20 * (-t449 + t272) + (-t262 + s26) * s46 - (-t218 + s25) 
     #* s56 + t324
      t471 = intHLs16s25s26s34s56x1121D4eps0()
      t472 = t58 * s34
      t473 = intHLs16s25s26s34s56x1321D6eps0()
      t474 = intHLs16s25s26s34s56x1131D4eps0()
      t475 = intHLs16s25s26s34s56x1221D4eps0()
      t476 = intHs16s25s26s34s56x1211D2eps0()
      t477 = intHLs160000x0122D4eps0()
      t478 = intHs160s26s34s56x1011D2eps1()
      t479 = intHLs160000x0121D2eps0()
      t454 = t10 * t213 - t218 - t454
      t218 = t26 * (t10 * (s45 + s56) + t218 + t70 - t359) - t456
      t480 = intHL0s25s260s56x1011D2eps1()
      t481 = intHL0s25s26s34s56x1110D2eps1()
      t482 = intHL0s25s260s56x1010D0eps1()
      t483 = intHL0s25s260s56x1020D2eps1()
      t484 = intHL0s25s26s34s56x1120D2eps0()
      t485 = intHL0s25s26s34s56x1210D2eps0()
      t486 = intHL0s25s26s34s56x1310D4eps0()
      t487 = intHs160s26s34s56x1013D4eps0()
      t488 = intHs160s26s34s56x1012D2eps0()
      t489 = intHs160s26s34s56x1021D2eps0()
      t490 = 0.1e1_dp / t45
      t491 = s34 * t243
      t492 = t86 * s34
      t493 = t490 * t86
      t494 = (t480 + t481) * t21
      t495 = t488 * t158
      t223 = t202 * (t133 * (-t10 * (-t332 + t268 + t223) + t404 + epinv
     # * (-t427 - t426 + t425 - t258 - t315) + (-t335 - t49 - s13 + s25)
     # * s25 + s26 * t422 - t422 * s34 + (-t10 * t423 + s25) * s56) + t1
     #42 * t242)
      t422 = 0.1e1_dp / ecossin ** 2.0_dp * propW34 * im
      t423 = epinv * t365
      t289 = t86 * (-t490 * (t307 * t311 + t312 * (s13 * (-t289 + t45) -
     # t304 * t314)) + t165 * t341 + t163 * (epinv * t214 + s25 - t293) 
     #+ t423 * t220) * t28
      t79 = t422 * (t24 * (t142 * t476 + t84 * (epinv * t142 - s13 - s16
     # - s26 + s34 + s56 - t236) + t86 * (epinv * (t217 * t99 + t226 * t
     #366) + s34 * (t248 * t489 * t158 + t294 * t297) + (t159 - t96 - t9
     #7) * t158 * t265 + t158 * (-t248 * t478 + t376 * t79) * t21)) + t8
     #6 * (t103 * (t14 * t218 - t17 * t454) + t105 * (t14 * (epinv * t21
     #8 + t456) + t249 * (-t291 * t454 + t456)) + t18 * (t158 * (t345 * 
     #t5 + (-t20 * (-s16 * t263 + (s12 + s13 + s16 + s26) * s34 + t268) 
     #+ (t13 * t263 + s26 + t352) * s26 + (-t318 + t319 - s46) * s46 + (
     #t10 * t261 + t13 * (s12 - s46) - t326) * s56 - t264 - t266 + s16 *
     #* 2.0_dp + t353 + t354) * (-t179 + t149) * t158) - t160 * t228) + 
     #t457 * t479) + t289)
      t99 = epinv * t104
      t218 = t101 * t408
      t73 = t109 * t446 + t110 * (t13 * (s56 * t229 + t460) + t20 * (s12
     # * t73 + s26 * t275 - s34 * t328 - t399 + t444) + s13 * (-t459 - t
     #57 + s25) - s25 * s56 + s26 * t229) + t158 * (t151 * t376 - t82 * 
     #(t256 * s16 - t269 + epinv * t265 + s25 * t253 + (t255 + t254) * s
     #26 - (t20 * (s16 - s34 - s13 + s45 + s25) + s12 + t40 + s46) * s46
     # + (t10 * t90 - t20 * t251 + s12 - t40) * s56 - t266 - t264 + t267
     #) - t375 * t248 - t287 * t248 * t406) + t218 * s34 * t294
      t261 = t122 * t292
      t263 = epinv * intHs160000x0111D0eps1()
      t41 = t18 * (t123 * (epinv * t467 + s16 - s34 - t16 + t260 + t464 
     #- t468) + t157 * (-t430 * s34 + t416 - t431) + t158 * (t4 * (s13 *
     # t290 + t246 * t288 + t288 * t49) - t167 * t345) + t228 * t94 + t2
     #61 * t294 + t263 * t142) + t24 * t73 + t28 * (-t111 * t341 - t143 
     #* t443 + t214 * (t172 - t99) - t463 * (epinv * t233 - s25 + t41 + 
     #t459 + t464 - t465 - t466 + t57) - t95 * (t20 * (-(s26 - s16 + s25
     #) * s45 + t259 - t268 - t332 + t396) + (t352 - t326) * s46 + (t20 
     #* (s12 + s16 + s13 - s45 - s25) - s46) * s56 - t472) - t469 * (epi
     #nv * t231 - t276 * t455 - s26 + t280 - t346 - t465) + t475 * (-s34
     # * t436 + (-t438 - s56) * s56 + t340 - t439 + t440 + t441) - t282 
     #* (t20 * (s13 * t283 - t286) - s25 * (t83 + t20 * (s45 * t9 + s46 
     #* t9) + t92 + t285 - t284)) * t490)
      t73 = t422 * t86
      t228 = t112 * (epinv * (-t20 * (-t333 + t397 + t432) - (t369 + s13
     #) * s25 - (t20 * (s12 - s13) + s16 + t16) * s26 + (-t459 - t61 + s
     #25) * s56 + t425 + t255 * s16) + t362 + t419 + t421 - t420) + t121
     # * (epinv * t470 + s16 + s26 - s34 - t16 - t278 + t280 - t468 + t4
     #9) + t467 * t85 + t470 * t7
      t246 = t115 * (-t10 * (-t364 - t445) + t269 + t325 + epinv * t446 
     #+ (-t370 + s26) * s26 - t347 * s34 + (-t351 + t350 + t371 + t70) *
     # s56) + t127 * (-t424 * s34 + t258 + t315 + t426 + t427)
      t258 = epinv * t452
      t213 = t18 * t228 + t24 * t246 + t28 * (-t106 * (epinv * t232 + 8.
     #0_dp * t213 - t346 - t466) - t117 * (t20 * (s46 * (t317 + s46) - t
     #449) - (-s46 * t20 + s25) * s56 + t324) - t118 * (t10 * (t442 * t9
     #3 + t268 + t279) - t281 + epinv * t443 + (t16 - t256) * s26 + (9.0
     #_dp * s25 + t10 * (s26 - s13 + s45 + s46)) * s56) - t2 * t232 - t4
     #34 * (t10 * t435 * t93 + t269 + t404 + epinv * (-t441 + s34 * t436
     # - (-t438 - s56) * s56 - t340 + t439 - t440) + (-t359 + t16 + s26)
     # * s34 + (t338 + t329 - t255 - t260 - s16) * s56 + t330 + t331 + t
     #332) + t44 * (epinv * (-t20 * (-s13 * (s16 - s25 - s56) + s26 * (s
     #16 - s34 - s45 + s25) + s56 * (s16 - s45 - s25) - t268 + t432 + t4
     #44 - t449 + t450) + s25 * s34 + s46 * t45) - t400 + t401) + t88 * 
     #(-t414 * s46 + t417 + t418) - t231 * t471) + t258 * t457
      t228 = t9 * t377
      t41 = t422 * (t213 * t86 + t28 * (t86 * (-t379 + t380 + t381 + t38
     #2) + t228) * t490 * t239) + t73 * t41 + t79 + t422 * (t28 * (t86 *
     # (-t1 * t214 - t114 * (t421 - t420 + epinv * (t414 * s46 - t417 - 
     #t418) + t362 + t419) + t215 * t220 - t227 * t233 + t25 * (t49 + s2
     #6 + s16 - s34 - s13 - s46 + s56 + s25) + t46 * (-epinv * (t20 * (s
     #25 * t411 + s46 * t410 + t272) + t235 * s25) + s25 * t142)) + (-t4
     #92 * (t484 + t485) - t238 - t241) * t490 * t239 + t493 * (t239 * (
     #t491 + t482 - t483) + t311 * t409) * t9 + t494 * t493 * t239) + (t
     #100 * (-t20 * t433 + epinv * (t430 * s34 - t416 + t431) - t257 * s
     #26 + s56 * t429 + t267 - t361) + t217 * t3 + t222 * t226 + t234 * 
     #(s13 + s26 - s56) + t495 * s34 * t248) * t86 * t24 + t223)
      t79 = t299 * s34
      t213 = t9 * t244
      t214 = t9 * t150
      t217 = t9 * t129
      t220 = t9 * t130
      t223 = t9 * t245
      t226 = t9 * t76
      t231 = t9 * t65
      t232 = -t128 * t294 + t86 * (t158 * (t265 * (-t214 + t29) + (-t217
     # + t80) * (t20 * t268 + s13 * t253 + (-t256 - t254) * s26 + (t316 
     #- s12 + s16 - s34 + s45) * s45 + (t318 + t251 + s46) * s46 + (-t13
     # * t90 - t319 + t57) * s56 + t320 + t321) - t220 * t274) - t237 * 
     #t6 + (t226 - t148) * (-t20 * (t93 * (-t302 + t360) + t364) - (-t25
     #5 + t60 + s26) * s26 - (t20 * t349 + s56 + t356 - t358) * s56 - t3
     #61 + t362) + t158 * (t223 - t487) * t248 * s34 + t231 * t306)
      t233 = t19 * t24
      t235 = t9 * t69
      t246 = -t235 + t145
      t248 = t9 * t72
      t253 = -t248 + t146
      t254 = s34 * intHLs16s25s26s34s56x1311D4eps0()
      t257 = t9 * t74
      t156 = t156 * t9
      t260 = (-t138 * (epinv * t270 + s26 + t221 - t225 + t70) - t270 * 
     #t373) * s34 - t137 * t337 - t140 * t306 + (t257 - t147) * (-t10 * 
     #t364 - t269 - (-t370 + s26) * s26 + (-t280 + t16 + s26) * s34 + (t
     #338 - t350 - t371 - t70) * s56 + t93 * (s34 * t38 - t10 * t360) - 
     #t325) + t81 * t274 * t158 + t156 * t237
      t264 = t295 * epinv
      t265 = -t131 * (-t20 * t429 * s56 + epinv * t323 + s26 * t447 - t4
     #47 * s34 + t10 * t433 + t325) - t294 * (t264 + t461)
      t219 = t73 * (t18 * t260 + t24 * t265 + t28 * (-t136 * t303 - t141
     # * t262 + t152 * (s25 * (-t255 + t219 + t70) - epinv * t262 + (t34
     #6 - s26) * s46 + (t16 + s46) * s56 + t460) + t246 * (t13 * (s46 * 
     #s56 + t279) + t20 * (-s13 * t275 + s26 * s45 + (s26 + s45 + s56) *
     # s56) + s25 * (t278 + t277 + s26) + s26 * s46 - t281) + (t296 * t9
     # - t462) * t294 + t458 * t334 + t253 * (-t20 * (-t332 + t320) + (-
     #s26 - t255 - t369 + s25) * s25 + (t346 - t230) * s45 - (-t346 + t2
     #30 + s26) * s46 + (t16 - t230 + s46) * s56) + t451 * (-epinv * (-t
     #20 * (s45 * t45 + t305 - t450) + s25 * (s16 - s34 - s46) + s26 * (
     #s16 - s46 + s25) + s56 * (s16 - s46 - s25) - t268 + t391 + t432) -
     # t20 * t402) + t254 * t300 + t372 * t9 * t142))
      t260 = t9 * t378
      t262 = t9 * t53
      t265 = t9 * t51
      t266 = t9 * t63
      t60 = t422 * (t28 * (t490 * (-t389 * t386 * t86 + t260 * t239) + t
     #86 * (-t192 + t135) * ((t16 + t385) * t387 - t256 * (t20 * t384 + 
     #s26))) - t86 * (t126 * t24 * t323 + epinv2 * (t14 * (-t13 * (t394 
     #+ t268) + t20 * (-t93 * (t40 + t70 + t229) + t332) + (-propZ25 * (
     #s56 * t395 + t259 - t332 + t396 + t397 - t398) + s16 + t336) * s25
     # + (s13 * t38 + s26 + t54 + t60) * s26 - (t256 + t326) * s34 + (t3
     #8 * t90 - t16 + t229 + t54 - t57) * s56) - t17 * (-s13 * t229 - (t
     #256 + s12) * s26 + (-t91 - s12 + s26 + s25 + s56) * s56 + t93 * t3
     #92 - t354 + t394)) * t158 * t490) + t202 * ((t262 - t155) * ((t10 
     #* t30 + t49) * s13 - t387 * (t388 - s25)) + (t265 - t154) * (s13 *
     # (-t16 + t388) - t387 * t405) + t266 * t337))
      t90 = t9 * intHLs160000x0113D4eps1() - intHLs160000x0113D4eps0()
      t256 = t9 * intHLs160000x0122D4eps1()
      t259 = intHLs160000x0111D2eps1() * t21
      t267 = s16 * (t256 - t477) - t383 - t259
      t144 = -t181 + t144
      t181 = t9 * intHLs16s25s26s34s56x1411D6eps1()
      t268 = t9 * intHLs16s25s26s34s56x1141D6eps1()
      t237 = t422 * (t86 * (t18 * ((s13 * t403 - t30 * t387 + t404) * (t
     #183 - t153) + (t185 - t139) * t142 * s34) + t28 * (t144 * (-s46 * 
     #t298 + t20 * s25 * (-s13 + t273) + t279) + (t268 - t412) * t402 + 
     #(t181 - t374) * t237 * s34)) + t18 * (t186 - t125) * t294 + t393 *
     # t86 * t267 * t250)
      t269 = 0.2e1_dp / 0.3e1_dp
      t270 = epinv * s25
      t272 = s26 + s16 - s34 + s56 - s25
      t273 = s25 * t101
      t274 = epinv * t30
      t277 = t176 + t151
      t278 = t390 * t407
      t279 = (-s34 * t406 - t78) * t9
      t280 = t313 * t4
      t281 = t287 * t124
      t286 = t12 * t222
      t287 = s34 * t297
      t288 = epinv * t366
      t289 = -t165 + t111
      t282 = t282 * t283
      t283 = t389 * t390
      t291 = t217 * t272
      t297 = t21 * t478
      t304 = t101 * t112
      t305 = t18 * t100
      t218 = t18 * (s34 * (t218 + t495) - t160 * t344) + t24 * (-t133 * 
     #(epinv * t405 + s16 + s26 - s34 + s56) + t158 * (t297 + t291) - t3
     #04 * t216) + t28 * (t310 * t289 + t434 * (-t274 + s34) + t490 * (s
     #34 * (-t307 + t484 + t485) + t21 * (-t480 - t481) + t9 * (-t482 + 
     #t483) - t282 - t283) - t163 * t101) + (t28 * t44 - t305) * t301 * 
     #epinv + epinv2 * (t14 * (s25 * (propZ25 * t216 + t10) + t275 * t38
     # + t230) + t17 * (t385 + t251)) * t490 * t158
      t230 = t82 * (t270 + s26 + s16 - s34 + s56)
      t251 = s34 * t5
      t275 = t80 * t272
      t306 = t106 * t292
      t311 = t86 * (t18 * (s34 * t489 + t230) + t24 * (t251 - t275)) * t
     #158
      t218 = t28 * (t86 * (t118 * (s26 + s56 + t273) - t258 + t306) - t2
     #28 * t490) + t86 * t218 + t86 * (t14 * (t22 * (t7 + t234) + t26 * 
     #t471 + t286) + t18 * (-t157 * t301 + t287 + t288) + t24 * (-t109 *
     # t367 - t110 * t216 + t115 * (-epinv * t367 + s34) - t121 * t9 - t
     #127 * t405 + t158 * (-t281 + t375) + t261) + t28 * (-t114 * t216 -
     # t30 * t475 + t301 * t95 - t469 * t9 + t490 * (t379 - t380 - t381 
     #- t382))) + t158 * (t18 * (t86 * (t149 + t279 + t280) - t159 + t96
     # + t97) + t86 * t277 * t24) + t28 * (t490 * (t238 + t241) + t86 * 
     #(-t479 - t172 - t25) + t143 + t86 * (s34 * t312 * t490 + t104) * e
     #pinv + t493 * (s34 * (t409 - t243) + t278) * t9) + t87 * (t81 + t4
     #61) + t169 * t180 + t311 + t202 * t94 * t344
      t228 = t125 - t186
      t243 = t228 * t24
      t261 = t352 * propW16
      t301 = t348 * propW16
      t311 = t14 * (t301 + t11 + t20) + t17 * (1.0_dp + t261)
      t315 = s25 + s26 - s34 + s56
      t316 = t38 * propW16
      t317 = epinv * t44
      t36 = (epinv * t36 + t113) * s34
      t113 = t103 * t311
      t311 = t105 * (epinv * t311 - t301 * t250)
      t319 = t7 + t8
      t320 = (t295 + t366) * epinv
      t321 = s34 * t408
      t323 = epinv2 * (t14 * (t316 * t30 - t10 + t11) + t17 * (t388 * pr
     #opW16 - 1.0_dp)) * t158
      t251 = t86 * (epinv * (s34 * t28 * t46 - t305 * t315) + t158 * (t1
     #8 * (-s34 * (t488 + t489) - t230) + t24 * (t275 - t297 - t251)) + 
     #t165 * t28 * t310) + t86 * (t14 * (-t12 * t461 - t22 * t319) + t18
     # * (-t157 * t315 - t320 - t321) + t24 * (-t115 * (epinv * t45 + s3
     #4) - t158 * t375 + t45 * (-t304 - t109 - t110) + t9 * (-t129 * t27
     #2 * t158 + t121 + t122)) + t28 * (s34 * t117 - t111 * t310 - t114 
     #* (epinv * t342 + s25 + s26 + s56) - t118 * (t83 + t313 + s26 + s5
     #6) - t143 * t448 - t342 * t88) + t323) + t158 * (t18 * (t86 * (-t1
     #49 - t279 - t280) + t159 - t96 - t97) + (-t151 - t176 + t281) * t8
     #6 * t24) - t81 * t87 + t86 * (-t14 * (t286 + t169) + t113 + t311) 
     #+ t86 * (t45 * (t95 + t317) - t171 + t25 - t99 + t36) * t28
      t223 = t487 - t223
      t279 = t248 - t146
      t235 = -t235 + t145
      t280 = (-t214 + t29) * t158
      t281 = propW16 * s16
      t286 = t281 * t250
      t304 = t286 * t86 * t90
      t191 = -t13 * t15 * (-t1 + t7 + t8 - t222 - t461) + t20 * (t24 * (
     #t280 + t126 + t127 + t128) + t86 * (t28 * (-s34 * t191 + t141 * t4
     #48 + t152 * (t83 + t313) - t193 * t45 + t216 * t279 + t390 * t235)
     # + s34 * (t158 * t223 + t138) * t24 - t220 * t18 + t212 * t14) + (
     #t133 + t19 + t131) * t24 * epinv) + t38 * ((s16 * (t477 - t479 - t
     #256 - t258) + t259 + t383) * t86 * t250 * propW16 + t15 * t2 - t18
     #9) - 12.0_dp * t304 - t251
      t193 = t20 + epinv
      t251 = s12 + s26 + s16 - s34 + s56
      t256 = t20 * t251
      t259 = t10 * t368
      t305 = t39 * t413 + t259 + t70 + t92
      t251 = t10 * t251
      t315 = t251 + s25
      t302 = t20 * (s12 + s26 + s16) - t302
      t324 = -t415 + s25
      t325 = -t256 + s25
      t328 = -t358 + t40 + t369 + t384
      t276 = t276 * t413 + t368 * t38 + s25 + t92
      t329 = t20 * t410 + t16
      t318 = t271 + t16 + t318
      t330 = -t20 * (s16 - s34 + s13 - s45 - s56) + s25 - t54
      t331 = t39 * s26
      t337 = t281 * t10
      t281 = 12.0_dp * t281
      t338 = -s13 + s25 + s26 + s45 + s56
      t341 = t310 * t475
      t342 = t463 + t469
      t344 = t310 * t434
      t122 = (t121 + t122) * t292
      t259 = t86 * (t14 * (-t12 * t3 - t26 * (t227 + t471)) + t18 * (-t1
     #58 * (t290 * t4 + t343 * t5) - t168 + t160 * t276) + t24 * (-t161 
     #- t122) + t28 * (epinv * (t329 * t46 - t344 + t365) + t342 * t9 + 
     #t215 - t306) + (t14 * (t281 + t11 + t20) + t17 * (1.0_dp + t337)) 
     #* (t258 + t479)) + t18 * (t86 * (t195 * t124 * t343 - t276 * t94) 
     #- t84) + t87 * (-t6 - t222 - t461) + t86 * (-t114 * (epinv * t330 
     #+ s25 + t57) + t117 * t329 - t118 * (epinv * t318 + s25 + t57) - t
     #143 * t318 - t330 * t88 + t35 * (t165 - t111) + t25 - t341) * t28 
     #+ t86 * (-t100 * (t20 * (-t101 * t363 + s26 - s34 + s56) + s25) - 
     #t112 * (-t10 * (s12 + s16 - s34) - t20 * t285 + t273) - t115 * (-t
     #20 * (-epinv * s26 + s25) - t259 + t270 + t351 - t92 - t331) + t17
     #7) * t24
      t273 = t56 * t58 + t59 * t62
      t276 = t179 - t149
      t285 = t58 * t135
      t62 = t62 * t136
      t290 = t45 * t451 + t79
      t230 = t86 * (-t230 * t18 + t275 * t24) * t158
      t36 = t18 * (t158 * (t492 * (t406 * t9 - t488 - t489) - t96 - t97)
     # + t86 * (-t140 * t302 - t147 * t305 - t155 * t315 + t156 - t321) 
     #- t86 * (t295 + t366) * epinv) + t86 * (-t28 * (t152 * (t107 * t33
     #8 + s25) + t254) + t113 + t311 + t323 - t23 * t85) + t24 * (t86 * 
     #(-t131 * (-t251 + t116) - t133 * (t270 - t256) - t154 * t325 - t15
     #8 * t277) - t126 - t127 - t128) + t86 * (t18 * (-t148 * t328 + t15
     #8 * t276) - t169 * t14) + t86 * (t273 * t9 + t45 * (-t474 + t95 + 
     #t317) + t36 - t99 - t285 - t62) * t28 + t195 * t98 * t18 + t86 * (
     #-t137 * t324 - t157 * t368 - t158 * (t375 + t297 + t291)) * t24 - 
     #t86 * t290 * t28 * epinv - t233 * t193 + t86 * (t18 * (t302 * t65 
     #+ t305 * t74 + t315 * t53 + t328 * t76) + t24 * (t324 * t63 + t325
     # * t51)) * t9 + t230
      t53 = (t46 - t114 - t118) * epinv
      t65 = epinv * t152 + t141
      t74 = s26 - s56 - s25
      t76 = s16 + s56
      t113 = t331 + t16 + t49 + t76
      t169 = t40 + t33 + s16 + s34 + s56
      t230 = t16 + t49 + s26 - s16 + s34 - s56
      t251 = t119 - s13 + s45 + s46 + s25
      t275 = t92 - t91 + t326
      t277 = t101 * s12
      t291 = t20 * (t277 + s26 + s25)
      t48 = t18 * (-t158 * t97 - t86 * (t108 * t121 + t110 * t252 + t123
     # * t292 + t133 * (t16 + t274 + t102) + t230 * t94 + t263)) + t24 *
     # (-t101 * t84 - t288 * t86) + t24 * (t86 * (-t109 * t58 - t158 * (
     #t151 + t175 + t176) - t166 + t167 - t178) - t476) - t7 * t87 + (t2
     #76 * t86 + t159 - t96) * t158 * t18 + t180 * (t22 * (-t4 - t5) + t
     #26 * (-t1 - t2)) + t86 * (epinv * (-t365 + t163 - t104 + t105) + t
     #172 + t174) * t28 + t86 * (t114 * (epinv * (-t92 + t91 + t413) + t
     #291) - t118 * (-t20 * (-epinv * t251 + s26) - t16) + t21 * t469 - 
     #t434 * (-t49 * t9 - s16 - s56 - t16 + t270 - t331) + t44 * (-s26 *
     # t193 + t284 - t83) - t46 * (-epinv * t48 + s25) - t74 * t95 + t88
     # * (-t92 + t247 + t413) + t9 * (t409 + t106) + t103 + t471) * t28 
     #- t180 * (t215 * t26 + t22 * t234)
      t83 = t137 + t154 + t19
      t27 = t86 * (t28 * (t113 * t458 - t136 * t216 - t141 * t275 - t143
     # * t251 - t152 * (epinv * t275 - s25 + s26 - s56) + t169 * t473 - 
     #t451 * (-t34 + t116) + t227 + t79) + (-t113 * t322 - t169 * t355 +
     # t372) * t28 * t9 - t28 * (t145 + t146) * t74 - t210 * t24 * t158 
     #+ t140 * t18 * t216 - t27 * t71) + t24 * (t86 * (t158 * (t32 * (t2
     #17 - t80) + t220 * t50) + t156 - t168 - t194 - t211) - t126 - t280
     #) + t28 * (t86 * (t9 * (t216 * t59 + t74 * (t69 + t72)) + t25) - t
     #474) + t202 * (t75 * (t257 - t147) + t77 * (t226 - t148) - t199 - 
     #t231 * t216) - t180 * (t3 + t6) * t22 + t9 * t18 * (t51 + t63) + t
     #23 * t83
      t12 = t10 * t86 * ((-t262 + t155) * t32 * t18 + t28 * t326 * (t192
     # - t135)) + t13 * (-t180 * t85 * t22 + t15 * (t1 + t2 - t4 - t5 + 
     #t215 - t222 - t234 + t307 + t312 + t242)) - t20 * t27 - t38 * (t86
     # * (t18 * (t187 + t188) + t28 * (s34 * (-t181 + t374) + t182 * t29
     #8 + t58 * (-t268 + t412))) + t15 * (-s25 * t83 + t3 + t6 - t71)) -
     # t86 * (t14 * (-t12 * t242 - t22 * t222 - t26 * (t307 + t312)) + t
     #18 * (-t112 * (epinv * t392 + t291) - t127 * (t49 + t30) + t160 * 
     #t230) + t28 * (t108 * t463 + t289 * (t346 + t40 + t49 + s56) - t47
     #5 * t50)) - t48 + 9.0_dp * t15 * t85
      t26 = s13 - s25 - s45 - s46 - s56
      t27 = -t20 * (s12 + s16 - s34 - s25) + s56 - t40
      t32 = s12 + s13 - s45 - s46 - s56 - s25
      t34 = s12 - s13 + s25 + s45 + s46 + s56
      t40 = s12 + s13 + s16 + s26 - s34 - s45 - s46
      t48 = t358 - t335 - t61 + s25
      t51 = t20 * t45
      t61 = t346 + t92 - s26
      t42 = t42 + t271 + t57
      t63 = 1.0_dp + t107
      t75 = t14 * (-t281 + t11 + t20)
      t77 = t17 * (1.0_dp - t337)
      t83 = t14 * (-t301 + t11 + t20) + t17 * (1.0_dp - t261)
      t87 = -t20 * t428 + s25
      t89 = t336 - t92 - t120 + t89
      t91 = t359 - t92 - t120 - s26
      t92 = t352 + t45
      t76 = t10 * t413 + t20 * t76 + s25
      t102 = t20 * (s26 + s45 + s46 + s56) + s25 - t359
      t94 = t94 * t43
      t44 = t73 * (t158 * (t18 * (-t82 * (-t26 * t9 + t277) - t32 * (t96
     # + t97)) + t24 * (t176 + t151)) + t24 * (t100 * (-t20 * (-t277 + t
     #274) + s25) + t109 * t76 + t115 * (epinv * t76 - s26 - s34 + s56 -
     # t54 - t70)) + t28 * (-t117 * t102 - t118 * (t10 * (epinv * s13 + 
     #s56) + t20 * (-epinv * (t119 + t93) + s26) + t339) - t44 * (epinv 
     #* t92 - t51) - t46 * (epinv * t102 - t70) - t95 * t92 + (s25 * t38
     # + t10 * (s34 + s56) - t335 - t49) * (t165 - t111) - t171 - t163 *
     # t21) - t160 * t18 * t43 + t149 * t18 * t40 * t158 ** 2.0_dp)
      t46 = t422 * (t18 * (t86 * (-t9 * (t204 + t197 + t198 + t64 + t67 
     #+ t203) + t199 * t63) + t128) + t24 * (t86 * (-t80 * t26 * t158 + 
     #t313 * t408) + t476) + t86 * (t490 * (-s34 * t485 + t307 * t310 + 
     #t312 * t314 + t282 + t283 - t379 + t380) + t9 * (t27 * t72 + t42 *
     # t56 + t48 * t59 + t61 * t69)) * t28)
      t11 = t73 * (t24 * (t126 * t403 + t127 * t30 + t131 * (epinv * t40
     #3 - t54 - t70)) + t28 * (-t135 * t42 - t145 * t61 - t146 * t27 - t
     #152 * (epinv * t91 - s26 + s56 + t16) + t490 * (t381 + t382)) + t1
     #33 * t18 * (-t274 + t33) - epinv2 * (t14 * (t316 * t45 + t11 + t20
     #) + t17 * (t51 * propW16 + 1.0_dp)) * t490)
      t11 = t11 + t46 + t422 * (t28 * (t490 * (t86 * (-s34 * t484 + t9 *
     # (t482 - t483) + t494) - t238 - t241) + t86 * (t423 - t2 - t25 + t
     #215) - t86 * t342 * t108) + t86 * (t263 * t18 + t479 * t83) + t86 
     #* (t158 * (t32 * (t214 - t29) + t34 * t81 + t217 * t26) - t156 + t
     #320) * t24) + t422 * (t18 * (t19 * t193 + t86 * (t209 + t242 + t20
     #0 + t201 + t206 + t207 + t208)) + t28 * (-t86 * (t136 * t48 + t141
     # * t91 + t143 * (t20 * t224 - t93)) + (t86 * (-t310 * t409 - t278 
     #+ t491) + t377) * t490 * t9) + t86 * (t24 * (t157 * (s12 - s16 - s
     #26 + s34 - s56) + t222 + t287 + t461 + t6 - t195 * t130 * t34) + t
     #258 * t83))
      t14 = t422 * (t28 * (t490 * (t492 * (-t213 + t486) + t240 - t260) 
     #+ t86 * (-t309 * t473 + t9 * (t309 * t355 - t296 - t372) - t227 - 
     #t344 + t462 - t471 + t71) + t86 * (t322 * t9 - t458) * (t51 + s34)
     #) + t202 * (s34 * t373 + t85))
      t11 = t10 * t422 * (t18 * (t492 * (-t185 + t139) + t125 - t186 + t
     #86 * (-t183 + t153) * t50) + t86 * t144 * t74 * t28) + t20 * t73 *
     # (t18 * t319 + t24 * (t168 + t3)) + t269 * (t44 + t422 * (t86 * (t
     #103 * (t77 + t75) + t105 * (epinv * (t77 + t75) + 12.0_dp * t286) 
     #+ t24 * (-t167 + t4 + t5)) + t202 * (t110 * t87 + t112 * (epinv * 
     #t87 + t16 + t54 + t57) + t123 * (epinv * t10 - t39) + t122 + t94 +
     # t195 * (-t78 * t40 * t158 + t32 * t98)) + t84 * t63 * t24 + (-t11
     #4 * (epinv * t89 - t10 * t229 + t16 - t57) - t88 * t89 + t1 - t172
     # + t99) * t86 * t28)) + 0.4e1_dp / 0.3e1_dp * t11 + 0.8e1_dp / 0.3
     #e1_dp * t14 - 0.16e2_dp / 0.3e1_dp * t73 * t28 * t290 - 16.0_dp * 
     #t286 * t73 * t90 - 8.0_dp * t73 * (t28 * (s34 * (-t181 + t374) + t
     #45 * (-t268 + t412)) + t267 * t250 * propW16)
      t14 = t54 + t293 + s25
      t17 = epinv * t58
      t26 = t293 + s25
      t6 = t6 - t85
      t14 = t18 * (-t112 * (t54 + t16 - t17 + t57) - t115 * (t54 + t70 -
     # t17 + s26 + s34 - s56) - t177 - t94) + t28 * (-epinv * (t452 + t3
     #65) - t106 * t21 - t114 * (-epinv * t14 + s25 + t57) - t118 * (-ep
     #inv * t26 + s25 + t57) + t14 * t88 + t143 * t26 - t463 * t9 - t173
     # + t174 + t341) + t18 * (t9 * (t37 * t43 - t124) - t161 - t162 + t
     #166 - t7 - t8) - t24 * (s25 * t84 + t168) + t28 * (epinv * (t344 +
     # t164) + t9 * (t31 * t35 - t469) - t2 - t215 + t227 + t25 + t471 -
     # t479) + t23 * (-t3 - t4 - t5)
      t4 = struc60PM * (-t10 * t28 * t65 + t20 * (t28 * (-t143 + t53 - t
     #88 + t117) + (-t9 * (t150 + t129 + t130) + t29 + t80 + t81) * t158
     # * t18) - t158 * (t24 * (t9 * (t82 + t98) - t96 - t97) + (t179 - t
     #149) * t158 * t18)) + struc9PM * (t10 * t28 * (t141 * t47 + t246 *
     # t326 + t253 * t66) - t13 * t249 * (t3 + t4 + t5) + t20 * (t9 * (-
     #t18 * t205 - t273 * t28) + t18 * (s25 * (-t126 - t127 - t128) - t1
     #31 * t132 - t133 * t134 + t52 * (t154 - t265) + t55 * (-t262 + t15
     #5) + t156 + t199 + t200 + t201 + t208 + t209) + t28 * (-t45 * t474
     # - t254 + t285 + t62) + t211 * t24 + t23 * t6 + t152 * t28 * (t107
     # * t47 + s25) - t196 * s25 * t18 - t290 * t28 * epinv) + t38 * (t1
     #8 * (s25 * t228 - t187 - t188) - t189 * s25 - t249 * t6) - t14)
      t3 = t422 * (t68 * struc10PM + t86 * t4 + struc5PM * (t10 * (t86 *
     # (t18 * (-t220 + t81) + t28 * (t141 * t338 + t235 * t326 - t279 * 
     #t66)) + (t492 * t223 - t214 + t29) * t158 * t24) + t13 * (t15 * (t
     #471 + t3 + t227) + t180 * (-t170 * t22 + t212)) - t20 * t36 - 24.0
     #_dp * t304 + t38 * (-t202 * t184 * (t256 + s25) + t15 * (t1 - t85)
     # + t492 * (t138 + t139 - t185) * t24 + t243 - t189) - t259 - 12.0_
     #dp * t86 * t250 * propW16 * t267 + 9.0_dp * t15 * (t2 - t7 - t8)) 
     #+ struc6PM * t12)
      t4 = (-t20 * t65 + t117 - t143 + t53 - t88) * t28 * struc61PM * t8
     #6
      t1 = t422 * (struc2PM * (-t13 * (t15 * (-t471 - t1 + t7 - t222 + t
     #234) - t180 * (t22 * t8 + t212)) + t20 * (t86 * (t24 * (s34 * (t10
     #1 * t138 + t158 * t487 + t140 - t231 + t373) + t216 * (t266 - t137
     #) + t272 * (-epinv * t131 - t126 + t265)) + t28 * (s34 * (-t486 * 
     #t490 + t136 - t190 + t299) + t216 * (t248 - t146) + t9 * (-t390 * 
     #t69 - t45 * t56 + t296)) - t262 * t18 * t30) + t24 * (t101 * t19 +
     # t158 * (-t9 * (t492 * t245 + t150) + t29) + t128 - t154 * t272 * 
     #t86) + t28 * (epinv * (-t105 * t86 + t451) + t490 * (t9 * (t492 * 
     #t244 + t378) - t240) + t86 * (t135 * t45 + t145 * t390 - t103 - t4
     #62) + t152 + t474) + t202 * (t155 * t30 - t220 + t264)) + t38 * (t
     #243 - t189) + t218 + 9.0_dp * t15 * (t2 - t8)) + t191 * (struc3PM 
     #+ struc8PM) + t4)
      result = t269 * t3 + (-t10 * t73 * t453 * t250 * t90 - t20 * t237 
     #+ t41 / 3.0_dp - (t60 + t219 + t422 * (t18 * t232 + t28 * (t490 * 
     #(-t239 * t240 + t86 * (t9 * t386 * t407 + (t213 - t486) * t239 * s
     #34)) + t86 * (-t142 * t71 - t474 * (-s46 * t326 + s56 * (s16 + s13
     # - s46 - s25 + t437 - s56) + t340 + t440 - t472 - t327 * s45) + t9
     # * (t303 * t59 - t322 * t334) + t79 * (epinv * t300 - s16 + s25 - 
     #s26 + s34 + t221 - t225 - t49)) + t86 * (t355 * t9 - t473) * (-t20
     # * (t93 * (t308 + t326) + t333) + (t49 + s16) * s25 - (-t49 - t356
     # - s26) * s26 - (-t255 - s16 - t357 + s34) * s34 - (t358 - t236 - 
     #s16 - s13 + s25 + s56) * s56 - t332)) + t233 * (epinv * t294 - s16
     # + s25 - s26 + s34 + t221 - t225 - t49))) * t269 + t422 * t202 * t
     #8 * t294) * struc1PM + t11 * struc7PM + 0.4e1_dp / 0.3e1_dp * t1

           ampNonresonantLightFullImC4PM = result
       end function ampNonresonantLightFullImC4PM

       function ampNonresonantLightFullImC4PP()
           implicit none
           complex(dp) :: ampNonresonantLightFullImC4PP

           complex(dp) :: result

      result = 0.0_dp

           ampNonresonantLightFullImC4PP = result
       end function ampNonresonantLightFullImC4PP

       function ampNonresonantLightFullImC7MM()
           implicit none
           complex(dp) :: ampNonresonantLightFullImC7MM

           complex(dp) :: result

      result = 0.0_dp

           ampNonresonantLightFullImC7MM = result
       end function ampNonresonantLightFullImC7MM

       function ampNonresonantLightFullImC7MP()
           implicit none
           complex(dp) :: ampNonresonantLightFullImC7MP

           complex(dp) :: result

      result = 0.0_dp

           ampNonresonantLightFullImC7MP = result
       end function ampNonresonantLightFullImC7MP

       function ampNonresonantLightFullImC7PM()
           implicit none
           complex(dp) :: ampNonresonantLightFullImC7PM

           complex(dp) :: result

      result = 0.0_dp

           ampNonresonantLightFullImC7PM = result
       end function ampNonresonantLightFullImC7PM

       function ampNonresonantLightFullImC7PP()
           implicit none
           complex(dp) :: ampNonresonantLightFullImC7PP

           complex(dp) :: result

      result = 0.0_dp

           ampNonresonantLightFullImC7PP = result
       end function ampNonresonantLightFullImC7PP

       function ampNonresonantLightFullMM()
           implicit none
           complex(dp) :: ampNonresonantLightFullMM
           complex(dp) ::  t1,t10,t100,t1000,t1001,t1002,t1003,t1004,t1005,t1006,t1007,t1008
           complex(dp) ::  t1009,t101,t1010,t1011,t1012,t1013,t1014,t1015,t1016,t1017,t1018,t1019
           complex(dp) ::  t102,t1020,t1021,t1022,t1023,t1024,t1025,t1026,t1027,t1028,t1029,t103
           complex(dp) ::  t1030,t1031,t1032,t1033,t1034,t1035,t1036,t1037,t1038,t1039,t104,t1040
           complex(dp) ::  t1041,t1042,t1043,t1044,t1045,t1046,t1047,t1048,t1049,t105,t1050,t1051
           complex(dp) ::  t1052,t1053,t1054,t1055,t1056,t1057,t1058,t1059,t106,t1060,t1061,t1062
           complex(dp) ::  t1063,t1064,t1065,t1066,t1067,t1068,t1069,t107,t1070,t1071,t1072,t1073
           complex(dp) ::  t1074,t1075,t1076,t1077,t1078,t1079,t108,t1080,t1081,t1082,t1083,t1084
           complex(dp) ::  t1085,t1086,t1087,t1088,t1089,t109,t1090,t1091,t1092,t1093,t1094,t1095
           complex(dp) ::  t1096,t1097,t1098,t1099,t11,t110,t1100,t1101,t1102,t1103,t1104,t1105
           complex(dp) ::  t1106,t1107,t1108,t1109,t111,t1110,t1111,t1112,t1113,t1114,t1115,t1116
           complex(dp) ::  t1117,t1118,t1119,t112,t1120,t1121,t1122,t1123,t1124,t1125,t1126,t1127
           complex(dp) ::  t1128,t1129,t113,t1130,t1131,t1132,t1133,t1134,t1135,t1136,t1137,t1138
           complex(dp) ::  t1139,t114,t1140,t1141,t1142,t1143,t1144,t1145,t115,t116,t117,t118
           complex(dp) ::  t119,t12,t120,t121,t122,t123,t124,t125,t126,t127,t128,t129
           complex(dp) ::  t13,t130,t131,t132,t133,t134,t135,t136,t137,t138,t139,t14
           complex(dp) ::  t140,t141,t142,t143,t144,t145,t146,t147,t148,t149,t15,t150
           complex(dp) ::  t151,t152,t153,t154,t155,t156,t157,t158,t159,t16,t160,t161
           complex(dp) ::  t162,t163,t164,t165,t166,t167,t168,t169,t17,t170,t171,t172
           complex(dp) ::  t173,t174,t175,t176,t177,t178,t179,t18,t180,t181,t182,t183
           complex(dp) ::  t184,t185,t186,t187,t188,t189,t19,t190,t191,t192,t193,t194
           complex(dp) ::  t195,t196,t197,t198,t199,t2,t20,t200,t201,t202,t203,t204
           complex(dp) ::  t205,t206,t207,t208,t209,t21,t210,t211,t212,t213,t214,t215
           complex(dp) ::  t216,t217,t218,t219,t22,t220,t221,t222,t223,t224,t225,t226
           complex(dp) ::  t227,t228,t229,t23,t230,t231,t232,t233,t234,t235,t236,t237
           complex(dp) ::  t238,t239,t24,t240,t241,t242,t243,t244,t245,t246,t247,t248
           complex(dp) ::  t249,t25,t250,t251,t252,t253,t254,t255,t256,t257,t258,t259
           complex(dp) ::  t26,t260,t261,t262,t263,t264,t265,t266,t267,t268,t269,t27
           complex(dp) ::  t270,t271,t272,t273,t274,t275,t276,t277,t278,t279,t28,t280
           complex(dp) ::  t281,t282,t283,t284,t285,t286,t287,t288,t289,t29,t290,t291
           complex(dp) ::  t292,t293,t294,t295,t296,t297,t298,t299,t3,t30,t300,t301
           complex(dp) ::  t302,t303,t304,t305,t306,t307,t308,t309,t31,t310,t311,t312
           complex(dp) ::  t313,t314,t315,t316,t317,t318,t319,t32,t320,t321,t322,t323
           complex(dp) ::  t324,t325,t326,t327,t328,t329,t33,t330,t331,t332,t333,t334
           complex(dp) ::  t335,t336,t337,t338,t339,t34,t340,t341,t342,t343,t344,t345
           complex(dp) ::  t346,t347,t348,t349,t35,t350,t351,t352,t353,t354,t355,t356
           complex(dp) ::  t357,t358,t359,t36,t360,t361,t362,t363,t364,t365,t366,t367
           complex(dp) ::  t368,t369,t37,t370,t371,t372,t373,t374,t375,t376,t377,t378
           complex(dp) ::  t379,t38,t380,t381,t382,t383,t384,t385,t386,t387,t388,t389
           complex(dp) ::  t39,t390,t391,t392,t393,t394,t395,t396,t397,t398,t399,t4
           complex(dp) ::  t40,t400,t401,t402,t403,t404,t405,t406,t407,t408,t409,t41
           complex(dp) ::  t410,t411,t412,t413,t414,t415,t416,t417,t418,t419,t42,t420
           complex(dp) ::  t421,t422,t423,t424,t425,t426,t427,t428,t429,t43,t430,t431
           complex(dp) ::  t432,t433,t434,t435,t436,t437,t438,t439,t44,t440,t441,t442
           complex(dp) ::  t443,t444,t445,t446,t447,t448,t449,t45,t450,t451,t452,t453
           complex(dp) ::  t454,t455,t456,t457,t458,t459,t46,t460,t461,t462,t463,t464
           complex(dp) ::  t465,t466,t467,t468,t469,t47,t470,t471,t472,t473,t474,t475
           complex(dp) ::  t476,t477,t478,t479,t48,t480,t481,t482,t483,t484,t485,t486
           complex(dp) ::  t487,t488,t489,t49,t490,t491,t492,t493,t494,t495,t496,t497
           complex(dp) ::  t498,t499,t5,t50,t500,t501,t502,t503,t504,t505,t506,t507
           complex(dp) ::  t508,t509,t51,t510,t511,t512,t513,t514,t515,t516,t517,t518
           complex(dp) ::  t519,t52,t520,t521,t522,t523,t524,t525,t526,t527,t528,t529
           complex(dp) ::  t53,t530,t531,t532,t533,t534,t535,t536,t537,t538,t539,t54
           complex(dp) ::  t540,t541,t542,t543,t544,t545,t546,t547,t548,t549,t55,t550
           complex(dp) ::  t551,t552,t553,t554,t555,t556,t557,t558,t559,t56,t560,t561
           complex(dp) ::  t562,t563,t564,t565,t566,t567,t568,t569,t57,t570,t571,t572
           complex(dp) ::  t573,t574,t575,t576,t577,t578,t579,t58,t580,t581,t582,t583
           complex(dp) ::  t584,t585,t586,t587,t588,t589,t59,t590,t591,t592,t593,t594
           complex(dp) ::  t595,t596,t597,t598,t599,t6,t60,t600,t601,t602,t603,t604
           complex(dp) ::  t605,t606,t607,t608,t609,t61,t610,t611,t612,t613,t614,t615
           complex(dp) ::  t616,t617,t618,t619,t62,t620,t621,t622,t623,t624,t625,t626
           complex(dp) ::  t627,t628,t629,t63,t630,t631,t632,t633,t634,t635,t636,t637
           complex(dp) ::  t638,t639,t64,t640,t641,t642,t643,t644,t645,t646,t647,t648
           complex(dp) ::  t649,t65,t650,t651,t652,t653,t654,t655,t656,t657,t658,t659
           complex(dp) ::  t66,t660,t661,t662,t663,t664,t665,t666,t667,t668,t669,t67
           complex(dp) ::  t670,t671,t672,t673,t674,t675,t676,t677,t678,t679,t68,t680
           complex(dp) ::  t681,t682,t683,t684,t685,t686,t687,t688,t689,t69,t690,t691
           complex(dp) ::  t692,t693,t694,t695,t696,t697,t698,t699,t7,t70,t700,t701
           complex(dp) ::  t702,t703,t704,t705,t706,t707,t708,t709,t71,t710,t711,t712
           complex(dp) ::  t713,t714,t715,t716,t717,t718,t719,t72,t720,t721,t722,t723
           complex(dp) ::  t724,t725,t726,t727,t728,t729,t73,t730,t731,t732,t733,t734
           complex(dp) ::  t735,t736,t737,t738,t739,t74,t740,t741,t742,t743,t744,t745
           complex(dp) ::  t746,t747,t748,t749,t75,t750,t751,t752,t753,t754,t755,t756
           complex(dp) ::  t757,t758,t759,t76,t760,t761,t762,t763,t764,t765,t766,t767
           complex(dp) ::  t768,t769,t77,t770,t771,t772,t773,t774,t775,t776,t777,t778
           complex(dp) ::  t779,t78,t780,t781,t782,t783,t784,t785,t786,t787,t788,t789
           complex(dp) ::  t79,t790,t791,t792,t793,t794,t795,t796,t797,t798,t799,t8
           complex(dp) ::  t80,t800,t801,t802,t803,t804,t805,t806,t807,t808,t809,t81
           complex(dp) ::  t810,t811,t812,t813,t814,t815,t816,t817,t818,t819,t82,t820
           complex(dp) ::  t821,t822,t823,t824,t825,t826,t827,t828,t829,t83,t830,t831
           complex(dp) ::  t832,t833,t834,t835,t836,t837,t838,t839,t84,t840,t841,t842
           complex(dp) ::  t843,t844,t845,t846,t847,t848,t849,t85,t850,t851,t852,t853
           complex(dp) ::  t854,t855,t856,t857,t858,t859,t86,t860,t861,t862,t863,t864
           complex(dp) ::  t865,t866,t867,t868,t869,t87,t870,t871,t872,t873,t874,t875
           complex(dp) ::  t876,t877,t878,t879,t88,t880,t881,t882,t883,t884,t885,t886
           complex(dp) ::  t887,t888,t889,t89,t890,t891,t892,t893,t894,t895,t896,t897
           complex(dp) ::  t898,t899,t9,t90,t900,t901,t902,t903,t904,t905,t906,t907
           complex(dp) ::  t908,t909,t91,t910,t911,t912,t913,t914,t915,t916,t917,t918
           complex(dp) ::  t919,t92,t920,t921,t922,t923,t924,t925,t926,t927,t928,t929
           complex(dp) ::  t93,t930,t931,t932,t933,t934,t935,t936,t937,t938,t939,t94
           complex(dp) ::  t940,t941,t942,t943,t944,t945,t946,t947,t948,t949,t95,t950
           complex(dp) ::  t951,t952,t953,t954,t955,t956,t957,t958,t959,t96,t960,t961
           complex(dp) ::  t962,t963,t964,t965,t966,t967,t968,t969,t97,t970,t971,t972
           complex(dp) ::  t973,t974,t975,t976,t977,t978,t979,t98,t980,t981,t982,t983
           complex(dp) ::  t984,t985,t986,t987,t988,t989,t99,t990,t991,t992,t993,t994
           complex(dp) ::  t995,t996,t997,t998,t999

           complex(dp) :: result

      t1 = intHs160s26s34s56x1031D4eps1()
      t2 = 1.0_dp - epinv
      t3 = 8.0_dp
      t4 = gw ** 2.0_dp
      t5 = gb ** 2.0_dp
      t6 = t5 ** 2.0_dp
      t7 = 9.0_dp * t4 ** 2.0_dp
      t8 = propZ25 * s25
      t9 = -t3 * t5 * t4 + t8 * (-t7 + t6)
      t10 = intHs16s25s26s34s56x1211D4eps0()
      t11 = -t9
      t12 = s16 + s26 - s34 - s45
      t13 = -s13 + s25 + s26 + s45 + s46 + s56
      t14 = s12 + s16 + s26 - s34 - s45
      t15 = intHLs16s25s26s34s56x1112D4eps0()
      t16 = 6.0_dp
      t17 = 4.0_dp
      t18 = t8 * t16
      t19 = t17 + t18
      t20 = t8 * t5
      t21 = t7 * t8
      t22 = t5 * (t19 * t4 + t20) + t21
      t23 = intHs16s25s26s34s56x1111D2eps0()
      t24 = s12 + s16 + s26 + s46 + s56
      t25 = intHs160s26s34s56x1020D2eps1()
      t26 = s12 + s16 + s26 - s34 - s45 + s56
      t27 = s16 + s26 - s34 + s56
      t28 = 2.0_dp
      t29 = -t28 + epinv
      t30 = intHs160s26s34s56x1013D4eps0()
      t31 = intHs160000x0211D2eps0()
      t32 = intHs160000x0211D2eps1()
      t33 = intHs16s25s26s34s56x1211D4eps1()
      t34 = intHLs16s25s26s34s56x1112D4eps1()
      t35 = t28 + epinv
      t36 = intHLs16s25s26s34s56x1211D2eps1()
      t37 = t2 * s12
      t38 = intHs160s26s34s56x1012D2eps1()
      t39 = intHs160s26s34s56x1013D4eps1()
      t40 = intHs16s25s26s34s56x1112D2eps1()
      t41 = epinv * t14
      t42 = t41 - s12
      t43 = intHL0s25s26s34s56x1210D2eps1()
      t44 = t28 * t5 * t4
      t19 = t8 * (t7 + t6) + t44 * t19 / 2.0_dp
      t45 = -s13 + s25 + s26 + s45 + s46
      t46 = t28 * s56
      t47 = -t46 - t45
      t48 = s34 + s45 + s46 + s56
      t49 = t13 * t48
      t50 = s12 * t47 + t49
      t51 = s25 + s26 + s56
      t52 = intHL0s25s26s34s56x1310D4eps1()
      t53 = intHs16s25s26s34s56x1120D2eps0()
      t54 = intHs16s25s26s34s56x1210D2eps0()
      t55 = intHs16s25s26s34s56x1310D4eps0()
      t56 = intHL0s25s260s56x1031D4eps0()
      t57 = intHs16s25s26s34s56x1130D4eps0()
      t58 = intHLs16s25s26s34s56x1111D2eps0()
      t59 = t28 * s45
      t60 = t59 + s25 - s12 - s16 + s34 + s56 - s13 + s46
      t61 = intHs16s25s26s34s56x1141D6eps0()
      t62 = s25 - s12 - s13
      t63 = -s13 + s25
      t64 = s16 - s34
      t65 = t28 * s26
      t66 = 3.0_dp
      t67 = t28 * (s26 + s45 - s13)
      t68 = t66 * s46
      t69 = t63 * s34
      t70 = s16 * t63
      t71 = s56 ** 2.0_dp
      t72 = s56 * t71
      t73 = t28 * t71
      t74 = intHs16s25s26s34s56x1311D4eps1()
      t75 = s12 ** 2.0_dp
      t76 = s12 * t75
      t77 = epinv * t13
      t78 = intHLs16s25s26s34s56x1113D4eps1()
      t79 = s34 + s45
      t80 = t28 * t79
      t81 = s16 + s26 - s56 - s46 - t80
      t82 = t41 * s12
      t83 = intHLs16s25s26s34s56x1112D2eps1()
      t84 = intHs16s25s26s34s56x1211D2eps1()
      t85 = t77 * t48
      t86 = intHLs16s25s26s34s56x1121D2eps0()
      t87 = s16 + s26 + s46
      t88 = -s45 - s56 - s46
      t89 = -s45 - s46
      t90 = t89 * s46
      t91 = s25 * t48
      t92 = s12 * t51
      t93 = s26 * s34
      t94 = s26 + s45 + s56
      t95 = t89 * s26
      t96 = t94 * s56
      t97 = t96 - t95
      t98 = t68 * s56
      t99 = intHs16s25s26s34s56x1221D4eps0()
      t100 = s25 - s16 - s26 + s34 - s56
      t101 = s16 + s26
      t102 = s34 + s45 - s56 - s46
      t103 = t28 * t101
      t104 = -t102 + t103
      t105 = s25 - s16 + s34 + s45
      t106 = s26 + s46 + s56
      t107 = t28 * t106
      t108 = s25 + s26 + s45 + s46 + s56
      t109 = s25 + s45 + s56
      t110 = s25 + s34
      t111 = s26 ** 2.0_dp
      t112 = s26 * t111
      t113 = s46 ** 2.0_dp
      t114 = s46 * t113
      t115 = t79 * s46
      t116 = t110 * s26
      t117 = s16 * t108
      t118 = t109 * s34
      t119 = s46 + s56
      t120 = s25 * t119
      t121 = s26 * t119
      t122 = s46 * s56
      t123 = -t89
      t124 = 5.0_dp
      t125 = t16 * s56
      t126 = t124 * s34
      t127 = t66 * s26
      t128 = t28 * t123
      t129 = t89 * s45
      t130 = s56 * t119
      t131 = s34 * t88
      t132 = t17 * t111
      t133 = t66 * s45
      t134 = -t124 * t131 - t28 * t130 + (t126 + t128 - t127 + s56 + s25
     #) * s25 + s26 * (-t123 * t17 - t125 + t126) - t129 - t132 - t133 *
     # s56
      t135 = s16 + s26 - s34 - s45 - s46 - s56
      t136 = s25 + s45
      t137 = s56 - s46
      t138 = t28 * t136
      t139 = s26 - s34 - s45
      t140 = t28 * t139
      t141 = -s26 - s45
      t142 = s34 * s56
      t143 = -s26 + s45
      t144 = s16 ** 2.0_dp
      t145 = s16 * t144
      t146 = s56 * t89
      t147 = t143 * s34
      t148 = s26 + s56
      t149 = s34 ** 2.0_dp
      t150 = s34 * t149
      t151 = s56 * t148
      t152 = t17 * s45
      t153 = t152 * s46
      t154 = -s26 + s34
      t155 = -s45 - s56
      t156 = -s34 - s56
      t157 = t155 * s26
      t158 = s34 * s45
      t159 = t156 * s45
      t160 = t154 * s56
      t161 = s26 * s45
      t162 = t144 * t108
      t163 = s45 * t71
      t164 = -s25 - s26 - s46
      t165 = s34 - s46
      t166 = -s34 - s46
      t167 = s45 ** 2.0_dp
      t168 = s25 * s26
      t169 = t168 * s46
      t170 = (-t149 - t113) * s45
      t171 = t149 * s56
      t172 = s45 * s56
      t173 = -t172 + t93
      t174 = s25 * s45
      t175 = t174 * s46
      t176 = t111 * s34
      t177 = t175 - t176
      t178 = s34 - s45
      t179 = s26 * s56
      t180 = -t113 - t71
      t181 = s25 - s45 + s56
      t182 = s25 ** 2.0_dp
      t183 = s45 * t137
      t184 = s25 * s46
      t185 = t135 * s13
      t186 = t111 * s46
      t187 = t180 * s56
      t188 = t12 * t13
      t189 = s12 - s34 - s45 + s56 + s46 + t103
      t190 = intHs16s25s26s34s56x1312D6eps0()
      t191 = t65 + t109
      t192 = s25 * s56
      t193 = s25 * t139
      t194 = t66 * s16
      t195 = t194 * s25
      t196 = s25 * (s26 - s34 - s45 - s46 - s56)
      t197 = s34 * s46
      t198 = t172 - t197
      t199 = t16 * s46
      t200 = t17 * s26
      t201 = t124 * s26
      t202 = t201 * s46
      t203 = (t200 + t199) * s56
      t204 = -s25 - s26 - s45
      t205 = s25 + s46
      t206 = t184 * s56
      t207 = t205 * s26
      t208 = -s16 - s26 + s34
      t209 = s25 * t113
      t210 = t141 * s46
      t211 = t136 * s26
      t212 = t211 * s46
      t213 = s25 - s46
      t214 = t213 * s46
      t215 = t214 + t192
      t216 = s25 * t149
      t217 = s34 - s56
      t218 = s26 * t101
      t219 = t217 * s46
      t220 = s16 * s56
      t221 = s34 * t113
      t222 = s26 * s46
      t223 = s26 + s46
      t224 = s56 * t223
      t225 = t16 * ((t224 + t113) * s56 + t186)
      t226 = 13.0_dp * s46
      t227 = s12 * s25
      t228 = ((t28 * t193 + (t191 - s13) * s13 + s46 * t109 - t111 + t11
     #3 + t192 - t179 + t195 + t227) * s12 + t124 * ((-t210 + t168) * s1
     #6 + (-s26 * t204 - t197) * s56 + t212) - t17 * ((s56 * t204 - t113
     # - t184) * s16 - t206 + t207 * s34) + t28 * (s25 * t144 + (s25 * t
     #208 - t142 + t71 - t93) * s45 + t71 * (s25 - s34 + s56) + t114 + t
     #209) - t66 * (((s16 + s56 + s46) * s34 - t111) * s25 + (t219 - t21
     #8) * s45 + (-t220 + t93) * s56 + t221) + 7.0_dp * s46 * (t220 + t2
     #22) - (t28 * t198 - t66 * t180 + s16 * (t46 - s25 + t68) - (-s26 -
     # s34 - s45 - s46) * s45 - t142 - t196 + t202 + t203 + t185) * s13 
     #- s45 * t215 + t111 * t64 - t167 * (-s25 - s16 + s46) + t112 + t21
     #6 + t225 + t226 * t179) * s12 + t188 * t106 * t104
      t229 = intHLs16s25s26s34s56x1213D6eps0()
      t230 = intHLs16s25s26s34s56x1221D4eps0()
      t231 = s16 + s26 + s46 + s56
      t232 = s25 + s26
      t233 = t124 * t232
      t234 = t17 * s56
      t235 = t28 * s46
      t236 = t148 * t17 + s16 + s25 + s34 + s45 + t235
      t237 = s34 + s45 + s46
      t238 = s16 * (t123 - t233 - t234)
      t239 = t237 * s46
      t240 = s26 * t232
      t241 = t232 * s34
      t242 = s26 - s46
      t243 = t242 * s56
      t244 = s45 * t148
      t245 = t16 * t244
      t246 = t124 * t241
      t247 = t66 * (t240 - t142)
      t248 = t152 * s25
      t249 = t28 * (t243 - t222)
      t250 = s25 - s45
      t251 = 7.0_dp * s26
      t252 = t154 * s34
      t253 = t119 * s46
      t254 = s26 + s34
      t255 = s26 * t148
      t256 = s45 * t223
      t257 = t16 * s45
      t258 = t257 * s56
      t259 = -t124 * t255 + t17 * t256 - t28 * (s46 * t254 + t144 - t158
     #) - s16 * (-t166 * t66 + t125 + t250 + t251) + t167 - t196 + t252 
     #+ t253 + t258
      t260 = 9.0_dp * s34
      t261 = t66 * t79
      t262 = t89 * s34
      t263 = s56 * t156
      t264 = t93 * t232
      t265 = s26 + s56
      t266 = t265 * s45
      t267 = s25 * s34
      t268 = s26 - s34 + s45
      t269 = t174 * t143
      t270 = t144 * (s25 + s26 - s45 - s46)
      t271 = s56 - s45
      t272 = t119 * s45
      t273 = t271 * s56
      t274 = s16 * (s25 * (-t65 + t119 + t260 + t152) - t16 * t146 - t17
     # * t263 - t28 * (t167 + t111 - t113) + t66 * (s26 * (t235 + s56 + 
     #t261) - t262))
      t275 = t88 * t149
      t276 = t93 * t119
      t277 = t179 * t141
      t278 = t205 * s56
      t279 = s45 * s46
      t280 = t279 * s56
      t281 = (t278 + t184) * s34
      t282 = s26 * t71
      t283 = 7.0_dp * t276
      t284 = t13 * t104
      t285 = intHLs16s25s26s34s56x1212D4eps0()
      t286 = t65 + s25 - s16 + s34 + s45 + s56
      t287 = s26 + s45 + s46 + s56
      t288 = 7.0_dp * s25
      t289 = -s34 + s45 + s56
      t290 = t156 * s26
      t291 = t155 * s34
      t292 = t289 * s46
      t293 = (t286 - s13) * s13
      t294 = s25 * t154
      t295 = t28 * t120
      t296 = t124 * t294
      t297 = t66 * t287
      t298 = t66 * t148
      t299 = t17 * t165
      t300 = -s26 - s45 - s46
      t301 = t300 * s45
      t302 = s46 * t165
      t303 = t124 * s46
      t304 = t17 * t217
      t305 = -t88
      t306 = 10.0_dp
      t307 = t16 * t119
      t308 = 7.0_dp * s45
      t309 = t306 * s34
      t310 = 12.0_dp * s26
      t311 = t3 * s56
      t312 = 9.0_dp * s46
      t313 = s45 - s46 - s56
      t314 = s25 + s56
      t315 = t313 * s25
      t316 = t315 * s34
      t317 = t93 * s45
      t318 = t279 * t217
      t319 = s25 - s34
      t320 = t319 * s26
      t265 = t265 * s34
      t321 = s25 * t167
      t322 = s45 * t110
      t323 = t165 * t167
      t324 = t149 * t119
      t325 = s26 * t149
      t326 = t188 * t231
      t327 = t28 * s25
      t328 = t327 * t76
      t175 = s12 * (s12 * (s16 * (t287 + t288) - t111 + t113 + t290 + t2
     #91 + t292 + t293 + t295 - t296 - t248) + s26 * (s56 * t226 - 9.0_d
     #p * t267) + t124 * (t111 * t314 + t316 + t161 * (-s25 + s56 + s46)
     #) + t16 * (s56 * (t224 + t320 + t113) + t222 * t232) + t17 * s46 *
     # (t192 - t197) - t28 * (s56 * t322 - t136 * t71 - t114 + t175 - t2
     #09 - t321 - t72) - t66 * (s34 * (-t267 + t71) + t111 * t178 + t317
     # + t318) - 7.0_dp * s46 * (-t222 + t265) + s13 * (t28 * (s34 * (s2
     #6 - s45 + s56) - t172) - t66 * (-t302 + t144 + t71) - t202 - t203 
     #+ s16 * (t59 + t299 - t298 + s25) - t149 + t196 + t301 - t185) + s
     #16 * (s16 * (t297 + t288) + s25 * (-t308 + t310 + t307 - t309) + s
     #26 * (-t17 * t178 + t311 + t312) - (t304 - t303 + s45) * t305 + t1
     #32) + t112 - t170 + t323 + t324 + t325) + t326 * t104 + t328
      t203 = intHs16s25s26s34s56x1212D4eps0()
      t329 = s26 - s34 + s46 + s56
      t330 = t66 * s34
      t331 = s34 * t119
      t332 = t136 * s34
      t333 = t250 * s56
      t334 = t250 * s46
      t335 = t28 * (s16 * t329 - t331) + s26 * (-t313 - t330 + s26) + t1
     #44 + t149 + t332 + t333 + t334 + t174
      t336 = -s25 - s34 - s56
      t337 = t336 * s26
      t338 = t217 * s56
      t339 = s34 * t223
      t340 = (t201 + t125) * s46
      t341 = s34 - s56 - s46
      t342 = 9.0_dp * s56
      t343 = t306 * s46
      t344 = t17 * t240
      t345 = t124 * s25
      t346 = s25 + s45
      t347 = t319 * s45
      t348 = (-s34 * t204 + t167) * s34
      t349 = t168 * s45
      t350 = t205 * s46
      t351 = s34 + s46 + s56
      t352 = s25 * t111
      t353 = t351 * t182
      t354 = 7.0_dp * s46
      t355 = 13.0_dp * s56
      t356 = intHs16s25s26s34s56x1222D6eps0()
      t357 = t3 * s25
      t358 = s25 + s45 + s46
      t359 = s13 ** 2.0_dp
      t360 = t358 * s46
      t361 = t16 * s34
      t362 = t107 + s16
      t363 = 9.0_dp * s25
      t364 = -t154
      t365 = 7.0_dp * t364
      t366 = s34 * t148
      t367 = 7.0_dp * s56
      t368 = -t16 * t302
      t369 = 17.0_dp * s26
      t370 = t306 * t119
      t371 = 13.0_dp * s34
      t372 = t17 * s34
      t373 = 16.0_dp * s56
      t374 = 18.0_dp * s46
      t375 = 13.0_dp * s45
      t376 = 11.0_dp * s45
      t377 = t124 * t111
      t378 = 9.0_dp * t279
      t379 = t168 * s34
      t380 = t113 * s56
      t381 = -s26 + s34
      t382 = s45 * t113
      t383 = t167 * s56
      t384 = s26 * t167
      t385 = s34 * t167
      t386 = t167 * s46
      t387 = t359 * t135
      t388 = s26 * t113
      t389 = s56 * t165
      t390 = s25 * t71
      t391 = 15.0_dp
      t392 = t306 * s26
      t393 = -t260 + t392
      t394 = s26 * (s45 * t393 + 30.0_dp * t122) - t124 * (t155 * t71 + 
     #t349) + t16 * t158 * (s25 - s46 - s56) + t17 * (t390 + t209 + t216
     # + t114 - t176) - t28 * (-t174 * t313 - t112 + t275 - t382 + t383 
     #+ t384 - t385 + t386 + t387) - t3 * (s25 * (t389 - t111 + t197) + 
     #t221) - 14.0_dp * s46 * (s26 * t381 - t71) + (s56 * (t179 + t161 -
     # t197) + t388) * t391 - 13.0_dp * t379 + 13.0_dp * t380 + 11.0_dp 
     #* s26 * (s46 * t346 + t192) - 12.0_dp * t179 * t154 + 7.0_dp * s56
     # * (t279 - t142) - s13 * (s56 * (t133 + t226 + t367) + t124 * (-t3
     #66 + t144) - t28 * (t196 + t129 - t149 - t111) + 11.0_dp * t121 + 
     #s16 * (t119 * t3 - t133 - t327 + t365) + t161 + t368 + t152 * s34)
     # + s16 * (s16 * (t124 * t305 + t127 + t363) + s56 * (20.0_dp * s46
     # + t376) - t180 * t306 + 7.0_dp * t131 + s25 * (-t308 - t371 + t36
     #9 + t370) + s26 * (t373 + t375 + t374 - t372) - t167 + t377 + t378
     #) + t325
      t395 = (-t148 + t327) * s12
      t278 = s12 * (s12 * (s25 * (t201 - t361) - t17 * (t255 + t174) + t
     #244 * t66 + t28 * (t262 + t278 - t359 + t360) + s13 * (t105 * t28 
     #+ s56 + t127) + s16 * (t128 - s26 + t357) - t142 + t395) + t394) +
     # t188 * t104 * t362
      t394 = intHs16s25s26s34s56x1212D4eps1()
      t286 = -t286
      t396 = (-7.0_dp + epinv) * s25
      t397 = t2 * t313
      t398 = t66 * epinv
      t399 = t398 * s34
      t400 = t2 * s26
      t401 = t341 * epinv
      t402 = -t111 - t144
      t403 = t402 * t2
      t404 = epinv * s34
      t405 = t404 * t119
      t406 = -t28 * t405 + epinv * (t334 + t333 + t332 + t149 + t174) + 
     #s16 * (-t28 * (t401 + t400) + t102) + s26 * (t397 + s34 - t399) + 
     #t272 + t331 + t403
      t407 = t2 * s25
      t408 = epinv * s45
      t409 = epinv * t165
      t410 = t400 + s56
      t411 = epinv * s56
      t412 = t28 * epinv
      t413 = -1.0_dp + t412
      t414 = t2 * s46
      t415 = t167 * t2
      t416 = t222 * t2
      t417 = t179 * t2
      t418 = t404 * t223
      t419 = t172 * t2
      t420 = t122 * t2
      t421 = t158 * epinv
      t422 = t142 * epinv
      t423 = t2 * t144
      t424 = t71 * t2
      t425 = t113 * t2
      t426 = -t16 * t420
      t427 = -t124 * t416
      t428 = -t66 + epinv
      t429 = epinv * t119
      t430 = t124 * epinv
      t431 = t17 - t430
      t432 = t431 * s34
      t433 = s46 - t411
      t434 = t199 * epinv
      t435 = epinv * t306
      t436 = t435 * s46
      t437 = t132 * t2
      t438 = -7.0_dp + t398
      t439 = t297 * t2
      t440 = epinv * s26
      t441 = -s34 - t440
      t442 = t2 * t111
      t443 = t161 * t2
      t444 = t441 * s25
      t445 = t174 * s34
      t446 = t2 * s45
      t447 = epinv * s25
      t448 = -t447 + t446
      t449 = -t446 + s34
      t450 = t447 * t113
      t451 = t414 + t400
      t452 = t400 + s25
      t453 = t184 * epinv
      t454 = s56 * t451
      t455 = -1.0_dp + epinv
      t456 = -t446 - s25
      t457 = epinv + 1.0_dp
      t458 = t158 * t2
      t459 = epinv * t149
      t460 = t174 * t457
      t461 = t2 * t72
      t462 = t2 * t114
      t463 = epinv * t182
      t464 = t2 * s13
      t465 = t464 * t135
      t466 = t112 * t2
      t431 = s13 * (-t17 * (-t418 - t417) - t28 * (s34 * (t404 + s26 - s
     #45 + s56) - t419) + t66 * (t423 + t422 - t421 + t425 + t424 - t197
     #) + s16 * (t124 * t409 - t17 * (s34 - s46 + t411) + t410 * t66 - t
     #407 + t408 - t59) - s25 * (t413 * t48 + t400) - s45 * (-t414 - t40
     #0) + t149 + t415 - t426 - t427 + t465) + s16 * (s16 * (s25 * t438 
     #- t439) + s25 * (-t124 * (s34 * t29 - t429) - t307 + t308 + t200 *
     # t428) + s26 * (-s45 * t431 - t311 + t432 - 9.0_dp * t433 + t436) 
     #+ t305 * (-t124 * t433 + s45 - t234 + t432 + t434) - t437) + s45 *
     # (-t425 + t463 - t149) + t167 * (t447 + t414 - s34) - t324 - t325 
     #- t466 + t447 * t111 + t353 * epinv
      t432 = 13.0_dp * t417 * s46
      t265 = s12 * (-s12 * (-epinv * (-t118 + t337 + t292 - t111 + t113)
     # - t2 * (s13 * (t286 + s13) - s16 * t287) - t111 + t113 + t290 + t
     #291 + t292 - t396 * s16 + t295 - t296 - t248) - t124 * (s26 * (s26
     # * (t404 + s25) - t174) + s46 * (t444 + t443) + s56 * (t442 + t444
     # + t443) + t445 + t404 * t113) - t16 * (s56 * (t454 + t425 + t320 
     #- t453) + t222 * t452) - t17 * (t206 - t221 + t404 * (t211 + t71))
     # + t28 * (epinv * (t349 + t348) + s46 * (t460 + t459) + s56 * (t46
     #0 + t459 + t458) + t456 * t71 - t209 - t321 - t461 - t462) + t66 *
     # (s34 * (-t271 * t447 - t267 + t71) + s46 * (s34 * t448 - t419) + 
     #t111 * t449 + t317 + t450 + t447 * t71) + 9.0_dp * s34 * (t411 * (
     #-s26 - s46) + t168) + 7.0_dp * s46 * (t222 * t455 + t265) + t431 -
     # t432 - t436 * t93) + t284 * t406 - t328
      t295 = intHLs16s25s26s34s56x1321D6eps0()
      t431 = t46 + s25 + s34 + s45 + s46 + t127
      t436 = s16 * t51
      t444 = s25 - s45 - s46
      t460 = t444 + t127 + t372
      t467 = t158 + t113 + t71
      t468 = t28 * t467
      t469 = 7.0_dp * s34
      t470 = t142 * t124
      t471 = t148 * t167
      t472 = s46 * t71
      t473 = -t122 - t149
      t474 = t144 * t51
      t475 = s56 * t141
      t476 = t179 * s45
      t477 = t271 * t111
      t478 = s16 * (-t17 * (-t292 + t158) + t28 * (t167 + t111 + t113 + 
     #t71) + s25 * (t107 - s45 - t469) + s26 * (-t59 + s46 + t234 - t469
     #) - t470)
      t479 = -t279 - t71
      t480 = t17 * (s34 * t479 - t155 * t222 + t280)
      t481 = t306 * t276
      t482 = t65 - t102
      t483 = t28 * s34
      t484 = -t305 + t483
      t485 = 9.0_dp * t119
      t486 = t28 * (s26 * (t305 - t126) + t122 + t149) + s34 * (-t485 + 
     #s45) + t113 - t167 + t71
      t487 = t330 * t106
      t488 = t487 * t482
      t489 = intHs16s25s26s34s56x1213D6eps0()
      t490 = t17 * s25
      t491 = t155 * s46
      t492 = s25 * t143
      t493 = t267 * t66
      t494 = t17 * (t179 - t197)
      t495 = t124 * s56
      t496 = t3 * s26
      t497 = t124 * s45
      t498 = t267 - t113
      t499 = t111 + t158
      t500 = t111 * s45
      t501 = s34 * t71
      t502 = -t350 - t192
      t503 = t111 * s56
      t504 = -s25 - s34 + s46
      t505 = t3 * s34
      t506 = -t143 - t483
      t507 = t16 * s26
      t508 = t124 * t119
      t509 = t28 * t119
      t510 = t506 * t106
      t511 = t510 * t482
      t512 = t28 * t145
      t513 = s16 * (s16 * (t508 - t330 - s45 + t507) + s26 * (-t17 * (t2
     #8 * t341 + s45) + t507) + s34 * (s45 - 7.0_dp * t119 + s34) - t509
     # * (t59 - t119)) + t511 + t512
      t214 = s12 * (s12 * (-t28 * t492 + s16 * (t287 + t490) - s34 * t28
     #7 - t111 + t113 + t120 - t179 + t293 - t491 - t493 + t227) - 7.0_d
     #p * s26 * t498 + t124 * (s26 * (t272 + t120) + s34 * t502 + t503) 
     #- t17 * (s46 * (-t192 + t158) + t176 + t317 + t501) + t28 * (s45 *
     # (-t168 + t71) + t314 * t71 + t114 + t209 + t216) + t66 * (s25 * t
     #499 - t165 * t172 + t500) + t225 - 9.0_dp * t197 * t148 - s13 * (-
     #t28 * t159 + t66 * (-t338 - t93 + t144 + t113) + t340 - s16 * (t30
     #4 + t136 - t303 - t127) + t149 - t196 - t301 + t494 + t185) - s16 
     #* (-s16 * (t345 + t297) - s25 * (-t133 + t496 + t307 - t469) - s26
     # * (t342 + t343 + t497 - t372) + t305 * (-t495 + t372 - t199) - t1
     #32) - s45 * (t214 - t149 + t192) - t167 * t504 + t112 + t324 + t32
     #5 + t179 * (t226 - t505)) + t13 * t513
      t225 = intHLs16s25s26s34s56x1312D6eps0()
      t297 = intHLs16s25s26s34s56x1231D6eps0()
      t304 = t109 + t127
      t317 = s46 * t136
      t514 = t319 * s56
      t515 = t66 * (-t244 + t436)
      t516 = t28 * t305
      t517 = -t516 + s34
      t518 = s34 * t178
      t519 = t93 * t3
      t520 = t28 * t517
      t521 = s26 - s56
      t522 = t66 * t521
      t523 = t17 * s46
      t524 = -t59 + t110 - t523 + t522
      t525 = t79 * s45
      t526 = t178 * s56
      t527 = t3 * s45
      t528 = t306 * s56
      t529 = t66 * t167
      t530 = t154 * s45
      t531 = t154 * t113
      t532 = t142 * s45
      t533 = -t149 - t167
      t534 = s16 * (s46 * (t527 + t528) - t124 * t180 + t28 * (-t526 + t
     #111) + s25 * (t65 + t508 - t372) - s26 * (t299 + s45 - t367) + t26
     #2 + t529)
      t276 = t3 * t276
      t535 = t16 * t122 * t119
      t390 = s12 * (s12 * (t92 - t28 * (t241 - t224 + t174) + s13 * (t30
     #4 - s13) + t113 + t240 + t71 + t317 + t514 + t515) - t124 * (s46 *
     # (t530 + t267) + t264 - t282 + t531 + t267 * s56) + t17 * s56 * (-
     #t142 + t184) + t28 * (t113 * t136 + t114 + t186 + t390 + t445 + t4
     #74 - t532 + t72) + t66 * (t119 * t168 + t280 + t471 + t503) - s13 
     #* (-s16 * t524 - s46 * (-t125 + t372) + t28 * (t179 + t279) - t66 
     #* (t526 + t95 + t111 - t113 - t71) - t196 + t525 - t93 + t185) - s
     #25 * t533 + t163 - t500 + t534 + t93 * t178 - t276 + t535 + t122 *
     # t393) - t13 * (s16 * (-t17 * (s26 * t305 + t122) - t28 * (-t167 +
     # t113 + t71) - t518 + 9.0_dp * t331 + t519 + t520 * s16) + t488)
      t526 = intHLs16s25s26s34s56x1222D6eps0()
      t536 = t148 + t327
      t537 = t16 * s25
      t538 = t66 * s56
      t539 = t138 + t201 + t538
      t540 = s16 * (t46 + t127 + t537)
      t541 = s25 + s45 + s46 + s56
      t542 = s25 * t79
      t543 = t66 * (t244 - t168)
      t544 = t3 * s46
      t545 = 9.0_dp * s26
      t546 = 11.0_dp * s56
      t547 = 12.0_dp * s46
      t548 = s26 * (t133 - t372 + t528 + t547)
      t549 = 16.0_dp * s46
      t550 = (t343 + t375) * s46
      t308 = (t308 + t549 + t125) * s56
      t551 = s25 - s56
      t552 = s45 - s46
      t553 = t552 * s46
      t554 = -t122 + t93
      t555 = -t222 + t142
      t556 = 22.0_dp
      t557 = 11.0_dp * s46
      t558 = t306 * t111
      t559 = t122 * t3
      t511 = t28 * (t511 + t145) + s16 * (s16 * (t485 - t330 + s45 + t49
     #6) - t17 * t180 + s26 * (t119 * t391 - t371 - t497) + t149 - t167 
     #- 7.0_dp * t272 + t558 + t559 - t371 * t119)
      t560 = t536 * s12
      t561 = s13 * t539 - t17 * t542 + t28 * (s46 * t541 + t192 - t359 -
     # t93) - t142 + t540 - t543 + t560
      t96 = s12 * (s12 * t561 + t124 * s56 * ((s26 + s46) * s45 - t142) 
     #+ t17 * (s25 * (t111 + t113 + t71) + t114 + t158 * t551) - t28 * (
     #s25 * (-s45 * t313 - t149) + s45 * (t553 - t172) - t112 - t384 + t
     #387 - t500) - t3 * (s46 * (t332 - t192 + t197) + s56 * (t267 - t11
     #1)) - t306 * (s56 * t554 - t186) - t66 * (t110 * t161 - t163 - t72
     #) - 9.0_dp * s26 * (s46 * (-s25 - s45) - t192 + t267 - t71) - 13.0
     #_dp * s46 * t555 - s13 * (s46 * (t545 + t546) + t124 * t96 - t28 *
     # (-s45 * t237 + t111 + t196) - t66 * (t142 - t161) + t368 - s16 * 
     #(s26 + s34 - s45 - t544 - t234 + t327 - s16) - t93) + s16 * (s16 *
     # (t305 + t345 + t127) + t550 + s25 * (t545 - t133 + t370 - t469) +
     # t131 + t548 + t377 + t529 + t308) + t325 - t93 * (t547 + t507) + 
     #t122 * (s26 * t556 + t557)) + t13 * t511
      t368 = -s26 + s34 - s46
      t370 = s16 * t305
      t387 = t154 * s46
      t561 = t28 * t244
      t562 = t49 * t27
      t563 = t28 * s16
      t564 = s12 + t563
      t565 = s16 + s25
      t566 = s12 - s13
      t567 = t28 * t565
      t568 = -s13 + s25 + s26
      t569 = t28 * t232
      t570 = t568 * t101
      t571 = -s16 + s25 + s34
      t572 = s25 + s26 + s34 - s13
      t573 = -t208
      t574 = t572 * t573
      t575 = -s12 - s13 + s25 + s26 + s34
      t576 = t59 + t575
      t577 = t66 * s12
      t578 = t577 + s13
      t579 = t28 * s12
      t580 = -s16 + t579
      t581 = -s13 + s16 + s25
      t582 = t124 * s12
      t583 = s25 - s34 - s13
      t584 = s16 + s45
      t585 = s12 * t566
      t586 = t578 * s16
      t587 = (-t28 * t584 + s12 + s13 - s25 + s34 - s46 - t127) * s46
      t588 = -t28 * t300
      t589 = -t59 - s25 + s16 - s34 - s56 + s13 - s46
      t590 = s12 * (t589 + s12) + t49
      t591 = t564 * s25
      t592 = s13 * s16
      t593 = t568 * s34
      t594 = t573 * s46
      t595 = t594 * t576
      t596 = t573 * t167
      t597 = t573 * t113
      t598 = propW16 * t51
      t599 = t27 * t590 * t598 - s12 * t111
      t600 = t194 + t579
      t601 = t66 * t565
      t602 = t66 * t232
      t603 = s16 - s26 - s34
      t604 = t603 - t327
      t605 = t17 * s12
      t606 = 7.0_dp * s12
      t607 = t101 + t605
      t608 = s12 * s13
      t609 = t28 * t75
      t610 = s12 - s45 - s46
      t611 = t21 * (-t28 * t599 - s12 * (-s26 * (t566 + t567) - t591 + t
     #592) - s34 * (-s12 * (-t569 + s13) + t570 - t593) - (s12 * t571 + 
     #t574) * s45 - s56 * (s56 * (t581 + t588 - t577 + s56) - s25 * t580
     # - s26 * (-t581 + t582 - s26) - s34 * (-t101 - t577 + s34) - s45 *
     # (-t583 - t127 - t563 - s45) - t585 - t586 - t587) - t595 - t596 -
     # t597)
      t612 = intHs16s25s26s34s56x1221D4eps1()
      t613 = t2 * s34
      t614 = t2 * s16
      t615 = s56 * t457
      t616 = epinv * s46
      t617 = t2 * s56
      t618 = -t617 + t616
      t619 = t457 * s46
      t620 = t618 * s26
      t621 = t457 * s45
      t622 = t457 * s34
      t623 = epinv * t111
      t624 = t267 * epinv
      t625 = t457 * t113
      t626 = epinv * t136
      t627 = t398 - t28
      t628 = epinv * t79
      t629 = t627 * s16
      t630 = epinv * t113
      t631 = s26 - s45 + s46
      t632 = t179 * epinv
      t633 = t631 * s34
      t634 = epinv * t71
      t635 = t143 * epinv
      t636 = t172 * epinv
      t637 = t122 * epinv
      t638 = t222 + t71
      t639 = s56 * t155
      t640 = t29 * s25
      t641 = epinv * t112
      t642 = t623 * s46
      t643 = t168 * epinv
      t644 = epinv * t72
      t131 = s26 * (-t463 + t167) + t275 + s13 * (t124 * t638 - t28 * (t
     #632 + t634 + t633 + t459 - t279) + t3 * t224 - t66 * (s34 * (t635 
     #+ s56) - t113 + t636 + t637) - s16 * (t629 - t28 * t626 - t66 * (t
     #400 - s34 + s46) - t125 + epinv * t137 + s25 + s45 - t404 * t124) 
     #+ s25 * (s26 * t413 - t28 * t628 + s34 + s45 + t414 + t617) - s45 
     #* (-t400 - s56) + t111 + t149 + t415 + t422 - t630 - t153 * epinv 
     #+ t465) + s16 * (-s16 * (-t66 * (epinv * t287 + t640) + t516) - t1
     #31 * t66 + t17 * t90 + t3 * t157 - epinv * t134 - s25 * (t508 - t2
     #60 + t310 - t257) + 7.0_dp * t639 - 7.0_dp * t222 - 11.0_dp * t122
     #) + s45 * (s25 * (t447 + s56 + s46) + s45 * (t447 + s56)) + t641 +
     # t642 + t404 * (s46 * (s25 + s45 - s46) + t182) - t643 * t271 + t4
     #11 * t113 + t644 - t158 * (s45 + t411)
      t463 = t111 - t197
      t645 = t184 * s34
      t646 = -t414 + s34
      t647 = t411 * t149
      t648 = t161 * epinv
      t649 = 18.0_dp * s26
      t395 = s25 * (-t126 - t152 + t507) + t28 * (s25 * (t619 + t615) + 
     #t620 + t122 * t457) - t66 * (-t244 + t111) + s13 * (-epinv * t106 
     #* t28 + s26 - s56 + t407 + t446 - t464 + t613 - t614) - s16 * (s26
     # * t35 + t396 - t414 - t446 + t615) + s26 * (t447 + t622) + s46 * 
     #(t621 - t613) + t457 * t71 + t422 - t458 + t623 + t624 + t625 + t3
     #95
      t396 = t77 * t12
      t650 = intHLs16s25s26s34s56x1321D6eps1()
      t651 = epinv * t51
      t652 = -s26 + s34 - s56
      t653 = t652 * s46
      t654 = t404 * (s26 + s46 + s56)
      t655 = t139 * s34
      t656 = -t28 * t655
      t657 = t66 * t654
      t658 = -t484 * t614 + t151 + t244 - t653 + t656 + t657 + t91
      t659 = t351 * s56
      t660 = t2 * t359
      t661 = s16 * s34
      t662 = t93 * t17
      t663 = -t404 - s56
      t664 = t279 * t2
      t665 = t178 * epinv
      t666 = s26 * t165
      t667 = t167 + t113
      t668 = 7.0_dp * t111
      t669 = -t93 - t113
      t670 = -s16 + s45
      t671 = -s26 - s56
      t672 = t232 * s56
      t673 = s16 * s46
      t674 = t89 * t71
      t675 = t64 * s46
      t676 = s16 * s26
      t467 = -epinv * (s34 * (s46 * (-t303 - t342) - t668) - t28 * (t113
     # * t141 + t141 * t71 + t148 * t158 - t325 - t474) + t66 * (t472 + 
     #t380 + t471) + t480 - t481 + s25 * (s34 * (-t508 + t133) - t28 * t
     #473 + s26 * (t119 - t469) + t113 + t167 + t172 + t279 + t71) + s46
     # * t667 - t275 + t385 + t476 + t477 + t478 + t72) + s13 * (t124 * 
     #t146 + t17 * (s45 * (-s26 - s34 + t411) - t616 * t217) + t28 * (ep
     #inv * t467 + t302 + t666) + t66 * (-s26 * t665 + s56 * t663 - t149
     # + t442 - t664) + s26 * (t616 + t407 - s56) - t415 + t459 + t640 *
     # t48 + t614 * t460 - t465) + s16 * (-t158 + t71) + s25 * (t675 + t
     #220 + t167) - t149 * t178 + t182 * t48 + t114 - t503 + t279 * t79 
     #- t676 * t271
      t207 = s12 * (s12 * (s12 * (t651 + s34) + t28 * (-t149 + t111 + t2
     #07 + t197) + t66 * (t179 + t661) + epinv * (t66 * (-t244 + t436 - 
     #t93) + s25 * (-t59 - t330 + s26) + s26 * t137 + t111 + t262 - t142
     # * t28) + s25 * t351 - s45 * t652 + t659 + t660 - t464 * t431 + t6
     #62) - t124 * (s34 * (-t111 - t673) + t349 - t674 - t192 * t89) + t
     #16 * s26 * t473 - t17 * (s46 * (s45 * t565 + t122 + t149) + t267 *
     # t670) + t28 * (s16 * (-t174 - t167 + t111 + t168 - t113) + s25 * 
     #(s26 * (s26 + s34 - s46) - t113) + s34 * t144 + s45 * t669 + s56 *
     # (-t158 - t168) + t112 - t385 - t72) - t3 * s26 * (t172 - t661) + 
     #t66 * (s34 * (s56 * t87 + t113) + s45 * (-t220 - t111) + s56 * (-t
     #672 - t673) + t149 * (-s25 - s16 - s56) - t388 - t471) + 7.0_dp * 
     #s46 * (s45 * t671 + t93) - t467) - t284 * t658
      t467 = intHLs16s25s26s34s56x1221D4eps1()
      t503 = -s26 - s34 - s56
      t677 = t503 * s45
      t678 = -s25 - s26 - s45 - s46
      t679 = -s26 + s45 + s46
      t680 = t679 * s46
      t681 = t678 * s34
      t682 = t28 * t51 * s26
      t683 = t124 * (-t95 + t172 + t71)
      t684 = t17 * t119
      t685 = -s34 + s45 + s46
      t686 = t178 * s46
      t687 = t149 * s46
      t688 = t199 * s56
      t689 = t136 * s56
      t690 = s56 * (t689 + t168 - t197) + t169
      t691 = -t71 + t197
      t692 = s16 + s56
      t693 = s26 * t691 - t279 * t692 - t380
      t694 = t161 * t119
      t695 = t136 * s45
      t696 = s16 * (-t695 - t222) + s34 * (t255 - t167) + t149 * t678 - 
     #t186
      t697 = -s16 * t71 + s25 * (s16 * t341 - t113) + s45 * (t661 - t111
     #) - t471 - t72 + t661 * t223
      t698 = s25 - s45 - s56
      t699 = (s25 + s16 - s56) * t111
      t220 = s25 * (s25 * t88 - t167) + s34 * (t317 + t220 - t182) + t11
     #3 * t89 + t112 - t171 + t699 + t676 * t698
      t700 = t391 * s26
      t701 = 9.0_dp * s16
      t702 = t122 * (t700 + t701)
      t703 = 7.0_dp * t148 * t174 + 7.0_dp * t388
      t704 = t306 * (t122 * t314 + t694)
      t211 = t124 * (s16 * (t172 + t113) - t221 - t174 * t166) + t16 * t
     #690 - t17 * t697 - t28 * t220 - t3 * t693 - t66 * t696 + epinv * (
     #s26 * (t519 + t688) - t124 * s45 * (t266 - t111) + t17 * (s56 * (t
     #686 + t179) + t186) - t28 * (s26 * (t553 + t158) + t552 * t71 + t2
     #70 + t279 * t685) + t66 * (-t277 + t221 - t325) + t283 + s25 * (s3
     #4 * (-t330 - t497 + t684) + s26 * (t59 + t119 + t505) - s45 * (t59
     # + t119)) + s56 * (-t149 - t263) + t274 + t380 - t385 - t687 + t15
     #8 * t156) - s13 * (t17 * (s34 * t552 + t179) - t28 * (t93 - t279) 
     #- t66 * (s25 * t155 + t111 - t113 - t149 - t184 - t267) + t559 - e
     #pinv * t259 - s16 * (-t59 + t299 + s25 + t522) - t142 + t167 - t16
     #8 + t683 + t465) - s34 * (t211 - t172 + t71) + t702 + t703 + t704
      t220 = intHLs16s25s26s34s56x1212D4eps1()
      t263 = t428 * s34
      t299 = t65 * t2
      t519 = -t617 - t414
      t522 = s34 * t79
      t553 = t28 * (t91 + t522) + s16 * (-t397 + t263 + t299) + s34 * (s
     #26 * t428 - t414 - t617) + s45 * t519 + t423 + t442 - t400 * t313
      t559 = t447 * t119
      t286 = t464 * t286
      t690 = t17 * epinv
      t693 = -t16 + t690
      t696 = t693 * s34
      t697 = -1.0_dp + t398
      t705 = 7.0_dp * epinv
      t706 = -t705 + t66
      t707 = t617 + t400
      t708 = t404 + t414
      t428 = s25 * (t428 * t48 + t400)
      t709 = t197 * t2
      t710 = -t616 - s34
      t711 = t404 + t407
      t712 = t93 * (-t447 + s46)
      t713 = 1.0_dp - epinv
      t714 = t455 * s46
      t715 = s34 + s45 - s46
      t323 = epinv * (-t170 + t323 + t324 + t325 + t112) + s13 * (-t17 *
     # (-t147 - t417) - t28 * (s56 * (t613 - t446) + t404 * t143) + t66 
     #* (t423 + t425 + t424 + t149 - t709) - t426 - t427 - s16 * (-t17 *
     # t708 + t28 * t446 - t66 * t707 + t361 + t407) + t415 + t443 - t45
     #9 + t664 - t428 + t465) + s16 * (-s16 * (s25 * t706 + t439) + s25 
     #* (s45 * (1.0_dp - t705) + t16 * (t429 + s34) + t17 * (s26 * t697 
     #- s46 - s56) - t404 * t306) + s26 * (-t693 * s34 - t2 * (t312 + t3
     #11 + t152)) - t305 * (-t2 * (-t303 + s45 - t234) + t696) - t437) +
     # s25 * (t331 - t111) - t112 + t279 * t715
      t437 = t447 - t446
      t439 = s25 + t408
      t323 = t124 * (s25 * (-t623 + t458 + t648) + s46 * (t624 + t443) +
     # s56 * (s26 * (t400 + t446) + t624) - t176) + t16 * (s56 * (s26 * 
     #t711 + t425 + t454) + t222 * (t407 + t400)) + t17 * (s25 * (s56 * 
     #(-t616 + s45 + s56) - t90) - t613 * t113) - t28 * (s25 * (s45 * t4
     #56 - t267) + s56 * (-epinv * t322 - t182) + t437 * t71 + t450 - t4
     #61 - t462 - t184 * t439) - t3 * s56 * (t93 - t184) - t66 * (s26 * 
     #(s34 * (t446 + s25 - s34) - t174) + s34 * (s45 * t710 - t2 * t267 
     #- t167 - t197) + s56 * (-t149 - t664) + t111 * (-t404 - t446) + t6
     #13 * t71) + t432 - 9.0_dp * t712 - 7.0_dp * s46 * (s26 * (-t404 + 
     #t714) + t142 * t713) - t323
      t456 = intHLs16s25s26s34s56x1312D6eps1()
      t716 = t693 / 2.0_dp
      t717 = t319 * s46
      t718 = s46 + s56
      t719 = t718 * epinv
      t397 = t28 * s34 * (t719 + s34 + s45) + s16 * (t28 * (-t429 + t400
     #) - s45 + s46 + s56 + t263 + t614) + s26 * (s34 * t716 - t397) + s
     #45 * (-t617 - t414 + s25) + t267 + t442 + t514 + t717
      t720 = -1.0_dp + t690
      t721 = -t17 + epinv
      t722 = t88 * epinv
      t723 = t447 * t143
      t724 = -t407 + t408
      t725 = t217 * epinv
      t726 = t17 * t413
      t727 = t133 * epinv
      t728 = t430 * s45
      t729 = 11.0_dp * epinv
      t730 = -s34 - t411
      t731 = -s34 - t408
      t732 = s25 - t408
      t733 = t411 * s25
      t734 = t404 * t111
      t735 = t617 + t414
      t736 = t350 + t192
      t737 = t617 * t111
      t738 = t447 * t149
      t739 = t142 * s46
      t740 = s34 - t408
      t741 = t448 * t111
      t454 = (s46 * (t414 + s25) + t454) * s56
      t438 = s16 * (-t66 * t722 + t126 - t684 + (-t66 + t430) * s25 + s2
     #6 * t438) - t16 * t616 * t89 - s56 * (s46 * (9.0_dp - t729) + t59)
     # + t124 * (-t155 * t411 - t113) - t17 * (s34 * (s34 - t722) + t71)
     # + t3 * t331 + s25 * (t16 * t429 - t684 + s26 * t726 + (-t705 + t1
     #6) * s34 + s45 - t727) - s26 * (-s26 * t693 + (-13.0_dp + t690) * 
     #s34 + t2 * (t342 + t343) + s45 - t728) + t158 - t68 * s45
      t693 = t504 * epinv
      t742 = t149 * s45
      t410 = -t119 * t182 + t150 + s13 * (-t17 * (s34 * (-t616 + s26 - s
     #45) - t417) - t28 * (s56 * t449 + t421) + t66 * (-s34 * t441 - t19
     #7 + t422 + t423 + t424 + t425) - t426 - t427 - s16 * (-t17 * (t725
     # + s46) - t410 * t66 + t361 + t59 - t724 + t616 * t124) - s25 * (t
     #29 * t48 + t400) + t415 + t443 - t459 + t664 + t465) + s16 * t438 
     #+ s34 * (t405 - t182) + s45 * (s34 * (t404 - s56) - t182 + t630) -
     # t352 - t466 - t693 * t167 + t459 * s26 - t742 + t279 * t552 + t17
     #4 * (-t429 - s45)
      t438 = t411 * t3
      t410 = s12 * (-s12 * (t447 - s16 + s34) - t28 * (-s34 * t341 - t72
     #3) + t66 * (t624 + t144) - t662 - s16 * (s25 * t720 + s26 * t721 +
     # s45 + t126 - t509 - t722) - s26 * (t617 - t404 + s25) - s46 * (ep
     #inv * t319 - t446 - t617) - t267 - t286 + t425 - t442 - t660 + t15
     #8 * t457 - t411 * t319) + t124 * (s26 * (s25 * t735 + s45 * t735) 
     #+ t737 + t404 * t736) + t16 * (s26 * (t149 + t416) + t454) - t17 *
     # (s34 * (s56 * t730 - t174) + s46 * (s34 * t731 + t733) + t221 - t
     #734 + t93 * t732) + t28 * (-t178 * t184 - t178 * t192 + t349 * t45
     #7 + t448 * t71 + t145 + t385 - t450 + t461 + t462 - t738) + t66 * 
     #(s25 * (s34 * t740 + t113) + t319 * t71 + t741 + t172 * t708) + t4
     #32 - 9.0_dp * s34 * (t616 * t671 + t179) - 7.0_dp * s26 * (-t624 -
     # t425 + t93) - 7.0_dp * t739 - t410 + t93 * (-t343 + t438)
      t443 = intHLs16s25s26s34s56x1222D6eps1()
      t448 = t66 * t119
      t671 = t17 * t79
      t735 = -t671 - t448 + s26
      t743 = t66 * (s34 - s46 + t440)
      t744 = t28 * t148
      t745 = t143 * s46
      t746 = s25 * t735
      t747 = t341 * s56
      t748 = (t17 * t429 + s25 - t614 - t665 + t743 - t744) * s16
      t749 = t387 + t93
      t750 = t66 * t749
      t751 = t28 * (epinv * (s26 * t143 + s56 * t143 + t745) + t149 + t1
     #58 + t179 - t279)
      t654 = t17 * t654
      t752 = t179 * t17
      t753 = t29 * s26
      t754 = s26 + t408
      t755 = t754 * s26
      t756 = t93 * epinv
      t757 = 9.0_dp * epinv
      t758 = (-t306 + t757) * s26
      t759 = t407 + s45
      t760 = epinv * t167
      t761 = t446 * t113
      t762 = t404 * t179
      t763 = t205 * s45
      t764 = t168 * t2
      t765 = t400 * t71
      t766 = t616 - s56
      t308 = s13 * (s45 * (-t372 - t367) + s46 * (s56 * (-12.0_dp + t729
     #) + t758) + t124 * (s56 * (-t400 + t408 + t411) - t161 + t197) - t
     #16 * (s46 * t708 + t71) + t28 * (t421 - t623 - t149 - t664 - t415)
     # + t66 * (-t422 + t755) + s16 * (-t716 * s25 + s46 * (epinv * t3 -
     # 7.0_dp) + s56 * (-t28 + t690) + t330 - t665 - t753 - t614) + s25 
     #* (-t716 * s26 + t28 * (t263 - t722) - t257 - t508) - t756) + s16 
     #* (s16 * (-epinv * (t345 + t127) + s25 + s26 + t414 + t446 - t615)
     # + s45 * (t557 + t200 + t367) + t28 * (t179 + t167) + t66 * (s34 *
     # t300 + t71) - epinv * (-s34 * t305 + t308 + t377 + t529 + t548 + 
     #t550) + s25 * (-t306 * t429 + t199 + t495 + t727 - (1.0_dp + t757)
     # * s26 - t706 * s34 + s45 - s25) + 9.0_dp * s46 * t223 + 12.0_dp *
     # t122) + s26 * (-t459 - t182) - t112 - t501 - t645 + t142 * (s25 -
     # s45 + s34)
      t508 = t172 * t223
      t548 = t279 * s34
      t550 = t110 * epinv
      t706 = t408 * t71
      t727 = t182 * s46
      t767 = t718 * s56
      t768 = epinv * t556
      t769 = s56 * (s25 * t391 + s26 * (23.0_dp - t768) - 11.0_dp * t616
     #) + 14.0_dp * t222
      t135 = s46 * t769 + t124 * (epinv * (t501 - t508) - t548) + t16 * 
     #(s45 * (t267 + t111) + t734) - t17 * (s25 * (t623 - t267) + s45 * 
     #(t404 * t551 - t182) + s56 * (s56 * (t447 - s56) - t111) - t321 + 
     #t450 - t462) + t28 * (s26 * (s34 * (s25 + s34 - s45 - s56) - t760)
     # + s46 * (s45 * t437 + t149) + t111 * (-s25 - s34 - t408) + t149 *
     # t759 + t385 - t641 + t761 + t636 * t250 - t660 * t135 - t447 * t1
     #67) - t3 * (s34 * (s46 * (-t616 - t626 + s56) - t733) - t349 + t41
     #1 * t111 + t192 * t766) + t306 * (s46 * (s26 * (s25 - t440) - t634
     #) + t762) - t66 * (s56 * (-t167 - t182) + t644 + t706 - t727 + t16
     #1 * (-t550 - s45)) - 9.0_dp * s56 * (-t764 - t763) - 9.0_dp * t712
     # + 9.0_dp * t765 - 9.0_dp * t222 * (t626 - s26) + 7.0_dp * s46 * (
     #t717 + t174) + 7.0_dp * t163 + 13.0_dp * t616 * t555 + 13.0_dp * t
     #694 + t308 + 12.0_dp * s46 * (t756 + t767)
      t135 = s12 * (s12 * (-t240 * t66 - t28 * (t90 + t660) - epinv * (t
     #28 * (s25 * (t119 - t80) - s46 * t88 - t93) - t543 + t540 - t142) 
     #- s56 * (s16 - s34 + s45 - s46) - t71 + t95 + t464 * t539 - t752 -
     # t560 * epinv) + t135) - t284 * (t747 + t748 + t746 - t111 - t751 
     #+ t750 - t654)
      t308 = intHLs16s25s26s34s56x1213D6eps1()
      t437 = s26 - s56 - s46
      t450 = -t437 + t80
      t539 = t119 + s26
      t540 = s16 * (t617 + t414 - t753 - s25 - s45 - t483 + s16) + s25 *
     # t450 + t139 * (t2 * t539 - s34)
      t543 = t327 * epinv
      t555 = t142 * t2
      t641 = t425 + t424
      t644 = t404 * t136
      t694 = t404 - s45
      t712 = t457 * s25
      t717 = t422 * s45
      t769 = t186 * t2
      t770 = t149 + t182
      t771 = -t447 - t414 + s34
      t215 = s13 * (-t17 * t417 + t28 * (s45 * (-t617 - s34) + t709) - t
     #66 * t641 + t426 + t427 + s16 * (t28 * (-t617 + s34) - t66 * t414 
     #- s26 + s45 - t640 - s16) + s25 * (t263 - t753 - t722 - t509 - t13
     #3) + s45 * (t404 - t414 - t400) - t252 - t415 + t555 - t465) + s16
     # * (s16 * (t108 - t543) - s25 * ((-1.0_dp + t430) * s26 - t627 * s
     #34 - t413 * (-t509 + s45) + s25) - s26 * (t124 * t766 - t133 * t2 
     #+ t17 * t411 - t199 + t483 + t753) - t305 * (epinv * (s45 + t523 +
     # t538) + t483 - t523 - t538)) + s26 * (t149 - t182) + s45 * (epinv
     # * t215 + t113 - t142 + t149) + s46 * t770 + s56 * t770 + t167 * t
     #771 + t466 + t734 + t216 * t2
      t191 = s12 * (s12 * (-epinv * (s25 * (t119 + t140) - t111 + t113 -
     # t179 - t491 + t195) - s26 * t314 - t111 + t113 - t491 - t660 + t4
     #64 * t191 - t447 * s12) - t124 * (s26 * (t559 - t417) + t161 * t51
     #9 + t197 * (t617 + s26)) + t16 * (t454 + t769) + t17 * (s26 * (t20
     #5 * t404 + t184 + t514) - t206 * epinv) - t28 * (s25 * (s45 * t694
     # - t542 + t630) + s46 * (-t425 + t158) + t71 * (t447 + t613 - t446
     #) + t176 - t461 - t717 + t161 * (-t712 + t613)) + t66 * (s25 * (s5
     #6 * (t404 + s45 + s56) + t158) + s46 * (t174 + t419 + t644) + t113
     # * (-t613 + s25) + t741 + t762) + t432 + t215 + 7.0_dp * t400 * t1
     #13) + t284 * t540
      t195 = intHLs16s25s26s34s56x1231D6eps1()
      t206 = -s34 + s45 - s46
      t215 = t28 * (t91 + t244 + t151 - t653) - t614 * t517 - t655 + t65
     #7
      t252 = -t414 - s34
      t426 = t457 * s26
      t432 = t404 * s46
      t454 = t432 + t179
      t461 = epinv * t16
      t491 = s34 * t539
      t514 = t122 * t17
      t517 = s34 + s45 + s56
      t539 = s16 * s25
      t657 = t172 * s16
      t129 = -epinv * (s56 * (s46 * t393 - t142 * t17) + t124 * (s26 * (
     #-t93 + t71) + s46 * (-t387 - t530)) + t28 * (t474 - t532 + t382 + 
     #t114 + t72 + t186) + t66 * (s56 * (-t129 + t111) + t384) - t276 + 
     #t535 + s25 * (t119 * t127 - t124 * t491 + t468 + t514 - t533) - s4
     #5 * (s26 * t254 - t71) + t325 + t534) + s13 * (-t17 * t454 + t28 *
     # (s45 * t252 + t632) + t66 * (s46 * t646 - t178 * t411 - t440 * t8
     #9 + t442 + t634) - t683 + s34 * (-t426 + t408) - t149 - t415 + t42
     #8 + t614 * t524 - t122 * (t3 - t461) - t465) + s34 * (s16 * (-s25 
     #- s26 - s45 + s56 - s46) + t167 - t184 + t672 - t681)
      t276 = t279 * s16
      t129 = s12 * (s12 * (-epinv * (-t28 * (-t224 + t93) + t515 + s25 *
     # (t106 - t80) + t111 - t338 - t90) - s56 * t206 - t660 - t90 + t95
     # + t464 * t304 - t682 - t651 * s12) + t124 * (s46 * (t174 + t673) 
     #+ t657) + t16 * (s56 * (t689 + t168) + t169) - t17 * (s46 * (t366 
     #- t539) + s56 * (-s56 * t692 - t539) - t113 * t319 - t471 - t500) 
     #+ t28 * (s25 * (s25 * t517 + t167) + s26 * (-s16 * t698 + t267) - 
     #t112 + t114 + t382 - t699 + t727) + t3 * (s56 * (-t90 + t179) + t2
     #76) + t66 * (s16 * t167 + s45 * (s34 * t213 + t539) + t186 + t222 
     #* s16) + t702 + t703 + t704 + t129) + t284 * t215
      t304 = intHLs16s25s26s34s56x1411D6eps1()
      t393 = t101 * t2
      t428 = t141 * epinv
      t465 = intHs16s25s26s34s56x1123D6eps1()
      t515 = s16 + s26 - s45
      t524 = t515 - t483
      t533 = t541 * s34
      t534 = s12 * t148
      t535 = s12 * (-t28 * (-s26 * t685 + t146) - t66 * t533 - s13 * (t2
     #31 - t330) + s16 * t358 + t113 + t168 + t317 + t672 + t71 - t534) 
     #+ t524 * t231 * t13
      t673 = intHL0s25s260s56x1022D4eps0()
      t683 = t139 * t148
      t685 = t148 * t287
      t692 = t46 + t123 + s25
      t450 = s12 * (s25 * t28 * t692 - s13 * t536 + t168 * t66 + t685) +
     # t13 * (-s25 * t450 + t436 + t683)
      t536 = intHL0s25s26s34s56x1120D2eps1()
      t699 = -s12 * t47 - t49
      t702 = intHs16s25s26s34s56x1113D4eps1()
      t703 = t231 + t579
      t704 = s25 + s56 - s13 + s46
      t727 = s16 - s34 + s45
      t734 = t28 * t704
      t741 = s12 * (t727 + t734 + t127 + s12) + t284
      t762 = intHs16s25s26s34s56x1132D6eps1()
      t772 = t28 * t87
      t773 = s25 * t148
      t326 = -s12 * (t17 * (t244 + t222) + t28 * (t773 + t360) + t66 * (
     #-t533 + t151) - s13 * (-t217 * t66 + t772) + s16 * (t28 * t358 + s
     #56) - t93 + t122 * t124) - t326 + t609 * t148
      t533 = intHL0s25s260s56x1021D2eps1()
      t774 = intHs160s26s34s56x1031D4eps0()
      t775 = intHs16s25s26s34s56x1114D6eps0()
      t776 = intHLs16s25s26s34s56x1411D6eps0()
      t777 = intHLs16s25s26s34s56x1311D4eps0()
      t778 = intHs160000x0111D0eps1()
      t779 = intHL0s25s26s34s56x1130D4eps1()
      t780 = intHL0s25s26s34s56x1120D2eps0()
      t781 = intHL0s25s26s34s56x1210D2eps0()
      t782 = intHL0s25s26s34s56x1310D4eps0()
      t783 = intHLs160000x0112D2eps1()
      t784 = t28 + t8
      t785 = t4 * t66 * t8 + t5 * t784
      t786 = intHLs16s25s26s34s56x1112D2eps0()
      t787 = intHLs16s25s26s34s56x1113D4eps0()
      t788 = intHs16s25s26s34s56x1312D6eps1()
      t789 = intHs16s25s26s34s56x1222D6eps1()
      t790 = intHs16s25s26s34s56x1213D6eps1()
      t791 = intHs16s25s26s34s56x1310D4eps1()
      t792 = intHL0s25s260s56x1031D4eps1()
      t793 = intHs16s25s26s34s56x1112D2eps0()
      t794 = intHs16s25s26s34s56x1130D4eps1()
      t795 = intHLs16s25s26s34s56x1211D2eps0()
      t796 = intHs16s25s26s34s56x1311D4eps0()
      t797 = intHLs16s25s26s34s56x1114D6eps1()
      t798 = intHs16s25s26s34s56x1211D2eps0()
      t799 = intHs160s26s34s56x1022D4eps0()
      t800 = s12 + s16 + s26 - s45
      t801 = s12 + s16
      t802 = (t101 - t483) * s45
      t803 = (t800 - t483) * s56
      t804 = s16 * t801
      t805 = (t564 + s26) * s26
      t806 = -s12 + s34
      t807 = s34 * t806
      t808 = -t101 * t330 + t28 * t807 - t802 + t803 + t804 + t805
      t809 = intHs16s25s26s34s56x1220D4eps0()
      t810 = -s13 + s16
      t811 = -s13 + s26
      t812 = s25 - s16 - s26 + s34
      t813 = (-t811 - t327) * s34
      t594 = s26 * (t810 + t327 + s26) - s45 * t812 + s56 * (s16 - s34 +
     # s45 - s13 + s46 + t569 + s56) + t591 - t592 + t594 + t813
      t814 = intHLs16s25s26s34s56x1141D6eps0()
      t815 = t232 - s13
      t816 = s16 * t232
      t817 = s12 * (-s13 + s25 + s26 + s45)
      t818 = s16 * s45
      t819 = s12 + s16 - s34
      t820 = (-t65 - t63) * s34
      t821 = t801 * s13
      t822 = t28 * t819
      t823 = t93 * t568
      t824 = s45 * (s16 * t815 + t818 - t93) + (-s13 * t64 - t28 * (t93 
     #- t818) - t332 + t816 + t817) * s46 + ((-s34 + t563) * s45 + t232 
     #* t801 + t820 - t821 + t822 * s46) * s56 + t113 * t819 + t71 * t81
     #9 - t823
      t825 = intHLs16s25s26s34s56x1132D6eps0()
      t826 = s16 - s26
      t827 = -t28 * t806
      t828 = s25 * t801
      t829 = s25 + s12 - s13
      t830 = t17 * t64
      t831 = t28 * t64
      t832 = (t581 + s26) * s26
      t833 = t70 + t832
      t820 = s26 * t833 + s45 * (s45 * t826 - s13 * t826 + s16 * (t65 + 
     #s25) - s26 * (s25 + t483)) + s46 * (s46 * (-t143 + t822) + t28 * (
     #-t821 + t828 + t820 + t111) + (t63 + t194 + t579) * s26 + s45 * (-
     #t568 + t827 + t194 - s45)) + s56 * (s56 * (s12 + s26 - s45 + t831)
     # + t28 * (t828 - t592 + t820 + t111) + (t829 + t194) * s26 + s45 *
     # (-t568 - t483 + t194 - s45) + (-t143 * t28 + t577 + t830) * s46 -
     # t608) - t823 * t28
      t823 = intHLs16s25s26s34s56x1123D6eps0()
      t834 = t28 * t223
      t835 = t14 * s46
      t836 = s12 * s56
      t817 = s26 * (s45 * (-s25 + s16 - s34 + s13 - s45) - t593 + t70) +
     # s46 * (t835 + t817) + t111 * t581 + t12 * (-s46 * (-s25 - s45 + s
     #13 - t65) + s56 * (-s13 + t136 + t834 + s56)) + t112 + t836 * t205
      t837 = intHs16s25s26s34s56x1321D6eps0()
      t838 = s12 + s13
      t839 = t838 + t563
      t840 = s13 + s25
      t841 = s13 * (-s13 + t565)
      t842 = t17 * s16
      t843 = t66 * t101
      t844 = t28 * t801
      t845 = t28 * t810
      t846 = s26 + s13
      t847 = t838 * s13
      t848 = t66 * t111
      t849 = s12 - s16
      t850 = t28 * t63
      t851 = t66 * s13
      t852 = t810 * s25
      t853 = t28 * t111
      t854 = s25 + s16 - s34
      t855 = s45 - s13
      t856 = t28 * t319
      t857 = s12 + s13 - s16
      t858 = t17 * t113
      t859 = t124 * t223
      t860 = (-s13 + s16 + s25 + s26) * s45
      t815 = t28 * t72 + s26 * (t842 * s25 + s26 * (-t566 + t327) + t579
     # * t840 + t841) + s34 * (-s13 * t815 - s25 * (t843 + t579) + t267)
     # + s45 * (-(t844 + s13) * s25 - (s25 - t845) * s26 - (t846 - t327)
     # * s34 + t847 + t848 + t860) + s46 * (s46 * (t28 * t855 + s46 + t1
     #27 + t854) + (-t849 + t850) * s26 + s45 * (t319 + t201 + t563 - t8
     #51 + s45) - t241 + t847 + t852 + t853) + s56 * (s56 * (t66 * t855 
     #+ t319 + t563 + t859) + s13 * t857 - (s12 - t567 + t851) * s26 + s
     #45 * (t66 * t810 + s45 + t319 + t507) + (t124 * t855 + t194 + t496
     # + t856) * s46 - t241 + t852 + t858 + t848) - t359 * t801 + t828 *
     # t839
      t861 = intHLs160000x0121D2eps0()
      t862 = intHLs160000x0122D4eps0()
      t863 = intHL0s25s26s34s56x1220D4eps0()
      t864 = s25 + s26 + s34 + s56
      t865 = s12 * s34
      t47 = t47 * t865 + t49 * t864
      t866 = intHL0s25s260s56x1012D2eps1()
      t867 = intHL0s25s260s56x1013D4eps1()
      t868 = intHL0s25s260s56x1010D0eps0()
      t869 = intHL0s25s260s56x1011D2eps0()
      t870 = intHL0s25s260s56x1020D2eps0()
      t871 = intHL0s25s26s34s56x1110D2eps0()
      t872 = intHs16s25s26s34s56x1113D4eps0()
      t873 = intHLs16s25s26s34s56x1311D4eps1()
      t874 = t29 * t101
      t875 = t429 + t874 + t80
      t876 = intHs160000x0111D0eps0()
      t877 = intHL0s25s26s34s56x1130D4eps0()
      t878 = intHs160s26s34s56x1022D4eps1()
      t879 = intHs16s25s26s34s56x1220D4eps1()
      t880 = intHs16s25s26s34s56x1321D6eps1()
      t881 = intHs16s25s26s34s56x1231D6eps1()
      t882 = t66 * s25
      t883 = s13 - t327
      t884 = t3 * s16
      t885 = -s16 + s25
      t886 = t28 * s13
      t887 = t838 * s16
      t888 = t28 * t144
      t889 = t16 * t111
      t890 = s13 + s34
      t891 = s25 - s16 - s13
      t892 = t16 * s16
      t893 = s16 + s46
      t894 = t306 * s45
      t895 = 11.0_dp * s16
      t896 = 7.0_dp * s16
      t897 = t124 * s13
      t898 = 11.0_dp * s26
      t899 = s56 * (t896 - t897 + t343 + t898 + t133 + s12 + t327 - t372
     #) + s26 * (t374 + t894) - t124 * (-t320 - t279) + t17 * s25 * t893
     # + t3 * (t676 + t113) - t66 * (s25 + s16 + s45) * s34 + t888 - s12
     # * (t242 * t28 + s12 - t133 - t891) - s13 * (t223 * t3 - s13 + t13
     #6 - t330 + t892) + t149 - t167 + 7.0_dp * t111 + 7.0_dp * t818 - 7
     #.0_dp * t197 + t895 * s46
      t320 = t17 * t72 + t28 * (-s13 * t144 + t114) + s13 * (-t887 - t60
     #8) + (t801 * (t842 + s13) + t75) * s25 + s26 * (-s26 * (t577 + s13
     # - t490) - s16 * (-t357 + t851) + s12 * (t882 + s13 - t563 - s12) 
     #+ s13 * t63) + s34 * (-s34 * t883 + s12 * (t846 - t882) + s13 * (t
     #65 - s25 + t194 + s13) - t537 * t101) + s45 * (s45 * (-s13 + s25 -
     # s26 + s34) - (s13 + t194 + t579) * s25 + (t577 + t884 - s13 - t32
     #7) * s26 + s34 * (t66 * t885 - s12 + s34 - t201 - t886) + t847 + t
     #887 + t888 + t668) + s46 * (s46 * (-t66 * t890 + s12 + t138 + t251
     # + t842) + t124 * t320 + t17 * s16 * t136 + t28 * (-t158 + t144) -
     # t66 * s34 * t565 + s12 * (s16 + s25 - s34 + s45) + s13 * (t178 * 
     #t28 + s13 - s25 - t194 - t201) + t149 + 7.0_dp * s26 * t584 + t889
     #) + s56 * t899
      t668 = intHs160000x0121D2eps1()
      t846 = t12 * (t13 + s12) + t75
      t899 = intHs160000x0112D2eps1()
      t900 = (-t12 - s12) * s12
      t901 = t49 + t900
      t902 = intHLs160000x0112D2eps0()
      t903 = intHLs160000x0122D4eps1()
      t904 = intHL0s25s26s34s56x1220D4eps1()
      t905 = intHL0s25s260s56x1011D2eps1()
      t906 = intHL0s25s26s34s56x1110D2eps1()
      t907 = intHL0s25s260s56x1010D0eps1()
      t908 = intHL0s25s260s56x1020D2eps1()
      t909 = intHL0s25s260s56x1022D4eps1()
      t910 = intHLs16s25s26s34s56x1211D4eps1()
      t911 = t29 * s13
      t912 = t28 * (t429 + s34 + s45) + t628 + t874
      t913 = t35 * s12
      t914 = intHs16s25s26s34s56x1112D4eps1()
      t915 = intHLs16s25s26s34s56x1121D2eps1()
      t916 = s16 + s25 + s26
      t917 = t168 * t28
      t918 = s12 - s13 + s16
      t919 = (t918 + t327 + s26) * s26
      t920 = s12 + s16 + s25 + s26 - s34
      t921 = t920 * s46
      t922 = t916 * s45
      t923 = t267 * t568
      t924 = (s25 * (t918 + s25) - s34 * t811 + (t28 * t916 - t806) * s4
     #5 - t821 + t919 + t921) * s46
      t816 = (-s13 * t916 + s25 * t110 + t111 + t816 + t917 + t922) * s4
     #5
      t678 = (s56 * (-t28 * t678 + s56 + t918) + s25 * (s34 + t918 + s25
     #) + s45 * (t602 - s13 + t563 + s45) + s46 * (t28 * (s12 + s16 + s4
     #5) + t602 - t890 + s46) - t821 + t919) * s56
      t925 = epinv * s12
      t926 = intHLs16s25s26s34s56x1141D6eps1()
      t927 = intHLs16s25s26s34s56x1132D6eps1()
      t928 = t101 + t882
      t929 = t28 * t566
      t930 = t28 * t811
      t931 = s12 + s26 - s34
      t932 = -t28 * t920 + s45
      t933 = t66 * t801
      t934 = t928 * s45
      t935 = t932 * s46
      t936 = intHLs16s25s26s34s56x1123D6eps1()
      t937 = -s13 + s25 + s34
      t938 = (t14 + t327) * s46
      t939 = (s12 + s16 + s25 + s26 - s34 - s45) * s46
      t940 = t28 * s25 * (s16 * (-s13 + s25) + (t581 + s26) * s26) - t66
     # * s25 * (t593 + t167) + epinv * t817 - s46 * (s25 * (s12 - s16 - 
     #s13 + s25) + s26 * (t918 + s26) - (t811 - t327) * s34 + s45 * (s12
     # + s16 - s34 + s13 + t882 - s45) - t821 + t939) - s56 * (s25 * (t3
     #30 + t152) - s25 * (s26 - t566 + t563 - s25) + t938 + t192) - t174
     # * (t66 * t937 + s26 - t563)
      t941 = intHLs16s25s26s34s56x1131D4eps0()
      t942 = s34 * t13
      t943 = s12 * (-s13 * t148 + s26 * t351 + t130 - t436 + t542 + t561
     # - t92) + t942 * t231
      t944 = intHs16s25s26s34s56x1123D6eps0()
      t945 = intHs16s25s26s34s56x1121D2eps0()
      t946 = -s26 - s34 - s45
      t947 = s34 * t141
      t948 = intHs16s25s26s34s56x1122D4eps0()
      t949 = s34 + s45 - s56
      t950 = s16 * t506
      t951 = t949 * s34
      t517 = t28 * (-s26 * (s45 + t483) - t172 + t950 + t951) - s46 * t5
     #17 + t111 + t144 - t71
      t952 = t101 - t372
      t953 = t108 * s26
      t954 = s34 * t205
      t955 = intHs16s25s26s34s56x1131D4eps0()
      t956 = s16 + s26 + s56
      t957 = t956 - t483
      t958 = t358 * s26
      t959 = intHLs16s25s26s34s56x1122D4eps0()
      t960 = s26 + t327
      t961 = s16 * t960
      t962 = t287 * s56
      t963 = -t482
      t964 = t139 * t106
      t965 = s12 * (t46 + s26)
      t966 = t12 * t231
      t967 = t966 * t108
      t968 = intHs16s25s26s34s56x1132D6eps0()
      t969 = intHs16s25s26s34s56x1121D2eps1()
      t970 = intHs16s25s26s34s56x1122D4eps1()
      t971 = t400 + t408
      t972 = t404 * t17
      t102 = -t28 * epinv * (s34 * (-s34 + s56 - s45) + t172) + s16 * (-
     #t28 * t971 + t102 - t972) - s26 * (epinv * (t59 + t372) - t102 + t
     #400) - s46 * (-t613 - t446 + t411) - s56 * (-s34 - s45 + t411) - t
     #423
      t973 = t106 * s25
      t974 = t404 + s46
      t975 = t756 * t16
      t976 = s12 * (-t124 * s56 * t974 - t17 * (s34 * (-t414 - t407 - s4
     #5 - s56) + t172 - t95) - t28 * (s34 * (-s26 + t408) - t459 - t90 +
     # t973) - t66 * t151 + epinv * (t953 + t172) + s13 * (-t17 * t613 +
     # t235 + t538 - t874) + s16 * (epinv * (t287 - t483) - s56 + t640 -
     # t128) - t975) + t102 * t13 - t609 * (t404 - s26 - s56)
      t977 = intHs16s25s26s34s56x1131D4eps1()
      t978 = t874 + t411
      t979 = -t613 + s46
      t980 = -s25 + t411
      t981 = intHLs16s25s26s34s56x1122D4eps1()
      t982 = -t698 + t235
      t983 = s12 + s45 + s56
      t735 = -s25 * t735
      t984 = t271 * s46
      t985 = t801 * s56
      t482 = s25 * (s56 * (-t527 - t354) - t257 * s46) - t17 * (s25 * (-
     #s45 * t204 + t71) + t222 * t314 + t267 * (s26 + s45 + s56 + s25)) 
     #+ t28 * (s25 * (s46 * (-s12 - s34) - t836) + s26 * ((-s12 - s16 + 
     #s34 - s45 - s46) * s46 - t71) + s46 * (-t675 + t158) + s56 * (-t81
     #8 + t197) - t186 + t674 + t90 * s12) - t66 * (s46 * (t985 + t182 +
     # t818) + s56 * t182 + t209 + t179 * t136) + epinv * (s12 * (-t28 *
     # (t962 + t542 + t161) + s26 * (s25 - s34 - s46) + t961 + t560) - s
     #13 * (s16 * (t482 + s16) + t964 - t965) + t967) + s13 * (s46 * t28
     # * t931 + s16 * t982 + s56 * t983 - t157 + t735 + t984) - s16 * (t
     #167 - t182) - s26 * (s45 * t584 - t182 - t539 + t985) - s46 * (s45
     # * t271 + t539) - s56 * (s56 * (s12 + s16 + s56) + t167) + t111 * 
     #t698 - t113 * t271
      t674 = intHLs16s25s26s34s56x1131D4eps1()
      t675 = t88 * s26
      t698 = intHs16s25s26s34s56x1231D6eps0()
      t985 = intHs160s26s34s56x1021D2eps1()
      t986 = t103 + s12
      t987 = t573 * s45
      t988 = t14 * s56
      t989 = intHs16s25s26s34s56x1120D2eps1()
      t990 = (-t108 + s13) * t27
      t991 = intHs160000x0121D2eps0()
      t992 = intHs160000x0112D2eps0()
      t993 = intHLs160000x0111D2eps0()
      t994 = intHL0s25s260s56x1012D2eps0()
      t995 = intHL0s25s260s56x1013D4eps0()
      t996 = intHL0s25s260s56x1021D2eps0()
      t997 = intHs16s25s26s34s56x1121D4eps0()
      t998 = -s13 - s16 + s25 + s34 + s46
      t999 = -t28 * t155
      t1000 = intHLs16s25s26s34s56x1121D4eps0()
      t1001 = s16 + s26 + s34 + s45
      t1002 = t1001 + t509
      t1003 = t806 * s56
      t1004 = t946 * s25
      t1005 = intHs16s25s26s34s56x1112D4eps0()
      t589 = s12 * (-t589 - s12) + t1002 * t13
      t1006 = intHLs16s25s26s34s56x1211D4eps0()
      t1007 = intHs16s25s26s34s56x1121D4eps1()
      t1008 = intHLs16s25s26s34s56x1121D4eps1()
      t1009 = -s25 - s26 + s34
      t1010 = t1009 * s26
      t1011 = s26 * t718
      t1012 = t29 * s16
      t1013 = 0.1e1_dp / t14
      t1014 = 0.1e1_dp / t51
      t189 = 0.1e1_dp / t189
      t1015 = 0.1e1_dp / s25
      t1016 = 0.1e1_dp / t27
      t1017 = 0.1e1_dp / s12
      t1018 = t905 + t906
      t1019 = -t907 + t908
      t1020 = t994 + t996
      t1021 = t780 + t781
      t1022 = 0.1e1_dp / ecossin ** 2.0_dp
      t1023 = 0.1e1_dp / gw ** 2.0_dp
      t1024 = t699 * t1014
      t1025 = t2 * t1014
      t1026 = (-t846 * t991 + t901 * t992 + t997 * (s12 * (-t998 - t999 
     #+ s12) + t188) - t1005 * t589 - t1007 * (-s12 * (t28 * (epinv * t1
     #55 + s25 + s46) + t133 + t234 + epinv * (-s25 - s34 - s46) + s26 +
     # s34 - t614 + t911) + t2 * (t188 + t75)) - t85 * t778) * t1015 * t
     #9
      t1027 = t2 * t789
      t1028 = t1027 * t11
      t514 = t19 * (t1015 * (t959 * (s12 * (-t28 * (t962 + t161) + s25 *
     # (-t80 + s26) + s26 * t166 + t961 + t560) + s13 * (s16 * (t963 - s
     #16) - t964 + t965) + t967) + t1000 * (-s13 * t1002 + t28 * (t192 +
     # t161 + t113 + t71 + t184) + t66 * (t272 + t121) + t1003 - t1004 +
     # t111 + t117 + t158 + t167 + t197 + t514 + t93) + t589 * t1006 + t
     #1008 * (t28 * (s45 * (s25 + s34 + s45 + t440) + s46 * (t447 - s26 
     #+ s34 + s45) + s56 * (t447 + s12 - s26 + s34 + s45) + t267 + t630 
     #+ t634 + t1010) + t66 * epinv * (s45 * t718 + t1011) + epinv * (-t
     #1004 + t1003 + t167 + t93 + t111 + t158 + t197) - s13 * t912 + t10
     #12 * t108 + t637 * t17)) + t1024 * t1020 + t1025 * t1015 * (t1019 
     #* t699 + t450 * t909 - t47 * t904) - t1024 * t1015 * t1018 * t29) 
     #+ t1026 + t1028 * t278 * t1015 * t189
      t560 = t19 * t1014
      t589 = t1022 * t1013
      t630 = t589 * t1023 * propW34
      t961 = intHs160s26s34s56x1011D2eps0() + intHs16s25s26s34s56x1110D2
     #eps0()
      t212 = t203 * (s12 * (s12 * (-t118 + t117 + t337 + t292 - t111 + t
     #113 + t293) + t124 * (-t179 * t204 - t176 + t212 - t221) + t16 * (
     #s56 * (t224 + t350) + t186) - t17 * s34 * (s26 * t346 + t71) + t28
     # * (s46 * (t149 + t174) + s56 * (t347 + t149) + t114 + t163 + t72 
     #+ t348 + t349) + t66 * (s25 * (t158 - t142 + t113 + t71 - t197) + 
     #s45 * (-t219 + t111) + t162) - 9.0_dp * t142 * t223 - s13 * (-t17 
     #* (t339 - t179) + t28 * (t91 + t172 + t149) + t66 * (-t338 + t158 
     #+ t144 + t113) - s16 * (t124 * t165 - t127 + t136 - t234) - t168 -
     # t301 + t340 + t185) + s16 * (t344 + s26 * (-t124 * t178 + t342 + 
     #t343) - t305 * (t124 * t217 - t199) - t345 * t341) - s45 * (-t182 
     #- t113) + t167 * t213 + t112 + t352 + t353 + t222 * (t354 - t309 +
     # t355)) + t284 * t335) + t394 * t265
      t265 = t11 * t356
      t108 = t1015 * (t1013 * (-t914 * (s12 * (-t28 * (s25 + s56 + s46 -
     # t408) - t133 + epinv * (s25 + s56 + s46) - s26 - t613 + t614 - t9
     #11 + t37) + t13 * (t28 * t429 - t446 - t448 - t613 + t874)) - t945
     # * (s12 * (-s46 * t946 + s56 * t254 + t122 * t28 + t351 * t63 + t1
     #13 + t71 - t947) + t562) - t948 * (s12 * (-s34 * (t495 + t507) - t
     #17 * t954 + t28 * t518 - s13 * t952 + s16 * (t108 - t483) + t172 +
     # t953) + t13 * t517 - t609 * s34) - t969 * (s12 * (t28 * (s46 * (s
     #26 + t411) + t151 + t244) + t98 - s13 * (t619 - t613 + s16 + s26 +
     # t411 + t46) + s16 * t541 + s25 * (-t613 + s26) + s46 * (t712 + t6
     #21 - t613 + t440) + s56 * (t712 - t613 + t440) - t458 + t625 + t63
     #4 + t756 - t534) + t85 * t27) - t970 * t976 - t189 * t212) - t942 
     #* t793) * t9
      t212 = t1015 * t1013
      t175 = t212 * (-t910 * (s12 * (t28 * (s25 - s34 + s56 + s46 - t408
     #) + t200 - epinv * (s25 + s34 + s56 + s46) + s16 * t35 + t911 + t9
     #13) - t13 * t912) + t981 * t482 + t189 * (t230 * (t28 * (-s34 * t2
     #31 * t284 + t51 * t76) - s12 * (s12 * ((-t236 + s13) * s13 + t158 
     #+ t71 + t238 + t239 - t247 - t249 + t245 + t246 + t248) - t124 * s
     #45 * (t267 - t111 + t266) + t17 * (-t280 + t282 + t281 + t186) + t
     #264 * t3 - t28 * (s45 * (t93 + t71) + s46 * (s45 * t268 - t71) + t
     #113 * t143 + t269 + t270) + t66 * (-t149 * t232 + t221 - t277) + s
     #13 * (t259 + t185) + s25 * (-t272 + t121) + s34 * (t273 - t167) - 
     #t187 + t274 + t275 + t283 + t199 * t179)) + t285 * t175)) * t19
      t259 = t2 * intHs16s25s26s34s56x1210D2eps1()
      t266 = intHs160s26s34s56x1020D2eps0() * t1015
      t270 = t942 * t795
      t274 = t1022 * t1023 * propW34
      t277 = t868 - t869 - t870 - t871
      t283 = t2 * t33
      t293 = s34 * t40
      t301 = t2 * t1016
      t337 = t301 * (-s34 * t26 * t38 + t594 * t879 - t808 * t878)
      t35 = t34 * t35
      t340 = t36 * (t37 + s16 + s26 - s34 - s45)
      t346 = t1015 * t1017
      t18 = t346 * (t13 * (-t22 * (s34 * t340 + t12 * t35) - t9 * t48 * 
     #t876) + t9 * (-t668 * t846 + t899 * t901) * epinv + epinv2 * (t5 *
     # (t4 * (t17 * (s12 * (s25 * t600 + s26 * (t601 - s13 + t579) - t59
     #2) + s34 * (s12 * (-t602 + s13) - t570 + t593) + (s12 * t604 - t57
     #4) * s45 + s56 * (s56 * (-t581 - t588 + t605 - s56) + (s13 + t605)
     # * s16 + (t577 - s16) * s25 + s26 * (-t581 + t606 - s26) + s34 * (
     #-t607 + s34) + s45 * (-s25 - s12 + s34 + s13 - t127 - t563 - s45) 
     #+ t587 - t608 + t609) - t595 - t596 - t597) - 12.0_dp * t599 - t18
     # * t27 * (-s12 * t568 + s45 * (t575 + s45) + s46 * (t576 + s46) + 
     #s56 * (-t28 * t610 + s56 + t572) + t593 + t598 * t590)) - t20 * (s
     #12 * (s13 * t27 + s56 * t368 - t347 - t370 + t387 - t561 - t71 + t
     #92) + t562)) + t611) * t1014 * t1016)
      t131 = t1015 * (t1017 * (t13 * (t12 * (t283 - t10) + t293 * t42) -
     # t612 * (s12 * (s12 * t395 + s56 * (s46 * (t544 + t649) + 11.0_dp 
     #* t179) + t124 * (s46 * (-t2 * t93 + t419) + t174 * t154 + t93 * t
     #155) + t16 * s26 * ((s25 + s26) * s46 + s56 * (epinv * (s34 - s45)
     # + s25) + t168) - t17 * (epinv * t177 + t165 * t192 + t501 + t645 
     #- t72) - t28 * (s45 * (s46 * t646 + t459) + s46 * (epinv * (t149 -
     # t168) - t350) + t167 * (-epinv * t166 - s25) + t71 * (-t409 - s25
     #) + t647 + t459 * t232) + t306 * s56 * (t122 + t161) - t66 * (s25 
     #* (s34 * (t635 - s34) + t636) + s56 * (s56 * (-t446 + t440) + t158
     # + t623) + t221 - t648 * t154 + t616 * (-t160 + t161)) - t131 + 7.
     #0_dp * s56 * t463 + 7.0_dp * t161 * t223 + 7.0_dp * t388 - 9.0_dp 
     #* t379) + t396 * t100 * t104) * t189) + t337) * t11
      t18 = t630 * (t1017 * (t1014 * (t1015 * (t277 * t699 + t47 * t863)
     # - t2 * t699 * t866) + (-t467 * (-t28 * (t651 * t76 + t284 * (s16 
     #* (t613 - s45 - s56 - s46) - t149 - t222 - t224 + t339 - t71 - t91
     # + t677 - t404 * t106)) + s12 * (s12 * (epinv * (t238 + t239 + t15
     #8 + t71 - t247 - t249 + t245 + t246 + t248) + s13 * (-epinv * t236
     # + s25 + s34 + s45 + s56 + t127 - t464) - s56 * t552 - t161 + t680
     # + t681 - t682) + t211)) - t526 * t96 - t220 * (s12 * (s12 * (t124
     # * t447 * t154 - t28 * t559 - epinv * (t291 + t290 + t292 - t111 +
     # t113) + s16 * (t2 * t287 + s25 - 7.0_dp * t447) - s25 * t254 - t1
     #11 + t113 + t290 + t291 + t292 - t660 - t286 + t248 * epinv) + t32
     #3) + t284 * t553 - t328 * epinv) + t443 * t135) * t1015 * t189) * 
     #t19 + t18 + t131)
      t47 = t1015 * (s34 * (intHs160s26s34s56x1012D2eps0() + intHs160s26
     #s34s56x1021D2eps0()) + t29 * (-intHs160s26s34s56x1011D2eps1() - in
     #tHs16s25s26s34s56x1110D2eps1()))
      t85 = t9 * t99 * (((t28 * (t120 + t121 + t122) + (-t105 - t107 + s
     #13) * s13 + t111 + t113 + t71 + t115 + t116 - t117 + t118) * s12 +
     # t124 * s46 * t173 + t16 * t179 * t178 - t17 * t177 + t28 * (t149 
     #* t164 + t165 * t71 + t166 * t167 + t169 + t170 - t171) + t66 * ((
     #t159 + t93) * s25 + (t157 + t158 - t71) * s26 + (t160 - t161) * s4
     #6 - t162 - t163) + (t28 * (t151 + t149) + t66 * (t147 - t146 + t14
     #4) - s16 * (t138 - t137 + t126 - t127) + s25 * (t119 - t140) - t14
     #1 * s45 + t113 - t142 + t153 + t185) * s13 + s16 * t134 + (t183 - 
     #t182 + t113 - t184) * s34 - t112 - t186 + t187 + t168 * t181 - t17
     #4 * t136) * s12 + t188 * t100 * t104) * t189 * t1015 - t11 * t84 *
     # (-s12 * t14 + t85)
      t85 = t630 * (t1017 * t85 + t19 * (t1017 * (t83 * (t49 - t82) + t1
     #014 * t1015 * (-t450 * t673 + t536 * (t404 * t699 + t49 * t51) - t
     #533 * (s12 * (s13 * (-t407 - s26 - s56) + s25 * (t2 * t692 - t753)
     # + t685) - t13 * (s25 * (t613 - s26 + s45 + t722) - t436 - t683)) 
     #+ t613 * t43 * t50)) + t86 * (t28 * t97 - s13 * (t46 + t87) - t88 
     #* s16 - t90 + t91 - t92 + t93 + t98) * t1015) + (t9 * (t1015 * t2 
     #* t25 - t53 - t54) + t47 * t11) * t1016 * t26 - t58 * t22 * t60 * 
     #t1015 + t346 * (t12 * (-t15 * t22 - t9 * (epinv * t32 + t31)) - t2
     #3 * t11 * t24) * t13)
      t88 = t2 * intHLs160000x0113D4eps1() - intHLs160000x0113D4eps0()
      t96 = t346 * t589 * propW16 * propW34 * s16 * t785
      t98 = epinv * intHLs160000x0121D2eps1()
      t105 = t29 * intHLs160000x0111D2eps1()
      t118 = s16 * (t2 * (-t783 + t903) + t861 - t862 + t902 + t98) - t9
     #93 - t105
      t120 = t2 * intHs16s25s26s34s56x1141D6eps1()
      t131 = t2 * intHs16s25s26s34s56x1114D6eps1()
      t134 = -t131 + t775
      t24 = t942 * t24
      t135 = t304 * t19 * s34
      t138 = t1015 * (-t120 + t61)
      t147 = t589 * propW34
      t24 = t147 * (t1023 * (t1017 * (t22 * (t1015 * (t814 * t824 - t926
     # * (-epinv * t824 + t678 + t816 + t923 + t924) - t24 * t776) + t18
     #8 * t797) - t135 * (s12 * (t464 - t617 - t414 - t407 - t428 - s16 
     #+ s34 - t65 - s12) - t13 * (-t429 + t393 - s34 - s45)) * t1015 + t
     #24 * t1015 * t134 * t9) + t138 * (s26 * t62 + (t65 + t64) * s45 + 
     #(t65 + s25 + s16 - s34 + s45 - s13 + s46) * s46 + (t67 + s25 - s12
     # + s16 - s34 + t68) * s56 - t69 + t70 + t73) * t11) + t346 * t118 
     #* t785 * t590 * propW16)
      t70 = t2 * t1017
      t153 = t2 * t794
      t162 = t153 * t13
      t170 = t50 * s34
      t97 = (-t535 * t944 + t955 * (s12 * (-t28 * (t954 - t179) - t66 * 
     #t366 - s13 * t957 + s16 * (s25 + s26 - s34 + s45 + s46 + s56) + t1
     #11 - t146 + t518 + t71 + t958 - t865) + t188 * t27) + t977 * (s12 
     #* (s12 * (-t404 + t744) - t17 * t97 + t28 * (s26 * t980 + s34 * (t
     #414 + t407 + s45 + s56) - t192 - t360) - t66 * t404 * t148 - t688 
     #+ epinv * (t518 - t146 + t958 + t111 + t71) - s13 * (-t28 * t979 -
     # t234 + t978) + s16 * (epinv * (s26 - s34 + s45 + s46 + s56) + t64
     #0 - t516)) + t396 * t27)) * t1015 * t9
      t177 = t2 * t779
      t187 = t2 * t792
      t188 = t1017 * t1013
      t211 = t346 * (t1013 * (t13 * (s34 * t703 * t777 + t177 * t48) + t
     #817 * t823 + t820 * t825 + t927 * (epinv * t820 + s25 * t833 + s45
     # * (s13 * t928 - t66 * s25 * (t254 + s25) - t218 - t934) + s46 * (
     #-t28 * (-t821 + t111 + t182) - (s16 + t929) * s25 - (t28 * t918 + 
     #t882) * s26 - (-t930 + s25) * s34 - s45 * (t28 * t931 + s13 - s45 
     #+ t194 + t490) + t935) + s56 * (s56 * (-t918 - t588 - t882 - s56) 
     #- t28 * s25 * (t566 + s25) - t493 - s45 * (t127 - s13 + t563 + t53
     #7 + s45) - s46 * (t855 + t345 + t933 - t483 + t200 + s46) + t821 -
     # t919) - t493 * t568) + t936 * t940) - t187 * t13) * t22
      t92 = t11 * (t190 * t228 + t214 * t489) + t19 * (t228 * t229 - t29
     #5 * (s12 * (s12 * (-t92 + t28 * (t174 + t142) - t66 * (-t241 - t24
     #4 + t436) - s13 * (t431 - s13) - s26 * (s25 + s56 - s46) - t111 - 
     #t262) + t124 * s34 * (s46 * (s25 + s46) + t192) + t28 * (s25 * t47
     #3 - s26 * (t149 + t113 + t71) + s45 * (t338 + t93 - t113) - t474) 
     #- t66 * (t445 + t472 + t380 + t471) + 7.0_dp * t264 - s13 * (s16 *
     # t460 - t17 * (-t146 - t197) + t66 * (t366 - t256 + t111) - t149 -
     # t167 + t196 - t222 - t468 - t185) - s25 * (-t475 - t210 + t167 + 
     #t113 + t71) - t114 + t275 - t385 - t386 - t72 - t476 - t477 - t478
     # - t480 + t481 + t312 * t142) + t13 * (-s16 * (-s16 * t28 * t484 +
     # t486) + t488)) + t214 * t225 + t297 * t390 + t650 * t207 - t456 *
     # (s12 * t410 + t284 * t397) - t308 * t191 - t195 * t129)
      t129 = t796 * t12
      t185 = t19 * t873
      t191 = t1015 * s34
      t48 = t274 * (t1013 * (t1017 * (t13 * (-t1015 * t22 * t48 * t877 +
     # t11 * t129) - t1024 * t867 * t2 * t19) + t191 * (t1017 * (t11 * t
     #741 * t872 - t185 * (s12 * (-t28 * (-epinv * t204 - s25 - t414 + t
     #464 - t617) + t127 + t727 + s12) - t13 * t875)) - t301 * t39 * t11
     # * t26) + t1015 * (t1017 * t92 - t11 * (t320 * t881 + t815 * t880)
     # * t2) * t189) - t22 * t787 + t1015 * (-t1 * t2 + t774) * t9)
      t92 = t191 * t9
      t82 = t274 * (t1013 * t11 * (t191 * t30 + t55) * t1016 * t26 + t18
     #8 * (-t92 * t702 * (epinv * t741 - t13 * t703) - t19 * t78 * (t13 
     #* t81 + t82) + t11 * t74 * (t12 * (t77 + s12) + t75)) + t1015 * (t
     #1013 * t11 * t57 + t1017 * t22 * t56) * t13 + t70 * t212 * (t11 * 
     #t326 * t762 + t170 * t52 * t560 + t465 * t535 * t9))
      t48 = t82 + t48 + t274 * (t211 + t188 * (t19 * (t1015 * (t674 * (e
     #pinv * t943 - t16 * t179 * t89 + t17 * (s16 * (-t146 + t279) - s26
     # * t479 + s56 * (s46 * t983 + t172)) - t28 * (-s12 * t113 - s16 * 
     #t667 - s26 * t667 + t111 * t155 - t72 + s12 * (t210 - t179 - t71) 
     #+ s13 * (s46 * (s12 + s26 - s34 + s56) + s56 * (s12 + s26 + s56) +
     # t244 + t370 + t91) + s16 * (t675 - t71) + s25 * (t146 * t66 + t28
     # * t479 - s12 * t119 - s26 * (t516 + s34) - t113 - t167 + t291 - t
     #370 - t91) + s56 * (t302 - t167) - t186 - t197 * t89 + t93 * s46))
     # + t941 * t943 - t170 * t782 * t1014) + t1024 * t995) + t97) + t10
     #13 * (t1015 * (t189 * (t320 * t698 + t815 * t837 - t70 * (t214 * t
     #790 + t228 * t788)) - t162 - t326 * t968 * t1017) - t301 * t791 * 
     #t26) * t11)
      t70 = 0.2e1_dp / 0.3e1_dp
      t82 = 0.2e1_dp / 0.9e1_dp
      t97 = 0.1e1_dp / 0.9e1_dp
      t18 = -(t85 + t18 + t274 * (t1017 * (t1013 * (t1015 * (-t265 * t27
     #8 * t189 + t22 * t915 * (t816 + t924 + t923 + t678 - t925 * (s26 *
     # t838 - (t65 + t565) * s45 - s46 * (s25 + s16 + s45 - s13 + t65 + 
     #s46) + (-s25 + s12 - s16 - t67 - t68) * s56 + t227 - t241 + t592 -
     # t73))) - t49 * t798 * t11) + t108 + t175) - t19 * t786 + t1013 * 
     #(t26 * (-t266 + t259) * t9 + t1015 * (-t26 * t961 - t594 * t809 + 
     #t799 * t808) * t11) * t1016 + t270 * t22 * t1015 * t1013) + t630 *
     # (t1017 * t514 + t1015 * (-t9 * t989 * (t26 * t447 + t990) + t11 *
     # t985 * (s34 * (epinv * t26 + s34 - t986) + t804 + t805 - t987 + t
     #988)) * t1016 - t560 * t1015 * t1017 * t1021 * t50 * s34)) * t97 +
     # 0.4e1_dp / 0.3e1_dp * t96 * t590 * t88 + t24 * t70 + t48 * t82
      t24 = s25 - s26 + s56 - s13 + s46 - t822 + t133
      t26 = s26 + s45 - s13 + s46
      t48 = s16 + t579
      t49 = s12 + s13 + s16
      t50 = s25 + s26 + s34
      t67 = s25 - s12 - s16 + s26 + s34 - s13 + s46 + t59
      t73 = s25 * t48
      t85 = t49 * s26
      t108 = (t569 + s16) * s34
      t155 = t50 * s46
      t170 = (-t67 - s56) * s56
      t175 = t569 * s45
      t196 = t568 * s25
      t204 = (t26 + t327 + s56) * s56
      t207 = t65 + t882
      t210 = (t577 + t563) * s25
      t211 = (s12 + s13 + s16 + s25) * s26
      t214 = (t882 + s26) * s34
      t228 = t207 * s45
      t236 = (-s12 - s16 + s25 + s34 + s45) * s46
      t238 = (s16 + t851) * s25
      t239 = s25 + s16 - s45
      t245 = s12 + s26
      t246 = t28 * t245
      t247 = s25 - s12 - s16 - s13 + s46 + t59 + t330
      t249 = (t63 + t844 + t127) * s34
      t256 = t28 * (t249 + t987)
      t264 = t17 * t149
      t275 = t565 * s26
      t278 = t565 * s16
      t284 = t579 * s16
      t286 = t136 - t822
      t287 = s25 + s26 - s13 + s46
      t292 = -t822 + t133 + t287
      t207 = t207 * s34
      t320 = t65 * t49
      t323 = t17 * s13
      t326 = s12 + s16 + s26 - s34
      t328 = t28 * t326
      t337 = t327 - s16 + s26 + s34 + s45
      t287 = t28 * t287 + t133 - t819
      t348 = (t563 + s25 + s12 + s13) * s26
      t351 = (t65 + s16 + t490) * s34
      t353 = (t127 + s16 + t490) * s45
      t360 = t28 * (t73 - t71)
      t366 = (s16 + t323) * s25
      t390 = (t101 + t490) * s45
      t395 = (t328 + t882 - s45) * s46
      t396 = t28 * t364
      t410 = s45 + s13
      t431 = s34 * t217
      t436 = -s45 - s56 + s13 - s46
      t450 = s25 + s26 + s46 + s56
      t458 = s12 * t450
      t460 = s13 * s26
      t468 = s25 - s16 + s56 - s13
      t473 = s12 - s34 - s46
      t474 = -t28 * t473 + t133 + t468
      t476 = t66 * t819
      t473 = t28 * t468 + t133 - t473
      t478 = t734 - s12 - s16 + s26 + s34 + t133
      t479 = t66 * t63
      t480 = s12 + s16 + s26
      t481 = -t142 * t66 + s16 * t480 + (-t65 - t479) * s34 - (t330 + s1
     #6) * s45 + (t800 - t330) * s46
      t482 = s25 + s12
      t488 = -s25 - s12 - s16
      t514 = t28 * t149
      t516 = s25 - s16 + s26
      t535 = t819 * s46
      t67 = (-t67 - s56) * s56
      t541 = s12 + s25 + s26
      t562 = t79 * t232
      t570 = (t45 + s56) * s56
      t572 = -s12 * t885 + s46 * t541 + t562 + t570
      t574 = s34 * (-s34 + t480)
      t326 = s46 * t326 + t570 + t574 + t987
      t570 = s12 + s25 + s26 - s34
      t575 = s25 * t101
      t576 = t364 * s45
      t581 = t570 * s46
      t587 = (-s13 + s25 + s26 - s34 + s45 + s46 + s56) * s56
      t588 = s12 * s16
      t589 = t28 * t13
      t590 = t28 * t12
      t593 = s12 - s34 - s45 - s46
      t594 = t79 * t960
      t48 = (-s46 + t48) * s25
      t595 = (-t488 + s26) * s26
      t596 = t48 + t595 - t594 + t988
      t574 = -t28 * t574 + t197 - t802 + t804 + t805 + t988
      t597 = epinv * t593
      t599 = t404 * t593
      t611 = t63 * s25
      t610 = t179 * t28 - s13 * t51 - s25 * (-t123 - t744 - s25) - s26 *
     # t300 - s34 * t610 - t146 + t149 + t71
      t625 = -t194 + t886
      t635 = t194 * t810
      t645 = s13 * t565
      t660 = s25 * (-t884 - t490) + t66 * (t645 - t144) + s12 * (-t601 +
     # s13) + t359
      t667 = s13 * (s16 + t327)
      t678 = s25 * t232
      t681 = t17 * t144
      t682 = t842 + t127 + t579
      t683 = s25 - s12 - s16
      t685 = (t605 + t851) * s16
      t688 = (t896 + t577 + s13 + t327) * s26
      t692 = s25 * (-s25 + t49)
      t699 = t124 * t144
      t703 = 7.0_dp * t101
      t722 = -t124 * t402
      t583 = t28 * (s25 * t583 + t75)
      t727 = t66 * (s16 + s45 - s13)
      t741 = t124 * t205
      t800 = s13 * (t103 + s25 - s34 - s45 + s46)
      t802 = t101 * s45
      t808 = 9.0_dp * s26 * t893 + 9.0_dp * t539
      t815 = t16 * (t111 + t184)
      t816 = t124 * s16
      t817 = t896 + s13
      t819 = 9.0_dp * s12
      t820 = (t816 + t606) * s16
      t824 = t66 * (t828 + t149)
      t833 = t3 * s12
      t846 = t16 * t101
      t855 = t124 * t801
      t901 = s34 * (t810 + t602) + s12 * (-t810 - t233) + s16 * (-t363 -
     # t310 - t563 + t851) - s26 * (t251 + t537)
      t911 = (s16 * (t819 + t817) + t609) * s25
      t912 = t592 * t578
      t918 = t851 * t144
      t824 = -t918 + s26 * (s26 * (t582 + t884 + s25 + t323 + s26) + s12
     # * (t345 + s13 + t701) + s16 * (t896 + t357) + t841 + t609) + s34 
     #* t901 + s45 * (s45 * (t569 + s34 - t194) - t17 * (t591 + t111) + 
     #t699 - (t933 + s13 + t490) * s26 - s34 * (t66 * t826 - s34 - t345 
     #+ t838) + t685) + s46 * (s46 * (-t330 + t844) + t824 + s26 * (t884
     # + s25 + t606 + t886 + s26) + (-t582 - t392 - t882 + s13 - t701) *
     # s34 - (-t1009 * t28 + t801) * s45 - t608 + t820 + t609) + s56 * (
     #s56 * (-t59 + t476 + s26) + s16 * (t896 + t819) + t28 * (t167 + t1
     #11 + t75) + t824 + (s25 + t833 + t895 + t886) * s26 + (-t101 * t30
     #6 + s13 - t582 - t882) * s34 - (t577 + t327 + t846 - t372) * s45 +
     # (-t59 + t855 + s26 - t361) * s46 - t608) + t911 - t912
      t901 = (-t345 - t842 - s12 - s26) * s26
      t919 = (t17 * t205 - t28 * t890 + t251 + t577 + t816) * s56
      t923 = (-t682 + s34) * s34
      t924 = t564 * t182
      t920 = t920 * t113
      t928 = t359 - t75
      t931 = t579 * s13
      t940 = t931 + t928
      t943 = -t577 + s13
      t953 = 7.0_dp * t144
      t954 = t66 * t482
      t958 = 13.0_dp * s16
      t960 = t28 * t110
      t962 = 7.0_dp * t801
      t964 = t17 * t364
      t967 = (s13 + t606) * s16
      t976 = (t577 + t842 - t886) * s25
      t983 = (-t964 - t816 - t579 + s45) * s45
      t1003 = t149 + t111
      t1004 = t182 + t75
      t1024 = t66 * t1003
      t1026 = t28 * t482
      t1029 = t12 + t1026
      t1030 = t16 * s12
      t1031 = t1030 + s13
      t1032 = t1031 * s16
      t1033 = (-t816 - t127 - s12 + t372 + s45) * s45
      t1034 = t66 * t149
      t1035 = t816 + s13
      t1036 = (s13 + t833) * s16
      t1037 = (t1030 + t884 + s13) * s25
      t1038 = t953 + t1036 + t1037 - t359 + t75
      t1039 = t306 * s25
      t1040 = t66 * (-t460 + t144)
      t1041 = t930 + s16 + t490
      t1042 = (t896 + t605) * s25
      t1043 = t810 * s16
      t49 = -t28 * (t1043 - t149 + t111) - (t345 + t892 + s12) * s26 - (
     #t28 * t49 + s26 - t537) * s34 - t1042
      t1044 = s25 + s12 - s34
      t1045 = t515 - t330 + t579
      t1046 = t17 * t482
      t1047 = (t954 + t892 + t886 + s26) * s26
      t1048 = (s26 + t960) * s45
      t1049 = (s16 * (t1035 + t833) + t609) * s25
      t1050 = (t1024 + t681 + (t962 + s13 + t327) * s26 - (t582 + t392 +
     # t884 + t490 - t886) * s34 + s45 * (-t1044 * t28 + s45 - t200 - t8
     #16) + t210 - t608 + t967 + t609) * s46
      t1051 = (t28 * (s46 * t1045 + t111) + t1034 + t681 + (t582 + t892 
     #+ s13 + t327) * s26 - (t545 + t884 + t1046 - t886) * s34 + s45 * (
     #-t856 - t816 - t127 - s12 + s45) + t1032 + t210 + t585) * s56
      t1052 = t1045 * t113
      t1045 = t1045 * t71
      t1053 = t512 + s16 * (-s16 * t943 - t940) + s26 * (t1047 + t1038) 
     #+ s34 * (s34 * t1041 + s12 * (-t127 - t845 - t537) + s16 * (-t1039
     # + t323) - t377 - 9.0_dp * t275 - t1040) + s45 * (t49 + t1048) + t
     #1049 + t1050 + t1051 + t1052 + t1045
      t1054 = (t896 + t605 + t507) * s34
      t1055 = (s16 * t943 + t940) * s16
      t1056 = t17 * t182
      t1057 = (t582 + s13) * s16
      t1058 = t28 * t182
      t1059 = s25 - s12
      t1060 = t566 * s16
      t1061 = t17 * t801
      t1062 = (t933 - t886) * s25
      t578 = (t578 + s16) * s16
      t1063 = t139 + t1026
      t1064 = t582 + t886
      t1065 = t17 * t101
      t1066 = s13 + t579
      t1067 = t66 * t144
      t1068 = t65 + t345 - s13
      t1069 = s13 * (s26 + t563)
      t1070 = t884 * s25
      t1071 = (t842 + t127 + t288) * s26
      t811 = -t811 - t882
      t515 = t515 + t827
      t1072 = (t954 + t842 + t886 + s26) * s26
      t1073 = t1066 * s16
      t1044 = (-t28 * ((t103 + t1044) * s45 - t149 - t75) - t66 * t402 +
     # (t892 + s13 + t327 + t606) * s26 + (-t251 - t816 - t882 + s13 - t
     #605) * s34 + t1032 + t167 + t210 - t608) * s46
      t1074 = (t28 * (s46 * t515 + t111 + t149) + t1067 + (t855 + s13 + 
     #t327) * s26 + (-t816 - t954 + s13 - t507) * s34 + s45 * (-t856 - t
     #842 - t127 - s12 + s45) + t1057 + t210 + t585) * s56
      t1075 = (-t144 - t1073) * s13
      t1037 = (t1072 + t1067 + t1057 + t1037 - t359 + t75) * s26
      t1076 = t515 * t113
      t515 = t515 * t71
      t1057 = -t923 - t1058 + (s16 + t886) * s25 - (t840 + t884 + t605) 
     #* s26 - t75 - t1057 - t699 - t848
      t573 = (t573 - t327) * s45
      t839 = t573 + t28 * (s25 * (t839 - s25) + s34 * (t65 - t1059 + t19
     #4 - s34) - t111) - t681 - (t896 - t882 + s12) * s26 - t1060
      t820 = -s13 * t66 * t885 + s25 * t28 * t801 + t1056 + t1072 + t820
     # - t928
      t1072 = t924 + t145
      t600 = t28 * t1072 - epinv * (s34 * (-s12 * t1068 - s34 * t811 + t
     #1069 - t1070 - t1071) + s45 * (-s45 * t604 + t28 * (t887 - t111) -
     # (t345 + s12 + t194) * s26 - s34 * (t838 - t345 + t563 - s34) - t1
     #042 + t144) + t1049 + t1044 + t1074 + t1075 + t1037 + t1076 + t515
     #) + s16 * (t586 - t847 + t75) - (s13 * t600 - s16 * t849) * s25 + 
     #s26 * t820 + s34 * t1057 + s45 * t839 + (t28 * (t149 + t578 + t182
     # + t75) + t848 + (t842 + s13 + t606 + t537) * s26 - (t201 + t1061 
     #+ t327) * s34 + s45 * (-t816 + t330 - t200 - t579 + s45) - t608 + 
     #t1062) * s46 + (t28 * (s46 * t1063 + t111 + t144 + t149 + t182) + 
     #s16 * t1064 + (t582 + s13 + t194 + t537) * s26 - (t1065 + t577 + t
     #327) * s34 + s45 * (-t364 * t66 - s12 + s45 - t816) + t1062 + t585
     #) * s56 + t113 * t1063 + t71 * t1063
      t820 = s12 - t886
      t839 = t820 * s16
      t847 = t851 * s12
      t1057 = -t847 - t928
      t1062 = 9.0_dp * t144
      t1063 = s12 * (t101 * t3 + t63)
      t1077 = 16.0_dp * s26
      t1078 = -s25 + s26
      t1079 = t306 * s16
      t1080 = 9.0_dp * t101
      t1081 = (-s13 + t194 + t579 + s25) * s25
      t1082 = (t892 + t833) * s16
      t1083 = t822 + s25 + s26 - s45
      t1084 = (t892 + t606) * s16
      t1085 = (t63 + t703 + t605) * s34
      t1086 = (t577 + t884 + s25 + t886 + s26) * s26
      t1087 = t1057 * s16
      t1088 = (-t801 * (s13 - t194) + t591) * s25
      t1089 = t566 * t66
      t1083 = -t28 * (t150 - t145) + s26 * (t1086 + s25 * (t842 + t327) 
     #+ s12 * (s25 + t701 + s12) + s13 * (-t565 - s13) + t1062) + s34 * 
     #(-s16 * (t1077 + t882) - t28 * (t168 + t75) + 7.0_dp * t402 + s13 
     #* (t232 + t194) - t182 - t1063 + t1085) + s45 * (-s45 * (s16 + s25
     # - s26 + s34) + s34 * (-t330 + t507) + t124 * s16 * t154 - t28 * s
     #26 * t1078 + s12 * (t330 + s25 + s16 - s26) + s13 * (t319 + t194) 
     #- t182) + (t1082 + t264 + t609 + (t1079 + t882 + s13 + t606) * s26
     # + (-t1080 - t1030 + s13 - t327) * s34 - t608 + t983 + t1081 + t84
     #8) * s46 + (t28 * (s46 * t1083 + t111) + t264 + (t582 + t882 + s13
     # + t701) * s26 + (-t582 + s13 - t496 - t701 - t327) * s34 + t1033 
     #+ t1081 + t585 + t1084) * s56 + t113 * t1083 + t71 * t1083 + t1087
     # + t1088 + t1089 * t144
      t1090 = s13 - t579
      t1091 = t306 * s12
      t1092 = (-t885 + s13) * s13
      t1093 = 12.0_dp * s25
      t1094 = s12 + t886
      t1095 = (s16 + t1094) * s16
      t1096 = s25 - s26 + s34
      t1097 = 7.0_dp * t93
      t1098 = t28 * (t149 + t75)
      t1099 = t17 * t566
      t1100 = t17 * t810
      t1101 = s34 - s13
      t1102 = t28 * t1101
      t1103 = (t582 + t884) * s25
      t1104 = (t17 * t565 - t361 + t545 + t582 - t886) * s46
      t1105 = t28 * (t149 - t167 + t113)
      t1106 = t66 * (t1060 - t608)
      t1107 = s25 + t194
      t1108 = (s13 + t1107) * s13
      t1109 = 16.0_dp * s16
      t1110 = t591 - s12 * (-t842 + s13) - (-t816 + s13) * s16
      t1111 = t144 * t566 + t112
      t209 = t28 * (t209 - t150 + t145) + t66 * t1111 + s25 * t1110 + s2
     #6 * (s26 * (t345 + t1079 + t605) + s25 * (t1079 + t327) + t1062 - 
     #s12 * (-s12 - t882 + s13 - t701) - t1108) + s34 * (t1085 - s26 * (
     #t1109 + t496) - t17 * t575 - t609 - t953 + s13 * (t65 + s25 + t194
     #) - t1063 - t182) + s45 * (-s45 * t1096 - t66 * (t149 - t592 + t11
     #1) - t888 - s25 * (s25 - t857) + (-t566 - t896) * s26 + (t577 - s1
     #3 + t846) * s34) + (t28 * (t804 + t149 + t111) - s25 * (-t892 - t5
     #77 + s13 - s25) + (t842 + t288 + t579) * s26 + (-t601 + s13 - t200
     # - t579) * s34 + (-t563 - s25 + s12 - s26 + s34) * s45 - t608) * s
     #46 + s56 * (s56 * (t250 + t328) + t377 + t264 + t1084 + s25 * (t84
     #2 - s13 + t579 + s25) + (t345 + t1030 + t895) * s26 + (-t582 - t10
     #80 + s13 - t327) * s34 + s45 * (-t964 - t892 - s12 + s45) + t395 +
     # t585) + t1087
      t328 = (t17 * (s16 - s34 + s46) + t251 + t577 + t850) * s56
      t570 = t570 * t113
      t964 = t28 * (t149 * t883 - t570 - t72) + epinv * t209 + s16 * (s1
     #6 * t1090 + t940) - (s16 * (t1091 + s13 + t701) + t609) * s25 + s2
     #6 * (s26 * (-t363 - s12 - t563 - s26) + s12 * (-t363 - t842 + s13)
     # - s16 * (s16 + 18.0_dp * s25) + t1092) + (s12 * (s16 + t537 - t88
     #6) - t101 * (-s26 + t323 - t1093)) * s34 + s45 * (s45 * (-t1096 * 
     #t28 + s16) + t28 * s34 * (t1078 * t66 - s34 + t563 + t838) - t132 
     #+ (t701 + t605) * s25 - (t577 + t842 + s13 - t357) * s26 - t1095) 
     #+ (s25 * (-t545 - t884 + t361) - t124 * t218 + t17 * (t661 + t174)
     # + s12 * (-t251 - t345 - t892 + t410 + t372) - s13 * (s16 + t483) 
     #+ s16 * t670 + t1097 - t1098) * s46 + s56 * (-t328 - t1106 - t1105
     # - (t363 + t892 + t1099) * s26 + (t251 + t1100 + t537 + t579) * s3
     #4 - (t1078 * t17 - t1102 + t933) * s45 - t144 - t1103 - t1104 - t8
     #89)
      t1062 = 11.0_dp * t144
      t1063 = s25 - s16 - s26
      t1085 = (t844 - s13 + s25) * s25
      t1096 = t579 + s25 + s16 + s26 - s34 - s45
      t1110 = t127 + t845 + t490
      t1111 = t124 * (-t592 + t111)
      t1112 = (t1079 + t537) * s26
      t930 = -t930 - t882 - s16
      t1113 = t28 * t838
      t1114 = -t149 + t111
      t1115 = (t816 + t579) * s25
      t1116 = t1114 * t28 - (s12 + t851) * s16 + (t816 + s12 + t327) * s
     #26 + (t1113 + s26 - t490 + t194) * s34 + t1115
      t1117 = -t516 - t483
      t827 = -t28 * (t828 + t75) - t1024 - t1082 - (t840 + t1079 + t606)
     # * s26 + (t582 + t392 + t882 + t701 - t886) * s34 - s45 * (s45 - t
     #816 - s25 - t200 - t827) + t608
      t1082 = t143 + t330 - t844
      t1084 = -t28 * (-s46 * t1082 + t111 + t828) - t1034 - t1084 - (t84
     #0 + t582 + t701) * s26 + (t1080 + t882 + t605 - t886) * s34 - s45 
     #* (-t816 - t482 - t127 + t483 + s45) - t585
      t1035 = (s16 * t1030 + s16 * t1035 + t75) * s25
      t1107 = (t1086 + s13 * t891 + s16 * (t701 + t537) + t75 + t577 * t
     #1107) * s26
      t1089 = (-t1089 - t563) * t144
      t682 = t28 * t682 * s34
      t1118 = s45 * (-t483 - t1063) - t28 * (t1060 + t111) - t264 - t699
     # - (s12 + t701 - t327) * s26 + (t251 + t1079 - s25 + t605) * s34 +
     # t692
      t1058 = t1086 + t1062 + t1058 + (11.0_dp * s12 + s13) * s16 + (t56
     #6 + t563) * s25 - t359 + t75
      t591 = epinv * (t1089 - t1057 * s16 - s34 * (-s12 * t1110 - s34 * 
     #t930 - t1040 - t1070 - t1111 - t1112) - s45 * (-s45 * t1117 - t111
     #6) + s46 * t827 + s56 * t1084 + t113 * t1082 + t71 * t1082 - t1035
     # - t1107) + s16 * (s16 * (t1030 - s13) + t28 * t585 - t359) + s25 
     #* (t591 + t284 + t1043 - t608) + s26 * t1058 + s34 * (t682 - t306 
     #* t804 - t377 - t609 - (t1109 + s13 + t606) * s26 - t611) + s45 * 
     #t1118 + (t28 * (-(s12 + t194 + t396) * s45 + t75) + t1024 + t699 +
     # (t884 + t882 + s13 + t606) * s26 - (t896 + t582 + s25 + t507) * s
     #34 + t1036 + t167 - t608 + t1085) * s46 + (t28 * (s46 * t1096 + t1
     #11) + t1034 + t699 + (t896 + t582 + t882 + s13) * s26 - (t896 + t2
     #01 + s25 + t605) * s34 + s45 * (-t892 - t127 - s12 + t372 + s45) +
     # t1085 + t585 + t967) * s56 + t113 * t1096 + t71 * t1096
      t1036 = t28 * t150
      t1043 = s12 * (t896 + t507)
      t1057 = 12.0_dp * s16
      t1058 = t16 * t144
      t1085 = (t201 + t1057) * s26
      t1086 = t66 * s16 * (t838 + s16)
      t1096 = s16 + s34
      t1118 = t842 * s12
      t1119 = (t577 - t886) * s16
      t1120 = t149 - t608
      t1121 = t17 * t232
      t1122 = t1121 + s16 - t886
      t1123 = (-t816 - t127) * s13
      t1124 = (t251 + t895) * s26
      t1125 = (t884 + t507) * s25
      t1126 = -t602 - s16 + t886
      t1127 = t28 * t885
      t831 = -t232 + t831
      t1128 = s16 + s25 + s26 - s34
      t1129 = t66 * (t149 + t144)
      t1130 = -t327 * t801 - t1129 - s26 * (t840 + t816 + t605 + s26) + 
     #(t392 + t884 + t882 + t605 - t886) * s34 + t1128 * s45 - t585 - t1
     #118
      t933 = t143 - t933 + t361
      t1131 = (-t482 - t842 - t127 + t483 + s45) * s45
      t1132 = (t840 + t582 + t884) * s26
      t1133 = (t392 + t882 + t701 + t605 - t886) * s34
      t1134 = t28 * (t828 + t111)
      t1135 = (-t816 - t1030) * s16
      t1136 = (t801 - t330) * s46
      t1137 = (s26 * (-t892 - t577 - s25 - t886 - s26) - t577 * (s25 + t
     #563) - s16 * (t816 + t537) + t1092 - t75) * s26
      t1026 = s56 * (-t28 * (s25 + s12 + s46) - t194 - t200 + t890 - s56
     #) - t28 * t1120 - t1067 - t132 - s25 * (s25 + t816 + t577 - s13) -
     # (t896 + t577 + t537 - t886) * s26 - (-t1065 - t1026 + s13) * s34 
     #- s45 * (-s25 + s16 + s34 + s13 + t246 - s45) - s46 * (s46 + t201 
     #+ t954 + t842 - s13 - t483) - t1119
      t1065 = s45 * (t1096 * t28 + t232) - t1086 + t264 - s25 * (t838 - 
     #s25) - (t840 + t194 + t579) * s26 - (t201 + t892 - s25 + t605) * s
     #34 - t853
      t1026 = epinv * (s34 * (s12 * t1122 + s34 * t1126 + t1123 + t1124 
     #+ t1125 + t888) + s45 * (s45 * t831 + t28 * (s34 * (t838 - t1127 -
     # s34) + t111) - t1086 + (s12 + t567) * s26 + t1115) + s46 * (-t113
     #6 + t1130) + s56 * (s46 * t933 + s56 * t1082 - t1034 - t1131 - t11
     #32 + t1133 - t1134 + t1135 - t585) - t1035 + t912 + t1137 + t918) 
     #- s26 * (s26 * (-t882 - t842 - s12 - s26) + s25 * (-t892 - t327) -
     # t1067 + s12 * (s13 - t327 - t194) + t1108) - s34 * (s25 * t916 + 
     #t1043 - t1054 + t1058 + t1085 + t609 - t645) - s45 * t1065 - s46 *
     # (-t921 - t66 * (t149 + t111 + t144) - (t892 + t1046) * s26 + (s25
     # + t846 + t605) * s34 + (-t330 + s16) * s45 - t1081 - t585 - t1118
     #) - s56 * t1026 + t1088 - t912 - t918 - t1036
      t50 = t28 * t50
      t1046 = t17 * t522
      t1065 = -s16 - s13
      t1081 = t1065 * s16
      t580 = -s25 * t580
      t1086 = t28 * t167
      t1057 = s34 * t1122 - s12 * (t201 + t537 - t886) - s16 * (t563 - t
     #323) - s26 * (t1057 + t357 - t851) - 9.0_dp * t539 - 9.0_dp * t111
      t344 = s16 * (-t838 ** 2.0_dp - t1073) + (s16 * (s13 + t194 + t606
     #) + t609) * s25 + s26 * (s26 * (s26 + t577 + t842 + s25 + t886) + 
     #t1067 + (t582 + t842 + s13) * s25 - t359 + t586 + t75) + s34 * t10
     #57 + s45 * (-s45 * (s16 - t50) - s16 * (t201 - s16) - t16 * s25 * 
     #t64 + t28 * (-s12 * (t254 + t327) + s34 * (s26 + s34 - s13) + t592
     #) - t66 * s16 * t806 - t344) + s46 * (-s46 * (t232 + t330) + t1034
     # - s26 * (t62 + s26) - (t896 + t954 + t392 - t886) * s34 - (-s12 +
     # s16 + s25 + s26 - s34) * s45 - t608 - t1081 - t580) + s56 * (s56 
     #* (-t28 * (s26 + s34 + s46) + s12 + s13 - s16 - s25 - s45 - s56) +
     # t1129 - (t1065 * t66 + s25 - t605) * s26 - (t392 + t884 + t882 - 
     #s13 + t605) * s34 - (t840 + t201 - t330 + t892 + t579) * s45 + s46
     # * (-t126 - t127 - t327 + s12 - s16 - s45 + s13 - s46) - t580 + t6
     #85 + t75 + t1086)
      t580 = 17.0_dp * s16
      t586 = 14.0_dp * s12
      t954 = t3 * t111
      t1057 = t124 * t149
      t1065 = 14.0_dp * s25
      t1073 = t28 * t165
      t1088 = t391 * s16
      t1108 = -s16 - t200 - t537 + t886
      t1118 = s16 * (-t124 * t838 - t194) + s26 * (t582 + s13 + t357 + t
     #895 + t496) + t28 * (s34 * (t28 * t826 - s34 - t345 + t838) + t104
     #2)
      t1129 = -t1121 - t483 + t194
      t1138 = (t586 + s13) * s16
      t1139 = (t582 + t563) * s25
      t1140 = t28 * (t167 - t608)
      t1080 = t17 * (t111 + t75) + t1057 + t1062 + (t586 + t580 + t851) 
     #* s26 - (t345 + t819 + t1088 + t369 - t886) * s34 - (-t126 + t1080
     # + t882 + t605) * s45 + t1138 + t1139 + t1140
      t1141 = t59 + t126 - t1061 - t1078
      t1142 = s16 - s46
      t1143 = t319 * s34
      t1144 = t17 * t167
      t568 = -t1142 * t28 + t17 * t806 + t133 + t568
      t841 = (s12 * (t577 + t1088 + t1039 + s13) + t28 * (s16 * (t896 + 
     #t537) + t841)) * s26
      t1064 = (s13 * t1064 - s16 * (s12 - t323)) * s16
      t1145 = t28 * (t111 * (t896 + s25 + t605 + t851 + s26) + t1049 + t
     #145)
      t233 = s56 * (t233 - t410 + t605 + t194 - t1073 + s56) + t1120 * t
     #66 - t16 * t402 + t1056 + (t701 - t323 + t606) * s25 + (t1065 - s1
     #3 + t833 + t895) * s26 + (-t703 - t345 + s13 - t605) * s34 + (t330
     # - t842 - s25 + s12 - s26 + s13) * s45 + s46 * (t3 * t541 - s13 + 
     #s46 - t126 - t133 + t892) + t75 + t884 * s12
      t365 = s46 * (t1128 * t66 - t59 + t605) + t1004 * t17 + t144 * t30
     #6 + t28 * (t967 + t167 - t608) + (t855 - t323) * s25 + (t1039 + t5
     #86 + s13 + t580) * s26 - (13.0_dp * t101 + t882 + t819) * s34 + (-
     #t365 - t1079 + s25 - t605) * s45 + t954 + t1057
      t402 = t682 - s16 * (t819 + t1079) - t1056 - t609 - t889 - (s13 + 
     #t580 + t327 + t606) * s26 + t366
      t541 = -s16 * t28 * t1090 + s13 * (-t605 - t886) + t75
      t682 = -s45 * (s16 + t483 + t490) - t17 * t770 + t28 * s16 * (-s16
     # + t1094) + (t17 * t838 + t896) * s25 + (t838 - t816 + t537) * s26
     # + (-t1059 * t17 + t507 + t701) * s34
      t206 = epinv * (s34 * (s34 * t1108 + s16 * (16.0_dp * s25 + t194) 
     #+ s26 * (t392 + t580 + t1093) + s12 * (t251 + t1039 + s16 - t886) 
     #- s13 * (t816 + s26)) + s45 * (s45 * t1129 + t1118) + s46 * (s46 *
     # t1141 - t1080) + s56 * (s56 * (t568 + s56) + s45 * (t958 - t469) 
     #+ s46 * (t260 - t892) - t124 * (-t1143 - t279) - t28 * (s25 * t114
     #2 + t111) + t391 * (s34 * t101 - t676) - t66 * (-t174 + t75) + s12
     # * (-t124 * t250 + t165 * t3 + s13 - 14.0_dp * s16 - t310) - s13 *
     # (t201 - t206 + t194) + t113 - 11.0_dp * t144 + 11.0_dp * t161 - t
     #1144) - t841 + t1064 - t1145) + s16 * t541 + (s16 * (-s13 * t16 + 
     #t577) - t605 * s13) * s25 + (t16 * t852 + t182 * t3 - t28 * t359 +
     # 12.0_dp * t144 + s12 * (s12 - s13 + t490 + t895)) * s26 + s34 * t
     #402 + s45 * t682 + s46 * t365 + s56 * t233
      t233 = (t816 + t882 + s13 + t579 + s26) * t111
      t14 = t28 * t14
      t359 = t77 - s12 - s16 - s26 + s34 + s45
      t365 = t41 - s25 - s26 - s45 - s56 + s13 - s46
      t402 = intHs16s25s26s34s56x1411D6eps0()
      t541 = s56 * (t65 + s25 + s16 - s34 + s45 - s13 + s46 + s56) + t81
     #2 * t89 - t592 + t813 + t828 + t832
      t676 = t1101 * s34
      t682 = -t1089 + s34 * (-s12 * t1110 - s34 * t930 - t1040 - t1070 -
     # t1111 - t1112) + s45 * (-s45 * t1117 - t1116) - s46 * t827 - s56 
     #* t1084 - t113 * t1082 - t71 * t1082 + t1035 + t1087 + t1107
      t254 = s34 * (-s12 * t1068 - s34 * t811 + t1069 - t1070 - t1071) +
     # s45 * (-s45 * t604 - s16 * (t127 + t288) - t28 * (s16 * t1101 + t
     #111) + t296 + s12 * (-t254 + t563 - t490) + t144 + t676) + t1037 +
     # t1044 + t1049 + t1074 + t1075 + t1076 + t515
      t296 = t1145 - s34 * (s34 * t1108 - t28 * (-t1103 + t608) + t1067 
     #+ t558 + (-t897 + s12) * s16 - (s13 - t580 - t606 - t1093) * s26) 
     #- s45 * (s45 * t1129 + t1118) - s46 * (s46 * t1141 - t1080) - s56 
     #* (s56 * (t568 + s56) - t66 * t75 - t1057 - t1062 - t1144 - t853 -
     # (t586 + t851) * s16 - (12.0_dp * s12 + t1088 + t897) * s26 - (-t1
     #01 * t391 + s13 - t345 - t833) * s34 + (t582 + t898 + t958 + t882 
     #+ s13 - t469) * s45 - s46 * (-t260 - t497 + t892 + s13 - t327 + t8
     #33 - s46) - t1139 + t608) - t1064 + t841
      t515 = s16 * (t1065 + s16)
      t558 = (t816 + t605) * s25
      t568 = (-s12 - t563 - t288 - s26) * s26
      t558 = s26 * (t568 - s12 * (-s13 + t563 + t357) + t1092 - t515) + 
     #s34 * (s34 * (-t882 + s13) - s12 * (-t345 + s13) - t101 * (-t363 -
     # s26 + t886)) + s45 * (t573 - t28 * (t887 + t111) - (t844 + s13 - 
     #t288) * s26 - s34 * (-t838 + t345 - t127 - t563 + s34) - t144 + t4
     #90 * t564) + s46 * (-t581 - (t605 + t194 + t537) * s26 - s34 * (-t
     #1121 - t844 + s13 + s34) + (t882 + s16) * s45 - t578 - t585 - t558
     # - t848) + s56 * (s56 * (-t28 * (s12 + s16 - s34 + s46) - t200 - t
     #63 - s56) + t931 - t132 - s16 * (-t1090 + s16) - (t577 + t842 + t5
     #37 - t886) * s26 - s34 * (-t1121 - s12 - t845 + s34) - s45 * (-s45
     # - t882 + s16 - s34 + s13 + t246) - s46 * (-t66 * t806 - s13 + s46
     # + t201 + t567) - t558) - t911 + t592 * (t1066 + s16)
      t564 = 14.0_dp * s26
      t567 = (t582 + t892) * s25
      t573 = 11.0_dp * s34
      t578 = s45 * (-s16 + t327 + t372) - t28 * ((t1113 + t843 - t490) *
     # s34 + t1115) + t264 + s16 * (t897 + s12 - s16) + (t838 - t816 - t
     #537) * s26
      t580 = s34 * (t65 + s16 - t323 + t537) + s12 * (-t843 + t323 - t35
     #7) + s13 * (t251 + t701) + s25 * (-t1109 - t564) + s26 * (-t507 - 
     #t895) - t1067
      t567 = t28 * (t1035 + t145 + t233) + s16 * (-t582 * s13 + s16 * (t
     #582 - t323) - t28 * t928) + (t28 * (s16 * (t816 + t357) - t1092) +
     # s12 * (t958 - s13 + t357 + s12)) * s26 + s34 * t580 + s45 * t578 
     #+ s46 * (s46 * (-t28 * t670 + s25 - t126 + t127 + t605) + t17 * t7
     #5 + t1057 + t1140 + t953 + t954 + (t1088 + t586 + s13 + t537) * s2
     #6 - (t819 + t958 + t369 - t323 + t288) * s34 - (t319 * t66 + t251 
     #+ t605 + t701) * s45 + t1138 + t567) + s56 * (s56 * (t201 + s25 - 
     #s45 - s13 + t235 + t1061 - t361 + s56) + t953 - t847 + (t1091 - s1
     #3) * s16 + (t958 - s13 + t833 + t537) * s26 - (-t897 + t700 + t103
     #0 + t958 + t288) * s34 + (-t816 + s12 - s26 + s34 + s13 - t882) * 
     #s45 + s46 * (t245 * t3 - s13 + s46 - t133 + t327 - t573 + t892) + 
     #t567 + t75 + t1057 + t889)
      t578 = t882 - t886
      t323 = t28 * (t570 + t72) + s16 * (-s16 * t1090 - t940) + (s12 * t
     #896 + s16 * t817 + t75) * s25 + s26 * (-t568 + s12 * (t842 - s13 +
     # t537) - t1092 + t515) + s34 * (s34 * t578 + s12 * (-s16 - t490 + 
     #t886) + t101 * (-t363 - s26 + t323)) + s45 * (s45 * (t885 - t396) 
     #- t28 * ((s12 + t194) * s25 + s34 * (t838 - t1127 + t127 - s34)) +
     # t132 + (-t345 + t577 + t842 + s13) * s26 + t1095) + (t1098 + t377
     # + s16 * (t1031 + s16) + (t816 + t357 + t606) * s26 - (t251 + t345
     # + t1061 - t886) * s34 - (t801 + t882) * s45 + t1042 - t608) * s46
     # + s56 * (t328 + t1106 + t1105 + (t892 + t1099 + t357) * s26 - (t2
     #51 + t345 + t1100 + t579) * s34 + (-t66 * t683 - t1102 + t200) * s
     #45 + t1042 + t1104 + t144 + t889)
      t328 = t822 + t127 + s25 - s45 + s56 - s13 + s46
      t363 = 12.0_dp * t148
      t515 = t2 * t762 - t968
      t568 = t295 * (-t918 - s34 * (s12 * t1122 + s34 * t1126 + t1123 + 
     #t1124 + t1125 + t888) - s45 * (s45 * t831 - s12 * (-s26 + t194 - t
     #960) - t17 * s34 * t885 + t28 * (s26 * t916 - t676) + t66 * t1081 
     #+ t345 * s16) - s46 * (-t1136 + t1130) - s56 * (s46 * t933 + s56 *
     # t1082 - t1034 - t1131 - t1132 + t1133 - t1134 + t1135 - t585) + t
     #1035 - t1137 - t912) + t297 * t824
      t249 = t1013 * (t189 * (-t323 * t698 - t489 * t682 + t558 * t837) 
     #- t944 * t481 + t955 * (-s16 * t488 - s26 * (-t482 - t563 - s26) +
     # t184 - t249 + t803 - t987 + t514) + t977 * (-t28 * (s25 * (s16 + 
     #t329) + s34 * t436 - s56 * t436 - t675) + epinv * (t28 * t431 + s1
     #6 * (t181 + t396 + s16) - s25 * t368 + s26 * t271 + s34 * t410 + t
     #111 - t172 - t93 * t66) + s12 * (-t235 + epinv * (t148 - t483) + t
     #1012)) + (t284 + (-t479 - s26) * s34 + (-t826 - t330) * s45 + (t23
     #9 - t330 + t246) * s46 + s56 * (t45 - t330 + s56) + t275 + t278) *
     # t515 - t702 * s34 * (t328 * epinv - t589)) + t774
      t329 = t995 * t22 * t1014
      t396 = t9 * t1016
      t52 = t2 * t52
      t436 = t52 - t782
      t323 = -t190 * t254 + t2 * (t254 * t788 + t323 * t881 - t558 * t88
     #0 + t682 * t790)
      t479 = t1025 * t867
      t482 = t2 * t1015
      t488 = t2 * t22
      t64 = t212 * (t11 * t323 + t22 * (-t195 * (-t28 * (t920 + t924 + t
     #72) + t918 + epinv * t824 + (s12 * t625 - t635) * s25 + s26 * (t90
     #1 + t660) + s34 * (t923 + s26 * (t884 + t127) + t28 * t678 - s12 *
     # (-t842 - t127 - s12) + t539 - t667 + t681) + s45 * (s45 * (t569 +
     # s34 + t194) - t28 * (t692 - t149) - t132 - (-t28 * t683 + s26) * 
     #s34 - t685 - t688 - t699) + (s26 * (t505 - t288) - t66 * t522 + s1
     #2 * (-t703 + t126 + t410 - t882) + s16 * (-t345 - t392 + s45 + t50
     #5) - t722 - t583) * s46 + s56 * (-t919 + t17 * (-s34 * t164 - t802
     #) + t28 * (t695 + t800 - t182 - t113) - t66 * s16 * t64 + s12 * (-
     #t727 - t741 + s34 - t200) - t518 - t808 - t815 - 11.0_dp * t168 - 
     #t896 * s46) + t912) - t225 * t682 - t229 * t254 + t308 * t600 + t4
     #56 * (t145 * t17 - t1036 + t591) + t650 * t1026)) * t189
      t1 = t274 * (t11 * (t1013 * (t13 * t796 + t74 * (t77 - t14)) + t48
     #2 * (t1013 * t465 * t481 - t1)) + t1013 * (t479 * t22 + t191 * (t5
     #60 * t436 + (t2 * t39 - t30) * t1016 * t11)) * t593 + t212 * (s34 
     #* t328 * t872 - t162) * t9 + t488 * t1015 * (-t1013 * t13 * t779 +
     # t792) + t64)
      t1 = t1 + t274 * (t1015 * (t11 * t249 + t22 * (t1013 * (t13 * t877
     # - t189 * t568 - t674 * (epinv * (-t170 - t73 - t85 + t108 + t155 
     #+ t175) + t28 * (t204 + t922 + t921 + t196)) - t823 * (t67 + t211 
     #- t214 - t228 - t236 + t210) - t825 * (t320 + (-t602 + s16) * s45 
     #- s46 * t286 + s56 * (-s56 - t292) - t207 + t210) + t927 * (s25 * 
     #(t65 + t882) + epinv * (-t320 + s46 * t286 - s56 * (-s56 - t292) -
     # t210 - t818 + t207 + t133 * t232) - t932 * s46 - s56 * (-s56 - t2
     #6 - t490) - t238 + t934) + t936 * (t66 * t109 * s25 + t917 + epinv
     # * (-t170 - t211 + t214 + t228 + t236 - t210) + t938 - t238) + t94
     #1 * (-t175 + t67 + t73 + t85 - t108 - t155)) - t56) + t57 * t9 * t
     #13 * t1013) + t1013 * (t396 * (-t2 * t791 + t55) - t329) * t593)
      t30 = t1015 * t9
      t39 = t871 * t1015
      t55 = t22 * t1014
      t53 = t1016 * (t11 * (-t266 - t53 - t54 + t259) - t30 * t961) + t5
     #5 * (t39 - t994 - t996)
      t54 = t30 * (-t31 - t991 - t992) * t60
      t26 = t630 * (t1015 * (t11 * (-t13 * t876 - t914 * (epinv * t478 -
     # t66 * t704 + s12 + s16 - s34 - t152 - t65) - t969 * (epinv * t326
     # + t575 + t576 + t581 + t587 + t588 - t69) + t970 * (t284 + epinv 
     #* (t256 - t264 + s56 * (s56 + t247) + t115 - t804 - t805) + (-t17 
     #* t63 - t65) * s34 + (-t826 - t372) * s45 + (t239 - t372 + t246) *
     # s46 + s56 * (t45 - t372 + s56) + t275 + t278) + t356 * t567 * t18
     #9) + t22 * (-t915 * (epinv * t572 + t196 + t204 + t921 + t922) + t
     #981 * (s25 * (t127 + t490) + epinv * (s46 * t337 + t287 * s56 - t1
     #44 - t348 + t351 + t353 - t360) - s56 * (-t26 - t345 - s56) - t366
     # + t390 + t395) - t230 * t344 * t189) + t9 * (-t23 * (t65 + s25 + 
     #s12 + s16 - s34 + s56 - s13 + s46) - t989 * (-t447 * t593 - t990) 
     #* t1016) - t36 * t19 * s34 * t359) + t593 * t53 + t54)
      t45 = -t1053 * t285 - t443 * (t1072 * t17 + t28 * (-t150 + t233) +
     # t206) + t467 * (t28 * (t920 - t150 + t924 + t72) - t918 - epinv *
     # t344 + (-s12 * t625 + t635) * s25 + s26 * (-t901 - t660) + s34 * 
     #(t1054 - t28 * (t678 + t75) - t1043 - t1058 - t1085 - t539 + t667)
     # + s45 * (-s45 * (t194 + t50) + t1114 * t17 + t28 * t692 + t699 + 
     #(t17 * t245 - t327 + t816) * s34 + t685 + t688) + (s26 * (-t260 + 
     #t288) + t583 + t722 - s12 * (-t703 + t410 - t882 + t361) + s16 * (
     #t345 - t260 + t392 - s45) + t1046) * s46 + s56 * (t919 + s12 * (t7
     #27 + t741 - t483 + t200) + s16 * (t354 + t194) + s26 * (11.0_dp * 
     #s25 - t126) - t17 * (s34 * (s25 + s16 + s46) - t802) - t28 * (t695
     # + t800 - t149 - t182 - t113) + t808 + t815) - t912) + t526 * t296
      t50 = t673 * t596
      t53 = t630 * t1015
      t54 = t533 * (s25 * (-t801 + t597) + t562 - t595 - t988)
      t64 = t220 * (-t28 * (-t150 + t924 + t145) + epinv * t1053 + ((-s1
     #6 + t851) * s16 - t579 * t810) * s25 + s26 * (-t1047 - t17 * s25 *
     # t565 - s13 * (-t882 + s16 - s13) - t75 - t953 - t579 * (t842 + s2
     #5)) + s34 * (-t1054 + s16 * (t896 + t833) + s26 * (t201 + t840 + t
     #958 + t606) + t28 * (t611 + t75)) + s45 * (s45 * (-s26 + t960) - s
     #26 * (t882 + t469) + t17 * t807 - t28 * (s25 * (-s25 + s12 - s34) 
     #- t111 - t144 + t645) + s12 * s26 + s16 * (-t882 - t505 + t507)) -
     # (t28 * t1004 + t681 + (t962 + s13 + t537) * s26 - (t582 + t327 + 
     #t846) * s34 - t608 + t967 + t976 + t983 + t1024) * s46 - (t28 * (s
     #46 * t1029 + t111 + t182) + t681 + (t16 * t565 + s13 + t582) * s26
     # - (t201 + t892 + t327 + t605) * s34 + t585 + t976 + t1032 + t1033
     # + t1034) * s56 - t113 * t1029 - t71 * t1029 + t1055) * t189
      t67 = t1025 * t866 * t593
      t38 = t274 * (t22 * (t1013 * (t1015 * (t1014 * (-t536 * (s26 * (-t
     #883 + s26) + s56 * (s56 + t569 + s45 - s13 + s46) - t232 * t89 + t
     #599 + t611) + t593 * (-t868 + t869 + t870) - t863 * t610 - t54 - t
     #2 * t909 * t596) - t270 - t64) + t67) + t786) + t212 * (s34 * (t40
     # * t9 * t365 + (-t38 * t396 + t43 * t560) * t593 * t2) - (t32 + t6
     #68 + t899) * t9 * t60 * epinv) + t212 * (t1016 * (-t574 * t799 - t
     #985 * (s34 * (-t986 + t597 + s34) + t804 + t805 - t987 + t988)) + 
     #t189 * (-t394 * (-t512 + epinv * t1083 + s26 * (-t1047 - t1038) + 
     #s34 * (-s34 * t1041 + t28 * ((t816 + t577) * s25 - t608 + t839) + 
     #t66 * ((t566 + t601) * s26 + t144) + t377) + s45 * (-t49 - t1048) 
     #- t1045 - t1049 - t1050 - t1051 - t1052 + t1055) - t612 * t964 - t
     #1027 * t567) - t77 * t778) * t11)
      t49 = epinv2 * t593
      t4 = t274 * (t212 * (-t49 * (t5 * (t4 * (s25 * (propZ25 * t16 * t2
     #7 + t3) + t830 + t363) - t20 * t571) + t21 * (t854 + t744)) * t101
     #6 + t488 * (t593 * (t907 - t908) + t610 * t904) - t593 * t1021 * t
     #19 * s34 + t22 * t593 * t1018 * t29) * t1014 + t1013 * (t1016 * t4
     #7 * t593 - t13 * t798 - t359 * t84) * t9 + (t1015 * t58 * t60 + t3
     #65 * t83) * t1013 * t22 - t191 * t793 * t11 + t301 * t212 * t11 * 
     #(t25 * t593 - t541 * t879 + t574 * t878))
      t5 = t274 * (-t22 * (t1013 * t365 * t78 + t787) + t212 * (t13 * t2
     #2 * t777 + t185 * t359) * s34)
      t13 = t120 - t61
      t20 = intHs16s25s26s34s56x1411D6eps1() * t2
      t21 = intHLs16s25s26s34s56x1114D6eps0() * t22
      t25 = t274 * (t1013 * (t22 * (t1015 * (t814 * (s45 * t516 - (s12 +
     # s16 - s34 - s45) * s56 + t241 - t828 - t85 - t535) + t926 * (-epi
     #nv * (-s34 * t450 - s45 * t51 + t117 + t458 + t460) + t196 + t204 
     #+ t921 + t922) + t942 * t776) - t365 * t797) + t359 * t135 * t1015
     # + t1015 * ((-t69 + t587 + t581 + t576 + t588 + t575) * t13 + t942
     # * t134) * t11) + t9 * (-t20 + t402) - t21)
      t47 = 0.1e1_dp / 0.3e1_dp
      t60 = 0.1e1_dp / 0.18e2_dp
      t1 = t1 * t97 - t82 * t5 - (t4 + t38 + t53 * (t11 * (-t10 * t24 - 
     #t1005 * t478 - t1007 * (epinv * t473 - t468 * t66 + s12 - s34 - s4
     #6 - t152) - t189 * (t1083 * t203 + t209 * t99) + t33 * (-epinv * t
     #24 - s13 + s25 + s46 + s56 + t152 - t476 - t65) - t945 * t326 - t9
     #48 * (-t256 - s56 * (s56 + t247) - t115 + t804 + t805 + t264) - t9
     #97 * t473 + t809 * t541 * t1016) + t19 * (-t1008 * (epinv * t474 +
     # t28 * (s25 - s16 + s45 + s56 - s13)) - t34 * (epinv * t24 + t589)
     #) + t22 * (-t1000 * t474 - t1006 * t478 - t15 * t24 + t189 * t45 -
     # t86 * t572 - t910 * (epinv * (-t603 + t734 + t133) - t913 - t590)
     # - t959 * (-s46 * t337 - t287 * s56 + t144 + t348 - t351 - t353 + 
     #t360) + t50 * t1014)) + t26) * t60 + t25 * t47
      t4 = s25 * t237
      t5 = s16 + s26 - s34 - s45 - s46
      t24 = s16 * (-t28 * t552 + s25 + t201 + t234 + t330)
      t25 = t66 * (-t244 + t339 + t111)
      t26 = t28 * (-t745 + t144)
      t38 = t464 * t5
      t45 = s34 * t110
      t61 = t28 * t104
      t64 = s26 * (t891 - s26) - (s25 - s26 - s13) * s34 + s46 * (-s46 -
     # t65 - s16 + s34 - s45 + s13) + s56 * (-s56 - t584 - t834) - t821 
     #+ t828 - t860
      t62 = t28 * (t144 + t113 + t71) + t848 + (-t62 + t816) * s26 - s34
     # * (-t62 + t200 + t194 - s34) + s45 * t937 + (t566 + t201 - t330 +
     # t842) * s46 + (t17 * t893 + s12 + t201 - t483) * s56 + t608 - t82
     #8 + t887
      t69 = t237 * t864 - t865
      t73 = t12 * t87
      t77 = s12 * (-t330 + s56 + t772) + t73
      t82 = -s13 + s25 - s46 + s56
      t85 = t213 * s45
      t108 = t314 * s34
      t5 = s13 * t5
      t109 = s12 * (-s12 * t82 + s46 * (t201 - t330) + t28 * (-s26 * t15
     #4 + t113 - t146) + s16 * (t65 + t523 - t181) - t168 + t85 + t108 +
     # t5) + t12 * t223 * t104
      t120 = t504 * s45
      t135 = s12 * (-s12 * (t28 * (s25 - s16 + s34 - s13 - s46) - t521) 
     #+ s26 * (-t260 + t507) + s46 * (t898 + t523 - t505) + t124 * (t122
     # + t144) + t28 * (t5 - t120 + t45 - t168) - t66 * t475 - s16 * (t2
     #8 * t551 + s45 - t343 + t469 - t898) - t161 - t338) + t12 * (t834 
     #+ s16) * t104
      t155 = s16 + s26 - s46
      t162 = -t80 + t155
      t164 = -s13 - s16 + s25 + s34 - s46 + s56
      t170 = t28 * (-t146 + t111 + t113)
      t175 = t164 * s12
      t181 = s12 * (t175 + t124 * t749 - t1067 + s16 * (-t201 + s25 + t3
     #72 - t199) + t120 + t168 - t45 - t5 - t170) - t104 * (s16 * (s16 +
     # t834 - s34) + t223 * t506)
      t182 = t124 * t364
      t185 = -t387
      t196 = t93 * t16
      t204 = t17 * t339
      t206 = t455 * s26
      t164 = t37 * t164
      t207 = t404 + t414 - t446
      t209 = t89 * s25
      t210 = (s25 + s34 + s45) * s34
      t211 = (-s46 + t693) * s45
      t214 = t404 - s45 + s46
      t228 = t400 + s46
      t233 = -t404 - t400
      t236 = -t414 - t446
      t238 = t236 * s56
      t239 = s25 + s45 - s56
      t100 = t12 * t100 * t104
      t6 = t44 * (t66 * (t598 * t784 + t8) + t28) + t8 * (t7 * (t28 * t5
     #98 + 1.0_dp) + t6)
      t7 = -s13 + s25 + s26 - s46 + s56
      t8 = -t128 + s34
      t44 = t330 * t223
      t245 = s16 * (t864 - t303 - t133)
      t246 = t197 * t124
      t247 = t66 * t418
      t249 = t408 * t148
      t254 = t980 * s46
      t256 = -t447 + s34
      t108 = s12 * (t37 * t7 - t124 * s46 * (t404 + s26) - t17 * (-t219 
     #+ t244) + t28 * (-t425 + t254 + t249) + t66 * (-t165 * t440 - t108
     # - t174) - s26 * (-t617 - t407 - s34) - s34 * (-epinv * t314 + s34
     #) - s45 * t256 - t2 * (t5 - t245) + t442) - t104 * (t28 * (t4 + t2
     #44 - t653 + t142) - t614 * t8 - t655 + t247)
      t259 = -s13 + s25 + s26 + s34 + s56
      t264 = -t123 + t483
      t266 = s16 * (t51 - t128 + t372)
      t270 = epinv * t148
      t275 = t148 * s46
      t278 = t94 * epinv
      t284 = s34 * t165
      t247 = s12 * (-t124 * t418 + t28 * (s45 * (t270 - s25) - t267 - t2
     #75) - t66 * (-t677 + t284) + t662 + s25 * t971 + s34 * (t626 - s56
     #) + s46 * (t278 - s25 - s45) + t2 * (t266 - t5) + t417 - t425 + t4
     #42 + t459 + t37 * t259) - t104 * (t237 * s56 - t264 * t614 - t197 
     #+ t247 + t4 + t656 - t95)
      t286 = t28 * (s25 - s13 - s46) + s26 + t538
      t287 = -t141
      t288 = -s26 + s34 + s46
      t292 = s34 * t87
      t296 = s12 * (s12 * (t744 + s25 + s16 + s34 - s13 + s46) + t25 + t
     #26 + t24 + t193 - t5 - t389 - t522 + t71 + t752) + t292 * t61
      t301 = t404 + s56 + s46
      t320 = t100 * epinv
      t323 = -t242 + t80
      t82 = s12 * (t28 * (s26 * (-t404 + s25) + t238 - t425 - t442) - t6
     #6 * (t542 - t93 - t709) + t427 - s16 * (t17 * t414 + t28 * t441 + 
     #s16 + t127 - t408 - t617 + t640) - s34 * (t617 - t447 + s45) - s46
     # * (-t446 + s25) - t149 - t38 + t723 + t37 * t82) - t104 * (s16 * 
     #(t414 - t753 - s25 - s45 - t483 + s16) + s25 * t323 + t139 * (t2 *
     # t223 - s34))
      t326 = -t124 + t690
      t328 = t250 * s26
      t183 = s12 * (-t37 * t286 + s46 * (t3 * t404 - t758) - t124 * (-s5
     #6 * t740 - t756) + t16 * t542 + t17 * (t425 + t122 + t290) + t28 *
     # (s34 * t256 + s45 * t771 + t38 - t623 + t643) + t66 * (-t288 * t4
     #11 + t184 - t328) + 7.0_dp * t198 - s16 * (-t614 - s46 * (9.0_dp -
     # t435) - s56 * t326 - t28 * (t447 + s26 + s45) + t66 * (-t428 + s2
     #5 + s34) - t404) + t111 - t424 - t648) - t104 * (-t17 * t418 - t28
     # * (epinv * (-s26 * t631 + t279) + t149 + t183 - t291) + t750 - s1
     #6 * (-t17 * t616 - s25 - s56 + t614 + t65 + t665 - t743) + s25 * (
     #-t671 + s26 - t68) - t111 + t243)
      t243 = t28 * t480
      t256 = -t243 + t715
      t290 = s12 + s16 + s26 + s46
      t291 = t87 + t579
      t337 = -t28 * s46 * t806 + s16 * (t235 + t287) + s26 * (-t483 - t5
     #52 + s26) - t279 - t836
      t344 = t28 * (t530 + t950 + t149) - t662 + s12 * t952 + t111 - t12
     #2 + t144 - t272 - t331
      t348 = s16 * t123
      t351 = t237 * s56
      t292 = t534 + t292
      t353 = s12 * t957
      t359 = t12 * t27
      t360 = t237 * t27
      t365 = -s12 * t166 + t360
      t360 = t360 * epinv
      t366 = s16 * (-s16 - t65 + t715) - t139 * t223 + t965
      t387 = s26 * t12 + t835 - t836
      t396 = (-t80 + s12 + s16 + s26 - s46) * s56
      t410 = t1001 + t235
      t418 = t235 + s12 + s16 + s26 + s34 + s45
      t428 = epinv * t418
      t450 = t2 * t465
      t87 = t274 * (t1013 * (t1015 * (t9 * (-s12 * t57 - t77 * t968 + (s
     #12 * (t87 - t330) + t524 * t87) * (t450 - t944)) + t19 * t941 * t2
     #92) - t74 * epinv * t9 * t12) + t22 * (t1015 * (t1013 * (t337 * t8
     #25 + t674 * (epinv * t292 + t28 * (t351 + t267 + t922 + t921)) + t
     #823 * t387 + t927 * (epinv * t337 + t396 - t493 + t575 - t934 + t9
     #35) + t936 * (epinv * t387 + t28 * t575 - t542 * t66 - t939 + t988
     #)) + t56) + t995 * t593 * t1013 * t1014 + t212 * (t291 * t777 + t8
     #73 * (epinv * t291 - t14)) * s34 + t212 * (t177 - t877) * t237) + 
     #t212 * (s34 * (-t256 * t872 - t702 * (epinv * t256 + s16 + s26 + s
     #46 + t579)) - t955 * (t353 + t359) - t977 * (s12 * (-t28 * (-t613 
     #+ s56 + s46) + t978) + t359 * epinv) + t37 * t189 * (-t62 * t881 +
     # t64 * t880)) * t11)
      t7 = t109 * t229 - t181 * t225 + t195 * t108 - t295 * (s12 * (s12 
     #* t259 - s46 * t94 + t124 * t339 - t113 + t193 + t255 + t266 - t5 
     #- t522 - t561) + t104 * (s16 * t264 + t44)) - t297 * (s12 * (s12 *
     # t7 - t28 * (t244 + t122 + t113) + t66 * t666 + t111 - t160 + t193
     # + t245 - t5 + t246) + t104 * (s16 * t8 + t44)) + t308 * t82 + t45
     #6 * (s12 * (t164 + s34 * (t523 + t507) + t124 * (s46 * t233 - t756
     #) + t28 * (-t442 + t238 - t425 - t542) - t66 * (t423 + t522) - epi
     #nv * (-t267 - t149) + s16 * (-t124 * t228 - t16 * t710 + s45 + t40
     #7 - t972) - t184 - t211 - t38 + t764) + t104 * (-t28 * (t616 + s34
     # + s45) * s34 - s16 * (t614 + t28 * (-t616 + t400) - t330 + t214) 
     #+ s26 * (s34 * (t66 - t412) - t414 + t446) - s46 * (-t446 + s25 - 
     #s34) - t442 - t542)) + t650 * t247
      t8 = t2 * t788
      t44 = t2 * t790
      t82 = t945 * t365
      t84 = epinv * t84
      t73 = t220 * (s12 * (t164 - t17 * s34 * (t440 + t714) + t28 * (s46
     # * (-t617 - s25) - t419 - t425 - t442) - t66 * (t423 + t174 + t210
     #) + t196 + t427 + s16 * (-t124 * t451 + t407 + t446 - t696) - t38 
     #+ t764 - t211 + t404 * t110) + t104 * (-t28 * (-t209 + t210) - s16
     # * (t614 + t299 - t330 + t207) + s26 * (t330 - t207) - s46 * (-t61
     #3 - t446) - t442)) - t230 * t296 + t285 * (s12 * (-t175 + t202 + t
     #1067 + t170 - s16 * (t136 - t859 + t372) - t120 - t168 + t45 + t5 
     #- t204) + t73 * t104) + t443 * t183 - t526 * (s12 * (-s12 * t286 +
     # s46 * (t545 + t523 - t505) - t124 * t173 + t28 * (t5 + t85 + t267
     # + t111 - t168) + t66 * t288 * s56 - s16 * (-t66 * t287 - s16 + s3
     #4 + t234 + t327 - t343) + t161 - t71) + t104 * (s16 * (-t178 + t12
     #7 + t523 + s16) + t834 * t506))
      t49 = t22 * (-t50 + t54) + t49 * t6
      t40 = t22 * (s12 * t795 - t340) + t11 * t40 * t42
      t42 = epinv * t237
      t35 = s34 * t40 + t1014 * t49 + t12 * (t11 * (epinv * (t32 + t668)
     # + t31 + t991) + t22 * (-t35 - t15) + t9 * (-t283 + t10)) + t189 *
     # (t19 * t73 + t9 * (t612 * (s12 * (s12 * (t464 + t619 - t407 + t61
     #4 - t613 + s26 + s56) + t124 * t301 * s26 - t17 * t554 - t28 * (s3
     #4 * (epinv * t239 + s56) - t113 + t459 + t623 - t71 - t643 + t408 
     #* t205) - t66 * (s26 * (-s26 + t411) + t197) - s16 * (t627 * s16 -
     # t124 * (t404 + t400) - t66 * t730 - t543 - t684 + t732) - s34 * (
     #-t616 - s25 - s45) + t149 + t38 + t492 - t634 - t637) + t320) + t9
     #9 * (s12 * (s12 * t998 - t28 * (s34 * t239 + t111 + t149 - t168 + 
     #t763) - t66 * (t179 + t144) + s16 * (-t182 + s45 + t327 - t538) - 
     #t130 + t197 - t5 + t201 * s34) + t100)) - t265 * t135) + t42 * (t1
     #1 * t778 + t899 * t9) - t37 * t879 * t9
      t37 = (-t2 * t43 + t780 + t781) * t19
      t40 = t536 * (t237 * t51 + t599)
      t24 = t19 * t467 * (t61 * (-s16 * (t613 - s45 - s46) + s34 * (-t40
     #0 + s34 + s56) + s46 * (-t613 + s26 + s56) - t677 + t4) + s12 * (s
     #12 * (epinv * (t744 + s25 + s16 + s34 + s46) - s25 - s26 - s34 + s
     #46 - s56 + t464) - t124 * t208 * s46 - t17 * (s26 * t178 + t146 + 
     #t661) + t28 * (t350 + t142) + t66 * (s45 * (s25 + s16 + s34) + t45
     #) + epinv * (t25 + t26 + t752 - t747 + t24 + t193 - t522) - s56 * 
     #t101 - t218 - t575 + t38))
      t5 = t1013 * t189 * (t9 * (t203 * (s12 * (t175 - t124 * t185 - t28
     # * (-s45 * t336 + t111 + t113 + t122 + t45) - t1067 + s16 * (-t182
     # + s25 - t199) + t168 - t334 - t5 + t196) - t104 * (t28 * (-s16 * 
     #t368 - t197) + s26 * (s26 - t330 - t552) + t144 + t149 + t174 + t3
     #32 + t334)) - t394 * (s12 * (t164 + t124 * s46 * (-t404 + t206) + 
     #t28 * (s45 * (-epinv * t336 - s56) - t420 - t425 - t442 + t459 + t
     #624) + t204 - t975 - s16 * (t124 * (t404 + t400 + s46) - t372 - t4
     #34 - t759) + s45 * (t414 - s25 - s34) - t38 - t45 + t453 + t764 - 
     #t1067 * t2) + t104 * (epinv * (t334 + t332 + t149 + t174) + s16 * 
     #(-t28 * (t400 + t409) + t715) + s26 * (-t399 - t414 + t446 + s34) 
     #+ t115 + t403 - t432 * t28))) + t1028 * t135 - t24)
      t24 = s12 * t231
      t5 = t274 * (t1015 * (t793 * t11 * s34 + t1013 * (t593 * (t22 * (t
     #1019 * t2 + t29 * (-t905 - t906) + t868 - t869 - t870) + t37 * s34
     #) + t40 * t19 + t488 * t909 * t596 + t19 * (-t2 * t904 + t863) * t
     #69) * t1014 + t5 + t24 * t86 * t22 * t1013) + t1013 * (t22 * t83 +
     # t798 * t9) * t237 - t67 * t22 * t1013)
      t5 = t5 + t53 * t35 + t630 * (t1015 * (-t19 * t366 * t959 + t23 * 
     #t290 * t9) + (t1005 * t418 + t1007 * (s12 * t29 + t12 * t2) + t876
     # * t237 + t914 * (-t243 + t428 - s34 - s45 - t68) + t948 * t344 + 
     #t969 * (s12 * (t619 - t613 + s16 + s26 + s56) + t360) + t970 * (-t
     #28 * (s26 * t801 + t588) + epinv * t344 + s34 * t607 + (-t12 - t57
     #9) * s46 - t111 - t144 + t802 - t836) + t997 * (s12 - s16 - s26 + 
     #s34 + s45) + t82) * t1015 * t11 + (-s12 * t58 + t1000 * t410 + t10
     #06 * t418 + t1008 * (epinv * t410 - t590) + t910 * (t428 - t14) + 
     #t915 * (t231 * t925 + t267 + t351 + t921 + t922) - t981 * (epinv *
     # t366 + t17 * t267 + t390 + t395 - t396 - t575)) * t1015 * t22 + t
     #9 * (t1015 * t992 + t84) * t237 + t55 * (-t39 + t994 + t996) * t59
     #3 + t30 * (t809 + t989) * s12)
      t5 = t47 * t147 * (t1023 * (t22 * (t1015 * (s34 * (t290 * t776 + t
     #304 * (epinv * t290 - s12 - s16 - s26 + s34 + s45)) - t814 * (-t93
     # + t535 + t818) + t926 * (-epinv * (-s46 * t806 + t348 - t93) + t2
     #67 + t921 + t922 + t351)) - t12 * t797) + t191 * t290 * t134 * t11
     # + t138 * (s16 + s26 - s34 + s46 + s56) * t9 * s12) + t1015 * t118
     # * t593 * t785 * propW16) + t60 * t5 + t70 * t147 * propW16 * s16 
     #* t785 * t593 * t1015 * t88 - t97 * (t274 * (-t22 * (t1013 * t162 
     #* t78 + t1015 * t187) + t1013 * (t482 * (s12 * t794 + t762 * t77) 
     #- t129) * t9 + t1013 * (-t488 * t867 + t191 * (-t52 + t782) * t19)
     # * t1014 * t593 + t212 * (t11 * (s12 * (t62 * t698 - t64 * t837) +
     # t109 * (-t8 + t190) + t181 * (t44 - t489)) + t19 * t7) * t189) + 
     #t87)
      t7 = s12 - s34 - s45 - s46 - s56
      t14 = s13 * (s16 + s26 + s34)
      t25 = t28 * t829
      t25 = t28 * (s16 * (t1060 + t828) - t114 - t72) + s26 * (s26 * (t8
     #50 + t577) + s12 * (t816 + t850 + s12) + t842 * t63) + s34 * (-s12
     # * t916 - t101 * t578 + t267) + s45 * (-s45 * (s16 + s26 - s34 - s
     #13) - s12 * t890 - t28 * (s16 * t110 + t111 + t116) + t14 - t144 +
     # t45 - t194 * s26) + s46 * (s46 * (-t59 - t201 + s12 + s34 + s13 -
     # t842) + t28 * (t585 + t828 - t111) + s16 * (t582 - s13 - s16) + (
     #t1030 - s13 + t327 - t194) * s26 + s34 * (-t101 - t850 - t577 + s3
     #4) + (-t816 - t25 + s34 - t507) * s45) + s56 * (s56 * (-t59 - t201
     # - t842 + s12 - t199) - t113 * t16 + t28 * (t828 - t608 - t144) - 
     #t848 + (-t816 + t850 + t605) * s26 + s34 * (-t101 - t25 + s34) + s
     #45 * (-t850 - t200 - t194 - s45) + (t890 - t392 - t884 - t152 + t5
     #79) * s46 + t1119 + t75) + t588 * t820
      t26 = t28 * t287
      t14 = s26 * (s26 * (s12 + t490 - t886) + s12 * (t345 + s16 - t886)
     # - t842 * t883) + (t101 * t886 - t882 * t986) * s34 + s45 * (s45 *
     # (-t103 + s25 + s34 + s13) - t124 * t575 + t66 * t1063 * s34 - s12
     # * (t1063 * t28 + t890) + s16 * t101 + t14 + t149) + s46 * (s46 * 
     #(t857 - t26 - s46) - t931 + s16 * (-t943 + s16) + (t810 + t882 + t
     #605) * s26 + s34 * (s34 - t28 * (s12 + s16 + s26 - s13) - t882) + 
     #s45 * (-t601 - t200 + t886 - s45) + t75 + t828 * t66) + s56 * (s56
     # * (s12 - s16 - s34 - t26 - t68 - s56) - t28 * (t167 + t608) + t66
     # * (t828 - t113) + s26 * (t882 - s16 + t929 - s26) + s34 * (-t28 *
     # (s16 + s26 - s13) - s12 - t882 + s34) + (t28 * (s12 - s26 + s13) 
     #- t1096 - t882) * s45 + (-t17 * t287 + t28 * t849 - t1101) * s46 +
     # t839) + t828 * (t842 + s12)
      t14 = t28 * (-t592 * t801 + t216) + t14
      t26 = t12 * t237
      t30 = s12 * (t12 * t66 + s56) - t26 + t609
      t35 = s12 * (s12 - t59 + t155) - t237 * t524
      t38 = t28 * t231
      t39 = t38 + s45
      t43 = s25 - s45 + s46
      t45 = s13 * t39
      t49 = t124 * t1003
      t50 = -t124 * t271
      t52 = t28 * t217
      t53 = t52 + t552
      t52 = t52 - t200 + t552
      t54 = t28 * (-t331 - t272 + t111 + t144)
      t62 = s25 * (t217 + t128)
      t63 = t306 * (-t172 + t111)
      t64 = s34 * t166
      t67 = (s16 * (-t50 - s46 + t496 - t309) + t66 * (-s26 * (-t17 * t2
     #71 + t469) + t71) + 12.0_dp * t522 + t686 + t167 + t62 + t63 - 13.
     #0_dp * t142) * s16
      t73 = s13 * (s16 * t52 + s26 * t53 + t253 + t525 - t54)
      t77 = t279 - t71
      t87 = t158 * t136
      t94 = -t149 - t71
      t108 = t124 * t172 * t119
      t109 = t17 * (s34 * t94 + t112)
      t110 = t51 + t235
      t116 = t124 * t51
      t118 = t3 * t121
      t120 = t17 * t389
      t127 = t17 * t51
      t128 = s45 + t127 + t68
      t129 = t279 + t142
      t134 = s25 * (-t497 + t448 + t496 - t361)
      t135 = t3 * (-t172 + t111 - t197) + t66 * (-t158 - t113 + t71) - t
     #1086 - 9.0_dp * t129 + s26 * (-t16 * t552 - t469 + t546) + t134
      t138 = -t279 + t168
      t146 = (t193 + t111) * s46
      t147 = t192 * t139
      t155 = (-t167 + t222) * s56
      t160 = t167 * t232
      t164 = t119 * t142 - t112 - t324 - t352 + t380 + t388
      t170 = t949 * t111
      t173 = 12.0_dp * t179 * t79
      t175 = t110 * s12
      t181 = t26 * t106
      t182 = s12 * (s12 * (-t124 * (t240 - t197) - t66 * (s56 * t314 - t
     #241 + t334) - s16 * (t59 + t116 + t199) - t113 + t45 - t677 - t118
     # + t120 + t174 * t28 - t175) + t16 * s34 * t138 + t124 * s45 * (s5
     #6 * (-s34 + s56) + t168) + t17 * t164 + t28 * (s46 * (-t90 - t167 
     #+ t71) - t216) - t66 * (s34 * (s45 * t232 + t93) + t282 + t146 + t
     #147) + s16 * (-s16 * t128 - t135) - t385 - t73 - t742 + t155 - t16
     #0 + 7.0_dp * t280 + 7.0_dp * t170 + t173 + t222 * (t260 + t376)) +
     # t181 * t104
      t183 = t110 * epinv
      t194 = t183 + s34
      t196 = t124 * t651
      t201 = t28 * t51
      t204 = s25 * t341
      t205 = -t204 + t71
      t39 = t464 * t39
      t207 = 9.0_dp * t222
      t52 = -s16 * t52 - s26 * t53 - t253 - t525 + t54
      t53 = -s26 * t726
      t54 = -s25 + s26 - s34 + s56
      t208 = t464 * t52
      t210 = s26 * (t111 + t71)
      t211 = s25 * (t157 + t172 + t459 - t111 + t142) + s34 * (t172 + t1
     #67) + t462 + t761 - t210 + t424 * s46
      t218 = t316 * epinv
      t238 = s26 * (s46 * (-t447 - t414 + s56) - t411 * t314) + s45 * (t
     #559 - t71) - t150 - t216 + t500 - t218
      t239 = (-t623 - t167) * s25
      t241 = epinv * (s46 * (s45 * t381 + t149) + s56 * (-t161 + t149))
      t243 = s34 * (-t755 - t158 - t71) + t349 * epinv
      t245 = s46 * (s34 * (12.0_dp * epinv - 14.0_dp) - t376) + 17.0_dp 
     #* t422
      t247 = s12 * (t194 - t235) + t124 * (-t1009 * t440 - t122) + t17 *
     # (-s26 * t289 - t172 + t254) - t28 * (s45 * (-t414 + t712 + s34) +
     # t192 + t240 + t71) + t3 * t440 * t119 + t66 * (epinv * t205 - t43
     #1) - s16 * (-t17 * t449 - t196 + t201 + t404 + 7.0_dp * t414) - s3
     #4 * t439 - t425 + t39 + 7.0_dp * t709 - t207 - t422 * t16
      t253 = t237 * t104
      t215 = t253 * t215
      t53 = s12 * (s12 * t247 - s26 * t245 - t124 * t243 - t16 * (s34 * 
     #(s25 * (-s46 + t440) + t634 + t637) + t169 - t471 + t279 * t980) -
     # t17 * (s56 * (-t425 + t149 + t111) + t706 + t239 - t623 * t631) +
     # t28 * t211 + t3 * t241 + t306 * t93 * (s56 - t440) - t66 * t238 -
     # epinv * (t167 * t54 - t742) - s16 * (s16 * (-t232 * t28 * t413 - 
     #t124 * t236 + t17 * t730 - t28 * t663) - t124 * t408 * t156 + s34 
     #* (s46 * (14.0_dp * epinv - 13.0_dp) + t469 + t411 * t391) - t16 *
     # t415 + t28 * (-t639 - t459) - t66 * (s46 * (t414 - s56) + t634) -
     # 9.0_dp * t142 - 9.0_dp * t664 + s25 * (t16 * t974 - t429 * t66 + 
     #s34 + t53 + t728 + t999) + s26 * (t53 + s45 * (t3 - epinv) + s56 *
     # (t16 - t729) + 9.0_dp * t710 + 12.0_dp * t974)) - t379 + t208 - 7
     #.0_dp * t111 * t433 - 7.0_dp * t400 * t149 - 7.0_dp * t687 + 7.0_d
     #p * t158 * (s25 - s46 + t411) + 9.0_dp * t122 * t79) + t215
      t156 = t714 + s34
      t169 = -t142 + t168
      t211 = s45 * (s25 * (epinv + 1.0_dp) + s34)
      t215 = t429 - s45
      t238 = (1.0_dp + t398) * s34
      t241 = s46 + t408
      t243 = 23.0_dp * epinv
      t245 = 19.0_dp
      t247 = 32.0_dp * epinv
      t254 = -27.0_dp
      t255 = 12.0_dp * s34
      t256 = (t404 + s26) * s56
      t259 = t222 * t79
      t264 = s16 * (s16 * (t28 * t720 * s25 + s34 * (t391 - t690) + s45 
     #* (t66 + t412) + t306 * t766 - t312 - 14.0_dp * t400 + 13.0_dp * t
     #411) + s34 * (-17.0_dp * s34 + 23.0_dp * s46) + s45 * (s46 * (21.0
     #_dp - t768) + t133 - t255) + s56 * (s45 * (-epinv * t245 + 11.0_dp
     #) + s46 * t326 + t404 * t254 + t435 * s56) - t16 * (t421 - t425 + 
     #t71 + t760) + t28 * (t111 * (-t3 + t729) + t459) - 24.0_dp * s34 *
     # t766 + s25 * (t16 * t215 + t17 * (s26 * t720 - t238) - t306 * t24
     #1 - s56) + s26 * (s56 * (-21.0_dp + t247) - t414 * t245 + (-t243 +
     # 31.0_dp) * s34 + (9.0_dp - 13.0_dp * epinv) * s45)) + s56 * (-s25
     # * t437 + t634)
      t265 = t122 + t111
      t266 = t455 * s56
      t194 = s12 * (-t194 * t28 + t523 + t956) + t124 * t71 - t16 * (-t3
     #41 * t447 - t149 + t664) + t17 * t211 + t28 * (s45 * (t404 + t615 
     #+ t426) + t267 - t39 + t425) - t306 * epinv * t169 + t66 * t773 + 
     #17.0_dp * t416 + 12.0_dp * t454 - 9.0_dp * s56 * t156 + 11.0_dp * 
     #s26 * t441 - 7.0_dp * t634 + 7.0_dp * t111 + 7.0_dp * t184 + s16 *
     # (-s16 * t721 + s25 * (-t435 + t66) + s45 * (-t430 + t28) - 12.0_d
     #p * t270 + t342 + 11.0_dp * t364 + t404 + 14.0_dp * t414) + t440 *
     # (-18.0_dp * s56 + t505) - t226 * s34
      t270 = -t237 * t104
      t272 = t270 * (t751 - t750 + t654 - t747 - t748 + t735 + t111)
      t146 = s12 * (s12 * t194 - (s34 * (s26 * (-t247 + 25.0_dp) + t226)
     # - 26.0_dp * t648 - t179 * t124) * s56 + t16 * (-t168 * t215 + t17
     #4 * t429 + t112 - t150 - t216 - t218 - t500) + t17 * (s34 * (-t167
     # + t168) - t462 - t738 - t761 + t415 * s46) + t245 * t623 * t217 -
     # t28 * (s56 * (s45 * t439 + t267) + t145 * t29 + t208 - t352 - t61
     #6 * t71 + t440 * t167 + t408 * (t322 + t149)) + t3 * (t532 + t239 
     #- t163) + t306 * (t146 + t440 * (t347 - t149 - t111 - t71) - t742)
     # + t391 * t408 * t265 + t66 * (t122 * (t400 - s56) - t471) - 14.0_
     #dp * s45 * (t267 + t256) + 12.0_dp * s34 * (s26 * (t447 + s45) + s
     #34 * (s46 * t713 - t411)) - 24.0_dp * t259 - 16.0_dp * t279 * (-t6
     #13 + s56) - 16.0_dp * t176 - t264 - 9.0_dp * t501 - 9.0_dp * t642 
     #+ 7.0_dp * t113 * (t266 + t206) + 11.0_dp * s56 * (s34 * (t719 + s
     #34) + t111 + t636) + 23.0_dp * epinv * t259 + t369 * t149) + t272
      t194 = t617 + t393
      t215 = epinv * t552
      t218 = t222 * epinv
      t226 = -t172 + t149
      t239 = -t167 - t142
      t247 = t217 * s45
      t264 = t267 * t79
      t272 = t315 + t158
      t280 = t267 * t119
      t249 = s12 * (-t28 * t974 - t194) + t124 * epinv * t1003 - t16 * (
     #t422 + t222) + t17 * (t421 - t179) + t28 * (s25 * (-t414 - s26 - s
     #56) + t389 + t279 * t457) + t66 * (t637 + t634 - t111 + t197 + t21
     #8) - s16 * (-t629 - epinv * (-t260 + t496 + t125) - s34 + t215 + t
     #538 - t640 + t859) - s25 * (-t278 - s34) - t249 + t39 - t425 - t43
     #2 - t71 - t947 + t440 * (-t260 + t311)
      t62 = s12 * (s12 * t249 + t124 * (s16 * t179 + t276 + t282 + t388)
     # + t16 * s46 * (t161 + t767) + t17 * (s16 * (-t168 + t113 + t71) +
     # t508) + t28 * (s16 * t272 + s26 * t272 + s45 * (t192 + t113 + t71
     # + t184) - t111 * t444 - t144 * t551 + t114 + t72 + t280) + t66 * 
     #(s26 * (-s16 * t89 + t179) + t657 + t267 * t101) + epinv * (-t28 *
     # (s46 * (t695 + t267 + t279) - t145 + t321) - t66 * t87 - t108 + t
     #109 - s25 * (t130 + t149) - s26 * (-s26 * (s46 + 7.0_dp * t271 - t
     #573) + t66 * t77 - t167 - 12.0_dp * t226 - t62 - t483 * (-t367 + t
     #257)) - s46 * (t64 + t659) + t383 - t548 + t67 + 7.0_dp * t149 * t
     #271 + t532 * t3) - s16 * t239 - s26 * t239 + s45 * t144 - s45 * (t
     #247 + t149) - s46 * (s34 * (-s16 - s26 + s34 + s45 + s56) - t144) 
     #- t171 + t208 - t221 - t264 + t122 * (t392 + t884)) + t320 * t237
      t89 = 14.0_dp * s46
      t239 = 23.0_dp * s34
      t249 = -t279 + t111
      t272 = t245 * s45
      t265 = t265 * s45
      t276 = 14.0_dp * s56
      t160 = s26 * (s26 * (t312 + t357) + s56 * (-32.0_dp * s34 - 26.0_d
     #p * s45 + t68)) + t16 * s25 * (-s34 * t718 + s45 * (s34 - s56 - s4
     #6) + t1011) - t17 * (s46 * (-t90 - t167) - t216) - t245 * t111 * t
     #217 + t28 * (-t472 + t383 + t742 + t160 + t73 + t145 + t385) + t30
     #6 * s26 * (s45 * (-s25 + s34) + t111 + t149 + t71) - t391 * t265 -
     # 23.0_dp * t259 - 12.0_dp * (-t331 + t168) * s34 - 11.0_dp * s56 *
     # (s56 * (s34 + s45) + t197) - 7.0_dp * t113 * t148 + s16 * (s16 * 
     #(t59 + t343 + t564 + t357 - t372 + t355) + s26 * (s46 * t245 + 32.
     #0_dp * s56 - t239 - t375) + s56 * (s34 * t254 - t272 + t523 + t528
     #) - t16 * (t525 + t113) + t28 * (t149 + t134) + t556 * t249 - 24.0
     #_dp * t197) + t72 + t158 * (t276 + t549)
      t245 = t28 * t76
      t254 = s16 + s26 - s34 - s45 + s56
      t259 = s12 * t254 - t26
      t278 = t457 * s16
      t162 = s12 * (s12 + t162) + t270
      t283 = t509 - s12 + s16 + s26 + s34 + s45
      t286 = s12 * (s12 + t590 + s56) - t26
      t288 = (s12 + s34 + s45 + s46) * t231
      t289 = s16 + s26 - s56
      t290 = s16 * (-s46 + t327)
      t291 = t139 * s46
      t140 = t137 + t140
      t292 = t882 * t139
      t181 = s12 * (-s12 * (t106 + t882) - s16 * (t744 + t882 + s46) - t
     #106 * t140 - t292) + t181
      t296 = s12 * s46
      t316 = -s16 * t982
      t185 = t28 * t185
      t320 = s16 * (t250 - t235)
      t326 = (t133 + t372) * s46
      t336 = -t28 * s46 * t242
      t298 = t327 * t75 + s12 * (t327 * t139 + s16 * (s16 - t715 + t298 
     #+ t327) + t106 * (-t79 + t744)) - t966 * t237
      t337 = t609 * s46
      t51 = s12 * (s12 * (t107 + t882) + s16 * (t51 * t66 + s45 + t235) 
     #+ t106 * (-t143 * t66 + s56 - t235 - t372) + t292) - t237 * (s16 *
     # (t509 + t287) + t510)
      t48 = -(s25 - s12 - s16 - s26 + s34 + s45) * s56 + t48 - t594 + t5
     #95
      t107 = t956 + t235
      t143 = -t973 + t93
      t287 = -t17 * t313 - t361 + t496
      t292 = 17.0_dp * s56
      t340 = t136 * t149
      t344 = -t111 + t142
      t345 = 16.0_dp * s34
      t104 = t26 * t104 * t362
      t94 = s12 * (s12 * (s16 * (t330 - t564 + s45 - t490 - t528 - t547)
     # - s26 * (t545 + t276) + s46 * (-t700 + t505 - t367) - t124 * (t14
     #4 + t71) + t129 * t16 + t17 * t143 - t28 * (t677 - t267 - t45 + t1
     #13)) - s46 * (s45 * (-21.0_dp * s26 + t309) - 11.0_dp * t93) - s56
     # * (s26 * (-s45 * t556 - t345 - t354) - t470 + t158 * t3) + t16 * 
     #(t379 + t472 - t687 - t171 - t112) - t17 * (s26 * (-t492 + t149 + 
     #t158) + s46 * (t193 + t167) - t114 + t147 - t382) - t28 * (s34 * (
     #t695 + t113) + t145 + t282 + t383 + t384 + t73 + t340) + t66 * s46
     # * t344 + 13.0_dp * t265 - 9.0_dp * s26 * t669 + 9.0_dp * s56 * (t
     #172 - t111 + t113) - s16 * (s16 * (t17 * t319 + t235 - t257 + t392
     # + t495) + s34 * (-t343 - t355) + s45 * (-t292 - t374) - s46 * (t5
     #44 + t125) - t28 * t94 + 14.0_dp * t111 + s25 * t287 + s26 * (t303
     # + t276 - t371 - t272)) + t72) + t104
      t104 = t66 * t552
      t265 = 12.0_dp * s56
      t272 = s25 * t287 / 2.0_dp
      t276 = t179 * s46
      t287 = t503 * t167
      t303 = t322 * s34
      t322 = t107 * s12
      t120 = s12 * (s12 * (t124 * t463 + t28 * t973 + t66 * (-t93 - t279
     # + t144 + t71) + t118 - t120 + s16 * (t856 + t354 + t496 + t125) +
     # t113 - t267 - t45 + t677 + t322) + t124 * (s56 * t226 + t687) + t
     #17 * (s26 * (t522 - t113) + t112 - t380 - t501) + t28 * (s25 * (t1
     #21 + t111) + s45 * (-t973 - t113) + s46 * t180 + t145 - t280 + t38
     #6) - t66 * (s26 * (t267 - t222 - t71) + t739) - t173 + s16 * (s16 
     #* (t495 - t104 + t496 + t327 - t372) + s46 * (-t260 - t523) + t66 
     #* t71 - 11.0_dp * t129 + t514 + t63 + s26 * (-t894 + t265 + t199 -
     # t573) - t122 - t167 + t272) + t216 + t221 + t73 - t276 - t287 + t
     #303 + 7.0_dp * t318 - 7.0_dp * t170 + t222 * (-t376 - t505) + t258
     # * s34) - t237 * t513
      t121 = epinv * t7
      t170 = t124 * t123
      t173 = t391 * s56
      t180 = t184 * t139
      t258 = -t484
      t318 = t487 * t963
      t54 = s12 * (s12 * (t175 + s56 * (t523 - t361) - t124 * t1010 - t2
     #8 * t763 + t66 * t205 + t118 + s16 * (t354 + t116 - s34 + t152) + 
     #t113 - t158 - t45 - 7.0_dp * t197) - t124 * t161 * t319 - t16 * (s
     #56 * (t115 + t142) + t379) + t17 * (t111 * t43 + t112 - t163 - t38
     #0) - t28 * (s46 * (-t90 + t71) - t216) - t3 * (s46 * (-t530 - t149
     #) + s56 * (-t149 + t161)) + t66 * (t445 + t282 - t388 + t180 + t14
     #7) - s16 * (-s16 * (t170 - t483 + t127) - s34 * (-t173 - t89) - t1
     #24 * t159 + t16 * t167 - t66 * (-t113 + t71) + t378 - t514 - t954 
     #- s26 * (t312 + s45 - t255 + t546) - t134) - s45 * (s45 * t54 - t1
     #49) + t73 + 7.0_dp * s56 * t499 + 7.0_dp * t325 + t93 * (-t392 - t
     #292 - t547)) - t237 * (t318 + s16 * (-t520 * s16 - t28 * (-t119 **
     # 2.0_dp + t167) - s34 * (-s34 + t485 + s45) + t200 * t258))
      t116 = t271 * t66 + s46
      t127 = t192 * t79
      t90 = (t193 + t90 + t167 - t142 - t71) * s46
      t134 = (s25 * t271 + t149 + t158 + t168) * s26
      t147 = (-t557 - t265) * s45
      t123 = t66 * t123
      t159 = t237 * t231
      t61 = -(s12 * (s12 * t956 - t17 * t244 + t28 * (-t672 - t317 + t14
     #4 - t168 - t71) - t66 * t275 - s16 * (-t1078 * t28 + t152 - t166) 
     #- s34 * (-s25 + s26 - s45 - s46) - t113 + t45) + t124 * (s56 * (-t
     #179 + t93 + t279) - t186) + t28 * (s26 * (-t262 - t240 + t172 - t1
     #49 + t315 - t113) + s45 * t736 + t281 + t386 - t472) + t66 * (s34 
     #* (t168 + t71) + s45 * t638 - t324 + t471) - s16 * (s16 * (t569 + 
     #t123 + t234) + s56 * (t354 + t495) - t17 * t344 + t28 * (t518 + t1
     #13) - t529 + s26 * (t544 - s34 + t152 + t528) - t115 + t272) - s34
     # * (s45 * t358 - t111 + t113) - t340 - t380 - t500 - t72 - t73 - t
     #179 * (t354 + t507)) * s12 + t159 * t61 * s34
      t72 = s25 - s34 + s45
      t175 = s25 + s26 - s34 + s56
      t205 = t137 * s46
      t123 = s12 * (s12 * (t296 + t28 * (s26 * t72 + t111 + t184 + t689 
     #+ t71) - t389 * t66 + t202 + t494 + s16 * (t201 + t133 - s34 + t52
     #3) + t113 - t332 - t45) - t124 * (-t93 * t154 - t279 * t652) - t16
     # * s34 * (s46 * (s26 - s34) - t142) - t17 * s56 * (s34 * (s56 - s4
     #5) - t111 + t161) + t28 * (s26 * (t240 - t347 + t192 - t167 + t71 
     #+ t184) + s45 * t502 - t280 - t380 - t383) - t66 * (t379 + t739 + 
     #t163 - t186) + s16 * (s16 * (t175 * t28 + t123) + s34 * (-t354 - t
     #342 + t330) + t17 * (-t167 + t111) + t28 * (t273 - t158) + s26 * (
     #t59 + t307 - t469) + t272 + t205 - t279 * t124) + s45 * (t522 - t1
     #11) + s46 * (s56 * t521 - t167) - t114 + t264 + t531 + t73 - t309 
     #* t179) - t237 * (s16 * (t563 * t258 + t486) + t318)
      t163 = t16 - t430
      t186 = s25 * (-s34 * t697 - t17 * t228 - t28 * (epinv * t313 + s56
     #))
      t202 = t279 * epinv
      t228 = s34 * (s56 * (-t617 - t408) - t167) + t180 - t383 + t737 + 
     #t161 * (-s45 + t411)
      t234 = t733 * t139
      t258 = s26 * (epinv * (-s45 * t72 + t184 + t71) + s26 * (-t407 - t
     #400) - t333) + t113 * (t617 - t408) + t127 - t548 - t411 * t167 + 
     #t234 - t282 - t453 * t79 + t382 + t321
      t262 = t404 * (s46 * t381 + t142)
      t156 = -t124 * (s26 * (-s34 * t214 + t202) - s34 * (t149 + t202) +
     # s56 * (-s46 * t740 + t149) + t111 * t979) + t16 * t262 - t17 * t2
     #28 + t28 * t258 + t306 * t93 * t730 + t66 * (s45 * (t424 + t267) +
     # s46 * (-t422 - t179 + t623 - t149) - t379 * epinv) + s16 * (s16 *
     # (t236 * t66 - t28 * (t2 * t232 + t404 + t617) + t372) - t124 * t1
     #56 * s45 + s34 * (s56 * (t3 - t757) - t309 - 7.0_dp * t616) - t17 
     #* (t442 - t415 - t197) - t28 * (s56 * (t617 + t408) + t421) + t66 
     #* (-t122 + t459) + s26 * ((-t705 + 9.0_dp) * s34 + t16 * t618 - t5
     #44 + t59 * t29) + s46 * (t414 + t411) + t186) + s45 * (s34 * (-s56
     # + t408) - t111 * t457 + t459) + s46 * (t632 + t424 + t415) + t113
     # * (t400 - t613) + t208 + t462 + t267 * (t622 + s26 + t408) + 7.0_
     #dp * t279 * t148 + 9.0_dp * s34 * (t179 + t158)
      t214 = t122 + t158
      t228 = t168 - t113
      t52 = s13 * t52
      t236 = t93 * t79
      t258 = t111 - t113
      t262 = -t172 + t222
      t264 = (-t167 + t71) * s46
      t271 = 13.0_dp * s26
      t273 = t271 + t292
      t280 = (-t617 + t404 - t400) * s45
      t281 = t3 * s26 * (t266 + t714)
      t282 = (t617 + t404) * s56
      t313 = s46 * (t424 - t415) + t385 + t462 + t738 + t761
      t315 = s56 * (-t527 + t1077) - t207 * epinv
      t318 = 12.0_dp * t179 * (s45 * t713 - t404)
      t319 = 11.0_dp * s26 * (t279 * t713 + t93)
      t183 = s12 * (-t235 + t183 - s16 - s26 + s34 - s56) + t124 * (s26 
     #* (t447 - t400) - t432) + t16 * t331 + t17 * (s56 * (-t404 - t414)
     # - t144) - t28 * t211 + t66 * (s25 * (-t401 - s46) - t149 - t424 +
     # t664 - t756) + t1097 + s16 * (epinv * (t59 + t199) + t196 - t250 
     #+ 7.0_dp * t341 - t545) + s25 * t503 + t39 - t425 - t280 + t281
      t128 = s12 * t183 + t315 * s34 + t124 * (s45 * (-t643 + t149 + t28
     #2) + t501 - t197 * t217) - t16 * t404 * t138 - t17 * (epinv * t164
     # + s46 * (-t275 - t294) + t112 + t145 - t321 + t174 * t242) + t28 
     #* t313 - t306 * t93 * (s34 + s45 - s46) - t66 * (s34 * (s34 * (-s2
     #5 - t440) - t149 + t168 - t408 * t232) + s46 * (epinv * (-t193 - t
     #111) + t111) - t234 + t765) + 7.0_dp * s34 * (t174 - t142) + 7.0_d
     #p * t111 * (-t617 - t404 + t446) + 7.0_dp * t279 * (t617 - s34) + 
     #s16 * (s16 * (epinv * t128 + t257 - t310 - t311 - t523 + t573) + s
     #56 * (t345 + t376) - t17 * (-t209 + t71) - t306 * (s34 * t715 - t2
     #79) - t66 * t498 - 12.0_dp * t111 + epinv * t135 + s26 * (s34 * t5
     #56 - t173 - t354 + t375) + s56 * t213 - t167) + t167 * (-t617 - t4
     #00 + t550) + t208 + t192 * t148 + t122 * t452 + t459 * s45 + t318 
     #+ t319
      t135 = t400 * t3
      t138 = t16 * (t617 + s46)
      t164 = t192 * t2
      t183 = -t442 + t459
      t196 = t200 * t2
      t209 = t167 * t457
      t211 = 9.0_dp * s45
      t232 = t17 * t425
      t234 = t66 * t424
      t242 = (t496 + t327) * t2
      t250 = t124 * (-t617 + s45)
      t215 = -t66 * t215
      t257 = (-t416 + t415) * s56
      t275 = t400 * t167
      t294 = t2 * t145
      t310 = s25 * (s26 * (t400 - s45) + s56 * (t400 - t613 - t446) + t2
     #91 + t760) + s46 * (-t424 - t142 + t415) + t236 - t462 - t761 + t2
     #94
      t311 = t446 * t111
      t313 = (-t442 + t664) * s56
      t315 = t174 * epinv
      t314 = epinv * t150 + s34 * (t616 * t314 + t172) + s46 * (t315 - t
     #420) + t466 + t222 * (-t447 - t414)
      t324 = t404 - t446
      t325 = 12.0_dp * t417 * s45
      t241 = -s12 * (-t28 * t708 - t194) - t124 * t183 - t17 * (s34 * t2
     #41 - t420) - t28 * (s26 * (-t407 + s34) - t184 - t164) - t66 * (s4
     #6 * (t447 + t446) + t142 - t423 - t424) - t281 + 7.0_dp * t405 - s
     #16 * (epinv * (t354 - t260) - t28 * (t407 - s34) + s45 - t135 - t1
     #38) + s45 * (-t617 - t447 - t400 - s34) - t267 - t39 + t425 + t435
     # * t93
      t126 = s12 * (s12 * t241 + t124 * (-s34 * t249 + t324 * t71) + t17
     # * t314 + t28 * t310 + t306 * t404 * t262 - t66 * (s26 * (-t424 + 
     #t267) + s34 * (-t338 - t197) - t769 - t404 * t174 - t404 * t167) +
     # s16 * (s16 * (s34 * (-t17 + t435) + s46 + t242 - t250 - t215) + s
     #46 * (-t126 - t211) + t28 * (s26 * (-t163 * s45 - s46 * t627 + t61
     #7 * t16) + t522) + t306 * (s45 * (-t617 - t404) + t442) + 11.0_dp 
     #* t616 * t79 + s25 * (-t28 * (s46 * t413 + s45 - t617) + t263 + t1
     #96) - t420 + t209 + t93 * (t243 - 9.0_dp) - t232 + t234 + t404 * (
     #t373 - t255) - t142 * t3) + s34 * (epinv * t228 + t113 + t167 + t1
     #74) + t149 * (t712 + s45) - t208 + t257 + t275 - 11.0_dp * t664 * 
     #s26 - 9.0_dp * s34 * (t256 + t202) + 7.0_dp * t459 * t552 - 7.0_dp
     # * t311 - 7.0_dp * t313 - 14.0_dp * t756 * t79 - t325 + t93 * (epi
     #nv * t273 - t199)) - t270 * t406
      t202 = t306 * (t442 - t419)
      t241 = s25 * (-t442 + t620) + s45 * (-t643 + t425 + t164 - t453) +
     # s46 * t641 + t167 * (-t414 + s25) - t294 + t267 * (t617 - t616)
      t135 = s12 * (-t28 * t414 + s34 - t393 - t617) + s34 * (t251 + t12
     #5) + t124 * (-t442 + t709) + t17 * (-t422 - t149 - t420) + t28 * s
     #25 * (t616 + t206 + t266) - t66 * (s45 * t646 + t184 + t423 + t424
     # + t756) + t281 + s16 * (-t28 * t711 + s45 - t135 - t138 - 7.0_dp 
     #* t710) - s45 * (-t617 + t404 - t400 + s25) + t39 - t425 - t624
      t59 = s12 * (s12 * t135 + t124 * (t459 * t119 + t446 * t71 + t150 
     #+ t501) - t16 * t158 * t433 - t17 * (s26 * (-t421 - t425 - t459 + 
     #t184) + s46 * (-t1143 - t174) - t385 + t466 + t122 * t252 + t404 *
     # t71) + t28 * t241 - t3 * s34 * (t142 + t218) + t66 * (s46 * (-t44
     #2 - t422) - t765 + t267 * (s45 - t440)) - 13.0_dp * t236 + t318 + 
     #t319 + s16 * (s16 * (s34 * (t306 - t690) - s46 + t215 - t242 + t25
     #0) - 11.0_dp * epinv * t129 + s34 * (t404 * t28 + t173 + t199) + t
     #232 - t234 - 13.0_dp * t522 + s26 * ((21.0_dp - t729) * s34 + s46 
     #* (-t17 + t461) - 12.0_dp * t617 + t59 * t163) + t186 - t209 + t42
     #0 - 9.0_dp * s46 * t694 - t202) + s34 * (t315 - t425 + t168 + t760
     #) + s56 * (t416 - t415) + t149 * (t712 + t408) + t208 - t275 + 7.0
     #_dp * s34 * (s46 * t754 - t623) + 7.0_dp * t311 + 7.0_dp * t313 + 
     #9.0_dp * t158 * t217 + t345 * t179) - t270 * t397
      t135 = -t404 - t407
      t138 = t404 * t79
      t173 = t158 * t213
      t139 = (s26 * t141 - s56 * t139) * s25
      t136 = s34 * (s26 * (t271 + t544 + t292) - t172 * t306) - t124 * (
     #s34 * (t623 - t167) + t71 * (-t446 - s34) + t197 * t740) - t16 * (
     #s46 * (t756 + t193) - t150) + t17 * (t113 * t707 + t321 - t466 + t
     #717 + t739) + t28 * (epinv * (t134 + t90 - t127) + t114 + t264 - t
     #294 + t382 + t139) + t66 * (s46 * t183 + t233 * t71 + t647 + t267 
     #* (s34 - t440)) - 16.0_dp * t236 + t325 - s16 * (-s16 * (s34 * (12
     #.0_dp - t690) + t124 * (-t617 + t446) - t242 - t414) - s34 * (s46 
     #* (7.0_dp - t430) - t438) - t17 * (-t400 * t116 + t425) - t28 * t1
     #38 + t202 + t234 - 9.0_dp * t664 - s25 * (t28 * (-t617 + t616 - t6
     #21) - t196 - t199 - t238) + t415 - t420 - t93 * (25.0_dp - t757) +
     # 16.0_dp * t951) - s34 * (-t315 + t425 + t168 - t760) + t208 - t25
     #7 - t275 + t459 * t136 + 7.0_dp * t313 + 7.0_dp * t311 + 7.0_dp * 
     #t173 + 11.0_dp * s45 * (t222 * t713 + t149) + 9.0_dp * t142 * t441
      t85 = -epinv * t61 + s12 * (s12 * (t124 * t214 + t149 * t16 - t17 
     #* (-s56 * t268 - t161 - t184) + t28 * (t240 + t85 + t192 + t71) + 
     #t207 - 7.0_dp * t339 + s16 * (-7.0_dp * t165 + t201 + t152) + t113
     # + t267 - t45) - t124 * s34 * (t167 + t71) + t16 * (s56 * (t149 - 
     #t167) - t150 + t180 + t161 * t178) - t17 * (-s56 * t258 + t321) - 
     #t28 * (t472 + t382 + t139 - t210 + t114) + t66 * (s45 * (t111 - t7
     #1) - t216 + t276 - t388) + s16 * (s16 * (t201 + t170 - t361) + s34
     # * (-t557 + t371 - t528) + t16 * s45 * t178 - t28 * t639 + t66 * t
     #205 + t132 - t378 + s25 * (s34 + t200 + t199 + t999) + s26 * (t527
     # - t371 + t125 + t547)) + s34 * (t168 + t113) - t52 - 7.0_dp * t11
     #1 * t165 - 7.0_dp * t173 + 7.0_dp * t687 - 11.0_dp * s45 * (t149 +
     # t222) - 11.0_dp * t179 * s34 - t93 * (-t371 + t547) - t122 * (t21
     #1 + t505))
      t132 = s12 * (-t235 - t261 + t289) + t1002 * t237 + t609
      t139 = t26 * t2
      t141 = s12 - s16 - s26 - s45
      t161 = t28 * t259
      t168 = t166 * s56
      t140 = t28 * (s34 * t237 - t172) + s16 * (t140 + s16) + s26 * (t13
     #7 - t80 + s26) + t168
      t170 = s12 * (s16 * t175 - t106 * t652 + t193 + t227) + t159 * s34
      t106 = s12 * (t458 + t106 * (-t679 - t483) + t117 + t193) - t237 *
     # (t370 - t491)
      t4 = -t197 * t28 + t244 + t348 + t4 - t680
      t117 = t237 * (t91 + t244 + t151 + t370 - t653)
      t151 = -t414 + s34 - t411
      t173 = t66 * t324
      t175 = (t296 + t4) * s12
      t8 = t22 * (t1015 * (t593 * (t187 - t56) + t288 * s34 * t777 * t10
     #13) + t479 * t593 * t7 * t1013) + t1013 * (t1015 * (t189 * (t182 *
     # t229 - t195 * t53 - t297 * t54) + t7 * (t436 * t1014 * t593 * s34
     # - t177 * t237) - t927 * (-t337 + epinv * t51 + s12 * (t336 + t326
     # + s25 * (t148 - t235 - t261) - t244 + t320) - t237 * (s25 * (-t50
     #9 + s26 - t261) + t157 + t316 + t639 - t984 - t185)) + t936 * (epi
     #nv * t181 + s12 * (t296 - s25 * (t744 - s46 - t261) - s46 * t323 -
     # t290) + t237 * (s25 * (t65 - t119 - t261) + t290 - t291))) + t78 
     #* (s12 * (t41 - t80 + s16 + s26 + s56 - s46) - t237 * t81)) * t19 
     #+ t1013 * (t1015 * (-t450 * t35 * t231 + (t8 - t190) * t189 * t182
     #) + t259 * t796 - t74 * (s12 * (-t426 - t278 + t622 + t621 - t411 
     #- s12) + t42 * t12)) * t9
      t42 = t346 * (t170 * t941 + t181 * t823 - t189 * (t120 * t225 + t1
     #23 * t295 + t308 * (s12 * t128 + t253 * t540) + t456 * t59 + t650 
     #* (s12 * (s12 * (-s12 * (t414 - s34) + t17 * (s34 * (-t616 + s26) 
     #+ s56 * (-t400 - s46) - t149) + t28 * (s26 * (epinv * t72 - s25) +
     # s56 * t724 - t424 - t442 + t453) + t66 * (s46 * (-s25 + s34 + t41
     #1) + t555 + t677) - s16 * (-t17 * t646 - t28 * (-t617 - t407 - t40
     #0) + t446 * t66 + t404) - t174 + t39 - t425 - t644 - t222 * t163) 
     #+ t156) + t253 * t658)) - t51 * t825 + t674 * (epinv * t170 - t28 
     #* (t175 - t117)) + t877 * t237 * t7 + t873 * s34 * (s12 * (s12 - t
     #235 - t261 + t426 + t278 + t429) + t237 * t875)) * t19
      t51 = t212 * t1017
      t13 = t1022 * propW34 * (t1023 * (t9 * (t212 * t13 * (s12 - s45 - 
     #s46 - s56) * t231 - t402 + t20) + t21 + t51 * (t22 * (t231 * t776 
     #+ t304 * (epinv * t231 - s12 - s16 - s26 + s34 + s45)) + t231 * (t
     #131 - t775) * t9) * t237 * s34 + t188 * (t1015 * (t106 * t814 + t9
     #26 * (epinv * t106 + t117 - t175)) + t797 * (s12 * (t41 + t254) - 
     #t26)) * t19) + propW16 * t51 * t785 * t593 * t7 * (s16 * (t2 * (t7
     #83 - t903) - t861 + t862 - t902 - t98) + t105 + t993))
      t20 = t212 * t19 * t189 * (t230 * t61 - t285 * (s12 * (s12 * (-t32
     #2 + t143 * t28 + t17 * t219 - t66 * (-t338 - t279 + t144) - t118 -
     # t377 + s16 * (-t856 + s45 - t496 - t307) - t113 + t267 + t45 - t6
     #77) - s26 * (s34 * (-t342 - t199) + t147) - t124 * (s45 * t691 - t
     #176) + t17 * (s56 * (-t158 + t113) - t112 + t388) - t28 * (t134 + 
     #t90 - t127 + t145) + t66 * (-s46 * t1003 + t154 * t71 - t171 + t37
     #9) - s16 * (s16 * (-t50 + s46 + t496 + t327 - t372) + s56 * (-t505
     # + t538) + t17 * (s26 * t116 - t113) + t28 * t522 - t122 + t167 - 
     #t246 + t272 + t63 - 9.0_dp * t93 - 9.0_dp * t279) - t216 - t221 + 
     #t276 + t287 - t303 - t73 - 7.0_dp * s56 * t249 + 7.0_dp * t500) + 
     #t253 * t966))
      t21 = t612 * t62 + t99 * (s12 * (s12 * (t353 + s34 * (-t545 + t152
     # - t125) + t28 * t317 + t66 * (t130 + t222 + t144) + s16 * (-t260 
     #+ t496 + t125 + t43) + t174 - t302 + t328 + t333 - t45 + t49 + t17
     #9 * t3) + s34 * (s56 * (t527 - t564) - 11.0_dp * t111) + t28 * (s4
     #6 * (t193 - t167 - t279) + t145 - t269) - t66 * (s26 * t77 + t87) 
     #+ s25 * (-t224 - t149 + t93 - t71) + s46 * (-t64 - t659 + t111) + 
     #t471 - t548 + t67 + t73 + 12.0_dp * s26 * (s45 * (s34 - s56) + t14
     #9) - 7.0_dp * s45 * t1003 + 7.0_dp * s56 * t1003 - t108 + t109) + 
     #t100 * t237)
      t4 = t1015 * (-t1008 * (epinv * t132 + t161) - t15 * t30 + t189 * 
     #(t146 * t443 - t526 * (t245 * t110 + s12 * (s12 * (s26 * (17.0_dp 
     #* s46 + t898) + s34 * (-t496 - t547) + s56 * (t312 + t649) - t16 *
     # (t204 + t279) + t28 * (t677 - t45 + t113) + t306 * t169 - t248 + 
     #7.0_dp * t71 + s16 * (s16 + t89 + t1039 + t497 - s34 + t363)) + t1
     #60) - t237 * t511)) + t298 * t959 - t34 * (epinv * t30 + t161) + t
     #915 * (s12 * (t121 * t231 + t296 + t4) - t117) + t981 * (-t337 + e
     #pinv * t298 + s12 * (t336 + t326 + s25 * (t148 - t671 - t68) - t24
     #4 + t320) - t237 * (t316 + t639 + t157 - t984 + t746 - t185)) - t1
     #025 * t904 * t7 * t69 - t36 * s34 * (s12 * (s12 - t80 - t414 + t42
     #6 + t278 + t411) - t26)) * t19
      t15 = t55 * t593 * (t482 * t909 * t48 + t1020 * t7)
      t4 = t9 * (t1015 * (t10 * t30 + t1007 * (s12 * (t101 * t716 - t28 
     #* t433 + t446 - t538 + t613 + t616) + t139)) + (-t1015 * t876 + t8
     #4) * t7 * t237) + t1015 * (epinv * (t259 * t668 + t286 * t32) + t1
     #89 * t21 + t948 * (s12 * (t353 - t198 * t66 + t1046 + s16 * (t65 +
     # t137 - t261 + s16) + s26 * (t137 - t261 + s26) - t279 + t168) - t
     #237 * t517) - t969 * (s12 * (-s12 * t101 + s46 * (-t141 + s46) + s
     #56 * (-t141 + t235 + s56) + t802) - t121 * t365) + t970 * (s12 * (
     #s12 * (-t28 * t301 + t978) + t113 * t28 + t138 * t17 - t66 * (-s46
     # * (t404 + s45) + s56 * (-t446 - s46)) - s16 * (t614 + t299 + t151
     # + t173) + s26 * (-t151 - t173) + s46 * t731 - s56 * (t616 + t622)
     # - t442 + t71) - t102 * t237) + t293 * (t41 * t593 + t24)) * t11 +
     # t4 + t15
      t10 = t23 * (-t900 - t159) + t237 * (t1005 * t283 + t7 * t992) + t
     #33 * (s12 * (s12 * t716 - t628 * t66 + t671 + t101 * (-t17 + t398)
     # - t617) + t139) - t997 * (-s12 * (t28 * t956 - t715) + t26) + (t2
     #45 * t107 - t94) * t189 * (-t1027 + t356)
      t15 = t203 * (s12 * (s12 * (s12 * (t956 - t1073) + t17 * t214 + t2
     #8 * t773 + t66 * (t334 + t144 + t71) + t118 - 7.0_dp * t331 + t49 
     #+ s16 * (t354 - t260 + t496 + t327 + t125) + s45 * (s25 - s26 - s5
     #6) + t113 - t45 - t93 * t306) - s26 * (s34 * t273 - t147) - t124 *
     # t71 * t79 + t17 * (s26 * t258 + t122 * t166 - t150 + t180) - t28 
     #* (s25 * (t172 - t179 + t167 - t111 + t142) + t114 - t145 + t382 +
     # t264) - t306 * s34 * t262 - t66 * (-s26 * t638 + t87) - s16 * (-s
     #16 * (t495 - t104 + t496 + t327 - t309) + s56 * (t345 - t538) - t3
     #06 * (t247 + t111) - 12.0_dp * t149 + t858 + s25 * (-t17 * t223 + 
     #s34 - t46) + s26 * (t239 + t894 - t265 - t199) + t122 + t167 + 11.
     #0_dp * t115) - s34 * t228 - t155 - t216 + t384 - t52 + 14.0_dp * t
     #236 + 7.0_dp * s46 * t226 + 7.0_dp * t477 - 7.0_dp * t742 + 9.0_dp
     # * s34 * t129) + t270 * t335) - t394 * t126
      t15 = t15 * t189 + t259 * t991 + t286 * t31 + t82 * t7 + t914 * t2
     #37 * (-epinv * t283 - s12 + s34 + s45 + t103 + t448)
      t6 = (t22 * (-t48 * t673 + t533 * (s25 * (-t801 + t121) + t562 - t
     #595 - t988)) + epinv2 * t6 * t7) * t1014 * t593
      t3 = t346 * t630 * (t11 * t15 + t19 * (-t1000 * t132 + t189 * (t22
     #0 * (s12 * (s12 * (s12 * (t28 * t646 - t194) + s34 * (t392 + t367)
     # + t124 * (-t442 - t158) - t16 * t284 + t17 * s46 * (-t404 - s25 +
     # t266) + t28 * (s25 * (t616 - s45) + s26 * t135 - t164) - t66 * (t
     #423 + t282 - t664) + t281 - s16 * (-t135 * t28 - t16 * t519 + t400
     # * t3 - t309 - t446) - t280 + t39 - t425 - t267 * t457) + t136) - 
     #t270 * t553) + t467 * (t28 * (t165 * t76 + t270 * (s16 * (s34 - s4
     #5 - s56 - s46) + s56 * t300 - t149 + t633 - t71 - t91 + t95)) - t8
     #5)) + t237 * (-t1006 * t283 - t910 * (-epinv * t946 + t28 * (t429 
     #- s26 + s34 + s45) + t1012 - t913)) + t7 * t1014 * (t69 * t863 + t
     #40)) + t9 * t10 + t6)
      t3 = t3 + t630 * t1017 * t4 + t274 * (t1017 * (-t92 * t793 * t593 
     #+ t1013 * (t237 * (t19 * t83 + t9 * (epinv * t1015 * (-t778 + t899
     #) + t798)) + (t22 * (t1015 * t277 + t2 * (t1015 * t1019 - t866) - 
     #t1015 * t1018 * t29) + t37 * t191) * t1014 * t593) * t7 + t20) + t
     #212 * (t11 * (-t2 * t879 + t809 + t989) + t22 * (-s34 * t795 + t7 
     #* t86 + t58)) * t231)
      t4 = t274 * t97 * (-t28 * (t22 * (-t2 * (t927 + t936 + t308 + t443
     # + t456) + t225 + t229 + t526 + t823 + t825) + t9 * (t2 * (t465 + 
     #t762 + t788 + t789 + t790) - t190 - t356 - t489 - t944 - t968)) + 
     #t11 * (-t2 * (t970 + t394) + t203 + t948) + t19 * (t959 + t285) - 
     #t2 * (t981 + t220) * t22) * struc11Step4 * t1015
      result = mytree * t18 + t1 * struc4Step5 + t5 * struc8Step5 + (-t4
     #7 * t13 + t60 * t3 + t70 * t96 * t593 * t7 * t88 + t97 * (t630 * (
     #-t329 * t593 * t7 * t1017 + (t1017 * (-t955 * (-s12 * t140 + t27 *
     # (t26 - t75)) - t977 * (s12 * (-s12 * (-t509 - t725 + t874) - epin
     #v * t140 - t38 * t305) + t360 * t12)) + t2 * t189 * (t14 * t880 + 
     #t25 * t881)) * t1015 * t11 + t42 + t346 * (s34 * (t162 * t872 + t7
     #02 * (epinv * t162 + t288)) + t231 * (t35 * t944 + t515 * (-t609 +
     # s12 * (t235 + t133 - t289) + t26)) + (-t44 + t489) * t189 * t120)
     # * t9) + t274 * (t1017 * t8 + t19 * t787 + t212 * (t189 * (-t14 * 
     #t837 - t25 * t698) + t231 * (-t153 + t57)) * t11))) * struc9Step5 
     #+ t4

           ampNonresonantLightFullMM = result
       end function ampNonresonantLightFullMM

       function ampNonresonantLightFullMP()
           implicit none
           complex(dp) :: ampNonresonantLightFullMP

           complex(dp) :: result

      result = 0.0_dp

           ampNonresonantLightFullMP = result
       end function ampNonresonantLightFullMP

       function ampNonresonantLightFullPM()
           implicit none
           complex(dp) :: ampNonresonantLightFullPM

           complex(dp) :: result

      result = 0.0_dp

           ampNonresonantLightFullPM = result
       end function ampNonresonantLightFullPM

       function ampNonresonantLightFullPP()
           implicit none
           complex(dp) :: ampNonresonantLightFullPP
           complex(dp) ::  t1,t10,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109
           complex(dp) ::  t11,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t12
           complex(dp) ::  t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t13,t130
           complex(dp) ::  t131,t132,t133,t134,t135,t136,t137,t138,t139,t14,t140,t141
           complex(dp) ::  t142,t143,t144,t145,t146,t147,t148,t149,t15,t150,t151,t152
           complex(dp) ::  t153,t154,t155,t156,t157,t158,t159,t16,t160,t161,t162,t163
           complex(dp) ::  t164,t165,t166,t167,t168,t169,t17,t170,t171,t172,t173,t174
           complex(dp) ::  t175,t176,t177,t178,t179,t18,t180,t181,t182,t183,t184,t185
           complex(dp) ::  t186,t187,t188,t189,t19,t190,t191,t192,t193,t194,t195,t196
           complex(dp) ::  t197,t198,t199,t2,t20,t200,t201,t202,t203,t204,t205,t206
           complex(dp) ::  t207,t208,t209,t21,t210,t211,t212,t213,t214,t215,t216,t217
           complex(dp) ::  t218,t219,t22,t220,t221,t222,t223,t224,t225,t226,t227,t228
           complex(dp) ::  t229,t23,t230,t231,t232,t233,t234,t235,t236,t237,t238,t239
           complex(dp) ::  t24,t240,t241,t242,t243,t244,t245,t246,t247,t248,t249,t25
           complex(dp) ::  t250,t251,t252,t253,t254,t255,t256,t257,t258,t259,t26,t260
           complex(dp) ::  t261,t262,t263,t264,t265,t266,t267,t268,t269,t27,t270,t271
           complex(dp) ::  t272,t273,t274,t275,t276,t277,t278,t279,t28,t280,t281,t282
           complex(dp) ::  t283,t284,t285,t286,t287,t288,t289,t29,t290,t291,t292,t293
           complex(dp) ::  t294,t295,t296,t297,t298,t299,t3,t30,t300,t301,t302,t303
           complex(dp) ::  t304,t305,t306,t307,t308,t309,t31,t310,t311,t312,t313,t314
           complex(dp) ::  t315,t316,t317,t318,t319,t32,t320,t321,t322,t323,t324,t325
           complex(dp) ::  t326,t327,t328,t329,t33,t330,t331,t332,t333,t334,t335,t336
           complex(dp) ::  t337,t338,t339,t34,t340,t341,t342,t343,t344,t345,t346,t347
           complex(dp) ::  t348,t349,t35,t350,t351,t352,t353,t354,t355,t356,t357,t358
           complex(dp) ::  t359,t36,t360,t361,t362,t363,t364,t365,t366,t367,t368,t369
           complex(dp) ::  t37,t370,t371,t372,t373,t374,t375,t376,t377,t378,t379,t38
           complex(dp) ::  t380,t381,t382,t383,t384,t385,t386,t387,t388,t389,t39,t390
           complex(dp) ::  t391,t392,t393,t394,t395,t396,t397,t398,t399,t4,t40,t400
           complex(dp) ::  t401,t402,t403,t404,t405,t406,t407,t408,t409,t41,t410,t411
           complex(dp) ::  t412,t413,t414,t415,t416,t417,t418,t419,t42,t420,t421,t422
           complex(dp) ::  t423,t424,t425,t426,t427,t428,t429,t43,t430,t431,t432,t433
           complex(dp) ::  t434,t435,t436,t437,t438,t439,t44,t440,t441,t442,t443,t444
           complex(dp) ::  t445,t446,t447,t448,t449,t45,t450,t451,t452,t453,t454,t455
           complex(dp) ::  t456,t457,t458,t459,t46,t460,t461,t462,t463,t464,t465,t466
           complex(dp) ::  t467,t468,t469,t47,t470,t471,t472,t473,t474,t475,t476,t477
           complex(dp) ::  t478,t479,t48,t480,t481,t482,t483,t484,t485,t486,t487,t488
           complex(dp) ::  t489,t49,t490,t491,t492,t493,t494,t495,t496,t497,t498,t499
           complex(dp) ::  t5,t50,t500,t501,t502,t503,t504,t505,t506,t507,t508,t509
           complex(dp) ::  t51,t510,t511,t512,t513,t514,t515,t516,t517,t518,t519,t52
           complex(dp) ::  t520,t521,t522,t523,t524,t525,t526,t527,t528,t529,t53,t530
           complex(dp) ::  t531,t532,t533,t534,t535,t536,t537,t538,t539,t54,t540,t541
           complex(dp) ::  t542,t543,t544,t545,t546,t547,t548,t549,t55,t550,t551,t552
           complex(dp) ::  t553,t554,t555,t556,t557,t558,t559,t56,t560,t561,t562,t563
           complex(dp) ::  t564,t565,t566,t567,t568,t569,t57,t570,t571,t572,t573,t574
           complex(dp) ::  t575,t576,t577,t578,t579,t58,t580,t581,t582,t583,t584,t585
           complex(dp) ::  t586,t587,t588,t589,t59,t590,t591,t592,t593,t594,t595,t596
           complex(dp) ::  t597,t598,t599,t6,t60,t600,t601,t602,t603,t604,t605,t606
           complex(dp) ::  t607,t608,t609,t61,t610,t611,t612,t613,t614,t615,t616,t617
           complex(dp) ::  t618,t619,t62,t620,t621,t622,t623,t624,t625,t626,t627,t628
           complex(dp) ::  t629,t63,t630,t631,t632,t633,t634,t635,t636,t637,t638,t639
           complex(dp) ::  t64,t640,t641,t642,t643,t644,t645,t646,t647,t648,t649,t65
           complex(dp) ::  t66,t67,t68,t69,t7,t70,t71,t72,t73,t74,t75,t76
           complex(dp) ::  t77,t78,t79,t8,t80,t81,t82,t83,t84,t85,t86,t87
           complex(dp) ::  t88,t89,t9,t90,t91,t92,t93,t94,t95,t96,t97,t98
           complex(dp) ::  t99

           complex(dp) :: result

      t1 = intHL0s25s26s34s56x1120D2eps0()
      t2 = 3.0_dp
      t3 = 2.0_dp
      t4 = t2 * s25
      t5 = t4 * propZ25
      t6 = -t3 + t5
      t7 = gw ** 2.0_dp
      t8 = gb ** 2.0_dp
      t9 = propZ25 * s25
      t10 = t9 * t8
      t11 = t7 * t6
      t12 = t11 + t10
      t13 = intHL0s25s26s34s56x1130D4eps0()
      t14 = intHL0s25s26s34s56x1210D2eps0()
      t15 = intHL0s25s26s34s56x1220D4eps0()
      t16 = intHL0s25s26s34s56x1310D4eps0()
      t17 = intHL0s25s26s34s56x1130D4eps1()
      t18 = 1.0_dp - epinv
      t19 = intHL0s25s26s34s56x1220D4eps1()
      t20 = intHL0s25s26s34s56x1310D4eps1()
      t21 = intHL0s25s26s34s56x1120D2eps1()
      t22 = intHL0s25s26s34s56x1210D2eps1()
      t23 = intHs16s25s26s34s56x1112D2eps0()
      t24 = 4.0_dp
      t25 = t24 - t5
      t26 = t7 * t25 + t10
      t27 = s25 + s26 + s56
      t28 = intHs16s25s26s34s56x1113D4eps1()
      t29 = intHs16s25s26s34s56x1114D6eps0()
      t30 = intHs16s25s26s34s56x1123D6eps0()
      t31 = intHs16s25s26s34s56x1213D6eps0()
      t32 = intHs16s25s26s34s56x1114D6eps1()
      t33 = intHs16s25s26s34s56x1123D6eps1()
      t34 = intHs16s25s26s34s56x1213D6eps1()
      t35 = intHs16s25s26s34s56x1112D2eps1()
      t36 = epinv * t35
      t37 = (t21 + t22) * t18
      t38 = (t23 + t36) * t26
      t39 = t17 + t19 + t20
      t40 = -t18 * t39 + t13 + t15 + t16
      t41 = (-t18 * (t33 + t34) + t28 + t30 + t31) * t26
      t42 = t32 * t18
      t43 = -t42 + t29
      t44 = t26 * t43
      t45 = 6.0_dp
      t46 = intHLs160000x0111D0eps0()
      t47 = intHLs16s25s26s34s56x1112D4eps0()
      t48 = intHs16s25s26s34s56x1110D2eps0()
      t49 = s16 + s26 - s34 + s56
      t50 = intHLs160000x0111D0eps1()
      t51 = intHLs16s25s26s34s56x1111D2eps0()
      t52 = intHs16s25s26s34s56x1212D4eps1()
      t53 = t7 * t25 + t10
      t54 = t3 * s13
      t55 = t54 + s12
      t56 = s13 + s16
      t57 = t24 * s12
      t58 = s13 + s25
      t59 = -s12 + s16
      t60 = 5.0_dp
      t61 = s25 - s12
      t62 = t24 * s16
      t63 = t60 * s26
      t64 = s25 - s26
      t65 = s16 - s34
      t66 = t2 * s12
      t67 = t3 * t65
      t68 = s13 - s46
      t69 = s16 + s26 - s34
      t70 = t3 * s12
      t71 = t24 * s13
      t72 = t2 * s16
      t73 = 10.0_dp
      t74 = t2 * t58
      t75 = s12 + s13
      t76 = 7.0_dp
      t77 = t76 * s25
      t78 = t2 * t75
      t79 = t3 * s16
      t80 = t73 * s26
      t81 = t2 * s46
      t82 = t76 * s12
      t83 = (s13 + t57) * s16
      t84 = s12 - s13
      t85 = s25 ** 2.0_dp
      t86 = s25 * t85
      t87 = s56 ** 2.0_dp
      t88 = s12 * t84
      t89 = s26 ** 2.0_dp
      t90 = s34 ** 2.0_dp
      t91 = t2 * t90
      t92 = t24 * (t88 + t87 - t85) + t60 * t89 + (t70 - t72 + t71) * s2
     #5 + (s12 * t73 + s16 - t74) * s26 + (t79 + s26 - t78 + t77) * s34 
     #+ (-t3 * t64 + s34 + t57) * s46 + (t81 + t82 + t79 - s34 + t80 - t
     #71) * s56 + t83 - t91
      t93 = t55 * s16
      t94 = s16 ** 2.0_dp
      t95 = s16 * t94
      t96 = t2 * t88
      t97 = s25 - s13 + s16 - s34 + s46 + s56 + s12
      t98 = t3 * s26
      t99 = t97 + t98
      t100 = t2 * epinv
      t101 = t3 - t100
      t102 = s25 - s13 + s16 - s34 + s46 + s56
      t103 = s34 + s56
      t104 = -t101
      t105 = t104 * s13
      t106 = epinv * s25
      t107 = t18 * s16
      t108 = epinv * s34
      t109 = epinv * s56
      t110 = epinv * s46
      t111 = t2 * (-t110 + s56)
      t112 = t24 * (-t109 + s26)
      t113 = t45 * s26
      t114 = t113 * epinv
      t115 = t24 * s25
      t116 = t2 * s34
      t117 = -s26 + s34
      t118 = s46 + s56
      t119 = s34 * s46
      t120 = t118 * s56
      t121 = s26 * s34
      t122 = s34 + s26
      t123 = s25 * s46
      t124 = t101 * s12
      t125 = (-t102 - t98) * t103
      t126 = t3 * epinv
      t127 = 1.0_dp - t126
      t128 = -t127
      t129 = t18 * s25
      t130 = t128 * s13
      t131 = s34 - s46 - s56
      t132 = t24 * epinv
      t133 = t132 * s26
      t134 = s16 + s25
      t135 = t3 * t134
      t136 = s25 + s26
      t137 = s34 - s56
      t138 = t2 * t136
      t139 = s34 - s26 - s56
      t140 = s25 + s56
      t141 = t127 * s12
      t142 = intHs16s25s26s34s56x1311D4eps1()
      t143 = s13 - s16
      t144 = t3 * t143
      t145 = t144 - t66
      t146 = t3 * s25
      t147 = t146 + s12
      t148 = s25 - s13 - s34 + s56
      t149 = s16 + s46
      t150 = t2 * s26
      t151 = t3 * t148
      t152 = s26 + s56
      t153 = s25 - s13 + s16 - s34 + s46
      t154 = t3 * s56
      t155 = (-t153 - t154 - t150 - s12) * s12
      t156 = (-t149 - t151 - t150) * t152 + t155
      t157 = s13 - s16 - s12
      t158 = t157 * s12
      t159 = intHs16s25s26s34s56x1131D4eps1()
      t160 = -s13 + s25
      t161 = -t160 - t79 - t150
      t162 = s25 - s13 + s46 + t67 + t150
      t163 = s25 - s34
      t164 = s16 + s26 + s56
      t165 = t2 * t163
      t166 = t3 * t163
      t167 = -t2 * (s13 - s16 - s46 - s56) + t113 + t166
      t168 = s13 * (t164 + t165)
      t169 = s16 * s56
      t170 = t163 * s46
      t171 = s16 * s26
      t172 = s26 * s56
      t173 = -s56 * t163
      t174 = s12 ** 2.0_dp
      t175 = t2 * (t174 - t90 - t85)
      t176 = t24 * t173
      t177 = -t60 * s26 * t163
      t178 = t3 * (-s16 * t163 - t172)
      t179 = t45 * s25
      t180 = t179 * s34
      t181 = t59 * s25
      t182 = (-s25 + s13 + s12 - t72) * s26
      t183 = (-s12 + s16 + s26 - s34) * s46
      t184 = s13 * s16
      t185 = (t161 + s34) * s34
      t186 = (t162 + s56) * s56
      t187 = s12 * s13
      t188 = t24 * t89
      t189 = intHL0s25s260s56x1010D0eps0()
      t190 = s12 + s16
      t191 = t3 * t160
      t192 = s25 - s13 + s16 - s34 + s12
      t193 = s12 + s13 + s16 + s26 - s45
      t194 = s13 ** 2.0_dp
      t195 = (t192 + t98) * s45
      t196 = -(-t190 - t54 + s25) * s25 - (-t190 * t2 + s25 - s26 - t54)
     # * s26 - (s16 - s26 + s12 - t191 + s34) * s34 + s46 * t193 - (t3 *
     # (s25 - s13 - s34) - s16 + s26 + s45 - s12 + s56) * s56 + t174 - t
     #194 + t94 - t195 + t70 * s16
      t197 = intHL0s25s260s56x1020D2eps0()
      t198 = intHs16s25s26s34s56x1211D4eps0()
      t199 = s12 + s16 + s26 + s46
      t200 = intHs16s25s26s34s56x1121D4eps1()
      t201 = t2 * t97
      t202 = intHs16s25s26s34s56x1211D4eps1()
      t203 = intHs16s25s26s34s56x1112D4eps0()
      t204 = s16 + s46 + s12
      t205 = t204 * t3 + t148 + t150
      t206 = intHs16s25s26s34s56x1111D2eps0()
      t151 = t204 + t151 + t150
      t207 = s13 - s26
      t208 = t3 * t59
      t209 = t3 * (s25 + s16 - s12)
      t210 = (t54 + s12 - t72) * s25
      t211 = t2 * t134
      t212 = intHs16s25s26s34s56x1113D4eps0()
      t213 = t24 * s26
      t214 = t2 * t143
      t215 = t88 - t90
      t216 = t24 * t85
      t217 = t3 * t75
      t218 = intHs16s25s26s34s56x1221D4eps0()
      t219 = 9.0_dp * s12
      t220 = t60 * s12
      t221 = t90 + t85
      t222 = t3 * t89
      t78 = -t2 * t221 + t24 * t88 + (-t79 + t78) * s25 + (t219 - s13 - 
     #t4) * s26 + (-t2 * (s13 - s26 + s12) + t79 + t179) * s34 + (-s25 +
     # s34 + s26 + t57) * s46 + (-t68 + t220 - t166 + t150 + s56) * s56 
     #+ t83 + t222
      t83 = intHs16s25s26s34s56x1122D4eps0()
      t223 = t2 * s13
      t224 = t76 * s13
      t225 = t45 * s16
      t226 = 13.0_dp
      t227 = t73 * s25
      t228 = 8.0_dp
      t229 = t228 * s26
      t230 = t226 * s25
      t231 = t60 * s16
      t232 = t45 * s13
      t233 = s16 - s26
      t234 = t3 * s34
      t235 = t45 * s12
      t236 = t60 * s34
      t237 = t76 * t85
      t238 = t215 * t45 + (t220 + t223 - s16) * s16 + (-t225 + t66 + t22
     #4) * s25 + (s12 * t226 + s26 - t227 - t72) * s26 + (-t232 + t231 +
     # t230 + t229 - t57) * s34 + (-t233 + t235 + t234 - t4) * s46 + (s1
     #2 * t228 + s46 * t3 + s56 + t150 - t179 + t236 - t56) * s56 - t237
      t239 = intHs16s25s26s34s56x1211D2eps0()
      t240 = t61 - t98
      t241 = (t143 - s25) * s25
      t242 = t3 * t87
      t243 = t2 * t89
      t244 = (-t54 + s16 + t66) * s26 + t240 * s34 - (s25 - s26 - s12) *
     # s46 + (-t3 * (s13 + s34 - s12) + s16 + s25 + s46 + t63) * s56 - t
     #158 + t241 + t243 + t242
      t245 = intHs16s25s26s34s56x1121D2eps0()
      t246 = t70 + s13
      t247 = -s16 + s25
      t248 = s26 + s12
      t249 = t2 * t248
      t250 = t3 * t248
      t251 = s13 - s16 + s34 - s46
      t252 = t3 * (t88 + t87) + t188 + (t246 - s25) * s25 + (t235 - t54 
     #+ s25 + t72) * s26 - (-t247 + t249) * s34 + (t65 + t250) * s46 + (
     #-t251 * t3 + s25 + t113 + t57) * s56 + t94 + t66 * s16
      t253 = intHs16s25s26s34s56x1212D4eps0()
      t254 = intHs16s25s26s34s56x1131D4eps0()
      t255 = intHs16s25s26s34s56x1112D4eps1()
      t256 = t24 * t97
      t257 = intHs16s25s26s34s56x1311D4eps0()
      t258 = intHs16s25s26s34s56x1221D4eps1()
      t259 = t3 * t68
      t260 = t24 * t65
      t261 = intHs16s25s26s34s56x1122D4eps1()
      t262 = t79 - t57
      t263 = -t160 * t2 + s12 - t229 - t231
      t264 = s16 + s26
      t265 = t3 * t264
      t266 = t265 - t116 - t57
      t267 = t3 * (s25 - s13 + s46 - s12) + t113 - t236 + t62
      t268 = s16 * (-s16 + t75)
      t269 = t3 * ((-t2 * t59 - t160) * s26 - t87 + t268)
      t270 = t24 * (t88 - t89)
      t271 = intHs16s25s26s34s56x1211D2eps1()
      t272 = s12 + s25
      t273 = t131 * t272
      t274 = (-t143 + t70 + s25) * s25
      t275 = t3 * t272
      t276 = t275 * s26
      t277 = intHs16s25s26s34s56x1121D2eps1()
      t278 = s12 + s34
      t279 = t3 * t278
      t280 = intHLs16s25s26s34s56x1212D4eps0()
      t281 = s16 + s26 - s34 - s45 - s46 + s56
      t282 = intHLs16s25s26s34s56x1113D4eps1()
      t283 = s25 + s45 + s46
      t284 = epinv * t283
      t285 = intHLs16s25s26s34s56x1212D4eps1()
      t286 = intHs16s25s26s34s56x1132D6eps0()
      t287 = t2 * t278
      t288 = t3 * t164
      t289 = -t288 + t287
      t290 = intHs16s25s26s34s56x1222D6eps0()
      t291 = t24 * t278
      t292 = t3 * (s25 - s26 - s56) - t72 + t291
      t293 = intHLs16s25s26s34s56x1114D6eps0()
      t294 = intHLs16s25s26s34s56x1123D6eps0()
      t295 = s26 - s56
      t296 = intHLs16s25s26s34s56x1211D2eps0()
      t297 = s12 + s34 - s56
      t298 = intHLs16s25s26s34s56x1312D6eps0()
      t299 = intHs16s25s26s34s56x1141D6eps0()
      t300 = -s12 + s16 + s26 - s34 + s56
      t301 = intHs16s25s26s34s56x1321D6eps0()
      t302 = -t300 + t146
      t303 = intHLs16s25s26s34s56x1222D6eps0()
      t304 = s25 + s26 + s34 - s56
      t305 = t304 + t70
      t306 = intHLs16s25s26s34s56x1121D2eps0()
      t307 = s12 + s25 + s26
      t308 = intHLs16s25s26s34s56x1132D6eps0()
      t309 = intHs16s25s26s34s56x1312D6eps0()
      t310 = s16 - s34 - s12
      t311 = -t310 + t146
      t312 = -t264 + t287
      t313 = intHLs16s25s26s34s56x1113D4eps0()
      t314 = intHLs16s25s26s34s56x1122D4eps0()
      t315 = intHLs16s25s26s34s56x1122D4eps1()
      t316 = intHLs16s25s26s34s56x1213D6eps0()
      t317 = -s16 + s25 + s34
      t318 = s56 - s12
      t319 = -t3 * t318 + t317
      t320 = intHs16s25s26s34s56x1231D6eps0()
      t321 = t3 * t300
      t322 = s25 - t321
      t323 = -t3 + epinv
      t324 = 1.0_dp - t9
      t325 = intHLs160000x0122D4eps0()
      t326 = intHs160000x0111D0eps0()
      t327 = -s13 + s25 + s26 - s34 + s56
      t328 = intHs160000x0112D2eps0()
      t329 = intHs160000x0121D2eps0()
      t148 = t148 * t2 + t204 + t213
      t204 = intHs160000x0211D2eps0()
      t330 = intHLs160000x0121D2eps1()
      t331 = intHLs16s25s26s34s56x1114D6eps1()
      t332 = intHLs16s25s26s34s56x1123D6eps1()
      t333 = intHLs16s25s26s34s56x1312D6eps1()
      t334 = intHLs16s25s26s34s56x1211D2eps1()
      t335 = intHs16s25s26s34s56x1141D6eps1()
      t336 = intHs16s25s26s34s56x1321D6eps1()
      t337 = intHLs16s25s26s34s56x1222D6eps1()
      t338 = intHs160s26s34s56x1020D2eps0()
      t339 = s12 + s13 + s16 + s26
      t340 = intHLs16s25s26s34s56x1112D2eps1()
      t341 = s46 + s56 + s45 - s12
      t342 = intHs160s26s34s56x1012D2eps0()
      t343 = intHs160s26s34s56x1013D4eps0()
      t344 = intHs160s26s34s56x1021D2eps0()
      t345 = intHs16s25s26s34s56x1130D4eps0()
      t346 = intHs16s25s26s34s56x1120D2eps1()
      t347 = s25 - s13 - s34 + s56 - s12
      t348 = intHs16s25s26s34s56x1120D2eps0()
      t349 = intHs16s25s26s34s56x1210D2eps0()
      t350 = intHs16s25s26s34s56x1310D4eps0()
      t351 = intHs160s26s34s56x1031D4eps0()
      t352 = intHs160s26s34s56x1022D4eps0()
      t353 = -t49 + t70
      t354 = intHs16s25s26s34s56x1220D4eps0()
      t355 = -s25 + s16 + s34 + s26 - s56 + t217
      t356 = intHs160s26s34s56x1021D2eps1()
      t357 = t18 * s12
      t358 = s25 + s26 + s34 + s56
      t359 = intHs16s25s26s34s56x1121D4eps0()
      t360 = intHLs16s25s26s34s56x1132D6eps1()
      t361 = intHLs16s25s26s34s56x1121D2eps1()
      t362 = intHs16s25s26s34s56x1312D6eps1()
      t363 = intHLs16s25s26s34s56x1213D6eps1()
      t364 = intHs16s25s26s34s56x1231D6eps1()
      t365 = intHs16s25s26s34s56x1132D6eps1()
      t366 = intHs16s25s26s34s56x1222D6eps1()
      t367 = intHLs160000x0122D4eps1()
      t368 = intHs16s25s26s34s56x1210D2eps1()
      t369 = intHs16s25s26s34s56x1310D4eps1()
      t370 = intHs160s26s34s56x1031D4eps1()
      t371 = intHs160s26s34s56x1022D4eps1()
      t372 = intHs16s25s26s34s56x1220D4eps1()
      t373 = intHs160000x0111D0eps1()
      t374 = intHs160000x0112D2eps1()
      t375 = intHL0s25s26s34s56x1110D2eps0()
      t376 = intHs160s26s34s56x1011D2eps0()
      t377 = intHs160000x0121D2eps1()
      t378 = intHs160000x0211D2eps1()
      t379 = intHs160s26s34s56x1020D2eps1()
      t380 = intHs160s26s34s56x1012D2eps1()
      t381 = intHs160s26s34s56x1013D4eps1()
      t382 = intHs16s25s26s34s56x1130D4eps1()
      t383 = intHLs160000x0111D2eps0()
      t384 = intHs16s25s26s34s56x1110D2eps1()
      t385 = intHLs160000x0112D2eps0()
      t386 = t225 * propW16
      t387 = t386 * t324
      t388 = t7 * (-t3 - t387 + t5) + t10
      t389 = intHLs160000x0112D2eps1()
      t390 = t387 * t7
      t391 = intHLs16s25s26s34s56x1112D2eps0()
      t392 = intHLs160000x0121D2eps0()
      t393 = intHLs16s25s26s34s56x1112D4eps1()
      t394 = -t2 + t126
      t395 = intHL0s25s260s56x1010D0eps1()
      t396 = intHL0s25s260s56x1020D2eps1()
      t397 = -t3 * t84
      t398 = t60 * s25
      t399 = (t54 - s16 - t398) * t247
      t400 = t2 * t247
      t401 = s12 * (-t161 + s12)
      t402 = t3 * t136
      t403 = t402 * t247
      t404 = t247 * t136
      t405 = s12 * t69
      t406 = (-t70 + s13) * s16
      t407 = -s12 * t84
      t408 = s12 * t94
      t409 = (-t59 - t54) * t85
      t410 = s12 * t90
      t411 = t3 * t86
      t412 = t60 * s13
      t413 = t76 * s16
      t414 = 9.0_dp * s25
      t226 = t226 * s16
      t415 = t2 * t64
      t416 = t2 * t94
      t417 = t228 * t85
      t418 = t228 * s16
      t419 = t223 * s16
      t420 = s13 + s25 - s46
      t421 = t27 * t99
      t422 = t421 * propW16 * s16 * t49
      t423 = -propZ25 * (t3 * (t422 + t86) + (s25 * (t54 - s16) - t407) 
     #* s16 + ((t70 + t400) * s26 + t70 * t134 - t399 + t88) * s26 - (t4
     #01 + t403) * s34 + (t404 + t405) * s46 + ((t420 + t234 - t150 - t6
     #2 - s56) * s56 + t3 * ((t84 - t79) * s25 - t89 - t94) + t216 + (t7
     #0 + s13 - t418 + t115) * s26 + (-s25 - s13 - s12 + t150 + t62 - s3
     #4) * s34 + (s25 + s34 - s26 + s12 - t79) * s46 + t88 + t419) * s56
     # + t408 + t409 + t410) * t2 + t417
      t424 = t421 * t386 * t49
      t425 = t10 * (-(s25 * (s16 - t397) - t407) * s16 - (-(-t70 + t400)
     # * s26 - (s13 - t62 - s12) * s12 + t399) * s26 + (t401 - t403) * s
     #34 + (t404 - t405) * s46 - ((t68 + t234 - t138 - s56) * s56 + t3 *
     # (s25 * t56 - t89) - t216 + (t190 * t3 + s13 - t179) * s26 - (t75 
     #- t138 + s34) * s34 - (-s12 + s25 + s26 - s34) * s46 + t88 - t406)
     # * s56 - t408 + t409 - t410 + t411)
      t426 = intHL0s25s26s34s56x1110D2eps1()
      t427 = intHs160s26s34s56x1011D2eps1()
      t428 = intHs16s25s26s34s56x1411D6eps0()
      t429 = 0.1e1_dp / s25
      t430 = 0.1e1_dp / t27
      t431 = 0.1e1_dp / t49
      t432 = 0.1e1_dp / t99
      t433 = t342 + t344
      t434 = 0.1e1_dp / ecossin ** 2.0_dp
      t435 = 0.1e1_dp / gw ** 2.0_dp
      t436 = t340 * epinv
      t437 = t393 * t394
      t438 = t352 * t353
      t439 = t356 * (-s16 + s34 - s26 - s56 + t357)
      t440 = t433 * s12
      t441 = t338 * t431
      t442 = t290 * t292
      t443 = t8 * t435
      t444 = t443 * t434
      t445 = t444 * t429 * propW34
      t446 = 0.1e1_dp / s16
      t447 = t348 + t349
      t448 = t1 + t14
      t449 = t430 ** 2.0_dp
      t450 = (t329 + t204) * t148
      t451 = t326 * t327
      t452 = epinv2 * t446
      t399 = t452 * (t7 * (t423 * s25 + t24 * (t89 * (s12 + t400) + t409
     #) + t3 * ((-s25 * (t79 + s12 - t71) - t407) * s16 + (s12 * (-t143 
     #+ t4 + s12) - t399 * t3) * s26 - (t404 * t24 + t401) * s34 + (t405
     # + t403) * s46 + ((t68 - t225 + t234 + t415 - s56) * s56 - (s12 - 
     #t412) * s16 + (t66 - t413 - t71) * s25 + (t70 + s13 - t226 + t414)
     # * s26 - (t75 - t225 + t415 + s34) * s34 + (s34 - s26 + s12 + t146
     # - t72) * s46 + t88 - t222 + t417 - t416) * s56 + t408 + t410) + t
     #424) + t425) * t430 * t431
      t401 = (-t395 + t396) * t18
      t403 = t401 * t449
      t404 = t337 * t18
      t405 = s34 * t430
      t407 = t53 * t431
      t238 = t445 * (t12 * (epinv * (t297 * t334 - t50) + t404 * t305 + 
     #t405 * t448) + t432 * (t26 * (t255 * (epinv * t205 - t229 - t256) 
     #- t271 * (-epinv * t244 - t158 - t273 + t274 + t276) + t450 - t451
     #) + t53 * (-t238 * t83 + t239 * t244 + t245 * t252 - t253 * t92 - 
     #t258 * (-t3 * ((t58 + t70 - t72) * s26 + t85 - t87 - t90 - t94) + 
     #t188 - t96 + epinv * t78 + (t54 - t220) * s25 + (-t113 + t54 + s12
     # - t62) * s34 + (-t3 * (s25 - s16 + s34 - s26) - t66) * s46 - (-t2
     #60 - t113 + s12 + t259) * s56 - t93) - t261 * (epinv * t238 + s25 
     #* t262 + s34 * t263 + s46 * t266 + s56 * t267 - t269 - t270 + t91)
     # + t277 * (epinv * t252 - t174 + t181 - t182 + t183 - t184 + t185 
     #+ t186 + t187 + t222 + t94)) - t399 + t403 * t196 * t12) + t407 * 
     #t447 * t347)
      t244 = t189 - t197
      t252 = t377 + t378
      t399 = t323 * t427
      t408 = t323 * t426
      t409 = s12 * t380
      t353 = t353 * t371
      t410 = t379 * t431
      t407 = t26 * (-t292 * t366 + t431 * (-t410 * t339 - t353 + t409)) 
     #+ (-s34 * t22 + t19 * t358) * t430 * t12 - t407 * (t347 * t368 + t
     #355 * t372)
      t415 = t385 * t388
      t423 = t361 * epinv
      t425 = t384 * t323
      t92 = t445 * (t18 * t407 + t432 * (t12 * (t196 * t244 * t449 + t32
     #7 * (t408 + t375) * t430) + t26 * (epinv * (t148 * t252 - t327 * t
     #373) - t35 * (epinv * (-t2 * s25 * t117 + t3 * s25 * t140 + s13 * 
     #(s34 - t135) + (-t137 * t3 + s16 + s46 + t138) * s16 + t139 * s34 
     #+ t123) + (t141 + t3 * (epinv * t131 + s26) - s34 + s46 + s56 + t1
     #07 + t129 + t130 - t133) * s12) + (-t399 - t376) * t431 * t151) - 
     #t52 * t53 * (t3 * (t94 + t90 + t89) + epinv * t92 - (-t56 + t57 + 
     #s25) * s25 - (-t59 * t60 + t58) * s26 + (-t61 + t54 - t63 - t62) *
     # s34 + (-t64 + t67 - t66) * s46 - (-t2 * t69 - s56 + t68 + t70) * 
     #s56 - t93 - t96)) + t415 + t423 * t12 * t307 - t425 * t26 * t431)
      t93 = t200 + t202
      t196 = t389 * (epinv * t388 + t390)
      t78 = t445 * (t12 * (-t46 + t51) + t432 * (t26 * (t151 * t206 + t1
     #98 * t199 + t203 * t205 - t23 * ((t211 - t57) * s26 - t3 * (t88 - 
     #t85) - (t55 - s16) * s16 - (-t207 + t208 + t4 - s34) * s34 + (t134
     # - t70) * s46 - (s34 - t209) * s56 - t210)) + t53 * ((epinv * t199
     # - t113 - t201) * t93 - t218 * t78)) + t196) + t92 + t238 + t445 *
     # (t12 * (-t280 * t281 + t283 * t314 - t285 * (epinv * t281 + s12 +
     # s25 + s26) + t296 * t297 - t303 * t305 + t306 * t307 + t315 * (-s
     #25 - s26 - s12 + t284) + t341 * (-t436 - t391) + t430 * (-t15 * t3
     #58 + t21 * (-s25 - s26 - s56 + t108)) - t437) + t26 * (t431 * (t34
     #6 * (epinv * t347 + s12 + s13 + s16 + s26) - t48 + t438 + t439 - t
     #440 + t441 * t339) + t442 + t199 * t359 * t432) + t354 * t53 * t35
     #5 * t431)
      t92 = t333 * t18
      t199 = t92 - t298
      t205 = t300 * t351
      t238 = s12 * t343
      t281 = t34 * t18
      t305 = t309 * t311
      t355 = t320 * t322
      t388 = t286 * t289
      t407 = t332 * t18
      t453 = s12 * t381
      t454 = t300 * t370
      t455 = t289 * t365
      t311 = t311 * t362
      t322 = t322 * t364
      t456 = epinv * t374
      t215 = t18 * (-t312 * t33 + t431 * (-t339 * t382 - t453 + t454) - 
     #t455 - t311 - t322) + t432 * (t28 * (epinv * ((-s34 * t76 + t113) 
     #* s25 + t2 * (s25 * s56 + t90) - t24 * (t121 - t85) - t3 * (s56 * 
     #t122 - t123) + s13 * (-t79 + s56 + t116 - t115) + (-t117 * t2 + s1
     #6 + t115 + t118) * s16 - t119 - t120) + (t3 * (s25 + s46 + t108 + 
     #t107) - s34 + t105 - t106 + t111 + t112 - t114 + t124) * s12 - t12
     #5) + t156 * t257 + t212 * (-t2 * t215 - (t143 * t24 + s12) * s25 +
     # (t45 * t61 + t72) * s26 + (t70 + t214 - t77 - t213) * s34 + (t65 
     #+ t146 - t66) * s46 - (t122 * t3 - s13 - s16 + s46 + s56 - t4 + t5
     #7) * s56 + t94 + t216 - t217 * s16) + t254 * (-s12 * t167 - t168 +
     # t169 + t170 + t171 - t175 - t176 - t177 - t178 - t180 + t87 + t89
     #) + t456 * t327)
      t131 = t445 * (t12 * (t18 * (t405 * t20 + t307 * t360 + t319 * t36
     #3) - t13 - t47) + t215 * t26 + t53 * (t432 * (-t142 * (-epinv * t1
     #56 + s25 * t145 + t131 * t147 - t3 * (s26 * t147 + t85) + t158) - 
     #t159 * (epinv * (s12 * t167 + t168 - t169 - t170 - t171 + t175 + t
     #176 + t177 + t178 + t180 - t87 - t89) + t3 * (t187 + t186 + t185 -
     # t184 + t183 + t181 - t182 + t94 - t174) + t188)) + t431 * t347 * 
     #(t18 * t369 - t350))) + t445 * (t12 * (t17 * t18 + t282 * (-s25 + 
     #s56 - s12 + t284) + t283 * t313 + (t295 + t275) * (t407 - t294) + 
     #t297 * t199 - t307 * t308 - t316 * t319 - t405 * t16) + t26 * (t30
     # * t312 + (s25 - s16 + s56 + t279) * (-t281 + t31) + t302 * (-t18 
     #* t336 + t301) + t431 * (t339 * t345 - t205 + t238) + t305 + t355 
     #+ t388 + t327 * t328 * t432))
      t156 = t331 * t18
      t167 = -t156 + t293
      t168 = t335 * t18
      t170 = t42 - t29
      t175 = (-t168 + t299) * t300
      t176 = t175 * t429
      t177 = intHs16s25s26s34s56x1411D6eps1() * t18
      t178 = t170 * (s12 + s34 + s56)
      t180 = t429 * t167
      t215 = t18 * t367
      t217 = t330 * epinv
      t284 = intHLs160000x0111D2eps1() * t323
      t297 = t434 * propW34 * t8
      t319 = -t18 * intHLs160000x0113D4eps1() + intHLs160000x0113D4eps0(
     #)
      t339 = 0.4e1_dp / 0.3e1_dp
      t347 = 0.2e1_dp / 0.3e1_dp
      t405 = 0.2e1_dp / 0.9e1_dp
      t457 = 0.1e1_dp / 0.9e1_dp
      t458 = t339 * t297 * propW16 * s16 * t324 * t429 * t319
      t459 = 1.0_dp + epinv
      t460 = t382 * t431
      t461 = t345 * t431
      t462 = t282 * epinv
      t463 = t366 * t18
      t464 = t26 * (t463 - t290)
      t465 = t18 * t372
      t466 = (t285 - t315 + t340) * epinv
      t467 = t277 * t459
      t468 = (t465 - t346 - t354) * t431
      t469 = t168 - t299
      t470 = t26 * t469
      t471 = t45 * t470
      t472 = t363 * t18
      t473 = t18 * t364
      t474 = t365 * t18
      t475 = -1.0_dp + epinv
      t476 = t475 * t336
      t477 = t475 * t362
      t478 = t475 * t33
      t479 = t475 * t34
      t480 = t475 * t365
      t481 = 1.0_dp - epinv
      t482 = t475 * t366
      t483 = t335 * t475 + t299
      t11 = t11 * t391
      t484 = t482 + t290
      t485 = t9 * (t484 * t8 + t7 * ((t285 - t315 + t340) * epinv + t277
     # * (1.0_dp + epinv) + t337 * t481 + t245 + t280 - t303 - t314))
      t486 = t475 * t363
      t335 = t9 * (t7 * (t333 * t475 + t258 + t261 + t298 + t30 + t301 +
     # t309 + t31 + t316 + t462 + t476 + t477 + t478 + t479 + t486 + t52
     #) + t8 * (t335 * t481 - t299))
      t487 = t9 * t7
      t488 = 12.0_dp
      t335 = t2 * t485 + t228 * t7 * (t479 + t30 + t31 + t52 + t258 + t2
     #61 + t301 + t309 + t476 + t477 + t478) + t24 * (t7 * (t316 + t298 
     #- t245 - t467 + t462 - t92 - t472) + t10 * (t320 + t286 + t159 - t
     #474 - t473)) - t3 * (t7 * (t280 - t303 + t337) + t12 * t313 + t7 *
     # ((-t337 + t340 - t315 + t285) * epinv - t314) + t10 * (t362 - t36
     #3 - t333 + t336 - t301 - t309 + t316 + t298 - t258 - t261 - t52 - 
     #t30 - t31 + t33 + t34) + t10 * epinv * (-t362 + t363 + t333 - t336
     # + t282 - t33 - t34)) - t45 * t335 + 16.0_dp * t7 * (t364 * t475 +
     # t159 + t286 + t320 + t480) + t488 * t7 * (t9 * (t364 * t481 + t36
     #5 * t481 - t159 - t286 - t320) + t290 + t482) + 24.0_dp * t7 * t48
     #3 - 18.0_dp * t487 * t483 + t10 * (t391 - t303 - t314 + t280 - t24
     #5 + t466 - t467 + t404) + t11 - 9.0_dp * t487 * t484
      t476 = (t18 * (t363 + t331) - t293 - t316) * t12
      t483 = (t18 * (t332 + t333) - t294 - t298) * t12
      t484 = t360 * t18
      t485 = t369 * t431
      t489 = (-t308 + t280 - t282 + t484) * t12
      t490 = t452 * t430
      t491 = (t334 + t340 + t361) * epinv
      t492 = t128 * t285
      t493 = -t404 + t303
      t494 = t493 * t12
      t495 = t2 * (t464 + t494) - t24 * (t26 * (-t18 * (t364 + t365) + t
     #159 + t286 + t320) + t483) - t3 * (t26 * (t18 * (-t362 - t336 - t3
     #3 - t34 - t485) + t350 * t431 + t258 + t261 + t30 + t301 + t309 + 
     #t31 + t52) + t489) - t45 * (-t470 + t476) - t12 * (t391 + t306 - t
     #315 + t296 + t492 + t491 + t490) - t26 * (t431 * (-epinv * t346 + 
     #t18 * (-t410 - t372 + t368) - t348 - t349 + t354 + t441) - t245 - 
     #t467)
      t496 = -t24 + t100
      t379 = t18 * t379
      t497 = t18 * t368
      t494 = -t2 * t494
      t73 = t228 * t26 * (t18 * (t33 + t364) - t30 - t320) - t24 * (t26 
     #* (-t18 * (t336 + t362) - t218 - t253 + t301 + t309) + t483) - t3 
     #* (t26 * (-t104 * t28 + t142 * t18 - t159 * t496 + t431 * (-t18 * 
     #(t370 + t371 + t372 + t381 + t382 + t369) + t343 + t345 + t350 + t
     #351 + t352 + t354) + t23 - t257) + t489) - t45 * (t26 * (t18 * (t2
     #61 - t32 - t34) - t212 - t254 + t29 + t31 - t83) + t476) + t464 * 
     #t76 + t470 * t488 + t73 * t26 * (-t286 + t474) - t12 * (t391 + t30
     #6 - t315 + t296 + t490) - t26 * (t497 * t431 + t239 + t245) - t12 
     #* (t492 + t491) - t26 * (t128 * t35 - t18 * t271 + t277 * t323 + t
     #431 * (t18 * (t380 + t346 + t356) + t431 * (-t379 + t338) - t342 -
     # t344 - t348 - t349) + (-t52 - t258) * (-t60 + t132)) - t494
      t470 = intHLs16s25s26s34s56x1231D6eps1()
      t476 = intHLs16s25s26s34s56x1321D6eps1()
      t483 = intHLs16s25s26s34s56x1221D4eps0()
      t489 = intHLs16s25s26s34s56x1231D6eps0()
      t492 = intHLs16s25s26s34s56x1321D6eps0()
      t498 = intHLs16s25s26s34s56x1221D4eps1()
      t499 = t24 * (-t286 - t159 + t474) * t26 + t3 * (t12 * (t18 * (t47
     #0 + t476) + t313 + t462 - t489 - t492) + t26 * (t18 * (t364 + t33 
     #- t460) - t261 - t30 - t320 + t461)) + t471 - t12 * ((-t315 + t340
     #) * epinv + t18 * (-t337 + t498) + t285 + t303 - t314 + t391 - t48
     #3) - t26 * (t459 * (t52 + t258 - t277) + t218 - t245 + t253 + t290
     # - t463 + t468)
      t500 = intHLs160000x0211D2eps0()
      t501 = intHLs16s25s26s34s56x1121D4eps0()
      t502 = intHLs16s25s26s34s56x1211D4eps0()
      t503 = intHL0s25s260s56x1012D2eps1()
      t504 = intHL0s25s260s56x1013D4eps1()
      t505 = intHL0s25s260s56x1031D4eps1()
      t506 = s13 + s25 + s26 + s56
      t507 = s25 + s26 - s34 + s56
      t508 = t97 + t150
      t509 = intHLs16s25s26s34s56x1141D6eps0()
      t510 = t152 + t74
      t511 = t3 * t152
      t74 = t511 + t74
      t512 = intHL0s25s260s56x1011D2eps0()
      t513 = intHLs160000x0211D2eps1()
      t514 = intHLs16s25s26s34s56x1121D4eps1()
      t515 = s13 - s34
      t516 = intHLs16s25s26s34s56x1311D4eps0()
      t517 = t18 * s13
      t518 = s34 * t459 + t517
      t519 = propW16 * s16
      t520 = t519 * t324
      t521 = t7 * (t2 * (t520 + t9) - t3) + t10
      t522 = s13 - s34 - s45 - s46 - s56
      t523 = 11.0_dp
      t524 = 9.0_dp * s26
      t525 = s45 + s12
      t526 = t488 * s13
      t527 = (t70 + s13 + t72) * s25
      t55 = (t55 + s16) * s16
      t528 = s46 ** 2.0_dp
      t529 = t84 * s13
      t230 = t3 * (t55 + t527 + t90) - t45 * (-t529 + t195 + t528) - t60
     # * t87 + t216 + t243 + (s13 * t523 + s12 + t230 + t413) * s26 - (t
     #24 * t56 + t179 + t524 + t70) * s34 - (s26 * t523 + t45 * t525 + t
     #146 + t260 - t526) * s46 - (s45 * t45 + t2 * t65 - t523 * t68 + s2
     #5 + t220 + t229) * s56
      t260 = t228 * s25
      t530 = t24 * t525
      t531 = (t66 + t62) * s25
      t532 = t195 + t528 + t87
      t533 = t24 * t532
      t534 = t3 * t90
      t535 = -t2 * (-t529 - t85) + (t232 + t72 + t260 - s26) * s26 - (t1
     #13 + t70 + s13 + t398 + t72) * s34 - (t530 - t234 + s25 + t72 - t2
     #24 + t229) * s46 - (s46 * t228 + s25 - t224 - t234 + t524 + t530 +
     # t72) * s56 + t55 + t531 - t533 + t534
      t536 = t76 * s26
      t537 = t76 * t207
      t538 = t76 * s46
      t539 = t195 + t528
      t540 = t90 + t89
      t530 = -t2 * (-t529 + t87 - t85) - t24 * t539 + t3 * t540 + (s12 +
     # t412 + t414 + t62) * s26 - (t70 + s13 + t398 + t72 + t536) * s34 
     #- (t530 - t234 + s25 + t72 - t537) * s46 + (-s45 * t24 + s34 + t23
     #2 - t538 - t63 - t66 - t79) * s56 + t531 + t55
      t531 = t3 * t58
      t541 = t3 * t56
      t542 = t2 * t525
      t543 = s13 - s26 - s46
      t544 = t3 * t85
      t545 = t2 * ((s16 + t531) * s26 - t195 - t528 + t529 - t87) - (s12
     # + t541 + t4 + t213 - s34) * s34 - (-t207 * t45 + s25 + t542 + t67
     #) * s46 - (-t45 * t543 + s25 + t542 + t67) * s56 + t527 + t55 + t5
     #44
      t546 = s25 + s16 + s12
      t547 = s25 + s16 - s34 + s45 + s12
      t548 = -s12 - s16 - s25 + s13
      t549 = t548 * s13
      t242 = -(t546 - t223) * s26 - s34 * t207 - (t547 - t54 + t150 + s4
     #6) * s46 - (-t2 * t68 + t3 * (s25 + s16 - s34 + s12) + s45 + t113)
     # * s56 - t195 - t549 - t243 - t242
      t550 = s16 + s25 - s34
      t551 = t24 * t207
      t552 = t3 * t525
      t553 = (t75 + s16) * s16
      t554 = (t75 + t79 + s25) * s25
      t555 = (t211 + t71) * s26 - t3 * (-t529 + t195 + t528 + t87) + (-t
     #75 - t135 - t150 + s34) * s34 - (t552 + t550 - t551) * s46 - (-t24
     # * t543 + t550 + t552) * s56 + t553 + t554
      t556 = t2 * s45
      t55 = -t2 * (-t529 + t195 + t528 - t89) - t3 * (t87 - t85) + (s12 
     #+ t412 + t77 + t62) * s26 + (-s12 - t541 - t63 - t4 + s34) * s34 -
     # (t542 - t232 + t67 + s25 + t63) * s46 - (-t60 * t68 + t250 + t556
     # + t65) * s56 + t527 + t55
      t232 = t79 + s12
      t527 = s12 * s16
      t81 = -t3 * (t195 + t528 - t89) + (t232 + s25) * s25 + (t134 * t24
     # + t75) * s26 + (-s12 - t135 - t213 + s34) * s34 - (-t2 * t207 + t
     #550 + t552) * s46 + (t3 * (s13 - s45) - t248 - t81 - s56) * s56 + 
     #t187 - t194 + t94 + t527
      t248 = s45 + s46 + s56
      t528 = -t24 * t248 + t412
      t542 = t2 * t248
      t552 = intHLs16s25s26s34s56x1211D4eps1()
      t557 = t412 - t542
      t558 = t3 * t248
      t559 = t223 - t558
      t560 = epinv * t522
      t561 = t3 * t547
      t562 = s46 + s26
      t192 = t192 * s13
      t563 = s16 - s34 + s45 + s12
      t564 = t2 * t563
      t565 = s16 - s34 + s12
      t566 = t228 * s13
      t567 = (t190 - t412 + s25) * s25
      t568 = t157 * s13
      t569 = -s25 + t71
      t570 = t3 * t207
      t571 = t3 * t543
      t572 = (t547 - t571 + s56) * s56
      t573 = (t547 - t570 + s46) * s46
      t574 = t71 * s26
      t575 = t190 - t54
      t576 = s25 - s13 + s26
      t577 = epinv * t242 + (t575 + s25) * s25 - s34 * t160 + (t3 * t576
     # + s46 + t563) * s46 + (t3 * (s25 - s13 + s46 + s26) + t563 + s56)
     # * s56 + t195 + t568 + t191 * s26
      t578 = t54 * s26
      t579 = (t547 + t98 - t223 + s46) * s46
      t547 = (t3 * t562 + s56 - t223 + t547) * s56
      t580 = s13 + s34
      t581 = t54 + s34
      t582 = t317 + t54
      t583 = s13 - s45 - s46 - s56
      t584 = s16 + s26 - s34 + s46 + s56
      t127 = t3 * t127 * s13
      t585 = t358 + t223
      t586 = s25 + s34
      t587 = s16 - s26 - s56
      t588 = t3 * t586
      t589 = t588 - t587 + t412
      t590 = t60 * epinv
      t591 = intHLs16s25s26s34s56x1131D4eps0()
      t592 = intHLs16s25s26s34s56x1311D4eps1()
      t593 = intHL0s25s260s56x1011D2eps1()
      t594 = intHL0s25s260s56x1012D2eps0()
      t595 = intHL0s25s260s56x1013D4eps0()
      t596 = intHL0s25s260s56x1021D2eps0()
      t597 = intHL0s25s260s56x1021D2eps1()
      t598 = t58 * t18
      t599 = intHLs16s25s26s34s56x1141D6eps1()
      t600 = intHL0s25s260s56x1031D4eps0()
      t601 = intHL0s25s260s56x1022D4eps0()
      t602 = t152 + t531
      t603 = t84 * t134
      t604 = s16 * s25
      t605 = s16 * t84
      t606 = (t79 + t150) * s16
      t607 = s25 * t136
      t608 = (-t117 + s16) * s16
      t609 = s25 * s34
      t610 = s16 + t84
      t611 = t94 * t610
      t612 = s16 * t90
      t613 = t79 - t4
      t614 = t3 * (t607 - t94)
      t615 = t402 * s25
      t616 = -t146 + s16
      t617 = s16 * s34
      t618 = (-s13 - s34 + s46) * s25
      t619 = (-s25 + s16 - s26) * s25
      t620 = s16 * t65
      t621 = (t157 - s25) * t85
      t422 = t422 * t3
      t603 = t10 * (s25 * t94 + ((s16 + t4) * s26 + t216 + t416 + t603 +
     # t604) * s26 - ((s25 + s16 + s26) * s25 + t605 + t606) * s34 + (t6
     #08 + t607) * s46 + (s56 * t134 + s46 * t134 - t3 * ((-s25 - s16 + 
     #s34 - s26) * s16 - t85) + t603 - t609 + t115 * s26) * s56 - t157 *
     # t85 + t86 + t611 + t612)
      t622 = epinv * s26
      t623 = s13 + s16 + s26 - s34 - s45 - s46
      t624 = intHL0s25s260s56x1022D4eps1()
      t625 = t584 + t71
      t626 = -epinv * t76 + t24
      t627 = intHLs16s25s26s34s56x1131D4eps1()
      t101 = t101 * s13
      t628 = t3 * t27
      t629 = t18 * t371
      t630 = -t629 + t352
      t631 = t376 * t431
      t632 = t18 * t624
      t433 = -t356 * (-s13 - s16 + s34 + s46 - s26 + s45 + t560) - t433 
     #* t522
      t81 = t431 * t433 + t432 * (-t218 * t530 + t23 * t555 - t230 * t83
     # + t239 * t242 + t245 * t81 - t253 * t535 - t258 * (t2 * t532 + t2
     #4 * t568 - t3 * t569 * s26 + epinv * t530 + s34 * t569 + (t564 + t
     #113 - t224 + t115) * s46 + (t45 * t562 + t115 - t224 + t564) * s56
     # + t567) - t261 * (epinv * t230 + t3 * ((t561 - t412 + t213) * s46
     # + (t24 * t562 - t412 + t561) * s56) - t45 * t192 + t533 - t526 * 
     #s26) + t271 * t577 + t277 * (epinv * t81 - t192 * t3 + t195 + t547
     # - t574 + t579) + t35 * (epinv * t555 - t192 + t195 + t572 + t573 
     #- t578) - t52 * (t2 * t539 - t24 * (-t568 - t87) + t222 + epinv * 
     #t535 + (-9.0_dp * s13 + t190 + t4) * s26 + (-t136 + t71) * s34 + (
     #t564 - t537 + t115) * s46 + (t24 * t565 + t398 + t524 + t538 + t55
     #6 - t566) * s56 + t567)) - t463 * t528
      t230 = t597 * (t598 + s26 + s56)
      t433 = t593 * t323
      t526 = (t594 + t596) * t58
      t530 = t598 * t503
      t532 = t296 * t515
      t533 = t334 * t518
      t397 = (t399 * s26 * t26 + t490 * (-t7 * (-t24 * (-t613 * t89 - t6
     #21) + t3 * (((t84 + s25 + s16) * s16 + t90) * s16 + (s16 * t84 + (
     #s16 + t397) * s25 - t417 + t416) * s26 - ((s25 - s13 + s12) * s16 
     #+ t72 * s26 - t614) * s34 + (t608 - t615) * s46 + (s56 * t616 + s1
     #2 * t616 - s16 * t420 + (t231 - t260) * s26 - t3 * (t618 + t617 - 
     #t94) - t216) * s56) + t424 + t5 * ((t605 + t619 + t606) * s34 + (-
     #t608 + t607) * s46 - (-s46 * t247 - t24 * s26 * t247 - t247 * t84 
     #- t3 * (-t620 + t85) + t609) * s56 - t247 * (-s26 * (t84 + t72 + t
     #115) - t87) - t94 * (t84 + s25) - t612 - t95 - t621 + t400 * t89 -
     # t422)) + t603)) * t432 * t431
      t67 = t445 * (t12 * (epinv * (t50 - t513) - t18 * t552 - t280 * t6
     #25 + t285 * (-epinv * t625 + s25 - s26 - s56 - t67 + t71) + t303 *
     # t589 - t314 * (s46 + t224) + t315 * (s13 * t626 - t110 + t115 + t
     #511) - t340 * (-epinv * t223 + s13 + s25 - t110) - t361 * (-epinv 
     #* s16 - epinv * t54 + s13 + s25 + s26 + s56) + t391 * (s46 + t223)
     # + t430 * (-t530 - t512) - t483 * (-s34 + t412) + t496 * (t393 + t
     #514) + t498 * ((t24 - t590) * s13 + s25 + s26 + s56 + t108 + t234)
     # + t532 - t533) + t290 * t26 * t528 + t397)
      t397 = t401 + t189 - t197
      t401 = t18 * t380
      t417 = epinv * t378
      t528 = t332 * t510
      t535 = t58 * t595
      t537 = t506 * t600
      t538 = t33 * t18
      t539 = t538 - t30
      t555 = t18 * t430
      t556 = t555 * (t504 * t58 + t505 * t506)
      t510 = t294 * t510
      t515 = t516 * t515
      t561 = t627 * (t628 + t101)
      t562 = t18 * t381
      t518 = t445 * (t12 * (-t282 * (-t110 + t127 + t146) - t298 * t580 
     #- t308 * t74 + t313 * (s46 + t71) - t316 * t582 - t489 * t585 - t4
     #92 * t581 - t518 * t592 - t556 - t510 + t515 - t561) + t26 * (-t28
     #6 * t557 - t320 * t559 + (-t542 + t71) * t539 + t583 * (-t301 - t3
     #09) + t562 * t522 * t431) + t389 * (-t2 * t520 * t7 + epinv * t521
     #) + t521 * (t217 + t385 + t392))
      t55 = t518 + t445 * (t12 * (t18 * (t333 * t580 + t360 * t74 + t363
     # * t582 + t470 * t585 + t476 * t581 + t528) + t430 * (t537 + t535)
     #) + t26 * (t18 * (t364 * t559 + t365 * t557 + t431 * (t370 * t623 
     #- t382 * t584) + t583 * (t336 + t362)) + t431 * (-t343 * t522 + t3
     #45 * t584 - t351 * t623) + t432 * (t142 * t577 + t159 * (epinv * t
     #55 - t192 * t24 + t3 * (t579 + t547 + t195) - t566 * s26) + t212 *
     # t545 + t242 * t257 + t254 * t55 + t28 * (epinv * t545 + t3 * (-t1
     #92 + t573 + t572 + t195) - t574) + t508 * (t417 + t204))))
      t192 = t156 - t293
      t242 = t599 * t18
      t518 = s13 * t591
      t521 = t297 * t429
      t168 = t521 * (t435 * (t12 * (t192 * t58 + t506 * (t242 - t509) + 
     #t518) + t26 * (t170 * t583 + (-t248 + t54) * (t168 - t299))) + (s1
     #6 * (t215 - t325) - t284 - t383) * t324 * propW16)
      t248 = t445 * t12
      t299 = t248 * (t501 + t47)
      t542 = 0.4e1_dp / 0.9e1_dp
      t545 = 0.1e1_dp / 0.3e1_dp
      t55 = t347 * t168 + t405 * t55 - t457 * (t445 * (t12 * (t46 - t51 
     #- t500 + t502) + t508 * (t12 * t397 * t432 * t449 * t507 - t12 * t
     #408 * t430 * t432) + t401 * t26 * t522 * t431) + t67 + t445 * (t12
     # * (t306 * (t54 + s16) + t430 * (-t601 * t602 - t230 - t433 + t526
     #)) + t26 * t81) + t445 * (t12 * (t430 * (-t375 * t432 * t508 + t60
     #2 * t632) - t404 * t589) + t26 * (t432 * (s26 * (t631 + t203 + t35
     #9) + t508 * (epinv * (t202 + t373) + t198 - t206 + t326) + (t97 + 
     #t213) * (-epinv * (t374 + t377) - t328 - t329) + (s25 - s13 + s16 
     #- s34 + s46 + s56 + s12 + t622 + t98) * (t200 + t255)) + (t465 - t
     #346 - t354) * t431 * t584 + (t3 * (s13 - s34 - s46 - s45) + s16 + 
     #s26 - s56) * t431 * t630))) + t542 * t445 * t26 * t583 * (t281 - t
     #31) - t545 * t299 - t458
      t67 = t33 + t34
      t81 = epinv * t627
      t97 = t67 * t18
      t168 = t17 * t430
      t192 = t192 * t12
      t299 = t339 * t445 * (t170 * t26 + t192) - t347 * t248 * t430 * (-
     #t18 * t19 + t15) - t405 * t445 * (t12 * (t430 * (t18 * (t22 - t503
     # - t597) - t21 * t323 - t1 - t14 + t594 + t596) + t459 * (t315 - t
     #334 + t498) - t296 - t306 + t314 + t483 - t423) + t26 * (-t36 - t2
     #3)) + t542 * t445 * (t12 * (t18 * (t333 + t337 + t360) + t430 * (t
     #18 * (t20 - t504 - t505 - t624) - t16 + t595 + t600 + t601) - t282
     # - t285 - t298 - t303 - t308 - t591 - t592 - t81) + t26 * (t97 - t
     #28 - t30 - t31)) - 0.8e1_dp / 0.9e1_dp * t248 * (t13 * t430 + t18 
     #* (-t168 - t332 - t363) + t294 + t316)
      t339 = t2 + epinv
      t465 = t339 * s13
      t508 = s25 + s16 - s34 + s26 + s56 - t54
      t547 = s13 - s16 - s26 + s34 - s56
      t557 = s13 - s16 + s25 - s26 + s34 - s56
      t559 = t3 * t580 - t164
      t564 = t7 * (-t3 + t387 + t5) + t10
      t567 = s13 - s16 + s25 + s34
      t568 = t3 + epinv
      t569 = t54 - s46
      t259 = s25 + s16 - s34 + s56 + s12 + t150 - t259
      t572 = -t54 + s12
      t574 = s25 + s46 + s12
      t577 = t3 * (s13 - s16 + s34)
      t579 = t3 * t94
      t236 = -t60 * t90 - t216 + (-t225 + s13 - t57) * s25 + (t54 - s12 
     #- t398 - t72 - s26) * s26 + (t264 * t76 + t220 + t414 - t54) * s34
     # + (s25 + s16 + s45 + s12 + t98 + s46) * s46 + (-t3 * (s16 - s46) 
     #+ s45 - t115 + t236) * s56 + t195 + t406 + t529 - t579
      t214 = (t272 - t214) * s26
      t406 = (s13 - s16 - s26 + s34) * s46
      t143 = t143 * s25
      t144 = (-t272 + t144 - t150 + s34) * s34
      t414 = (t572 + s16) * s16
      t581 = (t574 + t150 - t577 + s56) * s56
      t582 = s13 - s46 - s12
      t583 = t87 + t89
      t275 = -t2 * t583 + t534 + (s12 + t541 + s25) * s25 - (s16 - t223 
     #+ t275) * s26 - (t70 + s13 + t211) * s34 + (-t563 - t191 - t150 - 
     #s46) * s46 - (-t2 * t582 - s34 + s45 + t135 + t536) * s56 - t195 +
     # t529 + t553
      t534 = s26 + s34 + s46
      t541 = (t146 + s16) * s26
      t563 = t190 * s16
      t544 = (t70 + t72) * s25 + (-t24 * t264 + s13 - t398 - t66) * s34 
     #+ (-s25 + s13 - s16 - s45 - s12 - t98 - s46) * s46 + (-t3 * t534 -
     # s12 + s13 + s25 - s45 - s56) * s56 - t195 + t541 + t563 + t91 + t
     #544
      t584 = s12 + s16 + s25 + s26 - s34 + s56
      t585 = s12 + s13 + s16 + s25 - s34 - s46 + s56
      t589 = s13 - s56
      t598 = t3 * t589
      t603 = s25 + s46 + s26 + s45 - t598
      t606 = s25 + s34 + s26
      t608 = s46 + s45
      t621 = t3 * t608
      t589 = -t2 * t589 + t606 + t621
      t64 = (-s46 + t157 - s25) * s25 + (-s25 + s16 + s12 + s26) * s26 +
     # s34 * t64 + (t565 + t98 + s56) * s56
      t157 = s13 - s16 + s34 - s46 - s12
      t565 = t69 * s46
      t190 = (t190 + s25) * s25 + (-t84 - t79) * s26 - (t58 - t98) * s34
     # + (t157 - t150 - s56) * s56 + t184 - t565 - t222
      t135 = (t134 + t54) * s26 + (-s13 - s26 - s12 - t135 + s34) * s34 
     #+ (-t525 + t571 - s56) * s56 - t195 + t529 + t553 + t554 - t573
      t553 = t54 + s25
      t554 = s16 + t166
      t571 = t164 + t166
      t573 = s13 * t27
      t534 = s25 * t534 + s26 * t152 + (s16 + s26 + s56 + s12) * s34 + s
     #46 * t152 - t90 - t573
      t623 = s12 * s25
      t625 = -s25 + t72
      t633 = s16 * t134
      t634 = -s34 + s12
      t635 = s16 * t68
      t636 = t79 - s25
      t610 = t610 * s16
      t637 = t625 * s26
      t638 = s16 * t264
      t639 = s25 * t87
      t640 = s26 + s34 + s56
      t641 = s25 * s26
      t642 = s13 * t49
      t232 = t10 * (-s12 * t85 + (t610 + t565) * s16 + (t24 * t94 + t3 *
     # (-s25 * t247 + t527) - t623 - t184 + t637) * s26 + (-s26 * t62 - 
     #t3 * t633 - t605 + t607) * s34 + (t636 * s56 + t2 * t620 + t3 * (t
     #619 + t527) - s25 * t634 - t635 + t231 * s26) * s56 + t612 - t86 +
     # t604 * t232)
      t449 = t449 * (-t189 + t197)
      t619 = t218 * t554
      t643 = t334 * t459
      t644 = t430 * t12
      t645 = t644 * t448
      t560 = t445 * (t12 * (t430 * (t15 * t589 + t21 * (s25 + s46 + s26 
     #+ s45 + t560 - t598) + t597 * (-t106 + t517) + t601 * t553) + t643
     # * s13) + t432 * (t26 * (t135 * t23 + t190 * t245 - t203 * t259 - 
     #t239 * t64 - t328 * t585) + t543 * (t26 * (epinv * t373 - t206 + t
     #326) - t375 * t12 * t430) + t83 * t53 * t236 + t26 * (epinv * t93 
     #+ t198 + t359 - t619) * t584) + t463 * t26 * t557 + t645 * t522)
      t598 = t395 - t396
      t646 = t258 * t554
      t95 = t490 * (t7 * (t228 * s25 * (s56 * (s25 + s26) + t641) - t24 
     #* (-s25 * t89 + (t94 + t85) * s34 - t639 - t86 + t609 * s56 + t609
     # * s26) + t3 * (s12 * (-s16 * (t640 + t146 - s16) + t146 * t27) + 
     #s16 * (-s25 * t49 + s46 * t49 - t121 - t172 + t638 - t642 - t87 + 
     #t90)) + t424 - t5 * ((-t623 + t565 + t90) * s16 + (s25 * t147 + s2
     #6 * t134 - t184 + t579) * s26 - (t3 * t638 + t605 + t607) * s34 + 
     #((-s13 - s34 + s46 + s26 + s16) * s16 + s25 * (t634 + t402)) * s56
     # + t272 * t85 + t94 * t84 + t95 + t639 + t422)) - t232) * t431 + t
     #261 * t53 * (epinv * t236 - t143 + t144 + t214 + t222 - t406 + t41
     #4 - t529 + t581)
      t147 = t392 * t564
      t95 = t445 * (t12 * (epinv * t513 + t18 * (t430 * (t430 * t432 * t
     #534 * t598 - t22 * t522 + t503 * t58) + t552) + t430 * (-t526 + t5
     #12) + t500 - t501 - t502 + t430 * (-t426 * t432 * t543 + t593) * t
     #323) + t432 * t95 + t147 + t432 * (-t253 * t275 - t255 * (epinv * 
     #t259 - t113 - t201) - t52 * (epinv * t275 + t143 - t144 - t214 - t
     #222 + t406 - t414 + t529 - t581) + t584 * (-epinv * t646 + t399 * 
     #t431)) * t26)
      t143 = t248 * s13
      t144 = -t242 + t509
      t42 = -t42 + t29
      t201 = t521 * (t435 * (t12 * (s13 * t144 + t167 * t68) + t42 * t58
     #0 * t26) + (s16 * (t18 * (-t367 + t389) + t325 - t385) + t284 + t3
     #83) * t324 * propW16)
      t214 = t18 * t504
      t232 = t18 * t20
      t236 = (t332 + t363) * t18
      t259 = t444 * propW34
      t252 = t259 * (t26 * (t18 * (-t429 * (t33 * t559 + t34 * t567 + t3
     #65 * t547) - t362) + t309 + t429 * (t212 * t544 + t28 * (epinv * t
     #544 + (t546 - t570 - s34) * s34 + t118 * t580 - t549 + t578) + t58
     #4 * (epinv * (-t142 * t152 + t159 * t571) + t204 + t329)) * t432) 
     #+ t429 * (t430 * (-s13 * t600 - t13 * t603 + t522 * (-t232 + t16) 
     #+ t58 * (t214 - t595)) - t236 * t569) * t12) + t445 * (t12 * (s13 
     #* (t18 * (t430 * t505 - t333 - t476) + t298 + t492) + s34 * t516 +
     # t282 * t420 + t580 * t591 + t592 * (s13 + t108) + t627 * (s13 * t
     #568 + t108) + t569 * (t294 + t316) + t555 * t17 * t603) + t26 * (t
     #286 * t547 + t30 * t559 + t31 * t567) + t26 * (epinv * t252 - t152
     # * t257 + t254 * t571) * t432 * t584)
      t64 = -t347 * t201 - t405 * t252 + t457 * (t95 + t560 + t445 * (t1
     #2 * (s13 * t296 - s46 * t391 + t18 * (-t430 * (t19 * t589 + t553 *
     # t624) + t514) - t280 * t547 - t285 * (-t459 * t49 + t465) - t306 
     #* t508 - t314 * (-t49 + t54) - t315 * (-epinv * t49 + t459 * t54 +
     # s25) + t340 * (-t110 + s25) - t361 * (epinv * t508 - s13) - t483 
     #* (t234 + s13) - t498 * (s34 * t126 + t465) + t449 * t432 * t534) 
     #+ t26 * (-t290 * t557 + t432 * (epinv * (t135 * t35 + t190 * t277 
     #- t271 * t64 - t374 * t585) + t631 * t584)) + t217 * t564)) - t542
     # * t143 * (-t18 * (t360 + t470) + t308 + t489) - t545 * t143 * t49
     #3 - t458
      t95 = -t164 + t279
      t135 = t247 + t279
      t190 = -t164 + t287
      t201 = s25 + s26 + s45 + s46 + s56
      t252 = -s34 + s46 + s45
      t275 = s16 * t201
      t283 = t283 * s34
      t287 = t252 * s26
      t252 = t252 * s56
      t367 = t98 * s56
      t406 = t172 + t85
      t414 = s25 * t152
      t444 = s56 * t136
      t465 = t10 * propW34
      t493 = epinv2 * t8
      t508 = -s12 + s25 - s26 - s56
      t522 = epinv * t508
      t526 = -s25 + s16 + s26 + s56 + t70
      t529 = -s16 + s34 + s56 - t146 + t66
      t534 = t640 - t4 + t57
      t543 = -s16 + s25 - s26 + s34 - s56
      t544 = t3 * t543
      t546 = epinv * t99 + s13 - s25 - s26 - s45 - s46 - s56
      t547 = t122 + t154 - t4 + t57
      t549 = t3 * t317
      t553 = t136 + t154
      t557 = t2 * s56
      t559 = t606 + t557
      t560 = t18 * s34
      t565 = -s13 + s25 + s26 + s45 + s46 + s56
      t567 = s12 + s16 + s26 - s34 - s45
      t569 = -t146 + s34 + t66
      t570 = s16 - s56 - t235 - t234 + t115
      t578 = s25 + s16 - s56 - t70
      t579 = t369 + t372 + t382
      t381 = t278 * t381
      t580 = t30 * t190
      t581 = t31 * t135
      t584 = t302 * t336
      t585 = t33 * t190
      t589 = t34 * t135
      t603 = t454 * t431
      t302 = t301 * t302
      t631 = t278 * t343
      t634 = t565 * t18
      t638 = t431 ** 2.0_dp
      t639 = (t346 + t368) * t18
      t647 = t356 * (-s16 - s26 - s56 + t560 + t357)
      t434 = t434 * t429 * t435
      t275 = t434 * ((t12 * (epinv * (t340 * t567 - t565 * (t334 + t361)
     #) + t430 * (-t103 * t448 + t15 * t559 + t21 * (-epinv * t103 + s25
     # + s26 + t154)) + t449 * t193) + t26 * (t292 * (t463 - t290) + t43
     #1 * (t278 * (-t401 + t342 + t344) + t307 * (-t639 + t348 + t349) +
     # t95 * (t629 - t352) - t647) - t206 - t338 * t307 * t638)) * propW
     #34 * t8 - t493 * (t7 * propW34 * (t3 * (-t2 * t414 - t3 * t406 - t
     #252 - t275 + t283 - t287 + t642 - t87 - t89) + t5 * (t3 * (t444 + 
     #t641) + t252 + t275 - t283 + t287 - t642 + t85 + t87 + t89)) - t46
     #5 * (-t252 + t283 - t287 + t642 - t275 - t87 - t89 + t85 - t367)) 
     #* t446 * t430 * t431)
      t193 = t275 + t445 * (t12 * (t391 * t567 + t546 * (-t285 - t315) +
     # t565 * (-t296 - t306) + t99 * (-t280 - t314) + t555 * (t193 * t43
     #0 * t598 + t103 * t22 - t19 * t559)) + t26 * (t218 * t534 + t23 * 
     #t578 + t239 * t508 - t245 * t526 + t253 * t547 - t258 * (-epinv * 
     #t534 + t544 + t66) - t261 * (epinv * t570 - t288 + t291) + t271 * 
     #(s25 + s12 + t522) - t277 * (epinv * t526 - s12 + s16 + s26 - s34 
     #+ s56) + t35 * (epinv * t578 + s12 + s34) - t52 * (-epinv * t547 +
     # t549 + t66) - t83 * t570 + t379 * t307 * t638))
      t170 = t26 * (t170 * t278 * t429 + t176 + t177 - t428)
      t20 = -t347 * t259 * (t12 * t180 * t565 + t170) - t405 * (t445 * (
     #t12 * (t430 * (t103 * t16 + t13 * t553) + t634 * (-t333 - t337 - t
     #360)) + t26 * (t18 * (-t603 + t455 + t584 + t311 + t322 + t585 + t
     #589) - t198 - t203 - t305 - t355 - t359 - t388 - t302 - t631 * t43
     #1 - (t345 + t350 + t354) * t431 * t307)) + t445 * (t12 * (-t282 * 
     #t546 - t313 * t99 + t565 * (t298 + t303 + t308) - t555 * (t103 * t
     #20 + t17 * t553)) + t26 * (-t142 * (s12 + t522 + t146) - t159 * (-
     #epinv * t569 - t321) + t212 * t529 + t254 * t569 - t257 * t508 - t
     #28 * (-epinv * t529 + t279) + t323 * (-t200 - t202 - t255) + t431 
     #* (t18 * (t307 * t579 + t381) + t205) - t580 - t581))) - t457 * t1
     #93 - t542 * t248 * t565 * (-t236 + t294 + t316)
      t103 = -t2 * t27 * s34 + t3 * (t444 + t90) - s13 * (s25 + s16 + s2
     #6 + s56 - t234) - s16 * t139 - (s46 - s26 + s12) * s25 + t85 + t87
     # + t89
      t193 = s46 + s56 + s12
      t275 = s13 * t587 + t507 * (t193 + t98)
      t448 = s16 + s25 - s34 + s56
      t508 = s25 - s34 + s26
      t121 = t3 * (s56 * t508 - t609) - s13 * t448 + s16 * (t508 + t154)
     # + s56 * t193 - t121 + t641 + t85 + t90
      t508 = t24 * s34
      t522 = t136 * s46
      t59 = -t2 * t85 - t24 * t90 + (-t59 + t223) * s25 + (-t59 - t146 +
     # s26) * s26 + (t79 + t77 - t551) * s34 + (-s16 + s46 + s26 + s12 -
     # t4 + t508) * s56 + t522 + t54 * s16
      t526 = (s46 + s56 - t548 + t98 - s34) * s34
      t529 = t3 * t640
      t534 = s46 + s12
      t546 = t3 * ((s46 + s26 + s12) * s25 + s34 * t507) + s13 * (-t234 
     #+ s16) + (s25 + s34 - s26 - s56) * s16
      t444 = t3 * t444 - (-s16 + s34 - s46 - s12) * s25 - s26 * t117 - s
     #34 * s56 - t573 + t85 + t87 + t4 * s26
      t547 = s12 + s16 + s26 - s34 - s45 - s46 + s56
      t548 = -s12 + s25 + s45 + s46
      t551 = epinv * t548
      t553 = -s13 * t636 + (-t146 + t231) * s25 + t563
      t559 = s16 * s46
      t563 = t604 * t575
      t569 = (t79 + s13) * t85
      t570 = (s34 * t636 - t559) * t136
      t573 = 9.0_dp * s16
      t575 = s25 + s13 - s12
      t578 = (-t625 * s26 - t553) * s26
      t625 = t340 * t341
      t629 = t432 * t53
      t640 = t152 * t239
      t648 = t391 * t341
      t649 = t555 * t12 * t429
      t22 = t259 * (t26 * (t431 * (t429 * (t18 * (s34 * (-t371 + t380) -
     # t372 * t543) - t338 + t48 - t108 * t356) + t497) + t429 * (epinv 
     #* (t255 * t327 + t258 * t546 - t275 * t277) + t93 * (epinv * t151 
     #- s12 + s13 - s16 - s25 + s34 - s46 - s56 - t98) - t640 * t327) * 
     #t432) - t649 * (t21 + t22) * s56 + (-t136 * t306 - t137 * t296 + t
     #280 * t547 + t285 * (epinv * t547 + s25 + s26) + t303 * t304 - t31
     #4 * t548 - t315 * (-s25 - s26 + t551) + t437 + t648) * t429 * t12)
      t93 = t152 * t271
      t108 = t429 * t432
      t356 = t430 * t507
      t497 = t50 * epinv
      t507 = t432 * t430
      t84 = t490 * (t7 * (t24 * (t578 + t570 - t569 - t563 + t86) - t3 *
     # ((-s25 - s13 + s46 + s12 + t225 - t234 + t150 + s56) * s56 + (-t4
     #12 + t66) * s16 + (t54 + t573) * s25 + (t84 - t4 + t226) * s26 + (
     #t575 - t225 - t150 + s34) * s34 + (-t117 + t72) * s46 + t222 - t21
     #6 + t416) * s56 + t424 - t5 * (((t582 + t234 - t150 - t62 - s56) *
     # s56 - t3 * (t94 + t89 - t85) + (-t70 + t223) * s16 - (s13 + t231)
     # * s25 + (t575 - t418) * s26 + (t84 + t150 + t62 - s34) * s34 + (t
     #117 - t79) * s46) * s56 - t563 - t569 + t570 + t86 + t578 + t422))
     # + t10 * ((-t553 - t637) * s26 + ((-t582 + t166 + t150 + s56) * s5
     #6 + t3 * (t89 + t85) - (s13 + t72) * s25 - (-t84 + t79 - t4) * s26
     # - (t84 + t146 + t150 - s34) * s34 - s46 * t117 + t184) * s56 + t8
     #6 - t563 - t569 + t570)) * t431 * t432 - t415 - t196
      t84 = t259 * (t26 * (t431 * (t429 * ((t399 + t376) * t432 * t151 +
     # t425) - t348 - t349) + t108 * (t451 - t450) + t108 * (t148 * (-t3
     #77 - t378) + t327 * (-t93 + t373)) * epinv) + t429 * t84 + t429 * 
     #(t507 * t327 * (t356 * t397 - t375 - t408) + t46 - t51 + t497) * t
     #12)
      t117 = -t92 + t298
      t148 = -t472 + t316
      t196 = -t407 + t294
      t216 = t212 * t554
      t67 = t259 * (t429 * (t12 * t47 - t26 * t345) + t649 * (t17 + t19)
     # * s56 + t108 * t26 * (epinv * (-t103 * t159 + t142 * t444) + t327
     # * (-t456 - t328 - t216)) + t26 * (-s34 * t429 * t562 + t350) * t4
     #31) + t259 * (t26 * (t18 * (t429 * (s56 * t67 + t382) - t485) + t4
     #29 * (-s56 * (t30 + t31) + s34 * t343 * t431) + t108 * (-t103 * t2
     #54 + t257 * t444 + t28 * (-epinv * t327 * t554 + t137 * t99))) + t
     #429 * (t136 * (-t484 + t308) + t137 * t117 - t282 * (-s25 + s56 + 
     #t551) - t313 * t548 + (t317 - t154) * t148 + (t295 + t146) * t196 
     #+ (t232 - t13 - t15 - t16) * t430 * s56) * t12)
      t42 = t521 * (t435 * (t42 * t26 * s56 + t192 * (s25 - s56)) + (s16
     # * (-t217 + t325 - t392 - t215) + t284 + t383) * t324 * propW16)
      t22 = -t347 * t42 + t405 * t67 + t457 * (t84 + t22 + t445 * (t12 *
     # (epinv * (-t136 * t361 - t137 * t334 + t625) - t404 * t304) + t26
     # * (t431 * (s34 * (-t342 - t344 + t352) - t346 * (t106 + s16 - s34
     # + s26 + s56) + t354 * t543 + t379) + t432 * (t121 * t23 + t151 * 
     #(t198 + t359) + t218 * t546 - t245 * t275 + t327 * (t203 - t206) -
     # t35 * (-epinv * t121 + t526) + t52 * (epinv * (t2 * (t609 + t89) 
     #+ t3 * (t641 - t90 + t87) + s13 * (t134 - t529) + s25 * t193 + t53
     #4 * s26 + s56 * t534 + t617 - t85 + t63 * s56) + t421))) + t645 * 
     #s56 + t629 * (t253 * (-t3 * (t90 - t87) + t243 + (t75 - s25) * s25
     # + (t191 + s12) * s26 + (-t54 + s16 + t4) * s34 + (t574 - t54 + t6
     #3) * s56 + t184 + t522) + t261 * (epinv * t59 + t526) + t59 * t83)
     #)) - t458
      t42 = (t332 + t363) * t18 - t294 - t316
      t59 = t333 + t337 + t360
      t67 = t156 - t293
      t84 = t3 * (t18 * t59 - t282 - t298 - t303 - t308)
      t99 = t24 * t42
      t103 = t45 * t67
      t108 = -t2 + t132
      t121 = t503 * t18
      t134 = -t2 * (t391 + t436) + t228 * t313 + t24 * (t128 * t282 + t1
     #8 * (t470 + t363 + t337) - t303 - t316 + t483 - t489) + t280 * t60
     # + t3 * (t104 * t627 + t18 * (-t592 + t476 + t333) + t430 * (-t18 
     #* (t624 + t504 + t505) + t595 + t600 + t601) - t298 - t306 - t423 
     #- t492 + t516) + t314 * t76 + t45 * (t18 * (t599 + t360 + t331 + t
     #332) - t293 - t294 - t308 - t509 + t591) + t108 * t498 + t18 * t33
     #4 + t285 * (-t2 + t590) - t315 * t626 - t430 * (-t18 * t597 - t121
     # + t594 + t596) - t296 + t340 + t361
      t39 = t39 * t430
      t16 = (t13 + t15 + t16) * t430
      t137 = t430 * (-t1 - t14 - t452 + t37)
      t151 = t18 * (t315 + t285)
      t156 = (-t334 - t361) * epinv
      t67 = t67 * t12
      t191 = t45 * (t67 - t44)
      t192 = t475 * t332
      t193 = epinv * (t52 + t261)
      t226 = t481 * t28
      t232 = t475 * t331
      t243 = t475 * t32
      t259 = t52 + t261
      t275 = t18 * t285
      t295 = t18 * t315
      t304 = t18 * t362
      t327 = t12 * t282 + t7 * (t280 + t314) + t7 * (-t275 - t295 + t436
     #) + t10 * (epinv * t259 - t18 * t28 + t212 + t253 - t286 - t290 + 
     #t298 + t303 + t308 - t333 - t337 + t365 - t484 + t83) + t10 * ((-t
     #365 + t333 + t337) * epinv - t309 + t463 + t304)
      t11 = -t2 * t487 * (t285 * t481 + t315 * t481 + t261 - t280 - t314
     # - t436 + t52) + t228 * t7 * (t486 - t212 - t83 - t253 + t286 + t2
     #90 + t294 + t309 + t316 + t192 - t193 + t226 + t477 + t480 + t482)
     # - t24 * (t7 * (-t303 - t308 - t298 - t261 - t52 + t404 + t92 + t4
     #84) + t10 * (t316 + t294 - t30 - t31 + t538 + t281 - t407 - t472))
     # - t3 * t327 - t45 * (t487 * (epinv * (t333 + t337 + t365) + t290 
     #+ t303 + t309 + t477 + t482) + t9 * (t7 * (t360 * t475 - t193 - t2
     #12 + t226 - t253 + t286 + t298 + t308 - t333 - t337 - t365 - t83) 
     #+ t8 * (t32 * t481 + t232 - t29 + t293))) - t488 * t7 * (t331 * t4
     #81 - t293 + t9 * (t486 + t30 + t31 + t294 + t316 + t192 + t478 + t
     #479)) - 18.0_dp * t487 * (t243 + t29 + t293 + t232) + 24.0_dp * t7
     # * (t243 + t29) + 16.0_dp * t7 * (t479 + t30 + t31 + t478) + t10 *
     # (t391 + t314 + t280 + t261 + t52 - t275 - t295 + t436) + t11
      t29 = (-t18 * (t362 + t365 + t366) + t28 + t286 + t290 + t309) * t
     #26
      t32 = (t18 * t259 - t23 - t253 - t36 - t83) * t26
      t36 = t24 * t26 * (t97 - t30 - t31)
      t92 = -t45 * t44
      t97 = -t3 * (t12 * (t18 * (t360 + t363 + t332 + t333 + t337 - t39)
     # + t16 - t294 - t298 - t303 - t308 - t316) + t29) - t12 * (t430 * 
     #(-t1 - t14 + t37) - t151 + t280 + t314) - t32 + t36 + t92
      t17 = t24 * (t12 * (t18 * (t360 + t470) - t308 - t489) + t26 * (t5
     #38 - t30)) - t3 * (t12 * (-t306 + t314 + t316 + t294 + t298) + t26
     # * (t31 - t212 - t83 + t286) + t12 * (-t18 * (t476 + t363 + t332 +
     # t333) + t315 * t459 + t430 * (t18 * (t624 + t504 + t505 - t17) + 
     #t13 - t595 - t600 - t601) + t568 * t627 + t492 + t591 + t592) + (-
     #t365 + t28 - t34) * t18 * t26) + t45 * (t12 * (t242 - t509) - t44)
     # + t494 + t12 * (t339 * (-t285 - t498) + t361 * (1.0_dp + t126) + 
     #t430 * (t18 * (t597 + t503 + t19) - t15 - t21 - t594 - t596) - t28
     #0 + t296 - t391 - t436 - t483 + t643) + t26 * (t128 * t261 - t290 
     #+ t463 - t52)
      t126 = s46 + s45 - s12
      t192 = t3 * t126
      t193 = t192 + t317
      t226 = t543 + t70
      t232 = -t152 + t549 + t57
      t242 = s34 + s46 + s45 - s12
      t243 = s34 - s46 - s45 + s12 - t154
      t259 = -s56 * t228 - t126 * t76 + s25
      t275 = t24 * t126
      t293 = t45 * s56
      t295 = t24 * s56
      t126 = t2 * t126
      t327 = -s12 - s13 + s25 + s45 + s46 + s56
      t331 = -s12 + s25 + s26 + s45 + s46 + s56
      t339 = t3 * t548 + t152
      t366 = -t3 * (s13 - s46 - s45 + s12) - s16 + s25 - s26 + s34 + s56
      t397 = t153 + t511
      t407 = t3 * t153
      t408 = t407 + s12 + t557 + t213
      t415 = s12 + s13 + s16 + s26 - s34 - s45 - s46
      t416 = t2 * t608
      t104 = t104 * s12
      t128 = t128 * s12
      t418 = s13 * t616
      t421 = t572 * t85
      t422 = -t231 + t115
      t423 = (-t146 + t72) * s12
      t411 = t5 * (-(-(-t225 + t115) * s26 + s12 * t636 - (t179 - t54) *
     # s25 + t60 * t633 - t419) * s26 + ((t63 + t4) * s16 + t605 - t614)
     # * s34 + (s16 * (-t146 + s34 - t150 - s16) + t615) * s46 - (-(t4 -
     # t62) * s56 + (-s12 - t398 - t536) * s25 - t2 * t635 + t24 * t550 
     #* s16 + t3 * (-t618 + t527) + t80 * s16) * s56 + t421 - t611 - t61
     #2 - t604 * (-t54 + s12 + t72) + t411)
      t63 = -t7 * propW34 * (t3 * (-t24 * t86 + t3 * (-t422 * t89 - t421
     #) + (s25 * (t70 + t231 - t71) + t610 + t90) * s16 + (s13 * t422 + 
     #(-s25 * t488 + t573) * s25 + t76 * t94 + t423) * s26 + ((-t79 - t3
     #98 - t536) * s16 + t24 * t607 - t605) * s34 + (s16 * (-s34 + t63 +
     # t115 + s16) - t115 * t136) * s46 + ((-t179 + t413) * s56 + s25 * 
     #(-t227 + t413) + (17.0_dp * s16 - 14.0_dp * s25) * s26 - t24 * t61
     #8 + t45 * t620 - t60 * t635 + t423) * s56) + t411) + t465 * (t3 * 
     #(-t616 * t89 + t86) + (s25 * (t179 + s12) + t94 + t418 - t4 * s16)
     # * s26 - (t3 * (t607 + t94) + (-s25 - s13 + s26 + s12) * s16) * s3
     #4 + (-s16 * (t122 + t146 - s16) + t615) * s46 + (-t613 * s56 + (t7
     #7 - t62) * s26 + t3 * (-s16 - s34 + s46) * s25 + t60 * t85 + t418 
     #- t559 + t623) * s56 - t563 + t611 + t612 + t421)
      t86 = s56 * t60 + t153 * t2 + s12 + t113
      t72 = t221 * t3 + t24 * t583 - s16 * t75 + (t79 - t223) * s25 - (-
     #t179 - t79 + s12 + t71) * s26 + (-t113 - t79 + t223 - t115) * s34 
     #+ (s16 + s45 - t54 + t165 + t213 + s46) * s46 + (t163 * t45 - t24 
     #* t68 + t524 + t525 + t72) * s56 - t174 + t194 + t195
      t75 = s13 - s46 - s56
      t122 = t407 + t557 + t213
      t194 = t550 * t122
      t195 = s16 + s56
      t149 = -(t77 + t229) * s34 + t2 * ((s46 + s25) * s25 + t169) - t24
     # * ((s46 - s34) * s34 - t171 + t174) + t3 * s16 * t149 - t45 * (s3
     #4 * t195 - t641) + t60 * s25 * t195 + s12 * (t24 * t75 - t211 - t2
     #29 + t234) + s13 * (-t79 - t4 + t508)
      t169 = (t246 - s16) * s16
      t171 = -t2 * t264 + s12 - t115 + t54
      t56 = (t56 + t4) * s26
      t195 = (-t3 * (s25 - s34 - s12) - t233) * s46
      t211 = t557 * (-s25 - s13 + s34 + s46 + t70 + t150)
      t221 = (t54 + t66) * s16
      t219 = t219 * s26
      t227 = t24 * (t88 - t90 + t87 + t89 - t85)
      t145 = epinv * (-t2 * t87 - t188 + s26 * t145 - t240 * s34 + s46 *
     # t240 - (-t157 * t3 + s25 + t536) * s56 + t158 - t241) - t158 - t2
     #73 + t274 + t276
      t157 = -15.0_dp * s25
      t240 = t228 * s34
      t62 = t228 * t90 - t24 * t268 - t45 * t88 + t222 + t237 + (s16 * t
     #523 - t224 - t66) * s25 - (s12 * t488 + s13 + t157 - t573) * s26 +
     # (-s16 * t488 - 17.0_dp * s26 + t157 + t566 + t70) * s34 + (-t235 
     #+ s26 + t77 - t240 + t62) * s46 + (s25 * t523 - s34 * t488 + s56 +
     # t225 - t82 + t98) * s56
      t82 = epinv * t62 - s25 * t262 - s34 * t263 - s46 * t266 - s56 * t
     #267 + t269 + t270 - t91
      t91 = -s25 - s13 - s16 + t220
      t157 = -t207 + t209
      t207 = -s25 - s16 + s34 + s26 + t70
      t209 = -s25 - s13 - s16 + s34 + s46 + t249
      t237 = (s13 - t208 - s25) * s25
      t241 = t3 * (t88 + t89)
      t246 = t3 * t268
      t61 = -t2 * (t88 - t87 - t85) + t24 * t540 + (-t70 + t231 - t223) 
     #* s25 + (t225 - t54 - t220 + t260) * s26 + (-t225 + s12 - t80 - t7
     #7 + t71) * s34 + (t2 * t61 + t265 - t508) * s46 + (t136 * t76 - t3
     # * (s13 - s46 + s12) + t231 - t240) * s56 - t246
      t76 = -s16 + s34 - s46
      t77 = t24 * t163
      t80 = -t122 * t152 + t155
      t88 = t3 * t140
      t122 = t88 + t126 + s26
      t155 = t2 * t341 + t136
      t231 = -t2 + epinv
      t240 = -t241 - s26 * t91 + (-t157 + s34) * s34 - s46 * t207 + (-t2
     #09 - s56) * s56 - t169 - t237
      t249 = s25 - s16 + s34 + s56
      t260 = t249 + t70
      t153 = t153 + t250
      t75 = t198 * t408 + t200 * (-s12 * t231 + t102 * t2 + t109 + t113)
     # - t206 * t397 - t218 * t149 + t23 * (-t174 * t3 + s12 * (t3 * t75
     # - t213 - t550) + t194) + t245 * t240 - t258 * (-t3 * (t89 - t85) 
     #+ t96 + epinv * t149 + (-t54 + s16 + t220) * s25 + (s13 + t220 + t
     #400) * s26 + (-t58 + t208 + t150 - s34) * s34 + (-t69 + t146 + t66
     #) * s46 + (-t3 * t310 - s56 - t150 + t420) * s56 + t169) - t261 * 
     #t82 + t271 * t145 - t329 * t86 + t35 * ((t3 * (-epinv * t118 + s26
     #) + s46 + s56 + t107 + t129 + t130 - t133 - t560 + t141) * s12 + t
     #194 * epinv) - t52 * (epinv * (t3 * (s34 * t171 - t195 - t210 + t5
     #6 + t94) - t211 - t221 - t219 - t227) + (-s13 + t57 + s25) * s25 +
     # (t3 * (s16 + s26 - s12) - s13 - s34) * s34 + t118 * (t317 + t66) 
     #+ t169 + t98 * (t247 + t66) + t96)
      t53 = t445 * (t12 * (-t243 * t296 + t285 * (-t341 * t590 + t126 + 
     #t249) + t315 * (epinv * t259 + t24 * t341 + s25 + s26) - t361 * (-
     #epinv * (-t247 + t621 + t557) + s45 + s46 + s56 + t128) + t430 * (
     #t548 * t596 + t597 * (-s25 - s46 - s26 - s56 - s45 + s12 + t551)) 
     #+ (t606 + t557 + t275) * (-t404 + t303)) + t26 * (-t338 * t415 * t
     #638 + t232 * t463) + t432 * (t26 * t75 + t53 * (-t239 * ((t557 + t
     #536) * s56 - t3 * (t76 * s26 + s56 * t76) + t188 + (-t251 + t154 +
     # t150 + s12) * s12 + s13 * (-t511 + s25) + (-s16 + s34 - s46 + s56
     # - s25) * s25) - t277 * (epinv * (t241 + s26 * t91 + (t157 - s34) 
     #* s34 + s46 * t207 + (t209 + s56) * s56 + t169 + t237) - t174 + t1
     #81 - t182 + t183 - t184 + t185 + t186 + t187 + t94 + t222))))
      t8 = t434 * (t507 * t493 * t63 * t446 * t431 + (t12 * (t18 * (t430
     # * (-t136 * t19 + t339 * t624 - t503 * t548) - t514) - t306 * (-t1
     #92 + t247 - t557) + t314 * t259 - t334 * (epinv * t243 - s12 + s34
     # + s45 + s46) - t430 * (t433 + t512) + t483 * (-t69 - t293 - t275)
     # + t498 * (-epinv * (t24 * t608 + t293 + t69) + s12 * t108 + s34 +
     # t416 + t154) - t217 + t507 * (t244 * t430 * t72 - t375 * t408)) +
     # t26 * (-t232 * t290 + t431 * (t18 * (t327 * t368 - t366 * t372) +
     # t376 * t408 * t432) + t379 * t415 * t638 + (t202 * t408 + t373 * 
     #t397 - t377 * t86) * t432 * epinv)) * propW34 * t8)
      t19 = t26 * t431
      t63 = t629 * t212 * (-t24 * (t541 + t90 + t85) + t246 + t96 + (-t2
     #25 + s12 + t71) * s25 + (t136 * t228 - s12 + t225 - t71) * s34 + (
     #-t79 + t66 - t77) * s46 + (-t163 * t60 - t233 * t3 + s56 + t57 - t
     #68) * s56 + t113 * s12)
      t13 = t445 * (t12 * (t18 * (t122 * t332 + t155 * t360 + t193 * t36
     #3 + t341 * t476) + t430 * (-t13 * t136 + t331 * t600 + t548 * t595
     #) + t385) + t26 * (t18 * (-t260 * t34 - t272 * t336) + t431 * (-t3
     #27 * t350 + t345 * t415 - t205 + t238) + t432 * (t142 * (epinv * t
     #80 - t158 - t273 + t274 + t276) + t159 * (-t3 * t161 * s34 + epinv
     # * t61 + t3 * ((-t162 - s56) * s56 + t174 - t181 + t182 - t183 + t
     #184 - t187 - t90 - t94) - t188) + t204 * t408 + t254 * t61 + t257 
     #* t80 + t28 * (epinv * (-t173 * t60 - t228 * (s34 * t136 - t641) +
     # t24 * (-t119 + t123 + t90 + t85) + t3 * ((t118 + t98 + t165 + s16
     #) * s16 - t172) + s13 * (-t79 + s56 - t77) - t120) + (t3 * (s25 + 
     #s46 + t107) + t105 - t106 + t111 + t112 - t114 - t560 + t124) * s1
     #2 - t125) + t328 * t397)) - t63)
      t13 = t445 * (t26 * (t289 * (-t474 + t286) + t312 * (-t538 + t30) 
     #+ t417 * t408 * t432 + t18 * (t454 - t453) * t431) + t18 * (t136 *
     # t168 + t242 * t333) * t12) + t445 * (t12 * (epinv * t389 - t193 *
     # t316 - t242 * t298 - t282 * (-t132 * t341 + t140 + t192) - t294 *
     # t122 - t308 * t155 + t341 * (-t492 + t516) - t591 * (-t126 + s25 
     #- t295) + t592 * (epinv * t341 + s12 - s34 - s45 - s46) - t627 * (
     #epinv * (-t416 + s25 - t295) + t104 + t558) - t555 * (t331 * t505 
     #+ t504 * t548)) + t26 * (t226 * (-t473 + t320) + t260 * t31 + t272
     # * (-t304 + t301 + t309) + t18 * (t327 * t369 - t382 * t415) * t43
     #1 + t456 * t397 * t432)) + t13
      t60 = t445 * (t12 * (t144 * t341 + t167 * t548) + t26 * (t175 + t1
     #78))
      t61 = t248 * t341
      t63 = t61 * (-0.5e1_dp / 0.9e1_dp * t280 - 0.8e1_dp / 0.9e1_dp * t
     #313)
      t8 = t347 * t60 - t405 * t13 + t457 * (t445 * (t12 * (t393 * t496 
     #+ t46 + t497 + t501 - t51) + t19 * (t18 * (t353 - t409) + t425 + t
     #48)) + t8 + t53 + t445 * (t12 * (t430 * (t136 * (t15 + t21) - t339
     # * t601 + t548 * t594) - t392 + t625 * (-1.0_dp + t100)) + t26 * (
     #t431 * (t327 * (-t348 - t349) - t346 * (epinv * t327 + s12 + s13 +
     # s16 + s26 - s34 - s45 - s46) + t354 * t366 - t438 - t439 + t440) 
     #+ t432 * (-t153 * t203 - t255 * (epinv * t153 - t229 - t256) + t31
     #8 * t359 + t326 * t397 - t62 * t83)) + t403 * t12 * t432 * t72 + t
     #432 * (t26 * t427 * t431 - t426 * t644) * t408 * t323 + t629 * t25
     #3 * (-t3 * t171 * s34 - t3 * (-t210 + t56 - t195 + t94) + t211 + t
     #219 + t221 + t227))) - t542 * t61 * (t18 * t470 - t489) + t545 * t
     #248 * (t47 + t648) + t63
      t13 = t430 * t40
      t15 = -t3 * (t12 * (-t313 + t13 - t462) + t29) + t36 + t92 - t12 *
     # ((-t334 + t340 - t361) * epinv + t137 - t296 - t306 + t391) - t32
      t13 = -t3 * (t12 * (t18 * (-t363 - t332 + t282) + t13 + t294 - t31
     #3 + t316) + t29) + t191 + t36 - t12 * (-t306 - t296 + t156 + t137)
     # - t32
      t21 = t506 + t234
      t29 = t3 * t506 + s34
      t32 = (s25 + s34 + s45 + s46 - s16) * s16
      t36 = t608 * t69
      t51 = t152 - t235 + t79 + t77
      t53 = t164 - t66 + t166
      t56 = s25 + s16 - s34 + s46 + s56 + s12
      t60 = t3 * (s13 + s26) + t56
      t61 = t56 + t98 + t223
      t56 = t56 + t98 + t412
      t62 = t201 - t54
      t63 = s25 - s16 + s46 + s26 + s56 + s45 - t223
      t68 = t550 - t70
      t69 = s16 - t66 + t166
      t72 = t567 + t54
      t75 = -t108
      t76 = t75 * s12
      t80 = t587 + t166 - t57
      t82 = t163 - t70
      t86 = s12 + s26 + s56
      t90 = s46 + s26 + s56 + s45
      t91 = t2 * t152
      t92 = t608 + t91 + t54 + t115
      t96 = epinv * t86
      t75 = t75 * s13
      t100 = t152 * t24 + t3 * (s13 + s34 + s46 + s45) - s16 + t179
      t102 = t489 * t29
      t29 = t470 * t29
      t53 = t445 * (t12 * (t556 + t47 + t501 + t502) + (-t603 + t455 + t
     #584 + t311 + t322) * t18 * t26) + t445 * (t12 * (t18 * (-t360 * t9
     #2 - t29) + t430 * (-t537 - t535) + t561) + t26 * (t135 * (t281 - t
     #31) + t190 * t539 + t431 * (t278 * (t562 - t343) + t307 * (t18 * t
     #579 - t345 - t350 - t354)) - t302 - t305 - t355 - t388)) + t445 * 
     #(t12 * (-s13 * t516 + t117 * (t588 - s16 + s46 + s26 + s56 + s45) 
     #+ t148 * (t3 * t90 + t115 - t65) + t196 * (t91 + s13 + t398 + t621
     #) + t21 * (-t18 * t476 + t492) - t282 * (t60 * epinv - s13 - s26 -
     # s45 - s46 - s56 - t4) + t308 * t92 - t313 * t60 - t592 * (-t234 -
     # t517) + t102) + t26 * (-t142 * (s12 - t96 + t146) - t159 * (epinv
     # * t53 - t321) - t212 * t69 - t254 * t53 + t257 * t86 - t28 * (epi
     #nv * t69 + t279) + t205 * t431))
      t60 = s34 * (-t18 * intHLs16s25s26s34s56x1411D6eps1() + intHLs16s2
     #5s26s34s56x1411D6eps0())
      t69 = t144 * t506
      t91 = t200 + t202 + t255
      t32 = t445 * (t12 * (t280 * t61 + t430 * (t356 * (t18 * t598 - t18
     #9 + t197) + t375) + t404 * t100) + t26 * (epinv * (-t373 + t374 + 
     #t377 + t378) - t23 * t68 + t231 * t91 + t51 * t83 + t204 - t326 + 
     #t328 + t329) + t385 * t564 + t389 * (epinv * t564 - t390) - t490 *
     # (-t7 * (-t24 * t406 - t3 * (t252 - t283 + t287 + t32 - t642 + t87
     # + t89) + t45 * (t27 * t49 * t519 - t414) + t5 * (-(-s25 + s13 + s
     #16) * s16 - (-t146 + s13 - s26) * s26 - (s25 - s13 - s16 + s26) * 
     #s34 - (s13 + s34 - s46 - s45 - t402 - s56) * s56 + t85 + t36 - t62
     #8 * t519 * t49)) + t10 * (-t252 + t283 - t287 - t32 + t642 - t87 -
     # t89 + t85 - t367)) * t431)
      t85 = (t426 + t593) * t323
      t92 = (-t632 + t601) * t602
      t105 = t445 * (epinv * (t12 * t513 + t330 * t564) + t26 * (t431 * 
     #(t278 * (t401 - t342 - t344) + t307 * (t639 - t348 - t349) + t630 
     #* t95 + t647) + t307 * (-t379 + t338) * t638) - t437 * t12 + t644 
     #* (t58 * (-t594 - t596) + t230 + t512 + t85 + t92))
      t107 = t445 * (t12 * (t430 * t530 + t500) + t26 * (t292 * (-t463 +
     # t290) + t198 + t203 + t206 + t359) + t147)
      t4 = t107 + t105 + t32 + t445 * (t12 * (t285 * (epinv * t61 - s26 
     #- s45 - s46 - s56 - t4 - t577) + t296 * t62 - t303 * t100 + t306 *
     # t63 + t314 * t56 + t315 * (epinv * t56 - t2 * (s13 + s26 + s56) -
     # s45 - s46 - t398) + t334 * (epinv * t62 + s13 + s34) - t340 * (ep
     #inv * t72 - s13 - s25) + t361 * (epinv * t63 + s13 + s25 + s26 + s
     #56) - t391 * t72 - t394 * (t514 + t552) - t498 * (t3 * t358 + t75)
     #) + t26 * (t218 * (s16 + t166 - t57) + t239 * t86 - t245 * t82 + t
     #253 * t80 + t258 * (s16 * t323 - t3 * (-s25 * t459 + s26 + s56 - t
     #560) + t76) + t261 * (epinv * t51 - t288 + t291) + t271 * (-s25 - 
     #s12 + t96) - t277 * (epinv * t82 + s12 - s16 - s26 + s34 - s56) - 
     #t35 * (epinv * t68 + s12 + s34) + t52 * (epinv * t80 + t549 + t66)
     #))
      t32 = -s13 + s25 + s26 - s34 + s45 + s46 + s56
      t51 = s13 - s45 - s46
      t56 = s12 + s16 + s25 + s26 - s34 - s45 - s46
      t57 = t2 - t443
      t61 = t24 * t429
      t62 = propZ25 * t57 - t61
      t63 = t2 + t443
      t66 = t3 * t429
      t68 = propZ25 * t63 - t66
      t63 = -propZ25 * t63 + t66
      t61 = -propZ25 * t57 + t61
      t66 = propZ25 * (-t386 - t2 - t443) + t429 * (t386 + t3)
      t72 = propW16 * s34 + 1.0_dp
      t80 = t139 * t201
      t82 = -s13 + s56 + s12 + t402
      t86 = -s13 + s16 - s34 + s56 + s12 + t402
      t96 = -s13 + s16 - s34 + s12 + t154 + t138
      t9 = t428 * (-t57 * t9 + t24) + t435 * (t26 * (t278 * t429 * t43 +
     # t300 * t429 * t469 - t177) + t67 * t429 * t565) + (t254 + t83 + t
     #212) * t62 * s12 + t429 * (s16 * (t215 - t325) - t284 - t383) * t3
     #24 * propW16
      t43 = t429 * t435
      t36 = t43 * (t12 * ((t32 * t334 - t340 * t56 + t358 * t498 - t361 
     #* t51) * epinv + t280 * t82) + t431 * (t26 * (t18 * (-t226 * t372 
     #+ t272 * t368 + t278 * t380 - t371 * t95) + t346 * (-t106 - s16 + 
     #s34 - s26 - s56 + t357) + t376 + t647) + t490 * (t7 * (t3 * (-s16 
     #* (t2 * t27 + t608) + t642 + t80) - t5 * (s16 * (-t3 * (-propW16 *
     # t87 - propW16 * t89 + s25 * (propW16 * t139 + 1.0_dp) + s26 * (-p
     #ropW16 * t154 + t72) + s56 * t72) - t608) + t642 + t80 + t628 * pr
     #opW16 * t94) + t386 * t27 * t49) - t10 * ((-t160 - s26) * s26 + s3
     #4 * t576 + (-s25 + s13 + s34 - s46 - s45 - t98 - s56) * s56 + t184
     # - t36))) + t379 * t26 * t300 * t638)
      t32 = t43 * (t12 * (-t18 * (t393 + t514 + t552) + t285 * (epinv * 
     #t82 + s13 - s25 - s26 - s45 - s46 - s56) + t296 * t32 - t306 * t51
     # + t314 * t96 + t315 * (epinv * t96 + s13 - s25 - s26 - s45 - s46 
     #- s56) + t358 * t483 - t391 * t56) + t26 * (t245 * (t70 + s16) + t
     #258 * (t544 + t76) - t271 * (s25 + t357) + t277 * (s16 * t459 + s2
     #6 - s34 + s56 + t128) + t35 * (-s34 + t128) + t394 * t91 + t431 * 
     #t48 + t52 * (t549 + t76) + t442))
      t32 = t32 + epinv * (t330 * t66 + t50 * t68 + t513 * t63 + t62 * (
     #t374 + t377 + t378)) + t389 * (epinv * t66 - t387 * t429) + t500 *
     # t63 + t62 * (t328 + t329 + t204) + t66 * (t385 + t392) + t68 * (t
     #46 + t47 + t501 + t502) + t26 * t43 * (t431 * (t226 * t354 + t272 
     #* (-t348 - t349) + t278 * (-t342 - t344) + t323 * (t384 + t427) + 
     #t352 * t95 - t441 * t300) - t463 * t292) + t36 + t239 * s12 * t61
      t36 = s25 - s16 + s26 + s56 - t54
      t50 = t152 + t79 + t77
      t51 = -s46 + s26 + s56 + t146 - t224
      t54 = t24 * t58
      t56 = -s16 + t529 + t54
      t57 = s25 - s46
      t63 = t57 - t223
      t66 = s25 - s16 + s34 - s46 - t412
      t67 = t587 + t166
      t70 = -s13 - s34 + s46 + s26 + s45 + t88
      t72 = s25 * t27
      t76 = t448 + t98
      t77 = t448 + t150
      t5 = t445 * (t12 * (-t303 * t56 + t306 * t36 - t314 * t51 - t315 *
     # (epinv * t51 + t511 + t54) + t361 * (epinv * t36 + s13 + s25 + s2
     #6 + s56) + t430 * (t555 * t598 * t70 + t375) - t483 * (t358 - t71)
     # - t498 * (t358 * t568 + t75)) + t26 * (epinv * (-t200 - t202 - t2
     #55 - t373) + t253 * t67 + t431 * (t18 * (s26 * t368 - t372 * t77 -
     # t410 * t76) - t376) + t50 * t83 - t326 + t640) - t490 * (t7 * ((s
     #16 * t6 + (propZ25 * t152 * t45 - propZ25 * t116 - t24 + t5) * s25
     # - t152 * t45 + t234) * s16 + t72 * t25) - t10 * ((t586 - s16) * s
     #16 - t72)) * t431)
      t5 = t445 * (t12 * (t430 * (t58 * (t121 - t594 - t596) + t230 + t5
     #12 + t85 + t92) + t496 * (-t393 - t514) - t46 - t497 - t532 + t533
     #) + t26 * (-t431 * (t425 + t48) - t198 - t203 + t206 - t359)) + t5
     # + t445 * (t12 * (-t280 * t66 - t285 * (epinv * t66 + t223 + t549)
     # + t340 * (epinv * t63 + s13 + s25) + t391 * t63 - t496 * t552 + t
     #404 * t56) + t26 * (epinv * (t261 * t50 + t52 * t67 - t550 * (t35 
     #+ t277) + t646 + t93) + t431 * (-s26 * t447 + t346 * (s25 + s16 - 
     #s34 + s56 - t622 + t98) + t354 * t77 - t399) - t550 * (t23 + t245)
     # + t619 + t338 * t76 * t638) + t449 * t12 * t70)
      t6 = t248 * (t167 * t58 + t60 + t69)
      t5 = -t347 * t6 + t405 * (t445 * (t12 * (epinv * (t330 + t389 + t5
     #13) + t18 * (t21 * t476 + t29 + t528) + t430 * (t506 * (-t18 * t50
     #5 + t600) + t58 * (-t214 + t595)) + t74 * (t484 - t308) + t385 + t
     #392 + t500 - t510 + t515) + t26 * (-t18 * t431 * (s26 * t369 + t38
     #2 * t76) + t328 + t329 + t204)) + t445 * (t12 * (t199 * (s25 + s13
     # - s16 + t234) - t21 * t492 - t282 * (s25 * t568 - t110 + t127) - 
     #t313 * (t57 - t71) - t591 * (t27 - t223) - t592 * (s34 * t568 + t5
     #17) - t627 * (t27 * t568 + t101) + (-t65 + t531) * (t472 - t316) -
     # t102) + t26 * (epinv * (t28 * t554 + t374 + t377 + t378) + t152 *
     # (-epinv * t142 - t257) + t431 * (s26 * t350 + t345 * t76) + t571 
     #* (epinv * t159 + t254) + t216))) + t457 * t5 - t545 * t248 * (t47
     # + t501 + t502)
      t6 = struc4PP * (-t84 - t430 * (t430 * ((-t395 + t396) * t18 + t18
     #9 - t197) + t452) + t285 - t296 - t306 + t315 - t391 - t491 - t99 
     #- t103) - struc7PP * (-t3 * (t313 + t462) + (t340 - t315 - t285) *
     # epinv - t280 - t314 + t391) + struc9PP * t134 + (-t84 - t391 - t3
     #06 + t315 + t285 - t296 - t491 - t490 - t99 - t103) * (struc5PP + 
     #struc6PP)
      t1 = (-t3 * (t12 * t40 + t27 * t41) - t45 * t44 * t27 + t12 * (t1 
     #+ t14 - t37) + t38 * t27) * struc16PP * t430
      t1 = t445 * (t12 * t6 + struc11PP * t15 + struc12PP * t13 + struc1
     #3PP * (t24 * t12 * (t236 - t294 - t316) - t3 * (t12 * (t18 * (-t36
     #0 - t333 - t337 + t282 - t39) + t298 + t303 + t308 - t313 + t16) +
     # t41) - t12 * (-t306 - t314 - t280 - t296 + t156 + t151 + t137) + 
     #t38 + t191) - struc14PP * t11 + struc15PP * t97 + struc17PP * t17 
     #+ struc21PP * (t2 * t464 + t24 * t26 * (t18 * (t364 + t365) - t159
     # - t286 - t320) + t3 * (t12 * (-t18 * (t363 + t333) + t298 + t313 
     #+ t316 + t462) + t26 * (t18 * (t362 + t336 + t33 + t34 - t460) - t
     #258 - t261 - t30 - t301 - t309 - t31 - t52 + t461)) - t12 * (t391 
     #- t303 - t314 + t280 + t466 + t404) - t26 * (-t245 + t468 - t467) 
     #+ t471) - struc22PP * t335 + struc24PP * t495 + struc25PP * t73 + 
     #struc27PP * t499 + t1)
      t2 = t297 * t457 * (t24 * (t43 * t42 * t565 * t12 + s12 * t62 * (t
     #253 + t218)) + t3 * (t61 * (s12 * t23 + t198 + t203 + t359) + (t12
     # * (t27 * t81 + t59 * t634) + t26 * (t18 * (t300 * t431 * (t370 + 
     #t382) - t589 - t585 - t311 - t322 - t455 - t584) + t431 * (t272 * 
     #t350 - t300 * (t345 + t351) + t631))) * t429 * t435 + t257 * s12 *
     # t62 + t68 * s34 * (epinv * t592 + t516) + (t12 * (t27 * t591 + t2
     #82 * (epinv * t86 + s13 - s25 - s26 - s45 - s46 - s56) + t313 * t8
     #6 - t565 * (t298 + t303 + t308)) + t26 * (t142 * (t357 + t146) - t
     #159 * (t3 * t49 + t104) - t261 * (s16 + s26 + s56 + t104 - t234) -
     # t28 * (-t234 + t104) + t302 + t305 + t355 + t388 + t580 + t581) -
     # t19 * (t272 * t369 + t381) * t18) * t429 * t435) + t45 * t9 - t48
     #8 * t520 * t429 * t319 + t32) * struc35PP
      result = t405 * t1 + t5 * struc10PP + t299 * struc19PP + t64 * str
     #uc20PP + t55 * struc30PP + t22 * struc31PP + (-t131 * t405 - t78 *
     # t457 + t347 * t297 * (t435 * (t26 * (t178 * t429 + t176 + t177 - 
     #t428) + t180 * (s12 + s25 - s56) * t12) + t429 * (s16 * (t217 - t3
     #25 + t392 + t215) - t383 - t284) * t324 * propW16) - t458) * struc
     #32PP + t8 * struc33PP + (-t347 * t297 * (t435 * (t429 * (t167 * (t
     #90 + t146) - t518 + t60 + t69) * t12 + t170) + t429 * (s16 * (-t21
     #5 + t325) + t284 + t383) * t324 * propW16) - t405 * t53 + t457 * t
     #4 - t458 + t143 * t542 * t483) * struc34PP + t20 * struc36PP + t2

           ampNonresonantLightFullPP = result
       end function ampNonresonantLightFullPP

       function ampNonresonantLightFullReC4MM()
           implicit none
           complex(dp) :: ampNonresonantLightFullReC4MM

           complex(dp) :: result

      result = 0.0_dp

           ampNonresonantLightFullReC4MM = result
       end function ampNonresonantLightFullReC4MM

       function ampNonresonantLightFullReC4MP()
           implicit none
           complex(dp) :: ampNonresonantLightFullReC4MP
           complex(dp) ::  t1,t10,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109
           complex(dp) ::  t11,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t12
           complex(dp) ::  t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t13,t130
           complex(dp) ::  t131,t132,t133,t134,t135,t136,t137,t138,t139,t14,t140,t141
           complex(dp) ::  t142,t143,t144,t145,t146,t147,t148,t149,t15,t150,t151,t152
           complex(dp) ::  t153,t154,t155,t156,t157,t158,t159,t16,t160,t161,t162,t163
           complex(dp) ::  t164,t165,t166,t167,t168,t169,t17,t170,t171,t172,t173,t174
           complex(dp) ::  t175,t176,t177,t178,t179,t18,t180,t181,t182,t183,t184,t185
           complex(dp) ::  t186,t187,t188,t189,t19,t190,t191,t192,t193,t194,t195,t196
           complex(dp) ::  t197,t198,t199,t2,t20,t200,t201,t202,t203,t204,t205,t206
           complex(dp) ::  t207,t208,t209,t21,t210,t211,t212,t213,t214,t215,t216,t217
           complex(dp) ::  t218,t219,t22,t220,t221,t222,t223,t224,t225,t226,t227,t228
           complex(dp) ::  t229,t23,t230,t231,t232,t233,t234,t235,t236,t237,t238,t239
           complex(dp) ::  t24,t240,t241,t242,t243,t244,t245,t246,t247,t248,t249,t25
           complex(dp) ::  t250,t251,t252,t253,t254,t255,t256,t257,t258,t259,t26,t260
           complex(dp) ::  t261,t262,t263,t264,t265,t266,t267,t268,t269,t27,t270,t271
           complex(dp) ::  t272,t273,t274,t275,t276,t277,t278,t279,t28,t280,t281,t282
           complex(dp) ::  t283,t284,t285,t286,t287,t288,t289,t29,t290,t291,t292,t293
           complex(dp) ::  t294,t295,t296,t297,t298,t299,t3,t30,t300,t301,t302,t303
           complex(dp) ::  t304,t305,t306,t307,t308,t309,t31,t310,t311,t312,t313,t314
           complex(dp) ::  t315,t316,t317,t318,t319,t32,t320,t321,t322,t323,t324,t325
           complex(dp) ::  t326,t327,t328,t329,t33,t330,t331,t332,t333,t334,t335,t336
           complex(dp) ::  t337,t338,t339,t34,t340,t341,t342,t343,t344,t345,t346,t347
           complex(dp) ::  t348,t349,t35,t350,t351,t352,t353,t354,t355,t356,t357,t358
           complex(dp) ::  t359,t36,t360,t361,t362,t363,t364,t365,t366,t367,t368,t369
           complex(dp) ::  t37,t370,t371,t372,t373,t374,t375,t376,t377,t378,t379,t38
           complex(dp) ::  t380,t381,t382,t383,t384,t385,t386,t387,t388,t389,t39,t390
           complex(dp) ::  t391,t392,t393,t394,t395,t396,t397,t398,t399,t4,t40,t400
           complex(dp) ::  t401,t402,t403,t404,t405,t406,t407,t408,t409,t41,t410,t411
           complex(dp) ::  t412,t413,t414,t415,t416,t417,t418,t419,t42,t420,t421,t422
           complex(dp) ::  t423,t424,t425,t426,t427,t428,t429,t43,t430,t431,t432,t433
           complex(dp) ::  t434,t435,t436,t437,t438,t439,t44,t440,t441,t442,t443,t444
           complex(dp) ::  t445,t446,t447,t448,t449,t45,t450,t451,t452,t453,t454,t455
           complex(dp) ::  t456,t457,t458,t459,t46,t460,t461,t462,t463,t464,t465,t466
           complex(dp) ::  t467,t468,t469,t47,t470,t471,t472,t473,t474,t475,t476,t477
           complex(dp) ::  t478,t479,t48,t480,t481,t482,t483,t484,t485,t486,t487,t488
           complex(dp) ::  t489,t49,t490,t491,t492,t493,t494,t495,t496,t497,t498,t5
           complex(dp) ::  t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t6,t60
           complex(dp) ::  t61,t62,t63,t64,t65,t66,t67,t68,t69,t7,t70,t71
           complex(dp) ::  t72,t73,t74,t75,t76,t77,t78,t79,t8,t80,t81,t82
           complex(dp) ::  t83,t84,t85,t86,t87,t88,t89,t9,t90,t91,t92,t93
           complex(dp) ::  t94,t95,t96,t97,t98,t99

           complex(dp) :: result

      t1 = intHLs160000x0112D2eps0()
      t2 = intHLs16s25s26s34s56x1111D2eps0()
      t3 = intHs160000x0112D2eps0()
      t4 = intHs160s26s34s56x1021D2eps1()
      t5 = intHs160s26s34s56x1022D4eps0()
      t6 = intHs16s25s26s34s56x1112D4eps0()
      t7 = intHs16s25s26s34s56x1121D4eps0()
      t8 = intHs16s25s26s34s56x1211D4eps0()
      t9 = intHs16s25s26s34s56x1131D4eps0()
      t10 = 4.0_dp
      t11 = propZ25 * s25
      t12 = -t10 + t11
      t13 = intHs16s25s26s34s56x1221D4eps0()
      t14 = intHs16s25s26s34s56x1311D4eps0()
      t15 = intHs16s25s26s34s56x1311D4eps1()
      t16 = 2.0_dp
      t17 = -t16 + epinv
      t18 = 3.0_dp
      t19 = gw ** 2.0_dp
      t20 = gb ** 2.0_dp
      t21 = t18 * s25
      t22 = t21 * t19 * propZ25
      t23 = -t20 * t12
      t24 = t23 + t22
      t25 = 1.0_dp - epinv
      t26 = t20 * t12
      t27 = t26 - t22
      t28 = intHLs16s25s26s34s56x1112D4eps0()
      t29 = t16 + t11
      t30 = t20 * t29 + t22
      t31 = intHs16s25s26s34s56x1310D4eps0()
      t32 = s16 + s26 - s34 + s56
      t33 = intHLs16s25s26s34s56x1222D6eps1()
      t34 = s12 + s25
      t35 = t16 * t34
      t36 = t10 * s26
      t37 = t36 + t35 + s16
      t38 = intHLs16s25s26s34s56x1211D2eps1()
      t39 = 6.0_dp
      t40 = 5.0_dp
      t41 = t18 * s26
      t42 = t40 * s25
      t43 = t39 * s12
      t44 = t42 + t43 + s34 - s56 + t41
      t45 = intHLs16s25s26s34s56x1121D2eps1()
      t46 = s25 + s26 + s56
      t47 = intHLs16s25s26s34s56x1112D2eps1()
      t48 = -s13 + s26 + s45 + s46 + s56
      t49 = t16 * t48 + s25
      t50 = intHs16s25s26s34s56x1411D6eps0()
      t51 = t16 * s12
      t52 = t51 + s25
      t53 = intHs16s25s26s34s56x1321D6eps1()
      t54 = t51 + t21
      t55 = intHs16s25s26s34s56x1231D6eps1()
      t56 = t10 * s12
      t57 = t56 + t21
      t58 = intHLs16s25s26s34s56x1132D6eps1()
      t59 = t16 * s26
      t60 = t59 + s25
      t61 = intHLs16s25s26s34s56x1312D6eps1()
      t62 = s25 + s16
      t63 = t16 * (s12 + s26)
      t64 = t62 + t63
      t65 = intHs16s25s26s34s56x1312D6eps1()
      t66 = t21 + t63
      t67 = intHs16s25s26s34s56x1213D6eps1()
      t68 = s12 + s25 + s26
      t69 = t16 * t68 + s34
      t70 = s25 + s26 - s56
      t71 = intHLs16s25s26s34s56x1123D6eps1()
      t72 = t16 * s25
      t73 = t72 + t41 - s56
      t74 = intHLs16s25s26s34s56x1213D6eps1()
      t75 = t41 + t35 - s56
      t76 = intHs16s25s26s34s56x1132D6eps1()
      t77 = t72 + t56 + s26 + s34 - s56
      t78 = intHs16s25s26s34s56x1123D6eps1()
      t79 = t16 * (s12 + s34) + t70
      t80 = intHs16s25s26s34s56x1110D2eps1()
      t81 = intHs16s25s26s34s56x1220D4eps0()
      t82 = intHs16s25s26s34s56x1130D4eps0()
      t83 = intHs16s25s26s34s56x1120D2eps1()
      t84 = t25 * s25
      t85 = intHs160s26s34s56x1020D2eps1()
      t86 = intHs16s25s26s34s56x1220D4eps1()
      t87 = intHs16s25s26s34s56x1130D4eps1()
      t88 = intHs16s25s26s34s56x1120D2eps0()
      t89 = intHs16s25s26s34s56x1210D2eps0()
      t90 = intHs16s25s26s34s56x1210D2eps1()
      t91 = intHs16s25s26s34s56x1122D4eps1()
      t92 = t17 * s25
      t93 = t16 * epinv
      t94 = intHLs16s25s26s34s56x1112D2eps0()
      t95 = intHLs16s25s26s34s56x1122D4eps1()
      t96 = s13 - s45 - s46
      t97 = t16 * s56
      t98 = t97 - t96
      t99 = intHLs16s25s26s34s56x1222D6eps0()
      t100 = intHs16s25s26s34s56x1212D4eps1()
      t101 = intHLs16s25s26s34s56x1211D2eps0()
      t102 = intHLs16s25s26s34s56x1212D4eps1()
      t103 = t16 * t96
      t104 = t18 * s56
      t105 = t56 - t104 + s25 + s26 - s34 + t103
      t106 = epinv * t105
      t107 = intHLs16s25s26s34s56x1212D4eps0()
      t108 = intHs16s25s26s34s56x1222D6eps0()
      t109 = intHLs16s25s26s34s56x1121D2eps0()
      t110 = intHs160s26s34s56x1022D4eps1()
      t111 = intHs16s25s26s34s56x1121D2eps1()
      t112 = 1.0_dp + epinv
      t113 = t51 * t112
      t114 = intHs16s25s26s34s56x1121D4eps1()
      t115 = intHs16s25s26s34s56x1211D4eps1()
      t116 = intHs16s25s26s34s56x1112D4eps1()
      t117 = -t18 + t93
      t118 = intHLs160000x0111D0eps0()
      t119 = intHLs160000x0111D0eps1()
      t120 = intHLs160000x0112D2eps1()
      t121 = intHLs16s25s26s34s56x1112D4eps1()
      t122 = intHs16s25s26s34s56x1122D4eps0()
      t123 = intHs16s25s26s34s56x1212D4eps0()
      t124 = intHs160s26s34s56x1031D4eps1()
      t125 = intHs16s25s26s34s56x1121D2eps0()
      t126 = intHs16s25s26s34s56x1131D4eps1()
      t127 = -t56 + t92
      t128 = intHs16s25s26s34s56x1221D4eps1()
      t129 = -t51 + t92
      t130 = intHs16s25s26s34s56x1141D6eps0()
      t131 = intHs16s25s26s34s56x1321D6eps0()
      t132 = intHs16s25s26s34s56x1231D6eps0()
      t133 = intHLs16s25s26s34s56x1132D6eps0()
      t134 = intHLs16s25s26s34s56x1312D6eps0()
      t135 = intHs16s25s26s34s56x1312D6eps0()
      t136 = intHs16s25s26s34s56x1113D4eps1()
      t137 = intHs16s25s26s34s56x1114D6eps0()
      t138 = intHs16s25s26s34s56x1213D6eps0()
      t139 = intHLs16s25s26s34s56x1113D4eps0()
      t140 = -t104 - s26 + t103
      t141 = intHLs16s25s26s34s56x1122D4eps0()
      t142 = intHLs16s25s26s34s56x1114D6eps0()
      t143 = intHLs16s25s26s34s56x1123D6eps0()
      t144 = intHLs16s25s26s34s56x1213D6eps0()
      t145 = intHs16s25s26s34s56x1132D6eps0()
      t146 = intHs16s25s26s34s56x1123D6eps0()
      t147 = intHs16s25s26s34s56x1310D4eps1()
      t148 = intHs16s25s26s34s56x1110D2eps0()
      t149 = intHs160s26s34s56x1020D2eps0()
      t150 = intHLs16s25s26s34s56x1113D4eps1()
      t151 = intHs160s26s34s56x1031D4eps0()
      t152 = 0.1e1_dp / s25
      t153 = t39 * t19 * propZ25
      t154 = t20 * (-propZ25 * t16 + 8.0_dp * t152) + t153
      t155 = intHs16s25s26s34s56x1211D2eps1()
      t156 = 0.1e1_dp / t32
      t157 = t135 * t66
      t158 = t138 * t69
      t159 = t145 * t77
      t160 = t146 * t79
      t161 = t131 * t54
      t162 = t132 * t57
      t163 = s34 * t136
      t164 = t34 * t86
      t165 = t52 * t87
      t166 = t126 * t127
      t167 = t125 * s12
      t168 = t60 * t133
      t169 = t64 * t134
      t170 = t139 * t140
      t171 = -t141 * t98
      t172 = t2 * t29
      t173 = t34 * t81 + t52 * t82
      t174 = t9 + t13 + t14
      t175 = t58 * t60
      t64 = t61 * t64
      t54 = t53 * t54
      t57 = t55 * t57
      t176 = t65 * t66
      t177 = t67 * t69
      t178 = t76 * t77
      t179 = t78 * t79
      t180 = t25 * t147
      t181 = t15 * t17
      t182 = t26 * t174
      t183 = t152 * t25
      t184 = -t7 - t8
      t185 = t101 * s34
      t186 = t109 * t46
      t187 = t37 * t99
      t188 = t49 * t94
      t189 = t121 * t117
      t105 = t105 * t107
      t190 = t25 * t90
      t191 = t149 * t152
      t192 = t44 * t108
      t193 = t1 * t29
      t194 = (t122 + t123) * t60
      t195 = t91 * (t93 * s26 - s26 - s34 + s56 - t56 + t92)
      t196 = t148 * t156
      t197 = t17 * t80
      t198 = t114 + t115
      t199 = t116 * t117
      t200 = t25 * intHs16s25s26s34s56x1222D6eps1()
      t201 = t200 * t44
      t202 = t47 * t49
      t203 = s34 * t38
      t37 = t33 * t37
      t204 = t37 * t25
      t205 = t25 * t110
      t206 = intHs160000x0112D2eps1() * epinv
      t207 = t111 * (t113 + s25)
      t208 = (t83 * (t51 + t84) + t197) * t156
      t209 = t151 * t154
      t210 = t24 * (t152 * (t205 - t208 - t206 - t207) - t155) + t27 * (
     #t152 * (-t199 + t201) + t156 * (t183 * t85 - t88 - t89) - t152 * t
     #198 * t17) + t27 * (t156 * (-t191 + t190) - t192 * t152) + (t95 * 
     #(t16 * (epinv * t98 + s26) + s25) - t102 * (-t59 + t106 - s25) - t
     #185 - t186 + t187 - t188 + t189 - t105) * t152 * t30 + t152 * (-t1
     #00 * ((-t18 + epinv) * s25 - t56 - t59 * t25) - t194 - t195 - t196
     #) * t24 + t152 * (t12 * t184 - (-t3 - t4 - t5) * t12 - t193) * t20
     # - t209 + t152 * ((-t45 * t46 + t119 - t120 - t202 - t203) * epinv
     # + t118 - t204) * t30
      t211 = t25 * intHLs16s25s26s34s56x1114D6eps1()
      t212 = t211 - t142
      t213 = t25 * intHs16s25s26s34s56x1141D6eps1()
      t214 = -t213 + t130
      t215 = t25 * intHs16s25s26s34s56x1114D6eps1()
      t216 = intHs16s25s26s34s56x1411D6eps1() * t25
      t217 = (-t215 + t137) * s34
      t218 = t214 * t52
      t174 = t174 * s25
      t219 = t212 * t30
      t220 = t19 * propZ25
      t221 = intHL0s25s26s34s56x1110D2eps0()
      t222 = t46 - t103
      t223 = intHLs160000x0111D2eps0()
      t224 = t220 * s25
      t225 = t224 + t20
      t226 = s25 + s16 - s34
      t227 = t226 + t63
      t228 = t16 * (s12 - s45 - s46)
      t229 = s26 + s16 - s34 + s13 - s56 + t228
      t230 = intHLs16s25s26s34s56x1411D6eps0()
      t231 = intHs160s26s34s56x1011D2eps0()
      t232 = s26 + s16 - s34 - s56
      t233 = t16 * (s12 + s13 - s45 - s46)
      t234 = t232 + t233
      t235 = t16 * (s12 + s13 - s45) + s26 + s16 - s34 - s46 - s56
      t236 = intHLs16s25s26s34s56x1311D4eps0()
      t237 = t16 * s45
      t238 = t237 + s25 + s26 - s13 + s46 + s56
      t239 = intHL0s25s260s56x1013D4eps1()
      t240 = intHL0s25s260s56x1010D0eps0()
      t241 = intHL0s25s260s56x1011D2eps0()
      t242 = intHL0s25s260s56x1020D2eps0()
      t243 = t10 * s56
      t244 = s46 ** 2.0_dp
      t245 = s16 - s34
      t246 = s45 + s46 + s56
      t247 = s56 ** 2.0_dp
      t248 = s45 + s46
      t249 = t10 * s56 * t248 + t16 * (s45 * t245 + (s25 + s16 - s34 + s
     #45) * s46 + s56 * t226) + t18 * (s26 * t246 + t247) - s12 * (s25 +
     # s26 - s13 + s45 - s46 + s56) + (-t16 * (s16 - s34 + s45 + s46) - 
     #s25 - t41 - t243 + s13) * s13 + (s26 + s16 - s34 + s45) * s25 + s4
     #5 ** 2.0_dp + t244
      t250 = intHL0s25s260s56x1022D4eps1()
      t251 = s26 + s56
      t252 = t72 + t251
      t253 = t16 * t246
      t254 = s13 * t252
      t255 = t251 * t246
      t256 = (-t253 - s26 - s25) * s25 + t254 - t255
      t257 = s16 + s26 + s56
      t258 = t18 * s34
      t259 = -t258 + t257
      t260 = t16 * s34
      t261 = -t260 + t257
      t262 = t16 * t248
      t263 = t232 - t262
      t264 = s13 * t259 + t261 * t263 + t51 * t261
      t265 = intHs16s25s26s34s56x1112D2eps1()
      t266 = t16 * (s13 - s45 - s46 - s56)
      t267 = -t266 + s25
      t268 = intHL0s25s26s34s56x1220D4eps1()
      t269 = t46 + t260
      t270 = s25 + s26 + s34 + s56
      t271 = t46 + t262
      t272 = s13 * t269 - t270 * t271
      t273 = s12 + s13
      t274 = t16 * s16
      t275 = s25 - s12
      t276 = s16 - s34 - s13
      t277 = s25 - s12 + s16 - s34
      t278 = t16 * (s26 - s13 + s45)
      t279 = s25 + s26
      t280 = s13 - s45
      t281 = s34 ** 2.0_dp
      t282 = (t16 * t276 + s45 + t275 + t41) * s45
      t283 = t273 * s13
      t284 = t273 * s25
      t285 = t16 * (s16 * t273 - (s12 + s13 + s16 + s26) * s34 - t247) -
     # (-t18 * t273 - s26 - t274) * s26 - (t277 + t278 + s46) * s46 + (t
     #10 * t280 + t18 * (s12 - s46) - t279) * s56 + t281 + s16 ** 2.0_dp
     # - t282 - t283 + t284
      t286 = s12 + s13 + s16
      t287 = t16 * s13
      t288 = s12 - s26
      t280 = (-s34 + t286) * s25 + (s25 - s12 - s13) * s26 + (-t287 - s2
     #5 - s12 + s26 + s45) * s45 + (-t16 * t280 + s46 - t288) * s46 + (-
     #t288 - t103 + s56) * s56 + t283
      t289 = s25 - s13
      t290 = s26 - s13
      t291 = t16 * t290
      t292 = -s12 + s16 - s13
      t293 = s16 * s25
      t294 = t289 * s34
      t277 = t16 * t247 - t292 * s13 + (-t287 + t275) * s26 - (s12 - s16
     # + s34 - t291 - s45) * s45 - (-t226 - t278 - s46) * s46 + (-t18 * 
     #t96 + t277 + t59) * s56 + t293 - t294
      t295 = -s13 + s25 + s26 + s56
      t262 = t295 + t262
      t296 = intHs160s26s34s56x1012D2eps0()
      t297 = -t234
      t298 = intHs160s26s34s56x1021D2eps0()
      t299 = intHs160s26s34s56x1013D4eps0()
      t300 = intHLs160000x0122D4eps0()
      t301 = intHs160s26s34s56x1011D2eps1()
      t302 = intHL0s25s260s56x1011D2eps1()
      t303 = intHL0s25s26s34s56x1110D2eps1()
      t304 = intHL0s25s260s56x1010D0eps1()
      t305 = intHL0s25s260s56x1020D2eps1()
      t306 = intHL0s25s26s34s56x1120D2eps0()
      t307 = intHL0s25s26s34s56x1210D2eps0()
      t308 = intHL0s25s26s34s56x1310D4eps0()
      t309 = t18 * s16
      t310 = t40 * s13
      t311 = 7.0_dp
      t312 = t10 * s25
      t313 = -s13 * t311 + s26 + t312
      t314 = s25 - s26
      t315 = s26 + s45 + s46
      t316 = t18 * s13
      t317 = t10 * t315
      t318 = t311 * s34
      t319 = s25 ** 2.0_dp
      t320 = s12 * s13
      t321 = t18 * t247
      t322 = t16 * (t319 + t320) - (t309 + s13) * s25 - (-t310 + t62 + s
     #26) * s26 + t313 * s34 + (-s16 * t40 + s25 + t316 - t317 + t318) *
     # s56 + t248 * (-s16 * t39 + t10 * t314 + 8.0_dp * s34) - t321 + t3
     #10 * s16
      t323 = (-s12 - s16 - s25 - s26) * s26
      t324 = s13 * s16
      t325 = t286 * s25
      t326 = t60 * s34
      t290 = -t72 - t290
      t327 = (t274 + t56 + t310) * s25
      t237 = (-t237 + s25) * s46
      t328 = s13 ** 2.0_dp
      t329 = (t287 - s25 + s12 - s45 + s46 - s56) * s56
      t330 = (s12 + s13 + s16 + s25 + s26) * s26
      t331 = t279 * s34
      t332 = intHLs16s25s26s34s56x1321D6eps0()
      t333 = s25 + s16 - s13
      t334 = t16 * t275
      t335 = t10 * s34
      t336 = -t16 * (t248 * (t279 + t260) + t320) + (t51 + s16) * s25 - 
     #(-t51 - t333 - s26) * s26 - (t334 - t316 - s16 + s34) * s34 - (t33
     #5 + s25 - s16 - s13 - t228 + s56) * s56 - t324
      t337 = -s13 + s16
      t338 = t10 * s13
      t339 = s26 + s16
      t340 = s12 + s16
      t341 = s13 * t340
      t342 = t337 * s25
      t343 = (t59 - t338 + t21) * s34
      t315 = t16 * (-t248 * (t339 - t258) + t341) + (-t62 + t316 - s26) 
     #* s26 + (-t16 * t315 - s56 - t333 + t335) * s56 - t342 + t343
      t333 = intHLs16s25s26s34s56x1131D4eps0()
      t344 = t16 * t273
      t345 = s12 - s45
      t346 = t16 * t345
      t347 = (t51 + s25 + s16 + s13 + s26) * s26
      t348 = (t344 + s16) * s25
      t349 = t16 * t279
      t350 = intHLs16s25s26s34s56x1221D4eps0()
      t351 = -t72 - t41 + s13
      t352 = s25 + s26 + s34
      t353 = -s25 + s16 - s34 + s13 - s46 + t346
      t354 = t352 * s46
      t355 = t16 * (-s45 * t352 + t320)
      t356 = s26 - s34
      t357 = (-s12 + s13 - s16 - s25 - s26) * s26
      t291 = (-s25 - t291) * s34 - t16 * ((-t16 * t339 - s12 + s13 - s25
     # + s34 - s45 - s46 - s56) * s56 - t248 * (t274 + t356) - t293 + t3
     #20 + t357) - t338 * s16
      t293 = s13 + s16
      t273 = t10 * t273
      t358 = t16 * (s26 + s16 - s13 + s45 + s46)
      t359 = t248 * (t339 - t260)
      t281 = -t16 * (-t359 + t324 - t247 + t281) + (t293 + s26) * s25 + 
     #(t274 + t273 + t41) * s34 + (s25 + t358 - t258) * s56
      t360 = -t287 + s25
      t361 = s16 + s26 - s34
      t362 = s25 - s34
      t363 = t248 * t361
      t364 = t360 * s26
      t358 = -t16 * (-t363 + t324 - t247) + s25 * t293 + (t344 + t339 - 
     #s34) * s34 + (t362 + t358) * s56 + t364
      t365 = -t72 + t316 - s26
      t366 = s25 - s26 - s16 + s34
      t367 = t248 * t366
      t368 = (t18 * t289 + s16 + s26 + t51) * s26
      t276 = (t16 * (s12 + s45 + s46) + t18 * t276 + s25 + t36) * s56
      t369 = t16 * (t367 - t325 + t320)
      t370 = t18 * (t324 - t247)
      t371 = intHLs16s25s26s34s56x1231D6eps0()
      t372 = t51 + t293
      t373 = s25 + s45 + s46
      t288 = t16 * t288
      t374 = t10 * t373
      t375 = (-t51 - t337 - s26) * s26
      t376 = (-t372 + s25) * s25
      t377 = t16 * (t248 * (t349 + s34) + t320) + t321 + (t16 * t289 + s
     #26) * s34 + (t18 * (s34 - s13) - s16 - t288 + t374) * s56 + t324 +
     # t375 + t376
      t378 = s26 + s34
      t379 = t10 * t248
      t380 = t40 * s34
      t381 = t311 * s25
      t289 = t10 * (t72 + t378) * t248 + t16 * (s16 * t96 + t319 + t320)
     # + t321 - (t273 + t309) * s25 + (t10 * t289 + s26) * s34 + (-t18 *
     # t293 - t288 + t379 + t380 + t381) * s56 - t347
      t293 = t18 * t245
      t382 = t373 * s16
      t383 = s25 * t248
      t384 = s25 * s26
      t385 = s26 + s45 + s56
      t232 = t232 * s25
      t386 = (-t385 + s13) * s13
      t345 = t345 * s25
      t387 = -t16 * (-s46 * t385 - t244 - t345 + t386) + t232 + t21 * s1
     #3
      t388 = s26 * s46
      t389 = s46 * s56
      t232 = -t16 * (-s45 * s46 - t244 - t345 + t386) + t18 * (s13 * s25
     # + t389) + t232 + t388
      t345 = intHs160000x0211D2eps0()
      t386 = intHL0s25s260s56x1031D4eps0()
      t390 = s26 - s56
      t391 = s46 * t390
      t392 = -s16 + s25
      t393 = -t392 - t258
      t394 = s26 * s34
      t395 = t10 * s34 * t248 + t16 * (s34 + s45 + s46) * s25 + s13 * t3
     #93 - s16 * t271 + s25 * t46 + t394 + t258 * s56
      t396 = intHLs16s25s26s34s56x1211D4eps1()
      t397 = s13 - s56
      t398 = t16 * t245
      t399 = t18 * t397
      t400 = t40 * s46
      t401 = t39 * s45
      t314 = -t56 + t400 + t401 - t398 + t314 - t399
      t402 = t10 * t397
      t403 = t311 * s46
      t404 = 8.0_dp * s45
      t405 = t39 * t248
      t406 = t405 + t334 + s26 - s16 + s34 - t399
      t407 = 8.0_dp * t248
      t408 = intHLs16s25s26s34s56x1121D4eps1()
      t409 = s46 + s56
      t401 = t40 * t409 + s26 + t21 - t338 + t401
      t334 = s56 * t40 + t334 - t338 - t361 + t405
      t410 = t311 * s56
      t411 = s13 - s45 - s56
      t400 = -t39 * t411 + t21 + t400
      t412 = t310 * s25
      t385 = t16 * (-s13 * t251 + s26 * s45 + s56 * t385) + t18 * (t389 
     #+ t319) + s25 * (t410 + t405 + s26) + t388 - t412
      t388 = -s16 + s25 + s34
      t389 = t16 * t340
      t405 = s26 ** 2.0_dp
      t413 = t248 * t388
      t414 = -s34 * t365
      t344 = (t344 + t309 - s25) * s25
      t357 = t16 * (-(t398 - s25 + s12 + s26) * s56 + t320 + t357)
      t415 = t10 * (t413 + t324)
      t416 = t18 * s46
      t417 = (-s13 + s25 + s26) * s34
      t286 = t16 * (s13 * t286 + (s25 - s16 + s34 - s13) * s45 - t244 + 
     #t323 - t325 + t417) - (t16 * (s16 + s45) + t356) * s46 - (t16 * (s
     #12 + s26 + s16 - s34 + s13) + t416) * s56
      t418 = t338 - t21
      t419 = s26 + s16 - s34 - s13 + s45 + s46
      t420 = t51 * s13
      t421 = intHLs16s25s26s34s56x1221D4eps1()
      t352 = t10 * t352 * t248 + t321 + t420 + epinv * (-s34 * t351 + (-
     #t353 + s56) * s56 - t347 - t348 + t354 - t355) + (-t338 + t21 + s2
     #6) * s34 + (t380 + t374 - t316 - s16 - t288) * s56 + t324 + t375 +
     # t376
      t374 = t10 * (t248 * (t72 + s26) + t247 + t319) - t412 + epinv * (
     #-t16 * t290 * s45 - t16 * (-t331 + t329 + t330 + t244 - t328) + t2
     #37 - t327) + (-t287 + t21) * s26 + (9.0_dp * s25 + t10 * (s26 - s1
     #3 + s45 + s46)) * s56
      t375 = t16 * t419 + s25
      t363 = -t363 - t247 + t341
      t376 = t72 - t310 + s16
      t412 = t18 * t337
      t337 = t72 * t337
      t422 = t16 * t360
      t423 = t10 * t383
      t253 = t253 + s25
      t378 = t16 * (s34 * t246 + t319) + t423 + s13 * (-t21 + t398) - s1
     #6 * t253 + s25 * t378 + t21 * s56
      t424 = intHLs16s25s26s34s56x1141D6eps0()
      t425 = s13 * (t97 + s25)
      t426 = t46 * t253
      t427 = t425 - t426
      t428 = s12 * t46
      t429 = intHLs16s25s26s34s56x1131D4eps1()
      t430 = t16 * t32
      t431 = t430 + s25
      t432 = -t72 + t32
      t433 = s13 * (-t21 + t430) - t253 * t432
      t434 = s25 + s56
      t435 = t16 * t251
      t436 = t253 * (t21 + t435) - t287 * (t16 * t434 + s26)
      t437 = epinv * t262
      t438 = epinv * t267
      t439 = t59 + t226
      t440 = propW16 * t227
      t441 = t440 * t46 * t32
      t442 = s12 + s16 + s26 - s34
      t443 = intHs16s25s26s34s56x1211D2eps0()
      t444 = (t10 * t32 + t51) * s13 - t253 * (t430 - s25)
      t445 = intHLs160000x0121D2eps1()
      t446 = t440 * s16
      t447 = t16 * (t446 - s25)
      t448 = t10 * s45
      t449 = t446 * t39
      t450 = t20 * (t29 * (-t18 * t409 - s26 + t316 - t448 - t72) + t449
     #) + t22 * (t18 * (s13 - s46 - s56) - s26 + t447 - t448)
      t451 = intHLs160000x0121D2eps0()
      t447 = t10 * t411 - t416 + t447
      t338 = t29 * (t338 - t10 * (s45 + s56) - t416 - t72) + t449
      t452 = t18 * epinv
      t453 = intHL0s25s260s56x1012D2eps0()
      t454 = intHL0s25s260s56x1013D4eps0()
      t455 = intHL0s25s260s56x1021D2eps0()
      t456 = intHs160000x0111D0eps0()
      t457 = intHL0s25s26s34s56x1130D4eps0()
      t458 = -t16 * t411 + s25 + s46
      t459 = intHLs160000x0211D2eps0()
      t460 = t16 * (s12 + s13 - s56)
      t416 = t448 - t460 + t388 + t416
      t448 = t379 - t460 + t388
      t460 = intHs160000x0121D2eps0()
      t366 = -t51 + t379 + t366 - t399
      t461 = intHLs16s25s26s34s56x1211D4eps0()
      t462 = intHLs16s25s26s34s56x1121D4eps0()
      t260 = -t260 - t392
      t372 = t16 * t260 * t248 + (-t51 - s13 - s16) * s16 + s25 * t372 +
     # s26 * t392 + (-t16 * (s25 - s12 - s16) + s26 + t316 - s34) * s34 
     #+ s56 * t393
      t294 = t16 * (t324 + t294) - (s26 + t389 + t316 - s25) * s25 + (t3
     #12 - t398) * s45 - (-t312 + t398 + s26) * s46 + (t21 - t398 + s46)
     # * s56
      t339 = t10 * t341 - t321 - t337 - (t376 + s26) * s26 + (-t310 + t2
     #1 + s26) * s34 + (-t72 + t380 - t317 - t412) * s56 + t248 * (s34 *
     # t39 - t10 * t339)
      t380 = intHs16s25s26s34s56x1111D2eps0()
      t393 = intHL0s25s260s56x1021D2eps1()
      t463 = s26 + s56 + t84
      t464 = epinv * s56
      t465 = t25 * s26
      t466 = -s34 * t25 + s16 + s26 + s56
      t467 = t93 * s34
      t468 = t32 + t467
      t469 = intHs160000x0211D2eps1()
      t470 = -t10 + t452
      t471 = intHL0s25s260s56x1031D4eps1()
      t472 = intHs16s25s26s34s56x1112D2eps0()
      t473 = intHLs16s25s26s34s56x1311D4eps1()
      t474 = intHL0s25s26s34s56x1220D4eps0()
      t475 = intHL0s25s26s34s56x1120D2eps1()
      t476 = epinv * s34
      t477 = -t476 + s25 + s26 + s56
      t478 = -t16 * (s12 + s26 + s13 - s45 - s46 - s56) - t245
      t479 = intHLs16s25s26s34s56x1321D6eps1()
      t480 = intHLs16s25s26s34s56x1231D6eps1()
      t481 = intHs160000x0111D0eps1()
      t482 = intHL0s25s26s34s56x1130D4eps1()
      t483 = intHs16s25s26s34s56x1113D4eps0()
      t484 = intHLs160000x0211D2eps1()
      t485 = intHs160000x0121D2eps1()
      t486 = intHL0s25s26s34s56x1210D2eps1()
      t487 = intHL0s25s26s34s56x1310D4eps1()
      t488 = intHs160s26s34s56x1012D2eps1()
      t489 = intHs160s26s34s56x1013D4eps1()
      t490 = intHL0s25s260s56x1022D4eps0()
      t491 = 0.1e1_dp / t46
      t492 = 0.1e1_dp / ecossin ** 2.0_dp
      t493 = t152 * t24
      t358 = t493 * (-t126 * (t16 * t375 * s56 + epinv * t358 + s26 * t4
     #22 - t422 * s34 - t10 * t363 + t337) - t267 * t345 - t358 * t9)
      t422 = t492 * propW34
      t494 = t239 * t491
      t495 = t136 * t27
      t238 = t422 * (t152 * (t24 * (-t82 * t249 * t156 - t438 * t469) - 
     #t490 * t30 * t256 * t491 + t156 * (t277 * t81 + t280 * t31) * t27)
     # + t25 * (t30 * (t152 * (-t262 * t482 - t267 * t471 + t294 * t74 -
     # t336 * t479 + t377 * t480 + t385 * t71 + t436 * t58) - t494 * t22
     #2) + t493 * t87 * t249 * t156 + (t124 * t229 + t315 * t78 + t339 *
     # t76 + t378 * t65 + t395 * t67 + t433 * t53 + t444 * t55) * t152 *
     # t27) + t152 * (-t30 * (t236 * t238 + t473 * (epinv * t238 - s16 +
     # s25 - s26 + s34 + t379 - t399 - t51)) - t495 * (epinv * t478 + s2
     #6 + t379 - t399 + t72)) * s34)
      t249 = -t25 * t487 + t308
      t496 = t25 * t489
      t497 = t483 * s34
      t498 = s34 * t30
      t277 = t422 * (t27 * (-t14 * t267 - t15 * (-t51 + t379 + t438 + s2
     #5 - s26 - s16 + s34 - t399) + t152 * (t156 * (-t25 * (t147 * t280 
     #+ t277 * t86) + (t496 - t299) * t297 * s34) - t497 * t478)) + t498
     # * t152 * t249 * t491 * t222 - t183 * t61 * t30 * t372)
      t232 = t277 + t238 + t422 * (t30 * (t152 * (-t133 * t436 + t134 * 
     #t372 + t139 * t232 - t143 * t385 - t144 * t294 - t150 * (-t18 * s2
     #5 * t397 - epinv * t232 + t16 * t319 + t384 - t391 + t423) + t262 
     #* t457 + t267 * t386 + t332 * t336 + t333 * (-s46 * t279 + (-s25 +
     # s16 + s13 - s46 + t346 - s56) * s56 - t326 + t347 + t348 - t349 *
     # s45) - t371 * t377 - t429 * (-epinv * (-t16 * (s45 * t46 + t394 -
     # t428) + (s16 - s34 - s46) * s25 + (s25 + s16 - s46) * s26 + (-s25
     # + s16 - s46) * s56 - t247 + t254 + t405) - t16 * t427) + t250 * t
     #25 * t256 * t491) + t454 * t222 * t491) - (t131 * t433 + t132 * t4
     #44 + t135 * t378 + t138 * t395 + t145 * t339 + t146 * t315 + t151 
     #* t229) * t152 * t27 + t358)
      t238 = t100 * (t415 + t357 + epinv * (t16 * (t413 - t405 + t320) -
     # (t389 + s13) * s25 - (t16 * (s12 - s13) + s16 + t21) * s26 + (s25
     # - t293 - t63) * s56 + t414 + t316 * s16) + t343 - t344) + t108 * 
     #t322 + t114 * (epinv * t334 + s16 + s26 - s34 - t21 + t310 - t407 
     #- t410 + t51) + t116 * (epinv * t406 + s16 - s34 - t21 + t288 + t4
     #02 - t407) + t125 * t291 + t128 * (t10 * (t367 + t324 - t247) + ep
     #inv * (t369 - t368 + t414 - t276 + t370) - (t51 + t309 + s13 - s25
     #) * s25 + s26 * t418 - t418 * s34 + (-t10 * t419 + s25) * s56 + t4
     #20)
      t254 = -t111 * (epinv * t291 - t360 * s34 + s56 * t375 - t16 * t36
     #3 + t342 + t364) + t122 * t281 + t123 * (t16 * (s12 * t295 + s26 *
     # t251 - s34 * t373 + t382 - t383) + t18 * (s56 * t245 + t384) + s1
     #3 * (-t59 + s25 - t293) - s25 * s56 + s26 * t245) + t13 * (s34 * t
     #365 + t276 + t368 - t369 - t370) + t91 * (-t10 * (-t359 + t341) + 
     #t321 + epinv * t281 + (t376 + s26) * s26 - t313 * s34 + (t72 + t31
     #7 - t318 + t412) * s56 + t337)
      t256 = t422 * t152
      t276 = t25 * t33
      t277 = t85 * t156
      t281 = t25 * t156
      t288 = t112 * s34
      t275 = t3 * t448 + t366 * t460 + t380 * (s13 + s26 - s56) + t288 *
     # t265 * t267 - t83 * (-t321 + t287 * s16 + epinv * t280 - t292 * s
     #25 + (t316 - t275) * s26 + t360 * s34 - (t16 * (s25 + s16 - s34 - 
     #s13 + s45) + s12 + t41 + s46) * s46 + (t10 * t96 - t16 * t226 + s1
     #2 - t41) * s56 - t282 - t283) * t156
      t282 = epinv * t445
      t283 = t30 * t491
      t255 = t422 * (t152 * (t24 * t275 + t27 * (t262 * t456 + t334 * t7
     # + t406 * t6 - t281 * (t110 * t264 + t277 * t285)) + t30 * (t107 *
     # t286 + t118 * t458 - t2 * (t51 + s25 + s26 + s16 - s34 - s13 - s4
     #6 + s56) + t28 * t400 + t314 * t461 + t401 * t462 + t45 * (-epinv 
     #* (-t16 * (s13 * (s25 - s16 + s56) - s25 * s45 + (s25 + s16 - s34 
     #- s45) * s26 + (-s25 + s16 - s45) * s56 - t247 + t382 + t405 + t42
     #8) + s25 * s34 + s46 * t46) + t425 - t426) - t459 * t416 - t276 * 
     #t289 + t393 * (-t16 * (-s13 * t463 + t255) - s25 * (-t16 * (-s45 *
     # t25 - s46 * t25) + t104 + t465 - t464 + t84)) * t491) + t282 * t4
     #50) + t283 * (t453 + t455) * t222)
      t275 = t88 + t89
      t291 = t156 ** 2.0_dp
      t292 = t25 * intHL0s25s260s56x1012D2eps1()
      t294 = s34 * t152
      t295 = t115 * t470
      t313 = t472 * s34
      t263 = t422 * (t152 * (t24 * (t206 * t448 + t313 * t267) + t27 * (
     #-t200 * t322 + t295 * t267 + t437 * t481)) + t30 * (t491 * (t152 *
     # (t222 * (-t221 + t240 - t241 - t242) + t272 * t474 + t475 * (s13 
     #* (t46 - t467) - t271 * t477)) - t292 * t222) - t294 * (t101 * t45
     #8 + t38 * (epinv * t458 + s25 - t266)) + t152 * (t119 * t458 - t41
     #6 * t484) * epinv) + t191 * t27 * t285 * t291 + t152 * (t24 * (t14
     #8 * t235 + t231 * t234 - t275 * t280) + t27 * (t264 * t5 + t4 * (s
     #13 * t468 + t263 * t466 + t51 * t466))) * t156)
      t264 = t306 + t307
      t271 = t25 * t488
      t285 = t485 * epinv
      t315 = s34 * t486
      t62 = t152 * (t156 * epinv2 * (t20 * (-t16 * propW34 * (t11 * (-s5
     #6 * t442 + t323 + t324 - t325 + t413 + t417) + t16 * (-t248 * (t72
     # + t41 + t245) + t324) - t321 + (t273 + s16) * s25 + (s13 * t39 + 
     #s26 + t56 + t62) * s26 - (t287 + t279) * s34 + (t39 * t96 - t21 + 
     #t245 + t56 - t59) * s56) + t441 * t39 * propW34) - t153 * propW34 
     #* s25 * (s13 * t245 + (t287 + s12) * s26 + (-s25 + s12 - s26 + t10
     #3 - s56) * s56 - t248 * t439 + t284 - t441)) + t30 * propW34 * (t2
     #5 * (t222 * (-t315 - t304 + t305) - t268 * t272) + s34 * t222 * t2
     #64 - t222 * (t302 + t303) * t17)) * t491
      t62 = t152 * (t451 * propW34 * t450 + t1 * propW34 * (t20 * t338 +
     # t22 * t447) + t120 * (t20 * propW34 * (epinv * t338 - t449) + t22
     #4 * propW34 * (t452 * t447 - t449))) - t155 * propW34 * t27 * (t43
     #7 - s26 - s16 + s34 - s13 + s56 - t228) + t62 + propW34 * (t152 * 
     #(t156 * (t17 * (t234 * t301 + t235 * t80) + (-t271 + t296 + t298) 
     #* t297 * s34 + t190 * t280) + t285 * t366) + t262 * t443) * t24
      t43 = t492 * t62 + t263 + t255 + t256 * (t24 * t254 + t27 * t238 +
     # t30 * (t102 * (epinv * t286 + t343 - t344 + t357 + t415) + t109 *
     # (-t16 * ((s25 + s26 - s16) * s45 + t247 + t323 + t324 - t325) + (
     #t274 - t279) * s46 + (-t16 * (s25 - s12 - s16 - s13 + s45) - s46) 
     #* s56 - t326) + t121 * (epinv * t400 - t312 - t403 + 8.0_dp * t411
     #) - t141 * (t16 * (s45 * t290 + t244 - t328 + t329 + t330 - t331) 
     #- t237 + t327) - t47 * (-epinv * t387 + s25 * t262) + t94 * t387 +
     # t95 * t374 + t99 * t289 - t350 * (s34 * t351 + (t353 - s56) * s56
     # + t347 + t348 - t354 + t355) + t396 * (epinv * t314 - s25 + t293 
     #+ t402 - t403 - t404 + t43 + t59) + t408 * (epinv * t401 - t311 * 
     #t409 - s26 + t310 - t312 - t404) + t421 * t352))
      t62 = t25 * intHLs160000x0122D4eps1()
      t80 = intHLs160000x0111D2eps1() * t17
      t96 = (-t62 + t300) * s16 + t223 + t80
      t153 = t25 * intHLs16s25s26s34s56x1141D6eps1()
      t222 = t25 * intHLs16s25s26s34s56x1411D6eps1()
      t228 = t25 * intHLs160000x0113D4eps1() - intHLs160000x0113D4eps0()
      t234 = -s16 + s25 - s26 + s34 - s56
      t235 = epinv * s25
      t237 = s25 + s34
      t238 = epinv * t32
      t244 = t28 * t30
      t247 = t17 * t301
      t248 = t81 * t234
      t254 = t252 * t490
      t255 = t393 * t463
      t263 = epinv * t119
      t272 = t83 * (t235 + s26 + s16 - s34 + s56)
      t273 = (t296 + t298) * s34
      t280 = t493 * (s34 * t5 + t148 + t247 + t248)
      t284 = t172 * t20
      t286 = t47 * epinv
      t287 = epinv * t45
      t289 = t121 * t470
      t290 = t100 * t112
      t23 = t152 * (t24 * (-t122 * t260 + t388 * (t290 + t123) + t91 * (
     #-epinv * t260 + s34)) + t27 * ((t112 * t265 + t472) * s34 + t111 *
     # epinv * t392 - t200 * t261) + t30 * ((t112 * t38 - t107 + t286 + 
     #t94) * s34 - t102 * (t288 + s25 - s16) - t141 * t237 + t392 * (t28
     #7 + t109) - t118 + t282 - t289) + t23 * t380) + t156 * (t27 * (t15
     #2 * (t273 + t272 + t149) - t190) + t280 + epinv2 * (t20 * ((propZ2
     #5 * t388 - t10) * s25 - t251 * t39 - t398) - t22 * (t226 + t435)) 
     #* t152 * t491) + t284 * t152 + t152 * ((t491 * (-t306 - t307 + t47
     #4) + t101) * s34 + t25 * (t270 * t33 + t491 * (t304 - t305)) + t49
     #1 * (t17 * (t302 + t303) + t254 + t255) - t263) * t30 + t292 * t28
     #3 + t27 * (t88 + t89) * t156
      t110 = s34 * t110
      t226 = t86 * t234
      t251 = s34 * t488
      t250 = t250 * t252
      t260 = t197 * t156
      t288 = t12 * t460
      t293 = t99 * t270
      t296 = t31 * t156
      t297 = t25 * t87
      t298 = epinv * t469
      t234 = t152 * (t20 * (-t12 * t82 - t29 * t386) + t27 * (-t132 * t3
     #2 + t297 - t298)) + t24 * (t152 * ((-t138 - t483) * s34 - t135 * t
     #388 + t234 * (-t9 - t131)) - t14) + t30 * (t152 * (t133 * t46 + t1
     #43 * t252 + t144 * t388 + t150 * (t476 + s25)) + t25 * (t152 * t47
     #1 + t494) + t333 + t294 * (t473 + t134 + t139)) + (t152 * (-t126 *
     # t24 * t234 - t120 * t30) + t30 * t429) * epinv + t24 * (-t112 * t
     #15 - t296 + t294 * (-t112 * t136 - t156 * t299)) + t25 * (t24 * (t
     #147 * t156 + t152 * ((t156 * t489 + t67) * s34 + t234 * t53 + t388
     # * t65)) + t55 * t27 * t32 * t152 + t152 * ((t487 * t491 - t61) * 
     #s34 - t252 * t71 - t388 * t74 - t46 * t58) * t30) + t283 * (-t294 
     #* t308 - t454) - t152 * (t12 * t345 + t193) * t20
      t283 = t50 - t216
      t301 = t283 * t24
      t306 = t220 * (t1 - t345 + t386 - t82) + t301 + t219
      t307 = 9.0_dp * t220
      t23 = -t16 * t234 - t18 * (t152 * (t26 * t8 + t244) + t220 * (-t2 
     #- t7 - t380 + t460 - t451)) + t39 * t306 + t152 * (t24 * (-t128 * 
     #(epinv * t432 + s16 + s26 - s34 + s56) - t13 * t432 + t156 * t231)
     # + t30 * (t32 * t350 - t421 * (s34 - t238) + t491 * (t221 + t242) 
     #- t95 * (epinv * t237 + s25 + s26 + s56) - t293) + t108 * t27 * t2
     #61) + t152 * (t24 * (t295 + t260) + t125 * t392 * t27) + t30 * (-t
     #152 * t462 + t491 * (t152 * (-t240 + t241) - t453 - t455)) + t23 +
     # t183 * (-t114 * t24 + t156 * (t24 * (-t226 - t110) + t27 * (-t251
     # - t85)) + t30 * (t491 * ((t486 - t268) * s34 - t250) + t408)) + t
     #152 * (t27 * (s34 * t4 * t156 + t485) - t498 * t475 * t491) * epin
     #v + (-t7 * t12 + t29 * t451 + t288) * t152 * t20 + t307 * t8
      t234 = s25 + s26 - s34 + s56
      t237 = t39 * propW16
      t261 = t274 * propW16
      t295 = 1.0_dp + t261
      t306 = t237 * s16
      t308 = t16 + t306 + t11
      t272 = -t273 - t272 - t149
      t115 = t114 + t115
      t273 = t7 + t8
      t310 = epinv * t38
      t314 = s34 * t265
      t316 = t277 * t25
      t317 = t1 * t308
      t320 = epinv2 * (t20 * (t237 * t32 - t10 + t11) + t22 * (t430 * pr
     #opW16 - 1.0_dp)) * t156
      t308 = t120 * ((t20 * t308 + t22 * t295) * epinv - t306 * t225)
      t321 = t26 * t345
      t62 = ((t282 - t300 + t62 + t451) * s16 - t223 - t80) * t152 * t22
     #5 * propW16 + t220 * (t174 - t82) + t219
      t80 = t496 - t299
      t223 = t25 * t71
      t299 = -t223 + t143
      t300 = t25 * t74
      t322 = -t300 + t144
      t323 = -t25 * t61 + t134
      t324 = t25 * t58
      t325 = t24 * t156
      t252 = t152 * (t30 * (t323 * s34 - t139 * t362 + t150 * (t476 + t8
     #4) + t252 * t299 + t388 * t322 + t46 * (-t324 + t133) + t28) + s34
     # * (t156 * t80 - t136) * t24 + t297 * t27) + t26 * (-t152 * t82 + 
     #t13 + t14 + t9) - (t15 + t126 + t128) * t24 * epinv + t325 * (-t31
     # + t180)
      t295 = t1 * t295
      t326 = propW16 * s16
      t327 = t326 * t225
      t328 = -t327 * t152 * t228
      t62 = -t16 * t252 + t18 * t220 * (t295 + t2 + t7 + t8 - t345 - t46
     #0) + t39 * t62 - 12.0_dp * t328 - t152 * (t24 * (t156 * (t25 * (t2
     #26 + t110) - t231) - t290 * t46) + t27 * ((-epinv * t4 + t271) * t
     #156 * s34 - t298) - t308 - t321 + t293 * t30) - t152 * (t20 * (t12
     # * t273 - t288 - t317) + t24 * (t25 * t115 + t46 * (-t122 - t123) 
     #- t91 * (epinv * t46 + s34) - t260) + t27 * (-(t111 * t234 + t485)
     # * epinv - t314 + t316) + t30 * ((-t310 - t94) * s34 + t102 * (epi
     #nv * t257 + s25 + s26 + s56) + t107 * t257 + t95 * (t476 + s26 + s
     #56 + t84) - t287 * t46) - t320) + t152 * (t125 * t27 * t234 + t284
     #) - t156 * (t27 * (t152 * t272 + t190 - t88 - t89) - t280) - t152 
     #* ((-t286 - t101) * s34 - t141 * t362 + t118 - t186 + t189 + t263 
     #- t276 * t270) * t30
      t101 = t16 + epinv
      t110 = s12 + s26 + s16 - s34 + s56
      t186 = t16 * t110
      t234 = t10 * t340
      t252 = t356 * t40 + t104 + t234 + t72
      t110 = t10 * t110
      t257 = t110 + s25
      t271 = -s13 + s25 + s26 + s45 + s56
      t280 = t311 * t356 + t340 * t39 + s25 + t104
      t284 = t16 * (s26 - s13 + s45 + s56) + t21
      t258 = t16 * (s12 + s26 + s16) - t258
      t286 = -t16 * t442 + s25
      t288 = -t186 + s25
      t290 = t389 - t335 + t41 + t434
      t293 = t326 * t10
      t298 = 1.0_dp + t293
      t311 = -t16 * (s16 - s34 + s13 - s45 - s56) + s25 - t56
      t278 = t243 + t21 + t278
      t329 = t40 * s26
      t326 = 12.0_dp * t326
      t330 = t16 + t326 + t11
      t331 = t396 + t408
      t334 = t270 * t350
      t336 = epinv * t484
      t115 = t115 * t470
      t337 = t270 * t421
      t199 = t154 * t6 + t27 * (t152 * (t156 * (t259 * (t5 - t205) + t4 
     #* t468) - t200 * t280) + t155) + t152 * (t20 * (t12 * t3 + t29 * (
     #t2 + t459) + t330 * t451) + t24 * (t100 * (s25 * t112 - t10 * (s12
     # + s16 - s34) - t16 * t465) + t111 * (t16 * (t112 * t340 + s26 - s
     #34 + s56) + s25) + (t122 + t123) * t60 + t91 * (t16 * (epinv * s26
     # - s25) - t104 - t234 + t235 + t318 - t329) + t199 + t115) + t27 *
     # (t108 * t280 + t206) + t30 * (-t107 * t311 - t141 * t278 + t25 * 
     #t331 + t284 * t94 - t187 - t289 - t334 + t336 - t461) + t282 * (t2
     #0 * t330 + t22 * t298)) + t152 * ((t284 * t47 - t337) * epinv - t1
     #02 * (epinv * t311 + s25 + t59) - t95 * (epinv * t278 + s25 + t59)
     # + t204 - t462) * t30
      t204 = t26 * t184
      t78 = t25 * t78
      t234 = t156 * t27
      t259 = t15 * t101
      t278 = (-t175 - t64) * t25
      t280 = t26 * t152
      t226 = (-t27 * t469 + t498 * t473) * epinv + t25 * (t24 * (t226 * 
     #t156 + t286 * t65) + t27 * (t251 * t156 + t252 * t76 + t257 * t55 
     #+ t258 * t67)) - t308 - t321 - t325 * t231
      t110 = t152 * t226 + t152 * (t27 * (-t132 * t257 - t138 * t258 - t
     #145 * t252 - t146 * t290 + t25 * (t277 + t124) - t285) - t320) + t
     #24 * (-t152 * (t125 * t340 + t126 * (-t110 + t92) + t128 * (t235 -
     # t186) + t131 * t288 + t135 * t286 + t260) - t259) + t152 * (t156 
     #* (t24 * (-t247 - t248 - t148) + t27 * t272) + t27 * (t78 * t290 -
     # t314) + t30 * (t46 * (epinv * t429 - t109) + t118 + t168 + t169 -
     # t185 + t263) + t53 * t25 * t24 * t288 - t317 * t20) + t182 + t234
     # * (t190 - t88 - t89) + t152 * ((-t310 + t236) * s34 + t150 * (t93
     # * t271 + s25) + t46 * (t333 - t287) + t278) * t30 - t280 * (t151 
     #+ t460)
      t80 = -t10 * (t152 * (t27 * (t297 - t82) + t30 * (t139 * t271 + t2
     #79 * t299 + t322 * t68)) + t325 * (t294 * t80 + t180 - t31)) - t16
     # * t110 - t18 * (t152 * (-t204 + t244) + t220 * (-t298 * t451 - t2
     # + t3 - t459)) - 24.0_dp * t328 + t39 * (t220 * (t174 + t295 - t34
     #5 - t151 - t460) + t152 * t214 * (t186 + s25) * t27 + t493 * (t136
     # + t137 - t215) * s34 + t301 + t219) + t199 - 9.0_dp * t220 * t184
     # - 12.0_dp * t152 * t225 * propW16 * t96
      t110 = s25 - s26 + s56
      t136 = s16 + s56
      t97 = t97 + s25 - s13 + s45 + s46
      t186 = t104 + t279 - t103
      t199 = t51 + t312 + t41 + s56
      t214 = t112 * s12
      t226 = t16 * (t214 + s25 + s26)
      t49 = epinv * t49
      t231 = t30 * t152
      t67 = t25 * t67 - t138
      t138 = t25 * t480
      t244 = t473 * s34
      t247 = t152 * t20
      t108 = -t200 + t108
      t200 = epinv * t481
      t84 = t27 * (t114 * t117 + t116 * t470 + t128 * (t21 + t113 + t238
     #) + t13 * (t51 + t32) + t200) + t30 * (t45 * (-s26 * t101 + t464 -
     # t84) - t47 * (-t49 + s25) - t95 * (-t16 * (-epinv * t97 + s26) - 
     #t21)) + t285 * t24
      t33 = t152 * t84 + t209 + t154 * t3 + t24 * (t208 * t152 + t443) +
     # t24 * (t112 * t155 + t152 * (t122 * t60 + t195 + t196 - t205 + t2
     #07)) + t231 * ((-t484 - t119 + t120) * epinv - t118 + t185 + t188 
     #- t28 - t459 - t474 - t475) + t247 * (-(t4 + t5 + t380) * t12 + t1
     #93) + t234 * (t191 - t190 + t88 + t89) + t231 * (t102 * (epinv * (
     #-t104 + t356 + t103) + t226) + t107 * (-t104 + t356 + t233) + t109
     # * t110 + t117 * t396 + t17 * t408 + t199 * t99 + t25 * (-t199 * t
     #33 + t121 + t268) - t350 * t52 - t421 * (-t51 * t25 - s16 - s56 - 
     #t21 + t235 - t329) + t476 * t38) + t152 * (t100 * (epinv * t439 + 
     #t226) + t123 * t227 + t108 * (t51 + s26 - s16 + s34 - s56 + t21) -
     # t316) * t27 + t247 * (t12 * t456 - t12 * t460 + t29 * t462)
      t38 = t152 * (t27 * (t217 + t218) + t30 * (s34 * (-t222 + t230) + 
     #t212 * t390 + t60 * (-t153 + t424))) + t220 * (s25 * t9 + t2 + t46
     #1 - t7)
      t84 = t280 * t6 + t220 * (t1 + t4 + t5 + t380 - t456 + t460 + t462
     #)
      t34 = t152 * (t30 * t279 * (-t324 + t133) + (-t25 * t55 + t132) * 
     #t34 * t27)
      t15 = t10 * t34 - t16 * (t27 * (t388 * t152 * t67 - t131 - t135 - 
     #t15) + t30 * (t152 * (t110 * (t25 * (t71 + t74) - t143 - t144) + t
     #139 * t186 + t141 * t97 - t150 * (-epinv * t186 + s25 - s26 + s56)
     # - t323 * t388 + (t138 - t371) * (t51 + t329 + t21 + t136) + (t25 
     #* t479 - t332) * (t41 + s16 + s34 + s56 + t35) - t244) + t333) + t
     #20 * (t12 * (-t152 * t7 + t9) - t152 * (t2 + t461) * t29) + t24 * 
     #(-t152 * (t166 + t206 + t167) - t296) + t25 * (t24 * (t124 * t152 
     #+ t156 * (t152 * (t165 + t164) + t147)) + t27 * (t152 * (t179 + t1
     #78) + t53 + t65) - t231 * t482) - (t159 + t160 + t163) * t152 * t2
     #7 + t231 * (t429 * (-t36 + t92) + t457) - t493 * t173 * t156) + t1
     #8 * t84 + t39 * t38 + t33 - t307 * t6
      t29 = s12 - s13 + s25 + s45 + s46 + s56
      t33 = s12 + s13 + s16 + s26 - s34 - s45 - s46
      t34 = -s13 + s25 + s45 + s46 + s56
      t36 = -s12 - s13 + s25 + s45 + s46 + s56
      t38 = t20 * (t16 - t306 + t11) + t22 * (1.0_dp - t261)
      t53 = t20 * (-t16 + t326 - t11) + t22 * (-1.0_dp + t293)
      t41 = t16 * (s25 - s12 - s16 + s34) + s56 - t41
      t55 = t10 * t356 + t136 * t16 + s25
      t84 = t16 * t46
      t92 = t312 + t104 - s26
      t42 = t59 + t42 + t243
      t51 = s25 * t39 + t10 * (s34 + s56) - t309 - t51
      t97 = t274 + t46
      t101 = 1.0_dp + t93
      t103 = -t16 * t361 + s25
      t63 = t335 - t309 + s25 - t63
      t104 = t124 * t25
      t9 = t256 * (t24 * ((t469 + t485) * epinv + t125 * (s12 - s16 - s2
     #6 + s34 - s56) + t126 * (epinv * t431 - t56 - t72) + t13 * t32 + t
     #431 * t9 + t313 - t104) + t27 * (t128 * (t35 - t238) + t157 + t158
     # + t161 + t162 + t200) + t30 * (t117 * t331 + t134 * t63 + t170 + 
     #t171 - t336) + epinv2 * (t20 * (t237 * t46 + t11 + t16) + t22 * (t
     #84 * propW16 + 1.0_dp)) * t491 - t282 * t38)
      t11 = t256 * (t30 * (t133 * t42 + t143 * t92 + t144 * t41 + t150 *
     # (epinv * t140 - s26 + s56 + t21) + t25 * t491 * (t268 * t270 + t2
     #50 - t315)) - t281 * (t147 * t36 + t29 * t87 + t34 * t86) * t24 + 
     #t159 * t27)
      t13 = t422 * (t27 * (t294 * (-t215 + t137) + t50 - t216 + t152 * (
     #-t213 + t130) * t52) + t231 * (-t211 + t142) * t110)
      t22 = t256 * (t24 * (t111 * (t16 * (t214 - t238) + s25) + t122 * t
     #55 + t91 * (epinv * t55 - s26 - s34 + s56 - t56 - t72) - t205) + t
     #27 * (t100 * (epinv * t103 + t21 + t56 + t59) + t103 * t123 + t192
     # - t281 * (t277 * t33 + t36 * t90)) + t30 * (t102 * (-t10 * t245 +
     # t106 + t21 - t59) + t45 * (epinv * t97 - t84) + t47 * (-t72 + t49
     #) + t51 * t99 + t95 * (t16 * (-epinv * t98 + s26) + t243 + t381) +
     # t105 + t185 + t188 + t189 - t263)) + t422 * (t152 * (t1 * t53 + t
     #120 * (epinv * t53 - 12.0_dp * t327) + t30 * (t109 * t97 + t203 * 
     #t17 - t276 * t51 - t118)) + t24 * (t101 * t155 + t152 * (t156 * (t
     #83 * (t25 * t34 + t214) + t148 + t197) + t4 + t5)) + t152 * (t116 
     #* (epinv * t10 - t40) + t156 * (t149 * t33 * t156 + t275 * t36) + 
     #t115 - t201) * t27)
      t33 = t256 * (t30 * ((t222 - t230) * s34 + t46 * (t153 - t424)) + 
     #t96 * t225 * propW16)
      t35 = t422 * (t30 * (t152 * (t25 * (t269 * t479 - t471 - t482) - t
     #269 * t332 - t337 + t386 + t457 - t461 - t462) + t491 * (-t239 * t
     #25 + t294 * t249 + t454) + t152 * (t138 - t371) * (s34 + t84)) - t
     #152 * (t497 + t6) * t27)
      t40 = t256 * (t24 * (t206 + t3) + t27 * t273)
      t45 = t429 * t46 + t244
      t49 = 0.4e1_dp / 0.3e1_dp
      t9 = t10 * t13 + t16 * t40 + t49 * (t11 + t9 + t422 * (t152 * (t27
     # * (t79 * (-t78 + t146) + t456) - t38 * t451) + t30 * (t152 * (t2 
     #+ t28 - t459) + t491 * (t152 * (-t17 * t302 - t475 * t477 - t221 +
     # t240 - t241 - t242) - t292 + t453 + t455)) + t259 * t27 + t493 * 
     #(t156 * (t31 * t36 + t34 * t81) + t151 + t345 + t460)) + t422 * (t
     #24 * (t152 * (t82 * t29 * t156 + t476 * t265) + t443) + t27 * (t15
     #2 * (-t25 * (t178 + t54 + t57 + t176 + t177) + t163 * t101) + t14)
     # + t231 * (t25 * (-t41 * t74 - t42 * t58 + t491 * (-t304 + t305) -
     # t61 * t63 - t71 * t92) - t491 * (t17 * t303 + t270 * t474 + t254 
     #+ t255) + t491 * t264 * s34))) + 16.0_dp * t327 * t256 * t228 - 0.
     #8e1_dp / 0.3e1_dp * t35 + 0.16e2_dp / 0.3e1_dp * t256 * t30 * t45 
     #+ 0.2e1_dp / 0.3e1_dp * t22 - 8.0_dp * t33
      t11 = epinv * t60
      t13 = t266 + t56 + s25
      t14 = t266 + s25
      t22 = -t6 + t151
      t11 = t24 * (s25 * t155 - t100 * (t59 + t56 + t21 - t11) - t198 * 
     #t17 - t205 + t206 + t207 + t3 + t4 + t5) + t30 * (-t102 * (-epinv 
     #* t13 + s25 + t59) - t95 * (-epinv * t14 + s25 + t59) + t188 + t2 
     #- t28 - t451 - t459 + t461 + t462) + t24 * (-t44 * t108 - t116 * t
     #117 - t91 * (t72 + t56 - t11 + s26 + s34 - s56) + t194) + t30 * ((
     #t202 - t484 + t337 - t445) * epinv + t107 * t13 - t121 * t17 + t14
     # * t141 + t25 * (t37 - t396 - t408) - t187 + t334) - t204
      t11 = struc9MP * (t10 * t30 * (-t139 * t48 + t279 * (t223 - t143) 
     #+ t68 * (t300 - t144)) + t16 * (t24 * (t128 * t129 - t104 - t161 -
     # t162 + t167 + t174) + t24 * (t126 * t127 + t25 * (t57 + t54) + t6
     #6 * (t25 * t65 - t135) + t67 * t69 + t77 * (t25 * t76 - t145) + t7
     #9 * (t78 - t146) + t181 * s25) + t30 * (epinv * t45 + s34 * t236 -
     # t150 * (t93 * t48 + s25) + t333 * t46 - t168 - t169 - t278) + t49
     #5 * s34 - t26 * t22) + t18 * t224 * t184 + t39 * (t24 * (-s25 * t2
     #83 - t217 - t218) + t219 * s25 + t224 * t22) + t11) - struc60MP * 
     #t27 * (t25 * t85 + t32 * (-t16 * (-t25 * t86 - t180 - t297 + t31 +
     # t81 + t82) - t25 * t83 - t190 + t88 + t89) - t149) * t291
      t13 = t422 * (t23 * struc2MP + t62 * (struc3MP + struc8MP) + t231 
     #* (-t16 * (epinv * t150 + t139) + (t47 - t95 - t102) * epinv - t10
     #7 - t141 + t94) * struc61MP)
      t9 = (t10 * t256 * t446 * t225 * t228 - t16 * t422 * (t152 * (t27 
     #* ((s13 * t431 - t253 * t32 + t420) * (t213 - t130) + (t215 - t137
     #) * t262 * s34) + t30 * ((t211 - t142) * (-t16 * s25 * (s13 - t246
     #) - t391 + t319) + (-t153 + t424) * t427 + (-t222 + t230) * t229 *
     # s34)) + t440 * t152 * t96 * t225 + t27 * (t216 - t50) * t267) - 0
     #.2e1_dp / 0.3e1_dp * t232 + t43 / 3.0_dp + t256 * t8 * t27 * t267)
     # * struc1MP + t9 * struc7MP
      result = t49 * t13 + 0.2e1_dp / 0.3e1_dp * t422 * (t152 * t11 + (-
     #t16 * (t24 * (t156 * (-t152 * t173 + t180 - t31) - t181) + t152 * 
     #(t20 * (-t12 * t6 - t172) + t24 * (-t128 * t129 + t25 * (t156 * (t
     #165 + t164) + t124) - t166 - t167) + t27 * (-t163 - t157 - t158 - 
     #t159 - t160 - t161 - t162) + t30 * (t143 * t73 + t144 * t75 + t168
     # + t169 - t170 - t171)) + t182 + t183 * (t179 + t57 + t176 + t177 
     #+ t178 + t54) * t27 + t152 * (t150 * (-epinv * t140 + s25 + s26 - 
     #s56) - t25 * (t71 * t73 + t74 * t75 + t175 + t64) + t28) * t30) + 
     #t18 * t220 * (t1 + t3 + t4 + t5 - t7 - t8) + t39 * (t27 * (t217 * 
     #t152 + t218 * t152 - t216 + t50) + t219 * t152 * t70 + (t2 - t6 + 
     #t174) * propZ25 * t19) - t210) * struc10MP + t80 * struc5MP + t15 
     #* struc6MP) + t9

           ampNonresonantLightFullReC4MP = result
       end function ampNonresonantLightFullReC4MP

       function ampNonresonantLightFullReC4PM()
           implicit none
           complex(dp) :: ampNonresonantLightFullReC4PM
           complex(dp) ::  t1,t10,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109
           complex(dp) ::  t11,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t12
           complex(dp) ::  t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t13,t130
           complex(dp) ::  t131,t132,t133,t134,t135,t136,t137,t138,t139,t14,t140,t141
           complex(dp) ::  t142,t143,t144,t145,t146,t147,t148,t149,t15,t150,t151,t152
           complex(dp) ::  t153,t154,t155,t156,t157,t158,t159,t16,t160,t161,t162,t163
           complex(dp) ::  t164,t165,t166,t167,t168,t169,t17,t170,t171,t172,t173,t174
           complex(dp) ::  t175,t176,t177,t178,t179,t18,t180,t181,t182,t183,t184,t185
           complex(dp) ::  t186,t187,t188,t189,t19,t190,t191,t192,t193,t194,t195,t196
           complex(dp) ::  t197,t198,t199,t2,t20,t200,t201,t202,t203,t204,t205,t206
           complex(dp) ::  t207,t208,t209,t21,t210,t211,t212,t213,t214,t215,t216,t217
           complex(dp) ::  t218,t219,t22,t220,t221,t222,t223,t224,t225,t226,t227,t228
           complex(dp) ::  t229,t23,t230,t231,t232,t233,t234,t235,t236,t237,t238,t239
           complex(dp) ::  t24,t240,t241,t242,t243,t244,t245,t246,t247,t248,t249,t25
           complex(dp) ::  t250,t251,t252,t253,t254,t255,t256,t257,t258,t259,t26,t260
           complex(dp) ::  t261,t262,t263,t264,t265,t266,t267,t268,t269,t27,t270,t271
           complex(dp) ::  t272,t273,t274,t275,t276,t277,t278,t279,t28,t280,t281,t282
           complex(dp) ::  t283,t284,t285,t286,t287,t288,t289,t29,t290,t291,t292,t293
           complex(dp) ::  t294,t295,t296,t297,t298,t299,t3,t30,t300,t301,t302,t303
           complex(dp) ::  t304,t305,t306,t307,t308,t309,t31,t310,t311,t312,t313,t314
           complex(dp) ::  t315,t316,t317,t318,t319,t32,t320,t321,t322,t323,t324,t325
           complex(dp) ::  t326,t327,t328,t329,t33,t330,t331,t332,t333,t334,t335,t336
           complex(dp) ::  t337,t338,t339,t34,t340,t341,t342,t343,t344,t345,t346,t347
           complex(dp) ::  t348,t349,t35,t350,t351,t352,t353,t354,t355,t356,t357,t358
           complex(dp) ::  t359,t36,t360,t361,t362,t363,t364,t365,t366,t367,t368,t369
           complex(dp) ::  t37,t370,t371,t372,t373,t374,t375,t376,t377,t378,t379,t38
           complex(dp) ::  t380,t381,t382,t383,t384,t385,t386,t387,t388,t389,t39,t390
           complex(dp) ::  t391,t392,t393,t394,t395,t396,t397,t398,t399,t4,t40,t400
           complex(dp) ::  t401,t402,t403,t404,t405,t406,t407,t408,t409,t41,t410,t411
           complex(dp) ::  t412,t413,t414,t415,t416,t417,t418,t419,t42,t420,t421,t422
           complex(dp) ::  t423,t424,t425,t426,t427,t428,t429,t43,t430,t431,t432,t433
           complex(dp) ::  t434,t435,t436,t437,t438,t439,t44,t440,t441,t442,t443,t444
           complex(dp) ::  t445,t446,t447,t448,t449,t45,t450,t451,t452,t453,t454,t455
           complex(dp) ::  t456,t457,t458,t459,t46,t460,t461,t462,t463,t464,t465,t466
           complex(dp) ::  t467,t468,t469,t47,t470,t471,t472,t473,t474,t475,t476,t477
           complex(dp) ::  t478,t479,t48,t480,t481,t482,t483,t484,t485,t486,t487,t488
           complex(dp) ::  t489,t49,t490,t491,t492,t493,t494,t5,t50,t51,t52,t53
           complex(dp) ::  t54,t55,t56,t57,t58,t59,t6,t60,t61,t62,t63,t64
           complex(dp) ::  t65,t66,t67,t68,t69,t7,t70,t71,t72,t73,t74,t75
           complex(dp) ::  t76,t77,t78,t79,t8,t80,t81,t82,t83,t84,t85,t86
           complex(dp) ::  t87,t88,t89,t9,t90,t91,t92,t93,t94,t95,t96,t97
           complex(dp) ::  t98,t99

           complex(dp) :: result

      t1 = intHLs160000x0112D2eps0()
      t2 = intHLs16s25s26s34s56x1111D2eps0()
      t3 = intHs160000x0112D2eps0()
      t4 = intHs160s26s34s56x1021D2eps1()
      t5 = intHs160s26s34s56x1022D4eps0()
      t6 = intHs16s25s26s34s56x1112D4eps0()
      t7 = intHs16s25s26s34s56x1121D4eps0()
      t8 = intHs16s25s26s34s56x1211D4eps0()
      t9 = intHs16s25s26s34s56x1131D4eps0()
      t10 = 4.0_dp
      t11 = propZ25 * s25
      t12 = -t10 + t11
      t13 = intHs16s25s26s34s56x1221D4eps0()
      t14 = intHs16s25s26s34s56x1311D4eps0()
      t15 = intHs16s25s26s34s56x1311D4eps1()
      t16 = 2.0_dp
      t17 = -t16 + epinv
      t18 = 3.0_dp
      t19 = gw ** 2.0_dp
      t20 = gb ** 2.0_dp
      t21 = t18 * s25
      t22 = t21 * t19 * propZ25
      t23 = -t20 * t12
      t24 = t23 + t22
      t25 = 1.0_dp - epinv
      t26 = t20 * t12
      t27 = t26 - t22
      t28 = intHLs16s25s26s34s56x1112D4eps0()
      t29 = t11 + t16
      t30 = t20 * t29
      t31 = t30 + t22
      t32 = intHs16s25s26s34s56x1310D4eps0()
      t33 = s16 + s26 - s34 + s56
      t34 = intHLs16s25s26s34s56x1222D6eps1()
      t35 = s12 + s25
      t36 = t10 * s26
      t37 = t16 * t35
      t38 = s16 + t37 + t36
      t39 = intHLs16s25s26s34s56x1211D2eps1()
      t40 = intHs16s25s26s34s56x1222D6eps1()
      t41 = 6.0_dp
      t42 = 5.0_dp
      t43 = t41 * s12
      t44 = t42 * s25
      t45 = t18 * s26
      t46 = s34 - s56 + t43 + t45 + t44
      t47 = intHLs16s25s26s34s56x1121D2eps1()
      t48 = s25 + s26 + s56
      t49 = intHLs16s25s26s34s56x1112D2eps1()
      t50 = -s13 + s26 + s45 + s46 + s56
      t51 = t16 * t50 + s25
      t52 = intHs16s25s26s34s56x1411D6eps0()
      t53 = intHs160s26s34s56x1031D4eps1()
      t54 = intHs16s25s26s34s56x1121D2eps0()
      t55 = intHs16s25s26s34s56x1131D4eps1()
      t56 = t10 * s12
      t57 = t17 * s25
      t58 = t57 - t56
      t59 = intHs16s25s26s34s56x1221D4eps1()
      t60 = t16 * s12
      t61 = t57 - t60
      t62 = intHs16s25s26s34s56x1141D6eps0()
      t63 = t60 + s25
      t64 = intHs16s25s26s34s56x1321D6eps0()
      t65 = t21 + t60
      t66 = intHs16s25s26s34s56x1231D6eps0()
      t67 = t21 + t56
      t68 = intHLs16s25s26s34s56x1132D6eps0()
      t69 = t16 * s26
      t70 = t69 + s25
      t71 = intHLs16s25s26s34s56x1312D6eps0()
      t72 = s16 + s25
      t73 = t16 * (s26 + s12)
      t74 = t72 + t73
      t75 = intHs16s25s26s34s56x1312D6eps0()
      t76 = t21 + t73
      t77 = intHs16s25s26s34s56x1113D4eps1()
      t78 = intHs16s25s26s34s56x1114D6eps0()
      t79 = intHs16s25s26s34s56x1213D6eps0()
      t80 = s12 + s25 + s26
      t81 = t16 * t80 + s34
      t82 = intHLs16s25s26s34s56x1113D4eps0()
      t83 = s13 - s45 - s46
      t84 = t16 * t83
      t85 = t18 * s56
      t86 = -s26 - t85 + t84
      t87 = intHLs16s25s26s34s56x1122D4eps0()
      t88 = t16 * s56
      t89 = t83 - t88
      t90 = intHLs16s25s26s34s56x1114D6eps0()
      t91 = s25 + s26 - s56
      t92 = intHLs16s25s26s34s56x1123D6eps0()
      t93 = t16 * s25
      t94 = -s56 + t45 + t93
      t95 = intHLs16s25s26s34s56x1213D6eps0()
      t96 = -s56 + t37 + t45
      t97 = intHs16s25s26s34s56x1132D6eps0()
      t98 = s26 + s34 - s56 + t93 + t56
      t99 = intHs16s25s26s34s56x1123D6eps0()
      t100 = t16 * (s34 + s12) + t91
      t101 = intHs16s25s26s34s56x1310D4eps1()
      t102 = intHs16s25s26s34s56x1110D2eps0()
      t103 = intHs160s26s34s56x1020D2eps0()
      t104 = intHLs16s25s26s34s56x1113D4eps1()
      t105 = intHs16s25s26s34s56x1220D4eps1()
      t106 = intHs16s25s26s34s56x1130D4eps1()
      t107 = intHs160s26s34s56x1031D4eps0()
      t108 = 0.1e1_dp / s25
      t109 = t41 * t19 * propZ25
      t110 = t20 * (-propZ25 * t16 + 8.0_dp * t108) + t109
      t111 = intHs16s25s26s34s56x1211D2eps1()
      t112 = intHLs16s25s26s34s56x1132D6eps1()
      t113 = intHLs16s25s26s34s56x1312D6eps1()
      t114 = intHs16s25s26s34s56x1312D6eps1()
      t115 = intHs16s25s26s34s56x1213D6eps1()
      t116 = intHLs16s25s26s34s56x1123D6eps1()
      t117 = intHLs16s25s26s34s56x1213D6eps1()
      t118 = intHs16s25s26s34s56x1132D6eps1()
      t119 = intHs16s25s26s34s56x1123D6eps1()
      t120 = intHs16s25s26s34s56x1110D2eps1()
      t121 = intHs16s25s26s34s56x1220D4eps0()
      t122 = intHs16s25s26s34s56x1130D4eps0()
      t123 = intHs16s25s26s34s56x1120D2eps1()
      t124 = t25 * s25
      t125 = intHs160s26s34s56x1020D2eps1()
      t126 = intHs16s25s26s34s56x1321D6eps1()
      t127 = intHs16s25s26s34s56x1231D6eps1()
      t128 = intHs160s26s34s56x1022D4eps1()
      t129 = intHs160000x0112D2eps1()
      t130 = intHs16s25s26s34s56x1121D2eps1()
      t131 = 1.0_dp + epinv
      t132 = t60 * t131
      t133 = intHs16s25s26s34s56x1121D4eps1()
      t134 = intHs16s25s26s34s56x1211D4eps1()
      t135 = intHs16s25s26s34s56x1112D4eps1()
      t136 = t16 * epinv
      t137 = t136 - t18
      t138 = intHLs160000x0111D0eps0()
      t139 = intHLs160000x0111D0eps1()
      t140 = intHLs160000x0112D2eps1()
      t141 = intHLs16s25s26s34s56x1112D4eps1()
      t142 = intHs16s25s26s34s56x1122D4eps0()
      t143 = intHs16s25s26s34s56x1212D4eps0()
      t144 = intHLs16s25s26s34s56x1222D6eps0()
      t145 = intHs16s25s26s34s56x1212D4eps1()
      t146 = intHLs16s25s26s34s56x1211D2eps0()
      t147 = intHLs16s25s26s34s56x1212D4eps1()
      t148 = s26 - s34 + s25
      t149 = t148 - t85 + t84 + t56
      t150 = intHLs16s25s26s34s56x1212D4eps0()
      t151 = intHs16s25s26s34s56x1122D4eps1()
      t152 = intHLs16s25s26s34s56x1112D2eps0()
      t153 = intHLs16s25s26s34s56x1122D4eps1()
      t154 = intHs16s25s26s34s56x1222D6eps0()
      t155 = intHLs16s25s26s34s56x1121D2eps0()
      t156 = intHs16s25s26s34s56x1120D2eps0()
      t157 = intHs16s25s26s34s56x1210D2eps0()
      t158 = intHs16s25s26s34s56x1210D2eps1()
      t159 = 0.1e1_dp / t33
      t160 = s34 * t77
      t161 = t64 * t65
      t162 = t66 * t67
      t163 = t99 * t100
      t164 = t75 * t76
      t165 = t79 * t81
      t166 = t97 * t98
      t167 = t68 * t70
      t168 = t71 * t74
      t169 = t25 * (t159 * (t105 * t35 + t106 * t63) + t53)
      t170 = t55 * t58
      t171 = s12 * t54
      t172 = t2 * t29
      t173 = t25 * t101
      t174 = t9 + t13 + t14
      t175 = t70 * t112
      t74 = t74 * t113
      t176 = t126 * t65
      t177 = t127 * t67
      t178 = t114 * t76
      t81 = t115 * t81
      t98 = t118 * t98
      t100 = t119 * t100
      t179 = t159 * (-t108 * (t121 * t35 + t122 * t63) - t32 + t173)
      t180 = t15 * t17
      t181 = t26 * t174
      t182 = t108 * t25
      t183 = t155 * t48
      t184 = t51 * t152
      t185 = s34 * t146
      t186 = t137 * t141
      t187 = t38 * t144
      t188 = t135 * t137
      t189 = t154 * t46
      t190 = t158 * t25
      t191 = t190 * t159
      t192 = (t142 + t143) * t70
      t193 = t102 * t159
      t194 = t151 * (epinv * t69 - s26 - s34 + s56 - t56 + t57)
      t195 = -t125 * t182 + t156 + t157
      t196 = t17 * t120
      t197 = t123 * (t124 + t60) + t196
      t198 = t133 + t134
      t199 = t7 + t8
      t200 = s34 * t39
      t201 = t47 * t48
      t202 = t25 * t34
      t203 = t202 * t38
      t204 = t1 * t29
      t205 = epinv * t129
      t206 = t25 * t128
      t207 = t130 * (t132 + s25)
      t208 = t182 * t40
      t23 = t107 * t110 + t108 * (t20 * (t12 * t199 + t204) + t31 * ((t4
     #9 * t51 + t200 + t201) * epinv - t138 + t203)) + t24 * (t108 * (t1
     #59 * t197 + t205 - t206 + t207) + t111) + t27 * (t108 * t17 * t198
     # + t159 * t195 - t208 * t46) + t27 * (t108 * (t103 * t159 + t188 +
     # t189) - t191) + t23 * (t3 + t4 + t5) * t108 + t108 * ((-t139 + t1
     #40) * epinv + t147 * (epinv * t149 - s25 - t69) + t149 * t150 - t1
     #53 * (t16 * (-epinv * t89 + s26) + s25) + t183 + t184 + t185 - t18
     #6 - t187) * t31 + t108 * (t145 * ((-t18 + epinv) * s25 - t56 - t69
     # * t25) + t192 + t193 + t194) * t24
      t149 = t25 * intHs16s25s26s34s56x1141D6eps1()
      t209 = -t149 + t62
      t210 = t25 * intHLs16s25s26s34s56x1114D6eps1()
      t211 = t210 - t90
      t212 = t25 * intHs16s25s26s34s56x1114D6eps1()
      t213 = intHs16s25s26s34s56x1411D6eps1() * t25
      t214 = t209 * t63
      t215 = (-t212 + t78) * s34
      t174 = t174 * s25
      t216 = t211 * t31
      t217 = t19 * propZ25
      t218 = intHL0s25s26s34s56x1110D2eps0()
      t219 = t48 - t84
      t220 = intHLs160000x0111D2eps0()
      t221 = t217 * s25
      t222 = t20 + t221
      t223 = s16 - s34 + s25
      t224 = t73 + t223
      t225 = s26 + s16 - s34 + s13 - s56
      t226 = t16 * (s45 + s46 - s12)
      t227 = -t226 + t225
      t228 = intHLs16s25s26s34s56x1411D6eps0()
      t229 = intHs160s26s34s56x1011D2eps0()
      t230 = s26 + s16 - s34 - s56
      t231 = t16 * (s13 - s45 - s46 + s12)
      t232 = t231 + t230
      t233 = t16 * (s13 - s45 + s12) + s16 + s26 - s34 - s46 - s56
      t234 = intHLs16s25s26s34s56x1311D4eps0()
      t235 = s45 * t16 - s13 + s25 + s26 + s46 + s56
      t236 = intHL0s25s260s56x1010D0eps0()
      t237 = intHL0s25s260s56x1011D2eps0()
      t238 = intHL0s25s260s56x1020D2eps0()
      t239 = s25 - s12
      t240 = t16 * s13
      t241 = s13 - s25
      t242 = t16 * (s26 - s13)
      t243 = -s13 + s26 + s45
      t244 = t16 * t243
      t245 = s16 - s34 + s25 - s12
      t246 = s12 + s13 - s16
      t247 = s16 * s25
      t248 = t241 * s34
      t249 = s56 ** 2.0_dp
      t250 = t16 * t249 + t246 * s13 + (t239 - t240) * s26 + (s16 - s34 
     #- s12 + t242 + s45) * s45 + (t244 + t223 + s46) * s46 + (-t18 * t8
     #3 + t245 + t69) * s56 + t247 + t248
      t251 = s16 + s26 + s56
      t252 = t18 * s34
      t253 = t251 - t252
      t254 = t16 * s34
      t255 = t251 - t254
      t256 = s45 + s46
      t257 = t16 * t256
      t230 = -t257 + t230
      t258 = s13 * t253 + t230 * t255 + t255 * t60
      t259 = intHs16s25s26s34s56x1112D2eps1()
      t260 = t16 * (s13 - s45 - s46 - s56)
      t261 = -t260 + s25
      t262 = s12 + s13
      t263 = t16 * s16
      t264 = s16 - s34 - s13
      t265 = s25 + s26
      t266 = s13 - s45
      t267 = s34 ** 2.0_dp
      t268 = t262 * s13
      t269 = (t16 * t264 + s45 + t239 + t45) * s45
      t270 = t262 * s25
      t245 = -t16 * (-s16 * t262 + (s12 + s13 + s16 + s26) * s34 + t249)
     # - (-t18 * t262 - s26 - t263) * s26 - (t245 + t244 + s46) * s46 + 
     #(t10 * t266 - t18 * (s46 - s12) - t265) * s56 + s16 ** 2.0_dp + t2
     #67 - t268 - t269 + t270
      t271 = s12 + s13 + s16
      t272 = s26 - s12
      t266 = t16 * t266
      t273 = (-s34 + t271) * s25 - (s12 + s13 - s25) * s26 + (s26 - s25 
     #- s12 - t240 + s45) * s45 + (t272 - t266 + s46) * s46 + (t272 - t8
     #4 + s56) * s56 + t268
      t274 = intHL0s25s260s56x1020D2eps1()
      t275 = intHL0s25s26s34s56x1120D2eps0()
      t276 = intHL0s25s26s34s56x1210D2eps0()
      t277 = intHL0s25s26s34s56x1310D4eps0()
      t278 = -s13 + s25 + s26 + s56
      t279 = t278 + t257
      t280 = intHs160s26s34s56x1012D2eps0()
      t281 = -t232
      t282 = intHs160s26s34s56x1021D2eps0()
      t283 = intHs160s26s34s56x1013D4eps0()
      t284 = intHs160s26s34s56x1011D2eps1()
      t285 = intHL0s25s260s56x1011D2eps1()
      t286 = intHL0s25s26s34s56x1110D2eps1()
      t287 = intHL0s25s260s56x1010D0eps1()
      t288 = intHLs160000x0122D4eps0()
      t289 = intHL0s25s26s34s56x1220D4eps1()
      t290 = t48 + t254
      t291 = s25 + s26 + s34 + s56
      t292 = t257 + t48
      t293 = t10 * s56
      t294 = s45 + s46 + s56
      t295 = s16 - s34
      t296 = s46 ** 2.0_dp
      t297 = t10 * s56 * t256 + t16 * (s45 * t295 + (s16 - s34 + s45 + s
     #25) * s46 + s56 * t223) + t18 * (s26 * t294 + t249) - s12 * (s26 -
     # s13 + s45 - s46 + s56 + s25) + (-t16 * (s16 - s34 + s45 + s46) - 
     #s25 - t45 - t293 + s13) * s13 + (s26 + s16 - s34 + s45) * s25 + t2
     #96 + s45 ** 2.0_dp
      t298 = intHL0s25s260s56x1022D4eps1()
      t299 = s26 + s56
      t300 = t299 + t93
      t301 = t16 * t294
      t302 = s13 * t300
      t303 = t299 * t294
      t304 = (t301 + s26 + s25) * s25 - t302 + t303
      t305 = propW16 * s16
      t306 = t305 * t224
      t307 = s13 - s45 - s56
      t308 = t16 * (t306 - s25)
      t309 = t18 * s46
      t310 = t10 * t307 + t308 - t309
      t311 = t10 * s13
      t312 = t306 * t41
      t313 = t29 * (-t10 * (s45 + s56) - t309 - t93 + t311) + t312
      t314 = t18 * epinv
      t315 = intHLs16s25s26s34s56x1131D4eps1()
      t316 = s26 ** 2.0_dp
      t317 = s26 * s34
      t318 = s12 * t48
      t301 = t301 + s25
      t319 = s13 * (t88 + s25)
      t320 = t48 * t301
      t321 = t319 - t320
      t322 = t16 * t33
      t323 = t322 + s25
      t324 = t60 * s13
      t325 = t33 - t93
      t326 = s13 * (-t21 + t322) - t301 * t325
      t327 = s56 + s25
      t328 = t16 * t299
      t329 = t301 * (t21 + t328) - t240 * (t16 * t327 + s26)
      t330 = intHLs160000x0121D2eps1()
      t331 = t10 * s45
      t332 = s46 + s56
      t333 = t18 * s13
      t334 = t20 * (t29 * (-t18 * t332 - s26 - t331 + t333 - t93) + t312
     #)
      t308 = t22 * (t18 * (s13 - s46 - s56) - s26 + t308 - t331)
      t335 = s16 - s25
      t336 = t335 - t252
      t337 = intHLs16s25s26s34s56x1211D4eps1()
      t338 = s26 - s25
      t339 = s13 - s56
      t340 = t41 * s45
      t341 = t42 * s46
      t342 = t16 * t295
      t343 = t18 * t339
      t344 = -t343 - t338 - t342 + t341 + t340 - t56
      t345 = 7.0_dp
      t346 = 8.0_dp * s45
      t347 = t345 * s46
      t348 = t18 * t295
      t349 = t10 * t339
      t350 = t41 * t256
      t351 = t16 * t239
      t352 = -t343 + s26 - s16 + s34 + t351 + t350
      t353 = 8.0_dp * t256
      t272 = t16 * t272
      t354 = intHLs16s25s26s34s56x1121D4eps1()
      t340 = t332 * t42 + s26 + t21 - t311 + t340
      t355 = t10 * s25
      t356 = t42 * s13
      t357 = s16 + s26 - s34
      t358 = s56 * t42 - t311 + t350 + t351 - t357
      t359 = t345 * s56
      t341 = -t307 * t41 + t21 + t341
      t360 = s25 ** 2.0_dp
      t361 = t356 * s25
      t350 = t16 * (-s13 * t299 + s26 * s45 + (s26 + s45 + s56) * s56) +
     # t18 * (s46 * s56 + t360) + s25 * (s26 + t350 + t359) + s26 * s46 
     #- t361
      t362 = t16 * t262
      t363 = t18 * s16
      t364 = s16 - s25 - s34
      t365 = s12 + s16
      t366 = t16 * t365
      t367 = s26 - t333 + t93
      t368 = t256 * t364
      t369 = s12 * s13
      t370 = t367 * s34
      t371 = s13 * s16
      t372 = (-s12 + s13 - s16 - s25 - s26) * s26
      t373 = (-t362 - t363 + s25) * s25
      t374 = (t21 + t69 - t311) * s34
      t375 = t10 * (t371 - t368)
      t376 = t16 * (-(s26 - s25 + s12 + t342) * s56 + t372 + t369)
      t377 = s13 + s16
      t378 = t377 + t60
      t379 = s25 * s45
      t380 = (t378 + s26 - s34) * s25
      t381 = -t16 * (-t296 + t379) + (-t266 + s26) * s46 - (-t309 + s25)
     # * s56 + t380
      t382 = s25 * s26
      t383 = s26 - s13 + s45 + s56
      t384 = -s45 + s12
      t385 = s26 - s34
      t386 = -t16 * (s16 - s13 + s45) - t385
      t387 = s12 + s16 + s26 - s34
      t388 = t16 * t387
      t389 = s13 * t365
      t390 = (-s13 + s25 + s26) * s34
      t391 = (t388 + t309) * s56
      t392 = t16 * (-s25 * t365 - s45 * t364 - t296 + t372 + t389 + t390
     #)
      t393 = -t21 + t311
      t394 = s16 - s25 + s26 - s34
      t395 = s26 + s16 - s34 - s13 + s45 + s46
      t396 = t256 * t394
      t271 = t271 * s25
      t397 = (-t18 * t241 + s16 + s26 + t60) * s26
      t264 = (t16 * (s45 + s46 + s12) + t18 * t264 + s25 + t36) * s56
      t398 = t16 * (t369 - t396 - t271)
      t399 = t18 * (t371 - t249)
      t400 = intHLs16s25s26s34s56x1221D4eps1()
      t401 = -s13 + s16
      t402 = s25 + s26 + s34
      t403 = s45 + s46 + s25
      t404 = t42 * s34
      t405 = t10 * t403
      t406 = -s13 + t45 + t93
      t407 = -t16 * t384
      t408 = -s16 + s34 - s13 + s46 + s25 + t407
      t409 = (s16 + s13 + s25 + t60 + s26) * s26
      t410 = t402 * s46
      t411 = (t362 + s16) * s25
      t412 = t16 * (-s45 * t402 + t369)
      t413 = (-t378 + s25) * s25
      t414 = (-t401 - t60 - s26) * s26
      t415 = t18 * t249
      t416 = s26 + t93
      t417 = (s12 + s16 + s25 + s26) * s26
      t266 = -t16 * (-s34 * t265 - s45 * t416 + (s12 + s13 - s25 - s45 +
     # s46 - s56) * s56 + t296 + t417) - (t333 + t263 + t56) * s25 + (t2
     #66 + s25) * s46
      t418 = -t240 + s25
      t419 = t16 * t395 + s25
      t420 = s26 + s16
      t242 = (-t242 - s25) * s34 + t16 * ((t16 * t420 + s12 - s13 + s25 
     #- s34 + s45 + s46 + s56) * s56 + t256 * (t385 + t263) + t247 - t36
     #9 - t372) - t311 * s16
      t247 = t256 * t357
      t372 = t389 - t247 - t249
      t421 = t418 * s26
      t422 = t401 * s25
      t423 = s16 + t93 - t356
      t424 = s13 * t345 - s26 - t355
      t425 = s26 + s45 + s46
      t426 = t345 * s34
      t427 = t18 * t401
      t428 = t10 * t425
      t262 = t10 * t262
      t429 = t16 * (s26 + s16 - s13 + s45 + s46)
      t430 = t256 * (t420 - t254)
      t267 = -t16 * (t267 - t430 + t371 - t249) + (t377 + s26) * s25 + (
     #t45 + t262 + t263) * s34 + (t429 - t252 + s25) * s56
      t401 = t93 * t401
      t431 = t16 * t418
      t432 = s34 - s25
      t247 = -t16 * (-t247 + t371 - t249) + s25 * t377 - (-t420 - t362 +
     # s34) * s34 + (-t432 + t429) * s56 + t421
      t362 = s25 * t256
      t429 = s26 + s34
      t433 = t10 * t362 + t16 * (s34 * t294 + t360) + s13 * (-t21 + t342
     #) - s16 * t301 + t429 * s25 + t21 * s56
      t434 = intHLs16s25s26s34s56x1141D6eps0()
      t435 = t403 * s16
      t248 = t16 * (-t248 + t371) - (s26 + t366 + t333 - s25) * s25 + (-
     #t342 + t355) * s45 - (t342 + s26 - t355) * s46 + (t21 - t342 + s46
     #) * s56
      t436 = t10 * t389 - t401 - t415 - (t423 + s26) * s26 + (t21 + s26 
     #- t356) * s34 + (-t428 - t427 - t93 + t404) * s56 + t256 * (s34 * 
     #t41 - t10 * t420)
      t338 = t16 * (t369 + t360) - t415 - (t363 + s13) * s25 - (t72 - t3
     #56 + s26) * s26 - s34 * t424 + (-s16 * t42 + s25 + t333 + t426 - t
     #428) * s56 + t256 * (-s16 * t41 - t10 * t338 + 8.0_dp * s34) + t35
     #6 * s16
      t437 = t70 * s34
      t438 = intHLs16s25s26s34s56x1321D6eps0()
      t439 = s16 - s13 + s25
      t440 = t10 * s34
      t351 = -t16 * (t256 * (t265 + t254) + t369) + (s16 + t60) * s25 - 
     #(-t439 - t60 - s26) * s26 - (t351 - s16 - t333 + s34) * s34 - (-s1
     #6 - s13 + s25 + t440 + t226 + s56) * s56 - t371
      t420 = t16 * (-t256 * (t420 - t252) + t389) - (t72 - t333 + s26) *
     # s26 - (t16 * t425 + s56 + t439 - t440) * s56 + t374 - t422
      t425 = intHLs16s25s26s34s56x1131D4eps0()
      t439 = t16 * t265
      t441 = intHLs16s25s26s34s56x1221D4eps0()
      t442 = intHLs16s25s26s34s56x1231D6eps0()
      t443 = t10 * t256
      t444 = t345 * s25
      t377 = t10 * (t429 + t93) * t256 - t16 * (-s16 * t83 - t360 - t369
     #) + t415 - (t363 + t262) * s25 + (-t10 * t241 + s26) * s34 + (-t18
     # * t377 + t272 + t404 + t443 + t444) * s56 - t409
      t429 = intHs160000x0211D2eps0()
      t445 = intHL0s25s260s56x1031D4eps0()
      t446 = s26 - s56
      t447 = intHL0s25s260s56x1022D4eps0()
      t448 = intHL0s25s260s56x1021D2eps1()
      t449 = epinv * s56
      t450 = t25 * s26
      t451 = intHLs16s25s26s34s56x1231D6eps1()
      t452 = t25 * s34
      t453 = -t452 + s26 + s16 + s56
      t454 = t136 * s34
      t455 = t454 + t33
      t456 = intHs160000x0211D2eps1()
      t457 = t314 - t10
      t458 = intHL0s25s260s56x1031D4eps1()
      t459 = intHs16s25s26s34s56x1112D2eps0()
      t460 = intHLs16s25s26s34s56x1311D4eps1()
      t461 = intHL0s25s26s34s56x1220D4eps0()
      t462 = intHL0s25s26s34s56x1120D2eps1()
      t463 = epinv * s34
      t464 = -s26 - s56 + t463 - s25
      t465 = -t16 * t307 + s25 + s46
      t466 = -t16 * (s26 + s13 - s45 - s46 - s56 + s12) - t295
      t467 = intHLs16s25s26s34s56x1321D6eps1()
      t468 = intHs160000x0111D0eps1()
      t469 = intHL0s25s26s34s56x1130D4eps1()
      t470 = intHs16s25s26s34s56x1113D4eps0()
      t471 = intHLs160000x0211D2eps1()
      t472 = t16 * (s13 - s56 + s12)
      t309 = -t364 + t331 - t472 + t309
      t472 = -t364 - t472 + t443
      t473 = intHs160000x0121D2eps1()
      t394 = -t394 - t343 + t443 - t60
      t254 = t335 - t254
      t378 = t16 * t254 * t256 - (s13 + t60 + s16) * s16 + s25 * t378 - 
     #s26 * t335 - (-t16 * (s16 - s25 + s12) - s26 - t333 + s34) * s34 +
     # s56 * t336
      t474 = intHLs160000x0121D2eps0()
      t475 = intHs16s25s26s34s56x1111D2eps0()
      t476 = intHL0s25s260s56x1012D2eps0()
      t477 = intHL0s25s260s56x1013D4eps0()
      t478 = intHL0s25s260s56x1021D2eps0()
      t479 = intHs160000x0111D0eps0()
      t480 = intHL0s25s26s34s56x1130D4eps0()
      t481 = intHLs160000x0211D2eps0()
      t482 = intHs160000x0121D2eps0()
      t483 = intHLs16s25s26s34s56x1211D4eps0()
      t484 = intHLs16s25s26s34s56x1121D4eps0()
      t485 = intHL0s25s26s34s56x1210D2eps1()
      t486 = intHs160s26s34s56x1012D2eps1()
      t487 = intHs16s25s26s34s56x1211D2eps0()
      t488 = epinv * t279
      t489 = t69 + t223
      t490 = propW16 * t224
      t491 = t490 * t48 * t33
      t492 = 0.1e1_dp / t48
      t493 = 0.1e1_dp / ecossin ** 2.0_dp
      t454 = (t25 * t289 - t461) * (s13 * t290 - t291 * t292) - t462 * (
     #s13 * (-t454 + t48) + t292 * t464)
      t494 = t146 * t465 + t39 * (epinv * t465 + s25 - t260)
      t43 = (-t139 * t465 + t309 * t471) * epinv + s34 * t494 - t141 * (
     #epinv * t341 + 8.0_dp * t307 - t347 - t355) - t147 * (epinv * (s46
     # * t386 - t391 + t392) + t373 + t374 + t376 + t375) - t153 * (t10 
     #* (t256 * t416 + t249 + t360) - t361 + epinv * t266 + (t21 - t240)
     # * s26 + (9.0_dp * s25 + t10 * (s26 - s13 + s45 + s46)) * s56) - t
     #337 * (epinv * t344 - s25 - t346 - t347 + t348 + t349 + t43 + t69)
     # + t49 * (-epinv * (t16 * (s25 * t384 + s46 * t383 + t296) + t225 
     #* s25) + s25 * t279) - t354 * (epinv * t340 - t332 * t345 - s26 - 
     #t346 - t355 + t356) - t400 * (t10 * t402 * t256 + t324 + epinv * (
     #s34 * t406 - (-t408 - s56) * s56 - t409 + t410 - t411 - t412) + (t
     #21 + s26 - t311) * s34 + (t405 + t272 - s16 - t333 + t404) * s56 +
     # t371 + t413 + t414 + t415) + t454 * t492
      t225 = epinv * t330
      t296 = t493 * t108 * propW34
      t307 = t159 ** 2.0_dp
      t332 = t25 * intHL0s25s260s56x1012D2eps1()
      t346 = s34 * t459
      t347 = t25 * t40
      t361 = t134 * t457
      t384 = t103 * t108
      t402 = t384 * t27 * t245 * t307
      t404 = t493 * propW34
      t230 = t404 * (t108 * (t24 * (-t261 * t346 - t3 * t472 - t394 * t4
     #82) + t27 * (-t261 * t361 + t338 * t347)) + t31 * (t108 * (-t138 *
     # t465 + t309 * t481 - t344 * t483 + t492 * (t219 * (t218 - t236 + 
     #t237 + t238) - t448 * (-t16 * (s13 * (-t124 - s26 - s56) + t303) +
     # s25 * (t16 * (-s45 * t25 - s46 * t25) - t85 - t450 + t449 - t124)
     #)) + t202 * t377) + t332 * t219 * t492) + t108 * (t24 * (-t102 * t
     #233 + t123 * (t240 * s16 + epinv * t273 + t246 * s25 + (-t239 + t3
     #33) * s26 + t418 * s34 - (t16 * (s16 - s34 - s13 + s45 + s25) + s1
     #2 + t45 + s46) * s46 + (t10 * t83 - t16 * t223 + s12 - t45) * s56 
     #- t268 - t269 - t415) - t229 * t232 + t273 * (t156 + t157)) + t27 
     #* (-t258 * t5 - t4 * (s13 * t455 + t230 * t453 + t453 * t60))) * t
     #159 - t402)
      t239 = t476 + t478
      t246 = -t13 * (-t367 * s34 + t264 + t397 - t398 - t399) - t142 * t
     #267 - t143 * (t16 * (s12 * t278 + s26 * t299 - s34 * t403 - t362 +
     # t435) + t18 * (s56 * t295 + t382) + s13 * (-t348 - t69 + s25) - s
     #25 * s56 + s26 * t295) - t151 * (-t10 * (t389 - t430) + t415 + epi
     #nv * t267 + (t423 + s26) * s26 + s34 * t424 + (-t426 + t428 + t427
     # + t93) * s56 + t401) - t475 * (s13 + s26 - s56)
      t243 = t24 * t246 + t27 * (-t154 * t338 - t242 * t54 - t279 * t479
     # - t352 * t6 - t358 * t7) + t31 * (-t144 * t377 + t150 * (-s46 * t
     #386 + t391 - t392) - t152 * (-t16 * ((-t243 - s46) * s46 + t379) -
     # (-s46 * t16 + s25) * s56 + t380) - t155 * (-t16 * ((s26 - s16 + s
     #25) * s45 + t249 - t271 + t371 - t417) + (-t265 + t263) * s46 + (t
     #16 * (s16 + s13 - s45 - s25 + s12) - s46) * s56 - t437) + t2 * (s2
     #6 + s16 - s34 - s13 - s46 + s56 + s25 + t60) - t28 * t341 - t340 *
     # t484 - t47 * (-epinv * (s25 * s34 + s46 * t48 + t16 * (s13 * (s16
     # - s25 - s56) + (-s16 + s34 + s45 - s25) * s26 + (-s16 + s45 + s25
     #) * s56 + t249 - t316 - t318 + t379 - t435)) + t319 - t320) - t87 
     #* t266 + t441 * (-s34 * t406 + (-t408 - s56) * s56 + t409 - t410 +
     # t411 + t412))
      t246 = propW34 * t31
      t243 = t108 * (propW34 * t243 - t474 * propW34 * (t308 + t334)) - 
     #t246 * t239 * t492 * t219
      t266 = t25 * t486
      t267 = t125 * t159
      t268 = t182 * t159
      t269 = s34 * t485
      t278 = s34 * (t275 + t276)
      t303 = (t285 + t286) * t17
      t120 = propW34 * (t108 * (s34 * (-t131 * t259 * t261 + t159 * t281
     # * (t266 - t280 - t282)) - t159 * (t120 * t233 + t232 * t284) * t1
     #7 - t191 * t273) - t279 * t487) * t24
      t72 = -t108 * (t1 * propW34 * (t20 * t313 + t22 * t310) + t140 * (
     #t20 * propW34 * (epinv * t313 - t312) + t221 * propW34 * (t310 * t
     #314 - t312)) + epinv2 * (t20 * (-t16 * propW34 * (t11 * (-s56 * t3
     #87 - t271 - t368 + t371 + t390 - t417) + t16 * (-t256 * (t295 + t4
     #5 + t93) + t371) - t415 + (s16 + t262) * s25 + (s13 * t41 + s26 + 
     #t56 + t72) * s26 - (t265 + t240) * s34 + (t41 * t83 - t21 + t295 +
     # t56 - t69) * s56) + t491 * t41 * propW34) - t109 * propW34 * s25 
     #* (s13 * t295 + (t240 + s12) * s26 + (-s26 - s25 + s12 + t84 - s56
     #) * s56 - t256 * t489 + t270 - t491)) * t159 * t492) + t246 * t108
     # * t492 * (t25 * (t269 - t274 + t287) - t278 + t303) * t219 + t120
     # + propW34 * (t111 * (t488 - s26 - s16 + s34 - s13 + s56 + t226) +
     # t268 * (t128 * t258 + t267 * t245)) * t27
      t43 = t493 * t72 + t493 * t243 + t230 + t296 * (t24 * (-(t129 * t4
     #72 + t394 * t473) * epinv + t130 * (epinv * t242 - t418 * s34 + s5
     #6 * t419 - t16 * t372 + t421 + t422)) + t27 * (-t133 * (epinv * t3
     #58 + s16 + s26 - s34 - t21 - t353 + t356 - t359 + t60) - t135 * (e
     #pinv * t352 + s16 - s34 - t21 - t272 + t349 - t353) - t145 * (epin
     #v * (t16 * (-t316 + t369 - t368) - (t366 + s13) * s25 - (-t16 * (s
     #13 - s12) + s16 + t21) * s26 + (-t348 - t73 + s25) * s56 + t370 + 
     #t333 * s16) + t373 + t374 + t376 + t375) - t59 * (t10 * (t371 - t2
     #49 - t396) + t324 + epinv * (-t264 + t370 - t397 + t398 + t399) + 
     #(-t363 - s13 - t60 + s25) * s25 + s26 * t393 - t393 * s34 + (-t10 
     #* t395 + s25) * s56) - t488 * t468) + t31 * t43 - t225 * (t308 + t
     #334))
      t72 = t25 * t127
      t83 = -t72 + t66
      t109 = t460 * s34
      t120 = t104 * ((t331 - t333 + t93) * s25 - epinv * t381 + (-s26 + 
     #t355) * s46 + (t21 + s46) * s56 + t382) + t25 * (t113 * t378 - t11
     #6 * t350 - t117 * t248 + t279 * t469 + t351 * t467) + t315 * (-epi
     #nv * (-t16 * (s45 * t48 + t317 - t318) - (-s16 + s34 + s46) * s25 
     #- (-s16 + s46 - s25) * s26 - (-s16 + s46 + s25) * s56 - t249 + t30
     #2 + t316) - t16 * t321) + t68 * t329 + t92 * t350 - t304 * t447 * 
     #t492 + t109 * (epinv * t235 - s16 + s25 - s26 + s34 - t343 + t443 
     #- t60)
      t129 = -t25 * t451 + t442
      t191 = t25 * t458
      t226 = s34 * t234
      t230 = t112 * t25
      t232 = t25 * t53
      t233 = t126 * t25
      t240 = epinv * t456
      t242 = t25 * intHL0s25s260s56x1013D4eps1()
      t235 = t404 * (t108 * (t24 * (t247 * t9 + t261 * (t240 + t429) + t
     #55 * (s56 * t16 * t419 + epinv * t247 + s26 * t431 - s34 * t431 - 
     #t10 * t372 + t401)) + t27 * (t227 * (-t232 + t107) + t420 * t99 + 
     #t433 * t75 + t436 * t97 - t233 * t326) + t31 * (t248 * t95 + t261 
     #* (t191 - t445) - t279 * t480 - t351 * t438 - t381 * t82 - t425 * 
     #(-t439 * s45 - s46 * t265 + (s16 + s13 - s46 - s25 - t407 - s56) *
     # s56 + t409 + t411 - t437) - t71 * t378 + t129 * (t16 * (t256 * (t
     #439 + s34) + t369) + t415 + (-t16 * t241 + s26) * s34 + (t18 * (s3
     #4 - s13) - s16 + t272 + t405) * s56 + t371 + t413 + t414) + t226 *
     # t235 - t230 * t329)) + t31 * (t242 - t477) * t492 * t219)
      t241 = t25 * intHL0s25s26s34s56x1310D4eps1()
      t243 = t241 - t277
      t245 = t25 * intHs160s26s34s56x1013D4eps1()
      t246 = t25 * t298
      t120 = t404 * (t27 * (t14 * t261 + t15 * (epinv * t261 - s16 + s25
     # - s26 + s34 - t343 + t443 - t60)) + t108 * (t243 * t219 * s34 + t
     #246 * t304) * t492 * t31 + (t25 * (t27 * (t101 * t273 + t105 * t25
     #0) - t24 * t106 * t297) + t27 * (-t245 + t283) * t281 * s34) * t10
     #8 * t159) + t235 + t296 * (t27 * ((t466 * t470 + t77 * (epinv * t4
     #66 + s26 - t343 + t443 + t93)) * s34 + t159 * (-t121 * t250 - t273
     # * t32) - t25 * (t114 * t433 + t118 * t436 + t119 * t420) + t326 *
     # t64 + t83 * ((t10 * t33 + t60) * s13 - t301 * (t322 - s25)) + (-t
     #115 * t25 + t79) * (t10 * s34 * t256 + t16 * (s34 + s45 + s46) * s
     #25 + s13 * t336 - s16 * t292 + s25 * t48 + t317 + t252 * s56)) + t
     #31 * t120 + t122 * t24 * t297 * t159)
      t219 = t25 * intHLs16s25s26s34s56x1411D6eps1()
      t235 = t25 * intHLs16s25s26s34s56x1141D6eps1()
      t247 = t235 - t434
      t90 = -t210 + t90
      t210 = (t219 - t228) * s34
      t248 = t25 * intHLs160000x0122D4eps1()
      t249 = intHLs160000x0111D2eps1() * t17
      t250 = t25 * intHLs160000x0113D4eps1() - intHLs160000x0113D4eps0()
      t258 = 0.2e1_dp / 0.3e1_dp
      t264 = s26 + s16 - s34 + s56 - s25
      t270 = epinv * t33
      t271 = s25 * t131
      t272 = epinv * s25
      t273 = t144 * t291
      t281 = t141 * t457
      t292 = epinv * t473
      t297 = t12 * t482
      t302 = epinv * t4
      t304 = t108 * t24
      t308 = t30 * t108
      t309 = t25 * t105
      t310 = t300 * t447
      t312 = t448 * (t124 + s26 + s56)
      t313 = t108 * t492
      t275 = -t108 * t138 + t492 * (t108 * (t17 * t286 + t310 + t312) - 
     #t476 - t478) - t87 + t313 * (-t275 - t276 + t461) * s34 + t182 * t
     #492 * (-t274 + t287)
      t276 = (t280 + t282) * s34
      t286 = t123 * (s26 + s16 - s34 + s56 + t272)
      t314 = (t309 - t121) * t264
      t316 = s34 * t5
      t284 = t17 * t284
      t223 = t27 * (t108 * (t276 + t286) + t156 + t157) + t304 * (t314 +
     # t316 + t284) + t313 * epinv2 * (t20 * (-s25 * (propZ25 * t364 + t
     #10) - t299 * t41 - t342) - t22 * (t328 + t223))
      t299 = t24 * t145
      t313 = t27 * t130
      t34 = t131 * (s34 * (t259 * t27 + t31 * t39) - t299 * t364) + t25 
     #* (-t27 * (t255 * t40 + t267) + t31 * t34 * t291) - t31 * t484 + (
     #-t31 * t47 - t313) * t335 * epinv + t196 * t24 * t159
      t34 = t108 * t34 + t159 * t223 + t31 * t275 + t159 * (t27 * (t108 
     #* (s34 * (t302 - t266) + t103) - t190) + t304 * (-t452 * t128 + t1
     #02 + t229)) + t108 * (t20 * (-t12 * (t475 + t7) + t297) + t27 * (t
     #154 * t255 - t335 * t54 + t292 + t346)) + t31 * (t108 * (t147 * t3
     #64 - t153 * (t271 + s26 + s56) - t155 * t335 + t25 * t354 + t33 * 
     #t441 + t400 * (-s34 + t270) - t273 + t225 - t281) + t332 * t492) +
     # (-t13 * t325 - t142 * t254 - t143 * t364 + t151 * (-epinv * t254 
     #+ s34) - t59 * (epinv * t325 + s16 + s26 - s34 + s56)) * t108 * t2
     #4 + t304 * (-t133 * t25 + t361) + t308 * (t474 + t2) + t108 * (epi
     #nv * (-s34 * t462 * t492 - t139) + t492 * (t17 * t285 + t218 - t23
     #6 + t237 + t238) + t185 + t492 * (s34 * (-t289 + t485) - t298 * t3
     #00) * t25) * t31
      t40 = t25 * t106
      t191 = epinv * (-t304 * t55 * t264 - t31 * t315) + t20 * (t204 * t
     #108 - t12 * t14) + t24 * (t131 * t15 + t159 * t32 + (t131 * t77 + 
     #t159 * t283) * t108 * s34) + t25 * (t24 * (-t101 * t159 + t108 * (
     #-s34 * t115 + t114 * t364 + t126 * t264)) + (s34 * t113 + t112 * t
     #48 + t116 * t300 - t117 * t364) * t108 * t31 - t127 * t27 * t33 * 
     #t108) + t31 * (t492 * (s34 * t277 * t108 + t477) - t104 - t425) + 
     #t108 * (t20 * (t12 * (t429 + t122) + t29 * t445) + t27 * (t33 * t6
     #6 + t240 - t40)) + t31 * (t108 * (epinv * t140 + s34 * (-t241 * t4
     #92 - t460 - t71) - t300 * t92 + t364 * t95 - t48 * t68 - t191) - t
     #242 * t492) + t304 * (s34 * (-t245 * t159 + t470 + t79) + t264 * (
     #-t9 - t64) - t364 * t75)
      t223 = -t213 + t52
      t241 = t223 * t24
      t254 = t28 * t31
      t255 = 9.0_dp * t217
      t275 = t305 * t16
      t277 = 1.0_dp + t275
      t285 = t41 * propW16
      t298 = s25 + s26 - s34 + s56
      t317 = t305 * t41
      t318 = t317 + t11 + t16
      t245 = -t283 + t245
      t283 = t25 * t117
      t319 = t283 - t95
      t320 = t25 * t116
      t325 = -t320 + t92
      t326 = t24 * t159
      t328 = t7 + t8
      t329 = t1 * t318
      t322 = epinv2 * (t20 * (t285 * t33 - t10 + t11) + t22 * (t322 * pr
     #opW16 - 1.0_dp)) * t159
      t331 = epinv * t139
      t318 = t140 * (epinv * (t20 * t318 + t22 * t277) - t317 * t222)
      t333 = t108 * (t202 * t291 + t463 * t49) * t31
      t39 = t108 * (epinv * t313 * t298 + t299 * t48 * t131 + t318) + t1
     #59 * (t195 * t27 + t304 * (t309 * t264 + t196)) + t108 * (t20 * (t
     #12 * t429 + t172) + t240 * t27) + t159 * (t27 * (t108 * (s34 * (t2
     #80 + t282 - t266) + t103 + t286) - t190) + t304 * (-t121 * t264 + 
     #t102 + t229 + t284 + t316)) + t108 * (t20 * (-t12 * t328 + t297 + 
     #t329) + t24 * (t151 * (epinv * t48 + s34) + t25 * (-s34 * t128 * t
     #159 - t133 - t134) + t48 * (t142 + t143)) + t27 * (s34 * (t302 * t
     #159 + t259) + t298 * t54 + t292) + t31 * (s34 * t152 - t147 * (epi
     #nv * t251 + s25 + s26 + s56) - t150 * t251 - t153 * (t124 + s26 + 
     #s56 + t463) - t432 * t87 - t273) + t322) + t108 * (s34 * (epinv * 
     #t39 + t146) + t48 * (epinv * t47 + t155) - t138 - t186 - t331) * t
     #31 + t333
      t128 = t1 * t277
      t146 = t305 * t222
      t195 = t146 * t108 * t250
      t240 = t16 + epinv
      t251 = t305 * t10
      t264 = 1.0_dp + t251
      t252 = t16 * (s26 + s16 + s12) - t252
      t266 = -t388 + s25
      t273 = s26 + s16 - s34 + s56 + s12
      t277 = t10 * t273
      t273 = t16 * t273
      t280 = -t273 + s25
      t282 = t327 + t366 - t440 + t45
      t297 = -s13 + s25 + s26 + s45 + s56
      t298 = t10 * t365
      t299 = t385 * t42 + t298 + t85 + t93
      t302 = t277 + s25
      t309 = -t16 * (s16 - s34 + s13 - s45 - s56) + s25 - t56
      t313 = t345 * t385 + t365 * t41 + s25 + t85
      t316 = t16 * t383 + t21
      t244 = t21 + t244 + t293
      t327 = t42 * s26
      t305 = 12.0_dp * t305
      t333 = t305 + t11 + t16
      t334 = t337 + t354
      t335 = epinv * t471
      t336 = t291 * t400
      t38 = -t110 * t6 + t27 * (t108 * (-t159 * (t253 * t5 + t4 * t455) 
     #+ t347 * t313) - t111) + t108 * (t20 * (-t12 * t3 - t29 * t481 - t
     #333 * t474) + t24 * (-t130 * (t16 * (t131 * t365 + s26 - s34 + s56
     #) + s25) - t145 * (-t10 * (s16 - s34 + s12) - t16 * t450 + t271) -
     # t151 * (t16 * (epinv * s26 - s25) + t272 - t298 + t426 - t85 - t3
     #27) + t457 * (-t133 - t134) + t70 * (-t142 - t143) - t188) + t27 *
     # (t206 * t253 * t159 - t154 * t313 - t205) + t31 * (t147 * (epinv 
     #* t309 + s25 + t69) + t150 * t309 - t152 * t316 + t153 * (epinv * 
     #t244 + s25 + t69) + t244 * t87 - t25 * t334 + t291 * t441 + t281 -
     # t335) - t225 * (t20 * t333 + t22 * t264)) + t108 * (epinv * (-t31
     #6 * t49 + t336) + t38 * (-t202 + t144) + t483 + t484) * t31 - t308
     # * t2
      t188 = t26 * t199
      t244 = (t456 + t473) * epinv
      t253 = t48 * t425
      t271 = (t175 + t74) * t25
      t281 = t15 * t240
      t298 = t315 * t48
      t53 = t108 * (t24 * (-t25 * (t114 * t266 + t126 * t280) + t196 * t
     #159) + t318 + (-t298 - t109) * t31 * epinv - (t115 * t252 + t118 *
     # t299 + t119 * t282 + t127 * t302 + t267) * t25 * t27) + t159 * (t
     #27 * (t108 * (t276 + t103 + t286) + t156 + t157 - t190) + t304 * (
     #t314 + t102 + t284)) - t181 + t108 * (t20 * (t12 * (t107 + t429 + 
     #t482) + t329) + t24 * (t159 * t229 + t266 * t75 + t365 * t54 + t55
     # * (t57 - t277) + t59 * (-t273 + t272) + t64 * t280) + t27 * (s34 
     #* t259 + t25 * (-s34 * t486 * t159 - t53) + t252 * t79 + t282 * t9
     #9 + t299 * t97 + t302 * t66 + t244) + t31 * (-t104 * (t136 * t297 
     #+ s25) - t226 - t331 - t253) + t322) + t281 * t24 + t108 * (epinv 
     #* (t200 + t201) - t138 - t167 - t168 + t183 + t185 + t271) * t31
      t97 = s16 * (-t248 + t288) + t220 + t249
      t38 = -t10 * (t108 * (t27 * (t40 - t122) + t31 * (t265 * t325 + t2
     #97 * t82 - t319 * t80)) + t326 * (t245 * t108 * s34 + t173 - t32))
     # + t16 * t53 + t18 * (t108 * (-t254 - t188) + t217 * (t264 * t474 
     #+ t2 - t3 + t481)) + 24.0_dp * t195 + t41 * (t108 * t209 * (t273 +
     # s25) * t27 + t304 * (-t212 + t78 + t77) * s34 + t217 * (t128 - t4
     #82 - t429 + t174 - t107) + t241 + t216) - t38 - 12.0_dp * t108 * t
     #222 * propW16 * t97 + 9.0_dp * t217 * t199
      t53 = t25 * t125
      t99 = epinv * (-t147 - t153 + t49) - t150 + t152 - t87
      t118 = epinv * t104 + t82
      t119 = s26 - s56 - s25
      t125 = s16 + s56
      t127 = t21 + t125 + t60 + t327
      t183 = s16 + s34 + s56 + t37 + t45
      t201 = -s13 + s45 + s46 + s25 + t88
      t209 = t265 + t85 - t84
      t226 = t21 + s26 - s16 + s34 - s56 + t60
      t229 = s56 + t45 + t60 + t355
      t252 = t131 * s12
      t264 = t16 * (t252 + s26 + s25)
      t266 = epinv * t468
      t193 = t304 * (t142 * t70 + t193 + t194 - t206 + t207)
      t194 = t108 * (-t12 * (t475 + t4 + t5) + t204) * t20
      t29 = t27 * (t159 * (t384 - t190 + t156 + t157) - t208 * t226) + t
     #108 * (t20 * (t12 * t479 - t12 * t482 + t29 * t484) + t27 * (t13 *
     # (t33 + t60) + t133 * t137 + t135 * t457 + t143 * t224 + t145 * (e
     #pinv * t489 + t264) + t154 * t226 + t59 * (t21 + t132 + t270) + t2
     #66) + t31 * (-t119 * t155 + t137 * t337 + t144 * t229 + t147 * (ep
     #inv * (t385 - t85 + t84) + t264) + t150 * (t385 - t85 + t231) - t1
     #53 * (t16 * (epinv * t201 - s26) - t21) + t17 * t354 + t25 * (t141
     # + t289) - t400 * (-t60 * t25 - s16 - s56 - t21 + t272 - t327) - t
     #441 * t63 + t49 * (epinv * t51 - s25)) + t292 * t24) + t108 * (epi
     #nv * (-t471 - t139 + t140 + t200) + t47 * (-s26 * t240 - t124 + t4
     #49) - t138 + t184 + t185 - t28 - t202 * t229) * t31 + t193 + t194
      t35 = t108 * (t83 * t35 * t27 + t31 * t265 * (-t230 + t68))
      t83 = s25 * t9
      t84 = t108 * (t27 * (t215 + t214) + t31 * (s34 * (-t219 + t228) + 
     #t211 * t446 + t70 * (t434 - t235))) + t217 * (t483 + t2 - t7 + t83
     #)
      t132 = t26 * t6 * t108 + t217 * (t475 - t479 + t482 + t484 + t1 + 
     #t4 + t5)
      t29 = t10 * t35 - t16 * (t24 * (-t205 * t108 + t179) + t27 * (t25 
     #* (t108 * (-t115 * t364 + t100 + t98) + t114 + t126) - t15 - t64 -
     # t75) + t31 * (t108 * (t25 * (-t113 * t364 - t119 * (t116 + t117) 
     #+ t127 * t451 + t183 * t467) + t480) + t425) + t108 * (t24 * (-t17
     #1 + t169 - t170) + t27 * (t364 * t79 - t160 - t163 - t166) + t31 *
     # (t104 * (epinv * t209 - s25 + s26 - s56) - t25 * t469 + t315 * (t
     #57 - t36) + t364 * t71 - t438 * t183 - t442 * t127 + t82 * t209 + 
     #t87 * t201 + t119 * (t92 + t95) - t109) + t30 * (-t2 - t483)) + t2
     #6 * (-t108 * t7 + t9)) + t18 * t132 + t41 * t84 + t110 * (t3 + t10
     #7) + t24 * (t111 * t131 + t487) + t29 + t108 * (t197 * t24 - t53 *
     # t27) * t159 - (t481 + t462 + t461) * t108 * t31 - t255 * t6
      t30 = s13 - s45 - s46 - s56 - s25 + s12
      t35 = s13 - s45 - s46 - s56 - s25 - s12
      t36 = s13 - s25 - s45 - s46 - s56
      t57 = t48 + t263
      t45 = -t16 * (s16 - s34 - s25 + s12) + s56 - t45
      t79 = -t257 - s26 - t85 + t311
      t84 = t10 * t385 + t125 * t16 + s25
      t110 = t16 * t48
      t115 = -s26 + t85 + t355
      t44 = t293 + t69 + t44
      t125 = t16 * (s26 + s45 + s46 + s56) + s25 - t311
      t126 = t20 * (-t317 + t11 + t16) + t22 * (1.0_dp - t275)
      t127 = -t16 * t357 + s25
      t73 = -t363 + t440 - t73 + s25
      t85 = -t257 + t148 - t85 + t262
      t131 = 1.0_dp + t136
      t132 = t20 * (t305 - t11 - t16) + t22 * (-1.0_dp + t251)
      t47 = t296 * (t140 * (epinv * t132 - 12.0_dp * t146) + t159 * (t24
     # * (-t123 * (t25 * t36 - t252) + t102 + t196) + t30 * t27 * (t190 
     #- t156 - t157)) + t31 * (t147 * (epinv * t85 - t10 * t295 + t21 - 
     #t69) + t150 * t85 + t155 * t57 + t47 * (epinv * t57 - t110) + (s25
     # * t41 + t10 * (s34 + s56) - t363 - t60) * (-t202 + t144) - t138 +
     # t200 * t17) + t142 * t24 * t84 + t27 * (s12 + s13 + s16 + s26 - s
     #34 - s45 - s46) * (-t53 + t103) * t307 - t347 * t27 * t46)
      t57 = t296 * (t31 * (t25 * (-t112 * t44 - t113 * t73 - t115 * t116
     # - t117 * t45 + t492 * (t289 * t291 - t269)) + t492 * (-t291 * t46
     #1 + t462 * t464 + t278 - t310 - t312) - t303 * t492) - (t100 + t17
     #6 + t177 + t81 + t98) * t25 * t27 + t326 * (-t121 * t36 - t122 * t
     #35 - t30 * t32))
      t30 = t404 * (t108 * (-t126 * t474 + t27 * t479 + t31 * (t246 * t3
     #00 * t492 + t2 + t28 - t481)) + t268 * (t101 * t30 + t105 * t36 + 
     #t106 * t35) * t24 + t31 * t492 * t239)
      t11 = t30 + t57 + t296 * (t24 * (s34 * (epinv * t259 + t459) + t13
     # * t33 + t323 * t9 - t54 * (s26 + s16 - s34 + s56 - s12) + t55 * (
     #epinv * t323 - t56 - t93)) + t27 * (t59 * (t37 - t270) + t76 * (-t
     #114 * t25 + t75) + t165 + t160 * t131) + t31 * (t104 * (epinv * t7
     #9 - s26 + s56 + t21) + t45 * t95 + t68 * t44 + t71 * t73 + t79 * t
     #82 + t87 * (t16 * t339 - t256) + t92 * t115 + t492 * t25 * (t274 -
     # t287)) + epinv2 * (t20 * (t285 * t48 + t11 + t16) + t22 * (t110 *
     # propW16 + 1.0_dp)) * t492 - t225 * t126) + t404 * (t24 * (t108 * 
     #(-t232 + t107 + t429 + t482 + t244) + t487) + t27 * (t108 * (t166 
     #+ t161 + t162 + t163 + t266) + t14 + t281) + t31 * (t492 * (t108 *
     # (-t218 + t236 - t237 - t238) - t332) - t335 * t108 + t108 * t334 
     #* t137))
      t22 = t404 * (t31 * (t108 * (t25 * (-t290 * t467 + t458 + t469) + 
     #t290 * t438 + t336 - t445 - t480 + t483 + t484) + t492 * (t243 * t
     #108 * s34 + t242 - t477) + t108 * t129 * (s34 + t110)) + t108 * (s
     #34 * t470 + t6) * t27)
      t30 = t296 * (t31 * (t247 * t48 + t210) + t222 * t97 * propW16)
      t35 = t109 + t298
      t11 = -t10 * t404 * (t27 * (t108 * (t149 - t62) * t63 - t52 + t213
     # + t108 * (t212 - t78) * s34) + t108 * t90 * t119 * t31) + t16 * t
     #296 * (t24 * (t205 + t3) + t27 * t328) + t258 * (t47 + t404 * (t10
     #8 * (t1 * t132 + t24 * (t130 * (-t16 * (-t252 + t270) + s25) + t15
     #1 * (epinv * t84 - s26 - s34 + s56 - t56 - t93) - t206 + t4 + t5))
     # + t108 * (t127 * t143 + t135 * (epinv * t10 - t42) + t145 * (epin
     #v * t127 + t21 + t56 + t69) + t457 * (t133 + t134) + t189) * t27 +
     # (t125 * t152 + t153 * (t10 * (epinv * s13 + s56) - t16 * (epinv *
     # (t256 + t88) - s26) + t444) + t49 * (epinv * t125 - t93) + t185 +
     # t186 - t331) * t108 * t31 + t111 * t131 * t24)) + 0.4e1_dp / 0.3e
     #1_dp * t11 + 0.8e1_dp / 0.3e1_dp * t22 + 0.16e2_dp / 0.3e1_dp * t2
     #96 * t31 * t35 + 16.0_dp * t146 * t296 * t250 - 8.0_dp * t30
      t22 = t260 + s25 + t56
      t30 = epinv * t70
      t36 = t260 + s25
      t37 = t6 - t107
      t17 = t24 * (s25 * t111 - t135 * t137 - t145 * (t21 - t30 + t69 + 
     #t56) - t198 * t17 + t192 + t205 - t206 + t207) + t31 * (epinv * (-
     #t330 - t471) - t141 * t17 - t153 * (-epinv * t36 + s25 + t69) + t2
     #5 * (-t337 - t354) - t187) + t188 + t24 * (-t151 * (t93 - t30 + s2
     #6 + s34 - s56 + t56) + t46 * (t347 - t154) + t3 + t4 + t5) + t31 *
     # (-t147 * (-epinv * t22 + s25 + t69) + t150 * t22 + t291 * (epinv 
     #* t400 + t441) + t51 * (epinv * t49 + t152) + t87 * t36 + t2 + t20
     #3 - t28 - t474 - t481 + t483 + t484)
      t9 = t404 * (struc2PM * (t16 * t191 - t18 * (t108 * (t26 * t8 + t2
     #54) + t217 * (-t474 - t475 + t482 - t2 - t7)) + t41 * (t217 * (s25
     # * t14 + t1 - t122 - t429 + t445) + t241 + t216) + t34 + t255 * t8
     #) + (-t16 * (t108 * (t31 * (s34 * (-t113 * t25 + t71) + t104 * (t1
     #24 + t463) + t300 * t325 + t364 * t319 + t432 * t82 + t48 * (-t230
     # + t68) + t28) + t40 * t27 + s34 * (t245 * t159 - t77) * t24) + t2
     #6 * (-t108 * t122 + t13 + t14 + t9) - (t15 + t55 + t59) * t24 * ep
     #inv + t326 * (t173 - t32)) + t18 * t217 * (t128 + t2 + t7 + t8 - t
     #429 - t482) - t41 * ((s16 * (t288 - t474 - t225 - t248) + t220 + t
     #249) * t108 * t222 * propW16 + t217 * (-t174 + t122) - t216) + 12.
     #0_dp * t195 + t39) * (struc3PM + struc8PM) - (-t118 * t16 + t99) *
     # t31 * struc61PM * t108)
      t13 = struc60PM * (t10 * t31 * t118 - t16 * (t31 * t99 + (t25 * (t
     #101 + t105 + t106) - t121 - t122 - t32) * t159 * t27) - t159 * (t2
     #4 * (t25 * (t123 + t158) - t156 - t157) + (t53 - t103) * t159 * t2
     #7)) + struc9PM * (t10 * t31 * (t265 * (t320 - t92) - t50 * t82 + t
     #80 * (t283 - t95)) + t16 * (s34 * (t234 * t31 + t27 * t77) + t24 *
     # (-t166 - t232 - t164 - t165 + t171 + t83) + t24 * (s25 * (t13 + t
     #14 + t180) + t25 * (t81 + t98 + t100 + t178) + t55 * t58 + t59 * t
     #61 + t65 * (t233 - t64) + t67 * (-t66 + t72) - t163) + t31 * (epin
     #v * t35 - t104 * (t136 * t50 + s25) - t167 - t168 + t253 + t271) +
     # t26 * t37) - t18 * t221 * t199 - t41 * (t24 * (s25 * t223 + t214 
     #+ t215) + t221 * t37 - t216 * s25) + t17)
      t10 = (t10 * t296 * t306 * t222 * t250 + t120 * t258 + t16 * t404 
     #* (t108 * (t27 * ((s13 * t323 - t301 * t33 + t324) * (-t149 + t62)
     # + (-t212 + t78) * t279 * s34) + t31 * (t321 * t247 + t90 * (-t16 
     #* s25 * (-t294 + s13) - s46 * t446 + t360) + t210 * t227)) + t27 *
     # (-t213 + t52) * t261 + t490 * t108 * (s16 * (t248 - t288) - t220 
     #- t249) * t222) - t43 / 3.0_dp + t296 * t8 * t27 * t261) * struc1P
     #M + t11 * struc7PM
      result = t258 * t404 * (t108 * t13 + (-t16 * (t24 * (t179 - t180) 
     #+ t108 * (t20 * (-t12 * t6 - t172) + t24 * (-t59 * t61 + t169 - t1
     #70 - t171) + t27 * (-t166 - t160 - t161 - t162 - t163 - t164 - t16
     #5) + t31 * (-t82 * t86 - t87 * t89 + t92 * t94 + t95 * t96 + t167 
     #+ t168)) + t181 + t182 * (t81 + t98 + t100 + t176 + t177 + t178) *
     # t27 + t108 * (t104 * (-epinv * t86 + s25 + s26 - s56) - t25 * (t1
     #16 * t94 + t117 * t96 + t175 + t74) + t28) * t31) + t18 * t217 * (
     #t1 + t3 + t4 + t5 - t7 - t8) + t41 * (t27 * (t108 * t214 + t108 * 
     #t215 - t213 + t52) + (t174 + t2 - t6) * propZ25 * t19 + t216 * t10
     #8 * t91) + t23) * struc10PM + t38 * struc5PM + t29 * struc6PM) + 0
     #.4e1_dp / 0.3e1_dp * t9 + t10

           ampNonresonantLightFullReC4PM = result
       end function ampNonresonantLightFullReC4PM

       function ampNonresonantLightFullReC4PP()
           implicit none
           complex(dp) :: ampNonresonantLightFullReC4PP

           complex(dp) :: result

      result = 0.0_dp

           ampNonresonantLightFullReC4PP = result
       end function ampNonresonantLightFullReC4PP

       function ampNonresonantLightFullReC7MM()
           implicit none
           complex(dp) :: ampNonresonantLightFullReC7MM

           complex(dp) :: result

      result = 0.0_dp

           ampNonresonantLightFullReC7MM = result
       end function ampNonresonantLightFullReC7MM

       function ampNonresonantLightFullReC7MP()
           implicit none
           complex(dp) :: ampNonresonantLightFullReC7MP

           complex(dp) :: result

      result = 0.0_dp

           ampNonresonantLightFullReC7MP = result
       end function ampNonresonantLightFullReC7MP

       function ampNonresonantLightFullReC7PM()
           implicit none
           complex(dp) :: ampNonresonantLightFullReC7PM

           complex(dp) :: result

      result = 0.0_dp

           ampNonresonantLightFullReC7PM = result
       end function ampNonresonantLightFullReC7PM

       function ampNonresonantLightFullReC7PP()
           implicit none
           complex(dp) :: ampNonresonantLightFullReC7PP

           complex(dp) :: result

      result = 0.0_dp

           ampNonresonantLightFullReC7PP = result
       end function ampNonresonantLightFullReC7PP

